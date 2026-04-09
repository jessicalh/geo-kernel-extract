"""SecondaryLoader: shared protein iteration with ring-type tagging.

Every secondary analysis tool iterates proteins and needs:
  - matched atom indices (from WT-ALA delta)
  - per-atom bitmask of which ring types the atom "sees"
  - protein-level ring type inventory
  - element, distance to nearest removed ring

The loader yields these per protein so each tool picks what it needs.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mutation_set.config import Config, load_config
from mutation_set.dataset import list_proteins
from mutation_set.kernels import KernelLayout, assemble_kernels, normalize_kernels

# Ring type names in index order (matches nmr_extract.RingType enum)
RING_TYPE_NAMES = [
    "PHE", "TYR", "TRP_benzene", "TRP_pyrrole", "TRP_perimeter",
    "HIS", "HID", "HIE",
]
N_RING_TYPES = len(RING_TYPE_NAMES)

# Stratum definitions: each maps to a predicate on atom ring-type bitmask
# and protein ring-type set.  Built from config strata list.
_STRATUM_PREDICATES = {
    "hie_only":  lambda mask: (mask & (1 << 7)) != 0 and (mask & ~(1 << 7)) == 0,
    "phe_only":  lambda mask: (mask & (1 << 0)) != 0 and (mask & ~(1 << 0)) == 0,
    "tyr_only":  lambda mask: (mask & (1 << 1)) != 0 and (mask & ~(1 << 1)) == 0,
    "trp_only":  lambda mask: _is_trp_only(mask),
    "no_hie":    lambda mask: (mask & (1 << 7)) == 0 and mask != 0,
    "all":       lambda mask: True,
}


def _is_trp_only(mask: int) -> bool:
    """TRP-only: sees TRP_benzene, TRP_pyrrole, or TRP_perimeter, nothing else."""
    trp_bits = (1 << 2) | (1 << 3) | (1 << 4)
    return (mask & trp_bits) != 0 and (mask & ~trp_bits) == 0


@dataclass
class ProteinRecord:
    """Per-protein data yielded by the loader."""
    pid: str
    protein: object          # nmr_extract Protein
    matched_idx: np.ndarray  # indices of matched atoms
    element: np.ndarray      # element numbers for matched atoms
    ring_dist: np.ndarray    # distance to nearest removed ring
    atom_ring_masks: np.ndarray   # (M,) int — bitmask per matched atom
    protein_ring_types: set[int]  # which ring types exist in this protein
    has_mopac: bool


def _build_atom_ring_masks(protein, matched_idx: np.ndarray) -> np.ndarray:
    """Build per-atom bitmask of ring types seen, for matched atoms only."""
    M = len(matched_idx)
    masks = np.zeros(M, dtype=np.int32)
    rc = protein.ring_contributions
    if rc is None or rc.n_pairs == 0:
        return masks

    # Build lookup: atom_index → set of ring types
    atom_types: dict[int, int] = {}
    for i in range(rc.n_pairs):
        ai = int(rc.atom_index[i])
        rt = int(rc.ring_type[i])
        atom_types[ai] = atom_types.get(ai, 0) | (1 << rt)

    for j, ai in enumerate(matched_idx):
        masks[j] = atom_types.get(int(ai), 0)
    return masks


def _protein_ring_set(protein) -> set[int]:
    """Which ring types exist in this protein."""
    rg = protein.ring_geometry
    if rg is None or rg.n_rings == 0:
        return set()
    return set(int(rt) for rt in rg.ring_type)


def iter_proteins(cfg: Config, require_mopac: bool = False,
                  max_proteins: int = 0):
    """Yield ProteinRecord for each protein with delta data.

    Args:
        cfg: Config from calibration.toml.
        require_mopac: skip proteins without MOPAC data.
        max_proteins: 0 = all, >0 = limit (for quick testing).
    """
    from nmr_extract import load

    features_dir = cfg.paths.features
    proteins = list_proteins(features_dir)
    if max_proteins > 0:
        proteins = proteins[:max_proteins]

    for pid in proteins:
        try:
            p = load(features_dir / pid)
        except (FileNotFoundError, ValueError):
            continue

        if p.delta is None:
            continue

        has_mopac = p.mopac is not None
        if require_mopac and not has_mopac:
            continue

        matched_idx = np.where(p.delta.scalars.matched_mask)[0]
        if len(matched_idx) < cfg.data.min_matched_atoms:
            continue

        yield ProteinRecord(
            pid=pid,
            protein=p,
            matched_idx=matched_idx,
            element=p.element[matched_idx],
            ring_dist=p.delta.scalars.nearest_removed_ring_dist[matched_idx],
            atom_ring_masks=_build_atom_ring_masks(p, matched_idx),
            protein_ring_types=_protein_ring_set(p),
            has_mopac=has_mopac,
        )


def assign_stratum(atom_ring_mask: int, strata: list[str]) -> str | None:
    """Return the first matching stratum name for an atom, or None."""
    for s in strata:
        pred = _STRATUM_PREDICATES.get(s)
        if pred and pred(atom_ring_mask):
            return s
    return None


def nearest_ring_type(protein, atom_idx: int) -> int | None:
    """Ring type of the nearest ring to this atom, or None."""
    rc = protein.ring_contributions
    if rc is None:
        return None
    atom_rc = rc.for_atom(atom_idx)
    if atom_rc.n_pairs == 0:
        return None
    nearest = int(np.argmin(atom_rc.distance))
    return int(atom_rc.ring_type[nearest])


def ridge_fit(X: np.ndarray, y: np.ndarray,
              lam: float) -> tuple[np.ndarray, float]:
    """Ridge regression, returns (predictions, R²)."""
    XtX = X.T @ X + lam * np.eye(X.shape[1])
    w = np.linalg.solve(XtX, X.T @ y)
    pred = X @ w
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean(axis=0, keepdims=True)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
    return pred, r2


ELEMENT_NAMES = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}


def setup_sdk(cfg: Config):
    """Ensure nmr_extract SDK is importable."""
    sdk_path = str(cfg.paths.sdk)
    if sdk_path not in sys.path:
        sys.path.insert(0, sdk_path)
