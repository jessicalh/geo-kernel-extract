"""Kernel layout and assembly from nmr_extract SDK.

Each kernel is a T2 (L=2, 5 components) vector per atom.  The layout
is a named registry — no magic indices.  Assembly uses the SDK's
.as_block() and named accessors throughout.

Kernel groups:
    ring_type   — BS/HM/Disp/PQ per 8 ring types = 32
    bond_cat    — MC/MopacMC per 5 bond categories = 10
    total       — calculator totals = 7
    efg         — electric field gradient = 6
    per_ring    — individual ring BS/HM/Chi/PQ/HM_H/DispChi for top-K nearest = 6 × K
    TOTAL = 55 + 6K  (K from config, default 6 → 91 kernels)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from nmr_extract import RingType, N_RING_TYPES
from .config import Config


# ── Kernel names ─────────────────────────────────────────────────────

_RING_TYPE_NAMES = [rt.name for rt in RingType]

_CALC_PREFIXES = {
    "biot_savart":   "BS",
    "haigh_mallion": "HM",
    "dispersion":    "Disp",
    "pi_quadrupole": "PQ",
}

_BOND_CAT_NAMES = [
    "backbone_total", "sidechain_total", "aromatic_total",
    "CO_nearest", "CN_nearest",
]

_TOTAL_SHORT = {
    "mcconnell":              "MC_total",
    "coulomb":                "Coulomb_total",
    "hbond":                  "HBond_total",
    "ring_susceptibility":    "RingSusc_total",
    "pi_quadrupole.shielding": "PQ_total",
    "mopac.coulomb":          "MopacCoulomb_total",
    "mopac.mcconnell":        "MopacMC_total",
}

_EFG_SHORT = {
    "coulomb.efg_backbone":       "EFG_bb",
    "coulomb.efg_aromatic":       "EFG_aro",
    "mopac.coulomb.efg_backbone": "MopacEFG_bb",
    "mopac.coulomb.efg_aromatic": "MopacEFG_aro",
    "apbs.efg":                   "APBS_EFG",
    "delta.apbs.delta_efg":       "DeltaAPBS_EFG",
}

_PER_RING_CALC_SHORT = {
    "bs":       "BS",
    "hm":       "HM",
    "chi":      "Chi",
    "pq":       "PQ",
    "hm_H":     "HM_H",
    "disp_chi": "DispChi",  # computed: disp_scalar × chi.T2
}


@dataclass(frozen=True)
class KernelLayout:
    """Named kernel registry with index lookup."""
    names: list[str]
    n_kernels: int
    # Group boundaries for slicing
    ring_type_end: int
    bond_cat_end: int
    total_end: int
    efg_end: int
    per_ring_k: int            # number of nearest rings
    per_ring_calcs: list[str]  # calculator names for per-ring kernels

    @staticmethod
    def from_config(cfg: Config) -> KernelLayout:
        names = []
        # Ring-type: calculator × ring type
        for calc in cfg.kernels.ring_calculators:
            prefix = _CALC_PREFIXES[calc]
            for rt in _RING_TYPE_NAMES:
                names.append(f"{prefix}_{rt}")
        ring_type_end = len(names)

        # Bond-category: MC × 5 categories, MopacMC × 5
        for _ in cfg.kernels.bond_calculators:
            for cat in _BOND_CAT_NAMES:
                names.append(f"MC_{cat}")
        for _ in cfg.kernels.bond_calculators_mopac:
            for cat in _BOND_CAT_NAMES:
                names.append(f"MopacMC_{cat}")
        bond_cat_end = len(names)

        # Calculator totals
        for calc in cfg.kernels.total_calculators:
            names.append(_TOTAL_SHORT[calc])
        total_end = len(names)

        # EFG
        for calc in cfg.kernels.efg_calculators:
            names.append(_EFG_SHORT[calc])
        efg_end = len(names)

        # Per-ring: individual ring T2 for top-K nearest rings.
        # Each ring gets one kernel per calculator (BS, HM, Chi, PQ, HM_H).
        # Ordered by distance: ring0 is nearest, ring5 is 6th nearest.
        # The gating handles sparsity — atoms with fewer than K rings
        # have zero magnitude in the distant slots.
        per_ring_calcs = cfg.kernels.per_ring_calculators
        K = cfg.kernels.per_ring_k
        for ri in range(K):
            for calc in per_ring_calcs:
                prefix = _PER_RING_CALC_SHORT[calc]
                names.append(f"{prefix}_ring{ri}")

        return KernelLayout(
            names=names,
            n_kernels=len(names),
            ring_type_end=ring_type_end,
            bond_cat_end=bond_cat_end,
            total_end=total_end,
            efg_end=efg_end,
            per_ring_k=K,
            per_ring_calcs=per_ring_calcs,
        )

    def index(self, name: str) -> int:
        return self.names.index(name)


def _resolve(protein, dotpath: str):
    """Navigate p.mopac.coulomb.shielding etc."""
    obj = protein
    for part in dotpath.split("."):
        obj = getattr(obj, part, None)
        if obj is None:
            return None
    return obj


def assemble_kernels(protein, idx: np.ndarray, layout: KernelLayout
                     ) -> np.ndarray:
    """Build (M, n_kernels, 5) kernel T2 array using SDK accessors.

    Returns raw (unnormalized) kernels.
    """
    M = len(idx)
    K = np.zeros((M, layout.n_kernels, 5), dtype=np.float64)
    p = protein
    col = 0

    # Ring-type T2: each calculator provides (N, 8, 5) via .as_block()
    for calc_name in ["biot_savart", "haigh_mallion", "dispersion", "pi_quadrupole"]:
        group = getattr(p, calc_name)
        block = group.per_type_T2.as_block()[idx]  # (M, 8, 5)
        K[:, col:col + N_RING_TYPES, :] = block
        col += N_RING_TYPES

    # Bond-category T2: (N, 5, 5) via .as_block()
    K[:, col:col + 5, :] = p.mcconnell.category_T2.as_block()[idx]
    col += 5
    if p.mopac:
        K[:, col:col + 5, :] = p.mopac.mcconnell.category_T2.as_block()[idx]
    col += 5

    # Calculator totals
    for dotpath in ["mcconnell.shielding", "coulomb.shielding",
                    "hbond.shielding", "ring_susceptibility",
                    "pi_quadrupole.shielding"]:
        obj = _resolve(p, dotpath)
        K[:, col, :] = obj.T2[idx]
        col += 1
    for dotpath in ["mopac.coulomb.shielding", "mopac.mcconnell.shielding"]:
        obj = _resolve(p, dotpath)
        if obj is not None:
            K[:, col, :] = obj.T2[idx]
        col += 1

    # EFG T2
    for dotpath in ["coulomb.efg_backbone", "coulomb.efg_aromatic"]:
        K[:, col, :] = _resolve(p, dotpath).T2[idx]
        col += 1
    for dotpath in ["mopac.coulomb.efg_backbone", "mopac.coulomb.efg_aromatic"]:
        obj = _resolve(p, dotpath)
        if obj is not None:
            K[:, col, :] = obj.T2[idx]
        col += 1
    obj = _resolve(p, "apbs.efg")
    if obj is not None:
        K[:, col, :] = obj.T2[idx]
    col += 1
    obj = _resolve(p, "delta.apbs.delta_efg")
    if obj is not None:
        K[:, col, :] = obj.T2[idx]
    col += 1

    # Per-ring T2 RESIDUALS: individual ring T2 minus per-type mean.
    #
    # For single-ring-per-type proteins (the common case), the per-type
    # sum IS the single ring's contribution, so the residual is zero —
    # correctly indicating no new information beyond the summed kernel.
    # For multi-ring proteins, the residual captures the angular structure
    # unique to each ring that the type sum collapses.
    #
    # The gating handles zeros naturally: zero residual = zero magnitude
    # = gate suppresses the kernel.
    #
    # Map calculator names to per-type-sum column blocks (0-31):
    #   bs  → columns 0-7   (BS per-ring-type)
    #   hm  → columns 8-15  (HM per-ring-type)
    #   pq  → columns 24-31 (PQ per-ring-type)
    #   chi, hm_H, disp_chi → no per-type decomposition available, use raw
    _PERRING_CALC_TO_TYPE_BLOCK = {"bs": 0, "hm": 8, "pq": 24}

    rc = p.ring_contributions
    if rc is not None and rc.n_pairs > 0:
        # Count rings per type for this protein (for computing per-type mean)
        rg = p.ring_geometry
        rings_per_type = np.zeros(N_RING_TYPES, dtype=int)
        if rg is not None:
            for rt in rg.ring_type.astype(int):
                rings_per_type[rt] += 1

        for j in range(M):
            atom_rc = rc.for_atom(idx[j])
            if atom_rc.n_pairs == 0:
                continue
            order = np.argsort(atom_rc.distance)
            for ri_slot, ri in enumerate(order[:layout.per_ring_k]):
                ring_type = int(atom_rc.ring_type[ri])
                for ci, calc_name in enumerate(layout.per_ring_calcs):
                    kernel_col = col + ri_slot * len(layout.per_ring_calcs) + ci

                    if calc_name == "disp_chi":
                        # Dispersion-weighted Chi: disp_scalar × chi.T2
                        # scalar × L=2 = L=2, equivariant by construction.
                        # Encodes angular position weighted by vertex contact.
                        K[j, kernel_col, :] = (atom_rc.disp_scalar[ri]
                                               * atom_rc.chi.T2[ri])
                    elif calc_name in _PERRING_CALC_TO_TYPE_BLOCK:
                        # Subtract per-type mean: type_sum / n_rings_of_type
                        ring_t2 = getattr(atom_rc, calc_name).T2[ri]
                        type_col = _PERRING_CALC_TO_TYPE_BLOCK[calc_name] + ring_type
                        n_of_type = max(rings_per_type[ring_type], 1)
                        type_mean = K[j, type_col, :] / n_of_type
                        K[j, kernel_col, :] = ring_t2 - type_mean
                    else:
                        # Chi, hm_H: no per-type decomposition, use raw
                        K[j, kernel_col, :] = getattr(atom_rc, calc_name).T2[ri]
    col += layout.per_ring_k * len(layout.per_ring_calcs)

    assert col == layout.n_kernels, f"Built {col}, expected {layout.n_kernels}"
    return K


def normalize_kernels(kernels: np.ndarray, std_floor: float
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Per-protein normalization.  Returns (normalized, scales).

    Scales are the pre-normalization std per kernel — the magnitude
    information the MLP needs to bridge per-protein and global.
    """
    n_kernels = kernels.shape[1]
    scales = np.zeros(n_kernels, dtype=np.float64)
    for k in range(n_kernels):
        s = kernels[:, k, :].std()
        if s > std_floor:
            kernels[:, k, :] /= s
            scales[k] = s
    return kernels, scales
