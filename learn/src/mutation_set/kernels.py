"""Kernel layout and assembly from nmr_extract SDK.

Each kernel is a T2 (L=2, 5 components) vector per atom.  The layout
is a named registry — no magic indices.  Assembly uses the SDK's
.as_block() and named accessors throughout.

Kernel groups:
    ring_type   — BS/HM/Disp per 8 ring types = 24
    bond_cat    — MC/MopacMC per 5 bond categories = 10
    total       — calculator totals = 6
    efg         — electric field gradient = 6
    TOTAL = 46
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from nmr_extract import RingType, N_RING_TYPES
from .config import Config


# ── Kernel names ─────────────────────────────────────────────────────

_RING_TYPE_NAMES = [rt.name for rt in RingType]

_CALC_PREFIXES = {
    "biot_savart":  "BS",
    "haigh_mallion": "HM",
    "dispersion":    "Disp",
}

_BOND_CAT_NAMES = [
    "backbone_total", "sidechain_total", "aromatic_total",
    "CO_nearest", "CN_nearest",
]

_TOTAL_SHORT = {
    "mcconnell":          "MC_total",
    "coulomb":            "Coulomb_total",
    "hbond":              "HBond_total",
    "ring_susceptibility": "RingSusc_total",
    "mopac.coulomb":      "MopacCoulomb_total",
    "mopac.mcconnell":    "MopacMC_total",
}

_EFG_SHORT = {
    "coulomb.efg_backbone":       "EFG_bb",
    "coulomb.efg_aromatic":       "EFG_aro",
    "mopac.coulomb.efg_backbone": "MopacEFG_bb",
    "mopac.coulomb.efg_aromatic": "MopacEFG_aro",
    "apbs.efg":                   "APBS_EFG",
    "delta.apbs.delta_efg":       "DeltaAPBS_EFG",
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

        return KernelLayout(
            names=names,
            n_kernels=len(names),
            ring_type_end=ring_type_end,
            bond_cat_end=bond_cat_end,
            total_end=total_end,
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
    for calc_name in ["biot_savart", "haigh_mallion", "dispersion"]:
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
                    "hbond.shielding", "ring_susceptibility"]:
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
