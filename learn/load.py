"""
Data loading for learn/.

Loads extracted .npy features for WT/ALA pairs.
Handles the SphericalTensor packing: [T0, T1[3], T2[5]] = 9 doubles.
Provides T2 extraction and scaling.
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

FEATURES_BASE = Path(__file__).parent / "features"

def features_dir(run: str = "FirstExtraction") -> Path:
    """Return features directory for a named extraction run."""
    return FEATURES_BASE / run

# Default for backward compat — callers should use features_dir(run) explicitly
FEATURES_DIR = FEATURES_BASE / "FirstExtraction"

# SphericalTensor packing: [T0, T1[0], T1[1], T1[2], T2[0], T2[1], T2[2], T2[3], T2[4]]
T0_SLICE = slice(0, 1)     # 1 component
T1_SLICE = slice(1, 4)     # 3 components
T2_SLICE = slice(4, 9)     # 5 components

# Ring type indices (8 types, matching C++ RingTypeIndex enum)
RING_TYPES = ["PHE", "TYR", "TRP6", "TRP5", "TRP9", "HIE", "HID", "HIP"]

# Amino acid type codes (matching C++ AminoAcid enum, 0-19)
AA_NAMES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL"]


@dataclass
class Conformation:
    """All features for one conformation (WT or ALA)."""
    protein_id: str
    variant: str  # "wt" or "ala"
    path: Path

    # Identity
    pos: np.ndarray = None           # (N, 3)
    element: np.ndarray = None       # (N,) int
    residue_index: np.ndarray = None # (N,) int
    residue_type: np.ndarray = None  # (N,) int

    # Calculator shielding tensors as SphericalTensor [T0, T1[3], T2[5]]
    bs: np.ndarray = None            # (N, 9)
    hm: np.ndarray = None            # (N, 9)
    mc: np.ndarray = None            # (N, 9)
    ringchi: np.ndarray = None       # (N, 9)
    pq: np.ndarray = None            # (N, 9)
    disp: np.ndarray = None          # (N, 9)
    coulomb: np.ndarray = None       # (N, 9)
    hbond: np.ndarray = None         # (N, 9)
    orca: np.ndarray = None          # (N, 9) — DFT truth

    # Per-ring-type decompositions
    bs_type_T0: np.ndarray = None    # (N, 8)
    bs_type_T2: np.ndarray = None    # (N, 40) = 8 types × 5 T2 comps
    hm_type_T0: np.ndarray = None    # (N, 8)
    hm_type_T2: np.ndarray = None    # (N, 40)
    pq_type_T0: np.ndarray = None    # (N, 8)
    pq_type_T2: np.ndarray = None    # (N, 40)
    disp_type_T0: np.ndarray = None  # (N, 8)
    disp_type_T2: np.ndarray = None  # (N, 40)

    # McConnell bond-category T2
    mc_category_T2: np.ndarray = None  # (N, 25) = 5 categories × 5 T2

    # Coulomb decompositions
    coulomb_E: np.ndarray = None           # (N, 3)
    coulomb_efg_backbone: np.ndarray = None  # (N, 9)
    coulomb_efg_aromatic: np.ndarray = None  # (N, 9)

    # APBS solvated fields (may be absent)
    apbs_E: np.ndarray = None        # (N, 3) or None
    apbs_efg: np.ndarray = None      # (N, 9) or None

    # xTB charges (may be absent)
    xtb_charges: np.ndarray = None   # (N,) or None

    # Ring proximity counts
    ring_counts: np.ndarray = None   # (N, 4)

    @property
    def n_atoms(self):
        return self.pos.shape[0] if self.pos is not None else 0


def _load_optional(path: Path) -> Optional[np.ndarray]:
    if path.exists():
        return np.load(path)
    return None


def load_conformation(protein_id: str, variant: str,
                      features_dir: Path = FEATURES_DIR) -> Optional[Conformation]:
    """Load all features for one conformation."""
    d = features_dir / protein_id / variant
    if not (d / "pos.npy").exists():
        return None

    c = Conformation(protein_id=protein_id, variant=variant, path=d)

    # Required identity arrays
    c.pos = np.load(d / "pos.npy")
    c.element = np.load(d / "element.npy")
    c.residue_index = np.load(d / "residue_index.npy")
    c.residue_type = np.load(d / "residue_type.npy")

    # Calculator shielding (SphericalTensor format)
    c.bs = np.load(d / "bs_shielding.npy")
    c.hm = np.load(d / "hm_shielding.npy")
    c.mc = np.load(d / "mc_shielding.npy")
    c.ringchi = np.load(d / "ringchi_shielding.npy")
    c.pq = np.load(d / "pq_shielding.npy")
    c.disp = np.load(d / "disp_shielding.npy")
    c.coulomb = np.load(d / "coulomb_shielding.npy")
    c.hbond = np.load(d / "hbond_shielding.npy")
    c.orca = _load_optional(d / "orca_total.npy")

    # Per-ring-type decompositions
    c.bs_type_T0 = np.load(d / "bs_per_type_T0.npy")
    c.bs_type_T2 = np.load(d / "bs_per_type_T2.npy")
    c.hm_type_T0 = np.load(d / "hm_per_type_T0.npy")
    c.hm_type_T2 = np.load(d / "hm_per_type_T2.npy")
    c.pq_type_T0 = np.load(d / "pq_per_type_T0.npy")
    c.pq_type_T2 = np.load(d / "pq_per_type_T2.npy")
    c.disp_type_T0 = np.load(d / "disp_per_type_T0.npy")
    c.disp_type_T2 = np.load(d / "disp_per_type_T2.npy")

    # McConnell and Coulomb decompositions
    c.mc_category_T2 = np.load(d / "mc_category_T2.npy")
    c.coulomb_E = np.load(d / "coulomb_E.npy")
    c.coulomb_efg_backbone = np.load(d / "coulomb_efg_backbone.npy")
    c.coulomb_efg_aromatic = np.load(d / "coulomb_efg_aromatic.npy")

    # Ring proximity
    c.ring_counts = np.load(d / "bs_ring_counts.npy")

    # Optional arrays
    c.apbs_E = _load_optional(d / "apbs_E.npy")
    c.apbs_efg = _load_optional(d / "apbs_efg.npy")
    c.xtb_charges = _load_optional(d / "xtb_charges.npy")

    return c


def T2(arr: np.ndarray) -> np.ndarray:
    """Extract T2 components (last 5 columns) from SphericalTensor array."""
    return arr[:, T2_SLICE]

def T0(arr: np.ndarray) -> np.ndarray:
    """Extract T0 component (first column) from SphericalTensor array."""
    return arr[:, T0_SLICE]

def T1(arr: np.ndarray) -> np.ndarray:
    """Extract T1 components (columns 1-3) from SphericalTensor array."""
    return arr[:, T1_SLICE]


def classical_sum(c: Conformation) -> np.ndarray:
    """Sum of all 8 classical calculator shielding tensors. (N, 9)."""
    return c.bs + c.hm + c.mc + c.ringchi + c.pq + c.disp + c.coulomb + c.hbond


def residual(c: Conformation) -> np.ndarray:
    """DFT minus classical sum. (N, 9) SphericalTensor."""
    return c.orca - classical_sum(c)


def list_proteins(features_dir: Path = FEATURES_DIR):
    """List all protein IDs that have both wt and ala extracted."""
    proteins = []
    if not features_dir.exists():
        return proteins
    for d in sorted(features_dir.iterdir()):
        if not d.is_dir():
            continue
        if (d / "wt" / "pos.npy").exists() and (d / "ala" / "pos.npy").exists():
            proteins.append(d.name)
    return proteins


def load_pair(protein_id: str, features_dir: Path = FEATURES_DIR):
    """Load both WT and ALA conformations for a protein."""
    wt = load_conformation(protein_id, "wt", features_dir)
    ala = load_conformation(protein_id, "ala", features_dir)
    return wt, ala
