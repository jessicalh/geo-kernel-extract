"""
Typed wrappers for NMR shielding feature arrays.

Every NPY array exported by the C++ extractor gets a Python type that
enforces shape on construction and exposes named properties instead of
raw index arithmetic.  Type-distinct even when shapes coincide:
ShieldingTensor and EFGTensor are both (*, 9) but incompatible types.

SphericalTensor packing (9 components):
    [0]     T0     isotropic scalar
    [1:4]   T1     antisymmetric pseudovector  (x, y, z)
    [4:9]   T2     traceless symmetric tensor  (m = -2, -1, 0, +1, +2)
"""

from __future__ import annotations

from enum import IntEnum

import numpy as np


# ── Enums ────────────────────────────────────────────────────────────

class RingType(IntEnum):
    PHE          = 0   # Phe benzene
    TYR          = 1   # Tyr phenol
    TRP_benzene  = 2   # Trp 6-ring
    TRP_pyrrole  = 3   # Trp 5-ring
    TRP_perimeter = 4  # Trp indole 9-atom perimeter
    HIS          = 5   # His imidazole (neutral, N-delta protonated)
    HID          = 6   # Hid imidazole (doubly protonated)
    HIE          = 7   # Hie imidazole (neutral, N-epsilon protonated)

N_RING_TYPES = 8

class BondCategory(IntEnum):
    backbone_total  = 0   # BackboneOther + PeptideCO + PeptideCN
    sidechain_total = 1   # SidechainCO + SidechainOther
    aromatic_total  = 2   # Aromatic bonds
    CO_nearest      = 3   # Nearest PeptideCO bond
    CN_nearest      = 4   # Nearest PeptideCN bond

N_BOND_CATEGORIES = 5


# ── Base tensor types ────────────────────────────────────────────────

class SphericalTensor:
    """(*, 9) array: T0 + T1 + T2 in standard packing."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 9:
            raise ValueError(
                f"SphericalTensor needs last dim 9, got shape {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def T0(self) -> np.ndarray:
        """Isotropic component (*, 1)."""
        return self._data[..., 0:1]

    @property
    def T1(self) -> np.ndarray:
        """Antisymmetric pseudovector (*, 3)."""
        return self._data[..., 1:4]

    @property
    def T2(self) -> np.ndarray:
        """Traceless symmetric tensor (*, 5), m = -2..+2."""
        return self._data[..., 4:9]

    @property
    def isotropic(self) -> np.ndarray:
        """Scalar isotropic value (*,)."""
        return self._data[..., 0]

    @property
    def T2_magnitude(self) -> np.ndarray:
        """L2 norm of T2 components (*,)."""
        return np.linalg.norm(self._data[..., 4:9], axis=-1)

    @property
    def n_atoms(self) -> int:
        return self._data.shape[0]

    def __repr__(self) -> str:
        return f"{type(self).__name__}(shape={self._data.shape})"


class VectorField:
    """(*, 3) vector field — electric field, magnetic field, or position."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(
                f"VectorField needs last dim 3, got shape {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def x(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def y(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def z(self) -> np.ndarray:
        return self._data[..., 2]

    @property
    def magnitude(self) -> np.ndarray:
        return np.linalg.norm(self._data, axis=-1)

    @property
    def direction(self) -> np.ndarray:
        """Unit vector (*,3). Zero where magnitude is zero."""
        mag = self.magnitude[..., np.newaxis]
        return np.divide(self._data, mag, where=mag > 1e-12,
                         out=np.zeros_like(self._data))

    @property
    def n_atoms(self) -> int:
        return self._data.shape[0]

    def __repr__(self) -> str:
        return f"VectorField(shape={self._data.shape})"


# ── Typed tensor subtypes ────────────────────────────────────────────

class ShieldingTensor(SphericalTensor):
    """Shielding contribution from a calculator.  Units: ppm."""
    pass


class EFGTensor(SphericalTensor):
    """Electric field gradient in spherical tensor form.  Units: V/A^2."""
    pass


# ── Per-ring-type decompositions ─────────────────────────────────────

class PerRingTypeT0:
    """(*, 8) isotropic shielding decomposed by ring type."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_RING_TYPES:
            raise ValueError(
                f"PerRingTypeT0 needs last dim {N_RING_TYPES}, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def for_ring_type(self, rt: RingType | str) -> np.ndarray:
        """T0 contribution from one ring type (*,)."""
        idx = RingType[rt] if isinstance(rt, str) else int(rt)
        return self._data[..., idx]

    @property
    def total(self) -> np.ndarray:
        """Sum over all ring types (*,)."""
        return self._data.sum(axis=-1)

    def __repr__(self) -> str:
        return f"PerRingTypeT0(shape={self._data.shape})"


class PerRingTypeT2:
    """(*, 40) T2 shielding decomposed by ring type.  8 types x 5 T2 components."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_RING_TYPES * 5:
            raise ValueError(
                f"PerRingTypeT2 needs last dim {N_RING_TYPES * 5}, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def for_ring_type(self, rt: RingType | str) -> np.ndarray:
        """T2 components for one ring type (*, 5)."""
        idx = RingType[rt] if isinstance(rt, str) else int(rt)
        return self._data[..., idx * 5:(idx + 1) * 5]

    @property
    def total(self) -> np.ndarray:
        """Sum of T2 across all ring types (*, 5)."""
        reshaped = self._data.reshape(*self._data.shape[:-1], N_RING_TYPES, 5)
        return reshaped.sum(axis=-2)

    def magnitude_per_type(self) -> np.ndarray:
        """L2 norm of T2 per ring type (*, 8)."""
        reshaped = self._data.reshape(*self._data.shape[:-1], N_RING_TYPES, 5)
        return np.linalg.norm(reshaped, axis=-1)

    def __repr__(self) -> str:
        return f"PerRingTypeT2(shape={self._data.shape})"


# ── Per-bond-category decomposition ──────────────────────────────────

class PerBondCategoryT2:
    """(*, 25) McConnell T2 by bond category.  5 categories x 5 T2 components."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_BOND_CATEGORIES * 5:
            raise ValueError(
                f"PerBondCategoryT2 needs last dim {N_BOND_CATEGORIES * 5}, "
                f"got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def for_category(self, cat: BondCategory | str) -> np.ndarray:
        """T2 components for one bond category (*, 5)."""
        idx = BondCategory[cat] if isinstance(cat, str) else int(cat)
        return self._data[..., idx * 5:(idx + 1) * 5]

    @property
    def total(self) -> np.ndarray:
        """Sum of T2 across all categories (*, 5)."""
        reshaped = self._data.reshape(*self._data.shape[:-1], N_BOND_CATEGORIES, 5)
        return reshaped.sum(axis=-2)

    def magnitude_per_category(self) -> np.ndarray:
        """L2 norm of T2 per category (*, 5)."""
        reshaped = self._data.reshape(*self._data.shape[:-1], N_BOND_CATEGORIES, 5)
        return np.linalg.norm(reshaped, axis=-1)

    def __repr__(self) -> str:
        return f"PerBondCategoryT2(shape={self._data.shape})"


# ── Scalar feature types ─────────────────────────────────────────────

class McConnellScalars:
    """(*, 6) McConnell summary scalars."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 6:
            raise ValueError(f"McConnellScalars needs last dim 6, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def co_sum(self) -> np.ndarray:
        """Sum of McConnell f for PeptideCO bonds."""
        return self._data[..., 0]

    @property
    def cn_sum(self) -> np.ndarray:
        """Sum of McConnell f for PeptideCN bonds."""
        return self._data[..., 1]

    @property
    def sidechain_sum(self) -> np.ndarray:
        """Sum of McConnell f for SidechainCO bonds."""
        return self._data[..., 2]

    @property
    def aromatic_sum(self) -> np.ndarray:
        """Sum of McConnell f for Aromatic bonds."""
        return self._data[..., 3]

    @property
    def nearest_CO_dist(self) -> np.ndarray:
        """Distance to nearest PeptideCO bond (A).  Sentinel if absent."""
        return self._data[..., 4]

    @property
    def nearest_CN_dist(self) -> np.ndarray:
        """Distance to nearest PeptideCN bond (A).  Sentinel if absent."""
        return self._data[..., 5]


class CoulombScalars:
    """(*, 4) Coulomb E-field summary scalars.  Units: V/A."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 4:
            raise ValueError(f"CoulombScalars needs last dim 4, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def E_magnitude(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def E_bond_proj(self) -> np.ndarray:
        """E-field projected onto primary bond axis."""
        return self._data[..., 1]

    @property
    def E_backbone_frac(self) -> np.ndarray:
        """Backbone E component along total E direction."""
        return self._data[..., 2]

    @property
    def aromatic_E_magnitude(self) -> np.ndarray:
        return self._data[..., 3]


class HBondScalars:
    """(*, 3) hydrogen bond summary scalars."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(f"HBondScalars needs last dim 3, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def nearest_dist(self) -> np.ndarray:
        """Distance to nearest H-bond (A).  Sentinel if none."""
        return self._data[..., 0]

    @property
    def inv_d3(self) -> np.ndarray:
        """1/r^3 of nearest H-bond (A^-3)."""
        return self._data[..., 1]

    @property
    def count_within_3_5A(self) -> np.ndarray:
        """Number of H-bonds within 3.5 A."""
        return self._data[..., 2]


class MopacAtomScalars:
    """(*, 3) MOPAC per-atom electronic structure."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(f"MopacAtomScalars needs last dim 3, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def charge(self) -> np.ndarray:
        """Mulliken charge from PM7."""
        return self._data[..., 0]

    @property
    def s_pop(self) -> np.ndarray:
        """s orbital population."""
        return self._data[..., 1]

    @property
    def p_pop(self) -> np.ndarray:
        """p orbital population (0 for H atoms)."""
        return self._data[..., 2]


class RingCounts:
    """(*, 4) ring proximity counts at distance thresholds."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 4:
            raise ValueError(f"RingCounts needs last dim 4, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def within_3A(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def within_5A(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def within_8A(self) -> np.ndarray:
        return self._data[..., 2]

    @property
    def within_12A(self) -> np.ndarray:
        return self._data[..., 3]


class DsspBackbone:
    """(*, 5) DSSP secondary structure and backbone geometry."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 5:
            raise ValueError(f"DsspBackbone needs last dim 5, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def phi(self) -> np.ndarray:
        """Backbone phi angle (radians)."""
        return self._data[..., 0]

    @property
    def psi(self) -> np.ndarray:
        """Backbone psi angle (radians)."""
        return self._data[..., 1]

    @property
    def sasa(self) -> np.ndarray:
        """Solvent accessible surface area (A^2)."""
        return self._data[..., 2]

    @property
    def ss_helix(self) -> np.ndarray:
        """1 if helix (H/G/I), 0 otherwise."""
        return self._data[..., 3]

    @property
    def ss_sheet(self) -> np.ndarray:
        """1 if sheet (E/B), 0 otherwise."""
        return self._data[..., 4]


# ── Delta-specific types ─────────────────────────────────────────────

class DeltaScalars:
    """(*, 6) mutation delta metadata per atom."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 6:
            raise ValueError(f"DeltaScalars needs last dim 6, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def matched(self) -> np.ndarray:
        """1 if atom was matched WT->ALA, 0 otherwise."""
        return self._data[..., 0]

    @property
    def matched_mask(self) -> np.ndarray:
        """Boolean mask of matched atoms."""
        return self._data[..., 0] > 0.5

    @property
    def delta_T0(self) -> np.ndarray:
        """Delta isotropic shielding (ppm)."""
        return self._data[..., 1]

    @property
    def nearest_removed_ring_dist(self) -> np.ndarray:
        """Distance to nearest removed ring (A).  99.0 if none."""
        return self._data[..., 2]

    @property
    def delta_partial_charge(self) -> np.ndarray:
        """ff14SB partial charge delta."""
        return self._data[..., 3]

    @property
    def delta_mopac_charge(self) -> np.ndarray:
        """MOPAC Mulliken charge delta (0 if no MOPAC)."""
        return self._data[..., 4]

    @property
    def match_distance(self) -> np.ndarray:
        """Spatial distance between matched WT/ALA atoms (A)."""
        return self._data[..., 5]


class DeltaAPBS:
    """(*, 12) APBS electrostatics delta: delta_E(3) + delta_EFG(9)."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 12:
            raise ValueError(f"DeltaAPBS needs last dim 12, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def delta_E(self) -> VectorField:
        """WT-ALA electric field delta (V/A)."""
        return VectorField(self._data[..., :3])

    @property
    def delta_efg(self) -> EFGTensor:
        """WT-ALA EFG delta as SphericalTensor."""
        return EFGTensor(self._data[..., 3:])


class DeltaRingProximity:
    """(*, R*6) geometry of removed aromatic rings.

    R = number of removed rings (varies per protein).
    Per ring: [distance, z, rho, theta, mcconnell_factor, exp_decay].
    """
    __slots__ = ('_data', '_n_rings')
    COLS_PER_RING = 6

    def __init__(self, data: np.ndarray):
        if data.shape[-1] % self.COLS_PER_RING != 0:
            raise ValueError(
                f"DeltaRingProximity last dim must be multiple of "
                f"{self.COLS_PER_RING}, got {data.shape}")
        self._data = data
        self._n_rings = data.shape[-1] // self.COLS_PER_RING

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def n_removed_rings(self) -> int:
        return self._n_rings

    def ring(self, i: int) -> np.ndarray:
        """All 6 features for removed ring i: (*, 6)."""
        return self._data[..., i * 6:(i + 1) * 6]

    def distance(self, i: int) -> np.ndarray:
        """Distance to removed ring i (A)."""
        return self._data[..., i * 6]

    def z(self, i: int) -> np.ndarray:
        """Cylindrical z in ring frame (A)."""
        return self._data[..., i * 6 + 1]

    def rho(self, i: int) -> np.ndarray:
        """Cylindrical radial coordinate (A)."""
        return self._data[..., i * 6 + 2]

    def theta(self, i: int) -> np.ndarray:
        """Cylindrical angle (radians)."""
        return self._data[..., i * 6 + 3]

    def mcconnell_factor(self, i: int) -> np.ndarray:
        """(3cos^2(theta) - 1) / r^3 (A^-3)."""
        return self._data[..., i * 6 + 4]

    def exp_decay(self, i: int) -> np.ndarray:
        """Exponential decay factor with tau=4A."""
        return self._data[..., i * 6 + 5]

    def __repr__(self) -> str:
        return (f"DeltaRingProximity(shape={self._data.shape}, "
                f"n_rings={self._n_rings})")


# ── Sparse types ─────────────────────────────────────────────────────

class BondOrders:
    """(M, 3) sparse bond order edge list: [atom_i, atom_j, wiberg_order]."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(f"BondOrders needs shape (M, 3), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def atom_i(self) -> np.ndarray:
        """First atom index (int)."""
        return self._data[:, 0].astype(np.int32)

    @property
    def atom_j(self) -> np.ndarray:
        """Second atom index, j > i (int)."""
        return self._data[:, 1].astype(np.int32)

    @property
    def order(self) -> np.ndarray:
        """Wiberg bond order (float)."""
        return self._data[:, 2]

    @property
    def n_bonds(self) -> int:
        return self._data.shape[0]

    def for_atom(self, idx: int) -> np.ndarray:
        """Bond orders involving atom idx: (K, 2) of [other_atom, order]."""
        mask_i = self._data[:, 0].astype(int) == idx
        mask_j = self._data[:, 1].astype(int) == idx
        rows_i = np.column_stack([self._data[mask_i, 1], self._data[mask_i, 2]])
        rows_j = np.column_stack([self._data[mask_j, 0], self._data[mask_j, 2]])
        if rows_i.size == 0 and rows_j.size == 0:
            return np.empty((0, 2), dtype=np.float64)
        parts = [p for p in (rows_i, rows_j) if p.size > 0]
        return np.vstack(parts)

    def to_dense(self, n_atoms: int) -> np.ndarray:
        """Symmetric (N, N) bond order matrix."""
        mat = np.zeros((n_atoms, n_atoms), dtype=np.float64)
        i = self._data[:, 0].astype(int)
        j = self._data[:, 1].astype(int)
        mat[i, j] = self._data[:, 2]
        mat[j, i] = self._data[:, 2]
        return mat

    def __repr__(self) -> str:
        return f"BondOrders(n_bonds={self.n_bonds})"


# ── Global types ─────────────────────────────────────────────────────

class MopacGlobal:
    """(4,) graph-level MOPAC scalars."""
    __slots__ = ('_data',)

    def __init__(self, data: np.ndarray):
        if data.shape != (4,):
            raise ValueError(f"MopacGlobal needs shape (4,), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def heat_of_formation(self) -> float:
        """kcal/mol from PM7."""
        return float(self._data[0])

    def __repr__(self) -> str:
        return f"MopacGlobal(hof={self.heat_of_formation:.1f} kcal/mol)"
