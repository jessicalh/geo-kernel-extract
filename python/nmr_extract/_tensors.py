"""Tensor wrappers with e3nn Irreps.

Every wrapper holds a numpy array and exposes:
  .data    → np.ndarray (always available)
  .torch() → torch.Tensor
  .irreps  → e3nn.o3.Irreps instance

SphericalTensor packing (9 components):
    [0]     T0   isotropic scalar       (L=0, even parity)
    [1:4]   T1   antisymmetric vector   (L=1, odd parity)
    [4:9]   T2   traceless symmetric    (L=2, even parity, m=-2..+2)

T2 component ordering matches e3nn: m = -2, -1, 0, +1, +2.
Normalization is isometric (Frobenius-preserving), not orthonormal-on-sphere.
"""

from __future__ import annotations

import numpy as np
import torch
from e3nn.o3 import Irreps

from ._types import N_RING_TYPES, N_BOND_CATEGORIES, RingType, BondCategory


class SphericalTensor:
    """Irreducible spherical tensor: T0 (L=0) + T1 (L=1) + T2 (L=2).

    Shape ``(*, 9)`` packed as ``[T0, T1[3], T2[5]]``.

    T2 component ordering: m = -2, -1, 0, +1, +2 (matches e3nn).
    Normalization: isometric (Frobenius-preserving).

    Args:
        data: numpy array with last dimension 9.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("1x0e + 1x1o + 1x2e")

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 9:
            raise ValueError(
                f"{type(self).__name__}: last dim must be 9, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        """Raw numpy array ``(*, 9)``."""
        return self._data

    @property
    def irreps(self) -> Irreps:
        """``Irreps("1x0e+1x1o+1x2e")``."""
        return self.IRREPS

    def torch(self) -> torch.Tensor:
        """Zero-copy conversion to ``torch.Tensor``."""
        return torch.from_numpy(self._data)

    @property
    def T0(self) -> np.ndarray:
        """Isotropic component ``(*, 1)``.  L=0, even parity."""
        return self._data[..., 0:1]

    @property
    def T1(self) -> np.ndarray:
        """Antisymmetric pseudovector ``(*, 3)``.  L=1, odd parity."""
        return self._data[..., 1:4]

    @property
    def T2(self) -> np.ndarray:
        """Traceless symmetric tensor ``(*, 5)``.  L=2, even parity, m=-2..+2."""
        return self._data[..., 4:9]

    @property
    def isotropic(self) -> np.ndarray:
        """Scalar isotropic value ``(*,)`` (T0 squeezed)."""
        return self._data[..., 0]

    @property
    def T2_magnitude(self) -> np.ndarray:
        """L2 norm of T2 components ``(*,)``."""
        return np.linalg.norm(self._data[..., 4:9], axis=-1)

    @property
    def n_atoms(self) -> int:
        """Number of atoms (first dimension)."""
        return self._data.shape[0]

    def __repr__(self) -> str:
        return f"{type(self).__name__}(shape={self._data.shape}, irreps='{self.IRREPS}')"


class ShieldingTensor(SphericalTensor):
    """Shielding contribution in ppm.  Same structure as :class:`SphericalTensor`."""
    pass


class EFGTensor(SphericalTensor):
    """Electric field gradient in V/A^2.  Same structure as :class:`SphericalTensor`."""
    pass


class VectorField:
    """3D vector field (positions, E-fields, B-fields).

    Shape ``(*, 3)``.  Cartesian (x, y, z) components.

    Args:
        data: numpy array with last dimension 3.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("1x1o")

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(f"VectorField: last dim must be 3, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def irreps(self) -> Irreps:
        return self.IRREPS

    def torch(self) -> torch.Tensor:
        return torch.from_numpy(self._data)

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
    def n_atoms(self) -> int:
        return self._data.shape[0]

    def __repr__(self) -> str:
        return f"VectorField(shape={self._data.shape}, irreps='{self.IRREPS}')"


# ── Per-ring-type decompositions ────────────────────────────────────


class PerRingTypeT0:
    """Isotropic (T0) contribution decomposed by ring type.

    Shape ``(*, 8)`` -- one scalar per :class:`RingType`.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("8x0e")

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_RING_TYPES:
            raise ValueError(
                f"PerRingTypeT0: last dim must be {N_RING_TYPES}, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def irreps(self) -> Irreps:
        return self.IRREPS

    def for_type(self, rt: RingType) -> np.ndarray:
        return self._data[..., int(rt)]

    @property
    def total(self) -> np.ndarray:
        return self._data.sum(axis=-1)

    def __repr__(self) -> str:
        return f"PerRingTypeT0(shape={self._data.shape})"


class PerRingTypeT2:
    """T2 (L=2) contribution decomposed by ring type.

    Shape ``(*, 40)`` = 8 ring types x 5 T2 components.
    Use :meth:`as_block` for ``(*, 8, 5)`` view.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("8x2e")

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_RING_TYPES * 5:
            raise ValueError(
                f"PerRingTypeT2: last dim must be {N_RING_TYPES * 5}, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def irreps(self) -> Irreps:
        return self.IRREPS

    def for_type(self, rt: RingType) -> np.ndarray:
        i = int(rt)
        return self._data[..., i * 5:(i + 1) * 5]

    @property
    def total(self) -> np.ndarray:
        return self.as_block().sum(axis=-2)

    def as_block(self) -> np.ndarray:
        """Reshape to ``(N, 8, 5)``."""
        return self._data.reshape(*self._data.shape[:-1], N_RING_TYPES, 5)

    def __repr__(self) -> str:
        return f"PerRingTypeT2(shape={self._data.shape})"


# ── Per-bond-category decomposition ─────────────────────────────────


class PerBondCategoryT2:
    """T2 (L=2) McConnell contribution decomposed by bond category.

    Shape ``(*, 25)`` = 5 bond categories x 5 T2 components.
    Use :meth:`as_block` for ``(*, 5, 5)`` view.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("5x2e")

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_BOND_CATEGORIES * 5:
            raise ValueError(
                f"PerBondCategoryT2: last dim must be {N_BOND_CATEGORIES * 5}, "
                f"got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def irreps(self) -> Irreps:
        return self.IRREPS

    def for_category(self, cat: BondCategory) -> np.ndarray:
        i = int(cat)
        return self._data[..., i * 5:(i + 1) * 5]

    @property
    def total(self) -> np.ndarray:
        return self.as_block().sum(axis=-2)

    def as_block(self) -> np.ndarray:
        return self._data.reshape(*self._data.shape[:-1], N_BOND_CATEGORIES, 5)

    def __repr__(self) -> str:
        return f"PerBondCategoryT2(shape={self._data.shape})"


# ── Scalar feature types ────────────────────────────────────────────


class RingCounts:
    """(*, 4) ring proximity counts at 3/5/8/12 A."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 4:
            raise ValueError(f"RingCounts: last dim must be 4, got {data.shape}")
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

    def __repr__(self) -> str:
        return f"RingCounts(shape={self._data.shape})"


class McConnellScalars:
    """(*, 6) McConnell summary: CO/CN/sidechain/aromatic sums, nearest dists."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 6:
            raise ValueError(f"McConnellScalars: last dim must be 6, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def co_sum(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def cn_sum(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def sidechain_sum(self) -> np.ndarray:
        return self._data[..., 2]

    @property
    def aromatic_sum(self) -> np.ndarray:
        return self._data[..., 3]

    @property
    def nearest_CO_dist(self) -> np.ndarray:
        return self._data[..., 4]

    @property
    def nearest_CN_dist(self) -> np.ndarray:
        return self._data[..., 5]

    def __repr__(self) -> str:
        return f"McConnellScalars(shape={self._data.shape})"


class CoulombScalars:
    """(*, 4) Coulomb E-field summary.  Units: V/A."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 4:
            raise ValueError(f"CoulombScalars: last dim must be 4, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def E_magnitude(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def E_bond_proj(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def E_backbone_frac(self) -> np.ndarray:
        return self._data[..., 2]

    @property
    def aromatic_E_magnitude(self) -> np.ndarray:
        return self._data[..., 3]

    def __repr__(self) -> str:
        return f"CoulombScalars(shape={self._data.shape})"


class HBondScalars:
    """(*, 3) H-bond summary: nearest_dist, 1/r^3, count_within_3.5A."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(f"HBondScalars: last dim must be 3, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def nearest_dist(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def inv_d3(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def count_within_3_5A(self) -> np.ndarray:
        return self._data[..., 2]

    def __repr__(self) -> str:
        return f"HBondScalars(shape={self._data.shape})"


class DsspScalars:
    """(*, 5) DSSP: phi, psi, sasa, ss_helix, ss_sheet."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 5:
            raise ValueError(f"DsspScalars: last dim must be 5, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def phi(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def psi(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def sasa(self) -> np.ndarray:
        return self._data[..., 2]

    @property
    def ss_helix(self) -> np.ndarray:
        return self._data[..., 3]

    @property
    def ss_sheet(self) -> np.ndarray:
        return self._data[..., 4]

    def __repr__(self) -> str:
        return f"DsspScalars(shape={self._data.shape})"


class MopacScalars:
    """(*, 3) MOPAC per-atom: charge, s_pop, p_pop."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(f"MopacScalars: last dim must be 3, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def charge(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def s_pop(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def p_pop(self) -> np.ndarray:
        return self._data[..., 2]

    def __repr__(self) -> str:
        return f"MopacScalars(shape={self._data.shape})"


class MopacGlobal:
    """(4,) graph-level MOPAC scalars."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape != (4,):
            raise ValueError(f"MopacGlobal: expected shape (4,), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def heat_of_formation(self) -> float:
        return float(self._data[0])

    def __repr__(self) -> str:
        return f"MopacGlobal(hof={self.heat_of_formation:.1f} kcal/mol)"


class BondOrders:
    """(B, 3) sparse: [atom_i, atom_j, wiberg_order]."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(f"BondOrders: expected (B, 3), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def atom_i(self) -> np.ndarray:
        return self._data[:, 0].astype(np.intp)

    @property
    def atom_j(self) -> np.ndarray:
        return self._data[:, 1].astype(np.intp)

    @property
    def order(self) -> np.ndarray:
        return self._data[:, 2]

    @property
    def n_bonds(self) -> int:
        return self._data.shape[0]

    def to_dense(self, n_atoms: int) -> np.ndarray:
        mat = np.zeros((n_atoms, n_atoms), dtype=np.float64)
        i, j = self.atom_i, self.atom_j
        mat[i, j] = self._data[:, 2]
        mat[j, i] = self._data[:, 2]
        return mat

    def __repr__(self) -> str:
        return f"BondOrders(n_bonds={self.n_bonds})"


class DeltaScalars:
    """(*, 6) mutation delta metadata."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 6:
            raise ValueError(f"DeltaScalars: last dim must be 6, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def matched_mask(self) -> np.ndarray:
        return self._data[..., 0] > 0.5

    @property
    def delta_T0(self) -> np.ndarray:
        return self._data[..., 1]

    @property
    def nearest_removed_ring_dist(self) -> np.ndarray:
        return self._data[..., 2]

    @property
    def delta_partial_charge(self) -> np.ndarray:
        """ff14SB partial charge delta (WT - ALA)."""
        return self._data[..., 3]

    @property
    def delta_mopac_charge(self) -> np.ndarray:
        """MOPAC Mulliken charge delta (0 if no MOPAC)."""
        return self._data[..., 4]

    @property
    def match_distance(self) -> np.ndarray:
        """Spatial distance between matched WT/ALA atoms (A)."""
        return self._data[..., 5]

    def __repr__(self) -> str:
        return f"DeltaScalars(shape={self._data.shape})"


class DeltaAPBS:
    """(*, 12) APBS delta: delta_E(3) + delta_EFG(9)."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 12:
            raise ValueError(f"DeltaAPBS: last dim must be 12, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def delta_E(self) -> VectorField:
        return VectorField(self._data[..., :3])

    @property
    def delta_efg(self) -> EFGTensor:
        return EFGTensor(self._data[..., 3:])

    def __repr__(self) -> str:
        return f"DeltaAPBS(shape={self._data.shape})"


class DeltaRingProximity:
    """(*, R*6) removed ring geometry.  Per ring: [dist, z, rho, theta, mcconnell, exp_decay]."""

    __slots__ = ("_data", "_n_rings")
    COLS_PER_RING = 6

    def __init__(self, data: np.ndarray):
        if data.shape[-1] % self.COLS_PER_RING != 0:
            raise ValueError(
                f"DeltaRingProximity: last dim must be multiple of "
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
        return self._data[..., i * 6:(i + 1) * 6]

    def distance(self, i: int) -> np.ndarray:
        return self._data[..., i * 6]

    def __repr__(self) -> str:
        return f"DeltaRingProximity(shape={self._data.shape}, n_rings={self._n_rings})"
