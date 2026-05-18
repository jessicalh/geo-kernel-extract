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
    """Shielding contribution tensor.  Same structure as :class:`SphericalTensor`.

    Per-instance units are declared by the corresponding catalog entry
    (``ArraySpec.units``). DFT-derived shielding (orca_*, tripeptide_*,
    larsen_hbond_*_shielding) is in ppm. Classical-kernel-derived
    shielding (bs_*, hm_*, mc_*, pq_*, disp_*, hbond_*, ringchi_*,
    coulomb_shielding) is in the kernel's native unit (ppm·T/nA, Å⁻¹,
    Å⁻³, Å⁻⁵, Å⁻⁶, V/Å²) — calibration multiplies by the relevant
    parameter to map to ppm. See OBJECT_MODEL.md drift-table section.
    """
    pass


class EFGTensor:
    """Electric field gradient in V/A^2 — symmetric-traceless T2 only.

    Shape ``(*, 5)``. All EFGs in this codebase are constructed from
    symmetric outer-product physics — q·(3r⊗r/r⁵ − I/r³) for the
    Coulomb-family EFGs (water, Coulomb, MOPAC Coulomb, AIMNet2) and
    the Hessian of φ for APBS. After the explicit traceless projection
    each calculator performs, both T0 (trace) and T1 (antisymmetric
    pseudovector) are structural zeros. Only the symmetric-traceless
    T2 (5 real-spherical-tesseral components m=-2..+2) carries signal.

    Re-typed from a 9-component SphericalTensor subclass to a standalone
    5-component class on 2026-05-18 (codex review R2 M1 expansion). The
    old shape emitted 4 always-zero channels per atom; the new shape
    saves storage AND signals the physics correctly to e3nn / MACE /
    NequIP downstream consumers via the ``1x2e`` Irreps declaration.

    Args:
        data: numpy array with last dimension 5.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("1x2e")

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 5:
            raise ValueError(
                f"{type(self).__name__}: last dim must be 5 (T2-only), "
                f"got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        """Raw numpy array ``(*, 5)``. Components m=-2,-1,0,+1,+2."""
        return self._data

    @property
    def irreps(self) -> Irreps:
        """``Irreps("1x2e")`` — symmetric-traceless rank-2, parity-even."""
        return self.IRREPS

    def torch(self) -> torch.Tensor:
        """Zero-copy conversion to ``torch.Tensor``."""
        return torch.from_numpy(self._data)

    @property
    def T2(self) -> np.ndarray:
        """The 5-component T2 array ``(*, 5)``. Equivalent to ``.data``."""
        return self._data

    @property
    def T2_magnitude(self) -> np.ndarray:
        """L2 norm of T2 components ``(*,)``. Rotationally invariant."""
        return np.linalg.norm(self._data, axis=-1)

    @property
    def n_atoms(self) -> int:
        return self._data.shape[0]

    def __repr__(self) -> str:
        return f"{type(self).__name__}(shape={self._data.shape}, irreps='{self.IRREPS}')"


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
    """``(*, 4)`` MOPAC per-atom: charge, s_pop, p_pop, valency."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 4:
            raise ValueError(f"MopacScalars: last dim must be 4, got {data.shape}")
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

    @property
    def valency(self) -> np.ndarray:
        """Sum of Wiberg bond orders (QM coordination number)."""
        return self._data[..., 3]

    def __repr__(self) -> str:
        return f"MopacScalars(shape={self._data.shape})"


class MopacGlobal:
    """(4,) graph-level MOPAC scalars: [hof, dipole_x, dipole_y, dipole_z]."""

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

    @property
    def dipole(self) -> np.ndarray:
        """Molecular dipole vector (Debye), shape (3,)."""
        return self._data[1:4]

    @property
    def dipole_magnitude(self) -> float:
        return float(np.linalg.norm(self._data[1:4]))

    def __repr__(self) -> str:
        return (f"MopacGlobal(hof={self.heat_of_formation:.1f} kcal/mol, "
                f"|dipole|={self.dipole_magnitude:.2f} D)")


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
    """(*, 12) APBS delta: delta_E(3) + delta_EFG_full_sphericaltensor(9).

    The 9-component delta_EFG slice is the full T0+T1+T2 SphericalTensor
    packing produced by MutationDeltaResult — T0 and T1 are structurally
    zero (the difference of two symmetric-traceless matrices is itself
    symmetric-traceless), but the emitted shape stays 9 for backward
    compatibility with the mutation-delta consumer. The `.delta_efg`
    accessor returns a SphericalTensor (not EFGTensor) because the new
    EFGTensor class is T2-only 5-component. Downstream that wants just
    T2 should use ``apbs.delta_efg.T2`` to get the 5-component view.
    """

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
    def delta_efg(self) -> "SphericalTensor":
        """Returns the full 9-component packed SphericalTensor. T0 and T1
        slots are structurally zero but emitted for the mutation-delta
        path's backward-compat schema. Use `.T2` for the 5 nonzero
        components."""
        return SphericalTensor(self._data[..., 3:])

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


# ── AIMNet2 types ──────────────────────────────────────────────────────


class AIMNet2Charges:
    """(N,) per-atom Hirshfeld charges from AIMNet2 wB97M model."""

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def charges(self) -> np.ndarray:
        return self._data

    def __repr__(self) -> str:
        return f"AIMNet2Charges(n={len(self._data)})"


class AIMNet2AimEmbedding:
    """(N, 256) learned electronic structure embedding per atom.

    Geometry-dependent: changes per frame.  Encodes hybridisation,
    polarisability, conjugation, charge transfer.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim == 2 and data.shape[-1] != 256:
            raise ValueError(f"AIMNet2AimEmbedding: expected 256 dims, got {data.shape[-1]}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def __repr__(self) -> str:
        return f"AIMNet2AimEmbedding(shape={self._data.shape})"


class AIMNet2Polarisability:
    """(N, 3) per-atom charge-polarisation gradient via autograd.

    The vector field is dL/d(r_i) where L = sum_j q_j^2 over
    non-sentinel atoms (the L2-of-charges objective; sum(q) is
    constant under AIMNet2's charge-conservation projection so its
    gradient is ~0). Charge-weighted per-atom polarisability gradient.

    Companion scalar is the L2 norm of the vector, stored separately
    in `aimnet2_polarisability_scalar.npy`.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("1x1o")  # vector under SO(3); odd parity

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(
                f"AIMNet2Polarisability: expected last dim 3, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def irreps(self) -> Irreps:
        return self.IRREPS

    @property
    def vectors(self) -> np.ndarray:
        """Per-atom polarisability gradient vectors ``(N, 3)``."""
        return self._data

    @property
    def norms(self) -> np.ndarray:
        """Per-atom L2 norm ``(N,)``. Equal to ``aimnet2_polarisability_scalar.npy``
        up to floating-point precision; this property recomputes from the
        vectors so it's always consistent with the loaded vector field."""
        return np.linalg.norm(self._data, axis=-1)

    def __repr__(self) -> str:
        return f"AIMNet2Polarisability(shape={self._data.shape})"

