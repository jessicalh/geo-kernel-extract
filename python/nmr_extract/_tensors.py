"""Tensor wrappers for project-native NMR tensor payloads.

Every wrapper holds a numpy array and exposes:
  .data    → np.ndarray (always available)
  .torch() → torch.Tensor

Raw SphericalTensor/EFGTensor/T2 block payloads use the project basis defined
by the C++ SphericalTensor implementation. They are not e3nn tensors until an
explicit conversion method is called.

SphericalTensor packing (9 components):
    [0]     T0   isotropic scalar       (L=0, even parity)
    [1:4]   T1   Cartesian Levi-Civita dual x,y,z (L=1, even parity)
    [4:9]   T2   traceless symmetric    (L=2, even parity, m=-2..+2)

Canonical PackFull9 layout string:
    "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"

T2 component ordering is m = -2, -1, 0, +1, +2 in the project real-tesseral
isometric basis. Use project_t2_to_e3nn/project_full9_to_e3nn or wrapper
.to_e3nn() methods before passing data to e3nn Irreps consumers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from e3nn.o3 import Irreps

from ._types import N_RING_TYPES, N_BOND_CATEGORIES, RingType, BondCategory


PROJECT_FULL9_TENSOR_BASIS = "project_native_full9_spherical_tensor_v1"
PROJECT_T2_TENSOR_BASIS = "project_native_t2_isometric_real_tesseral_v1"
PROJECT_FULL9_COMPONENT_ORDER = (
    "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
)
PROJECT_T2_COMPONENT_ORDER = "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
PROJECT_FULL9_TENSOR_FRAME = "conformation_cartesian_xyz"
PROJECT_T2_TENSOR_FRAME = "conformation_cartesian_xyz"
PROJECT_E3NN_EXPORT_NOTE = (
    "raw project tensor; call to_e3nn()/to_e3nn_T2() or "
    "project_t2_to_e3nn() before using e3nn Irreps"
)
PACK_FULL9_IRREP_LAYOUT = PROJECT_FULL9_COMPONENT_ORDER

PROJECT_T2_TO_E3NN = np.array([
    [0.0, 0.0, 0.0, 1.0, 0.0],
    [1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, -0.5, 0.0, -np.sqrt(3.0) / 2.0],
    [0.0, 1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, np.sqrt(3.0) / 2.0, 0.0, -0.5],
], dtype=np.float64)


@dataclass(frozen=True)
class E3nnTensor:
    data: np.ndarray
    irreps: Irreps

    def torch(self) -> torch.Tensor:
        return torch.from_numpy(self.data)


def project_t2_to_e3nn(t2: np.ndarray) -> np.ndarray:
    """Convert project-native T2 coordinates to e3nn ``2e`` coordinates."""
    arr = np.asarray(t2)
    if arr.shape[-1] != 5:
        raise ValueError(f"project_t2_to_e3nn: last dim must be 5, got {arr.shape}")
    return arr @ PROJECT_T2_TO_E3NN.T


def project_full9_to_e3nn(full9: np.ndarray) -> np.ndarray:
    """Convert only the T2 block of PackFull9 project tensors to e3nn."""
    arr = np.asarray(full9)
    if arr.shape[-1] != 9:
        raise ValueError(f"project_full9_to_e3nn: last dim must be 9, got {arr.shape}")
    out = np.array(arr, dtype=np.result_type(arr.dtype, np.float64), copy=True)
    out[..., 4:9] = project_t2_to_e3nn(arr[..., 4:9])
    return out


class SphericalTensor:
    """Project-native spherical tensor: T0 + Cartesian T1 + real-tesseral T2.

    Shape ``(*, 9)`` packed as ``[T0, T1[3], T2[5]]``.

    T2 component ordering: m = -2, -1, 0, +1, +2 in the project basis.
    Normalization: isometric (Frobenius-preserving).

    Args:
        data: numpy array with last dimension 9.
    """

    __slots__ = ("_data",)
    TENSOR_BASIS = PROJECT_FULL9_TENSOR_BASIS
    COMPONENT_ORDER = PROJECT_FULL9_COMPONENT_ORDER
    TENSOR_FRAME = PROJECT_FULL9_TENSOR_FRAME
    E3NN_EXPORT = PROJECT_E3NN_EXPORT_NOTE

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
    def tensor_basis(self) -> str:
        return self.TENSOR_BASIS

    @property
    def component_order(self) -> str:
        return self.COMPONENT_ORDER

    @property
    def tensor_frame(self) -> str:
        return self.TENSOR_FRAME

    @property
    def e3nn_export(self) -> str:
        return self.E3NN_EXPORT

    def torch(self) -> torch.Tensor:
        """Zero-copy conversion to ``torch.Tensor``."""
        return torch.from_numpy(self._data)

    def to_e3nn(self) -> E3nnTensor:
        return E3nnTensor(
            project_full9_to_e3nn(self._data),
            Irreps("1x0e + 1x1e + 1x2e"),
        )

    def to_e3nn_T2(self) -> E3nnTensor:
        return E3nnTensor(project_t2_to_e3nn(self.T2), Irreps("1x2e"))

    @property
    def T0(self) -> np.ndarray:
        """Isotropic component ``(*, 1)``.  L=0, even parity."""
        return self._data[..., 0:1]

    @property
    def T1(self) -> np.ndarray:
        """Antisymmetric pseudovector ``(*, 3)``.  L=1, even parity."""
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
        return (
            f"{type(self).__name__}(shape={self._data.shape}, "
            f"basis='{self.TENSOR_BASIS}')"
        )


class ShieldingTensor(SphericalTensor):
    """Shielding contribution tensor.  Same structure as :class:`SphericalTensor`.

    Per-instance units are declared by the corresponding catalog entry
    (``ArraySpec.units``). DFT-derived shielding (orca_*,
    larsen_hbond_*_shielding) is in ppm. Classical-kernel-derived
    payloads use their catalog-declared native units (for example ppm·T/nA,
    Å⁻¹, Å⁻³, Å⁻⁵, or V/Å²). Some are geometry-only kernels with a deferred
    scale, not physical shielding in final units; consult
    ``ArraySpec.scaling_contract`` rather than inferring from this wrapper.
    """
    pass


class EFGTensor:
    """Electric field gradient in V/A^2 — symmetric-traceless T2 only.

    Shape ``(*, 5)``. All EFGs in this codebase are constructed from
    symmetric outer-product physics — q·(3r⊗r/r⁵ − I/r³) for the
    Coulomb-family EFGs (water, Coulomb, MOPAC Coulomb, AIMNet2) and
    the Hessian of φ for APBS. For finite analytic rows, T0 (trace) and T1
    (antisymmetric pseudovector) are structural zeros. The legacy force-field
    and MOPAC Coulomb paths apply an elementwise nonfinite replacement after
    traceless projection; an exceptional sanitized source matrix can regain a
    trace without a mask, although this wrapper serializes only its T2 block.
    See the per-array catalog validity contract.

    Re-typed from a 9-component SphericalTensor subclass to a standalone
    5-component class on 2026-05-18 (codex review R2 M1 expansion). The
    old shape emitted 4 always-zero channels per atom; the new shape
    saves storage and signals the physics correctly as a raw project
    T2 tensor. Call ``to_e3nn()`` before passing to e3nn / MACE /
    NequIP downstream consumers.

    Args:
        data: numpy array with last dimension 5.
    """

    __slots__ = ("_data",)
    TENSOR_BASIS = PROJECT_T2_TENSOR_BASIS
    COMPONENT_ORDER = PROJECT_T2_COMPONENT_ORDER
    TENSOR_FRAME = PROJECT_T2_TENSOR_FRAME
    E3NN_EXPORT = PROJECT_E3NN_EXPORT_NOTE

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 5:
            if data.shape[-1] == 9:
                raise ValueError(
                    f"{type(self).__name__}: last dim must be 5 (T2-only, "
                    f"post-2026-05-18 schema). Got {data.shape} — this "
                    f"looks like the pre-2026-05-18 9-component "
                    f"(T0+T1+T2) packing. Schema rev intentionally drops "
                    f"the structurally-zero T0 (trace) and T1 (antisymmetric "
                    f"pseudovector) channels for finite analytic EFG rows. To "
                    f"migrate older NPYs in place: `T2 = old[..., 4:9]` "
                    f"and load that. To re-extract: run with current build "
                    f"of WaterFieldResult / CoulombResult / MopacCoulombResult "
                    f"/ ApbsFieldResult / AIMNet2Result.")
            raise ValueError(
                f"{type(self).__name__}: last dim must be 5 (T2-only), "
                f"got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        """Raw numpy array ``(*, 5)``. Components m=-2,-1,0,+1,+2."""
        return self._data

    @property
    def tensor_basis(self) -> str:
        return self.TENSOR_BASIS

    @property
    def component_order(self) -> str:
        return self.COMPONENT_ORDER

    @property
    def tensor_frame(self) -> str:
        return self.TENSOR_FRAME

    @property
    def e3nn_export(self) -> str:
        return self.E3NN_EXPORT

    def torch(self) -> torch.Tensor:
        """Zero-copy conversion to ``torch.Tensor``."""
        return torch.from_numpy(self._data)

    def to_e3nn(self) -> E3nnTensor:
        return E3nnTensor(project_t2_to_e3nn(self._data), Irreps("1x2e"))

    def to_e3nn_T2(self) -> E3nnTensor:
        return self.to_e3nn()

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
        return (
            f"{type(self).__name__}(shape={self._data.shape}, "
            f"basis='{self.TENSOR_BASIS}')"
        )


class PositionField:
    """Cartesian positions with an affine SE(3) transformation law.

    Shape ``(*, 3)``.  Positions obey ``p' = R p + t`` and deliberately do
    not expose an e3nn irrep: treating them as a translation-free ``1o``
    vector field is a contract error.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 3:
            raise ValueError(
                f"PositionField: last dim must be 3, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

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
    def n_atoms(self) -> int:
        return self._data.shape[0]

    def __repr__(self) -> str:
        return f"PositionField(shape={self._data.shape}, law=\"p'=R p+t\")"


class VectorField:
    """3D polar vector field.

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


class QualifiedVectorField(VectorField):
    """Cartesian triples whose producer-level vector law is conditional.

    The stored component order is still ``x,y,z`` and the conditional
    physical representation is polar ``1o``.  However, a sanitizer,
    finite-grid estimator, learned charge source, or independently rerun
    calculator prevents the SDK from advertising an unconditional e3nn
    irrep for the end-to-end emitted feature. Read the catalog entry's
    ``transformation`` and ``validity`` before using the components.
    """

    IRREPS = None
    CONDITIONAL_IRREPS = Irreps("1x1o")

    @property
    def irreps(self):
        """No unconditional producer-level irrep is claimed."""
        return self.IRREPS

    @property
    def conditional_irreps(self) -> Irreps:
        """Polar irrep for the catalog entry's explicitly stated condition."""
        return self.CONDITIONAL_IRREPS

    def __repr__(self) -> str:
        return (
            f"QualifiedVectorField(shape={self._data.shape}, "
            "irreps=None, conditional_irreps='1x1o')"
        )


class MagneticVectorField(VectorField):
    """3D axial vector field for magnetic B-field diagnostics.

    Shape ``(*, 3)``.  Cartesian (x, y, z) components with even parity.
    """

    IRREPS = Irreps("1x1e")

    def __repr__(self) -> str:
        return f"MagneticVectorField(shape={self._data.shape}, irreps='{self.IRREPS}')"


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


class PerRingTypeT1:
    """T1 (L=1 axial) contribution decomposed by ring type.

    Shape ``(*, 24)`` = 8 ring types x 3 T1 components.
    Use :meth:`as_block` for ``(*, 8, 3)`` view.
    """

    __slots__ = ("_data",)
    IRREPS = Irreps("8x1e")

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_RING_TYPES * 3:
            raise ValueError(
                f"PerRingTypeT1: last dim must be {N_RING_TYPES * 3}, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def irreps(self) -> Irreps:
        return self.IRREPS

    def for_type(self, rt: RingType) -> np.ndarray:
        i = int(rt)
        return self._data[..., i * 3:(i + 1) * 3]

    @property
    def total(self) -> np.ndarray:
        return self.as_block().sum(axis=-2)

    def as_block(self) -> np.ndarray:
        """Reshape to ``(*, 8, 3)``."""
        return self._data.reshape(*self._data.shape[:-1], N_RING_TYPES, 3)

    def __repr__(self) -> str:
        return f"PerRingTypeT1(shape={self._data.shape})"


class PerRingTypeT2:
    """Project-native T2 contribution decomposed by ring type.

    Shape ``(*, 40)`` = 8 ring types x 5 T2 components.
    Use :meth:`as_block` for ``(*, 8, 5)`` view.
    """

    __slots__ = ("_data",)
    TENSOR_BASIS = PROJECT_T2_TENSOR_BASIS
    COMPONENT_ORDER = PROJECT_T2_COMPONENT_ORDER
    TENSOR_FRAME = PROJECT_T2_TENSOR_FRAME
    E3NN_EXPORT = PROJECT_E3NN_EXPORT_NOTE

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != N_RING_TYPES * 5:
            raise ValueError(
                f"PerRingTypeT2: last dim must be {N_RING_TYPES * 5}, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def tensor_basis(self) -> str:
        return self.TENSOR_BASIS

    @property
    def component_order(self) -> str:
        return self.COMPONENT_ORDER

    @property
    def tensor_frame(self) -> str:
        return self.TENSOR_FRAME

    @property
    def e3nn_export(self) -> str:
        return self.E3NN_EXPORT

    def for_type(self, rt: RingType) -> np.ndarray:
        i = int(rt)
        return self._data[..., i * 5:(i + 1) * 5]

    @property
    def total(self) -> np.ndarray:
        return self.as_block().sum(axis=-2)

    def as_block(self) -> np.ndarray:
        """Reshape to ``(N, 8, 5)``."""
        return self._data.reshape(*self._data.shape[:-1], N_RING_TYPES, 5)

    def to_e3nn(self) -> E3nnTensor:
        converted = project_t2_to_e3nn(self.as_block()).reshape(self._data.shape)
        return E3nnTensor(converted, Irreps(f"{N_RING_TYPES}x2e"))

    def __repr__(self) -> str:
        return f"PerRingTypeT2(shape={self._data.shape}, basis='{self.TENSOR_BASIS}')"


# ── Per-bond-category decomposition ─────────────────────────────────


class PerBondCategoryT2:
    """Project-native T2 McConnell contribution decomposed by bond category.

    Shape ``(*, 25)`` = 5 bond categories x 5 T2 components.
    Use :meth:`as_block` for ``(*, 5, 5)`` view.
    """

    __slots__ = ("_data",)
    TENSOR_BASIS = PROJECT_T2_TENSOR_BASIS
    COMPONENT_ORDER = PROJECT_T2_COMPONENT_ORDER
    TENSOR_FRAME = PROJECT_T2_TENSOR_FRAME
    E3NN_EXPORT = PROJECT_E3NN_EXPORT_NOTE

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
    def tensor_basis(self) -> str:
        return self.TENSOR_BASIS

    @property
    def component_order(self) -> str:
        return self.COMPONENT_ORDER

    @property
    def tensor_frame(self) -> str:
        return self.TENSOR_FRAME

    @property
    def e3nn_export(self) -> str:
        return self.E3NN_EXPORT

    def for_category(self, cat: BondCategory) -> np.ndarray:
        i = int(cat)
        return self._data[..., i * 5:(i + 1) * 5]

    @property
    def total(self) -> np.ndarray:
        return self.as_block().sum(axis=-2)

    def as_block(self) -> np.ndarray:
        return self._data.reshape(*self._data.shape[:-1], N_BOND_CATEGORIES, 5)

    def to_e3nn(self) -> E3nnTensor:
        converted = project_t2_to_e3nn(self.as_block()).reshape(self._data.shape)
        return E3nnTensor(converted, Irreps(f"{N_BOND_CATEGORIES}x2e"))

    def __repr__(self) -> str:
        return f"PerBondCategoryT2(shape={self._data.shape}, basis='{self.TENSOR_BASIS}')"


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


class McConnellNearFieldCounts:
    """(*, 2) McConnell near-field audit counts below 3 A.

    Columns are accepted and rejected source-target pair counts after the
    same filters used by the McConnell kernel evaluation path.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 2:
            raise ValueError(
                f"McConnellNearFieldCounts: last dim must be 2, got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def accepted_lt3A(self) -> np.ndarray:
        return self._data[..., 0]

    @property
    def rejected_lt3A(self) -> np.ndarray:
        return self._data[..., 1]

    def __repr__(self) -> str:
        return f"McConnellNearFieldCounts(shape={self._data.shape})"


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
    """(*, 4) Coulomb E-field summary. Units: V/A.

    ``E_bond_proj`` is defined only for hydrogens with a valid typed
    parent atom. It is NaN for non-hydrogen atoms and for parentless
    hydrogens. Column 2 is the signed projection of the backbone field onto
    the total-field direction in V/A; despite its historical ``*_frac`` name,
    it is not a ratio and is set to zero when the total field is near zero.
    """

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
        """Parent-to-H projection for supported H atoms; NaN otherwise."""
        return self._data[..., 1]

    @property
    def E_backbone_frac(self) -> np.ndarray:
        """Historical alias for :attr:`E_backbone_projection`; not a fraction."""
        return self._data[..., 2]

    @property
    def E_backbone_projection(self) -> np.ndarray:
        """Signed ``E_backbone · unit(E_total)`` in V/A."""
        return self._data[..., 2]

    @property
    def aromatic_E_magnitude(self) -> np.ndarray:
        return self._data[..., 3]

    def __repr__(self) -> str:
        return f"CoulombScalars(shape={self._data.shape})"


class HBondScalars:
    """(*, 4) H-bond summary with mixed units.

    Columns are nearest accepted distance (Å), its inverse cube (Å⁻³), count
    of accepted source evaluations within ``hbond_counting_radius`` (default
    3.5 Å), and ``Σ(3cos²θ−1)/r³`` (Å⁻³) over every accepted evaluation. The
    count therefore need not equal the number summed. A legacy zero
    distance/inverse cube means no accepted source when the companion
    ``hbond_flags`` availability is zero.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 4:
            raise ValueError(f"HBondScalars: last dim must be 4, got {data.shape}")
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

    @property
    def mcconnell_scalar(self) -> np.ndarray:
        return self._data[..., 3]

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
    """``(*, 4)`` legacy MOPAC printed-format per-atom projection.

    Columns are charge (six decimals), s and p populations (five decimals),
    and the sum of retained three-decimal orders from the reconstructed
    first-six compact bond row. The values are reconstructed from libmopac
    structs; this wrapper does not expose the separate full-precision arrays.
    """

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
        """Legacy sum of first-six/F6.3/>0.01 compact Wiberg orders."""
        return self._data[..., 3]

    def __repr__(self) -> str:
        return f"MopacScalars(shape={self._data.shape})"


class MopacGlobal:
    """``(4,)`` legacy printed MOPAC graph projection.

    The heat of formation is represented at five decimal places and the three
    dipole components at three. Unquantized binary64 API values are available
    in the direct MOPAC group ("full precision" means no compatibility decimal
    projection, not physical exactness). The molecular dipole is a polar vector;
    for a charged system it depends on the MOPAC input origin/translation.
    """

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


class MopacAtomPopulations:
    """``(N, 12)`` legacy printed-format population projection.

    Columns 0:5 are charge/electron-trace/s/p/d at 6/4/5/5/5 decimal
    places. p and d are NaN when their shells are not live. Column 10 is the
    three-decimal CSC-diagonal valency; column 11 is the sum of retained
    first-six/F6.3 compact bond orders. The legacy f-population column 5 and
    per-atom-dipole columns 6:10 have no struct source and are explicit NaN.
    All populated values are reconstructed from libmopac structs.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 12:
            raise ValueError(f"MopacAtomPopulations: expected (N, 12), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def net_charge(self) -> np.ndarray:
        return self._data[:, 0]

    @property
    def electron_density(self) -> np.ndarray:
        """Atomic electron population (AO-density trace)."""
        return self._data[:, 1]

    @property
    def s_population(self) -> np.ndarray:
        return self._data[:, 2]

    @property
    def p_population(self) -> np.ndarray:
        return self._data[:, 3]

    @property
    def d_population(self) -> np.ndarray:
        return self._data[:, 4]

    @property
    def f_population(self) -> np.ndarray:
        return self._data[:, 5]

    @property
    def dipole_contribution(self) -> np.ndarray:
        return self._data[:, 6:10]

    @property
    def mopac_valency(self) -> np.ndarray:
        return self._data[:, 10]

    @property
    def project_valency(self) -> np.ndarray:
        return self._data[:, 11]


class MopacAtomicOrbitalPopulations:
    """``(N, 9)`` legacy F10.5 AO-population projection.

    Values are reconstructed from the libmopac atom-AO-density diagonal and
    quantized to five decimal places. Only each atom's live AO width is
    populated; non-existent per-atom AO columns are NaN rather than defensive
    zeros. The direct MOPAC group retains the unquantized padded API density
    blocks.
    These diagonal-only s/px/py/pz/d values are not a closed O(3)
    representation because the AO-density coherences are omitted.
    """

    __slots__ = ("_data",)
    COLUMNS = ("s", "px", "py", "pz", "x2_minus_y2", "xz", "z2", "yz", "xy")
    diagnostic_frame_dependent = True
    model_scalar_source = False

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 9:
            raise ValueError(
                f"MopacAtomicOrbitalPopulations: expected (N, 9), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def s(self) -> np.ndarray:
        return self._data[:, 0]

    @property
    def p(self) -> np.ndarray:
        return self._data[:, 1:4]

    @property
    def d(self) -> np.ndarray:
        return self._data[:, 4:9]


class MopacAtomicOrbitalPopulationTotals:
    """``(N, 3)`` legacy AO-table-derived shell totals [s, p, d].

    Each AO diagonal is first projected through its F10.5 printed value and
    those projected entries are then summed, exactly matching the old parser;
    this is distinct from independently projecting the finished shell sum.
    p/d are NaN where the corresponding shell is not live; unquantized
    binary64 API-derived sums are available in the direct MOPAC group.
    """

    __slots__ = ("_data",)
    COLUMNS = ("s_total", "p_total", "d_total")
    diagnostic_frame_dependent = False
    model_scalar_source = True

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(
                f"MopacAtomicOrbitalPopulationTotals: expected (N, 3), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def s_total(self) -> np.ndarray:
        return self._data[:, 0]

    @property
    def p_total(self) -> np.ndarray:
        return self._data[:, 1]

    @property
    def d_total(self) -> np.ndarray:
        return self._data[:, 2]


class MopacUniqueBondOrders:
    """``(U, 8)`` legacy-compatible API-unique Wiberg projection.

    The a<b row set is the complete sparse API pair set, independent of the
    first-six compact table. Each symmetric CSC entry is first projected
    through F6.3 and max/mean are then formed as in the old parser. Columns 4, 5, and 6
    are NaN because the former ALLBONDS printed-entry count/indices have no
    struct source; column 7 is the topology-bond index in this API-unique row
    basis.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 8:
            raise ValueError(f"MopacUniqueBondOrders: expected (U, 8), got {data.shape}")
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
    def max_order(self) -> np.ndarray:
        return self._data[:, 2]

    @property
    def mean_order(self) -> np.ndarray:
        return self._data[:, 3]

    @property
    def topology_bond_index(self) -> np.ndarray:
        return self._data[:, 7]


class MopacTopologyBondOrdersFull:
    """``(B, 8)`` legacy-compatible topology/API-pair bridge.

    Present orders are reconstructed from the sparse CSC pair and quantized to
    three decimals; unique-pair indices use the complete API-unique row basis,
    not the first-six compact subset. The final ALLBONDS printed-entry-count
    column is text-only and always NaN. API absence includes Wiberg values
    omitted at MOPAC's 0.01 sparse threshold.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 8:
            raise ValueError(f"MopacTopologyBondOrdersFull: expected (B, 8), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def order(self) -> np.ndarray:
        return self._data[:, 3]

    @property
    def present(self) -> np.ndarray:
        return self._data[:, 4].astype(bool)


class BondOrders:
    """(K, 3) sparse pair rows: [atom_i, atom_j, order].

    ``K`` is the caller-defined serialized pair axis (for example compact or
    retained libmopac-unique), not necessarily the topology bond count.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(f"BondOrders: expected (K, 3), got {data.shape}")
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
        """Materialize a symmetric matrix with representational zero fill.

        For sparse MOPAC surfaces, an unlisted pair may have been censored by
        the API threshold and/or compact first-six projection. The dense zero
        is therefore not evidence that the physical Wiberg order was measured
        as zero.
        """
        mat = np.zeros((n_atoms, n_atoms), dtype=np.float64)
        i, j = self.atom_i, self.atom_j
        mat[i, j] = self._data[:, 2]
        mat[j, i] = self._data[:, 2]
        return mat

    def __repr__(self) -> str:
        return f"BondOrders(n_bonds={self.n_bonds})"


class DeltaScalars:
    """(*, 6) WT-minus-mutant metadata on dense WT atom rows.

    Columns are match mask, shielding T0 delta, nearest removed-ring distance,
    configured-force-field charge delta, compatibility MOPAC charge delta, and
    typed match distance. A matched row with no removed ring uses 99 Å.
    """

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
        """Configured ChargeSource partial-charge delta (WT - mutant)."""
        return self._data[..., 3]

    @property
    def delta_mopac_charge(self) -> np.ndarray:
        """Legacy F15.6 MOPAC Coulson charge delta.

        Zero is unmasked and ambiguous between unavailable MOPAC and a true
        zero compatibility-charge delta.
        """
        return self._data[..., 4]

    @property
    def match_distance(self) -> np.ndarray:
        """Spatial distance between matched WT/mutant atoms (Å)."""
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
    The entire NPY is absent when ``MutationDeltaResult`` has no paired APBS
    delta; absence is not represented by a zero-filled file.
    """

    __slots__ = ("_data",)
    efg_structural_zero_components = ("T0", "T1_x", "T1_y", "T1_z")
    tensor_basis = PROJECT_T2_TENSOR_BASIS

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

    @property
    def delta_efg_full9(self) -> "SphericalTensor":
        """Alias for the legacy full-9 EFG envelope."""
        return self.delta_efg

    @property
    def delta_efg_t2(self) -> EFGTensor:
        """Physical APBS delta EFG T2 slice from columns 7:12."""
        return EFGTensor(self._data[..., 7:12])

    def __repr__(self) -> str:
        return f"DeltaAPBS(shape={self._data.shape})"


class DeltaRingProximity:
    """(*, R*6) removed-ring geometry on dense WT atom rows.

    Per ring: ``[dist, z, rho, theta, mcconnell, exp_decay]``. The entire NPY
    is absent unless a mutation site removes at least one WT ring; unmatched
    rows in a present file are compatibility zeros gated by ``delta_scalars``.
    """

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
    """``(N,)`` per-atom ``charges`` output from the loaded TorchScript model.

    A conforming AIMNet2 model intends these to be Hirshfeld charges. The
    producer accepts a caller-selected ``.jpt`` and checks the runtime
    interface, shape and finiteness, but does not identify a particular
    functional or certify scalar/equivariant behaviour. The float64 NPY
    widens the Torch result; it does not mean the model ran in binary64.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim != 1:
            raise ValueError(
                f"AIMNet2Charges: expected (N,), got {data.shape}")
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
    """``(N, 256)`` raw float32 ``aim`` output from the loaded model.

    These are opaque learned channels. The producer does not attach a
    channel-level physical interpretation or guarantee a Cartesian, O(3),
    parity, or e3nn transformation law. Keep the 256-channel axis fixed; do
    not rotate it as spatial components. The values are geometry-dependent
    and finite-validated as part of the all-or-nothing main AIMNet2 result.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 256:
            raise ValueError(
                f"AIMNet2AimEmbedding: expected (N, 256), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def __repr__(self) -> str:
        return f"AIMNet2AimEmbedding(shape={self._data.shape})"


class AIMNet2ChargeResponseGradient:
    """(N, 3) AIMNet2 charge-response proxy/diagnostic via autograd.

    The vector field is dL/d(r_i) where L = sum_j q_j^2 over
    non-sentinel atoms (the L2-of-charges objective; for a conforming AIMNet2
    model, its charge-conservation projection makes sum(q) constant and that
    alternative gradient ~0). It is a sum-of-squared-charges sensitivity to
    atomic coordinates, NOT a Buckingham polarizability α = ∂μ/∂E
    and NOT an atom-resolved charge Jacobian.

    Companion scalar is the L2 norm of the vector, stored separately
    in `aimnet2_charge_response_gradient_scalar.npy`.

    Conditional law only: if the loaded model makes ``sum_j q_j^2`` an O(3)-
    invariant scalar, this coordinate gradient is a Cartesian polar vector.
    The producer accepts a caller-selected ``.jpt`` and validates autograd
    shape/finiteness, not that scalar-invariance property, so no unconditional
    e3nn irrep or exact float32 rerun covariance is advertised here.
    """

    __slots__ = ("_data",)
    IRREPS = None
    CONDITIONAL_IRREPS = Irreps("1x1o")

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(
                f"AIMNet2ChargeResponseGradient: expected (N, 3), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def irreps(self):
        """No unconditional irrep is producer-certified for an arbitrary model."""
        return self.IRREPS

    @property
    def conditional_irreps(self) -> Irreps:
        """Polar-vector irrep if the loaded model's charge objective is scalar."""
        return self.CONDITIONAL_IRREPS

    @property
    def vectors(self) -> np.ndarray:
        """Per-atom charge-response gradient vectors ``(N, 3)``."""
        return self._data

    @property
    def norms(self) -> np.ndarray:
        """Per-atom L2 norm ``(N,)``. Equal to ``aimnet2_charge_response_gradient_scalar.npy``
        up to floating-point precision; this property recomputes from the
        vectors so it's always consistent with the loaded vector field."""
        return np.linalg.norm(self._data, axis=-1)

    def __repr__(self) -> str:
        return f"AIMNet2ChargeResponseGradient(shape={self._data.shape})"
