"""Per-ring sparse contributions and ring geometry reference.

ring_contributions.npy — (P, 40) one row per union (atom, aromatic-ring)
neighbour row. Different calculators can contribute different subsets; an
untouched calculator block is its default zero and has no block-specific mask:
    [0]     atom_index
    [1]     ring_index
    [2]     ring_type           RingTypeIndex 0-7
    [3]     distance            A
    [4]     rho                 A
    [5]     z                   A (signed)
    [6]     theta               rad
    [7]     mcconnell_factor    (3cos^2 theta - 1) / r^3
    [8]     exp_decay           exp(-distance / 4.0)
    [9:18]  bs_G                SphericalTensor — BS shielding kernel
    [18:27] hm_H                SphericalTensor — HM raw integral (pure T2)
    [27:36] hm_G                SphericalTensor — HM unscaled shielding kernel G
                            (literature intensity is applied downstream)
    [36]    disp_scalar
    [37]    disp_contacts
    [38]    cos_phi             azimuthal angle cosine (relative to vertex 0)
    [39]    sin_phi             azimuthal angle sine (relative to vertex 0)

ring_geometry.npy — (R, 10) one row per aromatic ring:
    [0]     ring_index
    [1]     ring_type
    [2]     residue_index
    [3:6]   center (x, y, z)
    [6:9]   normal (nx, ny, nz)
    [9]     radius
"""

from __future__ import annotations

import numpy as np

from ._types import RingType
from ._tensors import SphericalTensor


class RingContributions:
    """Sparse union of per-(atom, aromatic-ring) neighbour rows.

    Shape ``(P, 40)``. Each row carries shared geometry plus calculator-owned
    ring-current/dispersion blocks. The axis is a union: a calculator may not
    have accepted a row that another calculator created, and its untouched
    zero is therefore absence-or-physical-zero without a per-calculator mask.

    Use :meth:`for_atom` and :meth:`for_ring_type` to filter rows.

    Args:
        data: numpy array of shape ``(P, 40)``.
    """

    __slots__ = ("_data",)
    COLS = 40

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != self.COLS:
            raise ValueError(
                f"RingContributions: expected (P, {self.COLS}), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def n_pairs(self) -> int:
        return self._data.shape[0]

    # ── Index columns ───────────────────────────────────────────────

    @property
    def atom_index(self) -> np.ndarray:
        return self._data[:, 0].astype(np.intp)

    @property
    def ring_index(self) -> np.ndarray:
        return self._data[:, 1].astype(np.intp)

    @property
    def ring_type(self) -> np.ndarray:
        return self._data[:, 2].astype(np.intp)

    # ── Geometry ────────────────────────────────────────────────────

    @property
    def distance(self) -> np.ndarray:
        return self._data[:, 3]

    @property
    def rho(self) -> np.ndarray:
        return self._data[:, 4]

    @property
    def z(self) -> np.ndarray:
        return self._data[:, 5]

    @property
    def theta(self) -> np.ndarray:
        return self._data[:, 6]

    @property
    def mcconnell_factor(self) -> np.ndarray:
        return self._data[:, 7]

    @property
    def ring_chi_scalar(self) -> np.ndarray:
        """Preferred geometry-scalar alias for legacy column 7.

        The separately emitted ``ringchi_scalar.npy`` carries the production
        RingSusceptibilityResult field on the same union axis, but its default
        zero is likewise ambiguous when that calculator did not accept a row.
        """
        return self._data[:, 7]

    @property
    def exp_decay(self) -> np.ndarray:
        return self._data[:, 8]

    # ── Physics kernels (SphericalTensor views) ─────────────────────

    @property
    def bs(self) -> SphericalTensor:
        """Biot-Savart shielding kernel G. Cols 9-17."""
        return SphericalTensor(self._data[:, 9:18])

    @property
    def hm_H(self) -> SphericalTensor:
        """Haigh-Mallion raw surface integral H (pure T2). Cols 18-26."""
        return SphericalTensor(self._data[:, 18:27])

    @property
    def hm(self) -> SphericalTensor:
        """Haigh-Mallion unscaled shielding kernel G; intensity is downstream. Cols 27-35."""
        return SphericalTensor(self._data[:, 27:36])

    @property
    def disp_scalar(self) -> np.ndarray:
        return self._data[:, 36]

    @property
    def disp_contacts(self) -> np.ndarray:
        return self._data[:, 37].astype(np.intp)

    @property
    def cos_phi(self) -> np.ndarray:
        """Cosine of azimuthal angle in ring plane (relative to vertex 0)."""
        return self._data[:, 38]

    @property
    def sin_phi(self) -> np.ndarray:
        """Sine of azimuthal angle in ring plane (relative to vertex 0)."""
        return self._data[:, 39]

    # ── Filtering ───────────────────────────────────────────────────

    def for_atom(self, idx: int) -> RingContributions:
        mask = self._data[:, 0].astype(np.intp) == idx
        return RingContributions(self._data[mask])

    def for_ring_type(self, rt: RingType) -> RingContributions:
        mask = self._data[:, 2].astype(np.intp) == int(rt)
        return RingContributions(self._data[mask])

    def __repr__(self) -> str:
        return f"RingContributions(n_pairs={self.n_pairs})"


class RingGeometry:
    """Per-aromatic-ring geometry reference table.

    Shape ``(R, 10)`` where R is ``Protein::RingCount()`` (aromatic only).
    One row per ring: index, type, residue, center, normal, radius.

    Args:
        data: numpy array of shape ``(R, 10)``.
    """

    __slots__ = ("_data",)
    COLS = 10

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != self.COLS:
            raise ValueError(
                f"RingGeometry: expected (R, {self.COLS}), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def n_rings(self) -> int:
        return self._data.shape[0]

    @property
    def ring_index(self) -> np.ndarray:
        return self._data[:, 0].astype(np.intp)

    @property
    def ring_type(self) -> np.ndarray:
        return self._data[:, 1].astype(np.intp)

    @property
    def residue_index(self) -> np.ndarray:
        return self._data[:, 2].astype(np.intp)

    @property
    def center(self) -> np.ndarray:
        return self._data[:, 3:6]

    @property
    def normal(self) -> np.ndarray:
        return self._data[:, 6:9]

    @property
    def radius(self) -> np.ndarray:
        return self._data[:, 9]

    def __repr__(self) -> str:
        return f"RingGeometry(n_rings={self.n_rings})"


class RingPairGeometry:
    """All ``i < j`` aromatic-ring geometry rows, shape ``(P, 13)``."""

    __slots__ = ("_data",)
    COLS = 13

    def __init__(self, data: np.ndarray):
        if data.ndim != 2 or data.shape[1] != self.COLS:
            raise ValueError(
                f"RingPairGeometry: expected (P, {self.COLS}), got {data.shape}")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def ring_a(self) -> np.ndarray:
        return self._data[:, 0].astype(np.intp)

    @property
    def ring_b(self) -> np.ndarray:
        return self._data[:, 1].astype(np.intp)

    @property
    def center_distance(self) -> np.ndarray:
        return self._data[:, 6]

    @property
    def normal_dot(self) -> np.ndarray:
        return self._data[:, 7]

    @property
    def normal_cross_magnitude(self) -> np.ndarray:
        return self._data[:, 8]

    @property
    def is_fused(self) -> np.ndarray:
        return self._data[:, 12] != 0

    def __repr__(self) -> str:
        return f"RingPairGeometry(n_pairs={self._data.shape[0]})"
