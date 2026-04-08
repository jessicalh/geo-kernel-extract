"""Per-ring sparse contributions and ring geometry reference.

ring_contributions.npy — (P, 57) one row per (atom, ring) pair:
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
    [27:36] hm_G                SphericalTensor — HM shielding kernel (intensity * H)
    [36:45] pq_G                SphericalTensor
    [45:54] chi_G               SphericalTensor
    [54]    disp_scalar
    [55]    disp_contacts
    [56]    gaussian_density

ring_geometry.npy — (R, 10) one row per ring:
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
    """Sparse per-(atom, ring) pair contributions.

    Shape ``(P, 57)`` where P = number of evaluated (atom, ring) pairs.
    Each row carries geometry, all ring current kernels as
    :class:`SphericalTensor` views, and dispersion scalars.

    Use :meth:`for_atom` and :meth:`for_ring_type` to filter rows.

    Args:
        data: numpy array of shape ``(P, 57)``.
    """

    __slots__ = ("_data",)
    COLS = 57

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
        """Haigh-Mallion shielding kernel G (intensity * H). Cols 27-35."""
        return SphericalTensor(self._data[:, 27:36])

    @property
    def pq(self) -> SphericalTensor:
        """Pi-quadrupole kernel. Cols 36-44."""
        return SphericalTensor(self._data[:, 36:45])

    @property
    def chi(self) -> SphericalTensor:
        """Ring susceptibility kernel. Cols 45-53."""
        return SphericalTensor(self._data[:, 45:54])

    @property
    def disp_scalar(self) -> np.ndarray:
        return self._data[:, 54]

    @property
    def disp_contacts(self) -> np.ndarray:
        return self._data[:, 55].astype(np.intp)

    @property
    def gaussian_density(self) -> np.ndarray:
        return self._data[:, 56]

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
    """Per-ring geometry reference table.

    Shape ``(R, 10)`` where R = number of rings in the protein.
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
