"""Trajectory data reader — loads the H5 master file from a GROMACS scan.

    from nmr_extract import load_trajectory
    traj = load_trajectory("output/trajectory.h5")

    traj.positions          # (T, N, 3) float64 — per-frame xyz
    traj.frame_times        # (T,) float64 — time in ps
    traj.n_atoms            # int
    traj.n_frames           # int

    # Rollup: legacy GromacsProtein per-atom Welford columns
    traj.rollup.bs_T0.mean  # (N,) ring current isotropic mean
    traj.rollup.bs_T0.std   # (N,) ring current isotropic std
    traj.rollup.names       # list of column names

    # Bonds: per-bond length statistics
    traj.bonds.atom_a       # (B,) atom indices
    traj.bonds.atom_b       # (B,) atom indices
    traj.bonds.length_mean  # (B,) mean bond length (A)
    traj.bonds.length_std   # (B,) bond length std

    # Welford H5 groups — written by *WelfordTrajectoryResult subclasses
    # (BS / HM / McConnell / Eeq / Sasa / HBondCount), 2026-05-17/18.
    # Each group is Optional — None when the corresponding C++ TR was
    # not attached for the run that produced this trajectory.h5.
    traj.welford.bs.t0.mean                  # (N,) ppm_T_per_nA
    traj.welford.bs.t1.mean                  # (N, 3) Cartesian Levi-Civita dual
    traj.welford.bs.t2.mean                  # (N, 5) real-spherical-tesseral m-basis
    traj.welford.bs.t0_dxdt.mean             # (N,) cadence-normalized rate
    traj.welford.bs.mean_dt_ps               # scalar — trajectory cadence
    traj.welford.bs.frame_index_range        # (first, last)
    traj.welford.eeq.charge.mean             # (N,) elementary_charge
    traj.welford.eeq.charge_dxdt.units       # "elementary_charge_per_ps"
    traj.welford.hbond_count.occupancy_fraction.mean  # (N,) ∈ [0,1]
    traj.welford.hbond_count.occupancy_fraction.units # "dimensionless"
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


# ─── Legacy GromacsProtein rollup ──────────────────────────────────


@dataclass(frozen=True)
class WelfordColumn:
    """Mean and std for one accumulated quantity across all atoms."""
    name: str
    mean: np.ndarray    # (N,)
    std: np.ndarray     # (N,)


class TrajectoryRollup:
    """Named access to per-atom Welford statistics.

    Access by attribute: rollup.bs_T0.mean, rollup.aimnet2_charge.std
    Access by index: rollup[0].mean (same as rollup.bs_T0.mean)
    """

    def __init__(self, names: list[str],
                 means: np.ndarray, stds: np.ndarray):
        self._names = names
        self._means = means   # (N, K)
        self._stds = stds     # (N, K)
        self._by_name: dict[str, int] = {n: i for i, n in enumerate(names)}

    @property
    def names(self) -> list[str]:
        return list(self._names)

    @property
    def n_columns(self) -> int:
        return len(self._names)

    def __getattr__(self, name: str) -> WelfordColumn:
        if name.startswith("_"):
            raise AttributeError(name)
        if name not in self._by_name:
            raise AttributeError(
                f"No rollup column '{name}'. "
                f"Available: {', '.join(self._names[:10])}...")
        k = self._by_name[name]
        return WelfordColumn(
            name=name,
            mean=self._means[:, k],
            std=self._stds[:, k],
        )

    def __getitem__(self, k: int) -> WelfordColumn:
        return WelfordColumn(
            name=self._names[k],
            mean=self._means[:, k],
            std=self._stds[:, k],
        )

    def as_block(self, columns: Optional[list[str]] = None) -> np.ndarray:
        """Return (N, K) mean array for selected columns (default: all)."""
        if columns is None:
            return self._means.copy()
        idxs = [self._by_name[c] for c in columns]
        return self._means[:, idxs]


@dataclass(frozen=True)
class BondRollup:
    """Per-bond length statistics across trajectory frames."""
    atom_a: np.ndarray       # (B,) uint64
    atom_b: np.ndarray       # (B,) uint64
    length_mean: np.ndarray  # (B,) float64
    length_std: np.ndarray   # (B,) float64

    @property
    def n_bonds(self) -> int:
        return len(self.atom_a)


# ─── Welford H5 groups (Phase 2b/C, 2026-05-17/18) ──────────────────
#
# One WelfordMoments per per-atom channel. Each Welford TR group
# exposes its channels as named attributes; shape is (N,) for scalar
# channels, (N, K) for per-component channels (T1 K=3, T2 K=5).


@dataclass(frozen=True)
class WelfordMoments:
    """7-stat block for one Welford channel.

    Mirrors the H5 datasets `<prefix>_{mean,m2,std,min,max,min_frame,
    max_frame}` plus the per-dataset `units` attribute. Per-component
    channels (T1, T2) carry (N, K) arrays; scalar channels carry (N,).
    """
    mean: np.ndarray
    m2: np.ndarray
    std: np.ndarray
    min: np.ndarray
    max: np.ndarray
    min_frame: np.ndarray   # uint frame index
    max_frame: np.ndarray   # uint frame index
    units: str              # base unit string from per-dataset H5 attr


def _read_moments(grp, prefix: str) -> WelfordMoments:
    """Read a 7-stat Welford block from an HDF5 group.

    Reads `<prefix>_{mean, m2, std, min, max, min_frame, max_frame}`
    and the `units` attribute from `<prefix>_mean` (authoritative
    per-dataset). Works for both 1D scalar channels and 2D
    per-component channels — numpy carries the shape through.
    """
    ds_mean = grp[f"{prefix}_mean"]
    units = ds_mean.attrs.get("units", "")
    if isinstance(units, bytes):
        units = units.decode()
    return WelfordMoments(
        mean=ds_mean[:],
        m2=grp[f"{prefix}_m2"][:],
        std=grp[f"{prefix}_std"][:],
        min=grp[f"{prefix}_min"][:],
        max=grp[f"{prefix}_max"][:],
        min_frame=grp[f"{prefix}_min_frame"][:],
        max_frame=grp[f"{prefix}_max_frame"][:],
        units=str(units),
    )


def _decode_attr(value):
    """Decode bytes → str (HDF5 round-trips string attrs as bytes)."""
    if isinstance(value, bytes):
        return value.decode()
    return value


def _group_units(grp) -> str:
    return str(_decode_attr(grp.attrs.get("units", "")))


def _group_mean_dt_ps(grp) -> float:
    return float(grp.attrs.get("mean_dt_ps", 0.0))


def _group_frame_index_range(grp) -> tuple[int, int]:
    arr = grp.attrs.get("frame_index_range", None)
    if arr is None:
        return (0, 0)
    arr = np.asarray(arr).ravel()
    return (int(arr[0]), int(arr[1]))


# ─── Tensor Welford groups (BS / HM / Mc) ───────────────────────────


@dataclass(frozen=True)
class _TensorWelfordGroup:
    """Common shape for BS / HM / McConnell Welford rollups."""
    # Scalar channels
    t0: WelfordMoments               # (N,)
    t2magnitude: WelfordMoments      # (N,)
    # Per-component channels
    t1: WelfordMoments               # (N, 3) — Cartesian Levi-Civita dual
    t2: WelfordMoments               # (N, 5) — real-spherical-tesseral m-basis
    # Delta channels (on T0)
    t0_delta: WelfordMoments
    t0_abs_delta: WelfordMoments
    t0_delta_squared: WelfordMoments
    t0_dxdt: WelfordMoments
    # Provenance
    t0_rms_delta: np.ndarray         # (N,) finalize-derived sqrt(<Δ²>)
    n_frames_per_atom: np.ndarray
    delta_n_per_atom: np.ndarray
    # Group-level attributes
    mean_dt_ps: float
    frame_index_range: tuple[int, int]
    irrep_layout_t1: str             # "v_x,v_y,v_z" — Cartesian
    irrep_layout_t2: str             # "m-2,m-1,m0,m+1,m+2"
    units: str                       # primary value-channel unit


@dataclass(frozen=True)
class BsWelfordGroup(_TensorWelfordGroup):
    """Per-atom Biot-Savart Welford rollup from /trajectory/bs_welford/."""


@dataclass(frozen=True)
class HmWelfordGroup(_TensorWelfordGroup):
    """Per-atom Haigh-Mallion Welford rollup from /trajectory/hm_welford/."""


@dataclass(frozen=True)
class McConnellWelfordGroup(_TensorWelfordGroup):
    """Per-atom McConnell Welford rollup from /trajectory/mc_welford/."""


def _load_tensor_welford(f, h5_path: str):
    """Read a tensor Welford group (BS / HM / Mc); return raw kwargs dict."""
    if h5_path not in f:
        return None
    grp = f[h5_path]
    return dict(
        t0=_read_moments(grp, "t0"),
        t2magnitude=_read_moments(grp, "t2magnitude"),
        t1=_read_moments(grp, "t1"),
        t2=_read_moments(grp, "t2"),
        t0_delta=_read_moments(grp, "t0_delta"),
        t0_abs_delta=_read_moments(grp, "t0_abs_delta"),
        t0_delta_squared=_read_moments(grp, "t0_delta_squared"),
        t0_dxdt=_read_moments(grp, "t0_dxdt"),
        t0_rms_delta=grp["t0_rms_delta"][:],
        n_frames_per_atom=grp["n_frames_per_atom"][:],
        delta_n_per_atom=grp["delta_n_per_atom"][:],
        mean_dt_ps=_group_mean_dt_ps(grp),
        frame_index_range=_group_frame_index_range(grp),
        irrep_layout_t1=str(_decode_attr(
            grp.attrs.get("irrep_layout_t1", ""))),
        irrep_layout_t2=str(_decode_attr(
            grp.attrs.get("irrep_layout_t2", ""))),
        units=_group_units(grp),
    )


def _load_bs_welford(f) -> Optional[BsWelfordGroup]:
    kwargs = _load_tensor_welford(f, "/trajectory/bs_welford")
    return None if kwargs is None else BsWelfordGroup(**kwargs)


def _load_hm_welford(f) -> Optional[HmWelfordGroup]:
    kwargs = _load_tensor_welford(f, "/trajectory/hm_welford")
    return None if kwargs is None else HmWelfordGroup(**kwargs)


def _load_mc_welford(f) -> Optional[McConnellWelfordGroup]:
    kwargs = _load_tensor_welford(f, "/trajectory/mc_welford")
    return None if kwargs is None else McConnellWelfordGroup(**kwargs)


# ─── Scalar Welford groups (Eeq / Sasa / HBondCount) ────────────────


def _load_scalar_welford_channels(grp, value_name: str) -> dict:
    """Read the 5 standard scalar Welford channels keyed on `value_name`.

    Per the C++ TR convention: `<value>`, `<value>_delta`,
    `<value>_abs_delta`, `<value>_delta_squared`, `<value>_dxdt`.
    """
    return {
        value_name: _read_moments(grp, value_name),
        f"{value_name}_delta": _read_moments(grp, f"{value_name}_delta"),
        f"{value_name}_abs_delta": _read_moments(grp, f"{value_name}_abs_delta"),
        f"{value_name}_delta_squared": _read_moments(grp, f"{value_name}_delta_squared"),
        f"{value_name}_dxdt": _read_moments(grp, f"{value_name}_dxdt"),
    }


@dataclass(frozen=True)
class EeqWelfordGroup:
    """Per-atom EEq charge Welford rollup from /trajectory/eeq_welford/."""
    charge: WelfordMoments
    charge_delta: WelfordMoments
    charge_abs_delta: WelfordMoments
    charge_delta_squared: WelfordMoments
    charge_dxdt: WelfordMoments
    rms_delta: np.ndarray
    n_frames_per_atom: np.ndarray
    delta_n_per_atom: np.ndarray
    mean_dt_ps: float
    frame_index_range: tuple[int, int]
    units: str


def _load_eeq_welford(f) -> Optional[EeqWelfordGroup]:
    if "/trajectory/eeq_welford" not in f:
        return None
    grp = f["/trajectory/eeq_welford"]
    chans = _load_scalar_welford_channels(grp, "charge")
    return EeqWelfordGroup(
        rms_delta=grp["rms_delta"][:],
        n_frames_per_atom=grp["n_frames_per_atom"][:],
        delta_n_per_atom=grp["delta_n_per_atom"][:],
        mean_dt_ps=_group_mean_dt_ps(grp),
        frame_index_range=_group_frame_index_range(grp),
        units=_group_units(grp),
        **chans,
    )


@dataclass(frozen=True)
class SasaWelfordGroup:
    """Per-atom SASA Welford rollup from /trajectory/sasa_welford/."""
    sasa: WelfordMoments
    sasa_delta: WelfordMoments
    sasa_abs_delta: WelfordMoments
    sasa_delta_squared: WelfordMoments
    sasa_dxdt: WelfordMoments
    rms_delta: np.ndarray
    n_frames_per_atom: np.ndarray
    delta_n_per_atom: np.ndarray
    mean_dt_ps: float
    frame_index_range: tuple[int, int]
    units: str


def _load_sasa_welford(f) -> Optional[SasaWelfordGroup]:
    if "/trajectory/sasa_welford" not in f:
        return None
    grp = f["/trajectory/sasa_welford"]
    chans = _load_scalar_welford_channels(grp, "sasa")
    return SasaWelfordGroup(
        rms_delta=grp["rms_delta"][:],
        n_frames_per_atom=grp["n_frames_per_atom"][:],
        delta_n_per_atom=grp["delta_n_per_atom"][:],
        mean_dt_ps=_group_mean_dt_ps(grp),
        frame_index_range=_group_frame_index_range(grp),
        units=_group_units(grp),
        **chans,
    )


@dataclass(frozen=True)
class HBondCountWelfordGroup:
    """Per-atom H-bond count Welford rollup from /trajectory/hbond_count_welford/.

    Adds `occupancy_fraction` channel (dimensionless ∈ [0,1]) on top of
    the standard scalar shape, plus the `source_radius_A` group
    attribute (default 3.5 Å, the donor-acceptor counting cutoff).
    """
    count: WelfordMoments
    count_delta: WelfordMoments
    count_abs_delta: WelfordMoments
    count_delta_squared: WelfordMoments
    count_dxdt: WelfordMoments
    occupancy_fraction: WelfordMoments
    rms_delta: np.ndarray
    n_frames_per_atom: np.ndarray
    delta_n_per_atom: np.ndarray
    mean_dt_ps: float
    frame_index_range: tuple[int, int]
    units: str
    source_radius_A: float


def _load_hbond_count_welford(f) -> Optional[HBondCountWelfordGroup]:
    if "/trajectory/hbond_count_welford" not in f:
        return None
    grp = f["/trajectory/hbond_count_welford"]
    chans = _load_scalar_welford_channels(grp, "count")
    return HBondCountWelfordGroup(
        occupancy_fraction=_read_moments(grp, "occupancy_fraction"),
        rms_delta=grp["rms_delta"][:],
        n_frames_per_atom=grp["n_frames_per_atom"][:],
        delta_n_per_atom=grp["delta_n_per_atom"][:],
        mean_dt_ps=_group_mean_dt_ps(grp),
        frame_index_range=_group_frame_index_range(grp),
        units=_group_units(grp),
        source_radius_A=float(grp.attrs.get("source_radius_A", 3.5)),
        **chans,
    )


# ─── Container for all six Welford groups ───────────────────────────


@dataclass(frozen=True)
class WelfordAccess:
    """Optional Welford H5 groups attached to TrajectoryData.

    Each field is None when the corresponding *WelfordTrajectoryResult
    was not attached during the C++ extraction run that produced this
    trajectory.h5. Field set is fixed; new Welford TRs require an SDK
    update (Export Everything Upstream — no implicit consumer).
    """
    bs: Optional[BsWelfordGroup] = None
    hm: Optional[HmWelfordGroup] = None
    mc: Optional[McConnellWelfordGroup] = None
    eeq: Optional[EeqWelfordGroup] = None
    sasa: Optional[SasaWelfordGroup] = None
    hbond_count: Optional[HBondCountWelfordGroup] = None


# ─── TrajectoryData and load_trajectory ─────────────────────────────


@dataclass
class TrajectoryData:
    """All trajectory-level data for one protein from the H5 master file.

    This is the training data for the waveform model, the analysis
    data for the ridge, and the movie data for VTK.
    """
    protein_id: str
    n_atoms: int
    n_frames: int

    # Per-frame positions: (T, N, 3) for VTK movies + waveform input
    positions: np.ndarray       # (T, N, 3) float64
    frame_times: np.ndarray     # (T,) float64

    # Per-atom rollup statistics (legacy GromacsProtein columns)
    rollup: TrajectoryRollup

    # Per-bond rollup statistics
    bonds: Optional[BondRollup]

    # Welford H5 groups (Phase 2b/C, 2026-05-17/18). Always present;
    # individual fields are None when the corresponding TR didn't run.
    welford: WelfordAccess = field(default_factory=WelfordAccess)


def load_trajectory(path: str | Path) -> TrajectoryData:
    """Load a trajectory H5 master file.

    The file is written by GromacsProtein::WriteH5 after a trajectory scan,
    optionally with /trajectory/<kind>_welford/ groups attached by the six
    *WelfordTrajectoryResult subclasses (2026-05-17/18).
    """
    import h5py

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Trajectory H5 not found: {path}")

    with h5py.File(path, "r") as f:
        protein_id = f.attrs.get("protein_id", path.stem)
        n_atoms = int(f.attrs["n_atoms"])
        n_frames = int(f.attrs["n_frames"])

        # Positions: stored as (T*N, 3), reshape to (T, N, 3)
        pos_raw = f["positions"][:]          # (T*N, 3)
        T_stored = int(f.attrs.get("positions_shape_T", n_frames))
        N_stored = int(f.attrs.get("positions_shape_N", n_atoms))
        positions = pos_raw.reshape(T_stored, N_stored, 3)

        frame_times = f["frame_times"][:]    # (T,)

        # Legacy GromacsProtein rollup
        means = f["rollup/mean"][:]          # (N, K)
        stds = f["rollup/std"][:]            # (N, K)
        names_raw = f["rollup/names"][:]
        names = [n.decode() if isinstance(n, bytes) else n for n in names_raw]

        # Bonds (optional)
        bonds = None
        if "bonds" in f:
            bonds = BondRollup(
                atom_a=f["bonds/atom_a"][:],
                atom_b=f["bonds/atom_b"][:],
                length_mean=f["bonds/length_mean"][:],
                length_std=f["bonds/length_std"][:],
            )

        # Welford H5 groups (each Optional — missing group → None)
        welford = WelfordAccess(
            bs=_load_bs_welford(f),
            hm=_load_hm_welford(f),
            mc=_load_mc_welford(f),
            eeq=_load_eeq_welford(f),
            sasa=_load_sasa_welford(f),
            hbond_count=_load_hbond_count_welford(f),
        )

    return TrajectoryData(
        protein_id=protein_id if isinstance(protein_id, str) else protein_id.decode(),
        n_atoms=n_atoms,
        n_frames=n_frames,
        positions=positions,
        frame_times=frame_times,
        rollup=TrajectoryRollup(names, means, stds),
        bonds=bonds,
        welford=welford,
    )
