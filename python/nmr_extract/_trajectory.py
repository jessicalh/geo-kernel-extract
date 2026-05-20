"""Trajectory data reader — loads either H5 schema produced by the C++ extractor.

    from nmr_extract import load_trajectory
    traj = load_trajectory("output/trajectory.h5")

    traj.positions          # (T, N, 3) float64 — normalized regardless of source schema
    traj.frame_times        # (T,) float64 — time in ps
    traj.n_atoms            # int
    traj.n_frames           # int

Schema tolerance: the reader auto-detects which schema the H5 file uses
and adapts. Codex 2026-05-18 caught the production-output incompatibility:

  - **Analysis schema (current `TrajectoryProtein::WriteH5`)**: root
    attrs `{protein_id, n_atoms, finalized}` only. Frame metadata at
    `/trajectory/frames/{time_ps, original_index}` (+ `n_frames` group
    attr). Positions atom-major at `/trajectory/positions/xyz`. No
    `/rollup` group — analysis mode replaces the legacy rollup with the
    six Welford H5 groups below. `traj.rollup is None`.
  - **Legacy ensemble (`GromacsProtein::WriteH5`)**: root attrs include
    `n_frames`, `positions_shape_T/N`. Positions frame-major at
    `/positions`. Frame times at `/frame_times`. Rollup at `/rollup/`.
    `traj.rollup` is a `TrajectoryRollup`.

Both schemas: positions are returned as `(T, N, 3)` after the reader
normalizes the on-disk layout.

    # Legacy rollup (Optional — only present in ensemble H5 files)
    if traj.rollup is not None:
        traj.rollup.bs_T0.mean  # (N,) ring current isotropic mean

    # Bonds: per-bond length statistics (legacy ensemble path only)
    if traj.bonds is not None:
        traj.bonds.length_mean

    # Welford H5 groups — written by *WelfordTrajectoryResult subclasses
    # (BS / HM / McConnell / Eeq / Sasa / HBondCount), 2026-05-17/18.
    # Each group is Optional — None when the corresponding C++ TR was
    # not attached for the run that produced this trajectory.h5.
    traj.welford.bs.t0.mean                  # (N,) ppm_T_per_nA
    traj.welford.bs.t0.units                 # "ppm_T_per_nA" (sample-channel)
    traj.welford.bs.t0.m2_units              # "ppm_T^2_per_nA^2" (Welford M2 accumulator)
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
    max_frame}` plus their per-dataset `units` attributes. Per-component
    channels (T1, T2) carry (N, K) arrays; scalar channels carry (N,).

    Units carry two distinct strings:
    - `units`: base unit of the sample (mean/std/min/max). For `t0` this
      is e.g. `"ppm_T_per_nA"`; for `charge` it's `"elementary_charge"`.
    - `m2_units`: unit of the Welford M2 accumulator (sum of (sample-mean)²),
      which has squared dimension. For `t0_m2` it's `"ppm_T^2_per_nA^2"`;
      for `charge_delta_squared_m2` it's `"elementary_charge^4"` (the
      sample is already squared, so M2 is base⁴). The C++ writer goes to
      explicit effort to emit honest per-dataset units; the SDK preserves
      both so downstream calibration math doesn't have to guess.

    Per-dataset units for `min_frame` / `max_frame` are always
    `"frame_index"` and not stored separately.
    """
    mean: np.ndarray
    m2: np.ndarray
    std: np.ndarray
    min: np.ndarray
    max: np.ndarray
    min_frame: np.ndarray   # uint frame index
    max_frame: np.ndarray   # uint frame index
    units: str              # base unit string from <prefix>_mean.attrs
    m2_units: str           # squared-units string from <prefix>_m2.attrs


def _read_moments(grp, prefix: str) -> WelfordMoments:
    """Read a 7-stat Welford block from an HDF5 group.

    Reads `<prefix>_{mean, m2, std, min, max, min_frame, max_frame}` and
    propagates the per-dataset `units` attribute from BOTH `<prefix>_mean`
    (→ `WelfordMoments.units`) and `<prefix>_m2` (→
    `WelfordMoments.m2_units`). C++ writes distinct unit strings for
    the squared accumulator vs the value channels (codex MEDIUM finding
    2026-05-18); collapsing to one would lose honesty. Works for both
    1D scalar channels and 2D per-component channels — numpy carries
    the shape through.
    """
    ds_mean = grp[f"{prefix}_mean"]
    ds_m2 = grp[f"{prefix}_m2"]
    units = ds_mean.attrs.get("units", "")
    m2_units = ds_m2.attrs.get("units", "")
    if isinstance(units, bytes):
        units = units.decode()
    if isinstance(m2_units, bytes):
        m2_units = m2_units.decode()
    return WelfordMoments(
        mean=ds_mean[:],
        m2=ds_m2[:],
        std=grp[f"{prefix}_std"][:],
        min=grp[f"{prefix}_min"][:],
        max=grp[f"{prefix}_max"][:],
        min_frame=grp[f"{prefix}_min_frame"][:],
        max_frame=grp[f"{prefix}_max_frame"][:],
        units=str(units),
        m2_units=str(m2_units),
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
    n_frames_per_atom: np.ndarray    # (N,) count of t0 samples
    delta_n_per_atom: np.ndarray     # (N,) count of t0_delta/abs/sq samples
    # Count of VALID-dt t0_dxdt samples (only frames with dt > MIN_DT_PS
    # contribute). Per codex 2026-05-18: distinct from delta_n because
    # zero-dt frames are skipped in the rate accumulator rather than
    # zero-filled. A well-formed trajectory has dxdt_n == delta_n on every
    # atom; mismatch flags frame-duplication or stride misconfiguration.
    dxdt_n_per_atom: np.ndarray      # (N,) count of t0_dxdt samples
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
        dxdt_n_per_atom=_read_dxdt_n_or_fallback(grp),
        mean_dt_ps=_group_mean_dt_ps(grp),
        frame_index_range=_group_frame_index_range(grp),
        irrep_layout_t1=str(_decode_attr(
            grp.attrs.get("irrep_layout_t1", ""))),
        irrep_layout_t2=str(_decode_attr(
            grp.attrs.get("irrep_layout_t2", ""))),
        units=_group_units(grp),
    )


def _read_dxdt_n_or_fallback(grp) -> np.ndarray:
    """Read `dxdt_n_per_atom` if present, else fall back to `delta_n_per_atom`.

    The `dxdt_n_per_atom` dataset is new (codex 2026-05-18 fix: zero-dt
    frames are skipped in the rate accumulator, so dxdt has its own
    counter). Older H5 files from before the fix don't have it; for those,
    return `delta_n_per_atom` so downstream consumers don't crash. The
    fallback IS the wrong number for files that had zero-dt frames, but
    those files also have the zero-fill bias the fix addresses — using
    delta_n_per_atom as a proxy is the least-broken default.
    """
    if "dxdt_n_per_atom" in grp:
        return grp["dxdt_n_per_atom"][:]
    return grp["delta_n_per_atom"][:]


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
    # Count of VALID-dt *_dxdt samples (codex 2026-05-18). Distinct from
    # delta_n_per_atom because zero-dt frames are skipped in the rate
    # accumulator rather than zero-filled.
    dxdt_n_per_atom: np.ndarray
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
        dxdt_n_per_atom=_read_dxdt_n_or_fallback(grp),
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
    # Count of VALID-dt *_dxdt samples (codex 2026-05-18). Distinct from
    # delta_n_per_atom because zero-dt frames are skipped in the rate
    # accumulator rather than zero-filled.
    dxdt_n_per_atom: np.ndarray
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
        dxdt_n_per_atom=_read_dxdt_n_or_fallback(grp),
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
    # Count of VALID-dt *_dxdt samples (codex 2026-05-18). Distinct from
    # delta_n_per_atom because zero-dt frames are skipped in the rate
    # accumulator rather than zero-filled.
    dxdt_n_per_atom: np.ndarray
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
        dxdt_n_per_atom=_read_dxdt_n_or_fallback(grp),
        mean_dt_ps=_group_mean_dt_ps(grp),
        frame_index_range=_group_frame_index_range(grp),
        units=_group_units(grp),
        source_radius_A=float(grp.attrs.get("source_radius_A", 3.5)),
        **chans,
    )


# ─── Energy time-series groups (system-scalar + per-atom) ──────────


@dataclass(frozen=True)
class GromacsEnergyTimeSeriesGroup:
    """Per-frame system-scalar timeline from /trajectory/gromacs_energy_time_series/.

    Shape: per-frame system scalars (T,). The two tensor channels (virial,
    pressure_tensor) carry their 3×3 layout as the second axis, (T, 9),
    with layout order XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ matching the source
    GromacsEnergy struct.

    Primary low-energy-state-filtering channel: `total_energy`. Selection
    of bottom-N% lowest-energy frames MUST mask NaN entries first — both
    source-absent frames and missing .edr columns populate NaN, and a
    naive `np.argsort(e)[:N]` is unsafe when the bottom-N% slice runs
    into the NaN region (numpy sorts NaN to the end, so the failure is
    silent: you get fewer-than-N valid frames or all-NaN selections).
    Safe pattern (R3 codex F1 2026-05-18):

        e = traj.energy.gromacs.total_energy
        valid = np.isfinite(e)
        n_select = int(0.1 * valid.sum())
        rank_in_valid = np.argsort(e[valid])[:n_select]
        selected = np.where(valid)[0][rank_in_valid]

    Provenance fields exposed (added 2026-05-18 R2 review): source-attached
    mask + count for the GromacsEnergyResult source-attached gate;
    `energy_frame_times_ps` carries the .edr time-stamp matched by
    EnergyAtTime() per frame — `energy_frame_times_ps - frame_times`
    is the snap distance (zero on exact match, nonzero when the .edr
    stride differs from the trajectory stride).
    """
    # Electrostatic (kJ/mol)
    coulomb_sr: np.ndarray       # (T,)
    coulomb_recip: np.ndarray
    coulomb_14: np.ndarray
    # Bonded (kJ/mol)
    bond: np.ndarray
    angle: np.ndarray
    urey_bradley: np.ndarray
    proper_dih: np.ndarray
    improper_dih: np.ndarray
    cmap_dih: np.ndarray
    # Van der Waals (kJ/mol)
    lj_sr: np.ndarray
    lj_14: np.ndarray
    disper_corr: np.ndarray
    # Thermodynamic state
    potential: np.ndarray        # kJ/mol
    kinetic: np.ndarray
    total_energy: np.ndarray
    enthalpy: np.ndarray
    temperature: np.ndarray      # K
    pressure: np.ndarray         # bar
    volume: np.ndarray           # nm^3
    density: np.ndarray          # kg/m^3
    # Box dimensions (nm)
    box_x: np.ndarray
    box_y: np.ndarray
    box_z: np.ndarray
    # Per-group temperatures (K)
    T_protein: np.ndarray
    T_non_protein: np.ndarray
    # 3×3 tensors as (T, 9)
    virial: np.ndarray           # (T, 9) kJ/mol
    pressure_tensor: np.ndarray  # (T, 9) bar
    # Frame indexing
    frame_indices: np.ndarray    # (T,) uint64
    frame_times: np.ndarray      # (T,) ps
    energy_frame_times_ps: np.ndarray  # (T,) ps — matched .edr row time
                                       # (NaN on source-absent frames)
    # Source-attached gate provenance (R2 review 2026-05-18)
    source_attached_per_frame: np.ndarray  # (T,) uint8 — 1=attached, 0=absent
    source_attached_count: int
    units: str                   # primary energy unit string


def _load_gromacs_energy_time_series(f) -> Optional[GromacsEnergyTimeSeriesGroup]:
    path = "/trajectory/gromacs_energy_time_series"
    if path not in f:
        return None
    g = f[path]
    return GromacsEnergyTimeSeriesGroup(
        coulomb_sr=g["coulomb_sr"][:],
        coulomb_recip=g["coulomb_recip"][:],
        coulomb_14=g["coulomb_14"][:],
        bond=g["bond"][:],
        angle=g["angle"][:],
        urey_bradley=g["urey_bradley"][:],
        proper_dih=g["proper_dih"][:],
        improper_dih=g["improper_dih"][:],
        cmap_dih=g["cmap_dih"][:],
        lj_sr=g["lj_sr"][:],
        lj_14=g["lj_14"][:],
        disper_corr=g["disper_corr"][:],
        potential=g["potential"][:],
        kinetic=g["kinetic"][:],
        total_energy=g["total_energy"][:],
        enthalpy=g["enthalpy"][:],
        temperature=g["temperature"][:],
        pressure=g["pressure"][:],
        volume=g["volume"][:],
        density=g["density"][:],
        box_x=g["box_x"][:],
        box_y=g["box_y"][:],
        box_z=g["box_z"][:],
        T_protein=g["T_protein"][:],
        T_non_protein=g["T_non_protein"][:],
        virial=g["virial"][:],
        pressure_tensor=g["pressure_tensor"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        energy_frame_times_ps=g["energy_frame_times_ps"][:]
            if "energy_frame_times_ps" in g else np.full(
                g["frame_times"].shape, np.nan, dtype=np.float64),
        source_attached_per_frame=g["source_attached_per_frame"][:]
            if "source_attached_per_frame" in g else np.ones(
                g["frame_times"].shape, dtype=np.uint8),
        source_attached_count=int(g.attrs.get(
            "source_attached_count", len(g["frame_times"]))),
        units=_group_units(g),
    )


@dataclass(frozen=True)
class BondedEnergyTimeSeriesGroup:
    """Per-atom 7-channel bonded-energy breakdown from /trajectory/bonded_energy_time_series/.

    Shape: each channel is (N, T) float64 in kJ/mol. The `total` channel
    is the running sum of the six interaction-type channels per atom per
    frame; emitted alongside per Export-Everything-Upstream.

    Split convention: GROMACS CHARMM36m interaction energies are split
    evenly among the 2..5 atoms participating in each interaction.
    """
    bond: np.ndarray             # (N, T) kJ/mol
    angle: np.ndarray
    urey_bradley: np.ndarray
    proper_dih: np.ndarray
    improper_dih: np.ndarray
    cmap_dih: np.ndarray         # name matches GromacsEnergy.cmap_dih
    total: np.ndarray
    frame_indices: np.ndarray    # (T,) uint64
    frame_times: np.ndarray      # (T,) ps
    # Source-attached gate provenance (R2 review 2026-05-18)
    source_attached_per_frame: np.ndarray  # (T,) uint8
    source_attached_count: int
    units: str
    split_convention: str
    split_convention_note: str   # "one of several valid attributions; ..."
    system_scope: str            # "protein_slice_only; sum != .edr whole-system term"


def _load_bonded_energy_time_series(f) -> Optional[BondedEnergyTimeSeriesGroup]:
    path = "/trajectory/bonded_energy_time_series"
    if path not in f:
        return None
    g = f[path]
    return BondedEnergyTimeSeriesGroup(
        bond=g["bond"][:],
        angle=g["angle"][:],
        urey_bradley=g["urey_bradley"][:],
        proper_dih=g["proper_dih"][:],
        improper_dih=g["improper_dih"][:],
        cmap_dih=g["cmap_dih"][:],
        total=g["total"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:]
            if "source_attached_per_frame" in g else np.ones(
                g["frame_times"].shape, dtype=np.uint8),
        source_attached_count=int(g.attrs.get(
            "source_attached_count", len(g["frame_times"]))),
        units=_group_units(g),
        split_convention=str(_decode_attr(
            g.attrs.get("split_convention", ""))),
        split_convention_note=str(_decode_attr(
            g.attrs.get("split_convention_note", ""))),
        system_scope=str(_decode_attr(
            g.attrs.get("system_scope", ""))),
    )


@dataclass(frozen=True)
class WaterFieldTimeSeriesGroup:
    """Per-atom water E-field + EFG timeline from /trajectory/water_field_time_series/.

    Shape: each channel is (N, T, ...) where the trailing axis carries
    the per-channel layout:

      efield, efield_first    (N, T, 3)  float64  V/Angstrom
      efg, efg_first          (N, T, 5)  float64  V/Angstrom²
                                         (real-spherical-tesseral T2 only)
      n_first, n_second       (N, T)     uint32   shell-occupancy count

    Water EFG carries no T0 (trace-zero by construction) and no T1
    (antisymmetric pseudovector vanishes for symmetric r⊗r contributions).
    Only T2 is emitted; consumers must NOT reshape to 9 components.

    "Absent, not faked" provenance: `source_attached_per_frame` is a
    per-frame bool mask. Frames where WaterFieldResult was not attached
    (no-solvent extraction) are NaN-filled on the float channels and
    carry `absent_sentinel = 0xFFFFFFFFu` on the uint32 shell-count
    datasets. When zero frames had the source attached, the entire H5
    group is absent — `load_trajectory(...).water_field.time_series` is
    `None` in that case.
    """
    efield: np.ndarray              # (N, T, 3) V/Angstrom
    efield_first: np.ndarray        # (N, T, 3)
    efg: np.ndarray                 # (N, T, 5) V/Angstrom² (T2 only)
    efg_first: np.ndarray           # (N, T, 5)
    # Shell-occupancy counts as RAW uint32 — absent-source frames carry
    # the sentinel `count_absent_sentinel` (0xFFFFFFFF). A naive
    # `np.mean(n_first)` over the raw array will see 4.29e9; downstream
    # code MUST mask first via `count_absent_sentinel` or call the
    # `n_first_float` / `n_second_float` accessors below which return
    # float arrays with NaN on absent frames.
    n_first: np.ndarray             # (N, T) uint32 — shell-occupancy
    n_second: np.ndarray            # (N, T) uint32
    count_absent_sentinel: int      # 0xFFFFFFFF — uint32 doesn't carry NaN
    frame_indices: np.ndarray       # (T,)  uint64
    frame_times: np.ndarray         # (T,)  ps
    source_attached_per_frame: np.ndarray  # (T,) uint8 — 1=attached, 0=absent
    source_attached_count: int
    # Layout strings (from H5 attrs)
    efield_layout: str              # "x,y,z"
    efield_parity: str              # "1o"
    efg_irrep_layout: str           # "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
    efg_parity: str                 # "2e"
    efield_cutoff_A: float          # 15.0
    n_first_cutoff_A: float         # 3.5
    n_second_cutoff_A: float        # 5.5

    @property
    def n_first_float(self) -> np.ndarray:
        """`n_first` promoted to float64 with NaN on source-absent frames.

        Use this instead of the raw uint32 `n_first` when computing means
        or feeding to a model — the uint32 sentinel `0xFFFFFFFF` would
        otherwise poison aggregations (4.29e9 waters per atom per absent
        frame). R3 codex F3 2026-05-18.
        """
        out = self.n_first.astype(np.float64)
        out[self.n_first == self.count_absent_sentinel] = np.nan
        return out

    @property
    def n_second_float(self) -> np.ndarray:
        """`n_second` promoted to float64 with NaN on source-absent frames."""
        out = self.n_second.astype(np.float64)
        out[self.n_second == self.count_absent_sentinel] = np.nan
        return out


def _load_water_field_time_series(f) -> Optional[WaterFieldTimeSeriesGroup]:
    path = "/trajectory/water_field_time_series"
    if path not in f:
        return None
    g = f[path]
    # absent_sentinel comes from the n_first dataset attr (set by C++);
    # fall back to the documented 0xFFFFFFFF if attr missing on older H5.
    if "absent_sentinel" in g["n_first"].attrs:
        sentinel = int(g["n_first"].attrs["absent_sentinel"])
    else:
        sentinel = 0xFFFFFFFF
    return WaterFieldTimeSeriesGroup(
        efield=g["efield"][:],
        efield_first=g["efield_first"][:],
        efg=g["efg"][:],
        efg_first=g["efg_first"][:],
        n_first=g["n_first"][:],
        n_second=g["n_second"][:],
        count_absent_sentinel=sentinel,
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        efield_layout=str(_decode_attr(g.attrs.get("efield_layout", ""))),
        efield_parity=str(_decode_attr(g.attrs.get("efield_parity", ""))),
        efg_irrep_layout=str(_decode_attr(g.attrs.get("efg_irrep_layout", ""))),
        efg_parity=str(_decode_attr(g.attrs.get("efg_parity", ""))),
        efield_cutoff_A=float(g.attrs.get("efield_cutoff_A", float("nan"))),
        n_first_cutoff_A=float(g.attrs.get("n_first_cutoff_A", float("nan"))),
        n_second_cutoff_A=float(g.attrs.get("n_second_cutoff_A", float("nan"))),
    )


@dataclass(frozen=True)
class WaterFieldWelfordGroup:
    """Per-atom Welford rollup of water E-field + EFG from /trajectory/water_field_welford/.

    Channels:
      efield_{x,y,z}, efield_magnitude              (N,)
      efield_first_{x,y,z}, efield_first_magnitude  (N,)
      efg_t2[5], efg_t2magnitude                    (N, 5) and (N,)
      efg_first_t2[5], efg_first_t2magnitude        (N, 5) and (N,)
      n_first, n_second                             (N,)
      Delta variants on 3 primary scalars: efield_magnitude, n_first, n_second
        (signed / abs / squared / dxdt / rms_delta)

    Like the TimeSeries TR, EFG has T0+T1 structural zeros — only T2.
    `source_attached_per_frame` mask is emitted alongside; Welford updates
    skip source-absent frames so the mean isn't biased toward zero.
    """
    # Per-component E-field
    efield_x: WelfordMoments
    efield_y: WelfordMoments
    efield_z: WelfordMoments
    efield_magnitude: WelfordMoments
    efield_first_x: WelfordMoments
    efield_first_y: WelfordMoments
    efield_first_z: WelfordMoments
    efield_first_magnitude: WelfordMoments
    # EFG T2 per-component + |T2|
    efg_t2: WelfordMoments          # (N, 5) per-component
    efg_t2magnitude: WelfordMoments
    efg_first_t2: WelfordMoments
    efg_first_t2magnitude: WelfordMoments
    # Shell occupancy counts
    n_first: WelfordMoments
    n_second: WelfordMoments
    # Delta variants on the 3 primary scalars
    efield_magnitude_delta: WelfordMoments
    efield_magnitude_abs_delta: WelfordMoments
    efield_magnitude_delta_squared: WelfordMoments
    efield_magnitude_dxdt: WelfordMoments
    efield_magnitude_rms_delta: np.ndarray   # (N,)
    n_first_delta: WelfordMoments
    n_first_abs_delta: WelfordMoments
    n_first_delta_squared: WelfordMoments
    n_first_dxdt: WelfordMoments
    n_first_rms_delta: np.ndarray
    n_second_delta: WelfordMoments
    n_second_abs_delta: WelfordMoments
    n_second_delta_squared: WelfordMoments
    n_second_dxdt: WelfordMoments
    n_second_rms_delta: np.ndarray
    # Provenance
    n_frames_per_atom: np.ndarray   # (N,)
    delta_n_per_atom: np.ndarray
    dxdt_n_per_atom: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    mean_dt_ps: float
    frame_index_range: tuple[int, int]


def _load_water_field_welford(f) -> Optional[WaterFieldWelfordGroup]:
    path = "/trajectory/water_field_welford"
    if path not in f:
        return None
    g = f[path]
    return WaterFieldWelfordGroup(
        efield_x=_read_moments(g, "efield_x"),
        efield_y=_read_moments(g, "efield_y"),
        efield_z=_read_moments(g, "efield_z"),
        efield_magnitude=_read_moments(g, "efield_magnitude"),
        efield_first_x=_read_moments(g, "efield_first_x"),
        efield_first_y=_read_moments(g, "efield_first_y"),
        efield_first_z=_read_moments(g, "efield_first_z"),
        efield_first_magnitude=_read_moments(g, "efield_first_magnitude"),
        efg_t2=_read_moments(g, "efg_t2"),
        efg_t2magnitude=_read_moments(g, "efg_t2magnitude"),
        efg_first_t2=_read_moments(g, "efg_first_t2"),
        efg_first_t2magnitude=_read_moments(g, "efg_first_t2magnitude"),
        n_first=_read_moments(g, "n_first"),
        n_second=_read_moments(g, "n_second"),
        efield_magnitude_delta=_read_moments(g, "efield_magnitude_delta"),
        efield_magnitude_abs_delta=_read_moments(g, "efield_magnitude_abs_delta"),
        efield_magnitude_delta_squared=_read_moments(g, "efield_magnitude_delta_squared"),
        efield_magnitude_dxdt=_read_moments(g, "efield_magnitude_dxdt"),
        efield_magnitude_rms_delta=g["efield_magnitude_rms_delta"][:],
        n_first_delta=_read_moments(g, "n_first_delta"),
        n_first_abs_delta=_read_moments(g, "n_first_abs_delta"),
        n_first_delta_squared=_read_moments(g, "n_first_delta_squared"),
        n_first_dxdt=_read_moments(g, "n_first_dxdt"),
        n_first_rms_delta=g["n_first_rms_delta"][:],
        n_second_delta=_read_moments(g, "n_second_delta"),
        n_second_abs_delta=_read_moments(g, "n_second_abs_delta"),
        n_second_delta_squared=_read_moments(g, "n_second_delta_squared"),
        n_second_dxdt=_read_moments(g, "n_second_dxdt"),
        n_second_rms_delta=g["n_second_rms_delta"][:],
        n_frames_per_atom=g["n_frames_per_atom"][:],
        delta_n_per_atom=g["delta_n_per_atom"][:],
        dxdt_n_per_atom=g["dxdt_n_per_atom"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        mean_dt_ps=_group_mean_dt_ps(g),
        frame_index_range=_group_frame_index_range(g),
    )


@dataclass(frozen=True)
class HydrationGeometryTimeSeriesGroup:
    """Per-atom SASA-normal water polarisation timeline from
    /trajectory/hydration_geometry_time_series/.

    Six channels: net water dipole vector + SASA outward normal + first-
    shell water count + the polarisation-signal trio (alignment,
    coherence, half_shell_asymmetry). The trio IS the polarisation
    signal — calibration weight on these channels is what we want to
    learn.

    Shape:
      dipole_vector, surface_normal  (N, T, 3)  float64
      first_shell_count               (N, T)    uint32 (sentinel
                                                 `0xFFFFFFFF` on
                                                 source-absent frames)
      half_shell_asymmetry, dipole_alignment, dipole_coherence
                                      (N, T)    float64 (NaN on absent)

    "Absent, not faked" provenance: `source_attached_per_frame` is a
    per-frame uint8 mask. When zero frames had the source attached, the
    H5 group is absent and `load_trajectory(...).hydration_geometry.
    time_series` is `None`.
    """
    dipole_vector: np.ndarray         # (N, T, 3) e·Å
    surface_normal: np.ndarray        # (N, T, 3) unit vector
    first_shell_count: np.ndarray     # (N, T) uint32 — sentinel on absent
    half_shell_asymmetry: np.ndarray  # (N, T)
    dipole_alignment: np.ndarray      # (N, T)  cos angle
    # dipole_coherence: source formula `|Σ d_i| / n_shell` (e·Å in
    # numerator, dimensionless in denominator) → e·Å, NOT a [0,1]
    # dimensionless order parameter. The H5 dataset attr `units` is
    # "e_Angstrom". R6 codex 2026-05-18.
    dipole_coherence: np.ndarray      # (N, T)  e·Å
    count_absent_sentinel: int        # 0xFFFFFFFF
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    dipole_vector_layout: str
    dipole_vector_parity: str
    surface_normal_layout: str
    surface_normal_parity: str
    # SASA-normal frame disambiguator — sibling HydrationShell TR uses
    # COM. Same dataset name (half_shell_asymmetry) means the same
    # physical question in different reference frames; never aggregate
    # cross-frame without checking this attr.
    reference_frame: str              # "SASA_normal"
    polarisation_signal_channels: str

    @property
    def first_shell_count_float(self) -> np.ndarray:
        """`first_shell_count` promoted to float64 with NaN on sentinel."""
        out = self.first_shell_count.astype(np.float64)
        out[self.first_shell_count == self.count_absent_sentinel] = np.nan
        return out


def _load_hydration_geometry_time_series(f) -> Optional[HydrationGeometryTimeSeriesGroup]:
    path = "/trajectory/hydration_geometry_time_series"
    if path not in f:
        return None
    g = f[path]
    sentinel = int(g["first_shell_count"].attrs.get(
        "absent_sentinel", 0xFFFFFFFF))
    return HydrationGeometryTimeSeriesGroup(
        dipole_vector=g["dipole_vector"][:],
        surface_normal=g["surface_normal"][:],
        first_shell_count=g["first_shell_count"][:],
        half_shell_asymmetry=g["half_shell_asymmetry"][:],
        dipole_alignment=g["dipole_alignment"][:],
        dipole_coherence=g["dipole_coherence"][:],
        count_absent_sentinel=sentinel,
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        dipole_vector_layout=str(_decode_attr(
            g.attrs.get("dipole_vector_layout", ""))),
        dipole_vector_parity=str(_decode_attr(
            g.attrs.get("dipole_vector_parity", ""))),
        surface_normal_layout=str(_decode_attr(
            g.attrs.get("surface_normal_layout", ""))),
        surface_normal_parity=str(_decode_attr(
            g.attrs.get("surface_normal_parity", ""))),
        reference_frame=str(_decode_attr(g.attrs.get("reference_frame", ""))),
        polarisation_signal_channels=str(_decode_attr(
            g.attrs.get("polarisation_signal_channels", ""))),
    )


@dataclass(frozen=True)
class HydrationGeometryWelfordGroup:
    """Per-atom Welford rollup of SASA-normal water polarisation features.

    Per-component Vec3 Welford on dipole_vector + surface_normal. Scalar
    Welford on the 4 polarisation scalars (asymmetry, alignment, coherence,
    shell_count) with signed/abs/sq/dxdt delta variants on each.
    """
    # Per-component Vec3
    dipole_vector_x: WelfordMoments
    dipole_vector_y: WelfordMoments
    dipole_vector_z: WelfordMoments
    dipole_magnitude: WelfordMoments
    surface_normal_x: WelfordMoments
    surface_normal_y: WelfordMoments
    surface_normal_z: WelfordMoments
    # Polarisation scalars
    half_shell_asymmetry: WelfordMoments
    dipole_alignment: WelfordMoments
    dipole_coherence: WelfordMoments
    shell_count: WelfordMoments
    # Delta variants on the 4 primary scalars
    half_shell_asymmetry_delta: WelfordMoments
    half_shell_asymmetry_abs_delta: WelfordMoments
    half_shell_asymmetry_delta_squared: WelfordMoments
    half_shell_asymmetry_dxdt: WelfordMoments
    half_shell_asymmetry_rms_delta: np.ndarray
    dipole_alignment_delta: WelfordMoments
    dipole_alignment_abs_delta: WelfordMoments
    dipole_alignment_delta_squared: WelfordMoments
    dipole_alignment_dxdt: WelfordMoments
    dipole_alignment_rms_delta: np.ndarray
    dipole_coherence_delta: WelfordMoments
    dipole_coherence_abs_delta: WelfordMoments
    dipole_coherence_delta_squared: WelfordMoments
    dipole_coherence_dxdt: WelfordMoments
    dipole_coherence_rms_delta: np.ndarray
    shell_count_delta: WelfordMoments
    shell_count_abs_delta: WelfordMoments
    shell_count_delta_squared: WelfordMoments
    shell_count_dxdt: WelfordMoments
    shell_count_rms_delta: np.ndarray
    # Provenance
    n_frames_per_atom: np.ndarray
    delta_n_per_atom: np.ndarray
    dxdt_n_per_atom: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    mean_dt_ps: float
    frame_index_range: tuple[int, int]
    # SASA-normal frame disambiguator (sibling HS Welford uses COM).
    reference_frame: str               # "SASA_normal"
    # Documents that 0.0 in dipole_alignment is bimodal: either a real
    # zero-cos-angle alignment or the |dipole_sum|<NEAR_ZERO_NORM
    # branch in HydrationGeometryResult.cpp where the alignment is
    # set to 0.0 for lack of a meaningful direction. Consumers reading
    # alignment_mean ≈ 0 should consult dipole_magnitude before
    # treating it as a polarisation signal.
    dipole_alignment_zero_sentinel: str


def _load_hydration_geometry_welford(f) -> Optional[HydrationGeometryWelfordGroup]:
    path = "/trajectory/hydration_geometry_welford"
    if path not in f:
        return None
    g = f[path]
    return HydrationGeometryWelfordGroup(
        dipole_vector_x=_read_moments(g, "dipole_vector_x"),
        dipole_vector_y=_read_moments(g, "dipole_vector_y"),
        dipole_vector_z=_read_moments(g, "dipole_vector_z"),
        dipole_magnitude=_read_moments(g, "dipole_magnitude"),
        surface_normal_x=_read_moments(g, "surface_normal_x"),
        surface_normal_y=_read_moments(g, "surface_normal_y"),
        surface_normal_z=_read_moments(g, "surface_normal_z"),
        half_shell_asymmetry=_read_moments(g, "half_shell_asymmetry"),
        dipole_alignment=_read_moments(g, "dipole_alignment"),
        dipole_coherence=_read_moments(g, "dipole_coherence"),
        shell_count=_read_moments(g, "shell_count"),
        half_shell_asymmetry_delta=_read_moments(g, "half_shell_asymmetry_delta"),
        half_shell_asymmetry_abs_delta=_read_moments(g, "half_shell_asymmetry_abs_delta"),
        half_shell_asymmetry_delta_squared=_read_moments(g, "half_shell_asymmetry_delta_squared"),
        half_shell_asymmetry_dxdt=_read_moments(g, "half_shell_asymmetry_dxdt"),
        half_shell_asymmetry_rms_delta=g["half_shell_asymmetry_rms_delta"][:],
        dipole_alignment_delta=_read_moments(g, "dipole_alignment_delta"),
        dipole_alignment_abs_delta=_read_moments(g, "dipole_alignment_abs_delta"),
        dipole_alignment_delta_squared=_read_moments(g, "dipole_alignment_delta_squared"),
        dipole_alignment_dxdt=_read_moments(g, "dipole_alignment_dxdt"),
        dipole_alignment_rms_delta=g["dipole_alignment_rms_delta"][:],
        dipole_coherence_delta=_read_moments(g, "dipole_coherence_delta"),
        dipole_coherence_abs_delta=_read_moments(g, "dipole_coherence_abs_delta"),
        dipole_coherence_delta_squared=_read_moments(g, "dipole_coherence_delta_squared"),
        dipole_coherence_dxdt=_read_moments(g, "dipole_coherence_dxdt"),
        dipole_coherence_rms_delta=g["dipole_coherence_rms_delta"][:],
        shell_count_delta=_read_moments(g, "shell_count_delta"),
        shell_count_abs_delta=_read_moments(g, "shell_count_abs_delta"),
        shell_count_delta_squared=_read_moments(g, "shell_count_delta_squared"),
        shell_count_dxdt=_read_moments(g, "shell_count_dxdt"),
        shell_count_rms_delta=g["shell_count_rms_delta"][:],
        n_frames_per_atom=g["n_frames_per_atom"][:],
        delta_n_per_atom=g["delta_n_per_atom"][:],
        dxdt_n_per_atom=g["dxdt_n_per_atom"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        mean_dt_ps=_group_mean_dt_ps(g),
        frame_index_range=_group_frame_index_range(g),
        reference_frame=str(_decode_attr(g.attrs.get("reference_frame", ""))),
        dipole_alignment_zero_sentinel=str(_decode_attr(
            g.attrs.get("dipole_alignment_zero_sentinel", ""))),
    )


@dataclass(frozen=True)
class HydrationGeometryAccess:
    """Hydration geometry TR family — TimeSeries + Welford rollup pair.
    Either slot is None when the corresponding C++ TR was not attached
    OR attached but had source_attached_count == 0 (no-solvent extraction).
    """
    time_series: Optional[HydrationGeometryTimeSeriesGroup] = None
    welford: Optional[HydrationGeometryWelfordGroup] = None


@dataclass(frozen=True)
class HydrationShellTimeSeriesGroup:
    """Per-atom COM-based hydration shell timeline from
    /trajectory/hydration_shell_time_series/.

    Four scalar channels (N, T):
      half_shell_asymmetry   COM-reference frame, fraction ∈ [0, 1]
      mean_water_dipole_cos  water orientation order parameter
      nearest_ion_distance   Angstrom; +inf sentinel when no ion within
                             ion_cutoff (default 20 Å). Use np.isfinite()
                             to filter buried atoms before aggregation.
      nearest_ion_charge     elementary charge; 0 when no ion in cutoff

    Sibling of HydrationGeometryTimeSeries (SASA-normal frame). Both
    methods accumulate per `feedback_methods_accumulate` — calibration
    weights them separately.
    """
    half_shell_asymmetry: np.ndarray   # (N, T)
    mean_water_dipole_cos: np.ndarray  # (N, T)
    nearest_ion_distance: np.ndarray   # (N, T) Å; +inf on no-ion-in-cutoff
    nearest_ion_charge: np.ndarray     # (N, T) e
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    reference_frame: str               # "COM"
    nearest_ion_distance_sentinel: str  # description string


def _load_hydration_shell_time_series(f) -> Optional[HydrationShellTimeSeriesGroup]:
    path = "/trajectory/hydration_shell_time_series"
    if path not in f:
        return None
    g = f[path]
    return HydrationShellTimeSeriesGroup(
        half_shell_asymmetry=g["half_shell_asymmetry"][:],
        mean_water_dipole_cos=g["mean_water_dipole_cos"][:],
        nearest_ion_distance=g["nearest_ion_distance"][:],
        nearest_ion_charge=g["nearest_ion_charge"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        reference_frame=str(_decode_attr(g.attrs.get("reference_frame", ""))),
        nearest_ion_distance_sentinel=str(_decode_attr(
            g.attrs.get("nearest_ion_distance_sentinel", ""))),
    )


@dataclass(frozen=True)
class HydrationShellWelfordGroup:
    """Per-atom Welford rollup of COM-based hydration shell features.

    `nearest_ion_distance` is on a CONDITIONAL Welford (R6 codex
    2026-05-18): only frames with a finite ion-in-cutoff distance
    contribute. Atoms with no ion observed in cutoff across any frame
    finalize to NaN (n=0). For the "is there an ion here?" question,
    use `ion_present_fraction.mean ∈ [0, 1]`. For the conditional
    distance distribution given there IS an ion, use
    `nearest_ion_distance.{mean, std, min, max}` directly — finite-only
    accumulation means no inf-poisoning.

    Per-atom provenance counters:
      n_ion_present_per_atom   divisor for `nearest_ion_distance` base
                               Welford (frames with finite distance).
      n_ion_delta_per_atom     divisor for delta variants of
                               `nearest_ion_distance` (consecutive
                               attached frames both finite).
      n_ion_dxdt_per_atom      divisor for `nearest_ion_distance_dxdt`
                               (n_ion_delta + dt > 0 gate).
    Other channels use delta_n_per_atom / dxdt_n_per_atom as before.
    """
    half_shell_asymmetry: WelfordMoments
    mean_water_dipole_cos: WelfordMoments
    nearest_ion_distance: WelfordMoments
    nearest_ion_charge: WelfordMoments
    ion_present_fraction: WelfordMoments
    # Delta variants on all 4 scalars
    half_shell_asymmetry_delta: WelfordMoments
    half_shell_asymmetry_abs_delta: WelfordMoments
    half_shell_asymmetry_delta_squared: WelfordMoments
    half_shell_asymmetry_dxdt: WelfordMoments
    half_shell_asymmetry_rms_delta: np.ndarray
    mean_water_dipole_cos_delta: WelfordMoments
    mean_water_dipole_cos_abs_delta: WelfordMoments
    mean_water_dipole_cos_delta_squared: WelfordMoments
    mean_water_dipole_cos_dxdt: WelfordMoments
    mean_water_dipole_cos_rms_delta: np.ndarray
    nearest_ion_distance_delta: WelfordMoments
    nearest_ion_distance_abs_delta: WelfordMoments
    nearest_ion_distance_delta_squared: WelfordMoments
    nearest_ion_distance_dxdt: WelfordMoments
    nearest_ion_distance_rms_delta: np.ndarray
    nearest_ion_charge_delta: WelfordMoments
    nearest_ion_charge_abs_delta: WelfordMoments
    nearest_ion_charge_delta_squared: WelfordMoments
    nearest_ion_charge_dxdt: WelfordMoments
    nearest_ion_charge_rms_delta: np.ndarray
    # Provenance — generic
    n_frames_per_atom: np.ndarray
    delta_n_per_atom: np.ndarray
    dxdt_n_per_atom: np.ndarray
    # Provenance — nearest_ion_distance conditional counters
    n_ion_present_per_atom: np.ndarray
    n_ion_delta_per_atom: np.ndarray
    n_ion_dxdt_per_atom: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    mean_dt_ps: float
    frame_index_range: tuple[int, int]
    reference_frame: str               # "COM"
    nearest_ion_distance_welford_policy: str


def _load_hydration_shell_welford(f) -> Optional[HydrationShellWelfordGroup]:
    path = "/trajectory/hydration_shell_welford"
    if path not in f:
        return None
    g = f[path]
    bases = ("half_shell_asymmetry", "mean_water_dipole_cos",
             "nearest_ion_distance", "nearest_ion_charge")

    def read_chs(base: str) -> dict:
        return {
            base: _read_moments(g, base),
            f"{base}_delta": _read_moments(g, f"{base}_delta"),
            f"{base}_abs_delta": _read_moments(g, f"{base}_abs_delta"),
            f"{base}_delta_squared": _read_moments(g, f"{base}_delta_squared"),
            f"{base}_dxdt": _read_moments(g, f"{base}_dxdt"),
            f"{base}_rms_delta": g[f"{base}_rms_delta"][:],
        }

    chs = {}
    for base in bases:
        chs.update(read_chs(base))

    return HydrationShellWelfordGroup(
        ion_present_fraction=_read_moments(g, "ion_present_fraction"),
        n_frames_per_atom=g["n_frames_per_atom"][:],
        delta_n_per_atom=g["delta_n_per_atom"][:],
        dxdt_n_per_atom=g["dxdt_n_per_atom"][:],
        n_ion_present_per_atom=g["n_ion_present_per_atom"][:],
        n_ion_delta_per_atom=g["n_ion_delta_per_atom"][:],
        n_ion_dxdt_per_atom=g["n_ion_dxdt_per_atom"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        mean_dt_ps=_group_mean_dt_ps(g),
        frame_index_range=_group_frame_index_range(g),
        reference_frame=str(_decode_attr(g.attrs.get("reference_frame", ""))),
        nearest_ion_distance_welford_policy=str(_decode_attr(
            g.attrs.get("nearest_ion_distance_welford_policy", ""))),
        **chs,
    )


@dataclass(frozen=True)
class HydrationShellAccess:
    """COM-based hydration shell TR family. Either slot is None when
    the C++ TR didn't run or had source_attached_count == 0."""
    time_series: Optional[HydrationShellTimeSeriesGroup] = None
    welford: Optional[HydrationShellWelfordGroup] = None


# ─── Dihedral time-series (per-residue, R rows × T frames) ──────────


@dataclass(frozen=True)
class DihedralTimeSeriesGroup:
    """Per-residue backbone + sidechain dihedral timelines from
    /trajectory/dihedral_time_series/. First per-residue (R, T) TR
    in the SDK (2026-05-19).

    Convention pins (also recorded as group attrs):
    - All angles in radians.
    - IUPAC signed dihedral atan2(y, x). Value range [-π, π] with
      discontinuity at ±π; consumers comparing across frames must
      use circular differences.
    - omega_deviation = WrapPi(omega - π), emitted for every well-
      defined peptide bond INCLUDING X→Pro (cis/trans isomerism at
      X-Pro is real signal, use `omega_is_xpro` mask to interpret).
    - rama_region uint8 legend: 0=unassigned, 1=αR, 2=β, 3=αL,
      4=PPII, 5=other. Resolution order αR → αL → PPII → β → other.
    - chi raw radians — PHE/TYR chi2 ring-flip + ASP/GLU carboxylate
      symmetry NOT applied; consumers needing mod-π apply themselves.
      TRP chi2 and HIS chi2 are NOT symmetric (CD1/CD2/ND1 chemically
      distinct).
    - Backbone connectivity is determined from the LegacyAmber bond
      graph (covalent C(i)-N(i+1) peptide bond), NOT from chain_id /
      sequence_number / terminal_state / insertion_code. Correct on
      antibody insertion-coded structures, engineered chimeras with
      non-monotonic numbering, and cyclic peptides.

    DSSP convention warning: `DsspResult.Phi/Psi` returned by the C++
    side via libdssp uses the NEGATED IUPAC convention. The values in
    THIS group use IUPAC directly (matching PlanarGeometryResult).
    Downstream code comparing this TR to DsspResult must negate one
    side. See `angle_convention` attr for the full note.

    Per-frame fields (R, T):
      phi              C(i-1)-N(i)-CA(i)-C(i); NaN at N-terminus and
                       non-covalent residue boundaries.
      psi              N(i)-CA(i)-C(i)-N(i+1); NaN at C-terminus and
                       non-covalent residue boundaries.
      omega            CA(i)-C(i)-N(i+1)-CA(i+1).
      omega_deviation  WrapPi(omega - π) ∈ [-π, π].
      chi              (R, T, 4) — NaN where chi[k] is not defined or
                       per-frame geometry is degenerate.
      rama_region      uint8 5-region literature grid (see legend).

    Static per-residue masks:
      chi_exists (R, 4) uint8 — chi[k] structurally cacheable; does
                                 NOT guarantee finite chi value at
                                 runtime, use isfinite(chi).
      omega_is_xpro / is_glycine / is_proline / is_pre_proline /
      residue_terminal_state (R,) uint8 — see source notes.

    Per-atom broadcast lookup:
      residue_index_per_atom (N,) int32 — atom_i → residue_i for
                                          SDK / viewer atom-axis
                                          broadcast (option C from
                                          design discussion).

    Frame metadata:
      frame_indices (T,) uint, frame_times (T,) float64 ps,
      source_attached_per_frame (T,) uint8 — trivially always 1
      (positions always present at tp.Seed); emitted for SDK
      uniformity with conditionally-attached-source TRs.
    """
    # Per-frame channels
    phi: np.ndarray              # (R, T) float64 radians
    psi: np.ndarray              # (R, T) float64 radians
    omega: np.ndarray            # (R, T) float64 radians
    omega_deviation: np.ndarray  # (R, T) float64 radians
    chi: np.ndarray              # (R, T, 4) float64 radians
    rama_region: np.ndarray      # (R, T) uint8
    # Static per-residue masks
    chi_exists: np.ndarray              # (R, 4) uint8
    omega_is_xpro: np.ndarray           # (R,) uint8
    is_glycine: np.ndarray              # (R,) uint8
    is_proline: np.ndarray              # (R,) uint8
    is_pre_proline: np.ndarray          # (R,) uint8
    residue_terminal_state: np.ndarray  # (R,) uint8
    chain_id_per_residue: np.ndarray    # (R,) variable-length strings
    # Per-atom lookup
    residue_index_per_atom: np.ndarray  # (N,) int32
    # Per-frame metadata
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    # Group provenance + convention pins
    angle_units: str            # "radians"
    periodicity: str            # "2pi"
    value_range: str            # "[-pi, pi] ..."
    angle_convention: str       # IUPAC convention + DSSP-divergence note
    chain_break_policy: str     # bond-graph derivation note
    omega_deviation_policy: str
    rama_region_legend: str
    rama_region_boundaries: str
    chi_symmetry_caveats: str
    residue_terminal_state_legend: str
    source_attached_policy: str  # "always_attached"
    chunking_policy: str


@dataclass(frozen=True)
class DihedralBinTransitionGroup:
    """Per-residue rotamer + Rama-region transition statistics from
    /trajectory/dihedral_bin_transition/. AV companion to
    DihedralTimeSeriesGroup (which carries per-frame raw angles).

    Bin labels match DihedralTimeSeriesGroup.rama_region verbatim:
      0=unassigned (phi or psi NaN at termini / non-bonded gaps),
      1=αR, 2=β, 3=αL, 4=PPII, 5=other. Boundaries follow
      Lovell-Richardson 2003 favored regions (Berkholz/Adzhubei tight
      PPII cone) — see DihedralTimeSeriesGroup.rama_region_boundaries
      attr for inclusive degree-coordinate ranges + resolution order.
    Chi rotamer bins (Lovell-Richardson 2003 convention, 2026-05-19
    fix):
      0=g+ (0 < chi < 120°), 1=trans (|chi| ≥ 120°), 2=g-
      (-120° < chi ≤ 0°), 255=unassigned (chi not defined for this AA
      or per-frame geometry degenerate). chi == 0 exactly lands in g-
      (rare in MD). Note: ChiRotamerSelectionTrajectoryResult (in the
      slated-for-removal ScanForDftPointSet config) uses a DIFFERENT
      strict-`>` convention that puts ±120° in gauche; this group
      uses the Lovell-Richardson trans-inclusive convention. The two
      are documented-divergent.

    Bin 0 (unassigned) IS populated in `backbone_bin_occupancy[:, 0]`
    (codex-review-2026-05-19 fix): every frame contributes to exactly
    one bin, so `sum(backbone_bin_occupancy[ri, :]) == n_frames` for
    all residues, and `n_frames_observed[ri] == sum(occupancy[ri, 1:])`.

    Transition gate: BOTH prev and curr frame must have an observed bin
    (non-unassigned) for a transition to count. Consecutive-frame walk;
    if intermediate frames are unobserved (NaN bin), the transition
    chain breaks and re-starts on the next observation.

    Stats per residue:
      backbone_transition_count   (R,)        uint32
      backbone_dominant_region    (R,)        uint8
      n_frames_observed           (R,)        uint32 (phi+psi both finite)
      backbone_bin_occupancy      (R, 6)      uint32 frame counts
      chi_transition_count        (R, 4)      uint32
      chi_dominant_rotamer        (R, 4)      uint8 (255 = no observation)
      chi_n_frames_observed       (R, 4)      uint32 (chi[k] finite)
      chi_rotamer_occupancy       (R, 4, 3)   uint32 frame counts

    Per-frame metadata:
      frame_indices (T,), frame_times (T,) ps,
      source_attached_per_frame (T,) uint8 (trivially all-1; positions
      always present so this TR has no conditional source).

    Backbone connectivity uses the canonical Protein::BackbonePredecessor
    / BackboneSuccessor walk (PeptideCN bond-graph filter) — same as
    DihedralTimeSeriesGroup. Cyclic / antibody insertion-coded
    structures get correct phi/psi via the wrap edge; the bin classification
    follows from the dihedral computation.
    """
    backbone_transition_count: np.ndarray  # (R,) uint32
    backbone_dominant_region: np.ndarray   # (R,) uint8
    n_frames_observed: np.ndarray          # (R,) uint32
    backbone_bin_occupancy: np.ndarray     # (R, 6) uint32
    chi_transition_count: np.ndarray       # (R, 4) uint32
    chi_dominant_rotamer: np.ndarray       # (R, 4) uint8
    chi_n_frames_observed: np.ndarray      # (R, 4) uint32
    chi_rotamer_occupancy: np.ndarray      # (R, 4, 3) uint32
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    # Convention pins
    backbone_bin_legend: str
    backbone_bin_boundaries: str
    chi_rotamer_legend: str
    transition_gate: str
    angle_convention: str
    source_attached_policy: str


@dataclass(frozen=True)
class Dssp8TimeSeriesGroup:
    """Per-residue per-frame DSSP 8-state SS + H-bond partners from
    /trajectory/dssp8_time_series/. Conditional-source TR: DsspResult
    attaches only when PerFrameRunOptions::skip_dssp is false.

    SS code: 0=H (alpha), 1=G (3_10), 2=I (pi), 3=E (extended),
             4=B (beta bridge), 5=T (turn), 6=S (bend), 7=C (coil);
             255 = no observation.
    H-bond partner: residue index of the partner, or -1 if no partner
                     at that slot. DSSP records up to 2 acceptors and
                     2 donors per residue.
    H-bond energy: kcal/mol; NaN if no partner at that slot.

    DSSP convention quirk: DsspResult.Phi/Psi (NOT in this group; see
    DihedralTimeSeriesGroup) returns NEGATED IUPAC phi/psi. This group
    does not mirror phi/psi/SASA from DSSP — those live in
    DihedralTimeSeriesGroup (IUPAC convention) and a future Sasa TR.
    """
    ss8_code: np.ndarray                   # (R, T) uint8
    hbond_acceptor_partner: np.ndarray     # (R, T, 2) int32
    hbond_acceptor_energy: np.ndarray      # (R, T, 2) float64
    hbond_donor_partner: np.ndarray        # (R, T, 2) int32
    hbond_donor_energy: np.ndarray         # (R, T, 2) float64
    residue_index_per_atom: np.ndarray     # (N,) int32
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    ss8_legend: str
    ss8_unassigned_sentinel: int           # 255
    hbond_partner_sentinel: str
    hbond_energy_units: str                # "kcal/mol"
    source: str
    source_attached_policy: str            # "conditional"


def _load_dssp8_time_series(f) -> Optional[Dssp8TimeSeriesGroup]:
    path = "/trajectory/dssp8_time_series"
    if path not in f:
        return None
    g = f[path]
    def _attr(name: str) -> str:
        return str(_decode_attr(g.attrs.get(name, "")))
    return Dssp8TimeSeriesGroup(
        ss8_code=g["ss8_code"][:],
        hbond_acceptor_partner=g["hbond_acceptor_partner"][:],
        hbond_acceptor_energy=g["hbond_acceptor_energy"][:],
        hbond_donor_partner=g["hbond_donor_partner"][:],
        hbond_donor_energy=g["hbond_donor_energy"][:],
        residue_index_per_atom=g["residue_index_per_atom"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        ss8_legend=_attr("ss8_legend"),
        ss8_unassigned_sentinel=int(g.attrs.get("ss8_unassigned_sentinel", 255)),
        hbond_partner_sentinel=_attr("hbond_partner_sentinel"),
        hbond_energy_units=_attr("hbond_energy_units"),
        source=_attr("source"),
        source_attached_policy=_attr("source_attached_policy"),
    )


@dataclass(frozen=True)
class Dssp8TransitionGroup:
    """Per-residue DSSP 8-state transition statistics from
    /trajectory/dssp8_transition/. AV companion to Dssp8TimeSeriesGroup.

    Stats per residue:
      ss8_transition_count   (R,) uint32
      ss8_dominant           (R,) uint8 (255 if no observation)
      n_frames_observed      (R,) uint32

    Per-residue per state (R, 8):
      ss8_occupancy          uint32 frame count per state

    Per-residue transition matrix (R, 8, 8):
      ss8_transition_matrix  uint32, M[ri, prev, curr] = count of
        consecutive observed-pair transitions (prev -> curr). Diagonal
        is identically zero (self-transitions not counted).

    Transition gate: BOTH prev and curr observed for the consecutive
    pair to count. Source-absent frames break the transition chain.
    """
    ss8_transition_count: np.ndarray       # (R,) uint32
    ss8_dominant: np.ndarray               # (R,) uint8
    n_frames_observed: np.ndarray          # (R,) uint32
    ss8_occupancy: np.ndarray              # (R, 8) uint32
    ss8_transition_matrix: np.ndarray      # (R, 8, 8) uint32
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    ss8_legend: str
    transition_matrix_layout: str
    transition_gate: str
    source: str
    source_attached_policy: str


def _load_dssp8_transition(f) -> Optional[Dssp8TransitionGroup]:
    path = "/trajectory/dssp8_transition"
    if path not in f:
        return None
    g = f[path]
    def _attr(name: str) -> str:
        return str(_decode_attr(g.attrs.get(name, "")))
    return Dssp8TransitionGroup(
        ss8_transition_count=g["ss8_transition_count"][:],
        ss8_dominant=g["ss8_dominant"][:],
        n_frames_observed=g["n_frames_observed"][:],
        ss8_occupancy=g["ss8_occupancy"][:],
        ss8_transition_matrix=g["ss8_transition_matrix"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        ss8_legend=_attr("ss8_legend"),
        transition_matrix_layout=_attr("transition_matrix_layout"),
        transition_gate=_attr("transition_gate"),
        source=_attr("source"),
        source_attached_policy=_attr("source_attached_policy"),
    )


@dataclass(frozen=True)
class Dssp8Access:
    """DSSP 8-state TR family. Either slot is None when the C++ TR
    didn't run for the extraction (skip_dssp=true) or had
    source_attached_count == 0."""
    time_series: Optional[Dssp8TimeSeriesGroup] = None
    transitions: Optional[Dssp8TransitionGroup] = None


@dataclass(frozen=True)
class RingPuckerTimeSeriesGroup:
    """Per-ring per-frame Cremer-Pople saturated-ring pucker (Q, θ) +
    per-aromatic-ring χ₂ from /trajectory/ring_pucker_time_series/.

    Two axes:
      pucker_Q     (S, T) float64 Å — saturated rings only (Pro
                                       pyrrolidine + any future saturated)
      pucker_theta (S, T) float64 degrees [0, 360) — saturated rings
      aromatic_chi2 (A, T) float64 radians — aromatic rings

    Per-ring static metadata for the SDK / viewer broadcast:
      saturated_parent_residue_index (S,) int32
      aromatic_parent_residue_index  (A,) int32

    Source: PlanarGeometryResult (conditionally attached when
    LegacyAmber substrate is populated). NaN within attached frames
    indicates per-ring degenerate geometry (incomplete ring or
    collinear vertices).

    Pucker theta mod 72° gives envelope vs twist endo/exo classification
    for 5-rings (Cremer-Pople 1975).
    Aromatic chi2: matches DihedralTimeSeriesGroup.chi[1] of the parent
    residue (PHE/TYR ring-flip canonical observable per Akke & Weininger
    2023). Per-frame value is instantaneous — flip kinetics not
    measurable from one frame.

    Slot fields are None when the ring count on that axis is zero
    (e.g., a protein with no aromatic residues would have A=0 and no
    aromatic_chi2 dataset emitted; the field is None here).
    """
    pucker_Q: Optional[np.ndarray]              # (S, T) or None
    pucker_theta: Optional[np.ndarray]          # (S, T) or None
    aromatic_chi2: Optional[np.ndarray]         # (A, T) or None
    saturated_parent_residue_index: np.ndarray  # (S,) int32 (may be empty)
    aromatic_parent_residue_index: np.ndarray   # (A,) int32 (may be empty)
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    source_attached_count: int
    n_saturated_rings: int
    n_aromatic_rings: int
    pucker_convention: str
    pucker_Q_units: str                         # "Angstrom"
    pucker_theta_units: str                     # "degrees"
    aromatic_chi2_units: str                    # "radians"
    aromatic_chi2_convention: str
    source: str
    source_attached_policy: str


@dataclass(frozen=True)
class JCouplingTimeSeriesGroup:
    """Per-residue per-frame Karplus ³J observables from
    /trajectory/j_coupling_time_series/.

    Nine datasets across eight channel families (R, T) float64 in Hz;
    thin Karplus transform of the persisted phi/chi1 timeline. NaN
    within frames indicates the channel is structurally absent for this
    residue (PRO has no HN; GLY has no Cβ; GLY/ALA have no chi1;
    SER/CYS/THR have non-carbon chi1 terminal so J(N,Cγ) / J(C',Cγ)
    NaN; N-terminus has no C'(prev) so J(Hα,C') NaN; ALA's methyl Cβ
    is deliberately excluded from J(Hα,Hβ)).

      J_HN_Halpha          ³J(HN, Hα) via H-N-CA-HA dihedral; phi
                           observable. Vuister & Bax 1993 JACS 115:7772
                           (DOI 10.1021/ja00070a024).
      J_HN_Halpha_Vogeli   Same atomic dihedral, alternate Karplus
                           parametrization. Vögeli, Ying, Grishaev &
                           Bax 2007 JACS 129:9377 (DOI 10.1021/
                           ja070324o), Table 1 "rigid" row. Methods-
                           accumulate alternate (both channels stay).
      J_HN_Cbeta           ³J(HN, Cβ) via H-N-CA-CB dihedral; orthogonal
                           phi observable. Wang & Bax 1996 JACS 118:2483
                           NMR/X-ray refined fit row 3 (DOI 10.1021/
                           ja9535524).
      J_HN_Cprime          ³J(HN, C') via H-N-CA-C dihedral; phi
                           observable. Wang & Bax 1996 Table 1 ROW 4
                           (θ=0°, NMR/X-ray refined fit, A=4.32,
                           B=+0.84, C=0.00). B is POSITIVE; J can be
                           slightly negative (~-0.04 Hz minimum).
                           Row-mapping fixed 2026-05-20 per codex F1:
                           prior bundle attributed row 2 values
                           (3.75, +2.19, 1.28) to this channel.
      J_Halpha_Cprime      ³J(Hα, C') via HA-CA-N-C'(prev) dihedral;
                           phi observable. The 3-bond path crosses
                           the peptide bond at the previous residue's
                           C'; rotation is around the N-CA axis (phi
                           axis), per Vuister teaching lecture sect
                           6.1 + Vogeli 2007 page 9384. Wang & Bax
                           1996 Table 1 ROW 2 (θ=-60°, A=3.75,
                           B=+2.19, C=1.28). B is POSITIVE. NaN at
                           N-terminus (no C'(prev)). Atom path +
                           row-mapping fixed 2026-05-20 per codex F1
                           + F2: prior bundle used HA-CA-C-N(next)
                           (psi axis, NOT phi) with row 4 values.
      J_N_Cgamma           ³J(N, Cγ) via N-CA-CB-CG (= chi1); chi1
                           rotamer. Pérez et al. 2001 JACS 123:7081
                           (DOI 10.1021/ja003724j), Table 2 consensus
                           row.
      J_Cprime_Cgamma      ³J(C', Cγ) via C-CA-CB-CG; chi1 rotamer
                           with C' leading. Pérez 2001 Table 2
                           consensus row, A=2.31, B=-0.87, C=0.55
                           (byte-verified 2026-05-19 — earlier circulated
                           value (1.74, -0.57, 0.25) was wrong).
      J_Halpha_Hbeta2      ³J(Hα, Hβ2) via HA-CA-CB-HB2; chi1 rotamer.
                           Pérez 2001 Table 2 consensus row. Methylene
                           pro-R/pro-S assignment follows IUPAC HB2/HB3
                           atom-name convention. Ile/Val/Thr methine
                           Hβ emitted in this slot.
      J_Halpha_Hbeta3      As J_Halpha_Hbeta2 but for HB3; for the
                           prochiral methylene pair. Ile/Val/Thr methine
                           Hβ mirrored in this slot (same value).
                           Gly/Ala: NaN (no methylene Hβ).

    Karplus form (post-codex-F6 + project-sign repair, 2026-05-20):
      Backbone channels (HN-Hα, HN-Cβ, HN-C', Hα-C'): J = A·cos²(φ +
      θ_offset) + B·cos(φ + θ_offset) + C, where φ is canonical
      Ramachandran C(prev)-N-CA-C (IUPAC signed atan2, radians), and
      θ_offset is stored in the project phi convention (opposite sign
      from Wang-Bax / Vögeli printed phi):
        HN-Hα: θ_offset = +π/3
        HN-Cβ: θ_offset = -π/3
        HN-C': θ_offset = 0
        Hα-C': θ_offset = +π/3
      Chi1 channels (J(N,Cγ), J(C',Cγ), J(Hα,Hβ{2,3})): J = A·cos²(α)
      + B·cos(α) + C, where α is the actual 4-atom atomic dihedral
      computed directly from positions (Pérez 2001 Table 2 footnote c
      form -- the per-coupling (A, B, C) internalize the Cα-substituent
      offset).

      This is the same sign convention used by DihedralTimeSeries
      (`phi_DSSP = -phi_IUPAC`). The LiteratureAnchoredProbeOn1UBQ
      test is the executable guard for the project-sign mapping.

    Coefficients per channel (byte-verified 2026-05-19):
      J_HN_Halpha          A=6.51,  B=-1.76, C=1.60 (range [1.48, 9.87])
      J_HN_Halpha_Vogeli   A=7.97,  B=-1.26, C=0.63 (range [0.58, 9.86])
      J_HN_Cbeta           A=3.39,  B=-0.94, C=0.07 (range [0.005, 4.40])
      J_HN_Cprime          A=4.32,  B=+0.84, C=0.00 (range [-0.04, 5.16])
                                    (Wang-Bax row 4; fixed 2026-05-20)
      J_Halpha_Cprime      A=3.75,  B=+2.19, C=1.28 (range [0.96, 7.22])
                                    (Wang-Bax row 2; fixed 2026-05-20)
      J_N_Cgamma           A=1.29,  B=-0.49, C=0.37 (range [0.32, 2.15])
      J_Cprime_Cgamma      A=2.31,  B=-0.87, C=0.55 (range [0.47, 3.73])
      J_Halpha_Hbeta{2,3}  A=7.23,  B=-1.37, C=2.22 (range [2.16, 10.82])

    All numerics live in src/PhysicalConstants.h and are byte-verified
    against the source PDFs in references/.

    Static per-residue masks (R,) uint8:
      J_HN_Halpha_exists      1 if C(prev) + H + N + CA + C + HA all
                              cached (gates both J_HN_Halpha and
                              J_HN_Halpha_Vogeli).
      J_HN_Cbeta_exists       1 if C(prev) + H + N + CA + C + CB all
                              cached (PRO=0, GLY=0).
      J_HN_Cprime_exists      1 if C(prev) + H + N + CA + C all cached
                              (PRO=0).
      J_Halpha_Cprime_exists  1 if C(prev) + N + CA + C + HA all
                              cached. The physical atom path is
                              HA-CA-N-C'(prev), but the emitted value
                              is phi-derived and therefore also needs
                              current-residue C. N-terminus=0; codex
                              F2/F4 2026-05-20 flipped this from
                              C-terminus.
      J_chi1_exists           1 if chi1 atoms valid (residue has chi1
                              defined; GLY/ALA → 0). Necessary but
                              NOT sufficient for J_N_Cgamma /
                              J_Cprime_Cgamma (those further require
                              Element::C at chi[0].a[3]).
      J_N_Cgamma_exists       1 if J_chi1_exists AND chi[0].a[3]
                              element is Carbon (SER/CYS/THR → 0
                              because chi1 terminal is OG/SG/OG1).
      J_Cprime_Cgamma_exists  1 if J_N_Cgamma_exists AND C/CA cached.
      J_Halpha_Hbeta_exists   1 if HA + CA + CB + (HB2 or HB3 or HB
                              methine) cached (GLY=0, ALA=0).

    Per-atom lookup (N,) int32:
      residue_index_per_atom  atom_i → residue_i broadcast for the SDK /
                              viewer (atom-axis access to the per-residue
                              J channels).

    GLY caveat: GLY uses Residue.HA which is HA2 by Residue.h convention;
    HA3 is not separately measured here. The Vuister-Bax 1993 fit DID
    include glycine ³J values, so the published (A, B, C) absorb the
    pro-R/pro-S averaging error. ALA Hβ channels are NaN by design
    (methyl ≠ methylene observable). Consumers needing strict pro-R/
    pro-S resolution should compute directly from the per-atom indices.

    Source: positions + Residue backbone-cache (H, N, CA, HA, CB, C) +
    chi1 atom indices + per-residue Hβ atoms (looked up by IUPAC name)
    + C'(prev) (resolved via Protein::BackbonePredecessor bond-graph
    query, NOT ri-1 / chain_id-equality which is a banned adjacency
    anti-pattern). No source ConformationResult dependency (positions
    present from tp.Seed; source_attached_per_frame trivially all-1
    for SDK uniformity under the OBJECT_MODEL "Conditional-attach TR"
    canonical statement).
    """
    J_HN_Halpha: np.ndarray                  # (R, T) Hz
    J_HN_Halpha_Vogeli: np.ndarray           # (R, T) Hz
    J_HN_Cbeta: np.ndarray                   # (R, T) Hz
    J_HN_Cprime: np.ndarray                  # (R, T) Hz
    J_Halpha_Cprime: np.ndarray              # (R, T) Hz
    J_N_Cgamma: np.ndarray                   # (R, T) Hz
    J_Cprime_Cgamma: np.ndarray              # (R, T) Hz
    J_Halpha_Hbeta2: np.ndarray              # (R, T) Hz
    J_Halpha_Hbeta3: np.ndarray              # (R, T) Hz
    J_HN_Halpha_exists: np.ndarray           # (R,) uint8
    J_HN_Cbeta_exists: np.ndarray            # (R,) uint8
    J_HN_Cprime_exists: np.ndarray           # (R,) uint8
    J_Halpha_Cprime_exists: np.ndarray       # (R,) uint8
    J_chi1_exists: np.ndarray                # (R,) uint8
    J_N_Cgamma_exists: np.ndarray            # (R,) uint8
    J_Cprime_Cgamma_exists: np.ndarray       # (R,) uint8
    J_Halpha_Hbeta_exists: np.ndarray        # (R,) uint8
    residue_index_per_atom: np.ndarray       # (N,) int32
    frame_indices: np.ndarray
    frame_times: np.ndarray
    source_attached_per_frame: np.ndarray
    n_residues: int
    n_atoms: int
    n_frames: int
    karplus_form: str
    J_HN_Halpha_coefficients: str
    J_HN_Halpha_Vogeli_coefficients: str
    J_HN_Cbeta_coefficients: str
    J_HN_Cprime_coefficients: str
    J_Halpha_Cprime_coefficients: str
    J_N_Cgamma_coefficients: str
    J_Cprime_Cgamma_coefficients: str
    J_Halpha_Hbeta_coefficients: str
    dihedral_convention: str
    GLY_caveat: str
    units: str                                # "Hz"
    absent_sentinel: str                      # "NaN"
    source: str
    source_attached_policy: str


def _load_j_coupling_time_series(f) -> Optional[JCouplingTimeSeriesGroup]:
    path = "/trajectory/j_coupling_time_series"
    if path not in f:
        return None
    g = f[path]
    def _attr(name: str) -> str:
        return str(_decode_attr(g.attrs.get(name, "")))
    return JCouplingTimeSeriesGroup(
        J_HN_Halpha=g["J_HN_Halpha"][:],
        J_HN_Halpha_Vogeli=g["J_HN_Halpha_Vogeli"][:],
        J_HN_Cbeta=g["J_HN_Cbeta"][:],
        J_HN_Cprime=g["J_HN_Cprime"][:],
        J_Halpha_Cprime=g["J_Halpha_Cprime"][:],
        J_N_Cgamma=g["J_N_Cgamma"][:],
        J_Cprime_Cgamma=g["J_Cprime_Cgamma"][:],
        J_Halpha_Hbeta2=g["J_Halpha_Hbeta2"][:],
        J_Halpha_Hbeta3=g["J_Halpha_Hbeta3"][:],
        J_HN_Halpha_exists=g["J_HN_Halpha_exists"][:],
        J_HN_Cbeta_exists=g["J_HN_Cbeta_exists"][:],
        J_HN_Cprime_exists=g["J_HN_Cprime_exists"][:],
        J_Halpha_Cprime_exists=g["J_Halpha_Cprime_exists"][:],
        J_chi1_exists=g["J_chi1_exists"][:],
        J_N_Cgamma_exists=g["J_N_Cgamma_exists"][:],
        J_Cprime_Cgamma_exists=g["J_Cprime_Cgamma_exists"][:],
        J_Halpha_Hbeta_exists=g["J_Halpha_Hbeta_exists"][:],
        residue_index_per_atom=g["residue_index_per_atom"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        n_residues=int(g.attrs["n_residues"]),
        n_atoms=int(g.attrs["n_atoms"]),
        n_frames=int(g.attrs["n_frames"]),
        karplus_form=_attr("karplus_form"),
        J_HN_Halpha_coefficients=_attr("J_HN_Halpha_coefficients"),
        J_HN_Halpha_Vogeli_coefficients=_attr("J_HN_Halpha_Vogeli_coefficients"),
        J_HN_Cbeta_coefficients=_attr("J_HN_Cbeta_coefficients"),
        J_HN_Cprime_coefficients=_attr("J_HN_Cprime_coefficients"),
        J_Halpha_Cprime_coefficients=_attr("J_Halpha_Cprime_coefficients"),
        J_N_Cgamma_coefficients=_attr("J_N_Cgamma_coefficients"),
        J_Cprime_Cgamma_coefficients=_attr("J_Cprime_Cgamma_coefficients"),
        J_Halpha_Hbeta_coefficients=_attr("J_Halpha_Hbeta_coefficients"),
        dihedral_convention=_attr("dihedral_convention"),
        GLY_caveat=_attr("GLY_caveat"),
        units=_attr("units"),
        absent_sentinel=_attr("absent_sentinel"),
        source=_attr("source"),
        source_attached_policy=_attr("source_attached_policy"),
    )


def _load_ring_pucker_time_series(f) -> Optional[RingPuckerTimeSeriesGroup]:
    path = "/trajectory/ring_pucker_time_series"
    if path not in f:
        return None
    g = f[path]
    def _attr(name: str) -> str:
        return str(_decode_attr(g.attrs.get(name, "")))
    S = int(g.attrs["n_saturated_rings"])
    A = int(g.attrs["n_aromatic_rings"])
    return RingPuckerTimeSeriesGroup(
        pucker_Q=g["pucker_Q"][:] if "pucker_Q" in g else None,
        pucker_theta=g["pucker_theta"][:] if "pucker_theta" in g else None,
        aromatic_chi2=g["aromatic_chi2"][:] if "aromatic_chi2" in g else None,
        saturated_parent_residue_index=g["saturated_parent_residue_index"][:],
        aromatic_parent_residue_index=g["aromatic_parent_residue_index"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        source_attached_count=int(g.attrs["source_attached_count"]),
        n_saturated_rings=S,
        n_aromatic_rings=A,
        pucker_convention=_attr("pucker_convention"),
        pucker_Q_units=_attr("pucker_Q_units"),
        pucker_theta_units=_attr("pucker_theta_units"),
        aromatic_chi2_units=_attr("aromatic_chi2_units"),
        aromatic_chi2_convention=_attr("aromatic_chi2_convention"),
        source=_attr("source"),
        source_attached_policy=_attr("source_attached_policy"),
    )


def _load_dihedral_bin_transition(f) -> Optional[DihedralBinTransitionGroup]:
    path = "/trajectory/dihedral_bin_transition"
    if path not in f:
        return None
    g = f[path]
    def _attr(name: str) -> str:
        return str(_decode_attr(g.attrs.get(name, "")))
    return DihedralBinTransitionGroup(
        backbone_transition_count=g["backbone_transition_count"][:],
        backbone_dominant_region=g["backbone_dominant_region"][:],
        n_frames_observed=g["n_frames_observed"][:],
        backbone_bin_occupancy=g["backbone_bin_occupancy"][:],
        chi_transition_count=g["chi_transition_count"][:],
        chi_dominant_rotamer=g["chi_dominant_rotamer"][:],
        chi_n_frames_observed=g["chi_n_frames_observed"][:],
        chi_rotamer_occupancy=g["chi_rotamer_occupancy"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        backbone_bin_legend=_attr("backbone_bin_legend"),
        backbone_bin_boundaries=_attr("backbone_bin_boundaries"),
        chi_rotamer_legend=_attr("chi_rotamer_legend"),
        transition_gate=_attr("transition_gate"),
        angle_convention=_attr("angle_convention"),
        source_attached_policy=_attr("source_attached_policy"),
    )


def _load_dihedral_time_series(f) -> Optional[DihedralTimeSeriesGroup]:
    path = "/trajectory/dihedral_time_series"
    if path not in f:
        return None
    g = f[path]

    def _attr(name: str) -> str:
        return str(_decode_attr(g.attrs.get(name, "")))

    # chain_id_per_residue is variable-length strings — h5py returns
    # bytes; decode to str for ergonomics.
    raw_chain_ids = g["chain_id_per_residue"][:]
    chain_ids = np.array(
        [c.decode("utf-8", errors="replace") if isinstance(c, (bytes, bytearray))
         else str(c) for c in raw_chain_ids],
        dtype=object)

    return DihedralTimeSeriesGroup(
        phi=g["phi"][:],
        psi=g["psi"][:],
        omega=g["omega"][:],
        omega_deviation=g["omega_deviation"][:],
        chi=g["chi"][:],
        rama_region=g["rama_region"][:],
        chi_exists=g["chi_exists"][:],
        omega_is_xpro=g["omega_is_xpro"][:],
        is_glycine=g["is_glycine"][:],
        is_proline=g["is_proline"][:],
        is_pre_proline=g["is_pre_proline"][:],
        residue_terminal_state=g["residue_terminal_state"][:],
        chain_id_per_residue=chain_ids,
        residue_index_per_atom=g["residue_index_per_atom"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        source_attached_per_frame=g["source_attached_per_frame"][:],
        angle_units=_attr("angle_units"),
        periodicity=_attr("periodicity"),
        value_range=_attr("value_range"),
        angle_convention=_attr("angle_convention"),
        chain_break_policy=_attr("chain_break_policy"),
        omega_deviation_policy=_attr("omega_deviation_policy"),
        rama_region_legend=_attr("rama_region_legend"),
        rama_region_boundaries=_attr("rama_region_boundaries"),
        chi_symmetry_caveats=_attr("chi_symmetry_caveats"),
        residue_terminal_state_legend=_attr("residue_terminal_state_legend"),
        source_attached_policy=_attr("source_attached_policy"),
        chunking_policy=_attr("chunking_policy"),
    )


@dataclass(frozen=True)
class WaterFieldAccess:
    """Water-field TR family — TimeSeries + Welford rollup pair.

    Either slot is None when the corresponding C++ TR was not attached
    (or attached but source-absent in 0/T frames). Mirrors the energy
    pair's EnergyAccess shape.
    """
    time_series: Optional[WaterFieldTimeSeriesGroup] = None
    welford: Optional[WaterFieldWelfordGroup] = None


@dataclass(frozen=True)
class EnergyAccess:
    """Energy time-series groups attached to TrajectoryData.

    Either slot is None when the corresponding C++ TR didn't run for
    the extraction that produced this trajectory.h5. The `gromacs` slot
    is the gate for low-energy-state ML model framings; the `bonded`
    slot adds per-atom strain features for the same frames.
    """
    gromacs: Optional[GromacsEnergyTimeSeriesGroup] = None
    bonded: Optional[BondedEnergyTimeSeriesGroup] = None


# ─── Container for all six Welford groups ───────────────────────────


@dataclass(frozen=True)
class WelfordAccess:
    """Kernel-class Welford H5 groups (BS / HM / McConnell / Eeq / Sasa /
    HBondCount). These were the first Welford rollup TRs landed (Phase 2b,
    2026-05-17) — kernel-level statistics with a uniform shape.

    Each field is None when the corresponding *WelfordTrajectoryResult
    was not attached during the C++ extraction run that produced this
    trajectory.h5. NEW Welford TRs that pair with a TimeSeries sibling
    route through their source-family Access container instead (e.g.
    `traj.water_field.welford`, `traj.hydration_geometry.welford`,
    `traj.hydration_shell.welford`) — those are domain-namespaced rather
    than added to this kernel-Welford grab-bag. Pattern picked up
    2026-05-18 with the F7 water SDK landing; the kernel-Welford set in
    this dataclass is fixed and historical.
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

    Two H5 schemas exist:
    - **Analysis (current)** written by `TrajectoryProtein::WriteH5` +
      per-TR `WriteH5Group(file)`. Root attrs `{protein_id, n_atoms,
      finalized}`. Frame metadata under `/trajectory/frames/`. Positions
      under `/trajectory/positions/xyz` (atom-major). Each attached TR
      writes its own `/trajectory/<group>/`. No `/rollup` group.
    - **Legacy ensemble** written by `GromacsProtein::WriteH5` (old
      ensemble accumulation path). Root attrs include `n_frames` /
      `positions_shape_{T,N}`. Positions at `/positions` (frame-major).
      Frame times at `/frame_times`. Legacy rollup at `/rollup/`.

    `load_trajectory` detects which schema is present and adapts; the
    fields below normalize on the analysis schema's representation
    (positions in `(T, N, 3)` regardless of source). Fields absent in
    a given schema are `None` or empty arrays — check `traj.rollup is
    not None` before using the legacy rollup.
    """
    protein_id: str
    n_atoms: int
    n_frames: int

    # Per-frame positions: (T, N, 3) for VTK movies + waveform input
    positions: np.ndarray       # (T, N, 3) float64
    frame_times: np.ndarray     # (T,) float64

    # Per-atom rollup statistics (legacy GromacsProtein columns). None
    # for analysis-schema H5 files — those replace the rollup notion
    # with the Welford H5 groups below.
    rollup: Optional[TrajectoryRollup] = None

    # Per-bond rollup statistics (legacy ensemble path only).
    bonds: Optional[BondRollup] = None

    # Welford H5 groups (Phase 2b/C, 2026-05-17/18). Always present;
    # individual fields are None when the corresponding TR didn't run.
    welford: WelfordAccess = field(default_factory=WelfordAccess)

    # Per-frame timeline groups (2026-05-18 batch onward). Always present
    # as a container; individual fields are None when the corresponding
    # TR didn't run. Unblocks the (a) all-frames and (b) low-energy-state
    # ML model framings on the 676-trajectory fleet — see
    # `project_fleet_676_2026-05-18`.
    energy: EnergyAccess = field(default_factory=EnergyAccess)

    # Water field TR pair (TimeSeries + Welford rollup).
    water_field: WaterFieldAccess = field(default_factory=WaterFieldAccess)

    # Hydration geometry TR pair (SASA-normal water polarisation).
    hydration_geometry: HydrationGeometryAccess = field(
        default_factory=HydrationGeometryAccess)

    # Hydration shell TR pair (COM-based water shell features).
    hydration_shell: HydrationShellAccess = field(
        default_factory=HydrationShellAccess)

    # Per-residue dihedral timeline (first (R, T) TR — 2026-05-19).
    # None when the C++ TR didn't run for the extraction that produced
    # this trajectory.h5. Movie-target consumer: broadcasts to atom axis
    # via `dihedrals.residue_index_per_atom` at render time.
    dihedrals: Optional[DihedralTimeSeriesGroup] = None

    # Per-residue rotamer + Rama-region transition statistics (AV
    # companion to dihedrals; 2026-05-19).
    dihedral_bin_transitions: Optional[DihedralBinTransitionGroup] = None

    # DSSP 8-state TR family (TS + transition pair; 2026-05-19).
    dssp8: Dssp8Access = field(default_factory=Dssp8Access)

    # Per-ring (saturated + aromatic) pucker timeline (2026-05-19).
    ring_pucker: Optional[RingPuckerTimeSeriesGroup] = None

    # Per-residue Karplus ³J observables (2026-05-19).
    j_coupling: Optional[JCouplingTimeSeriesGroup] = None


def _read_frame_metadata(f, n_atoms_hint: int):
    """Return (n_frames, frame_times) from either H5 schema.

    Analysis schema (current C++ TrajectoryProtein): `/trajectory/frames/`
    group with `n_frames` group-attr + `time_ps` dataset.
    Legacy schema (GromacsProtein ensemble): root attr `n_frames` +
    `/frame_times` dataset.
    Returns `(0, empty array)` when neither is present (e.g., bare H5
    fixture used for testing one specific group in isolation).
    """
    if "trajectory" in f and "frames" in f["trajectory"]:
        frames_grp = f["/trajectory/frames"]
        n_frames = int(frames_grp.attrs.get("n_frames", 0))
        if "time_ps" in frames_grp:
            return n_frames, frames_grp["time_ps"][:]
        return n_frames, np.array([], dtype=np.float64)
    if "n_frames" in f.attrs and "frame_times" in f:
        return int(f.attrs["n_frames"]), f["frame_times"][:]
    return 0, np.array([], dtype=np.float64)


def _read_positions(f, n_atoms_hint: int, n_frames_hint: int) -> np.ndarray:
    """Return (T, N, 3) positions from either H5 schema.

    Analysis schema: `/trajectory/positions/xyz` written atom-major as
    `(N*T, 3)` by `PositionsTimeSeriesTrajectoryResult::WriteH5Group`.
    Legacy schema: `/positions` written frame-major as `(T*N, 3)` by
    `GromacsProtein::WriteH5`. Both normalized to `(T, N, 3)` here.
    """
    if "trajectory" in f and "positions" in f["trajectory"] \
            and "xyz" in f["/trajectory/positions"]:
        pos_grp = f["/trajectory/positions"]
        N = int(pos_grp.attrs.get("n_atoms", n_atoms_hint))
        T = int(pos_grp.attrs.get("n_frames", n_frames_hint))
        pos_raw = pos_grp["xyz"][:]
        # Atom-major (N*T, 3) → (N, T, 3) → (T, N, 3)
        return pos_raw.reshape(N, T, 3).transpose(1, 0, 2)
    if "positions" in f:
        pos_raw = f["positions"][:]
        T = int(f.attrs.get("positions_shape_T", n_frames_hint))
        N = int(f.attrs.get("positions_shape_N", n_atoms_hint))
        return pos_raw.reshape(T, N, 3)
    return np.empty((0, n_atoms_hint, 3), dtype=np.float64)


def _read_legacy_rollup(f) -> Optional[TrajectoryRollup]:
    """Return legacy `/rollup/` TrajectoryRollup or None.

    Only present in GromacsProtein ensemble H5 files. Analysis H5 files
    written by the current TrajectoryProtein writer replace this notion
    with the per-TR Welford groups; for those, return None and let the
    caller use `traj.welford.<kind>` instead.
    """
    if "rollup" not in f or "mean" not in f["rollup"]:
        return None
    means = f["rollup/mean"][:]
    stds = f["rollup/std"][:]
    names_raw = f["rollup/names"][:]
    names = [_decode_attr(n) for n in names_raw]
    return TrajectoryRollup(names, means, stds)


def load_trajectory(path: str | Path) -> TrajectoryData:
    """Load a trajectory H5 master file.

    Supports both H5 schemas (see `TrajectoryData` docstring):
    - Analysis (current `TrajectoryProtein` + per-TR `WriteH5Group`)
    - Legacy ensemble (`GromacsProtein::WriteH5`)

    Schema detection is automatic; fields absent in one schema are
    `None` or empty in the returned `TrajectoryData`. The six Welford
    H5 groups are loaded independently of which schema's frame
    metadata is present — they live at `/trajectory/<kind>_welford/`
    in both.
    """
    import h5py

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Trajectory H5 not found: {path}")

    with h5py.File(path, "r") as f:
        protein_id = _decode_attr(f.attrs.get("protein_id", path.stem))
        n_atoms = int(f.attrs["n_atoms"])

        n_frames, frame_times = _read_frame_metadata(f, n_atoms)
        positions = _read_positions(f, n_atoms, n_frames)
        rollup = _read_legacy_rollup(f)

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

        # Per-frame time-series groups (each Optional)
        energy = EnergyAccess(
            gromacs=_load_gromacs_energy_time_series(f),
            bonded=_load_bonded_energy_time_series(f),
        )

        # Water field TR pair (TimeSeries + Welford)
        water_field = WaterFieldAccess(
            time_series=_load_water_field_time_series(f),
            welford=_load_water_field_welford(f),
        )

        # Hydration geometry TR pair
        hydration_geometry = HydrationGeometryAccess(
            time_series=_load_hydration_geometry_time_series(f),
            welford=_load_hydration_geometry_welford(f),
        )

        # Hydration shell TR pair (COM-based)
        hydration_shell = HydrationShellAccess(
            time_series=_load_hydration_shell_time_series(f),
            welford=_load_hydration_shell_welford(f),
        )

        # Per-residue dihedral timeline (first (R, T) TR — 2026-05-19)
        dihedrals = _load_dihedral_time_series(f)
        dihedral_bin_transitions = _load_dihedral_bin_transition(f)
        dssp8 = Dssp8Access(
            time_series=_load_dssp8_time_series(f),
            transitions=_load_dssp8_transition(f),
        )
        ring_pucker = _load_ring_pucker_time_series(f)
        j_coupling = _load_j_coupling_time_series(f)

    return TrajectoryData(
        protein_id=protein_id,
        n_atoms=n_atoms,
        n_frames=n_frames,
        positions=positions,
        frame_times=frame_times,
        rollup=rollup,
        bonds=bonds,
        welford=welford,
        energy=energy,
        water_field=water_field,
        hydration_geometry=hydration_geometry,
        hydration_shell=hydration_shell,
        dihedrals=dihedrals,
        dihedral_bin_transitions=dihedral_bin_transitions,
        dssp8=dssp8,
        ring_pucker=ring_pucker,
        j_coupling=j_coupling,
    )
