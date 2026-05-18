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
    of bottom-N% lowest-energy frames is a numpy one-liner downstream:

        e = traj.energy.gromacs.total_energy
        low_idx = np.argsort(e)[:int(0.1 * len(e))]
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
    cmap: np.ndarray
    total: np.ndarray
    frame_indices: np.ndarray    # (T,) uint64
    frame_times: np.ndarray      # (T,) ps
    units: str
    split_convention: str


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
        cmap=g["cmap"][:],
        total=g["total"][:],
        frame_indices=g["frame_indices"][:],
        frame_times=g["frame_times"][:],
        units=_group_units(g),
        split_convention=str(_decode_attr(
            g.attrs.get("split_convention", ""))),
    )


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
    )
