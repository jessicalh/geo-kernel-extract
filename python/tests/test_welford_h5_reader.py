"""SDK round-trip tests for the six Welford H5 groups.

Verifies that:
  - load_trajectory() exposes /trajectory/<kind>_welford/ as
    WelfordAccess.{bs, hm, mc, eeq, sasa, hbond_count},
  - each per-channel block surfaces as a WelfordMoments with the
    seven (mean, m2, std, min, max, min_frame, max_frame) arrays,
  - per-component channels (T1[3], T2[5]) carry the right inner axis,
  - per-dataset `units` H5 attributes propagate to WelfordMoments.units,
  - group-level attributes (mean_dt_ps, frame_index_range, irrep_layout_*)
    survive the trip,
  - missing groups return None (no exception),
  - legacy aliases (t2mag_*, eeq's bare 'mean'/'std', etc.) are NOT
    exposed on the typed wrapper — canonical names only.

Synthesises trajectory.h5 from scratch using h5py — does not depend on
a real C++ extraction. The schema mirrors *WelfordTrajectoryResult::
WriteH5Group exactly. If the C++ schema drifts, these tests will catch it
when run against a real H5.
"""

import h5py
import numpy as np
import pytest

from nmr_extract import (
    load_trajectory,
    BsWelfordGroup,
    HmWelfordGroup,
    McConnellWelfordGroup,
    EeqWelfordGroup,
    SasaWelfordGroup,
    HBondCountWelfordGroup,
    WelfordAccess,
    WelfordMoments,
)


# Tiny atom count — small enough to hand-verify on test failure.
N_ATOMS = 10


# ─── Synthesis helpers ──────────────────────────────────────────────


def _make_moments_1d(n_atoms: int, seed: int) -> dict:
    """Build the seven 1D arrays of a Welford 7-stat block."""
    rng = np.random.default_rng(seed)
    return {
        "mean":      rng.normal(size=n_atoms).astype(np.float64),
        "m2":        rng.uniform(0.0, 1.0, size=n_atoms).astype(np.float64),
        "std":       rng.uniform(0.0, 0.1, size=n_atoms).astype(np.float64),
        "min":       rng.normal(-2.0, 0.1, size=n_atoms).astype(np.float64),
        "max":       rng.normal( 2.0, 0.1, size=n_atoms).astype(np.float64),
        "min_frame": rng.integers(0, 100, size=n_atoms).astype(np.uint64),
        "max_frame": rng.integers(0, 100, size=n_atoms).astype(np.uint64),
    }


def _make_moments_2d(n_atoms: int, k: int, seed: int) -> dict:
    """Build a per-component (N, K) Welford 7-stat block."""
    rng = np.random.default_rng(seed)
    return {
        "mean":      rng.normal(size=(n_atoms, k)).astype(np.float64),
        "m2":        rng.uniform(0.0, 1.0, size=(n_atoms, k)).astype(np.float64),
        "std":       rng.uniform(0.0, 0.1, size=(n_atoms, k)).astype(np.float64),
        "min":       rng.normal(-2.0, 0.1, size=(n_atoms, k)).astype(np.float64),
        "max":       rng.normal( 2.0, 0.1, size=(n_atoms, k)).astype(np.float64),
        "min_frame": rng.integers(0, 100, size=(n_atoms, k)).astype(np.uint64),
        "max_frame": rng.integers(0, 100, size=(n_atoms, k)).astype(np.uint64),
    }


def _write_moments_block(grp, prefix: str, arrays: dict,
                          base_units: str, m2_units: str) -> None:
    """Emit the seven `<prefix>_<stat>` datasets with units attrs."""
    for stat in ("mean", "m2", "std", "min", "max", "min_frame", "max_frame"):
        ds = grp.create_dataset(f"{prefix}_{stat}", data=arrays[stat])
        if stat == "m2":
            ds.attrs["units"] = m2_units
        elif stat in ("min_frame", "max_frame"):
            ds.attrs["units"] = "frame_index"
        else:
            ds.attrs["units"] = base_units


def _write_required_traj_root(f: h5py.File, n_atoms: int, n_frames: int) -> None:
    """Emit the analysis-schema frame metadata so load_trajectory has
    coherent (T, N, 3) positions + (T,) frame_times before reaching the
    welford groups.

    Mirrors what the current C++ writer produces:
      - Root attrs: protein_id, n_atoms, finalized only
      - /trajectory/frames/{time_ps, original_index} + n_frames group attr
      - /trajectory/positions/xyz atom-major (N*T, 3) +
        n_atoms / n_frames / result_name / finalized group attrs

    No `/rollup` group — that's a legacy GromacsProtein artifact, absent
    in the analysis-mode H5 the SDK is the consumer for. Tests that
    specifically need to exercise the legacy-rollup load path can use
    `_write_legacy_rollup_root` instead.
    """
    f.attrs["protein_id"] = "synthetic_test"
    f.attrs["n_atoms"]    = n_atoms
    f.attrs["finalized"]  = True

    frames_grp = f.create_group("/trajectory/frames")
    frames_grp.attrs["n_frames"] = np.uint64(n_frames)
    frames_grp.create_dataset(
        "time_ps", data=np.linspace(0.0, 1.0, n_frames, dtype=np.float64))
    frames_grp.create_dataset(
        "original_index", data=np.arange(n_frames, dtype=np.uint64))

    pos_grp = f.create_group("/trajectory/positions")
    pos_grp.attrs["result_name"] = "PositionsTimeSeriesTrajectoryResult"
    pos_grp.attrs["n_atoms"]     = np.uint64(n_atoms)
    pos_grp.attrs["n_frames"]    = np.uint64(n_frames)
    pos_grp.attrs["finalized"]   = True
    # Atom-major (N*T, 3) — matches PositionsTimeSeriesTrajectoryResult.cpp
    # layout: for atom i, the T frames are contiguous, then the next atom.
    pos_grp.create_dataset(
        "xyz", data=np.zeros((n_atoms * n_frames, 3), dtype=np.float64))


def _write_legacy_rollup_root(f: h5py.File, n_atoms: int, n_frames: int) -> None:
    """Emit the legacy GromacsProtein ensemble schema (pre-2026-05-17
    trajectory pipeline). Kept so the SDK's backward-compat path can be
    exercised end-to-end — current production output uses the
    `_write_required_traj_root` schema above."""
    f.attrs["n_atoms"]               = n_atoms
    f.attrs["n_frames"]              = n_frames
    f.attrs["positions_shape_T"]     = n_frames
    f.attrs["positions_shape_N"]     = n_atoms
    f.create_dataset("positions",
                     data=np.zeros((n_frames * n_atoms, 3), dtype=np.float64))
    f.create_dataset("frame_times",
                     data=np.linspace(0.0, 1.0, n_frames, dtype=np.float64))
    f.create_dataset("rollup/mean", data=np.zeros((n_atoms, 1), dtype=np.float64))
    f.create_dataset("rollup/std",  data=np.zeros((n_atoms, 1), dtype=np.float64))
    f.create_dataset("rollup/names", data=np.array([b"dummy"], dtype="S16"))


def _write_bs_welford(f: h5py.File, n_atoms: int) -> dict:
    """Synthesise /trajectory/bs_welford/ and return the source dict
    for cross-check on read-back."""
    grp = f.create_group("/trajectory/bs_welford")
    grp.attrs["result_name"]     = "BsWelfordTrajectoryResult"
    grp.attrs["n_frames"]        = np.uint64(750)
    grp.attrs["finalized"]       = True
    grp.attrs["ddof"]            = np.int32(1)
    grp.attrs["mean_dt_ps"]      = 20.0
    grp.attrs["frame_index_range"] = np.array([0, 749], dtype=np.uint64)
    grp.attrs["irrep_layout_t1"] = "v_x,v_y,v_z"
    grp.attrs["irrep_layout_t2"] = "m-2,m-1,m0,m+1,m+2"
    grp.attrs["units"]           = "ppm_T_per_nA"

    base = "ppm_T_per_nA"
    sq   = "ppm_T^2_per_nA^2"
    rate = "ppm_T_per_nA_per_ps"
    rate_sq = "ppm_T^2_per_nA^2_per_ps^2"

    # Scalar channels
    t0          = _make_moments_1d(n_atoms, seed=1)
    t2magnitude = _make_moments_1d(n_atoms, seed=2)
    _write_moments_block(grp, "t0",          t0,          base, sq)
    _write_moments_block(grp, "t2magnitude", t2magnitude, base, sq)

    # Per-component channels
    t1 = _make_moments_2d(n_atoms, 3, seed=3)
    t2 = _make_moments_2d(n_atoms, 5, seed=4)
    _write_moments_block(grp, "t1", t1, base, sq)
    _write_moments_block(grp, "t2", t2, base, sq)

    # Delta channels (T0)
    t0_delta         = _make_moments_1d(n_atoms, seed=5)
    t0_abs_delta     = _make_moments_1d(n_atoms, seed=6)
    t0_delta_squared = _make_moments_1d(n_atoms, seed=7)
    t0_dxdt          = _make_moments_1d(n_atoms, seed=8)
    _write_moments_block(grp, "t0_delta",         t0_delta,         base, sq)
    _write_moments_block(grp, "t0_abs_delta",     t0_abs_delta,     base, sq)
    _write_moments_block(grp, "t0_delta_squared", t0_delta_squared, sq,   "ppm_T^4_per_nA^4")
    _write_moments_block(grp, "t0_dxdt",          t0_dxdt,          rate, rate_sq)

    # Provenance
    # sqrt of a synthetic-Welford mean is allowed to be NaN; the
    # readback path only cares about shape + dtype + units propagation.
    with np.errstate(invalid="ignore"):
        t0_rms_delta = np.sqrt(t0_delta_squared["mean"]).astype(np.float64)
    nframes = np.full(n_atoms, 750, dtype=np.uint64)
    deltan  = np.full(n_atoms, 749, dtype=np.uint64)
    # dxdt_n equals delta_n on a well-formed (no zero-dt) trajectory.
    dxdtn  = np.full(n_atoms, 749, dtype=np.uint64)
    rms_ds = grp.create_dataset("t0_rms_delta",      data=t0_rms_delta)
    rms_ds.attrs["units"] = base
    nf_ds = grp.create_dataset("n_frames_per_atom", data=nframes)
    nf_ds.attrs["units"] = "frame_count"
    dxn_ds = grp.create_dataset("dxdt_n_per_atom", data=dxdtn)
    dxn_ds.attrs["units"] = "frame_count"
    dn_ds = grp.create_dataset("delta_n_per_atom",  data=deltan)
    dn_ds.attrs["units"] = "frame_count"

    # Legacy aliases (deprecated; reader should ignore)
    for name, src in (
        ("t2mag_mean", t2magnitude["mean"]),
        ("t2mag_std",  t2magnitude["std"]),
        ("t2mag_min",  t2magnitude["min"]),
        ("t2mag_max",  t2magnitude["max"]),
    ):
        ds = grp.create_dataset(name, data=src)
        ds.attrs["units"]          = base
        ds.attrs["deprecated_use"] = f"canonical name is t2magnitude_{name.split('_', 1)[1]}"

    return {
        "t0": t0, "t2magnitude": t2magnitude, "t1": t1, "t2": t2,
        "t0_delta": t0_delta, "t0_abs_delta": t0_abs_delta,
        "t0_delta_squared": t0_delta_squared, "t0_dxdt": t0_dxdt,
        "t0_rms_delta": t0_rms_delta,
        "n_frames_per_atom": nframes,
        "delta_n_per_atom": deltan,
    }


def _write_eeq_welford(f: h5py.File, n_atoms: int) -> dict:
    grp = f.create_group("/trajectory/eeq_welford")
    grp.attrs["result_name"]      = "EeqWelfordTrajectoryResult"
    grp.attrs["n_frames"]         = np.uint64(750)
    grp.attrs["finalized"]        = True
    grp.attrs["ddof"]             = np.int32(1)
    grp.attrs["mean_dt_ps"]       = 20.0
    grp.attrs["frame_index_range"] = np.array([0, 749], dtype=np.uint64)
    grp.attrs["units"]            = "elementary_charge"

    base = "elementary_charge"
    sq   = "elementary_charge^2"
    rate = "elementary_charge_per_ps"
    rate_sq = "elementary_charge^2_per_ps^2"

    src = {}
    for i, name in enumerate(("charge", "charge_delta", "charge_abs_delta",
                              "charge_delta_squared", "charge_dxdt")):
        src[name] = _make_moments_1d(n_atoms, seed=100 + i)
    _write_moments_block(grp, "charge",               src["charge"],               base, sq)
    _write_moments_block(grp, "charge_delta",         src["charge_delta"],         base, sq)
    _write_moments_block(grp, "charge_abs_delta",     src["charge_abs_delta"],     base, sq)
    _write_moments_block(grp, "charge_delta_squared", src["charge_delta_squared"], sq,   "elementary_charge^4")
    _write_moments_block(grp, "charge_dxdt",          src["charge_dxdt"],          rate, rate_sq)

    with np.errstate(invalid="ignore"):
        rms = np.sqrt(src["charge_delta_squared"]["mean"]).astype(np.float64)
    grp.create_dataset("rms_delta", data=rms).attrs["units"] = base
    grp.create_dataset("n_frames_per_atom",
                       data=np.full(n_atoms, 750, dtype=np.uint64)).attrs["units"] = "frame_count"
    grp.create_dataset("delta_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    # dxdt_n equals delta_n on a well-formed (no zero-dt) trajectory.
    grp.create_dataset("dxdt_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    return src


def _write_tensor_welford(f: h5py.File, group_path: str,
                          base: str, sq: str, rate: str, rate_sq: str,
                          delta_sq_unit: str, delta_sq_m2_unit: str,
                          n_atoms: int) -> dict:
    """Shared writer for HM and Mc Welford groups (same shape as BS)."""
    grp = f.create_group(group_path)
    grp.attrs["result_name"]     = group_path.rsplit("/", 1)[-1]
    grp.attrs["n_frames"]        = np.uint64(750)
    grp.attrs["finalized"]       = True
    grp.attrs["ddof"]            = np.int32(1)
    grp.attrs["mean_dt_ps"]      = 20.0
    grp.attrs["frame_index_range"] = np.array([0, 749], dtype=np.uint64)
    grp.attrs["irrep_layout_t1"] = "v_x,v_y,v_z"
    grp.attrs["irrep_layout_t2"] = "m-2,m-1,m0,m+1,m+2"
    grp.attrs["units"]           = base

    seed_base = abs(hash(group_path)) % 10000
    t0          = _make_moments_1d(n_atoms, seed=seed_base + 1)
    t2magnitude = _make_moments_1d(n_atoms, seed=seed_base + 2)
    t1 = _make_moments_2d(n_atoms, 3, seed=seed_base + 3)
    t2 = _make_moments_2d(n_atoms, 5, seed=seed_base + 4)
    _write_moments_block(grp, "t0",          t0,          base, sq)
    _write_moments_block(grp, "t2magnitude", t2magnitude, base, sq)
    _write_moments_block(grp, "t1", t1, base, sq)
    _write_moments_block(grp, "t2", t2, base, sq)

    t0_delta         = _make_moments_1d(n_atoms, seed=seed_base + 5)
    t0_abs_delta     = _make_moments_1d(n_atoms, seed=seed_base + 6)
    t0_delta_squared = _make_moments_1d(n_atoms, seed=seed_base + 7)
    t0_dxdt          = _make_moments_1d(n_atoms, seed=seed_base + 8)
    _write_moments_block(grp, "t0_delta",         t0_delta,         base, sq)
    _write_moments_block(grp, "t0_abs_delta",     t0_abs_delta,     base, sq)
    _write_moments_block(grp, "t0_delta_squared", t0_delta_squared,
                          delta_sq_unit, delta_sq_m2_unit)
    _write_moments_block(grp, "t0_dxdt",          t0_dxdt,          rate, rate_sq)

    with np.errstate(invalid="ignore"):
        rms = np.sqrt(t0_delta_squared["mean"]).astype(np.float64)
    grp.create_dataset("t0_rms_delta", data=rms).attrs["units"] = base
    grp.create_dataset("n_frames_per_atom",
                       data=np.full(n_atoms, 750, dtype=np.uint64)).attrs["units"] = "frame_count"
    grp.create_dataset("delta_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    # dxdt_n equals delta_n on a well-formed (no zero-dt) trajectory.
    grp.create_dataset("dxdt_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    return {
        "t0": t0, "t2magnitude": t2magnitude, "t1": t1, "t2": t2,
        "t0_delta": t0_delta, "t0_abs_delta": t0_abs_delta,
        "t0_delta_squared": t0_delta_squared, "t0_dxdt": t0_dxdt,
    }


def _write_sasa_welford(f: h5py.File, n_atoms: int) -> dict:
    grp = f.create_group("/trajectory/sasa_welford")
    grp.attrs["result_name"]      = "SasaWelfordTrajectoryResult"
    grp.attrs["n_frames"]         = np.uint64(750)
    grp.attrs["finalized"]        = True
    grp.attrs["ddof"]             = np.int32(1)
    grp.attrs["mean_dt_ps"]       = 20.0
    grp.attrs["frame_index_range"] = np.array([0, 749], dtype=np.uint64)
    grp.attrs["units"]            = "Angstrom^2"

    src = {}
    for i, name in enumerate(("sasa", "sasa_delta", "sasa_abs_delta",
                              "sasa_delta_squared", "sasa_dxdt")):
        src[name] = _make_moments_1d(n_atoms, seed=300 + i)
    _write_moments_block(grp, "sasa",               src["sasa"],               "Angstrom^2",        "Angstrom^4")
    _write_moments_block(grp, "sasa_delta",         src["sasa_delta"],         "Angstrom^2",        "Angstrom^4")
    _write_moments_block(grp, "sasa_abs_delta",     src["sasa_abs_delta"],     "Angstrom^2",        "Angstrom^4")
    _write_moments_block(grp, "sasa_delta_squared", src["sasa_delta_squared"], "Angstrom^4",        "Angstrom^8")
    _write_moments_block(grp, "sasa_dxdt",          src["sasa_dxdt"],          "Angstrom^2_per_ps", "Angstrom^4_per_ps^2")

    with np.errstate(invalid="ignore"):
        rms = np.sqrt(src["sasa_delta_squared"]["mean"]).astype(np.float64)
    grp.create_dataset("rms_delta", data=rms).attrs["units"] = "Angstrom^2"
    grp.create_dataset("n_frames_per_atom",
                       data=np.full(n_atoms, 750, dtype=np.uint64)).attrs["units"] = "frame_count"
    grp.create_dataset("delta_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    # dxdt_n equals delta_n on a well-formed (no zero-dt) trajectory.
    grp.create_dataset("dxdt_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    return src


def _write_hbond_count_welford(f: h5py.File, n_atoms: int) -> dict:
    grp = f.create_group("/trajectory/hbond_count_welford")
    grp.attrs["result_name"]      = "HBondCountWelfordTrajectoryResult"
    grp.attrs["n_frames"]         = np.uint64(750)
    grp.attrs["finalized"]        = True
    grp.attrs["ddof"]             = np.int32(1)
    grp.attrs["mean_dt_ps"]       = 20.0
    grp.attrs["frame_index_range"] = np.array([0, 749], dtype=np.uint64)
    grp.attrs["units"]            = "pairs"
    grp.attrs["source_radius_A"]  = 3.5

    src = {}
    for i, name in enumerate(("count", "count_delta", "count_abs_delta",
                              "count_delta_squared", "count_dxdt")):
        src[name] = _make_moments_1d(n_atoms, seed=200 + i)
    _write_moments_block(grp, "count",               src["count"],               "pairs",        "pairs^2")
    _write_moments_block(grp, "count_delta",         src["count_delta"],         "pairs",        "pairs^2")
    _write_moments_block(grp, "count_abs_delta",     src["count_abs_delta"],     "pairs",        "pairs^2")
    _write_moments_block(grp, "count_delta_squared", src["count_delta_squared"], "pairs^2",      "pairs^4")
    _write_moments_block(grp, "count_dxdt",          src["count_dxdt"],          "pairs_per_ps", "pairs^2_per_ps^2")

    src["occupancy_fraction"] = _make_moments_1d(n_atoms, seed=250)
    # Clamp the synthetic data to [0, 1] to mirror the C++ contract.
    for stat in ("mean", "min", "max"):
        src["occupancy_fraction"][stat] = np.clip(src["occupancy_fraction"][stat], 0.0, 1.0)
    _write_moments_block(grp, "occupancy_fraction", src["occupancy_fraction"],
                         "dimensionless", "dimensionless")

    with np.errstate(invalid="ignore"):
        rms = np.sqrt(src["count_delta_squared"]["mean"]).astype(np.float64)
    grp.create_dataset("rms_delta", data=rms).attrs["units"] = "pairs"
    grp.create_dataset("n_frames_per_atom",
                       data=np.full(n_atoms, 750, dtype=np.uint64)).attrs["units"] = "frame_count"
    grp.create_dataset("delta_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    # dxdt_n equals delta_n on a well-formed (no zero-dt) trajectory.
    grp.create_dataset("dxdt_n_per_atom",
                       data=np.full(n_atoms, 749, dtype=np.uint64)).attrs["units"] = "frame_count"
    return src


# ─── Fixtures ───────────────────────────────────────────────────────


@pytest.fixture
def h5_with_bs_welford(tmp_path):
    """Trajectory H5 containing only the BS Welford group."""
    path = tmp_path / "trajectory.h5"
    with h5py.File(path, "w") as f:
        _write_required_traj_root(f, N_ATOMS, n_frames=4)
        src = _write_bs_welford(f, N_ATOMS)
    return path, src


@pytest.fixture
def h5_with_bs_eeq_hbond(tmp_path):
    """BS + EEq + HBondCount; HM/Mc/SASA intentionally omitted."""
    path = tmp_path / "trajectory.h5"
    with h5py.File(path, "w") as f:
        _write_required_traj_root(f, N_ATOMS, n_frames=4)
        bs_src    = _write_bs_welford(f, N_ATOMS)
        eeq_src   = _write_eeq_welford(f, N_ATOMS)
        hbond_src = _write_hbond_count_welford(f, N_ATOMS)
    return path, {"bs": bs_src, "eeq": eeq_src, "hbond_count": hbond_src}


@pytest.fixture
def h5_with_all_six(tmp_path):
    """All six Welford groups present — exercises every loader."""
    path = tmp_path / "trajectory.h5"
    with h5py.File(path, "w") as f:
        _write_required_traj_root(f, N_ATOMS, n_frames=4)
        _write_bs_welford(f, N_ATOMS)
        _write_tensor_welford(
            f, "/trajectory/hm_welford",
            base="Angstrom^-1", sq="Angstrom^-2",
            rate="Angstrom^-1_per_ps", rate_sq="Angstrom^-2_per_ps^2",
            delta_sq_unit="Angstrom^-2", delta_sq_m2_unit="Angstrom^-4",
            n_atoms=N_ATOMS)
        _write_tensor_welford(
            f, "/trajectory/mc_welford",
            base="Angstrom^-3", sq="Angstrom^-6",
            rate="Angstrom^-3_per_ps", rate_sq="Angstrom^-6_per_ps^2",
            delta_sq_unit="Angstrom^-6", delta_sq_m2_unit="Angstrom^-12",
            n_atoms=N_ATOMS)
        _write_eeq_welford(f, N_ATOMS)
        _write_sasa_welford(f, N_ATOMS)
        _write_hbond_count_welford(f, N_ATOMS)
    return path


# ─── Test cases ─────────────────────────────────────────────────────


class TestBsWelfordReadback:
    """(a) Synthetic fixture round-trip for BS Welford group."""

    def test_welford_access_attached(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        assert isinstance(traj.welford, WelfordAccess)
        assert isinstance(traj.welford.bs, BsWelfordGroup)

    def test_scalar_channel_shape(self, h5_with_bs_welford):
        path, src = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.welford.bs.t0.mean.shape == (N_ATOMS,)
        assert traj.welford.bs.t2magnitude.mean.shape == (N_ATOMS,)
        np.testing.assert_array_equal(traj.welford.bs.t0.mean, src["t0"]["mean"])
        np.testing.assert_array_equal(
            traj.welford.bs.t2magnitude.std, src["t2magnitude"]["std"])

    def test_per_component_t1_shape(self, h5_with_bs_welford):
        path, src = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.welford.bs.t1.mean.shape == (N_ATOMS, 3)
        np.testing.assert_array_equal(traj.welford.bs.t1.mean, src["t1"]["mean"])

    def test_per_component_t2_shape(self, h5_with_bs_welford):
        path, src = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.welford.bs.t2.mean.shape == (N_ATOMS, 5)
        np.testing.assert_array_equal(traj.welford.bs.t2.std, src["t2"]["std"])

    def test_delta_channels_present(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        for chan_name in ("t0_delta", "t0_abs_delta",
                          "t0_delta_squared", "t0_dxdt"):
            chan = getattr(bs, chan_name)
            assert isinstance(chan, WelfordMoments)
            assert chan.mean.shape == (N_ATOMS,)

    def test_units_attribute_per_channel(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.t0.units == "ppm_T_per_nA"
        assert bs.t1.units == "ppm_T_per_nA"
        assert bs.t2.units == "ppm_T_per_nA"
        # dxdt is rate-shaped — different unit string
        assert bs.t0_dxdt.units == "ppm_T_per_nA_per_ps"
        # delta_squared values are in squared units
        assert bs.t0_delta_squared.units == "ppm_T^2_per_nA^2"

    def test_group_level_attributes(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.mean_dt_ps == 20.0
        assert bs.frame_index_range == (0, 749)
        assert bs.irrep_layout_t1 == "v_x,v_y,v_z"
        assert bs.irrep_layout_t2 == "m-2,m-1,m0,m+1,m+2"
        assert bs.units == "ppm_T_per_nA"

    def test_provenance_arrays(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.t0_rms_delta.shape == (N_ATOMS,)
        assert bs.n_frames_per_atom.shape == (N_ATOMS,)
        assert bs.delta_n_per_atom.shape == (N_ATOMS,)
        np.testing.assert_array_equal(
            bs.n_frames_per_atom, np.full(N_ATOMS, 750, dtype=np.uint64))


class TestMissingGroupResilience:
    """(b) Groups that aren't in the H5 surface as None, no error."""

    def test_bs_present_but_hm_missing(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.welford.bs is not None
        assert traj.welford.hm is None
        assert traj.welford.mc is None
        assert traj.welford.eeq is None
        assert traj.welford.sasa is None
        assert traj.welford.hbond_count is None

    def test_partial_group_set(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        assert isinstance(traj.welford.bs, BsWelfordGroup)
        assert isinstance(traj.welford.eeq, EeqWelfordGroup)
        assert isinstance(traj.welford.hbond_count, HBondCountWelfordGroup)
        assert traj.welford.hm is None
        assert traj.welford.mc is None
        assert traj.welford.sasa is None


class TestLegacyAliasNotExposed:
    """(c) Legacy t2mag_* and friends are not surfaced on the typed wrapper.

    The C++ writer emits them for backward compatibility with pre-Phase-2b
    H5 consumers. The SDK only exposes the canonical names.
    """

    def test_no_t2mag_attribute_on_bs_group(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        # Canonical channel is t2magnitude
        assert hasattr(bs, "t2magnitude")
        # The dataclass field set is fixed (frozen). Confirm by direct
        # attribute access — should raise.
        for legacy_name in ("t2mag_mean", "t2mag_std", "t2mag_min", "t2mag_max"):
            with pytest.raises(AttributeError):
                getattr(bs, legacy_name)


class TestEeqScalarGroup:
    """(d) EEq is the canonical scalar (non-tensor, non-occupancy) shape."""

    def test_charge_channel_shape(self, h5_with_bs_eeq_hbond):
        path, src = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        eeq = traj.welford.eeq
        assert eeq.charge.mean.shape == (N_ATOMS,)
        np.testing.assert_array_equal(eeq.charge.mean,
                                      src["eeq"]["charge"]["mean"])

    def test_charge_dxdt_units(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        assert traj.welford.eeq.charge_dxdt.units == "elementary_charge_per_ps"

    def test_eeq_lacks_tensor_channels(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        # EEq is dataclass-frozen; should have no t1/t2/t2magnitude
        for tensor_name in ("t1", "t2", "t2magnitude"):
            with pytest.raises(AttributeError):
                getattr(traj.welford.eeq, tensor_name)


class TestHBondCountOccupancyFraction:
    """(e) HBondCount adds occupancy_fraction; units must be per-dataset."""

    def test_occupancy_fraction_is_dimensionless(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        assert traj.welford.hbond_count.occupancy_fraction.units == "dimensionless"

    def test_count_units_separate_from_occupancy(self, h5_with_bs_eeq_hbond):
        """Per-dataset units, not group-level: count_mean's units differ
        from occupancy_fraction_mean's units within the same H5 group."""
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        hc = traj.welford.hbond_count
        assert hc.count.units == "pairs"
        assert hc.occupancy_fraction.units == "dimensionless"
        # Group-level units attribute is whichever the writer set —
        # the typed wrapper preserves it on the group, but each channel
        # carries its own.
        assert hc.units == "pairs"

    def test_occupancy_fraction_shape(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        of = traj.welford.hbond_count.occupancy_fraction
        assert of.mean.shape == (N_ATOMS,)
        assert of.min_frame.shape == (N_ATOMS,)

    def test_source_radius_A_attr(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        assert traj.welford.hbond_count.source_radius_A == pytest.approx(3.5)


class TestAllSixGroups:
    """Smoke coverage for HM / Mc / SASA loaders (BS / EEq / HBondCount
    covered above)."""

    def test_all_six_present(self, h5_with_all_six):
        traj = load_trajectory(h5_with_all_six)
        assert isinstance(traj.welford.bs, BsWelfordGroup)
        assert isinstance(traj.welford.hm, HmWelfordGroup)
        assert isinstance(traj.welford.mc, McConnellWelfordGroup)
        assert isinstance(traj.welford.eeq, EeqWelfordGroup)
        assert isinstance(traj.welford.sasa, SasaWelfordGroup)
        assert isinstance(traj.welford.hbond_count, HBondCountWelfordGroup)

    def test_hm_tensor_shape_and_units(self, h5_with_all_six):
        traj = load_trajectory(h5_with_all_six)
        hm = traj.welford.hm
        assert hm.t1.mean.shape == (N_ATOMS, 3)
        assert hm.t2.mean.shape == (N_ATOMS, 5)
        assert hm.t0.units == "Angstrom^-1"
        assert hm.t0_dxdt.units == "Angstrom^-1_per_ps"

    def test_mc_tensor_shape_and_units(self, h5_with_all_six):
        traj = load_trajectory(h5_with_all_six)
        mc = traj.welford.mc
        assert mc.t2.mean.shape == (N_ATOMS, 5)
        assert mc.t0.units == "Angstrom^-3"

    def test_sasa_scalar_shape_and_units(self, h5_with_all_six):
        traj = load_trajectory(h5_with_all_six)
        sasa = traj.welford.sasa
        assert sasa.sasa.mean.shape == (N_ATOMS,)
        assert sasa.sasa.units == "Angstrom^2"
        assert sasa.sasa_dxdt.units == "Angstrom^2_per_ps"


class TestTrajectoryWithoutWelfordGroups:
    """Pre-Phase-2b trajectory.h5 (no /trajectory/<kind>_welford/ groups)
    still loads — WelfordAccess is empty but present."""

    def test_no_welford_groups(self, tmp_path):
        path = tmp_path / "trajectory.h5"
        with h5py.File(path, "w") as f:
            _write_required_traj_root(f, N_ATOMS, n_frames=4)
        traj = load_trajectory(path)
        assert isinstance(traj.welford, WelfordAccess)
        assert traj.welford.bs is None
        assert traj.welford.hm is None
        assert traj.welford.mc is None
        assert traj.welford.eeq is None
        assert traj.welford.sasa is None
        assert traj.welford.hbond_count is None


# ─── M2 units honesty (codex MEDIUM finding 2026-05-18) ─────────────


class TestPerDatasetM2Units:
    """C++ writes distinct unit strings for `<prefix>_m2` (squared
    accumulator) vs `<prefix>_mean` (sample). The SDK preserves both
    via `WelfordMoments.units` + `WelfordMoments.m2_units`."""

    def test_bs_t0_m2_units_squared(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.t0.units == "ppm_T_per_nA"
        assert bs.t0.m2_units == "ppm_T^2_per_nA^2"

    def test_bs_t1_m2_units_squared(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.t1.units == "ppm_T_per_nA"
        assert bs.t1.m2_units == "ppm_T^2_per_nA^2"

    def test_bs_t0_delta_squared_units_pre_squared(self, h5_with_bs_welford):
        """The *_delta_squared channel stores Δ² samples in base² units
        (e.g. ppm_T^2_per_nA^2). Its m2 accumulator is THEN squared
        again — base⁴ — because Welford M2 = sum (sample-mean)²."""
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.t0_delta_squared.units == "ppm_T^2_per_nA^2"
        assert bs.t0_delta_squared.m2_units == "ppm_T^4_per_nA^4"

    def test_bs_t0_dxdt_rate_units(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.t0_dxdt.units == "ppm_T_per_nA_per_ps"
        assert bs.t0_dxdt.m2_units == "ppm_T^2_per_nA^2_per_ps^2"

    def test_hbond_occupancy_fraction_dimensionless(self, h5_with_bs_eeq_hbond):
        """occupancy_fraction's mean and m2 are BOTH dimensionless —
        the indicator is unitless so its square is too."""
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        hb = traj.welford.hbond_count
        assert hb.occupancy_fraction.units == "dimensionless"
        assert hb.occupancy_fraction.m2_units == "dimensionless"


# ─── Production-schema (analysis H5) round-trip ─────────────────────


class TestProductionSchema:
    """Verify the new schema written by current C++ TrajectoryProtein
    loads: /trajectory/frames/time_ps + /trajectory/positions/xyz +
    root attrs of just protein_id/n_atoms/finalized, with NO /rollup."""

    def test_protein_id_and_n_atoms_from_root_attrs(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.protein_id == "synthetic_test"
        assert traj.n_atoms == N_ATOMS

    def test_n_frames_from_trajectory_frames_group(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.n_frames == 4  # set by _write_required_traj_root
        assert traj.frame_times.shape == (4,)
        # Monotone non-decreasing — sanity, not physics
        assert np.all(np.diff(traj.frame_times) >= 0.0)

    def test_positions_shape_T_N_3(self, h5_with_bs_welford):
        """Production schema: atom-major (N*T, 3) on disk → (T, N, 3)
        in-memory after the reader's transpose."""
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.positions.shape == (4, N_ATOMS, 3)

    def test_no_rollup_in_production_schema(self, h5_with_bs_welford):
        """Production trajectory.h5 has no /rollup — analysis mode
        replaces it with the Welford H5 groups."""
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        assert traj.rollup is None


# ─── Legacy GromacsProtein schema backward-compat ──────────────────


class TestDxdtNPerAtom:
    """`dxdt_n_per_atom` is a new H5 dataset (codex 2026-05-18 fix). It
    tracks the count of VALID-dt samples that contributed to the dxdt
    accumulator — distinct from `delta_n_per_atom` because zero-dt
    frames are skipped rather than zero-filled.

    On a well-formed fixture trajectory, dxdt_n == delta_n at every
    atom (no duplicated timestamps). For older H5 files that predate
    the fix, the SDK falls back to `delta_n_per_atom` so reads don't
    crash."""

    def test_bs_dxdt_n_present_and_equals_delta_n(self, h5_with_bs_welford):
        path, _ = h5_with_bs_welford
        traj = load_trajectory(path)
        bs = traj.welford.bs
        assert bs.dxdt_n_per_atom.shape == (N_ATOMS,)
        np.testing.assert_array_equal(bs.dxdt_n_per_atom, bs.delta_n_per_atom)

    def test_eeq_dxdt_n_present(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        eeq = traj.welford.eeq
        assert eeq.dxdt_n_per_atom.shape == (N_ATOMS,)
        np.testing.assert_array_equal(eeq.dxdt_n_per_atom, eeq.delta_n_per_atom)

    def test_hbond_dxdt_n_present(self, h5_with_bs_eeq_hbond):
        path, _ = h5_with_bs_eeq_hbond
        traj = load_trajectory(path)
        hb = traj.welford.hbond_count
        assert hb.dxdt_n_per_atom.shape == (N_ATOMS,)

    def test_legacy_h5_without_dxdt_n_falls_back_to_delta_n(self, tmp_path):
        """Older H5 files (pre-codex-fix) don't have dxdt_n_per_atom.
        Reader falls back to delta_n_per_atom and the dataclass still
        constructs."""
        path = tmp_path / "trajectory.h5"
        with h5py.File(path, "w") as f:
            _write_required_traj_root(f, N_ATOMS, n_frames=4)
            _write_bs_welford(f, N_ATOMS)
            # Drop the new dataset to simulate a pre-fix H5 file.
            del f["/trajectory/bs_welford/dxdt_n_per_atom"]
        traj = load_trajectory(path)
        bs = traj.welford.bs
        # Fallback: dxdt_n_per_atom is None of delta_n_per_atom (same array).
        np.testing.assert_array_equal(bs.dxdt_n_per_atom, bs.delta_n_per_atom)


class TestLegacyRollupSchema:
    """Verify the legacy GromacsProtein ensemble schema still loads.
    `/positions` + `/frame_times` + `/rollup/{mean,std,names}` at root,
    root attrs include n_frames + positions_shape_T/N."""

    def test_legacy_rollup_populated(self, tmp_path):
        path = tmp_path / "legacy.h5"
        with h5py.File(path, "w") as f:
            _write_legacy_rollup_root(f, N_ATOMS, n_frames=4)
        traj = load_trajectory(path)
        assert traj.rollup is not None
        assert traj.rollup.n_columns == 1
        assert traj.n_atoms == N_ATOMS
        assert traj.n_frames == 4
        assert traj.frame_times.shape == (4,)
        assert traj.positions.shape == (4, N_ATOMS, 3)

    def test_legacy_schema_no_welford_groups(self, tmp_path):
        path = tmp_path / "legacy.h5"
        with h5py.File(path, "w") as f:
            _write_legacy_rollup_root(f, N_ATOMS, n_frames=4)
        traj = load_trajectory(path)
        assert traj.welford.bs is None
        assert traj.welford.hm is None
        assert traj.welford.mc is None
        assert traj.welford.eeq is None
        assert traj.welford.sasa is None
        assert traj.welford.hbond_count is None
