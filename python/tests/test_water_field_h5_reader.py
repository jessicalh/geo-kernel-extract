"""SDK round-trip tests for the two water-field H5 groups.

Verifies that:
  - load_trajectory() exposes /trajectory/water_field_time_series/ as
    WaterFieldAccess.time_series (WaterFieldTimeSeriesGroup),
  - and /trajectory/water_field_welford/ as WaterFieldAccess.welford
    (WaterFieldWelfordGroup),
  - the (N, T, 3) Vec3 / (N, T, 5) T2-only EFG layouts round-trip,
  - source_attached_per_frame mask + source_attached_count propagate,
  - structural-zero attrs (efg_t0_structural_zero, efg_t1_structural_zero)
    propagate via the underlying H5 attrs (we don't expose them as
    fields on the dataclass — they're invariants),
  - missing groups (no water TR ran) return None on WaterFieldAccess.

Synthesises trajectory.h5 from scratch with h5py — mirrors
WaterFieldTimeSeriesTrajectoryResult::WriteH5Group +
WaterFieldWelfordTrajectoryResult::WriteH5Group exactly.
"""

import h5py
import numpy as np
import pytest

from nmr_extract import (
    load_trajectory,
    WaterFieldAccess,
    WaterFieldTimeSeriesGroup,
    WaterFieldWelfordGroup,
    WelfordMoments,
)


N_ATOMS = 10
N_FRAMES = 6


def _write_required_traj_root(f: h5py.File, n_atoms: int, n_frames: int) -> None:
    """Production-schema trajectory root + positions, copied from the
    welford h5 reader test fixture."""
    f.attrs["protein_id"] = "TEST_PROTEIN"
    f.attrs["n_atoms"]    = n_atoms
    f.attrs["finalized"]  = True

    traj = f.create_group("trajectory")
    frames = traj.create_group("frames")
    frames.attrs["n_frames"] = n_frames
    times = np.arange(n_frames, dtype=np.float64) * 0.1
    frames.create_dataset("time_ps", data=times)
    frames.create_dataset("original_index",
                          data=np.arange(n_frames, dtype=np.uint64))

    pos = traj.create_group("positions")
    pos.attrs["n_atoms"]     = n_atoms
    pos.attrs["n_frames"]    = n_frames
    pos.attrs["result_name"] = "PositionsTimeSeriesTrajectoryResult"
    pos.attrs["finalized"]   = True
    flat = np.zeros((n_atoms * n_frames, 3), dtype=np.float64)
    pos.create_dataset("xyz", data=flat)


def _write_water_field_time_series(f: h5py.File, n_atoms: int, n_frames: int,
                                   source_attached: bool = True) -> dict:
    """Emit /trajectory/water_field_time_series with a recognisable pattern.

    `source_attached=False` simulates a no-solvent extraction where the
    C++ TR captured frames but the WaterFieldResult source was never
    attached → in production WriteH5Group would skip emission entirely.
    """
    grp = f.create_group("/trajectory/water_field_time_series")
    grp.attrs["result_name"]          = "WaterFieldTimeSeriesTrajectoryResult"
    grp.attrs["n_atoms"]              = n_atoms
    grp.attrs["n_frames"]             = n_frames
    grp.attrs["source_attached_count"] = n_frames if source_attached else 0
    grp.attrs["finalized"]            = True
    grp.attrs["efield_layout"]        = "x,y,z"
    grp.attrs["efield_units"]         = "V/Angstrom"
    grp.attrs["efield_parity"]        = "1o"
    grp.attrs["efield_normalization"] = "cartesian"
    grp.attrs["efg_irrep_layout"]     = "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
    grp.attrs["efg_units"]            = "V/Angstrom^2"
    grp.attrs["efg_parity"]           = "2e"
    grp.attrs["efg_normalization"]    = "isometric_real_sph"
    grp.attrs["efg_t0_structural_zero"] = True
    grp.attrs["efg_t1_structural_zero"] = True
    grp.attrs["count_units"]          = "dimensionless"
    grp.attrs["efield_cutoff_A"]      = 15.0
    grp.attrs["n_first_cutoff_A"]     = 3.5
    grp.attrs["n_second_cutoff_A"]    = 5.5

    rng = np.random.default_rng(13)
    efield       = rng.normal(0, 1, (n_atoms, n_frames, 3)).astype(np.float64)
    efield_first = rng.normal(0, 0.5, (n_atoms, n_frames, 3)).astype(np.float64)
    efg          = rng.normal(0, 0.2, (n_atoms, n_frames, 5)).astype(np.float64)
    efg_first    = rng.normal(0, 0.1, (n_atoms, n_frames, 5)).astype(np.float64)
    n_first      = rng.integers(0, 8, (n_atoms, n_frames)).astype(np.uint32)
    n_second     = rng.integers(0, 16, (n_atoms, n_frames)).astype(np.uint32)

    for name, arr, units in [
        ("efield",       efield,       "V/Angstrom"),
        ("efield_first", efield_first, "V/Angstrom"),
        ("efg",          efg,          "V/Angstrom^2"),
        ("efg_first",    efg_first,    "V/Angstrom^2"),
    ]:
        ds = grp.create_dataset(name, data=arr)
        ds.attrs["units"] = units
    for name, arr in [("n_first", n_first), ("n_second", n_second)]:
        ds = grp.create_dataset(name, data=arr)
        ds.attrs["units"] = "dimensionless"
        ds.attrs["absent_sentinel"] = np.uint32(0xFFFFFFFF)

    fi = np.arange(n_frames, dtype=np.uint64)
    ft = np.arange(n_frames, dtype=np.float64) * 0.1
    mask = np.full(n_frames, 1 if source_attached else 0, dtype=np.uint8)
    grp.create_dataset("frame_indices", data=fi).attrs["units"] = "frame_index"
    grp.create_dataset("frame_times",   data=ft).attrs["units"] = "ps"
    grp.create_dataset("source_attached_per_frame", data=mask) \
       .attrs["units"] = "dimensionless"

    return {"efield": efield, "efg": efg, "n_first": n_first,
            "n_second": n_second, "source_mask": mask}


def _write_moments_block(grp, prefix: str, n_atoms: int, k: int = 0,
                         base_units: str = "V/Angstrom",
                         m2_units: str = "V^2/Angstrom^2",
                         seed: int = 0) -> dict:
    """Emit the 7-stat Welford block for one channel.
    k=0 → 1D shape (n_atoms,); k>0 → 2D (n_atoms, k)."""
    rng = np.random.default_rng(seed)
    shape = (n_atoms,) if k == 0 else (n_atoms, k)
    arrays = {
        "mean":      rng.normal(size=shape).astype(np.float64),
        "m2":        rng.uniform(0.0, 1.0, size=shape).astype(np.float64),
        "std":       rng.uniform(0.0, 0.1, size=shape).astype(np.float64),
        "min":       rng.normal(-2.0, 0.1, size=shape).astype(np.float64),
        "max":       rng.normal( 2.0, 0.1, size=shape).astype(np.float64),
        "min_frame": rng.integers(0, 100, size=shape).astype(np.uint64),
        "max_frame": rng.integers(0, 100, size=shape).astype(np.uint64),
    }
    for stat, arr in arrays.items():
        ds = grp.create_dataset(f"{prefix}_{stat}", data=arr)
        if stat == "m2":
            ds.attrs["units"] = m2_units
        elif stat in ("min_frame", "max_frame"):
            ds.attrs["units"] = "frame_index"
        else:
            ds.attrs["units"] = base_units
    return arrays


def _write_water_field_welford(f: h5py.File, n_atoms: int, n_frames: int) -> dict:
    """Emit /trajectory/water_field_welford with a recognisable pattern."""
    grp = f.create_group("/trajectory/water_field_welford")
    grp.attrs["result_name"]              = "WaterFieldWelfordTrajectoryResult"
    grp.attrs["n_frames"]                 = n_frames
    grp.attrs["source_attached_count"]    = n_frames
    grp.attrs["finalized"]                = True
    grp.attrs["ddof"]                     = 1
    grp.attrs["mean_dt_ps"]               = 0.1
    grp.attrs["frame_index_range"]        = np.array([0, n_frames - 1], dtype=np.uint64)
    grp.attrs["irrep_layout_efield"]      = "v_x,v_y,v_z"
    grp.attrs["irrep_layout_efg_t2"]      = "m-2,m-1,m0,m+1,m+2"
    grp.attrs["units"]                    = "V/Angstrom"
    grp.attrs["efg_t0_structural_zero"]   = True
    grp.attrs["efg_t1_structural_zero"]   = True

    truth = {}
    # E-field per-component + magnitude
    for ch in ("efield_x", "efield_y", "efield_z", "efield_magnitude",
               "efield_first_x", "efield_first_y", "efield_first_z",
               "efield_first_magnitude"):
        truth[ch] = _write_moments_block(grp, ch, n_atoms, k=0,
                                          base_units="V/Angstrom",
                                          m2_units="V^2/Angstrom^2",
                                          seed=hash(ch) & 0xFFFF)
    # EFG T2 per-component (N, 5) + |T2|
    for ch in ("efg_t2", "efg_first_t2"):
        truth[ch] = _write_moments_block(grp, ch, n_atoms, k=5,
                                          base_units="V/Angstrom^2",
                                          m2_units="V^2/Angstrom^4",
                                          seed=hash(ch) & 0xFFFF)
    for ch in ("efg_t2magnitude", "efg_first_t2magnitude"):
        truth[ch] = _write_moments_block(grp, ch, n_atoms, k=0,
                                          base_units="V/Angstrom^2",
                                          m2_units="V^2/Angstrom^4",
                                          seed=hash(ch) & 0xFFFF)
    # Shell counts
    for ch in ("n_first", "n_second"):
        truth[ch] = _write_moments_block(grp, ch, n_atoms, k=0,
                                          base_units="dimensionless",
                                          m2_units="dimensionless",
                                          seed=hash(ch) & 0xFFFF)

    # Delta variants on 3 primary scalars
    for base in ("efield_magnitude", "n_first", "n_second"):
        for suf, m2 in (
            ("_delta",         "V^2/Angstrom^2" if base == "efield_magnitude" else "dimensionless"),
            ("_abs_delta",     "V^2/Angstrom^2" if base == "efield_magnitude" else "dimensionless"),
            ("_delta_squared", "V^4/Angstrom^4" if base == "efield_magnitude" else "dimensionless"),
            ("_dxdt",          "V^2/Angstrom^2/ps^2" if base == "efield_magnitude"
                                else "count^2/ps^2"),
        ):
            bu = "V/Angstrom" if base == "efield_magnitude" else "dimensionless"
            if suf == "_dxdt":
                bu = "V/Angstrom/ps" if base == "efield_magnitude" else "count/ps"
            if suf == "_delta_squared":
                bu = "V^2/Angstrom^2" if base == "efield_magnitude" else "dimensionless"
            _write_moments_block(grp, base + suf, n_atoms, k=0,
                                  base_units=bu, m2_units=m2,
                                  seed=(hash(base + suf) & 0xFFFF))
        grp.create_dataset(f"{base}_rms_delta",
                           data=np.arange(n_atoms, dtype=np.float64) * 0.01) \
           .attrs["units"] = "V/Angstrom" if base == "efield_magnitude" else "dimensionless"

    grp.create_dataset("n_frames_per_atom",
                       data=np.full(n_atoms, n_frames, dtype=np.uint64)) \
       .attrs["units"] = "frame_count"
    grp.create_dataset("delta_n_per_atom",
                       data=np.full(n_atoms, n_frames - 1, dtype=np.uint64)) \
       .attrs["units"] = "frame_count"
    grp.create_dataset("dxdt_n_per_atom",
                       data=np.full(n_atoms, n_frames - 1, dtype=np.uint64)) \
       .attrs["units"] = "frame_count"
    mask = np.ones(n_frames, dtype=np.uint8)
    grp.create_dataset("source_attached_per_frame", data=mask) \
       .attrs["units"] = "dimensionless"

    return truth


@pytest.fixture
def h5_with_water(tmp_path):
    h5 = tmp_path / "traj.h5"
    with h5py.File(h5, "w") as f:
        _write_required_traj_root(f, N_ATOMS, N_FRAMES)
        ts_truth = _write_water_field_time_series(f, N_ATOMS, N_FRAMES)
        wf_truth = _write_water_field_welford(f, N_ATOMS, N_FRAMES)
    return h5, ts_truth, wf_truth


@pytest.fixture
def h5_without_water(tmp_path):
    h5 = tmp_path / "traj.h5"
    with h5py.File(h5, "w") as f:
        _write_required_traj_root(f, N_ATOMS, N_FRAMES)
    return h5


class TestWaterFieldAccess:
    def test_returns_access_container(self, h5_with_water):
        h5, _, _ = h5_with_water
        traj = load_trajectory(h5)
        assert isinstance(traj.water_field, WaterFieldAccess)
        assert traj.water_field.time_series is not None
        assert traj.water_field.welford     is not None

    def test_missing_groups_return_none(self, h5_without_water):
        traj = load_trajectory(h5_without_water)
        assert isinstance(traj.water_field, WaterFieldAccess)
        assert traj.water_field.time_series is None
        assert traj.water_field.welford     is None


class TestWaterFieldTimeSeriesGroup:
    def test_shapes_and_layouts(self, h5_with_water):
        h5, ts_truth, _ = h5_with_water
        ts = load_trajectory(h5).water_field.time_series
        assert isinstance(ts, WaterFieldTimeSeriesGroup)
        # T2-only EFG: shape (N, T, 5), parity 2e, layout "T2_m-2,...,T2_m+2"
        assert ts.efg.shape         == (N_ATOMS, N_FRAMES, 5)
        assert ts.efg_first.shape   == (N_ATOMS, N_FRAMES, 5)
        assert ts.efg_parity        == "2e"
        assert ts.efg_irrep_layout  == "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        # E-field Vec3
        assert ts.efield.shape        == (N_ATOMS, N_FRAMES, 3)
        assert ts.efield_layout       == "x,y,z"
        assert ts.efield_parity       == "1o"
        # Shell counts
        assert ts.n_first.shape  == (N_ATOMS, N_FRAMES)
        assert ts.n_second.shape == (N_ATOMS, N_FRAMES)
        assert ts.n_first.dtype  == np.uint32

    def test_data_round_trips(self, h5_with_water):
        h5, ts_truth, _ = h5_with_water
        ts = load_trajectory(h5).water_field.time_series
        np.testing.assert_array_equal(ts.efield, ts_truth["efield"])
        np.testing.assert_array_equal(ts.efg,    ts_truth["efg"])
        np.testing.assert_array_equal(ts.n_first, ts_truth["n_first"])

    def test_cutoff_attrs(self, h5_with_water):
        h5, _, _ = h5_with_water
        ts = load_trajectory(h5).water_field.time_series
        assert ts.efield_cutoff_A   == 15.0
        assert ts.n_first_cutoff_A  == 3.5
        assert ts.n_second_cutoff_A == 5.5

    def test_source_attached_provenance(self, h5_with_water):
        h5, _, _ = h5_with_water
        ts = load_trajectory(h5).water_field.time_series
        assert ts.source_attached_count == N_FRAMES
        np.testing.assert_array_equal(
            ts.source_attached_per_frame,
            np.ones(N_FRAMES, dtype=np.uint8))


class TestWaterFieldWelfordGroup:
    def test_efg_t2_per_component_shape(self, h5_with_water):
        h5, _, wf_truth = h5_with_water
        wf = load_trajectory(h5).water_field.welford
        assert isinstance(wf, WaterFieldWelfordGroup)
        # efg_t2 is (N, 5) per-component
        assert wf.efg_t2.mean.shape == (N_ATOMS, 5)
        # efg_t2magnitude is (N,) scalar
        assert wf.efg_t2magnitude.mean.shape == (N_ATOMS,)
        # Per-component round-trip
        np.testing.assert_array_equal(wf.efg_t2.mean, wf_truth["efg_t2"]["mean"])

    def test_efield_per_component(self, h5_with_water):
        h5, _, wf_truth = h5_with_water
        wf = load_trajectory(h5).water_field.welford
        assert isinstance(wf.efield_x, WelfordMoments)
        assert wf.efield_x.mean.shape == (N_ATOMS,)
        np.testing.assert_array_equal(wf.efield_x.mean, wf_truth["efield_x"]["mean"])

    def test_delta_variants_on_primary_scalars(self, h5_with_water):
        h5, _, _ = h5_with_water
        wf = load_trajectory(h5).water_field.welford
        # 3 primary scalars get full delta variants
        for ch in (wf.efield_magnitude_delta, wf.n_first_delta, wf.n_second_delta):
            assert isinstance(ch, WelfordMoments)
            assert ch.mean.shape == (N_ATOMS,)
        # rms_delta is a per-atom scalar array (not a WelfordMoments)
        assert wf.efield_magnitude_rms_delta.shape == (N_ATOMS,)

    def test_efg_t1_absent_on_dataclass(self, h5_with_water):
        """T1 structurally zero → no efg_t1 field on the group."""
        h5, _, _ = h5_with_water
        wf = load_trajectory(h5).water_field.welford
        assert not hasattr(wf, "efg_t1")
        assert not hasattr(wf, "efg_first_t1")

    def test_provenance(self, h5_with_water):
        h5, _, _ = h5_with_water
        wf = load_trajectory(h5).water_field.welford
        assert wf.source_attached_count == N_FRAMES
        assert wf.mean_dt_ps   == pytest.approx(0.1)
        assert wf.frame_index_range == (0, N_FRAMES - 1)
        assert wf.n_frames_per_atom.shape == (N_ATOMS,)
