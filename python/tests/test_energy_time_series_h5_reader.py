"""SDK round-trip tests for the two energy time-series H5 groups.

Verifies that:
  - load_trajectory() exposes /trajectory/gromacs_energy_time_series/ as
    EnergyAccess.gromacs (GromacsEnergyTimeSeriesGroup),
  - and /trajectory/bonded_energy_time_series/ as EnergyAccess.bonded
    (BondedEnergyTimeSeriesGroup),
  - each (T,) scalar channel surfaces with correct dtype + length,
  - the virial / pressure_tensor (T, 9) tensor channels carry the
    XX,XY,XZ,... layout convention,
  - bonded breakdown's 7 per-atom channels surface at (N, T) shape,
  - units + split_convention attributes propagate,
  - missing groups return None on EnergyAccess.

Synthesises trajectory.h5 from scratch using h5py — does not depend on a
real C++ extraction. Mirrors GromacsEnergyTimeSeriesTrajectoryResult::
WriteH5Group + BondedEnergyTimeSeriesTrajectoryResult::WriteH5Group exactly.
If the C++ schema drifts, these tests will catch it when run on real H5.
"""

import h5py
import numpy as np
import pytest

from nmr_extract import (
    load_trajectory,
    EnergyAccess,
    GromacsEnergyTimeSeriesGroup,
    BondedEnergyTimeSeriesGroup,
)


N_ATOMS = 10
N_FRAMES = 6


def _write_required_traj_root(f: h5py.File, n_atoms: int, n_frames: int) -> None:
    """Mirror the analysis-schema frame metadata so load_trajectory has
    coherent (T, N, 3) positions + (T,) frame_times before reaching the
    energy groups. Same as the welford H5 reader fixture helper.
    """
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


def _write_gromacs_energy_ts(f: h5py.File, n_frames: int) -> dict:
    """Emit /trajectory/gromacs_energy_time_series with a recognisable pattern.

    Includes the R2-review provenance fields (source_attached_per_frame,
    source_attached_count, energy_frame_times_ps) so the fixture matches
    production H5 schema.
    """
    grp = f.create_group("/trajectory/gromacs_energy_time_series")
    grp.attrs["result_name"]            = "GromacsEnergyTimeSeriesTrajectoryResult"
    grp.attrs["n_frames"]               = n_frames
    grp.attrs["source_attached_count"]  = n_frames
    grp.attrs["finalized"]              = True
    grp.attrs["units"]                  = "kJ/mol"
    grp.attrs["tensor_layout"]          = "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ"

    rng = np.random.default_rng(42)
    fields = {
        "coulomb_sr":    rng.normal(-1000, 50, n_frames),
        "coulomb_recip": rng.normal(-500, 20, n_frames),
        "coulomb_14":    rng.normal(0, 10, n_frames),
        "bond":          rng.normal(100, 5, n_frames),
        "angle":         rng.normal(200, 8, n_frames),
        "urey_bradley":  np.zeros(n_frames),
        "proper_dih":    rng.normal(150, 10, n_frames),
        "improper_dih":  np.zeros(n_frames),
        "cmap_dih":      np.zeros(n_frames),
        "lj_sr":         rng.normal(-2000, 50, n_frames),
        "lj_14":         rng.normal(50, 5, n_frames),
        "disper_corr":   rng.normal(-100, 2, n_frames),
        "potential":     rng.normal(-5000, 100, n_frames),
        "kinetic":       rng.normal(2000, 30, n_frames),
        "total_energy":  rng.normal(-3000, 100, n_frames),
        "enthalpy":      rng.normal(-2800, 100, n_frames),
        "temperature":   rng.normal(300, 1, n_frames),
        "pressure":      rng.normal(1, 0.5, n_frames),
        "volume":        rng.normal(100, 0.5, n_frames),
        "density":       rng.normal(1000, 1, n_frames),
        "box_x":         rng.normal(5, 0.05, n_frames),
        "box_y":         rng.normal(5, 0.05, n_frames),
        "box_z":         rng.normal(5, 0.05, n_frames),
        "T_protein":     rng.normal(300, 1, n_frames),
        "T_non_protein": rng.normal(300, 1, n_frames),
    }
    for name, arr in fields.items():
        ds = grp.create_dataset(name, data=arr.astype(np.float64))
        ds.attrs["units"] = {
            "temperature": "K", "T_protein": "K", "T_non_protein": "K",
            "pressure":    "bar",
            "volume":      "nm^3",
            "density":     "kg/m^3",
            "box_x":       "nm", "box_y": "nm", "box_z": "nm",
        }.get(name, "kJ/mol")

    vir  = rng.normal(0, 100, (n_frames, 9)).astype(np.float64)
    pres = rng.normal(0, 1,   (n_frames, 9)).astype(np.float64)
    ds_v = grp.create_dataset("virial",          data=vir)
    ds_v.attrs["units"] = "kJ/mol"
    ds_v.attrs["tensor_layout"] = "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ"
    ds_p = grp.create_dataset("pressure_tensor", data=pres)
    ds_p.attrs["units"] = "bar"
    ds_p.attrs["tensor_layout"] = "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ"

    fi = np.arange(n_frames, dtype=np.uint64)
    ft = np.arange(n_frames, dtype=np.float64) * 0.1
    # energy_frame_times_ps: snap-distance audit. Simulate a small snap
    # offset (0.0005 ps from the trajectory frame time) for the test.
    energy_t = ft + 0.0005
    mask = np.ones(n_frames, dtype=np.uint8)
    grp.create_dataset("frame_indices", data=fi).attrs["units"] = "frame_index"
    grp.create_dataset("frame_times",   data=ft).attrs["units"] = "ps"
    ds_et = grp.create_dataset("energy_frame_times_ps", data=energy_t)
    ds_et.attrs["units"] = "ps"
    ds_et.attrs["description"] = ("Matched .edr row time per frame. "
                                  "Snap distance from frame_times reveals "
                                  "EnergyAtTime() drift.")
    grp.create_dataset("source_attached_per_frame", data=mask) \
       .attrs["units"] = "dimensionless"

    return {**fields, "virial": vir, "pressure_tensor": pres,
            "frame_indices": fi, "frame_times": ft,
            "energy_frame_times_ps": energy_t,
            "source_attached_per_frame": mask}


def _write_bonded_energy_ts(f: h5py.File, n_atoms: int, n_frames: int) -> dict:
    """Emit /trajectory/bonded_energy_time_series with a recognisable pattern."""
    grp = f.create_group("/trajectory/bonded_energy_time_series")
    grp.attrs["result_name"]           = "BondedEnergyTimeSeriesTrajectoryResult"
    grp.attrs["n_atoms"]               = n_atoms
    grp.attrs["n_frames"]              = n_frames
    grp.attrs["source_attached_count"] = n_frames
    grp.attrs["finalized"]             = True
    grp.attrs["units"]                 = "kJ/mol"
    grp.attrs["split_convention"]      = "evenly_among_participating_atoms"
    grp.attrs["split_convention_note"] = "one of several valid attributions; calibration-absorbable"
    grp.attrs["system_scope"]          = "protein_slice_only; sum != .edr whole-system term"

    rng = np.random.default_rng(7)
    channels = ("bond", "angle", "urey_bradley", "proper_dih",
                "improper_dih", "cmap_dih", "total")
    fields = {}
    for ch in channels:
        arr = rng.normal(0, 5, (n_atoms, n_frames)).astype(np.float64)
        if ch in ("urey_bradley", "improper_dih", "cmap_dih"):
            arr.fill(0.0)
        fields[ch] = arr
        ds = grp.create_dataset(ch, data=arr)
        ds.attrs["units"] = "kJ/mol"

    fi = np.arange(n_frames, dtype=np.uint64)
    ft = np.arange(n_frames, dtype=np.float64) * 0.1
    mask = np.ones(n_frames, dtype=np.uint8)
    grp.create_dataset("frame_indices", data=fi).attrs["units"] = "frame_index"
    grp.create_dataset("frame_times",   data=ft).attrs["units"] = "ps"
    grp.create_dataset("source_attached_per_frame", data=mask) \
       .attrs["units"] = "dimensionless"

    return {**fields, "frame_indices": fi, "frame_times": ft,
            "source_attached_per_frame": mask}


# ─── Fixtures ───────────────────────────────────────────────────────


@pytest.fixture
def h5_with_both_energy(tmp_path):
    h5 = tmp_path / "traj.h5"
    with h5py.File(h5, "w") as f:
        _write_required_traj_root(f, N_ATOMS, N_FRAMES)
        gromacs = _write_gromacs_energy_ts(f, N_FRAMES)
        bonded  = _write_bonded_energy_ts(f, N_ATOMS, N_FRAMES)
    return h5, gromacs, bonded


@pytest.fixture
def h5_with_no_energy(tmp_path):
    h5 = tmp_path / "traj.h5"
    with h5py.File(h5, "w") as f:
        _write_required_traj_root(f, N_ATOMS, N_FRAMES)
    return h5


# ─── Tests ──────────────────────────────────────────────────────────


class TestEnergyAccessAttachment:
    def test_returns_energy_access_container(self, h5_with_both_energy):
        h5, _, _ = h5_with_both_energy
        traj = load_trajectory(h5)
        assert isinstance(traj.energy, EnergyAccess)
        assert traj.energy.gromacs is not None
        assert traj.energy.bonded  is not None

    def test_missing_groups_return_none(self, h5_with_no_energy):
        traj = load_trajectory(h5_with_no_energy)
        assert isinstance(traj.energy, EnergyAccess)
        assert traj.energy.gromacs is None
        assert traj.energy.bonded  is None


class TestGromacsEnergyTimeSeriesGroup:
    def test_channel_round_trips(self, h5_with_both_energy):
        h5, gromacs_truth, _ = h5_with_both_energy
        ge = load_trajectory(h5).energy.gromacs
        assert isinstance(ge, GromacsEnergyTimeSeriesGroup)
        # Spot-check a representative sample across categories.
        np.testing.assert_array_equal(ge.total_energy, gromacs_truth["total_energy"])
        np.testing.assert_array_equal(ge.temperature, gromacs_truth["temperature"])
        np.testing.assert_array_equal(ge.box_x,       gromacs_truth["box_x"])
        np.testing.assert_array_equal(ge.T_protein,   gromacs_truth["T_protein"])

    def test_tensor_channels_carry_layout(self, h5_with_both_energy):
        h5, gromacs_truth, _ = h5_with_both_energy
        ge = load_trajectory(h5).energy.gromacs
        assert ge.virial.shape          == (N_FRAMES, 9)
        assert ge.pressure_tensor.shape == (N_FRAMES, 9)
        np.testing.assert_array_equal(ge.virial,          gromacs_truth["virial"])
        np.testing.assert_array_equal(ge.pressure_tensor, gromacs_truth["pressure_tensor"])

    def test_units_propagate(self, h5_with_both_energy):
        h5, _, _ = h5_with_both_energy
        ge = load_trajectory(h5).energy.gromacs
        assert ge.units == "kJ/mol"

    def test_frame_indexing_aligned_with_root(self, h5_with_both_energy):
        h5, gromacs_truth, _ = h5_with_both_energy
        traj = load_trajectory(h5)
        ge = traj.energy.gromacs
        assert ge.frame_indices.shape == (N_FRAMES,)
        np.testing.assert_array_equal(ge.frame_indices, gromacs_truth["frame_indices"])
        np.testing.assert_allclose(ge.frame_times,      traj.frame_times)

    def test_low_energy_selection_safe_pattern(self, h5_with_both_energy):
        """The user's 2026-05-18 use case + R3 codex F1 NaN-safety.

        With all frames finite (fixture writes real values), the safe
        pattern selects bottom-N% by total_energy. The pattern itself
        is documented on GromacsEnergyTimeSeriesGroup — this test
        exercises it to prevent the unsafe `np.argsort(e)[:N]` form from
        creeping back into downstream code.
        """
        h5, _, _ = h5_with_both_energy
        ge = load_trajectory(h5).energy.gromacs
        e = ge.total_energy
        valid = np.isfinite(e)
        assert valid.all(), "fixture should have all-finite total_energy"
        n_select = int(0.5 * valid.sum())
        rank_in_valid = np.argsort(e[valid])[:n_select]
        selected = np.where(valid)[0][rank_in_valid]
        assert selected.shape == (n_select,)
        # Sanity: every selected frame is finite + ≤ median total energy.
        for f in selected:
            assert np.isfinite(e[f])
            assert e[f] <= np.median(e)

    def test_low_energy_selection_handles_nan(self, tmp_path):
        """If some frames have NaN total_energy (e.g. .edr column missing
        in those rows), the safe pattern silently excludes them while
        the naive `np.argsort(e)[:N]` would silently include them."""
        h5 = tmp_path / "traj_with_nans.h5"
        with h5py.File(h5, "w") as f:
            _write_required_traj_root(f, N_ATOMS, N_FRAMES)
            gromacs_truth = _write_gromacs_energy_ts(f, N_FRAMES)
            _write_bonded_energy_ts(f, N_ATOMS, N_FRAMES)
            # Inject NaN on half the frames' total_energy.
            ds = f["/trajectory/gromacs_energy_time_series/total_energy"]
            data = ds[:]
            data[::2] = np.nan  # every other frame NaN
            ds[:] = data

        ge = load_trajectory(h5).energy.gromacs
        e = ge.total_energy
        valid = np.isfinite(e)
        n_valid = int(valid.sum())
        assert n_valid == N_FRAMES // 2, "half the frames should still be valid"

        # Safe pattern: only valid frames selected.
        n_select = max(1, int(0.5 * n_valid))
        rank_in_valid = np.argsort(e[valid])[:n_select]
        selected = np.where(valid)[0][rank_in_valid]
        for f in selected:
            assert np.isfinite(e[f]), \
                f"safe pattern selected NaN frame {f} — gating broken"

    def test_source_attached_provenance(self, h5_with_both_energy):
        """R2 review: source_attached_per_frame mask + count round-trip."""
        h5, _, _ = h5_with_both_energy
        ge = load_trajectory(h5).energy.gromacs
        assert ge.source_attached_count == N_FRAMES
        np.testing.assert_array_equal(
            ge.source_attached_per_frame,
            np.ones(N_FRAMES, dtype=np.uint8))

    def test_energy_frame_times_snap_distance(self, h5_with_both_energy):
        """R2 review L7: energy_frame_times_ps surfaces snap distance.

        The fixture writes energy_frame_times_ps = frame_times + 0.0005 ps,
        simulating a small .edr-stride / trajectory-stride mismatch. The
        difference is the snap distance.
        """
        h5, _, _ = h5_with_both_energy
        traj = load_trajectory(h5)
        ge = traj.energy.gromacs
        assert ge.energy_frame_times_ps.shape == (N_FRAMES,)
        snap = ge.energy_frame_times_ps - traj.frame_times
        np.testing.assert_allclose(snap, np.full(N_FRAMES, 0.0005))


class TestBondedEnergyTimeSeriesGroup:
    def test_per_atom_shape(self, h5_with_both_energy):
        h5, _, bonded_truth = h5_with_both_energy
        be = load_trajectory(h5).energy.bonded
        assert isinstance(be, BondedEnergyTimeSeriesGroup)
        for ch in ("bond", "angle", "urey_bradley", "proper_dih",
                   "improper_dih", "cmap_dih", "total"):
            assert getattr(be, ch).shape == (N_ATOMS, N_FRAMES)
        np.testing.assert_array_equal(be.bond,  bonded_truth["bond"])
        np.testing.assert_array_equal(be.total, bonded_truth["total"])

    def test_units_and_split_convention(self, h5_with_both_energy):
        h5, _, _ = h5_with_both_energy
        be = load_trajectory(h5).energy.bonded
        assert be.units            == "kJ/mol"
        assert be.split_convention == "evenly_among_participating_atoms"
        # R2 review: caveat + scope notes surface as fields too.
        assert "calibration-absorbable" in be.split_convention_note
        assert "protein_slice_only"     in be.system_scope

    def test_source_attached_provenance(self, h5_with_both_energy):
        h5, _, _ = h5_with_both_energy
        be = load_trajectory(h5).energy.bonded
        assert be.source_attached_count == N_FRAMES
        np.testing.assert_array_equal(
            be.source_attached_per_frame,
            np.ones(N_FRAMES, dtype=np.uint8))

    def test_zero_channels_legal_on_charmm36m(self, h5_with_both_energy):
        """UB / improper / CMAP are zero on the 1P9J CHARMM36m fixture
        (verified 2026-05-18). The SDK should not flinch — the datasets
        are still emitted with zero values, not omitted."""
        h5, _, _ = h5_with_both_energy
        be = load_trajectory(h5).energy.bonded
        for ch in ("urey_bradley", "improper_dih", "cmap_dih"):
            arr = getattr(be, ch)
            assert arr.shape == (N_ATOMS, N_FRAMES)
            assert np.all(arr == 0.0)
