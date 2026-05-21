"""SDK round-trip tests for /trajectory/ring_neighbourhood_trajectory_stats/.

Verifies that:
  - load_trajectory() exposes the group as
    TrajectoryData.ring_neighbourhood_trajectory_stats
    (RingNeighbourhoodTrajectoryStatsGroup);
  - the 4D `geometry` (N, T, R_max, 4) dataset round-trips with correct
    shape/dtype + channel-wise accessor helpers (distance(), rho(), z(),
    in_plane_angle());
  - static `ring_membership_per_atom` (N, R_max) int32 round-trips
    correctly including -1 sentinel for unfilled slots;
  - `live_slots_mask()` returns the right boolean mask;
  - convention attrs (channel_layout, units, source_attached_policy,
    static_snapshot_origin, etc.) propagate verbatim;
  - missing group (no aromatic-ring/atom pairs OR TR not registered)
    returns None.

Synthesises trajectory.h5 from scratch with h5py — mirrors
RingNeighbourhoodTrajectoryStats::WriteH5Group exactly.
"""

import h5py
import numpy as np
import pytest

from nmr_extract import (
    load_trajectory,
    RingNeighbourhoodTrajectoryStatsGroup,
)


N_ATOMS    = 12
N_FRAMES   = 4
R_MAX      = 3


def _write_required_traj_root(f: h5py.File) -> None:
    f.attrs["protein_id"] = "TEST_PROTEIN"
    f.attrs["n_atoms"]    = N_ATOMS
    f.attrs["finalized"]  = True

    traj = f.create_group("trajectory")
    frames = traj.create_group("frames")
    frames.attrs["n_frames"] = N_FRAMES
    times = np.arange(N_FRAMES, dtype=np.float64) * 0.1
    frames.create_dataset("time_ps", data=times)
    frames.create_dataset("original_index",
                          data=np.arange(N_FRAMES, dtype=np.uint64))

    pos = traj.create_group("positions")
    pos.attrs["n_atoms"]     = N_ATOMS
    pos.attrs["n_frames"]    = N_FRAMES
    pos.attrs["result_name"] = "PositionsTimeSeriesTrajectoryResult"
    pos.attrs["finalized"]   = True
    pos.create_dataset("xyz",
                       data=np.zeros((N_ATOMS * N_FRAMES, 3),
                                      dtype=np.float64))


def _write_ring_neighbourhood(f: h5py.File) -> dict:
    """Emit /trajectory/ring_neighbourhood_trajectory_stats with a
    recognisable pattern. Atoms 0..7 carry 2 rings (slots 0, 1, sentinel
    in slot 2); atoms 8..11 carry 1 ring (sentinel in slots 1, 2).
    """
    grp = f.create_group("/trajectory/ring_neighbourhood_trajectory_stats")

    grp.attrs["result_name"] = "RingNeighbourhoodTrajectoryStats"
    grp.attrs["n_atoms"]   = N_ATOMS
    grp.attrs["n_frames"]  = N_FRAMES
    grp.attrs["r_per_atom_max"] = R_MAX
    grp.attrs["finalized"] = True
    grp.attrs["ring_current_spatial_cutoff_A"] = 15.0
    grp.attrs["channel_layout"] = "distance,rho,z,in_plane_angle"
    grp.attrs["units"] = "Angstrom,Angstrom,Angstrom,radians"
    grp.attrs["in_plane_angle_range"] = "[0, 2*pi) ; NaN on ring axis"
    grp.attrs["z_sign_convention"] = "z > 0 = same side as ring normal"
    grp.attrs["nan_semantics"] = (
        "padded slots NaN; in-plane angle NaN at singular cases")
    grp.attrs["aromatic_only"] = "true"
    grp.attrs["source_attached_policy"] = "always_attached"
    grp.attrs["static_snapshot_origin"] = "frame_0"

    # ring_membership_per_atom (N, R_max) int32 with -1 sentinel
    membership = np.full((N_ATOMS, R_MAX), -1, dtype=np.int32)
    # Atoms 0..7 see rings 0 and 1
    for i in range(8):
        membership[i, 0] = 0
        membership[i, 1] = 1
    # Atoms 8..11 see ring 0 only
    for i in range(8, N_ATOMS):
        membership[i, 0] = 0
    grp.create_dataset("ring_membership_per_atom", data=membership)

    # geometry (N, T, R_max, 4) float64
    # Live slots: distance = atom_idx * 0.1 + ring_slot * 0.01 + frame * 0.001
    # rho, z, in_plane_angle = related patterns
    # Padded slots: NaN
    geom = np.full((N_ATOMS, N_FRAMES, R_MAX, 4), np.nan, dtype=np.float64)
    for i in range(N_ATOMS):
        for t in range(N_FRAMES):
            for r in range(R_MAX):
                if membership[i, r] == -1:
                    continue  # leave as NaN
                base = i * 0.1 + r * 0.01 + t * 0.001
                geom[i, t, r, 0] = base                # distance
                geom[i, t, r, 1] = base * 0.5          # rho
                geom[i, t, r, 2] = base * 0.5          # z (so dist^2 = rho^2 + z^2)
                geom[i, t, r, 3] = (r * 1.0) % (2.0 * np.pi)  # in_plane_angle
    grp.create_dataset("geometry", data=geom)

    # Per-frame metadata
    frame_indices = np.arange(N_FRAMES, dtype=np.uint64)
    frame_times = np.arange(N_FRAMES, dtype=np.float64) * 0.1
    source_attached = np.ones((N_FRAMES,), dtype=np.uint8)
    grp.create_dataset("frame_indices", data=frame_indices)
    grp.create_dataset("frame_times", data=frame_times)
    grp.create_dataset("source_attached_per_frame", data=source_attached)

    return {
        "membership": membership,
        "geometry": geom,
        "frame_indices": frame_indices,
        "frame_times": frame_times,
    }


def test_ring_neighbourhood_round_trip(tmp_path):
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        expected = _write_ring_neighbourhood(f)

    traj = load_trajectory(h5_path)
    assert traj.ring_neighbourhood_trajectory_stats is not None
    g = traj.ring_neighbourhood_trajectory_stats
    assert isinstance(g, RingNeighbourhoodTrajectoryStatsGroup)

    # Shape + dtype contracts
    assert g.geometry.shape == (N_ATOMS, N_FRAMES, R_MAX, 4)
    assert g.geometry.dtype == np.float64
    assert g.ring_membership_per_atom.shape == (N_ATOMS, R_MAX)
    assert g.ring_membership_per_atom.dtype == np.int32

    # Membership round-trip
    np.testing.assert_array_equal(
        g.ring_membership_per_atom, expected["membership"])

    # Geometry round-trip (NaN-aware: equal_nan via numpy.testing)
    np.testing.assert_array_equal(
        g.geometry, expected["geometry"], strict=False)


def test_ring_neighbourhood_channel_accessors(tmp_path):
    """The .distance() / .rho() / .z() / .in_plane_angle() helpers
    expose the right channel slices."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        _write_ring_neighbourhood(f)
    traj = load_trajectory(h5_path)
    g = traj.ring_neighbourhood_trajectory_stats

    assert g.distance().shape == (N_ATOMS, N_FRAMES, R_MAX)
    assert g.rho().shape == (N_ATOMS, N_FRAMES, R_MAX)
    assert g.z().shape == (N_ATOMS, N_FRAMES, R_MAX)
    assert g.in_plane_angle().shape == (N_ATOMS, N_FRAMES, R_MAX)

    # Channel 0/1/2/3 vs full geometry
    np.testing.assert_array_equal(g.distance(), g.geometry[..., 0],
                                   strict=False)
    np.testing.assert_array_equal(g.rho(), g.geometry[..., 1], strict=False)
    np.testing.assert_array_equal(g.z(), g.geometry[..., 2], strict=False)
    np.testing.assert_array_equal(g.in_plane_angle(),
                                   g.geometry[..., 3], strict=False)


def test_ring_neighbourhood_live_mask(tmp_path):
    """live_slots_mask() identifies populated (atom, ring) pairs."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        _write_ring_neighbourhood(f)
    traj = load_trajectory(h5_path)
    g = traj.ring_neighbourhood_trajectory_stats

    mask = g.live_slots_mask()
    assert mask.shape == (N_ATOMS, R_MAX)
    assert mask.dtype == bool
    # Atoms 0..7 have 2 live slots each = 16 live; atoms 8..11 have 1 = 4
    # Total live = 20
    assert mask.sum() == 20
    # Padded slots (membership == -1) MUST be False
    assert not mask[0, 2]  # atom 0 slot 2 is padded
    assert not mask[8, 1]  # atom 8 slot 1 is padded
    # Live slots MUST be True
    assert mask[0, 0]
    assert mask[8, 0]


def test_ring_neighbourhood_padded_slots_are_nan(tmp_path):
    """Geometry values in padded (sentinel) slots must be NaN across all
    frames + channels. Reader contract: use isfinite() to filter."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        _write_ring_neighbourhood(f)
    traj = load_trajectory(h5_path)
    g = traj.ring_neighbourhood_trajectory_stats

    mask = g.live_slots_mask()
    for i in range(N_ATOMS):
        for r in range(R_MAX):
            if mask[i, r]:
                continue
            for t in range(N_FRAMES):
                for ch in range(4):
                    assert np.isnan(g.geometry[i, t, r, ch]), (
                        f"padded slot ({i}, {t}, {r}, {ch}) must be NaN, "
                        f"got {g.geometry[i, t, r, ch]}")


def test_ring_neighbourhood_convention_attrs(tmp_path):
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        _write_ring_neighbourhood(f)
    traj = load_trajectory(h5_path)
    g = traj.ring_neighbourhood_trajectory_stats

    assert g.channel_layout == "distance,rho,z,in_plane_angle"
    assert "Angstrom" in g.units
    assert "radians" in g.units
    assert g.aromatic_only.startswith("true")
    assert g.source_attached_policy.startswith("always_attached")
    assert g.static_snapshot_origin.startswith("frame_0")
    assert g.ring_current_spatial_cutoff_A == 15.0
    assert g.n_atoms == N_ATOMS
    assert g.n_frames == N_FRAMES
    assert g.r_per_atom_max == R_MAX


def test_ring_neighbourhood_absent(tmp_path):
    """No group emitted => SDK returns None (graceful absence handling).
    Matches the C++ side's R_max == 0 skip-the-group contract."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        # Deliberately do NOT emit the ring_neighbourhood group
    traj = load_trajectory(h5_path)
    assert traj.ring_neighbourhood_trajectory_stats is None
