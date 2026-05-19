"""SDK round-trip tests for /trajectory/dihedral_bin_transition/."""

import h5py
import numpy as np
import pytest

from nmr_extract import (
    load_trajectory,
    DihedralBinTransitionGroup,
)


N_ATOMS    = 12
N_RESIDUES = 4
N_FRAMES   = 5


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
                       data=np.zeros((N_ATOMS * N_FRAMES, 3), dtype=np.float64))


def _write_dihedral_bin_transition(f: h5py.File) -> dict:
    grp = f.create_group("/trajectory/dihedral_bin_transition")
    grp.attrs["result_name"] = "DihedralBinTransitionTrajectoryResult"
    grp.attrs["n_residues"]  = N_RESIDUES
    grp.attrs["n_atoms"]     = N_ATOMS
    grp.attrs["n_frames"]    = N_FRAMES
    grp.attrs["finalized"]   = True
    grp.attrs["backbone_bin_count"] = 6
    grp.attrs["chi_count"]          = 4
    grp.attrs["chi_bin_count"]      = 3
    grp.attrs["backbone_bin_legend"] = (
        "0=unassigned, 1=alphaR, 2=beta, 3=alphaL, 4=PPII, 5=other.")
    grp.attrs["backbone_bin_boundaries"] = (
        "Resolution order: alphaR -> alphaL -> PPII -> beta -> other")
    grp.attrs["chi_rotamer_legend"] = (
        "0=gplus, 1=trans, 2=gminus, 255=unassigned")
    grp.attrs["transition_gate"] = (
        "Both prev and curr frame must have an observed bin")
    grp.attrs["angle_convention"] = (
        "IUPAC signed dihedral; BackbonePredecessor/Successor bond-graph")
    grp.attrs["residue_axis"] = "protein_residue_index"
    grp.attrs["source"] = "positions + AminoAcidType.chi_angles"
    grp.attrs["source_attached_policy"] = "always_attached"

    bb_trans = np.array([0, 2, 1, 0], dtype=np.uint32)
    bb_dom   = np.array([0, 1, 1, 5], dtype=np.uint8)  # 0 = unassigned for N-term
    n_obs    = np.array([0, N_FRAMES, N_FRAMES, N_FRAMES], dtype=np.uint32)
    grp.create_dataset("backbone_transition_count", data=bb_trans)
    grp.create_dataset("backbone_dominant_region", data=bb_dom)
    grp.create_dataset("n_frames_observed", data=n_obs)

    bb_occ = np.zeros((N_RESIDUES, 6), dtype=np.uint32)
    bb_occ[1, 1] = 3; bb_occ[1, 2] = 2  # res 1 mostly alphaR with some beta
    bb_occ[2, 1] = 5                      # res 2 all alphaR
    bb_occ[3, 5] = 5                      # res 3 all other
    grp.create_dataset("backbone_bin_occupancy", data=bb_occ)

    chi_trans = np.array([[0, 0, 0, 0],
                          [1, 0, 0, 0],
                          [0, 0, 0, 0],
                          [0, 0, 0, 0]], dtype=np.uint32)
    chi_dom = np.array([[255, 255, 255, 255],
                        [  0,   1, 255, 255],
                        [  0,   0,   0, 255],
                        [255, 255, 255, 255]], dtype=np.uint8)
    chi_n_obs = np.array([[0, 0, 0, 0],
                          [5, 5, 0, 0],
                          [5, 5, 5, 0],
                          [0, 0, 0, 0]], dtype=np.uint32)
    grp.create_dataset("chi_transition_count", data=chi_trans)
    grp.create_dataset("chi_dominant_rotamer", data=chi_dom)
    grp.create_dataset("chi_n_frames_observed", data=chi_n_obs)

    chi_occ = np.zeros((N_RESIDUES, 4, 3), dtype=np.uint32)
    chi_occ[1, 0, 0] = 3; chi_occ[1, 0, 1] = 2  # chi1 of res1 g+/t mix
    chi_occ[1, 1, 1] = 5                          # chi2 of res1 all t
    grp.create_dataset("chi_rotamer_occupancy", data=chi_occ)

    grp.create_dataset("frame_indices",
                       data=np.arange(N_FRAMES, dtype=np.uint64))
    grp.create_dataset("frame_times",
                       data=np.arange(N_FRAMES, dtype=np.float64) * 0.1)
    grp.create_dataset("source_attached_per_frame",
                       data=np.ones((N_FRAMES,), dtype=np.uint8))

    return {"bb_trans": bb_trans, "bb_dom": bb_dom, "n_obs": n_obs,
            "bb_occ": bb_occ, "chi_trans": chi_trans, "chi_dom": chi_dom,
            "chi_n_obs": chi_n_obs, "chi_occ": chi_occ}


def test_bin_transition_round_trip(tmp_path):
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        expected = _write_dihedral_bin_transition(f)

    traj = load_trajectory(h5_path)
    assert traj.dihedral_bin_transitions is not None
    g = traj.dihedral_bin_transitions
    assert isinstance(g, DihedralBinTransitionGroup)

    np.testing.assert_array_equal(g.backbone_transition_count, expected["bb_trans"])
    np.testing.assert_array_equal(g.backbone_dominant_region, expected["bb_dom"])
    np.testing.assert_array_equal(g.n_frames_observed, expected["n_obs"])
    np.testing.assert_array_equal(g.backbone_bin_occupancy, expected["bb_occ"])
    np.testing.assert_array_equal(g.chi_transition_count, expected["chi_trans"])
    np.testing.assert_array_equal(g.chi_dominant_rotamer, expected["chi_dom"])
    np.testing.assert_array_equal(g.chi_n_frames_observed, expected["chi_n_obs"])
    np.testing.assert_array_equal(g.chi_rotamer_occupancy, expected["chi_occ"])

    # Bin-occupancy invariant: per-residue sum across non-unassigned
    # bins equals n_frames_observed.
    for ri in range(N_RESIDUES):
        sum_observed = g.backbone_bin_occupancy[ri, 1:].sum()
        assert sum_observed == g.n_frames_observed[ri]


def test_bin_transition_attrs(tmp_path):
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        _write_dihedral_bin_transition(f)
    g = load_trajectory(h5_path).dihedral_bin_transitions
    assert "alphaR" in g.backbone_bin_legend
    assert "Resolution order" in g.backbone_bin_boundaries
    assert "gplus" in g.chi_rotamer_legend
    assert "Both prev and curr" in g.transition_gate
    assert "BackbonePredecessor" in g.angle_convention
    assert g.source_attached_policy == "always_attached"


def test_bin_transition_missing_group(tmp_path):
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
    traj = load_trajectory(h5_path)
    assert traj.dihedral_bin_transitions is None
