import h5py
import numpy as np

from nmr_extract import load_trajectory


def test_larsen_six_group_loader_and_source_count(tmp_path):
    path = tmp_path / "larsen.h5"
    mask = np.array([1, 0, 1], dtype=np.uint8)
    frame_indices = np.array([0, 2, 4], dtype=np.uint64)
    frame_times = np.array([0.0, 2.0, 4.0])
    tensor_paths = (
        "larsen_hbond_1pHB_shielding_time_series",
        "larsen_hbond_2pHB_shielding_time_series",
        "larsen_hbond_1pHaB_shielding_time_series",
        "larsen_hbond_2pHaB_shielding_time_series",
    )
    with h5py.File(path, "w") as f:
        f.attrs["protein_id"] = "synthetic"
        f.attrs["n_atoms"] = 2
        f.attrs["finalized"] = True
        trajectory = f.create_group("trajectory")
        frames = trajectory.create_group("frames")
        frames.attrs["n_frames"] = 3
        frames.create_dataset("time_ps", data=frame_times)
        frames.create_dataset("original_index", data=frame_indices)
        positions = trajectory.create_group("positions")
        positions.attrs["n_atoms"] = 2
        positions.attrs["n_frames"] = 3
        positions.attrs["result_name"] = \
            "PositionsTimeSeriesTrajectoryResult"
        positions.attrs["finalized"] = True
        positions.create_dataset("xyz", data=np.zeros((2, 3, 3)))
        positions.create_dataset("frame_indices", data=frame_indices)
        positions.create_dataset("frame_times", data=frame_times)
        for name in tensor_paths:
            g = f.create_group(f"trajectory/{name}")
            g.create_dataset("xyz", data=np.zeros((2, 3, 9)))
            g.create_dataset("frame_indices", data=frame_indices)
            g.create_dataset("frame_times", data=frame_times)
            g.create_dataset("source_attached_per_frame", data=mask)
            g.attrs["source_attached_count"] = np.uint64(2)
            g.attrs["source_attached_policy"] = (
                "conditional_larsen_grid_source")
            g.attrs["atom_axis"] = "protein_atom_index"
            g.attrs["frame_axis"] = "trajectory_frame_row"
            g.attrs["irrep_layout"] = (
                "PackFull9: [T0, T1_cartesian_xyz, "
                "T2_real_tesseral_m-2..m+2]")
        for name, dataset, dtype in (
            ("larsen_hbond_water_term_time_series", "water_term", float),
            ("larsen_hbond_count_time_series", "count", np.int32),
        ):
            g = f.create_group(f"trajectory/{name}")
            g.create_dataset(dataset, data=np.zeros((2, 3), dtype=dtype))
            g.create_dataset("frame_indices", data=frame_indices)
            g.create_dataset("frame_times", data=frame_times)
            g.create_dataset("source_attached_per_frame", data=mask)
            g.attrs["source_attached_count"] = np.uint64(2)
            g.attrs["source_attached_policy"] = (
                "conditional_larsen_grid_source")
            g.attrs["atom_axis"] = "protein_atom_index"
            g.attrs["frame_axis"] = "trajectory_frame_row"

    traj = load_trajectory(path)
    assert traj.larsen_hbond is not None
    assert traj.larsen_hbond.hbond_1pHB.xyz.shape == (2, 3, 9)
    assert traj.larsen_hbond.hbond_2pHaB.source_attached_count == 2
    assert traj.larsen_hbond.hbond_2pHaB.source_attached_count == int(mask.sum())
    assert traj.larsen_hbond.water_term.values.shape == (2, 3)
    assert traj.larsen_hbond.count.values.dtype == np.int32
    assert traj.larsen_hbond.count.atom_axis == "protein_atom_index"
