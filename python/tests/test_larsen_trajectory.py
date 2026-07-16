from dataclasses import fields

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
    tensor_attrs = {
        "n_atoms": 2,
        "n_frames": 3,
        "finalized": True,
        "source_attached_count": np.uint64(2),
        "source_attached_policy": "conditional_larsen_grid_source",
        "atom_axis": "protein_atom_index",
        "frame_axis": "trajectory_frame_row",
        "irrep_layout": (
            "PackFull9: [T0, T1_cartesian_xyz, "
            "T2_real_tesseral_m-2..m+2]"),
        "normalization": "isometric_real_sph",
        "parity": "mixed",
        "coordinate_frame": "conformation_cartesian_xyz",
        "tensor_basis": "project_native_full9_spherical_tensor_v1",
        "tensor_component_order": (
            "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"),
        "tensor_frame": "conformation_cartesian_xyz",
        "tensor_t1_semantics": (
            "Cartesian Levi-Civita dual x,y,z (not real-Y1m): "
            "a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2); "
            "axial a'=det(R) R a; generically nonzero"),
        "tensor_t1_structural_zero": False,
        "tensor_structural_zero_components": "none",
        "e3nn_export": (
            "explicit project-basis to e3nn conversion required before use"),
        "normalization_scope": (
            "xyz tensor payload: T2 uses isometric real-tesseral "
            "normalization; T1 uses the tensor_t1_semantics convention"),
        "transformation": (
            "even_rank2 under proper rotations: T'=R T R^T; signed-rho "
            "DFT-grid lookup is chirality-conditioned and has no "
            "improper-transform contract"),
        "units": "ppm",
    }
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
            g.attrs["result_name"] = f"Synthetic{name}Result"
            for attr_name, value in tensor_attrs.items():
                g.attrs[attr_name] = value
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

    expected_fields = {
        "xyz", "frame_indices", "frame_times",
        "source_attached_per_frame", "result_name", *tensor_attrs,
    }
    for name, timeline in zip(tensor_paths, (
            traj.larsen_hbond.hbond_1pHB,
            traj.larsen_hbond.hbond_2pHB,
            traj.larsen_hbond.hbond_1pHaB,
            traj.larsen_hbond.hbond_2pHaB,
    )):
        assert {item.name for item in fields(timeline)} == expected_fields
        np.testing.assert_array_equal(timeline.frame_indices, frame_indices)
        np.testing.assert_array_equal(timeline.frame_times, frame_times)
        np.testing.assert_array_equal(timeline.source_attached_per_frame, mask)
        assert timeline.result_name == f"Synthetic{name}Result"
        for attr_name, value in tensor_attrs.items():
            assert getattr(timeline, attr_name) == value
