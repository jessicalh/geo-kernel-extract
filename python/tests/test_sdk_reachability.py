"""Round-trip coverage for producer data that was previously unreachable."""

from dataclasses import fields

import h5py
import numpy as np

from nmr_extract import (
    AIMNet2ChargeTimeSeriesGroup,
    BondLengthStatsGroup,
    BsShieldingTimeSeriesGroup,
    BsT0AutocorrelationGroup,
    DsspGroup,
    GeometryGroup,
    HmShieldingTimeSeriesGroup,
    LocalBackboneGeometryNpyGroup,
    McConnellShieldingTimeSeriesGroup,
    PositionsTimeSeriesGroup,
    SasaTimeSeriesGroup,
    TrajectoryAtomsGroup,
    VectorField,
    load,
    load_trajectory,
)
from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
)


def _save(path, stem, data):
    data = np.asarray(data)
    np.save(path / f"{stem}.npy", data)
    return data


def test_all_previously_dropped_npy_arrays_reach_their_owner(tmp_path):
    n_atoms = 4
    n_residues = 2
    n_bonds = 2
    write_required_sdk_npys(
        tmp_path, n_atoms=n_atoms, n_residues=n_residues, n_bonds=n_bonds)
    write_minimal_topology_sidecar(
        tmp_path, n_atoms=n_atoms, n_residues=n_residues, n_bonds=n_bonds)

    expected = {
        "bond_length": _save(tmp_path, "bond_length", [1.1, 1.4]),
        "bond_direction": _save(
            tmp_path, "bond_direction", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        "bond_geometry_valid": _save(
            tmp_path, "bond_geometry_valid", np.array([1, 0], dtype=np.uint8)),
        "tau_N_CA_C": _save(tmp_path, "tau_N_CA_C", [0.11, 0.12]),
        "tau_N_CA_C_valid": _save(
            tmp_path, "tau_N_CA_C_valid", np.array([1, 0], dtype=np.uint8)),
        "angle_N_CA_CB": _save(tmp_path, "angle_N_CA_CB", [0.21, 0.22]),
        "angle_N_CA_CB_valid": _save(
            tmp_path, "angle_N_CA_CB_valid",
            np.array([1, 1], dtype=np.uint8)),
        "angle_CB_CA_C": _save(tmp_path, "angle_CB_CA_C", [0.31, 0.32]),
        "angle_CB_CA_C_valid": _save(
            tmp_path, "angle_CB_CA_C_valid",
            np.array([0, 1], dtype=np.uint8)),
        "angle_Cprev_N_CA": _save(
            tmp_path, "angle_Cprev_N_CA", [0.41, 0.42]),
        "angle_Cprev_N_CA_valid": _save(
            tmp_path, "angle_Cprev_N_CA_valid",
            np.array([1, 0], dtype=np.uint8)),
        "angle_CA_C_Nnext": _save(
            tmp_path, "angle_CA_C_Nnext", [0.51, 0.52]),
        "angle_CA_C_Nnext_valid": _save(
            tmp_path, "angle_CA_C_Nnext_valid",
            np.array([0, 1], dtype=np.uint8)),
        "cb_deviation": _save(tmp_path, "cb_deviation", [0.61, 0.62]),
        "cb_deviation_valid": _save(
            tmp_path, "cb_deviation_valid", np.array([1, 1], dtype=np.uint8)),
        "cb_residual_vector": _save(
            tmp_path, "cb_residual_vector",
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
        "cb_residual_vector_valid": _save(
            tmp_path, "cb_residual_vector_valid",
            np.array([1, 0], dtype=np.uint8)),
        "hbond_pairs_index": _save(
            tmp_path, "hbond_pairs_index",
            np.arange(18, dtype=np.int32).reshape(3, 6)),
        "hbond_pairs_geometry": _save(
            tmp_path, "hbond_pairs_geometry",
            np.arange(15, dtype=np.float64).reshape(3, 5) / 10.0),
        "hbond_pairs_angle_valid": _save(
            tmp_path, "hbond_pairs_angle_valid",
            np.array([1, 0, 1], dtype=np.uint8)),
        "dssp_torsion_angle": _save(
            tmp_path, "dssp_torsion_angle",
            np.arange(12, dtype=np.float64).reshape(2, 6) / 10.0),
        "dssp_torsion_sin": _save(
            tmp_path, "dssp_torsion_sin",
            np.arange(12, dtype=np.float64).reshape(2, 6) / 20.0),
        "dssp_torsion_cos": _save(
            tmp_path, "dssp_torsion_cos",
            np.arange(12, dtype=np.float64).reshape(2, 6) / 30.0),
        "dssp_torsion_valid": _save(
            tmp_path, "dssp_torsion_valid",
            np.array([[1, 1, 1, 0, 0, 0], [1, 1, 1, 1, 0, 0]],
                     dtype=np.uint8)),
        "dssp_hbond_partner_residue_index": _save(
            tmp_path, "dssp_hbond_partner_residue_index",
            np.arange(16, dtype=np.int32).reshape(4, 4)),
        "omega_sin": _save(tmp_path, "omega_sin", [0.71, 0.72]),
        "omega_cos": _save(tmp_path, "omega_cos", [0.81, 0.82]),
        "omega_valid": _save(
            tmp_path, "omega_valid", np.array([1, 0], dtype=np.uint8)),
    }

    protein = load(tmp_path)

    assert isinstance(protein.geometry, GeometryGroup)
    np.testing.assert_array_equal(
        protein.geometry.bond_length, expected["bond_length"])
    assert isinstance(protein.geometry.bond_direction, VectorField)
    np.testing.assert_array_equal(
        protein.geometry.bond_direction.data, expected["bond_direction"])
    np.testing.assert_array_equal(
        protein.geometry.bond_geometry_valid,
        expected["bond_geometry_valid"])

    local = protein.local_backbone_geometry
    assert isinstance(local, LocalBackboneGeometryNpyGroup)
    for stem in (
        "tau_N_CA_C", "tau_N_CA_C_valid",
        "angle_N_CA_CB", "angle_N_CA_CB_valid",
        "angle_CB_CA_C", "angle_CB_CA_C_valid",
        "angle_Cprev_N_CA", "angle_Cprev_N_CA_valid",
        "angle_CA_C_Nnext", "angle_CA_C_Nnext_valid",
        "cb_deviation", "cb_deviation_valid", "cb_residual_vector",
        "cb_residual_vector_valid",
    ):
        np.testing.assert_array_equal(getattr(local, stem), expected[stem])

    for stem in (
        "hbond_pairs_index", "hbond_pairs_geometry",
        "hbond_pairs_angle_valid",
    ):
        np.testing.assert_array_equal(
            getattr(protein.hbond, stem.removeprefix("hbond_")),
            expected[stem])

    assert isinstance(protein.dssp, DsspGroup)
    for stem in (
        "dssp_torsion_angle", "dssp_torsion_sin", "dssp_torsion_cos",
        "dssp_torsion_valid", "dssp_hbond_partner_residue_index",
    ):
        np.testing.assert_array_equal(
            getattr(protein.dssp, stem.removeprefix("dssp_")),
            expected[stem])

    for stem in ("omega_sin", "omega_cos", "omega_valid"):
        np.testing.assert_array_equal(
            getattr(protein.planar_geometry, stem),
            expected[stem])


def _write_full9_group(trajectory, name, values, frame_indices, frame_times,
                       units, mcconnell=False):
    group = trajectory.create_group(name)
    attrs = {
        "result_name": f"{name}Result",
        "n_atoms": values.shape[0],
        "n_frames": values.shape[1],
        "finalized": True,
        "tensor_basis": "project_native_full9_spherical_tensor_v1",
        "tensor_component_order": (
            "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"),
        "tensor_frame": "conformation_cartesian_xyz",
        "parity": "0e+1e+2e",
        "tensor_parity": "even",
        "tensor_transformation": "even_rank2: T'=R T R^T",
        "tensor_t1_semantics": "producer tensor T1 semantics",
        "tensor_t1_structural_zero": False,
        "tensor_structural_zero_components": "none",
        "e3nn_export": "explicit conversion required",
        "normalization": "isometric_real_sph",
        "normalization_scope": "producer normalization scope",
        "units": units,
    }
    if mcconnell:
        attrs["coordinate_frame"] = "conformation_cartesian_xyz"
    else:
        attrs["irrep_layout"] = attrs["tensor_component_order"]
    for key, value in attrs.items():
        group.attrs[key] = value
    group.create_dataset("xyz", data=values)
    group.create_dataset("frame_indices", data=frame_indices)
    group.create_dataset("frame_times", data=frame_times)
    return {
        "xyz": values,
        "frame_indices": frame_indices,
        "frame_times": frame_times,
        **attrs,
    }


def _assert_complete_surface(actual, expected):
    assert {field.name for field in fields(actual)} == set(expected)
    for name, value in expected.items():
        observed = getattr(actual, name)
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(observed, value)
        else:
            assert observed == value


def test_all_previously_unread_live_h5_groups_are_reachable(tmp_path):
    n_atoms = 2
    n_frames = 3
    root_indices = np.array([10, 20, 30], dtype=np.uint64)
    root_times = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    result_indices = np.array([101, 102, 103], dtype=np.uint64)
    result_times = np.array([1.1, 2.1, 3.1], dtype=np.float64)
    path = tmp_path / "trajectory.h5"

    expected = {}
    with h5py.File(path, "w") as h5:
        h5.attrs["protein_id"] = "reachability"
        h5.attrs["n_atoms"] = n_atoms
        h5.attrs["finalized"] = True

        atoms = h5.create_group("atoms")
        expected["atoms"] = {
            "element": np.array([1, 6], dtype=np.int32),
            "residue_index": np.array([0, 1], dtype=np.uint64),
            "pdb_atom_name": np.array([b"H", b"CA"]),
        }
        atoms.create_dataset("element", data=expected["atoms"]["element"])
        atoms.create_dataset(
            "residue_index", data=expected["atoms"]["residue_index"])
        atoms.create_dataset(
            "pdb_atom_name", data=np.array(["H", "CA"],
                                             dtype=h5py.string_dtype()))

        trajectory = h5.create_group("trajectory")
        provenance = {
            "xtc_path": "input.xtc",
            "tpr_path": "input.tpr",
            "edr_path": "input.edr",
            "configuration": "complete",
            "window_mode": "local",
            "window_start": 7,
            "window_len": 3,
        }
        for key, value in provenance.items():
            trajectory.attrs[key] = value

        frames_group = trajectory.create_group("frames")
        frames_group.attrs["n_frames"] = n_frames
        frames_group.create_dataset("time_ps", data=root_times)
        frames_group.create_dataset("original_index", data=root_indices)

        positions = trajectory.create_group("positions")
        positions_xyz = np.arange(
            n_atoms * n_frames * 3, dtype=np.float64).reshape(
                n_atoms, n_frames, 3)
        position_attrs = {
            "result_name": "PositionsTimeSeriesTrajectoryResult",
            "n_atoms": n_atoms,
            "n_frames": n_frames,
            "finalized": True,
        }
        for key, value in position_attrs.items():
            positions.attrs[key] = value
        xyz = positions.create_dataset("xyz", data=positions_xyz)
        xyz_attrs = {
            "physical_class": "cartesian_position",
            "coordinate_frame": "conformation_cartesian_xyz",
            "parity": "odd",
            "transformation": "affine_position: p'=R p+t",
            "units": "Angstrom",
        }
        for key, value in xyz_attrs.items():
            xyz.attrs[key] = value
        positions.create_dataset("frame_indices", data=result_indices)
        positions.create_dataset("frame_times", data=result_times)
        expected["positions"] = {
            "xyz": positions_xyz,
            "frame_indices": result_indices,
            "frame_times": result_times,
            **position_attrs,
            **xyz_attrs,
        }

        expected["bs"] = _write_full9_group(
            trajectory, "bs_shielding_time_series",
            np.arange(54, dtype=np.float64).reshape(2, 3, 9),
            result_indices + 10, result_times + 10.0, "ppm_T_per_nA")
        expected["hm"] = _write_full9_group(
            trajectory, "hm_shielding_time_series",
            np.arange(54, 108, dtype=np.float64).reshape(2, 3, 9),
            result_indices + 20, result_times + 20.0, "Angstrom^-1")
        expected["mc"] = _write_full9_group(
            trajectory, "mc_shielding_time_series",
            np.arange(108, 162, dtype=np.float64).reshape(2, 3, 9),
            result_indices + 30, result_times + 30.0, "Angstrom^-3",
            mcconnell=True)

        charge = trajectory.create_group("aimnet2_charge_time_series")
        charge_values = np.arange(6, dtype=np.float64).reshape(2, 3) / 10.0
        expected["charge"] = {
            "charge": charge_values,
            "frame_indices": result_indices + 40,
            "frame_times": result_times + 40.0,
            "result_name": "AIMNet2ChargeTimeSeriesTrajectoryResult",
            "n_atoms": n_atoms,
            "n_frames": n_frames,
            "finalized": True,
            "irrep_layout": "T0",
            "parity": "0e",
            "units": "elementary_charge",
        }
        for key, value in expected["charge"].items():
            if not isinstance(value, np.ndarray) and key != "charge":
                charge.attrs[key] = value
        charge.create_dataset("charge", data=charge_values)
        charge.create_dataset(
            "frame_indices", data=expected["charge"]["frame_indices"])
        charge.create_dataset(
            "frame_times", data=expected["charge"]["frame_times"])

        sasa = trajectory.create_group("sasa_time_series")
        sasa_values = np.arange(6, dtype=np.float64).reshape(2, 3) + 0.5
        expected["sasa"] = {
            "sasa": sasa_values,
            "frame_indices": result_indices + 50,
            "frame_times": result_times + 50.0,
            "result_name": "SasaTimeSeriesTrajectoryResult",
            "n_atoms": n_atoms,
            "n_frames": n_frames,
            "finalized": True,
            "irrep_layout": "raw_scalar_no_exact_o3_irrep",
            "parity": "mixed",
            "units": "Angstrom^2",
            "coordinate_frame": "lab_fixed_fibonacci_sampling_grid",
            "transformation": "producer finite-grid transformation",
            "directional_metadata_scope": "sasa dataset only",
            "sasa_units": "Angstrom^2",
        }
        for key, value in expected["sasa"].items():
            if (not isinstance(value, np.ndarray)
                    and key not in {"sasa_units"}):
                sasa.attrs[key] = value
        sasa_dataset = sasa.create_dataset("sasa", data=sasa_values)
        sasa_dataset.attrs["units"] = expected["sasa"]["sasa_units"]
        sasa.create_dataset(
            "frame_indices", data=expected["sasa"]["frame_indices"])
        sasa.create_dataset(
            "frame_times", data=expected["sasa"]["frame_times"])

        bonds = trajectory.create_group("bond_length_stats")
        expected["bonds"] = {
            "length_mean": np.array([1.0, 1.1]),
            "length_std": np.array([0.01, 0.02]),
            "length_min": np.array([0.9, 1.0]),
            "length_max": np.array([1.2, 1.3]),
            "length_delta_mean": np.array([0.1, 0.2]),
            "length_delta_std": np.array([0.03, 0.04]),
            "atom_a": np.array([0, 0], dtype=np.uint64),
            "atom_b": np.array([1, 1], dtype=np.uint64),
            "order": np.array([1, 2], dtype=np.int8),
            "category": np.array([3, 4], dtype=np.int8),
            "result_name": "BondLengthStatsTrajectoryResult",
            "n_bonds": 2,
            "n_frames": n_frames,
            "finalized": True,
            "units": "Angstrom",
        }
        for key, value in expected["bonds"].items():
            if isinstance(value, np.ndarray):
                bonds.create_dataset(key, data=value)
            else:
                bonds.attrs[key] = value

        autocorrelation = trajectory.create_group("bs_t0_autocorrelation")
        expected["autocorrelation"] = {
            "rho": np.array([[1.0, 0.5], [1.0, 0.25]]),
            "lag_frames": np.array([0, 1], dtype=np.uint64),
            "lag_times_ps": np.array([0.0, 2.0]),
            "result_name": "BsT0AutocorrelationTrajectoryResult",
            "n_atoms": n_atoms,
            "n_lags": 2,
            "n_frames": n_frames,
            "finalized": True,
            "sample_interval_ps": 2.0,
            "units": "dimensionless",
            "estimator": "biased",
            "mean_convention": "full_range",
        }
        for key, value in expected["autocorrelation"].items():
            if isinstance(value, np.ndarray):
                autocorrelation.create_dataset(key, data=value)
            else:
                autocorrelation.attrs[key] = value

    result = load_trajectory(path)

    np.testing.assert_array_equal(result.frame_times, root_times)
    np.testing.assert_array_equal(result.frame_indices, root_indices)
    np.testing.assert_array_equal(result.positions, positions_xyz.transpose(1, 0, 2))
    assert result.finalized is True
    for key, value in provenance.items():
        assert getattr(result, key) == value

    assert isinstance(result.atoms, TrajectoryAtomsGroup)
    _assert_complete_surface(result.atoms, expected["atoms"])
    assert isinstance(result.positions_time_series, PositionsTimeSeriesGroup)
    _assert_complete_surface(result.positions_time_series, expected["positions"])
    assert isinstance(result.bs_shielding_time_series,
                      BsShieldingTimeSeriesGroup)
    _assert_complete_surface(result.bs_shielding_time_series, expected["bs"])
    assert isinstance(result.hm_shielding_time_series,
                      HmShieldingTimeSeriesGroup)
    _assert_complete_surface(result.hm_shielding_time_series, expected["hm"])
    assert isinstance(result.mc_shielding_time_series,
                      McConnellShieldingTimeSeriesGroup)
    _assert_complete_surface(result.mc_shielding_time_series, expected["mc"])
    assert isinstance(result.aimnet2_charge_time_series,
                      AIMNet2ChargeTimeSeriesGroup)
    _assert_complete_surface(
        result.aimnet2_charge_time_series, expected["charge"])
    assert isinstance(result.sasa_time_series, SasaTimeSeriesGroup)
    _assert_complete_surface(result.sasa_time_series, expected["sasa"])
    assert isinstance(result.bond_length_stats, BondLengthStatsGroup)
    _assert_complete_surface(result.bond_length_stats, expected["bonds"])
    assert isinstance(result.bs_t0_autocorrelation, BsT0AutocorrelationGroup)
    _assert_complete_surface(
        result.bs_t0_autocorrelation, expected["autocorrelation"])
