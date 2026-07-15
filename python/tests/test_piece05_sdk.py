"""Synthetic SDK coverage for the piece-05 electrostatic schema."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from nmr_extract import (
    AIMNet2AimProjectionWelfordGroup,
    AIMNet2ChargeResponseGradient,
    AIMNet2ChargeResponseGradientTimeSeriesGroup,
    ApbsEfgTimeSeriesGroup,
    ApbsEfieldTimeSeriesGroup,
    CATALOG,
    EFGTensor,
    ShieldingTensor,
    VectorField,
    load,
    load_trajectory,
)
from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
)


N_ATOMS = 4
N_FRAMES = 3
BASIS_ID = "splitmix64_0xA17E20260708_achlioptas_32x256_element_HCNOS"
M7_SENTINEL_MAX_RANK = 17
M7_SENTINEL_RANK3_POLICY = "legacy_sdk_probe"


def _save(out_dir, stem: str, shape, dtype=np.float64) -> np.ndarray:
    data = np.arange(np.prod(shape), dtype=dtype).reshape(shape)
    np.save(out_dir / f"{stem}.npy", data)
    return data


def test_piece05_catalog_contract():
    aim_required = {
        "aimnet2_efg_sidechain": 5,
        "aimnet2_E": 3,
        "aimnet2_E_backbone": 3,
        "aimnet2_E_sidechain": 3,
        "aimnet2_E_aromatic": 3,
        "aimnet2_energy_mlp": None,
        "aimnet2_energy_shifted_local": None,
        "aimnet2_energy_terms": 6,
        "aimnet2_d3_e_disp_atom": None,
        "aimnet2_d3_cn": None,
        "aimnet2_d3_c6_stats": 3,
        "aimnet2_aim_projection": 32,
    }
    for stem, cols in aim_required.items():
        assert CATALOG[stem].required is True
        assert CATALOG[stem].cols == cols

    assert CATALOG["aimnet2_energy_terms"].native_axis == "protein"
    assert BASIS_ID in CATALOG["aimnet2_aim_projection"].description
    assert "aimnet_extra" not in CATALOG

    apbs = {
        "apbs_phi": None,
        "apbs_E_clamp_mask": None,
        "apbs_E_clamp_scale": None,
        "apbs_nonfinite_sanitizer_mask": None,
        "apbs_E_total_diagnostic": 3,
        "apbs_efg_total_diagnostic": 5,
    }
    for stem, cols in apbs.items():
        assert CATALOG[stem].required is False
        assert CATALOG[stem].cols == cols
    assert "reaction" in CATALOG["apbs_E"].description.lower()
    assert "reaction" in CATALOG["apbs_efg"].description.lower()

    assert "coulomb_E_solvent" not in CATALOG
    assert "coulomb_efg_solvent" not in CATALOG
    aromatic_projection_doc = CATALOG["coulomb_aromatic_E_proj"].description
    assert "parent-to-H projection" in aromatic_projection_doc
    assert "NaN for non-H or parentless atoms" in aromatic_projection_doc

    assert CATALOG["eeq_chi_eff"].group == "eeq"
    assert CATALOG["eeq_hardness"].cols == 2
    eeq_coulomb_cols = {
        "eeq_coulomb_efg": 9,
        "eeq_coulomb_E": 3,
        "eeq_coulomb_E_backbone": 3,
        "eeq_coulomb_E_sidechain": 3,
        "eeq_coulomb_E_aromatic": 3,
        "eeq_coulomb_efg_backbone": 5,
        "eeq_coulomb_efg_sidechain": 5,
        "eeq_coulomb_efg_aromatic": 5,
        "eeq_coulomb_scalars": 4,
        "eeq_coulomb_aromatic_E_proj": None,
        "eeq_coulomb_aromatic_n_src": None,
    }
    for stem, cols in eeq_coulomb_cols.items():
        assert CATALOG[stem].group == "eeq_coulomb"
        assert CATALOG[stem].cols == cols


def test_static_loader_wires_piece05_groups(tmp_path):
    write_required_sdk_npys(tmp_path, N_ATOMS, n_residues=1)
    write_minimal_topology_sidecar(tmp_path, N_ATOMS, n_residues=1)

    # Optional canonical APBS surfaces.
    _save(tmp_path, "apbs_E", (N_ATOMS, 3))
    _save(tmp_path, "apbs_efg", (N_ATOMS, 5))
    _save(tmp_path, "apbs_phi", (N_ATOMS,))
    _save(tmp_path, "apbs_E_clamp_mask", (N_ATOMS,), np.uint8)
    _save(tmp_path, "apbs_E_clamp_scale", (N_ATOMS,))
    _save(tmp_path, "apbs_nonfinite_sanitizer_mask", (N_ATOMS,), np.uint8)
    _save(tmp_path, "apbs_E_total_diagnostic", (N_ATOMS, 3))
    _save(tmp_path, "apbs_efg_total_diagnostic", (N_ATOMS, 5))

    # Optional EEQ solver diagnostics and its Coulomb calculator.
    _save(tmp_path, "eeq_charges", (N_ATOMS,))
    _save(tmp_path, "eeq_cn", (N_ATOMS,))
    _save(tmp_path, "eeq_chi_eff", (N_ATOMS,))
    _save(tmp_path, "eeq_hardness", (N_ATOMS, 2))
    _save(tmp_path, "eeq_coulomb_efg", (N_ATOMS, 9))
    for stem in (
        "eeq_coulomb_E",
        "eeq_coulomb_E_backbone",
        "eeq_coulomb_E_sidechain",
        "eeq_coulomb_E_aromatic",
    ):
        _save(tmp_path, stem, (N_ATOMS, 3))
    for stem in (
        "eeq_coulomb_efg_backbone",
        "eeq_coulomb_efg_sidechain",
        "eeq_coulomb_efg_aromatic",
    ):
        _save(tmp_path, stem, (N_ATOMS, 5))
    _save(tmp_path, "eeq_coulomb_scalars", (N_ATOMS, 4))
    _save(tmp_path, "eeq_coulomb_aromatic_E_proj", (N_ATOMS,))
    _save(tmp_path, "eeq_coulomb_aromatic_n_src", (N_ATOMS,), np.int32)

    protein = load(tmp_path)

    assert not hasattr(protein.coulomb, "E_solvent")
    assert not hasattr(protein.coulomb, "efg_solvent")
    assert isinstance(protein.apbs.E, VectorField)
    assert isinstance(protein.apbs.efg, EFGTensor)
    assert protein.apbs.phi.shape == (N_ATOMS,)
    assert protein.apbs.E_clamp_mask.dtype == np.uint8
    assert protein.apbs.nonfinite_sanitizer_mask.dtype == np.uint8
    assert isinstance(protein.apbs.E_total_diagnostic, VectorField)
    assert isinstance(protein.apbs.efg_total_diagnostic, EFGTensor)

    assert protein.aimnet2.efg_sidechain.data.shape == (N_ATOMS, 5)
    assert protein.aimnet2.E.data.shape == (N_ATOMS, 3)
    assert protein.aimnet2.energy_terms.shape == (1, 6)
    assert protein.aimnet2.aim_projection.shape == (N_ATOMS, 32)
    assert protein.aimnet2.aim_projection.dtype == np.float32

    assert protein.eeq.chi_eff.shape == (N_ATOMS,)
    assert protein.eeq.hardness.shape == (N_ATOMS, 2)
    assert isinstance(protein.eeq.coulomb.efg, ShieldingTensor)
    assert isinstance(protein.eeq.coulomb.E, VectorField)
    assert isinstance(protein.eeq.coulomb.efg_sidechain, EFGTensor)


def _write_trajectory_root(f: h5py.File) -> None:
    f.attrs["protein_id"] = "piece05"
    f.attrs["n_atoms"] = N_ATOMS
    f.attrs["finalized"] = True
    traj = f.create_group("trajectory")
    frames = traj.create_group("frames")
    frames.attrs["n_frames"] = N_FRAMES
    frames.create_dataset(
        "time_ps", data=np.arange(N_FRAMES, dtype=np.float64) * 0.5)
    frames.create_dataset(
        "original_index", data=np.arange(N_FRAMES, dtype=np.uint64))
    pos = traj.create_group("positions")
    pos.attrs["n_atoms"] = N_ATOMS
    pos.attrs["n_frames"] = N_FRAMES
    pos.attrs["result_name"] = "PositionsTimeSeriesTrajectoryResult"
    pos.attrs["finalized"] = True
    pos.create_dataset(
        "xyz", data=np.zeros((N_ATOMS, N_FRAMES, 3), dtype=np.float64))
    pos.create_dataset(
        "frame_indices", data=np.arange(N_FRAMES, dtype=np.uint64))
    pos.create_dataset(
        "frame_times", data=np.arange(N_FRAMES, dtype=np.float64) * 0.5)


def _write_apbs_common_attrs(g) -> None:
    g.attrs["field_quantity"] = \
        "reaction_field_total_minus_homogeneous_vacuum_reference"
    g.attrs["reference_dielectric"] = 1.0
    g.attrs["reference_ionic_strength_M"] = 0.0
    g.attrs["reference_mobile_ion_count"] = 0
    g.attrs["reference_subtracts"] = "all_solute_charges"
    g.attrs["diagnostic_total_unclamped"] = True
    g.attrs["apbs_grid_mode"] = "single_manual"
    g.attrs["apbs_manual_grid_padding_A"] = 70.0
    g.attrs["apbs_manual_grid_min_dim_A"] = 70.0
    g.attrs["apbs_temperature_K"] = 310.0
    g.attrs["apbs_thermal_voltage_V"] = 0.02672
    g.attrs["max_potential_derivative_rank"] = 2
    g.attrs["higher_derivatives_present"] = False
    g.attrs["rank3_policy"] = "not_emitted_no_local_frame"
    g.attrs["source"] = "ApbsFieldResult"
    g.attrs["source_attached_policy"] = \
        "always_attached_with_defensive_absence_mask"


def _write_apbs_grid_datasets(g) -> None:
    dims = np.tile(np.array([97, 97, 97], dtype=np.uint64), (N_FRAMES, 1))
    lengths = np.tile(np.array([70.0, 71.0, 72.0]), (N_FRAMES, 1))
    origin = np.tile(np.array([-35.0, -35.5, -36.0]), (N_FRAMES, 1))
    spacing = lengths / 96.0
    g.create_dataset("apbs_grid_dims_per_frame", data=dims)
    g.create_dataset("apbs_grid_lengths_A_per_frame", data=lengths)
    g.create_dataset("apbs_grid_origin_A_per_frame", data=origin)
    g.create_dataset("apbs_grid_spacing_A_per_frame", data=spacing)


def _write_apbs_groups(f: h5py.File, *, legacy: bool = False) -> None:
    fi = np.arange(N_FRAMES, dtype=np.uint64)
    ft = np.arange(N_FRAMES, dtype=np.float64) * 0.5

    e = f.create_group("/trajectory/apbs_efield_time_series")
    e.attrs["n_atoms"] = N_ATOMS
    e.attrs["n_frames"] = N_FRAMES
    e.attrs["irrep_layout"] = "x,y,z"
    e.attrs["normalization"] = "cartesian"
    e.attrs["parity"] = "1o"
    e.attrs["units"] = "V/Angstrom"
    e.create_dataset(
        "xyz", data=np.arange(N_ATOMS * N_FRAMES * 3, dtype=np.float64)
        .reshape(N_ATOMS, N_FRAMES, 3))
    e.create_dataset("frame_indices", data=fi)
    e.create_dataset("frame_times", data=ft)
    if not legacy:
        _write_apbs_common_attrs(e)
        e.attrs["efield_clamp_units"] = "V/Angstrom"
        e.attrs["efield_clamp_threshold"] = 100.0
        e.create_dataset(
            "source_attached_per_frame",
            data=np.array([1, 0, 1], dtype=np.uint8))
        e.create_dataset(
            "clamp_mask",
            data=np.zeros((N_ATOMS, N_FRAMES), dtype=np.uint8))
        clamp_scale = np.ones((N_ATOMS, N_FRAMES), dtype=np.float64)
        clamp_scale[:, 1] = np.nan
        e.create_dataset("clamp_scale", data=clamp_scale)
        _write_apbs_grid_datasets(e)

    efg = f.create_group("/trajectory/apbs_efg_time_series")
    efg.attrs["n_atoms"] = N_ATOMS
    efg.attrs["n_frames"] = N_FRAMES
    efg.attrs["t2_basis"] = "project_native_t2_isometric_real_tesseral_v1"
    efg.attrs["t2_component_order"] = \
        "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
    efg.attrs["t2_frame"] = "cartesian_xyz_emitted_frame"
    efg.attrs["t2_parity"] = "even"
    efg.attrs["units"] = "V/Angstrom^2"
    efg.create_dataset(
        "t2", data=np.arange(N_ATOMS * N_FRAMES * 5, dtype=np.float64)
        .reshape(N_ATOMS, N_FRAMES, 5))
    efg.create_dataset("frame_indices", data=fi)
    efg.create_dataset("frame_times", data=ft)
    if not legacy:
        _write_apbs_common_attrs(efg)
        efg.create_dataset(
            "source_attached_per_frame",
            data=np.array([1, 0, 1], dtype=np.uint8))
        _write_apbs_grid_datasets(efg)


def _write_projection_group(f: h5py.File) -> None:
    g = f.create_group("/trajectory/aimnet2_aim_projection_welford")
    g.attrs["n_atoms"] = N_ATOMS
    g.attrs["n_frames"] = N_FRAMES
    g.attrs["source_attached_count"] = N_FRAMES
    g.attrs["projection_dim"] = 32
    g.attrs["basis_id"] = BASIS_ID
    g.attrs["units"] = "dimensionless"
    g.attrs["irrep_layout"] = "feature_vector"
    g.attrs["parity"] = "0e"
    g.attrs["source"] = \
        "AIMNet2Result.aimnet2_aim projected by committed element basis"
    g.attrs["source_attached_policy"] = \
        "always_attached_with_skip_on_absent_source"
    g.attrs["finalized"] = True
    shape = (N_ATOMS, 32)
    for name in (
        "projection_mean", "projection_m2", "projection_std",
        "projection_min", "projection_max",
    ):
        g.create_dataset(name, data=np.zeros(shape, dtype=np.float64))
    for name in ("projection_min_frame", "projection_max_frame"):
        g.create_dataset(name, data=np.zeros(shape, dtype=np.uint64))
    g.create_dataset(
        "n_frames_per_atom",
        data=np.full(N_ATOMS, N_FRAMES, dtype=np.uint64))
    g.create_dataset(
        "frame_indices", data=np.arange(N_FRAMES, dtype=np.uint64))
    g.create_dataset(
        "frame_times", data=np.arange(N_FRAMES, dtype=np.float64) * 0.5)
    g.create_dataset(
        "source_attached_per_frame",
        data=np.ones(N_FRAMES, dtype=np.uint8))


def _write_mopac_coulomb_group(f: h5py.File) -> None:
    def write_group(path: str) -> h5py.Group:
        g = f.create_group(path)
        g.attrs["n_atoms"] = N_ATOMS
        g.attrs["n_frames"] = N_FRAMES
        g.attrs["source_attached_count"] = N_FRAMES
        g.create_dataset(
            "t2",
            data=np.zeros((N_ATOMS, N_FRAMES, 5), dtype=np.float64))
        g.create_dataset(
            "frame_indices", data=np.arange(N_FRAMES, dtype=np.uint64))
        g.create_dataset(
            "frame_times", data=np.arange(N_FRAMES, dtype=np.float64) * 0.5)
        g.create_dataset(
            "source_attached_per_frame",
            data=np.ones(N_FRAMES, dtype=np.uint8))
        return g

    # Mirror the C++ writer's two-group shape: canonical owns the payload;
    # M7 derivative-policy attributes live only on the historical alias.
    write_group("/trajectory/mopac_coulomb_efg_time_series")
    legacy = write_group(
        "/trajectory/mopac_coulomb_shielding_time_series")
    # Deliberately non-default sentinels make it impossible for SDK defaults
    # to masquerade as a successful read from the legacy group.
    legacy.attrs["max_potential_derivative_rank"] = M7_SENTINEL_MAX_RANK
    legacy.attrs["higher_derivatives_present"] = True
    legacy.attrs["rank3_policy"] = M7_SENTINEL_RANK3_POLICY


def _write_mopac_mc_group(f: h5py.File, *, legacy: bool) -> None:
    g = f.create_group("/trajectory/mopac_mc_shielding_time_series")
    g.attrs["n_atoms"] = N_ATOMS
    g.attrs["n_frames"] = N_FRAMES
    g.attrs["source_attached_count"] = N_FRAMES
    g.attrs["normalization"] = "isometric_real_sph"
    g.attrs["units"] = "Angstrom^-3"
    g.attrs["source"] = \
        "MopacMcConnellResult.mopac_mc_shielding_contribution"
    g.attrs["source_attached_policy"] = "conditional"
    if legacy:
        g.attrs["irrep_layout"] = \
            "0e,1e_x,1e_y,1e_z,2e_m-2..+2"
        g.attrs["parity"] = "0e+1e+2e"
    else:
        g.attrs["tensor_basis"] = \
            "project_native_full9_spherical_tensor_v1"
        g.attrs["tensor_component_order"] = \
            "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        g.attrs["tensor_frame"] = "conformation_cartesian_xyz"
        g.attrs["tensor_parity"] = "even"
        g.attrs["e3nn_export"] = \
            "explicit project-basis to e3nn conversion required before use"
    g.create_dataset(
        "xyz",
        data=np.zeros((N_ATOMS, N_FRAMES, 9), dtype=np.float64))
    g.create_dataset(
        "frame_indices", data=np.arange(N_FRAMES, dtype=np.uint64))
    g.create_dataset(
        "frame_times", data=np.arange(N_FRAMES, dtype=np.float64) * 0.5)
    g.create_dataset(
        "source_attached_per_frame",
        data=np.ones(N_FRAMES, dtype=np.uint8))


def test_trajectory_piece05_round_trip(tmp_path):
    h5 = tmp_path / "piece05.h5"
    with h5py.File(h5, "w") as f:
        _write_trajectory_root(f)
        _write_apbs_groups(f)
        _write_projection_group(f)
        _write_mopac_coulomb_group(f)

    traj = load_trajectory(h5)
    e = traj.apbs_efield
    assert isinstance(e, ApbsEfieldTimeSeriesGroup)
    assert e.xyz.shape == (N_ATOMS, N_FRAMES, 3)
    assert e.clamp_mask.dtype == np.uint8
    assert np.isnan(e.clamp_scale[:, 1]).all()
    assert e.field_quantity == \
        "reaction_field_total_minus_homogeneous_vacuum_reference"
    assert e.reference_subtracts == "all_solute_charges"
    assert e.apbs_grid_dims_per_frame.shape == (N_FRAMES, 3)
    assert e.apbs_temperature_K == 310.0
    assert e.max_potential_derivative_rank == 2
    assert e.higher_derivatives_present is False
    assert e.rank3_policy == "not_emitted_no_local_frame"

    efg = traj.apbs_efg
    assert isinstance(efg, ApbsEfgTimeSeriesGroup)
    assert efg.t2.shape == (N_ATOMS, N_FRAMES, 5)
    assert efg.apbs_grid_spacing_A_per_frame.shape == (N_FRAMES, 3)
    assert efg.reference_dielectric == 1.0
    assert efg.max_potential_derivative_rank == 2

    projection = traj.aimnet2_aim_projection_welford
    assert isinstance(projection, AIMNet2AimProjectionWelfordGroup)
    assert projection.projection_mean.shape == (N_ATOMS, 32)
    assert projection.projection_dim == 32
    assert projection.basis_id == BASIS_ID
    np.testing.assert_array_equal(
        projection.n_frames_per_atom,
        np.full(N_ATOMS, N_FRAMES, dtype=np.uint64))

    mopac = traj.mopac_coulomb_efg_time_series
    assert mopac.max_potential_derivative_rank == M7_SENTINEL_MAX_RANK
    assert mopac.higher_derivatives_present is True
    assert mopac.rank3_policy == M7_SENTINEL_RANK3_POLICY
    assert traj.mopac_coulomb_shielding_time_series is mopac


def test_apbs_old_h5_compatibility_synthesizes_documented_defaults(tmp_path):
    h5 = tmp_path / "legacy_apbs.h5"
    with h5py.File(h5, "w") as f:
        _write_trajectory_root(f)
        _write_apbs_groups(f, legacy=True)

    traj = load_trajectory(h5)
    e = traj.apbs_efield
    np.testing.assert_array_equal(
        e.source_attached_per_frame,
        np.ones(N_FRAMES, dtype=np.uint8))
    np.testing.assert_array_equal(
        e.clamp_mask,
        np.zeros((N_ATOMS, N_FRAMES), dtype=np.uint8))
    np.testing.assert_array_equal(
        e.clamp_scale,
        np.ones((N_ATOMS, N_FRAMES), dtype=np.float64))
    assert e.apbs_grid_dims_per_frame.shape == (0, 3)
    assert e.apbs_grid_dims_per_frame.dtype == np.uint64
    assert e.apbs_grid_lengths_A_per_frame.shape == (0, 3)
    assert np.isnan(e.apbs_temperature_K)
    assert e.field_quantity == ""
    assert e.diagnostic_total_unclamped is False
    assert e.max_potential_derivative_rank == 2
    assert e.rank3_policy == "not_emitted_no_local_frame"

    efg = traj.apbs_efg
    np.testing.assert_array_equal(
        efg.source_attached_per_frame,
        np.ones(N_FRAMES, dtype=np.uint8))
    assert efg.apbs_grid_origin_A_per_frame.shape == (0, 3)
    assert np.isnan(efg.apbs_thermal_voltage_V)
    assert efg.higher_derivatives_present is False


def test_mopac_mc_project_metadata_round_trip(tmp_path):
    h5 = tmp_path / "mopac_mc_project_metadata.h5"
    with h5py.File(h5, "w") as f:
        _write_trajectory_root(f)
        _write_mopac_mc_group(f, legacy=False)

    group = load_trajectory(h5).mopac_mc_shielding_time_series
    assert group.xyz.shape == (N_ATOMS, N_FRAMES, 9)
    assert group.tensor_basis == \
        "project_native_full9_spherical_tensor_v1"
    assert group.tensor_component_order == \
        "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
    assert group.tensor_frame == "conformation_cartesian_xyz"
    assert group.tensor_parity == "even"
    assert group.e3nn_export == \
        "explicit project-basis to e3nn conversion required before use"
    assert group.legacy_irrep_layout == ""
    assert group.legacy_parity == ""
    with pytest.warns(DeprecationWarning, match="legacy_irrep_layout"):
        assert group.irrep_layout == ""
    with pytest.warns(DeprecationWarning, match="legacy_parity"):
        assert group.parity == ""


def test_mopac_mc_legacy_metadata_gets_honest_authoritative_defaults(tmp_path):
    h5 = tmp_path / "mopac_mc_legacy_metadata.h5"
    with h5py.File(h5, "w") as f:
        _write_trajectory_root(f)
        _write_mopac_mc_group(f, legacy=True)

    group = load_trajectory(h5).mopac_mc_shielding_time_series
    assert group.tensor_basis == \
        "project_native_full9_spherical_tensor_v1"
    assert group.tensor_component_order == \
        "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
    assert group.tensor_frame == "conformation_cartesian_xyz"
    assert group.tensor_parity == "even"
    assert "project_full9_to_e3nn(xyz)" in group.e3nn_export
    assert "project_t2_to_e3nn" not in group.e3nn_export
    assert "to_e3nn()" not in group.e3nn_export
    assert group.legacy_irrep_layout == \
        "0e,1e_x,1e_y,1e_z,2e_m-2..+2"
    assert group.legacy_parity == "0e+1e+2e"
    with pytest.warns(DeprecationWarning, match="legacy_irrep_layout"):
        assert group.irrep_layout == group.legacy_irrep_layout
    with pytest.warns(DeprecationWarning, match="legacy_parity"):
        assert group.parity == group.legacy_parity


def test_charge_response_labels_are_unambiguous():
    for stem in (
        "aimnet2_charge_response_gradient",
        "aimnet2_charge_response_gradient_scalar",
    ):
        description = CATALOG[stem].description
        assert "proxy/diagnostic" in description
        assert "NOT Buckingham" in description

    tensor_doc = AIMNet2ChargeResponseGradient.__doc__
    trajectory_doc = AIMNet2ChargeResponseGradientTimeSeriesGroup.__doc__
    assert "proxy/diagnostic" in tensor_doc
    assert "NOT a Buckingham" in tensor_doc
    assert "proxy/diagnostic" in trajectory_doc
    assert "NOT Buckingham" in trajectory_doc
