"""SDK round-trip tests for the topology sidecar.

Verifies that ``nmr_extract.load()`` exposes the new ``Protein.topology``
group (TopologyGroup) with ``bonds``, ``rings``, ``ring_membership``,
and parsed ``manifest`` attributes, and that the catalog declares
``native_axis`` / ``mechanism`` / ``irreps`` / ``units`` /
``sign_convention`` / ``tensor_rank`` / ``parity`` for every entry.
"""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
import pytest

from nmr_extract import (
    CATALOG,
    PositionField,
    ShieldingTensor,
    Residues,
    Bonds,
    Rings,
    RingMembership,
    ExtractionManifest,
    TopologyGroup,
    load,
)
from nmr_extract._catalog import ALLOWED_NATIVE_AXES
from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
    _RESIDUES_DTYPE,
    _BONDS_DTYPE,
    _RINGS_DTYPE,
    _RING_MEMBERSHIP_DTYPE,
)


N_ATOMS = 16


_DIRECTIONAL_STEMS = {
    # Identity, sparse geometry, local frames, and positions.
    "pos", "ring_contributions", "ring_direction_to_center",
    "ring_geometry", "ring_pair_geometry", "bond_direction",
    "bond_geometry_valid", "cb_deviation", "cb_deviation_valid",
    "cb_residual_vector", "cb_residual_vector_valid", "spatial_neighbors",
    "hbond_nearest_dir", "hbond_flags", "hbond_pairs_geometry",
    "hbond_pairs_index", "hbond_pairs_angle_valid",
    "atom_sasa", "atom_sasa_fraction",
    "sasa_normal", "water_polarization",
    "water_hbond_candidates", "water_hbond_counts", "water_hbond_nearest",
    "sidechain_co_frame", "sidechain_co_source_bonds",
    "sidechain_co_frame_quality",
    "piquad_local_tensor", "piquad_local_T2", "piquad_local_frame",
    "piquad_local_geometry", "piquad_quad_scalar",
    # Ring-current vectors/tensors.
    "bs_shielding", "bs_per_type_T1", "bs_per_type_T2", "bs_total_B",
    "bs_ring_B_field", "bs_ring_B_cylindrical", "hm_shielding",
    "hm_per_type_T1", "hm_per_type_T2", "hm_ring_B_field",
    # McConnell and typed side-chain carbonyl outputs.
    "mc_peptide_co_fixed", "mc_peptide_co_bo", "mc_peptide_co_rhombic",
    "mc_peptide_cn_fixed", "mc_peptide_cn_bo",
    "mc_backbone_other_fixed", "mc_backbone_other_bo",
    "mc_sidechain_co_fixed", "mc_sidechain_co_bo",
    "mc_sidechain_other_fixed", "mc_sidechain_other_bo",
    "mc_disulfide_fixed", "mc_disulfide_bo",
    "mc_aromatic_fixed", "mc_aromatic_bo",
    "mc_backbone_xh_fixed", "mc_backbone_xh_bo",
    "mc_sidechain_xh_fixed", "mc_sidechain_xh_bo",
    "mc_s_h_fixed", "mc_s_h_bo", "mc_nearest_co_dir",
    "mc_nearest_co_midpoint", "mc_nearest_co_T2", "mc_nearest_cn_T2",
    "mc_nearest_co_bond_index", "mc_nearest_cn_bond_index",
    "mopac_mc_co_sum", "mopac_mc_cn_sum", "mopac_mc_sidechain_sum",
    "mopac_mc_aromatic_sum", "mopac_mc_co_nearest",
    "mopac_mc_nearest_co_dist", "mopac_mc_nearest_cn_dist",
    "mopac_mc_nearest_co_T2", "mopac_mc_nearest_cn_T2",
    "mopac_mc_backbone_total", "mopac_mc_sidechain_total",
    "mopac_mc_aromatic_total", "mopac_mc_shielding",
    "sidechain_co_fixed_T2", "sidechain_co_bo_T2",
    # Coulomb/APBS/water/AIMNet2 vectors and EFGs.
    "coulomb_efg", "coulomb_efg_t2", "coulomb_E",
    "coulomb_E_backbone", "coulomb_E_sidechain", "coulomb_E_aromatic",
    "coulomb_efg_backbone", "coulomb_efg_sidechain",
    "coulomb_efg_aromatic",
    "mopac_coulomb_efg", "mopac_coulomb_E",
    "mopac_coulomb_E_backbone", "mopac_coulomb_E_sidechain",
    "mopac_coulomb_E_aromatic", "mopac_coulomb_efg_backbone",
    "mopac_coulomb_efg_sidechain", "mopac_coulomb_efg_aromatic",
    "mopac_coulomb_aromatic_E_proj", "mopac_coulomb_aromatic_n_src",
    "mopac_coulomb_E_clamp_mask", "mopac_coulomb_E_clamp_scale",
    "eeq_coulomb_efg", "eeq_coulomb_E", "eeq_coulomb_E_backbone",
    "eeq_coulomb_E_sidechain", "eeq_coulomb_E_aromatic",
    "eeq_coulomb_efg_backbone", "eeq_coulomb_efg_sidechain",
    "eeq_coulomb_efg_aromatic", "apbs_E", "apbs_efg", "apbs_phi",
    "apbs_E_clamp_mask", "apbs_E_clamp_scale",
    "apbs_nonfinite_sanitizer_mask",
    "apbs_E_total_diagnostic", "apbs_efg_total_diagnostic",
    "water_efield", "water_efield_first", "water_efg",
    "water_efg_first", "water_shell_counts", "water_efield_clamp_mask",
    "water_efield_clamp_scale", "water_efield_first_clamp_mask",
    "water_efield_first_clamp_scale", "aimnet2_E", "aimnet2_E_backbone",
    "aimnet2_E_sidechain", "aimnet2_E_aromatic", "aimnet2_efg",
    "aimnet2_efg_aromatic", "aimnet2_efg_backbone",
    "aimnet2_efg_sidechain", "aimnet2_charge_response_gradient",
    "aimnet2_charge_response_gradient_scalar",
    # External quantum/reference and mutation tensors.
    "orca_total", "orca_diamagnetic", "orca_paramagnetic",
    "larsen_hbond_shielding", "larsen_hbond_1pHB_shielding",
    "larsen_hbond_2pHB_shielding", "larsen_hbond_1pHaB_shielding",
    "larsen_hbond_2pHaB_shielding",
    "larsen_hbond_diagnostic_CB_shielding", "delta_shielding",
    "wt_shielding_diamagnetic", "wt_shielding_paramagnetic",
    "mut_shielding_diamagnetic", "mut_shielding_paramagnetic",
    "delta_shielding_diamagnetic", "delta_shielding_paramagnetic",
    "delta_apbs", "delta_ring_proximity", "delta_scalars", "delta_graph",
    "mutation_atom_map",
    # External-frame and reflection-sensitive diagnostics.
    "mopac_global", "mopac_atom_populations",
    "mopac_atomic_orbital_populations",
    "mopac_atomic_orbital_population_totals", "gromacs_energy",
    "gromacs_box",
    "mopac_dipole_debye", "mopac_dipole_point_charge_debye",
    "mopac_dipole_hybridization_debye", "mopac_atom_ao_density",
    "mopac_atomic_orbital_populations_full_precision",
    "mopac_bond_ao_density_directed", "mopac_bond_ao_density",
    "mopac_lmo_occupied_coefficients", "mopac_lmo_virtual_coefficients",
    "bonded_energy",
    "dssp_observed", "dssp_backbone", "dssp_ss8", "dssp_ppii",
    "dssp_hbond_energy", "dssp_hbond_partner_residue_index", "dssp_chi",
    "dssp_torsion_angle", "dssp_torsion_sin", "dssp_torsion_cos",
    "dssp_torsion_valid", "omega_actual", "omega_valid",
    "omega_deviation", "omega_sin", "omega_cos", "aromatic_chi2",
    "pucker_Q", "pucker_theta", "larsen_hbond_pairs_geometry",
    "larsen_hbond_pairs_isotropic", "larsen_hbond_pairs_index",
    "larsen_hbond_count", "larsen_hbond_water_term",
    "larsen_corner_imputed",
    "larsen_imputed_pair_count", "larsen_sidechain_carbonyl_pair_count",
    "larsen_hbond_pairs",
    "larsen_sidechain_donor_candidates", "larsen_sidechain_donor_atoms",
}


# ── Catalog metadata invariants ────────────────────────────────────


class TestCatalogMetadata:
    """Resolves OI-016: every catalog entry declares axis + irreps + units."""

    def test_every_entry_declares_native_axis(self):
        bad = [k for k, s in CATALOG.items() if not s.native_axis]
        assert not bad, f"Entries missing native_axis: {bad}"

    def test_every_entry_declares_mechanism(self):
        bad = [k for k, s in CATALOG.items() if not s.mechanism]
        assert not bad, f"Entries missing mechanism: {bad}"

    def test_every_entry_declares_units_and_truth_contract(self):
        for field in ("units", "coordinate_frame", "transformation", "validity"):
            bad = [k for k, s in CATALOG.items() if not getattr(s, field)]
            assert not bad, f"Entries missing {field}: {bad}"
        vague = [k for k, s in CATALOG.items() if s.units == "mixed"]
        assert not vague, f"Entries with unqualified mixed units: {vague}"

    def test_geometry_kernels_declare_deferred_physical_scaling(self):
        stems = {
            "ring_contributions", "ringchi_scalar", "ringchi_per_type_T0",
            "disp_per_type_T0", "aromatic_r6_proximity_per_type_T0",
            "bs_total_B", "bs_ring_B_field", "bs_ring_B_cylindrical",
            "hm_shielding", "hm_per_type_T0", "hm_per_type_T1",
            "hm_per_type_T2", "hm_ring_B_field", "hm_welford",
            "mc_welford", "sidechain_co_fixed_T2", "sidechain_co_bo_T2",
            "sidechain_co_scalar_audit",
        }
        stems.update(
            stem for stem, spec in CATALOG.items()
            if spec.group == "mcconnell" and spec.tensor_rank == 2
        )
        stems.update(
            stem for stem, spec in CATALOG.items()
            if spec.group == "pi_quadrupole" and
            ("kernel" in spec.description.lower() or
             "tensor" in spec.description.lower())
        )
        missing = [stem for stem in sorted(stems)
                   if not CATALOG[stem].scaling_contract]
        assert not missing, f"Geometry kernels missing scaling contract: {missing}"

    def test_every_entry_has_valid_parity(self):
        bad = [
            k for k, s in CATALOG.items()
            if s.parity not in ("even", "odd", "mixed")
        ]
        assert not bad, f"Entries with invalid parity: {bad}"
        assert {
            k for k, s in CATALOG.items() if s.parity == "mixed"
        } == {
            "atom_sasa", "atom_sasa_fraction", "apbs_E_clamp_mask",
            "apbs_E_clamp_scale", "apbs_nonfinite_sanitizer_mask", "apbs_phi",
            "bonded_energy", "bs_ring_B_cylindrical", "cb_deviation",
            "cb_residual_vector", "delta_apbs", "delta_ring_proximity",
            "dssp_backbone", "dssp_chi", "dssp_ppii",
            "gromacs_energy", "hbond_pairs_geometry", "larsen_corner_imputed",
            "larsen_hbond_pairs", "larsen_hbond_pairs_geometry",
            "larsen_hbond_pairs_isotropic",
            "larsen_imputed_pair_count", "larsen_sidechain_donor_candidates", "mopac_atom_ao_density",
            "mopac_bond_ao_density", "mopac_bond_ao_density_directed", "mopac_global",
            "mopac_lmo_occupied_coefficients", "mopac_lmo_occupied_coefficient_storage_native", "mopac_lmo_virtual_coefficients",
            "mopac_lmo_virtual_coefficient_storage_native", "piquad_local_T2", "piquad_local_frame",
            "piquad_local_geometry", "piquad_local_tensor", "pucker_theta",
            "ring_contributions", "ring_geometry", "ring_pair_geometry",
            "sasa_normal", "sasa_welford", "sidechain_co_frame",
            "spatial_neighbors", "water_hbond_candidates", "water_polarization",
        }

    def test_tensor_rank_is_in_0_1_2(self):
        bad = [k for k, s in CATALOG.items() if s.tensor_rank not in (0, 1, 2)]
        assert not bad, f"Entries with invalid tensor_rank: {bad}"

    def test_topology_sidecar_entries_registered(self):
        for stem in ("residues", "bonds", "rings", "ring_membership"):
            assert stem in CATALOG, f"CATALOG missing topology entry {stem!r}"
            assert CATALOG[stem].required, (
                f"topology sidecar entry {stem!r} must be required")
            assert CATALOG[stem].mechanism == "topology"

    def test_axis_value_set(self):
        """Native axes are drawn from a closed set declared by the contract."""
        seen = {s.native_axis for s in CATALOG.values()}
        unknown = seen - ALLOWED_NATIVE_AXES
        assert not unknown, f"Unknown native_axis values: {unknown}"

    def test_shielding_entries_carry_sign_convention(self):
        """Every shielding tensor entry declares its sign convention."""
        for stem, spec in CATALOG.items():
            if (
                "shielding" in stem
                and spec.tensor_rank == 2
                and spec.mechanism != "electrostatic_efg"
            ):
                assert spec.sign_convention, (
                    f"shielding entry {stem!r} missing sign_convention")

    def test_every_inventoried_directional_name_has_explicit_contract(self):
        missing_names = _DIRECTIONAL_STEMS - CATALOG.keys()
        assert not missing_names, f"Directional NPYs absent from catalog: {missing_names}"

        missing_fields = {
            stem: tuple(
                field for field in
                ("coordinate_frame", "transformation", "validity")
                if not getattr(CATALOG[stem], field)
            )
            for stem in sorted(_DIRECTIONAL_STEMS)
        }
        missing_fields = {
            stem: fields for stem, fields in missing_fields.items() if fields
        }
        assert not missing_fields, (
            f"Directional NPYs missing freeze-contract fields: {missing_fields}")

    @pytest.mark.parametrize("stem", [
        "bond_direction", "ring_direction_to_center", "mc_nearest_co_dir",
        "hbond_nearest_dir",
        "eeq_coulomb_E", "eeq_coulomb_E_backbone",
        "eeq_coulomb_E_sidechain", "eeq_coulomb_E_aromatic",
        "water_efield", "water_efield_first",
    ])
    def test_exact_global_polar_vector_law(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == "polar_vector: v'=R v"

    @pytest.mark.parametrize("stem", [
        "aimnet2_E", "aimnet2_E_backbone", "aimnet2_E_sidechain",
        "aimnet2_E_aromatic", "aimnet2_charge_response_gradient",
    ])
    def test_aimnet_vectors_are_polar(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.irreps == "1o"
        assert spec.parity == "odd"
        assert spec.tensor_rank == 1

    @pytest.mark.parametrize("stem", [
        "coulomb_E", "coulomb_E_backbone", "coulomb_E_sidechain",
        "coulomb_E_aromatic", "mopac_coulomb_E",
        "mopac_coulomb_E_backbone", "mopac_coulomb_E_sidechain",
        "mopac_coulomb_E_aromatic",
    ])
    def test_coulomb_family_E_is_polar(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == "polar_vector: v'=R v"
        assert spec.irreps == "1o"
        assert spec.parity == "odd"

    @pytest.mark.parametrize("stem", [
        "mopac_coulomb_aromatic_E_proj",
        "mopac_coulomb_aromatic_n_src",
        "mopac_coulomb_E_clamp_mask",
        "mopac_coulomb_E_clamp_scale",
        "mopac_mc_co_sum", "mopac_mc_cn_sum",
        "mopac_mc_sidechain_sum", "mopac_mc_aromatic_sum",
        "mopac_mc_co_nearest", "mopac_mc_nearest_co_dist",
        "mopac_mc_nearest_cn_dist",
    ])
    def test_new_mopac_scalar_surfaces_are_intrinsic(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "intrinsic_geometry"
        assert spec.irreps == "0e"
        assert spec.parity == "even"
        assert spec.tensor_rank == 0
        assert "O(3)-invariant" in spec.transformation

    @pytest.mark.parametrize("stem", [
        "mc_nearest_co_bond_index", "mc_nearest_cn_bond_index",
    ])
    def test_mcconnell_nearest_bond_indices_are_geometry_owned(self, stem):
        spec = CATALOG[stem]
        assert spec.group == "mcconnell"
        assert spec.coordinate_frame == "intrinsic_geometry"
        assert spec.units == "index"
        assert spec.irreps == "0e"
        assert spec.parity == "even"
        assert spec.tensor_rank == 0
        assert "solely by accepted source-midpoint geometry" in (
            spec.transformation)
        assert "-1 means no accepted source" in spec.validity
        assert "mopac_topology_bond_orders_full.npy" in spec.validity


    @pytest.mark.parametrize("stem", [
        "apbs_E_clamp_mask", "apbs_E_clamp_scale",
    ])
    def test_apbs_clamp_diagnostics_declare_finite_grid_scalar_law(
            self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == (
            "lab_fixed_apbs_finite_difference_grid")
        assert spec.irreps == ""
        assert spec.parity == "mixed"
        assert spec.tensor_rank == 0
        assert "continuum rotation/translation/reflection-invariant scalar" in (
            spec.transformation)
        assert "no exact O(3) law" in spec.transformation
        assert "ApbsFieldResult is absent if APBS fails" in spec.validity

    def test_apbs_sanitizer_declares_finite_grid_outcome_law(self):
        spec = CATALOG["apbs_nonfinite_sanitizer_mask"]
        assert spec.coordinate_frame == (
            "lab_fixed_apbs_finite_difference_grid")
        assert spec.irreps == ""
        assert spec.parity == "mixed"
        assert spec.tensor_rank == 0
        assert (
            "continuum rotation/translation/reflection-invariant "
            "finite-value diagnostic"
        ) in spec.transformation
        assert "no exact O(3) outcome law" in spec.transformation
        assert "bit0 reaction E" in spec.validity
        assert "bit3 total EFG" in spec.validity
        assert "zero means no sanitizer fired" in spec.validity

    def test_hydration_contract_separates_sampled_normal_descendants(self):
        spec = CATALOG["water_polarization"]
        assert "cols0:3 net-water-dipole polar vector" in spec.transformation
        assert "cols3:6 outward SASA normal is continuum polar" in (
            spec.transformation)
        assert "cols6:8 are continuum invariant" in spec.transformation
        assert "live finite-grid normal" in spec.transformation
        assert "cols8:10 are exact O(3)-invariant scalars" in (
            spec.transformation)

    @pytest.mark.parametrize("stem,fragment", [
        ("bond_geometry_valid", "row-aligned bond direction"),
        ("cb_deviation_valid", "chiral-L-CB output"),
        ("cb_residual_vector_valid", "chiral-L-CB output"),
        ("hbond_flags", "atom classifications"),
        ("hbond_pairs_index", "sparse donor/acceptor"),
        ("hbond_pairs_angle_valid", "angle-availability mask"),
        ("omega_valid", "signed peptide omega"),
        ("sidechain_co_source_bonds", "sparse typed bond/atom/residue"),
        ("sidechain_co_frame_quality", "raw normal norm"),
        ("water_shell_counts", "first/second-shell water counts"),
        ("water_efield_clamp_mask", "threshold diagnostic"),
        ("water_efield_clamp_scale", "scalar scale"),
        ("water_efield_first_clamp_mask", "threshold diagnostic"),
        ("water_efield_first_clamp_scale", "scalar scale"),
        ("water_hbond_counts", "candidate/pass counts"),
        ("water_hbond_nearest", "mixed scalar/identity row"),
        ("delta_scalars", "joint rigid transform of WT and mutant"),
        ("delta_graph", "graph row"),
        ("larsen_sidechain_donor_atoms", "typed donor identity"),
    ])
    def test_directional_companions_declare_exact_invariant_law(
            self, stem, fragment):
        spec = CATALOG[stem]
        expected_frame = (
            "shared_wt_mut_intrinsic"
            if stem in {"delta_scalars", "delta_graph"}
            else "intrinsic_geometry"
        )
        assert spec.coordinate_frame == expected_frame
        assert spec.irreps == "0e"
        assert spec.parity == "even"
        assert spec.tensor_rank == 0
        assert "exact O(3)-invariant" in spec.transformation
        assert fragment in spec.transformation
        assert spec.validity

    @pytest.mark.parametrize("stem", [
        "bs_total_B", "bs_ring_B_field", "hm_ring_B_field",
        "bs_per_type_T1", "hm_per_type_T1",
    ])
    def test_exact_global_axial_vector_law(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == "axial_vector: a'=det(R) R a"

    @pytest.mark.parametrize("stem", [
        "bs_shielding", "hm_shielding", "mc_peptide_co_fixed",
        "mc_peptide_co_rhombic", "mc_peptide_cn_fixed",
        "mc_backbone_other_fixed", "mc_sidechain_co_fixed",
        "mc_sidechain_other_fixed", "mc_disulfide_fixed",
        "mc_aromatic_fixed", "mc_aromatic_bo",
        "mc_backbone_xh_fixed", "mc_sidechain_xh_fixed", "mc_s_h_fixed",
        "mc_nearest_co_T2", "mc_nearest_cn_T2",
        "mopac_mc_nearest_co_T2", "mopac_mc_nearest_cn_T2",
        "mopac_mc_backbone_total", "mopac_mc_sidechain_total",
        "mopac_mc_aromatic_total", "mopac_mc_shielding",
        "sidechain_co_fixed_T2", "eeq_coulomb_efg",
    ])
    def test_exact_global_full_rank2_law(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == "even_rank2: T'=R T R^T"

    @pytest.mark.parametrize("stem", [
        "mc_peptide_co_bo", "mc_peptide_cn_bo", "mc_backbone_other_bo",
        "mc_sidechain_co_bo", "mc_sidechain_other_bo", "mc_disulfide_bo",
        "mc_aromatic_bo",
        "mc_backbone_xh_bo", "mc_sidechain_xh_bo", "mc_s_h_bo",
        "sidechain_co_bo_T2",
    ])
    def test_mopac_weighted_rank2_matches_its_fixed_sibling(self, stem):
        # A scalar Wiberg weight cannot change a tensor's transform law:
        # the _bo channel must carry the same law as its _fixed twin.
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == "even_rank2: T'=R T R^T"
        assert spec.irreps == CATALOG["mc_peptide_co_fixed"].irreps

    @pytest.mark.parametrize("stem", [
        "bs_per_type_T2", "hm_per_type_T2", "eeq_coulomb_efg_backbone",
        "eeq_coulomb_efg_sidechain", "eeq_coulomb_efg_aromatic",
        "water_efg", "water_efg_first",
    ])
    def test_exact_global_native_t2_law(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == (
            "even_rank2_native_T2: reconstruct Cartesian T, apply "
            "T'=R T R^T, then decompose in project-native T2 basis"
        )

    @pytest.mark.parametrize("stem", [
        "aimnet2_efg", "aimnet2_efg_aromatic",
        "aimnet2_efg_backbone", "aimnet2_efg_sidechain",
    ])
    def test_aimnet_native_t2_is_even_rank2(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert "project-native T2 basis" in spec.transformation
        assert spec.parity == "even"
        assert spec.tensor_rank == 2

    @pytest.mark.parametrize("stem", [
        "coulomb_efg", "mopac_coulomb_efg", "coulomb_efg_t2",
        "coulomb_efg_backbone", "coulomb_efg_sidechain",
        "coulomb_efg_aromatic", "mopac_coulomb_efg_backbone",
        "mopac_coulomb_efg_sidechain", "mopac_coulomb_efg_aromatic",
    ])
    def test_coulomb_family_efg_is_even_rank2(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert "T'=R T R^T" in spec.transformation
        assert spec.structural_zero_components == "T0,T1_x,T1_y,T1_z"

    @pytest.mark.parametrize("stem", [
        "apbs_E", "apbs_E_total_diagnostic",
        "apbs_efg", "apbs_efg_total_diagnostic",
    ])
    def test_apbs_finite_grid_physical_laws_are_qualified(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation.startswith("continuum ")
        assert "no exact O(3) law" in spec.transformation
        assert "finite-grid envelope" in spec.transformation

    def test_exact_affine_and_local_component_contracts(self):
        assert CATALOG["pos"].transformation == "cartesian_position: p'=R p+t"
        assert CATALOG["pos"].wrapper is PositionField
        assert CATALOG["mc_nearest_co_midpoint"].transformation == (
            "cartesian_position: p'=R p+t")
        assert CATALOG["mc_nearest_co_midpoint"].wrapper is PositionField
        assert CATALOG["bs_ring_B_cylindrical"].coordinate_frame == (
            "ring_cylindrical_components")
        assert CATALOG["bs_ring_B_cylindrical"].structural_zero_components == (
            "B_phi")
        assert CATALOG["bs_ring_B_cylindrical"].wrapper is np.ndarray
        assert not CATALOG["bs_ring_B_cylindrical"].irreps
        assert CATALOG["piquad_local_tensor"].coordinate_frame == (
            "ring_local_vertex0_gauge")
        assert CATALOG["piquad_local_T2"].coordinate_frame == (
            "ring_local_vertex0_gauge")
        assert "coefficients are invariant under a global proper rotation" in (
            CATALOG["piquad_local_T2"].transformation)
        assert not CATALOG["piquad_local_tensor"].irreps
        assert not CATALOG["piquad_local_tensor"].e3nn_export
        assert CATALOG["piquad_local_tensor"].wrapper is np.ndarray
        assert "col6 is the acute" in CATALOG["ring_contributions"].transformation
        assert "combined tensor-evaluation validity" in (
            CATALOG["piquad_local_geometry"].validity)
        ring_geometry = CATALOG["ring_geometry"]
        assert "ordinary nondegenerate ring" in ring_geometry.transformation
        assert "underdetermined SVD normal" in ring_geometry.transformation
        assert "lab-basis-dependent unit normal" in ring_geometry.validity
        ring = CATALOG["ring_contributions"]
        assert ring.tensor_rank == 0
        assert "cols9:18,18:27,27:36" in ring.tensor_basis
        assert "T0,T1_x,T1_y,T1_z" in ring.tensor_component_order
        assert "cols18:22" in ring.structural_zero_components
        assert not ring.e3nn_export
        scalar = CATALOG["piquad_quad_scalar"]
        assert scalar.coordinate_frame == "intrinsic_geometry"
        assert scalar.transformation.startswith("rotation_invariant scalar")
        assert "no dedicated validity mask" in scalar.validity

    @pytest.mark.parametrize("stem", [
        "larsen_hbond_shielding", "larsen_hbond_1pHB_shielding",
        "larsen_hbond_2pHB_shielding", "larsen_hbond_1pHaB_shielding",
        "larsen_hbond_2pHaB_shielding",
        "larsen_hbond_diagnostic_CB_shielding",
    ])
    def test_procs15_hbond_tensors_keep_type_and_proper_rotation_caveat(
            self, stem):
        spec = CATALOG[stem]
        assert "under proper rotations" in spec.transformation
        assert "no improper-transform contract" in spec.transformation
        assert spec.wrapper is ShieldingTensor
        assert spec.parity == "even"
        assert spec.irreps == CATALOG["bs_shielding"].irreps
        assert spec.e3nn_export == CATALOG["bs_shielding"].e3nn_export

    def test_ideal_l_cb_outputs_have_only_a_proper_rotation_contract(self):
        deviation = CATALOG["cb_deviation"]
        assert deviation.coordinate_frame == "intrinsic_chiral_lookup"
        assert deviation.transformation.startswith(
            "rotation-invariant under proper rotations")
        assert "no improper-transform contract" in deviation.transformation
        assert deviation.parity == "mixed"
        assert deviation.tensor_rank == 0
        assert not deviation.irreps

        residual = CATALOG["cb_residual_vector"]
        assert residual.coordinate_frame == "conformation_cartesian_xyz"
        assert residual.transformation.startswith(
            "polar displacement under proper rotations")
        assert "no single improper-transform parity" in residual.transformation
        assert residual.parity == "mixed"
        assert not residual.irreps

    @pytest.mark.parametrize(
        "stem", ["atom_sasa", "atom_sasa_fraction", "sasa_normal"])
    def test_finite_sasa_stencil_does_not_claim_exact_o3_irreps(self, stem):
        spec = CATALOG[stem]
        assert "finite" in spec.transformation
        assert "no exact O(3) law" in spec.transformation
        assert "recorded" in spec.transformation
        assert spec.parity == "mixed"
        assert not spec.irreps
        assert not spec.e3nn_export

    def test_exact_mixed_block_contracts(self):
        expected_fragments = {
            "ring_contributions": "col5 pseudoscalar",
            "ring_geometry": "cols3:6 affine position",
            "ring_pair_geometry": "cols9:11 are pseudoscalars",
            "spatial_neighbors": "cols2:5 polar unit vector",
            "hbond_pairs_geometry": "cols2:5 polar H-to-O unit vector",
            "piquad_local_frame": "x and y are polar",
            "piquad_local_geometry": "col4 cos_theta is pseudoscalar",
            "sidechain_co_frame": "cols0:3 affine position",
            "water_polarization": "cols0:3 net-water-dipole polar vector",
            "water_hbond_candidates": "cols9:12 water-O affine position",
            "delta_apbs": "cols0:3 polar delta-E",
            "delta_ring_proximity": "z pseudoscalar",
            "mopac_global": "cols1:4 molecular dipole polar vector",
            "mopac_atom_ao_density": "reducible AO-basis matrix",
            "mopac_bond_ao_density": "reducible interatomic AO block",
            "mopac_lmo_occupied_coefficients": "no fixed rowwise O(3)/e3nn law",
            "bonded_energy": "no unconditional improper-transform law",
            "gromacs_energy": "cols23:32 virial and cols32:41 pressure",
            "dssp_backbone": "phi/psi cols0:2 are signed dihedral",
            "dssp_ss8": "physical O(3)-invariant libdssp eight-class one-hot",
            "dssp_chi": "cos/exists invariant; sin is pseudoscalar",
            "larsen_hbond_pairs_geometry": "col2 signed rho",
            "larsen_hbond_pairs": "col18 signed-rho pseudoscalar",
            "larsen_sidechain_donor_candidates": "col10 signed rho",
        }
        for stem, fragment in expected_fragments.items():
            assert fragment in CATALOG[stem].transformation, stem

    def test_dssp_static_boundary_and_companion_contracts_are_explicit(self):
        for stem in ("dssp_backbone", "dssp_ss8", "dssp_ppii"):
            assert "0.001" in CATALOG[stem].transformation, stem
            assert "PDB" in CATALOG[stem].transformation, stem
        energy = CATALOG["dssp_hbond_energy"]
        assert "5e-3 kcal/mol" in energy.transformation
        assert "dssp_hbond_partner_residue_index.npy" in energy.validity
        partner = CATALOG["dssp_hbond_partner_residue_index"]
        assert "different boundary partner" in partner.transformation
        assert "-1 means no mapped partner" in partner.validity
        for stem in ("dssp_observed", "dssp_torsion_valid"):
            assert CATALOG[stem].parity == "even"
            assert "exact rotation/translation/reflection-invariant" in (
                CATALOG[stem].transformation)

    @pytest.mark.parametrize("stem", [
        "dssp_torsion_angle", "dssp_torsion_sin", "omega_actual",
        "omega_deviation", "omega_sin", "aromatic_chi2",
    ])
    def test_homogeneous_signed_geometry_is_explicitly_0o(self, stem):
        spec = CATALOG[stem]
        assert spec.irreps == "0o"
        assert spec.parity == "odd"

    def test_mutation_and_mopac_row_axes_are_the_serialized_axes(self):
        mutation_stems = {
            "delta_shielding", "delta_scalars", "mutation_atom_map",
            "delta_graph", "delta_apbs",
            "delta_ring_proximity", "wt_shielding_diamagnetic",
            "wt_shielding_paramagnetic", "mut_shielding_diamagnetic",
            "mut_shielding_paramagnetic", "delta_shielding_diamagnetic",
            "delta_shielding_paramagnetic",
        }
        assert all(CATALOG[stem].native_axis == "atom"
                   for stem in mutation_stems)
        assert CATALOG["mopac_atomic_orbital_populations"].native_axis == "atom"
        assert CATALOG[
            "mopac_atomic_orbital_population_totals"].native_axis == "atom"

    def test_directional_array_specs_asdict_json_round_trip(self):
        for stem in sorted(_DIRECTIONAL_STEMS):
            payload = asdict(CATALOG[stem])
            # `wrapper` is the one pre-existing non-JSON dataclass member;
            # encode its stable producer class name while exercising every
            # metadata field, including the three freeze-contract additions.
            payload["wrapper"] = payload["wrapper"].__name__
            decoded = json.loads(json.dumps(payload, allow_nan=False))
            assert decoded["stem"] == stem
            assert decoded["coordinate_frame"] == CATALOG[stem].coordinate_frame
            assert decoded["transformation"] == CATALOG[stem].transformation
            assert decoded["validity"] == CATALOG[stem].validity


# ── TopologyGroup load + wrappers ──────────────────────────────────


def _required_identity_npys(out_dir, n_atoms):
    """Minimal required-by-catalog NPYs the loader checks for."""
    write_required_sdk_npys(out_dir, n_atoms)


def _required_calculator_npys(out_dir, n_atoms):
    write_required_sdk_npys(out_dir, n_atoms)


def test_trajectory_frame_loads_parent_invariants_and_full_manifest(tmp_path):
    n_atoms = 8
    n_residues = 2
    n_bonds = 3
    n_aromatic_rings = 2
    n_saturated_rings = 1
    npys_dir = tmp_path / "npys"
    frame_dir = npys_dir / "frame_000123"

    write_required_sdk_npys(
        frame_dir,
        n_atoms=n_atoms,
        n_residues=n_residues,
        n_bonds=n_bonds,
        n_aromatic_rings=n_aromatic_rings,
        n_saturated_rings=n_saturated_rings,
    )
    write_minimal_topology_sidecar(
        npys_dir,
        n_atoms=n_atoms,
        n_residues=n_residues,
        n_bonds=n_bonds,
        n_aromatic_rings=n_aromatic_rings,
        n_saturated_rings=n_saturated_rings,
        protein_id="parent_full_manifest",
    )

    frame_pos = (
        np.arange(n_atoms * 3, dtype=np.float64).reshape(n_atoms, 3) + 0.125
    )
    np.save(frame_dir / "pos.npy", frame_pos)

    category_info = np.zeros(
        n_atoms,
        dtype=[("atom_index", "<i4"), ("parent_marker", "<i4")],
    )
    category_info["atom_index"] = np.arange(n_atoms, dtype=np.int32)
    category_info["parent_marker"] = 2718
    np.save(npys_dir / "atoms_category_info.npy", category_info)

    child_manifest = {
        "schema_version": "1.0",
        "extractor": "nmr_extract",
        "feature_metadata": {
            "mcconnell": {"manifest_owner": "frame_child"},
        },
    }
    (frame_dir / "extraction_manifest.json").write_text(
        json.dumps(child_manifest)
    )

    assert not (npys_dir / "pos.npy").exists()
    assert not (frame_dir / "atoms_category_info.npy").exists()
    for stem in ("residues", "bonds", "rings", "ring_membership"):
        assert not (frame_dir / f"{stem}.npy").exists()
    assert "axis_sizes" not in child_manifest

    protein = load(frame_dir)

    parent_manifest = json.loads(
        (npys_dir / "extraction_manifest.json").read_text()
    )
    parent_category_info = np.load(
        npys_dir / "atoms_category_info.npy", allow_pickle=False
    )

    np.testing.assert_array_equal(protein.pos.data, frame_pos)
    assert protein.n_atoms == parent_manifest["axis_sizes"]["atom"]
    assert protein.category_info.n_atoms == len(parent_category_info)
    np.testing.assert_array_equal(
        protein.category_info.data, parent_category_info
    )
    for stem in ("residues", "bonds", "rings", "ring_membership"):
        parent_table = np.load(npys_dir / f"{stem}.npy", allow_pickle=False)
        typed_table = getattr(protein.topology, stem)
        assert len(typed_table) == len(parent_table)
        np.testing.assert_array_equal(typed_table.data, parent_table)
    assert protein.topology.manifest.data == parent_manifest
    assert protein.topology.manifest.protein_id == "parent_full_manifest"


@pytest.mark.parametrize(
    "missing_stem",
    (
        "atoms_category_info",
        "residues",
        "bonds",
        "rings",
        "ring_membership",
    ),
)
def test_trajectory_frame_missing_parent_topology_names_exact_path(
    tmp_path, missing_stem
):
    npys_dir = tmp_path / "npys"
    frame_dir = npys_dir / "frame_000123"
    write_required_sdk_npys(frame_dir, n_atoms=4, n_residues=1)
    write_minimal_topology_sidecar(npys_dir, n_atoms=4, n_residues=1)
    np.save(
        npys_dir / "atoms_category_info.npy",
        np.zeros(4, dtype=[("atom_index", "<i4")]),
    )
    missing_path = npys_dir / f"{missing_stem}.npy"
    missing_path.unlink()

    with pytest.raises(FileNotFoundError) as exc_info:
        load(frame_dir)

    assert str(missing_path.resolve()) in str(exc_info.value)


def test_numeric_frame_directory_uses_parent_layout_by_convention(tmp_path):
    frame_dir = tmp_path / "frame_000123"
    write_minimal_topology_sidecar(
        tmp_path,
        n_atoms=4,
        n_residues=1,
        protein_id="parent_manifest",
    )
    np.save(
        tmp_path / "atoms_category_info.npy",
        np.zeros(4, dtype=[("atom_index", "<i4")]),
    )
    write_required_sdk_npys(frame_dir, n_atoms=4, n_residues=1)
    write_minimal_topology_sidecar(
        frame_dir,
        n_atoms=4,
        n_residues=1,
        protein_id="child_manifest_must_not_replace_parent",
    )

    protein = load(frame_dir)

    assert protein.n_atoms == 4
    assert protein.topology.residues.n_residues == 1
    assert protein.topology.manifest.protein_id == "parent_manifest"


class TestTopologyLoad:

    @pytest.fixture
    def fake_extraction(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # 4 residues, 4 atoms each. Two synthetic bonds, one aromatic
        # ring with 6 members, one saturated ring with 5 members.
        n_residues = N_ATOMS // 4
        residues = np.zeros(n_residues, dtype=_RESIDUES_DTYPE)
        for i in range(n_residues):
            residues[i]["residue_index"] = i
            residues[i]["residue_number"] = i + 1
            residues[i]["atom_count"] = 4
            residues[i]["prev_residue_index"] = i - 1 if i > 0 else -1
            residues[i]["next_residue_index"] = i + 1 if i + 1 < n_residues else -1
        np.save(tmp_path / "residues.npy", residues)

        bonds = np.zeros(2, dtype=_BONDS_DTYPE)
        bonds[0] = (0, 0, 1, 0, 0, 0, 0, 0, 1)   # bond_index, a, b, order=Single, cat=PeptideCO, rotate, arom, peptide, backbone
        bonds[1] = (1, 1, 2, 4, 1, 0, 0, 1, 1)   # peptide bond
        np.save(tmp_path / "bonds.npy", bonds)

        rings = np.zeros(2, dtype=_RINGS_DTYPE)
        rings[0] = (0, 0, 0, 6, 0, 0, 1, 100, -1)  # aromatic PheBenzene
        rings[1] = (1, 1, 8, 5, 0, 0, 2, 200, -1)  # saturated Pro
        np.save(tmp_path / "rings.npy", rings)

        memb = np.zeros(11, dtype=_RING_MEMBERSHIP_DTYPE)
        for k in range(6):
            memb[k] = (0, k, k, 1, 0, 0)
        for k in range(5):
            memb[6 + k] = (1, 6 + k, k, 1, 0, 0)
        np.save(tmp_path / "ring_membership.npy", memb)

        manifest = {
            "schema_version": "1.0",
            "extractor": "nmr_extract",
            "generated_at_utc": "2026-05-13T00:00:00Z",
            "protein_id": "test_fixture",
            "topology": {"has_atom_semantic": True,
                         "has_ff_atom_types": False,
                         "has_ff_mass": False},
            "axis_sizes": {
                "atom": N_ATOMS, "residue": n_residues, "bond": 2,
                "aromatic_ring": 1, "saturated_ring": 1, "ring": 2,
                "ring_membership": 11,
            },
            "axis_alignment": {"atom": "rows align with atom_index"},
        }
        (tmp_path / "extraction_manifest.json").write_text(json.dumps(manifest))
        return tmp_path

    def test_residues_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        r = p.topology.residues
        assert isinstance(r, Residues)
        assert r.n_residues == 4
        # Total atom_count must sum to atom axis size (validation invariant).
        assert int(r.atom_count.sum()) == N_ATOMS
        # First residue has prev = -1 (chain start).
        assert r.prev_residue_index[0] == -1
        # Last residue has next = -1 (chain end).
        assert r.next_residue_index[-1] == -1

    def test_protein_topology_is_a_topology_group(self, fake_extraction):
        p = load(fake_extraction)
        assert isinstance(p.topology, TopologyGroup)

    def test_bonds_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        b = p.topology.bonds
        assert isinstance(b, Bonds)
        assert b.n_bonds == 2
        assert b.atom_index_a[0] == 0 and b.atom_index_b[0] == 1
        assert b.atom_index_a[1] == 1 and b.atom_index_b[1] == 2
        assert b.is_peptide[0] == False  # bond 0 is Single
        assert b.is_peptide[1] == True   # bond 1 is Peptide
        assert b.is_backbone[0] == True
        assert b.is_backbone[1] == True

    def test_rings_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        r = p.topology.rings
        assert isinstance(r, Rings)
        assert r.n_rings == 2
        assert r.is_aromatic[0]
        assert r.is_saturated[1]
        assert r.atom_count[0] == 6
        assert r.atom_count[1] == 5
        assert r.parent_residue_number[0] == 100
        assert r.parent_residue_number[1] == 200
        assert r.fused_partner_ring_id[0] == -1

    def test_ring_membership_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        m = p.topology.ring_membership
        assert isinstance(m, RingMembership)
        assert m.n_rows == 11
        # First six rows reference ring_id 0 (aromatic).
        assert np.all(m.ring_id[:6] == 0)
        # Next five reference ring 1 (saturated).
        assert np.all(m.ring_id[6:] == 1)
        # All rows are vertices in our model.
        assert np.all(m.is_vertex)

    def test_manifest_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        man = p.topology.manifest
        assert isinstance(man, ExtractionManifest)
        assert man.schema_version == "1.0"
        assert man.protein_id == "test_fixture"
        assert man.has_atom_semantic()
        assert man.axis_size("atom") == N_ATOMS
        assert man.axis_size("bond") == 2
        assert man.axis_size("ring") == 2
        assert man.axis_size("ring_membership") == 11


class TestStrictLoad:
    """Loader rejects extractions without the topology sidecar."""

    def test_missing_bonds_npy_fails(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # No residues / bonds / rings / ring_membership / manifest.
        with pytest.raises(FileNotFoundError):
            load(tmp_path)

    def test_missing_manifest_fails(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # residues/bonds/rings/ring_membership present; manifest absent.
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        (tmp_path / "extraction_manifest.json").unlink()
        with pytest.raises(FileNotFoundError, match="extraction_manifest"):
            load(tmp_path)


class TestValidationInvariants:
    """SDK enforces the codex first-pass validation gates at load time."""

    def test_axis_size_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # Write a topology sidecar with a manifest that lies about axis sizes.
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # Re-emit a corrupted manifest claiming a different atom count.
        manifest_path = tmp_path / "extraction_manifest.json"
        man = json.loads(manifest_path.read_text())
        man["axis_sizes"]["atom"] = N_ATOMS + 7   # lie
        manifest_path.write_text(json.dumps(man))
        with pytest.raises(ValueError, match="atom axis"):
            load(tmp_path)

    def test_bond_endpoint_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS, n_bonds=2)
        # Corrupt bonds.npy with an out-of-range endpoint.
        bonds = np.load(tmp_path / "bonds.npy")
        bonds[0]["atom_index_a"] = N_ATOMS + 100  # invalid
        np.save(tmp_path / "bonds.npy", bonds)
        with pytest.raises(ValueError, match="bonds.npy"):
            load(tmp_path)

    def test_ring_membership_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=1)
        # Membership rows are empty by default; need to add a corrupt one.
        memb = np.zeros(1, dtype=_RING_MEMBERSHIP_DTYPE)
        memb[0]["ring_id"] = 99  # invalid (only 1 ring exists)
        memb[0]["atom_index"] = 0
        memb[0]["is_vertex"] = 1
        np.save(tmp_path / "ring_membership.npy", memb)
        man_path = tmp_path / "extraction_manifest.json"
        man = json.loads(man_path.read_text())
        man["axis_sizes"]["ring_membership"] = 1
        man_path.write_text(json.dumps(man))
        with pytest.raises(ValueError, match="ring_id"):
            load(tmp_path)

    def test_residue_atom_count_sum_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # Corrupt residues.npy so atom_count sum != N_ATOMS.
        residues = np.load(tmp_path / "residues.npy")
        residues[0]["atom_count"] = 999
        np.save(tmp_path / "residues.npy", residues)
        with pytest.raises(ValueError, match="atom_count sums"):
            load(tmp_path)

    def test_residue_row_identity_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        residues = np.load(tmp_path / "residues.npy")
        residues[1]["residue_index"] = 99
        np.save(tmp_path / "residues.npy", residues)
        with pytest.raises(ValueError, match="residue_index"):
            load(tmp_path)

    def test_bond_row_identity_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS, n_bonds=2)
        bonds = np.load(tmp_path / "bonds.npy")
        bonds[1]["bond_index"] = 99
        np.save(tmp_path / "bonds.npy", bonds)
        with pytest.raises(ValueError, match="bond_index"):
            load(tmp_path)

    def test_ring_row_identity_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=2)
        rings = np.load(tmp_path / "rings.npy")
        rings[1]["ring_id"] = 99
        np.save(tmp_path / "rings.npy", rings)
        with pytest.raises(ValueError, match="ring_id"):
            load(tmp_path)

    def test_atom_residue_counts_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # Keep all residue references in range but redistribute one atom
        # onto residue 0; residues.npy atom_count no longer matches.
        residue_index = np.load(tmp_path / "residue_index.npy")
        residue_index[-1] = 0
        np.save(tmp_path / "residue_index.npy", residue_index)
        with pytest.raises(ValueError, match="residue_index.npy counts"):
            load(tmp_path)

    def test_atom_residue_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        residue_index = np.load(tmp_path / "residue_index.npy")
        residue_index[0] = 99
        np.save(tmp_path / "residue_index.npy", residue_index)
        with pytest.raises(ValueError, match="residue_index.npy references"):
            load(tmp_path)

    def test_ring_parent_residue_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=1)
        rings = np.load(tmp_path / "rings.npy")
        rings[0]["parent_residue_index"] = 99
        np.save(tmp_path / "rings.npy", rings)
        with pytest.raises(ValueError, match="parent_residue_index"):
            load(tmp_path)

    def test_ring_native_axis_noncontiguous_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=2)
        rings = np.load(tmp_path / "rings.npy")
        rings[1]["native_axis_index"] = 3
        np.save(tmp_path / "rings.npy", rings)
        with pytest.raises(ValueError, match="native_axis_index"):
            load(tmp_path)

    def test_ring_atom_count_membership_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=1)
        rings = np.load(tmp_path / "rings.npy")
        rings[0]["atom_count"] = 2
        np.save(tmp_path / "rings.npy", rings)

        memb = np.zeros(1, dtype=_RING_MEMBERSHIP_DTYPE)
        memb[0]["ring_id"] = 0
        memb[0]["atom_index"] = 0
        memb[0]["is_vertex"] = 1
        np.save(tmp_path / "ring_membership.npy", memb)
        man_path = tmp_path / "extraction_manifest.json"
        man = json.loads(man_path.read_text())
        man["axis_sizes"]["ring_membership"] = 1
        man_path.write_text(json.dumps(man))

        with pytest.raises(ValueError, match="atom_count does not match"):
            load(tmp_path)


# The 6 new fields on atoms_category_info are exercised end-to-end on
# the C++ side at tests/test_category_info_projection.cpp::WriteFeaturesEmitsNpy
# (header check) and the round-trip from real extraction is implicit:
# accessing ``info.chain_id`` on an NPY without the column raises
# ValueError from numpy structured-array access (loud fail), which is
# the contract.
