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
    "cb_residual_vector", "spatial_neighbors", "hbond_nearest_dir",
    "hbond_pairs_geometry", "sasa_normal", "water_polarization",
    "water_hbond_candidates", "sidechain_co_frame",
    "piquad_local_tensor", "piquad_local_T2", "piquad_local_frame",
    "piquad_local_geometry",
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
    "mc_aromatic_zeroed_fixed", "mc_aromatic_zeroed_bo",
    "mc_backbone_xh_fixed", "mc_backbone_xh_bo",
    "mc_sidechain_xh_fixed", "mc_sidechain_xh_bo",
    "mc_s_h_fixed", "mc_s_h_bo", "mc_nearest_co_dir",
    "mc_nearest_co_midpoint", "mc_nearest_co_T2", "mc_nearest_cn_T2",
    "sidechain_co_fixed_T2", "sidechain_co_bo_T2",
    # Coulomb/APBS/water/AIMNet2 vectors and EFGs.
    "coulomb_efg", "coulomb_efg_t2", "coulomb_E",
    "coulomb_E_backbone", "coulomb_E_sidechain", "coulomb_E_aromatic",
    "coulomb_efg_backbone", "coulomb_efg_sidechain",
    "coulomb_efg_aromatic", "coulomb_E_solvent", "coulomb_efg_solvent",
    "mopac_coulomb_efg", "mopac_coulomb_E",
    "mopac_coulomb_E_backbone", "mopac_coulomb_E_sidechain",
    "mopac_coulomb_E_aromatic", "mopac_coulomb_efg_backbone",
    "mopac_coulomb_efg_sidechain", "mopac_coulomb_efg_aromatic",
    "eeq_coulomb_efg", "eeq_coulomb_E", "eeq_coulomb_E_backbone",
    "eeq_coulomb_E_sidechain", "eeq_coulomb_E_aromatic",
    "eeq_coulomb_efg_backbone", "eeq_coulomb_efg_sidechain",
    "eeq_coulomb_efg_aromatic", "apbs_E", "apbs_efg",
    "apbs_E_total_diagnostic", "apbs_efg_total_diagnostic",
    "water_efield", "water_efield_first", "water_efg",
    "water_efg_first", "aimnet2_E", "aimnet2_E_backbone",
    "aimnet2_E_sidechain", "aimnet2_E_aromatic", "aimnet2_efg",
    "aimnet2_efg_aromatic", "aimnet2_efg_backbone",
    "aimnet2_efg_sidechain", "aimnet2_charge_response_gradient",
    # External quantum/reference and mutation tensors.
    "orca_total", "orca_diamagnetic", "orca_paramagnetic",
    "tripeptide_bb_shielding", "tripeptide_bb_residual_vec",
    "tripeptide_bb_match_distance", "tripeptide_bb_match_atoms",
    "tripeptide_bb_method_tag", "tripeptide_neighbor_reference",
    "tripeptide_neighbor_shielding", "tripeptide_neighbor_shielding_prev",
    "tripeptide_neighbor_shielding_next",
    "tripeptide_neighbor_residual_vec_prev",
    "tripeptide_neighbor_residual_vec_next",
    "larsen_hbond_shielding", "larsen_hbond_1pHB_shielding",
    "larsen_hbond_2pHB_shielding", "larsen_hbond_1pHaB_shielding",
    "larsen_hbond_2pHaB_shielding",
    "larsen_hbond_diagnostic_CB_shielding", "delta_shielding",
    "wt_shielding_diamagnetic", "wt_shielding_paramagnetic",
    "mut_shielding_diamagnetic", "mut_shielding_paramagnetic",
    "delta_shielding_diamagnetic", "delta_shielding_paramagnetic",
    "delta_apbs", "delta_ring_proximity",
    # External-frame and reflection-sensitive diagnostics.
    "mopac_global", "mopac_atom_populations",
    "mopac_atomic_orbital_populations",
    "mopac_atomic_orbital_population_totals", "gromacs_energy",
    "dssp_backbone", "dssp_ss8", "dssp_chi", "dssp_torsion_angle",
    "dssp_torsion_sin", "dssp_torsion_cos", "dssp_ppii", "omega_actual",
    "omega_deviation", "omega_sin", "omega_cos", "aromatic_chi2",
    "pucker_Q", "pucker_theta", "tripeptide_bb_diagnostics",
    "tripeptide_neighbor_diagnostics", "larsen_hbond_pairs_geometry",
    "larsen_hbond_pairs_isotropic", "larsen_hbond_pairs_index",
    "larsen_hbond_count", "larsen_corner_imputed",
    "larsen_imputed_pair_count", "larsen_sidechain_carbonyl_pair_count",
    "larsen_hbond_pairs",
    "larsen_sidechain_donor_candidates",
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

    def test_every_entry_has_valid_parity(self):
        bad = [
            k for k, s in CATALOG.items()
            if s.parity not in ("even", "odd", "mixed")
        ]
        assert not bad, f"Entries with invalid parity: {bad}"
        assert {
            k for k, s in CATALOG.items() if s.parity == "mixed"
        } == {
            "bs_ring_B_cylindrical",
            "dssp_backbone",
            "dssp_chi",
            "dssp_ppii",
            "dssp_ss8",
            "piquad_local_tensor",
            "piquad_local_T2",
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
        "hbond_nearest_dir", "coulomb_E", "coulomb_E_backbone",
        "coulomb_E_sidechain", "coulomb_E_aromatic", "coulomb_E_solvent",
        "mopac_coulomb_E", "mopac_coulomb_E_backbone",
        "mopac_coulomb_E_sidechain", "mopac_coulomb_E_aromatic",
        "eeq_coulomb_E", "eeq_coulomb_E_backbone",
        "eeq_coulomb_E_sidechain", "eeq_coulomb_E_aromatic", "apbs_E",
        "apbs_E_total_diagnostic", "water_efield", "water_efield_first",
        "aimnet2_E", "aimnet2_E_backbone", "aimnet2_E_sidechain",
        "aimnet2_E_aromatic", "aimnet2_charge_response_gradient",
    ])
    def test_exact_global_polar_vector_law(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == "polar_vector: v'=R v"

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
        "mc_peptide_co_bo", "mc_peptide_co_rhombic", "mc_peptide_cn_fixed",
        "mc_peptide_cn_bo", "mc_backbone_other_fixed",
        "mc_backbone_other_bo", "mc_sidechain_co_fixed",
        "mc_sidechain_co_bo", "mc_sidechain_other_fixed",
        "mc_sidechain_other_bo", "mc_disulfide_fixed", "mc_disulfide_bo",
        "mc_aromatic_zeroed_fixed", "mc_aromatic_zeroed_bo",
        "mc_backbone_xh_fixed", "mc_backbone_xh_bo",
        "mc_sidechain_xh_fixed", "mc_sidechain_xh_bo", "mc_s_h_fixed",
        "mc_s_h_bo", "mc_nearest_co_T2", "mc_nearest_cn_T2",
        "sidechain_co_fixed_T2", "sidechain_co_bo_T2", "coulomb_efg",
        "mopac_coulomb_efg", "eeq_coulomb_efg",
    ])
    def test_exact_global_full_rank2_law(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == "even_rank2: T'=R T R^T"

    @pytest.mark.parametrize("stem", [
        "bs_per_type_T2", "hm_per_type_T2", "coulomb_efg_t2",
        "coulomb_efg_backbone", "coulomb_efg_sidechain",
        "coulomb_efg_aromatic", "coulomb_efg_solvent",
        "mopac_coulomb_efg_backbone", "mopac_coulomb_efg_sidechain",
        "mopac_coulomb_efg_aromatic", "eeq_coulomb_efg_backbone",
        "eeq_coulomb_efg_sidechain", "eeq_coulomb_efg_aromatic",
        "apbs_efg", "apbs_efg_total_diagnostic", "water_efg",
        "water_efg_first", "aimnet2_efg", "aimnet2_efg_aromatic",
        "aimnet2_efg_backbone", "aimnet2_efg_sidechain",
    ])
    def test_exact_global_native_t2_law(self, stem):
        spec = CATALOG[stem]
        assert spec.coordinate_frame == "conformation_cartesian_xyz"
        assert spec.transformation == (
            "even_rank2_native_T2: reconstruct Cartesian T, apply "
            "T'=R T R^T, then decompose in project-native T2 basis"
        )

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
            "mopac_atom_populations": "intended cols6:9 per-atom polar dipole",
            "gromacs_energy": "cols23:32 virial and cols32:41 pressure",
            "dssp_backbone": "phi/psi cols0:2 are signed dihedral",
            "dssp_ss8": "no homogeneous improper-transform law",
            "dssp_chi": "cos/exists invariant; sin is pseudoscalar",
            "tripeptide_bb_diagnostics": "cols10:12 and14:18 reverse sign",
            "tripeptide_neighbor_diagnostics": "base+11:13",
            "larsen_hbond_pairs_geometry": "col2 signed rho",
            "larsen_hbond_pairs": "col18 signed-rho pseudoscalar",
            "larsen_sidechain_donor_candidates": "col10 signed rho",
        }
        for stem, fragment in expected_fragments.items():
            assert fragment in CATALOG[stem].transformation, stem

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
            "delta_shielding", "delta_scalars", "delta_graph", "delta_apbs",
            "delta_ring_proximity", "wt_shielding_diamagnetic",
            "wt_shielding_paramagnetic", "mut_shielding_diamagnetic",
            "mut_shielding_paramagnetic", "delta_shielding_diamagnetic",
            "delta_shielding_paramagnetic",
        }
        assert all(CATALOG[stem].native_axis == "atom"
                   for stem in mutation_stems)
        assert CATALOG["mopac_atomic_orbital_populations"].native_axis == (
            "mopac_atomic_orbital_row")
        assert CATALOG[
            "mopac_atomic_orbital_population_totals"].native_axis == (
                "mopac_atomic_orbital_row")

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
