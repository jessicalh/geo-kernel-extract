# CODEX Wire-Up Coverage Report

Generated from the emitted `catalog_coverage.json` manifests, not from the snapshot or previous row output.

## Outputs

- 1P9J trajectory output: `/shared/2026Thesis/shielding-calcsets/data/row-design/full_wireup_20260611/1p9j`
- 720 static output: `/shared/2026Thesis/shielding-calcsets/data/row-design/full_wireup_20260611/720`
- Exact per-component counts are in each dataset `catalog_coverage.json` and the flattened `catalog_sidecar_support.csv`; this report uses compact notation where all components share the same count.

## Scope

- Scoped producer catalog fields: 122
- Excluded by contract: Larsen/Larsen H-bond, tripeptide, mutation deltas, rediscover, legacy duplicate groups, and dead imposed-geometry tensors absent from the generated catalog.
- 1P9J rows: 635,346 = 846 atoms x 751 DFT frames over the 1501-frame trajectory stride-2 mapping.
- 720 rows: 475,116 across 720 discovered static pose directories.

## Manifest Summary

- 1P9J: field_count=122; atom_row_sidecar=59, native_axis_sidecar=19, row_columns=39, structured_native_sidecar=5; row_reductions=19
- 720: field_count=122; atom_row_sidecar=47, atom_row_sidecar_absent=12, native_axis_sidecar=18, native_axis_sidecar_absent=1, row_columns=39, structured_native_sidecar=5; row_reductions=18

## Coverage Table

| Field | Group | Axis | 1P9J representation | 1P9J counts | 1P9J reduction | 720 representation | 720 counts | 720 reduction |
|---|---|---|---|---:|---|---|---:|---|
| pos | identity | atom | atom_row_sidecar: row_design_field_pos.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_pos.npy | 475116 x3 |  |
| element | identity | atom | row columns: catalog_element | 635346 |  | row columns: catalog_element | 475116 |  |
| residue_index | identity | atom | row columns: catalog_residue_index | 635346 |  | row columns: catalog_residue_index | 475116 |  |
| residue_type | identity | atom | row columns: catalog_residue_type | 635346 |  | row columns: catalog_residue_type | 475116 |  |
| ring_contributions | identity | ring_contribution_pair | native_axis_sidecar: row_design_field_ring_contributions_native.npy | 3804062 x40 | row_design_field_ring_contributions_reduction.npy (401 comps; min=0 max=635346; common: 635283 x80, 0 x80, 227217 x40, 608958 x40, 256575 x40; exact in manifest/CSV) | native_axis_sidecar: row_design_field_ring_contributions_native.npy | 1025753 x40 | row_design_field_ring_contributions_reduction.npy (401 comps; min=0 max=475116; common: 385542 x80, 0 x80, 204161 x40, 199743 x40, 79020 x40; exact in manifest/CSV) |
| ring_direction_to_center | identity | ring_contribution_pair | native_axis_sidecar: row_design_field_ring_direction_to_center_native.npy | 3804062 x3 | row_design_field_ring_direction_to_center_reduction.npy (31 comps; min=0 max=635346; common: 635283 x6, 0 x6, 227217 x3, 608958 x3, 256575 x3; exact in manifest/CSV) | native_axis_sidecar: row_design_field_ring_direction_to_center_native.npy | 1025753 x3 | row_design_field_ring_direction_to_center_reduction.npy (31 comps; min=0 max=475116; common: 385542 x6, 0 x6, 204161 x3, 199743 x3, 79020 x3; exact in manifest/CSV) |
| ring_geometry | identity | aromatic_ring | native_axis_sidecar: row_design_field_ring_geometry_native.npy | 11265 x10 | row_design_field_ring_geometry_reduction.npy (51819 x10) | native_axis_sidecar: row_design_field_ring_geometry_native.npy | 3015 x10 | row_design_field_ring_geometry_reduction.npy (15034 x10) |
| enrichment_role | enrichment | atom | row columns: catalog_enrichment_role | 635346 |  | row columns: catalog_enrichment_role | 475116 |  |
| enrichment_hybridisation | enrichment | atom | row columns: catalog_enrichment_hybridisation | 635346 |  | row columns: catalog_enrichment_hybridisation | 475116 |  |
| enrichment_flags | enrichment | atom | row columns (8): catalog_enrichment_flags_0 ... catalog_enrichment_flags_7 | 635346 x8 |  | row columns (8): catalog_enrichment_flags_0 ... catalog_enrichment_flags_7 | 475116 x8 |  |
| ff_partial_charge | charge_assignment | atom | row columns: catalog_ff_partial_charge | 635346 |  | row columns: catalog_ff_partial_charge | 475116 |  |
| ff_pb_radius | charge_assignment | atom | row columns: catalog_ff_pb_radius | 635346 |  | row columns: catalog_ff_pb_radius | 475116 |  |
| bs_shielding | biot_savart | atom | atom_row_sidecar: row_design_field_bs_shielding.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_bs_shielding.npy | 475116 x9 |  |
| bs_per_type_T0 | biot_savart | atom | row columns (8): catalog_bs_per_type_T0_0 ... catalog_bs_per_type_T0_7 | 635346 x8 |  | row columns (8): catalog_bs_per_type_T0_0 ... catalog_bs_per_type_T0_7 | 475116 x8 |  |
| bs_per_type_T1 | biot_savart | atom | atom_row_sidecar: row_design_field_bs_per_type_T1.npy | 635346 x24 |  | atom_row_sidecar: row_design_field_bs_per_type_T1.npy | 475116 x24 |  |
| bs_per_type_T2 | biot_savart | atom | atom_row_sidecar: row_design_field_bs_per_type_T2.npy | 635346 x40 |  | atom_row_sidecar: row_design_field_bs_per_type_T2.npy | 475116 x40 |  |
| bs_total_B | biot_savart | atom | atom_row_sidecar: row_design_field_bs_total_B.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_bs_total_B.npy | 475116 x3 |  |
| bs_ring_B_field | biot_savart | ring_contribution_pair | native_axis_sidecar: row_design_field_bs_ring_B_field_native.npy | 3804062 x3 | row_design_field_bs_ring_B_field_reduction.npy (31 comps; min=0 max=635346; common: 635283 x6, 0 x6, 227217 x3, 608958 x3, 256575 x3; exact in manifest/CSV) | native_axis_sidecar: row_design_field_bs_ring_B_field_native.npy | 1025753 x3 | row_design_field_bs_ring_B_field_reduction.npy (31 comps; min=0 max=475116; common: 385542 x6, 0 x6, 204161 x3, 199743 x3, 79020 x3; exact in manifest/CSV) |
| bs_ring_B_cylindrical | biot_savart | ring_contribution_pair | native_axis_sidecar: row_design_field_bs_ring_B_cylindrical_native.npy | 3804062 x3 | row_design_field_bs_ring_B_cylindrical_reduction.npy (31 comps; min=0 max=635346; common: 635283 x6, 0 x6, 227217 x3, 608958 x3, 256575 x3; exact in manifest/CSV) | native_axis_sidecar: row_design_field_bs_ring_B_cylindrical_native.npy | 1025753 x3 | row_design_field_bs_ring_B_cylindrical_reduction.npy (31 comps; min=0 max=475116; common: 385542 x6, 0 x6, 204161 x3, 199743 x3, 79020 x3; exact in manifest/CSV) |
| bs_ring_counts | biot_savart | atom | row columns: catalog_bs_ring_counts_0, catalog_bs_ring_counts_1, catalog_bs_ring_counts_2, catalog_bs_ring_counts_3 | 635346 x4 |  | row columns: catalog_bs_ring_counts_0, catalog_bs_ring_counts_1, catalog_bs_ring_counts_2, catalog_bs_ring_counts_3 | 475116 x4 |  |
| hm_shielding | haigh_mallion | atom | atom_row_sidecar: row_design_field_hm_shielding.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_hm_shielding.npy | 475116 x9 |  |
| hm_per_type_T0 | haigh_mallion | atom | row columns (8): catalog_hm_per_type_T0_0 ... catalog_hm_per_type_T0_7 | 635346 x8 |  | row columns (8): catalog_hm_per_type_T0_0 ... catalog_hm_per_type_T0_7 | 475116 x8 |  |
| hm_per_type_T1 | haigh_mallion | atom | atom_row_sidecar: row_design_field_hm_per_type_T1.npy | 635346 x24 |  | atom_row_sidecar: row_design_field_hm_per_type_T1.npy | 475116 x24 |  |
| hm_per_type_T2 | haigh_mallion | atom | atom_row_sidecar: row_design_field_hm_per_type_T2.npy | 635346 x40 |  | atom_row_sidecar: row_design_field_hm_per_type_T2.npy | 475116 x40 |  |
| hm_ring_B_field | haigh_mallion | ring_contribution_pair | native_axis_sidecar: row_design_field_hm_ring_B_field_native.npy | 3804062 x3 | row_design_field_hm_ring_B_field_reduction.npy (31 comps; min=0 max=635346; common: 635283 x6, 0 x6, 227217 x3, 608958 x3, 256575 x3; exact in manifest/CSV) | native_axis_sidecar: row_design_field_hm_ring_B_field_native.npy | 1025753 x3 | row_design_field_hm_ring_B_field_reduction.npy (31 comps; min=0 max=475116; common: 385542 x6, 0 x6, 204161 x3, 199743 x3, 79020 x3; exact in manifest/CSV) |
| pq_per_type_T0 | pi_quadrupole | atom | row columns (8): catalog_pq_per_type_T0_0 ... catalog_pq_per_type_T0_7 | 635346 x8 |  | row columns (8): catalog_pq_per_type_T0_0 ... catalog_pq_per_type_T0_7 | 475116 x8 |  |
| piquad_quad_scalar | pi_quadrupole | ring_contribution_pair | native_axis_sidecar: row_design_field_piquad_quad_scalar_native.npy | 3804062 | row_design_field_piquad_quad_scalar_reduction.npy ([635346, 635283, 635283, 227217, 608958, 256575, 249788, 254329, 544987, 0, 0]) | native_axis_sidecar: row_design_field_piquad_quad_scalar_native.npy | 1025753 | row_design_field_piquad_quad_scalar_reduction.npy ([475116, 385542, 385542, 204161, 199743, 79020, 80757, 78994, 0, 0, 150868]) |
| disp_per_type_T0 | dispersion | atom | row columns (8): catalog_disp_per_type_T0_0 ... catalog_disp_per_type_T0_7 | 635346 x8 |  | row columns (8): catalog_disp_per_type_T0_0 ... catalog_disp_per_type_T0_7 | 475116 x8 |  |
| mc_peptide_co_fixed | mcconnell | atom | atom_row_sidecar: row_design_field_mc_peptide_co_fixed.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_peptide_co_fixed.npy | 475116 x9 |  |
| mc_peptide_co_bo | mcconnell | atom | atom_row_sidecar: row_design_field_mc_peptide_co_bo.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_peptide_co_bo.npy | 475116 x9 |  |
| mc_peptide_co_rhombic | mcconnell | atom | atom_row_sidecar: row_design_field_mc_peptide_co_rhombic.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_peptide_co_rhombic.npy | 475116 x9 |  |
| mc_peptide_cn_fixed | mcconnell | atom | atom_row_sidecar: row_design_field_mc_peptide_cn_fixed.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_peptide_cn_fixed.npy | 475116 x9 |  |
| mc_peptide_cn_bo | mcconnell | atom | atom_row_sidecar: row_design_field_mc_peptide_cn_bo.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_peptide_cn_bo.npy | 475116 x9 |  |
| mc_backbone_other_fixed | mcconnell | atom | atom_row_sidecar: row_design_field_mc_backbone_other_fixed.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_backbone_other_fixed.npy | 475116 x9 |  |
| mc_backbone_other_bo | mcconnell | atom | atom_row_sidecar: row_design_field_mc_backbone_other_bo.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_backbone_other_bo.npy | 475116 x9 |  |
| mc_sidechain_co_fixed | mcconnell | atom | atom_row_sidecar: row_design_field_mc_sidechain_co_fixed.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_sidechain_co_fixed.npy | 475116 x9 |  |
| mc_sidechain_co_bo | mcconnell | atom | atom_row_sidecar: row_design_field_mc_sidechain_co_bo.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_sidechain_co_bo.npy | 475116 x9 |  |
| mc_sidechain_other_fixed | mcconnell | atom | atom_row_sidecar: row_design_field_mc_sidechain_other_fixed.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_sidechain_other_fixed.npy | 475116 x9 |  |
| mc_sidechain_other_bo | mcconnell | atom | atom_row_sidecar: row_design_field_mc_sidechain_other_bo.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_sidechain_other_bo.npy | 475116 x9 |  |
| mc_disulfide_fixed | mcconnell | atom | atom_row_sidecar: row_design_field_mc_disulfide_fixed.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_disulfide_fixed.npy | 475116 x9 |  |
| mc_disulfide_bo | mcconnell | atom | atom_row_sidecar: row_design_field_mc_disulfide_bo.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_disulfide_bo.npy | 475116 x9 |  |
| mc_aromatic_zeroed_fixed | mcconnell | atom | atom_row_sidecar: row_design_field_mc_aromatic_zeroed_fixed.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_aromatic_zeroed_fixed.npy | 475116 x9 |  |
| mc_aromatic_zeroed_bo | mcconnell | atom | atom_row_sidecar: row_design_field_mc_aromatic_zeroed_bo.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_aromatic_zeroed_bo.npy | 475116 x9 |  |
| mc_nearfield_counts | mcconnell | atom | row columns: catalog_mc_nearfield_counts_0, catalog_mc_nearfield_counts_1 | 635346 x2 |  | row columns: catalog_mc_nearfield_counts_0, catalog_mc_nearfield_counts_1 | 475116 x2 |  |
| mc_nearest_co_dir | mcconnell | atom | atom_row_sidecar: row_design_field_mc_nearest_co_dir.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_mc_nearest_co_dir.npy | 475116 x3 |  |
| mc_nearest_co_midpoint | mcconnell | atom | atom_row_sidecar: row_design_field_mc_nearest_co_midpoint.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_mc_nearest_co_midpoint.npy | 475116 x3 |  |
| mc_nearest_co_T2 | mcconnell | atom | atom_row_sidecar: row_design_field_mc_nearest_co_T2.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_nearest_co_T2.npy | 475116 x9 |  |
| mc_nearest_cn_T2 | mcconnell | atom | atom_row_sidecar: row_design_field_mc_nearest_cn_T2.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mc_nearest_cn_T2.npy | 475116 x9 |  |
| coulomb_efg | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_efg.npy | 635346 x9 |  | atom_row_sidecar_absent | 0 x9 |  |
| coulomb_E | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_E.npy | 635346 x3 |  | atom_row_sidecar_absent | 0 x3 |  |
| coulomb_E_backbone | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_E_backbone.npy | 635346 x3 |  | atom_row_sidecar_absent | 0 x3 |  |
| coulomb_E_sidechain | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_E_sidechain.npy | 635346 x3 |  | atom_row_sidecar_absent | 0 x3 |  |
| coulomb_E_aromatic | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_E_aromatic.npy | 635346 x3 |  | atom_row_sidecar_absent | 0 x3 |  |
| coulomb_efg_backbone | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_efg_backbone.npy | 635346 x5 |  | atom_row_sidecar_absent | 0 x5 |  |
| coulomb_efg_sidechain | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_efg_sidechain.npy | 635346 x5 |  | atom_row_sidecar_absent | 0 x5 |  |
| coulomb_efg_aromatic | coulomb | atom | atom_row_sidecar: row_design_field_coulomb_efg_aromatic.npy | 635346 x5 |  | atom_row_sidecar_absent | 0 x5 |  |
| coulomb_scalars | coulomb | atom | row columns: catalog_coulomb_scalars_0, catalog_coulomb_scalars_1, catalog_coulomb_scalars_2, catalog_coulomb_scalars_3 | 635346 x4 |  | row columns: catalog_coulomb_scalars_0, catalog_coulomb_scalars_1, catalog_coulomb_scalars_2, catalog_coulomb_scalars_3 | 0 x4 |  |
| coulomb_aromatic_E_proj | coulomb | atom | row columns: catalog_coulomb_aromatic_E_proj | 635346 |  | row columns: catalog_coulomb_aromatic_E_proj | 0 |  |
| coulomb_aromatic_n_src | coulomb | atom | row columns: catalog_coulomb_aromatic_n_src | 635346 |  | row columns: catalog_coulomb_aromatic_n_src | 0 |  |
| hbond_scalars | hbond | atom | row columns: catalog_hbond_scalars_0, catalog_hbond_scalars_1, catalog_hbond_scalars_2, catalog_hbond_scalars_3 | 635346 x4 |  | row columns: catalog_hbond_scalars_0, catalog_hbond_scalars_1, catalog_hbond_scalars_2, catalog_hbond_scalars_3 | 475116 x4 |  |
| hbond_nearest_dir | hbond | atom | atom_row_sidecar: row_design_field_hbond_nearest_dir.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_hbond_nearest_dir.npy | 475116 x3 |  |
| hbond_flags | hbond | atom | row columns: catalog_hbond_flags_0, catalog_hbond_flags_1, catalog_hbond_flags_2 | 635346 x3 |  | row columns: catalog_hbond_flags_0, catalog_hbond_flags_1, catalog_hbond_flags_2 | 475116 x3 |  |
| dssp_backbone | dssp | atom | row columns (5): catalog_dssp_backbone_0 ... catalog_dssp_backbone_4 | 635346 x5 |  | row columns (5): catalog_dssp_backbone_0 ... catalog_dssp_backbone_4 | 475116 x5 |  |
| dssp_ss8 | dssp | atom | row columns (8): catalog_dssp_ss8_0 ... catalog_dssp_ss8_7 | 635346 x8 |  | row columns (8): catalog_dssp_ss8_0 ... catalog_dssp_ss8_7 | 475116 x8 |  |
| dssp_hbond_energy | dssp | atom | row columns: catalog_dssp_hbond_energy_0, catalog_dssp_hbond_energy_1, catalog_dssp_hbond_energy_2, catalog_dssp_hbond_energy_3 | 635346 x4 |  | row columns: catalog_dssp_hbond_energy_0, catalog_dssp_hbond_energy_1, catalog_dssp_hbond_energy_2, catalog_dssp_hbond_energy_3 | 475116 x4 |  |
| dssp_chi | dssp | atom | row columns (12): catalog_dssp_chi_0 ... catalog_dssp_chi_11 | 635346 x12 |  | row columns (12): catalog_dssp_chi_0 ... catalog_dssp_chi_11 | 475116 x12 |  |
| atom_sasa | sasa | atom | row columns: catalog_atom_sasa | 635346 |  | row columns: catalog_atom_sasa | 475116 |  |
| sasa_normal | sasa | atom | atom_row_sidecar: row_design_field_sasa_normal.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_sasa_normal.npy | 475116 x3 |  |
| water_efield | water_field | atom | atom_row_sidecar: row_design_field_water_efield.npy | 635346 x3 |  | atom_row_sidecar_absent | 0 x3 |  |
| water_efield_first | water_field | atom | atom_row_sidecar: row_design_field_water_efield_first.npy | 635346 x3 |  | atom_row_sidecar_absent | 0 x3 |  |
| water_efg | water_field | atom | atom_row_sidecar: row_design_field_water_efg.npy | 635346 x5 |  | atom_row_sidecar_absent | 0 x5 |  |
| water_efg_first | water_field | atom | atom_row_sidecar: row_design_field_water_efg_first.npy | 635346 x5 |  | atom_row_sidecar_absent | 0 x5 |  |
| water_shell_counts | water_field | atom | row columns: catalog_water_shell_counts_0, catalog_water_shell_counts_1 | 635346 x2 |  | row columns: catalog_water_shell_counts_0, catalog_water_shell_counts_1 | 0 x2 |  |
| hydration_shell | hydration | atom | row columns: catalog_hydration_shell_0, catalog_hydration_shell_1, catalog_hydration_shell_2, catalog_hydration_shell_3 | [635346, 635346, 501941, 635346] |  | row columns: catalog_hydration_shell_0, catalog_hydration_shell_1, catalog_hydration_shell_2, catalog_hydration_shell_3 | 0 x4 |  |
| water_polarization | water_polarization | atom | row columns (10): catalog_water_polarization_0 ... catalog_water_polarization_9 | 635346 x10 |  | row columns (10): catalog_water_polarization_0 ... catalog_water_polarization_9 | 0 x10 |  |
| eeq_charges | eeq | atom | row columns: catalog_eeq_charges | 635346 |  | row columns: catalog_eeq_charges | 475116 |  |
| eeq_cn | eeq | atom | row columns: catalog_eeq_cn | 635346 |  | row columns: catalog_eeq_cn | 475116 |  |
| gromacs_energy | gromacs | protein | native_axis_sidecar: row_design_field_gromacs_energy_native.npy | 43 comps; min=0 max=751; common: 751 x40, 0 x3; exact in manifest/CSV | row_design_field_gromacs_energy_reduction.npy (43 comps; min=0 max=635346; common: 635346 x40, 0 x3; exact in manifest/CSV) | native_axis_sidecar_absent | 0 x43 |  |
| bonded_energy | bonded | atom | row columns (7): catalog_bonded_energy_0 ... catalog_bonded_energy_6 | 635346 x7 |  | row columns (7): catalog_bonded_energy_0 ... catalog_bonded_energy_6 | 0 x7 |  |
| mopac_charges | mopac_core | atom | row columns: catalog_mopac_charges | 635346 |  | row columns: catalog_mopac_charges | 475116 |  |
| mopac_scalars | mopac_core | atom | row columns: catalog_mopac_scalars_0, catalog_mopac_scalars_1, catalog_mopac_scalars_2, catalog_mopac_scalars_3 | 635346 x4 |  | row columns: catalog_mopac_scalars_0, catalog_mopac_scalars_1, catalog_mopac_scalars_2, catalog_mopac_scalars_3 | 475116 x4 |  |
| mopac_bond_orders | mopac_core | bond | native_axis_sidecar: row_design_field_mopac_bond_orders_native.npy | 1072765 x3 | row_design_field_mopac_bond_orders_reduction.npy (635346 x4) | native_axis_sidecar: row_design_field_mopac_bond_orders_native.npy | 765120 x3 | row_design_field_mopac_bond_orders_reduction.npy (475116 x4) |
| mopac_bond_neighbors | mopac_core | mopac_bond_neighbor_pair | native_axis_sidecar: row_design_field_mopac_bond_neighbors_native.npy | 2145530 x4 | row_design_field_mopac_bond_neighbors_reduction.npy (635346 x13) | native_axis_sidecar: row_design_field_mopac_bond_neighbors_native.npy | 1530240 x4 | row_design_field_mopac_bond_neighbors_reduction.npy (475116 x13) |
| mopac_global | mopac_core | protein | native_axis_sidecar: row_design_field_mopac_global_native.npy | 751 x4 | row_design_field_mopac_global_reduction.npy (635346 x4) | native_axis_sidecar: row_design_field_mopac_global_native.npy | 720 x4 | row_design_field_mopac_global_reduction.npy (475116 x4) |
| mopac_atom_populations | mopac_core | atom | row columns (12): catalog_mopac_atom_populations_0 ... catalog_mopac_atom_populations_11 | [635346, 635346, 635346, 331942, 5257, 0, 0, 0, 0, 0, 635346, 635346] |  | row columns (12): catalog_mopac_atom_populations_0 ... catalog_mopac_atom_populations_11 | [475116, 475116, 475116, 233480, 1796, 0, 0, 0, 0, 0, 475116, 475116] |  |
| mopac_atomic_orbital_populations | mopac_core | atom | row columns (9): catalog_mopac_atomic_orbital_populations_0 ... catalog_mopac_atomic_orbital_populations_8 | [635346, 331942, 331942, 331942, 5257, 5257, 5257, 5257, 5257] |  | row columns (9): catalog_mopac_atomic_orbital_populations_0 ... catalog_mopac_atomic_orbital_populations_8 | [475116, 233480, 233480, 233480, 1796, 1796, 1796, 1796, 1796] |  |
| mopac_bond_valencies | mopac_core | atom | row columns: catalog_mopac_bond_valencies | 635346 |  | row columns: catalog_mopac_bond_valencies | 475116 |  |
| mopac_bond_orders_unique | mopac_core | mopac_unique_pair | native_axis_sidecar: row_design_field_mopac_bond_orders_unique_native.npy | [2948278, 2948278, 2948278, 2948278, 2948278, 2948278, 2948278, 647362] | row_design_field_mopac_bond_orders_unique_reduction.npy (635346 x9) | native_axis_sidecar: row_design_field_mopac_bond_orders_unique_native.npy | [2076213, 2076213, 2076213, 2076213, 2076213, 2076213, 2076213, 478093] | row_design_field_mopac_bond_orders_unique_reduction.npy (475116 x9) |
| mopac_topology_bond_orders_full | mopac_core | bond | native_axis_sidecar: row_design_field_mopac_topology_bond_orders_full_native.npy | 647362 x8 | row_design_field_mopac_topology_bond_orders_full_reduction.npy (635346 x9) | native_axis_sidecar: row_design_field_mopac_topology_bond_orders_full_native.npy | 478093 x8 | row_design_field_mopac_topology_bond_orders_full_reduction.npy ([475116, 475112, 475112, 475112, 475112, 475112, 475112, 475112, 475112]) |
| mopac_coulomb_efg | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_efg.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_mopac_coulomb_efg.npy | 475116 x9 |  |
| mopac_coulomb_E | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_E.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_mopac_coulomb_E.npy | 475116 x3 |  |
| mopac_coulomb_E_backbone | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_E_backbone.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_mopac_coulomb_E_backbone.npy | 475116 x3 |  |
| mopac_coulomb_E_sidechain | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_E_sidechain.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_mopac_coulomb_E_sidechain.npy | 475116 x3 |  |
| mopac_coulomb_E_aromatic | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_E_aromatic.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_mopac_coulomb_E_aromatic.npy | 475116 x3 |  |
| mopac_coulomb_efg_backbone | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_efg_backbone.npy | 635346 x5 |  | atom_row_sidecar: row_design_field_mopac_coulomb_efg_backbone.npy | 475116 x5 |  |
| mopac_coulomb_efg_sidechain | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_efg_sidechain.npy | 635346 x5 |  | atom_row_sidecar: row_design_field_mopac_coulomb_efg_sidechain.npy | 475116 x5 |  |
| mopac_coulomb_efg_aromatic | mopac_coulomb | atom | atom_row_sidecar: row_design_field_mopac_coulomb_efg_aromatic.npy | 635346 x5 |  | atom_row_sidecar: row_design_field_mopac_coulomb_efg_aromatic.npy | 475116 x5 |  |
| mopac_coulomb_scalars | mopac_coulomb | atom | row columns: catalog_mopac_coulomb_scalars_0, catalog_mopac_coulomb_scalars_1, catalog_mopac_coulomb_scalars_2, catalog_mopac_coulomb_scalars_3 | 635346 x4 |  | row columns: catalog_mopac_coulomb_scalars_0, catalog_mopac_coulomb_scalars_1, catalog_mopac_coulomb_scalars_2, catalog_mopac_coulomb_scalars_3 | 475116 x4 |  |
| apbs_E | apbs | atom | atom_row_sidecar: row_design_field_apbs_E.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_apbs_E.npy | 475116 x3 |  |
| apbs_efg | apbs | atom | atom_row_sidecar: row_design_field_apbs_efg.npy | 635346 x5 |  | atom_row_sidecar: row_design_field_apbs_efg.npy | 475116 x5 |  |
| orca_total | orca | atom | atom_row_sidecar: row_design_field_orca_total.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_orca_total.npy | 475116 x9 |  |
| orca_diamagnetic | orca | atom | atom_row_sidecar: row_design_field_orca_diamagnetic.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_orca_diamagnetic.npy | 475116 x9 |  |
| orca_paramagnetic | orca | atom | atom_row_sidecar: row_design_field_orca_paramagnetic.npy | 635346 x9 |  | atom_row_sidecar: row_design_field_orca_paramagnetic.npy | 475116 x9 |  |
| atoms_category_info | identity | atom | structured_native_sidecar: row_design_field_atoms_category_info_native.npy + index row_design_field_atoms_category_info_index.csv | 846 |  | structured_native_sidecar: row_design_field_atoms_category_info_native.npy + index row_design_field_atoms_category_info_index.csv | 475116 |  |
| aimnet2_charges | aimnet2 | atom | row columns: catalog_aimnet2_charges | 635346 |  | row columns: catalog_aimnet2_charges | 475116 |  |
| aimnet2_aim | aimnet2 | atom | atom_row_sidecar: row_design_field_aimnet2_aim.npy | 635346 x256 |  | atom_row_sidecar: row_design_field_aimnet2_aim.npy | 475116 x256 |  |
| aimnet2_efg | aimnet2 | atom | atom_row_sidecar: row_design_field_aimnet2_efg.npy | 635346 x5 |  | atom_row_sidecar: row_design_field_aimnet2_efg.npy | 475116 x5 |  |
| aimnet2_efg_aromatic | aimnet2 | atom | atom_row_sidecar: row_design_field_aimnet2_efg_aromatic.npy | 635346 x5 |  | atom_row_sidecar: row_design_field_aimnet2_efg_aromatic.npy | 475116 x5 |  |
| aimnet2_efg_backbone | aimnet2 | atom | atom_row_sidecar: row_design_field_aimnet2_efg_backbone.npy | 635346 x5 |  | atom_row_sidecar: row_design_field_aimnet2_efg_backbone.npy | 475116 x5 |  |
| aimnet2_charge_response_gradient | aimnet2 | atom | atom_row_sidecar: row_design_field_aimnet2_charge_response_gradient.npy | 635346 x3 |  | atom_row_sidecar: row_design_field_aimnet2_charge_response_gradient.npy | 475116 x3 |  |
| aimnet2_charge_response_gradient_scalar | aimnet2 | atom | row columns: catalog_aimnet2_charge_response_gradient_scalar | 635346 |  | row columns: catalog_aimnet2_charge_response_gradient_scalar | 475116 |  |
| pyramidalization | planar_geometry | atom | row columns: catalog_pyramidalization | 635346 |  | row columns: catalog_pyramidalization | 475116 |  |
| omega_actual | planar_geometry | residue | native_axis_sidecar: row_design_field_omega_actual_native.npy | 39803 | row_design_field_omega_actual_reduction.npy (620326) | native_axis_sidecar: row_design_field_omega_actual_native.npy | 29224 | row_design_field_omega_actual_reduction.npy (463775) |
| omega_deviation | planar_geometry | residue | native_axis_sidecar: row_design_field_omega_deviation_native.npy | 39803 | row_design_field_omega_deviation_reduction.npy (620326) | native_axis_sidecar: row_design_field_omega_deviation_native.npy | 29224 | row_design_field_omega_deviation_reduction.npy (463775) |
| aromatic_chi2 | planar_geometry | aromatic_ring | native_axis_sidecar: row_design_field_aromatic_chi2_native.npy | 11265 | row_design_field_aromatic_chi2_reduction.npy (51819) | native_axis_sidecar: row_design_field_aromatic_chi2_native.npy | 3015 | row_design_field_aromatic_chi2_reduction.npy (15034) |
| pucker_Q | planar_geometry | saturated_ring | native_axis_sidecar: row_design_field_pucker_Q_native.npy | 751 | row_design_field_pucker_Q_reduction.npy (3755) | native_axis_sidecar: row_design_field_pucker_Q_native.npy | 965 | row_design_field_pucker_Q_reduction.npy (4825) |
| pucker_theta | planar_geometry | saturated_ring | native_axis_sidecar: row_design_field_pucker_theta_native.npy | 751 | row_design_field_pucker_theta_reduction.npy (3755) | native_axis_sidecar: row_design_field_pucker_theta_native.npy | 965 | row_design_field_pucker_theta_reduction.npy (4825) |
| omega_is_xpro | planar_geometry | residue | native_axis_sidecar: row_design_field_omega_is_xpro_native.npy | 40554 | row_design_field_omega_is_xpro_reduction.npy (635346) | native_axis_sidecar: row_design_field_omega_is_xpro_native.npy | 29944 | row_design_field_omega_is_xpro_reduction.npy (475116) |
| residues | topology | residue | structured_native_sidecar: row_design_field_residues_native.npy + index row_design_field_residues_index.csv | 54 |  | structured_native_sidecar: row_design_field_residues_native.npy + index row_design_field_residues_index.csv | 29944 |  |
| bonds | topology | bond | structured_native_sidecar: row_design_field_bonds_native.npy + index row_design_field_bonds_index.csv | 862 |  | structured_native_sidecar: row_design_field_bonds_native.npy + index row_design_field_bonds_index.csv | 478093 |  |
| rings | topology | ring | structured_native_sidecar: row_design_field_rings_native.npy + index row_design_field_rings_index.csv | 16 |  | structured_native_sidecar: row_design_field_rings_native.npy + index row_design_field_rings_index.csv | 3980 |  |
| ring_membership | topology | ring_membership | structured_native_sidecar: row_design_field_ring_membership_native.npy + index row_design_field_ring_membership_index.csv | 96 |  | structured_native_sidecar: row_design_field_ring_membership_native.npy + index row_design_field_ring_membership_index.csv | 22862 |  |

## Dataset-Specific Absence

720 has explicit zero-count manifest entries for fields absent from the static producer arrays:

- `coulomb_efg`: atom_row_sidecar_absent, counts 0 x9
- `coulomb_E`: atom_row_sidecar_absent, counts 0 x3
- `coulomb_E_backbone`: atom_row_sidecar_absent, counts 0 x3
- `coulomb_E_sidechain`: atom_row_sidecar_absent, counts 0 x3
- `coulomb_E_aromatic`: atom_row_sidecar_absent, counts 0 x3
- `coulomb_efg_backbone`: atom_row_sidecar_absent, counts 0 x5
- `coulomb_efg_sidechain`: atom_row_sidecar_absent, counts 0 x5
- `coulomb_efg_aromatic`: atom_row_sidecar_absent, counts 0 x5
- `water_efield`: atom_row_sidecar_absent, counts 0 x3
- `water_efield_first`: atom_row_sidecar_absent, counts 0 x3
- `water_efg`: atom_row_sidecar_absent, counts 0 x5
- `water_efg_first`: atom_row_sidecar_absent, counts 0 x5
- `gromacs_energy`: native_axis_sidecar_absent, counts 0 x43

1P9J has no scoped field represented as absent.

## Reduction Policies

- `ring_contribution_pair`: per-(atom,ring): pair count, sum over rings, nearest-ring value by ring_contributions distance, and ring-type sum bins
- `aromatic_ring`: primary ring-membership policy: first topology membership matching the catalog native axis
- `protein`: broadcast protein-axis values to every emitted atom row
- `bond`: generic incident endpoint sum for bond-like native axis; review downstream semantics before modelling
- `mopac_bond_neighbor_pair`: mopac_bond_neighbors incident aggregate: count, sum, mean, max by atom endpoint
- `mopac_unique_pair`: generic incident endpoint sum for bond-like native axis; review downstream semantics before modelling
- `residue`: broadcast residue-axis values by atom.residue_index
- `saturated_ring`: primary ring-membership policy: first topology membership matching the catalog native axis

Generic reductions flagged for downstream review before modelling: `mopac_bond_orders`, `mopac_bond_orders_unique`, `mopac_topology_bond_orders_full`.

