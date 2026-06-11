# CODEX Spine Wire-up Report

## Probe Result

| Dataset | Manifest | Fields | Failed | Result | Flat coverage compared |
| --- | --- | ---: | ---: | --- | --- |
| 720 static pose P84337 | `/tmp/codex-spine-probe-static/spine_reachability_static.json` | 122 | 66 | FAIL | false |
| 1P9J trajectory | `/tmp/codex-spine-probe-1p9j/spine_reachability_trajectory.json` | 122 | 81 | FAIL | false |

The acceptance gate did not pass: the static probe still has 66 failed scoped fields and the 1P9J trajectory probe has 81 failed scoped fields. The failures are explicit resident-source absences or shape/topology absence reasons; no field is treated as present by default.

Because the gate did not pass, no cleanup deletions were performed. `RowDesignCatalogCoverage`, the flat emitter path, the 19GB sidecar dump, and the three old-sink sidecars remain untouched.

## Implementation Notes

- The FieldKind spine is driven from `ScopedProducerCatalog()` and uses `FieldSpec.cols` for fixed component width.
- Providers are explicit per field: static producer array, dense H5 time-series, sparse DFT by original frame, typed topology, or dataset absent.
- Static producer reads use resident native rows, not atom-major `rowFor()` for non-atom fields.
- Structured topology fields are sourced through `QtProtein`/`QtTopology`; numeric `value()` is not used for structured rows.
- Dense 1P9J mappings include the typed DSSP SS8/H-bond/chi broadcasts and `water_polarization`; `dssp_backbone` remains absent because this H5 lacks the DSSP per-residue SASA component required for its full five-column representation.

## Scoped Field Matrix

| Field | Group | Axis | Static provider | Static extent | Static result | 1P9J provider | 1P9J extent | 1P9J result |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `pos` | identity | atom | static-producer-array | 412r x 1f x 3c | PASS (1236 reads) | dense-h5-time-series | 846r x 1501f x 3c | PASS (3809538 reads) |
| `element` | identity | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | typed-topology | 846r x 1f x 1c | PASS (846 reads) |
| `residue_index` | identity | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | typed-topology | 846r x 1f x 1c | PASS (846 reads) |
| `residue_type` | identity | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | typed-topology | 846r x 1f x 1c | PASS (846 reads) |
| `ring_contributions` | identity | ring_contribution_pair | dataset-absent | 0r x 0f x 40c | FAIL: malformed-shape: calibration/features/FinalExtraction/P84337/ring_contributions.npy has 58 columns; catalog ring_contributions expects 40 | dataset-absent | 0r x 0f x 40c | FAIL: not-produced-in-dataset |
| `ring_direction_to_center` | identity | ring_contribution_pair | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `ring_geometry` | identity | aromatic_ring | static-producer-array | 1r x 1f x 10c | PASS (10 reads) | dataset-absent | 0r x 0f x 10c | FAIL: not-produced-in-dataset |
| `enrichment_role` | enrichment | atom | dataset-absent | 0r x 0f x 0c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `enrichment_hybridisation` | enrichment | atom | dataset-absent | 0r x 0f x 0c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `enrichment_flags` | enrichment | atom | dataset-absent | 0r x 0f x 8c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 8c | FAIL: not-produced-in-dataset |
| `ff_partial_charge` | charge_assignment | atom | typed-topology | 412r x 0f x 1c | FAIL: topology-mismatch | typed-topology | 846r x 1f x 1c | PASS (846 reads) |
| `ff_pb_radius` | charge_assignment | atom | dataset-absent | 0r x 0f x 0c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `bs_shielding` | biot_savart | atom | static-producer-array | 412r x 1f x 9c | PASS (3708 reads) | dense-h5-time-series | 846r x 1501f x 9c | PASS (11428614 reads) |
| `bs_per_type_T0` | biot_savart | atom | static-producer-array | 412r x 1f x 8c | PASS (3296 reads) | dataset-absent | 0r x 0f x 8c | FAIL: not-produced-in-dataset |
| `bs_per_type_T1` | biot_savart | atom | dataset-absent | 0r x 0f x 24c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 24c | FAIL: not-produced-in-dataset |
| `bs_per_type_T2` | biot_savart | atom | static-producer-array | 412r x 1f x 40c | PASS (16480 reads) | dataset-absent | 0r x 0f x 40c | FAIL: not-produced-in-dataset |
| `bs_total_B` | biot_savart | atom | static-producer-array | 412r x 1f x 3c | PASS (1236 reads) | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `bs_ring_B_field` | biot_savart | ring_contribution_pair | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `bs_ring_B_cylindrical` | biot_savart | ring_contribution_pair | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `bs_ring_counts` | biot_savart | atom | static-producer-array | 412r x 1f x 4c | PASS (1648 reads) | dataset-absent | 0r x 0f x 4c | FAIL: not-produced-in-dataset |
| `hm_shielding` | haigh_mallion | atom | static-producer-array | 412r x 1f x 9c | PASS (3708 reads) | dense-h5-time-series | 846r x 1501f x 9c | PASS (11428614 reads) |
| `hm_per_type_T0` | haigh_mallion | atom | static-producer-array | 412r x 1f x 8c | PASS (3296 reads) | dataset-absent | 0r x 0f x 8c | FAIL: not-produced-in-dataset |
| `hm_per_type_T1` | haigh_mallion | atom | dataset-absent | 0r x 0f x 24c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 24c | FAIL: not-produced-in-dataset |
| `hm_per_type_T2` | haigh_mallion | atom | static-producer-array | 412r x 1f x 40c | PASS (16480 reads) | dataset-absent | 0r x 0f x 40c | FAIL: not-produced-in-dataset |
| `hm_ring_B_field` | haigh_mallion | ring_contribution_pair | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `pq_per_type_T0` | pi_quadrupole | atom | static-producer-array | 412r x 1f x 8c | PASS (3296 reads) | dataset-absent | 0r x 0f x 8c | FAIL: not-produced-in-dataset |
| `piquad_quad_scalar` | pi_quadrupole | ring_contribution_pair | dataset-absent | 0r x 0f x 0c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `disp_per_type_T0` | dispersion | atom | static-producer-array | 412r x 1f x 8c | PASS (3296 reads) | dataset-absent | 0r x 0f x 8c | FAIL: not-produced-in-dataset |
| `mc_peptide_co_fixed` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_peptide_co_bo` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_peptide_co_rhombic` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_peptide_cn_fixed` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_peptide_cn_bo` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_backbone_other_fixed` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_backbone_other_bo` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_sidechain_co_fixed` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_sidechain_co_bo` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_sidechain_other_fixed` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_sidechain_other_bo` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_disulfide_fixed` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_disulfide_bo` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_aromatic_zeroed_fixed` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_aromatic_zeroed_bo` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_nearfield_counts` | mcconnell | atom | dataset-absent | 0r x 0f x 2c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 2c | FAIL: not-produced-in-dataset |
| `mc_nearest_co_dir` | mcconnell | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `mc_nearest_co_midpoint` | mcconnell | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `mc_nearest_co_T2` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mc_nearest_cn_T2` | mcconnell | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `coulomb_efg` | coulomb | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `coulomb_E` | coulomb | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `coulomb_E_backbone` | coulomb | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `coulomb_E_sidechain` | coulomb | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `coulomb_E_aromatic` | coulomb | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `coulomb_efg_backbone` | coulomb | atom | dataset-absent | 0r x 0f x 5c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `coulomb_efg_sidechain` | coulomb | atom | dataset-absent | 0r x 0f x 5c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `coulomb_efg_aromatic` | coulomb | atom | dataset-absent | 0r x 0f x 5c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `coulomb_scalars` | coulomb | atom | dataset-absent | 0r x 0f x 4c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 4c | FAIL: not-produced-in-dataset |
| `coulomb_aromatic_E_proj` | coulomb | atom | dataset-absent | 0r x 0f x 0c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `coulomb_aromatic_n_src` | coulomb | atom | dataset-absent | 0r x 0f x 0c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `hbond_scalars` | hbond | atom | static-producer-array | 412r x 1f x 4c | PASS (1648 reads) | dataset-absent | 0r x 0f x 4c | FAIL: not-produced-in-dataset |
| `hbond_nearest_dir` | hbond | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `hbond_flags` | hbond | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `dssp_backbone` | dssp | atom | static-producer-array | 412r x 1f x 5c | PASS (2060 reads) | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `dssp_ss8` | dssp | atom | static-producer-array | 412r x 1f x 8c | PASS (3296 reads) | dense-h5-time-series | 846r x 1501f x 8c | PASS (10158768 reads) |
| `dssp_hbond_energy` | dssp | atom | static-producer-array | 412r x 1f x 4c | PASS (1648 reads) | dense-h5-time-series | 846r x 1501f x 4c | PASS (5079384 reads) |
| `dssp_chi` | dssp | atom | static-producer-array | 412r x 1f x 12c | PASS (4944 reads) | dense-h5-time-series | 846r x 1501f x 12c | PASS (15238152 reads) |
| `atom_sasa` | sasa | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | dense-h5-time-series | 846r x 1501f x 1c | PASS (1269846 reads) |
| `sasa_normal` | sasa | atom | static-producer-array | 412r x 1f x 3c | PASS (1236 reads) | dense-h5-time-series | 846r x 1501f x 3c | PASS (3809538 reads) |
| `water_efield` | water_field | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 3c | PASS (3809538 reads) |
| `water_efield_first` | water_field | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 3c | PASS (3809538 reads) |
| `water_efg` | water_field | atom | dataset-absent | 0r x 0f x 5c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 5c | PASS (6349230 reads) |
| `water_efg_first` | water_field | atom | dataset-absent | 0r x 0f x 5c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 5c | PASS (6349230 reads) |
| `water_shell_counts` | water_field | atom | dataset-absent | 0r x 0f x 2c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 2c | PASS (2539692 reads) |
| `hydration_shell` | hydration | atom | dataset-absent | 0r x 0f x 4c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 4c | PASS (5079384 reads) |
| `water_polarization` | water_polarization | atom | dataset-absent | 0r x 0f x 10c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 10c | PASS (12698460 reads) |
| `eeq_charges` | eeq | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `eeq_cn` | eeq | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `gromacs_energy` | gromacs | protein | dataset-absent | 0r x 0f x 43c | FAIL: producer-array-absent | dense-h5-time-series | 1r x 1501f x 43c | PASS (64543 reads) |
| `bonded_energy` | bonded | atom | dataset-absent | 0r x 0f x 7c | FAIL: producer-array-absent | dense-h5-time-series | 846r x 1501f x 7c | PASS (8888922 reads) |
| `mopac_charges` | mopac_core | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `mopac_scalars` | mopac_core | atom | static-producer-array | 412r x 1f x 4c | PASS (1648 reads) | dataset-absent | 0r x 0f x 4c | FAIL: not-produced-in-dataset |
| `mopac_bond_orders` | mopac_core | bond | static-producer-array | 333r x 1f x 3c | PASS (999 reads) | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `mopac_bond_neighbors` | mopac_core | mopac_bond_neighbor_pair | dataset-absent | 0r x 0f x 4c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 4c | FAIL: not-produced-in-dataset |
| `mopac_global` | mopac_core | protein | static-producer-array | 1r x 1f x 4c | PASS (4 reads) | dataset-absent | 0r x 0f x 4c | FAIL: not-produced-in-dataset |
| `mopac_atom_populations` | mopac_core | atom | dataset-absent | 0r x 0f x 12c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 12c | FAIL: not-produced-in-dataset |
| `mopac_atomic_orbital_populations` | mopac_core | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mopac_bond_valencies` | mopac_core | atom | dataset-absent | 0r x 0f x 0c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `mopac_bond_orders_unique` | mopac_core | mopac_unique_pair | dataset-absent | 0r x 0f x 8c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 8c | FAIL: not-produced-in-dataset |
| `mopac_topology_bond_orders_full` | mopac_core | bond | dataset-absent | 0r x 0f x 8c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 8c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_efg` | mopac_coulomb | atom | dataset-absent | 0r x 0f x 9c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 9c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_E` | mopac_coulomb | atom | static-producer-array | 412r x 1f x 3c | PASS (1236 reads) | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_E_backbone` | mopac_coulomb | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_E_sidechain` | mopac_coulomb | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_E_aromatic` | mopac_coulomb | atom | dataset-absent | 0r x 0f x 3c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 3c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_efg_backbone` | mopac_coulomb | atom | static-producer-array | 412r x 1f x 5c | PASS (2060 reads) | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_efg_sidechain` | mopac_coulomb | atom | dataset-absent | 0r x 0f x 5c | FAIL: producer-array-absent | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_efg_aromatic` | mopac_coulomb | atom | static-producer-array | 412r x 1f x 5c | PASS (2060 reads) | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `mopac_coulomb_scalars` | mopac_coulomb | atom | static-producer-array | 412r x 1f x 4c | PASS (1648 reads) | dataset-absent | 0r x 0f x 4c | FAIL: not-produced-in-dataset |
| `apbs_E` | apbs | atom | static-producer-array | 412r x 1f x 3c | PASS (1236 reads) | dense-h5-time-series | 846r x 1501f x 3c | PASS (3809538 reads) |
| `apbs_efg` | apbs | atom | static-producer-array | 412r x 1f x 5c | PASS (2060 reads) | dense-h5-time-series | 846r x 1501f x 5c | PASS (6349230 reads) |
| `orca_total` | orca | atom | static-producer-array | 412r x 1f x 9c | PASS (3708 reads) | sparse-dft-by-original | 846r x 751f x 9c | PASS (5718114 reads) |
| `orca_diamagnetic` | orca | atom | static-producer-array | 412r x 1f x 9c | PASS (3708 reads) | sparse-dft-by-original | 846r x 751f x 9c | PASS (5718114 reads) |
| `orca_paramagnetic` | orca | atom | static-producer-array | 412r x 1f x 9c | PASS (3708 reads) | sparse-dft-by-original | 846r x 751f x 9c | PASS (5718114 reads) |
| `atoms_category_info` | identity | atom | typed-topology | 412r x 1f x structuredc | PASS (412 reads) | typed-topology | 846r x 1f x structuredc | PASS (846 reads) |
| `aimnet2_charges` | aimnet2 | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | dense-h5-time-series | 846r x 1501f x 1c | PASS (1269846 reads) |
| `aimnet2_aim` | aimnet2 | atom | static-producer-array | 412r x 1f x 256c | PASS (105472 reads) | dense-h5-time-series | 846r x 1501f x 256c | PASS (325080576 reads) |
| `aimnet2_efg` | aimnet2 | atom | static-producer-array | 412r x 1f x 5c | PASS (2060 reads) | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `aimnet2_efg_aromatic` | aimnet2 | atom | static-producer-array | 412r x 1f x 5c | PASS (2060 reads) | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `aimnet2_efg_backbone` | aimnet2 | atom | static-producer-array | 412r x 1f x 5c | PASS (2060 reads) | dataset-absent | 0r x 0f x 5c | FAIL: not-produced-in-dataset |
| `aimnet2_charge_response_gradient` | aimnet2 | atom | static-producer-array | 412r x 1f x 3c | PASS (1236 reads) | dense-h5-time-series | 846r x 1501f x 3c | PASS (3809538 reads) |
| `aimnet2_charge_response_gradient_scalar` | aimnet2 | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | dense-h5-time-series | 846r x 1501f x 1c | PASS (1269846 reads) |
| `pyramidalization` | planar_geometry | atom | static-producer-array | 412r x 1f x 1c | PASS (412 reads) | dataset-absent | 0r x 0f x 0c | FAIL: not-produced-in-dataset |
| `omega_actual` | planar_geometry | residue | static-producer-array | 24r x 1f x 1c | PASS (24 reads) | dense-h5-time-series | 54r x 1501f x 1c | PASS (81054 reads) |
| `omega_deviation` | planar_geometry | residue | static-producer-array | 24r x 1f x 1c | PASS (24 reads) | dense-h5-time-series | 54r x 1501f x 1c | PASS (81054 reads) |
| `aromatic_chi2` | planar_geometry | aromatic_ring | static-producer-array | 1r x 1f x 1c | PASS (1 reads) | dense-h5-time-series | 15r x 1501f x 1c | PASS (22515 reads) |
| `pucker_Q` | planar_geometry | saturated_ring | static-producer-array | 0r x 1f x 1c | PASS (0 reads) | dense-h5-time-series | 1r x 1501f x 1c | PASS (1501 reads) |
| `pucker_theta` | planar_geometry | saturated_ring | static-producer-array | 0r x 1f x 1c | PASS (0 reads) | dense-h5-time-series | 1r x 1501f x 1c | PASS (1501 reads) |
| `omega_is_xpro` | planar_geometry | residue | static-producer-array | 24r x 1f x 1c | PASS (24 reads) | dense-h5-time-series | 54r x 1f x 1c | PASS (54 reads) |
| `residues` | topology | residue | typed-topology | 24r x 1f x structuredc | PASS (24 reads) | typed-topology | 54r x 1f x structuredc | PASS (54 reads) |
| `bonds` | topology | bond | typed-topology | 412r x 1f x structuredc | PASS (412 reads) | typed-topology | 862r x 1f x structuredc | PASS (862 reads) |
| `rings` | topology | ring | typed-topology | 1r x 1f x structuredc | PASS (1 reads) | typed-topology | 16r x 1f x structuredc | PASS (16 reads) |
| `ring_membership` | topology | ring_membership | typed-topology | 5r x 1f x structuredc | PASS (5 reads) | typed-topology | 96r x 1f x structuredc | PASS (96 reads) |

## Deleted Artifacts

None. The acceptance gate failed, so the flat path and sidecars were preserved.
