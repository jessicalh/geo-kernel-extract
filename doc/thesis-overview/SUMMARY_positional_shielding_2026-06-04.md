# Positional Geometric-Kernel Shielding Summary, 2026-06-04

## Path aliases

- `LEARN=/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn`
- `S1=/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage1-topology-workup-20260514`
- `AA=/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage1-atom-audit`
- `OPUS=/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage1-topology-workup-20260514/derived-opus`
- `STAGE2=/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage2`

## Scope / Language Discipline

- The calculators output geometric kernels; ORCA provides shielding targets. All `R^2` values below are kernel-vs-ORCA shielding correlations or diagnostic fits, not pointwise agreement to literature shifts and not calibrated shielding laws. Source: `LEARN/ARGUMENT.md:5-9`, `S1/docs/STAGE1_COMPENDIUM.md:148-160`.
- The current clean topology workup uses `calibration/features/Stage1BMRB_20260513_topology`, 720 OK protein directories, 425,599 matched atom rows, 788 feature columns, and 23 target columns. Source: `S1/README.md:3-24`, `S1/docs/STAGE1_COMPENDIUM.md:259-266`.
- The main positional/local-geometry cluster is `S1/derived/local_geometry_cloud_probe_20260518_002`, plus the dihedral, response-shape, directionality, torsion-axis, and `derived-opus` v2 chart outputs. Source: `S1/docs/STAGE1_COMPENDIUM.md:136-140`, `S1/docs/SCRIPT_PROVENANCE.md:74-93`, `OPUS/README.md:55-64`.
- Older atom-audit results are baseline/framing: useful for naming the questions, but the topology-audited `20260513` export supersedes its regex/display-label debt. Source: `AA/STAGE1_720_BASELINE_2026-05-10.md:9-11`, `AA/STAGE1_720_BASELINE_2026-05-10.md:444-459`, `S1/README.md:16-24`.

## 1. Direct Questions

### Core Stage 1 Atom/Kernel Questions

- Do named atom/context strata, not just elements, carry distinct kernel-vs-ORCA shielding signal? Source: `AA/PHYSICS_OVERLAP_DELIVERABLE.md:19-35`, `AA/PHYSICS_OVERLAP_DELIVERABLE.md:173-184`.
- For each supported named atom stratum, which geometric-kernel blocks correlate with `total_iso`, `dia_iso`, `para_iso`, `total_t2_norm`, `dia_t2_norm`, and `para_t2_norm` ORCA targets? Source: `S1/R/01_block_ridge.R:23-33`, `S1/R/01_block_ridge.R:64-100`.
- Are high diagnostic `R^2` values support/capacity artifacts, or do they survive protein-level blocking and resampling? Source: `AA/STAGE1_720_FINDINGS_2026-05-09.md:45-71`, `S1/R/07_mechanism_blocked_cv.R:3-9`, `S1/R/21_blocked_stability.R:3-8`.
- Which source blocks add independent signal after other blocks are present, and which are redundant/proxy blocks? Source: `S1/R/02_mechanism_independence.R` as described by `S1/docs/SCRIPT_PROVENANCE.md:15-17`; pairwise definitions from `S1/R/12_atom_pairwise_independence.R` are summarized in `S1/docs/ATOM_SIGNAL_QUESTION_TABLES.md:18-35`.
- How many overlapping source dimensions are active in each atom/target/normalization stratum, without collapsing the result to a single "3-dimensional" story? Source: `S1/R/03_dimension_basis.R:60-73`, `S1/R/03_dimension_basis.R:102-170`, `S1/docs/ATOM_SIGNAL_QUESTION_TABLES.md:30-35`.

### Positional / Geometry Questions

- How do backbone phi/psi/omega/chi/SASA/DSSP descriptors correlate with raw ORCA shielding summaries in regular proteins, separate from mutation deltas? Source: `S1/derived/dihedral_backbone_analysis_20260518_002/DIHEDRAL_BACKBONE_ANALYSIS_NOTES.md:7-25`, `S1/R/39_dihedral_backbone_analysis.R` row counts in `S1/docs/SECONDARY_EXPLORATION_STATUS.md:307-318`.
- Over supported phi/psi bins, do raw atom-native calculator-kernel intensities co-vary with ORCA shielding target surfaces? Source: `S1/derived/dihedral_conjunction_probe_20260518_001/CONJUNCTION_PROBE_NOTES.md:7-18`.
- Are negative source-target surface correlations real anti-phase signal rather than failure? Source: `S1/derived/directionality_phase_audit_20260518_003/DIRECTIONALITY_PHASE_AUDIT_NOTES.md:7-16`.
- In the residue-local N-CA-C frame, what nearest-neighbour clouds differ between high and low within-protein ORCA shielding states? Source: `S1/derived/local_geometry_cloud_probe_20260518_002/LOCAL_GEOMETRY_CLOUD_NOTES.md:9-27`, `S1/R/45_local_geometry_cloud_probe.R:226-287`.
- Around phi and psi torsion axes, which neighbour classes at which angles are associated with higher or lower within-protein ORCA shielding z for a target backbone atom? Source: `OPUS/CHART_DESIGN.md:7-30`, `OPUS/derived/torsion_axis_v2_20260520_001/TORSION_AXIS_V2_NOTES.md:8-33`.
- In pair-derived torsion path models, do target-level contributor densities explain ORCA shielding beyond geometry-only covariates, and which terms have stored coefficients? Source: `S1/derived/torsion_pair_path_rollup_20260520_003/TORSION_PAIR_PATH_ROLLUP_NOTES.md:6-22`.

### Stage 2 Trajectory Scaffolds

- 1P9J trajectory scripts are scaffolds for dynamic follow-up, not completed positional results in this cluster: sigma variance vs `(1 - NH-S^2)`, sigma essential dynamics, fast-sigma vs Rex spatial colocation, and block-averaged SEM convergence. Source: `STAGE2/sigma_variance_vs_nh_s2.py:12-32`, `STAGE2/sigma_essential_dynamics.py:19-35`, `STAGE2/markwick_overlay_1p9j.py:5-20`, `STAGE2/block_avg_convergence.py:14-29`.

## 2. Methods / Math

### Export Contract / Strata

- Python/R must consume typed facts exported by C++ through `nmr_extract.load()` and avoid chemistry inference from filenames, regexes, display labels, or array order. Source: `S1/README.md:16-24`.
- Core support gate for named atom categories: `n_atoms >= 80` and `n_proteins >= 20`; target IDs are `total_iso`, `dia_iso`, `para_iso`, `total_t2_norm`, `dia_t2_norm`, `para_t2_norm`. Source: `S1/R/01_block_ridge.R:12-15`, `S1/R/01_block_ridge.R:23-26`, `S1/R/01_block_ridge.R:35-47`.
- Main row-count contract: 720 proteins; 425,599 matched atom rows; 475,116 topology atoms; 138,744 block-ridge rows; 275,232 dimension-source rows; 35,136 dimension-target rows; 58,560 per-atom top contributor rows; 667,584 per-atom dimension proxy rows; 1,464,000 blocked-stability samples. Source: `S1/docs/STAGE1_COMPENDIUM.md:259-290`.

### Scaling

- Global z scaling, per column:
  - `z_ij = (x_ij - mean_j) / sd_j` when at least 3 finite values and `sd_j > 0`; otherwise output 0 for that column. Source: `S1/R/matrix_common.R:47-62`.
- Within-protein scaling:
  - for protein `p`, `z_ij^(p) = (x_ij - mean_{p,j}) / sd_{p,j}` using the same finite-value and nonzero-SD rules; protein groups with fewer than 3 rows output 0. Source: `S1/R/matrix_common.R:64-73`.
- In R/36, "raw" is not unscaled physical units; it is global column centering/scaling via `safe_scale`. Source: `S1/R/36_backbone_atom_multisource_blocks.R:223-229`, audit caveat in `S1/docs/SECONDARY_EXPLORATION_STATUS.md:254-263`.

### Ridge Signal-Capture Diagnostics

- Ridge prediction:
  - `beta_hat_lambda = (X^T X + lambda I)^(-1) X^T y`
  - `y_hat = X beta_hat_lambda`
  - default `lambda = 0.01` in core block-ridge work. Source: `S1/R/matrix_common.R:75-88`, `S1/R/01_block_ridge.R:15`.
- Diagnostic fit metric:
  - `R^2 = 1 - sum_i (y_i - yhat_i)^2 / sum_i (y_i - mean(y))^2`
  - finite paired values only; fewer than 3 finite paired values returns `NA`. Source: `S1/R/matrix_common.R:90-115`.
- Block-ridge rows fit one block `b=(mechanism, group, component_kind)` at a time inside one named atom stratum, normalization, and ORCA target. Source: `S1/R/01_block_ridge.R:28-33`, `S1/R/01_block_ridge.R:64-100`.
- Protein-blocked CV holds out whole proteins in cyclic folds over sorted protein IDs; raw scaling is fit on train and applied to test; within-protein scaling is computed within each protein because folds are protein-disjoint. Source: `S1/R/07_mechanism_blocked_cv.R:86-137`.
- Blocked-stability resampling samples whole proteins without replacement, with `sample_fraction=0.80`, `n_resamples=25`, refits selected top blocks, and records sample `R^2` and sample rank. Source: `S1/R/21_blocked_stability.R:15-29`, `S1/R/21_blocked_stability.R:105-111`, `S1/R/21_blocked_stability.R:243-339`.

### Independence / Overlap

- Mechanism unique contribution:
  - `R2_full = R^2(X_all, y)`
  - `R2_without_m = R^2(X_all without mechanism m, y)`
  - `unique_delta_R2(m) = R2_full - R2_without_m`
  - `redundancy_gap(m) = R2_marginal(m) - unique_delta_R2(m)`. Source: `S1/R/02_mechanism_independence.R` summarized in `S1/docs/SCRIPT_PROVENANCE.md:15-17`.
- Pairwise block diagnostics:
  - `pair_R2 = R^2([pred_a, pred_b], y)`
  - `overlap = R2_a + R2_b - pair_R2`
  - `unique_a_given_b = pair_R2 - R2_b`
  - `unique_b_given_a = pair_R2 - R2_a`
  - `gain_over_best = pair_R2 - max(R2_a, R2_b)`. Source: `S1/R/12_atom_pairwise_independence.R` output definitions in `S1/docs/ATOM_SIGNAL_QUESTION_TABLES.md:18-35`.
- Positive participation of top-20 contributors:
  - `D_eff = (sum_i max(r_i,0))^2 / sum_i max(r_i,0)^2`
  - This is an effective positive-contributor count, not an SVD dimension. Source: `S1/R/13_atom_signal_question_tables.R:57-64`.

### SVD / Dimension Basis

- For scaled `X`, compute truncated SVD `X = U D V^T`, with `k = min(max_pcs, n-1, p)` and default `max_pcs=12`. Source: `S1/R/03_dimension_basis.R:9-16`, `S1/R/03_dimension_basis.R:60-73`.
- Dimension scores are `U_k D_k`; variance fraction is `d_j^2 / sum(X^2)`; source loading mass for a block is `sum_{features in block} v_{fj}^2 / sum_f v_{fj}^2`; target `R^2` is ridge from the first `u` PC scores to target `y`. Source: `S1/R/03_dimension_basis.R:102-170`.
- Same-dimension proxy rows pair the top source on a dimension with secondary sources on that same dimension; `secondary_to_top_loading_ratio = secondary_loading_l2_fraction / top_loading_l2_fraction`; relation is `same_source_type` or `cross_source_type`. Source: `S1/R/13_atom_signal_question_tables.R:282-343`.

### Backbone Atom Multi-Source Models

- Full model rows fit all 788 feature columns for each backbone atom role, normalization, and target. Source: `S1/R/36_backbone_atom_multisource_blocks.R:413-434`.
- Source-block rows store `single_r2`, `drop_source_r2`, `drop_delta_r2 = full_r2 - drop_source_r2`, fractions of full `R^2`, and prediction SD; coefficient vectors are not persisted in these R/36 CSVs. Source: `S1/R/36_backbone_atom_multisource_blocks.R:258-295`.
- Shared prediction space standardizes each source-block prediction vector, computes its correlation matrix, eigenvalues, `first_pc_variance_fraction`, `n_pcs_80pct`, and `full_minus_best_single_r2`. Source: `S1/R/36_backbone_atom_multisource_blocks.R:487-528`.
- Top-four commonality fits all nonempty subsets of the top source blocks and applies inclusion-exclusion to partition `R^2` into unique/shared/suppressor components. Source: `S1/R/36_backbone_atom_multisource_blocks.R:562-638`.

### Dihedral / Placement Surface Probes

- The dihedral sidecar is not mutation-delta; it reads all atom rows plus raw per-protein ORCA total/diamagnetic/paramagnetic arrays. Source: `S1/derived/dihedral_backbone_analysis_20260518_002/DIHEDRAL_BACKBONE_ANALYSIS_NOTES.md:7-15`.
- Phi and psi are residue-level DSSP values broadcast to atoms; first residues are excluded for phi, last residues excluded for psi; angular values are wrapped to `(-pi, pi]`; smooth grids mask low effective sample count. Source: `S1/derived/dihedral_backbone_analysis_20260518_002/DIHEDRAL_BACKBONE_ANALYSIS_NOTES.md:16-25`.
- Main feature sets include `phi_only`, `psi_only`, `phi_psi`, `phi_psi_interaction`, `omega`, `dssp_sasa_ss`, chi sets, and `all_local_dihedral`; fits use the same ridge diagnostic and `lambda=0.01`. Source: `S1/R/39_dihedral_backbone_analysis.R` row counts and run identity in `S1/docs/SECONDARY_EXPLORATION_STATUS.md:307-318`.
- Conjunction source intensity per atom-native raw calculator block:
  - scale block columns in the atom role by z.
  - `block_intensity_i = sqrt(mean_j z_{ij}^2)`
  - `block_signed_mean_i = mean_j z_{ij}`.
  - Placement arrays (`dssp_backbone`, `dssp_chi`, `dssp_ss8`, omega) and non-atom-native projected arrays are excluded. Source: `S1/derived/dihedral_conjunction_probe_20260518_001/CONJUNCTION_PROBE_NOTES.md:7-18`.
- Response-shape probe compares source and target surfaces with weighted correlations, rank correlations, z-cosine, z-RMSE, high-high overlap, adjacent-bin gradients, region contrasts, and theme-pair covariation. Source: `S1/derived/dihedral_response_shape_probe_20260518_002/RESPONSE_SHAPE_PROBE_NOTES.md:8-23`.
- Directionality audit classifies positive, negative, and near-zero direction separately; negative correlation is retained as anti-phase signal. Source: `S1/derived/directionality_phase_audit_20260518_003/DIRECTIONALITY_PHASE_AUDIT_NOTES.md:7-16`.

### Local Geometry Cloud / Rank Composition

- Local geometry cloud parameters: cutoff `5.50 A`, max neighbours `24`, voxel bin `1.0 A`, radial bin `0.5 A`; nearest neighbours require distance `>0.35 A` and `<=5.5 A`. Source: `S1/R/45_local_geometry_cloud_probe.R:55-58`, `S1/R/45_local_geometry_cloud_probe.R:386-397`.
- Target roles are `bb_N`, `bb_CA`, `bb_C`, `bb_O`, `bb_H`, `bb_HA`; targets are `total_iso`, `total_t2_norm`. Source: `S1/R/45_local_geometry_cloud_probe.R:226-230`.
- Target within-protein z is computed within `(protein_id, atom_role_label, target_id)`. High and low classes use `target_z >= q75` and `target_z <= q25`; middle values remain in continuous signal outputs. Source: `S1/R/45_local_geometry_cloud_probe.R:256-287`, caveat in `S1/derived/local_geometry_cloud_probe_20260518_002/LOCAL_GEOMETRY_CLOUD_NOTES.md:18-27`.
- Residue-local frame:
  - target atom at origin.
  - `e_x = unit(CA - N)`.
  - `e_y = unit((C - CA) - projection_on_e_x)`.
  - `e_z = cross(e_x, e_y)`.
  - neighbour local coordinates are dot products against `(e_x,e_y,e_z)`. Source: `S1/R/45_local_geometry_cloud_probe.R:336-424`.
- Continuous outputs aggregate `mean_target_z`, `mean_abs_target_z`, target value, distance by voxel, XY, XZ, radial shell, and nearest-neighbour rank. Source: `S1/R/45_local_geometry_cloud_probe.R:448-529`.
- High/low density deltas:
  - `P_high(b) = n_pairs_high(b) / sum_b n_pairs_high(b)` within `(target_atom_role,target_id,high)`.
  - `P_low(b) = n_pairs_low(b) / sum_b n_pairs_low(b)`.
  - `delta_high_low(b) = P_high(b) - P_low(b)`.
  - `log2_high_low = log2((P_high+1e-9)/(P_low+1e-9))`. Source: `S1/R/45_local_geometry_cloud_probe.R:758-817`.
- Geometry-cloud rank composition computes element probabilities by nearest-neighbour rank for high and low shielding classes, currently plotted for `bb_CA` and `bb_O`, target IDs `total_iso` and `total_t2_norm`, ranks `1-12`, elements `C/N/O/H`. Source: `S1/R/45_local_geometry_cloud_probe.R:1034-1072`.

### Torsion-Axis v2 / Angle-Contributor Charts

- Torsion-axis v2 stores one row per `(target_atom_role, target_id, placement_region, torsion, radius_bin, z_along_bin, theta_bin, contributor_class)`, with cylindrical coordinates around the phi or psi axis. Source: `OPUS/derived/torsion_axis_v2_20260520_001/TORSION_AXIS_V2_NOTES.md:8-18`.
- Phi frame: origin `N(i)`, z-axis `N(i)->CA(i)`, x-reference `C(i-1)`; psi frame: origin `CA(i)`, z-axis `CA(i)->C(i)`, x-reference `N(i)`. Source: `OPUS/CHART_DESIGN.md:17-24`.
- v2 fixes v1 by keeping full per-class breakdown instead of one plurality label, storing target-level sums/SE/t statistics, cylindrical bins, scaffold flags, and 0.50 A x 0.50 A x 5 degree bins. Source: `OPUS/derived/torsion_axis_v2_20260520_001/TORSION_AXIS_V2_NOTES.md:34-47`.
- Cell statistics:
  - `mean_target_z = sum_target_z / n_target_instances`.
  - `var_target_z` and `se_target_z` are meaningful only when `n_target_instances >= 2`.
  - `t_target_z = mean_target_z / se_target_z`.
  - v2 chart significance gate is `|t| >= 2.00`; cells below the gate still render faded in the chart. Source: `OPUS/derived/torsion_axis_v2_20260520_001/TORSION_AXIS_V2_NOTES.md:19-33`, `OPUS/derived/angle_contributor_chart_v2_20260521_008/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md:7-10`, `OPUS/derived/angle_contributor_chart_v2_20260521_008/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md:47-78`.
- Chart reading:
  - fill = signed statistic / sign of within-protein z direction.
  - size/count evidence is unique target instances, not raw pair count.
  - cumulative signed evidence around the ring distinguishes sinusoidal closed-ring from asymmetric broken-ring phase. Source: `OPUS/derived/angle_contributor_chart_v2_20260521_008/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md:11-39`, `OPUS/derived/angle_contributor_chart_v2_20260521_008/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md:63-104`.

### Torsion Pair Path Models

- Pair path rows keep target-neighbour provenance, then aggregate to target-level contributor density features before fitting, so one target shielding value is not counted once per neighbour. Source: `S1/derived/torsion_pair_path_rollup_20260520_003/TORSION_PAIR_PATH_ROLLUP_NOTES.md:10-14`.
- Model form per scene:
  - `y_z = G gamma + C beta + eps`
  - `G` includes phase/dihedral/SASA/SS geometry covariates.
  - `C` includes top contributor density features.
  - reported fits include geometry-only, contributors-only, and geometry-plus-contributors; ridge `lambda=1.0` in R/58. Source: `S1/R/58_torsion_pair_path_analysis.R` summarized by `S1/docs/SCRIPT_PROVENANCE.md:74-93` and `S1/derived/torsion_pair_path_rollup_20260520_003/TORSION_PAIR_PATH_ROLLUP_NOTES.md:16-22`.
- Term CSVs persist `full_model_beta_ols`, `full_model_beta_ridge`, `unique_r2_drop_ols`, and `unique_r2_drop_ridge`; these are model coefficients/diagnostics, not calibrated physical constants. Source: `S1/derived/torsion_pair_path_rollup_20260520_003/path_rollup_top_terms_by_atom_target.csv`.

## 3. Concrete Results

### 3.1 Row-Count Contract

| Analysis block | Count(s) | Source |
|---|---:|---|
| Feature export | 720 OK directories; 425,599 matched atom rows; 788 feature columns; 23 target columns | `S1/docs/STAGE1_COMPENDIUM.md:259-266` |
| Core ridge/dimension work | 138,744 block-ridge rows; 275,232 dimension-source rows; 35,136 dimension-target rows | `S1/docs/STAGE1_COMPENDIUM.md:267-271` |
| Per-atom contributor/pair/proxy | 58,560 top-contributor rows; 556,320 pairwise rows; 667,584 dimension proxy-pair rows | `S1/docs/STAGE1_COMPENDIUM.md:272-275` |
| Blocked stability | 2,928 atom groups; 1,464,000 blocked-stability samples; 0 reproduction check failures | `S1/docs/STAGE1_COMPENDIUM.md:276-290` |
| Dihedral sidecar | 45,036 dihedral fit rows; 41,472 smooth cells; 31,080 supported smooth cells | `S1/docs/STAGE1_COMPENDIUM.md:291-296` |
| Placement/response/directionality | 3,360 conjunction rows; 3,360 response-shape rows; 3,360 directionality rows; 186 moderate/strong negative examples | `S1/docs/STAGE1_COMPENDIUM.md:297-335` |
| Local geometry cloud | 167,706 target backbone atoms; 8,049,472 target-neighbour pairs; 20,838 nearest-rank signal rows; 5,045 nearest-rank delta rows | `S1/docs/STAGE1_COMPENDIUM.md:336-399` |
| Torsion-axis v2 | 18,285,392 pair rows; 1,017,768 `(cell,class)` rows; 294,377 cells with `|t| >= 2`; 72 scenes; 8/8 checks | `OPUS/README.md:55-56` |
| v2 bb_N/total_iso chart | 6 scenes; 3,673 significant cells; cumulative sums alpha_R phi/psi `+177/-113`, extended beta `+277/-17`, left alpha `+5019/+3591` | `OPUS/README.md:57-58`, detailed CSV below |

Caption: these counts prove table coverage and support, not physical correctness by themselves.

### 3.2 Backbone Atom Full-Feature Ridge, Total Targets

Every number in this table is from `S1/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_full_model_ridge.csv`; method source: `S1/R/36_backbone_atom_multisource_blocks.R:413-434`.

| Atom role | Normalization | Target | n rows | n proteins | Full diagnostic R^2 |
|---|---|---|---:|---:|---:|
| bb_C | raw/global-z | total_iso | 27,475 | 720 | 0.551 |
| bb_C | raw/global-z | total_t2_norm | 27,475 | 720 | 0.780 |
| bb_C | within_protein_z | total_iso | 27,475 | 720 | 0.469 |
| bb_C | within_protein_z | total_t2_norm | 27,475 | 720 | 0.834 |
| bb_CA | raw/global-z | total_iso | 27,475 | 720 | 0.615 |
| bb_CA | raw/global-z | total_t2_norm | 27,475 | 720 | 0.896 |
| bb_CA | within_protein_z | total_iso | 27,475 | 720 | 0.619 |
| bb_CA | within_protein_z | total_t2_norm | 27,475 | 720 | 0.878 |
| bb_H | raw/global-z | total_iso | 25,812 | 720 | 0.967 |
| bb_H | raw/global-z | total_t2_norm | 25,812 | 720 | 0.919 |
| bb_H | within_protein_z | total_iso | 25,812 | 720 | 0.925 |
| bb_H | within_protein_z | total_t2_norm | 25,812 | 720 | 0.887 |
| bb_HA | raw/global-z | total_iso | 24,960 | 720 | 0.976 |
| bb_HA | raw/global-z | total_t2_norm | 24,960 | 720 | 0.923 |
| bb_HA | within_protein_z | total_iso | 24,960 | 720 | 0.950 |
| bb_HA | within_protein_z | total_t2_norm | 24,960 | 720 | 0.902 |
| bb_N | raw/global-z | total_iso | 27,475 | 720 | 0.725 |
| bb_N | raw/global-z | total_t2_norm | 27,475 | 720 | 0.867 |
| bb_N | within_protein_z | total_iso | 27,475 | 720 | 0.536 |
| bb_N | within_protein_z | total_t2_norm | 27,475 | 720 | 0.789 |
| bb_O | raw/global-z | total_iso | 27,475 | 720 | 0.490 |
| bb_O | raw/global-z | total_t2_norm | 27,475 | 720 | 0.798 |
| bb_O | within_protein_z | total_iso | 27,475 | 720 | 0.374 |
| bb_O | within_protein_z | total_t2_norm | 27,475 | 720 | 0.775 |

Critical caption:
- Hydrogen and HA show very high in-sample signal capture, especially `total_iso`; this is not a calibration law and not an out-of-sample claim.
- Carbonyl O remains lower for `total_iso`, especially within-protein z, but has strong `total_t2_norm` signal.
- Nitrogen is not one thing: backbone-N total-iso is lower under within-protein z than raw/global-z, while total-T2 remains high.

### 3.3 Shared Prediction Space / Effective Source Complexity

Every number in this table is from `S1/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_shared_prediction_space.csv`; method source: `S1/R/36_backbone_atom_multisource_blocks.R:487-528`.

| Atom role | Target view | n sources | median `n_pcs_80pct` over 12 target-normalization cells | min-max `n_pcs_80pct` | median abs source-prediction corr, total targets | Meaning |
|---|---:|---:|---:|---:|---:|---|
| bb_C | all 6 targets x 2 norms | 18 | 9 | 8-9 | 0.070-0.094 | Many weakly correlated source prediction directions; no single "3-D" summary. |
| bb_CA | all 6 targets x 2 norms | 18 | 9 | 8-10 | 0.063-0.074 | Broad multi-source structure with high full model R^2. |
| bb_H | all 6 targets x 2 norms | 18 | 8 | 7-9 | 0.070-0.099 | High signal but still multiple prediction axes. |
| bb_HA | all 6 targets x 2 norms | 18 | 8 | 8-10 | 0.068-0.110 | High signal with nontrivial axis mixture. |
| bb_N | all 6 targets x 2 norms | 18 | 8 | 8-9 | 0.076-0.109 | Backbone-N remains multi-mechanism; do not merge with sidechain N. |
| bb_O | all 6 targets x 2 norms | 17-18 | 8 | 7-9 | 0.094-0.152 | Lower iso R^2 but still multiple axes; one hbond group is all-zero in R/38. |

Separate top-20 positive participation, not SVD dimension:

Every number in this table is from `S1/derived/tables/ml_atom_signal_shape_numeric.csv`; participation formula source: `S1/R/13_atom_signal_question_tables.R:57-64`.

| Atom role | raw total_iso median participation | raw total_t2_norm median participation | within total_iso median participation | within total_t2_norm median participation |
|---|---:|---:|---:|---:|
| bb_C | 13.15 | 11.48 | 11.21 | 8.27 |
| bb_CA | 11.63 | 12.89 | 10.14 | 8.46 |
| bb_H | 12.08 | 14.11 | 8.88 | 8.30 |
| bb_HA | 11.66 | 13.91 | 9.43 | 8.41 |
| bb_N | 10.82 | 12.59 | 11.39 | 8.59 |
| bb_O | 10.83 | 9.98 | 10.10 | 8.28 |

Critical caption:
- The files I found do not contain the exact phrase or table `H~20, C~6, N~3, O~12`; the audited current tables expose `n_pcs_80pct` and top-20 positive participation separately. Source checked by targeted `rg` over `S1/docs`, `S1/R`, and `AA`.
- Do not average these into one dimension count. The story is atom-role/target/normalization-specific.

### 3.4 Backbone Atom Signal / Stability / Proxy Overview

Every number in this table is from `S1/derived/backbone_atom_overlap_proxy_20260517_003/bb_atom_type_overview.csv`; method source: `S1/R/35_backbone_atom_overlap_proxy.R` row-count contract in `S1/docs/STAGE1_COMPENDIUM.md:307-318`.

| Atom role | Contributor rows | Stable independent rows | Weak stable rows | High-signal unstable rows | Median row R^2 | q90 row R^2 | Median unique | q90 unique | Stable pair-overlap rows | Proxy rows | Proxy rows loading >=0.25 | Proxy rows both top3>=0.50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bb_C | 3,840 | 571 | 4 | 1,337 | 0.099 | 0.480 | 0.094 | 0.439 | 239 | 43,776 | 4,530 | 2,127 |
| bb_CA | 3,840 | 573 | 5 | 1,353 | 0.100 | 0.526 | 0.095 | 0.497 | 175 | 43,776 | 4,326 | 2,058 |
| bb_H | 3,600 | 535 | 0 | 1,372 | 0.110 | 0.603 | 0.099 | 0.547 | 139 | 41,040 | 9,264 | 1,042 |
| bb_HA | 3,600 | 537 | 3 | 1,426 | 0.114 | 0.624 | 0.108 | 0.570 | 130 | 41,040 | 6,036 | 1,447 |
| bb_N | 3,840 | 569 | 0 | 1,443 | 0.109 | 0.477 | 0.100 | 0.428 | 194 | 43,776 | 5,058 | 1,974 |
| bb_O | 3,840 | 569 | 8 | 1,157 | 0.086 | 0.456 | 0.083 | 0.426 | 215 | 43,776 | 5,292 | 2,149 |

Critical caption:
- Stable-independent means `top3_rate >= 0.50` and `unique >= 0.01` in the readout thresholds, not a physical law. Source: `S1/R/26_blocked_stability_readout.R` thresholds summarized in `S1/docs/SECONDARY_EXPLORATION_STATUS.md:170-177`.
- `median_top3_rate` in this overview is 0 for all six atom roles because it is a median over all contributor rows, many of which are not top-stable; use the stable row counts above for thresholded readout.

### 3.5 Top Source Blocks by Atom Role

Every number in this table is from `S1/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv`; source-block method source: `S1/R/36_backbone_atom_multisource_blocks.R:258-295`.

| Atom role | Source block | Physics theme | n features | median single R^2 | q90 single R^2 | max single R^2 | q90 drop-delta R^2 | cells single R^2 >=0.20 |
|---|---|---|---:|---:|---:|---:|---:|---:|
| bb_C | ring_current | ring_response | 121 | 0.444 | 0.661 | 0.721 | 0.107 | 12 |
| bb_C | efg_coulomb | electrostatic_efg | 34 | 0.277 | 0.648 | 0.657 | 0.016 | 6 |
| bb_C | mopac_efg | electrostatic_efg | 34 | 0.254 | 0.623 | 0.647 | 0.022 | 6 |
| bb_CA | ring_current | ring_response | 121 | 0.514 | 0.762 | 0.780 | 0.107 | 12 |
| bb_CA | efg_coulomb | electrostatic_efg | 34 | 0.258 | 0.640 | 0.709 | 0.015 | 6 |
| bb_CA | mopac_efg | electrostatic_efg | 34 | 0.240 | 0.621 | 0.702 | 0.024 | 6 |
| bb_H | ring_current | ring_response | 121 | 0.565 | 0.894 | 0.915 | 0.108 | 12 |
| bb_H | mopac_bond_anisotropy | bond_anisotropy | 40 | 0.102 | 0.819 | 0.883 | 0.005 | 3 |
| bb_H | bond_anisotropy | bond_anisotropy | 40 | 0.101 | 0.818 | 0.884 | 0.005 | 3 |
| bb_HA | ring_current | ring_response | 121 | 0.594 | 0.924 | 0.951 | 0.098 | 12 |
| bb_HA | bond_anisotropy | bond_anisotropy | 40 | 0.131 | 0.843 | 0.900 | 0.003 | 4 |
| bb_HA | mopac_bond_anisotropy | bond_anisotropy | 40 | 0.134 | 0.840 | 0.897 | 0.003 | 4 |
| bb_N | ring_current | ring_response | 121 | 0.456 | 0.644 | 0.653 | 0.099 | 12 |
| bb_N | efg_coulomb | electrostatic_efg | 34 | 0.312 | 0.613 | 0.627 | 0.015 | 6 |
| bb_N | mopac_efg | electrostatic_efg | 34 | 0.302 | 0.591 | 0.602 | 0.017 | 6 |
| bb_O | efg_coulomb | electrostatic_efg | 34 | 0.292 | 0.641 | 0.650 | 0.015 | 6 |
| bb_O | ring_current | ring_response | 121 | 0.383 | 0.637 | 0.641 | 0.092 | 11 |
| bb_O | mopac_efg | electrostatic_efg | 34 | 0.277 | 0.617 | 0.623 | 0.020 | 6 |

Critical caption:
- `drop_delta_R^2` measures signal lost when a source block is removed from the all-feature model; it is not a physical coefficient.
- Ring-current blocks are strong across all backbone roles, but electrostatic EFG blocks are also large for C/N/O; the story is not "ring only".

### 3.6 Dihedral Feature-Set Ridge, Backbone Roles

Every number in this table is from `S1/derived/dihedral_backbone_analysis_20260518_002/dihedral_backbone_ranked_effects.csv`; method/caveat source: `S1/derived/dihedral_backbone_analysis_20260518_002/DIHEDRAL_BACKBONE_ANALYSIS_NOTES.md:7-25`.

| Atom role | Normalization | Target | n rows | phi-only R^2 | psi-only R^2 | phi+psi R^2 | phi+psi+interaction R^2 | DSSP/SASA/SS R^2 | all-local-dihedral R^2 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| bb_CA | raw_global_z | total_iso | 29,944 | 0.293 | 0.079 | 0.302 | 0.342 | 0.194 | 0.694 |
| bb_O | raw_global_z | total_iso | 29,944 | 0.118 | 0.384 | 0.395 | 0.455 | 0.295 | 0.502 |
| bb_C | raw_global_z | total_iso | 29,944 | 0.225 | 0.130 | 0.273 | 0.306 | 0.217 | 0.450 |
| bb_N | raw_global_z | total_iso | 29,944 | 0.154 | 0.007 | 0.208 | 0.232 | 0.066 | 0.446 |
| bb_HA | raw_global_z | total_iso | 27,429 | 0.223 | 0.094 | 0.238 | 0.254 | 0.201 | 0.382 |
| bb_H | raw_global_z | total_iso | 28,268 | 0.011 | 0.037 | 0.039 | 0.070 | 0.210 | 0.264 |
| bb_CA | raw_global_z | total_t2_norm | 29,944 | 0.215 | 0.052 | 0.227 | 0.321 | 0.199 | 0.580 |
| bb_N | raw_global_z | total_t2_norm | 29,944 | 0.083 | 0.021 | 0.113 | 0.145 | 0.050 | 0.395 |
| bb_O | raw_global_z | total_t2_norm | 29,944 | 0.040 | 0.095 | 0.140 | 0.180 | 0.079 | 0.363 |
| bb_HA | raw_global_z | total_t2_norm | 27,429 | 0.020 | 0.044 | 0.083 | 0.099 | 0.054 | 0.202 |
| bb_C | raw_global_z | total_t2_norm | 29,944 | 0.014 | 0.066 | 0.114 | 0.126 | 0.005 | 0.173 |
| bb_H | raw_global_z | total_t2_norm | 28,268 | 0.005 | 0.006 | 0.014 | 0.027 | 0.046 | 0.091 |
| bb_CA | within_protein_z | total_iso | 29,944 | 0.243 | 0.046 | 0.247 | 0.295 | 0.158 | 0.649 |
| bb_N | within_protein_z | total_iso | 29,944 | 0.153 | 0.011 | 0.195 | 0.221 | 0.069 | 0.460 |
| bb_O | within_protein_z | total_iso | 29,944 | 0.080 | 0.304 | 0.312 | 0.373 | 0.199 | 0.429 |
| bb_C | within_protein_z | total_iso | 29,944 | 0.172 | 0.056 | 0.206 | 0.233 | 0.143 | 0.399 |
| bb_HA | within_protein_z | total_iso | 27,429 | 0.170 | 0.051 | 0.187 | 0.196 | 0.126 | 0.335 |
| bb_H | within_protein_z | total_iso | 28,268 | 0.011 | 0.029 | 0.030 | 0.064 | 0.177 | 0.230 |
| bb_CA | within_protein_z | total_t2_norm | 29,944 | 0.180 | 0.020 | 0.190 | 0.295 | 0.150 | 0.547 |
| bb_N | within_protein_z | total_t2_norm | 29,944 | 0.080 | 0.029 | 0.105 | 0.140 | 0.057 | 0.411 |
| bb_O | within_protein_z | total_t2_norm | 29,944 | 0.028 | 0.058 | 0.102 | 0.150 | 0.044 | 0.318 |
| bb_HA | within_protein_z | total_t2_norm | 27,429 | 0.038 | 0.038 | 0.077 | 0.098 | 0.053 | 0.204 |
| bb_C | within_protein_z | total_t2_norm | 29,944 | 0.015 | 0.066 | 0.107 | 0.118 | 0.006 | 0.157 |
| bb_H | within_protein_z | total_t2_norm | 28,268 | 0.008 | 0.006 | 0.021 | 0.031 | 0.046 | 0.098 |

Critical caption:
- Dihedral geometry carries substantial within-protein signal for bb_CA, bb_N, bb_O, and bb_C, but this is a diagnostic correlation to ORCA targets in regular proteins.
- bb_O total-iso is psi-led; bb_N total-iso is phi-led; bb_H is much weaker in pure phi/psi but has DSSP/SASA/SS signal. Do not average across atom roles.

### 3.7 Placement Surface / Response Shape / Directionality

Pattern counts from `S1/derived/dihedral_conjunction_probe_20260518_001/conjunction_bin_alignment_summary.csv`, `S1/derived/dihedral_response_shape_probe_20260518_002/response_shape_similarity.csv`, and `S1/derived/directionality_phase_audit_20260518_003/directionality_phase_summary.csv`; method sources: `CONJUNCTION_PROBE_NOTES.md:7-18`, `RESPONSE_SHAPE_PROBE_NOTES.md:8-23`, `DIRECTIONALITY_PHASE_AUDIT_NOTES.md:7-16`.

| Probe | Rows | Pattern counts |
|---|---:|---|
| Conjunction bin alignment | 3,360 | weak_or_mixed 2,166; signed_covariation 548; placement_modulated_distinct 303; moderate_abs_covariation 260; strong_abs_covariation 83 |
| Response shape similarity | 3,360 | weak_or_mixed_shape 1,604; moderate_intensity_shape 803; moderate_signed_shape 409; strong_signed_shape 261; source_modulated_distinct_shape 200; strong_intensity_shape 83 |
| Directionality phase audit | 3,360 | low_visible 1,167; magnitude_visible 676; weak_or_mixed 637; signed_visible 564; magnitude_and_signed_visible 316 |

Selected high-strength directionality rows, every number from `S1/derived/directionality_phase_audit_20260518_003/directionality_phase_summary.csv`:

| Atom | Target | Norm | Source/theme | n bins | Magnitude direction/strength | Signed direction/strength | Signed tail phase | Relation type |
|---|---|---|---|---:|---|---|---|---|
| bb_C | dia_iso | raw_global_z | electronic_embedding | 27 | positive 0.506 | negative 0.861 | antiphase_tail_support | magnitude_and_signed_visible |
| bb_HA | dia_iso | raw_global_z | hbond_solvation | 25 | negative 0.193 | negative 0.897 | antiphase_tail_support | signed_visible |
| bb_C | para_iso | raw_global_z | electronic_embedding | 27 | positive 0.474 | positive 0.860 | aligned_tail_support | magnitude_and_signed_visible |
| bb_CA | dia_iso | raw_global_z | polarisability | 27 | positive 0.849 | negative 0.617 | antiphase_tail_support | magnitude_and_signed_visible |
| bb_CA | dia_iso | raw_global_z | hbond_solvation | 27 | near_zero 0.082 | negative 0.818 | antiphase_tail_support | signed_visible |
| bb_C | para_iso | raw_global_z | hbond_solvation | 27 | negative 0.294 | positive 0.811 | aligned_tail_support | signed_visible |

Critical caption:
- Negative signed correlations are retained as anti-phase signal; they are not discarded as weak rows.
- Co-high overlap answers aligned high-high questions only; anti-phase tails must be read with the anti-phase columns.

### 3.8 Local Geometry Cloud, Radial Signal

Every number in this table is from `S1/derived/local_geometry_cloud_probe_20260518_002/geometry_cloud_radial_signal_summary.csv`; method source: `S1/R/45_local_geometry_cloud_probe.R:448-529`.

Rows shown require at least 5,000 neighbour pairs and are sorted by `max_abs_bin_target_z`; sparse-bin risk still remains, so use `n_pairs`, `n_target_instances`, and `rms_bin_target_z` together.

| Target atom | Target | Region | Neighbour | Relation | radial bins | n pairs | n target instances | signed mean target z | RMS bin target z | max abs bin target z |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|
| bb_HA | total_t2_norm | extended_beta_box | C | same_residue | 10 | 57,554 | 40,769 | -0.033 | 0.148 | 6.109 |
| bb_HA | total_t2_norm | extended_beta_box | N | next_residue | 6 | 12,356 | 12,356 | 0.011 | 0.183 | 6.109 |
| bb_C | total_iso | extended_beta_box | H | same_residue | 9 | 61,676 | 39,084 | 0.057 | 0.077 | 6.071 |
| bb_HA | total_iso | extended_beta_box | C | same_residue | 10 | 57,554 | 40,769 | -0.200 | 0.235 | 5.933 |
| bb_HA | total_iso | extended_beta_box | N | next_residue | 6 | 12,356 | 12,356 | -0.236 | 0.307 | 5.933 |
| bb_C | total_iso | alpha_R_box | H | next_residue | 8 | 38,993 | 36,814 | -0.099 | 0.250 | 4.987 |
| bb_C | total_t2_norm | alpha_R_box | H | next_residue | 8 | 38,993 | 36,814 | -0.183 | 0.190 | 4.878 |
| bb_N | total_iso | alpha_R_box | H | previous_residue | 7 | 46,668 | 23,914 | -0.008 | 0.131 | 4.608 |
| bb_CA | total_iso | alpha_R_box | H | same_residue | 9 | 84,157 | 45,350 | -0.438 | 0.513 | 4.461 |
| bb_CA | total_t2_norm | alpha_R_box | H | same_residue | 9 | 84,157 | 45,350 | -0.315 | 0.354 | 3.787 |

Critical caption:
- This table says supported local neighbour shells have nonzero signed within-protein ORCA shielding-z structure.
- It does not say the neighbour causes the shielding state; the cloud is an indicative geometry diagnostic.

### 3.9 Local Geometry Cloud, Nearest-Rank High-Low Delta

Every number in this table is from `S1/derived/local_geometry_cloud_probe_20260518_002/geometry_cloud_nearest_rank_delta.csv`; method source: `S1/R/45_local_geometry_cloud_probe.R:758-817`.

Rows shown require `n_pairs_total >= 1,000` and are sorted by absolute `probability_delta_high_low`.

| Target atom | Target | Rank | Neighbour | Relation | n high | n low | probability delta high-low | log2 high/low | n total |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|
| bb_H | total_iso | 2 | C | previous_residue | 6,704 | 718 | 0.036 | 3.224 | 7,422 |
| bb_H | total_iso | 3 | C | previous_residue | 187 | 6,088 | -0.036 | -5.024 | 6,275 |
| bb_O | total_t2_norm | 2 | H | distant_sequence | 383 | 6,330 | -0.035 | -4.047 | 6,713 |
| bb_O | total_t2_norm | 2 | N | next_residue | 6,404 | 503 | 0.035 | 3.671 | 6,907 |
| bb_H | total_iso | 2 | O | distant_sequence | 1 | 5,488 | -0.033 | -12.420 | 5,489 |
| bb_O | total_t2_norm | 3 | N | next_residue | 715 | 6,365 | -0.033 | -3.154 | 7,080 |
| bb_O | total_t2_norm | 3 | C | same_residue | 5,827 | 272 | 0.032 | 4.421 | 6,099 |
| bb_H | total_iso | 3 | C | same_residue | 4,699 | 84 | 0.028 | 5.806 | 4,783 |
| bb_CA | total_iso | 14 | H | same_residue | 1,184 | 5,697 | -0.026 | -2.266 | 6,881 |
| bb_H | total_t2_norm | 2 | C | previous_residue | 2,465 | 6,588 | -0.025 | -1.418 | 9,053 |
| bb_H | total_t2_norm | 3 | C | previous_residue | 4,353 | 302 | 0.024 | 3.849 | 4,655 |
| bb_CA | total_iso | 7 | O | same_residue | 148 | 4,255 | -0.024 | -4.845 | 4,403 |

Critical caption:
- The `<=q25` low and `>=q75` high class definitions are inclusive, not strict `<`/`>`; source: `S1/R/45_local_geometry_cloud_probe.R:273-276`.
- Rank deltas are class-conditional probability differences, not raw count differences.

### 3.10 Torsion Pair Path Rollup

Every number in this table is from `S1/derived/torsion_pair_path_rollup_20260520_003/path_rollup_atom_target_summary.csv`; rollup method source: `S1/derived/torsion_pair_path_rollup_20260520_003/TORSION_PAIR_PATH_ROLLUP_NOTES.md:6-22`.

| Atom | Target | Scenes | median n targets | median full ridge R^2 | max full ridge R^2 | median geometry-only R^2 | median contributor-only R^2 | median unique contributor R^2 | median overlap R^2 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| bb_C | total_iso | 6 | 12,791 | 0.631 | 0.644 | 0.251 | 0.601 | 0.378 | 0.218 |
| bb_CA | total_iso | 6 | 12,791 | 0.785 | 0.893 | 0.153 | 0.751 | 0.632 | 0.096 |
| bb_H | total_iso | 6 | 12,304 | 0.552 | 0.635 | 0.153 | 0.544 | 0.436 | 0.131 |
| bb_HA | total_iso | 6 | 12,333 | 0.448 | 0.605 | 0.141 | 0.429 | 0.367 | 0.122 |
| bb_N | total_iso | 6 | 12,791 | 0.612 | 0.621 | 0.231 | 0.565 | 0.388 | 0.178 |
| bb_O | total_iso | 6 | 12,791 | 0.576 | 0.631 | 0.295 | 0.531 | 0.305 | 0.240 |
| bb_C | total_t2_norm | 6 | 12,791 | 0.411 | 0.427 | 0.149 | 0.355 | 0.225 | 0.084 |
| bb_CA | total_t2_norm | 6 | 12,791 | 0.689 | 0.710 | 0.524 | 0.522 | 0.178 | 0.215 |
| bb_H | total_t2_norm | 6 | 12,304 | 0.449 | 0.705 | 0.101 | 0.435 | 0.359 | 0.078 |
| bb_HA | total_t2_norm | 6 | 12,333 | 0.432 | 0.449 | 0.145 | 0.307 | 0.198 | 0.108 |
| bb_N | total_t2_norm | 6 | 12,791 | 0.445 | 0.465 | 0.200 | 0.332 | 0.265 | 0.099 |
| bb_O | total_t2_norm | 6 | 12,791 | 0.615 | 0.705 | 0.185 | 0.580 | 0.475 | 0.176 |

Selected coefficient-bearing terms, every number from `S1/derived/torsion_pair_path_rollup_20260520_003/path_rollup_top_terms_by_atom_target.csv`:

| Atom | Target | Level | Scene | Contributor | n targets | corr with shielding | single R^2 | ridge beta | unique ridge R^2 drop |
|---|---|---|---|---|---:|---:|---:|---:|---:|
| bb_H | total_iso | atom | extended beta / PSI | O / carbonyl O / distant sequence | 12,304 | -0.731 | 0.535 | -0.650 | 0.242 |
| bb_H | total_iso | atom | extended beta / PHI | O / carbonyl O / distant sequence | 12,304 | -0.731 | 0.535 | -0.648 | 0.241 |
| bb_H | total_iso | atom | alpha-R / PSI | O / carbonyl O / distant sequence | 13,038 | -0.611 | 0.373 | -0.578 | 0.161 |
| bb_H | total_iso | atom | alpha-R / PHI | O / carbonyl O / distant sequence | 13,038 | -0.611 | 0.373 | -0.572 | 0.157 |
| bb_H | total_iso | atom | left-alpha / PHI | O / carbonyl O / distant sequence | 1,418 | -0.549 | 0.302 | -0.363 | 0.067 |
| bb_H | total_iso | atom | left-alpha / PSI | O / carbonyl O / distant sequence | 1,418 | -0.549 | 0.302 | -0.360 | 0.066 |
| bb_C | total_t2_norm | residue_atom | extended beta / PHI | PRO:N / backbone N / next residue | 12,791 | -0.453 | 0.205 | -0.472 | 0.197 |
| bb_C | total_t2_norm | residue_atom | extended beta / PSI | PRO:N / backbone N / next residue | 12,791 | -0.453 | 0.205 | -0.470 | 0.195 |

Critical caption:
- These `beta` values are ridge model coefficients after feature scaling in a covariance/path model. They are not physical calibration constants in ppm or a causal mediation effect.

### 3.11 Torsion-Axis v2 bb_N/total_iso Chart

Every number in this table is from `OPUS/derived/angle_contributor_chart_v2_20260521_008/angle_contributor_v2_chart_manifest.csv`; chart method source: `OPUS/derived/angle_contributor_chart_v2_20260521_008/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md:1-39`.

| Scene | significant cells | cells total | cumulative sum | cumulative max abs | significant n<5 | significant n>=25 | median n among significant |
|---|---:|---:|---:|---:|---:|---:|---:|
| bb_N total_iso alpha_R phi | 680 | 1,440 | 176.6 | 191.8 | 43 | 514 | 71 |
| bb_N total_iso alpha_R psi | 429 | 1,440 | -113.1 | 110.9 | 57 | 284 | 78 |
| bb_N total_iso extended_beta phi | 762 | 1,440 | 276.7 | 197.3 | 32 | 621 | 113 |
| bb_N total_iso extended_beta psi | 531 | 1,440 | -16.5 | 198.6 | 67 | 336 | 74 |
| bb_N total_iso left_alpha phi | 733 | 1,440 | 5,019 | 378.8 | 122 | 309 | 18 |
| bb_N total_iso left_alpha psi | 538 | 1,440 | 3,591 | 400.5 | 111 | 212 | 16 |

Top rendered contributors by scene, every number from `OPUS/derived/angle_contributor_chart_v2_20260521_008/angle_contributor_v2_rank.csv`:

| Scene | Rank | Contributor class | target instances | pairs | signed evidence sum | abs evidence sum | significant bins |
|---|---:|---|---:|---:|---:|---:|---:|
| alpha_R phi | 1 | sidechain H, previous residue | 27,489 | 27,489 | -18.78 | 303.7 | 49 |
| alpha_R phi | 2 | backbone N, previous residue | 13,485 | 13,485 | -17.25 | 296.9 | 51 |
| alpha_R phi | 3 | alpha H, previous residue | 12,555 | 12,555 | -191.8 | 290.6 | 55 |
| alpha_R psi | 1 | sidechain H, same residue | 44,113 | 44,113 | -110.9 | 208.2 | 33 |
| alpha_R psi | 2 | carbonyl O, distant sequence | 8,279 | 8,279 | 37.38 | 171.5 | 30 |
| alpha_R psi | 3 | amide H, near sequence | 15,988 | 15,988 | 17.31 | 169.7 | 29 |
| extended_beta phi | 1 | alpha H, previous residue | 11,486 | 11,486 | -197.3 | 335.9 | 53 |
| extended_beta phi | 2 | backbone N, previous residue | 12,791 | 12,791 | 169.7 | 326.7 | 58 |
| extended_beta phi | 3 | sidechain H, same residue | 48,211 | 48,211 | 9.13 | 323.3 | 49 |
| extended_beta psi | 1 | sidechain H, same residue | 48,211 | 48,211 | -198.6 | 239.3 | 38 |
| extended_beta psi | 2 | carbonyl O, distant sequence | 4,810 | 4,810 | 77.18 | 194.8 | 34 |
| extended_beta psi | 3 | amide H, distant sequence | 5,282 | 5,282 | -46.17 | 174.3 | 28 |
| left_alpha phi | 1 | sidechain H, distant sequence | 916 | 916 | 378.8 | 378.8 | 63 |
| left_alpha phi | 2 | alpha C, same residue | 1,418 | 1,418 | 375.5 | 375.5 | 60 |
| left_alpha phi | 3 | sidechain H, previous residue | 3,219 | 3,219 | 368.3 | 368.3 | 50 |
| left_alpha psi | 1 | sidechain H, distant sequence | 916 | 916 | 400.5 | 400.5 | 62 |
| left_alpha psi | 2 | carbonyl C, same residue | 1,418 | 1,418 | 379.9 | 379.9 | 62 |
| left_alpha psi | 3 | sidechain H, near sequence | 705 | 705 | 275.7 | 275.7 | 40 |

Critical caption:
- v2 chart evidence is a sign/rank/phase visualization; absolute cumulative values scale with theta-bin width and are not physical magnitudes. Source: `OPUS/derived/angle_contributor_chart_v2_20260521_008/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md:27-39`.
- Left-alpha scenes have smaller median significant-cell support (`16-18`) than alpha-R or extended-beta scenes (`71-113`), so left-alpha visual patterns are stronger in cumulative signed evidence but rest on thinner support.

### 3.12 Stage 2 Trajectory Scripts

| Script | Question | Status in file | Source |
|---|---|---|---|
| `sigma_variance_vs_nh_s2.py` | Compare per-atom `sigma_T0` trajectory variance `(ppm)^2` to `(1 - NH-S^2)` mobility proxy for BMRB 5801 / 1P9J. | Emits per-atom variance; residue/BMRB bridge is TODO. | `STAGE2/sigma_variance_vs_nh_s2.py:12-32`, `STAGE2/sigma_variance_vs_nh_s2.py:116-121` |
| `sigma_essential_dynamics.py` | SVD/PCA of per-frame sigma tensor stack to test low-rank sigma manifold. | Placeholder output; per-calc TS wiring is scaffold. | `STAGE2/sigma_essential_dynamics.py:19-35`, `STAGE2/sigma_essential_dynamics.py:72-99` |
| `markwick_overlay_1p9j.py` | Spatial colocation of fast-timescale sigma fluctuation with micro/millisecond Rex; expected null due timescale mismatch. | Scaffold; per-calc shielding TS read not wired. | `STAGE2/markwick_overlay_1p9j.py:5-20`, `STAGE2/markwick_overlay_1p9j.py:103-143` |
| `block_avg_convergence.py` | Block-averaged SEM convergence for per-atom T0 timelines. | Implemented for named TS channels present in H5. | `STAGE2/block_avg_convergence.py:14-29`, `STAGE2/block_avg_convergence.py:105-139` |

## Gaps / What I Could Not Find

- I did not find a final physical calibration-law coefficient table in the surveyed `20260514` positional cluster. R/36 stores diagnostic ridge `R^2`, source drop deltas, prediction-space correlations, and commonality, but not source coefficient vectors. Source checked: `S1/derived/backbone_atom_multisource_blocks_20260517_003/*.csv`; method confirms no coefficient persistence in `S1/R/36_backbone_atom_multisource_blocks.R:258-295`.
- Coefficient-bearing positional rows were found in the torsion pair path rollup (`full_model_beta_ols`, `full_model_beta_ridge`), but those are scaled covariance/path model coefficients, not calibrated physical constants. Source: `S1/derived/torsion_pair_path_rollup_20260520_003/path_rollup_top_terms_by_atom_target.csv`.
- I did not find the exact "H~20, C~6, N~3, O~12" effective-dimension numbers in targeted searches over `S1/docs`, `S1/R`, and `AA`. Current audited files provide `n_pcs_80pct` and positive top-20 participation; both are reported above and should not be collapsed into one number.
- Stage 2 contains trajectory-analysis scaffolds, not completed 1P9J numerical results in the files read here. Source: `STAGE2/sigma_variance_vs_nh_s2.py:116-121`, `STAGE2/sigma_essential_dynamics.py:91-99`, `STAGE2/markwick_overlay_1p9j.py:136-143`.
- I did not read broad recursive `derived/`, `runs/`, `figures/`, or `png/` trees. I used named docs/scripts and named CSVs only, per the non-interactive constraint.
