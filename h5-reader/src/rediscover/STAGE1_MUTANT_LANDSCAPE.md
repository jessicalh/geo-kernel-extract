# Stage-1 Mutant Landscape, Recent-Only

This is a read-only landscape map of the recent Stage-1 mutant-analysis work. It is not a conclusion, not a ranking, and not a recommendation to merge any Stage-1 result into the current rediscover engine.

Source root used below:

`LEARN=/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn`

Hard scope:

- Read only files with `mtime >= 2026-05-12`.
- Do not treat older work as resolved background unless a recent file itself summarized it.
- Do not investigate the IUPAC/topology episode. I noted the matching-policy provenance where it was part of recent run accounting, but did not spelunk it.
- Do not read the whole slow tree. I read recent docs, run ledgers, current summary tables, targeted CSV rollups, and selected notes. Large matrices, PDFs, images, and most per-file raw outputs were not opened.

## One-Screen Map

The recent Stage-1 mutant work is best understood as a sequence of progressively stricter views over the same broad experimental substrate: many feature families, many atom identities, many target components, two normalization regimes, and several forms of blocking/stability analysis. The central pattern is not "mechanism X wins." The central pattern is that visibility changes with the question being asked:

- Raw/global magnitude views often keep ring response, AIMNet2/charge embeddings, and high-R2 apparent structure visible.
- Within-protein-Z views often reveal or amplify angular/local/tensor structure, EFG-like proxies, and some geometry/secondary-structure effects.
- Atom-type and named-atom views refine element-level stories. In particular, "nitrogen is hard" was later refined into "backbone N is hard in some contexts; sidechain N can be excellent."
- Stability, blocked-CV, pair-overlap, and multi-source drop-one tests repeatedly changed the interpretation of high single-source R2.
- Several apparent signals became "real but conditional," "visible but proxy-like," or "present in one target/norm/atom stratum and absent in another."

The current state is therefore a plural map:

- Ring response is broad, high-reach, tensor-specific in many tables, and especially important in raw/global views, but it is also overlapped with ring dispersion, ring EFG/pi-quadrupole, secondary structure, and AIMNet2-like embeddings.
- Charge/electrostatic features are a mixture: AIMNet2 charge embeddings are powerful but proxy/learned; simple EEQ/MOPAC core charges are often tiny in multi-source backbone blocks; Coulomb/MOPAC EFG are visible in some Stage-1 tables; APBS EFG is mostly not load-bearing in the targeted backbone multi-source read.
- McConnell/bond-anisotropy is not dead, but Stage-1 treats it as a conditional, overlapping mechanism. It appears in feature definitions, weak/independence/proxy tables, and MOPAC-vs-classical variants, but in the May 17 backbone multi-source block its drop-one contribution is small compared with ring response and charge/EFG themes.
- EFG/field-like material is a family, not one thing. Stage-1 includes Coulomb E-field, Coulomb EFG, MOPAC E-field/EFG, AIMNet2 EFG, APBS E-field/EFG, and APBS mutation-delta features. Apparent EFG visibility can be tensor-specific and normalisation-sensitive, but blocked-CV and multi-source checks caution against reading every high EFG row as a robust physical law.
- H-bond and solvation are mostly weak/sparse in the Stage-1 mutant tables, with one important hard-null: `bb_O/hbond_grid` was diagnosed as structurally all-zero in the May 17 backbone multi-source diagnostic.
- The geometry/torsion sidecars moved the work from feature leaderboards toward placement/phase maps. They are regular-protein/raw-shielding sidecars, not mutation-delta mechanism decompositions, and their notes explicitly ask readers not to over-read them causally.

## Date-Ordered Supersession Trail

### 2026-05-14: export and first full tables

The recent workup starts from a large validated export: 720 mutation directories, 425,599 matched atom rows, 788 feature columns, and 23 target columns. The compendium frames the project as "what shielding signal is visible after topology-aware atom identity" rather than a global leaderboard. It explicitly tracks exact named atom, normalization, target, and calculator/source rows. Citation: `LEARN/stage1-topology-workup-20260514/docs/STAGE1_COMPENDIUM.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/docs/SCRIPT_PROVENANCE.md` (mtime 2026-05-20).

The May 14 core tables include:

- block ridge and mechanism views,
- mechanism independence,
- dimension basis/source loading,
- mechanism blocked CV,
- residue/topology coverage,
- feature catalog and source dictionary.

Important early tension: apparent mechanism strength and blocked/stable mechanism strength are not the same thing. For example, mechanism independence tables show strong marginal and unique charge/ring patterns, while blocked-CV summaries show some EFG-like apparent rows collapsing under protein blocking. Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/mechanism_independence_summary.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/tables/mechanism_blocked_cv_summary.csv` (mtime 2026-05-14).

### 2026-05-15: secondary explorations and adversarial audit

The next layer asked whether first-pass summaries were robust to confounds:

- normalization contrast: raw/global vs within-protein-Z,
- weak independent signal,
- source overlap graphs,
- tensor/component split,
- high-R2 audit,
- atom role stratification,
- ML packaging views.

Citation: `LEARN/stage1-topology-workup-20260514/docs/SECONDARY_EXPLORATIONS.md` (mtime 2026-05-15).

The adversarial audit found and fixed several interpretation-risk bugs:

- source-family mapping had treated DSSP inconsistently in some scripts,
- tensor split used the wrong absolute-magnitude calculation,
- family-node named-atom counts were overcounted,
- empty threshold pairs were hidden,
- the high-R2 smoke pass was biased and later superseded by balanced full25 blocked-stability chunks,
- append guards and shape/range checks were added.

This matters because some attractive early reads were later downgraded from "signal" to "needs stability/support" or "view artifact." Citation: `LEARN/stage1-topology-workup-20260514/docs/ADVERSARIAL_AUDIT.md` (mtime 2026-05-15).

### 2026-05-15 to 2026-05-16: atom detail and blocked-stability full25

Atom-detail views extended the work beyond element pooling. The table docs warn that downstream views should not refit and should preserve named atom, normalization, target, and source. They also record singleton/invalid named-atom support groups, e.g. MET H1/H2/H3 and GLY O'' in within-protein singleton named-atom groups. Citation: `LEARN/stage1-topology-workup-20260514/docs/ATOM_SIGNAL_QUESTION_TABLES.md` (mtime 2026-05-15).

The balanced full25 blocked-stability rollup then superseded the earlier high-R2 smoke pass. It contains 58,560 contributor rows, 2,928 atom groups, 58,560 pair rows, and 667,584 dimension-proxy rows, with 65 readout checks passing. It reports:

- 8,695 stable independent contributors,
- 47 weak-stable independent contributors,
- 23,460 high observed single-signal rows with unstable top3 rank,
- 2,446 stable pair-overlap rows.

The important change is epistemic: high observed R2 became only one axis. Stable top3, unique contribution, support warnings, pair gain, and blocked structure became equally important. Citation: `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/blocked_stability_readout.md` (mtime 2026-05-16); `LEARN/stage1-topology-workup-20260514/docs/SECONDARY_EXPLORATION_STATUS.md` (mtime 2026-05-20).

### 2026-05-17: signal visibility, relationship synthesis, results packet

The current signal-visibility packet refined the question again. It did not just ask "which source has high R2"; it tracked row reach, exact source, stable source counts, high-unstable rows, low-correlation rows, raw-vs-Z normalization contrast, tensor specificity, and source relationships.

Key current packet outputs include:

- `results_v1_20260517_008` as the current narrowed results packet after several superseded versions.
- `relationship_report_003` as the current relationship synthesis.
- `signal_visibility_atlas_20260516_003` as the current atlas after row-count tightening.

Citation: `LEARN/stage1-topology-workup-20260514/docs/SECONDARY_EXPLORATION_STATUS.md` (mtime 2026-05-20); `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/RESULTS_PACKET_NOTES.md` (mtime 2026-05-17).

The source-family summary is broad rather than reductive:

- Ring has very large row reach, many exact-source rows, and many high rows, but also many unstable/high rows.
- AIMNet2 has high median and stable counts, but is a learned/proxy ceiling, not a physical mechanism.
- EFG has a meaningful stable footprint and strong tensor specificity in some summaries.
- MOPAC and MOPAC bond-anisotropy are present but mixed.
- H-bond, solvation, and simple charges are sparse/weak in several aggregate views.

Citation: `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_source_family_summary.csv` (mtime 2026-05-17).

### 2026-05-17 to 2026-05-18: backbone atom multi-source blocks and hard nulls

The backbone multi-source blocks were a major refinement because they performed true ridge refits from NPY inputs and then asked conditional questions such as "what happens when this theme is dropped from the full block?" Their notes warn that raw is global-Z, drop-one deltas are conditional, and prediction correlations omit constants. Citation: `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/BB_ATOM_MULTISOURCE_BLOCKS_NOTES.md` (mtime 2026-05-18).

This block makes several nuances concrete:

- Ring response has the largest median single-R2 and largest conditional drop-one deltas across the six backbone atom roles in this block.
- Electrostatic EFG has moderate single-R2 but smaller conditional drop deltas.
- Charge embedding has moderate single-R2 and modest conditional deltas, largely tied to AIMNet2 in source-block detail.
- Bond anisotropy has modest single-R2 and very small drop deltas in this particular backbone multi-source block.
- H-bond/solvation is weak by both single and drop-one measures.

For example, median drop-delta R2 by theme was roughly:

| Atom role | Ring response | Charge embedding | Electrostatic EFG | Bond anisotropy | Hbond/solvation |
|---|---:|---:|---:|---:|---:|
| bb_C | 0.102 | 0.042 | 0.035 | 0.005 | 0.001 |
| bb_CA | 0.085 | 0.029 | 0.035 | 0.004 | 0.001 |
| bb_H | 0.077 | 0.021 | 0.034 | 0.004 | 0.001 |
| bb_HA | 0.068 | 0.021 | 0.034 | 0.004 | 0.001 |
| bb_N | 0.104 | 0.038 | 0.030 | 0.004 | 0.002 |
| bb_O | 0.098 | 0.039 | 0.029 | 0.005 | 0.001 |

Citation: `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_physics_theme_summary.csv` (mtime 2026-05-18).

The source-block detail matters for APBS and EFG:

- `efg_coulomb` and `mopac_efg` have moderate single-R2 in the backbone multi-source block.
- `efg_apbs` is tiny in this block, with median single-R2 around 0.0007 to 0.0025 across bb atoms and zero drop-delta threshold counts.
- AIMNet2 EFG is also small in this source-block table, despite AIMNet2 as a broader embedding being strong elsewhere.

Citation: `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18).

The hard null was also made explicit: `bb_O/hbond_grid` produced constant predictions because it had zero raw nonzero columns, zero scaled nonzero columns, rank 0, and prediction SD 0. Citation: `LEARN/stage1-topology-workup-20260514/derived/backbone_constant_prediction_diagnostic_20260517_002/CONSTANT_PREDICTION_DIAGNOSTIC_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/docs/SECONDARY_EXPLORATION_STATUS.md` (mtime 2026-05-20).

### 2026-05-18 to 2026-05-20: geometry, directionality, torsion, and phase sidecars

The later sidecars shifted from mutation-delta mechanism scoring toward regular-protein geometry and phase questions. They are explicitly sidecars, not final mechanism decompositions.

Important constraints from the notes:

- Dihedral/backbone analysis uses regular-protein raw shielding, not mutation-delta.
- Smooth surface masks must not be read as evidence in low-support angular regions.
- Conjunction probes compare placement and calculator-block intensity but do not residualize placement and do not claim causal attribution.
- Directionality/phase audit treats negative correlation as signal direction, not failure.
- Rotational geometry model notes describe possible local-frame and spherical/Wigner/cyclic bases, but do not claim a new law.
- Local-geometry clouds are visual and target-centered; high-low density differences are visualization aids, not gates.
- Torsion-pair path analysis is covariance/provenance, not causal mediation.

Citations: `LEARN/stage1-topology-workup-20260514/derived/dihedral_backbone_analysis_20260518_002/DIHEDRAL_BACKBONE_ANALYSIS_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/dihedral_conjunction_probe_20260518_001/CONJUNCTION_PROBE_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/directionality_phase_audit_20260518_003/DIRECTIONALITY_PHASE_AUDIT_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/directionality_phase_audit_20260518_003/ROTATIONAL_GEOMETRY_MODEL_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/local_geometry_cloud_probe_20260518_002/LOCAL_GEOMETRY_CLOUD_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/torsion_pair_path_rollup_20260520_003/TORSION_PAIR_PATH_ROLLUP_NOTES.md` (mtime 2026-05-20).

The open-ended thread here is important: the Stage-1 work increasingly treats geometry as the selection variable. It is not enough to name a residue or IUPAC category. The later tools ask where a contributor sits in a supported local frame, around a torsion axis, or inside a Ramachandran-supported placement region.

### 2026-05-20 to 2026-05-22: derived-opus v2 and chart substrate repair

The Opus chart work is a clear supersession example. The v1 chart was described as "viva winner" style, but the newer README lists six substrate weaknesses:

1. dominant-label plurality collapsed contributor mixtures,
2. signed density used pair-count sample size over target-level estimator,
3. no per-cell evidence weight,
4. Cartesian voxelization broke rotational symmetry,
5. scaffold atoms were not flagged,
6. resolution was coarse.

The v2 pipeline replaced this with cylindrical cells, actual contributor class, unique target counts, t-stat gates, scaffold flags, 5-degree bins, and evidence panels. Citation: `LEARN/stage1-topology-workup-20260514/derived-opus/README.md` (mtime 2026-05-21); `LEARN/stage1-topology-workup-20260514/derived-opus/derived/angle_contributor_chart_v2_20260522_009/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md` (mtime 2026-05-22).

The v2 chart notes explicitly constrain the claim:

- load-bearing reads are sign of per-cell within-protein z, phase structure around torsion axis, and rank ordering across scenes;
- the chart does not claim physical magnitudes;
- signed-t scales with bin width;
- rows, not single bubbles, should be read;
- evidence-distribution panels expose low-n bins.

Citation: `LEARN/stage1-topology-workup-20260514/derived-opus/derived/angle_contributor_chart_v2_20260522_009/ANGLE_CONTRIBUTOR_V2_CHART_NOTES.md` (mtime 2026-05-22).

The latest 2026-05-22 chart manifests preserve atom-specific phase complexity:

- bb_N default/focus page: left-alpha phi/psi cumulative sums are very large positive; alpha-R and extended-beta are much smaller and mixed/near-zero in sign.
- bb_H: positive in alpha-R, extended-beta, and left-alpha scenes.
- bb_HA: positive in alpha-R and left-alpha, negative in extended-beta.
- bb_CA: negative in alpha-R, strongly positive in extended-beta and left-alpha.
- bb_C: positive across the shown total_iso scenes, especially left-alpha.
- bb_O: negative in alpha-R and left-alpha, positive in extended-beta.

Citations: `LEARN/stage1-topology-workup-20260514/derived-opus/derived/angle_contributor_chart_v2_20260522_009/angle_contributor_v2_chart_manifest.csv` (mtime 2026-05-22); `LEARN/stage1-topology-workup-20260514/derived-opus/derived/angle_contributor_chart_v2_full_atoms_20260522_009/*/angle_contributor_v2_chart_manifest.csv` (mtime 2026-05-22); `LEARN/stage1-topology-workup-20260514/derived-opus/CHART_DESIGN.md` (mtime 2026-05-21).

This is not a mechanism conclusion. It is a later geometry/phase view that complicates any single "ring effect" or "backbone atom effect" story.

### 2026-05-25: physics bridge, actual-physics overview, and kernel definitions

Two recent files summarize/refine older physics interpretations:

- The bridge maps 55 kernels to 8 physics groups and an element-specific predictive-dimensionality story.
- The actual-physics overview reframes element-level interpretations using atom-type calibration.

Citations: `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25); `LEARN/src/actual_physics/OVERVIEW.md` (mtime 2026-05-25).

The bridge says raw kernel space can look universally low-dimensional/magnitude-dominated, while normalization strips magnitude and reveals element-specific angular structure. It gives an element story:

- H: ring current dominant in normalized physics bridge, with high normalized R2 and ring-current contribution.
- C: lower-dimensional but EFG-dominated in that bridge; charge-polarization gap is specifically a carbon phenomenon.
- N: blurred/weak in element-pooled views, with near-field exceptions.
- O: richer/high-dimensional, with dispersion becoming much more visible after normalization.

Citation: `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25).

The actual-physics overview then refines the nitrogen story: sidechain N can be excellent while backbone N is hard. Citation: `LEARN/src/actual_physics/OVERVIEW.md` (mtime 2026-05-25).

The May 25 kernel definitions show the Stage-1 feature universe was deliberately broad: ring current, ring susceptibility, ring dispersion, pi-quadrupole/ring EFG, bond anisotropy, MOPAC bond anisotropy, Coulomb/MOPAC/APBS/AIMNet2 EFG and fields, hbond, Larsen hbond grid, solvation/water field, AIMNet2 dyadic lifts, per-ring residuals, radial RBF partitions, angular partitions, and mutation scalar features. Citations: `LEARN/src/mutation_set/kernels.py` (mtime 2026-05-25); `LEARN/src/mutation_set/scalars.py` (mtime 2026-05-25); `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14).

## Mechanism Landscape

### Ring response

Ring response is one of the broadest and most repeatedly visible families, but it should not be collapsed into one law.

What Stage-1 tested under ring-like headings:

- Biot-Savart ring current,
- Haigh-Mallion ring current,
- ring susceptibility,
- ring dispersion,
- pi-quadrupole/ring EFG,
- per-ring residualized variants,
- radial RBF partitions,
- angular axial/equatorial partitions,
- tensor and scalar decompositions.

Citations: `LEARN/src/mutation_set/kernels.py` (mtime 2026-05-25); `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14).

Where it carries:

- It has the largest row reach and many stable/high rows in the May 17 results packet.
- In backbone multi-source blocks, ring response has the highest median single-R2 and largest median conditional drop deltas across bb_C, bb_CA, bb_H, bb_HA, bb_N, and bb_O.
- Raw/global views often favor ring because magnitude/proximity is part of the signal.
- Tensor/component summaries show ring can be strongly tensor-specific, not only scalar isotropic.
- Pair-overlap tables show ring current often overlaps with ring dispersion, ring EFG, secondary structure, and AIMNet2-like embeddings.

Citations: `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_source_family_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_physics_theme_summary.csv` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_tensor_component_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_stable_pair_overlap.csv` (mtime 2026-05-16).

Where it fails or becomes ambiguous:

- Many high ring rows are high observed but unstable in blocked/top3 checks.
- Ring strength changes under normalization. The May 17 normalization summary shows many ring rows raw-only and a negative median raw-to-Zp delta, which means within-protein-Z can suppress the magnitude piece that made ring obvious.
- Ring-current, dispersion, ring EFG, and geometry can co-vary, making source-family attribution difficult.

Citations: `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_normalization_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_high_signal_unstable_rank.csv` (mtime 2026-05-16); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_dimension_proxy_stability_summary.csv` (mtime 2026-05-16).

Bridge to current rediscover:

- The current rediscover clean Pople ring law is compatible with Stage-1 ring being broad and real, but Stage-1 warns that mutation-delta ring visibility can include magnitude, tensor, dispersion, ring EFG, and geometry overlap. The current de-circularised ring law should not be overextended into every Stage-1 ring row.

### Bond anisotropy / McConnell

Stage-1 includes both classical McConnell and MOPAC-derived McConnell/bond anisotropy variants.

What was tested:

- `mc_scalars`, `mc_category_T2`, `mc_shielding`,
- `mopac_mc_scalars`, `mopac_mc_category_T2`, `mopac_mc_shielding`,
- scalar and tensor forms,
- exact source, source-family, weak-independent, proxy, and multi-source conditional views.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14); `LEARN/src/mutation_set/kernels.py` (mtime 2026-05-25).

Where it appears:

- It is present in the dimension source mixture tables. Top-mechanism counts include bond anisotropy especially for hydrogens and sidechain hydrogens, with smaller counts for backbone roles.
- It appears in weak-independent and source-overlap/proxy tables, including McConnell vs MOPAC-McConnell pairs and McConnell overlap with ring/EFG families.
- In the backbone multi-source block, bond anisotropy has moderate median single-R2 for C/CA/H/HA but small conditional drop deltas across all six backbone roles.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/core_dimension_source_mixture_by_atom.csv` (mtime 2026-05-15); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/weak_independent_summary.csv` (mtime 2026-05-15); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_dimension_proxy_scope_summary.csv` (mtime 2026-05-16); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18).

Where it fails or becomes ambiguous:

- Conditional drop-one effect in the May 17 backbone block is small, which means bond anisotropy is often absorbed by other blocks in that design.
- It is not clear from Stage-1 alone whether McConnell is a standalone observable law, a joint-model contributor, a geometry-restricted mechanism, or partly a proxy for neighboring electronic structure.
- MOPAC-McConnell variants do not simply rescue all cases; they track classical McConnell in some overlap tables and diverge in others.

Bridge to current rediscover:

- Current rediscover's fair valid-source McConnell standalone failure is not contradicted by Stage-1. Stage-1's McConnell evidence is from mutation-delta/static/all-protein contexts and mixed source families; it suggests McConnell belongs in joint/ensemble or geometry-controlled tests, not as a single-mechanism standalone expectation.
- The Stage-1 kernel/filter design also supports the current need to avoid self/near-field source leakage. `bones/agent-brief.md` describes source filters such as self-source, sequential exclusion, and near-field filters. Citation: `LEARN/bones/agent-brief.md` (mtime 2026-05-25).

### Charge and electrostatic embedding

Stage-1 separates several things that are easy to conflate:

- AIMNet2 charge/electronic embeddings,
- EEQ charges,
- MOPAC core charge,
- Coulomb E-field and shielding,
- Coulomb EFG,
- MOPAC Coulomb E-field/EFG,
- APBS E-field/EFG,
- mutation-delta APBS features.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14); `LEARN/src/mutation_set/scalars.py` (mtime 2026-05-25).

Where it carries:

- Early mechanism-independence tables show charges with strong marginal and unique patterns in some within-protein-Z targets.
- Dimension source mixture tables are dominated by `charges|aimnet2|scalar` as a top loading across many named atom/dimension rows.
- Backbone multi-source theme summaries show charge embedding with moderate median single-R2 and modest conditional drop deltas.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/mechanism_independence_summary.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/tables/core_dimension_source_mixture_by_atom.csv` (mtime 2026-05-15); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_physics_theme_summary.csv` (mtime 2026-05-18).

Where it becomes ambiguous:

- In the May 17 backbone source-block summary, `charges_eeq` and `mopac_core_charge` are tiny compared with `aimnet2` charge embedding.
- AIMNet2 charge-like features are powerful but learned/proxy-like; they should not be read as simple q/r3 physics.
- Charge and EFG families overlap heavily in source-proxy tables.

Citations: `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_dimension_proxy_scope_summary.csv` (mtime 2026-05-16).

Bridge to current rediscover:

- Current rediscover's q/r3 backbone static-T2 law is a more specific mathematical extraction than Stage-1's broad "charge embedding" family. Stage-1 supports keeping charge/electrostatic mechanisms alive, but it does not by itself distinguish simple q/r3 from learned AIMNet2 or EFG-adjacent proxies.

### EFG, field, Buckingham-like material, and APBS

Stage-1's EFG/field family is broad and not identical to the current rediscover Buckingham/APBS question.

What was tested:

- `coulomb_E`, `coulomb_efg_*`, `coulomb_shielding`,
- `mopac_coulomb_E`, `mopac_coulomb_efg_*`, `mopac_coulomb_shielding`,
- `aimnet2_efg`,
- `apbs_E`, `apbs_efg`,
- `delta_apbs`,
- water E-field/EFG and solvation features.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14).

Where it appears:

- EFG has meaningful row reach and stable counts in the May 17 source-family results.
- EFG is strongly tensor-specific in the tensor/component summary.
- Coulomb and MOPAC EFG have moderate single-R2 in the May 17 backbone multi-source source-block summary.
- In dimension source mixtures, electrostatic EFG appears as a top mechanism in some atom roles, especially backbone hydrogens and some C/CA rows.

Citations: `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_source_family_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_tensor_component_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/tables/core_dimension_source_mixture_by_atom.csv` (mtime 2026-05-15).

Where it fails or is contested:

- Mechanism blocked-CV showed an electrostatic/EFG family with high apparent R2 in some axes but low median CV, a warning that apparent EFG signal can be fragile under protein blocking.
- APBS specifically is not load-bearing in the May 17 backbone multi-source block: `efg_apbs` has tiny median single-R2 and no meaningful drop-delta threshold counts across bb atom roles.
- The source-proxy scope summary shows many APBS-related proxy pairs have very low fraction with both stability, e.g. charge/AIMNet2 vs APBS, Coulomb EFG vs APBS, MOPAC EFG vs APBS.
- Stage-1 APBS is folded into `electrostatic_efg`; it is not a clean test of current APBS radii, Coulomb prefactors, Buckingham axis projection, or current local-frame EFG formulation.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/mechanism_blocked_cv_summary.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_dimension_proxy_scope_summary.csv` (mtime 2026-05-16).

Bridge to current rediscover:

- Stage-1 suggests "EFG-like" and "field-like" signals are not uniformly absent, but APBS itself was not load-bearing in the targeted backbone multi-source block. This lines up with the current rediscover finding that APBS placeholder radii were not the main issue, but it does not solve the current field/EFG math question.
- The main warning is frame and confounder control: Stage-1 EFG visibility often lives in tensor, normalization, and source-overlap contexts. Current local-frame EFG failures should be treated as a separate mathematical/mechanistic problem rather than dismissed by Stage-1 aggregate EFG rows.

### Dispersion

Dispersion appears most clearly as ring dispersion and is part of the ring-response neighborhood, but the recent bridge gives it a distinct oxygen-heavy role.

Where it appears:

- Ring dispersion appears in feature catalog and kernel definitions as scalar, tensor, and shielding forms.
- It is present in stable pair-overlap and weak-independent views.
- The physics bridge says oxygen dispersion changes strongly after normalization, going from a weak raw signal to a much more important normalized contributor.
- In the backbone multi-source source-block table, ring dispersion has moderate median single-R2 and nontrivial drop counts for some atom roles, especially bb_N and bb_O.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_stable_pair_overlap.csv` (mtime 2026-05-16); `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18).

Open ambiguity:

- Dispersion is not cleanly isolated from ring proximity, aromatic geometry, ring current, ring EFG/pi-quadrupole, and atom-type response.
- The oxygen dispersion story is strong in the bridge but should be checked against the later atom-type and multi-source tables rather than imported as a universal law.

Bridge to current rediscover:

- If the current all-category emitter adds or prioritizes dispersion, Stage-1 suggests oxygen and normalization-controlled contexts are the first places to look, but not as a single winner.

### H-bond and solvation

The H-bond family includes physical hbond scalars/shielding, Larsen-grid hbond contributions, water term, and solvation/water-field features.

Where it appears:

- It is part of the broad feature catalog and source dictionary.
- Hbond kernel has some weak/stable rows in the blocked-stability readout.
- Secondary structure/DSSP and hbond-energy views can have weak independent signal, especially as geometry/proxy markers.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_weak_stable_independent_contributors.csv` (mtime 2026-05-16).

Where it fails:

- Hbond/solvation has sparse stability in the May 17 source-family summary.
- In the backbone multi-source block, hbond/solvation has low median single-R2 and negligible drop deltas.
- `bb_O/hbond_grid` was structurally all-zero in the constant-prediction diagnostic.

Citations: `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_source_family_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_physics_theme_summary.csv` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/backbone_constant_prediction_diagnostic_20260517_002/CONSTANT_PREDICTION_DIAGNOSTIC_NOTES.md` (mtime 2026-05-18).

Bridge to current rediscover:

- H-bond features should be emitted/audited, but Stage-1 warns that absence can be structural rather than mechanistic. Before interpreting a null, check nonzero columns, rank, and support.

### Quadrupole / ring EFG / pi-quadrupole

This family sits between ring response and electrostatic tensor physics.

Where it appears:

- Feature catalog includes pi-quadrupole scalar and T2 per ring type plus shielding.
- Pair-overlap and dimension-proxy tables show ring EFG/pi-quadrupole overlapping with ring current, dispersion, and EFG.
- It contributes to the bridge's N/O complexity story but rarely stands alone in aggregate summaries.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_dimension_proxy_scope_summary.csv` (mtime 2026-05-16); `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25).

Open ambiguity:

- It may be a real tensor/source component, a ring-proximity proxy, or both depending on atom, target, and normalization.

### AIMNet2 and learned/proxy ceilings

AIMNet2 is one of the most visible families in Stage-1, but its interpretation is deliberately different from the classical mechanisms.

Where it carries:

- May 17 source-family summary has high median R2 and many stable rows for AIMNet2.
- Blocked-stability top-family rows often put AIMNet2 at the top for several axes.
- Dimension source mixtures show high AIMNet2 loading fractions across many atom roles.
- Actual-physics overview reports strong atom-type calibration values, including sidechain N.

Citations: `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_source_family_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_source_family_by_surface_norm.csv` (mtime 2026-05-16); `LEARN/stage1-topology-workup-20260514/derived/tables/core_dimension_source_mixture_by_atom.csv` (mtime 2026-05-15); `LEARN/src/actual_physics/OVERVIEW.md` (mtime 2026-05-25).

Interpretation:

- AIMNet2 is a ceiling/proxy, not a recovered physical law. It can tell us signal is learnable and can expose missing variables, but it cannot by itself say whether charge, EFG, bond anisotropy, dispersion, or geometry is the physical cause.

Bridge to current rediscover:

- This agrees with current rediscover's AIMNet2 role as learnable ceiling. Use it to set expectations and identify residual learnability, not to replace mechanism recovery.

### APBS specifically

APBS appears in Stage-1 as solvated E-field, solvated EFG, and mutation-delta APBS E+EFG. Citation: `LEARN/stage1-topology-workup-20260514/derived/tables/feature_catalog.csv` (mtime 2026-05-14).

What the recent Stage-1 work suggests:

- APBS is present in the feature universe and source-proxy tables.
- In targeted backbone multi-source blocks, `efg_apbs` is small and non-load-bearing relative to Coulomb and MOPAC EFG.
- APBS has low co-stable fractions in several source-proxy rows, suggesting it is not a strong substitute for Coulomb/MOPAC EFG in those Stage-1 views.

Citations: `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/readout_dimension_proxy_scope_summary.csv` (mtime 2026-05-16).

What it does not settle:

- It does not settle current APBS placeholder radii, field prefactor, axis projection, or local-frame EFG math.
- It does not prove APBS cannot matter in other atom/geometry/target regimes, because the APBS read above is from specific Stage-1 tables and source-block design.

## Stratification Findings

### Element pooling can obscure atom-type reality

The older-looking element bridge gives H/C/N/O narratives, but recent atom-type material complicates them. The most explicit refinement is nitrogen:

- Element-pooled bridge: N looks blurred/weak, no group dominates strongly, near-field exceptions exist.
- Actual-physics overview: "Nitrogen hard" was wrong as a blanket statement; backbone N is hard, sidechain N can be excellent.
- Backbone atom-type overview: bb_N's aggregate row/stability distribution is not uniquely catastrophic compared with bb_C/bb_CA/bb_H/bb_HA/bb_O, but source/block behavior differs.

Citations: `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25); `LEARN/src/actual_physics/OVERVIEW.md` (mtime 2026-05-25); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_overlap_proxy_20260517_003/bb_atom_type_overview.csv` (mtime 2026-05-17).

### Named atom and residue views matter

The Stage-1 table design repeatedly preserves exact named atom, residue, backbone role, surface/region, target component, normalization, and source identity. This is not bookkeeping noise. It changes interpretation:

- singletons and low-support atom groups can create attractive but invalid rows,
- backbone and sidechain atoms have different mechanism mixes,
- residue-specific named atoms can have high signal for reasons that do not transfer to the element.

Citations: `LEARN/stage1-topology-workup-20260514/docs/ATOM_SIGNAL_QUESTION_TABLES.md` (mtime 2026-05-15); `LEARN/stage1-topology-workup-20260514/derived/tables/core_dimension_source_mixture_by_atom.csv` (mtime 2026-05-15).

### Raw vs within-protein-Z is not a cosmetic choice

Normalization changes the question:

- raw/global-Z preserves magnitude/protein-size/proximity effects,
- within-protein-Z emphasizes within-protein angular/local variation,
- ring often loses when moving to within-protein-Z,
- EFG/AIMNet2/secondary-structure rows can gain under within-protein-Z,
- some bond-anisotropy and APBS/EFG rows appear as Zp-only in normalization-contrast tables.

Citations: `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_normalization_summary.csv` (mtime 2026-05-17); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/norm_contrast_summary.csv` (mtime 2026-05-15); `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25).

### Geometry stratification is becoming the preferred sample axis

The later sidecars increasingly treat geometry as the unit of interpretation:

- phi/psi and DSSP placement,
- local target-centered frames,
- torsion-axis angular placement,
- scaffold/forced-geometry flags,
- evidence-weighted per-cell views.

The chart design explicitly says some angles are forced by geometry: theta for C(i) is phi and theta for N(i+1) is psi. This is exactly the kind of confound that residue/IUPAC names can obscure. Citation: `LEARN/stage1-topology-workup-20260514/derived-opus/CHART_DESIGN.md` (mtime 2026-05-21).

Bridge to current rediscover:

- This strongly supports selecting current rediscover scenarios by geometry, not by IUPAC/residue name. It also supports a distinction between maths-methods that extract signal from all data and sample-methods that choose clean geometric regions.

## Maths and Fitting Explorations

The recent work did not use one fit style. It used a family of mathematical probes:

- ridge fits over block/source features,
- mechanism independence and unique delta-R2,
- protein-blocked CV,
- balanced blocked-stability full25,
- top3 stability and support warnings,
- pairwise source overlap/gain,
- commonality and drop-one multi-source refits,
- tensor/component specificity,
- dimension/PCA-like source loading and incremental dimension rows,
- weak-independent low-R2 contributor retention,
- high-R2 audits,
- normalization contrast.

Citations: `LEARN/stage1-topology-workup-20260514/docs/STAGE1_COMPENDIUM.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/docs/SECONDARY_EXPLORATION_STATUS.md` (mtime 2026-05-20); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/BB_ATOM_MULTISOURCE_BLOCKS_NOTES.md` (mtime 2026-05-18).

Several methodological lessons recur:

- High R2 is not enough. It can be unstable, support-limited, or proxy-driven.
- Drop-one deltas are conditional on the full model and should not be read as marginal mechanism strength.
- Prediction correlations compare fitted predictions and omit constants; they are not raw feature correlations.
- Weak low-R2 rows can still be independently stable and should not be discarded automatically.
- Tensor-specific response can be real even if scalar summaries look modest.
- Protein blocking can overturn apparent rows.
- Dimension-proxy rows are not the same scope as contributor stability rows.

Citations: `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/blocked_stability_readout.md` (mtime 2026-05-16); `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/BB_ATOM_MULTISOURCE_BLOCKS_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/high_r2_audit_summary.csv` (mtime 2026-05-15).

Bridge to current rediscover:

- For current McConnell and field/EFG, Stage-1 argues against only asking "does standalone source X fit?" A fair maths-methods layer should include partial/joint/drop-one/controlled fits, with explicit support and overlap diagnostics.
- For current all-atom/all-category emit, Stage-1 argues for preserving exact source, source family, tensor component, target, normalization, atom role, named atom, and support metadata so later analyses do not have to reconstruct provenance.

## Dimension and Complexity

There is a real evolution and tension here.

The physics bridge gives a compact element-specific effective-dimension story:

- raw kernel space looks universally 3-dimensional and magnitude/protein-size confounded,
- normalized H is high-dimensional/ring-heavy,
- C is around 6-dimensional and EFG-heavy in that bridge,
- N is around 3-dimensional/blurred with near-field exceptions,
- O is around 12-dimensional with dispersion emerging strongly after normalization.

Citation: `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25).

The May 14/15 topology workup dimension tables complicate that:

- named-atom dimension target R2 often continues improving through 12 dimensions;
- the "best" row in the `dimension_target_r2.csv` sample is usually at 12 dimensions across element/norm/target combinations;
- dimension source mixtures are often dominated by AIMNet2 charge embeddings as top loading, with ring, EFG, bond-anisotropy, and dispersion appearing as secondary or later-dimension components;
- top-loading counts differ by atom role and surface.

Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/dimension_target_r2.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/tables/plot_dimension_target_summary.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/tables/core_dimension_source_mixture_by_atom.csv` (mtime 2026-05-15); `LEARN/stage1-topology-workup-20260514/derived/tables/core_dimension_incremental_independence.csv` (mtime 2026-05-15).

How to hold the ambiguity:

- The bridge may be describing a physics-grounded kernel grouping and effective predictive structure.
- The topology workup dimension tables may be showing named-atom/all-feature/proxy-loaded complexity, where AIMNet2 and source overlap inflate the dimensional landscape.
- These are not simply contradictory; they answer different questions. But they should not be collapsed into one clean H/C/N/O dimension count.

Bridge to current rediscover:

- Effective dimension should be reported with the fitting universe and stratification attached. For current all-category emit, a per-element dimension number without atom-role/source-family context is likely misleading.

## Hard Cases, Nulls, and Resolutions

Resolved or refined failures:

- DSSP source-family mapping inconsistency was fixed in the adversarial audit. Citation: `LEARN/stage1-topology-workup-20260514/docs/ADVERSARIAL_AUDIT.md` (mtime 2026-05-15).
- Tensor split magnitude bug was fixed. Citation: same audit.
- Family-node named-atom overcount was fixed. Citation: same audit.
- Hidden empty threshold pairs were made visible. Citation: same audit.
- High-R2 smoke bias was superseded by balanced full25 blocked stability. Citation: `LEARN/stage1-topology-workup-20260514/docs/SECONDARY_EXPLORATION_STATUS.md` (mtime 2026-05-20).
- Several results packet versions were superseded by later current packets: atlas `_003`, relationship report `_003`, results_v1 `_008`, backbone multi-source `_003`, multi-source audit `_002`. Citation: same status file.
- `bb_O/hbond_grid` was resolved as a structural all-zero/constant-prediction case, not an ambiguous model failure. Citation: `LEARN/stage1-topology-workup-20260514/derived/backbone_constant_prediction_diagnostic_20260517_002/CONSTANT_PREDICTION_DIAGNOSTIC_NOTES.md` (mtime 2026-05-18).
- Opus v1 chart substrate weaknesses were replaced by v2 chart design. Citation: `LEARN/stage1-topology-workup-20260514/derived-opus/README.md` (mtime 2026-05-21).

Still nuanced or unresolved:

- High observed signal with unstable rank remains common. It was not "fixed away"; it is part of the landscape.
- Weak but stable independent rows remain. They are not high-R2 results, but they can point to real minor mechanisms.
- EFG visibility versus local-frame/current rediscover EFG failure is unresolved because the Stage-1 and rediscover formulations differ.
- McConnell visibility in mutation-delta/static all-protein contexts versus current standalone local T2 failure is unresolved and should be tested jointly/geometrically.
- Geometry/phase charts are suggestive but not causal mechanism decompositions.

## Confounders Identified

The recent work repeatedly exposes these confounders:

- Normalization: raw magnitude vs within-protein local/angular variation.
- Element pooling: especially nitrogen; sidechain and backbone differ.
- Atom-type/named-atom support: singleton, low-support, and invalid groups can mislead.
- Source overlap: ring, dispersion, ring EFG, secondary structure, AIMNet2, charge, EFG, and MOPAC blocks overlap.
- Learned/proxy features: AIMNet2 can dominate without identifying the physical law.
- Tensor/component labels: tensor target bugs and component-specificity matter.
- High-R2 support: high rows can be unstable or support-limited.
- Geometry forcing: torsion-axis angles can be chemically/geometrically forced rather than independent.
- Sign/phase: negative correlation can be a directional signal, not a failure.
- Constant feature blocks: nulls can be structural matrix nulls, not mechanism nulls.
- Matching/topology policy: present in provenance, not investigated here per guardrail.

Citations: `LEARN/stage1-topology-workup-20260514/docs/ADVERSARIAL_AUDIT.md` (mtime 2026-05-15); `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/readout_002/blocked_stability_readout.md` (mtime 2026-05-16); `LEARN/stage1-topology-workup-20260514/derived-opus/CHART_DESIGN.md` (mtime 2026-05-21); `LEARN/stage1-topology-workup-20260514/derived/directionality_phase_audit_20260518_003/DIRECTIONALITY_PHASE_AUDIT_NOTES.md` (mtime 2026-05-18); `LEARN/stage1-topology-workup-20260514/derived/backbone_constant_prediction_diagnostic_20260517_002/CONSTANT_PREDICTION_DIAGNOSTIC_NOTES.md` (mtime 2026-05-18).

## Surprises and Reversals

- "Nitrogen is hard" was refined into a stratification artifact: sidechain N can be excellent while backbone N remains difficult. Citation: `LEARN/src/actual_physics/OVERVIEW.md` (mtime 2026-05-25).
- Oxygen dispersion became much more important after normalization in the bridge. Citation: `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25).
- EFG can look strong in apparent/tensor contexts and weak under blocked/local/current contexts. This is a tension, not a resolved contradiction. Citations: `LEARN/stage1-topology-workup-20260514/derived/tables/mechanism_blocked_cv_summary.csv` (mtime 2026-05-14); `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/results_v1_tensor_component_summary.csv` (mtime 2026-05-17).
- APBS was present as a field/EFG source but tiny in the backbone multi-source source-block table. Citation: `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/bb_atom_source_block_summary.csv` (mtime 2026-05-18).
- High-R2 smoke was explicitly downgraded by support/stability audits. Citation: `LEARN/stage1-topology-workup-20260514/docs/ADVERSARIAL_AUDIT.md` (mtime 2026-05-15).
- Negative phase/correlation was reframed as directional signal rather than failure. Citation: `LEARN/stage1-topology-workup-20260514/derived/directionality_phase_audit_20260518_003/DIRECTIONALITY_PHASE_AUDIT_NOTES.md` (mtime 2026-05-18).
- Opus v1 visual "winner" style was replaced because its evidence substrate collapsed contributor mixtures and used the wrong sample-weight intuition. Citation: `LEARN/stage1-topology-workup-20260514/derived-opus/README.md` (mtime 2026-05-21).
- bb_O is a phase dissenter in the latest total_iso Opus full-atom chart: negative in alpha-R and left-alpha while other atoms often show positive left-alpha. Citation: `LEARN/stage1-topology-workup-20260514/derived-opus/derived/angle_contributor_chart_v2_full_atoms_20260522_009/bb_O__total_iso/angle_contributor_v2_chart_manifest.csv` (mtime 2026-05-22).

## Open / Contested Threads

These should remain open rather than prematurely resolved:

1. McConnell status:
   - Stage-1 mutation-delta/static evidence says bond anisotropy can be present.
   - May 17 backbone multi-source drop-one says it is small in that conditional block.
   - Current rediscover standalone valid-source McConnell fails.
   - The unresolved question is whether McConnell is a joint/geometry-restricted contributor, a mutation-delta-specific effect, or still partially contaminated by source-definition choices in some views.

2. Field/EFG status:
   - Stage-1 EFG family has real visibility in tensor/source-family contexts.
   - APBS-specific backbone multi-source signal is tiny.
   - Protein-blocked CV can kill apparent EFG rows.
   - Current rediscover local-frame EFG mostly fails.
   - The unresolved question is formulation and confounding: frame, target, atom stratum, APBS vs Coulomb/MOPAC, and field-vs-EFG distinction.

3. Ring law vs ring ecosystem:
   - Current rediscover has a clean Pople ring law.
   - Stage-1 ring family includes ring current, susceptibility, dispersion, pi-quadrupole, tensor specificity, raw magnitude, and geometry overlap.
   - The unresolved question is how much of Stage-1 ring visibility is classical current versus adjacent ring physics and geometry.

4. Charge law vs charge embedding:
   - Current rediscover q/r3 is specific and clean in some backbone static-T2 contexts.
   - Stage-1 charge visibility often comes through AIMNet2 embeddings and broad electrostatic blocks.
   - The unresolved question is how to separate simple classical charge kernels from learned electronic proxies in all-atom/all-category data.

5. Dimension story:
   - Physics bridge gives compact element-specific effective dimensions.
   - Named-atom topology workup shows 12-dim improvements and AIMNet2-heavy source loadings.
   - The unresolved question is which dimensionality definition is appropriate for current mechanism recovery.

6. Geometry charts:
   - Later Opus charts show strong atom-specific phase patterns.
   - Their own notes disallow causal/magnitude claims.
   - The unresolved question is how to turn geometry/phase maps into controlled sample-methods without smuggling in forced geometry as chemistry.

7. H-bond:
   - Sparse/weak in many aggregate views.
   - Some weak-stable rows exist.
   - Some nulls are structural all-zero cases.
   - The unresolved question is whether H-bond matters in narrow, geometry-correct subsamples rather than broad all-atom mutation-delta tables.

## Bridge to Current Rediscover Open Problems

### McConnell

Stage-1 suggests:

- Keep McConnell in the joint/ensemble layer, not only standalone.
- Compare classical and MOPAC McConnell variants.
- Add source filters and near-field/self exclusions as first-class provenance.
- Use geometry-selected sample-methods for clean tests, and maths-methods for confounder-controlled all-data extraction.

Stage-1 warns:

- A standalone failure does not erase mutation-delta/static evidence.
- A mutation-delta/static signal does not prove standalone local T2 recovery.
- Drop-one smallness in a full model can mean overlap, not absence.

### Field / EFG / APBS

Stage-1 suggests:

- Separate APBS E-field/EFG from Coulomb EFG and MOPAC EFG in current outputs.
- Keep tensor component, frame, axis, and target metadata explicit.
- Use blocked/stability and source-overlap diagnostics before interpreting high EFG rows.

Stage-1 warns:

- APBS itself was not load-bearing in the May 17 backbone multi-source block.
- Stage-1 EFG visibility is not a direct fix for current local-frame EFG failure.
- The current issue is likely mechanism/math/formulation, not just APBS radii.

### All-Atom / All-Category Emitter

Stage-1 suggests the emitter should preserve:

- exact source block and source family,
- mechanism label and component kind,
- scalar/vector/tensor decomposition,
- target component and tensor summary,
- atom role, named atom, residue, surface/region,
- normalization axis,
- support counts, nonzero/rank diagnostics,
- exclusion/filter provenance,
- raw and within-protein variants.

Stage-1 warns:

- If these fields are collapsed, later analyses will recreate the same confounds: element pooling, source overlap, high-R2 over-read, and structural null ambiguity.

### Geometric Stratification

Stage-1 suggests:

- Choose scenarios by geometry and support, not IUPAC name alone.
- Use local frames, torsion-axis placement, radial shells, angular bins, and scaffold/forced-geometry flags.
- Treat negative phase as directional signal when appropriate.

Stage-1 warns:

- Geometry can be forced by backbone topology; mark this instead of interpreting it as independent chemistry.
- Smooth surfaces and chart bins require evidence-weighting and low-n display.

### Maths-Methods vs Sample-Methods

Stage-1 supports both:

- Maths-methods: ridge, partial/drop-one, source overlap, blocked CV, tensor specificity, normalization contrast, dimension source mixtures.
- Sample-methods: clean geometry-selected regions, supported local frames, torsion/placement regions, source-filtered neighborhoods.

Stage-1 warns:

- Neither mode replaces the other. Clean samples can become too narrow or support-limited; all-data maths can become proxy-driven or confounded.

### Three-Tool Stool: ORCA / MOPAC / APBS

Stage-1 suggests:

- MOPAC is useful as electronic/EFG/bond-order contrast, especially where classical sources overlap.
- APBS should be kept, but not assumed load-bearing.
- ORCA/DFT targets remain the anchor for mechanism recovery.

Stage-1 warns:

- APBS, MOPAC, Coulomb, and AIMNet2 are not interchangeable "electrostatics." They have different stability and proxy behavior.
- MOPAC EFG and Coulomb EFG can be visible where APBS EFG is tiny.

## What I Could Not Reach or Did Not Read

Read limits were intentional because the spinner is large and the user requested a strict recent-only top filter.

Not read:

- Any file with `mtime < 2026-05-12`.
- The full 42+ GB `derived/` tree and 1.8+ GB `derived-opus/` tree.
- Large binary matrices, NPY payloads, PDFs, PNGs, and most generated figures.
- The IUPAC/topology investigation or matching-policy details, beyond noting that matching policy appears in recent provenance.
- Every individual per-atom/per-residue CSV. I read summary tables and targeted rows for the requested mechanisms.
- `pooled_control_v2.pdf`; I found it but did not open it because it is a PDF/image artifact and there were no adjacent recent text notes in the targeted pass.

Reached recent clusters:

- `LEARN/stage1-topology-workup-20260514/docs/` (mtime 2026-05-15 to 2026-05-20).
- `LEARN/stage1-topology-workup-20260514/derived/tables/` summary/core tables (mtime 2026-05-14 to 2026-05-15).
- `LEARN/stage1-topology-workup-20260514/derived/secondary/tables/blocked_stability_full25_rollup_20260516_clean001/` readouts (mtime 2026-05-16).
- `LEARN/stage1-topology-workup-20260514/derived/signal_visibility_results_v1_20260517_008/` summary packet (mtime 2026-05-17).
- `LEARN/stage1-topology-workup-20260514/derived/backbone_atom_multisource_blocks_20260517_003/` and audits/diagnostics (mtime 2026-05-18).
- `LEARN/stage1-topology-workup-20260514/derived/*geometry*`, `*dihedral*`, `*torsion*`, and `*directionality*` notes (mtime 2026-05-18 to 2026-05-20).
- `LEARN/stage1-topology-workup-20260514/derived-opus/` README/design/current v2 chart notes and manifests (mtime 2026-05-21 to 2026-05-22).
- `LEARN/stage1-mutations/notes/physics_analysis_bridge.md` (mtime 2026-05-25).
- `LEARN/src/actual_physics/OVERVIEW.md` (mtime 2026-05-25).
- `LEARN/src/mutation_set/kernels.py` and `scalars.py` (mtime 2026-05-25).
- `LEARN/bones/` planning/filter notes relevant to geometry and source filtering (mtime 2026-05-25).

## Bottom Line, Without Collapsing the Landscape

The recent Stage-1 mutant work is a map of mechanism visibility under changing stratifications and increasingly strict audits. It does not license a single winner. The most useful transfer to current rediscover is methodological:

- preserve exact provenance,
- separate learned ceilings from physical laws,
- use geometry as a sample axis,
- use joint/partial maths for overlapping mechanisms,
- read APBS/EFG/MOPAC/Coulomb as distinct tools,
- keep tensor/frame/normalization/support metadata attached to every result,
- report unresolved tensions explicitly.

The plurality is the result.
