# All-Atom Charge-Complete Fit Spec - 2026-06-03

Status: **landed and superseded in part** (trued 2026-06-04): Piece 1 emitted charge scalars at `4bb9a0198594773c016160270b6a872b2f9b23fb`; Piece 2 fit/partition landed at `5b5525b45f363e7c4a05cfff6a02a19d2d8c15ab` and the alpha fix at `44e786c094c4f9cd1ae4ebe451454fc3e31406d6`. Build3/Build4/Stage2 and the maths audit now own the current fit story; use `STATE.md`/`NOW.md` before quoting priority or between-axis claims from this spec.
This document specifies one small emit
extension plus the first charge-complete all-atoms joint fit and
partition-by-condition analysis. It extends
`PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md` and targets the v1
`per_atom_substrate` shape committed at `00ec16876d530bfa9ed5a58be83976559eef9dd8`.

Reference v1 run:

```text
/tmp/rediscover-runs/2026-06-03-per-atom-substrate-v1-fixed
558,360 rows = 846 atoms x 660 DFT frames
row_id == frame_slot * n_atoms + atom_index
```

This spec is not a request to build or run anything. The later implementation
has two internal checkpoints, and the second checkpoint must not consume output
until the first one passes.

## Purpose And Scope

The chunk has two pieces:

1. **Small emit extension, C++ spine.** Add two raw charge source scalars to
   `per_atom_substrate`: `ff14sb_charge` and
   `mopac_welford_mean_charge`, with present flags and metadata. This is only
   an additive source-coverage extension. It must not reshape, reorder, or
   reinterpret the existing target, mechanism, conditioning, embedding, audit,
   or named-query outputs.

2. **All-atoms joint fit and partition, Python edge.** On the charge-complete
   emitted substrate, fit DFT T2 in the frozen 2e basis using held-out
   block-CV, then slice held-out predictions by emitted input-side conditions.
   Python reads the emitted rows, sidecars, and column specs only. It does not
   open `trajectory.h5`, ORCA outputs, older run directories, or pair/source
   dumps.

The fit reads one substrate directory. The old `/tmp/combined-mopac-layer3`
multi-directory merge is the pattern for feature tiering, target conversion,
LOAO and frame-block scoring, CIs, and audit outputs, but the new fit must not
repeat its join strategy. The new join is row-aligned inside the emitted
`per_atom_substrate`.

## Gate Sequence

Run gates in this order:

```text
Piece 1 emit-extension gates
  -> charge-complete substrate is accepted
  -> Piece 2 fit and partition gates
```

If any Piece 1 gate fails, the fit is blocked. The analysis must report that it
has no charge-complete substrate rather than silently falling back to v1 or
joining external charge files.

## Piece 1: Small Emit Extension

### Required Columns

Add four row-aligned columns to `per_atom_substrate_rows.csv`, appended after
the existing v1 columns:

```text
ff14sb_charge
ff14sb_charge_present
mopac_welford_mean_charge
mopac_welford_mean_charge_present
```

Appending preserves the ordinal meaning of every existing v1 column. Existing
sidecars and query outputs keep their current shapes.

The scalar values are repeated on every DFT frame row for the same atom because
both sources are atom-axis/static in the current catalog. A missing value is
encoded as `NaN` with its present flag set to `0`.

### Source Contract

| emitted column | C++ source | resident backing source | unit | frame dependence | present flag |
|---|---|---|---|---|---|
| `ff14sb_charge` | `Catalog::value(ArrayId::Ff14sbCharge, atom, row)` | `QtAtom::partialCharge` when `QtAtom::hasPartialCharge` | `e` | static atom axis | `Catalog::present(Ff14sbCharge, atom, row)` and finite value |
| `mopac_welford_mean_charge` | `Catalog::value(ArrayId::MopacChargeWelfordMean, atom, row)` | `h5->mopacChargeWelford()->value[atom].mean` | `e` | static atom axis | `Catalog::present(MopacChargeWelfordMean, atom, row)` and finite value |

Important boundary: do not use `ArrayId::MopacCharge`. The catalog marks it
absent because no per-frame MOPAC charge trajectory exists. The only MOPAC raw
charge source for this chunk is the atom-axis Welford mean.

### Metadata

`Catalog.{h,cpp}` already has the needed source ArraySpecs:

| ArrayId | catalog name | rank | axes | residence | unit | available when |
|---|---|---|---|---|---|---|
| `Ff14sbCharge` | `ff14sb_charge` | scalar | atom, no frame | `StaticTopol` | `e` | topology atoms have partial charges |
| `MopacChargeWelfordMean` | `mopac_charge_welford_mean` | scalar | atom, no frame | `StaticTopol` | `e` | H5 has MOPAC charge Welford moments and the atom has `n_frames_per_atom > 0` |

The emit must add column-spec metadata for the new row columns. Recommended
entries in `per_atom_substrate_column_specs.json`:

| array | column | units | irreps | mechanism | is_feature | native_axis |
|---|---|---|---|---|---:|---|
| `per_atom_substrate_rows` | `ff14sb_charge` | `e` | `0e` | `charges` | true | `rediscover_target_row` |
| `per_atom_substrate_rows` | `ff14sb_charge_present` | | `0e` | `provenance_qc` | false | `rediscover_target_row` |
| `per_atom_substrate_rows` | `mopac_welford_mean_charge` | `e` | `0e` | `charges` | true | `rediscover_target_row` |
| `per_atom_substrate_rows` | `mopac_welford_mean_charge_present` | | `0e` | `provenance_qc` | false | `rediscover_target_row` |

No new physics feature sidecar is required for this chunk. If the Python loader
later needs an explicit catalog entry for row CSV metadata, add it as a loader
description of `per_atom_substrate_rows`, not as a replacement for these
columns.

### Emit-Extension Checks

The emit-extension checkpoint passes only if all checks below pass.

**R and uniqueness**

- `rows == n_atoms * n_dft_frames`.
- Expected 1P9J values remain `846`, `660`, and `558,360`.
- `row_id` is dense and equals `frame_slot * n_atoms + atom_index`.
- `(atom_index, original_frame_index)` is unique.
- Every NPY sidecar first dimension equals the row count.

**Charge scalar coverage**

- New columns exist in `per_atom_substrate_rows.csv`.
- Present flags are `0` or `1`.
- Present flag `1` implies finite scalar value.
- Present flag `0` implies `NaN` scalar value.
- For each atom, `ff14sb_charge` is constant over the 660 frame rows where
  present.
- For each atom, `mopac_welford_mean_charge` is constant over the 660 frame
  rows where present.
- The manifest support block reports
  `ff14sb_charge_present_rows`,
  `mopac_welford_mean_charge_present_rows`, and
  `charge_complete_rows`, where charge-complete means finite
  `ff14sb_charge`, `mopac_welford_mean_charge`, and `aimnet2_charge`.

**Oracle parity unchanged**

- All pre-existing NPY sidecars from v1 are byte-identical against the v1 run
  when emitted from the same input and settings:
  `per_atom_substrate_target_T2.npy`,
  `per_atom_substrate_target_T0.npy`,
  `per_atom_substrate_features_classical.npy`,
  `per_atom_substrate_features_conditioning.npy`,
  `per_atom_substrate_driver_modulation_by_atom.npy`,
  `per_atom_substrate_backbone_audit.npy`, and
  `per_atom_substrate_aimnet2_embedding.npy`.
- `per_atom_substrate_rows.csv` is byte-identical to v1 after projecting away
  the four appended charge columns.
- `per_atom_substrate_ring_identity.csv` and default `query_results/` products
  are unchanged unless review explicitly moves charge-source diagnostics into a
  named query.
- DFT target parity stays exact under the existing frame-alignment diagnostic.

**Backbone regression unchanged**

- Selecting the six backbone strata (`N`, `CA`, `C`, `O`, `HN`, `HA`) from the
  new substrate and ignoring the four new charge columns reproduces the
  existing broad-backbone compatibility checks.
- `per_atom_substrate_backbone_audit.npy` continues to match the broad-backbone
  local charge, field, and McConnell audit payloads for the same
  `(atom_index, original_frame_index)` rows and cutoffs.
- The previous backbone combined-fit pattern must not change its held-out
  scores beyond ordinary floating-point tolerance when it is run on the
  unchanged existing feature surfaces. The new scalar columns are not allowed to
  perturb old reducers.

**Charge positive-control readiness**

- The charge source agreement edge is directly evaluable from one emitted
  substrate:
  `formal_charge`, `ff14sb_charge`, `mopac_welford_mean_charge`, and
  `aimnet2_charge`.
- No H5 read, external charge table, or multi-directory join is needed for that
  positive control.

## Piece 2: All-Atoms Joint Fit

### Inputs

The analysis input is the accepted charge-complete substrate directory. Required
files:

```text
per_atom_substrate_rows.csv
per_atom_substrate_manifest.json
per_atom_substrate_column_specs.json
per_atom_substrate_target_T2.npy
per_atom_substrate_features_classical.npy
per_atom_substrate_features_conditioning.npy
per_atom_substrate_driver_modulation_by_atom.npy
per_atom_substrate_aimnet2_embedding.npy
```

`per_atom_substrate_target_T0.npy`, `per_atom_substrate_backbone_audit.npy`,
`per_atom_substrate_ring_identity.csv`, and named query outputs may be read for
audit or display, but they are not required for the T2 fit.

### Join Contract

The fit joins only within the emitted substrate:

- CSV row `i` joins to sidecar row `i`.
- `driver_modulation_by_atom` joins by `atom_index`.
- `aimnet2_embedding` joins by row id.
- Feature names and blocks come from `per_atom_substrate_column_specs.json`.

No join to `/tmp/combined-mopac-layer3`, broad-backbone directories, MOPAC
directories, all-atom-equivariant directories, `trajectory.h5`, ORCA outputs,
or raw pair/source rows is permitted.

### Target

The target is DFT total T2:

```text
y_library = per_atom_substrate_target_T2.npy[:, 0:5]
y_2e = y_library @ change_of_basis.get_C().T
```

The C++ library order is `[xy, yz, zz, xz, xx-yy]`. The fitter reports R2 on
the 5-component 2e vector after the frozen `get_C()` conversion. It does not
score only `|T2|`.

Every T2 feature block emitted in the same library T2 basis is transformed with
the same `get_C()` matrix before fitting. Scalar and vector features are not
projected in Python.

### Feature Tiers

Primary model is a global all-atoms multi-output ridge fit. Scores are reported
globally and sliced by emitted atom type. The global fit is what makes the
between-axis use all 846 atoms instead of the roughly 54 atoms available in each
old backbone stratum.

Primary atom type for score tables: emitted `stratum`. Also report `role` and
`ff_atom_type_ord` coverage so review can decide whether a later table should
use a deeper identity rollup.

Predeclared tiers:

| tier | included feature blocks |
|---|---|
| `classical_mechanisms_combined` | `ring_jb_T0/T2`, `charge_q_over_r3_T2`, `mc_lit_T0/T2_valid`, `mopac_coulomb_shielding_T2`, `mopac_mc_shielding_T2`, `ff14sb_field_x/y/z/mag`, `apbs_E_x/y/z/mag`, `apbs_efg_T2` |
| `plus_AIMNet2` | classical tier plus `aimnet2_charge`, `aimnet2_crg_scalar`, `aimnet2_crg_x/y/z`, and `per_atom_substrate_aimnet2_embedding.npy` |
| `all` | `plus_AIMNet2` plus `formal_charge`, `ff14sb_charge`, and `mopac_welford_mean_charge` |

Conditioning columns are slicers by default, not predictive covariates. This
keeps the response curves interpretable. A later sensitivity may add an
`all_plus_conditioners` tier, but that is outside this chunk unless review
explicitly pulls it in.

Embedding default: use the full 256-dimensional f32 sidecar under ridge for the
global all-atoms fit. If a partition-specific refit is added later and becomes
feature-wide relative to atom count, use unsupervised PCA fit on training rows
only and record the transform in the audit.

### Ridge And Preprocessing

- Multi-output ridge predicts all five T2 components together.
- Add an intercept.
- Fit imputation means and feature standardization on training rows only.
- Fill missing feature values with the training mean for that feature.
- Default ridge alpha is `10.0` for continuity with
  `/tmp/combined-mopac-layer3`.
- If alpha selection is changed, choose it by training-only inner CV from a
  predeclared grid and report the grid. The held-out test fold must not choose
  alpha.

### Held-Out Axes

Report held-out test R2, never in-sample R2.

**Between atoms: leave-atoms-out**

- Unit of independence is `atom_index`.
- For each held-out atom, train on all other atoms.
- The between score uses atom means across frames: one `xbar` and one `ybar`
  per atom, with held-out predictions for held-out atom means.
- Primary global between score uses all 846 atoms.
- Per-atom-type table rows slice held-out atom predictions by `stratum` and
  report the number of atoms in that slice.
- CI is jackknife over atoms.

**Within atoms: held-out frame block**

- Center `x` and `y` within each atom using training frames only.
- Hold out a contiguous block of original frames, default 20 percent.
- Purge at least one adjacent frame on both sides of the held-out block.
- Train on centered training rows and score centered held-out rows.
- Report cross-split lag-1 adjacency count; expected target is zero after purge.
- CI is jackknife over atoms.
- Report AR(1)-aware effective N by T2 component:

```text
N_eff_atom_component = n * (1 - rho1) / (1 + rho1)
```

Clamp `rho1` to the open interval `(-0.999, 0.999)` and cap `N_eff` to
`[1, n]`. Sum over atoms for each component, then report min, median, and max
over the five components plus median lag-1 rho.

### Score Outputs

Write these artifacts:

```text
allatom_fit_score_table.csv
allatom_fit_score_long.csv
allatom_fit_report.md
join_coverage.csv
feature_blocks_used.json
run_audit.json
```

Required columns in `allatom_fit_score_table.csv`:

```text
target
atom_type_axis
atom_type
tier
tier_label
rows
n_atoms_between_global
n_atoms_in_slice
n_original_frames
n_features
ridge_alpha
variance_share_between
variance_share_within
within_N_eff_min
within_N_eff_median
within_N_eff_max
median_lag1_rho
thin_flag
p_gt_atoms_flag
between_LOAO_test_R2
between_LOAO_test_R2_jackknife_se
between_LOAO_test_R2_ci95_low
between_LOAO_test_R2_ci95_high
within_frameblock_test_R2
within_frameblock_test_R2_jackknife_se
within_frameblock_test_R2_ci95_low
within_frameblock_test_R2_ci95_high
within_split_strategy
test_frames
purged_train_frames
cross_split_lag1_pairs
feature_blocks
```

`thin_flag` is set when a reported atom-type slice has fewer than 10 atoms or
when any tier has too few held-out atoms for a stable jackknife. `p_gt_atoms`
is set when the feature count is greater than or equal to the number of atoms
used for that score slice.

## Partition By Condition

Partitioning uses emitted input-side columns only. It may derive algebraic
conditioners from emitted input-side columns, but every derived formula must be
listed in `run_audit.json`. Partitioning must not use DFT target values,
held-out residuals, or fitted coefficients to define bins.

### Predeclared Condition Families

| family | emitted inputs |
|---|---|
| atom identity | `stratum`, `role`, `element_ord`, `ff_atom_type_ord`, `formal_charge`, `planar_group_ord`, `polar_h_kind_ord`, `aromatic`, `is_exchangeable`, residue and locant fields |
| geometry and isolation | `nearest_ring_r`, `nearest_charge_r`, `nearest_bond_midpoint_r`, `nearest_heavy_atom_r`, ring/charge/bond counts, self-or-bonded counts |
| driver magnitude | `abs_ring_jb_T2`, `abs_charge_T2`, `abs_mc_lit_T2`, `abs_mopac_coulomb_T2`, `abs_mopac_mc_T2`, `abs_apbs_efg_T2`, `abs_apbs_E`, `abs_ff14sb_E`, `abs_aimnet2_crg` |
| driver modulation | `sd_*_by_atom` columns from `per_atom_substrate_driver_modulation_by_atom.npy` |
| support and source availability | `ring_present`, `charge_present`, `mc_lit_valid_present`, `ff14sb_field_present`, `apbs_E_present`, `apbs_efg_present`, `mopac_coulomb_shielding_present`, `mopac_mc_shielding_present`, `aimnet2_*_present`, new charge present flags |
| charge-source agreement | emitted `formal_charge`, `ff14sb_charge`, `mopac_welford_mean_charge`, `aimnet2_charge`, plus predeclared absolute differences and sign-agreement indicators |

Continuous conditioners use predeclared bins:

- geometry distances: existing threshold bands where available (`4A`, `6A`,
  `8A`, `10A`) plus quantiles for nearest-distance columns;
- driver magnitudes and modulation: quintiles within the scored atom-type
  slice, with bin edges written to the audit;
- charge-source agreement: quantiles of absolute pairwise charge differences
  and categorical sign agreement/disagreement.

Categorical conditioners use emitted categories directly, with sparse levels
rolled up before scoring. Minimum recommended reporting support is 10 atoms and
500 held-out rows per curve point; thinner points remain in CSV with a thin
flag but are not eligible for the favorable-case shortlist.

### Response Curves

Primary curves slice held-out predictions from the predeclared CV folds. They
do not refit per bin. For each conditioner, atom type, tier, and CV axis,
report:

```text
condition_family
condition_name
bin_label
bin_low
bin_high
atom_type
tier
axis
rows
n_atoms
N_eff_median
heldout_R2
heldout_R2_ci95_low
heldout_R2_ci95_high
delta_R2_vs_classical
delta_R2_vs_previous_tier
thin_flag
```

Curves can rise or fall. The report should read the shape directly: monotone
rise, monotone fall, U-shape, threshold behavior, flat response, or unstable
thin response.

### Emergent Favorable Cases

After all predeclared curves are written, select a small favorable-case table
from the same curve outputs. A case is eligible only if:

- it is from a predeclared condition family and binning rule;
- it has enough atoms and N_eff for the reported axis;
- its held-out R2 or tier delta is positive with a positive CI lower bound, or
  it is one of the top ranked bins by held-out R2 with the uncertainty clearly
  reported;
- it was not selected using DFT residuals, target magnitude, or manual
  post-hoc chemistry labels.

Recommended shortlist size: 10 to 20 rows. Include unfavorable or falling
curves in the response-curve outputs even when they are not shortlisted.

Required artifacts:

```text
partition_response_curves.csv
partition_response_curves_long.csv
partition_favorable_cases.csv
partition_report.md
```

## Fit-Stage Checks

The fit checkpoint passes only if all checks below pass.

**Input acceptance**

- Manifest says `relationship = per_atom_substrate` and
  `relationship_kind = per_atom_aggregate`.
- Shape is `846 x 660` for the reference 1P9J run unless review intentionally
  points at another run.
- New charge columns and flags exist.
- `charge_complete_rows` is reported and matches the filter used by the fit.
- Every selected feature column is found through emitted column specs or row
  schema metadata.

**No external merge**

- The audit states that Python read only the emitted substrate directory.
- The audit states that Python did not open H5, ORCA, older broad-backbone,
  MOPAC, all-atom-equivariant, or source/pair dump directories.
- Sidecar joins are by row id; driver modulation joins by atom id.

**Basis and target**

- `abs(C.T @ C - I).max()` for `change_of_basis.get_C()` is below `1e-12`.
- Target sidecar shape is `(R, 5)`.
- All T2 blocks used in the fit have five components and are transformed
  consistently.

**CV integrity**

- Between-axis held-out atoms are never in the training set.
- Within-axis held-out frame block is not used for preprocessing, imputation,
  standardization, alpha selection, or atom-centering means.
- Purged adjacent frames are absent from train and test.
- Reported scores are held-out R2.
- `N_eff` and lag-1 rho are reported for every atom-type/tier row.

**Partition integrity**

- Every condition column is emitted or algebraically derived from emitted
  input-side columns.
- Bin definitions are written before ranking favorable cases.
- No condition uses DFT target, residual, or fitted coefficient values.
- Thin bins are flagged and excluded from the favorable shortlist.

## Expected Relationships To Assess

These are expectations to test, not acceptance criteria.

| relationship | expected direction | probability before fit |
|---|---|---:|
| All-atoms global between-axis score is more stable than the old per-backbone-stratum between score | narrower CIs because the held-out atom axis has up to 846 atoms | 0.70 |
| `plus_AIMNet2` improves within-axis recovery in chemically specific and high-modulation rows | positive delta most visible in charge-response and embedding-sensitive partitions | 0.65 |
| `all` improves charge-sensitive partitions after adding raw FF14SB and MOPAC Welford charges | strongest in charged heavy, polar H, carbonyl/amide sidechains, and high charge-disagreement bins | 0.55 |
| APBS EFG and MOPAC Coulomb shielding remain moderately aligned after sign/scale differences | charge/electrostatic partitions show better recovery than pooled all-atoms | 0.65 |
| Ring Johnson-Bovey aggregate recovers best near aromatic H and near-ring geometry bins | response rises as nearest ring distance falls and ring support/dominance rises | 0.60 |
| McConnell aggregate recovers best with valid bond support and low self/bonded contamination | response improves in HN/peptide-plane and high valid-bond-support bins | 0.55 |
| FF14SB, MOPAC Welford, and AIMNet2 raw charges agree unevenly | agreement is better within neutral/aliphatic contexts and worse in polar/charged contexts | 0.55 |
| APBS field and FF14SB vacuum field remain divergent in screened or polar environments | field response curves may fall or change sign in polar/charged bins | 0.60 |

Uncertainty is highest for thin atom types such as sulfur and narrow deep
identity partitions. Those should be reported with atom counts and N_eff rather
than promoted by point estimate alone.

## Open Questions For Review

1. **Row CSV versus sidecar for raw charge scalars.** This spec chooses appended
   row columns because the values are scalar, static, and useful for filters.
   Review can choose a tiny row-aligned sidecar instead if loader ergonomics
   demand it, but existing v1 sidecars should remain unchanged.

2. **Python catalog representation for row columns.** `_catalog.py` currently
   describes NPY sidecars, while per-column metadata lives in
   `per_atom_substrate_column_specs.json`. Decide whether row CSV columns need a
   formal Python `ArraySpec` entry or whether column specs are the source of
   truth.

3. **Tier semantics.** This spec keeps conditioning variables out of the three
   primary tiers. Decide whether a fourth `all_plus_conditioners` sensitivity is
   worth adding later.

4. **Embedding handling.** Default is full 256-d embedding under ridge for the
   global all-atoms fit. Review whether to require training-only PCA for the
   primary run for continuity with `/tmp/combined-mopac-layer3`.

5. **Atom-type axis.** Primary table uses emitted `stratum`. Review whether the
   first report should also include `role`, `ff_atom_type_ord`, or a deep
   identity rollup as first-class score axes.

6. **Favorable-case threshold.** The shortlist rule uses held-out R2 or tier
   delta with uncertainty and minimum support. Review the exact minimum atom
   count, row count, and N_eff thresholds before execution.

7. **Charge-complete strictness.** This spec blocks the fit if the charge source
   columns are missing and reports charge-complete rows explicitly. Review
   whether partial charge coverage should be allowed for non-charge tiers or
   whether all tiers should use the same charge-complete row set for clean
   comparisons.
