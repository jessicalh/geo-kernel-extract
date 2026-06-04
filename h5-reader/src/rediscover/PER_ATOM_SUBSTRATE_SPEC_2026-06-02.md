# Per-Atom Aggregate Substrate Spec - 2026-06-02

Status: **landed as v1 and extended**; trued 2026-06-04. v1 substrate landed at `00ec16876d530bfa9ed5a58be83976559eef9dd8`, later extended through Build1, Build4, Stage2, and the maths-audit fixes. Use `STATE.md`/`NOW.md` for current limitations: 1P9J between/LOAO positives are retracted, 720-WT owns transfer, and clean verdict re-runs are parked.
This document specifies the lean #58
node-store realization that an 846-atom joint fit plus partition-by-condition
analysis consumes.

## Purpose And Scope

The substrate is one row per `(atom, DFT-present frame)` carrying identity,
DFT target, mechanism aggregates, and input-side conditioning variables. It is
the small per-atom aggregate face of the node-store contract in
`NODE_STORE_CONTRACT_2026-06-02.md`: C++ owns association, source selection,
kernel construction, and reduction; Python reads named aggregate arrays for
fitting, CV, plots, and the frozen `get_C()` basis conversion.

The substrate is explicitly **not**:

- a per-pair/source dump;
- a Python re-sum of source rows;
- a second Python physics model;
- a replacement for the transient pair-index query in the parent contract.

It enables a joint all-atoms design matrix with about 558,360 rows for 1P9J
(`846 atoms * 660 DFT-labeled frames`), plus conditioning partitions that slice
only on inputs: geometry, source dominance, neighbor support, driver modulation,
and self/bonded status. The current `/tmp/combined-mopac-layer3` run is the
immediate consumer prototype: it joined emitted rows by `(atom_index,
original_index)`, used the frozen `get_C()` conversion, and combined ring,
charge, McConnell, Buckingham/APBS, MOPAC, and AIMNet2 blocks. The #58 substrate
removes the multi-directory merge and source-row reassembly from that workflow.

## Primary Row Axis

Canonical row key:

```text
row_id = dense row ordinal in emission order
atom_index = QtProtein atom row
frame_slot = ordinal among run.frameMap.dftRows()
h5_row = trajectory H5 frame row
original_frame_index = frameMap.originalIndex(h5_row)
```

Emission order should match `RelationshipEngine::RunTraversal`: DFT row outer,
atom inner. For a dense all-atom run, `row_id == frame_slot * n_atoms +
atom_index` must hold. Consumers must still join by `(atom_index,
original_frame_index)` when comparing to older substrates, because the
MOPAC-layer3 audit showed `h5_row` can differ across emitted runs while the
physical frame key is stable.

Primary tensor frame:

- Default #58 all-atom aggregate tensors are in the H5/ORCA-aligned molecular
  frame used by `AllAtomEquivariant`. The manifest must retain the existing DFT
  frame-alignment diagnostic from `OutputManifest`.
- Backbone-local sidecars may be emitted as compatibility/audit payloads, but
  they are not the primary all-atom joint-fit target.
- T2 storage follows the existing C++ library isometric ordering
  `[xy, yz, zz, xz, xx-yy]`; the analysis basis is obtained by multiplying by
  the frozen `change_of_basis.get_C()` matrix. ArraySpec metadata still declares
  these 5 columns as `1x2e`.

## Files And Row Alignment

Preferred materialization under a managed `tools/rediscover_run.py` run:

```text
per_atom_substrate_rows.csv
per_atom_substrate_manifest.json
per_atom_substrate_column_specs.json
per_atom_substrate_target_T2.npy
per_atom_substrate_target_T0.npy                  optional
per_atom_substrate_features_classical.npy
per_atom_substrate_features_conditioning.npy
per_atom_substrate_aimnet2_embedding.npy          optional/separable
```

Row contract:

- CSV row `i` has `row_id == i`.
- Every NPY sidecar has first dimension `R == len(per_atom_substrate_rows.csv)`.
- Sidecar row `i` belongs to CSV row `i`.
- Missing unavailable mechanisms are encoded as `present=0` in CSV plus `NaN`
  payloads in sidecars, not as dropped rows.
- `per_atom_substrate_aimnet2_embedding.npy` is separable by default. Its
  absence must not invalidate the classical substrate.

Equivalent node-store form is allowed if it preserves the same metadata:

```text
atoms/frame_slots/...
arrays shaped (n_frames, n_atoms, components)
row_id(frame_slot, atom_index) = frame_slot * n_atoms + atom_index
```

The CSV-plus-NPY form is preferred for the first #58 delivery because it matches
the existing rediscover manifest and ArraySpec discipline.

## Schema: Identity And Target

Identity columns live in `per_atom_substrate_rows.csv`.

| column | dtype | unit | irreps | mechanism | sign | source |
|---|---:|---|---|---|---|---|
| row_id | int64 | | | topology | | new ordinal |
| atom_index | int32 | | | topology | | `QtAtom::atomIndex` |
| h5_row | int32 | | | topology | | `RunTraversal` row |
| frame_slot | int32 | | | topology | | ordinal of DFT row |
| original_frame_index | int32 | | | topology | | `frameMap.originalIndex` |
| time_ps | float64 | ps | 0e | topology | | trajectory time |
| element_ord | int16 | | 0e categorical | topology | | `QtAtom::element` |
| ff_atom_type_ord | int16 | | 0e categorical | topology | | `QtAtom::ffAtomType` |
| atom_name | string | | | topology | | BMRB display projection |
| residue_index | int32 | | | topology | | `QtAtom::residueIndex` |
| residue_number | int32 | | | topology | | `QtResidue::address` |
| amino_acid_ord | int16 | | 0e categorical | topology | | `QtResidue::aminoAcid` |
| backbone_role_ord | int8 | | 0e categorical | topology | | `QtAtom::backboneRole` |
| locant_ord | int8 | | 0e categorical | topology | | `QtAtom::locant` |
| planar_group_ord | int8 | | 0e categorical | topology | | `QtAtom::planarGroup` |
| polar_h_kind_ord | int8 | | 0e categorical | topology | | `QtAtom::polarH` |
| ring_position_primary_ord | int8 | | 0e categorical | topology | | `QtAtom::ringPositionPrimary` |
| role | string/int16 | | 0e categorical | topology | | derived typed role |
| stratum | string/int16 | | 0e categorical | topology | | partition label |

Proposed `role` labels:

```text
backbone_N, backbone_CA, backbone_C, backbone_O, backbone_HN,
backbone_HA, aromatic_H, polar_sidechain_H, aliphatic_H,
aromatic_heavy, sidechain_heavy, terminal_or_cap, other
```

Proposed `stratum` labels start with the previous six backbone strata
`N, CA, C, O, HN, HA` and add all-atom partitions:

```text
aromatic_H, polar_H, aliphatic_H, aromatic_heavy, charged_heavy,
carbonyl_sidechain, amide_sidechain, sulfur, other_heavy, other
```

Target payloads:

| column or sidecar | dtype | shape | unit | irreps | mechanism | sign |
|---|---:|---:|---|---|---|---|
| dft_present | int8 | CSV | | | quantum_reference | |
| dft_T2_0..4 | float64 | `(R,5)` | ppm | 1x2e | quantum_reference | shielding convention from `_catalog.py`: `sigma_ab = -dB_sec_a/dB0_b` |
| dft_T0 | float64 | `(R,)` optional | ppm | 1x0e | quantum_reference | same |
| dft_raw_3x3 | float64 | optional `(R,3,3)` | ppm | 0e+1e+2e | quantum_reference | same |

The required target for the joint fit is `dft_T2_0..4`. `dft_T0` is optional
and must not be required for the T2 pipeline.

## Schema: Mechanism Aggregate Features

Classical feature sidecars should use grouped arrays plus explicit column specs
instead of hundreds of CSV float columns. CSV carries present flags, support
counts, and scalar conditioning values.

### Required aggregate feature groups

| group | columns | unit | irreps | mechanism tag | sign convention | status |
|---|---|---|---|---|---|---|
| ring_jb | `ring_jb_T0`, `ring_jb_T2_0..4`, `ring_jb_present`, `ring_n`, `ring_valid_n` | ppm | 1x0e + 1x2e | ring_current | Johnson-Bovey shielding, same shielding sign as target | **Partial**: source-level `jb_T0`/`jb_T2_local_*` are already computed by `RingCurrentNeighborhood` and `ComposedRelationships::ringAttacher`; aggregate sum is not yet emitted. Add a C++ fold over attached ring slots. |
| charge_q_over_r3 | `charge_q_over_r3_T2_0..4`, `charge_present`, `charge_n`, `charge_cutoff_A` | current broad uses CoulombKe-scaled charge EFG-like T2; document as V/A^2-equivalent until final unit name is blessed | 1x2e | electrostatic_efg | EFG sign from `q*(3 rr/r5 - I/r3)` with target-to-source displacement; no shielding sign unless converted | **Exists for backbone** as `BroadBackbone::chargeKernelT2FromSources` and `broad_backbone_aggregated_charge_literature_kernel_T2.npy`; must be generalized to all atoms. |
| mc_lit_valid | `mc_lit_T0_valid`, `mc_lit_T2_valid_0..4`, `mc_lit_valid_present`, `bond_n_valid`, `mc_near_field_ratio` | ppm | 1x0e + 1x2e | bond_anisotropy | McConnell literature-scaled shielding, target shielding sign | **Exists for backbone** in `BroadBackboneSink` valid-source columns; source kernel from `McConnellSourceLiteratureKernelLocal`. Needs all-atom frame decision and generalized carrier. |
| mopac_coulomb_shielding | `mopac_coulomb_shielding_T2_0..4`, `mopac_coulomb_shielding_present` | ppm | 1x2e | electrostatic_efg | already a shielding T2; not raw EFG | **Exists all-atom** in `AllAtomEquivariantSink` and `Catalog::MopacCoulombShielding`. |
| mopac_mc_shielding | `mopac_mc_shielding_T2_0..4`, `mopac_mc_shielding_present` | ppm | 1x2e | bond_anisotropy | shielding sign | **Exists all-atom** in `AllAtomEquivariantSink` and `Catalog::MopacMcShielding`. |
| ff14sb_field | `ff14sb_field_x..z`, `ff14sb_field_mag`, `ff14sb_field_present`, `charge_n` | e/A^2 or CoulombKe-scaled field; choose one and record | 1x1o | charges | field sign `E = sum q_i * (r_atom - r_i)/r^3`; existing broad local field uses `-q*disp/r^3` where `disp = r_i-r_atom` | **Exists for backbone** in `ReduceBroadBackboneSources`; must be generalized to all atoms. |
| buckingham_apbs_field | `apbs_E_x..z`, `apbs_E_mag`, `apbs_E_present` | V/A | 1x1o | electrostatic_efg | APBS electric field as loaded from H5 | **Exists all-atom** in `AllAtomEquivariantSink` and backbone-local in `BuckinghamEfieldSink`. |
| apbs_efg | `apbs_efg_T2_0..4`, `apbs_efg_present` | V/A^2 | 1x2e | electrostatic_efg | APBS EFG T2, symmetric traceless | **Exists all-atom** in `AllAtomEquivariantSink`; local backbone version in `EfgFeatureSink`. |
| aimnet2_charge_crg | `aimnet2_charge`, `aimnet2_crg_scalar`, `aimnet2_crg_x..z`, present flags | e, e^2/A | 1x0e + 1x0e + 1x1o | aimnet2 | CRG is labelled charge-response-gradient, not polarizability | **Exists all-atom** in `AllAtomEquivariantSink`; backbone-local in `Aimnet2FeatureSink`. |
| aimnet2_embedding | `embedding_000..255` sidecar | dimensionless | 256x0e | aimnet2 | separable learned embedding | **Exists all-atom** as `all_atom_equivariant_aimnet2_embedding.npy`; store as optional f32 sidecar. |

### Recommended compact feature arrays

`per_atom_substrate_features_classical.npy` should be row-major float32 or
float64 with column metadata in `per_atom_substrate_column_specs.json`. Initial
column order:

```text
ring_jb_T0
ring_jb_T2_0..4
charge_q_over_r3_T2_0..4
mc_lit_T0_valid
mc_lit_T2_valid_0..4
mopac_coulomb_shielding_T2_0..4
mopac_mc_shielding_T2_0..4
ff14sb_field_x, ff14sb_field_y, ff14sb_field_z, ff14sb_field_mag
apbs_E_x, apbs_E_y, apbs_E_z, apbs_E_mag
apbs_efg_T2_0..4
aimnet2_charge
aimnet2_crg_scalar
aimnet2_crg_x, aimnet2_crg_y, aimnet2_crg_z
```

This gives 47 non-conditioning float feature columns before optional aliases or
support scalars. The previous `/tmp/combined-mopac-layer3` fit used 44
classical columns for the backbone-only design and 81 columns after adding
32 embedding PCs plus AIMNet2 charge/CRG; the #58 all-atom version should keep
the full 256-d embedding in a separate f32 sidecar and let the fitter choose
PCA or direct regularization.

## Schema: Conditioning Variables

Conditioning variables are input-side partition controls. They must be computed
from topology, coordinates, source support, and emitted mechanism features only.
They must never use DFT target values or residuals.

Proposed required conditioning columns:

| group | columns | unit | irreps | mechanism | status |
|---|---|---|---|---|---|
| geometry nearest | `nearest_ring_r`, `nearest_charge_r`, `nearest_bond_midpoint_r`, `nearest_heavy_atom_r` | A | 0e | geometry | Existing KD trees can answer these; add min-distance reducers. |
| neighbor counts | `ring_count_4A`, `ring_count_6A`, `ring_count_8A`, `charge_count_4A`, `charge_count_6A`, `charge_count_10A`, `bond_count_4A`, `bond_count_8A`, `bond_count_10A` | count | 0e | geometry | Existing `SpatialIndexSet::near/range`; add count reducer. |
| dominance fractions | `ring_abs_T2_dom_frac`, `charge_abs_T2_dom_frac`, `mc_abs_T2_dom_frac`, `mopac_abs_T2_dom_frac`, `aimnet2_embedding_norm_rank_frac` | dimensionless | 0e | conditioning | Needs C++ aggregate support for top-source contribution or feature-block magnitude. |
| support validity | `ring_present`, `charge_present`, `mc_lit_valid_present`, `apbs_present`, `mopac_present`, `aimnet2_present` | bool | 0e | provenance_qc | Present flags already exist in current sinks. |
| self/bonded | `ring_self_or_bonded_n`, `bond_self_or_bonded_n`, `charge_excluded_same_residue_n`, `has_self_or_bonded_driver` | count/bool | 0e | topology | Ring and bond flags exist on source slots; aggregate counts need #58 fold. |
| driver magnitude per row | `abs_ring_jb_T2`, `abs_charge_T2`, `abs_mc_lit_T2`, `abs_mopac_coulomb_T2`, `abs_mopac_mc_T2`, `abs_apbs_efg_T2`, `abs_apbs_E`, `abs_ff14sb_E`, `abs_aimnet2_crg` | mechanism units | 0e | conditioning | Norms of emitted input-side feature groups. |
| driver modulation per atom | `sd_ring_jb_T2_by_atom`, `sd_charge_T2_by_atom`, `sd_mc_lit_T2_by_atom`, `sd_mopac_coulomb_T2_by_atom`, `sd_mopac_mc_T2_by_atom`, `sd_apbs_efg_T2_by_atom`, `sd_apbs_E_by_atom`, `sd_ff14sb_E_by_atom`, `sd_aimnet2_crg_by_atom` | mechanism units | 0e | conditioning | C++ or post-emission reducer over input features only; repeat per row or store atom-axis sidecar. |

Open to vet: dominance fractions require either a top-source accumulator inside
the reducer or a named pair-index query. They should not be computed from a
resident pair dump.

## Emission Plan Grounded In Existing Code

The #58 emitter should be a new case, tentatively `--case per_atom_substrate`,
implemented as a target-axis carrier over `RelationshipEngine::RunTraversal`
or the typed record-shape overload already used by `AllAtomEquivariant`.

Existing code to reuse:

- `RelationshipEngine::RunTraversal` is the canonical DFT-row-outer,
  atom-inner walk and already supplies atom, H5 row, original index, target, and
  attached sources.
- `Catalog.{h,cpp}` already exposes positions, APBS E-field/EFG, AIMNet2
  charge/CRG/embedding, MOPAC-Coulomb shielding, MOPAC-Mc shielding, static
  FF14SB charges, and MOPAC Welford mean charge. It also marks per-frame MOPAC
  charges absent.
- `AllAtomEquivariant.{h,cpp}` already defines the all-atom stratum, lab-frame
  target, all-atom source selectors, and all-atom target-axis sidecars. It
  should provide the identity/target/direct-H5 feature pattern, but #58 must
  not keep its source CSV as the main product.
- `BroadBackbone.{h,cpp}` and `BroadBackboneSink` already contain the C++
  reducers for ring/bond/charge source aggregation, FF14SB field, valid-source
  McConnell, and charge EFG-like T2. These reducers are currently backbone-local
  and must be generalized to all atoms or explicitly scoped.
- `RingCurrentNeighborhood` and `ComposedRelationships::ringAttacher` already
  compute source-level Johnson-Bovey `jb_T0` and `jb_T2`. #58 needs a C++ sum
  over those attached source kernels.
- `EfgFeature`, `BuckinghamEfield`, and `Aimnet2Feature` prove the no-source
  per-atom feature carrier pattern through the unified traversal.
- `OutputManifest` already records row-alignment contracts, support counts,
  normalization, and sidecar file lists. #58 should extend that manifest rather
  than invent another output registry.
- The Python ArraySpec discipline in
  `/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py` defines the
  required metadata shape: `native_axis`, `irreps`, `units`,
  `sign_convention`, `tensor_rank`, `parity`, `mechanism`, and `is_feature`.

Reducer status summary:

| mechanism | already emits target-axis aggregate? | required #58 action |
|---|---|---|
| DFT T2 target | yes, all-atom in `AllAtomEquivariantSink` | reuse row-aligned target sidecar and `get_C()` metadata |
| APBS E-field/EFG | yes, all-atom in `AllAtomEquivariantSink` | reuse, normalize metadata irreps/parity |
| AIMNet2 charge/CRG/embedding | yes, all-atom in `AllAtomEquivariantSink` | reuse, keep embedding separable |
| MOPAC Coulomb shielding T2 | yes, all-atom in `AllAtomEquivariantSink` | reuse; label as EFG-derived shielding, not raw EFG |
| MOPAC McConnell shielding T2 | yes, all-atom in `AllAtomEquivariantSink` | reuse |
| FF14SB field | backbone aggregate only in `BroadBackbone` | generalize reducer to all atoms |
| charge q/r3 T2 | backbone aggregate only in `BroadBackbone::chargeKernelT2FromSources` | generalize reducer to all atoms |
| McConnell valid-source literature kernel | backbone aggregate only in `BroadBackboneSink` | generalize, with frame choice vetted |
| ring Johnson-Bovey sum | source-level only in ring-current path | add C++ aggregate fold; no new physics |
| geometry/conditioning dominance | not emitted as target-axis aggregate | add C++ aggregate counters/top-source stats or defer to pair-index query |

## Size Estimate

Rows:

```text
R = 846 atoms * 660 DFT frames = 558,360 rows
```

Binary payload estimates:

| payload | columns | dtype | estimate |
|---|---:|---:|---:|
| target T2 + T0 + direct classical features | about 55 columns | float64 | about 246 MB |
| conditioning features | about 30 columns | float32/float64 mixed | about 70-135 MB |
| identity/present/support columns | about 35 columns | int/float in CSV or compact sidecar | about 40-120 MB depending on CSV text |
| AIMNet2 embedding | 256 columns | float32 | about 572 MB decimal |

Expected lean size:

- Without embedding: about 350-500 MB for binary sidecars plus compact CSV.
- With embedding sidecar: about 0.9 GB binary-equivalent, still under 1 GB if
  embedding is f32 and wide float payloads stay in NPY sidecars rather than CSV.
- A float64 embedding would add about 1.14 GB by itself and is out of scope for
  the lean form.

This is two orders of magnitude smaller than the old resident pair expansion
because source payloads are reduced once per target row and the pair index stays
a pointer/query structure. Any full pair materialization remains a named
transient under `tools/rediscover_run.py` and is deleted by the run framework's
drop-old substrate policy.

## Fit And Partition Use

The all-atom substrate should support:

- a joint T2 fit across all 846 atoms, with `role` and `stratum` as model
  covariates or partition keys;
- leave-atom-out and grouped-by-role tests where the between-axis has many more
  atoms than the six-stratum backbone fit;
- condition partitions such as isolated-ring, charge-dominated, McConnell-valid,
  MOPAC-present, APBS-present, AIMNet2-present, high-driver-modulation, and
  low-self/bonded-support;
- block ablations using the exact feature groups listed above.

Expected relationships to assess:

| feature block | expected strongest partitions | probability before fit |
|---|---|---:|
| MOPAC-Coulomb shielding T2 | heavy atoms and charge-dominated local environments | high, about 0.70 |
| AIMNet2 embedding/CRG | within-atom modulation and chemically specific residual structure | medium-high, about 0.65 |
| ring Johnson-Bovey aggregate | aromatic H, aromatic-heavy neighbors, near-ring partitions | medium-high, about 0.60 |
| McConnell valid-source aggregate | HN, peptide-plane atoms, high valid-bond-support rows | medium, about 0.55 |
| charge q/r3 T2 aggregate | charged/polar neighborhoods and long-range electrostatic partitions | medium, about 0.50 |
| APBS EFG and APBS/Buckingham field | solvent/continuum-sensitive rows; useful mainly as partitioned support | medium-low, about 0.40 |

The determinability question improves materially versus the current backbone
fit: the full all-atom between-axis has up to 846 atoms. The non-embedding
feature count should stay below about 90. The full 256-d embedding makes the
raw feature count roughly 300-350, still below 846 globally but too wide for
some narrow partitions; those partitions should use ridge, PCA, or exclude the
embedding block.

## Validation Gates

The vet/build phase should require these checks before analysis consumes #58:

- `manifest.json` reports `relationship = per_atom_substrate`,
  `relationship_kind = per_atom_aggregate`, `normalization = raw_lab_frame`
  or the vetted replacement string, and row alignment for every sidecar.
- `R == n_atoms * n_dft_frames` for dense all-atom runs.
- `(atom_index, original_frame_index)` is unique.
- DFT target sidecar rows match CSV `dft_present` flags.
- Present flags match finite sidecar payloads.
- Existing all-atom direct feature sidecars match the corresponding #58 columns
  cell-for-cell for APBS, AIMNet2, and MOPAC target-axis payloads.
- Backbone rows can reproduce the existing broad-backbone aggregate values
  within tolerance when selecting the six backbone strata and the same cutoffs.
- No `*_sources.csv` pair/source dump is emitted as part of the default #58
  substrate. Optional pair diagnostics must be named query outputs and managed
  as transient substrate.

## Open Questions For Vet

1. **Conditioning set.** Which proposed conditioning variables are required for
   the first fit, and which can wait for the pair-index query surface? Dominance
   fractions are useful but require top-source tracking or a query.

2. **Embedding placement.** Default recommendation is a separable f32 sidecar.
   Vet whether the run framework should drop/keep it independently from the
   classical substrate.

3. **Frame choice for generalized reducers.** Primary all-atom tensors should be
   lab/equivariant to match the node-store contract. Backbone-local compatibility
   sidecars are useful but should not drive the all-atom schema unless explicitly
   blessed.

4. **Per-mechanism aggregate choice.** For ring, use the Johnson-Bovey
   `jb_T0/jb_T2` source kernel sum, not only the old dipolar sum. For charge,
   use the existing `chargeKernelT2FromSources`-style q/r3 EFG-like T2, not the
   older `charge_dipole` `mu` alone. Confirm these are the blessed reducers for
   the first all-atom joint fit.

5. **Feature count versus atom count.** Global all-atom fits are likely
   determinate with 846 atoms, but narrow partitions may have fewer atoms than
   the full feature count, especially with the 256-d embedding. Decide the rule:
   ridge-only, PCA embedding, or drop embedding for thin partitions.

6. **Reducers not yet producing exactly this aggregate.** Ring JB aggregate,
   all-atom charge q/r3 aggregate, all-atom FF14SB field, all-atom valid-source
   McConnell, and dominance/count conditioning require C++ aggregate folds even
   though their source kernels/selectors already exist. This is implementation
   work, not a Python recompute.

7. **MOPAC charge source.** `Catalog::MopacCharge` is absent for per-frame
   charges; only `MopacChargeWelfordMean` exists as a static atom-axis charge.
   The first #58 substrate should not advertise per-frame MOPAC charge unless a
   trajectory time series is added upstream.

8. **ArraySpec additions.** The parent Python catalog does not yet describe all
   broad-backbone literature-kernel sidecars or the new #58 sidecars. The build
   phase must add ArraySpec entries with units, irreps, sign conventions,
   native axes, mechanisms, dtypes, and feature eligibility before consumers use
   the files.

## Pair-Index / Deep-Pairwise Query Surface

Status: **design carried forward**: the resident pointer/index idea landed as in-memory indexes plus named query outputs, not as a resident pair dump.
The per-atom aggregate substrate above is the reduced fit substrate. The
pair-index is the unreduced, queryable pointer surface. Both halves are part of
the same design and should be vetted together before either is built.

### Purpose

The pair-index is a queryable pointer index over the same source associations
that `RelationshipEngine::RunTraversal` already generates and the same
per-frame KD trees that `SpatialIndexSet` already owns. It keeps the deep
pairwise view available for Stage-1, de-circ, law-hunter, and reader navigation
without restoring the old source-payload dump.

A pair is a pointer tuple:

```text
(target_idx, source_idx, source_kind)
```

where `target_idx` identifies the target atom/frame case and `source_idx` plus
`source_kind` point back into the resident node store or source cloud:

```text
source_kind = atom_charge_site | ring_center | bond_midpoint |
              all_bond_midpoint | target_feature
```

Expected v1 mapping:

- atom charge sites point to atom nodes, with the charge source selected by the
  query (`ff14sb`, `aimnet2`, `mopac_welford_mean` when present);
- ring centers point to ring nodes / `CloudKind::RingCenters`;
- bond midpoints point to bond nodes / `CloudKind::BondMidpoints` or
  `CloudKind::AllBondMidpoints`, depending on whether the query wants only
  anisotropic McConnell bonds or every topology bond;
- target features such as APBS E-field, APBS EFG, AIMNet2 atom features, and
  MOPAC target-axis tensors are self-target pointers, not neighbor rows.

No source payload is copied into the pair. Per-pair quantities are computed
only when a query asks for them:

- `disp`, `r`, `inv_r3`, `cos_theta`, and `dipolar` use the existing C++ verb
  path (`verbs::pos`, `verbs::ringGeom`, `verbs::displacement`, bond-axis
  geometry, and charge-site displacement);
- per-pair T2/T0 kernels use the same source attachers / literature-kernel
  producers already used by `RunTraversal`, `RingCurrentNeighborhood`,
  `BroadBackbone`, and `AllAtomEquivariant`;
- `heavyParent`, ring self/bonded, bond endpoint-self, McConnell near-field,
  and producer-valid flags are pointer attributes computed with the same C++
  topology and geometry rules, not Python proxies.

Uncertainty to vet: the exact `source_kind` enum names can change during the
build. The contract is the pointer semantics and C++ dereference behavior, not
the spelling above.

### Lazy Versus Held

Default recommendation for v1: **lazy pair-index queries**.

Lazy mode regenerates adjacency from the per-frame KD trees and the
relationship selectors at query time. Storage is approximately zero beyond the
resident node store, `SpatialIndexSet`, topology, catalog, and frame map already
needed by the model. This is the lowest-risk first build because it surfaces
existing behavior and avoids committing to resident adjacency layout before the
v1 query set is known.

Optional held mode stores compact adjacency:

```text
target_idx, source_idx, source_kind, pointer_flags
```

For 1P9J, the parent node-store contract estimates about 30-60M pointer pairs.
At compact integer widths, held adjacency is expected to be about 0.3-0.5 GB.
That is acceptable as an optional acceleration/cache, but it should not be the
default unless lazy query time becomes a real bottleneck.

Expected relationship to assess:

| mode | expected v1 outcome | probability before build |
|---|---|---:|
| lazy from KD trees | enough for named query outputs and vetting, with minimal disk risk | high, about 0.75 |
| held adjacency | useful for repeated interactive navigation or many hunter scans | medium, about 0.55 |
| default pair dump | unnecessary for #58 and likely to recreate the old disk problem | low, about 0.10 |

### Vectorized Query API

The live API, when it exists, must be vectorized: one query returns arrays.
There must never be a Python call per pair.

Canonical query shape:

```text
pairs(stratum | atom_set | frame_range | case_filter, source_kinds, cutoff_policy)
    -> named arrays
```

Minimum returned pointer arrays:

```text
target_row_id
target_atom_index
frame_slot
h5_row
original_frame_index
source_kind
source_id
source_atom_index              optional / -1 for non-atom sources
source_cloud_index             optional / -1 when not KD-backed
source_category_ord
pointer_flags
```

Optional lazy-computed arrays, selected by query:

```text
disp_x, disp_y, disp_z
r
inv_r3
cos_theta
dipolar_3cos2m1_over_r3
kernel_T0
kernel_T2_0..4
source_value
q_over_r3
self_or_bonded
producer_valid
near_field
```

Named v1 query families to vet:

| query | answer | main consumer |
|---|---|---|
| `pairs(stratum, frame_range, source_kinds)` | pointer arrays plus requested lazy geometry/kernel columns | Stage-1/de-circ small pairwise views |
| `pairs(atom_set, frame_range, source_kinds)` | same, with explicit atom selection | targeted audits and reader navigation |
| `top_sources(stratum_or_atom_set, frame_range, mechanism, metric, k)` | top-k source pointers and contribution magnitudes per target row | dominance analysis and case inspection |
| `dominance_fractions(stratum_or_atom_set, frame_range, mechanism, metric)` | per-target `max_abs_contribution / sum_abs_contribution`, optionally top-k share | conditioning variables deferred by the aggregate spec |
| `hunter(law, constraints)` | least-confounded cases ranked by isolation, support, source dominance, and input modulation | #55 case selection |
| `reader_pairs(atom, frame_window, source_kinds)` | navigable pair rows for one atom over a small frame window | reader strips and explanatory views |
| `source_conditioning(stratum_or_atom_set, frame_range, mechanism)` | source-level support counts, distance bands, self/bonded counts, and source-category mix | all-atoms partition conditioning |

Dominance metrics should be input-side only. Examples:

```text
abs(kernel_T2) contribution
abs(kernel_T0) contribution
abs(q_over_r3)
abs(dipolar * source_value)
abs(source tensor T2 magnitude)
```

The query must not rank cases by DFT residuals or target values. It can join to
target identity for display, but selection metrics are source geometry,
provenance, support, and emitted input-side features.

Expected relationships to assess:

| query result | expected use | probability before build |
|---|---|---:|
| dominance fractions from top-source queries | stronger conditioning than raw neighbor counts for ring, charge, and McConnell partitions | high, about 0.70 |
| hunter least-confounded cases | finds cleaner per-law examples than global strata alone | medium-high, about 0.65 |
| reader atom/frame-window pairs | makes strips inspectable without carrying a resident pair table | high, about 0.75 |
| source-level conditioning for all-atoms partitions | improves partition diagnosis when aggregate features are ambiguous | medium, about 0.55 |

### Interim Delivery: Named Query Output

The parent contract is delivery-agnostic and the live pybind11 binding is
deferred. For the first #58 vet/build pass, the pair-index can be delivered as
a **small named query result** emitted by the run framework.

Interim form:

```text
query_results/
  <query_name>_manifest.json
  <query_name>_rows.csv                  identity, source ids, flags, categories
  <query_name>_geometry.npy              optional selected geometry columns
  <query_name>_kernel_T2.npy             optional selected T2 kernel columns
  <query_name>_kernel_T0.npy             optional selected T0 kernel column
  <query_name>_contribution.npy          optional scalar contribution metric
```

The manifest records:

```text
query_name
relationship_or_mechanism
source_kinds
atom_set_or_stratum
frame_range
cutoff_policy
lazy_columns_computed
units
irreps
sign_conventions
row_count
source = pair_index_query
delivery = transient_named_query_output
```

The run framework may drop old query outputs under the same managed-output
policy used for other rediscover products. Query outputs are for specific
analysis questions and should normally be KB-MB, not GB-TB. A larger flat
materialization is still allowed as a named one-off tool when a job explicitly
needs it, but it is not the default pair-index delivery.

### Boundary And Disk Rules

The default #58 build must not emit a resident `*_sources.csv` or full pair dump.
The pair-index default is lazy from KD trees plus topology/catalog dereference.
The aggregate substrate default is one row per `(atom, DFT-present frame)`.

Python reads query outputs and aggregate arrays. Python does not recompute
geometry, fields, projections, per-pair kernels, self/bonded status, or
dominance fractions from raw coordinates. Those are C++ query results because
they depend on the same verbs, source selectors, filters, and attachers as the
model.

Disk expectation:

| product | expected size | default? |
|---|---:|---|
| per-atom aggregate substrate without embedding | about 350-500 MB | yes |
| optional AIMNet2 embedding sidecar | about 572 MB f32 | separable |
| lazy pair-index resident storage | approximately 0 | yes |
| held pair-index adjacency | about 0.3-0.5 GB | optional |
| named pair query output | KB-MB for v1 queries | yes, transient |
| full flat pair materialization | tens of GB for broad global jobs | one-off only |
| unbounded pair/source dump | can drift toward the old 68 GB problem or worse | no |

There is no 1 TB path in the default design. If a proposed query would generate
that scale, it is not a v1 named query; it needs a scoped materialization plan,
a drop-old policy, and an explicit consumer.

### What This Unblocks

The pair-index supplies the things the aggregate substrate intentionally does
not store:

- dominance fractions: the aggregate spec deferred
  `ring_abs_T2_dom_frac`, `charge_abs_T2_dom_frac`,
  `mc_abs_T2_dom_frac`, and related top-source measures here;
- #55 hunter case selection: least-confounded examples can be selected by
  source isolation, dominance, support, and input modulation without scanning a
  raw pair dump;
- reader strips: an `(atom, frame_window)` view can show the source pointers and
  lazy geometry/kernels for navigable cases;
- all-atoms partition conditioning: source-level mix, distance bands, and
  self/bonded support can condition aggregate fits without becoming aggregate
  feature columns by default;
- de-circ and Stage-1 audits: small pairwise views remain available as named
  products, but they no longer define the resident substrate.

Expected relationship to assess: using pair-index dominance and source-support
queries should make the per-atom aggregate fits easier to interpret, especially
where several mechanisms have similar aggregate magnitudes. Probability before
build: medium-high, about 0.65. The uncertainty is strongest for narrow
all-atom partitions where source support may be too thin for stable conclusions.

### Relationship To Per-Atom Aggregates

The two halves are complementary:

```text
per-atom aggregates = reduced fit substrate
pair-index          = unreduced queryable pointers
together            = #58
```

The aggregate substrate answers:

```text
For this target atom/frame, what are the reduced mechanism features, target,
identity, and conditioning variables for a joint all-atoms fit?
```

The pair-index answers:

```text
For this target atom/frame or selected case set, which source pointers produced
the support, what were their lazy geometry/kernel values, and which sources
dominated?
```

The aggregate substrate should not grow source rows to answer pairwise
questions. The pair-index should not become a second aggregate substrate. When a
query computes a reducer-like answer, it should be named as a query result and
kept transient unless it graduates into the vetted aggregate schema.

### Open Questions For Vet: Pair-Index

1. **Lazy versus held for v1.** Recommendation is lazy default, held adjacency
   optional after timing real named queries. Vet whether any first-pass consumer
   needs repeated scans badly enough to justify held adjacency immediately.

2. **Python reach before pybind11.** Recommendation is named-query emission
   through the run framework now, live binding later. Vet whether the first v1
   query set is small enough for this, or whether the binding should be pulled
   forward.

3. **V1 query set.** Recommended minimum is
   `pairs(...)`, `top_sources(...)`, `dominance_fractions(...)`, and
   `reader_pairs(...)`. Vet whether `hunter(...)` is v1 or a second pass built
   on the first named outputs.

4. **Contribution metrics.** Decide per mechanism whether dominance is based on
   T2 magnitude, T0 magnitude, `q_over_r3`, dipolar intensity, or another
   source-side scalar. Record units and sign/magnitude choices in the query
   manifest.

5. **Cutoff policy.** Use the same relationship cutoffs as the aggregate
   reducers unless a query explicitly asks for a diagnostic cutoff. Diagnostic
   cutoffs must be named in the manifest so they do not silently redefine the
   fit substrate.

6. **Pointer flags.** Bless the v1 flag set: `self_or_bonded`,
   `producer_valid`, `near_field`, `present`, and any source-kind-specific
   validity bits. These should come from C++ producer rules.

7. **Graduation rule.** If a pair-index query becomes a repeated input feature,
   decide whether it stays a transient query result or graduates into the
   per-atom aggregate schema with ArraySpec metadata and validation gates.

## Deep Per-Residue Identity (from our own topology)

Status: **partly landed**: perceived-topology identity is emitted in the substrate; deeper categorical-engine work is parked.
This extends the identity/stratum surface
above and the pair-index surface below it.

Boundary: identity is surfaced from the topology the model already perceived:
`atoms_category_info.npy`, `residues.npy`, `rings.npy`,
`ring_membership.npy`, and the typed reader objects `QtAtom`, `QtResidue`,
`QtRing`, `QtRingMembership`, and `QtTopology`. The substrate reads perceived
identity. It must not re-perceive atoms or residues from external IUPAC
naming/ordering tables, and it must not use atom-name order as chemistry.
`iupac_atom_name` is allowed only as the already-emitted name projection on
`QtAtomNames`, useful for display and joins; it is not a source of chemical
classification.

### Purpose

The first identity table keeps the all-atom substrate lean, but it still mostly
names the six backbone strata plus broad sidechain buckets. The deep identity
surface should make every atom row partitionable by the chemistry substrate the
model already holds: residue type, locant, branch, diastereotopic index,
planarity, polar-H kind, ring position, charge state, exchangeability,
equivalence class, and ring membership.

This serves two consumers:

- all-atoms partitioning, where sidechain and ring-aware strata prevent
  chemically different atoms from being pooled only as `sidechain_heavy` or
  `aromatic_H`;
- the ring-current habitat, where per-ring type, ring position,
  vertex-versus-substituent status, and fused-ring context are high-value
  conditioning variables for Johnson-Bovey aggregate fits, top-source queries,
  and hunter case selection.

### Atom And Residue Identity Fields

These fields are static on the atom/residue axes and are repeated on the
per-row CSV for `(atom, DFT-present frame)` unless the vet/build phase chooses
a compact atom-axis identity sidecar with manifest-declared row joins. Either
delivery is acceptable if consumers can join by `atom_index` without
recomputing identity.

| proposed substrate field | already in topology sidecar? | typed source/accessor | #58 exposure status | use |
|---|---|---|---|---|
| `iupac_atom_name` | yes: `atoms_category_info.iupac_atom_name` | `QtProtein::atomNames(i).iupac` / `QtAtomNames::iupac` | new explicit exposure; current spec only has generic `atom_name` as display projection | display, audit, ML/external shift joins; never classification |
| `locant_ord` | yes: `atoms_category_info.locant` | `QtAtom::locant` | already in first identity table; keep | sidechain position partition |
| `branch_outer_ord`, `branch_inner_ord` | yes: `atoms_category_info.branch_outer`, `branch_inner` | `QtAtom::branch.outer`, `QtAtom::branch.inner` | new substrate exposure | branch-aware sidechain split |
| `di_index_ord` | yes: `atoms_category_info.di_index` | `QtAtom::diIndex` | new substrate exposure | prochiral/paired-H split |
| `planar_group_ord` | yes: `atoms_category_info.planar_group` | `QtAtom::planarGroup` | already in first identity table; keep | amide, carboxylate, guanidinium, imidazole, aromatic partitions |
| `polar_h_kind_ord` | yes: `atoms_category_info.polar_h_kind` | `QtAtom::polarH` | already in first identity table; keep | exchangeable/polar-H partitions |
| `prochiral_ord` | yes: `atoms_category_info.prochiral` | `QtAtom::prochiral` | new substrate exposure | ProR/ProS conditioning |
| `backbone_role_ord` | yes: `atoms_category_info.backbone_role` | `QtAtom::backboneRole` | already in first identity table; keep | six backbone strata and sidechain exclusion |
| `amino_acid_ord` | yes: `residues.residue_type` and replicated as `atoms_category_info.residue_type` | `QtResidue::aminoAcid` via `QtAtom::residueIndex` | already in first identity table; keep | residue-aware sidechain strata |
| `equivalence_class` | yes: `atoms_category_info.equivalence_class` | `QtAtom::equivalenceClass` | new substrate exposure | equivalent atom grouping and grouped CV |
| `aromatic` | yes: `atoms_category_info.aromatic` | `QtAtom::aromatic` | new substrate exposure; ring rows also expose aromaticity through `QtRing::ringKind` | aromatic atom split |
| `formal_charge` | yes: `atoms_category_info.formal_charge` | `QtAtom::formalCharge` | new substrate exposure | charged sidechain partitions |
| `is_exchangeable` | yes: `atoms_category_info.is_exchangeable` | `QtAtom::isExchangeable` | new substrate exposure | exchangeable-H and H-bond habitat partitions |

`ff_atom_type_ord` remains useful and is already in the first identity table.
It is a typed AMBER force-field surface (`QtAtom::ffAtomType`), not an atom-name
string. It can support existing predicates such as
`QtAtom::IsAromaticRingHydrogen()` without re-reading naming tables.

### Ring Topology Distinctions

The ring identity fields should be available for two related cases:

- target atom identity: which ring(s) the target atom itself belongs to, and
  which ring position labels the perceived atom substrate assigned;
- ring-current source conditioning: which source ring type/kind/fusion context
  generated a contribution, including top-source and hunter query outputs.

The current topology can expose the following without new chemistry
perception:

| proposed substrate/query field | already in topology sidecar? | typed source/accessor | #58 exposure status | use |
|---|---|---|---|---|
| `ring_type_index_ord` | yes: `rings.ring_type_index` | `QtRing::TypeIndex()`; display only via `QtRing::TypeName()` (`PHE`, `TYR`, `TRP6`, `TRP5`, `TRP9`, `HIS`, `HID`, `HIE`, `PRO`) | new exposure on target-owned ring slots and ring-source query rows | split PHE/TYR/TRP/HIS/PRO ring habitats |
| `ring_kind_ord` | yes: `rings.ring_kind` | `QtRing::ringKind`, `QtRing::IsAromatic()` | new exposure | aromatic versus saturated ring conditioning |
| `fused_partner_ring_id` | yes: `rings.fused_partner_ring_id` | `QtRing::fusedPartnerRingId`, `QtRing::IsFused()` | new exposure | TRP fused-ring split and bridge/fusion context |
| `ring_atom_count` | yes: `rings.atom_count` | `QtRing::atomIndices.size()` after membership join; type cross-check via `QtRing::RingSizeValue()` | new exposure | five-, six-, and nine-membered ring split |
| `ring_atom_order` | yes: `ring_membership.ring_atom_order` | `QtRingMembership::ringAtomOrder` from `QtTopology::ringMembershipsForAtom()` / `ringMembershipsForRing()` | new exposure | canonical ring-walk position |
| `is_vertex` | yes: `ring_membership.is_vertex` | `QtRingMembership::isVertex` | new exposure; 1P9J fixture currently has all membership rows as vertices | vertex/substituent partition |
| `is_substituent` | yes: `ring_membership.is_substituent` | `QtRingMembership::isSubstituent` | new exposure; current fixture has all zero and the code comments reserve it for future extension | future substituent split without schema churn |
| `ring_position_primary_ord` | yes: `atoms_category_info.ring_position_primary` | `QtAtom::ringPositionPrimary` | already in first identity table; keep | atom-side ring position label |
| `ring_position_secondary_ord` | yes: `atoms_category_info.ring_position_secondary` | `QtAtom::ringPositionSecondary` | new exposure | fused/perimeter secondary label |
| `aromatic` | yes: atom flag in `atoms_category_info.aromatic`; ring aromaticity through `rings.ring_kind` | `QtAtom::aromatic`; `QtRing::IsAromatic()` | expose both atom and ring contexts where relevant | ring atom and ring-source conditioning |

Recommended target-row columns for atom-owned ring context:

```text
ring_membership_n
ring0_id, ring0_type_index_ord, ring0_kind_ord, ring0_atom_count,
ring0_fused_partner_ring_id, ring0_atom_order, ring0_is_vertex,
ring0_is_substituent
ring1_id, ring1_type_index_ord, ring1_kind_ord, ring1_atom_count,
ring1_fused_partner_ring_id, ring1_atom_order, ring1_is_vertex,
ring1_is_substituent
ring2_id, ring2_type_index_ord, ring2_kind_ord, ring2_atom_count,
ring2_fused_partner_ring_id, ring2_atom_order, ring2_is_vertex,
ring2_is_substituent
```

Use `-1` ordinals/ids and `0` booleans for absent slots. The `QtTopology`
comment says atoms may belong to up to 2-3 rings for fused TRP bridgeheads.
Uncertainty to vet: whether fixed `ring0..ring2` CSV columns are sufficient for
all supported proteins, or whether this should be a compact atom-axis
`per_atom_ring_identity` sidecar keyed by `(atom_index, membership_slot)`.
The important rule is not to collapse a multi-ring atom to one ring.

Recommended ring-source conditioning columns for aggregate rows or named
pair-index queries:

```text
ring_source_dominant_type_index_ord
ring_source_dominant_kind_ord
ring_source_dominant_atom_count
ring_source_dominant_fused_partner_ring_id
ring_source_dominant_is_fused
ring_source_type_support_counts_0..7       aromatic ring current types
ring_source_type_support_counts_8          saturated PRO ring, expected zero JB
```

These are source-side conditioning variables. They are computed by the C++
aggregate fold or pair-index query from `QtRing` and contribution magnitudes,
not from DFT residuals.

### Deepened Role And Stratum Labels

The existing role and stratum proposals stay in place. Deep identity extends
them by adding ring-aware and sidechain-aware refinements. The labels are
derived labels over the already-perceived typed fields above; they are not a
new perception pass.

Proposed role extensions:

```text
aromatic_H.PHE, aromatic_H.TYR, aromatic_H.TRP6, aromatic_H.TRP5,
aromatic_H.TRP9, aromatic_H.HIS, aromatic_H.HID, aromatic_H.HIE

aromatic_H.<ring_type>.<ring_position_primary>
aromatic_heavy.vertex.<ring_type>.<ring_position_primary>
aromatic_heavy.fused.<ring_type>.<ring_position_primary>.<ring_position_secondary>
ring_substituent.<ring_type>

sidechain_H.aliphatic.<amino_acid>.<locant>.<branch_outer>.<di_index>
sidechain_H.polar.<polar_h_kind>
sidechain_H.exchangeable.<polar_h_kind>
sidechain_heavy.<amino_acid>.<locant>.<planar_group>
sidechain_heavy.charged.<formal_charge>.<planar_group>
sidechain_heavy.aromatic.<ring_type>.<ring_position_primary>
sidechain_heavy.sulfur.<amino_acid>.<locant>
```

Proposed stratum extensions:

```text
aromatic_H_by_ring_type
aromatic_H_by_ring_type_position
aromatic_heavy_vertex
aromatic_heavy_fused
ring_substituent
ring_bridge_or_fusion

sidechain_heavy_by_residue_locant
sidechain_H_by_residue_locant_branch
sidechain_H_prochiral_di
polar_H_by_kind
exchangeable_H_by_kind
charged_heavy_by_formal_charge_planar_group
amide_sidechain_by_atom_identity
carboxylate_by_atom_identity
guanidinium_by_atom_identity
imidazole_by_atom_identity
sulfur_by_residue_locant
```

The partition code may choose coarser labels when a partition is too thin.
For example, `aromatic_H.TRP5.Heteroatom_NH` may be useful for inspection but
too narrow for a global coefficient table; it can roll up to
`aromatic_H.TRP5` or `aromatic_H_by_ring_type`.

### What This Serves

Expected relationships to assess:

| conditioning/partition | expected use | probability before fit/query |
|---|---|---:|
| `aromatic_H.<ring_type>.<ring_position>` plus source `ring_type_index_ord` | cleaner ring-current aggregate response and cleaner hunter ring habitats than one pooled `aromatic_H` stratum | high, about 0.75 |
| fused TRP slots using `ring_type_index_ord`, `ring_position_secondary_ord`, and `fused_partner_ring_id` | separates benzene, pyrrole, and perimeter contexts for target atoms and top ring sources | medium-high, about 0.65 |
| `sidechain_heavy.<amino_acid>.<locant>.<planar_group>` | sharper all-atom partitions for MOPAC, APBS/EFG, charge, and AIMNet2 residual structure than `sidechain_heavy` alone | medium-high, about 0.65 |
| `polar_H_by_kind` and `exchangeable_H_by_kind` | better conditioning for electrostatic, H-bond, and AIMNet2 charge-response groups | medium, about 0.55 |
| `branch_outer/inner`, `di_index`, `prochiral`, and `equivalence_class` | distinguishes paired hydrogens and grouped atoms for leave-group-out CV and modulation tests | medium, about 0.50 |
| `formal_charge` plus `planar_group` | separates carboxylate, guanidinium, termini, and neutral polar groups for charge-driven partitions | medium-high, about 0.60 |

For the ring-current aggregate specifically, the clean habitat is not just
`near_ring` or `aromatic_H`. It is the joint condition:

```text
target atom perceived identity
+ target atom ring position/membership when present
+ source ring type/kind/fusion context
+ source dominance/support
+ ring geometry columns (r, rho, z, in-plane angle, cos_theta)
```

The per-atom aggregate substrate should carry the stable identity and compact
ring-source conditioning. The pair-index should answer detailed source-ring
questions, including top source rings and hunter cases, without emitting a
resident source dump.

### Catalog And Validation Notes

`_catalog.py` already declares the topology sidecars:
`atoms_category_info` on native axis `atom`, `rings` on native axis `ring`, and
`ring_membership` on native axis `ring_membership`; all are mechanism
`topology` and not ML features. If the build phase emits a new compact
atom-axis ring identity sidecar instead of CSV slot columns, it needs its own
ArraySpec entry with `native_axis="atom"` or an explicit
`atom_ring_membership` axis name.

Validation gates to add when this is built:

- every deep identity column must match the typed object loaded from the same
  sidecar row: `QtAtom`, `QtResidue`, `QtRing`, or `QtRingMembership`;
- `iupac_atom_name` must match `QtAtomNames::iupac` and must not participate in
  role/stratum derivation;
- `ring_membership_n` and `ring0..ringN` slots must match
  `QtTopology::ringMembershipsForAtom(atom_index)`;
- `ring_atom_order` slots must preserve `ring_membership.npy` canonical walk
  order;
- `ring_position_primary_ord` and `ring_position_secondary_ord` must match
  `QtAtom` fields, even when a fixed ring-slot layout is used;
- `is_substituent` must be surfaced even when all rows are zero in the current
  fixture, so future perceived substituents do not require a schema rethink;
- role and stratum labels must be reproducible from emitted typed columns
  alone, without name-table lookups.

Open vet questions:

1. **Ring-slot layout.** Fixed `ring0..ring2` columns are simple and match the
   current reader comment, but a compact atom-axis membership sidecar may be
   cleaner if future topologies can exceed three memberships.

2. **String placement.** `iupac_atom_name` is useful in CSV for audit and joins,
   but it is a projection. Vet whether it belongs in the main rows CSV or a
   parallel label sidecar.

3. **Label granularity.** The deep role vocabulary can become sparse. Vet the
   first-pass roll-up rules for coefficient tables, hunter filters, and plots.

4. **Dominant source ring identity.** The aggregate schema can carry dominant
   source-ring conditioning if the reducer tracks top source contribution.
   Otherwise it should first appear as a named pair-index query result and
   graduate only after repeated use.
