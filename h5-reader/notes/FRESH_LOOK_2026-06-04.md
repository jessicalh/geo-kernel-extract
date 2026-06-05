# Fresh-look readiness review - 2026-06-04

Scope: read-only readiness review for the next session: integrate `wt-720-build` into
`h5-reader-pysr-spike`, run the 720 static ingest on the frozen May 13/14 data, build
the combined 1P9J+720 Model-1 e3nn shielding emulator, and continue reader UI work
(radius view, infelicities, Linux AppImage).

Read: today's checkpoint/handoff, the 720 recipe/data note, the static-ingest and
learning specs, the advisor precis copies, the UI/stabilisation notes, and the built
720 ingest source in `/shared/2026Thesis/nmr-720`. I did not build, run tests, run the
720 ingest, attempt a merge, or change code. This note is the only intended write.

## Verdict

**Not ready as-is for the whole next-session plan. Fix these first.**

The direction is good: the frozen 720 raw per-protein NPY directories are usable, the
built worktree has a real static-cohort emitter, and the reader UI has a solid model
spine. The blockers are not conceptual. They are contract mismatches a tired session
could paper over:

- The frozen 720 data is old-schema raw NPY output with no `.LGS` files. The built
  static ingest currently requires exact per-protein `.LGS` single-pose manifests,
  DFT `meta_json` entries, and current field widths.
- The built static target path is **raw ORCA text via `.LGS` DFT meta**, while the
  frozen NPY root also carries `orca_total.npy`/dia/para. Those are loaded as optional
  snapshot fields, but they are **not** the current `BuildTarget` source.
- Current EFG fields are 5-component T2 arrays. The frozen EFG-like arrays predate
  the EFG schema revision and must be converted/trimmed before strict loading.
- The existing e3nn script is a 1P9J ring-source, local-frame, train-atom-centered
  interpolation script. Model 1 needs a combined corpus with **raw absolute DFT
  shielding**, not mutation deltas and not per-atom-centered modulation.
- The built static emitter writes lean per-atom aggregate sidecars and no default
  per-source table. A source-sum e3nn needs either a C++ source-tensor/chewer path or
  an explicit decision to train from aggregate sidecars instead.

So: the next session can start, but should not jump straight to "merge, run 720, train
combined model" without a short version-surgery and model-contract pass.

## Readiness By Track

### 720 frozen data

The frozen raw data is usable. `LEARN_720_RECIPE_AND_DATA_2026-06-04.md` is right:
`/shared/2026Thesis/nmr-shielding/calibration/features/Stage1BMRB_20260513_topology`
has 720 complete per-protein directories with topology sidecars, positions, MOPAC,
ring/H-bond/etc. arrays, and absolute ORCA arrays:

- `orca_total.npy`, `orca_diamagnetic.npy`, `orca_paramagnetic.npy`
- `wt_shielding_diamagnetic.npy`, `wt_shielding_paramagnetic.npy`
- `delta_*` arrays from the mutation workflow
- topology sidecars: `atoms_category_info.npy`, `residues.npy`, `bonds.npy`,
  `rings.npy`, `ring_membership.npy`, `extraction_manifest.json`

Do **not** use the old Stage-1 workup's `feature_matrix.npy` or `target_matrix.npy`
for Model 1. That derived matrix is mutation-delta/matched-row material and excludes
mutation-site atoms. Model 1 wants all atoms and absolute WT DFT shielding.

### Built static ingest

The worktree implementation is substantial and mostly aligned with the spec:

- `FrameNpyLoader::LoadSnapshotExactFields` is a strict expected-file loader and does
  not enumerate `*.npy`.
- `RunLoader::LoadStaticPoseExact` requires exact `.LGS` paths for cohort entries.
- `RunStaticCohortSubstrateEmit` writes row-aligned `static_rows.csv`,
  `static_target_total_raw.npy`, dia/para raw targets, T0/T2 targets, local-frame T2
  where valid, invariants, grouped sidecars, manifests, and a disk budget gate.
- The static manifest records `dft_target_source = raw_orca_out via .LGS DFT meta`,
  `source_pair_default_emit=false`, `python_second_model_allowed=false`, and the
  T2-orientation policy.

The problem is drop-in readiness against the frozen May data:

- The frozen per-protein directories have no `.LGS`. The static cohort reader requires
  `proteins[].producer_lgs` pointing to an exact `.LGS` file.
- `DftShieldingLoader` requires `.LGS -> dft.frames[].meta_json -> files.out_primary`.
  The frozen directories have ORCA NPYs; the corresponding ORCA text lives under
  `calibration/{ID}/{ID}_WT_nmr.out`.
- `RunData::staticPoseExpectedFields()` currently marks several post-schema fields as
  required and validates catalog column counts. Old EFG-width arrays will fail this
  path unless adapted.
- The static preflight does not enforce exactly 720 proteins; it accepts any non-empty
  unique cohort. The spec says the 720 run should validate exactly 720.
- `static_manifest.json` sets `stage2_between_axis_statistics_ready=true` and
  `stage3_model1_e3nn_training_corpus_ready=true` even though oracle parity is recorded
  as not run by the static emit. Treat those flags as aspirational until the actual
  old-schema run and gate pass.
- The single-pose static path does not hydrate FF14SB charges from the frozen
  `.prmtop`; trajectory charge hydration is wired through `trajectory.topology_top`.
  If FF14SB charge/field slabs matter, the adapter must add a charge source or record
  FF14SB absent honestly.

### Merge readiness

I did not attempt a merge. From source inspection, `wt-720-build` touches shared
reader/rediscover surfaces:

- `h5-reader/src/io/FrameNpyLoader.*`
- `h5-reader/src/io/QtProteinLoader.*`
- `h5-reader/src/model/QtProtein.h`
- `h5-reader/src/rediscover/{RunData,Catalog,PerAtomSubstrate,main_extract}.*`

Most reader risk looks additive: normal reader `LoadFromManifest(single_pose)` still
uses lenient `LoadPose`; strict loading is used by rediscover static paths. `QtProtein`
only adds `RunLoader` as a friend. Still, this is shared IO/model code, so a post-merge
reader build and a real 1P9J open are not optional.

The nanoflann issue is mostly provenance, not merge content. The 720 worktree still has
an untracked symlink at `h5-reader/extern/nanoflann/nanoflann.hpp`, while the main tree
has a real vendored header. Verify the merged build uses the real header/CMake include
path and does not depend on the worktree symlink.

### Combined Model-1

The scientific target is clear: **raw absolute DFT shielding**, total tensor primary,
dia/para diagnostic, no mutation delta.

The implementation path is not yet clear enough to call ready:

- Existing `analysis/equiv_t2_e3nn.py` consumes `ring_current_sources.csv` plus
  `rediscover_ring_current_sources_target_local_T2.npy`.
- It filters to local-frame-valid ring through-space source rows.
- It de-means target and prediction per atom using train frames.
- It reports 1P9J within-axis modulation recovery.

That is fine as an interpolation v0. It is not the combined 1P9J+720 shielding
emulator yet. For the combined model, fix the contract first:

- Target: train on absolute total shielding (`static_target_total_raw` and/or
  `static_target_T2`/T0), not `target_matrix.npy`, not `delta_*`, and not only
  train-atom-centered modulation.
- Orientation: choose raw molecular/lab-frame equivariant training, or local-frame
  training with explicit `frame_valid`/variant policy. Do not pool arbitrary lab-frame
  T2 components in a scalar/non-equivariant way.
- Static rows have one frame per protein. The 1P9J per-atom de-meaning protocol does
  not transfer directly to 720 statics.
- Source inputs: either add a bounded C++ source-tensor/chewer export for static rows,
  or intentionally train a smaller first model from aggregate sidecars. Do not pretend
  the current static aggregate files are the same input object as `ring_current_sources.csv`.
- Splits: use protein/static-pose held-outs for the 720 leg and blocked/purged frames
  for 1P9J; all preprocessing must be train-only.

e3nn remains the right library choice. The blocker is the data contract, not MACE vs
e3nn.

### Reader UI

The reader is real enough for UI work, but not desk-ready in the requested sense:

- Opening 1P9J is likely if `.LGS` and runtime paths are valid, but runtime was not
  verified here.
- The infelicities list is credible: debug-flavoured toolbar, "Instrument" label,
  no obvious selection Clear, per-frame overlay log noise, default `npy:dssp_chi`,
  and `QtFrame::eeqCharge` returning placeholder zero.
- Radius/local isolation is not built. The right path is still a separate
  `displayMolecule_`, source/display index map, `AtomFilter::WithinRadiusOf`,
  isolate toggle/radius control/clear, and source-index-based picking/overlays.
- Stabilisation infrastructure is useful but not the radius view. Plane lock is
  camera-only; `TransformedConformation` is the display-space seam.
- Linux is source-buildable, not installable. There are no install rules, CPack,
  AppImage flow, `$ORIGIN` RPATH strategy, or clean-machine smoke process.

### Docs drift

The handoff now exists and correctly points to this file first. One doc drift remains:
there are duplicate learning/precis docs. The source-tree
`h5-reader/src/rediscover/SPEC_LEARNING_MODEL_2026-06-04.md` is updated; the
`doc/thesis-overview/SPEC_LEARNING_MODEL_2026-06-04.md` copy still contains stale
status text saying companion specs were not present. The thesis-overview precis body is
mostly aligned, but its status block is older than the source-tree precis. The next
session should treat the handoff's path choices as authoritative or de-duplicate.

## Version-Surgery Specifics

This is the exact reconciliation needed before running the frozen 720 ingest.

1. **Create a current static input surface.**
   The built static cohort path wants a JSON cohort manifest with unique protein IDs
   and exact `.LGS` paths (`producer_lgs`/`lgs_path`/`path`). The frozen data has only
   raw per-protein NPY directories. Either:

   - synthesize per-protein `.LGS` single-pose manifests pointing at each frozen
     feature directory and `extraction_manifest.json`; or
   - add an explicit old-schema loader that accepts raw `pose_dir` entries while still
     recording equivalent manifest provenance.

   Also add the missing exact-720 preflight: 720 unique protein IDs, no duplicates, no
   silent 719/721 success.

2. **Choose and wire the DFT target source.**
   The preferred low-risk route is to synthesize `.LGS` DFT entries with
   `frame_index=0` and `meta_json` files whose `files.out_primary` points to
   `calibration/{ID}/{ID}_WT_nmr.out`. That keeps the current `DftShieldingLoader`,
   raw ORCA parser, atom-count validation, dia+para check, and coordinate-alignment
   diagnostics.

   If using the frozen `orca_total.npy`/dia/para arrays instead, then the adapter must
   inject those arrays into `DftFrameSet`/`BuildTarget` as the canonical target. Merely
   loading `OrcaTotal` as an optional `FrameNpyLoader` field is not enough. Validate
   atom count and `total == dia + para` componentwise, and record that the target came
   from producer ORCA NPY rather than raw ORCA text.

3. **Convert old EFG schema to current EFG schema.**
   Current `QtFieldCatalog.gen.h` expects 5-column T2 EFG fields for:

   - `apbs_efg`
   - `aimnet2_efg`
   - `aimnet2_efg_aromatic`
   - `aimnet2_efg_backbone`
   - `mopac_coulomb_efg_backbone`
   - `mopac_coulomb_efg_aromatic`
   - any carried water/coulomb EFG fields (`water_efg`, `water_efg_first`,
     `coulomb_efg_*`)

   The frozen data predates the May 18 EFG revision and carries old 9-column EFG-like
   arrays. Adapter rule:

   - if old 9-column EFG is old spherical layout (`T0`, `T1[3]`, `T2[5]`), verify
     T0/T1 are structurally zero or negligible and write columns 4..8 as the current
     5-column T2 array;
   - if any old 9-column EFG is actually raw 3x3, decompose with the same C++ basis
     (`DecomposeLibrary`) and write the T2 five-vector;
   - record the conversion in static manifest provenance.

   Do **not** apply this trim to shielding arrays that are supposed to remain 9-column
   `0e+1e+2e` fields (`bs_shielding`, `mc_shielding`, `hbond_shielding`,
   `mopac_coulomb_shielding`, `mopac_mc_shielding`, `orca_*`, etc.). The current
   static catalog intentionally reads their T2 leg from columns 4..8.

4. **Handle AIMNet2 naming honestly.**
   Frozen data has `aimnet2_charges.npy`, `aimnet2_aim.npy`, `aimnet2_efg*.npy`, plus
   `aimnet2_polarisability.npy` and `aimnet2_polarisability_scalar.npy`.
   Current code treats charge-response-gradient fields as optional:
   `aimnet2_charge_response_gradient.npy` and
   `aimnet2_charge_response_gradient_scalar.npy`.

   Do not rename polarizability to CRG unless the lead explicitly vets that physics
   equivalence. Leave CRG absent with support flags, and optionally preserve
   polarizability as old-schema provenance/optional future feature.

5. **Decide FF14SB charge hydration.**
   Frozen calibration inputs have `{ID}_WT.prmtop`, not the trajectory
   `topology_top` path used by `LoadFf14sbChargesFromTopol`. The static loader will
   otherwise emit absent FF14SB charge support. Either add a vetted `.prmtop`/topology
   charge hydration step for static poses, or explicitly accept/record FF14SB absent
   and rely on MOPAC/AIMNet2/EEQ charge slabs.

6. **Keep provenance blunt.**
   The emitted manifest should say this is `old_schema_stage1_topology_20260513`
   adapted to current static schema. It should record:

   - old feature root path and digest/provenance;
   - `.LGS`/meta synthesis if used;
   - DFT target source (`raw_orca_out` vs producer `orca_*.npy`);
   - EFG conversion mode per field;
   - CRG absent / polarizability old-only status;
   - FF14SB charge source or absence;
   - memory strategy;
   - thin-support threshold (= 10 unless changed before run);
   - oracle parity status.

## Fix-These-First Verdict

ready y/n: **N**.

top 3 must-fix-before-next-session items:

1. Build the old-schema adapter/staging manifest before the 720 run: exact 720 cohort, per-protein `.LGS`/DFT meta or an explicit old-schema loader, current EFG 5-column conversion, target provenance, and FF14SB/AIMNet2 handling.
2. Fix the combined Model-1 training contract: absolute raw total DFT shielding target, static-compatible preprocessing/splits, raw/equivariant vs local-frame orientation policy, and a source-tensor/chewer or explicit aggregate-sidecar input path.
3. Treat the 720 merge as unproven until the lead does the actual merge check and then rebuilds/runs the reader and extractor from the merged tree, with the real nanoflann header and no symlink dependency.

version-surgery specifics:

- Use frozen raw per-protein NPY directories, not the old Stage-1 `target_matrix.npy`.
- Synthesize/adapter-provide current static manifests because frozen dirs have no `.LGS`.
- Feed Model 1 absolute WT DFT shielding: `orca_total`/dia/para or raw ORCA WT text, never `delta_*`.
- Convert old 9-column EFG surfaces to current 5-column T2 EFG; leave true 9-column shielding surfaces alone.
- Keep AIMNet2 CRG absent unless vetted; do not silently map polarizability to CRG.
- Hydrate FF14SB charges from a vetted static topology source or record them absent.
- Record all old-schema conversions and do not mark corpus/statistics ready until the actual run and oracle/post-merge gates pass.
