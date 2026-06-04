# SPEC: 1P9J equivariant interpolation graph scoping - 2026-06-04

Status: scoping/build-plan document; trued 2026-06-04. The companion
`SPEC_LEARNING_MODEL_2026-06-04.md` exists, and the bounded interpolation build
later landed in `BUILD_INTERP_RESULT_2026-06-04.md` with clean held-out metrics
(T2 component R2 +0.4776, |T2| r +0.8188, T0 modulation R2 +0.4537).

## Verdict

**GO for the bounded 1P9J within-axis advisor graph.**

The realistic deliverable is a simple e3nn equivariant interpolation harness that predicts held-out DFT shielding modulation on 1P9J and emits one advisor-facing graph. It is close because the real e3nn T2 source-sum model already exists, the clean split/centering protocol now exists, the 1P9J substrate exists, and the recorded e3nn result already reproduces the old T2 signal: frame-split T2 R2 about 0.466 and |T2| r about 0.757.

**NO-GO for pretending this is the full Stage-3 all-atom source-level chewer.** The current live all-atom Build1 substrate is aggregate-only, with no resident source CSV. A full all-atom variable-source e3nn predictor needs the chewer/per-source carrier work or a bounded source export. That is Stage-3 build work, not the advisor graph.

The graph can still truthfully say: "we are heading toward the Stage-3 equivariant model." It must not say: transferability, between-axis recovery, BMRB validation, process model, or destination.

## Grounding Read

Primary local grounding:

- `CODEX_BRIEF_SPEC_INTERP_SCOPING_2026-06-04.md`: this task and hard rules.
- `STATE.md` and `NOW.md`: current Stage-2 status, true-LOAO correction, live substrate, and reporting standards.
- `analysis/equiv_t2_e3nn.py`: ring-current e3nn source-sum exemplar.
- `analysis/equiv_t2_backbone_e3nn.py`: heterogeneous source-kind e3nn exemplar.
- `analysis/equiv_t2_efg_e3nn.py`: clean Schur/e3nn-style EFG protocol template.
- `analysis/e3nn_protocol.py`: blocked/purged split, train-only centering, train-only normalization.
- `analysis/change_of_basis.py`: the only allowed library-T2 to e3nn-2e basis map.
- `analysis/PATTERNS.md`: Python consumes emitted substrate only; equivariant means e3nn.
- `analysis/FINDINGS.md` and `MODEL_PLACEMENT_PROPOSAL.md`: recorded ring/T2 results and model-boundary rationale.
- `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1`: live aggregate substrate with T2/T0 targets.
- `/tmp/rediscover-runs/2026-06-04-true-loao`: true between-axis correction.
- `../app` and `../model`: reader signal/display surface for optional overlay.

Could not ground directly in this pass:

- Named memories from the brief were not directly accessible; their facts were grounded through `STATE.md`, `NOW.md`, `analysis/FINDINGS.md`, and the brief itself.
- `SPEC_LEARNING_MODEL_2026-06-04.md` now exists; use it with
  `BUILD_INTERP_RESULT_2026-06-04.md` for the Stage-3/down-payment status.
- The historical source-e3nn run directory `/tmp/rdc-composed` was recorded in `STATE.md`, but I did not find a current `ring_current_sources.csv` or `broad_backbone_sources.csv` under `/tmp` in this read-only pass. The build must preflight this.

## What This Is

The deliverable is a **1P9J within-axis e3nn interpolation tool**, with the advisor graph as the first-class output.

It predicts DFT shielding across the trajectory's sampled geometries:

- primary target: total shielding T2, five components, never scalar-collapsed;
- companion target: sigma_iso/T0, alongside the T2 result;
- split: held-out frame block with purged neighboring frames;
- axis: within-atom modulation on one protein;
- inputs: emitted C++ substrate only;
- model: e3nn equivariant predictor extended from the existing exemplar, not a hand-rolled second model.

The honest reading:

- This is a geometry-sampler interpolation: each frame is an instantaneous geometry -> shielding sample.
- This is within-axis only: it tests held-out geometries for the same 1P9J protein, not transfer to new proteins or new atom environments.
- This is correlate-not-match: the graph should show recovery of DFT modulation, not claim chemical-shift accuracy.
- This is direction, not destination: a visible down-payment toward Stage 3, not the finished Stage-3 conditioned GNN.
- This is one protein, thin in the ring/source subsets, and must support-flag atom counts and effective-N.

## What Already Exists

### The e3nn model spine exists

`analysis/equiv_t2_e3nn.py` is already the correct model class for the graphable source-level result:

- e3nn `o3.spherical_harmonics` builds the `2e` angular channels.
- e3nn `TensorProduct` / `FullyConnectedTensorProduct` handles the `1o x 1o -> 2e` cross term.
- an invariant radial MLP gates the angular terms;
- source contributions are scatter-added to `(atom, frame)` groups;
- target T2 is the C++-emitted local-frame/library-basis sidecar mapped through `change_of_basis`;
- there is no Python tensor projection fallback.

`STATE.md` records that this e3nn rebuild passed the parity gate: frame-split T2 R2 about 0.466 and |T2| r about 0.757, reproducing the retired hand-rolled result. It also reports leave-atoms-out T2 about 0.370, explicitly thin and not the gate.

### The clean protocol exists

`analysis/e3nn_protocol.py` now provides:

- contiguous held-out frame block;
- neighboring-frame purge;
- train-only per-atom target centering;
- train-only feature normalization;
- explicit train/test masks so purged rows cannot leak into training.

`POSTMORTEM_E3NN_PROTOCOL_FIX_2026-06-04.md` says the code paths were aligned to this protocol, but the e3nn re-run after the protocol fix is held for the lead. So the graph build should expect a small metric move, not assume the old clean-vs-leaky verdict has been re-run.

### The 1P9J aggregate substrate exists

The live Build1 substrate is present at:

```text
/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1
```

It contains:

- `per_atom_substrate_rows.csv`
- `per_atom_substrate_manifest.json`
- `per_atom_substrate_column_specs.json`
- `per_atom_substrate_target_T2.npy`
- `per_atom_substrate_target_T0.npy`
- `per_atom_substrate_features_classical.npy`
- `per_atom_substrate_features_conditioning.npy`
- `per_atom_substrate_features_ring_paths.npy`
- `per_atom_substrate_features_method_paths.npy`
- `per_atom_substrate_aimnet2_embedding.npy`

Its manifests report:

- 846 atoms;
- 660 DFT frames;
- 558,360 aggregate rows;
- T2 components frame-aligned;
- `target_T0` present as `dft_sigma_iso`;
- sidecar row alignment by `row_id == frame_slot * n_atoms + atom_index`.

This is enough for aggregate/T0 companion checks and for a weaker aggregate-shadow e3nn fallback. It is not, by itself, the source-row input for the true source-sum e3nn model.

### The Stage-2 within-axis evidence is strong enough for context

Current within-axis facts:

- unified D_ab-sum within-frameblock R2 about 0.432 on 25 atoms;
- charge q/r3 within R2 about 0.278;
- ring-current within R2 about 0.275;
- e3nn source-sum ring T2 R2 about 0.466 and |T2| r about 0.757;
- ring scalar held-out-atom/coupled within-modulation result: k about 21 ppm*A^3, coupled within-atom R2 about 0.62.

These numbers make the graph plausible. They do not make a transferability claim.

### The between-axis correction is load-bearing

`POSTMORTEM_TRUE_LOAO_2026-06-04.md`, `STATE.md`, and `NOW.md` retract the prior Stage-2 LOAO/between readings:

- charge true-LOAO R2 about 0.036;
- ring true-LOAO R2 about -1.04;
- unified true-LOAO overfit/null;
- one-line read: true between-atom recovery on 1P9J is approximately null.

Therefore this spec must not quote the old Stage-2 LOAO numbers as between recovery. Any held-out-atom result in the advisor graph is a within-modulation audit, not between-axis transferability.

## Feasibility

### Recommendation

Build the graph as **Track A: source-e3nn v0**, with the current e3nn ring/source model as the primary path.

This is the closest honest graph because it is already the real equivariant source-sum model:

```text
emitted source geometry -> e3nn angular/tensor-product channels -> scatter-pool -> T2 (+ T0 head) -> held-out DFT graph
```

The build delta is packaging:

- locate or regenerate the bounded source substrate used by `equiv_t2_e3nn.py`;
- reuse the e3nn class and shared protocol;
- add a T0 scalar head inside the same model;
- write prediction/metric artifacts;
- render the advisor graph.

If source rows are unavailable and the lead does not approve bounded re-emission, then Track A must stop before training. In that case the only graphable fallback is Track B below, which is useful but weaker.

### Fallback

Track B is an **aggregate-shadow e3nn v0** on the Build1 aggregate substrate. It would consume emitted irreps and scalar conditioners from `per_atom_substrate_*` sidecars, not source rows.

That fallback can still be e3nn and T2-valued, but it is not the full variable-source Stage-3 architecture. Label it:

```text
aggregate-shadow e3nn v0; emitted source sums, not full source chewer
```

Use Track B only if Track A's source input cannot be materialized within the lead-approved build loop. Track B should not be sold as "the" Stage-3 model; it is a readable direction graph from emitted shadows.

### Hidden work assessment

Very close for the advisor graph:

- model code: present;
- clean protocol: present;
- basis map: present;
- prior metrics: present;
- T0 target: present in the Build1 aggregate substrate;
- graph/rendering: small analysis-output delta.

Not close for full all-atom Stage-3:

- no all-atom source-level chewer is present in the live Build1 aggregate substrate;
- no Python protein/spatial graph may be constructed to fake it;
- source-level all-atom model should wait for the chewer/per-source carrier.

## Model Scope

### Primary model

Extend the existing e3nn exemplar. Do not author a second model.

Minimum primary model:

- keep the current T2 `2e` source-sum path from `analysis/equiv_t2_e3nn.py`;
- add a companion `0e`/scalar T0 output head inside the same module;
- keep train-atom centering for both T2 and T0;
- report T2 as primary and T0 as companion.

The T0 head must be part of the same e3nn/DeepSets source-sum model. Do not create a separate sklearn/ridge scalar fit for the advisor headline.

For display, the model can report both:

- modulation prediction: target and prediction centered by train-atom mean;
- baseline-restored prediction: add the train-atom DFT mean back for an intuitive predicted-vs-DFT scatter.

The primary R2 should be the modulation R2. The baseline-restored scatter is for interpretation only, because absolute sigma_iso can be dominated by per-atom chemical baseline.

### Input boundary

Allowed:

- emitted CSV/NPY substrate files;
- emitted source vectors, source normals/axes, invariant scalars, source/group IDs;
- emitted target T2/T0 sidecars;
- `change_of_basis.get_C()`.

Forbidden:

- opening `trajectory.h5` in Python;
- recomputing Pople, charge, field, McConnell, or any source geometry in Python;
- projecting raw tensors into T2 in Python;
- building a Python protein/spatial graph, neighbor list, or secondary aggregate;
- selecting rows or atoms by DFT residual/fit quality;
- using old LOAO numbers as between-axis evidence.

## Held-Out Protocol

Use the shared clean e3nn protocol:

1. Split by a contiguous held-out frame block, default 20 percent.
2. Purge at least one neighboring frame on each side of the held-out block.
3. Center T2 and T0 per atom using training rows only.
4. Normalize invariant/source features using training source rows only.
5. Train only on `g_tr`.
6. Score only on `g_te`.
7. Keep purged rows out of both train and test metrics.
8. Report leave-atoms-out only as an audit, not as the graph gate, because source support is thin on 1P9J.

The graph caption must state "held-out frames, train-atom baseline removed/restored as labelled." If a baseline-restored panel is shown, it must also show the modulation R2.

## Advisor Graph Design

The graph is the deliverable. Prefer one file:

```text
interp_1p9j_advisor_graph.png
```

and optionally a PDF:

```text
interp_1p9j_advisor_graph.pdf
```

### Main layout

Use a 2 x 2 figure.

Panel A: **T2 held-out component scatter**

- x: held-out DFT T2 components, de-meaned by train-atom mean;
- y: e3nn predicted T2 components, same centering;
- one point per component per held-out `(atom, frame)`;
- y=x line;
- annotate 5-component R2;
- color by atom role or source subset; keep thin subsets visible.

Panel B: **T2 magnitude recovery**

- x: held-out DFT `|T2|` modulation;
- y: predicted `|T2|` modulation;
- annotate correlation r and optional magnitude R2;
- this is rotation-invariant and advisor-readable.

Panel C: **sigma_iso/T0 companion scatter**

- primary x/y: held-out DFT vs predicted T0 modulation;
- annotate T0 modulation R2;
- optional inset: baseline-restored sigma_iso scatter with train-atom means added back;
- caption says the T0 head is companion, not the T2 thesis target.

Panel D: **direction caption / model sketch**

Compact text plus a small flow:

```text
1P9J geometry sampler -> emitted C++ substrate -> e3nn source-sum predictor -> held-out DFT
```

Required caption language:

```text
Within-axis 1P9J interpolation only. The trajectory samples instantaneous geometries; this is not a dynamics/process model. The graph shows correlate-not-match recovery of held-out DFT modulation and is a direction signal toward the Stage-3 equivariant model, not a transferability or BMRB validation claim.
```

### Optional secondary plots

Only if time remains after the main graph:

- per-atom mini time-series for 2 or 3 coupled atoms, held-out block shaded;
- residual histogram for T2 components;
- support table: atoms, frames, source rows, N_eff, train/test frame counts;
- baseline-restored absolute sigma_iso panel.

These are optional. Do not let them displace the main predicted-vs-DFT scatter.

## Optional Reader Overlay

Reader crossover is optional and not part of the GO gate.

The reader already has:

- ORCA DFT descriptors in `TrajectorySignalCatalog.cpp`;
- `DftShieldingStore` sampling total/dia/para T0 and T2 magnitude by atom/frame;
- dashboard strips and static atom-color modes.

What is missing is a predicted-shielding artifact source. A demo overlay would need a small loader/descriptor path for the prediction artifact, then the dashboard could compare predicted vs ORCA DFT strips or color residuals. This is useful for a Track-D demo, but it is UI/loader work and should not be mixed into the first advisor graph gate.

## Build Sketch For The Gated Loop

### Output directory

Use a fresh, explicit output directory under `/tmp/rediscover-runs`, for example:

```text
/tmp/rediscover-runs/2026-06-04-interp-1p9j-e3nn
```

Drop-old only by explicit full path. Never glob. Never write bulky artifacts under `/shared`.

### Step 0: preflight gates

Read only at preflight:

- `analysis/PATTERNS.md`
- `analysis/e3nn_protocol.py`
- `analysis/change_of_basis.py`
- `analysis/equiv_t2_e3nn.py`
- input substrate manifest(s)

Hard preflight checks:

- disk: `/tmp` free >= 20 GiB;
- disk: `/tmp/rediscover-runs` total < 15 GiB before writing, keeping Build1 and Build4;
- `change_of_basis.get_C()` orthogonality max error < 1e-12;
- source-mode input files exist, or stop for lead approval to re-emit:
  - source CSV with emitted `disp_local`, `source_normal_local` or parity-safe axis fields, invariant scalars, group IDs;
  - emitted target local T2 sidecar;
  - T0 target column or sidecar;
- aggregate Build1 files exist if T0 companion or fallback is used:
  - `per_atom_substrate_rows.csv`;
  - `per_atom_substrate_target_T2.npy`;
  - `per_atom_substrate_target_T0.npy`;
  - `per_atom_substrate_column_specs.json`;
- no code path opens `trajectory.h5`;
- no code path computes forbidden kernel/projection expressions.

If source-mode files are absent and the lead approves bounded re-emission, the re-emission is a separate gated action in the build loop. It is not Python physics, but it is still a build step and must be explicitly logged.

### Step 1: extend the e3nn harness, not the model family

Implement by factoring or extending the existing exemplar, not by starting a new model:

- reuse the e3nn T2 source-sum module;
- reuse `e3nn_protocol` split/centering/normalization;
- add a `0e` T0 companion output head to the same module;
- keep T2 loss primary;
- scale T2/T0 losses from train rows only so T0 cannot dominate the T2 fit;
- preserve `--cross exact|learnable|both` behavior if using the ring source path;
- write all metrics for exact and learnable if both are run, but choose the graph model by held-out T2 R2.

Do not add a second scalar model. Do not add a Python source aggregator.

### Step 2: train and score held-out frames

Run with:

- fixed random seed already used by the exemplar;
- blocked/purged split;
- train-only feature normalization;
- train-only per-atom target centering;
- T2 5-component R2;
- `|T2|` correlation;
- T0 modulation R2;
- optional leave-atoms-out audit, reported separately and support-flagged.

Expected source-mode sanity range:

- If the same canonical ring source substrate is used, T2 R2 should be near the recorded 0.466 and |T2| r near 0.757.
- A clean-protocol re-run may move the number; a large collapse must be investigated before graphing.
- A crude graph is still useful if the held-out T2 recovery is positive and the caption is honest.

Hard graph gate:

- prediction rows non-empty;
- held-out rows non-empty;
- no train/test leakage in split audit;
- T2 R2 finite;
- graph file exists and is non-empty;
- metrics JSON and prediction CSV row counts match.

Interpretation gate:

- T2 R2 >= 0.25: advisor graph is strong enough for "direction" label.
- 0.05 <= T2 R2 < 0.25: graph is weak but can be shown as a crude early direction if labelled plainly.
- T2 R2 < 0.05 or negative: NO-GO for advisor graph unless the lead explicitly wants a negative-control slide.

### Step 3: write artifacts

Required artifacts:

```text
interp_1p9j_predictions.csv
interp_1p9j_metrics.json
interp_1p9j_advisor_graph.png
interp_1p9j_run_audit.json
```

Recommended fields in `interp_1p9j_predictions.csv`:

- `atom_index`
- `h5_row`
- `original_frame_index`
- `split` (`train`, `test`, `purged`)
- `target_T2_0..4_centered`
- `pred_T2_0..4_centered`
- `target_T0_centered`
- `pred_T0_centered`
- optional baseline-restored T2/T0 columns;
- atom role/stratum/source-support fields if emitted.

Recommended fields in `interp_1p9j_metrics.json`:

- model path and mode: source-e3nn or aggregate-shadow e3nn;
- input directory and required files;
- row counts, atoms, frames, source rows;
- split strategy, test frames, purged frames, cross-split lag-1 adjacency count;
- T2 component R2;
- `|T2|` r;
- T0 modulation R2;
- optional baseline-restored R2;
- optional leave-atoms-out audit;
- support flags;
- disk audit;
- no-H5/no-recompute boundary audit.

### Step 4: SETI report

The post-run summary should be numbers and shapes first:

- T2 scatter shape and R2;
- `|T2|` scatter/correlation;
- T0 companion R2;
- support counts;
- whether source-mode or aggregate-shadow mode was used;
- exact caveat label.

Allowed bucket labels:

- "direction signal";
- "weak direction signal";
- "negative/NO-GO graph";
- "source-e3nn v0";
- "aggregate-shadow e3nn v0".

Avoid:

- "validated";
- "transferable";
- "between recovery";
- "BMRB-ready";
- "process model";
- "solves Stage 3".

## What Proves It Worked

The build worked if:

1. The run consumes only emitted substrate files.
2. `change_of_basis` is the only basis conversion.
3. The split audit shows no held-out leakage.
4. The advisor graph exists and has non-empty held-out predicted-vs-DFT scatter.
5. T2 5-component R2 and `|T2|` r are reported.
6. sigma_iso/T0 companion R2 is reported or the missing reason is fail-loud and grounded.
7. The caption says within-axis, geometry-sampler, correlate-not-match, direction-not-destination.
8. The post-run report states whether the input was true source-e3nn or aggregate-shadow fallback.

## Scope Boundaries For The Lead

GO:

- one 1P9J graph;
- held-out frame interpolation;
- T2 primary, T0 companion;
- existing e3nn model family;
- realistic labels.

Not in this build:

- full all-atom Stage-3 chewer;
- 720-WT transfer;
- BMRB training or validation;
- reader overlay unless separately approved;
- re-opening Stage-2 true-LOAO;
- new Python physics or source aggregation.

## Final GO/NO-GO

**GO** for a bounded source-e3nn v0 advisor graph, because the e3nn predictor and clean protocol already exist and previous 1P9J held-out T2 recovery is strong enough for a "direction" plot.

**Conditional GO** for an aggregate-shadow fallback if source rows cannot be materialized in the lead-approved loop. It must be labelled weaker.

**NO-GO** for claiming the current live aggregate substrate alone proves the full all-atom Stage-3 source-level model. That requires the chewer/per-source carrier and belongs to the later Stage-3 build.
