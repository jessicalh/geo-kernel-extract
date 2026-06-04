# SPEC maths-fix and review cycle - 2026-06-04

Status: DRAFT stash for lead plan-vet. This document is a librarian/assembler pass only.
It consolidates the already-recorded audit outcomes, the two parked verdict re-runs, the
lead+Claude maths agenda, and the review loop. It makes no maths judgment, re-opens no
disposition, and resolves no open question.

Lane rule: every unresolved item below is explicitly for lead+Claude - not specified here.
Codex may grind future run mechanics only after lead go. The lead owns all git.

Grounding limits of this assembly pass:
- No git was run. Branch and commit/hash-like identifiers are copied from `NOW.md`,
  `STATE.md`, and the postmortems, not independently verified here.
- No fit, e3nn, build, extractor, or re-run was run. Existing run outputs were not
  regenerated.
- The documented `/tmp` e3nn substrate candidates were not live-grounded in this pass;
  the future runner must verify required sidecars before firing.
- The named memories were not accessible as files here; this spec grounds their content
  only where the local docs cite them.

## Part 1 - Parked verdict re-runs

### 1A. LOAO re-fit on the corrected true-between path

Purpose: confirm the retraction using the corrected true between-atom path. This is not a
new between claim. It is a confirmation-of-retraction check.

Sources:
- `MATHS_AUDIT_CHECKLIST_2026-06-04.md`: issue #2 disposition.
- `POSTMORTEM_TRUE_LOAO_2026-06-04.md`: validated true-LOAO diagnostic.
- `POSTMORTEM_LOAO_FIX_2026-06-04.md`: code-only main-path fix; full re-fit held.
- `analysis/stage2_law_fits.py`: `--true-loao`, `stage22_evaluate_loao`,
  `stage23_eval_probability`, and shared `true_loao_*` atom-mean machinery.

Substrate it reads:
- `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`
- Required files are the Build4 per-atom substrate CSV/NPY sidecars and CaseHunter
  manifests, including `per_atom_substrate_column_specs.json`,
  `per_atom_substrate_rows.csv`, `per_atom_substrate_target_T2`, the feature sidecars,
  dominance sidecars, and dominance bins.
- Forbidden inputs, per the run-audit pattern: `trajectory.h5`, ORCA outputs,
  per-source dumps, and older rediscover fitting directories.

Exact command set, subject to lead go:

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover

ls \
  /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_column_specs.json \
  /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_rows.csv \
  /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_target_T2.npy

/usr/bin/python3 analysis/stage2_law_fits.py \
  --substrate-dir /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4 \
  --out-dir /tmp/rediscover-runs/2026-06-04-true-loao-refit \
  --true-loao \
  --true-loao-permutations 1000
```

Production-table refresh, only if lead explicitly wants the main Stage 2 tables regenerated
after the code-only fix:

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover

ls \
  /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_column_specs.json \
  /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_rows.csv \
  /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_target_T2.npy

/usr/bin/python3 analysis/stage2_law_fits.py \
  --substrate-dir /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4 \
  --out-dir /tmp/rediscover-runs/2026-06-04-stage2-fits-loao-refit
```

Whether the diagnostic-only command is sufficient, or the full production-table refresh
must also be run, is for lead+Claude - not specified here.
If any 1A preflight `ls` fails, stop and ask for lead direction; do not silently
re-extract or substitute another substrate.

Held-out gate:
- TRUE LOAO trains on training-atom means only.
- Held-out atom means are transformed with training-atom feature and target statistics.
- Scores are on physical held-out atom means.
- The own-heldout-atom centering path is retired from the main LOAO/between path and kept
  only as an explicit historical diagnostic.

Anti-circular gate:
- Selection remains from input-side CaseHunter/dominance strata, not DFT recovery.
- The run reads Build4 substrate sidecars only; no trajectory, ORCA, per-source dumps,
  or old fitting directories.
- `change_of_basis.get_C()` orthogonality is recorded before fitting; current documents
  record `|C.T C - I|max=1.11e-16`.

Disk gate:
- Output must stay under `/tmp/rediscover-runs`; never under `/shared`.
- `analysis/stage2_law_fits.py` disk guard requires >=20 GiB filesystem free and
  `/tmp/rediscover-runs` <15 GiB before write.
- Use a fresh lead-approved output directory unless the lead explicitly accepts
  drop-existing behavior.

SETI gate:
- Report numbers, null positions, support flags, and artifacts. Do not editorialize.
- Treat this as confirmation of the retraction, not as a new 1P9J between result.
- If the rerun differs materially from the recorded diagnostic, stop and hand to
  lead+Claude - not specified here.

Expected verdict, copied as the held expectation from the audit docs:
- 1P9J true between-atom recovery remains approximately null: charge about `0.036`
  (`0.03`-class), ring about `-1.0`, unified negative/extreme overfit with null-class
  position. This is confirmation-of-retraction only, not a result to quote as between.
- Every prior LOAO/between number remains retracted. Any future change to that closed
  disposition is outside this spec and for lead+Claude - not specified here.
- Between / transferability / combine-depth remains deferred to the 720-WT instrument.
  Any future change to that closed disposition is outside this spec and for lead+Claude -
  not specified here.

Output that proves it:
- Diagnostic command: `true_loao_recovery.csv`, `true_loao_unified_drop_one.csv`,
  `true_loao_atom_predictions.csv`, `true_loao_null_scores.csv`,
  `true_loao_run_audit.json`, and the generated true-LOAO postmortem.
- Production refresh command: `stage2_kernel_summary.csv`,
  `stage2_path_coefficients.csv`, `stage2_unified_dab_summary.csv`,
  `stage2_unified_dab_intensities.csv`, `stage2_run_audit.json`, plus the convention
  ledger.
- Proof condition: the output tables use true atom-mean LOAO, not the historical
  own-heldout-atom centering diagnostic.

### 1B. e3nn clean-protocol re-run

Purpose: obtain the actual clean-protocol number for the within combine's "three paths
agree" statement after ring/broad e3nn were aligned to the shared blocked/purged EFG
protocol. This is the held verdict for issue #3.

Sources:
- `POSTMORTEM_E3NN_PROTOCOL_FIX_2026-06-04.md`: code-only fix, no e3nn re-run.
- `analysis/e3nn_protocol.py`: shared split, train-only centering, train-only
  normalization.
- `analysis/equiv_t2_e3nn.py`: ring-current e3nn fitter.
- `analysis/equiv_t2_backbone_e3nn.py`: broad-backbone e3nn fitter.
- `analysis/ENV.md`: system python and CUDA library path.
- `STATE.md` 2026-06-01 e3nn sections: prior ring e3nn baseline and broad-backbone
  corrected-axis context.

Substrates it reads:
- Ring e3nn consumes an extractor output directory containing:
  `ring_current_sources.csv` and required target sidecar
  `rediscover_ring_current_sources_target_local_T2.npy`.
  `analysis/ENV.md` documents `/tmp/rdc-composed` and `/tmp/rediscover-rebuild-npy`
  as historical/candidate directories, but current existence was not grounded here.
- Broad e3nn consumes an extractor output directory containing:
  `broad_backbone_aggregated.csv`, `broad_backbone_sources.csv`, and required target
  sidecar `broad_backbone_aggregated_target_local_T2.npy`.
  The docs identify `/tmp/rdc-broad-backbone-axes` and `/tmp/rdc-broad-backbone` as
  historical/candidate directories, but current existence was not grounded here.
- If sidecars are absent, whether to re-extract is for lead+Claude - not specified here.

Exact command set, subject to lead go and sidecar verification:

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/analysis

NV=$HOME/.local/lib/python3.12/site-packages/nvidia
TORCHLIB=$HOME/.local/lib/python3.12/site-packages/torch/lib
export LD_LIBRARY_PATH="$NV/cu13/lib:$NV/cudnn/lib:$NV/nccl/lib:$NV/cusparselt/lib:$NV/nvshmem/lib:$TORCHLIB:$LD_LIBRARY_PATH"

/usr/bin/python3 test_change_of_basis.py

RING_DIR=/tmp/rdc-composed
ls \
  "$RING_DIR"/ring_current_sources.csv \
  "$RING_DIR"/rediscover_ring_current_sources_target_local_T2.npy
/usr/bin/python3 equiv_t2_e3nn.py "$RING_DIR" --cross both --epochs 4000 --lr 3e-3

BROAD_DIR=/tmp/rdc-broad-backbone-axes
ls \
  "$BROAD_DIR"/broad_backbone_aggregated.csv \
  "$BROAD_DIR"/broad_backbone_sources.csv \
  "$BROAD_DIR"/broad_backbone_aggregated_target_local_T2.npy
/usr/bin/python3 equiv_t2_backbone_e3nn.py "$BROAD_DIR" --with-axes --epochs 4000 --lr 3e-3
```

Note for future runner: `LD_LIBRARY_PATH` is copied from `analysis/ENV.md`, but the
future runner must verify the exact installed CUDA wheel directories before firing.
If any preflight `ls` fails, stop and ask for lead direction; do not silently re-extract.
If the clean rerun is intended to compare directly to the old leaky command line, the
old leaky logs/outputs must be identified and cited; their exact location was not
grounded in this assembly pass. That comparison is for lead+Claude - not specified here.

Held-out gate:
- Ring/broad frame-split gate must be blocked/purged, using `make_split_masks`.
- Train rows only are used for target centering and feature normalization.
- Forward-time prediction centering must pass `center_mask=g_tr`.
- Leave-atoms-out is reported or opt-in depending on script; it is not the primary
  gate for the e3nn clean protocol.

Anti-circular gate:
- Inputs are C++-emitted CSV/NPY substrates only.
- Target sidecars are required and row-aligned; no Python projection fallback.
- `change_of_basis` is the frozen basis map, not a rederived analysis projection.
- No DFT target or residual participates in source selection.

Disk gate:
- e3nn scripts print to stdout and do not write large model artifacts by default.
- If logs are captured, put them under a lead-approved `/tmp/rediscover-runs/...`
  run directory, not `/shared`.
- Do not re-extract missing sidecars without a separate lead-vetted run plan.

SETI gate:
- Record the printed frame-split T2 `R2`, `|T2| r`, effective N/stratum support,
  split strategy, purged frames, and `cross_split_lag1_pairs`.
- Do not convert clean-vs-leaky differences into a maths verdict. The verdict is
  for lead+Claude - not specified here.

Expected verdict:
- A clean-protocol number exists for ring and broad e3nn, replacing the previously
  suspect all-group de-mean/random-split number.
- Whether the clean number preserves, weakens, or kills the "three paths agree"
  statement is for lead+Claude - not specified here.

Output that proves it:
- `test_change_of_basis.py` passes under the same env.
- Ring log prints `split=blocked` behavior through the shared protocol and final
  `frame-split: T2 5-vec R2(test)=... |T2| corr(test) r=...`, plus the final gate line.
- Broad log prints per-stratum `split=blocked`, `purged_frames`, `cross_lag1_pairs`,
  per-stratum `T2 R2`, `|T2| r`, and final broad-backbone summary.
- Structural proof remains the postmortem: old random split calls and `g_tr = ~g_te`
  are gone from ring/broad paths; train-only centering and `center_mask=g_tr` are in use.

## Part 2 - Lead+Claude maths agenda

Do not answer any item in this part. Every item is for lead+Claude - not specified here.

Axis framing for the discussion, copied from the recorded dispositions and not a new maths
call: the within/frame axis is the current claim-bearing instrument; the between/atom
LOAO axis is retracted and may be discussed only as a case-study/diagnostic on 1P9J;
between / transferability remains deferred to the 720-WT instrument.

### Validity gates: could the combine be a mirage?

- Regularization plus effective DOF for the 26-term unified combine:
  for lead+Claude - not specified here.
- The right null for a combine, including structured null versus shuffle-target null:
  for lead+Claude - not specified here.
- Basis-invariance, including `change_of_basis` and per-type sums:
  for lead+Claude - not specified here.
- Held-out recovery, including which held-out axis is claim-bearing after true LOAO:
  for lead+Claude - not specified here.

Collinearity is not a validity gate in the assembled agenda. It describes the combine.
The local docs frame collinearity as expected/possibly confirmatory for
shadows-of-one-object: recovery is the robust column-space claim; per-shadow drop-one
attribution is the soft part. How to word or grade that is for lead+Claude - not
specified here.

### Attribution and reporting

- Per-shadow drop-one marginals such as field/McConnell drops are attribution, not
  the core recovery claim: for lead+Claude - not specified here.
- Report recovery as the robust claim and attribution as soft: for lead+Claude - not
  specified here.
- Statistical position plus determinability, not "survives/overfit": for lead+Claude -
  not specified here.
- Probability where N earns it; case-study where N is tiny: for lead+Claude - not
  specified here.

### Physics-architecture fold-ins

- Grade against the fixed-eigenstructure null: the parameter-free
  `(3cos^2 theta - 1)` angular shape with node at `cos^2 theta = 1/3`, not a generic
  fit family. How to apply it to the combine is for lead+Claude - not specified here.
- Define path agreement as coefficient agreement in physical units on a dominance-clean
  stratum, not "the three R2s are close." The equivariant leg of that
  coefficient-agreement gate depends on 1B actually running and producing the clean
  coefficient; the exact gate is for lead+Claude - not specified here.
- Maintain the convention ledger: `ringchi` opposite-sign/different convention; WaterField
  `+Hessian` versus APBS/MOPAC `-Hessian`; Larsen ppm tensors separate from geometric
  H-bond shadow; ring current excluded from the unified symmetric `D_ab` sum. Any
  mathematical consequence is for lead+Claude - not specified here.
- Tied angular form with per-mechanism radial/intensity is a fold-in from
  `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md`; whether to make it part of this
  maths cycle is for lead+Claude - not specified here.

### Standing constraints while discussing the maths

- Settled within-axis results stand: charge q/r^3 within `0.28/z451`, ring current
  within `0.28/z155`, unified combine within `0.43/z263`, with field+McConnell not
  charge-in-a-coat.
- Between/LOAO numbers are retracted, not provisional.
- 1P9J is the within instrument; the 720-WT owns between. Any future change to that
  closed disposition is outside this spec and for lead+Claude - not specified here.
- Do not quote any between/LOAO number as a result.

## Part 3 - Review-cycle protocol

Future loop checklist:

- Read first: `NOW.md`, the top Stage 2 block of `STATE.md`,
  `MATHS_AUDIT_CHECKLIST_2026-06-04.md`, `POSTMORTEM_MATHS_WALK_2026-06-04.md`,
  the four resolution postmortems, `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md`,
  and this spec.
- Prove comprehension to the lead before firing: settled within results, retracted
  between/LOAO numbers, #1 benign, #2 fixed/retracted, #3 code-fixed/re-run held,
  and the split between validity and attribution.
- Plan-vet risky loops before firing. The lead owns the go/no-go and all git.
- Preflight substrates and sidecars. If an expected `/tmp` substrate is missing, stop
  and ask for lead direction; do not silently re-extract.
- Preflight disk: lean `/tmp/rediscover-runs` budget, no result data under `/shared`,
  and no uncontrolled large artifacts.
- Preflight env for e3nn: system `/usr/bin/python3`, pinned torch/e3nn versions, and
  `LD_LIBRARY_PATH` per `analysis/ENV.md`.
- Run only the lead-approved command set. Capture stdout/stderr logs if a future
  postmortem needs exact printed gates.
- Preserve held-out discipline: true held-out atom means for LOAO; blocked/purged
  held-out frame gates for e3nn; train-only preprocessing everywhere.
- Preserve anti-circular discipline: selection from input-side conditions only; no
  DFT target/residual/score in selection predicates.
- Preserve SETI discipline: numbers, support, gates, and artifacts first; no
  editorial verdict from codex.
- After new physics-math or verdict runs, run adversarial review before any claim is
  promoted. The review must check: no maths call made by codex, no re-opened
  disposition, artifacts match commands, held-out/anti-circular/disk/SETI gates, and
  wording does not overclaim.
- Optional second adversarial pass over now-fixed code: inspect
  `analysis/stage2_law_fits.py`, `analysis/e3nn_protocol.py`,
  `analysis/equiv_t2_e3nn.py`, `analysis/equiv_t2_backbone_e3nn.py`, and the EFG
  template path. This is optional and for lead+Claude - not specified here.
- Commit only green, and only by the lead. Agents do not merge, switch, rebase, PR,
  reset, or commit.

## Index - decided dispositions and citations

Branch discipline, as recorded:
- Branch `h5-reader-pysr-spike`; never merge/switch/rebase/PR. Lead owns all git.
  Source: `CODEX_BRIEF_SPEC_MATHS_REVIEW_2026-06-04.md`, `NOW.md`,
  `NEXT_SESSION_PROMPT.md`.

Stage 2 settled within context, as recorded:
- Stage 2 `ecbddd1`: fitter decomposed, per-mechanism law fits, unified `D_ab` sum.
  Source: `NOW.md`, `STATE.md`.
- Stage 2.1 `f92109c`: happy-spot sweep and frame ablation. Source: `NOW.md`,
  `STATE.md`.
- Stage 2.2 `943915f`: unified vet, real combine not charge-in-a-coat on within.
  Source: `NOW.md`, `STATE.md`.
- Stage 2.3 `bp5cixi7k`: probability close. Later audit retracts LOAO/between
  readings but leaves within-axis charge/ring/unified standing. Source: `NOW.md`,
  `STATE.md`, `MATHS_AUDIT_CHECKLIST_2026-06-04.md`.

Maths walk:
- Maths-walk `bmgflkzi8`: read-only adversarial walk found three concrete issues.
  Source: `NOW.md`, `POSTMORTEM_MATHS_WALK_2026-06-04.md`.

Issue #1 disposition:
- `POSTMORTEM_DIAPARA_CHECK_2026-06-04.md`: dia+para=total holds at all emitted
  components to ORCA print-rounding; split targets sound; total target was not
  implicated. The checklist marks this `RESOLVED-BENIGN 2026-06-04`.

Issue #2 disposition:
- `POSTMORTEM_TRUE_LOAO_2026-06-04.md` / `b6e4d2e`, as recorded: true-LOAO run shows
  1P9J true between-atom recovery is approximately null; prior LOAO/between numbers
  were mislabeled within-modulation and are retracted; between / transferability /
  combine-depth defers to 720-WT.
- `POSTMORTEM_LOAO_FIX_2026-06-04.md` / code-fix landing `bh5f0e7ve`, as recorded:
  main LOAO/between path now uses shared `true_loao_*` atom-mean machinery; full
  re-fit was not run and was held per lead.

Issue #3 disposition:
- `POSTMORTEM_E3NN_PROTOCOL_FIX_2026-06-04.md` / `b804wd9rr`, as recorded: ring/broad
  e3nn aligned to clean EFG protocol through shared `analysis/e3nn_protocol.py`;
  blocked/purged split, train-only centering/normalization, and `center_mask=g_tr`
  are in place. Code-only fix; e3nn re-run not run; clean-vs-leaky verdict pending.

Physics-architecture fold-ins:
- `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md`: architecture stands; fold-ins are
  fixed-eigenstructure null, coefficient agreement in physical units, convention
  ledger, and tied angular/per-mechanism radial framing. All mathematical grading
  remains for lead+Claude - not specified here.

Current audit state:
- `NOW.md` records: audit code-complete; #1 benign, #2 and #3 fixed; parked verdict
  re-runs are LOAO re-fit and e3nn clean re-run; do not quote between/LOAO numbers.
