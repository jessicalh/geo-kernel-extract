# Opus adversarial review — SPEC_MATHS_FIX_AND_REVIEW_2026-06-04.md

> **Historical review record — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Reviewer: opus adversarial pass. Scope: lane-check, Part-2 agenda PREP (no resolution),
Part-1 runnability sharpen, faithfulness. DOCS-ONLY; no code/git/build/fit run. Read-only
everywhere else. I make and resolve **no maths call** — same lane rule binds me.

Grounding limits of THIS pass:
- The Bash filesystem-listing call (`ls` of the `/tmp` substrate + e3nn candidate dirs) was
  DENIED by the sandbox here — the same restriction ENV.md documents (this env denies
  `import torch`/`numpy` and the compiler). So I could NOT live-ground `/tmp` existence either;
  I verified entry points, flags, gates, and cite-accuracy by reading the code and the
  postmortems. The spec's own "not grounded here" hedge (spec L16-17, L142-145) is therefore
  faithful, and ENV.md independently confirms the sidecar-presence risk (ENV.md L48-51).
- No git run; commit-hash-like tokens are taken as-recorded, not verified (matches spec L11-13).

---

## 1. LANE CHECK — does the spec CONSOLIDATE without DECIDING?

**Verdict: lane-clean (y).** I read the spec adversarially for any maths validity/attribution
call, any clean-vs-leaky verdict, any changed/re-opened disposition. I found none. Every open
item carries the explicit "for lead+Claude - not specified here" tag, and the count is high
(the phrase appears on essentially every unresolved line in Parts 1-2). Codex's claim that it
left all such calls to lead+Claude holds against the text.

Specific lane-discipline I checked and confirmed clean:
- Part 1A expected verdict (spec L97-105) is stated as **"confirmation-of-retraction, not a new
  between claim"** (L25-26, L93-94, L101-102) — it copies the already-recorded null, it does not
  re-derive or re-bless it. Correct lane.
- Part 1B expected verdict (spec L196-200) explicitly defers the clean-vs-leaky meaning:
  "Whether the clean number preserves, weakens, or kills 'three paths agree' is for lead+Claude"
  (L199-200). It does NOT assert leakage was real or that the number will move. Correct lane.
- Part 2 (spec L211-267) transcribes the validity/attribution split verbatim from the checklist
  and answers nothing. "Collinearity is not a validity gate ... It describes the combine"
  (L226-230) is a faithful transcription of the checklist's own framing (checklist L40-44), not a
  new judgment.
- Standing constraints (spec L259-267) restate settled-within / retracted-between as already
  decided; they introduce no new number.

**Two soft lane-adjacent notes (NOT drift, but worth the lead's eye):**
- (a) Spec L98-100 reproduces the retracted between numbers (`charge ~0.036, ring ~-1.0,
  unified negative/extreme overfit`) as the **held expectation for the re-run**. This is legitimate
  (it is the documented expectation a confirmation run must reproduce), and the spec fences it
  hard as not-a-result four times. But note the project rule "do not quote any between/LOAO number
  as a result" (spec L267): the numbers appear here only as *expected-diagnostic-output*, which is
  the one allowed use. I judge this in-lane, but flag it so the lead confirms the framing is the
  intended one. No change required from a lane standpoint.
- (b) Spec L97 calls the expectation "copied as the held expectation from the audit docs." That is
  the correct provenance; it is not the spec endorsing the number.

No line-referenced drift to report.

---

## 2. PART-2 AGENDA PREP (structure + evidence + options; RESOLVE NOTHING)

I answer none of these. For each open item I (i) sharpen the question, (ii) gather and cite the
local evidence already on disk, and (iii) lay out the option space — so the live lead+Claude
discussion starts one read in. Every item below remains **for lead+Claude**.

### VALIDITY — could the combine be a mirage?

**V1. Regularization + effective DOF over the 26 terms.**
- Sharpened question: with a ridge penalty, the *effective* degrees of freedom is
  `tr(X(XᵀX+λI)⁻¹Xᵀ)`, not the nominal 26 — so "26 terms / 25 atoms" overstates the free
  parameters at the operating λ. The question is: at the within-combine's selected α, what is the
  effective DOF, and is the within R²=0.43 reported at a fixed α or an inner-CV α?
- Local evidence: term count is exactly 26 (verified: `UNIFIED_SPECS`,
  `analysis/stage2_law_fits.py:139-148` — charge(1) + mc_category 0-4(5) + mopac backbone/aromatic(2)
  + water_efg reconciled(1) + hbond_geometric(1) + pq_type 0-7(8) + disp_type 0-7(8) = 26). The
  TRUE-LOAO atom-axis has N=25 atoms (`POSTMORTEM_TRUE_LOAO` table) — so on the BETWEEN axis 26 terms
  over 25 atom-means is p≥N (the postmortem's own `thin_TRUE_LOAO_p_ge_atoms` support flag).
  STATE.md L38 separately records the within combine "0.43 stable under shrinkage" and a "0.10
  shrunk" LOAO. The all-atoms-fit handoff (STATE.md L128-130) already flagged a **fixed-alpha
  under-regularization caveat** ("fixed alpha=10 on ~80 features ... likely under-regularized ...
  plausibly OUR regularization artifact") — directly relevant precedent that the lead should fold in.
- Option space (for lead+Claude, not chosen here): (a) report effective DOF at the operating α
  alongside nominal 26; (b) distinguish the WITHIN axis (many rows, ~3.9k) from the BETWEEN axis
  (25 atom-means, where p≥N is the real exposure) and apply the DOF concern only where it bites;
  (c) decide whether the within 0.43 needs an inner-CV-α restatement or is already on one.

**V2. The right NULL for a combine — structured vs shuffle-target.**
- Sharpened question: a shuffle-target null breaks ALL structure (it tests "is there any signal");
  a structured null (e.g. phase/block-preserving, or shuffling target atom-means across atoms while
  keeping features and folds) tests the claim-bearing thing — "is the combine's recovery more than
  what the column-space geometry alone buys." Which null does each axis need?
- Local evidence: the code already implements **axis-appropriate nulls** — Stage 2.3 "nulls shuffle
  the axis-appropriate target structure: within atom for frame modulation, across atoms for LOAO"
  (`POSTMORTEM_MATHS_WALK` S3, `stage2_law_fits.py:2085/2129/2171`). TRUE-LOAO uses "1000 deranged
  shuffles of target atom-means across atoms; structural features and LOAO folds unchanged"
  (`POSTMORTEM_TRUE_LOAO` L5). So a structured-ish null already exists for the between axis. The open
  question the checklist raises (checklist L38) is whether the COMBINE specifically wants a different
  null than the per-kernel ones.
- Option space: (a) keep the across-atom derangement as the between null and add a within-combine
  structured null; (b) consider a block/phase null on the within axis given AR(1) ρ≈0.53 (NOW.md L92,
  "effective ≪ raw") so the within null respects autocorrelation; (c) decide whether shuffle-target
  is reported as a floor (necessary-not-sufficient) with the structured null as the real gate.

**V3. Basis-invariance — change_of_basis / per-type sums.**
- Sharpened question: is the combine recovery invariant under the frozen library↔e3nn basis map
  (it should be, if it is a genuine column-space projection), and do the per-type partial sums
  behave consistently under that map?
- Local evidence (strong, already on disk): the map is "frozen, orthogonal, round-tripped, and
  Wigner-D tested against rotated Cartesian symmetric-traceless tensors" (`POSTMORTEM_MATHS_WALK` S4,
  `change_of_basis.py:162`, `test_change_of_basis.py:56/65`); the fit hard-fails if
  `|CᵀC − I|` is not tiny (`stage2_law_fits.py:722`), and the recorded value is
  `|CᵀC − I|max = 1.11e-16` (`POSTMORTEM_TRUE_LOAO` L3; spec L83). So orthogonality is a positive
  control already passing.
- Option space: (a) treat basis-invariance as already-established by the frozen/Wigner-tested map
  and the orthogonality gate, and merely cite it; (b) if a stronger demonstration is wanted, run the
  combine in both bases and show R² identity (a cheap confirmatory, not a new claim) — but whether
  that is worth a run is for lead+Claude.

**V4. Held-out recovery — which axis is claim-bearing after true LOAO.**
- Sharpened question: true-LOAO retired the between axis on 1P9J (it is ~null / case-study-N), so the
  claim-bearing held-out axis for the combine is now the **within/frame** held-out (blocked frame
  folds), not the atom-mean LOAO. The lead needs to confirm: is "held-out recovery" for the combine
  reported on the within/frame axis, with the atom-LOAO demoted to a documented case-study?
- Local evidence: within combine 0.43 (z263) is the within/frame held-out number (STATE.md L49);
  the between/LOAO is retracted to ~null (`POSTMORTEM_TRUE_LOAO`; NOW.md L60-64). The within evaluator
  uses blocked frame folds with purged rows (`POSTMORTEM_MATHS_WALK` S3,
  `stage2_law_fits.py:1463`, `allatom_fit_common.py:1423/1502`).
- Option space: (a) declare within/frame held-out as the claim-bearing axis for the combine and
  carry the atom-LOAO only as the retracted case-study; (b) defer any between-combine depth statement
  entirely to the 720-WT (the standing disposition, spec L266).

### ATTRIBUTION — describe, do not dismiss

**A1. Collinearity DESCRIBES, drop-one is SOFT.**
- Sharpened question: how to *word and grade* the per-shadow drop-one marginals so they read as a
  fingerprint of a shared source (confirmatory for shadows-of-one-object), not as instability that
  undercuts the combine.
- Local evidence (note a real tension the lead must reconcile): the drop-one marginals are NOT
  stable in sign/scale across cuts. Stage 2.2 within drop-one: "mopac_field +0.194, mc_cat4 +0.171,
  charge adds little" (STATE.md L37). Stage 2.3 within drop-one: "field +0.198" (STATE.md L53). But
  the TRUE-LOAO between drop-one is wildly inflated and overfit: "mopac_field_backbone +333.631;
  mc_category_4 +153.668; mc_category_0 +100.840; charge_total +52.528" (`POSTMORTEM_TRUE_LOAO` L12)
  — i.e. on the between axis the marginals are uninterpretable (p≥N overfit). So the "soft" framing
  must be axis-aware: drop-one is informative-but-soft on WITHIN, uninterpretable on the retracted
  BETWEEN. The checklist's framing ("recovery robust, attribution soft", checklist L42-44) is the
  right shape; the lead's wording should separate the two axes.
- Option space: (a) report drop-one only on the within axis, explicitly labeled soft/collinear; (b)
  show the collinearity structure directly (condition number / pairwise block correlation of the
  field & McConnell blocks) as the "fingerprint" rather than leaning on drop-one deltas; (c) state
  the one-object prediction (collinearity EXPECTED) so a reviewer does not read it as a defect
  (NOW.md L100-112 already argues this — the lead has the language).

### PHYSICS FOLD-INS

**P1. Fixed `(3cos²−1)` eigenstructure null.**
- Sharpened question: grade statistical position against the parameter-free D_ab angular shape
  (eigenvalues +2/ρ³,−1/ρ³,−1/ρ³; node at cos²θ=1/3), not a generic suggestable-fit family — i.e.
  the null/baseline is "the fixed magic-angle shape," and the question is whether the fit beats /
  matches THAT, with the free-5 fit as the alternative.
- Local evidence: the code already has the machinery — `fit_free5`/`free5_predict` and a
  `free5_beats_fixed_shape` / `fixed-shape` stat_position branch (`stage2_law_fits.py:279-292,
  545-563`). The physics doc says PySR already rediscovers `0.224≈1/3` as a positive control
  (`PHYSICS_ARCHITECTURE_UNIFICATION` L21-26, citing `FINDINGS.md:46`). So the fixed-shape-vs-free-5
  comparison is implemented; the open question is how to *apply it to the COMBINE* (the combine is a
  26-term sum, not a single angular law).
- Option space (for lead+Claude): (a) apply the fixed-eigenstructure null per-mechanism (where a
  single angular shape is meaningful) and report the combine separately; (b) define a combine-level
  "does it beat the per-mechanism fixed-shape baseline sum" comparison. Which is appropriate is a
  maths call — not made here.

**P2. Coefficient agreement in PHYSICAL UNITS (path agreement).**
- Sharpened question: "three paths agree" must mean the ridge/PySR/equivariant estimates of the SAME
  ppm-response coefficient agree on a dominance-clean stratum **in physical units**, not "the three
  R²s are close" (checklist L48, physics doc L27-30). The open question is the exact coefficient and
  the exact agreement tolerance on which stratum.
- Local evidence: per-kernel coefficients with CIs already exist — charge `+9.30 [8.10,10.51]`,
  field `−0.83 [−1.23,−0.43]`, ring `+0.69 [0.21,1.17]` (STATE.md L13-18). The convention ledger is
  emitted with the run (`stage2_law_fits.py:2867-2877`). The e3nn re-run (Part 1B) is the path whose
  clean coefficient is currently missing — so P2 is *gated on Part 1B actually running*. The lead
  should note: the coefficient-agreement gate cannot be evaluated for the equivariant path until the
  clean e3nn number exists.
- Option space: (a) define the agreement gate as overlapping CIs on the per-mechanism coefficient on
  a named dominance-clean stratum; (b) defer the equivariant leg of "three paths agree" until 1B
  produces the clean coefficient.

**P3. The convention ledger (the sign/convention ledger).**
- Sharpened question: confirm the three documented conventions are correctly applied in the fit so
  they do not resurface as false validation-fails: (i) ringχ opposite-sign/different convention →
  excluded from path agreement; (ii) WaterField `+Hessian` vs APBS/MOPAC `−Hessian` → water flipped
  before unified use; (iii) Larsen ppm tensors kept separate from the geometric H-bond shadow; ring
  current excluded from the symmetric D_ab sum.
- Local evidence (all confirmed in code, cite-accurate): water sign flip is live —
  `ShadowSpec("water_field_efg_reconciled", ..., "water_efg", sign=-1.0)`
  (`stage2_law_fits.py:144`); the producer-side convention reconciliation is explicit (APBS Hessian
  at `ApbsEfgTimeSeriesTrajectoryResult.cpp:160`, Water same form at `WaterFieldResult.cpp:125`,
  flip at `stage2_law_fits.py:145` per `POSTMORTEM_MATHS_WALK` S2). Ring excluded + Larsen separate
  is emitted in `convention_ledger.md` (`stage2_law_fits.py:2873-2876`) and asserted in S5 of the
  maths-walk. The ringχ-opposite finding is the data result jb/bs/hm agree +0.994 / ringχ opposite
  (physics doc L64-65; STATE.md cited as L82-87 in the checklist).
- Option space: this is largely *confirm-and-cite* rather than decide; the only live maths question
  is whether the WaterField/APBS sign reconciliation needs a positive-control check before the joint
  field fit (physics doc L31-34 says "reconcile before the joint field fit"). Whether to spend that
  check is for lead+Claude.

**P4. Tied angular form + per-mechanism radial.**
- This is flagged by the physics doc (L35-38) as a parameter-reduction fold-in (one M_ab angular
  machine, per-mechanism intensity). The spec correctly parks the question of *whether to make it
  part of this maths cycle* (spec L255-257). Local evidence: the e3nn exemplar already does
  shared-angular/per-type-radial (physics doc L13-17, "the strongest endorsement"). No PREP beyond
  noting it is a design fold-in, not a validity gate — keep it off the validity list.

**Agenda-prep usability assessment:** the prep above is one read from a live discussion — each item
has the question sharpened, the on-disk evidence cited to file:line, and the option space laid out,
with **nothing resolved**. The single most valuable cross-cutting observation for the lead: **almost
every validity/attribution concern is AXIS-DEPENDENT** — the WITHIN/frame axis is the claim-bearing
one (many rows, blocked folds, 0.43) and the BETWEEN/atom axis is retracted/case-study (25 atoms,
p≥N, overfit drop-ones). Framing the whole agenda as "within = the instrument, between = deferred to
720-WT" collapses several of the open items. (Stated as a framing the evidence supports; the maths
calls themselves remain for lead+Claude.)

---

## 3. PART-1 SHARPEN — are the two parked re-runs actually runnable?

### 1A. LOAO re-fit (true-between) — RUNNABLE as written, two sharpenings.

Entry point + flags VERIFIED against code:
- `--substrate-dir`, `--out-dir`, `--true-loao`, `--true-loao-permutations` all exist
  (`stage2_law_fits.py:153-160`). `--true-loao` returns early after writing only the
  `true_loao_*` artifacts (`main()` L2850-2854) — so the spec's diagnostic-command output list
  (spec L108-110) is **correct and complete**.
- The production-refresh command (default, no stage flag) falls through to L2855+ and writes
  `stage2_kernel_summary.csv`, `stage2_path_coefficients.csv`, `stage2_unified_dab_summary.csv`,
  `stage2_unified_dab_intensities.csv`, `stage2_run_audit.json`, plus `convention_ledger.md`
  (`main()` L2858-2879) — the spec's production output list (spec L111-114) is **correct**; it
  could also name `convention_ledger.md` explicitly (it says "plus the convention ledger" at
  L114, which is fine).
- The default run now genuinely uses the true-LOAO atom-mean machinery in the kernel/stage22/
  stage23 paths (`true_loao_single_feature_fit` at L501/911; `true_loao_atom_mean_arrays` +
  `true_loao_score_atom_means` at L647-648/960-961/1619-1620) — confirming the LOAO_FIX postmortem
  and the spec's "proof condition: output tables use true atom-mean LOAO" (spec L116). Good.

Gates VERIFIED:
- Disk guard requires ≥20 GiB free and <15 GiB in `/tmp/rediscover-runs` (`MIN_DISK_FREE_BYTES`
  / `MAX_REDISCOVER_BYTES`, L44-45; checks L194-202) — spec L86-87 is **exact**.
- `--true-loao` (and the stage flags) run with `drop_existing=True` (L2832); the bare production
  command runs with `drop_existing=False` and `mkdir(exist_ok=True)` (L203) — so a **fresh**
  production out-dir (the spec uses `2026-06-04-stage2-fits-loao-refit`, a new name) writes cleanly
  without deleting anything. The spec's "use a fresh lead-approved output directory unless the lead
  explicitly accepts drop-existing" (spec L88-89) is the correct and safe instruction.

Expected verdict HONEST: spec L97-105 copies charge≈0.036, ring≈−1.0, unified negative/extreme
overfit / null-class — matches `POSTMORTEM_TRUE_LOAO` (charge 0.036, ring −1.041, unified −104.579,
p0.698/z0.060) exactly. Framed as confirmation-of-retraction. Honest.

**Sharpening 1A-i (preflight the substrate):** the spec lists the Build4 sidecars (spec L37-42) but
the command set does not include an explicit preflight check. The disk guard does NOT check substrate
existence — `load_stage2_inputs` will fail at `load_json(substrate_dir / "per_atom_substrate_column_specs.json")`
(L689) if the dir is absent. Recommend the spec add a one-line preflight `ls` of
`per_atom_substrate_column_specs.json` + `per_atom_substrate_rows.csv` before firing (the spec's own
Part-3 checklist already says "preflight substrates," spec L281-282 — fold that into 1A's command
block so a cold runner does it inline). I could not `ls` it here (sandbox denial), so this preflight
matters more, not less.

**Sharpening 1A-ii (state the determinism/permutation seed):** `--true-loao-permutations 1000` matches
the floor (`n_perm = max(1000, ...)`, L2720). Fine as-is; optional nicety is to record that the null
is 1000 deranged across-atom shuffles (spec already implies this at L98-99 via the expected null
position). No change required.

### 1B. e3nn clean-protocol re-run — RUNNABLE in form, BUT gated on un-grounded sidecars.

Entry points + flags VERIFIED against code:
- `equiv_t2_e3nn.py`: positional `out_dir` (nargs="?"), `--cross both`, `--epochs 4000`, `--lr 3e-3`
  (L67-71) — the spec's `equiv_t2_e3nn.py "$RING_DIR" --cross both --epochs 4000 --lr 3e-3`
  (spec L159) is **valid**.
- `equiv_t2_backbone_e3nn.py`: positional `out_dir`, `--with-axes` (L128), `--epochs 4000`,
  `--lr 3e-3` (L137-138) — the spec's
  `equiv_t2_backbone_e3nn.py "$BROAD_DIR" --with-axes --epochs 4000 --lr 3e-3` (spec L162) is
  **valid**. Note the broad script's `--cross` defaults to `exact` (L123), and the spec passes none,
  so it uses `exact` — which matches the script's own docstring example (script L70). Consistent.
- `test_change_of_basis.py` is the correct gate per ENV.md (L32-38). The spec runs it first
  (spec L156), matching the documented order.
- Required inputs are fail-loud: ring needs `ring_current_sources.csv` +
  `rediscover_ring_current_sources_target_local_T2.npy` (`equiv_t2_e3nn.py:77/89-95`, `sys.exit` if
  absent); broad needs `broad_backbone_sources.csv` + `broad_backbone_aggregated.csv` +
  `broad_backbone_aggregated_target_local_T2.npy` (`equiv_t2_backbone_e3nn.py:169-184`, `sys.exit` if
  absent). The spec's required-file list (spec L135-145) is **cite-accurate**, and its anti-circular
  gate "Target sidecars are required and row-aligned; no Python projection fallback" (spec L181) is
  **correct** — the scripts have no fallback.

**The one real runnability gap (the spec already flags it, honestly):** the candidate dirs
`/tmp/rdc-composed`, `/tmp/rediscover-rebuild-npy` (ring) and `/tmp/rdc-broad-backbone-axes`,
`/tmp/rdc-broad-backbone` (broad) are **not grounded** — the spec says so plainly (spec L137-145,
L165-169), and ENV.md independently warns that the canonical `/tmp/rediscover-out-v2` has only CSVs,
that the NPY sidecar is REQUIRED, and that one must point at a dir that also has the sidecar or
re-extract (ENV.md L45-51). **I could not `ls` these dirs (sandbox denied the listing), so I cannot
upgrade the spec's hedge to a confirmation.** This is correctly handled: Part 1B is "runnable
*subject to sidecar verification*" (spec L147), and the spec defers "whether to re-extract if absent"
to lead+Claude (spec L146). That is the right disposition — do not let a runner silently re-extract
(NOW.md / the "no file discovery" rule).

**Sharpening 1B-i (make the preflight a hard inline step, not prose):** the spec's command block sets
`RING_DIR`/`BROAD_DIR` and fires; the sidecar-existence check lives only in the surrounding prose.
Recommend the spec promote it into the command block as an explicit `ls "$RING_DIR"/ring_current_sources.csv
"$RING_DIR"/rediscover_ring_current_sources_target_local_T2.npy` (and the broad triple) **before** the
python call, with "STOP and ask lead if absent" — so a cold runner cannot skip it. The scripts fail
loud anyway, but a pre-python `ls` distinguishes "wrong dir" from "torch import segfault" cleanly.

**Sharpening 1B-ii (the LD_LIBRARY_PATH honesty note is right; add the verify step inline):** the spec
copies the `LD_LIBRARY_PATH` from ENV.md and tells the runner to verify the installed CUDA wheel dirs
(spec L165-166). ENV.md confirms the exact same export (ENV.md L23-25) and that `import torch`
segfaults without it (ENV.md L17-21, memory `feedback_cuda_ld_path`). Good. Optional: have the runner
confirm `test_change_of_basis.py` passes (spec already requires this, L156/L204) as the de-facto
"torch imports cleanly" smoke before the two fits — which it does. No change needed beyond 1B-i.

**Sharpening 1B-iii (the clean-vs-leaky comparison target is ungrounded — the spec says so):** the spec
notes "if the clean rerun is intended to compare to the old leaky command line, the old leaky
logs/outputs must be identified and cited; their location was not grounded" (spec L167-169). This is
honest and correct — the comparison baseline is the gap. The lead should note: Part 1B can produce the
CLEAN number standalone, but the "clean-vs-leaky delta" needs the old leaky artifact located first.
That is correctly deferred (spec L169).

Expected verdict HONEST: spec L196-200 says only "a clean-protocol number exists ... whether it
preserves/weakens/kills 'three paths agree' is for lead+Claude." That is the correct non-committal
shape — it does not pre-judge that leakage was real. Honest.

---

## 4. FAITHFULNESS — is the disposition index cite-accurate?

I checked the spec's Index (spec L306-356) and Part-1 sources against the postmortems, the checklist,
NOW.md, and STATE.md. **Cite-accurate overall;** the numbers, dispositions, and file attributions
match. Specifics:

- Issue #1 (spec L328-331): "dia+para=total holds at all emitted components (T0/T1/T2) to ORCA
  print-rounding; RESOLVED-BENIGN" — matches `POSTMORTEM_DIAPARA_CHECK` (T0 max 1.0e-3, T2 max
  1.632993e-3 = isometric projection of 1e-3, inverse 3×3 max 1.0e-3; "split SOUND") and checklist
  L12-16. **Faithful.**
- Issue #2 (spec L333-340): true-LOAO null (charge 0.036, ring −1.0, unified −105) + retraction +
  720-WT deferral + code-fix-held — matches `POSTMORTEM_TRUE_LOAO` and `POSTMORTEM_LOAO_FIX` ("FULL
  RE-FIT NOT RUN — held per lead"). **Faithful.**
- Issue #3 (spec L342-346): ring/broad aligned to clean EFG protocol via shared `e3nn_protocol.py`,
  blocked/purged + train-only centering + `center_mask=g_tr`, code-only, re-run held — matches
  `POSTMORTEM_E3NN_PROTOCOL_FIX` line-for-line. **Faithful.**
- Settled-within (spec L259-264, L313-322): charge within 0.28/z451, ring 0.28/z155, unified 0.43/z263,
  field+McConnell not charge-in-a-coat — matches STATE.md L48-55 and NOW.md L59-79. **Faithful.**
- 26-term combine (spec L217): **independently verified** against `UNIFIED_SPECS`
  (`stage2_law_fits.py:139-148`) = exactly 26 terms. **Accurate.**
- `|CᵀC − I|max = 1.11e-16` (spec L83): matches `POSTMORTEM_TRUE_LOAO` L3. **Accurate.**

**Faithfulness flags (minor; none are spec errors, but the lead should know):**
- (f1) The spec, the postmortems, NOW.md, and STATE.md all carry **base-36-style commit-like tokens**
  (`bp5cixi7k`, `bmgflkzi8`, `b804wd9rr`, `bh5f0e7ve`, `bp5cixi7k`) alongside conventional short
  hashes (`ecbddd1`, `f92109c`, `943915f`, `b6e4d2e`). The spec explicitly says these are
  copied-not-verified (spec L11-13, L308-311) and that no git was run. I also ran no git. **So these
  identifiers are UNVERIFIED by both author and reviewer** — the spec is honest about this, but the
  lead is the only one who can confirm they resolve. Not a defect; a known open.
- (f2) The maths-walk's convention-ledger cites (`POSTMORTEM_MATHS_WALK` S5 at
  `stage2_law_fits.py:2580/2582/2583`) point at lines that now hold `true_loao_permutation_null`
  machinery, not ledger text — the ledger is now emitted at L2867-2877. This is **staleness inside
  the maths-walk postmortem**, NOT a spec error (the spec does not reproduce those line numbers). Flag
  it only so a future reader of the maths-walk is not confused; the substantive claim (ring/Larsen
  excluded) is still true and emitted.
- (f3) Spec L99/L101 says unified true-LOAO is "−105 ... with null-class position"; the postmortem
  has unified −104.579, pct 30.2, z0.060. Rounding to −105 and "null-class" (pct 30) is **faithful**.
  No overstatement.

I found **no misquote or overstatement that inflates a result.** Where the spec rounds, it rounds
conservatively; where it is uncertain, it says so.

---

## Bottom line

The spec is a faithful, lane-disciplined librarian pass. It consolidates without deciding, its two
parked re-runs are runnable (1A fully; 1B in form, gated on un-grounded sidecars it honestly flags),
and its disposition index is cite-accurate. The highest-value sharpenings are inlining the substrate/
sidecar preflights into the command blocks (both 1A and 1B), and giving the lead the cross-cutting
framing that the entire Part-2 agenda is axis-dependent (within = instrument, between = deferred). The
maths calls themselves remain, correctly, for lead+Claude.
