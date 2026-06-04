# Codex brief — SPEC the maths-fix + review cycle (consolidate, do NOT decide)

> **Historical session brief — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Status: **DRAFT — pending lead plan-vet, not yet fired.**

You own the grind; the lead vets + judges + owns ALL git. This is a **librarian/assembler** pass: stash
the maths-audit outcome + the parked re-runs + the review protocol into ONE executable spec the lead
shelves ~2 days and resumes from cold. **The maths discussion itself is lead + Claude, NOT codex** — your
job is to CONSOLIDATE what is already decided and STRUCTURE the mechanics, never to make or re-open a
maths call. Honest stakes: a faithful stash now means we resume the hard-maths chase losslessly instead
of re-walking the audit.

## CONTEXT (read first)
- `MATHS_AUDIT_CHECKLIST_2026-06-04.md` (the closed agenda + dispositions), `POSTMORTEM_MATHS_WALK_2026-06-04.md`
  (the 3 issues), and the four resolution postmortems: `POSTMORTEM_DIAPARA_CHECK_2026-06-04.md` (#1 benign),
  `POSTMORTEM_TRUE_LOAO_2026-06-04.md` + `POSTMORTEM_LOAO_FIX_2026-06-04.md` (#2), `POSTMORTEM_E3NN_PROTOCOL_FIX_2026-06-04.md` (#3).
- `NOW.md` (the COMBINE NOTE + "Current step") + top of `STATE.md` (the STAGE 2 block) — the settled within
  results the audit does NOT threaten, and the retracted between/LOAO numbers.
- `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md` — the fold-ins (fixed-eigenstructure null; coefficient
  agreement in physical units; the convention ledger).
- The code the re-runs touch (read, do not run): `analysis/stage2_law_fits.py` (the LOAO/true-between path,
  ~lines 657-665/1678-1682/2066 + the `true_loao_*` machinery), `analysis/e3nn_protocol.py` (the shared
  clean protocol), `analysis/equiv_t2_e3nn.py` + `analysis/equiv_t2_backbone_e3nn.py` (the now-aligned paths).
- Memories: `feedback_applied_maths_over_methodology_caveats`, `feedback_law_as_statistical_position`,
  `feedback_two_path_validation`, `feedback_adversarial_review_physics`, `feedback_plan_vet_before_firing_risky_loops`.

## WHAT THE SPEC MUST ASSEMBLE (three parts + an index)

**Part 1 — the two parked verdict re-runs, as codex-ready sub-briefs (mechanical, runnable later).**
For each: the exact script/entry point + flags, the substrate it reads, the held-out / anti-circular / disk
/ SETI gates, the expected verdict, and what output proves it.
- (a) **LOAO re-fit** on the corrected true-between path (reuses the validated `true_loao_*` machinery).
  Expected: 1P9J true between-atom recovery ~null (charge ~0.036, ring ~−1.0, unified overfit) — the 720-WT
  owns between. State it as confirmation-of-retraction, not a new claim.
- (b) **e3nn clean-protocol re-run** (the held one) — the actual clean-vs-leaky number for the combine's
  "three paths agree," now that ring/broad e3nn use the shared blocked/purged `e3nn_protocol.py`.

**Part 2 — the lead+Claude maths agenda (LIST it, do NOT resolve it; mark every item "for lead+Claude").**
Faithfully transcribe the open agenda, split as the checklist already splits it:
- VALIDITY gates (could the combine be a mirage): regularization + effective DOF (26 terms), the right NULL
  for a combine (structured vs shuffle-target), basis-invariance (`change_of_basis` / per-type sums),
  held-out recovery. **Collinearity is NOT a validity gate** — it DESCRIBES the combine (expected/confirmatory
  for shadows-of-one-object); report recovery as the robust claim, attribution as soft.
- The physics-architecture fold-ins: grade against the fixed `(3cos²−1)` eigenstructure null (node at
  cos²θ=1/3); path-agreement = COEFFICIENT agreement in PHYSICAL UNITS on a dominance-clean stratum, not
  "the three R²s are close"; the convention ledger (ringχ opposite-sign; WaterField +Hessian vs APBS −Hessian).
- Do NOT answer any of these. They are the live discussion the lead and Claude hold together.

**Part 3 — the review-cycle protocol.** Adversarial review after new physics-math; the optional SECOND
adversarial pass over the now-fixed code; gate discipline (commit only green; human-in-loop each loop;
plan-vet risky loops before firing; lead owns ALL git). Make it a checklist a future loop runs.

**Index — the decided dispositions, verbatim + cited** (commit hashes + postmortem files) so a cold resume
trusts the stash without re-deriving.

## After you: an opus agent adversarially reviews this spec
An opus agent will (1) verify you stayed in your lane — that NO maths call or disposition was made or
re-opened here — and flag anything that drifted; (2) PREP (structure + evidence-gather + frame the
options for) the Part-2 lead+Claude agenda WITHOUT resolving it; (3) sharpen the Part-1 sub-briefs.
Write the stash to survive that review: faithful, cited, lane-disciplined.

## HARD RULES
- **DOCS ONLY.** Write exactly ONE file:
  `/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/SPEC_MATHS_FIX_AND_REVIEW_2026-06-04.md`.
  Change NO code, run NO build / fit / e3nn / re-run, run NO `git`. Everything else READ-ONLY.
- **CONSOLIDATE, NEVER DECIDE.** Make no maths judgment; re-open no disposition; resolve no open question.
  Every unresolved item is explicitly "for lead+Claude — not specified here."
- Settled within results stand (charge within 0.28/z451; ring 0.28/z155; unified combine within 0.43/z263,
  field+McConnell not charge-in-a-coat). Between/LOAO numbers are RETRACTED, not provisional. Quote no
  between number as a result.
- Branch `h5-reader-pysr-spike` — never merge/switch/rebase/PR. **Lead owns ALL git.** Truthful, cited, no overclaim.
