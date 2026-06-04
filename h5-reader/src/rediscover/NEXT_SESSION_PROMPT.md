# Continuation prompt — rediscover (next session)

Paste to a fresh session. Branch `h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR. The LEAD owns ALL git.

You are resuming the nmr-shielding "rediscover" work. The MATHS AUDIT is **CODE-COMPLETE** (2026-06-04): the
numbers were chased, three concrete issues found + dispositioned, the docs made honest. The case study now stands
cleanly — **1P9J is a WITHIN instrument**: charge/ring/the combine recover on the within axis (the cookies); the
between axis was a mislabel (retracted → 720-WT). The spine now is the parked **VERDICT RE-RUNS** (confirm the
within numbers) + the secondary tracks — NOT re-doing the audit.

## STEP 0 — ORIENT (read, in order)
- `src/rediscover/INDEX.md` — the cold-start map + minimal reading path. Follow it.
- `NOW.md` (live marker) → `STATE.md` top-to-first-`---` (current synthesis; below the divider = archive) →
  `MATHS_AUDIT_CHECKLIST_2026-06-04.md` (the audit agenda) → `POSTMORTEM_MATHS_WALK_2026-06-04.md`.
- Memories: `feedback_law_as_statistical_position`, `feedback_transparent_cutoffs`,
  `feedback_applied_maths_over_methodology_caveats`, `feedback_correlate_not_match`,
  `feedback_run_the_algorithm_get_a_cookie`, `project_rediscover_state`, `project_unified_stats_engine`,
  `project_stage3_equivariant_gnn`, `feedback_token_economy_codex_codes`, `feedback_all_statistics_minimize_python_15gb`.

## STEP 1 — PROVE COMPREHENSION + THINK OUT LOUD (a gate, not a formality). Tell the lead back:
- **SETTLED & STANDING** (survived the audit): 1P9J = WITHIN instrument; charge q/r³ (within 0.28, z451) + ring
  current (within 0.28, z155) indicative+determinable; the unified combine's WITHIN recovery (0.43, z263) is a
  REAL combine (MOPAC-field + McConnell, charge≈0), not charge-in-a-coat. field/mc/hbond standalone ~null.
- **RESOLVED by the audit:** #1 dia/para split benign; **#2 the "LOAO" was within-modulation not between (true-LOAO:
  1P9J between ~null) → between numbers RETRACTED, code fixed, between defers ENTIRELY to the 720-WT**; #3 ring/broad
  e3nn leakage code-fixed to the clean EFG protocol. Two VERDICT RE-RUNS remain parked (STEP 2).
- **THE FRAME:** statistical position + determinability, never survives/overfit; probability where N earns it,
  case-study where tiny; 1P9J = WITHIN instrument; ~0.03 ≈ null, ~0.2 = something-or-trash by determinability;
  collinearity DESCRIBES the combine (expected for shadows-of-one-object), NOT a dismissal.
- FLAG anything off/stale. WAIT for the lead's go.

## STEP 2 — THE PARKED VERDICT RE-RUNS (the spine; on the lead's go)
The audit is code-complete; what remains are the two "are the numbers real, finally" confirmation runs (likely
one overnight batch): (a) the **LOAO re-fit** — now correct, though 1P9J's between is ~null regardless (720-WT
owns between); (b) the **e3nn clean-protocol re-run** — the actual clean-vs-leaky number for the combine's
"three paths agree," the one survivor the #3 fix leaves pending. Keep VALIDITY (DOF / null / basis-invariance /
held-out) separate from ATTRIBUTION (collinearity — describe, don't dismiss); fold in the physics-architecture
sharpeners (fixed `(3cos²−1)` null; coefficient-agreement-in-physical-units). Optionally a **second adversarial
pass** over the now-fixed code (lead floated it). Dispositions/maths = lead + Claude, NOT codex; the re-runs are
codex grind, gated, lead owns git. Closed dispositions are in `MATHS_AUDIT_CHECKLIST_2026-06-04.md`.

## STEP 3 — secondary tracks (on the lead's steer, AFTER the audit gates)
- **720-WT statics pilot, B-path** — now the ONLY between instrument: the durable static-pose/mutant-tree
  ingest (a `StaticPoseConformation` path in the rediscover extractor, 1P9J-trajectory oracle-parity gate).
  Uses the EXISTING 720 DFTs (same r²SCAN, absolute σ, `/shared/2026Thesis/consolidated/`, NO new ORCA). Also
  the seed of the unified stats engine ([[project_unified_stats_engine]]).
- **Stage 3 "ring-toss" predictor** — a proper DESCRIPTION of the equivariant-conditioned-GNN architecture, via
  agent + review ([[project_stage3_equivariant_gnn]]): feature pile (dihedrals+motion, Larsen DFTs untried,
  AIMNet2, MOPAC, OF3 embedding+OF3, geometry, all calculator outputs); train on DFT (the last stable ground),
  shoot at BMRB in the dark (the ensemble our short MD never samples), transfer on the 720-WT. Prediction not
  explanation; R² IS the metric; needs the chewer (per-source→GPU). Ground the physics FIRST, then toss rings.
- **ubiquitin-50fr** — cheap 2nd dynamic protein (1ns@20ps ≈ 50fr keeps the main signal; frame-ablation green).

## DISCIPLINE + CONSTRAINTS (load-bearing)
- codex grinds; the lead vets + judges; the lead's CONTEXT is irreplaceable — don't spend it on the
  build/gate/commit grind. Human-in-loop each loop; plan-vet risky loops BEFORE firing; adversarial review
  after new physics-math; gates are the truth, commit only green.
- **LEAD owns ALL git** (agents edit/build only — never commit/reset/rebase). Branch never merge/switch/PR.
- C++ filters over the in-memory indexes / minimize Python (the fitter is a decomposed ~4k-line pile — do not
  regrow it). Lean DISK < 15 G, generous RAM (swap OK). nmr_extract extractions SACRED. Never open
  trajectory.h5 in Python; frozen `get_C()`. Do NOT investigate the IUPAC-topology episode. ORCA = the one
  Trp-cage run only (not yet spent).
- **The maths discussion is lead + Claude, NOT codex.** Codex/agents for the build / gate / audit-walk grind only.
