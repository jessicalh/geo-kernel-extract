# Continuation prompt — rediscover (next session)

Paste to a fresh session. Branch `h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR. The LEAD owns ALL git.

You are resuming the nmr-shielding "rediscover" work. The grounded physics + equivariant case study has neat
numbers for a student-scale exploration — but **"neat" is not "real" yet, and the spine of this session is the
careful MATHS AUDIT: are the numbers real?** Do not chase new results until it closes.

## STEP 0 — ORIENT (read, in order)
- `src/rediscover/INDEX.md` — the cold-start map + minimal reading path. Follow it.
- `NOW.md` (live marker) → `STATE.md` top-to-first-`---` (current synthesis; below the divider = archive) →
  `MATHS_AUDIT_CHECKLIST_2026-06-04.md` (the audit agenda) → `POSTMORTEM_MATHS_WALK_2026-06-04.md`.
- Memories: `feedback_law_as_statistical_position`, `feedback_transparent_cutoffs`,
  `feedback_applied_maths_over_methodology_caveats`, `feedback_correlate_not_match`,
  `feedback_run_the_algorithm_get_a_cookie`, `project_rediscover_state`, `project_unified_stats_engine`,
  `project_stage3_equivariant_gnn`, `feedback_token_economy_codex_codes`, `feedback_all_statistics_minimize_python_15gb`.

## STEP 1 — PROVE COMPREHENSION + THINK OUT LOUD (a gate, not a formality). Tell the lead back:
- **SETTLED** (the audit does NOT threaten): charge q/r³ (within 0.28, z451) + ring current (within 0.28,
  z155) indicative+determinable on the WITHIN axis; the unified combine's WITHIN recovery (0.43, z263) is a
  REAL combine carried by MOPAC-field + McConnell (charge≈0), not charge-in-a-coat.
- **PROVISIONAL** (the audit decides): every LOAO/between number + the combine's DEPTH/attribution claim.
  **CODE-CONFIRMED: the "LOAO" centers by the held-out atom's own mean → WITHIN-modulation, NOT between-atom
  recovery → 1P9J has NO clean between axis; the 720-WT is the ONLY between instrument.**
- **THE FRAME:** statistical position + determinability, never survives/overfit; probability where N earns it,
  case-study where tiny (ring-5); 1P9J = WITHIN instrument; ~0.03 ≈ null, ~0.2 = something-or-trash by
  determinability; collinearity DESCRIBES the combine (expected for shadows-of-one-object), NOT a dismissal.
- FLAG anything off/stale. WAIT for the lead's go.

## STEP 2 — THE MATHS AUDIT (the spine; once the lead confirms)
Walk `MATHS_AUDIT_CHECKLIST` **together — lead + Claude, NOT codex.** The 3 issues, **#2 first** (decide whether
a true LOAO/atom-MEAN test is even meaningful on one protein, or between is deferred entirely to the 720-WT);
then #1 (DFT dia/para T0-only validation) and #3 (ring/broad e3nn de-mean/unpurged-split leakage → align to the
clean EFG e3nn path). Keep VALIDITY (DOF / null / basis-invariance / held-out) separate from ATTRIBUTION
(collinearity-affected — describe, don't dismiss). Fold in the physics-architecture sharpeners (fixed
`(3cos²−1)` eigenstructure null; coefficient-agreement-in-physical-units). Decide what is REAL, what stays
provisional, what to fix. **A SECOND audit pass after**, if it surfaces more. Any fixes: lead-vetted, gated,
commit-only-green, lead owns git.

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
