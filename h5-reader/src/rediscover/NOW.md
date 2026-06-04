# NOW — the live process marker (rediscover, 1P9J)

**This is the MOVING MARKER: where we are *right now*. Update it every loop.** It is short by
design — `STATE.md` holds the detail/history; this is the "you are HERE" pointer. Branch
`h5-reader-pysr-spike` (NEVER merge/switch).

## The goal (unchanged)
Recover the expected CLASSICAL NMR shielding relationships from DFT on 1P9J as physics
**explanation** (not prediction); grade by **statistical position**. Arc: signal → equation
calibrations → ensemble → equivariant transferability pilot.

## The process we are in
The **partition arc**: fit-all (engine) → PARTITION by input-side conditions → favourable
clean exemplars EMERGE → laws fit on those → between-calculator network → equations table →
statistical-position grading.

## Loop ledger (this session, newest last)
- Loop 1 `4bb9a01` — charge-scalars emit (ff14sb + mopac welford). gated.
- Loop 2 `5b5525b` / 2b `44e786c` — all-atoms ridge fit + alpha-selection + AIMNet2 tier switch.
- Channel-completion `b583d7c` (Loop 3) + `9965c2b`/`1d0f56a`/`5a39288` (Loop 3b) — full
  cheap-wire menu + T0/T1/T2/dia/para targets.
- Build 1 `a353104` — C++ partition-filter tool (isolation primitives + bin-ids + CaseHunter).
- Build-1 fix `d9ba53d` — charge self/same-residue exclusion (H1) + real anti-circular guard
  (M1) + units (M2).
- Build 2 `a2d69cb` — the partition. N total-T2 0.63/0.81 (win); para-T2/T1/H-bond-lift negative;
  H strata ANTI-PREDICT (−54) — diagnosed as the **global-fit-sliced, NOT per-type** artifact.
- Build 3 `d35d7ec` — per-type fit alongside global + dia-T2 + dominance curves. **Per-type
  RESCUES H** (HN −54→+0.58 within; all H strata recover within 0.31–0.72) AND improves the
  heavy strata (O within 0.18→0.72, aromatic_heavy + N up). Per-type is the fit architecture.
- Doc truth-pass `43aa5bf` — grinder trued the docs (lead-committed; restore point `ad02f875`).
- Build 4 `fba5cd3` — dominance handle extended to the law-bearing FIVE: field (MOPAC-Coulomb) +
  H-bond get dominant_fraction/gap + C++ quintile bin-ids + CaseHunter habitats; charge widened
  (`charge_wide`); bins C++-side for all 5. Gated green (31/31 parity, 5e-13 two-path field/hbond,
  7/7 tests), build4 substrate (build1 kept). **C++ dominance arc complete.**
- Stage 2 `ecbddd1` — fitter DECOMPOSED (reproduce-gate byte-for-byte) + per-mechanism law fits +
  unified D_ab-sum. **charge + MOPAC-field = recovered laws; unified through-space R² 0.43;** ring
  real-but-thin (5 atoms); McConnell/H-bond standalone-null. Recorded in STATE.md (STAGE 2 block).
- Stage 2.1 `f92109c` — happy-spot sweep + frame-ablation: signal pops with MODULATION not isolation
  ("loud not isolated"); field+unified pop, nulls NOT rescued, support thin→720-WT; **ubiquitin-50fr keeps
  the signal** (cheap 2nd protein viable).
- Stage 2.2 `943915f` — unified vet: **REAL combine, NOT charge-in-a-coat** (carried by MOPAC-field +
  McConnell, charge ≈0); within stable, LOAO modest+N-limited (report by PROBABILITY, not overfit-dismiss).
- Stage 2.3 `bp5cixi7k` — probability close DONE (lead-verified): charge/ring/unified = indicative+determinable
  on WITHIN; between indeterminate→720-WT; **field standalone ~null** (0.03-class, top combine contributor);
  mc/hbond ~null. within cutoff-robust, between fragile.
- Maths-walk `bmgflkzi8` (read-only, overnight; `POSTMORTEM_MATHS_WALK_2026-06-04.md`) — adversarial walk
  nmr_extract→e3nn. Maths MOSTLY SOUND (strong SUPPORTS: producer formulae/conventions, the T2 basis +
  frozen/Wigner-tested get_C, train-only CV, axis-appropriate nulls, Larsen+ring excluded from the combine).
  **3 concrete issues for the MATHS TABLE — verify together, do NOT fix yet:** (1) DFT dia+para validated on
  T0 only, not T1/T2 (affects the dia/para SPLIT targets, not the total combine); (2) **the "LOAO" path
  centers by the held-out atom's OWN mean → measures WITHIN-modulation, NOT between-atom recovery** — if real,
  the WITHIN results stand but every LOAO/BETWEEN number is a mislabeled within-number (reframes the
  between-axis story + the 720-WT-as-between framing; also reinforces "1P9J is a within instrument"). THE
  consequential one. (3) ring/broad e3nn: all-group de-mean + unpurged random frame splits → possible leakage
  (the EFG e3nn path is clean = the template to align to).

## Current step
**THE GATE — maths audit IN PROGRESS: #1 + #2 closed, #3 open.** The WITHIN-axis results STAND (the audit does
NOT threaten them): charge q/r³ (within 0.28, z451), ring current (within 0.28, z155), the unified combine's
WITHIN recovery (0.43, z263, field+McConnell, not charge-in-a-coat). **#2 RESOLVED (true-LOAO `b6e4d2e`): 1P9J's
TRUE between-atom recovery is ~null — charge 0.036 (0.03-class), ring −1.0, unified −105 (overfit, null p0.70).
Every prior "LOAO/between" number (charge 0.38 / ring 0.17 / unified 0.26) was mislabeled within-modulation →
RETRACTED. Between / transferability / combine-DEPTH defers ENTIRELY to the 720-WT (1P9J has NO clean between
axis).** Code fix landing (`bh5f0e7ve`): the LOAO/between path → true between-atom recovery; full re-fit HELD per
lead. #1 RESOLVED-BENIGN (dia/para split sound). **#3 OPEN:** ring/broad e3nn leakage — bears on the WITHIN
combine's "3 paths agree," the next audit item. Do NOT quote any between/LOAO number.

**Stage 2 + ALL follow-ups landed (`ecbddd1` → `f92109c` → `943915f` → 2.3 probability close).**
RESULTS (probability + determinability; **1P9J = WITHIN axis ONLY — true-LOAO killed the between axis here**):
**charge q/r³ = indicative+determinable on WITHIN** (0.28, z451); **ring current = indicative+determinable on
WITHIN** (0.28, z155); **unified D_ab-sum = a REAL combine on WITHIN** — NOT charge-in-a-coat, carried by
**MOPAC-field + McConnell** (charge≈0), within 0.43 z263. **field standalone = ~null** (0.03-class; coeff
nonzero, recovery null-class; its value is the TOP combine contributor, not a standalone law). **mc + hbond
standalone = ~null.** **ALL between/LOAO numbers RETRACTED** (true-LOAO: 1P9J between ~null → the 720-WT owns
between + the combine's depth). Happy-spot: signal pops with **modulation, not isolation** ("loud not
isolated"); nulls not rescued. Frame: recovery survives ~50 frames → ubiquitin-50fr viable.

**REPORTING STANDARD (now law):** statistical position + determinability, NEVER survives/overfit. Probability
ONLY where N earns it (the within/frame axis, or many-atom sets); **CASE STUDY where N is tiny** (few-atom
between, e.g. ring-5 — report the found fit honestly; squeeze probability from the within/frame axis where the
between can't). Scale: ~0.03≈nothing; ~0.2 = something-or-trash decided by determinability; higher = clearly
something. [[feedback_law_as_statistical_position]] [[feedback_transparent_cutoffs]]. **1P9J is a WITHIN
instrument** (one protein, deep trajectory → probability lives on the within/modulation axis); its
between/LOAO is thin-BY-CONSTRUCTION = case-study → the **720-WT is the between-axis probability instrument**
(many proteins). Lead 1P9J verdicts on WITHIN; treat every LOAO/between figure as suggestive-pending-720-WT.

NEXT (open, lead to steer): (0) **land + verify + commit the probability close** (`bp5cixi7k`) — read tiny-N
between p-values as case-study, within p-values as real probability. (1) **720-WT statics pilot, B-path** — the
durable static-pose/mutant-tree ingest = the cure for the thin between-axis (many atoms → between case-studies
become probabilities) AND the seed of the unified stats engine ([[project_unified_stats_engine]]); same r²SCAN
+ same `.out`, absolute σ present, ~1–2 day adapter. (2) **ubiquitin-50fr** — cheap 2nd dynamic protein, frame
result says viable. (3) **McConnell/H-bond → joint/ensemble fit** (their home). Detail: STATE.md STAGE 2 block
+ Stage 2.1–2.3 follow-ups.
SIDE-FRAME: `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md` (vindicates architecture; fold-ins applied in stage 2).

## Live decisions / through-lines (the marker's real content)
- **THE COMBINE, IF IT HOLDS, IS THE DEEPEST RESULT — bigger than any standalone** (lead, 2026-06-04): the
  calculators-as-shadows-of-one-object claim made EMPIRICAL — the unified D_ab-sum recovers (within 0.43)
  where the individual shadows can't isolate; it carries what no part does alone. **BUT confidence is PENDING
  a JOINT MATHS discussion — lead + Claude, together, NOT codex.** Agenda, splitting VALIDITY from
  ATTRIBUTION (lead's correction 2026-06-04): **VALIDITY gates** (could the combine be a mirage) =
  regularization + effective DOF (26 terms), the right NULL for a combine (structured vs shuffle-target),
  basis-invariance (change_of_basis / per-type sums), held-out recovery — **collinearity is NOT on this
  list.** **ATTRIBUTION / DESCRIBE** (not a dismissal): collinearity among the shadow blocks is EXPECTED +
  possibly CONFIRMATORY for shadows-of-one-object — recovery is the projection onto the spanned column space,
  INVARIANT to collinearity; collinearity only destabilizes the per-shadow drop-one marginals. So report the
  RECOVERY as the robust claim, DESCRIBE the collinearity (fingerprint of the shared source), treat
  per-shadow attribution as the soft part. Dismissing the combine for collinearity = dismissing it for
  behaving as the one-object hypothesis predicts.
  **Do NOT over-claim the combine until we have talked the maths through.** Held provisional, deepest-if-holds.
- **Laws come from dominance-isolated clean exemplars, fit per-mechanism, UN-SHRUNK.** The
  per-type fit is a separate MODEL layer that may shrink but never owns the laws. **Dominance is
  the gate** between fine-clean (law) and coarse-mixed (model) — it protects the exemplars from
  the coarse wash. Bring dominance in from the start, not after.
- **Fit must be per-type, not pooled-then-sliced** (Build 2's H anti-prediction = the bill for
  pooling incompatible atom-type physics; Stage-1's per-element lesson).
- **dia-T2 is due diligence** for "signal is in total" (para-T2 ≈0; gauge-dependent split).
- **Categorical-engine / deep per-residue-position type granularity = PARKED** (real thesis
  material, not now; protein-limited on 1P9J; the transferability pilot's payoff).
- Disciplines: C++ filters over the in-memory indexes / minimize Python (no gazelle); lean DISK
  < 15 GB, generous RAM (swap OK); ORCA = exactly one more run, on Trp-cage (carries the
  DFT-reference extras); plan-vet risky loops before firing; adversarial review after new
  physics-math; gates are the truth, commit only green; nmr_extract extractions SACRED.

## Pointers
`STATE.md` (full detail/history), `PARTITION_FILTER_DESIGN.md` (durable filter architecture), the named memories (`feedback_*`/`project_*`).
