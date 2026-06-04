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

## Current step
**Stage 2 landed (`ecbddd1`) — the finish.** Fitter decomposed (reproduce-gate byte-for-byte), then
per-mechanism law fits on the dominance-clean exemplars + the unified D_ab-sum. RESULTS: **charge q/r³ =
recovered-law** (9.30, LOAO 0.38, paths agree); **field/MOPAC-Coulomb = recovered-law** (−0.83, weak but
nonzero, PySR agrees — MOPAC>Amber vindicated); **ring = form-recovered-scale-fitted** (real but THIN, 5
atoms); **McConnell + H-bond = can't-make-it-work standalone** (CIs span 0; joint-fit territory). **Unified
D_ab-sum recovers through-space total-T2** (within 0.43 / LOAO 0.26) — "calculators as shadows" combine
fits; intensities real but NOT literature-clean (correlate-not-match held). 3 paths: ridge + PySR +
equivariant-Schur (full e3nn-per-source still the deferred chewer). ONE protein, thin strata → nulls PROVISIONAL.
NEXT (open, lead to steer): (1) **happy-spot sweep** — response curves recovery-vs-cleanliness (dominance /
isolation / + a geometric-noise axis we lack — CaseHunter gates isolation/motion/quiet, none = low geom
noise), strict→loose: does it POP toward clean (tests "noisy geometry limits visibility")? rescue the nulls?
cheap, existing substrate. (2) **720-WT statics pilot** — same r²SCAN + same `.out` files (CONFIRMED,
absolute σ present); lots of rings → fattens ring's thin between-axis + cross-protein charge; needs the
substrate emitted on the 720 WTs (bounded static run). (3) frame-count ablation (cheap, within-axis).
(4) McConnell/H-bond → joint/ensemble fit. Full detail: STATE.md STAGE 2 block.
SIDE-FRAME: `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md` (vindicates architecture; fold-ins applied in stage 2).

## Live decisions / through-lines (the marker's real content)
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
