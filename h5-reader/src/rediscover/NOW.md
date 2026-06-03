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

## Current step
**Build 3 landed.** Verdict: per-type is the right fit (rescues H, improves heavy);
per-type-WITHIN is trustworthy + where the wins are; per-type-BETWEEN is thin (40–69 atoms,
p≥atoms) → the static axis wants the hierarchical/interaction model (846 + per-type coeffs) or
more proteins. **Channel = total-T2** (dia/para each recover worse — gauge-dependent dilution;
total is gauge-invariant). **H-bond→HN reopened** (HN recoverable per-type). Dominance curves in.
NEXT: lead check-in (docs restore point) → doc-truth-pass grinder → then between-calculator
network + equations table, built on the per-type-within foundation, dominance-gated.

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
`STATE.md` (full detail/history), the partition-filter design (designer's plan — *not yet saved
as a doc; gap*), the named memories (`feedback_*`/`project_*`).
