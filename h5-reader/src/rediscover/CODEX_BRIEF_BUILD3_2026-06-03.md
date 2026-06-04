# Codex brief — Build 3: fit-architecture loop (per-type fit alongside global + dia-T2 + dominance curve)

> **Historical run brief — not current truth (trued 2026-06-04).** Session
> provenance only; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Status: **landed** at `d35d7ec46e2c6b2064ab90228b6a997e4714e11a`; results are in `POSTMORTEM_BUILD3.md`.

You own the grind; the owner vets. Build 2 fit ONE global ridge across all 846 atoms and
SLICED it by stratum — which anti-predicts the H strata (HN between R² −54, robust CI, not
no-data) because one coefficient set is mis-applied across atom types with incompatible T2
physics. Build 3 fixes the fit to do BOTH — the all-atoms global fit (determinability) AND a
genuine per-type fit (Stage-1-style, each stratum its own coefficients) — and compares them.
Plus dia-T2 (channel due-diligence) and the dominance response curve. Python edge; extend the
committed `src/rediscover/analysis/allatom_fit_piece2.py`. Input substrate:
`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1` (already has the `dom_*` scalars).

## DISCIPLINES (hard)

- **Minimize Python — read the reduced substrate ONLY** (columns/sidecars/bin-ids + the `dom_*`
  scalar conditioners). No trajectory.h5 / per-source / older dirs. No physics re-derivation.
  Per-type fit, dia-T2, and the dominance binning are STATISTICS on already-C++-computed
  columns — allowed; do not re-derive any kernel/field/charge value.
- **DISK GUARD:** `df` before write; abort+report if free < 20 G; lean ~tens-of-MB result dir,
  drop-old by explicit FULL PATH only (never regex); total rediscover < 15 G; never `/shared`.
- Held-out everywhere, train-only PCA, same charge-complete rows, frozen `get_C()`,
  anti-circular (bins from INPUT-side only), SETI (numbers + shapes, NO verdicts).

## A. PER-TYPE fit alongside the global fit — the "both" (the #1 thing)

Add a PER-TYPE fit: for each stratum, a SEPARATE ridge fit on ONLY that stratum's atoms —
between = LOAO within the stratum's atoms (atom-means), within = frame-block on the stratum's
rows, alpha-selected per stratum. Report each stratum's between/within held-out R² from the
per-type fit **ALONGSIDE** the existing global-sliced R² (Build 2's), in the same score table,
with a `fit_scope` column (`global_sliced` vs `per_type`) and a per-stratum delta.

Headline questions to answer (as numbers, not verdicts):
- Does per-type fitting **rescue or clarify the H strata** — does HN go from −54 (global,
  anti-prediction) to recovery-or-honest-near-zero (per-type)? Same for HA/aromatic_H/polar_H.
- Do the strong strata (N, aromatic_heavy, O) hold under per-type?
- What does per-type's thin between-axis (~54 atoms/stratum) cost vs global's 846 — i.e., the
  determinability tradeoff. **Support-flag the thin per-type between-axis (atom-count/N_eff).**

(The within axis — many frames per atom — is the cleaner per-type comparison; the H rescue
should be most visible there.) Note for a future loop: the eventual unification is one model
with type×feature interactions (846-atom determinability + per-type coefficients); this loop
does the straightforward literal both — global + per-type side-by-side — to settle the question.

## B. dia-T2 target — the channel due-diligence

Add `target_dia_T2` as a fourth target (total-T2 / dia-T2 / para-T2 / T1), fit under both
scopes. Build 2 showed para-T2 recovers ≈0 while total-T2 N is 0.63. dia-T2 tells us whether
the recoverable signal lives in the diamagnetic channel (where through-space neighbor effects
would sit). Report dia-T2 per-stratum alongside total/para. (Note honestly: dia/para split is
gauge-origin-dependent — report, don't over-interpret the mechanism.)

## C. Dominance response curve

Bin the EXISTING `dom_*` scalar conditioners (`dominant_fraction_{ring,charge,mc}`) by quantile
(a lean Python quantile on the already-C++-computed scalar — NOT physics-derivation) and add
the dominance response curves to the partition (recovery vs dominance bin, per mechanism ×
stratum × target). **FLAG in the postmortem:** for consistency the `dom_*` bin-id belongs in
C++ alongside the other bin-ids — it should ride the next substrate re-emit; we bin the C++
scalar in Python here only to avoid re-emitting 3 GB for one column.

## Artifacts + output

Extend the score table with `fit_scope` (global_sliced / per_type) + the per-stratum deltas;
add dia-T2 rows; add dominance response curves to the partition outputs. Results to a fresh
drop-old dir `/tmp/rediscover-runs/2026-06-03-build3-fit-arch/`. Commit the SCRIPT atomically
on `h5-reader-pysr-spike` — NEVER merge/switch. Heavy result data not committed.

OUTPUT — TINY VERDICT: write `src/rediscover/POSTMORTEM_BUILD3.md` (≤55 lines): the
per-type-vs-global per-stratum table for total-T2 (esp. **did per-type rescue HN/H?** — the
headline); the dia-T2 vs total/para channel read per stratum; the dominance response-curve
shapes; the determinability cost of per-type (thin-between support flags); fit-stage +
partition-integrity checks; disk-guard; commit hash + run dir. Print ONLY a ≤12-line summary +
that path; DO NOT echo diffs/tables to stdout.
