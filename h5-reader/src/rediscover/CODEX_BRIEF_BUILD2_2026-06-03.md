# Codex brief — Build 2: the Python partition (the deliverable)

Status: **landed** at `a2d69cbb829dcc4901486c1162567fc99ba38b95`; retained as the execution brief.

You own the grind; the owner vets. Build 2 is the PARTITION — the product of the whole
effort: fit all atoms (the engine), then partition held-out recovery by input-side
conditions so the favourable cases EMERGE, per mechanism × stratum × target. Extend the
committed fitter `src/rediscover/analysis/allatom_fit_piece2.py` (the alpha-selection +
tier-switch version). Input substrate (the FIXED one): `/tmp/rediscover-runs/2026-06-03-
per-atom-substrate-build1` (charge axis repaired, anti-circular guard real).

## DISCIPLINES (hard)

- **Minimize Python — read the reduced substrate ONLY.** Read emitted columns/sidecars/
  column_specs + the C++ bin-id columns. Do NOT open trajectory.h5, ORCA, per-source, or
  older dirs. Do NOT re-derive/re-aggregate physics or re-compute bins in Python — the C++
  spine already reduced to (atom,frame) and emitted the bin-ids; bin by LOOKUP on those.
  Python does only: the ridge/CV solve, per-bin recovery on the bin-ids, plots. No gazelle.
- **DISK GUARD (insurance — Build 2 writes only ~tens of MB):** `df` before writing; abort+
  report if free < 20 G; write the result dir to a fresh drop-old dir; deletes by explicit
  FULL PATH only (never regex); total rediscover < 15 G; never write/delete under `/shared`.
  No redundant in-RAM copies of the big NPYs.
- **Held-out everywhere** (between-LOAO atom-means, within-frameblock purged, AR(1) N_eff),
  train-only PCA on the embedding, same charge-complete rows, frozen `get_C()`, anti-circular
  (bins from INPUT-side only — never DFT target/residual/coef). **SETI: numbers + curve
  shapes + the favourable shortlist; NO verdicts** (no "law" grading — that's a later step).

## A. Targets — fit three, primary stays total-T2

`total-T2` (primary), `para-T2` (cleaner — chemistry not diluted by the additive
diamagnetic part), `T1` (field-linear DIAGNOSTIC, NMR-silent — fit specifically vs the
vector field/charge mechanisms; report labeled, not a shift claim).

## B. Tiers — the new-mechanism lift, per stratum

Per target, per-stratum (N/CA/C/O/HN/HA + sidechain strata): `classical` →
`classical_plus_newmech` (add Larsen per-class hbond, pq, disp, the new field/EFG paths) →
`plus_AIMNet2` → `all`. **Headline: does adding the new classical mechanisms — H-bond
especially — lift the under-explained strata (HN, O)?** Report `delta_R2_vs_classical` per
stratum.

## C. THE PARTITION (the centerpiece)

Slice held-out recovery by the predeclared input-side condition families — using the
**C++ bin-id columns** (incl. the now-FIXED charge isolation: `dom_charge`, `gap_charge`,
`nearest_charge_r` are non-degenerate) — atom identity, geometry/isolation, driver
magnitude, driver modulation, support, charge-agreement. For each (mechanism × stratum ×
target × condition): **response curves** (recovery vs condition-bin), read the SHAPE
(monotone rise / fall / U / threshold / flat). The favourable cases EMERGE as the favourable
partition. This is categorical insight into the KIND of quantum effect a nearby-atom
environment creates — source-kind (mechanism × type) × response-kind (which of T2/para/T1
recovers).

**SUPPORT-FLAGGING (load-bearing — 1P9J is one protein):** every category/bin reports
**atom-count + N_eff**; fine categories thin between-atom are FLAGGED, never promoted on a
point estimate. Report the trend where the support is; flag where it's thin.

## D. Favourable cases ↔ the hunter

Intersect the high-recovery partitions with the **CaseHunter candidates** (the navigable
clean-habitat cases already emitted at `equations/<mechanism>/cases_manifest.csv`) → the
favourable-case shortlist (emergent-honest, not hunted). Emit the ranked shortlist as the
navigable manifest the reader strips will load.

## OUT OF SCOPE (deferred to step 2 — do NOT fold in)

The ring-current cross-method agreement (bs/hm/jb) and the field divergence study are the
between-calculator NETWORK = step 2, a separate loop. Build 2 is the partition only.

## Artifacts + output

Score tables (per target × tier × stratum, held-out R² + CIs + deltas), partition response
curves (per mechanism × stratum × target × condition, with atom-count/N_eff + shape),
favourable-case shortlist, run_audit. Results to a fresh drop-old dir
`/tmp/rediscover-runs/2026-06-03-build2-partition/`. Commit the analysis SCRIPT atomically on
`h5-reader-pysr-spike` — NEVER merge/switch. Heavy result data not committed.

OUTPUT — TINY VERDICT: write `src/rediscover/POSTMORTEM_BUILD2.md` (≤50 lines): per-target ×
per-stratum headline R² + the new-mechanism lift on HN/O; the partition response-curve shapes
per mechanism × condition (which rise / fall / flat, support-flagged); the emergent
favourable-case shortlist (the habitats); confirmation the charge partition is now
non-degenerate (multi-bin); fit-stage + partition-integrity checks; commit hash + run dir.
Print ONLY a ≤12-line summary + that path; DO NOT echo diffs/tables to stdout.
