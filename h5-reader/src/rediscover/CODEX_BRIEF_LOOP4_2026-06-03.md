# Codex brief — Loop 4: refit on the COMPLETE substrate + free ring-current validation + field divergence

You own the grind; the owner vets. Loop 4 is the stats edge on the now-complete substrate
(`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-piece3b-final`, ~3 GB, 558,360 rows).
Extend the committed fit script `src/rediscover/analysis/allatom_fit_piece2.py` (the
alpha-selection + tier-switch version, `44e786c`). This is Python at the EDGE — the heavy
reduction is already done in the C++ substrate; Python only fits / correlates / partitions.

## DISCIPLINES (hard)

- **Minimize Python; read the reduced substrate ONLY.** Read the emitted columns/sidecars/
  column_specs from the piece3b-final directory. Do NOT open `trajectory.h5`, ORCA, older
  run dirs, or any per-source/pair data. Do NOT re-derive or re-aggregate physics in Python —
  the C++ spine already reduced everything to (atom,frame). Python does the fit (ridge/CV),
  the cross-method correlations, and the partition slicing — nothing heavier. No gazelle.
- **Lean outputs**, drop-old, total rediscover-runs stays < 15 GB.
- **Held-out everywhere** (between-LOAO over atom means, within-frameblock purged, AR(1)
  N_eff), train-only PCA on the embedding, same charge-complete row set across tiers,
  frozen `get_C()`. Anti-circular partitions (input-side bins only). SETI — numbers + curve
  shapes, NO verdicts.

## A. Targets — fit three, primary stays total-T2

1. **total-T2** (the thesis target; as before).
2. **para-T2** (`target_para_T2`) — the paramagnetic part, where the chemistry lives; dia is
   near-additive-atomic. Expect a cleaner per-mechanism story than total. Report side-by-side
   with total-T2.
3. **T1** (`target_T1`) — the field-linear DIAGNOSTIC (NMR-silent, not a shift claim). Fit it
   SPECIFICALLY against the vector field/charge mechanisms (apbs_E, ff14sb_field,
   mopac_coulomb_E, aimnet2 CRG): does the linear-field response show up in T1 where T2 (the
   E² / through-space home) misses it? Report as a diagnostic, labeled.

## B. Tiers — measure the new mechanisms' lift

Per target, predeclared tiers, scored per-stratum (N/CA/C/O/HN/HA + sidechain strata):
`classical` (the original set) → `classical_plus_newmech` (add Larsen per-class hbond, pq,
disp, the new field/EFG sources) → `plus_AIMNet2` → `all`. **Headline question: does adding
the new classical mechanisms — H-bond especially — lift the under-explained strata (HN, O)?**
Report `delta_R2_vs_classical` per stratum. Conditioning columns stay slicers, not features.

## C. Free ring-current cross-method validation (NO DFT, all frames)

The clean comparable set is **bs / hm / jb** — emit their pairwise correlation + slope on
`*_T2` (expect ≈1; report the structured divergence — near-field, fused rings — as the
physics). **ringχ is OPPOSITE convention** (anti-correlates bs/hm at −0.72; bare
susceptibility χ/r³, no shielding minus-sign, Å⁻³ not ppm·T/nA): sign-flip AND rescale it
before any comparison, OR report it as a separate-convention path — NEVER naive-slope ringχ
vs bs/hm (it would fake a disagreement). State which you did.

## D. Field divergence study (NOT an agreement test)

FF14SB / APBS / MOPAC / AIMNet2 / water field+EFG diverge by construction (vacuum vs
solvated vs polarizable). Do NOT score them as a slope≈1 agreement. Instead: (i) the CLEAN
within-method definitional check — confirm `EFG ≈ ∇(field)` and `field ≈ Coulomb-sum(charges)`
within a single method (a break = pipeline bug); (ii) report the cross-method divergence
(screening/polarization) as a structured study, not a pass/fail.

## E. DSSP donor/acceptor — select by atom role

The DSSP raw backup broadcasts donor AND acceptor records onto all four peptide-plane atoms.
Select by role: donor record → HN/N atoms, acceptor record → O/C′ atoms (use the emitted
`dssp_chem_donor_flag`/`acceptor_flag` to gate). Use as a partition conditioner.

## Artifacts + output

Score tables (per target × tier × stratum, held-out R² + CIs + deltas), the ring-current
validation table, the field definitional+divergence report, partition response curves +
favourable-case shortlist. Results to a fresh drop-old dir
`/tmp/rediscover-runs/2026-06-03-loop4-refit/`. Commit the script (analysis code) atomically
on `h5-reader-pysr-spike` — NEVER merge/switch. Heavy result data not committed.

OUTPUT — TINY VERDICT: write `src/rediscover/POSTMORTEM_LOOP4.md` (≤60 lines): per-target ×
per-stratum headline R² + the new-mechanism lift on HN/O; the bs/hm/jb agreement + how ringχ
was handled; the field definitional-check result; the T1-vs-field-mechanism diagnostic
result; emergent favourable partitions; commit hash + run dir. Print ONLY a ≤12-line summary
+ that path; DO NOT echo diffs/tables to stdout.
