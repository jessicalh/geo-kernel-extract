# Codex brief — Loop 2b: alpha-selection + analysis-side AIMNet2 switch

> **Historical run brief — not current truth (trued 2026-06-04).** Session
> provenance only; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Status: **landed** at `44e786c094c4f9cd1ae4ebe451454fc3e31406d6`; retained as the execution brief.

You own the grind; the owner vets. Loop 2b settles a confound from Loop 2: the
`+AIMNet2` and `all` tiers HURT the between-axis under fixed `alpha=10` (global
between +0.104 → +0.058 → +0.053). On ~80 features over only 846 atom-means, a fixed
alpha carried over from an older, lower-dimensional fit is likely under-regularized
— so the negative delta is plausibly OUR regularization choice, not a property of
AIMNet2. This loop makes the test fair. Two analysis-side changes to the committed
`src/rediscover/analysis/allatom_fit_piece2.py` (commit `5b5525b`). Keep it focused.

## HARD GUARDRAIL — read this first

AIMNet2 is **always-on in the producer and the emitted substrate** — non-negotiable,
hardware-required, the owner's standing rule across many sessions. The substrate
ALWAYS carries the AIMNet2 columns. Therefore:

- **DO NOT** add any `--no-aimnet2` flag, `WITH_GPU=OFF`/build option, aimnet2 test
  label, or any mechanism that avoids COMPUTING AIMNet2 anywhere in the producer,
  build, or test machinery. Do not touch the C++ / extractor / CMake / ctest at all.
- "Switchable" here means **only** analysis-side: the FIT's feature-tier selection —
  whether the Python fit *consumes* the always-present AIMNet2 columns. Express it as
  feature-tier / feature-block selection (the existing `classical` vs `plus_AIMNet2`
  vs `all` structure made into a clean controlled comparison), NOT as a flag named to
  disable AIMNet2. The substrate is untouched; only fit consumption toggles.

## Change 1 — train-only inner-CV alpha selection (settles the confound)

For each tier × held-out axis (between-LOAO and within-frameblock), select the ridge
alpha by **training-only inner CV** from a predeclared grid
`[0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5]`. The held-out TEST fold must NOT participate
in alpha selection (no leakage). Record the selected alpha per tier × axis in
`run_audit.json` and add it to the score table. This is the spec's sanctioned
alpha-selection path (`ALLATOM_FIT_SPEC_2026-06-03.md` "Ridge And Preprocessing").
Keep `alpha=10` runnable as a labelled baseline for continuity, but the primary
reported scores use selected alpha.

## Change 2 — AIMNet2 as a clean analysis-side feature-tier switch

Make the controlled comparison `classical` (AIMNet2 features OFF) vs `plus_AIMNet2`
(ON) a clean parameter of the analysis (e.g. a `--tiers` / feature-block selector),
so the owner can run the AIMNet2-on vs AIMNet2-off comparison directly under selected
alpha. Train-only PCA on the embedding stays when AIMNet2 is on. No producer flag, no
`--no-aimnet2`; the substrate always has the columns.

## Re-run + the decisive question

Re-run the fit + partition with selected alpha across the three tiers (same
charge-complete rows, train-only PCA, all Loop-2 disciplines). Report — as numbers,
not verdicts (SETI) — whether the `+AIMNet2` between-axis delta is still
negative/zero or flips POSITIVE once alpha is selected per tier. Report the selected
alphas alongside the deltas so the comparison is legible.

## Disciplines (unchanged from Loop 2)

- Read the emitted substrate directory ONLY
  (`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-charge-scalars-piece1-final`).
  NEVER open `trajectory.h5`, ORCA, older/MOPAC/all-atom-equivariant dirs, or pair
  dumps. Audit states this.
- No Python physics recompute. Frozen `get_C()` (`abs(C.TC - I).max() < 1e-12`).
  Score the 5-component 2e T2, not `|T2|`.
- Anti-circular partition (input-side bins only; never DFT target/residual/coef).
- Held-out everywhere; report held-out R² only; AR(1) N_eff per component.
- Fit-stage checks ALL pass, PLUS a new one: **alpha-selection integrity** — the
  held-out test fold did not choose alpha; selected alphas are recorded. If any check
  fails, STOP and report; do not present scores as valid.

## Environment / artifacts / commit

Analysis venv (`src/rediscover/analysis/venv`). Extend `allatom_fit_piece2.py` (or a
clearly-named sibling) — committed code. Results to a fresh drop-old dir, e.g.
`/tmp/rediscover-runs/2026-06-03-allatom-fit-piece2b-alpha/` (ephemeral; do not commit
heavy data). Branch `h5-reader-pysr-spike` — **NEVER merge/switch/rebase/PR/checkout**.
One atomic commit of the script change. Message: describe train-only alpha selection
+ the analysis-side AIMNet2 feature-tier switch.

## Report back — drift + postmortem (plain, no editorializing)

- What landed vs asked; any drift.
- Each fit-stage check (incl. alpha-selection integrity): pass/fail.
- The **selected alphas** per tier × axis.
- Per-tier × per-stratum between/within held-out R² under selected alpha, and the
  tier deltas — explicitly: did the `+AIMNet2` between delta flip sign vs Loop 2's
  fixed-alpha result?
- Commit hash; results-dir path.
