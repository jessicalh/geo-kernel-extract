# Codex brief — Stage 2: the law fits (finish what we started)

Status: **DRAFT — pending lead plan-vet, not yet fired.**

You own the grind; the lead vets + judges + owns ALL git. This is the PAYOFF: at the
dominance-clean audible places (Build 4's handles), fit each law-bearing kernel's law the
multiple ways we have, and say plainly — **is the kernel not-crap, what is the coefficient, its
statistical position, do the paths agree** — plus the NEW unified fit the physics-architecture
review handed us. Python analysis on the build4 substrate. Read first: `NOW.md`,
`PARTITION_FILTER_DESIGN.md`, `PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md`,
`analysis/FINDINGS.md`, and the existing fitter `analysis/allatom_fit_piece2.py` + the exemplar
template (`analysis/equiv_t2_*e3nn*.py`, `analysis/change_of_basis.py`).

## SCOPE
Substrate (READ-ONLY): `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4` — has
`dominant_fraction_{ring,charge,mc,field,hbond}` + C++ quintile bins + CaseHunter habitats +
the per-type shadow sums. Five law-bearing kernels: ring (Pople/J-B), charge (q/r³), McConnell
(bond Δχ), field (Buckingham, MOPAC-Coulomb), H-bond (Larsen). Target: **total-T2** (5-component,
frozen `get_C`; NEVER scalar-collapse).

## HARD RULES (the good rules)
- **Python ONLY fits the emitted substrate.** NO `trajectory.h5`, NO per-source dump, NO
  re-deriving any kernel/field/charge value in Python (rollback offense). frozen `get_C()`
  (`abs(C.TC−I).max()<1e-12`); score the 5-component T2, not `|T2|`.
- **DON'T GROW THE PYTHON PILE — the decompose REDUCES it.** Reuse the existing exemplar template
  (e3nn equivariant, PySR distill, ridge, change_of_basis); do NOT author a new equivariant model.
- **NO GIT** (run no git command). Leave gated-green uncommitted; lead commits. Branch
  `h5-reader-pysr-spike` — never merge/switch/rebase/PR.
- **DISK GUARD:** `df` before write, abort if free < 20 G; results to a fresh drop-old dir;
  deletes by EXPLICIT FULL PATH only (never a glob); total rediscover < 15 G; never `/shared`.
  KEEP build4 + build1.
- **Held-out everywhere** (between-LOAO atom-means, within-frameblock, AR(1) N_eff).
  **Anti-circular:** SELECT exemplars on INPUT-side dominance/geometry only (high
  `dominant_fraction_<mech>` + CaseHunter habitat); TEST on DFT total-T2 — NEVER select by DFT
  fit quality. **SETI:** numbers + curve shapes + bucket verdicts (recovered-law /
  form-recovered-scale-fitted / can't-make-it-work-for-now), support-flagged (atom-count/N_eff);
  never over-claim a thin stratum. **statistical position, never "law"-gate.**

## WRITE BUDGET — semi-deterministic; ESTIMATE before writing (size-gate)

The deliberate outputs (score tables, response curves, coefficients, cases manifests) are small and
bounded — anchors: Build 2 = 53 MB, Build 3 = 165 MB; stage 2 ≈ a few hundred MB. The ONLY path that
can balloon is PER-SOURCE / pairwise materialization (the boa-constrictor: the ~68 GB
all-atom-equivariant emit, the 1.7–5.5 GB backbone per-source CSVs). The equivariant path needs
per-source displacement vectors; build4 is `(atom,frame)`-REDUCED, so per-source must be materialized.
Therefore:
- **SCOPE per-source to the DOMINANCE-CLEAN EXEMPLARS ONLY** — the small selected subset (the point of
  the dominance gate). NEVER materialize per-source for the full 558,360 rows × all mechanisms; that is
  the tens-to-hundreds-of-GB path and is forbidden.
- **ESTIMATE before any per-source write:** print `n_exemplar_windows × sources/window × floats × 8 B`.
  If any single estimate trends past ~2 GB, STOP and report — do not write it. (Per-source for the
  clean exemplars should be MB-to-low-GB.)
- **Per-source intermediates are TRANSIENT** — drop-old by full path the moment the fit consumes them;
  only the reduced fit outputs (coefficients, tables, manifests) persist. Total transient + output stays
  inside rediscover < 15 G AND the 20 G-free abort.
- If a path genuinely cannot be scoped to the clean exemplars without a full-substrate per-source dump,
  STOP and report it as a design question — do NOT dump. (The equivariant per-source may need a small
  scoped C++ emit later; ridge + PySR on the already-emitted per-type reduced sums carry chunk 2 meanwhile.)

## Chunk 1 — DECOMPOSE the fitter (debt; gated by EXACT reproduction)
Split `allatom_fit_piece2.py` (~4.1k lines) into clean modules (substrate-read / basis+get_C /
partition+dominance / per-type ridge / fit-paths / reporting). **GATE: the decomposed pipeline
reproduces the Build-3 score table NUMBER-FOR-NUMBER** (per-type vs global-sliced; total/dia/para;
the dominance response curves). If ANY number moves, STOP — the refactor changed behaviour; do not
proceed to the fits. Only chunk-1-green unlocks chunks 2–3.

## Chunk 2 — per-mechanism dominance-clean law fits (the core)
For EACH of the five kernels, on its dominance-clean audible exemplars (input-side selected):
- Fit the kernel's expected law, **UN-SHRUNK**, on those exemplars.
- The PATHS where each applies: equivariant-T2 (e3nn), PySR (symbolic), ridge.
- Report per kernel: **(a) not-crap?** — does a real fit emerge anywhere (coefficient
  significantly non-zero on the clean stratum); **(b) the COEFFICIENT** in physical units;
  **(c) STATISTICAL POSITION** graded against the **fixed-eigenstructure NULL** — the
  parameter-free `(3cos²−1)` magic-angle shape (`cos²θ=1/3` node), not a generic fit family
  (`PHYSICS_ARCHITECTURE_UNIFICATION` fold-in); **(d) PATH-AGREEMENT in PHYSICAL UNITS** — do
  equivariant/PySR/ridge recover the SAME coefficient on the clean stratum (NOT "the R²s are close").
- **CONVENTION LEDGER** (fold-in): reconcile WaterField `+Hessian` vs APBS `−Hessian` sign before
  the field fit; treat ringχ as a different-convention path — NEVER naive-slope ringχ vs bs/hm
  (false validation-fail).

## Chunk 3 — the unified D_ab-sum fit (the NEW fourth path) + cross-kernel + examples
Every D_ab-family mechanism is ONE angular kernel `(3cos²−1)/rⁿ` on a different source-type. Fit
the UNIFIED form on the dominance-clean (through-space-dominated) atoms, on the EMITTED per-type
shadow sums (`charge_T2`, `mc_lit_T2`, field, hbond, pq/disp per-type) — NO per-source dump:
  `total-T2 ≈ Σ_type intensity_type · [per-type D_ab-shadow sum]`,
fitting per-type **INTENSITIES** (shared angular form across the D_ab family; per-type radial
stays — charge r⁻³, disp r⁻⁶).
- **RING IS SEPARATE:** keep bs/hm/jb as their own current-loop block (rank-≤1 `G=−V_a n_b`, NOT
  the symmetric D_ab) with its Pople intensity — do NOT fold ring into the shared-angular constraint.
- The local/bonded part is the **RESIDUAL** (AIMNet2's domain) — do NOT force the unified sum to
  explain it; the clean exemplars are where through-space dominates.
- Report: (a) does the unified sum recover total-T2 on the through-space-dominated atoms (held-out);
  (b) the recovered per-type **intensities + their statistical position (the fit-space described)**,
  and COMPARE them to literature as a REFERENCE only (sign, order of magnitude, ballpark) — NOT a
  pass/fail. With literature constants as inputs and a full numerical-methods stack (MD, DFT,
  kernels, fit) between, the recovered scale will NOT equal literature, and that is EXPECTED — the
  `form-recovered, scale-fitted` bucket, a real result. The deliverable is the RELATIONSHIP (does
  the unified form correlate / hold) and the fit-space, not a literature match (correlate-not-match;
  literature is sprinkles, not the bar); (c) how it compares to the per-mechanism fits + per-type ridge.
- **GUARD:** CONSTRAINED LINEAR fit — NOT a grand pointwise fit to total σ, NOT a single
  equivariant net over all calculators. Disagreeing shadows are DIAGNOSTIC, not averaged
  (correlate-not-match). Larsen ppm-tensor ≠ the geometric hbond shadow — do NOT double-count
  (don't treat the ppm tensor as another un-calibrated shadow).
- **CROSS-KERNEL:** the comparison table (each kernel's standing).
- **NAVIGABLE EXAMPLES:** extend the per-kernel `equations/<mech>/cases_manifest.csv` (the
  dominance-clean exemplars used) for the reader strips.

## Gates (ALL green per chunk; lead commits on green)
- Chunk 1: Build-3 numbers reproduced EXACTLY.
- frozen `get_C`; held-out everywhere; anti-circular (no DFT in selection); same charge-complete rows.
- Disk guard; lean result dir; never `/shared`.
- Integrity: input-acceptance, basis/targets, CV, partition.

## Resilience
Gate + bank each chunk before the next (chunk 1 must be green before 2–3). Short runs + retry if
the codex context-image bug bites. Results to a fresh drop-old dir
`/tmp/rediscover-runs/2026-06-04-stage2-fits/`.

## OUTPUT — TINY VERDICT
Write `src/rediscover/POSTMORTEM_STAGE2.md` (≤60 lines): the decompose module list + reproduce-gate
result; per-kernel (not-crap? coefficient? statistical position vs the physical null? paths agree in
physical units?); the unified D_ab-sum result (recovers through-space total-T2? intensities vs
literature?); the cross-kernel table; disk; run dir. Print ONLY a ≤12-line summary + that path.
DO NOT echo diffs/tables to stdout.
