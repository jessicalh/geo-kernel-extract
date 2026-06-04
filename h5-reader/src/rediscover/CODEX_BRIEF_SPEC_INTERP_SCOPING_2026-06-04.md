# Codex brief — SCOPE + build-plan the 1P9J interpolation tool (the equivariant model, simply)

> **Historical session brief — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Status: **DRAFT — pending lead plan-vet, not yet fired.**

You own the grind; the lead vets + judges + owns ALL git. Deliverable: ONE spec that both SCOPES feasibility
and, if GO, DOUBLES as the build plan — so the moment the lead green-lights, the build runs from this same
doc (no second briefing round). **No build in this pass — spec only.**

## THE LEAD'S FRAMING (load-bearing — design to it)
- The interpolation tool, **even built badly, is a huge win**: it shows the advisor — on a GRAPH he can see —
  that we are *heading toward* the Stage-3 maths, even if we are not there yet. **The graph is the deliverable.**
- **It is (mostly) the equivariant model itself**, not a throwaway scalar fit. Build the real e3nn equivariant
  predictor (the Stage-3 architecture, in simple form) applied to **1P9J within-axis** — reuse/EXTEND the
  existing `equiv_t2_*e3nn*` exemplar; **never author a hand-rolled second model**
  (`feedback_no_python_physics_except_labeled_integrity_test`).
- **We are very close already.** Existing results are most of an interpolator: equivariant T2 modulation
  R²≈0.44 (held-out frames, 5-comp), |T2| r≈0.75, and the through-space leave-atoms-out k≈21 ppm·Å³ →
  R²≈0.62 (coupled). The lift is packaging + the graph, not new physics.
- **Be realistic.** It is a *simple, not great* interpolator: within-axis only; the trajectory is a geometry
  sampler (instantaneous map, not a process); correlate-not-match (`feedback_run_the_algorithm_get_a_cookie`,
  `feedback_correlate_not_match`). Label it honestly for the advisor: **direction, not destination.**

## CONTEXT (read first)
- Memory `project_rediscover_state` (the equiv-T2 + leave-atoms-out numbers above) + `project_stage3_equivariant_gnn`
  + `project_thesis_reporting_arc`.
- Top of `STATE.md` (STAGE 2 block — unified within R²≈0.43, field+McConnell combine) + `NOW.md`.
- The e3nn EXEMPLAR to extend (do NOT re-author the model): `analysis/equiv_t2_e3nn.py`,
  `analysis/equiv_t2_efg_e3nn.py`, `analysis/change_of_basis.py`, `analysis/e3nn_protocol.py` (the clean
  blocked/purged protocol — use it), `analysis/FINDINGS.md`; the decomposed fitter `analysis/allatom_fit_common.py`
  + `analysis/stage2_law_fits.py`.
- The live 1P9J substrate path (see `STATE.md` for the current `/tmp/rediscover-runs/...build*` dir).
- `SPEC_LEARNING_MODEL_2026-06-04.md` (this interpolator is its 1P9J v0 / down-payment) — keep them consistent.
- Reader display: skim `src/app` + `src/model` for any predicted-vs-DFT overlay surface (a possible Track-D
  crossover for the demo).

## WHAT THE SPEC MUST COVER
1. **What it is:** the simple equivariant predictor — T2-valued (σ_iso/T0 alongside) — interpolating 1P9J DFT
   shielding across the trajectory's sampled geometries, within-axis only. The honest caveats UP FRONT.
2. **Feasibility / closeness:** trace exactly what already exists vs the small delta to a runnable
   interpolator. What is reused (the exemplar net, the clean protocol, the substrate) vs added (the held-out
   1P9J interpolation harness, the graph). State the delta honestly — confirm "very close" or flag hidden work.
3. **The graph(s) — design them explicitly** (this is the advisor-facing payoff): predicted-vs-DFT held-out
   scatter (T2 recovery + σ_iso), the recovery R², and an honest "heading-toward-Stage-3" caption. Optionally
   the reader overlay (predicted-vs-DFT in the viewer) — scope as OPTIONAL, cross-ref the UI work.
4. **GO / NO-GO + build plan:** the bar is LOW (even a crude graph is a win) but the LABEL must be realistic.
   If GO, give the build sketch ready for a gated codex loop: held-out protocol, the R² + graph outputs,
   disk/SETI gates, what proves it works. Scope-only here; the build is a separate gated loop on the lead's go.

## After you: an opus agent adversarially reviews this spec
An opus agent will pressure-test the feasibility honesty (is it really close, or is there hidden work?),
guard the within-axis / direction-not-destination framing against overclaim, sharpen the graph design, and
confirm the GO/NO-GO bar + realistic label. Write it to survive that.

## HARD RULES
- **DOCS ONLY (this pass).** Write exactly ONE file:
  `/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/SPEC_INTERP_1P9J_SCOPING_2026-06-04.md`.
  Change NO code, run NO build / training / fit, run NO `git`. Everything else READ-ONLY.
- **Extend the e3nn exemplar; never a hand-rolled second model.** T2 never scalar-collapsed. Python consumes
  the emitted substrate only (never open `trajectory.h5`). Don't regrow the fitter pile — reuse it. **No
  protein/spatial model in Python — not even a secondary aggregate; spatial work + the protein model live in C++.**
- Realistic labeling is mandatory: within-axis, geometry-sampler, correlate-not-match, direction-not-destination.
- Branch `h5-reader-pysr-spike` — never merge/switch/rebase/PR. **Lead owns ALL git.** Truthful, cited, no overclaim.
