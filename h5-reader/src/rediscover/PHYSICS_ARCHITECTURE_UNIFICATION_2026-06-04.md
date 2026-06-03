# Physics-architecture unification review (2026-06-04)

Review of `doc/calculators/physics-architecture.tex` against the rediscover fitting maths, by a
grounded agent (read NOW/STATE/PARTITION_FILTER_DESIGN/FINDINGS + the fitter / change-of-basis /
equivariant exemplar + `PerAtomSubstrate.cpp`, then the .tex). Captured as a frame for a later
session. **The architecture stands** — these are analysis-side fold-ins, not a redirect.

## Framing (lead-hand)

The .tex is not a new idea to chase — it is our own maths handed back. Its central claim (one
kernel `D_ab = ∂_a∂_b(1/ρ) = (3ρ̂ρ̂−δ)/ρ³`; calculators as projections onto different source
axes/lists) is literally what the emit computes (charge `q·D_ab` at `PerAtomSubstrate.cpp:646`;
McConnell/hbond `dipolar·(ĥĥ−I/3)` at `:808`; all → `DecomposeLibrary` T0/T1/T2 →
`change_of_basis`). It independently re-derives our single most important choice — **shared
angular form, per-type radial** (the e3nn exemplar) — the strongest endorsement in the document.
It does not change the architecture; it sharpens a few analysis-side things and hands us a better
null.

## FOLD INTO STAGE 2 (cheap, free, sharpen the deliverable)

- **Fixed-eigenstructure NULL (the best free win).** `D_ab` eigenvalues `+2/ρ³,−1/ρ³,−1/ρ³`,
  traceless, magic-angle node at `cos²θ=1/3` — a parameter-free angular prior. Grade statistical
  position against THIS fixed `(3cos²−1)` shape (a named physical null), not a generic
  suggestable-fit range. PySR already rediscovers `0.224≈1/3` (`FINDINGS.md:46`) — the node is
  real (positive control). [.tex §3]
- **Path-agreement in PHYSICAL UNITS.** equivariant/PySR/ridge all estimate the SAME response
  coefficient that turns the geometric shadow into ppm. The agreement gate is *coefficient
  agreement on a dominance-clean stratum, in physical units* — NOT "the three R²s are close."
  This is "the coefficient is the deliverable" made into the gate. [.tex §1, §6]
- **Convention ledger.** The .tex documents signs we've been bitten by: ring minus sign (= our
  ringχ-opposite finding, `STATE.md:82-87`), WaterField `+Hessian` vs APBS `−Hessian` (a real
  sign divergence between two field blocks). Reconcile WaterField/APBS sign before the joint
  field fit; treat ringχ as a different-convention path. Fold into analysis notes so they don't
  resurface as false validation-fails. [.tex §6, §7]
- **Tied angular, per-mechanism radial in ridge/PySR too.** McConnell/RingSusceptibility/HBond
  are one `M_ab` machine on different axes; fit ONE angular law with per-mechanism intensity,
  not five independent angular shapes. Parameter reduction, grounded (the e3nn path already does
  this). [.tex §5]

## DEFER (the one new structural idea — data-driven, not pre-emptive)

- **Axis-angle separation in the dominance gate.** Because ring/McConnell/hbond share the angular
  machine on different axes, magnitude-dominance alone may not separate their LAWS when two
  co-located sources share an axis direction (the tensor analogue of the differencing-lockstep
  problem, `FINDINGS.md:213`). Fix = a cheap C++ axis-angle separation scalar (same family as
  `gap_to_2nd_*`). DEFER: proceed on BUILD 4 magnitude-dominance; add this only if stage-2 fits
  show ring/mc/hbond contaminating each other. [.tex §5, §11]

## GUARD (disciplines, reinforced)

- Do NOT read "one hidden object" as "one grand pointwise fit to total σ." The .tex itself says
  correlate-not-match; disagreeing shadows are DIAGNOSTIC, not to be averaged. Conflicts with the
  per-type evidence + minority-drowning (`STATE.md:357`). [.tex §11]
- Larsen ppm-tensor vs geometric H-bond shadow: we carry BOTH under similar names (geometric
  `dipolar·(ĥĥ−I/3)` at `:808` vs the real Larsen per-class ppm tensors). The joint fit must not
  double-count them or treat the ppm tensor as another un-calibrated shadow. Labeling hazard. [.tex §1, §8]

## WRITE-UP NOTES (no action)

- The architecture PERMITS T1 (asymmetric `M_ab` term) but the DFT reference doesn't populate it
  through-space (T1≈0, `STATE.md:142`) — note it; the .tex slightly over-sells a channel the data
  leaves empty.
- BiotSavart/Haigh-Mallion are rank-≤1 outer products (`G=−V_a n_b`), NOT the symmetric-traceless
  `D_ab` shadow — a different object. The .tex doesn't pick which ring calculator is cleanest; our
  data did (jb/bs/hm agree +0.994, ringχ opposite).

## MSc realism

Every fold-in is analysis-side, on fits we already run; no new producer code, no grand model.
Explicitly do NOT attempt a single equivariant net over all calculators — the .tex's "one
network" is a picture; our per-mechanism dominance-clean fits are the deliverable-sized version.
