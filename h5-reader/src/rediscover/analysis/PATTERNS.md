# Rediscover analysis — controlling patterns (the Python consumer side)

Read this FIRST, and front-load it in every analysis/fitter agent brief. This is
the constitution for the Python consumer surface, the analogue of the C++ library's
`PATTERNS.md` / `OBJECT_MODEL.md` / `CONSTITUTION.md`. It is the *push* that keeps
this surface in a good basin; once it's second nature the work is just productive.
The `#30` tripwire test is the backstop that catches a drift — not the lead.

## The boundary (load-bearing)

1. **The model is the spine, and it lives in C++.** The C++ producer computes AND
   EMITS the physics — kernels, fields, sums, the spherical (T0/T1/T2) decomposition,
   local frames, the DFT target. Python never recomputes any of it.
   (`feedback_model_is_spine`.)
2. **Python only fits + reads emitted features.** The fitter (ridge / e3nn / PySR)
   consumes the emitted CSV/NPY substrate: it READS emitted columns and NPYs; it does
   not re-derive physics. (`feedback_no_parallel_h5_in_python` — and never open
   `trajectory.h5` in Python; the reader owns H5.)
3. **Equivariant = e3nn.** The equivariant model uses real e3nn (`o3.spherical_harmonics`,
   tensor products, Wigner-D). NEVER hand-roll spherical harmonics or the projection in
   numpy/torch. (The deleted `equiv_t2.py`/`lib_T2` was the offense this doc exists to
   prevent.)
4. **One frozen change-of-basis.** The library↔e3nn `2e` map is the pinned constant in
   `change_of_basis.py` (`get_C()`). Reuse it; never re-derive it in a model path or
   inline a projection.

## The rules that follow

5. **NO recompute.** No `(3cos²−1)/r³`, no `q·d/r³`, no tensor projection, in any
   fitter/analysis path. The ONE exception: a **labeled pinned fixture test** (e.g.
   `test_change_of_basis.py`) — a one-shot equivalence assert on fixed inputs,
   justified by a separate source, never reusable recompute code.
   (`feedback_no_python_physics_except_labeled_integrity_test`.)
6. **Flag, don't hand-roll.** If a fit needs a quantity the substrate doesn't emit,
   FLAG it — the C++ emit extends (a codex task); you do NOT recompute it in Python.
   (The orientation-vector gap → an additive `BroadBackboneSink` emit, not a Python
   normal.)
7. **Report effective N; don't oversell.** Per-stratum, correlate-not-match; flag thin
   strata rather than force-fitting a number. (`feedback_seti`, `feedback_correlate_not_match`.)
8. **The exemplar is the template — debt in it propagates.** Catch technical debt
   where you are, in the moment; don't leave it for a cleanup that the velocity will
   outrun. (`feedback_catch_debt_in_the_moment`.)

## Why a doc and not just a test

Claude+codex generating code is stochastic descent on a high-dimensional landscape
full of local minima. A controlling patterns doc is the push that finds and holds a
good basin; the tripwire test only catches a fall after it starts. An agent that
reads this and works within it stays productive without anyone re-stating the rules.
That is the goal: pushy now, second nature soon.
