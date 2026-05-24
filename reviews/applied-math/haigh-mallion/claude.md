# haigh-mallion — claude review (readability focus)

- **Targets:** src/HaighMallionResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

# Adversarial readability review — HaighMallionResult

## Verdict

**HaighMallionResult.h** — Reads cleanly. The header doc-comment is the model for how this kind of physics file should open: it states the two-step construction (raw integral H → shielding kernel G), the symmetry/trace properties, units, the quadrature rule with citation, and the parallel to BiotSavart, all before any code. A chemist who knows ring-current theory can follow the through-line top to bottom on one pass.

**HaighMallionResult.cpp** — Mostly coherent and well-signposted; the banner comments above each static function carry the math (kernel form, subdivision logic, fan triangulation) so the story holds. The one real friction is `Compute`, where the genuine three-step physics (integrate → contract with normal → outer-product to rank-1) is interleaved with verbose `GeometryChoice::Record` bookkeeping blocks that visually dominate and break the narrative; the reader must skip past provenance scaffolding to find the four lines that do the physics. A documentation discrepancy on the G sign (header says `G = n_b·(H·n)_a`, cpp says `G = -n_b·V_a`) is the only thing that would stop a careful reader.

## 1. Coherent story / readability

- HaighMallionResult.cpp:234–344 — the physics (Steps 1–3 at 281–294) is buried among three multi-line `choices.Record(...)` lambda blocks (239–246, 258–266, 270–279); as written, a reader must scan past ~40 lines of provenance bookkeeping to locate the four lines of actual computation — add a `// --- physics ---` signpost at 281 to mark where the math resumes.
- HaighMallionResult.cpp:296–323 — the find-or-create `RingNeighbourhood` loop plus polar-coordinate setup (z/rho/theta) is inline in the hot loop and reads as unrelated geometry housekeeping wedged between G-construction and storage; a one-line signpost `// locate/create per-ring record` at 296 would keep the through-line.
- HaighMallionResult.cpp:312–319 — `z`, `rho`, `theta` are computed and stored but the cylindrical-frame convention (z along normal, rho in-plane) is only inferable from the math; add `// cylindrical coords about ring axis` so the reader needn't reverse-engineer it.

## 2. Naming carries meaning

- HaighMallionResult.cpp:90–96 — `rS`, `rho`, `rho3`, `rho5` are clear from the banner; `r` (the field point, param at :81) is the one terse name — fine given the comment, but consider `r` → matches header's "field point" only after reading line 73.
- HaighMallionResult.cpp:310 — `Vec3 d` for the atom→center displacement is bare; `d_plane` at :314 is good — name it `disp` or `to_center` to match the existing `direction_to_center`.
- HaighMallionResult.cpp:401 `PackST_HM` / out — `out` ordering (T0 at [0], T1 at [1..3], T2 at [4..8]) is the SphericalTensor 9-vector layout; name is fine, the layout is obvious from the loop.

## 3. Visible math structure (grouping)

- HaighMallionResult.cpp:281–294 — the Step 1/2/3 comments are exactly the right grouping; this is the clearest part of `Compute`. No change.
- HaighMallionResult.cpp:99–103 — the `K_ab` double loop fuses weight, area, and the two kernel terms into one expression; acceptable since the banner (68–77) names it, and the inline `// K_ab = ...` at :98 restates the formula. Clean.
- HaighMallionResult.cpp:336–341 — per-type accumulation is clearly a distinct stage but rides directly after storage with no break; the existing `// Per-type T0 and T2 sums` comment carries it. Fine.

## 4. Function / method naming

- HaighMallionResult.cpp:79 `AccumulateTensor` — accumulates into `H` but the name doesn't say *which* tensor or that it's per-triangle; `AccumulateTriangleIntegral` would say what it computes. Minor.
- HaighMallionResult.cpp:158 `SurfaceIntegral` — good: says quantity, banner says it returns H (symmetric/traceless/A^-1).
- HaighMallionResult.cpp:401 `PackST_HM` — the `_HM` suffix is noise (it's a file-local static packing any SphericalTensor); not misleading, just redundant. Leave it.

## 5. Comments as signposts

- HaighMallionResult.cpp:244 — `// TODO: HM sampler needs AccumulateTensor refactor` is process/history prose that doesn't help a reader of the math; either make it actionable or drop it.
- HaighMallionResult.cpp:287–290 — the minus-sign justification is good and grounded (cites the σ = -dB/dB₀ convention and the BS parallel); keep.
- HaighMallionResult.cpp:326–330 — the trailing `// Decompose(H): T0=0, T1=0` / `// Decompose(G): T0, T1, T2` comments are excellent terse signposts on a non-obvious property. Keep.
- HaighMallionResult.h:18–19 — header says the b-index of G carries `n_b` and `G_ab = n_b·(H·n)_a` with no minus sign, but the cpp constructs `G(a,b) = -n_b·V_a` (HaighMallionResult.cpp:294, :387). Stale/incomplete comment — update the header to show the leading minus and cite the same convention note as cpp:287–290.

## 6. Correctness (secondary)

- HaighMallionResult.cpp:316 — `theta = std::atan2(d_plane.norm(), std::abs(z))` uses `d_plane.norm()` where `rho` (computed on the prior line) is the same value; harmless duplication, but `std::abs(z)` folds both sides of the ring plane onto θ∈[0,π/2] — check this is the intended convention (loses which face of the ring the atom is on; `z` retains sign, `theta` does not).
- HaighMallionResult.cpp:337 — `ti >= 0 && ti < 8` silently drops rings whose `TypeIndexAsInt()` ≥ 8 from the per-type sums while still adding them to `G_total`; check that ring type count ≤ 8 is guaranteed, else per-type and total diverge silently — consider logging the dropped case.
- HaighMallionResult.cpp:377 vs :230 — `SampleShieldingAt` guards with `distance < geom.radius` (skip inside ring footprint) but `Compute` does not apply that same radius guard, relying instead on the filter set (:213–216); confirm the filter's near-field exclusion covers the inside-the-ring case so the two paths agree on the same physics the banner at :362 claims ("Same physics as Compute()").
- HaighMallionResult.cpp:93 — singularity guard `continue`s on a single near-coincident quadrature point while still accumulating the other 6; near the surface this silently biases H rather than refining — note the adaptive subdivision (125–148) is the intended mitigation, so likely fine, but the dropped-point bias is invisible to a reader. Check threshold vs. L2 subdivision interplay.
- HaighMallionResult.cpp:171 — fan triangulation pairs `geom.center` with consecutive `verts[i], verts[j]`; vertex winding determines triangle normal orientation, but `AccumulateTensor` uses `.norm()` of the cross product (unsigned area at :85), so orientation cannot flip the integrand sign — correct and robust to winding. No issue.
