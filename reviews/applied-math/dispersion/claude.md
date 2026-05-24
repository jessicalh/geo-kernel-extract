# dispersion — claude review (readability focus)

- **Targets:** src/DispersionResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**DispersionResult.h** — Reads cleanly. The header comment states the kernel, the tracelessness proof, the filters, and the units up front, so a physicist meets the math before the class. The interface is small and the names are honest. No trouble.

**DispersionResult.cpp** — The math core (`DispSwitchingFunction`, `ComputeDispVertex`) is genuinely clear and well-signposted. The story breaks in `Compute`, where the actual physics — find rings, exclude, sum vertices, decompose — is interleaved with three verbose `GeometryChoice::Record` provenance blocks and a hand-rolled RingNeighbourhood find-or-create, so the through-line is buried under bookkeeping. The dominant readability cost is the repeated `CalculatorConfig::Get("long_string_key")` calls inlined into expressions and log strings, which fuse a config lookup with the math and bloat several lines past readability.

---

## 1. Coherent story / readability (primary)

- DispersionResult.cpp:66-67 — `rc2`/`rs2` each built from two inline `CalculatorConfig::Get("...cutoff")` / `Get("...onset...")` calls, so a single squaring reads as a 4-call expression — as written, a reader must visually de-duplicate the repeated config key to see it is just `Rc²` and `Rs²` — hoist the two distances to named locals (`r_cut`, `r_switch`) at function top, then `rc2 = r_cut*r_cut`.
- DispersionResult.cpp:108 — guard fuses two config lookups, two comparisons, and an early-return on one line — split into two lines (singularity guard; cutoff guard) so each physical reason reads on its own.
- DispersionResult.cpp:218-346 — the per-ring loop body is ~130 lines where the three `choices.Record(...)` lambda blocks (near-field, through-bond, taper, noise-floor) outweigh the four lines of actual physics — the dispersion-summation story is hard to trace through the provenance scaffolding — consider a one-line signpost before each physics step (`// --- physics: sum kernel over vertices ---`) so the eye can skip the Record blocks.
- DispersionResult.cpp:302-327 — find-or-create RingNeighbourhood plus the z/rho/theta geometry derivation is a 25-line digression inside the main loop unrelated to dispersion — a `// locate/create ring neighbourhood record` signpost would mark it as bookkeeping, not kernel math.
- DispersionResult.cpp:337,339 — bare `8` and `5` (ring-type count, T2 component count) as loop/index bounds — magic numbers; name them (`kRingTypeCount`, `kT2Components`) or comment inline.
- DispersionResult.cpp:404,418 — `N * 40` and the `t*5 + c` indexing encode 8×5 flattening with no nearby reminder that 40 = 8 ring types × 5 T2 components — add `// 8 ring types × 5 T2`.

## 2. Naming carries meaning

- DispersionResult.cpp:110 — `double S` is the switching value; single letter, but the header/comment establish `S(r)` as the switching function so it is defensible — acceptable, no change.
- DispersionResult.cpp:272,296,332 — `s_ring` / `disp_scalar` is the isotropic 1/r⁶ sum; fine — but `K_ring`/`K_total`/`disp_total` use bare `K` which the comments tie to the kernel tensor; acceptable given the header defines `K_ab`.
- DispersionResult.cpp:392 — `PackST_D` — the `_D` suffix (Dispersion) is cryptic and the name does not say it flattens a SphericalTensor to 9 doubles T0|T1|T2 — rename `PackSphericalTensor9` or comment the layout.
- DispersionResult.cpp:315,320 — `d` and `d_plane` reused inside the neighbourhood block shadow the kernel's `d = r_atom - r_vertex` meaning from the header; here `d` is atom-to-ring-center — fine locally but a one-word note (`// d: atom → ring center`) prevents conflating the two `d`s.

## 3. Visible math structure (grouping)

- DispersionResult.cpp:121-124 — the `K_ab` double loop is the one place the tensor is actually built; clear and matches the header formula — clean.
- DispersionResult.cpp:270-298 — vertex-sum stage is well isolated and reads top-to-bottom — clean.
- DispersionResult.cpp:329-341 — store-results and per-type-accumulation are two distinct stages run back-to-back with no blank-line/label separation from the find-or-create above — a blank line + `// accumulate per ring-type` already-present comment is fine; minor.

## 4. Function / method naming

- DispersionResult.cpp:62 — `DispSwitchingFunction(double r)` returns dimensionless S(r); name says what it computes — clear.
- DispersionResult.cpp:101 — `ComputeDispVertex` returns the per-vertex K + scalar; clear.
- DispersionResult.cpp:136 — `BondedToVertices` returns the exclusion set; clear and self-documenting.
- DispersionResult.h:59 — `SampleShieldingAt` returns a kernel tensor, not shielding (the header is explicit the system outputs kernels not shielding) — name overstates; consider `SampleKernelAt` to match the project's kernel/shielding distinction. (Naming convention, not a bug.)

## 5. Comments as signposts

- DispersionResult.cpp:26-49 and 51-61 — two long block comments both explain the same R_SWITCH/R_CUT taper rationale (4.3/5.0, MD fluctuation, <0.1% truncation) back-to-back — redundant; the second block restates the first — collapse to one.
- DispersionResult.cpp:45-48 — comment cites `R_switch = 4.3 A` / `R_cut = 5.0 A` as literals, but the code reads these from `CalculatorConfig` (learnable) — risk of going stale if the TOML changes — phrase as "default 4.3/5.0 A, see calculator config".
- DispersionResult.cpp:111 — `// below switching threshold` is a good terse signpost — keep.
- DispersionResult.cpp:282 — `// Only fires when r is in the taper range but S < 1e-15` references `1e-15` literal while the code uses `dispersion_switching_noise_floor` config — minor stale-number risk; say "below noise floor".
- DispersionResult.h:34 — `Unit C6 = 1 — parameter is learnable per ring type` is a clear, grounded one-liner — good.

## 6. Correctness (secondary)

- DispersionResult.cpp:108 vs 381 — both paths guard `r > cutoff`, but `Compute` (108) also rejects `r < singularity_guard` while `SampleShieldingAt` (381) applies the same guard — consistent; no issue.
- DispersionResult.cpp:283 — the noise-floor provenance fires on `r > onset && r < cutoff`, but `ComputeDispVertex` returns invalid when `S < noise_floor`; if `r` is just below `onset`, `S≈1` and the vertex is valid, so the branch correctly excludes only true taper-range floors — looks consistent.
- DispersionResult.cpp:337 — `ti >= 0 && ti < 8` silently drops contributions from ring types with index ≥ 8 (no log) while `disp_total` (343) still includes them — check whether a ring type can have `TypeIndexAsInt() >= 8`; if so the per-type NPY and the total disagree silently.
- DispersionResult.cpp:247 — through-bond exclusion skips the *entire* ring if the atom is bonded to *any* vertex; consistent with the header's stated intent — no bug, convention is documented.
- DispersionResult.cpp:374-375 — `SampleShieldingAt` rejects `ring_dist < radius` (inside ring) and `ring_dist < singularity_guard`, but unlike `Compute` it applies no `DipolarNearFieldFilter` — check the two paths are intended to differ; the header says the filter is "same physics as all other ring calculators," so the grid path omitting it may be deliberate (raw field sampling) but is worth confirming.
