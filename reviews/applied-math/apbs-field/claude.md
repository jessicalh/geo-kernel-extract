# apbs-field — claude review (readability focus)

- **Targets:** src/ApbsFieldResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**ApbsFieldResult.h** — Clean. The header comment tells the whole story up front: dependency, what is stored, that the grid comes from a PB solve in kT/e, how E and EFG are extracted, the no-fallback policy, and the unit conversion with its constant. A domain reader follows it top to bottom without effort.

**ApbsFieldResult.cpp** — Mostly coherent and reads as a clear three-act story (grid utilities → solve+extract → query/write). The trilinear interpolation and central-difference stages are well-signposted with grounded comments. Two friction points: the header comment block still advertises a "vacuum Coulomb fallback" that the code explicitly does not do, and the per-atom sanitise/convert/store loop fuses several distinct steps with no internal signposts. Naming is generally good; the worst offenders are the C-bridge `Arr` arrays and a couple of unexplained magic numbers.

---

## 1. Coherent story / readability

- ApbsFieldResult.cpp:284–286 — Section banner says "try APBS, fall back to vacuum Coulomb on failure" but `Compute` returns `nullptr` with an explicit "No fallback" message; the banner contradicts the code and the header policy — change banner to "try APBS; nullptr on failure (no fallback)".
- ApbsFieldResult.cpp:242–273 — One loop body fuses four conceptually distinct steps (NaN/inf sanitise, magnitude clamp, unit conversion, store + decompose) with only a bare `// Sanitise` label — add 2–4 word signposts: `// clamp runaway field`, `// kT/e → V/A`, `// store Mat3 + T2`.
- ApbsFieldResult.cpp:160–179 — Grid-sizing comment lists four facts but the code computes only three of them, and `grid_dim = 161` sits orphaned below the block it belongs to — move `grid_dim` declaration up under the comment so the "161 per axis" line and the variable are adjacent.
- ApbsFieldResult.cpp:243–260 — E-field gets clamp+zero-on-NaN but EFG gets only zero-on-NaN (no magnitude clamp); the asymmetry is intentional-looking but unexplained — one-line note why EFG needs no sanity limit (derives from already-clamped E? check).

## 2. Naming carries meaning

- ApbsFieldResult.cpp:140,164–167 — `xArr/yArr/zArr` read as type-noise, not quantity; these are atom coordinate components in Å — rename `x_coords`/`y_coords`/`z_coords` (units in the nearby comment).
- ApbsFieldResult.cpp:164–171 — `lo`/`hi`/`extent` are fine, but `lo` is seeded from atom 0 then re-min'd including atom 0; harmless, but `lo`/`hi` would read clearer as `bbox_min`/`bbox_max`.
- ApbsFieldResult.cpp:81–90 — In `FieldGradientFromGrid`, index `j` is the differentiation direction and `i` the field component; the `EFG(i,j)` convention is correct but the loop ordering (outer `j`, inner `i`) inverts reader expectation — a one-line `// EFG(i,j) = dE_i/dr_j` at line 82 pins it.
- ApbsFieldResult.cpp:35 — `frac` is good; the per-axis `fx/fy/fz` (the within-cell offsets) and `ix/iy/iz` (cell indices) are standard and clear — no change.

## 3. Visible math structure (grouping)

- ApbsFieldResult.cpp:57–64 — The 8-term trilinear sum is dense but correctly laid out one corner per line with aligned weights; the shape is legible — no change.
- ApbsFieldResult.cpp:81–119 — `FieldGradientFromGrid` cleanly sequences build → symmetrize → de-trace, each with a comment header; this is the model the rest of the file should follow — no change.
- ApbsFieldResult.cpp:236–273 — The extraction loop is the one block where stages are not visibly grouped (see 1.b); the math (`E *= KT_OVER_E_298K; EFG *= KT_OVER_E_298K`) is correct but buried among sanitise/store noise.

## 4. Function / method naming

- ApbsFieldResult.cpp:68,81 — `ElectricFieldFromGrid` / `FieldGradientFromGrid` say exactly what they return — clean.
- ApbsFieldResult.cpp:126 — `ComputeViaApbs` returning `bool` (success) while mutating `conf` in place is clear from context — clean.
- ApbsFieldResult.h:46–48 — `ElectricFieldAt`/`FieldGradientAt`/`FieldGradientSphericalAt` are unambiguous — clean.

## 5. Comments as signposts

- ApbsFieldResult.cpp:93–104 — The symmetrization comment is excellent science but verbose (12 lines) and carries process history ("Per R4 codex review 2026-05-18") — keep the physics, cut the review-provenance line; the "why symmetrize" is the load-bearing part.
- ApbsFieldResult.cpp:39–40 — `floor() not static_cast<int>()` comment is exactly the right kind of grounded gotcha note — keep.
- ApbsFieldResult.cpp:107–114 — Traceless-projection comment is clear and correctly motivated (self-potential delta-function artifact) — keep; could tighten but not required.
- ApbsFieldResult.cpp:242 — `// Sanitise` is too terse for the four-step block it heads (see 1.b).
- ApbsFieldResult.cpp:262–264 — Unit-conversion comment duplicates the header verbatim; one is enough — trim inline to `// kT/(e·A) → V/A; kT/(e·A^2) → V/A^2`.
- ApbsFieldResult.cpp:347–350 — WriteFeatures T2-only comment is accurate and well-placed; the "line 102" cross-reference is now stale (symmetrization is at line 105) — drop the line number, say "symmetrized + traceless above".

## 6. Correctness (secondary)

- ApbsFieldResult.cpp:179,182–185 — Magic numbers `grid_dim = 161`, `pdie 4.0`, `sdie 78.54`, `ionic_strength 0.15` are commented but `40.0`/`30.0` padding at lines 175–176 and the `2026-05-18` schema constants are not parameterised; values look standard (check 78.54 vs 78.4 for water at 298 K, and that 0.15 M matches the ionic-strength units the bridge expects).
- ApbsFieldResult.cpp:45–47 — Interpolate returns 0.0 for out-of-bounds points; the central-difference stencils at lines 71–73 and 84–86 step ±spacing, so an atom within one cell of the grire edge silently yields a one-sided (biased) derivative against a zero — likely harmless given 20 Å padding, but worth a guard-comment that padding guarantees interior stencils.
- ApbsFieldResult.cpp:251–254 — Field clamp rescales E to `APBS_SANITY_LIMIT` but leaves the corresponding EFG unclamped, so a clamped atom keeps an unscaled (possibly huge) gradient — confirm this is intended (clamp is a guard, not a physical rescale of the pair).
- ApbsFieldResult.cpp:288–311 — `Compute` returns a default-constructed result (not nullptr) when `ChargeAssignmentResult` is absent (line 292) but nullptr when APBS fails — two different "no data" signals; check downstream consumers distinguish them correctly.
- ApbsFieldResult.h:43 vs .cpp:43 — `Interpolate` divides by `spacing(d)` with no zero-guard; a degenerate axis (spacing 0) would NaN, but the bridge sets spacing, so low risk — no action unless bridge can emit zero spacing.
