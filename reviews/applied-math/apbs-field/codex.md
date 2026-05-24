# apbs-field — codex review (readability focus)

- **Targets:** src/ApbsFieldResult.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 10,614
- **Brief:** `../_brief.md`

---

**Verdict**

`src/ApbsFieldResult.h` mostly tells a coherent story: it gives the reader the physical quantity, dependency, APBS path, units, and no-fallback policy before the API. The only real drag is that it front-loads more implementation detail than a header needs, but it is understandable on one pass.

`src/ApbsFieldResult.cpp` has the right computational through-line, but it does not consistently read as a clean narrative. The main path is recognizable: prepare atom arrays → size grid → solve APBS → cache grid → finite-difference fields → symmetrize/project → store/export. The story breaks where comments become history/provenance prose, where stale section labels contradict the implementation, and where compact math names like `E`, `EFG`, `i/j/a/b/d` force the reader to re-derive roles from loops.

**1. Coherent Story / Readability**

`src/ApbsFieldResult.h:3` — clear file-level story — no change needed.

`src/ApbsFieldResult.cpp:21` — “from old ApbsSolver” is process history, not part of the physical story — replace with `// Grid interpolation`.

`src/ApbsFieldResult.cpp:68` — `ElectricFieldFromGrid` is readable, but the reader must infer units from far away — add a terse signpost such as `// central difference gradient`.

`src/ApbsFieldResult.cpp:81` — `FieldGradientFromGrid` computes derivative-of-field, then symmetrizes, then removes trace; this is the core math but the stages are buried inside one function — keep computation unchanged, but visually group as `// finite differences`, `// symmetric projection`, `// traceless projection`.

`src/ApbsFieldResult.cpp:93` — the long symmetrization comment interrupts the narrative with schema history and review provenance — reduce to a grounded comment, e.g. `// symmetric projection`; if needed, one extra sentence: “Finite differences can leave antisymmetric interpolation noise.”

`src/ApbsFieldResult.cpp:107` — this comment is scientifically useful but too long for the local flow — keep the reason, shorten to `// traceless projection` plus one sentence: “Subtract the discretized self-potential trace.”

`src/ApbsFieldResult.cpp:139` — “Separate x, y, z arrays for the C bridge” is clear — no change needed.

`src/ApbsFieldResult.cpp:160` — grid sizing is one of the clearer blocks; the physical assumptions are visible before the numbers — no change needed.

`src/ApbsFieldResult.cpp:242` — “Sanitise” hides two distinct behaviors: zero non-finite fields and clamp extreme field magnitude — split labels into `// finite-value guard` and `// field magnitude cap`.

`src/ApbsFieldResult.cpp:284` — “fall back to vacuum Coulomb on failure” contradicts the implementation and the header — replace with `// Main compute: APBS solve only`.

`src/ApbsFieldResult.cpp:327` — the WriteFeatures section comment is partly schema archaeology — replace with `// Feature export`.

**2. Naming Carries Meaning**

`src/ApbsFieldResult.cpp:69` — `E` is conventional but still terse; acceptable locally, but `electric_field` would read better in storage/sanitization code.

`src/ApbsFieldResult.cpp:82` — `EFG` is domain-recognizable, but `field_gradient` or `efg_tensor` would better signal it is a `Mat3`.

`src/ApbsFieldResult.cpp:83` / `89` — `j` and `i` require the reader to decode column-as-coordinate and row-as-field-component — rename mentally/locally to `coord_axis` and `field_axis`.

`src/ApbsFieldResult.cpp:130` — `non_authoritative_radii` is clear — no change needed.

`src/ApbsFieldResult.cpp:140` — `xArr`, `yArr`, `zArr` are C-bridge-shaped names; understandable, though `x_coords`, `y_coords`, `z_coords` would be clearer.

`src/ApbsFieldResult.cpp:164` — `lo` / `hi` are common but terse — `min_corner` / `max_corner` would better tell the geometry story.

`src/ApbsFieldResult.cpp:179` — `grid_dim` is clear enough, but since it is one axis dimension, `grid_points_per_axis` would be more precise.

`src/ApbsFieldResult.cpp:197` — `gridResult` mixes camelCase in otherwise snake-ish locals — not a readability bug, but `apbs_grid` would carry meaning better.

`src/ApbsFieldResult.cpp:335` — `N` is conventional but less readable than `atom_count`.

**3. Visible Math Structure**

`src/ApbsFieldResult.cpp:57` — trilinear interpolation is compressed into one eight-term return; mathematically standard, but a reader must verify all corner weights manually — add `// trilinear interpolation`, or name corner values/weights if this block is often audited.

`src/ApbsFieldResult.cpp:81` — finite-difference derivative, symmetrization, and traceless projection are the right sequence, but the function name only announces the first part — either rename to include projection or add block labels making the full pipeline visible.

`src/ApbsFieldResult.cpp:93-116` — the important structure is “symmetrize, then de-trace”; the current comment makes the schema consequence more prominent than the math — foreground the two projections.

`src/ApbsFieldResult.cpp:251-254` — field magnitude clamping has a magic threshold name but no local physical meaning — add `// sanity clamp` or make sure `APBS_SANITY_LIMIT` is documented where defined.

`src/ApbsFieldResult.cpp:262-266` — unit conversion is clear and well placed immediately before storage — no change needed.

**4. Function / Method Naming**

`src/ApbsFieldResult.h:46` — `ElectricFieldAt` is clear, though it does not expose units — acceptable because the header documents units.

`src/ApbsFieldResult.h:47` — `FieldGradientAt` is slightly generic; `ElectricFieldGradientAt` or `EfgTensorAt` would say more precisely what is returned.

`src/ApbsFieldResult.h:48` — `FieldGradientSphericalAt` is understandable, but `ElectricFieldGradientSphericalAt` would match the physical quantity.

`src/ApbsFieldResult.cpp:68` — `ElectricFieldFromGrid` is good.

`src/ApbsFieldResult.cpp:81` — `FieldGradientFromGrid` hides that the returned tensor has been symmetrized and traceless-projected — clearer: `ProjectedFieldGradientFromGrid` or `TracelessEfgFromGrid`.

`src/ApbsFieldResult.cpp:126` — `ComputeViaApbs` is clear.

`src/ApbsFieldResult.cpp:333` — `WriteFeatures` is framework-generic and acceptable.

**5. Comments As Signposts**

`src/ApbsFieldResult.h:42` — “Factory” is fine but slightly generic — `// compute APBS fields`.

`src/ApbsFieldResult.cpp:8` — useful bridge signpost — no change needed.

`src/ApbsFieldResult.cpp:21-23` — decorative banner plus history — replace with `// Grid interpolation`.

`src/ApbsFieldResult.cpp:74` — good terse math comment — no change needed.

`src/ApbsFieldResult.cpp:93-104` — verbose, historical, and review-process prose — replace with `// symmetric projection`.

`src/ApbsFieldResult.cpp:107-114` — too much local exposition — replace with `// traceless projection`; keep one compact reason if desired.

`src/ApbsFieldResult.cpp:196` — “Call the C bridge” is useful — no change needed.

`src/ApbsFieldResult.cpp:242` — “Sanitise” is vague — replace with `// finite-value guard`.

`src/ApbsFieldResult.cpp:268` — “both Mat3 AND SphericalTensor” over-emphasizes representation — replace with `// store tensors`.

`src/ApbsFieldResult.cpp:284-286` — stale and wrong — replace with `// Main compute`.

`src/ApbsFieldResult.cpp:347-350` — comment references stale line number and schema history — replace with `// T2 feature export`; remove the line-number reference.

**6. Correctness**

`src/ApbsFieldResult.cpp:164` — `xArr[0]` is read before any empty-conformation guard — if `AtomCount()` can be zero, this is an out-of-bounds bug; add an early `n_atoms == 0` guard or confirm the framework forbids empty conformations.

`src/ApbsFieldResult.cpp:291-293` — returns an `ApbsFieldResult` with `conf_ == nullptr` when `ChargeAssignmentResult` is missing; query methods then dereference `conf_` — if direct `Compute()` calls are possible, return `nullptr` or make the missing dependency impossible here. If the result scheduler guarantees dependencies, document that assumption.

`src/ApbsFieldResult.cpp:284` — not a runtime bug, but the stale “fall back” comment is correctness-adjacent because it tells future maintainers the opposite of the code.

No sign/frame/unit bug is obvious from the inlined code alone.
