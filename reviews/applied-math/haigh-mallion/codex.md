# haigh-mallion — codex review (readability focus)

- **Targets:** src/HaighMallionResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 20,145
- **Brief:** `../_brief.md` (readability-first)
- **See also:** `codex-correctness.md` (earlier correctness-first pass)

---

Verdict: `HaighMallionResult.h` mostly tells the physical story up front, but it over-compresses the tensor convention at the most important point: the sign and outer-product order do not match the `.cpp`. A reader can follow the intent, but would have to stop and reconcile `G_ab = n_b * (H.n)_a`, `G = n (x) V`, and the implemented `G(a,b) = -n_b * V_a`.

Verdict: `HaighMallionResult.cpp` has a generally coherent top-down narrative: quadrature rule → triangle integral → adaptive refinement → ring surface integral → atom/ring accumulation → feature output. The main readability breaks are convention density around `H`, `V`, and `G`, a few cryptic scalar names in the geometry block, and duplicated construction of `G` where the sign convention is easy to miss.

**1. Coherent Story / Readability**
`src/HaighMallionResult.h:15` — header says `G_ab = n_b * (H . n)_a`, but implementation uses `G_ab = -n_b * V_a` — include the minus sign here or explicitly say the sign is applied in `Compute`.

`src/HaighMallionResult.h:25` — `"G = n (x) V"` conflicts with the indexed form and implementation; as written, a reader must infer which side is rows vs columns — write `G_ab = -V_a n_b` instead.

`src/HaighMallionResult.cpp:281` — the story is clear here, but `H`, `V`, and `G` are still dense symbolic names at the critical physical transition — rename or annotate as `surface_integral_H`, `effective_field_V`, `shielding_kernel_G`.

`src/HaighMallionResult.cpp:296` — “Find or create RingNeighbourhood” interrupts the math narrative between constructing `G` and storing/decomposing it — add a short label before line 325 like `// store ring contribution`, or move the neighbourhood bookkeeping under a clearer persistence section.

`src/HaighMallionResult.cpp:310` — the local geometry block introduces `d`, `z`, `rho`, `theta` without a signpost — as written, a reader must decode that this is ring-frame cylindrical geometry. Suggest `// ring-frame coordinates`.

`src/HaighMallionResult.cpp:360` — `SampleShieldingAt` says “Same physics as Compute” but silently applies different filters: no bonded exclusion, skips `distance < geom.radius`, and no atom-specific filtering — say `// grid-point exclusions` near lines 375–378.

**2. Naming Carries Meaning**
`src/HaighMallionResult.cpp:35` — `TriQuadPoint` is good; `lambda` and `weight` are clear for quadrature.

`src/HaighMallionResult.cpp:83` — `H` is mathematically meaningful but underspecified in code-local scope — suggest `surface_H` or `surface_integral_H`.

`src/HaighMallionResult.cpp:90` — `rS` is readable only after the comment — suggest `surface_point`.

`src/HaighMallionResult.cpp:91` — `rho` is domain-appropriate, but its frame/sign matters — suggest `field_minus_surface`.

`src/HaighMallionResult.cpp:285` — `V` is too small for an important physical quantity — suggest `effective_field`.

`src/HaighMallionResult.cpp:291` — `G` is tolerable after the nearby comment, but `shielding_kernel` would carry the returned quantity better.

`src/HaighMallionResult.cpp:310` — `d` hides direction and origin — suggest `center_to_atom` or `atom_from_center`.

`src/HaighMallionResult.cpp:315` — `rho` reuses the earlier kernel displacement name for cylindrical radius — suggest `radial_distance`.

`src/HaighMallionResult.cpp:401` — `PackST_HM` is terse and implementation-shaped — suggest `PackHaighMallionTensor`.

**3. Visible Math Structure**
`src/HaighMallionResult.cpp:98` — kernel evaluation fuses weight, area, dipolar tensor, and accumulation in one expression — split into `kernel_ab` and `quadrature_weight` if touched.

`src/HaighMallionResult.cpp:131` — adaptive threshold selection is clear enough, but the two thresholds would read better as named locals before the branch — no computation change needed.

`src/HaighMallionResult.cpp:281` — the key physical sequence is well grouped: `H = SurfaceIntegral`, `V = H * normal`, `G = outer product`. This is the strongest narrative part of the file.

`src/HaighMallionResult.cpp:325` — storage and spherical decomposition are visibly grouped and readable.

`src/HaighMallionResult.cpp:380` — `SampleShieldingAt` duplicates the `H → V → G` math; readability would benefit from the same exact signpost wording as `Compute`.

**4. Function / Method Naming**
`src/HaighMallionResult.h:50` — `Compute` is conventional for result classes and acceptable.

`src/HaighMallionResult.h:54` — `SampleShieldingAt` is clear, though it returns a decomposed `SphericalTensor`; `SampleShieldingTensorAt` would be more explicit.

`src/HaighMallionResult.cpp:44` — `Gauss7` is understandable but terse — suggest `TriangleGauss7Rule`.

`src/HaighMallionResult.cpp:79` — `AccumulateTensor` does not say which tensor or domain — suggest `AccumulateTriangleSurfaceIntegral`.

`src/HaighMallionResult.cpp:125` — `AccumulateAdaptive` is too generic — suggest `AccumulateAdaptiveTriangleIntegral`.

`src/HaighMallionResult.cpp:158` — `SurfaceIntegral` is clear in context, but `ComputeSurfaceIntegralH` would make the returned quantity explicit.

**5. Comments As Signposts**
`src/HaighMallionResult.h:3` — header block is useful but long; it reads more like design notes than a quick contract — keep the equation/units/convention, move empirical/process prose out if this is meant as API documentation.

`src/HaighMallionResult.h:23` — “This parallels BiotSavartResult…” is useful context but not necessary for reading this class — consider shortening to `// Biot-Savart analogue`.

`src/HaighMallionResult.cpp:27` — banner comment is readable, but large banners dominate the file — terse replacement: `// triangle quadrature`.

`src/HaighMallionResult.cpp:68` — good signpost, but verbose — terse replacement: `// dipolar surface integral`.

`src/HaighMallionResult.cpp:108` — good signpost — terse replacement: `// adaptive subdivision`.

`src/HaighMallionResult.cpp:238` — `// ---- GeometryChoice: surface integral ----` is process-oriented and partly opaque — terse replacement: `// record source ring`.

`src/HaighMallionResult.cpp:250` — `// Apply filter` restates code — terse replacement: `// pair filters`.

`src/HaighMallionResult.cpp:269` — good signpost — terse replacement: `// adaptive refinement`.

`src/HaighMallionResult.cpp:287` — this is an important convention comment and should stay, but align it with the header exactly.

`src/HaighMallionResult.cpp:396` — feature layout comment is clear and useful.

**6. Correctness**
`src/HaighMallionResult.h:15` / `src/HaighMallionResult.cpp:287` — possible documentation/convention mismatch, not enough to call a code bug: header omits the minus sign while implementation applies it. Check against `BiotSavartResult` and Boyd & Skrynnikov convention.

`src/HaighMallionResult.h:25` / `src/HaighMallionResult.cpp:294` — possible outer-product wording mismatch: `G = n (x) V` suggests `n_a V_b`, while code implements `-V_a n_b`. Check intended tensor index order.

`src/HaighMallionResult.cpp:316` — `theta = atan2(d_plane.norm(), abs(z))` discards side-of-plane sign. This may be intentional cylindrical geometry; I would only check whether downstream consumers expect signed polar angle.

No definite arithmetic, bounds, or accumulation bug is visible from the inlined text alone.
