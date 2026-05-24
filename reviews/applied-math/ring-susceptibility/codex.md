# ring-susceptibility — codex review (readability focus)

- **Targets:** src/RingSusceptibilityResult.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 9,033
- **Brief:** `../_brief.md`

---

**Verdict**

`src/RingSusceptibilityResult.h` mostly tells a coherent story: it names the physical model, gives the tensor form, states units, and warns that the tensor is asymmetric and non-traceless. The only weak spot is that the public API names do not fully say whether they return the full asymmetric susceptibility tensor or only a decomposed shielding feature.

`src/RingSusceptibilityResult.cpp` is scientifically legible in `ComputeRingChiKernel`, but the through-line breaks in `Compute`: the physical calculation is buried inside bookkeeping, filtering, logging, neighbourhood mutation, and feature storage. A reader can reconstruct it, but as written they must hold kernel quantities, filter semantics, ring-neighbour state, and tensor accumulation in their head at once.

**1. Coherent Story / Readability**

`src/RingSusceptibilityResult.cpp:45` — `RingChiKernelResult` appears after a long derivation comment, but its fields do not map cleanly back to the named terms in the comment — suggestion: name the fields as the outputs in the derivation, e.g. `full_tensor_over_r3`, `dipolar_kernel_over_r3`, `scalar_kernel`.

`src/RingSusceptibilityResult.cpp:61-75` — setup is readable, but the story jumps from displacement to scalar without labeling the geometric quantities — suggestion: add terse block labels such as `// displacement geometry` and `// scalar kernel`.

`src/RingSusceptibilityResult.cpp:77-91` — the two tensor builds are clear only if the reader keeps the formula above in view — suggestion: split with labels `// dipolar kernel` and `// full susceptibility tensor`; consider named identity term or `delta_ab`.

`src/RingSusceptibilityResult.cpp:124-128` — this is process prose, not a signpost, and it interrupts the main computation before the atom loop — suggestion: reduce to `// exclusion filters` plus only the boundary detail if it is genuinely needed.

`src/RingSusceptibilityResult.cpp:141-219` — the main loop combines six concerns: atom iteration, spatial lookup, filter setup, geometry-choice audit logging, neighbourhood creation, kernel storage, and accumulation — suggestion: add visible block labels in sequence: `// nearby rings`, `// filter pair`, `// ring neighbourhood`, `// store kernel`, `// accumulate atom total`.

`src/RingSusceptibilityResult.cpp:176-205` — neighbourhood lookup/creation dominates the middle of the physical loop — as written, a reader must pause the tensor story to reverse-engineer data ownership — suggestion: add a short signpost comment saying this is attaching per-ring features to the atom, not part of the kernel math.

`src/RingSusceptibilityResult.cpp:193-201` — “Cylindrical coordinates in ring frame” is underspecified because `theta` is not the usual azimuthal cylindrical angle — as written, a reader must infer that this is a folded polar angle from the ring normal — suggestion: rename/comment as `// folded polar coordinates` or name `theta` more specifically.

`src/RingSusceptibilityResult.cpp:231-249` — `SampleShieldingAt` repeats filtering and kernel accumulation, but the relationship to `Compute` is not stated — suggestion: add `// grid-point accumulation` and label the three cutoffs.

`src/RingSusceptibilityResult.cpp:253-256` — `PackST_RS` is opaque and appears abruptly — suggestion: rename or comment `// feature layout`.

`src/RingSusceptibilityResult.h:3-20` — clean overall; the header gives a useful physical narrative.

**2. Naming Carries Meaning**

`src/RingSusceptibilityResult.cpp:46` — `M_over_r3` preserves formula notation but not meaning — suggestion: `full_chi_tensor_over_r3` or `shielding_kernel_tensor`.

`src/RingSusceptibilityResult.cpp:47` — `K` is too compressed outside the derivation — suggestion: `dipolar_kernel`.

`src/RingSusceptibilityResult.cpp:48` — `f` is too compressed — suggestion: `scalar_kernel`.

`src/RingSusceptibilityResult.cpp:61` — `d` is acceptable in math, but costs readability in production code — suggestion: `ring_to_atom`.

`src/RingSusceptibilityResult.cpp:69` — `d_hat` is recognizable but still local-symbol-heavy — suggestion: `ring_to_atom_hat`.

`src/RingSusceptibilityResult.cpp:141-143` — `ai`, `ca` are compact but force decoding — suggestion: `atom_index`, `atom`.

`src/RingSusceptibilityResult.cpp:150-152` — `ri`, `geom` are tolerable in tight loops, but this loop is long — suggestion: `ring_index`, `ring_geom`.

`src/RingSusceptibilityResult.cpp:179` — `rn` is cryptic across a 27-line block — suggestion: `ring_neighbour`.

`src/RingSusceptibilityResult.cpp:191` — `direction_to_center = kernel.direction` appears semantically reversed, because `kernel.direction` is “from ring center to atom” — suggestion: rename one side or add a contract comment. A domain reader will otherwise question the sign.

`src/RingSusceptibilityResult.cpp:198` — `theta` sounds like azimuth in cylindrical coordinates, but the expression is inclination from the normal folded by `abs(z)` — suggestion: `polar_angle_to_ring_normal` or `folded_theta`.

`src/RingSusceptibilityResult.cpp:253` — `PackST_RS` is internal shorthand — suggestion: `PackRingChiSphericalTensor`.

**3. Visible Math Structure**

`src/RingSusceptibilityResult.cpp:74-91` — the scalar, dipolar tensor, and full tensor are computed in sequence, which is good; the grouping would be stronger with names matching the physical outputs — suggestion: keep the order, improve labels and field names.

`src/RingSusceptibilityResult.cpp:87-91` — the full tensor expression fuses asymmetric term, ring-normal term, dipolar term, identity term, and `1/r3` scaling — suggestion: introduce named local terms inside the loop, e.g. `asymmetric_term`, `normal_term`, `dipolar_term`, without changing arithmetic.

`src/RingSusceptibilityResult.cpp:154-163` — kernel computation happens before filtering because distance is needed, but that dependency is not visible — suggestion: comment `// compute distance for filters`.

`src/RingSusceptibilityResult.cpp:207-218` — decomposition and accumulation are adjacent and clear enough; `chi_spherical` and `ringchi_shielding_contribution` carry useful meaning.

`src/RingSusceptibilityResult.h:12-20` — math structure is unusually visible for a header; no major issue.

**4. Function / Method Naming**

`src/RingSusceptibilityResult.cpp:54` — `ComputeRingChiKernel` is mostly clear, but “Chi” may be too insider-shorthand for a reader scanning names — suggestion: `ComputeRingSusceptibilityKernel`.

`src/RingSusceptibilityResult.h:43` — `SampleShieldingAt` says where but not that it returns spherical decomposition of the accumulated full tensor — suggestion: `SampleShieldingTensorAt` or `SampleRingChiShieldingAt`.

`src/RingSusceptibilityResult.cpp:253` — `PackST_RS` does not say what layout it writes — suggestion: `PackSphericalTensor9`.

`src/RingSusceptibilityResult.h:39-41` — `Compute` is conventional for this result type; clean.

`src/RingSusceptibilityResult.h:46` — `WriteFeatures` is conventional; clean.

**5. Comments As Signposts**

`src/RingSusceptibilityResult.cpp:25-43` — the derivation comment is useful, but visually heavy and partly duplicates the header — suggestion: keep the formula and outputs; remove banner weight or reduce prose.

`src/RingSusceptibilityResult.cpp:74` — good signpost, but verbose — suggestion: `// scalar kernel`.

`src/RingSusceptibilityResult.cpp:77` — good signpost — suggestion: `// dipolar kernel`.

`src/RingSusceptibilityResult.cpp:83-84` — useful formula, but it restates the code directly after the derivation block — suggestion: `// full susceptibility tensor`.

`src/RingSusceptibilityResult.cpp:124-128` — too explanatory for an in-loop setup comment — suggestion: `// exclusion filters`.

`src/RingSusceptibilityResult.cpp:145` — good signpost — suggestion: `// nearby rings`.

`src/RingSusceptibilityResult.cpp:157` — good, but slightly wordy — suggestion: `// filter pair`.

`src/RingSusceptibilityResult.cpp:164` — banner-style and process-oriented — suggestion: `// record exclusion`.

`src/RingSusceptibilityResult.cpp:176-178` — useful context, but too long — suggestion: `// ring-neighbour record`.

`src/RingSusceptibilityResult.cpp:193` — misleading if `theta` is not azimuth — suggestion: `// ring-frame coordinates` or `// folded polar coordinates`.

`src/RingSusceptibilityResult.cpp:207` — good signpost — suggestion: `// store ring kernel`.

`src/RingSusceptibilityResult.cpp:212` — good signpost — suggestion: `// accumulate atom tensor`.

`src/RingSusceptibilityResult.cpp:217` — good signpost — suggestion: `// decompose atom total`.

**6. Correctness**

No confirmed correctness bug from the inlined text alone.

`src/RingSusceptibilityResult.cpp:191` — convention risk: `direction_to_center` receives a vector documented as ring center to atom — suggestion: check the `RingNeighbourhood` contract before trusting the sign.

`src/RingSusceptibilityResult.cpp:198` — convention risk: `theta = atan2(rho, abs(z))` is a folded polar inclination, not cylindrical azimuth — suggestion: check downstream expectation for `theta`; if intended, rename/comment it so it does not look like an angle bug.
