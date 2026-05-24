# biot-savart — codex review (readability focus)

- **Targets:** src/BiotSavartResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 20,558
- **Brief:** `../_brief.md` (readability-first)

---

**Verdict**

`src/BiotSavartResult.h` mostly tells the physics story, but it over-compresses the sign/kernel distinction and leaves the public API ambiguous about whether it returns a geometric kernel or an intensity-scaled shielding tensor.

`src/BiotSavartResult.cpp` has a coherent high-level order: segment field → Johnson-Bovey loop → per-atom accumulation → sampling → export. The main readability break is inside `Compute()`, where physics, filtering, geometry-choice logging, neighborhood feature construction, tensor accumulation, and export-side counters are interleaved; a reader can follow it, but only by carrying several parallel responsibilities at once.

**1. Coherent Story / Readability**

`src/BiotSavartResult.h:10` — says `G_ab = n_b * B_a * PPM_FACTOR`, then later line 23 says `G_ab = -n_b * B_a * PPM_FACTOR`; the story contradicts itself before the sign convention arrives — make line 10 include the minus sign or explicitly say “before shielding sign convention.”

`src/BiotSavartResult.h:54` — `SampleShieldingAt` says “Returns the shielding tensor,” but the implementation appears to return the unit-current geometric kernel sum, not intensity-scaled shielding — say “returns spherical tensor decomposition of the summed unit-current G kernel” if that is intended.

`src/BiotSavartResult.cpp:145` — the filter setup comment is long and implementation-specific before the main physics loop starts; as written, a reader must switch from physics to policy plumbing before reaching the computation — replace with a short signpost like `// atom-ring filters`.

`src/BiotSavartResult.cpp:176` — the `GeometryChoice` recording block interrupts the ring-current computation before distance/filter/B/G are visible; as written, a reader must skip 20 lines of logging/sampling machinery to recover the physical path — add a terse signpost such as `// provenance sample` or move the reader-facing comment to say this block records one sampler per ring.

`src/BiotSavartResult.cpp:219` — this is where the physical through-line resumes, but it is buried after logging/filter code — a short block label like `// unit-current field` would help.

`src/BiotSavartResult.cpp:237` — “Find or create RingNeighbourhood” begins a large feature-population block in the middle of tensor computation; as written, a reader must track that the physics result is already computed and this is storage/feature geometry — label it `// ring-neighbour features`.

`src/BiotSavartResult.cpp:254` — cylindrical coordinate construction is readable but dense enough that its role is not immediately clear — use `// target in ring frame`.

`src/BiotSavartResult.cpp:280` — storing `G`, spherical decomposition, and `B` is clear.

`src/BiotSavartResult.cpp:297` — accumulation is clear.

`src/BiotSavartResult.cpp:319` — proximity counts are easy to read, but they count over all `ca.ring_neighbours`, not just rings accepted in this invocation if the atom already had neighbours — if accumulation across earlier calculators is intended, say so; otherwise this is worth checking.

**2. Naming Carries Meaning**

`src/BiotSavartResult.cpp:47` — `dl_m` is acceptable because the preceding formula defines it.

`src/BiotSavartResult.cpp:48` — `dA_m` and `dB_m` are defined in the comment, but still easy to confuse with magnetic field `B`; suggest `start_to_point_m` and `end_to_point_m` if touching this code.

`src/BiotSavartResult.cpp:59` — `factor` hides physical content — suggest `biot_savart_scale`.

`src/BiotSavartResult.cpp:60` — `dotTerm` names syntax, not quantity — suggest `endpoint_projection_term`.

`src/BiotSavartResult.cpp:92` — `B` is fine in the local physics block.

`src/BiotSavartResult.cpp:167` — `G_total` is clear enough.

`src/BiotSavartResult.cpp:237` — `rn` is compact but local; acceptable after the comment, though `ring_neighbour` would reduce decoding cost.

`src/BiotSavartResult.cpp:251` — `d` is too generic in a geometry-heavy block — suggest `center_to_atom`.

`src/BiotSavartResult.cpp:256` — `d_plane` is understandable but could be clearer — suggest `in_plane_offset`.

`src/BiotSavartResult.cpp:433` — `PackST` is terse and codebase-ish; suggest `PackSphericalTensor`.

**3. Visible Math Structure**

`src/BiotSavartResult.cpp:87` — unit conversion is well grouped.

`src/BiotSavartResult.cpp:96` — upper/lower loop construction is clear.

`src/BiotSavartResult.cpp:225` — tensor construction is visibly separated and explained well.

`src/BiotSavartResult.cpp:251` — ring-frame coordinate derivation mixes `z`, radial projection, polar angle, and azimuth in one block — keep computation unchanged but split with labels: `// axial/radial coordinates` and `// azimuth reference`.

`src/BiotSavartResult.cpp:285` — B-field cylindrical projection repeats the ring-frame offset calculation from lines 251–258; as written, a reader must verify the repeated math is intentionally the same — add a signpost like `// project B into ring frame`.

`src/BiotSavartResult.cpp:312` — `total_G_spherical` and `bs_shielding_contribution` are both decompositions of related totals; as written, a reader must infer why both exist — add a terse distinction comment or confirm whether both are intentionally retained.

`src/BiotSavartResult.cpp:356` and `src/BiotSavartResult.cpp:386` — sampling methods repeat the same ring acceptance logic; not asking for abstraction, but the narrative would benefit from matching comments so readers know the duplicated guards are deliberate.

**4. Function / Method Naming**

`src/BiotSavartResult.h:57` — `SampleShieldingAt` is likely misleading if it returns the unit-current kernel tensor rather than an intensity-scaled shielding tensor — suggest `SampleShieldingKernelAt` or document the unit-current return explicitly.

`src/BiotSavartResult.cpp:43` — `WireSegmentField` is clear; units are carried by parameter suffixes.

`src/BiotSavartResult.cpp:77` — `JohnsonBoveyField` is clear, though `JohnsonBoveyBField` would make the returned quantity explicit.

`src/BiotSavartResult.cpp:433` — `PackST` is too compressed for a file aimed at scientific readability — suggest `PackSphericalTensor`.

**5. Comments as Signposts**

`src/BiotSavartResult.h:3` — useful file-level signpost.

`src/BiotSavartResult.h:18` — citation/comment is useful, but dense; acceptable because it anchors convention.

`src/BiotSavartResult.cpp:26` — section header is clear, though visually heavy.

`src/BiotSavartResult.cpp:87` — good signpost: keep.

`src/BiotSavartResult.cpp:96` — “Upper loop” / “Lower loop” are good.

`src/BiotSavartResult.cpp:176` — `// ---- GeometryChoice: ring current ----` is process-oriented; replace with `// provenance sampler`.

`src/BiotSavartResult.cpp:200` — good signpost.

`src/BiotSavartResult.cpp:207` — `// ---- GeometryChoice: near-field exclusion ----` is process-oriented; replace with `// record exclusion`.

`src/BiotSavartResult.cpp:225` — useful, but too verbose through line 231 given the header already explains the sign; replace with `// shielding sign convention` plus one formula line if needed.

`src/BiotSavartResult.cpp:237` — replace with `// ring-neighbour features`.

`src/BiotSavartResult.cpp:254` — replace with `// target in ring frame`.

`src/BiotSavartResult.cpp:263` — good signpost, but lines 264–265 are explanatory prose; could be shortened to `// asymmetric ring side`.

`src/BiotSavartResult.cpp:327` — replace `// ---- GeometryChoice: ring shells ----` with `// record ring shells`.

`src/BiotSavartResult.cpp:348` — comment says “No filter set” then “DipolarNearFieldFilter still applied”; that reads contradictory even if meaning is “manual near-field guard” — rewrite as “No atom/topology filters; keep distance validity guards.”

**6. Correctness**

`src/BiotSavartResult.h:10` / `src/BiotSavartResult.h:23` — real documentation inconsistency: the sign of `G_ab` differs. The implementation uses the negative sign.

`src/BiotSavartResult.h:54` / `src/BiotSavartResult.cpp:386` — possible API contract bug: `SampleShieldingAt` claims shielding but appears to sum unit-current `G` kernels without multiplying by `ring.Intensity()`. If this is intentional, the name/comment should say kernel; if not, this is a physics-scale bug.

`src/BiotSavartResult.cpp:321` — possible accumulation hazard: ring proximity counters are incremented from `ca.ring_neighbours` without visible reset in this file. If `Compute()` can run more than once on the same conformation, counts may double. Check initialization/reset ownership before treating this as a bug.

`src/BiotSavartResult.cpp:258` — `theta = atan2(d_plane.norm(), abs(z))` discards the sign of `z`. That may be intentional for a polar angle folded across the ring plane, but the name `theta` alone does not reveal the convention. Check downstream expectations before changing.
