# pi-quadrupole — codex review (readability focus)

- **Targets:** src/PiQuadrupoleResult.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 9,263
- **Brief:** `../_brief.md`

---

**Verdict**

`src/PiQuadrupoleResult.h` mostly tells a coherent story. The top-level derivation is unusually helpful: it establishes source, tensor, scalar side-effect, properties, and units before the API. The main weak spot is that the public method names/comments blur “geometric kernel,” “shielding,” and “features,” so the API surface is less clear than the derivation.

`src/PiQuadrupoleResult.cpp` is scientifically traceable but not fully narrative. The kernel function has the right order, but the tensor build collapses five physical terms into one dense assignment, so a reader must keep the paper formula, coefficient signs, powers of `r`, and index roles in their head at once. The compute loop also mixes filtering, geometry-choice logging, neighbourhood creation, kernel storage, per-type accumulation, and total accumulation without enough short signposts.

**1. Coherent Story / Readability**

`src/PiQuadrupoleResult.h:64` — “Grid sampling: evaluate PQ EFG kernel” followed by `SampleShieldingAt` makes the story shift from EFG kernel to shielding without saying whether parameters are applied — rename/comment as “sample accumulated PQ EFG tensor” or similar.

`src/PiQuadrupoleResult.cpp:88-97` — the five-term EFG formula is implemented as one fused assignment — split into named intermediates matching the displayed formula, e.g. `axial_dd_term`, `cross_nd_term`, `isotropic_dd_term`, `normal_normal_term`, `trace_term`. As written, a reader must re-derive the correspondence term by term.

`src/PiQuadrupoleResult.cpp:130-134` — the filter comment is long and process-like for a block whose purpose is simple — replace with a terse signpost such as `// near-field exclusions`, with maybe one short note if the quadrupole-specific rationale is essential.

`src/PiQuadrupoleResult.cpp:181-205` — finding/creating `RingNeighbourhood` interrupts the physics flow and hides that this is metadata for later consumers — add a block label like `// ring-neighbour metadata`.

`src/PiQuadrupoleResult.cpp:196-201` — local cylindrical coordinates are computed without saying what frame they are in — add `// ring-frame coordinates`. As written, a reader must infer that `z`, `rho`, and `theta` are relative to the ring normal/plane.

`src/PiQuadrupoleResult.cpp:212-218` — per-type scalar/T2 accumulation appears before total tensor accumulation, but the narrative does not explain that this is feature bookkeeping, not part of the core tensor sum — add `// per-ring-type features`.

`src/PiQuadrupoleResult.cpp:239-257` — `SampleShieldingAt` repeats sampling/filter logic with different exclusions from `Compute`, but the story does not say why grid sampling uses only distance/radius/cutoff and not the full filter set — add a signpost or rationale if intentional.

**2. Naming Carries Meaning**

`src/PiQuadrupoleResult.cpp:47` — `G` is mathematically defensible but opaque outside the derivation — consider `efg_kernel`.

`src/PiQuadrupoleResult.cpp:48` — `scalar` is too generic for a separate physical effect — consider `buckingham_a_kernel` or `quad_field_scalar`.

`src/PiQuadrupoleResult.cpp:61` — `d` is conventional in the displayed formula, but code readability would improve with `atom_from_ring`.

`src/PiQuadrupoleResult.cpp:74` — `dn` requires remembering “d dot n” — consider `height_along_normal`.

`src/PiQuadrupoleResult.cpp:88` — `diag_term` is a matrix-location name, not a physical/mathematical one — consider `delta_trace_term` or `kronecker_delta_term`.

`src/PiQuadrupoleResult.cpp:145` — `total_pairs` actually counts accepted atom-ring kernel evaluations — consider `accepted_pairs`.

`src/PiQuadrupoleResult.cpp:148` — `ca` is too compressed for a central mutable object — consider `atom` or `conf_atom`.

`src/PiQuadrupoleResult.cpp:213` — `ti` is readable only locally — consider `ring_type_index`.

**3. Visible Math Structure**

`src/PiQuadrupoleResult.cpp:69-79` — scalar setup is clear enough; clean on structure.

`src/PiQuadrupoleResult.cpp:81-99` — tensor construction should visibly mirror the five displayed mathematical terms — split the assignment into one named value per line-group, then sum them.

`src/PiQuadrupoleResult.cpp:225-226` — “Store accumulated shielding contribution (pure T2)” hides the visible stage “decompose total EFG tensor” — comment/name should say `// decompose total EFG`.

`src/PiQuadrupoleResult.cpp:261-265` — packing `T0/T1/T2` is mechanically clear; clean.

**4. Function / Method Naming**

`src/PiQuadrupoleResult.h:65` / `src/PiQuadrupoleResult.cpp:239` — `SampleShieldingAt` sounds like final shielding, but returns a decomposed geometric EFG sum without coefficients — clearer: `SamplePiQuadTensorAt`, `SamplePiQuadKernelAt`, or `SamplePqEfgAt`.

`src/PiQuadrupoleResult.cpp:54` — `ComputePiQuadKernel` is clear enough, though `ComputePiQuadEfgKernel` would better state the returned quantity.

`src/PiQuadrupoleResult.cpp:261` — `PackST_PQ` is compressed and file-local — clearer: `PackPiQuadSphericalTensor`.

**5. Comments As Signposts**

`src/PiQuadrupoleResult.h:3-41` — the derivation comment is long, but it is grounded and useful; I would keep it.

`src/PiQuadrupoleResult.cpp:25-44` — duplicates the header derivation almost verbatim — consider replacing with `// axial quadrupole EFG kernel` plus a pointer to the header derivation.

`src/PiQuadrupoleResult.cpp:78` — comment restates the formula immediately below — acceptable, but could be shorter: `// Buckingham A kernel`.

`src/PiQuadrupoleResult.cpp:81-87` — useful as a formula signpost, but only if the implementation below is split to match it.

`src/PiQuadrupoleResult.cpp:162` — `// Apply filter: source extent = ring diameter` is good.

`src/PiQuadrupoleResult.cpp:169` — `// ---- GeometryChoice: filter exclusion ----` is visually noisy — replace with `// record exclusion`.

`src/PiQuadrupoleResult.cpp:181` — `// Find or create RingNeighbourhood` is clear but mechanical — better signpost: `// ring-neighbour metadata`.

`src/PiQuadrupoleResult.cpp:207` — good concise signpost.

`src/PiQuadrupoleResult.cpp:212` — good concise signpost, though `// per-ring-type features` carries more meaning.

`src/PiQuadrupoleResult.cpp:220` — good concise signpost.

**6. Correctness**

No definite correctness bug from the inlined text alone.

`src/PiQuadrupoleResult.cpp:130-138` vs `239-257` — grid sampling does not use the same filter set as atom computation. This may be intentional for arbitrary points, but it is worth checking because it changes which rings contribute near bonded/topological exclusions.
