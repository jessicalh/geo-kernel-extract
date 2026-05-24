# dispersion — codex review (readability focus)

- **Targets:** src/DispersionResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 11,099
- **Brief:** `../_brief.md` (readability-first)

---

**Verdict**

`src/DispersionResult.h`: Mostly coherent. The header gives a clear physical story before the API, but it is doing too much: it repeats implementation details, filter policy, and validation rationale that belong closer to the computation.

`src/DispersionResult.cpp`: The computation is readable in broad outline: switch function → per-vertex kernel → bonded exclusions → atom/ring accumulation → feature writing. The main readability cost is compression inside `Compute`: physics, geometry-choice logging, filtering, accumulation, and storage are interleaved, so the mathematical through-line gets interrupted even though the algorithm itself is understandable.

**1. Coherent Story / Readability**

`src/DispersionResult.h:3` — header comment is a full algorithm note, not just an interface signpost. A reader must separate public contract from implementation policy. Suggest trimming header to purpose, output units, and stored quantities; leave filters/switching details in `.cpp`.

`src/DispersionResult.cpp:26` — switching comment is coherent but overlong and partly repeated by lines 51–61. Suggest one compact formula block plus terse labels: `// switching polynomial`, `// cutoff rationale`.

`src/DispersionResult.cpp:184` — ring-level filter comment is useful, but it interrupts the compute narrative before any accumulation starts. Suggest shorter signpost: `// ring near-field filter`.

`src/DispersionResult.cpp:226-241` — filtering and `GeometryChoice` recording are fused. As written, a reader must mentally skip logging scaffolding to recover the physical step. Suggest a signpost comment before the branch and a shorter logging label, e.g. `// reject near-field atom`.

`src/DispersionResult.cpp:244-258` — same issue for through-bond exclusion: the science is clear, but the logging block dominates the branch. Suggest `// reject bonded atom` and keep the detailed reason inside the recorded metadata.

`src/DispersionResult.cpp:270-298` — the actual core computation is easy to miss because it is shorter than the surrounding bookkeeping. Suggest a stronger signpost: `// vertex kernel sum`.

`src/DispersionResult.cpp:302-327` — “Find or create RingNeighbourhood” is procedural, but this block also computes ring-relative geometry. As written, a reader must infer why `z`, `rho`, and `theta` are being set here. Suggest splitting with comments: `// neighbour record` and `// ring-frame coordinates`.

`src/DispersionResult.cpp:329-343` — storage, spherical decomposition, per-type accumulation, and tensor total are packed together. Suggest labels: `// store ring result`, `// per-type sums`, `// atom total`.

`src/DispersionResult.cpp:363-388` — `SampleShieldingAt` is readable but less coherent than `Compute`: it applies a different-looking ring filter without explaining whether this is intentionally equivalent. Suggest a signpost: `// grid-sampling guard`.

**2. Naming Carries Meaning**

`src/DispersionResult.cpp:62` — `DispSwitchingFunction` is understandable but abbreviated. Suggest `DispersionSwitchingFunction`.

`src/DispersionResult.cpp:66-70` — `rc2`, `rs2`, `num`, `den` are compact math names. A reader must map them back to cutoff/onset. Suggest `cutoff_sq`, `switch_sq`, `numerator`, `denominator`.

`src/DispersionResult.cpp:94-98` — `DispVertexResult::K` is mathematically recognizable, but not self-describing. Suggest `kernel_tensor`.

`src/DispersionResult.cpp:110` — `S` is formula-faithful but easy to lose in code. Suggest `switch_weight`.

`src/DispersionResult.cpp:113` — `d` needs its direction/frame in the name. Suggest `atom_minus_vertex`.

`src/DispersionResult.cpp:271-273` — `K_ring`, `s_ring` are partly clear, partly compressed. Suggest `ring_kernel_tensor`, `ring_scalar_sum`.

`src/DispersionResult.cpp:279` — `vr` is opaque. Suggest `vertex_result`.

`src/DispersionResult.cpp:303` — `rn` is local but cryptic in a dense block. Suggest `ring_neighbour`.

`src/DispersionResult.cpp:315-323` — `d`, `z`, `d_plane` need frame meaning. Suggest `center_to_atom`, `normal_offset`, `in_plane_offset`.

`src/DispersionResult.cpp:392` — `PackST_D` is very compressed. Suggest `PackSphericalTensorFloat64`.

**3. Visible Math Structure**

`src/DispersionResult.cpp:118-124` — scalar and tensor kernel are computed directly, but the inverse powers are not grouped by meaning. Suggest named intermediates `inv_r6`, `inv_r8`, or `radial_scalar` if consistent with surrounding code.

`src/DispersionResult.cpp:120-124` — tensor expression fuses dyadic term, trace term, and switch weight. Suggest named terms: `dyadic_term`, `isotropic_term`, then assign `switch_weight * (...)`.

`src/DispersionResult.cpp:319-323` — ring-frame coordinate derivation is visible but unlabeled. Suggest `// ring-frame coordinates`.

`src/DispersionResult.cpp:331` — decomposition appears at storage time without a math-stage label. Suggest `// spherical decomposition`.

**4. Function / Method Naming**

`src/DispersionResult.h:55` — `Compute` is conventional for result classes and acceptable.

`src/DispersionResult.h:59` — `SampleShieldingAt` returns only the dispersion spherical tensor, not total shielding. Suggest `SampleDispersionAt`.

`src/DispersionResult.cpp:101` — `ComputeDispVertex` is understandable but abbreviated. Suggest `ComputeDispersionVertexKernel`.

`src/DispersionResult.cpp:136` — `BondedToVertices` is clear enough, though `AtomsBondedToRingVertices` would be more explicit.

`src/DispersionResult.cpp:392` — `PackST_D` does not tell a human what format is being packed. Suggest `PackSphericalTensorFloat64`.

**5. Comments As Signposts**

`src/DispersionResult.cpp:26-49` — too verbose for a local signpost. Replace the local label with `// switching function`; keep only formula/reference if needed.

`src/DispersionResult.cpp:75-92` — useful, but repeats the header. Replace block label with `// per-vertex kernel`.

`src/DispersionResult.cpp:131-134` — clear; could be shorter: `// bonded vertex set`.

`src/DispersionResult.cpp:152-161` — useful overview, but repeats nearby comments. Keep as `// atom-ring accumulation`.

`src/DispersionResult.cpp:260` — “dispersion taper” is really recording metadata, not applying the taper. Suggest `// record taper parameters`.

`src/DispersionResult.cpp:281-282` — comment is helpful but too process-specific. Suggest `// switching floor`.

`src/DispersionResult.cpp:302` — procedural comment. Suggest `// neighbour record`.

`src/DispersionResult.cpp:335` — good terse signpost: `// Per-type accumulation`.

**6. Correctness**

No definite correctness bug from the inlined text alone.

One thing to check, not a claimed bug: `src/DispersionResult.cpp:372-376` uses a simpler grid-sampling guard than `Compute`’s `DipolarNearFieldFilter` with `source_extent = 2.0 * geom.radius`. If `SampleShieldingAt` is intended to match production atom evaluation semantics, verify that this difference is intentional.
