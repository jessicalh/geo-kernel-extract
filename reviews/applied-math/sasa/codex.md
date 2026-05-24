# sasa — codex review (readability focus)

- **Targets:** src/SasaResult.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 6,952
- **Brief:** `../_brief.md`

---

**Verdict**

`src/SasaResult.h`: Mostly coherent. The header tells the SASA story clearly, but it promises stored result-like access while the class actually reads values back from `ProteinConformation`; that ownership model is not visible here.

`src/SasaResult.cpp`: The main algorithm is followable top to bottom: generate sphere points → expand atom radius → test occlusion → compute SASA and normal → write outputs. It is not symbol soup, but several small compressions make the reader infer physical meaning from variable names and comments, especially around the neighbor cutoff and surface-normal construction.

**1. Coherent Story / Readability**

`src/SasaResult.h:38` — `AtomSASA` / `AllSASA` imply `SasaResult` owns SASA values, but the only member is `conf_`; as written, a reader must discover later that results are stored on atoms, not in this result object — suggestion: add a terse comment near `conf_`, e.g. `// values stored on atoms`.

`src/SasaResult.cpp:62-65` — the neighbor search logic is physically sensible but compressed: `max_vdw`, `search_radius`, and the extra `probe_radius` require mental reconstruction of “all atoms whose expanded probe spheres could cover this test sphere” — suggestion: rename `max_vdw` to `max_occluder_vdw_radius` and add `// occluder search radius`.

`src/SasaResult.cpp:69-85` — the inner loop is readable, but the narrative jumps from sampling points to accumulating normals without a visible label for the exposed-point stage — suggestion: add a short block label before line 69, e.g. `// exposed surface samples`.

`src/SasaResult.cpp:88-97` — SASA assignment and normal normalization are close enough, but they are two conceptual outputs fused into one unlabeled end block — suggestion: split with comments `// accessible area` before line 88 and `// exposed normal` before line 91.

`src/SasaResult.cpp:100-107` — `radii_source` is recorded as numeric `0.0` with unit `"Bondi_1964"`; as written, a reader must infer this is metadata forced through a numeric API — suggestion: add `// metadata-only choice` or, if supported elsewhere, use a string metadata field.

**2. Naming Carries Meaning**

`src/SasaResult.cpp:19` — `n` is acceptable locally, but `FibonacciSphere(int n)` hides that this is a point count — suggestion: `point_count`.

`src/SasaResult.cpp:49` — `unit_sphere` sounds like a geometric object, but it holds sample directions — suggestion: `sphere_directions` or `sample_directions`.

`src/SasaResult.cpp:57` — `i` is fine as a loop index, but downstream names would read better if the atom role were explicit in derived quantities — suggestion: keep `i`, but prefer `atom_pos`, `expanded_atom_radius` over `pos_i`, `r_i`.

`src/SasaResult.cpp:63` — `max_vdw` lacks role and unit/frame — suggestion: `max_occluder_vdw_radius`.

`src/SasaResult.cpp:67` — `exposed` is understandable but slightly underspecified — suggestion: `exposed_count` or `exposed_samples`.

`src/SasaResult.cpp:72` — `occluded` is clear.

`src/SasaResult.cpp:118` — `AtomSASA` returns Angstroms squared but the name does not say units — acceptable in this code style; clearer would be `AtomSasaAngstrom2` only if the project uses unit-bearing names elsewhere.

**3. Visible Math Structure**

`src/SasaResult.cpp:23-27` — the Fibonacci lattice formula is compact and recognizable, but the vertical coordinate and azimuth are not named; as written, a reader must parse spherical-coordinate construction from the trig — suggestion: name intermediates such as `z`, `radius_xy`, and `azimuth`.

`src/SasaResult.cpp:70-77` — test-point generation, occluder radius expansion, distance calculation, and occlusion decision are clear but tightly packed — suggestion: add one signpost comment `// occlusion test`; no algorithm change needed.

`src/SasaResult.cpp:89` — area computation fuses sphere area, exposed fraction, and assignment — suggestion: introduce `exposed_fraction` as a named intermediate if this code is meant to be read by non-maintainers.

`src/SasaResult.cpp:93-97` — normalization thresholding is structurally clear.

**4. Function / Method Naming**

`src/SasaResult.h:33` — `Compute` is conventional for this result framework, but by itself does not say it mutates atom fields — suggestion: if names can change, `ComputeAtomSasa` would be clearer; otherwise add a comment that outputs are stored on atoms.

`src/SasaResult.cpp:15` — `BondiRadius` is clear enough, though `BondiVdwRadius` already carries the fuller meaning — suggestion: the wrapper adds little narrative value; if kept, `BondiVdwRadiusForElement`.

`src/SasaResult.cpp:19` — `FibonacciSphere` is broadly understandable, but it returns directions rather than a sphere — suggestion: `FibonacciSphereDirections`.

`src/SasaResult.h:38-39` — `AtomSASA` and `AllSASA` are clear at call sites.

**5. Comments As Signposts**

`src/SasaResult.h:3-19` — clean and useful; it gives method, parameters, and outputs without overexplaining.

`src/SasaResult.cpp:54` — “GeometryChoice: record the parameters used” is process prose rather than a math signpost — suggestion: `// parameter provenance`.

`src/SasaResult.cpp:62` — useful, but slightly wordy — suggestion: `// occluder search radius`.

`src/SasaResult.cpp:91-92` — good scientific signpost; could be terser — suggestion: `// exposed normal`.

`src/SasaResult.cpp:100` — “Record a single GeometryChoice summarising...” is verbose relative to the block — suggestion: `// SASA provenance`.

`src/SasaResult.cpp:123` — useful because the storage model is surprising; keep it, or move the same idea to the header too.

`src/SasaResult.cpp:135` — restates filename/shape already visible in adjacent code, but harmless — suggestion: keep if this is the project’s feature-writing convention.

`src/SasaResult.cpp:141` — useful output signpost; acceptable.

**6. Correctness**

`src/SasaResult.cpp:47-49,89` — if `sasa_n_points` can be `0` or negative, the code can divide by zero or create an invalid vector size; I do not know whether config validation rejects this elsewhere — suggestion: check config validation requires `sasa_n_points > 0`.

`src/SasaResult.cpp:63` — comment says sulfur is the largest Bondi radius in the table; I do not know whether `Element::S` is still maximal in `PhysicalConstants.h` — suggestion: check the table if more elements have been added.

No definite physics/sign/frame bug is visible from the inlined code.
