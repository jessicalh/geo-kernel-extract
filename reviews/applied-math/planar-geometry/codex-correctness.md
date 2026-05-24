# planar-geometry — codex review

- **Targets:** src/PlanarGeometryResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 41216
- **Brief:** `../_brief.md`

---

**Algorithm Correctness**
src/PlanarGeometryResult.cpp:40 — `Dihedral` has no guard for zero/near-zero plane normals; collinear triples can make an undefined dihedral return `0` via `atan2(0,0)` — guard `b2`, `n1`, and `n2`; return `kNaN`.

src/PlanarGeometryResult.cpp:187 — sorting pyramidalization neighbours by atom index makes the sign depend on numbering, not a topology/improper order; this conflicts with the stated CHARMM/AMBER sign convention — use a canonical chemical order or document the sign as arbitrary.

src/PlanarGeometryResult.cpp:75 — degenerate neighbour plane returns `0.0`, indistinguishable from exact planarity — return a missing/degenerate sentinel if the ABI allows, or label this as degenerate-zero behavior.

src/PlanarGeometryResult.cpp:35 — I don’t know: dihedral sign should be checked against the library’s canonical dihedral/AMBER convention; this implementation uses the `b1` arm, while some standard formulas use `-b1` — check one four-point sign fixture.

src/PlanarGeometryResult.h — no additional correctness findings.

**Variable Naming Clarity**
src/PlanarGeometryResult.cpp:70 — `A/B/C/D/G/n` obscure central atom, neighbour plane, centroid, and normal roles — `central_pos`, `neighbor*_pos`, `neighbor_centroid`, `plane_normal`.

src/PlanarGeometryResult.cpp:142 — `cs/sn` hide the Cremer-Pople mode projection — `q2_cos_sum`, `q2_sin_sum`.

src/PlanarGeometryResult.cpp:264 — `nb` hides that these are atom indices — `neighbor_indices`.

src/PlanarGeometryResult.cpp:348 — `ri` is a ring index after earlier meaning residue index — `arom_ring_i`.

src/PlanarGeometryResult.cpp:377 — `ri` repeats the same ambiguity for saturated rings — `sat_ring_i`.

src/PlanarGeometryResult.h — variable names are clear.

**Grouping Of Math Operations**
src/PlanarGeometryResult.cpp:35 — dihedral setup, plane normals, and signed-angle construction are packed together without visible math stages — group as bond vectors, plane normals, signed angle.

src/PlanarGeometryResult.cpp:142 — Cremer-Pople projection is hidden in a terse accumulator loop — make the unscaled Q2 cosine/sine projection a named step.

src/PlanarGeometryResult.h — clean on math grouping.

**Method / Function Naming**
src/PlanarGeometryResult.cpp:70 — `Pyramidalization` sounds like a general chemistry quantity or angle, but returns signed displacement in Å — `PyramidalizationDisplacement`.

src/PlanarGeometryResult.cpp:100 — `CremerPople5Ring` is acceptable but slightly underspecified — `CremerPople5RingQTheta`.

src/PlanarGeometryResult.h — public method names are clear within the result interface.

**Comments**
src/PlanarGeometryResult.h:2 — file header is process/history-heavy — `// planar geometry features`.

src/PlanarGeometryResult.h:6 — cross-result block is dated process detail — `// cross-result readers`.

src/PlanarGeometryResult.h:25 — convention block is too large for inline comments — `// geometry conventions`.

src/PlanarGeometryResult.cpp:30 — decorative banner and explanation are too verbose — `// signed dihedral`.

src/PlanarGeometryResult.cpp:49 — wrap comment is history-heavy — `// wrap to pi`.

src/PlanarGeometryResult.cpp:60 — pyramidalization comment is verbose, and `A == G` is not required for planarity — `// signed plane distance`.

src/PlanarGeometryResult.cpp:107 — CP normal comment mixes equations, rationale, and review history — `// CP mean normal`.

src/PlanarGeometryResult.cpp:153 — degeneracy comment is citation/history-heavy — `// planar-ring guard`.

src/PlanarGeometryResult.cpp:168 — neighbour helper comment overexplains implementation detail — `// neighbor endpoints`.

src/PlanarGeometryResult.cpp:187 — sorting comment claims chemistry portability despite atom-index-defined sign — `// atom-index order`.

src/PlanarGeometryResult.cpp:278 — omega block is long policy/history text — `// peptide omega`.

src/PlanarGeometryResult.cpp:465 — mask write comment overexplains consumer interpretation — `// X-Pro mask`.
