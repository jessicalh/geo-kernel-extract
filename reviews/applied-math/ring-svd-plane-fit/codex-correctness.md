# ring-svd-plane-fit — codex review

- **Targets:** src/Ring.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 15254
- **Brief:** `../_brief.md`

---

**1. Algorithm Correctness**
src/Ring.cpp:8 — no concrete correctness issue found in the centroid/SVD/radius computation as written. I don’t know from these files alone whether `atom_indices` is always cyclic; check ring construction because lines 34–39 rely on the first three vertices for normal orientation.

**2. Variable Naming Clarity**
src/Ring.h:40 — `total_B_at_center` gives quantity/location but not units/frame — suggest `total_field_center`.
src/Ring.h:42 — `total_G_T0_diagnostic` is opaque; unclear what `G` and `T0` mean — suggest `total_G_T0` or `G_T0_sum` only if those symbols are documented elsewhere.
src/Ring.cpp:27 — `coords` hides that these are centered vertex coordinates — suggest `centered_vertices`.

**3. Grouping Of Math Operations**
src/Ring.cpp:17 — centroid step is clear.
src/Ring.cpp:23 — SVD plane-fit step is clear.
src/Ring.cpp:34 — orientation step is clear.
src/Ring.cpp:42 — radius step is clear.

**4. Method / Function Naming**
src/Ring.h:75 — `ComputeGeometry` is accurate for center, normal, radius, vertices.
src/Ring.h:221 — `CreateRing` is accurate.

**5. Comments**
src/Ring.h:3 — verbose file banner comment — suggest `// ring types`.
src/Ring.h:5 — explanatory/stale-prone design prose — suggest remove or `// type constants`.
src/Ring.h:9 — stale: says 8 types, but code defines 9 including `ProPyrrolidineRing` — suggest `// ring classes`.
src/Ring.h:178 — verbose chemistry explanation block — suggest `// saturated ring`.
src/Ring.cpp:23 — comment is longer than needed but grounded — suggest `// best-fit plane`.
src/Ring.cpp:34 — clear but can be shorter — suggest `// orient normal`.
src/Ring.cpp:65 — verbose failure-policy comment — suggest `// invalid enum guard`.
