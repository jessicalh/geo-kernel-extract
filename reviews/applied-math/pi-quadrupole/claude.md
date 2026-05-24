# pi-quadrupole — claude review (readability focus)

- **Targets:** src/PiQuadrupoleResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**PiQuadrupoleResult.h** — Reads cleanly. The header comment is an exemplary signpost: it states the physical model, gives the closed-form kernel, names every symbol (`d`, `r`, `dn`), states the three structural properties (symmetric / traceless / pure-T2), and flags the units. A physicist can read it once and know exactly what the file computes and why. No findings of substance.

**PiQuadrupoleResult.cpp** — Mostly coherent and follows the header's story. `ComputePiQuadKernel` is the core and is well-signposted; the kernel formula is restated adjacent to the loop that implements it. The story does break in two places: the `RingNeighbourhood` find-or-create block inside `Compute` mixes bookkeeping with a second, unexplained set of cylindrical-coordinate computations (`z`, `rho`, `theta`) that have no comment tying them to the physics; and a cluster of magic numbers (`8`, `5`, `40`, `9`) carry array-layout meaning that the reader must infer. Correct, but the cylindrical block and the literal dimensions are the spots a fresh reader stalls.

---

## Findings

### 1. Coherent story / readability

- `PiQuadrupoleResult.cpp:181-205` — the find-or-create `RingNeighbourhood` block interrupts the quadrupole-kernel story with unrelated cylindrical-coordinate setup (`z`/`rho`/`theta`); as written, a reader must guess these feed *other* ring calculators, not this one — add a one-line signpost `// shared ring-neighbour geometry (used by other calcs)`.
- `PiQuadrupoleResult.cpp:201` — `theta = atan2(d_plane.norm(), std::abs(z))` silently folds both half-spaces onto a single polar angle via `abs(z)`; the intent (axial symmetry) is invisible — note `// fold to [0,pi/2] (axial symmetry)`.
- `PiQuadrupoleResult.cpp:214` `if (ti >= 0 && ti < 8)` and `:216` `for c < 5` — bare `8` (ring types) and `5` (T2 components) read as magic; name or comment them so the per-type accumulation shape is self-evident.

### 2. Naming carries meaning

- `PiQuadrupoleResult.cpp:74` — `dn` is the height above the ring plane (d·n); the inline comment saves it, but `d_dot_n` or `height` would carry the meaning without the comment crutch.
- `PiQuadrupoleResult.cpp:88` — `diag_term` is fine but slightly understates it; `// isotropic delta_ab piece` already in the formula comment covers it — acceptable as-is.
- `PiQuadrupoleResult.cpp:182,189-204` — `rn` / `new_rn` are terse but the type `RingNeighbourhood` is right there; acceptable.
- `PiQuadrupoleResult.cpp:269` — `N` for atom count is clear in context; fine.

### 3. Visible math structure (grouping)

- `PiQuadrupoleResult.cpp:90-99` — the kernel is a single fused 5-term expression per (a,b); the adjacent comment block (`:83-87`) maps it well, so the fusion is legible. Clean.
- `PiQuadrupoleResult.cpp:69-72` — the `r2/r5/r7/r9` power ladder is grouped and ordered; reads as a deliberate precompute step. Clean.
- `PiQuadrupoleResult.cpp:196-201` — cylindrical projection (`d_vec`, `z`, `d_plane`, `rho`, `theta`) is an unnamed sub-step embedded mid-bookkeeping; grouping it under a one-line label (see 1) would make the shape obvious.

### 4. Function / method naming

- `PiQuadrupoleResult.cpp:261` — `PackST_PQ` packs a SphericalTensor into a flat 9-double layout (T0,T1,T2); name is adequate but the `_PQ` suffix implies PQ-specificity when the packing is generic — minor, leave it.
- `ComputePiQuadKernel`, `Compute`, `SampleShieldingAt`, `WriteFeatures` — all say what they compute and return; clean.

### 5. Comments as signposts

- `PiQuadrupoleResult.cpp:25-44` — the duplicated kernel/properties block restates the header verbatim; acceptable as a local reference, but consider trimming to `// kernel G_ab = T_abcd n_c n_d — see header` to avoid two copies drifting.
- `PiQuadrupoleResult.cpp:130-134` — good, grounded signpost (why the topology filter matters more for quadrupole than dipole).
- `PiQuadrupoleResult.cpp:196-201` — non-obvious cylindrical block is unlabeled (see 1/3); needs the 2-4 word signpost.
- `PiQuadrupoleResult.cpp:42` — `// Traceless: yes (Laplace; verified numerically in tests)` is a good grounded signpost.

### 6. Correctness (secondary)

- `PiQuadrupoleResult.cpp:248-251` — `SampleShieldingAt` guards with `distance < geom.radius` (skip if inside ring radius) but `Compute` uses the full filter set (`MinDistance` + `DipolarNearField` + `RingBondedExclusion`); the two paths apply *different* near-field criteria, so grid samples and per-atom values won't agree near rings — confirm this divergence is intended.
- `PiQuadrupoleResult.cpp:64` — `ComputePiQuadKernel` returns a default (zero) result on the singularity guard but still sets nothing on `distance`/`direction`; the caller at `:164` reads `kernel.distance` (0.0) into the filter context — check that a zeroed distance can't pass `DipolarNearFieldFilter` and produce a spurious zero accumulation. Likely benign (G is zero) but worth confirming.
- `PiQuadrupoleResult.cpp:214` — `ti < 8` silently drops ring types with `TypeIndexAsInt() >= 8` from per-type accumulation while still adding them to `G_total` (`:221`); if any ring type index can exceed 7, the per-type NPYs and the total tensor disagree — confirm 8 is the true type ceiling.
- `PiQuadrupoleResult.cpp:208-209` — `rn->quad_tensor`/`quad_spherical` are overwritten (not accumulated) per ring; correct since one `rn` is unique per ring index, but note the find-or-create at `:183-188` relies on that uniqueness — fine as written.
