# ring-susceptibility — claude review (readability focus)

- **Targets:** src/RingSusceptibilityResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**RingSusceptibilityResult.h** — Clear. The header comment lays out the tensor formula, its symmetry properties (asymmetric, non-traceless), and the McConnell-with-n̂ provenance before any code, so a physicist reads the physics before the API. The five-method class is self-explanatory.

**RingSusceptibilityResult.cpp** — Reads well as a coherent story: the kernel-math block carries an exemplary banner comment that names all three tensor terms and what irrep each contributes, then `Compute` walks atom → nearby rings → filter → store → accumulate in sensible order. The math is grouped into named stages (`f`, `K`, `M_over_r3`) matching the header. The only genuine friction is the cylindrical-coordinate `theta` computation, which is opaque and unlabeled, and a couple of cryptic local/helper names. No symbol soup; the through-line is followable in one pass.

---

## Findings

### 1. Coherent story / readability
- cpp:198 — `theta = atan2(d_plane.norm(), abs(z))` recomputes `d_plane.norm()` (already `rho` on line 197) and the `abs(z)` folds both hemispheres to one angle without explanation — as written, a reader must guess this is a polar tilt angle measured from the ring axis, sign-collapsed. — Reuse `rho` and add `// polar angle from ring axis (hemisphere-folded)`.
- cpp:90 — `(a == b ? 1.0 : 0.0)` is the Kronecker δ inline; appears twice (also :81) and reads as a magic ternary mid-expression. — Acceptable given the banner, but a `// δ_ab` trailing comment on each closes the gap.

### 2. Naming carries meaning
- cpp:46 `M_over_r3` — good, units-bearing. cpp:48 `f` — matches header "susceptibility scalar" and is documented; acceptable. cpp:47 `K` — single letter, but the inline comment names it; fine.
- cpp:179 `rn` — bare abbreviation for the `RingNeighbourhood*` reused across ~30 lines; mildly taxing. — Rename `neighbourhood` or `ring_nb`.
- cpp:61,194 `d` — reused as two different displacement vectors in the same function scope region (kernel helper line 61 vs Compute line 194); each is local so no bug, but the repeat invites confusion. — In Compute, name it `to_atom` to distinguish from the kernel's `d`.

### 3. Visible math structure (grouping)
- Clean. The three derived quantities are sequenced f → K → M_over_r3 with section comments matching the header's term-by-term breakdown; the build/decompose split (Mat3 accumulate then `SphericalTensor::Decompose`) is obvious at cpp:213/218.
- cpp:193-201 — the cylindrical-frame block is a named step ("Cylindrical coordinates in ring frame") but the z/rho/theta trio is computed inline; the `theta` line is the one unnamed-intermediate snag (see axis 1).

### 4. Function / method naming
- `ComputeRingChiKernel` says what it computes and returns (the kernel struct); good. `SampleShieldingAt` / `WriteFeatures` clear.
- cpp:253 `PackST_RS` — the `_RS` suffix is a file-disambiguation hack and `PackST` abbreviates SphericalTensor; reads as machine shorthand. — `PackSphericalTensor` (file-static, no collision risk).

### 5. Comments as signposts
- Strong overall: the cpp:25-43 banner is a model signpost (names each term and its irrep). cpp:74,77,83 inline labels are terse and correct.
- cpp:124-128 — the filter-set rationale comment is good and load-bearing (explains why topology check supplements distance at the boundary).
- cpp:176-178 — "may already have an entry from another calculator (BiotSavart)" is a useful cross-calculator signpost; keep.
- No verbose/AI-style or stale comments found.

### 6. Correctness (secondary)
- cpp:241-243 vs filter path — `SampleShieldingAt` reimplements the accept logic as three inline distance checks (singularity, `< radius`, `> cutoff`) instead of the `KernelFilterSet` used in `Compute`; check that this intentional simplification matches the production kernel (no `RingBondedExclusionFilter` applies at arbitrary grid points, so likely fine, but the `< radius` reject is a different boundary than `DipolarNearFieldFilter` uses).
- cpp:208-210 — `chi_tensor`/`chi_spherical`/`chi_scalar` are overwritten per matching ring without accumulation, while `M_total` accumulates; confirm `RingNeighbourhood` is intended to hold per-ring (not summed) values — appears intentional (per-neighbour record vs per-atom sum) but worth a one-line note.
- cpp:72 `cos_theta = d_hat.dot(ring_normal)` — assumes `ring_normal` is unit-length; check `RingGeometry::normal` is normalized at construction (no guard here).
- No sign/transpose/bounds issues spotted; the singularity guard (cpp:64) precedes all division and the `r==0` early return zeroes cleanly.
