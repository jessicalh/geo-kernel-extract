# bonded-energy — codex review (readability focus)

- **Targets:** src/BondedEnergyResult.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 17,441
- **Brief:** `../_brief.md`

---

**Verdict**

`BondedEnergyResult.h` mostly tells a coherent story: purpose, inputs, outputs, and units are visible. The main readability break is the compressed `atoms[5]` / `p[3]` encoding, which forces the reader to carry a type-dependent schema in their head.

`BondedEnergyResult.cpp` has a clear top-level shape: energy kernels, per-frame accumulation, feature output. The story breaks inside the dihedral/CMAP math and in `Compute`, where meaningful physical quantities collapse into `p0`, `p1`, `ix`, `p[]`, `fi`, `wj`, and `target`.

**1. Coherent Story / Readability**

- `src/BondedEnergyResult.h:42` — `atoms[5]` is compact but schema-dependent. As written, a reader must remember the atom arity and ordering for each interaction type. Suggest a short arity/order comment per type, or named local unpacking at use sites.
- `src/BondedEnergyResult.h:46` — `p[3]` is the biggest story break in the header. As written, every use of `ix.p[0]` requires jumping back here. Suggest `params[3]` at minimum, plus named constants or local aliases like `equilibrium_distance_nm`, `force_constant`, `multiplicity`.
- `src/BondedEnergyResult.cpp:45` — `DihedralAngle` is mathematically dense. As written, a reader must reconstruct “bond vectors → plane normals → signed angle” from `b1`, `b2`, `b3`, `n1`, `n2`, `m`. Suggest signpost comments: `// bond vectors`, `// plane normals`, `// signed angle`.
- `src/BondedEnergyResult.cpp:91` — CMAP has a good overall sequence, but `dx`, `fi`, `fj`, `wi`, `wj` hide units and axes. Suggest names like `grid_step_rad`, `phi_grid`, `psi_grid`, `phi_weight`, `psi_weight`.
- `src/BondedEnergyResult.cpp:136` — `ix` makes the main loop read mechanically rather than physically. Suggest `interaction`.
- `src/BondedEnergyResult.cpp:144` — `energy` plus `target` is workable but abstract. Suggest `interaction_energy` and `energy_channel` or `per_atom_energy`.

**2. Naming Carries Meaning**

- `src/BondedEnergyResult.h:64` — `cmap_grid_spacing` sounds like a physical spacing, but it is used as a point count. Suggest `cmap_grid_size` or `cmap_point_count`.
- `src/BondedEnergyResult.h:87` — `UBEnergy` relies on abbreviation. Suggest `UreyBradleyEnergy`.
- `src/BondedEnergyResult.h:90` — `ProperDihEnergy` / `ImproperDihEnergy` are understandable but clipped. Suggest `ProperDihedralEnergy` and `ImproperDihedralEnergy`.
- `src/BondedEnergyResult.cpp:133` — `count_ub` is less legible than the log output. Suggest `count_urey_bradley`.
- `src/BondedEnergyResult.cpp:196` — `cmap_idx` is clear enough locally; no issue.

**3. Visible Math Structure**

- `src/BondedEnergyResult.cpp:48` — the dihedral computation needs visible grouping. Suggest separating with comments or blank-line stages: bond vectors, plane normals, signed angle.
- `src/BondedEnergyResult.cpp:95` — CMAP mapping and interpolation are structurally present but compressed. Suggest named intermediate stages: `// angle grid coordinates`, `// periodic cell`, `// bilinear weights`.
- `src/BondedEnergyResult.cpp:147` — each switch case hides parameters behind `ix.p[]`. Suggest local unpacking per case, e.g. `r0_nm`, `k_bond`, `theta0_rad`, `multiplicity`, before calling the evaluator.

**4. Function / Method Naming**

- `src/BondedEnergyResult.cpp:38` — `EvalUB` is terse. Suggest `EvalUreyBradleyEnergy`.
- `src/BondedEnergyResult.cpp:65` — `EvalProperDih` says the type but not clearly that it returns energy. Suggest `EvalProperDihedralEnergy`.
- `src/BondedEnergyResult.cpp:73` — same for `EvalImproperDih`; suggest `EvalImproperDihedralEnergy`.
- `src/BondedEnergyResult.cpp:86` — `EvalCMAP` is acceptable but could be clearer as `EvalCmapEnergy`.
- `src/BondedEnergyResult.h:79` — `Compute` is vague in isolation, but acceptable as a `ConformationResult` factory if this is the local convention.

**5. Comments As Signposts**

- `src/BondedEnergyResult.cpp:12` — section banner is clear but decorative. If touched, use `// energy kernels`.
- `src/BondedEnergyResult.cpp:77` — comment is useful. Terse replacement: `// periodic difference`.
- `src/BondedEnergyResult.cpp:94` — useful but slightly wordy. Terse replacement: `// angle grid coordinates`.
- `src/BondedEnergyResult.cpp:99` — “sufficient — CMAP grids are smooth” is a justification, not a signpost, and may invite scientific debate. Terse replacement: `// bilinear CMAP`.
- `src/BondedEnergyResult.cpp:137` — restates the guard. Terse replacement: `// index guard`.
- `src/BondedEnergyResult.cpp:216` — clear, but could be tighter: `// per-atom split`.
- `src/BondedEnergyResult.cpp:222` — clear; no change needed.

**6. Correctness**

- `src/BondedEnergyResult.cpp:139` — the index guard validates only `ix.n_atoms`, but the switch later reads fixed arity by type. If malformed parameters can reach this point, a dihedral with `n_atoms < 4` or CMAP with `n_atoms < 5` can use default atom slots and divide by the wrong arity. Suggest asserting or guarding type-specific arity.
- `src/BondedEnergyResult.cpp:89` — CMAP checks `grid.empty()` and `spacing < 2`, but not `grid.size() >= spacing * spacing`. If malformed grids are possible, this can index out of bounds. Suggest an invariant assertion or guard.
- `src/BondedEnergyResult.cpp:79` — `M_PI` is non-standard C++ unless provided by the platform. If portability matters, check this build configuration.
- `src/BondedEnergyResult.h:73` / `src/BondedEnergyResult.h:79` — the header uses `std::type_index` and `std::unique_ptr` without directly including `<typeindex>` and `<memory>`. This may compile transitively, but include hygiene is fragile.
- No definite physics/sign-convention bug from the inline text. I would specifically regression-check `DihedralAngle` against the GROMACS convention, but I am not calling it wrong from this excerpt.
