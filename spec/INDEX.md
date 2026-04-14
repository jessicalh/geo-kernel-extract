# Document Index — NMR Shielding Tensor Prediction System

**READ THIS FIRST** when starting a new context.

## Reading Order

### Tier 1 — Start here (what the system is)

1. **doc/ARCHITECTURE.md** — structure, provenance paths, pipeline,
   kernel equations, diagrams. Points at code with line numbers.
2. **PATTERNS.md** — what holds the system together, what breaks it.
   Anti-patterns, numerical stability, C++ rules, lessons learned.

### Tier 2 — Open while working (field-level detail)

3. **OBJECT_MODEL.md** — every class, property, type, unit. The
   canonical concrete specification. Code against this.
4. **EXTRACTION_ORDER.md** — dependency graph, feature manifest
   (189 features), what each ConformationResult computes.
5. **CALCULATOR_PARAMETER_API.md** — 93 tuneable parameters with
   equations, defaults, literature references.

### Tier 3 — Design decisions (ask before diving in)

These documents are deep. Reading them burns significant context
for significant gain — they ground you in the physics and sign
conventions that every calculator depends on. Ask the user whether
the task requires this depth. The user is aware of these tiers and
may direct you to dive deep on one or more — if so, read thoroughly.
Skimming these is worse than not reading them.

6. **spec/CONSTITUTION.md** — supreme constraint. Sign conventions,
   inviolable rules, anti-simplification rules, minimum output
   contracts per calculator.
7. **spec/MATHS_GOALS.md** — the three mathematical pillars: T0
   parameter corrections, T2 angular residual, classical calculation
   correctness. Validation checklist.
8. **GEOMETRIC_KERNEL_CATALOGUE.md** — full derivations for the 8
   classical + 2 MOPAC-derived calculators. McConnell tensor from
   first principles.

### When implementing a calculator

Read an existing calculator first (e.g., McConnellResult.cpp).
You must:
- Produce full rank-2 tensor output (Mat3 + SphericalTensor)
- Hold a KernelFilterSet and call AcceptAll() before every evaluation
- Record GeometryChoices via GeometryChoiceBuilder during Compute(),
  showing every inclusion, exclusion, and triggered event with
  entities, roles, outcomes, and named numbers with units
- Document the unit chain in comments
- Verify sign with a known physical scenario

See PATTERNS.md "Numerical Stability" section.

### When adding or changing WriteFeatures output

Validate with `./smoke_tests` (~80s), not individual calculator tests.
See **spec/TEST_FRAMEWORK.md** for the full test tier decision tree.

Any new NPY file must be registered in the Python reader SDK so
consumers can load it.  See **spec/EXTRACTION_SDK.md** for the
design and **python/API.md** for the full column-level reference.

The reader lives at `python/nmr_extract/`.  The format contract is
`python/nmr_extract/_catalog.py` — one `ArraySpec` per NPY file.
If you add a column to `ring_contributions.npy` or create a new
array, update `_catalog.py` and the corresponding wrapper class,
then run `python -m pytest python/tests/` to verify.

## All Documents

### Top-level
- **OBJECT_MODEL.md** — canonical concrete specification
- **EXTRACTION_ORDER.md** — dependency graph + feature manifest
- **CALCULATOR_PARAMETER_API.md** — model-to-calculator parameter interface
- **PATTERNS.md** — what holds the system together (living doc)
- **GEOMETRIC_KERNEL_CATALOGUE.md** — calculator physics, derivations

### doc/
- **ARCHITECTURE.md** — system map with diagrams
- **diagrams/** — mermaid sources + rendered SVG/PNG
- **generated/doxygen/** — full class UML, collaboration graphs

### spec/
- **INDEX.md** — this file
- **CONSTITUTION.md** — supreme constraint
- **MATHS_GOALS.md** — mathematical validation plan
- **USE_CASES.md** — the 5 use cases (--pdb, --orca, --mutant, --trajectory, --trajectory --analysis)
- **GEOMETRY_CHOICE_BRIEF.md** — GeometryChoice recording spec
- **DEPENDENCIES.md** — external library list
- **DIRECTORY_SET.md** — directory structure (historical)
- **TEST_FRAMEWORK.md** — test tiers, what to run when, smoke test design
- **ENSEMBLE_MODEL.md** — GromacsProtein trajectory-based model (revised 2026-04-12, replaces old EnsembleConformation design)
- **TRAJECTORY_EXTRACTION.md** — full-system trajectory path: explicit solvent E-field, water packing, ion field calculators (design, 2026-04-12)
- **TIMING_4876_ATOMS.md** — measured per-calculator times on 4876-atom protein (2026-04-12)
- **OUTSTANDING_GROMACS_PATH.md** — open items for trajectory extraction path (GeometryChoice, KernelFilterSet, TOML, SDK tests)
- **POLARISABILITY_ROADMAP_2026-04-13.md** — charge polarisation proxy approaches: SASA normal, HydrationGeometryResult, water-embedded AIMNet2, EEQ, E-field variance. What was rejected and why.
- **ANALYSIS_TRAJECTORY_2026-04-14.md** — analysis trajectory mode: exhaustive per-frame H5, per-ring K=6, ridge+MLP predictions, projections. 10-protein workspace design (tentative, evolving).
- **meta-docs-review/** — 2026-04-03 documentation audit artifacts

### ui/
- **CLAUDE.md** — UI-specific rules, library API quick-reference,
  current state, what you own vs read-only
- **UI_ROADMAP.md** — forward-looking per-calculator visualization plan
  (stable — nothing here is needed for correctness)
- **GEOMETRY_CHOICE_UI_NOTES_20260407.md** — library-to-UI handoff
  for GeometryChoice display (what each calculator records, UI ideas)
- **src/REST_INTERFACE_SPEC.md** — all 20+ REST commands with parameters
- **doc/generated/doxygen/** — viewer class diagrams, collaboration graphs

### python/ — Extraction Reader SDK
- **spec/EXTRACTION_SDK.md** — design intent, trust boundary, what it covers
- **python/API.md** — consumer reference: every type, field, column
- **python/nmr_extract/_catalog.py** — format contract (one entry per NPY)
- **python/doc/** — Sphinx autodoc source (build with `sphinx-build`)
- **python/tests/test_load.py** — 79 tests against real extractions

### Historical (learn/bones/)
Session notes, early specs, analysis passes, resolved feedback.
Reference only.

### Source
- **src/** — 54 headers, ~14.5K lines
- **ui/src/** — 14 viewer source files (Qt6/VTK)
- **tests/** — GTest suite (287 tests, 40 test files, ~2.5 hours with MOPAC)
- **extern/** — header-only libraries (nanoflann, sphericart)

### External
- **/shared/2026Thesis/consolidated/** — 723 WT+ALA protein pairs
- **tests/data/fleet** — GROMACS CHARMM36m ensemble

## Current Status (2026-04-13)

19 calculators (8 classical, 2 MOPAC-derived, AIMNet2, SASA,
WaterField, HydrationShell, GromacsEnergy, HydrationGeometry, EEQ,
plus foundation: Geometry, SpatialIndex, Enrichment). DSSP extended
with 8-class SS, H-bond energies, chi1-4. Calibration settled at
R²=0.818 (per-element ridge, 55 kernels).

Two-pass trajectory streaming: GromacsProtein (adapter +
accumulators) + GromacsFrameHandler (streaming XTC reader, PBC fix,
frame lifecycle). Tested end-to-end on 1ZR7_6721 (479 atoms, 3525
water, 25 ions). 74 NPY arrays registered in SDK catalog (added water_polarization,
eeq_charges, eeq_cn).

**H5 master file (2026-04-13):** HighFive integration complete.
WriteH5 produces one file per protein with per-frame positions
(T, N, 3), 48-column Welford rollup (mean + std per atom), per-bond
length statistics, and frame times. SDK `load_trajectory()` reads
it. WriteCatalog produces the same rollup as CSV. Accumulators
cover all 6 classical kernel T0+|T2|, AIMNet2 charge + EFG, APBS,
water E-field (magnitude + 3 components), hydration geometry, DSSP
phi/psi + chi1-4 + H-bond energy, bond angles, SASA + surface normal,
4 frame-to-frame derivative trackers, and transition counters (chi
rotamers, SS changes). Hard-fail error handling on AIMNet2, WaterField,
HydrationShell, GromacsEnergy.

**Solvent calculator compliance (2026-04-13):** WaterFieldResult and
HydrationShellResult now under the TOML + GeometryChoice umbrella.
WaterFieldResult has KernelFilterSet (MinDistanceFilter). 4 new TOML
parameters registered. SasaResult extended with surface normal output
(sasa_normal.npy). nmr_extract reads UDP logging config from
~/.nmr_tools.toml [logging] section.

### Next: polarisability calculators

See **POLARISABILITY_ROADMAP_2026-04-13.md** for 5 planned approaches
to charge polarisation. DONE: SASA surface normal (item 1),
HydrationGeometryResult (item 2, writes water_polarization.npy),
EEQ calculator (item 4, writes eeq_charges.npy + eeq_cn.npy).
Remaining: water-embedded AIMNet2 (item 3), E-field variance
routing (item 5).

**Next work:** Analysis trajectory mode (`--trajectory --analysis`).
See **spec/ANALYSIS_TRAJECTORY_2026-04-14.md** for the tentative
design: exhaustive per-frame H5 with ridge+MLP predictions and
projections onto Stage 1 eigenvectors. 10-protein workspace at
`/shared/2026Thesis/fleet_calibration/`.
