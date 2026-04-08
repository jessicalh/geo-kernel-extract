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
   classical calculators. McConnell tensor from first principles.

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
- **USE_CASES.md** — the 4 use cases (--pdb, --orca, --mutant, --fleet)
- **GEOMETRY_CHOICE_BRIEF.md** — GeometryChoice recording spec
- **DEPENDENCIES.md** — external library list
- **DIRECTORY_SET.md** — directory structure (historical)
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

## Current Status (2026-04-07)

10 calculators (8 classical + 2 MOPAC-derived), MOPAC electronic
structure, three builders, OperationRunner, CLI entry point.
287 tests pass, 0 failures.

### Next: calibration pipeline

~80 calculator parameters externalised to TOML, swept against
DFT WT-ALA deltas across 723 proteins. The measure is T2 R² —
how well calibrated classical kernels reproduce the angular
structure of DFT shielding. This is what makes the system a
calibrated scientific instrument.

Feature extraction is distributed: each ConformationResult implements
WriteFeatures(), and ConformationResult::WriteAllFeatures() traverses
all attached results. This exports the full tensor data (geometric
kernels + DFT reference when available) as NPY arrays for the
calibration pipeline.
