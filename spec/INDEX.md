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

### Historical (learn/bones/)
Session notes, early specs, analysis passes, resolved feedback.
Reference only.

### Source
- **src/** — 54 headers, ~14.5K lines
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
