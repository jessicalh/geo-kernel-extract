# Document Index — NMR Shielding Tensor Prediction System

**READ THIS FIRST** when starting a new context.

## Authority Hierarchy

1. **spec/CONSTITUTION.md** — supreme constraint on principles,
   sign conventions, anti-simplification rules, inviolable rules.
2. **OBJECT_MODEL.md** — canonical concrete specification.
   Every class, property, type, unit, dependency. When the constitution
   and object model disagree on concrete types, OBJECT_MODEL.md wins.
3. **EXTRACTION_ORDER.md** — what each ConformationResult computes,
   dependency graph, feature manifest (189 features).
4. **CALCULATOR_PARAMETER_API.md** — the interface between the
   ParameterCorrectionResult (e3nn model) and the 8 classical calculators.
   93 tuneable parameters with equations, defaults, references.
5. **spec/LAYER0_PLAN.md** — implementation passes. What gets built
   in what order.
6. **PATTERNS.md** — how to implement things. Anti-patterns. C++ rules.
   Updated after each correction pass. Living document.
7. **spec/MATHS_GOALS.md** — mathematical validation plan. Three pillars:
   T0 corrections, T2 residual analysis, classical calculation correctness.
8. **GEOMETRIC_KERNEL_CATALOGUE.md** — mathematical foundation for
   the 8 classical calculators. McConnell full tensor derivation. Build
   order and verification strategy.

## Reading Order for a New Context

1. This file (spec/INDEX.md)
2. Session state in Claude memory (captures latest project state)
3. spec/CONSTITUTION.md — skim for principles, read sign convention
   and inviolable rules carefully
4. OBJECT_MODEL.md — the concrete types. This is what you code against.
5. spec/LAYER0_PLAN.md — where we are in implementation
6. PATTERNS.md — how to do things safely, what gets code rejected
7. spec/MATHS_GOALS.md — what must be mathematically correct
8. GEOMETRIC_KERNEL_CATALOGUE.md — calculator physics and build order
9. APPLIED_MATHS_FIXES.md — numerical stability notes per calculator

Do NOT read FEEDBACK.md first — most items are resolved.

## All Documents (in this repo)

### Top-level working documents
- **OBJECT_MODEL.md** — canonical concrete specification
- **EXTRACTION_ORDER.md** — dependency graph + feature manifest
- **CALCULATOR_PARAMETER_API.md** — model-to-calculator parameter interface
- **PATTERNS.md** — agent implementation guide (living doc)
- **GEOMETRIC_KERNEL_CATALOGUE.md** — calculator physics, McConnell derivation
- **APPLIED_MATHS_FIXES.md** — numerical stability pad per calculator
- **FEEDBACK.md** — resolved contradictions (historical)

### spec/ (design principles + plans)
- **INDEX.md** — THIS FILE
- **CONSTITUTION.md** — supreme constraint (~1200 lines)
- **LAYER0_PLAN.md** — implementation passes
- **MATHS_GOALS.md** — mathematical validation plan (3 pillars)
- **USE_CASES.md** — the 4 use cases (agreed 2026-04-03)
- **BUILDER_AND_PIPELINE_SPEC.md** — builder + OperationRunner design
- **ASSESSMENT_20260401.md** — honest assessment at end of Layer 0
- **PROTONATION_DESIGN_HISTORY.md** — design vs reality for protonation
- **SESSION_20260403_FLEET_LOADER.md** — fleet loader session notes
- **meta-docs-review/** — documentation audit artifacts (2026-04-03).
  For the human author, not for working agents. Contains critiques,
  proposed edits, and a Gemini review triage.
- **CATEGORICAL_ANALYSIS_PASS1.md** — 470 ops from v1 extractor
- **DESIGN_REVIEW_PRELIM.md** — two-pass review (historical)
- **DESIGN_BRIEF_DRAFT_1.md** — original design brief (historical)
- **DEPENDENCIES.md** — external library list
- **REWRITE_DECISIONS.md** — why we rewrote (historical)
- **UNIFICATION_BUGS.md** — bugs found in the old project
- **DIRECTORY_SET.md** — directory structure (historical)
- **DEPENDENCY_ANALYSIS.md** — dependency analysis (historical)

### Source and tests
- **src/** — source code (~90 files, 14.5K lines)
- **tests/** — GTest suite (283 tests, 40 test files)
- **extern/** — header-only libraries (nanoflann, sphericart)
- **CMakeLists.txt** — build system (links libgromacs, OpenBabel,
  cifpp, libdssp, APBS bridge, Eigen3)

### External (not in this repo)
- **/mnt/extfast/2026thesis/biot-savart/** — old project (REFERENCE ONLY)
- **/mnt/extfast/2026thesis/consolidated/** — 734 WT+ALA protein pairs
- **/mnt/extfast/fleet_results/** — GROMACS CHARMM36m fleet data

## Current Status (2026-04-06)

**Layer 0 complete + 10 calculators (8 classical + 2 MOPAC-derived) +
MOPAC electronic structure + builder unification + CLI entry point.**
260 tests (3 KaML skipped).

### What exists in code

**Three builders** returning BuildResult (protein + charges + net_charge):
- BuildFromPdb: protonates with linked reducelib (Richardson lab),
  assigns ff14SB charges. For bare PDB input.
- BuildFromOrca: charges from AMBER prmtop. For DFT data.
- BuildFromGromacs: charges from GROMACS TPR via libgromacs.
  For MD ensemble data.

**OperationRunner**: single home for all ordered sequences.
Run (standard), RunMutantComparison, RunEnsemble. Replaces
both Pipeline and CalculationRunner. MOPAC and APBS run
automatically when charges are available.

**nmr_extract**: CLI entry point for all 4 use cases
(--pdb, --orca, --mutant, --fleet). See spec/USE_CASES.md.

**8 classical geometric kernel calculators**, batch-validated on
723 proteins: McConnell, Coulomb, RingSusceptibility, HBond,
BiotSavart, HaighMallion, PiQuadrupole, Dispersion. Full
Mat3+SphericalTensor output. T2 independence verified.

**2 MOPAC-derived calculators** using conformation electronic
structure from MopacResult: MopacCoulombResult (EFG from QM
charges), MopacMcConnellResult (bond anisotropy weighted by
Wiberg bond order). Both produce independent T2 angular features.

**MOPAC electronic structure** (MopacResult): PM7+MOZYME
semiempirical calculation provides per-atom Mulliken charges,
s/p orbital populations, and per-bond Wiberg bond orders for
every conformation. Runs as a conformation electronic structure
precondition before calculators that depend on it.

**Infrastructure**: CovalentTopology (geometry-to-topology boundary),
WriteFeatures (46 NPY arrays per conformation), KernelFilterSet
with 4 filters (DipolarNearField, SelfSource, SequentialExclusion,
RingBondedExclusion), typed conformation hierarchy (Crystal,
Prediction, MDFrame, Derived).

**External tools**: DSSP (libdssp), APBS (PB solve, first-class),
MOPAC (PM7+MOZYME, conformation electronic structure), reduce
(linked C++ library for protonation), OpenBabel (bond perception).

### What does NOT yet exist in code

- ParameterCorrectionResult (e3nn model)
- FeatureExtractionResult (feature manifest as ConformationResult)
- PredictionResult (ML inference)
- Copy-and-modify / rebuild pattern for re-protonation
- Filter TOML configuration
- Batch test migration to OperationRunner (4 files, manual chains)

### Doc coherency notes (2026-04-03)

Critique documents exist for the three foundational docs:
- spec/CONSTITUTION_CRITIQUE.md
- spec/OBJECT_MODEL_CRITIQUE.md
- spec/PATTERNS_CRITIQUE.md

These identify stale claims, missing descriptions, and aspirational
sections that need collaborative review before editing the originals.
The physics (geometric kernels, sign conventions, minimum representation
contract, filter framework) is accurate. The plumbing descriptions
(copy semantics, framework stores, named accessors) need updating.
