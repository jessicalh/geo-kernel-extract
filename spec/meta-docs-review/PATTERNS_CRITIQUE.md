# Patterns Critique — Code vs Document

**Date**: 2026-04-03
**Purpose**: Audit of PATTERNS.md against actual code. Same categories
as the other critiques.

---

## Lessons 1-5, 7-12: Core Patterns

**ALL TRUE.** Every core pattern is accurately reflected in code:

- Private ConformationAtom construction (friend ProteinConformation)
- Template Result<T>() access (the ONLY template pattern)
- Singleton guarantee (one result per type)
- Dependency declaration via Dependencies() returning type_index vector
- Compute() factory returning unique_ptr
- Fields typed by functional analysis (grouped by result type on
  ConformationAtom, not by physics category)
- Diagnostic error messages with state values
- Return codes, not exceptions (RunResult with error string)
- Equations in comments, Eigen in code
- Both representations always present (Mat3 + SphericalTensor on 50+
  fields — perfectly consistent)
- Ring type virtual interface (Intensity(), JBLobeOffset(), etc.)

These are the backbone of the document and they're all correct.

---

## Lesson 6: Shielding Contribution Contract

**PARTIALLY STALE.** The lesson says:

> "At the end of Attach(), if the ParameterCorrectionResult is present,
> the calculator calls SubtractCalculatorContribution with its shielding
> output."

ParameterCorrectionResult does not exist. SubtractCalculatorContribution
does not exist. The shielding_contribution fields ARE written by every
calculator (bs_shielding_contribution, mc_shielding_contribution, etc.)
and these are real, tested, and output to NPY. But the residual
subtraction mechanism is aspirational.

**Recommendation**: Keep the shielding contribution storage part (it's
true and load-bearing for WriteFeatures). Mark the residual subtraction
as "planned for ParameterCorrectionResult integration."

---

## Lessons 13-19: Protonation, Charges, Loaders

**TRUE.** Protonation flows through Protein, typed charge sources,
variant index contract, element verification at loading boundaries,
each protein complete on its own, FinalizeConstruction order,
full McConnell tensor.

---

## Lessons 20-24: Calculator Physics

**ALL TRUE and critical.** These lessons document the physics decisions
that made all 8 calculators work:

- Lesson 20: Calculators compute GEOMETRIC KERNELS, not parameterized
  output. G_tensor is the kernel; sigma = I * G is the physics.
- Lesson 21: 6-step calculator validation (formula → batch → magnitude
  → per-source → DFT proximity → T2 independence)
- Lesson 22: Sign convention sigma = I * G (G has minus sign)
- Lesson 23: BS-HM T2 redundancy (cos = 0.999)
- Lesson 24: Fused ring representation (TRP5+TRP6 vs TRP9)

---

## Anti-Patterns Section

**TRUE.** All listed anti-patterns are correctly avoided in code:
- No string traversals in calculators
- No cifpp/CCD outside loading boundary
- No adapter/wrapper/bridge classes
- No template metaprogramming beyond Result<T>()
- No shadow data structures
- No premature optimization (SIMD, pools, alignas)
- No exception hierarchies
- No utility/helpers/common namespaces

---

## C++ Rules Section

**TRUE.** C++17, unique_ptr for polymorphic, raw pointers for
non-owning back-refs, PascalCase types/methods, snake_case_ private,
UPPER_CASE constants.

---

## What's MISSING from PATTERNS.md

### 1. Pipeline vs CalculationRunner distinction

Two orchestration paths exist and agents need to know which to use:

- **Pipeline.h/RunClassicalCalculators()**: Flexible, takes
  PipelineOptions (charge source, run_xtb, run_apbs, skip_dssp).
  Used by GromacsEnsembleLoader's RunAllFrames() and test_full_pipeline.
- **CalculationRunner.h**: Named rigid sequences (RunPdbAnalysis,
  RunSingleDft, RunMutantComparison). Returns RunResult with error.
  Used by batch tests and ORCA loading path.

Pipeline is the convenience wrapper. CalculationRunner is the canonical
sequences for specific analysis patterns. Both respect dependency order.

### 2. WriteFeatures pattern

Each ConformationResult implements WriteFeatures() to serialize its own
contribution to NPY files. ConformationResult::WriteAllFeatures()
traverses all attached results plus writes identity arrays (positions,
elements, residues). This is the distributed feature output mechanism —
no centralized FeatureExtractionResult.

### 3. CovalentTopology as geometry-to-topology boundary

CovalentTopology::Resolve() takes atoms, rings, residues, and ONE set
of positions. Returns a self-contained topology object with bonds,
per-atom bond indices, and H-parent assignments. The positions are used
only for covalent radius checks, then discarded. The result is pure
topology that Protein delegates through.

### 4. GromacsEnsembleLoader authority separation

TPR is single authority for topology (atom names, residue names,
charges, elements). Pose PDBs provide positions only. NamingRegistry
translates CHARMM→canonical at the loading boundary. Charges returned
separately — caller constructs appropriate ChargeSource.

### 5. RingBondedExclusionFilter

The fourth concrete KernelEvaluationFilter, added during the calculator
audit (2026-04-02). Topological check: ring vertices + their bonded
neighbours are excluded from through-space kernel evaluation. Used by
BS, HM, PQ, RS. Fixed the HM fused ring artifact (1.127→1.000).
Distinct from DipolarNearFieldFilter (distance-based) — this is
topology-based.

### 6. CHARMM switching function in DispersionResult

Smooth C1-continuous taper (Brooks 1983) between R_switch=4.3A and
R_cut=5.0A. NOT a hard cutoff despite what APPLIED_MATHS_FIXES.md
says. Prevents feature discontinuities across MD ensembles.

### 7. Per-calculator shielding_contribution fields

Every calculator writes a `*_shielding_contribution` SphericalTensor
on ConformationAtom. These are what WriteFeatures outputs to NPY files.
They represent the geometric kernel decomposed into T0+T1+T2, before
any intensity/parameter multiplication. This is the primary output of
the calculator layer.
