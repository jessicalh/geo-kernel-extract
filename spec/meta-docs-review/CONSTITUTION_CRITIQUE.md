# Constitution Critique — Code vs Document

**Date**: 2026-04-03
**Purpose**: Line-by-line audit of CONSTITUTION.md against the actual
C++ code. For each section, one of:
- **TRUE**: The code does what the doc says.
- **STALE**: The doc describes something that was planned, partially
  built, or superseded by a different design.
- **ASPIRATIONAL**: The doc describes future work that is still wanted
  but not yet built.
- **WRONG**: The doc contradicts what the code actually does.

The constitution created the coherence that made 8 calculators and 250
tests work. The goal is to make it *true* again, not to weaken it.
Principles that guided good code stay. Claims about nonexistent code go.

---

## Hardware Context (line 17)

**TRUE.** 128GB, no premature optimisation, 15A neighbourhoods stored.
The code does exactly this (SpatialIndexResult stores all neighbours
within 15A on every atom).

---

## The Object Model

### ProteinBuildContext (line 26)

**MOSTLY TRUE.** ProteinBuildContext.h exists with pdb_source,
protonation_tool, pH, force_field, prmtop_path, tleap_script_path,
crystal_resolution, organism, deposition_date, stripped, assumptions.
Has Clone() method.

**STALE detail**: "A copy may keep its build context or take a new
one" — Protein is non-copyable. Build context cloning exists but the
copy-and-modify pattern it serves does not.

### Protonation State (line 49)

**WRONG**: "One protonation state per ProteinConformation." The code
puts protonation on Protein (it determines the atom list). Different
protonation = different Protein = different atom indices. The
PROTONATION_DESIGN_HISTORY.md documents this explicitly: "The
constitution text is stale."

**Recommendation**: Change to "One protonation state per Protein.
Different protonation = different Protein instance (rebuild, not copy)."

### Protein (line 56)

**TRUE.** Protein holds sequence, residues, atoms, rings, bonds, build
context. Does not hold geometry or computed properties. Owns
conformations via factory methods.

**STALE detail in API example (line 64)**:
- `protein.MDFrames()` — does not exist. Actual: `MDFrameAt(i)`,
  `MDFrameCount()`.
- `protein.Predictions()` — does not exist. Actual: `PredictionAt(i)`,
  `PredictionCount()`.
- `protein.Conformations()` — does not exist as a typed collection
  accessor. Actual: `ConformationAt(i)`, `ConformationCount()`.
- `protein.CrystalConformation()` — actual name is `CrystalConf()`.

### ProteinConformation hierarchy (line 78)

**PARTIALLY TRUE.** The base class design is exactly as described:
positions const, back-pointer to Protein, ConformationResult singletons,
template access.

**STALE hierarchy**: Constitution describes:
```
ProteinConformation
+-- ExperimentalConformation
|   +-- CrystalConformation
|   +-- NMRConformation
+-- ComputedConformation
|   +-- PredictionConformation
|   +-- MinimisedConformation
|   +-- MDFrameConformation
+-- DerivedConformation
```

Code has (flat, no intermediate classes):
```
ProteinConformation
+-- CrystalConformation
+-- PredictionConformation
+-- MDFrameConformation
+-- DerivedConformation
```

NMRConformation, MinimisedConformation, ExperimentalConformation,
ComputedConformation do not exist. The intermediate grouping classes
were never needed — the flat hierarchy works.

**STALE detail**: "confidence_per_residue" on PredictionConformation —
actual field is a single `double confidence_`.

**Copy policy** (line 136): GeometryOnly and Full copy — **NOT
IMPLEMENTED.** Protein is non-copyable (deleted copy/move). No copy
mechanism exists. See Protonation Design History.

### Residue (line 139)

**TRUE** for the core: type, variant, sequence position, chain.
**STALE detail**: "Backbone geometry (phi, psi, omega)" — phi/psi come
from DsspResult (per-conformation), not stored on Residue directly.
Omega and chi angles: Residue.h has chi_atom_indices but computation
is in DsspResult. This is a placement question, not a bug.

### Bond (line 149)

**TRUE.** Bond.h has order, category, atom_a, atom_b. Categories match
exactly. Bond types are typed, not string-identified.

**STALE detail**: "length, direction, midpoint" — these are per-
conformation (computed by GeometryResult, stored on ProteinConformation
as bond_lengths/bond_directions/bond_midpoints), not on Bond itself.
The Bond is topology; geometry is per-conformation. This is correct
design but the constitution implies they're on the Bond object.

### Atom (line 157)

**TRUE in spirit, STALE in detail.** The code splits atom data between:
- Atom.h: element, pdb_atom_name, residue_index, bond_indices,
  parent_atom_index (identity only, no position)
- ConformationAtom.h: all computed properties (position, charges,
  tensors, ring neighbourhoods, etc.)

The constitution describes this split correctly in the ConformationAtom
section (line 933) but the "Atom (in a ProteinConformation)" heading
conflates the two. A reader might think there's one Atom object with
everything on it.

### Ring (line 178)

**TRUE.** Ring.h matches exactly: class hierarchy, virtual properties,
vertex indices, parent residue, fused partner. All 8 types implemented.

---

## Units and Representations (line 210)

**TRUE.** All units match the code. Angstroms, radians, elementary
charges, Mat3 + SphericalTensor stored together.

---

## Sign Convention (line 227)

**TRUE.** G_ab = -n_b * B_a * PPM_FACTOR. Verified by test. Consistent
across all calculators.

---

## ConformationResult (line 242)

**TRUE** for the core mechanism: named singletons, template access,
dependency checking, permanent attachment.

**STALE details**:
- "Throws if not yet computed" (line 257) — actually aborts with
  fprintf(stderr, "FATAL:...") + std::abort(). Not an exception.
  This is deliberate (PATTERNS.md: "return codes, not exceptions").
- `conformation.ResultNames()` (line 319) — does not exist.
- `conformation.BiotSavart()` convenience wrappers (line 263) — do
  not exist. All access is via Result<T>() template.

### Known result types (line 349)

**STALE entries**:
- `ParameterCorrectionResult` — does NOT exist as a .cpp/.h file.
  ConformationAtom.h has placeholder fields (correction_predicted_delta
  etc.) but no result class.
- `FeatureExtractionResult` — does NOT exist. Feature output is
  distributed: each ConformationResult implements WriteFeatures().
  ConformationResult::WriteAllFeatures() traverses all attached results.
- `PredictionResult` — does NOT exist. ConformationAtom.h has
  placeholder fields (predicted_T0, predicted_T2, confidence, tier).

**All other result types are ACCURATE** — GeometryResult through
RingSusceptibilityResult all exist with the listed dependencies.

---

## Extraction Order (line 335)

**TRUE in principle**: dependency graph drives order. The Pipeline.cpp
implements this as a hardcoded-but-correct tier-0/tier-1 sequence, not
a dynamic resolver. The CalculationRunner.cpp provides alternative
orderings for specific use cases. Both respect dependencies.

**ASPIRATIONAL**: "The pipeline resolves the order" dynamically (line
340). In practice the order is written out explicitly in Pipeline.cpp
and CalculationRunner.cpp. The dependency checking at AttachResult time
catches violations but doesn't auto-resolve order.

---

## Layer 0 Foundation (line 384)

**MOSTLY TRUE.** PDB loading, bond detection, ring detection,
protonation detection, DSSP, charge assignment, APBS, ORCA delta
pipeline — all working and tested.

**STALE items in Layer 0 acceptance criteria (line 480)**:
- "PDB loads and produces a Protein with a CrystalConformation" —
  PdbFileReader now produces AddConformation() (base class), not
  CrystalConformation. This was a deliberate change: bare PDB lacks
  crystal metadata. SESSION_20260403_FLEET_LOADER.md documents this.
- "PROPKA runs and produces pKa predictions" — PropkaProtonator.cpp
  exists but is NOT connected to protein construction. It can produce
  pKa values but the builder that would USE them doesn't exist.
- "The copy-and-modify pattern works" — NOT IMPLEMENTED. Protein is
  non-copyable. See Protonation Design History.
- "MD frames load as MDFrameConformations" — TRUE now (via
  GromacsEnsembleLoader). Constitution says "OpenMM or AMBER" but
  actual implementation is GROMACS/libgromacs.

---

## Inviolable Rules (line 500)

**Rule 1** (object model): TRUE.
**Rule 2** (no string parsing in calculators): TRUE.
**Rule 3** (geometry in conformations): TRUE.
**Rule 4** (nanoflann for spatial queries): TRUE.
**Rule 5** (sphericart for decomposition): TRUE.
**Rule 6** (sign convention): TRUE.
**Rule 7** (features subclassed): **STALE.** No Feature class hierarchy
exists. Features are output via WriteFeatures() per ConformationResult.
The 189-feature manifest in EXTRACTION_ORDER.md is implemented as NPY
arrays written by each calculator's WriteFeatures method, not as
subclassed Feature objects.
**Rule 8** (ring types are classes): TRUE.
**Rule 9** (one-way, accumulative): TRUE.
**Rule 10** (precompute for next extractor): TRUE.
**Rule 11** (terrible categoricals allowed): TRUE.

---

## External Libraries (line 537)

**STALE entries**:
- "e3nn: equivariant neural network (Python)" — not integrated.
- "LibTorch: TorchScript inference (C++)" — not integrated.

**MISSING**:
- libgromacs (GROMACS 2026.0) — used by GromacsEnsembleLoader.

---

## Copy Semantics (line 551)

**ENTIRELY STALE.** The entire section (lines 551-635) describes a
copy-and-modify pattern that does not exist. Protein is non-copyable.
The protonation design history documents why: different protonation
means different atom lists, which means different flat indices, which
means a REBUILD not a copy.

The principle (same geometry, different protonation → different
physics → delta is the experiment) is still correct and still wanted.
The mechanism is wrong. Should describe the REBUILD pattern instead.

---

## Minimum Representation Contract (line 639)

**TRUE.** Every calculator produces the full representation described.
Mat3 + SphericalTensor for all tensors. Per-ring and per-bond
attribution. This section is one of the most valuable in the
constitution and the code matches it exactly.

---

## Anti-Simplification Rule (line 714)

**TRUE.** The code stores full tensors everywhere. No calculator
produces scalar-only output.

---

## Core Principle: If You Computed It, You Stored It (line 747)

**TRUE.** All computed quantities are stored on ConformationAtom or
RingNeighbourhood/BondNeighbourhood structures.

---

## Enforcement: The Framework Stores (line 846)

**WRONG.** The constitution describes:
```
conformation.StoreRingContribution(atomIndex, ringIndex, result);
// Writes to atom's ring neighbourhood AND emits UDP log
```

This API does not exist. Calculators write directly to ConformationAtom
fields:
```
auto& ca = conf.MutableAtomAt(ai);
ca.ring_neighbours.push_back(rn);
ca.per_type_G_T0_sum[ti] += rn.G_spherical.T0;
```

The typed store interface with automatic logging and overwrite
detection was never built. What exists instead:
- OperationLog.h: a general-purpose UDP/stderr logger with channels.
  Calculators call OperationLog::Info() at macro level (start/finish
  of computation), not per-property-store level.
- The singleton guarantee (one writer per field) is enforced by
  convention and the dependency graph, not by a store framework.

The PRINCIPLE (one writer per field, dependency ordering prevents
overwrites) is upheld. The MECHANISM described is fiction.

---

## UDP Logging (line 971)

**PARTIALLY TRUE.** OperationLog.h exists with UDP support, channel
bitmask, JSON format, scoped logging. Calculators do log.

**WRONG**: The per-property-store logging described (lines 977-1010)
with automatic overwrite detection and read-before-write errors does
not exist. The actual logging is at the operation level ("BiotSavart
starting on protein X with N atoms"), not at the property level
("stored G_tensor on atom 42 from ring PHE-7").

---

## Seven Query Patterns (line 1027)

**PARTIALLY IMPLEMENTED:**
1. Nearest N rings by distance — YES via SpatialIndexResult
   RingsWithinRadius.
2. All rings filtered by type — YES via ProteinConformation
   rings_by_type map.
3. Ring pairs by mutual geometry — YES via ring_pairs vector.
4. Atoms filtered by role near a point — PARTIAL. KD-tree gives atoms
   near point. Role filtering is manual intersection (no combined
   query method).
5. Bonds filtered by category near a point — PARTIAL. Same pattern.
6. Graph neighbours by bond count — YES via MolecularGraphResult BFS.
7. Atoms filtered by element within radius — PARTIAL. Same pattern
   as #4.

**STALE detail**: "Query result enrichment" (line 1074) with automatic
direction/distance/cylindrical-coordinate attachment — this happens
inside each calculator's Compute() method, not as a general query
enrichment framework. The calculators compute cylindrical coordinates
themselves and store them on RingNeighbourhood. Correct behaviour,
different mechanism.

---

## Atom Classification / Role (line 1100)

**TRUE.** AtomRole enum, set at enrichment, const thereafter. All 17
roles exist. EnrichmentResult assigns them using typed properties only.

---

## Bond Category Taxonomy (line 1157)

**TRUE.** All 8 categories exist and match.

---

## Ring Classification (line 1179)

**TRUE.** All 8 ring types as classes with virtual properties. Matches
code exactly.

---

## Summary: What Needs Changing

### Sections that are TRUE and load-bearing (DO NOT TOUCH):
- Units and Representations
- Sign Convention
- Minimum Representation Contract + Anti-Simplification Rule
- "If You Computed It, You Stored It"
- Ring Classification
- Bond Category Taxonomy
- Atom Role Classification
- Inviolable Rules 1-6, 8-11
- ConformationResult core mechanism (singletons, template access,
  dependencies)

### Sections that need factual corrections:
- Protonation State: "per ProteinConformation" → "per Protein"
- ProteinConformation hierarchy: remove intermediate classes, add note
  that flat hierarchy was deliberate
- API examples: fix accessor names to match actual code
- Known result types: mark ParameterCorrectionResult,
  FeatureExtractionResult, PredictionResult as "planned, not yet
  implemented"
- Copy semantics: rewrite to describe rebuild pattern
- "Framework stores" section: rewrite to describe actual mechanism
  (direct writes, singleton guarantee by convention)
- UDP logging: describe actual operation-level logging, not
  per-property-store fiction
- Layer 0 acceptance: update to reflect actual state (PDB → base
  conformation, GROMACS not OpenMM, copy-and-modify not implemented)
- External libraries: add libgromacs, mark e3nn/LibTorch as future
- Inviolable Rule 7: rewrite to describe WriteFeatures pattern
- Query patterns: describe what's implemented vs aspirational

### Sections that should be REMOVED or moved to history:
- None. Every section has a principle worth preserving. The stale parts
  need correction, not deletion.
