# Layer 0 Implementation Plan

**All source code goes in nmr-shielding/src/.**
**Every agent gets: CONSTITUTION.md, OBJECT_MODEL.md, PATTERNS.md, and their pass spec.**
**After every implementation pass, a correction pass reviews against the constitution and updates PATTERNS.md.**

## Progress (2026-04-03)

Pass 0 through Pass 2 (and their correction passes) are COMPLETE.
All 8 classical calculators implemented and batch-validated.
Fleet loader (GROMACS) operational. Feature output (NPY) working.
283 tests passing (+ 3 KaML skipped on ARM64).

**Done (Layer 0):** Core object model, PDB loading, GeometryResult,
DemoResult, OpenBabel bonds, DSSP, APBS (real PB), xTB (real binary),
spatial index, enrichment (16 atom roles), molecular graph (BFS),
protonation detection, ToolPaths, OperationLog, NamingRegistry,
IupacAtomIdentity, Protein::FinalizeConstruction pattern.

**Done (2026-04-01):** Protonation pipeline (ProtonationState,
PropkaProtonator, KamlProtonator, Protonator interface,
Henderson-Hasselbalch), typed ChargeSource hierarchy
(ParamFile/Prmtop/GmxTpr/Stub), ForceField enum, OrcaRunLoader
(prmtop+XYZ→Protein), OrcaShieldingResult (dia+para+total tensors
with element-verified ordering).

**Done (2026-04-02, early):** GEOMETRIC_KERNEL_CATALOGUE.md (mathematical
foundation for all 8 calculators, McConnell full tensor derivation).
MutationDeltaResult (rich WT-ALA comparison as ConformationResult on
WT conformation, with ring proximity, cylindrical coords, exponential
decay, distance decay curve). McConnellResult (first classical
calculator, full asymmetric tensor M with T0+T1+T2, batch-validated
on 465 protein pairs with 288,741 matched atoms, 0 failures).
APPLIED_MATHS_FIXES.md. Spec docs consolidated into repo under spec/.

**Done (2026-04-02, late):** CoulombResult (vacuum E-field V/A and EFG
V/A² from partial charges, decomposed backbone/sidechain/aromatic,
solvent = APBS − vacuum). RingSusceptibilityResult (full McConnell
tensor from ring centers, T0=f at 3e-15 across 914K pairs, per-ring-type
validated). HBondResult (full McConnell tensor from DSSP H-bond
partners). KernelEvaluationFilter framework (ABC + DipolarNearFieldFilter,
SelfSourceFilter, SequentialExclusionFilter, KernelFilterSet on all
calculators, per-filter rejection logging). T2 independence verified
across all 6 calculator pairs (mean |cos| 0.38-0.47, random = 0.36).
Unit fixes: APBS converted to V/A, vacuum fallback removed, ke and
kT/e as named constants. PATTERNS.md lesson 21 (full analytical
validation process).

**Done (2026-04-02, session 2):** BiotSavartResult (Johnson-Bovey
double-loop, wire segments in SI, G = -n⊗B×PPM rank-1 tensor).
HaighMallionResult (fan triangulation, 7-point Gaussian quadrature,
adaptive subdivision; TWO tensors: raw integral H symmetric traceless
A^-1, full kernel G = -n⊗(H·n) rank-1). Sign convention verified
analytically: sigma = I × G gives +1.40 ppm at 3A above PHE (Case
1995). Batch-validated on 465 pairs, 906K atom-ring pairs, 0 failures.
BS-HM T2 redundancy confirmed (cos = 0.999). T2 independent of all
4 prior calculators (0.38-0.62). DFT proximity 7.8x. TRP fused ring
superposition: both additive (BS 1.000, HM 1.000 after
RingBondedExclusionFilter; previously HM 1.127 was an artifact). Boyd &
Skrynnikov JACS 2002 tensor construction. 249 tests.

**Done (2026-04-02, session 3):** PiQuadrupoleResult (Stone T-tensor
EFG from point axial quadrupole, corrected formula from catalogue —
original violated Laplace, new is traceless at 5e-15 across 915K pairs).
DispersionResult (London 1/r^6 kernel summed over ring vertices, distance
filter [1.5, 5.0] A, 462K vertex contacts across 468 proteins).
Both batch-validated on 468 proteins, 0 failures. T2 independence:
PQ vs McConnell 0.38, PQ vs BS 0.57, Disp vs BS 0.71. Per-type
validated: 5-membered rings (TRP5, HIE) produce stronger PQ signal
than 6-membered (expected from 1/r^5 decay). 265 tests, 8/8 calculators.

**Done (2026-04-02, session 3, cont'd):** Prmtop recovery — regenerated
prmtop for 253 proteins using tLeap ff14SB. Training set expanded to
723 clean protein pairs. Calculator audit: RingBondedExclusionFilter
added to BS, HM, RS, PQ. All 8 calculators re-validated on 720-723
proteins. 266 tests.

**Done (2026-04-02 to 2026-04-03):** WriteFeatures on all 8 calculators
+ ORCA: ~29 NPY arrays per conformation. Pipeline.cpp
(RunClassicalCalculators with PipelineOptions for charges, xTB, APBS).
CovalentTopology extracted as geometry→topology boundary object.
Conformation type hierarchy: PredictionConformation, MDFrameConformation
added. Orthogonal access: Conformation() for primary, typed accessors
for metadata consumers. Loader honesty: OrcaRunLoader→AddPrediction,
PdbFileReader→AddConformation (no fake crystal metadata).

**Done (2026-04-03):** GromacsEnsembleLoader — reads GROMACS TPR via
libgromacs for topology authority, pose PDBs for positions.
PreloadedChargeSource wraps TPR charges. NamingRegistry CHARMM→canonical
translation. 10 MDFrameConformations per protein with walker, time_ps,
weight, rmsd, rg. Full pipeline per frame verified: 14 results × 10
frames = 290 NPY arrays. 283 tests pass, baselines byte-identical.

**Done (2026-04-03, late):** Builder + OperationRunner unification.
BuildResult as common return type from all 3 builders (PDB, ORCA,
GROMACS). OperationRunner replaces both Pipeline and CalculationRunner
as single home for all ordered sequences. reduce (Richardson lab)
linked as C++ library for PDB protonation — PDB-standard H naming,
HIS tautomer assignment. BuildFromPdb protonates bare PDB, assigns
ff14SB charges, computes net formal charge. xTB and APBS now
first-class citizens (run automatically when charges available).
nmr_extract CLI: --pdb, --orca, --mutant, --fleet. 263 tests pass.

**Next:** ParameterCorrectionResult (e3nn model integration).
Batch test migration to OperationRunner (4 files).
StubChargeSource removal. APBS C bridge global state fix.

---

## Pass 0: Working Protein + DemoResult (the Rosetta Stone)

**Status**: COMPLETE (2026-04-01)

**Goal**: Prove the architecture end-to-end with a real protein and
real data flowing through the ConformationResult framework. This pass
produces the example code that every subsequent agent copies.

**Scope**:

1. **Build system**: CMakeLists.txt with modular structure. Link Eigen3,
   cifpp, sphericart headers from extern/. GTest for testing. Must compile
   on Linux ARM64 (DGX Spark, GCC 13).

2. **Core types**: Vec3 (Eigen::Vector3d), Mat3 (Eigen::Matrix3d),
   SphericalTensor with WORKING Decompose() and Reconstruct() via
   sphericart. Unit tested: decompose a known Mat3, reconstruct it,
   verify roundtrip. T0, T1, T2 components verified against hand
   calculation.

3. **Naming and identity foundations**: Adapt from old code. These are
   the typed side of the string boundary. cifpp provides CCD data at
   construction time; these tables ARE the runtime authority.

   Adapt from biot-savart/src/Sequence/:
   - **AminoAcid enum** + three-letter code parsing + IsAromatic
   - **AminoAcidType table**: the single authority for amino acid
     chemistry. Every atom, ring, chi angle, protonation variant for
     all 20 amino acids. Ring detection reads from this.
   - **IupacAtomIdentity**: IUPAC ground truth for atom role, bonded
     partners, backbone/sidechain classification. The VALIDATOR.

   Adapt from biot-savart/src/Registry/:
   - **NamingRegistry**: the gatekeeper for residue and atom name
     translation between tool contexts (PDB, AMBER, CHARMM, ORCA,
     OpenMM, etc.). All translations happen here, once, at the tool
     boundary. After translation, typed objects carry properties.

   **The rule**: after PDB loading, objects answer questions about
   themselves. `atom.CovalentRadius()`, `residue.IsAromatic()`,
   `ring.NitrogenCount()`. Not library calls with strings as
   parameters. The information was resolved at construction time.
   See PATTERNS.md "Objects answer questions about themselves."

4. **PDB reader**: Adapt the existing cifpp-based reader from
   biot-savart/src/Protein/PdbFileReader.cpp to produce new typed
   objects. This is adaptation, not rewrite — the old reader handles
   alternate conformations, insertion codes, multi-model, hydrogen
   naming, GROMACS quirk filtering, CCD-based hydrogen parent
   assignment, and sentinel-based naming convention auto-detection.

   cifpp is a construction-time dependency. The reader calls cifpp,
   the NamingRegistry, and the AminoAcidType table to produce:
   - Protein with residues, atoms (Element enum, AtomRole from
     IupacAtomIdentity, not string matching)
   - Bond connectivity (via OpenBabel if available, or deferred)
   - One CrystalConformation with const positions and ConformationAtom
     per atom (position set, all computed fields default-initialised)
   - Ring detection from AminoAcidType ring definitions + atom
     presence (HIS/HID/HIE distinguished by which nitrogen carries H)
   - Protonation state detected from hydrogen atoms present

5. **Object model skeleton**: Protein, ProteinConformation (base +
   CrystalConformation at minimum), ProteinBuildContext, Residue, Atom,
   Bond, Ring type hierarchy (all 8 types with const virtual properties).
   ConformationAtom with ALL fields from OBJECT_MODEL (default-initialised).
   ConformationResult ABC. Template Result<T>(), HasResult<T>(),
   AttachResult<T>(), AllResults(). Conformation back-pointer to Protein.

6. **DemoResult**: A trivial ConformationResult that exercises the full
   pattern. Computes something real but simple — e.g., per-atom distance
   to nearest ring center (uses Eigen, ring geometry, stores a double
   and a Vec3 direction per atom). Declares a dependency on
   GeometryResult. Has a Compute() factory. Stores data in its own
   typed struct. Query methods that answer questions ("what's the
   nearest ring to this atom?").

7. **GeometryResult**: Ring geometry (center, normal from SVD, radius,
   vertex positions), bond geometry (midpoint, length, direction).
   Pre-built collections (rings_by_type, bonds_by_category). Ring pair
   properties. This is needed by DemoResult and by everything after.
   Uses Eigen for SVD. Uses sphericart for nothing yet — but
   SphericalTensor is used in the DemoResult to decompose a trivial
   test tensor, proving the pipeline works.

8. **Two conformations**: Load 1UBQ. Create CrystalConformation from
   the PDB coordinates. Create a second conformation (e.g., jitter the
   positions slightly, or load a second model if multi-model available,
   or use AddDerived). Attach GeometryResult and DemoResult to BOTH.
   Verify they are independent — changing one doesn't affect the other.

9. **Example calling code**: A main() or test that demonstrates the
   complete pattern:
   ```
   auto protein = LoadProtein("1UBQ.pdb");
   auto& conf1 = protein.CrystalConformation();
   conf1.AttachResult(GeometryResult::Compute(conf1));
   conf1.AttachResult(DemoResult::Compute(conf1));
   auto& demo = conf1.Result<DemoResult>();
   double d = demo.NearestRingDistance(42);
   ```
   THIS CODE IS THE PATTERN. Every subsequent agent reads it.

**Deliverables**:
- Compiling project with modular CMake
- SphericalTensor with tested Decompose/Reconstruct
- PDB reader producing typed Protein from 1UBQ.pdb
- Ring type hierarchy with 8 types, correct properties
- ConformationResult framework with template access
- GeometryResult: real ring/bond geometry
- DemoResult: full pattern end-to-end
- Two conformations with independent results
- Example code demonstrating the pattern
- Unit tests for all of the above

**Agent gets**: Constitution, Object Model, this pass spec, the old
PDB reader source for adaptation, 1UBQ.pdb, sphericart headers,
nanoflann headers

**What this pass does NOT include**: DSSP, charges, APBS, xTB,
spatial index (nanoflann), McConnell, Biot-Savart, features, ML.
Just the skeleton that everything hangs on, proven working.

---

## Pass 0R: Review + Correction

Check:
- Does it compile and all tests pass?
- Is the template Result<T>() mechanism correct and clean?
- Are ring types classes with virtual const properties, not enum lookups?
- Are positions const after construction?
- Does the conformation back-pointer work?
- Are the two conformations truly independent?
- Is SphericalTensor Decompose/Reconstruct consistent with sphericart?
- Is the PDB reader free of string-based atom identity in the NEW code?
  (Old reader may use strings internally during parsing — that's OK.
  The OUTPUT must be typed objects.)
- Is the example code clean enough to be the pattern for all agents?
- Create initial PATTERNS.md from the working code

---

## Pass 1: Build System + External Tool Integration

**Status**: not started

**Goal**: Every external tool callable, linked, and tested. The working
demo from Pass 0 stays working throughout.

**Scope**:

1. **OpenBabel integration**: Bond perception, hybridisation assignment.
   Adapt from old code. Test: detect bonds in 1UBQ, verify count and
   types. Feed into Bond objects with BondCategory enum.

2. **libdssp integration**: Link and call. Test: compute DSSP for 1UBQ,
   verify secondary structure assignments match known values.

3. **APBS bridge**: Adapt from old code. Test: run APBS on 1UBQ with
   test charges, get back E-field vectors. Verify known values.

4. **xTB invocation**: Call xtb binary, parse output. Test: run on a
   small fragment, get charges and HOMO-LUMO gap.

5. **nanoflann spatial index**: Build KD-tree from atom positions.
   Test: nearest-neighbour queries return correct atoms and distances
   for known geometry.

6. **PROPKA integration**: Call propka, parse pKa predictions. Test:
   run on 1UBQ, get predicted pKa for HIS residues.

7. **Modular CMake**: Each tool integration is its own target.
   Adding a new tool = new subdirectory + one add_subdirectory line.

**Deliverables**:
- All external tools callable and tested
- Bond detection via OpenBabel producing typed Bond objects
- nanoflann KD-tree building and querying
- Pass 0 demo still works
- Integration test: load 1UBQ, detect bonds, run DSSP, build
  spatial index, all in sequence

**Agent gets**: Constitution, Object Model, Patterns (from 0R),
Pass 0 code, old tool integration source for adaptation,
1UBQ.pdb, tool binaries/libraries

---

## Pass 1R: Review + Correction

Check:
- Do all tool invocations work on DGX Spark (ARM64)?
- Are tool results returned as typed data, not strings?
- Is each tool's output consistent with the old code's output for 1UBQ?
- Is the CMake structure modular?
- Update PATTERNS.md with tool integration patterns

---

## Pass 2: First Real ConformationResults

**Status**: not started

**Goal**: Replace DemoResult with real science. The first real
ConformationResult types that wrap the tools from Pass 1.

**Scope**:

1. **EnrichmentResult**: Atom role assignment (16 roles from bond
   connectivity + residue type), hybridisation from OpenBabel,
   categorical booleans. Pre-built atoms_by_role collection.

2. **DsspResult**: Wraps libdssp output. Per-residue SS, phi, psi,
   SASA, H-bond partners. Physics query methods.

3. **ChargeAssignmentResult**: Load ff14SB charges. Per-atom partial
   charge and VdW radius. Physics query methods.

4. **SpatialIndexResult**: Wraps nanoflann. Three KD-trees (atoms,
   ring centers, bond midpoints). 15A neighbour lists with stored
   distances and directions. Depends on GeometryResult.

5. **ApbsFieldResult**: Wraps APBS output. Per-atom solvated E-field
   (Vec3) and EFG (Mat3 + SphericalTensor). Depends on
   ChargeAssignmentResult.

6. **XtbChargeResult**: Wraps xTB output. Per-atom Mulliken charges,
   HOMO-LUMO gap, polarisability tensor. Per-bond Wiberg bond orders.

7. **Dependency checking proven**: Attach APBS before Charges → error.
   Attach DSSP (no deps) → works. The framework catches mistakes.

**Deliverables**:
- Six ConformationResult types, all following the DemoResult pattern
- All attach to 1UBQ CrystalConformation and pass tests
- Dependency checking works for real dependencies
- DemoResult can be retired or kept as documentation
- Integration test: full enrichment pipeline on 1UBQ

**Agent gets**: Constitution, Object Model, Patterns, Pass 0+1 code,
1UBQ.pdb

---

## Pass 2R: Review + Correction

Check:
- Are results named singletons? (Result<DsspResult>() not results[2])
- Does dependency checking work for real cases?
- Are APBS EFG tensors stored as BOTH Mat3 and SphericalTensor?
- Are xTB charges stored alongside ff14SB charges (both available)?
- Is UDP logging automatic on every store?
- Does the seven query patterns requirement hold?
- Update PATTERNS.md

---

## Pass 3: MD Trajectory Loading + ORCA Delta Loading

**Status**: not started

**Goal**: Reading conformations in bulk and reading DFT comparison data.

**Scope**:

1. **Trajectory reader**: Read GROMACS XTC/TRR or PDB multi-model,
   create MDFrameConformations via protein factory methods.
   Frame metadata: frame_index, time_picoseconds, walker_index,
   boltzmann_weight.

2. **ORCA delta reader**: Load WT and ALA ORCA output files, compute
   delta tensors (WT - ALA shielding), match atoms between structures
   (position-based, not name-based). Store as OrcaShieldingResult.

3. **Delta decomposition**: Decompose delta tensors via sphericart.

4. **Bulk enrichment**: Run the Pass 2 pipeline on multiple
   conformations. Prove the system handles 50+ frames.

**Deliverables**:
- MDFrameConformation creation from trajectory frames
- OrcaShieldingResult with delta tensors
- Atom matching between WT and ALA structures
- Unit tests verifying T0 and T2 match known reference values
- Bulk enrichment test: load trajectory, enrich 50 frames

**Agent gets**: Constitution, Object Model, Patterns, Pass 0-2 code,
test trajectory, test ORCA files

---

## Pass 3R: Review + Correction

Check:
- Are MD frames typed as MDFrameConformation with proper metadata?
- Does the protein's conformation list handle 500 frames?
- Is the ORCA delta decomposition using correct normalization?
- Is atom matching position-based, not name-based?
- Update PATTERNS.md

---

## Pass 4: ParameterCorrectionResult Training (Python)

**Status**: not started

**Goal**: Train the e3nn model that provides corrected parameters to
classical calculators.

**Scope**:

1. **Training data preparation**: Extract features from Pass 0-3
   results (APBS E-field, xTB charges, DSSP, DFT shielding) into
   tensors suitable for e3nn.

2. **Model architecture**: e3nn with full spherical tensor output
   (L=0+L=1+L=2). Additional L=0 output heads for 93 global
   parameter corrections.

3. **Training**: On ~553 protein WT-ALA mutant set. Validated on
   held-out proteins.

4. **Export**: TorchScript checkpoint for C++ inference.

**This model is NOT an NMR predictor.** It learns corrected values for
classical physics parameters from environment data + DFT.

**Deliverables**:
- Training script reading from Pass 0-3 outputs
- Trained checkpoint + TorchScript export
- Validation report: residuals, parameter correction values,
  comparison with literature defaults

**Agent gets**: Constitution, Calculator Parameter API, training data

---

## Pass 4R: Review

Check:
- Does the model capture environment effects?
- Are residuals largest near aromatic rings and charged residues?
- Are the 93 parameter corrections physically reasonable?
- Is the T2 prediction meaningful (not just noise)?

---

## Pass 5: ParameterCorrectionResult C++ Integration

**Status**: not started

**Goal**: Load the trained model as a ConformationResult type.
Implement the query API and residual tracking.

**Scope**:

1. **ParameterCorrectionResult**: Load TorchScript model via LibTorch.
   Per-atom tensor predictions, per-calculator parameter corrections,
   residual that updates as classical results attach.

2. **Two-path architecture**: Every subsequent calculator runs with
   default AND corrected parameters. Delta logged.

3. **SubtractCalculatorContribution**: Working residual update.
   Mock calculator test: attach a fake calculator contribution,
   verify residual decreases.

**Deliverables**:
- ParameterCorrectionResult as a ConformationResult
- Query API working: PredictedDelta(), ResidualT0(), corrections
- Mock calculator residual test
- Integration test: full pipeline through model inference

**Agent gets**: Constitution, Object Model, Calculator Parameter API,
Patterns, Pass 0-4 code, trained checkpoint

---

## Pass 5R: Review + Correction

Check:
- Can a mock calculator attach and see the residual decrease?
- Is the per-irrep breakdown correct (T0, T1, T2 independently)?
- Are parameter corrections global, not per-atom?
- Does the two-path architecture log deltas correctly?
- Update PATTERNS.md

---

## After Layer 0

The system can:
- Load a PDB into a typed Protein with typed conformations
- Enrich with roles, hybridisation, categoricals
- Run DSSP, assign charges, build spatial index
- Run APBS for solvated E-field and EFG tensors
- Compute xTB charges, HOMO-LUMO gaps, polarisabilities
- Load hundreds of MD frames
- Load and decompose ORCA WT/ALA deltas
- Query corrected parameters for all 8 classical calculators
- Track residuals per irrep level as classical results attach
- Copy a protein with new protonation and re-enrich

Physics agents (Layer 1+) get this as their starting point.
They focus on calculators and features, not plumbing.

The DemoResult pattern and PATTERNS.md tell them exactly how.
