# Comment-phase arc — src/ comment truthfulness (working tracker)

Transient working artifact for the `src/` comment-truthfulness phase of the
campaign in memory `project_doc_squareaway_2026-05-26`. **Delete when the phase
ends** and the single source of truth lands (`project_docs_single_source_of_truth`).
Not a permanent design doc.

## The bar — every comment, every group

CODE IS GROUND TRUTH. Verify each comment against the `.cpp` implementation —
not the header comment, not a doc.

1. **Truthful.** Every comment is true of the code. Fix a found lie in the same
   pass — never defer a lie.
2. **The comment serves THIS code, never the domain.** It earns its place only
   where the code needs help being understood. Never a supplemental textbook —
   no re-deriving physics/math, no domain tutorial. That belongs in the spec/docs.
3. **Pipeline / project narration goes away.** A comment that tours the
   program's global operation — naming other components and the order they run
   ("Seed does that on first frame", "that is Trajectory::Run's job",
   "Phase 6 …") — is project textbook; delete it (that architecture lives in
   OBJECT_MODEL / PATTERNS, not the source). A **precondition** is the opposite
   and survives: it constrains how THIS code is used ("call Seed() before
   reading the Protein", "caller owns the returned pointer", "positions must be
   PBC-whole"). Keep a precondition only if real and verified, phrased as a
   requirement on this code's use — never as a tour of the pipeline. This is the
   single sharpest filter; the trajectory groups (G10–G14) are thick with
   pipeline narration and lose most of it.
4. **Inline comments are the highest risk — default delete.** Keep one only for
   durable, non-obvious, code-local utility: a unit, a gotcha, a contract, a
   non-obvious *why*.
5. **Kill darlings.** Cut AI handwaving, motivational/marketing filler, narration
   that just restates adjacent code, stale "future"/changelog notes, and
   over-long blocks that are a maintenance/risk surface. Say less, not more;
   grounded plain language. When in doubt, cut.
6. **Never touch code.** Comment text only. If a comment can't be made true
   without a code change, delete the comment and log the code issue.

## Method per group

1. **Inventory.** Enumerate every file in the group (regenerated from `find src`
   at group start so nothing slips) as a checkbox, each with the docs to
   cross-check it against. The group ends only when all boxes are ticked.
2. **Lead (Claude)** edits comments to the bar — surgical, comment-only.
3. While editing, cross-check code+comment against the now-truthful doc set and
   **log every doc(`.md`)-vs-code mismatch** to my per-agent file
   `doc-discrepancies/<UTC-stamp>-<slug>-disc.md` (NEVER a shared file). One entry
   per line: `doc:line | doc claim | code truth | src file:line | DOC|CODE`. A
   later aggregation pass dedupes all `*-disc.md` into one doc-cleanup.
4. **Follow (codex, reasoning = xhigh)** examines the diff: accuracy vs source,
   clarity, unkilled darlings / broad-op / textbook, and a CODE-TOUCHED check
   (the diff must be comment-only).
5. Act on findings (re-verify vs source) → commit the group atomically → update
   the baton + tick boxes here.

Completeness check: the union of the group membership rules below must equal
`find src -maxdepth 2 \( -name '*.h' -o -name '*.cpp' \)` minus the exclusions.
Verify no file is unassigned before declaring the phase done.

## Excluded

- `src/generated/` (`LegacyAmberSemanticTables.{h,cpp}`) — generator output;
  hand-edits are futile. Skip.

## Doc abbreviations

OM = OBJECT_MODEL.md · CON = spec/CONSTITUTION.md · PAT = PATTERNS.md ·
ARCH = doc/ARCHITECTURE.md · CAT = GEOMETRIC_KERNEL_CATALOGUE.md ·
MG = spec/MATHS_GOALS.md · AM = spec/APPLIED_MATHEMATICS.md ·
PF = spec/PHYSICS_FOUNDATIONS.md · API = python/API.md + python/nmr_extract/_catalog.py.

## The arc  (core model → calculators outward → job-running → core code)

Status: [ ] todo · [~] in progress · [x] done (commit)

- [~] **G0 — verified stragglers** (never-defer; spans groups, done first):
  BuildResult.h, TrajectoryProtein.h, ConformationAtom.h, OrcaShieldingResult.h,
  MopacResult.h. Docs: OM, CON.
- [ ] **G1 — core object model & Result base.** Protein, ProteinConformation,
  ConformationAtom (structure/identity comments only — per-calculator field
  blocks are verified within each calculator's group), ConformationResult, Atom,
  AtomEvent, Residue, Bond, Ring, BuildResult, ProteinBuildContext,
  ProteinTopology, GeometryResult, GeometryChoice, CalculatorContract,
  CalculatorConfig, RecordBag, DenseBuffer, SelectionRecord, errors. Docs: OM, CON, PAT.
- [ ] **G2 — core types, enums, constants.** Types, SemanticEnums, AminoAcidType,
  PhysicalConstants, SolventEnvironment. Docs: OM, CON, PF.
- [ ] **G3 — topology & perception.** CovalentTopology, RingTopology,
  MolecularGraphResult, LegacyAmberTopology, LarsenResidue, NamingRegistry. Docs: OM, CON, PAT.
- [ ] **G4 — charge / force field / protonation / readers.** ChargeSource,
  AmberChargeResolver, AmberPreparedChargeSource, AmberLeapInput,
  ForceFieldChargeTable, ChargeAssignmentResult, ProtonationState,
  ProtonationDetectionResult, ReduceProtonation, PdbFileReader,
  GromacsToAmberReadbackBlock. Docs: CON, ARCH.
- [ ] **G5 — calculators: ring current.** BiotSavartResult, HaighMallionResult,
  McConnellResult, RingSusceptibilityResult. Docs: CAT, CON, MG, AM, PF, OM, API.
- [ ] **G6 — calculators: electrostatics.** CoulombResult, ApbsFieldResult,
  apbs_bridge, PiQuadrupoleResult, EeqResult, MopacCoulombResult. Docs: CAT, CON, MG, AM, PF, OM, API.
- [ ] **G7 — calculators: H-bond / tripeptide / Larsen.** HBondResult,
  LarsenHBondShieldingResult, LarsenHBondGrid, PlanarGeometryResult,
  TripeptideBackboneShieldingResult, TripeptideNeighborShieldingResult,
  TripeptideDftTable, TripeptidePoseAssembler. Docs: CAT, CON, OM, API.
- [ ] **G8 — calculators: QM / charge-derived / dispersion.** MopacResult,
  MopacMcConnellResult, AIMNet2Result, AIMNet2ChargeResponseGradientResult,
  EnrichmentResult, DispersionResult. Docs: CAT, CON, OM (MOPAC §2253), API.
- [ ] **G9 — calculators: solvent / structure / reference / delta.**
  WaterFieldResult, HydrationShellResult, HydrationGeometryResult, SasaResult,
  DsspResult, SpatialIndexResult, GromacsEnergyResult, BondedEnergyResult,
  OrcaShieldingResult, MutationDeltaResult, DemoResult. Docs: CAT, CON, OM, API.
- [ ] **G10 — trajectory core machinery.** Trajectory, TrajectoryProtein,
  TrajectoryAtom, TrajectoryResult, TrajectoryMoments, FullSystemReader,
  GromacsFrameHandler, GromacsFramePullResult, KernelEvaluationFilter,
  FrameNpyEmitter, FramePdbEmitter, RunConfiguration. Docs: OM, PAT §§13-18, ARCH.
- [ ] **G11 — trajectory results: ring-current + electrostatics.** `Bs*`, `Hm*`,
  `McConnell{Welford,ShieldingTimeSeries}`, `RingSusceptibility*TimeSeries`,
  `Apbs*`, `Eeq*`, `MopacCoulomb*TimeSeries`, `MopacMcConnell*TimeSeries`,
  `PiQuadrupole*TimeSeries`. Docs: OM, API, CAT.
- [ ] **G12 — trajectory results: H-bond + tripeptide.** `HBondCountWelford`,
  `HBondShieldingTimeSeries`, `LarsenHBond*`, `Tripeptide{Backbone,Neighbor}*`
  (TimeSeries/Welford/MethodTag/ResidualVec). Docs: OM, API, CAT.
- [ ] **G13 — trajectory results: QM/charge + solvent.** `AIMNet2*` (TimeSeries/
  Welford/Embedding), `Mopac{BondOrder,Charge}Welford`, `MopacVsFf14Sb*`,
  `Dispersion*TimeSeries`, `WaterField*`, `Hydration{Shell,Geometry}*`, `Sasa*`. Docs: OM, API.
- [ ] **G14 — trajectory results: structure / geometry / selection.** `Dihedral*`,
  `Dssp8*`, `Rmsd*`, `RingPucker*`, `RingNeighbourhoodTrajectoryStats`,
  `BondLengthStats`, `ChiRotamerSelection`, `PositionsTimeSeries`,
  `DftPoseCoordinator`, `JCouplingTimeSeries`, `GromacsEnergyTimeSeries`,
  `BondedEnergyTimeSeries`. Docs: OM, API, ARCH.
- [ ] **G15 — job-running / orchestration / CLI.** OperationRunner, Session,
  nmr_extract.cpp, Cli/* (CommonOptions, FrameEmission, ModeSpec, Parse,
  PrintUsage, Validate), OrcaRunLoader. Docs: ARCH, CLAUDE.md (5-mode spec), CON.
- [ ] **G16 — core infrastructure / IO.** RuntimeEnvironment, OperationLog,
  NpyWriter, TopologySidecar, CategoryInfoProjection. Docs: ARCH, OM, API.
