# Constitution: NMR Shielding Prediction System

This document is the supreme constraint on all implementation work.
Every agent receives this at the top of every instruction. Violation
of any principle here means the work is WRONG and must be redone.
The agent must know this. To violate it is anathema.

To create alternate object definitions, dive into string processing,
or bypass the typed object model is forbidden. But the GOAL that
motivated such thinking should not be abandoned -- the objects should
simply be used as intended, extended if necessary.

---

## Hardware Context

128 GB RAM minimum. Avoid premature optimisation. If storing precomputed
distances, properties, tensors on every atom within 15A makes the code
simpler and the extractors faster, do it. Memory is not the constraint.
Correctness and accessibility of information to extractors IS.

---

## The Object Model (sacrosanct)

### ProteinBuildContext

How this protein instance was built. PDB source, protonation tool +
version + pH, force field, what was stripped, what was assumed, crystal
resolution, organism, deposition date. Friend to the Protein. Travels
with the protein instance. A copy may keep its build context or take
a new one.

Tool configuration (where binaries live, default parameters) is global,
not per-protein. ProteinBuildContext records what happened, not where
the tools are.

**Why this contract exists (for agents creating new build contexts):**
Every protein instance needs a build context because it is the mutable
part of a protein definition. When we copy a protein to test different
protonation, the build context tells us what changed. When we compare
WT and mutant, it tells us which is which. When we serialise results,
it is the provenance record. Future work you can't see will need to
rely on this. Provide the basics (source, pH, force field) so
downstream consumers don't invent their own tracking.

### Protonation State

The protonation state for a ProteinConformation of a given protein.
Determines which titratable groups are charged, which histidine
tautomer exists, which cysteines are bonded. Derived from build
context + structure. One protonation state per ProteinConformation.

### Protein

Sequence and science data. Amino acid types, chain addressing (chain
is an addressing scheme, not an entity), canonical atom templates.
Does NOT hold geometry. Does NOT hold computed properties. The protein
is what the molecule IS, independent of where its atoms are.

The protein owns its conformation list. Conformations are created
through typed factory methods on the protein:

```
protein.AddCrystalConformation(positions, metadata)
protein.AddMDFrame(positions, metadata)
protein.AddPrediction(positions, metadata)
protein.AddDerived(parent, description)

protein.CrystalConformation()    // exactly one or throws
protein.MDFrames()               // typed collection
protein.Predictions()            // typed collection
protein.Conformations()          // all, heterogeneous
```

No agent creates a ProteinConformation directly. No loose conformations.
**Exception:** trajectory streaming (GromacsFrameHandler) creates
free-standing conformations via the public constructor. These point at
the Protein for topology but are never added to its conformations\_
vector. They live for one frame and die. See spec/ENSEMBLE\_MODEL.md.

### ProteinConformation (typed hierarchy)

A conformation without its protein is meaningless. One geometric
instance. Outside of a ProteinConformation, no geometric computation
is allowed. Holds exactly one set of atomic positions.

**Positions are const after construction.** New geometry = new
ProteinConformation. Spatial index is never stale.

ALL computed properties live on the ProteinConformation as typed
ConformationResult objects, or on its atoms and rings (stored there
by those result objects). This is where ALL work happens. Every
extractor operates within a ProteinConformation. Systems of spatial
query and property access accumulate over extraction passes within
the ProteinConformation.

**Class hierarchy, self-type-reporting:**

```
ProteinConformation (base -- NOT abstract, fully functional)
+-- ExperimentalConformation
|   +-- CrystalConformation
|   |     resolution_angstroms, r_factor, temperature_kelvin, pdb_id
|   +-- NMRConformation
|         ensemble_member, restraint_count
+-- ComputedConformation
|   +-- PredictionConformation
|   |     tool (AlphaFold/OpenFold3/ESMFold), confidence_per_residue
|   +-- MinimisedConformation
|   |     parent, method, force_field, energy, converged
|   +-- MDFrameConformation
|         frame_index, time_picoseconds, walker_index, boltzmann_weight
+-- DerivedConformation
      parent, derivation_description, properties_invalidated
```

The base class does the heavy lifting: holds positions (const),
provides spatial queries, accumulates ConformationResult objects,
supports the full extraction pipeline. Subclasses add source-specific
metadata.

A ProteinConformation always holds a valid pointer to its owning
Protein. This is not optional. Conformations are never orphaned --
the Protein's lifetime always exceeds its conformations' lifetimes.
Do not add reference counting, weak pointers, or null checks for
this pointer. The design guarantees it is valid. Do not defend
against scenarios that cannot happen in this system.

**Why this contract exists (for agents creating new conformation types):**
Future agents will attach ConformationResults you don't know about yet.
They need spatial queries, atom access, and result attachment to work
identically on your conformation type as on every other type. Your
metadata is what makes your type WORTH having — it distinguishes your
conformation from others in the protein's list. The base class machinery
is what makes it USABLE by everyone who comes after you. You are building
for people you haven't met.

Copy policy:
- **GeometryOnly**: copy positions, clear all results
- **Full**: copy everything (only valid if build context unchanged)

### Residue

All properties of the amino acid at a sequence position. All possible
angles and chemical properties as an amino acid. Backbone geometry
(phi, psi, omega), sidechain geometry (chi angles), secondary structure
assignment (from DSSP). Type, protonation variant, sequence position,
chain address (not chain entity).

### Bond

All bonds, distinguishing those which change with protonation state
and those which do not. All chemical and mathematical properties:
bond order, bond category (PeptideCO, PeptideCN, BackboneOther,
SidechainCO, Aromatic, Disulfide, SidechainOther, Unknown),
start atom, end atom, length, direction, midpoint. Bond types are
typed, not string-identified.

### Atom (in a ProteinConformation)

Each atom in a ProteinConformation knows:
- Element, role (backbone, sidechain, hydrogen, heavy), partial charge,
  van der Waals radius, hybridisation (from OpenBabel at enrichment)
- Distance and direction to all other atoms within 15 Angstroms,
  with ACTUAL distances, stored and accessible
- Assigned bonds (typed)
- Every accumulated property from prior extraction passes:
  field values (as Mat3 AND as sphericart-decomposed components),
  B-field vectors, E-field vectors, scalar contributions per calculator
- NMR-specific categorical properties set at enrichment:
  is_backbone, is_amide_H, is_alpha_H, is_methyl, is_aromatic_residue,
  parent_hybridisation, H-bond donor/acceptor role
  (from the categorical analysis: at minimum the 19 classification
  schemes from the categorical analysis (now in learn/bones/), expressed as
  typed properties, not strings)
- If a future extractor needs a property not yet present, it is added
  to the atom's property set properly, with name AND description from
  the tool or calculation that produced it

### Ring

Typed by class hierarchy. A PheBenzeneRing, a HisImidazoleRing,
a TrpPyrroleRing are different types with different properties
(intensity, lobe offset, nitrogen count, aromaticity, size).
Ring types are classes, not enum values looked up in a table.

A ring has: vertices (atom indices), geometry (center, normal, radius),
type properties, parent residue.

**Why this contract exists (for agents creating new ring types):**
The ring current calculators evaluate every ring the same way: vertices,
center, normal, intensity. Your ring type's specific physics (how much
current, what lobe offset, whether nitrogen is present) is what makes
the Biot-Savart result DIFFERENT for your type vs others. If you don't
provide it through the typed interface, someone will hardcode it as a
string comparison. Later, when we refine intensities from mutant data
or add a new ring type, the typed interface means the refinement
propagates everywhere automatically.

Rings accumulate properties over extraction passes. After each full
pass over the ring's atoms, a ring-level property update runs.
If a geometric or field property can be precomputed for the ring
as a whole (e.g., ring-frame B-field components, the geometric kernel
summed across ring vertices), it is computed and stored on the ring
at the end of that pass, ready for the next extractor.

It is NOT forbidden to traverse the atoms of a specific ring to
harvest properties for calculation. In fact, it must be encouraged.
Ring-specific atomic traversal is how ring properties are built.

---

## Units and Representations (explicit contract)

- All distances: double, in Angstroms
- All directions: Vec3 (Eigen::Vector3d), normalised where documented
- All positions: Vec3, in Angstroms, in the ProteinConformation's coordinate frame
- All field tensors: Mat3 (Eigen::Matrix3d), in Angstrom-based physical units
- All spherical decompositions: sphericart components in canonical normalization
  (set once when a tensor is produced, stored alongside the Mat3)
- All angles: double, in radians
- All charges: double, in elementary charge units (e)

When an extractor receives a distance, it is Angstroms. When it
receives a tensor, it is Mat3 AND sphericart components. No conversion
at point of use. No ambiguity about units. No undocumented scaling.

---

## Sign Convention (chosen once, used everywhere)

NMR convention: shielding is positive. sigma > 0 = diamagnetic shielding.
Ring current above an aromatic ring produces positive sigma (the atom
is more shielded). In the ring plane, sigma is negative (deshielded).

The Biot-Savart geometric kernel G_ab follows this convention:
G_ab = -n_b * B_a * PPM_FACTOR (negative, matching NMR literature).

All calculators, all features, all models use this convention.
Verified by test: a proton 3 Angstroms above a PHE ring on the
normal axis has sigma > 0.

---

## ConformationResult: Named Singletons with Dependency Graph

Computed results are NOT fields on a mega-struct. They are separate
typed objects that attach to a ProteinConformation and never detach.

Each result type is a **named singleton** on the ProteinConformation.
Accessed by template, not by iteration:

```
auto& dssp = conformation.Result<DsspResult>();
auto& apbs = conformation.Result<ApbsFieldResult>();
auto& bs = conformation.Result<BiotSavartResult>();
auto& coulomb = conformation.Result<CoulombResult>();
```

The template accessor returns the real typed object. Throws if not
yet computed. The code reads like English:
  conformation.Result<BiotSavartResult>().TensorAtAtom(42, ring3)

Optional named convenience wrappers exist as sugar:
```
auto& bs = conformation.BiotSavart();       // calls Result<BiotSavartResult>()
```
These are optional and not the mechanism. See "Result Access: Template
Mechanism" below.

Each result type declares its dependencies:

```
class McConnellResult : ConformationResult {
    requires: {SpatialIndexResult}
    ...
};
```

At attach time, dependencies are checked. Missing dependency =
immediate logged error, not a runtime surprise.

**Once attached, permanent.** A ProteinConformation that has a
DsspResult keeps it for its lifetime. Results accumulate. Nothing
is removed.

When a result is attached, it stores properties on the atoms and
rings it computed for. The BiotSavartResult stores per-atom per-ring
tensors. The ring summaries update as part of attachment. The next
result type finds the data already on the atoms.

**Results ARE the calculation.** Their methods ARE the physics queries:

```
bs.SumT0ByRingType(atomIdx, PheBenzene)
hm.NearestRingContribution(atomIdx)
coulomb.BackboneField(atomIdx)
coulomb.SidechainField(atomIdx)
```

Not getters returning raw data. Methods answering physics questions.

### Result Access: Template Mechanism

Result access uses C++ templates for type safety and extensibility.
Adding a new ConformationResult type does NOT require modifying
ProteinConformation. The agent writes the result class; access works
automatically.

    // Type-safe access (throws if not attached)
    auto& bs = conformation.Result<BiotSavartResult>();

    // Safe existence check
    if (conformation.HasResult<MopacResult>()) { ... }

    // Iteration over all attached results
    for (auto& [type_id, result] : conformation.AllResults()) {
        std::cout << result->Name() << std::endl;
    }

    // What's been computed?
    auto names = conformation.ResultNames();

Internal storage: unordered_map<type_index, unique_ptr<ConformationResult>>.
The template methods perform the cast. Named convenience wrappers
(conformation.BiotSavart()) are optional sugar on top, not the
mechanism.

**Why templates, not named accessors as the mechanism:**
Named accessors require modifying ProteinConformation every time a new
result type is added. In agent-driven development, this is a reliable
source of errors: the agent writes a perfect result class but forgets
to add the accessor, or misspells it, or puts it in the wrong file.
Templates make the pattern self-contained within the result class.

---

## Extraction Order: Dependency-Graph-Driven

The manual phase ordering is replaced by the dependency graph. The
pipeline resolves the order:

1. Check which results are attached
2. Find results whose dependencies are all satisfied
3. Compute and attach them
4. Repeat until all desired results are attached

The order emerges from the science (can't run Coulomb before charges)
not from a numbered list. An agent adding a new result type declares
dependencies and the pipeline slots it in.

### Known result types (with dependencies)

```
GeometryResult                requires: nothing
DsspResult                    requires: nothing
ChargeAssignmentResult        requires: nothing
EnrichmentResult              requires: nothing
OrcaShieldingResult           requires: nothing (loaded from files)
SpatialIndexResult            requires: GeometryResult
ApbsFieldResult               requires: ChargeAssignmentResult
MopacResult                   requires: nothing (external tool, runs early as
                                conformation electronic structure precondition)
MolecularGraphResult          requires: SpatialIndexResult
BiotSavartResult              requires: SpatialIndexResult, GeometryResult
HaighMallionResult            requires: SpatialIndexResult, GeometryResult
McConnellResult               requires: SpatialIndexResult, GeometryResult
CoulombResult                 requires: ChargeAssignmentResult, SpatialIndexResult
MopacCoulombResult            requires: MopacResult, SpatialIndexResult
MopacMcConnellResult          requires: MopacResult, SpatialIndexResult, GeometryResult
HBondResult                   requires: DsspResult, SpatialIndexResult
DispersionResult              requires: SpatialIndexResult, GeometryResult
PiQuadrupoleResult            requires: SpatialIndexResult, GeometryResult
RingSusceptibilityResult      requires: SpatialIndexResult, GeometryResult

(Feature extraction is distributed: each ConformationResult implements
WriteFeatures(). See ConformationResult.h.)
```

### What each result stores on atoms/rings

Each result type, when attached, populates properties on the atoms
and rings of the ProteinConformation. The next result type in the
resolved order finds these properties already present. See the
per-calculator minimum output section for exactly what each result
stores.

---

## Layer 0: Foundation

Everything in this layer must be working and tested BEFORE any feature
extraction agents run. This is the infrastructure that all result
types depend on. If Layer 0 is broken, nothing above it can be trusted.

### PDB loading and object construction

- Parse PDB/mmCIF via cifpp
- Build Protein (sequence, residue types, canonical atom templates)
- Build first ProteinConformation from PDB coordinates
  (CrystalConformation with resolution, R-factor, temperature)
- The Protein owns the ProteinConformation. Created via factory method.

### Bond detection and ring detection

- Bond detection via OpenBabel bond perception, producing typed Bond
  objects (bond order, bond category: PeptideCO, PeptideCN,
  BackboneOther, SidechainCO, Aromatic, Disulfide, SidechainOther)
- Ring detection from residue type + atom presence, producing typed
  Ring objects using the class hierarchy (PheBenzeneRing,
  HisImidazoleRing, HidImidazoleRing, HieImidazoleRing,
  TrpPyrroleRing, TrpBenzeneRing, PhenolRing, IndolePerimeter)
- Ring geometry: center, normal (from SVD), radius, vertex positions

### Protonation state detection

- Detect protonation state from hydrogen atoms present in the structure
- HIS/HID/HIE distinguished by which nitrogen carries hydrogen
- Titratable group states (ASP/GLU protonation, LYS neutralisation,
  CYS disulfide bonding) detected from atom presence
- The protonation state determines ring types (HIS variants) and
  charge assignments

### External tool integration

All external tool invocations must be working and tested:

- **DSSP**: secondary structure, phi/psi, SASA, H-bond partners.
  Invoked via libdssp. Results stored as DsspResult on the
  ProteinConformation.
- **Charge assignment**: tleap + ff14SB force field. Produces partial
  charges on all atoms. Results stored as ChargeAssignmentResult.
- **APBS**: Poisson-Boltzmann solvation. Requires charges first.
  Produces solvated E-field (Vec3) and EFG (Mat3 + sphericart) at
  each atom. Results stored as ApbsFieldResult.
- **OpenBabel**: bond perception, hybridisation assignment. Used during
  bond detection and enrichment.

### DFT mutation comparison pipeline

The alanine mutation delta-shielding pipeline must be working:

- Loading wild-type ORCA output files (shielding tensors per atom)
- Loading ALA mutant ORCA output files (shielding tensors per atom)
- Atom matching between WT and ALA structures (identifying
  corresponding atoms despite different residue types)
- Computing delta tensors: WT shielding minus ALA shielding per atom
- Results stored as OrcaShieldingResult on the ProteinConformation
- This is the ground truth that the ML model learns to predict

### Protonation tools

At least two protonation prediction tools, tested and integrated:

- **PROPKA**: pKa prediction from structure. Given a PDB, returns
  predicted pKa values for all titratable residues.
- **KaML** (or equivalent second tool): ML-based pKa prediction.
  Provides an independent estimate for comparison.

The **copy-and-modify pattern** for generating conformations with
different protonation states must be tested:

1. Take an existing Protein with its ProteinConformations
2. Copy it (proper copy constructor)
3. Apply a new ProtonationState to the copy (e.g., different pH
   from PROPKA predictions)
4. Re-detect ring types (HIS variant may change)
5. Re-run charge assignment (charges change with protonation)
6. Re-run enrichment and extraction on the copy
7. The delta between original and copy IS the experiment

This is how pH scanning works: same protein, same geometry,
different pH -> different protonation -> different charges ->
different Coulomb field -> different predictions.

### MD frame loading

- Reading trajectory frames (e.g., from OpenMM or AMBER) into typed
  MDFrameConformation objects via the protein's factory method
- Each frame carries: frame_index, time_picoseconds, walker_index,
  boltzmann_weight
- Positions are const on each MDFrameConformation
- Ensemble analysis: one protein, many MDFrameConformations, each
  independently enriched and extracted

### Layer 0 acceptance criteria

All of the above must have passing unit tests before any agent
proceeds to implement result types (Layer 2 and beyond). Specifically:

- PDB loads and produces a Protein with a CrystalConformation
- Bonds are detected and typed correctly
- Rings are detected, typed by class hierarchy, and have correct geometry
- Protonation state is detected correctly for HIS/HID/HIE
- DSSP runs and produces secondary structure assignments
- Charge assignment runs and produces partial charges
- APBS runs and produces solvated fields
- ORCA outputs load and delta tensors are computed correctly
- PROPKA runs and produces pKa predictions
- The copy-and-modify pattern works: copy protein, change protonation,
  verify ring types and charges update
- MD frames load as MDFrameConformations with correct metadata

---

## Inviolable Rules

1. **No extractor or calculator may violate or alter the object model.**
   If you need information, query the objects. If the object doesn't
   have the property you need, ADD it to the object properly.

2. **No string parsing of atom names, residue names, or element symbols
   in calculators, condition sets, or feature extractors.** Typed
   properties, set once at load/enrichment time.

3. **All geometry lives in ProteinConformations.** No geometric computation
   outside a ProteinConformation.

4. **All spatial queries through nanoflann.** No linear scans over atoms.

5. **All tensor decomposition through sphericart.** One normalization.
   One library. Both sides.

6. **One sign convention (documented above), used everywhere.**

7. **Features are subclassed.** Adding a feature = writing a class.

8. **Ring types are classes.** Type hierarchy captures the science.

9. **Extraction is one-way and accumulative.** Each result adds to
   what prior results stored. Nothing is overwritten or removed.

10. **If a property can be precomputed for the next extractor, it is.**
    After each pass over atoms or rings, run a property update.
    The goal: the next extractor finds what it needs already there.

11. **Terrible categorical ideas are allowed.** A condition set can be
    crude or sophisticated. Both must use typed object properties.
    The crude version gets us running. The robust version comes later.

---

## External Libraries

- **Eigen3**: Vec3, Mat3, SVD
- **nanoflann**: spatial KD-tree, neighbour queries
- **sphericart**: spherical harmonics, canonical normalization
- **cifpp**: PDB/mmCIF parsing, CCD
- **libdssp**: secondary structure, H-bonds
- **OpenBabel**: bond perception, hybridisation
- **APBS**: Poisson-Boltzmann solvation
- **e3nn**: equivariant neural network (Python-side calibration model,
  learn/c_equivariant/ — not linked into C++)

---

## Protein Instance Copy Semantics

A Protein is an instantiated object with both static properties
(sequence, amino acid types, canonical atom templates) and dynamic
properties that depend on the instance context (protonation state,
build context, charge assignments).

### The copy-and-modify pattern

It must be possible, tested, and routine to:

1. Take an existing Protein instance with its ProteinConformations
2. Copy it (proper copy constructor)
3. Apply a new ProtonationState to the copy
4. Apply a new ProteinBuildContext to the copy (or keep the original)
5. Copy ProteinConformations from the original to the new instance,
   preserving geometry and any properties that remain valid
6. Re-run enrichment and extraction on the copy with the new
   foundational properties

This is how mutant analysis works: same protein, same geometry,
different protonation -> different ring types -> different features
-> different predictions. The delta between the original and the
copy IS the experiment.

This is how pH scanning works: same protein, same geometry,
different pH -> different protonation -> different charges ->
different Coulomb field.

This is how ensemble analysis works: one protein, many
ProteinConformations, each with its own enrichment.

### Copy policies

**GeometryOnly**: copy positions, clear all ConformationResult objects.
Use when foundational properties (protonation, charges) have changed
and all computed results must be recomputed.

**Full**: copy everything including all attached ConformationResult
objects. Only valid if ProteinBuildContext is unchanged. Use for
creating working copies where results are still valid.

### Categories of validity on copy

When copying a protein + ProteinConformation with modified foundational
properties, some computed properties remain valid and some must
be recomputed:

**Always valid on GeometryOnly copy (geometry-only):**
- Atom positions (these ARE the ProteinConformation)
- Bond connectivity (topology doesn't change with protonation
  for most bonds -- but disulfides and protonation-dependent
  bonds DO change)
- Neighbour lists and distances (purely geometric)
- DSSP secondary structure (backbone geometry unchanged)

**Invalid on copy (depend on foundational properties):**
- Partial charges (depend on protonation + force field)
- APBS E-field (depends on charges)
- Ring type classification (depends on protonation for HIS)
- All calculator ConformationResult objects (depend on charges,
  ring types, etc.)
- All features (derived from results)
- All predictions (derived from features)

**Conditionally valid (depends on what changed):**
- Hybridisation (usually unchanged unless bonds change)
- H-bond classification (may change if protonation changes donors)
- Categorical atom properties (may change with protonation)

The copy constructor must accept parameters that specify WHAT
changed, so it can mark the right properties as invalid and
preserve the rest. This is not an optimisation -- it is correctness.
A ProteinConformation copied with new charges must NOT retain the
old APBS field.

### Unit tested

The copy-and-modify pattern must have explicit tests:
- Copy protein, change protonation, verify ring types update
- Copy protein, change charges, verify Coulomb features change
- Copy ProteinConformation from protein A to protein B, verify
  positions transfer but field results do not
- Verify that features on the copy are different from the original
  when foundational properties differ

---

## Minimum Representation Contract

Every calculator and enrichment step must produce AT LEAST the
representations specified below. An implementation that produces
a correct scalar but omits the tensor is WRONG. An implementation
that uses a geometric shortcut instead of the full linear algebra
is WRONG, even if the scalar result matches.

The model and downstream extractors need the FULL representation.
Physics models need tensors. Equivariant models need spherical
components. The viewer needs vectors and matrices. Collapsing to
scalars destroys information that cannot be recovered.

### Per-calculator minimum output (at each evaluated atom)

**Biot-Savart ring current (BiotSavartResult):**
- Mat3: full geometric kernel G_ab (rank-1 but stored as full matrix)
- sphericart: T0, T1[3], T2[5] decomposition of G_ab
- Vec3: magnetic field B at the point
- Per-ring attribution: which ring produced which contribution
  (not just the sum)
- Ring-frame cylindrical B-field components: B_n, B_rho, B_phi
  (for the nearest ring, stored on the atom)

**Haigh-Mallion surface integral (HaighMallionResult):**
- Mat3: full rank-2 tensor (NOT rank-1)
- sphericart: T0, T1[3], T2[5]
- Vec3: effective B-field (from the B-field integral, separate
  from the tensor integral)
- Per-ring attribution

**McConnell bond anisotropy (McConnellResult):**
- Mat3: dipolar tensor per bond
- sphericart: T0, T1[3], T2[5]
- Per-bond-category subtotals (PeptideCO, PeptideCN, BackboneOther,
  SidechainCO, Aromatic, Disulfide, SidechainOther) as separate
  Mat3 + sphericart each
- Scalar McConnell factor (3cos^2 theta - 1)/r^3 is derived FROM the
  tensor trace, not computed instead of the tensor

**Coulomb EFG (CoulombResult):**
- Vec3: electric field E at the point
- Mat3: electric field gradient tensor V_ab
- sphericart: T0, T1[3], T2[5] of V_ab
- Decomposed: backbone vs sidechain vs solvent contributions
  (separate Vec3 and Mat3 for each)

**Pi-quadrupole (PiQuadrupoleResult):**
- Mat3: quadrupole EFG tensor
- sphericart: T0, T1[3], T2[5]
- Per-ring attribution

**Ring susceptibility anisotropy (RingSusceptibilityResult):**
- Mat3: dipolar tensor from ring center
- sphericart: T0, T1[3], T2[5]
- Per-ring attribution

**London dispersion (DispersionResult):**
- Mat3: anisotropic dispersion tensor (per ring, summed over vertices)
- double: isotropic 1/r^6 contribution (NOT omitted)
- sphericart: T0, T1[3], T2[5]
- Scalar: contact count, sum of 1/r^6

**H-bond dipolar (HBondResult):**
- Mat3: dipolar tensor to H-bond partner(s)
- sphericart: T0, T1[3], T2[5]
- double: isotropic cos^2 theta / r^3 contribution
- Geometry: distance, D-H-A angle, donor/acceptor classification

**APBS solvated field (ApbsFieldResult):**
- Vec3: solvated E-field
- Mat3: solvated EFG (field gradient)
- sphericart: T0, T1[3], T2[5] of the EFG
- These are FROM APBS, not recomputed. The authoritative solvated answer.

### Anti-simplification rule

If an agent produces an implementation where a calculator returns
only scalar T0 values, or only the diagonal of a tensor, or uses
`3*cos^2 theta - 1` instead of computing the full dipolar tensor and
extracting the trace: that implementation is WRONG.

The full tensor is the minimum. Scalars are derived from tensors.
Tensors are decomposed via sphericart. Both forms are stored.
The next extractor, the model, and the viewer all need the full
representation. Optimising it away is destroying information.

If an agent proposes simplifying the representation for performance
or clarity, the answer is NO. 128 GB of RAM. Store everything.
Clarity comes from the object model, not from reducing what we compute.

---

## Agent Anti-Pattern Enforcement

Every agent instruction must include:

"Your output will be reviewed against the minimum representation
contract. If any calculator produces less than the specified
minimum (e.g., scalar instead of tensor, sum instead of per-source
decomposition), the work will be rejected and redone. If you
propose simplifying or omitting a representation, explain why
in a comment and leave the full representation in place. The
reviewer will decide whether the simplification is acceptable.
Do not decide for them."

---

## Core Principle: If You Computed It, You Stored It

No extractor may compute a derived quantity (distance, angle, tensor,
scalar, field value) and hold it only as a local variable. If the
quantity was computed in the context of an atom, it is stored on that
atom. If it was computed in the context of a ring, it is stored on
that ring. The next extractor finds it there.

If you need a quantity, you retrieve it from the object that owns it.
You do not recompute it. If it's not there, that's a dependency
ordering error (a missing ConformationResult dependency), not a
missing getter you work around.

---

## Atom Property Model

### Static core (const after construction, set at load time)
- Element (typed enum, not string)
- PDB atom name (string, but ONLY for display/serialisation -- never
  used for identity in extractors or calculators)
- Residue membership (index into protein's residue list)
- Bond list (indices into ProteinConformation's bond list)

### Enrichment properties (set once, append-only)
Added by ConformationResult objects during attachment. Once set, never
overwritten or removed within a ProteinConformation. Each property
records its source result type.

- Role: backbone, sidechain, hydrogen, heavy (from enrichment)
- Categorical: is_amide_H, is_alpha_H, is_methyl, is_hbond_donor,
  is_hbond_acceptor, is_aromatic_residue, etc. (from enrichment)
- Partial charge, VdW radius (from ChargeAssignmentResult)
- Hybridisation (from OpenBabel at enrichment)
- DSSP: secondary structure code, phi, psi, SASA (from DsspResult,
  via residue)
- Spatial neighbourhood (from SpatialIndexResult, via nanoflann):
  all atoms within 15A, each with stored distance and direction Vec3

### Ring neighbourhood (structured, per nearby ring)
For each aromatic ring within range of this atom, a structure:
- Ring reference (index + type)
- Distance to ring center (double, Angstroms)
- Direction to ring center (Vec3, normalised)
- Cylindrical coordinates in ring frame (rho, z, theta)
- Per-calculator contributions from this ring:
  - G tensor (Mat3) + sphericart decomposition (T0, T1[3], T2[5])
  - B-field (Vec3) + cylindrical components (B_n, B_rho, B_phi)
  Built by BiotSavartResult et al. Exported via WriteFeatures().

### Bond neighbourhood (structured, per nearby bond)
For each bond within range:
- Bond reference (index + category)
- Distance to bond midpoint
- McConnell geometry (angle theta, factor (3cos^2 theta - 1)/r^3)
- Per-calculator contribution:
  - Dipolar tensor (Mat3) + sphericart decomposition
  Built by McConnellResult. Exported via WriteFeatures().

### Field accumulation (from calculator ConformationResults)
- Per-calculator total tensor: Mat3 + sphericart decomposition
- Total field (sum of all calculators): Mat3 + sphericart
- APBS solvated E-field (Vec3) and EFG (Mat3 + sphericart)
- Coulomb vacuum E-field (Vec3) + decomposed (backbone/sidechain)
- H-bond geometry and tensor

### Feature export
All per-atom tensor data is exported via WriteFeatures() for the
calibration pipeline and upstream e3nn prediction model.

---

## Ring Property Model

### Static core (from ring detection)
- Type (class -- PheBenzeneRing, HisImidazoleRing, etc.)
- Vertex atom indices
- Parent residue index
- Fused partner index (for TRP)

### Geometry (from ProteinConformation positions)
- Center (Vec3)
- Normal (Vec3, from SVD)
- Radius (double)
- Vertex positions (vector of Vec3)

### Accumulated properties (updated after each atom pass)
After each ConformationResult attachment that traverses the ring's
atoms, a ring-level property update runs. Properties that can be
precomputed as ring totals are stored here:

- Total B-field at ring center (Vec3)
- Ring current intensity used (double, from ring type)
- Sum of G_T0 contributions to all nearby atoms (for diagnostics)
- Any ring-frame summary statistics needed by the next result type

---

## Enforcement: The Framework Stores, Not the Extractor

An extractor does not directly write `atom.someProperty = value`.
The ConformationResult, during attachment, stores properties on
atoms and rings through a typed store interface:

```
// The result computes and stores
Result result = myCalculation(atom, ring, conformation);

// The framework stores it on the atom and logs it
conformation.StoreRingContribution(atomIndex, ringIndex, result);
// This writes to atom's ring neighbourhood structure
// AND emits UDP log: "BiotSavartResult -> atom 42 x ring PHE-7"
```

The store is:
- Typed (you can't store a double where a Mat3 goes)
- Logged (every store emits to UDP log automatically)
- Append-only (writing to a slot that's already filled is an error,
  logged as WARNING -- it means dependency ordering is wrong or a
  result is being computed twice)

The retrieval is:
- Typed (you get back what was stored, not a void*)
- Checked (retrieving from an empty slot is an error, logged as
  ERROR -- it means a dependency wasn't met)

This means: the UDP log is a complete record of what was computed,
when, by whom, for which atom/ring. When something goes wrong,
the log tells you which result stored the wrong value, not just
that the final feature is wrong.

---

## Static vs Dynamic Object Properties

Every object in the model has two kinds of properties:

### Static: determined by what the object IS
- An atom's element is static. It doesn't change between results.
- A ring's type (PheBenzene, HisImidazole) is static after
  protonation state is applied.
- A residue's amino acid type is static.
- A bond's category (PeptideCO, Aromatic, etc.) is static after
  bond detection.

Static properties are set at construction or enrichment and are
const thereafter. They travel with copies. They are not
result-dependent.

### Dynamic: determined by computation within a ProteinConformation
- An atom's partial charge is dynamic (depends on force field + protonation)
- An atom's ring neighbourhood tensors are dynamic (computed by
  ConformationResult objects)
- A ring's accumulated field properties are dynamic

Dynamic properties are set by ConformationResult attachments and are
ProteinConformation-specific. They do NOT travel with GeometryOnly
copies. They travel with Full copies only if the ProteinBuildContext
is unchanged.

This distinction matters for the copy-and-modify pattern:
when you copy a protein with a new protonation state, static
properties transfer, dynamic properties are invalidated and
must be recomputed.

### The Placement Rule (for agents)

- **Protein**: identity and topology. Element, residue type, bond
  connectivity, ring membership, sequence. Never changes with geometry.
  Never modified by a ConformationResult.
- **ProteinBuildContext**: provenance. How this protein instance was
  built. Applies to all conformations of this protein.
- **ProteinConformation**: positions and ALL computed data. Everything
  that depends on geometry lives here, attached as ConformationResult
  objects. This is where work happens.
- **EnrichmentResult** properties (role, hybridisation, categoricals)
  depend on topology and element, not geometry. They are the same across
  conformations with the same protonation. They are computed per
  conformation (simple, correct) not shared (complex, premature).

If you are unsure where something goes: if it changes when you move
an atom, it goes on ProteinConformation. If it doesn't, it goes on
Protein. If it describes how the Protein was prepared, it goes on
ProteinBuildContext.

### ConformationAtom: the per-atom computed data store

Each ProteinConformation holds a vector<ConformationAtom> parallel to
the Protein's atom list. ConformationAtom has typed, named fields for
every per-atom result from every ConformationResult type. The fields
are declared upfront with default values (zero, Vec3::Zero(),
SphericalTensor{}, etc.). ConformationResult objects fill them in
during Attach().

ConformationAtom is the "wide table" — one row per atom, one column
per computed property. Feature extraction reads across all columns.
The result objects own the dependency logic and per-source detail
(per-ring, per-bond attributions). ConformationAtom owns the per-atom
totals and summaries.

**Why this design:** The field names and types ARE the spec. An agent
implementing BiotSavartResult can see exactly which fields to fill in
by reading ConformationAtom.h. A SphericalTensor field cannot be
filled with a scalar. The type system enforces full tensor output.

**Singleton guarantee:** Each ConformationResult type is a singleton
on a ProteinConformation. BiotSavartResult is computed once, attached
once, and never replaced. No other result type writes to its fields.
If you need different parameters, that's a different conformation
(copy-and-modify), not a second BiotSavartResult on the same
conformation. This is why field accumulation on ConformationAtom is
safe — each field has exactly one writer, guaranteed by the type
system.

**Adding new fields:** If a new ConformationResult type produces
per-atom data, add new fields to ConformationAtom. Do not arbitrarily
subdivide by physics categories — keep typed by the functional
analysis performed. BiotSavartResult fields, McConnellResult fields,
CoulombResult fields. Not ring_current_scalars, ring_current_tensors,
electrostatic_scalars. The result type is the organising principle.

---

## UDP Logging (mandatory, automatic)

Every property store operation emits a UDP log entry. This is not
optional. It is not per-extractor. The framework does it.

Log entry format:
```json
{
  "result_type": "BiotSavartResult",
  "target": "atom",
  "index": 42,
  "source": "ring PHE-7 (index 3)",
  "property": "G_tensor",
  "T0": -0.0234,
  "stored": true
}
```

If a property is stored twice (overwrite attempt):
```json
{
  "result_type": "McConnellResult",
  "target": "atom",
  "index": 42,
  "property": "G_tensor",
  "WARNING": "slot already filled by BiotSavartResult",
  "stored": false
}
```

If a property is read before being stored:
```json
{
  "result_type": "CoulombResult",
  "target": "atom",
  "index": 42,
  "property": "apbs_efield",
  "ERROR": "property not set -- ApbsFieldResult not attached?",
  "retrieved": null
}
```

This logging is the crash diagnosis infrastructure for the extraction
pipeline. When features are wrong, the log tells you exactly where
the wrong value entered the system.

---

## Core Principle: Extractors Ask Physics Questions

Extractors do not iterate raw atom arrays. They ask physics questions
of the ProteinConformation and get back typed, enriched objects with
everything already computed. The ProteinConformation speaks physics.

---

## Seven Query Patterns

These are the queries extractors actually need, derived from
the categorical analysis of existing extraction code. The
ProteinConformation must provide all seven efficiently.

### 1. Nearest N rings by distance from a point
Returns: ring objects sorted by distance, each carrying geometry
(center, normal, radius, vertices), type properties, direction
and distance FROM the query point (precomputed by the query,
not by the caller), and all accumulated field properties.

### 2. All rings filtered by type
Returns: prebuilt collection. "All PHE rings" is a list maintained
by the ProteinConformation, not a filter applied at query time.
Combinable with spatial queries: "PHE rings within 10A" intersects
the type list with the KD-tree result.

### 3. Ring pairs by mutual geometry
Returns: pairs of rings with their mutual angle (n1.n2),
cross product magnitude (|n1 x n2|), center-to-center distance,
and type pair. For stacking analysis, interference patterns,
and multi-ring effects.

### 4. Atoms filtered by role near a point
Returns: atoms matching a role filter (BackboneHeavy, AmideH,
MethylH, etc.) within radius R of a point. Each atom carries
all its properties. The role IS the filter, not a post-query check.

### 5. Bonds filtered by category near a point
Returns: bonds matching a category (PeptideCO, PeptideCN,
SidechainCO, Aromatic, Disulfide) within radius R. Each bond
carries its typed properties, start/end atoms, geometry.

### 6. Graph neighbours by bond count
Returns: atoms reachable within N bonds of a source atom,
with the bond count (BFS distance) and the path. For
through-bond features: e^{-lambda*n} decay, conjugation chains,
electronegativity sums along the path.

### 7. Atoms filtered by element within radius
Returns: atoms of a specific element(s) within radius R.
For Coulomb sums (all charged atoms), dispersion contacts
(all heavy atoms), and similar calculations.

---

## Query Result Enrichment

When a spatial query returns an object, the result includes the
spatial relationship to the query point. This is computed once
by the query, not recomputed by the extractor.

A ring returned by "nearest rings from point P" carries:
- Everything the ring already has (geometry, type, accumulated fields)
- PLUS: distance from P to ring center (double)
- PLUS: direction from P to ring center (Vec3, normalised)
- PLUS: cylindrical coordinates of P in ring frame (rho, z, theta)
- PLUS: angle between (P - center) and ring normal (double)

An atom returned by "atoms near point P" carries:
- Everything the atom already has (properties, fields)
- PLUS: distance from P (double)
- PLUS: direction from P (Vec3, normalised)

This means: an extractor that queries for nearby rings and then
needs the angle to the ring normal DOES NOT compute it. It's
already on the result. The "oh I can calculate that too" property
emerges from the query returning everything the spatial relationship
implies.

---

## Atom Classification: Role, Not Subclass

Atoms are NOT subclassed by element in the deep hierarchy sense.
Element is a property (enum). What matters for extractors is ROLE:
the atom's function in the protein's chemistry and its relevance
to NMR physics.

### Role taxonomy (set at enrichment, const thereafter)

**Heavy backbone:**
- BackboneN: peptide nitrogen
- BackboneCA: alpha carbon
- BackboneC: carbonyl carbon
- BackboneO: carbonyl oxygen

**Heavy sidechain:**
- SidechainC: sidechain carbon (further classified by hybridisation)
- SidechainN: sidechain nitrogen (charged or neutral, donor or not)
- SidechainO: sidechain oxygen (charged or neutral, acceptor or not)
- SidechainS: sulfur (cysteine, methionine)
- AromaticC: carbon in an aromatic ring
- AromaticN: nitrogen in an aromatic ring (HIS, TRP)

**Hydrogen (by environment, which determines NMR response):**
- AmideH: bonded to backbone peptide N
- AlphaH: bonded to CA
- MethylH: bonded to terminal CH3 carbon
- AromaticH: bonded to aromatic ring carbon
- HydroxylH: bonded to OH (SER, THR, TYR)
- OtherH: all other hydrogens

The role is determined at enrichment from:
- Element (from PDB)
- Bond connectivity (from bond detection)
- Residue type (from sequence)
- Position in residue (backbone vs sidechain, from atom template)

The role is typed, not a string. Extractors filter by role.
"Give me all AmideH atoms within 6A" is a first-class query.

### Why role, not element subclass

Element determines physics (nuclear charge, covalent radius,
electronegativity). Role determines NMR context (how this atom
responds to ring currents, what Buckingham sensitivity it has,
whether it's an H-bond donor). Both are needed. Element is a
property. Role is the classification extractors filter on.

An AromaticH on a PHE ring and an AmideH on the adjacent backbone
are both hydrogen. They have the same element. But they respond
to the ring current completely differently: the AromaticH is IN
the ring current, the AmideH is a PROBE of the ring current.
The role captures this.

---

## Bond Category Taxonomy

Bond types are typed, not string-identified. The category is set at
bond detection time and is const thereafter.

- **PeptideCO**: backbone C=O
- **PeptideCN**: backbone C-N
- **BackboneOther**: N-CA, CA-C, CA-CB
- **SidechainCO**: sidechain C=O (ASN, ASP, GLN, GLU)
- **Aromatic**: within aromatic rings
- **Disulfide**: S-S bonds
- **SidechainOther**: all other sidechain bonds
- **Unknown**: unclassified

The finer granularity is required because McConnell bond anisotropy
depends on bond order -- sidechain C=O has much stronger anisotropy
than sidechain C-C. Coarser groupings can be derived at query time.

These are the categories used for McConnellResult per-bond-category
subtotals and for the bond spatial query filter.

---

## Ring Classification: Type IS the Class

Ring types are classes with physics properties baked in:

```
Ring (base)
+-- SixMemberedRing
|   +-- BenzeneRing (PHE)        -- I=-12.0, offset=0.64, N=0, full aromatic
|   +-- PhenolRing (TYR)         -- I=-11.28, offset=0.64, N=0, full aromatic
|   +-- TrpBenzeneRing (TRP 6)   -- I=-12.48, offset=0.64, N=0, full aromatic
+-- FiveMemberedRing
|   +-- TrpPyrroleRing (TRP 5)   -- I=-6.72, offset=0.52, N=1, reduced aromatic
|   +-- HisImidazoleRing (HIS)   -- I=-5.16, offset=0.50, N=2, weak aromatic
|   +-- HidImidazoleRing (HID)   -- I=-5.16, offset=0.50, N=2, weak aromatic
|   +-- HieImidazoleRing (HIE)   -- I=-5.16, offset=0.50, N=2, weak aromatic
+-- FusedRing
    +-- IndolePerimeter (TRP 9)  -- I=-19.2, offset=0.60, N=1, full aromatic
```

Each class provides: intensity, lobe offset, nitrogen count,
aromaticity level, ring size, vertex count. These are const
properties of the type, not looked up in a table.

The ProteinConformation maintains per-type collections:
- allBenzeneRings, allPhenolRings, allImidazoleRings, etc.
These are prebuilt at ring detection time. Spatial queries can
intersect with them: "BenzeneRings within 10A of point P."

### Ring-to-ring relationships (precomputed)

For every pair of rings in the protein, precompute and store:
- Center-to-center distance
- Normal-normal dot product (n1.n2)
- Normal-normal cross product magnitude (|n1 x n2|)
- Whether they are fused partners (TRP 5+6)

This is O(R^2) where R is the number of rings, typically < 20.
Trivial memory. Enormous value for stacking analysis and
multi-ring interference features.

---

## Collection Patterns

The ProteinConformation maintains prebuilt collections for the
patterns extractors need. These are built once during enrichment
and available to all subsequent ConformationResult attachments.

### By type
- Rings by type class (see above)
- Atoms by role (AmideH list, AromaticC list, etc.)
- Bonds by category (PeptideCO list, Aromatic list, etc.)
- Residues by amino acid type (all PHE residues, all HIS, etc.)

### By spatial region (from nanoflann)
- Atoms within R of any query point (KD-tree, on demand)
- Rings within R of any query point (separate KD-tree on ring centers)
- Bonds within R of any query point (KD-tree on bond midpoints)

### Intersection
"PHE rings within 10A" = type collection intersected with spatial query.
The ProteinConformation provides this as a single query, not two
queries that the extractor intersects manually.

---

## Accumulation After Ring Traversal

After each ConformationResult that traverses ring atoms:

1. Each atom that was visited has new properties stored on it
2. The ring runs a post-pass update that reads from its member
   atoms and computes ring-level summaries
3. Ring-to-ring pair properties are updated if the new data
   affects them (e.g., if we now know the B-field at each ring
   center from the other rings, update the mutual field properties)

The post-pass update is part of the result attachment, not a separate
step. When BiotSavartResult finishes computing contributions for
all atoms of ring PHE-7, it immediately updates PHE-7's ring-level
properties before moving to the next ring.

---

## Extensibility

### Adding a new ProteinConformation type
1. Define subclass with type-specific metadata
2. Add factory method on Protein
3. Done

### Adding a new ConformationResult type
1. Define subclass of ConformationResult
2. Declare dependencies (other result types that must be attached first)
3. Implement computation (what it stores on atoms/rings)
4. Done -- access via conformation.Result<YourType>() works automatically

The base ProteinConformation provides all the machinery. New types
inherit it. The goal: a feature extraction expert grabs a fresh
notebook (new conformation type or result type) and all the pencils
and calculators are already there.

---

## Implementation Process: Incremental Agent Passes

We do NOT implement everything in one shot. Agents run in layers:

**Layer 0** (Foundation): External tools and infrastructure.
PDB loading, bond/ring detection, protonation detection, DSSP,
charge assignment, APBS, ORCA delta pipeline, protonation tools,
MD frame loading. ALL external tool invocations working and tested
before any feature agent touches anything. See the Layer 0 section
above for the complete specification.

**Layer 1**: Core object model. Protein, ProteinConformation,
ProteinBuildContext, Atom, Residue, Bond, Ring hierarchy. No
physics. Just the structural objects with typed properties,
spatial queries, and the ConformationResult attachment framework.
Unit tested.

**Layer 2**: First result types. SpatialIndexResult, DsspResult,
ChargeAssignmentResult. These are the foundation that everything
else depends on. Attach to ProteinConformations from Layer 1.
Unit tested.

**Layer 3**: Ring current calculators. BiotSavartResult,
HaighMallionResult, PiQuadrupoleResult, RingSusceptibilityResult,
DispersionResult. Each is a separate agent pass. Each is reviewed
against the constitution before the next runs. Each stores full
tensors + sphericart on atoms. Unit tested.

**Layer 4**: Bond and atom calculators. McConnellResult,
CoulombResult, HBondResult, ApbsFieldResult. Same pattern.

**Layer 5**: Feature extraction. Reads from all previous results.
Produces the 189-feature manifest.

**Layer 6**: ML prediction and heuristic classification.

After EACH layer:
- Review the agent's output against the constitution
- Check: did the agent use the object model or bypass it?
- Check: are tensors full (Mat3 + sphericart) or simplified?
- Check: are results stored on atoms/rings via the framework?
- Fix violations before proceeding to the next layer
- The object model may need corrections from reality the agent
  found -- that's OK, but the correction goes into the model, not
  into a hack in the agent's code

This way complexity is introduced incrementally. Each layer is
reviewable. Violations are caught before they compound. The object
model evolves from contact with reality but stays coherent.

---

## Design Goal: Library as Language

The object model should feel like a language for expressing
physics questions, not a burden of getters and setters.

Good: "for each AmideH near this PHE ring, what's the B-field
component along the ring normal?" -- expressed as a query returning
enriched results, one dot product per result.

Bad: "get atom list, filter by name == 'H', filter by residue
has backbone N in bond list, compute distance to ring center,
check against threshold, get ring normal from ring geometry,
compute B-field from stored tensor..." -- expressed as a manual
traversal through raw data structures.

The difference is whether the object model has already answered
the structural questions (role, spatial relationship, ring
geometry) so the extractor only asks the physics question.

---

## What This Document Does NOT Specify

- Class names, file layout, API signatures (design phase)
- Exact feature list (feature specification, informed by categorical analysis)
- Training hyperparameters, model architecture (ML phase)
- Viewer implementation (UI phase)

These follow FROM the object model and dependency-driven extraction order.

---

## Thesis Hypothesis (what this system is FOR)

This system is NOT an NMR chemical shift predictor competing with
empirical methods. It is an exploration of whether classical
electromagnetic calculations (Biot-Savart ring currents, Haigh-Mallion
surface integrals, McConnell bond anisotropy, Buckingham E-field, etc.)
can be refined by training equivariant tensor product models against
DFT shielding tensors, and whether the decomposed per-calculator
results serve as physics-grounded features for full tensor
understanding of shielding perturbations.

The hypothesis: classical physics provides the geometry and the
decomposition. DFT provides the ground truth. Equivariant ML learns
the corrections and the confidence. Mutant DFT isolates which
classical calculator captures which physical mechanism. The end
product is not a single prediction but a decomposed, attributed,
confidence-gated set of physics estimates that a downstream model
(training against RefDB experimental NMR) can use as input features.

## Foundational Data Sources

### ProCS15 Tripeptide Database (1.9 million DFT calculations)
Reference shielding values for every tripeptide context. This IS the
baseline: what the local electronic environment contributes to
shielding, independent of ring currents. Our ring current calculators
measure the PERTURBATION from this baseline. The tripeptide data
must be a queryable resource, accessible to any extractor that
needs "what would this atom's shielding be without ring current
effects?"

The existing code has TripeptideDatabase, TripeptideAssembler, and
TripeptideConformation. These must be preserved and integrated as
a ConformationResult type (TripeptideResult) or a protein-level
data source.

### ORCA r2SCAN DFT Calculations
WT and mutant shielding tensors. These are full r2SCAN calculations,
not quick estimates. They provide ground truth for:
- WT-ALA deltas (isolating ring current contribution)
- Charge-flip mutant deltas (isolating Buckingham contribution)
- Ring modification deltas (isolating intensity from geometry)
- Salt-bridge breaker deltas (isolating solvation contribution)

### MD Conformational Ensembles
Every protein at the end of the day will have 5 × 10ns PLUMED
replicas. These are real trajectories that the full extraction
pipeline must process. This is not future work. It is the actual
data pipeline. Layer 0 must handle loading hundreds of
MDFrameConformations per protein efficiently.

### RefDB Experimental NMR
The upstream model trains against this. We provide decomposed features.
We do not train against RefDB directly in this system.

## Mutant Comparison: First-Class Operation

Classical features must be extracted on BOTH WT and mutant structures.
The delta in classical calculations (BS_WT - BS_ALA) tells us what
our classical model PREDICTS the mutation should do. The delta in DFT
(DFT_WT - DFT_ALA) tells us what ACTUALLY happened. The residual
(DFT_delta - classical_delta) is what classical misses -- and that
is exactly what the ML should learn.

Current pipeline extracts features from WT only and trains the ML to
predict the full DFT delta. The ML has to learn both the physics that
classical gets right AND the corrections. If we give it the classical
delta as a feature, it only has to learn the corrections.

For each mutant pair:
- WT protein instance → full extraction (all ConformationResults)
- Mutant protein instance → full extraction (same pipeline)
- Delta features: WT_feature - mutant_feature, per atom, per calculator
- These delta features ARE training inputs alongside the raw WT features
- The model sees: "classical says X should change, DFT says Y changed,
  learn Y - X"

This applies to all mutant categories:
- WT-ALA: ring current isolation
- ASP→ASN: Buckingham isolation (classical Coulomb delta vs DFT delta)
- TYR→PHE: ring intensity (classical BS delta vs DFT delta)
- LYS→MET: solvation (classical APBS delta vs DFT delta)

The ConformationResult framework handles this naturally: two Protein
instances, each with their conformations and results, queried and
differenced at the feature level.

## MOPAC Conformation Electronic Structure

MOPAC PM7+MOZYME computes QM-derived electronic properties for every
conformation as a ConformationResult (MopacResult). These are quantum-
informed charges and bond orders that capture polarisation and electron
sharing that the force field cannot represent. MOPAC runs as a
conformation electronic structure precondition — early in the pipeline,
before the geometric kernel calculators that depend on its output.

MOPAC PM7+MOZYME runs in ~45 seconds on 889-atom proteins via
subprocess (`/home/jessica/micromamba/envs/mm/bin/mopac`).

Both ff14SB charges (from ChargeAssignmentResult) and MOPAC charges
are available. Feature extractors and downstream calculators can use
either or both. MopacCoulombResult and MopacMcConnellResult depend
on MopacResult for QM charges and bond orders respectively.

### MopacResult

A ConformationResult that stores PM7+MOZYME computed properties.
Requires: nothing (external tool invocation, runs early).
Named accessor: conformation.Result<MopacResult>()

Stores per atom:
- Mulliken charge (double, elementary charge units)
- s-orbital population (double)
- p-orbital population (double)
- Valency — sum of Wiberg bond orders (double)
- Bond neighbours with continuous Wiberg bond orders (vector)

Stores per bond:
- Wiberg bond order (double, 0.01–3.0, continuous)
- O(1) lookup by atom pair: BondOrder(atom_a, atom_b)
- Topology bridge: TopologyBondOrders() parallel to protein.Bonds()

Stores per molecule:
- Heat of formation (double, kcal/mol)
- Dipole moment (Vec3)

## Mutation Delta: ConformationResult on the WT

MutationDeltaResult is a ConformationResult that attaches to the WT
conformation. The mutant conformation is an extra input to Compute,
not a dependency. This follows the ConformationResult pattern: the
delta is a property of THIS conformation, learned by comparing it
to a mutant.

```
auto delta = MutationDeltaResult::Compute(wt_conf, mutant_conf);
wt_conf.AttachResult(std::move(delta));
auto& d = wt_conf.Result<MutationDeltaResult>();
```

### What it stores (internally, not on ConformationAtom)
- Atom correspondence: position-based matching with element filter
  and bijection enforcement (0.5A tolerance)
- Per matched atom: DFT shielding delta (Mat3 + SphericalTensor),
  APBS field delta, charge deltas, DSSP delta, graph delta
- Ring proximity: cylindrical coordinates of each matched atom
  relative to each removed ring (z, rho, theta, McConnell factor,
  exponential decay)
- Mutation sites: which residues changed, what the WT/mutant types
  were, which rings were removed
- Summary statistics by element, by distance bin, backbone/sidechain

### Precondition checks
- WT conformation must have OrcaShieldingResult (formal dependency)
- Mutant conformation must have OrcaShieldingResult (runtime check)
- Both proteins must have the same residue count

### Future mutation types

Currently only MutationDeltaResult exists (handles WT-ALA and any
other mutation type). Future mutation categories (charge-flip,
ring-modification, solvation) may be implemented as separate
ConformationResult types (e.g. ChargeMutationDeltaResult) if they
need different atom matching or different stored data. The singleton
guarantee (one result per type) naturally supports this: a single
conformation could have both a MutationDeltaResult (vs ALA) and a
ChargeMutationDeltaResult (vs ASN) attached simultaneously.

## Feature Bar: Irrep-Level Breakdown per Source

This is the vision. Before any classical calculator runs, each atom
in a conformation already has rich tensor data at every irrep level.
Classical calculators must IMPROVE on what this baseline provides.
If a calculator doesn't improve the prediction beyond what's already
captured, it isn't adding physics.

### Per-atom feature sources and their irrep levels

**ORCA r2SCAN (non-mutant protein baseline DFT):**
- L=0: isotropic shielding (trace/3)
- L=1: antisymmetric part (T1, 3 components)
- L=2: traceless symmetric anisotropy (T2, 5 components)
- Diamagnetic and paramagnetic contributions separately
- Full Mat3 + sphericart decomposition
- This is a FEATURE INPUT, not a training target
- RefDB proteins may use r2SCAN-3c (faster, different basis set)
- Method recorded in ProteinBuildContext; model sees which level

**APBS solvated field:**
- L=0: E-field magnitude, EFG trace
- L=1 (1e): E-field vector (3 components)
- L=2 (2e): EFG tensor (5 components via sphericart)

**MOPAC PM7 (per atom):**
- L=0: Mulliken charge, s/p orbital populations, Wiberg bond orders

**DSSP (per residue, projected to atoms):**
- L=0: secondary structure one-hot, phi, psi, SASA

**Mutant DFT comparison (WT minus mutant, per atom):**
- L=0: delta isotropic
- L=1: delta antisymmetric
- L=2: delta anisotropy
- Full tensor comparison, same DFT level on both sides
- Existing mutants: r2SCAN. These are ground truth. Not rerun.

### Environment tools set the floor

The external tools (APBS, MOPAC, DSSP, ORCA DFT) provide tensor data
at multiple irrep levels. The calibration pipeline trains on this data.
Its prediction already incorporates everything these tools provide,
including their angular structure.

A classical calculator earns its place by explaining T2 structure
BEYOND what the environment data already captures. The measure: does
adding a calculator's T2 to the feature set reduce the residual
against DFT deltas? If not — the calculator's angular physics didn't
tell the calibration model anything APBS EFG and MOPAC charges hadn't
already provided.

This is the true bar: not "does the calculator produce T2?" (trivially
yes — the equations produce full tensors) but "does the calculator's T2
REDUCE the residual beyond what environment tensors already explain?"

Environment tensor data by tool:
- **APBS**: E-field (L=1, Vec3), EFG (L=2, Mat3 + sphericart)
- **MOPAC**: charges, orbital populations, bond orders (L=0 scalars)
- **ORCA DFT**: full shielding tensor (L=0+L=1+L=2, feature input)
- **DSSP**: L=0 only (SS, phi, psi, SASA)
- **Charge assignment**: L=0 only (partial charges)

The L=2 data from APBS means the model already sees angular
environment structure. A classical calculator must add angular physics
that these tools CANNOT capture — ring-specific geometry, bond-specific
anisotropy, distance-dependent dispersion — to justify its existence
at L=2.

### The bar for classical calculators

Each classical calculator must produce FULL TENSOR output (Mat3 +
sphericart, all irrep levels). Its contribution is measured by:
does adding this calculator's output as a feature improve the model's
prediction of the mutant delta, beyond what APBS + MOPAC + DSSP + DFT
baseline already capture?

If yes: the calculator adds physics the environment data missed.
If no: the calculator is redundant with the raw field data.

This is how we justify each classical computation scientifically,
not just numerically. The improvement must be in the TENSOR
structure (L=2), not just the scalar (L=0).

### DFT method mixing

Mutant deltas: r2SCAN (existing, not rerun, ground truth)
RefDB baselines: r2SCAN-3c (3-5x faster, feasible for 700 proteins)

Mixing is OK because:
- Mutant deltas are WITHIN-METHOD (WT and ALA both r2SCAN)
- Basis set errors cancel in the delta
- RefDB baselines are INPUT FEATURES, not training targets
- Model knows which method via ProteinBuildContext metadata
- NEVER subtract across methods (r2SCAN-3c minus r2SCAN = garbage)

## Calibration: How Parameters Are Tuned

The ~93 tuneable calculator parameters (ring current intensities, bond
anisotropies, Buckingham coefficients — see CALCULATOR_PARAMETER_API.md)
are calibrated by the external Python e3nn model (learn/c_equivariant/)
against DFT WT-ALA delta tensors. Calibrated values enter the C++
system as TOML configuration overriding literature defaults.

The C++ system is a transparent instrument: it computes geometric
kernels, loads DFT reference tensors when available, and exports both
via WriteFeatures() as NPY arrays. The calibration comparison happens
in the Python pipeline.

### The T2 residual: where classical angular physics breaks down

After classical calculators run with calibrated parameters, the L=2
residual (DFT delta T2 minus classical kernel T2) reveals where the
classical equations get the ANGULAR structure wrong, not just the
magnitude. This is a first-class thesis result:

- T2 residual small → classical angular physics is correct. The
  equations get the shape right; they just needed magnitude tuning.
- T2 residual systematically nonzero near specific ring types → the
  Johnson-Bovey or Haigh-Mallion angular approximation breaks down
  for those types.
- T2 residual large at scattered atoms → data quality issue (bad DFT
  convergence, atom matching error, geometry artefact).
- T2 residual large far from rings → missing physics mechanism.

The T2 residual map across proteins, per ring type, per calculator,
is a diagnostic no one has produced before. It tells us not just
WHETHER classical models work, but WHERE and WHY they fail angularly.
