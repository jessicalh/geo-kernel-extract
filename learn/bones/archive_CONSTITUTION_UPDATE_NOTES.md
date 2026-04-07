# Constitution Update Notes

These corrections supersede corresponding sections in CONSTITUTION.md
and OBJECT_MODEL.md. The next agent pass must integrate these.

---

## 1. Environment → ProteinBuildContext

The polymorphic Environment hierarchy is eliminated. Replace with:

**ProteinBuildContext**: How this protein instance was built. PDB source,
protonation tool + version + pH, force field, what was stripped, what
was assumed, crystal resolution, organism, deposition date. Friend to
the Protein. Travels with the protein instance. A copy may keep its
build context or take a new one.

**Tool configuration** (where binaries live, default parameters): global,
not per-protein.

**Computed results**: live on conformations as typed result objects (see §4).

---

## 2. Conformation → ProteinConformation (typed hierarchy)

A conformation without its protein is meaningless. Renamed to
ProteinConformation.

**Positions are const after construction.** New geometry = new
ProteinConformation. Spatial index is never stale.

**Class hierarchy, self-type-reporting:**

```
ProteinConformation (base — NOT abstract, fully functional)
├── ExperimentalConformation
│   ├── CrystalConformation
│   │     resolution_angstroms, r_factor, temperature_kelvin, pdb_id
│   └── NMRConformation
│         ensemble_member, restraint_count
├── ComputedConformation
│   ├── PredictionConformation
│   │     tool (AlphaFold/OpenFold3/ESMFold), confidence_per_residue
│   ├── MinimisedConformation
│   │     parent, method, force_field, energy, converged
│   └── MDFrameConformation
│         frame_index, time_picoseconds, walker_index, boltzmann_weight
└── DerivedConformation
      parent, derivation_description, properties_invalidated
```

The base class does the heavy lifting: holds positions (const),
provides spatial queries, accumulates result objects, supports the
full extraction pipeline. Subclasses add source-specific metadata.

---

## 3. Protein owns the conformation list

The protein creates conformations through typed factory methods.
No agent creates a ProteinConformation directly. No loose conformations.

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

Copy policy:
- GeometryOnly: copy positions, clear all results
- Full: copy everything (only valid if build context unchanged)

---

## 4. ConformationResult: named singletons with dependency graph

Computed results are NOT fields on a mega-struct. They are separate
typed objects that attach to a conformation and never detach.

Each result type is a **named singleton** on the conformation.
Accessed by name, not by iteration:

```
auto& dssp = conformation.Dssp();           // DsspResult&
auto& apbs = conformation.ApbsField();      // ApbsFieldResult&
auto& bs = conformation.BiotSavart();       // BiotSavartResult&
auto& coulomb = conformation.Coulomb();     // CoulombResult&
```

Each accessor returns the real typed object. Throws if not yet
computed. The code reads like English:
  conformation.BiotSavart().TensorAtAtom(42, ring3)

Each result type declares its dependencies:

```
class McConnellResult : ConformationResult {
    requires: {SpatialIndexResult}
    ...
};
```

At attach time, dependencies are checked. Missing dependency =
immediate logged error, not a runtime surprise.

**Once attached, permanent.** A conformation that has a DsspResult
keeps it for its lifetime. Results accumulate. Nothing is removed.

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

Known result types (with dependencies):

```
DsspResult                requires: nothing
ChargeAssignmentResult    requires: nothing
SpatialIndexResult        requires: nothing
ApbsFieldResult           requires: ChargeAssignmentResult
MolecularGraphResult      requires: SpatialIndexResult
BiotSavartResult          requires: SpatialIndexResult
HaighMallionResult        requires: SpatialIndexResult
McConnellResult           requires: SpatialIndexResult
CoulombResult             requires: ChargeAssignmentResult, SpatialIndexResult
HBondResult               requires: DsspResult
DispersionResult          requires: SpatialIndexResult
PiQuadrupoleResult        requires: SpatialIndexResult
RingSusceptibilityResult  requires: SpatialIndexResult
OrcaShieldingResult       requires: nothing (loaded from files)
FeatureExtractionResult   requires: all physics results
PredictionResult          requires: FeatureExtractionResult
```

---

## 5. Extraction order derived from dependencies

The manual phase ordering (Phase 0-6) is replaced by the
dependency graph. The pipeline resolves the order:

1. Check which results are attached
2. Find results whose dependencies are all satisfied
3. Compute and attach them
4. Repeat until all desired results are attached

The order emerges from the science (can't run Coulomb before
charges) not from a numbered list. An agent adding a new result
type declares dependencies and the pipeline slots it in.

The EXTRACTION_ORDER.md document becomes: here are the known result
types, their dependencies, and what each stores on atoms/rings.
The order is implied.

---

## 6. Conformation types are extensible

Adding a new conformation type:
1. Define subclass with type-specific metadata
2. Add factory method on Protein
3. Done

Adding a new result type:
1. Define subclass of ConformationResult
2. Declare dependencies
3. Implement computation
4. Add named accessor on ProteinConformation
5. Done

The base ProteinConformation provides all the machinery. New types
inherit it. The goal: a feature extraction expert grabs a fresh
notebook (new conformation type or result type) and all the pencils
and calculators are already there.

---

## 7. Implementation process: incremental agent passes

We do NOT implement everything in one shot. Agents run in layers:

**Layer 1**: Core object model. Protein, ProteinConformation,
ProteinBuildContext, Atom, Residue, Bond, Ring hierarchy. No
physics. Just the structural objects with typed properties,
spatial queries, and the result attachment framework. Unit tested.

**Layer 2**: First result types. SpatialIndexResult, DsspResult,
ChargeAssignmentResult. These are the foundation that everything
else depends on. Attach to conformations from Layer 1. Unit tested.

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
- The object model may need corrections from reality the agent found
  -- that's OK, but the correction goes into the model, not into a
  hack in the agent's code

This way complexity is introduced incrementally. Each layer is
reviewable. Violations are caught before they compound. The object
model evolves from contact with reality but stays coherent.
