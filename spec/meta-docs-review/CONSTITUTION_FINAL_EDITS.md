# Constitution: Final Replacement Text

These are the replacement sections for CONSTITUTION.md, written in
the constitution's own voice. True statements only. No TBDs, no
"planned," no aspirational language. When these are applied, the
constitution is true again.

The reasoning for each edit is in CONSTITUTION_PROPOSED_EDITS.md
(the earlier draft with notes). This file is the clean text.

---

## Replace: Protonation State (lines 49-52)

### Protonation State

Protonation determines the atom list. Which titratable groups are
charged, which histidine tautomer exists, which cysteines form
disulfide bonds -- these decisions change which atoms are present.
Adding or removing hydrogen atoms changes the flat index array that
every ring vertex, bond endpoint, ConformationAtom, and spatial
neighbour list is indexed into.

**Different protonation = different Protein.** Not a different
conformation of the same protein. A ProteinConformation holds
geometry for a fixed atom list. Many conformations share the same
Protein (same protonation, different poses). Two protonation states
require two Protein instances.

Protonation state is recorded in the ProteinBuildContext (which tool,
which pH, which decisions) and in the Residue variant indices (which
variant was applied to each titratable residue).

---

## Replace: ProteinConformation hierarchy (lines 94-112)

**Class hierarchy:**

```
ProteinConformation (base -- NOT abstract, fully functional)
+-- CrystalConformation
|     resolution_angstroms, r_factor, temperature_kelvin, pdb_id
+-- PredictionConformation
|     method, confidence
+-- MDFrameConformation
|     walker, time_picoseconds, boltzmann_weight, rmsd_nm, rg_nm
+-- DerivedConformation
      derivation_description
```

The hierarchy is flat. The base class does the heavy lifting: holds
positions (const), provides ConformationResult attachment and template
access, holds the ConformationAtom wide table. Subclasses add
provenance metadata only. Future types are direct subclasses of
ProteinConformation, following the same pattern.

**Primary access**: `protein.Conformation()` returns the first
conformation as ProteinConformation& regardless of subtype. This is
what calculators and tests use. Typed accessors (`CrystalConf()`,
`MDFrameAt()`, `PredictionAt()`) exist for metadata consumers.

---

## Insert after ProteinConformation hierarchy: Protein-Conformation
   Access Path

### Accessing the Protein from a Conformation

A ProteinConformation holds a valid pointer to its owning Protein.
This is the designed access path for identity information:

```cpp
// Identity (element, bonds, residue type) -- from Protein
const auto& protein = conf.ProteinRef();
Element elem = protein.AtomAt(i).element;
const Residue& res = protein.ResidueAt(protein.AtomAt(i).residue_index);
const Ring& ring = protein.RingAt(ri);

// Geometry and computed fields -- from ProteinConformation
Vec3 pos = conf.PositionAt(i);
const auto& ca = conf.AtomAt(i);  // ConformationAtom
double charge = ca.partial_charge;
```

The Protein is non-copyable and non-movable. It lives on the heap
via unique_ptr. Conformations are owned by the Protein and cannot
outlive it. The raw pointer is safe because the lifetime guarantee
is structural.

**Do not** duplicate identity onto ConformationAtom. Element, residue
type, bond connectivity, ring membership live on Protein and are
accessed through ProteinRef(). This is the separation between what
the molecule IS (Protein) and where its atoms ARE (ProteinConformation).

---

## Replace: Copy Semantics section (lines 551-635)

## Protein Rebuild Pattern

### The principle

Analysing the same protein geometry under different chemical
conditions is a first-class operation:

- Different protonation (pH scanning)
- Different mutation (WT vs ALA)
- Different charge model (ff14SB vs xTB vs CHARMM36m)

The delta between two conditions IS the experiment. WT protein with
full extraction. ALA mutant with full extraction. Delta in DFT
shielding is the training target. Delta in classical calculations is
what the model learns to predict.

### The mechanism: rebuild, not copy

Different protonation changes the atom list. Different atom list
means different flat indices. The correct mechanism is a clean
rebuild from symbolic data:

1. Take heavy-atom identity from the source (residue types, atom
   names, sequence)
2. Apply new protonation decisions
3. Produce the new atom list (with hydrogens placed by the
   protonation tool)
4. Build a new Protein with fresh indices via FinalizeConstruction
5. Create conformations using heavy-atom position transfer (symbolic
   matching by residue + atom name bridges old and new index spaces)
6. Run the full extraction pipeline

MutationDeltaResult demonstrates this: two separate Proteins (WT and
mutant), atom matching by position + element, delta tensors computed.
The three loaders (PdbFileReader, OrcaRunLoader, GromacsEnsembleLoader)
each build fresh Proteins from different source formats, converging
on the same FinalizeConstruction pipeline.

### Why not copy

The two indexing systems (symbolic identity and geometric flat array)
are fused at FinalizeConstruction. There is no safe way to patch one
without rebuilding the other. See spec/PROTONATION_DESIGN_HISTORY.md.

---

## Replace: Framework Stores section (lines 846-877) + UDP Logging
   section (lines 971-1015)

## Enforcement: Structural Guarantees and Tracing

### The singleton guarantee

Each ConformationResult type is a singleton on a ProteinConformation.
BiotSavartResult is computed once, attached once, never replaced.
Each ConformationAtom field has exactly one writer, determined by
which ConformationResult type owns it. AttachResult checks:

- Not null
- Not already attached (singleton violation: logged error, rejected)
- All dependencies present (missing dependency: logged error, rejected)

This prevents overwrites and ordering errors structurally. The
scenario where two results write the same field cannot arise.

### Filter rejection logging

Every KernelFilterSet logs per-evaluation rejections with specific
context:

```
RingBondedExclusionFilter: atom=47 is bonded to ring 3
DipolarNearFieldFilter: atom=123 distance=0.82A < threshold=1.15A
  (source_extent=2.30A)
SelfSourceFilter: atom=15 IS source_atom=15
```

Each record includes: filter name, atom index, source indices,
distance, source extent, human-readable reason. Per-filter rejection
counts are reported at computation end. This is the "what happened to
atom 42" diagnostic for kernel evaluations.

### Operation logging

OperationLog (OperationLog.h) provides structured JSON logging over
UDP or stderr with 22 named channels. Every result attachment,
external tool invocation, and computation boundary emits a log entry.
Scoped logging tracks elapsed time automatically. Channel bitmask
controls verbosity. Warnings and errors always emit regardless of
channel mask.

---

## Replace: Known result types dependency graph (lines 349-371)

### Known result types (with dependencies)

```
GeometryResult                requires: nothing
DsspResult                    requires: nothing
ChargeAssignmentResult        requires: nothing
EnrichmentResult              requires: nothing
OrcaShieldingResult           requires: nothing (loaded from files)
SpatialIndexResult            requires: GeometryResult
ApbsFieldResult               requires: ChargeAssignmentResult
XtbChargeResult               requires: nothing (external tool)
MolecularGraphResult          requires: SpatialIndexResult
BiotSavartResult              requires: SpatialIndexResult, GeometryResult
HaighMallionResult            requires: SpatialIndexResult, GeometryResult
McConnellResult               requires: SpatialIndexResult, GeometryResult
CoulombResult                 requires: ChargeAssignmentResult, SpatialIndexResult
HBondResult                   requires: DsspResult, SpatialIndexResult
DispersionResult              requires: SpatialIndexResult, GeometryResult
PiQuadrupoleResult            requires: SpatialIndexResult, GeometryResult
RingSusceptibilityResult      requires: SpatialIndexResult, GeometryResult
MutationDeltaResult           requires: OrcaShieldingResult (on both WT and mutant)
ProtonationDetectionResult    requires: nothing
```

The dependency graph drives execution order. Each result type declares
its dependencies via `Dependencies()`. AttachResult checks them at
attach time. Missing dependency = immediate diagnostic, not a runtime
surprise.

---

## Replace: Inviolable Rule 7 (line 524)

7. **Feature output is per-result.** Each ConformationResult knows
   what it computed and writes its own features to disk via
   WriteFeatures(). ConformationResult::WriteAllFeatures() traverses
   all attached results. No centralised feature extraction step
   rewrites what the calculators already stored.

---

## Replace: External Libraries (lines 537-548)

## External Libraries

- **Eigen3**: Vec3, Mat3, SVD
- **nanoflann**: spatial KD-tree, neighbour queries
- **sphericart**: spherical harmonics, canonical normalization
- **cifpp**: PDB/mmCIF parsing
- **libdssp**: secondary structure, H-bonds
- **OpenBabel**: bond perception, hybridisation
- **APBS**: Poisson-Boltzmann solvation (C library bridge)
- **libgromacs**: GROMACS TPR reading (topology, charges, elements)
