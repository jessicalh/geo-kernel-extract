# Object Model: NMR Shielding Tensor Prediction System

This is the concrete object model. Every class, every property, every type,
every unit, every dependency. Implementation agents code from this.
If a property is not in this document, it does not exist in the code.

**NOTE (2026-04-03): This document is not aspirational. The library is
mature for this functionality — 4 builders, 4 use cases,
293 tests passing. Read the code before assuming any section is current.**

**NOTE (2026-04-13):** SasaResult, AIMNet2Result, WaterFieldResult,
HydrationShellResult, GromacsEnergyResult, and the extended DsspResult
(8-class SS, H-bond energies, chi1-4) are now documented in
"Complete property flow" below.

**Trajectory-scope entities** (TrajectoryProtein, TrajectoryAtom,
TrajectoryResult, Trajectory, Session, RunConfiguration, RecordBag,
DenseBuffer) are documented below after the ConformationResult section.
The old Gromacs* classes (GromacsProtein, GromacsProteinAtom,
GromacsRunContext, AnalysisWriter) are in `learn/bones/` — their
replacement landed as the trajectory-scope object model. Design
working-notes and pending appendices (NmrAtomIdentity, full
TrajectoryResult catalog, H5 metadata schema) are in
`spec/pending_include_trajectory_scope_2026-04-22.md`; that file is
not authoritative for anything already landed — the code + the
trajectory-scope section below win.

**Copy-and-modify pattern: SUPERSEDED.** The copy-and-modify sections in
this document are design history. The pattern was originally proposed for
pH scanning and re-protonation but was never needed for the actual use
cases (A-D in USE_CASES.md). The system creates one Protein per input
source with N conformations, each independently enriched. There is no
protein copying. Do not implement copy-and-modify. The sections remain
as design history documenting how the Protein/ProteinConformation
separation evolved — the distinction they enforced is real and load-bearing,
the mechanism described is not.

**Calibration pipeline.** The ~93 tuneable calculator parameters are
calibrated externally by the Python e3nn model (learn/c_equivariant/)
against DFT WT-ALA delta tensors. Calibrated values enter the C++
system as TOML configuration overriding literature defaults. There is
no embedded model in the C++ library.

Aligned with CONSTITUTION.md (2026-03-29 revision). Key changes from prior
version: Environment replaced by ProteinBuildContext, flat Conformation
replaced by typed ProteinConformation hierarchy, manual phases replaced by
dependency-graph-driven ConformationResult singletons, results accessed by
name with physics query methods.

---

## Core Types (from Eigen, used everywhere)

```
Vec3 = Eigen::Vector3d     // 3-component vector
Mat3 = Eigen::Matrix3d     // 3x3 matrix
```

### SphericalTensor
Irreducible decomposition of a 3x3 tensor via sphericart.

| Property | Type | Description |
|----------|------|-------------|
| T0 | double | Isotropic component (trace / 3) |
| T1 | array<double, 3> | Antisymmetric (pseudovector) components |
| T2 | array<double, 5> | Traceless symmetric components (m = -2..+2) |

Static methods:
- `Decompose(Mat3) -> SphericalTensor`: via sphericart, canonical normalization
- `Reconstruct() -> Mat3`: inverse

Both Mat3 AND SphericalTensor are stored. No conversion at point of use.

### FieldValue
A calculator result attributed to a specific source.

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| tensor | Mat3 | depends on calculator | Full 3x3 tensor |
| spherical | SphericalTensor | same | Sphericart decomposition |
| source_calculator | CalculatorId enum | - | Which calculator |
| source_index | size_t | - | Index of source ring/bond/atom |

Note: no string identifier for the source. CalculatorId enum is
sufficient. Strings for identity are forbidden (see PATTERNS.md).

### KernelEvaluationFilter (ABC) and KernelFilterSet

Every calculator that evaluates a geometric kernel holds a
KernelFilterSet: an ordered list of KernelEvaluationFilter objects.
Before computing a kernel at a field point, the calculator assembles a
KernelEvaluationContext from the evaluation geometry and passes it to
the filter set. All filters must accept for the kernel to proceed.

Filters exist for first-principles physics reasons, not for empirical
tuning. The test for a filter: the calculator works but not as well
without it. A filter that changes nothing is dead code.

```
KernelEvaluationContext
    distance           double     field point to source center (A)
    source_extent      double     spatial extent of source (A)
    sequence_separation int       residue seq distance (-1 if N/A)
    atom_index         size_t     field point atom
    source_atom_a      size_t     first source atom (SIZE_MAX if N/A)
    source_atom_b      size_t     second source atom (SIZE_MAX if N/A)

KernelEvaluationFilter (ABC)
    Accept(ctx) -> bool           should this evaluation proceed?
    Name() -> const char*         identifier for logging and TOML
    Description() -> const char*  physics reason this filter exists
    RejectReason(ctx) -> string   what/where/which for this rejection
                                  (specific values from ctx, not generic)

KernelFilterSet
    Add(filter)                   add a filter (ownership transferred)
    AcceptAll(ctx) -> bool        all must accept; tracks rejections
    ReportRejections() -> string  per-filter rejection counts
    SetLogRejections(bool)        enable per-rejection detail logging
    RejectionReasons() -> vector<string>  human-readable what/where/which
```

Concrete filters (KernelEvaluationFilter.h):

| Filter | Criterion | Physics |
|--------|-----------|---------|
| DipolarNearFieldFilter | distance > source_extent/2 | Multipole expansion invalid inside source distribution |
| SelfSourceFilter | atom_index ∉ {source_a, source_b} | Field undefined at source itself |
| SequentialExclusionFilter | seq_separation >= min_sep | Through-bond coupling not modelled by through-space kernel |

Filter sets per calculator:

| Calculator | Filters | Source extent |
|------------|---------|---------------|
| McConnell | SelfSource + DipolarNearField | bond length |
| CoulombResult | SelfSource (i≠j only) | 0 (point source) |
| RingSusceptibility | DipolarNearField | ring diameter |
| HBondResult | SelfSource + DipolarNearField | N...O distance |
| BiotSavartResult | DipolarNearField | ring diameter |
| HaighMallionResult | DipolarNearField | ring diameter |

Filter sets are configurable from TOML. Rejection counts are logged
per filter per calculator invocation. Per-rejection detail (which atom,
which source, at what distance) is available when diagnostic logging
is enabled.

---

## Element (enum)

```
enum class Element { H, C, N, O, S, Unknown };
```

Properties (compile-time, per element):
| Element | Atomic number | Covalent radius (A) | Electronegativity (Pauling) |
|---------|---------------|--------------------|-----------------------------|
| H | 1 | 0.31 | 2.20 |
| C | 6 | 0.76 | 2.55 |
| N | 7 | 0.71 | 3.04 |
| O | 8 | 0.66 | 3.44 |
| S | 16 | 1.05 | 2.58 |

---

## Hybridisation (enum)

```
enum class Hybridisation { sp, sp2, sp3, Unassigned };
```

Set by OpenBabel during enrichment. Const after enrichment.

---

## AtomRole (enum)

```
enum class AtomRole {
    // Heavy backbone
    BackboneN,      // peptide nitrogen
    BackboneCA,     // alpha carbon
    BackboneC,      // carbonyl carbon
    BackboneO,      // carbonyl oxygen

    // Heavy sidechain
    SidechainC,     // sidechain carbon (non-aromatic)
    SidechainN,     // sidechain nitrogen
    SidechainO,     // sidechain oxygen
    SidechainS,     // sulfur (CYS, MET)
    AromaticC,      // carbon in an aromatic ring
    AromaticN,      // nitrogen in an aromatic ring (HIS, TRP)

    // Hydrogen (by NMR-relevant environment)
    AmideH,         // bonded to backbone peptide N
    AlphaH,         // bonded to CA
    MethylH,        // bonded to terminal CH3 carbon
    AromaticH,      // bonded to aromatic ring carbon
    HydroxylH,      // bonded to OH (SER, THR, TYR)
    OtherH,         // all other hydrogens

    Unknown
};
```

### Assignment criteria (from enrichment, no string parsing)

| Role | Element | Bond connectivity | Residue context |
|------|---------|-------------------|-----------------|
| BackboneN | N | bonded to CA of same residue and C of previous residue | backbone position from Residue.N index |
| BackboneCA | C | bonded to backbone N and backbone C | backbone position from Residue.CA index |
| BackboneC | C | bonded to backbone O (double bond) and backbone N+1 | backbone position from Residue.C index |
| BackboneO | O | bonded to backbone C (double bond) | backbone position from Residue.O index |
| SidechainC | C | not backbone, not in aromatic ring | all other carbons |
| SidechainN | N | not backbone, not in aromatic ring | sidechain nitrogens |
| SidechainO | O | not backbone | sidechain oxygens |
| SidechainS | S | any | CYS SG, MET SD |
| AromaticC | C | member of an AromaticRing vertex set | detected from ring membership |
| AromaticN | N | member of an AromaticRing vertex set | HIS ND1/NE2, TRP NE1 |
| AmideH | H | bonded to BackboneN | Residue.H index (not PRO) |
| AlphaH | H | bonded to BackboneCA | Residue.HA index (or HA2 for GLY) |
| MethylH | H | bonded to a terminal sp3 carbon with 3 H neighbours | parent C has exactly 3 H bonds and 1 C/S bond |
| AromaticH | H | bonded to AromaticC | parent is ring member |
| HydroxylH | H | bonded to O that is bonded to C | SER OG, THR OG1, TYR OH |
| OtherH | H | none of the above | remaining hydrogens |

---

## BondCategory (enum)

```
enum class BondCategory {
    PeptideCO,      // backbone C=O (C to O, double bond)
    PeptideCN,      // backbone C-N (C to N+1, partial double)
    BackboneOther,  // N-CA, CA-C, CA-CB
    SidechainCO,    // sidechain C=O (ASN, ASP, GLN, GLU)
    Aromatic,       // within aromatic ring
    Disulfide,      // S-S crosslink
    SidechainOther, // all other sidechain bonds
    Unknown
};
```

Note: The constitution specifies five categories (PeptideCO, PeptideCN,
Sidechain, Aromatic, Disulfide). We retain the finer granularity
(BackboneOther, SidechainCO, SidechainOther) because McConnell bond
anisotropy calculations require distinguishing double-bond sidechain
carbonyls from single-bond sidechain connections, and backbone N-CA/CA-C
bonds from peptide bonds. The coarser groupings can be derived from
these at query time.

### Assignment criteria (from bond detection)

| Category | Atom A role | Atom B role | Bond order |
|----------|-------------|-------------|------------|
| PeptideCO | BackboneC | BackboneO | Double |
| PeptideCN | BackboneC (residue i) | BackboneN (residue i+1) | Peptide |
| BackboneOther | any backbone | any backbone | Single |
| SidechainCO | SidechainC | SidechainO | Double |
| Aromatic | AromaticC/AromaticN | AromaticC/AromaticN | Aromatic |
| Disulfide | SidechainS | SidechainS | Single (crosslink) |
| SidechainOther | any sidechain | any sidechain | any |

---

## HeuristicTier (enum)

```
enum class HeuristicTier { REPORT, PASS, SILENT };
```

---

## ProteinBuildContext

How this protein instance was built. Immutable once constructed.
Records provenance (what happened), not experimental conditions
(which live on the ProteinConformation subtype metadata).

**Why this contract exists (for agents creating new build contexts):**
Every protein instance needs a build context because it is the mutable
part of a protein definition. When we copy a protein to test different
protonation, the build context tells us what changed. When we compare
WT and mutant, it tells us which is which. When we serialise results,
it is the provenance record. Future work you can't see will need to
rely on this. Provide the basics so downstream consumers don't invent
their own tracking.

### Static properties (const after construction)

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| pdb_source | string | - | PDB ID or file path |
| deposition_date | string | - | PDB deposition date (if known) |
| organism | string | - | Source organism (if known) |
| crystal_resolution | double | Angstroms | NaN if not crystallographic |
| protonation_tool | string | - | "PROPKA 3.5.1", "KaML 1.0", "Manual" |
| protonation_pH | double | - | pH used for protonation prediction |
| force_field | string | - | "ff14SB", "CHARMM36m", etc. |
| stripped | vector<string> | - | What was removed ("waters", "heteroatoms") |
| assumptions | vector<string> | - | What was assumed ("missing loops rebuilt") |

### Query methods
- `Clone() -> unique_ptr<ProteinBuildContext>`
- `Describe() -> string`: human-readable provenance summary

### Copy semantics
Deep copy. All properties preserved.

---

## ProtonationState

Per-protein protonation decisions. Value type. Created by a Protonator
(PROPKA, KaML, tleap default) and consumed by a topology builder to
produce a Protein with the correct hydrogen atoms. Changing protonation
means a new Protein (copy-and-modify), not mutation of an existing one.

### Static properties (const after construction)

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| name | string | - | e.g. "propka_pH7.0" |
| pH | double | - | NaN for default |
| tool | ProtonationTool enum | - | PROPKA, KaML, TLeap, Manual |
| tool_version | string | - | "PROPKA 3.5.1", "KaML-CBTrees" |
| residues | vector<ResidueProtonation> | - | Per-residue decisions |

### ResidueProtonation (value type)

| Property | Type | Description |
|----------|------|-------------|
| residue_index | size_t | Into protein's residue list |
| amino_acid | AminoAcid enum | Typed, not string |
| variant_index | int | Into AminoAcidType::variants, -1 = default |
| pKa | double | Predicted pKa (NaN if not computed) |
| is_charged | bool | Carries formal charge in this state |

### Query methods
- `ForResidue(size_t) -> const ResidueProtonation*`
- `NetDecisionCharge() -> int` (only counts explicit decisions)
- `NetChargeForProtein(vector<AminoAcid>) -> int` (full protein including non-titratable)
- `Describe() -> string`

### Copy semantics
Value type. Full copy.

## ChargeSource (typed hierarchy)

Where per-atom charges come from. The force field determines charges,
naming, and VdW parameters. Each source is a distinct type, not a
string-dispatched function.

```
ChargeSource (abstract)
├── ParamFileChargeSource   — ff14SB from flat parameter file (fallback)
├── PrmtopChargeSource      — ff14SB/ff19SB from AMBER prmtop (authoritative)
├── GmxTprChargeSource      — CHARMM36m from GROMACS .tpr (fleet data)
└── StubChargeSource        — uniform test charges
```

ChargeAssignmentResult::Compute(conf, ChargeSource&) is the typed factory.
The param-file and stub convenience factories delegate to it.

### ForceField (enum)

```
enum class ForceField { Amber_ff14SB, Amber_ff19SB, CHARMM36m, Unknown };
```

Recorded in ProteinBuildContext as provenance.

---

## Protein

Sequence and science data. What the molecule IS, independent of geometry.
Does NOT hold positions. Does NOT hold computed properties.

The protein owns its conformation list. Conformations are created
through typed factory methods on the protein:

```
protein.AddCrystalConformation(positions, metadata) -> CrystalConformation&
protein.AddNMRConformation(positions, metadata) -> NMRConformation&
protein.AddMDFrame(positions, metadata) -> MDFrameConformation&
protein.AddPrediction(positions, metadata) -> PredictionConformation&
protein.AddMinimised(positions, metadata) -> MinimisedConformation&
protein.AddDerived(parent, description) -> DerivedConformation&

protein.CrystalConformation()    // exactly one or throws
protein.NMRConformations()       // typed collection
protein.MDFrames()               // typed collection
protein.Predictions()            // typed collection
protein.Conformations()          // all, heterogeneous
```

No agent creates a ProteinConformation directly. No loose conformations.

### Static properties (const after construction)

| Property | Type | Unit | Set at |
|----------|------|------|--------|
| residues | vector<Residue> | - | construction |
| atoms | vector<unique_ptr<Atom>> | - | construction |
| aromatic_rings | vector<Ring> | - | construction (topology only) |
| covalent_bonds | vector<Bond> | - | construction |
| build_context | unique_ptr<ProteinBuildContext> | - | construction |
| conformations | vector<unique_ptr<ProteinConformation>> | - | via factory methods |

### Query methods
- `AtomAt(size_t) -> const Atom&`
- `ResidueAt(size_t) -> const Residue&`
- `RingAt(size_t) -> const Ring&`
- `BondAt(size_t) -> const Bond&`
- `AtomCount() -> size_t`
- `ResidueCount() -> size_t`
- `RingCount() -> size_t`
- `BondCount() -> size_t`

### Construction: FinalizeConstruction

Every loader must call `FinalizeConstruction(positions)` after adding
all atoms and residues. This single call performs, in order:
1. `CacheResidueBackboneIndices()` — backbone N/CA/C/O/H/HA indices
2. `DetectAromaticRings()` — from AminoAcidType + atom presence
3. `DetectCovalentBonds(positions)` — via OpenBabel bond perception

Order matters: rings need backbone cache (to know residue types),
bonds need rings (to classify aromatic bonds). Call BEFORE creating
any ProteinConformation.

### Ownership and lifetime
Protein is non-movable and non-copyable. It lives on the heap via
unique_ptr, created by LoadProtein() which returns unique_ptr<Protein>
in a LoadResult. Conformations hold a raw Protein* back-pointer that
is valid for the Protein's lifetime. Non-movable means the pointer
never dangles. No move constructors, no move assignment.

For the copy-and-modify pattern (pH scanning, mutant comparison):
use a Copy() factory method that creates a NEW Protein on the heap
with modified protonation/build context:
- `Copy(new_protonation)`: new protein, applies new protonation,
  invalidates ring types for titratable residues, invalidates charges
- `Copy(new_build_context)`: new protein, applies new build context
- Geometry-only properties (positions, distances) transfer on copy.
  Charge-dependent properties (APBS, Coulomb, features) invalidated.

### Relationship to other objects
- Owns Residues, Atoms, Rings, Bonds, ProteinBuildContext, ProteinConformations
- Does NOT own computed properties (those live on ProteinConformations
  as ConformationResult objects)

---

## Residue

One amino acid at a sequence position.

### Static properties (const after construction)

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| type | AminoAcid enum | - | Standard amino acid type |
| sequence_number | int | - | PDB sequence number |
| chain_id | string | - | Chain addressing (not entity) |
| insertion_code | string | - | PDB column 27 (e.g. "A", "B") |
| atom_indices | vector<size_t> | - | Into protein atom list |
| protonation_variant_index | int | - | -1 if not titratable |
| N | size_t | - | Backbone N atom index, NONE if absent |
| CA | size_t | - | Alpha carbon atom index |
| C | size_t | - | Carbonyl carbon atom index |
| O | size_t | - | Carbonyl oxygen atom index |
| H | size_t | - | Amide H atom index, NONE for PRO |
| HA | size_t | - | Alpha H atom index (HA2 for GLY) |
| CB | size_t | - | Beta carbon atom index, NONE for GLY |
| chi[4] | ChiAtoms | - | Dihedral atom indices |

DSSP properties are conformation-dependent. They live on the
DsspResult ConformationResult, accessed via conformation.Dssp().

### Query methods
- `IsAromatic() -> bool` (from AminoAcidType)
- `IsTitratable() -> bool`
- `HasAmideH() -> bool`
- `ChiAngleCount() -> int`
- `SequenceAddress() -> (chain_id, sequence_number, insertion_code)`

### Copy semantics
Value type within Protein. Copies with Protein.

---

## Atom

Each atom in the protein. Identity only -- no position.
Position lives in a ProteinConformation.

### Static properties (const after construction)

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| element | Element enum | - | Atomic element |
| pdb_atom_name | string | - | Display/serialisation ONLY |
| residue_index | size_t | - | Into protein's residue list |
| bond_indices | vector<size_t> | - | Into protein's bond list |

Element-derived properties (non-virtual, from free functions in Types.h):
| Property | Type | Unit | Description |
|----------|------|------|-------------|
| CovalentRadius() | double | Angstroms | From element |
| Electronegativity() | double | Pauling | From element |
| IsHBondDonorElement() | bool | - | N, O return true |
| IsHBondAcceptorElement() | bool | - | N, O return true |
| parent_atom_index | size_t | - | Nearest bonded heavy atom (SIZE_MAX if not H) |

Note: there is NO Atom subclass hierarchy. One flat Atom class.
CovalentRadius() etc. dispatch on the element enum via free functions.
parent_atom_index is on the base class, SIZE_MAX for non-hydrogen atoms.
This was a deliberate simplification from Pass 1: the old subclasses
(HydrogenAtom, CarbonAtom, etc.) added virtual dispatch for no benefit.

### Enrichment properties (set once, append-only)

Set by enrichment ConformationResult objects during attachment.
Once set, never overwritten or removed within a ProteinConformation.
Each property records its source result type.

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| role | AtomRole enum | - | EnrichmentResult | NMR-relevant classification |
| hybridisation | Hybridisation enum | - | EnrichmentResult | From OpenBabel |
| is_backbone | bool | - | EnrichmentResult | role in backbone set |
| is_amide_H | bool | - | EnrichmentResult | role == AmideH |
| is_alpha_H | bool | - | EnrichmentResult | role == AlphaH |
| is_methyl | bool | - | EnrichmentResult | role == MethylH |
| is_aromatic_H | bool | - | EnrichmentResult | role == AromaticH |
| is_on_aromatic_residue | bool | - | EnrichmentResult | residue.IsAromatic() |
| is_hbond_donor | bool | - | EnrichmentResult | H bonded to N or O |
| is_hbond_acceptor | bool | - | EnrichmentResult | N or O with lone pair |
| parent_is_sp2 | bool | - | EnrichmentResult | For H: parent hybridisation == sp2 |

### ConformationAtom: the per-atom computed data store

Each ProteinConformation holds a `vector<ConformationAtom>` parallel to
the Protein's atom list. ConformationAtom has typed, named fields for
every per-atom result. Fields are declared upfront with default values.
ConformationResult objects fill them in during Attach().

**Construction and ownership:**

ConformationAtom has a PRIVATE constructor. Only ProteinConformation can
create them. The vector is built once in the ProteinConformation
constructor from the Protein's atom count and the provided positions,
and is never resized. You cannot construct a ConformationAtom outside
of a ProteinConformation. The compiler enforces this.

```cpp
class ConformationAtom {
    friend class ProteinConformation;
public:
    Vec3 Position() const { return position_; }

    // Computed fields below — written by ConformationResult::Attach(),
    // read by feature extraction and query methods. Public because the
    // singleton guarantee means each field has exactly one writer.
    AtomRole role = AtomRole::Unknown;
    // ... all fields listed in tables below ...

private:
    explicit ConformationAtom(Vec3 pos) : position_(pos) {}
    const Vec3 position_;  // set once at construction, never changed
};
```

**Identity goes through Protein, not ConformationAtom.** Element, bonds,
residue membership, ring membership — these are identity properties on
the Protein's Atom objects. ConformationAtom does NOT duplicate them.
Query identity through the ProteinConformation's back-pointer to Protein:

```cpp
// Identity (from Protein — never changes with geometry)
const auto& identity = conf.Protein().AtomAt(42);
Element elem = identity.element;
for (size_t bi : identity.bond_indices) {
    const auto& bond = conf.Protein().BondAt(bi);
}

// Computed data (from ConformationAtom — geometry-dependent)
const auto& ca = conf.AtomAt(42);
double charge = ca.partial_charge;
Vec3 B = ca.total_B_field;
SphericalTensor bs = ca.bs_shielding_contribution;
```

**T2 completeness rule:** Every field that stores a scalar (T0) quantity
derived from a tensor MUST have a companion SphericalTensor field storing
the full decomposition. If a calculator produces a Mat3, the
SphericalTensor is stored alongside it. No data is left on the floor.
Upstream models will explore this data for physics work — they need full
tensor information at every irrep level (L=0, L=1, L=2), not just
scalars. A ConformationAtom with scalar-only fields is INCOMPLETE.

**Adding new fields:** If a new layer's ConformationResult produces
per-atom data, you may add new fields to ConformationAtom. Do not
arbitrarily subdivide by physics categories — keep typed by the
functional analysis performed. BiotSavartResult fields, McConnellResult
fields, CoulombResult fields. The result type is the organising
principle.

#### Core geometry

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| position | Vec3 | Angstroms | (construction) | In conformation coordinate frame |

#### Charges and radii (from ChargeAssignmentResult, MopacResult)

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| partial_charge | double | elementary charge (e) | ChargeAssignmentResult | ff14SB force field |
| vdw_radius | double | Angstroms | ChargeAssignmentResult | ff14SB force field |
| mopac_charge | double | elementary charge (e) | MopacResult | PM7 Mulliken charge |
| mopac_s_pop | double | electrons | MopacResult | s-orbital population |
| mopac_p_pop | double | electrons | MopacResult | p-orbital population |
| mopac_valency | double | dimensionless | MopacResult | Sum of Wiberg bond orders |
| mopac_bond_neighbours | vector&lt;MopacBondNeighbour&gt; | - | MopacResult | Sorted descending by order |

### Spatial neighbourhood (on ProteinConformation)

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| spatial_neighbours | vector<AtomNeighbour> | - | All atoms within 15A |

AtomNeighbour:
| Field | Type | Unit |
|-------|------|------|
| atom_index | size_t | - |
| distance | double | Angstroms |
| direction | Vec3 | normalised |

### Ring neighbourhood (structured per nearby ring)

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| ring_neighbours | vector<RingNeighbourhood> | - | BiotSavartResult et al. | Per nearby ring |

See RingNeighbourhood class below.

### Bond neighbourhood (structured per nearby bond)

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| bond_neighbours | vector<BondNeighbourhood> | - | McConnellResult | Per nearby bond |

See BondNeighbourhood class below.

### Field accumulation (from ConformationResult objects)

Ring current totals (populated by BiotSavartResult, HaighMallionResult):
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| total_B_field | Vec3 | dimensionless | BiotSavartResult |
| total_G_tensor | Mat3 | dimensionless | BiotSavartResult |
| total_G_spherical | SphericalTensor | dimensionless | BiotSavartResult |
| per_type_G_T0_sum | array<double, 8> | dimensionless | BiotSavartResult |
| per_type_G_T2_sum | array<array<double,5>, 8> | dimensionless | BiotSavartResult |
| per_type_hm_T0_sum | array<double, 8> | dimensionless | HaighMallionResult |
| per_type_hm_T2_sum | array<array<double,5>, 8> | dimensionless | HaighMallionResult |
| hm_shielding_contribution | SphericalTensor | ppm | HaighMallionResult | total shielding from HM |
| n_rings_within_3A | int | - | BiotSavartResult |
| n_rings_within_5A | int | - | BiotSavartResult |
| n_rings_within_8A | int | - | BiotSavartResult |
| n_rings_within_12A | int | - | BiotSavartResult |
| mean_ring_distance | double | Angstroms | BiotSavartResult |
| nearest_ring_atom_distance | double | Angstroms | BiotSavartResult |
| G_iso_exp_sum | double | dimensionless | BiotSavartResult |
| G_T2_exp_sum | array<double, 5> | dimensionless | BiotSavartResult |
| G_iso_var_8A | double | dimensionless | BiotSavartResult |
| bs_shielding_contribution | SphericalTensor | ppm | BiotSavartResult | total shielding from all rings (intensity-weighted sum of G_spherical over all rings) |

Bond anisotropy totals (populated by McConnellResult):
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| mcconnell_co_sum | double | Angstrom^-3 | McConnellResult |
| mcconnell_cn_sum | double | Angstrom^-3 | McConnellResult |
| mcconnell_sidechain_sum | double | Angstrom^-3 | McConnellResult |
| mcconnell_aromatic_sum | double | Angstrom^-3 | McConnellResult |
| mcconnell_co_nearest | double | Angstrom^-3 | McConnellResult |
| nearest_CO_midpoint | Vec3 | Angstroms | McConnellResult |
| nearest_CO_dist | double | Angstroms | McConnellResult |
| nearest_CN_dist | double | Angstroms | McConnellResult |
| T2_CO_nearest | SphericalTensor | Angstrom^-3 | McConnellResult |
| T2_CN_nearest | SphericalTensor | Angstrom^-3 | McConnellResult |
| T2_backbone_total | SphericalTensor | Angstrom^-3 | McConnellResult |
| T2_sidechain_total | SphericalTensor | Angstrom^-3 | McConnellResult |
| T2_aromatic_total | SphericalTensor | Angstrom^-3 | McConnellResult |
| dir_nearest_CO | Vec3 | normalised | McConnellResult |
| mc_shielding_contribution | SphericalTensor | Angstrom^-3 | McConnellResult | total shielding from bond anisotropy — sum of FULL McConnell tensors M_ab/r^3 (asymmetric, non-traceless, T0+T1+T2). See GEOMETRIC_KERNEL_CATALOGUE.md for derivation. |

Coulomb field totals (populated by CoulombResult):
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| coulomb_E_total | Vec3 | V/Angstrom | CoulombResult |
| coulomb_E_backbone | Vec3 | V/Angstrom | CoulombResult |
| coulomb_E_sidechain | Vec3 | V/Angstrom | CoulombResult |
| coulomb_E_aromatic | Vec3 | V/Angstrom | CoulombResult |
| coulomb_EFG_total | Mat3 | V/Angstrom^2 | CoulombResult |
| coulomb_EFG_total_spherical | SphericalTensor | V/Angstrom^2 | CoulombResult |
| coulomb_EFG_backbone | Mat3 | V/Angstrom^2 | CoulombResult |
| coulomb_EFG_backbone_spherical | SphericalTensor | V/Angstrom^2 | CoulombResult |
| coulomb_EFG_aromatic | Mat3 | V/Angstrom^2 | CoulombResult |
| coulomb_EFG_aromatic_spherical | SphericalTensor | V/Angstrom^2 | CoulombResult |
| coulomb_E_solvent | Vec3 | V/Angstrom | CoulombResult (= apbs - vacuum) |
| coulomb_EFG_solvent | Mat3 | V/Angstrom^2 | CoulombResult (= apbs - vacuum) |
| coulomb_E_magnitude | double | V/Angstrom | CoulombResult |
| coulomb_E_bond_proj | double | V/Angstrom | CoulombResult |
| coulomb_E_backbone_frac | double | dimensionless | CoulombResult |
| aromatic_E_magnitude | double | V/Angstrom | CoulombResult |
| aromatic_E_bond_proj | double | V/Angstrom | CoulombResult |
| aromatic_n_sidechain_atoms | int | - | CoulombResult |
| coulomb_shielding_contribution | SphericalTensor | ppm | CoulombResult | Buckingham shielding from E-field (A*E_z + B*E_z^2 for T0, gamma*EFG for T2) |

Note: coulomb_E_ring_proj (E_total . nearest_ring_normal) is computed
at feature extraction time, not in CoulombResult, because it depends
on ring geometry data from BiotSavartResult.

APBS fields (populated by ApbsFieldResult):
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| apbs_efield | Vec3 | V/Angstrom | ApbsFieldResult |
| apbs_efg | Mat3 | V/Angstrom^2 | ApbsFieldResult |
| apbs_efg_spherical | SphericalTensor | V/Angstrom^2 | ApbsFieldResult |

H-bond properties (populated by HBondResult):
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| hbond_nearest_dist | double | Angstroms | HBondResult |
| hbond_nearest_dir | Vec3 | normalised | HBondResult |
| hbond_nearest_tensor | Mat3 | Angstrom^-3 | HBondResult |
| hbond_nearest_spherical | SphericalTensor | Angstrom^-3 | HBondResult |
| hbond_inv_d3 | double | Angstrom^-3 | HBondResult |
| hbond_is_backbone | bool | - | HBondResult |
| hbond_count_within_3_5A | int | - | HBondResult |
| hbond_is_donor | bool | - | HBondResult |
| hbond_is_acceptor | bool | - | HBondResult |
| hbond_shielding_contribution | SphericalTensor | ppm | HBondResult | H-bond dipolar shielding |

Ring-based shielding contributions (per-atom totals, populated by respective results):
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| piquad_shielding_contribution | SphericalTensor | ppm | PiQuadrupoleResult | quadrupole shielding contribution |
| ringchi_shielding_contribution | SphericalTensor | ppm | RingSusceptibilityResult | ring susceptibility shielding |
| disp_shielding_contribution | SphericalTensor | ppm | DispersionResult | dispersion shielding contribution |

Graph topology (populated by MolecularGraphResult):
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| graph_dist_ring | int | bonds | MolecularGraphResult |
| graph_dist_N | int | bonds | MolecularGraphResult |
| graph_dist_O | int | bonds | MolecularGraphResult |
| eneg_sum_1 | double | Pauling units | MolecularGraphResult |
| eneg_sum_2 | double | Pauling units | MolecularGraphResult |
| n_pi_bonds_3 | int | - | MolecularGraphResult |
| is_conjugated | bool | - | MolecularGraphResult |
| bfs_to_nearest_ring_atom | int | bonds | MolecularGraphResult |
| bfs_decay | double | dimensionless | MolecularGraphResult |

ORCA DFT shielding (populated by OrcaShieldingResult):

Each protein gets its own OrcaShieldingResult on its own conformation.
WT and mutant are separate Proteins. Comparison (delta computation) is
done by MutantProteinConformationComparison, not stored on ConformationAtom.

| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| orca_shielding_total | Mat3 | ppm | OrcaShieldingResult |
| orca_shielding_total_spherical | SphericalTensor | ppm | OrcaShieldingResult |
| orca_shielding_diamagnetic | Mat3 | ppm | OrcaShieldingResult |
| orca_shielding_diamagnetic_spherical | SphericalTensor | ppm | OrcaShieldingResult |
| orca_shielding_paramagnetic | Mat3 | ppm | OrcaShieldingResult |
| orca_shielding_paramagnetic_spherical | SphericalTensor | ppm | OrcaShieldingResult |
| has_orca_shielding | bool | - | OrcaShieldingResult |

Note: diamagnetic and paramagnetic stored separately (not just total).
The dia/para decomposition is physics — the paramagnetic contribution
correlates with ring current effects. Both are exported as features
via WriteFeatures() for the calibration pipeline.

### Query methods
- `IsHydrogen() -> bool`: element == Element::H
- `IsHeavy() -> bool`: element != Element::H
- `IsBackbone() -> bool`: is_backbone enrichment property
- `ParentAtom() -> size_t`: for H atoms, parent_atom_index
- `BondDirection(conformation) -> Vec3`: for H atoms, direction from parent to H

### Copy semantics
- Static properties: always preserved
- Enrichment properties: preserved if what-changed does not affect them
  (e.g., hybridisation preserved unless bonds change)
- Conformation-dependent properties: invalidated on copy with new
  foundational properties. Specifically:
  - Change protonation -> invalidate: charges, APBS, Coulomb, ring types,
    all ConformationResult fields, features, predictions
  - Change charges -> invalidate: APBS, Coulomb, features, predictions
  - Change geometry -> invalidate: ALL dynamic properties
  - Same protein same geometry -> preserve: spatial neighbours, DSSP

---

## Bond

Typed covalent bond between two atoms.

### Static properties (const after construction)

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| atom_index_a | size_t | - | First atom in protein's atom list |
| atom_index_b | size_t | - | Second atom in protein's atom list |
| order | BondOrder enum | - | Single, Double, Triple, Aromatic, Peptide |
| category | BondCategory enum | - | From bond detection |
| is_rotatable | bool | - | From bond detection |

### Geometry methods (take positions, conformation-dependent)
- `Midpoint(positions) -> Vec3`: 0.5 * (A + B), Angstroms
- `Length(positions) -> double`: |B - A|, Angstroms
- `Direction(positions) -> Vec3`: normalised (B - A)

### Query methods
- `IsPeptideBond() -> bool`: order == Peptide
- `IsPeptideCO() -> bool`: category == PeptideCO
- `IsBackbone() -> bool`: category in {PeptideCO, PeptideCN, BackboneOther}
- `IsAromatic() -> bool`: category == Aromatic

### Copy semantics
Value type. Always preserved on copy (topology does not change with
protonation for most bonds). Exception: disulfide bonds may change
if CYS protonation changes to CYX.

---

## Ring (class hierarchy)

Ring types ARE classes with physics properties baked in. Each type provides
const properties derived from the ring's identity, not from a lookup table.

**Why this contract exists (for agents creating new ring types):**
The ring current calculators evaluate every ring the same way: vertices,
center, normal, intensity. Your ring type's specific physics (how much
current, what lobe offset, whether nitrogen is present) is what makes
the Biot-Savart result DIFFERENT for your type vs others. If you don't
provide it through the typed interface, someone will hardcode it as a
string comparison. Later, when we refine intensities from mutant data
or add a new ring type, the typed interface means the refinement
propagates everywhere automatically.

### Base class: Ring

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| atom_indices | vector<size_t> | - | Vertex atoms in protein's atom list |
| type_index | RingTypeIndex enum | - | Specific type identifier |
| parent_residue_index | size_t | - | Into protein's residue list |
| parent_residue_number | int | - | PDB sequence number (display only) |
| fused_partner_index | size_t | - | SIZE_MAX if not fused |

Virtual const properties (overridden by each type class):
| Property | Type | Unit | Description |
|----------|------|------|-------------|
| intensity | double | nA/T | DFT-calibrated effective ring current intensity |
| literature_intensity | double | nA/T | Giessner-Prettre 1969 literature value |
| jb_lobe_offset | double | Angstroms | Johnson-Bovey pi-electron lobe offset |
| nitrogen_count | int | - | Number of nitrogen atoms in ring |
| aromaticity | RingAromaticity enum | - | Full, Reduced, Weak |
| ring_size | int | - | Number of vertex atoms (5, 6, or 9) |
| type_name | string | - | "PHE", "TYR", etc. |

### Ring type class hierarchy

```
Ring (abstract base)
|
+-- SixMemberedRing (abstract, ring_size = 6)
|   |
|   +-- PheBenzeneRing
|   |     intensity = -12.0 nA/T
|   |     literature_intensity = -12.0 nA/T
|   |     jb_lobe_offset = 0.64 A
|   |     nitrogen_count = 0
|   |     aromaticity = Full
|   |     ring_size = 6
|   |     type_index = PheBenzene
|   |
|   +-- TyrPhenolRing
|   |     intensity = -11.28 nA/T
|   |     literature_intensity = -11.28 nA/T
|   |     jb_lobe_offset = 0.64 A
|   |     nitrogen_count = 0
|   |     aromaticity = Full
|   |     ring_size = 6
|   |     type_index = TyrPhenol
|   |
|   +-- TrpBenzeneRing
|         intensity = -12.48 nA/T
|         literature_intensity = -12.48 nA/T
|         jb_lobe_offset = 0.64 A
|         nitrogen_count = 0
|         aromaticity = Full
|         ring_size = 6
|         type_index = TrpBenzene
|
+-- FiveMemberedRing (abstract, ring_size = 5)
|   |
|   +-- TrpPyrroleRing
|   |     intensity = -6.72 nA/T
|   |     literature_intensity = -6.72 nA/T
|   |     jb_lobe_offset = 0.52 A
|   |     nitrogen_count = 1
|   |     aromaticity = Reduced
|   |     ring_size = 5
|   |     type_index = TrpPyrrole
|   |
|   +-- HisImidazoleRing
|   |     intensity = -5.16 nA/T
|   |     literature_intensity = -5.16 nA/T
|   |     jb_lobe_offset = 0.50 A
|   |     nitrogen_count = 2
|   |     aromaticity = Weak
|   |     ring_size = 5
|   |     type_index = HisImidazole
|   |     NOTE: unspecified tautomer. Both ND1 and NE2 may carry H.
|   |     Used when protonation state is unknown or ambiguous.
|   |
|   +-- HidImidazoleRing
|   |     intensity = -5.16 nA/T
|   |     literature_intensity = -5.16 nA/T
|   |     jb_lobe_offset = 0.50 A
|   |     nitrogen_count = 2
|   |     aromaticity = Weak
|   |     ring_size = 5
|   |     type_index = HidImidazole
|   |     NOTE: N-delta protonated tautomer. ND1 carries H, NE2 is
|   |     the unprotonated nitrogen. Sparse training data -- type
|   |     exists in the science even if the training set exercises it
|   |     infrequently. DO NOT OMIT.
|   |
|   +-- HieImidazoleRing
|         intensity = -5.16 nA/T
|         literature_intensity = -5.16 nA/T
|         jb_lobe_offset = 0.50 A
|         nitrogen_count = 2
|         aromaticity = Weak
|         ring_size = 5
|         type_index = HieImidazole
|         NOTE: N-epsilon protonated tautomer. NE2 carries H, ND1 is
|         the unprotonated nitrogen. Sparse training data -- type
|         exists in the science. DO NOT OMIT.
|
+-- FusedRing (abstract)
    |
    +-- IndolePerimeterRing
          intensity = -19.2 nA/T
          literature_intensity = -19.2 nA/T
          jb_lobe_offset = 0.60 A
          nitrogen_count = 1
          aromaticity = Full
          ring_size = 9
          type_index = TrpPerimeter
          NOTE: The 9-atom TRP indole perimeter. This is the physical
          ring current path. TRP5+TRP6 are the sub-rings; TRP9 is
          the perimeter that carries the actual current. v1 did not
          have this type (7 types). Rewrite has 8.
```

Rings accumulate properties over extraction passes. After each full
pass over the ring's atoms, a ring-level property update runs.
If a geometric or field property can be precomputed for the ring
as a whole, it is computed and stored on the ring at the end of that
pass, ready for the next ConformationResult.

It is NOT forbidden to traverse the atoms of a specific ring to
harvest properties for calculation. Ring-specific atomic traversal
is how ring properties are built.

### RingTypeIndex (enum)

```
enum class RingTypeIndex {
    PheBenzene    = 0,
    TyrPhenol     = 1,
    TrpBenzene    = 2,
    TrpPyrrole    = 3,
    TrpPerimeter  = 4,
    HisImidazole  = 5,
    HidImidazole  = 6,
    HieImidazole  = 7,
    Count         = 8
};
```

### Ring geometry (conformation-dependent, on ProteinConformation)

| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| center | Vec3 | Angstroms | GeometryResult |
| normal | Vec3 | normalised | GeometryResult (from SVD of vertex positions) |
| radius | double | Angstroms | GeometryResult |
| vertices | vector<Vec3> | Angstroms | GeometryResult |

### Accumulated ring properties (dynamic, set by ConformationResult post-pass updates)

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| total_B_at_center | Vec3 | dimensionless | BiotSavartResult | Sum of B from all OTHER rings at this ring's center |
| intensity_used | double | nA/T | BiotSavartResult | The ring type's effective intensity |
| total_G_T0_diagnostic | double | dimensionless | BiotSavartResult | Sum of |G_T0| over all nearby atoms |
| mutual_B_from | map<size_t, Vec3> | dimensionless | BiotSavartResult | B-field at this ring's center from ring[key] |

### Query methods
- `IsFused() -> bool`: fused_partner_index != SIZE_MAX
- `Size() -> RingSize`: FiveMembered or SixMembered
- `Aromaticity() -> RingAromaticity`: from type class
- `NitrogenCount() -> int`: from type class
- `TypeIndexAsInt() -> int`: for one-hot encoding
- `ComputeGeometry(positions) -> Ring::Geometry`

### Copy semantics
- Static properties (type, vertices, parent): always preserved
- Ring type: INVALIDATED if protonation changes HIS tautomer
  (HIS -> HID or HIE changes the ring type class)
- Geometry: INVALIDATED if conformation changes
- Accumulated properties: INVALIDATED on any foundational property change

---

## ProteinConformation (typed hierarchy)

A conformation without its protein is meaningless. One geometric
instance. Outside of a ProteinConformation, no geometric computation
is allowed. Holds exactly one set of atomic positions.

A ProteinConformation always holds a valid pointer to its owning
Protein. This pointer is set at construction and never null.
Conformations are never orphaned. Do not add reference counting,
weak pointers, or null checks. The Protein's lifetime always
exceeds its conformations'.

**Positions are const after construction.** New geometry = new
ProteinConformation. Spatial index is never stale.

ALL computed properties live on the ProteinConformation as typed
ConformationResult objects, or on its atoms and rings (stored there
by those result objects). This is where ALL work happens. Every
extractor operates within a ProteinConformation.

**Why this contract exists (for agents creating new conformation types):**
Future agents will attach ConformationResults you don't know about yet.
They need spatial queries, atom access, and result attachment to work
identically on your conformation type as on every other type. Your
metadata is what makes your type WORTH having -- it distinguishes your
conformation from others in the protein's list. The base class machinery
is what makes it USABLE by everyone who comes after you. You are building
for people you haven't met.

### Class hierarchy (self-type-reporting)

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

### Core geometry (const after construction)

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| atom_positions | vector<Vec3> | Angstroms | One position per atom in protein |
| protonation | ProtonationState | - | Protonation state for this conformation |

### ConformationResult storage

Results are named singletons on the ProteinConformation. Accessed by
name, not by iteration. Each accessor returns the real typed object.
Throws if not yet computed.

```
auto& dssp = conformation.Dssp();           // DsspResult&
auto& apbs = conformation.ApbsField();      // ApbsFieldResult&
auto& bs = conformation.BiotSavart();       // BiotSavartResult&
auto& coulomb = conformation.Coulomb();     // CoulombResult&
```

Results have physics query methods, not raw data getters:

```
bs.SumT0ByRingType(atomIdx, PheBenzene)
hm.NearestRingContribution(atomIdx)
coulomb.BackboneField(atomIdx)
coulomb.SidechainField(atomIdx)
```

### Result access mechanism

Results are accessed via C++ templates. Adding a new ConformationResult
does NOT require modifying ProteinConformation.

    auto& bs = conformation.Result<BiotSavartResult>();
    if (conformation.HasResult<MopacResult>()) { ... }
    for (auto& [type_id, result] : conformation.AllResults()) { ... }

The named accessors below are convenience wrappers, not the mechanism:

### Named result accessors

| Accessor | Returns | Description |
|----------|---------|-------------|
| `Dssp()` | `DsspResult&` | Secondary structure, phi/psi, SASA, H-bonds |
| `ChargeAssignment()` | `ChargeAssignmentResult&` | Partial charges, VdW radii |
| `SpatialIndex()` | `SpatialIndexResult&` | KD-trees, neighbour lists |
| `GeometryData()` | `GeometryResult&` | Ring/bond geometry, collections |
| `Enrichment()` | `EnrichmentResult&` | Atom roles, categoricals |
| `ApbsField()` | `ApbsFieldResult&` | Solvated E-field and EFG |
| `MolecularGraph()` | `MolecularGraphResult&` | BFS topology, graph features |
| `BiotSavart()` | `BiotSavartResult&` | Ring current geometric kernels |
| `HaighMallion()` | `HaighMallionResult&` | Surface integral tensors |
| `McConnell()` | `McConnellResult&` | Bond anisotropy tensors |
| `Coulomb()` | `CoulombResult&` | Coulomb E-field and EFG |
| `HBond()` | `HBondResult&` | H-bond geometry and tensors |
| `Dispersion()` | `DispersionResult&` | London dispersion tensors |
| `PiQuadrupole()` | `PiQuadrupoleResult&` | Quadrupole EFG |
| `RingSusceptibility()` | `RingSusceptibilityResult&` | Ring susceptibility anisotropy |
| `OrcaShielding()` | `OrcaShieldingResult&` | DFT shielding tensors |
| `Result<MopacResult>()` | `MopacResult&` | MOPAC charges, orbital pops, bond orders |
| `Result<MopacCoulombResult>()` | `MopacCoulombResult&` | Coulomb EFG from QM charges |
| `Result<MopacMcConnellResult>()` | `MopacMcConnellResult&` | Bond-order-weighted anisotropy |

### Result attachment

```
conformation.AttachResult(unique_ptr<ConformationResult> result);
```

At attach time:
1. Dependencies are checked. Missing dependency = immediate logged error.
2. The result computes its values and stores properties on atoms/rings.
3. Once attached, permanent. Nothing is removed.
4. Results accumulate. The next result type finds prior data already there.

### Derived geometry (cached by GeometryResult)

| Property | Type | Unit | Source |
|----------|------|------|--------|
| ring_geometries | vector<Ring::Geometry> | Angstroms | GeometryResult |
| bond_lengths | vector<double> | Angstroms | GeometryResult |
| bond_directions | vector<Vec3> | normalised | GeometryResult |
| bond_midpoints | vector<Vec3> | Angstroms | GeometryResult |
| bounding_min | Vec3 | Angstroms | GeometryResult |
| bounding_max | Vec3 | Angstroms | GeometryResult |
| center_of_geometry | Vec3 | Angstroms | GeometryResult |
| radius_of_gyration | double | Angstroms | GeometryResult |

### Spatial indices (built by SpatialIndexResult)

| Property | Type | Source |
|----------|------|--------|
| atom_kd_tree | nanoflann::KDTreeSingleIndexAdaptor | SpatialIndexResult |
| ring_center_kd_tree | nanoflann::KDTreeSingleIndexAdaptor | SpatialIndexResult |
| bond_midpoint_kd_tree | nanoflann::KDTreeSingleIndexAdaptor | SpatialIndexResult |

### Pre-built collections (built by GeometryResult)

| Property | Type | Source |
|----------|------|--------|
| atoms_by_role | map<AtomRole, vector<size_t>> | EnrichmentResult |
| rings_by_type | map<RingTypeIndex, vector<size_t>> | GeometryResult |
| bonds_by_category | map<BondCategory, vector<size_t>> | GeometryResult |
| residues_by_type | map<AminoAcid, vector<size_t>> | GeometryResult |

### Ring pair properties (built by GeometryResult)

| Property | Type | Source |
|----------|------|--------|
| ring_pairs | vector<RingPair> | GeometryResult |

RingPair:
| Field | Type | Unit |
|-------|------|------|
| ring_a | size_t | ring index |
| ring_b | size_t | ring index |
| center_distance | double | Angstroms |
| normal_dot | double | dimensionless |
| normal_cross_mag | double | dimensionless |
| is_fused | bool | - |

### Seven query patterns (Constitution requirement)

1. **Nearest N rings by distance from point**
   `NearestRings(Vec3 point, int n) -> vector<RingQueryResult>`
   Returns: ring objects sorted by distance, each carrying geometry,
   type properties, direction and distance FROM query point, cylindrical
   coordinates of point in ring frame.

2. **All rings filtered by type**
   `RingsByType(RingTypeIndex type) -> const vector<size_t>&`
   Returns: pre-built collection.

3. **Ring pairs by mutual geometry**
   `RingPairs() -> const vector<RingPair>&`
   Returns: pre-built pairs.

4. **Atoms filtered by role near point**
   `AtomsByRole(AtomRole role, Vec3 point, double radius) -> vector<AtomQueryResult>`
   Returns: intersection of role collection and spatial query.

5. **Bonds filtered by category near point**
   `BondsByCategory(BondCategory cat, Vec3 point, double radius) -> vector<BondQueryResult>`
   Returns: intersection of category collection and spatial query.

6. **Graph neighbours by bond count**
   `GraphNeighbours(size_t atom, int max_bonds) -> vector<GraphNeighbour>`
   Returns: atoms reachable within N bonds, with bond count and path.

7. **Atoms filtered by element within radius**
   `AtomsByElement(Element elem, Vec3 point, double radius) -> vector<AtomQueryResult>`

### Query result types

RingQueryResult:
| Field | Type | Unit | Description |
|-------|------|------|-------------|
| ring_index | size_t | - | Into protein's ring list |
| ring | const Ring& | - | Full ring object |
| geometry | const Ring::Geometry& | - | This conformation's geometry |
| distance | double | Angstroms | From query point to ring center |
| direction | Vec3 | normalised | From query point to ring center |
| rho | double | Angstroms | Cylindrical radial distance in ring frame |
| z | double | Angstroms | Signed height above ring plane |
| theta | double | radians | Angle from ring normal |
| ring_neighbourhood | const RingNeighbourhood* | - | If accumulated, pointer to stored data |

AtomQueryResult:
| Field | Type | Unit | Description |
|-------|------|------|-------------|
| atom_index | size_t | - | Into protein's atom list |
| atom | const Atom& | - | Full atom object |
| distance | double | Angstroms | From query point |
| direction | Vec3 | normalised | From query point |

BondQueryResult:
| Field | Type | Unit | Description |
|-------|------|------|-------------|
| bond_index | size_t | - | Into protein's bond list |
| bond | const Bond& | - | Full bond object |
| distance | double | Angstroms | From query point to bond midpoint |
| direction | Vec3 | normalised | From query point to bond midpoint |

GraphNeighbour:
| Field | Type | Description |
|-------|------|-------------|
| atom_index | size_t | Target atom |
| bond_count | int | BFS distance |
| path | vector<size_t> | Bond indices along path |

### Copy policies

**GeometryOnly**: copy positions, clear all ConformationResult objects.
Use when foundational properties (protonation, charges) have changed
and all computed results must be recomputed.

**Full**: copy everything including all attached ConformationResult
objects. Only valid if ProteinBuildContext is unchanged. Use for
creating working copies where results are still valid.

### Copy semantics
- Atom positions: always preserved (these ARE the conformation)
- Bond connectivity / neighbour lists / distances: preserved (purely geometric)
- DSSP: preserved (backbone geometry unchanged)
- Partial charges: INVALIDATED if protonation changes
- APBS fields: INVALIDATED if charges change
- All ConformationResult objects: INVALIDATED if charges or ring types change
- Features: INVALIDATED if any ConformationResult changes
- Predictions: INVALIDATED if features change

---

## ConformationResult (ABC)

Base class for all computed results that attach to a ProteinConformation.
Each result type is a named singleton. Accessed by name, checked at
attach time.

```
class ConformationResult {
public:
    virtual ~ConformationResult() = default;
    virtual string Name() const = 0;
    virtual vector<type_index> Dependencies() const = 0;
};
```

Note: the ABC does NOT have an Attach() method. Computation happens
in the static Compute() factory on each subclass. The framework calls
AttachResult() on the ProteinConformation, which checks dependencies
via the type_index vector. Dependencies are compile-time type-safe,
not string-matched.

**WHY this contract exists:** Results accumulate on a ProteinConformation
over the extraction pipeline. Each result declares what it needs (other
results that must already be attached) and what it provides (properties
on atoms/rings). The dependency graph replaces manual phase ordering.
The attach method is where computation happens: the result reads from
prior results and stores its own properties. Once attached, permanent.

### Known result types (with dependencies)

```
GeometryResult                requires: nothing
DsspResult                    requires: nothing
ChargeAssignmentResult        requires: nothing
MopacResult                   requires: nothing (conformation electronic structure,
                                runs early as precondition for MOPAC-derived calculators)
EnrichmentResult              requires: nothing
SpatialIndexResult            requires: GeometryResult
ApbsFieldResult               requires: ChargeAssignmentResult
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
OrcaShieldingResult           requires: nothing (loaded from files)
SasaResult                    requires: SpatialIndexResult
AIMNet2Result                 requires: SpatialIndexResult (model loaded externally)
WaterFieldResult              requires: SpatialIndexResult + SolventEnvironment (runtime)
HydrationShellResult          requires: SpatialIndexResult + SolventEnvironment (runtime)
GromacsEnergyResult           requires: nothing (reads from .edr file at runtime)

(Feature extraction is distributed: each ConformationResult implements
WriteFeatures(). ConformationResult::WriteAllFeatures() traverses all
attached results and exports NPY arrays for the calibration pipeline.)
```

### What each result stores on atoms/rings

Each result type, when attached, populates properties on the atoms
and rings of the ProteinConformation. The next result type in the
resolved order finds these properties already present.

See the per-result sections below and the Atom field accumulation
tables above for exactly what each result stores.

---

## Trajectory-scope entities

**Current-truth note:** A library-reference restatement of this
material appears at the end of this file as
`# Addition — Trajectory-scope entities (2026-04-24, library-reference form)`.
The addition reflects the code in the current tree; this section is
retained for design history. Where the two disagree, the addition is
authoritative.

Trajectory-scope classes are modular accumulators on top of the
per-frame ConformationResult pipeline — not replacements for it.
The static-analysis path (Protein + ProteinConformation +
ConformationResult + OperationRunner::Run) is unchanged and
authoritative; trajectory-scope types extend the object model to
cover multi-frame runs without changing any of the above.

The organizing discipline (typed objects, private constructors,
one writer per field, dependencies declared via type_index,
self-serialising results) transfers directly from conformation
scope. The structural differences are multi-phase compute (per-frame
+ end-of-stream), internal accumulator state on Results rather than
on the per-atom store, and richer serialisation patterns (dense
buffers for time series, selection bags for event streams). See
PATTERNS.md §§13-18 for the disciplines; see below for the concrete
types.

### TrajectoryProtein

A protein in trajectory context. Wraps a `Protein` (identity +
topology), adds the per-atom trajectory-scope running buffer and the
attached TrajectoryResults. Replaces the old `GromacsProtein` role
per spec/pending_include_trajectory_scope_2026-04-22.md §3.

| Property | Type | Description |
|----------|------|-------------|
| `protein_` | `unique_ptr<Protein>` | Owned Protein (identity + topology) |
| `charges_` | `unique_ptr<ChargeSource>` | Owned charge source from TPR |
| `sys_reader_` | `FullSystemReader` | TPR topology handler, borrowed by GromacsFrameHandler |
| `bonded_params_` | `BondedParameters` | From TPR parse, consumed via RunOptions.bonded_params per frame |
| `atoms_` | `vector<TrajectoryAtom>` | Per-atom running buffer, parallel to Protein.atoms |
| `results_` | `unordered_map<type_index, unique_ptr<TrajectoryResult>>` | Attached Results, singleton-per-type |
| `results_attach_order_` | `vector<TrajectoryResult*>` | Attach order = dispatch order |
| `dense_buffers_` | `unordered_map<type_index, unique_ptr<DenseBufferBase>>` | Per-Result-owner dense buffers transferred at Finalize |
| `finalized_` | `bool` | True after FinalizeAllResults |

Key operations:

- `BuildFromTrajectory(dir)` — parses md.tpr, builds Protein + charges.
  Does NOT finalize (bond detection needs first-frame geometry).
- `Seed(positions, time_ps)` — runs `Protein::FinalizeConstruction`,
  adds canonical conformation via `AddMDFrame`, allocates
  `TrajectoryAtoms`. Called once by `Trajectory::Run` Phase 2. Must
  precede TrajectoryResult attach.
- `CanonicalConformation()` — returns conf0; same object
  `Protein::Conformation()` returns in static paths. Permanent on
  `Protein.conformations_` with its ConformationResults.
- `TickConformation(positions)` — creates an ephemeral
  `ProteinConformation` pointing at the wrapped Protein for topology;
  lifetime is the caller's iteration (frames 1..N).
- `AttachResult(unique_ptr<TrajectoryResult>)` — singleton-per-type
  check, appends to results_attach_order_. Dependency validation is
  at Trajectory::Run Phase 4, not here.
- `Result<T>()` / `HasResult<T>()` / `AllResults()` / `ResultsInAttachOrder()` — typed access.
- `DispatchCompute(conf, traj, frame_idx, time_ps)` — iterates
  attached Results, calls each one's Compute.
- `FinalizeAllResults(traj)` — iterates, calls Finalize, sets
  finalized_.
- `AdoptDenseBuffer<T>(buffer, owner_type)` /
  `GetDenseBuffer<T>(owner_type)` — Result ownership transfer for
  per-atom × stride dense data.
- `WriteH5(file)` / `WriteFeatures(dir)` — traverses attached
  Results, each emits its own group / NPY files.

Identity flows through the wrapped Protein:

```cpp
const TrajectoryAtom& ta = tp.AtomAt(42);
const Atom& identity = tp.ProteinRef().AtomAt(42);  // element, bonds, residue
```

### TrajectoryAtom

Per-atom trajectory-scope data store. Private constructor (only
`TrajectoryProtein` constructs via friend access); one instance per
Protein atom, constructed once by `TrajectoryProtein::Seed` and
never resized.

Three coexisting field shapes (see PATTERNS.md §13):

**Typed accumulator fields.** One writer per field. Current fields
are BS shielding Welford state written by
`BsWelfordTrajectoryResult`:

| Field | Type | Writer | Description |
|-------|------|--------|-------------|
| `bs_t0_mean`, `bs_t0_m2`, `bs_t0_std`, `bs_t0_min`, `bs_t0_max` | `double` | BsWelford | Welford state for BS T0 shielding; std is Finalize-only |
| `bs_t0_min_frame`, `bs_t0_max_frame` | `size_t` | BsWelford | Frame indices of extrema |
| `bs_n_frames` | `size_t` | BsWelford | Samples accumulated; used as denominator |
| `bs_t2mag_{mean,m2,std,min,max,min_frame,max_frame}` | `double` / `size_t` | BsWelford | Same for \|T2\| magnitude |
| `bs_t0_delta_{mean,m2,std,min,max,n}` | `double` / `size_t` | BsWelford | Frame-to-frame T0 delta |

Additional fields are added as `*Welford`-style TrajectoryResults
land.

**Typed struct vectors for known-shape per-source data.** Aspirational;
no current fields. Parallel to `ConformationAtom::ring_neighbours` but
holding trajectory-averaged stats (e.g. future
`RingNeighbourhoodTrajectoryStats` — per atom, per nearby ring, the
trajectory-averaged geometry + per-calculator tensor statistics).

**Per-atom event bag.** `RecordBag<AtomEvent> events` — open-shape
event bag for per-atom events emitted by scan-mode detectors and
lifetime/transition accumulators (rotamer transitions, H-bond
form/break, anomaly markers).

No identity duplication: element, residue, bonds, ring membership
come from the wrapped Protein via `tp.ProteinRef().AtomAt(i)`.

### TrajectoryResult (ABC)

Base class for per-trajectory modular calculators. Parallel to
`ConformationResult` at conformation scope but with multi-phase
compute.

Virtual methods:

| Method | Signature | Role |
|--------|-----------|------|
| `Name()` | `-> string` | Human-readable name for logs / attach diagnostics |
| `Dependencies()` | `-> vector<type_index>` | TrajectoryResult OR ConformationResult types that must be attached/running. Validated at Trajectory::Run Phase 4 |
| `Compute(conf, tp, traj, frame_idx, time_ps)` | `-> void` | Per-frame work. Reads conf and prior TrajectoryAtom state; writes owned fields; pushes to atom events or run-scope selections |
| `Finalize(tp, traj)` | `-> void` | End-of-stream. Converts accumulator state (M2 → std); transfers DenseBuffer ownership to tp; derives final fields |
| `WriteFeatures(tp, output_dir)` | `-> int` | Optional. NPY emission |
| `WriteH5Group(tp, file)` | `-> void` | Optional. H5 group emission |

`Compute` is virtual (dispatched per-frame per-Result); `Name` and
`Dependencies` are pure virtual; the rest are default-implementable.
Results construct via static `Create(tp)` factories that return
unique_ptr; the factory allocates internal per-atom buffers sized
from `tp.AtomCount()` or similar (PATTERNS.md §15: Seed precedes
Attach, so factories see a finalized Protein).

### Known TrajectoryResult types (as of landing)

| Type | Scope | Dependencies | Lifecycle | Output |
|------|-------|--------------|-----------|--------|
| `BsWelfordTrajectoryResult` | per-atom | `BiotSavartResult` | AV | TrajectoryAtom fields (bs_t0_\*, bs_t2mag_\*, bs_t0_delta_\*) + `/trajectory/bs_welford/` |
| `BsShieldingTimeSeriesTrajectoryResult` | per-atom | `BiotSavartResult` | FO | `DenseBuffer<SphericalTensor>` + `/trajectory/bs_shielding_time_series/` (N, T, 9) with irrep_layout / normalization / parity attrs |
| `BsAnomalousAtomMarkerTrajectoryResult` | per-atom | `BsWelfordTrajectoryResult`, `BiotSavartResult` | AV | per-atom events bag with kinds `BsAnomalyHighT0`, `BsAnomalyLowT0` |
| `BsT0AutocorrelationTrajectoryResult` | per-atom × lag | `BiotSavartResult` | FO | `DenseBuffer<double>` + `/trajectory/bs_t0_autocorrelation/` (N, N_LAGS=120) with estimator / mean_convention attrs |
| `BondLengthStatsTrajectoryResult` | per-bond | (none) | AV | internal `vector<PerBondWelford>` + `/trajectory/bond_length_stats/` |
| `PositionsTimeSeriesTrajectoryResult` | per-atom | (none) | FO | `DenseBuffer<Vec3>` + `/trajectory/positions/` (N, T, 3) |
| `ChiRotamerSelectionTrajectoryResult` | per-residue (emits) | (none) | AV | `traj.MutableSelections()` push per transition detected |

Lifecycle: **AV** = always-valid mid-stream (Compute updates
output fields in place each frame); **FO** = Finalize-only (Compute
appends to internal buffers; output materialises at Finalize).

Additional types are spec'd in
`spec/pending_include_trajectory_scope_2026-04-22.md` Appendix F;
they land as the time-series, Welford, and scan-mode families fill in.

### Trajectory

Process entity representing one traversal of an XTC + TPR + EDR
source. Holds both process-role state (handler, env) during Run and
record-role state (frame times, frame indices, selection bag, source
paths) after.

| Property | Type | Description |
|----------|------|-------------|
| `xtc_path_`, `tpr_path_`, `edr_path_` | `filesystem::path` | Source paths |
| `edr_frames_` | `vector<GromacsEnergy>` | Preloaded at construction for O(log T) time lookup |
| `env_` | `TrajectoryEnv` | Single-slot per-frame environment (solvent, current_energy pointer, frame idx, time) — populated by handler / Trajectory::Run each frame |
| `frame_count_`, `frame_times_`, `frame_indices_` | `size_t` / vectors | Growing record during Run |
| `selections_` | `RecordBag<SelectionRecord>` | Run-scope event bag; Results push directly |
| `output_dir_` | `filesystem::path` | Record field for downstream writers |
| `handler_` | `unique_ptr<GromacsFrameHandler>` | During-Run only |
| `state_` | `enum State` | Constructed / Running / Complete |

`Run(tp, config, session, extras, output_dir)` drives the traversal
in 8 named phases:

1. Open handler (mount XTC + build PBC fixer).
2. Read frame 0 + `tp.Seed` (finalize Protein + create conf0 +
   init TrajectoryAtoms).
3. Attach TrajectoryResults from `config.TrajectoryResultFactories()`
   + caller extras.
4. Validate dependencies and caller-supplied resources (e.g.
   `session.Aimnet2Model()` if `config.RequiresAimnet2()`).
5. Build per-frame `RunOptions` template from config + tp + session.
6. Frame 0: populate env, run OperationRunner on canonical conf,
   tp.DispatchCompute, record.
7. Per-frame loop: skip stride-1, read next, tickify, update env,
   OperationRunner, tp.DispatchCompute, record.
8. `tp.FinalizeAllResults(traj)`. Selections were pushed directly
   during Compute/Finalize by the Results themselves; no collection
   sweep here.

Returns `Status` (`errors.h`); `kOk` on success.

`WriteH5(file)` emits the `/trajectory/` group with source paths,
frame metadata, and one selection sub-group per kind present in the
bag. The per-atom and per-Result outputs come from
`TrajectoryProtein::WriteH5`.

### Session

Process-wide resources wrapper. Owns the AIMNet2 model lifetime and
orchestrates `RuntimeEnvironment::Load`, `OperationLog::LoadChannelConfig`,
and session-start logging.

| Property | Type | Description |
|----------|------|-------------|
| `aimnet2_model_` | `unique_ptr<AIMNet2Model>` | Loaded once, passed to `Trajectory::Run` via the config's `RequiresAimnet2` flag + session access |
| `last_error_` | `string` | Error from the most recent failed Load call |

Methods: `LoadFromToml()`, `LoadAimnet2Model(path)`, `Aimnet2Model()`,
`HasAimnet2Model()`, `LastError()`. Returns `Status` from load
methods.

### RunConfiguration

Typed description of a trajectory-run shape. Three named static
factories (matching Use Case E/F shapes + the DFT scan loop):

- **`ScanForDftPointSet()`** — *Slated for removal; a later doc pass
  will handle the replacement.* Cheap per-frame set (no MOPAC, no
  APBS, no Coulomb; Geometry / SpatialIndex / Enrichment / DSSP /
  BiotSavart / SASA). Attaches `BsWelfordTrajectoryResult`
  (placeholder; scan-mode emitters + DftPoseCoordinator land with the
  scan family). For choosing DFT pose frames from MD.
- **`PerFrameExtractionSet()`** — production canonical. Full
  classical stack every frame, MOPAC skipped (sparse-frame only),
  Coulomb skipped (APBS supersedes at N > 1000 atoms), AIMNet2
  required. Stride 2 by default (25 ns × 1250 frames → 625 sampled).
  Attaches the current exemplar Results (BsWelford, BsShieldingTimeSeries,
  BsAnomalousAtomMarker, BsT0Autocorrelation, BondLengthStats,
  PositionsTimeSeries). Expands as the catalog fills.
- **`FullFatFrameExtraction()`** — PerFrameExtractionSet + MOPAC on
  a selected frame subset (DFT pose set, μs harvester checkpoints).

Each configuration records: per-frame `RunOptions` base (with skip
flags), required `ConformationResult` types (for Phase 4 dependency
validation), `TrajectoryResult` factory list (attach order = dispatch
order), stride, and a `RequiresAimnet2` flag.

### RecordBag<Record>

Template container for typed event streams. Used at two scopes:

- `RecordBag<SelectionRecord>` on `Trajectory` — run-scope
  frame-level events (rotamer transitions, pose candidates, RMSD
  spikes).
- `RecordBag<AtomEvent>` on each `TrajectoryAtom` — atom-scope
  events attributed to a specific atom (anomaly flags, H-bond
  form/break, ring flips).

The record type must expose `kind` (type_index), `frame_idx`
(size_t), `time_ps` (double) as explicit top-level fields —
windowed queries are direct reads, not string lookups.

Primitives: `Push`, `All`, `Count`, `Kinds`. Affordance queries:
`ByKind<T>()`, `ByKindSinceFrame<T>(from)`,
`ByKindSinceTime<T>(since_ps)`, `MostRecent<T>()`, `CountByKind<T>()`.

```cpp
struct SelectionRecord {
    std::type_index kind;
    std::size_t frame_idx = 0;
    double time_ps = 0.0;
    std::string reason;
    std::map<std::string, std::string> metadata;
};

struct AtomEvent {
    std::type_index emitter;         // which Result pushed
    std::type_index kind;            // emitter-specific discriminator
    std::size_t frame_idx = 0;
    double time_ps = 0.0;
    std::map<std::string, std::string> metadata;
};
```

`AtomEvent` carries a separate `emitter` alongside `kind` because
one emitter typically produces several semantic kinds at atom scope
(`BsAnomalousAtomMarker` emits both `BsAnomalyHighT0` and
`BsAnomalyLowT0`).

### DenseBuffer<T>

Template for contiguous per-atom × stride storage of Results' FO
outputs. Owned by `TrajectoryProtein` after a Result transfers
ownership at Finalize, keyed by the owning Result's `type_index`.
Atom-major layout: per-atom slice is contiguous.

Typical payloads and emission shapes:

| T | sizeof(T) | Example Result | H5 shape |
|---|-----------|----------------|----------|
| `double` | 8 B | `BsT0Autocorrelation` | `(N, stride)` |
| `Vec3` | 24 B | `PositionsTimeSeries` | `(N, stride, 3)` |
| `SphericalTensor` | 72 B | `BsShieldingTimeSeries` | `(N, stride, 9)` |
| `Mat3` | 72 B | future `*EFGTimeSeries` | `(N, stride, 3, 3)` |

H5 emission builds an explicit flat `vector<double>` via named
component access — `.x()/.y()/.z()`, `.T0/.T1[k]/.T2[k]`,
`m(i,j)` — and writes with a typed `DataSpace`. No reinterpret_cast
on the raw buffer; no assumption about Eigen / struct packing
layout. Every future `DenseBuffer<T>` emitter follows this pattern.

See PATTERNS.md §§13-18 for the disciplines; for the catalog of
pending TrajectoryResults + the deferred H5 metadata schema see the
Extended section at the end of this file.

---

## GeometryChoice (TBD — aspirational, next after calculator tuning)

Every calculator makes implicit geometric choices: distance cutoffs,
near-field exclusion zones, sequential separation thresholds, ring
current intensities, lobe offsets. A GeometryChoice is a runtime record
of one such decision — what was decided, on which model objects, with
what value, for what physics reason. Created by calculators, deposited
on the conformation, drawable in the viewer.

See spec/GEOMETRY_CHOICE_BRIEF.md for the full design and per-calculator
manifest. Prior CalculationArea type vocabulary is in learn/bones/.

---

## Calculator Shielding Contribution Contract

Every classical calculator ConformationResult (BiotSavartResult,
HaighMallionResult, McConnellResult, CoulombResult, PiQuadrupoleResult,
RingSusceptibilityResult, DispersionResult, HBondResult,
MopacCoulombResult, MopacMcConnellResult) must produce
TWO kinds of output:

### 1. Geometric output (for features and reuse)
The raw geometric kernel or field in the calculator's natural units
(dimensionless for ring current kernels, Angstrom^-3 for dipolar tensors,
V/Angstrom for E-fields, etc.). Stored per atom per source (ring/bond).
This is what the feature extractor reads. It can be reused with different
parameters without recomputing geometry.

### 2. Shielding contribution (for calibration and upstream models, in ppm)
The total shielding contribution from this calculator, per atom, as a
SphericalTensor in ppm. This is the geometric output multiplied by the
appropriate parameter (intensity, anisotropy, Buckingham coefficient, etc.).
Both representations are exported via WriteFeatures() as NPY arrays — the
calibration pipeline uses them to tune parameters against DFT deltas, and
the upstream e3nn prediction model consumes them as physics-grounded
tensor features at all irrep levels (T0, T1, T2):

- BiotSavart: sigma_atom = Sum_rings[ I_type * G_tensor ]
- HaighMallion: sigma_atom = Sum_rings[ J_type * HM_tensor ]
- McConnell: sigma_atom = Sum_bonds[ Dchi_cat * dipolar_tensor ]
- Coulomb: sigma_atom = A_elem * E_z + B_elem * E_z^2 (T0) + gamma_elem * EFG (T2)
- PiQuadrupole: sigma_atom = Sum_rings[ Q_type * quad_tensor ]
- RingSusceptibility: sigma_atom = Sum_rings[ Dchi_ring_type * chi_tensor ]
- Dispersion: sigma_atom = Sum_rings[ alpha_disp_elem * disp_tensor ]
- HBond: sigma_atom = eta_elem * hbond_tensor

### What this enables

The shielding contribution is stored per atom as:
| Property | Type | Unit | Source result |
|----------|------|------|---------------|
| shielding_contribution | SphericalTensor | ppm | each calculator |

This is SEPARATE from the geometric output. Both are stored. The
geometric output feeds features; the shielding contribution provides
the parameterised tensor for calibration comparison and for the
upstream e3nn prediction model that consumes full tensor features.

After all calculators attach, the T2 residual between classical
kernel output and DFT delta tensors shows where the angular physics
fails. This comparison happens in the calibration pipeline (Python),
which reads both sets of NPY arrays.

### Per-calculator shielding contribution

| Calculator | Formula | Parameter source |
|------------|---------|-----------------|
| BiotSavart | I_type * G_T0/T1/T2 | Ring type defaults or TOML-calibrated values |
| HaighMallion | J_type * HM_T0/T1/T2 | Ring type defaults or TOML-calibrated values |
| McConnell | Dchi_cat * dipolar_T0/T1/T2 | Bond category defaults or TOML-calibrated values |
| Coulomb | A*E + B*E^2 (T0), gamma*EFG (T2) | Element defaults or TOML-calibrated values |
| PiQuadrupole | Q_type * quad_T0/T1/T2 | Ring type defaults or TOML-calibrated values |
| RingSusceptibility | Dchi_ring * chi_T0/T1/T2 | Ring type defaults or TOML-calibrated values |
| Dispersion | alpha_elem * disp_T0/T1/T2 | Element defaults or TOML-calibrated values |
| HBond | eta_elem * hbond_T0/T1/T2 | Element defaults or TOML-calibrated values |

---

### Unit Test Contract for Classical Calculators

Every classical calculator must have unit tests that verify:

1. **Full tensor output with DEFAULT parameters against known analytical values.**
   Not just T0. Specific T2 components must match analytical predictions
   for a known test geometry (e.g., proton at specific position relative
   to a PHE ring). This proves the physics is correct independent of
   any model.

2. **Shielding contribution in ppm.** The intensity-weighted / parameter-weighted
   total must match the analytical prediction.

3. **Two-path consistency.** Default and corrected paths produce the same
   output when corrected parameters equal default parameters.

A test that only checks T0 allows the agent to skip T2 computation.
A test that checks specific T2[0..4] values forces full tensor implementation.

---

## Per-Result Minimum Output (Constitution contract)

**BiotSavartResult:**
- Mat3: full geometric kernel G_ab per atom per ring
- sphericart: T0, T1[3], T2[5] decomposition of G_ab
- Vec3: magnetic field B at the point per ring
- Per-ring attribution
- Ring-frame cylindrical B-field: B_n, B_rho, B_phi
- Post-pass: total B, per-type sums, ring counts, ring distances

**HaighMallionResult:**
- Mat3: raw surface integral H per atom per ring (symmetric, traceless, A^-1)
- SphericalTensor of H: pure T2 (T0=0, T1=0 by construction)
- Vec3: effective B-field V = H·n per ring (A^-1)
- Full shielding kernel G = -n⊗V accumulated into hm_shielding_contribution
  (rank-1, same structure as BS, sign convention: sigma = I × G)
- Per-ring attribution
- Post-pass: per-type sums

**McConnellResult:**
- Mat3: dipolar tensor per bond per atom
- sphericart: T0, T1[3], T2[5]
- Per-bond-category subtotals (peptide CO, CN, sidechain, aromatic)
  as separate Mat3 + sphericart each
- Scalar factor derived FROM tensor trace, not computed instead

**CoulombResult:**
- Vec3: E-field at each atom
- Mat3: EFG tensor V_ab at each atom
- sphericart: T0, T1[3], T2[5] of V_ab
- Decomposed: backbone vs sidechain vs solvent contributions

**PiQuadrupoleResult:**
- Mat3: quadrupole EFG tensor per atom per ring
- sphericart: T0, T1[3], T2[5]
- Per-ring attribution

**RingSusceptibilityResult:**
- Mat3: dipolar tensor from ring center per atom per ring
- sphericart: T0, T1[3], T2[5]
- Per-ring attribution

**DispersionResult:**
- Mat3: anisotropic dispersion tensor per ring (summed over vertices)
- double: isotropic 1/r^6 contribution
- sphericart: T0, T1[3], T2[5]
- Contact count per ring

**HBondResult:**
- Mat3: dipolar tensor to H-bond partner(s)
- sphericart: T0, T1[3], T2[5]
- double: isotropic cos^2 theta / r^3 contribution
- Geometry: distance, D-H-A angle, donor/acceptor classification

**ApbsFieldResult:**
- Vec3: solvated E-field from APBS
- Mat3: solvated EFG from APBS
- sphericart: T0, T1[3], T2[5] of EFG

**OrcaShieldingResult:**
- Per-atom shielding tensors: diamagnetic, paramagnetic, total
  (Mat3 + sphericart each)
- Element-verified atom ordering (ORCA nucleus element must match protein atom)
- One result per protein per conformation. WT and mutant are separate Proteins.
- Mutant comparison (atom matching, delta computation) is a separate
  static operation (MutantProteinConformationComparison), not part of this result.

---

## RingNeighbourhood

Per-atom, per-ring structured result. One entry per aromatic ring
within range of an atom. Built by ring current ConformationResult objects.

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| ring_index | size_t | - | BiotSavartResult | Into protein's ring list |
| ring_type | RingTypeIndex | - | BiotSavartResult | Type of the source ring |
| distance_to_center | double | Angstroms | BiotSavartResult | |r - center| |
| direction_to_center | Vec3 | normalised | BiotSavartResult | (center - r) / |center - r| |
| rho | double | Angstroms | BiotSavartResult | Cylindrical radial in ring frame |
| z | double | Angstroms | BiotSavartResult | Signed height above ring plane |
| theta | double | radians | BiotSavartResult | atan2(rho, z) |
| G_tensor | Mat3 | dimensionless | BiotSavartResult | BS geometric kernel G=-n⊗B×PPM (rank-1) |
| G_spherical | SphericalTensor | dimensionless | BiotSavartResult | Decomposition of G |
| B_field | Vec3 | Tesla | BiotSavartResult | JB B-field from unit current (1 nA) |
| B_cylindrical | Vec3 | Tesla | BiotSavartResult | (B_rho, 0, B_z) in ring frame |
| hm_tensor | Mat3 | Angstrom^-1 | HaighMallionResult | Raw surface integral H (symmetric, traceless) |
| hm_spherical | SphericalTensor | Angstrom^-1 | HaighMallionResult | Decomposition of H (pure T2) |
| hm_B_field | Vec3 | Angstrom^-1 | HaighMallionResult | Effective B-field V = H·n |
| quad_tensor | Mat3 | Angstrom^-4 | PiQuadrupoleResult | Quadrupole EFG |
| quad_spherical | SphericalTensor | Angstrom^-4 | PiQuadrupoleResult | |
| quad_scalar | double | Angstrom^-4 | PiQuadrupoleResult | (3cos^2 theta - 1)/r^4 |
| chi_tensor | Mat3 | Angstrom^-3 | RingSusceptibilityResult | Ring susceptibility dipolar tensor |
| chi_spherical | SphericalTensor | Angstrom^-3 | RingSusceptibilityResult | |
| chi_scalar | double | Angstrom^-3 | RingSusceptibilityResult | (3cos^2 theta - 1)/r^3 |
| disp_tensor | Mat3 | Angstrom^-6 | DispersionResult | London dispersion tensor |
| disp_spherical | SphericalTensor | Angstrom^-6 | DispersionResult | |
| disp_scalar | double | Angstrom^-6 | DispersionResult | Sum of 1/r^6 over vertices |
| disp_contacts | int | - | DispersionResult | Count of vertex contacts in range |
| gaussian_density | double | dimensionless | BiotSavartResult | Learned per-type spatial envelope |

---

## BondNeighbourhood

Per-atom, per-bond structured result. One entry per bond within range
of an atom. Built by McConnellResult.

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| bond_index | size_t | - | McConnellResult | Into protein's bond list |
| bond_category | BondCategory | - | McConnellResult | Category of source bond |
| distance_to_midpoint | double | Angstroms | McConnellResult | |r - midpoint| |
| direction_to_midpoint | Vec3 | normalised | McConnellResult | (midpoint - r) / d |
| dipolar_tensor | Mat3 | Angstrom^-3 | McConnellResult | Full dipolar tensor |
| dipolar_spherical | SphericalTensor | Angstrom^-3 | McConnellResult | Decomposition |
| mcconnell_scalar | double | Angstrom^-3 | McConnellResult | Derived from tensor trace |

---

## AtomNeighbourhood

Per-atom, per-atom spatial relationship. Built by SpatialIndexResult.

| Property | Type | Unit | Source result | Description |
|----------|------|------|---------------|-------------|
| atom_index | size_t | - | SpatialIndexResult | Into protein's atom list |
| distance | double | Angstroms | SpatialIndexResult | |r_a - r_b| |
| direction | Vec3 | normalised | SpatialIndexResult | (r_b - r_a) / d |

---

## SpatialIndex

Wrapper around nanoflann KD-tree. Three instances on each
ProteinConformation, built by SpatialIndexResult.

### Atom spatial index
- Built from: atom_positions (vector<Vec3>)
- Query: `RadiusSearch(Vec3 point, double radius) -> vector<pair<size_t, double>>`
- Query: `KNearestSearch(Vec3 point, int k) -> vector<pair<size_t, double>>`

### Ring center spatial index
- Built from: ring_geometries[i].center for all rings
- Same query interface

### Bond midpoint spatial index
- Built from: bond_midpoints for all bonds
- Same query interface

All spatial queries through these indices. No linear scans over atoms.
Constitution Rule 4.

---

## Feature (base class)

Each feature is a subclass. Adding a feature = writing a class.

| Property | Type | Description |
|----------|------|-------------|
| name | string | Unique identifier for manifest |
| irrep | IrrepType enum | L0, L1e, L1o, L2e |
| components | int | 1 (L0), 3 (L1), 5 (L2) |

### IrrepType enum
```
enum class IrrepType { L0, L1e, L1o, L2e };
```

### Virtual method
```
virtual void Compute(size_t atom_index,
                     const ProteinConformation& conf,
                     FeatureOutput& out) const = 0;
```

Features ONLY READ from ConformationResult objects on the
ProteinConformation. They do NOT compute physics. They do NOT access
raw positions or charges directly. They read what prior
ConformationResult objects stored.

### FeatureOutput
```
struct FeatureOutput {
    vector<double> L0;          // scalar features
    vector<Vec3> L1;            // vector features
    vector<IrrepType> L1_parity; // e or o per L1 feature
    vector<array<double,5>> L2; // tensor features
};
```

---

## Framework Store Interface

ConformationResult objects store properties on atoms and rings through
the ProteinConformation's typed store interface:

```
conformation.StoreRingContribution(atom_index, ring_index, result);
conformation.StoreBondContribution(atom_index, bond_index, result);
conformation.StoreAtomProperty(atom_index, property_name, value);
```

Each store operation:
- Is typed (cannot store double where Mat3 goes)
- Is logged (UDP log entry emitted automatically)
- Is append-only (writing to a filled slot is WARNING)

Each retrieval:
- Is typed (returns correct type)
- Is checked (empty slot is ERROR)

### Enforcement: The Framework Stores, Not the Extractor

An extractor does not directly write `atom.someProperty = value`.
The ConformationResult, during attachment, stores properties on
atoms and rings through the typed store interface.

---

## Static vs Dynamic Object Properties

### Static: determined by what the object IS
- An atom's element is static.
- A ring's type (PheBenzene, HisImidazole) is static after
  protonation state is applied.
- A residue's amino acid type is static.
- A bond's category is static after bond detection.

Static properties are set at construction or enrichment and are
const thereafter. They travel with copies.

### Dynamic: determined by computation within a ProteinConformation
- An atom's partial charge is dynamic (from ChargeAssignmentResult)
- An atom's ring neighbourhood tensors are dynamic (from
  BiotSavartResult et al.)
- A ring's accumulated field properties are dynamic

Dynamic properties are set by ConformationResult attachments and are
ProteinConformation-specific. They do NOT travel with GeometryOnly
copies. They travel with Full copies only if the ProteinBuildContext
is unchanged.

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

---

## Complete property flow: which ConformationResult sets what

### Construction (PDB loading)
Protein: residues, atoms, rings (topology), bonds (topology), build context
ProteinConformation: positions, protonation state

### GeometryResult (requires: nothing)
Ring geometry (center, normal, radius, vertices), bond geometry (midpoint,
length, direction), bounding box, center of geometry, radius of gyration,
pre-built collections (rings_by_type, bonds_by_category, residues_by_type),
ring pair properties.

### DsspResult (requires: nothing)
Per-residue: secondary structure (1-char, H/G/I/E/B/T/S/C), phi (radians),
psi (radians), SASA (Angstroms^2), H-bond acceptors[2] with energy
(kcal/mol), H-bond donors[2] with energy (kcal/mol).

**WriteFeatures output:**
- `dssp_backbone.npy` (N, 5): per-atom backbone geometry
- `dssp_ss8.npy` (N, 8): 8-class secondary structure one-hot
  (H=α-helix, G=3₁₀-helix, I=π-helix, E=strand, B=isolated bridge,
  T=turn, S=bend, C=coil)
- `dssp_hbond_energy.npy` (N, 4): H-bond energies (acc0, acc1, don0, don1)
- `dssp_chi.npy` (N, 12): chi1-4 dihedral angles as (cos, sin, exists) × 4

### ChargeAssignmentResult (requires: nothing)
Per-atom: partial charge (e), VdW radius (Angstroms).

### MopacResult

**Requires:** nothing (external tool, conformation electronic structure precondition)
**Accessor:** `conformation.Result<MopacResult>()`

#### What it computes
PM7+MOZYME semiempirical electronic structure for the full protein.
Calls MOPAC as subprocess (~45s for 889 atoms). Provides conformation-
dependent charges and bond orders that downstream calculators draw on.

#### Input
Atom positions and elements from ProteinConformation, net charge.

#### What it stores (per atom on ConformationAtom)
- mopac_charge: double, Mulliken charge (elementary charge units)
- mopac_s_pop: double, s-orbital population (electrons)
- mopac_p_pop: double, p-orbital population (electrons)
- mopac_valency: double, sum of Wiberg bond orders
- mopac_bond_neighbours: vector of MopacBondNeighbour, sorted descending by order

#### What it stores (on the result object)
- Per-bond Wiberg bond orders: O(1) lookup via BondOrder(atom_a, atom_b)
- Topology bridge: TopologyBondOrders() parallel to protein.Bonds()
- Heat of formation (kcal/mol), dipole moment (Vec3)

#### Physics query methods
- `ChargeAt(atom_index) -> double`
- `SPopAt(atom_index) -> double`
- `PPopAt(atom_index) -> double`
- `ValencyAt(atom_index) -> double`
- `BondOrder(atom_a, atom_b) -> double` (O(1), symmetric, 0.0 if no bond)
- `TopologyBondOrder(bond_index) -> double` (parallel to protein.Bonds())
- `HeatOfFormation() -> double`

### MopacCoulombResult

**Requires:** MopacResult, SpatialIndexResult
**Accessor:** `conformation.Result<MopacCoulombResult>()`

Coulomb EFG from MOPAC QM charges. Same kernel as CoulombResult
(Calculator 4 in GEOMETRIC_KERNEL_CATALOGUE.md), different charge source.
See Calculator 9 in the catalogue for full description.

**TBD:** ConformationAtom field table, per-atom storage details,
WriteFeatures output manifest. Fields exist in code (ConformationAtom.h,
mopac_coulomb_ prefix). See src/MopacCoulombResult.cpp for implementation.

### MopacMcConnellResult

**Requires:** MopacResult, SpatialIndexResult, GeometryResult
**Accessor:** `conformation.Result<MopacMcConnellResult>()`

McConnell bond anisotropy tensor weighted by MOPAC Wiberg bond order.
Same kernel as McConnellResult (Calculator 3), each bond weighted by
measured electron sharing. See Calculator 10 in the catalogue.

**TBD:** ConformationAtom field table, per-atom storage details,
WriteFeatures output manifest. Fields exist in code (ConformationAtom.h,
mopac_mc_ prefix). See src/MopacMcConnellResult.cpp for implementation.

### EnrichmentResult (requires: nothing)
Per-atom: role (AtomRole), hybridisation (Hybridisation), categorical
booleans (is_backbone, is_amide_H, is_alpha_H, is_methyl, is_aromatic_H,
is_on_aromatic_residue, is_hbond_donor, is_hbond_acceptor, parent_is_sp2).
Pre-built collection: atoms_by_role.

### SpatialIndexResult (requires: GeometryResult)
KD-trees (atom positions, ring centers, bond midpoints), neighbour lists
(all atoms within 15A, with stored distance and direction).

### ApbsFieldResult (requires: ChargeAssignmentResult)
Per-atom: APBS solvated E-field (Vec3), APBS EFG (Mat3 + SphericalTensor).

### MolecularGraphResult (requires: SpatialIndexResult)
Per-atom: graph_dist_ring, graph_dist_N, graph_dist_O, eneg_sum_1,
eneg_sum_2, n_pi_bonds_3, is_conjugated, bfs_to_nearest_ring_atom,
bfs_decay.

### BiotSavartResult (requires: SpatialIndexResult, GeometryResult)
Per-atom per-ring (RingNeighbourhood): G_tensor, G_spherical, B_field,
B_cylindrical, gaussian_density.
Per-atom totals: total_B_field, total_G_tensor, total_G_spherical,
per_type_G_T0_sum, per_type_G_T2_sum, ring counts, ring distances,
exp-weighted sums, variance.
Per-ring: total_B_at_center, intensity_used, diagnostics, mutual B.
Per-atom shielding: bs_shielding_contribution (SphericalTensor, ppm).

### HaighMallionResult (requires: SpatialIndexResult, GeometryResult)
Per-atom per-ring (RingNeighbourhood): hm_tensor, hm_spherical, hm_B_field.
Per-atom totals: per_type_hm_T0_sum, per_type_hm_T2_sum.
Per-atom shielding: hm_shielding_contribution (SphericalTensor, ppm).

### PiQuadrupoleResult (requires: SpatialIndexResult, GeometryResult)
Per-atom per-ring (RingNeighbourhood): quad_tensor, quad_spherical,
quad_scalar.
Per-atom shielding: piquad_shielding_contribution (SphericalTensor, ppm).

### RingSusceptibilityResult (requires: SpatialIndexResult, GeometryResult)
Per-atom per-ring (RingNeighbourhood): chi_tensor, chi_spherical,
chi_scalar.
Per-atom shielding: ringchi_shielding_contribution (SphericalTensor, ppm).

### DispersionResult (requires: SpatialIndexResult, GeometryResult)
Per-atom per-ring (RingNeighbourhood): disp_tensor, disp_spherical,
disp_scalar, disp_contacts.
Per-atom shielding: disp_shielding_contribution (SphericalTensor, ppm).

### McConnellResult (requires: SpatialIndexResult, GeometryResult)
Per-atom per-bond (BondNeighbourhood): dipolar_tensor, dipolar_spherical,
mcconnell_scalar.
Per-atom totals: mcconnell_co_sum, mcconnell_cn_sum,
mcconnell_sidechain_sum, mcconnell_aromatic_sum, mcconnell_co_nearest,
T2_CO_nearest, T2_CN_nearest, T2_backbone_total, T2_sidechain_total,
T2_aromatic_total, nearest_CO_midpoint, nearest_CO_dist, nearest_CN_dist,
dir_nearest_CO.
Per-atom shielding: mc_shielding_contribution (SphericalTensor, ppm).

### CoulombResult (requires: ChargeAssignmentResult, SpatialIndexResult)
Per-atom: coulomb_E_total, _backbone, _sidechain, _aromatic,
coulomb_EFG_total (+ spherical), _backbone (+ spherical), _aromatic
(+ spherical), _solvent (derived: APBS - vacuum), coulomb_E_magnitude,
_bond_proj, _backbone_frac, aromatic_E_magnitude, _bond_proj,
aromatic_n_sidechain_atoms.
Per-atom shielding: coulomb_shielding_contribution (SphericalTensor, ppm).

### HBondResult (requires: DsspResult, SpatialIndexResult)
Per-atom: hbond_nearest_dist, _dir, _tensor, _spherical, _inv_d3,
_is_backbone, _count_within_3_5A, _is_donor, _is_acceptor.
Per-atom shielding: hbond_shielding_contribution (SphericalTensor, ppm).

### OrcaShieldingResult (requires: nothing)
Per-atom: diamagnetic + paramagnetic + total shielding tensors (Mat3 +
SphericalTensor each). Parsed from ORCA NMR output with element-verified
atom ordering. Each protein gets its own result — WT and mutant are
separate Proteins. Comparison is MutantProteinConformationComparison.

### SasaResult (requires: SpatialIndexResult)
Per-atom Shrake-Rupley solvent-accessible surface area using a Fibonacci
sphere (~92 points, r_vdW + r_probe). Each point checked for occlusion
by nearby atoms via SpatialIndexResult. SASA = fraction_exposed × sphere_area.

**Parameters:** `sasa_probe_radius` (default 1.4 Å), `sasa_n_points` (default 92).
**Physics query:** `AtomSASA(i) -> double`, `AllSASA() -> vector<double>`.
**WriteFeatures output:** `atom_sasa.npy` (N,) float64, Angstroms² per atom.

### AIMNet2Result (requires: SpatialIndexResult; model loaded externally)
Neural network charges and EFG via AIMNet2 (libtorch, CUDA mandatory).
Produces per-atom Hirshfeld charges, 256-dim AIM embedding, and Coulomb
EFG decomposed by source using the same dipolar kernel as CoulombResult.
The model (.jpt) is loaded once at startup and shared across conformations.
Failure is non-silent — returns nullptr and logs an error.

**CUDA required.** No CPU fallback. Unknown elements → nullptr.

**charge_sensitivity:** computed via autograd (d(charges)/d(positions)) only
if `aimnet2_sensitivity_mode = "autograd"` in TOML. Default: not computed.
The perturbation approach was removed — it produced conformation-specific
values that don't generalise.

**WriteFeatures output:**
- `aimnet2_charges.npy` (N,): Hirshfeld charges (elementary charge)
- `aimnet2_aim.npy` (N, 256): electronic embedding
- `aimnet2_efg.npy` (N, 9): Coulomb EFG total (SphericalTensor)
- `aimnet2_efg_aromatic.npy` (N, 9): Coulomb EFG from aromatic atoms
- `aimnet2_efg_backbone.npy` (N, 9): Coulomb EFG from backbone atoms

### WaterFieldResult (requires: SpatialIndexResult + SolventEnvironment)
Per-atom explicit-water E-field and EFG. Same Coulomb kernel as
CoulombResult but summing over water charges (TIP3P: O = -0.834e,
H = +0.417e). Builds its own spatial index over water oxygen positions.
Available in trajectory mode only (requires SolventEnvironment with water
positions for the frame).

**What it captures that APBS doesn't:** water orientation fluctuations,
bridging water, structural water, cavity effects.

**WriteFeatures output:**
- `water_efield.npy` (N, 3): total E-field from all water within cutoff (V/Å)
- `water_efield_first.npy` (N, 3): first-shell E-field (< 3.5 Å from O)
- `water_efg.npy` (N, 9): total EFG (SphericalTensor)
- `water_efg_first.npy` (N, 9): first-shell EFG
- `water_shell_counts.npy` (N, 2): [n_first, n_second] water O counts

### HydrationShellResult (requires: SpatialIndexResult + SolventEnvironment)
Per-atom hydration shell geometry from explicit water positions. Characterises
the local dielectric environment that no geometry-only kernel can see.

**What it computes:**
- Half-shell asymmetry: fraction of first-shell waters on the solvent-exposed
  side vs the protein interior (UCSB method).
- Mean water dipole orientation: cos(angle) between water dipole and
  atom→water O vector, averaged over first shell (water order parameter).
- Nearest ion distance and charge.

**WriteFeatures output:**
- `hydration_shell.npy` (N, 4): [asymmetry, dipole_cos, ion_dist, ion_charge]

### GromacsEnergyResult (requires: nothing; reads .edr file at runtime)
Per-frame aggregate energy terms from the GROMACS simulation. Reads the .edr
energy file and finds the frame nearest to the conformation's time (ps).
Whole-system quantities (not per-atom).

**GromacsEnergy struct fields:**
| Field | Unit | Description |
|-------|------|-------------|
| time_ps | ps | Frame time |
| coulomb_sr | kJ/mol | PME real-space Coulomb |
| coulomb_recip | kJ/mol | PME reciprocal-space Coulomb |
| coulomb_14 | kJ/mol | 1-4 intramolecular Coulomb |
| lj_sr | kJ/mol | Lennard-Jones short-range |
| potential | kJ/mol | Total potential energy |
| temperature | K | System temperature |
| pressure | bar | System pressure |
| volume | nm³ | Box volume |

`CoulombTotal() = coulomb_sr + coulomb_recip` (excludes 1-4).

**WriteFeatures output:**
- `gromacs_energy.npy` (1, 9): [t, Coul_SR, Coul_recip, Coul_14, LJ, pot, T, P, V]

### Calibration Pipeline (external — not a ConformationResult)

The ~93 tuneable calculator parameters are calibrated by the Python
e3nn model (learn/c_equivariant/) against DFT WT-ALA delta tensors.
The flow is: C++ geometric kernels → NPY features (via WriteFeatures)
→ Python e3nn calibration → TOML parameters → C++ (next run).

The calibrated values override literature defaults for ring current
intensities, bond anisotropies, Buckingham coefficients, and other
calculator-specific parameters. See CALCULATOR_PARAMETER_API.md for
the full parameter set (93 parameters with equations and references).

Feature extraction is distributed: each ConformationResult writes its
own NPY arrays via WriteFeatures(). ConformationResult::WriteAllFeatures()
traverses all attached results. There is no centralized
FeatureExtractionResult — the pattern is simpler and keeps each
result responsible for its own output.

---

## Extended: Trajectory-scope detail (catalog + pending H5 schema)

Material below is true and planned but digressive for the core
model section above. Living checklists and pending schemas.

### TrajectoryResult catalog (living)

Grouped by configuration. ✓ = landed; ⏳ = pending. Source of truth
for the pending rows is
`spec/pending_include_trajectory_scope_2026-04-22.md` Appendix F.

**`PerFrameExtractionSet` (production canonical):**

| ✓/⏳ | Type | Source ConformationResult | Lifecycle | Emission |
|-----|------|---------------------------|-----------|----------|
| ✓ | `PositionsTimeSeriesTrajectoryResult` | (positions) | FO | `DenseBuffer<Vec3>` → `/trajectory/positions/` |
| ✓ | `BsWelfordTrajectoryResult` | BiotSavart | AV | TrajectoryAtom Welford fields + `/trajectory/bs_welford/` |
| ✓ | `BsShieldingTimeSeriesTrajectoryResult` | BiotSavart | FO | `DenseBuffer<SphericalTensor>` → `/trajectory/bs_shielding_time_series/` |
| ✓ | `BsAnomalousAtomMarkerTrajectoryResult` | BsWelford + BiotSavart | AV | per-atom events bag |
| ✓ | `BsT0AutocorrelationTrajectoryResult` | BiotSavart | FO | `DenseBuffer<double>` → `/trajectory/bs_t0_autocorrelation/` |
| ✓ | `BondLengthStatsTrajectoryResult` | (positions) | AV | internal per-bond Welford → `/trajectory/bond_length_stats/` |
| ⏳ | `HmWelfordTrajectoryResult` | HaighMallion | AV | Welford fields + group |
| ⏳ | `HmShieldingTimeSeriesTrajectoryResult` | HaighMallion | FO | `DenseBuffer<SphericalTensor>` |
| ⏳ | `McConnellWelfordTrajectoryResult` | McConnell | AV | Welford fields + per-category sums |
| ⏳ | `McConnellShieldingTimeSeriesTrajectoryResult` | McConnell | FO | `DenseBuffer<SphericalTensor>` |
| ⏳ | `CoulombFieldTimeSeriesTrajectoryResult` | Coulomb | FO | E-field Vec3 + EFG Mat3 dense |
| ⏳ | `ApbsFieldTimeSeriesTrajectoryResult` | ApbsField | FO | `DenseBuffer<Vec3>` + `DenseBuffer<Mat3>` |
| ⏳ | `WaterEnvironmentTimeSeriesTrajectoryResult` | WaterField | FO | multiple dense buffers |
| ⏳ | `HydrationShellTimeSeriesTrajectoryResult` | HydrationShell | FO | |
| ⏳ | `HydrationGeometryTimeSeriesTrajectoryResult` | HydrationGeometry | FO | |
| ⏳ | `AIMNet2ChargeTimeSeriesTrajectoryResult` | AIMNet2 | FO | per-atom charges over time |
| ⏳ | `AIMNet2EmbeddingTimeSeriesTrajectoryResult` | AIMNet2 | FO | 256-dim per atom per frame — optional, large |
| ⏳ | `EeqChargeWelfordTrajectoryResult` | Eeq | AV | |
| ⏳ | `SasaTimeSeriesTrajectoryResult` | Sasa | FO | |
| ⏳ | `SasaWelfordTrajectoryResult` | Sasa | AV | |
| ⏳ | `HBondTimeSeriesTrajectoryResult` | HBond | FO | |
| ⏳ | `HBondCountWelfordTrajectoryResult` | HBond | AV | |
| ⏳ | `DihedralTimeSeriesTrajectoryResult` | Dssp | FO | per-residue dihedrals |
| ⏳ | `DihedralBinTransitionTrajectoryResult` | Dssp | AV | per-residue transition counters |
| ⏳ | `Dssp8TimeSeriesTrajectoryResult` | Dssp | FO | per-residue SS code + H-bond energies |
| ⏳ | `Dssp8TransitionTrajectoryResult` | Dssp | AV | per-residue SS transitions |
| ⏳ | `BondedEnergyTimeSeriesTrajectoryResult` | BondedEnergy | FO | per-atom energy decomposition |
| ⏳ | `GromacsEnergyTimeSeriesTrajectoryResult` | GromacsEnergy | FO | per-frame aggregate (no atom axis) |
| ⏳ | `PiQuadrupoleShieldingTimeSeriesTrajectoryResult` | PiQuadrupole | FO | |
| ⏳ | `RingSusceptibilityShieldingTimeSeriesTrajectoryResult` | RingSusceptibility | FO | |
| ⏳ | `DispersionShieldingTimeSeriesTrajectoryResult` | Dispersion | FO | |
| ⏳ | `RingNeighbourhoodTrajectoryStats` | multiple ring calculators | FO | rich per-atom-per-ring struct vectors (Pattern A) |

**`ScanForDftPointSet` (scan mode, for DFT pose selection):** *Slated
for removal; a later doc pass will handle the replacement.*

| ✓/⏳ | Type | Lifecycle | Emission |
|-----|------|-----------|----------|
| ✓ | `ChiRotamerSelectionTrajectoryResult` | AV | run-scope SelectionBag push per transition |
| ⏳ | `RmsdTrackingTrajectoryResult` | AV | per-frame RMSD vs reference frame |
| ⏳ | `RmsdSpikeSelectionTrajectoryResult` | AV | SelectionBag push on threshold crossing |
| ⏳ | `Dssp8TransitionTrajectoryResult` | AV | (shared with PerFrameExtractionSet) |
| ⏳ | `DftPoseCoordinatorTrajectoryResult` | FO | reads other emitters' SelectionRecords at Finalize, deduplicates, pushes reduced set back under its own kind |

**`FullFatFrameExtraction` (selected-frame MOPAC):**

| ✓/⏳ | Type | Lifecycle | Emission |
|-----|------|-----------|----------|
| ⏳ | `MopacChargeWelfordTrajectoryResult` | AV (sparse) | |
| ⏳ | `MopacBondOrderWelfordTrajectoryResult` | AV (sparse) | internal per-bond (not TrajectoryBond) |
| ⏳ | `MopacCoulombShieldingTimeSeriesTrajectoryResult` | FO (sparse) | |
| ⏳ | `MopacMcConnellShieldingTimeSeriesTrajectoryResult` | FO (sparse) | |
| ⏳ | `MopacVsFf14SbReconciliationTrajectoryResult` | FO | per-atom \|cos\| of MOPAC-Coulomb T2 vs DFT-delta T2 |

### H5 metadata schema (pending / partial)

The design pass spec'd a full `/metadata/source/` + `/atoms/identity/`
schema (see pending-include file §7). What's currently emitted:

- Per-Result group attributes: `result_name`, `n_frames`, `finalized`,
  and Result-specific semantics (e.g. `irrep_layout`, `normalization`,
  `parity`, `units`, `estimator`, `mean_convention`). These ARE
  schema contracts — downstream Python consumers parse them verbatim.
- `/atoms/` group with `element`, `residue_index`, `pdb_atom_name`
  emitted by `TrajectoryProtein::WriteH5`. Minimal passthrough from
  `Protein`.
- `/trajectory/source/` group attributes: `xtc_path`, `tpr_path`,
  `edr_path`, `configuration`.
- `/trajectory/frames/` datasets: `time_ps`, `original_index`.
- `/trajectory/selections/<kind>/` per-kind sub-groups walked via
  `selections_.Kinds()`.

Pending full schema (awaits `NmrAtomIdentity` landing):

- `/metadata/source/` with protein provenance attrs (pdb_source,
  deposition_date, crystal_resolution, protonation_tool, force_field,
  stripped, assumptions, extractor_version, extraction_date) from
  `ProteinBuildContext`.
- `/atoms/identity/` with the typed NmrAtomIdentity fields
  (iupac_atom_position, nmr_class, locant, methyl_group,
  chi_participation, ring_atom_role, ring_membership, residue_category,
  hydropathy, etc.) for BMRB / RefDB binding downstream.

Emission is currently sufficient for the Python SDK and viewer to
consume the implemented handlers' output. The full schema is a
pre-fleet-activation bundle flagged PROPOSAL-PENDING-USER-REVIEW in
the pending-include file.

### Deferred work pointers

- `NmrAtomIdentity` on Protein (§2 of pending-include, Appendix A for
  the enum generator) — the typed identity layer for external-data
  binding. Proposal-pending user review.
- Full TrajectoryResult catalog (Appendix F above) — ~25 classes to
  land, each cloning one of the seven canonical shapes.
- `TopologyTrajectoryResult` or equivalent — centralized topology
  (bonds, rings, parent atoms) emission so bond-scope Results stop
  duplicating topology in their own H5 groups.
- Byte-parity validation pass against archived pre-refactor output
  across 10 calibration proteins, gated pre-fleet-activation.
- NVRTC rpath fix (see memory entry `reference_nvrtc_rpath_fix`)
  so `TrajectoryRunDrivesLoop` runs without external
  `LD_LIBRARY_PATH`.


---

# Addition — Trajectory-scope entities (2026-04-24, library-reference form)

This section is an addition, not a replacement for the earlier
`## Trajectory-scope entities` and `## Extended: Trajectory-scope
detail` sections above. It describes the trajectory-scope object
model against the current tree (2026-04-24) in library-manual
form. Where this section disagrees with the earlier sections, trust
this section as current truth; the earlier sections are retained for
design history. Known trajectory-handling errors in the earlier
sections are flagged separately for review.

## Scope

Trajectory-scope classes are modular accumulators built on top of the
per-frame `ConformationResult` pipeline. The static-analysis path
(`Protein` → `ProteinConformation` → `ConformationResult` →
`OperationRunner::Run`) is unchanged and authoritative.

- `Protein` — the invariant protein.
- `ProteinConformation` — the recorded and calculated geometry of that
  protein (positions recorded; ring geometries, bond geometries,
  spatial indices, `ConformationResult`s, `ConformationAtom` fields
  all calculated from the positions).
- `TrajectoryProtein` — the concept of a protein in the context of a
  trajectory. Contains the invariant `Protein`, basic geometric
  properties calculated on conformation 0, the per-atom `TrajectoryAtom`
  buffer holding the geometric record over time, and the attached
  `TrajectoryResult`s holding the kernel record over time.
- `Trajectory` — the process that drives one traversal of a source
  (XTC + TPR + EDR). Process during `Run`; record after (frame times,
  frame indices, run-scope selection stream, last-frame env).

---

## TrajectoryProtein

Wraps a `Protein`, seats conformation 0 on it, allocates and owns the
parallel `TrajectoryAtom` vector, and holds the attached
`TrajectoryResult`s and any `DenseBuffer<T>`s they transfer at
`Finalize`. Non-copyable.

### Fields

| Field                   | Type                                                      | Description                                       |
|-------------------------|-----------------------------------------------------------|---------------------------------------------------|
| `protein_`              | `unique_ptr<Protein>`                                     | Wrapped invariant protein                         |
| `charges_`              | `unique_ptr<ChargeSource>`                                | TPR-derived charge source                         |
| `sys_reader_`           | `FullSystemReader`                                        | Topology + frame-split helper, borrowed by handler|
| `bonded_params_`        | `BondedParameters`                                        | Force-field bonded terms from TPR                 |
| `atoms_`                | `vector<TrajectoryAtom>`                                  | Per-atom store, parallel to `Protein.Atoms()`     |
| `results_`              | `unordered_map<type_index, unique_ptr<TrajectoryResult>>` | One TR per type                                   |
| `results_attach_order_` | `vector<TrajectoryResult*>`                               | Attach order = dispatch order                     |
| `dense_buffers_`        | `unordered_map<type_index, unique_ptr<DenseBufferBase>>`  | Keyed by owning TR's type                         |
| `finalized_`            | `bool`                                                    | Set by `FinalizeAllResults`                       |

### Build and seed

| Method                                       | Role                                                                                   |
|----------------------------------------------|----------------------------------------------------------------------------------------|
| `BuildFromTrajectory(dir_path)`              | Parse `dir_path/md.tpr`; build Protein + charges + bonded params. Does not seat conf0. |
| `Seed(first_frame_positions, time_ps)`       | `Protein::FinalizeConstruction` (bonds + rings), seat conf0 via `AddMDFrame`, allocate `TrajectoryAtom`s. |
| `Error()`                                    | Last-error string from `BuildFromTrajectory`.                                          |

### Conformation access

| Method                                       | Role                                                                     |
|----------------------------------------------|--------------------------------------------------------------------------|
| `CanonicalConformation()`                    | `protein_->ConformationAt(0)`. Same object `Protein::Conformation()` returns. |
| `TickConformation(positions)`                | Ephemeral per-frame `unique_ptr<ProteinConformation>` pointing at the wrapped Protein. |
| `ProteinRef()`                               | The wrapped Protein.                                                     |
| `AtomCount()`, `AtomAt(i)`, `MutableAtomAt(i)`, `Atoms()` | Per-atom access.                                            |

### Results

| Method                                       | Role                                                                     |
|----------------------------------------------|--------------------------------------------------------------------------|
| `AttachResult(unique_ptr<TR>)`               | Singleton-per-type check; returns false if already attached. Dependency validation is in `Trajectory::Run` Phase 4, not here. |
| `Result<T>()`, `HasResult<T>()`, `AllResults()`, `ResultsInAttachOrder()` | Typed access; attach order = dispatch order. |
| `DispatchCompute(conf, traj, frame_idx, time_ps)` | Iterates `results_attach_order_`, calls each TR's `Compute`.        |
| `FinalizeAllResults(traj)`                   | Iterates, calls `Finalize`; sets `finalized_`.                           |
| `AdoptDenseBuffer<T>(buffer, owner_type)`    | Move ownership of a `DenseBuffer<T>` into `dense_buffers_`, keyed by owner TR's `type_index`. |
| `GetDenseBuffer<T>(owner_type)`              | Retrieve; nullptr on missing or element-type mismatch.                   |

### Serialisation

| Method                                       | Role                                                                     |
|----------------------------------------------|--------------------------------------------------------------------------|
| `WriteH5(file)`                              | File-root attributes (`protein_id`, `n_atoms`, `finalized`); `/atoms/{element, residue_index, pdb_atom_name}`; delegates to each TR's `WriteH5Group`. |
| `WriteFeatures(dir)`                         | Delegates NPY emission per TR.                                           |

Identity goes through the wrapped Protein, not duplicated here:

```cpp
const TrajectoryAtom& ta = tp.AtomAt(42);              // trajectory record
const Atom& a           = tp.ProteinRef().AtomAt(42);  // element, bonds, residue
```

---

## TrajectoryAtom

Per-atom cell of the trajectory buffer. Private constructor; only
`TrajectoryProtein` constructs via `friend`. One instance per wrapped
Protein atom, allocated once by `Seed`, never resized.

Per-atom trajectory data lives in two shapes in this struct:

**Event bag** — `RecordBag<AtomEvent> events`. Open-ended per-atom
record stream tagged by `(emitter, kind, frame_idx, time_ps)`.
Queried by kind (`events.ByKind<T>()`, `events.ByKindSinceFrame<T>(...)`,
`events.CountByKind<T>()`). What an emitter puts in the bag is the
emitter's choice; later TRs consume by kind. The per-atom events bag
is not currently H5-emitted — events are accessible in-process for
cross-TR reads and test assertions.

**Finalized rollup fields** — typed fields written by one TR each,
for rolled-up statistics where a record per frame would be wasteful.
The BS Welford set written by `BsWelfordTrajectoryResult`:

| Field                                              | Type       | Writer                        |
|----------------------------------------------------|------------|-------------------------------|
| `bs_t0_mean`, `bs_t0_m2`, `bs_t0_std`              | `double`   | `BsWelfordTrajectoryResult`   |
| `bs_t0_min`, `bs_t0_max`                           | `double`   | `BsWelfordTrajectoryResult`   |
| `bs_t0_min_frame`, `bs_t0_max_frame`               | `size_t`   | `BsWelfordTrajectoryResult`   |
| `bs_n_frames`                                      | `size_t`   | `BsWelfordTrajectoryResult`   |
| `bs_t2mag_{mean, m2, std, min, max, min_frame, max_frame}` | `double`/`size_t` | `BsWelfordTrajectoryResult` |
| `bs_t0_delta_{mean, m2, std, min, max, n}`         | `double`/`size_t` | `BsWelfordTrajectoryResult` |

`_m2` fields carry the running Welford sum-of-squared-deviations;
`_std` is populated at `Finalize` from `m2 / (n - 1)`. `_std` is
undefined mid-stream.

Accumulator implementation objects (`Welford` structs, rolling
windows, full-history buffers) live inside the owning TR, not on
`TrajectoryAtom`. Identity (element, residue, bonds, ring
membership) is not duplicated here — reach through
`tp.ProteinRef().AtomAt(i)`.

---

## TrajectoryResult (ABC)

Base class for per-trajectory modular calculators. Parallel to
`ConformationResult` at conformation scope.

### Virtual interface

| Method                                                      | Required       | Role                                                                          |
|-------------------------------------------------------------|----------------|-------------------------------------------------------------------------------|
| `Name() -> string`                                          | pure           | For logs and attach diagnostics.                                              |
| `Dependencies() -> vector<type_index>`                      | pure           | TR or ConformationResult types that must precede. Validated in Phase 4.       |
| `Compute(conf, tp, traj, frame_idx, time_ps)`               | pure           | Per-frame work.                                                               |
| `Finalize(tp, traj)`                                        | default no-op  | End-of-stream synthesis.                                                      |
| `WriteFeatures(tp, dir) -> int`                             | default 0      | NPY emission.                                                                 |
| `WriteH5Group(tp, file)`                                    | default no-op  | H5 group emission.                                                            |

TRs construct via static `Create(const TrajectoryProtein&)` factories
returning `unique_ptr<TrajectoryResult>`. Factories run after
`tp.Seed`, so they may size internal buffers from `tp.AtomCount()`,
`tp.ProteinRef().BondCount()`, or `tp.ProteinRef().RingCount()`.

During `Compute` a TR reads the current `ProteinConformation`
(its `ConformationResult` output and `ConformationAtom` fields) and
writes into one or more of:

- its own internal accumulator state (private to the TR);
- `tp.MutableAtomAt(i)` fields it owns (one writer per field);
- `tp.MutableAtomAt(i).events` (per-atom event bag);
- `traj.MutableSelections()` (run-scope selection bag);
- internal per-bond / per-ring / per-residue state on itself.

At `Finalize` a TR converts accumulator state to output (Welford m2
to std, full history to ACF, counter to frequency), transfers a
`DenseBuffer<T>` to TP via `AdoptDenseBuffer`, or reads the
selection bag and pushes a reduced set under its own kind.

### Lifecycle shapes

**AV (always valid mid-stream).** `Compute` updates output fields in
place each frame. `Finalize` at most converts running variance to
std. Mid-run snapshots of `tp.AtomAt(i)` are meaningful.

**FO (Finalize only).** `Compute` appends to an internal buffer; the
output materialises at `Finalize` (dense buffer transfer, ACF from
full history, FFT). Mid-run snapshots are not meaningful.

### Current TRs

The TRs below illustrate the patterns in `PATTERNS.md §§13-18`. They
are worked examples of the available lifecycles and output shapes;
they are NOT a prescribed template set. New calculators clone the
most relevant exemplar and adapt — see `PATTERNS.md §17` on
duplication over chaining.

| Class                                           | Lifecycle      | `Dependencies()`                                         | Output                                                               |
|-------------------------------------------------|----------------|----------------------------------------------------------|----------------------------------------------------------------------|
| `BsWelfordTrajectoryResult`                     | AV             | `BiotSavartResult`                                       | `TrajectoryAtom` rollup fields + `/trajectory/bs_welford/`           |
| `BsShieldingTimeSeriesTrajectoryResult`         | FO             | `BiotSavartResult`                                       | `DenseBuffer<SphericalTensor>` → `/trajectory/bs_shielding_time_series/{xyz, frame_indices, frame_times}` |
| `BsAnomalousAtomMarkerTrajectoryResult`         | AV (emitter)   | `BsWelfordTrajectoryResult`, `BiotSavartResult`          | Per-atom events pushed to `ta.events`, kinds `BsAnomalyHighT0` / `BsAnomalyLowT0`. No H5 group of its own (events bag not H5-emitted yet). |
| `BsT0AutocorrelationTrajectoryResult`           | FO             | `BiotSavartResult`                                       | `DenseBuffer<double>` → `/trajectory/bs_t0_autocorrelation/{rho, lag_frames, lag_times_ps}` |
| `BondLengthStatsTrajectoryResult`               | AV             | none (reads `Bond` indices + positions)                  | Internal per-bond Welford → `/trajectory/bond_length_stats/`         |
| `PositionsTimeSeriesTrajectoryResult`           | FO             | none (reads positions)                                   | `DenseBuffer<Vec3>` → `/trajectory/positions/{xyz, frame_indices, frame_times}` |

Current `RunConfiguration` attach status:

| Configuration           | Attached TRs                                                                                   |
|-------------------------|-----------------------------------------------------------------------------------------------|
| `ScanForDftPointSet`    | `BsWelfordTrajectoryResult`                                                                   |
| `PerFrameExtractionSet` | All six in the table above                                                                     |
| `FullFatFrameExtraction`| Same as `PerFrameExtractionSet`                                                               |

---

## Trajectory

Process entity representing one traversal of a source (XTC + TPR +
EDR). Non-copyable. Holds source paths, preloaded EDR frames, the
single-slot per-frame environment, the run record, and during `Run`
the handler.

### Fields

| Field                                         | Type                              | Description                                            |
|-----------------------------------------------|-----------------------------------|--------------------------------------------------------|
| `xtc_path_`, `tpr_path_`, `edr_path_`         | `filesystem::path`                | Source paths.                                          |
| `edr_frames_`                                 | `vector<GromacsEnergy>`           | Preloaded at construction; `EnergyAtTime(t)` is O(log T). |
| `env_`                                        | `TrajectoryEnv`                   | Single-slot per-frame environment.                     |
| `frame_count_`, `frame_times_`, `frame_indices_` | `size_t` / `vector<double>` / `vector<size_t>` | Run record.                           |
| `selections_`                                 | `RecordBag<SelectionRecord>`      | Run-scope event stream.                                |
| `output_dir_`                                 | `filesystem::path`                | Recorded for post-run writers.                         |
| `handler_`                                    | `unique_ptr<GromacsFrameHandler>` | During-Run only; released at end of Phase 8.           |
| `state_`                                      | enum                              | `Constructed` → `Running` → `Complete`.                |

Constructor preloads EDR (`LoadEdr`) so frame-energy lookup is
available from Phase 6 onwards.

### Run (eight phases)

`Status Run(tp, config, session, extras = {}, output_dir = {})`:

1. **Open handler.** `handler_ = make_unique<GromacsFrameHandler>(tp)`; `handler_->Open(xtc_path, tpr_path)` mounts the stream and builds the PBC fixer.
2. **Read frame 0 and seed.** `handler_->ReadNextFrame()`; `tp.Seed(handler_->ProteinPositions(), handler_->Time())`.
3. **Attach TrajectoryResults.** Iterate `config.TrajectoryResultFactories()` then caller `extras`. Attach order is dispatch order.
4. **Validate dependencies and resources.** For each attached TR, every `type_index` in `Dependencies()` must be satisfied by another attached TR OR `config.RequiresConformationResult(t)`. If `config.RequiresAimnet2()`, `session.HasAimnet2Model()` must be true.
5. **Build base `RunOptions`.** From `config.PerFrameRunOptions()` + `tp.Charges()` + `tp.BondedParams()` + `session.Aimnet2Model()`.
6. **Frame 0.** Populate `env_` from handler + EDR; `OperationRunner::Run(tp.CanonicalConformation(), frame_opts)`; `tp.DispatchCompute(conf0, *this, 0, time)`; record.
7. **Per-frame loop.** `handler_->Skip()` (stride − 1) times; `handler_->ReadNextFrame()`; update `env_`; `tick = tp.TickConformation(positions)`; `OperationRunner::Run(*tick, frame_opts)`; `tp.DispatchCompute(*tick, *this, handler_->Index(), handler_->Time())`; record.
8. **Finalize.** `tp.FinalizeAllResults(*this)`. Selections were pushed during Compute / Finalize — no sweep here. `state_ = Complete`; release `handler_`.

### Return codes

From `errors.h`: `kOk`, `kXtcOpenFailed`, `kFrameReadFailed`,
`kAttachRejectedSingleton`, `kAttachDependencyUnmet`,
`kConfigRequiresAimnet2`, `kCalculatorPipelineFailed`. Every non-zero
return is paired with an `OperationLog::Error` diagnostic at the
failure site.

### Post-Run accessors

`IsComplete()`, `FrameCount()`, `TotalTimePs()`, `FrameTimes()`,
`FrameIndices()`, `Selections()`, `EnergyAtTime(t)`, `Env()`,
`XtcPath()`, `TprPath()`, `EdrPath()`, `OutputDir()`.

### H5 emission

`WriteH5(file)`:

- `/trajectory/source/` — attributes `xtc_path`, `tpr_path`, `edr_path`, `configuration`.
- `/trajectory/frames/` — datasets `time_ps` (T,), `original_index` (T,); attribute `n_frames`.
- `/trajectory/selections/<kind>/` — one group per kind in `selections_.Kinds()`, with `frame_idx`, `time_ps`, `reason` datasets and `n_records` attribute. The group path uses `kind.name()` (mangled C++ type name — compiler-dependent but stable within a build).

---

## TrajectoryEnv

Single-slot per-frame environment stash owned by `Trajectory`.
`Trajectory::Run` writes it at the top of each frame from the handler
and the preloaded EDR; calculators read it via `RunOptions` pointers,
TRs read it via `traj.Env()`. Overwritten each frame.

| Field                  | Type                    | Populated from                                  |
|------------------------|-------------------------|-------------------------------------------------|
| `solvent`              | `SolventEnvironment`    | `handler_->Solvent()`                           |
| `current_energy`       | `const GromacsEnergy*`  | `EnergyAtTime(handler_->Time())`; may be null   |
| `current_frame_idx`    | `size_t`                | `handler_->Index()` (0 on frame 0)              |
| `current_frame_time`   | `double`                | `handler_->Time()`                              |

A TR that needs cross-frame environment (running water-dipole
statistics, bridging-water histograms) keeps its own per-frame
buffer; `env_` is strictly this-frame.

---

## Session

Process-wide resources. Outlives any single `Trajectory`. Loaded
from TOML; carries large models (AIMNet2) loaded once and passed
into `Trajectory::Run` via the `session` parameter.

| Field             | Type                        | Description                                      |
|-------------------|-----------------------------|--------------------------------------------------|
| `aimnet2_model_`  | `unique_ptr<AIMNet2Model>`  | Loaded via `LoadAimnet2Model(path)`.             |
| `last_error_`     | `string`                    | Diagnostic from most recent failed Load.         |

Methods: `LoadFromToml()`, `LoadAimnet2Model(path)`, `Aimnet2Model()`,
`HasAimnet2Model()`, `LastError()`. Load methods return `Status`.

---

## RunConfiguration

Typed description of a trajectory-run shape. Three named static
factories.

### Fields

| Field                         | Type                                | Description                                          |
|-------------------------------|-------------------------------------|------------------------------------------------------|
| `name_`                       | `string`                            | For logs.                                            |
| `per_frame_opts_`             | `RunOptions`                        | Base `RunOptions` for `OperationRunner::Run`.        |
| `traj_factories_`             | `vector<TrajectoryResultFactory>`   | Attach order = dispatch order.                       |
| `required_conf_result_types_` | `unordered_set<type_index>`         | Checked in Phase 4 against each TR's `Dependencies()`. |
| `requires_aimnet2_`           | `bool`                              | Phase 4 enforces session has model loaded.           |
| `stride_`                     | `size_t`                            | Process every N-th frame (default 1).                |

### Factories

| Factory                     | Per-frame Conformation set                                                                                                                                                                                                                                                                                | Stride | Requires AIMNet2 | Use                                                        |
|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------|------------------|------------------------------------------------------------|
| `ScanForDftPointSet()` *(slated for removal)* | `GeometryResult`, `SpatialIndexResult`, `EnrichmentResult`, `DsspResult`, `BiotSavartResult`, `SasaResult`. MOPAC / APBS / Coulomb skipped.                                                                                                                                                                | 1      | no               | Rotamer / RMSD / bin-crossing detection for DFT pose selection. A later doc pass will handle the replacement. |
| `PerFrameExtractionSet()`   | `GeometryResult`, `SpatialIndexResult`, `EnrichmentResult`, `DsspResult`, `ChargeAssignmentResult`, `ApbsFieldResult`, `BiotSavartResult`, `HaighMallionResult`, `McConnellResult`, `RingSusceptibilityResult`, `PiQuadrupoleResult`, `DispersionResult`, `HBondResult`, `SasaResult`, `EeqResult`, `AIMNet2Result`, `WaterFieldResult`, `HydrationShellResult`, `HydrationGeometryResult`, `GromacsEnergyResult`, `BondedEnergyResult` (21 types). MOPAC skipped; vacuum Coulomb skipped (APBS supersedes). | 2      | yes              | Production canonical for 685-protein fleet.                |
| `FullFatFrameExtraction()`  | `PerFrameExtractionSet` with `skip_mopac = false`. MOPAC-family ConformationResult types are not yet in `required_conf_result_types_` — see `spec/pending_decisions_20260423.md` item 3.                                                                                                                  | 2      | yes              | Selected-frame MOPAC (DFT pose set, harvester checkpoints). |

Mutators for factory authoring: `SetName`, `MutablePerFrameRunOptions`,
`AddTrajectoryResultFactory`, `RequireConformationResult`,
`SetRequiresAimnet2`, `SetStride`.

---

## RecordBag&lt;Record&gt;

Typed container for a stream of tagged records. Used at two scopes:

- **Run-scope**: `Trajectory::selections_` as `RecordBag<SelectionRecord>`.
- **Atom-scope**: `TrajectoryAtom::events` as `RecordBag<AtomEvent>`.

`Record` must expose `kind` (`type_index`), `frame_idx` (`size_t`),
and `time_ps` (`double`) as top-level fields so windowed queries are
direct reads, not metadata lookups.

### API

| Method                                               | Role                                                             |
|------------------------------------------------------|------------------------------------------------------------------|
| `Push(record)`                                       | Append.                                                          |
| `All()`                                              | All records in push order.                                       |
| `Count()`                                            | Total count.                                                     |
| `Kinds()`                                            | Distinct kinds present, in first-seen order.                     |
| `ByKind(kind)` / `ByKind<T>()`                       | Records of one kind, in push order.                              |
| `ByKindSinceFrame(kind, from)` / `ByKindSinceFrame<T>(from)` | Records with `frame_idx >= from`.                          |
| `ByKindSinceTime(kind, t_ps)` / `ByKindSinceTime<T>(t_ps)`   | Records with `time_ps >= t_ps`.                            |
| `MostRecent(kind)` / `MostRecent<T>()`               | Most recent record of kind; nullptr if none.                     |
| `CountByKind(kind)` / `CountByKind<T>()`             | Count of one kind.                                               |

---

## SelectionRecord

Run-scope record type.

| Field        | Type                     | Note                              |
|--------------|--------------------------|-----------------------------------|
| `kind`       | `type_index`             | Emitting TR's `type_index`.       |
| `frame_idx`  | `size_t`                 |                                   |
| `time_ps`    | `double`                 |                                   |
| `reason`     | `string`                 | Human-readable cause.             |
| `metadata`   | `map<string, string>`    | Emitter-specific extras.          |

---

## AtomEvent

Atom-scope record type.

| Field        | Type                     | Note                                                 |
|--------------|--------------------------|------------------------------------------------------|
| `emitter`    | `type_index`             | TR that pushed.                                      |
| `kind`       | `type_index`             | Emitter-specific discriminator.                      |
| `frame_idx`  | `size_t`                 |                                                      |
| `time_ps`    | `double`                 |                                                      |
| `metadata`   | `map<string, string>`    | Emitter-specific extras.                             |

`emitter` and `kind` are separate because one emitter typically
pushes several kinds (e.g. `BsAnomalousAtomMarkerTrajectoryResult`
pushes `BsAnomalyHighT0` and `BsAnomalyLowT0`). Queries on the bag
discriminate on `kind`; `emitter` allows filtering by source when
relevant.

---

## DenseBuffer&lt;T&gt;

Contiguous per-atom × stride storage for FO `TrajectoryResult`
outputs. Atom-major layout: `storage[atom_idx * stride + offset]`.
Owned by `TrajectoryProtein` after a TR transfers ownership at
`Finalize`, keyed by the owning TR's `type_index`.

### Base and template

```cpp
class DenseBufferBase {
public:
    virtual ~DenseBufferBase() = default;
    virtual size_t AtomCount() const = 0;
    virtual size_t StridePerAtom() const = 0;
    virtual size_t ElementSizeBytes() const = 0;
};

template <typename T>
class DenseBuffer : public DenseBufferBase {
public:
    DenseBuffer(size_t atom_count, size_t stride_per_atom);
    T&       At(size_t atom_idx, size_t offset);
    const T& At(size_t atom_idx, size_t offset) const;
    T*       AtomSlicePtr(size_t atom_idx);
    const T* AtomSlicePtr(size_t atom_idx) const;
    T*       RawData();
    const T* RawData() const;
    size_t   TotalElementCount() const;
    // AtomCount() / StridePerAtom() / ElementSizeBytes() override the base.
};
```

Element types observed in current code: `double` (BsT0Autocorrelation
ρ(k)), `Vec3` (PositionsTimeSeries), `SphericalTensor`
(BsShieldingTimeSeries). `Mat3` reserved for future EFG
time-series TRs.

H5 emission flattens via explicit named-component access
(`.x() / .y() / .z()`, `.T0 / .T1[k] / .T2[k]`, `m(i,j)`), never
`reinterpret_cast` on the raw buffer. Each emitting TR writes its
own group and sets its own schema attributes — e.g.
`BsShieldingTimeSeriesTrajectoryResult` emits `irrep_layout`,
`normalization`, `parity`, `units`;
`BsT0AutocorrelationTrajectoryResult` emits `estimator`,
`mean_convention`.

---

## GromacsFrameHandler

Format-specific XTC/TPR reader. Mounts the stream, builds the PBC
fixer from the TPR, reads frames on demand. Pure reader: does not
create conformations, does not run calculators, does not write to
`traj.env_`, does not know `TrajectoryResult`s.

Additional trajectory formats would be sibling reader classes at
this layer, not a virtual-base hierarchy.

### Fields

| Field                 | Type                          | Description                                         |
|-----------------------|-------------------------------|-----------------------------------------------------|
| `tp_`                 | `TrajectoryProtein&`          | Borrowed; used for topology access via `SysReader`. |
| `reader_`             | `XtcStreamReader`             | XTC cursor.                                         |
| `wholer_`             | `unique_ptr<MoleculeWholer>`  | PBC fixer, built from TPR.                          |
| `protein_positions_`  | `vector<Vec3>`                | Last-read frame, protein slice.                     |
| `solvent_`            | `SolventEnvironment`          | Last-read frame, solvent slice.                     |
| `last_frame_time_`    | `double`                      |                                                     |
| `current_index_`      | `size_t`                      | XTC position; valid iff `has_read_`.                |
| `has_read_`           | `bool`                        |                                                     |

### Operations

| Method                                    | Role                                                                 |
|-------------------------------------------|----------------------------------------------------------------------|
| `Open(xtc_path, tpr_path)`                | Mount XTC, build wholer from TPR, sanity-check atom counts against `tp.SysReader().Topology().protein_count`. Does not read a frame. |
| `ReadNextFrame()`                         | Read one XTC frame, PBC-fix protein slice, split into `protein_positions_` + `solvent_`. Returns false at EOF. |
| `Skip()`                                  | Advance one frame without extracting; accessors return stale data afterward. |
| `Reopen()`                                | Reset cursor to start (legacy multi-pass flows).                     |
| `ProteinPositions()`, `Solvent()`, `Index()`, `Time()`, `HasRead()`, `error()` | Accessors.                         |

---

## H5 file layout

The file emitted after a successful `Run` has this top-level shape:

```
/                                        root attributes from TrajectoryProtein::WriteH5:
                                           protein_id, n_atoms, finalized
/atoms/
    element                              (N,) int
    residue_index                        (N,) uint64
    pdb_atom_name                        (N,) string

/trajectory/source/                      attributes: xtc_path, tpr_path, edr_path, configuration
/trajectory/frames/
    time_ps                              (T,) float64
    original_index                       (T,) uint64                         attr: n_frames
/trajectory/selections/<kind>/           one group per kind in selections_.Kinds()
    frame_idx                            (R,) size_t                         attr: n_records
    time_ps                              (R,) float64
    reason                               (R,) string

/trajectory/<result_name>/               one group per attached TR that implements WriteH5Group
    (TR-specific datasets + schema attributes)
```

Per-TR group names and attributes are owned by each TR's
`WriteH5Group`; downstream consumers parse them verbatim. There is
no central H5 schema module — the layout is the composition of what
each TR emits.

### Planned additions

Aspirational. Not yet in the tree.

- **`/header/` group.** File-level provenance: schema version,
  library version, run configuration name, source provenance from
  `ProteinBuildContext` (PDB source, build date, protonation tool,
  force field, assumptions). Root attributes currently emitted by
  `TrajectoryProtein::WriteH5` (`protein_id`, `n_atoms`, `finalized`)
  are candidates to migrate into this group when it lands.

---

END trajectory-scope replacement.
