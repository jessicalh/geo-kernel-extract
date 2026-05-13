# LegacyAmberTopology Matrix Pre-Spec

Date: 2026-04-28

Status: setup document for the calculator topology/charge migration.

Related planning/design record:
`spec/plan/openai-5.5-strong-architecture-layout.md`

## Purpose

This work is too large to migrate by feel. Before implementation starts,
the project needs a calculator-by-calculator matrix that makes the
existing topology and charge dependencies explicit.

The matrix converts the migration into bounded sessions:

```text
enumerate current access
derive LegacyAmberTopology accessors
derive calculator contracts
record gotchas before code changes
pin the acceptance test for each calculator row
then migrate one bounded slice at a time
```

The matrix is not architecture commentary. It is a working artifact used
to decide what `ProteinTopology`, `LegacyAmberTopology`,
`ForceFieldChargeTable`, and the calculator contracts must expose.

## Round-Robin Rule

This pre-spec is intended for a round-robin review process:

```text
Reviewer A: Opus session
Reviewer B: Opus session
Reviewer C: OpenAI 5.5 session
```

Each stage should end with agreement on the artifact before the next
stage begins. No implementation sweep starts until the matrix and gotcha
list are accepted.

## Matrix Deliverables

The final matrix must produce four concrete outputs.

### 1. Exhaustive Accessor List

The union of all calculator reads becomes the candidate accessor set for
the concrete `LegacyAmberTopology` and associated charge entities.

The matrix should also confirm the small `ProteinTopology` ABC surface.
That base surface should stay minimal. It names the topology family; it
does not absorb every concrete topology accessor.

Current topology information is spread across:

```text
Protein
Residue
Ring
Bond
Atom
CovalentTopology
ConformationAtom
TrajectoryProtein / TrajectoryAtom
ConformationResult-owned structures
TrajectoryResult-owned structures
```

The matrix decides what the concrete topology contract exposes. It
should expose what calculators actually need, not an imagined complete
object model.

### 2. Per-Calculator Contract Declaration

Each row derives a contract from what the calculator reads.

Examples:

```cpp
using Contract =
    CalculatorContract<LegacyAmberTopology, NoChargeTable>;

using Contract =
    CalculatorContract<LegacyAmberTopology, ForceFieldChargeTable>;

using Contract =
    CalculatorContract<LegacyAmberTopology, MopacChargeTable>;
```

The contract is not a design choice made per calculator during
implementation. It is derived from the matrix.

### 3. Gotcha List

Each calculator row includes concrete migration gotchas. These are not
generic cautions. They are code-path facts that affect the migration.

Use the gotcha taxonomy below.

### 4. Acceptance Test

Every row names the test that proves that calculator migrated correctly.

Default acceptance:

```text
existing blessed NPY output for that calculator is unchanged
```

Trajectory default acceptance:

```text
existing trajectory NPY/H5 output is unchanged for the row's owned fields,
plus existing trajectory discipline tests still pass
```

Do not create a bespoke test plan for every row by default. Prefer the
existing NPY/H5/existing GTest surface. If a row has no meaningful
existing acceptance test, mark that fact in the matrix and carry it as a
known untestable spot.

## Matrix Columns

The matrix should use these columns. Additional local notes are fine, but
these are the required fields.

```text
Result
Scope
Source files
Production path
Reads Protein atom identity
Reads Residue symbolic slots
Reads residue variant/protonation state
Reads Ring topology
Reads Bond topology
Reads Atom bond_indices / parent_atom_index
Reads ConformationAtom enrichment projection
Reads ForceFieldChargeTable
Reads calculated charges
Reads downstream topology-derived structures
Writes topology-derived structures
Writes charge projection/cache fields
Post-construction NamingRegistry/cifpp/string identity use
Required topology concrete class
Derived contract
Needed ProteinTopology ABC accessors
Needed LegacyAmberTopology accessors
Needed charge-table accessors
Gotchas
Gotcha resolutions
Acceptance test
Migration band
Owner/reviewer
Status
```

## Accessor Buckets

Use these buckets while filling matrix cells. The final accessor list
should be the union of concrete needs found in the rows.

### ProteinTopology ABC

The base class is deliberately small. Candidate accessors:

```text
ProteinTopology::Kind()
ProteinTopology::Name()
ProteinTopology::AtomCount()
ProteinTopology::ResidueCount()
ProteinTopology::IupacOrNull()
```

Do not add ring, bond, residue-slot, or charge accessors to the ABC
unless the matrix shows a real cross-topology need. Existing calculators
should bind to `LegacyAmberTopology`, not to the base class.

### Protein Atom Identity

```text
AtomCount()
AtomAt(i)
Atom::element
Atom::residue_index
Atom::pdb_atom_name / legacy_amber_atom_name
```

### Residue Symbolic Topology

```text
ResidueCount()
ResidueAt(i)
Residue::type
Residue::atom_indices
Residue::sequence_number
Residue::chain_id
Residue::insertion_code
Residue::protonation_variant_index
Residue::N / CA / C / O / H / HA / CB
Residue::chi[0..3]
```

### Ring Topology

```text
RingCount()
RingAt(i)
Ring::atom_indices
Ring::type_index
Ring::parent_residue_index
Ring::parent_residue_number
Ring::fused_partner_index
Ring::Intensity()
Ring::LiteratureIntensity()
Ring::JBLobeOffset()
Ring::NitrogenCount()
Ring::Aromaticity()
Ring::RingSizeValue()
Ring::TypeName()
```

### Bond Topology

```text
BondCount()
BondAt(i)
Bonds()
Bond::atom_index_a
Bond::atom_index_b
Bond::order
Bond::category
Bond::is_rotatable
CovalentTopology::BondIndicesFor(atom_index)
CovalentTopology::HydrogenParentOf(atom_index)
Atom::bond_indices
Atom::parent_atom_index
```

### Legacy Enrichment Projection

```text
ConformationAtom::role
ConformationAtom::hybridisation
ConformationAtom::is_backbone
ConformationAtom::is_amide_H
ConformationAtom::is_alpha_H
ConformationAtom::is_methyl
ConformationAtom::is_aromatic_H
ConformationAtom::is_on_aromatic_residue
ConformationAtom::is_hbond_donor
ConformationAtom::is_hbond_acceptor
ConformationAtom::parent_is_sp2
```

### Force-Field Charge Table

```text
ForceFieldChargeTable::ChargeAt(atom_index)
ForceFieldChargeTable::RadiusAt(atom_index)
ForceFieldChargeTable::TotalCharge()
ForceFieldChargeTable::SourceForceField()
ForceFieldChargeTable::SourceDescription()
```

Current compatibility reads to map:

```text
ConformationAtom::partial_charge
ConformationAtom::vdw_radius
ChargeAssignmentResult::ChargeAt(atom_index)
ChargeAssignmentResult::RadiusAt(atom_index)
ChargeAssignmentResult::TotalCharge()
```

### Calculated Charges

```text
MopacResult::ChargeAt(atom_index)
ConformationAtom::mopac_charge
ConformationAtom::eeq_charge
ConformationAtom::aimnet2_charge
```

These derive future result-owned charge contracts:

```text
MopacChargeTable
EeqChargeTable
Aimnet2ChargeTable
```

Only introduce the common wrapper when a migration row actually benefits
from it.

### Downstream Topology-Derived Structures

Track these separately. They are not source topology, but their shape is
downstream of the topology contract.

```text
ProteinConformation::ring_geometries
ProteinConformation::ring_pairs
ProteinConformation::bond_lengths
ProteinConformation::bond_directions
ProteinConformation::bond_midpoints
ProteinConformation::rings_by_type
ProteinConformation::bonds_by_category
ProteinConformation::residues_by_type
ConformationAtom::ring_neighbours
ConformationAtom::bond_neighbours
ConformationAtom::mopac_bond_neighbours
```

If these structures need contract names, add them after the matrix shows
the need.

## Gotcha Taxonomy

Use these exact labels in the matrix.

```text
NONE
POST_CONSTRUCTION_STRING_IDENTITY
NAMING_REGISTRY_AFTER_CONSTRUCTION
CIFPP_AFTER_CONSTRUCTION
PROTONATION_VARIANT_UNSET_SEMANTICS
RING_TYPING_CHANGES_MATH
BOND_CATEGORY_CHANGES_MATH
CONFORMATION_TOPOLOGY_PROJECTION
TRAJECTORY_SCOPE_CROSS_READ
DOWNSTREAM_STRUCTURE_CONTRACT
MUTATION_MATCHING_IDENTITY
CHARGE_SOURCE_MATCHING
FORCE_FIELD_SPECIFIC_MATH
H5_OR_NPY_SCHEMA_TOUCH
NO_BLESSED_TEST
```

Definitions:

```text
POST_CONSTRUCTION_STRING_IDENTITY
    Reads atom/residue names after topology resolution to decide identity.

NAMING_REGISTRY_AFTER_CONSTRUCTION
    Calls NamingRegistry after construction for calculator identity.

CIFPP_AFTER_CONSTRUCTION
    Uses cifpp after construction for calculator identity.

PROTONATION_VARIANT_UNSET_SEMANTICS
    Depends on protonation_variant_index < 0 fallback behavior.

RING_TYPING_CHANGES_MATH
    Relies on ring type classes/properties where retagging changes output.

BOND_CATEGORY_CHANGES_MATH
    Relies on BondCategory values where recategorising changes output.

CONFORMATION_TOPOLOGY_PROJECTION
    Reads topology-like facts from ConformationAtom fields.

TRAJECTORY_SCOPE_CROSS_READ
    Reads conformation-scope outputs and trajectory-scope accumulators in
    one result or migration slice.

DOWNSTREAM_STRUCTURE_CONTRACT
    Reads or writes RingNeighbourhood, BondNeighbourhood, ring geometry,
    bond geometry, or similar topology-derived structures.

MUTATION_MATCHING_IDENTITY
    Compares atoms/proteins across WT/mutant or equivalent structures.

CHARGE_SOURCE_MATCHING
    Assumes force-field charges match atom order or legacy atom names.

FORCE_FIELD_SPECIFIC_MATH
    Uses parameters or categories that are not portable across force fields.

H5_OR_NPY_SCHEMA_TOUCH
    Migration changes feature output schema or array construction.

NO_BLESSED_TEST
    No existing test proves equivalence for this row.

KNOWN_UNTESTABLE_WINDOW
    Row can be migrated only inside a named partial-test window.
```

## Seed Calculator Inventory

This is the starting row inventory from `src/` on 2026-04-28. The matrix
fill pass must verify it against code, `OperationRunner`, and
`RunConfiguration`.

### Conformation Scope

| Result | Source files | Production path |
|---|---|---|
| `GeometryResult` | `GeometryResult.{h,cpp}` | Standard + trajectory per-frame |
| `SpatialIndexResult` | `SpatialIndexResult.{h,cpp}` | Standard + trajectory per-frame |
| `EnrichmentResult` | `EnrichmentResult.{h,cpp}` | Standard + trajectory per-frame |
| `DsspResult` | `DsspResult.{h,cpp}` | Standard + trajectory per-frame unless skipped |
| `ChargeAssignmentResult` | `ChargeAssignmentResult.{h,cpp}` | Standard when charge source present |
| `ProtonationDetectionResult` | `ProtonationDetectionResult.{h,cpp}` | Construction/protonation path |
| `MopacResult` | `MopacResult.{h,cpp}` | Full-fat / MOPAC-enabled path |
| `ApbsFieldResult` | `ApbsFieldResult.{h,cpp}` | Standard when charges present and APBS not skipped |
| `BiotSavartResult` | `BiotSavartResult.{h,cpp}` | Standard + trajectory per-frame |
| `HaighMallionResult` | `HaighMallionResult.{h,cpp}` | Standard + trajectory per-frame |
| `McConnellResult` | `McConnellResult.{h,cpp}` | Standard + trajectory per-frame |
| `RingSusceptibilityResult` | `RingSusceptibilityResult.{h,cpp}` | Standard + trajectory per-frame |
| `PiQuadrupoleResult` | `PiQuadrupoleResult.{h,cpp}` | Standard + trajectory per-frame |
| `DispersionResult` | `DispersionResult.{h,cpp}` | Standard + trajectory per-frame |
| `CoulombResult` | `CoulombResult.{h,cpp}` | Standard when charges present and vacuum Coulomb not skipped |
| `MopacCoulombResult` | `MopacCoulombResult.{h,cpp}` | MOPAC-enabled path |
| `MopacMcConnellResult` | `MopacMcConnellResult.{h,cpp}` | MOPAC-enabled path |
| `HBondResult` | `HBondResult.{h,cpp}` | Standard when DSSP present |
| `SasaResult` | `SasaResult.{h,cpp}` | Standard + trajectory per-frame |
| `EeqResult` | `EeqResult.{h,cpp}` | Standard + trajectory per-frame |
| `AIMNet2Result` | `AIMNet2Result.{h,cpp}` | When AIMNet2 model is loaded; mandatory in per-frame trajectory set |
| `WaterFieldResult` | `WaterFieldResult.{h,cpp}` | Explicit-solvent trajectory path |
| `HydrationShellResult` | `HydrationShellResult.{h,cpp}` | Explicit-solvent trajectory path |
| `HydrationGeometryResult` | `HydrationGeometryResult.{h,cpp}` | Explicit-solvent trajectory path |
| `GromacsEnergyResult` | `GromacsEnergyResult.{h,cpp}` | Trajectory frame energy path |
| `BondedEnergyResult` | `BondedEnergyResult.{h,cpp}` | Trajectory bonded-parameter path |
| `OrcaShieldingResult` | `OrcaShieldingResult.{h,cpp}` | ORCA comparison path |
| `MutationDeltaResult` | `MutationDeltaResult.{h,cpp}` | WT/mutant comparison path |
| `MolecularGraphResult` | `MolecularGraphResult.{h,cpp}` | Available result; verify production use |
| `DemoResult` | `DemoResult.h` / `DemoResult.cpp` | Demo/test only; verify whether to migrate |

### Trajectory Scope

| Result | Source files | Production path |
|---|---|---|
| `BsWelfordTrajectoryResult` | `BsWelfordTrajectoryResult.{h,cpp}` | All trajectory configurations |
| `BsShieldingTimeSeriesTrajectoryResult` | `BsShieldingTimeSeriesTrajectoryResult.{h,cpp}` | Per-frame extraction set |
| `BsAnomalousAtomMarkerTrajectoryResult` | `BsAnomalousAtomMarkerTrajectoryResult.{h,cpp}` | Per-frame extraction set |
| `BsT0AutocorrelationTrajectoryResult` | `BsT0AutocorrelationTrajectoryResult.{h,cpp}` | Per-frame extraction set |
| `BondLengthStatsTrajectoryResult` | `BondLengthStatsTrajectoryResult.{h,cpp}` | Per-frame extraction set |
| `PositionsTimeSeriesTrajectoryResult` | `PositionsTimeSeriesTrajectoryResult.{h,cpp}` | Per-frame extraction set |
| `ChiRotamerSelectionTrajectoryResult` | `ChiRotamerSelectionTrajectoryResult.{h,cpp}` | Available/pending selection path |

## Row Template

Use this template for each calculator row during the fill pass.

```text
Result:
Scope:
Source files:
Production path:

Reads Protein atom identity:
Reads Residue symbolic slots:
Reads residue variant/protonation state:
Reads Ring topology:
Reads Bond topology:
Reads Atom bond_indices / parent_atom_index:
Reads ConformationAtom enrichment projection:
Reads ForceFieldChargeTable:
Reads calculated charges:
Reads downstream topology-derived structures:

Writes topology-derived structures:
Writes charge projection/cache fields:
Post-construction NamingRegistry/cifpp/string identity use:

Derived contract:
Needed ProteinTopology ABC accessors:
Needed LegacyAmberTopology accessors:
Needed charge-table accessors:
Gotchas:
Gotcha resolutions:
Acceptance test:
Migration band:
Owner/reviewer:
Status:
```

Allowed matrix cell values:

```text
NO
READ
WRITE
READ_WRITE
VERIFY
```

Use `VERIFY` only for the first enumeration pass. Consolidation must turn
every `VERIFY` into a concrete value.

## Migration Bands

Assign each row one migration band after consolidation.

```text
BAND_A_CLEAN
    Reads only clear topology/charge accessors. No gotchas.

BAND_B_ACCESSOR_FIRST
    Needs a missing LegacyAmberTopology or charge-table accessor before the
    calculator can migrate cleanly.

BAND_C_GOTCHA_FIRST
    Has a named gotcha that must be resolved before migration.

BAND_D_SCHEMA_OR_TEST_FIRST
    Needs blessed NPY/H5 coverage or schema decision before migration.

BAND_E_DEFER
    Not part of this sweep.
```

## Acceptance Test Rules

Each row must name one acceptance test before implementation.

Preferred acceptance surfaces:

```text
existing blessed NPY equivalence
existing trajectory NPY/H5 equivalence
existing GTest result-specific equivalence
existing MutationDelta equivalence
existing source-boundary load test
```

If no test exists, the row gets:

```text
Gotcha: NO_BLESSED_TEST
Gotcha: KNOWN_UNTESTABLE_WINDOW if it must move before coverage exists
Migration band: BAND_D_SCHEMA_OR_TEST_FIRST or explicit inclusion in the
known untestable window
```

Adding new tests is allowed only when the row genuinely cannot be judged
through existing NPY/H5/GTest surfaces and the migration would otherwise
be blind in a way that matters for the sweep. The default response to a
gap is to name the gap, not to expand the work.

## Known Resolution Candidates

These are concrete policies surfaced during round-robin review. They are
not substitutes for the matrix rows; they are starting candidates for the
`Gotcha resolutions` field.

### `ChargeAssignmentResult`

Likely gotchas:

```text
CHARGE_SOURCE_MATCHING
POST_CONSTRUCTION_STRING_IDENTITY
```

Candidate resolution:

```text
ForceFieldChargeTable is built during Protein construction, or in the
loader before the loader returns. ChargeAssignmentResult exposes/projects
that table; it does not become a second owner.

On `(residue, legacy_amber_atom_name)` lookup miss against the selected
force-field charge table, fail loud with the offending residue/name pair.
Do not preserve the silent per-atom 0.0 charge fallback as normal
behavior.

Before flipping this discipline, Stage 0 audits blessed proteins for
current fallback hits.
```

### `ProtonationDetectionResult`

Likely gotchas:

```text
PROTONATION_VARIANT_UNSET_SEMANTICS
POST_CONSTRUCTION_STRING_IDENTITY
```

Candidate resolution:

```text
The matrix row must decide whether this remains a ConformationResult,
becomes a construction helper, or is removed after loaders populate
ResidueChargeState / protonation_variant_index directly.

The desired endpoint is no post-construction string-dispatch calculator
identity path for protonation state. Variant state should be resolved
before or during topology construction, then read as typed residue
state.
```

### Atom Name Field Schedule

Likely gotcha:

```text
POST_CONSTRUCTION_STRING_IDENTITY
H5_OR_NPY_SCHEMA_TOUCH
```

Candidate decision for Stage 4:

```text
Either rename `Atom::pdb_atom_name` to `legacy_amber_atom_name` before
the calculator sweep, or explicitly keep the old field name through the
sweep and rename after. Do not leave the schedule implicit.
```

## Stage Gates

### Stage 0: Pre-Flight Audit

Before enumeration/migration, identify current blessed-data dependence
on behavior that will become fail-loud or contract-explicit.

Minimum audit:

```text
Walk blessed proteins under tests/golden/.
For each, run the existing charge-loading path with diagnostics on the
silent fallback path.
Count atoms whose charge/radius currently lands via fallback rather than
an explicit force-field table hit.
```

Output:

```text
list of blessed proteins with non-zero silent-fallback hit count
per-hit disposition: fix input, document as expected, or defer row
pre-flight risk register for behavior that will flip from silent to loud
```

This does not require inventing new per-calculator tests. It establishes
whether existing blessed artifacts depend on behavior the migration is
about to remove.

### Stage 1: Enumeration

Parallel review is allowed. Each reviewer owns a disjoint subset of
results and fills row templates from source.

Output:

```text
one row per result
all access cells filled READ/WRITE/READ_WRITE/NO/VERIFY
gotchas named
initial gotcha resolutions proposed
initial acceptance test named or NO_BLESSED_TEST marked
```

### Stage 2: Consolidation

Resolve duplicated language, remove `VERIFY`, and produce the union
accessor list.

Output:

```text
complete matrix
ProteinTopology ABC accessor list
LegacyAmberTopology accessor union
ForceFieldChargeTable accessor union
calculated-charge contract list
gotcha list by result
gotcha-resolution list by result
```

### Stage 3: Round-Robin Agreement

Two Opus reviewers and one OpenAI 5.5 reviewer read the consolidated
matrix. Each reviewer signs off or names specific row changes.

Output:

```text
accepted matrix
accepted ProteinTopology ABC shape
accepted accessor list
accepted calculator contracts
accepted migration bands
accepted gotcha resolutions
```

### Stage 4: Implementation Planning

Turn the migration bands into bounded implementation sessions.

Output:

```text
LegacyAmberTopology landing order
ProteinTopology ABC landing order
ForceFieldChargeTable landing order
calculator migration order
known untestable window
tests that must stay green during the window
tests expected to return green at the end
explicit decision for `pdb_atom_name` rename timing
```

### Stage 5: Migration Sweep

Each calculator migration is one bounded feature slice where possible.
If a row expands beyond one session, stop and either add the missing
accessor/contract explicitly or move the row to a later band.

## Untestable Window

The matrix should name the untestable window before implementation
starts.

The goal is not "tests are green after every keystroke." The goal is:

```text
the partial-test window is known
the tests that remain meaningful are named
the return-to-green point is named
the rows inside the window are known
```

Recovery means rollback to the known-working state. The project already
proved this path by rolling back the failed naming/topology attempt on
2026-04-28. Tests and NPY equivalence are signals for accepting the new
state; they are not the only recovery mechanism.

If the migration becomes architecturally incoherent, rollback is the
decision. Passing partial tests must not be used to keep an incoherent
topology/charge model alive. The goal is not to nurse a bad migration
through the suite; it is to either land the model cleanly or return to
the known-working state.
