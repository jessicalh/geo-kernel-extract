# OpenAI 5.5 Strong Architecture Layout

Date: 2026-04-28

Status: planning/design record for the topology/charge contract cleanup.

This document captures the current agreed shape for making the
extractor's topology and charge contracts explicit. It is written for
external AI agents and for future implementation sessions.

## Core Decision

Existing calculators are not topology-portable. They were built against
the current residue naming, residue slot, ring typing, bond category,
hydrogen parent, and charge conventions. That contract should be named
directly:

```text
LegacyAmberTopology
```

`LegacyAmberTopology` is the explicit calculator topology contract for
the existing calculators and the current NPY compatibility target. It is
not "the biological topology" and not a universal protein topology.

`CovalentTopology` remains the covalent bond graph component. It should
not be treated as the whole topology.

Loaded force-field partial charges and radii should also become a named
entity:

```text
ForceFieldChargeTable
```

Calculated charges from MOPAC, EEQ, AIMNet2, or future calculators are
result-owned charge tables, not the same object as load-time force-field
charges.

IUPAC identity is not itself the topology. It is a named atom-identity
mapping/projection owned by a concrete topology.

## Modeling Boundary

This is a modeling boundary, not a hostile-input boundary. The extractor
is run by a small set of known users, against curated proteins, with
tests and human review enforcing release discipline.

The purpose of `LegacyAmberTopology` is to name the language the current
calculators speak. It is not to protect the code from arbitrary future
coders or unsupported molecule classes.

Free contract-preserving checks are welcome. Broad defensive layers that
add indirection without naming a real extractor failure mode are not the
goal of this design.

Once a calculator has explicitly bound to its topology and charge
contracts, it may copy/cache whatever atom lists, ring lists, charge
arrays, or category vectors it wants at the start of `Compute()`. Those
local copies are below the modeling concern. The contract matters at the
entry point.

## Name Set

Use these names in new C++ API and comments:

```text
LegacyAmberTopology
    Existing calculator topology contract.

CovalentTopology
    Geometry-derived covalent bond graph component.

ForceFieldChargeTable
    Load-time force-field partial charges and radii.

CalculatedChargeTable
    Result-owned calculated partial charges, when a common wrapper is useful.

NetCharge
    Whole-system integer or near-integer total charge used by MOPAC and
    loader/run configuration.

ResidueChargeState
    Protonation/variant state: HID, HIE, HIP, ASH, GLH, CYX, CYM, LYN,
    ARN, TYM, etc. This can initially remain represented by the existing
    `protonation_variant_index`.
```

Avoid these names in new architecture:

```text
FormalCharge
    Too ambiguous for this project and not what the extractor stores.

Standard
    Already overloaded between PDB-ish, AMBER-ish, IUPAC-ish, and tool
    translation language.

Topology
    Too generic for new public API. Existing `Protein::Topology()` can
    remain as compatibility, but new code should name the contract.
```

If multiple names for an atom are retained, keep them as named
projections by source or purpose:

```text
source_atom_name
legacy_amber_atom_name
iupac_atom_name
bmrb_atom_name
```

Do not introduce a new generic `name`, `atom_name`, or `canonical_name`
field as an internal authority.

## Result And Table Naming

The suffix `Result` already means something specific in this codebase:

```text
ConformationResult
    science/computation output attached to one ProteinConformation

TrajectoryResult
    science/computation output accumulated over a trajectory
```

Do not use `Result` for the charge entity.

Use `ChargeTable` for readable per-atom charge values:

```text
ForceFieldChargeTable
    loaded force-field partial charges and radii

MopacChargeTable
    MOPAC-derived Mulliken charges, if/when factored out

EeqChargeTable
    EEQ-derived charges, if/when factored out

Aimnet2ChargeTable
    AIMNet2 Hirshfeld charges, if/when factored out
```

The rule:

```text
*Result      = run-attached science output
*ChargeTable = readable per-atom charge values
```

`ChargeAssignmentResult` may remain as a compatibility/result attachment
because it participates in the existing `ConformationResult` pipeline.
The authoritative loaded charge entity should be `ForceFieldChargeTable`.

## IUPAC Mapping

IUPAC should not be introduced as a replacement topology.

The cleaner model is:

```text
each concrete ProteinTopology may own an IUPAC atom-identity mapping
```

For the legacy calculators, that means:

```cpp
class LegacyAmberTopology : public ProteinTopology {
public:
    const IupacAtomMap& Iupac() const;
};
```

The topology answers structural and calculator-language questions:

```text
What atoms exist?
Which residues own them?
Which bonds, rings, parents, roles, and categories exist?
Which calculator language is valid?
```

The IUPAC mapping answers naming and atom-identity questions:

```text
For atom index i, what is its IUPAC/NMR identity?
For residue r and IUPAC locant L, which atom index is that?
How should this atom align to BMRB/RefDB/output naming?
```

Suggested shape:

```cpp
class IupacAtomMap {
public:
    const IupacAtomId& Atom(size_t atom_index) const;

    std::optional<size_t> FindAtom(
        size_t residue_index,
        IupacLocant locant) const;
};

struct IupacAtomId {
    AminoAcid residue_type;
    IupacLocant locant;
    NmrClass nmr_class;
    RingAtomRole ring_atom_role;
    MethylGroup methyl_group;
    std::string iupac_name;
};
```

The string field is for output/provenance. The typed fields are the
identity.

This keeps IUPAC useful where it belongs:

- BMRB/RefDB matching
- exported atom identity
- per-residue atom-stratified analysis
- future calculators that explicitly declare an IUPAC-aware topology

It also prevents IUPAC from becoming a new global canonical-name layer.
The mapping is interpreted in the context of the concrete topology that
created it.

`LegacyAmberTopology` can own an `IupacAtomMap`. A future richer
`ProteinTopology` can also own an `IupacAtomMap`. Those maps may be built
with different structural knowledge, but each remains attached to the
topology contract that produced it.

The rule:

```text
IUPAC is an atom-identity projection inside each concrete ProteinTopology.
It is not a replacement topology and not a global canonical name layer.
```

## Future Topologies

The desired long-term shape is one topology family rooted at:

```text
ProteinTopology
```

`ProteinTopology` is an abstract base class. `LegacyAmberTopology` is a
concrete topology contract in that family, not a parallel fork. A future
richer topology may hold IUPAC-aware residue roles, BMRB alignment,
ring-atom roles, pro-R/pro-S distinctions, and other topology layers,
but it should still descend from `ProteinTopology`.

Calculators are written against a concrete topology contract from the
start. They are not meant to be fiddled into a different topology later.
For AI-written calculators, this is especially important: the topology
name is part of the calculator's language.

The rule is:

```text
all topologies descend from ProteinTopology
each calculator declares one concrete topology contract
old calculators declare LegacyAmberTopology
future calculators declare the richer concrete topology they were built for
```

Do not silently change the meaning of an existing calculator by changing
the ring types, bond categories, atom names, or residue slots supplied
under its existing topology name.

## Required ABC Object Model

The architecture should use an explicit abstract base class for topology:

```cpp
enum class ProteinTopologyKind {
    LegacyAmber,
    // Future concrete topology kinds go here.
};

class ProteinTopology {
public:
    virtual ~ProteinTopology() = default;

    virtual ProteinTopologyKind Kind() const = 0;
    virtual std::string_view Name() const = 0;
    virtual size_t AtomCount() const = 0;
    virtual size_t ResidueCount() const = 0;

    virtual const IupacAtomMap* IupacOrNull() const { return nullptr; }
};
```

The base class should stay small. It names the topology family and gives
shared lifecycle/identity hooks. It should not become a giant virtual
interface containing every possible ring, bond, residue, IUPAC, and
charge accessor.

Concrete topology classes expose their own language:

```cpp
class LegacyAmberTopology final : public ProteinTopology {
public:
    ProteinTopologyKind Kind() const override {
        return ProteinTopologyKind::LegacyAmber;
    }
    std::string_view Name() const override { return "LegacyAmberTopology"; }

    size_t AtomCount() const override;
    size_t ResidueCount() const override;

    const LegacyAmberResidueTopology& ResidueAt(size_t residue_index) const;
    const LegacyAmberAtomTopology& AtomAt(size_t atom_index) const;

    size_t RingCount() const;
    const Ring& RingAt(size_t ring_index) const;

    const CovalentTopology& Bonds() const;
    const IupacAtomMap& Iupac() const;
    const IupacAtomMap* IupacOrNull() const override { return &Iupac(); }
};
```

Calculators bind to the concrete topology class, not to the abstract
base:

```cpp
using Contract =
    CalculatorContract<LegacyAmberTopology, ForceFieldChargeTable>;
```

This keeps the model strict. The ABC gives the project one topology
family. The concrete topology gives each calculator the exact language it
was written against.

Do not use runtime polymorphism to make calculators topology-agnostic.
There is no planned switch-up where an old calculator is handed a richer
topology and expected to keep meaning the same thing. A calculator's
declared topology is part of its model.

`CovalentTopology` is not a `ProteinTopology`. It is a first-class value
object owned by a concrete `ProteinTopology`.

`ForceFieldChargeTable` is not a `ProteinTopology`. It is a first-class
load-time charge/radius object owned by `Protein` and consumed through a
calculator charge contract.

## LegacyAmberTopology

`LegacyAmberTopology` is the full current calculator language. It
includes the following.

### 1. Legacy Atom And Residue Naming Contract

After construction, the atom-name string used by old calculators is a
legacy calculator name, not arbitrary source text.

Current compatibility field:

```cpp
Atom::pdb_atom_name
```

Future clearer name, if introduced:

```cpp
legacy_amber_atom_name
```

The important point is not the field name. The important point is that
all string-dependent legacy sites read one named convention. That
convention is the one used by:

- `AminoAcidType` atom/ring/chi tables
- ff14SB residue/atom charge keys
- current residue backbone slot detection
- current ring detection
- current protonation variant detection

It is AMBER/PDB-derived calculator language. It is not BMRB/IUPAC atom
identity.

Loader-specific source names may still exist for provenance or output,
but they must not silently become calculator identity.

### 2. Residue Symbolic Topology

The contract includes the residue index language currently cached on
`Residue`:

```cpp
Residue::atom_indices
Residue::N
Residue::CA
Residue::C
Residue::O
Residue::H
Residue::HA
Residue::CB
Residue::chi[4]
Residue::protonation_variant_index
Residue::type
Residue::sequence_number
Residue::chain_id
Residue::insertion_code
```

The backbone and chi fields are atom indices, not strings. The string
work belongs at construction/resolution time.

`protonation_variant_index` remains load-bearing because it selects
histidine tautomer typing and ff14SB residue-name lookup. Its ordering
contract remains:

```text
HIS: 0=HID, 1=HIE, 2=HIP
ASP: 0=ASH
GLU: 0=GLH
CYS: 0=CYX, 1=CYM
LYS: 0=LYN
ARG: 0=ARN
TYR: 0=TYM
```

### 3. Ring Topology

Ring typing is calculator math. It is not cosmetic naming.

`LegacyAmberTopology` includes:

```cpp
Ring list
Ring atom indices
RingTypeIndex
Ring subclass/type object
parent_residue_index
parent_residue_number
fused_partner_index
```

It also includes the type-specific ring properties used by calculators:

```cpp
Intensity()
LiteratureIntensity()
JBLobeOffset()
NitrogenCount()
Aromaticity()
RingSizeValue()
TypeName()
```

Changing ring identity or ring typing changes the mathematical result.
That should be treated as a new topology contract or a new calculator,
not as the same calculator with different symbols.

### 4. Covalent Bond Topology

`CovalentTopology` remains a good object name for the bond graph.

It includes:

```cpp
std::vector<Bond> bonds
std::vector<std::vector<size_t>> bond_indices_by_atom
std::vector<size_t> hydrogen_parent_by_atom
```

Each `Bond` carries:

```cpp
atom_index_a
atom_index_b
BondOrder
BondCategory
is_rotatable
```

The current categories are part of the legacy contract:

```text
PeptideCO
PeptideCN
BackboneOther
SidechainCO
Aromatic
Disulfide
SidechainOther
Unknown
```

The key modeling point:

```text
CovalentTopology is the bond graph component inside LegacyAmberTopology.
It is not an atom identity authority.
```

It is geometry-derived, but it depends on legacy symbolic inputs:
residue slots and ring identities. The AMBER/name contamination is
upstream of `CovalentTopology`, not primarily inside it.

### 5. Legacy Enrichment Projection

The existing enrichment fields are derived, but calculators consume them
as part of the legacy language. They should be treated as a projection
owned by or explicitly derived from `LegacyAmberTopology`.

Current fields:

```cpp
AtomRole role
Hybridisation hybridisation
bool is_backbone
bool is_amide_H
bool is_alpha_H
bool is_methyl
bool is_aromatic_H
bool is_on_aromatic_residue
bool is_hbond_donor
bool is_hbond_acceptor
bool parent_is_sp2
```

The current implementation writes these onto `ConformationAtom` through
`EnrichmentResult`. That can remain as the compatibility/output surface.
The conceptual authority is the legacy topology contract plus the
geometric/bond information it resolves.

## ForceFieldChargeTable

Loaded partial charges and radii need a first-class readable object.
They should not exist only as incidental writes to
`ConformationAtom::partial_charge` and `ConformationAtom::vdw_radius`.

`ForceFieldChargeTable` is a clear object-model peer to topology, not a
field bundle and not a calculator side effect. It records the loaded
partial-charge/radius model for this protein instance.

Use:

```cpp
class ForceFieldChargeTable {
public:
    ForceField SourceForceField() const;
    ChargeModelKind Kind() const;
    const std::string& SourceDescription() const;

    size_t AtomCount() const;
    double ChargeAt(size_t atom_index) const;
    double RadiusAt(size_t atom_index) const;
    double TotalCharge() const;

    const std::vector<AtomChargeRadius>& Values() const;
};
```

Suggested source kind enum:

```cpp
enum class ChargeModelKind {
    AmberPrmtop,
    GromacsTpr,
    Ff14SBParamFile,
    Preloaded,
    Unknown
};
```

`ForceFieldChargeTable` should be owned by `Protein` because these charges
arrive with the loaded protein/topology.

Suggested API:

```cpp
class Protein {
public:
    const ProteinTopology& TopologyBase() const;

    template<class TopologyT>
    const TopologyT& TopologyAs() const;

    const LegacyAmberTopology& LegacyAmber() const;

    bool HasForceFieldCharges() const;
    const ForceFieldChargeTable& ForceFieldCharges() const;
};
```

`ChargeAssignmentResult` then becomes a projection/attachment result. It
does not become a second owner of the loaded charge model; it exposes the
`Protein`'s `ForceFieldChargeTable` and writes compatibility/output fields
onto `ConformationAtom`.

```cpp
class ChargeAssignmentResult : public ConformationResult {
public:
    const ForceFieldChargeTable& ChargeTable() const;
};
```

It may continue to write:

```cpp
ConformationAtom::partial_charge
ConformationAtom::vdw_radius
```

Those fields are useful output/cache columns. They are not the authority.

Do not call this `LegacyAmberChargeTable`. The topology contract may be
`LegacyAmberTopology`, but the loaded charges may come from AMBER prmtop,
ff14SB params, CHARMM36m TPR, or another force-field source that matches
atom order and loader construction.

## Calculated Charge Tables

Calculated charges are result-owned.

Current examples:

```cpp
MopacResult
    Mulliken charges

EeqResult
    EEQ charges

AIMNet2Result
    AIMNet2 Hirshfeld charges
```

These should not be folded into `ForceFieldChargeTable`. They are
conformation-dependent results.

If a common interface is useful later, use the concept name:

```cpp
CalculatedChargeTable
```

But the current result-specific APIs are acceptable:

```cpp
conf.Result<MopacResult>().ChargeAt(i)
conf.AtomAt(i).eeq_charge
conf.AtomAt(i).aimnet2_charge
```

Future cleanup can move each to:

```cpp
conf.Result<MopacResult>().ChargeTable().ChargeAt(i)
conf.Result<EeqResult>().ChargeTable().ChargeAt(i)
conf.Result<AIMNet2Result>().ChargeTable().ChargeAt(i)
```

That is a consistency cleanup, not a prerequisite for the topology fix.

## Protein API Convention

New code should make the topology contract visible.

Suggested ownership names:

```cpp
class Protein {
private:
    std::unique_ptr<ProteinTopology> protein_topology_;
    std::unique_ptr<ForceFieldChargeTable> force_field_charges_;

    // Temporary compatibility storage only, if needed during migration.
    std::unique_ptr<CovalentTopology> legacy_bond_topology_projection_;
};
```

Do not keep using a member named `topology_` for the covalent graph once
`ProteinTopology` lands. That name is the ambiguity this cleanup is
removing. If a standalone covalent graph member remains temporarily, call
it `legacy_bond_topology_projection_`, `bond_topology_`, or another name
that says it is the bond graph, not the protein topology.

Suggested shape:

```cpp
class Protein {
public:
    const ProteinTopology& TopologyBase() const;

    template<class TopologyT>
    const TopologyT& TopologyAs() const;

    const LegacyAmberTopology& LegacyAmber() const;

    bool HasForceFieldCharges() const;
    const ForceFieldChargeTable& ForceFieldCharges() const;

    // Compatibility projections for existing calculators.
    size_t RingCount() const;
    const Ring& RingAt(size_t i) const;

    size_t BondCount() const;
    const Bond& BondAt(size_t i) const;
    const CovalentTopology& BondTopology() const;

    // Existing compatibility name, if retained.
    const CovalentTopology& Topology() const;
};
```

`TopologyAs<TopologyT>()` is a contract assertion, not a negotiation
mechanism. If the stored topology is not `TopologyT`, it should fail with
a clear programmer-error diagnostic. Implementation can use a
kind-guarded cast or `dynamic_cast`; the important behavior is that a
calculator asking for the wrong concrete topology does not silently run.

Prefer this in new code:

```cpp
const auto& topo = conf.ProteinRef().TopologyAs<LegacyAmberTopology>();
```

This compatibility helper is also acceptable when the concrete contract
is already clear:

```cpp
const auto& topo = conf.ProteinRef().LegacyAmber();
```

Avoid this in new calculator code:

```cpp
const auto& topo = conf.ProteinRef().Topology();
```

`Protein::RingAt`, `Protein::BondAt`, `Atom::bond_indices`, and
`Atom::parent_atom_index` can remain as compatibility projections.
They are convenient and widely used. They should not obscure the
contract name in new calculator work.

## Calculator Contract Convention

Use narrow template enforcement. Do not template `Protein`,
`ProteinConformation`, `ConformationResult`, or `TrajectoryProtein`.

New calculators declare their contract:

```cpp
struct NoChargeTable {};

template<class TopologyT, class ChargeTableT = NoChargeTable>
struct CalculatorContract {
    using Topology = TopologyT;
    using ChargeTable = ChargeTableT;
};
```

Example:

```cpp
class CoulombResult : public ConformationResult {
public:
    using Contract =
        CalculatorContract<LegacyAmberTopology, ForceFieldChargeTable>;
};
```

A ring calculator that does not read charges:

```cpp
class BiotSavartResult : public ConformationResult {
public:
    using Contract =
        CalculatorContract<LegacyAmberTopology>;
};
```

A MOPAC charge calculator:

```cpp
class MopacCoulombResult : public ConformationResult {
public:
    using Contract =
        CalculatorContract<LegacyAmberTopology, MopacChargeTable>;
};
```

The two template axes are intentionally orthogonal:

```text
topology contract
charge-table contract
```

That orthogonality matters because this project is expected to need at
least one more concrete topology before the maths work is done. A future
maths-aware topology should be a new `ProteinTopology` implementation,
not a reason to blur charge access, result ownership, and topology
identity.

Then use helper functions to make the contract visible at call sites:

```cpp
template<class CalculatorT>
const auto& RequiredTopology(const ProteinConformation& conf) {
    using TopologyT = typename CalculatorT::Contract::Topology;

    static_assert(std::is_base_of_v<ProteinTopology, TopologyT>);
    return conf.ProteinRef().TopologyAs<TopologyT>();
}

template<class CalculatorT>
const auto& RequiredChargeTable(const ProteinConformation& conf) {
    using ChargeTableT = typename CalculatorT::Contract::ChargeTable;

    if constexpr (std::is_same_v<ChargeTableT, ForceFieldChargeTable>) {
        return conf.ProteinRef().ForceFieldCharges();
    } else if constexpr (std::is_same_v<ChargeTableT, MopacChargeTable>) {
        return conf.Result<MopacResult>().ChargeTable();
    } else {
        static_assert(!std::is_same_v<ChargeTableT, ChargeTableT>,
                      "No charge accessor registered for calculator contract");
    }
}
```

Inside a calculator:

```cpp
const auto& topo = RequiredTopology<CoulombResult>(conf);
const auto& charges = RequiredChargeTable<CoulombResult>(conf);
```

This is the intended enforcement boundary. The core object model remains
concrete; the template layer names calculator contracts and routes access
to the right entity.

## Old Bond Topology

The old bond topology becomes:

```text
CovalentTopology, the bond graph component of LegacyAmberTopology.
```

It should keep doing the thing it does well:

- detect covalent bonds from one coordinate set
- store bond list
- store per-atom bond indices
- assign hydrogen parents
- classify current legacy bond categories

It should not become a source of atom naming, residue identity, IUPAC
identity, or future topology semantics.

The public API should move toward:

```cpp
protein.LegacyAmber().BondTopology()
protein.BondTopology()
```

`protein.Topology()` can remain temporarily as an old synonym for the
bond graph, but new code should avoid that name.

## Construction Shape

Conceptual target:

```cpp
void Protein::FinalizeConstruction(const std::vector<Vec3>& positions,
                                    double bond_tolerance) {
    protein_topology_ = LegacyAmberTopology::Resolve(
        atoms_, residues_, positions, bond_tolerance);

    const auto& legacy = TopologyAs<LegacyAmberTopology>();

    // Compatibility projections.
    rings_ = legacy.ProjectRings();
    legacy_bond_topology_projection_ = legacy.ProjectBondTopology();

    for (size_t i = 0; i < atoms_.size(); ++i) {
        atoms_[i]->bond_indices =
            legacy.Bonds().BondIndicesFor(i);
        atoms_[i]->parent_atom_index =
            legacy.Bonds().HydrogenParentOf(i);
    }
}
```

The implementation does not need to land in this exact internal order.
The important thing is that `LegacyAmberTopology` becomes the named
authority, while existing fields remain projections for old code.

## Migration Shape

The migration sequence is:

1. Add `ForceFieldChargeTable` and make `ChargeAssignmentResult` expose
   it while continuing to populate `ConformationAtom::partial_charge`
   and `ConformationAtom::vdw_radius`.

2. Add `LegacyAmberTopology` as a named object around the current
   residue slot, ring, and `CovalentTopology` resolution.

3. Add `Protein::LegacyAmber()`, `Protein::BondTopology()`, and
   `Protein::ForceFieldCharges()`.

4. Keep `Protein::RingAt`, `Protein::BondAt`, `Atom::bond_indices`,
   and `Atom::parent_atom_index` as compatibility projections.

5. Add calculator contract declarations and helper accessors.

6. Migrate calculators opportunistically from direct compatibility
   fields to:

   ```cpp
   RequiredTopology<ThisCalculator>(conf)
   RequiredChargeTable<ThisCalculator>(conf)
   ```

7. Keep NPY output compatibility as the behavioral reference for the
   existing calculators.

The blessed NPY shape is the compatibility endpoint for old calculators.
If the old calculators are moved behind `LegacyAmberTopology` and the
NPYs remain equivalent, the migration has preserved their contract.

## Project-Specific Failure Modes

These are the concrete problems this design is meant to prevent.

### AMBER Name Leakage

`Atom::pdb_atom_name` has been acting as a pseudo-canonical identity.
That is acceptable only if the convention is named and enforced at the
loader/construction boundary. It is not acceptable as an unnamed general
atom identity.

### Charge Source Mismatch

Force-field charges are parallel to atom indices and/or legacy atom
names. A charge table is valid only for the protein/topology instance it
was loaded with. `ForceFieldChargeTable` makes that relationship explicit.

### Ring Typing Changes Calculator Math

Ring type changes alter parameters and therefore the result. A new ring
typing scheme is not just a new spelling of the same math.

### Bond Categories Are Legacy Semantics

`BondCategory::Aromatic`, `PeptideCO`, `PeptideCN`, etc. are the current
calculator categories. A future topology with different aromaticity or
bond semantics must not silently reuse those categories as if nothing
changed.

### ConformationAtom Fields Are Projections

`ConformationAtom::partial_charge`, role booleans, hybridisation, and
related fields are useful output/cache fields. They should not be the
only authority for topology or charge identity.

### Generic `Topology` Hides The Contract

`Protein::Topology()` currently returns the covalent bond graph. New code
should not use that name for the full calculator contract.

## Final Rule

Calculators declare a topology contract and, if needed, a charge-table
contract.

```text
LegacyAmberTopology + ForceFieldChargeTable
LegacyAmberTopology + MopacChargeTable
LegacyAmberTopology + NoChargeTable
```

All topology contracts should belong to the `ProteinTopology` family.
The object model stays concrete. The template layer exists only to make
each calculator's chosen topology and charge contracts explicit and
mechanically hard to mix.
