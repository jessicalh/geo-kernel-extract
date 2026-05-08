# Legacy Amber Implementation Brief

Date: 2026-04-29

Status: active implementation brief for the first topology/charge landing.

Audience: external AI implementors and reviewers. This is the short
contract to read before touching code. The longer architecture record is
`openai-5.5-strong-architecture-layout.md`.

If this brief conflicts with looser transitional language in older
planning documents, this brief is the implementation rule for the first
landing.

## Goal

Land the smallest real object model that makes the existing Amber-based
calculator language explicit.

The migration is not a general topology platform. It is a disciplined
cleanup of the current model so future calculator and comparison work
does not keep rediscovering atom and residue identity through strings.

The required shape is:

```text
ProteinTopology
LegacyAmberTopology
CovalentTopology
ForceFieldChargeTable
CalculatorContract<LegacyAmberTopology, ChargeTableT>
```

`LegacyAmberTopology` is the concrete topology contract for the existing
calculators and the current NPY compatibility target.

## String Barrier

This is the first rule of the implementation.

```text
LegacyAmberTopology is the string barrier.
```

All loader/input/source naming, including CHARMM-style names, PDB atom
names, AMBER residue variants, compatibility display names, and
force-field atom names, is resolved before or during construction of
`LegacyAmberTopology` and its associated prepared state.

After that boundary:

```text
strings are display/provenance only
calculators do not perform string-based atom or residue identity lookup
calculators consume typed/indexed topology facts
```

If a calculator asks "what atom is this?" by inspecting a name string,
the barrier failed. If it asks `LegacyAmberTopology` for a typed/indexed
fact, the barrier is holding.

Any new calculator-side string-dispatch identity path is a fatal model
violation and grounds for rollback of that implementation session.

## Authority Map

After the first landing, authority should be readable from the object
names:

```text
Protein
    owns the protein instance and the shared atom/residue index space

ProteinTopology
    abstract base class for concrete topology contracts

LegacyAmberTopology
    owns the current Amber-derived calculator interpretation

CovalentTopology
    owns the named covalent bond graph component used by the topology

ForceFieldChargeTable
    owns loaded force-field partial charges and radii

ResidueChargeState / protonation_variant_index
    owns resolved protonation/variant state for this model

ChargeAssignmentResult
    temporary compatibility projection over ForceFieldChargeTable

ProtonationDetectionResult
    temporary compatibility projection/report over resolved residue state
```

Existing compatibility fields may remain, but they are projections:

```text
ConformationAtom::partial_charge
ConformationAtom::vdw_radius
Atom::bond_indices
Atom::parent_atom_index
Protein::RingAt / Protein::BondAt / Protein::Topology()
```

They must not become competing sources of truth.

The covalent-topology ownership rule is:

```text
CovalentTopology is not a ProteinTopology.
CovalentTopology is a value object owned by LegacyAmberTopology.
Protein::BondTopology() is a compatibility/read path to that component.
Protein::Topology(), if retained, is only an old synonym for BondTopology().
```

## Construction And Result Phases

The implementation must keep these phases distinct.

Preparation/construction may consume strings:

```text
input atom names
CHARMM/AMBER naming
residue variant names
force-field charge names
loader-specific source names
```

Preparation/construction produces typed/indexed state:

```text
LegacyAmberTopology
CovalentTopology
ResidueChargeState / protonation_variant_index
ForceFieldChargeTable
optional topology-owned IUPAC atom mapping hook
```

The calculator/result phase consumes prepared state:

```text
results may assign prepared topology/state facts into result/output surfaces
results may cache atom lists or arrays from the declared contract
results may not interpret names to decide atom or residue identity
```

There is no temporary class of "bad calculators" that is allowed to do
name lookup after topology construction.

## Thin Compatibility Results

`ChargeAssignmentResult` and `ProtonationDetectionResult` remain only as
thin compatibility adapters while the old result/dependency/output
surface is migrated.

While they remain `ConformationResult`s, they are still calculators in
the framework sense. Their special status is not that they are exempt
from the model. Their special status is that they are projections over
already prepared state.

Endpoint: once callers, outputs, and tests read the prepared state
directly, these compatibility results either go away or remain only as
boring reports. They must not accumulate new authority during the
transition.

### ChargeAssignmentResult

Allowed responsibilities:

```text
require ForceFieldChargeTable already exists
expose the table through the old result surface
copy/project charges and radii into old conformation atom fields
report assigned count, total charge, source, and related diagnostics
satisfy existing ConformationResult dependencies during migration
```

Forbidden responsibilities:

```text
own the authoritative loaded charge model
perform residue/atom string lookup after LegacyAmberTopology exists
silently invent zero charges for failed identity matches
become a second charge-assignment authority
```

### ProtonationDetectionResult

Allowed responsibilities:

```text
require resolved residue/protonation state already exists
report variant names/indexes through the old result surface
satisfy existing tests and dependency plumbing during migration
```

Forbidden responsibilities:

```text
independently decide residue identity after LegacyAmberTopology exists
perform post-construction string-dispatch protonation lookup
make protein validity depend on arbitrary calculator order
become a second protonation/variant authority
```

If either result still performs identity lookup from residue or atom name
strings during `Compute()`, the string barrier failed.

## CHARMM Translation

CHARMM translation is in scope because it is a correctness repair, not a
future-general topology feature.

The narrow rule is:

```text
CHARMM input identity -> explicit translation -> LegacyAmberTopology
calculators consume LegacyAmberTopology
```

Do not introduce `CHARMMTopology` as a peer runtime topology in this
landing. Do not make existing calculators portable between CHARMM and
Amber. Translate ingress identity into the explicit legacy Amber
contract, preserving source/provenance names for audit and output.

## IUPAC And Mutant Matching Door

Full IUPAC work is not part of this first landing, but the landing must
not close the door.

The shared protein atom index remains the identity spine:

```text
Protein atom index
LegacyAmberTopology atom index
future IUPAC/BMRB-facing atom identity
NPY/output rows
mutant comparison rows
```

The required door is:

```text
each concrete ProteinTopology may own an IUPAC atom-identity mapping
```

For now this can be a hook or stub. It must not become a generic
canonical-name layer. The typed identity fields are the future authority;
strings are for output/provenance.

Later mutant comparison should be able to match:

```text
residue position and residue identity context
backbone/common atoms by typed atom identity
sidechain atoms by residue-specific typed identity where meaningful
mutated sidechain atoms as unmatched/new/deleted by explicit policy
```

The first landing succeeds only if it keeps this later work possible
without adding calculator-local naming logic.

## Calculator Discipline

Existing calculators are not topology-portable. They bind to the
concrete topology they were written for:

```cpp
using Contract =
    CalculatorContract<LegacyAmberTopology, ChargeTableT>;
```

The contract is checked at the calculator entry point. After that, a
calculator may copy/cache atom lists, charge arrays, category vectors,
ring lists, or parent arrays for local work.

That local copying is acceptable. Local identity interpretation is not.

Rule:

```text
If a calculator needs identity, bond, ring, charge, or protonation data,
it gets that fact through the named topology/state object required by
its contract.
```

Do not add new direct string matching in calculator code. Do not add a
new generic `name`, `atom_name`, or `canonical_name` as internal
authority.

## First Landing Slice

The first implementation slice should create the small real architecture
landing, not the whole sweep.

Must exist:

```text
ProteinTopology
LegacyAmberTopology
CovalentTopology
ForceFieldChargeTable
CalculatorContract
```

Must be demonstrated:

```text
one real calculator or plumbing path obtains topology/charge through the
new contract

ChargeAssignmentResult and ProtonationDetectionResult remain old-surface
compatibility projections only

old output semantics are not intentionally changed

no new calculator-local atom/residue identity interpretation appears
```

Must not include:

```text
broad calculator sweep
general topology platform
CHARMM peer topology
full IUPAC naming framework
portable calculator redesign
intentional NPY drift unless agreed before the slice
```

## Acceptance Surface

NPY equivalence is the primary behavioral contract for this migration.

Default rule:

```text
existing NPY outputs should remain shape/value equivalent unless a
specific variation is proposed before implementation of that slice and
explicitly agreed
```

Practical gates:

```text
baseline gate:
    about 5 representative WT / non-mutated / mutant comparison proteins

full confidence gate:
    all 723 curated proteins before blessing the migration as non-breaking
```

Existing GTests, H5 checks, and logs are useful signals. They do not
override unexpected NPY drift.

## Working Protocol

This work may pass through bounded disruption. It must not pass through
confused disruption.

Each implementor session should be able to say:

```text
what slice was attempted
what object/model piece was put in place
what remains broken, if anything
whether the breakage is expected and bounded
whether NPY behavior is expected to remain equivalent
whether the string barrier still holds
```

Watcher review is part of the protocol. A watcher can say the slice
looks consistent with the agreed model, needs a bounded follow-up, or
should be rolled back. A passing partial test subset is not enough to
accept a slice that violates the object model.

End-of-session choices:

```text
continue
fix next
roll back
```

An accepted slice may become the new known-good point. A rejected slice
is fixed inside the bounded migration window or rolled back cleanly.

Rollback is the recovery mechanism if the model becomes incoherent. A
passing partial test subset must not be used to keep an incoherent
topology/charge model alive.

Hard rollback triggers:

```text
LegacyAmberTopology fails to hold the string barrier
new calculator-local string identity lookup appears
ChargeAssignmentResult becomes an independent charge authority
ProtonationDetectionResult becomes an independent identity authority
unexpected NPY drift is waved through without agreement
the slice turns into a general topology platform
the session ends with confused brokenness rather than bounded disruption
```

## Short Doctrine

```text
Small architecture, strict names.

Make Amber explicit.
Translate CHARMM into that explicit Amber contract.
Keep charge/protonation results as thin projections.
Use NPY equivalence as the end-of-window acceptance target.
Allow bounded disruption.
Reject confused disruption.
```
