# Pre-IUPAC Cruft Map

Date: 2026-04-29

Status: active planning document for the cleanup slice before topology-owned
IUPAC identity work.

Audience: external AI implementors and reviewers.

**2026-04-29 evening note on landing status.** The "first topology
landing" referenced below is **landed in the working tree but not yet
committed** (master HEAD remains `f8729e3`). The AMBER charge slice
(steps 1–6 in `spec/plan/amber-implementation-plan-2026-04-29.md`) sits
on top of that landing, also working-tree-only and uncommitted. Read
the central plan doc for the full working-tree state — this cruft map's
"the first topology landing put the basic objects in place" is correct
in the working tree but not yet on master.

This document controls the next cleanup slice. It exists because the first
`ProteinTopology` / `LegacyAmberTopology` landing put the basic objects in
place, but several old compatibility surfaces can still invite the wrong
model if a later session reads the code casually.

## Goal

Defend the current object model before adding direct IUPAC identity.

The first topology landing made the Amber-derived calculator language
explicit. This cleanup keeps that work from being diluted by old names,
silent fallbacks, or misplaced output concerns.

`ProteinTopology` is deliberately light. Its job is not to become a
large framework or a universal chemistry engine. Its job is to prevent
string assumptions from spreading past construction, and to give the
project a coherent place to name the typed/indexed facts that calculators
already rely on.

Because this project is developed across context boundaries, the code
must carry more of the model in its object names and ownership structure
than a single long-lived human session would need. If an implementation
requires remembering a conversation to understand where identity lives,
the model is not explicit enough.

The target model remains:

```text
Protein
    owns the protein instance and atom/residue index space

ProteinTopology
    abstract base for concrete topology contracts

LegacyAmberTopology
    the concrete topology contract for existing Amber-derived calculators
    and the string barrier for the current model

CovalentTopology
    the covalent bond graph component owned by LegacyAmberTopology

ForceFieldChargeTable
    the loaded force-field charge/radius table owned by Protein

ChargeAssignmentResult / ProtonationDetectionResult
    old result-surface projections over already prepared state
```

## Scope

This slice applies to:

```text
Protein
LegacyAmberTopology
ForceFieldChargeTable
ChargeSource / ff14SB charge loading
ChargeAssignmentResult as a compatibility projection
ProtonationDetectionResult as a compatibility projection
conformation-scope calculator contracts
```

This slice does not apply to:

```text
TrajectoryProtein
trajectory H5 /atoms export
trajectory result storage layout
full IUPAC atom identity implementation
mutant matching rewrite
calculator sweep
general topology platform
```

`TrajectoryProtein` may eventually export richer atom identity once the
protein-side model is stable. It is not the authority for this work and
must not be used to shape the topology model in this slice.

## Assumption Boundary

The project currently makes an Amber-derived assumption for the existing
calculator set. That assumption is allowed only when it is explicit:

```text
LegacyAmberTopology means the existing Amber-derived calculator language.
```

The goal is not to make calculators topology-portable. The goal is to
make their current topology language named, inspectable, and hard to
confuse with IUPAC, BMRB, CHARMM, display names, or trajectory export
names.

The goal is also not to make topology heavy. Add topology state only when
it removes an implicit string/name assumption, clarifies ownership, or
gives future code a typed/indexed fact it would otherwise rediscover
locally.

The recovered covalent-topology plan is:

```text
CovalentTopology is not a ProteinTopology.
CovalentTopology is owned by LegacyAmberTopology.
Protein::BondTopology() is a compatibility/read path to that component.
Protein::Topology(), if present, is only an old covalent-topology synonym.
```

The future IUPAC door is topology-owned symbolic identity:

```text
LegacyAmberTopology and future sister ProteinTopology objects may own
typed/symbolic atom-identity mappings internally.
Those mappings project to strings only for output, provenance, or
diagnostics.
This cleanup slice must not replace that future work with calculator-local
string matching.
```

The current shared identity spine remains:

```text
Protein atom index
LegacyAmberTopology atom index
ForceFieldChargeTable row
ConformationAtom row
NPY output row
future topology-owned IUPAC map row
```

## String Barrier Rule

`LegacyAmberTopology` is still the string barrier.

Preparation/construction may consume source names:

```text
PDB/mmCIF atom names
AMBER atom names
CHARMM atom names
AMBER residue variants
force-field charge table keys
loader-specific source names
```

After construction/preparation, calculators consume typed or indexed
facts. They do not rediscover residue or atom identity by name.

String use is acceptable in this slice only when it is clearly one of:

```text
loader/construction translation into typed/indexed state
force-field parameter-table lookup before ForceFieldChargeTable exists
display/provenance/output
test diagnostics
```

String use is not acceptable when it is:

```text
calculator-local atom identity logic
calculator-local residue identity logic
a new generic canonical-name layer
a workaround for missing topology facts
```

## Pinned Cruft Items

### 1. ff14SB Charge Helpers Live On A Result

`ChargeAssignmentResult` still exposes ff14SB parameter-file parsing and
variant-name mapping. That makes the compatibility result look like a
charge authority.

Required direction:

```text
trap ff14SB parameter parsing and Amber residue-name mapping inside
ParamFileChargeSource or an equivalent construction-boundary charge source
keep ChargeAssignmentResult as a projection over ForceFieldChargeTable
leave compatibility wrappers only if needed for tests or old call sites
```

Do not create a reusable Amber/string helper for this. A reusable helper
looks like a new identity service and weakens the barrier. The flat
ff14SB parameter file is string-keyed, so the source adapter may use
strings internally, but the only model object it should produce is the
prepared charge table.

### 2. Silent Charge Misses Are Not A Model

The flat ff14SB parameter-file path historically allowed missing atom
matches to become zero charge plus element-radius fallback. That may
preserve old output shape, but it must not be silent.

Required direction:

```text
ForceFieldChargeTable must be able to distinguish matched charge rows
from fallback rows
missing ff14SB entries must be reported explicitly
future policy can decide whether known misses are fatal
```

NPY equivalence matters. Therefore this cleanup may first make fallback
explicit without changing numeric output, then a later agreed slice can
turn selected misses into hard failures after the known cases are
classified.

### 3. `Protein::Topology()` Is A Naming Trap

`Protein::Topology()` historically meant the covalent bond graph. The new
topology model has `TopologyBase()`, `TopologyAs<T>()`, and
`LegacyAmber()`.

Required direction:

```text
new code must not call Protein::Topology()
use Protein::BondTopology() for the covalent bond graph
use Protein::LegacyAmber() for the concrete calculator topology
use Protein::TopologyBase() only when intentionally working with the ABC
```

If no source call sites require the compatibility alias, remove it. If it
must temporarily remain, mark it as deprecated and do not add new uses.

### 4. Construction-Time String Use Must Stay Quarantined

`Protein::ResolveProtonationStates`, `DetectAromaticRings`, and
`CacheResidueBackboneIndices` still read atom-name strings at construction
time. This is allowed because they translate source names into typed
indices and typed topology facts.

Required direction:

```text
do not copy those string checks into calculators
do not treat those functions as general naming utilities
future cleanup may move their outputs deeper into LegacyAmberTopology
```

### 5. Charge Projection Fields Are Not Authority

`ConformationAtom::partial_charge` and `ConformationAtom::vdw_radius`
remain compatibility projections consumed by existing APBS/Coulomb code.

Required direction:

```text
ForceFieldChargeTable is the owner of loaded force-field charges/radii
ChargeAssignmentResult copies values into old fields for compatibility
new charge-aware calculators should declare the charge-table contract
```

### 6. NamingRegistry Is Not The IUPAC Solution

The existing `NamingRegistry` contains useful evidence about CHARMM/IUPAC
gaps and old H5 naming problems. That evidence should inform the future
typed IUPAC map.

Required direction:

```text
do not expand NamingRegistry into an internal canonical-name authority
do not patch IUPAC semantics by adding calculator-local name translation
use NamingRegistry only at input/tool/output boundaries
```

Open clarity issue before code: `NamingRegistry.h` says unknown names
are rejected, but `TranslateAtomName()` currently returns the original
atom name unchanged when no rule exists. Decide explicitly whether this
is acceptable pass-through for same-convention names, or whether
boundary callers need a separate validation step / stricter translation
mode. Do not silently change this behavior inside an unrelated cleanup.

## Work Order

Do these before direct IUPAC implementation:

```text
1. Write this scope/cruft map and make it the entry point for the slice.
2. Trap ff14SB charge parsing/mapping authority inside the charge-source
   construction boundary, not on ChargeAssignmentResult and not in a
   reusable string helper.
3. Make flat-parameter charge misses explicit in the prepared charge table.
4. Remove or deprecate Protein::Topology() as a covalent-topology alias.
5. Rebuild and run focused tests.
6. Re-check that no source calculator gained string identity logic.
```

Do not do these in this slice:

```text
define IupacAtomId enums
populate IupacAtomMap
rewrite MutationDeltaResult matching
change TrajectoryProtein export
sweep calculators to new topology accessors
```

## Acceptance

A cleanup slice is accepted only if it can say:

```text
which cruft item was fixed
which compatibility surface remains and why
whether numeric NPY output is expected to stay equivalent
whether the string barrier still holds
whether any charge miss is explicit rather than silent
which tests were run
```

Rollback triggers:

```text
direct IUPAC work starts inside this cleanup slice
TrajectoryProtein becomes the authority for protein identity
ChargeAssignmentResult gains new charge authority
new calculator-local string identity lookup appears
unexpected NPY drift is waved through without agreement
```
