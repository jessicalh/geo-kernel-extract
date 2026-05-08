# Current Topology Anchor

Date: 2026-04-29

Status: active context-recovery anchor.

Audience: external AI implementors and reviewers.

If a future session sees mixed planning history, read this first. This
document is the compact statement of the current topology decision.

## Core Decision

The active model is:

```text
Protein
    owns the protein instance and the shared atom/residue index space
    owns the attached ProteinTopology instance
    owns the loaded ForceFieldChargeTable when charges are prepared

ProteinTopology
    small abstract base class for concrete topology contracts

LegacyAmberTopology
    concrete ProteinTopology for the existing Amber-derived calculator
    language
    string barrier output for the current model
    owner of the covalent bond graph component

CovalentTopology
    value object owned by LegacyAmberTopology
    stores the covalent bond graph, bond categories, per-atom bond lists,
    and hydrogen parent mapping

ForceFieldChargeTable
    load-time force-field charge/radius table owned by Protein
    peer prepared state, not a ProteinTopology
```

The important recovered plan is:

```text
CovalentTopology is not a ProteinTopology.
CovalentTopology is the bond graph component inside LegacyAmberTopology.
Protein::BondTopology() is a compatibility/read path to that component.
Protein::Topology() is only an old covalent-topology synonym if it remains.
```

Do not shrink the model back to loose mappings. Do not widen it into a
general topology platform. The current target is a small ABC plus one
explicit concrete topology contract for the calculators that already
exist.

## IUPAC Door

The smaller topology plan still leaves room for topology-owned symbolic
IUPAC identity:

```text
LegacyAmberTopology, and future sister ProteinTopology objects, may own
their own typed/symbolic atom-identity mappings.
Those mappings are internal model facts, not canonical strings.
Strings are projected from those mappings only for output, provenance,
or diagnostics.
```

This slice does not define or populate IUPAC identity. It only keeps the
object model clear enough that IUPAC can later be added orthogonally
inside the relevant concrete topology instead of by calculator-local
string matching.

## In-Between State

The code may temporarily expose compatibility surfaces on `Protein` and
`Atom`:

```text
Protein::RingAt
Protein::BondAt
Protein::BondTopology
Protein::Topology
Atom::bond_indices
Atom::parent_atom_index
```

These are projections or access paths over prepared topology state. They
must not become competing sources of truth.

The desired ownership pressure is:

```text
LegacyAmberTopology owns the typed/indexed topology facts.
Protein exposes compatibility accessors while old calculators still need
them.
Calculators should move toward explicit LegacyAmberTopology contracts.
```

## String Barrier

`LegacyAmberTopology` is the string barrier.

Source names may be consumed while constructing prepared state:

```text
PDB/mmCIF atom names
CHARMM atom names
AMBER atom names
AMBER residue variant names
flat force-field parameter keys
loader-specific source labels
```

After construction/preparation:

```text
calculators consume typed/indexed topology facts
calculators do not rediscover atom or residue identity by string
strings are display, provenance, diagnostics, or boundary translation
```

The flat ff14SB parameter file is string-keyed. That string work belongs
inside `ParamFileChargeSource` or an equivalent construction-boundary
charge source. The model object that crosses the boundary is the prepared
`ForceFieldChargeTable`.

## Immediate Pre-IUPAC Work

Before direct IUPAC identity work, the cleanup pass should:

```text
make ff14SB charge/PB-radius misses explicit preparation failures
keep flat parameter parsing inside the charge-source boundary
keep ChargeAssignmentResult as a projection over ForceFieldChargeTable
preserve NPY binary identity unless drift is explicitly agreed
leave TrajectoryProtein out of scope
leave mutant matching and IUPAC enum work out of scope
```

## Active Interruption: AMBER Ends And APBS Radii

This work temporarily sits above the topology/IUPAC stack:

```text
1. Repair AMBER terminal/end-state handling for the flat ff14SB path.
2. Repair APBS radius semantics so the APBS bridge receives PB radii,
   not LJ radii described as APBS radii.
3. Quarantine current GROMACS/CHARMM PB-radius placeholders with an
   explicit non-authoritative status and TODO.
4. Return to the topology/IUPAC work only after the AMBER/APBS repair is
   understood and tested.
```

AMBER policy:

```text
ff14SB flat-table rows are generated from AmberTools tleap with
leaprc.protein.ff14SB and set default PBRadii mbondi2.
The table carries explicit INTERNAL, NTERM, and CTERM keys.
Missing AMBER charge/PB-radius rows are build/preparation failures.
There is no numeric element fallback in the AMBER path.
```

GROMACS/CHARMM policy for this interruption:

```text
TPR paths currently provide real charges but not a defensible PB-radius
set for APBS.
Those radii are quarantined as compatibility placeholders.
APBS may run with a warning on this path, but those outputs are not the
standard APBS result until CHARMM/GROMACS PB radii are implemented.
The CHARMM conversion is a must-do follow-up after the AMBER end-state
and APBS wrapping work lands.
The follow-up should start from the GROMACS API surfaces already used by
the project, not assume that text dumps or hand translation are the only
available route. If loaded GROMACS topology/type/radius data exposes the
needed force-field parameters, use that before inventing a separate path.
```

## Context Reset Handoff

If context resets before this interruption is closed, resume here:

```text
Do not return to IUPAC/topology migration first.
First finish the AMBER terminal/protonation charge audit and any needed
AmberTools/prepared-topology generation path.
Then finish the GROMACS/CHARMM PB-radius path, or prove the exact
source of truth needed to finish it.
Then return to the topology/IUPAC stack.
```

The current code state after the AMBER/APBS repair is:

```text
AMBER flat-table charges and PB radii are authoritative when a row exists.
AMBER flat-table misses are failures, not element fallback.
PRMTOP charges and radii are authoritative when CHARGE and RADII arrays
are present and complete.
GROMACS/CHARMM TPR charges are real, but PB radii are quarantined
compatibility placeholders.
APBS warns when non-authoritative radii are present.
```

The GROMACS starting points are:

```text
src/GromacsEnsembleLoader.cpp
    Reads TPR state/topology through the GROMACS API.
    Currently copies tpr_atoms.atom[ai].q into ForceFieldChargeTable.
    Currently assigns kCompatibilityPlaceholderPbRadiusAngstrom.

src/FullSystemReader.cpp
    Reads gmx_mtop_t through read_tpx_state.
    Already includes force-field/topology surfaces such as
    forcefieldparameters, idef, ifunc, and mtop utilities.
    This is the best first place to inspect whether the GROMACS API
    exposes enough type/parameter information for a defensible radius
    assignment.
```

The goal of the next GROMACS/CHARMM slice is:

```text
Every protein atom loaded from a CHARMM/GROMACS TPR has a real charge
and a defensibly sourced PB/APBS radius.
ForceFieldChargeTable rows on that path are marked authoritative.
The placeholder status disappears for solved paths.
APBS no longer warns for solved GROMACS/CHARMM inputs.
No calculator performs atom or residue identity discovery by string.
```

Do not silently convert CHARMM/GROMACS LJ parameters into APBS radii.
If the GROMACS API does not expose PB/APBS radii directly, decide and
document a defensible CHARMM-to-APBS radius source before coding the
translation. Likely candidates to research are CHARMM/PDB2PQR/APBS
radius tables or a generated reference table from a standard toolchain.
The source must be named in the data file or generator. A header-only
set of imagined radii is not acceptable.

This is charge/radius preparation state, not a new topology decision.
Keep the prepared object as `ForceFieldChargeTable`. Keep
`ChargeAssignmentResult` as an output/provenance projection over that
prepared state.

## AMBER Terminal Variant Audit

The current ff14SB table is generated from AmberTools with
`leaprc.protein.ff14SB` and `mbondi2`. That is correct for the rows it
contains, but AmberTools does not necessarily provide every terminal
form of every protonation variant in the standard libraries.

The known local shape is:

```text
amino12.lib contains internal ASH, CYM, GLH, LYN.
aminont12.lib/aminoct12.lib contain terminal forms for common standard
residues plus variants such as ASP/GLU/CYX/HID/HIE/HIP.
The standard terminal libraries do not appear to contain NASH/CASH,
NCYM/CCYM, NGLH/CGLH, or NLYN/CLYN.
```

This is not a harmless caveat. The project has enough real inputs to
measure the issue:

```text
Audit terminal residues across the high-resolution test residues and
the full curated protein set.
Count terminal residue names, protonation/variant names, source labels,
and pH/end-tag cases where available.
Report unsupported ff14SB terminal protonation states explicitly.
```

If such cases appear, do not fall back to internal residue rows. The
next decision is a preparation/calculator policy decision, for example:

```text
require a PRMTOP/authoritative prepared topology for that input,
generate a custom AMBER library/table through a documented AmberTools
route,
standardize the input chemistry before extraction,
or mark the charge/radius preparation unsupported for that case.
```

Before implementing custom handling, research how AMBER users normally
handle terminal/protonation variants or residues not supplied by the
standard force-field libraries. Capture the chosen source path in the
generator or data file. This is follow-up research after the reset
handoff is anchored, not an excuse to add code-local fallbacks.

Preliminary AMBER practice check:

```text
LEaP builds a complete topology from residue libraries plus parameter
files; residue libraries carry atom names, atom types, connectivity,
and charges, while parameter/frcmod files supply force-field terms.
If LEaP does not know the residue or atom names, the system is not ready
until that is fixed.
For non-standard residues, the Amber tutorials point users toward
creating a residue library/prepi/mol2 and any needed frcmod file, then
loading those into LEaP and saving a PRMTOP.
For amino-acid-like residues in chains, prepgen uses explicit head/tail,
main-chain, omitted-atom, neighboring-type, and net-charge declarations.
```

Sources to keep with this decision:

```text
https://ambermd.org/tutorials/pengfei/index.php
https://ambermd.org/tutorials/basic/tutorial5/index.php
https://ambermd.org/tutorials/basic/tutorial4b/index.php
https://ambermd.org/antechamber/ac.html
```

The project policy should follow that shape: use an authoritative
AmberTools-generated prepared topology/table, or cleanly fail and mark
the case unsupported. Do not invent terminal-variant charge/radius rows
inside extractor code.

## Acceptance Before Returning To IUPAC

At minimum, before the interruption is considered closed:

```text
cmake --build build -j2
./build/unit_tests
./build/structure_tests --gtest_filter=ChargeFF14SBTest.*:ApbsFF14SBTest.RealChargesProduceNonZeroFields:PrmtopChargeTest.*:OrcaRunTest.PrmtopChargesAreCorrect
./build/trajectory_tests --gtest_filter=FleetLoaderTest.ChargesReturned
```

For the GROMACS/CHARMM radius slice, acceptance also requires:

```text
the placeholder PB-radius warning is absent on solved paths,
the source of radii is visible in code/data provenance,
APBS sees positive PB radii for every protein atom,
and NPY binary identity is checked before blessing any drift.
```

If NPY binary identity changes, characterize the exact output drift
before blessing. The current expected drift class is APBS or
charge/radius-dependent output caused by replacing placeholders or
fixing terminal-state charge/radius assignment.

## Read Order

For the current migration slice:

```text
1. spec/plan/current-topology-anchor-2026-04-29.md
2. spec/plan/amber-terminal-charge-generation-2026-04-29.md
3. spec/plan/legacy-amber-implementation-brief-2026-04-29.md
4. spec/plan/pre-iupac-cruft-map-2026-04-29.md
5. spec/plan/amber-implementation-plan-2026-04-29.md
   ← TODAY'S CENTRAL CAPTURE-OF-DECISIONS. Locks O2/O3/O4,
   contains the crystal projection rule, the substrate-vs-
   conformation split, the drift policy, and the post-slice
   sequencing (PHASE 0 → PHASE 1 (N1.A-G) → PHASE 2 → OpenBabel
   exit → N4). Six implementation steps GREEN in the working
   tree, 62/62 tests passing, uncommitted.
6. spec/plan/openai-5.5-strong-architecture-layout.md
```

Older round-robin review documents explain how the decision was reached.
They are review history, not active implementation instructions.
