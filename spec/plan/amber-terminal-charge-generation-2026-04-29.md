# AMBER Terminal Charge Generation Note

Date: 2026-04-29

Status: active handoff note, before GROMACS/CHARMM PB-radius work.

Audience: future implementors and reviewers after context reset.

## Scope

This note is about AMBER-side charge/PB-radius correctness for the
existing static PDB/ORCA paths.

This is not the GROMACS/CHARMM radius work. Any observations about
GROMACS residue-variant mapping are parked until the AMBER terminal
charge problem is closed.

## Verified Current Fact

`data/ff14sb_params.dat` was regenerated from:

```text
AmberTools tleap
source leaprc.protein.ff14SB
set default PBRadii mbondi2
```

The regenerated table matched the checked-in table byte-for-byte. The
flat ff14SB data source is therefore reproducible and has clear
provenance for the rows it contains.

## Problem

The flat table contains explicit rows for:

```text
INTERNAL residues
NTERM residues supplied by the standard ff14SB terminal libraries
CTERM residues supplied by the standard ff14SB terminal libraries
```

AmberTools does not appear to supply every terminal form of every
protonation variant in the standard ff14SB libraries. In particular,
local inspection found internal variants for:

```text
ASH
CYM
GLH
LYN
```

but not standard terminal units named:

```text
NASH / CASH
NCYM / CCYM
NGLH / CGLH
NLYN / CLYN
```

This means terminal protonation/variant cases can be real chemistry
cases that the current generated flat table cannot represent.

## Hard Rule

Do not hand-code chemistry to fill these rows.

No extractor-local table of imagined terminal charges, radii, or atom
renames should be added. If the extractor has no authoritative prepared
source for a terminal/protonation case, it must fail clearly or route
through a documented generation path.

This is methods-facing. The charge/radius source must be something the
project can name and defend.

## Preferred Engineering Route

The best implementation is proper C++ library integration with an
authoritative AMBER topology reader.

There is no bias toward shelling out or text parsing when a usable C++
library is available. If the current local AmberTools installation lacks
headers or a clean development surface, getting or building the needed
library/source is part of the work. Do not downgrade the design merely
because the first installed package is incomplete.

The priority order is:

```text
1. Use a proper C++ AMBER topology reader/library to consume prepared
   PRMTOP data and expose charges, PB radii, residue/atom names, and
   provenance.
2. If the authoritative AMBER interface for generation is LEaP, invoke
   LEaP as the prepared-topology generator, then consume its PRMTOP
   through the C++ reader path.
3. Use direct documented PRMTOP parsing only as a bounded fallback if no
   practical C++ reader can be obtained in time.
4. Never replace missing chemistry with extractor-local constants or
   hard-coded charge/radius rows.
```

Verification byte comparisons are acceptable as tests. They are not the
model. The model is: external/prepared AMBER chemistry in, typed
`ForceFieldChargeTable` out.

## Architecture Decision

This is not a data-audit problem and not a flat-table completion
problem. It is an AMBER topology-preparation problem.

The extractor must have a path whose methods-facing description is:

```text
Protein/conformation was prepared with AmberTools LEaP under
leaprc.protein.ff14SB, then AMBER PRMTOP charges and PB radii were read
back onto the extractor's protein atom indices.
```

The flat `ParamFileChargeSource` remains a compatibility and regression
path for cases already represented in `data/ff14sb_params.dat`. It is
not the authority for awkward termini, disulfides, caps, or missing
force-field units.

The authority for awkward AMBER cases is a generated AMBER topology:

```text
typed Protein + ProteinConformation
    -> AmberTopologyPreparation boundary
    -> LEaP script/PDB/sequence/cap/disulfide decisions
    -> PRMTOP/INPCRD
    -> PRMTOP charge/radius reader
    -> ForceFieldChargeTable on Protein
    -> ChargeAssignmentResult projection on Conformation
```

AMBER-facing strings are allowed inside `AmberTopologyPreparation`.
They must not leak into calculators. The output of the boundary is a
typed/indexed charge table and provenance.

The provenance is part of the result, not commentary. Record at least:

```text
force field name/version family
tool and executable path used for generation
leaprc files sourced
PBRadii setting
AMBER library/frcmod/custom files used
LEaP input script path or retained script text
generated PRMTOP/INPCRD path when retained
preparation policy selected
per-chain terminal policy decisions
per-residue AMBER unit selected for preparation
cap residues inserted and which protein residues they attach to
explicit disulfide bond commands emitted
any atom-name mapping that was not one-to-one by exact name
```

The extractor should be able to explain whether an atom's charge came
from stock ff14SB terminal libraries, stock internal libraries, capped
fragment preparation, an existing PRMTOP, or a custom AMBER library.
That source accounting matters as much as the numeric charge.

## AMBER Practice To Encode

AMBER/LEaP supplies several distinct preparation mechanisms that must
stay distinct in this model:

```text
Standard protein termini:
    Use the N* and C* terminal library units supplied by
    aminont12.lib and aminoct12.lib through leaprc.protein.ff14SB.

Explicit sequence control:
    Use LEaP sequence construction when charges/radii can be generated
    from typed residue state alone. Use loadPdbUsingSeq only when
    coordinate ingestion through PDB is required. For charge/PB-radius
    assignment, sequence is usually cleaner because the extractor, not
    raw PDB residue strings, owns the typed residue decisions.

Histidine:
    Use HID, HIE, or HIP. HIS is an input convenience alias, not the
    internal scientific state.

Cysteine/disulfide:
    Free cysteine is CYS. Disulfide cysteine is CYX and requires an
    explicit LEaP bond command between the two SG atoms. Do not model a
    disulfide by bonding CYS sulfur atoms.

Terminal/capped fragments:
    AMBER supports capped peptide models such as ACE ... NME/NHE. This
    is the correct prepared-topology route when the intended model is a
    capped experimental fragment or when a terminal protonation variant
    must be represented as an internal AMBER residue by explicit policy.

Unsupported stock ff14SB terminal variants:
    Local ff14SB libraries provide internal ASH, CYM, GLH, and LYN, but
    not NASH/CASH, NCYM/CCYM, NGLH/CGLH, or NLYN/CLYN. ARN and TYM are
    not present in the local ff14SB amino12 library. These must not be
    silently mapped to ASP/CYS/GLU/LYS/ARG/TYR.
```

If a terminal unsupported variant is encountered, there are only three
valid outcomes:

```text
1. The preparation policy marks the chain end as capped/cut-fragment
   chemistry, so LEaP prepares the original residue as an internal
   variant with ACE/NME/NHE caps and the extractor maps only the
   original protein atoms back.
2. A custom AMBER library/frcmod route is provided and recorded in
   provenance, then LEaP prepares the topology from that route.
3. Preparation fails clearly with the missing AMBER unit and residue.
```

The first route is a modeling decision, not a fallback. It must be
configured or selected by a named policy so the methods text can say
what was done.

## Acceptable Generation Paths

If unsupported terminal/protonation cases occur and must be supported,
acceptable paths are prepared-toolchain paths:

```text
Use an existing PRMTOP when one is available and authoritative for the
same atom ordering.

Generate a PRMTOP through AmberTools/LEaP using a documented residue
library and parameter route, then read CHARGE and RADII from that PRMTOP.

Generate an extension to the flat table only from AmberTools-produced
artifacts, with the generator recording the exact source libraries,
commands, terminal state, residue name, atom names, charges, radii, and
any custom library/frcmod files used.
```

The implementation may cache generated rows, but the cache is not the
authority. The authority is the recorded AmberTools/prepared-topology
source.

The clean target shape is:

```text
Protein/Conformation input
    -> AMBER preparation/generation boundary when needed
    -> prepared PRMTOP/inpcrd or equivalent AMBER artifact
    -> C++ AMBER topology reader
    -> ForceFieldChargeTable
    -> ChargeAssignmentResult projection
```

This keeps the chemistry outside calculators and outside C++ constants,
while still keeping the extractor's internal model typed and inspectable.

## Immediate Implementation Slice

Implement the smallest vertical slice of `AmberTopologyPreparation`:

```text
Inputs:
    Protein
    ProteinConformation
    ForceField = Amber_ff14SB
    AmberPreparationPolicy

Policy:
    UseStockTermini
    UseCappedFragmentsForUnsupportedTerminalVariants
    FailOnUnsupportedTerminalVariants

Preparation:
    Write a LEaP input script in a temporary work directory.
    Source leaprc.protein.ff14SB.
    Set default PBRadii mbondi2.
    Select AMBER units from typed residue state, not raw PDB strings.
    Prefer LEaP sequence construction for charge/PB-radius topology.
    Use loadPdbUsingSeq only if coordinate ingestion is needed.
    Add explicit SG-SG bond commands for CYX disulfides.
    saveamberparm to PRMTOP/INPCRD.

Readback:
    Read CHARGE and RADII from the generated PRMTOP.
    Map generated AMBER atoms back to the extractor atom indices.
    Ignore cap atoms when caps were policy-selected.
    Fail on ambiguous or missing atom mapping.
```

The direct PRMTOP parser is acceptable for this slice because the
PRMTOP format fields needed here are CHARGE, RADII, ATOM_NAME,
RESIDUE_LABEL, and RESIDUE_POINTER. C++ cpptraj linkage remains a good
future replacement if the project accepts the dependency/licensing and
build work.

LEaP sequence probes on this machine succeeded for the intended first
coverage cases:

```text
sequence { NALA ALA CALA }
sequence { ACE ASH NME }
sequence { ACE CYM NME }
sequence { ACE GLH NME }
sequence { ACE LYN NME }
sequence { NCYX CCYX } plus bond mol.1.SG mol.2.SG
```

These are the first tests for `AmberTopologyPreparation`; they exercise
AMBER-supported stock termini, capped internal variants, and CYX
disulfide handling without adding extractor-local chemistry.

Do not spend another session proving the current flat table is
byte-identical. That question is answered. The next proof is that LEaP
topology preparation succeeds or fails correctly for the AMBER cases.

## Conversion Semantics

AMBER preparation is allowed because the target is explicitly AMBER
ff14SB charges and PB radii.

The later CHARMM/GROMACS slice must not silently change the meaning of
trajectory frames by converting them through AMBER. If an external
conversion path is used there, the conversion must be the best available
toolchain for coverage and must record exactly what changed:

```text
source force field
target force field or radius model
conversion tool and version
input topology/trajectory files
atom/residue mapping policy
names/variants changed during conversion
charges changed during conversion
radii model applied after conversion
any residues dropped, capped, patched, or retyped
```

For CHARMM trajectories, preserving source semantics is the default.
Conversion is a named modeling decision, not a plumbing detail.

## AMBER Practice Reference

AMBER practice for chemistry outside the standard libraries is to define
the missing residue/unit and parameters, load them into LEaP, check the
unit, and save a prepared topology. Relevant Amber documentation already
recorded in the current topology anchor:

```text
https://ambermd.org/Manuals.php
https://ambermd.org/tutorials/pengfei/index.php
https://ambermd.org/tutorials/basic/tutorial0/
https://ambermd.org/tutorials/basic/tutorial5/index.php
https://ambermd.org/tutorials/basic/tutorial4b/index.php
https://ambermd.org/antechamber/ac.html
```

The project should follow that shape: prepared topology from a
documented AMBER toolchain, not extractor-local chemistry. For ordinary
stock proteins this means LEaP plus stock libraries. For non-standard
chemistry this means custom AMBER library/frcmod inputs. For cut
fragments this can mean explicit ACE/NME/NHE capping when that is the
selected model.

## Implementation Shape

Keep this as charge/radius preparation state:

```text
ForceFieldChargeTable remains the prepared object on Protein.
ChargeAssignmentResult remains a projection over ForceFieldChargeTable.
ParamFileChargeSource remains the flat ff14SB adapter.
PrmtopChargeSource remains the authoritative prepared-topology adapter.
```

If a generated terminal-variant supplement is added, it belongs at the
charge-source/generator boundary, not in calculators.

The stronger implementation is a new preparation source rather than a
larger flat-file source:

```text
AmberTopologyPreparation
AmberPreparationPolicy
AmberPreparedTopology
AmberPreparedChargeSource or PrmtopChargeSource over generated PRMTOP
```

`AmberTopologyPreparation` owns all AMBER naming decisions. Existing
`Protein`, `Residue`, `AminoAcidType::variants`, terminal state, and
protein topology/disulfide information remain the typed inputs.

## Known Mapping Bug To Check In This Slice

The AMBER flat-table lookup must respect the shared
`AminoAcidType::variants` index contract.

Known candidate:

```text
CYS variant 0 = CYX
CYS variant 1 = CYM
```

The flat ff14SB lookup must not silently map `CYM` to default `CYS`.
Similar checks should cover `ARN` and `TYM` if those variants can be
produced by the current protonation pipeline.

## Parked GROMACS Notes

These are preserved so they are not lost, but they are not the active
slice. Do not start here until the AMBER terminal charge generation
work is closed.

Observed during the charge audit:

```text
FullSystemReader::VariantFromCharmmResidueName appears to use its own
0/1/2 convention for several residues instead of the shared
AminoAcidType::variants contract.

The shared contract is:
    default state = -1
    HIS: 0 HID, 1 HIE, 2 HIP
    ASP: 0 ASH
    GLU: 0 GLH
    CYS: 0 CYX, 1 CYM
    LYS: 0 LYN
    ARG: 0 ARN
    TYR: 0 TYM

Any CHARMM/GROMACS mapping must map into that contract, not a local
encoding.
```

A second GROMACS note:

```text
GromacsEnsembleLoader has a newer registry-key based mapper that is
closer to the shared contract, but it should later be checked for
variant names it may miss, including AMBER-form names such as ASH,
GLH, CYM, CYX, LYN, ARN, and TYM.
```

Charge-source note for later:

```text
The old GmxTprChargeSource shells through `gmx dump`; the current loader
paths already read TPR data through linked GROMACS C++ APIs and preload
charges. Future GROMACS/CHARMM work should prefer C++ API linkage where
available, not text dump parsing.
```

## Acceptance

Before returning to topology/IUPAC or moving to CHARMM/GROMACS:

```text
LEaP-generated PRMTOP preparation works for a normal stock protein case.
LEaP-generated PRMTOP preparation works for stock terminal HID/HIE/HIP
and CYS/CYX cases.
CYX disulfide preparation uses CYX plus explicit SG-SG bond commands.
Internal ASH/CYM/GLH/LYN are covered by capped-fragment policy tests.
Terminal unsupported variants fail clearly under FailOnUnsupported policy.
Terminal unsupported variants are handled through explicit capping policy
or custom AMBER library/frcmod policy, not canonical fallback.
The flat-table variant mapping respects the AminoAcidType variant-index
contract.
Existing AMBER charge tests pass.
NPY binary identity is checked before blessing any output drift.
```

The later corpus scan is still useful, but it is not the start of this
slice. It should run after this AMBER preparation boundary exists, so
exceptions are classified against real supported policies rather than
against the current flat table.
