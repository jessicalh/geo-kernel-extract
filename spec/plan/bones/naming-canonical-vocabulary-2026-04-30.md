# Naming Canonical Vocabulary

Date: 2026-04-30

Status: active spec for the trajectory-load naming-canonicalization
discipline. Companion to `spec/plan/amber-implementation-plan-2026-04-29.md`
(the AMBER charge slice). Establishes the typed boundary by which
FF-port-prepared trajectories enter the project's typed Protein /
LegacyAmberTopology surface with canonical names and chemistry
recovered from upstream-authority files, not from any consume-side
inference or in-code chemistry tables.

Audience: external AI implementors and reviewers. Read this before
writing any naming-translation code in the trajectory loader path.

## Scope

This spec covers naming canonicalisation on the consume side of every
load path the project supports today:

```text
--orca / --mutant   (upstream tleap-prepared PRMTOP)
--pdb               (raw PDB → reduce → ff14SB)
--protonated-pdb    (user-attested protonated PDB → ff14SB)
--trajectory        (GROMACS-prepared MD; AMBER ff14SB only as of 2026-04-30)
```

Out of scope:

```text
CHARMM trajectories — the legacy 1A6J_5789 fleet fixture and its 5
    pre-existing FleetLoaderTest failures are out of scope. PHASE 0
    of the post-AMBER-slice queue replaces this fixture with AMBER
    trajectories (1P9J_5801 + 1Z9B_6577); CHARMM legacy is retired.

Force fields other than AMBER ff14SB. Future ports are typed-contract
    extensions; this spec does not commit to their convention details.

Modified residues, glycosylation, non-protein ligands, isopeptide
    crosslinks. The verdict-loud-fail path from the AMBER charge
    slice rejects these at the resolver boundary; canonicalisation
    is for accepted inputs only.
```

## Core principle

**Residue labels are advisory; atom-and-bond identity is the
load-bearing canonical fact.** Chemistry decisions (protonation
state, disulfide pairing, terminal capping) are made by the
upstream authority that prepared the input — tleap for ORCA,
GROMACS pdb2gmx for trajectories, our own AmberPreparedChargeSource
for the runtime-prepared variant. The consume side reads, walks,
audits, and applies; it never re-decides.

The convergence point: **PDB Chemical Component Dictionary (CCD)
is the FF-agnostic intermediate vocabulary**. Every modern
structural-biology toolkit (OpenFF Toolkit, ParmEd, gemmi, MDAnalysis,
RDKit) treats CCD-canonical names as the canonical surface; the
project does the same. CCD encodes: residue 3-letter codes (HIS,
CYS, ASP, ...), per-residue atom names per IUPAC 1979/1998
(N, CA, C, O, ND1, NE2, HD1, HE2, ...), bond lists, formal charges,
stereo configurations.

CCD does NOT enumerate FF-specific tautomer labels: there is no
HID / HIE / HIP / HISH / HSD / HSE / HSP CCD entry. Histidine is
a single CCD residue (HIS); the protonation tautomer is encoded by
which atoms are present, not by a name suffix.

## Per-path canonical authority

Each load path has exactly one upstream-authority file that records
the canonical chemistry decision. The loader reads that file; it
does not re-derive.

```text
--orca / --mutant
    Authority:   the upstream PRMTOP (path in OrcaRunFiles.prmtop_path)
    Format:      AMBER prmtop, written by AmberTools tleap.
    Names:       canonical AMBER (HID, HIE, HIP, CYX, ASH, GLH, LYN,
                 CYM, etc.) per ff14SB residue templates.
    Verification: not needed; tleap's output IS canonical AMBER by
                 construction. The project has run thousands of these
                 through the calibration corpus; format is trusted.

--pdb / --protonated-pdb (flat ff14SB table path)
    Authority:   the input PDB itself
    Format:      PDB-deposit canonical (HIS, CYS, ASP, etc.) per CCD.
    Names:       PDB CCD per atom (HD1, HE2, OXT, etc.).
    Verification: covered by the AMBER charge slice's flat-table
                 verdict (AnalyzeFlatTableCoverage).

--pdb / --protonated-pdb (AmberPreparedChargeSource path)
    Authority:   our generated PRMTOP from tleap
    Format:      AMBER prmtop produced by our AmberPreparedChargeSource
                 from typed Protein state.
    Names:       canonical AMBER per ff14SB residue templates.
    Verification: typed SG-SG distance (we own the chemistry decision
                 here per AmberTools tutorial methodology); recorded
                 in AmberPreparedChargeSource::Describe().

--trajectory
    Authority A: prep_run_*/topol.top "; residue N <name> rtp <rtp> q <q>"
                 comment lines.
                 GROMACS pdb2gmx writes one per residue; the rtp field
                 is canonical AMBER (HID / HIE / HIP / CYX / ASH / GLH /
                 LYN). Disulfide chemistry: CYS-paired residues show
                 ".name=CYS rtp=CYX" with HG absent and SG-SG bond
                 added.
    Authority B: prep_run_*/sources/structure_extracted_*.pdb (or
                 pdb_companion_*.pdb if structure_extract was not
                 run).
                 Deposit-canonical PDB; atom names per CCD/IUPAC for
                 atom-name canonicalisation (HB2/HB3 reference for
                 the amber14sb-port HB1/HB2 convention; CD1 reference
                 for ILE CD; etc.).
    Verification: cross-check the two authorities — every residue and
                 every heavy atom in topol.top has a counterpart in
                 the structure_extracted PDB; the per-residue chemistry
                 decision (rtp field) plus the heavy-atom topology
                 must be consistent. Logged.
```

## What GROMACS pdb2gmx actually does to the model

The trajectory path is the most complex; this is the load-bearing
detail:

```text
1. pdb2gmx -ff amber14sb reads input.pdb (deposit-canonical names).
2. Per-residue, looks up rtp template via aminoacids.rtp directly OR
   via aminoacids.r2b alias mapping (HISH→HIP rtp, HISD→HID rtp,
   ASPH→ASH rtp, LYSN→LYN rtp, CYS2→CYX rtp, etc.).
3. -his / -asp / -glu / -lys flags drive interactive tautomer/state
   selection. Hardcoded enums in pdb2gmx.cpp emit literal labels:
       HistidineStates → "HISD" "HISE" "HISH" "HIS1"
       AspartateStates → "ASP" "ASPH"
       (etc.)
   These are written into pdba->resinfo[].name in the SAME pass that
   sets pdba->resinfo[].rtp = chosen template.
4. specbond.cpp scans SG-SG distances. For pairs ≤ 0.2 nm, rewrites
   pdba->resinfo[].rtp = "CYS2" (which the .r2b table resolves to
   the CYX rtp template). pdba->resinfo[].name STAYS "CYS".
5. Topology is built from the .rtp template (canonical AMBER chemistry).
   Hydrogens generated per aminoacids.hdb (amber14sb HB1/HB2 numbering
   for β-methylenes; ILE CD; etc.).
6. topol.top written: per-residue comment line is
       ; residue N <.name> rtp <.rtp> q <q>
   The atom records carry .name. The TPR carries only .name (the .rtp
   comment is text-only).
```

Three distinct effects in one pass:

```text
A. Real chemistry change: atom-set delta.
   - Disulfide CYS loses HG (CYX has no sulfhydryl proton).
   - SG partial charge changes from -0.31 e (CYS) to -0.08 e (CYX).
   - Inter-residue SG-SG bond added to topology.
   - HIS gains HD1, HE2, or both per tautomer choice.

B. Real chemistry change: charge.
   - HIP carries +1; HID/HIE neutral.
   - ASP -1, ASH neutral. Etc.

C. Just labels: .name field divergence from CCD canonical.
   - HISH/HISD/HISE in .name (FF-port labels).
   - HB1/HB2 vs HB2/HB3 in atom records (amber14sb hdb numbering).
   - ILE CD vs CD1.
   - OXT vs OC1/OC2.
```

A and B are real chemistry made by GROMACS's specbond + rtp lookup.
The decisions are recorded in the topol.top rtp comment line.
C is a naming convention shift recorded in the .name field and
recoverable by structural alignment against the canonical reference.

## The walk procedure

One pass, per load:

```text
Step 1. Locate the canonical reference for this load.
    --orca / --mutant:        upstream PRMTOP at OrcaRunFiles.prmtop_path
    --pdb / --protonated-pdb: the input PDB
    --trajectory:             prep_run_*/topol.top (Authority A) and
                              sources/structure_extracted_*.pdb (Authority B)

    If a trajectory load doesn't have these on disk, the load fails
    with a typed error; we don't synthesise. The project owns the
    prep pipeline; if a future MD lacks the canonical reference,
    fix the prep pipeline. (User direction 2026-04-30: "There is
    no epistemic risk because we handle the steps.")

Step 2. Walk every residue.

    For each residue R in the FF-port form:
        - For trajectory loads: parse the topol.top "; residue N name
          rtp rtp" comment line. The rtp field is the canonical AMBER
          name (HID / HIE / HIP / CYX / ASH / GLH / LYN / CYM / etc.).
          The .name field is the FF-port label.
        - For non-trajectory loads: the input file's residue name is
          the canonical authority directly.

        Emit:
            R.canonical_residue_name = (rtp from topol.top, or .name
                                        from PRMTOP/PDB)
            R.ff_port_label          = (.name from topol.top, or
                                        same-as-canonical for other paths)
            R.protonation_variant    = derived from canonical name +
                                        atom presence (HIP / HID /
                                        HIE / CYX / CYS / CYM / ASH /
                                        GLH / LYN / TYM / ARN / etc.)

        If trajectory: cross-check against
        sources/structure_extracted_*.pdb at the same residue index.
        Log discrepancies.

Step 3. Walk every atom of every residue.

    For each atom A in R (FF-port form):
        - Match A against the corresponding atom in the canonical
          reference (Authority B for trajectories; same-file for
          others) by:
              (residue index, element, bonded-neighbor signature)
          Heavy atoms first; hydrogens by typed parent + AminoAcidType
          template.
        - Three outcomes:
            MATCHED     A.canonical_name = reference's name; if it
                        differs from A.ff_port_name, log the rename.
            ADDED       A is in the FF-port form but not in the
                        reference. Typed reason: protonation H added
                        by pdb2gmx (HD1/HE2/HG/etc.). Resolved via
                        AminoAcidType template; canonical_name comes
                        from the template's atom slot.
            STRIPPED    A is in the reference but not in the FF-port
                        form. Typed reason: chemistry change during
                        prep (HG missing on disulfide-paired CYS,
                        HE2 missing on HID, etc.). Confirmed against
                        the canonical residue rtp.

        After-walk: every FF-port atom has a canonical name. The
        protonation_variant_index is set on R from the typed
        AminoAcidType variants table.

Step 4. Chemistry decisions are READ, not made.

    The disulfide bond between two CYS residues (say residue 8 and
    residue 22 in 1P9J_5801) is recognised because:
        - topol.top rtp field for both is "CYX"
        - both have HG absent
        - the [ bonds ] section of topol.top contains an SG-SG pair

    The loader does NOT measure SG-SG distance to decide. GROMACS
    already decided via specbond.dat at prep time; the rtp record
    is GROMACS's typed self-report of that decision.

    Same for protonation: rtp HIP at residue 4 is GROMACS's record
    of the user's interactive choice (or geometric classifier
    decision). The loader applies it; it does not re-derive.

    The AmberPreparedChargeSource path is the ONE place we make
    chemistry decisions on the consume side, and that's because
    we ARE the upstream authority on that path (we run tleap from
    typed Protein state). The disulfide-distance methodology there
    is per AmberTools tutorial citation and is documented in
    spec/plan/amber-implementation-plan-2026-04-29.md.

Step 5. Emit the canonicalisation log.

    Per-load JSON record dropped alongside other run logs in the
    prep dir (or wherever the existing per-run log convention
    places them — see "Log location" below).

    Schema (see "Canonicalisation log schema" below).

Step 6. Apply canonical names.

    The typed Protein has every residue and every atom labeled with
    its canonical name. ProteinBuildContext carries a pointer to
    the canonicalisation log path. Calculators read canonical names;
    output writers project from them.
```

## Water and ion canonicalisation

Water and counter-ions are non-protein chemistry that lives alongside
the protein in solvated MD trajectories. They need the same
canonicalisation discipline, but the structural-match step degenerates:
each water molecule is uniform (3 atoms in TIP3P/SPC; 4 atoms in
TIP4P-family; 5 in TIP5P), and ions are single atoms. Per-molecule
walks are unnecessary; an aggregate canonicalisation entry per
species is sufficient.

### Authority for water and ion chemistry

Authority for water: the `#include "<ff>.ff/<watermodel>.itp"` line
in `prep_run_*/topol.top`, plus the included .itp file itself.
The .itp declares the moleculetype (residue label, per-atom names,
charges, masses, bonds, settles constraints).

For the project's current AMBER fixtures (1P9J_5801, 1Z9B_6577):

```text
File:           /home/jessica/gromacs/src/gromacs-2026.0/share/top/amber14sb.ff/tip3p.itp
Residue label:  SOL
Atom names:     OW (oxygen), HW1, HW2 (two hydrogens)
Atom types:     OW_tip3p, HW_tip3p
Charges:        OW = -0.834e, HW1 = +0.417e, HW2 = +0.417e (canonical TIP3P)
Geometry:       Settles constraint (rigid; OH bond 0.09572 nm; HOH angle 104.52°)
```

Authority for ions: the `aminoacids.rtp` (or a per-ion .itp file
included by topol.top) declares each ion species. For monovalent
ions in our fixtures (`NA`, `K`, `CL`), GROMACS uses the PDB CCD
3-letter codes directly as both residue label and atom name.

### Canonicalisation translation

Water (TIP3P, amber14sb-port):

```text
GROMACS              PDB CCD (HOH)             Source
SOL    (residue)  →  HOH                       amber14sb.ff/tip3p.itp
                                                 + xlateat.dat (HOH O→OW)
OW     (atom)     →  O                          (same)
HW1    (atom)     →  H1                         (same)
HW2    (atom)     →  H2                         (same)
```

Other water models are minor extensions of this table, citable to
their respective .itp files in the same `<ff>.ff/` directory:

```text
TIP4P / TIP4P-Ew:  +1 atom MW (virtual M-site)  →  CCD has no MW; emit as
                                                    "MW" with virtual_site=true
                                                    flag in canonicalisation log
TIP5P:             +2 atoms LP1, LP2 (virtual)  →  same shape
SPC / SPC/E:       same OW/HW1/HW2 names         →  same translation as TIP3P
OPC / OPC3:        same shape as TIP3P / TIP4P  →  same approach
```

Ions:

```text
GROMACS              PDB CCD                    Source
NA     (residue   →  NA  (sodium ion)           CCD (already canonical)
       and atom)
K      (same)     →  K   (potassium ion)        CCD (already canonical)
CL     (same)     →  CL  (chloride ion)         CCD (already canonical)
```

For the ions in our current fixtures, no translation is required.
GROMACS amber14sb writes CCD-canonical ion names directly. Future
ions (Mg, Ca, Zn, Mn, etc.) need verification when first encountered;
most are CCD-canonical but a few may have AMBER-port specifics
(e.g., MG vs MG2, CA vs CA2 — note the residue collision with
calcium-the-ion vs CA-the-protein-Cα).

### Virtual sites (flagged for future, NOT in current scope)

TIP4P-family water models have a 4th atom (M-site, MW, or EP) that
carries charge but has no mass; its position is constrained
geometrically from O, H1, H2. Same for TIP5P (two lone-pair virtual
sites L, LP1/LP2). For E-field, dipole, and APBS calculations,
calculators must know:

- Which atoms are virtual sites (no mass, no thermal motion).
- That charge lives on the M-site, not O, in TIP4P models.
- That virtual-site positions are deterministic from the real atoms.

Our current trajectories (1P9J, 1Z9B) are TIP3P (3-atom water,
no virtual sites). When TIP4P+ trajectories enter scope, the
canonicalisation log gets a per-water-model `virtual_sites: [...]`
entry naming the M-site name and its constrained-from atoms. The
typed Protein gains a per-atom `is_virtual_site: bool` field.

This is documented here so it isn't lost when needed; not
implemented now.

### Charges and dipoles

Water atom charges are read directly from the TPR via libgromacs
(already done by `FullSystemReader::ReadTopology`; the
`tpr_atoms.atom[ai].q` field gives the per-atom partial charge,
which for waters comes from the .itp moleculetype). For 1P9J and
1Z9B these are the canonical TIP3P values (-0.834, +0.417, +0.417);
recorded in the `Topology::water_O_q` / `water_H_q` summary on the
existing `FullSystemReader` output.

Calculators that compute water dipoles, polarization, water-field
contributions, hydration-shell geometry, etc. read positions and
charges from typed Protein/ProteinConformation. The naming
canonicalisation ensures these calculators see canonical CCD names
(`HOH`, `O`, `H1`, `H2`) regardless of the FF-port labels in the
input. Existing water-field / hydration calculators that match by
atom name on the input side need a one-time review: if they read
`OW` literally (FF-port label), they need to either (a) read the
canonical name `O` post-canonicalisation, or (b) be updated as
part of the Phase 2 calculator-by-calculator migration to consume
typed atom-role flags rather than name strings. The user's
direction on Phase 2 (substrate first, calculators later, on each
calculator's own schedule) governs this.

### Aggregate canonicalisation log entries

Per-water-molecule log entries would be enormous (23933 × per-residue
JSON entries for 1P9J alone). Water and ions get aggregate
canonicalisation log entries instead:

```jsonc
{
  "non_protein_canonicalisation": {
    "water": {
      "model": "TIP3P",
      "model_source": "prep_run_*/topol.top:#include amber14sb.ff/tip3p.itp",
      "ff_port_residue_label": "SOL",
      "ccd_residue_label": "HOH",
      "molecule_count": 23933,
      "atoms_per_molecule": 3,
      "virtual_sites": [],
      "atom_renames": [
        {"ff_port": "OW",  "canonical": "O",  "reason": "TIP3P GROMACS amber14sb-port convention"},
        {"ff_port": "HW1", "canonical": "H1", "reason": "(same)"},
        {"ff_port": "HW2", "canonical": "H2", "reason": "(same)"}
      ],
      "charges": {"O": -0.834, "H1": 0.417, "H2": 0.417},
      "geometry_constraint": "settles (rigid OH=0.09572nm, HOH=104.52deg)"
    },
    "ions": [
      {
        "ff_port_residue_label": "NA",
        "ccd_residue_label": "NA",
        "molecule_count": 28,
        "rename_required": false,
        "charge": +1.0
      },
      {
        "ff_port_residue_label": "CL",
        "ccd_residue_label": "CL",
        "molecule_count": 27,
        "rename_required": false,
        "charge": -1.0
      }
    ]
  }
}
```

This sits alongside the per-residue protein canonicalisation array
in the same JSON log file.

## Canonicalisation log schema

JSON. One file per load. Schema:

```jsonc
{
  "schema_version": "1",
  "created": "2026-04-30T...",
  "protein_id": "1P9J_5801",
  "load_path": "trajectory",          // or "orca", "mutant", "pdb",
                                      //   "protonated-pdb",
                                      //   "amber-prepared"
  "authority_files": {
    "topol_top": "prep_run_*/topol.top",
    "reference_pdb": "sources/structure_extracted_*.pdb",
    "tpr": "prep_run_*/production_*/production.tpr",
    "prmtop": null                    // or path for orca/mutant
  },
  "residues": [
    {
      "index": 4,
      "sequence_number": 4,
      "chain_id": "A",
      "ff_port_label": "HISH",
      "canonical_residue_name": "HIP",
      "ccd_residue_name": "HIS",
      "protonation_variant": "doubly",
      "protonation_variant_index": 2,
      "rtp_source_line": "prep_run_*/topol.top:77",
      "atom_renames": [],
      "atom_added": ["HD1", "HE2"],
      "atom_stripped": [],
      "decision_summary": "HIS at residue 4 received HD1 and HE2 via pdb2gmx -his choice 2 (HISH); rtp HIP confirmed in topol.top; charge +1.0"
    },
    {
      "index": 8,
      "sequence_number": 8,
      "chain_id": "A",
      "ff_port_label": "CYS",
      "canonical_residue_name": "CYX",
      "ccd_residue_name": "CYS",
      "protonation_variant": "disulfide",
      "protonation_variant_index": 0,
      "rtp_source_line": "prep_run_*/topol.top:185",
      "atom_renames": [
        {"ff_port": "HB1", "canonical": "HB2", "reason": "amber14sb hdb β-methylene numbering"},
        {"ff_port": "HB2", "canonical": "HB3", "reason": "amber14sb hdb β-methylene numbering"}
      ],
      "atom_added": [],
      "atom_stripped": ["HG"],
      "atom_stripped_reason": "disulfide formation; specbond.dat geometric SG-SG ≤ 0.2 nm at prep time; CYS rtp rewritten to CYX",
      "disulfide_partner_residue_index": 22,
      "decision_summary": "CYS at residue 8 paired in disulfide with residue 22 via GROMACS specbond at prep time; HG absent; SG charge -0.08e per CYX rtp"
    }
    // ... one entry per residue
  ],
  "summary": {
    "residues_total": 54,
    "residues_with_canonical_renames": 3,
    "residues_with_atom_changes": 9,
    "disulfide_pairs": 3,
    "histidine_protonation_distribution": {"HID": 0, "HIE": 0, "HIP": 3},
    "atoms_total": 846,
    "atoms_renamed_label_only": 142,
    "atoms_added_during_prep": 6,
    "atoms_stripped_during_prep": 6,
    "errors": []
  }
}
```

## Log location

Per user direction 2026-04-30: logs go where other per-run logs go,
findable but not cruft. Concretely, the canonicalisation log is
dropped at:

```text
trajectory:        prep_run_*/canonicalization_record.json
                   (alongside topol.top, mdout.mdp, step_*.stdout)
orca / mutant:     <orca_run_dir>/canonicalization_record.json
pdb / protonated:  <output_dir>/canonicalization_record.json
amber-prepared:    <work_dir>/canonicalization_record.json
                   (alongside the generated leap.in / prep.prmtop /
                    tleap.log already written by AmberPreparedChargeSource)
```

`ProteinBuildContext` gains a typed field `canonicalization_log_path`
(string) for downstream readers (methods text generators, H5
attribute writers, audit reviewers).

## What we do NOT do

Explicit non-decisions, with rationale:

```text
1. We do NOT measure SG-SG distances on the consume side for
   trajectory loads. GROMACS's specbond.dat already decided at prep
   time; the rtp comment in topol.top is the typed record. Distance
   measurement on consume would be a second authority on the same
   chemistry decision (PATTERNS.md anti-pattern: do not be a second
   authority on already-decided typed state).

   Exception: AmberPreparedChargeSource (the runtime-tleap path) DOES
   measure SG-SG distance because WE are the upstream authority on
   that path. The methodology is the AmberTools tutorial canonical;
   recorded in spec/plan/amber-implementation-plan-2026-04-29.md.

2. We do NOT parse GROMACS .r2b / .hdb / xlateat.dat at runtime, nor
   do we generate a project data file from them. The walk against
   topol.top rtp comments + structural match against canonical
   reference PDB obviates both. Fewer moving parts; less coupling
   to a specific GROMACS install version.

3. We do NOT add an in-code translation table beyond the small typed
   alias map already in NamingRegistry::to_canonical_. The typed
   discriminator hint (HISH→HIS during ingest) is fine; the actual
   canonicalisation happens by walk + match.

4. We do NOT use InterMol or ParmEd. Neither provides FF-agnostic
   naming canonicalisation; both preserve format-specific names.
   (Verified 2026-04-30 against InterMol intermol/atom.py and
   ParmEd documentation.)

5. We do NOT validate against PDB CCD as the primary authority.
   CCD is GLOBAL canonical; the local canonical (the project's prep
   dir) is more authoritative for the specific protein we're
   loading because IT is what was simulated. CCD becomes a periodic
   spot-check (run gemmi validation against the canonical reference
   PDB during fleet curation; not on every trajectory load).

6. We do NOT trust the .name field in the TPR for canonical residue
   identity. .name is the FF-port label (HISH, HISD, etc.); it does
   not propagate the rtp information. topol.top's rtp comment is the
   field that carries canonical AMBER on disk.

7. We do NOT make chemistry decisions. We read, walk, audit, and
   apply. The upstream authority owns the chemistry; we own the
   record-keeping.
```

## Discipline rule

```text
Every claim in the canonicalisation log must point at an
upstream-authority file. No claim is generated from consume-side
chemistry inference.

Methods text reads from the canonicalisation log. Reviewers can
spot-check any residue's claim by opening the cited file.

Chemistry decisions are READ from typed upstream records (PRMTOP,
topol.top rtp comment, input PDB) — never re-derived on the
consume side.
```

## Verified findings (2026-04-30)

Three reads against the actual files for 1P9J_5801 and 1Z9B_6577
confirm the architecture is grounded in what's on disk:

```text
Reference PDB sources/structure_extracted_1P9J_chain_A_model_0.pdb:
    Residues:    plain HIS, CYS, ASP, GLU, etc. (deposit-canonical)
    Atoms:       N, CA, C, O, CB, HB2, HB3, ND1, NE2, HD1, HE2, etc.
                 (PDB CCD / IUPAC 1979/1998 canonical)
    H1/H2/H3 at N-terminus (canonical capped)
    No protonation suffixes; no FF-port labels.

topol.top "; residue N name rtp rtp q charge" lines:
    Format consistent across both proteins.
    Regex: ^; residue\s+(\d+)\s+(\S+)\s+rtp\s+(\S+)\s+q\s+(\S+)$

    1P9J: 3 HIS → all HISH/HIP (user pipeline answered "2" thrice).
          6 CYS → all CYS/CYX (3 disulfide bridges; residues
                  8/16/22/33/35/44).
    1Z9B: 1 HIS → HISH/HIP (residue 75).
          1 HIS → HISE/HIE  (residue 121).
          (Mixed protonation; rtp field captures both choices.)

CYS atom-set delta (residue 8 of 1P9J):
    Reference:  C, CA, CB, H, HA, HB2, HB3, N, O, SG
    Trajectory: C, CA, CB, H, HA, HB1, HB2, N, O, SG
    Differences:
        - HB2/HB3 → HB1/HB2 (amber14sb hdb β-methylene numbering;
          label change only)
        - HG present in reference, absent in trajectory (real
          chemistry change; CYX rtp from specbond)
    Both detectable via the walk; both citable.
```

These three reads are the empirical foundation. Any future
canonicalisation log can be hand-spot-checked against them.

## Implementation order

```text
Step P0.NC1   Define typed CanonicalizationRecord (or similar) on
              ProteinBuildContext, plus a canonicalization log writer.

Step P0.NC2   Implement the topol.top rtp comment parser. Pure C++;
              regex-based; ~30 lines. Extracts (residue index,
              ff_port_name, rtp, charge) per residue. Loaded by
              the trajectory loader alongside the TPR.

Step P0.NC3   Implement the structural-match against the canonical
              reference PDB. Read the reference via PdbFileReader
              (existing); match heavy atoms by element +
              bonded-neighbor; H atoms by typed parent + AminoAcidType
              template. ~80 lines. Used only on the trajectory path
              (other paths' canonical authority is the input itself).

Step P0.NC4   Wire the walk into the trajectory loader. After
              FullSystemReader::ReadTopology and BuildProtein,
              run the walk; emit the log; apply canonical names.

Step P0.NC5   Tests: 1P9J + 1Z9B walks produce correct logs; each
              disulfide pair detected; each protonation variant
              correctly typed; HG missing on disulfide CYS recorded
              with the cited specbond.dat reason; β-methylene
              renames recorded with the amber14sb hdb reason.

Step P0.NC6   Extend log format to other load paths:
              ORCA / mutant: trivial log (canonical authority = PRMTOP;
                  no renames; chemistry from PRMTOP CHARGE/RADII).
              PDB / protonated-PDB (flat-table path): trivial (canonical
                  authority = input PDB).
              AmberPreparedChargeSource: log records WHAT WE DID
                  (capping policy, disulfide detection, atom mapping)
                  in the same format.

Step P0.NC7   Water and ion canonicalisation.
              Water:
                - Detect water model from topol.top #include line and
                  per-water atom count.
                - Apply the 4-row aggregate translation table
                  (SOL/OW/HW1/HW2 → HOH/O/H1/H2 for TIP3P).
                - Record water charges (already read from TPR by
                  FullSystemReader) into the canonicalisation log.
                - Emit aggregate water entry in the log (single block
                  for all 23933+ waters; not per-molecule).
              Ions:
                - Confirm GROMACS-port name matches CCD; emit aggregate
                  per-species block in the log; rename only if needed
                  (currently never for our monovalent ions).
              Virtual sites: out of scope; flagged as future when
                TIP4P+ trajectories enter scope.
              Calculators that match water atoms by FF-port name
              ("OW") instead of canonical ("O") get reviewed under
              Phase 2's per-calculator migration; not gated on this
              step.
```

NC1–NC5 land per-step like the AMBER charge slice. NC6 extends
the log surface to non-trajectory load paths; it doesn't change the
typed object model. NC7 adds non-protein chemistry (water, ions,
future virtual sites) to the same log surface.

## Sources

External (FF-agnostic vocabulary):

- [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)
- Westbrook, J. D. et al. (2015). The chemical component dictionary:
  complete descriptions of constituent molecules in experimentally
  determined 3D macromolecules in the Protein Data Bank.
  *Bioinformatics* 31(8):1274.
- [OpenFF Toolkit Topology — protein assignment via CCD chemistry,
  not residue labels](https://docs.openforcefield.org/projects/toolkit/en/latest/api/generated/openff.toolkit.topology.Topology.html)
- IUPAC-IUB Commission on Macromolecular Nomenclature (1979/1998).
  Atom names for amino acids; project ingestion via gemmi planned
  in spec/plan/amber-implementation-plan-2026-04-29.md N1.B.

External (GROMACS pdb2gmx behaviour, cited but not parsed at runtime):

- [GROMACS pdb2gmx documentation](https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html)
- [pdb2gmx input files reference](https://manual.gromacs.org/current/reference-manual/topologies/pdb2gmx-input-files.html)
- Filesystem on this machine:
    `/home/jessica/gromacs/src/gromacs-2026.0/share/top/amber14sb.ff/aminoacids.r2b`
    `/home/jessica/gromacs/src/gromacs-2026.0/share/top/amber14sb.ff/aminoacids.hdb`
    `/home/jessica/gromacs/src/gromacs-2026.0/share/top/xlateat.dat`
    `/home/jessica/gromacs/src/gromacs-2026.0/share/top/specbond.dat`
    `/home/jessica/gromacs/src/gromacs-2026.0/src/gromacs/gmxpreprocess/{pdb2gmx,pdb2top,hizzie,specbond,xlate}.cpp`

External (negative findings — what does NOT solve naming
canonicalisation, in case future readers are tempted):

- [InterMol intermol/atom.py — single .name field, no
  canonicalisation](https://github.com/shirtsgroup/InterMol/blob/master/intermol/atom.py)
- [ParmEd documentation — preserves format-specific names; does
  not translate](https://parmed.github.io/ParmEd/html/index.html)

In-tree:

- `spec/plan/amber-implementation-plan-2026-04-29.md` — AMBER
  charge slice (companion). The crystal projection rule, the
  substrate-vs-conformation split, the AMBER-as-standard lock.
- `spec/plan/current-topology-anchor-2026-04-29.md` — anchor for
  topology + charge work.
- `src/NamingRegistry.cpp` — existing typed alias map; extends only
  with the small typed-discriminator hints already required.
- `tests/data/fleet_amber/{1P9J_5801,1Z9B_6577}/` — empirical
  foundation; verified 2026-04-30.

## Doctrine

```text
Naming canonicalisation is structural alignment against an
upstream-authority record, not table lookup against a documentation
standard.

The upstream authority is local to each protein's prep tree.
The canonical AMBER record on disk for trajectories is the
topol.top rtp comment line.
The canonical name reference for atoms is the deposit-form PDB in
sources/.

The consume side reads, walks, audits, applies. It does not decide.
Chemistry decisions live with their upstream authorities.

The canonicalisation log is the methodology trail. Methods text
reads from it; reviewers spot-check via cited files.
```
