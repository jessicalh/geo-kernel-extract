# GromacsToAmberReadbackBlock — design

**Date:** 2026-05-02
**Status:** Design only. No implementation yet. Compiler-trace shape; not live state.

This document records the agreed shape of the import-time readback
block for AMBER-via-GROMACS trajectory loads. The block reads back
what GROMACS decided about chemistry (during pdb2gmx prep) and
applies those decisions as typed facts on the existing object model;
it does not re-decide and it does not stay resident.

It supersedes both the "NamingRegistry expansion for HISH" approach
and the prior session's "long-lived canonicalisation state on a typed
object" approach. Naming-canonicalisation is one strand of the work
the block does (HISH → HID where the rtp says so); the rest is
adopting GROMACS's chemistry decisions across atom-set, bonded list,
charges, terminal templates, and disulfide pairing.

Core constraints:

- It is a **compiler trace**, not live state. Calculators never see it.
  See memory `feedback_readback_block_is_a_compiler_trace`.
- The merge is **selective authority**, not "revert to original PDB"
  and not "alias HISH → HIS." Different facts have different
  authorities; the record records which authority was applied where.
- Compiled facts land on the existing typed object model. Nothing new
  stays resident.

## Inputs (authority files)

| Authority | Source | Provides |
|---|---|---|
| TPR | `production.tpr` (libgromacs `read_tpx_state`) | atom/residue order, charges, atom set per residue, bonded interactions, exclusions, FF-numerical invariants |
| topol.top rtp comment | `prep_run_*/topol.top` "; residue N <name> rtp <rtp> q <q>" lines | **canonical AMBER residue label** (HID/HIE/HIP/CYX/CYM/ASH/GLH/LYN/...) when GROMACS rewrote `.name` to an FF-port label (HISH/HISD/HISE) |
| Extracted PDB | `prep_run_*/sources/structure_extracted_*.pdb` (or `pdb_companion_*.pdb` if structure_extract did not run) | deposit-canonical IUPAC atom names where GROMACS renamed (HB1/HB2 ↔ HB2/HB3, ILE CD ↔ CD1, OXT ↔ OC1/OC2) |
| decisions.json | `prep_run_*/decisions.json` | prep-time chemistry decisions (titration, ion strategy, pH, capping policy) — provenance only |

## Selective-authority table (the merge logic)

| Fact | Authority | Notes |
|---|---|---|
| atom / residue order | TPR | what GROMACS wrote is what we walk |
| atom set membership (which atoms exist on a residue) | TPR | CYX has no HG; HIP has both HD1+HE2; ASH has HD2 |
| force-field partial charges | TPR | per-atom q from `tpr_atoms.atom[ai].q` |
| disulfide pairing (CYS → CYX) | TPR | SG-SG bond present in topology + atom set delta (HG removed) |
| explicit prep protonation (HID/HIE/HIP/CYX/CYM/ASH/GLH/LYN) | topol.top rtp comment line | the rtp field is the canonical AMBER name even when `.name` is the FF-port label |
| intermediate residue spelling (HISH / HISD / HISE / etc.) | topol.top rtp + atom-set | resolved against canonical via rtp; never aliased into NamingRegistry |
| atom-name spelling (HB1/HB2/HB3, ILE CD vs CD1, OC1/OC2 vs OXT) | extracted PDB | structural alignment by atom-set against deposit-canonical PDB |
| canonical typed identity | internal vocabulary (AminoAcidType + Element + AtomRole enums) | always the destination; the record explains how we got here |

## Outputs (compiled facts on existing typed objects)

The record applies these typed fields and then is discarded:

```text
Residue.protonation_variant_index  HID/HIE/HIP/CYX/CYM/ASH/GLH/LYN/...
                                   matched per AminoAcidType.h canonical
                                   indices (asserted by ValidateVariantIndices).
                                   Resolved by VariantIndexFromForceFieldName
                                   given the rtp-field name (NOT the .name
                                   FF-port label).

Residue.terminal_state             N/C/N+C/Internal (existing, set by
                                   Protein::ResolveResidueTerminalStates;
                                   could be over-ridden by the record if
                                   the rtp ever specifies cap residues
                                   like ACE/NME).

Atom::pdb_atom_name (canonical)    after merging extracted-PDB names back
                                   onto FF-port-renamed atoms (HB1/HB2/HB3
                                   ↔ HB2/HB3, etc.)

Protein atom set                   implicit; only present atoms exist.
                                   No marker needed — CYX has no HG
                                   ConformationAtom because pdb2gmx didn't
                                   write one.

CovalentTopology bonds             SG-SG inter-residue bonds tagged
                                   BondCategory::Disulfide from prep
                                   authority, not re-decided from
                                   distance heuristic.

LegacyAmberTopology FF data        mass / atom-type / ptype / atomtype-
                                   string / exclusions / fudge_qq / atnr
                                   / num_non_perturbed — already plumbed
                                   via LegacyAmberInvariants value-pack
                                   (commits landed 2026-05-02).
```

## Audit emission (the trace, materialized to disk)

Emit one JSON file alongside the topology:

```text
prep_run_*/.../gromacs_to_amber_readback_block.json
```

Contents (sketch):

```json
{
  "schema_version": 1,
  "load_path": "amber-trajectory",
  "authority_files": {
    "tpr": "production.tpr (sha256: ...)",
    "topol_top": "../topol.top (sha256: ...)",
    "extracted_pdb": "../sources/structure_extracted_<id>_chain_A_model_0.pdb (sha256: ...)",
    "decisions_json": "../decisions.json (sha256: ...)"
  },
  "residues": [
    {
      "index": 3,
      "tpr_name": "HISH",
      "rtp": "HIP",
      "canonical": "HIS",
      "variant_index": 2,
      "terminal_state": "internal",
      "source_line": "; residue   3 HISH rtp HIP q +1"
    },
    {
      "index": 7,
      "tpr_name": "CYS",
      "rtp": "CYX",
      "canonical": "CYS",
      "variant_index": 0,
      "atom_set_delta": ["HG: stripped (CYX has no sulfhydryl proton)"],
      "disulfide_partner_index": 41,
      "source_line": "; residue   7 CYS rtp CYX q  0"
    }
  ],
  "atom_renames": [
    {
      "atom_index": 23,
      "ff_port_name": "HB1",
      "canonical_name": "HB2",
      "match_method": "atom-set-aligned-against-extracted-pdb"
    }
  ],
  "warnings": [],
  "summary": {
    "n_residues": 56,
    "n_protonation_variants_applied": 4,
    "n_disulfide_bonds": 1,
    "n_atom_renames": 12
  }
}
```

The JSON is for debugging, methods text, and reviewer spot-check. NOT
for any calculator to consume at runtime. The file is regeneratable
from the authority files, so it does not need to be persisted across
runs (though emitting it at extraction time is cheap and convenient).

## What does NOT live on the record

These are explicitly out of scope (waters/ions are not in the typed
substrate per the 2026-04-30 walk-back):

- water moltype names / SETTLE geometry
- ion atom positions or names
- virtual-site geometry for water M-sites
- LJ pair table (FF-wide nonbonded numerics; not topology)

## Where the record lives in code (when implemented)

- A struct in the loader's stack frame (e.g. inside
  `FullSystemReader::BuildProtein` or a helper called from it).
- Parsers for topol.top rtp comments + extracted PDB live as
  free functions in the loader module.
- The merge logic is a single function: `BuildGromacsToAmberReadbackBlock(
  tpr_data, topol_path, extracted_pdb_path, decisions_path)
  -> GromacsToAmberReadbackBlock`.
- Application: `ApplyGromacsToAmberReadbackBlock(record, protein)` walks
  the typed Protein (Residues, Atoms, CovalentTopology) and writes
  the typed facts. Returns void; logs warnings.
- Emission: `EmitGromacsToAmberReadbackBlock(record, output_path)` writes
  the JSON.
- Disposal: the record goes out of scope at the end of `BuildProtein`.
- Calculators NEVER see the struct — there is no accessor for it on
  Protein, ProteinConformation, LegacyAmberTopology, or any other
  long-lived object.

## Why this is not "NamingRegistry expansion"

The NamingRegistry approach (add HISH/HISD/HISE aliases to the global
canonicalization table) was rejected for three reasons:

1. **It hides what GROMACS decided.** A residue labelled HISH could
   have rtp=HIP (charged) or rtp=HID/HIE if a tool drift renamed it.
   Aliasing flattens that signal.

2. **The decision is in the rtp, not in the name.** Reading rtp
   directly is the structurally correct path. The .name field is a
   transport artifact of pdb2gmx's amber14sb port.

3. **It pollutes the stable vocabulary** for all paths, including
   ones where HISH never appears (PDB load, ORCA load).

The record is the right shape because it merges from typed authorities
with explicit provenance, then disappears.

## Implementation order (when this is the active task)

1. Implement `topol.top` rtp comment parser (free function, returns a
   per-residue map).
2. Implement extracted-PDB reader (atom names by residue index, for
   the structural alignment pass).
3. Build `GromacsToAmberReadbackBlock` in `FullSystemReader::BuildProtein`
   from those inputs + the parsed TPR.
4. Apply the record to the typed Protein before
   `Protein::FinalizeConstruction` returns.
5. Emit the JSON audit file.
6. Wire `AmberTrajectoryFixtureTest` to verify HISH-class names
   resolve through the record (not aliased through NamingRegistry).
7. Verify the existing `BuildsAmberProtein_1P9J` / `_1Z9B` tests
   pass.

## Companion memory entries

- `feedback_readback_block_is_a_compiler_trace` —
  the compiler-trace rule.
- `feedback_capture_at_the_boundary` — the capture rule that
  populates `LegacyAmberInvariants`.
- `feedback_no_attach_lifecycle_for_invariant_data` — the rule
  that makes `LegacyAmberTopology` plain-fields rather than
  Optional/Attach.
- `project_charmm_retired_amber_only_2026-05-02` — the project
  state this design fits into.

## What this document is NOT

Not a plan-of-record for sequencing other work. Not an implementation
schedule. Not a user-facing brief. The record exists IF AND WHEN this
becomes the active task; until then, the test failures
(`BuildsAmberProtein_1P9J` / `_1Z9B`) are the gate that signals the
work is needed. They are pre-existing and not regressions from the
2026-05-02 cleanup pass.
