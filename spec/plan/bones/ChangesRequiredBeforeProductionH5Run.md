# Changes Required Before Production H5 Run

This document is the **single coordination point** for changes to the
`nmr-extract` pipeline that MUST be applied before running production
H5 extraction over the 685-protein fleet (scheduled ~2026-04-27).

## Scope

Entries here come from two sources, both active through the calibration
phase:

1. **Discoveries from the 10-protein calibration work itself** — e.g.,
   the 2026-04-20 BMRB/RefDB experimental-shift forensics that surfaced
   the ILE naming gap below.
2. **Discoveries that emerge while running statistical analyses on the
   10-protein calibration data** — upcoming, and *expected*. The
   mutations work covered most of this territory and a bit beyond, so
   we are close-to-right, but there is always something. The chances
   this ILE gap is the only thing we find are nil. Every subsequent
   find must be carefully folded in here before the fleet extraction,
   not patched in-place in application code.

Both sources feed the same list. This file is the living record.

## Discipline

- Naming changes go to the registry (`src/NamingRegistry.cpp`), never
  as ad-hoc translation rules in application code.
- Per-protein mapping files are acceptable *only* when per-protein
  provenance truly varies. If a change is uniform across the fleet
  (as the ILE gap below is), it belongs in the library as a canonical
  rule, not as per-protein data.
- No change listed here is activated incrementally. The fleet-wide
  vetting produces one coordinated library update, applied before
  production extraction begins. Activating pieces as they are found
  risks a partial fix that silently hides a different mismatch behind
  the accepted correction.
- Every entry below must record: discovery source, proposed code change
  (with pointer into `src/`), activation criteria, and — when applied —
  the activation date and commit SHA.

---

## 2026-04-20 — NamingRegistry: multiple CHARMM↔IUPAC coverage gaps plus one wildcard-β-methylene bug

**Status:** DEFERRED. Activate after fleet-wide BMRB/RefDB vetting.

**Source of discovery.** 10-protein experimental-shift forensics at
`h5-reader/notes/nmr_forensics/`. Empirical probe via
`_probe_naming_conventions.py` established BMRB and RefDB as uniform
IUPAC conventions across the 10 (expanded for BMRB, pseudo-collapsed
for RefDB). The packager precheck (`pack_experimental_shifts.py`
check 5) then tried to resolve every BMRB row to a CHARMM atom_index
via the H5 `/atoms/atom_name` dataset and enumerated every miss on
1DV0 — 8 categories of naming divergence surfaced. See `SUMMARY.md`
in the forensics directory for full tallies.

**Finding.** The H5 `/atoms/atom_name` dataset is built by
`AnalysisWriter` calling `NamingRegistry::TranslateAtomName` on each
CHARMM atom name. The registry's current rules cover only three cases:
backbone H↔HN, β-methylene HB1/HB2↔HB2/HB3, γ-methylene HG1/HG2↔HG2/HG3.
These don't cover the full CHARMM↔IUPAC divergences, and one of them
fires where it shouldn't.

### Missing rules — the H5 carries CHARMM names for these atoms

| Position | CHARMM name | IUPAC name | Residues affected |
|---|---|---|---|
| ILE γ-carbon | `CD` | `CD1` | ILE |
| ILE δ-methyl | `HD1`/`HD2`/`HD3` | `HD11`/`HD12`/`HD13` | ILE |
| ILE γ1-methylene | `HG11`/`HG12` | `HG12`/`HG13` | ILE |
| δ-methylene | `HD1`/`HD2` | `HD2`/`HD3` | ARG, LYS, PRO |
| ε-methylene | `HE1`/`HE2` | `HE2`/`HE3` | LYS |
| α-methylene | `HA1`/`HA2` | `HA2`/`HA3` | GLY |

CHARMM uses 1-based numbering consistently on methylenes; IUPAC uses
2/3-based. The existing β-methylene and γ-methylene rules handle
those two positions (with the ALA exception noted below); δ, ε, α,
and the ILE-specific locant + three-deep methyl numbering are
uncovered.

### Library bug — the H5 has duplicate atom_names on ALA β-methyl

CHARMM ALA β-methyl has three H atoms named HB1, HB2, HB3. The
registry's wildcard β-methylene rule (residue = `"*"`) applies to all
residues, but ALA's β position carries a 3-H methyl, not a 2-H
methylene. The translation under the current wildcard produces:

    CHARMM HB1  →  IUPAC HB2      (rule HB1→HB2 fires)
    CHARMM HB2  →  IUPAC HB3      (rule HB2→HB3 fires)
    CHARMM HB3  →  IUPAC HB3      (no rule, passes through)

The H5 `/atoms/atom_name` dataset therefore carries **two atoms named
HB3** and one atom named HB2 at every ALA residue. The underlying
atom data (positions, charges, shielding tensors, bond topology) is
correct — three distinct atom indices with parent = CB. Only the
`atom_name` display label is corrupted.

The library-side fix: narrow the wildcard β-methylene rule to exclude
ALA (or more defensibly: make it residue-set-specific, listing every
residue that has a β-methylene — ASP, ASN, CYS, GLU, GLN, HIS, LEU,
PHE, PRO, SER, TRP, TYR, LYS, ARG, MET). ALA β-methyl uses
HB1/HB2/HB3 identically in CHARMM and IUPAC; no translation is needed.

### Consequence (pre-existing H5 files on the 10-protein set)

All eight categories above surface as CHARMM (or corrupted-IUPAC for
ALA) names in the already-generated H5 files under
`fleet_calibration-{working,stats,backup}`. Downstream consumers that
assume uniform IUPAC — e.g., BMRB/RefDB shift binding against
`/atoms/atom_name` — fail across all these positions. Regenerating
the existing H5 files is not on the table.

For the 10-protein experimental-shifts work (2026-04-20), every case
is handled in-place by two typed Python data constants in
`h5-reader/notes/nmr_forensics/pack_experimental_shifts.py`:

- `H5_NAME_FROM_IUPAC` — keyed on (residue_type, IUPAC atom name),
  returns H5 atom name. Covers every missing-rule case in the table
  above. Direct dict lookup, no string parsing.
- `BMRB_METHYL_VIA_PARENT` — keyed on (residue_type, BMRB atom name)
  for ALA specifically. Binds the three BMRB rows (HB1/HB2/HB3) via
  `/topology/parent_atom_index` — the three CHARMM atoms that have
  CB as parent — and flags the binding `AMBIGUOUS_METHYL_ORDER`.
  Methyl Hs are chemically indistinguishable, so no physical
  information is lost; the ambiguity is in which BMRB row-number
  corresponds to which atom_index, which was never meaningful.

### Proposed code change

In `NamingRegistry::InitialiseAtomNameRules()`:

- **Narrow the wildcard β-methylene rule** to exclude ALA (or become
  residue-set-specific). Currently the `"*"` wildcard over-fires on
  ALA's 3-H methyl.
- **Add rules for missing cases**: bidirectional CHARMM↔Standard pairs
  for each divergence — ILE δ-methyl (6 rules), ILE γ-carbon CD↔CD1
  (2 rules), ILE γ1-methylene (4 rules), δ-methylene for ARG/LYS/PRO
  (4 rules × 3 residues = 12 but can use wildcard if all three
  behave the same — verify), LYS ε-methylene (4 rules), GLY
  α-methylene (4 rules).

Exact lines are present in `src/NamingRegistry.cpp` today as
commented-out code with the full rationale inline. Search for
`// ILE-specific CHARMM ↔ Standard rules` and the surrounding block.

### Activation criteria

1. Fleet of 685 proteins lands.
2. Re-run `_probe_naming_conventions.py` over all 685 × 2 = 1370
   deposition records. Confirm convention uniformity fleet-wide.
3. Run `pack_experimental_shifts.py` precheck 5 on one protein from
   the fleet to enumerate any remaining H5-vs-IUPAC gaps.
4. **If these eight categories are the only gaps:** uncomment the
   added rules AND apply the ALA β-methylene narrowing in
   `src/NamingRegistry.cpp`; run unit tests.
5. **If additional gaps surface:** add them as new sub-sections under
   this 2026-04-20 entry (this discovery's section accumulates all
   findings from the same audit pass), then activate together as a
   single coordinated change. No partial activation.
6. Re-run calibration-protein H5 regeneration on one protein first;
   verify `/atoms/atom_name` now carries IUPAC for every one of the
   eight positions — spot-check ALA β-methyl has three distinct names
   (HB1/HB2/HB3), ILE has CD1/HD11/HD12/HD13, LYS has HD2/HD3 and
   HE2/HE3, GLY has HA2/HA3.
7. Begin fleet-wide production H5 extraction.
8. Record activation date and commit SHA in the section header
   (flip status from DEFERRED to APPLIED).

### Consequences for already-generated calibration H5s

None. Regeneration is not in scope. The 10-protein local constants
in `pack_experimental_shifts.py` stay in place and continue to bind
shifts against the CHARMM-named atoms in the existing H5s. The
module carries a comment documenting that the constants are
short-term and become unnecessary when the library fix activates
AND the 10 proteins are re-extracted (neither currently planned).
