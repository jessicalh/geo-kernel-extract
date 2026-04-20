# NMR experimental-shift forensics for the 10 calibration proteins

Per-protein, per-database audit of BMRB (raw) and RefDB (corrected)
experimental chemical shifts against the CHARMM topology used for MD
+ analysis trajectory extraction.

**Read first:** `PRINCIPLES.md`. The seven principles governing this
work are the spec; everything else is implementation.

## Layout

```
h5-reader/notes/nmr_forensics/
├── README.md                    (this file — orientation)
├── PRINCIPLES.md                (discipline for upstream-data forensics)
├── SUMMARY.md                   (state + per-protein findings)
├── _starter_translation.toml    (modern-IUPAC reference; NOT operative)
├── audit_nmr.py                 (audit tool)
├── 1DV0_4757.bmrb.toml          (reviewed exemplar)
├── 1DV0_4757.refdb.toml         (reviewed exemplar)
├── 1DV0_4757.audit.json         (tool output, machine)
└── 1DV0_4757.audit.md           (tool output, human)
```

The other 9 calibration proteins have no per-protein tables yet; the
tool refuses to run on them until each has its own reviewed tables
derived per the 7-step checklist in `PRINCIPLES.md`.

## Data source

`/shared/2026Thesis/thesis-md-data-access-restricted/`:
- `chemical_shifts_685.csv` — raw BMRB depositions (713k rows, 8
  cols: bmrb_id, residue_seq, residue_name, atom_name, atom_type,
  shift_value, shift_error, ambiguity_code).
- `refdb_corrected_shifts_685.csv` — RefDB re-referenced subset
  (374k rows, 7 cols; includes `original_value` and
  `corrected_value`).

## Running the audit

```bash
python3 audit_nmr.py /shared/2026Thesis/fleet_calibration-working/<pdb>_<bmrb> <bmrb_id>
```

Refuses if `<pdb>_<bmrb>.bmrb.toml` or `<pdb>_<bmrb>.refdb.toml` are
absent. Derive from `_starter_translation.toml` (copy + add
provenance header) before running.

Outputs:
- `<pdb>_<bmrb>.audit.json` — complete machine-readable record.
- `<pdb>_<bmrb>.audit.md` — per-db summaries + every row + CHARMM
  atom match status.

## Tool features

- **Per-protein, per-database translation tables.** No cross-entry
  assumption.
- **`[residue_seq_offset]` support.** Per-database offset applied
  before topology lookup. 1DV0's -2 offset is the worked example.
- **Huxley discipline.** Every anomaly complains loudly; no silent
  translations, aliases, or drops.
- **Full classification vocabulary:** `MATCHED`,
  `MATCHED_WITH_ANOMALIES`, `ONE_TO_MANY` (pseudo-atoms),
  `STEREO_AMBIGUOUS` (ambig code ≥ 2), `UNMATCHED_ROW`.
- **Full anomaly vocabulary:** `RESIDUE_SEQ_OUT_OF_RANGE`,
  `RESIDUE_NAME_MISMATCH`, `RESIDUE_ALIAS_APPLIED`
  (HSP/HSD/HSE → HIS etc.), `NAME_TRANSLATION_APPLIED`,
  `UNMATCHED_NAME`, `ONE_TO_MANY`, `STEREO_AMBIGUOUS`,
  `NUCLEUS_MISMATCH`, `REFDB_CORRECTED_VALUE_DIFFERS`, and others.

## Per-protein table schema

```toml
# Provenance block — what, why, when
# Protein: <pdb>_<bmrb>
# Database: bmrb | refdb
# Derived: YYYY-MM-DD from _starter_translation.toml
# Reviewer: <name>
# Deposition metadata consulted: <yes/no, sources>
# Findings:
#   - ...

[residue_seq_offset]
# Optional; omit if no offset. Applied to every row's residue_seq
# before topology lookup.
offset = -2

[backbone]
# Universal CHARMM → database atom-name map (applied to all residues)
HN = "H"

[residue.GLY]
# Per-residue map
HA1 = "HA2"
HA2 = "HA3"

[residue.ILE]
# CHARMM 'CD' → BMRB 'CD1' (naming divergence verified 2026-04-19
# during 1HS5 first pass)
CD = "CD1"
# ...

[residue_name_aliases]
# CHARMM protonation variants → database canonical name
HSP = "HIS"
HSD = "HIS"
HSE = "HIS"
```

## Policies explicitly deferred

The audit does NOT decide:

- Pseudo-atom expansion (BMRB methyl `HG2` → three CHARMM Hs):
  expand? average? distribute?
- Stereo-ambiguity (ambig ≥ 2 pairs): swap? average? dual-carry?
- Unmatched atoms: NaN? mean-impute? drop?
- RefDB vs BMRB when they diverge: prefer corrected, or carry both?
- HIS protonation resolution when BMRB is silent.
- N-terminal NH3+ vs single-H depositions (1DV0, 1HS5 pattern).

Each is a downstream analysis question the API must preserve
inputs for, not pre-answer.

## Getting started from cold — the mechanics

A successor session. You've read `PRINCIPLES.md`. You want to audit
the next protein. This is the concrete "how":

### Where the data is

- **CHARMM topology (trusted):** `reference.pdb` in each protein dir
  at `/shared/2026Thesis/fleet_calibration-working/<pdb>_<bmrb>/`.
- **Raw BMRB shifts:** one combined CSV at
  `/shared/2026Thesis/thesis-md-data-access-restricted/chemical_shifts_685.csv`.
  Filter rows by `bmrb_id` column.
- **RefDB re-referenced shifts:** one combined CSV at
  `/shared/2026Thesis/thesis-md-data-access-restricted/refdb_corrected_shifts_685.csv`.
  Filter rows by `bmrb_id`.
- **BMRB deposition metadata (external, required for provenance):**
  `https://bmrb.io/data_library/summary/?bmrbId=<id>` —
  deposition date, spectrometer conditions, pH, temperature,
  referencing compound, depositor, construct notes.

### Concrete first 30 minutes on a new protein

Example: taking on 1HS5/4934.

1. **Read BMRB metadata.** Open `bmrb.io/data_library/summary/?bmrbId=4934`.
   Note: deposition date, pH, temperature, construct description,
   referencing standard, any notes on N-terminal / C-terminal
   chemistry, protonation state of histidines if mentioned.
2. **Derive the two per-protein tables.** Copy `_starter_translation.toml`
   → `1HS5_4934.bmrb.toml` and `1HS5_4934.refdb.toml`. Edit each to:
   - Replace the starter's narrative-scope comment with a short
     narrative block specific to this protein.
   - Add a `[provenance]` structured section with all six required
     keys (see `1DV0_4757.bmrb.toml` exemplar). Schema validator in
     `audit_nmr.py` refuses to run without this.
3. **Run the audit.** `python3 audit_nmr.py /shared/2026Thesis/fleet_calibration-working/1HS5_4934 4934`.
   Produces `1HS5_4934.audit.json` and `.md` alongside the tables.
4. **Read the audit output line-by-line** (not just the summary).
   Every `UNMATCHED_ROW`, every `RESIDUE_NAME_MISMATCH`, every
   `REFDB_CORRECTED_VALUE_DIFFERS` row is a finding to annotate.
5. **Apply PRINCIPLES §5 + §6:** write an AI observation per anomaly
   class; classify each as adjudication / escalation / deferred
   downstream / document-and-carry. If adjudication, update the
   TOML (offset, translation rule, alias) and re-run.
6. **Append findings to `SUMMARY.md`** under a new `## 1HS5/4934`
   section. Cross-reference the TOML provenance block and the
   audit JSON/MD.
7. **Update `TRACKING.md`** — mark the protein's status as
   `reviewed`.

### What to expect from the audit output

- **Status counts** at the top — quick glance at the distribution
  of `MATCHED` / `MATCHED_WITH_ANOMALIES` / `ONE_TO_MANY` /
  `STEREO_AMBIGUOUS` / `UNMATCHED_ROW` per database.
- **Anomaly counts** — every anomaly classification tallied.
- **CHARMM coverage by element** — how many atoms got matched per
  database.
- **Every row tabulated** — per BMRB and per RefDB row, status +
  matched CHARMM atom + anomaly list.
- **Every CHARMM atom tabulated** — matched by BMRB, matched by
  RefDB, or unmatched (absence ≠ anomaly per se; note it).

## Downstream consumption — how forensics outputs connect to the Python API

*This section preserves the design intent; the Python CLI spec
(`PYTHON_CLI_SPEC.md`, pending) will operationalise it.*

Each audited protein produces:

- `<pdb>_<bmrb>.bmrb.toml` + `.refdb.toml` — translation tables
  + provenance + adjudications (offset, aliases, etc.).
- `<pdb>_<bmrb>.audit.json` — every BMRB row + every RefDB row +
  every CHARMM atom, each with status and anomaly metadata.
- `<pdb>_<bmrb>.audit.md` — human-reviewable counterpart.

The eventual Python extension module (`h5reader`) will consume
these so that every CHARMM atom in the typed object graph carries
an optional `experimental_shift` record:

```
atom.experimental_shift = ExperimentalShift(
    bmrb_raw = 8.29,                       # from CSV, unchanged
    refdb_corrected = 8.31,                # from RefDB CSV
    refdb_delta = +0.02,                   # corrected - raw
    audit_status = "MATCHED_WITH_ANOMALIES",
    anomalies = ["NAME_TRANSLATION_APPLIED: ..."],
    provenance = "1DV0_4757.audit.json#/bmrb_audit/42",
)
```

The Python layer does not decide pseudo-atom expansion, RefDB
preference, or HIS protonation — it presents both values with
audit context; downstream analysis decides per
`PRINCIPLES.md` §6.

**An atom with no experimental shift** is just `atom.experimental_shift
is None`. Absence is a legitimate state; no synthesis, no imputation
at the API layer.

**An atom with `document-and-carry` anomalies** (e.g. 1DV0 residue 1
N-terminal GLN H) has `experimental_shift.audit_status = "UNMATCHED_ROW"`
and its anomalies attached. The API exposes it; the analyst decides
whether to use it.

The audit JSON's per-row indexing gives stable provenance anchors.
Downstream queries like "what's this atom's RefDB-corrected shift?"
return the value with a provenance pointer the reviewer can trace
to the source row.

## Next work (2026-04-19)

Per-protein tables for the 9 remaining proteins. Status tracked in
`TRACKING.md`.

First-pass hypotheses from the pre-refactor audit (see `SUMMARY.md`)
are starting observations — verify each during review; do not adopt
as findings without running the audit under the new per-protein
table discipline.
