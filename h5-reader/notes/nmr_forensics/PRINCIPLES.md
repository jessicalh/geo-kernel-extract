# Principles for upstream-data forensics

Codified 2026-04-19 during the BMRB/RefDB audit of the 10 calibration
proteins. Generalises to any adjudication between our clean internal
pipeline (CHARMM/ff14SB → MD → analysis H5) and upstream-data sources
of uncertain provenance (BMRB, RefDB, PDB headers, published tables,
consolidated training corpora, third-party chemical-shift
depositions).

## Why these principles exist

Upstream data is not like our internal data. Our pipeline is one
rigorous path, authored by one person, with one set of conventions.
BMRB is 40 years of community depositions from thousands of labs,
multiple IUPAC revisions, software-specific variants, inconsistent
stereo-assignment, ambiguous protonation states, construct-numbering
offsets, silent re-referencing. A user who curated 2400 candidate
proteins down to 685 survivors reports that even the survivors have
per-entry anomalies.

A successor coming to this work without these principles will
pattern-match — applying a translation "that usually works" — and
silently miss real findings. The principles are the discipline that
keeps findings loud.

## The seven principles

### 1. Be empirical

Record what the upstream data says, per row, per atom, per residue.
Do not infer. Do not pattern-match. Do not backfill from "what it
should be." If an entry claims atom X exists at residue Y, record
that claim; do not evaluate it against your expectation until the
adjudication step, which is a separate, explicit act.

### 2. Apply no assumptions

No assumption that entries share conventions, even within the same
database. No assumption that a modern IUPAC translation applies to
pre-1998 depositions. No assumption that RefDB corrections are
monotonic with BMRB originals. No assumption that identical atom
names across entries mean identical physics. No assumption that
residue numbers align across depositors.

When in doubt, the default answer is "the audit preserves the data;
the human adjudicates downstream." The audit never decides.

### 3. Every protein is its own case

Ten proteins means ten distinct forensics passes. The pattern that
worked for protein A may not work for protein B even if they share
a residue sequence. 1DV0's +2 offset came from a construct
convention; other proteins have other constructs and no offset.
Whatever is true about how BMRB entry 192 represents 1CBH is not
automatically true about how BMRB entry 4757 represents 1DV0.

There is no shortcut. No "apply the same table to all." Every
protein requires its own deep read.

### 4. Records are orthogonal across databases

BMRB raw and RefDB corrected are independent depositions of the
same underlying experiment. BMRB does not automatically provide
what RefDB expects; RefDB does not automatically inherit BMRB's
atom-naming scheme. One protein has TWO records — one per
database — each with its own per-protein translation table.

Treating them as interchangeable silently hides re-referencing
signals (e.g. 1DV0's 125 RefDB corrections, 1HS5's 96) and
naming-convention divergences (e.g. RefDB's atom-name
normalisation when BMRB's deposition used pseudo-atoms).

### 5. Intelligent AI observations on every anomaly

Every anomaly the audit flags carries both a mechanical
classification (e.g. `UNMATCHED_NAME`, `RESIDUE_NAME_MISMATCH`,
`REFDB_CORRECTED_VALUE_DIFFERS`) and — where a reasoning agent
examined it — a written observation attaching the finding to what
the data is likely telling us.

Example: 1DV0's 140 UNMATCHED + 406 NAME_MISMATCH firings are
mechanically "unmatched rows." The AI observation is: "BMRB used
full-construct residue numbering; PDB was re-numbered 1-based; the
sequences align cleanly under a -2 offset. Origin: HHR23A UBA
domain deposition predates modern IUPAC guidance on PDB-matched
numbering."

The AI observation is not an adjudication. It is a hypothesis about
what the data is showing, offered for the human reviewer's decision.

### 6. Sane resolution

Every anomaly resolves to exactly one of:

1. **An adjudication** — a decision encoded in the per-protein
   translation table (e.g. `[residue_seq_offset]` for 1DV0). The
   TOML's provenance block records what was decided and why.
2. **An escalation** — requires information the audit cannot
   access (e.g. 1HS5's ASP-1 N-terminal modification requires BMRB
   deposition metadata). The provenance block records what to
   check and where.
3. **A deferred downstream policy** — not audit-scope (e.g.
   pseudo-atom expansion, stereo-swap, unmatched-handling). The
   audit preserves both values/alternatives; downstream analysis
   decides.
4. **A document-and-carry** — the anomaly persists because the
   data is genuinely heterogeneous (e.g. the BMRB N-terminal H
   pattern seen in both 1DV0 and 1HS5 will not resolve; it's how
   BMRB represents that physics). The audit keeps the anomaly in
   the output for every subsequent run.

No silent drops. No "fix" that hides the signal.

### 7. Recorded

Every decision, hypothesis, adjudication, deferral, and escalation
lives in a durable artifact: the per-protein TOML provenance
block, the audit JSON/MD, the `SUMMARY.md`, the `SESSION_<date>.md`.
A successor opening this directory cold can trace every choice to
a written reason.

No decision in conversation or in memory alone. Memory is volatile
and session-scoped; the repository is the durable record.

## Checklist for a new protein audit

When adding an N+1'th protein to this discipline:

1. Read the upstream entry's header/metadata (deposition date,
   software used, construct description). Before running the tool.
2. Derive `<pdb>_<bmrb>.bmrb.toml` and `<pdb>_<bmrb>.refdb.toml`
   from `_starter_translation.toml`. Immediately update the
   provenance block to record what you read in step 1.
3. Run the audit. Read the full output (not just the summary).
4. For every anomaly: apply principle 5 (write an AI observation)
   then principle 6 (classify the resolution). For any anomaly
   that becomes an adjudication, update the per-protein TOML.
5. Re-run the audit after any TOML change. Verify the adjudicated
   anomalies no longer fire (or fire with a documented
   document-and-carry status).
6. Append to `SUMMARY.md` with the findings + AI observations for
   the new protein.
7. Remove the `REVIEW REQUIRED` provenance from the per-protein
   TOML header only after steps 1-6 are complete.

## Generalisation

These principles apply to any reconciliation between a rigorous
internal representation and an uncertain external deposition:

- PDB header vs MD topology.
- Force-field parameter tables vs published references.
- Calibration DFT vs experimental shielding references.
- Consolidated training data vs per-protein provenance.
- ML training corpora assembled from multiple sources.

In every case: treat external data with utmost suspicion, keep
per-source records orthogonal, flag every inconsistency, surface
AI observations explicitly, resolve sanely, record everywhere.
