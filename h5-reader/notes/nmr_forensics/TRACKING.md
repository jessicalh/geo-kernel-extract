# NMR forensics — per-protein tracking

Explicit per-protein status for the 10 calibration proteins. Update
this file after each session's work.

**Methodology:** `PRINCIPLES.md` 7-step checklist (§125–144). Each
step must complete before the next; no skipping.

**Status vocabulary:**
- `not started` — no per-protein tables exist.
- `in progress <step>` — tables exist; reviewer is mid-checklist.
- `blocked <reason>` — requires external input (BMRB metadata, user
  adjudication, etc.).
- `reviewed` — all 7 steps complete; `REVIEW REQUIRED` header removed;
  findings appended to `SUMMARY.md`.

## Tracking table

| Protein | BMRB id | Status | Step 1 metadata | Step 2 tables | Step 3 audit | Step 4 observations | Step 5 re-run | Step 6 SUMMARY | Step 7 clean | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| 1B1V | 4292 | **reviewed** | ✓ PSP1 cytokine PDB; tag_alignment | ✓ offset=0; starter sufficed | ✓ | ✓ silent N-term; 43% stereo-ambig | ✓ | ✓ | ✓ | Pure-H deposition; 0 RefDB corrections; 1 disulfide; clean |
| 1BWX | 3449 | **reviewed** | ✓ Human PTH PDB; tag_alignment | ✓ offset=0; HSP alias works | ✓ | ✓ 3 HSP all untestable; CHARMM 1-39 vs BMRB 1-34 (5-res silent gap) | ✓ | ✓ | ✓ | Most HSP-heavy (3); no CYS, no PRO; flexible hormone |
| 1CBH | 192 | **reviewed** | ✓ Cellobiohydrolase I CBM PDB; tag_alignment | ✓ offset=0; HSP alias works | ✓ | ✓ HSP-4 untestable; ambig codes deposited | ✓ | ✓ | ✓ | Oldest entry (1989, BMRB 192); 0 UNMATCHED; THR-1 silent N-term |
| 1DV0 | 4757 | **reviewed** | ✓ reference.pdb TITLE + seq match | ✓ offset applied | ✓ | ✓ | ✓ 1 doc-and-carry | ✓ | ✓ | +2 offset resolved; 1 N-term residual; 125 RefDB corrections |
| 1G26 | 4656 | **reviewed** | ✓ Granulin A PDB; tag_alignment | ✓ offset=0; HSP alias works | ✓ | ✓ HSP untestable from NMR; 0 ambig codes deposited | ✓ | ✓ | ✓ | First HSP exemplar; 0 UNMATCHED; silent N-term VAL-1 |
| 1HA9 | 5028 | **reviewed** | ✓ MCoTI-II macrocyclic knottin PDB; tag_alignment | ✓ offset=0; cyclic-derived topology noted | ✓ | ✓ NEW N-term strategy: cyclic-derived non-zwitterionic | ✓ | ✓ | ✓ | 0 UNMATCHED (cyclic origin); 0 corrections; H-only |
| 1HD6 | 4820 | **reviewed** | ✓ Pheromone Er-22 PDB; tag_alignment | ✓ offset=0; starter sufficed | ✓ | ✓ silent N-term omission strategy noted | ✓ | ✓ | ✓ | Disulfide-rich (6 CYS); H+Cα coverage; 0 UNMATCHED; 54 STEREO_AMBIGUOUS (all ambig=2); BMRB-CSV-corrected verified |
| 1HS5 | 4934 | **reviewed** | ✓ PDB COMPND p53 frag 324-357; tag_alignment | ✓ offset=0; starter sufficed | ✓ | ✓ | ✓ no TOML edits drove adjudication | ✓ | ✓ | ASP-1 NH3+ confirmed; 96 RefDB corrections; ARG-19 ambig=4 noted; cross-cutting BMRB-CSV-is-corrected finding |
| 1I2V | 4976 | **reviewed** | ✓ Heliomicin defensin PDB; tag_alignment | ✓ offset=0; HSP alias works | ✓ | ✓ HSP-31 untestable; 0 ambig codes (1G26 pattern) | ✓ | ✓ | ✓ | 0 UNMATCHED; ASP-1 silent N-term; 6 CYS |
| 1I8X | 4351 | **reviewed** | ✓ Carp Granulin-1 PDB; tag_alignment | ✓ offset=0; HSP alias works | ✓ | ✓ HSP-3 untestable; 53 corrections, 53/53 match corrected | ✓ | ✓ | ✓ | H+C deposition; granulin pair to 1G26 |

## Next up

**All 10 proteins reviewed.** No remaining audits in scope.

If new proteins are added to the calibration set (beyond the current 10),
follow the PRINCIPLES.md 7-step checklist: read deposition metadata,
derive `<pdb>_<bmrb>.bmrb.toml` + `.refdb.toml` from
`_starter_translation.toml`, run `audit_nmr.py`, write AI observations,
classify resolutions, update SUMMARY.md, mark this row reviewed.

## Audit reruns

When the starter translation table (`_starter_translation.toml`) or
the audit tool (`audit_nmr.py`) changes, every reviewed protein's
audit must be re-run and outputs re-inspected for new findings. Log
reruns here with date + trigger.

_No reruns of 1DV0 needed yet — its audit still reflects the current
tool + starter state as of 2026-04-19. 1HS5 was reviewed against the
same tool + starter on 2026-04-19; its provenance edits (post-audit
findings_summary) are narrative-only and the audit was re-run after
the edits to confirm clean parse._
