# NMR forensics — state and findings summary

Audit state as of 2026-04-19:

- **Tool:** `audit_nmr.py` — per-protein per-database audit, consumes
  `[residue_seq_offset]` to apply construct-numbering corrections.
  Refuses to run without per-protein tables (no silent fallback).
- **Starter:** `_starter_translation.toml` — reference; modern IUPAC
  conventions. Not operative.
- **Reviewed:** all 10 proteins (1DV0/4757, 1HS5/4934, 1HD6/4820,
  1B1V/4292, 1HA9/5028, 1G26/4656, 1CBH/192, 1I2V/4976, 1I8X/4351,
  1BWX/3449). Per-protein tables, audits, and provenance complete.
- **Principles codified in:** `PRINCIPLES.md`.

## Cross-cutting CORRECTION — BMRB CSV is post-RefDB-correction (2026-04-19)

Established during the 1HS5 audit and back-verified against 1DV0:
`chemical_shifts_685.csv`'s `shift_value` column is the **RefDB-corrected
value** wherever RefDB provides one, not the BMRB raw deposition value.

| BMRB id | RefDB rows where corrected ≠ original | BMRB CSV matches `corrected` | matches `original` | matches neither |
|---|---|---|---|---|
| 4757 (1DV0) | 125 | 125 | 0 | 0 |
| 4934 (1HS5) | 96 | 96 | 0 | 0 |
| 4820 (1HD6) | 37 | 37 | 0 | 0 |
| 4292 (1B1V) | 0 (no corrections) | — | — | — |
| 5028 (1HA9) | 0 (no corrections) | — | — | — |
| 4656 (1G26) | 0 (no corrections) | — | — | — |
| 192 (1CBH) | 0 (no corrections) | — | — | — |
| 4976 (1I2V) | 0 (no corrections) | — | — | — |
| 4351 (1I8X) | 53 | 53 | 0 | 0 |
| 3449 (1BWX) | 0 (no corrections) | — | — | — |

The 1DV0 SUMMARY's earlier wording "prefer RefDB `corrected_value` over
BMRB raw `shift_value`" is misleading. There is no BMRB raw to prefer —
the CSV already carries the corrected number. RefDB's standalone value is
its `original_value` column (the change-log of which rows were
re-referenced and by how much), not a separate truth source on
`shift_value`. To be re-verified for the remaining 8 proteins as each is
audited.

## 1DV0/4757 — the exemplar (reviewed)

**Protein:** DNA repair protein HHR23A, N-terminal UBA domain (45
residues).

**Sources of record (for traceability):**
- TOML provenance (BMRB): [`1DV0_4757.bmrb.toml`](1DV0_4757.bmrb.toml) — `[provenance]` + `[residue_seq_offset]` sections.
- TOML provenance (RefDB): [`1DV0_4757.refdb.toml`](1DV0_4757.refdb.toml) — `[provenance]` + `[residue_seq_offset]` sections.
- Audit outputs: [`1DV0_4757.audit.json`](1DV0_4757.audit.json), [`1DV0_4757.audit.md`](1DV0_4757.audit.md).

**Adjudication applied:** `[residue_seq_offset] offset = -2` on both
BMRB and RefDB tables. BMRB deposition used full-construct numbering;
our PDB was re-numbered 1-based after MD prep. Our residues 1..45 ≡
BMRB residues 3..47, identical sequence.

**After offset — clean state:**

| Database | Rows | Matched (incl. anomalies) | Unmatched | Coverage |
|---|---|---|---|---|
| BMRB | 434 | 433 | 1 | H: 308/349, C: 154/230, N: 41/58 |
| RefDB | 213 | 212 | 1 | H: 87/349, C: 84/230, N: 41/58 |

**AI observations:**

- The remaining 1 unmatched row (identical in BMRB and RefDB) is
  residue 3 (CHARMM seq 1) GLN backbone `H` at 8.86 ppm. Our CHARMM
  residue 1 has zwitterionic NH3+ (H1/H2/H3), no single backbone H.
  This is a **fundamental sample-vs-force-field chemistry mismatch**,
  not a translation gap. Same pattern anticipated at 1HS5's residue 1
  ASP. Policy: document-and-carry; downstream analysis decides how
  to reconcile (drop, map to pseudo-average of H1/H2/H3, or apply a
  pre-analysis N-terminal correction).
- **RefDB re-referencing is substantial.** 125 of 213 RefDB rows
  (59%) have `corrected_value ≠ original_value`. This is the highest
  re-referencing rate observed across the 10 (first-pass data).
  Downstream comparison to our DFT/kernel predictions must use
  RefDB's `corrected_value`, not BMRB's raw `shift_value`.
- **Heteronuclear coverage is strong.** 41/58 N (71%), 154/230 C
  (67% via BMRB, 36% via RefDB). Among the 10, 1DV0 is one of the
  two proteins (with 1HS5) carrying full H+C+N assignment depth.
- **No HSP residues.** No protonation adjudication required.
- **ONE_TO_MANY pseudos (24)** for methyl groups (THR HG2, VAL HG1/
  HG2, LEU HD1/HD2, etc.). Expansion policy deferred to downstream.
- **STEREO_AMBIGUOUS (4)** — low for a protein of this size;
  depositors stereo-assigned most β-methylenes.

## 1HS5/4934 — second exemplar (reviewed)

**Protein:** designed p53 tetramerisation-domain dimer fragment (residues
324-357 of human p53; PDB header: TITLE "NMR SOLUTION STRUCTURE OF
DESIGNED P53 DIMER", deposited 2000-12-22, ENGINEERED YES, MUTATION YES).
Chain A used for MD; CHARMM 34 residues, BMRB 33.

**Sources of record:**
- TOML provenance (BMRB): [`1HS5_4934.bmrb.toml`](1HS5_4934.bmrb.toml).
- TOML provenance (RefDB): [`1HS5_4934.refdb.toml`](1HS5_4934.refdb.toml).
- Audit outputs: [`1HS5_4934.audit.json`](1HS5_4934.audit.json), [`1HS5_4934.audit.md`](1HS5_4934.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0` (identity
numbering). Starter translation table sufficed.

**Audit outcome (after run):**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 336 | 313 | 12 ONE_TO_MANY, 10 STEREO_AMBIGUOUS, 1 UNMATCHED_ROW | C 95/178, H 227/278, N 37/55 |
| RefDB | 159 | 158 | 1 UNMATCHED_ROW | C 63/178, H 62/278, N 33/55 |

**AI observations:**

- **N-terminal NH3+ pattern confirmed.** ASP-1 BMRB row `(H, 8.29 ppm)`
  is the 1 UNMATCHED. CHARMM ASP-1 has H1/H2/H3 (zwitterionic NH3+) and
  no single backbone H. Identical mechanism to 1DV0/GLN-3. Document-and-
  carry; the audit keeps the row visible on every subsequent run.
- **No HSP residues.** No protonation adjudication required for 1HS5.
- **ARG-19 ambiguity_code = 4** on HB3 + HG3 — atypically high (versus
  the routine ambig=2 for sidechain amide pairs). Code 4 indicates
  intra-residue ambiguity (depositor unsure which proton is which); a
  stronger downstream-policy flag than ambig=2.
- **Methyl pseudos at 12 atoms** (THR-6 HG2; LEU-7/25/27 HD1+HD2; ILE-9
  HD1+HG2; ALA-24/30/32 HB) — ONE_TO_MANY, expansion deferred.
- **CHARMM C-term LYS-34 is BMRB-silent.** Topology has 34 residues; BMRB
  4934's coverage stops at GLY-33. Not an anomaly — there's nothing to
  match against and nothing fires.
- **N coverage strong (37/55, 67%).** Same H+C+N depth as 1DV0; together
  these are the two N-bearing proteins among the 10.
- **Re-referencing magnitude.** 96 corrected rows. Modest C shifts
  (Δ ≈ +0.46 ppm on backbone C at residue 1) and larger N shifts
  (Δ ≈ −2.23 ppm on backbone N at residues 1 and 2). Suggests an
  N-channel reference change applied to the deposition.

## 1HD6/4820 — third reviewed protein

**Protein:** Pheromone Er-22 — 37-residue disulfide-rich peptide from the
ciliate *Euplotes raikovi* (PDB 1HD6, deposited 2000-11-09; BMRB title
"1H and 13Ca Chemical Shift Assignments for the Pheromone Er-22"). Single
chain A, identity numbering, no engineered mutation, no expression tag.

**Sources of record:**
- TOML provenance (BMRB): [`1HD6_4820.bmrb.toml`](1HD6_4820.bmrb.toml).
- TOML provenance (RefDB): [`1HD6_4820.refdb.toml`](1HD6_4820.refdb.toml).
- Audit outputs: [`1HD6_4820.audit.json`](1HD6_4820.audit.json), [`1HD6_4820.audit.md`](1HD6_4820.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`. Starter
translation table sufficed.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 240 | 165 | 21 ONE_TO_MANY, 54 STEREO_AMBIGUOUS, 0 UNMATCHED | C 37/166 (Cα-only), H 245/256 (96%), N 0/42 |
| RefDB | 105 | 105 | 0 UNMATCHED | C 37/166, H 68/256, N 0/42 |

**AI observations:**

- **Zero UNMATCHED — silent N-terminal handling.** ASP-1 BMRB rows are
  CA, HA, HB2, HB3 only. The depositor simply omitted the backbone amide
  H rather than depositing a single-H value at an NH3+ N-terminus. Same
  underlying physics as 1DV0/1HS5 (zwitterionic NH3+ in solution), opposite
  depositor strategy (silent omission vs. confident single-H deposit). The
  audit catches the 1DV0/1HS5 case as UNMATCHED; the 1HD6 case as
  invisible. Worth noting downstream: BMRB N-terminal H presence is *not*
  evidence of single-H chemistry — depositor convention varies.
- **Disulfide-rich peptide, CHARMM CYS naming uniform.** 6 CYS (residues
  3, 10, 15, 18, 24, 32) likely form 3 disulfide bridges. CHARMM topology
  names all six as CYS (no CYX); no audit anomaly fires; no per-protein
  alias needed. (Re-examine if any HG (free thiol) rows fire UNMATCHED —
  none did.)
- **C coverage is Cα-only.** 37 C rows = 1 per residue (Cα). BMRB title
  is explicit: "1H and 13Ca". No Cβ/Cγ/etc. assignments deposited.
- **H coverage exceptional (245/256 = 96%).** Includes ONE_TO_MANY methyl
  pseudo expansions (1 BMRB pseudo → 3 CHARMM Hs). Of all 10 calibration
  proteins surveyed so far, 1HD6 has the highest fractional H coverage.
- **STEREO_AMBIGUOUS very high (54 rows, all ambig=2).** Driven by
  systematically-unassigned β-methylene stereo pairs across the
  disulfide-rich peptide. No high-ambiguity (code-4) outliers like
  1HS5/ARG-19. Defer to downstream policy; record the count.
- **Re-referencing modest.** 37 RefDB corrections; cross-cutting BMRB-CSV
  pattern holds (37/37 → corrected).

## 1B1V/4292 — fourth reviewed protein

**Protein:** PSP1, plasmatocyte-spreading peptide from *Pseudoplusia
includens* (insect cytokine, 23-residue short peptide). PDB 1B1V deposited
1998-11-23; SYNONYM PSP1, ENGINEERED YES, no mutation. Single chain A.

**Sources of record:**
- TOML provenance (BMRB): [`1B1V_4292.bmrb.toml`](1B1V_4292.bmrb.toml).
- TOML provenance (RefDB): [`1B1V_4292.refdb.toml`](1B1V_4292.refdb.toml).
- Audit outputs: [`1B1V_4292.audit.json`](1B1V_4292.audit.json), [`1B1V_4292.audit.md`](1B1V_4292.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`. Starter
table sufficed.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 129 | 67 | 7 ONE_TO_MANY, 55 STEREO_AMBIGUOUS, 0 UNMATCHED | H 143/161 (89%); C/N/O/S 0 |
| RefDB | 40 | 40 | 0 UNMATCHED | H 40/161; C/N/O/S 0 |

**AI observations:**

- **Pure proton-only deposition.** BMRB title is "NMR Structure of PSP1…";
  no C, no N, no S rows. Coverage of CHARMM heavy atoms is 0/n by design;
  H coverage 89% is excellent for a short peptide.
- **Silent N-terminal handling (1HD6 strategy).** GLU-1 BMRB rows are HA,
  HB3, HG2, HG3 only — depositor omitted backbone amide H rather than
  deposit a single-H value at NH3+. Same physics as 1DV0/1HS5/1HD6,
  consistent with 1HD6's strategy. Two strategies coexist across our 4
  reviewed proteins (1DV0+1HS5 single-H deposit; 1HD6+1B1V silent omit).
- **STEREO_AMBIGUOUS proportion very high (55/129 = 43%).** All ambig=2.
  Highest fractional stereo-ambiguity rate of the 4 reviewed proteins
  (1DV0 4/434 = 1%, 1HS5 10/336 = 3%, 1HD6 54/240 = 23%, 1B1V 43%).
  Pattern suggests short-peptide H-only depositions from late-1990s era
  systematically left β/γ methylene stereo pairs unassigned.
- **2 CYS, 1 disulfide expected** (residues 7, 19). CHARMM names CYS;
  no anomaly fires. Consistent with 1HD6 disulfide handling.
- **No re-referencing.** RefDB has 40 rows, 0 with corrected ≠ original.
  Clean entry — RefDB adds no information beyond a smaller H subset of
  BMRB.

## 1HA9/5028 — fifth reviewed protein

**Protein:** MCoTI-II (Momordica cochinchinensis Trypsin Inhibitor II), a
**macrocyclic knottin peptide** — naturally cyclic at the backbone. The
1HA9 linear PDB representation starts at a permutation point of the cyclic
backbone. PDB deposited 2001-04-02; SYNONYM MCOTI-II, OTHER_DETAILS SQUASH
TRYPSIN INHIBITOR. 34 residues, 6 CYS (3-disulfide cystine knot).

**Sources of record:**
- TOML provenance (BMRB): [`1HA9_5028.bmrb.toml`](1HA9_5028.bmrb.toml).
- TOML provenance (RefDB): [`1HA9_5028.refdb.toml`](1HA9_5028.refdb.toml).
- Audit outputs: [`1HA9_5028.audit.json`](1HA9_5028.audit.json), [`1HA9_5028.audit.md`](1HA9_5028.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`. Starter
table sufficed.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 186 | 129 | 12 ONE_TO_MANY, 45 STEREO_AMBIGUOUS, 0 UNMATCHED | H 210/227 (92%); C/N/O/S 0 |
| RefDB | 59 | 59 | 0 UNMATCHED | H 59/227; C/N/O/S 0 |

**AI observations:**

- **Cyclic-derived topology with non-zwitterionic termini** (NEW
  N-terminal strategy, third one observed). CHARMM SER-1 has a single HN
  (not the linear-peptide H1/H2/H3 NH3+) and CHARMM GLY-34 has a single O
  (not OT1/OT2 carboxylate). The MD prep treated MCoTI-II as a linear
  permutation of the cyclic backbone with mid-chain chemistry at both ends.
  Consequence: BMRB SER-1 backbone H 8.85 ppm matches cleanly to CHARMM
  HN — neither single-H-deposit anomaly (1DV0/1HS5) nor silent-omit
  (1HD6/1B1V) applies. **Three N-terminal strategies now observed across
  the 5 reviewed proteins.**
- **Pure proton deposition.** BMRB title is "Solution Structure of
  MCOTI-II"; 186 rows all H, no C, no N. H coverage exceptional (92%).
- **6 CYS handled uniformly** by CHARMM CYS naming (no CYX); no audit
  anomaly. Cystine-knot disulfide topology is in the structure, not the
  atom-name scheme.
- **STEREO_AMBIGUOUS 45 rows (24%), all ambig=2.** Intermediate rate;
  consistent with disulfide-rich peptide H-only deposition.
- **No re-referencing** (0 RefDB corrections). RefDB adds no information.

## 1G26/4656 — sixth reviewed protein (first HSP exemplar)

**Protein:** Granulin A N-terminal subdomain — 31-residue designed
disulfide-stabilised β-hairpin (cytokine-family, PDB FRAGMENT "N-TERMINAL
DOMAIN (RESIDUES 1-31)", ENGINEERED YES, MUTATION YES; SYNONYM HGA).
PDB 1G26 deposited 2000-10-17.

**Sources of record:**
- TOML provenance (BMRB): [`1G26_4656.bmrb.toml`](1G26_4656.bmrb.toml).
- TOML provenance (RefDB): [`1G26_4656.refdb.toml`](1G26_4656.refdb.toml).
- Audit outputs: [`1G26_4656.audit.json`](1G26_4656.audit.json), [`1G26_4656.audit.md`](1G26_4656.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`. Starter
table sufficed; HSP→HIS alias resolved correctly.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 166 | 153 | 13 ONE_TO_MANY, 0 STEREO_AMBIGUOUS, 0 UNMATCHED | H 192/211 (91%); C/N/O/S 0 |
| RefDB | 55 | 55 | 0 UNMATCHED | H 55/211; C/N/O/S 0 |

**AI observations:**

- **First HSP exemplar — alias machinery worked, protonation untestable.**
  CHARMM has HSP at residue 3 (doubly-protonated; +1 charge; ring atoms
  HD1, HD2, HE1, HE2). All 6 BMRB HIS-3 rows (H, HA, HB2, HB3, HD2, HE1)
  bound correctly to CHARMM HSP atoms via the `[residue_name_aliases]
  HSP = "HIS"` rule, firing 6 RESIDUE_ALIAS_APPLIED. **However, BMRB did
  NOT deposit HD1 (Nδ1-H) or HE2 (Nε2-H)** — exactly the two ring N-Hs
  that would distinguish HSP/HSD/HSE protonation. CHARMM's HSP choice is
  a prep-time decision with no NMR evidence either way; document-and-carry,
  downstream policy decides if/how to treat the protonation as a free
  variable.
- **Ambiguity codes entirely absent.** All 166 BMRB rows have empty
  `ambiguity_code` field — depositor did not deposit ambig codes at all.
  Consequently 0 STEREO_AMBIGUOUS fires; this is NOT evidence of clean
  stereo assignment, it is evidence that the depositor did not annotate
  ambiguity. Distinct from 1HD6/1B1V which deposited ambig=2 systematically.
- **VAL-1 silent N-term (1HD6/1B1V family).** BMRB row 1: HA, HB, HG1, HG2;
  no backbone H. CHARMM has H1/H2/H3 (NH3+); nothing fires.
- **6 CYS at 4, 10, 16, 17, 26, 27 — likely 3 disulfide bridges** stabilising
  the β-hairpin stack. CHARMM CYS uniform; no anomaly.
- **No re-referencing** (0 RefDB corrections).

## 1CBH/192 — seventh reviewed protein (oldest BMRB deposition)

**Protein:** Cellobiohydrolase I C-terminal domain — cellulose-binding
module from *Trichoderma reesei*. 36-residue, ENGINEERED YES, EC 3.2.1.91
(cellulase). PDB 1CBH deposited 1989-05-30 — **oldest deposition in the
calibration set**. BMRB id 192 — single-digit, indicating very early BMRB
entry.

**Sources of record:**
- TOML provenance (BMRB): [`1CBH_192.bmrb.toml`](1CBH_192.bmrb.toml).
- TOML provenance (RefDB): [`1CBH_192.refdb.toml`](1CBH_192.refdb.toml).
- Audit outputs: [`1CBH_192.audit.json`](1CBH_192.audit.json), [`1CBH_192.audit.md`](1CBH_192.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`. Starter
table sufficed.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 188 | 147 | 15 ONE_TO_MANY, 26 STEREO_AMBIGUOUS, 0 UNMATCHED | H 218/236 (92%); C/N/O/S 0 |
| RefDB | 63 | 63 | 0 UNMATCHED | H 63/236; C/N/O/S 0 |

**AI observations:**

- **Oldest BMRB entry — clean despite age.** Despite 1989 deposition era,
  the entry conforms to modern conventions for atom naming (post-2002
  IUPAC matching the starter): backbone H, methylene HB2/HB3, methyl
  pseudos. No legacy convention surprises caught by the audit. Either
  the entry was updated post-deposition or the depositor anticipated
  modern conventions.
- **HSP-4: same untestable-protonation pattern as 1G26-HSP-3.** BMRB
  deposits H, HA, HB2, HB3, HD2 (6.82 ppm), HE1 (8.56 ppm); ring N-Hs
  (HD1, HE2) absent. CHARMM HSP choice unsupported by NMR data.
- **THR-1 silent N-term.** BMRB rows HA, HB, HG2 only; CHARMM has H1/H2/H3
  NH3+. No anomaly fires — joins 1HD6/1B1V/1G26 silent-omit family.
- **Stereo ambiguity 26/188 = 14%, all ambig=2.** Depositor DID deposit
  ambig codes (156 code-1, 32 code-2 across 188 rows). Distinct from
  1G26 which deposited none.
- **4 CYS** at 8, 19, 25, 35 — likely 2 disulfide bridges in the
  cellulose-binding wedge fold. CHARMM CYS uniform; no anomaly.
- **No re-referencing** (0 RefDB corrections).

## 1I2V/4976 — eighth reviewed protein (second HSP exemplar with no ambig codes)

**Protein:** Heliomicin-LL — antifungal/antibacterial MUTANT of heliomicin
(insect defensin antimicrobial peptide). 44-residue, ENGINEERED YES,
MUTATION YES; PDB 1I2V deposited 2001-02-12.

**Sources of record:**
- TOML provenance (BMRB): [`1I2V_4976.bmrb.toml`](1I2V_4976.bmrb.toml).
- TOML provenance (RefDB): [`1I2V_4976.refdb.toml`](1I2V_4976.refdb.toml).
- Audit outputs: [`1I2V_4976.audit.json`](1I2V_4976.audit.json), [`1I2V_4976.audit.md`](1I2V_4976.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 237 | 219 | 18 ONE_TO_MANY, 0 STEREO_AMBIGUOUS, 0 UNMATCHED | H 273/296 (92%); C/N/O/S 0 |
| RefDB | 80 | 80 | 0 UNMATCHED | H 80/296; C/N/O/S 0 |

**AI observations:**

- **HSP-31 untestable** (third HSP exemplar with the same pattern). HD2
  (7.08 ppm), HE1 (8.49 ppm) deposited; HD1, HE2 absent. Cumulative
  pattern across 1G26/1CBH/1I2V: HIS protonation never resolvable from
  this depositor convention.
- **Ambig codes entirely absent** (237 rows, no codes). 0
  STEREO_AMBIGUOUS — silently ambiguous, not cleanly assigned. Joins
  1G26 in this depositor failure mode.
- **ASP-1 silent N-term** (1HD6/1B1V/1G26/1CBH/1I2V — five-protein
  silent-omit cluster).
- **6 CYS at 7, 18, 22, 32, 40, 42** — likely 3 disulfide bridges in the
  cysteine-stabilised α/β defensin fold. CHARMM CYS uniform.
- **No re-referencing** (0 RefDB corrections).

## 1I8X/4351 — ninth reviewed protein (HSP + moderate re-referencing)

**Protein:** Carp Granulin-1 N-terminal 30-residue fragment — CG1 1-30
peptide (cytokine granulin family, structurally analogous β-hairpin
stack to 1G26). Determined via semi-automatic ARIA NOE assignment.
PDB 1I8X deposited 2001-03-16.

**Sources of record:**
- TOML provenance (BMRB): [`1I8X_4351.bmrb.toml`](1I8X_4351.bmrb.toml).
- TOML provenance (RefDB): [`1I8X_4351.refdb.toml`](1I8X_4351.refdb.toml).
- Audit outputs: [`1I8X_4351.audit.json`](1I8X_4351.audit.json), [`1I8X_4351.audit.md`](1I8X_4351.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 204 | 151 | 15 ONE_TO_MANY, 38 STEREO_AMBIGUOUS, 0 UNMATCHED | H 181/200 (91%); C 53/141 (38%); N 0/33 |
| RefDB | 104 | 104 | 0 UNMATCHED | C 53/141; H 51/200 |

**AI observations:**

- **HSP-3 untestable** (granulin pair to 1G26 — both have HSP at residue
  3). HD2 (7.30 ppm), HE1 (8.64 ppm); HD1 and HE2 absent. Adds CA (56.48)
  and CB (30.08) shifts not seen in the other HSPs. 8 RESIDUE_ALIAS_APPLIED
  fired (6 H + 2 C rows, all on HIS-3).
- **First H+C deposition** in the latter half (since 1HD6 in the first
  half). 53 of 141 CHARMM C atoms covered (≈ Cα + Cβ for each residue).
- **VAL-1 silent N-term** (1HD6/1B1V/1G26/1CBH/1I2V/1I8X — six-protein
  silent-omit cluster).
- **STEREO_AMBIGUOUS 38 rows (19%, all ambig=2).** Depositor deposited
  ambig codes properly (166 code-1, 38 code-2 across 204 rows).
- **4 CYS at 4, 10, 16, 26** — likely 2 disulfide bridges.
- **BMRB-CSV-is-corrected pattern verified at scale.** Of 53 re-referenced
  rows, 53/53 BMRB shift_value match RefDB corrected (0 match original).
  Fourth confirmation (1DV0=125/125, 1HS5=96/96, 1HD6=37/37, 1I8X=53/53)
  — pattern is now robust.

## 1BWX/3449 — tenth reviewed protein (most HSP-heavy, MD-vs-BMRB length gap)

**Protein:** Human Parathyroid Hormone fragment PTH(1-39) — flexible
peptide hormone, no disulfide bridges. PDB 1BWX deposited 1998-09-29;
ENGINEERED YES.

**Sources of record:**
- TOML provenance (BMRB): [`1BWX_3449.bmrb.toml`](1BWX_3449.bmrb.toml).
- TOML provenance (RefDB): [`1BWX_3449.refdb.toml`](1BWX_3449.refdb.toml).
- Audit outputs: [`1BWX_3449.audit.json`](1BWX_3449.audit.json), [`1BWX_3449.audit.md`](1BWX_3449.audit.md).

**Adjudication applied:** none — `[residue_seq_offset] offset = 0`.

**Audit outcome:**

| Database | Rows | Matched (incl. anomalies) | Other status | Coverage |
|---|---|---|---|---|
| BMRB | 225 | 131 | 20 ONE_TO_MANY, 74 STEREO_AMBIGUOUS, 0 UNMATCHED | H 265/330 (80%); C/N/O/S 0 |
| RefDB | 66 | 66 | 0 UNMATCHED | H 66/330; C/N/O/S 0 |

**AI observations:**

- **MD topology longer than BMRB deposition.** PDB FRAGMENT "RESIDUES
  1-39" → CHARMM 39 residues; BMRB 3449 title is "PTH(1-34)" → BMRB
  rows only on residues 1-34. CHARMM residues 35-39 (VAL ALA LEU GLY
  ALA) are BMRB-silent: no rows, no anomaly, just 5 residues with no
  shifts to bind. Same pattern as 1HS5 (1-residue gap), here 5 residues.
  Downstream H/C predictions on 35-39 will lack experimental anchors.
- **Three HSP residues (9, 14, 32) — all untestable.** 18
  RESIDUE_ALIAS_APPLIED total (6 rows × 3 HSP). Each HSP gets H, HA,
  HB2, HB3, HD2 (~7.2-7.3), HE1 (~8.6); ring N-Hs (HD1, HE2) absent
  for all three. CHARMM HSP choice is a prep-time decision unsupported
  by NMR data on all three. Cumulatively across 5 HSP exemplars
  (1G26/1CBH/1I2V/1I8X/1BWX), every HSP follows the same
  HD2+HE1-only pattern. The depositor convention is universal for
  this set.
- **No CYS, no PRO** — first reviewed protein without disulfide-rich
  topology. Distinct from the 6 disulfide-rich peptides reviewed
  (1HD6, 1B1V, 1HA9, 1G26, 1CBH, 1I2V, 1I8X all have 2-6 CYS). PTH
  is a flexible soluble hormone.
- **STEREO_AMBIGUOUS 74 rows (33%, all ambig=2).** Highest absolute
  count of any reviewed protein; second-highest fractional rate after
  1B1V (43%). Depositor put ambig codes properly.
- **SER-1 silent N-term** (1HD6/1B1V/1G26/1CBH/1I2V/1I8X/1BWX —
  seven-protein silent-omit cluster).
- **No re-referencing** (0 RefDB corrections).



First-pass audit outputs were deleted in the mid-clean; the
per-protein tables still need to be derived. These observations
are from first-pass data and `SESSION_2026-04-19.md`; they MUST be
re-verified after each per-protein table is derived and the audit
re-run. Treat as hypotheses, not findings.

## Cross-cutting patterns (first-pass)

To verify per protein during review:

1. **CHARMM `HN` → BMRB `H`** universal for backbone amide. Encoded
   in starter as `[backbone] HN = "H"`; verify per-entry.
2. **Gly `HA1/HA2` → BMRB `HA2/HA3`** universal. Verify per-entry.
3. **Generic methylene `HB1/HB2`, `HG1/HG2`, `HD1/HD2`, `HE1/HE2`
   → BMRB 2/3 naming** (per-residue in starter). Verify per-entry
   that BMRB didn't use 1/2 (older convention) or primed names
   (HB, HB').
4. **CHARMM ILE `CD` → BMRB `CD1`** (no-locant vs locant divergence).
   Caught during 1HS5's first pass before it silently broke C-shift
   comparison. Verify any ILE-containing entry.
5. **Methyl pseudo-atoms** (THR HG2, VAL HG1/HG2, LEU HD1/HD2, ILE
   HG2+HD1, ALA HB, MET HE, LYS HZ): BMRB typically collapses to
   one pseudo; CHARMM has three explicit Hs. Flagged as
   ONE_TO_MANY; expansion deferred.
6. **HSP → HIS alias** fires `RESIDUE_ALIAS_APPLIED` on every HSP
   residue's rows. Protonation state not determinable from BMRB
   alone; requires deposition metadata check.
7. **N-terminal handling — three strategies observed across 5 reviewed
   proteins.** Both BMRB depositor convention AND CHARMM topology
   preparation contribute:
   - **Single-H deposit + CHARMM NH3+:** 1DV0 (GLN-3≡CHARMM-1), 1HS5
     (ASP-1). BMRB deposits a single `H`; CHARMM has H1/H2/H3. Mismatch
     fires UNMATCHED document-and-carry.
   - **Silent omit + CHARMM NH3+:** 1HD6 (ASP-1), 1B1V (GLU-1). BMRB
     omits backbone H entirely; CHARMM still has H1/H2/H3 NH3+. No
     anomaly fires (nothing to bind).
   - **Cyclic-derived non-zwitterionic CHARMM:** 1HA9 (MCoTI-II SER-1).
     MD prep treated cyclic peptide as linear with mid-chain chemistry
     at both ends — single HN at SER-1, single O at GLY-34. BMRB single
     `H` deposit matches cleanly; no anomaly. Distinct from a "linear
     peptide that happens to be acetylated/amidated" — there are no
     blocking groups, just no zwitterion.
   Implication: BMRB N-terminal H presence/absence is *not* a chemistry
   signal. The CHARMM N-terminal hydrogen count depends on the prep
   choice (zwitterionic vs cyclic-origin), independent of the deposit.
   Downstream must read both the BMRB row AND the CHARMM topology to
   know what is being matched.
8. **ARG ambiguity_code = 4 on β/γ methylenes** (1HS5/ARG-19): higher
   than the routine ambig=2 for amide NH2 pairs. Means depositor was
   unsure which intra-residue proton is which (vs 2 = stereo-pair
   unassigned). Worth distinguishing in any downstream swap policy.
9. **Disulfide-rich peptides keep CHARMM CYS naming.** 1HD6 has 6 CYS
   (likely 3 disulfide bridges); 1B1V has 2 CYS (likely 1 bridge). CHARMM
   topology files name all CYS as CYS regardless of disulfide state (no
   CYX as in AMBER); BMRB matches identity. No alias needed.
10. **STEREO_AMBIGUOUS proportion is a depositor-discipline signal,
    not a chemistry signal.** Reviewed-protein rates: 1DV0 1% (45 aa,
    2000), 1HS5 3% (34 aa, 2001), 1HA9 24% (34 aa, 2001), 1HD6 23%
    (37 aa, 2000), 1B1V 43% (23 aa, 1998), **1G26 0% (31 aa, 2000)**.
    All non-zero rates are ambig=2; 1G26 has *no ambiguity codes
    deposited at all* (empty field on every row). 0 STEREO_AMBIGUOUS
    in 1G26 is NOT evidence of clean stereo assignment — it is evidence
    that the depositor did not annotate ambiguity. Two failure modes
    silently coexist: ambig deposited honestly (1HD6/1B1V/1HA9) vs ambig
    not deposited at all (1G26). Downstream cannot infer assignment
    confidence from ambig-code absence.
11. **HSP/HSD/HSE protonation almost never testable from NMR alone.**
    1G26's HIS-3 deposition includes ring CARBON Hs (HD2, HE1) but
    omits ring NITROGEN Hs (HD1, HE2). The N-Hs are exchangeable in
    solution and rarely resolved at neutral pH. Without HD1 or HE2
    shifts, HSP (doubly protonated) vs HSD (Nδ1-H only) vs HSE (Nε2-H
    only) cannot be distinguished. CHARMM prep-time choice (HSP) is
    untested. Record per protein; defer adjudication to downstream
    sensitivity analysis (treat HIS protonation as free or fix by
    pH/structure).

## Next work

**All 10 calibration proteins reviewed (2026-04-19 → 2026-04-20).** Each
has `<pdb>_<bmrb>.bmrb.toml`, `<pdb>_<bmrb>.refdb.toml`, audit JSON, audit
markdown, and a section in this SUMMARY. TRACKING.md table marks each
reviewed.

Implementation-side: no further tool changes are known to be needed. The
starter `_starter_translation.toml` covered every name-mapping case
encountered; only one adjudicating TOML edit was required across the 10
proteins (1DV0's `[residue_seq_offset] offset = -2`).

The audit outputs are now the durable record. When downstream analysis
needs to reason about a protein's BMRB/RefDB binding to CHARMM atoms,
read its `<pdb>_<bmrb>.audit.md` first.
