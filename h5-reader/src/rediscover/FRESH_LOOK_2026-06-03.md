# Fresh-Look Validation — rediscover analysis program (2026-06-03)

Read-only safety-net review for the next session. No code or commits changed; this is
the only file written. Branch `h5-reader-pysr-spike`, HEAD `00ec168`. I verified the
2026-06-03 SESSION HANDOFF block in `STATE.md:7-46` against the committed code, the
result files in `/tmp`, the contract/spec docs, and the named memory files — trusting
nothing, re-deriving where I could. Framing throughout: "expected relationship +
probability + fit," never "law." Cited as `file:line`.

Bottom line: the arc is **solid and internally coherent**, the v1 substrate is real and
gated, and almost every cited number checks out against a result file. Two things to fix
before proceeding: (1) one IMMEDIATE-STATS number (the MOPAC-field r-values) has **no
locatable backing result file** in this run set; (2) the `ALLATOM_FIT_SPEC` has now
**landed** (the handoff still says "in flight") and is the next human-in-loop review.

---

## SOLID (verified true vs code / commits / result files)

**v1 substrate committed at `00ec168`.** `git show --stat 00ec168` confirms: adds
`PerAtomSubstrate.{h,cpp}` (1658-line cpp), wires `main_extract.cpp`, and extends
`python/nmr_extract/_catalog.py` (+32). Files match `STATE.md:10`.

**558,360 rows = 846 × 660.** Arithmetic checks (`846*660 = 558360`). The emitted run
`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-v1-fixed` manifest declares
`n_atoms: 846`, `n_dft_frames: 660`, `normalization: raw_lab_frame`. Row alignment
`row_id == frame_slot * n_atoms + atom_index` is declared in
`PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md:45,93` and `ALLATOM_FIT_SPEC_2026-06-03.md:14`,
and the positive-control report states it held on the emitted file
(`/tmp/per_atom_positive_control_analysis/report.md:4`). Lean: the run is ~1.05 GB
(571 MB of that is the separable AIMNet2 embedding), consistent with the contract's
~650 MB / ~80-MB-without-embedding sizing (`NODE_STORE_CONTRACT_2026-06-02.md:105-112`).

**The four reducer-generalizations are in the v1 code, in lab frame.**
- ring Johnson-Bovey: `JohnsonBoveySourceUnitKernelLocal` + emitted `ring_jb_T0`,
  `ring_jb_T2_0..4` (`PerAtomSubstrate.cpp:473,919-922`).
- charge q/r³: `charge_q_over_r3` + `charge_q_over_r3_T2_0..4`
  (`PerAtomSubstrate.cpp:525,923-925`).
- FF14SB field: same included-charge source set as q/r³ (`PerAtomSubstrate.cpp:690`).
- valid-source McConnell: `mc_lit_valid` with the producer self/near-field filter
  (`mc_source_is_self_or_bonded`, `PerAtomSubstrate.cpp:578,582,588,595`).
- frame stamped `raw_lab_frame` (`PerAtomSubstrate.cpp:1304,1576`).

**The backbone-reproduces-broad gate has its raw material in v1.** The substrate emits a
`per_atom_substrate_backbone_audit` array carrying `broad_charge_literature_kernel_T2_*`,
`broad_field_local_*`, `broad_mc_lit_T2_valid_*` (`PerAtomSubstrate.cpp:971-984,1269`;
present in the run's `per_atom_substrate_column_specs.json`). These are the
broad-backbone-local-frame columns whose purpose is exactly "validation can compare
current all-atom backbone rows against the existing broad substrate"
(`PerAtomSubstrate.cpp:714-715`). The "machine precision" claim in `STATE.md:10` is a
gate the build session ran; the audit columns exist to support it, but the pass/fail
record itself is not a committed artifact in this tree (re-run the gate if you need the
number — see VERIFY FIRST).

**Combined-score per atom type (all-tier).** Handoff `STATE.md:15` says N 0.73, C 0.74,
O 0.67, CA 0.54, HN 0.53, HA 0.49. These match the `within_frameblock_test_R2` of the
"all emitted" tier in `/tmp/combined-mopac-layer3/combined_score_table.csv` (N 0.7327,
C 0.7361, O 0.6688, CA 0.5393, HN 0.5315, HA 0.4910). Accurate to rounding.

**AIMNet2 lift +0.03–0.10, biggest at CA.** Within-frameblock classical→all: N +0.052,
CA +0.095, C +0.065, O +0.036, HN +0.058, HA +0.037 (same CSV). Range and the CA peak
are correct.

**Positive-control network (the stool holds).** `report.md` confirms:
APBS-EFG↔MOPAC-Coulomb-σ sign-flip r=0.702 (handoff 0.70), chain ff14sb_E→charge_T2
0.944 (0.94), apbs_E→apbs_efg 0.760 (0.76), apbs_efg→mopac-σ 0.549 (0.55), FF14SB-vs-APBS
field divergence as expected screening
(`/tmp/per_atom_positive_control_analysis/report.md:13,23-28,12`). All match.

**Per-mechanism standing (corrected-backbone capstone).** `CORRECTED_BACKBONE_SNAPSHOT.md`
confirms: ring γ_bare −11.326±3.668 nA/T, Pople/JB pooled γ 0.923±0.531
(`:8,28`); charge q/r³ N γ 9.701±0.227, absT2 r 0.776 (`:30`); McConnell valid-source
can't-make-it-work-for-now (`:29`); AIMNet2 ceiling C-within 0.717 / N-within 0.587
(`:48,34`); EFG local-frame ≈0 (`:33,40,47`). The charge PySR R²=0.775 figure is
recorded in `STATE.md:18,106` and the capstone narrative; it is the strongest symbolic
fit claim and is consistent with the standing.

**HN field "between R² 0.48".** The capstone records FF14SB HN between R² 0.48
(`STATE.md:18,110`); the parallel `static_environment_calibration.md:33` (a different,
mopac-field-dir run) shows HN ff14sb_charge_field_sigma between_LOAO_R2 0.5195. Same
mechanism, two runs, ~0.48 vs ~0.52 — within run-to-run noise, not a contradiction, but
note the two sources differ.

**The analysis arc is grounded in memory and consistent across docs.** The FIT-ALL-THEN-
PARTITION reframe in `STATE.md:27-31` (step 1) matches `project_law_example_hunter.md`
verbatim (the 2026-06-02 REFRAME, `:167-192`: fit across all ~846 atoms, between-axis
determinable because 846 ≫ features, partition on input conditions, happy cases emerge).
Step 4 (grade by statistical position) matches `feedback_law_as_statistical_position.md`
(PySR Pareto / Bayesian-IC / null; position not pass/fail). Steps 2-3 (between-calculator
network, `equations/<mechanism>/` pre-registration) match `project_between_calculator_network`
and the `project_law_example_hunter` output-structure block (`:109-121`). All seven named
POINTERS memory files exist.

**The workflow loop is stated consistently.** `STATE.md:33-39` — spec→vet→execute→
drift-assessment+postmortem; human-in-loop each loop; codex grinds / lead vets spec +
judges postmortem / lead context is irreplaceable; breakout at internal-check boundaries;
"expected relationship + probability + fit," never "law," never "ceremony." This matches
`feedback_token_economy_codex_codes`, `feedback_law_as_statistical_position`, and
`feedback_dont_ai_the_definition_of_a_law`. I found no stale "law-gating" or "model-not-
law stamp" framing left in the handoff block (the older `combined_score_report.md:1`
still opens "Model-not-law:", but that is an older result file's header, not the governing
arc — harmless, not load-bearing).

---

## DRIFTED / STALE / NEEDS A NOTE (not errors of substance, but flag them)

**1. `ALLATOM_FIT_SPEC_2026-06-03.md` has LANDED — the handoff still says "IN FLIGHT."**
`STATE.md:12` reads "IN FLIGHT: codex drafting `ALLATOM_FIT_SPEC_2026-06-03.md`." The
file now exists, dated 2026-06-03, **status "draft for review only"**
(`ALLATOM_FIT_SPEC_2026-06-03.md:3`). So the drafting is done; the actual next action is
the human-in-loop **spec review**, not waiting for it to land. The spec is well-formed: it
correctly scopes Piece 1 (the small charge-scalars emit that closes the positive-control
gap — `ff14sb_charge`, `mopac_welford_mean_charge`, `:25-30,67-89`) and Piece 2 (the
all-atoms joint fit + partition on the v1 substrate, reading one directory, NOT repeating
the `/tmp/combined-mopac-layer3` multi-directory join — `:32-43,45-57`). It targets v1 at
`00ec168` explicitly. This is the right next chunk and matches the arc.

**2. The IMMEDIATE-STATS "combined score" is from the OLD multi-directory workflow, not
the v1 substrate.** `/tmp/combined-mopac-layer3/run_audit.json` shows that run joined
THREE directories — a 500-frame capstone (`/tmp/rediscover-corrected-backbone-snapshot-1p9j`),
a 660-frame MOPAC dir (`/tmp/rediscover-task1-mopac-field`), and an h5 intersection of
250 — by `(atom_index, original_index)`, landing on 500 joined frames and the
54-atom-per-stratum backbone slices (`run_audit.json:2-6,63-67,97-102`). The v1 substrate
spec itself says this is "the immediate consumer prototype" that v1 is meant to **replace**
(`PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md:26-30`). So the headline combined-score and
per-mechanism standing numbers were produced by the workflow the v1 substrate retires —
they are valid as the prior reading, but the all-atoms fit on the actual 558,360-row v1
substrate is the IN-FLIGHT/next step and may move them. The handoff is honest about this
(it labels the all-atoms fit "the fix" for the unstable between-axis, `STATE.md:15`), but
a reader could mistake the combined-score for a v1-substrate result. It is not.

**3. "Combined score" is the WITHIN axis; the BETWEEN axis is weak/unstable.** The
handoff numbers are `within_frameblock_test_R2`. The companion `between_LOAO_test_R2` in
the same table is weak-to-negative at 54-atom strata (N classical 0.095, C classical
−0.161, O 0.070; `combined_score_table.csv`). The handoff says this plainly
("Between-axis was unstable at 54-atom strata → the 846-atom all-atoms fit is the fix",
`STATE.md:15`), so this is not a hidden problem — but the word "combined score" alone
omits that it is the within-axis read. The next session should keep the two axes labeled
distinctly when re-running on v1.

**4. STATE.md top line is dated "2026-05-31, EOD" while the live handoff is 2026-06-03.**
`STATE.md:1` vs `:7`. Cosmetic — the file accretes handoffs newest-first — but the stale
top-of-file date can mislead a tired reader about freshness.

---

## REAL GAP / COULD-NOT-VERIFY

**The MOPAC-field r-values have no locatable backing result file.** `STATE.md:16` states
"MOPAC-Coulomb-EFG carries the field signal — N r 0.505 (R² 0.33), C 0.579 (0.24)
strongest; HN 0.36; CA/HA weak; O absent." I grepped every result file under
`src/rediscover/`, `/tmp/rediscover-task1-mopac-field/`, `/tmp/combined-mopac-layer3/`,
and `/tmp/per_atom_positive_control_analysis/` for `0.505` and `0.579`. The ONLY hits are
`STATE.md` itself and `APBS_VALUES_VERIFICATION.md` (and the latter is a false match — a
PubMed Central URL `PMC5600505`, `:35`). The result files I did find tell a different
picture for the BETWEEN axis: `static_environment_calibration.md:47-52` shows
`mopac_coulomb_shielding_T2` between-axis r NEGATIVE for N/C/O (−0.107, −0.677, −0.613),
and `efg_localframe_audit.md` shows local-frame EFG ≈0. So the 0.505/0.579 are plausibly
**within-axis** field-carries-signal reads from a run whose report did not survive into
this tree (the MOPAC-field dir's `variance_decomposition`/`efg_localframe_audit` don't
contain them either). This is the one IMMEDIATE-STATS number a fresh reader cannot
reproduce — treat it as unverified until re-derived. Everything else in the stats block
reproduced.

---

## VERIFY FIRST (before proceeding with the all-atoms fit chunk)

1. **Review `ALLATOM_FIT_SPEC_2026-06-03.md` (human-in-loop) — it is the real next gate,
   not "in flight."** It is a clean two-piece chunk with the right boundary: Piece 1
   closes the charge-scalars gap the positive-control flagged (v1 emits only the AIMNet2
   scalar charge; FF14SB partial charge and MOPAC Welford-mean charge are not emitted —
   `summary.json:30-37`, `report.md:11`); Piece 2 is the FIT-ALL-THEN-PARTITION on v1.
   Confirm the gate sequence (`:45-57`) and that Piece 2 must report "no charge-complete
   substrate" rather than silently re-joining external files (the no-multi-dir-merge rule).

2. **Re-derive the MOPAC-field r-values (the 0.505/0.579 gap above) before citing them.**
   They are the basis for "the field signal is the MOPAC leg, not APBS" (`STATE.md:16,66`).
   The positive-control's `abs_apbs_E→abs_mopac_coulomb_T2 r=0.538` and the negative
   between-axis mopac rows suggest the within-axis is where this lives — confirm on v1.

3. **The APBS-EFG↔MOPAC sign-flip is a CONVENTION, confirm not a bug.** The positive
   control reads raw comp r=−0.702 and calls it "expected because APBS is raw EFG and
   MOPAC column is shielding-sign T2" (`report.md:13`). v1 records `sign_convention`
   per-column (`PerAtomSubstrate.cpp:1206`). Before the network edge is graded, confirm
   the sign relationship is the documented EFG-vs-shielding convention, not an unintended
   flip — the handoff itself flags this as open (`STATE.md:11`).

4. **Re-run the backbone-reproduces-broad gate if you need its number.** The audit columns
   are in v1 (`per_atom_substrate_backbone_audit`), but the machine-precision pass is a gate
   the build session ran, not a committed artifact here. Oracle parity is the explicit
   command in `NODE_STORE_CONTRACT_2026-06-02.md:142-147` (`--case all` covers ring+mc
   only; not broad/all-atom/efg/buckingham/aimnet2 — re-run scoped ctest for those).

5. **Keep within and between axes labeled distinctly** when the v1 all-atoms fit lands —
   the 846-atom fit is specifically meant to make the between (static) axis determinable
   where the 54-atom strata could not (`project_law_example_hunter.md:176-181`).

---

## ARC COHERENCE — no contradiction found

Fit → partition → network → equations → statistical-position is consistent across
`STATE.md:27-31`, `NODE_STORE_CONTRACT_2026-06-02.md` (the resident node store + lazy
pair-index, materialization-as-a-verb), `PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md` (the lean
per-atom aggregate face), `ALLATOM_FIT_SPEC_2026-06-03.md` (the charge-complete fit
chunk), and the memory files. The discipline (model-is-the-spine, no Python physics, T2
preserved end-to-end, no factories, nmr_extract /shared extractions sacred and only
rediscover substrate cleanup-eligible) is stated identically in the contract
(`:5-8,35-38`) and `STATE.md:41-44`. The one substantive caution is the MOPAC-field
number gap; the structural cautions are that the headline stats are pre-v1-workflow and
the next spec has already landed for review.
