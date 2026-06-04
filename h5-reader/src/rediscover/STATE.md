# Rediscover — current state (2026-06-04)

## STAGE 2 DONE (2026-06-04, `ecbddd1`) — per-mechanism law fits + unified D_ab-sum (fitter decomposed)

Chunk 1: the ~4k-line fitter DECOMPOSED (`allatom_fit_common/build3/legacy` + `stage2_law_fits`;
`piece2`=wrapper); **reproduce-gate GREEN — Build-3 score/partition artifacts BYTE-FOR-BYTE.** frozen
get_C (1.11e-16), 5-comp T2, anti-circular (clean strata by CaseHunter + input-side dominance);
conventions applied (ringχ excluded, WaterEFG→−Hessian, Larsen-ppm separate). Substrate build4; run
`/tmp/rediscover-runs/2026-06-04-stage2-fits` (165M). Gates lead-verified (no per-source dump,
Python-only no producer/C++/git, disk 6.7G/15G).

PER-KERNEL (clean strata; **statistical-position + bucket, NOT literature-match**):
- **charge q/r³ = RECOVERED-LAW** — coeff +9.30 [8.10,10.51]; within 0.28 / **LOAO 0.38** (strongest on
  the BETWEEN/static axis); paths agree; fixed D_ab shape supported. The cookie.
- **field (MOPAC-Coulomb) = RECOVERED-LAW** — coeff −0.83 [−1.23,−0.43]; weak (~0.03) but NONZERO, PySR
  agrees. Vindicates MOPAC>Amber.
- **ring (Pople) = form-recovered-scale-fitted** — coeff +0.69 [0.21,1.17]; within 0.28 / LOAO 0.17; REAL
  but THIN (5 atoms — aromatic-ring-facing H).
- **McConnell = can't-make-it-work standalone** (CI spans 0; home = joint fit).
- **H-bond = can't-make-it-work this cut** (CI spans 0; geometric + Larsen don't rescue).
- **UNIFIED D_ab-SUM** (through-space atoms, 25 atoms / 3.9k rows): **within R² 0.432 / LOAO 0.258** — the
  "calculators as shadows" COMBINE genuinely recovers through-space total-T2. Intensities charge +15.5,
  MOPAC-field −8.8 well-determined; mc/pq/disp huge/weak. Recovery REAL, calibration NOT literature-clean
  (correlate-not-match held; disagreeing shadows diagnostic, not averaged).
- 3 PATHS: ridge + PySR + equivariant-Schur (closed-form on reduced sums — NOT the full e3nn-per-source,
  which stays the deferred chewer); agree where the law holds.

### Stage 2.1–2.3 follow-ups (2026-06-04)
- **2.1 happy-spot sweep + frame-ablation** (`f92109c`): the signal pops with **MODULATION (driver worked),
  NOT isolation** — "loud, not isolated" (lead). Only **field** (modulation) + **unified** (dominance+modulation)
  rise toward clean; charge/ring/mc/hbond flat-or-fall; McConnell/H-bond NOT rescued at the clean end. Support
  THIN everywhere (one protein) → the sweep is power-starved → 720-WT gives it teeth. Frame-ablation: recovery
  survives at FEW frames (knees ring 20 / field 50 / unified 50 / charge-within 200; charge LOAO 0.40 holds at
  50fr) → **ubiquitin 1ns@20ps (50fr) keeps the main signal — cheap 2nd protein viable.**
- **2.2 unified vet** (`943915f`): **REAL combine, NOT charge-in-a-coat** — charge-alone within ≈0.007; the
  within recovery is carried by **MOPAC-field + McConnell** (drop-one: mopac_field +0.194, mc_cat4 +0.171),
  charge adds little. within 0.43 stable under shrinkage; **between-atom LOAO 0.26 is modest + N-limited**
  (26 terms / 25 atoms → 0.10 shrunk) — NOT "overfit-dismissed": a real-but-thin between signal, reportable
  with PROBABILITY framing; the 25-atom thinness → 720-WT is the cure.
- **REPORTING STANDARD (lead, 2026-06-04, now law):** grade by STATISTICAL POSITION + DETERMINABILITY, never
  binary survives/overfit. Two questions: (1) fairly-indicative-AND-not-an-artifact (highly reportable, even
  necessary, WITH caveats + probability framing); (2) robust-in-context. Scale: ~0.03 ≈ nothing; ~0.2 =
  potentially something OR trash, decided by whether we can DETERMINE what drives it; higher = clearly
  something. [[feedback_law_as_statistical_position]] [[feedback_transparent_cutoffs]].
- **2.3 probability close DONE** (`bp5cixi7k`; lead-verified, committed): permutation-null (1000×) +
  determinability + lead-scale, every recovery. VERDICT (within = 1P9J's probability axis; between =
  case-study → 720-WT): **charge within indicative+determinable** (R²0.28, z451); **charge LOAO indicative**
  (0.38, z33); **ring within indicative+determinable** (0.28, z155 — ring's probability lives on WITHIN; its
  LOAO 0.17/z3.2 is the thin case-study); **unified within indicative+determinable** (0.43, z263 —
  field+McConnell mixture); **unified LOAO above-null but INDETERMINATE** (0.26, z2.0 — atom-axis attribution
  undetermined → needs 720-WT). **field standalone = ~null (0.03-class)** — CORRECTS the Stage-2
  "field=recovered-law": coefficient nonzero but recovery null-class; field's real value is the TOP
  contributor to the COMBINE (drop-one +0.198), not a standalone law. **mc + hbond standalone = ~null**
  (CI spans 0). Cutoff: within p<0.001 + robust across the 0.5/0.7 sweep; LOAO fragile (p→0.12 at one cut) =
  the thin between axis. [[feedback_law_as_statistical_position]] [[feedback_transparent_cutoffs]].

CAVEATS: ONE protein; thin clean strata (ring 5, unified 25) → "1P9J across its structures," no population
inference; the two nulls are PROVISIONAL (this cut/data), not earned.

NEXT (open, lead to steer):
1. **HAPPY-SPOT SWEEP (the core-assumption test).** Stage 2 fit AT a clean threshold but did NOT sweep.
   The hunt proper = response curves recovery-vs-cleanliness (dominance / isolation / **a geometric-noise
   criterion we don't yet have** — CaseHunter gates isolation/motion/quiet, none = "low geometric noise";
   motion even pulls the other way), strict→loose: does recovery POP toward the cleanest spots ("noisy
   geometry limits visibility")? Could RESCUE McConnell/H-bond at THEIR happy spots. Cheap, existing substrate.
2. **720-WT STATICS PILOT.** Same r²SCAN + same `.out` files (CONFIRMED — absolute σ present). Lots of rings
   → fattens ring's thin BETWEEN axis + cross-protein charge validation (charge's LOAO 0.38 is its strongest
   axis = the static one the 720-WT set is made of). Needs the rediscover substrate emitted on the 720 WT
   structures (bounded static run). The transferability/statics pilot (arc layer 4).
3. **Frame-count ablation** — recovery vs n_frames (within-axis only; AR(1) ρ≈0.53 → effective ≪ raw). Cheap.
   **Doubles as the go/no-go for a fast 2nd-protein-WITH-dynamics run** (lead idea 2026-06-04): ubiquitin,
   one choice ns @ 20 ps ≈ 50 frames, DFT-cheap. 1P9J is already ~20 ps stride (751 fr / 15 ns) → subsampling
   1P9J to 50 frames SIMULATES it. If recovery survives at ~50 frames, ubiquitin-50 is viable (within-axis,
   2nd protein, +a few rings) and complements the static 720-WT. ORCA-budget: weigh vs the one-more-run-
   Trp-cage earmark ([[project_orca_budget_one_more_run_trpcage]]) — lead's call; ubiquitin-50 is a fraction
   of a full trajectory's cost.
4. **McConnell/H-bond → joint/ensemble fit** (their home, per the deep-audit).

---

## SESSION HANDOFF (2026-06-03 LATE) — all-atoms fit through Build 3 LANDED

ALLATOM_FIT_SPEC vetted + re-chunked to TWO human-in-loop loops (codex had set one
chunk for convenience; split to build-then-analysis). Both ran clean.

- **Loop 1 (Piece 1, C++ emit, `4bb9a01`):** `ff14sb_charge` +
  `mopac_welford_mean_charge` (+present) appended to `per_atom_substrate`,
  append-only. 5 emit gates green (NPY byte-parity 0-diff, value-identical existing
  CSV cols, backbone-audit byte-identical, tests 6/6). Judged-good drift: a 2-line
  additive `ReadScalarWelford` `n_per_atom` fallback in `QtTrajectoryH5.cpp` surfaces
  existing-but-unread MOPAC Welford charge (static-workaround pattern; values sane —
  87.7% FF14SB sign-agreement, divergence on near-zero aliphatics; also a latent
  reader fix). `charge_complete_rows = 558360` (100%). Charge-agreement positive
  control now live from one substrate.
- **Loop 2 (Piece 2, Python fit+partition, `5b5525b`):** all-atoms ridge on DFT T2
  (frozen get_C, 5-comp), 3 tiers on the SAME charge-complete rows, train-only PCA
  (32 comp) on the embedding. ALL fit-stage checks pass. Results EPHEMERAL:
  `/tmp/rediscover-runs/2026-06-03-allatom-fit-piece2/`; script committed.

RESULTS (held-out R², facts not verdicts):
- **Between/static axis became DETERMINABLE — the reframe's payoff.** N between
  **+0.60** (was +0.095 at 54-atom strata), within +0.69; aromatic_heavy between
  +0.61. Global classical: between +0.10, within +0.13.
- **AIMNet2 + charge tiers do NOT lift; they HURT between** (global between
  +0.104→+0.058→+0.053; within flat ~0.125). MOVED from the old 54-atom prior
  (+0.03–0.10 AIMNet2 within-lift). **CAVEAT: fixed alpha=10 on ~80 features over
  846 atom-means is likely under-regularized for the higher tiers → the negative
  delta is plausibly OUR regularization artifact, not "AIMNet2 carries no T2 signal";
  old win was on σ_iso (T0), here we fit T2.** OPEN: alpha-selection (inner-CV)
  sensitivity on the higher tiers before any AIMNet2-on-T2 conclusion
  (attribute-our-own-mistakes-first).
- **Emergent favourable partitions (anti-circular, the hunter payoff):** N at high
  |ring_jb_T2| Q5 between **0.81**; aromatic_heavy nearest-ring Q1 (closest) 0.79;
  N nearest-ring Q2 0.785 — rising with driver exercise.
- AR(1) N_eff ~37.5K (median lag-1 ρ 0.53). Thin: amide_sidechain 9 atoms, sulfur 7.

- **Loop 2b (alpha-selection + analysis-side AIMNet2 switch, `44e786c`):** train-only
  inner-CV alpha selection settled the confound — selected alpha **1000–10000** (we WERE
  under-regularized at 10). But correcting it did NOT rescue AIMNet2: global between
  `+AIMNet2` delta −0.046 → **−0.014** (shrank ~70%, still negative); within ~0. In every
  cleanly-fit stratum (N between 0.50 / within 0.69; aromatic_heavy 0.41/0.46; O 0.07/0.14)
  AIMNet2 adds ≈nothing. **AIMNet2 does not make it on T2.** CAVEAT (not asserted): its
  "magic" was on **σ_iso/T0**; coherent hypothesis = AIMNet2 carries the isotropic/local
  part, the T2 anisotropy is geometric and the classical through-space kernels already get
  it — a **T0 companion fit** would confirm the split. (H strata go wildly negative on the
  between axis — thin/LOAO-unstable + aromatic-H self-ring excluded; not interpretable.)
  "switchable" = ANALYSIS-side feature-tier toggle ONLY; producer/substrate AIMNet2 stays
  always-on (no flag/build/label) — [[feedback_aimnet2_required_no_weasel]] holds, refined.

CADENCE (2026-06-03, to protect the lead's finite context across the remaining loops):
codex writes a ≤40-line POSTMORTEM + prints only a short summary (NO diff dumps to
stdout); batch related work into one gated loop; checkpoint STATE richly (the WHY, not
just numbers) each loop so a fresh session resumes losslessly. Human-in-loop each loop
still holds.

CHANNEL-COMPLETION EMIT DONE — substrate is now COMPLETE (2026-06-03):
- **Loop 3 (`b583d7c`)** + **Loop 3b (3 gated chunks `9965c2b`/`1d0f56a`/`5a39288`,
  `CODEX_BRIEF_PIECE3B`)**: emitted the full cheap-wire menu, C++ spine only (0 .py;
  guard-confirmed), additive, oracle-parity byte-identical, ~3.0 GB (431 new cols).
  Landed: `T0/T1/T2/dia/para` targets (**T1 frame-VERIFIED** ~1e-4°, a completeness
  DIAGNOSTIC not a shift claim); 4 ring-current paths (bs/hm/ringχ/jb + per-type, total_B,
  counts); pq/disp per-type; McConnell per-category (mc + mopac + bond_orders); the full
  field/EFG agreement set (mopac field+EFG, aimnet2 EFG, water EFG); eeq_cn; H-bond
  best+backup (Larsen per-class 1pHB/2pHB/1pHaB/2pHaB donor/acceptor-resolved + hbond_scalars
  + DSSP chemical-flag & raw backup); conditioners (ss8/chi/omega/pyramidalization/ring_geom).
  fail-loud-AND-LOCATE worked: "absent" hbond geometry located in hbond_scalars + DSSP.
  Historical run: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-piece3b-final`, later drop-old replaced by the Build1 substrate.
- FLAGS: stray commit `91db21c` (doc/calculators/CONTINUE_PROMPT.md) interleaved — NOT 3b
  (likely a concurrent docs-effort codex; two writers on the repo). `/tmp/rediscover-runs`
  was brought back under budget by keeping the Build1 substrate and analysis result dirs.
- GATE for all this: emit keyed (atom,frame)-only, C++ derives/reduces in memory, lean DISK
  (<15 GB output) but generous RAM (swap OK — patient, not sloppy). [[feedback_all_statistics_minimize_python_15gb]]

ADVERSARIAL REVIEW of the 3b emit (opus, 2026-06-03) — **substrate SAFE TO FIT, no
critical/high.** Verified: T0/T1/T2 round-trip; dia+para=total (1e-3 = ORCA print round);
per-type sums → aggregates at machine precision; EFG set shares units/convention/frame
(aimnet2_efg vs mopac_coulomb_efg +0.967); mopac_bond_orders 52% = clean topology split
(zero partial); hbond_scalars order proven (s₀=nearest_dist, |s₁−1/s₀³|=5.5e-17); frame +
offset integrity clean. HANDLE IN NEXT ANALYSIS:
- **MEDIUM — ringχ is opposite-convention.** `ringchi_shielding` anti-correlates bs/hm
  (−0.72): producer decomposes bare susceptibility χ/r³ WITHOUT the shielding minus-sign,
  and it's in Å⁻³ not ppm·T/nA. bs/hm/jb are the clean comparable set (+0.994). Sign-flip +
  rescale ringχ before any cross-method slope, OR report it as a different-convention path —
  NEVER naive-slope ringχ vs bs/hm (else a false "validation fail"). Producer-side; handle
  in analysis (no re-emit). (Optional future: relabel ringχ sign_convention at the producer.)
- **LOW** — `target_T1`/`bs_total_B` stored as Cartesian antisymmetric vector but labeled
  irrep `1o` (harmless; conventions-note). DSSP raw backup broadcasts donor+acceptor onto all
  4 peptide-plane atoms → the fitter must select by atom role (donor→HN/N, acceptor→O/C′).

BUILD 1 DONE (`a353104`, 2026-06-03) — the C++ PARTITION-FILTER tool. Pure C++ (CaseHunter.cpp
+ PerAtomSubstrate.cpp; 0 .py). Isolation primitives `gap_to_2nd_{ring,charge,bond}_r` +
`dominant_fraction_{ring,charge,mc}` (lean per-(atom,frame) scalar; query_frame_slots stayed
1 — no boa-constrictor) + bin-id columns + the `CaseHunter` (typed-habitat × frame-windows
over rowPairContributions; isolation∧motion∧quiet; DFT recovery MEASURED-not-selected-on —
anti-circular). Gates: oracle byte-parity (append-only 26→32 conditioning cols); **dominance
two-path cross-check PASS to ~5e-13** (new per-(atom,frame) dominant_fraction == 1-frame
named-query); hunter 24 candidates each ring/charge/mc, deterministic. Disk-guard honored
(full-path deletes only; drop-old replaced piece3b-final). Substrate is now
`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1`. (ORCA scratch self-cleans ~30 G
per round — transient, not accumulating.)

GRANULARITY NOTE (2026-06-03, Jessica's question before Build 2): the partition is categorical
insight into the KIND of quantum effect nearby atoms create. We HAVE the conditioner
granularity on both axes — source-kind (mechanism × type × geometry/orientation × isolation ×
ss8) and response-kind (T0/T1/T2/dia/para = the effect-type axis; dia/para+T1 are what let us
categorize the kind, not just magnitude). The LIMIT is DATA, not granularity: 1P9J is one
protein → fine categories are RICH within-atom (modulation over 660 frames) but THIN
between-atom (≈24 clean cases/mechanism, 16 rings). Robust between-atom categorical statistics
need the 2nd case (Trp-cage) / the 720-backbone transferability pilot. Build 2 must
support-flag (atom-count/N_eff per category), never over-claim a thin fine-categorical trend.

BUILD-1 ADVERSARIAL REVIEW (opus, 2026-06-03) — **anti-circularity HOLDS by construction**
(every selection/ranking input traced; none touches DFT; dft_recovery_R2 measured-after-only).
Bulk clean (score, windows, ring/bond gaps, dominance, bins, manifest, index-discipline). THREE
fixes found → FIX LOOP fired (`CODEX_BRIEF_BUILD1FIX`):
- **H1 (HIGH):** the ChargeSites cloud includes the target atom itself → `nearest_charge_r`
  collapses, `gap_to_2nd_charge` pinned ~1 Å, `bin_nearest/gap_charge` constant-0 (zero-info
  partition axes; charge→N is a headline). Pre-existing bug Build 1 exposed; dominance path
  already excludes self+same-residue, so gap/dominance disagree. Fix: exclude self+same-residue
  in gap/nearest for ChargeSites (re-values those columns — intentional re-baseline).
- **M1 (MEDIUM):** the `anti_circular_assertion` is a vacuous no-op (`dftTouchedDuringSelection`
  never set → hardcoded true, can't go red). Invariant holds today but the guard is fake. Fix:
  a REAL guard (DFT-read during selection throws) + a test that a deliberate leak makes it fire.
- **M2 (MEDIUM):** `"_r"` substring mislabels counts/magnitudes/fractions as units=Å. Fix:
  endsWith/explicit.
FIX LOOP DONE (`d9ba53d`, 2026-06-03) — all three fixed, surgically, pure C++:
- H1: charge bins now NON-DEGENERATE (`bin_nearest_charge_r` {0:511k,1:45k,2:1899,3:3};
  `bin_gap_to_2nd_charge_r` {0:273k,1:133k,2:145k,3:7277}; `nearest_charge_r` 1.23–8.23 Å).
  Charge→N is a usable partition axis again. gap/dominance now share the self+same-residue-
  excluded charge neighbour set.
- M1: `DftFrameSet` scoped selection-read guard — DFT read during CaseHunter selection THROWS;
  regression test proves it fires; `anti_circular_assertion: true` now backed by a real guard.
- M2: units corrected (only charge_excluded_same_residue_n, abs_ring_jb_T2, dominant_fraction_ring).
- Surgical parity: ONLY the charge-distance columns + bins + units metadata + charge cases-manifest
  changed; DFT targets + all other conditioners/mechanisms/bins/queries byte-identical. Disk-guard
  honored (full-path drop-old). Substrate (correct) = `…-build1`.

BUILD 2 DONE (`a2d69cb`) — the partition. Integrity clean (pure Python, lean 53M, all CV/
partition checks). **N total-T2 between 0.63 / within 0.81** (the win) + favourable habitats
(MC-magnitude monotone-rise on N; modulation thresholds) + 60 navigable hunter-intersected cases.
THREE reframing negatives: para-T2 ≈0 everywhere (overturns "para is cleaner" — gauge-dependent +
likely local/out-of-through-space-reach → total-T2 is the target); T1 ≈0/neg (field-linear not a
recoverable channel); new mechanisms (incl. H-bond) don't lift HN/O. Charge bins confirmed multi-bin
(H1 fix held). Gap: dominance not binned in the partition.

KEY DIAGNOSIS (2026-06-03, VERIFIED IN CODE) — the H-strata wild-negative (HN −54, HA −16,
aromatic_H −24, …) is **ANTI-PREDICTION, NOT no-data**: counts normal-to-large (HN 52, aliphatic_H
213 > N's 54), jackknife CIs tight & far-below-zero, variance-share normal. CAUSE: the fit is ONE
**GLOBAL ridge over all 846 atoms, SLICED by stratum** (`evaluate_between_tier`/`evaluate_within_tier`
train on ALL atoms; `build_score_rows` only slices) — **NOT per-type**. One coefficient set is
mis-applied to H (incompatible T2 physics: tiny proton CSA, aromatic-H self-ring excluded, H-bond/
local-dominated). Stage-1's per-element fitting was dropped when the all-atoms fit replaced it;
"condition on atom type" became SLICING not FITTING. **⇒ "H-bond doesn't lift HN" is CONFOUNDED —
not judgeable until per-type.** Fix = do BOTH (all-atoms determinability + per-type coefficients).
[[feedback_stage1_lessons]] [[feedback_no_simplification]]

BUILD 3 DONE (`d35d7ec`, 2026-06-03) — fit-architecture loop. Integrity clean (pure Python,
gates pass, anti-circular, disk-light). FINDINGS (total-T2, tier=all, per-type vs global-sliced,
between/within):
- **Per-type RESCUES the H strata** (the pooled-fit anti-prediction confirmed as OUR artifact):
  HN −54/−77 → −0.05/**+0.578**; polar_H → −0.05/+0.717; aromatic_H → −0.05/+0.499; HA → +0.05/+0.307.
- **Per-type also improves the heavy strata:** O 0.08/0.18 → 0.20/**0.72**; aromatic_heavy 0.44/0.49
  → 0.64/0.67; N 0.63/0.81 → 0.52/0.91 (N pays a small BETWEEN cost from thinness).
- **Determinability tension:** per-type BETWEEN uses 40–69 atoms/stratum (p≥atoms, thin-flagged) vs
  global's 846. So per-type-WITHIN (many frames) is the trustworthy read + where the wins are;
  per-type-BETWEEN is thin/overfit → static axis wants the hierarchical/interaction model or more
  proteins. "Both" is complementary: per-type-within = per-stratum truth; global-between = determinable static.
- **Channel SETTLED: total-T2 is the target.** total beats BOTH dia and para in every stratum
  (HN total 0.58 vs dia/para 0.25; O 0.72 vs 0.35/0.43) — total is gauge-invariant, dia/para
  gauge-dependent (noise cancels in total). dia-T2 due-diligence done; "para cleaner"/"dia channel"
  both refuted. Don't fit the split.
- **H-bond→HN REOPENED** — HN recoverable per-type (within 0.58); the earlier "H-bond doesn't lift
  HN" was the pooled-fit confound. Judge H-bond fairly in the per-type fit.
- Dominance response curves landed (14k rows, the isolation axis); C++ dom bin-id flagged for next emit.

NEXT: lead check-in (rediscover docs → restore point) → doc-truth-pass grinder (`CODEX_BRIEF_DOC_TRUTH_PASS`,
lead owns git, agent edits-only) → then between-calculator network (step 2) + equations table (step 3),
on the per-type-within + dominance-gated + total-T2 foundation.

## SUPERSEDED SESSION SNAPSHOT (2026-06-03)

The v1 handoff body was removed because Build 1/2/3 superseded it; durable discipline retained: spec -> vet -> execute -> drift/postmortem, human-in-loop, C++ spine owns physics, Python only fits emitted substrate, nmr_extract outputs sacred.

---

## SESSION HANDOFF (2026-06-02, EOD) — all-atom emit + three-tool stool

DONE this session:
- **Multi-stride `.LGS`** (`7ac1b24`): `1p9j-calibration-dense-mopac-live-orca.LGS` — dense MOPAC
  all-frames extraction + LIVE 660 ORCA DFTs (incl 910/932; the holes were a CONSOLIDATION lag, the
  jobs exist in `/shared/2026Thesis/1p9j-orcas/jobs`), original-index join formalized in the schema.
- **All-atom equivariant emit** (`2ecbe59`): the law-discovery foundation. 846 atoms × 660 DFT rows,
  all 9 ring types + all 8 BondCategory (`AllBondMidpoints`), RAW lab-frame equivariant geometry
  (disp_*/orientation_* vectors, q_over_r3, APBS/AIMNet2 payloads) — **NO imposed per-atom local frame**
  (e3nn is equivariant; frames from physics+KD-tree not the backbone test). Targets DFT raw/T2/σ_iso in
  the ORCA-aligned molecular frame. ArraySpecs carry irreps/units/sign/tensor_rank/parity/mechanism.
  Gates green (ctest 5/5, oracle PASS). The later large all-atom emit (~68 GB after MOPAC) was deleted;
  do not treat any `/tmp/rdc-all-atom-*` path as current. The build target is
  `NODE_STORE_CONTRACT_2026-06-02.md`, not a resident raw dump.
- **APBS field/EFG: NOT the radii** (A/B: real-vs-placeholder modest, `APBS_RADII_AB_WORKAROUND.md`),
  **NOT the vendor** (our wiring; Stage-1 used real prmtop radii) → back to mechanism / a MATHS-METHODS
  question. Stage-1 mining: APBS-EFG tiny, **Coulomb/MOPAC EFG moderate → the field signal is the MOPAC
  leg**, not APBS.
- **Stage-1 mutant landscape** (`STAGE1_MUTANT_LANDSCAPE.md`): "visibility changes with the question
  asked" — validates maths-vs-sample, geometric-not-IUPAC, preserve-provenance, AIMNet2-as-ceiling.
  McConnell→joint confirmed; self/bonded fix validated.
- **Run-framework** (`tools/rediscover_run.py`, `RUN_FRAMEWORK.md`): standard root + drop-old + dry-run;
  nmr_extract-guard patch firing. Cleaned 33 GB of superseded rediscover substrate. **nmr_extract
  outputs are SACRED (16h atomic) — never delete** (`feedback_static_workaround_not_producer_redo`).

REVIEWER #52 DONE (`ALL_ATOM_EMIT_REVIEW.md`): model discipline CLEAN (typed spine, no string-dispatch,
no-imposed-frame real). Current truth: #29 is **PARTIAL**, not done. `RunTraversal`/`PerRecordSink`
now owns ring/mc/charge, broad, and the no-source feature runners (`efg`,
`buckingham_efield`, `aimnet2_features`). Still outside it is the rich all-atom target/source
carrier hand walk at `AllAtomEquivariant.cpp:409`. MOPAC landed in `dca30b8` by extending the
existing all-atom path, not by adding a fourth MOPAC runner. Keep that physics wiring, but fold it
with the all-atom carrier when #29 is finished. Provenance gaps remain: atom-role/named-atom into
NPY sidecars, support/N_eff/nonzero-rank to manifest, `normalization=raw` tag. Framework guard DONE
(refuses /shared + extraction-signature; tested).

NEXT: build to **`NODE_STORE_CONTRACT_2026-06-02.md`**. The real #29 is a typed relationship index
over **all** record shapes: source-sum carriers (ring/mc/charge + broad), no-source per-target
feature carriers, and the remaining richer all-atom raw-source carrier. No new sibling drivers.
MOPAC is no longer pending; its selectors/attachers move with the all-atom fold. Then the named,
described outputs feed the **law-discovery maths model** at the Python edge
(e3nn/PySR/ridge/CV/plots/frozen `get_C()`), not a second Python protein model. larsen postponed
(grounded, `LARSEN_GROUNDING.md`). Only drop rediscover substrate; /shared extractions remain sacred.

## CORRECTED BACKBONE CAPSTONE (2026-06-02) — consolidated corrected snapshot

Full per-stratum × per-mechanism run on the FULLY-CORRECTED substrate (10 Å cutoff recorded,
valid-source McConnell, literature-scaled ring `jb_T*`, between-led statics; all corrections audited
active; gates green). Doc: `CORRECTED_BACKBONE_SNAPSHOT.md`. Charts/CSVs under
`/tmp/rediscover-corrected-backbone-snapshot-1p9j/` — **EPHEMERAL** (copy if keeping; the doc is
committed). Verdict buckets only: recovered law / form-recovered-scale-fitted / can't-make-it-work-for-now
(never "null"). Standing:
- **Ring = recovered law** — γ_bare −11.3 nA/T = Pople six-membered intensity; γ_lit ~0.77–0.92
  (compatible with 1). NOTE: the scalar aromatic-H LOAO 0.62 is NOT re-derived in this T2 snapshot
  (and not claimed); per-backbone ring rows are weak because backbone atoms aren't ring-facing.
- **Charge q/r³ = strongest backbone static T2** — N-between |T2| r 0.776 (γ 9.70±0.23), C-within 0.672
  (γ 10.22±0.29); strongest symbolic fit (PySR R² 0.775; fixed r⁻³ comparator ≤0.98).
  form-recovered-scale-fitted (per-stratum scales vary — O γ=158).
- **McConnell valid-source = can't-make-it-work-for-now ALONE** (corrected substrate confirms; high
  in-sample between-c, no out-of-sample convergent constant). Real test = the joint/ensemble fit (layer 3).
- **HN field = the one clear field-like static read** — FF14SB between R² 0.48 (A −22.4±6.4),
  Buckingham 0.40; other strata can't-work.
- **AIMNet2 = strong learnable ceiling** (not a law) — C-within R² 0.72, N-within 0.59, lifts vs
  physics N 7.7 / C 3.5; the embedding carries the local residual.
- **EFG local-frame = can't-make-it-work-for-now** (~0; lab-frame confound stays killed).
Equation fits: charge cleanest (r⁻³ recovered, PySR 0.775); bond moderate; ring weak on backbone (home
is aromatic H). Cookies to stand on: ring (law) + charge (backbone static T2). McConnell → joint fit.

## UPDATE 2026-06-02 (LATE) — ring de-circularisation CONFIRMED; emit-surface audit; McConnell Δχ next

Three landed today (Python-only or research; NO producer/library change):

**Ring de-circularisation CONFIRMED — the units "failure" was the recovered law** (`b79da2f`,
`analysis/ring_literature_decirc.py` + `RING_LITERATURE_DECIRC.md`). The old bare-kernel γ=−11.3 is
recovered as an INTENSITY: all-valid within-T0 γ_bare = −11.33 ± 3.67 nA/T, dead in the Pople
six-membered range (−11…−12.5). Literature-scaled (`jb_T*`) γ_lit compatible with 1 on BOTH axes —
within (modulation) and between (static across-atom, intercept ≈24 ppm = aromatic-H baseline):
T0 within 0.77±0.47 (r=0.61), T0 between 0.64±0.52 (r=0.72), T2 within 1.05±1.13, T2 between
0.92±0.53. CAVEATS (honest): n=1 protein; the ring signal CONCENTRATES (effective signal-N ~2.2
atoms even in the 41-atom pool — all rows thin); per-type all thin (TYR + TRP-perimeter compatible
with table; PHE modest scale-excess 1.38; HIS γ_lit 4.1 genuinely off — the protonation-dependent
imidazole, where the Pople table is weakest). Ring moves "form-recovered, scale-fitted" →
**de-circularised recovered law (six-membered), limited one-protein confidence**. The BIG QUESTION
DISSOLVES for ring: the residual is THIN-N (candidate 4), not a physics failure. The 750-DFT set
tightens effective-N directly (script reruns unchanged).

**Emit-surface audit done** (`EMIT_SURFACE_AUDIT.md`, opus). Confirms the "too much stuff" smell:
every SourceSum relationship emits the blessed kernel AND parallel geometric scalar sums for Python
reassembly (`look03_coefficient.py:86` literally says the weighted aggregate "is a C++ reducer to
add, not a Python re-sum"). Buckets A(keep blessed)/B(raw input)/C(cut reassembly-scratch: ring+mc
`sum_dipolar_*`, `field_z`/`field_mag`)/D(consolidate dupes: DFT target 3×/row, 5× writeNpyF64).
KEY: `bondKernelT2FromSources` ALREADY builds the McConnell tensor C++-side → the blessed Δχ emit is
an EXTENSION, not new code. Smell is CONVERGING not spreading. Stage cuts after replacements (ring
`sum_dipolar_*` consumed by ~10 scripts); DRY → #29.

**McConnell Δχ literature done** (`MCCONNELL_DCHI_LITERATURE.md`, codex — web verified working; the
earlier codex "failure" was over-broad-brief context overflow, NOT broken web). Unit chain (ties to
ring's −n·B sign): `σ_ppm = −0.5535·q·f` (scalar; q=Δχ in 10⁻⁶ cm³/mol, f=(3cos²θ−1)/r³); tensor
prefactor 1.6605 before the 1/3. Carbonyl Δχ DISPUTED ~×6 across sources. Recommended first set
(lead to confirm): peptide C=O +2.41 / C-N −5.42 (Williamson-Asakura, the protein lineage); sidechain
C=O +2.41 (single-family) or +6.34 (Abraham amide); **aromatic = 0 — RING already carries the π
current, an aromatic McConnell term DOUBLE-COUNTS the ring**.

**McConnell cell DONE (#36, `3d8086d`, `analysis/MCCONNELL_LITERATURE_DECIRC.md`).** Emitted the
blessed `mc_lit_T2_local_*` (per-category WA Δχ on the rediscover bond model via
`bondKernelT2FromSources`; aromatic=0 verified contributes exactly 0; T0≈0 traceless as designed;
ring+mc oracle `--case all` exit 0, ctest pass, no Python recompute). McConnell is NOT the clean ring story:
- **Form is real at backbone C** — |T2| r 0.75, comp_r −0.62, and it SURVIVES leave-atoms-out (LOAO
  R² 0.38 within / 0.31 between — the ONLY stratum with real out-of-sample T2 signal). N/CA weaker;
  HN (expected to lead — amide H sees peptide C=O/C–N) is WEAK (|T2| r 0.27, LOAO ~0).
- **Provisional WA Δχ SCALE does NOT de-circularise** — γ_lit scatters −2.1…+38 across strata, lands
  at 1 nowhere meaningful → **form-recovered, scale-fitted, NOT de-circularised** (expected: Δχ are
  provisional/web-cited, carbonyl disputed ~×6; the scale question is MOOT until the primary-refs debt
  is paid — `references/incoming/TECH_DEBT_mcconnell_dchi_primary_refs_2026-06-02.md`).
- **C puzzle / OPEN THREAD:** C's coupling is robust but NEGATIVE and ~2× literature (γ −2.14 ± 0.022,
  very tight). Likely the carbonyl C sits ON its own C=O bond → far-field bond-anisotropy formula
  weakest exactly there (near-field, like ring's in-plane breakdown). CHECK whether the McConnell emit
  excludes self/bonded contributions; if not, C's signal is probably its own carbonyl near-field, not
  a recovery. A thread, not a win.
- The table's "zero-circularity recovered law" labels (N/O between) are SE artifacts (γ 2.6 / 38±30,
  trivially "compatible with 1") — ignore; nothing de-circularises cleanly.
Ring stays the clean de-circularised law; McConnell = form present (esp C, pending the self-bond
check), scale not yet. The emit-surface C3 fix landed (blessed `mc_lit_*` added; `sum_dipolar_*` kept).

**McConnell Δχ calibration — DID NOT CONVERGE (2026-06-02, `b69134d`, `analysis/MCCONNELL_DCHI_CALIBRATION.md`).**
ORCA-free layer-2 calibration (ORCA busy with the production 750; deferred the first-principles route).
Per-category q scatters −40…+329 across strata, SEs ≈ the values, R²≈0 except C. NO coherent
DFT-calibrated Δχ — **not a cookie** (`feedback_run_the_algorithm_get_a_cookie`). A couple of strata
(N/O between) drift near Abraham-amide magnitude (q_CO≈6, q_CN≈−14) but that's a cloud, not
convergence. C is the only robust fit (q_CO −5.73±0.15, R² 0.38, |T2| r 0.76) but **sign-flipped +
wrong magnitude; cause UNKNOWN** — carbonyl-C-on-its-own-bond near-field is the leading SUSPICION
(where a far-field point-dipole must fail), NOT established. McConnell PARKED: its transferable Δχ
awaits the deferred first-principles project — three routes found 2026-06-02: ORCA magnetizability
(real toggle, `%eprnmr` "General magnetic properties: Magnetizability NO" — just off, not absent),
the NICS-fit (ghost grid + the proven dft-ex1 `! rSCAN def2-SVP NMR` format → `$SCF_Chemical_Shift`),
or LeanSCF's own `cpscf_solve` B-field response (`/shared/dft-ex1`). Ring/BS/HM are the clear visible signal.

**VET (2026-06-02, `brgjygpmm`, `MCCONNELL_PIPELINE_VET.md`) FOUND THE BUG — not the mechanism.**
Rediscover's broad/standalone McConnell source sum INCLUDES the target atom's OWN amide + near-field
bonds; the producer EXCLUDES them (`SelfSourceFilter` + `DipolarNearFieldFilter` ratio 0.5,
`../src/KernelEvaluationFilter.h`). Rediscover rejects only `r≤1e-6` (`BroadBackbone.cpp:342-345`;
comment "no self/bonded concept for bonds"). Backbone N/C/O ARE the amide atoms → dominated by their
OWN bond at ~0.5 Å (far-field point-dipole invalid); the C sign-flip = its own C=O/C–N. Ring works
precisely BECAUSE it excludes `is_self_or_bonded` (`ComposedRelationships.cpp:132-177`). CLEAN: frame,
T2 basis, units/sign, categorization (aromatic q=0). Minor: cutoff mismatch (rediscover 8 Å vs
producer 10 Å). Stage-1 confirms the mechanism is REAL but tests a DIFFERENT quantity (720-protein
WT-ALA mutation-delta static ridge vs 1-protein per-frame local-T2) — so the fix gives a FAIR test,
not a guaranteed cookie. The "not a cookie" was a BUG, Stage-1-prior vindicated again
(`feedback_stage1_prior_is_real_signal`, F1 pattern). **McConnell UN-PARKED. FIX (pending lead go):**
mirror producer source filters (self + near-field) + align cutoff + keep all-source columns for
diagnostics → re-emit 1P9J → re-run decirc/calibration with before/after C-stratum audit. C++ +
re-emit + re-run (oracle parity + ctest gated; no ORCA).

**FIX DONE + FAIR TEST (2026-06-02, `1ceab65`): McConnell does NOT converge — and the C signal was the
ARTIFACT.** Additive fix landed: per-source `mc_source_is_self_or_bonded` flag + `mc_lit_T2_local_valid_*`;
all-source columns byte-identical (163K aggregate + 13.5M source rows) → ring+mc oracle `--case all` exit 0,
ctest 4/4; standalone/oracle path untouched; consumers gained `--mc-source-mode valid|all` (default
valid). Removing C's own bonds (53,500 self-sources = 27k own-C=O + 26.5k own-C–N) COLLAPSES the C
stratum: |T2| r 0.755→0.092, R² 0.382→0.0008, q_CO −5.73±0.15 → +3.63±9.77 (noise). So the strong C
result was ENTIRELY own-bond near-field, NOT a far-field McConnell law — a would-be false finding the
bug was hiding. Fair-test per-stratum valid |T2| r is weak-modest (strongest O-between ~0.66; CA ~0.40);
per-category Δχ still scatters with huge SEs — NO coherent convergent constant. **CAN'T MAKE IT WORK FOR NOW** (this test / this data — bug fixed,
fair test given; provisional, NOT a definitive null — revisit with the 750, the MOPAC bond-order
variant, target isolation, or a subtler mistake the deep audit #40 may surface). McConnell's
Stage-1 reality stands (different quantity: 720-protein WT-ALA mutation-delta static ridge, untouched).
The valid-source filter is now in place for future McConnell + the MOPAC bond-order variant. Cookie
framing: ran it fair, no cookie. Bug-catch value: the self/bonded poison affected ALL McConnell-consuming
analyses; the C artifact would have been a false finding.

**DEEP AUDIT (2026-06-02, opus, `MCCONNELL_DEEP_AUDIT.md`, 11 candidates ranked; #40) — null NOT
earned; it points to the JOINT/ENSEMBLE fit.** Deepest (C1, HIGH): the de-circ/calibration correlate
a McConnell-ONLY predictor against the FULL DFT total T2 — McConnell is a MINORITY contributor drowned
in non-McConnell variance. Ring de-circularised STANDALONE only because ring-facing H is ring-DOMINATED;
the backbone T2 is multi-mechanism, so McConnell's fair test is the JOINT multi-mechanism design matrix
= the ensemble (arc layer 3), NOT standalone. So: "can't make McConnell-ALONE work for now"; the right
home is the joint fit. Loose ends — **RESOLVED** (recorded here + in memory so they stop resurfacing and burning tokens —
`feedback_record_resolutions_durably`): **C2** fixed in #42 — McConnell/anisotropic-bond cutoff default
is now 10 Å (= producer; `main_extract.cpp:186`, `AllAtomEquivariant.h:19`); the ring-centre 8 Å is the
INTENTIONAL aromatic-neighbourhood convention, not a truncation bug. **C5** between-axis-led report
landed in #42. **C6** the QtBond A/B sign worry is NOT a consumed-T2 sign bug — #41
(`MCCONNELL_LOOSE_ENDS_VET.md`) proved the McConnell M tensor is EVEN under A/B swap: cosθ=d̂·b̂ co-flips
with the bond axis, so term-1 cosθ·(d̂⊗b̂) is invariant and terms 2–3 are quadratic; the deep-audit's
term-1-symmetric-part claim forgot cosθ co-flips. LIVE RESIDUAL (feeds the equivariant emit, NOT a
McConnell bug): raw `bond_axis_local_*` is index-oriented (producer min/max), so as a POLAR 1o vector
its per-bond sign is arbitrary — the all-atom equivariant emit must not ingest it as a meaningful odd
vector without chemical orientation (memory `project_mcconnell_bsign_resolved`). **SDK** consumed
`mc_lit_*`/MOPAC columns → close in the #51 provenance pass. CLEAN: T2 basis, trace order/sign, DFT
target frame/basis, the self-source fix.
**METHOD INSIGHT (general):** standalone de-circularisation is fair only for a DOMINATED stratum (ring);
minority contributors must be tested in the JOINT/ensemble design matrix, not alone.

## OPEN (2026-06-02): the weak DYNAMIC field σ-response — concrete, testable; do NOT read as a real null, and do NOT excuse with "treatment mismatch"

Jessica's check: over 15 ns backbones *do* swing through varying fields (input varies: |EFG|
per-atom frame-std ≈0.20). What's weak is the within-axis Buckingham σ_iso RESPONSE (R²≈0.10).
**"Treatment mismatch" (CPCM) is NOT the explanation — it's too global/convenient, and wrong
for the within axis** (the within field-fluctuation is the protein's own charge motion, which
FF14SB/APBS capture; CPCM is a continuum reaction field ~set by geometry, not a separate
within-fluctuating source — a between/absolute footnote at most). The real, specific fact:
**the within σ IS learnable (AIMNet2 within HN 0.74) but the FIELD explains almost none of it
(0.10).** Three falsifiable candidates: (1) **units/prefactor bug** — FF14SB field missing the
Coulomb prefactor (audit); `A·E_proj + B·|E|²` mixes E and |E|² at different powers of the
missing factor → mis-conditioned → test: fix + re-fit; (2) **projection axis** — test |E| vs
E_proj; (3) **field is a genuinely weak within-driver** (within σ is local-geometry-driven,
embedding-captured). **Decisive diagnostic: the raw frame-to-frame correlation of ΔField vs Δσ,
per atom, before any fit** — present ⇒ fit/units bug (1,2), fixable; absent ⇒ field really is a
weak within-driver (3). Run that, don't hand-wave.

## OPEN — THE BIG QUESTION (2026-06-02): why doesn't a STRONG relationship de-circularise un-fitted?

Ring current is strong (form |r|=0.65; fitted universal k≈21 → held-out R²=0.62) yet the
**un-fitted literature kernel does NOT predict DFT** (literature-fixed de-circularisation:
ring-T0 γ=−11.3, McConnell-T0 γ=−4.75 — sign-flipped, magnitude-off; both bucket
"form-recovered, scale-fitted"). For a relationship this strong, "why doesn't it fit
un-fitted?" is the real physics question — but **only answerable after the little stuff is
mowed.** Candidate confounds to clear FIRST:
1. ~~units/sign~~ **RESOLVED (2026-06-02, `UNITS_AND_ISSUES_AUDIT.md`, bl7pvjdhf): a units/
   scaling + LABELING bug, NOT a physics failure.** The de-circularisation read producer BARE
   kernels (geometric unit-current BS, unit-Δχ McConnell), not literature-ppm predictions —
   so ring γ=−11.3 absorbs the omitted `LiteratureIntensity`, McConnell γ=−4.75 the Δχ×unit
   prefactor. The negative SIGN is the CORRECT ring-current convention (`G=−n·B·PPM_FACTOR`,
   negative diamagnetic intensity → shielding above ring), NOT a σ-vs-δ flip (target is ORCA σ).
   **FIX → re-run:** ring uses `jb_T*` ppm cols or BS×LiteratureIntensity; McConnell emits/applies
   a named Δχ+unit prefactor; then γ should → ≈1 = genuinely de-circularised. **So ring is LIKELY
   a clean recovered law pending the corrected test.** Plus a SYSTEMIC units-mislabeling cascade
   (BS H5 attr, `Catalog` KernelBs/Mc, `BareKernelColumns`, broad sidecar names all falsely "ppm")
   — full map + 5-item next-fix-set in `UNITS_AND_ISSUES_AUDIT.md`.
2. ~~in-plane point-dipole breakdown~~ **MOWED (2026-06-02, J-B cell `165ed08`): NOT the answer.**
   In-plane, BOTH point-dipole (T0 R²=0.02) and J-B (T0 R²=0.01, T2 comp r=0.01) are ~0 — the
   in-plane band has ~no DFT ring-modulation (self-ring near-constant CV 0.03, de-meaned away);
   the FINDINGS R²=0.67 was form-reconstruction of the producer kernel, not DFT signal. The
   misfit is NOT an in-plane-form gap. (Byproduct: the J-B/BS kernel is now emitted spine-side,
   `jb_T0`/`jb_T2_local_*`, wired through both traversals, oracle-PASS — `JOHNSON_BOVEY_REGION_RECOVERY.md`.)
3. **per-ring-type intensity** — one constant can't serve PHE/TYR/TRP/HIS if intensities differ.
4. **thin N** (~3–7 coupled aromatic H; γ SE ±3.6) — maybe unresolved, not failed.
Deepest residual once (1)–(4) are gone: **the classical kernel is a MODULATION (Δσ) model,
not absolute-σ** — a literature constant predicts the frame-to-frame modulation, not absolute
shielding, so an absolute un-fitted fit *should* miss. SEQUENCE: mow (1)–(4), THEN the residual
is the answerable question. Do NOT chase the big question with confounds live (confound-chasing).

## UPDATE 2026-06-01 (LATE) — applied-maths cleanup arc + static calibration; AIMNet2 NEXT

A full applied-maths audit + cleanup of the rediscovery analysis, then the between-axis
static calibration. Two independent audits (Opus agent + codex) on the original code →
verified fixes → capstone re-run + reusable charts → the calibration table.

**Methods error found + fixed (the big one) — the variance axis.** The uniform per-atom
de-mean answered only the WITHIN-atom (dynamic) axis and *deleted* the BETWEEN-atom
(static-environment) signal — exactly where the electrostatic mechanisms live (Jessica's
Stage-1 mutant intuition). Built `analysis/variance_decomposition.py` (between = LOAO
whole-atom, within = train-only-centred frame split, AR(1) N_eff, variance shares) —
Opus-audit-verified. The between axis uncovered real static signal that was hidden:
**charge→N T2 between |T2| r = 0.78**; Buckingham / charge-field HN σ_iso between R² ≈ 0.40–0.48.

**Two real bugs fixed:**
- **q/r² → q/r³** charge distillation (wrong radial power masked signal) → charge-N
  R²=0.981. Commit `f299a01`.
- **efg lab-frame rotation confound** — feature+target were lab-frame; de-meaning a
  *tumbling* tensor left a co-movement confound. Fixed by **local-frame EFG emit**
  (`fb90bbd`). Honest result: the efg "O 0.34" was almost entirely the confound — clean
  local-frame efg ≈ 0 across the backbone. **EFG carries ~nothing for the backbone T2**;
  the cleanup caught a would-be-false claim.

**F1 false alarm cleared (Jessica's Stage-1 prior vindicated):** the bond-anisotropy T2 is
the project's CANONICAL McConnell `K·χ` tensor (library-identical `McConnellResult.cpp:90`,
derived in `GEOMETRIC_KERNEL_CATALOGUE.md`), NOT home-rolled. Disclose the traceless-χ (PCS)
convention. Opus-1's "non-standard" finding was wrong; codex's decline was right.

**Static calibration (`d7c63e0`; `analysis/static_environment_calibration.py` +
`STATIC_ENVIRONMENT_CALIBRATION.md`):** recovered the calibration COEFFICIENT per stratum
on the between axis + within-protein jackknife uncertainty + literature comparison.
**VERDICT BUCKET (Jessica): "form-recovered, scale-fitted" for now** — kernel forms recover
(charge-N |T2| r 0.78; McConnell form), but coefficients/scales are fitted, NOT yet
de-circularised to first-principles literature values (γ≠1, e.g. McConnell-N γ≈26; the
γ-vs-units interpretation is the open thread). The clean un-fitted de-circularised law
remains ring-current Pople (k≈21, the aromatic-H result).

**Confidence scope (Jessica):** n=1 PROTEIN but many structures within it (≈500 frames,
≈50 atom-environments/stratum) → LIMITED within-protein confidence metrics apply
(jackknife-over-atoms, block-bootstrap-over-frames, autocorr-aware), scoped "this protein
across its structures"; NO population (across-proteins) inference. Jackknife SEs show which
coefficients are resolved (Buckingham HN A=−15.6±2.8) vs not (N=−100±244).

**THE THESIS REPORTING ARC (memory `project_thesis_reporting_arc` — the steering doc):**
(1) signal [done — variance decomposition] → (2) equation calibrations [done —
form-recovered, scale-fitted; all-with-something + within-protein uncertainty] →
(3) ensemble model [the good calibrated mechanisms; AIMNet2 embedding plugs in HERE IF
magic] → (4) equivariant transferability pilot [720 WT backbones, statics as comparison;
data-gated]. The fitting IS instrument calibration — the coefficient is the deliverable,
not the R².

**Detailed record on disk:** `APPLIED_MATHS_AUDIT.md`, `APPLIED_MATHS_AUDIT_codex.md`,
`FIXES_AUDIT_opus.md`, `VARIANCE_DECOMPOSITION_METHOD.md`, `STATIC_ENVIRONMENT_CALIBRATION.md`,
refreshed `BACKBONE_LAW_EVIDENCE.md` / `EFG_ARC_EVIDENCE.md`, capstone charts
(`analysis/rediscover_capstone_charts.py`, reusable on the 750). Commits `f299a01`,
`fb90bbd`, `d7c63e0`.

**NEXT (this session): AIMNet2 embedding-ceiling = thesis-arc layer 3.** Emit the 256-d
embedding (`ArrayId::Aimnet2Embedding`; it's a fail-loud stub like efg/buckingham were →
emit-then-fit) + the ceiling fit: physics-only R² vs physics+embedding gap, per atom-type —
does a learned local rep capture what the kernels can't, esp. **Cα** (the atom no
through-space/field kernel touches)? Black-box → report as the labeled **learnable
ceiling**, not a recovered law; "if magic, we take it." 1P9J 750-DFT full set lands ~3 days
→ re-fire the capstone charts + calibration + variance decomposition on it.

## UPDATE 2026-06-01 (EOD) — full arc closed; differencing PARKED; WORK CATALOG defined

The backbone exemplar demonstrates the full **B→A→de-circularised** arc (commits through
184a1ee). De-circularised: the UN-FITTED literature-coefficient kernel T2 predicts DFT for
N (component r 0.69), O (0.53), C (0.51); CA/HN/HA need the fit; ring null. Honest
multi-mechanism mixture.

**#35 differencing system-ID: PARKED.** Implemented + run (`analysis/differencing_system_id.py`;
scope `DIFFERENCING_EMIT_SCOPE.md` — NO C++ emit needed, catalog complete). The smoothness gate
FAILS at the current every-other-frame DFT sampling (Δσ ~90% noise, lag-1 ρ≈0.05); a faint
bond-led CA T2 response (r≈0.36) is consistent with the static story but sub-gate. Needs a DENSE
consecutive-frame DFT burst (NOT more every-other frames) + more contributors. NOTE: the 1P9J
750-DFT set lands in FULL in ~3 days (every-other; 5×4-core×64GB×~month) — that enriches the
STATIC analyses (re-run the exemplar on the full 750), NOT differencing's spacing. The script is
reusable when a dense consecutive-frame campaign lands.

**WORK CATALOG defined: `WORK_CATALOG.md`** — the overall work breakdown (Done / #29 engine /
item-2 calculators / #30 guard / parked #35 / the data dependency / future axes), structured for
feature-by-feature codex delegation against the complete emit-catalog + the exemplar template.

## UPDATE 2026-06-01 (latest) — BACKBONE equivariant-T2 EXEMPLAR landed (2 runs, 8 strata)

The exemplar for the rest of the backbone work. `analysis/equiv_t2_backbone_e3nn.py`
(Opus-authored, codex-built/run/debt-caught, lead-verified): heterogeneous-source
equivariant model — `o3.spherical_harmonics("2e")` of each source's direction,
per-source-TYPE radial MLP (ring/bond/charge), scatter-pooled per atom, target via the
frozen `get_C()`. Predicts per-atom DFT **T2** across the 8 strata (N/CA/C/O/HN/HA/HA2/HA3).

Emit-fix: `BroadBackboneSink.cpp` ONLY — appended `source_normal_local_*` (rings) +
`bond_axis_local_*` (bonds) per source (data was already in `SourceSlot`; mirrors
RecordSink; charge → zero sentinel). Ring/mc oracle intact by construction (no other
C++ touched). The first `--with-axes` run on `/tmp/rdc-broad-backbone-axes` is now
marked TAINTED because it consumed ring normals and index-oriented bond axes through a
polar `1o` cross path.

Parity re-check (`eb2bf03`, same substrate, corrected consumer): `--with-axes` now uses
only parity/sign-safe `l=2` axis terms (`Y2(axis)`, equivalent to traceless axis⊗axis);
`disp_local_*` remains the only polar `1o` input. Corrected frame-split T2 R²:
  disp-only → corrected axes:  N .573→.630  CA .330→.418  C .563→.574  O .587→.666
                               HN .686→.712  HA .598→.648  HA2 .792→.823(thin,4at)  HA3 .884→.911(thin,4at)
Verdict: the axes lift **SHRANK** versus the tainted cross-path run, but did not vanish;
it remains positive in every stratum. Solid on the 54/52/50-atom strata; HA2/HA3 thin
(4 coupled atoms — correlate-not-match, flagged).

Discipline: no recompute (lead re-grep clean), frozen-C reused, e3nn (no hand-roll),
C++ untouched except the broad sink. Axis columns are not consumed as polar `1o`. Debt caught
(cached Y2 features, LOAO opt-in, dead-code removed). The Python-consumer discipline is now a
controlling doc: `analysis/PATTERNS.md` (front-load in every fitter brief). NEXT backbone
work copies this exemplar.

## UPDATE 2026-06-01 — Python physics RETIRED; equiv-T2 rebuilt on e3nn (authored)

Per `MODEL_PLACEMENT_PROPOSAL.md` + the lead's decisions. The Python end-runs are
gone; the equivariant fitter is e3nn on the C++ export. The C++ substrate was NOT
touched (no recompile; the carrier `compact_npy` change is a separate codex task).

- **Deleted** `analysis/equiv_t2.py` (its numpy `lib_T2` was a byte-for-byte
  clone of `DecomposeLibrary` — the projection end-run).
- **Rebuilt** as `analysis/equiv_t2_e3nn.py`: `o3.spherical_harmonics("2e")` for
  Y2(r̂)/Y2(n̂) + a `1o⊗1o→2e` cross term (BOTH `--cross exact` fixed-path and
  `--cross learnable` FullyConnectedTensorProduct; `--cross both` reports both and
  picks the better by frame-split R² — decision 4) + invariant radial MLP +
  scatter-pool. Consumes the C++-emitted per-source `disp_local`/
  `source_normal_local` + the REQUIRED emitted target NPY
  `rediscover_ring_current_sources_target_local_T2.npy` (fail-loud if absent; NO
  numpy projection fallback). Target mapped library-basis→e3nn-2e by a PINNED
  constant.
- **`analysis/change_of_basis.py` + `analysis/test_change_of_basis.py`**: the 5×5
  library-T2↔e3nn-2e change-of-basis is the ONLY surviving "recompute," written as
  a fixture-pinned TEST — orthogonality/round-trip + a Wigner-D equivariance
  round-trip vs the C++ library tensor (handles e3nn's y-z-x convention). The
  derivation uses e3nn's own `ReducedTensorProducts("ij=ji", i="1o")` so features
  and target share the e3nn-2e basis (verified: e3nn 1o is the plain 3-vector,
  no hidden permutation — `_reduce.py:126-129`).
- **Pople-comparison recompute arrays DELETED** (lead decision 3, no kept "labeled
  integrity test"): `sumpool_kernel.py`, `refine_kernel.py`, `pysr_distill.py`,
  `sumpool_t0.py`, `sumpool_mcconnell.py` keep their pooling FITS but no longer
  build `(3cos²−1)/r³`; comparisons read the C++-emitted `dipolar`/`bare_T0`, and
  PySR carries the symbolic "form fell out" claim. `look03_coefficient.py` reads
  the emitted `sum_dipolar_producer_valid` aggregate (no `Σ intensity·dipolar`
  re-sum; intensity-weighted aggregate is a future C++ reducer).
  `look_charge_dipole.py` lost the `Σ q·d/r³` field recompute — fits only the
  emitted `mu_*` (the field is the C++ `buckingham_efield`/APBS relationship, a
  fail-loud stub on this branch). `look01` repointed to the emitted aggregate col.
- **GREP PROOF**: no `(3cos²−1)/r³`, `/r**3`, `q·d/r³`, or `lib_T2` projection
  arithmetic remains in any analysis `.py` outside the pinned change-of-basis
  test. (Run from `analysis/`: see the grep block in the handoff report.)
- **Env** (decision 1): system `/usr/bin/python3` has torch 2.11.0+cu130 + e3nn
  0.6.0. PINNED in `analysis/requirements-e3nn.txt`; `LD_LIBRARY_PATH` gotcha +
  run commands in `analysis/ENV.md`. rediscover NPYs NOT added to the SDK catalog
  (decision 2 — never-merge spike).
- **RUN + GATE: PASSED (lead, 2026-06-01).** Ran in system python (torch
  2.11+cu130 / e3nn 0.6.0, `LD_LIBRARY_PATH` per `ENV.md`) on `/tmp/rdc-composed`.
  Change-of-basis: all 3 checks pass — fixed the derivation to float64 (e3nn
  defaults float32, which under-precised C to ~1e-7 and tripped the test's strict
  thresholds; the matrix was always exactly orthogonal). `_C_FROZEN` frozen as the
  exact 5×5 (0,±1,±0.5,±√3/2); `get_C()` loads it with NO e3nn (orthogonal to
  1.1e-16). e3nn fit `--cross both`: **frame-split T2 R²=0.466** (baseline 0.467),
  **|T2| r=0.757** (baseline 0.756) — reproduces the retired hand-rolled result to
  noise, confirming the equivariant signal is REAL, not a hand-rolling artifact.
  `exact` cross-term beat `learnable` (angular structure fixed, not fitted);
  leave-atoms-out 0.370 reported, not gated (thin ~7 coupled aromatic H). Offense
  grep-clean (no physics recompute outside the labeled change-of-basis test).
  Committed on `h5-reader-pysr-spike` (never merged).

## UPDATE 2026-06-01 (later) — broad-backbone BUILT + GATED (commit 35f3768)

Big-one #1 done: `broad_backbone` composed THROUGH the engine (Claude-authored,
codex-built+gated, lead-verified). All 6 backbone frame classes resolve **0%
invalid** (N/CA/C/O=54, HN=52, HA=58); new N-frame convention z=unit(CA−N), x in
the peptide plane. Heterogeneous selectors [rings/bonds/charge-field via the
GENERAL KD backends] — breadth from the selector list, not a procedural walk.
Carrier target-repeat fix held (27-col source rows, no `dft_*` columns). Tests
12/12 (√6 + frames + rotation-equivariance + reducer-sum). Ring/mc byte-parity
oracle STILL exit-0 (lead independently re-ran — NO regression); commit touched
only the new broad files + additive `LocalFrameBasis` (existing builders
untouched) + `main_extract` + tests + `_catalog.py`.

Per-atom-type σ_iso R² (correlate-not-match, first-pass, features
[ring_sum_dipolar, bond_sum_dipolar, field_z, field_mag]): HN 0.45, N 0.45, C
0.31, O 0.20, HA 0.24, **CA 0.055** (weakest — CA is local-bonding-dominated, so
through-space mechanisms barely touch it; diagnostic about kernel completeness,
NOT a physical conclusion). Cutoff sweep 6/10/12 Å **flat** (field
short-range-saturated, not truncation-starved — revises the charge_dipole
"sweep it" hypothesis).

SCALE FINDING (→ #29 input): per-source CSV is 1.7–5.5 GB across the sweep (charge
= 8.1M rows at 6 Å). The target-repeat fix prevented worse, but per-source CSV
doesn't scale for the charge mechanism → NPY-for-source-rows (or don't emit
per-source charge rows) belongs in the totality design.

ENGINE FINDING (→ #29): broad needed a SIBLING runner (`RunBroadBackbone`) because
`RunRelationship` hardcodes the ring/mc scalar-sum sink. DECIDED (Jessica): fix
the engine, but as a REVIEWED TOTALITY design — the fold/sink/carrier across ALL
9 relationship shapes, minimal+clarifying abstraction, NOT whack-a-mole. #29,
blocked until the design + review eyes land. Sibling runner kept meanwhile;
oracle untouched.

## UPDATE 2026-06-01 — functional API BUILT + byte-parity-VERIFIED (commit 99cdc85)

BUILT + GATED (codex build, lead independent re-verify). Compiles + ctest pass —
one real fix: Qt's `slots` macro collided with `verbs::slots` → renamed
`verbs::ringSlots` in the NEW files only. The composed engine reproduces the
procedural oracle **byte-for-byte** — all 4 CSVs (141000/20500/812205/26000 rows)
+ 12 sidecar NPYs identical (independently re-run by the lead into a fresh dir,
GATE: PASS). Physics held by identity: ring k=21.1 / coupled R²=0.616, equiv-T2
basis 4.88e-8 / R²=0.467 / |T2| 0.756, mc R²=0.549 / readout 0.919; DFT frame rot
~9e-5°. Commit 99cdc85 touched ONLY the new API files + `main_extract` +
`CMakeLists` (cells/spine/RecordSink/library untouched — frozen-oracle rule held).
The functional API is **now the proven default path**; broad-backbone (#26) is
unblocked. Compile-knob answer: a Claude Agent-tool subagent STILL can't reach the
compiler here (even `which cmake` denied) — author-with-Claude / build-with-codex
is the working split.

The original authoring note follows (now superseded by the BUILT status above):

### (superseded) functional API authored, UNBUILT

Authored the composable functional API that SURFACE_DESIGN specifies, replacing
the three monolithic procedural walks as the DEFAULT path (cells kept as the
reference oracle). New files: `Verbs.{h,cpp}` (Layer-1 primitive verbs, thin
over the spine), `Relationship.{h,cpp}` (Layer-2 curried-closure combinators +
the named `Relationship` bundle + `atomsWhere`/`slotsBackend`/`nearBackend`
builders), `RelationshipEngine.{h,cpp}` (Layer-3 one pure iterated-closure
loop), `ComposedRelationships.{h,cpp}` (ring_current + mcconnell rebuilt as
composed `Relationship`s). `main_extract.cpp` gained `--engine
{composed|procedural}` (default composed; procedural runs the oracle cells for
the parity diff). `analysis/oracle_parity.py` is the gate runner (diffs
composed vs procedural CSVs+NPYs byte-for-byte).

**The compiler was NOT reachable** from the authoring agent's sandbox (Bash
denied `cmake --build`, even `which cmake`) — answers the open knob question:
agent Bash here cannot compile. So the API is authored + rigorously
self-reviewed, NOT built. **codex takes the compile + the oracle gate** — full
inventory, self-review, and the exact build/test/gate commands are in
`HANDOFF_TO_CODEX.md`. Topology reused not regenerated; reader owns H5; GUI
untouched; library not linked; never merged.

## NEXT (read `BROAD_BACKBONE_NEXT.md`) — do the BROAD case, not more narrow cells

Decided 2026-05-31: STOP building narrow single-stratum × single-mechanism cells
("single-bond-thingies"). The next move is **every backbone atom, all mechanisms
in one heterogeneous bundle** — which (a) is the thesis target (backbone shifts),
(b) stress-tests whether the architecture composes the GENERAL analysis or has
baked narrowness, (c) forces the unbuilt backbone frames (Cα/CO/HA/N). Full plan +
the charge_dipole carry-overs (field-not-μ, cutoff sweep, AIMNet2/APBS, carrier
target-repeat fix) are in `BROAD_BACKBONE_NEXT.md`. charge_dipole cell committed
`3103e73` (μ null, field carries: Buckingham field_z r=0.46, LOAO R²=0.21).

## UPDATE 2026-05-31 (late) — multi-scenario surface BUILT, faithful-rebuild gate PASSED

The general surface is implemented and validated (built by **codex** — codex has
full unsandboxed agency via `codex exec --dangerously-bypass-approvals-and-sandbox`
from the lead; Claude Agent-tool subagents are sandbox-blocked from compile, a
wrong-knob misconfig, not a hard limit — see `reference_subagent_build_agency`).
Committed `1bd61aa` on `h5-reader-pysr-spike` (NEVER MERGED).

- **Spine built** (`src/rediscover/`): `Catalog`, `TemporalIndex`, `TypedAtomIndex`
  (scoped `select`/`selectUnique`, IUPAC-trap-safe — the positional `front()` anchor
  fallback is REMOVED), `SpatialIndexSet` (per-cloud KD trees, near + range/annulus),
  `RingGeometryCache`, `ChargeStore` + FF14SB read from `topol.top` inline `[atoms]`
  (typed resnr/order cross-check, no glob), `ResidentIndexes`, `OutputManifest`,
  `AnalysisBody`. Per-relationship schema + `relationship_kind` + T2 sidecar NPYs
  documented in `python/nmr_extract/_catalog.py`. No `PbcCellSeries` (PBC=None).
- **ring_current + mcconnell ported** to the Body/catalog/index surface and
  **reproduce the ORACLE from a fresh rebuild** (not the one-off output):
  ring k=21.1, coupled within-atom R²=0.616, equivariant T2 R²=0.468, |T2| r=0.757,
  basis 4.88e-8, frame rot ~1e-4°; mcconnell scalar R²=0.55, kernel readout
  r=0.918/R²=0.843. `ctest h5reader_rediscover_tests` passes; GUI untouched.
- **The 7 others fail loud** (exit 2, ValidateScenario): buckingham_efield, efg,
  charge_dipole, charge_quadrupole, larsen_hbond, charge_response_gradient,
  aimnet2_embedding. MOPAC charge source = absent→loud (per-frame data lands AM);
  AIMNet2 recognized but multipole reducers intentionally unrunnable.

**Next (build, when data/decisions land):** wire the 7 fail-loud stubs — charge
multipoles once charges flow (FF14SB done, MOPAC AM, AIMNet2 charge), the
per-atom-feature items (efg/efield/CRG/embedding) against the carrier, larsen
once its detection/classifier decisions are made. Equivariant-T2 path is proven on
ring; extend to the new items per the frame resolution. The one emergent issue
codex hit + handled: topol↔model atom-name aliasing during charge load (typed
residue/order cross-checks, not positional).

## UPDATE 2026-05-31 (evening) — canonical re-run done + scalar fit landed

The value-affecting Codex fixes are IN, built clean (lead session), re-run on
1P9J, and verified. Canonical output: `/tmp/rediscover-out-v2/` (same shape:
141000/20500 ring, 812205/26000 mc). Fixes that landed:
- **Ring identity + self/bonded exclusion** — `ring_index` + `is_self_or_bonded`
  per source; aggregated split `sum_dipolar_all` vs `sum_dipolar_producer_valid`
  + `n_sources_valid`. Self/fused detected by identity (own-ring set + shared-atom
  overlap), NOT a distance proxy; self-rings verified at r≈2.49 Å, 21.6% of rows.
- **Aromatic-H frame anchor → typed CG/CD2** (the topology-convention time bomb):
  `typedRingAnchor` picks the unique γ-carbon (or δ-carbon CD2 for TRP-benzene) by
  typed `Locant`, no name strings, no positional first-atom. `frame_anchor_atom_index`
  emitted, 100% populated (13 distinct anchors).
- **Ring normal flipped to canonical traversal** `(v1−v0)×(v2−v0)` (local to the
  rediscover code; `FitRingGeometry`/GUI untouched).
- **McConnell**: `bond_axis_local` (unit axis in local frame, verified |·|=1) +
  `bond_atom_a/b` endpoints; `cutoff_A` recorded (CLI `--mc-cutoff`, default 8.0,
  the conventions' aromatic value — producer's exact MC cutoff still unverified).
- nanoflann `[[nodiscard]]` captured.

**Scalar fit ran and is credible (see `analysis/FINDINGS.md`).** Pipeline:
sum-pooling NN learns the per-source function, PySR distils it. Recovering the
Pople form from the producer kernel `bare_T0` is partly CIRCULAR (reverse-engineers
the producer's own Giessner-Prettre formula) — demoted from the headline. The
NON-circular, instantaneous result (ring current is a static map; the trajectory
is a geometry sampler, not a process): a universal coefficient k≈21 ppm·Å³ predicts
the within-atom shielding of HELD-OUT ATOMS (leave-atoms-out, autocorrelation-free)
at R²=0.62 on the ~7 coupled atoms, against independent DFT, on the identity-clean
`sum_dipolar_producer_valid`. Thin (one protein, ~7 coupled aromatic H) but real.

**T2 Cartesian-frame caveat — RESOLVED (2026-05-31 eve), the disciplined way.**
The check lives in the READER, not a Python h5py hack (the reader owns H5): the
ORCA parser now additively retains the `CARTESIAN COORDINATES (ANGSTROEM)` block
(`DftAtomShielding::orca_coord`), and `ExtractionSupport::CheckDftFrameAlignment`
Kabsch-compares the ORCA-tensor frame to the H5 positions the extractor already
holds. Verdict on 1P9J (500 frames × 846 atoms): rotation mean 8.9e-5°, max
2.4e-4°; RMSD 0.0005 Å. The ORCA tensors ARE in the H5 frame — **T2 components
are valid as emitted, no rotation needed** (the sub-millidegree residual is
float-print precision). The equivariant T2 model (task #23) is UNBLOCKED.

**Still open:** (b) DFT
validation (element-equality, raw total≈dia+para) — DEFERRED (defensive, not
fit-affecting). (c) the literature-coefficient-FIXED test (un-fitted constant →
DFT) — the final de-circularising check, not yet run. (d) McConnell producer-kernel
reconstruction gap (R²≈0.55 — likely fuller anisotropy than one bond-axis angle).

## The extractor works end-to-end (verified in the lead session)

`h5reader-rediscover` (branch `h5-reader-pysr-spike`) **compiles, links,
and runs on 1P9J**, and the basis fixture passes (`√6`).

- Build: `cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests`
- Test: `ctest --test-dir build/linux-gcc -R h5reader_rediscover_tests`
- Run: `build/linux-gcc/h5reader_extract --run /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft --out <outdir> --case all`
- Verified output: 4 CSVs — `ring_current_{sources,aggregated}.csv`
  (141000 / 20500 rows; 41 aromatic ring-facing H × 500 DFT frames),
  `mcconnell_{sources,aggregated}.csv` (812205 / 26000; 52 backbone HN ×
  500). 846 atoms, 751 frames, 500 DFT rows, 0 gaps. Two row kinds per
  case (un-summed per-source + aggregated summed-feature) sharing the
  identity + DFT-target columns.

**Verified:** builds, runs, well-formed tables, basis `√6`.
**NOT yet verified:** the physics correctness of the *values* (Codex is
checking — see below).

## Fix already applied

`RingCurrentNeighborhood.cpp` + `McConnellNeighborhood.cpp` had the wrong
include `"../model/QtTrajectoryH5.h"` → corrected to `"../io/QtTrajectoryH5.h"`
(that header lives in `src/io/`). This was the only compile error.

## Two confirmed must-fix caveats (Jessica confirmed both real)

1. **Aromatic-H local-frame anchor — topology-convention violation, the
   project's classic time bomb.** The frame x-axis is anchored on the ring's
   FIRST canonical-walk atom (positional) instead of the chemistry-typed
   anchor the conventions require: CG for PHE/TYR/TRP-pyrrole/HID/HIE/HIP,
   CD2 for TRP-benzene (per `spec/substrate_conventions_2026-05-30.md`). This
   is identity-from-position — the banned pattern (`feedback_identity_from_chemistry_not_position`,
   the IUPAC-revert episode). FIX: anchor x on centroid→(typed atom);
   **re-read the exact per-residue anchors, do not guess.** Invariants
   (z = ring normal, r, cosθ = z/r) are unaffected; only the azimuthal x/y
   frame is wrong. File: `src/rediscover/LocalFrameBasis.{h,cpp}`.
2. **T2 Cartesian frame.** ORCA DFT tensors are in the ORCA-input frame; H5
   kernels in the MD frame. If those orientations differ, T2 *components*
   don't correspond — comparison is meaningless. T0 (iso) and |T2|
   (magnitude) are rotation-invariant and safe. RESOLVE empirically:
   compare ORCA-input Cartesian coords (from a `.out`) to H5 positions at
   that frame. Same → T2 stands; rotated → rotate the DFT tensor into the
   H5 frame via the shared atoms, or restrict T2 to invariants + flag
   components unverified. **Do not assume.**

## In flight: Codex critique

A Codex critique of the WORKING code + output is running (background,
started ~16:30). Output log: `/tmp/claude-1000/-shared-2026Thesis-nmr-shielding/37ab6452-233d-4c70-9077-8b3440b92984/tasks/blu78kz6j.output`
(huge log — grep near the end for the verdict). It assesses the physics of
the values, both caveats above (incl. running the frame-match check),
match-to-design, the additive edits, and the nanoflann `radiusSearch`
`[[nodiscard]]` warning. Read it, fold it in, then apply the two fixes +
whatever else it surfaces.

## Agent build-agency (flagged twice)

Subagents' Bash is sandboxed in this environment and denied the compiler/
CMake even with their override; only the lead session's
`dangerouslyDisableSandbox` works. So agents could WRITE but not BUILD —
which is why the include bug shipped blind. Agents have had build agency
before; this is a sandbox/settings restriction to RESTORE so the dozen
incoming agents can compile what they write. NOT yet fixed. Interim:
lead/Codex builds; agents write + review.

## Next actions (rough order)

1. Read Codex's report (`blu78kz6j`); apply its fixes.
2. Fix caveat 1 (anchor → typed CG/CD2, re-reading the convention) — committed.
3. Resolve caveat 2 (frame check → rotate or flag).
4. Restore agent build-agency (settings) for the incoming agents.
5. Run a fitter on the substrate (task #22) ONLY after correctness is
   settled. The fitter is OPEN — ridge / scalar SR / equivariant SR /
   equivariant sum-pooling; NOT decided, NOT PySR-only, GNN/equivariant
   not foreclosed.

## Durable discipline (unchanged)

Substrate is method-agnostic; fitter open. Truthful docs, no hyperbole.
Reuse the reader; additive edits only; GUI untouched. Experimental
one-shot branch (no integration). Guardrail memory:
`feedback_build_inmemory_export_dont_relitigate`.

## Codex critique findings (landed, blu78kz6j)

**Confirmed correct:** builds + test pass; scalar physics right — ring
`cosθ=z/r` + `(3cos²θ−1)/r³`; McConnell bond-axis `(3cos²θ−1)/r³` (not /3r³);
`σ_iso=trace/3`; library basis ordering + `√6` fixture; intensity from
`LiteratureIntensity`; per-source rows sum exactly to the aggregate. No
`inf`; `nan` only in 30 ring rows (`ring_in_plane_angle`, near-axis, expected).

**Real problems to fix (source semantics — the scalar math itself is sound):**
1. **Self/bonded ring included** (biggest new physics finding): ring-current
   uses every H5 ring slot incl. the H's OWN parent ring; the producer
   excludes self/bonded. Fix: emit `ring_index` + an `is_self_or_bonded`
   flag; split `sum_dipolar_all` vs `sum_dipolar_producer_valid`.
2. **Ring source identity missing** — add `ring_index` to `SourceSlot` + CSV
   (McConnell emits `bond_index`; ring doesn't).
3. **Aromatic-frame normal sign not canonicalized** — SVD normal not flipped
   to the ring-traversal convention (library `Ring.cpp:32` does). Scalar OK;
   `disp_local` / local-tensor / azimuthal components can randomly flip.
   Fix: flip normal per `(v1−v0)×(v2−v0)`.
4. **Aromatic-H anchor (caveat 1, confirmed)** — `ring.atomIndices.front()`
   not the typed CG/CD2; works today only by coincidence. Fix: typed anchor
   + emit `frame_anchor_atom_index`.
5. **McConnell source rows lack the bond-axis vector in the local frame** —
   an equivariant/tensor fit can't reconstruct the McConnell tensor from
   scalar+midpoint. Fix: add `bond_axis_local_{x,y,z}` + endpoint indices.
6. **McConnell cutoff hard-coded 8.0 Å** — producer uses a configurable
   cutoff (commonly 10 Å); design says don't hide cutoffs. Fix: record/require
   it; use 10 Å when comparing to the producer.
7. **T2 frame unverified (caveat 2, confirmed)** — fix: mark T2 columns
   frame-unverified, OR Kabsch-verify ORCA-input vs H5 frame and rotate the
   tensors back if ORCA reoriented.
8. **DFT validation weaker than spec** — `DftShieldingLoader` doesn't check
   parsed-element == protein-atom element, nor raw 3×3 `total≈dia+para` (only
   iso). Fix: add both.
9. **nanoflann nodiscard** (`FrameSpatialIndex.cpp:59`) — cosmetic (`hits`
   filled by reference); tidy anyway.

**DECIDED — C–H only for now** (Jessica, 2026-05-31). Honest caveat she
raised: the accumulated narrowing choices — single protein, DFT subset,
C–H only, cutoffs — may starve the fit of statistical depth; flagged, and
we follow it to the end regardless (report effective N alongside any fit).
`IsAromaticRingHydrogen` (HA/H4/H5) = aromatic
**C–H only**; it excludes N-bound aromatic/exchangeable H (TRP `HE1`,
protonated HIS N–H). Right if the stratum is aromatic C–H; if it should be
"all aromatic-ring-attached H," change to parent-heavy-atom-in-aromatic-ring.

Codex bottom line: core formulas, basis, tensor retention, row shape, and
sampled magnitudes are sound; the real risks are source semantics
(self/bonded rings, unverified tensor frame, missing orientation/provenance
columns, normal/anchor conventions).
