# Thesis discussion — Part 1: MSc framing and the set of wins

**Date:** 2026-04-30
**Status:** Discussion record, captured straight. Part 1 of a multi-part conversation record. Not a plan-of-record. Not a substitute for `md-rerun-685-discussion-priors-2026-04-30.md` which captures the MD/MDP priors separately.

---

## Section 1 — What this document is

This is a record of the MSc framing the discussion landed on at the end of 2026-04-30, and the structure of work that follows from it, captured without ranking, prescription, or "first move" filtering. The next conversational moves should pick up from these items, not from a curated subset.

When the assistant proposed candidates A-F (narrow niche choices) in the prior turn, it was scoping wrong. The user's reframing pulled the conversation back. This document captures the post-reframing picture as raw as practical.

---

## Section 2 — The MSc framing (user's words)

> "this is a 1 year MSc program and we do not need to make large claims. We need to show a competent, small, early stage research project with a good lit review which addresses physics and biology in the context of ML for NMR."

> "The dept head told me EVERY MSc student in this program leaves with a long list of what their investigation would have done next. The fact we found and incorporated these things is gold. What makes us keep the gold is showing some effect."

> "If we can find some small thing the latest papers left out to touch on, nifty. If we root it just in anisotropy and a way to use it in this, that is fine too."

> "This discussion is the tidal zone where we can creatively look at the literature and ask what we can meet, before we do our final big MD run and our final deep pass on the calculators."

> "The MD time series is a must per my advisor."

Implications captured straight:

- A competent, small, early-stage MSc research project is the deliverable shape.
- Large claims are not needed and should not be forced.
- "Calculators that incorporate real NMR findings are gold" — for methodology, for discussion, for the lit review chapter.
- The fact of having found and incorporated literature findings into the calculators is itself part of what makes the work valuable.
- "Showing some effect" of those incorporated findings is what keeps the value defensible.
- Leaving with a long list of what investigation would have done next is normal and expected for this program.
- The MD time series is a must — non-negotiable per advisor.
- The GNN model on top is the stated intent ("at least a GNN model on top of that").
- The lit review is part of the work; addresses physics + biology + ML for NMR.

---

## Section 3 — What we already have to lean on

Captured as inventory, not curated by current usefulness:

- **OpenFold3 local + embeddings on the DGX Spark.**
- **RefDB and BMRB** as experimental anchors.
- **Three 5090 machines with heft** for MD + ML work.
- **The T2 (rank-2 tensor) extractor** and the maths and ideas behind it. The thesis-load-bearing claim that other ML-for-CS predictors don't share.
- **The ability to quickly build new calculators.** The kernel architecture is in place; adding a new physics-decomposed kernel is bounded work.
- **Stage 1 mutant calibration corpus.** 720 proteins, 446K atoms, 55 kernels, R² = 0.818 with per-element / per-atom-type stratification. The "calculators see real things" foundation.
- **AMBER charge slice.** Six steps green in working tree (62/62 acceptance tests passing); awaiting commit. The AMBER cleanup the conversation references.
- **IUPAC topology layer.** AtomTopology + IupacAtomName + AtomReference shipped 2026-04-26 (commits fc76f47 / 145d2cc / 6448272). With outstanding debts: pro-R/pro-S verification, variant overrides, DFT-compare completeness (MutationDeltaResult).
- **260-DFT calibration set.** 10 proteins × 26 frames at 1 ns intervals from the 2026-04-18 batch.
- **685-protein MD coverage fleet** (existing trajectories that the re-run will replace).
- **2-protein dense validation set** (1Z9B + 1P9J in test fixtures; per memory entry `project_dense_validation_2proteins_20260424.md`).
- **PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md** running translation file (paper → candidate calculator).

---

## Section 4 — Structure of the set of wins (user's articulation)

> "Clean up the amber stuff as we are doing now, and walk through the calculators and make them better. Use the IUPAC naming output to feed better stats and a chapter from the mutant work. That is one set of wins before we get to the MD (or rather, while we work through the MD and the models)."

Captured as the user described it:

```
WORKING NOW              AMBER charge slice cleanup (in progress, 6 steps green,
                         awaiting commit)

WORK PASS                Walk through the calculators and make them better.
                         Incorporate the literature findings the cache + 2026
                         survey surfaced. Show some effect.

CHAPTER                  Use the IUPAC naming output to feed better stats and
                         a chapter from the mutant work.

PARALLEL / LATER         The MD re-run.
                         The model (GNN on top).

LEAVE WITH               Long list of what the investigation would have done
                         next. (Per dept-head expectation.)
```

The key sequencing claim: **AMBER cleanup + calculator pass + IUPAC chapter is one set of wins that lands before / while the MD re-run and the model work proceed.** Not a strict ordering; the calculator pass and the IUPAC chapter can run while the MD planning + execution is in motion.

---

## Section 5 — Calculator pass items (literature findings → kernel implications)

Flat list. Each item is: the kernel(s) it touches, the literature finding, the kernel implication. No ranking, no "first move" recommendation, no cost ordering. Pick from this list during the conversation that follows.

These items come from the Track C local-cache scan (full report on disk in conversation tool-results) and Track B Feb-Apr 2026 scan, both run 2026-04-30.

### 5.1 Ring-normal stability (Sahakyan-Vendruscolo 2013, JPCB 117 1989-1998)

- **Finding:** Three-atom ring-normal definitions are unstable to MD out-of-plane fluctuations. Their fix: compute two independent ring normals (e.g., atoms 1-2-3 and atoms 4-5-6 of a six-membered ring) and average.
- **Kernels affected:** `BiotSavartResult`, `HaighMallionResult` — both compute ring-current contributions and need ring normals.
- **Possible kernel implication:** Replace single 3-atom normal with two-normal averaging across all aromatic ring kernels. Expose Δθ between the two normals as a per-frame diagnostic.
- **Cache filename:** `references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-{1-4}.txt`.

### 5.2 Single-loop vs two-loop ring current (Agarwal et al. 1977, Can J Chem 55 2575-2581)

- **Finding:** [10]-paracyclophane benchmark: single-loop S=0 fits experiment with r=0.9928 versus two-loop S=1.28 Å with r=0.8883.
- **Kernels affected:** `BiotSavartResult` (loop separation parameter).
- **Possible kernel implication:** Audit current loop separation default. Document with citation. Re-derive constants against the Agarwal benchmark explicitly.
- **Cache filename:** `references-text/agarwal-1977-ring-currents-local-anisotropy-paracyclophane-text-{1-3}.txt`.

### 5.3 H-bond angle θ over distance d (Yi-McDermott 2024, J Phys Chem Lett 15 2270-2278)

- **Finding:** AF-QM/MM on DHFR-trimethoprim, 1000 ns. Hydrogen-bond ANGLE θ correlates strongly with ψ and shifts; H-bond DISTANCE only weakly. Backbone amide ¹⁵N excursions track angle, not distance.
- **Kernels affected:** Any existing kernel that uses H-bond geometry. (Verify: which?)
- **Possible kernel implication:** Add per-residue H-bond angle as an input feature. If no current kernel uses H-bond geometry, this is a candidate new kernel.
- **Cache filename:** `references-text/yi-mcdermott-2024-temperature-shifts-conformational-dynamics-text-{1-3}.txt`.

### 5.4 Buckingham σ^(1)·E coefficient (Buckingham 1960, Can J Chem 38 300-307)

- **Finding:** Theoretical foundation of σ^(1)·E linear electric-field shielding. Δσ ≈ −2×10⁻¹² E_z − 10⁻¹⁸ E². At E ≈ 7 Å from a unit charge, the linear term is ~0.2 ppm — about 20× the quadratic term. Specific coefficients: −3.0 ppm/au for nucleic acids per Case 1995, −3.4 ppm/au for proteins.
- **Kernels affected:** `ApbsField` and any Coulomb-from-water-positions kernel.
- **Possible kernel implication:** Sanity-check our existing EF kernel coefficients against these analytical values at known geometries. Document the σ^(1) coefficient explicitly in kernel comments. If kernel under- or over-shoots the expected magnitude, that's a finding and a fix.
- **Cache filename:** `references-text/buckingham-1960-chemical-shifts-polar-groups-text-{1-3}.txt`.

### 5.5 Per-element CSA expectations as kernel validation (Boyd-Skrynnikov 2002, JACS 124 1832-1833 + Sahakyan-Vendruscolo 2013)

- **Findings:**
  - Boyd-Skrynnikov 2002: Closed-form σ_xz formula extending Johnson-Bovey σ_zz to the full rank-2 ring-current tensor. At G42 of fibronectin type-2 module, ring-current-dominant N-H⋯π hydrogen bond gives σ_rc CSA contribution of 16.6 ppm.
  - Sahakyan-Vendruscolo 2013: For ¹⁷O, EF sensitivity reaches 80 ppm.
- **Kernels affected:** All T2-emitting kernels (`BiotSavartResult`, `HaighMallionResult`, EF kernel).
- **Possible kernel implication:** Add geometric test cases to the unit test suite that assert tensor magnitude reaches these values at appropriate geometries. Catches T2 normalization or sign bugs that would otherwise ship into the re-run.
- **Cache filenames:** `references-text/boyd-skrynnikov-2002-ring-current-chemical-shielding-anisotropy-text-1.txt`; `references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-{1-4}.txt`.

### 5.6 EF dominance for ¹³C / ¹⁵N / ¹⁷O (Sahakyan-Vendruscolo 2013)

- **Finding:** Direct analysis of Pople, Waugh-Fessenden-Johnson-Bovey, Haigh-Mallion ring-current contributions and electric-field effects on RNA base shifts. For ¹³C, EF dominates over RC (R = 0.702 vs 0.257). For ¹H, RC+EF reaches R = 0.961 (rms 0.079 ppm). For ¹⁷O, EF sensitivity up to 80 ppm.
- **Kernels affected:** EF kernel weighting; the relative contribution of RC vs EF kernels in the calibration ridge.
- **Possible kernel implication:** Audit the relative magnitude of ring-current vs EF kernel contributions per-element. If our calibration up-weights RC for heavy nuclei, we may be missing where the signal lives. Stage-1 atom-type stratification may already reveal this; worth a cross-check.

### 5.7 CH3Shift kernel decomposition for methyl side chains (Sahakyan et al. 2011, JBNMR 50 331-346)

- **Finding:** δ = δ_rc^rot + Δδ_dih + Δδ_ring + Δδ_ma + Δδ_EF + Δδ_dist. Four physical contributions plus a phenomenological distance-based term, all polynomial functions of interatomic distances within a 6.5 Å sphere (active region) plus a 1.8 Å neutral region. Differentiable so usable as restraints. RMSD 0.133-0.198 ppm on Ala/Thr/Val/Leu/Ile methyl ¹H.
- **Kernels affected:** Methyl-specific calculators (currently absent or under-served).
- **Possible kernel implication:** Add a methyl-specific kernel that mirrors CH3Shift's decomposition. Direct ancestor of our kernel approach; coverage gap currently. Chapter material on its own.
- **Cache filename:** `references-text/sahakyan-2011-methyl-chemical-shifts-proteins-text-{1-6}.txt`.

### 5.8 Static a^MD vs averaged a^Xray calibration discipline (Li-Brüschweiler 2012, JBNMR 54 257-265)

- **Finding:** PPM predictor. The load-bearing equation: δ_exp = Σ a_j^(k) · ⟨f_j^(k)(r_n)⟩ + δ_0^(k). The fitted `a` parameters obtained from average X-ray structures vs from MD ensembles are NOT equivalent — even when ⟨r_n⟩_MD = r_Xray — because dynamic averaging gets absorbed into a_j^Xray.
- **Kernels affected:** Not a kernel change. Calibration architecture documentation in `learn/`.
- **Possible implication:** Document the static-a^MD calibration discipline explicitly in `learn/CLAUDE.md` and code comments. Methods-chapter requirement; ensures the calibration architecture is defensible to a reviewer.
- **Cache filename:** `references-text/li-bruschweiler-2012-ppm-shift-predictor-ensembles-text-{1-3}.txt`.

### 5.9 Block-averaged convergence reporting (Markwick / de Gortari / Robustelli, multiple)

- **Finding:** The cache's older lineage uses block averaging to defend convergence claims. Without it, our averaging claims do not survive review.
- **Kernels affected:** Trajectory-scope analysis infrastructure.
- **Possible implication:** Implement `BlockAveragedConvergenceResult` (already in `PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` Idea 1) before the MD re-run, so the analysis is ready when the trajectories land.

### 5.10 Crankshaft cadence-aliasing characterization (Yi-McDermott 2024)

- **Finding:** Picosecond-timescale crankshaft motion (ψ_i / φ_{i+1} anticorrelated, residence times 0.23-0.89 ps) produces 25 ppm ¹⁵N excursions. At 10 ps cadence the motion is aliased.
- **Kernels affected:** No kernel change. Diagnostic / methods-chapter work.
- **Possible implication:** Run a dense-burst diagnostic on a 2-protein subset alongside the main re-run; quantify how much ¹⁵N variance the standard 10 ps cadence loses. Frame as a methods contribution.

### 5.11 Live items that may not be calculator pass but bear on it

- **Cache hygiene:** `references-text/bratholm-2017-procs15-qm-chemical-shielding-refinement-text-{1-7}.txt` and corresponding summary file actually contain Venetos 2023 silicate text, not Bratholm 2017. PDF on disk; needs re-ingestion via `scripts/references/ingest_pdf.sh`.
- **EDMD lineage as thesis context:** Cavalli 2007 CHESHIRE → Robustelli 2009 → Robustelli-Kohlhoff-Cavalli-Vendruscolo 2010 → Robustelli 2012 → Papaleo 2018 → Gadanecz 2026. The lit review chapter sits in this lineage; the thesis can position relative to it.

---

## Section 6 — IUPAC chapter substance

User's framing:

> "Use the IUPAC naming output to feed better stats and a chapter from the mutant work."

Material that follows:

- **Re-stratify the existing 720-protein / 446K-atom mutant calibration by IUPAC atom name** instead of (or alongside) AMBER atom name. Should give finer per-atom-class signal than AMBER stratification.
- **BMRB uses IUPAC natively** — can cross-reference our calibration set against BMRB chemical shifts for atoms we have both for. Validation against experiment, not just DFT.
- **pro-R / pro-S, methylene / methyl, sidechain-by-residue-type stratifications** all become tractable through `AtomReference` and `IupacAtomName`.
- **Atom matching across WT/mutant pairs** becomes more robust. Current matching uses some combination of position + element + topology; IUPAC name + AtomReference is more direct.

Outstanding debts the IUPAC pass should clear (per memory entry `project_iupac_topology_landed_20260426.md`):

- **pro-R / pro-S verification.** Geometric assignment exists but unverified.
- **Variant overrides.** Non-standard residues need the override pathway.
- **DFT-compare completeness.** `MutationDeltaResult` at `src/MutationDeltaResult` currently stores only `delta_shielding` total. Per memory entry `project_dft_compare_calculator_completeness.md`: must add WT original + mutant compared + diamagnetic + paramagnetic, all six tensors plus three deltas; can adopt `AtomReference` for matching.

The mutant validation gate (per memory entry `project_mutant_validation_gate_20260426.md`) is two FULL runs against all 723 WT+ALA pairs — catches drift, coverage gaps, real-data wrongness that smoke tests miss.

The chapter that emerges from this work writes itself given:
- Stage 1 atom-type stratification reveals (sidechain N R²=0.887 vs backbone N R²=0.387 — element-pooling artifact)
- IUPAC re-stratification sharpens those reveals
- BMRB cross-reference validates against experiment
- The DFT-compare completeness work gives full tensor data for a tensor-aware narrative

---

## Section 7 — The "leave with a list" expectation

> "The dept head told me EVERY MSc student in this program leaves with a long list of what their investigation would have done next."

The list-of-next-things is part of the deliverable, not a sign of incompleteness. Items that would naturally land on it given the current scope:

- Crankshaft dense-burst diagnostic at sub-picosecond cadence on a multi-protein subset.
- Full S² / Lipari-Szabo relaxation analysis from MD (the Lai-Brooks 10-20 replicas regime).
- LEGOLAS comparison: predict CS via LEGOLAS on our trajectories alongside our kernel-decomposed predictions.
- ShiftML3-style equivariant ensemble model for full T2 prediction.
- Cross-FF comparison (ff19SB / OPC vs ff14SB / TIP3P) on a single-protein subset.
- More proteins, deeper sampling, longer trajectories.
- ML-on-MD-trajectories for relaxation observable back-calculation per Lesovoy 2025.
- Replica-averaged metainference per Papaleo 2018 if we ever want to pull experimental restraints into MD.
- AMP-BMS/MM 2026 anchored sub-trajectory for DFT-quality MD on a small subset.
- Müntener-style high-throughput de novo BMRB expansion as the calibration data evolves.

This list is itself part of what the lit review and discussion chapters argue from.

---

## Section 8 — What this document is NOT

- Not a plan-of-record. The decisions in here are not committed; they are the post-reframing picture.
- Not a substitute for `md-rerun-685-discussion-priors-2026-04-30.md`, which captures the MD/MDP/replica/cadence priors separately.
- Not a curated subset. Items in §5 are surfaced raw; rankings, "first move" recommendations, cost orderings have been deliberately removed.
- Not a complete capture of the conversation. This is Part 1. Subsequent parts will capture decisions as they land, calculator pass progress, IUPAC chapter work, and the MD re-run plan as it locks.

---

## Document hygiene

When decisions land, capture them in a separate plan-of-record document with this one cross-referenced as the conversation that produced them. Until then, this document is the picture as the user articulated it on 2026-04-30, with the assistant's interpretive shading deliberately backed off so the items are readable on their own.
