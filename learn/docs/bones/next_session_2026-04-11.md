# Next Session Prompt

Read spec/INDEX.md first.  Then read learn/src/actual_physics/OVERVIEW.md
— it describes the physics calibration pipeline built on 2026-04-10.
Then read learn/docs/twenty_eight_realities_2026-04-10.md for the full
validation document.

## Where we are

The calibration is done.  55 core geometric kernels (rank-2 traceless
symmetric tensors), per-protein normalised, per-element ridge
regression with two fair categorical scalars (kernel scale factors,
mutation type identity).  Weighted per-element R² = 0.818 across
69,080 atoms in 110 proteins.  No MLP, no gating, no hidden layers.
The ridge coefficients ARE the calibrated physical constants.

28 NMR physics realities verified analytically and tensor-directly,
each grounded in a specific literature citation.  Seven publication-
quality R figures produced.  LaTeX section drafted.

## Key findings from this session

1. **Per-element physics is real.**  Ring current dominates H (R²=0.949),
   EFG dominates C (R²=0.691), N is multi-mechanism (R²=0.704),
   O is mixed (R²=0.632).  Pooled analysis hides all of this.

2. **The MLP was an engineering detour.**  Per-element ridge (0.818)
   beats the gated MLP (0.61).  Four tiny per-element models
   (hidden=4) also beat pooled hidden=64.

3. **Kernel distance dependence at 3 significant figures.**
   BS slope=-3.04 (theory -3), PQ=-5.05 (theory -5), HM=-3.06.
   HM/BS magnitude ratio 1.009 ± 0.037.

4. **Valence, bond order, and dipole add zero.**  The kernels encode
   everything these scalars know.

5. **Mutation type matters, especially for heavy atoms.**  +0.22 R²
   for nitrogen.  This is Case 1995's ring-type-specific intensity
   factors, extended to T2.

## What to do next

### Immediate (before thesis writing)

1. **Restart the AzimuthalExtraction** (613 proteins remaining).
   Run in tmux:
   ```bash
   cd /shared/2026Thesis/nmr-shielding
   nohup python3 -u learn/extract.py --run AzimuthalExtraction --resume \
     > calibration/features/AzimuthalExtraction/extract_stdout.log 2>&1 &
   ```
   When complete, rerun all actual_physics/ scripts for final 723-protein
   numbers.

2. **Steal from analyze.py:**
   - Section 7 (OLS vs ridge per element) — quantifies whether more
     signal is in-span or out-of-span for each element
   - Section 6 (leave-one-out ablation per element) — which kernels
     are individually indispensable vs redundant

3. **Pin down the thesis calibration table:**
   55 kernels × 4 elements × ~5 mutation types = ridge coefficient
   table with confidence intervals.  This is the deliverable.

### For the thesis chapter

4. **The 28 realities are Section 4.x** — "The Tensor Decomposition
   of the Shielding Perturbation: 723 Mutants."  Each reality is
   one paragraph + figure/table.  Realities 1-12 come from the
   existing secondary analysis tools (divergence, ring_strata,
   kernel_structure, ablation).  Realities 13-28 come from
   actual_physics/ scripts.

5. **The progressive calibration is Section 4.y** — shows what
   physics helps (element, mutation type, kernel scales) and what
   doesn't (valence, bond order, dipole).

6. **Limitations section is honest** — per-protein ceiling (0.81),
   HIE modeling (0.062), near-field accuracy, no structural
   relaxation.

### Do NOT do

- Do not try to improve R².  0.818 is the answer.
- Do not add features.  Valence/bond/dipole are zero.
- Do not train larger models.  Ridge IS the model.
- Do not say "wait for more data."  110 proteins already show
  the physics; 723 will sharpen the numbers, not change them.

## Files to know

- `learn/src/actual_physics/` — 7 analysis scripts + OVERVIEW.md
- `learn/R/twenty_eight_realities.R` — 7 publication figures
- `learn/docs/twenty_eight_realities_2026-04-10.md` — full validation
- `learn/docs/element_physics_2026-04-10.md` — element decomposition
- `learn/docs/calibrated_weights_2026-04-10.md` — weight vectors
- `learn/docs/realities_latex.tex` — thesis LaTeX
- `references/` — 4 PDFs (Boyd 2002, Sahakyan 2013, Case 1995,
  Buckingham 1960) + ANNOTATED_BIBLIOGRAPHY.md (90+ papers)
