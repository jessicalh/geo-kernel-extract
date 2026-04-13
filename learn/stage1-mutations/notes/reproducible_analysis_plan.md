# Reproducible Analysis Plan

2026-04-13.  Turning the discovery process into a defensible
methodology on the final 723-protein data.

---

## The problem

The calibration was built through successive analysis: forward
selection, ablation, "does this kernel help?", element stratification
discovered partway through.  The conclusions are correct but the
process is not reproducible — it was problem-solving, not methodology.

A thesis examiner can ask: "you selected 55 kernels by iterative
analysis on 110 proteins.  How do I know this isn't overfitted to
your discovery path?"

The answer: a principled pipeline that starts from everything and
lets the data decide.  Run it on the full 723-protein Stage1Results
extraction.  If it arrives at the same conclusions, the iterative
discovery was correct.  If it arrives at different conclusions, we
deal with that honestly.

---

## Phase 1: Full-space characterisation (no selection)

All kernels.  All atoms.  All elements.  No pre-filtering.

### What to compute

- **Eigenspectrum** of the full kernel T2 space, per element.
  Raw and per-protein-normalised.  How many effective dimensions?
  Do eigenvalues come in blocks of 5 (one per L=2 component)?

- **Inter-kernel cosine matrix.**  For each pair of kernel families,
  mean |cos| across all atoms.  This is the empirical independence
  structure.  Visualise as a heatmap per element.

- **Per-kernel magnitude distributions.**  Median, IQR, fraction
  nonzero.  Which kernels are active at which distances?

- **PCA vs PLS.**  PCA finds directions of maximum variance.  PLS
  (partial least squares) finds directions of maximum covariance
  with the target.  If they agree, the dominant signal is also the
  predictive signal.  If they disagree (as they do for H), the
  prediction lives in the angular fine structure.

### Tools

- Python: assemble full kernel arrays via existing SDK + kernels.py.
  Export per-element design matrices as CSV/feather.
- R: prcomp() for PCA, pls::plsr() for PLS.  Compare loading
  structures.

### Output

- Per-element eigenspectrum plot (raw vs normalised)
- Inter-kernel cosine heatmap per element
- PCA vs PLS dimension comparison (loadings on physics groups)

---

## Phase 2: Data-driven selection

Let the data choose which kernel families matter.  Do not use
forward selection (greedy, path-dependent).  Use methods that
examine the full space simultaneously.

### Group LASSO

R package glmnet with group penalties.  Define groups by physics
family:

| Group | Kernels | Physics |
|-------|---------|---------|
| BS | BS_PHE, BS_TYR, ... (8) | Biot-Savart ring current |
| HM | HM_PHE, HM_TYR, ... (8) | Haigh-Mallion ring current |
| Disp | Disp_PHE, ... (8) | London dispersion |
| PQ | PQ_PHE, ... (8) | Pi-quadrupole EFG |
| MC | MC_backbone, ... (5) | McConnell bond anisotropy |
| MopacMC | MopacMC_backbone, ... (5) | MOPAC bond-order-weighted MC |
| Totals | MC_total, Coulomb_total, ... (7) | Calculator totals |
| EFG | EFG_bb, EFG_aro, ... (6) | Electric field gradient |

The LASSO path (lambda vs group inclusion) shows the order in which
physics families enter as regularisation decreases.  This is the
data-driven version of forward selection — but simultaneous, not
greedy.

Compare the LASSO entry order to our forward selection order.
Agreement validates the forward selection.  Disagreement is a
finding.

### Elastic net path

LASSO is L1 (sparse).  Ridge is L2 (dense).  Elastic net
interpolates.  The optimal alpha (L1/L2 ratio) tells you whether
the kernel space is truly sparse (a few families dominate) or dense
(many families contribute).  Per element.

### Cross-validated R-squared curve

For each level of regularisation, 5-fold CV R-squared.  The curve
shape tells you: how much signal is in the top families, and when
do you start overfitting?

### Tools

- R: glmnet::cv.glmnet() with group argument.  grpreg package for
  group LASSO specifically.
- Python: sklearn.linear_model.ElasticNetCV for quick check.

### Output

- LASSO path plot per element (lambda vs group coefficient)
- Elastic net optimal alpha per element
- CV R-squared curve per element

---

## Phase 3: Nonlinear check

The calibration uses linear ridge.  Is linearity sufficient, or
does the kernel space contain nonlinear signal that ridge misses?

### Kernel ridge regression

sklearn.kernel_ridge.KernelRidge with RBF kernel.  Compare to
linear ridge on the same CV splits.  The gap (if any) is the
nonlinear signal.

### Random forest feature importance

sklearn.ensemble.RandomForestRegressor or R randomForest.
Permutation importance per kernel.  Does the importance ranking
match the ridge coefficient ranking?  If yes, the linear
decomposition captures the full structure.  If specific kernels
gain importance in the nonlinear model, their functional form may
be inadequate (the kernel computes the wrong geometry).

### What we expect

The MLP (hidden=64) got 0.61 pooled.  Per-element ridge gets 0.818
weighted.  The MLP is worse despite being nonlinear.  We expect
nonlinear methods to confirm that ridge is sufficient — that the
kernel T2 space is genuinely linear in its relationship to the DFT
target.

If nonlinear methods DO find signal for specific elements (likely
N, which has 5 competing mechanisms), that's an honest finding:
"the kernel interactions for nitrogen are not purely additive."

### Tools

- Python: sklearn for kernel ridge and random forest.
- R: randomForest, kernlab::ksvm for kernel methods.

### Output

- Linear vs nonlinear R-squared per element (table)
- Random forest importance vs ridge coefficient correlation
- If nonlinear wins somewhere: which kernels, which elements

---

## Phase 4: Stability

### Leave-protein-out cross-validation

Not leave-atom-out (atoms within a protein are correlated).
Leave-one-protein-out (LOPOCV) or leave-10%-proteins-out.
Per-element R-squared distribution across folds.

### Bootstrap confidence intervals

Resample proteins (not atoms) with replacement, refit ridge,
collect coefficient distributions.  95% CI on each kernel's
ridge coefficient per element.  Kernels whose CI includes zero
are not significantly contributing.

### Regularisation sensitivity

Ridge R-squared as a function of lambda.  Is there a plateau
(robust) or a knife-edge (fragile)?  Per element.

### Tools

- R: boot package for bootstrap, custom LOPOCV loop.
- Python: sklearn cross_val_score with GroupKFold (groups = proteins).

### Output

- LOPOCV R-squared distribution per element (violin/box plot)
- Bootstrap CI table: kernel × element × (coeff, lower, upper)
- Lambda sensitivity curve per element

---

## Practical workflow

1. **Python export script** (new, in src/actual_physics/):
   Load all 723 proteins via SDK.  Assemble full kernel arrays.
   Export per-element CSV files with columns:
   protein_id, atom_index, kernel_1_c1..c5, ..., kernel_55_c1..c5,
   target_c1..c5, distance, element.
   Also export the kernel name list and physics group assignments.

2. **R analysis script** (new, in R/):
   Reads the CSVs.  Runs all four phases.  Produces figures and
   tables.  Self-contained — does not call Python.

3. **Comparison script** (Python):
   Reads the R output and the existing actual_physics/ output.
   Produces a concordance table: does the principled analysis agree
   with the iterative analysis?

---

## What this gives the thesis

- "We performed group LASSO, elastic net, PLS, and nonlinear
  comparison on the full 55-kernel space across 723 proteins."
- "The data-driven selection confirms the per-element physics
  decomposition found through iterative analysis."
- "Linear ridge is sufficient: nonlinear methods add [nothing /
  X for nitrogen]."
- "Leave-protein-out CV gives R² = X ± Y per element.  Bootstrap
  CIs on the calibrated coefficients are [table]."
- Every number has a proper statistical pedigree.
