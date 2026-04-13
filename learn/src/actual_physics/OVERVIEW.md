# actual_physics/ — Kernel Validation and Calibration

Physics-grounded analysis of the geometric kernels.  Everything here
answers one question: **do the kernels see what NMR theory says they
should see?**

No MLP, no gating, no hidden layers.  Ridge regression and tensor
geometry.  The outputs are tables of calibrated physical constants,
not model predictions.

## Scripts (in order of use)

### 1. element_physics.py — Per-element kernel decomposition

Tests the literature prediction: ring current dominates H, EFG
dominates C (Boyd & Skrynnikov 2002; Sahakyan & Vendruscolo 2013).
Stratifies 91 kernels by element and physics group.  Forward
selection per element shows which kernels matter where.

**Key finding:** H needs 2 kernels for 90% of R².  N needs 15.

Output: `output/secondary/element_physics/`

### 2. twenty_realities.py — Analytical reality tests (13-20)

Eight analytical tests against known NMR physics: distance falloff,
angular peaking, BS-HM agreement, Buckingham sign, mechanical mutant
controls, element weighting.  Ridge-based.

Output: `output/secondary/twenty_realities/`

### 3. tensor_realities.py — Tensor-direct reality tests (T21-T28)

Eight tests on the raw T2 vectors — no fitting.  Distance slopes
(BS -3.04, PQ -5.05), BS-HM magnitude ratio (1.009), angular
dependence, inter-kernel independence.  Three-significant-figure
agreement with multipolar theory.

Output: `output/secondary/tensor_realities/`

### 4. export_for_r.py — Per-atom data export for R figures

Exports 69K-row CSV with per-atom tensor geometry (dist, theta,
rho, z, kernel magnitudes, cosine similarities) for the R plotting
script at `learn/R/twenty_eight_realities.R`.

Output: `output/actual_physics/atom_tensor_data.csv`

### 5. clean_calibration.py — 55-kernel per-element ridge

Strips to 55 core kernels (ring-type + bond-category + totals +
EFG).  No per-ring residuals, no RBF, no angular basis.  Per-element
ridge with forward selection and cumulative R².

**Key finding:** Weighted per-element R² = 0.690 with 55 kernels,
no model.

Output: `output/actual_physics/clean_calibration/`

### 6. physics_calibration.py — Progressive fair-scalar calibration

Adds per-protein normalization and tests each scalar interaction
independently: kernel scales, mutation type (which residue removed),
MOPAC valence, bond order, molecular dipole.

**Key finding:** Scales (+0.06) and mutation type (+0.12-0.22) help.
Valence, bond order, dipole add exactly zero.  Weighted R² = 0.818
with fair physics only.

Output: `output/actual_physics/physics_calibration/`

### 7. per_element_calibration.py — MLP comparison (historical)

Trains tiny per-element models (hidden=4) to compare gated MLP
against ridge.  Shows the MLP adds nothing for H, modest gains for
C, and the pooled model (0.61) is worse than per-element ridge
(0.718).

Output: `output/actual_physics/calibration/`

## The calibration result

55 core kernels, per-protein normalization, element stratification,
mutation type identity.  Ridge regression.

| Element | Base | + scales | + mutation type | Fair set |
|---------|------|----------|----------------|----------|
| H       | 0.875 | +0.063  | +0.027         | **0.949** |
| C       | 0.520 | +0.069  | +0.121         | **0.691** |
| N       | 0.450 | +0.059  | +0.220         | **0.704** |
| O       | 0.407 | +0.048  | +0.195         | **0.632** |
| Weighted | 0.684 |         |                | **0.818** |

No MLP needed.  The ridge coefficients are the physical constants.

## Related files

- `learn/R/twenty_eight_realities.R` — 7 publication figures
- `learn/docs/twenty_eight_realities_2026-04-10.md` — full writeup
- `learn/docs/element_physics_2026-04-10.md` — element analysis
- `learn/docs/calibrated_weights_2026-04-10.md` — weight interpretation
- `learn/docs/realities_latex.tex` — thesis LaTeX section

## Running everything

```bash
cd learn/src

# Element physics + forward selection
python3 -c "
from mutation_set.config import load_config
from secondary.loader import setup_sdk
from actual_physics.element_physics import run
cfg = load_config('calibration.toml'); setup_sdk(cfg); run(cfg)
"

# 28 realities (analytical + tensor)
python3 -c "
from mutation_set.config import load_config
from secondary.loader import setup_sdk
from actual_physics.twenty_realities import run
from actual_physics.tensor_realities import run as run_t
cfg = load_config('calibration.toml'); setup_sdk(cfg)
run(cfg); run_t(cfg)
"

# Clean calibration (55 core kernels, per-element ridge)
python3 -c "
from mutation_set.config import load_config
import sys; cfg = load_config('calibration.toml')
sys.path.insert(0, str(cfg.paths.sdk))
from actual_physics.clean_calibration import run
run(cfg)
"

# Physics calibration (progressive: +scales, +mutation type, +valence, +dipole)
python3 -c "
from mutation_set.config import load_config
import sys; cfg = load_config('calibration.toml')
sys.path.insert(0, str(cfg.paths.sdk))
from actual_physics.physics_calibration import run
run(cfg)
"

# R figures
cd learn/R && Rscript twenty_eight_realities.R
```

All scripts have been run on 720 proteins (Stage1Results extraction,
446K atoms).  Output in src/output/actual_physics/.  Additional
scripts added in 2026-04-13: full_space_analysis.py, orca_dia_para.py,
completeness_checks.py.  Full reproduction via
stage1-mutations/analysis/run_all.sh.
