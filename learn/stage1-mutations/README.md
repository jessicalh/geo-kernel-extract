# stage1-mutations/ — Thesis Chapter: Mutation Calibration

Established 2026-04-12.  Analysis on 720 proteins, 446,006 atoms
from Stage1Results extraction.

## What this is

The thesis chapter answering: what does an always-T2 kernel feature
extraction tell us about the physics of shielding?

Read docs/stage1_plan.md for the full plan.
Read docs/INDEX.md for the document index.

## Reproducing

```bash
cd learn/src
bash ../stage1-mutations/analysis/run_all.sh    # all analyses
python3 ../stage1-mutations/analysis/verify_numbers.py  # check
cd learn
Rscript stage1-mutations/analysis/stage1_figures.R  # figures
```

## Contents

```
analysis/
  run_all.sh           — runs all 11+ analysis scripts
  verify_numbers.py    — checks cited numbers against JSON output
  stage1_figures.R     — publication figures from CSVs

notes/
  dimension_inventory_raw.md    — full inventory of dimensions
  argument_chain.md             — the thesis claim and its grounding
  literature_grounding.md       — what's cited, extended, or novel
  dia_para_decomposition.md     — WT-ALA delta of DFT dia/para
  new_arrays_assessment.md      — new Stage1Results arrays tested
  dimensionality_honest.md      — H=20, C=6, N=3, O=12
  physics_group_grounding.md    — each group's physics, per element
  normalisation_physics.md      — what normalisation does and doesn't
  dimension_tree.md             — full tree with cosines and citations
  physics_analysis_bridge.md    — bridging table: groups × elements
  reproducible_analysis_plan.md — 4-phase defensible methodology
  completeness_and_nonlinear.md — LPOCV, bootstrap, RF nonlinear check

figures/                        — PDFs from stage1_figures.R
references/                     — papers (to be populated)
```

## Key numbers (720 proteins)

| Element | n | R² (norm) | LPOCV R² | Dims | RF delta |
|---------|--------|---------|---------|------|---------|
| H | 230,135 | 0.856 | 0.844 | 20 | +0.002 |
| C | 133,488 | 0.484 | 0.456 | 6 | **+0.128** |
| N | 39,954 | 0.267 | 0.172 | 3 | **+0.169** |
| O | 42,429 | 0.304 | 0.213 | 12 | +0.013 |

RF nonlinear signal follows paramagnetic ordering (Saito 2010):
N > C >> O > H.  All 20/20 bootstrap stable.

## Do not modify

This directory is frozen for the thesis.  New analysis goes in
new directories.  The numbers here must match verify_numbers.py.
