# Compute Budget

## Extraction cost

| Scenario | Frames | Time/frame | Total CPU-time | Wall (3 machines) |
|---|---|---|---|---|
| 6 initial-samples per protein | 4,110 | 1.3s | 1.5 hrs | 30 min |
| 50-frame subsample per protein | 34,250 | 1.3s | 12.4 hrs | 4 hrs |
| Full trajectory (5005 frames) | 3,428,425 | 1.3s | 1,238 hrs | 17 days |

The 6-pose pilot is nearly free. A 50-frame subsample (every 100th
frame) gives ensemble statistics at 1% of the full cost. Full
trajectory extraction is a significant commitment.

No additional charge calculations needed. The EFG geometry dominance
finding (EFG_GEOMETRY_DOMINANCE.md) shows topology charges carry the
same directional signal as MOPAC (cos 0.9296 vs 0.9323). The
geometry-only extraction has everything we need.

## Analysis memory budget (128 GB machine)

| Data | Size | Notes |
|---|---|---|
| Ensemble-averaged T0, 685 proteins | ~500 MB | Mean + var + skew per atom per kernel |
| Centred RefDB shifts | ~50 MB | 713K matched shifts |
| Per-ring contributions (6 poses) | ~2 GB | Sparse, ~7K pairs x 6 poses x 685 proteins |
| Per-ring contributions (50 frames) | ~15 GB | Sparse |
| Kernel covariance matrices | ~4 GB | 685 proteins x ~1500 atoms x 20x20 |
| Working headroom | ~100 GB | For R/Python intermediate objects |

6-pose and 50-frame analyses fit comfortably. Full-trajectory
per-atom covariance would need streaming accumulation (the
EnsembleConformation design) -- never hold all frames in memory.

## GPU budget (4x RTX 5090)

The analytical layers don't need GPU. The GNN does. These run in
parallel -- analytical work on CPU while GNN trains on GPU. No
competition for resources.

## R machine (2 TB, short loan)

If available, use for:
- Massive cross-protein PCA (Layer 7) with all atoms in memory
- Publication-quality figures with full data density
- Any analysis that needs >128 GB working set

Plan the analysis to run unattended. Prepare scripts, ship data,
run, collect results. Don't iterate on the big machine.

## Timeline sketch (speculative)

| Week | Analytical | GNN |
|---|---|---|
| 1-2 | Layer 0-2: census, centering, pilot on 6 poses | Piece 1-3: types, registry, graph builder |
| 3-4 | Layer 3-4: T0 correlation, distance stratification | Piece 4-5: equivariance tests, layers |
| 5-6 | Extract 50-frame subsample, Layer 5: covariance | Piece 6-7: loss, evaluation |
| 7-8 | Layer 6-7: per-ring tracking, emergent clusters | Piece 8: training loop, first runs |
| 9-10 | Full trajectory extraction (background) | GNN training on 6-pose data |
| 11-12 | Final figures, cross-validation of findings | GNN on full ensemble, final results |
| 13 | Thesis writing | Thesis writing |

This is speculative. Reality will reorder it. The key constraint
is: analytical findings by week 8 so the GNN has design guidance.
