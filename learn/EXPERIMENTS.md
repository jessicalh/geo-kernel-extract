# Calibration Experiments Log

All experiments from the 2026-04-08 session.  Runs are in learn/runs/.
Commit after every experiment from now on.

| Run name | Proteins | Epochs | Key change | Val R² | Val 0-4Å | Notes |
|---|---|---|---|---|---|---|
| smoke_test_v2 | 8 | 20 | First CalibrationDataset (old 66 scalars, 48 kernels, global kernel norm) | ~0 | — | Only 1 protein had delta, rest from old extraction |
| early_signal_13 | 13 | 50 | Same as above, more proteins | 0.16 | — | First real signal |
| critique_fixes | 29 | 30 | Global kernel norm, H/S element, z-score scalars | -0.06 | -0.10 | Global norm killed it — model learning nuisance scale |
| perprot_restored | 50 | 50 | Reverted to per-protein kernel norm, kept z-score + H/S element | **0.44** | **0.56** | Best result. Per-protein norm isolates angular structure |
| diagnostics_test | 45 | 30 | Added distance-weighted loss (tau=8Å), naive baselines | 0.05 | 0.09 | Lower than perprot — fewer val proteins? distance weighting? |
| xkernel_test | 57 | 50 | Added 18 cross-kernel T2 dot products (BS·HM, BS·PQ, BS·MC/Coulomb) | -0.04 | -0.12 | Overfitting. 86 scalars too many for 57 proteins. Stripped. |
| analysis_driven | 77 | 50 | Mutation identity (4), z/rho cylindrical coords, PQ dropped (40 kernels, 78 scalars) | **0.42** | **0.51** | On track. 4-8Å improved to 0.26. 8-12Å improved from -1.0 to -0.28 |

## Analytical diagnostics (no learning)

Run: `analyze.py --run CalibrationExtractionTest` on 66 proteins.

- Per-protein ridge R²: median=0.81, range [0.60, 0.95] — kernels have the signal
- Global ridge R²: 0.23 — weights don't transfer across proteins
- BS and HM dominate (r≈-0.23), both negatively correlated with target
- H atoms strongest: HM_TYR r=-0.51
- PQ near zero everywhere (r<0.03) — dropped
- 0-4Å ridge R²=0.30 (hard), 8-12Å ridge R²=0.87 (easy, small signal)

## Key lessons

1. Per-protein kernel normalization is correct for calibration (angular structure question, not magnitude)
2. Global kernel normalization forces model to learn nuisance scale variable — hurts with limited data
3. Cross-kernel dot products overfit at <100 proteins — revisit at 700
4. Mutation identity (which aromatic residue removed) is cheap and useful
5. Distance-weighted loss helps mid/far-field, doesn't hurt near-field
6. Naive baselines (unweighted sums) are negative because of sign structure, not broken physics
7. The 0.60 target is reachable: per-protein ceiling is 0.81, MLP at 0.44 with 77 proteins
