# McConnell literature-scaled de-circularisation

> **Evidence archive — not current truth (trued 2026-06-04).** Preserve as
> analysis provenance; check `../NOW.md` and corrected `../STATE.md` before
> quoting quantitative claims.

Read-only analysis of emitted `mc_lit_T0` and `mc_lit_T2_local_*` columns from `broad_backbone_aggregated.csv`. Python only correlates emitted columns against emitted DFT targets; it does not open `trajectory.h5`, apply Delta-chi, rebuild tensors, or call the frozen change-of-basis helper.

Delta-chi values are provisional single-family Williamson-Asakura values from `src/rediscover/MCCONNELL_DCHI_LITERATURE.md`: peptide C=O +2.41, peptide C-N -5.42, sidechain C=O +2.41, aromatic 0 in 10^-6 cm^3/mol. Aromatic is zero because RING carries the pi current.

Input out-dir: `/tmp/rediscover-mc-lit-fresh`

CSV artifact: `/tmp/rediscover-mc-lit-decirc/mcconnell_literature_decirc.csv`

Audit artifact: `/tmp/rediscover-mc-lit-decirc/mcconnell_literature_decirc_audit.json`

## T2 lead

| stratum | axis | rows | atoms | N_eff(atom) | N_eff(AR1 frames) | gamma_lit +/- SE | component_r | absT2_r | LOAO R2 | bucket |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| HN | within | 26000 | 52 | 46.7356 | 2.2127e+04 | -0.3495 +/- 0.3107 | -0.0299 | 0.2748 | -0.0005 | form-recovered-scale-fitted |
| HN | between | 52 | 52 | 46.7356 | 2.2127e+04 | -1.1534 +/- 0.7969 | 0.0781 | 0.0943 | -0.0390 | form-recovered-scale-fitted |
| N | within | 27000 | 54 | 52.6814 | 2.5922e+04 | 0.2773 +/- 0.0312 | 0.0258 | 0.1724 | 0.0006 | form-recovered-scale-fitted |
| N | between | 54 | 54 | 52.6814 | 2.5922e+04 | 2.6041 +/- 0.9647 | 0.6533 | 0.9138 | 0.0640 | zero-circularity recovered law |
| CA | within | 27000 | 54 | 27.1593 | 1.9400e+04 | 1.3567 +/- 0.5556 | 0.0342 | 0.4066 | 0.0008 | zero-circularity recovered law |
| CA | between | 54 | 54 | 27.1593 | 1.9400e+04 | 2.9034 +/- 0.6368 | 0.1163 | 0.1810 | -0.0007 | form-recovered-scale-fitted |
| C | within | 27000 | 54 | 52.7654 | 2.5926e+04 | -2.1368 +/- 0.0225 | -0.6176 | 0.7549 | 0.3813 | form-recovered-scale-fitted |
| C | between | 54 | 54 | 52.7654 | 2.5926e+04 | -1.9164 +/- 0.4109 | -0.3987 | 0.5090 | 0.3097 | form-recovered-scale-fitted |
| O | within | 27000 | 54 | 53.6927 | 2.6176e+04 | 2.4929 +/- 0.2294 | 0.0343 | 0.0615 | 0.0012 | form-recovered-scale-fitted |
| O | between | 54 | 54 | 53.6927 | 2.6176e+04 | 38.0253 +/- 30.3759 | 0.5374 | 0.3486 | -0.0432 | zero-circularity recovered law |
| HA | within | 29000 | 58 | 37.7413 | 2.2477e+04 | 0.3478 +/- 0.2521 | 0.0310 | 0.2838 | -4.9731e-05 | form-recovered-scale-fitted |
| HA | between | 58 | 58 | 37.7413 | 2.2477e+04 | 1.8259 +/- 0.3792 | -0.0309 | -0.1931 | 0.0299 | form-recovered-scale-fitted |

## T0 trace audit

| stratum | axis | rows | atoms | gamma_lit +/- SE | scalar_r | max_abs_mc_lit_T0_ppm | bucket |
|---|---:|---:|---:|---:|---:|---:|---|
| HN | within | 26000 | 52 | NA +/- NA | -0.0049 | 2.149612e-16 | insufficient |
| HN | between | 52 | 52 | NA +/- NA | -0.1547 | 2.149612e-16 | insufficient |
| N | within | 27000 | 54 | NA +/- NA | 0.0054 | 4.774501e-15 | insufficient |
| N | between | 54 | 54 | NA +/- NA | -0.1391 | 4.774501e-15 | insufficient |
| CA | within | 27000 | 54 | NA +/- NA | -0.0018 | 1.282973e-16 | insufficient |
| CA | between | 54 | 54 | NA +/- NA | -0.0044 | 1.282973e-16 | insufficient |
| C | within | 27000 | 54 | NA +/- NA | -0.0054 | 5.935501e-15 | insufficient |
| C | between | 54 | 54 | NA +/- NA | -0.2037 | 5.935501e-15 | insufficient |
| O | within | 27000 | 54 | NA +/- NA | -0.0030 | 3.025213e-15 | insufficient |
| O | between | 54 | 54 | NA +/- NA | -0.0535 | 3.025213e-15 | insufficient |
| HA | within | 29000 | 58 | NA +/- NA | 0.0043 | 9.916836e-17 | insufficient |
| HA | between | 58 | 58 | NA +/- NA | 0.1118 | 9.916836e-17 | insufficient |

## Notes

- `mc_lit_T0` is a trace audit channel; the emitted PCS tensor is traceless, so gamma is usually undefined.
- `gamma_lit` is a scale diagnostic only. Correlations are unfitted component and magnitude comparisons.
- Confidence is within-protein only: leave-one-atom jackknife SE plus AR(1)-deflated frame N_eff. Thin strata are flagged, not forced.
- Verdict buckets are mechanical diagnostics; interpretation is reserved for the lead.
