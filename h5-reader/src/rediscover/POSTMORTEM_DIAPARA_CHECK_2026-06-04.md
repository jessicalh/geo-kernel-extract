# Dia/Para Full-Tensor Check - 2026-06-04
Substrate: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`; rows=558360; specs=`per_atom_substrate_column_specs.json`.
Sidecars: total=`per_atom_substrate_target_T*`, dia=`per_atom_substrate_target_dia_T*`, para=`per_atom_substrate_target_para_T*`.
Residual: `abs(dia + para - total)` in ppm; distribution is elementwise over emitted T components.

| group | shape | median | p99 | max |
|---|---:|---:|---:|---:|
| T0 | (558360,) | 3.333333e-4 | 6.666667e-4 | 1.000000e-3 |
| T1 | (558360,3) | 2.342571e-14 | 1.000000e-3 | 1.000000e-3 |
| T2 | (558360,5) | 5.952183e-14 | 1.414214e-3 | 1.632993e-3 |

Row-max distributions: T0 median/p99/max=3.333333e-4/6.666667e-4/1.000000e-3; T1=5.000000e-4/1.000000e-3/1.000000e-3; T2=7.071068e-4/1.414214e-3/1.632993e-3.
Inverse 3x3 tensor residual from T0/T1/T2: median/p99/max=2.845131e-14/1.000000e-3/1.000000e-3; all nine Cartesian component maxima are 1.000000e-3.
Verdict: split SOUND. T2's 1.632993e-3 is the expected isometric-basis projection of ORCA 1e-3 printed tensor rounding, not component drift.
Affected scope: none for build4 dia/para split sidecar targets; TOTAL target remains raw total and was not implicated.
Scratch evidence: `/tmp/rediscover-runs/2026-06-04-diapara-check/diapara_full_tensor_check.json`.
