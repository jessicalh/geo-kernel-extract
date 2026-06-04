# McConnell Delta-Chi Calibration

> **Evidence archive — not current truth (trued 2026-06-04).** Preserve as
> analysis provenance; check `../NOW.md` and corrected `../STATE.md` before
> quoting quantitative claims.

Status: ORCA-free layer-2 calibration from emitted substrate only. The calibrated Delta-chi values below are provisional DFT-calibrated coefficients with within-protein confidence, not a de-circularisation claim.

Input out-dir: `/tmp/rediscover-mc-lit-fresh`
Single-gamma CSV: `/tmp/rediscover-mc-lit-decirc/mcconnell_literature_decirc.csv`
Artifacts: `/tmp/rediscover-mc-dchi-calibration`

Relation used: the C++ McConnell source scales the local geometric tensor as `T2_ppm = -(1e24/N_A) * q / 3 * T2_geom`. With `q = Delta-chi/(10^-6 cm^3 mol^-1)`, `1e24/N_A = 1.660539067`, so `beta_geom = -0.553513022 * q` and `q_cal = -beta_geom / 0.553513022`.

The broad source CSV does not carry aggregate per-category unweighted columns. It does carry emitted source tensor rows with `bond_category`; this report source-sums those emitted tensor components by `row_id` and category, then removes the exact C++ WA scalar. It does not evaluate a distance, angle, tensor projection, H5 read, ORCA job, or emitter.

## Calibrated Delta-Chi Lead

| stratum | axis | q_CO | q_CN | q_sidechain_CO | CO/CN | absT2 r | R2 | LOAO R2 | N_eff atom | nearest CO / CN / SC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| HN | within | 6.665 +/- 2.295 | 6.348 +/- 1.451 | 15.082 +/- 2.286 | 1.050 | 0.399 | 0.031 | 0.028 | 44.938 | Abraham aliphatic amide / ApSimon / Schneider carbonyl |
| HN | between | 0.926 +/- 6.145 | 8.527 +/- 5.949 | 12.503 +/- 19.385 | 0.109 | 0.486 | 0.014 | -0.056 | 44.938 | WA / ApSimon / ApSimon carbonyl |
| N | within | -40.076 +/- 6.887 | -1.320 +/- 0.180 | -18.770 +/- 17.613 | 30.359 | 0.511 | 0.006 | 0.005 | 52.685 | WA / ApSimon / WA peptide-like |
| N | between | 37.602 +/- 25.975 | -14.178 +/- 6.169 | -81.115 +/- 116.072 | -2.652 | 0.830 | 0.122 | 0.052 | 52.685 | Schneider / Abraham aliphatic amide / WA peptide-like |
| CA | within | -13.017 +/- 6.987 | -13.849 +/- 3.845 | 28.361 +/- 8.582 | 0.940 | 0.396 | 0.007 | 0.005 | 32.327 | WA / Abraham aliphatic amide / Schneider carbonyl |
| CA | between | -17.758 +/- 13.021 | -29.504 +/- 8.210 | 16.609 +/- 59.360 | 0.602 | 0.169 | 0.068 | 0.004 | 32.327 | WA / Abraham aliphatic amide / Schneider carbonyl |
| C | within | -5.730 +/- 0.146 | 11.549 +/- 0.126 | 9.081 +/- 6.445 | -0.496 | 0.755 | 0.382 | 0.382 | 52.972 | WA / ApSimon / Abraham aliphatic amide |
| C | between | -6.154 +/- 9.455 | 10.387 +/- 2.018 | -39.110 +/- 61.965 | -0.592 | 0.499 | 0.355 | 0.301 | 52.972 | WA / ApSimon / WA peptide-like |
| O | within | 4.561 +/- 0.537 | -74.853 +/- 10.122 | 14.583 +/- 43.155 | -0.061 | 0.458 | 0.003 | 0.003 | 53.720 | Abraham aliphatic amide / Abraham aliphatic amide / Schneider carbonyl |
| O | between | 329.066 +/- 294.918 | -35.548 +/- 42.459 | 144.721 +/- 231.329 | -9.257 | 0.415 | 0.173 | -0.095 | 53.720 | Schneider / Abraham aliphatic amide / Schneider carbonyl |
| HA | within | 6.854 +/- 1.147 | 0.681 +/- 1.213 | 12.103 +/- 2.902 | 10.064 | 0.258 | 0.026 | 0.023 | 42.134 | Abraham aliphatic amide / ApSimon / ApSimon carbonyl |
| HA | between | 4.960 +/- 3.541 | -9.793 +/- 2.824 | 1.008 +/- 12.983 | -0.506 | 0.419 | 0.069 | 0.014 | 42.134 | Abraham aliphatic amide / Schneider / WA peptide-like |

WA assumes CO/CN = -0.445. The fitted CO/CN ratios vary by stratum and axis; rows with sign flips or very large SEs should be treated as diagnostics rather than transferable constants.

## Fit Diagnostics

| stratum | axis | rows | atoms | active | rank | component r | absT2 r | R2 | LOAO R2 | AR1 N_eff | intercept norm | thin |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| HN | within | 26000 | 52 | 52 | 3 | 0.177 | 0.399 | 0.031 | 0.028 | 2.177e+04 | NA | ok |
| HN | between | 52 | 52 | 52 | 3 | 0.812 | 0.486 | 0.014 | -0.056 | 2.177e+04 | 7.181 | ok |
| N | within | 27000 | 54 | 54 | 3 | 0.077 | 0.511 | 0.006 | 0.005 | 2.597e+04 | NA | ok |
| N | between | 54 | 54 | 54 | 3 | 0.963 | 0.830 | 0.122 | 0.052 | 2.597e+04 | 90.789 | ok |
| CA | within | 27000 | 54 | 54 | 3 | 0.082 | 0.396 | 0.007 | 0.005 | 2.013e+04 | NA | ok |
| CA | between | 54 | 54 | 54 | 3 | 0.842 | 0.169 | 0.068 | 0.004 | 2.013e+04 | 25.500 | ok |
| C | within | 27000 | 54 | 54 | 3 | 0.618 | 0.755 | 0.382 | 0.382 | 2.618e+04 | NA | ok |
| C | between | 54 | 54 | 54 | 3 | 0.995 | 0.499 | 0.355 | 0.301 | 2.618e+04 | 90.858 | ok |
| O | within | 27000 | 54 | 54 | 3 | 0.053 | 0.458 | 0.003 | 0.003 | 2.621e+04 | NA | ok |
| O | between | 54 | 54 | 54 | 3 | 0.994 | 0.415 | 0.173 | -0.095 | 2.621e+04 | 1650.184 | ok |
| HA | within | 29000 | 58 | 58 | 3 | 0.160 | 0.258 | 0.026 | 0.023 | 2.356e+04 | NA | ok |
| HA | between | 58 | 58 | 58 | 3 | 0.625 | 0.419 | 0.069 | 0.014 | 2.356e+04 | 3.968 | ok |

## Step 1 Implied Readout

Arithmetic from the prior single-gamma WA run: `q_CO = gamma_lit * 2.41`, `q_CN = gamma_lit * -5.42`.

| stratum | axis | gamma_lit | implied CO | CO spread | implied CN | CN spread |
|---|---:|---:|---:|---|---:|---|
| HN | within | -0.349 +/- 0.311 | -0.842 +/- 0.749 | sign-flip; nearest WA | 1.894 +/- 1.684 | sign-flip; nearest ApSimon |
| HN | between | -1.153 +/- 0.797 | -2.780 +/- 1.921 | sign-flip; nearest WA | 6.251 +/- 4.319 | sign-flip; nearest ApSimon |
| N | within | 0.277 +/- 0.031 | 0.668 +/- 0.075 | below spread; nearest WA | -1.503 +/- 0.169 | above spread; nearest ApSimon |
| N | between | 2.604 +/- 0.965 | 6.276 +/- 2.325 | inside spread; nearest Abraham aliphatic amide | -14.114 +/- 5.229 | inside spread; nearest Abraham aliphatic amide |
| CA | within | 1.357 +/- 0.556 | 3.270 +/- 1.339 | inside spread; nearest WA | -7.353 +/- 3.012 | inside spread; nearest Schneider |
| CA | between | 2.903 +/- 0.637 | 6.997 +/- 1.535 | inside spread; nearest Abraham aliphatic amide | -15.736 +/- 3.451 | below spread; nearest Abraham aliphatic amide |
| C | within | -2.137 +/- 0.023 | -5.150 +/- 0.054 | sign-flip; nearest WA | 11.582 +/- 0.122 | sign-flip; nearest ApSimon |
| C | between | -1.916 +/- 0.411 | -4.619 +/- 0.990 | sign-flip; nearest WA | 10.387 +/- 2.227 | sign-flip; nearest ApSimon |
| O | within | 2.493 +/- 0.229 | 6.008 +/- 0.553 | inside spread; nearest Abraham aliphatic amide | -13.511 +/- 1.243 | inside spread; nearest Abraham aliphatic amide |
| O | between | 38.025 +/- 30.376 | 91.641 +/- 73.206 | above spread; nearest Schneider | -206.097 +/- 164.638 | below spread; nearest Abraham aliphatic amide |
| HA | within | 0.348 +/- 0.252 | 0.838 +/- 0.607 | below spread; nearest WA | -1.885 +/- 1.366 | above spread; nearest ApSimon |
| HA | between | 1.826 +/- 0.379 | 4.400 +/- 0.914 | inside spread; nearest Abraham aliphatic amide | -9.896 +/- 2.055 | inside spread; nearest Schneider |

Step-1 disagreement flags: N between implies CO about +6.3 and CN about -14.1, closest to Abraham aliphatic amide. C within/between sign-flips both categories (CO about -5, CN about +10 to +12), consistent with the carbonyl-C near-field/self-bond warning rather than a transferable far-field Delta-chi. O between is far above the literature spread and weakly predictive out-of-sample. HN/HA within are low-scale diagnostics.

## Literature Spread Used

| category | WA | Abraham aliphatic amide | ApSimon | Schneider |
|---|---:|---:|---:|---:|
| PeptideCO | 2.41 | 6.34 | 12.65 | 14.45 |
| PeptideCN | -5.42 | -14.25 | -3.61 | -7.23 |
| SidechainCO | 2.41 | 6.34 | 12.65 | 14.45 |

## Self Audit

- Aggregated source-sum vs emitted aggregate `mc_lit_T2_local_*`: RMS=1.81596365e-08, max_abs=1.04399998e-07. This checks the category source-sum is the emitted McConnell tensor, not a rebuilt one.
- Aromatic category excluded: rows=1753541, emitted abs T2 sum=0.00000000. RING owns the pi current.
- Manual conversion example: N between PeptideCO beta=-20.81325721; q=-beta/0.55351302=37.60210937; reported q=37.60210937; abs_diff=0.00000000.
- No ORCA, no C++ change, no re-emit, no `trajectory.h5`, and no coordinate/tensor physics recompute in Python.

## Rerun

`python3 src/rediscover/analysis/mcconnell_dchi_calibration.py --out-dir <new-broad-emitted-dir> --decirc-csv <new-single-gamma-csv>`
