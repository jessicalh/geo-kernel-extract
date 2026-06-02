# Ring Literature-Scaled De-Circularisation

Substrate: `/tmp/rediscover-jb-parity-v2/composed`
CSV artifact: `/tmp/rediscover-ring-literature-decirc/ring_literature_decirc.csv`
Audit artifact: `/tmp/rediscover-ring-literature-decirc/ring_literature_decirc_audit.json`

Read-only Python analysis of the aromatic ring-facing H stratum. Inputs are the emitted `ring_current_sources.csv` columns `jb_T*` and `jb_unit_T*`, plus the emitted local DFT T2 sidecar. No C++ was changed, no data was re-emitted, and no H5 file is read.

Axes: `within` is per-atom de-meaned with no intercept. `between` first averages each atom over frames, then fits an intercept plus one shared gamma across atoms. T2 between uses a five-component intercept vector and one shared gamma; the table reports that intercept vector norm.

Thin rows: `thin_flag=thin` means atom-signal participation N_eff < 5; these rows are diagnostics, not strong per-type claims. Confidence is delete-atom jackknife within this one protein; AR(1) frame N_eff is reported as a dependence diagnostic, not a population sample size.

## Headline

Confirmed unit finding: the prior bare-kernel T0 gamma is recovered as an intensity, not as a failed dimensionless coefficient. The corrected literature-scaled rows are compatible with gamma=1 on the all-valid stratum by delete-atom jackknife, but the per-ring-type rows are thin and must not be oversold.

Per-ring-type lead read, within T0:

| ring_type | atom_signal_Neff | gamma_lit | gamma_bare nA/T | Pople nA/T | verdict |
| --- | --- | --- | --- | --- | --- |
| PHE | 2.4962 | 1.3836 +/- 0.0494 | -16.6033 +/- 0.5924 | -12.00 | form-recovered-scale-fitted / form-recovered-scale-fitted |
| TYR | 2.7558 | 1.1841 +/- 0.1378 | -13.3571 +/- 1.5543 | -11.28 | zero-circularity recovered law / zero-circularity recovered law |
| TRP-6 | 1.0399 | 1.1576 +/- 1.8095 | -14.4468 +/- 22.5831 | -12.48 | zero-circularity recovered law / zero-circularity recovered law |
| TRP-5 | 1.1799 | 1.4934 +/- 2.4216 | -10.0359 +/- 16.2730 | -6.72 | zero-circularity recovered law / zero-circularity recovered law |
| TRP-perimeter | 1.1486 | 0.8533 +/- 0.1650 | -16.3825 +/- 3.1680 | -19.20 | zero-circularity recovered law / zero-circularity recovered law |
| HIS | 3.1526 | 4.1279 +/- 0.3339 | -21.2999 +/- 1.7231 | -5.16 | form-recovered-scale-fitted / form-recovered-scale-fitted |

All per-type within-T0 rows are below the effective-atom threshold, so their fitted scales are diagnostics. The six-membered pooled/all-valid row is the robust unit check because the aromatic H stratum is dominated by six-membered-ring modulation.

All-valid stratum cross-checks:

- T0 within all-valid: gamma_lit=0.7722 +/- 0.4732, gamma_bare=-11.3259 +/- 3.6683 nA/T vs Pople -12.50..-11.00; intercept=NA; lit bucket=zero-circularity recovered law; bare bucket=zero-circularity recovered law.
- T0 between all-valid: gamma_lit=0.6423 +/- 0.5206, gamma_bare=-9.5946 +/- 4.5340 nA/T vs Pople -12.50..-11.00; intercept=24.2516; lit bucket=zero-circularity recovered law; bare bucket=zero-circularity recovered law.
- T2 within all-valid: gamma_lit=1.0462 +/- 1.1322, gamma_bare=-15.9796 +/- 10.5704 nA/T vs Pople -12.50..-11.00; intercept=NA; lit bucket=zero-circularity recovered law; bare bucket=zero-circularity recovered law.
- T2 between all-valid: gamma_lit=0.9230 +/- 0.5306, gamma_bare=-13.6405 +/- 4.6707 nA/T vs Pople -12.50..-11.00; intercept=2.3974; lit bucket=zero-circularity recovered law; bare bucket=zero-circularity recovered law.

Plain comparison: all-valid within T0 gamma_bare=-11.33 nA/T lands in the Pople six-membered range (-11 to -12.5). TYR and the TRP perimeter are compatible with their table intensities; PHE and HIS within rows are form-recovered-scale-fitted; TRP-5 is too thin to claim even though the broad SE overlaps its pyrrole value.

## Results

| ring_type | axis | target | atoms_active | atom_signal_neff | ar1_frame_neff | thin_flag | gamma_lit | gamma_bare_nA_per_T | Pople_nA_per_T | intercept | unfitted_metric | LOAO_R2_lit | lit_bucket | bare_bucket |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_valid | within | T0 | 41 | 2.2431 | 1.0215e+04 | thin | 0.7722 +/- 0.4732 | -11.3259 +/- 3.6683 | -12.50..-11.00 | NA | r=0.6088 | 0.1740 | zero-circularity recovered law | zero-circularity recovered law |
| all_valid | between | T0 | 41 | 2.2431 | 1.0215e+04 | thin | 0.6423 +/- 0.5206 | -9.5946 +/- 4.5340 | -12.50..-11.00 | 24.2516 | r=0.7190 | 0.0725 | zero-circularity recovered law | zero-circularity recovered law |
| all_valid | within | T2 | 41 | 1.6106 | 8584.5270 | thin | 1.0462 +/- 1.1322 | -15.9796 +/- 10.5704 | -12.50..-11.00 | NA | comp_r=0.4731, |T2|r=0.5668 | -0.1092 | zero-circularity recovered law | zero-circularity recovered law |
| all_valid | between | T2 | 41 | 1.6106 | 8584.5270 | thin | 0.9230 +/- 0.5306 | -13.6405 +/- 4.6707 | -12.50..-11.00 | 2.3974 | comp_r=0.2388, |T2|r=0.3632 | 0.0123 | zero-circularity recovered law | zero-circularity recovered law |
| PHE | within | T0 | 12 | 2.4962 | 2178.6993 | thin | 1.3836 +/- 0.0494 | -16.6033 +/- 0.5924 | -12.00 | NA | r=0.4764 | 0.2265 | form-recovered-scale-fitted | form-recovered-scale-fitted |
| PHE | between | T0 | 12 | 2.4962 | 2178.6993 | thin | 1.7658 +/- 0.9235 | -21.1897 +/- 11.0826 | -12.00 | 24.3254 | r=0.2838 | -0.0017 | zero-circularity recovered law | zero-circularity recovered law |
| PHE | within | T2 | 12 | 2.7058 | 1623.2912 | thin | 2.5541 +/- 0.2548 | -30.6495 +/- 3.0573 | -12.00 | NA | comp_r=0.3445, |T2|r=0.4338 | 0.1168 | form-recovered-scale-fitted | form-recovered-scale-fitted |
| PHE | between | T2 | 12 | 2.7058 | 1623.2912 | thin | 2.1218 +/- 1.6991 | -25.4614 +/- 20.3886 | -12.00 | 2.3445 | comp_r=0.0846, |T2|r=0.1840 | -0.0492 | zero-circularity recovered law | zero-circularity recovered law |
| TYR | within | T0 | 41 | 2.7558 | 1.0605e+04 | thin | 1.1841 +/- 0.1378 | -13.3571 +/- 1.5543 | -11.28 | NA | r=0.3032 | 0.0899 | zero-circularity recovered law | zero-circularity recovered law |
| TYR | between | T0 | 41 | 2.7558 | 1.0605e+04 | thin | 1.2229 +/- 0.2211 | -13.7944 +/- 2.4938 | -11.28 | 24.3252 | r=0.4090 | 0.1148 | zero-circularity recovered law | zero-circularity recovered law |
| TYR | within | T2 | 41 | 3.0531 | 8778.0917 | thin | 2.3053 +/- 0.4919 | -26.0040 +/- 5.5490 | -11.28 | NA | comp_r=0.3043, |T2|r=0.2835 | 0.0857 | form-recovered-scale-fitted | form-recovered-scale-fitted |
| TYR | between | T2 | 41 | 3.0531 | 8778.0917 | thin | 1.6186 +/- 1.0040 | -18.2576 +/- 11.3254 | -11.28 | 2.3384 | comp_r=0.1500, |T2|r=0.1243 | -0.0315 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-6 | within | T0 | 19 | 1.0399 | 4709.1427 | thin | 1.1576 +/- 1.8095 | -14.4468 +/- 22.5831 | -12.48 | NA | r=0.3003 | -0.1456 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-6 | between | T0 | 19 | 1.0399 | 4709.1427 | thin | 1.7997 +/- 0.6377 | -22.4597 +/- 7.9583 | -12.48 | 24.3126 | r=0.5813 | 0.2567 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-6 | within | T2 | 19 | 1.0464 | 4286.3549 | thin | 2.0062 +/- 1.6867 | -25.0377 +/- 21.0496 | -12.48 | NA | comp_r=0.2766, |T2|r=0.3636 | 0.0185 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-6 | between | T2 | 19 | 1.0464 | 4286.3549 | thin | 2.4495 +/- 0.1952 | -30.5699 +/- 2.4366 | -12.48 | 2.3249 | comp_r=0.1993, |T2|r=0.2835 | 0.0120 | form-recovered-scale-fitted | form-recovered-scale-fitted |
| TRP-5 | within | T0 | 19 | 1.1799 | 4878.1236 | thin | 1.4934 +/- 2.4216 | -10.0359 +/- 16.2730 | -6.72 | NA | r=0.1154 | -0.0263 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-5 | between | T0 | 19 | 1.1799 | 4878.1236 | thin | 9.0307 +/- 1.5809 | -60.6865 +/- 10.6239 | -6.72 | 24.3184 | r=0.5889 | 0.3002 | form-recovered-scale-fitted | form-recovered-scale-fitted |
| TRP-5 | within | T2 | 19 | 1.1557 | 4466.0238 | thin | 5.2729 +/- 2.1732 | -35.4336 +/- 14.6040 | -6.72 | NA | comp_r=0.1949, |T2|r=0.3310 | 0.0307 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-5 | between | T2 | 19 | 1.1557 | 4466.0238 | thin | 6.1992 +/- 3.8865 | -41.6588 +/- 26.1172 | -6.72 | 2.3349 | comp_r=0.1066, |T2|r=0.3109 | -0.0387 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-perimeter | within | T0 | 19 | 1.1486 | 4785.5860 | thin | 0.8533 +/- 0.1650 | -16.3825 +/- 3.1680 | -19.20 | NA | r=0.3469 | 0.1153 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-perimeter | between | T0 | 19 | 1.1486 | 4785.5860 | thin | 0.8672 +/- 0.1758 | -16.6511 +/- 3.3758 | -19.20 | 24.3132 | r=0.5856 | 0.2932 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-perimeter | within | T2 | 19 | 1.0943 | 4377.6800 | thin | 1.1702 +/- 0.3410 | -22.4677 +/- 6.5464 | -19.20 | NA | comp_r=0.3028, |T2|r=0.3675 | 0.0831 | zero-circularity recovered law | zero-circularity recovered law |
| TRP-perimeter | between | T2 | 19 | 1.0943 | 4377.6800 | thin | 1.1187 +/- 0.4540 | -21.4789 +/- 8.7164 | -19.20 | 2.3352 | comp_r=0.1820, |T2|r=0.2937 | -0.0047 | zero-circularity recovered law | zero-circularity recovered law |
| HIS | within | T0 | 38 | 3.1526 | 7973.4484 | thin | 4.1279 +/- 0.3339 | -21.2999 +/- 1.7231 | -5.16 | NA | r=0.2069 | 0.0423 | form-recovered-scale-fitted | form-recovered-scale-fitted |
| HIS | between | T0 | 38 | 3.1526 | 7973.4484 | thin | 2.5775 +/- 6.2485 | -13.2999 +/- 32.2424 | -5.16 | 24.3709 | r=0.0476 | -0.0743 | zero-circularity recovered law | zero-circularity recovered law |
| HIS | within | T2 | 38 | 4.8401 | 7102.2268 | thin | 9.1337 +/- 2.5702 | -47.1298 +/- 13.2622 | -5.16 | NA | comp_r=0.1825, |T2|r=0.0935 | 0.0288 | form-recovered-scale-fitted | form-recovered-scale-fitted |
| HIS | between | T2 | 38 | 4.8401 | 7102.2268 | thin | 4.6076 +/- 3.9024 | -23.7752 +/- 20.1363 | -5.16 | 2.3549 | comp_r=-0.0513, |T2|r=-0.1707 | -0.0520 | zero-circularity recovered law | zero-circularity recovered law |

## Pople Intensities

| ring_type | source | literature intensity nA/T |
| --- | --- | --- |
| PHE | benzene | -12.00 |
| TYR | phenol | -11.28 |
| TRP-6 | tryptophan 6-ring | -12.48 |
| TRP-5 | tryptophan pyrrole | -6.72 |
| TRP-perimeter | tryptophan indole perimeter | -19.20 |
| HIS | imidazole | -5.16 |

## Verdict Buckets

- `zero-circularity recovered law`: literature-scaled gamma is compatible with 1 within about two delete-atom jackknife SEs, and/or bare gamma is compatible with the relevant Pople intensity or all-valid six-membered dominated range.
- `form-recovered-scale-fitted`: the emitted kernel shape correlates, but the reported scale is not compatible with the fixed literature scale by that mechanical rule.
- Scientific interpretation remains reserved for the lead: this is one protein, with correlated frames, and thin ring-type rows are explicitly flagged.

## Self Audit

- Manual all-valid within T0 gamma from emitted `jb_T0` and DFT: numerator=1942.53662255, denominator=2515.48443653, manual=0.77223162, function=0.77223162, abs_diff=0.00000000.
- One emitted source-row scale check: atom=56, h5_row=0, ring_type=0, jb_T0=-0.0052658600, jb_unit_T0*intensity=-0.0052658600, abs_diff=3.9999999840e-12.
- Grouped source `jb_unit_T0` vs aggregate bare `bare_T0`: r=0.9825, RMS=0.01107287, max_abs=0.03541865 ppm_T_per_nA.
- Unfitted metrics (`r`, `component_r`, `|T2|r`) are computed directly from emitted kernels and DFT targets. Gamma and LOAO R2 are separate scale diagnostics and are not used to produce those correlations.
- The script reads CSV/NPY substrate only. It does not open H5 and does not evaluate any ring-current formula or tensor projection.

## Rerun

For the 750-DFT substrate, rerun with `python3 src/rediscover/analysis/ring_literature_decirc.py --out-dir <new-emitted-dir>`; no row count is hard-coded.
