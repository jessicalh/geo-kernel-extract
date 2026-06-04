# Johnson-Bovey Region Recovery

> **Evidence archive — not current truth (trued 2026-06-04).** Preserve as
> analysis provenance; check `../NOW.md` and corrected `../STATE.md` before
> quoting quantitative claims.

Substrate: `/tmp/rediscover-jb-parity-v2/composed`
CSV artifact: `/tmp/rediscover-jb-region-recovery/johnson_bovey_region_recovery.csv`

Read-only analysis of the aromatic ring-facing H stratum. Point-dipole rows use the emitted scalar `sum_dipolar` source column; J-B rows use the emitted source-level `jb_T0` / `jb_T2_local_*` columns computed in C++ from the Johnson-Bovey double-loop with literature G-P intensities.

Geometry regions are parameterized and source-level: `all_valid` excludes self/bonded rings; `in_plane_band` is `|cos_theta| <= 0.2`; `near_field_r_lt4` is `r < 4 A`; `axial_far` is `r >= 4 A` and `|cos_theta| >= 0.7`. The H5 `ring_in_plane_angle` is retained as emitted geometry but is an in-plane azimuth, so region membership uses `cos_theta` and `r`.

All scores are within-atom centered. SEs are delete-atom jackknife within this one protein. `loao_gamma_scaled_R2` fits only a scalar gamma on all other atoms before scoring each held-out atom; direct correlation columns are unfitted.

## Results

| region | model | target | atoms_total | atom_signal_neff | source_rows | pearson_r | pearson_R2 | pearson_R2_jackknife_se | component_r | component_r_jackknife_se | absT2_r | absT2_r_jackknife_se | gamma | gamma_jackknife_se | loao_gamma_scaled_R2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_valid | point_dipole | T0 | 41 | 3.4902 | 110500 | 0.6493 | 0.4215 | 0.0761 | NA | NA | NA | NA | 18.0779 | 5.5618 | 0.3587 |
| all_valid | johnson_bovey | T0 | 41 | 2.2431 | 110500 | 0.6088 | 0.3706 | 0.0684 | NA | NA | NA | NA | 0.7722 | 0.4732 | 0.1740 |
| all_valid | johnson_bovey | T2 | 41 | 1.6106 | 110500 | NA | NA | NA | 0.4731 | 0.0403 | 0.5668 | 0.0286 | 1.0462 | 1.1322 | -0.1092 |
| in_plane_band | point_dipole | T0 | 41 | 10.2875 | 18368 | 0.1401 | 0.0196 | 0.0070 | NA | NA | NA | NA | 25.1525 | 7.2191 | 0.0165 |
| in_plane_band | johnson_bovey | T0 | 41 | 6.3808 | 18368 | 0.1127 | 0.0127 | 0.0046 | NA | NA | NA | NA | 0.9035 | 0.3325 | 0.0095 |
| in_plane_band | johnson_bovey | T2 | 41 | 6.7119 | 18368 | NA | NA | NA | 0.0095 | 0.0111 | -0.1257 | 0.0613 | 0.2019 | 0.2441 | -0.0002 |
| near_field_r_lt4 | point_dipole | T0 | 41 | 2.7041 | 3181 | 0.5912 | 0.3496 | 0.0751 | NA | NA | NA | NA | 16.4623 | 6.8359 | 0.2593 |
| near_field_r_lt4 | johnson_bovey | T0 | 41 | 1.8228 | 3181 | 0.5554 | 0.3084 | 0.0633 | NA | NA | NA | NA | 0.6827 | 0.5718 | 0.0201 |
| near_field_r_lt4 | johnson_bovey | T2 | 41 | 1.4472 | 3181 | NA | NA | NA | 0.4374 | 0.0396 | 0.5689 | 0.0301 | 0.9879 | 1.2645 | -0.1930 |
| axial_far | point_dipole | T0 | 41 | 11.7398 | 42397 | 0.0863 | 0.0074 | 0.0256 | NA | NA | NA | NA | 9.5759 | 11.9440 | -0.0140 |
| axial_far | johnson_bovey | T0 | 41 | 6.3111 | 42397 | 0.0294 | 0.0009 | 0.0218 | NA | NA | NA | NA | 0.1574 | 0.7946 | -0.0364 |
| axial_far | johnson_bovey | T2 | 41 | 10.6763 | 42397 | NA | NA | NA | 0.1177 | 0.0369 | 0.1199 | 0.0570 | 1.2955 | 0.5249 | 0.0096 |

## Interpretation

- In-plane band: point-dipole T0 R2=0.0196 +/- 0.0070; fixed J-B T0 R2=0.0127 +/- 0.0046; fixed J-B T2 component r=0.0095 +/- 0.0111. This does not show J-B recovering the in-plane DFT modulation that the scalar point-dipole misses.
- r < 4 A near-field diagnostic: point-dipole T0 R2=0.3496 +/- 0.0751; fixed J-B T0 R2=0.3084 +/- 0.0633; fixed J-B T2 component r=0.4374 +/- 0.0396, |T2| r=0.5689 +/- 0.0301. This region is thin, with atom-signal N_eff=1.4472.
- The all-valid rows still carry the main ring signal, especially J-B T2 magnitude, but the in-plane source partition alone is near-null against DFT in this single-protein substrate. Verdict remains reserved for the lead.

## Self Audit

- Source-level emitted unit-current J-B T0 valid-sum vs emitted H5 `bare_T0`: r=0.9825, RMS=0.0111 ppm_T_per_nA, max_abs=0.0354 ppm_T_per_nA.
- Manual emitted-row sum check for one in-plane atom/frame: atom=72, h5_row=0, source_rows=1, manual_jb_T0=-0.00852413, grouped_jb_T0=-0.00852413, abs_diff=0.00000000.
- Python did not open trajectory.h5, evaluate Biot-Savart/J-B, compute `(3cos^2-1)/r^3`, or project tensors; it summed emitted source columns and read emitted target sidecars.
- Point-dipole T2 is not reported because `sum_dipolar` is a scalar T0 proxy, not an emitted tensor law.

## De-Circularised Read

The de-circularised test is the direct J-B correlation/gamma rows: the kernel uses fixed literature G-P intensities in C++ and no fitted scale for the reported correlation. Gamma and LOAO scaled R2 are diagnostics, not the verdict. Scientific verdict remains reserved for the lead.
