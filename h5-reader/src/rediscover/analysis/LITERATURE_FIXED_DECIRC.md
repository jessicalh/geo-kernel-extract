# Literature-Coefficient-Fixed De-Circularisation

Read-only Python analysis of the already-emitted ring-current and McConnell substrate.
Substrate: `/tmp/rdc-ring-mc-capstone`
CSV artifact: `/tmp/rediscover-literature-fixed-decirc/literature_fixed_decirc.csv`

Unit handling: emitted fixed kernels and DFT targets are compared in matched ppm units after within-atom centering. The reported `gamma` is therefore a dimensionless post-hoc multiplier. Correlation columns are un-fitted; `gamma_scaled_R2` is only the scale diagnostic.

Thin-stratum handling: `atoms_total` is the emitted home-stratum atom count. `atom_signal_neff` is a participation-ratio effective atom count from emitted fixed-kernel modulation, and `atom_signal_top90_count` is the number of atoms carrying 90% of that modulation.

Parameterized for the 750-frame substrate: rerun with `--out-dir` pointing at the later emitted ring/mc directory; no row count is hard-coded.

## Results

| mechanism | target | protocol | atoms_total | atom_signal_neff | atom_signal_top90_count | frame_effective_n_sum | target_median_lag1 | pearson_r | pearson_R2 | component_r | absT2_r | gamma | gamma_jackknife_se | gamma_scaled_R2 | verdict_bucket |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ring_current | T0 | frame_split | 41 | 3.2646 | 7 | 1.5412e+04 | 0.1130 | -0.6539 | 0.4275 | NA | NA | -11.3509 | 3.6314 | 0.4272 | form-recovered-scale-fitted (gamma not approx 1) |
| ring_current | T0 | leave_atoms_out | 41 | 3.2646 | 7 | 1.5412e+04 | 0.1130 | -0.6504 | 0.4231 | NA | NA | -11.3042 | 3.6476 | 0.3546 | form-recovered-scale-fitted (gamma not approx 1) |
| ring_current | T2 | frame_split | 41 | 2.0496 | 7 | 1.1927e+04 | 0.2414 | NA | NA | 0.0397 | 0.4704 | 0.8676 | 0.6141 | 0.0009 | zero-circularity recovered law (gamma approx 1) |
| ring_current | T2 | leave_atoms_out | 41 | 2.0496 | 7 | 1.1927e+04 | 0.2414 | NA | NA | 0.0382 | 0.4615 | 0.8901 | 0.5202 | 0.0006 | zero-circularity recovered law (gamma approx 1) |
| mcconnell | T0 | frame_split | 52 | 43.8006 | 42 | 2.0788e+04 | 0.0920 | -0.4040 | 0.1632 | NA | NA | -4.7778 | 0.4034 | 0.1631 | form-recovered-scale-fitted (gamma not approx 1) |
| mcconnell | T0 | leave_atoms_out | 52 | 43.8006 | 42 | 2.0788e+04 | 0.0920 | -0.4085 | 0.1669 | NA | NA | -4.7514 | 0.4094 | 0.1644 | form-recovered-scale-fitted (gamma not approx 1) |
| mcconnell | T2 | frame_split | 52 | 47.1025 | 44 | 1.7014e+04 | 0.1840 | NA | NA | 0.0054 | 0.0881 | -0.0017 | 0.0544 | -0.0002 | form-recovered-scale-fitted (gamma not approx 1) |
| mcconnell | T2 | leave_atoms_out | 52 | 47.1025 | 44 | 1.7014e+04 | 0.1840 | NA | NA | 0.0014 | 0.0904 | 0.0092 | 0.0527 | -0.0001 | form-recovered-scale-fitted (gamma not approx 1) |

## Self Audit

- No fitting is used for `pearson_r`, `pearson_R2`, `component_r`, or `absT2_r`; those correlate the fixed emitted kernel directly against DFT on the scoring rows.
- `gamma` is computed separately as `(K dot DFT) / (K dot K)` on centered emitted arrays. The leave-atoms-out `gamma_scaled_R2` refits gamma on all other atoms for each held-out atom.
- T2 sidecars and targets are mapped with the frozen `change_of_basis.get_C()` matrix; no basis is re-derived in this script.
- Manual gamma check ring_current T0: numerator=-151.49455245, denominator=13.40165595, manual=-11.30416667, function=-11.30416667, abs_diff=0.00000000.
- Manual gamma check mcconnell T0: numerator=-520.90829503, denominator=109.63270302, manual=-4.75139516, function=-4.75139516, abs_diff=0.00000000.

## Verdict Buckets

The bucket rule is mechanical and deliberately reserved: if `gamma = 1` lies within about two delete-atom jackknife SEs, the row is tagged `zero-circularity recovered law`; otherwise it is tagged `form-recovered-scale-fitted`. The scientific call remains for the lead.
