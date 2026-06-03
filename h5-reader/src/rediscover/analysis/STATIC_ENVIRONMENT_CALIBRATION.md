> SUPERSEDED / pre-units-audit (2026-06-02): see `../UNITS_AND_ISSUES_AUDIT.md` + `../STATE.md`. The `mcconnell_T2_fixed` and `ring_current_T2_fixed` rows carry `calibration_flag = calibrated-to-physics`, `literature_value = 1.0000`, "fixed emitted literature-kernel multiplier" — this is FALSE: the script labels bare emitted kernels as fixed-literature with coefficient 1 (`static_environment_calibration.py:493-523`), but the emit provides BARE kernels (ring unit-current `ppm_T_per_nA`, McConnell unit-Δχ `Angstrom^-3`), not literature-ppm predictions. So those γ are unitful fitted multipliers (e.g. McConnell-N γ≈26, McConnell-HN γ=−11.0, ring-N γ=−391), NOT dimensionless-≈1 de-circularised checks; "calibrated-to-physics" is mislabeled. Treat the verdict bucket as "form-recovered, scale-fitted" (per STATE). `between_LOAO_absT2_r` correlations are scale-invariant and stand.

# Static Environment Calibration

Between-axis recovery from per-atom means. Python reads only emitted CSV/NPY substrate,
uses the frozen `change_of_basis.get_C()` T2 map, and reports correlations rather than
claiming exact agreement. Verdict reserved.

## Calibration Table

| mechanism | target | stratum | n_atoms_between | primary_coef_name | primary_coef | primary_se | secondary_coef_name | secondary_coef | secondary_se | coefficient_units | calibration_metric_name | calibration_metric | between_LOAO_R2 | between_LOAO_absT2_r | literature_value | literature_units | recovered_literature_ratio | calibration_flag |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| apbs_charge_EFG_T2 | T2 | N | 54 | gamma | -15.2073 | 2.1234 |  | NA | NA | dimensionless fitted multiplier of emitted APBS EFG T2 | between_LOAO_absT2_r | 0.5172 | -18.1916 | 0.5172 | NA |  | NA | empirical/form-only |
| apbs_charge_EFG_T2 | T2 | CA | 54 | gamma | -1.0776 | 0.6683 |  | NA | NA | dimensionless fitted multiplier of emitted APBS EFG T2 | between_LOAO_absT2_r | 0.1380 | -3.0223 | 0.1380 | NA |  | NA | empirical/form-only |
| apbs_charge_EFG_T2 | T2 | C | 54 | gamma | -22.5717 | 1.2735 |  | NA | NA | dimensionless fitted multiplier of emitted APBS EFG T2 | between_LOAO_absT2_r | -0.4678 | -68.5026 | -0.4678 | NA |  | NA | empirical/form-only |
| apbs_charge_EFG_T2 | T2 | O | 54 | gamma | 174.2099 | 4.9138 |  | NA | NA | dimensionless fitted multiplier of emitted APBS EFG T2 | between_LOAO_absT2_r | -0.0455 | -71.7750 | -0.0455 | NA |  | NA | empirical/form-only |
| apbs_charge_EFG_T2 | T2 | HN | 52 | gamma | 2.0060 | 0.3138 |  | NA | NA | dimensionless fitted multiplier of emitted APBS EFG T2 | between_LOAO_absT2_r | 0.0855 | -3.0374 | 0.0855 | NA |  | NA | empirical/form-only |
| apbs_charge_EFG_T2 | T2 | HA | 58 | gamma | 1.9715 | 0.4704 |  | NA | NA | dimensionless fitted multiplier of emitted APBS EFG T2 | between_LOAO_absT2_r | -0.1036 | -0.5803 | -0.1036 | NA |  | NA | empirical/form-only |
| buckingham_efield_sigma | sigma_iso | N | 54 | A_Eproj | -100.7241 | 243.6616 | B_Emag_sq | 42.8825 | 138.0670 | ppm per emitted electric-field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.9020 | -0.9020 | NA | NA |  | NA | empirical/form-only |
| buckingham_efield_sigma | sigma_iso | CA | 54 | A_Eproj | -7.5810 | 5.9361 | B_Emag_sq | 8.7948 | 3.8167 | ppm per emitted electric-field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.0367 | -0.0367 | NA | NA |  | NA | empirical/form-only |
| buckingham_efield_sigma | sigma_iso | C | 54 | A_Eproj | -0.7191 | 9.5340 | B_Emag_sq | 0.6377 | 4.3584 | ppm per emitted electric-field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.4262 | -0.4262 | NA | NA |  | NA | empirical/form-only |
| buckingham_efield_sigma | sigma_iso | O | 54 | A_Eproj | 263.4324 | 86.6885 | B_Emag_sq | 116.1675 | 41.4251 | ppm per emitted electric-field unit; ppm per emitted field-unit squared | between_LOAO_R2 | 0.0911 | 0.0911 | NA | NA |  | NA | empirical/form-only |
| buckingham_efield_sigma | sigma_iso | HN | 52 | A_Eproj | -15.6245 | 2.8050 | B_Emag_sq | -9.3380 | 2.0473 | ppm per emitted electric-field unit; ppm per emitted field-unit squared | between_LOAO_R2 | 0.3972 | 0.3972 | NA | NA |  | NA | empirical/form-only |
| buckingham_efield_sigma | sigma_iso | HA | 58 | A_Eproj | 0.0170 | 0.8436 | B_Emag_sq | -9.2918e-04 | 1.9538 | ppm per emitted electric-field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.3214 | -0.3214 | NA | NA |  | NA | empirical/form-only |
| ff14sb_charge_EFG_T2 | T2 | N | 54 | gamma | 9.7006 | 0.2274 |  | NA | NA | dimensionless fitted multiplier of emitted charge EFG-like T2 | between_LOAO_absT2_r | 0.7764 | -10.2360 | 0.7764 | NA |  | NA | empirical/form-only |
| ff14sb_charge_EFG_T2 | T2 | CA | 54 | gamma | 5.2530 | 0.5692 |  | NA | NA | dimensionless fitted multiplier of emitted charge EFG-like T2 | between_LOAO_absT2_r | -0.0357 | -2.8568 | -0.0357 | NA |  | NA | empirical/form-only |
| ff14sb_charge_EFG_T2 | T2 | C | 54 | gamma | 10.2167 | 0.2915 |  | NA | NA | dimensionless fitted multiplier of emitted charge EFG-like T2 | between_LOAO_absT2_r | 0.1499 | -54.7388 | 0.1499 | NA |  | NA | empirical/form-only |
| ff14sb_charge_EFG_T2 | T2 | O | 54 | gamma | 157.8795 | 13.4092 |  | NA | NA | dimensionless fitted multiplier of emitted charge EFG-like T2 | between_LOAO_absT2_r | 0.2441 | -94.8070 | 0.2441 | NA |  | NA | empirical/form-only |
| ff14sb_charge_EFG_T2 | T2 | HN | 52 | gamma | 0.8029 | 0.1894 |  | NA | NA | dimensionless fitted multiplier of emitted charge EFG-like T2 | between_LOAO_absT2_r | -0.3547 | -3.1980 | -0.3547 | NA |  | NA | empirical/form-only |
| ff14sb_charge_EFG_T2 | T2 | HA | 58 | gamma | 1.8190 | 0.2487 |  | NA | NA | dimensionless fitted multiplier of emitted charge EFG-like T2 | between_LOAO_absT2_r | -0.3384 | -0.4877 | -0.3384 | NA |  | NA | empirical/form-only |
| ff14sb_charge_field_sigma | sigma_iso | N | 54 | A_field_z | -332.7077 | 349.6537 | B_field_mag_sq | 201.9640 | 295.3977 | ppm per emitted FF14SB field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.2135 | -0.2135 | NA | NA |  | NA | empirical/form-only |
| ff14sb_charge_field_sigma | sigma_iso | CA | 54 | A_field_z | 7.2267 | 34.5612 | B_field_mag_sq | -55.9875 | 231.0676 | ppm per emitted FF14SB field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.1062 | -0.1062 | NA | NA |  | NA | empirical/form-only |
| ff14sb_charge_field_sigma | sigma_iso | C | 54 | A_field_z | 2.7095 | 14.3506 | B_field_mag_sq | 23.6332 | 28.5962 | ppm per emitted FF14SB field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.0981 | -0.0981 | NA | NA |  | NA | empirical/form-only |
| ff14sb_charge_field_sigma | sigma_iso | O | 54 | A_field_z | -91.2127 | 267.5936 | B_field_mag_sq | 1466.8372 | 1494.7917 | ppm per emitted FF14SB field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.0815 | -0.0815 | NA | NA |  | NA | empirical/form-only |
| ff14sb_charge_field_sigma | sigma_iso | HN | 52 | A_field_z | -22.4314 | 6.4192 | B_field_mag_sq | 38.8038 | 27.3706 | ppm per emitted FF14SB field unit; ppm per emitted field-unit squared | between_LOAO_R2 | 0.4834 | 0.4834 | NA | NA |  | NA | empirical/form-only |
| ff14sb_charge_field_sigma | sigma_iso | HA | 58 | A_field_z | -6.3120 | 3.9215 | B_field_mag_sq | -27.8795 | 50.0176 | ppm per emitted FF14SB field unit; ppm per emitted field-unit squared | between_LOAO_R2 | -0.0179 | -0.0179 | NA | NA |  | NA | empirical/form-only |
| mcconnell_T2_fixed | T2 | N | 54 | gamma | 25.9586 | 5.7059 |  | NA | NA | dimensionless multiplier of emitted fixed McConnell T2 | between_LOAO_absT2_r | 0.4428 | -18.3383 | 0.4428 | 1.0000 | fixed emitted literature-kernel multiplier | 25.9586 | calibrated-to-physics |
| mcconnell_T2_fixed | T2 | CA | 54 | gamma | 18.7569 | 3.8247 |  | NA | NA | dimensionless multiplier of emitted fixed McConnell T2 | between_LOAO_absT2_r | -0.0942 | -2.8350 | -0.0942 | 1.0000 | fixed emitted literature-kernel multiplier | 18.7569 | calibrated-to-physics |
| mcconnell_T2_fixed | T2 | C | 54 | gamma | 121.2482 | 2.4535 |  | NA | NA | dimensionless multiplier of emitted fixed McConnell T2 | between_LOAO_absT2_r | 0.3069 | -39.5032 | 0.3069 | 1.0000 | fixed emitted literature-kernel multiplier | 121.2482 | calibrated-to-physics |
| mcconnell_T2_fixed | T2 | O | 54 | gamma | 419.0313 | 31.3223 |  | NA | NA | dimensionless multiplier of emitted fixed McConnell T2 | between_LOAO_absT2_r | 0.1275 | -95.0805 | 0.1275 | 1.0000 | fixed emitted literature-kernel multiplier | 419.0313 | calibrated-to-physics |
| mcconnell_T2_fixed | T2 | HN | 52 | gamma | -11.0278 | 0.4753 |  | NA | NA | dimensionless multiplier of emitted fixed McConnell T2 | between_LOAO_absT2_r | 0.2029 | -1.3353 | 0.2029 | 1.0000 | fixed emitted literature-kernel multiplier | -11.0278 | calibrated-to-physics |
| mcconnell_T2_fixed | T2 | HA | 58 | gamma | -4.3548 | 0.3904 |  | NA | NA | dimensionless multiplier of emitted fixed McConnell T2 | between_LOAO_absT2_r | 0.0530 | -0.1492 | 0.0530 | 1.0000 | fixed emitted literature-kernel multiplier | -4.3548 | calibrated-to-physics |
| mcconnell_scalar_bare | sigma_iso | N | 54 | k_bare | -12.6987 | 26.5053 |  | NA | NA | ppm A^3 per emitted bare bond-anisotropy sum | between_LOAO_R2 | -2.7744 | -2.7744 | NA | NA |  | NA | empirical/form-only |
| mcconnell_scalar_bare | sigma_iso | CA | 54 | k_bare | -5.1812 | 14.9502 |  | NA | NA | ppm A^3 per emitted bare bond-anisotropy sum | between_LOAO_R2 | -0.0626 | -0.0626 | NA | NA |  | NA | empirical/form-only |
| mcconnell_scalar_bare | sigma_iso | C | 54 | k_bare | 0.3919 | 2.1406 |  | NA | NA | ppm A^3 per emitted bare bond-anisotropy sum | between_LOAO_R2 | -1.8538 | -1.8538 | NA | NA |  | NA | empirical/form-only |
| mcconnell_scalar_bare | sigma_iso | O | 54 | k_bare | -16.7395 | 61.0746 |  | NA | NA | ppm A^3 per emitted bare bond-anisotropy sum | between_LOAO_R2 | -0.4984 | -0.4984 | NA | NA |  | NA | empirical/form-only |
| mcconnell_scalar_bare | sigma_iso | HN | 52 | k_bare | -4.7744 | 1.6485 |  | NA | NA | ppm A^3 per emitted bare bond-anisotropy sum | between_LOAO_R2 | 0.1308 | 0.1308 | NA | NA |  | NA | empirical/form-only |
| mcconnell_scalar_bare | sigma_iso | HA | 58 | k_bare | 1.3103 | 0.8050 |  | NA | NA | ppm A^3 per emitted bare bond-anisotropy sum | between_LOAO_R2 | -0.0229 | -0.0229 | NA | NA |  | NA | empirical/form-only |
| ring_current_T2_fixed | T2 | N | 54 | gamma | -391.7037 | 110.6952 |  | NA | NA | dimensionless multiplier of emitted fixed ring-current T2 | between_LOAO_absT2_r | 0.2213 | -18.0393 | 0.2213 | 1.0000 | fixed emitted literature-kernel multiplier | -391.7037 | calibrated-to-physics |
| ring_current_T2_fixed | T2 | CA | 54 | gamma | -6.5543 | 19.4399 |  | NA | NA | dimensionless multiplier of emitted fixed ring-current T2 | between_LOAO_absT2_r | -0.1324 | -3.0370 | -0.1324 | 1.0000 | fixed emitted literature-kernel multiplier | -6.5543 | calibrated-to-physics |
| ring_current_T2_fixed | T2 | C | 54 | gamma | -263.8020 | 98.9507 |  | NA | NA | dimensionless multiplier of emitted fixed ring-current T2 | between_LOAO_absT2_r | -0.0778 | -76.9954 | -0.0778 | 1.0000 | fixed emitted literature-kernel multiplier | -263.8020 | calibrated-to-physics |
| ring_current_T2_fixed | T2 | O | 54 | gamma | -1569.4129 | 333.4185 |  | NA | NA | dimensionless multiplier of emitted fixed ring-current T2 | between_LOAO_absT2_r | -0.1550 | -107.6639 | -0.1550 | 1.0000 | fixed emitted literature-kernel multiplier | -1569.4129 | calibrated-to-physics |
| ring_current_T2_fixed | T2 | HN | 52 | gamma | 18.3875 | 10.1299 |  | NA | NA | dimensionless multiplier of emitted fixed ring-current T2 | between_LOAO_absT2_r | -0.1537 | -3.3119 | -0.1537 | 1.0000 | fixed emitted literature-kernel multiplier | 18.3875 | calibrated-to-physics |
| ring_current_T2_fixed | T2 | HA | 58 | gamma | -19.5352 | 4.4818 |  | NA | NA | dimensionless multiplier of emitted fixed ring-current T2 | between_LOAO_absT2_r | 0.2277 | -0.5035 | 0.2277 | 1.0000 | fixed emitted literature-kernel multiplier | -19.5352 | calibrated-to-physics |
| ring_current_scalar_bare | sigma_iso | N | 54 | k_bare | 34.6589 | 278.2302 |  | NA | NA | ppm A^3 per emitted bare ring sum | between_LOAO_R2 | -2.2320 | -2.2320 | NA | 21.0000 | ppm A^3 for the canonical Pople-style ring-current control | 1.6504 | empirical/form-only |
| ring_current_scalar_bare | sigma_iso | CA | 54 | k_bare | 14.2576 | 37.9824 |  | NA | NA | ppm A^3 per emitted bare ring sum | between_LOAO_R2 | -0.3291 | -0.3291 | NA | 21.0000 | ppm A^3 for the canonical Pople-style ring-current control | 0.6789 | empirical/form-only |
| ring_current_scalar_bare | sigma_iso | C | 54 | k_bare | 28.7167 | 22.1680 |  | NA | NA | ppm A^3 per emitted bare ring sum | between_LOAO_R2 | -0.0396 | -0.0396 | NA | 21.0000 | ppm A^3 for the canonical Pople-style ring-current control | 1.3675 | empirical/form-only |
| ring_current_scalar_bare | sigma_iso | O | 54 | k_bare | 181.1786 | 222.9615 |  | NA | NA | ppm A^3 per emitted bare ring sum | between_LOAO_R2 | -0.1188 | -0.1188 | NA | 21.0000 | ppm A^3 for the canonical Pople-style ring-current control | 8.6276 | empirical/form-only |
| ring_current_scalar_bare | sigma_iso | HN | 52 | k_bare | 2.7096 | 11.5176 |  | NA | NA | ppm A^3 per emitted bare ring sum | between_LOAO_R2 | -0.0901 | -0.0901 | NA | 21.0000 | ppm A^3 for the canonical Pople-style ring-current control | 0.1290 | empirical/form-only |
| ring_current_scalar_bare | sigma_iso | HA | 58 | k_bare | 12.2267 | 3.8638 |  | NA | NA | ppm A^3 per emitted bare ring sum | between_LOAO_R2 | 0.1088 | 0.1088 | NA | 21.0000 | ppm A^3 for the canonical Pople-style ring-current control | 0.5822 | empirical/form-only |

## Literature And De-Circularising Notes

- Buckingham A is a real shielding-polarizability coefficient. The sourced literature supports the
  expansion and amide dipole-shielding polarizabilities, but this emitted APBS scalar field does
  not give a defensible stratum/unit-matched scalar A to plug without inventing a number.
- Ring-current scalar rows are bare broad-substrate sums. The canonical ring-current control
  recovered k near 21 ppm A^3 elsewhere; the fixed T2 rows are the direct sidecar comparison here.
- McConnell scalar rows are untyped broad bond sums. The fixed T2 rows compare to the emitted
  McConnell-style sidecar, with the producer-branch caveat in the notes column.
- Charge field and charge EFG rows are form-only empirical calibrations; no clean universal
  coefficient was found for these emitted sidecars.

## Self Audit

- Manual OLS re-derivation: buckingham_efield_sigma N primary=-100.7241, script=-100.7241, abs_diff=3.041123e-11.
- DOI audit completed with live DOI/Crossref resolution for all DOI strings cited in this report.
- Grep audit passed for H5 access and Python-side physics reconstruction patterns.

## Source Notes

- Buckingham 1960: DOI https://doi.org/10.1139/v60-040.
- Boyd et al. 2003: DOI https://doi.org/10.1021/ja034855y.
- Pople 1956: DOI https://doi.org/10.1063/1.1742701.
- McConnell 1957: DOI https://doi.org/10.1063/1.1743676.
