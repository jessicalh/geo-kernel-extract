# Corrected Backbone Snapshot

Fresh corrected 1P9J broad-backbone snapshot from `/tmp/rediscover-corrected-backbone-snapshot-1p9j`. Python analysis reads emitted CSV/NPY substrate only; no Python H5 read, no ORCA run, no `trajectory.h5` write, and no C++ edit was made in this pass.

## Lead Read

- Ring source law: standalone Pople/JB all-valid T2 between gamma=0.923 +/- 0.531, component_r=0.239, |T2|r=0.363, LOAO_R2=0.012; bucket `recovered law`. The fresh T2-between LOAO value is not the expected 0.62, so this report does not claim that number.
- Ring unit check: all-valid T0 within gamma=0.772 +/- 0.473; bare intensity gamma=-11.326 +/- 3.668 nA/T, matching the six-membered Pople scale within this one-protein jackknife.
- McConnell: valid-source/self-or-bonded filtering is active; alone it still lands in `can't-make-it-work-for-now`, with useful signal mainly as a joint-fit component.
- Charge q/r3: calibrated as `form-recovered-scale-fitted`; it is the strongest broad static T2 mechanism for several strata, especially N and C/within-frame reads.
- Field/Buckingham: scalar field calibrations are uneven; HN is the one clear field-like static read, while most other strata remain in `can't-make-it-work-for-now`.
- EFG local-frame: local-frame audit and nonlinear EFG distill stay near zero or negative R2; bucket `can't-make-it-work-for-now`.
- AIMNet2: useful learnable ceiling only, not a physical law; reported as `form-recovered-scale-fitted` ceiling rows with ridge jackknife diagnostics.

## Correction Audit

- 10 A broad bond cutoff is recorded in `broad_backbone_aggregated.csv` as `bond_cutoff_A=10`; sampled `ring_cutoff_A=8` and `mc_source_near_field_ratio=0.5`.
- `mc_source_is_self_or_bonded` is present in `broad_backbone_sources.csv`; valid-source aggregate columns `mc_lit_T2_local_valid_0..4` are present and the McConnell analyses were run with `--mc-source-mode valid`.
- Standalone `ring_current_sources.csv` contains `jb_T0`, `jb_T2_local_0..4`, and matching `jb_unit_T*` columns, so the literature-scaled Pople/JB correction is active in the ring check.
- Between-axis-led reads and thin/effective-N flags are present in the McConnell, ring, and variance analysis outputs.

## Per-Stratum Mechanism Table

All signals are correlate-not-match diagnostics. T2 rows show both between-atom and within-frame axes. Static scalar-field rows show between-led LOAO R2 and within-frame R2 because scalar sigma has no T2 component axis. Coefficient SEs are within-protein delete-atom jackknife SEs where emitted by the analysis script.

| stratum | mechanism | signal | equation / coefficient +/- SE | effective N | verdict bucket |
| --- | --- | --- | --- | --- | --- |
| N | ring current | T2 between c=0.205, absT2=0.221; within c=0.073, absT2=0.060 | T2 = gamma*K_ring; broad gamma=-391.704 +/- 110.695; Pople/JB pooled gamma=0.923 +/- 0.531 | atoms=54; AR1 Neff_med=2.220e+04 | recovered law |
| N | McConnell valid-source | T2 between c=0.958, absT2=0.186; within c=0.088, absT2=0.415 | T2 = q_CO*G_CO + q_CN*G_CN + q_SC*G_SC; q_CO=16.174 +/- 22.868, q_CN=-30.164 +/- 18.950, q_SC=-68.219 +/- 107.092 | atoms=54; atom_Neff=48.894; AR1=2.129e+04 | can't-make-it-work-for-now |
| N | charge q/r3 T2 | T2 between c=0.658, absT2=0.776; within c=0.123, absT2=0.229 | T2 = gamma*K_charge_q/r3; gamma=9.701 +/- 0.227 | atoms=54; AR1 Neff_med=2.220e+04 | form-recovered-scale-fitted |
| N | FF14SB static field | sigma between R2=-0.213; within R2=0.058 | sigma = A*E_z + B*E_mag^2; A=-332.708 +/- 349.654, B=201.964 +/- 295.398 | atoms=54; AR1 Neff_med=2.327e+04 | can't-make-it-work-for-now |
| N | APBS/Buckingham field | sigma between R2=-0.902; within R2=0.005 | sigma = A*E_proj + B*E_mag^2; A=-100.724 +/- 243.662, B=42.883 +/- 138.067 | atoms=54; AR1 Neff_med=2.327e+04 | can't-make-it-work-for-now |
| N | APBS EFG local-frame | T2 between c=-0.059, absT2=0.517; within c=0.168, absT2=0.082; local audit c=0.153, absT2=0.102 | T2 = gamma*K_APBS_EFG_local; gamma=-15.207 +/- 2.123; nonlinear local R2=0.023 | atoms=54; AR1 Neff_med=2.220e+04 | can't-make-it-work-for-now |
| N | AIMNet2 ceiling | T2 between R2=-1.399; within R2=0.587; lift vs physics=7.742/0.072 | ridge(physics + AIMNet2 charge + CRG + embedding PCs); alpha=10; ceiling diagnostic, not a physical coefficient | atoms=54; rows=27000; atom jackknife SE in CSV | form-recovered-scale-fitted |
| CA | ring current | T2 between c=-0.034, absT2=-0.132; within c=0.095, absT2=0.069 | T2 = gamma*K_ring; broad gamma=-6.554 +/- 19.440; Pople/JB pooled gamma=0.923 +/- 0.531 | atoms=54; AR1 Neff_med=2.080e+04 | recovered law |
| CA | McConnell valid-source | T2 between c=0.842, absT2=0.164; within c=0.083, absT2=0.393 | T2 = q_CO*G_CO + q_CN*G_CN + q_SC*G_SC; q_CO=-17.244 +/- 13.151, q_CN=-29.332 +/- 8.208, q_SC=15.515 +/- 59.442 | atoms=54; atom_Neff=32.337; AR1=2.018e+04 | can't-make-it-work-for-now |
| CA | charge q/r3 T2 | T2 between c=0.414, absT2=-0.036; within c=0.171, absT2=0.281 | T2 = gamma*K_charge_q/r3; gamma=5.253 +/- 0.569 | atoms=54; AR1 Neff_med=2.080e+04 | form-recovered-scale-fitted |
| CA | FF14SB static field | sigma between R2=-0.106; within R2=0.005 | sigma = A*E_z + B*E_mag^2; A=7.227 +/- 34.561, B=-55.988 +/- 231.068 | atoms=54; AR1 Neff_med=2.435e+04 | can't-make-it-work-for-now |
| CA | APBS/Buckingham field | sigma between R2=-0.037; within R2=0.035 | sigma = A*E_proj + B*E_mag^2; A=-7.581 +/- 5.936, B=8.795 +/- 3.817 | atoms=54; AR1 Neff_med=2.435e+04 | can't-make-it-work-for-now |
| CA | APBS EFG local-frame | T2 between c=-0.260, absT2=0.138; within c=0.078, absT2=0.033; local audit c=0.056, absT2=0.059 | T2 = gamma*K_APBS_EFG_local; gamma=-1.078 +/- 0.668; nonlinear local R2=0.005 | atoms=54; AR1 Neff_med=2.080e+04 | can't-make-it-work-for-now |
| CA | AIMNet2 ceiling | T2 between R2=0.082; within R2=0.374; lift vs physics=0.447/0.137 | ridge(physics + AIMNet2 charge + CRG + embedding PCs); alpha=10; ceiling diagnostic, not a physical coefficient | atoms=54; rows=27000; atom jackknife SE in CSV | form-recovered-scale-fitted |
| C | ring current | T2 between c=0.168, absT2=-0.078; within c=0.055, absT2=0.059 | T2 = gamma*K_ring; broad gamma=-263.802 +/- 98.951; Pople/JB pooled gamma=0.923 +/- 0.531 | atoms=54; AR1 Neff_med=2.361e+04 | recovered law |
| C | McConnell valid-source | T2 between c=0.992, absT2=0.156; within c=0.030, absT2=0.090 | T2 = q_CO*G_CO + q_CN*G_CN + q_SC*G_SC; q_CO=-26.012 +/- 14.402, q_CN=-2.515 +/- 17.270, q_SC=-29.588 +/- 72.738 | atoms=54; atom_Neff=43.149; AR1=1.956e+04 | can't-make-it-work-for-now |
| C | charge q/r3 T2 | T2 between c=0.659, absT2=0.150; within c=0.544, absT2=0.672 | T2 = gamma*K_charge_q/r3; gamma=10.217 +/- 0.292 | atoms=54; AR1 Neff_med=2.361e+04 | form-recovered-scale-fitted |
| C | FF14SB static field | sigma between R2=-0.098; within R2=0.052 | sigma = A*E_z + B*E_mag^2; A=2.709 +/- 14.351, B=23.633 +/- 28.596 | atoms=54; AR1 Neff_med=2.582e+04 | can't-make-it-work-for-now |
| C | APBS/Buckingham field | sigma between R2=-0.426; within R2=0.015 | sigma = A*E_proj + B*E_mag^2; A=-0.719 +/- 9.534, B=0.638 +/- 4.358 | atoms=54; AR1 Neff_med=2.582e+04 | can't-make-it-work-for-now |
| C | APBS EFG local-frame | T2 between c=0.524, absT2=-0.468; within c=0.122, absT2=0.042; local audit c=-0.121, absT2=0.046 | T2 = gamma*K_APBS_EFG_local; gamma=-22.572 +/- 1.274; nonlinear local R2=0.017 | atoms=54; AR1 Neff_med=2.361e+04 | can't-make-it-work-for-now |
| C | AIMNet2 ceiling | T2 between R2=-0.377; within R2=0.717; lift vs physics=3.483/0.067 | ridge(physics + AIMNet2 charge + CRG + embedding PCs); alpha=10; ceiling diagnostic, not a physical coefficient | atoms=54; rows=27000; atom jackknife SE in CSV | form-recovered-scale-fitted |
| O | ring current | T2 between c=0.212, absT2=-0.155; within c=0.031, absT2=0.037 | T2 = gamma*K_ring; broad gamma=-1.569e+03 +/- 333.419; Pople/JB pooled gamma=0.923 +/- 0.531 | atoms=54; AR1 Neff_med=2.248e+04 | recovered law |
| O | McConnell valid-source | T2 between c=0.993, absT2=0.015; within c=0.047, absT2=0.533 | T2 = q_CO*G_CO + q_CN*G_CN + q_SC*G_SC; q_CO=124.718 +/- 56.874, q_CN=-25.887 +/- 36.210, q_SC=14.114 +/- 179.204 | atoms=54; atom_Neff=50.168; AR1=2.134e+04 | can't-make-it-work-for-now |
| O | charge q/r3 T2 | T2 between c=0.358, absT2=0.244; within c=0.062, absT2=0.090 | T2 = gamma*K_charge_q/r3; gamma=157.880 +/- 13.409 | atoms=54; AR1 Neff_med=2.248e+04 | form-recovered-scale-fitted |
| O | FF14SB static field | sigma between R2=-0.082; within R2=0.026 | sigma = A*E_z + B*E_mag^2; A=-91.213 +/- 267.594, B=1.467e+03 +/- 1.495e+03 | atoms=54; AR1 Neff_med=2.385e+04 | can't-make-it-work-for-now |
| O | APBS/Buckingham field | sigma between R2=0.091; within R2=0.003 | sigma = A*E_proj + B*E_mag^2; A=263.432 +/- 86.688, B=116.167 +/- 41.425 | atoms=54; AR1 Neff_med=2.385e+04 | can't-make-it-work-for-now |
| O | APBS EFG local-frame | T2 between c=0.769, absT2=-0.045; within c=0.024, absT2=0.008; local audit c=-0.019, absT2=0.023 | T2 = gamma*K_APBS_EFG_local; gamma=174.210 +/- 4.914; nonlinear local R2=-5.040e-04 | atoms=54; AR1 Neff_med=2.248e+04 | can't-make-it-work-for-now |
| O | AIMNet2 ceiling | T2 between R2=-0.063; within R2=0.419; lift vs physics=0.613/0.089 | ridge(physics + AIMNet2 charge + CRG + embedding PCs); alpha=10; ceiling diagnostic, not a physical coefficient | atoms=54; rows=27000; atom jackknife SE in CSV | form-recovered-scale-fitted |
| HN | ring current | T2 between c=0.079, absT2=-0.154; within c=0.205, absT2=0.241 | T2 = gamma*K_ring; broad gamma=18.388 +/- 10.130; Pople/JB pooled gamma=0.923 +/- 0.531 | atoms=52; AR1 Neff_med=1.788e+04 | recovered law |
| HN | McConnell valid-source | T2 between c=0.812, absT2=0.482; within c=0.178, absT2=0.397 | T2 = q_CO*G_CO + q_CN*G_CN + q_SC*G_SC; q_CO=1.148 +/- 6.154, q_CN=8.617 +/- 5.913, q_SC=12.243 +/- 19.621 | atoms=52; atom_Neff=44.973; AR1=2.175e+04 | can't-make-it-work-for-now |
| HN | charge q/r3 T2 | T2 between c=-0.002, absT2=-0.355; within c=0.278, absT2=0.364 | T2 = gamma*K_charge_q/r3; gamma=0.803 +/- 0.189 | atoms=52; AR1 Neff_med=1.788e+04 | form-recovered-scale-fitted |
| HN | FF14SB static field | sigma between R2=0.483; within R2=0.207 | sigma = A*E_z + B*E_mag^2; A=-22.431 +/- 6.419, B=38.804 +/- 27.371 | atoms=52; AR1 Neff_med=2.079e+04 | form-recovered-scale-fitted |
| HN | APBS/Buckingham field | sigma between R2=0.397; within R2=0.104 | sigma = A*E_proj + B*E_mag^2; A=-15.625 +/- 2.805, B=-9.338 +/- 2.047 | atoms=52; AR1 Neff_med=2.079e+04 | form-recovered-scale-fitted |
| HN | APBS EFG local-frame | T2 between c=0.469, absT2=0.086; within c=0.116, absT2=0.202; local audit c=0.120, absT2=0.229 | T2 = gamma*K_APBS_EFG_local; gamma=2.006 +/- 0.314; nonlinear local R2=0.024 | atoms=52; AR1 Neff_med=1.788e+04 | can't-make-it-work-for-now |
| HN | AIMNet2 ceiling | T2 between R2=0.198; within R2=0.460; lift vs physics=0.695/0.100 | ridge(physics + AIMNet2 charge + CRG + embedding PCs); alpha=10; ceiling diagnostic, not a physical coefficient | atoms=52; rows=26000; atom jackknife SE in CSV | form-recovered-scale-fitted |
| HA | ring current | T2 between c=0.303, absT2=0.228; within c=0.411, absT2=0.469 | T2 = gamma*K_ring; broad gamma=-19.535 +/- 4.482; Pople/JB pooled gamma=0.923 +/- 0.531 | atoms=58; AR1 Neff_med=1.698e+04 | recovered law |
| HA | McConnell valid-source | T2 between c=0.625, absT2=0.414; within c=0.161, absT2=0.258 | T2 = q_CO*G_CO + q_CN*G_CN + q_SC*G_SC; q_CO=5.151 +/- 3.528, q_CN=-9.796 +/- 2.832, q_SC=0.397 +/- 13.035 | atoms=58; atom_Neff=42.120; AR1=2.354e+04 | can't-make-it-work-for-now |
| HA | charge q/r3 T2 | T2 between c=0.199, absT2=-0.338; within c=0.092, absT2=0.272 | T2 = gamma*K_charge_q/r3; gamma=1.819 +/- 0.249 | atoms=58; AR1 Neff_med=1.698e+04 | form-recovered-scale-fitted |
| HA | FF14SB static field | sigma between R2=-0.018; within R2=0.032 | sigma = A*E_z + B*E_mag^2; A=-6.312 +/- 3.922, B=-27.879 +/- 50.018 | atoms=58; AR1 Neff_med=2.098e+04 | can't-make-it-work-for-now |
| HA | APBS/Buckingham field | sigma between R2=-0.321; within R2=0.021 | sigma = A*E_proj + B*E_mag^2; A=0.017 +/- 0.844, B=-9.292e-04 +/- 1.954 | atoms=58; AR1 Neff_med=2.098e+04 | can't-make-it-work-for-now |
| HA | APBS EFG local-frame | T2 between c=0.001, absT2=-0.104; within c=0.013, absT2=0.144; local audit c=0.010, absT2=0.151 | T2 = gamma*K_APBS_EFG_local; gamma=1.971 +/- 0.470; nonlinear local R2=-0.005 | atoms=58; AR1 Neff_med=1.698e+04 | can't-make-it-work-for-now |
| HA | AIMNet2 ceiling | T2 between R2=0.016; within R2=0.401; lift vs physics=0.645/0.063 | ridge(physics + AIMNet2 charge + CRG + embedding PCs); alpha=10; ceiling diagnostic, not a physical coefficient | atoms=58; rows=25000; atom jackknife SE in CSV | form-recovered-scale-fitted |

## Equation Fits

Backbone distillation evidence was regenerated from the corrected broad substrate with six strata, default 4000-epoch fits, and atom-split validation. The symbolic PySR pass was run from `src/rediscover/analysis/venv` with per-category bond fits enabled.

| mechanism | radial form read | learned radial fit | fixed comparator |
| --- | --- | --- | --- |
| ring | learned log-power range -2.67..-1.78 | learned radial R2 range 0.065..0.212 | fixed r^-3 comparator R2 range 0.001..0.245 |
| bond | learned log-power range -2.86..-2.08 | learned radial R2 range 0.515..0.612 | fixed r^-3 comparator R2 range 0.172..0.979 |
| charge | learned log-power range -4.11..-2.12 | learned radial R2 range 0.100..0.443 | fixed r^-3 comparator R2 range 0.029..0.981 |

| PySR label | variables | best equation | R2 | sample rows |
| --- | --- | --- | --- | --- |
| ring | r,cos_theta,ring_intensity | square(cos_theta) / (r / (ring_intensity + 14.574982)) | 0.0646 | 6000 |
| bond_pooled | r,cos_theta,bond_category | (-0.098248325 / (square((r + -2.186969) - cos_theta) + 0.0073036416)) / 2.700291 | 0.2618 | 6000 |
| bond_cat0 | r,cos_theta | (square(cos_theta) + -0.39135703) / (r + -1.946662) | 0.0458 | 6000 |
| bond_cat1 | r,cos_theta | square(square(cos_theta / r)) * -83.444695 | 0.0260 | 6000 |
| bond_cat3 | r,cos_theta | square(square(square((r / (square(r - (cos_theta * 3.5472)) + 2.454844)) / 2.0084963) * 2.194088)) | 0.1252 | 6000 |
| bond_cat4 | r,cos_theta | square(1.8261545 / r) | 0.0633 | 6000 |
| charge | r,source_q_e | square(square(((3.15573 - r) / (source_q_e + 0.46519893)) / square(r))) | 0.7749 | 6000 |

Layer-2 calibrated coefficients are the per-stratum coefficients in the main table: Pople/JB ring gamma from `jb_T*`, McConnell valid-source category q values, charge q/r3 gamma, APBS local-frame EFG gamma, FF14SB field A/B, and APBS/Buckingham field A/B. The PySR equations are learned-gate diagnostics, not replacements for those calibrated physical readouts.

## Artifact Paths

- broad substrate: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/broad_backbone`
- standalone ring substrate: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/ring_current`
- EFG substrate: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/efg`
- Buckingham substrate: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/buckingham_efield`
- AIMNet2 substrate: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/aimnet2_features`
- ring de-circ CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/ring_literature_decirc/ring_literature_decirc.csv`
- McConnell valid de-circ CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/mcconnell_literature_decirc/mcconnell_literature_decirc.csv`
- McConnell dchi CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/mcconnell_dchi_calibration/mcconnell_dchi_calibration.csv`
- static calibration CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/static_environment_calibration/static_environment_calibration.csv`
- variance decomposition CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/variance_decomposition/variance_decomposition.csv`
- EFG local audit CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/efg_localframe_audit/efg_localframe_audit.csv`
- EFG distill CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/efg_distill/efg_distill_summary.csv`
- AIMNet2 ensemble CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/aimnet2_ensemble/aimnet2_ceiling_ensemble.csv`
- backbone radial fits CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/backbone_law_evidence/radial_fit_summary.csv`
- PySR summary CSV: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/backbone_law_evidence/pysr_summary_gate_chosen.csv`

## Chart Paths

- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/variance_shares_sigma_iso.png`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/variance_shares_sigma_iso.svg`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/variance_shares_T2.png`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/variance_shares_T2.svg`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/decircularising_correlations.png`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/decircularising_correlations.svg`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/t2_fit_metrics.png`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/t2_fit_metrics.svg`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/efg_localframe_correction.png`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/efg_localframe_correction.svg`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/efg_localframe_fit.png`
- `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/capstone_charts/efg_localframe_fit.svg`

## Discipline Notes

- Per-stratum results are shown for N, CA, C, O, HN, and HA throughout; HA2/HA3 are not promoted into the capstone table.
- Confidence is within this single protein: delete-atom jackknife for fitted coefficients and atom-level metrics, plus AR(1) effective-frame diagnostics where the source scripts emit them.
- Static mechanisms are between-axis-led in the coefficient readout. Within-frame values are diagnostics, not population confidence.
- This report uses the allowed verdict buckets only: `recovered law`, `form-recovered-scale-fitted`, and `can't-make-it-work-for-now`.
- Future 750-DFT reruns are parameterized by passing the new emitted out-dir paths to the same analysis scripts; no hard-coded row count is required in the analysis consumers.
