# TRUE LOAO Postmortem - 2026-06-04
Run dir: `/tmp/rediscover-runs/2026-06-04-true-loao`
Substrate: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`; Build4 CSV/NPY only; no extraction/DFT; frozen `get_C`, `|C.T C-I|max=1.11e-16`.
TRUE LOAO: train on training-atom means only; transform held-out atom means with training-atom feature/target stats; score physical held-out atom means.
Null: 1000 deranged shuffles of target atom-means across atoms; structural features and LOAO folds unchanged.
Support note: charge N is the actual Stage2 charge-wide subset tied to old 0.38; the older 54-atom backbone stratum is a different comparison.
| kernel | N atoms | old mislabeled within-LOAO R2 | TRUE-LOAO R2 | delta | null p/z/pct | determinability | lead-scale | support |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- | --- |
| charge | 13 | 0.377 | 0.036 | -0.341 | p 0.001; z 5.700; pct 100.000 | single q/r3; Stage2 charge-wide subset, N-limited case-study | ~null / 0.03-class (case-study N) | thin_atoms_case_study |
| ring | 5 | 0.165 | -1.041 | -1.206 | p 0.940; z -1.180; pct 6.000 | single JB ring shadow; N=5 case-study, probability discrete | negative/null; no atom-mean recovery (case-study N) | thin_atoms_case_study |
| unified_Dab_sum | 25 | 0.258 | -104.579 | -104.837 | p 0.698; z 0.060; pct 30.200 | unified drop-one: mopac_field_backbone +333.631; mc_category_4 +153.668; mc_category_0 +100.840; charge_total +52.528 | negative/null; no atom-mean recovery (case-study N) | thin_TRUE_LOAO_p_ge_atoms;thin_atoms_case_study |
Unified drop-one determinability: mopac_field_backbone +333.631; mc_category_4 +153.668; mc_category_0 +100.840; charge_total +52.528.
One-line read: TRUE between-atom recovery on 1P9J is ~null; between-axis probability should move to the 720-WT.
Difference is explicit: the old number is within-atom modulation because the held-out atom was centered by its own mean; TRUE LOAO is between-atom atom-mean recovery.
Artifacts: `true_loao_recovery.csv`, `true_loao_unified_drop_one.csv`, `true_loao_atom_predictions.csv`, `true_loao_null_scores.csv`, `true_loao_run_audit.json`.
Disk: `/tmp/rediscover-runs` 6.7G before write (<15G); output drop-old=True; build4+build1 kept.
