# Stage 2.2 Postmortem - 2026-06-04

> **Historical run record — not current truth (trued 2026-06-04).** Preserve as
> provenance. Positive 1P9J LOAO/between claims here are superseded by
> `POSTMORTEM_TRUE_LOAO_2026-06-04.md`; check `NOW.md`/`STATE.md` before
> quoting.

Run dir: `/tmp/rediscover-runs/2026-06-04-stage2_2-unified-vet`
Substrate: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`
Scope: Build4 CSV/NPY sidecars only; no extraction/DFT/per-source; frozen `get_C`, `|C.T C-I|max=1.11e-16`; five-component total-T2.
Selection: D_ab CaseHunter atoms + ring dominance<0.7 + max D_ab dominance>=0.5; input-side dominance only; 25 atoms / 3903 rows / 26 unified terms.
Overfit vet: OLS 0.432/0.258 -> ridge ALL within 0.435 (JK 0.104,0.766; BB 0.397,0.464); LOAO 0.101 (JK 0.021,0.180; BB 0.094,0.104).
Shrinkage: within alpha 100, effDOF 19.3/26 (+intercept 20.3); within N_eff 515.8.
LOAO shrinkage: alpha mode/min/max 100000/100000/100000; effDOF median 3.1/26; support `thin_LOAO_p_ge_atoms;thin_atoms`.
Bucket: overfit=within-survives; LOAO-shrinks-to-thin-positive; old LOAO 0.26 does not survive as-is and should not outrank within.
Charge table: charge-alone within 0.007, LOAO 0.037; minus-charge within 0.408, LOAO 0.081; all within 0.435, LOAO 0.101.
Deltas ALL-charge: within +0.428, LOAO +0.064; ALL-minus-charge: within +0.027, LOAO +0.020; bucket=REAL combine on within; LOAO support thin.
Drop-one losses (within/LOAO, positive means term carries recovery): mopac_field_backbone +0.194/+0.041; mc_category_4 +0.171/+0.029; mc_category_0 +0.056/+0.012; charge_total +0.027/+0.020.
Drop-one removals that improve within (negative loss = noisy/overfit term): mc_category_1 -0.018/+0.003; mc_category_3 -0.007/+0.001; mc_category_2 -0.006/-0.000.
Verdict: real combine, not charge-in-a-coat; within survives shrinkage, LOAO is overfit-sensitive/thin-positive.
Artifacts: `stage2_2_charge_ablation_recovery.csv`, `stage2_2_drop_one_ablation.csv`, `stage2_2_unified_all_ridge_intensities.csv`, `stage2_2_run_audit.json`.
Disk: `/tmp/rediscover-runs` 6.7G before write (<15G); output drop-old=False; build4+build1 kept.
