# Build 2 Postmortem

Run dir: `/tmp/rediscover-runs/2026-06-03-build2-partition`.
Script commit: `a2d69cbb829dcc4901486c1162567fc99ba38b95`.
Substrate: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1`.
Read scope: manifest/specs/rows, classical/hbond/method-path reduced sidecars, partition bin IDs, embedding, total-T2/para-T2/T1 targets, CaseHunter manifests only.
Disk guard: PASS; free before write 128.3 GiB, `/tmp/rediscover-runs` 3.16 GiB, result dir 53 MiB.
Checks: input, basis/targets, CV, no-external-merge, partition integrity all PASS.

All-tier held-out R2 B/W core strata, total-T2: ALL .096/.144; N .628/.809; CA -1.188/-2.599; C -.164/-.336; O .081/.182; HN -54.380/-76.872; HA -16.406/-24.738.
All-tier held-out R2 B/W core strata, para-T2: ALL .001/.075; N -.010/.034; CA -.003/.008; C -.012/.036; O -.003/.084; HN -.013/.002; HA -.006/.046.
All-tier held-out R2 B/W core strata, T1 diagnostic: ALL -.047/.006; N -.017/.029; CA -.196/-.045; C -.045/-.131; O -.090/-.002; HN -7.129/-1.946; HA -1.511/-.769.
Fine sidechain support: amide_sidechain n=9 and sulfur n=7 are thin; retained with flags, not promoted.

+newmech lift vs classical on HN/O, total-T2: HN -14.575/-28.352; O +.011/+.030.
+newmech lift vs classical on HN/O, para-T2: HN +.003/+.086; O -.002/+.030.
+newmech lift vs classical on HN/O, T1 diagnostic: HN -4.055/-.456; O -.013/+.001.
Headline read: new mechanisms lift O within-frame for total/para and para-HN within-frame; they do not rescue total-T2 HN. T1 is diagnostic only.

Partition curves: 64,116 rows; every bin carries rows, atom count, N_eff, support_flag, and response_shape.
Supported shape read, total-T2: charge/field/mc/ring driver modulation is mostly threshold or U; MC magnitude on N shows monotone rise; charge magnitude on N often monotone fall; ring/charge geometry has threshold/rise/fall but many thin outer bins.
Supported shape read, para-T2: charge/field/mc driver magnitude is mostly fall or threshold; ring geometry is mainly threshold; HN/O gains are within-frame, not between-atom.
Supported shape read, T1 diagnostic: field/charge/mc driver modulation is threshold or U; MC magnitude includes monotone rise; support and charge-agreement fine categories are mostly unstable-thin.
Support policy: thin_atoms/thin_rows/thin_N_eff bins are flagged and excluded from the favorable shortlist.

Charge partition fixed: `bin_nearest_charge_r` bins [0,1,2,3], `bin_gap_to_2nd_charge_r` bins [0,1,2,3], `bin_abs_charge_T2` and `bin_sd_charge_T2` bins [0,1,2,3,4].
Note: this substrate has no separate emitted `dom_charge` bin column; dominant-charge scalar exists in Build 1 conditioning but Build 2 did not re-bin it in Python.

Favorable partitions: top supported habitats are total-T2 N strata, especially MC magnitude Q4/Q5 monotone-rise, charge modulation Q2 threshold, field modulation Q4 threshold, aimnet2 modulation Q1/Q2 threshold, and ring modulation Q1/Q4 threshold.
CaseHunter intersection: 60 navigable rows in `partition_casehunter_shortlist.csv`.
Top intersected habitats: charge/N atom 287 windows 877-884, 157-164, 0-4 via `sd_charge_T2_by_atom` Q2; charge/N atom 330 window 37-44 via `sd_mopac_coulomb_T2_by_atom` Q2 and `abs_charge_T2` Q2.
Additional intersected charge/N atoms include 503, 414, 489, 592, and 513 across the emitted windows; ranking is by partition held-out R2 then CaseHunter score.

Artifacts: `allatom_fit_score_table.csv`, `partition_response_curves.csv`, `partition_favorable_partitions.csv`, `partition_casehunter_shortlist.csv`, plots, and `run_audit.json`.
Deferred: ring/field cross-method validations and field divergence network remain step 2.
SETI: numbers and curve shapes only; no verdicts.
