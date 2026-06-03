# Build 3 postmortem - fit architecture loop
Run dir: `/tmp/rediscover-runs/2026-06-03-build3-fit-arch/`
Script commit: `d35d7ec46e2c6b2064ab90228b6a997e4714e11a` on `h5-reader-pysr-spike`; script only committed, heavy results not committed.

## Total-T2, tier=all, per-type vs global-sliced
| stratum | n_type | global_b | per_b | d_b | global_w | per_w | d_w | support |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| HN | 52 | -54.380 | -0.051 | +54.329 | -76.872 | 0.578 | +77.450 | thin_between<100,p>=atoms |
| HA | 58 | -16.406 | 0.047 | +16.453 | -24.738 | 0.307 | +25.045 | thin_between<100,p>=atoms |
| aromatic_H | 41 | -23.711 | -0.049 | +23.663 | -29.372 | 0.499 | +29.872 | thin_between<100,p>=atoms |
| polar_H | 40 | -12.398 | -0.048 | +12.350 | -10.839 | 0.717 | +11.556 | thin_between<100,p>=atoms |
| N | 54 | 0.628 | 0.519 | -0.108 | 0.809 | 0.911 | +0.102 | thin_between<100,p>=atoms |
| aromatic_heavy | 69 | 0.443 | 0.639 | +0.196 | 0.491 | 0.670 | +0.179 | thin_between<100,p>=atoms |
| O | 54 | 0.081 | 0.198 | +0.117 | 0.182 | 0.718 | +0.536 | thin_between<100,p>=atoms |

H strata numeric read: global anti-prediction collapses to near-zero between and positive within under per-type.
Strong-strata read: N/aromatic_heavy/O retain positive per-type within; N pays between cost vs global, aromatic_heavy/O do not.
Determinability cost: per-type between uses 7-213 atoms by stratum, 40-69 for the headline rows, vs 846 global atoms; amide_sidechain=9 and sulfur=7 are also `thin`.

## Channel read, tier=all, within R2
| stratum | scope | total | dia | para |
| --- | --- | ---: | ---: | ---: |
| HN | per_type | 0.578 | 0.246 | 0.245 |
| HA | per_type | 0.307 | 0.195 | 0.195 |
| aromatic_H | per_type | 0.499 | 0.199 | 0.197 |
| polar_H | per_type | 0.717 | 0.325 | 0.323 |
| N | per_type | 0.911 | 0.204 | 0.232 |
| aromatic_heavy | per_type | 0.670 | 0.135 | 0.184 |
| O | per_type | 0.718 | 0.350 | 0.430 |
Gauge-dependent split: report only; total/dia/para R2 are not additive mechanism attribution.

## Dominance response curves
Dominance rows: 14,080; targets=total/dia/para/T1, scopes=global_sliced/per_type, mechanisms=ring/charge/mc.
total-T2 all/within global shapes: charge U18 flat5 fall9 rise6 threshold17 thin26; mc flat5 fall10 rise4 threshold25 thin26; ring U19 fall5 rise5 threshold30 thin20.
total-T2 all/within per-type shapes: charge U11 fall14 rise6 threshold24 thin26; mc U10 flat5 fall5 rise5 threshold19 thin26; ring U9 flat5 rise5 threshold40 thin20.
C++ next emit flag: `dominant_fraction_{ring,charge,mc}` bin-id should ride the substrate; Build 3 used Python quantiles on C++ scalars only to avoid re-emitting.

## Integrity
Score rows: 496; partition rows: 185,056; dominance bins are input-side only; no target/residual/coef conditioners.
Fit checks passed: input_acceptance, no_external_merge, basis_and_targets, cv_integrity, partition_integrity.
Held-out everywhere; same 558,360 charge-complete rows; train-only PCA; frozen `change_of_basis.get_C()`; no H5/per-source/older-dir reads.
Disk guard: df before write passed with 163.6 GiB free, `/tmp/rediscover-runs` 3.21 GiB before write; result dir 165M, total rediscover 3.4G; never `/shared`.
SETI: numbers and shapes only, no mechanism verdict.
