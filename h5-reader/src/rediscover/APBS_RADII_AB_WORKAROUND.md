# APBS Static PB-Radii A/B and Workaround Status

Date: 2026-06-02

Branch: `h5-reader-pysr-spike`

Scope: ORCA-free APBS-only test. This does not rewire or redo trajectory
producer extraction. It tests a workaround-stage substitution of static real PB
radii into APBS while keeping the same 1P9J trajectory positions and TPR stored
charges.

## Context

The 1P9J trajectory TPR path installs the compatibility PB-radius placeholder
`1.5 A` for every protein atom. APBS consumes `partial_charge` and `pb_radius`
directly, so the placeholder changes the APBS dielectric boundary, even though
it does not change positions or charges.

The static real-radius source used here is:

- `../data/ff14sb_params.dat`
- header provenance: `AmberTools tleap`, `leaprc.protein.ff14SB`
- PB-radius model: `mbondi2`
- table column: `PB_RADIUS(Angstrom)`

No local 1P9J-specific prmtop was found. The local prmtops found by search are
unrelated calibration/test systems, so the cleanest available static source for
1P9J is the flat ff14SB parameter table. Mapping coverage was complete:
`internal_matches=808 terminal_matches=38 atoms=846`.

Charges were not replaced by flat-param charges. Each A/B branch uses the same
TPR stored charges. Only the PB radii differ:

- placeholder branch: TPR charges + uniform `1.5 A` PB radii
- real branch: TPR charges + ff14SB flat `pb_radius` values

Radius distribution for the real branch:

| source | n | min | max | mean | median abs delta vs 1.5 A | counts |
|---|---:|---:|---:|---:|---:|---|
| ff14SB flat `pb_radius` | 846 | 1.2 | 1.8 | 1.44001 | 0.2 | `1.2:319;1.3:85;1.5:81;1.55:73;1.7:281;1.8:7` |

## Harness

Added and built:

- `../tools/apbs_radii_ab_workaround.cpp`
- CMake target: `apbs_radii_ab_workaround`

Build command:

```bash
cmake --build ../build --target apbs_radii_ab_workaround -j2
```

A/B command:

```bash
../build/apbs_radii_ab_workaround \
  --trajectory-dir ../tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z \
  --h5 ../molprobity_runs/1P9J_5801/extract/trajectory.h5 \
  --params ../data/ff14sb_params.dat \
  --out /tmp/apbs-radii-ab-workaround \
  --frames 0,250,501
```

The harness constructs APBS inputs at analysis/workaround time, installs the
intended charge/radius table on each seeded `TrajectoryProtein`, and solves only
`ChargeAssignmentResult` + `ApbsFieldResult`. It does not use ORCA and does not
rewrite producer trajectory extraction.

The APBS logs confirm the intended A/B:

- placeholder solves warn: 846 atoms on non-authoritative placeholder PB radii
- real solves use: `ff14SB:tpr_charges+flat_pb_radius:../data/ff14sb_params.dat`

Sample output files:

- `/tmp/apbs-radii-ab-workaround/ab_metrics.csv`
- `/tmp/apbs-radii-ab-workaround/radii_summary.csv`
- `/tmp/apbs-radii-ab-workaround/placeholder_1p5A/frame_*/apbs_E.npy`
- `/tmp/apbs-radii-ab-workaround/placeholder_1p5A/frame_*/apbs_efg.npy`
- `/tmp/apbs-radii-ab-workaround/real_static_radii/frame_*/apbs_E.npy`
- `/tmp/apbs-radii-ab-workaround/real_static_radii/frame_*/apbs_efg.npy`

## A/B Frames

| h5_row | original_index | time_ps |
|---:|---:|---:|
| 0 | 0 | 0 |
| 250 | 500 | 5000 |
| 501 | 1002 | 10020 |

Surface/buried stratification uses per-frame atom SASA from `SasaResult` on the
same positions. `buried_sasa_q1` is the lowest-SASA quartile, and
`surface_sasa_q4` is the highest-SASA quartile.

## A/B Metrics

All atoms:

| metric | average across 3 frames | range |
|---|---:|---:|
| E-vector cosine median | 0.999727 | 0.999700 to 0.999770 |
| E-vector cosine p05 | 0.983208 | 0.982315 to 0.983816 |
| E magnitude ratio median, real/placeholder | 0.988844 | 0.987758 to 0.990006 |
| E magnitude ratio p05 | 0.884617 | 0.879477 to 0.893990 |
| E magnitude ratio p95 | 1.037590 | 1.034790 to 1.041560 |
| EFG T2 component correlation | 0.999369 | 0.999336 to 0.999389 |
| EFG `|T2|` correlation | 0.999532 | 0.999468 to 0.999614 |
| EFG magnitude ratio median, real/placeholder | 1.000680 | 1.000580 to 1.000760 |
| EFG magnitude ratio p05 | 0.947750 | 0.945806 to 0.951023 |
| EFG magnitude ratio p95 | 1.050430 | 1.047120 to 1.054650 |

Buried SASA Q1:

| metric | average across 3 frames | range |
|---|---:|---:|
| E-vector cosine median | 0.999868 | 0.999835 to 0.999904 |
| E-vector cosine p05 | 0.993634 | 0.992150 to 0.995236 |
| E magnitude ratio median, real/placeholder | 0.995088 | 0.993824 to 0.996034 |
| E magnitude ratio p05 | 0.937100 | 0.934649 to 0.941262 |
| E magnitude ratio p95 | 1.030500 | 1.027660 to 1.032010 |
| EFG T2 component correlation | 0.999847 | 0.999810 to 0.999880 |
| EFG `|T2|` correlation | 0.999881 | 0.999862 to 0.999905 |
| EFG magnitude ratio median, real/placeholder | 0.999582 | 0.999195 to 1.000060 |
| EFG magnitude ratio p05 | 0.975621 | 0.968490 to 0.981408 |
| EFG magnitude ratio p95 | 1.022340 | 1.018360 to 1.026800 |

Surface SASA Q4:

| metric | average across 3 frames | range |
|---|---:|---:|
| E-vector cosine median | 0.999341 | 0.999200 to 0.999430 |
| E-vector cosine p05 | 0.965367 | 0.952093 to 0.972401 |
| E magnitude ratio median, real/placeholder | 0.978762 | 0.977137 to 0.980224 |
| E magnitude ratio p05 | 0.843819 | 0.829833 to 0.856395 |
| E magnitude ratio p95 | 1.064510 | 1.037210 to 1.080700 |
| EFG T2 component correlation | 0.998646 | 0.998576 to 0.998714 |
| EFG `|T2|` correlation | 0.999025 | 0.998901 to 0.999254 |
| EFG magnitude ratio median, real/placeholder | 1.005960 | 1.001510 to 1.012070 |
| EFG magnitude ratio p05 | 0.911631 | 0.900371 to 0.933511 |
| EFG magnitude ratio p95 | 1.074860 | 1.061010 to 1.085920 |

## Verdict

Engineering verdict from this A/B: the placeholder PB radii are a modest APBS
perturbation, not a material field/EFG change for these 1P9J frames.

The field directions are nearly unchanged, all-atom median field magnitudes
shift by about `1-2%`, and EFG T2 correlations remain about `0.9993-0.9999`.
Surface atoms show wider tails, especially field magnitude p05/p95, but the
central APBS field/EFG behavior stays very close to the placeholder run.

Scientific interpretation is reserved for the lead. The narrow conclusion here
is that the APBS "can't-make-work" issue is not explained by the uniform
`1.5 A` PB-radius placeholder alone.

## Workaround-Fixed Field/EFG

Not applied.

The requested workaround-fixed field/EFG sidecar is conditional on the A/B being
material. Because this run classifies the placeholder effect as modest, I did
not re-emit the full rediscover-frame APBS field/EFG set.

The harness already has a conditional sidecar path for future use:

```bash
../build/apbs_radii_ab_workaround \
  --trajectory-dir ../tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z \
  --h5 ../molprobity_runs/1P9J_5801/extract/trajectory.h5 \
  --params ../data/ff14sb_params.dat \
  --out /tmp/apbs-radii-ab-workaround \
  --frames 0,250,501 \
  --emit-rows-csv /tmp/rdc-buckingham-capstone/buckingham_efield_aggregated.csv
```

If run, that emits `/tmp/apbs-radii-ab-workaround/real_static_radii_apbs.h5`
with `apbs_E` and `apbs_efg` for the unique `h5_row` values from the rediscover
CSV, using the documented workaround provenance.
