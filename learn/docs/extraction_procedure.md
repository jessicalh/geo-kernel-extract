# Running a Feature Extraction

## What a "run" is

A run is a named extraction: one version of the extractor binary,
one set of calculator parameters, applied to the 723 WT-ALA mutant
pairs.  Each run gets its own directory under
`calibration/features/{run_name}/`.  Never overwrite a run.  If the
extractor changes, make a new run with a new name.

Previous runs:

| Run | When | Notes |
|-----|------|-------|
| FirstExtraction | early | Initial test, old extractor |
| CalibrationExtractionTest | early | 66 proteins, pre-SDK |
| GatedCalibration | 2026-04 | Full 723, pre-azimuthal |
| AzimuthalExtraction | 2026-04-09 | Current.  112/723 done.  Added cos/sin phi, per-ring azimuthal. |

## Pre-flight

1. **Build the extractor.**  From repo root:
   ```bash
   cd /shared/2026Thesis/nmr-shielding
   cmake --build build --target nmr_extract -j$(nproc)
   ```
   The binary is `build/nmr_extract`.  extract.py checks it exists.

2. **Verify calibration data.**  723 protein directories should exist:
   ```bash
   ls calibration/ | wc -l   # expect ~725 (723 + features/ + maybe a few meta)
   ```
   Each has `{ID}_WT.xyz`, `{ID}_WT.prmtop`, `{ID}_ALA.xyz`,
   `{ID}_ALA.prmtop`, and optionally `{ID}_WT_nmr.out`.

3. **Choose a run name.**  Add it to `KNOWN_RUNS` in `extract.py`
   to suppress the warning, or just accept the warning.

## Running

From the repo root (not learn/):

```bash
# Dry run — shows what would execute, extracts nothing
python learn/extract.py --run NewRunName --dry-run

# Full run — sequential, one protein at a time
# MOPAC is ~10 min/protein.  723 proteins ≈ 120 hours.
python learn/extract.py --run NewRunName

# Resume after interruption — skips proteins with pos.npy already present
python learn/extract.py --run NewRunName --resume

# Single protein (for testing)
python learn/extract.py --run NewRunName --protein P84477

# With calculator parameter overrides
python learn/extract.py --run NewRunName --config data/calculator_params.toml
```

For long runs, use tmux or nohup:

```bash
nohup python -u learn/extract.py --run NewRunName --resume \
  > calibration/features/NewRunName/extract_stdout.log 2>&1 &
```

The `-u` flag is important — unbuffered output so the log is readable
while the extraction runs.

## What it produces

Per protein: a directory `calibration/features/{run}/{ID}/` containing
NPY arrays (currently ~60 files per protein).  The exact set depends
on the extractor version — new calculators add new arrays.

Per run: `extract_log.jsonl` in the run directory.  One JSON line per
protein with `protein_id`, `ok`, `error`, `elapsed`, `n_arrays`.

## Monitoring

```bash
# How many done?
ls calibration/features/NewRunName/ | wc -l

# Any failures?
grep '"ok": false' calibration/features/NewRunName/extract_log.jsonl

# Tail the log
tail -f calibration/features/NewRunName/extract_stdout.log
```

## Pointing the analysis at the new run

After extraction completes, update `learn/src/calibration.toml`:

```toml
[paths]
features = "/shared/2026Thesis/nmr-shielding/calibration/features/NewRunName"
```

Then re-run all analysis scripts per `src/actual_physics/OVERVIEW.md`.
The scripts read from `cfg.paths.features` — nothing else needs to
change.

## Cost

MOPAC dominates extraction time.  Approximate per-protein costs:

| Calculator | Time |
|-----------|------|
| Geometry + SpatialIndex | < 1s |
| Coulomb + McConnell + HBond | < 1s |
| BiotSavart + HaighMallion | < 1s |
| PiQuadrupole + Dispersion + RingSusc | < 1s |
| DSSP + SASA | < 1s |
| MOPAC (PM7/MOZYME) | ~8-12 min |
| AIMNet2 | ~0.2s (GPU) |
| APBS | ~4s |

Total per protein: ~10 min (dominated by MOPAC).
Total for 723 proteins: ~120 hours sequential.

If MOPAC is not needed (geometry-only analysis), pass a config that
disables it.  The geometry-only basis analysis in
`src/actual_physics/geometry_only_basis.py` shows that 44/55 core
kernels survive without MOPAC, with the full→geo-only R² gap being
only 0.007 for H.

## Never overwrite

Extraction cost is real.  A completed run is an asset.  If you need
to change something:
- Changed extractor binary → new run name
- Changed calculator params → new run name
- Need to add proteins → `--resume` on same run
- Something went wrong → fix and `--resume`, don't delete
