# Rediscover Run Framework

This is a script-level convention for rediscover test runs. It does not change
C++ code, does not rebuild anything, and does not run ORCA.

## Root And Layout

Use one managed root for rediscover runs:

```bash
export REDISCOVER_RUN_ROOT=/tmp/rediscover-runs
```

Default root: `/tmp/rediscover-runs`.

Each run lives at:

```text
<root>/<run-name>/
```

`<run-name>` must be a single path component using letters, digits, `.`, `_`,
and `-`. Prefer names that carry the case and date, for example
`2026-06-02-broad-backbone-750`.

The C++ extractor writes emitted substrate files at the top of the run
directory. The current extractor CLI uses `--out`, so the wrapper invokes it as:

```bash
tools/rediscover_run.py emit 2026-06-02-broad-backbone-750 \
  --run /path/to/calcset-or-run.LGS \
  --bin build/linux-gcc/h5reader_extract \
  --case broad_backbone
```

If a future extractor binary exposes `--out-dir`, pass
`--extract-out-option --out-dir`. No source rebuild is required by the wrapper.

Analysis artifacts are kept under:

```text
<root>/<run-name>/analysis/<analysis-name>/
```

For analysis scripts that support the usual options, the wrapper injects:

```text
--out-dir <root>/<run-name>
--artifact-dir <root>/<run-name>/analysis/<script-stem>
--report-md <root>/<run-name>/analysis/<script-stem>/<script-stem>.md
```

Example:

```bash
tools/rediscover_run.py analyze 2026-06-02-broad-backbone-750 \
  src/rediscover/analysis/mcconnell_literature_decirc.py
```

## Substrate Vs Results

Substrate is the large emitted producer surface: per-run CSVs and NPY/NPZ
sidecars such as `*_sources.csv`, `*_aggregated.csv`, and emitted tensor
payloads. The cleaner identifies substrate from the producer
`<run>/manifest.json` and the runner's `<run>/.rediscover-run.json`.

Results are the small analysis outputs under `<run>/analysis/`: de-circ CSVs,
calibration JSON, charts, and generated markdown reports. Cleanup leaves these
in place.

The runner also keeps:

- `<run>/.rediscover-run.json` — wrapper metadata and cleanup history.
- `<run>/manifest.json` — producer substrate manifest.
- `<root>/manifest.json` — small index of managed runs.

## Keep Policy

Default substrate retention is:

```text
REDISCOVER_KEEP_SUBSTRATE=1
```

That means only the newest completed, non-recent run keeps its substrate by
default. Older managed runs keep their analysis results but lose superseded
substrate. Override per command with `--keep-substrate N`.

Analysis outputs are retained indefinitely by this framework.

## Drop-Old-On-Rerun

After a successful `emit`, the wrapper refreshes the root manifest and then
automatically drops superseded substrate using the keep policy. This is the only
automatic deletion path, and it only targets managed run directories under the
configured root.

Safety rules:

- Only directories with `.rediscover-run.json` are managed cleanup candidates.
- Only files recorded as substrate by the producer manifest or runner metadata
  are deleted.
- Before deleting any path, the cleaner verifies that the resolved path is
  under the managed run root, is not under `/shared`, and is not an
  `nmr_extract` extraction or inside one.
- `nmr_extract` extractions (16-hour atomic inputs: `trajectory.h5` plus
  `extraction_manifest.json`, `npys/frame_*`, and `pdbs/`) are categorically
  off the cleanup table; this framework manages only rediscover-engine
  substrate.
- Runs with status `running` are skipped.
- Runs touched within `REDISCOVER_ACTIVE_MINUTES` are skipped.
- Default active window: `720` minutes.
- Failed runs are skipped.
- The current run emitted by the command is protected by the keep policy and by
  its fresh mtime.

Existing scattered `/tmp` directories are not auto-adopted and are not touched.
This avoids disturbing a currently-running ad-hoc output directory.

## Manifest

Refresh and print the managed-run index:

```bash
tools/rediscover_run.py manifest
```

The index is written to:

```text
<root>/manifest.json
```

It records run name, timestamps, status, total size, substrate size,
analysis size, and whether substrate is still present.

## Clean

On-demand cleanup is dry-run by default:

```bash
tools/rediscover_run.py clean
```

Actually delete superseded substrate:

```bash
tools/rediscover_run.py clean --force
```

Keep two latest substrate sets instead of one:

```bash
tools/rediscover_run.py clean --keep-substrate 2 --force
```

For manually inspected managed runs that lack producer metadata, the cleaner can
fall back to top-level `*.csv`, `*.csv.gz`, `*.npy`, and `*.npz` files:

```bash
tools/rediscover_run.py clean --include-patterns
```

Use `--include-patterns --force` only after checking the dry-run output.
