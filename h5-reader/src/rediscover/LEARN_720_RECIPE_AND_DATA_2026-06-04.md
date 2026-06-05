# Learn 720 Recipe And Data - 2026-06-04

Read-only spinner pass over:

```text
/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn
```

Built on the earlier survey:

```text
/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/LEARN_STAGE1_UNDERSTANDING_2026-06-04.md
```

No extraction, build, copy, or git action was run.  Only this report was
written.

## Executive Answer

**Recipe found:** yes.

The original 720 producer extraction recipe is the `learn/extract.py` batch
runner, documented in `learn/docs/extraction_procedure.md`.  It invokes
`build/nmr_extract --mutant` once per WT/ALA calibration pair, using the
default calculator config unless overridden.  For the Stage-1 BMRB/topology
runs, `--stage1-audit` is the important flag: it requires WT and ALA ORCA
`*_nmr.out` files and validates the Stage-1/sidecar outputs after each protein.

**Old data usable:** yes, with one important boundary.

The complete frozen per-protein NPY substrate is present at:

```text
/shared/2026Thesis/nmr-shielding/calibration/features/Stage1BMRB_20260513_topology
```

It has 720 complete protein directories, 100 NPY arrays per complete protein,
MOPAC arrays, topology sidecars, positions, ORCA absolute shielding arrays, and
WT-ALA delta arrays.  It is stale relative to the current producer/schema work,
but it is not missing the data needed for first-pass Model-1 training.

Do **not** use the old topology workup's `feature_matrix.npy` /
`target_matrix.npy` as the Model-1 all-atom absolute-shielding corpus.  That
matrix export is Stage-1 mutation-delta/matched-row material:
425,599 matched atom rows, delta target columns, and mutation-site residues
excluded.  Use the raw per-protein NPY directories, or a new static ingest
adapter over them.

**Recommendation:** use the old frozen per-protein data and skip the producer
re-run for the first combined 1P9J + 720 Model-1 corpus.  Re-run only if the
lead requires exact current producer schema/.LGS manifests/T2-only EFG arrays
with no adapter layer.

## Recipe: How To Re-run The 720 Extraction

### Entry Point

Canonical runner:

```text
/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/extract.py
```

Procedure doc:

```text
/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/docs/extraction_procedure.md
```

From repo root:

```bash
cd /shared/2026Thesis/nmr-shielding
cmake --build build --target nmr_extract -j$(nproc)

python3 learn/extract.py \
  --run NewRunName \
  --resume \
  --stage1-audit
```

The Stage-1 topology run's log did not echo the exact shell command, but its
record matches that form:

```text
WARNING: 'Stage1BMRB_20260513_topology' not in known runs: [...]
Proteins: 723  Jobs: 722  Skipped: 1
...
Done: 719 OK, 3 failed, 28742s total
```

The actual run had one pre-existing skipped output, 719 new OK jobs, and 3
expected missing-ORCA-input failures.  The final validated directory count is
720 OK.

### What `extract.py` Runs

For each protein, `extract.py` constructs:

```text
build/nmr_extract --mutant \
  --wt  calibration/{ID}/{ID}_WT \
  --ala calibration/{ID}/{ID}_ALA \
  --output calibration/features/{run}/{ID} \
  --config data/calculator_params.toml
```

With `--stage1-audit`, it also:

- requires WT and ALA ORCA `*_nmr.out` inputs before starting the job;
- validates required Stage-1 arrays after extraction;
- validates topology sidecar files added in the 2026-05-13 topology pass.

Relevant required outputs in `learn/extract.py` include:

```text
atoms_category_info.npy
delta_shielding_diamagnetic.npy
delta_shielding_paramagnetic.npy
wt_shielding_diamagnetic.npy
wt_shielding_paramagnetic.npy
mut_shielding_diamagnetic.npy
mut_shielding_paramagnetic.npy
aimnet2_polarisability.npy
aimnet2_polarisability_scalar.npy
pyramidalization.npy
omega_actual.npy
omega_deviation.npy
omega_is_xpro.npy
aromatic_chi2.npy
pucker_Q.npy
pucker_theta.npy
residues.npy
bonds.npy
rings.npy
ring_membership.npy
extraction_manifest.json
```

### Inputs

Calibration pair root:

```text
/shared/2026Thesis/nmr-shielding/calibration/{ID}
```

Expected files per pair:

```text
{ID}_WT.xyz
{ID}_WT.prmtop
{ID}_WT_nmr.out
{ID}_ALA.xyz
{ID}_ALA.prmtop
{ID}_ALA_nmr.out
```

The 720 set is the 723 calibration pairs minus three expected missing ORCA
shielding inputs:

```text
A0A075FQU3  missing WT nmr.out
A0A7C4ZM98  missing ALA nmr.out
A0A7J2L4W1  missing WT nmr.out
```

### MOPAC State

MOPAC is ON for this recipe when using the default calculator config:

```text
/shared/2026Thesis/nmr-shielding/data/calculator_params.toml
```

The extraction procedure doc says MOPAC dominates runtime and describes
PM7/MOZYME as the expensive leg.  The actual 2026-05-13/14 topology run wrote
MOPAC arrays for all 720 complete proteins:

```text
mopac_charges.npy
mopac_bond_orders.npy
mopac_coulomb_E.npy
mopac_coulomb_efg_aromatic.npy
mopac_coulomb_efg_backbone.npy
mopac_coulomb_scalars.npy
mopac_coulomb_shielding.npy
mopac_global.npy
mopac_mc_category_T2.npy
mopac_mc_scalars.npy
mopac_mc_shielding.npy
mopac_scalars.npy
```

All of those files are present in 720/720 complete directories.

### Output Layout

Per-run root:

```text
/shared/2026Thesis/nmr-shielding/calibration/features/{run}
```

Per-protein output:

```text
/shared/2026Thesis/nmr-shielding/calibration/features/{run}/{ID}
```

Run-level logs:

```text
extract_log.jsonl
extract_background.log
extract_background.pid
```

Per-protein outputs are flat NPY/JSON directories.  A representative complete
directory contains 100 NPY arrays plus `extraction_manifest.json`.

## Old 720 Data State

### Raw Per-protein Producer Output

Canonical current clean 720 feature run:

```text
/shared/2026Thesis/nmr-shielding/calibration/features/Stage1BMRB_20260513_topology
```

Dates:

```text
feature root mtime:       2026-05-14 08:07:06 +0100
extract log mtime:        2026-05-14 08:08:11 +0100
first sampled manifest:   generated_at_utc 2026-05-13T23:07:47Z
last sampled output:      2026-05-14 08:08:10 +0100
```

No files under that feature root were newer than 2026-05-15 in the metadata
check.  Treat it as frozen.

Completeness checks:

```text
720 protein directories
720 "ok": true lines in extract_log.jsonl
3 "ok": false lines, all expected missing ORCA inputs
720 directories with 100 NPY files each
720 pos.npy
720 delta_shielding.npy
720 delta_shielding_diamagnetic.npy
720 delta_shielding_paramagnetic.npy
720 ring_contributions.npy
720 ring_geometry.npy
720 atoms_category_info.npy
720 residues.npy
720 bonds.npy
720 rings.npy
720 ring_membership.npy
720 extraction_manifest.json
720 MOPAC core/Coulomb/McConnell arrays
720 ORCA absolute arrays: orca_total.npy, orca_diamagnetic.npy, orca_paramagnetic.npy
```

Representative sample from `A0A062V9G2`:

```text
pos.npy:                         (732, 3)
orca_total.npy:                  (732, 9)
orca_diamagnetic.npy:            (732, 9)
orca_paramagnetic.npy:           (732, 9)
delta_shielding.npy:             (732, 9)
ring_contributions.npy:          (3057, 59)
ring_geometry.npy:               (7, 10)
mopac_charges.npy:               (732,)
mopac_bond_orders.npy:           (714, 3)
mopac_coulomb_E.npy:             (732, 3)
mopac_coulomb_efg_aromatic.npy:  (732, 9)
mopac_coulomb_efg_backbone.npy:  (732, 9)
mopac_mc_category_T2.npy:        (732, 25)
mopac_mc_shielding.npy:          (732, 9)
mopac_scalars.npy:               (732, 4)
residues.npy:                    48 structured rows
bonds.npy:                       739 structured rows
rings.npy:                       9 structured rows
ring_membership.npy:             52 structured rows
```

This is the data to use if we skip the re-run.

### Derived Topology Workup Output

Current workup root:

```text
/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage1-topology-workup-20260514
```

Main compendium:

```text
/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage1-topology-workup-20260514/docs/STAGE1_COMPENDIUM.md
```

Validated table/matrix export:

```text
/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage1-topology-workup-20260514/derived/tables
```

Important rows:

```text
run_validation.csv:       720 OK directories
feature_catalog.csv:      109 feature-catalog rows
feature_matrix.npy:       (425599, 788) float32
target_matrix.npy:        (425599, 23) float32
atom_context.csv:         matched WT atom rows
atom_context_all.csv:     all topology atom rows
```

The `matrix_export_summary.json` reports:

```text
proteins:           720
matched_atom_rows:  425,599
feature_columns:    788
target_columns:     23
```

The `STAGE1_COMPENDIUM.md` reports the wider row contract:

```text
720 extractor output directories validated as OK
425,599 matched atom rows
788 feature columns
23 target columns
475,116 total topology atom rows
49,517 unmatched mutation-residue atoms
172 reproduction-audit checks, 0 failures
```

Boundary for Model-1:

- `feature_matrix.npy` / `target_matrix.npy` are good evidence that the raw
  producer output is readable and feature-rich.
- They are **not** the desired all-atom absolute-shielding Model-1 corpus.
- `target_columns.csv` shows all 23 targets come from `delta_shielding`.
- `export_matrix.py` explicitly keeps only matched WT rows because those are
  the Stage-1 shielding-delta target rows.
- `export_context.py` keeps all atoms separately because the matched matrix
  excludes mutation-site residue atoms.

For Model-1, the useful substrate is the raw per-protein NPY data, especially:

```text
pos.npy
atoms_category_info.npy
residues.npy
bonds.npy
rings.npy
ring_membership.npy
orca_total.npy
orca_diamagnetic.npy
orca_paramagnetic.npy
wt_shielding_diamagnetic.npy
wt_shielding_paramagnetic.npy
all MOPAC arrays
all kernel/shadow arrays
```

## How Stale Is It?

The old data is frozen at the 2026-05-13/14 topology-sidecar producer state.
As of 2026-06-04, it is about three weeks old.

Staleness is real but mostly schema/producer-contract staleness, not absence of
the science payload:

- The raw data has old 9-component 3x3 tensor arrays for EFG/shielding-style
  outputs.  The current repo later changed EFG schema: commit `b580ebe`
  on 2026-05-18 says "EFG schema rev: all 5 calculators drop
  structurally-zero T0+T1 channels."  A current strict loader may therefore
  expect T2-only EFG surfaces or current metadata.
- AIMNet2 naming/metadata evolved after this run.  The old data has
  `aimnet2_polarisability*`; current docs/code discuss AIMNet2 charge response
  gradient naming.
- The current rediscover/static-ingest design now wants a lean per-(protein,
  atom) static substrate, `.LGS`/manifest discipline, exact required-field
  loading, T2 orientation handling, and grouped sidecars.  That was designed
  after this Stage-1 topology run.
- The current static ingest wants all atoms and absolute WT shielding target.
  The old raw NPYs have those.  The old Stage-1 derived matrix does not.
- The current repo has many later trajectory/result and comment-truthfulness
  changes.  Those do not by themselves prove the May 14 static WT data is
  scientifically invalid.

Practical interpretation:

- **Usable as raw data:** yes.
- **Drop-in current-schema producer output:** no.
- **Drop-in Stage-3 all-atom Model-1 matrix:** no.
- **Good enough to avoid a costly MOPAC producer re-run for first-pass
  training:** yes, if the static ingest/exporter is allowed to adapt old NPY
  schema to the current Model-1 substrate contract.

## Recommendation

Use the old frozen per-protein NPY data and skip the 720 producer re-run for
the first combined 1P9J + 720 Stage-3 Model-1 corpus.

Reasoning:

1. The costly substrate is present and complete: 720/720 validated directories,
   all required topology sidecars, positions, MOPAC fields, kernel arrays, and
   absolute ORCA shielding arrays.
2. The missing 3 proteins are not extraction failures; they lack required ORCA
   shielding inputs, and all docs already treat the clean set as 720.
3. A re-run would mainly buy schema freshness, not new DFT targets or an
   obviously missing MOPAC payload.
4. The existing Stage-1 matrix export proves the producer output can be read
   reproducibly, but it should not be reused as the training target because it
   is mutation-delta/matched-row data.
5. The right low-cost path is a current static ingest/export pass over the old
   per-protein NPY directories: all atoms, absolute `orca_total` or vetted
   WT dia/para/total target, current feature naming, current T2 conversion, and
   explicit old-schema provenance in the emitted manifest.

Tradeoff:

- **Skip re-run:** saves producer/MOPAC time and avoids touching the expensive
  archive; accepts adapter/schema-conversion risk.
- **Re-run:** cleaner current-schema output and fewer adapter caveats; costs a
  full 720 producer pass, creates large new data, and still uses the same old
  r2SCAN ORCA target files.

I would re-run only if one of these is a hard requirement:

- the static ingest must consume only current `.LGS`/single-pose producer
  manifests and exact current field catalog metadata;
- the lead refuses an old-schema adapter/converter;
- a later audit finds a real calculator bug in the 2026-05-13/14 arrays, not
  just metadata drift;
- the Model-1 substrate must include new post-2026-05-14 fields that are not
  reconstructible from the old NPYs.

## Tight Summary

Recipe found: **yes** -
`/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/extract.py` plus
`/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/docs/extraction_procedure.md`;
run from `/shared/2026Thesis/nmr-shielding` with
`python3 learn/extract.py --run <new-run> --resume --stage1-audit`, which
drives `build/nmr_extract --mutant` over the 723 WT/ALA calibration pairs with
MOPAC ON under the default config.

Old data usable: **yes, raw per-protein NPYs are usable** -
`/shared/2026Thesis/nmr-shielding/calibration/features/Stage1BMRB_20260513_topology`
has 720 complete frozen outputs from 2026-05-13/14, including MOPAC fields and
absolute ORCA DFT shielding arrays.  It is about three weeks stale and predates
current EFG/schema/static-ingest conventions; the old derived Stage-1
`target_matrix.npy` is **not** the Model-1 target because it is matched
mutation-delta data.

Recommendation: **skip the producer re-run for first-pass Model-1 training**;
consume the frozen per-protein NPYs through a current static ingest/exporter
that records old-schema provenance and converts to the current all-atom,
absolute-shielding substrate.  Re-run only for exact current producer schema or
if a real calculator bug is found.
