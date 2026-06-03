# Codex brief — Piece 1: charge-scalars emit extension (`per_atom_substrate`)

Status: **landed** at `4bb9a0198594773c016160270b6a872b2f9b23fb`; retained as the execution brief.

You own the grind here; the lead vets and judges. This is the first build loop of
the all-atoms charge-complete fit. The stakes: the all-atoms fit's charge-source
positive control and its `all` feature tier both need two raw charge scalars on
the substrate. The change is small and additive — but it touches the **producer
substrate**, so **byte-parity against the v1 run is sacred**. Keep this run
SHORT: there is a known codex context-image platform bug on long runs; if you
crash mid-run, that is the platform bug, not your work, and the lead will re-fire
a shorter scope. Don't pad.

## Scope — Piece 1 ONLY

Implement ONLY "Piece 1: Small Emit Extension" from
`src/rediscover/ALLATOM_FIT_SPEC_2026-06-03.md`. Do **NOT** implement Piece 2 (the
Python fit / partition). Do **NOT** reshape, reorder, or reinterpret any existing
v1 output — this is purely additive source-coverage.

## What to add (`PerAtomSubstrate.{h,cpp}`)

Append four row-aligned columns to `per_atom_substrate_rows.csv`, AFTER all
existing v1 columns (appending preserves every existing column's ordinal meaning):

```
ff14sb_charge
ff14sb_charge_present
mopac_welford_mean_charge
mopac_welford_mean_charge_present
```

Exact source contract:

- `ff14sb_charge` ← `Catalog::value(ArrayId::Ff14sbCharge, atom, row)`; present ⇔
  `Catalog::present(ArrayId::Ff14sbCharge, atom, row)` AND value finite. Static
  atom-axis: repeat the same value on all 660 frame rows for that atom.
- `mopac_welford_mean_charge` ← `Catalog::value(ArrayId::MopacChargeWelfordMean,
  atom, row)`; present ⇔ `Catalog::present(...)` AND finite. Static atom-axis,
  repeated over the atom's 660 frame rows.
- **DO NOT use `ArrayId::MopacCharge`** (the catalog marks it absent — no
  per-frame MOPAC charge trajectory exists). The only MOPAC raw charge source for
  this chunk is the atom-axis Welford mean.
- Missing value ⇒ emit `NaN` with present flag `0`. Present flag is `0` or `1` only.

**First**, verify `ArrayId::Ff14sbCharge` and `ArrayId::MopacChargeWelfordMean`
exist in `Catalog.{h,cpp}` with the residence/units the spec table claims
(`ALLATOM_FIT_SPEC_2026-06-03.md:93-98`). If either is absent or named
differently, **STOP and report** — do not invent or rename an ArrayId.

Column-spec metadata: add 4 entries to `per_atom_substrate_column_specs.json`
exactly per the spec table (`:103-108`): array `per_atom_substrate_rows`; units
`e` for the two values, blank for the present flags; irreps `0e`; mechanism
`charges` for the values, `provenance_qc` for the flags; `is_feature` true for
values, false for flags; `native_axis` `rediscover_target_row`.

Manifest support block: add `ff14sb_charge_present_rows`,
`mopac_welford_mean_charge_present_rows`, and `charge_complete_rows` (rows with
finite `ff14sb_charge` AND `mopac_welford_mean_charge` AND `aimnet2_charge`).

SDK: per spec open-question Q2, `per_atom_substrate_column_specs.json` is the
source of truth for row-CSV column metadata. Add a `python/nmr_extract/_catalog.py`
loader description ONLY if genuinely required; do not reshape or replace existing
catalog entries. If unsure, leave `_catalog.py` untouched and flag it in the
postmortem.

## Reproduce v1 input + settings (so byte-parity is meaningful)

Derive the EXACT input calcset/`.LGS` path and config (ring/bond/charge cutoffs,
`mc_near_field_ratio`, align) from the v1 reference run's provenance:
`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-v1-fixed/.rediscover-run.json`
(+ its `per_atom_substrate_manifest.json`). Use the run framework
(`tools/rediscover_run.py`, `RUN_FRAMEWORK.md`): a fresh standard root, drop-old,
and the **nmr_extract guard MUST be active** (it refuses `/shared` paths and the
extraction signature — `/shared` extractions are SACRED, never write/delete
there). Emit the new charge-complete run to that fresh root.

The v1 dir `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-v1-fixed` is the
**read-only parity baseline** — never modify or delete it.

Build (you are unsandboxed): `cmake --build build/linux-gcc --target
h5reader_extract h5reader_rediscover_tests`. Run the emit via the
`per_atom_substrate` case.

## Emit-extension gates — ALL must pass before any commit

(From the spec's "Emit-Extension Checks". Run each; report each.)

1. **R & uniqueness** — rows == n_atoms·n_dft_frames; values 846 / 660 / 558,360;
   `row_id` dense and == `frame_slot*n_atoms+atom_index`; `(atom_index,
   original_frame_index)` unique; every NPY sidecar first dim == row count.
2. **Charge scalar coverage** — new cols exist; present flags 0/1; present⇒finite,
   absent⇒NaN; per atom, `ff14sb_charge` constant over its 660 present frame rows;
   same for `mopac_welford_mean_charge`; manifest support block reports the 3 counts.
3. **Existing output UNCHANGED (parity)** — the 7 pre-existing v1 NPY sidecars
   (`target_T2`, `target_T0`, `features_classical`, `features_conditioning`,
   `driver_modulation_by_atom`, `backbone_audit`, `aimnet2_embedding`) must be
   **byte-identical** to the v1 run — they are not rewritten by this change, so
   `cmp`/sha256 each (cheap and load-bearing; this is where the real rigor lives).
   For `rows.csv`, do NOT attempt raw byte-identity (appended columns lengthen
   every line by construction): parse BOTH CSVs and assert the existing columns
   are **value-identical row-for-row**, and that the only delta is the 4 new
   columns, in order, at the end. `ring_identity.csv` + default `query_results/`
   unchanged. DFT-target parity exact under the existing frame-alignment
   diagnostic. Determinism note: if any pre-existing NPY differs, FIRST confirm
   the v1 input + settings were reproduced exactly — a non-substantive diff almost
   always means a settings/path mismatch, not a code defect — before treating it
   as a perturbation.
4. **Backbone regression UNCHANGED** — the 6 backbone strata (N/CA/C/O/HN/HA),
   ignoring the 4 new cols, reproduce the broad-backbone compatibility checks;
   `backbone_audit.npy` matches; the new columns perturb no old reducer's held-out
   scores beyond floating-point tolerance.
5. **Charge positive-control readiness** — `formal_charge`, `ff14sb_charge`,
   `mopac_welford_mean_charge`, `aimnet2_charge` all evaluable from the one
   emitted substrate (no H5 read, external table, or multi-dir join).

If ANY gate fails: STOP, do NOT commit, report exactly which gate and the diff. A
failed parity gate means the emit perturbed existing output — a real defect, not a
tolerance issue.

## Commit (only if all gates green)

Branch is `h5-reader-pysr-spike`. **NEVER merge / switch / rebase / PR / checkout
another branch.** One atomic commit: the `PerAtomSubstrate` emit change + the
column-spec/manifest additions (+ a `_catalog.py` loader line only if required).
Message: describe the additive charge-scalars emit and that v1 byte-parity held.
Do not commit emitted run data (it lives in the drop-old run root).

## Report back — drift-assessment + postmortem

End with a short, plain section (no editorializing — the lead judges):
- What landed vs what the Piece-1 spec asked. Any drift?
- Each gate's result (pass/fail + key numbers: both `*_present_rows`,
  `charge_complete_rows`, parity byte-diff = 0).
- Commit hash if committed, or the blocking failure.
- Anything you had to decide that the spec left open.
