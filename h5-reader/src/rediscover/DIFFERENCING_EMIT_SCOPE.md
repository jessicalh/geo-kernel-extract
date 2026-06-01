# Differencing Emit Scope (#35)

Goal: support frame-to-frame differencing as a system-ID probe: delta shielding
versus delta geometry per mechanism, frame-to-frame. This is a Python-side data
transform over emitted columns, not a physics recompute. Python may subtract
emitted values across frames, but must not rebuild kernels, projections, fields,
or source contributions that the C++ substrate did not emit.

## Conclusion

No required C++ emit addition for the first differencing analysis.

The current broad-backbone substrate already emits the pieces needed to build
consecutive DFT-frame pairs per target atom, join source rows back to their
target frame, and align physical sources across adjacent frames:

- `broad_backbone_aggregated.csv`: one `(atom, frame)` row with `row_id`,
  `atom_index`, `h5_row`, `original_index`, `time_ps`, emitted aggregate
  mechanism features, and the DFT target.
- `broad_backbone_sources.csv`: one `(atom, frame, source)` row with `row_id`,
  target atom identity, mechanism, emitted source geometry, and stable physical
  source identity columns.

Any pair index, delta time, source-pair key, churn flag, or finite-difference
column should be derived in Python from this substrate.

## Frame Pairing

Differencing needs consecutive DFT-frame pairs per target atom. The current
aggregated emit is sufficient:

1. Load `broad_backbone_aggregated.csv`.
2. Keep rows with `dft_present == 1` and `frame_valid == 1`.
3. Group by `atom_index`.
4. Sort each group by `h5_row` (or `time_ps`; `h5_row` is the extractor's sorted
   DFT-row traversal key).
5. Pair adjacent rows in that sorted DFT subset.
6. Record derived `row_id_t`, `row_id_tp1`, `h5_row_t`, `h5_row_tp1`,
   `original_index_t`, `original_index_tp1`, and `dt_ps`.

These are consecutive DFT samples, not necessarily adjacent H5 rows. `time_ps`
already gives the actual interval if the analysis wants a rate
`delta / dt_ps` rather than an unscaled finite difference.

No C++ pair table is required. A future additive `dft_frame_ordinal` on the
aggregated row would be an ergonomic convenience only, not a blocker.

## Source Alignment

Per-source frame-to-frame alignment is possible from the current source emit
after joining source rows to the aggregated row on `row_id`.

Use stable physical keys, never source-row order:

| mechanism | stable source key | current emit status |
|---|---|---|
| ring | `(atom_index, "ring", ring_index)` | `ring_index` is the absolute topology ring id; stable across frames. |
| bond | `(atom_index, "bond", bond_index)` | `bond_index` is emitted; `bond_atom_a/b` are also emitted as a backup provenance key. |
| charge | `(atom_index, "charge", source_atom_index)` plus `charge_source` from the joined aggregate row | `source_atom_index` is the absolute source atom id; stable across frames. |

KD membership can churn at the cutoff boundary. That does not make identities
unalignable; it means the Python analysis must classify each adjacent pair:

- present in both frames: compute per-source finite differences of emitted
  source geometry columns.
- present only at `t` or only at `t+1`: count/report source churn separately,
  and do not silently interpret cutoff entry/exit as a smooth geometric delta.
- aggregate mechanism columns may still be differenced directly; those deltas
  include any membership churn by construction.

For frozen-membership ring rows, align by `ring_index`, not by slot position.
For broad-backbone KD ring rows, `ring_index` is also the stable physical key.

## What To Difference

Target:

- scalar: `dft_sigma_iso`
- tensor option: emitted DFT T2 components (`dft_total_T2_*`) or the target-local
  T2 sidecar, keyed by aggregated row order, if the tensor analysis is used

Aggregated mechanism features:

- ring: `ring_sum_dipolar`, `ring_n`, emitted ring literature-kernel T2 columns
- bond: `bond_sum_dipolar`, `bond_n`, emitted bond literature-kernel T2 columns
- charge: `field_local_*`, `field_z`, `field_mag`, `mu_local_*`, emitted charge
  literature-kernel T2 columns
- combined: emitted `literature_kernel_T2_*`

Per-source geometry:

- common: `disp_local_x/y/z`, `r`, `cos_theta`, `dipolar`
- ring: `source_normal_local_x/y/z`, `ring_index`, ring type/intensity metadata
- bond: `bond_axis_local_x/y/z`, `bond_index`, bond endpoints/category metadata
- charge: `disp_local_x/y/z`, `r`, `source_q_e`, `source_atom_index`

Subtracting any of these emitted values across aligned frames is allowed. Do not
derive an un-emitted field, EFG, projection, or kernel in Python. If a later
analysis specifically needs per-source charge field or EFG contribution deltas,
that would be a separate C++ emit request; it is not needed for the first
frame-to-frame geometry differencing probe.

## Analysis Plan

Build two Python tables:

1. `atom_pairs`: adjacent DFT-frame pairs from the aggregated rows, with target
   deltas and aggregate mechanism deltas.
2. `source_pairs`: source rows joined to `atom_pairs` by `(row_id_t,
   row_id_tp1)` and the stable mechanism-specific source key, with per-source
   geometry deltas and churn flags.

Then evaluate:

- delta target versus aggregate delta mechanism features, per backbone stratum.
- delta target versus aligned per-source delta geometry, grouped by mechanism.
- same analysis for windows `w = 1, 3, 5` DFT samples if the lag-1 gate says
  single-step differences are too noisy.
- source churn rates per mechanism/cutoff, so cutoff-boundary effects are visible.

The C++ `WindowFn`/delta combinator remains the right producer-side abstraction
if differencing is ever promoted into the engine, but #35 does not need that.

## Autocorrelation Guard

Consecutive DFT samples are correlated, so differencing needs the same gate
started in `analysis/diag_differencing.py`:

- lag-1 autocorrelation of the target per atom and of emitted source/mechanism
  feature series.
- noise cost estimate from `Var(diff) / (2 Var(level)) ~= 1 - autocorr`.
- decorrelation check: compare mechanism/source feature collinearity in levels
  versus differences, including `w = 1, 3, 5` windows.
- report effective sample size or, at minimum, flag thin/highly correlated strata.

Differencing is a technique for probing response. It is not source identity, and
it should not be accepted as helpful unless the smoothness and decorrelation
checks pass.

## Scientific Framing

This is a geometric-perturbation probe: how emitted shielding targets change as
the emitted geometry changes between nearby sampled conformations. It is not a
claim about relaxation dynamics or causal time evolution. The trajectory supplies
nearby geometries; the system-ID claim is the local response relation
`delta geometry -> delta shielding`.

## Emit Catalog Addition

Required emit changes: none.

The follow-up implementation should consume the existing broad-backbone
aggregated/source CSVs and sidecars, then construct pair and source-delta tables
in Python. If ergonomics later justify an additive helper, the only reasonable
candidate is `dft_frame_ordinal` on aggregated rows; it must not duplicate target
data onto source rows or require C++ to compute finite differences.

## Open Questions

- Should the first report use scalar `dft_sigma_iso`, local T2, or both? The
  substrate supports both, but the scalar probe is the smallest first gate.
- How much KD source churn occurs at the current cutoffs, and does it dominate
  aggregate deltas for any mechanism?
- For charge, are emitted aggregate field deltas enough for the first mechanism
  probe, or does a later per-source charge-field contribution emit become useful?
