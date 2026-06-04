> **Historical run record — not current truth (trued 2026-06-04).** Preserve as
> provenance; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Build1Fix Postmortem - 2026-06-03

Branch: h5-reader-pysr-spike; no switch, merge, or rebase performed.
Commit: `d9ba53dc0a7929eb4d4ce8185f6d9e6ff868c374`.
Run path: /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1

H1 fixed: CloudKind::ChargeSites distance/gap now exclude target atom,
zero/self distance, and same residue, matching the charge contribution source set.
Intentional re-baseline columns: nearest_charge_r, gap_to_2nd_charge_r,
bin_nearest_charge_r_4_6_8_10, bin_gap_to_2nd_charge_r_4_6_8_10.

M1 fixed: DftFrameSet has a scoped selection read guard.
Any DFT target read during CaseHunter selection throws runtime_error.
Added guard-fires regression test; normal CaseHunter path passes clean.

M2 fixed: conditioning unit labeling is explicit/_r-suffix based.
Counts, magnitudes, and fractions are not mislabeled A.
Unit metadata changed only for charge_excluded_same_residue_n,
abs_ring_jb_T2, and dominant_fraction_ring.

Build/test gates: cmake build of h5reader_extract and h5reader_rediscover_tests passed.
ctest -R h5reader_rediscover_tests passed after replacement.
Final emit: atoms=846, dft_frames=660, rows=558360.
CaseHunter anti_circular_assertion=true in final manifest.

Charge bin distribution after H1:
bin_nearest_charge_r_4_6_8_10 = {0:511367, 1:45091, 2:1899, 3:3}.
bin_gap_to_2nd_charge_r_4_6_8_10 = {0:273102, 1:133088, 2:144893, 3:7277}.
nearest_charge_r range = 1.2307903305559433..8.23339418172693 A.
gap_to_2nd_charge_r range = 1.2342593791458967e-07..1.3941541438454306 A.

Oracle parity: no missing or extra emitted files after metadata copy.
Changed files versus Build1 were exactly:
equations/charge/cases_manifest.csv, per_atom_substrate_column_specs.json,
per_atom_substrate_features_conditioning.npy, per_atom_substrate_manifest.json,
per_atom_substrate_partition_bins.npy.
All other sidecars, DFT targets, query results, and non-charge equation outputs were byte-identical.
All conditioning columns except nearest_charge_r and gap_to_2nd_charge_r were equal.
All partition columns except the two charge-distance bins were equal.
Rows, atoms, frame count, row alignment, R outputs, and uniqueness sidecars were unchanged.

Disk guard: df before emit had >=161G free on /tmp, above the 20G abort threshold.
Old Build1 was replaced only by explicit full paths under /tmp/rediscover-runs.
No /shared deletion or emitted-output write was performed.
Final /tmp/rediscover-runs size: 3.2G, below the 15G guard.
