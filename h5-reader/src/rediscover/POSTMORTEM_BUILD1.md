# Build 1 postmortem

> **Historical run record — not current truth (trued 2026-06-04).** Preserve as
> provenance; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Branch: h5-reader-pysr-spike.
Run dir: /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1.
Commit hash: `a3531040655955a4d393e314c5a69de51b1c54bf`.

Disk guard: PASS. Free was >=160G after drop-old before primary emit; final free 156G.
Rediscover output: 3.2G total, 3.1G Build1, under the 15G budget.
Deleted exact paths: /tmp/rediscover-runs/2026-06-03-per-atom-substrate-piece3b-final.
Deleted exact paths: /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1-detcheck.
Deleted exact paths: /tmp/rediscover-runs/2026-06-03-per-atom-substrate-oracle-head.
Deleted exact paths: /tmp/oracle-head-src and /tmp/oracle-head-build.
No write/delete under /shared.

Build/test: PASS, h5reader_extract + h5reader_rediscover_tests; ctest passed.
Oracle parity: PASS. 26 old NPY/CSV files byte-identical.
Append-only conditioning: PASS, old (558360,26) equals new first 26 cols of (558360,32).
Backbone audit and DFT T0/T1/T2 targets: PASS, byte-identical to HEAD oracle emit.
Dominance two-path: PASS; max abs diffs ring/charge/mc = 4.995e-13/4.996e-13/4.996e-13.
Hunter: PASS; candidates ring=24, charge=24, mc=24.
Hunter deterministic: PASS; 38 emitted files matched byte-for-byte in detcheck.
Anti-circularity: PASS; selection predicate is typed_habitat && isolation && motion && quiet.
DFT recovery is measured only as dft_recovery_R2 in cases_manifest.csv after selection.
R/uniqueness: PASS; 846 atoms x 660 DFT rows = 558360, dense row_id, unique atom/frame.

New conditioning ranges:
gap_ring 1.919e-07..10.0467; gap_charge 0.959992..1.91988; gap_bond 1.621e-07..5.08524.
dom_ring 0.208778..1; dom_charge 0.0398653..1; dom_mc 0.0627701..0.924099.
Partition bins: PASS; shape (558360,25), values -1..4, 25 manifest columns.

C++ spine: PASS. Emitted values are computed in PerAtomSubstrate/CaseHunter C++; Python was validation only.
query_frame_slots remains 1; no per-source full-frame dump added.
