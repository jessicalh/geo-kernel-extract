# Codex brief — Build 1: the C++ partition-filter tool (isolation primitives + bin-ids + CaseHunter)

> **Historical run brief — not current truth (trued 2026-06-04).** Session
> provenance only; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Status: **landed** at `a3531040655955a4d393e314c5a69de51b1c54bf`; fixed by Build1Fix.

You own the grind; the owner vets. This is the C++ half of the PARTITION step — the
**filters live in C++ over the resident indexes**, not in Python. Design reference:
the partition-filter design just produced (the owner does NOT lean Python — keep selection/
filtering/binning in C++; Python is only the later ridge/CV fit). De-phased deliberately:
this is ONE build, not six. Read `PerAtomSubstrate.cpp` (the existing emit + the
`rowPairContributions`/`dominance` machinery) + `SpatialIndexSet`/`TemporalIndex`/
`TypedAtomIndex` headers first.

## HARD RULES

- **C++ spine ONLY. Zero Python computes/derives any emitted value.** Mirror the existing
  `PerAtomSubstrate.cpp` mechanism/conditioner blocks (typed, no string dispatch, ArraySpec
  per output, present-flags, metadata). Producer/CMake/ctest untouched beyond this rediscover
  code. Read existing H5/NPY + the live model/indexes only; no re-extraction.
- **DISK GUARD (don't crash the machine — single filesystem, ORCA resident + writing):**
  run `df` before any emit; if free < **20 G**, ABORT and report (do NOT write into a
  near-full disk). **drop-old / clean up as you go:** the new substrate emit REPLACES the
  prior `2026-06-03-per-atom-substrate-piece3b-final`; you MAY also remove your OWN
  no-longer-needed intermediate / superseded run dirs to reach the result and stay under
  budget. Total rediscover output stays **< 15 G**. **When you delete ANYTHING, use the
  explicit FULL PATH — NEVER a regex/glob** (no `rm /tmp/rdc-*`-style wildcards; name each
  directory you remove, and print it). NEVER write or delete under `/shared`.
- **Lean / keyed (atom,frame) only.** No per-source/pairwise dump. Critically: do **NOT**
  raise `query_frame_slots` to 660 to get full-frame dominance — that re-emits the heavy
  pairs CSV (the boa-constrictor). Compute dominant-fraction as a SEPARATE lean scalar
  reduction emitting only the per-(atom,frame) scalar.

## Piece A — isolation primitives (in `conditioningFeatures`, `PerAtomSubstrate.cpp:1703`)

Add, row-aligned per-(atom,frame), append-only:
- `gap_to_2nd_{ring,charge,bond}_r` — sort the first two distances from the existing KD
  `near()` result (`nearestDistance`, `:1709`); gap = r(2nd) − r(1st).
- `dominant_fraction_{ring,charge,mc}` — `max|contribution| / Σ|contribution|` over that
  mechanism's pairs, reusing the existing `dominance` lambda (`:2755`) over
  `rowPairContributions`, as a lean per-(atom,frame) scalar (NOT via the heavy named query).

## Piece B — bin-id columns (C++ binning, so Python never derives bins)

Emit per-(atom,frame) bin-id columns for the predeclared condition families (the
`ALLATOM_FIT_SPEC` families: geometry distance bands 4/6/8/10 Å fixed-edge; driver-magnitude
& modulation quintiles — compute the quintile edges C++-side in one pass, record edges in the
manifest). Python will bin by lookup, not by deriving.

## Piece C — `CaseHunter` (new `CaseHunter.{h,cpp}`)

Per mechanism: candidate set = typed habitat (`TypedAtomIndex::select`) × frame-windows
(`TemporalIndex::range`). Per candidate, over `rowPairContributions` across the window,
compute the THREE input-side cleanliness metrics: isolation (`dominant_fraction ≥ θ_dom`
AND `gap_to_2nd ≥ θ_gap`), motion (within-window variance of the mechanism's input driver),
quiet-confounders (within-window variance of the summed OTHER mechanisms); score =
var(driver)/var(others). Rank; **MEASURE** DFT recovery per candidate (the curve y-axis) but
**NEVER select on it** (anti-circular — selection predicate touches input geometry only).
Deterministic (same thresholds → identical output).

## Piece D — cases manifest (the hunter's named product)

Emit `equations/<mechanism>/cases_manifest.csv` (or json): top-N navigable addresses
`(protein, atom_index, frame_window_begin, frame_window_end, mechanism, strip_metric_ids)` +
the per-case input-side cleanliness metrics + the measured DFT recovery. Pointers into the
same (atom,frame) addressing as the existing pair index. This is what the reader strips load.

## Gates (ALL pass before commit)

- **Disk guard** passed (free ≥ 20 G at emit; drop-old; total < 15 G).
- **Oracle parity:** every pre-existing NPY/CSV byte-identical (append-only); `backbone_audit`
  byte-identical; DFT target parity exact.
- **Dominance two-path cross-check:** the new per-(atom,frame) `dominant_fraction` equals the
  existing 1-frame named-query `dominance` value on frame-slot 0 (the two paths agree).
- **Hunter:** produces > 0 candidates for at least the ring and charge habitats; deterministic
  re-run identical; anti-circularity assertion (no DFT field read inside the selection predicate).
- R & uniqueness (846/660/558,360; dense row_id; new columns finite where present).

## Build / commit / output

`cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests`
(unsandboxed). Branch `h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR/checkout. One
atomic commit. Write `src/rediscover/POSTMORTEM_BUILD1.md` (≤50 lines): disk-guard result
(free before/after), each gate pass/fail, the dominance two-path agreement, per-mechanism
candidate counts, new-column ranges, commit hash, run dir. Print ONLY a ≤10-line summary +
that path; DO NOT echo diffs to stdout.
