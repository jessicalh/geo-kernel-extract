# Codex brief — Build 1 fix loop (H1 + M1 + M2 from the adversarial review)

Status: **landed** at `d9ba53dc0a7929eb4d4ce8185f6d9e6ff868c374`; retained as the fix-loop brief.

You own the grind; the owner vets. Three targeted fixes to the just-landed partition-filter
code (commit `a353104`), from the adversarial review, BEFORE the Build-2 fit/partition
consumes it. C++ spine only — zero Python computes/derives any emitted value. This re-emits
the substrate (the charge-distance columns intentionally change — the canonical re-run
pattern: re-baseline the affected columns, byte-parity everything else).

## H1 (HIGH) — charge cloud self-inclusion makes the charge isolation/distance columns degenerate

Root cause: the self-exclusion branch in `gapToSecondDistance` / `nearestDistance`
(`PerAtomSubstrate.cpp:178-208`) only fires for `CloudKind::Atoms`; the `ChargeSites` cloud
includes the target atom itself (and same-residue atoms), so for a charge query the nearest
"source" is the atom at r=0 → `nearest_charge_r` collapses, `gap_to_2nd_charge_r` is pinned
~1 Å, and `bin_nearest_charge_r`/`bin_gap_to_2nd_charge_r` are constant-0 (zero-information).

Fix: for `CloudKind::ChargeSites`, exclude the target atom AND same-residue atoms in BOTH
`gapToSecondDistance` and `nearestDistance` — match the dominance path's neighbour set
(`makeChargeContribution` rejects `sourceAtom == targetAtom` and same-residue,
`PerAtomSubstrate.cpp:622-624`). Generalize the self-test (skip a source whose `entity_index`
== the target atom, or `r < epsilon`) rather than hard-coding Atoms-only.

This re-values `nearest_charge_r`, `gap_to_2nd_charge_r`, and the dependent charge bins —
INTENTIONAL. They are supposed to change; re-baseline them.

## M1 (MEDIUM) — make the anti-circularity guard REAL (it is currently a vacuous no-op)

`CaseHunter.cpp:343` declares `bool dftTouchedDuringSelection = false;` and never sets it
true, so the throw at `:377-378` is dead and `anti_circular_assertion: true` is a hardcoded
green light that can't go red. The invariant holds today by construction — but the guard
protecting the program's load-bearing rule must be able to FIRE.

Fix: make it a genuine guard — during the CaseHunter SELECTION phase (`evaluateInputCandidate`
+ ranking), instrument the DFT-target read path so ANY read of `run.dft` / the DFT target
(orca total/T0/T1/T2/dia/para) THROWS (or sets the flag that throws). Cleanest is a scoped
"in-selection" mode on the `Body`/DFT accessor that asserts-on-DFT-read; pick what's least
invasive. **REQUIRED: add a test that DELIBERATELY reads the DFT target inside the selection
path and confirms the guard THROWS** — prove it can go red. The normal run must still pass
(selection touches only geometry / FF14SB charge / literature kernels, never DFT). After this,
`anti_circular_assertion: true` means the real guard ran and stayed clean.

## M2 (MEDIUM) — units mislabel from over-broad `"_r"` substring

The column-spec writer (`PerAtomSubstrate.cpp` diff line ~493) uses
`name.contains("_r") ? "A" : ""`, which matches "_ring"/"_residue" etc. — mislabeling
`charge_excluded_same_residue_n` (count), `abs_ring_jb_T2` (tensor magnitude),
`dominant_fraction_ring` (dimensionless) as Å. Fix: match the actual distance suffix
(`name.endsWith("_r")`) and/or set units explicitly per column. Only true distance columns
carry `units = A`.

## L1/L2 (LOW, one-liners) — optional but cheap

Note in the manifest/doc: the selection "motion" window covers all window rows while
`dft_recovery_R2` is over the DFT-present rows only (different frame sets, both correct);
`input_frames` = window size, not the finite-isolation count.

## Gates (ALL pass before commit)

- **H1:** charge-distance columns re-valued (no self r=0; `gap_to_2nd_charge_r` no longer
  pinned ~1 Å — report its new range); `bin_nearest_charge_r` + `bin_gap_to_2nd_charge_r`
  now populate **> 1 bin** (report the bin distribution — must be non-degenerate); the charge
  `gap`/`nearest` neighbour set now matches the dominance neighbour set (consistency check).
- **M1:** the deliberate-DFT-read-in-selection test makes the guard **THROW** (demonstrated);
  the normal run passes with the real guard clean.
- **M2:** `per_atom_substrate_column_specs.json` units corrected (the three named columns no
  longer `A`; only distance columns `A`).
- **Oracle parity:** every data column byte-identical to the `build1` run EXCEPT the
  intentionally-changed charge-distance columns (`nearest_charge_r`, `gap_to_2nd_charge_r`,
  dependent charge bins) — list exactly which changed; `column_specs.json` units intentionally
  corrected. DFT targets / backbone-audit / all mechanism columns byte-identical.
- R/uniqueness unchanged (846/660/558,360).
- **DISK GUARD:** `df` before emit, ABORT+report if free < 20 G; drop-old (new emit REPLACES
  `…-build1`); deletes by explicit FULL PATH only, never regex; total < 15 G; never write/
  delete under `/shared`.

## Build / commit / output

`cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests`
(unsandboxed). Branch `h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR/checkout. One
atomic commit. Write `src/rediscover/POSTMORTEM_BUILD1FIX.md` (≤45 lines): each fix done;
the guard-fires test result; the new charge-bin distribution (non-degeneracy proof); the
gap/dominance neighbour-set consistency; oracle parity (byte-identical set + the intentional
charge re-baseline + units); disk-guard result; commit hash; run dir. Print ONLY a ≤10-line
summary + that path; DO NOT echo diffs to stdout.
