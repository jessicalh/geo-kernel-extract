# Codex brief — Loop 3 (Piece 3): channel-completion emit (missing classical mechanisms)

> **Historical run brief — not current truth (trued 2026-06-04).** Session
> provenance only; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Status: **landed** at `b583d7c87af2291ef7c025c828be60f2de0ecff9`; retained as the execution brief.

You own the grind; the owner vets. Goal: emit into `per_atom_substrate` the classical
mechanism channels the H5 ALREADY carries but the rediscover `Catalog` does not read,
so the joint fit and the between-calculator network see the full mechanism set. The
standout is **H-bond** — it targets the under-explained HN / O strata directly. All
additive; the producer is untouched; you are wiring through EXISTING H5 data.

## HARD RULE — emitter discipline, ZERO Python (read this twice)

**Every emitted value is computed in the C++ SPINE. Under NO circumstances does Python
compute, derive, transform, or touch any emitted value — not now, not "just for this
one channel," not as a "quick helper."** Python's only role anywhere downstream is to
READ the emitted substrate later. This is model-is-the-spine + no-piracy +
keep-Python-complexity-knocked-back; a Python physics/derivation path in the emit is a
rollback offense. If you feel reached for Python to produce a column, STOP and report
instead.

Follow the emitter rules, mirroring the existing exemplars in `PerAtomSubstrate.cpp`:

- Extend `Catalog` with one ArrayId per new channel, wired to the H5 buffers
  `QtTrajectoryH5` ALREADY loads for the reader — `QtFrame::hbondShielding` /
  `pqShielding` / `dispShielding` / `waterEfield` / `sasa` / `eeqCharge` /
  `eeqCoordinationNumber` / `hmShielding` / `rsShielding` read these today. **REUSE
  those buffers.** Add a `QtTrajectoryH5` loader ONLY if a buffer genuinely isn't loaded
  (mirror the Loop-1 2-line `n_per_atom` fallback; minimal, additive).
- Emit through the EXISTING `per_atom_substrate` C++ carrier (carrier on
  `RunTraversal`, NOT a sibling driver). Mirror EXACTLY how the existing mechanism
  blocks (`ring_jb`, `charge_q_over_r3`, `mc_lit`, `ff14sb_field`, `apbs_efg`) are
  aggregated and written: typed code (no string dispatch), ArraySpec per output,
  present-flags, units/irreps/mechanism metadata in `per_atom_substrate_column_specs.json`.
  Do for the new channels EXACTLY what the file already does for ring/charge/mc.
- Read existing H5 slabs ONLY. NO producer change, NO new physics computation, NO
  re-running `nmr_extract`. These calculators already ran; wire their outputs through.
- No file discovery / fail-loud: if a channel's slab is ABSENT in this run, REPORT it
  and emit present=0 / NaN. Do NOT fabricate, glob, or try-and-fail.

## The channels (mirror each on the existing block of the same shape)

Per-frame T2 shielding (mirror `ring_jb_T2` / `mc_lit_T2`):
- **Larsen H-bond** shielding T2 + geometry the H5 carries (nearest dist scalar, nearest
  dir vec3, count, donor/acceptor flags) — mechanism `hbond`.
- **π-quadrupole** shielding T2 — mechanism `pi_quadrupole`.
- **dispersion** shielding T2 — mechanism `dispersion`.
- **Haigh-Mallion** shielding T2 AND **ring-susceptibility (ringchi)** shielding T2 — the
  producer's alternate ring kernels (H5 has bs/hm/rs; we emit only bs/jb today) —
  mechanism `ring_current_alt` (for the ring agreement edge).

Per-frame vector/scalar context (mirror `ff14sb_field` / `apbs_E`):
- **water/hydration**: water E-field vec3 + n-first/second shell counts + half-shell
  asymmetry + dipole cos — mechanism `water_field`.
- **SASA** scalar + SASA normal vec3 — mechanism `sasa`.

Per-atom static scalars (mirror Loop-1 `ff14sb_charge` / `mopac_welford_mean_charge` —
repeat over the atom's frames):
- **EEQ charge + EEQ coordination number** — mechanism `eeq`.

T2 frame: emit in the SAME raw lab frame as the existing T2 blocks (the H5 kernels were
verified in the H5 frame ≈ lab, ~1e-4°). Mirror the existing T2 emit convention exactly;
impose no new frame.

## Gates (same discipline as Loop 1 — ALL pass before commit)

- R & uniqueness: 846 / 660 / 558,360; dense `row_id`; `(atom_index,
  original_frame_index)` unique; sidecar first-dims match.
- New-channel coverage: present flags 0/1; present⇒finite; axis correct (EEQ constant
  over an atom's 660 frames; per-frame channels vary); report each channel's
  present-count + a value range.
- Oracle parity UNCHANGED: every pre-existing NPY sidecar byte-identical to the prior
  committed substrate; existing CSV columns value-identical row-for-row; `ring_identity`
  + default `query_results/` unchanged; DFT-target parity exact.
- Backbone regression UNCHANGED: `backbone_audit.npy` byte-identical; new blocks perturb
  no old reducer beyond fp tolerance.

If ANY gate fails: STOP, do not commit, report which.

## Build / commit

`cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests`
(unsandboxed). Emit to a fresh drop-old run dir. Branch `h5-reader-pysr-spike` —
**NEVER merge/switch/rebase/PR/checkout.** One atomic commit. `_catalog.py`:
`column_specs.json` is source of truth (Q2); touch it only if required, and flag it.

## OUTPUT — TINY VERDICT (new cadence, important)

Write `src/rediscover/POSTMORTEM_PIECE3.md`, **≤40 lines**: drift (what landed vs asked;
any loader added; any channel slab absent in the H5); each gate pass/fail; per-channel
present-counts + ranges; commit hash; run-dir path. Then print ONLY a ≤10-line summary
plus that path. **DO NOT echo diffs to stdout** — diffs live in git (`git show --stat`).
Keep stdout small; the owner's reading budget is the scarce resource.
