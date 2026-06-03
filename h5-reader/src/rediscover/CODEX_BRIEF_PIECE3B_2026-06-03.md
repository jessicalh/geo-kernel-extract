# Codex brief — Loop 3b: the full channel-completion + T1/dia/para targets + comparison-method paths

Status: **landed** in chunks `9965c2bd1f38f4c138e7dc3836206b45f3e91525`, `1d0f56aac4c8a53310b1214663bd65aa8f8fc0b3`, and `5a39288b30b069bac31e25c38d1c6df5c981d91f`.

You own the grind; the owner vets. This is the FAT, deliberate emit: wire every
cheap-available channel + the full shielding-target decomposition into
`per_atom_substrate`, additive on commit `b583d7c` (Loop 3). The point is to never have
to walk this path again — codex compute is cheap, a re-walk is not. Read
`ALLATOM_FIT_SPEC_2026-06-03.md` + `PerAtomSubstrate.cpp` (the existing mechanism blocks
are your template) first.

## THE GATE — derivation in C++ memory; only reduced (atom,frame) aggregates cross the edge

This is the discipline that keeps us off a terabyte. The node-store contract, as one rule:

- **The emit is keyed (atom, frame) ONLY.** Everything written is a per-(atom,frame)
  REDUCED AGGREGATE. NO per-source rows, NO pairwise dump — those are the TB path (the
  68 GB we already lived). The per-source intermediate is computed, folded, and DISCARDED
  per atom-frame in C++ memory; it never reaches disk.
- **C++ derives + reduces in memory; Python does only statistics.** Every reduction —
  sum-over-neighbours, per-ring-type / per-bond-category aggregation, the cross-method
  columns — happens in the C++ spine BEFORE the edge. Python receives reduced aggregates
  and only fits / correlates / partitions. Python never derives, never opens the H5,
  never sees a per-source row. (If Python would need per-source to aggregate, you've been
  forced to dump per-source → TB. Keep the aggregation in C++.)
- **SIZE-GATE (pre-commit check):** before any chunk commits, ESTIMATE the flat footprint
  (rows × new float64 columns × 8 B; anchor: v1 ≈ 1.05 GB, ~4.5 MB/column at 558,360 rows).
  Expected total after 3b ≈ **2–3 GB** (drop-old). If any estimate trends toward tens of
  GB, STOP and report — that means a per-source/pairwise axis crept in. GB is fine; a
  second key axis is the alarm.
- Per-type / per-category breakdowns ARE allowed flat (they are (atom,frame) aggregates,
  ~1 GB) — they are reduced in C++ memory, emitted as columns, NOT via per-source rows.

## HARD RULES (read twice)

- **ZERO Python on any emitted value.** Everything computed in the C++ spine, mirroring
  the existing `PerAtomSubstrate.cpp` blocks (typed, no string dispatch, ArraySpec per
  output, present-flags, units/irreps/mechanism metadata). A Python emit path is a
  rollback offense. Producer/extractor/CMake/ctest untouched. Read existing H5/NPY only;
  NO re-extraction.
- **Fail-loud-AND-LOCATE.** For EACH item below, first LOCATE where it lives in the
  producer output — read `python/nmr_extract/_catalog.py` for the array name/group/shape;
  the datum may be packed in a multi-component array (`hbond_scalars`), per-residue
  (DSSP), or under a non-obvious name. Only after locating do you wire it. If — after
  locating — it is genuinely absent from THIS run's H5/NPY, emit present=0/NaN and record
  the item + the reason in the manifest. "Absent under the canonical name" is NOT
  "absent" — locate first. (Locate ≠ glob: read the documented catalog, don't try-and-fail.)
- **Best + backup.** Prefer the producer's own physical assignment over a hand-rolled
  mapping. Where an assignment is genuinely uncertain, emit BOTH the best form and a
  backup form so the fit decides (the existing `sum_all`/`sum_valid` pattern).

## A. Shielding target decomposition — add T1, dia, para

Currently the substrate emits target T0 + T2. Add, from the raw `orca_total` (and
`orca_diamagnetic` / `orca_paramagnetic`) 3×3 — all recoverable, NO new ORCA run:

- **`target_T1`** (rank-1 / antisymmetric, 3 components, irreps `1o`). NMR-silent for
  shifts → label it a kernel-completeness DIAGNOSTIC, not a shift target. **Verify its
  frame** the way T2 was verified (Kabsch ORCA-frame vs H5 positions; the existing
  `buckingham_efield_target_T1_unverified` is unverified — confirm before trusting). If
  the frame check fails, emit it flagged `t1_frame_unverified`, do not silently trust.
- **`target_para_T0/T1/T2`** and **`target_dia_T0/T1/T2`** — the diamagnetic/paramagnetic
  split of the target. The paramagnetic part is where the chemistry lives; emitting it
  lets the fit target para separately from the near-additive dia.
- Keep existing `target_T2` / `target_T0` byte-identical (append the new target sidecars).

## B. Comparison-method paths — the FREE validation (emit ALL paths per quantity)

The owner's priority: emit every independent path to each quantity so cross-method
agreement is evaluable with NO DFT, on all frames. Mirror the existing T2 block emit.

- **Ring current — 4 paths:** `bs_shielding` (+ `bs_per_type_T0/T2`, `bs_total_B`,
  `bs_ring_counts`), `hm_shielding` (+ `hm_per_type_T0/T2`), `ringchi_shielding`, and the
  existing Johnson-Bovey fold (already emitted). Emit each as its own labeled mechanism so
  bs/hm/jb can be correlated as the comparable shielding set; ringchi is opposite-convention and must not be naive-slope-scored.
- **π-quadrupole / dispersion per-type:** `pq_per_type_T0/T2`, `disp_per_type_T0/T2`.
- **McConnell paths + per-category:** `mc_category_T2`, `mc_scalars`,
  `mopac_mc_category_T2`, `mopac_mc_scalars`, `mopac_bond_orders` (mechanistic bond-order→Δχ).
- **Field / EFG paths (definitional check + divergence STUDY — NOT a clean agreement edge):**
  `mopac_coulomb_E`, `mopac_coulomb_efg_backbone/aromatic`, `mopac_coulomb_scalars`,
  `aimnet2_efg` (+ `aimnet2_efg_aromatic/backbone`), `water_efg`, `water_efield_first`.
  (FF14SB field, APBS E/EFG, mopac_coulomb_shielding already emitted.) **Analysis framing
  (record for the next analysis, do not treat field cross-method as slope≈1):** unlike ring current,
  these will NOT agree — FF14SB(vacuum) / APBS(solvated) / MOPAC / AIMNet2(gas-phase) /
  explicit-water diverge by construction (screening, polarization). The CLEAN field
  validation is WITHIN-method definitional — `EFG = ∇(field)`, `field = Coulomb-sum(charges)`
  for the SAME method (a break = pipeline bug). Cross-method comparison is a
  structured-divergence study (the screening/polarization IS the signal), NOT a positive control.
- **Charge paths:** `eeq_cn` (coordination — the missing one), `mopac_scalars`.

## C. H-bond — best + backup (the donor/acceptor question, done right)

- **Best (producer's physical assignment):** the four Larsen per-class shielding T2 arrays
  `larsen_hbond_1pHB / 2pHB / 1pHaB / 2pHaB` (+ `larsen_hbond_water_term`) — each already
  assigned to the correct atom per Larsen 2015 Table 2 (donor/acceptor + amide/Hα resolved).
- **Geometry:** the full `hbond_scalars` (nearest_dist, 1/r³, count, McConnell Σ) — count
  already emitted; add nearest_dist + the rest.
- **Backup conditioner (uncertain mapping → emit both forms):** DSSP donor/acceptor —
  (i) the chemical per-atom flag (donor→HN/N, acceptor→O/C′) AND (ii) the raw per-residue
  donor/acceptor energy+count (`dssp_hbond_energy`, `hbond_donor/acceptor_partner/energy`)
  carried onto the peptide-plane atoms. Let the fit decide which encoding carries signal.

## D. Conditioners (input-side partition axes)

`dssp_ss8` (secondary structure), `dssp_chi`, `omega_actual`, `pyramidalization`,
`ring_geometry`. Slicers, not predictive features (consistent with existing conditioning policy).

## Excluded (do NOT emit — record as intentional)

`coulomb_*` (retired vacuum Coulomb), `delta_* / wt_ / mut_` (Stage-1 mutant mode),
`mopac_global`, `tripeptide_bb_shielding` (ProCS15 — semi-circular as a feature;
reference-only, not a feature).

## Resilience (this is a big run — bank progress against the codex context-image bug)

Implement + gate + commit in THREE atomic family-chunks so a crash never loses banked
work: (1) targets A + ring-current paths B; (2) the rest of B (field/EFG/charge/McConnell
paths); (3) C (hbond best+backup) + D (conditioners). Each chunk: oracle-parity +
backbone-regression + per-channel coverage gates BEFORE its commit. Branch
`h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR/checkout.

## Gates (every chunk, ALL pass before that chunk commits)

- R & uniqueness: 846 / 660 / 558,360; dense `row_id`; unique `(atom_index,
  original_frame_index)`; sidecar first-dims match.
- New-channel coverage: present flags 0/1; present⇒finite; axis correct (static repeats
  over frames, per-frame channels vary); report each new channel's present-count + range +
  (if absent) the located-but-absent reason.
- Oracle parity UNCHANGED: every pre-existing NPY sidecar byte-identical; existing CSV
  columns value-identical row-for-row; `ring_identity` + default `query_results/`
  unchanged; DFT-target T0/T2 parity exact.
- Backbone regression UNCHANGED: `backbone_audit.npy` byte-identical.
- T1 frame: report the Kabsch frame-check result; flag if unverified.

## OUTPUT — TINY VERDICT

Write `src/rediscover/POSTMORTEM_PIECE3B.md`, **≤60 lines** (it's a big emit): per-chunk
commit hashes; each gate pass/fail; per-channel present-counts + ranges + any
located-but-absent items with reasons; the T1 frame-check verdict; run-dir path. Print
ONLY a ≤12-line summary + that path. DO NOT echo diffs to stdout.
