# Session 2026-05-10 (Part 2) — Tripeptide BB + Neighbor calculators land

This is the continuation writeup from the 2026-05-10 session that
opened with calculator-queue + data-backend setup (part 1 writeup:
`spec/plan/session_2026-05-10_calculator_queue_and_data_backend.md`)
and continued through the night of 2026-05-10 / morning of 2026-05-11
to land the BB + Neighbor tripeptide calculators end-to-end.

## What landed

Three calculators + their substrate layer + smoke tests + the
infrastructure plumbing.

### Source files added

- `src/TripeptideDftTable.{h,cpp}` — libpq-backed loader for the
  ProCS15 tripeptide DFT data from the local `tensorcs15` Postgres
  replica. Forward-declares `PGconn` to keep `<libpq-fe.h>` out of
  callers' transitive includes. Parses geometry + tensor JSONB
  columns. Identifies central residue backbone atoms (N/CA/C/O) by
  the gotham element-pattern heuristic. Exposes
  `QueryNearest(letter, phi, psi, chi1..chi4, n_chi_axes=-1)`. The
  `n_chi_axes` override is the chi-fallback knob: callers walk from
  the residue's natural chi depth down to phi/psi-only on misses,
  dropping chi columns from the SQL `WHERE` clause. (The gotham
  reference impl zero-fills chi3/chi4 and misses K/R/E/M/Q on
  chi-specific lookups; our extension uses actual SQL column drop —
  fixed the ARG/LYS miss-then-fallback-to-phi-only chain that the
  first build exhibited.)

- `src/TripeptidePoseAssembler.{h,cpp}` — the validation layer
  between DFT-record atoms and protein-residue atoms. Two independent
  paths (DFT canonical ordering + typed substrate role cross-check)
  gate atom emission. Records Vec3 residual per atom as an ML
  feature, NOT a stratification stat. See
  `feedback_two_path_validation` and `feedback_residual_as_ml_feature`
  for the disciplines this layer enforces.

- `src/TripeptideBackboneShieldingResult.{h,cpp}` — σ_BB^i
  `ConformationResult`. Per-residue: query DB with chi-fallback,
  assemble pose via `TripeptidePoseSide::Central`, emit per-atom
  rotated tensors + residual vectors. Stores on `ConformationAtom`:
  `tripeptide_bb_shielding_tensor` + `_spherical` + `_match_distance`
  + `_residual_vec` + `_has_match` + `_method_tag`.

- `src/TripeptideNeighborShieldingResult.{h,cpp}` — Δσ_BB^{i±1} per
  Larsen 2015 Eq 3. **Paper-verified cap-side reading.** Per residue
  i: query (i-1)'s tripeptide, assemble C-terminal ALA cap; query
  (i+1)'s tripeptide, assemble N-terminal ALA cap; subtract AAA-at-
  standard-angles reference; accumulate at residue i's atoms.
  Per-direction residual vectors stored separately
  (`tripeptide_neighbor_residual_vec_prev`, `_next`) so the ML model
  can attend to each side's alignment quality independently.

- `tests/test_tripeptide_backbone_shielding.cpp` — BB smoke on
  1UBQ_pm6dh3plus.pdb. Substrate cross-check + residual bounds.

- `tests/test_tripeptide_neighbor_shielding.cpp` — Neighbor smoke on
  1UBQ_pm6dh3plus.pdb. Both frame_types observed; per-direction
  residual sanity bounds.

### Source files modified

- `src/ConformationAtom.h` — added per-atom tripeptide BB +
  neighbor storage fields.
- `src/RuntimeEnvironment.{h,cpp}` — section-aware TOML parser
  (tracks `[databases]` section). New `TensorCs15Dsn()` accessor.
- `src/Session.{h,cpp}` — owns `unique_ptr<TripeptideDftTable>`;
  `LoadTripeptideDftTable()` opens libpq once and holds for Session
  lifetime. Empty DSN = table left null, calculators return nullptr
  at Compute.
- `src/errors.h` — `kSessionTripeptideDbLoadFailed = 0x0005`.
- `CMakeLists.txt` — `find_package(PostgreSQL REQUIRED)` +
  `PostgreSQL::PostgreSQL` link to `nmr_shielding`. New sources +
  tests added.

## Smoke results on 1UBQ_pm6dh3plus.pdb

This is Larsen's PM6-D3H+ optimised geometry — the published-RMSD
validation target from Larsen 2015. Pulled into
`/mnt/expansion/larsen_archive/structures/` during part 1 of this
session.

**BB calculator:**

- 74/76 residues matched (terminal pair skipped — no full φ/ψ).
- 1033/1232 atoms assigned.
- Mean Kabsch BB RMSD 0.55 Å, max 1.48 Å (one outlier).
- No backbone atoms with residual > 3 Å (after the SER OG↔O bug fix).
- ASA SER residues correctly hit `frame_type='orca_input_orientation'`
  (project SER PBE regen); the rest hit `'gaussian_standard_orientation'`
  (Larsen OPBE).

**Neighbor calculator:**

- 75/76 residues received ≥1 neighbor contribution.
- 516 per-atom Δσ accumulations applied.
- Both OPBE and ORCA-PBE frame_types observed in `prev_frame_type` /
  `next_frame_type` per-residue match records (SER as i±1 neighbor
  hits the PBE rows).
- Per-direction residual vectors populated separately.

## Bugs caught + fixed during this session

1. **Wrong Larsen Eq 3 reading site (memory entry was wrong).**
   Initial design (memory entry `project_larsen_neighbor_axa_reuse`)
   said "read at the central residue's atoms in the (i-1)
   tripeptide". This was wrong. Reading Larsen 2015 pp. 2–4 directly
   (the Cβ/Val worked example pins it) shows the actual reading is
   at the **flanking ALA cap atoms** — C-terminal ALA for the i-1
   contribution, N-terminal ALA for i+1. The cap represents the
   residue position in the protein corresponding to i.

   The memory entry has been corrected. The calculator was
   implemented against the paper, not the wrong memory.

2. **Sidechain re-rotation rotating the backbone C/O.** The first
   port of gotham's sidechain re-rotation used the heuristic "rotate
   atoms from CB onward, skip the last 2 which are usually backbone
   C/O". For SER (and any residue where Gaussian's ordering
   interleaves backbone with sidechain — SER's order is `N H CA HA
   CB **C** HB HB OG **O** HG`), the "last 2" are O and HG, not C
   and O. The backbone C got silently rotated as sidechain, producing
   3–5 Å residuals on BB atoms.

   Fix: identify central C and O explicitly via
   `IdentifyCentralBackbone` and skip those specific indices in the
   sidechain rotation loop (regardless of position).

3. **SER OG ↔ backbone O swap.** Greedy element + nearest-distance
   matching saw two O atoms in SER (sidechain OG at position 8,
   backbone O at position 9) and silently swapped them when the
   aligned OG ended up closer to the protein's BB O. The substrate
   cross-check only fired for atoms matching the typed Residue cache
   slots (res.N, res.CA, etc.) — for sidechain atoms it was a no-op.

   Fix: **eager typed binding** of backbone N/CA/C/O to the typed
   `Residue` cache slots BEFORE any element+nearest-distance search.
   Backbone role-pinning at both ends. The remaining greedy matching
   runs only on sidechain atoms (where the swap risk is residue-
   internal sidechain heteroatoms — still not perfectly clean, but
   the backbone is now airtight).

4. **Chi-fallback losing K/R/E/M/Q hits.** Initial port passed
   `chi3=0, chi4=0` to `QueryNearest` on fallback, but those rows
   exist with non-zero grid values for chi3/chi4. The query
   `WHERE chi3=0 AND chi4=0` matched no rows. Fix: `n_chi_axes`
   override parameter that drops chi columns from the WHERE clause
   on fallback. After fix: 74/76 matched (was 51/76).

5. **The wrong "central atoms" reading for Δσ_BB^{i±1}.** Already
   covered as #1 above but worth re-emphasising: the memory entry
   was authoritative-looking until verified against the paper.

## Architectural decisions made

### Two-path validation discipline

Cross-substrate matching (DFT records ↔ protein residues) must be
typed at both ends. The matching layer:

1. Identifies cap atoms by canonical Gaussian/ORCA ordering on the
   DFT side (build-time or element-pattern walk).
2. Maps each cap slot to a typed `Residue` cache slot on the protein
   side.
3. Cross-validates: protein atom's typed `BackboneRole`/`Locant`
   must match the expected role for the cap slot. Mismatch = fail
   loud.

Per the new memory entry `feedback_two_path_validation`. The current
implementation has the protein side fully typed and the DFT side
typed only on BB N/CA/C/O (gotham heuristic). The proper fix — fully
typed DFT side — is the queued typed-topology refactor.

### Residual as Vec3 ML feature

Per the new memory entry `feedback_residual_as_ml_feature`. Direction
+ magnitude both matter to the upstream model. Emit per-atom as
(N, 3) NPY alongside the tensor. Don't gate emission on residual
magnitude.

### Frame_type passthrough

The DB row's `frame_type` column (`'gaussian_standard_orientation'`
for OPBE, `'orca_input_orientation'` for SER PBE) is stashed per
record + propagated to per-atom `tripeptide_bb_method_tag` (byte:
0=none, 1=OPBE, 2=PBE) for downstream calibration to route methods
separately. Per the existing `project_serine_pbe_discontinuity`
discipline.

## Queued work

### Typed tripeptide topology refactor

The sidechain matching beyond Cβ is still heuristic (element +
nearest-distance with cross-check when matched-atom happens to be a
typed slot). To make matching role↔role at both ends, the DFT side
needs its own typed substrate. The Plan agent designed it; user
added pose-metadata enrichment (governing dihedral + rotation axis
+ bond parent per atom for residual decomposition).

**Full spec:** `spec/plan/typed-tripeptide-topology-design-2026-05-10.md`.

**Memory entry:** `project_typed_tripeptide_topology_design`.

**Cost:** 5–6 hours, two sessions.

### Run on 1P9J trajectory

The smoke tests run on 1UBQ_pm6dh3plus.pdb (static, Larsen's
geometry). The actual Stage 2 study system is 1P9J — 750 MD frames
at 20 ps stride. Run both calculators against the trajectory frames
to surface bugs under MD geometry variance (substantially more
chi-grid mismatch than the static crystal/PM6 geometry).

### HBondHα implementation

The third queued slice from the calculator queue (per
`project_post_csa_calculator_queue`). Kernel × η_Hα form per
`project_hbond_halpha_design`. Larsen's HB grid data is in the
`hydrogenbondnmrlogs.tar` ERDA archive we pulled in part 1 but
haven't ingested.

### Loth 2005 PDF acquisition

For A.2 AmideTensorGeometryResult. Still held pending the paper
acquisition.

## Memory entries written/updated this session

**Corrected:**
- `project_larsen_neighbor_axa_reuse` — fixed the wrong "central
  atoms" reading; pinned to the paper's Cβ/Val example.
- `project_post_csa_calculator_queue` — marked BB + Neighbor LANDED;
  HBondHα remains queued.
- `reference_gotham_assembler` — noted the port + four improvements
  (chi3/chi4, frame_type, two-path validation, Vec3 residual).

**New:**
- `feedback_two_path_validation` — discipline.
- `feedback_residual_as_ml_feature` — discipline.
- `project_tripeptide_calculators_landed` — what shipped this session.
- `project_typed_tripeptide_topology_design` — queued architectural
  work.

`MEMORY.md` index updated with the new entries.

## Build state

- `nmr_shielding` library: builds clean.
- `nmr_extract` binary: builds + links clean.
- `unit_tests`: 114/114 pass — no regression.
- `structure_tests`: TripeptideBackboneShieldingTest + 
  TripeptideNeighborShieldingTest both pass.
- All compiler diagnostics during the session were clangd-only
  (the `compile_commands.json` doesn't see Eigen include path);
  the actual GCC build is clean.

## Things flagged for next session

1. **1P9J trajectory smoke** (most important; the actual study
   system).
2. **Typed tripeptide topology refactor** (Phase 1: psql dump →
   TOML → generator extension).
3. **HBondHα implementation.**
4. **A.2 Amide / Loth 2005 PDF acquisition.**

## Provenance

- Session opened 2026-05-10 with calculator-queue/data-backend setup
  (part 1).
- Tripeptide port work began with the post-compact prompt's pointer
  to the gotham assembler.
- User intervention on residual-as-ML-feature + two-path-validation
  + tripeptide-model-builder-as-future-work all shaped this session's
  architectural decisions.
- Plan agent invoked for the typed-topology design; full output
  captured in `spec/plan/typed-tripeptide-topology-design-2026-05-10.md`.
- Doc revision pass at end of session (this writeup + memory entry
  updates + CLAUDE.md update + post-compact prompt) — explicitly
  requested by user to keep the durable record honest with what
  landed.
