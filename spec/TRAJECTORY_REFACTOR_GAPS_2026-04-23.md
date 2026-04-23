# Trajectory Refactor — Open Gaps (2026-04-23)

Landing session: replaced `GromacsProtein`/`GromacsProteinAtom`/
`GromacsRunContext`/`GromacsFinalResult` with the
`TrajectoryProtein`/`TrajectoryAtom`/`TrajectoryResult`/`Trajectory`/
`RunConfiguration`/`RunContext`/`DenseBuffer`/
`SelectionEmittingTrajectoryResult` + `BsWelfordTrajectoryResult`
pattern. Classes fully defined, library + driver + tests compile, two
of three streaming tests pass, third (end-to-end PerFrameExtractionSet
with AIMNet2 + APBS) running overnight.

This document inventories every known gap left by the landing, ordered
by consequence. Low-hanging fruit is marked **[LHF]**.

## Blocking gaps (must close before fleet can run)

### G1. Stride support in `Trajectory::Run` **[LHF]**

Production `RunAnalysis` used `STRIDE=2` over ~1200-frame trajectories
to harvest ~600 frames per protein. `Trajectory::Run`'s Phase 4 loop
currently iterates every frame (stride=1). The 685-protein fleet
extraction targets the stride=2 cadence.

**Shape of fix**: `RunContext` gains `size_t stride_ = 1` field with
setter `SetStride(size_t s)`. `Trajectory::Run` Phase 4 reads
`ctx.Stride()`; if > 1, calls `handler_->Skip()` (stride−1) times
between each `handler_->Next()` dispatch. `frame_indices_` records
the original XTC index so downstream consumers see the gap.

**Effort**: ~20 lines. Next session's first item.

### G2. MOPAC skip flag not honoured in per-frame run (!) **[LHF if reproduces]**

Log from 2026-04-23 test run shows `MopacResult::Compute [BEGIN]` at
frame 0 despite `PerFrameExtractionSet::per_frame_opts_.skip_mopac =
true`. Either:
  a. `OperationRunner::Run` doesn't check `skip_mopac` for frame 0
  b. The flag is being overwritten between `RunContext::PerFrameRunOptions()`
     and `OperationRunner::Run` somewhere along the handler path
  c. Missing check in my `GromacsFrameHandler::Open` init_opts path

**Shape of fix**: investigate `OperationRunner::Run` path for
`skip_mopac` honour. If the flag is correctly checked in steady state
but not frame 0, patch. If flag is lost in the options plumbing, trace
the assignment.

**Effort**: 30 min read + fix. Next session second item.

### G3. Per-frame stride / selected-frame mechanism for FullFatFrameExtraction

`FullFatFrameExtraction` is meant for MOPAC on selected frames only
(DFT pose set, μs harvester checkpoints). Neither `Trajectory::Run`
nor `RunContext` currently supports a selected-frame filter.

**Shape of fix** (two options per WIP §6):
  a. Pre-filter XTC externally; run FullFat on the filtered file.
  b. `RunContext::SetSelectedFrameIndices(std::set<size_t>)`;
     `GromacsFrameHandler::Next()` skips non-selected frames.

Option (b) is cleaner but requires scan-pass to produce the index set.
Deferred until SelectionEmittingTrajectoryResult subclasses exist
(G4).

### G4. SelectionEmittingTrajectoryResult has zero implementers

The mixin interface + `Trajectory::Run` collection loop exist, but no
concrete TrajectoryResult implements it. That means `ScanForDftPointSet`
cannot emit DFT pose selections yet — the whole scan+extract two-pass
use case is dormant.

**Shape of fix**: implement `ChiRotamerSelectionTrajectoryResult`
(per WIP Appendix F ScanForDftPointSet row), `RmsdTrackingTrajectoryResult`,
`RmsdSpikeSelectionTrajectoryResult`, `DftPoseCoordinatorTrajectoryResult`.
Each is one subclass with `Compute` detecting events + `SelectionRecords()`
accessor returning the accumulated records.

**Effort**: 1 session. Not blocking fleet if scan+extract is skipped.

### G5. Byte-parity vs pre-refactor trajectory output not verified

Only `BsWelfordTrajectoryResult` is migrated. The other ~40 Welford
fields that the old `GromacsProteinAtom` tracked (per-atom mc_T0,
hm_T0, pq_T0, disp_T0, aimnet2_charge, water_emag, water_n_first,
eeq_charge, sasa, chi1–4_cos, ss8_transitions, phi_cos, psi_cos,
half_shell, dipole_cos, nearest_ion_dist, wp_dipole_x/y/z,
wp_normal_x/y/z, wp_asymmetry, wp_alignment, wp_coherence,
wp_first_shell_n, dssp_hbond_energy, mean_bond_angle_cos, bond
length Welford, etc.) no longer accumulate.

**Shape of fix**: follow-up sessions add one TrajectoryResult per
source ConformationResult per WIP Appendix F:

  - `McWelfordTrajectoryResult` (copy of BsWelford for mc_T0, mc_T2mag)
  - `HmWelfordTrajectoryResult`
  - `RingSusceptibilityWelfordTrajectoryResult`
  - `PiQuadrupoleWelfordTrajectoryResult`
  - `DispersionWelfordTrajectoryResult`
  - `HBondWelfordTrajectoryResult`
  - `ApbsFieldWelfordTrajectoryResult`
  - `AIMNet2WelfordTrajectoryResult` (charge + EFG)
  - `SasaWelfordTrajectoryResult`
  - `WaterEnvironmentWelfordTrajectoryResult`
  - `HydrationWelfordTrajectoryResult`
  - `EeqWelfordTrajectoryResult`
  - `DsspDynamicsWelfordTrajectoryResult` (phi/psi/chi cos +
    transition counters)
  - `BondLengthStatsTrajectoryResult` (per-bond length Welford, Result-
    internal per §3 option (b))
  - Position time series (`PositionsTimeSeriesTrajectoryResult` with
    DenseBuffer<Vec3>) for movie + RMSF

**Effort**: ~5–8 sessions. This is the remainder of WIP 1b.

## Non-blocking architectural gaps

### G6. AnalysisWriter hasn't dissolved yet

`RunAnalysis` in `nmr_extract.cpp` still uses `AnalysisWriter` as its
output surface. Per WIP §7 + Appendix H sub-checkpoint 1e,
AnalysisWriter should migrate into a family of `*TimeSeriesTrajectoryResult`
classes each owning a `DenseBuffer` + emitting its own H5 group via
`WriteH5Group`. Dissolution waits for those classes to land.

**Effort**: ~1 session after G5 is substantially done.

### G7. Positions time series not emitted in trajectory H5

`GromacsProtein` used to store raw `frame_positions_` (T×N×3 doubles)
and emit them at `/positions` in the trajectory H5. `TrajectoryProtein`
doesn't. The replacement is `PositionsTimeSeriesTrajectoryResult`
owning a `DenseBuffer<Vec3>` (part of G5). Until then, movies +
RMSF-from-positions depend on the old H5 shape.

**Effort**: ~2 hours of one G5 session.

### G8. Per-atom identity not emitted beyond element + residue_index

`TrajectoryProtein::WriteH5` emits `/atoms/element`,
`/atoms/residue_index`, `/atoms/pdb_atom_name`. The richer per-atom
typed identity (NmrAtomIdentity: `IupacAtomPosition`, `NmrClass`,
`MethylGroup`, etc.) is explicitly deferred to post-rollup per WIP
Appendix H deferred-stage.

**No action**: deferred per WIP, not a mistake.

### G9. Two-pass scan+extract pattern is dormant in the driver

`nmr_extract.cpp RunTrajectory` is now single-pass
(PerFrameExtractionSet) instead of the old pass1 scan + pass2 extract
on selected frames. Returns when G4 lands.

**No action**: documented in code comment.

## Quality-of-life gaps

### G10. `smoke_tests` BINARY DIFF regression (pre-session 2026-04-23) **[LHF triage]**

Pre-dates trajectory refactor. Both `SmokeTest.NoDft` and
`SmokeTest.WithDft` fail binary comparison against blessed baselines
on ~30 NPY arrays (`bs_per_type_T0.npy`, `mc_shielding.npy`,
`coulomb_*`, `mopac_*`, `pq_*`, `ring_*`, `eeq_*`, `hbond_*`,
`disp_*`, `dssp_chi.npy`). Sizes match, content differs.

**Investigation path**:
  a. Check blessed-baseline timestamps (`tests/golden/blessed/smoke/
     {nodft,withdft}/`) vs `git log` on each affected calculator.
  b. If blessed pre-dates a recent calculator change, re-bless after
     confirming output is correct.
  c. If blessed post-dates calculator changes, the differences are
     non-determinism or a real regression — use `h5diff`/`numpy
     .allclose` to quantify.

**Effort**: 1–2 hours diagnostic. Tracked in memory entry
`project_smoke_test_regression_20260423.md`.

### G11. Clangd/editor indexing doesn't find Eigen **[LHF dev UX]**

clangd reports cascading errors on every Trajectory-family header
because it can't find `Eigen/Dense`. CMake finds it fine. The project
CLAUDE.md mentions pointing clangd via `--compile-commands-dir=build/
<preset>`. Per-preset symlinks are forbidden, but a
`.clangd` file in `src/` could add the Eigen include path.

**Effort**: 10 minutes. Not blocking anything.

### G12. `GromacsFrameHandler::AddFramePath` + `frame_paths_` tracking dropped

The old `GromacsProtein::frame_paths_` accumulated per-frame NPY
output directory paths. My refactor didn't port this to any of the new
classes — `handler.Next(opts, output_dir, ...)` still creates the
frame_NNNN/ directories and writes NPYs, but nobody records the paths.
Low consequence — downstream consumers discover frames by iterating
`output_dir/frame_*` anyway.

**Effort**: zero (just document). Or 30 min if a TrajectoryResult
wants the paths.

### G13. Library is not `#include "GromacsRunContext.h"` clean in comments

Several source files have stale comments referencing `GromacsRunContext`
/ `GromacsProtein` / `AccumulateFrame`:
  - `src/JobSpec.h:32`
  - `src/OperationRunner.h:62, 74`
  - `src/GromacsEnergyResult.h:90`, `src/GromacsEnergyResult.cpp:9`
  - `src/BondedEnergyResult.h:57`

Mechanical sweep to update comments to the new class names.

**Effort**: 20 min. Low value but tidies the codebase.

## Inventory summary

| ID | Gap | LHF? | Blocking? | Session estimate |
|---|---|---|---|---|
| G1 | Stride support | yes | yes | 0.1 session |
| G2 | skip_mopac honoured | yes | yes | 0.1 session |
| G3 | Selected-frame FullFat | no | partial | 0.3 session |
| G4 | Selection mixin impls | no | no (but scan+extract dead) | 1 session |
| G5 | Other ~40 Welfords | no | yes (byte-parity) | 5–8 sessions |
| G6 | AnalysisWriter dissolve | no | no | 1 session after G5 |
| G7 | Positions time series | no | no | rolled into G5 |
| G8 | NmrAtomIdentity emission | no | no (WIP deferral) | separate stage |
| G9 | Two-pass pattern | no | no | rolled into G4 |
| G10 | smoke_tests bless | yes | no | 0.1 session diagnostic |
| G11 | clangd Eigen path | yes | no | 0.1 session |
| G12 | frame_paths tracking | no | no | 0.1 session if needed |
| G13 | Stale comments sweep | yes | no | 0.1 session |

## Rough sequencing

**Immediate (next session, ~1 sitting):** G1 stride + G2 MOPAC flag +
G11 clangd + G13 comments sweep. All LHF, close them together. Leaves
the trajectory architecture in a clean code-complete state before the
long grind on G5.

**Short-term (2–3 sessions after):** G10 smoke diagnostic (re-bless or
quantify non-determinism), G4 selection-emitting subclasses for scan
mode — unlocks G3 / G9 two-pass pattern.

**Long-term (5–8 sessions):** G5 Welford migrations, one source
ConformationResult at a time. G7 rolls in naturally. G6 dissolves
AnalysisWriter after.
