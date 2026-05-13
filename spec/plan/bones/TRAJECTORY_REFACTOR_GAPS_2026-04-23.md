# Trajectory Refactor — Open Gaps (historical, 2026-04-23)

**HISTORICAL RECORD.** Session-end gap inventory from the 2026-04-23
trajectory-scope refactor landing. Individual gap statuses have
moved since this file was written. **Do not use this file to plan
current work.** For current truth, use:

- `spec/pending_decisions_20260423.md` — items that need user choice
- `spec/doc_wrongness_20260423.md` — current doc/code contradictions
- The code in `src/` — authoritative source of what is implemented
- `OBJECT_MODEL.md` trajectory-scope entities + `PATTERNS.md §§13-18`
  — current shape and patterns

## Known drift from this file's content (at minimum)

- **G1 (stride support)**: implemented. Lives on `RunConfiguration`
  (`SetStride` / `Stride`), consumed in `Trajectory::Run` Phase 7 via
  `handler_->Skip()` between dispatches. Not on `RunContext` — that
  class was dissolved.
- **G3 (selected-frame mechanism)**: `RunContext` no longer exists.
  Selected-frame filtering is not currently implemented; if and when
  it returns, the surface will live on `RunConfiguration` or directly
  in `Trajectory::Run`, not on a separate `RunContext`.
- **G4 (selection-emitting Results)**: mixin interface +
  `dynamic_cast` collection approach was abandoned. Direct
  `traj.MutableSelections().Push(...)` during a TR's own `Compute` /
  `Finalize` replaces it. `ChiRotamerSelectionTrajectoryResult`
  landed as a concrete emitter but is not currently attached by any
  `RunConfiguration` factory — open decision tracked in
  `spec/pending_decisions_20260423.md`.
- **G5 (other Welford migrations)**: partially landed —
  `BsWelford`, `BsShieldingTimeSeries`, `BsAnomalousAtomMarker`,
  `BsT0Autocorrelation`, `BondLengthStats`, `PositionsTimeSeries`.
  The remaining classes in the original list (`HmWelford`,
  `McConnellWelford`, `ApbsFieldWelford`, `SasaWelford`, etc.) are
  still follow-up work.
- **G13 (stale-comments sweep)**: the specific file:line references
  (`src/JobSpec.h:32`, `src/OperationRunner.h:62, 74`,
  `src/GromacsEnergyResult.h:90`, `src/GromacsEnergyResult.cpp:9`,
  `src/BondedEnergyResult.h:57`) no longer contain the flagged
  `GromacsRunContext`/`GromacsProtein`/`AccumulateFrame` references
  as of 2026-04-23.

Other gap statuses (G2, G6-G12) have not been individually
re-audited; trust the code over any claim here.

## Eight-phase Trajectory::Run — phase numbers in the text below are stale

The original text below describes a 5-phase Run with stride work
"in Phase 4." The current `Trajectory::Run` is eight-phase:

1. Open handler.
2. Read frame 0 + `tp.Seed`.
3. Attach TrajectoryResults.
4. Validate dependencies + session resources.
5. Build per-frame `RunOptions` base.
6. Frame 0 dispatch.
7. Per-frame loop (stride consumed here).
8. `FinalizeAllResults`.

See `OBJECT_MODEL.md` trajectory-scope Run description and
`src/Trajectory.cpp` for the authoritative phase list.

---

[Session-end gap content below this line is preserved for history.
Treat it as a 2026-04-23 snapshot, not a current work queue.]

## Original blocking gaps (at 2026-04-23 landing)

- **G1** Stride support in `Trajectory::Run` — CLOSED (see drift note above).
- **G2** MOPAC `skip_mopac` not honoured at frame 0 — status not
  re-verified post-landing; treat as unverified.
- **G3** Per-frame stride / selected-frame for `FullFatFrameExtraction`
  — still open; mechanism shape has moved (see drift note).
- **G4** `SelectionEmittingTrajectoryResult` implementers — mechanism
  changed (see drift note); `ChiRotamerSelection` landed but
  unattached.
- **G5** Other ~40 Welford migrations — partial (see drift note).

## Original non-blocking architectural gaps

- **G6** AnalysisWriter dissolution — still follow-up.
- **G7** Positions time series in trajectory H5 — landed
  (`PositionsTimeSeriesTrajectoryResult`).
- **G8** Per-atom identity emission — still deferred;
  `TrajectoryProtein::WriteH5` emits element + residue_index +
  pdb_atom_name only. The NmrAtomIdentity-on-Atom design originally
  proposed for this gap was **superseded 2026-04-28** by the
  `LegacyAmberTopology` + calculator-contract architecture
  (`spec/plan/openai-5.5-strong-architecture-layout.md` and memory entry
  `project_proteintopology_architecture` — typed contract attached
  to Protein with typed semantic fields absorbed in). When G8 is
  closed, it will be via `LegacyAmberTopology`'s output projection
  surface at the H5 boundary (AMBER strings native, IUPAC and BMRB
  as derived projections).
- **G9** Two-pass scan+extract pattern — unchanged status: dormant in
  driver, returns when scan-emitter TRs are attached.

## Original quality-of-life gaps

- **G10** `smoke_tests` BINARY DIFF regression — still open; tracked
  in memory entry `project_smoke_test_regression_20260423`.
- **G11** clangd / Eigen include path — dev UX; unchanged.
- **G12** `frame_paths_` tracking — dropped; downstream consumers
  iterate `output_dir/frame_*`.
- **G13** `#include "GromacsRunContext.h"` stale-comments sweep —
  named files no longer contain matches (see drift note).
