# Trajectory Refactor — Landing State (historical, 2026-04-23)

**HISTORICAL RECORD.** Session-end note from the trajectory-scope
refactor landing. The current shape of the trajectory-scope object
model and its patterns lives in `OBJECT_MODEL.md` (trajectory-scope
entities section) and `PATTERNS.md §§13-18`; when those disagree with
anything here, they are the truth. Class names, phase counts, and
file pointers in this file reflect the 2026-04-23 session-end tree
at the moment of writing, not the current tree — specifically,
`RunContext.{h,cpp}` and `SelectionEmittingTrajectoryResult.h` were
dissolved after this session, the mixin-plus-`dynamic_cast`-sweep
approach was abandoned in favour of direct `RecordBag::Push`,
`Trajectory::Run` is eight phases not five, and
`spec/WIP_OBJECT_MODEL.md` was renamed to
`spec/pending_include_trajectory_scope_2026-04-22.md`.

## What the session landed

Trajectory-scope framework replacing `GromacsProtein` /
`GromacsProteinAtom` / `GromacsRunContext` / `GromacsFinalResult`
(all moved to `learn/bones/` with a `GromacsProtein_MIGRATION_NOTE.md`).
New classes: `TrajectoryProtein`, `TrajectoryAtom`,
`TrajectoryResult`, `Trajectory`, `RunConfiguration`, `DenseBuffer`,
plus `TrajectoryMoments.h` helpers and the first concrete TR,
`BsWelfordTrajectoryResult`. `GromacsFrameHandler` kept in its
pure-reader role (format-specific name preserved deliberately).

`nmr_extract.cpp::RunTrajectory` rewired to `Trajectory::Run`.
`RunAnalysis` left on the legacy AnalysisWriter path; its dissolution
into per-Result `WriteH5Group` emitters was named as follow-up work.

At session-end: `unit_tests` 82/82, three of the rewritten
`gromacs_streaming_tests` green
(`TrajectoryBuildAndScan`, `BsWelfordAttachAndFinalize`,
`TrajectoryStrideSkipsFramesCorrectly`). A fourth streaming test
(`TrajectoryRunDrivesLoop`) required external `LD_LIBRARY_PATH` for
PyTorch's bundled CUDA runtime; NVRTC rpath fix was deferred.

## Sandbox exercise

A clean-room sandbox at
`analysis-speculative/trajectory_framework_sketch/` implemented the
framework from spec alone, without reading `src/`. Convergence with
the `src/` implementation confirmed several shapes as
pattern-essential rather than residue: attach-map + attach-order
duality on TrajectoryProtein, singleton-per-class AttachResult,
private-ctor TrajectoryAtom with typed-output fields only,
`DenseBuffer<T>` + ownership transfer at Finalize, and the
`Create(tp) → Compute → Finalize` TR shape with `Compute` updating
TrajectoryAtom fields per atom per frame.

The sandbox also surfaced contradictions in the spec itself
(exceptions-vs-return-codes between PATTERNS.md §9 and WIP §5,
`MeanAtAtom`/`StdAtAtom` convenience methods that could not be
honestly implemented without a `tp` parameter). The former was
settled by making `Trajectory::Run` return `Status`; the latter by
removing the convenience methods and reading TrajectoryAtom fields
directly.

## End-of-session cleanups

- `bool dispatch` parameter removed from `GromacsFrameHandler::Next`;
  handler is now a pure reader (`Open` / `ReadNextFrame` / `Skip` /
  `Reopen`) and dispatch happens in `Trajectory::Run`.
- `MeanAtAtom` / `StdAtAtom` removed from
  `BsWelfordTrajectoryResult`; callers read `tp.AtomAt(i).bs_t0_mean`
  etc. directly.
- `DispatchCompute`, `FinalizeAllResults`, `WriteH5`, `WriteFeatures`
  preserved on `TrajectoryProtein` as named operations on the entity
  (not wrappers to be inlined).

## Open items

For current open items, pending user decisions, and known-wrong doc
content, see:

- `spec/pending_decisions_20260423.md` — items requiring user choice
  (AllWelfords revival, `ChiRotamerSelection` attachment,
  `FullFatFrameExtraction` MOPAC deps)
- `spec/doc_wrongness_20260423.md` — observational audit of doc/code
  contradictions
- `OBJECT_MODEL.md` + `PATTERNS.md §§13-18` — authoritative current
  object-model shape and patterns
- Memory entry `reference_nvrtc_rpath_fix` — cause + CMake fix for
  the NVRTC bundled-library rpath issue

## Session provenance

- Pre-session tag: `master @ 83b62b6` (spec: trajectory-scope WIP
  object model + 2026-04-22 design pass).
- Landing session date: 2026-04-23.
- Landing state at time of writing: uncommitted; tree dirty.
