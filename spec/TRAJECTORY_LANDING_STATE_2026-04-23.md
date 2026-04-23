# Trajectory Refactor — Landing State (2026-04-23)

Companion to `spec/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md`. Written at the
end of the landing session. Captures what lives in the tree, what the
clean-room sandbox exercise revealed, what is residue vs pattern, and what
the next session starts from.

## What is in the tree

Trajectory-scope framework under `src/`:

- `TrajectoryResult.{h,cpp}` + `SelectionEmittingTrajectoryResult.h`
- `TrajectoryAtom.h` (BS fields only; anti-pattern warning in header block)
- `TrajectoryProtein.{h,cpp}`
- `DenseBuffer.h`
- `BsWelfordTrajectoryResult.{h,cpp}` — first and only concrete
- `Trajectory.{h,cpp}` — 5-phase Run (with stride support)
- `RunConfiguration.{h,cpp}` — three named static factories
- `RunContext.{h,cpp}` — class, SetStride / SetAimnet2Model / AttachExtra
- `TrajectoryMoments.h` — free functions `MomentsUpdate`, `MinMaxUpdate`,
  `MomentsStd`; factors Welford ritual out of Compute bodies
- `GromacsFrameHandler.{h,cpp}` — takes TrajectoryProtein now; still
  GROMACS-named per §9 format-specific

Old `GromacsProtein.{h,cpp}`, `GromacsProteinAtom.h`, `GromacsRunContext.{h,
cpp}`, `GromacsFinalResult.{h,cpp}` moved to `learn/bones/` with
`GromacsProtein_MIGRATION_NOTE.md`.

`nmr_extract.cpp`: RunTrajectory uses Trajectory::Run; RunAnalysis still
uses GromacsFrameHandler + AnalysisWriter directly (explicitly unchanged
this session — AnalysisWriter dissolution is separate work).

`tests/test_gromacs_streaming.cpp` rewritten: three tests pass —
TrajectoryBuildAndScan, BsWelfordAttachAndFinalize, TrajectoryStrideSkipsFramesCorrectly.
A fourth test TrajectoryRunDrivesLoop exists but requires external
`LD_LIBRARY_PATH` to load PyTorch's bundled CUDA runtime; see NVRTC note.

Build: `make -j nmr_shielding gromacs_streaming_tests` clean. `unit_tests`
82/82. Three streaming tests pass with LD_LIBRARY_PATH export.

## What the sandbox exercise revealed

A clean-room sandbox at `analysis-speculative/trajectory_framework_sketch/`
implemented the framework from spec only, reading order CONSTITUTION →
PATTERNS → OBJECT_MODEL → USE_CASES → WIP_OBJECT_MODEL, with no reference
to `src/`. Sorting the sandbox against `src/` against the spec:

### Core pattern — converged between the independent implementations

Both agents landed on the same shape for these. They are pattern-essential,
not residue:

- Attach-map + attach-order-vector duality on TrajectoryProtein. The
  tension ("dispatch needs order, callers need typed-lookup") is inherent
  to what the spec asks; both implementations hit it and stored both.
- TrajectoryResult base with four virtuals + optional serialization hooks.
- Singleton-per-class AttachResult.
- 5-phase Trajectory::Run, phases 3 and 4 iterating
  `tp.ResultsInAttachOrder()` explicitly for Compute.
- SelectionEmittingTrajectoryResult as mixin, collected via `dynamic_cast`
  in Phase 5.
- Private ctor + `friend TrajectoryProtein` on TrajectoryAtom; typed-output
  fields only; no Welford/DeltaTracker instances anywhere in the struct.
- `DenseBuffer<T>` template; AdoptDenseBuffer ownership transfer.
- `BsWelfordTrajectoryResult` pattern: `Create(tp)` factory, `Compute`
  updates TrajectoryAtom fields per-atom-per-frame, `Finalize` converts M2
  to std.

### Src/ residue — in `src/` but not in the sandbox; not pattern-essential

- **`bool dispatch` parameter on `GromacsFrameHandler::Next`** — CLEANED
  UP end-of-session. Handler.Next now takes only `(opts, output_dir)`
  and never dispatches TrajectoryResults (role leak eliminated).
  Trajectory::Run Phase 3 + Phase 4 call `tp.DispatchCompute(conf, idx,
  time)` after the per-frame pipeline returns. Callers outside
  Trajectory::Run (e.g. tests that exercise Compute directly) also call
  `tp.DispatchCompute` — it is the named operation the trajectory
  protein performs on itself per frame, not a wrapper to be inlined.
- **Per-frame NPY write inside `ProcessFrame`** — NOT CLEANED this
  session. Handler still writes `frame_NNNN/*.npy` when given a non-empty
  output_dir. Belongs in Trajectory::Run or a writer. Deferred to the
  AnalysisWriter dissolution work, which owns the writer-layer decision.
- **Attach-time partial dep check comment** — NOT CLEANED. AttachResult
  does singleton only (correct per §12 RESOLVED) but the comment still
  suggests attach-time pre-check is possible. Stale; strip in next
  session's read-and-cut pass.
- **`MeanAtAtom` / `StdAtAtom` on `BsWelfordTrajectoryResult`** — CLEANED
  UP end-of-session. Methods removed from both header and cpp. Callers
  read `tp.AtomAt(i).bs_t0_mean` etc. directly, matching the private-ctor
  + typed-fields pattern on TrajectoryAtom.

### Configure/run tangles — spec defensible, but defensibility is the weakest link

Both implementations reach for these because the spec's §6 pseudocode
shows the shape. Neither agent questioned the premise:

- **Factory lambda list on RunConfiguration**. Both have
  `std::vector<std::function<std::unique_ptr<TrajectoryResult>(const TrajectoryProtein&)>>`.
  The lambdas defer construction only because RunConfiguration is built
  standalone, before tp exists. A free-function model
  (`AttachPerFrameExtractionSet(tp, ctx)`) would eliminate the list.
  Spec commits to RunConfiguration as typed-first-class; spec does not
  argue WHY standalone. Worth examining.
- **EDR plumbing via handler back-pointer**. My `BindTrajectory(traj)`
  passes handler a pointer so it can call `traj.EnergyAtTime(time)` before
  each `OperationRunner::Run`. Spec §5 EDR-timing subsection sanctions
  this. Sandbox elides the decision by stubbing the handler. Asymmetric
  ownership the rest of the object model doesn't need; could be solved
  differently, but spec-committed.
- **RunContext existing as a distinct entity from Trajectory**. Both
  keep it; spec describes it separately; the reason for the split is
  not examined in spec text.

### Spec-internal contradictions surfaced by the sandbox

These are real issues in the spec itself, not implementation choices:

- **Exceptions vs return codes**. PATTERNS.md §9 / "Exception hierarchies"
  forbid `throw` in calculator/result/pipeline code. WIP §5's concrete
  `Trajectory::Run` body throws `std::runtime_error` and
  `std::logic_error`. Both implementations follow each section's own code
  (AttachResult returns bool per §3, Trajectory::Run throws per §5) but
  the contradiction is unresolved at the spec level.
- **`MeanAtAtom` / `StdAtAtom` convenience methods on TrajectoryResult
  subclasses**. §4's worked example shows these with no tp parameter and
  no back-pointer. There is no honest implementation. Both agents kept
  the spec signature and returned 0.0 with a source note. Methods should
  be removed from the spec or the spec should show how they're supposed
  to work.

## Out of this session's scope

Listed in `TRAJECTORY_REFACTOR_GAPS_2026-04-23.md`. Summary:

- G4 (selection-emitting concrete Results) — the scan+extract two-pass
  workflow is dormant until these land.
- G5 (~40 other Welford migrations per Appendix F catalogue) — the grind.
- G6 (AnalysisWriter dissolution) — THIS IS THE LOAD-BEARING GAP. Current
  `RunAnalysis` uses AnalysisWriter as its output surface, bypassing
  TrajectoryProtein::WriteAllFeatures entirely. Production analysis is
  on the old pattern. Spec-session meeting-of-minds was that AnalysisWriter
  should dissolve into a family of `*TimeSeriesTrajectoryResult` classes,
  each owning a DenseBuffer and emitting via `WriteH5Group` — exactly the
  BsWelford pattern at scale. This session did not land that work. It
  must land before `RunAnalysis` is considered refactor-complete.
- NVRTC bundled library rpath — the CUDA JIT dependency on
  `libnvrtc-builtins.so.13.0` still requires external LD_LIBRARY_PATH
  export. Documented in memory. Real fix (patchelf-in-CMake or CTest
  ENVIRONMENT wrapper) still deferred.
- Smoke test BINARY DIFF regression — pre-dates this session, not
  investigated.

## Reading order for the next session

Start of next session should load, in order:

1. `spec/TRAJECTORY_LANDING_STATE_2026-04-23.md` (this file)
2. `spec/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md` (the gap inventory)
3. `analysis-speculative/trajectory_framework_sketch/README.md` (sandbox
   read — which tangles are core vs configure/run)
4. `spec/WIP_OBJECT_MODEL.md` §0 non-regression contract, §3 anti-patterns
5. Only then: whichever src/ file is being worked on

If the next session starts with AnalysisWriter dissolution, also read
`spec/ANALYSIS_TRAJECTORY_2026-04-14.md` to see what H5 groups
AnalysisWriter currently emits — each becomes one `*TimeSeriesTrajectoryResult`.

## Spec issues to resolve before more code is written

From the sandbox exercise:

1. Pick one: `throw` in Run, bool in AttachResult, or consistent return
   codes throughout. Either way, update PATTERNS.md + WIP §5 so they
   agree.
2. Remove `MeanAtAtom` / `StdAtAtom` from WIP §4. User directive
   end-of-session: "they can go away". Sandbox and src/ both agree they
   don't work without tp threading; the convenience was false. Spec §4
   should drop them from the BsWelford worked example; callers read
   TrajectoryAtom fields directly.
3. Decide whether RunConfiguration-as-standalone-typed-object is the
   right shape, or whether free-function attach-sequences would do.
   Affects the factory list question.

These are spec decisions, not implementation. Ideally settled before G5
grinding starts — the 40 Welford migrations will each need to know the
answer to (1) and (2).

## Session provenance

- Pre-session tag: `master @ 83b62b6` (spec: trajectory-scope WIP object
  model + 2026-04-22 design pass)
- Landing session date: 2026-04-23
- Landing state: uncommitted; tree dirty
- End-of-session cleanup: `bool dispatch` removed from handler.Next
  (role leak fix — handler no longer dispatches);
  `tp.DispatchCompute(conf, idx, time)` kept as the named operation
  on TrajectoryProtein, called from Trajectory::Run Phase 3 + 4;
  `MeanAtAtom`/`StdAtAtom` removed from BsWelfordTrajectoryResult
  (per user directive — false convenience).
- End-of-session course-correction: mid-cleanup I also deleted
  `DispatchCompute` on the "one-line wrapper" reasoning and replaced
  it with inline for loops at each call site. User called this out as
  mistaking the object model's organising principle for mechanical
  wrapping. Named operations on entities ARE the organisation; they
  hold the shape the model is trying to lay down. `DispatchCompute`
  restored. Same principle preserves `FinalizeAllResults`, `WriteH5`,
  `WriteFeatures` on TrajectoryProtein — each is a named operation,
  not a wrapper.
- Test status at landing: `unit_tests` 82/82,
  `gromacs_streaming_tests` 3/3 fast tests (Build+Scan,
  BsWelfordAttachAndFinalize, TrajectoryStrideSkipsFramesCorrectly)
  green after dispatch cleanup. `TrajectoryRunDrivesLoop` (long AIMNet2
  full-run) not re-run post-cleanup but the code paths it exercises are
  unchanged by the cleanup.
