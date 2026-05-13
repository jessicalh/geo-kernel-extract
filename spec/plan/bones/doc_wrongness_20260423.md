# Documentation wrongness findings — 2026-04-23

This is an observational audit only. No replacement text is proposed; each entry pairs a claim in docs/comments with what the code actually does, citing evidence by `path:line`.

## Findings

### 1. `TrajectoryProtein::AttachResult` docstring claims dependency check that does not exist

- **Where**: `src/TrajectoryProtein.h:141-146`, `src/TrajectoryProtein.cpp:132-142`
- **Claim**: Header docstring says `AttachResult` "Checks singleton (one per type) and dependencies (every declared type_index must already be attached as another TrajectoryResult OR must be a ConformationResult type — the latter is validated by Trajectory::Run against the RunConfiguration's per-frame set)." The cpp comment block directly above the implementation repeats the dependency-check claim.
- **Reality**: `src/TrajectoryProtein.cpp:143-164` performs only the singleton check and the `results_attach_order_.push_back`. No `Dependencies()` call, no loop over declared deps, no diagnostic about missing deps. The actual dep check happens in `Trajectory::Run` Phase 4 at `src/Trajectory.cpp:134-143`.
- **Note**: `OBJECT_MODEL.md:1429-1431` correctly documents this (`"Dependency validation is at Trajectory::Run Phase 4, not here."`). The header/cpp comment blocks are inconsistent with OBJECT_MODEL and with the implementation.

### 2. `TrajectoryResult::Dependencies` docstring claims attach-time check

- **Where**: `src/TrajectoryResult.h:50-55`
- **Claim**: "Attach-time check fires if any declared dependency is missing."
- **Reality**: As in Finding 1, `TrajectoryProtein::AttachResult` does not check dependencies. Check fires in `Trajectory::Run` Phase 4 at `src/Trajectory.cpp:134-143`.

### 3. `RunConfiguration.h` docstring claims "this session populates only BsWelfordTrajectoryResult factory"

- **Where**: `src/RunConfiguration.h:26-30`
- **Claim**: "This session populates only the BsWelfordTrajectoryResult factory in PerFrameExtractionSet (and in ScanForDftPointSet, since BS runs per-frame there too). The other TrajectoryResult factories remain empty placeholders until follow-up sessions migrate their fields."
- **Reality**: `src/RunConfiguration.cpp:169-192` attaches six factories in `PerFrameExtractionSet`: `BsWelford`, `BsShieldingTimeSeries`, `BsAnomalousAtomMarker`, `BsT0Autocorrelation`, `BondLengthStats`, `PositionsTimeSeries`.
- **Note**: Session-narrative phrasing ("this session") is also flagged in Finding 22.

### 4. `RunConfiguration.cpp` "two worked examples" comment contradicts the factory list it sits in front of

- **Where**: `src/RunConfiguration.cpp:157-168`
- **Claim**: "The two worked examples currently exercise the two canonical shapes: BsWelfordTrajectoryResult … and PositionsTimeSeriesTrajectoryResult …. The rest of Appendix F lands in follow-up sessions — one AddTrajectoryResultFactory line per class."
- **Reality**: Six factories are registered immediately after the comment at `src/RunConfiguration.cpp:169-192`.

### 5. `RunConfiguration.h` claims AIMNet2 check is Phase 2, actually Phase 4

- **Where**: `src/RunConfiguration.h:71-81`
- **Claim**: "True means `Trajectory::Run` Phase 2 will throw if the corresponding pointer in the resolved RunOptions is null."
- **Reality**: `src/Trajectory.cpp:145-150` does the `RequiresAimnet2` check in Phase 4 (not Phase 2), and returns `kConfigRequiresAimnet2` (a `Status`, not an exception). Phase 2 at `src/Trajectory.cpp:86-98` is reading frame 0 + `tp.Seed`.

### 6. `RunConfiguration.h` describes `ScanForDftPointSet` as emitting selections; actual factory attaches only `BsWelfordTrajectoryResult`

- **Where**: `src/RunConfiguration.h:6-7`
- **Claim**: "ScanForDftPointSet — cheap per-frame, selection-emitting TrajectoryResults for DFT pose choice"
- **Reality**: `src/RunConfiguration.cpp:54-86` attaches only `BsWelfordTrajectoryResult` (not a selection emitter). Comment at cpp lines 73-79 acknowledges "This session: only BsWelford" and lists rotamer / rmsd / pose-coordinator as "follow-up sessions populate the rest."

### 7. `nmr_extract.cpp` comment describes Run as 5-phase with `RunContext`

- **Where**: `src/nmr_extract.cpp:227-230`
- **Claim**: "Trajectory owns file paths, EDR, and drives the 5-phase Run(). RunConfiguration describes what the run does; RunContext holds the caller's choices."
- **Reality**: `src/Trajectory.cpp:50` and `src/Trajectory.h:14-44` say eight phases. No `RunContext` exists in `src/`; the files are in `learn/bones/` (confirmed absent from `src/` at this date).

### 8. `RunConfiguration.cpp` references `RunContext::SetSelectedFrameIndices`

- **Where**: `src/RunConfiguration.cpp:205`
- **Claim**: "ScanForDftPointSet run, or RunContext::SetSelectedFrameIndices in a future iteration."
- **Reality**: `RunContext` is not in `src/`; the class was dissolved (see Finding 7).

### 9. `nmr_extract.cpp` references `RunContext` as holding caller choices

- **Where**: `src/nmr_extract.cpp:230`
- **Claim**: "RunConfiguration describes what the run does; RunContext holds the caller's choices."
- **Reality**: `Trajectory` owns frame record, `RunConfiguration` owns run shape, `Session` owns runtime resources. No `RunContext` class.

### 10. Stale `WIP_OBJECT_MODEL.md` file references throughout `src/`

- **Where** (multiple):
  - `src/TrajectoryProtein.h:7` ("per spec/WIP_OBJECT_MODEL.md §3")
  - `src/TrajectoryAtom.h:34` ("see spec/WIP_OBJECT_MODEL.md §3")
  - `src/TrajectoryResult.h:24, 26` (§3 + §4)
  - `src/RunConfiguration.h:30` ("§6 + Appendix F")
  - `src/BsWelfordTrajectoryResult.h:16` ("§4 worked example")
  - `src/BondLengthStatsTrajectoryResult.h:7` ("§3 option (b)")
  - `src/DenseBuffer.h:18, 20` (Appendix C, §4)
- **Claim**: Each reference cites `spec/WIP_OBJECT_MODEL.md`.
- **Reality**: File was renamed to `spec/pending_include_trajectory_scope_2026-04-22.md` per `spec/INDEX.md:105`. No `WIP_OBJECT_MODEL.md` exists at `spec/` or root.

### 11. `FinalizeProtein()` referenced after rename to `Seed()`

- **Where**:
  - `src/TrajectoryProtein.h:66-67` ("FinalizeProtein() completes construction from first frame positions")
  - `src/TrajectoryProtein.cpp:23-24` ("The protein is not finalized here — FinalizeProtein must be called by the frame handler…")
  - `spec/ENSEMBLE_MODEL.md:49, 59, 148, 272` ("FinalizeProtein()…")
- **Claim**: Method named `FinalizeProtein()` is the completion-from-first-frame call.
- **Reality**: `src/TrajectoryProtein.cpp:69-89` defines `Seed()`; no `FinalizeProtein` method exists on `TrajectoryProtein`. OBJECT_MODEL.md:1419 correctly names `Seed`.

### 12. `TrajectoryAtom.h` "Three coexisting shapes" claims Pattern A is present

- **Where**: `src/TrajectoryAtom.h:41-64`
- **Claim**: "Three coexisting shapes hold per-atom trajectory data: Pattern A — typed struct vectors for known-shape per-source data. …Lands as classes like RingNeighbourhoodTrajectoryStats as the catalog fills in. Pattern B — typed accumulator fields, one writer per field…. Pattern C — the per-atom event bag."
- **Reality**: The struct body at `src/TrajectoryAtom.h:77-123` contains only Pattern B (bs_t0_* / bs_t2mag_* / bs_t0_delta_* fields at lines 86-111) and Pattern C (`RecordBag<AtomEvent> events` at line 119). No Pattern A fields are declared.
- **Note**: `OBJECT_MODEL.md:1474-1478` correctly labels Pattern A as "Aspirational; no current fields." The header itself does not say Pattern A is aspirational.

### 13. `Trajectory.h` ctor comment places EDR preload at "Phase 3+"; body says preload happens before all phases

- **Where**: `src/Trajectory.h:41-44` vs `src/Trajectory.cpp:35-43`
- **Claim**: `Trajectory.h:41-44` — "EDR frames are preloaded in the constructor so Phase 6 already has them available." (consistent with body.) `src/Trajectory.cpp:35-36` comment inside ctor says "Preload EDR now; available through Phase 3+ for per-frame energy lookup via EnergyAtTime."
- **Reality**: EDR is consumed in Phase 6 (`src/Trajectory.cpp:169`) and Phase 7 (line 210). Phase 3 attaches `TrajectoryResult`s (`src/Trajectory.cpp:101-122`) and does not touch `EnergyAtTime`. "Phase 3+" is accurate in the sense "any time from Phase 3 onward", but the header/body pair disagree on which phase first uses the data (6 vs 3).

### 14. `ScanForDftPointSet` factory comment claims "BiotSavart runs because ring proximity matters for selection" but required ConformationResult list contradicts the "cheap set" framing

- **Where**: `src/RunConfiguration.cpp:43-86`
- **Claim**: Class-level comment (lines 45-52) says "Cheap per-frame: no MOPAC, no vacuum Coulomb, no APBS, no AIMNet2. … BiotSavart runs because ring proximity matters for selection; BsWelford aggregates it."
- **Reality**: `src/RunConfiguration.cpp:66-71` requires `GeometryResult`, `SpatialIndexResult`, `EnrichmentResult`, `DsspResult`, `BiotSavartResult`, `SasaResult`. HM and McConnell, listed in the "Tier 1 — cheap set (keep Geometry / SpatialIndex / Enrichment / DSSP / BS / HM / McConnell / SASA)" narrative in `OBJECT_MODEL.md:1596-1600`, are NOT required here. Cross-doc inconsistency between cpp factory and OBJECT_MODEL description of the same config.

### 15. `spec/TRAJECTORY_LANDING_STATE_2026-04-23.md` references `SelectionEmittingTrajectoryResult.h`, `RunContext.{h,cpp}`, 5-phase Run, phases 3+4 for compute

- **Where**: `spec/TRAJECTORY_LANDING_STATE_2026-04-23.md:12, 17, 19, 45, 58, 60-61, 74, 110, 161`
- **Claims** (various):
  - Line 12: "`SelectionEmittingTrajectoryResult.h`" in the tree
  - Line 17: "5-phase Run (with stride support)"
  - Line 19: "`RunContext.{h,cpp}` — class, SetStride / SetAimnet2Model / AttachExtra"
  - Line 45: reading order includes `WIP_OBJECT_MODEL`
  - Line 58: "5-phase Trajectory::Run, phases 3 and 4 iterating `tp.ResultsInAttachOrder()` explicitly for Compute."
  - Line 60-61: "SelectionEmittingTrajectoryResult as mixin, collected via `dynamic_cast` in Phase 5."
  - Line 74: "Trajectory::Run Phase 3 + Phase 4 call `tp.DispatchCompute(conf, idx, time)` after the per-frame pipeline…"
  - Line 110: "RunContext existing as a distinct entity from Trajectory"
  - Line 161: references `spec/WIP_OBJECT_MODEL.md` §0 / §3
- **Reality**: `SelectionEmittingTrajectoryResult.h` and `RunContext.{h,cpp}` do not exist in `src/`. `Trajectory::Run` is eight phases (`src/Trajectory.cpp:50`, `src/Trajectory.h:14-44`). `DispatchCompute` is called in Phase 6 (line 185) and Phase 7 (line 226). No mixin, no `dynamic_cast` sweep — `src/Trajectory.cpp:241-243` says "Selections are emitted directly … No dynamic_cast sweep here." The `WIP_OBJECT_MODEL.md` file was renamed (see Finding 10).
- **Note**: This file's surface-level purpose is historical ("Landing State"), but it is referred to as useful reading from current docs (`spec/INDEX.md:106`) and is dated 2026-04-23, matching today; every stale reference inside it will read as current to a reader who trusts the date.

### 16. `spec/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md` leads with stale class list and a stride "gap" that is closed in code

- **Where**: `spec/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md:3-7, 17-30, 49-75, 198-209, 64-77`
- **Claims**:
  - Lines 3-7 intro names `RunContext` and `SelectionEmittingTrajectoryResult` as part of the landed pattern.
  - Lines 17-30 (G1): "Production `RunAnalysis` used `STRIDE=2` … `Trajectory::Run`'s Phase 4 loop currently iterates every frame (stride=1)."
  - Lines 24-25 (G1 shape-of-fix): "`RunContext` gains `size_t stride_ = 1` field with setter `SetStride(size_t s)`."
  - Lines 49-62 (G3): `RunContext::SetSelectedFrameIndices` pattern.
  - Lines 64-77 (G4): "SelectionEmittingTrajectoryResult has zero implementers."
  - Lines 198-209 (G13): "Library is not `#include \"GromacsRunContext.h\"` clean in comments … stale comments at `src/JobSpec.h:32`, `src/OperationRunner.h:62, 74`, `src/GromacsEnergyResult.h:90`, `src/GromacsEnergyResult.cpp:9`, `src/BondedEnergyResult.h:57`."
- **Reality**:
  - Stride IS implemented: `src/RunConfiguration.h:88-89` (`SetStride`/`Stride`) on `RunConfiguration` (not on `RunContext`), consumed at `src/Trajectory.cpp:198-204`. G1 as written is stale.
  - `RunContext` absent from `src/` (see Finding 7).
  - `SelectionEmittingTrajectoryResult` absent. `ChiRotamerSelectionTrajectoryResult` (`src/ChiRotamerSelectionTrajectoryResult.{h,cpp}`) exists and pushes selections directly via `traj.MutableSelections()` (no mixin). `OBJECT_MODEL.md:1522, 2398` lists it.
  - G13-listed stale comments: grep of `src/JobSpec.h`, `src/OperationRunner.h`, `src/GromacsEnergyResult.h`, `src/GromacsEnergyResult.cpp`, `src/BondedEnergyResult.h` returns zero `GromacsRunContext`/`GromacsProtein`/`AccumulateFrame` hits. G13's specific line references no longer match.

### 17. `spec/INDEX.md` file/test counts and date mismatch

- **Where**: `spec/INDEX.md:131, 133, 140`
- **Claim**: Line 131: "src/ — 74 headers, 62 .cpp files." Line 133: "tests/ — GTest suite (287 tests, 40 test files, …)." Line 140: "Current Status (2026-04-15)".
- **Reality**: `src/*.h` = 88 files, `src/*.cpp` = 70. `tests/test_*.cpp` = 47 files (49 `.cpp` total). Today is 2026-04-23; the "Current Status (2026-04-15)" date header is stale by itself, independent of whether content is still true.

### 18. `spec/INDEX.md` Current Status section and "Two-pass trajectory streaming" description references dissolved classes

- **Where**: `spec/INDEX.md:156-159`
- **Claim**: "Two-pass trajectory streaming: GromacsProtein (adapter + accumulators) + GromacsFrameHandler (streaming XTC reader, PBC fix, frame lifecycle) + GromacsRunContext (bonded params, EDR energy frames, cursor state)."
- **Reality**: `GromacsProtein` and `GromacsRunContext` are in `learn/bones/`, no longer in `src/`. The live design is `TrajectoryProtein` + `GromacsFrameHandler` + `Trajectory` (not `GromacsRunContext`). `RunAnalysis` itself is stubbed (`src/nmr_extract.cpp:301-310`).

### 19. `spec/INDEX.md` lists `ENSEMBLE_MODEL.md` as the trajectory model with a mid-April revision date

- **Where**: `spec/INDEX.md:94`
- **Claim**: "ENSEMBLE_MODEL.md — GromacsProtein trajectory-based model (revised 2026-04-12, replaces old EnsembleConformation design)"
- **Reality**: `spec/ENSEMBLE_MODEL.md:1-26` header itself says the document is partially historical and the live model is in `OBJECT_MODEL.md` + `PATTERNS.md §§13-18`. The line in INDEX.md still points at it as the trajectory model without flagging the partial obsolescence inside.

### 20. `spec/ENSEMBLE_MODEL.md` body still describes `GromacsProtein` / `AccumulateFrame` / `SelectFrames` / `WriteCatalog` / `AllWelfords` as live API

- **Where**: `spec/ENSEMBLE_MODEL.md:30-346` (continuous body of the document)
- **Claim**: The "Architecture: GromacsProtein pattern" section and following sections describe the build paths, accumulation methods, two-pass scan/extract, `GromacsFinalResult`, and `AllWelfords()` as the working design.
- **Reality**: `GromacsProtein.{h,cpp}`, `GromacsProteinAtom.h`, `GromacsRunContext.{h,cpp}`, `GromacsFinalResult.{h,cpp}`, `AnalysisWriter.{h,cpp}` are all in `learn/bones/`. `ENSEMBLE_MODEL.md:3-15` flags this at the top, but the body is not annotated or struck through.

### 21. `PATTERNS.md` "AllWelfords pattern" section describes live behaviour of moved classes

- **Where**: `PATTERNS.md:1024-1048`
- **Claim**: "GromacsProteinAtom has ~45 Welford accumulators. WriteCatalog (CSV) and WriteH5 both need the column list. … The fix: `AllWelfords()` returns a `vector<NamedWelford>` with `{name, pointer-to-Welford}` pairs…. This is not ideal — the accumulation in AccumulateFrame is still manual…. The AccumulateFrame code is write-once (change it when you add a new calculator)…"
- **Reality**: `GromacsProteinAtom`, `AllWelfords`, `AccumulateFrame`, `WriteCatalog` all in `learn/bones/`; `TrajectoryAtom` replaces `GromacsProteinAtom` with one-writer-per-field discipline and no central enumeration. The pattern is presented as a lesson learned without flagging that the underlying classes are no longer live.

### 22. Session-narrative comments in source code

- **Where**:
  - `src/RunConfiguration.h:26` ("This session populates only…")
  - `src/RunConfiguration.h:77-79` ("See user feedback 2026-04-23…")
  - `src/RunConfiguration.cpp:73` ("This session: only BsWelford.")
  - `src/RunConfiguration.cpp:74` ("Follow-up sessions populate the rest per Appendix F:")
  - `src/RunConfiguration.cpp:163` ("The rest of Appendix F lands in follow-up sessions…")
  - `src/RunConfiguration.cpp:217, 221` ("follow-up session", "in this session")
  - `src/JobSpec.h:30` ("Fleet removed 2026-04-12.")
  - `src/JobSpec.cpp:121` ("--fleet mode removed (2026-04-12).")
  - `src/nmr_extract.cpp:219` ("RunFleet (--fleet mode) removed 2026-04-12.")
  - `src/nmr_extract.cpp:232-234` ("This session: single-pass Trajectory::Run with PerFrameExtractionSet. The former two-pass scan+extract pattern comes back when scan-mode TrajectoryResults land and push to traj.MutableSelections() (future session).")
  - `src/nmr_extract.cpp:283-287` ("STUBBED 2026-04-23. The prior implementation …")
  - `src/RuntimeEnvironment.h:15` ("Live surface (2026-04-05):")
  - `src/ChargeAssignmentResult.h:43` / `.cpp:186` ("No-arg Compute removed 2026-04-03.")
  - `src/ChargeSource.h:162` / `.cpp:290` ("StubChargeSource — REMOVED 2026-04-03.")
  - `src/OrcaRunLoader.h:17` ("As of 2026-04-02, prmtop regenerated via tLeap …")
  - `src/NamingRegistry.cpp:155-250` (multiple date-tagged "DISCOVERED 2026-04-20, DEFERRED until fleet-wide vetting" / "scheduled to land around 2026-04-27" / "memory:project_charmm_iupac_gaps_2026-04-20" blocks)
- **Claim**: Each of these comments describes a specific session, a scheduled future date, or a "this session" scope.
- **Reality**: Flagged as stale by nature per the sweep's method ("session-narrative comments rot by their nature").

### 23. Session-narrative in design docs

- **Where**:
  - `spec/DIRECTORY_SET.md:5` ("Last known good state before this session's changes.")
  - `spec/TRAJECTORY_EXTRACTION.md:39` ("Built this session")
  - `spec/SESSION_STATE_20260412.md:3` ("What was built this session")
  - `spec/PLANNED_CALCULATORS_2026-04-22.md:289` ("next session…")
  - `spec/MICROSECOND_MD_HARVESTER_2026-04-22.md:8` ("Not for the next session. The next session is Session 1 drafting…")
  - `spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md:6` ("The next session should re-read, not rubber-stamp.")
  - `spec/TRAJECTORY_LANDING_STATE_2026-04-23.md:6, 31, 131, 150, 153, 155, 164` ("this session", "next session")
  - `spec/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md:231` ("Immediate (next session, ~1 sitting)")
- **Claim**: Each reference scopes content to a specific session.
- **Reality**: Flagged as stale by nature.

### 24. `OBJECT_MODEL.md` claim about where ConformationResult computation happens

- **Where**: `OBJECT_MODEL.md:1319-1323`, `OBJECT_MODEL.md:1325-1330`
- **Claim**: "The framework calls AttachResult() on the ProteinConformation, which checks dependencies via the type_index vector." (correct.) Lines 1329-1330 add: "The attach method is where computation happens: the result reads from prior results and stores its own properties."
- **Reality**: `src/ProteinConformation.cpp:22-59` (`AttachResult`) does singleton check + dependency check + store — no computation. Computation is done in each subclass's static `Compute(conf, …)` factory before `AttachResult` is called (`PATTERNS.md:105-111`).
- **Note**: Two lines apart in the same section disagree. Line 1319-1321 frames attach as dependency check; line 1329 frames attach as where computation happens.

### 25. `CLAUDE.md` "Current state (as of 2026-04-17)" is stale

- **Where**: `CLAUDE.md:273-320`
- **Claim**: Section is dated 2026-04-17 and reports fleet progress ("200 of 685 proteins H5-complete"), viewer + reader session dates, calibration DFT completion 2026-04-18, etc.
- **Reality**: Today is 2026-04-23. The 2026-04-17 label is by itself stale; some dated facts within (calibration DFT 2026-04-18) post-date the header. Header date does not match its contents.

### 26. `OBJECT_MODEL.md` "as of landing" section header

- **Where**: `OBJECT_MODEL.md:1512`
- **Claim**: "Known TrajectoryResult types (as of landing)" introduces a table.
- **Reality**: Phrase is session-narrative by nature ("landing" refers to a specific prior session). Flagged per the sweep method, not per content discrepancy — the table rows themselves match `src/` implementations.

### 27. `spec/INDEX.md` externs list omits HighFive and xdrfile

- **Where**: `spec/INDEX.md:134`
- **Claim**: "extern/ — header-only libraries (nanoflann, sphericart)"
- **Reality**: `ls extern/` shows `HighFive/`, `nanoflann.hpp`, `sphericart*.{h,hpp}`, `xdrfile/`, plus `cuda_*`, `macros.hpp`, `templates*.hpp`.

### 28. `Trajectory.cpp` comment describing EDR reader provenance

- **Where**: `src/Trajectory.cpp:13`, `src/Trajectory.cpp:339`
- **Claim**: "GROMACS EDR reader (copied layout from former GromacsRunContext)." / "LoadEdr (copied from the dissolved GromacsRunContext::LoadEdr)"
- **Reality**: This is a factually-correct provenance comment (GromacsRunContext is dissolved; the layout was copied from it). Not a wrongness. Flagged as session-narrative adjacent: it relies on a reader knowing what the "former" and "dissolved" references mean.
- **Note**: Kept for completeness; a reviewer may choose to treat this as acceptable provenance documentation rather than rot.

## Unverified

- `spec/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md` G2 claim that MOPAC `skip_mopac` is not honoured at frame 0: requires running the test to verify, beyond static-read scope.
- `spec/INDEX.md:133` "287 tests" — the test count depends on runtime enumeration; test file count verified (49 .cpp, 47 test_*.cpp), individual GTest counts not enumerated.
- `OBJECT_MODEL.md:1512` "Known TrajectoryResult types (as of landing)": table contents (dependencies, lifecycle labels, emission paths) spot-checked against `src/Bs*TrajectoryResult.h` headers for `BsWelfordTrajectoryResult`, `BsShieldingTimeSeries…`, `BondLengthStatsTrajectoryResult`; the rest would require opening each cpp.
- The trajectory-scope section text in `OBJECT_MODEL.md:1377-1611` sometimes refers to "Phase 4 dependency validation" which matches the cpp; whether other phase-number citations are all consistent with an eight-phase labelling was not exhaustively cross-checked.

## Verified clean

- `src/Trajectory.cpp` code body — `Trajectory::Run` is implemented in eight named phases; no `throw` statements anywhere in the file (consistent with `PATTERNS.md` "return status, not exceptions").
- `src/TrajectoryProtein.cpp:69-89` `Seed` — implementation matches docstring at `src/TrajectoryProtein.h:72-94` (calls `FinalizeConstruction`, `AddMDFrame`, `InitTrajectoryAtoms`; does NOT run per-frame calculators).
- `src/ProteinConformation.cpp:22-59` `AttachResult` — implementation matches `PATTERNS.md:91-101` Pattern 4 (singleton check + dep check with diagnostic).
- `OBJECT_MODEL.md:1429-1431` description of `TrajectoryProtein::AttachResult` correctly states dep validation is at Phase 4, not at attach time (OBJECT_MODEL is right; header/cpp docstrings in src are wrong — see Finding 1).
- `OBJECT_MODEL.md:1474-1478` Pattern A labelled "Aspirational; no current fields" — correct against `src/TrajectoryAtom.h` body (see Finding 12 — inconsistency is in the header, not in OBJECT_MODEL).
- `OBJECT_MODEL.md:1550-1567` Run described in 8 named phases — matches `src/Trajectory.cpp:50` and `src/Trajectory.h:14-44`.
- `PATTERNS.md:242-258` §15 "Trajectory::Run is the 8-phase orchestrator" — matches code.
- `PATTERNS.md:197-217` §13 "TrajectoryAtom: private construction, three coexisting shapes" — matches `src/TrajectoryAtom.h` (private ctor + friend + three-shape framing, Pattern A aspirational).
- `src/BsWelfordTrajectoryResult.h` CROSS-RESULT READ (writer side) block — references `BsAnomalousAtomMarkerTrajectoryResult` as a reader; file exists at `src/BsAnomalousAtomMarkerTrajectoryResult.{h,cpp}`.
- `src/RunConfiguration.cpp:111-195` `PerFrameExtractionSet` required ConformationResult set — each `typeid(X)` resolves to a header actually `#include`d in RunConfiguration.cpp (GeometryResult, SpatialIndexResult, EnrichmentResult, DsspResult, ChargeAssignmentResult, ApbsFieldResult, BiotSavartResult, HaighMallionResult, McConnellResult, RingSusceptibilityResult, PiQuadrupoleResult, DispersionResult, HBondResult, SasaResult, EeqResult, AIMNet2Result, WaterFieldResult, HydrationShellResult, HydrationGeometryResult, GromacsEnergyResult, BondedEnergyResult).
- `src/Trajectory.cpp:241-243` comment "Selections are emitted directly … No dynamic_cast sweep here" — consistent with cpp code (no dynamic_cast, no mixin collection pass).
- `src/GromacsFrameHandler.cpp:20-95` — `Open` does NOT read a frame (explicit `has_read_ = false` on line 45); matches `Trajectory.h:77` "No frame is read here".
- `spec/INDEX.md:105` description of `pending_include_trajectory_scope_2026-04-22.md` as not authoritative for landed work — consistent with the renaming narrative; stale content inside that file is expected per its own status.
