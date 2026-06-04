# Reader robustness backlog (2026-05-30 audit)

Status trued 2026-06-04 against `UI_STATE_OVERVIEW_2026-06-04.md`.

| Item | Current status |
|---|---|
| Reader-local H5 exception boundaries / try-catch | done for the current reader-local pattern; verify legacy readers before adding new work |
| Required `/trajectory/positions` validation | done; malformed required positions hard-error |
| Slider scrub reflection / pile-up handling | done |
| QSettings persistence | done |
| Optional TR polish | verify before acting |
| Windows crash capture and residual render probe | open |

Settled fix list from the Codex + Claude general-purpose audits run
2026-05-30. Both agents read the qt6-cpp skill and audited the reader
against six robustness dimensions with an "appropriateness lens" (no
dogmatic decoupling; only flag what will actively hurt in 6 months).

The audit conversation and per-dimension findings are in the prior
session transcript. This doc is the actionable output.

## What the audits confirmed is fine as-is

- **Reload / hot-swap**: `ReaderMainWindow::onOpenDirectory`
  (`src/app/ReaderMainWindow.cpp:674`) launches a fresh process via
  `QProcess::startDetached`. Deliberate choice that sidesteps in-place
  rebuild entirely. Do not refactor to support in-window reload unless
  someone explicitly asks for it.
- **Architectural complexity**: `ReaderMainWindow`, `Conformation`,
  `DashboardDisplayController` all earn their keep. ReaderMainWindow
  is wiring, not god class. Conformation is deliberate facade
  (`Conformation.h:11-22` documents why it stays thin). DashboardDisplayController
  is large but follows clean per-TR builder dispatch.
- **Connection lifecycle**: dialogs and docks are construct-once +
  show/hide. Qt's parent-tree auto-disconnect handles cleanup. No
  ACONNECT-without-disconnect leaks under current usage.
- **Long-running stability**: memory profile is stable. REST server
  lifecycle clean. Snapshot residency bounded. 24h idle uptime fine.
- **Settings persistence absence**: stale. `QSettings` persistence is
  present as of the June 4 UI survey; keep only as historical rationale.

## Fixes to land

### 1. Wrap unwrapped `Read*` functions in `QtTrajectoryH5.cpp` — status: verify/done-current-pattern

**What**: ~25 `Read*TimeSeries` / `ReadSpecial*` functions currently lack
`try/catch` blocks. Only the 5 newest readers
(`ReadKernelCoherence`, `ReadDihedralAutocorrelation`,
`ReadReorientationalDynamics`, `ReadKernelDynamics`,
`ReadIRedOrderParameters`) have the proper pattern.

**Why**: A HighFive exception inside any unwrapped reader propagates
up through the constructor and fails the WHOLE load — not just that
group. Sparse-tolerance is theatre for these: works for absent groups,
breaks for truncated/corrupt ones. Real failure mode on flaky network
mounts or 40+ GB H5 files.

**How**: Mechanical change. The infrastructure exists:
- `LogReadException(group_path, kind, e)` helper
- noexcept handler pattern (see the 5 wrapped readers as template)
Wrap each `Read*` body in:
```cpp
try {
    // existing read code
} catch (const HighFive::Exception& e) {
    LogReadException(...);
    out.reset();
} catch (const std::exception& e) {
    // ...
}
```

**Where**: `src/io/QtTrajectoryH5.cpp` (25 functions).

**Effort**: ~half day mechanical.

---

### 2. Make missing/malformed `/trajectory/positions` a hard load error — status: done

**What**: Currently `/trajectory/positions` is documented as
always-present but `TrajectoryConformation::atomPosition`
(`src/model/TrajectoryConformation.cpp:54`) returns zero when
positions are absent. `MoleculeScene::Build` and `setFrame` then
render those zeros as if they were real coordinates — atoms collapsed
at origin, silently wrong data.

**Why**: This is the worst error-handling mode: silent wrong-data
display, not visible failure. A partial or truncated H5 opens
successfully and shows a physically false structure. Users have no
indication anything is wrong.

**How**: At `QtTrajectoryH5` construction, validate:
- positions buffer exists (non-null)
- positions buffer frame count matches the trajectory's declared
  frame count
- positions buffer atom count matches the topology's atom count

On any mismatch: throw with a clear message that propagates to the
loader as a hard load failure (rc=3 from main, `ErrorBus` warning).

**Where**: `src/io/QtTrajectoryH5.cpp` constructor validation block,
likely after the positions read.

**Effort**: ~half day with tests.

---

### 3. Fix slider-drag pile-up in `DashboardDisplayController::extendToFrame` — status: done

**What**: `DashboardDisplayController::extendToFrame`
(`src/app/DashboardDisplayController.cpp:2449`) currently iterates
from `startFrame` to `frame`, calling `Conformation::requestSnapshot`
for each gap frame for any series with `needsFrameSnapshot=true`.
Dragging the slider end-to-end can stack up 750 synchronous NPY
directory reads on the GUI thread inside one `valueChanged` cascade.

**Why**: Multi-second GUI freeze on rapid slider drag, which is normal
interactive use. Letting intermediates be `Pending` gaps is a UX
compromise (placeholder flicker) not a real fix.

**How**: Use the proper Qt slider-drag pattern. Distinguish between:

- **Human scrub during drag** (`QSlider::isSliderDown() == true`):
  only update cheap things — playhead position, time cursor on strip,
  molecular scene atom positions (in-memory, ~5ms). NO snapshot
  fetches, NO `extendToFrame` walking.
- **Slider release** (`QSlider::sliderReleased` signal): run
  `extendToFrame` once for the final frame; fetch snapshot async if
  needed.
- **Programmatic playback** (timer-driven, `isSliderDown() == false`):
  run as today; the per-tick gap is always 1 frame so iteration is
  cheap.

`QtPlaybackController` needs to expose `isSliderDown()` state (or the
slider's signal directly); `DashboardDisplayController::extendToFrame`
bails early if `playback.isSliderDown()` is true.

**Where**: `src/app/QtPlaybackController.{h,cpp}` + branch in
`src/app/DashboardDisplayController.cpp::extendToFrame`.

**Effort**: ~30-50 lines of change; ~half day with manual verification.

---

### 4. Normalize optional TR error handling — status: verify

**What**: Some optional TR readers catch and warn internally (the 5
new ones); others rely on exceptions propagating to the top-level
loader. Malformed optional TR data degrades silently in some cases
and fails the whole trajectory load in others. Inconsistent.

**Why**: Surprising and inconsistent behavior. Either an optional TR
that fails to read should always degrade to null (the user sees
"this signal not available"), or it should always fail load — but
not both depending on which reader hits the error.

**How**: Adopting #1 (wrap all readers) accomplishes this. After #1,
every `Read*` failure becomes "this TR group is null, log a warning."
The loader's top-level try/catch becomes the safety net for truly
unrecoverable cases (out-of-memory, file disappeared mid-read).

**Where**: Falls out of fix #1.

**Effort**: zero additional (subsumed by fix #1).

---

### 5. Wire Windows `CrashHandler` — status: open

**What**: `CrashHandler.cpp:170-173` Windows path is a stub. Linux
minidump works; Windows produces nothing on crash.

**Why**: Already tracked in `POLISH_BACKLOG.md` item 9. First Windows
crash without diagnostics is a debugging dead-end. Critical before
any adviser distribution targets Windows.

**How**: Implement Windows minidump via `MiniDumpWriteDump` API
(SEH-handled, writes a .dmp file with companion .txt for object
census + connection audit + log tail). The qt6-cpp skill at
`/home/jessica/.claude/skills/qt6-cpp/references/crash-diagnosis.md`
documents the pattern; templates at
`/home/jessica/.claude/skills/qt6-cpp/templates/crash-handler.h` and
`.cpp` are the starting point.

**Where**: `src/diagnostics/CrashHandler.{h,cpp}` Windows branch.

**Effort**: ~1 day including testing on a Windows build.

---

### 6. Resolve the `MoleculeScene::setFrame` rendering diagnostic probe — status: open

**What**: `MoleculeScene.cpp` around line 501-527 has a diagnostic
probe (`mapper_->Modified()` after the atom-position loop) that was
added to investigate intermittent end-of-trajectory atom-render drop.
The comment frames it as a probe, not a settled fix. `feedback_vtk_bounds_cache`
memory entry suggests adjacent VTK pipeline confusion.

**Why**: "For now" code that risks becoming permanent while the visual
bug remains ambiguous. Each frame change pays for an extra
`Modified()` call whose causal role hasn't been confirmed. If it's
actually fixing the bug, document it as the settled answer. If it
isn't, remove it.

**How**:
- Reproduce the original render-drop bug on a known trajectory with
  the probe disabled
- If still happens: the probe is doing something real; document it
  as the fix with a comment explaining why `Modified()` is needed
  here vs the normal VTK invalidation chain
- If no longer happens: the probe is dead code; remove it

**Where**: `src/app/MoleculeScene.cpp` lines 501-527 (probe area) +
related `feedback_vtk_bounds_cache` workaround at lines 434-445.

**Effort**: ~half day investigation + decision.

---

### Judgment call: minimal versioned settings persistence — status: done

**What**: No app settings persistence today. Window layout, recent
files, dashboard/panel configurations reset on every restart.

**The disagreement**: Codex recommended adding it ("medium user-pain
risk over six months"); Claude said fine as-is ("do not pre-emptively
add"). User decides.

**If we do it**:
- Use `QMainWindow::saveState(version)` + `restoreState(blob, version)`
- Version-tag the blob; discard on version mismatch (schema-evolution
  safe pattern)
- Persist via `QSettings` for window geometry + recent files; custom
  serialization for dashboard signal/panel state via `QSaveFile` if
  not using `QSettings`
- Tolerant restore semantics: unknown/missing docks ignored, not
  errors

**Where**: New `src/app/AppSettings.{h,cpp}` plus
`ReaderMainWindow::saveState`/`restoreState` overrides.

**Effort**: ~1-2 days including the panel-config serialization.

**Defer until** the user feels the friction of rebuilding the layout
each restart.

---

## Suggested work order

These are independent; can land in any order, by any session.
Suggested dependency-aware sequence:

1. **Fix #2 (positions hard requirement)** — important for data
   correctness; cheap; ship first.
2. **Fix #3 (slider drag pile-up)** — important for interactive use;
   cheap; user-visible quality improvement.
3. **Fix #1 (wrap Read* functions)** — important for partial-data
   tolerance; mechanical; can be done in parallel with #2/#3 by a
   different session. Fix #4 falls out automatically.
4. **Fix #5 (Windows CrashHandler)** — before Windows distribution
   pass begins. Self-contained; can be done in any session that has
   a Windows build environment.
5. **Fix #6 (resolve rendering probe)** — important for thesis demos
   where rendering correctness shows; investigation rather than
   coding; do when there's bandwidth for the investigation.
6. **Judgment call (settings persistence)** — only if/when user
   feels the friction.

## What is NOT in this backlog (do not re-litigate)

- Multi-protein simultaneous support (user explicitly removed from scope)
- Architectural decoupling of ReaderMainWindow / Conformation /
  DashboardDisplayController (both audits agreed: appropriate as-is)
- Strip UI internals (off-limits)
- Stage 2 chewer additions (covered separately in
  `analysis-speculative/reader_as_platform_2026-05-29.md` +
  `analysis-speculative/stage2_pysr_campaign_2026-05-29.md`)

## Provenance

- Audit prompt: see prior session transcript
- Codex output: `/tmp/codex-reader-robustness-output.txt` (transient)
- Claude agent output: `/tmp/claude-1000/-shared-2026Thesis-nmr-shielding/.../tasks/` (transient)
- Audit conversation: prior session 2026-05-30, with explicit
  "appropriateness lens" framing rejecting dogmatic patterns
- Both agents read `/home/jessica/.claude/skills/qt6-cpp/SKILL.md` as
  required pre-reading
