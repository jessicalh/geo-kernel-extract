# Codex — STEP 2.3: availability + renderability as first-class roles on DashboardSignalModel

## FIRST — read the qt6-cpp skill before touching any code
Read, in this order, and treat them as the law for this work (you have full filesystem access — open them):
- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/model-view.md`  ← custom models, roles, `dataChanged`, delegates
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md` ← signal/slot, ownership, startup/shutdown

This is `QAbstractListModel` role work — the model-view reference is directly on point. Honor it:
`CENSUS_REGISTER` in any new QObject ctor, `ACONNECT` (not raw `connect`), `ASSERT_THREAD(this)` on
thread-sensitive methods, **no new `QTimer`**, full error handling at any external boundary, portable Qt
(Windows/macOS/Linux first-class).

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`, on commit
`8007e2e` + the uncommitted 2.2 capability table). Scope: `h5-reader/src` ONLY. Never link the
`nmr_shielding` library, never write H5. **Verify from the code before changing anything, ever.**

## Why this matters
We're making `DashboardSignalModel` the single source of truth for selected metrics (the next step gives the
mutations a named controller). Before that, the model must actually *carry* the two facts views keep
re-deriving ad hoc: **is this metric's data available**, and **is each of its display modes renderable**.
Today neither is model state — availability is recomputed in scattered places and a selected metric whose
descriptor goes missing is *silently dropped* mid-rebuild. You're making the model hold the truth so the
controller and every view read one place.

## Read first (confirm before building)
- `DashboardSignalModel` — the rows + existing roles (`src/model/DashboardSignalModel.h`,
  `src/model/DashboardSignalModel.cpp`). One row per selected metric (`DashboardSignal`).
- `TrajectoryFieldAvailability` (`src/model/TrajectoryFieldAvailability.h`) — the loaded-run classifier
  (`Absent / NoFramePayload / AllMissing / AllZeroStructural / AllZeroObserved / Available`) +
  `allowsDescriptor()` / `canSampleDescriptor()`. How `TrajectorySignalCatalog` already consults it
  (`fieldAvailability()`).
- The 2.2 capability table `model::DisplayModeCapabilityFor(mode)` (`src/model/DisplayModeCapability.h`) —
  the single answer for `hasVisibleSurface` / `buildsPanelWidget` / `emitsPanelRef`.
- The silent-drop site: `DashboardDisplayController` rebuild, where a selected signal whose descriptor isn't
  found in the filtered catalog is skipped with `if (!descriptor) continue;` — find it and confirm.

## Deliverable (mostly additive model state)
1. **Availability on each selected row.** Add a role/accessor on `DashboardSignalModel` for the metric's
   availability state + reason, resolved from `TrajectoryFieldAvailability` for that descriptor. The model
   needs a handle to availability (inject the same `shared_ptr<TrajectoryFieldAvailability>` the catalog
   gets, via a setter wired in `ReaderMainWindow` where the catalog/inspector are wired — read that wiring).
   Recompute the row's availability on a `refreshAvailability()` call (e.g. after a run loads). Emit
   `dataChanged` for the availability role when it changes.
2. **Per-mode renderability** comes from `DisplayModeCapabilityFor` — expose, per row, which of its
   `displayModeIds` are renderable (have a surface or build a panel) vs ref-only/none. A role or a small
   accessor (`renderableModeCount(row)` / `modeRenderability(row)`); don't duplicate the table.
3. **No silent vanish.** A selected metric whose descriptor can't be resolved must carry an "unavailable"
   state and stay in the model (so views can show it greyed/with a reason), instead of being dropped. Fix the
   `if (!descriptor) continue;` drop to mark-unavailable-and-skip-rendering, not erase. (This is the one
   intended behavior change: unavailable shows as unavailable rather than disappearing.)
4. **Counts the model answers** (cheap accessors, derived from rows): `selectedCount()`,
   `renderableSelectedCount()`, `unavailableCount()`, `noRendererCount()`.
5. **Single source for the instrument.** Repoint `/dashboard/state`'s per-metric `availability` and per-mode
   flags at these new model accessors (STEP 1 computed availability in RestServer ad hoc; now it reads the
   model). The JSON shape stays the same.

Keep it additive: existing add/remove/mode/toggle behavior is unchanged except the no-silent-vanish fix.
Do NOT introduce the named controller yet (that's the next step) and do NOT change the dialog's gating logic
(2.1 owns that). Availability *recording* here; availability *policy* (refuse-add vs mark) is the controller's
call later — for now record + expose, don't hard-refuse adds.

## Behaviour check (state it in notes)
Additive roles + counts shouldn't change any existing `/dashboard/state` value for a healthy run; the only
new observable is that a metric with a missing descriptor reports `availability != Available` and stays
listed instead of vanishing. I will verify by driving the instrument before/after.

## Output
Implement, build green (`cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target
h5reader -j$(nproc)`). **Do NOT launch, do NOT git.** Append a short section to
`notes/INSTRUMENT_DASHBOARD_2026-06-05.md` (new roles/accessors + the no-vanish fix), and note anything in the
code that contradicted this brief (code wins).
