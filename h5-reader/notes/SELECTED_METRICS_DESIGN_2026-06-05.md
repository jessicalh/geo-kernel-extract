# Selected Metrics Design Review - 2026-06-05

Scope: read-only design review of selected metric tracking, dashboard panel visibility, data availability, and the existing H5 -> typed model -> views flow. Source scope is `h5-reader/src`. The qt6-cpp model/view and architecture references were read first; the review treats the Windows/macOS/Linux Qt source as first-class portable code.

## Verdict

The internals are steerable, but the current model is not coherent enough to hand off as-is.

The nearest existing source of truth for selected metrics is `DashboardSignalModel::signals_`, exposed as `activeSignals()` in `src/model/DashboardSignalModel.h:44` and stored at `src/model/DashboardSignalModel.h:82`. That should become the single selected-metrics list.

Today, however, rendering and panel membership also depend on `DashboardPanelModel::panels_[].displays`, defined in `src/model/DashboardPanelModel.h:28` and stored at `src/model/DashboardPanelModel.h:102`. The display controller renders only the intersection of active signals and active-panel display refs (`src/app/DashboardDisplayController.cpp:1529`, `src/app/DashboardDisplayController.cpp:1540`, `src/app/DashboardDisplayController.cpp:1575`, `src/app/DashboardDisplayController.cpp:1583`, `src/app/DashboardDisplayController.cpp:2453`). That makes "metric is selected" and "metric is displayed in this active panel" separate, partially synchronized facts.

The rot is not a Qt begin/end model bug. `DashboardSignalModel` and `DashboardPanelModel` mostly obey Qt model mutation contracts. The problem is architectural: a widget-side dialog and an anonymous main-window coordinator are trying to keep two model containers in sync, while dock visibility, selector enabled state, renderer support, and data availability are all derived through separate ad-hoc tests.

## Current Tracking

### Selected signal list

`DashboardSignalModel` is a `QAbstractListModel` with one row per selected dashboard signal (`src/model/DashboardSignalModel.h:11`). Each row stores a `DashboardSignal` containing a `DisplaySignalBinding`, label, `displayModeIds`, enabled flag, axis, anchor, and value shape (`src/model/DashboardSignal.h:261`, `src/model/DashboardSignal.h:271`).

Rows are added by `DashboardSignalModel::addSignal()` with proper `beginInsertRows()` / `endInsertRows()` and `signalAdded()` emission (`src/model/DashboardSignalModel.cpp:233`). Descriptor-based addition creates a `DashboardSignal`, normalizes `displayModeIds`, and chooses the first display mode as the binding's primary mode (`src/model/DashboardSignalModel.cpp:246`, `src/model/DashboardSignalModel.cpp:264`, `src/model/DashboardSignalModel.cpp:267`).

Rows are removed by `removeRows()`, `removeSignal()`, `removeSignalAt()`, and `clear()` (`src/model/DashboardSignalModel.cpp:159`, `src/model/DashboardSignalModel.cpp:273`, `src/model/DashboardSignalModel.cpp:277`, `src/model/DashboardSignalModel.cpp:288`). Mode toggles mutate the same signal row via `setDisplayModes()`, `addDisplayMode()`, `removeDisplayMode()`, and `toggleDisplayMode()` (`src/model/DashboardSignalModel.cpp:361`, `src/model/DashboardSignalModel.cpp:382`, `src/model/DashboardSignalModel.cpp:398`, `src/model/DashboardSignalModel.cpp:420`).

This model is the correct place for the single selected-metrics list, but it is not currently treated as authoritative by itself.

### Panel display refs

`DashboardPanelModel` is a second `QAbstractListModel`, with one row per dashboard panel (`src/model/DashboardPanelModel.h:38`). Each `DashboardPanel` owns a vector of `DashboardDisplayRef` entries (`src/model/DashboardPanelModel.h:14`, `src/model/DashboardPanelModel.h:28`). A display ref is `(signalId, displayModeId, channelId)`, so it is not a metric row; it is a placement/rendering reference.

`DashboardPanelModel` always creates at least one panel through `ensureOnePanel()` (`src/model/DashboardPanelModel.cpp:111`, `src/model/DashboardPanelModel.cpp:450`). If the last panel is removed, `removePanelAt()` does not remove the panel row; it clears that panel's display refs and emits `displayRefRemoved`, `displayRefsChanged`, and `panelRemoved` (`src/model/DashboardPanelModel.cpp:249`). Multiple-panel removal removes a row and emits removed refs for that panel (`src/model/DashboardPanelModel.cpp:263`, `src/model/DashboardPanelModel.cpp:274`).

Display refs are created by `DisplayRefsForSignal()` (`src/model/DashboardPanelModel.cpp:74`). Strip modes become channel refs; implemented static panel modes become a sentinel `channelId="panel"` (`src/model/DashboardPanelModel.cpp:81`, `src/model/DashboardPanelModel.cpp:86`). This function also contains a duplicate panel-mode predicate that must be kept aligned with `DashboardDisplayController.cpp` (`src/model/DashboardPanelModel.cpp:53`).

### Synchronization loop

`ReaderMainWindow` installs an anonymous `DashboardSignalPanelCoordinator` (`src/app/ReaderMainWindow.cpp:99`). It wires signal removal to panel-ref removal (`src/app/ReaderMainWindow.cpp:111`, `src/app/ReaderMainWindow.cpp:124`) and panel-ref removal back to signal removal when the signal has no remaining refs (`src/app/ReaderMainWindow.cpp:115`, `src/app/ReaderMainWindow.cpp:133`, `src/app/ReaderMainWindow.cpp:139`).

This is the only intentional cleanup loop. It explains the current behavior:

- Removing a signal shrinks panel refs.
- Closing a panel shrinks selected signals only if that close removes the last display ref for a signal.
- Removing one display mode can remove the whole signal if that was its last remaining display ref.

The loop is useful, but it is the wrong layer. It is hidden in `ReaderMainWindow.cpp`, not part of either model's contract, so callers cannot reason about selected metrics by reading the selected-metrics model alone.

## Current Events

### Add metric

The dialog add path is two-step:

1. `SignalDisplayDialog::onAddSelected()` gathers checked mode boxes, validates each mode through `TrajectorySignalCatalog::canBind()`, and builds a `displayModes` list (`src/app/SignalDisplayDialog.cpp:1120`, `src/app/SignalDisplayDialog.cpp:1124`, `src/app/SignalDisplayDialog.cpp:1131`).
2. It calls `DashboardSignalModel::addSignal()` (`src/app/SignalDisplayDialog.cpp:1154`) and then creates matching panel refs with `DisplayRefsForSignal()` / `DashboardPanelModel::addDisplayRef()` (`src/app/SignalDisplayDialog.cpp:1161`, `src/app/SignalDisplayDialog.cpp:1166`).

The ordering matters. `signalAdded` fires before panel refs exist, so the display controller rebuilds once with no refs and then relies on `displayRefsChanged` to rebuild again. The source comments call out that exact race (`src/app/DashboardDisplayController.cpp:1229`).

### Remove metric

The active Remove button removes the signal row (`src/app/SignalDisplayDialog.cpp:1313`, `src/app/SignalDisplayDialog.cpp:1321`). The coordinator then removes all panel refs for that signal (`src/app/ReaderMainWindow.cpp:124`, `src/app/ReaderMainWindow.cpp:129`).

### Toggle display mode

Active mode checkboxes first mutate the signal row (`src/app/SignalDisplayDialog.cpp:1290`) and then add/remove panel refs in the active panel (`src/app/SignalDisplayDialog.cpp:1295`, `src/app/SignalDisplayDialog.cpp:1300`, `src/app/SignalDisplayDialog.cpp:1303`, `src/app/SignalDisplayDialog.cpp:1305`).

This can drift if the panel operation fails or has no matching refs. The selected signal row and the placement refs are separate writes, with no transaction-like API.

### Close panel

Closing a panel tab calls `DashboardPanelModel::removePanelAt(row)` (`src/app/DashboardStripDock.cpp:352`, `src/app/DashboardStripDock.cpp:356`). Ref removals emitted from `DashboardPanelModel` then drive selected-signal pruning through the coordinator if the signal refcount reaches zero (`src/app/ReaderMainWindow.cpp:133`, `src/app/ReaderMainWindow.cpp:139`).

This mostly matches the lead model, "closing a panel shrinks selected metrics," but only indirectly and only by reference count. The selected list is not the authority for panel state; panel refs can delete selected signals.

## Why The Selector Looks Greyed Out

There are several distinct grey/disabled states that the UI currently collapses into one user experience.

### Toolbar and dock metric buttons are focus-gated only

The top-level `Metrics...` action starts disabled (`src/app/ReaderMainWindow.cpp:831`) and is enabled only when `AtomSelection` has focus (`src/app/ReaderMainWindow.cpp:328`, `src/app/ReaderMainWindow.cpp:330`). `onOpenSignalDisplays()` also refuses to open without focus (`src/app/ReaderMainWindow.cpp:1192`, `src/app/ReaderMainWindow.cpp:1196`).

The dashboard dock's internal `Metrics...` button is the same: constructed disabled (`src/app/DashboardStripDock.cpp:91`) and enabled only when selection has focus (`src/app/DashboardStripDock.cpp:191`, `src/app/DashboardStripDock.cpp:193`).

None of this derives from selected metrics, data availability, or whether the current selected metric list is non-empty. It is purely an anchor-focus gate.

### Candidate mode boxes are current-row and renderer-gated

After `refreshCatalog()`, the dialog resets the descriptor table, clears candidate selection, and calls `onCandidateSelectionChanged()` (`src/app/SignalDisplayDialog.cpp:958`, `src/app/SignalDisplayDialog.cpp:965`, `src/app/SignalDisplayDialog.cpp:969`). There is no code that selects the first descriptor row after the reset. With no current candidate row, `record` is null and every candidate mode checkbox is disabled (`src/app/SignalDisplayDialog.cpp:1060`, `src/app/SignalDisplayDialog.cpp:1064`, `src/app/SignalDisplayDialog.cpp:1069`, `src/app/SignalDisplayDialog.cpp:1072`).

Even when a descriptor row is current, a checkbox is enabled only if the mode exists and `modeHasVisibleSurface()` returns true (`src/app/SignalDisplayDialog.cpp:1069`, `src/app/SignalDisplayDialog.cpp:1070`). `modeHasVisibleSurface()` currently includes `strip.*`, `static.bar.sequence`, `static.curve.lag.animated`, `static.chord.coupling`, and `static.fixed_freq` (`src/app/SignalDisplayDialog.cpp:218`). It excludes `static.table`, `static.atomColor`, and `static.tensor`, even if the descriptor data exists.

The tooltip says "This display mode does not have an implemented visible renderer" (`src/app/SignalDisplayDialog.cpp:1074`). That is a renderer-capability statement, not a data-availability statement. The UI makes it feel like "nothing is selectable," but the causes can be no current row, no anchor, no implemented renderer, or unavailable data.

### Add button is not list-state-derived

The Add button is enabled only when there is an active model, an anchor candidate, a current descriptor row, and at least one enabled+checked mode checkbox (`src/app/SignalDisplayDialog.cpp:1091`, `src/app/SignalDisplayDialog.cpp:1103`). It does not reflect "the selected metrics list has entries"; it reflects whether the current widget state can create another row.

### Active mode boxes are active-row and catalog-gated

The active-side mode boxes are enabled only if an active signal row is selected, the descriptor is found in the catalog, the mode has a visible surface, and the descriptor's supported modes contain that kind (`src/app/SignalDisplayDialog.cpp:1250`, `src/app/SignalDisplayDialog.cpp:1253`, `src/app/SignalDisplayDialog.cpp:1255`, `src/app/SignalDisplayDialog.cpp:1267`, `src/app/SignalDisplayDialog.cpp:1270`).

There is no auto-selection of the first active signal when the dialog opens with a non-empty selected list. Therefore the active controls can be grey even when `DashboardSignalModel::activeSignals()` is non-empty.

## Panel/Dock Visibility Today

All property docks start hidden, including the dashboard strip dock (`src/app/ReaderMainWindow.cpp:391`, `src/app/ReaderMainWindow.cpp:396`, `src/app/ReaderMainWindow.cpp:398`). The only general visible path is the `QDockWidget::toggleViewAction()` menu/toolbar path (`src/app/ReaderMainWindow.cpp:400`, `src/app/ReaderMainWindow.cpp:415`, `src/app/ReaderMainWindow.cpp:418`, `src/app/ReaderMainWindow.cpp:428`).

Adding a metric does not show the dock. The add path ends after updating models and selecting the active signal in the dialog (`src/app/SignalDisplayDialog.cpp:1179`); it does not call `DashboardStripDock::show()`, `ReaderMainWindow::revealDockQueued()`, or toggle the dock action.

The display controller rebuilds tracks and owned panels (`src/app/DashboardDisplayController.cpp:1504`, `src/app/DashboardDisplayController.cpp:1685`, `src/app/DashboardDisplayController.cpp:1695`) and the dock updates the stack widget/status (`src/app/DashboardStripDock.cpp:281`, `src/app/DashboardStripDock.cpp:297`, `src/app/DashboardStripDock.cpp:303`). Neither layer derives dock visibility from selected metric count.

So the current behavior diverges from the desired rule:

- Desired: dashboard/panel is up by default whenever selected metrics are non-empty, unless the user explicitly disabled it.
- Current: dashboard dock is hidden by default and remains hidden after selected metrics become non-empty unless the user explicitly opens it.

## Data Availability

The availability model is conceptually right and should be the gate.

`TrajectoryFieldAvailability` classifies descriptor data as `Absent`, `NoFramePayload`, `AllMissing`, `AllZeroStructural`, `AllZeroObserved`, or `Available` (`src/model/TrajectoryFieldAvailability.h:24`). `allowsDescriptor()` returns true only for visible states (`src/model/TrajectoryFieldAvailability.h:128`, `src/model/TrajectoryFieldAvailability.h:142`). `AllZeroObserved` is intentionally visible; `AllZeroStructural` is not (`src/model/TrajectoryFieldAvailability.h:142`, `src/model/TrajectoryFieldAvailability.h:145`, `src/model/TrajectoryFieldAvailability.h:150`).

`ReaderMainWindow` builds availability after loading, then installs it on the `TrajectorySignalCatalog` and the atom inspector (`src/app/ReaderMainWindow.cpp:290`, `src/app/ReaderMainWindow.cpp:291`, `src/app/ReaderMainWindow.cpp:294`, `src/app/ReaderMainWindow.cpp:295`).

`TrajectorySignalCatalog` filters its visible descriptor list through `availability_->allowsDescriptor()` (`src/model/TrajectorySignalCatalog.cpp:1176`, `src/model/TrajectorySignalCatalog.cpp:1182`, `src/model/TrajectorySignalCatalog.cpp:1186`). The dialog refreshes from `catalog->descriptorList()` (`src/app/SignalDisplayDialog.cpp:958`), so the picker table is indirectly availability-gated. `canSample()` also checks `availability_->canSampleDescriptor()` for temporal strip sampling (`src/model/TrajectorySignalCatalog.cpp:1282`, `src/model/TrajectorySignalCatalog.cpp:1286`).

The gap: `DashboardSignalModel` and `DashboardPanelModel` do not know availability. `DashboardSignalModel::addSignal()` accepts any descriptor object it is handed (`src/model/DashboardSignalModel.cpp:246`) and no model role records availability. The display controller looks up descriptors in the filtered catalog (`src/app/DashboardDisplayController.cpp:1562`); if a selected signal's descriptor is not in that filtered catalog, it silently does not render (`src/app/DashboardDisplayController.cpp:1563`, `src/app/DashboardDisplayController.cpp:1564`).

Availability and renderer capability are currently conflated in the UI. Availability decides whether a descriptor should be offered. `modeHasVisibleSurface()` decides whether a display mode has an implemented renderer (`src/app/SignalDisplayDialog.cpp:218`). Those are orthogonal:

- Available data + no renderer: descriptor can exist, but a particular display mode cannot produce a visible surface yet.
- Unavailable data + implemented renderer: descriptor should not be offered or should be shown disabled with an availability reason.
- Available data + implemented renderer: selectable and renderable.

The current grey checkboxes explain only the renderer axis, not the data axis, and the user cannot see which axis blocked the action.

## Existing Model-Is-Spine Flow

The code already has a strong spine and the selected-metrics design should fit it.

Load starts at `QtProteinLoader::LoadRunPath()` (`src/io/QtProteinLoader.h:82`), called by `main_reader.cpp` (`src/main_reader.cpp:191`, `src/main_reader.cpp:193`). The loader returns a typed `QtLoadResult` containing `QtProtein` and a `Conformation` subclass (`src/io/QtProteinLoader.h:29`, `src/io/QtProteinLoader.h:30`, `src/io/QtProteinLoader.h:31`). For trajectory runs, `LoadTrajectory()` reads topology sidecar data, constructs `QtTrajectoryH5`, builds `QtProtein`, then builds `TrajectoryConformation` (`src/io/QtProteinLoader.cpp:85`, `src/io/QtProteinLoader.cpp:92`, `src/io/QtProteinLoader.cpp:121`, `src/io/QtProteinLoader.cpp:134`).

`TrajectoryConformation` owns the dense H5 layer and exposes the common `Conformation` seam (`src/model/TrajectoryConformation.h:1`, `src/model/TrajectoryConformation.h:54`, `src/model/TrajectoryConformation.h:55`). `Conformation` owns the frame/snapshot seam and emits `snapshotReady()` when frame-local data becomes resident (`src/model/Conformation.h:82`, `src/model/Conformation.h:89`, `src/model/Conformation.cpp:31`, `src/model/Conformation.cpp:56`).

`ReaderMainWindow` then wires typed model objects into views: scene, picker, inspector, selection model, signal catalog, dashboard signal model, panel model, dialog, and dashboard dock (`src/app/ReaderMainWindow.cpp:198`, `src/app/ReaderMainWindow.cpp:270`, `src/app/ReaderMainWindow.cpp:288`, `src/app/ReaderMainWindow.cpp:290`, `src/app/ReaderMainWindow.cpp:296`, `src/app/ReaderMainWindow.cpp:297`, `src/app/ReaderMainWindow.cpp:298`, `src/app/ReaderMainWindow.cpp:379`).

The dashboard signal -> panel -> strip flow is therefore a branch of the same spine:

H5 / frame snapshots -> typed descriptors and samplers -> `TrajectorySignalCatalog` / `DashboardSignalModel` -> `DashboardDisplayController` -> `DashboardStripDock` / `StripStackWidget`.

The selected-metrics list belongs in that model branch. It should not be bolted onto `SignalDisplayDialog`, `DashboardStripDock`, or `StripStackWidget`.

## Proposed Clean Model

### 1. Promote `DashboardSignalModel` to the single selected-metrics source of truth

Keep one row per selected metric in `DashboardSignalModel`. Treat `activeSignals()` as selected metrics, not just "active signals."

Add explicit model roles or accessors for:

- `SelectedMetricId` / existing UUID
- descriptor id, concept key, source kind, anchor, follows-focus
- selected display modes
- enabled flag
- availability state and reason
- renderability state per display mode

The model should be able to answer:

- selected metric count
- renderable selected display count
- selected-but-unavailable count
- selected-but-no-renderer count

Qt model/view implication: these are model facts and roles. Dialogs, docks, toolbar actions, and status labels should subscribe to `rowsInserted`, `rowsRemoved`, `modelReset`, and `dataChanged`, not maintain their own count.

### 2. Demote `DashboardPanelModel` to placement/layout

`DashboardPanelModel` should not be a co-owner of selected metrics. It should store layout: panel tabs, active panel, and per-panel placement refs for selected metric display modes.

Allowed model relationship:

- `DashboardSignalModel` says what is selected.
- `DashboardPanelModel` says where selected display refs are placed.
- A single model-layer coordinator/controller owns compound mutations so the two models do not drift.

The current anonymous `DashboardSignalPanelCoordinator` should become a named model-layer class, for example `DashboardSelectionController`, with methods such as:

- `addMetric(descriptor, anchor, modes, targetPanel)`
- `removeMetric(signalId)`
- `setMetricMode(signalId, mode, enabled, targetPanel)`
- `removePanel(panelId, PanelRemovalPolicy::RemoveReferencedMetrics)`
- `clearAllMetrics()`

The dialog and dock should call these methods, rather than writing `DashboardSignalModel` and `DashboardPanelModel` separately.

### 3. Define selected-list events

Events that grow the selected list:

- User adds an available descriptor with at least one selected display mode.
- Programmatic/test API adds the same through the same controller.

Events that shrink the selected list:

- User removes an active metric.
- User disables the last selected display mode for a metric.
- User closes a panel and the design policy is "closing a panel removes its referenced metrics."
- Availability changes to unavailable, if the product decision is prune; otherwise mark unavailable and stop rendering.

Do not let a widget-side refcount callback be the only place that enforces this. The selected-list model should emit normal Qt model changes for every grow/shrink event.

### 4. Derive panel/dock visibility from selected metrics plus user intent

Use one application-level visibility rule:

`dashboardDockVisible = selectedMetricCount > 0 && !dashboardDockUserDisabled`

When `selectedMetricCount` changes from 0 to non-zero, call the existing `revealDockQueued(dashboardStripDock_)` unless the user explicitly hid/disabled the dock. When the user toggles the dock closed, set `dashboardDockUserDisabled = true`. When the user explicitly opens it, clear that flag. This is portable Qt state, not platform-specific window behavior.

The dock's contents should derive from renderable selected displays:

- If selected metrics exist and renderable displays exist, show panels/strips.
- If selected metrics exist but no selected display has a renderer, keep the selected list visible in the selector and show a dashboard status explaining "selected metrics have no implemented dashboard surface."
- If selected metrics exist but data later becomes unavailable, show the selected list with unavailable state and no renderer output.

Do not use dock visibility as the source of truth for whether metrics exist.

### 5. Separate data availability from renderer capability

`TrajectoryFieldAvailability` gates what data exists. It belongs at descriptor-selection time:

- The picker candidate model should offer only descriptors allowed by `TrajectoryFieldAvailability::allowsDescriptor()`, or show unavailable descriptors disabled with a reason if discoverability is desired.
- `DashboardSignalModel` should validate adds against catalog/availability so no code path can add an unavailable descriptor silently.
- Existing selected rows should carry availability state so views can show "unavailable" rather than disappearing when `catalog->findDescriptor()` returns null.

Renderer capability is a different model axis:

- Replace local `modeHasVisibleSurface()` with shared display-mode capability metadata, probably in `DashboardDisplayController`/catalog-adjacent model code, not buried in `SignalDisplayDialog`.
- Replace duplicate panel-mode predicates in `DashboardPanelModel.cpp` and `DashboardDisplayController.cpp` with one shared function or table (`src/model/DashboardPanelModel.cpp:53`, `src/app/DashboardDisplayController.cpp:59`).
- Display modes with no renderer should be marked "not renderable yet," not treated as unavailable data.

### 6. Selector should reflect reality

The dialog should not present an all-grey control surface when model state exists.

Concrete fixes:

- After `refreshCatalog()` and after anchor/filter changes, if the candidate proxy has rows and no current candidate row, select row 0.
- When the active signal model has rows and no current active row, select row 0.
- Enable/disable candidate mode boxes based on three explicit facts: descriptor available, mode supported by descriptor, mode renderable. Show different reasons for each.
- The Add button should explain which required fact is missing: anchor, descriptor row, renderable mode, or availability.
- The active-side controls should reflect selected metrics from `DashboardSignalModel` even if the selected metric currently has no renderable display.

## Internals Assessment

`DashboardSignalModel`: structurally acceptable Qt model, but too permissive. It lacks catalog/availability validation, and mode mutation emits only row changes; it does not guarantee corresponding panel placement changes.

`DashboardPanelModel`: structurally acceptable layout model, but it is currently treated as a second owner of metric existence. `DisplayRefsForSignal()` is useful, but the duplicate `isPanelDisplayMode()` predicate is a maintenance hazard.

`DashboardDisplayController`: good place for rendering derivation, not for selection truth. It correctly rebuilds from model signals, but it renders only the active-panel intersection, so its status cannot be the selected-list truth. It also silently drops selected signals whose descriptors are absent from the filtered catalog.

`DashboardStripDock`: view/dock layer only. It should not decide what selected metrics exist. It currently gates the metric picker button by focus only and never derives dock visibility from selected metric count.

`StripStackWidget`: coherent renderer surface. It already documents that data ownership is outside the widget (`src/app/StripStackWidget.h:1`). Keep it that way.

`SignalDisplayDialog`: currently does too much coordination. It writes both selected signal rows and panel refs, owns local renderer-capability rules, and bases enabled state on current widget selection rather than model state. It should become a view/editor over the selected-metrics model and a caller of a model-layer selection controller.

## Fix Plan

1. Establish a named model-layer coordinator/controller for selected metric mutations. Move the anonymous cleanup loop out of `ReaderMainWindow.cpp`.

2. Make `DashboardSignalModel` the selected-metrics source of truth. Add count/availability/renderability roles or helper queries. Ensure every grow/shrink event mutates this model first and emits Qt model changes.

3. Turn `DashboardPanelModel` into placement state. Keep panel tabs and refs, but stop treating refcount callbacks as the hidden authority for selected metric existence.

4. Add a shared display-mode capability table. Use it from the dialog, panel ref generation, and display controller. Remove duplicated panel-mode lists.

5. Push `TrajectoryFieldAvailability` through the whole selection path. Catalog filtering is a good start, but model adds and existing selected rows also need availability awareness.

6. Implement dock visibility derivation in `ReaderMainWindow` or a small UI state controller: selected non-empty reveals dashboard unless user disabled it; empty selection can hide or show an empty state depending on product choice.

7. Fix selector reflection: auto-select first candidate/active row when appropriate, and split disabled reasons into "no anchor," "data unavailable," "mode unsupported," and "mode has no renderer."

This is not a full dashboard rebuild. The existing typed data flow and the two Qt models can be kept. But the selected metric list must become a real model contract, and the current two-model/widget-coordinator behavior should be steered into a single selected-metrics source of truth plus derived layout/rendering views.
