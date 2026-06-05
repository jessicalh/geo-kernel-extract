# UI Interaction Paths Review - 2026-06-05

Read-only review of `h5-reader/src` interaction paths. No build, launch, git, or source edits were performed.

## 1. Metric -> Panel -> Visible Dock: Root Cause And Fix

### Finding

The live bug is real and code-level: adding a metric can successfully create the signal, panel display refs, strip tracks, and/or owned panels, but nothing shows the `DashboardStripDock`. The dock is explicitly hidden at startup, and neither opening the Metrics dialog nor committing Add calls `show()`, `setVisible(true)`, `raise()`, `toggleViewAction()`, or any equivalent dock reveal path.

The user sees nothing because the rendered result is placed into a hidden dock.

### Flow Map

1. Startup hides all property docks.
   - `kSettingsVersion` is bumped for the hidden-docks layout policy at `src/app/ReaderMainWindow.cpp:88` and `src/app/ReaderMainWindow.cpp:92`.
   - Inspector, Selection, and Strip are tabified together, then all three are hidden: `src/app/ReaderMainWindow.cpp:372`, `src/app/ReaderMainWindow.cpp:387`, `src/app/ReaderMainWindow.cpp:388`, `src/app/ReaderMainWindow.cpp:396`, `src/app/ReaderMainWindow.cpp:397`, `src/app/ReaderMainWindow.cpp:398`.
   - `restoreAllSettings()` runs after dock creation; a missing/version-mismatched state leaves those hidden defaults intact: `src/app/ReaderMainWindow.cpp:511`, `src/app/ReaderMainWindow.cpp:515`.

2. The Metrics action only opens the dialog.
   - Toolbar action is disabled until there is selection focus: `src/app/ReaderMainWindow.cpp:820`, `src/app/ReaderMainWindow.cpp:821`, `src/app/ReaderMainWindow.cpp:328`, `src/app/ReaderMainWindow.cpp:340`.
   - Opening the dialog refreshes catalog and shows/raises only the dialog: `src/app/ReaderMainWindow.cpp:1173`, `src/app/ReaderMainWindow.cpp:1181`, `src/app/ReaderMainWindow.cpp:1182`, `src/app/ReaderMainWindow.cpp:1183`, `src/app/ReaderMainWindow.cpp:1184`.
   - The dock's own Metrics button is also only a request to open that dialog: `src/app/DashboardStripDock.cpp:126`, `src/app/DashboardStripDock.cpp:127`, wired to the same slot at `src/app/ReaderMainWindow.cpp:446`, `src/app/ReaderMainWindow.cpp:447`.

3. The dialog is intentionally renderer-agnostic.
   - Its header says it only updates dashboard models and does not create strip widgets/tables/overlays: `src/app/SignalDisplayDialog.h:1`, `src/app/SignalDisplayDialog.h:3`, `src/app/SignalDisplayDialog.h:5`.

4. Add creates model state.
   - Add requires an anchor, selected descriptor, and at least one checked display mode: `src/app/SignalDisplayDialog.cpp:1083`, `src/app/SignalDisplayDialog.cpp:1095`.
   - The default mode selection tries Strip first, then the first supported visual mode: `src/app/SignalDisplayDialog.cpp:294`, `src/app/SignalDisplayDialog.cpp:302`, `src/app/SignalDisplayDialog.cpp:310`, `src/app/SignalDisplayDialog.cpp:1058`, `src/app/SignalDisplayDialog.cpp:1065`, `src/app/SignalDisplayDialog.cpp:1071`, `src/app/SignalDisplayDialog.cpp:1076`.
   - `onAddSelected()` calls `DashboardSignalModel::addSignal(...)`: `src/app/SignalDisplayDialog.cpp:1146`, `src/app/SignalDisplayDialog.cpp:1147`, `src/app/SignalDisplayDialog.cpp:1152`.
   - `DashboardSignalModel::addSignal()` normalizes, inserts the row, emits `signalAdded`, and returns the UUID: `src/model/DashboardSignalModel.cpp:233`, `src/model/DashboardSignalModel.cpp:240`, `src/model/DashboardSignalModel.cpp:242`, `src/model/DashboardSignalModel.cpp:243`.

5. Add creates panel refs for supported strip/panel modes.
   - `DisplayRefsForSignal()` creates refs for static panel modes and `strip.*` modes: `src/model/DashboardPanelModel.cpp:74`, `src/model/DashboardPanelModel.cpp:86`, `src/model/DashboardPanelModel.cpp:87`, `src/model/DashboardPanelModel.cpp:90`, `src/model/DashboardPanelModel.cpp:94`, `src/model/DashboardPanelModel.cpp:99`, `src/model/DashboardPanelModel.cpp:101`.
   - The dialog adds those refs and activates the target panel: `src/app/SignalDisplayDialog.cpp:1155`, `src/app/SignalDisplayDialog.cpp:1158`, `src/app/SignalDisplayDialog.cpp:1161`.
   - `DashboardPanelModel::addDisplayRef()` emits `displayRefAdded` and `displayRefsChanged`: `src/model/DashboardPanelModel.cpp:332`, `src/model/DashboardPanelModel.cpp:338`, `src/model/DashboardPanelModel.cpp:340`, `src/model/DashboardPanelModel.cpp:341`.

6. Controller rebuild now happens correctly after refs exist.
   - The controller still rebuilds on `signalAdded`, but also rebuilds on `displayRefsChanged` specifically because add is a two-step signal/ref operation: `src/app/DashboardDisplayController.cpp:1209`, `src/app/DashboardDisplayController.cpp:1210`, `src/app/DashboardDisplayController.cpp:1229`, `src/app/DashboardDisplayController.cpp:1238`, `src/app/DashboardDisplayController.cpp:1239`.
   - `rebuild()` filters static panels by active panel refs and builds owned panel widgets for implemented modes: `src/app/DashboardDisplayController.cpp:1504`, `src/app/DashboardDisplayController.cpp:1575`, `src/app/DashboardDisplayController.cpp:1583`, `src/app/DashboardDisplayController.cpp:1588`, `src/app/DashboardDisplayController.cpp:1592`, `src/app/DashboardDisplayController.cpp:1640`, `src/app/DashboardDisplayController.cpp:1685`, `src/app/DashboardDisplayController.cpp:1696`, `src/app/DashboardDisplayController.cpp:1697`.
   - Temporal strip modes build `ChannelBuffer` tracks: `src/app/DashboardDisplayController.cpp:1665`, `src/app/DashboardDisplayController.cpp:1667`, `src/app/DashboardDisplayController.cpp:1668`, `src/app/DashboardDisplayController.cpp:2302`, `src/app/DashboardDisplayController.cpp:2310`, `src/app/DashboardDisplayController.cpp:2345`.

7. The dock receives rendered content, but only inside the hidden dock.
   - `ownedPanelsChanged` moves owned panel widgets into `StripStackWidget`: `src/app/DashboardStripDock.cpp:112`, `src/app/DashboardStripDock.cpp:119`, `src/app/DashboardStripDock.cpp:122`, `src/app/StripStackWidget.cpp:441`, `src/app/StripStackWidget.cpp:452`, `src/app/StripStackWidget.cpp:455`.
   - `stripTracksChanged` pushes temporal tracks into `StripStackWidget`: `src/app/DashboardStripDock.cpp:280`, `src/app/DashboardStripDock.cpp:286`, `src/app/DashboardStripDock.cpp:296`, `src/app/DashboardStripDock.cpp:301`.
   - `StripStackWidget::paintEvent()` paints panels only when the widget is visible in the dock hierarchy: `src/app/StripStackWidget.cpp:562`, `src/app/StripStackWidget.cpp:568`, `src/app/StripStackWidget.cpp:589`, `src/app/StripStackWidget.cpp:592`.

8. No reveal exists on Add.
   - `onAddSelected()` ends after setting status and selecting the active row; there is no dock reveal call: `src/app/SignalDisplayDialog.cpp:1163`, `src/app/SignalDisplayDialog.cpp:1171`, `src/app/SignalDisplayDialog.cpp:1179`.
   - Search confirms the only explicit dialog show is for the Metrics dialog, not the strip dock: `src/app/ReaderMainWindow.cpp:1182`.

### Concrete Fix

Add one UI-side reveal helper in `ReaderMainWindow`, and call it after a successful metric display add or display-mode enable:

```cpp
void ReaderMainWindow::showDashboardStripDock()
{
    if (!dashboardStripDock_) return;
    dashboardStripDock_->setVisible(true);
    resizeDocks({dashboardStripDock_}, {360}, Qt::Horizontal);
    dashboardStripDock_->raise();
    QTimer::singleShot(0, this, [this]() {
        if (!dashboardStripDock_) return;
        resizeDocks({dashboardStripDock_}, {360}, Qt::Horizontal);
        dashboardStripDock_->raise();
    });
}
```

Best trigger points:

- Preferred: have `SignalDisplayDialog` emit a `displayAdded(panelId, signalId, addedRefs)` signal from `onAddSelected()` after `addedRefs > 0`, and connect that to `showDashboardStripDock()`.
- Also connect display-mode enable paths that add refs (`SignalDisplayDialog::onActiveModeToggled()` or `DashboardPanelModel::displayRefAdded`) so enabling a mode later is visible too.
- If using `DashboardPanelModel::displayRefAdded`, keep the reveal in `ReaderMainWindow`, not the model; the model should remain UI-neutral.

This fixes the live bug for supported strip/static-panel display modes. It also gives a single place to harden the all-tabified/all-hidden Qt case with a queued resize/raise.

### Additional Metric-Path Dead Ends

There are display modes that the dialog can accept but that do not currently produce a visible surface:

- `static.tensor` gets a `"panel"` ref in `DisplayRefsForSignal()`, but `DashboardDisplayController::isPanelMode()` excludes it: `src/model/DashboardPanelModel.cpp:67`, `src/model/DashboardPanelModel.cpp:71`, `src/app/DashboardDisplayController.cpp:63`, `src/app/DashboardDisplayController.cpp:68`. The controller explicitly says the tensor-glyph trigger is intentionally omitted until a future explicit gesture: `src/app/DashboardDisplayController.cpp:1645`, `src/app/DashboardDisplayController.cpp:1655`, `src/app/DashboardDisplayController.cpp:1690`, `src/app/DashboardDisplayController.cpp:1693`.
- Other static modes such as `static.table`, `static.atomColor`, `static.scalar`, `static.geometry`, `static.system`, and `static.category` are advertised in descriptors but are neither `strip.*` nor in `isPanelDisplayMode()`, so `DisplayRefsForSignal()` skips them and `addedRefs` can be `0`: `src/model/DashboardPanelModel.cpp:86`, `src/model/DashboardPanelModel.cpp:90`, `src/app/SignalDisplayDialog.cpp:1163`, `src/app/SignalDisplayDialog.cpp:1169`.

Concrete fix for those: do not offer modes with no visible implementation in the Add row, or route them to a visible placeholder/table/color surface. `static.tensor` needs either a real reveal gesture or should not be presented as an Add-visible mode.

## 2. Dock / Panel Visibility

### Flow Map

- `View -> Panels` menu is created in `buildUi()`: `src/app/ReaderMainWindow.cpp:699`, `src/app/ReaderMainWindow.cpp:700`, `src/app/ReaderMainWindow.cpp:701`.
- The three dock widgets are tabified together: `src/app/ReaderMainWindow.cpp:372`, `src/app/ReaderMainWindow.cpp:387`.
- All three start hidden: `src/app/ReaderMainWindow.cpp:391`, `src/app/ReaderMainWindow.cpp:396`, `src/app/ReaderMainWindow.cpp:398`.
- Each dock contributes its own `QDockWidget::toggleViewAction()` to the menu: `src/app/ReaderMainWindow.cpp:403`, `src/app/ReaderMainWindow.cpp:415`, `src/app/ReaderMainWindow.cpp:425`.
- The Tools toolbar's Panels button points at the same menu/actions: `src/app/ReaderMainWindow.cpp:428`, `src/app/ReaderMainWindow.cpp:430`, `src/app/ReaderMainWindow.cpp:434`, `src/app/ReaderMainWindow.cpp:435`.
- On `visibilityChanged(true)`, the app resizes the dock to 360 and raises it: `src/app/ReaderMainWindow.cpp:418`, `src/app/ReaderMainWindow.cpp:420`, `src/app/ReaderMainWindow.cpp:422`, `src/app/ReaderMainWindow.cpp:423`.

### Break / Tangle

The visibility recovery path only runs when the user explicitly toggles a panel. It is not called by the Add flow.

From code alone, the intended toggle path should resize and raise a toggled dock. What cannot be settled without runtime is whether Qt's all-tabified-then-all-hidden state defeats `setVisible(true)` or leaves the dock with unusable width. The code has the right `visibilityChanged(true)` mitigation, but only a runtime check can confirm it:

- Trigger View -> Panels -> Strip from all docks hidden.
- Assert `dashboardStripDock_->isVisible() == true`.
- Assert width is approximately 360 or otherwise usable.
- Assert the Strip tab is the raised tab.

If that fails, put the queued resize/raise helper above behind both the user toggle and the metric Add reveal.

The REST/harness path preserves prior visibility, so it is not a recovery path from startup-hidden docks: hide stashes `isVisible()`, restore re-applies the stashed value: `src/app/ReaderMainWindow.cpp:595`, `src/app/ReaderMainWindow.cpp:610`, `src/app/ReaderMainWindow.cpp:611`, `src/app/ReaderMainWindow.cpp:623`, `src/app/ReaderMainWindow.cpp:627`.

## 3. Pick -> Selection -> Inspector

### Flow Map

- `QtAtomPicker` is installed as an event filter on the VTK widget: `src/app/QtAtomPicker.cpp:33`, `src/app/QtAtomPicker.cpp:48`.
- Double-click is consumed and converted to a pick: `src/app/QtAtomPicker.cpp:55`, `src/app/QtAtomPicker.cpp:57`, `src/app/QtAtomPicker.cpp:59`, `src/app/QtAtomPicker.cpp:60`.
- Picking ray-tests all atoms at the current frame through the transformed conformation: `src/app/QtAtomPicker.cpp:94`, `src/app/QtAtomPicker.cpp:101`, `src/app/QtAtomPicker.cpp:103`, `src/app/QtAtomPicker.cpp:116`, `src/app/QtAtomPicker.cpp:125`.
- `ReaderMainWindow` wires `atomPicked` to `AtomSelection::applyPick`, scene reveal clearing, and a picker render tag: `src/app/ReaderMainWindow.cpp:308`, `src/app/ReaderMainWindow.cpp:311`, `src/app/ReaderMainWindow.cpp:318`, `src/app/ReaderMainWindow.cpp:321`.
- Plain pick replaces the selection and emits `changed()` plus `focusChanged(atomIdx)`: `src/model/AtomSelection.cpp:131`, `src/model/AtomSelection.cpp:144`, `src/model/AtomSelection.cpp:150`, `src/model/AtomSelection.cpp:154`, `src/model/AtomSelection.cpp:155`.
- Shift-pick toggles membership and emits the same fanout or `cleared()`: `src/model/AtomSelection.cpp:159`, `src/model/AtomSelection.cpp:170`, `src/model/AtomSelection.cpp:171`, `src/model/AtomSelection.cpp:201`, `src/model/AtomSelection.cpp:202`.
- Inspector is retargeted on focus and rebuilt: `src/app/ReaderMainWindow.cpp:324`, `src/app/ReaderMainWindow.cpp:325`, `src/app/QtAtomInspectorDock.cpp:212`, `src/app/QtAtomInspectorDock.cpp:218`, `src/app/QtAtomInspectorDock.cpp:245`, `src/app/QtAtomInspectorDock.cpp:251`.
- Selection dock is a model-bound display surface: `src/app/SelectionDock.cpp:65`, `src/app/SelectionDock.cpp:69`, `src/app/SelectionDock.cpp:71`, `src/app/SelectionDock.cpp:77`.
- Measurement overlay is wired to the same `changed()` signal and scene refresh: `src/app/ReaderMainWindow.cpp:343`, `src/app/ReaderMainWindow.cpp:345`, `src/app/ReaderMainWindow.cpp:348`, `src/app/ReaderMainWindow.cpp:363`.
- Measurement overlay is visible by default and shows selected slots as sphere actors: `src/app/MeasurementOverlay.h:115`, `src/app/MeasurementOverlay.cpp:242`, `src/app/MeasurementOverlay.cpp:250`, `src/app/MeasurementOverlay.cpp:262`, `src/app/MoleculeScene.cpp:424`, `src/app/MoleculeScene.cpp:429`.

### Break / Tangle

Picking does not reveal the Inspector or Selection docks. This is explicit in the startup-hidden comment: picking updates Inspector contents but does not force-show the dock: `src/app/ReaderMainWindow.cpp:391`, `src/app/ReaderMainWindow.cpp:393`.

The path is not completely silent because the measurement overlay should show in the viewport after a successful pick. However, the detailed inspector path is hidden by design. If the lead's rule is that a successful pick should surface the inspector content, then this path needs the same style of UI reveal helper for `inspectorDock_`.

Concrete fix: on first successful pick or on `AtomSelection::focusChanged`, reveal Inspector if it is currently hidden, or show an unobtrusive status/toolbar indicator that the Inspector has updated. Keep the measurement overlay as the viewport-level confirmation.

## 4. Camera

### Flow Map

- Camera input filter is installed after the picker, with double-click left to the picker: `src/app/ReaderMainWindow.cpp:261`, `src/app/ReaderMainWindow.cpp:264`, `src/app/CameraInputFilter.cpp:43`, `src/app/CameraInputFilter.cpp:49`.
- Mouse drag and wheel become `CameraGesture`s: `src/app/CameraInputFilter.cpp:95`, `src/app/CameraInputFilter.cpp:103`, `src/app/CameraInputFilter.cpp:108`, `src/app/CameraInputFilter.cpp:119`, `src/app/CameraInputFilter.cpp:127`, `src/app/CameraInputFilter.cpp:144`, `src/app/CameraInputFilter.cpp:152`.
- Gesture path applies deltas, syncs clipping, and requests render immediately: `src/app/CameraComposer.cpp:722`, `src/app/CameraComposer.cpp:785`, `src/app/CameraComposer.cpp:844`, `src/app/CameraInputFilter.cpp:133`, `src/app/CameraInputFilter.cpp:135`, `src/app/CameraInputFilter.cpp:153`, `src/app/CameraInputFilter.cpp:155`.
- Frame tick path pushes atom positions, writes camera state, updates overlays, syncs clipping, and schedules render: `src/app/MoleculeScene.cpp:377`, `src/app/MoleculeScene.cpp:393`, `src/app/MoleculeScene.cpp:410`, `src/app/MoleculeScene.cpp:416`, `src/app/MoleculeScene.cpp:424`, `src/app/MoleculeScene.cpp:438`, `src/app/MoleculeScene.cpp:441`.
- Clipping sync uses cached padded per-frame bounds, not stale VTK actor bounds: `src/app/MoleculeScene.cpp:365`, `src/app/MoleculeScene.cpp:373`, `src/app/MoleculeScene.cpp:302`, `src/app/MoleculeScene.cpp:306`.

### Break / Tangle

Gesture camera control is visible immediately.

Toolbar camera modes are inconsistent:

- Plane lock writes the camera, syncs clipping, and renders immediately: `src/app/MoleculeScene.cpp:469`, `src/app/MoleculeScene.cpp:474`, `src/app/MoleculeScene.cpp:487`, `src/app/MoleculeScene.cpp:488`, `src/app/MoleculeScene.cpp:489`.
- Focus and Newman actions only call `CameraComposer::setMode(...)`: `src/app/ReaderMainWindow.cpp:997`, `src/app/ReaderMainWindow.cpp:1014`, `src/app/ReaderMainWindow.cpp:1017`, `src/app/ReaderMainWindow.cpp:1027`.
- `CameraComposer::setMode()` captures state and emits `modeChanged()`, but it does not write the camera or request render: `src/app/CameraComposer.cpp:46`, `src/app/CameraComposer.cpp:103`, `src/app/CameraComposer.cpp:109`.

Concrete fix: after Focus/Newman `setMode()`, mirror the plane-lock immediate apply path: `composer->write(currentFrame)`, `scene_->syncCameraClippingRange()`, `scene_->requestRender(...)`, then update action state. Otherwise a paused user can click Focus/Newman and see no camera change until the next frame tick or gesture.

## 5. Transform

### Flow Map

- The loaded conformation is wrapped in `TransformedConformation`; startup defaults to backbone fit if at least three backbone atoms exist: `src/app/ReaderMainWindow.cpp:163`, `src/app/ReaderMainWindow.cpp:168`, `src/app/ReaderMainWindow.cpp:169`, `src/app/ReaderMainWindow.cpp:172`, `src/app/ReaderMainWindow.cpp:173`.
- Fallback to all-atom fit is explicit when backbone fit is underdetermined: `src/app/ReaderMainWindow.cpp:174`, `src/app/ReaderMainWindow.cpp:177`, `src/app/ReaderMainWindow.cpp:180`.
- Scene, picker, inspector, metric dialog, and dashboard all receive the transformed conformation pointer: `src/app/ReaderMainWindow.cpp:199`, `src/app/ReaderMainWindow.cpp:256`, `src/app/ReaderMainWindow.cpp:258`, `src/app/ReaderMainWindow.cpp:271`, `src/app/ReaderMainWindow.cpp:302`, `src/app/ReaderMainWindow.cpp:380`.
- The toolbar action toggles all-atom vs backbone fit: `src/app/ReaderMainWindow.cpp:800`, `src/app/ReaderMainWindow.cpp:805`, `src/app/ReaderMainWindow.cpp:1040`, `src/app/ReaderMainWindow.cpp:1046`, `src/app/ReaderMainWindow.cpp:1059`.
- `TransformedConformation::atomPosition()` computes/caches a per-frame transform and applies `R * raw + T`: `src/model/TransformedConformation.cpp:51`, `src/model/TransformedConformation.cpp:58`, `src/model/TransformedConformation.cpp:62`, `src/model/TransformedConformation.cpp:63`.
- `setMode()` clears the transform cache, rebuilds reference positions, and emits `transformChanged()`: `src/model/TransformedConformation.cpp:76`, `src/model/TransformedConformation.cpp:89`, `src/model/TransformedConformation.cpp:95`, `src/model/TransformedConformation.cpp:125`.
- Transform changes force a current-frame refresh, making scene/overlays/picker-space positions update without waiting for playback: `src/app/ReaderMainWindow.cpp:183`, `src/app/ReaderMainWindow.cpp:190`, `src/app/MoleculeScene.cpp:337`, `src/app/MoleculeScene.cpp:344`.

### Break / Tangle

No transform visibility break found in the source path. The main caveat is an intentional split: snapshots are forwarded raw, while `atomPosition()` is transformed: `src/model/TransformedConformation.cpp:66`, `src/model/TransformedConformation.cpp:72`, `src/model/TransformedConformation.cpp:73`. Consumers that mix snapshot fields with positions need to keep that distinction in mind.

## 6. Overlay Toggles

### Flow Map

- Toolbar creates Ribbon, Rings, Butterfly, and B-field actions; Ribbon/Rings default on, Butterfly/B-field default off: `src/app/ReaderMainWindow.cpp:831`, `src/app/ReaderMainWindow.cpp:833`, `src/app/ReaderMainWindow.cpp:837`, `src/app/ReaderMainWindow.cpp:839`, `src/app/ReaderMainWindow.cpp:843`, `src/app/ReaderMainWindow.cpp:845`, `src/app/ReaderMainWindow.cpp:850`, `src/app/ReaderMainWindow.cpp:852`.
- Connections are deferred until scene/playback exist: `src/app/ReaderMainWindow.cpp:857`, `src/app/ReaderMainWindow.cpp:862`, `src/app/ReaderMainWindow.cpp:869`.
- Ribbon and Rings just set actor visibility and request render: `src/app/ReaderMainWindow.cpp:874`, `src/app/ReaderMainWindow.cpp:877`, `src/app/ReaderMainWindow.cpp:878`, `src/app/ReaderMainWindow.cpp:880`, `src/app/ReaderMainWindow.cpp:883`, `src/app/ReaderMainWindow.cpp:884`.
- Butterfly and B-field set visibility and call `refreshCurrentFrame()` when turning on so skipped expensive data updates run before render: `src/app/ReaderMainWindow.cpp:891`, `src/app/ReaderMainWindow.cpp:894`, `src/app/ReaderMainWindow.cpp:895`, `src/app/ReaderMainWindow.cpp:903`, `src/app/ReaderMainWindow.cpp:906`, `src/app/ReaderMainWindow.cpp:907`.
- MoleculeScene owns and builds all overlays: `src/app/MoleculeScene.cpp:244`, `src/app/MoleculeScene.cpp:252`, `src/app/MoleculeScene.cpp:254`, `src/app/MoleculeScene.cpp:257`, `src/app/MoleculeScene.cpp:259`, `src/app/MoleculeScene.cpp:262`.
- Per-frame fanout updates every overlay then renders once: `src/app/MoleculeScene.cpp:424`, `src/app/MoleculeScene.cpp:425`, `src/app/MoleculeScene.cpp:428`, `src/app/MoleculeScene.cpp:441`.

### Break / Tangle

No UI-code break found for the toolbar overlay toggles.

The expensive overlays rely on their owner path to force a same-frame refresh:

- `QtFieldGridOverlay::setFrame()` skips work while hidden: `src/app/QtFieldGridOverlay.cpp:228`, `src/app/QtFieldGridOverlay.cpp:233`.
- `QtBFieldStreamOverlay::setFrame()` also skips work while hidden: `src/app/QtBFieldStreamOverlay.cpp:254`, `src/app/QtBFieldStreamOverlay.cpp:258`.
- Their `setVisible(true)` methods only flip visibility/state; they do not self-refresh: `src/app/QtFieldGridOverlay.cpp:293`, `src/app/QtFieldGridOverlay.cpp:299`, `src/app/QtBFieldStreamOverlay.cpp:291`, `src/app/QtBFieldStreamOverlay.cpp:295`.

That is fine for the current toolbar path because `ReaderMainWindow` calls `refreshCurrentFrame()` on enable. If future callers toggle these overlays directly, they must also refresh the scene or the actors can become visible before current-frame data is populated.

## Prioritized Fix List

1. Reveal Strip dock after successful metric display add.
   - Add a `ReaderMainWindow::showDashboardStripDock()` helper that `setVisible(true)`, resizes to 360, raises, and repeats resize/raise on a zero-delay queued callback.
   - Call it from a new `SignalDisplayDialog` success signal or from `DashboardPanelModel::displayRefAdded`.
   - Include active-mode toggles that add refs, not just first Add.

2. Runtime-check all-tabified/all-hidden dock recovery.
   - From all docks hidden, trigger View -> Panels -> Strip and assert visible, raised, and usable width.
   - If Qt defeats direct visibility, route all dock reveals through the queued helper and use that for toolbar/menu toggles too.

3. Stop accepting display modes that cannot surface visibly.
   - Hide/disable `static.table`, `static.atomColor`, `static.scalar`, etc. until they have actual surfaces.
   - Either implement a `static.tensor` reveal gesture or stop presenting it as a visible Add mode.

4. Make Focus/Newman camera actions visible immediately.
   - After `setMode()`, call `write(currentFrame)`, sync clipping, and request render, matching plane lock.

5. Decide whether pick should reveal Inspector.
   - The viewport marker already gives visible feedback, but the inspector content stays hidden by design. If advisor-facing behavior expects details to appear, reveal Inspector on first focus pick.

6. Keep expensive overlay refresh discipline centralized.
   - Current toolbar path is correct. Future direct callers of `setVisible(true)` for Butterfly/B-field must also force `scene_->refreshCurrentFrame()`.
