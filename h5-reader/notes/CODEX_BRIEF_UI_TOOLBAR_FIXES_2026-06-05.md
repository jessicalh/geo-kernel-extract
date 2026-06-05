# Codex — three reader UI-polish fixes (lead's live test feedback, 2026-06-05)

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch h5-reader-pysr-spike). Edit only
`src/app/ReaderMainWindow.cpp` + `src/app/ReaderMainWindow.h`. No git. qt6 discipline (ACONNECT, ASSERT_THREAD;
no new persistent QTimer). Build to verify: `cmake --build build/linux-rwdi --target h5reader -j$(nproc)`. Do
NOT launch the GUI — Claude relaunches. The lead is testing live; the three diagnoses below are confirmed at
the cited lines.

## 1. Split the single toolbar into TWO rows (playback vs tools)
`buildToolbar()` (line 742) builds ONE toolbar `tb` (`addToolBar("Playback")`, line 743) holding everything.
The lead wants the playback bar and the action toolbar **separate**.
- Keep the **Playback** toolbar = Play/Pause + Step back/fwd + frame slider + fps (lines 754-777).
- After the fps block (~line 779), do `addToolBarBreak();` then create a SECOND toolbar and **reassign `tb` to
  it** (`tb = addToolBar("Tools"); tb->setObjectName("ToolsToolbar"); tb->setMovable(false);
  tb->setFont(toolbarFont);`) and store it in a new member `toolsToolbar_`. Because the rest of `buildToolbar`
  already uses `tb->...`, every subsequent item (camera-mode group Focus/Newman/Plane/Free, the "All-atom fit"
  toggle, the hidden Harness-marker, "Metrics...", the overlay toggles) then lands on the Tools row with no
  further changes. (Leave `playbackToolbar_` pointing at the first toolbar.)
- The **Panels** button is added later in the ctor body (~lines 466-473) to `playbackToolbar_`; point it at
  `toolsToolbar_` instead so it sits with the other actions.

## 2. Panels don't appear when toggled
The three docks start hidden (lines 441-443) and were given `minimumSize(0,0)` on the dock AND inner widget
(loop ~lines 404-410, the "thin tab-only strip" override); with the `kSettingsVersion` bump discarding saved
width, toggling a panel visible shows it at ~0 px → "nothing happens."
- **Remove the `minimumSize(0,0)` override loop** so docks keep a natural minimum width.
- AND make a toggled-on panel actually appear usefully: wire each dock's `toggleViewAction()` `toggled(bool)`
  (or the dock's `visibilityChanged(bool)`) so that on becoming visible the window does
  `resizeDocks({dock}, {360}, Qt::Horizontal);` + `dock->raise();`. Verify: click **Panels → Inspector** and
  the Inspector dock appears at a usable width.

## 3. "Free" appears twice
`buildToolbar()` adds a checkable **Free** camera-mode radio (lines 820-827, checked by default) AND a
`cameraStateLabel_` QLabel showing "Free" immediately beside it (lines 829-831); the label is also driven with
"Locked: …"/"Follow: …" text in `updateCameraModeActions` (~lines 1003-1020).
- **Remove `cameraStateLabel_`** entirely: its creation (829-831), every `cameraStateLabel_->setText(...)` use
  in `updateCameraModeActions`, and the member in `ReaderMainWindow.h`. **Keep** the camera-mode radio group
  (Focus/Newman/Plane/Free) — `Free` is the functional "release the lock" radio in the exclusive group and
  must stay; the active mode is already shown by which radio is checked. This removes the duplicate "Free".

## 4. No clean "no metric selected" initial state
The reader **seeds a default metric on load** — `addInitialGenericDashboardSignal(..., "npy:dssp_chi", ...)`
(helper ~lines 99-119; the call is ~line 326 in the ctor body). So the metric box always opens with that
default showing/selected; there is no clean unselected state. The lead wants the reader to open with **no
metric selected** — empty dashboard, nothing pre-chosen in the picker.
- **Remove the default-metric seed**: delete the `addInitialGenericDashboardSignal(...)` call (and the helper
  function + the `DashboardSignalPanelCoordinator` seed path if they become unused). The surgery only *gated*
  this seed on availability; the lead wants it gone outright.
- Confirm the metric-selection box / dashboard strip opens with **nothing selected** (no default row, no
  pre-checked item). If a combo/list defaults to the first entry, set it to an explicit "none" initial state.

## Boundaries note
Issues 1-3 are confined to `ReaderMainWindow.{cpp,h}`. Issue 4's seed removal is also in `ReaderMainWindow`,
but if clearing the picker's default selection needs a small touch in the dashboard model / `SignalDisplayDialog`,
that is allowed — keep it minimal and report it. Still: no git, build to verify, do not launch.

## Report
Per fix: what changed + how to test. Build result. Anything you couldn't do in the allowed files.
