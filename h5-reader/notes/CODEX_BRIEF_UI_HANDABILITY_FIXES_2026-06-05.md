# Codex — reader handability fixes (lead-revised, opus-verified)

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch h5-reader-pysr-spike). Edit `h5-reader/src`
(app + the dialog + catalog + the dock widget). No git. Build to verify: `cmake --build build/linux-rwdi
--target h5reader -j$(nproc)`. Do NOT launch — Claude relaunches and runtime-verifies. Context: the reader
"mostly works but is too funky to hand off"; the opus-verified UI-paths map (`notes/UI_PATHS_REVIEW_2026-06-05.md`)
found the breaks. **The lead revised the panels fix: do NOT auto-pop the dock — instead make the dock SELECT
(the Panels toggle) reliably work.**

## 1. Make the dock-select reliably show the dock (NOT auto-pop on Add)
The Add flow renders the panel into `DashboardStripDock`, which starts hidden — that is FINE; the user shows
it via Panels→Strip. The real bug is that the toggle does not surface a USABLE dock. Root cause (opus
verification): `DashboardStripDock` sets `setSizePolicy(Expanding, Ignored)` with **no `setMinimumWidth`**
(only `setMinimumHeight(64)`), so when toggled visible from the all-tabified-then-all-hidden state it collapses
to ~0 px width → invisible. The round-1 `visibilityChanged(true)` → resize-to-360 + raise mitigation is not
reliably landing.
- **Give the property docks a sensible `setMinimumWidth` (~240-280 px)** so they cannot collapse to 0.
- **Make the Panels-toggle reveal reliable:** route dock reveals through a queued helper — `setVisible(true)`
  then `QTimer::singleShot(0, …)` → `resizeDocks({dock},{360},Qt::Horizontal)` + `raise()` — used by the
  View→Panels menu and the Tools-toolbar Panels button. Handle the all-tabified-then-all-hidden Qt quirk so
  toggling one dock actually shows it raised at usable width.
- **Do NOT auto-reveal the dock on metric Add.** The panel renders into the dock; the user reveals it via
  Panels. Goal: from all-hidden, click Panels→Strip (or Inspector/Selection) → the dock appears, raised,
  ~360 wide, painting.

## 2. Disable the show-nothing mode checkboxes (co-equal with #1)
Opus confirmed: `static.table` is offered + ENABLED on nearly every descriptor and renders nothing;
`static.atomColor` / `static.tensor` / `static.scalar` / `static.geometry` / `static.system` /
`static.category` likewise have no visible surface (addedRefs=0, or excluded by `isPanelMode`, or the
tensor-glyph trigger is intentionally omitted). So a user can check a mode, Add, and see nothing even with the
dock shown.
- In `SignalDisplayDialog::onCandidateSelectionChanged` (the `control.box->setEnabled(supported)` path, ~lines
  1061-1064), treat the modes with NO visible implementation as **not supported** — do not enable their
  checkboxes. Keep enabled only the modes that actually render: the `strip.*` modes plus any static mode with
  a real surface. Adjust the default-check fallback (~lines 1071-1079) so it never lands on a dead-end mode.
- Verify a normal strip-capable descriptor still has an enabled, working mode (a field whose only option is a
  dead-end mode should honestly show NO enabled mode rather than a lying checkbox).

## 3. Focus / Newman camera — apply on click
Plane-lock applies immediately (`MoleculeScene.cpp:487-489`: `write` + `syncCameraClippingRange` +
`requestRender`). Focus (`ReaderMainWindow.cpp:1014`) and Newman (`:1027`) call only `composer->setMode(...)`,
which emits `modeChanged` but does NOT write/render → on a paused frame they do nothing visible.
- After `setMode` in `onFocusCameraTriggered` / `onNewmanProjectionTriggered`, mirror plane-lock:
  `composer->write(currentFrame)` + `scene_->syncCameraClippingRange()` + `scene_->requestRender(...)`, then
  update action state.

## Do NOT
- Do NOT auto-reveal the inspector on pick (deferred — the measurement-overlay sphere already gives feedback).
- Do NOT auto-pop the Strip dock on Add (see #1).

## Boundaries
`h5-reader/src` (app + dialog + catalog + dock). No git. qt6 discipline (ACONNECT, ASSERT_THREAD; the
`singleShot(0)` queued calls are one-shots, allowed — no new PERSISTENT QTimer). Build to verify. Do NOT
launch. Report per-fix + how to test each.
