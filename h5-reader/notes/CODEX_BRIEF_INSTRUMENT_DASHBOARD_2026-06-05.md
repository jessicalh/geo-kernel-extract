# Codex — instrument the dashboard selection path (REST + structured logging)

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`). **Portable Qt6
(Windows/macOS/Linux) — the cross-platform source is first-class, treat it as real.** Scope: `h5-reader/src`
ONLY. Never read/link the `nmr_shielding` library, never write H5, never trigger extraction.

**Verify from the code before changing anything, ever.** This brief cites `file:line`, but the code is the
truth — if a cite is stale, follow the code and tell me what you found.

## Why this matters (the stakes, honestly)
We're about to do a careful, staged rebuild of the selected-metrics model (the "add a metric → see a panel"
path that currently fails silently). To do that without four human round-trips at the screen, we need to
**drive and observe** the metric/panel UI headlessly. You're building that instrument first. It is purely
additive observability — **no behavior change to any existing UI path** — and it becomes the harness every
later step is verified against. Get this right and the rest goes fast and safe.

## Read these first (the real paths — start here, confirm before you build)
- **RestServer** (`src/app/RestServer.{h,cpp}`): a `QHttpServer` on loopback. Handlers run on the GUI thread
  (`ASSERT_THREAD(this)`) so they read/mutate models directly. JSON via `QJsonDocument/Object/Array` only.
  Helpers: `jsonResponse()`, `errorResponse(msg, status)`. It already holds
  `QPointer<model::DashboardSignalModel> signalModel_`, `QPointer<model::DashboardPanelModel> panelModel_`,
  `QPointer<ReaderMainWindow> readerWindow_`. Existing route **GET `/dashboard/signals`**
  (`RestServer.cpp:1046`) already serializes the selected list — model your new code on it.
- **The real ADD path** — read it and mirror it EXACTLY; do NOT invent a parallel add:
  `SignalDisplayDialog::onAddSelected` (`src/app/SignalDisplayDialog.cpp:1106`) →
  `model::DashboardSignalModel::addSignal(...)` (`src/model/DashboardSignalModel.cpp:233`) →
  `model::DisplayRefsForSignal(...)` (`src/model/DashboardPanelModel.cpp:74`) →
  `DashboardPanelModel::addDisplayRef(...)`. **Order matters**: `signalAdded` fires first (controller
  rebuilds with no refs yet), then `displayRefsChanged` fires and the controller rebuilds again — this
  two-step is documented at `DashboardDisplayController.cpp:1229-1239`. Your REST add must reproduce that
  order so the rebuild lands. (When a later step introduces a named `DashboardSelectionController`, this REST
  handler becomes its second caller — so faithfulness now pays off then.)
- **The cleanup cascade**: anonymous `DashboardSignalPanelCoordinator` (`src/app/ReaderMainWindow.cpp:102-146`)
  — `signalRemoved → removeDisplayRefsForSignal`; `displayRefRemoved → removeSignal` when refcount==0. REST
  remove should call `removeSignal` and let the coordinator cascade (don't duplicate the cascade).
- **Dock reveal**: `ReaderMainWindow::revealDockQueued(QDockWidget*)` (`src/app/ReaderMainWindow.cpp:633`) +
  `dashboardStripDock_`. Docks start hidden (`ReaderMainWindow.cpp:396-398`).
- **Logging**: structured logger → stderr + UDP 9997; per-module `Q_LOGGING_CATEGORY` is the idiom
  (`Q_LOGGING_CATEGORY(cRest, "h5reader.rest")` at `RestServer.cpp:57`).

## Deliverable 1 — REST routes (the observable surface)
Add to `RestServer::registerRoutes()` (+ the `RestServer.h` doc-comment route table). Each handler:
`ASSERT_THREAD(this)`, guard null models with `errorResponse(..., ServiceUnavailable)`, QJson only.

- **GET `/dashboard/state`** — the single source of observable truth:
  ```json
  { "selected": [ { "id","descriptor_id","concept_key","label","anchor":{...},"enabled",
                    "availability":"<state-or-unknown>",
                    "modes":[ {"mode","enabled","is_panel_ref","renderable_panel","has_visible_surface"} ] } ],
    "selected_count": N,
    "dock":   { "visible": bool, "width": int, "raised": bool },
    "render": { "owned_panel_count": int, "strip_track_count": int } }
  ```
  - `availability`: resolve each selected descriptor against the loaded-run availability layer
    (`TrajectoryFieldAvailability` — reuse how `TrajectorySignalCatalog` consults `availability_`; do NOT
    re-derive). If you can't reach availability cleanly from RestServer yet, return `"unknown"` with a
    `// TODO(availability-step)` — do NOT fabricate a state.
  - per-mode flags — report ALL THREE truthfully (they disagree on purpose, see the note below; exposing the
    disagreement IS the instrument's job — do NOT reconcile here): `is_panel_ref` = the panel-model panel-ref
    predicate; `renderable_panel` = the controller's `isPanelMode`; `has_visible_surface` = the dialog's
    `modeHasVisibleSurface`.
  - `render` counts: read how `DashboardStripDock` receives `ownedPanelsChanged` / `stripTracksChanged` and
    where the live counts sit (`StripStackWidget` / `DashboardDisplayController`). Expose a `const` getter if
    none exists. These two ints let a headless caller assert "a panel actually got built," without pixels.
- **POST `/dashboard/metric`** `{"descriptor_id":"...","anchor":{"atom":int,"frame":int}|{"follows_focus":true},"modes":["strip....",...]}`
  — resolve the descriptor via the catalog (read the catalog's lookup), build the anchor, then reproduce the
  dialog add (`addSignal` → `DisplayRefsForSignal` → `addDisplayRef`, that order). Return
  `{"id":uuid,"added_refs":int}`. Descriptor not found / not available → `errorResponse` with the reason.
- **POST `/dashboard/metric/remove`** `{"id":uuid}` — `removeSignal(id)`; coordinator cascades. Return
  `{"removed":true}`.
- **POST `/dashboard/metric/mode`** `{"id":uuid,"mode":"...","enabled":bool}` — mirror
  `SignalDisplayDialog::onActiveModeToggled` (read it): mutate the signal row's modes, then add/remove the
  panel refs in the active panel. Return the new mode list.
- **POST `/dashboard/dock`** `{"visible":bool}` — show via `revealDockQueued(dashboardStripDock_)`, hide via
  `setVisible(false)`. Return `{"visible":bool,"width":int}`.

## Deliverable 2 — structured logging (the transition trace)
Add `Q_LOGGING_CATEGORY(cDash, "h5reader.dashboard")` (declare where the other categories live; pick the
right home — likely `DashboardDisplayController.cpp`). Emit ONE greppable, structured line per transition
(`qCInfo(cDash).noquote() << "key=val ..."`), hooked **where the transition actually occurs** (model
signals / coordinator / controller rebuild tail / dock reveal) so **both UI-driven and REST-driven** changes
log — NOT in the REST layer:
- signal added (id, descriptor_id, modes, anchor) / signal removed (id)
- display mode toggled (id, mode, enabled)
- panel refs changed (signal id, ref_count)
- controller rebuild outcome (owned_panel_count, strip_track_count, active_panel)
- dock visibility changed (visible, width)

A log tail must be able to prove the chain: `add X → selected=[X] → refs=N → rebuild panels=1 tracks=0 →
dock visible width=360`.

## A real subtlety — don't paper over it
Three predicates classify a display mode and they DISAGREE (verified in code):
- `modeHasVisibleSurface` (`SignalDisplayDialog.cpp:218`) — dialog "can the user pick it" gate
- `isPanelMode` (`DashboardDisplayController.cpp:63`) — controller "build an `AbstractStripPanel` widget"
- `isPanelDisplayMode` (`DashboardPanelModel.cpp:61`) — panel-model "emit a `panel` sentinel ref"

`static.tensor`: ref=yes, panel-widget=no, surface=no (intentional — deferred scene overlay).
`static.spectrum.power`: ref=yes, panel-widget=**YES**, surface=**NO** (a renderable panel the dialog greys).
For THIS step just **expose all three truthfully** in `/dashboard/state`. Reconciling them is a later step;
reporting the disagreement is the point.

## Discipline
`CENSUS_REGISTER` in any new QObject ctor; `ACONNECT` (not raw `connect`) for any new connection;
`ASSERT_THREAD(this)` on thread-sensitive methods; **no new `QTimer`** (QHttpServer + the existing
`revealDockQueued` singleShot are fine — don't add polling); full error handling at the JSON boundary;
portable (no platform-specific calls); match the file's existing idioms.

## Output
Implement, then build green: `cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi
--target h5reader -j$(nproc)`. **Do NOT launch, do NOT git.** Then write a short
`notes/INSTRUMENT_DASHBOARD_2026-06-05.md`: routes added (with example `curl`), log lines (with example
tail), and anything in the code that contradicted this brief (code wins — tell me).
