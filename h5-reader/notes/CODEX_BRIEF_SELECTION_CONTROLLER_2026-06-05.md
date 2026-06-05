# Codex — STEP 2.5a: the named DashboardSelectionController (single owner of mutations)

## FIRST — read the qt6-cpp skill before touching any code (full filesystem access — open them)
- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md` ← ownership, signal/slot, object lifetime, the "one owner" principle
- `/home/jessica/.claude/skills/qt6-cpp/references/model-view.md`  ← model mutation contracts, dataChanged

Honor it: `CENSUS_REGISTER` in the new QObject ctor, `ACONNECT` (never raw `connect`), `ASSERT_THREAD(this)`
on mutators, **no new `QTimer`**, portable Qt (Windows/macOS/Linux first-class). This is the highest-stakes
change in the series — a coordination refactor — so the architecture reference is the lens.

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`; floor `8007e2e` +
uncommitted 2.2/2.3). Scope: `h5-reader/src` ONLY. Never link `nmr_shielding`, never write H5.
**Verify from the code before changing anything, ever.** This refactor will be opus-adversarially reviewed
and instrument-verified before it's trusted — so make it faithful and obvious, not clever.

## Why this matters
Two models — `DashboardSignalModel` (what's selected) and `DashboardPanelModel` (where it's placed) — are
kept in sync by an **anonymous coordinator buried in `ReaderMainWindow.cpp`**, and callers write BOTH models
directly in two-step sequences that can drift. That hidden two-owner coupling is the rot behind "add a metric,
nothing happens." This step gives those compound mutations **one named owner** so every caller goes through a
single, testable door. **It is a pure structural refactor: behavior must be preserved exactly** — same adds,
same removes, same cascade, same rendering. No new UX yet (the dock reveal/reset rules are the NEXT step).

## Read first and confirm (cite what you find)
- The anonymous `DashboardSignalPanelCoordinator` — `src/app/ReaderMainWindow.cpp` (~lines 99-147): it wires
  `signalRemoved → panelModel.removeDisplayRefsForSignal` and `displayRefRemoved → signalModel.removeSignal`
  when the signal's ref count hits 0, with a `signalsBeingRemoved_` reentrancy guard. **This logic is what
  moves into the controller — preserve its exact semantics, including the guard.**
- The dialog's compound writes (the two-step pattern): `SignalDisplayDialog::onAddSelected` (addSignal →
  `DisplayRefsForSignal` → `addDisplayRef`), `onRemoveActive`, `onActiveModeToggled` (mutate modes → add/remove
  refs). The ordering matters — `signalAdded` fires before refs exist, then `displayRefsChanged` drives the
  controller's rebuild (`DashboardDisplayController` documents this two-step at the rebuild wiring).
- The dock's panel close: `DashboardStripDock` → `DashboardPanelModel::removePanelAt`.
- The REST routes I built that mutate selection: `/dashboard/metric`, `/dashboard/metric/remove`,
  `/dashboard/metric/mode` in `src/app/RestServer.cpp` — they currently replicate the dialog's two-step.
- `DashboardSignalModel` (now the source of truth, carries availability/renderability roles after 2.3) and
  `DashboardPanelModel` (placement) and `model::DisplayRefsForSignal`.

## Deliverable — `DashboardSelectionController`
1. **New class** `DashboardSelectionController : public QObject` (place it per the architecture reference —
   `src/app` or `src/model`, your judgment; it holds non-owning `QPointer`s to both models and is owned by
   `ReaderMainWindow`). It is the single owner of compound selected-metric mutations. Verbs:
   - `QUuid addMetric(const SignalDescriptor&, const SignalAnchor&, const QStringList& modes, panelId/target)`
     — does addSignal → `DisplayRefsForSignal` → `addDisplayRef` **in one call**, preserving the existing
     order so the controller's rebuild still fires (signalAdded then displayRefsChanged). Returns the new id.
   - `bool removeMetric(const QUuid& signalId)` — removes the signal; cascades ref removal.
   - `bool setMetricMode(const QUuid& signalId, const QString& mode, bool enabled, panelId/target)` — mutates
     the signal's modes AND the placement refs together (the dialog's `onActiveModeToggled` body).
   - `int removePanel(const QUuid& panelId, PanelRemovalPolicy policy)` with
     `enum class PanelRemovalPolicy { RemoveReferencedMetrics }` (the only policy for now) — removes the panel
     and the metrics it referenced, via the cascade.
   - `void clearAllMetrics()` — empties the selected list + refs (used by the next step's run-load reset).
2. **Absorb the anonymous coordinator.** Delete `DashboardSignalPanelCoordinator` from `ReaderMainWindow.cpp`;
   move its `signalRemoved`/`displayRefRemoved` cascade + the `signalsBeingRemoved_` guard INTO the controller
   (so the cascade still happens for any removal path). Instantiate the controller where the coordinator was.
3. **Route every caller through the controller — no caller writes both models directly anymore:**
   - `SignalDisplayDialog`: `onAddSelected → controller->addMetric(...)`; `onRemoveActive → removeMetric`;
     `onActiveModeToggled → setMetricMode`. Inject the controller pointer into the dialog the same way the
     models are injected.
   - `DashboardStripDock` close-panel → `controller->removePanel(panelId, RemoveReferencedMetrics)`.
   - `RestServer` `/dashboard/metric*` → call the controller verbs (it becomes the controller's second caller,
     which is the point — one door for UI and REST).
4. **A count signal for the next step:** emit `selectedCountChanged(int)` (and/or `selectionChanged()`) from
   the controller whenever the selected count changes. (2.5b wires dock visibility to it; don't implement dock
   visibility here.)
5. **Instrument trace:** emit a `h5reader.dashboard` (`cDash`) line per verb
   (`event=selection_add/remove/set_mode/remove_panel/clear id=... count=...`) so the controller's actions are
   greppable alongside the existing model-level lines.

## Preserve exactly (this is the whole risk)
- Same add/remove/mode/close-panel **outcomes** as today; the `DashboardDisplayController` rebuild must still
  fire on the same signals (verify the rebuild triggers survive the routing change).
- The removal cascade semantics + reentrancy guard — byte-for-byte behavior, just relocated.
- The 2.2 capability-table asymmetries and 2.3 availability roles — untouched.
- Do NOT add dock-visibility behavior, run-load reset, or any UX change here. Pure structural move.

## Output
Implement, build green (`cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target
h5reader -j$(nproc)`). **Do NOT launch, do NOT git.** Append a `notes/INSTRUMENT_DASHBOARD_2026-06-05.md`
section (the controller's API + which callers now route through it + the deleted coordinator), and note
anything in the code that contradicted this brief (code wins — say what you found).
