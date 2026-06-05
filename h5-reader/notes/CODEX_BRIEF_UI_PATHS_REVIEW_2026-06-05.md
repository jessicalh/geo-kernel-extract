# Codex — comprehensive UI interaction-paths review (READ-ONLY)

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch h5-reader-pysr-spike). **READ-ONLY: map +
diagnose. Do NOT edit, build, git, or launch.** The lead has now twice seen panels fail to appear — "I select
a field, I add it, and alas no panel" — and a targeted fix did not take. Stop patching blind: map the UI
**interaction paths** end-to-end (user action → signal/slot → model → view → visible result), cite `file:line`,
and find where each one breaks, tangles, or dead-ends. The lead's principle: the reader is the advisor-facing
vetting surface; a path that produces no visible result is broken.

## Paths to map (priority order)

1. **PRIORITY — Metric → panel → visible dock.** The live bug: add a metric, see no panel. Trace the FULL
   chain: `SignalDisplayDialog` (select a descriptor → the Add action) → `DashboardSignalModel::addSignal` →
   `DashboardPanelModel` display refs → `DashboardStripDock` + `DashboardDisplayController` building the strip
   → the **`DashboardStripDock` visibility**. The Strip dock now starts HIDDEN (the docks-hidden-on-startup
   change) and the saved layout was discarded (kSettingsVersion bump). Questions to answer concretely:
   - When a signal is added (or when the Metrics dialog opens / commits), does ANYTHING show or raise the
     Strip dock? If not, that is the bug — the panel renders into a hidden dock.
   - Does the recent resize-to-360 + raise-on-show wiring actually make a toggled dock appear, or does the
     all-tabified-then-all-hidden state (a known Qt quirk) defeat `setVisible(true)`?
   - Is the panel even being created/populated (DashboardPanelModel / DisplayController), or does the add
     short-circuit? Pin EXACTLY why the user sees nothing after Add, and the concrete fix.

2. **Dock / panel visibility.** The three docks (Inspector/Selection/Strip) start hidden + tabified together;
   View→Panels menu + the Tools-toolbar Panels button (shared `toggleViewAction`s); the resize/raise-on-show.
   Trace the tabify + hidden + toggle interaction; does toggling actually yield a usable-width visible dock?

3. **Pick → selection → inspector.** Double-click → `QtAtomPicker` → `AtomSelection` → inspector dock (hidden)
   + measurement overlay. Does a pick surface anything the user can see while the inspector dock is hidden?

4. **Camera.** Gesture → `CameraInputFilter` → `CameraComposer` → render + the new clipping sync.

5. **Transform.** Backbone-fit default → `TransformedConformation::atomPosition` → scene / overlays / picker.

6. **Overlay toggles.** Ribbon / Rings / Butterfly / B-field → scene overlays → render.

## Output
Write `notes/UI_PATHS_REVIEW_2026-06-05.md`: per path, a short flow map (the chain, cited) + the break/tangle
+ the concrete fix. **Lead with the metric→panel→dock root cause + fix.** Then a prioritized fix list.
READ-ONLY — diagnose only; we implement from the map.

## Boundaries
No edit / build / git / launch. `h5-reader/src` scope. Cite `file:line`. Be concrete — the actual break, not
"it's complex." If a path can't be settled from code alone, say what a runtime check would resolve it.
