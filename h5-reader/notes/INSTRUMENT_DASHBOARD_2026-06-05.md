# Dashboard Selection Instrumentation - 2026-06-05

## Routes

All examples assume the app was started with `--rest` and the scraped port is in `$PORT`.

```bash
curl -s "http://127.0.0.1:$PORT/dashboard/state" | jq .
```

Returns selected metrics plus dock/render state:

```json
{
  "selected": [
    {
      "id": "uuid",
      "descriptor_id": "h5:kernel_dynamics_psd",
      "concept_key": "kernel_dynamics.psd",
      "label": "Kernel power spectrum (PSD) - Atom 0",
      "anchor": {"kind": "atom", "atom": 0},
      "enabled": true,
      "availability": "Available",
      "modes": [
        {
          "mode": "static.spectrum.power",
          "enabled": true,
          "is_panel_ref": true,
          "renderable_panel": true,
          "has_visible_surface": false
        }
      ]
    }
  ],
  "selected_count": 1,
  "dock": {"visible": true, "width": 360, "raised": true},
  "render": {"owned_panel_count": 1, "strip_track_count": 0}
}
```

```bash
curl -s -X POST "http://127.0.0.1:$PORT/dashboard/metric" \
  -H 'content-type: application/json' \
  -d '{"descriptor_id":"h5:kernel_dynamics_psd","anchor":{"atom":0,"frame":0},"modes":["static.spectrum.power"]}' | jq .
```

```bash
curl -s -X POST "http://127.0.0.1:$PORT/dashboard/metric" \
  -H 'content-type: application/json' \
  -d '{"descriptor_id":"h5:kernel_dynamics_psd","anchor":{"follows_focus":true},"modes":["static.spectrum.power"]}' | jq .
```

```bash
curl -s -X POST "http://127.0.0.1:$PORT/dashboard/metric/mode" \
  -H 'content-type: application/json' \
  -d '{"id":"<uuid>","mode":"static.spectrum.power","enabled":false}' | jq .
```

```bash
curl -s -X POST "http://127.0.0.1:$PORT/dashboard/metric/remove" \
  -H 'content-type: application/json' \
  -d '{"id":"<uuid>"}' | jq .
```

```bash
curl -s -X POST "http://127.0.0.1:$PORT/dashboard/dock" \
  -H 'content-type: application/json' \
  -d '{"visible":true}' | jq .
```

## Log Tail

Category: `h5reader.dashboard`. The logs are emitted from the model/controller/window transition sites, not from REST handlers, so UI-driven changes produce the same trace.

Example chain:

```text
event=signal_added id=<uuid> descriptor_id=h5:kernel_dynamics_psd concept_key=kernel_dynamics.psd modes=[static.spectrum.power] anchor=atom:0 selected=[<uuid>]
event=controller_rebuild owned_panel_count=0 strip_track_count=0 active_panel=<panel_uuid>
event=panel_refs_changed signal_id=<uuid> ref_count=1 panel_id=<panel_uuid> mode=static.spectrum.power channel=panel
event=controller_rebuild owned_panel_count=1 strip_track_count=0 active_panel=<panel_uuid>
event=dock_visibility_changed visible=1 width=360
```

Other transition lines:

```text
event=display_mode_toggled id=<uuid> mode=static.spectrum.power enabled=0
event=signal_removed id=<uuid> selected=[]
```

## Code Findings

- Availability was reachable cleanly through `TrajectorySignalCatalog::fieldAvailability()`, so `/dashboard/state` reports the catalog availability state and only returns `"unknown"` if no availability object is wired.
- The brief's `anchor.frame` field has no corresponding member in the current `model::SignalAnchor`; REST validates the optional frame range but does not store it.
- The three display-mode predicates intentionally disagree and are now exposed from their original homes. For example, `static.spectrum.power` is panel-ref/renderable but not dialog-visible; `static.tensor` is panel-ref only.
- `ReaderMainWindow::revealDockQueued()` is private, so REST uses a narrow public `setDashboardDockVisible()` wrapper that calls the same queued reveal path for show and `setVisible(false)` for hide.

Build verification:

```bash
cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target h5reader -j$(nproc)
```

Result: succeeded.

## Metric Picker Probe

The Metric Picker selector is observable without a human at the screen:

```bash
curl -s "http://127.0.0.1:$PORT/dashboard/picker" | jq .
```

```bash
curl -s -X POST "http://127.0.0.1:$PORT/dashboard/picker/open" \
  -H 'content-type: application/json' \
  -d '{"atom":0}' | jq .
```

`POST /dashboard/picker/open` accepts an optional `atom` integer. When supplied, the route applies the same selection pick path used by `/selection/pick`, then opens through the Metrics action gate; if `AtomSelection` still has no focus, it returns the gate reason instead of forcing the dialog open. The successful response and `GET /dashboard/picker` both return:

```json
{ "open": true, "candidate_row": 0, "anchor": {}, "modes": [] }
```

Selector reflection also emits one `h5reader.dashboard` line whenever the candidate selector refreshes:

```bash
tail -f "$LOG" | rg 'event=picker_selector_state'
```

```text
event=picker_selector_state current_row=0 enabled=[strip] disabled=[spectrum:This descriptor does not offer that display mode.,...]
```

## STEP 2.2 Shared Display Mode Capability Table

- Added header-only `src/model/DisplayModeCapability.h`; `DisplayModeCapabilityFor()` is now the single truth table for picker visibility, static panel construction, and panel-ref emission.
- Behavior proof from the current predicates before editing: `strip.*` remains visible-only; `static.bar.sequence`, `static.curve.lag.animated`, `static.chord.coupling`, and `static.fixed_freq` remain true for all three flags; `static.spectrum.power` remains dialog-greyed but panel-renderable/ref-emitting; `static.tensor` remains scene-overlay ref-only; unknown modes remain false for all flags.
- Repointed the dialog picker and active-mode gates, plus the `DashboardModeHasVisibleSurface()` compatibility wrapper, to `hasVisibleSurface`.
- Repointed `DashboardDisplayController`'s static panel build gate, plus `IsRenderableDashboardPanelMode()`, to `buildsPanelWidget`.
- Repointed `DashboardPanelModel`'s sentinel panel-ref gate, plus `IsPanelDisplayMode()`, to `emitsPanelRef`.
- Repointed `/dashboard/state` mode flags to one `DisplayModeCapabilityFor()` result per mode; `/dashboard/picker` reason logic keeps the same strings, now fed by the shared visible-surface flag.
- Code contradicted the brief: none found.

## STEP 2.3 DashboardSignalModel Availability and Renderability Roles

- Added `DashboardSignalModel` row availability state as first-class model state. New roles: `availabilityState`, `availability`, and `availabilityReason`; new accessors: `setFieldAvailability()`, `refreshAvailability()`, `availabilityAt()`, `availabilityName()`, `availabilityReason()`, and `isVisibleAvailable()`.
- Wired `ReaderMainWindow` to inject the same `shared_ptr<const TrajectoryFieldAvailability>` into `DashboardSignalModel` that already feeds `TrajectorySignalCatalog` and `QtAtomInspectorDock`.
- Added per-selected-row renderability from the STEP 2.2 table without duplicating the table. New role/accessors: `modeRenderability`, `renderableModeCount()`, `ModeRenderabilityFor()`, and count accessors `selectedCount()`, `renderableSelectedCount()`, `unavailableCount()`, `noRendererCount()`.
- Repointed `/dashboard/state` per-metric `availability` to `DashboardSignalModel::availabilityName(row)`, per-mode flags to `DashboardSignalModel::ModeRenderabilityFor(mode)`, and `selected_count` to `DashboardSignalModel::selectedCount()`. The JSON field names stay unchanged.
- Fixed the controller rebuild path so unavailable selected rows are skipped for rendering via `DashboardSignalModel::isVisibleAvailable(row)` but remain selected rows with an availability reason. A descriptor that is still unresolved after that is logged as `event=signal_descriptor_unresolved ... action=skip_render` instead of disappearing silently.
- Behavior check: healthy-run `/dashboard/state` values should remain unchanged for selected metrics and mode flags; the intended new observable is that a selected metric filtered out by availability stays listed and reports a non-visible availability state/reason instead of vanishing from render rebuilds.
- Build verification: `cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target h5reader -j$(nproc)` succeeded.
- Code contradicted the brief: `/dashboard/state` was already listing selected rows from `DashboardSignalModel`; the silent drop was specifically in `DashboardDisplayController::rebuild()` render construction, not in the REST selected-array enumeration.

## STEP 2.5a DashboardSelectionController Mutation Owner

- Added `src/app/DashboardSelectionController.{h,cpp}` as the named owner for compound selected-metric mutations. It is a `QObject`, calls `CENSUS_REGISTER(this)`, stores non-owning `QPointer`s to `DashboardSignalModel` and `DashboardPanelModel`, and is owned by `ReaderMainWindow`.
- Public mutation API:
  - `addMetric(descriptor, anchor, modes, PanelTarget, followsFocus, label, addedRefs)` preserves the existing order: `DashboardSignalModel::addSignal()` first, then `DisplayRefsForSignal()`, then `DashboardPanelModel::addDisplayRef()`.
  - `removeMetric(signalId)` removes the selected signal and lets the signal-removal cascade remove refs.
  - `setMetricMode(signalId, mode, enabled[, PanelTarget])` preserves the existing mode order: mutate signal modes first, then add/remove panel refs.
  - `removePanel(panelId, PanelRemovalPolicy::RemoveReferencedMetrics)` removes the panel and prunes only metrics whose ref count reaches zero through the display-ref cascade.
  - `clearAllMetrics()` removes panel refs and selected signals for the next run-load reset step.
- Deleted the anonymous `DashboardSignalPanelCoordinator` from `ReaderMainWindow.cpp`. Its `signalRemoved -> removeDisplayRefsForSignal()` and `displayRefRemoved -> removeSignal()` cascade moved into `DashboardSelectionController`, including the `signalsBeingRemoved_` reentrancy guard.
- Routed callers through the controller:
  - `SignalDisplayDialog::onAddSelected`, `onRemoveActive`, and `onActiveModeToggled`.
  - `DashboardStripDock::onPanelTabCloseRequested`.
  - REST routes `POST /dashboard/metric`, `POST /dashboard/metric/remove`, and `POST /dashboard/metric/mode`.
- `DashboardDisplayController` rebuild triggers are preserved. It still rebuilds on `signalAdded` before refs exist and again on `displayRefsChanged` after refs are added; the comment now names `DashboardSelectionController` as the two-step source.
- New controller count signals: `selectedCountChanged(int)` and `selectionChanged()`, emitted only when `DashboardSignalModel::selectedCount()` changes.
- New `h5reader.dashboard` verb traces:
  - `event=selection_add id=<uuid> count=<n> panel_id=<uuid> added_refs=<n> modes=[...]`
  - `event=selection_remove id=<uuid> count=<n> removed=<0|1>`
  - `event=selection_set_mode id=<uuid> count=<n> mode=<mode> enabled=<0|1> refs_changed=<n>`
  - `event=selection_remove_panel id=<panel_uuid> count=<n> removed_metrics=<n> removed=<0|1>`
  - `event=selection_clear id=all count=<n> refs_removed=<n>`
- Build verification: `cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target h5reader -j$(nproc)` succeeded.
- Code contradicted the brief: no behavioral contradiction found. Implementation detail from the code: `setMetricMode()` needs catalog access because `DashboardSignalModel` rows store descriptor ids and modes, but not the channel descriptors required by `DisplayRefsForSignal()`.

## STEP 2.5b Dock Reveal-on-Add and Run-Load Reset

- `ReaderMainWindow` now consumes `DashboardSelectionController::selectedCountChanged(int)`.
  It tracks the previous selected count and only reveals the dashboard strip when the count increases and
  `dashboardStripDock_` is hidden. The reveal path is still `revealDockQueued(dashboardStripDock_)`, preserving
  the existing queued resize/raise behavior. New log line:
  `event=dock_reveal_on_add count=<n>`.
- Run-load clean slate is enforced after `restoreAllSettings()` in `ReaderMainWindow::resetDashboardStateForRunLoad()`.
  That keeps window geometry/dock layout restoration separate, then calls `DashboardSelectionController::clearAllMetrics()`,
  hides `dashboardStripDock_`, and resets the local selected-count baseline. New log line:
  `event=selection_reset_on_load count=0 dock_visible=0`.
- Added `POST /dashboard/picker/add`. It calls the dialog's real `SignalDisplayDialog::onAddSelected()` slot,
  the same slot wired to the Add Signal button, using the current picker row, anchor row, checked mode boxes,
  and panel target widgets. The response returns the resulting picker state, `selected_count`, and dashboard
  dock `{visible,width,raised}`.
- Source check for the nondeterministic `h5:dssp8_transition` seed: current `src/` has no dashboard-selection
  restore from `QSettings`, no `.LGS`/session restore, and no `ReaderMainWindow` hard-coded add. The only
  current `h5:dssp8_transition` occurrence in product code is descriptor registration in
  `TrajectorySignalCatalog.cpp`; it is catalog-visible, not selected. The old default selected-metric source
  appears to be retired code documented in stale REST tests/notes (`npy:dssp_chi`), not present in the current
  `src/` tree. Because no current source exists to suppress directly, the reset is the post-load defensive
  clear/hide described above.
- Added `GET /dashboard/strip/series` to dump active strip `ChannelBuffer` values/valid arrays; source check:
  `stripTracks()` is reconstructed from controller `series_`, not stored as a separate live list.
- 2026-06-05: Added `snapshotReady(frame)`-driven strip backfill for frame-NPY snapshot series. `SignalBuffer::backfill()` overwrites placeholder gaps in place, preserving contiguous history and valid-only y-range growth; `DashboardDisplayController` now connects `Conformation::snapshotReady` with `ACONNECT` and emits `stripTracksChanged()` after a backfill. Build verification: `cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target h5reader -j$(nproc)` succeeded.
- 2026-06-05: Reverted the moot `DashboardDisplayController` snapshot-backfill path (`onSnapshotReady`, the controller `snapshotReady` connection, `SignalBuffer::backfill()`, and `ChannelBuffer::set()`) and preserved strip history across `rebuild()` by moving old buffers into matching new series keyed by `(signal.id, channel.id, displayModeId)`; new metrics still start empty and `event=strip_sample` diagnostics remain. Build verification: `cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target h5reader -j$(nproc)` succeeded. Code check: in this checkout `activePanelChanged` calls `refreshPanelVisibility()`, not `rebuild()`.
