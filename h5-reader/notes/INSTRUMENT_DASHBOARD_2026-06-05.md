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
