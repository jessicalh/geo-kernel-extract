# Codex — STEP 2.2: one shared display-mode capability table (behavior-preserving)

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`). **Portable Qt6
(Windows/macOS/Linux) — cross-platform source is first-class.** Scope: `h5-reader/src` ONLY. Never link the
`nmr_shielding` library, never write H5.

**Verify from the code before changing anything, ever.** Cites are against the current working tree (commit
`8007e2e`); if one is stale, follow the code and tell me.

## Why this matters
Three separate predicates each classify a display mode, and they disagree. They're the seam the whole
selected-metrics rebuild rests on — so before the single-source-of-truth controller lands, we collapse them
into ONE shared table that is the single answer to "what can this display mode do." **This is a pure
refactor: behavior must be byte-identical afterward.** No mode becomes newly pickable or newly greyed. We're
removing duplication and a drift hazard, not changing what the user sees.

## The three predicates today (read all three, confirm the truth table)
- `modeHasVisibleSurface(m)` — dialog "can the user pick it" gate (`src/app/SignalDisplayDialog.cpp`, the
  global `DashboardModeHasVisibleSurface`).
- `isPanelMode(m)` — controller "build an `AbstractStripPanel` widget via setOwnedPanels"
  (`src/app/DashboardDisplayController.cpp:63`).
- `isPanelDisplayMode(m)` — panel-model "emit a `panel` sentinel ref" (`src/model/DashboardPanelModel.cpp:61`).

The EXACT current truth (reproduce it precisely — these two disagreements are deliberate and must survive):

| mode                        | hasVisibleSurface | buildsPanelWidget | emitsPanelRef |
|-----------------------------|:---:|:---:|:---:|
| `strip.*` (prefix)          |  T  |  F  |  F  |
| `static.bar.sequence`       |  T  |  T  |  T  |
| `static.spectrum.power`     | **F** |  T  |  T  |  ← renderable panel the dialog greys (PRESERVE as-is)
| `static.curve.lag.animated` |  T  |  T  |  T  |
| `static.chord.coupling`     |  T  |  T  |  T  |
| `static.fixed_freq`         |  T  |  T  |  T  |
| `static.tensor`             | **F** | **F** |  T  |  ← scene-overlay; ref-only (PRESERVE the asymmetry)
| anything else               |  F  |  F  |  F  |

(`static.spectrum.power` staying greyed in the dialog is INTENTIONAL for this step — whether to make it
pickable is a separate deliberate call the lead will make later. Do not change it here.)

## Deliverable
1. **New shared header** `src/model/DisplayModeCapability.h` (model-adjacent, not in the dialog):
   ```cpp
   namespace h5reader::model {
   struct DisplayModeCapability {
       bool hasVisibleSurface = false;  // dialog: offer as a pickable display mode
       bool buildsPanelWidget = false;  // controller: build an AbstractStripPanel via setOwnedPanels
       bool emitsPanelRef     = false;  // panel-model: track with a "panel" sentinel display ref
   };
   DisplayModeCapability DisplayModeCapabilityFor(const QString& mode);
   }
   ```
   Implement `DisplayModeCapabilityFor` to return exactly the truth table above (a single `.cpp` — e.g.
   `src/model/DisplayModeCapability.cpp` — added to the model sources in `CMakeLists.txt`, OR header-only
   `inline` if you prefer no new TU; your call, keep it clean).
2. **Replace the three predicates with table lookups** at all call sites:
   - dialog `DashboardModeHasVisibleSurface(m)` → `DisplayModeCapabilityFor(m).hasVisibleSurface` (keep the
     existing free function as a one-line wrapper if other code/REST calls it, or repoint callers — your
     judgment, but no behavior change).
   - controller `isPanelMode(m)` → `…buildsPanelWidget` (`DashboardDisplayController.cpp`).
   - panel-model `isPanelDisplayMode(m)` → `…emitsPanelRef` (`DashboardPanelModel.cpp`).
   Delete the now-duplicated literal mode lists. Keep the comments that explain WHY tensor/spectrum disagree —
   move them onto the table.
3. The STEP-1 instrument already reports these three flags per mode in `/dashboard/state`; repoint it (and the
   `/dashboard/picker` reason logic) at the shared table too, so there's exactly one source.

## Behavior-preserving proof (state it in your notes)
The point is zero behavior change. Confirm by reasoning that every mode's three booleans are unchanged from
the table above. (I will verify at runtime by diffing `/dashboard/state` per-mode flags and the
`/dashboard/picker` enabled/reason output before vs after — they must be identical.)

## Discipline
`CENSUS_REGISTER`/`ACONNECT`/`ASSERT_THREAD` where applicable (this is mostly free functions, so likely
none); no new `QTimer`; portable; match existing idioms. If adding a `.cpp`, wire it into the h5-reader
`CMakeLists.txt` model sources.

## Output
Implement, build green (`cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target
h5reader -j$(nproc)`). **Do NOT launch, do NOT git.** Append a short section to
`notes/INSTRUMENT_DASHBOARD_2026-06-05.md` (the new header + the call sites repointed), and note anything in
the code that contradicted this brief (code wins).
