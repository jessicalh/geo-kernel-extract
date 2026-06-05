# Desk-ready punch list — 2026-06-05 (lead, mouse-driven)

The selected-metrics coherence rebuild (2.1–2.5b: capability table, availability roles, single-source-of-truth
`DashboardSelectionController`, reveal-on-add, run-load reset) is DONE and verified at the plumbing level
(add → controller → panel; opus-clean cascade). **But that was plumbing.** Driving by mouse, the lead surfaced
the layer that actually blocks handing this to an advisor. None of it is where selection-state lives; it's
what the UI *means and shows*. This is the real desk-ready list.

> Verification caveat: a REST-driven run on `:1` got tangled with the lead's mouse (same window) — discarded,
> no conclusions. Going forward: isolate the test display or hand off explicitly so my drives never collide.

## Week sequencing (lead, 2026-06-05) — the controlling order
1. **Make the strips work** — DONE + committed `141e7b6` (B2). Lead mouse-confirming on `:1`.
2. **Camera + "Kabsch gyroscope" + hard backbone modes** — the 3D/stabilization pass. Maps to A3 (Focus/Free)
   + D1 (reveal steals camera) + the transform/Kabsch behaviour. "Gyroscope" = the fit may leave a residual
   spin/precession (investigate before touching). "Hard backbone modes" = a stricter backbone lock. NEXT.
3. **Remove cruft + things that almost work** — A1/A2 (toolbar chrome, meaningless buttons), B1 (hide
   render-nothing modes), the dead `SignalBuffer::isValidAt`, the `static.tensor` deferred glyph, etc.
4. **Back to fields + display modes** — C1 (field glossary) + the display-mode semantics. LAST.

## A — Chrome wrong / meaningless
- **A1. Playback bar and button bar not separated.** (A round-1 "split into Playback + Tools rows" supposedly
  landed — verify why the lead still sees them merged: regression? not visually distinct? built binary stale?)
- **A2. Meaningless buttons along the top** ("the biggies"). Inventory every toolbar action; cut/rename cruft
  (candidates: `Instrument` harness preset, duplicate Free, camera-mode buttons with no effect).
- **A3. Focus vs Free camera — no visible difference.** What is Focus supposed to do; why identical to Free.

## B — Metric display is semantically hollow (the core)
- **B1. Display modes have no display impact.** strip / glyph / spectrum / table — prominent checkboxes, but
  picking different ones changes nothing. Most modes don't render (capability table: most are
  `buildsPanelWidget=false`/`hasVisibleSurface=false`). Fix direction: **don't OFFER modes that render nothing**
  (hide, not grey — the lead's "don't show options with nothing in them"); keep only modes with a real surface.
- **B2. A "strip" shows one dot for the whole run.** Either a static field mis-cast as a time series, or the
  per-frame sampling emits one point (broken). [Investigation Q1.] Verdict pending: sampling bug vs
  playback-tied sampling vs static-field-shouldn't-be-a-strip.
- **B3. Fields give no logical response and can't be told apart by meaning.** (Overlaps B1/C1.)

## C — Can't tell what you're selecting
- **C1. No field identity / no dictionary.** Only the last column names the field; nothing explains what it
  *means*. Need a per-metric glossary (descriptor → human meaning + units) surfaced at selection time.
  (Format/approach is a "discuss first" — design the glossary source before building.)

## D — Selecting disrupts the 3D view
- **D1. Selection kicks into "video modes" that lose the stabilization.** The fit-reference/backbone
  stabilization (TransformedConformation) gets thrown away → view looks random. [Investigation Q2.] Likely a
  coupling where adding a strip / "Follow current frame" starts playback or bypasses `atomPosition()` so the
  molecule tumbles. Verdict pending.

## Verdicts (opus investigation 2026-06-05, code-only)

- **B2 — strip one dot. RESOLVED 2026-06-05** (continuity test, H5-vs-NPY). Root cause was a THIRD, deeper bug
  than the two below: `TransformedConformation::loadSnapshot` (`:73`) returned `inner_->snapshot(frame)` — the
  resident *accessor* (non-null only if that frame was already loaded) — instead of *triggering* the inner's
  synchronous load. The display renders from positions (not snapshots), so the inner only ever held the
  startup frame → every NPY strip sample read that stale slot (1/16). H5 strips (in-memory dense buffers) were
  unaffected (16/16) — which is how H5-vs-NPY isolated it. FIX: `loadSnapshot` now calls
  `inner_->requestSnapshot(frame)` before forwarding. PLUS strip history now PERSISTS across `rebuild()` (was
  `series_ = std::move(next)` with empty buffers, discarding history on every add/remove/focus change; now
  buffers carry over by (signal,channel,mode)). Verified: NPY z-series matches the NPY ground truth 16/16.
  The (a)/(c) below were real secondary observations but not the cause of the empty strips.
- **B2 (original verdict) — strip one dot. TWO compounding causes.**
  - (a) *By design, all backends*: strip series fill **lazily as playback advances** — `extendToFrame` samples
    only `[lastFrame()+1 .. current]`; "the PAST stays, the FUTURE is never written" (`StripChartChannel.h:16`).
    A reader **paused at frame 0 shows exactly one point**, for any field. Fills as you play/scrub. → DESIGN
    CHOICE: pre-fill `[0..N-1]` on add for the *cheap dense-H5* sources (in-memory; the scrub-deferral guard
    exists for expensive snapshot/NPY sources, not these). Trades "honest no-future strip" for a filled chart.
  - (c) *Wrong display model for this field*: `h5:dssp8_transition` is a per-RUN aggregate **transition matrix**
    — `ss8_transition_count` shape `(R,)`, one cumulative count per residue; its sampler **ignores `frame`**
    (`DashboardDisplayController.cpp:988`), so every frame's value is identical → flat line even fully sampled.
    It should render as a **matrix/heatmap panel, never a strip**. The genuine per-frame DSSP field is
    `h5:dssp8_time_series`. → FIX: drop `strip.*` from transition-matrix/EventRecord descriptors; remove
    `dssp8_transition`/`dihedral_bin_transition` from `kDensePaths`. (Lead's NPY hypothesis: directionally
    right that per-frame data is present-but-unseen for *temporal* fields, BUT this field isn't NPY and isn't
    temporal — it's an H5 aggregate mis-cast. NPY sources are deferred-for-cost = the "fancy scheme," not the
    culprit here.)
- **D1 — selection disrupts view. NOT adding a metric; molecule stabilization NEVER lost.** The
  `TransformedConformation` Kabsch fit is always-on through playback (molecule never tumbles). What re-aims the
  CAMERA is a deliberate left-click on a strip's **"Reveal" hotspot** → `focusCameraOnReveal` → a per-frame
  `Subset/Atom/Dihedral` camera lock that re-aims every frame ("video"). Easy to trigger by a stray click. →
  FIX: decouple reveal-highlight from camera re-aim; gate the `focusCameraOnReveal` (`MoleculeScene.cpp:533`)
  behind a modifier, or make it one-shot (aim once, drop back to Free). UX CHOICE.
- **A3 — Focus vs Free.** Genuinely different modes; difference only visible **in motion**. Focus inherits the
  current sight + distance, so paused it looks identical to Free; diverges on playback (Focus keeps the
  residue backbone plane centered/face-on). → FIX: snap Focus to a canonical view (don't inherit sight,
  `CameraComposer.cpp:190-199`) so it visibly reorients; + show a "locked" toolbar indicator
  (`updateCameraModeActions` leaves Atom/Bond/Subset unchecked, `:1061-1066`).
- **A1/A2/Q4 — toolbar.** Already TWO separate toolbars on separate rows (`PlaybackToolbar` + `ToolsToolbar`,
  `addToolBarBreak` at `:827`). "Not separated" = **visual perception** (no gap; Tools bar opens with a leading
  separator that runs into Playback). → FIX (cosmetic): spacer/heading, or drop the leading separator (`:834`).
  The "meaningless buttons" = the **camera-mode group** (Focus/Newman/Plane lock/Free): 3 of 4 disabled most of
  the time (need 1-focus / 4-sel / 3-sel), the enabled one (Focus) invisible paused. → FIX: hide-until-usable
  or fold into a single "Camera" popup; remove the (already-hidden) `instrumentAction_`. No duplicate Free, no
  visible Instrument button.
- **B1 — modes no impact.** Most display modes have `buildsPanelWidget=false`/`hasVisibleSurface=false`
  (capability table) → render nothing. → FIX: stop OFFERING render-nothing modes (hide, not grey).
- **C1 — glossary.** Still open; needs a field-meaning source (descriptor → human meaning + units). Discuss
  the source before building.

**One real "wrong display model" flag:** B2(c) only. Everything else is small-fix / config / cosmetic / UX-choice.

## Testing infra
- Isolated **Xvfb `:99`** is up (software GL) — my REST drives + screenshots run there, never on the lead's
  `:1`. Real-path verification (the `/dashboard/picker/add` → genuine `onAddSelected` probe) + play-and-watch
  per backend happen here. **Nothing is "done" until the lead confirms by mouse on `:1`.**

## Status
- "All of these things are addressable." (lead) — week-long deliverability push; advisor-shareable is the bar.
- Verdicts in hand. Choices owed to the lead: B2(a) strip pre-fill (honest vs filled), D1 reveal camera-grab
  (opt-in/one-shot), C1 glossary source. Quick wins ready: A1 (cosmetic bar gap), B1 (hide dead modes),
  A2 (camera-group cleanup). Deeper: B2(c) (matrix panel + catalog fix), A3 (Focus snap), D1.
