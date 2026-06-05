# Codex mission — execute the UI coherence surgery (the full cut, one pass)

Working root: `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`; do not switch). You
own the grind; Claude reviews the diff; the lead tests + owns ALL git. Honest stakes: this is the surgery that
turns the reader from accreted-ad-hoc into a coherent advisor-grade vetting tool. You already mapped it — now
execute your own prioritized list with the lead's calls folded in. ONE coherent cut: Claude rebuilds +
relaunches ONCE for the lead to test the whole thing.

## Your spec is your own review
`notes/UI_COHERENCE_REVIEW_2026-06-05.md` — its "Prioritized Surgery List" + per-concern fixes (cited
`file:line`) ARE the implementation spec. Also read `notes/UI_SURGERY_2026-06-05.md` (the lead's decisions).
The lead's organizing principle: **no buttons or fields for what we don't need.**

## The cut (priority order)

1. **Camera clipping on free-camera gestures (the zoom cutoff).** Your recommended fix: cache the padded
   per-frame bounds in `MoleculeScene` during `setFrame`; add `MoleculeScene::syncCameraClippingRange()` that
   calls `ResetCameraClippingRange(cachedPaddedBounds_)` when bounds are valid; call it after EVERY camera
   write — `CameraInputFilter` gesture (after `applyGesture`, before `requestRender`), `ResetCamera`, REST
   camera set/clear — and collapse the `lockCameraToSelectionPlane` special-case into the same helper.

2. **Loaded-run field availability (the BIG piece — get it right).** Build the shared availability layer you
   designed: `TrajectoryFieldAvailability` (near `model`/`io`, NOT widget code) built after `QtTrajectoryH5`
   load + frame-NPY discovery, classifying each field {Absent, NoFramePayload, AllMissing, AllZeroStructural,
   AllZeroObserved, Available}. Consume it in BOTH:
   - **inspector** — group adders return `bool`; a group is created only if ≥1 real row was added; dash-only
     groups disappear.
   - **metric catalog / `SignalDisplayDialog`** — expose only descriptors/channels whose availability is not
     Absent / AllMissing / AllZeroStructural; `canSample()` consults availability.
   - **gate** the default `npy:dssp_chi` seed on availability (else remove it).
   ONE shared layer, not per-panel filters (else inspector/dialog/strip keep disagreeing). This is the lead's
   "show nothing if nothing exists; find the empty stuff on load."

3. **Pane recovery obvious.** Surface the existing dock `toggleViewAction`s in a proper **View → Panels** menu
   + a clearer toolbar grouping. Share the same `QAction`s between menu + toolbar so checked state stays right.

4. **Transform UI → ONE switch (lead's call).** Collapse the four-option transform menu to a single switch
   toggling **backbone-fit ↔ all-atom-fit**. REMOVE Identity + Center-COM entirely. **Default on load =
   backbone-fit** (`FitSubset` over the typed backbone subset, reference frame 0) — set it as the
   `TransformedConformation`'s initial mode so the reader opens stationary. The switch flips to all-atom
   `FitReference`. Keep the shared Kabsch + degeneracy policy. ("Center COM" is gone, so that mislabel
   resolves itself.)

5. **Cut / quarantine the accreted cruft (your "Cut or fold in" list).** Hide/rename the **"Instrument"**
   toolbar action (keep the REST hook); remove dormant **`QtSelectionOverlay`** (or fold any useful marker
   into `MeasurementOverlay`); move **`/proc/self/statm`** out of the render path; demote the per-frame
   **`qCInfo` overlay logs** to debug; fix the stale comments (`installRenderTimer`; the render-verb comment
   in `ReaderMainWindow`); fix **`QtFieldGridOverlay::setMode`** to recompute the CURRENT frame, not frame 0;
   store the typed mapper in **`QtBFieldStreamOverlay`** (drop the per-update `dynamic_cast`); settle the
   **PROBE** comment (document as the fix or remove). Camera checked-state (Atom/Bond/Subset with no checked
   action): add a visible "Locked/Follow" indicator if quick, else log it for a follow-up.

## qt6 discipline (non-negotiable)
`CENSUS_REGISTER` in any new QObject; `ACONNECT` for connections; `ASSERT_THREAD` atop VTK-mutating methods
(GUI-thread); NO new persistent `QTimer`; renders ONLY via `MoleculeScene::requestRender`; UDP-log meaningful
state changes. The availability layer is model code (plain types; QObject/signals only if it truly needs them).

## Boundaries
- **`h5-reader/src` only.** Do NOT touch the extractor (`/shared/2026Thesis/nmr-shielding/src`), the model's
  H5/calculator loading semantics, or the rediscover code. No new physics.
- **No git** — the lead owns ALL git. Edit + build.
- **Build to verify:** `cmake --build build/linux-rwdi --target h5reader -j$(nproc)` must succeed. **Do NOT
  launch the GUI** — Claude rebuilds + relaunches for the lead.
- Keep `TransformedConformation` as the display seam; don't bypass it. Model-is-spine.

## Report
Per surgery item: files changed + what + why (cited). The availability layer's design as built (the
classification + where it plugs in). Build result. How to test each (zoom no longer clips; opens stationary
on backbone-fit + the switch flips to all-atom; empty groups/channels gone; View→Panels works; Instrument +
cruft gone). Any risks or pieces left for a follow-up.
