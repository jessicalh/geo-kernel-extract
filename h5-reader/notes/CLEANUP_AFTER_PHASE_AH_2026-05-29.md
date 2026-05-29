# Cleanup punch list after Phases A-H (2026-05-29)

Phases A-H landed the 5 new TR surfaces (iRED, KernelDynamics,
ReorientationalDynamics, DihedralAutocorrelation, KernelCoherence) end-to-end
into h5-reader. All three ctest targets are green:

```
h5reader_model_tests    — catalog + anchor + cascade + lockstep regression
h5reader_scene_math_tests
h5reader_rest_smoke
```

The implementation passed adversarial review by an Anthropic agent + a Codex
agent (sandbox-bypassed). The reviews surfaced real bugs that were tagged
NOW (must-fix-before-Phase-D) vs LATER (Phase H polish). NOW items #1, #3,
#5 (from the first round) were folded into Phases D/E/F/G as templates so
the clones didn't carry the defects. The items below are what remained
after that inline-application — they were deliberately deferred per the
user's "wrap through phases completely before going into long fix mode"
directive.

Plan file: `/home/jessica/.claude/plans/validated-questing-metcalfe.md`
Scope doc: `h5-reader/notes/SCOPE_NEW_TRS_2026-05-29.md`

## NOW — must fix this session (or before in-app use)

These are bugs that prevent the new descriptors from actually working in
the running app even though tests pass.

### NOW-1: `takeOwnedPanels` drain-on-tick

**Where:** `h5-reader/src/app/DashboardDisplayController.cpp` (the
`takeOwnedPanels()` method moves `ownedPanels_` OUT) +
`h5-reader/src/app/DashboardStripDock.cpp:277` (calls
`controller_->takeOwnedPanels()` from `refreshTracks()` which is
connected to `stripTracksChanged()`) +
`h5-reader/src/app/DashboardDisplayController.cpp:1205` (`setFrame()`
emits `stripTracksChanged()` on every frame tick).

**Symptom:** First frame change after adding any new-TR signal blanks
the owned panels from the strip widget. iRED/KernelDynamics/etc.
display flashes once then vanishes.

**Fix:** Separate signal `ownedPanelsChanged` emitted ONLY from
`rebuild()` (not from `setFrame()`). Dock connects to it independently
of `stripTracksChanged`. `takeOwnedPanels()` stays as-is.

```cpp
// In controller .h:
signals:
    void stripTracksChanged();       // existing
    void ownedPanelsChanged();       // new

// In rebuild():
    ...
    ownedPanels_ = std::move(nextPanels);
    ...
    emit stripTracksChanged();
    emit ownedPanelsChanged();       // new

// In setFrame() / refreshPanelVisibility():
    // emit stripTracksChanged() — keep as-is, owned panels NOT touched

// In dock setSignalModels():
    ACONNECT(controller_, &DashboardDisplayController::ownedPanelsChanged,
             this, [this]() { stackWidget_->setOwnedPanels(controller_->takeOwnedPanels()); });
```

### NOW-2: SignalDisplayDialog doesn't know about static-panel modes

**Where:** `h5-reader/src/app/SignalDisplayDialog.cpp` — search for
mode-kind registration. The 4 new mode strings are not recognized:
`static.bar.sequence`, `static.spectrum.power`,
`static.curve.lag.animated`, `static.chord.coupling`.

**Symptom:** User opens "Add Signal", picks iRED, the dialog can't
select the bar mode and falls back to `static.table`. No panel renders.
End-to-end UI flow is broken for the 5 new TRs.

**Fix:** Add a new mode-kind constant + per-mode-kind dialog branch.
Look for how `static.table` is registered and clone the pattern for
each of the 4 new modes. Each new mode should preserve its mode-id
into the signal's `displayModeIds`.

Bonus: filter the dialog's descriptor list to also surface
`samplingStatus == Valid` static descriptors (currently the dialog
probably filters by `descriptor.temporal`, which would exclude the
new Valid-bypassed iRED/etc.).

### NOW-3: SequenceBarPanel multi-row overdraw on shared residue

**Where:** `h5-reader/src/app/SequenceBarPanel.cpp:131-158` — bar x
position is computed from `row.residue_index` only.

**Symptom:** Reorient's S²/τ_e/R1/R2/NOE panels have 3 rows per residue
(NH, CaHa, CO when all are present). Current code draws all three at
the same x; only the last drawn is visible. The bar chart is
silently wrong.

**Fix:** When `row.kind > 0`, offset the bar's x within the residue
slot by a fraction of `slotWidth` per kind. E.g. for 3 kinds:
NH at slot center, CaHa at center+0.33×slot, CO at center+0.66×slot
(or center ± offset based on kind).

```cpp
// Inside the bar-paint loop:
const double slotWidth = plot.width() / std::max(1.0, xSpan + 1.0);
const double barWidth = std::max(2.0, slotWidth * 0.75 / 3.0);  // narrower
const double kindOffset = (row.kind == 0)
    ? 0.0
    : (static_cast<double>(row.kind) - 2.0) * barWidth * 1.1;
const double xCenter = plot.left()
    + (row.residue_index - residueMin_) / std::max(1.0, xSpan) * plot.width()
    + kindOffset;
```

Same fix needed for `rowAtX()` (the click hit-test) so clicking
selects the right kind.

### NOW-4: Owned panels emit for all active signals regardless of active tab

**Where:** `h5-reader/src/app/DashboardDisplayController.cpp` —
`rebuild()` iterates `activeModel_->activeSignals()` without filtering
through `panelModel_`. The temporal-strip path uses
`seriesIsVisibleInActivePanel()` to filter; the owned-panel path
doesn't.

**Symptom:** Multiple dashboard tabs share the same owned-panel set;
adding iRED to tab A also shows it on tab B.

**Fix:** In the rebuild loop's panel-mode branch, check
`panelModel_->isSignalReferencedInActivePanel(signal.id)` (or whatever
the existing check is named — see how `seriesIsVisibleInActivePanel`
works for strip series) and skip if not.

Also: build a `DashboardDisplayRef` for the static panel modes the
same way strip modes do, so `removeDisplayRefsForSignal` cleanup
also drops the owned panel reference.

## LATER — defer to a new plan

**Next session pick-up:** these 10 items are organised into a
4-phase implementation plan at
`h5-reader/notes/PLAN_LATER_ITEMS_2026-05-29.md`.

These don't block use; they're polish/extension items.

### LATER-1: SceneRevealOverlay BondVector lookup helper

`h5-reader/src/app/SceneRevealOverlay.cpp` currently has duplicated
logic in `atomsForBinding` and `reveal()` for walking the iRED and
Reorient identity tables. Extract a free function
`lookupBondVector(h5, residue, kind) → optional<{tail, head}>` and
have both call sites consume it. Saves ~20 lines + makes Phase 2 (if
a third BondVector TR ever lands) trivial.

### LATER-2: HighFive narrowing/exception robustness

`h5-reader/src/io/QtTrajectoryH5.cpp` reader functions catch nothing.
If H5 has unexpected dtype (int64 where int32 expected, etc.),
HighFive throws and the whole trajectory load aborts. Wrap each
`Read*` reader body in try/catch, log the exception with the
specific group_path, leave the buffer null. The existing `WarnGroupAbsent`/
`WarnShapeMismatch` pattern is the destination.

Per Codex: int32 identity columns (`residue_index`, `n_atom`, etc.)
narrowly match the producer's dtype but if a future producer ever
emits int64, current reads silently narrow. Add a dtype assertion
or read into `int64_t` and cast down with range-check.

### LATER-3: REST `anchorToJson` BondVector arm

`h5-reader/src/app/RestServer.cpp:83` — `anchorToJson()` has no
BondVectorAnchor case. `/dashboard/signals` returns `{}` for the
anchor field on any iRED/Reorient signal. Add a `bond_vector` JSON
kind with `residue` and `kind` int fields.

### LATER-4: Mat3 tensor-glyph rendering

`h5-reader/src/app/SceneRevealOverlay.cpp` currently highlights bond
vector endpoints. Reorient's `bond_orientation_tensor` (per-vector
Mat3) is loaded into the buffer but not yet rendered as an ellipsoid
glyph attached to the bond. Add a `vtkPolyData` ellipsoid actor
parameterized by the tensor's eigendecomposition, lit and translucent
so it overlays the bond line.

### LATER-5: FixedFreqBlock J(ω) panel

Reorient's `spectral_density_j` is (V, 5) — J at 5 KTB Larmor
frequencies. Currently no dedicated panel; descriptor uses
`static.bar.sequence` placeholder. Build a small line-plot panel that
shows J(ω) at the 5 frequencies as a (log-x, log-y) scatter+line
plot, one per selected bond vector. Reuse PowerSpectrumPanel
mechanics with 5 discrete x positions instead of a dense F grid.

### LATER-6: Chi[0..3] dihedral autocorrelation

Phase F landed phi/psi only. The producer emits `chi_acf` (R, 4, L)
and `chi_corr_time_ps` (R, 4) for the 4 chi torsions per residue.
Add 8 more descriptors (4 chi corr_time scalars + 4 chi ACF curves)
following the existing phi/psi pattern. Or composite them under one
descriptor with PerClassBlock × 4 channels.

### LATER-7: Unify the three-bucket panel render order

`h5-reader/src/app/StripStackWidget.cpp` paints in three sequential
blocks: tracks_, then spectrumTracks_, then ownedPanels_. Phase B's
shortcut. For heterogeneous interleaving (a researcher comparing a
SequenceBarPanel against a TemporalStripPanel for the same residue),
the order should be controller-determined. Either:
- Unify all three into one `panels_` vector with the order the
  controller specifies, OR
- Add an explicit `panelOrderHint` per panel that the widget sorts on.

### LATER-8: Drop `mutable` on `ownedPanels_`

`h5-reader/src/app/DashboardDisplayController.h:183` — `mutable
std::vector<std::unique_ptr<AbstractStripPanel>> ownedPanels_`. The
only mutating user is `takeOwnedPanels()` which is non-const. The
`mutable` is unused. Drop it.

### LATER-9: SequenceBarPanel multi-channel overlay

Currently each panel renders one scalar across residues. Multi-overlay
(e.g. S² and τ_e on the same axis with twin-y) would let the user
compare relaxation observables in one chart. Defer until users ask.

### LATER-10: Status text for static-only signals

`DashboardDisplayController::updateStatusText()` (line ~1380) reports
"No active strip display modes" when only owned-panel signals are
active. Track a separate `activeOwnedPanelCount_` and report
combined.

## After-cleanup baseline check

Before declaring cleanup done, this must still pass:

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-debug
cmake --build . -j
ctest --output-on-failure
# 3/3 must pass
```

And the manual sanity (after-NOW-1+2 fixed): launch h5reader,
open Add Signal dialog, pick `h5:ired_s2`, confirm `static.bar.sequence`
appears in the mode picker, add the signal, observe the SequenceBarPanel
render across residues, then step the playhead and confirm the panel
STAYS rendered (the drain-on-tick bug is fixed).

## File inventory (untracked, ready to commit after cleanup)

New files (12):
- `h5-reader/src/app/AbstractStripPanel.{h,cpp}` — Phase B chassis
- `h5-reader/src/app/SequenceBarPanel.{h,cpp}` — Phase C
- `h5-reader/src/app/PowerSpectrumPanel.{h,cpp}` — Phase D
- `h5-reader/src/app/LagDecayPanel.{h,cpp}` — Phase D
- `h5-reader/src/app/ChordCouplingPanel.{h,cpp}` — Phase G
- `h5-reader/src/model/QtBondVectorBuffers.h` — Phase C/E
- `h5-reader/src/model/QtPerAtomChannelBuffers.h` — Phase D/G
- `h5-reader/notes/SCOPE_NEW_TRS_2026-05-29.md` — input scope doc
- `h5-reader/notes/CLEANUP_AFTER_PHASE_AH_2026-05-29.md` — this file

Modified files (15):
- `h5-reader/CMakeLists.txt`
- `h5-reader/src/app/DashboardDisplayController.{h,cpp}`
- `h5-reader/src/app/DashboardStripDock.cpp`
- `h5-reader/src/app/SceneRevealOverlay.cpp`
- `h5-reader/src/app/StripStackWidget.{h,cpp}`
- `h5-reader/src/io/QtTrajectoryH5.{h,cpp}`
- `h5-reader/src/model/DashboardSignal.{h,cpp}`
- `h5-reader/src/model/QtPerResidueBuffers.h`
- `h5-reader/src/model/SignalDictionary.h`
- `h5-reader/src/model/TrajectorySignalCatalog.cpp`
- `h5-reader/tests/dashboard_model_tests.cpp`
