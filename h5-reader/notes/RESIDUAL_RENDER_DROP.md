# Residual atom-render drop — what we know, what we've tried, what's left

**TL;DR.** One intermittent end-of-trajectory frame where molecule spheres
vanish while overlays (ring polygons, ribbon, streamlines) remain drawn.
The headline VTK bounds-cache bug is fixed (feedback_vtk_bounds_cache).
One sighting since the fix. Multiple plausible remaining causes, ordered
below by likelihood. Observability first — don't paper over the probe
added in step 7+ without understanding why it was needed.

**Guiding principle (from this session).** *Race conditions always
point to something deeper that needs love and clarity.* If a `Modified()`
added as a diagnostic probe makes the symptom disappear, the symptom is
not fixed — the true cause (broken Modified() chain, GPU upload ordering,
etc.) is still there waiting to surface under a different workload.

## What's already in place

- Per-frame bounds computed in one pass in `MoleculeScene::setFrame`
  and passed to `renderer_->ResetCameraClippingRange(padded)` with 5 Å
  padding. **Owns** the clipping-plane problem.
- Duplicate zero-arg `ResetCameraClippingRange()` removed from the
  camera-follow block (`MoleculeScene.cpp:~282`, was line 289).
- `mapper_->Modified()` added after `molecule_->Modified()` as a
  **diagnostic probe** (not a settled fix). Documented inline so the
  next reader knows it may be load-bearing or may be noise.
- Per-50-frame snapshot logging: RSS, actor count, actor bounds,
  actor visibility. Catches a repeat of the bounds-cache pattern or a
  new progressive-render degradation immediately.

## Four hypotheses, ordered by likelihood

### (A) Duplicate `ResetCameraClippingRange()` — DONE

The old zero-arg call at camera-follow time used `vtkActor::GetBounds`
whose cache does not invalidate on `SetAtomPosition + Modified`. On a
trajectory loop-wrap the delta is large; near/far lagged the true
bounds and spheres fell behind near while overlays with different
depth treatment stayed visible.

**Status.** Removed. The explicit-bounds call after overlay updates is
the single owner of clipping-range sync.

### (B) `Modified()` propagation through the composite mapper

`vtkOpenGLMoleculeMapper` builds internal sphere- and cylinder-imposter
mappers from the vtkMolecule. In principle `molecule_->Modified()`
propagates through the input connection and triggers re-upload on both
inner mappers. In practice there is at least one observed drop that
looks like a missed re-upload.

**Probe in place.** `mapper_->Modified()` added at `MoleculeScene.cpp`
after `molecule_->Modified()`. See inline comment.

**Experiment.**

1. Run a full-loop playback three times. If the drop never repeats,
   the composite chain was the cause and the probe is load-bearing.
2. If the drop does recur, the probe is noise — remove it and pursue
   (C) and (D).

**If probe is load-bearing — the real fix.**

Options, from least-to-most invasive:

- Trust the probe and document it as a known VTK quirk. Cheapest,
  least satisfying — "we add a Modified() because VTK doesn't always
  propagate one" is not a root-cause story.
- Replace `vtkOpenGLMoleculeMapper` with our own typed sphere and
  cylinder mappers driven off `vtkPoints` from QtFrame. Most invasive,
  but eliminates the composite-chain category of bugs entirely and
  lets us attach per-atom glyphs / colour bubbles for the time-series
  illustrator in one pipeline instead of layered overlays. See
  `TIME_SERIES_EXPANSION.md` for the overlap — this may be the right
  call when the next session starts.
- Call `molecule_->GetPoints()->Modified()` in addition to
  `molecule_->Modified()`. Middle ground: forces the specific upstream
  object the mapper's compute chain actually watches. Cheap experiment.

### (C) FXAA + translucency pass-ordering at end-of-trajectory

We deliberately use FXAA (cheap, AMD-friendly) instead of depth peeling
(qt6-cpp skill's default). FXAA runs after opaque + translucent passes;
if a translucent overlay (ring polygon edges, butterfly isosurfaces
with opacity < 1) fires a pass reorder at the exact frame the molecule
spheres draw, the spheres could be clipped by the stencil mask carried
into the post-pass.

**Experiment.**
1. Run the playback with `fieldGrid_->setVisible(false)` and
   `bfieldStream_->setVisible(false)` but rings + ribbon on. If the
   drop vanishes, a translucent overlay is implicated.
2. If isolated, disable FXAA (`renderer_->SetUseFXAA(false)` in
   `MoleculeScene::Build`) and retry. Comparison pass → order vs FXAA
   itself.

### (D) Picker / playback race on the last rendered frame

The user observed the last drop while clicking to stop playback near
end-of-trajectory. `QtPlaybackController::stop()` kills the timer;
`QtAtomPicker::handleMouseClick` runs a ray cast on the current frame.
Both touch VTK state on the GUI thread, but the click event can interrupt
a `setFrame` half-way through overlay fan-out. We issue exactly one
`renderWindow_->Render()` at the end of `setFrame`; a picker-triggered
`requestRender()` firing between molecule update and overlay update
would render spheres at frame t-1 while ribbon is at frame t.

**Experiment.**

1. Add a boolean `frameUpdateInProgress_` in `MoleculeScene`, asserted
   at `setFrame` start and cleared after Render. Have `requestRender()`
   short-circuit when asserted (no render, just log). If the drop
   disappears, we have the culprit and need to queue or drop the
   conflicting request.
2. Alternative diagnostic: log every entry into `setFrame` and
   `requestRender` with the current frame index and a monotonic call
   counter. The UDP trace should show the interleaving at the drop.

## What NOT to do

- Do not add a QTimer on frame rendering to "just rerender at the end
  to be safe." That papers over (D) and violates the Qt-discipline
  memory (feedback_qt_discipline).
- Do not `renderWindow_->Render()` inside overlay `setFrame`
  implementations. The single-render contract documented in
  `MoleculeScene.h` is load-bearing for (C).
- Do not remove the per-50-frame snapshot logging after viva. It cost
  nothing and it caught the bounds-cache bug the first time.

## How to instrument further if a drop recurs

In `MoleculeScene::setFrame`, immediately before `renderWindow_->Render()`:

```cpp
// Drop-forensics: one-line per frame with the state that matters.
qCInfo(cScene).noquote().nospace()
    << "frame=" << t
    << " mol.vis=" << (actor_->GetVisibility() ? 1 : 0)
    << " mol.bounds=[" << bounds[0] << "," << bounds[1] << ","
                        << bounds[2] << "," << bounds[3] << ","
                        << bounds[4] << "," << bounds[5] << "]"
    << " overlays.count=" << renderer_->GetActors()->GetNumberOfItems();
```

Then `tail -f` the UDP log during a playback and grep for the frame at
which the drop was observed. Compare bounds, visibility, and actor
count against known-good frames.

**The signal, not the noise.** Observability discipline: log enough
that the next incident shows the cause in one line of UDP. Speculating
beats guessing; data beats speculating.
