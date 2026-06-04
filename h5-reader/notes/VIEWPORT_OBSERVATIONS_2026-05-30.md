# Viewport observations — 2026-05-30

> **Status note (post-refactor):** the design landed at
> `spec/viewport_pipeline_2026-05-30.md` and the implementation landed
> in the same commit as this update. The post-refactor numbers live in
> `tests/scripts/HARNESS_BASELINE_PIPELINE_2026-05-30.md`. The sections
> below are preserved as historical record of the pre-refactor state —
> the cited file:line locations refer to the old `MoleculeScene.cpp`
> structure (centroid-delta block, internal `CameraPlaneLock` struct,
> single-layer renderer) that the refactor removed. Use the spec doc
> as the current contract; this doc as the why-we-did-it record.
>
> **Follow-up landed:** the `writeSubset` rotation half + the
> `POST /camera/focus_atom` convenience endpoint + small fixes from
> the Codex + Opus evaluative passes (camera-mode round-trip,
> QuietTrackballStyle::OnChar no-op, TransformedConformation CENSUS,
> DownAxis policy viewUp recompute, Picker render-source wiring,
> ReaderMainWindow.h docstring, pan scaling, writePlane policy
> assert, modeChanged toolbar sync). Two canonical recipes for the
> "focus atom + surroundings legible" goal are documented in
> `tests/scripts/HARNESS_BASELINE_PIPELINE_2026-05-30.md`.
>
> **June 4 truth note:** camera-mode sections are partly superseded by
> built atom/bond/dihedral/subset composition and plane lock. The local
> display-isolation observation remains current: no pick atom + radius +
> hide-nonlocal-atoms/residues UI exists yet.

Senior-engineer reading notes after a day with `h5-reader/src/app/`.
Output of the first pass on the viewport refactor scope named in
`QT_PRIMITIVES_ALIGNMENT_2026-05-30.md` (Group B + viewer animation
refactor) and `ROBUSTNESS_BACKLOG_2026-05-30.md` item 6 (rendering
probe). Engagement with the material, not a plan — that's the next
conversation.

The user's correction during this session, which anchors everything
below: *the current code doesn't really work, it all bounces around,
and both Claude and Codex have said "you can't have that" — but it's
VTK and of course you can; the current framework is the limiter, not
the toolkit.* Saved as `feedback_vtk_capable_dont_underclaim`. The
rest of this doc is what "the framework is the limiter" means
concretely.

## 1. Why it bounces (specific causes, by file:line)

The current frame-advance hot path is `MoleculeScene::setFrame`
(`MoleculeScene.cpp:418-592`). Three classes of camera state coexist
inside it; they interact badly. Each by itself is fine — the failure
mode is the composition.

### 1.1 Centroid-delta translation accumulates noise and fights everything else

`MoleculeScene.cpp:446-461` computes the per-frame centroid as the
mean of all atom positions and stores it in `lastCentroid_`.
`MoleculeScene.cpp:476-497` translates the camera by `(new_centroid
- last_centroid)` every frame, except when the explicit plane lock
owns the camera. This is the bouncing.

Three things compound:

- **Vibrational noise** is on the order of 0.1-0.3 Å rms per atom at
  300 K. For a 1500-atom protein the centroid is steadier, but its
  per-frame jitter is still tens to hundreds of milliangstroms — and
  it's correlated with the user's eye-direction expectation because
  the camera shifts in the SAME plane as the molecule's perceived
  motion. Vibrations the eye would smooth out become visible camera
  micro-jitter.
- **PBC wraparound**. When the protein crosses the simulation box
  boundary in the GROMACS run, the H5 positions wrap. A wrapped
  atom contributes a teleport-sized delta to the centroid for one
  frame. The camera follows it. The next frame the centroid is back
  but the camera is several Å away from where it was. Repeat.
- **Floating-point drift over hundreds of frames**. Even without PBC,
  delta accumulation is the wrong primitive. After 750 frames of
  1P9J the camera has been `+=`d 749 times. The composed transform
  is no longer what an absolute calculation would have produced.

This is the centroid-follow that was *useful* for the 1B1V_4292 diffusion
case the comment at `MoleculeScene.cpp:434-445` documents — when the
protein drifts across the simulation box over 25 ns, you need the
camera to stay with it or the molecule renders off-screen. But the
*implementation* — delta translation off the previous camera state
— is fragile. The right primitive is absolute: `SetFocalPoint(F)`
+ `SetPosition(F - V*D)` where F, V, D are computed each frame from
a target spec (atom, centroid, bond midpoint, dihedral central bond
midpoint), not derived by delta from the previous camera.

### 1.2 Camera input via Qt eventFilter, not VTK interactor

The h5-reader's Qt+VTK precedent for mouse handling is the picker:
`QtAtomPicker.cpp:48` installs `eventFilter` on the
`QVTKOpenGLNativeWidget` and intercepts `MouseButtonDblClick` before
VTK sees it. The same pattern is the primary path for camera
drag/wheel.

Today `MoleculeScene.cpp:88-91` installs the stock
`vtkInteractorStyleTrackballCamera`. Each manipulation method ends
in `rwi->Render()` — five sites
(`vtkInteractorStyleTrackballCamera.cxx:272, 297, 348, 393, 451`)
plus the wheel-dolly path at `:201-232`. These are independent
render triggers from the playback timer; ordering between them on
the same Qt tick is unspecified. That's the user-fights-camera
shape.

Qt-aware primary path (extends `QtAtomPicker` pattern):

1. Install an eventFilter on `vtkWidget_` that catches
   `MouseButtonPress / MouseMove / MouseButtonRelease / Wheel` for
   camera-drag intent. Camera math (`Azimuth/Elevation/Pan/Dolly`)
   moves to free helper functions on a `CameraGesture` value type
   — pulled out of the stock trackball, not invented. The filter
   schedules a render through the single coalesce path (section 5b)
   and returns `true` so VTK never sees the event.
2. Leave the stock trackball installed for events the eventFilter
   doesn't intercept (3D mice, hover, future touch). Call
   `interactor_->EnableRenderOff()` — the trackball's `rwi->Render()`
   hits the gate at `vtkRenderWindowInteractor.cxx:168-177` and
   no-ops. Camera mutation still works for the non-intercepted
   path; render scheduling stays ours.
3. Call `interactorStyle_->AutoAdjustCameraClippingRangeOff()`. The
   stock trackball's no-arg `ResetCameraClippingRange()`
   (`vtkInteractorStyleTrackballCamera.cxx:264`) hits the
   stale-bounds path documented at `feedback_vtk_bounds_cache`.
   Disable; rely on the per-frame explicit-bounds call at
   `MoleculeScene.cpp:547-553`.

Second-order benefit: gesture-aware locks. With camera input
flowing through Qt-side handlers, a 1-atom lock can intercept
rotation gestures and re-interpret them as "rotate around the
locked atom" — recompute the rotation center per frame, then call
the helper. The trackball today spins around a fixed focal-point
and a lock has to fight it; the eventFilter path doesn't.

### 1.3 The only mode that doesn't bounce is the 3-atom plane lock

`MoleculeScene.cpp:282-415` is the working precedent. The shape:

- `lockCameraToSelectionPlane(atoms)` captures the user's current
  camera state (view-up direction, distance, normal sign convention)
  into a `CameraPlaneLock` struct, expressed in the LOCAL frame of
  the 3 atoms (`localViewUp` is stored as (x, y, z) components in
  the plane's basis, not in world space).
- `applyCameraPlaneLock(frame)` recomputes the world-space camera
  from scratch each frame: builds the plane's orthonormal basis from
  the current frame's atom positions, transforms `localViewUp` back
  to world via the basis, computes the view direction with a sign-
  flip continuity guard (lines 390-395, `lastDirection.dot < 0` →
  flip `normalSign`), and writes `SetFocalPoint` + `SetPosition` +
  `SetViewUp` + `OrthogonalizeViewUp`.
- `setFrame` calls `applyCameraPlaneLock` BEFORE the centroid-delta
  branch and short-circuits the latter when the lock is active
  (`MoleculeScene.cpp:463-471`).

This is the sophisticated VTK pattern at the cost of ~130 lines. The
camera is absolute-written from a typed lock state, the normal-sign
continuity handles ring-flip-style degeneracies, and the basis is
recomputed per frame so user-driven view-up changes survive across
playback. It just doesn't bounce.

The lesson: **everything else should follow the same shape**, not
the centroid-delta shape.

## 2. What's missing (mapped to the user's stated wants)

The user named four capabilities in the kickoff prompt and the
follow-up: atom-locked centering, atom+bond centering, dihedral
sight-down, distance-filtered display. AtomSelection caps at 4 atoms
(`AtomSelection.h:66`, `kMaxAtoms = 4`) and already classifies the
selection as `Distance` (2) / `Angle` (3) / `Dihedral` (4) via
`GeometryKind` (`AtomSelection.h:56`). So *the data side* already
knows about all four cardinalities and the corresponding geometric
meanings. What's missing is one CameraLock-per-cardinality and the
generalised plane-lock pattern applied to each.

### 2.1 Atom-only lock (1 atom, missing)

Should be the simplest mode: pick an atom, F = atom position each
frame, V and D inherited from current camera at lock acquisition (so
the view doesn't teleport when the lock turns on). Re-applied per
frame from the absolute atom position, not delta-translated. No
basis, no sign-flip — just `SetFocalPoint(atom_pos_t)` +
`SetPosition(atom_pos_t - V*D)`.

This is what the user means when they ask for "the atom stays
centered through playback." It does NOT exist today; the closest is
`QtSelectionOverlay` (`QtSelectionOverlay.h:1-67`, marked dormant)
which used to move a sphere to the picked atom but never affected
the camera.

### 2.2 Atom+bond lock (2 atoms, missing)

Two atoms define a bond axis B = (atom2 - atom1).normalized(). F =
midpoint(atom1, atom2) (or atom1 if the user wants atom1 specifically
centered — design choice). V = inherited from current camera at lock
acquisition, projected to be perpendicular to B if the user wants
the bond visible (otherwise V = B for "look down the bond" which is
the bond-as-sight-axis variant). ViewUp = B - (B·V)V normalised, so
the bond axis is horizontal in screen space.

This is the camera-stays-aligned-with-the-bond mode. None of this
exists today.

### 2.3 Dihedral sight-down lock (4 atoms, exists ONE-SHOT, missing AS LOCK)

`MoleculeScene.cpp:624-694` has the right math in `focusCameraOnReveal`:
when the binding's anchor is an `AtomTupleAnchor` with ≥4 atoms, the
function sights down the axis (atom3 - atom2) (the central bond of
the 1-2-3-4 torsion) and sets the camera position behind atom2 along
that axis. ViewUp is preserved-then-orthogonalised against the new
view direction. This is the Newman-projection view — standard
chemistry.

But it's a ONE-SHOT call from `revealBinding` (`MoleculeScene.cpp:594-613`).
The camera is set once when reveal fires, then the next setFrame
delta-translates it. No per-frame re-application. As soon as the
user advances the frame, the bond rotates relative to the camera
and the Newman projection is lost.

The fix is to lift this math into a `CameraDihedralLock` with the
same per-frame re-apply discipline as the plane lock. The basis here
is straightforward: V = (atom3 - atom2).normalized() (sign chosen to
respect the user's current view direction at lock acquisition,
continuity-guarded between frames the same way the plane lock does),
F = (atom2 + atom3) * 0.5, ViewUp = (atom1 - atom2) projected
perpendicular to V (or atom4 - atom3, design choice — they should
agree at zero dihedral and differ proportional to the torsion at
nonzero, which is exactly what you want to see).

### 2.4 Distance-filtered atom display (missing entirely)

`vtkOpenGLMoleculeMapper` has no per-atom visibility API. Confirmed
against the VTK source: ghost arrays (`AllocateAtomGhostArray`) are
honoured by atom rendering but NOT by bond rendering, which would
leave bonds dangling. The clean approach is a separate "display
molecule" — a subset of the source `vtkMolecule` rebuilt on filter
change.

Sketch (NOT spec yet — for the design conversation):

- Keep `molecule_` as authoritative (all atoms, all bonds).
- Add `displayMolecule_` plus a forward index map `displayIdx[source]
  → display_idx_or_NPOS`. On filter change (new cutoff, new focus
  atom), walk source atoms, decide which ones survive, build the
  display molecule with AppendAtom/AppendBond preserving identity-
  to-source via a per-display-atom scalar array. Re-point the mapper
  with `SetInputData(displayMolecule_)`.
- On per-frame setFrame, only `SetAtomPosition` for the kept display
  atoms (cheap, ~5000 calls max).
- AppendAtom + AppendBond rebuild costs are O(N+B) but with constant
  factors well under 1 ms for 5000 atoms. The expensive bit is the
  GPU re-upload on the first render after rebuild — unavoidable.
- For the distance-from-moving-atom case (filter centred on the
  currently-locked atom), refresh the cutoff per frame too. At
  1500 atoms × O(N) pairwise distance, the per-frame filter is well
  under a millisecond; the GPU re-upload on subset change is the
  cost. Cache the previous frame's subset; only rebuild + re-upload
  when the subset membership actually changes.

The `Filtered` (or `Lensed`?) view layer is conceptually parallel
to the `CameraLock` layer — both are policies on the same underlying
data, both refresh per frame, both should be first-class objects in
the scene rather than functions buried in setFrame.

## 3. The architectural shape this points at

Not a plan, but the conversation seed:

Today: `setFrame` updates atom positions, runs centroid-delta
camera translation, conditionally hands the camera to the
plane-lock, fans `setFrame` to overlays, computes bounds, calls
`Modified()` twice, `ResetCameraClippingRange`, `Render`, logs.

Wanted: `setFrame` still drives everything from the one timer's
`frameChanged`, but the camera comes from a typed lock object
that writes the camera absolute each frame. Centroid-delta
accumulation goes away. The 3-atom plane lock pattern generalises:

- **`CameraLock`** — typed object, variants Centroid / Atom / Bond /
  Dihedral / Plane / Free. Each holds lock-local state (the
  equivalent of today's `CameraPlaneLock` struct at
  `MoleculeScene.cpp:148-154`). Each frame: compute focal point F,
  view direction V, view-up U, distance D from the conformation at
  the current frame; call `SetFocalPoint(F)` + `SetPosition(F -
  V*D)` + `SetViewUp(U)` + `OrthogonalizeViewUp()`. No delta
  translation. The plane-lock's per-frame `applyCameraPlaneLock`
  (`MoleculeScene.cpp:362-416`) is the template for each variant.
- **`AtomFilter`** — typed object, variants All / WithinRadiusOf(atom,
  R) / ResidueSubset. Maintains a display `vtkMolecule` that is a
  subset of the source. On filter change: rebuild via AppendAtom /
  AppendBond, `mapper->SetInputData(displayMolecule_)`. On frame
  change: `SetAtomPosition` on the kept atoms only. The 5000-atom
  rebuild cost is well under a millisecond; the expensive part is
  the first-render GPU re-upload after subset change.
- **Input**: Qt eventFilter on `vtkWidget_` is the primary camera
  drag/wheel path (section 1.2; matches `QtAtomPicker` precedent).
  Stock trackball stays installed for events the filter doesn't
  intercept, gated via `EnableRenderOff()` +
  `AutoAdjustCameraClippingRangeOff()`. Camera math
  (Azimuth/Elevation/Pan/Dolly) lives in free helper functions, not
  inside an interactor subclass.

Smoothing (`vtkCameraInterpolator` + `vtkQuaternionInterpolator`)
is future work, not on the path. The current bouncing is
centroid-delta accumulation + clipping-range chatter, not
interpolator drift — confirm with the harness (section 5b) before
adding smoothing infrastructure. Worth noting: `vtkTransform::SetupCamera`
internally orthogonalises the modelview matrix via cross products
regardless of stored ViewUp orthogonality, so the visible projection
is fine even with non-orthogonal stored ViewUp. The ViewUp-drift
failure mode is the cross-product gradient flipping sign when ViewUp
approaches parallelism with the view direction — a sudden ~180° spin
of the side-axis, not a continuous wobble. Relevant only if/when we
add interpolation between sparse keyframes.

`PlaneFrameMath.h` (1-79) is the test model — pure-function math
headers extracted from the scene so the camera geometry is
unit-testable without a VTK renderer or a loaded conformation. Every
new lock variant gets the same treatment.

## 4. The paint-some-frames bug — diagnostic angle pivots

`MoleculeScene.cpp:501-528` is the probe (`mapper_->Modified()` after
the atom-position loop). Source comments call it a probe, not a fix.

Tracing the VTK pipeline against `/home/jessica/builds/VTK-9.5.2/`:
the rebuild gate at `vtkMoleculeMapper.cxx:320-330` is sound.
`vtkMolecule::SetAtomPosition` bumps MTime per call
(`vtkMolecule.cxx:184-197`); molecule MTime exceeds AtomGlyphPolyData
MTime by N (atom count) at the end of the position loop;
`UpdateAtomGlyphPolyData` rebuilds; `Initialize()` bumps the polydata
MTime (`vtkDataObject.cxx:134-151`); the OpenGL mapper's VBO gate at
`vtkOpenGLPolyDataMapper.cxx:3742-3767` sees the bump via
`CurrentInput->GetMTime()` and re-uploads. The probe `mapper_->Modified()`
is redundant via the third clause `mapper->GetMTime() > polydata MTime`.

More plausible cause: clipping-range chatter. The camera position
drifts via the centroid-delta path at `MoleculeScene.cpp:476-497`.
When the protein diffuses or the user dollies in, the camera can
land inside (or very near) the protein bounds. `ResetCameraClippingRange(padded)`
at `MoleculeScene.cpp:547-553` then computes near/far from the
per-frame bounds + the camera position; the near plane can end up
very close to imposter quads at the front of the molecule.
`vtkOpenGLSphereMapper.cxx:102` has `if (d < 0.0) discard` in the
sphere intersection fragment shader — fragments near the near plane
fall through this discard. Atoms thin out at the front;
end-of-trajectory frames look like atoms dropped.

Diagnostic hook for the "did the rebuild actually fire" question:
`vtkOpenGLMoleculeMapper::GetFastAtomMapper()`
(`vtkOpenGLMoleculeMapper.h:40`) is public. The FastAtomMapper's
input is the AtomGlyphPolyData via the trivial producer
(`vtkOpenGLMoleculeMapper.cxx:53`). Attach a `vtkCommand::ModifiedEvent`
observer to `fastMapper->GetInputDataObject(0, 0)`. UDP-log per
fire. If rebuilds fire on every frame, the gate works — focus on
the clipping range / fragment discard. If a frame is missing a
rebuild, the gate has a bug we haven't found.

Experiments via the harness (section 5b): mark a front-facing atom
magenta; snapshot per frame across the end-of-trajectory range;
locate the magenta blob's pixel count per frame. Drop in pixel
count = fragments discarded. Correlate with logged near-plane
values from the per-frame bounds calc.

Two cheap probes worth trying before deep instrumentation:

- **`molecule_->GetPoints()->Modified()`** after the SetAtomPosition
  loop in `MoleculeScene::setFrame`. `vtkMolecule::SetAtomPosition`
  bumps the molecule's MTime but not the underlying `vtkPoints`
  MTime explicitly. The `vtkOpenGLPolyDataMapper` VBO gate checks
  several MTimes; if the points' MTime is the one a particular
  code path consults, the gate misses. If end-of-trajectory drop
  vanishes, this is the bug. One line.

- **`SetForceTranslucent(true)` on the BS/HM isosurface actor**.
  `vtkOpenGLRenderer.cxx:520-690` polls each actor's
  `HasTranslucentPolygonalGeometry()` per-frame to dispatch
  opaque-vs-translucent passes. For the isosurfaces, per-cell
  opacity computed near boundary values (0.0/1.0) can return
  different classifications frame-to-frame. Pinning translucent
  stabilises the classification and may resolve the BS/HM stacking
  bug (`POLISH_BACKLOG.md` item 0) as a side effect.

Per `feedback_race_condition_deeper_issue`: the probe stays
documented as an open question rather than declared a fix; this
doc IS that documentation.

## 5. Coupling with the in-flight Group B refactor

The Group B items in `QT_PRIMITIVES_ALIGNMENT_2026-05-30.md` (5 + 6)
migrate `Conformation::requestSnapshot` / `snapshotReady` and
`DftShieldingStore::requestFrame` / `frameReady` to
`QFuture<...>` + `QFutureWatcher`. The viewer animation refactor is
named as the natural co-landing because both touch the same code
paths.

Concretely: this viewport refactor will repeatedly visit
`MoleculeScene::setFrame` and the per-frame fan-out from
`QtPlaybackController::frameChanged`. The current snapshot-subscriber
chain (which lives in `DashboardDisplayController`, the inspector
dock, etc., per `ReaderMainWindow.cpp:245-247`) does NOT touch
`MoleculeScene` directly — scene gets positions through
`Conformation::atomPosition` synchronously, not through snapshots.
So the QFuture migration on snapshots and the camera/scene refactor
are *parallel work in adjacent files*, not one work that must
happen together. Item 5 + 6 can land separately. They're paired in
the QT_PRIMITIVES doc because the SAME engineer should know both
when modifying either, but they're not architecturally entangled.

The slider-drag pile-up fix already partially landed via
`DashboardDisplayController::scrubActive_` (set on
`QSlider::sliderPressed/Released`, wired at
`ReaderMainWindow.cpp:347-358`). Snapshot fetches now defer during
drag. This is fine; doesn't need to be revisited as part of the
viewport refactor.

## 5b. Instrumentation-first (added 2026-05-30 mid-session)

User correction during this session: *"I think we should add a REST
service to snapshot frames, and the ability to also specify a
selection colour of choice and to select multiple atoms via REST. I
suspect this will require trying and rejecting hypotheses and that
will be a very different matter if you have that plus proper UDP
logging set up."*

This reframes the next session's agenda. Before any of the lock
taxonomy / interactor subclass / smoothing work, build the
experimental harness. The deep VTK source trace produced a wall of
hypotheses (paint-cycle inversion via `QVTKRenderWindowAdapter.cxx:241-299`
where `renderWindow->Render()` only schedules the paint via
`ParentWidget->update()`; `vtkInteractorStyleTrackballCamera` issuing
its own `rwi->Render()` at five sites including line 272, fighting
our QTimer-driven setFrame; `vtkCameraInterpolator.cxx:435-470`
setting ViewUp without renormalising or orthogonalising;
`vtkPoints::SetPoint` not always bumping the points array's own MTime
so the mapper's MTime gate misses; `vtkRenderer::ResetCameraClippingRange`
no-arg form hitting the stale `vtkActor::GetBounds()` cache);
acting on those one-at-a-time without an experimental harness
produces a multi-session loop of "this should fix it" then "actually
that wasn't it." See [[feedback-instrument-before-refactor-qtvtk]].

**What needs to be wired** (REST endpoints, in priority order):

1. **`POST /screenshot {"target": "scene"}` — already exists** at
   `RestServer.cpp:472-491`; `captureScenePng` at `:126-148` uses
   `vtkWindowToImageFilter` + `vtkPNGWriter` to a
   `vtkUnsignedCharArray`. Two extensions needed for the harness:
   - Query parameter `?force_render=true|false`. The default for
     `vtkWindowToImageFilter` is `ShouldRerenderOn` (forces a fresh
     render before reading pixels). `force_render=false` calls
     `w2i->ShouldRerenderOff()` before `Update()` so the snapshot
     reads whatever is in the FBO — the right mode for the
     paint-cycle-inversion experiment.
   - `ASSERT_THREAD(this)` at the top of `captureScenePng`.
     `QHttpServer` lambdas run on the main thread by default, but
     this is the GL-touching path; verify explicitly.
2. **`POST /selection/atoms`** with body `{"atoms": [int, ...]}` —
   bulk-set the AtomSelection (today only `/selection/pick` exists,
   which simulates one atom + modifiers at a time; for scripted
   experiments we want one call that sets the whole 4-tuple). Must
   round-trip through AtomSelection's typed mutators so the
   GeometryKind classification + downstream signals fire identically
   to the GUI gesture path.
3. **`POST /selection/color`** with body `{"slot": int, "rgb":
   [r,g,b]}` — overrides the fixed Okabe-Ito palette in
   `AtomSelection::SlotColorRgb` for a given slot. Need to add a
   per-instance override map to AtomSelection (the static palette
   stays as the default). MeasurementOverlay reads the colour at
   sphere-build time; on override, the model emits a "colours
   changed" signal that triggers overlay re-colour. Use case: mark
   "the atom I expect to bounce" in bright red, snapshot per-frame,
   visually confirm whether it moves on-screen between frames where
   it shouldn't.
4. **UDP discipline audit + render-source tagging.** Per
   `feedback_udp_logging` + the chewer zero-write hot-path note,
   verify the per-frame setFrame log path is UDP-only. Stderr
   restricted to startup/shutdown/error.

   The render-source tag is the highest-value diagnostic add: the
   EndEvent observer at `MoleculeScene.cpp:101-111` already exists
   and logs render time. Set a thread-local `lastRenderSource_`
   (one of `timer/picker/rest/plane-lock/interactor/paint`) inside
   the new `requestRender()` coalesce path (item 5 below) before
   the queued invocation runs; the EndEvent observer reads it and
   UDP-logs `source=...` alongside the render time. Python harness
   then correlates `/frame/set` requests against actual render
   sources — if a render fires with `source=interactor` inside a
   `source=timer` interval, the trackball is still firing renders
   despite `EnableRenderOff()`.

5. **Render coalesce primitive.** Replace the four direct
   `renderWindow_->Render()` call sites
   (`MoleculeScene.cpp:237, 242, 337, 555`) with a coalesced
   `requestRender()`. `widget_->update()` is the Qt-aware render
   verb (schedules a paint event which goes through the adapter's
   `paint()` → `iren->Render()` chain at
   `QVTKRenderWindowAdapter.cxx:241-266`); coalescing via
   `QMetaObject::invokeMethod` with `Qt::QueuedConnection` drains
   multiple per-tick callers to one render:

   ```cpp
   void MoleculeScene::requestRender(RenderSource src) {
       ASSERT_THREAD(this);
       lastRenderSource_ = src;   // read by EndEvent observer
       if (renderPending_) return;
       renderPending_ = true;
       QMetaObject::invokeMethod(this, [this] {
           renderPending_ = false;
           if (vtkWidget_) vtkWidget_->update();
       }, Qt::QueuedConnection);
   }
   ```

   No `QTimer::singleShot(0)` — queued invocation drains on the
   next event-loop tick without adding a timer. The one application
   timer stays the playback timer.

**What this enables — experimental design per hypothesis**:

- **Paint-cycle inversion**: `POST /frame/set` then immediately
  `GET /scene/snapshot` (snap-A), then `GET /scene/snapshot?wait_ms=100`
  (snap-B). If snap-A != snap-B, the synchronous Render didn't reach
  the user-visible framebuffer until the next paint event — confirms
  the QVTKRenderWindowAdapter inversion documented at
  `QVTKRenderWindowAdapter.cxx:241-299`. Fix candidate: drive renders
  via `widget->update()` not `renderWindow_->Render()`.
- **Trackball-vs-timer race**: while playback is running, drive
  synthetic mouse-move events via QTest or by exposing a
  `POST /scene/synthetic_mouse_move` endpoint. Snapshot per frame +
  count UDP "render triggered by" log lines per frame. > 1 render
  per frame confirms the trackball's `rwi->Render()` at
  `vtkInteractorStyleTrackballCamera.cxx:272` competing with our
  timer. Fix candidate: `interactor_->EnableRenderOff()` + funnel
  all renders through `widget->update()`.
- **AtomGlyphPolyData rebuild miss** (end-of-trajectory drop): add a
  ModifiedEvent observer on the molecule mapper's internal
  AtomGlyphPolyData (need a friend-access shim or a temporary
  subclass; `vtkMoleculeMapper` exposes `GetAtomGlyphPolyData()` per
  9.5 source). UDP-log every rebuild. Playback frame 0..749 and
  count rebuilds. If end-of-trajectory frames don't log a rebuild,
  the gate failed — confirms our hypothesis-1 from the trace agent.
  Fix candidate: also bump `molecule_->GetPoints()->Modified()`
  before `molecule_->Modified()`.
- **ViewUp drift in dihedral lock**: pin a 4-atom selection, enable
  dihedral lock (once it exists), use `/selection/color` to mark
  atom 1 bright red. Snapshot every frame across playback. In
  ImageMagick / PIL, locate the red blob's screen-space position
  per frame. If the position drifts predictably with frame index,
  the ViewUp interpolator isn't renormalising — confirms
  `vtkCameraInterpolator.cxx:435-470`. Fix candidate: skip the
  interpolator for sustained-orientation locks; compute view-up
  per frame from the live atom geometry + `OrthogonalizeViewUp`.
- **Centroid-delta accumulation**: enable centroid follow, mark
  atom 0 red, snapshot every 10 frames across 750 frames. Track
  the red blob's screen position. Drift over playback ≠ noise
  confirms accumulation. Fix candidate: replace delta translation
  with absolute `SetFocalPoint(centroid_t)`.

**The marker-colour design is load-bearing.** User correction
2026-05-30 mid-session: *"You will need a special colour to be able
to see the selection and agreement not to give up and say oh we
can't see and go back to 5-minute-a-round user sessions. The only
way near-pixel but not pixel agreement is going to work is if you
have a quickly discernible set of fixed points via colour of atom
selection for this work."*

Pixel-perfect snapshot comparison is doomed for our renderer. FXAA
blurs edges (`MoleculeScene.cpp:80-81`: FXAA on, depth peeling off);
`vtkOpenGLSphereMapper`'s imposter shader rasterises spheres with
sub-pixel coverage in its fragment shader; GPU drivers can produce
indeterministic edge pixels frame-to-frame on the SAME camera state.
Two snapshots of "the same frame" will diff in noise even when
nothing is wrong.

**The harness has to be centroid-of-blob analysis, not pixel diff.**
Mark target atoms with a HIGH-CONTRAST colour that does NOT appear
anywhere else in the rendered scene — not in any element's CPK
colour, not in the background, not in any ribbon / ring / overlay.
Then find the centroid of each marker blob per snapshot via
connected-component analysis (PIL / OpenCV cv2.connectedComponents
or numpy `label` from scipy.ndimage). The marker's screen-space
position is sub-pixel accurate even when the surrounding render is
noisy.

**Specific colour design** (CPK collisions to avoid: O=red, N=blue,
H=white, C=black/dark-grey, S=yellow, P=orange, Cl=green, F=cyan,
Br=brown, I=purple-ish, Fe=orange-brown). Safe choices:

| slot | hex       | rgb            | name              | rationale |
|------|-----------|----------------|-------------------|-----------|
| 0    | `#FF00FF` | (255, 0, 255)  | pure magenta      | not in CPK; max distance from all CPK in HSL; eye-catching |
| 1    | `#00FF7F` | (0, 255, 127)  | spring green      | distinct from chlorine's pure green (#00FF00); high V/S |
| 2    | `#FF1493` | (255, 20, 147) | deep pink         | distinct from magenta but in similar hue family for "this is a marker" cue |
| 3    | `#9D00FF` | (157, 0, 255)  | vivid violet      | distinct from purple iodine; high contrast against most backgrounds |

Markers render SOLID (opacity=1.0), NOT translucent. Today
`MeasurementOverlay` uses 0.5 opacity, 0.85Å radius
(`MeasurementOverlay.cpp:21-22` constants), and lives in the main
renderer — both too weak for blob detection AND occludable by
geometry. For instrument mode: opacity=1.0, ~1.5Å radius, **and
move the actors to a second renderer at `SetLayer(1)`**
(`renderWindow_->SetNumberOfLayers(2)`; overlay renderer shares
the active camera with the main renderer for screen-space
alignment). Layered renderers composite back-to-front per layer
with depth reset between, so markers paint on top regardless of
whether the underlying atom is buried in the protein interior. The
overlay renderer can also disable FXAA to keep marker edges crisp
for centroid analysis. This supersedes the depth-test-off-on-property
approach (cleaner composition with the BS/HM translucent pass).

`MeasurementOverlay` is the right home for the implementation —
inject the overlay renderer into its constructor instead of the
main renderer for instrument mode. Per-actor color, opacity, and
radius become setters on the overlay; today they are
anonymous-namespace constants.

**Instrument-mode preset** (suggested REST endpoint):
`POST /selection/instrument` with `{enabled: true}` toggles all slot
colours to the table above AND forces opacity=1.0 AND boots
visibility on the MeasurementOverlay. Reversible via `{enabled: false}`
(restores Okabe-Ito + default opacity). This is a single-button
"make markers find-able" mode. The per-slot `/selection/color`
override stays for finer control.

**Workflow commitment** (the user explicitly asked for this in
writing): the next session does NOT hit a snag with the renderer and
fall back to 5-minute-round user check-ins. The harness produces
deterministic results; hypotheses are accepted or rejected by the
data; rejected hypotheses lead to the next hypothesis without
asking permission. The user is consulted on **design choices**
(should the lock release on selection change? should dihedral
sight-down preserve user rotation?), not on **whether the previous
experiment worked**. The harness is the answer to "did it work."
This commitment is binding for the viewport-refactor work and
documented in `feedback_instrument_before_refactor_qtvtk` memory
entry.

**Python harness sketch** (for the design conversation):

```python
# scripts/viewport_probe.py — drives experiments against running reader
import requests, hashlib
from pathlib import Path

BASE = "http://127.0.0.1:9988"

def screenshot(force_render: bool):
    r = requests.post(f"{BASE}/screenshot",
                      json={"target": "scene"},
                      params={"force_render": "true" if force_render else "false"})
    return r.content

def probe_paint_inversion(frame):
    requests.post(f"{BASE}/frame/set", json={"frame": frame})
    # snap_a: read FBO as-is (post-frame-set, before any new render)
    snap_a = screenshot(force_render=False)
    # snap_b: force a fresh render then read
    snap_b = screenshot(force_render=True)
    return hashlib.md5(snap_a).hexdigest() != hashlib.md5(snap_b).hexdigest()

# Run across frames; high rate = paint events not draining before snap_a fires.
inversions = sum(probe_paint_inversion(f) for f in range(0, 750, 25))
print(f"paint-cycle inversion rate: {inversions/30:.0%}")
```

The harness lands in `h5-reader/tests/scripts/viewport_probe.py`
(new dir; companion to existing `tests/rest/` pytest suite). Lives
alongside the reader so the experiments are reproducible per build.

**Existing RestServer surface area** (`src/app/RestServer.cpp`,
494 lines) confirms the foundation:

- `/protein/atoms`, `/selection`, `/selection/pick`, `/selection/clear`
- `/frame/current`, `/frame/set` (auto-pauses playback — REST
  experiments are always paused, no timer race)
- `/plane-lock`, `/plane-lock/enable`, `/plane-lock/disable`
- `/positions`, `/dashboard/signals`, `/shutdown`
- **`/screenshot`** — `{target:"scene"}` or `{target:"window"}`,
  returns PNG bytes via `vtkWindowToImageFilter` /
  `QWidget::grab()`

Adding `/selection/atoms` + `/selection/color` + `?force_render`
query on `/screenshot` fits the existing route + `jsonResponse` +
`ASSERT_THREAD` pattern. No architectural deviation; the
colour-override hookup in AtomSelection is ~30 lines, the
two-renderer overlay-layer change in MeasurementOverlay is ~20
lines (constructor takes overlay renderer instead of main), the
render-coalesce primitive is ~20 lines in MoleculeScene + reorder
of currentFrame_ assignment (`MoleculeScene.cpp:557` moves before
the requestRender call so the EndEvent observer reads the
correct frame). Plus the shutdown reorder
(`vtkWidget_->setRenderWindow(nullptr)` before drop of
`renderWindow_` smart pointer, see Section 6). Half a day total.

## 6. Qt+VTK proper-use summary

From the Qt+VTK specialist pass 2026-05-30 (agent loaded with the
qt6-cpp skill, per `feedback_agents_qtvtk_load_skill`). Today vs
proper Qt+VTK for this setup:

| Topic | Today | Proper Qt+VTK |
|---|---|---|
| Camera input | Trackball renders, app re-renders, they fight | Qt `eventFilter` on widget (`QtAtomPicker` precedent); trackball `EnableRenderOff` + `AutoAdjustCameraClippingRangeOff`; camera math as free helpers |
| Render trigger | `renderWindow_->Render()` at 4 sites | `vtkWidget_->update()` via coalesced `requestRender()`; one drain via `QMetaObject::invokeMethod(... Qt::QueuedConnection)` |
| Camera per frame | Centroid-delta `+=` on focal-point (accumulates) | Absolute `SetFocalPoint/Position/ViewUp` from typed `CameraLock` (matches plane-lock precedent) |
| Markers | `MeasurementOverlay` in main renderer, opacity 0.5, occludable | Overlay renderer at `SetLayer(1)`, opacity 1.0, CPK-distinct palette, depth-occlusion-immune |
| Translucency | FXAA on, peeling off, classification flicker possible | Depth peeling with runtime fallback for AMD iGPU; `SetForceTranslucent(true)` on BS/HM to pin classification |
| Shutdown | `renderWindow_->Finalize()` while adapter still attached | `vtkWidget_->setRenderWindow(nullptr)` first; adapter destructor calls Finalize correctly; then drop the smart pointer |
| Diagnostic | Per-frame log, snapshot every 50 | EndEvent observer carries `source=` thread-local; per-render UDP correlatable to REST / timer / picker |

The two load-bearing pieces — Qt `eventFilter` as primary
camera-input path and `widget_->update()` as the only render verb
— recognise that `QVTKOpenGLNativeWidget` is a `QOpenGLWidget`
first and a VTK adapter second. The h5-reader's `QtAtomPicker`
proves the pattern works in this codebase. The rest of the table
follows.

### 6.1 Shutdown sequence (concrete)

`ReaderMainWindow::shutdown()` today (`ReaderMainWindow.cpp:449-471`)
stops timers, then calls `renderWindow_->Finalize()`. The order is
wrong: the adapter is still attached, and any queued paint event
between Finalize and widget destruction hits dead GL resources.
Right sequence:

```cpp
void shutdown() {
    if (restServer_) restServer_->stopListening();   // /shutdown fires here
    for (auto* t : findChildren<QTimer*>()) t->stop();
    if (vtkWidget_) {
        vtkWidget_->setRenderWindow(static_cast<vtkGenericOpenGLRenderWindow*>(nullptr));
        // adapter dtor calls RenderWindow->Finalize() internally
    }
    // smart pointer release happens via normal destruction now
}
```

The explicit `renderWindow_->Finalize()` is redundant and
out-of-order; `setRenderWindow(nullptr)` does the Finalize via the
adapter destructor (`QVTKRenderWindowAdapter.cxx:150-166`) under a
`makeCurrent()` and detaches in the right order.

## 7. Open questions for the design conversation

These belong in the next session, not this observation pass:

1. **Lock target source of truth**. AtomSelection is the natural
   producer (already has 4-atom cap, GeometryKind classification). Does
   the lock READ from AtomSelection (selection drives lock), or does
   the lock have its own pinned target independent of selection (lock
   stays even when selection changes)? The plane lock today reads
   from AtomSelection and releases when selection changes
   (`ReaderMainWindow.cpp:305-311`). Sensible default; user behaviour
   on the new locks will tell us if that's the right release semantic.
2. **User-delta persistence**. When the user rotates around a locked
   bond axis, does that rotation persist across lock target changes
   (pin a new atom but keep the user's rotation)? Across mode changes
   (atom-lock → bond-lock)? Across playback start/stop? Defaults
   matter; user feel will sort them out.
3. **Smoothing default**. Off by default (snap-per-frame), on by
   default (visible quality difference), or per-mode default (snap for
   atom-track, smooth for bond/dihedral)? Snap-by-default is the safe
   start — adds option, no surprise.
4. **Distance filter and lock target**. If the user has an atom lock
   AND a distance filter, does the filter centre on the lock target by
   default? Probably yes; user's mental model is "show me what's near
   this atom".
5. **Camera-relative bond up-axis convention**. For bond locks, screen-
   horizontal vs screen-vertical for the bond, default to which? For
   dihedrals, atom1 above or atom4 above when looking down 2-3? Newman
   projection convention is atom1 toward you (in front); we should
   probably match.
6. **Plane lock release on selection change** — the existing behaviour
   (drops lock when selection changes, `ReaderMainWindow.cpp:305-311`)
   is right for plane but may not be right for the new locks. A bond
   lock should probably survive the user picking a third atom (because
   the third atom might be the start of building a dihedral selection).

## 8. Citations index

Concrete pointers for the design conversation, grouped by purpose:

- **Working precedent for per-frame absolute camera write**: `src/app/MoleculeScene.cpp:282-415` (`lockCameraToSelectionPlane` + `applyCameraPlaneLock`), `src/app/PlaneFrameMath.h:1-79`.
- **The bouncing centroid-delta translation**: `src/app/MoleculeScene.cpp:446-497`.
- **Dihedral one-shot math worth lifting**: `src/app/MoleculeScene.cpp:624-694` (`focusCameraOnReveal`).
- **Interactor style install point**: `src/app/MoleculeScene.cpp:88-91`.
- **Paint-some-frames probe**: `src/app/MoleculeScene.cpp:501-528`, with bounds-cache workaround context at `src/app/MoleculeScene.cpp:434-445` and `src/app/MoleculeScene.cpp:547-553`.
- **Per-frame fan-out**: `src/app/MoleculeScene.cpp:530-555`, callers `src/app/ReaderMainWindow.cpp:194-199`.
- **Playback timer**: `src/app/QtPlaybackController.cpp:1-100` (the one legitimate QTimer).
- **Selection model**: `src/model/AtomSelection.h:59-138`, geometry-kind classification at `AtomSelection.h:56`.
- **Picker**: `src/app/QtAtomPicker.cpp:1-128` (double-click ray-cast, 2 Å tolerance).
- **Plane-lock UI wiring**: `src/app/ReaderMainWindow.cpp:179-186`, `ReaderMainWindow.cpp:532-538`, `ReaderMainWindow.cpp:662-680`.
- **Existing reveal flow** (dashboard → scene): `src/app/ReaderMainWindow.cpp:335-336`, `src/app/MoleculeScene.cpp:594-613`.
- **VTK-specific findings** (per research pass 2026-05-30):
  - `vtkCameraInterpolator` + `vtkQuaternionInterpolator` (SLERP'd orientation) for smooth interpolation
  - No per-atom visibility on `vtkOpenGLMoleculeMapper`; ghost arrays orphan bonds; subset display molecule is the path
  - `vtkInteractorStyleTrackballCamera` subclass + override `Rotate/Pan/Dolly` for user-delta composition
  - `vtkCamera::SetPosition/SetFocalPoint/SetViewUp` + `OrthogonalizeViewUp` for the absolute-write
  - `vtkMolecule::SetAtomPosition` calls `Modified()` per-call (verified against VTK source); no known VTK 9.x bug matches the end-of-trajectory drop, pointing at a race/transition rather than invalidation

## 9. What this doesn't cover

- Concrete API names and signatures for the new collaborator classes
  — design conversation deliverable, not this one.
- Implementation order or commit slicing — once the design lands.
- Migration of the existing 3-atom plane lock into the new
  `CameraLock` taxonomy — yes it should happen, but the existing API
  + UX is settled, no rush.
- The BS/HM stacking issue (`POLISH_BACKLOG.md` item 0) — separate
  rendering bug in a different overlay layer. Doesn't block the
  viewport refactor; possibly shares root cause with the paint-some-
  frames bug if both are translucency-pass-ordering, but that's a
  later investigation.
- Cross-platform behaviour (Windows MSVC, macOS Clang) — current
  code is platform-conditional only in the diagnostics layer; the
  viewport refactor shouldn't introduce new platform divergence.
  `vtkInteractorStyleTrackballCamera` subclassing is portable.
