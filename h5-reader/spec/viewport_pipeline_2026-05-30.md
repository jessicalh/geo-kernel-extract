# Viewport pipeline — settled design (2026-05-30)

Implementation contract for the h5-reader viewport refactor. Names each
stage of the rendering data flow, fixes each stage's single
responsibility, defines the interface between stages, and locks in the
lock-as-fit consolidation that the empirical numbers (see
`tests/scripts/HARNESS_BASELINE_TRANSFORM_2026-05-30.md`) demand.

This doc is the contract. The next session(s) implement against it. It
lives forever as the architectural record. Code sketches in this doc are
signatures and struct shapes; full implementations belong in the work
that lands the changes.

Sibling references (load if you are about to implement):
- Empirical motivation: `tests/scripts/HARNESS_BASELINE_2026-05-30.md`,
  `tests/scripts/HARNESS_BASELINE_TRANSFORM_2026-05-30.md`
- Observations: `notes/VIEWPORT_OBSERVATIONS_2026-05-30.md`
- Idiom source-of-truth: `qt6-cpp` skill — `architecture.md` (timer +
  signal/slot patterns), `crash-diagnosis.md` (CENSUS / ACONNECT /
  QPointer), `3d-vtk.md` (depth peeling, smart pointers, render-window
  lifecycle)
- Today's working precedent: `src/app/MoleculeScene.cpp:282-415`
  (`lockCameraToSelectionPlane` + `applyCameraPlaneLock`)
- Qt eventFilter precedent in this codebase: `src/app/QtAtomPicker.cpp:48`
- Decorator precedent: `src/model/TransformedConformation.cpp:1-279`

---

## 0. Headline design lens — the lock IS the fit

A "camera lock" in this codebase is not a separate primitive that
composes on top of `TransformedConformation`. The transform layer
(`TransformedConformation`) applies a per-frame rigid-body Kabsch fit on
a subset of atoms and rewrites positions. A camera lock does the same
thing — Kabsch fit on a subset of atoms — but uses the fit to compute an
absolute camera write each frame (`SetFocalPoint`/`SetPosition`/
`SetViewUp`).

The 3-atom plane lock at `MoleculeScene.cpp:282-415` is the working
precedent. The four other locks the user wants (atom-only, bond,
dihedral sight-down, backbone) are the same primitive on different atom
subsets, with an optional small camera-orientation override
(`OrientationPolicy`) on top.

The empirical evidence: `transform_plus_lock` median 39 px / max 795 px
vs `with_plane_lock` alone 1.5 px (`HARNESS_BASELINE_TRANSFORM_2026-05-30.md`)
shows the two layers interfere catastrophically at the tails. The same
fit running twice (once writing positions, once writing the camera) on
overlapping atom subsets disagrees on which frame is "the reference".
The design rejects two-layer composition and unifies: one fit per
`CameraMode`, applied either to positions, the camera, or both, with a
single shared reference frame.

---

## 1. Pipeline overview

### 1.1 Stage list

End-to-end data flow from a timer tick or a Qt mouse event to
user-visible pixels. Each stage has one responsibility, named below.
Crossing a stage boundary is an interface call, not a method on the same
object reaching into a different concern.

| # | Stage                       | Owner class                                 | Lives on    |
|---|-----------------------------|---------------------------------------------|-------------|
| 1 | Tick source                 | `QtPlaybackController`                       | GUI thread  |
| 2 | Position source             | `Conformation` (`TransformedConformation`)   | GUI thread  |
| 3 | Camera composer             | `CameraComposer` + `CameraMode` + `OrientationPolicy` | GUI thread |
| 4 | Scene mutator               | `MoleculeScene::setFrame`                    | GUI thread  |
| 5 | Render scheduler            | `MoleculeScene::requestRender`               | GUI thread  |
| 6 | Qt paint                    | `QVTKOpenGLNativeWidget::paintGL` (Qt-owned) | GUI thread  |
| 7 | VTK render                  | `vtkRenderWindow::Render` (VTK-owned)        | GUI thread  |
| 8 | End-event correlation       | `MoleculeScene` EndEvent observer            | GUI thread  |

There is **one** path for user-input-driven frames (the Qt eventFilter
described in Section 5) and **one** path for timer-driven frames; both
funnel into stages 5-8.

### 1.2 Timer-driven path (the playback case)

```
[Stage 1] QtPlaybackController::advance
        emit frameChanged(t)
              │
              ▼
[Stage 4] MoleculeScene::setFrame(t)
   ├── reads Stage 2: TransformedConformation::atomPosition(t, i) for i in atoms
   ├── consults Stage 3: CameraComposer::write(t, cameraMode, orientationPolicy, userDelta)
   ├── for each Overlay: overlay->setFrame(t)
   └── calls Stage 5: requestRender(RenderSource::Timer)
              │
              ▼
[Stage 5] coalesce: set lastRenderSource_, schedule one widget->update()
              │
              ▼
[Stage 6] Qt fires paintGL on the widget
              │
              ▼
[Stage 7] QVTKRenderWindowAdapter::paint() → iren->Render() → renderWindow->Render()
              │
              ▼
[Stage 8] EndEvent observer logs "render | source=timer | ms=<n> | frame=<t>"
```

### 1.3 User-input-driven path (camera gesture)

```
[Qt event] mouseMove on QVTKOpenGLNativeWidget
              │
              ▼
[Stage X] CameraInputFilter::eventFilter
   ├── classifies gesture (drag / wheel / shift+drag = pan)
   ├── computes free-helper math (Azimuth/Elevation/Pan/Dolly on a CameraGesture)
   ├── if a lock is active: mutates userDelta_ on the active CameraLock
   │                       (the lock then re-applies on the NEXT setFrame; the
   │                        gesture's intent is preserved)
   ├── else (free camera): writes camera state directly via free helpers
   └── calls Stage 5: requestRender(RenderSource::CameraInput)
              │
              ▼ (same Stage 5-8 as above)
```

### 1.4 What this is NOT

- **Not two render verbs.** `widget_->update()` is the only render verb
  in application code. Direct `renderWindow_->Render()` is removed from
  `MoleculeScene` (it remains an implementation detail of
  `QVTKRenderWindowAdapter::paint` per `QVTKRenderWindowAdapter.cxx:241-266`).
- **Not centroid-delta camera translation.** Stage 3 absolute-writes the
  camera every frame. The `lastCentroid_` / `haveLastCentroid_` /
  delta-translation code at `MoleculeScene.cpp:476-497` is removed.
- **Not two layers of fit (transform + lock).** The lock IS the fit. The
  `TransformedConformation` Kabsch path is reused by the camera composer
  for `CameraMode::Subset(atoms)` — same Kabsch, same reference, written
  to the camera instead of positions.

---

## 2. Per-stage detail

### 2.1 Stage 1 — Tick source

**Owner:** `QtPlaybackController` (already exists, no changes)
**File:** `src/app/QtPlaybackController.{h,cpp}`

**Responsibility:** the one application-owned `QTimer`. Emits
`frameChanged(int t)` on every tick, scrub, step, or setFrame. No
business logic; all consumers subscribe.

**Input contract:** none on the timer path; for direct user input it
takes `setFrame(int)`, `play()`, `pause()`, `stepForward()`,
`stepBackward()`, `setFps(int)` slots.

**Output contract:** `frameChanged(int t)` signal. `0 <= t <
frameCount`. Emitted on the GUI thread (the timer's thread).

**Discipline:** `CENSUS_REGISTER`, `ACONNECT` on the internal QTimer's
`timeout`. ASSERT_THREAD on every public mutator.

**Why unchanged:** the one-timer discipline is sound. The bouncing is
not a timer-source problem; it's a camera-writer problem. Section 4
(`MoleculeScene::setFrame`) is where the changes land.

---

### 2.2 Stage 2 — Position source

**Owner:** `Conformation` base; in practice `TransformedConformation`
wrapping `TrajectoryConformation` or `SingleConformation`
**File:** `src/model/Conformation.{h,cpp}`,
`src/model/TransformedConformation.{h,cpp}` (already exist)

**Responsibility:** answer "what is atom `i`'s position at frame `t`?"
through the polymorphic `atomPosition(frame, atomIdx)` seam. Optionally
apply a rigid-body Kabsch fit on a subset of atoms (Identity /
CenterCom / FitReference / FitSubset).

**Input contract:** `(std::size_t frame, std::size_t atomIdx)`; frame in
`[0, frameCount())`, atomIdx in `[0, protein->atomCount())`.

**Output contract:** `model::Vec3` in Ångströms, in **display space**
(after the transform). Pure: no side effects, idempotent for fixed
inputs and mode.

**Discipline:** `CENSUS_REGISTER` on the wrapper; ASSERT_THREAD on
`setMode`. `transformChanged` signal emitted exactly once per `setMode`
call so consumers (Stage 4) re-evaluate.

**Why this stage exists separately:** it's the decorator over the raw H5
seam. Per the no-pluggable-interfaces rule, the inner `Conformation*` is
a typed concrete pointer, not a `std::function<>` or factory. The
wrapper is itself a `Conformation` so consumers don't know whether the
transform is active — polymorphism does the rest.

**Key observation for the design downstream:** when `CameraMode::Subset`
is active, the `CameraComposer` (Stage 3) reads its reference positions
from the **inner** conformation (raw, pre-transform) so the lock and the
transform agree on what "frame 0" means. See Section 3.

---

### 2.3 Stage 3 — Camera composer

**Owner:** new `CameraComposer` class
**Files:** `src/app/CameraComposer.{h,cpp}` (new),
`src/app/CameraComposerMath.h` (new, header-only),
`src/app/FitTargetMath.h` (new, header-only)

**Responsibility:** given the current frame, the active `CameraMode`,
the active `OrientationPolicy`, and the accumulated user-gesture delta,
produce the absolute camera state (`focalPoint`, `position`, `viewUp`)
and write it via `SetFocalPoint` / `SetPosition` / `SetViewUp` +
`OrthogonalizeViewUp`. No delta accumulation. No `Render` call.

**Input contract:**

```cpp
// src/app/CameraComposer.h
namespace h5reader::app {

class CameraComposer final : public QObject {
    Q_OBJECT
public:
    CameraComposer(vtkSmartPointer<vtkRenderer> renderer,
                   const model::QtProtein* protein,
                   model::Conformation* conformation,
                   QObject* parent = nullptr);

    // Switch lock target. `mode` carries the atom subset and the lock
    // kind; `policy` the orientation override. Bumps the lock's
    // reference frame to `currentFrame` so subsequent applies are
    // continuous with the user's current view. Idempotent on equal
    // (mode, policy).
    void setMode(CameraMode mode, OrientationPolicy policy,
                 std::size_t currentFrame);

    // Absolute camera write for frame `t`. Returns false if the mode's
    // atom subset became degenerate at this frame (collinear, missing
    // atoms, etc); caller (Stage 4) decides whether to fall back to
    // CameraMode::Free or surface the failure.
    [[nodiscard]] bool write(std::size_t t);

    // User-gesture deltas (called from the Qt eventFilter on
    // QVTKOpenGLNativeWidget). The delta is accumulated INTO the
    // active CameraMode's state (e.g. for CameraMode::Atom, a rotate
    // gesture sets userAzimuth_/userElevation_ on the Atom variant).
    // On the next write(t), the per-frame fit is applied AND the user
    // delta is composed on top. For CameraMode::Free, the delta is
    // applied directly to the camera state.
    void applyGesture(const CameraGesture& g);

    // Read-only inspectors for REST / tests.
    CameraMode mode() const;
    OrientationPolicy policy() const;

signals:
    void modeChanged();

private:
    // ... (state including last applied transform, reference positions,
    //      sign continuity state, accumulated user delta per mode)
};

}  // namespace h5reader::app
```

**Output contract:** writes to `renderer_->GetActiveCamera()`'s
`FocalPoint`, `Position`, `ViewUp` (with `OrthogonalizeViewUp` after).
Bumps the camera's `Modified()` MTime, which propagates to the next
render. Returns `bool` for degenerate-input signaling.

**Discipline:** `CENSUS_REGISTER` in ctor; ASSERT_THREAD on `setMode`,
`write`, `applyGesture`. UDP log line per `setMode`:
`camera | mode=<name> | atoms=[i,j,...] | policy=<name>`. UDP log line
per `applyGesture` at DEBUG level (cheap to filter).

**VTK-side concerns:** writes via the absolute setters
(`vtkCamera::SetFocalPoint` etc, source at `vtkCamera.cxx:286-353`). Per
the source, `SetViewUp` normalises but does not orthogonalise — we call
`OrthogonalizeViewUp` explicitly (separate method per
`vtkCamera.cxx:561-569`). After all three setters,
`Renderer::ResetCameraClippingRange(double bounds[6])` is **not** called
here; that responsibility belongs to Stage 4 (which already computes
per-frame bounds in one pass through atom positions per
`MoleculeScene.cpp:447-459`).

**Diagnostic hook:** thread-local `lastRenderSource_` is owned by Stage
5 and set immediately before the `widget_->update()` enqueue. The
EndEvent observer (Stage 8) reads it and UDP-logs
`source=camera-write | mode=<name> | frame=<t>` if the render was
triggered by a setMode call. Inter-stage correlation is by `frame=`
field.

#### 2.3.1 `CameraMode` — the typed lock target

```cpp
// src/app/CameraMode.h
namespace h5reader::app {

// Sum type for the camera's lock target. Variants carry the atoms that
// define the lock; the composer reads positions for those atoms each
// frame and computes the absolute camera write.
//
// "The lock IS the fit": each variant corresponds to a fit subset
// (the atom list) AND a built-in orientation policy default. The
// OrientationPolicy parameter to CameraComposer::setMode is the
// optional override; passing OrientationPolicy::Default uses the
// variant's natural default.
struct CameraMode {
    enum class Kind {
        Free,                 // no lock; camera state owned by user gestures
        Atom,                 // 1 atom — F = atom_t, V/D inherited from gesture state
        Bond,                 // 2 atoms — F = midpoint, axis = atom2 - atom1
        Dihedral,             // 4 atoms — F = midpoint(b,c), axis = (c - b)
        Plane,                // 3 atoms — F = centroid, basis = (b-a)×(c-a)
        Subset                // N atoms — F = subset centroid, R from Kabsch fit on subset
    };
    Kind kind = Kind::Free;
    std::vector<std::size_t> atoms;   // .empty() iff Free; .size() implies kind

    // Identity (for cache key + log).
    bool operator==(const CameraMode& other) const = default;
};

// Helpers (free functions, header-only):
inline std::size_t ExpectedAtomCount(CameraMode::Kind k);  // 0/1/2/4/3/-1
inline const char* NameFor(CameraMode::Kind k);             // "free"/"atom"/...

// Construction helpers — typed constructors so callers don't pass the
// wrong cardinality. Each variant rejects (via assertion + log) inputs
// that don't match the kind's expected atom count.
inline CameraMode FreeMode();
inline CameraMode AtomMode(std::size_t a);
inline CameraMode BondMode(std::size_t a, std::size_t b);
inline CameraMode DihedralMode(std::size_t a, std::size_t b,
                                std::size_t c, std::size_t d);
inline CameraMode PlaneMode(std::size_t a, std::size_t b, std::size_t c);
inline CameraMode SubsetMode(std::vector<std::size_t> atoms);

}  // namespace h5reader::app
```

#### 2.3.2 `OrientationPolicy` — orientation override on top of the fit

```cpp
// src/app/OrientationPolicy.h
namespace h5reader::app {

// What ViewUp should be at each frame, GIVEN the per-frame fit's
// origin and (where defined) axis. Most variants are no-ops on the
// fit's natural orientation; non-default variants override the
// ViewUp + axis-of-sight that the fit would otherwise inherit from
// user gesture state.
struct OrientationPolicy {
    enum class Kind {
        Default,              // use the CameraMode's natural default
        Free,                 // V/D inherited from gesture state (no override)
        PerpendicularToBond,  // V = an axis perpendicular to the bond axis
        DownAxis,             // V is what the user had; D is the bond axis
                              // (Newman projection sight-down; used by Dihedral)
        PerpendicularToPlane, // V is the plane normal; D is in-plane
    };
    Kind kind = Kind::Default;
    // For DownAxis: which two atoms define the axis. For Dihedral mode
    // the natural default is (atoms[1], atoms[2]); the override can pick
    // any two of the four.
    std::array<std::size_t, 2> axisAtoms{0, 0};
};

inline const char* NameFor(OrientationPolicy::Kind k);
inline OrientationPolicy DefaultPolicy();
inline OrientationPolicy DownAxisPolicy(std::size_t a, std::size_t b);
inline OrientationPolicy PerpToBondPolicy();
inline OrientationPolicy PerpToPlanePolicy();

}  // namespace h5reader::app
```

#### 2.3.3 Per-`CameraMode` write semantics

| Kind     | Reference per frame | View direction default          | View-up default                   |
|----------|---------------------|---------------------------------|-----------------------------------|
| Free     | none                | inherited from gesture state    | inherited                         |
| Atom     | F = atom_t          | inherited from gesture state    | inherited                         |
| Bond     | F = midpoint(a,b)_t | perpendicular to bond axis      | bond axis (horizontal in screen)  |
| Dihedral | F = midpoint(b,c)_t | down (c-b) axis                 | (a-b) projected ⊥ (c-b)            |
| Plane    | F = centroid(a,b,c)_t | plane normal (sign-continuity) | (b-a) projected ⊥ normal           |
| Subset   | F = centroid(subset)_t | inherited (Kabsch removes rotation, so inherited makes sense) | inherited |

For Bond / Dihedral / Plane: the table is the *default* the variant
implies. The `OrientationPolicy` override can change which axis is
"down" and which is "up" (e.g. `DownAxis(a, d)` instead of
`DownAxis(b, c)` on a Dihedral).

For Subset: this is the bridge between the transform layer and the
camera. `CameraMode::Subset(backbone_atoms)` writes the camera such that
the backbone-stabilised view is in the same frame as the transformed
positions — the same Kabsch, applied to the camera instead of the
positions. This replaces the `transform_plus_lock` composition outright:
flip the transform to Identity and put the same atom subset in
`CameraMode::Subset` and the result is one-fit, one-reference, no tail
interference.

#### 2.3.4 `FitTargetMath.h` — pure-function math for unit testing

```cpp
// src/app/FitTargetMath.h
namespace h5reader::math {

// Compute the per-frame fit anchor for each CameraMode kind. Pure
// functions so each shape is unit-testable without a renderer or a
// loaded conformation. Each returns nullopt for degenerate inputs.

struct AtomAnchor {
    model::Vec3 focal;        // F: where the camera looks at
    std::optional<model::Vec3> axis;       // for Bond/Dihedral/Plane: the natural sight axis
    std::optional<math::PlaneFrame> frame; // for Plane: full orthonormal basis (re-uses PlaneFrameMath)
};

// Atom: F = positions[0]; axis = nullopt; frame = nullopt.
std::optional<AtomAnchor> ComputeAtomAnchor(
    const std::array<model::Vec3, 1>& positions);

// Bond: F = midpoint; axis = (b - a).normalized(); frame = nullopt.
std::optional<AtomAnchor> ComputeBondAnchor(
    const std::array<model::Vec3, 2>& positions);

// Dihedral: F = midpoint(b, c); axis = (c - b).normalized();
// frame = (axis, viewUp from (a-b) projected ⊥ axis, n × axis).
std::optional<AtomAnchor> ComputeDihedralAnchor(
    const std::array<model::Vec3, 4>& positions);

// Plane: thin wrapper around math::computePlaneFrame in
// PlaneFrameMath.h — exists so callers can route through the same
// ComputeXAnchor name in the CameraComposer.
std::optional<AtomAnchor> ComputePlaneAnchor(
    const std::array<model::Vec3, 3>& positions);

// Subset: Kabsch fit at the current frame against reference positions.
// reference must be the same atoms at lock-acquisition (referenceFrame_).
// Returns a Transform3D the composer uses to rotate the inherited
// camera direction into the subset's stabilised frame.
struct Transform3D { model::Mat3 R = model::Mat3::Identity(); model::Vec3 T = model::Vec3::Zero(); };
std::optional<Transform3D> ComputeSubsetTransform(
    const std::vector<model::Vec3>& current,
    const std::vector<model::Vec3>& reference);

}  // namespace h5reader::math
```

`PlaneFrameMath.h:chooseContinuousNormal` is the precedent for the
sign-continuity guard (lifted from the existing plane lock at
`MoleculeScene.cpp:390-395`). All variants that have a "natural" axis
(Plane, Dihedral) reuse it.

`ComputeSubsetTransform` is the same Kabsch math
`TransformedConformation::KabschFit` already implements
(`TransformedConformation.cpp:228-276`). The header here exposes the
free function so `CameraComposer` calls it directly without depending on
the decorator type; the existing static method on
`TransformedConformation` becomes a thin wrapper around this free
function in the same commit.

---

### 2.4 Stage 4 — Scene mutator

**Owner:** `MoleculeScene::setFrame` (existing; trimmed)
**File:** `src/app/MoleculeScene.{h,cpp}`

**Responsibility:** for frame `t`, push atom positions into the
`vtkMolecule`, compute per-frame bounds, ask Stage 3 to write the
camera, fan `setFrame` to overlays, set the clipping range, schedule one
render via Stage 5.

**Input contract:** `void setFrame(int t)` — the `frameChanged` slot.
`t` in `[0, frameCount)` (out-of-range is a guarded early-return per
`MoleculeScene.cpp:422`).

**Output contract:** per-frame VTK state mutations (molecule positions,
overlay backing data, clipping range, camera). Exactly one render
scheduled at the end via Stage 5. `currentFrame_ = t` BEFORE the
schedule (so the EndEvent observer in Stage 8 reads the correct frame).

**Discipline:** `ASSERT_THREAD(this)` first line (already there per
`MoleculeScene.cpp:419`). UDP DEBUG log per setFrame
(`frame <t> applied | <ms> ms` per `MoleculeScene.cpp:570-571`). The
every-50-frames snapshot stays at DEBUG (`MoleculeScene.cpp:572-591`) —
it caught the bounds-cache bug once; cheap to leave.

**VTK-side concerns:**

1. `vtkMolecule::SetAtomPosition` per atom; molecule's MTime bumps
   per call (`vtkMolecule.cxx:184-197`). Add an explicit
   `molecule_->GetPoints()->Modified()` AFTER the loop — this is the
   cheap probe from `notes/POLISH_BACKLOG.md` item 0 and
   `VIEWPORT_OBSERVATIONS §4`. The VTK source confirms `SetAtomPosition`
   calls `Points->SetPoint(...)` (which does **not** bump the points
   array's MTime explicitly), then `this->Modified()` on the molecule.
   `vtkOpenGLPolyDataMapper`'s VBO gate consults several MTimes; if a
   particular code path looks at the points' MTime, the gate misses.
   One line, header-confirmed pathology, no risk.

2. `mapper_->Modified()` at `MoleculeScene.cpp:528` becomes redundant
   after (1). Remove it. Keep the comment as a brief history note
   pointing here.

3. Per-frame bounds computation (already in the loop at
   `MoleculeScene.cpp:447-459`) stays. Pass the padded bounds to
   `Renderer::ResetCameraClippingRange(double bounds[6])` per
   `vtkRenderer.cxx:1062` — this is the form that bypasses the
   `vtkActor::GetBounds()` cache per `feedback_vtk_bounds_cache`.

4. The centroid-delta camera translation at `MoleculeScene.cpp:476-497`
   is **removed**. The `lastCentroid_` / `haveLastCentroid_` fields are
   removed. The camera is written by Stage 3.

5. The plane-lock-specific block at `MoleculeScene.cpp:463-471` is
   **removed**. Its logic moves into Stage 3 (`CameraMode::Plane`); the
   lock-as-fit consolidation in Section 3 below describes the bridge.

6. `lockCameraToSelectionPlane` / `clearCameraPlaneLock` /
   `applyCameraPlaneLock` / `cameraPlaneLock_` / the `CameraPlaneLock`
   struct (`MoleculeScene.h:148-154`,
   `MoleculeScene.cpp:282-415`) become thin compatibility shims atop
   `CameraComposer::setMode(PlaneMode(a, b, c))`. The public API stays
   for the REST surface (`/plane-lock/*`) and the toolbar action; the
   internals delegate. See Section 3 below for the exact shim shape.

7. `focusCameraOnReveal` (`MoleculeScene.cpp:624-694`) becomes a thin
   shim atop `CameraComposer::setMode(DihedralMode(...) +
   DownAxisPolicy(b, c))` when the binding is an `AtomTupleAnchor` with
   ≥ 4 atoms; otherwise it falls back to `CameraMode::Atom(centroid-of-
   atoms)`. The math at `MoleculeScene.cpp:658-695` moves into
   `FitTargetMath::ComputeDihedralAnchor`.

**Diagnostic hook:** at the end of `setFrame`, BEFORE the render
schedule, the UDP log line is:
`scene | frame=<t> | atoms=<N> | bounds=[xmin,xmax][ymin,ymax][zmin,zmax] | overlays=<count>`.
The render is then enqueued via Stage 5 with
`RenderSource::Timer`.

**Slimmed setFrame skeleton:**

```cpp
void MoleculeScene::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!molecule_ || !protein_ || !conformation_) return;
    if (t == currentFrame_) return;
    if (t < 0 || std::size_t(t) >= conformation_->frameCount()) return;

    QElapsedTimer timer; timer.start();

    // 1. Position push: one pass through atoms updates molecule
    //    positions and accumulates per-frame bounds.
    double bounds[6] = {+1e30, -1e30, +1e30, -1e30, +1e30, -1e30};
    PushAtomPositions(t, bounds);  // helper, was the loop in setFrame

    // 2. Modified bumps — molecule + points (the second is the probe).
    molecule_->Modified();
    molecule_->GetPoints()->Modified();   // PROBE; see §6 below

    // 3. Stage 3 writes the camera (absolute). No-op if mode==Free.
    cameraComposer_->write(std::size_t(t));

    // 4. Fan to overlays — each updates its own backing data.
    if (ribbon_)        ribbon_->setFrame(t);
    if (ringPolygons_)  ringPolygons_->setFrame(t);
    if (fieldGrid_)     fieldGrid_->setFrame(t);
    if (bfieldStream_)  bfieldStream_->setFrame(t);
    if (selection_)     selection_->setFrame(t);
    if (measurement_)   measurement_->setFrame(t);
    if (reveal_)        reveal_->setFrame(t);

    // 5. Clipping range from per-frame bounds (NOT the no-arg form;
    //    see feedback_vtk_bounds_cache).
    constexpr double pad = 5.0;
    double padded[6] = {bounds[0]-pad, bounds[1]+pad,
                        bounds[2]-pad, bounds[3]+pad,
                        bounds[4]-pad, bounds[5]+pad};
    renderer_->ResetCameraClippingRange(padded);

    // 6. Diagnostic UDP line + currentFrame_ update BEFORE the render
    //    schedule so the EndEvent observer (Stage 8) reads correctly.
    currentFrame_ = t;
    qCDebug(cScene) << "scene | frame=" << t
                    << "| dt_ms=" << timer.elapsed();

    // 7. Schedule one render via Stage 5.
    requestRender(RenderSource::Timer);
}
```

The setFrame becomes ~30 lines (was ~175). The shed weight is camera
math (Stage 3), the redundant `mapper_->Modified()` (folded into
`molecule_->GetPoints()->Modified()`), the centroid-delta block, and
the plane-lock special-case branch.

---

### 2.5 Stage 5 — Render scheduler

**Owner:** `MoleculeScene::requestRender` (new shape; replaces today's
`requestRender` + the four direct `renderWindow_->Render()` sites at
`MoleculeScene.cpp:237, 242, 337, 555`)
**File:** `src/app/MoleculeScene.{h,cpp}`

**Responsibility:** coalesce multiple per-tick render requests into one
draw. Set the thread-local `lastRenderSource_` so Stage 8 can correlate
the resulting render to its source. The render verb is
`vtkWidget_->update()` (per the `qt6-cpp` 3d-vtk reference and the
adapter source at `QVTKRenderWindowAdapter.cxx:241-266` — `paint()`
calls `iren->Render()`, which calls `renderWindow_->Render()` when
`EnableRender` is true).

**Input contract:**

```cpp
enum class RenderSource {
    Timer,        // QtPlaybackController::frameChanged → setFrame
    CameraInput,  // Qt eventFilter gesture
    Picker,       // double-click pick triggered an overlay update
    Rest,         // REST handler mutated state and needs a redraw
    Overlay,      // overlay visibility toggle from the UI
    Reveal,       // dashboard reveal binding
    External,     // explicit consumer request that doesn't fit above
};

void MoleculeScene::requestRender(RenderSource src);
```

**Output contract:** at most one `widget_->update()` per Qt event-loop
tick, regardless of how many `requestRender` calls fired. After the
update is enqueued, `renderPending_` is true; the queued lambda flips
it back to false before calling `widget_->update()`. No direct
`renderWindow_->Render()` in application code.

**Discipline:** `ASSERT_THREAD(this)` first line. The thread-local
`lastRenderSource_` is intentionally a non-atomic plain field — the
reader is GUI-thread-only; cross-thread render scheduling is not on the
path and would be a bug if it happened (caught by the ASSERT). UDP log
at DEBUG: `render-scheduled | source=<name>`.

```cpp
void MoleculeScene::requestRender(RenderSource src) {
    ASSERT_THREAD(this);
    lastRenderSource_ = src;  // EndEvent observer (Stage 8) reads this
    if (renderPending_ || !vtkWidget_) return;
    renderPending_ = true;
    QPointer<MoleculeScene> self(this);
    QMetaObject::invokeMethod(this, [self]() {
        if (!self) return;
        self->renderPending_ = false;
        if (self->vtkWidget_) self->vtkWidget_->update();
    }, Qt::QueuedConnection);
}
```

No `QTimer::singleShot(0, …)` — `QMetaObject::invokeMethod(...,
Qt::QueuedConnection)` is the queued-tick primitive and does NOT count
as a timer per the qt6-cpp skill's timer-discipline section
(`architecture.md` §4).

**Why this shape:** the qt6-cpp 3d-vtk reference makes
`vtkWidget_->update()` the only proper render verb in a Qt+VTK app.
`renderWindow_->Render()` skips the adapter's FBO management and the
paint-event chain. With one scheduler funnelling everything, Stage 8 has
a single correlation point.

**Stores `vtkWidget_` reference:** `MoleculeScene` today does not hold
`QVTKOpenGLNativeWidget*`. The refactor passes the widget into the
ctor (replacing the bare `vtkSmartPointer<vtkGenericOpenGLRenderWindow>`
parameter at `MoleculeScene.h:80`) so the scheduler can call
`widget_->update()`. `MoleculeScene` already obtains its
`vtkGenericOpenGLRenderWindow` from `vtkWidget_->renderWindow()`.

---

### 2.6 Stage 6 — Qt paint

**Owner:** `QVTKOpenGLNativeWidget::paintGL` (VTK-owned, unchanged)
**Source:** `/home/jessica/builds/VTK-9.5.2/GUISupport/Qt/QVTKOpenGLNativeWidget.cxx:274-305`

**Responsibility:** Qt fires `paintGL` when the queued
`widget_->update()` from Stage 5 drains. `paintGL` resets GL state to
VTK's expectations and delegates to the adapter's `paint()` (source at
`QVTKRenderWindowAdapter.cxx:241-266`).

**Input contract:** none from application code; Qt's compositor handles
this stage.

**Output contract:** `adapter->paint()` returns; widget is blitted to
its default framebuffer.

**Discipline:** none from application code. Per the qt6-cpp 3d-vtk
reference, this is a contract VTK owns; we trust it.

**Why this stage is named:** to make explicit that we don't bypass it.
`renderWindow_->Render()` skips the paint-event chain and the FBO blit;
`widget_->update()` goes through it. The acceptance gates' "force a
fresh render then read" path (`vtkWindowToImageFilter` with
`ShouldRerenderOn`) lives in the screenshot REST handler and uses VTK's
own render path; the per-frame app path uses Qt's.

---

### 2.7 Stage 7 — VTK render

**Owner:** `vtkRenderWindow::Render` (VTK-owned)
**Source path:** `QVTKRenderWindowAdapter::paint` → `iren->Render`
(`vtkRenderWindowInteractor.cxx:168-177`) → `RenderWindow->Render`

**Responsibility:** walk the renderer collection, run each renderer's
pipeline (atom imposters, then overlays, then markers on the overlay
layer per Section 6), composite, frame.

**Input contract:** all `vtkSmartPointer`s on the renderer's actor
collection are non-null; mappers' inputs are valid `vtkDataObject`s;
camera state has been written (Stage 3) and clipping range has been
reset (Stage 4).

**Output contract:** GL framebuffer for the widget is updated;
`vtkCommand::EndEvent` fires on the render window after the frame is
done.

**Discipline:** none from application code. VTK owns this stage; we
respect the contract.

**Important constraint:** `interactor_->EnableRenderOff()` (called once
at MoleculeScene construction; see Section 5.2 below) gates
`vtkRenderWindowInteractor::Render` per `vtkRenderWindowInteractor.cxx:170`.
With this gate active, the stock trackball's `rwi->Render()` calls
(five sites in `vtkInteractorStyleTrackballCamera.cxx` at lines 272,
297, 348, 393, 451) become no-ops at the render level — the trackball
can still mutate camera state for events the Qt eventFilter doesn't
intercept, but it cannot trigger renders. The render schedule remains
ours.

**But:** `iren->Render()` is also called by the adapter's `paint()`
(`QVTKRenderWindowAdapter.cxx:256-263`). With `EnableRender=false`,
`iren->Render()` becomes a no-op — but the adapter's paint flow also
falls back to `renderWindow_->Render()` directly when there's no
interactor (`QVTKRenderWindowAdapter.cxx:260-263`). The flow with
interactor + EnableRender=false: adapter calls `iren->Render()`, the
gate blocks the inner `renderWindow_->Render()`, but `iren->Render()`
still fires `RenderEvent` (per `vtkRenderWindowInteractor.cxx:176`); the
adapter's paint chain then completes and the widget is blitted.

The catch: we need the paint chain to call `renderWindow_->Render()`
exactly once per `widget_->update()`. Reading the VTK source carefully:
when `iren->Render()` no-ops, the adapter's `paint()` does NOT
separately call `renderWindow_->Render()` — it expects the interactor's
Render to have done it. So with `EnableRenderOff()`, the
adapter+interactor combination produces no actual render on `paintGL`.

**Resolution:** install the interactor (per the current pattern at
`MoleculeScene.cpp:88-91`) so the adapter's paint chain HAS an
interactor object, but call `iren->EnableRenderOn()` — keep the
interactor's Render path active. Disable RENDERS coming from the
trackball by a different mechanism: a custom interactor style that
inherits from `vtkInteractorStyleTrackballCamera` but overrides
`Rotate/Pan/Dolly/Spin/EnvironmentRotate` to no-op. With the eventFilter
(Section 5) intercepting mouse events before VTK sees them, the
trackball's manipulators never fire — so the override is belt-and-
suspenders for VTK-internal triggers (3D mouse, touch, etc.) that the
eventFilter doesn't cover.

```cpp
// src/app/QuietTrackballStyle.h
namespace h5reader::app {

// Subclass of vtkInteractorStyleTrackballCamera that no-ops every
// camera manipulator. Installed on the render window's interactor so
// the adapter has a valid style (paint chain expects one), but
// 3D-mouse / touch / unforeseen-input paths never mutate the camera
// behind the application's back.
//
// The Qt eventFilter on QVTKOpenGLNativeWidget intercepts mouse/wheel
// events for the camera-input path (Section 5). This subclass exists
// for the events the filter does NOT intercept.
class QuietTrackballStyle : public vtkInteractorStyleTrackballCamera {
public:
    static QuietTrackballStyle* New();
    vtkTypeMacro(QuietTrackballStyle, vtkInteractorStyleTrackballCamera);

    void Rotate() override {}
    void Spin() override {}
    void Pan() override {}
    void Dolly() override {}
    void EnvironmentRotate() override {}
    // OnMouseWheelForward/Backward also no-op (they call Dolly internally
    // by default).
    void OnMouseWheelForward() override {}
    void OnMouseWheelBackward() override {}

protected:
    QuietTrackballStyle() = default;
    ~QuietTrackballStyle() override = default;
};

}  // namespace h5reader::app
```

This subclass is the **only** new VTK subclass introduced by the
refactor. It's a narrow override of behaviour that's overtly the wrong
default for this app (paint flow expects an interactor style; the stock
style fires renders we can't gate cleanly). Per the
no-pluggable-interfaces rule, it's a concrete typed class installed
once at construction, not a factory.

`AutoAdjustCameraClippingRangeOff()` is **also** called on the style
(belt-and-suspenders against the no-arg `ResetCameraClippingRange()` at
`vtkInteractorStyleTrackballCamera.cxx:264`) — even though all camera
manipulators are overridden to no-op, this protects against future VTK
code paths that might call the auto-adjust outside the manipulator
methods. One line; effectively free.

---

### 2.8 Stage 8 — End-event correlation

**Owner:** `MoleculeScene` EndEvent observer (already exists at
`MoleculeScene.cpp:101-111`; extended)
**File:** `src/app/MoleculeScene.{h,cpp}`

**Responsibility:** after each `vtkRenderWindow` render
(`vtkCommand::EndEvent` fires per the GenericOpenGL window's
`Frame()` at `vtkGenericOpenGLRenderWindow.cxx:73-77`), log a structured
UDP line correlating the render to its source and timing.

**Input contract:** reads `MoleculeScene::lastRenderSource_` and
`MoleculeScene::currentFrame_` (both GUI-thread only; ASSERT_THREAD
guarded).

**Output contract:** one UDP log line per render:
```
render | source=<timer|camera-input|picker|rest|overlay|reveal|external> |
         frame=<t> | ms=<elapsed> | mode=<camera mode name>
```

**Discipline:** the observer fires on the GUI thread (VTK invokes
observers from the render call). UDP-only; not stderr (the path is hot:
every render goes through here). Per `feedback_udp_logging` and
`feedback_capture_at_the_boundary`.

**Why this stage is named:** every other diagnostic in the system
ultimately routes back to "the render happened, and here's why". The
EndEvent is the natural boundary at which to capture that. Stage 5
writes a hint; Stage 8 reads it.

```cpp
// In MoleculeScene ctor (replaces the lambda at MoleculeScene.cpp:101-111):
auto cb = vtkSmartPointer<vtkCallbackCommand>::New();
cb->SetClientData(this);
cb->SetCallback([](vtkObject*, unsigned long, void* cd, void*) {
    auto* self = static_cast<MoleculeScene*>(cd);
    if (!self || !self->renderer_) return;
    const double ms = self->renderer_->GetLastRenderTimeInSeconds() * 1000.0;
    const char* src = NameFor(self->lastRenderSource_);
    const char* mode = NameFor(self->cameraComposer_->mode().kind);
    qCInfo(cScene).noquote() << "render | source=" << src
                              << "| frame=" << self->currentFrame_
                              << "| ms=" << QString::number(ms, 'f', 1)
                              << "| mode=" << mode;
});
renderWindow_->AddObserver(vtkCommand::EndEvent, cb);
```

---

## 3. The lock-as-fit consolidation

The core architectural move. Today the 3-atom plane lock at
`MoleculeScene.cpp:282-415` is a working precedent — sub-pixel-floor
camera writes per frame from an absolute computation off a typed lock
state. The other locks (atom / bond / dihedral / backbone) need to land
in the same shape; the `transform_plus_lock` interference shows what
happens when they don't.

### 3.1 The existing plane lock as the template

`lockCameraToSelectionPlane(atoms)` captures gesture state into a
`CameraPlaneLock` struct (atoms, localViewUp, normalSign, distance,
lastDirection); `applyCameraPlaneLock(frame)` rebuilds the plane basis
from the current frame's positions, applies the continuity guard, and
writes the camera absolutely. `setFrame` calls `applyCameraPlaneLock`
before the centroid-delta block and short-circuits the latter when the
lock is active.

In the new design, all of that math moves into `CameraComposer` +
`FitTargetMath`. The `MoleculeScene` keeps the public API for backward
compatibility with the REST surface and the toolbar:

```cpp
// src/app/MoleculeScene.h — shims, signatures unchanged
bool lockCameraToSelectionPlane(const std::vector<std::size_t>& atoms);
void clearCameraPlaneLock();
bool isCameraPlaneLocked() const;
std::vector<std::size_t> cameraPlaneLockAtoms() const;
```

Implementation becomes:

```cpp
// src/app/MoleculeScene.cpp — shim implementations
bool MoleculeScene::lockCameraToSelectionPlane(const std::vector<std::size_t>& atoms) {
    ASSERT_THREAD(this);
    if (atoms.size() != 3) return false;
    const std::size_t frame = currentFrame_ >= 0 ? std::size_t(currentFrame_) : 0;
    cameraComposer_->setMode(PlaneMode(atoms[0], atoms[1], atoms[2]),
                              DefaultPolicy(),
                              frame);
    requestRender(RenderSource::External);
    emit cameraPlaneLockChanged(true);
    return true;
}

void MoleculeScene::clearCameraPlaneLock() {
    ASSERT_THREAD(this);
    if (cameraComposer_->mode().kind != CameraMode::Kind::Plane) return;
    cameraComposer_->setMode(FreeMode(), DefaultPolicy(),
                              std::size_t(currentFrame_));
    requestRender(RenderSource::External);
    emit cameraPlaneLockChanged(false);
}

bool MoleculeScene::isCameraPlaneLocked() const {
    return cameraComposer_->mode().kind == CameraMode::Kind::Plane;
}

std::vector<std::size_t> MoleculeScene::cameraPlaneLockAtoms() const {
    const auto& m = cameraComposer_->mode();
    if (m.kind != CameraMode::Kind::Plane) return {};
    return m.atoms;
}
```

The signatures don't change; the toolbar action at
`ReaderMainWindow.cpp:722-740` and the REST endpoints
(`RestServer.cpp` `/plane-lock/*`) keep working with no edits.

### 3.2 The dihedral one-shot

`MoleculeScene::focusCameraOnReveal` at
`MoleculeScene.cpp:624-694` is a one-shot dihedral sight-down today
(camera is set once when a reveal binding fires, then the next setFrame
delta-translates it). The math at lines 658-695 is correct — `axis = (c
- b).normalized()`, `F = (b + c) * 0.5`, ViewUp preserved-then-
orthogonalised — but it's not re-applied per frame.

In the new design this becomes:

```cpp
void MoleculeScene::focusCameraOnReveal(const model::SignalBinding& binding,
                                         const std::vector<std::size_t>& atoms,
                                         int frame) {
    if (!renderer_ || !protein_ || !conformation_ || atoms.empty()) return;

    const auto* tuple = std::get_if<model::AtomTupleAnchor>(&binding.anchor);
    const bool canSightDown = tuple && tuple->atoms.size() >= 4;

    if (canSightDown) {
        // The exact math at MoleculeScene.cpp:658-695 lives in
        // FitTargetMath::ComputeDihedralAnchor now; this call wires
        // the lock so subsequent setFrame() re-apply the sight-down
        // every frame, not just this one.
        cameraComposer_->setMode(
            DihedralMode(tuple->atoms[0], tuple->atoms[1],
                          tuple->atoms[2], tuple->atoms[3]),
            DownAxisPolicy(tuple->atoms[1], tuple->atoms[2]),
            std::size_t(frame));
    } else {
        // Single-atom or arbitrary-tuple reveal: F = centroid of atoms,
        // V/D inherited from current camera. Atom lock; if the centroid
        // is what we want, package the atoms as Subset with the
        // centroid-only fit (Subset's Kabsch on N<3 falls back to
        // identity rotation with translation = centroid_t - centroid_ref,
        // i.e. pure follow-the-centroid).
        cameraComposer_->setMode(
            SubsetMode(atoms),
            DefaultPolicy(),
            std::size_t(frame));
    }
    // No requestRender here: revealBinding caller already issues one.
}
```

Per-frame re-application is automatic: Stage 4's `setFrame` calls
`cameraComposer_->write(t)` on every frame, which applies the dihedral
math (or any other mode) absolute-write for the current frame's atom
positions. The reveal is no longer a one-shot; it's a sustained lock
that releases when the reveal clears.

### 3.3 Per-frame composition with user-gesture delta

The user can drag-rotate while a lock is active. The new design holds
the gesture's accumulated delta as state on the active `CameraMode`
variant (or on the `CameraComposer` itself), and re-applies it AFTER the
per-frame fit. So:

- Each frame, `CameraComposer::write(t)`:
  1. Reads the relevant atoms' positions at `t` through Stage 2
     (`Conformation::atomPosition`).
  2. Calls the appropriate `FitTargetMath::Compute*Anchor` and
     `ComputeSubsetTransform` to derive the fit's natural camera state
     for this frame (focal point, view direction, view-up, distance).
  3. Applies the accumulated user-gesture delta on top
     (`Azimuth/Elevation/Roll/Dolly` of the camera around the fit's
     focal point). The delta is in the fit's local frame, so a rotation
     applied at frame `t-1` produces the same rotation relative to the
     atoms at frame `t`.
  4. Writes `SetFocalPoint/SetPosition/SetViewUp + OrthogonalizeViewUp`.

- When the user issues a gesture (Stage 5 path), `applyGesture` updates
  the delta on the active mode variant. The next setFrame's
  `cameraComposer_->write(t)` re-applies fit + delta together. The
  gesture is preserved across playback frames.

For the existing plane lock, the gesture state is approximately the
existing `localViewUp + distance` fields on `CameraPlaneLock`
(`MoleculeScene.h:148-154`). Generalised to the other modes: each
variant has a small gesture-state struct (e.g. `AtomLockGestureState`
holds `Azimuth/Elevation/Roll/Distance`; `DihedralLockGestureState`
holds those plus an `axisDirectionSign` for the continuity guard).

### 3.4 Why this resolves the `transform_plus_lock` interference

The interference was: `TransformedConformation` does Kabsch fit on
backbone, writes positions. The plane lock then re-fit from the
already-transformed atoms 12, 13, 14 — and the lock captured its
reference at lock-acquisition time, which was AFTER the transform's
first frame. The two layers' "frame 0" disagreed; the continuity guard
at `MoleculeScene.cpp:390-395` fired spuriously when the transformed
basis crossed the lock's stored reference direction.

The new design unifies: the user picks ONE policy. Either:

- **Transform layer Identity, CameraMode::Plane(12,13,14)** — same as
  today's working `with_plane_lock` case (median 1.5 px). The camera is
  written from a fit at the same atoms each frame.
- **Transform layer Identity, CameraMode::Subset(backbone_atoms)** —
  same Kabsch as today's `fit_subset` transform layer, but applied to
  the camera. The marker rides the stabilised backbone-frame view
  without a separate transform writing positions. The composition is
  one fit.
- **Transform layer FitSubset(backbone), CameraMode::Free** — today's
  `transform_only` (median 67 px). Positions written, camera not. Best
  current option for users who want the protein "rest-framed" on
  screen.
- **Transform layer FitSubset(backbone), CameraMode::Atom(14)** — both
  layers active, but the camera lock is just translation (Atom mode is
  F = atom_t, no rotation), so there's no Kabsch-vs-Kabsch fight.
  Single-atom centring on top of backbone-stabilisation is the
  best-of-both-worlds combination this design enables that today's
  combination cannot.

The acceptance gate for `transform_plus_lock` (median ≤ 5 px AND p95 ≤
10 px AND max ≤ 30 px per `HARNESS_BASELINE_TRANSFORM_2026-05-30.md`)
is met by either of:

- The new design's `FitSubset + Atom(14)` combination — the camera lock
  is translation-only (no rotation), so it composes additively with the
  transform's rotation.
- The new design's `Identity + Subset(backbone)` combination — one
  Kabsch, one reference, one write.

The harness experiment for the acceptance gate uses the second option
(no transform; camera does the work); the first option exists for users
who want both layers visible on the screen (e.g. backbone-rest-frame for
the protein, atom-centred for the marker).

---

## 4. The Qt+VTK boundary

The mouse-handling path. The h5-reader precedent
(`QtAtomPicker.cpp:48`) installs an `eventFilter` on
`QVTKOpenGLNativeWidget` and intercepts events before VTK sees them. The
camera-input path follows the same shape, extended to mouse-move /
wheel / press / release for camera gestures.

### 4.1 `CameraInputFilter` — the eventFilter

```cpp
// src/app/CameraInputFilter.h
namespace h5reader::app {

class CameraInputFilter final : public QObject {
    Q_OBJECT
public:
    CameraInputFilter(QVTKOpenGLNativeWidget* widget,
                      MoleculeScene* scene,
                      CameraComposer* composer,
                      QObject* parent = nullptr);
    ~CameraInputFilter() override;

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

private:
    void handleMouseDown(QMouseEvent* me);
    void handleMouseMove(QMouseEvent* me);
    void handleMouseUp(QMouseEvent* me);
    void handleWheel(QWheelEvent* we);

    QPointer<QVTKOpenGLNativeWidget> widget_;
    QPointer<MoleculeScene>          scene_;
    QPointer<CameraComposer>         composer_;

    enum class Gesture { None, Rotate, Pan, Roll };
    Gesture activeGesture_ = Gesture::None;
    QPointF lastPos_;
    Qt::KeyboardModifiers lastMods_ = Qt::NoModifier;
};

}  // namespace h5reader::app
```

**Event interception policy:** the filter catches and returns `true`
(consuming the event so VTK does not see it) for:

- `MouseButtonPress` with left/middle/right button, when not over the
  QVTKOpenGLNativeWidget's child widgets (overlay docks etc.). On press,
  records `lastPos_`, `lastMods_`, starts a `Gesture::Rotate` (left),
  `Pan` (middle), `Roll` (right). Shift modifier promotes left to Pan.
- `MouseMove` when `activeGesture_ != None`. Computes delta against
  `lastPos_`, calls `composer_->applyGesture(...)`, calls
  `scene_->requestRender(RenderSource::CameraInput)`. Returns true.
- `MouseButtonRelease` when `activeGesture_ != None`. Sets
  `activeGesture_ = None`. Returns true.
- `Wheel` events. Computes dolly factor from `we->angleDelta().y()`,
  calls `composer_->applyGesture(...)`. Returns true.
- `MouseButtonDblClick` is **NOT intercepted** — falls through to
  `QtAtomPicker` (the dbl-click filter at `QtAtomPicker.cpp:55-62`).
  The picker installs its own eventFilter; Qt event filters run in
  reverse install order, so installing `CameraInputFilter` AFTER
  `QtAtomPicker` lets the picker see double-clicks first.

Event filters that the filter does NOT intercept:

- `MouseButtonDblClick` — picker owns it.
- `KeyPress/KeyRelease` — Qt shortcuts / dock focus / etc.
- `Enter/Leave/HoverMove` — Qt cursor management.
- `MouseButtonPress` while the user is dragging a dock divider or
  resizing the main window — the press happens on a different widget,
  not on `vtkWidget_`.

### 4.2 Trackball gating

At `MoleculeScene` construction (replacing the stock trackball install
at `MoleculeScene.cpp:88-91`):

```cpp
if (auto* iren = renderWindow_->GetInteractor()) {
    vtkNew<QuietTrackballStyle> style;     // no-op overrides per §2.7
    style->AutoAdjustCameraClippingRangeOff();  // per §1.2 lesson
    iren->SetInteractorStyle(style);
    // iren->EnableRender stays ON (default). The adapter's paint chain
    // needs iren->Render() to NOT be a no-op so the render actually
    // fires when paintGL drains. QuietTrackballStyle no-ops the
    // mutators, so trackball-internal Render calls (at .cxx:272 et al)
    // simply don't happen.
}
```

This composition replaces the trio of "EnableRenderOff + AutoAdjust...Off
+ keep stock trackball" that the observations doc considered. The
practical issue: `EnableRenderOff()` also no-ops the adapter's paint
flow because the adapter's `paint()` only calls `iren->Render()` (no
fallback when an interactor is present, per
`QVTKRenderWindowAdapter.cxx:256-263`). With `QuietTrackballStyle`,
EnableRender stays on, paintGL renders normally, and the trackball
simply does nothing on its own initiative.

### 4.3 Wiring in `ReaderMainWindow`

```cpp
// src/app/ReaderMainWindow.cpp — added after picker_ construction
cameraInputFilter_ = new CameraInputFilter(vtkWidget_, scene_,
                                            scene_->cameraComposer(),
                                            this);
// Picker is installed first (existing line); CameraInputFilter is
// installed after, so Qt's filter chain calls CameraInputFilter
// FIRST. CameraInputFilter passes MouseButtonDblClick through (does
// not intercept), and the picker handles it.
```

### 4.4 Shutdown reorder

Today `ReaderMainWindow::shutdown` at
`ReaderMainWindow.cpp:509-531` stops timers then calls
`renderWindow_->Finalize()`. The order is wrong: the
`QVTKRenderWindowAdapter` is still attached, and any queued paint event
between Finalize and widget destruction touches dead GL resources.

Per the `qt6-cpp` 3d-vtk reference and
`QVTKOpenGLNativeWidget.cxx:90-135` + `QVTKRenderWindowAdapter.cxx:150-166`
(the adapter's destructor calls `RenderWindow->Finalize()` under
`makeCurrent()`, in the right order):

```cpp
void ReaderMainWindow::shutdown() {
    ASSERT_THREAD(this);
    if (shutdownDone_) return;
    shutdownDone_ = true;

    qCInfo(cWindow).noquote() << "shutdown entered";

    // 1. Stop the REST server synchronously. The /shutdown endpoint
    //    fires from a request, so we want the server to drain before
    //    we tear timers down.
    if (restServer_) restServer_->stopListening();

    // 2. Stop every timer owned by us or our children.
    for (auto* timer : findChildren<QTimer*>()) {
        if (timer->isActive()) timer->stop();
    }

    // 3. Detach the render window from the widget BEFORE dropping the
    //    smart pointer. The widget's setRenderWindow(nullptr) calls
    //    makeCurrent() and resets the adapter; the adapter's
    //    destructor (per QVTKRenderWindowAdapter.cxx:150-166) calls
    //    renderWindow->Finalize() in the right context.
    if (vtkWidget_) {
        vtkWidget_->setRenderWindow(static_cast<vtkGenericOpenGLRenderWindow*>(nullptr));
    }

    // 4. The renderWindow_ smart pointer drops as part of normal
    //    destruction now; Finalize already ran in step 3.

    qCInfo(cWindow).noquote() << "shutdown done";
}
```

The explicit `renderWindow_->Finalize()` at the old
`ReaderMainWindow.cpp:527` is removed — `setRenderWindow(nullptr)` does
it via the adapter in the right context. The `RestServer::stopListening`
synchronous call before timers stop is the second part of the reorder
(today `/shutdown` may race with timer shutdown if it fires from a
request after the playback timer has stopped).

---

## 5. The overlay-layer renderer

Markers (today: `MeasurementOverlay`'s ≤4 spheres + connecting line;
potentially the picker's selection sphere if it gains one) live in a
second `vtkRenderer` at `SetLayer(1)`, sharing the active camera with
the main renderer. The harness's marker-blob analysis needs the markers
to be visible regardless of depth occlusion — at the bottom of the
trajectory's drift sweep, the marker can be behind the protein interior;
the blob detector loses signal when the marker is occluded by imposter
spheres.

### 5.1 Two-renderer composition

```cpp
// src/app/MoleculeScene.cpp — ctor changes
renderer_ = vtkSmartPointer<vtkRenderer>::New();
renderer_->SetLayer(0);              // main scene layer
renderer_->SetBackground(1.0, 1.0, 1.0);
renderer_->SetUseFXAA(true);
// SetUseDepthPeeling stays 0 for now per the original comment
// (works on this hardware; flip if translucency artifacts on AMD); see §6.

overlayRenderer_ = vtkSmartPointer<vtkRenderer>::New();
overlayRenderer_->SetLayer(1);       // overlay layer — depth-occlusion-immune
overlayRenderer_->SetActiveCamera(renderer_->GetActiveCamera());  // shared camera
overlayRenderer_->SetBackground(0.0, 0.0, 0.0);
overlayRenderer_->SetBackgroundAlpha(0.0);  // transparent — composites over main
overlayRenderer_->SetUseFXAA(false);  // crisp marker edges for blob detection
overlayRenderer_->SetUseDepthPeeling(0);
overlayRenderer_->InteractiveOff();   // do not receive interaction events

renderWindow_->SetNumberOfLayers(2);
renderWindow_->AddRenderer(renderer_);
renderWindow_->AddRenderer(overlayRenderer_);
```

Per VTK's layered composition: each renderer paints in order
(layer 0, then layer 1) with depth reset between layers. Actors on
layer 1 paint on top of layer 0 regardless of their world-space z. The
shared camera means screen-space positions of markers and atoms remain
aligned.

### 5.2 Which actors move to overlay layer

Today's main-renderer actors that stay on layer 0:
- `actor_` — the `vtkMolecule` actor
- `ribbon_`'s ribbon actor
- `ringPolygons_`'s ring polygon actors + normal arrows
- `fieldGrid_`'s BS/HM isosurfaces (Butterfly)
- `bfieldStream_`'s streamlines
- `selection_`'s old single-atom yellow highlight (dormant per
  `MoleculeScene.cpp:174`)

Actors that move to layer 1:
- `measurement_` — the ≤4 spheres + connecting line. Instrument-mode
  markers must be findable by the harness.
- `reveal_`'s sphere highlights + connecting lines + tensor ellipsoid
  (the latter is opacity 0.4 per `SceneRevealOverlay.cpp:31`). The
  tensor glyph is a "should I see this" piece of UI; layered keeps it
  visible.

The picker (`QtAtomPicker`) doesn't render anything today, so no
actor moves.

### 5.3 Wiring the two-renderer constructor

`MeasurementOverlay::MeasurementOverlay` (and `SceneRevealOverlay`) take
both renderers in their constructor:

```cpp
// src/app/MeasurementOverlay.h — updated signature
MeasurementOverlay(vtkSmartPointer<vtkRenderer> mainRenderer,
                   vtkSmartPointer<vtkRenderer> overlayRenderer,
                   vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
                   QObject* parent = nullptr);
```

Internally the overlay adds its actors to `overlayRenderer_`. `mainRenderer`
remains in the signature for symmetry with overlays that paint on
layer 0 (and for a future opt-in that lets a marker overlay paint on
either layer); not used by the current implementation.

### 5.4 Depth peeling — runtime detection per
`project_h5reader_target_hardware`

The main renderer's depth peeling stays off for now (per
`MoleculeScene.cpp:80-81`'s comment "translucency artifacts on AMD
hardware"). If/when the BS/HM stacking bug at
`POLISH_BACKLOG.md` item 0 is investigated, the fix surface includes:

1. Try `SetForceTranslucent(true)` on the BS/HM isosurface actors (the
   second cheap probe from Section 6 below).
2. If that doesn't resolve it, turn depth peeling on for the main
   renderer behind a runtime hardware check:
   ```cpp
   const QString glRenderer = QString::fromUtf8(renderWindow_->ReportCapabilities());
   const bool useDepthPeeling = !glRenderer.contains("AMD", Qt::CaseInsensitive)
                                && !glRenderer.contains("Radeon", Qt::CaseInsensitive);
   if (useDepthPeeling) {
       renderer_->SetUseDepthPeeling(1);
       renderer_->SetMaximumNumberOfPeels(8);
       renderer_->SetOcclusionRatio(0.0);
       renderWindow_->SetAlphaBitPlanes(1);   // already on
       renderWindow_->SetMultiSamples(0);     // already on
   }
   ```

The overlay renderer never uses depth peeling — markers are solid
(opacity 1.0 in instrument mode); peeling is unnecessary and would add
cost.

### 5.5 BS/HM stacking bug

`POLISH_BACKLOG.md` item 0 (BS/HM isosurfaces persist frame-to-frame)
likely shares root cause with the paint-some-frames bug (Section 6).
The two-renderer split does not solve it; the cheap probe
`SetForceTranslucent(true)` on the BS/HM actors does. If the probe
clears the BS/HM stacking, the design's pipeline correctness story is
complete; if it doesn't, the BS/HM stacking is a separate translucency
classification bug to chase in a different session — not blocking the
viewport pipeline land.

---

## 6. The two cheap probes

### 6.1 `molecule_->GetPoints()->Modified()`

Located AFTER the SetAtomPosition loop in `MoleculeScene::setFrame` (the
new slim version per Section 2.4). Existing `mapper_->Modified()` at
`MoleculeScene.cpp:528` is removed.

**The VTK source basis:** `vtkMolecule::SetAtomPosition` calls
`Points->SetPoint(...)` and then `this->Modified()` on the molecule
(`vtkMolecule.cxx:184-197`). `Points->SetPoint` does NOT bump
the points array's own MTime explicitly — only the molecule's MTime
goes up.

The `vtkOpenGLPolyDataMapper`'s VBO gate at
`vtkOpenGLPolyDataMapper.cxx:3742-3767` checks several MTimes against
its cached VBO version; one of those checks (per the VTK source the
observations agent traced) consults the input data's `Points` MTime via
the trivial-producer chain inside `vtkOpenGLMoleculeMapper`. With the
molecule's MTime up but the points' MTime stale, that particular
fragment of the gate can miss when other conditions on its OR chain
fall out.

**The probe:** explicitly bump the points' MTime. One line.

**Expected effect:** end-of-trajectory atom-render drop disappears. If
it does, the cause is the missed points-MTime check. The fix stays in
the codebase as the settled answer (with the comment updated to point
at the VTK source line for posterity).

**If it doesn't:** the probe is removed. The investigation moves to
`notes/RESIDUAL_RENDER_DROP.md` for next-session attention. Per
`feedback_race_condition_deeper_issue`, the probe is documented as a
probe in the comment header until we've confirmed the cause.

### 6.2 `SetForceTranslucent(true)` on BS/HM actors

Located in `QtFieldGridOverlay::Build` and `QtBFieldStreamOverlay::Build`
on the actors created there. The actor property change is one line per
actor.

**The VTK source basis:** `vtkOpenGLRenderer::Render` at
`vtkOpenGLRenderer.cxx:520-690` polls each actor's
`HasTranslucentPolygonalGeometry()` per render to dispatch the
opaque-vs-translucent passes. The BS/HM isosurfaces' per-cell opacity
is computed by the upstream contouring filter; near boundary values
(0.0 / 1.0), classification can return different values frame-to-frame
on the same actor. With unstable classification, the actor flickers
between the opaque pass and the translucent pass; if the order in the
translucent pass varies frame-to-frame (depth-sorting changes), the
visible "stacking" bug per POLISH_BACKLOG.md item 0 appears.

**The probe:** force the classification to translucent with
`SetForceTranslucent(true)`. Always dispatched to the translucent pass
regardless of per-cell opacity values.

**Expected effect:** BS/HM stacking stabilises. Side effect: the
paint-some-frames bug may also clear if its root was BS/HM-translucency
ordering leaking into the atom-imposter pass ordering — the bug
manifests at end-of-trajectory where the BS/HM isosurfaces have
typically equilibrated to certain ring positions, which is consistent
with stable BS/HM helping.

**If it doesn't fix BS/HM stacking:** the actors with that property
return to the default; the BS/HM stacking is a deeper bug to chase
separately. The probe is removed.

### 6.3 Why these two stay in the design doc

Both probes are one-line VTK property changes with a precise hypothesis
backed by VTK source citations. They cost nothing if they don't help
(remove the line and the comment); if they do help, they're documented
fixes with VTK source pointers in the comment for any future engineer
who wonders "why this line". The harness measures the effect
deterministically (paint-some-frames bug + BS/HM stacking are both
visible in the viewport_probe.py marker-area metric — when atoms drop,
the marker_area also drops because the imposter-pass dropout takes the
marker's neighbour atoms with it).

---

## 7. Tests

### 7.1 Pure-function math units (extends `tests/scene_math_tests`)

The existing `tests/scene_math_tests.cpp` covers
`PlaneFrameMath::computePlaneFrame` + `chooseContinuousNormal`. The new
math headers extend the same target:

- **`FitTargetMath::ComputeAtomAnchor`** — trivially `F = positions[0]`,
  no degenerate case beyond input-shape guarantees. Test cases: a
  zero-vector position (returns the zero focal), a non-trivial position
  (returns it).
- **`FitTargetMath::ComputeBondAnchor`** — F = midpoint, axis = (b-a)/|b-a|.
  Test cases: unit-length bond along each cardinal axis, a tilted bond,
  degenerate zero-length bond returns nullopt.
- **`FitTargetMath::ComputeDihedralAnchor`** — Newman-projection
  geometry. Test cases: a canonical trans-180° dihedral on the XY plane
  (axis = +z, viewUp = (1, 0, 0) projected ⊥ z), a gauche dihedral
  (axis = +z, viewUp tilted), zero-length central bond returns nullopt,
  a degenerate (a, b, c, d) where (a - b) is parallel to (c - b)
  returns nullopt for the viewUp branch with a documented fallback.
- **`FitTargetMath::ComputePlaneAnchor`** — wraps
  `math::computePlaneFrame`; the wrapping is mechanical. One test
  case asserts the wrapping forwards to the underlying frame.
- **`FitTargetMath::ComputeSubsetTransform`** — Kabsch. Test cases:
  identity (current == reference returns R = I, T = 0), a known
  rotation (current = reference rotated 90° about z, returns R that
  rotates back), a reflection (det(V*U^T) < 0; returns proper rotation
  via determinant correction), N < 3 returns nullopt.
- **`CameraComposer` per-mode write semantics** (table in §2.3.3). Each
  mode gets one test that constructs a synthetic `CameraComposer` with
  a stub `Conformation`, calls `write(t)` for a known frame, and reads
  back the camera state via the renderer's `GetActiveCamera` setters.
  Tests run against a real `vtkRenderer` (cheap; no rendering required
  — `SetPosition/SetFocalPoint/SetViewUp` are pure state setters).

### 7.2 Existing harness experiments (gating)

`tests/scripts/viewport_probe.py` against the 1P9J fixture:

- **with_plane_lock (median ≤ 5 px AND p95 ≤ 15 px AND max ≤ 30 px)**.
  Today's baseline is median 1.5 / p95 24 / max 507 — the max comes
  from the normal-sign continuity flip handled correctly by
  `chooseContinuousNormal`. The new design preserves that math; the
  max should stay around 500 (one expected flip in 200 frames) but the
  median + p95 should improve marginally (no centroid-delta accumulation
  noise). Conservative gate: **with_plane_lock median ≤ 5 px AND p95 ≤
  30 px AND max ≤ 800 px** — i.e. don't regress today's numbers.
- **without_plane_lock (median ≤ 13 px)** per
  `HARNESS_BASELINE_TRANSFORM_2026-05-30.md`. Today's baseline is
  median 225 px (the centroid-delta bug). The new design removes
  centroid-delta entirely. Free camera with the user's initial
  ResetCamera state should drift only with the protein's actual COM
  motion in display space — at 1P9J's scale, ≤ 5 px / 200 frames.
- **transform_plus_lock (median ≤ 5 px AND p95 ≤ 10 px AND max ≤
  30 px)** per the same baseline doc. The new design's path here is
  `TransformedConformation::Mode::FitSubset(backbone) +
  CameraComposer::AtomMode(14)`. Single Kabsch on positions, atom
  centring on the camera; no Kabsch-vs-Kabsch interference.

### 7.3 New harness experiments (exercise lock-as-fit)

Each per-mode experiment uses the same per-frame screenshot + marker
centroid drift methodology. The harness adds a `--camera-mode` argument
that maps to a setMode REST call.

- **`atom_lock` (mode=Atom, atom=14, frames 200)** — marker should
  stay within 1-2 px of frame 0's centroid (atom 14 has thermal
  vibration only; camera follows it absolutely).
- **`bond_lock` (mode=Bond, atoms=13,14)** — marker should stay within
  the same range; the midpoint between two atoms is steadier than a
  single atom.
- **`dihedral_lock` (mode=Dihedral, atoms=12,13,14,15)** — Newman
  projection; marker drift dominated by torsion of atom 12 / 15 around
  the (13, 14) axis. Hard to set a numerical gate without an atomistic
  baseline; the visual gate is "atom 14 stays at the centre of the
  view, and atoms 12 / 15 visibly orbit it across the trajectory."
- **`subset_lock_backbone` (mode=Subset, atoms=backbone)** — equivalent
  to today's `transform_only`'s Kabsch but applied to the camera.
  Marker drift target: same as `transform_only`'s ~67 px (it's still
  atom-vibration noise within the backbone-rest-frame; the lock can't
  remove that without a transform).
- **`transform_subset_lock_atom` (transform=FitSubset(backbone),
  mode=Atom(14))** — the gate experiment for the
  `transform_plus_lock ≤ 5 / 10 / 30` acceptance. Position-side does
  the backbone rigid-body removal; camera-side does atom centring.
  No Kabsch interference.

### 7.4 Adversarial review

Per `feedback_codex_review_invocation` and
`feedback_adversarial_review_physics` (extended to scene-pipeline math):
after `CameraComposer` + `FitTargetMath` land (before the consolidation
shim wiring lands), spawn a Codex pass on the math headers + the
composer's per-frame write. The math is pure-function and unit-tested,
but sign conventions (axis direction, normal continuity) are exactly the
class of bug the adversarial review catches that smoke tests miss.

---

## 8. Open questions

Honest list. The design conversation that spawned this doc named some
of these; some are this doc's own un-decided edges. Naming them here
is the contract that the implementation session pauses to ask before
making the call.

1. **Lock release semantics when selection changes.** The plane lock
   today releases when `AtomSelection::changed` fires
   (`ReaderMainWindow.cpp:323-329`). The new design's atom / bond /
   dihedral locks: do they release on a selection change too? Argument
   for yes: the user picked the lock target from the selection, so a
   selection change means a new target intent. Argument for no: bond
   lock should survive picking a third atom (the third atom might be
   the start of a dihedral selection). **Proposed default**: release
   for Plane (status quo) and Subset; preserve for Atom, Bond, Dihedral
   (the user can pick more atoms to refine the selection without
   losing the camera lock). User confirms or overrides in the
   implementation session.

2. **User-delta persistence across mode changes.** When the user
   rotates around a locked atom (Atom mode), then switches to Bond
   mode, does the gesture's azimuth/elevation carry over? **Proposed
   default**: per-mode storage; switching modes resets the delta to
   the new mode's natural default. Reasoning: each lock has its own
   semantically meaningful "natural orientation" (Newman for dihedral,
   bond-horizontal for bond, free for atom), and inheriting the
   previous mode's delta would put the camera in an arbitrary state.

3. **Default mode at startup.** Today: ResetCamera is called once
   from `MoleculeScene::ResetCamera`. The new design's natural default
   is `CameraMode::Free` (no lock). **Proposed default**: Free at
   startup; the user opts into a lock via the toolbar or REST.
   Alternative considered + rejected: `CameraMode::Subset(backbone)`
   as the default to "stabilise the protein on screen" — rejected
   because it changes the user's expectation of a static initial view.

4. **Default OrientationPolicy per CameraMode variant.** The table at
   §2.3.3 picks one default per mode. The Plane variant's default
   uses the plane normal as the view axis (sight-perpendicular to the
   plane). For Bond, the default puts the bond horizontal in screen
   space. These are the "obvious" choices but are not the only
   reasonable ones; alternative considered: Plane's default could be
   "look at the plane edge-on" (axis in-plane). **Proposed**: ship the
   table at §2.3.3 as the default; revisit if a user-feedback session
   says otherwise.

5. **Does `SingleConformation` need any of this?** Single-pose
   conformations have one frame; per-frame re-application of a lock
   is a no-op. **Proposed**: the same `CameraComposer` API handles
   both; `Conformation::frameCount() == 1` means the composer runs
   once at setFrame(0) and never again. No code branch on conformation
   kind in the composer.

6. **Touch and 3D-mouse input.** The `QuietTrackballStyle` no-ops the
   standard mouse manipulators; the eventFilter catches Qt-level mouse
   events. Touch events on a touchscreen reach VTK via different paths
   that the quiet-style overrides cover (the style's `OnTouch*` methods
   are also no-opped if VTK 9.5 exposes them; otherwise they fall
   through to the base class no-op of the gesture-specific
   manipulators). **Open**: do we explicitly wire touch through the
   eventFilter? **Proposed**: no for now — the immediate target
   hardware (Linux 5090, GB10 Spark ARM64, Mac M3, Windows 8060S per
   `project_h5reader_target_hardware`) is desktop / laptop with mouse;
   touch can be wired later when a touch user shows up.

7. **What happens to today's `CameraPlaneLock` state struct in
   `MoleculeScene`.** The shim API stays (Section 3.1); the struct
   becomes redundant once the math moves to `CameraComposer`. **Proposed**:
   remove the struct in the same commit that moves the math, even
   though the public API is unchanged. The internal field
   `cameraPlaneLock_` (`MoleculeScene.h:194`) goes too.

8. **Snapshot at scrub start.** Today the slider-drag pile-up at
   `DashboardDisplayController::scrubActive_` defers snapshot fetches
   during drag (`ReaderMainWindow.cpp:365-376`). The viewport refactor
   does not touch this path, but the scrub-active flag could also be
   consulted by the EndEvent observer (Section 8) to throttle UDP log
   lines during scrub — a 30-fps drag pumps 30 log lines per second
   which may saturate the UDP listener. **Proposed**: leave it; the
   UDP listener is local and never the bottleneck in practice. If it
   becomes one, downgrade the per-render log line to DEBUG during
   scrub.

---

## 9. What this design rejects, and why

Explicit anti-shapes. Each is something the design conversation or the
empirical numbers ruled out; documenting them here prevents
re-litigation in a future session.

- **"Lock vs transform" as composing layers** — rejected per the
  harness numbers (`transform_plus_lock` median 39 px / max 795 px vs
  `with_plane_lock` 1.5 px). Two Kabsch fits on overlapping atom
  subsets disagree on which frame is "the reference"; the tail blows
  up. The unification (one fit, written to positions OR the camera OR
  composed via translation-only Atom lock) is the resolution.

- **VTK trackball subclass with `Rotate/Pan/Dolly` override that
  computes "rotate around the locked atom"** — rejected per the Qt+VTK
  boundary lesson. The `QtAtomPicker` precedent at `QtAtomPicker.cpp:48`
  proves Qt eventFilter on the widget is the in-house pattern. The
  trackball subclass would need to know about the active lock target,
  reach into application state, recompute renders — all things the
  eventFilter does naturally, with less VTK-side coupling. The
  `QuietTrackballStyle` (§2.7) is a much smaller subclass (no-op
  overrides only) — it exists to give the adapter's paint chain a
  valid style, not to do camera math.

- **`widget->update()` as the only render verb without renaming
  `Render()` patterns elsewhere** — kept. `widget->update()` IS the
  verb; paint events trigger VTK render via the adapter; the prior
  observations doc framed it slightly differently in one section
  (`VIEWPORT_OBSERVATIONS §5b`'s code sketch) but the right reading is:
  Stage 5's `requestRender` enqueues `widget->update()`; Qt fires
  paintGL; adapter's paint chain calls `iren->Render()`;
  `iren->Render()` calls `renderWindow->Render()` (the gate at
  `vtkRenderWindowInteractor.cxx:170` is OPEN — we keep `EnableRender`
  on); the EndEvent observer fires; Stage 8 logs.

- **Centroid-delta camera translation** — rejected. Three failure modes
  (vibrational noise correlated with eye motion; PBC wraparound; FP
  accumulation across hundreds of frames) per `VIEWPORT_OBSERVATIONS §1.1`,
  validated by the harness (`without_plane_lock` median 225 px). The
  new design absolute-writes the camera every frame from a typed
  CameraMode; the lastCentroid_ / haveLastCentroid_ fields are removed.

- **Adding a second QTimer for "camera smoothing at independent fps"**
  — rejected per `feedback_qt_discipline` and the `qt6-cpp`
  architecture reference's §4. One application timer; the rendering
  fans out from `frameChanged`. Smoothing (interpolation between
  sparse keyframes) is future work; the current bouncing is fit
  precision, not interpolator drift (per the observation that
  `vtkTransform::SetupCamera` internally orthogonalises modelview
  regardless of stored ViewUp — visible projection is fine even with
  non-orthogonal stored ViewUp).

- **A factory or abstract `CameraLock` base class with virtual
  `apply(frame)`** — rejected per `feedback_no_abstractions` and the
  no-pluggable-interfaces project rule. `CameraMode` is a typed sum;
  `CameraComposer::write(t)` dispatches via a switch on the kind. Each
  case calls the appropriate `FitTargetMath::Compute*Anchor` free
  function. No virtual call, no factory, no plugin point.

- **Per-atom visibility on `vtkOpenGLMoleculeMapper` for the distance
  filter use case** — confirmed unavailable in the observations doc
  trace; ghost arrays apply to atoms but not bonds, leaving bonds
  dangling. The clean approach is a separate "display molecule" subset
  (`VIEWPORT_OBSERVATIONS §2.4`). This is **not in this design** — the
  distance filter is a separate feature, not a viewport-pipeline
  concern. It will reuse this design's `CameraMode::Atom` (filter
  centred on the locked atom) when it lands, but the filter itself
  needs its own design pass.

- **Rendering off the GUI thread** — rejected per the qt6-cpp 3d-vtk
  reference §7. All VTK rendering happens on the GUI thread; cross-
  thread render scheduling is a bug (caught by ASSERT_THREAD). Worker
  threads in this codebase (the snapshot loader, the chewer) deliver
  results via signal/slot back to the GUI thread, which then mutates
  VTK state.

- **Removing the every-50-frames diagnostic snapshot block at
  `MoleculeScene.cpp:572-591`** — kept. It catches the bounds-cache
  bug class; cheap; UDP-only; off by default at DEBUG. The new
  per-render UDP log line (Stage 8) replaces the per-frame "render <ms>
  ms" line but the per-50-frames snapshot stays.

---

## Citations summary

For the implementation session — the file:line pointers that ground
each section in either the current h5-reader source or the VTK source.

**h5-reader source (paths under `/shared/2026Thesis/nmr-shielding/h5-reader/`):**
- `src/app/MoleculeScene.cpp:88-91` — stock trackball install (§2.7)
- `src/app/MoleculeScene.cpp:101-111` — EndEvent observer (§2.8)
- `src/app/MoleculeScene.cpp:282-415` — plane lock precedent (§3.1)
- `src/app/MoleculeScene.cpp:418-592` — setFrame to be slimmed (§2.4)
- `src/app/MoleculeScene.cpp:447-459` — per-frame bounds (§2.4)
- `src/app/MoleculeScene.cpp:463-471` — plane-lock branch to remove (§2.4)
- `src/app/MoleculeScene.cpp:476-497` — centroid-delta to remove (§2.4)
- `src/app/MoleculeScene.cpp:501-528` — paint-some-frames probe (§6.1)
- `src/app/MoleculeScene.cpp:530-538` — overlay fan-out (§2.4)
- `src/app/MoleculeScene.cpp:547-553` — explicit-bounds clipping range (§2.4)
- `src/app/MoleculeScene.cpp:624-694` — dihedral one-shot to lift (§3.2)
- `src/app/QtAtomPicker.cpp:48` — eventFilter precedent (§4.1)
- `src/app/QtPlaybackController.{h,cpp}` — the one timer (§2.1)
- `src/app/ReaderMainWindow.cpp:185-189` — TransformedConformation
  construction wiring (§3.4)
- `src/app/ReaderMainWindow.cpp:323-329` — selection-change-releases-lock
  (Open question 1)
- `src/app/ReaderMainWindow.cpp:509-531` — shutdown to reorder (§4.4)
- `src/app/ReaderMainWindow.cpp:722-740` — plane-lock toolbar wiring (§3.1)
- `src/app/MeasurementOverlay.cpp:21-22` — sphere radius/opacity constants (§5.3)
- `src/app/MeasurementOverlay.cpp:34-41` — instrument-mode palette (§5.3)
- `src/app/PlaneFrameMath.h:43-77` — math pattern precedent (§2.3.4)
- `src/model/TransformedConformation.cpp:228-276` — Kabsch fit math (§2.3.4)
- `src/model/AtomSelection.h:56-66` — GeometryKind classification (§1.1)
- `src/model/Conformation.h:73` — atomPosition seam (§2.2)
- `tests/scripts/viewport_probe.py` — the harness (§7)
- `tests/scene_math_tests.cpp` — existing math unit test target (§7.1)

**VTK 9.5.2 source (paths under `/home/jessica/builds/VTK-9.5.2/`):**
- `Common/DataModel/vtkMolecule.cxx:184-197` — SetAtomPosition bumps
  molecule MTime, not points MTime (§2.4, §6.1)
- `Domains/Chemistry/vtkMoleculeMapper.cxx:320-351` — UpdateGlyphPolyData
  gate (§2.4)
- `Domains/Chemistry/vtkOpenGLMoleculeMapper.{h,cxx}` — fast atom mapper
  accessor (§6.1 — used by the residual investigation hypothesis if
  the probe doesn't fix it)
- `GUISupport/Qt/QVTKOpenGLNativeWidget.cxx:90-135` — setRenderWindow
  lifecycle (§4.4)
- `GUISupport/Qt/QVTKOpenGLNativeWidget.cxx:274-305` — paintGL flow (§2.6)
- `GUISupport/Qt/QVTKRenderWindowAdapter.cxx:140-166` — adapter
  destructor calls Finalize (§4.4)
- `GUISupport/Qt/QVTKRenderWindowAdapter.cxx:241-266` — paint() → iren->Render() (§2.5, §2.6, §7)
- `Interaction/Style/vtkInteractorStyleTrackballCamera.cxx:240-273` —
  Rotate (§2.7)
- `Interaction/Style/vtkInteractorStyleTrackballCamera.cxx:264` —
  AutoAdjust ResetCameraClippingRange (§2.7, §4.2)
- `Interaction/Style/vtkInteractorStyleTrackballCamera.cxx:272, 297, 348, 393, 451` —
  five rwi->Render() sites (§2.7)
- `IO/Image/vtkWindowToImageFilter.cxx:57-63` — ShouldRerender (§7 — used
  by the screenshot REST handler)
- `Rendering/Core/vtkCamera.cxx:321-353` — SetViewUp normalises but does
  NOT orthogonalise (§2.3)
- `Rendering/Core/vtkCamera.cxx:360-387` — ComputeViewTransform (§2.3)
- `Rendering/Core/vtkCamera.cxx:561-569` — OrthogonalizeViewUp is
  separate (§2.3)
- `Rendering/Core/vtkRenderWindowInteractor.cxx:168-177` — EnableRender
  gate (§2.7)
- `Rendering/Core/vtkRenderer.cxx:1050-1068` —
  ResetCameraClippingRange forms (§2.4)
- `Rendering/OpenGL2/vtkGenericOpenGLRenderWindow.cxx:73-77` — Frame()
  fires WindowFrameEvent (§2.8)
- `Rendering/OpenGL2/vtkOpenGLPolyDataMapper.cxx:3742-3767` — VBO
  rebuild gate (§6.1)
- `Rendering/OpenGL2/vtkOpenGLRenderer.cxx:520-690` — render-pass
  dispatch (§6.2)
- `Rendering/OpenGL2/vtkOpenGLSphereMapper.cxx:102` — fragment discard
  (mentioned by the observations doc as part of the clipping-range
  hypothesis trail; not in this design's direct path)
