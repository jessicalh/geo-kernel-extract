# Stabilisation Feature Evaluation - 2026-06-04

Scope: read-only code discovery and evaluation. I did not run git or build/test
commands. This note evaluates the real implementation against the stale design
note and against the requested "local interaction view" goal.

## Verdict

The stale note is stale, but not because its design was implemented exactly.
There is now a real stabilisation system, and it is useful, but it is a
two-layer stabilisation mechanism:

1. A display-space rigid transform layer:
   `src/model/TransformedConformation.h:1-18`,
   `src/model/TransformedConformation.cpp:51-66`.
2. A camera-only lock/composition layer:
   `src/app/CameraComposer.h:1-16`,
   `src/app/CameraMode.h:1-30`,
   `src/app/CameraComposer.cpp:277-289`.

This is not yet the user-requested local interaction view. It does not let the
user pick an atom, set a radius, hide everything else, and see a clear
stabilised residue-neighbourhood. It also is not the stale design note's exact
"three selected atoms define a local display frame" mode, because the Kabsch
path explicitly rejects coplanar subsets, and any exact three-atom subset is
coplanar (`src/app/FitTargetMath.h:142-165`,
`src/app/FitTargetMath.h:211-220`,
`tests/scene_math_tests.cpp:741-755`).

Recommendation: build on the good infrastructure, but do not treat the current
plane lock as the local interaction view. Add a new explicit Local Interaction
View state that owns focus atom, radius, visible residue/atom set, transform
choice, and camera choice together.

## What The Stale Note Intended

`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md` says the feature was "not
implemented yet" (`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md:5-7`).
Its intended feature was:

- choose exactly three atoms and build a local reference frame
  (`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md:11-21`);
- stronger than camera plane lock, which was only meant to keep atoms in sight
  (`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md:17-21`);
- implement the stabilised view as a display-space transform so molecule
  geometry, overlays, picking, reveal, measurement, and camera agree
  (`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md:54-68`);
- expose it as "Stabilized view" or "Local frame" only for exactly three atoms
  (`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md:74-86`);
- keep camera plane lock separate from stabilised local view
  (`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md:80-86`).

That note's separation is still conceptually right: camera plane lock and
display-space stabilisation are different promises. The code now has both
kinds of infrastructure, but not the exact local-frame feature described there.

## Real Mechanism

### 1. Display-space transform layer

The real display transform is `TransformedConformation`, a decorator over
`model::Conformation`. Its only decorated virtual is `atomPosition(frame, atom)`:
in non-identity modes it computes or reads a cached rigid transform and returns
`R * raw + T` (`src/model/TransformedConformation.cpp:51-66`). Its header
explicitly says consumers holding a `Conformation*` - scene, picker, overlays,
REST `/positions` - see the wrapped coordinates (`src/model/TransformedConformation.h:13-18`).

Supported modes are:

- `Identity`: raw positions (`src/model/TransformedConformation.h:57-62`);
- `CenterCom`: translate all atoms by equal-weight center of mass
  (`src/model/TransformedConformation.cpp:157-168`);
- `FitReference`: Kabsch fit all atoms to a reference frame
  (`src/model/TransformedConformation.cpp:170-181`);
- `FitSubset`: Kabsch fit a subset, currently used by UI as backbone fit
  (`src/model/TransformedConformation.cpp:183-205`).

`KabschFit()` delegates to the shared `math::ComputeSubsetTransform()` and, on
degenerate input, freezes rotation and only aligns centroids
(`src/model/TransformedConformation.cpp:227-254`). This is a sensible fallback,
but it also means the advertised "FitSubset" is not "always locked"; it can
silently degrade to translation-only when the subset is rank-deficient.

Important boundary: `loadSnapshot()` deliberately forwards the raw snapshot and
does not rewrite the snapshot `Pos` column
(`src/model/TransformedConformation.cpp:68-76`). Most viewport code uses
`atomPosition()` and therefore sees transformed coordinates. Consumers that
read snapshot arrays directly still see source data. That is mostly correct,
but it is a place where "display position" and "raw calculator data" coexist.

### 2. Camera lock/composer layer

The camera side is `CameraComposer`: one owner of absolute per-frame camera
writes, replacing old delta camera paths (`src/app/CameraComposer.h:1-16`).
`CameraMode` defines typed locks: Free, Atom, Bond, Dihedral, Plane, and Subset
(`src/app/CameraMode.h:41-49`). The camera mode documentation calls Plane a
three-atom camera plane lock and Subset a Kabsch camera fit bridge to
`TransformedConformation::FitSubset` (`src/app/CameraMode.h:15-20`).

On every frame, `CameraComposer::write()` dispatches by active mode
(`src/app/CameraComposer.cpp:277-289`). It reads positions through its
`Conformation*`, which is the transformed wrapper when the scene is wired that
way (`src/app/CameraComposer.cpp:270-273`).

Plane mode is camera-only. It:

- reads three displayed atom positions (`src/app/CameraComposer.cpp:524-528`);
- computes a three-atom plane frame;
- applies a normal sign-continuity guard so the camera does not flip sides
  (`src/app/CameraComposer.cpp:532-546`);
- sets focal point to the plane centroid and camera position along the normal
  (`src/app/CameraComposer.cpp:561-564`).

It does not transform molecule coordinates. It does not hide atoms. It is a
view orientation/centering aid.

Subset camera mode is a camera stabiliser, not a display transform. It captures
reference subset positions and camera-relative pose at mode acquisition
(`src/app/CameraComposer.cpp:214-247`), then per frame Kabsch-fits current subset
to reference and rotates the captured camera pose by `R^T`
(`src/app/CameraComposer.cpp:567-652`). If the Kabsch fit is degenerate, it
freezes the camera write for that frame (`src/app/CameraComposer.cpp:593-599`).

User camera input is routed through `CameraInputFilter` into composer gestures
instead of VTK's stock trackball (`src/app/CameraInputFilter.h:1-27`,
`src/app/CameraInputFilter.cpp:43-75`,
`src/app/CameraInputFilter.cpp:95-151`). `QuietTrackballStyle` no-ops VTK camera
manipulators so the composer remains the owner of camera mutation
(`src/app/QuietTrackballStyle.h:1-25`,
`src/app/QuietTrackballStyle.h:39-50`). This is good infrastructure.

### 3. Wiring

`ReaderMainWindow` creates one `TransformedConformation` around the loaded
conformation and connects `transformChanged` to `scene_->refreshCurrentFrame()`
(`src/app/ReaderMainWindow.cpp:184-195`). The scene is then built against the
transformed wrapper (`src/app/ReaderMainWindow.cpp:197-205`).

The scene builds one full `vtkMolecule` actor from `conformation.atomPosition()`
at frame 0 (`src/app/MoleculeScene.cpp:192-221`), then updates every atom's
position from `conformation_->atomPosition()` per frame
(`src/app/MoleculeScene.cpp:336-352`). It creates the composer with the same
conformation pointer (`src/app/MoleculeScene.cpp:225-233`), and fans the same
conformation to ribbon, rings, field grids, B-field streams, selection,
measurement, and reveal overlays (`src/app/MoleculeScene.cpp:235-275`).

This largely satisfies the stale note's "one display-space coordinate system"
guardrail:

- picking raycasts against `conformation_->atomPosition()`
  (`src/app/QtAtomPicker.cpp:91-125`);
- measurement markers and polylines use `atomPosition()`
  (`src/app/MeasurementOverlay.cpp:223-289`);
- reveal markers/lines use `atomPosition()`
  (`src/app/SceneRevealOverlay.cpp:448-491`);
- the atom inspector shows display `xyz` via `atomPosition()`, while detail data
  still comes from snapshots (`src/app/QtAtomInspectorDock.cpp:287-296`);
- REST `/positions` explicitly returns active display-frame coordinates
  (`src/app/RestServer.cpp:970-1009`).

### 4. User-visible controls today

Toolbar camera actions expose Focus, Newman, Plane lock, and Free
(`src/app/ReaderMainWindow.cpp:737-783`). Plane lock requires exactly three
selected atoms and is implemented by the scene's compatibility shim
(`src/app/MoleculeScene.cpp:468-530`).

The transform menu exposes only:

- Identity;
- Center COM;
- Fit reference;
- Fit backbone.

See `src/app/ReaderMainWindow.cpp:787-825`. `FitSubset` in the UI always means
"backbone subset against reference frame 0"
(`src/app/ReaderMainWindow.cpp:1060-1074`). There is no UI action for "fit
selection", no local-frame toggle, no radius, and no isolation.

REST exposes more than the UI:

- `/transform` supports identity, center_com, fit_reference, and fit_subset,
  including explicit `subset_atoms` or `backbone_only`
  (`src/app/RestServer.cpp:454-564`);
- `/plane-lock` and `/plane-lock/enable` are old plane-lock shims
  (`src/app/RestServer.cpp:566-610`);
- `/camera/mode` supports free, atom, bond, dihedral, plane, and subset
  (`src/app/RestServer.cpp:612-825`).

So the feature is scriptable, but not exposed as a cohesive local interaction
view.

## State

### Works

- Global/backbone display stabilisation exists and is wired through the main
  viewport path.
- Camera plane lock exists and is tested end-to-end over REST. The test asserts
  active state, focal point tracking the plane centroid, camera direction
  parallel to the plane normal, direction continuity, and nonblank screenshots
  (`tests/rest/test_camera_plane_lock.py:1-16`,
  `tests/rest/test_camera_plane_lock.py:85-123`,
  `tests/rest/test_camera_plane_lock.py:126-164`).
- The historical harness says plane lock and transform both contributed:
  camera-only plane lock had median drift 1.52 px, raw/no-lock had 225.32 px,
  transform-only backbone fit had 67.26 px, and transform+lock had 39.01 px
  (`tests/scripts/HARNESS_BASELINE_TRANSFORM_2026-05-30.md:107-186`). Treat this
  as historical evidence, not a fresh run for this note.

### Partial

- Transform and camera locks are separate switches. There is no single feature
  that decides "for this focus atom and radius, use this atom set for display
  isolation, this transform subset, and this camera anchor".
- The code has a camera `SubsetMode`, but the normal UI does not expose
  arbitrary subset camera stabilisation.
- The code has display `FitSubset`, but the normal UI exposes only backbone
  subset. REST can send explicit subset atoms.
- A three-atom plane lock exists, but it is camera-only.
- A three-atom display Kabsch local frame, as suggested by the old note, is not
  actually possible with the current shared Kabsch policy because coplanar
  subsets return `nullopt`.

### Missing for the user's goal

The requested goal is: pick an atom, set a radius, see only that atom plus the
residues within that radius, and have a stabilised/zoomed view that makes the
interacting residues clear.

Missing pieces:

- no viewport visibility mask for atoms/residues;
- no separate local `vtkMolecule` actor;
- no "whole residues within radius" scene mode;
- no camera zoom/framing policy tied to the visible local set;
- no unified UI/REST command for focus atom + radius + isolation +
  stabilisation;
- no transform subset policy based on the radius-local atom/residue set.

`MoleculeScene` currently owns one full `vtkMolecule` actor with all atoms and
bonds (`src/app/MoleculeScene.cpp:192-221`) and only moves all atoms per frame
(`src/app/MoleculeScene.cpp:336-352`). There is no per-atom visibility or local
actor path in that code.

`AtomSelection` is intentionally capped at four atoms for distance/angle/dihedral
measurements (`src/model/AtomSelection.h:63-66`), so it cannot be the model for
"all residues in a radius".

## Existing Radius Logic

There is useful nearby-radius logic, but it is not a viewport isolation feature.
`NearbySignalModel` is explicitly "discovery state, not atom selection"
(`src/app/NearbySignalModel.h:1-5`). It stores a radius and focus anchor
(`src/app/NearbySignalModel.h:72-80`, `src/app/NearbySignalModel.h:90-95`).
The radius is clamped to 1.0-30.0 Angstrom
(`src/app/NearbySignalModel.cpp:224-230`). `rebuild()` computes focus position,
atom distances, residue best distances, bonds, rings, ring memberships, and
bond-vector candidates using `conformation_->atomPosition()`
(`src/app/NearbySignalModel.cpp:322-505`).

`SignalDisplayDialog` exposes that as a "Selection context" candidate table with
Live and Radius controls (`src/app/SignalDisplayDialog.cpp:670-701`,
`src/app/SignalDisplayDialog.cpp:811-813`,
`src/app/SignalDisplayDialog.cpp:1026-1029`).

This is salvageable for membership computation: it already answers "what atoms,
residues, bonds, rings, and signal anchors are near this focus atom in display
coordinates?" But it is currently a dashboard binding picker, not a scene
visibility controller.

## Cruft And Confusion

### Cut or retire: dormant `QtSelectionOverlay`

`QtSelectionOverlay` says it listens for atom picks and frame changes
(`src/app/QtSelectionOverlay.h:1-8`), but it is not connected to the current
selection flow. The current multi-atom `MeasurementOverlay` is documented as its
successor (`src/app/MeasurementOverlay.h:1-15`). `MoleculeScene.h` labels it
"dormant; superseded by measurement_" (`src/app/MoleculeScene.h:225-230`).

Despite that, the scene still constructs and builds it
(`src/app/MoleculeScene.cpp:257-260`) and calls `selection_->setFrame(t)` every
frame (`src/app/MoleculeScene.cpp:408-415`). This is real cruft. It should be
deleted or at least not constructed/called in normal scene setup.

### Design mismatch: exact three-atom Kabsch local frame

The stale design note names "Kabsch fit over the three selected atoms" as an
acceptable first rule (`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md:43-49`).
The current shared Kabsch helper rejects rank-deficient inputs including
coplanar subsets (`src/app/FitTargetMath.h:142-165`,
`src/app/FitTargetMath.h:211-220`). The test suite asserts coplanar input
returns nullopt (`tests/scene_math_tests.cpp:741-755`).

That policy is defensible for full 3D Kabsch camera stability, but it means the
old "exactly three atoms" stabilised display Kabsch idea cannot work as written.
If a future local-frame feature wants exactly three atoms, it needs a separate
plane-frame/wedge stabiliser, not the current rank-3 Kabsch helper.

### Duplicate or misleading plane-frame continuity

`PlaneFrameMath` still exposes `chooseContinuousNormal()`
(`src/app/PlaneFrameMath.h:14-24`, `src/app/PlaneFrameMath.h:74-77`), but
`CameraComposer::writePlane()` implements its own normal continuity state
directly (`src/app/CameraComposer.cpp:532-546`). The helper is still covered by
tests (`tests/scene_math_tests.cpp:242-275`), so this is not a blind delete, but
production ownership is confusing: plane-frame math and camera composer each
carry part of the same story.

### Misleading comments and names

- `CameraAnchorHelper.h` claims a residue N/CA/C plane lock makes the immediate
  neighborhood coherent and that the focus atom appears at the plane centroid
  (`src/app/CameraAnchorHelper.h:1-8`). The implementation's own comment is more
  accurate: the focus atom sits within a few Angstrom of that centroid
  (`src/app/CameraAnchorHelper.cpp:51-56`). This helper derives a camera plane
  or dihedral anchor; it does not implement a local interaction view.
- `ReaderMainWindow.h` says the transform menu is "Identity / Center COM / Fit
  backbone / Fit selection" (`src/app/ReaderMainWindow.h:231-232`), but the UI
  exposes Identity, Center COM, Fit reference, Fit backbone
  (`src/app/ReaderMainWindow.cpp:795-805`). There is no Fit selection action.
- `TransformedConformation.h` says `KabschFit` returns identity for `N < 3`
  (`src/model/TransformedConformation.h:125-130`), but implementation delegates
  to `ComputeSubsetTransform()` and then applies a centroid-translation fallback
  for any degenerate/nullopt case (`src/model/TransformedConformation.cpp:231-249`).

### REST validation gap

`POST /transform` parses `reference_frame` and passes it to
`TransformedConformation::setMode()` without checking it against
`frameCount()` (`src/app/RestServer.cpp:503-505`,
`src/app/RestServer.cpp:559`). `TrajectoryConformation::atomPosition()` reads
`positions()->at(atomIdx, frame)` without an explicit frame guard
(`src/model/TrajectoryConformation.cpp:36-39`). UI uses reference frame 0, so
this is mostly a REST robustness issue today.

`TransformedConformation::setMode()` also skips out-of-range subset atoms when
building `referencePositions_`, but keeps the original `subsetAtoms_`
(`src/model/TransformedConformation.cpp:103-112`). `computeTransform()` then
clips by the smaller size and slices the reference vector by `current.size()`
(`src/model/TransformedConformation.cpp:188-204`). REST validates subset atoms,
so this is not a likely user-visible bug in current normal use, but the helper's
internal contract is brittle.

### AI-cruft mixed with gifts

Several comments contain session-history language such as "Codex finding",
"agent decision", "implementation prompt", and memory-entry names. Examples:

- `src/app/CameraComposer.cpp:68-77`;
- `src/app/CameraComposer.cpp:589-607`;
- `src/model/TransformedConformation.cpp:210-226`;
- `src/app/FitTargetMath.h:142-164`;
- `src/app/ReaderMainWindow.cpp:184-190`.

Some of these comments preserve useful bug provenance. The cruft is that product
code now reads like an implementation transcript. Keep the technical facts
(why freezing on rank-deficient Kabsch exists, why VTK trackball is quiet), but
collapse the session provenance into concise comments or notes.

### Keep, do not cut as cruft: compatibility shims

The `/plane-lock` REST routes and `MoleculeScene::lockCameraToSelectionPlane()`
are compatibility shims, but they are still useful and tested
(`src/app/RestServer.cpp:566-610`,
`src/app/MoleculeScene.cpp:468-530`,
`tests/rest/test_camera_plane_lock.py:71-123`). Keep them while adding the new
local-view feature.

## Salvageable Parts

Keep and build on:

- `TransformedConformation`: it is the right display-space seam. Most viewport
  consumers already agree on its coordinates.
- `CameraComposer`, `CameraInputFilter`, and `QuietTrackballStyle`: they provide
  the correct "absolute camera write" ownership model.
- `FitTargetMath` and `PlaneFrameMath`: keep the pure math, but separate
  rank-3 Kabsch from any future rank-2/three-atom plane stabiliser.
- `NearbySignalModel` logic: extract or reuse the focus-radius distance and
  residue membership computation, but do not overload it as scene state.
- `MeasurementOverlay` and `SceneRevealOverlay`: use them for focus/local
  highlights instead of resurrecting `QtSelectionOverlay`.
- REST routes: useful for demo scripting, especially `/transform`,
  `/camera/mode`, `/camera/focus_atom`, `/positions`, and `/screenshot`.

## Gap To The Goal

The existing stabilisation gets partway to "clear interactions" only in the
camera/coordinate sense. It can reduce whole-protein rigid drift and keep a
target plane or atom family in view. It does not solve visual clutter, which is
the main thing the requested radius view needs.

For the user's goal, the primary missing feature is not more Kabsch math. It is
local visibility:

- focus atom;
- radius;
- all residues with at least one atom within radius;
- only those residues' atoms and bonds visible;
- camera zoom/framing around that visible set;
- optional stabilisation so the local set remains steady enough to compare
  frames.

The current full-molecule actor path makes everything visible all the time.
Without isolation, even a stabilised camera can still leave the user guessing
which residues matter.

## Cleanest Path For A Desk-ready Demo Tomorrow

Build on the current infrastructure, but create a new explicit local-view mode.

Suggested minimal path:

1. Add a Local Interaction View state outside `AtomSelection`: focus atom,
   radius, visible residues, visible atoms, reference frame, and enabled flag.
   Do not use `AtomSelection`; it is capped at four atoms by design.

2. Reuse/extract `NearbySignalModel` distance logic to compute membership:
   find all atoms within radius of the focus atom, lift to whole residues, then
   include every atom and bond internal to those residues. For "clear
   interactions", showing whole residues is more useful than showing only
   radius-hit atoms.

3. Add a local molecule rendering path. For a fast demo, prefer a separate
   `vtkMolecule` actor containing only visible atoms/bonds and hide the global
   molecule actor while local view is active. That is simpler and safer than
   trying to force per-atom visibility into `vtkOpenGLMoleculeMapper`.

4. Frame the local set explicitly. Use the visible atom bounds/centroid to set
   camera distance and clipping. Existing `CameraComposer` can own the write,
   but the feature should set one coherent camera mode rather than asking the
   user to combine toolbar actions manually.

5. Stabilise conservatively for the demo:
   - Use `TransformedConformation::FitSubset` on a non-coplanar subset when
     available: local visible backbone atoms plus nearby atoms, or fallback to
     whole backbone.
   - Use the current frame as the reference when entering local view, not always
     frame 0.
   - Do not attempt exact three-atom display Kabsch tomorrow; the current math
     rejects it by design.
   - Plane lock can be used as an orientation aid, but not as the feature.

6. Filter or suppress nonlocal overlays while local view is active. Measurement
   and reveal overlays can highlight focus/interacting atoms. Full-protein
   ribbon/ring/field overlays will likely reintroduce clutter unless filtered
   to visible residues or hidden.

7. After the demo, refactor:
   - retire `QtSelectionOverlay`;
   - clean stale comments/names;
   - formalise the local-view controller;
   - decide whether there is a separate rank-2/three-atom plane-frame
     stabiliser, distinct from rank-3 Kabsch.

## Final Recommendation

Build on the existing stabilisation, but do not just polish the old plane lock.
The good core is real: display-space `TransformedConformation` plus
composer-owned camera locks. The cruft is mostly around stale UI promises,
dormant overlay code, session-transcript comments, and the unresolved mismatch
between "three-atom local frame" and rank-3 Kabsch.

For the requested local interaction view, the clean path is:

1. implement radius-based residue isolation as a new scene mode;
2. render the isolated set in a separate local molecule actor;
3. use existing transform/camera infrastructure to stabilise and zoom that
   local actor;
4. keep plane lock as an optional orientation aid, not the definition of the
   feature.
