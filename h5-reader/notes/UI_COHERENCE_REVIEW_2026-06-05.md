# UI Coherence Review - 2026-06-05

Scope: read-only code-level review of the current `h5-reader-pysr-spike`
working tree, focused on the lead's five UI/3D coherence concerns. I did not
build, run, switch branches, use git, or touch the running reader. This report
is the only requested filesystem output.

## Top-line verdict

The app has the right architectural direction for molecular display coherence:
display transforms are centralized at the `Conformation::atomPosition()` seam,
the scene mostly renders through that seam, and the inspector/catalog code has
started to distinguish absent data from zero-valued data. But two user-visible
gaps remain high priority:

1. Camera clipping is still stale on free-camera zoom and drag. The scene
   updates clipping when the frame changes or plane lock is enabled, but not
   after ordinary camera gestures. This explains the reported zoom-in cutoff.
2. Field availability is not a first-class loaded-data concept. The inspector
   can show dash-only groups, and the metric picker/catalog can offer channels
   that are implemented in code but absent, all-gap, or structurally empty in
   the loaded run.

The pane toggles are not literally missing in code. They are added via
`QDockWidget::toggleViewAction()`, but they are weakly discoverable because they
live as text actions at the far end of a dense playback toolbar and have no
normal View menu or panel control.

## 1. Display Transforms

### What the code does

`TransformedConformation` is a decorator over `Conformation` that applies a
per-frame rigid transform only at `atomPosition(frame, atom)`.

- The header explicitly says all virtuals delegate to the inner conformation
  except `atomPosition`, and that `MoleculeScene`, overlays, picker, and REST
  `/positions` see the wrapper through polymorphism
  (`src/model/TransformedConformation.h:1`, `src/model/TransformedConformation.h:13`).
- `atomPosition()` returns raw inner positions in Identity mode; otherwise it
  computes/caches a frame transform and returns `R * raw + T`
  (`src/model/TransformedConformation.cpp:51`).
- Snapshots are deliberately not decorated: `loadSnapshot()` forwards the
  inner snapshot unchanged so snapshot `Pos` and calculator NPYs stay in the
  source frame (`src/model/TransformedConformation.cpp:68`).
- `ReaderMainWindow` wraps the loaded conformation once, connects
  `transformChanged` to `scene_->refreshCurrentFrame()`, and builds the scene
  against the wrapper (`src/app/ReaderMainWindow.cpp:184`,
  `src/app/ReaderMainWindow.cpp:191`, `src/app/ReaderMainWindow.cpp:203`).
- `/positions` also explicitly reads through the wrapper and documents that it
  returns the active display frame, not raw H5 positions
  (`src/app/RestServer.cpp:970`, `src/app/RestServer.cpp:980`,
  `src/app/RestServer.cpp:999`).

The available modes are Identity, CenterCom, FitReference, and FitSubset
(`src/model/TransformedConformation.h:57`). The toolbar exposes them as
Identity, Center COM, Fit reference, and Fit backbone
(`src/app/ReaderMainWindow.cpp:803`). The UI FitSubset action hard-codes
reference frame 0 and the typed backbone subset
(`src/app/ReaderMainWindow.cpp:1069`).

### Coherence verdict

The transform layer is mostly coherent and should not be ripped out. It is the
right place to remove rigid-body drift for display because every renderer-like
consumer already asks for `Conformation::atomPosition()`.

The source/display contract is also mostly coherent, but it must be stated
plainly in the UI and docs:

- H5 trajectory positions, snapshots, and calculator payloads remain in the
  source frame.
- Scene positions, overlays that use `Conformation::atomPosition()`, picker
  math, and REST `/positions` are in the active display frame.

That split is defensible, but it can surprise a user who expects `/positions`
or inspector position values to be raw H5 data. The code comments already say
this; the UI should expose "display frame" vs "source frame" where it matters.

### Issues and redundancy

- "Center COM" is not center of mass. The implementation uses an equal-weight
  centroid over all atoms (`src/model/TransformedConformation.cpp:157`). Rename
  it to "Center centroid" or implement mass weighting.
- "Fit reference" and "Fit backbone" are not equally useful UI choices. The
  all-atom fit removes global translation/rotation using every atom, which can
  smear sidechain/internal motion into the fit. Backbone fit is usually the
  better presentation default for proteins because it removes rigid-body drift
  while leaving sidechain motion readable. Keep all-atom fit only as an
  advanced/debug option unless there is a real advisor workflow for it.
- The toolbar transform menu does not let the user choose the reference frame.
  It always uses frame 0 for every UI fit path (`src/app/ReaderMainWindow.cpp:1077`).
  REST supports `reference_frame` (`src/app/RestServer.cpp:503`), so this is a
  UI under-specification rather than a model limitation.
- FitSubset handles invalid/out-of-range subset atoms in a brittle way inside
  the model: reference caching skips invalid atoms, while current-frame
  sampling iterates `subsetAtoms_` and then clips the reference to
  `current.size()` (`src/model/TransformedConformation.cpp:103`,
  `src/model/TransformedConformation.cpp:183`). REST validates subset indices
  (`src/app/RestServer.cpp:540`), and the toolbar uses a typed backbone subset,
  so this is low practical risk, but the model itself is not a clean public
  boundary.
- The Kabsch path is shared and documented, which is good. Degenerate fits
  return `nullopt` from `ComputeSubsetTransform()` on underdetermined or
  rank-deficient input (`src/app/FitTargetMath.h:130`,
  `src/app/FitTargetMath.h:185`). `TransformedConformation::KabschFit()` then
  freezes rotation and still aligns centroids when it has matched pairs
  (`src/model/TransformedConformation.cpp:210`,
  `src/model/TransformedConformation.cpp:238`). That is a sensible display
  policy, but the header comment saying FitSubset is simply a Kabsch fit hides
  this fallback behavior.
- No PBC unwrap is implemented in this layer by design
  (`src/model/TransformedConformation.h:25`). The fit modes can stabilize an
  already-whole trajectory, but they cannot make broken molecules whole.

### Recommended action

Keep the transform decorator. Make Fit backbone the prominent stabilizing mode,
rename or fix Center COM, add a visible reference-frame choice if fit modes are
meant for users, and label positions as "display frame" when the active
transform is not Identity. Consider making all-atom Fit reference advanced or
REST-only unless the lead wants it in the primary UI.

## 2. Camera Clipping Cutoff On Zoom

### What the code does

The scene disables VTK's automatic interactor clipping adjustment:

- `style->AutoAdjustCameraClippingRangeOff()` is called in the scene setup
  (`src/app/MoleculeScene.cpp:129`).

The scene does compute real atom bounds on frame changes:

- `PushAtomPositions()` writes each transformed atom position into the VTK
  molecule and accumulates bounds from the same transformed coordinates
  (`src/app/MoleculeScene.cpp:336`).
- `setFrame()` then pads those bounds and calls
  `renderer_->ResetCameraClippingRange(padded)`
  (`src/app/MoleculeScene.cpp:417`, `src/app/MoleculeScene.cpp:428`).

Plane-lock enable has an extra clipping reset:

- After `composer_->write(frame)`, `lockCameraToSelectionPlane()` recomputes
  current bounds and resets clipping because the camera moved and stale
  clipping may cut into the molecule (`src/app/MoleculeScene.cpp:481`,
  `src/app/MoleculeScene.cpp:487`, `src/app/MoleculeScene.cpp:499`).

Free-camera gestures do not do this:

- Mouse drag and wheel call `composer_->applyGesture()` and then only
  `scene_->requestRender(RenderSource::CameraInput)`
  (`src/app/CameraInputFilter.cpp:95`, `src/app/CameraInputFilter.cpp:141`).
- `requestRender()` schedules a VTK render but does not reset clipping
  (`src/app/MoleculeScene.cpp:299`).
- In Free mode, `CameraComposer::applyGesture()` directly writes the live
  camera position/focal point/view-up, including dolly zoom
  (`src/app/CameraComposer.cpp:785`, `src/app/CameraComposer.cpp:840`,
  `src/app/CameraComposer.cpp:844`).

### Root cause

The zoom cutoff is very likely stale near/far clipping after free-camera dolly.

The app has intentionally disabled the interactor's automatic clipping-range
maintenance, but it only resets clipping when the frame changes or when plane
lock is enabled. A free-mode wheel zoom moves the camera closer immediately and
then renders using the previous near/far planes. Once the camera crosses close
enough to the molecule, the old near plane can slice through the atom actors and
box-like overlay geometry. The code in `lockCameraToSelectionPlane()` already
documents this exact failure mode for camera moves.

### Recommended fix

Add a single scene-owned clipping sync path and call it after every camera write.

Suggested shape:

1. Store the last padded molecule bounds in `MoleculeScene` during `setFrame()`.
   This avoids calling `PushAtomPositions()` again on every wheel tick.
2. Add `MoleculeScene::syncCameraClippingRange()` that calls
   `renderer_->ResetCameraClippingRange(lastPaddedBounds_)` when the scene has
   valid bounds.
3. After `CameraComposer::applyGesture()` in `CameraInputFilter`, call
   `scene_->syncCameraClippingRange()` before `requestRender()`.
4. Also call the same helper after `ResetCamera()`, REST camera set/clear, and
   any future camera-write site. The current plane-lock special case can then
   collapse into the same helper.

This is a timing/staleness issue, not a fundamental problem with using padded
atom bounds. If the overlays can extend farther than the 5 Angstrom pad, make
the cached bounds overlay-aware later, but the immediate zoom cutoff does not
need a new camera system.

## 3. Empty Fields And All-Zero Channels

### What the inspector does

The inspector uses optional-aware adders, but they render missing values as
dash rows (`src/app/QtAtomInspectorDock.cpp:111`). Many top-level groups are
created unconditionally after a snapshot is present:

- Ring current (`src/app/QtAtomInspectorDock.cpp:313`)
- Quadrupole / Dispersion (`src/app/QtAtomInspectorDock.cpp:326`)
- Bond anisotropy (`src/app/QtAtomInspectorDock.cpp:333`)
- Electrostatics (`src/app/QtAtomInspectorDock.cpp:350`)
- H-bond (`src/app/QtAtomInspectorDock.cpp:362`)
- SASA (`src/app/QtAtomInspectorDock.cpp:374`)
- Water (`src/app/QtAtomInspectorDock.cpp:382`)
- Charges (`src/app/QtAtomInspectorDock.cpp:404`)
- Planar geometry (`src/app/QtAtomInspectorDock.cpp:425`)

Some groups are already conditionally hidden when no meaningful data exists,
for example DSSP, energy, MOPAC, ORCA, tripeptide, and Larsen groups
(`src/app/QtAtomInspectorDock.cpp:415`, `src/app/QtAtomInspectorDock.cpp:439`,
`src/app/QtAtomInspectorDock.cpp:457`, `src/app/QtAtomInspectorDock.cpp:473`).

So the inspector is half-correct: absence is represented, but the UI can still
show categories full of "missing" rows. The lead's "show the user nothing if
nothing exists" concern is valid.

### What the dashboard/catalog does

The metric catalog mostly knows whether a sampler is implemented in code, not
whether data exists in the loaded run.

- `hasImplementedTemporalSampler()` returns true from a static set of storage
  paths for dense H5 and returns true for all `FrameNpySnapshot` descriptors
  (`src/model/TrajectorySignalCatalog.cpp:206`,
  `src/model/TrajectorySignalCatalog.cpp:278`).
- `makeDescriptor()` turns that into `SampleStatus::Valid`
  (`src/model/TrajectorySignalCatalog.cpp:341`).
- `SignalDisplayDialog::refreshCatalog()` loads the full descriptor list into
  the candidate table without applying loaded-run availability
  (`src/app/SignalDisplayDialog.cpp:950`).
- `TrajectorySignalCatalog::filterDescriptors()` and `canSample()` only check
  descriptor metadata and `samplingStatus`, not actual loaded data
  (`src/model/TrajectorySignalCatalog.cpp:1188`,
  `src/model/TrajectorySignalCatalog.cpp:1264`).

The sampling path can represent gaps after a signal is chosen:

- `DashboardDisplayController::extendToFrame()` appends the sampler result or
  a gap (`src/app/DashboardDisplayController.cpp:2473`).
- `SignalBuffer::append()` converts gap samples into null optional values
  (`src/model/SignalTimeSeries.cpp:33`).
- `StripStackWidget` ignores invalid/non-finite samples while computing ranges
  and painting (`src/app/StripStackWidget.cpp:40`,
  `src/app/StripStackWidget.cpp:122`).

That makes the strip chart tolerant of absent data, but tolerance after
selection is not enough. The picker can still offer options that will produce
empty strips.

There is also old API surface where absence is collapsed to zeros:

- `QtFrame` documents an "absent-is-zero" backward-compat policy
  (`src/model/QtFrame.cpp:1`).
- Several accessors return zero or placeholder values for data not in the new
  format, such as total B field, McConnell scalar helpers, water fields, and
  EEQ charge (`src/model/QtFrame.cpp:103`, `src/model/QtFrame.cpp:201`,
  `src/model/QtFrame.cpp:311`, `src/model/QtFrame.cpp:336`).

The current inspector mostly avoids `QtFrame`, but any future UI fed from this
legacy surface will reintroduce fake zero values unless it also carries
availability metadata.

### Recommended availability model

Add a loaded-run availability layer and use it in both the inspector and the
signal catalog.

The minimum useful model:

- `Absent`: source object/dataset is not present.
- `NoFramePayload`: frame-local NPY snapshot is not sampled for this frame.
- `AllMissing`: present flags or source-attached masks say no valid rows.
- `AllZeroStructural`: mathematically fixed/structural zero, not a useful user
  channel.
- `AllZeroObserved`: data exists but all observed values are zero. This should
  be shown only when zero is scientifically meaningful for that field.
- `Available`: at least one meaningful finite sample exists.

Implementation points:

- Build `TrajectoryFieldAvailability` after `QtTrajectoryH5` load and after
  frame-NPY sample discovery. It should live near `model`/`io`, not in widget
  code.
- For dense H5 descriptors, resolve `storagePath` to the typed time-series
  pointer, check dimensions, check source-attached/present masks when available,
  and scan finite values. Use descriptor semantics to classify structural zeros.
- For frame NPY descriptors, scan the sampled snapshot group views or a cached
  snapshot schema summary. If loading every sampled frame is expensive, build a
  lazy cache on first opening the metric dialog, but still present it as
  loaded-run availability, not as per-widget guesses.
- Change inspector adders to return `bool` and only create a group if at least
  one child row was added. Dash-only groups should disappear.
- Change `TrajectorySignalCatalog` or `SignalDisplayDialog` to expose only
  descriptors/channels whose availability is not `Absent`, `AllMissing`, or
  `AllZeroStructural`. `canSample()` should also consult availability.
- Remove the default startup signal
  `npy:dssp_chi` unless the descriptor is available in the loaded run
  (`src/app/ReaderMainWindow.cpp:306`).

This should be treated as one shared availability layer, not ad-hoc filters in
each UI panel. Otherwise the inspector, metric dialog, and strip chart will keep
disagreeing.

## 4. Pane Toggle Buttons

### What the code does

The docks are created and then hidden at launch:

- Inspector dock is created earlier; selection and dashboard strip docks are
  created and tabified with it (`src/app/ReaderMainWindow.cpp:373`,
  `src/app/ReaderMainWindow.cpp:381`).
- All three are set invisible on launch
  (`src/app/ReaderMainWindow.cpp:412`, `src/app/ReaderMainWindow.cpp:417`).

The toggles do exist:

- After dock creation, `ReaderMainWindow` appends each dock's
  `toggleViewAction()` to the playback toolbar and relabels them "Inspector",
  "Selection", and "Strip" (`src/app/ReaderMainWindow.cpp:421`,
  `src/app/ReaderMainWindow.cpp:438`).

There is also a REST harness route to hide/restore all docks:

- `/docks/visible` calls `ReaderMainWindow::setDocksVisible()`
  (`src/app/RestServer.cpp:427`).
- The restore path returns each dock to its stashed prior visibility, so a dock
  that was already hidden stays hidden (`src/app/ReaderMainWindow.cpp:602`,
  `src/app/ReaderMainWindow.cpp:626`).

### Coherence verdict

The lead's concern is directionally right but needs precision: the pane toggles
are present in code, but they are not strong UI affordances.

Problems:

- They are appended to a dense playback toolbar after play controls, a 400px
  slider, camera actions, transform, Instrument, Metrics, and overlay toggles
  (`src/app/ReaderMainWindow.cpp:707`, `src/app/ReaderMainWindow.cpp:734`,
  `src/app/ReaderMainWindow.cpp:756`, `src/app/ReaderMainWindow.cpp:848`,
  `src/app/ReaderMainWindow.cpp:856`). On narrower windows they can be pushed
  into toolbar overflow.
- There is no View menu fallback using the same actions.
- The text actions are not visually grouped as "Panels", so a user can miss
  how to reopen a hidden inspector after picking an atom.
- The docks start hidden and picking an atom intentionally does not force-show
  the inspector (`src/app/ReaderMainWindow.cpp:412`), so discoverability matters
  more than usual.

### Recommended action

Keep `QDockWidget::toggleViewAction()` because it gives correct two-way Qt
binding. Move the actions into a proper `View -> Panels` menu, add a compact
toolbar panel control or visible panel group, and consider icons or a single
"Panels" menu button if toolbar width is tight. The same `QAction`s should be
shared between menu and toolbar so checked state stays correct.

## 5. Accreted Ad-hoc Cruft

This section separates code to keep, code to settle, and code to remove or
fold into a cleaner abstraction.

### Keep, but clean up

- `TransformedConformation`: keep. It is the correct architectural seam.
  Clean up labels and UI semantics as described above.
- `QDockWidget::toggleViewAction()` panel binding: keep. Add a real View menu
  and clearer panel affordance.
- Shared Kabsch implementation in `FitTargetMath.h`: keep. It has one
  degeneracy policy used by both transform and camera code
  (`src/app/FitTargetMath.h:130`, `src/model/TransformedConformation.cpp:210`).

### Must fix

- Free-camera clipping reset is missing. Add a scene-owned clipping sync helper
  and call it after camera gestures and camera REST changes
  (`src/app/CameraInputFilter.cpp:141`, `src/app/CameraComposer.cpp:785`,
  `src/app/MoleculeScene.cpp:417`).
- Field availability is not centralized. Add loaded-run availability and make
  the inspector/catalog consume it.
- The default dashboard seed signal is hard-coded to `npy:dssp_chi`
  (`src/app/ReaderMainWindow.cpp:306`). Gate it on availability or remove it.

### Cut or fold in

- `QtSelectionOverlay` appears superseded. It is still built by `MoleculeScene`
  (`src/app/MoleculeScene.cpp:257`), while `MeasurementOverlay` is wired to the
  actual `AtomSelection` model (`src/app/ReaderMainWindow.cpp:350`). The header
  comment says `MoleculeScene` keeps it dormant/superseded
  (`src/app/MoleculeScene.h:229`), and `MeasurementOverlay.h` calls itself the
  deliberate successor. Remove `QtSelectionOverlay` or fold any remaining
  useful marker behavior into `MeasurementOverlay`.
- The Linux-only `/proc/self/statm` diagnostic lives in the render-frame path
  (`src/app/MoleculeScene.cpp:446`). Move it to a diagnostics helper or remove
  it after the render-drop investigation is closed.
- The `PROBE` comment around points `Modified()` is no longer a settled-code
  comment (`src/app/MoleculeScene.cpp:372`). Either document it as the permanent
  VTK invalidation fix or remove it if the investigation ended elsewhere.
- `ReaderMainWindow`'s constructor/toolbar wiring uses zero-delay
  `QTimer::singleShot()` to connect actions after scene/playback exist
  (`src/app/ReaderMainWindow.cpp:885`). This is workable but fragile. Prefer
  constructing the toolbar after dependencies exist, or split action creation
  from dependency wiring in an explicit `wireToolbarActions()` call.
- The toolbar has an "Instrument" user-facing action that is described as a
  harness marker preset (`src/app/ReaderMainWindow.cpp:836`) and REST has a
  harness route for the same mode (`src/app/RestServer.cpp:390`). Keep the REST
  hook if the harness needs it, but hide or rename the toolbar action for
  normal UI.
- Camera checked state is incomplete for Atom/Bond/Subset modes. The code
  knowingly leaves no toolbar action checked for those modes
  (`src/app/ReaderMainWindow.cpp:1003`). Add a visible "Locked" or "Follow"
  state, or expose those modes in the UI if they are real user-facing modes.
- Overlay timing logs are still `qCInfo` every visible frame in both field
  overlays (`src/app/QtFieldGridOverlay.cpp:255`,
  `src/app/QtBFieldStreamOverlay.cpp:288`). Demote to debug or guard behind a
  diagnostics toggle after profiling.
- `QtFieldGridOverlay::setMode()` recomputes frame 0 instead of the current
  frame when visible (`src/app/QtFieldGridOverlay.cpp:263`). That is a
  correctness smell if mode changes can happen away from frame 0.
- `QtBFieldStreamOverlay` dynamic-casts its mapper during every update
  (`src/app/QtBFieldStreamOverlay.cpp:249`). Store the mapper in the ring state
  if this path stays.
- `QtBFieldStreamOverlay` comments reference
  `MoleculeScene::installRenderTimer()` (`src/app/QtBFieldStreamOverlay.cpp:263`),
  but the searched source does not show that as a current scene method. Update
  or remove stale comments.
- `ReaderMainWindow` comments still say the only render verb is widget update
  (`src/app/ReaderMainWindow.cpp:197`), while `MoleculeScene::requestRender()`
  now explicitly calls `iren->Render()`/`renderWindow_->Render()`
  (`src/app/MoleculeScene.cpp:308`). Fix the comment to match the current
  render contract.
- `QtFrame` absent-is-zero accessors should not feed user display without
  availability metadata (`src/model/QtFrame.cpp:1`,
  `src/model/QtFrame.cpp:336`). Mark as legacy, replace with optional/group
  views, or confine it to compatibility-only code.

## Prioritized Surgery List

1. Fix camera clipping on free camera gestures.
   Add scene-owned cached padded bounds and a clipping sync helper; call it
   after all camera writes before render. This should directly address the
   zoom cutoff.

2. Add loaded-run field availability.
   Use it to hide dash-only inspector groups, filter metric descriptors and
   channels, gate default dashboard signals, and prevent all-gap strips from
   being offered as if useful.

3. Make panel recovery obvious.
   Reuse the existing `toggleViewAction()` actions in a `View -> Panels` menu
   and a clearer toolbar/panel control. Do not rely only on far-right toolbar
   text actions.

4. Simplify transform UI semantics.
   Rename or fix Center COM, make Fit backbone the primary stabilizing mode,
   and add reference-frame control if fit modes are intended for normal users.

5. Remove or quarantine harness/probe code from the visible UI.
   Start with `Instrument`, dormant `QtSelectionOverlay`, `/proc/self/statm`,
   stale render comments, frame-info overlay logs, and the frame-0 recompute
   in `QtFieldGridOverlay::setMode()`.

## Answer To The Five Lead Concerns

1. Display transforms: the core architecture is sound. It is one wrapper at the
   `atomPosition()` seam. The main risks are mislabeled Center COM, too-prominent
   all-atom Fit reference, no UI reference-frame control, and unclear
   source-frame vs display-frame labeling.
2. Camera clipping cutoff: yes, there is a concrete bug path. Free camera zoom
   changes the camera but does not reset clipping after automatic clipping was
   disabled.
3. Empty fields shown: yes. Inspector groups can be dash-only, and metric
   descriptors are filtered by implemented sampler, not by data actually present
   and meaningful in the loaded run.
4. Missing pane-toggle buttons: not literally missing, but under-exposed.
   Toolbar toggles exist; they need a proper View/Panels affordance.
5. Accreted ad-hoc cruft: yes. There is enough harness/probe/legacy code in
   the UI path to justify cleanup after the two real correctness fixes above.
