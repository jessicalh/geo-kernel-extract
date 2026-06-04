# UI State Overview - 2026-06-04

Read-only survey of the `h5-reader` Qt/VTK UI side, grounded in:

- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/3d-vtk.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/cmake.md`
- `CLAUDE.md`

No build, test, git command, or code edit was run for this survey. The only
write is this note.

## Executive State

The real UI state is substantially beyond a prototype. The reader now has a
typed `.LGS` load path, a `QtProtein` / `Conformation` model spine, a VTK
molecule scene, playback, atom picking, ordered selection, camera composer
modes, transform modes, ribbon/ring/field overlays, inspector and dashboard
docks, and a loopback REST harness. The camera refactor described in the
late-May viewport notes is largely built: Qt owns input filtering, VTK trackball
is quiet, camera writes are absolute per frame, and reveal/selection camera
modes are typed.

The main missing UI capability for the three stated goals is not camera focus
but local display isolation. The code can focus or stabilize the camera around
selected atoms; it does not yet hide nonlocal atoms/residues or present a
radius control. Linux is source-buildable by documented presets, but not
installable as an adviser-ready artifact: there are no `install()` rules,
package presets, CPack config, AppImage path, or `linuxdeployqt` bundle.

## DOC State

### Designed And Built

`README.md` is now closer to the current implementation than the older design
notes. It describes a standalone Qt6/VTK reader that opens a calcset directory
or `.LGS`, dispatches to trajectory or single-pose conformations, supports the
GUI open action, and documents Linux/macOS/Windows build presets
(`README.md:9`, `README.md:196`, `README.md:205`). It also records platform
status: Linux x86_64 is verified, Linux ARM64/macOS/Windows presets exist but
are not verified (`README.md:318` onward).

`notes/SCOPE.md` still captures the intended application boundary: standalone,
read-only, no dependency on the main `nmr_shielding` library, a typed object
model, 3D playback/static pose, ordered selection up to four atoms, and the
field-grid physics exception. That matches the current code at a high level,
but some data-loading details predate the `.LGS` manifest path.

The late-May dashboard and panel plans are partly built. Current source
contains the panel stack, `DashboardStripDock`, `DashboardDisplayController`,
`SequenceBarPanel`, `PowerSpectrumPanel`, `LagDecayPanel`,
`ChordCouplingPanel`, `FixedFreqPanel`, and scene reveal tensor glyph plumbing
(`ReaderMainWindow.cpp:436`, `ReaderMainWindow.cpp:442`). This means several
"later" items in `notes/PLAN_LATER_ITEMS_2026-05-29.md` have landed:
`lookupBondVector`, Reorient tensor reveal, fixed-frequency panels, sequence
bar, and owned-panel status accounting.

`notes/ROBUSTNESS_BACKLOG_2026-05-30.md` is also partly stale in a good way.
The H5 reader now has file-local exception handling helpers and per-reader
try/catch boundaries (`src/io/QtTrajectoryH5.cpp:61`, `src/io/QtTrajectoryH5.cpp:70`,
`src/io/QtTrajectoryH5.cpp:200`), and required positions validation throws
instead of silently fabricating missing coordinates (`src/io/QtTrajectoryH5.cpp:88`).
The slider scrub backlog item is reflected in the current window wiring
(`ReaderMainWindow.cpp:448`). QSettings persistence is present
(`ReaderMainWindow.cpp:509`).

### Designed But Not Built

The radius/local-isolation view is designed but not implemented. The clearest
record is `notes/VIEWPORT_OBSERVATIONS_2026-05-30.md`, which explicitly calls
distance-filtered atom display "missing entirely" and recommends a
`displayMolecule_`, a source-to-display index map, and an `AtomFilter`
variant such as `WithinRadiusOf(atom, R)`. The same note warns against trying
to use per-atom visibility on `vtkOpenGLMoleculeMapper` because bonds become
orphaned. A source search found no current `displayMolecule_`,
`WithinRadiusOf`, local radius filter, isolate mode, or neighborhood endpoint
in `src/app`.

`notes/STABILIZED_LOCAL_VIEW_DESIGN_2026-05-29.md` is a separate camera/data
transform design: build a local reference frame from selected atoms and keep
that local context visually steady. The current `TransformedConformation` and
`CameraComposer` cover much of the stabilization side, but not the display
isolation side.

`notes/POLISH_BACKLOG.md` remains the best list of UI polish debt. Some items
are built or partially built, but the record still correctly names open areas:
field-grid extent/isovalue controls, glossary/help text for fields, per-metric
atom colouring, inspector in-place update polish, replacing the dynamic mapper
cast in the B-field overlay, cross-platform pass, Windows minidump, and the
residual render-drop investigation.

### Deferred Or Out Of Scope

`notes/VISION_AND_PROGRESS.md` is useful as a historical progress ledger but
predates current source in several places. It still names deliberately deferred
features: Python model mechanism, movie/export paths, GNN/network streaming,
and broader microsecond-tier ambitions. Those should not be treated as desk-
ready requirements.

`notes/BUILD_LAYOUT_PLAN_2026-05-23.md` explicitly defers packaging:
`packagePresets` "lands later" (`notes/BUILD_LAYOUT_PLAN_2026-05-23.md:68`) and
CPack/installer generation is out of scope until cross-platform builds are
green (`notes/BUILD_LAYOUT_PLAN_2026-05-23.md:206`). That matches current code.

## CODE State And Architecture

### Startup And Load

`src/main_reader.cpp` mostly follows the Qt/VTK startup discipline. It forces
desktop OpenGL before `QApplication` (`src/main_reader.cpp:67`), installs the
VTK default surface format before `QApplication` (`src/main_reader.cpp:71`),
creates the app (`src/main_reader.cpp:74`), installs the structured logger
before normal app logs (`src/main_reader.cpp:83`), sets VTK SMP to STDThread
(`src/main_reader.cpp:132`), warms `ErrorBus` / `ObjectCensus`, installs signal
handlers (`src/main_reader.cpp:142`), loads through
`QtProteinLoader::LoadRunPath` (`src/main_reader.cpp:193`), and then constructs
`ReaderMainWindow` (`src/main_reader.cpp:204`).

Load is now `.LGS`-disciplined. `QtProteinLoader.cpp` documents one entry path:
manifest first, no fallback sniffing (`src/io/QtProteinLoader.cpp:11`), and the
trajectory path constructs `QtTrajectoryH5` inside a try/catch at the external
library boundary (`src/io/QtProteinLoader.cpp:92`). `LoadFromManifest` dispatches
trajectory, single pose, and mutant pair modes (`src/io/QtProteinLoader.cpp:210`).

### Model And IO Spine

The model spine is coherent:

- `QtProtein` owns atoms, atom names, residues, topology, rings, and sidecar
  identity.
- `Conformation` abstracts positions and snapshots.
- `TrajectoryConformation` reads animated positions and per-frame NPY snapshots.
- `SingleConformation` presents one pose.
- `TransformedConformation` wraps another conformation and applies runtime
  rigid transforms through the same `atomPosition` seam
  (`src/model/TransformedConformation.cpp:24`, `src/model/TransformedConformation.cpp:51`).
- `QtFrame` is the per-frame atom view surface used by inspector/dashboard
  consumers.

The transform layer is important for UI behavior. It supports identity, center
of mass, all-atom fit, and backbone subset fit (`src/model/TransformedConformation.cpp:78`,
`src/model/TransformedConformation.cpp:153`). It forwards full snapshots
unchanged so calculator data and transformed display positions do not diverge
silently (`src/model/TransformedConformation.cpp:68`).

The H5 seam is now closer to the requested Qt discipline than older notes
suggest. `QtTrajectoryH5.cpp` keeps HighFive local to the `.cpp`, wraps reader
failures into warnings/null optional buffers, and hard-errors required position
shape mismatches (`src/io/QtTrajectoryH5.cpp:1`, `src/io/QtTrajectoryH5.cpp:88`).

### Window, Controllers, And Docks

`ReaderMainWindow` is the UI coordinator. Its constructor registers with the
census (`src/app/ReaderMainWindow.cpp:175`), builds UI/toolbars/status, wraps
the loaded conformation in `TransformedConformation`
(`src/app/ReaderMainWindow.cpp:184`), builds `MoleculeScene`
(`src/app/ReaderMainWindow.cpp:203`), creates playback/time viewport
controllers (`src/app/ReaderMainWindow.cpp:219`), wires playback to scene and
docks (`src/app/ReaderMainWindow.cpp:223`), creates the picker
(`src/app/ReaderMainWindow.cpp:261`), creates the Qt camera input filter
(`src/app/ReaderMainWindow.cpp:269`), and installs inspector, selection, and
dashboard docks (`src/app/ReaderMainWindow.cpp:275`, `src/app/ReaderMainWindow.cpp:376`,
`src/app/ReaderMainWindow.cpp:386`).

The window also exposes the main demo controls:

- Playback, frame slider, and FPS (`src/app/ReaderMainWindow.cpp:710`,
  `src/app/ReaderMainWindow.cpp:725`, `src/app/ReaderMainWindow.cpp:731`).
- Camera actions: Focus, Newman, Plane lock, Free
  (`src/app/ReaderMainWindow.cpp:747`, `src/app/ReaderMainWindow.cpp:757`,
  `src/app/ReaderMainWindow.cpp:767`, `src/app/ReaderMainWindow.cpp:776`).
- Transform menu: Identity, Center COM, Fit reference, Fit backbone
  (`src/app/ReaderMainWindow.cpp:790`, `src/app/ReaderMainWindow.cpp:794`).
- Instrument marker toggle (`src/app/ReaderMainWindow.cpp:827`).
- Metrics dialog action (`src/app/ReaderMainWindow.cpp:839`).
- Overlay toggles: Ribbon, Rings, Butterfly, B-field
  (`src/app/ReaderMainWindow.cpp:850`, `src/app/ReaderMainWindow.cpp:856`,
  `src/app/ReaderMainWindow.cpp:862`, `src/app/ReaderMainWindow.cpp:869`).

Shutdown is explicit and Qt/VTK-aware: it asserts GUI thread, deletes the REST
server synchronously, stops child timers, and detaches the render window from
the widget before GL teardown (`src/app/ReaderMainWindow.cpp:633`,
`src/app/ReaderMainWindow.cpp:646`, `src/app/ReaderMainWindow.cpp:658`,
`src/app/ReaderMainWindow.cpp:671`).

### Scene, VTK Pipeline, And Overlays

`MoleculeScene` owns the render-side architecture. It sets up two renderers,
main and overlay, sharing one camera (`src/app/MoleculeScene.cpp:91`,
`src/app/MoleculeScene.cpp:111`, `src/app/MoleculeScene.cpp:115`), configures a
two-layer render window with alpha bit planes and no MSAA
(`src/app/MoleculeScene.cpp:122`), installs `QuietTrackballStyle`
(`src/app/MoleculeScene.cpp:129`), and logs one render EndEvent line with source,
frame, render time, and camera mode (`src/app/MoleculeScene.cpp:140`).

Build constructs a single source `vtkMolecule`, appends all atoms and static
bonds (`src/app/MoleculeScene.cpp:192`, `src/app/MoleculeScene.cpp:197`,
`src/app/MoleculeScene.cpp:205`), feeds it to `vtkOpenGLMoleculeMapper`
(`src/app/MoleculeScene.cpp:214`), and then builds overlays in place:

- `QtBackboneRibbonOverlay` (`src/app/MoleculeScene.cpp:237`)
- `QtRingPolygonOverlay` (`src/app/MoleculeScene.cpp:242`)
- `QtFieldGridOverlay` (`src/app/MoleculeScene.cpp:247`)
- `QtBFieldStreamOverlay` (`src/app/MoleculeScene.cpp:252`)
- `QtSelectionOverlay` (`src/app/MoleculeScene.cpp:257`)
- `MeasurementOverlay` on the overlay renderer (`src/app/MoleculeScene.cpp:267`)
- `SceneRevealOverlay` on the overlay renderer (`src/app/MoleculeScene.cpp:272`)

Frame updates are one fan-out point. `setFrame` pushes atom positions, bumps
VTK modification time, writes camera state through the composer, updates every
overlay, resets clipping bounds, and schedules one render
(`src/app/MoleculeScene.cpp:354`, `src/app/MoleculeScene.cpp:370`,
`src/app/MoleculeScene.cpp:386`, `src/app/MoleculeScene.cpp:399`,
`src/app/MoleculeScene.cpp:408`, `src/app/MoleculeScene.cpp:428`,
`src/app/MoleculeScene.cpp:431`).

Render scheduling is coalesced and GUI-threaded through a queued Qt invocation
(`src/app/MoleculeScene.cpp:299`). The implementation currently calls
`iren->Render()` or `renderWindow_->Render()` inside that queued lambda
(`src/app/MoleculeScene.cpp:314`), because comments say `widget->update()` alone
blits a stale FBO (`src/app/MoleculeScene.cpp:308`). That is not a raw VTK
event-loop takeover, but it is a deliberate departure from the earlier comment
in `ReaderMainWindow` that says `vtkWidget_->update()` is the only render verb
(`src/app/ReaderMainWindow.cpp:197`).

Overlay state:

- Ribbon and ring polygons are ordinary VTK actor/pipeline overlays and are on
  by default.
- Butterfly field grid and B-field streamlines are off by default and
  expensive-when-visible (`src/app/ReaderMainWindow.cpp:862`,
  `src/app/ReaderMainWindow.cpp:869`).
- Field grid exposes mode/threshold/opacity setters, but the toolbar only wires
  visibility (`src/app/QtFieldGridOverlay.cpp:263`,
  `src/app/QtFieldGridOverlay.cpp:282`, `src/app/QtFieldGridOverlay.cpp:288`).
- B-field streamlines still use a dynamic cast to the mapper
  (`src/app/QtBFieldStreamOverlay.cpp:249`), matching an open polish item.

### Camera

The camera architecture is one of the strongest parts of the current UI.
`CameraMode` defines typed modes: Free, Atom, Bond, Dihedral, Plane, and Subset.
`CameraComposer` owns absolute per-frame camera writes and mode changes. `Plane`
releases on selection change; Atom/Bond/Dihedral/Subset are sustained. The scene
uses this for explicit toolbar actions and for dashboard reveal focus
(`src/app/MoleculeScene.cpp:562`, `src/app/MoleculeScene.cpp:581`,
`src/app/MoleculeScene.cpp:587`, `src/app/MoleculeScene.cpp:590`,
`src/app/MoleculeScene.cpp:593`).

`CameraInputFilter` is a Qt event filter on the QVTK widget. That matches the
skill guidance: Qt owns interaction; VTK trackball is present only to satisfy
adapter plumbing and is quieted by `QuietTrackballStyle`. This avoids the old
"VTK trackball fights app camera" problem recorded in the viewport notes.

The camera path should be considered built, not just designed. What is missing
is UI affordance clarity for modes reached indirectly through reveal: the
toolbar has only Focus/Newman/Plane/Free, while reveal can set Atom/Bond/Subset.
That can leave the camera in a real sustained mode with no matching toolbar
action checked.

### Selection, Picking, Inspector, Dashboard, REST

`QtAtomPicker` uses a Qt event filter on double-click, DPR-aware coordinates,
and a nearest-atom ray cast. It emits picks and stays intentionally dumb; the
`AtomSelection` model interprets plain versus Shift modifiers
(`src/app/ReaderMainWindow.cpp:253`, `src/app/ReaderMainWindow.cpp:315`).

`AtomSelection` is the single source of selection truth. It drives the inspector
focus, measurement overlay, metric action enablement, plane-lock release, and
scene refresh (`src/app/ReaderMainWindow.cpp:331`, `src/app/ReaderMainWindow.cpp:350`,
`src/app/ReaderMainWindow.cpp:355`). The selection cap of four atoms supports
distance/angle/dihedral operation but is not a general set-selection model.

`MeasurementOverlay` is the active selected-atom visual layer. `QtSelectionOverlay`
is still built, but it is effectively legacy/dormant now that the measurement
overlay owns selected atom markers.

The inspector and dashboard are wired well enough for a desk demo:
`QtAtomInspectorDock` receives focus atom and current frame; `DashboardStripDock`
receives model/context/selection/time viewport; the dashboard can reveal scene
bindings and has a metrics dialog path. The code also seeds an initial generic
`npy:dssp_chi` dashboard signal (`src/app/ReaderMainWindow.cpp:306`), which is
useful for smoke visibility but may be surprising in a polished workflow.

`RestServer` is a loopback GUI-thread test surface. It exposes health, frame,
selection, instrument mode, dock visibility, transform, plane lock, camera mode,
camera focus, positions, dashboard signals, screenshot, and shutdown. It does
not expose a radius/local-isolation view.

## Qt Discipline Assessment

### Follows The Discipline

- `CENSUS_REGISTER(this)` is present in core QObjects reviewed:
  `ReaderMainWindow`, `MoleculeScene`, `TransformedConformation`, camera and
  overlay classes.
- `ASSERT_THREAD(this)` is used at thread-sensitive UI/model mutation points,
  including `ReaderMainWindow::shutdown`, `MoleculeScene::Build`,
  `MoleculeScene::setFrame`, overlay `setFrame` methods, and transform mode
  changes.
- Most signal wiring uses `ACONNECT`, including transform, playback, picker,
  selection, docks, and toolbar actions.
- Startup follows the VTK/Qt ordering: `QSurfaceFormat` before `QApplication`,
  `QVTKOpenGLNativeWidget`, `vtkGenericOpenGLRenderWindow`, and
  `vtk_module_autoinit` (`CMakeLists.txt:278`).
- VTK work stays on the GUI side; overlays update backing data and do not
  schedule their own renders.
- The scene uses VTK smart pointers and a source/filter/mapper/actor shape.
  The ribbon overlay uses a proper VTK filter path; field/stream overlays own
  mutable data objects and force expensive pipeline updates only when visible.
- UDP/structured logging is installed before ordinary app logs
  (`src/main_reader.cpp:83`), and render logs carry a source tag.
- Persistent `QTimer` ownership appears concentrated in `QtPlaybackController`.
  The zero-delay `QTimer::singleShot` calls are single-shots, so they do not
  violate the one-persistent-timer rule.
- HighFive exceptions are caught at the loader/reader boundary, and optional
  malformed TR groups degrade to warnings/null buffers instead of unwinding into
  Qt's event loop.
- CMake has platform helper modules and preset indirection, matching the skill's
  deployment/build hygiene for multi-platform Qt/VTK projects.

### Violations, Risks, And Drift

- Raw connect: `main_reader.cpp` uses `QObject::connect` for `aboutToQuit` to
  `ReaderMainWindow::shutdown` (`src/main_reader.cpp:205`). Project discipline
  says use `ACONNECT`.
- Linux-only file read in app code: `MoleculeScene::setFrame` reads
  `/proc/self/statm` directly (`src/app/MoleculeScene.cpp:448`). That violates
  the project's cross-platform rule against Linux-only system calls outside
  diagnostics/platform modules. It also contradicts the polish backlog's "avoid
  Linux-only system calls outside diagnostics" guidance.
- Render-path comment drift: `ReaderMainWindow` says the scheduler can call
  `vtkWidget_->update()` as the only render verb (`src/app/ReaderMainWindow.cpp:197`),
  but `MoleculeScene::requestRender` now calls `iren->Render()` / `Render()`
  inside a queued lambda (`src/app/MoleculeScene.cpp:314`). The current behavior
  may be correct, but the documentation trail needs consolidation.
- VTK render-drop probe still in production path: `points->Modified()` is
  labelled a probe pending confirmation (`src/app/MoleculeScene.cpp:372`).
  Either settle it as the fix or retire it after the render-drop note is closed.
- Dormant overlay: `QtSelectionOverlay` is still built
  (`src/app/MoleculeScene.cpp:257`) even though picker selection is wired to
  `MeasurementOverlay`. This is harmless but confusing cruft.
- B-field mapper update still uses `dynamic_cast<vtkPolyDataMapper*>`
  (`src/app/QtBFieldStreamOverlay.cpp:249`), exactly the typed-mapper cleanup
  item in the polish backlog.
- Stale B-field comment references `MoleculeScene::installRenderTimer()`
  (`src/app/QtBFieldStreamOverlay.cpp:263`), but the current scene uses an
  EndEvent observer instead.
- Field-grid `setMode` recomputes frame zero when visible
  (`src/app/QtFieldGridOverlay.cpp:263`, `src/app/QtFieldGridOverlay.cpp:270`).
  There is no current UI control that trips this, but it is wrong for a future
  mode selector unless current frame is passed or stored.
- Placeholder model value: `QtFrame::eeqCharge` returns `0.0` regardless of the
  Welford buffer (`src/model/QtFrame.cpp:336`, `src/model/QtFrame.cpp:342`).
  That can mislead inspector/dashboard users if EEQ fields are shown.
- Retired registry files remain on disk and are called out by CMake as
  vestigial (`CMakeLists.txt:138`). That is acceptable but should not be
  rediscovered as active architecture.
- Toolbar wiring defers several connections with a zero-delay single-shot
  because recipients are constructed later (`src/app/ReaderMainWindow.cpp:847`,
  `src/app/ReaderMainWindow.cpp:881`). It is allowed by the timer rule, but it
  is a construction-order smell and makes toolbar behavior harder to audit.

## Cruft Register

- `src/app/QtSelectionOverlay.*`: legacy single-picked-atom overlay. It is still
  constructed and frame-updated, but the selection UX now runs through
  `AtomSelection` + `MeasurementOverlay`.
- `src/app/MoleculeScene.cpp:372`: render-drop "PROBE" comment and explicit
  points MTime bump. Needs resolution into a settled fix or removal.
- `src/app/MoleculeScene.cpp:448`: `/proc/self/statm` in app rendering code.
  Move behind diagnostics/platform abstraction or delete.
- `src/app/QtBFieldStreamOverlay.cpp:249`: dynamic mapper cast; store the typed
  mapper in the ring state.
- `src/app/QtBFieldStreamOverlay.cpp:263`: stale `installRenderTimer` comment.
- `src/app/QtFieldGridOverlay.cpp:270`: visible mode switch recomputes frame 0.
- `src/model/QtFrame.cpp:342`: EEQ placeholder returns zero even when the
  Welford buffer exists.
- `CMakeLists.txt:138`: retired `QtNamingRegistry` files remain on disk.
- `src/app/ReaderMainWindow.cpp:306`: default `npy:dssp_chi` dashboard signal is
  useful for smoke state but may be demo clutter unless made intentional.
- `src/app/ReaderMainWindow.cpp:829`: "Instrument" is a harness concept exposed
  as a toolbar label; useful internally, awkward as user-facing UI.

## Against The Three Goals

### 1. Open 1P9J, Work, And Not Suck Tomorrow

Likely state: if the 1P9J calcset has a valid `.LGS` and the host has the
documented Qt/VTK/HDF5 runtime paths, the reader should open and present a
usable scene. The code path for `.LGS` loading, typed topology, trajectory
positions, playback, picker, selection, inspector, transform, camera modes,
and dashboards is built. I did not run the app, so this is a code/doc
assessment, not a runtime verification.

The operational infelicities most likely to grate:

- The toolbar is dense and text-heavy: Focus, Newman, Plane lock, Free,
  Transform, Instrument, Metrics, Ribbon, Rings, Butterfly, B-field, plus
  playback controls. It is functional, but it reads as a debug/operator toolbar
  more than a polished viva surface.
- Selection affordance is underexplained in operation. Double-click selects;
  Shift-double-click adds/toggles; the Selection dock lists atoms, but there is
  no obvious Clear button or click-to-focus list affordance.
- "Instrument" is exposed as a main action even though it is a harness/marker
  preset. For a user, "Highlight" or a visual focus toggle would be clearer.
- Reveal-driven Atom/Bond/Subset camera modes can be active without a matching
  checked toolbar action. The camera can feel locked while the toolbar only
  shows none of the expected radio buttons.
- Field overlays are binary toggles. There is no visible mode/threshold/extent
  control for the butterfly, despite setters existing in code.
- Butterfly and B-field overlays log at info level for every visible frame
  (`src/app/QtFieldGridOverlay.cpp:255`, `src/app/QtBFieldStreamOverlay.cpp:288`).
  That is good for profiling but noisy for a clean demo.
- The default `npy:dssp_chi` dashboard signal may be surprising if the first
  screen should be quiet.
- Inspector/dashboard values still have placeholder risk for some fields
  (`QtFrame::eeqCharge`) and older backlog items around in-place update polish.
- There is no local radius view, which is the most important missing
  "look here, ignore the rest" interaction for explaining 1P9J.

Desk-ready tomorrow means the safest path is to demo existing strengths:
load `.LGS`, apply `Fit backbone` or `Center COM` if the molecule drifts,
use atom picking plus measurement markers, use Focus/Newman/Plane lock for
geometry, keep expensive overlays off until needed, and avoid promising radius
isolation unless implemented in a focused patch.

### 2. Radius / Local-Isolation View

Current code has camera localization, not molecule isolation. The relevant
built pieces are:

- Focus atom and ordered selection via `AtomSelection`.
- Position seam through `Conformation::atomPosition`.
- Camera modes Atom/Bond/Dihedral/Plane/Subset.
- `TransformedConformation` for global/stabilized coordinate transforms.
- Measurement/reveal overlays that operate on source atom indices.

The missing piece is a display filter for the molecule actor. Because
`vtkOpenGLMoleculeMapper` does not offer a clean per-atom visibility path
without bond artifacts, the path recorded in the notes is still the right one:

1. Add an `AtomFilter` state to `MoleculeScene`: `All` and
   `WithinRadiusOf(focusAtom, radius, residueExpanded=true)`.
2. Add a `displayMolecule_` separate from the source molecule plus maps:
   `sourceAtom -> displayAtom or npos` and optionally `displayAtom -> sourceAtom`.
3. On filter changes, compute membership from the current transformed
   conformation. For the requested "residues within R" behavior, include every
   residue that has any atom within radius of the picked atom, plus the picked
   atom's residue unconditionally.
4. Rebuild `displayMolecule_` with `AppendAtom` for included atoms and
   `AppendBond` only when both endpoints are included. Point the existing
   molecule mapper at `displayMolecule_`.
5. On each frame, update only included display atom positions. If membership is
   dynamic around the moving focus atom, recompute membership and rebuild only
   when the included set changes.
6. Keep selection/measurement/reveal overlays source-index based. They can
   continue to draw highlighted atoms even when the molecule actor is filtered,
   as long as overlay position lookup still uses source indices.
7. Add UI controls close to selection: an isolate toggle, a radius double
   spinbox, and a clear action. Enable the toggle only when there is a focus
   atom. Use icons/tooltips where possible rather than another long toolbar
   word.
8. Add REST support for harness/demo repeatability, for example
   `POST /local-view { atom, radius, residue_expanded }` and
   `POST /local-view/clear`.

This is a narrow scene/window feature, not a model rewrite. It should not
touch H5 loading, dashboard sampling, or the calculator source model. The
risky part is preserving picker/selection semantics while the displayed
molecule no longer has source atom indices; the clean mitigation is to leave
picker ray-cast over source positions, not over mapper atom ids, or maintain
the display/source index map explicitly.

### 3. Linux Installability

Linux buildability exists; Linux installability does not.

What exists:

- Linux presets for debug, RelWithDebInfo, and release
  (`CMakePresets.json:29`, `CMakePresets.json:35`, `CMakePresets.json:41`).
- Documented Linux dependencies and source-build command
  (`README.md:71`).
- Linux platform module for HDF5 include-order and compile settings.
- README troubleshooting acknowledges VTK runtime library path problems:
  missing VTK runtime libs can produce a black window, and users may need
  `LD_LIBRARY_PATH` for `~/VTK/lib` (`README.md:278`, `README.md:285`).

What is missing:

- No `install()` rules in `CMakeLists.txt`.
- No `packagePresets` in `CMakePresets.json`.
- No CPack config.
- No AppImage flow.
- No `linuxdeployqt` or equivalent Qt plugin/library bundling flow.
- No RPATH/install-name strategy for a self-contained Linux directory.
- No clean-machine smoke doc for an adviser binary.

The qt6-cpp skill is Windows-targeted for deployment details, so it cannot be
treated as sufficient for Linux packaging. For Linux, the practical next step
is a Linux-specific package track:

1. Add `install(TARGETS h5reader RUNTIME DESTINATION bin)` plus resource/plugin
   destinations if needed.
2. Set runtime search paths so bundled libraries resolve from `$ORIGIN` and
   `$ORIGIN/../lib`.
3. Choose one artifact path for tomorrow-level use: AppImage is simplest for a
   "download and run" adviser binary; CPack TGZ/DEB is cleaner if the target
   machine is known.
4. Bundle Qt platform plugins, Qt shared libraries, VTK shared libraries, HDF5,
   and any non-system runtime libraries. `linuxdeployqt` can do much of the Qt
   plugin work, but VTK and HDF5 still need explicit verification.
5. Test on a clean Ubuntu 24.04 machine or container with GUI/OpenGL access.
   The current README source-build path is not a substitute for this.

## Bottom Line

The UI is a real Qt/VTK application with a good model spine, a disciplined scene
pipeline, a mostly modern camera architecture, and useful docks/overlays. It is
not yet a polished local-explanation instrument: interaction labels are still
debug-flavored, some stale/dormant code remains, field overlays lack controls,
and the radius isolation feature is absent. The fastest radius path is a
scene-level display molecule plus an `AtomFilter::WithinRadiusOf` state driven
by the current selection focus and a small UI/REST surface. Linux is currently
source-buildable, not installable; packaging needs an explicit Linux artifact
plan outside the Windows-oriented skill discipline.
