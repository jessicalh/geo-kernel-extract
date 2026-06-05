# Visualisation Audit And Must-Haves - 2026-06-05

Read-only audit and proposal. Source scope was `h5-reader/src` only, plus the prior inventory note as a map. No source code was changed, no git commands were run, no H5 write path was touched, and no extraction/build/run was triggered.

Prior map checked: [notes/DISPLAY_MODE_INVENTORY_2026-06-05.md](DISPLAY_MODE_INVENTORY_2026-06-05.md).

External VTK references consulted for the VTK-specific rendering/pipeline statements:

- QVTK/Qt integration: https://vtk.org/doc/nightly/html/classQVTKOpenGLNativeWidget.html
- OpenGL render window: https://vtk.org/doc/nightly/html/classvtkGenericOpenGLRenderWindow.html
- Renderer/render-window layering and render calls: https://vtk.org/doc/nightly/html/classvtkRenderer.html and https://vtk.org/doc/nightly/html/classvtkRenderWindow.html
- Actor visibility/mapper use: https://vtk.org/doc/nightly/html/classvtkActor.html
- Glyph/transform pipeline: https://vtk.org/doc/nightly/html/classvtkTransformPolyDataFilter.html
- Ribbon, contour, streamline filters: https://vtk.org/doc/nightly/html/classvtkProteinRibbonFilter.html, https://vtk.org/doc/nightly/html/classvtkContourFilter.html, https://vtk.org/doc/nightly/html/classvtkStreamTracer.html

## Executive Finding

The lead's read is materially right. There is a small drawing-widget hierarchy after dispatch, but there is not a real typed visualisation class system. The product-level identity of a visualisation is mostly a `QString` mode id carried through catalog descriptors, dialog checkboxes, panel refs, and controller dispatch. The typed classes exist at the end of the pipeline as panels and overlays; they do not own whether a mode is offerable, what data it needs, or how it renders.

What reliably renders today:

- Temporal strip panels for `strip.*` modes when the controller can sample the signal.
- Several hard-coded static dashboard panels: sequence bars, lag curves, chord coupling, fixed-frequency curves, and one power-spectrum builder that is not normally selectable through the visible-mode gate.
- Main 3D molecular scene and scene overlays: ribbon, rings, field contours, B-field streamlines, measurement overlay, reveal overlay.

What is hollow or advertised-dead:

- `static.table`, `static.atomColor`, `static.tensor`, `static.scalar`, `static.efg`, `static.vector`, `static.vectorGlyph`, `static.category`, `static.per-class`, `static.rollup`, `static.event`, `static.system`, `static.embedding`, `static.geometry`, `static.topology`.
- `strip.spectrum` as an actual spectrum strip: the code has a spectrum panel class, but the live dock path clears spectrum tracks and sends all `strip.*` modes through temporal strips.
- Several component/category/event modes that technically draw pixels but do not provide the control or semantics their names imply.

## 1. Class System - Is There One?

Verdict: no, not at the visualisation-type level. There is a reusable panel base class and several concrete panel classes, but the system that decides what a visualisation is, whether it can be offered, and how it is routed is string dispatch.

### What Is Typed

There is a real panel base class for stackable 2D panels. `AbstractStripPanel` defines common panel behaviour and pure virtual paint/reveal hooks ([src/app/AbstractStripPanel.h:145](../src/app/AbstractStripPanel.h:145), [src/app/AbstractStripPanel.h:150](../src/app/AbstractStripPanel.h:150), [src/app/AbstractStripPanel.h:155](../src/app/AbstractStripPanel.h:155)). Shared painting helpers live under that base, including header/reveal UI drawing ([src/app/AbstractStripPanel.cpp:160](../src/app/AbstractStripPanel.cpp:160), [src/app/AbstractStripPanel.cpp:131](../src/app/AbstractStripPanel.cpp:131)).

Concrete subclasses exist:

- `SequenceBarPanel final : public AbstractStripPanel` ([src/app/SequenceBarPanel.h:43](../src/app/SequenceBarPanel.h:43)), with clickable reveal binding from sequence bars ([src/app/SequenceBarPanel.cpp:124](../src/app/SequenceBarPanel.cpp:124)).
- `PowerSpectrumPanel final : public AbstractStripPanel` ([src/app/PowerSpectrumPanel.h:24](../src/app/PowerSpectrumPanel.h:24)).
- `LagDecayPanel final : public AbstractStripPanel` ([src/app/LagDecayPanel.h:26](../src/app/LagDecayPanel.h:26)).
- `ChordCouplingPanel final : public AbstractStripPanel` ([src/app/ChordCouplingPanel.h:21](../src/app/ChordCouplingPanel.h:21)).
- `FixedFreqPanel final : public AbstractStripPanel` ([src/app/FixedFreqPanel.h:27](../src/app/FixedFreqPanel.h:27)).
- Private stack panels `TemporalStripPanel` and `SpectrumStripPanel` also subclass `AbstractStripPanel` ([src/app/StripStackWidget.cpp:60](../src/app/StripStackWidget.cpp:60), [src/app/StripStackWidget.cpp:288](../src/app/StripStackWidget.cpp:288)).

There are also typed scene overlay classes. `MoleculeScene` owns typed overlay objects for ribbon, rings, field grid, B-field streams, measurement, and reveal ([src/app/MoleculeScene.h:154](../src/app/MoleculeScene.h:154), [src/app/MoleculeScene.h:228](../src/app/MoleculeScene.h:228)). It builds them explicitly ([src/app/MoleculeScene.cpp:244](../src/app/MoleculeScene.cpp:244), [src/app/MoleculeScene.cpp:249](../src/app/MoleculeScene.cpp:249), [src/app/MoleculeScene.cpp:254](../src/app/MoleculeScene.cpp:254), [src/app/MoleculeScene.cpp:259](../src/app/MoleculeScene.cpp:259), [src/app/MoleculeScene.cpp:269](../src/app/MoleculeScene.cpp:269), [src/app/MoleculeScene.cpp:274](../src/app/MoleculeScene.cpp:274)).

### What Is Not Typed

The visualisation identity is not a typed object. It is a mode string.

The model stores display modes as strings: `DisplaySignalBinding` has `QString displayModeId` ([src/model/DashboardSignal.h:266](../src/model/DashboardSignal.h:266)), and `DashboardSignal` has `QStringList displayModeIds` ([src/model/DashboardSignal.h:275](../src/model/DashboardSignal.h:275)). `AllDisplayModes` concatenates descriptor `temporalModes` and `staticModes` strings ([src/model/DashboardSignal.cpp:285](../src/model/DashboardSignal.cpp:285), [src/model/DashboardSignal.cpp:292](../src/model/DashboardSignal.cpp:292)), and `SupportsDisplayMode` checks those same string lists ([src/model/DashboardSignal.cpp:299](../src/model/DashboardSignal.cpp:299)).

The capability table is also string keyed. Every `strip.` mode is marked visible by prefix, not by a strip visualisation class ([src/model/DisplayModeCapability.h:25](../src/model/DisplayModeCapability.h:25)). Only a small set of static string ids have capability rows ([src/model/DisplayModeCapability.h:29](../src/model/DisplayModeCapability.h:29), [src/model/DisplayModeCapability.h:41](../src/model/DisplayModeCapability.h:41)); unknown `static.` modes return no capability ([src/model/DisplayModeCapability.h:47](../src/model/DisplayModeCapability.h:47)).

The dialog also uses string mapping and string matching. UI kinds such as Strip, Spectrum, Table, Color map, Glyph, Bar, Curve, Chord, and Fixed freq are an enum in the dialog ([src/app/SignalDisplayDialog.cpp:56](../src/app/SignalDisplayDialog.cpp:56), [src/app/SignalDisplayDialog.cpp:70](../src/app/SignalDisplayDialog.cpp:70)), but each maps back to a canonical mode string like `strip.scalar`, `strip.spectrum`, `static.table`, `static.atomColor`, or `static.tensor` ([src/app/SignalDisplayDialog.cpp:156](../src/app/SignalDisplayDialog.cpp:156), [src/app/SignalDisplayDialog.cpp:178](../src/app/SignalDisplayDialog.cpp:178)). Matching a descriptor to a UI kind is string-prefix/string-equality logic ([src/app/SignalDisplayDialog.cpp:180](../src/app/SignalDisplayDialog.cpp:180), [src/app/SignalDisplayDialog.cpp:207](../src/app/SignalDisplayDialog.cpp:207)).

The panel model routes display references by string. `DisplayRefsForSignal` calls `DisplayModeCapabilityFor(modeId)` and emits a panel ref when the string's capability says `emitsPanelRef` ([src/model/DashboardPanelModel.cpp:89](../src/model/DashboardPanelModel.cpp:89), [src/model/DashboardPanelModel.cpp:96](../src/model/DashboardPanelModel.cpp:96)). Strip channel routing is duplicated string logic with `canonicalModeChannel` and `modeWantsChannel` ([src/model/DashboardPanelModel.cpp:21](../src/model/DashboardPanelModel.cpp:21), [src/model/DashboardPanelModel.cpp:31](../src/model/DashboardPanelModel.cpp:31)).

The display controller dispatches by mode string and storage path. It defines `isStripMode` as `startsWith("strip.")` ([src/app/DashboardDisplayController.cpp:47](../src/app/DashboardDisplayController.cpp:47)), duplicates channel routing for `strip.tensor.*` and `strip.vector.*` ([src/app/DashboardDisplayController.cpp:75](../src/app/DashboardDisplayController.cpp:75), [src/app/DashboardDisplayController.cpp:90](../src/app/DashboardDisplayController.cpp:90)), gates static panel construction through the capability table ([src/app/DashboardDisplayController.cpp:1605](../src/app/DashboardDisplayController.cpp:1605), [src/app/DashboardDisplayController.cpp:1607](../src/app/DashboardDisplayController.cpp:1607)), then uses explicit `storagePath` plus mode-id branches for each static panel ([src/app/DashboardDisplayController.cpp:1621](../src/app/DashboardDisplayController.cpp:1621), [src/app/DashboardDisplayController.cpp:1670](../src/app/DashboardDisplayController.cpp:1670)).

### How A Metric And Mode Become Pixels Today

1. The catalog creates a `SignalDescriptor` with mode strings. Helper lists produce tensor, scalar, vector, category, event, rollup, and static mode ids ([src/model/TrajectorySignalCatalog.cpp:84](../src/model/TrajectorySignalCatalog.cpp:84), [src/model/TrajectorySignalCatalog.cpp:204](../src/model/TrajectorySignalCatalog.cpp:204)). `makeDescriptor` stores those lists and derives temporal/static flags from list emptiness ([src/model/TrajectorySignalCatalog.cpp:299](../src/model/TrajectorySignalCatalog.cpp:299), [src/model/TrajectorySignalCatalog.cpp:332](../src/model/TrajectorySignalCatalog.cpp:332)).
2. The dialog presents those modes through kind checkboxes and filters. Descriptor rows summarize the available modes and expose raw mode ids in the tooltip ([src/app/SignalDisplayDialog.cpp:437](../src/app/SignalDisplayDialog.cpp:437), [src/app/SignalDisplayDialog.cpp:465](../src/app/SignalDisplayDialog.cpp:465)). Candidate checkboxes are created for all dialog mode kinds ([src/app/SignalDisplayDialog.cpp:801](../src/app/SignalDisplayDialog.cpp:801), [src/app/SignalDisplayDialog.cpp:824](../src/app/SignalDisplayDialog.cpp:824)).
3. The dialog enables a candidate mode only when a descriptor supports the mode string, data is available, and `DisplayModeCapabilityFor(modeId).hasVisibleSurface` is true ([src/app/SignalDisplayDialog.cpp:1191](../src/app/SignalDisplayDialog.cpp:1191), [src/app/SignalDisplayDialog.cpp:1208](../src/app/SignalDisplayDialog.cpp:1208)). Active signal checkboxes use the same capability gate ([src/app/SignalDisplayDialog.cpp:1399](../src/app/SignalDisplayDialog.cpp:1399), [src/app/SignalDisplayDialog.cpp:1416](../src/app/SignalDisplayDialog.cpp:1416)).
4. Adding or toggling a metric sends string mode ids to `DashboardSelectionController`, which stores them in the signal model and asks the panel model for display refs ([src/app/DashboardSelectionController.cpp:58](../src/app/DashboardSelectionController.cpp:58), [src/app/DashboardSelectionController.cpp:87](../src/app/DashboardSelectionController.cpp:87), [src/app/DashboardSelectionController.cpp:120](../src/app/DashboardSelectionController.cpp:120), [src/app/DashboardSelectionController.cpp:164](../src/app/DashboardSelectionController.cpp:164)).
5. `DashboardSignalModel::setData(DisplayModeRole)` mutates the signal's string mode list ([src/model/DashboardSignalModel.cpp:263](../src/model/DashboardSignalModel.cpp:263), [src/model/DashboardSignalModel.cpp:276](../src/model/DashboardSignalModel.cpp:276)).
6. On rebuild, panel-capable static mode ids become concrete `AbstractStripPanel` subclasses through hard-coded controller builders ([src/app/DashboardDisplayController.cpp:1576](../src/app/DashboardDisplayController.cpp:1576), [src/app/DashboardDisplayController.cpp:1670](../src/app/DashboardDisplayController.cpp:1670)).
7. Any signal with a `strip.*` mode goes into generic track construction ([src/app/DashboardDisplayController.cpp:1695](../src/app/DashboardDisplayController.cpp:1695), [src/app/DashboardDisplayController.cpp:1699](../src/app/DashboardDisplayController.cpp:1699)). `buildGenericTracks` loops `signal.displayModeIds`, filters by `isStripMode`, chooses channels, and creates active strip series ([src/app/DashboardDisplayController.cpp:2368](../src/app/DashboardDisplayController.cpp:2368), [src/app/DashboardDisplayController.cpp:2414](../src/app/DashboardDisplayController.cpp:2414)).
8. `DashboardStripDock::refreshTracks` converts controller tracks to `StripStackWidget::Track` objects and calls `setTracks` ([src/app/DashboardStripDock.cpp:297](../src/app/DashboardStripDock.cpp:297), [src/app/DashboardStripDock.cpp:312](../src/app/DashboardStripDock.cpp:312)). `StripStackWidget::setTracks` instantiates `TemporalStripPanel` for each track ([src/app/StripStackWidget.cpp:401](../src/app/StripStackWidget.cpp:401), [src/app/StripStackWidget.cpp:418](../src/app/StripStackWidget.cpp:418)).
9. Static owned panels are pushed from controller to dock to stack; the stack appends them after strips ([src/app/DashboardStripDock.cpp:121](../src/app/DashboardStripDock.cpp:121), [src/app/DashboardStripDock.cpp:125](../src/app/DashboardStripDock.cpp:125), [src/app/StripStackWidget.cpp:441](../src/app/StripStackWidget.cpp:441), [src/app/StripStackWidget.cpp:456](../src/app/StripStackWidget.cpp:456)).
10. Reveal clicks can enter the 3D scene through the dock/main-window connection to `MoleculeScene::revealBinding` ([src/app/ReaderMainWindow.cpp:416](../src/app/ReaderMainWindow.cpp:416), [src/app/MoleculeScene.cpp:519](../src/app/MoleculeScene.cpp:519), [src/app/MoleculeScene.cpp:538](../src/app/MoleculeScene.cpp:538)).

### OO Discipline Judgment

The project has good typed discipline in places: QObject overlays are explicit owned objects, panel classes are concrete, and reveal bindings are typed enough to move between panels and the scene. But visualisation identity itself is not OO. It is a string protocol spread across catalog helpers, capability rows, dialog mode-kind heuristics, panel-ref generation, and controller branch dispatch.

That is why hollow modes can exist. A catalog helper can advertise `static.table` or `static.atomColor` without any visualisation object behind it ([src/model/TrajectorySignalCatalog.cpp:94](../src/model/TrajectorySignalCatalog.cpp:94), [src/model/TrajectorySignalCatalog.cpp:100](../src/model/TrajectorySignalCatalog.cpp:100)), and the only thing preventing the dialog from enabling it is a separate capability table that returns empty for unknown static strings ([src/model/DisplayModeCapability.h:47](../src/model/DisplayModeCapability.h:47)).

## 2. What Actually Renders

### Main 3D Scene

The app uses Qt/VTK integration correctly at the top level: it sets the QVTK default surface format before `QApplication` ([src/main_reader.cpp:69](../src/main_reader.cpp:69), [src/main_reader.cpp:74](../src/main_reader.cpp:74)), creates a `QVTKOpenGLNativeWidget` and a `vtkGenericOpenGLRenderWindow` ([src/app/ReaderMainWindow.cpp:767](../src/app/ReaderMainWindow.cpp:767), [src/app/ReaderMainWindow.cpp:771](../src/app/ReaderMainWindow.cpp:771)), and `MoleculeScene` uses a main renderer plus an overlay renderer that share a camera ([src/app/MoleculeScene.cpp:102](../src/app/MoleculeScene.cpp:102), [src/app/MoleculeScene.cpp:124](../src/app/MoleculeScene.cpp:124)). The VTK class references for `QVTKOpenGLNativeWidget`, `vtkGenericOpenGLRenderWindow`, `vtkRenderer`, and `vtkRenderWindow` support this render-window/renderer setup.

### Temporal Strip Panels

What draws: `TemporalStripPanel::paint` draws the header, grid, time ticks, selected range, y labels, trace, border, and cursor ([src/app/StripStackWidget.cpp:86](../src/app/StripStackWidget.cpp:86), [src/app/StripStackWidget.cpp:119](../src/app/StripStackWidget.cpp:119)). It draws direct or decimated traces ([src/app/StripStackWidget.cpp:128](../src/app/StripStackWidget.cpp:128), [src/app/StripStackWidget.cpp:280](../src/app/StripStackWidget.cpp:280)).

What reaches it: all active `strip.*` modes are considered visibly renderable by the capability table ([src/model/DisplayModeCapability.h:25](../src/model/DisplayModeCapability.h:25), [src/model/DisplayModeCapability.h:28](../src/model/DisplayModeCapability.h:28)), and the controller routes any signal with a strip mode through `buildGenericTracks` ([src/app/DashboardDisplayController.cpp:1695](../src/app/DashboardDisplayController.cpp:1695), [src/app/DashboardDisplayController.cpp:1699](../src/app/DashboardDisplayController.cpp:1699)). The dock converts them to temporal tracks ([src/app/DashboardStripDock.cpp:297](../src/app/DashboardStripDock.cpp:297), [src/app/DashboardStripDock.cpp:312](../src/app/DashboardStripDock.cpp:312)), and `setTracks` creates `TemporalStripPanel` instances ([src/app/StripStackWidget.cpp:401](../src/app/StripStackWidget.cpp:401), [src/app/StripStackWidget.cpp:418](../src/app/StripStackWidget.cpp:418)).

Strip modes that can therefore draw pixels when sampling succeeds:

- Tensor strips: `strip.tensor.T0`, `strip.tensor.T1`, `strip.tensor.T2`, `strip.tensor.component` are advertised in tensor helpers ([src/model/TrajectorySignalCatalog.cpp:84](../src/model/TrajectorySignalCatalog.cpp:84), [src/model/TrajectorySignalCatalog.cpp:91](../src/model/TrajectorySignalCatalog.cpp:91)). T0/T1/T2 sampling exists for tensor rows ([src/app/DashboardDisplayController.cpp:292](../src/app/DashboardDisplayController.cpp:292), [src/app/DashboardDisplayController.cpp:308](../src/app/DashboardDisplayController.cpp:308)). `strip.tensor.component` is not a real component selector; for tensor sampling it returns one component path, and for T2 vectors it falls back to the first value ([src/app/DashboardDisplayController.cpp:305](../src/app/DashboardDisplayController.cpp:305), [src/app/DashboardDisplayController.cpp:314](../src/app/DashboardDisplayController.cpp:314)).
- Vector strips: `strip.vector.component` and `strip.vector.magnitude` are advertised ([src/model/TrajectorySignalCatalog.cpp:135](../src/model/TrajectorySignalCatalog.cpp:135), [src/model/TrajectorySignalCatalog.cpp:140](../src/model/TrajectorySignalCatalog.cpp:140)); vector sampling covers x/y/z/magnitude ([src/app/DashboardDisplayController.cpp:319](../src/app/DashboardDisplayController.cpp:319), [src/app/DashboardDisplayController.cpp:331](../src/app/DashboardDisplayController.cpp:331)).
- Scalar/count/category/per-class/rollup/event strips: helpers advertise them ([src/model/TrajectorySignalCatalog.cpp:120](../src/model/TrajectorySignalCatalog.cpp:120), [src/model/TrajectorySignalCatalog.cpp:204](../src/model/TrajectorySignalCatalog.cpp:204)); generic NPY sampling uses named x/y/z/magnitude/T0/T1/T2 fields or the first row value ([src/app/DashboardDisplayController.cpp:242](../src/app/DashboardDisplayController.cpp:242), [src/app/DashboardDisplayController.cpp:273](../src/app/DashboardDisplayController.cpp:273)).
- `strip.spectrum` is advertised in many helper lists ([src/model/TrajectorySignalCatalog.cpp:90](../src/model/TrajectorySignalCatalog.cpp:90), [src/model/TrajectorySignalCatalog.cpp:124](../src/model/TrajectorySignalCatalog.cpp:124), [src/model/TrajectorySignalCatalog.cpp:140](../src/model/TrajectorySignalCatalog.cpp:140), [src/model/TrajectorySignalCatalog.cpp:170](../src/model/TrajectorySignalCatalog.cpp:170)), but it goes through the same temporal-strip path unless removed or implemented as a real spectrum.

### Spectrum Strip Panel

There is a `SpectrumStripPanel` class and it paints a spectrum trace ([src/app/StripStackWidget.cpp:288](../src/app/StripStackWidget.cpp:288), [src/app/StripStackWidget.cpp:319](../src/app/StripStackWidget.cpp:319), [src/app/StripStackWidget.cpp:359](../src/app/StripStackWidget.cpp:359), [src/app/StripStackWidget.cpp:382](../src/app/StripStackWidget.cpp:382)). `StripStackWidget::setSpectrumTracks` can instantiate it ([src/app/StripStackWidget.cpp:420](../src/app/StripStackWidget.cpp:420), [src/app/StripStackWidget.cpp:439](../src/app/StripStackWidget.cpp:439)).

But the normal dock path currently clears spectrum tracks every refresh ([src/app/DashboardStripDock.cpp:313](../src/app/DashboardStripDock.cpp:313), [src/app/DashboardStripDock.cpp:315](../src/app/DashboardStripDock.cpp:315)). Therefore `strip.spectrum` is advertised but does not currently render as a spectrum strip through the dashboard.

### Static Dashboard Panels

The static panels that actually build are hard-coded in `DashboardDisplayController::rebuild`. The common gate is `DisplayModeCapabilityFor(modeId).buildsPanelWidget` ([src/app/DashboardDisplayController.cpp:1605](../src/app/DashboardDisplayController.cpp:1605), [src/app/DashboardDisplayController.cpp:1607](../src/app/DashboardDisplayController.cpp:1607)).

Confirmed rendering paths:

- `static.bar.sequence` -> `SequenceBarPanel`. Branches exist for iRED order parameters, reorientational dynamics, and dihedral autocorrelation ([src/app/DashboardDisplayController.cpp:1621](../src/app/DashboardDisplayController.cpp:1621), [src/app/DashboardDisplayController.cpp:1625](../src/app/DashboardDisplayController.cpp:1625), [src/app/DashboardDisplayController.cpp:1634](../src/app/DashboardDisplayController.cpp:1634), [src/app/DashboardDisplayController.cpp:1643](../src/app/DashboardDisplayController.cpp:1643), [src/app/DashboardDisplayController.cpp:1648](../src/app/DashboardDisplayController.cpp:1648), [src/app/DashboardDisplayController.cpp:1654](../src/app/DashboardDisplayController.cpp:1654)). `SequenceBarPanel` paints per-row bars and overlays ([src/app/SequenceBarPanel.cpp:135](../src/app/SequenceBarPanel.cpp:135), [src/app/SequenceBarPanel.cpp:326](../src/app/SequenceBarPanel.cpp:326)).
- `static.curve.lag.animated` -> `LagDecayPanel`. Branches exist for kernel ACF, reorientational dynamics, and dihedral autocorrelation ([src/app/DashboardDisplayController.cpp:1630](../src/app/DashboardDisplayController.cpp:1630), [src/app/DashboardDisplayController.cpp:1633](../src/app/DashboardDisplayController.cpp:1633), [src/app/DashboardDisplayController.cpp:1644](../src/app/DashboardDisplayController.cpp:1644), [src/app/DashboardDisplayController.cpp:1647](../src/app/DashboardDisplayController.cpp:1647), [src/app/DashboardDisplayController.cpp:1655](../src/app/DashboardDisplayController.cpp:1655), [src/app/DashboardDisplayController.cpp:1661](../src/app/DashboardDisplayController.cpp:1661)). `LagDecayPanel` paints lag curves and an animated cursor ([src/app/LagDecayPanel.cpp:54](../src/app/LagDecayPanel.cpp:54), [src/app/LagDecayPanel.cpp:150](../src/app/LagDecayPanel.cpp:150)).
- `static.chord.coupling` -> `ChordCouplingPanel` for kernel coherence ([src/app/DashboardDisplayController.cpp:1662](../src/app/DashboardDisplayController.cpp:1662), [src/app/DashboardDisplayController.cpp:1665](../src/app/DashboardDisplayController.cpp:1665)). The panel draws nodes and coupling arcs ([src/app/ChordCouplingPanel.cpp:31](../src/app/ChordCouplingPanel.cpp:31), [src/app/ChordCouplingPanel.cpp:116](../src/app/ChordCouplingPanel.cpp:116)).
- `static.fixed_freq` -> `FixedFreqPanel` for reorientational dynamics ([src/app/DashboardDisplayController.cpp:1666](../src/app/DashboardDisplayController.cpp:1666), [src/app/DashboardDisplayController.cpp:1670](../src/app/DashboardDisplayController.cpp:1670)). The panel filters finite samples and draws a log-frequency curve with markers ([src/app/FixedFreqPanel.cpp:54](../src/app/FixedFreqPanel.cpp:54), [src/app/FixedFreqPanel.cpp:179](../src/app/FixedFreqPanel.cpp:179)).
- `static.spectrum.power` -> `PowerSpectrumPanel` for kernel dynamics PSD ([src/app/DashboardDisplayController.cpp:1626](../src/app/DashboardDisplayController.cpp:1626), [src/app/DashboardDisplayController.cpp:1629](../src/app/DashboardDisplayController.cpp:1629)); `PowerSpectrumPanel` draws PSD polylines and a legend ([src/app/PowerSpectrumPanel.cpp:49](../src/app/PowerSpectrumPanel.cpp:49), [src/app/PowerSpectrumPanel.cpp:139](../src/app/PowerSpectrumPanel.cpp:139)). However, the capability table marks `static.spectrum.power` as not having a visible surface while still building a panel/ref ([src/model/DisplayModeCapability.h:31](../src/model/DisplayModeCapability.h:31), [src/model/DisplayModeCapability.h:33](../src/model/DisplayModeCapability.h:33)), so the normal dialog gate does not make it user-selectable.

### Static Tensor Glyph

This is the important half-built case. `static.tensor` is in the capability table, but it is invisible, does not build a panel, and only advertises a panel ref trigger ([src/model/DisplayModeCapability.h:37](../src/model/DisplayModeCapability.h:37), [src/model/DisplayModeCapability.h:40](../src/model/DisplayModeCapability.h:40)). The controller explicitly omits the trigger and says there is no UI gesture yet ([src/app/DashboardDisplayController.cpp:1675](../src/app/DashboardDisplayController.cpp:1675), [src/app/DashboardDisplayController.cpp:1693](../src/app/DashboardDisplayController.cpp:1693), [src/app/DashboardDisplayController.cpp:1748](../src/app/DashboardDisplayController.cpp:1748), [src/app/DashboardDisplayController.cpp:1751](../src/app/DashboardDisplayController.cpp:1751)).

The scene-side glyph machinery exists. `SceneRevealOverlay` exposes `revealTensor` ([src/app/SceneRevealOverlay.h:51](../src/app/SceneRevealOverlay.h:51), [src/app/SceneRevealOverlay.h:66](../src/app/SceneRevealOverlay.h:66)), builds a tensor sphere/transform/filter/mapper/actor pipeline ([src/app/SceneRevealOverlay.cpp:138](../src/app/SceneRevealOverlay.cpp:138), [src/app/SceneRevealOverlay.cpp:169](../src/app/SceneRevealOverlay.cpp:169)), and decomposes a tensor into transform/visibility in `applyTensorFrame` ([src/app/SceneRevealOverlay.cpp:401](../src/app/SceneRevealOverlay.cpp:401), [src/app/SceneRevealOverlay.cpp:445](../src/app/SceneRevealOverlay.cpp:445)). This aligns with VTK's actor/mapper and transform-filter pipeline model, but the dashboard does not currently invoke it for `static.tensor`.

There are stale comments that imply this should be wired: the display controller header says `static.tensor` should reveal a tensor glyph if an overlay exists ([src/app/DashboardDisplayController.h:117](../src/app/DashboardDisplayController.h:117), [src/app/DashboardDisplayController.h:122](../src/app/DashboardDisplayController.h:122)), and `ReaderMainWindow` comments say `static.tensor` fires an ellipsoid glyph ([src/app/ReaderMainWindow.cpp:418](../src/app/ReaderMainWindow.cpp:418), [src/app/ReaderMainWindow.cpp:421](../src/app/ReaderMainWindow.cpp:421)). The implementation contradicts those comments.

### Scene Overlays

These are real, separate from the dashboard display-mode system:

- Ribbon overlay: toolbar action exists and is checked by default ([src/app/ReaderMainWindow.cpp:918](../src/app/ReaderMainWindow.cpp:918), [src/app/ReaderMainWindow.cpp:923](../src/app/ReaderMainWindow.cpp:923)). `QtBackboneRibbonOverlay::Build` creates VTK points/arrays, a `vtkProteinRibbonFilter`, mapper, and actor, then adds it to the renderer ([src/app/QtBackboneRibbonOverlay.cpp:172](../src/app/QtBackboneRibbonOverlay.cpp:172), [src/app/QtBackboneRibbonOverlay.cpp:272](../src/app/QtBackboneRibbonOverlay.cpp:272)). It updates on frame changes ([src/app/QtBackboneRibbonOverlay.cpp:392](../src/app/QtBackboneRibbonOverlay.cpp:392), [src/app/QtBackboneRibbonOverlay.cpp:400](../src/app/QtBackboneRibbonOverlay.cpp:400)). VTK's `vtkProteinRibbonFilter` docs cover the filter role in that pipeline.
- Ring overlay: toolbar action exists and is checked by default ([src/app/ReaderMainWindow.cpp:924](../src/app/ReaderMainWindow.cpp:924), [src/app/ReaderMainWindow.cpp:928](../src/app/ReaderMainWindow.cpp:928)). `QtRingPolygonOverlay` builds polygon and arrow actors and adds them to the renderer ([src/app/QtRingPolygonOverlay.cpp:98](../src/app/QtRingPolygonOverlay.cpp:98), [src/app/QtRingPolygonOverlay.cpp:156](../src/app/QtRingPolygonOverlay.cpp:156)).
- Butterfly/field-grid overlay: toolbar action exists, off by default and labelled expensive ([src/app/ReaderMainWindow.cpp:930](../src/app/ReaderMainWindow.cpp:930), [src/app/ReaderMainWindow.cpp:935](../src/app/ReaderMainWindow.cpp:935)). The overlay builds contour actors from image data and `vtkContourFilter` ([src/app/QtFieldGridOverlay.h:5](../src/app/QtFieldGridOverlay.h:5), [src/app/QtFieldGridOverlay.cpp:52](../src/app/QtFieldGridOverlay.cpp:52), [src/app/QtFieldGridOverlay.cpp:125](../src/app/QtFieldGridOverlay.cpp:125)). VTK's `vtkContourFilter` docs cover this contour pipeline.
- B-field stream overlay: toolbar action exists, off by default and labelled expensive ([src/app/ReaderMainWindow.cpp:937](../src/app/ReaderMainWindow.cpp:937), [src/app/ReaderMainWindow.cpp:942](../src/app/ReaderMainWindow.cpp:942)). The overlay builds a structured grid, stream tracer, tube filter, mapper, and actor ([src/app/QtBFieldStreamOverlay.h:7](../src/app/QtBFieldStreamOverlay.h:7), [src/app/QtBFieldStreamOverlay.cpp:73](../src/app/QtBFieldStreamOverlay.cpp:73), [src/app/QtBFieldStreamOverlay.cpp:157](../src/app/QtBFieldStreamOverlay.cpp:157)). VTK's `vtkStreamTracer` docs cover the streamline part.
- Measurement overlay: `MeasurementOverlay` is the active successor to the dormant `QtSelectionOverlay` ([src/app/MeasurementOverlay.h:13](../src/app/MeasurementOverlay.h:13), [src/app/MeasurementOverlay.h:14](../src/app/MeasurementOverlay.h:14)). It builds sphere and line actors ([src/app/MeasurementOverlay.cpp:82](../src/app/MeasurementOverlay.cpp:82), [src/app/MeasurementOverlay.cpp:141](../src/app/MeasurementOverlay.cpp:141)) and updates selected atoms/lines on selection and frame changes ([src/app/MeasurementOverlay.cpp:222](../src/app/MeasurementOverlay.cpp:222), [src/app/MeasurementOverlay.cpp:288](../src/app/MeasurementOverlay.cpp:288)). `ReaderMainWindow` wires selection changes to it ([src/app/ReaderMainWindow.cpp:302](../src/app/ReaderMainWindow.cpp:302), [src/app/ReaderMainWindow.cpp:306](../src/app/ReaderMainWindow.cpp:306)).
- Scene reveal overlay: dashboard reveal requests call `MoleculeScene::revealBinding` ([src/app/ReaderMainWindow.cpp:416](../src/app/ReaderMainWindow.cpp:416), [src/app/ReaderMainWindow.cpp:417](../src/app/ReaderMainWindow.cpp:417)), and `MoleculeScene` can focus the camera on revealed atoms, residues, bonds, tuples, or ring/bond-vector anchors ([src/app/MoleculeScene.cpp:549](../src/app/MoleculeScene.cpp:549), [src/app/MoleculeScene.cpp:584](../src/app/MoleculeScene.cpp:584)).

## 3. Clean Class-System Proposal

Goal: keep the screen small and desk-ready, but make hollow modes structurally un-offerable. A visualisation should only appear in UI if there is a typed class behind it that can say "I support this descriptor and I can build or update the surface."

### Minimal Shape

Introduce a small registered hierarchy, not a broad framework:

```cpp
enum class VisualizationType {
    TemporalStrip,
    TensorGlyph,
    AtomColor,
    SequenceBar,
    LagCurve,
    ChordCoupling,
    FixedFrequency
};

class VisualizationDefinition {
public:
    virtual ~VisualizationDefinition() = default;
    virtual VisualizationType type() const = 0;
    virtual QString label() const = 0;
    virtual bool supports(const SignalDescriptor& descriptor) const = 0;
    virtual bool isAvailable(const VisualizationContext& context,
                             const SignalDescriptor& descriptor) const = 0;
    virtual DisplaySurface surface() const = 0;
};
```

Then split construction by surface rather than mode string:

- `StripVisualizationDefinition` returns typed `StripSeriesRequest` objects, with tensor/vector component policy as enums instead of `strip.tensor.T2` text.
- `PanelVisualizationDefinition` returns `std::unique_ptr<AbstractStripPanel>`.
- `SceneVisualizationDefinition` updates or reveals through a typed overlay object.

The existing `AbstractStripPanel` hierarchy can stay. It is useful for drawing and reveal behaviour ([src/app/AbstractStripPanel.h:145](../src/app/AbstractStripPanel.h:145), [src/app/AbstractStripPanel.h:180](../src/app/AbstractStripPanel.h:180)). The change is above it: the selector and controller should not decide visualisation identity with string prefixes.

### Registry As The Offerability Gate

Replace or wrap `DisplayModeCapabilityFor` with a registry of concrete definitions. Today the table marks every `strip.` string visible by prefix and separately lists a few static strings ([src/model/DisplayModeCapability.h:25](../src/model/DisplayModeCapability.h:25), [src/model/DisplayModeCapability.h:41](../src/model/DisplayModeCapability.h:41)). That is what allows catalog/UI drift. In the proposed model:

- `TrajectorySignalCatalog` describes data shape, storage, anchors, temporal/static nature, tensor/vector/scalar category. It does not hand-author UI mode strings.
- `VisualizationRegistry::supportedTypes(descriptor)` asks registered definitions whether they support the descriptor.
- `SignalDisplayDialog` builds filters and checkboxes from registered definitions, replacing `allModeKinds`, `canonicalModeId`, and `modeMatchesKind` ([src/app/SignalDisplayDialog.cpp:156](../src/app/SignalDisplayDialog.cpp:156), [src/app/SignalDisplayDialog.cpp:207](../src/app/SignalDisplayDialog.cpp:207), [src/app/SignalDisplayDialog.cpp:360](../src/app/SignalDisplayDialog.cpp:360)).
- The selection model stores typed selections: `VisualizationSelection { VisualizationType type; QVariant options; }`, not `QString displayModeId` ([src/model/DashboardSignal.h:266](../src/model/DashboardSignal.h:266), [src/model/DashboardSignal.h:275](../src/model/DashboardSignal.h:275)).
- `DashboardPanelModel::DisplayRefsForSignal` becomes a registry query instead of mode/channel string logic ([src/model/DashboardPanelModel.cpp:89](../src/model/DashboardPanelModel.cpp:89), [src/model/DashboardPanelModel.cpp:123](../src/model/DashboardPanelModel.cpp:123)).
- `DashboardDisplayController::rebuild` dispatches to definitions by typed surface. The current hard-coded branches from `static.bar.sequence` through `static.fixed_freq` become definitions registered only if the lead keeps them ([src/app/DashboardDisplayController.cpp:1621](../src/app/DashboardDisplayController.cpp:1621), [src/app/DashboardDisplayController.cpp:1670](../src/app/DashboardDisplayController.cpp:1670)).

This can be introduced incrementally. A first pass can keep string ids internally but require a startup validation rule: every catalog mode id must be registered, and every registered mode must have either a strip sampler, a panel builder, or a scene overlay handler. That one rule would have caught `static.table`, `static.atomColor`, and `strip.spectrum` immediately.

### Suggested Definitions For The Basic Screen

Keep the initial registry small:

- `TemporalStripVisualization`: owns the generic strip path and typed component policy. It covers scalar, vector magnitude/component, tensor T0/T2, count, and a deliberately chosen subset of category/event if the lead wants those plotted numerically.
- `TensorGlyphVisualization`: owns T2/rank-2 tensor reveal and scene overlay behaviour. It can reuse `SceneRevealOverlay::revealTensor` initially ([src/app/SceneRevealOverlay.h:51](../src/app/SceneRevealOverlay.h:51), [src/app/SceneRevealOverlay.cpp:376](../src/app/SceneRevealOverlay.cpp:376)).
- `AtomColorVisualization`: owns scalar-to-atom/residue colouring if chosen. This replaces hollow `static.atomColor`.
- `SequenceBarVisualization`: owns sequence/residue summary bars if chosen. It reuses `SequenceBarPanel`, but the storage-path branching moves behind the definition.

Optional definitions can be left out of the first registry: lag curves, chord coupling, fixed-frequency plots, and power spectra. If they are not registered, they cannot appear in the selector.

## 4. Must-Have Proposal For Lead Decision

The reader's job is vetting H5 fields before they become thesis claims. The must-have set should be visualisations that expose field correctness, spatial locality, and tensor structure. Below is a decision set, not a decision on behalf of the lead.

### Have: Strips

What it shows: temporal strip charts for sampled metrics at the current selection, including scalar values, vector magnitude/components, and tensor-derived T0/T2 values. The controller already samples tensor T0/T1/T2 and vector magnitude/components ([src/app/DashboardDisplayController.cpp:292](../src/app/DashboardDisplayController.cpp:292), [src/app/DashboardDisplayController.cpp:331](../src/app/DashboardDisplayController.cpp:331)), and the stack already paints temporal traces ([src/app/StripStackWidget.cpp:86](../src/app/StripStackWidget.cpp:86), [src/app/StripStackWidget.cpp:119](../src/app/StripStackWidget.cpp:119)).

Why must-have for vetting: strips are the fastest way to catch bad frames, discontinuities, missing channels, impossible magnitudes, and temporal disagreement between related fields before a number appears in a thesis table.

Data needed: existing dense H5/NPY/ORCA sampling plans and `ChannelBuffer` track path. `buildGenericTracks` and `extendToFrame` already own most of this path ([src/app/DashboardDisplayController.cpp:2368](../src/app/DashboardDisplayController.cpp:2368), [src/app/DashboardDisplayController.cpp:2414](../src/app/DashboardDisplayController.cpp:2414), [src/app/DashboardDisplayController.cpp:2539](../src/app/DashboardDisplayController.cpp:2539), [src/app/DashboardDisplayController.cpp:2575](../src/app/DashboardDisplayController.cpp:2575)).

Rough build cost: low. Keep and tighten. Remove or fix fake `strip.spectrum`; decide whether `strip.tensor.component` deserves a component selector or should be cut; keep T2 visible because tensor anisotropy is central.

### Candidate 1: Tensor Glyph Overlay

What it shows: a rank-2 tensor glyph at the relevant atom, residue, or bond-vector anchor, ideally preserving tensor shape/orientation rather than reducing the field to a scalar strip. This is the visual home for T2/tensor structure.

Why must-have for vetting: the rank-2 tensor is the thesis-bearing object. A T2 strip can say "large anisotropy here"; a glyph can show whether the tensor orientation and shape are physically plausible in the molecular frame.

Data needed: a full tensor or a defensible reconstruction/eigendecomposition policy for each supported descriptor. The scene already has a tensor-glyph pipeline and tensor math hook ([src/app/SceneRevealOverlay.cpp:138](../src/app/SceneRevealOverlay.cpp:138), [src/app/SceneRevealOverlay.cpp:401](../src/app/SceneRevealOverlay.cpp:401)), but `static.tensor` is not wired from the dashboard ([src/app/DashboardDisplayController.cpp:1675](../src/app/DashboardDisplayController.cpp:1675)).

Rough build cost: medium if limited to the reorient/tensor cases already close to `SceneRevealOverlay::revealTensor`; medium-high if generalized across all shielding/EFG tensor descriptors with UI selection, legend, scale controls, and missing-data handling.

### Candidate 2: Atom Or Residue Colour Map

What it shows: colour the molecular atoms/residues by the selected metric at the current frame or selected summary value. It should support scalar values and tensor-derived T2 magnitude at minimum, with a legend and explicit missing-data colour.

Why must-have for vetting: spatial outliers are hard to see in a strip stack. Per-metric colouring lets the reviewer see whether a suspicious value is localised to an atom/residue/structural region or is a global artefact.

Data needed: per-atom or per-residue values at the current frame, anchor mapping from descriptors, and an overlay or mapper path that can update colours on frame changes. The catalog already advertises `static.atomColor` widely ([src/model/TrajectorySignalCatalog.cpp:94](../src/model/TrajectorySignalCatalog.cpp:94), [src/model/TrajectorySignalCatalog.cpp:127](../src/model/TrajectorySignalCatalog.cpp:127), [src/model/TrajectorySignalCatalog.cpp:173](../src/model/TrajectorySignalCatalog.cpp:173)), but there is no corresponding capability/renderer row ([src/model/DisplayModeCapability.h:29](../src/model/DisplayModeCapability.h:29)).

Rough build cost: medium. The sampling and anchor logic already exists in pieces, but the VTK/Qt colour surface, legend, and update contract are new work.

### Candidate 3: Sequence Or Structural Summary Bar

What it shows: compact bars across residue index, bond-vector index, or another ordered structural axis. Existing examples include S2/correlation time/reorient metrics and dihedral summaries.

Why must-have for vetting: many NMR and structural questions are "where along the sequence or bond-vector set is this field wrong?" A sequence/structural bar is denser than many strips and supports reveal-to-scene interaction.

Data needed: static scalar rows with anchor ids and labels. The existing `SequenceBarPanel` already paints this and returns reveal bindings on click ([src/app/SequenceBarPanel.cpp:124](../src/app/SequenceBarPanel.cpp:124), [src/app/SequenceBarPanel.cpp:135](../src/app/SequenceBarPanel.cpp:135)). Current controller support is storage-path specific ([src/app/DashboardDisplayController.cpp:1771](../src/app/DashboardDisplayController.cpp:1771), [src/app/DashboardDisplayController.cpp:1891](../src/app/DashboardDisplayController.cpp:1891), [src/app/DashboardDisplayController.cpp:2154](../src/app/DashboardDisplayController.cpp:2154)).

Rough build cost: low-medium. Reuse the existing panel, but move the storage-path knowledge into one registered `SequenceBarVisualization` and cut unrelated static placeholders.

### Optional Instead Of Candidate 2 Or 3: Advanced Dynamics Panels

Power spectrum, lag decay, chord coupling, and fixed-frequency panels are real or nearly real, but they are specialised dynamics diagnostics rather than general field-vetting surfaces. If the lead keeps them, they should be explicit advanced visualisations, not incidental static modes. Their concrete panel branches are already present ([src/app/DashboardDisplayController.cpp:1626](../src/app/DashboardDisplayController.cpp:1626)).

### Minimalist Option: Strips Only

The lead can choose strips only for the first clean pass. In that case, cut all metric-mode static options from the selector and keep only structural scene overlays that are independently useful. This produces the smallest honest screen, but it leaves the rank-2 tensor visually under-exposed unless T2 strips are considered enough.

## 5. Remove List

This list is mechanical: mode or option, why it is hollow/dead/misleading, and where it is offered.

### `strip.spectrum`

Why cut or implement: it is advertised as a spectrum, but normal dashboard refresh clears spectrum tracks ([src/app/DashboardStripDock.cpp:313](../src/app/DashboardStripDock.cpp:313)). The existing `SpectrumStripPanel` is not reached through normal metric selection ([src/app/StripStackWidget.cpp:420](../src/app/StripStackWidget.cpp:420)).

Where it is spattered:

- Catalog helper lists: tensor, EFG, scalar, vector, and per-class modes include `strip.spectrum` ([src/model/TrajectorySignalCatalog.cpp:84](../src/model/TrajectorySignalCatalog.cpp:84), [src/model/TrajectorySignalCatalog.cpp:103](../src/model/TrajectorySignalCatalog.cpp:103), [src/model/TrajectorySignalCatalog.cpp:120](../src/model/TrajectorySignalCatalog.cpp:120), [src/model/TrajectorySignalCatalog.cpp:135](../src/model/TrajectorySignalCatalog.cpp:135), [src/model/TrajectorySignalCatalog.cpp:165](../src/model/TrajectorySignalCatalog.cpp:165)).
- ORCA total shielding includes `strip.spectrum` ([src/model/TrajectorySignalCatalog.cpp:1078](../src/model/TrajectorySignalCatalog.cpp:1078)).
- Dialog kind/label/canonical id expose "Spectrum" and map it to `strip.spectrum` ([src/app/SignalDisplayDialog.cpp:56](../src/app/SignalDisplayDialog.cpp:56), [src/app/SignalDisplayDialog.cpp:108](../src/app/SignalDisplayDialog.cpp:108), [src/app/SignalDisplayDialog.cpp:156](../src/app/SignalDisplayDialog.cpp:156)).
- Dialog mode filter, descriptor summary/tooltip, add checkboxes, and active checkboxes all surface it when present ([src/app/SignalDisplayDialog.cpp:437](../src/app/SignalDisplayDialog.cpp:437), [src/app/SignalDisplayDialog.cpp:801](../src/app/SignalDisplayDialog.cpp:801), [src/app/SignalDisplayDialog.cpp:843](../src/app/SignalDisplayDialog.cpp:843), [src/app/SignalDisplayDialog.cpp:874](../src/app/SignalDisplayDialog.cpp:874)).

### `static.table`

Why cut or implement: there is no dashboard table renderer behind the mode. Unknown/static unlisted modes return no capability ([src/model/DisplayModeCapability.h:47](../src/model/DisplayModeCapability.h:47)), and the controller has no table branch in its static panel dispatch ([src/app/DashboardDisplayController.cpp:1621](../src/app/DashboardDisplayController.cpp:1621)).

Where it is spattered:

- Static helper lists across tensor, EFG, scalar, vector, category, per-class, rollup, and event descriptors include `static.table` ([src/model/TrajectorySignalCatalog.cpp:94](../src/model/TrajectorySignalCatalog.cpp:94)).
- Many descriptor groups add it explicitly, including kernel coherence, dihedral, reorient, iRED, topology, geometry, and selection-style descriptors ([src/model/TrajectorySignalCatalog.cpp:613](../src/model/TrajectorySignalCatalog.cpp:613), [src/model/TrajectorySignalCatalog.cpp:636](../src/model/TrajectorySignalCatalog.cpp:636), [src/model/TrajectorySignalCatalog.cpp:757](../src/model/TrajectorySignalCatalog.cpp:757), [src/model/TrajectorySignalCatalog.cpp:1104](../src/model/TrajectorySignalCatalog.cpp:1104)).
- Dialog exposes a Table kind and maps it to `static.table` ([src/app/SignalDisplayDialog.cpp:56](../src/app/SignalDisplayDialog.cpp:56), [src/app/SignalDisplayDialog.cpp:108](../src/app/SignalDisplayDialog.cpp:108), [src/app/SignalDisplayDialog.cpp:156](../src/app/SignalDisplayDialog.cpp:156)).
- Dialog checkboxes and filters include the Table option even though the capability gate disables it ([src/app/SignalDisplayDialog.cpp:360](../src/app/SignalDisplayDialog.cpp:360), [src/app/SignalDisplayDialog.cpp:801](../src/app/SignalDisplayDialog.cpp:801), [src/app/SignalDisplayDialog.cpp:843](../src/app/SignalDisplayDialog.cpp:843)).

### `static.atomColor` / Color Map

Why cut or implement: this is a good candidate must-have, but today it has no renderer or capability row. It is therefore advertised-dead.

Where it is spattered:

- Catalog helper lists advertise it for tensor, EFG, scalar, and per-class static modes ([src/model/TrajectorySignalCatalog.cpp:94](../src/model/TrajectorySignalCatalog.cpp:94), [src/model/TrajectorySignalCatalog.cpp:111](../src/model/TrajectorySignalCatalog.cpp:111), [src/model/TrajectorySignalCatalog.cpp:127](../src/model/TrajectorySignalCatalog.cpp:127), [src/model/TrajectorySignalCatalog.cpp:173](../src/model/TrajectorySignalCatalog.cpp:173)).
- ORCA total shielding advertises it ([src/model/TrajectorySignalCatalog.cpp:1078](../src/model/TrajectorySignalCatalog.cpp:1078)).
- Dialog exposes "Color map" and maps it to `static.atomColor` ([src/app/SignalDisplayDialog.cpp:56](../src/app/SignalDisplayDialog.cpp:56), [src/app/SignalDisplayDialog.cpp:108](../src/app/SignalDisplayDialog.cpp:108), [src/app/SignalDisplayDialog.cpp:156](../src/app/SignalDisplayDialog.cpp:156)).
- Dialog checkboxes/filter/summary surface it like other modes ([src/app/SignalDisplayDialog.cpp:437](../src/app/SignalDisplayDialog.cpp:437), [src/app/SignalDisplayDialog.cpp:801](../src/app/SignalDisplayDialog.cpp:801), [src/app/SignalDisplayDialog.cpp:843](../src/app/SignalDisplayDialog.cpp:843), [src/app/SignalDisplayDialog.cpp:874](../src/app/SignalDisplayDialog.cpp:874)).

### `static.tensor` / Glyph Overlay

Why cut or implement: tensor glyph code exists in `SceneRevealOverlay`, but the metric-mode trigger is not wired. Keeping it advertised without a class/gesture is misleading, especially because T2/tensor is central.

Where it is spattered:

- Catalog helper lists and descriptors advertise `static.tensor` for tensor/EFG/reorient/ORCA fields ([src/model/TrajectorySignalCatalog.cpp:94](../src/model/TrajectorySignalCatalog.cpp:94), [src/model/TrajectorySignalCatalog.cpp:111](../src/model/TrajectorySignalCatalog.cpp:111), [src/model/TrajectorySignalCatalog.cpp:815](../src/model/TrajectorySignalCatalog.cpp:815), [src/model/TrajectorySignalCatalog.cpp:1078](../src/model/TrajectorySignalCatalog.cpp:1078)).
- Capability table includes `static.tensor`, but marks it invisible and no-panel ([src/model/DisplayModeCapability.h:37](../src/model/DisplayModeCapability.h:37)).
- Dialog exposes "Glyph" and maps it to `static.tensor` ([src/app/SignalDisplayDialog.cpp:56](../src/app/SignalDisplayDialog.cpp:56), [src/app/SignalDisplayDialog.cpp:108](../src/app/SignalDisplayDialog.cpp:108), [src/app/SignalDisplayDialog.cpp:156](../src/app/SignalDisplayDialog.cpp:156)).
- Controller comments explicitly say the trigger is omitted ([src/app/DashboardDisplayController.cpp:1675](../src/app/DashboardDisplayController.cpp:1675), [src/app/DashboardDisplayController.cpp:1748](../src/app/DashboardDisplayController.cpp:1748)).
- Stale comments in the controller header and main window still imply this should work ([src/app/DashboardDisplayController.h:117](../src/app/DashboardDisplayController.h:117), [src/app/ReaderMainWindow.cpp:418](../src/app/ReaderMainWindow.cpp:418)).

### Other Hollow Static Modes

Why cut or map to real definitions: these are catalog/static mode strings with no capability rows and no controller panel or overlay branch. They should disappear unless one is deliberately promoted into the typed registry.

Modes:

- `static.scalar`
- `static.efg`
- `static.vector`
- `static.vectorGlyph`
- `static.category`
- `static.per-class`
- `static.rollup`
- `static.event`
- `static.system`
- `static.embedding`
- `static.geometry`
- `static.topology`

Where they are spattered:

- Static helper lists advertise many of them ([src/model/TrajectorySignalCatalog.cpp:111](../src/model/TrajectorySignalCatalog.cpp:111)).
- Topology, geometry, frame NPY, embedding/system, EFG, vector, category, event, and ORCA descriptor blocks include these static ids ([src/model/TrajectorySignalCatalog.cpp:900](../src/model/TrajectorySignalCatalog.cpp:900), [src/model/TrajectorySignalCatalog.cpp:985](../src/model/TrajectorySignalCatalog.cpp:985), [src/model/TrajectorySignalCatalog.cpp:1078](../src/model/TrajectorySignalCatalog.cpp:1078)).
- Dialog summary/tooltip surfaces raw mode ids even when no renderer exists ([src/app/SignalDisplayDialog.cpp:437](../src/app/SignalDisplayDialog.cpp:437)).
- Dialog mode filtering maps broad UI kinds onto string groups, so hollow modes still affect discovery even when disabled ([src/app/SignalDisplayDialog.cpp:613](../src/app/SignalDisplayDialog.cpp:613)).

### Misleading Or Under-Specified Strip Modes

These draw pixels through generic temporal strips, but their names imply semantics the UI does not provide:

- `strip.tensor.component`: advertised in tensor and EFG helper lists ([src/model/TrajectorySignalCatalog.cpp:84](../src/model/TrajectorySignalCatalog.cpp:84), [src/model/TrajectorySignalCatalog.cpp:103](../src/model/TrajectorySignalCatalog.cpp:103)); sampling takes a fixed/first component path without a component chooser ([src/app/DashboardDisplayController.cpp:305](../src/app/DashboardDisplayController.cpp:305), [src/app/DashboardDisplayController.cpp:314](../src/app/DashboardDisplayController.cpp:314)).
- `strip.embedding.component`: advertised for embeddings ([src/model/TrajectorySignalCatalog.cpp:442](../src/model/TrajectorySignalCatalog.cpp:442)); the visible UI has no component selector, so the plotted component is not user-controlled.
- `strip.event` and `strip.category`: advertised in event/category helpers ([src/model/TrajectorySignalCatalog.cpp:151](../src/model/TrajectorySignalCatalog.cpp:151), [src/model/TrajectorySignalCatalog.cpp:194](../src/model/TrajectorySignalCatalog.cpp:194)); the strip renderer is numeric temporal paint, not a categorical/event lane ([src/app/StripStackWidget.cpp:86](../src/app/StripStackWidget.cpp:86)).
- RMSD advertises `strip.event`, which is semantically suspicious for a scalar time series ([src/model/TrajectorySignalCatalog.cpp:463](../src/model/TrajectorySignalCatalog.cpp:463)).

Mechanical cut: remove these from catalog helpers/descriptors unless a typed strip option exists for their semantics.

### UI Options To Remove Or Hide With The Hollow Modes

- Signal dialog mode kinds: Spectrum, Table, Color map, Glyph should be removed from the visible checkbox/filter set unless each has a registered visualisation ([src/app/SignalDisplayDialog.cpp:56](../src/app/SignalDisplayDialog.cpp:56), [src/app/SignalDisplayDialog.cpp:360](../src/app/SignalDisplayDialog.cpp:360), [src/app/SignalDisplayDialog.cpp:801](../src/app/SignalDisplayDialog.cpp:801), [src/app/SignalDisplayDialog.cpp:843](../src/app/SignalDisplayDialog.cpp:843)).
- Descriptor mode summary/tooltip should stop listing raw hollow mode ids once descriptors no longer carry them ([src/app/SignalDisplayDialog.cpp:437](../src/app/SignalDisplayDialog.cpp:437)).
- The "Metrics..." entry is valid only if the metric selector is limited to real visualisations. It is offered from the dashboard dock and main toolbar ([src/app/DashboardStripDock.cpp:92](../src/app/DashboardStripDock.cpp:92), [src/app/ReaderMainWindow.cpp:907](../src/app/ReaderMainWindow.cpp:907)).
- Expensive structural overlay toggles Butterfly and B-field are real, but they are not metric visualisations. If the lead wants a basic vetting screen, keep them out of the metric decision and consider hiding them from the default toolbar ([src/app/ReaderMainWindow.cpp:930](../src/app/ReaderMainWindow.cpp:930), [src/app/ReaderMainWindow.cpp:978](../src/app/ReaderMainWindow.cpp:978)).

## DECISIONS OWED TO THE LEAD

1. Must-have set: choose `Strips only`, or `Strips + Tensor Glyph`, or `Strips + Tensor Glyph + Atom Colour`, or `Strips + Tensor Glyph + Atom Colour + Sequence/Structural Bar`.
2. Keep/cut boundary: decide whether advanced dynamics panels (`static.curve.lag.animated`, `static.chord.coupling`, `static.fixed_freq`, `static.spectrum.power`) are part of the basic reader, advanced tools, or removed from the mode selector.
3. Tensor policy: decide whether T2/rank-2 tensor visibility is satisfied by strips for now, or whether tensor glyphs are mandatory before the reader is considered honest.
4. Class-system go/no-go: approve a small typed `VisualizationDefinition` registry so no mode can be advertised without a registered renderer, or accept a narrower cleanup that only deletes hollow strings from the current string-dispatch system.
5. Spectrum policy: either implement a real spectrum path using `SpectrumStripPanel`, or remove `strip.spectrum` everywhere it is offered.
