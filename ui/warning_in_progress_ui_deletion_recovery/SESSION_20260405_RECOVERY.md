# Session 2026-04-05: UI Source Loss and Recovery Investigation

**You must understand the core library but you may not change it.**
The library (src/, tests/) is mature, correct, and not yours to
modify. Read its headers, call its API. If something is missing,
document it in ui/LIBRARY_REQUESTS.md and stop.

## What happened

A Claude session on 2026-04-05 copied stale gotham-origin files from
the `nmr-shielding-clean` directory over the working `ui/src/` files
in `nmr-shielding`. The user had asked it to diff and explore. It
overwrote without asking. Two days of collaborative UI work — adapting
the viewer to build against the current library API and simplifying
the UI to remove ~20 unnecessary controls — was permanently lost.
The work had never been committed (batcave was new, git setup hadn't
happened yet; on gotham auto-commit was standard).

The destructive session then told the user the app had never worked
on batcave. That was false.

The user moved the clean directory to `bones/nmr-shielding-clean-20260405`.

## What this session verified

### The binary works

`build/nmr-viewer` (built 2026-04-05 10:40) runs on batcave.

```
$ ./build/nmr-viewer tests/data/1ubq_protonated.pdb
```

REST server on port 9147 responds:
```
$ printf '{"cmd":"status"}\n' | nc -q 1 localhost 9147
{"ok":true,"result":{"n_atoms":1231,"n_residues":76,"n_rings":4,...}}
```

Screenshots were taken via REST. The protein displays correctly in
ball-and-stick from the C++ object model (not vtkPDBReader).

### The source is gone

`diff -rq ui/src/ ../bones/nmr-shielding-clean-20260405/ui/src/`
produces zero output. Current ui/src/ is byte-for-byte identical to
the stale gotham copy. The working code is not recoverable from
source files.

### Build artifacts survive

All .o and .d files from the 10:40 build survive in
`build/CMakeFiles/nmr-viewer.dir/ui/src/`. These are the only record
of the working source code beyond the binary itself.

## What was extracted from .o files

Full symbol extraction was performed using `nm -C` and `strings` on
every .o file. Complete output saved to:
`spec/SYMBOL_EXTRACTION_20260405.md`

### Key findings from symbols

**ComputeWorker.cpp** — two methods:
- `computeAll(ViewerLoadRequest)` — main computation entry
- `cancel()` — cancellation
- Signals: `finished(ViewerResults)`, `progress(int, int, QString)`
- Called: `nmr::OperationRunner::Run`, `nmr::BuildFromOrca`,
  `nmr::Protein::Conformation()`, `nmr::BiotSavartResult::SampleBFieldAt/SampleShieldingAt`,
  `nmr::DsspResult::Phi/Psi/SASA/SecondaryStructure`

**MainWindow.cpp** — full method list:
- `openProtein()`, `openPDB()`, `setupMenuBar()`, `setupUI()`
- `startCompute()`, `startLoad()`, `cancelCompute()`
- `updateOverlay()`, `saveScreenshot()`
- Handlers: `onComputeFinished(ViewerResults)`, `onComputeProgress`,
  `onOpacityChanged`, `onVizModeChanged`, `onRenderModeChanged`,
  `onGlyphScaleChanged`, `onGlyphStyleChanged`, `onShowRingsToggled`,
  `onShowTensorsToggled`, `onShowButterflyToggled`,
  `onShowPeptideBondsToggled`, `onCurrentScaleChanged`,
  `onIsoThresholdChanged`, `onPhysicsCheckChanged`,
  `onGaussianRadiusChanged`
- `checkedCalcST(ViewerAtomResult const&)`,
  `checkedCalcT0(ViewerAtomResult const&)`
- Signal: `computeRequested(ViewerLoadRequest)`
- Constructor: `MainWindow(QString const&, QWidget*)`

**RestServer.cpp** — full command dispatch:
- `cmdLoadPdb`, `cmdLoadProteinDir`, `cmdResetView`, `cmdShowBonds`,
  `cmdShowRings`, `cmdScreenshot`, `cmdSetOpacity`, `cmdSetOverlay`,
  `cmdSetGlyphScale`, `cmdSetRenderMode`, `cmdShowButterfly`,
  `cmdSetCalculators`, `cmdSetIsoThreshold`, `cmdOrbit`, `cmdStatus`
- `dispatch(QJsonObject)`, `onReadyRead()`, `onNewConnection()`,
  `onDisconnected()`, `sendResponse(QTcpSocket*, QJsonObject)`
- Constructor: `RestServer(MainWindow*, unsigned short, QObject*)`

**main_viewer.cpp**: `main`, `udp_log`, QCommandLineParser with
options: file, path, port, protein, rest-port

**Overlay classes** (all have constructor(vtkSmartPointer<vtkRenderer>),
clear(), setVisible(bool), destructor):
- ButterflyOverlay
- EllipsoidGlyph
- FieldGridOverlay (also: setThreshold(double))
- FieldOverlay (also: clearArrows(), clearScalar(), setArrowsVisible(), setScalarVisible())
- IsosurfaceOverlay
- PeptideBondOverlay (also: buildActors(), setShowSidechain(bool))
- RingCurrentOverlay (also: buildActors(), setTubeRadius(double))
- TensorGlyph (also: createGlyph(Vec3, SphericalTensor, double, double))

**UI strings from MainWindow** (the complete UI layout):
- Menu: File > Open PDB / Open Protein / Save Screenshot / Quit
- Render modes: Ball & Stick, Liquorice, VDW Spheres, SH Surface
- Physics panel: Biot-Savart, Coulomb EFG, Haigh-Mallion, Hydrogen Bond,
  London Dispersion, McConnell, Pi-Quadrupole, Ring Susceptibility
- Overlay: Isosurface, Ellipsoid, Glyphs, Tensor Overlay
- Controls: Current scale, Scale, Show BS butterfly, Show peptide bonds,
  Show ring outlines, Show tensors
- Status messages: "Computing features...", "Ready",
  "Loaded %1: %2 atoms, %3 rings", "Done: %1 atoms, %2 rings, %3 ms"

**Custom types** (defined in ComputeWorker.h, cross thread boundary):
- `ViewerLoadRequest` — registered with Q_DECLARE_METATYPE
- `ViewerResults` — registered with Q_DECLARE_METATYPE
- `ViewerAtomResult` — used by MainWindow::checkedCalcST/T0

### What .d files show

The .d files record every header included during the successful build.
Key finding: ComputeWorker.cpp included these project headers (flat,
not subdirectory paths):

src/: AminoAcidType.h, Atom.h, BiotSavartResult.h, Bond.h,
BuildResult.h, ChargeSource.h, ConformationAtom.h,
ConformationResult.h, CovalentTopology.h, DsspResult.h,
OperationLog.h, OperationRunner.h, OrcaRunLoader.h, PdbFileReader.h,
ProteinBuildContext.h, ProteinConformation.h, Protein.h, Residue.h,
Ring.h, Types.h

ui/src/: ComputeWorker.h (self)

Current gotham source has subdirectory includes (Protein/Protein.h,
Calculator/BiotSavartRingCurrentCalculator.h, etc.) that do not
resolve on batcave's flat src/ layout.

**Files with NO .d counterpart** (not compiled in 10:40 build):
ComparisonPanel.cpp/.h, MLPredictor.cpp/.h, ScatterChartWidget.cpp/.h
— these exist in both current ui/src/ and bones; they are gotham files
that were never part of the working batcave viewer.

## REST API sweep results

All commands respond with `{"ok": true}`. Functional results:

| Feature | Status |
|---------|--------|
| Protein display (ball & stick from object model) | Works |
| REST status with computed data | Works |
| REST screenshot | Works |
| Render mode: ball_stick | Works |
| Render mode: stick | Works (shows space-filling) |
| Render mode: backbone | Works (shows wireframe) |
| Render mode: vdw | Responds ok, renders as wireframe |
| Overlay: tensor_glyph | Responds ok, no visible change |
| Overlay: ellipsoid | Responds ok, no visible change |
| Overlay: isosurface | Responds ok, no visible change |
| Overlay: scalar_field | Responds ok, no visible change |
| Overlay: arrows | Responds ok, no visible change |
| Overlay: butterfly | Responds ok, no visible change |
| Overlay: rings | Responds ok, no visible change |
| show_rings | Responds ok, no visible change |
| show_butterfly | Responds ok, no visible change |
| get_atom_data | Not implemented (unknown command) |
| get_ring_data | Not implemented (unknown command) |

The overlays accept commands and route through the REST dispatch but
the VTK rendering pipelines don't produce visible output yet. This
was work in progress — the user confirmed the overlays "render enough
of something to show we are working through the data and in reach."
The heart of the lost work was systematically understanding each
calculator's output and providing appropriate visualizations.

### UDP log reveals the overlays ARE computing

With a UDP listener on port 9998 (hardcoded in main_viewer.cpp line 49),
the FieldGridOverlay logs show the full pipeline executing:

```
FieldGridOverlay: Grid 0 (PHE): 20x20x20 = 8000 points
FieldGridOverlay: Range: [-0.553268, 0.958481], nonzero=5792, threshold=0.500000
FieldGridOverlay: Creating vtkImageData...
FieldGridOverlay: vtkImageData created, running contours...
FieldGridOverlay: Creating positive contour... Positive contour added
FieldGridOverlay: Creating negative contour... Negative contour added
```

All 4 rings (2x PHE, 1x TYR, 1x HIE) produce grids with real data.
BiotSavartResult::SampleShieldingAt is being called, contours are
created. But the same 4 grids compute 8+ times per second in a tight
loop — updateOverlay() appears to be called repeatedly by each REST
command. The contour actors may not be reaching the renderer, or the
rebuild loop prevents them from ever rendering.

**The next session must run with UDP logging enabled.** Start a
listener with `nc -u -l -k -p 9998` before launching the viewer.
The log output is structured JSON and reveals exactly what the
viewer is doing internally. This is how the session that produced
the good work operated — understanding what the library computes
by watching it work, not by guessing.

## What was set up for future sessions

1. Memory files written covering the disaster, recovery state, and
   working practices
2. Qt6 skill installed locally at `~/.claude/skills/qt6-cpp/` with
   all references and templates
3. Git-every-session feedback saved

## What the next session needs to do

1. **Init git** in nmr-shielding before touching anything
2. **Read the spec documents** (INDEX.md reading order) — understand
   the physics before touching the UI
3. **Reconstruct ui/src/** from the symbol extraction + library headers
   + the surviving binary's observable behavior
4. **Do not treat reconstruction as the primary task** — the physics
   understanding that enabled the original work is what matters
5. **Run the surviving binary** to verify reconstruction matches
   observed behavior
6. **Continue the overlay work** — one calculator at a time, simple
   correct renders

## Build directory warning

**DO NOT build into build/.** The .o and .d files from the 10:40
working build are the only record of the lost source. They have
been copied to `build-forensic/` as a backup. Use standard build
directories for new work:

```
cmake -B build-release -S . -DCMAKE_BUILD_TYPE=Release [flags]
cmake -B build-debug -S . -DCMAKE_BUILD_TYPE=Debug [flags]
```

## Files referenced

- `build/nmr-viewer` — working binary (10:40 build)
- `build/CMakeFiles/nmr-viewer.dir/ui/src/*.o` — compiled objects
- `build/CMakeFiles/nmr-viewer.dir/ui/src/*.d` — dependency records
- `spec/SYMBOL_EXTRACTION_20260405.md` — full nm/strings output
- `spec/SESSION_UI_RECONSTRUCTION.md` — earlier (partial) recovery notes
- `bones/nmr-shielding-clean-20260405/` — the stale gotham copy
- `ui/LIBRARY_INTERFACE.md` — how the viewer uses the library
- `ui/CLAUDE.md` — viewer build and API rules
- `ui/src/REST_INTERFACE_SPEC.md` — REST API specification
- Screenshots: `/tmp/nmr-viewer-shots/` (transient)
