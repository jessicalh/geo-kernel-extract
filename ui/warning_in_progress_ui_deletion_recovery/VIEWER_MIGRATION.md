# Viewer Migration: bs-viewer -> nmr-shielding viewer

**Purpose**: Rewire the existing Qt6+VTK molecular viewer from the old
biot-savart library (v2) to the new nmr-shielding ConformationResult
architecture. The viewer is working, pretty, and feature-rich. The job
is plumbing, not UI design.

---

## What exists (old project)

**Location**: `/mnt/extfast/2026thesis/biot-savart/src/Viewer/`

**17 source files**, one build target `bs-viewer`:

| File | What it does |
|------|-------------|
| main_viewer.cpp | Entry point. QApplication, CLI args, REST server, auto-load. |
| MainWindow.h/cpp | Qt main window. VTK renderer, all UI controls, 8 calculator checkboxes. ~800 lines. |
| ComputeWorker.h/cpp | **THE BRIDGE** to library. Loads protein, runs all 8 calculators per atom, produces ViewerResults. ~420 lines. Only file that includes library headers beyond Core/Types.h. |
| RestServer.h/cpp | TCP JSON REST API on port 9147. Load, camera, overlays, screenshots, turntable, sweep. ~300 lines. Full spec in REST_INTERFACE_SPEC.md. |
| ButterflyOverlay.h/cpp | B-field streamlines around rings (RK4-5, seeded at 1.5x ring radius). The signature visualization. |
| TensorGlyph.h/cpp | Per-atom spherical harmonic surface glyphs for T2 tensors. Red=deshielded, blue=shielded. |
| EllipsoidGlyph.h/cpp | Batched ellipsoid glyphs via vtkTensorGlyph. Fast path (one VTK actor for all atoms). |
| FieldOverlay.h/cpp | Scalar field (colored spheres) + E-field arrows at atoms. |
| FieldGridOverlay.h/cpp | Structured 3D grid isosurfaces from ViewerFieldGrid data. |
| IsosurfaceOverlay.h/cpp | Isosurface contours from T0 field grids. |
| RingCurrentOverlay.h/cpp | Aromatic ring outlines with normal arrows. |
| PeptideBondOverlay.h/cpp | Peptide bond visualization. |
| ComparisonPanel.h/cpp | Side panel for WT vs ALA DFT comparison display. |
| ScatterChartWidget.h/cpp | Qt Charts scatter plot (classical vs DFT). |
| Colormap.h | Moreland diverging blue-white-red colormap. |
| MLPredictor.h/cpp | Stub for ML predictions from .hpred files. |

**Build dependencies**: Qt6 (Widgets, OpenGLWidgets, Network), VTK 9.5
(IOChemistry, DomainsChemistry, DomainsChemistryOpenGL2, FiltersSources,
FiltersCore, FiltersGeneral, RenderingOpenGL2, RenderingCore,
RenderingAnnotation, GUISupportQt, InteractionStyle, ImagingHybrid,
FiltersFlowPaths, CommonColor, RenderingContext2D,
RenderingContextOpenGL2, ChartsCore, ViewsContext2D).

Links against: `nmr_v2` (old library), Qt6, VTK.

---

## What needs to change

### The ONLY file that touches the library is ComputeWorker

ComputeWorker.cpp is 420 lines. It includes the old library's headers
and calls the old API:

```
Old API (v2):                          New API (nmr-shielding):
─────────                              ──────────────────────
Protein::FromPdbFile()                 Same (Protein::FromPdbFile)
protein.ConformationAt(0)              protein.CrystalConformation()
conf.AtomPositions[]                   conf.AtomAt(i).position
conf.SetPartialCharges(charges)        ChargeAssignmentResult::Compute(conf)
conf.RingGeometries()                  conf.Result<GeometryResult>()
conf.EnsureBondGeometry()              (automatic from GeometryResult)
DsspCalculator::Calculate()            DsspResult::Compute(conf)
BiotSavartRCCalc::ShieldingTensorAtPoint()  conf.Result<BiotSavartResult>()
 ... (per atom, per calculator)        (results already on ConformationAtom)
```

The old viewer computes per-atom by calling each calculator's
`ShieldingTensorAtPoint(pt, protein, conf)` in a loop over atoms.
The new library pre-computes everything via ConformationResults that
write to ConformationAtom fields. The new ComputeWorker just reads
those fields.

### Migration strategy

1. **Copy all 17 files to `ui/src/`** (pure Qt/VTK code, no library changes)
2. **Rewrite ComputeWorker.cpp** (~420 lines) to use new API:
   - Build protein and conformation the new way
   - Attach results via dependency graph (or call a RunAllCalculators helper)
   - Read ConformationAtom fields instead of calling calculators per-point
   - Populate the same ViewerResults/ViewerAtomResult structs
3. **Update ComputeWorker.h** — remove old library includes, add new ones
4. **Everything else stays the same** — MainWindow, overlays, REST server,
   glyphs all work with ViewerResults, never touch the library directly

### What to watch for

- **ViewerAtomResult.element is a string** ("H", "C", etc.). The new
  library uses Element enum. The conversion is trivial but must happen
  in ComputeWorker.
- **ViewerAtomResult.atomName is a string**. Used for display only.
  Available from Protein's atom identity.
- **ViewerRingResult.ringType is a string** ("phe", "tyr", "trp6").
  Used for UI display. Get from ring.TypeName() or similar.
- **HeuristicClassifier** — exists in old code, may not exist in new.
  The tier system (REPORT/PASS/SILENT) is a viewer concern, not library.
  Keep it in the viewer or migrate as a simple helper.
- **Field grid computation** (Phase 8) evaluates BS at 20^3 grid points
  per ring. In the new library, BiotSavartResult writes to atoms, not
  arbitrary grid points. Either: (a) keep a ShieldingTensorAtPoint()
  method on BiotSavartResult for grid evaluation, or (b) compute the
  grid differently. The butterfly and isosurface overlays need this.
- **Charges**: old code loads from params file with string-matched
  residue/atom names. New code uses typed ChargeAssignmentResult from
  prmtop. The viewer needs a protein dir with prmtop, not a separate
  params file.

### What NOT to change

- MainWindow.cpp (800 lines of pure Qt/VTK — works as-is)
- All overlay classes (pure VTK — work with ViewerResults)
- RestServer (pure Qt networking — works with MainWindow slots)
- Colormap.h (standalone utility)
- ScatterChartWidget (Qt Charts — works with ViewerResults)

---

## ViewerResults data contract

The viewer's internal data model (ViewerResults) is the interface
between ComputeWorker and everything else. All 16 other files depend
on these structs. They are defined in ComputeWorker.h.

**Per atom** (ViewerAtomResult):
- Position, element, name, residue info (strings for display)
- Per-calculator T0 (8 doubles: bs, hm, mc, ld, ce, pq, rsa, hb)
- Total T0
- Full SphericalTensors (bs, hm, mc, total — for glyphs)
- B-field vector (for arrows)
- E-field vector (for arrows)
- DFT comparison (optional)
- Heuristic tier
- DSSP (secondary structure, phi, psi, SASA)

**Per ring** (ViewerRingResult):
- Center, normal, radius, intensity, residue, type string

**Per bond** (ViewerBondResult):
- Start, end positions, isPeptide, isSidechain

**Butterfly/field grids** (ViewerButterflyData, ViewerFieldGrid):
- 3D structured grids with B-field vectors or T0 scalars

The new ComputeWorker must populate the SAME structs. The rest of
the viewer doesn't care how they were computed.

---

## New ConformationAtom fields to map

**Verified against src/ConformationAtom.h (2026-04-02).**

After all ConformationResults are attached, each ConformationAtom has:

| ConformationAtom field | ViewerAtomResult target |
|------------------------|----------------------|
| .bs_shielding_contribution | bsST (direct), bsT0 (= .T0) |
| .hm_shielding_contribution | hmST (direct), hmT0 (= .T0) |
| .mc_shielding_contribution | mcST (direct), mcT0 (= .T0) |
| .disp_shielding_contribution | ldT0 (= .T0) |
| .coulomb_shielding_contribution | ceT0 (= .T0) |
| .coulomb_E_total | ceField |
| .piquad_shielding_contribution | pqT0 (= .T0) |
| .ringchi_shielding_contribution | rsaT0 (= .T0) |
| .hbond_shielding_contribution | hbT0 (= .T0) |
| .total_B_field | bsField |
| .total_G_tensor | totalST (decompose) |
| .Position() | position |
| .role | (convert to string for element/atomName) |
| .is_backbone | isBackbone |
| .hbond_is_donor | (available) |
| .hbond_is_acceptor | (available) |

Identity fields (element, atomName, residueName, residueNumber) come
from the Protein's atom/residue objects, not from ConformationAtom.

The total_G_tensor is the BS total only. For the viewer's totalST
(sum of all calculators), sum the 8 shielding_contribution tensors.

Note: old code has `totalST = Decompose(bsT + hmT + mcT + ldT + ceT + pqT + rsaT + hbT)`.
New code: sum the 8 `*_shielding_contribution` SphericalTensors, or
reconstruct Mat3s from each, sum, and decompose. Either works.

---

## Build integration

Add to nmr-shielding/CMakeLists.txt (conditional, like the old project):

```cmake
find_package(Qt6 COMPONENTS Widgets OpenGLWidgets Network QUIET)
find_package(VTK 9.5 COMPONENTS ... QUIET)

if(Qt6_FOUND AND VTK_FOUND)
    set(CMAKE_AUTOMOC ON)
    add_executable(nmr-viewer
        ui/src/main_viewer.cpp
        ui/src/MainWindow.cpp
        ui/src/ComputeWorker.cpp
        ... (all 17 files)
    )
    target_link_libraries(nmr-viewer PRIVATE
        nmr_shielding   # the new library target
        Qt6::Widgets Qt6::OpenGLWidgets Qt6::Network
        ${VTK_LIBRARIES}
    )
    vtk_module_autoinit(TARGETS nmr-viewer MODULES ${VTK_LIBRARIES})
endif()
```

The library builds as a static library (`nmr_shielding`). The viewer
links against it. No circular dependencies.

---

## Grid-point evaluation problem

The old viewer evaluates calculators at ARBITRARY 3D points (not just
atom positions) for butterfly streamlines and isosurface grids. The new
library's ConformationResult architecture computes at atom positions only.

Options:
1. **Add a public evaluate-at-point method** to BiotSavartResult (or a
   standalone function). This is the cleanest — the physics is the same,
   we just need to evaluate the kernel at a non-atom point.
2. **Pre-compute grids in a dedicated ConformationResult**. Heavy but
   consistent with the architecture.
3. **Keep the old per-point calculator call pattern** in a viewer-only
   helper that reuses the library's ring geometry but doesn't go through
   ConformationResult. Pragmatic but breaks the "one writer per field" rule.

Recommend option 1. The butterfly and isosurface are validation tools,
not features. A public method that evaluates the BS kernel at an arbitrary
point is honest and useful.

---

## REST API

The REST spec (REST_INTERFACE_SPEC.md) is viewer-internal. It doesn't
touch the library. Copy it as-is. The only commands that indirectly
touch the library are `load_pdb` and `load_protein_dir`, which call
through to ComputeWorker.

Future: the REST API could become an MCP server for Claude Code
integration. The JSON command structure is already MCP-compatible.

---

## Priority for first working build

1. Copy files to ui/src/
2. Rewrite ComputeWorker.cpp (the one file that matters)
3. Add CMake target
4. Build and verify: load a protein, see the molecule
5. Verify: overlay modes work (tensor glyphs, scalar field)
6. Verify: butterfly contours render around rings
7. Verify: REST API responds to status/screenshot commands
