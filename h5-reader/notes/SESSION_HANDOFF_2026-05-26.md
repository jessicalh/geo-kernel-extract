# h5-reader — session handoff (2026-05-26, killer-app increment 1 LANDED)

A long, productive session: the killer-app **design pass is settled** (pressure-tested with
the user) and **increment 1 — multi-atom selection + colour-coded spheres — is built, green,
and validated live on `:1`** (the user drove it: a 4-atom set via plain + Shift picks,
focus-tracked through 607 frames of playback, clean `rc=0`). The next mountain is **increment 2
— the measurement itself** (the graph paper made live).

---

## ⇩ PASTE PROMPT for the next session (copy this whole block)

```
Resume the h5-reader. Killer-app INCREMENT 1 is LANDED + validated and the design is settled —
this is a continue, not a restart. Deep-dive first (the project's standing rule); the living
state is the memory store + this handoff + the code. The OLDER notes/ design docs (REVIEW_BRIEF_*,
H5_READER_REWRITE_DESIGN, FEATURE_PLAN, VISION_AND_PROGRESS) are HISTORY from before the rewrite.

UNDERSTAND FIRST (read, in this order):
1. notes/SESSION_HANDOFF_2026-05-26.md (this) — landed state + the increment-2 plan.
2. Memory: project_h5reader_killer_app_multiatom_compare_20260526 (THE doc — vision + the
   SETTLED DESIGN + INCREMENT 1 LANDED + the Franklin "graph paper" thesis line),
   project_h5reader_rewire_landed_20260526 (the Conformation-base foundation),
   project_h5reader_buffering_decision_20260526 (why Conformation is shaped this way — do NOT
   relitigate), project_h5reader_display_surface_20260526 + project_h5reader_tier_mirror_audit_20260526
   (the later scope-views arc + the one-source-per-role tier rule).
3. Code: src/model/AtomSelection.{h,cpp} (the selected group, as a QAbstractListModel),
   src/model/ConformationGeometry.{h,cpp} (where increment 2's Distance/Angle/Dihedral go,
   beside RingGeometryAt), src/app/MeasurementOverlay.{h,cpp} (the ≤4 spheres — grows the
   lines + labels), src/app/SelectionDock.{h,cpp}, and the wiring block in
   src/app/ReaderMainWindow.cpp.
Qt/VTK discipline: CENSUS_REGISTER every QObject ctor, ACONNECT every connection, ASSERT_THREAD
on GUI methods, QPointer<T> for cross-object QObject* members.

STATE (verified 2026-05-26): the design pass settled two spine decisions — (1) ONE unified,
focus-tracked AtomSelection (plain double-click = replace the focus; Shift+double-click = toggle
the atom in the ≤4 ordered set), and (2) the user chose the FULL QT MODEL: AtomSelection IS the
QAbstractListModel (the group = the model, one source of truth), exposing rows via roles for
views AND typed atoms()/focus()/geometryKind()/SlotColor accessors for the renderer. Increment 1
LANDED: AtomSelection + MeasurementOverlay (≤4 Okabe-Ito colour-coded spheres, repositioned per
frame from atomPosition) + SelectionDock (the QListView panel, the model's first view); the
picker grew atomPicked(idx, Qt::KeyboardModifiers); ReaderMainWindow rewired picker → AtomSelection
→ {focusChanged: inspector + time-series; cleared: their clearSelection; changed: overlay +
scene refresh}. Built green on linux-rwdi (no warnings); validated live on :1 (clean rc=0).
QtSelectionOverlay (the old single yellow sphere) is left CONSTRUCTED-BUT-DORMANT — retire it
when increment 2 lands, but ASK before deleting. RUN ON A REAL DISPLAY (DISPLAY=:1); offscreen
segfaults (no FBO). clangd "QObject/Eigen/VTK not found" cascade is editor noise — the cmake
build is truth.

NEXT = INCREMENT 2 (the measurement — the graph paper made live):
1. ConformationGeometry gains Distance (Å), AngleDegrees (vertex = middle atom, [0,180]), and
   DihedralDegrees (signed Blondel-Karplus atan2, (−180,180]); GeometryMeasurement {kind, value,
   valid} + Measure(conf, frame, atoms) dispatching on count. All from Conformation::atomPosition
   → works for a pose and a trajectory.
2. MeasurementOverlay grows a connecting POLYLINE through the ordered atoms + a
   vtkBillboardTextActor3D VALUE LABEL — anchored at the midpoint (distance) / vertex (angle) /
   central-bond midpoint (dihedral) — both re-read every frame so the number + geometry HOLD as
   the molecule rotates. Add RenderingFreeType to the VTK COMPONENTS in cmake/DiscoverDeps.cmake
   for the label (RenderingAnnotation is already present and pulls it transitively, but make it
   explicit).
3. DEFER the prettier angle-arc / dihedral-wedge glyph — connecting lines + label first.
THESIS HOOK (free here): click ω's four backbone atoms and the measured dihedral should match the
inspector's extracted omega_actual at that frame — measured-from-positions vs extracted-field, the
reader's whole vetting purpose, in one gesture. Handle the ±180° wrap in any over-time plot.
LATER ARCS (NOT now): inc-3 derived geometric series in the dock; inc-4 comparison across the set;
inc-5 tabbed scope views + the data dictionary (display_surface memory). Bulk de-Qt rename of the
~115 non-QObject Qt* types is still pending (separate deliberate pass; new types are born honest).

DEFERRED BUG (orthogonal to the killer app): the molecule's ball-and-stick atoms intermittently
DROP in some playback frames, leaving only the overlays. This is the PRE-EXISTING
RESIDUAL_RENDER_DROP (a vtkOpenGLMoleculeMapper VBO re-upload miss; notes/RESIDUAL_RENDER_DROP.md +
the probe comment at MoleculeScene.cpp:316) — NOT caused by increment 1 (which only adds actors +
a model, never touches the molecule mapper). User said "figure out later." Worth fixing before
viva. NOT on this side: the catalog drift (hbond_scalars 3-vs-4, gromacs_energy 43-vs-42) lives in
python/nmr_extract/_catalog.py.

RULES: work only in h5-reader/; do NOT commit (no git ownership here; the user manages git).
atomPosition is the one position seam. New/rewritten types are born honest (no Qt prefix). Writer
is definitive for NPY/H5 shape. Ask before deleting/overwriting (incl. retiring QtSelectionOverlay).
Run on a real display.
```

---

## Where we are (2026-05-26) — landed

- **Reader foundation** (earlier today): the `Conformation` base + `TrajectoryConformation` /
  `SingleConformation`; the one position seam `atomPosition`; both run shapes render on `:1`.
  Memory `project_h5reader_rewire_landed_20260526`.
- **Killer-app design settled**: ONE unified, focus-tracked `AtomSelection`; the user chose the
  **full Qt model** (`QAbstractListModel`) for the selected group; heart-first sequencing.
  Memory `project_h5reader_killer_app_multiatom_compare_20260526` (THE doc — vision + settled
  design + increment-1 landed + the Franklin line).
- **Increment 1 LANDED + validated**: `AtomSelection` (the `QAbstractListModel`),
  `MeasurementOverlay` (≤4 Okabe-Ito spheres), `SelectionDock` (the panel); picker →
  `atomPicked(idx, modifiers)`; `ReaderMainWindow` rewired. Built green (no warnings); driven
  live on `:1` (4-atom set, 607 frames of focus-tracked playback, clean `rc=0`).

## What's next

1. **Increment 2 — the measurement** (detail in the paste prompt): `Distance` / `Angle` /
   `Dihedral` + `GeometryMeasurement` + `Measure()` in `ConformationGeometry`;
   `MeasurementOverlay` grows a connecting polyline + a `vtkBillboardTextActor3D` value label
   held through rotation. Thesis hook: measured ω vs extracted `omega_actual`.
2. Then inc-3 (geometric series in the dock), inc-4 (comparison across the set), inc-5 (tabbed
   scope views + data dictionary — memory `project_h5reader_display_surface_20260526`). Bulk
   de-`Qt` rename still pending.

## Deferred bug (NOT the killer app)

The molecule's ball-and-stick atoms intermittently drop in some playback frames, leaving only
the overlays — the **pre-existing `RESIDUAL_RENDER_DROP`** (`vtkOpenGLMoleculeMapper` VBO
re-upload miss; `notes/RESIDUAL_RENDER_DROP.md` + the probe comment at `MoleculeScene.cpp:316`),
not caused by increment 1. User: "figure out later." Worth fixing before viva.

## Doc note

`README.md` + `SCOPE.md` updated 2026-05-26 to mention the `AtomSelection` model + the
multi-atom selection. `CLAUDE.md` rules are unchanged (the killer app is additive; scope +
discipline hold). The older `notes/` design docs (REVIEW_BRIEF_*, H5_READER_REWRITE_DESIGN,
FEATURE_PLAN, VISION_AND_PROGRESS) are point-in-time history — the memories + this handoff +
the code are the living state.
