# Polish backlog — what's left between now and viva-ready

**Framing.** Jessica: *this is a thesis-defense-quality tool for advisor
review of protein analysis before write-up. Surfacing complex data is
the fastest way to avoid folly.* The reader's role is not just viewing —
it is the editing / vetting surface for the H5 fields that will be
cited in the thesis argument. A clear per-metric glossary + appropriate
tensor / metric / atom colouring lets an advisor spot which fields are
compromised or restrictive before the model gets trained on them.

Work items below are ordered by their contribution to that framing.
"Next session" means the session after the current one, which the user
already earmarked for the time-series illustrator expansion.

---

## Thesis-defense critical (the actual content of the reader)

### 0. BS / HM frame-to-frame stacking (user-observed 2026-04-17)

During playback, the BS/HM isosurfaces visually **stack** — new frames'
surfaces draw with a few extra frames / extents of persistence still
visible. Reported as "very visible" with slow FPS.

**Code review — where stacking is NOT coming from.**
- `QtFieldGridOverlay` updates a single `vtkImageData` per ring in
  place: `SetOrigin`, `SetSpacing`, per-voxel `scalars->SetValue`, then
  `scalars->Modified()` + `imageData->Modified()`. One contour filter
  per ring per sign. No per-frame actor creation.
- Actor count stable at 63 across frames (verified in snapshot log),
  so this is NOT accumulated actors.

**Most likely mechanism.**
- `SetForceTranslucent(true)` + `SetOpacity(~0.3)` on the field grid
  actors, combined with FXAA (no depth peeling), gives us
  order-independent-translucency-without-actually-being-OIT. When the
  new frame's isosurface draws with alpha < 1, and the previous
  frame's fragments weren't fully overwritten by the opaque pass,
  translucent pixels from the prior frame leak through.
- `vtkContourFilter` output polygon count / vertex positions change
  slightly frame-to-frame (ring centre moves, scalars change). The
  changed polygons draw in slightly different pixels; pixels from the
  previous frame that no longer have a polygon *in the current frame*
  may retain their old colour due to translucency sort.

**Experiments, re-ordered 2026-04-17 given current hardware (5090).**

1. **Enable depth peeling for one run.** Depth peeling is effectively
   free on the 5090 — the qt6-cpp skill's FXAA-by-default guidance
   is calibrated for the Windows 8060S iGPU, not this machine. In
   `MoleculeScene::Build`, call `renderer_->SetUseDepthPeeling(true)`
   + `SetMaximumNumberOfPeels(4)` + `SetOcclusionRatio(0.0)`. If
   stacking disappears, the mechanism is translucency-without-OIT and
   the fix becomes a **runtime branch** (see
   `project_h5reader_target_hardware`): depth peeling on NVIDIA
   discrete + Apple M3 + Blackwell GB10, FXAA-only on AMD integrated.

2. **Force alpha=1 on the field grid for one run.** Set
   `rg.actorShielded->GetProperty()->SetOpacity(1.0)` (same for
   deshielded), `SetForceTranslucent(false)`. Confirms translucency is
   the mechanism (stacking should disappear). Diagnostic, not a fix —
   solid butterfly obscures molecule.

3. **Explicit pipeline Update() before Render.** In
   `QtFieldGridOverlay::setFrame`, call `rg.contourShielded->Update();
   rg.contourDeshielded->Update();` inside the ring loop. Rules out
   "contour output was stale on the first Render" as an alternative
   explanation. If stacking remains after (1) and (2) both pass, this
   deeper change is next.

Same architectural family as hypothesis (C) in `RESIDUAL_RENDER_DROP.md`
(FXAA + translucency pass ordering). If (0) and the residual drop
share a root cause, a single fix (depth peeling OR explicit Update())
addresses both — worth checking.

**Next-session priority.** Before extending the illustrator to per-atom
glyphs / colour bubbles (also translucent), understand THIS behaviour.
A per-atom colour overlay would inherit the same mechanism unless we
fix the pass-ordering story here first.

### 1. BS / HM butterfly extent controls

Currently the butterfly isosurfaces render with a fixed isovalue and
grid extent. A reviewer asking "what would this lobe look like at
half the threshold?" has no UI answer.

**What to add.**
- Toolbar spin-boxes or a dock slider for isovalue.
- Grid extent (e.g., 6–15 Å box half-width around the centroid) as a
  second slider.
- Debounced re-evaluation via `QtBiotSavartCalc` /
  `QtHaighMallionCalc`; kernel eval is not cheap, so a 250 ms debounce
  avoids burning cycles on slider drag.

**Where.** `src/app/QtFieldGridOverlay.{h,cpp}` already owns the
grid and mapper; add the isovalue/extent parameters there, a small
`QtFieldGridOverlayControls` dock or inline toolbar.

### 2. Per-field glossary / help

Every scalar the time-series dock surfaces, every label in the atom
inspector, needs a one-paragraph description sourceable from the H5
writer's own documentation. This is the "advisor-can-spot-compromised-
fields" feature; without it the advisor has to go to source to know
what `hbond_nearest_distance` means in the extractor.

**Motivating observation (2026-04-17):** user noticed many `T0`
scalars read zero for a picked atom. Physically this falls into
three buckets, and the glossary must distinguish them:

1. **Traceless by construction — T0 ≡ 0 always.**
   - `APBS EFG T0`, `coulomb EFG T0`, any EFG field — ∇·E = 0 in
     charge-free vacuum (Poisson). T0 is numerical noise.
   - `pq_shielding T0` if `pq` is a quadrupole-moment contribution —
     the moment tensor is traceless by definition.
2. **Local-only — T0 = 0 at atoms outside the cutoff.**
   - `bs_shielding T0`, `hm_shielding T0`, `rs_shielding T0` —
     ring-current, ~5–7 Å range from ring plane.
   - `mc_shielding T0` — McConnell bond anisotropy, local to polar
     bonds (C=O, C=N).
   - `hbond_shielding T0` — only atoms in an H-bond see nonzero.
3. **Expected nonzero everywhere.**
   - `coulomb_shielding T0`, `aimnet2Shielding T0`,
     `disp_shielding T0`. If zero at many atoms, that's a BUG to flag
     to the extractor maintainer.

Glossary MUST show which bucket each field belongs to, because
"T0 = 0 here" is three totally different stories.

**What to add.**
- Extend the `Desc` table in `QtAtomTimeSeriesDock.cpp` with a `help`
  field (one paragraph of plain English).
- Same table pattern in `QtAtomInspectorDock` for its field rows.
- Tooltips on labels, AND a help dock / modal tab that prints the full
  glossary with group-by-group grouping, so the advisor can read it
  without mousing over every field.
- The glossary itself probably lives in a single `docs/FIELD_GLOSSARY.md`
  that both the reader and the nmr_extract SDK refer to. That file is
  OUT OF SCOPE of `h5-reader/` but the reader should display the
  contents from a data-driven table inline — NOT by reading a .md file
  at runtime (that's fragile). Structured constants in a .cpp.

**Scope.** ~186 H5 datasets × ~35 currently exposed = the atom-scoped
minimum. Next session expands this alongside the scope-selector UI.
The glossary is the *source of truth* for which fields are in which
dock and how they are displayed.

### 3. Tensor display on atoms

T2 is sacred (feedback_t2_sacred). The 5-component spherical tensor
traceless symmetric piece IS the thesis argument. Right now the reader
surfaces `tensor.xx`, `tensor.xy`, etc., as scalars in the time-series
dock. That is not "displaying the tensor" — it is displaying its
components one at a time.

**What to add.**
- Per-atom tensor glyph in the 3D view: principal-axes ellipsoid or
  hedgehog oriented by the eigenvector frame, coloured by anisotropy.
- In the inspector dock, a compact matrix view of the T2 components
  with eigenvalues / eigenvectors shown inline.
- In the time-series dock, a tensor-specific plot mode: three
  eigenvalue traces + anisotropy / asymmetry in one panel.

**Where.** `src/model/Types.h` already has `SphericalTensor` with the
right fields. The overlay is new (`src/app/QtTensorGlyphOverlay.{h,cpp}`
— duplicate-and-adapt `QtRingPolygonOverlay` pattern). The time-series
dock extension is a new tab-row type in the Desc table.

### 4. Per-metric colouring on atoms (the "colour bubble" representation)

Called out in `TIME_SERIES_EXPANSION.md` axis 1. The whole-protein
scan-at-a-glance view: every atom rendered with its sphere coloured
and/or sized by a chosen scalar at the current frame.

**What to add.**
- New overlay `QtAtomColourOverlay` that modulates the existing
  `vtkOpenGLMoleculeMapper`'s per-atom scalar array. vtkMolecule
  supports per-atom scalars via `SetAtomicNumberArrayName` and LUT
  coupling; that pattern should carry a per-frame float array.
- Toolbar: a single dropdown "colour by [element | scalar]" with the
  same Desc table as time-series (identical glossary).

**Coupling to work item 3.** Tensor colouring is a natural case of
per-metric colouring (anisotropy → colour). Design the overlay so
"colour by scalar" and "colour by tensor-scalar-derived quantity" share
the same LUT plumbing.

---

## Deferred from code review (thesis-polish, not thesis-critical)

### 5. Inspector in-place update instead of `rebuild()`

`QtAtomInspectorDock::rebuild()` calls `clear()` + repopulates on every
pick OR frame change. Current scale is fine; when the time-series
expansion adds hundreds of rows, this will hot-spot.

**What to add.**
- Build the tree once per `setContext` / `setPickedAtom`.
- On `setCurrentFrame(t)`, iterate existing items and `setText(1, …)`
  only on the value column.
- Assert-only on topology changes (atom picked changed).

Est. 1–2 h. Matters most after work item 2 lands (many more rows).

### 6. Picker → selection direct wire

Currently when the picker emits `atomPicked`, the selection overlay
picks up the new atom, and `refreshCurrentFrame()` fires a full
overlay re-run (including kernel evals for field grid / streamlines).
The selection itself does not need those evals; only the selection
overlay needed the position update.

**What to add.**
- Separate path: `MoleculeScene::setSelection(size_t atomIdx)` that
  updates `selection_` only and issues a single `requestRender()` —
  skipping ribbon / rings / grid / streamlines.
- Keep `refreshCurrentFrame()` for the legitimate cases (visibility
  toggles where a hidden overlay comes back on).

Saves kernel evals on every pick. Small, clean.

### 7. Typed `vtkPolyDataMapper` member in `QtBFieldStreamOverlay`

`QtBFieldStreamOverlay` currently does `dynamic_cast<vtkPolyDataMapper*>`
on the actor's mapper to grab it for property changes. The mapper is
known at construction; store a `vtkSmartPointer<vtkPolyDataMapper>` as
a member and drop the cast.

Est. 10 min. No runtime change, just cleaner.

---

## Sequencing — the next 2–3 sessions before cross-platform

Jessica's plan (this session):

1. **Next session — time-series illustrator expansion.** Items 2, 3, 4
   land here (per-field glossary, tensor glyph/display, per-metric
   colouring bubbles). `TIME_SERIES_EXPANSION.md` is the scope map.
2. **Session after — controls + graph reps for matrix values.** Item 1
   (BS / HM extent controls), item 3 refinement (best graph
   representation for T2 tensor time evolution — three eigenvalue
   traces? anisotropy-vs-asymmetry scatter? Wigner-D coefficient
   heatmap? decide here).
3. **Third session — thorough review pass.** Items 5, 6, 7 land
   (inspector in-place update, picker direct wire, typed mapper
   member). Full review against the qt6-cpp skill references again.
   Known issues catalogue clean, UDP log clean on playback, crash
   path exercised.

**Only after all three: cross-platform build pass.**

### 8. Build on macOS and Windows — do NOT do during Linux sessions

The shared directory `/shared/2026Thesis/nmr-shielding/h5-reader` is
accessible from separate macOS and Windows Qt-Creator-or-CMake
sessions. Plan: once sessions 1–3 above are done and Jessica is happy
with the Linux reader, she will run separate sessions on Windows and
macOS that build against the same source tree. Linux-session work
should prepare for that — keep `#ifdef` platform gates explicit, keep
`CMakePresets.json` honest, avoid Linux-only system calls outside the
crash-diagnosis layer — but should NOT attempt to build on the other
two platforms from here.

What to do now in Linux sessions to make the cross-platform pass cheap:

- Any new file added to `CMakeLists.txt`: make sure the path and case
  are Windows-correct (no mixed case that collides with case-insensitive
  filesystems).
- Any direct `#include` path: use forward slashes.
- Any POSIX-ism (`/proc/self/statm` RSS probe, `<sys/socket.h>`,
  `<fcntl.h>`, signal self-pipe): keep it inside `#ifdef __linux__` or
  `#ifdef __unix__` with a documented no-op on Windows.
- File paths: `QStandardPaths` / `QDir` only, never literal POSIX paths.

### 9. Wire Windows MiniDumpWriteDump — during the cross-platform pass

`src/diagnostics/CrashHandler.cpp:171-174` is a stub on Windows. Qt6-cpp
skill's #1 principle is "when it crashes at 2am, can you read the
dump?" Wire this during the Windows build pass (session 4+), not
before — it is Windows-only and needs testing on a real Windows crash.

### 10. Residual render drop

Tracked in `RESIDUAL_RENDER_DROP.md`. Four hypotheses, ordered by
likelihood. Two cheap probes in place (duplicate-ResetCameraClippingRange
removed; `mapper_->Modified()` added with diagnostic comment). Next
time it recurs, follow the instrumentation playbook in that file.

---

## Principles guiding this backlog

- **Surface complex data to avoid folly.** The reader's purpose is to
  let the advisor see into the H5 well enough to spot what's wrong
  before it becomes a thesis claim.
- **Race conditions = deeper issue.** A "fix" that makes a race go
  away without explaining why it was racing in the first place is a
  band-aid. Probes are fine; undocumented band-aids are not.
- **Glossary is the source of truth.** Field meaning is authored once
  (in the reader's data tables) and referenced everywhere (tooltip,
  help dock, inspector row, time-series legend).
- **T2 is sacred.** Don't simplify it to a scalar anywhere — not in
  the inspector, not in the time series, not in the colouring. A T2
  colour mapping collapses 5 components to 1 by explicit design
  choice, not by UI convenience.
