# h5-reader — Feature Plan

Per-feature planning for the work that turns the reader into the
thesis-data visibility surface. Paired with
`H5_FIELD_GLOSSARY.md`, which tells you *what* to show; this
document tells you *how* to build the renderers and docks that
show it.

**Date:** 2026-04-17. **Status:** working plan, **tentative**.

**Why tentative.** The 10-protein calibration DFTs currently
running against `/shared/2026Thesis/fleet_calibration-working/`
(project state `project_three_stages`, Stage 2) are the first
real test of which kernels carry signal against experiment at
the per-frame level. Their outcome will re-rank the features
in this plan: if water EFG or AIMNet2 EFG prove load-bearing,
their overlays and docks move up; if the ridge residual on
per-frame DFT is dominated by a single physics group, the
"predictions stacked-bar decomposition" dock (F14) becomes the
headline feature rather than a specialist view; if the MLP
nonlinear signal is confined to a narrow subset of kernels, the
55-kernel heat strip (F15) may shrink. Infrastructure I1–I8 is
largely DFT-insensitive and starts now; features that depend
on specific calibration outcomes will firm up once the DFTs
complete and the Phase 2 analysis lands. Re-rank this plan at
that point; the tracks and the refinement cadence stay.

## Intent

The glossary assigns every H5 field a primary presentation
modality across seven modalities (tensor glyph, colour bubble,
SAS surface, ribbon thickness, volumetric isosurface, n-D dock,
scalar line chart). The existing reader covers modality 7
(scalar line chart) and fragments of modality 6 (the BS/HM
butterfly and ring B-field streamline are in scene, not in a
dock, but count as "abstract view of a physical field"). The
other five modalities are missing.

**Approach: shared primitives first, features second.** Colour
bubble, SAS surface, and ribbon thickness all share a scalar
selector + colormap. Superquadric and Haeberlen tensor glyphs
share an eigendecomposition. Every dock that is not atom-scoped
shares a scope enum and picker-routing logic. The plan builds
those shared foundations once, then features become mechanical
table entries against them.

## How to read this document

Each section has the same shape:

- **Goal** — one sentence.
- **Unlocks** — which H5 fields this surfaces (from the
  glossary's Index). For infrastructure items, which features
  it enables.
- **Classes / files** — class names to create, pattern to copy,
  file locations under `src/`.
- **Data path** — what reads the H5, what transforms, what
  renders.
- **UI** — where it lives in the window, control surface.
- **Diagnostics** — UDP log category per
  `feedback_udp_logging`, sanity checks.
- **Depends on** — infrastructure or earlier features.
- **Effort** — S (≤ 1 session), M (2–3 sessions), L (3–5
  sessions). Sessions are ~4-hour Claude sessions.
- **Open questions** — what to decide in-session or punt.

Infrastructure and features are both numbered. Infrastructure
uses `I#`, features use `F#`.

**Tracks.** Work splits into four tracks that can progress in
parallel once the infra layer is in place:

- **Track A — In-scene overlays.** Renderers that draw on or
  around the molecule. `F1–F8`.
- **Track B — Abstract docks.** Separate windows with charts or
  projections. `F9–F18`.
- **Track C — Navigation + inspector.** Interactive glossary
  and inspector small-multiples. `F19–F20`. The glossary here is
  not a passive doc viewer but a first-class navigational
  surface: search, current-value panel, cross-links to the
  overlay / dock / inspector row that shows each field.
- **Track R — Refinement.** One dedicated session after every
  2–4 feature sessions to harmonise controls, colours,
  defaults, error states, and cross-platform appearance. `R1–R5`.
  See the Refinement cadence section below for the template.

Cross-platform bar (Linux MSVC / macOS Clang / Windows MSVC)
carries through per `project_h5_reader_scope`. All Qt-managed
file I/O, `QUdpSocket` networking, no platform `#ifdef` outside
the diagnostics layer.

---

# Summary table

| ID  | Feature                                 | Unlocks (fields)         | Track | Effort | Depends on |
|-----|------------------------------------------|--------------------------|-------|--------|------------|
| I1  | Scalar field registry + colormap utility | ~60 per-atom scalars     | infra | M | — |
| I2  | Tensor machinery (eigen, Haeberlen, SQ)  | ~25 tensor fields        | infra | M | — |
| I3  | Vector machinery (arrow norm + presets)  | ~10 vector fields        | infra | S | — |
| I4  | Per-residue aggregator                   | derived per-residue      | infra | S | — |
| I5  | SAS mesh primitive                       | SAS surface colouring    | infra | M | — |
| I6  | Scope framework (Atom/Residue/Ring/Frame)| all docks                | infra | M | — |
| I7  | Field-glossary lookup (from this repo)   | tooltips, dropdown labels| infra | S | — |
| I8  | Per-overlay UDP diagnostics pattern      | debuggability            | infra | S | — |
| F1  | Atom colour bubble                       | ~30 scalars              | A     | S | I1, I8 |
| F2  | SAS surface colouring                    | Mahalanobis + residuals  | A     | M | I1, I5 |
| F3  | Ribbon thickness modulation              | hbond_energy + derived   | A     | M | I1, I4 |
| F4  | Superquadric tensor glyph (shielding)    | 13 shielding tensors     | A     | M | I2 |
| F5  | Haeberlen principal-axis glyph (EFG)     | 9 EFG fields             | A     | S | I2, F4 |
| F6  | Vector arrow overlay                     | 10 vector fields         | A     | S | I3 |
| F7  | Relationship line segment                | 2 direction fields       | A     | S | F6 |
| F8  | Protein-wide \|B\| volumetric isosurface | total_B_field            | A     | M | BS evaluator|
| F9  | Time-series scope expansion              | residue/ring/system scope| B     | M | I6 |
| F10 | Ramachandran dock                        | phi, psi                 | B     | S | I6 |
| F11 | DSSP residue × time heatmap              | ss8                      | B     | S | I6 |
| F12 | Per-ring-type bar chart                  | 4 per-type T0 fields     | B     | S | I6 |
| F13 | K=6 nearest-ring bar chart               | 7 per-ring fields        | B     | M | I6 |
| F14 | Predictions stacked-bar decomposition    | raw_*/norm_* groups      | B     | M | I6 |
| F15 | MLP 55-kernel heat strip                 | mlp_kernel_weights       | B     | S | I6, I7 |
| F16 | AIMNet2 embedding projection (UMAP/PCA)  | aim, subspace_coords     | B     | L | I6 |
| F17 | System-energy dock                       | /energy/* (42 lines)     | B     | M | I6 |
| F18 | 3×3 tensor matrix view (virial, pressure)| virial, pressure_tensor  | B     | S | I6 |
| F19 | Inspector small-multiples for T2 stacks  | per-category tensor      | C     | M | I2 |
| F20 | Glossary Navigator (schema nav + search) | every H5 field (nav hub) | C     | M | I1, I6, I7 |
| R1  | Refinement — post-infra (defaults, chrome)| consistency             | R     | S | after infra |
| R2  | Refinement — post-scalar + tensor families| polish                  | R     | S | after F1–F5 |
| R3  | Refinement — post-glossary + vectors + vol| polish                  | R     | S | after F6–F8, F20 |
| R4  | Refinement — post-simple docks            | polish                  | R     | S | after F9–F12, F18 |
| R5  | Refinement — final (incl. F19, tooltips)  | polish                  | R     | M | after F13–F17 |

---

# Track 0 — Shared infrastructure

Build these first. Every feature past this point assumes they
exist.

## I1 — Scalar field registry + colormap utility

**Goal.** One table that names every per-atom, per-residue, and
per-frame scalar field in the H5, with the UI label, units,
default colormap and range, and the function that extracts the
value from `QtFrame` or `QtConformation`. Feeds three overlays
(F1, F2, F3) and the scope-expanded time-series (F9).

**Unlocks.** F1 (colour bubble), F2 (SAS surface), F3 (ribbon
thickness), F9 (scope expansion), F15 / F19 (label lookups),
F20 (tooltip content).

**Classes / files.**
- `src/model/ScalarField.h` — `struct ScalarField { QString id;
  QString label; QString units; Scope scope; ColormapPreset cmap;
  QPair<double,double> range; ExtractFn fn; }`. `ExtractFn` is a
  `std::function<double(const QtFrame&, size_t)>` so the index
  argument is atom / residue / ring depending on scope.
- `src/model/ScalarRegistry.h/.cpp` — singleton-ish registry
  populated at startup via free `registerScalars()` functions
  per group (one per H5 group, mirror the existing
  `QtAtomTimeSeriesDock::descAt()` pattern so migration is a
  mechanical refactor).
- `src/model/Colormap.h` — extend current `MakeDivergingCTF`
  with `MakeSequentialCTF(viridis / plasma / inferno)` and a
  `ApplyPresetToMapper` helper.

**Data path.** Registry built at file-load from
`QtConformation`; read via `ScalarRegistry::find(id)` or
`ScalarRegistry::all(scope)`. Consumers bind to `QComboBox` via
`scalarRegistry.toComboModel(Scope::Atom)`.

**UI.** Not user-visible directly. Backs the dropdowns in F1,
F2, F3, F9.

**Diagnostics.** `h5reader.scalar` category. Log field count at
startup, log invalid id lookups. Dump full registry on UDP
command (future admin hook).

**Open questions.** Do derived fields (ridge residual, per-
residue sums) live in the same registry or a sibling? Default
yes, flagged with `source = Derived` so inspectors and docs
handle them consistently.

## I2 — Tensor machinery

**Goal.** A small library for converting the `(9,)` spherical-
tensor layout into the renderings the glyphs need: 3×3 Cartesian
matrix, eigenvalues / eigenvectors, Haeberlen ordering, η, span,
skew, and the superquadric shape parameters (c_l, c_p from
Kindlmann γ=2.5 convention).

**Unlocks.** F4 (superquadric glyph), F5 (Haeberlen glyph), F19
(inspector small-multiples), eigenvalue time-series entries
added to F9.

**Classes / files.**
- `src/model/SphericalTensor.h/.cpp` — already used in the
  inspector as a container; extend with:
  - `Mat3 toCartesian(const double* st9) const` — reconstruct
    the 3×3 from T0 + T1 + T2 using the T0/T1/T2 to Cartesian
    transform (documented in the calculators; codify here).
  - `Eigendecomp eig(const Mat3&)` — Jacobi or wrap Eigen.
    Eigen is header-only, already a plausible dep.
  - `HaeberlenOrdering haeberlen(const Eigendecomp&)` —
    `|V_zz| ≥ |V_yy| ≥ |V_xx|`, η = (V_xx − V_yy) / V_zz.
  - `SuperquadricParams kindlmann(const Eigendecomp&, double
    gamma = 2.5)` — c_l, c_p parameters for
    `vtkSuperquadricSource` per the 2004 paper.
- Extract T0 / T1 / T2 helpers already exist; keep the irrep
  split explicit in the API — no opaque 9-vectors past this
  file.

**Data path.** Overlay reads `(T, N, 9)` slab from `QtFrame`
via existing accessor, converts to `Mat3`, eigen-decomposes,
builds glyph parameters, caches per-frame. No per-glyph re-
evaluation of the underlying physics (that stays in H5).

**UI.** None directly. The glyphs and the inspector consume it.

**Diagnostics.** `h5reader.tensor` category. Log any atom where
`|Tr(V)| > 1e-8` for a dataset expected to be traceless (EFGs);
advisers will spot drift this way.

**Open questions.** Vendor Eigen (MIT) or hand-roll a 3×3 Jacobi?
3×3 real symmetric Jacobi is ~80 lines and avoids a dep; prefer
hand-rolled. Full Mat3 (asymmetric) eigendecomp for shielding
tensors is harder but feasible via Schur; for the glyphs we
only need the symmetric T2 part so Jacobi suffices.

## I3 — Vector machinery

**Goal.** Standard per-atom arrow rendering: magnitude → arrow
length (with cap and log option), direction from the Vec3,
colour preset per source field. Shared by F6 / F7.

**Unlocks.** F6 (vector arrow overlay), F7 (relationship line
segment). Also supports any future per-atom gradient fields.

**Classes / files.**
- `src/app/ArrowRenderer.h/.cpp` — wraps `vtkGlyph3D` +
  `vtkArrowSource`. Input: per-atom Vec3 array, colour-by
  choice (magnitude scalar OR a fixed preset colour).
- `src/app/VectorFieldPresets.h` — named preset colours by
  source: `BackboneGrey`, `SidechainTan`, `AromaticPurple`,
  `SolventCyan`, `WaterBlue`, `BFieldGold`, etc. One source,
  one colour across the whole app.
- `QtFrame` additions where needed (most Vec3 accessors exist;
  audit for the full list in the glossary's Vector-arrow Index
  section).

**Data path.** Overlay reads `(T, N, 3)` slab, feeds into
`ArrowRenderer` which owns its `vtkPolyData + vtkGlyph3D +
vtkActor` chain. Arrow scale adaptive: `scale = min(|v| / max|v|,
cap_ratio)` with `cap_ratio` UI-tunable.

**UI.** Selector in the overlay dropdown. Legend actor with
colour + scale example.

**Diagnostics.** `h5reader.arrow` category. Log max |v| per
frame so scale sanity is visible on UDP.

**Open questions.** Should arrows render on all atoms or only
selected / filtered subset? At 5000 atoms with 10 vector fields
the clutter is real. Default: render only on picked atom + its
first-shell neighbours; have a "show all" toggle for protein-
scale views.

## I4 — Per-residue aggregator

**Goal.** Given a per-atom array and an aggregation rule
(sum / mean / max / min), produce a per-residue array suitable
for ribbon thickness and per-residue time series. Shared
between F3, F9, and any derived per-residue scalar in F20.

**Unlocks.** F3 (ribbon thickness), F9 scope expansion to
residues, per-residue derived scalars in the SAS surface overlay.

**Classes / files.**
- `src/model/ResidueAggregator.h/.cpp` — takes
  `const QtConformation&` plus a scalar id and aggregation rule,
  returns `std::vector<double>` of length R. Caches last result.
- `Aggregation` enum: Sum, Mean, Max, Min, RMSD (future).

**Data path.** Reads per-atom slab, reads `/atoms/residue_index`,
folds. O(N) per frame. Cacheable per (field_id, aggregation,
frame) if profiling says so.

**UI.** Consumed via a "per-residue summary" picker when a
scope-Residue Desc entry is a derived aggregation.

**Diagnostics.** `h5reader.aggregator` category, log cache
misses per frame.

**Open questions.** Expose aggregation rule in UI, or curate one
default per field? Default curate in the registry (the glossary
already picks: `bonded_energy/total` wants sum, Mahalanobis
wants mean, etc.).

## I5 — SAS mesh primitive

**Goal.** Build a solvent-accessible surface mesh of the
protein at each frame and expose it as a `vtkPolyData` with per-
vertex scalar mapping.

**Unlocks.** F2 (SAS surface colouring) and any future surface-
based visualisation.

**Classes / files.**
- `src/app/QtSasMesh.h/.cpp` — owns the build + cache:
  - Path A: `vtkMolecularSurfaceFilter` (Connolly / SES).
  - Path B: marching cubes over a Shrake-Rupley density grid
    built to match `SasaResult` (`SasaResult::ProbeRadius`
    matches `vtkMarchingCubes`).
- `src/app/ScalarToVertex.h` — maps per-atom scalar → per-vertex
  scalar by nearest-atom or inverse-distance-weighted.

**Data path.** Build once at file-load for a reference frame,
per-frame rebuild or deformation per frame (benchmark: 2000
atoms should mesh in ~100 ms via marching cubes at 1 Å grid).
Scalar mapping: iterate mesh vertices, find nearest atom via
`vtkOctreePointLocator`, sample the scalar.

**UI.** Consumed by F2.

**Diagnostics.** `h5reader.sasmesh` category. Log mesh vertex
count per frame, build time per frame, nearest-atom lookup
misses.

**Open questions.** Per-frame vs every-N-frames vs reference-only
re-mesh. Per-frame is right for accuracy; every-5-frames is a
plausible fallback if the cost is too high. Decide after first
benchmark on a 2000-atom protein.

## I6 — Scope framework

**Goal.** Let UI controls and Desc-table entries tag themselves
with the scope of the data they represent (Atom, Residue, Ring,
System), and route picker events / index arguments accordingly.

**Unlocks.** Every dock past F9. Existing `QtAtomTimeSeriesDock`
refactored to use the framework rather than its private Desc
table.

**Classes / files.**
- `src/model/Scope.h` — `enum class Scope { Atom, Residue, Ring,
  Frame, System };` and `struct ScopedIndex { Scope s; size_t
  idx; };`.
- `src/app/QtResiduePicker.h/.cpp`, `src/app/QtRingPicker.h/.cpp`
  — new pickers emitting `residuePicked(size_t)`,
  `ringPicked(size_t)`. Ray-cast hit-test against ribbon /
  ring polygons respectively.
- Refactor `ScalarField` in I1 to carry a Scope tag.
- Update the existing picker to `QtAtomPicker` explicitly and
  document that it owns atom-scope events only.

**Data path.** Pickers emit `ScopedIndex`; any dock can subscribe
to the right channel. Cross-scope queries (e.g. "click a residue,
highlight all its atoms") go via the aggregator.

**UI.** Scope picker chips or a scope radio in each dock header.

**Diagnostics.** `h5reader.scope` category. Log every pick with
scope + index + position for UDP-tailable debugging.

**Open questions.** Does the existing `atomPicked(size_t)` signal
keep its name (back-compat) or migrate to the `ScopedIndex`
channel? Keep the name, add `ScopedIndex` channel alongside, so
the refactor is non-breaking.

## I7 — Field glossary lookup

**Goal.** Parse `H5_FIELD_GLOSSARY.md` (or a derived JSON) into
an in-memory lookup so tooltips, inspector headers, and dropdown
labels pull descriptions from the one authoritative source.

**Unlocks.** F20 (tooltips), consistent UI labels across F1–F19.

**Classes / files.**
- `src/model/FieldGlossary.h/.cpp` — loads at startup. Input
  options:
  - Parse the markdown directly (fragile; depends on doc
    structure).
  - Generate `doc/field_glossary.json` at build time via a
    Python extractor script (`scripts/build_glossary.py`) and
    ship the JSON alongside the binary.
  - Hand-maintain a compact `.toml` with a `[[field]]` entry
    per dataset.
- Prefer the JSON-from-markdown extractor: the glossary stays
  the single source of truth, CI regenerates the JSON.

**Data path.** On load, read JSON into
`std::unordered_map<std::string, FieldDoc>`. Lookups by H5 path
(`"/projections/mahalanobis"`).

**UI.** Consumed via tooltip on any control that binds to a
scalar field; dedicated help dock in F20.

**Diagnostics.** `h5reader.glossary` category. Log
missing-field-description warnings.

**Open questions.** Extractor script in `scripts/` or inside
`h5-reader/`? Put it in `h5-reader/scripts/` since h5-reader
owns the consumer.

## I8 — Per-overlay UDP diagnostics pattern

**Goal.** Codify the diagnostic log lines every overlay should
emit so advisers can debug from the UDP stream without opening
the source. Per `feedback_udp_logging` and
`feedback_qt_discipline`, UDP is the primary debug channel.

**Unlocks.** Cross-cutting reliability.

**Pattern.** Each overlay class declares a category
(`h5reader.overlay.<name>`) and emits at minimum:

- at construction: `CENSUS_REGISTER`, log "created".
- at `Build`: atom / ring / residue count, field(s) bound.
- at `setFrame`: one line with frame index + summary statistics
  of what rendered (e.g. min/max of the bound scalar).
- at `setVisible`: log on/off transitions.
- at boundary failures: `try/catch` at HighFive seams, log
  specific field + frame + error, continue with partial state.

**No per-render-call spam.** Per-frame statistics go at DEBUG
level; errors at WARNING / ERROR.

---

# Track A — In-scene overlays

Order: F1 → F2 → F3 (scalar family), then F4 → F5 (tensor
family), then F6 → F7 (vectors), then F8. F1–F3 share the
scalar registry so they land nicely together.

## F1 — Atom colour bubble

**Goal.** Paint the atom spheres of the `vtkMolecule` with a
scalar selected from the registry, with a legend.

**Unlocks.** 30+ fields in the "Atom colour bubble" Index
section of the glossary.

**Classes / files.**
- `src/app/QtAtomColourOverlay.h/.cpp` — overlay contract
  (`Build` / `setFrame` / `setVisible`).
- Reuses `vtkMolecule::GetAtomicNumberArrayName` and a custom
  `vtkFloatArray` of per-atom scalars; point the mapper's
  lookup table at the array.
- `vtkScalarBarActor` for the legend.

**Data path.** On `setFrame(t)`: fetch the selected
`ScalarField`, iterate `0..N`, fill the float array from
`field.fn(frame, i)`, set the array on the mapper, mark it
modified. No re-build needed per frame.

**UI.** Toolbar action "Colour by…" with a dropdown of all
atom-scope scalars from I1; colormap + range picker in the
overlay options (mid-priority; range can default from the
registry).

**Diagnostics.** `h5reader.overlay.atomcolour` per I8.

**Depends on.** I1, I8.

**Effort.** S.

**Open questions.** Element colouring is the default; "off"
restores it. Should multiple scalars be paintable simultaneously
(alpha-blended)? No for v1; toggle one at a time.

## F2 — SAS surface colouring

**Goal.** Draw the SAS mesh (from I5), paint it by a scalar
from I1.

**Unlocks.** Mahalanobis, residual_*, hydration alignment /
coherence, bonded_energy/total painted as a surface.

**Classes / files.**
- `src/app/QtSasSurfaceOverlay.h/.cpp` — overlay contract.
- Consumes I5's `QtSasMesh` for geometry; `ScalarToVertex` for
  painting.

**Data path.** On `setFrame(t)`: ask `QtSasMesh` for the current
frame's mesh, ask `ScalarToVertex` to populate per-vertex
scalars from the selected field, push to the mapper.

**UI.** Toolbar toggle "Surface"; scalar dropdown reuses F1's
registry-driven selector (shared widget class).

**Diagnostics.** `h5reader.overlay.sassurface` per I8. Log
mesh vertex count, scalar min/max, time to map.

**Depends on.** I1, I5.

**Effort.** M (mesh build + interpolation + integration).

**Open questions.** Transparency: opaque vs 30% translucent
with molecule showing through. Default translucent so atoms
remain legible.

## F3 — Ribbon thickness modulation

**Goal.** Drive the existing ribbon's per-CA width from a
per-residue scalar (either a direct residue-scope field or a
per-atom field aggregated via I4).

**Unlocks.** `/dssp/hbond_energy` + derived per-residue sums
(bonded_energy total, Mahalanobis mean, CSP once wired).

**Classes / files.**
- `src/app/QtRibbonThicknessOverlay.h/.cpp` — feeds a per-
  residue width array into the existing
  `QtBackboneRibbonOverlay`.
- Widen `QtBackboneRibbonOverlay` to accept an optional width
  array and a width range `[0.3×, 2.0×]` around the default.

**Data path.** On `setFrame(t)`: fetch residue-scope field OR
compute per-residue aggregate via I4; normalise to width
range; push to ribbon filter; rerun the filter (~2-5 ms per
frame for 300 residues, already the ribbon's per-frame cost).

**UI.** "Ribbon width by…" dropdown in the ribbon section of
the toolbar.

**Diagnostics.** `h5reader.overlay.ribbon` per I8.

**Depends on.** I1, I4.

**Effort.** M.

**Open questions.** Colour and thickness as independent channels
or coupled? Independent — that's the point of the sausage
representation. Keep the existing SS-colour default.

## F4 — Superquadric tensor glyph (asymmetric shielding)

**Goal.** Render a Kindlmann superquadric at each atom for an
asymmetric rank-2 shielding tensor, shape from T2 eigenstructure,
surface colour from T0, optional T1 arrow.

**Unlocks.** `bs_shielding`, `hm_shielding`, `rs_shielding`,
`mc_shielding`, `pq_shielding`, `disp_shielding`,
`hbond_shielding`, `coulomb_shielding`, `aimnet2_shielding`,
`predictions/raw_T2`, `predictions/mlp_T2` (+ the McConnell /
Coulomb sub-decompositions once the tensor selector exposes
them).

**Classes / files.**
- `src/app/QtTensorGlyphOverlay.h/.cpp` — overlay contract,
  mode enum (Superquadric / HaeberlenAxis).
- `vtkSuperquadricSource` per atom OR a single
  `vtkTensorGlyph` with `vtkSuperquadricSource` as the input.
  VTK has `vtkTensorGlyph` that accepts symmetric tensors; use
  it for the T2 path, overlay a T1 arrow on top.
- T2 → `vtkDoubleArray(9)` per atom in Cartesian. I2 does the
  conversion.

**Data path.** On `setFrame(t)`: for every atom, read `(9,)`
SphericalTensor from `QtFrame`, convert to `Mat3` via I2,
write to a `vtkDoubleArray` with `NumberOfComponents = 9`;
hand to `vtkTensorGlyph`; colour by T0 via scalar array.

**UI.** "Tensor glyph" toolbar toggle with a "Tensor source"
dropdown (the 11 shielding fields + their decompositions).
Scale cap slider (per `feedback_surface_complex_data`, the
tensor magnitudes span orders of magnitude; global cap
default at 1.5×σ of the frame's magnitudes).

**Diagnostics.** `h5reader.overlay.tensorglyph` per I8. Log
max |T2| per frame, any atom with NaN / Inf (should never
happen per 2026-04-15 validation).

**Depends on.** I2.

**Effort.** M. The math is all in I2; the overlay is a wiring
job. Main cost is curating glyph scale across the protein and
across the 11 source fields.

**Open questions.** Render on all atoms or filter to "non-
negligible" atoms? At a full protein all-atom render is visual
noise. Default: show glyphs only on atoms where `|T0| > ε` for
the selected field (ε = protein-median / 10). Toggle for
"show all". Also: T1 rendering — small opaque arrow through
glyph, toggled off by default (T1 is typically small; surface
when an adviser asks).

## F5 — Haeberlen principal-axis glyph (traceless EFG)

**Goal.** For the traceless symmetric EFG fields, render the
Haeberlen principal-axis picture: V_zz double-ended arrow + an
oriented disk flattened by η.

**Unlocks.** `coulomb_total`, `coulomb_backbone`,
`coulomb_aromatic`, `apbs_efg`, `aimnet2_total`,
`aimnet2_backbone`, `aimnet2_aromatic`, `water/efg`,
`water/efg_first`.

**Classes / files.**
- Reuses `QtTensorGlyphOverlay` from F4 with the HaeberlenAxis
  mode selected.
- Compose from a `vtkArrowSource` (V_zz) + `vtkDiskSource`
  (η disk) per atom; or a small custom actor per atom.

**Data path.** On `setFrame(t)`: per atom, compute Haeberlen
ordering from I2, extract V_zz direction and η, build arrow +
disk actors with matching orientation.

**UI.** "EFG glyph" toolbar toggle, same overlay dropdown as F4
(with a modes sub-selection). Show η in a per-atom tooltip.

**Diagnostics.** `h5reader.overlay.haeberlen` per I8. Assert
`|Tr| < 1e-8` for all rendered atoms; warn otherwise.

**Depends on.** I2, F4 (class chassis).

**Effort.** S (extends F4).

**Open questions.** Colour mapping: by |V_zz| (magnitude) or
by η (asymmetry) or both (thickness / colour)? Default: colour
by |V_zz|, opacity by η so high-η disks appear more present.

## F6 — Vector arrow overlay

**Goal.** Per-atom arrow rendering for E-fields, B-field,
dipoles, surface normals.

**Unlocks.** `/efg/E_total`, `E_backbone`, `E_sidechain`,
`E_aromatic`, `E_solvent`, `apbs_efield`, `/water/efield`,
`efield_first`, `dipole_vector`, `surface_normal`, `/sasa/normal`.

**Classes / files.**
- `src/app/QtAtomArrowOverlay.h/.cpp` — overlay contract.
- Uses I3's `ArrowRenderer` and presets.

**Data path.** Standard per-atom Vec3 pipeline via
`vtkGlyph3D`. On `setFrame(t)`: fetch `(N,3)` slab, push to
renderer.

**UI.** "Arrows" toolbar action with dropdown over the ten
vector fields. Render only on picked atom and first-shell
neighbours by default (clutter control); toggle for "show all".

**Diagnostics.** `h5reader.overlay.arrow` per I8. Log max |v|
per field per frame.

**Depends on.** I3.

**Effort.** S.

**Open questions.** Multiple arrow fields simultaneously?
Yes — different colour presets per field so stacking is
legible (adviser wants to see E_vacuum vs E_APBS side by side).

## F7 — Relationship line segment

**Goal.** When an atom is picked, draw a thin line from the
atom to a target point (nearest-C=O midpoint, nearest H-bond
partner position) without magnitude encoding.

**Unlocks.** `/bond_aniso/dir_nearest_CO` + `nearest_CO_midpoint`,
`/hbond/nearest_dir` + partner position reconstruction.

**Classes / files.**
- `src/app/QtRelationshipOverlay.h/.cpp` — selection-gated,
  renders only on picked-atom subscription.
- Share `vtkLineSource` pattern.

**Data path.** Listen on `atomPicked(size_t)` + `frameChanged`,
read the chosen relationship field(s), draw one line per active
relationship.

**UI.** Lives in the Inspector dock as a toggle "show
relationships" with checkboxes for "nearest C=O" and "nearest
H-bond".

**Diagnostics.** `h5reader.overlay.relationship` per I8. Log
pick + target per event.

**Depends on.** F6 plumbing (line-source utility).

**Effort.** S.

**Open questions.** Extend to "nearest N pi-bond neighbours"
etc.? Punt until the inspector shows it being wanted.

## F8 — Protein-wide |B| volumetric isosurface

**Goal.** Sum the per-ring B-field evaluators over the whole
protein onto a wide grid; contour at ±0.1 mT; render alongside
the existing per-ring butterfly.

**Unlocks.** `/ring_current/total_B_field` (as a protein-wide
field view rather than per-atom arrows).

**Classes / files.**
- `src/app/QtProteinFieldOverlay.h/.cpp` — overlay contract;
  same pattern as `QtFieldGridOverlay`, wider grid, summed
  sources.
- Reuses `QtBiotSavartCalc::EvaluateBField` from the existing
  per-ring butterfly.

**Data path.** At `Build`: determine protein bounding box,
grid resolution (default 64³; UI-tunable). On `setFrame(t)`:
for every ring, evaluate B on the shared grid, accumulate into
a single `vtkImageData` scalar field; `vtkContourFilter` at
the adviser-chosen isovalues; two actors (shielded / deshielded).

**UI.** Toolbar toggle "Protein B-field" alongside the existing
"Butterfly" toggle; isovalue slider (shared control); grid
resolution slider.

**Diagnostics.** `h5reader.overlay.proteinfield` per I8. Log
grid build time, accumulated |B| statistics, ring count summed.

**Depends on.** Existing `QtBiotSavartCalc`.

**Effort.** M. Main cost: per-frame summation over 20 rings on
64³ = 5.2M voxel-evaluations. Benchmark and cache per-ring base
grids (translate + rotate into common frame) if > 300 ms/frame.

**Open questions.** Include HM field, just BS, or their sum?
Match the existing butterfly default (mode selector
BS / HM / Sum). Per-atom `total_B_field` values in the H5 are
already summed; use them to validate the volumetric sum at atom
positions ("the grid agrees with what the atoms see").

---

# Track B — Abstract docks

Order: I6 (scope framework) gates every dock; build it first,
then F9 → F11 in quick succession (straightforward per-residue
views), F10 (Ramachandran) and F12–F13 (bar charts share
pattern), then F14–F18 (more complex).

## F9 — Time-series scope expansion

**Goal.** Extend `QtAtomTimeSeriesDock` to accept per-residue,
per-ring, and system-scope scalars through the unified scope
framework (I6). Keep the atom-scope behaviour as a special case.

**Unlocks.** Residue-scoped fields (`/dihedrals/*`,
`/dssp/hbond_energy` as a chart complement to the ribbon
thickness), ring-scoped fields (`/ring_geometry/data` columns),
system-scope fields (temperature, pressure, volume, density
before F17's dedicated dock replaces them).

**Classes / files.**
- Refactor `src/app/QtAtomTimeSeriesDock` → rename to
  `QtScopedTimeSeriesDock` (keep `QtAtomTimeSeriesDock` as a
  thin alias for now to avoid breaking connections).
- Scope-picker `QComboBox` in the header (Atom / Residue /
  Ring / System). Field dropdown repopulates on scope change.
- Index wired from the appropriate picker (I6) for Atom /
  Residue / Ring; System ignores picker.

**Data path.** `ExtractFn` signature from I1 takes a
`ScopedIndex`, so the same registry drives every scope.

**UI.** Scope chips at the top of the dock. Chart axis label
updates per field.

**Diagnostics.** `h5reader.dock.timeseries` per I8.

**Depends on.** I1, I6.

**Effort.** M. Mostly a refactor of the existing dock + one
new piece (System scope).

**Open questions.** Multi-field overlay on the same chart? Pin
set of `(scope, id, index)` tuples; each gets a coloured line.
Plan this for v2 after single-scope works.

## F10 — Ramachandran dock

**Goal.** Joint-distribution scatter of (φ, ψ) per residue per
frame, current-frame highlight, optional colour by DSSP class.

**Unlocks.** `/dihedrals/phi`, `/dihedrals/psi`,
`/dihedrals/omega` (secondary trace panel below main scatter).

**Classes / files.**
- `src/app/QtRamachandranDock.h/.cpp` — owns a `QChart` with
  a `QScatterSeries` per DSSP class (8 series, toggleable via
  legend).
- Optional overlay of classical Ramachandran contours (α,
  β, L_α) as `QSplineSeries` on the same chart for reference.

**Data path.** On `Build`: load `/dihedrals/phi`, `psi` entirely
(typical (T,R) = (600, 300) = 180k points, fine for Qt Charts
at sub-sampled render). On `setFrame(t)`: rebuild the
current-frame highlight series (one point per residue at the
current frame) — fast.

**UI.** Dock widget with chart; legend; φ-range / ψ-range
defaults ±180°. Click a residue highlight → `residuePicked`
broadcast via I6 (bi-directional linkage with ribbon).

**Diagnostics.** `h5reader.dock.ramachandran` per I8.

**Depends on.** I6.

**Effort.** S.

**Open questions.** Auto-filter out GLY (distinct allowed region)
and PRO (restricted φ)? Include them with distinctive marker
shapes; it's informative, not a problem.

## F11 — DSSP residue × time heatmap

**Goal.** A 2D grid, one column per residue, one row per frame,
cell coloured by 8-class DSSP code.

**Unlocks.** `/dssp/ss8` as a time story; secondary
`/dssp/hbond_energy` view with diverging colormap.

**Classes / files.**
- `src/app/QtDsspHeatmapDock.h/.cpp` — `QGraphicsView` with a
  coloured grid; frame cursor as a horizontal line.
- Categorical colormap using the 8 DSSP-standard colours
  (α red, 3₁₀ cyan, π magenta, strand yellow, bridge tan,
  turn teal, bend green, coil white).

**Data path.** On `Build`: load `(T, R)` ss8 array, pre-render
the grid to a `QPixmap` (fast to display even at 1500 × 300).
On `setFrame(t)`: move the cursor.

**UI.** Dock widget with the grid + a colour legend +
mode-toggle ("SS class" / "hbond energy"). Click a cell →
`residuePicked` + `frameChanged` via I6.

**Diagnostics.** `h5reader.dock.dsspheat` per I8.

**Depends on.** I6.

**Effort.** S.

**Open questions.** None — this is well-established and the
data is small.

## F12 — Per-ring-type bar chart

**Goal.** For a selected atom and current frame, show eight
horizontal bars (one per RingTypeIndex) coloured by type,
length by the per-type T0 value for a selected calculator.

**Unlocks.** `/ring_current/bs_T0_per_type`,
`/ring_current/hm_T0_per_type`, `/quadrupole/pq_T0_per_type`,
`/dispersion/disp_T0_per_type`.

**Classes / files.**
- `src/app/QtPerTypeBarDock.h/.cpp` — `QChart` with
  `QHorizontalBarSeries`.
- Calculator selector in header.

**Data path.** On `atomPicked` + `setFrame(t)`: fetch
`(N, 8)` slab's row for selected atom, push to bar series.
Update ring-type colours from ring class hierarchy.

**UI.** Dock widget; calculator dropdown (BS / HM / PQ / Disp);
legend maps colour → type.

**Diagnostics.** `h5reader.dock.perype` per I8.

**Depends on.** I6 (atom picker).

**Effort.** S.

**Open questions.** Include a "weighted by ridge coefficient"
mode so the bars reflect predicted contribution rather than raw
kernel? Defer until predictions are confirmed reliable; raw
first.

## F13 — K=6 nearest-ring bar chart

**Goal.** For a selected atom, six horizontal bars, one per
nearest ring, ordered by distance, coloured by ring type, length
by the selected per-ring quantity.

**Unlocks.** `/per_ring/bs_T2`, `hm_T2`, `chi_T2`, `pq_T2`,
`hm_H_T2`, `disp_scalar`, `disp_T2`.

**Classes / files.**
- `src/app/QtPerRingBarDock.h/.cpp` — same `QHorizontalBarSeries`
  pattern as F12.
- Selector for field: `(calculator × scalar-or-|T2|)`.

**Data path.** On `atomPicked` + `setFrame(t)`: fetch
`(N, K, ...)` slab, order by `/per_ring/geometry` distance
column (from `geometry_fields`), push to bars.

**UI.** Dock widget; field dropdown. Hover a bar → highlights
the ring in the 3D scene and opens a tooltip with full per-ring
geometry.

**Diagnostics.** `h5reader.dock.perring` per I8.

**Depends on.** I6 (atom picker), `geometry_fields` interpretation.

**Effort.** M (the hover → highlight wiring is non-trivial).

**Open questions.** Show ring k index or geometric distance as
the y-axis tick? Show both (`"#0 · 3.2 Å · PHE1"`).

## F14 — Predictions stacked-bar decomposition

**Goal.** For a selected atom and current frame, render the
additive decomposition of the ridge / MLP prediction: six core
groups (ring_current, efg, bond_aniso, quadrupole, dispersion,
hbond) + three new-feature groups (water, charges, sasa),
stacked positive / negative, for T0 or |T2|.

**Unlocks.** `/predictions/raw_{ring_current, efg, bond_aniso,
quadrupole, dispersion, hbond, water, charges, sasa}` and the
`norm_*` counterparts; `/predictions/mlp_T0`, `mlp_T2`.

**Classes / files.**
- `src/app/QtPredictionsDock.h/.cpp` — stacked
  `QHorizontalBarSeries` with per-group colour.
- Selector for: raw vs norm; T0 vs |T2|; ridge vs MLP.

**Data path.** On `atomPicked` + `setFrame(t)`: read each group
slab's row for selected atom, compute |T2| if T2 mode, push to
stacked bars. Annotate total = sum.

**UI.** Dock widget with mode toggles. Per-group colour legend.
A second row "predicted total" + "measured sum" pair of bars so
reviewers see the delta (that's the ridge residual per atom
per frame).

**Diagnostics.** `h5reader.dock.predictions` per I8.

**Depends on.** I6.

**Effort.** M.

**Open questions.** Bar orientation (horizontal stacked is
easier for labels, vertical is easier for grouping multiple
atoms side-by-side). Start horizontal; switch if multi-atom
view is added.

## F15 — MLP 55-kernel heat strip

**Goal.** Horizontal heat-strip of `/predictions/mlp_kernel_weights`
for the selected atom, one cell per kernel, coloured by gate
weight. Kernel labels from the kernel manifest.

**Unlocks.** `/predictions/mlp_kernel_weights`,
`/projections/mlp_gate_weights` (same data path).

**Classes / files.**
- `src/app/QtKernelWeightsDock.h/.cpp` — custom widget with a
  `QGraphicsView` grid (55 cells × 1 row) + labels.
- Kernel manifest loader (JSON alongside the model artifacts);
  falls back to `"kernel_0"…"kernel_54"` if not present.

**Data path.** On `atomPicked` + `setFrame(t)`: fetch `(N, 55)`
slab's row, colour each cell via sequential viridis.

**UI.** Dock widget. Tooltip on each kernel cell pulls from I7
glossary (kernel name, what it means, typical value range).

**Diagnostics.** `h5reader.dock.kernelweights` per I8.

**Depends on.** I6 (atom picker), I7 (kernel descriptions).

**Effort.** S.

**Open questions.** Include a protein-wide reference row ("mean
across all atoms at this frame") for comparison? Yes — second
row below the per-atom row, same colormap.

## F16 — AIMNet2 embedding projection (UMAP / PCA)

**Goal.** 2D projection of the 256-D AIMNet2 embedding across
(atom × frame), coloured by element / residue type / DSSP, with
the currently picked atom's trajectory traced.

**Unlocks.** `/aimnet2_embedding/aim`; also `/projections/
subspace_coords` as a second mode.

**Classes / files.**
- `src/app/QtEmbeddingProjectionDock.h/.cpp` — `QChartView`
  scatter.
- UMAP or PCA reducer:
  - Minimal: principal components via eigen on the
    covariance (fast, deterministic, good enough for a
    reference tool). Use I2's eigendecomposition.
  - Ambitious: UMAP via `umappp` (header-only C++ UMAP). Slower
    to build (~10 s), more faithful for non-linear structure.
  - Start with PCA; add UMAP as a second mode once PCA is wired.

**Data path.** On file-load, sample ~10⁵ points from (atom ×
frame), fit PCA (or UMAP), project all ~10⁵ at display time and
the currently-picked atom's trajectory on top. Cache projection
coefficients so `setFrame` only moves the highlight.

**UI.** Dock; element-colour toggle; projection-mode selector
(PCA axes 1-2, 2-3, etc.; UMAP). Click a point → `atomPicked`
+ `frameChanged` via I6.

**Diagnostics.** `h5reader.dock.embedding` per I8. Log fit time,
variance-explained per PCA axis.

**Depends on.** I6.

**Effort.** L (UMAP integration is the long pole).

**Open questions.** Element-pooled vs per-element projection?
Per-element is what the Stage 1 PCA uses; do per-element as the
default (separate scatters per H / C / N / O). Share axes
across elements if an adviser asks.

## F17 — System-energy dock

**Goal.** Multi-panel line chart of the ~20 scalar `/energy/*`
time series. Frame cursor shared with the main playback.

**Unlocks.** `/energy/potential`, `kinetic`, `enthalpy`,
`temperature`, `pressure`, `volume`, `density`, `coulomb_sr`,
`coulomb_recip`, `bond`, `angle`, `urey_bradley`, `proper_dih`,
`improper_dih`, `cmap_dih`, `lj_sr`, `T_protein`,
`T_non_protein`, `box.x`, `box.y`, `box.z`.

**Classes / files.**
- `src/app/QtSystemEnergyDock.h/.cpp` — multiple
  `QChart` stacked vertically in a single `QScrollArea`.
- Group the panels: Thermodynamic (T, P, V, ρ), Coulomb
  (sr+recip), Bonded (6 terms), LJ, Box (3 edges), Enthalpy.

**Data path.** On `Build`: load every `(T,)` dataset in
`/energy`. On `setFrame(t)`: move the vertical cursor on every
panel.

**UI.** Dock widget; panel visibility toggles in the header so
reviewers can focus on a subset.

**Diagnostics.** `h5reader.dock.energy` per I8. On `Build`,
validate `T_protein` + `T_non_protein` matches global temperature
(within 1 K); warn on drift.

**Depends on.** I6 (System scope).

**Effort.** M.

**Open questions.** Y-axis sharing across panels? No — each
panel autoscales; reviewers want visible variation per term.

## F18 — 3×3 tensor matrix view (virial, pressure)

**Goal.** Render the 3×3 `virial` and `pressure_tensor` as a
matrix view with a line chart per diagonal element and a
secondary panel for off-diagonal entries (equilibrated systems
have ≈0 off-diagonal).

**Unlocks.** `/energy/virial`, `/energy/pressure_tensor`.

**Classes / files.**
- `src/app/QtMatrixTensorDock.h/.cpp` — `QGridLayout` of 9
  mini-charts + a cursor line.

**Data path.** On `Build`: load `(T, 9)` arrays, split into 9
`(T,)` series per the `virial_layout` attribute. On
`setFrame(t)`: cursor.

**UI.** Dock widget.

**Diagnostics.** `h5reader.dock.matrix` per I8. Warn on
sustained non-zero off-diagonal (> 10 bar mean).

**Depends on.** I6 (System scope).

**Effort.** S.

**Open questions.** Include virial separately from
pressure_tensor or tabbed? Two tabs in one dock.

---

# Track C — Inspector + tooltips

## F19 — Inspector small-multiples for T2 stacks

**Goal.** In the atom inspector, for any atom with a
per-category tensor decomposition (McConnell backbone /
sidechain / aromatic / CO_nearest / CN_nearest; Coulomb total /
backbone / aromatic; per-ring-type stacks for BS / HM / PQ /
Disp), render tiny superquadric glyphs side-by-side under the
group node so the adviser sees the decomposition at a glance.

**Unlocks.** Already-loaded T2 per-category / per-type fields
become interpretable in one view instead of requiring an overlay
toggle.

**Classes / files.**
- `src/app/QtMiniGlyphWidget.h/.cpp` — thumbnail-sized VTK
  render widget (or CPU-rasterised superquadric preview for
  performance).
- Extend `QtAtomInspectorDock` to host these widgets as custom
  tree-item payloads on the T2-decomposition group nodes.

**Data path.** On `setPickedAtom` / `setFrame`: inspector asks
I2 for each sub-tensor's glyph parameters, hands them to
`QtMiniGlyphWidget`. Pure display.

**UI.** Inside the existing inspector tree. No new top-level
docks.

**Diagnostics.** `h5reader.inspector.miniglyph` per I8.

**Depends on.** I2, F4 glyph machinery.

**Effort.** M.

**Open questions.** CPU-rasterised preview vs live VTK widget
per mini-glyph: at 10+ mini-glyphs visible, per-widget VTK
actors can be slow. Prototype with VTK, fall back to a raster
if profiling warrants.

## F20 — Glossary Navigator

**Goal.** The glossary as a first-class navigational surface,
not a passive help file. Advisers open one dock, browse every
H5 field, search / filter, see the field's current value at the
current atom and frame, and jump with one click to the overlay
or dock that surfaces it. The glossary becomes the primary
index into the data — a reviewer who does not know the codebase
can open the reader cold and use the navigator as their starting
point.

**Unlocks.** Every H5 field becomes discoverable without
docs. Right-click "show in glossary" from any dropdown /
inspector row / chart title closes the loop in the other
direction (UI element → glossary entry). Tooltips on the
dropdowns consume the same lookup.

**Classes / files.**
- `src/app/QtGlossaryDock.h/.cpp` — the dock. Three panes:
  - **Schema tree** (left): H5 groups as folders, fields as
    leaves. Same structure as the glossary markdown.
  - **Search / filter bar** (top): live filter by name, scope,
    units, modality, "has overlay", "has dock". Incremental;
    empty filter shows full tree.
  - **Detail panel** (right): for the selected field:
    - Name, H5 path, shape, dtype, units, scope.
    - One-paragraph physical meaning (from the glossary).
    - Provenance (calculator + file:line).
    - **Primary modality** + **useful modalities** as chips.
    - **Current value** at current atom × current frame (where
      defined), with units; updates on frame advance and on
      atom / residue / ring pick via I6.
    - **Show me** buttons — one per modality the field has.
      Click "Show on atom" → enables the colour-bubble overlay
      bound to this field. Click "Plot" → adds this field to
      the scope-expanded time-series dock. Click "Glyph" → turns
      on the tensor-glyph overlay with this field selected.
      Click "Paint surface" → SAS surface with this scalar.
      The catalogue of buttons comes from the field's modality
      assignments in the glossary.
    - **Related fields** — cross-links, e.g. `bs_shielding` →
      `bs_T0_per_type`, `bs_T2_per_type`; `coulomb_total` →
      `coulomb_backbone`, `coulomb_aromatic`; every shielding
      field → `predictions/raw_T2`.
- `src/model/GlossaryNavigatorModel.h/.cpp` —
  `QAbstractItemModel` wrapping the `FieldGlossary` JSON from
  I7 so the tree and the search filter are both cheap.
- `src/app/ShowMeRouter.h/.cpp` — takes a `(field_id, modality)`
  pair, dispatches to the right overlay or dock. Singleton-ish,
  registered by overlays / docks at startup so the router knows
  what's available.
- Reverse wiring: every dropdown entry, every inspector row,
  every chart title gets a context menu "Show in glossary" that
  calls `QtGlossaryDock::openField(id)`.

**Data path.** On file-load: `FieldGlossary` (from I7) provides
the entire field catalogue; `GlossaryNavigatorModel` wraps it;
the dock renders the tree. On `frameChanged` / atom pick / etc.:
detail panel refreshes the current-value display only for the
currently-selected field (O(1) per update). "Show me" buttons
call into the `ShowMeRouter` which in turn drives overlay /
dock settings.

**UI.** Default-docked on the right, tab-grouped with the atom
inspector. Keyboard shortcut **Ctrl+G** opens / raises it.
Toolbar "Glossary" button. The schema tree stays visible; the
detail panel is the interactive part. No modal dialogs —
navigation is one panel.

**Diagnostics.** `h5reader.dock.glossary` per I8. Log every
"Show me" dispatch with target overlay / dock so an adviser can
replay a session from the UDP stream. Log missing-description
lookups (should be zero; CI assertion that every
`ScalarRegistry` id has a `FieldGlossary` entry).

**Depends on.** I7 (glossary JSON), I1 (scalar registry for
current values and "Show me" routing), I6 (scope framework
for atom / residue / ring / frame-scoped current values).

**Effort.** M. The model + dock is one session; the ShowMeRouter
+ reverse-wiring ("show in glossary" context menus on every
bound widget) is a second. Worth the time — this is the feature
that lets the reader be self-explanatory to a professor opening
it cold.

**Open questions.** The glossary markdown is the authority;
changes there should propagate at build time via the extractor
script. What happens when a field is in the H5 but not in the
glossary? Schema tree shows it with a warning badge; detail
panel shows "no description — see `H5_FIELD_GLOSSARY.md`".
Drives contributors to keep the glossary complete.

**Tooltip discipline (incorporated, not a separate feature).**
Every widget bound to a field — dropdown entry, inspector row,
chart title, bar-segment label — calls
`setToolTip(FieldGlossary::tooltip(id))` at construction. The
first-line summary from the glossary entry is the tooltip. No
bespoke tooltip strings; the glossary is the single source. This
is part of every overlay / dock session (F1 onward) and audited
in R5.

---

# Refinement cadence (R1–R5)

Feature sessions add surface area; refinement sessions tighten
it. **One dedicated session after every 2–4 feature sessions.**
Without this cadence the app accumulates inconsistent chrome,
ad-hoc colour choices, and half-wired tooltips — then shipping
to an adviser requires a big-bang polish pass that takes longer
than the refinements would have.

## What a refinement session does

Not new features. Each refinement pass walks a checklist across
every surface touched since the last refinement. Scope-limited:
one session, one checklist.

**Chrome and layout.**
- Toolbar: icon style consistent, ordering logical (view /
  overlays / docks / help), separators in the right places.
- Dock defaults: which side, which tab group, which visible at
  startup. An adviser who opens the reader cold sees a
  sensible layout, not every dock stacked.
- Menus: right-click menus consistent across widgets. "Show in
  glossary" available on every field-bound widget.
- Keyboard accelerators: Ctrl+G glossary, Space play/pause,
  ←/→ step frame, F toggle frame-time display, etc. Documented
  in help.

**Colour and typography.**
- Colormaps: diverging Moreland for signed fields, sequential
  viridis for non-negative, categorical for enums (DSSP,
  BondCategory). Every colour-mapped overlay has a
  `vtkScalarBarActor` or `QtLegend` widget. No bare colours.
- Vector presets (I3): backbone grey, sidechain tan, aromatic
  purple, solvent cyan, water blue, B-field gold — one
  source, one colour, app-wide.
- Font hierarchy: dock titles, chart titles, axis labels, tick
  labels. Qt defaults are fine; check proportions only.

**Defaults and empty states.**
- Every overlay's "off by default" unless the user enabled it
  last session (persist in `QSettings`).
- Every dock shows a meaningful empty state when no atom /
  residue / ring picked (not a blank chart with "No data").
- Every dropdown has a sensible first entry.

**Tooltips and field descriptions.**
- Every field-bound widget calls `FieldGlossary::tooltip(id)`
  at construction. Audit: grep the codebase, find any widget
  that binds to a field_id without a tooltip call.
- Context menu "Show in glossary" on every such widget.

**Error and loading states.**
- HighFive seam: every `try / catch` logs specific field +
  frame + path at boundary failure, UI shows a non-intrusive
  "couldn't load `/path/to/field`" banner, the rest of the app
  keeps working with partial state.
- Long-running operations (SAS mesh build, UMAP fit): progress
  indicator, cancellation where possible, non-blocking UI.

**Frame-transition smoothness.**
- No flicker when toggling an overlay.
- Frame-advance under the target (50 ms for most overlays,
  100 ms with tensor glyph + volumetric isosurface both on).
- Measure. `QElapsedTimer` around `setFrame` at DEBUG level on
  UDP. Benchmark on 2000-atom fixture.

**Cross-platform.**
- Linux is the development platform; Mac and Windows appearance
  checks are deferred per `project_h5reader_target_hardware`.
  Each refinement session notes any Mac / Windows-affecting
  changes so the platform sweep has a clean list.

**UDP log audit.**
- Every overlay / dock emits its `h5reader.<scope>.<name>`
  category per I8.
- No per-render-call log spam at INFO; move to DEBUG.
- Errors at WARNING / ERROR with actionable context.

## R1–R5 scope

- **R1** (post-infra). No feature surface to polish yet; the
  session audits the chrome: toolbar skeleton, menu structure,
  default dock placement, keyboard map, `QSettings` persistence
  keys. Establishes the base look every feature inherits.
- **R2** (post F1–F3, F4–F5). Scalar and tensor families. Audit:
  colour-map consistency across colour bubble / SAS surface /
  ribbon / glyphs; legend positioning; scale caps; tooltip
  completeness on the ~50 fields these surface.
- **R3** (post F6–F7, F8, F20). Vectors, volumetric, glossary
  navigator. Audit: vector-preset palette, volumetric isovalue
  slider ergonomics, glossary "Show me" routing to every
  overlay + dock added so far, "Show in glossary" context menus
  on every field-bound widget.
- **R4** (post F9–F12, F18). Simple docks. Audit: dock default
  positions, scope-picker consistency across docks,
  frame-cursor synchronisation, chart styling uniformity.
- **R5** (final, includes F19 + tooltip audit). Inspector
  small-multiples land here. Final tooltip pass over every
  bound widget. Performance benchmark on the 2000-atom
  fixture. Prep for cross-platform pass (which runs in a
  separate session from a Mac / Windows dev environment).

---

# Sequencing and dependencies

Session order, interleaved with refinements.

1. **Infra bundle 1** (parallelisable): I1, I2, I3, I4, I6, I8.
   3–4 sessions.
2. **Infra bundle 2**: I5, I7. 2 sessions.
3. **R1** — post-infra chrome. 1 session.
4. **Scalar family**: F1 → F2 → F3. 2–3 sessions.
5. **Tensor family**: F4 → F5. 2 sessions.
6. **R2** — post-scalar + tensor. 1 session.
7. **Vector family**: F6 → F7. 1–2 sessions.
8. **Volumetric**: F8. 2 sessions.
9. **Glossary Navigator**: F20. 2 sessions.
10. **R3** — post-navigator. 1 session.
11. **Dock family — simple**: F9, F10, F11, F12, F18.
    3–4 sessions (F9 refactors the existing time-series dock;
    F10–F12, F18 share the `QChart` bar pattern).
12. **R4** — post-simple docks. 1 session.
13. **Dock family — complex**: F13, F14, F15, F17. 4 sessions.
14. **R5** — final polish + F19 inspector small-multiples +
    tooltip audit + 2000-atom benchmark. 2 sessions.
15. **Optional**: F16 AIMNet2 embedding projection. 3 sessions
    (defer past v1 if the PCA-only fallback in I2 is enough
    for the first adviser pass).

**Total with refinements.** ~24 sessions to the v1 feature set,
~27 with F16 ambitious UMAP. Refinement sessions are ~20% of
the budget — appropriate for a reference implementation that
professors will run cold.

**Minimum viable adviser demo.** I1 + I2 + I6 + F1 + F3 + F4 +
F10 + F11 + F20 + R1 + R2 + R3. Twelve sessions. Gets colour
bubbles, sausage ribbon, tensor glyphs, Ramachandran, DSSP
heatmap, glossary navigator with tooltips, and three rounds of
polish — a first cut that already looks like a thesis tool and
already explains itself to a reviewer.

---

# Cross-cutting discipline

Every feature session observes these, not because they're
optional polish but because the existing viewer rebuild
(three-restart cycle in `ui/`) proved they are load-bearing.

- **Invoke the `qt6-cpp` skill at session start** — the Qt /
  VTK patterns vary enough that codifying the pattern before
  writing code is always the faster path.
- **`CENSUS_REGISTER` every `QObject`, `ACONNECT` every signal
  / slot, `ASSERT_THREAD` every GUI-affinitised method** per
  `feedback_qt_discipline`.
- **UDP category per overlay / dock** per I8.
- **HighFive boundary is always wrapped in try / catch** with
  a specific-context log and graceful degradation. See
  `feedback_qt_citizen` on continuing with partial state.
- **File I/O via `QFile` / `QDir` / `QStandardPaths`** —
  cross-platform bar per `project_h5_reader_scope`. No POSIX
  shortcuts.
- **No `QTimer` on interactive controls** except
  `QtPlaybackController`. For animation, listen on its
  `frameChanged`.
- **No kernel re-evaluation that the H5 already answered.** The
  F8 volumetric B-field is the explicit exception (the H5 has
  per-atom `total_B_field` but not the protein-wide volumetric
  grid; re-eval uses the same evaluator as the existing
  butterfly). Every other overlay reads H5, transforms,
  renders — never computes physics.
- **Field names on screen match H5 paths** per
  `feedback_ui_data`. The glossary's "surfacing discipline"
  rules are contractual.

Performance targets (informational, refine per
machine once benchmarked on 2000-atom × 600-frame fixture):

- Frame advance: end-to-end cost ≤ 50 ms (20 fps feels
  responsive for scrubbing). With F4 + F8 both on, ≤ 100 ms
  is acceptable.
- File load: ≤ 15 s for a 2 GB file. SAS mesh pre-compute
  (I5) within this budget.
- Dock startup: ≤ 2 s after file-load.
- Memory: peak ≤ 12 GB on a 5090; embeddings dominate and the
  PCA fit holds.
