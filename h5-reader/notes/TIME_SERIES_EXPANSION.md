# Time-series illustrator — scope for the next session

**Context.** Step 7 landed the scaffold: `QtAtomTimeSeriesDock` with a
QComboBox of ~35 per-atom per-frame scalars plotted against frame
index via Qt6 Charts, with a vertical cursor tracking the current
frame. Tabified into the right dock area alongside the atom inspector.

The reader's reason for existing is to let a chemist or adviser read
protein dynamics through these time series. What's here is the smallest
possible thing; the real feature is bigger.

## Scope the next session should take

Jessica's direction (this session, 2026-04-17):

> When we wrap, which looks to be soon, the next session will be
> responsible for looking at whether time series data for an atom is
> glyphs, colour bubbles, strip charts or all three, and what values we
> have in the h5 and also which apply protein or residue or ring-wide.
> This will be our primary illustrator for time series data in protein
> context and we have the basics to do something really nice here. It
> is 100s of fields so next session work, but that is why we need this.

Four design axes for the next session:

### 1. Representation: strip chart vs glyph vs colour bubble vs all three

The current dock is a single strip chart. Real use cases probably want
several simultaneous representations:

- **Strip chart** (current): one scalar vs frame, cursor at current
  frame. Best for trend / variability.
- **Glyph on the atom**: per-frame glyph at the atom position in the 3D
  view, sized/oriented by a scalar. Best for spatial context.
- **Colour bubble**: atom rendered with size/colour modulated by the
  scalar at the current frame. Best for whole-protein scan at a glance.
- **Residue strip / ring strip**: one strip per residue (or ring)
  stacked in a residue-indexed axis. Best for "where along the chain
  is this signal active?"

Probably all three, with UI for choosing which to show.

### 2. Data scope: atom / residue / ring / protein / system

The H5 carries per-atom, per-residue, per-ring, and per-frame data:

- Per-atom, per-frame: ~30 scalars + ~20 vectors + ~15 SphericalTensors
  surfaced in step 6/7. Many more exist (per-type T0 per ring, K=6
  nearest rings, AIMNet2 embedding 256-D, etc.) — a hundred-plus fields
  total.
- Per-residue, per-frame: `dihedrals/{phi, psi, omega, chi1..4, chi*_{cos,sin}}`,
  `dssp/{ss8, hbond_energy}`, per-atom fields broadcast (any summary
  the extractor records per-residue).
- Per-ring, per-frame: `ring_geometry/data[T, n_rings, 7]`, `per_ring/*`
  (K=6 nearest rings per atom, with bs/hm/chi/pq/hm_H/disp T2).
- Per-frame system-level: `energy/*` (42 EDR terms: temperature,
  pressure, virial, box, LJ, Coulomb SR/recip, bonded decomposition).

That's multiplicity × hundreds of fields. The UI needs to let a user
pick the scope first, then a field, then a representation.

### 3. Cross-atom comparison

A single atom is limiting. The next session should probably support:

- Pin multiple atoms and overlay their time series on one chart.
- Pick a residue, show all its atoms as stacked strips.
- Pick a ring, show per-ring-vertex strips.

### 4. Aggregations

- Protein-wide mean / min / max of a scalar per frame.
- Standard deviation / histogram across atoms per frame.
- Rolling average (smoothing noisy trajectories).

## Inventory: what's in the H5

Check `fileformat/analysis_file.h` for the full 186-dataset schema. The
~35 scalars we surface in step 7 are a curated minimum; any field in
`analysis_file.h`'s struct can be added to the `descs()` table in
`QtAtomTimeSeriesDock.cpp` with one extra entry. That file's `Desc`
table is deliberately data-driven so the next session can expand it
mechanically.

For non-atom-scoped data (per-residue DSSP, per-ring geometry,
per-frame energies) the `Desc` table needs extending to carry a
`scope` tag plus a different index type (ring index, residue index,
etc.) — that's the first architectural change.

## What's in place right now

- `QtFrame` per-atom slab accessors (bsShielding, hmShielding, mcShielding,
  coulombShielding, apbsEfg, etc.) — see `model/QtFrame.h`.
- `QtFrame::ringGeometry(ringIdx)` + `ringVertices(ringIdx)` — ring-scoped.
- `QtFrame::dsspCode(residueIdx)` — residue-scoped. One sample among
  many we haven't yet exposed.
- Atom picker + selection overlay + atom inspector dock + time-series
  dock, all driven by one `QtAtomPicker::atomPicked(size_t)` signal +
  `QtPlaybackController::frameChanged(int)`.
- Shared `Desc` pattern in `QtAtomTimeSeriesDock.cpp` — extend the
  table for atom-scoped scalars; duplicate-and-adapt for residue/ring
  scopes.

Next session starts with the representation matrix and the scope-
selector UI, then expands the descriptor tables. Don't redesign what
landed in step 7; extend it.
