# UI Roadmap: Protein Tensor Viewer

**Written:** 2026-04-06, session 5 (MOPAC UI integration).

## What the viewer IS

A verification and exploration tool for the library's geometric kernel
computations. Every value displayed comes from the library objects
directly — the viewer is a window, not a calculator. When you click an
atom, you see what the library computed. When you toggle a visualization,
you see the spatial structure of a specific calculation. Nothing is
recomputed or re-derived in the viewer.

## What the viewer is NOT

- A prediction tool (no heuristic tiers, no "how close to DFT" overlays)
- An adapter layer (no ViewerResults, no intermediate structs)
- A second source of truth (if the viewer shows a bond midpoint,
  it reads `conf.bond_midpoints[i]`, it does not compute 0.5*(a+b))

## Current state (session 5)

**Working:**
- Ball-and-stick rendering with double bonds from topology
- Atom inspector: double-click → full object model including all 10
  calculator contributions, MOPAC electronic data (charges, orbital
  populations, valency, bond neighbours with Wiberg orders), vector
  fields (ff14SB and MOPAC Coulomb side by side), McConnell breakdowns
  (unweighted and bond-order-weighted)
- Bond inspector: shift+double-click → bond identity, MOPAC Wiberg
  order, endpoint charges, McConnell contribution
- Bond order color overlay: tubes colored by Wiberg order on a
  blue→red colormap (single→double)
- Operations log panel: UDP listener streaming library log in real time
- REST server with port fallback (tries 9147-9156)
- REST responds immediately during computation (no blocking spin loop)

**Picking:**
- Screen-space projection picking (WorldToDisplay for every atom/bond,
  nearest in 2D with depth tiebreaker). No VTK pickers — they don't
  work reliably with vtkMoleculeMapper's impostor rendering.

**Known issues:**
- Overlay modes (Heuristic, Classical, DFT Delta, Residual) and the
  8 physics checkboxes are from the prediction-focused era. They
  predate the mature kernel catalogue and MOPAC integration. They
  work but are not the right abstraction.
- Eigenvector/tensor glyph display is broken
- Isosurface overlay needs rethinking
- Crash on reload still present (command-line only loading)

## Vision: calculator-centric visualization

The sidebar reorganizes around calculators, not prediction modes.
Each calculator gets a section with 2-3 meaningful visualizations
of what it computed. The atom and bond inspectors show everything
per-object. The 3D view shows the spatial structure of each
calculation.

### Per-calculator visualizations (2-3 each)

**Ring current (Biot-Savart, Haigh-Mallion):**
1. Conformal isosurface butterfly — the B-field/shielding field
   sampled on a 3D grid around each ring, rendered as isosurfaces.
   Already half-built (ComputeWorker Phase 3/4 grid sampling).
   Should be a first-class display, not buried under "heuristic."
2. Ring outlines with current direction arrows (existing overlay,
   needs polish)
3. Per-atom G tensor eigenvectors at REPORT-tier atoms near rings

**McConnell bond anisotropy (unweighted + MOPAC-weighted):**
1. Tensor glyphs at atoms showing dipolar angular pattern
2. Bond coloring by McConnell scalar (which bonds dominate shielding)
3. Comparison mode: unweighted vs bond-order-weighted side by side
   to see where MOPAC changes the angular structure

**Coulomb EFG (ff14SB + MOPAC charges):**
1. E-field arrows at atoms (direction and magnitude)
2. EFG tensor glyphs (the traceless symmetric T2 pattern)
3. Charge delta surface: color atoms by (MOPAC - ff14SB) charge
   to show where conformation shifts electrostatics

**Pi-Quadrupole:**
1. EFG isosurface around each ring (1/r^5 decay pattern)
2. Per-atom scalar colored by quadrupole contribution

**Ring Susceptibility:**
1. Tensor glyphs (same kernel as McConnell but from ring center)

**London Dispersion:**
1. Contact surface: highlight atom-ring vertex pairs within cutoff

**H-Bond:**
1. H-bond vectors with strength coloring
2. Dipolar tensor glyphs at donor/acceptor atoms

**MOPAC electronic structure (not a calculator, but first-class data):**
1. Bond order coloring (implemented this session)
2. Charge comparison surface (MOPAC vs ff14SB delta)
3. Orbital population display (s/p ratio as hybridization indicator)

### CalculationAreas (future)

When CalculationArea objects land on ProteinConformation, each one
describes a spatial region where a specific calculation was applied
with its cutoff. The viewer displays these as:

1. Region list panel — all CalculationAreas with metadata
2. Click a region → highlight constituent atoms, show what
   calculation was applied and its parameters
3. 3D wireframe/convex hull for each region, color-coded by
   calculator type
4. Toggle per-calculator to show/hide its regions

These are self-documenting objects: the conformation knows where
each calculation reached and what it did. The viewer just shows them.

### Sidebar structure (target)

Replace the current overlay mode + physics checkboxes with:

```
[Rendering]
  Ball & Stick / VDW / Liquorice

[Ring Currents]
  □ Ring outlines
  □ BS butterfly isosurface
  □ HM isosurface
  Iso threshold: [___] A

[Bond Anisotropy]
  □ McConnell tensor glyphs
  □ MOPAC McConnell glyphs
  □ Bond McConnell scalar coloring

[Electric Fields]
  □ Coulomb E-field arrows
  □ MOPAC Coulomb E-field arrows
  □ EFG tensor glyphs
  □ Charge delta surface

[MOPAC Electronic]
  □ Bond order coloring
  □ Charge delta (PM7 - ff14SB)

[Regions] (when CalculationAreas arrive)
  □ Show calculation regions
  Per-calculator toggles
```

Each checkbox toggles ONE visualization of ONE calculator's output.
No modal overlay switching. No physics checkboxes that modify a
shared overlay. Each visualization reads directly from the library
and renders independently.

### Tensor display (ongoing study)

Good ways to show per-atom symmetric traceless tensors (T2, 5
components) are an active area. Options being explored:

- Ellipsoid glyphs (eigenvalue decomposition → axis lengths)
- Spherical harmonic surfaces (direct T2 → Y_2^m rendering)
- Superquadric glyphs
- Color-mapped eigenvector arrows

The eigenvector display was working but is currently broken. Fixing
it is a prerequisite for the per-calculator tensor visualizations.

## Implementation order

1. Fix picking (this session — screen-space projection)
2. Fix tensor glyph display
3. BS/HM butterfly as first-class isosurface
4. McConnell tensor glyphs with per-category coloring
5. E-field arrows for Coulomb and MOPAC Coulomb
6. Sidebar reorganization (calculators, not prediction modes)
7. CalculationArea display (when they land in the library)
