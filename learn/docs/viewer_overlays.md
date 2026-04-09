# Viewer Overlay Ideas for Secondary Analysis

The Qt/VTK viewer renders molecules with conformal tensor butterflies
and ring geometry.  These overlays add secondary analysis results as
translucent per-atom shells or halos, not atom recolouring.  Each
overlay answers a specific question about the calibration physics.

Data source: `output/secondary/viewer/{protein_id}.json`, exported
by `python -m secondary viewer_export --protein {ID}`.  The JSON
includes reality-check fields (element, residue_index, position) so
the viewer can verify referential integrity before applying overlays.

---

## 1. Polarisation disagreement halo

**Question:** Where does the quantum charge model see something
different from the force field?

**Data field:** `divergence_cosine` (1.0 = agree, 0 = orthogonal,
-1 = opposite).

**Overlay:** Translucent sphere at each atom.  Radius proportional
to `divergence_mag` (how big the difference is).  Colour mapped from
cosine: blue (agree, >0.8) → transparent (neutral, ~0.5) → red
(disagree, <0.2).  Most atoms are faintly blue or invisible.  The
few red halos show where polarisation changes the angular pattern.

**Thesis value:** On a protein with HIE, the red halos should
cluster near the imidazole — showing that the electrostatic
environment near histidine is where force field and QM diverge most.
On a protein with only PHE/TYR, the red halos should be sparse and
distant.

---

## 2. Ring-type zone overlay

**Question:** Which ring type dominates each atom's environment?

**Data field:** `stratum` (hie_only, phe_only, tyr_only, trp_only,
no_hie) and `nearest_ring_type`.

**Overlay:** Translucent shell coloured by stratum, using the
standard ring-type palette from common.R (orange = PHE, blue = TYR,
green = TRP, grey = HIE).  Radius scaled by 1/ring_dist so nearby
atoms have larger shells.  The ring itself is rendered as a filled
translucent disc with the same colour.

**Thesis value:** Shows the spatial "territory" of each ring type.
For a TRP residue, you see a large green zone extending 8-10A.  For
HIE, a smaller grey zone.  Overlapping zones (atoms that see
multiple ring types) get blended colours.  This is the spatial
version of the strata_atom_counts bar chart.

---

## 3. Residual magnitude shell

**Question:** Where on the molecule does the model fail?

**Data fields:** `self_fit_residual` (4 ring-type kernels only) and
`ridge_residual_mag` (all 91 kernels).

**Overlay — self-fit:** Translucent shell with radius proportional
to self_fit_residual magnitude.  Colour: uniform warm tone (amber).
Large shells = the ring current model doesn't explain this atom.
Small or absent = the model works here.

**Overlay — ridge:** Same idea, cooler tone (teal).  The difference
between the two overlays shows what the non-ring kernels (EFG, bond
anisotropy) rescue.  Atoms with large amber but small teal shells
are where EFG picks up what ring current misses.

**Thesis value:** Near HIE, almost every atom has a large amber
shell (self-fit R² = 0.062) but many have smaller teal shells
(ridge R²_all = 0.24 for hie_only).  The viewer shows the spatial
pattern of this rescue — likely concentrated on atoms along the
charge redistribution axis of the imidazole.

---

## 4. Target magnitude halo

**Question:** Where is the DFT signal we're trying to explain?

**Data field:** `target_t2_mag` (DFT WT-ALA delta T2 magnitude
in ppm).

**Overlay:** Translucent shell, radius proportional to T2 magnitude.
Colour: purple gradient (faint below 1 ppm, vivid above 2 ppm).
This is the ground truth — what the aromatic residue does to the
shielding tensor angular structure at each atom.

**Thesis value:** The target magnitude map shows the "reach" of each
ring type.  TRP produces vivid halos out to 8A+.  HIE produces
fainter halos, concentrated close to the ring.  Combined with
overlay 3, you see: large vivid purple (strong DFT signal) with
large amber overlay (model can't explain it) = the hard atoms.

---

## 5. Combined: the diagnostic stack

For a single protein (e.g., one with both TRP and HIE):

- Layer 1 (back): Ring-type zone (coloured translucent regions)
- Layer 2: Target magnitude (purple halos showing where DFT signal is)
- Layer 3 (front): Residual shells (amber/teal showing where model fails)

With VTK opacity control, the viewer can fade layers in/out.  A
poster figure would show the same protein three ways:
panel A = zone overlay (what's there),
panel B = target magnitude (what DFT sees),
panel C = residual overlay (what we can't explain).

---

## 6. Residual tensor butterflies

**Question:** What angular pattern can't the model explain?

**Data fields:** `ridge_residual_t2` (5 T2 components) and
`target_t2` (5 T2 components of DFT delta).

**Overlay:** Conformal tensor butterflies — the same glyph pipeline
the viewer already uses for calculator T2 output.  Two modes:

- **Target butterfly:** rendered from `target_t2`.  Shows the DFT
  ground truth angular pattern at each atom.  This is what the
  aromatic residue does to the shielding tensor.

- **Residual butterfly:** rendered from `ridge_residual_t2`.  Shows
  what the 91-kernel ridge CANNOT explain.  Where the model works,
  this butterfly is small.  Where it fails, the butterfly is large
  and points in the unexplained direction.

**Thesis value:** Side-by-side panels on the same protein: panel A
shows the target butterflies (what DFT sees), panel B shows the
residual butterflies (what's left after calibration).  The residual
butterflies near TRP should be small and randomly oriented.  Near
HIE they should be larger and systematically oriented — pointing
toward the imidazole nitrogen positions, revealing the physics the
ring current model misses.

This is the most powerful figure in the thesis: it shows the angular
structure of the model's ignorance, not just its magnitude.

---

## 7. Before/after φ comparison

After implementing the azimuthal angle (spec/AZIMUTHAL_ANGLE.md),
re-run viewer_export on the same protein.  The self_fit_residual
values will change for atoms near HIE (hypothesis: they shrink).
Overlay the before and after residual shells side by side, same
protein, same viewpoint.  The difference is the visual proof that
φ helps.

---

## Implementation notes

The viewer already renders per-atom translucent overlays for the
tensor butterflies.  These overlays use the same VTK pipeline:
vtkSphereSource per atom, opacity set by the annotation value,
colour from a lookup table.  The JSON annotations map directly to
the per-atom scalar arrays that VTK actors consume.

Reality check on load:
```
for each atom_index in annotations:
    assert viewer_element[atom_index] == annotation["element_number"]
    assert |viewer_position[atom_index] - annotation["position"]| < 0.01 A
```

If any check fails, refuse to render — the extraction and the
viewer are looking at different proteins or different conformations.
