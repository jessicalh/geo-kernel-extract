# H5 Field Glossary — what is in the analysis file, what it means, and how it should be shown

**Status:** living document, written 2026-04-17. Authoritative schema:
`fileformat/analysis_file.h`. Authoritative writer:
`src/AnalysisWriter.cpp`. If the schema header disagrees with this
document the header wins — file a bug against the glossary.

**Target audience.** Professors, advisers, collaborators reviewing the
trajectory data through `h5-reader`. The document assumes graduate-
level familiarity with NMR, magnetic shielding, electrostatics, and
MD simulation, but does not assume familiarity with this codebase. One
line per field gives the physics; the presentation recommendation
says how the reader should surface it so a reviewer can vet it on
first sight.

**Why a glossary.** The reader is a visibility surface. The whole
point of the analysis H5 is that it is exhaustive and self-contained
— everything the conformation knows is in it — so an adviser should
be able to open any atom and see every physical quantity the kernels
say about it, in a representation that matches the quantity's
geometry. A scalar belongs on a line chart. A vector belongs on an
arrow. A traceless rank-2 tensor belongs on an oriented glyph. A
256-dimensional embedding belongs on a projection. A strip chart
cannot carry a tensor honestly; a glyph cannot carry a time series
honestly. Assigning the right modality per field is the difference
between "readable" and "vettable".

---

## How to read this document

### Shape symbols

- `T` — number of sampled frames in the trajectory.
- `N` — number of protein atoms.
- `R` — number of residues.
- `B` — number of covalent bonds.
- `n_rings` — number of aromatic rings in the protein.
- `K` — the per-atom ring stack size, `K = 6` (six nearest rings).
- Array shapes are written in row-major (C-order) as `(T, N, 3)` etc.

### SphericalTensor layout (the (…,9) datasets)

Every `(…, 9)` dataset is a `SphericalTensor` decomposition of a
rank-2 tensor:

| Index | Irrep | Components | Physical content |
|-------|-------|------------|------------------|
| 0     | T0    | 1 (scalar) | isotropic / trace part |
| 1–3   | T1    | 3 (pseudovector) | antisymmetric part |
| 4–8   | T2    | 5 (m=−2,−1,0,+1,+2) | traceless symmetric part |

The H5 writer asserts `layout = "T0,T1[3],T2[5]"` as a dataset
attribute on every spherical-tensor dataset. T2 preservation is
load-bearing — the glossary flags tensor fields whose T2 is a thesis
claim, not decoration.

### Sign conventions

- Shielding tensor σ defined as σ_ab = −∂B_a^sec / ∂B_{0,b}. A
  diamagnetic current (I < 0) above a ring produces σ > 0 (shielded).
  Verified against Case 1995 for a proton 3 Å above PHE.
- EFG V_ab = second derivative of the Coulomb potential; traceless
  by Laplace. Haeberlen ordering |V_zz| ≥ |V_yy| ≥ |V_xx| and
  asymmetry η = (V_xx − V_yy) / V_zz ∈ [0, 1] are the canonical
  reporting form.
- Positions in Å, charges in e, energy in kJ/mol except DSSP H-bond
  energy in kcal/mol, angles in radians, temperature in K, pressure
  in bar.

### The presentation modalities

For every field the glossary marks one modality as **primary**
(the best-case defensible display for advisers) and flags the
others as **useful** (a supporting view worth offering) or
**N/A** (genuinely wrong for this data shape). The modalities
progress from atom-local through surface-local, residue-local,
protein-wide, and abstract.

1. **Tensor glyph on atom.** A 3D shape drawn at the atom's
   position that carries the tensor's eigenstructure. The physics
   literature's answer is **superquadric glyphs** (Kindlmann 2004;
   Schultz & Kindlmann 2010): unlike ellipsoids they do not suffer
   profile ambiguity, they encode sign of eigenvalues, and they
   interpolate prolate / oblate / spherical / cuboid over
   anisotropy. For rank-2 *asymmetric* shielding tensors (BS / HM /
   McConnell / ring-chi / H-bond) the T2 irrep is rendered as a
   superquadric with T0 colour on the surface; T1 may ride as a
   small arrow through the glyph. For traceless symmetric EFGs
   (Coulomb, APBS, water, AIMNet2) Haeberlen principal axes V_zz
   arrow + oriented disk for η is the cleanest reading.
2. **Colour bubble on atom.** The atom sphere recoloured by a
   scalar, with a legend and scale. For one-dimensional
   quantities that vary over the protein (T0 of a shielding
   tensor, partial charge, SASA as a quick-look). Must be paired
   with a colour-bar actor; diverging Moreland for signed fields,
   sequential viridis for non-negative fields. When the
   protein-wide spatial pattern is the story rather than per-atom
   numbers, prefer the SAS surface modality below — a continuous
   coloured mesh reads the pattern where a sea of coloured spheres
   degrades.
3. **SAS surface colouring.** The solvent-accessible surface of
   the protein (Shrake-Rupley mesh, `vtkMolecularSurfaceFilter`
   or marching cubes over a vdW density grid) painted with a
   per-atom or per-surface-point scalar. The PyMOL / Chimera /
   ChimeraX idiom every structural biologist reads fluently.
   Natural targets: Mahalanobis outlier-distance, ridge residual,
   hydration-dipole alignment and coherence, per-residue bonded
   strain summed onto surface atoms.
4. **Ribbon thickness modulation.** The backbone ribbon
   (`vtkProteinRibbonFilter`) with a per-residue scalar driving
   local width — the "sausage" representation from chemical-shift
   perturbation figures. The ribbon is drawn anyway; widening
   strong-signal regions is the cleanest visual hierarchy from a
   distance. Thickness and colour can carry different channels
   (e.g. thickness = Mahalanobis, colour = secondary structure).
   Natural targets: DSSP H-bond strength, per-residue bonded
   energy, per-residue Mahalanobis, per-residue chemical shift
   perturbation.
5. **Volumetric scalar isosurface.** A 3D scalar field sampled
   on a grid around the whole protein, contoured at one or more
   values — the cumulative environment picture. The reader
   already renders ring-local BS/HM butterflies and per-ring
   B-field streamlines; the protein-wide version sums those grids
   over every ring. Natural target: |total_B_field| isosurface at
   e.g. ±0.1 mT. Translucent, toggleable; belongs alongside the
   existing butterfly / streamline controls.
6. **3D or n-D graph.** Anything that does not live on the
   molecule itself: 256-D AIMNet2 embeddings (UMAP / PCA
   projection), dihedral joint distributions (Ramachandran
   scatter), per-ring stacks (K=6 bar chart by distance),
   frame-wise trace through PCA subspace, 8-class DSSP residue ×
   time heatmap. These need their own dock, not a strip chart.
7. **Scalar line chart.** One scalar vs frame with the cursor at
   the current frame — `QtAtomTimeSeriesDock` today. Right for
   trend and variability; wrong for any quantity with direction.

**When 3D relief is wrong.** Relief helps exactly when the data
is spatially situated on the molecule. It hurts for: Ramachandran
density (2D + contour lines beats a 3D bump); DSSP heatmaps
(categorical, no magnitude to raise); the 55-slot MLP
kernel-weight strip (too many spikes, self-occlusion);
system-energy line panels (time should stay linear). Don't 3D-ify
these just because the canvas can.

### What already exists in the reader

The recommendations below call for primitives, not features. The
reader today (2026-04-17) provides:

- scalar line chart for ~35 per-atom scalars (`QtAtomTimeSeriesDock`)
- atom inspector with per-frame rebuild of a tree of all per-atom
  fields, including T2 decompositions under a group node
- atom picker, selection overlay, backbone ribbon, ring polygons
- Biot-Savart / Haigh-Mallion butterfly isosurfaces (re-evaluated
  from ring geometry, not from the H5)
- B-field streamlines around rings
- `QtPlaybackController` for frame scrub
- structured UDP log on port 9997

It does **not** yet have atom colour overlays, atom-centred tensor
glyphs, E-field / B-field per-atom arrows, ring-indexed K=6
charts, a Ramachandran dock, a UMAP/PCA dock for the AIMNet2
embedding, or a per-residue heat-strip. For each H5 field, the
"Primary" recommendation below names the modality an adviser
should see; "Available today" flags whether that modality is
already wired up or whether adding it is a mechanical table entry
(line-chart / inspector), a new overlay class (colour bubble,
tensor glyph, vector arrow — duplicate the `QtRingPolygonOverlay`
pattern), or a new dock (Ramachandran, embedding projector,
K=6 ring chart).

---

# /meta — trajectory temporal metadata

Scope: one-per-frame, with protein-level attributes.

- **`frame_times (T,) float64, ps.** Simulation time at each
  sampled frame, read from GROMACS EDR. This is the horizontal
  axis on every time-series chart and the canonical way to talk
  about "when" an event happens.
  *Primary:* scalar line chart axis (every other field is plotted
  against this). *Useful:* frame-to-time conversion label in the
  status bar for `QtPlaybackController`. *N/A:* glyph / bubble /
  projection.
- **`frame_indices (T,) int32.** Original XTC frame index for
  each sampled frame (accounts for the `stride` attribute). Lets
  an adviser correlate the H5 to the underlying trajectory file.
  *Primary:* status-bar readout. *N/A:* glyph / bubble / projection
  / chart. Advisers should not look at frame indices in a chart;
  they want ps.
- **Attributes:** `protein_id` (string, e.g. `1CBH_192`),
  `n_atoms`, `n_frames`, `n_residues`, `stride` (every *stride*-th
  XTC frame is kept). Displayed in the title bar and the File dock.

---

# /atoms — static per-atom properties

Scope: per-atom, frame-invariant. These drive identity, filtering,
and colouring but they do not change over time. Typed enums are
documented in `h5-reader/src/model/Types.h`.

- **`element (N,) int32, atomic number.** 1=H, 6=C, 7=N, 8=O,
  16=S. Sets the element-dependent physics (Buckingham A/B, γ,
  magnetogyric ratio) and drives the default molecule colouring.
  *Primary:* element-colour on the `vtkMolecule` (already the
  default). *Useful:* filter dropdown for every other overlay
  ("only show H atoms"). *N/A:* chart, glyph, projection —
  element is identity, not a quantity.
- **`residue_index (N,) int32.** 0-based residue ordinal. Maps
  atom→residue for all per-residue data.
  *Primary:* atom inspector header ("ALA 42, chain A").
  *N/A* elsewhere.
- **`atom_name (N,) string.** PDB label ("CA", "CB", "N", "HB2").
  Human identification. Shown in inspector title and selection
  overlay label.
- **`atom_role (N,) int32 enum.** `AtomRole` — BackboneN,
  BackboneCA, BackboneC, BackboneO, sidechain heavy, hydrogens.
  *Primary:* inspector row + filter ("only backbone"). *Useful:*
  stratified time-series ("compare all backbone CA" — multi-atom
  overlay, requires scope expansion per
  `notes/TIME_SERIES_EXPANSION.md`).
- **`hybridisation (N,) int32 enum.** sp / sp² / sp³ / aromatic
  / unassigned. Inspector row; feeds bond-category logic but is
  itself a static label. *N/A* glyph / chart.
- **`n_bonded (N,) int32.** Covalent degree. Inspector row;
  sanity check for missing/extra bonds. *N/A* elsewhere.
- **`graph_dist_ring (N,) int32.** Bond-hops to the nearest
  aromatic ring (0 = on-ring, 1 = bonded-to-ring, …). Classic
  NMR shielding covariate.
  *Primary:* inspector row. *Useful:* atom colour bubble (sequential
  viridis; buried-in-ring atoms bright) for a ring-context
  overview. *N/A:* chart, glyph, projection.
- **`is_backbone (N,), is_conjugated (N,) int8.** Topological
  flags. Inspector rows; drive ribbon rendering. *N/A* chart.
- **`is_amide_H`, `is_alpha_H`, `is_methyl`, `is_aromatic_H`,
  `is_on_aromatic_residue`, `is_hbond_donor`, `is_hbond_acceptor`,
  `parent_is_sp2` (N,) int8.** Enrichment flags from the molecular
  graph. Per feedback_ui_data, expose all of them in the
  inspector; do not paper over them with a summary. They are the
  NMR-interesting partition of the atom set and should also appear
  as **filter chips** so a reviewer can say "only amide Hs" and
  see every overlay recompute.
- **`graph_dist_N`, `graph_dist_O (N,) int32.** Bond-hops to the
  nearest N / O heavy atom. Inspector rows.
  *Useful:* atom colour bubble for "polar-group nearness" (sequential
  colormap; cuts off above 5 hops). *N/A:* chart, glyph.
- **`eneg_sum_1`, `eneg_sum_2 (N,) float64.** Pauling
  electronegativity summed over first / second coordination shells.
  An integrated chemical environment polarity scalar.
  *Primary:* atom colour bubble (sequential viridis). *Useful:*
  inspector row; N=2 chart comparing static eneg against dynamic
  shielding components.
- **`n_pi_bonds_3 (N,) int32.** Count of π bonds within 3 bonds.
  Static delocalisation proxy. Inspector row + filter.
- **`bfs_to_nearest_ring (N,) int32.** Breadth-first-search distance
  to nearest ring (may differ from `graph_dist_ring` in fused-ring
  topology). Inspector row.
- **`bfs_decay (N,) float64.** exp(−bfs / scale), scale from TOML.
  Smooth ring-proximity weight used as a feature.
  *Primary:* atom colour bubble (sequential, 0 to 1). *Useful:*
  inspector row.
- **`partial_charge (N,) float64, e.** ff14SB / CHARMM static
  partial charge. The Coulomb EFG is built from this.
  *Primary:* atom colour bubble (diverging Moreland, symmetric
  about 0; −0.8 to +0.4 typical). *Useful:* inspector row.
- **`vdw_radius (N,) float64, Å.** Bondi / Cordero vdW radius.
  Used by SASA and dispersion. Inspector row.

---

# /residues — static per-residue properties

Scope: per-residue, frame-invariant.

- **`residue_name (R,) string, three-letter code.** "ALA", "GLY",
  "PRO". Drives ribbon colour and inspector title.
  *Primary:* ribbon colour-by-residue (already wired via
  `AminoAcid` enum dispatch). *Useful:* residue picker dropdown.
- **`residue_number (R,) int32.** PDB numbering (not necessarily
  contiguous). Inspector header + x-axis on any per-residue
  strip / heatmap.
- **`chain_id (R,) string.** Chain identifier; inspector header.

---

# /topology — bond graph and ring membership (CSR)

Scope: static, loaded once per file.

- **`bond_atoms (B, 2) int32.** Atom-pair list of bonds.
  *Primary:* molecular bond rendering (already via
  `vtkMolecule::AppendBond`). *Useful:* bond-centric inspector
  when a bond is picked (future; today the reader picks atoms,
  not bonds).
- **`bond_category (B,) int32 enum.** `BondCategory` — PeptideCO,
  PeptideCN, BackboneOther, SidechainCO, Aromatic, Disulfide,
  SidechainOther, Unknown. Drives McConnell `Δχ` per category.
  *Primary:* bond colour on the molecule (needs a legend dock).
  *Useful:* filter chip. Advisers asking "are the aromatic bonds
  correctly categorised?" need to see this on the 3D scene.
- **`bond_order (B,) int32 enum.** Single, Double, Triple,
  Aromatic, Peptide, Unknown. Drives line style on the molecule
  (double = two lines, aromatic = dashed ring). Already wired.
- **`parent_atom_index (N,) int32.** For every atom that is an H,
  the index of the heavy atom it is bonded to; −1 otherwise.
  Used to build X-H bond axes for Buckingham projection.
  Inspector row.
- **`ring_type (n_rings,) int32 enum.** `RingTypeIndex` — PHE, TYR,
  TRP6, TRP5, TRP9, HIS, HID, HIE. Drives ring-current intensity
  per ring. *Primary:* ring polygon colour (already wired via
  `QtRing` class hierarchy). *Useful:* per-ring inspector field
  when a ring is selected (future).
- **`ring_residue (n_rings,) int32.** Parent residue for each ring.
  Inspector row on the ring.
- **`ring_fused_partner (n_rings,) int32.** −1 if none; otherwise
  the paired ring index (TRP5 ↔ TRP6 share an edge).
  Inspector row. *Useful:* ring polygon line style change for
  fused rings so an adviser sees TRP's three-ring structure
  immediately.
- **`ring_offsets (n_rings+1,) int32.** CSR offsets into
  `ring_atom_indices`. Internal plumbing; no display.
- **`ring_atom_indices (Σ ring sizes,) int32.** Flat atom indices
  for all rings, in polygon order. Used by ring polygon, BS/HM
  butterfly, quadrature. *Primary:* ring polygon vertices
  (already wired). *N/A* elsewhere.

---

# /positions — per-frame coordinates

- **`xyz (T, N, 3) float64, Å.** Cartesian positions after PBC
  unwrap. Every per-frame overlay reads this.
  *Primary:* the molecule. `MoleculeScene::setFrame` drives
  everything. *Useful:* per-atom scalar line chart of x / y / z
  is already in the Desc table (not the best use of a chart slot,
  but legal). *Primary for a different view:* centroid trace or
  RMSD plot (system-wide), if added.

---

# /ring_current — Biot-Savart (BS) + Haigh-Mallion (HM) + ring-susceptibility (RS)

The ring-current group is the **largest and most-visible of the
shielding contributions** — dominates H R² in Stage 1. The 8-slot
per-ring-type axis is ordered by `RingTypeIndex`.

## Total shielding tensors (the user-facing T2 story)

- **`bs_shielding (T, N, 9) float64, ppm.** Full BS shielding
  tensor summed over all rings. Asymmetric (G = n ⊗ B), so T0,
  T1, T2 are all nonzero. T0 is the classical Johnson-Bovey
  isotropic shift; T2 is the angular residual that the thesis
  cares about.
  *Primary:* **superquadric tensor glyph at atom**, T2 shape,
  T0 surface colour (diverging Moreland; red = deshielded, blue
  = shielded). This is the thesis's money shot: an adviser should
  see every H atom within 6 Å of a ring wearing an oriented glyph.
  *Useful:* T0 and |T2| already in the time-series strip chart
  (extend to eigenvalue traces δ_ZZ, δ_YY, δ_XX once glyph lands).
  Atom colour bubble by T0 is a light-weight global scan.
  *Provenance:* `BiotSavartResult::SampleShieldingAt`.
- **`hm_shielding (T, N, 9) float64, ppm.** Same structure,
  surface-integral model (Haigh-Mallion).
  *Primary:* **tensor glyph** next to the BS glyph (toggle BS /
  HM / sum). Advisers comparing BS vs HM at close range need the
  same primitive with a model-selector.
  *Provenance:* `HaighMallionResult`.
- **`rs_shielding (T, N, 9) float64, ppm.** Ring-susceptibility
  contribution (bulk ring treated as point dipole at centre, full
  McConnell form). Asymmetric, all irreps nonzero.
  *Primary:* tensor glyph (same renderer; different model
  selector). Typically small compared to BS/HM.

## Per-ring-type stacks (the per-type decomposition)

- **`bs_T0_per_type (T, N, 8), hm_T0_per_type (T, N, 8) float64,
  ppm.** Isotropic contribution from each of 8 ring types. Most
  slots are zero at most atoms (H's near a PHE see only the PHE
  column). These are the learnable per-type coefficients from
  Stage 1 calibration.
  *Primary:* **3D/n-D chart** — a per-atom horizontal bar chart
  with 8 bars, one per ring type, coloured by type. Dock pattern
  `QtAtomRingTypeBarChart`, one entry per atom per frame, updates
  on frame advance. *Useful:* time-series for a single
  (atom × type) slot, extending the Desc table. *N/A:* glyph (no
  orientation), colour bubble (one scalar per atom, but loses the
  type decomposition).
- **`bs_T2_per_type (T, N, 8, 5), hm_T2_per_type (T, N, 8, 5)
  float64, ppm.** T2 per ring type. 40 scalars per atom per
  frame.
  *Primary:* **superquadric glyph per ring type** — 8 small glyphs
  stacked by ring type in the inspector dock, one per nonzero
  type. The inspector should render each as its own tiny tensor
  glyph widget with `Δχ / Ω / κ` printed alongside, so an adviser
  can see "PHE2 contributes this anisotropy, TYR0 contributes
  that". A full 3D glyph per type per atom is visual overload;
  keep that restricted to the total.

## Ring-proximity scalars

- **`n_rings_3A (T, N), n_rings_5A, n_rings_8A, n_rings_12A
  int16.** Ring centre count within R Å.
  *Primary:* scalar line chart (already in `QtAtomTimeSeriesDock`).
  *Useful:* atom colour bubble for `n_rings_5A` — "which atoms
  are in the aromatic cluster?". Integer colormap (0 / 1 / 2 / 3+).
- **`mean_ring_dist (T, N), nearest_ring_atom (T, N) float64, Å.**
  Spatial-context scalars.
  *Primary:* scalar line chart. *Useful:* atom colour bubble.
- **`G_iso_exp_sum (T, N) float64, ppm.** Exponentially weighted
  sum of isotropic ring kernels.
  *Primary:* scalar line chart. *Useful:* colour bubble.
- **`G_T2_exp_sum (T, N, 5) float64, ppm.** Exponentially
  weighted T2 (5 components).
  *Primary:* derived scalar |T2| on the line chart; full 5-vector
  in the inspector with a mini-glyph. *Useful:* tensor-glyph
  rendering using the 5-component T2 directly; a side-by-side
  small-multiples panel "all five channels of this atom" works
  well for advisers wanting to see angular-pattern evolution.
- **`G_iso_var_8A (T, N) float64, ppm².** Variance of ring-kernel
  contributions within 8 Å — a measure of cancellation vs
  reinforcement between nearby rings.
  *Primary:* scalar line chart. *Useful:* colour bubble with a
  log scale (the variance can span four orders of magnitude).

## Vector fields

- **`total_B_field (T, N, 3) float64, Tesla.** Summed secondary
  magnetic field from all ring currents at each atom.
  *Primary:* **protein-wide |B| volumetric isosurface** —
  reconstruct a 3D grid of |B(r)| from the same machinery that
  drives the per-ring butterfly (sum of `QtBiotSavartCalc` over
  every ring), contour at ±0.1 mT with translucent shielded /
  deshielded surfaces. This is the textbook cumulative-magnetic-
  environment picture and the one an NMR reviewer expects to see
  once the per-ring butterflies are in place. *Useful:* per-atom
  vector arrow when a single atom is picked (local direction +
  magnitude); |B| scalar line chart for time evolution.
  *Provenance:* sum over `BiotSavartResult` per-ring B evaluations.

---

# /efg — electric field gradient and derived shielding

Pure-electric tensors (Coulomb, APBS, water, AIMNet2 derived from
charges) are **symmetric and traceless** — pure T2, no T0, no T1.
Haeberlen ordering |V_zz| ≥ |V_yy| ≥ |V_xx|, asymmetry
η = (V_xx − V_yy) / V_zz. This is the canonical EFG reporting
convention and should drive the glyph.

## Classical Coulomb (ff14SB / CHARMM charges)

- **`coulomb_total (T, N, 9) float64, V/Å².** Total Coulomb EFG
  from all protein charges.
  *Primary:* **Haeberlen principal-axis glyph** — V_zz arrow
  through the atom, flattened disk perpendicular with aspect
  ratio set by η (η=0 is a circle, η=1 is a line). Colour by
  |V_zz|.
  *Useful:* atom colour bubble by η (viridis, 0→1). Scalar line
  chart of |T2|. *N/A:* colour bubble of V_zz alone would lose
  the sign.
  *Provenance:* `CoulombResult`, dipolar kernel summed over ff14SB
  charges.
- **`coulomb_backbone (T, N, 9), coulomb_aromatic (T, N, 9)
  float64, V/Å².** EFG decomposed by source subset.
  *Primary:* selector on the Haeberlen glyph ("total / backbone /
  aromatic / sidechain" toggle). *Useful:* side-by-side small
  multiples of three glyphs per atom in the inspector so advisers
  see how backbone and aromatic contributions compose the total.
- **`coulomb_shielding (T, N, 9) float64, ppm.** Shielding from
  applying the Buckingham A·E + γ·V mapping. Has T0 (from A·E +
  B·E² scalar), T2 from γ·V. NOT pure T2.
  *Primary:* **tensor glyph** (superquadric; same renderer as
  BS/HM). *Useful:* T0 colour bubble + |T2| line chart.

## E-field vectors

- **`E_total (T, N, 3), E_backbone, E_sidechain, E_aromatic,
  E_solvent (T, N, 3) float64, V/Å.** E-field decomposed by
  source. `E_solvent = apbs_efield − E_total` (what solvation
  adds).
  *Primary:* **vector arrow at atom** for E_total; secondary
  arrows for the decomposition pieces on hover or in the
  inspector with small-multiple arrows. Arrow colour by source
  (by convention: backbone = grey, sidechain = tan, aromatic =
  purple, solvent = cyan, total = black).
  *Useful:* |E| line chart (already in Desc table).
- **`E_magnitude (T, N) float64, V/Å.** |E_total|.
  *Primary:* scalar line chart + colour bubble. The colour bubble
  is especially useful for spotting strong-field regions in one
  glance.
- **`E_bond_proj (T, N) float64, V/Å.** Projection of E along
  the primary bond axis (the Buckingham σ_iso input).
  *Primary:* scalar line chart. *Useful:* colour bubble
  (diverging, signed).
- **`E_backbone_frac (T, N) float64, ∈ [0, 1].** Fraction of
  |E_total| supplied by backbone charges.
  *Primary:* scalar line chart. *Useful:* colour bubble
  (sequential viridis, 0 to 1).

## APBS solvated

- **`apbs_efg (T, N, 9) float64, V/Å².** APBS-derived EFG (full
  Poisson-Boltzmann solve with dielectric screening).
  *Primary:* Haeberlen principal-axis glyph (same renderer as
  `coulomb_total`). *Useful:* side-by-side comparison with
  `coulomb_total` to show solvent screening (APBS magnitudes
  typically smaller than vacuum Coulomb near surface atoms).
  *Provenance:* `ApbsEfgResult`; APBS replaces vacuum Coulomb at
  ~4 s vs 25 s per frame on 4876-atom proteins.
- **`apbs_efield (T, N, 3) float64, V/Å.** APBS E-field vector.
  *Primary:* vector arrow at atom, in a distinct colour from
  vacuum E_total so advisers see the screening direction.

## AIMNet2 charge-derived

- **`aimnet2_total (T, N, 9), aimnet2_backbone, aimnet2_aromatic,
  aimnet2_shielding (T, N, 9) float64.** EFG and shielding
  derived from AIMNet2 Hirshfeld charges — same kernel as
  `coulomb_total` but charges are geometry-responsive (from
  `AimNet2Result`). `aimnet2_shielding` is the Buckingham-mapped
  shielding including T0 (same shape as `coulomb_shielding`).
  *Primary:* Haeberlen glyph for EFG fields; superquadric
  shielding glyph for `aimnet2_shielding`. Same renderers as the
  Coulomb family; toggle the charge source. *Useful:*
  `Δ = AIMNet2 − Coulomb` as an atom colour bubble would show
  the polarisation signal directly; that is an inspector-level
  derived view, not a raw H5 field.

---

# /bond_aniso — McConnell bond magnetic anisotropy

McConnell is the **dominant T2 contributor** per
`GEOMETRIC_KERNEL_CATALOGUE.md` — hundreds of bonds, 1/r³ decay,
large Δχ for C=O. Full formula per bond (asymmetric,
non-traceless, T0+T1+T2 nonzero).

- **`mc_shielding (T, N, 9) float64, ppm.** Summed McConnell
  shielding over all bonds.
  *Primary:* **superquadric tensor glyph at atom**. T0 surface
  colour, T2 shape, T1 as small perpendicular arrow. This is the
  thesis's second-most-important glyph after ring current.
  *Useful:* T0 + |T2| line chart (already in Desc table).
  *Provenance:* `McConnellResult`.
- **`T2_backbone, T2_sidechain, T2_aromatic, T2_CO_nearest,
  T2_CN_nearest (T, N, 9) float64, ppm.** Shielding decomposed
  by bond category plus the two "single nearest bond" tensors
  (one C=O, one C-N).
  *Primary:* **selector on the McConnell glyph** — adviser picks
  total / backbone / sidechain / aromatic / nearest-CO /
  nearest-CN. All five rendered simultaneously in the inspector
  as small-multiple glyphs so the adviser sees which category
  contributes what shape.
  *Useful:* per-category T0 already implicit in `co_sum` etc.
- **`co_sum, cn_sum, sidechain_sum, aromatic_sum (T, N) float64,
  ppm.** McConnell scalar sums per category (the (3cos²θ−1)/r³
  sum, the T0 contribution before Δχ multiplication).
  *Primary:* scalar line chart. *Useful:* atom colour bubble,
  category-selectable (diverging Moreland).
- **`co_nearest (T, N) float64, ppm.** McConnell scalar from the
  single nearest C=O.
  *Primary:* scalar line chart. *Useful:* colour bubble.
- **`nearest_CO_dist, nearest_CN_dist (T, N) float64, Å.**
  Distance to the nearest bond of that category.
  *Primary:* scalar line chart (H-bond-like). *Useful:* colour
  bubble with a sequential colormap capped at 6 Å.
- **`nearest_CO_midpoint (T, N, 3), dir_nearest_CO (T, N, 3)
  float64, Å / unit vector.** Midpoint position of the nearest
  C=O and unit direction from atom to that midpoint.
  *Primary:* **render a thin line segment from the atom to the
  C=O midpoint** in the 3D scene — toggleable overlay, only drawn
  when an atom is selected. Geometric context for the nearest-CO
  tensor, and the exact input to the nearest-CO glyph.
  *Useful:* inspector row with the 3-vector plus distance.

---

# /quadrupole — π-electron quadrupole EFG

Aromatic rings carry an axial quadrupole moment; creates an EFG
via Stone's T-tensor formalism. Traceless, symmetric, pure T2.

- **`pq_shielding (T, N, 9) float64, ppm.** Total
  π-quadrupole-derived shielding.
  *Primary:* tensor glyph at atom (same superquadric renderer).
  Typically small magnitude relative to ring current, so the
  glyph should auto-scale per-channel.
  *Useful:* T0 + |T2| line chart.
  *Provenance:* `PointQuadrupoleResult`.
- **`pq_T0_per_type (T, N, 8), pq_T2_per_type (T, N, 8, 5)
  float64, ppm.** Per-ring-type decomposition.
  *Primary:* per-type horizontal bar chart for T0 (same pattern
  as ring-current per-type bars); inspector-side tiny glyphs for
  T2 per type.

---

# /dispersion — van der Waals C6 anisotropy

Small anisotropic contribution from 1/r⁸ dispersion over ring
vertices. Dominates H R² slightly (0.387 in Stage 1) though most
of that is scalar, not T2.

- **`disp_shielding (T, N, 9) float64, ppm.** Summed dispersion
  shielding.
  *Primary:* tensor glyph (same renderer, auto-scaled).
  *Useful:* T0 + |T2| line chart.
  *Provenance:* `DispersionResult`.
- **`disp_T0_per_type (T, N, 8), disp_T2_per_type (T, N, 8, 5)
  float64, ppm.** Per-ring-type decomposition.
  *Primary:* per-type bar chart for T0; inspector per-type
  glyphs for T2.

---

# /hbond — H-bond dipolar shielding

Weak in Stage 1 (R²≈0.02 for H) but present. Same McConnell form
as bond anisotropy, with the D-H…A direction replacing the bond
axis.

- **`hbond_shielding (T, N, 9) float64, ppm.** Summed H-bond
  contribution (asymmetric, all irreps).
  *Primary:* tensor glyph. *Useful:* T0 line chart.
  *Provenance:* `HBondResult`.
- **`nearest_spherical (T, N, 9) float64, ppm.** Tensor from the
  single strongest H-bond partner.
  *Primary:* tensor glyph in the inspector, as a small-multiple
  next to the total hbond glyph so advisers see whether a single
  partner dominates.
- **`nearest_dist (T, N) float64, Å.** Distance to nearest partner.
  *Primary:* scalar line chart. *Useful:* colour bubble, sequential
  colormap capped at 3.5 Å.
- **`nearest_dir (T, N, 3) float64, unit vector.** Direction
  to nearest partner.
  *Primary:* **line segment or small arrow** drawn from atom to
  partner when the atom is selected. Same pattern as
  `dir_nearest_CO`. Only rendered on selection.
- **`inv_d3 (T, N) float64, Å⁻³.** 1/r³ at the nearest-partner
  distance — the geometric factor in the Cornilescu-Bax form.
  *Primary:* scalar line chart. *Useful:* colour bubble.
- **`count_3_5A (T, N) int16.** Number of H-bond candidates
  in 3–5 Å.
  *Primary:* scalar line chart (integer step colormap).
- **`is_donor (T, N), is_acceptor (T, N), is_backbone (T, N)
  int8.** Per-frame flags (protonation could change these in
  principle; rarely does in a standard MD).
  *Primary:* boolean state in inspector. *Useful:* filter chips
  ("only donors").

---

# /sasa — solvent-accessible surface area

- **`sasa (T, N) float64, Å².** Shrake-Rupley SASA.
  *Primary:* **atom colour bubble**, sequential viridis, 0–80 Å².
  Classic protein-scale display. *Useful:* scalar line chart.
  *Provenance:* `SasaResult`.
- **`normal (T, N, 3) float64, unit vector.** Outward surface
  normal (only meaningful for exposed atoms; buried atoms may
  have degenerate normals).
  *Primary:* **vector arrow** on exposed atoms only (filter on
  `sasa > 0.5`). Arrow length fixed, colour by `sasa_asymmetry`
  from /water (surface directionality).

---

# /water — explicit solvent fields (TIP3P)

Water E-field and EFG computed directly from water atom charges;
plus hydration-geometry scalars. The hydration story is the newest
signal in the schema (Stage 2).

- **`efield (T, N, 3) float64, V/Å.** E-field from all waters
  within cutoff.
  *Primary:* vector arrow (third colour alongside protein E and
  APBS E in the inspector).
  *Useful:* |E| scalar line chart.
- **`efg (T, N, 9) float64, V/Å².** EFG from all waters,
  traceless symmetric.
  *Primary:* Haeberlen principal-axis glyph (same renderer as
  `coulomb_total`).
- **`efield_first (T, N, 3), efg_first (T, N, 9) float64.**
  First-shell-only versions (≤ 3.5 Å).
  *Primary:* same renderers, model selector "all waters / first
  shell". `efield − efield_first ≈ bulk` contribution.
- **`n_first (T, N), n_second (T, N) int16.** Shell counts.
  *Primary:* scalar line chart. *Useful:* integer colour bubble
  (0 / 1 / 2 / 3+). An adviser's first question is "how many
  waters are in contact?".
- **`half_shell_asymmetry (T, N) float64, ∈ [−1, 1].** COM-based
  hydration asymmetry.
  *Primary:* scalar line chart. *Useful:* colour bubble
  (diverging Moreland, symmetric about 0).
  *Provenance:* `HydrationShellResult`.
- **`dipole_cos (T, N) float64, ∈ [−1, 1].** Mean cos of water
  dipole orientation relative to atom→water direction.
  *Primary:* scalar line chart. *Useful:* diverging colour bubble.
- **`nearest_ion_dist, nearest_ion_charge (T, N) float64, Å / e.**
  Closest counter-ion + its formal charge.
  *Primary:* scalar line chart (dist). *Useful:* colour bubble
  capped at 5 Å; glyph at ion position would require rendering
  ion atoms which is not yet wired.
- **`dipole_vector (T, N, 3) float64.** Net dipole direction
  from first-shell waters.
  *Primary:* vector arrow at atom (small, fourth colour).
- **`surface_normal (T, N, 3) float64, unit vector.** SASA-derived
  outward normal used by `HydrationGeometryResult`. Same as
  `/sasa/normal` for most atoms; overlap intentional, separate
  storage documents the provenance.
- **`sasa_asymmetry (T, N) float64.** SASA-frame hydration
  asymmetry (distinct from `half_shell_asymmetry` — different
  reference direction).
  *Primary:* scalar line chart. *Useful:* colour bubble.
- **`sasa_dipole_align (T, N) float64, ∈ [−1, 1].** Alignment of
  average water dipole with SASA normal.
  *Primary:* **SAS surface colouring** (diverging Moreland) —
  ordered-hydration patches stand out as surface regions.
  *Useful:* scalar line chart; atom colour bubble as a fallback.
- **`sasa_dipole_cohere (T, N) float64, ∈ [0, 1].** Coherence
  magnitude of the water-dipole directions.
  *Primary:* **SAS surface colouring** (sequential viridis) —
  ice-like vs bulk-like hydration visible at a glance on the
  surface. *Useful:* scalar line chart. High coherence = ordered
  hydration shell.
- **`sasa_first_shell_n (T, N) int16.** First-shell count from
  the SASA-normal frame; differs slightly from `n_first`.

---

# /charges — geometry-responsive charge models

- **`aimnet2_charge (T, N) float64, e.** AIMNet2 Hirshfeld
  charges. Responsive to conformation.
  *Primary:* **atom colour bubble**, diverging Moreland,
  symmetric about 0. Advisers reading polarisation signal want
  this on the protein, not on a chart.
  *Useful:* scalar line chart (already in Desc table).
  *Provenance:* `AimNet2Result`; float32 in embedding but float64
  for per-atom charge.
- **`eeq_charge (T, N) float64, e.** Caldeweyher EEQ
  (electronegativity equilibration) charge.
  *Primary:* colour bubble. *Useful:* line chart + `aimnet2_charge
  − eeq_charge` as a derived diagnostic.
- **`eeq_cn (T, N) float64, dimensionless.** EEQ coordination
  number (soft covalent count).
  *Primary:* scalar line chart. *Useful:* colour bubble.

---

# /aimnet2_embedding — 256-D electronic embedding

- **`aim (T, N, 256) float32, dimensionless.** AIMNet2 hidden
  layer, 256 dimensions per atom per frame. Encodes the local
  electronic structure (orbital energies, occupations, charge
  response) that AIMNet2 uses to predict charges. Not
  human-readable per-component.
  *Primary:* **projection dock** — 2D UMAP (or PCA) of the
  embedding across (atom × frame), coloured by element. A trace
  overlay highlights the currently-picked atom's trajectory
  through embedding space. Use UMAP fit once at file-load
  (10-20 s on 600 × 300 × 256 = 45M points is marginal; subsample
  to ~10⁵ for fit and project the rest). A 2D canvas dock is the
  right shape — not a line chart and not a glyph.
  *Useful:* per-atom embedding **norm** (‖aim‖₂) as a scalar line
  chart; per-atom cosine with a chosen reference atom as a time
  series. *N/A:* glyph (no rank-2 structure), colour bubble on
  the raw 256 values (meaningless).
  *Provenance:* `AimNet2Result`. Stored float32 (feedback entry
  `feedback_embedding_float32`).

---

# /per_ring — K=6 nearest rings, per atom (ring-stratified kernels)

For every atom at every frame, the six nearest ring centres are
captured with their geometry and the per-ring contribution from
five calculators. This is the inspect-the-individual-ring view
that the summed shielding tensors lose.

- **`geometry (T, N, K, 6) float64, mixed units.** Six geometric
  scalars per ring for each of the K=6 nearest rings. The six
  field names are listed in `geometry_fields`.
  *Primary:* **ring-stratified table in the inspector** — six
  rows, one per nearest ring, with ring type chip (coloured by
  `ring_type` below) and the six scalars as columns. A small
  inset cartoon of the six rings at their true geometry relative
  to the picked atom would be an ambitious best-case.
  *Useful:* time series of the "nearest" slot's distance (k=0) is
  already useful on the scalar chart.
  *N/A:* glyph on atom (six structures at once would overwhelm).
- **`geometry_fields (6,) string.** Column labels for `geometry`.
  Read at file-load; drives the inspector column headers. Must
  not be ignored.
- **`ring_type (T, N, K) int8.** Ring type index for each of
  the K nearest rings. Drives ring-type colour chip in the
  inspector.
- **`bs_T2 (T, N, K, 5) float64, ppm.** Per-ring T2
  contribution to the BS summed shielding. Summing over K=6
  recovers most of `bs_shielding`'s T2 for atoms near the
  protein surface.
  *Primary:* **per-atom K=6 bar chart** in a dedicated dock —
  horizontal bars, one per nearest ring, length = |T2| for that
  ring, colour = ring type, order = distance ascending. Adviser
  sees at a glance "ring 0 dominates, rings 2-5 small". *Useful:*
  tiny per-ring glyphs in the inspector (the same small-multiple
  pattern as the per-type stacks).
- **`hm_T2 (T, N, K, 5), chi_T2 (T, N, K, 5), pq_T2 (T, N, K, 5),
  hm_H_T2 (T, N, K, 5) float64, ppm.** Same K=6 T2 stacks for
  Haigh-Mallion, ring-susceptibility, π-quadrupole, and
  H-atom-specific Haigh-Mallion.
  *Primary:* same K=6 bar-chart dock with a calculator-selector.
  Advisers comparing BS vs HM per-ring get this view for free.
- **`disp_scalar (T, N, K) float64, ppm.** Dispersion scalar per
  nearest ring.
  *Primary:* K=6 bar chart (scalar height, no T2 shape).
- **`disp_T2 (T, N, K, 5) float64, ppm.** Dispersion T2 per
  nearest ring.
  *Primary:* K=6 bar-chart dock entry (same pattern).

---

# /ring_geometry — per-ring, per-frame geometry

- **`data (T, n_rings, 7) float64.** Seven geometric scalars per
  ring per frame. Field names in `fields`. Typically
  [centre(3), normal(3), radius(1)] — whatever `fields` says at
  load time.
  *Primary:* the ring polygon overlay and BS/HM butterfly
  isosurfaces already read this. *Useful:* ring-scoped
  time-series dock (new scope in `QtAtomTimeSeriesDock`) for
  things like "radius of PHE1 over time" — scope expansion per
  `notes/TIME_SERIES_EXPANSION.md`.
- **`fields (7,) string.** Column labels. Drives the ring-scoped
  inspector.

---

# /bonded_energy — per-atom CHARMM36m bonded strain

Per-atom energy split evenly among participating atoms. Sums
exactly to the whole-system EDR terms (charge-conservation
equivalent check, validated 2026-04-15).

- **`bond (T, N) float64, kJ/mol.** Harmonic bond stretching.
  *Primary:* scalar line chart + atom colour bubble (log-scale
  sequential).
  *Useful:* outliers (>10 kJ/mol) flagged in inspector.
- **`angle (T, N) float64, kJ/mol.** Harmonic angle bending.
  *Primary:* scalar line chart + colour bubble.
- **`urey_bradley (T, N) float64, kJ/mol.** 1-3 harmonic
  restoring.
  *Primary:* scalar line chart.
- **`proper_dih (T, N) float64, kJ/mol.** Periodic proper
  dihedral.
  *Primary:* scalar line chart. *Useful:* colour bubble mapped
  onto the φ/ψ picture from /dihedrals.
- **`improper_dih (T, N) float64, kJ/mol.** Planarity restoring
  on sp² centres.
  *Primary:* scalar line chart. *Useful:* colour bubble flagging
  aromatic ring warping.
- **`cmap (T, N) float64, kJ/mol.** CHARMM map correction.
  *Primary:* scalar line chart. *Useful:* colour bubble.
- **`total (T, N) float64, kJ/mol.** Sum of the above.
  *Primary:* atom colour bubble (sequential viridis or log scale).
  Strain map over the protein at a glance.

Provenance: `BondedEnergyResult` evaluates CHARMM36m energy
functions from positions + TPR parameters; splits evenly among
participants.

---

# /energy — whole-system per-frame thermodynamic state

42 EDR terms; every one of these is **scope=System, not per-atom**.
The right dock is a single system-energy dock with all terms as
line charts against frame_time, not a table on every atom.

- **`coulomb_sr (T,), coulomb_recip (T,) float64, kJ/mol.** PME
  real-space and reciprocal-space Coulomb.
  *Primary:* line chart in the system-energy dock.
- **`bond, angle, urey_bradley, proper_dih, improper_dih, cmap_dih
  (T,) float64, kJ/mol.** Whole-system bonded terms.
  *Primary:* line chart. *Useful:* cross-check against
  per-atom-sum of `/bonded_energy` (should match to machine
  precision per 2026-04-15 validation).
- **`lj_sr (T,) float64, kJ/mol.** Short-range Lennard-Jones.
  *Primary:* line chart.
- **`potential, kinetic, enthalpy (T,) float64, kJ/mol.**
  *Primary:* line chart (on one axis so advisers see partitioning).
- **`temperature (T,) float64, K.** Kinetic temperature.
  *Primary:* line chart, reference line at thermostat setpoint.
- **`pressure (T,) float64, bar.** Virial-derived pressure.
  *Primary:* line chart, reference line at NPT setpoint.
- **`volume (T,) float64, nm³.** Box volume.
  *Primary:* line chart. *Useful:* detects NPT breathing vs NVT.
- **`density (T,) float64, kg/m³.** Mass density.
  *Primary:* line chart.
- **`box (T, 3) float64, nm.** Box edges.
  *Primary:* three line charts stacked (L_x, L_y, L_z).
- **`virial (T, 9) float64, kJ/mol.** 3×3 virial tensor flattened
  as XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ (layout confirmed by
  `virial_layout`).
  *Primary:* 3×3 matrix view in a dock; diagonal elements as
  line charts, off-diagonal as a secondary panel (equilibrated
  systems have ≈0 off-diagonal).
- **`virial_layout (9,) string.** Labels — drives the 3×3 view
  headers.
- **`pressure_tensor (T, 9) float64, bar.** Same 3×3 layout as
  virial. Diagonal → pressure; off-diagonal → shear stress.
  *Primary:* same 3×3 matrix view as virial.
- **`T_protein (T,), T_non_protein (T,) float64, K.** Per-group
  kinetic temperature.
  *Primary:* two lines on the temperature chart; divergence flags
  thermostat-coupling issues.

---

# /dihedrals — per-residue backbone and sidechain angles

All angles in radians, written by `DsspResult`. φ, ψ, ω are
backbone; χ1–χ4 are sidechain. cos/sin pairs are stored for
circular-statistics work.

- **`phi (T, R), psi (T, R) float64, rad.** Backbone dihedrals.
  *Primary:* **Ramachandran scatter dock** — one point per
  (residue, frame) in (φ, ψ) space, coloured by DSSP secondary
  structure or by residue type. Current-frame points highlighted.
  This is the canonical protein-dynamics view and advisers will
  expect to see it. *Useful:* per-residue time-series of φ alone
  or ψ alone (residue-scoped Desc table, scope expansion).
  *N/A:* glyph or colour bubble (angles on atoms, not residues,
  lose residue identity).
- **`omega (T, R) float64, rad.** Peptide-bond planarity
  (CA-C-N+1-CA+1). 180° for trans (normal), 0° for cis (rare,
  often before PRO).
  *Primary:* per-residue line chart with a reference at 180°.
  *Useful:* flag |ω − 180°| > 20° as twisted-peptide outliers
  (none were found in the 2026-04-15 validation, a sanity check
  the reader should make visible).
- **`chi1 (T, R), chi2, chi3, chi4 (T, R) float64, rad.**
  Sidechain dihedrals. NaN for residues without that dihedral
  (e.g. GLY has no χ1).
  *Primary:* **circular histogram** per residue (dock) — a polar
  plot of χ1 over frames shows rotamer occupancy immediately.
  One polar plot per selected residue, cycled via the residue
  picker. *Useful:* per-residue line chart (caveats around
  wrap-around at ±π; use cos/sin trick).
- **`chi1_cos, chi1_sin, chi2_cos, chi2_sin, chi3_cos, chi3_sin,
  chi4_cos, chi4_sin (T, R) float64, ∈ [−1, 1].** Trig pairs.
  *Primary:* enabler for circular mean / variance; probably not
  shown directly. An adviser cares about the angle, not its
  sine.
  *Useful:* 2D (cos χ, sin χ) scatter on the unit circle —
  equivalent-looking to a radial histogram, better for dense data.

---

# /dssp — secondary structure and H-bond assignment

- **`ss8 (T, R) int8.** 8-class DSSP code per residue per frame.
  Encoding in the group attribute `ss8_encoding`:
  `0=H(α), 1=G(3₁₀), 2=I(π), 3=E(strand), 4=B(bridge), 5=T(turn),
  6=S(bend), 7=C(coil)`.
  *Primary:* **per-residue heatmap over time** — residue on one
  axis, frame_time on the other, cell coloured by 8-class code.
  Adviser sees secondary-structure evolution immediately (helix
  unfolding, β-strand forming). DSSP already drives ribbon colour
  per frame; the heatmap adds the time dimension.
  *Useful:* fraction-helix / fraction-strand as a scalar line
  chart.
- **`hbond_energy (T, R) float64, kcal/mol.** Strongest acceptor
  energy at this residue (Kabsch-Sander electrostatic model).
  *Primary:* **ribbon thickness modulation** — sausage
  representation with ribbon width ∝ |hbond_energy|. Reviewers
  see H-bond strength in structural context directly on the
  protein. *Useful:* residue × frame heatmap (diverging Moreland)
  for the time story; per-residue line chart on an individual
  residue selection.

---

# /predictions — ridge and MLP shielding predictions

Output of the calibrated models (ridge from Stage 1, gated MLP for
N and C). Every field here is **a model prediction**, not a
physical observable. The thesis's discipline is that these are
diagnostics for whether the kernel set carries the signal; they
are not shielding.

- **`raw_T0 (T, N) float64, ppm.** Ridge prediction on
  unnormalised kernels. Isotropic.
  *Primary:* scalar line chart paired with the sum of measured
  T0 contributions. The **difference** is the ridge residual,
  which is visually more informative than either alone.
  *Useful:* colour bubble (diverging Moreland).
- **`raw_T2 (T, N, 5) float64, ppm.** Ridge prediction T2 (5
  components).
  *Primary:* tensor glyph using `raw_T2` as the T2 part (same
  renderer as measured shielding). Colour by a neutral/prediction
  fourth colour to distinguish from measured. *Useful:* |T2|
  scalar line chart.
- **`norm_T0 (T, N), norm_T2 (T, N, 5) float64, ppm.** Ridge
  predictions on per-frame-per-protein normalised kernels.
  *Primary / useful:* same modalities as `raw_*`, on a second
  toggle.
- **`raw_ring_current, raw_efg, raw_bond_aniso, raw_quadrupole,
  raw_dispersion, raw_hbond (T, N, 5) float64, ppm.** Per-group
  additive contributions to the ridge prediction.
  *Primary:* **stacked bar chart per atom per frame** — six
  bars side-by-side, positive and negative stacked, showing
  cancellation vs reinforcement. This is the thesis's "additive
  decomposition" picture; it belongs in a dedicated dock.
  *Useful:* per-group |T2| scalar line charts stacked vertically.
- **`norm_ring_current, norm_efg, norm_bond_aniso, norm_quadrupole,
  norm_dispersion, norm_hbond (T, N, 5) float64, ppm.** Same
  decomposition on normalised kernels.
  *Primary / useful:* same modalities as the raw counterparts,
  selector toggle.
- **`raw_water, raw_charges, raw_sasa (T, N, 5) float64, ppm.**
  New-feature groups; zero weights in Stage 1 but ride along so
  the raw kernel contribution is visible. Phase 2 calibration
  may assign non-zero weights.
  *Primary:* stacked bar chart alongside the core groups; zero
  contribution is itself informative.
- **`norm_water, norm_charges, norm_sasa (T, N, 5) float64, ppm.**
  Normalised versions.
- **`mlp_T0 (T, N), mlp_T2 (T, N, 5) float64, ppm.** Gated-MLP
  nonlinear prediction. For N and C, captures the +0.13 / +0.17
  R² above ridge.
  *Primary:* **ridge-MLP difference** rendered as a tensor glyph
  (δ_ridge, δ_mlp, δ_diff), showing where the nonlinear signal
  lives. Per-element selector. *Useful:* scalar line chart of
  `mlp_T0 − raw_T0` as "nonlinear residual".
- **`mlp_kernel_weights (T, N, 55) float64, dimensionless.**
  Per-kernel gating weights the MLP chose. 55 kernels total.
  *Primary:* **per-atom horizontal heat-strip** — one column per
  kernel (labels from the kernel manifest), colour = gate weight.
  Dock pattern: same as the K=6 bar chart but with 55 slots. The
  adviser sees which kernels the MLP is turning on or off at
  this atom at this frame. *Useful:* protein-wide atlas
  (per-element mean weight per kernel) as a reference.

---

# /projections — mathematical projections and distance metrics

- **`subspace_coords (T, N, 20) float64, dimensionless.**
  Per-element PCA subspace coordinates of the normalised-kernel
  vector. Padded to 20; actual dimension per element in
  `subspace_n_dims`.
  *Primary:* **2D projection dock** — pick two coordinates (via
  dropdown, colour by third if desired), show the atom's
  trajectory through subspace. One scatter plot per element so
  the element-specific subspaces (H=20, C=6, N=3, O=12) are
  respected. *Useful:* leading-coordinate line chart per atom;
  colour bubble by leading coord magnitude.
- **`subspace_n_dims (N,) int8.** Actual dims per atom's
  element. Static; shown in inspector + used to mask padding.
- **`residual_water_efg, residual_water_first, residual_aimnet2_efg,
  residual_hydration, residual_charges (T, N) float64, ppm.**
  Projection of a new feature onto the ridge residual direction
  — how much of the unexplained variance does this feature
  explain?
  *Primary:* **SAS surface colouring** (diverging Moreland,
  centred at 0) — where on the protein does this new feature fill
  the residual? Positive patches are calibration candidates.
  *Useful:* scalar line chart to watch excursions; atom colour
  bubble fallback.
- **`mahalanobis (T, N) float64, dimensionless σ.** Mahalanobis
  distance of the atom's kernel vector from the 720-protein
  calibration distribution (per-element mean + inverse covariance).
  *Primary:* **SAS surface colouring**, sequential colormap capped
  at a 3 σ outlier threshold — "which regions is this trajectory
  exploring uncalibrated geometry for?" reads as a continuous
  surface pattern far better than a sea of coloured atoms, and
  those surface regions are the DFT-candidate targets. *Useful:*
  ribbon thickness modulation if the question is per-residue;
  atom colour bubble as an atom-count alternative; scalar line
  chart to watch excursions over time.
- **`mlp_gate_weights (T, N, 55) float64.** Copy of
  `/predictions/mlp_kernel_weights` for consumer convenience.
  Same display as there.

---

# /events — reserved for Phase 3

Schema open until exploratory analysis reveals which discrete
events (rotamer flips, H-bond switches, water-shell exchanges,
charge jumps) matter. No fields defined yet; the reader should
treat `/events` as optional and silently absent today.

---

# Visualisation reference

This section names the four primitives the glossary's "Primary"
calls refer to. Each is a class the reader either has or should
add to serve the recommendations. Existing classes are credited;
new ones are sketched against the overlay contract
(`MoleculeScene` → `setFrame`).

## Superquadric tensor glyph on atom

**What.** A superquadric surface rendered at the atom position,
oriented by the T2 eigenvectors, sized by eigenvalues, coloured
by T0 (for shielding) or η (for EFG). Replaces the ambiguous
ellipsoid glyph; replaces cuboid-style cubes which mis-orient on
rotationally symmetric tensors (Kindlmann & Westin, 2004;
Schultz & Kindlmann, TVCG 2010).

**Why.** Ellipsoids look the same in profile for shapes of
different anisotropy. Superquadrics do not — oblate vs prolate
vs spherical are visually distinct, sign of eigenvalue is
distinct, and the glyph is symmetric under the tensor's symmetries.
Every symmetric second-order tensor visualisation since 2010 has
moved to superquadrics; we should not regress.

**How to implement.** VTK has `vtkTensorGlyph`; with the
superquadric source from
`vtkSuperquadricTensorGlyphFilter` (third-party, small, MIT-licensed)
or implemented directly via the Kindlmann convention (γ parameter
tuned to 2.5 for anisotropy discrimination). Overlay class
`QtTensorGlyphOverlay`, pattern `QtRingPolygonOverlay`. Input:
per-atom (T2 components, T0 for colour, scale cap). One actor
per atom for the full protein is feasible (N ≤ 5000).

**For asymmetric shielding tensors (BS / HM / McConnell / RS /
HBond).** T2 drives the glyph shape; T0 drives colour; T1 is
rendered as a small perpendicular arrow through the glyph,
toggleable. The asymmetric piece is small in magnitude for most
atoms but the T1 arrow makes it visible when it matters.

**For symmetric traceless EFGs (Coulomb / APBS / water /
AIMNet2).** Haeberlen axes: V_zz as a long double-ended arrow
through the atom, V_yy / V_xx as orthogonal short arrows or an
oriented disk with aspect ratio = (1−η). This is the standard
NMR crystallography rendering and the one an adviser with an
NMR background will recognise on sight.

## Vector arrow on atom

**What.** An arrow glyph at the atom position. Length ∝ magnitude
(clamped at a maximum for legibility), colour by a per-field hue.
For E-field, B-field, dipole, nearest-partner directions.

**How.** `vtkGlyph3D` with `vtkArrowSource`, coloured via a
`vtkColorTransferFunction` on magnitude. Overlay class
`QtAtomArrowOverlay`, selector for which vector field to show;
per-atom filter so only picked/hovered atoms render the arrow if
the protein is large.

## Atom colour bubble

**What.** The `vtkMolecule` atoms recoloured from the default
element colour to a scalar colormap plus a `vtkScalarBarActor`
legend. Per-atom scalar array drives the lookup.

**How.** Extend `MoleculeScene` so its mapper can consume a
per-atom scalar array. Overlay class `QtAtomColourOverlay` holds
a dropdown: element | any per-atom scalar from the H5 | derived
diagnostic (ridge residual, AIMNet2 − ff14SB charge, Mahalanobis).
Colormap defaults: viridis for non-negative, Moreland diverging
for signed. Honour the standard scales (e.g. SASA on [0, 80]
Å²; Mahalanobis on [0, 5] σ).

## SAS surface colouring

**What.** A solvent-accessible-surface mesh around the protein,
coloured by a per-atom or per-surface-point scalar. A continuous
coloured surface that reads protein-wide patterns at a glance.

**Why.** Atom colour bubbles work on small proteins with clear
contrasts; they degrade fast at 5000 atoms or when the signal
is subtle. A continuous SAS mesh carries the protein-wide story
in the form reviewers already read from PyMOL / Chimera /
ChimeraX figures. For Mahalanobis distance, ridge residual,
hydration asymmetry, and per-residue strain summed onto surface
atoms, the surface is the right canvas.

**How to implement.** VTK has two paths:
`vtkMolecularSurfaceFilter` (Connolly / SES-style, slower) or
`vtkMarchingCubes` over a Shrake-Rupley density grid (faster,
matches `SasaResult`). Per-frame remeshing is affordable (~100 ms
for a 2000-atom protein); per-atom scalar mapping via a
`vtkFloatArray` on the mesh points, interpolated to surface
vertices by nearest atom or inverse-distance-weighted. Overlay
class `QtSasSurfaceOverlay`, scalar-selector dropdown shared with
the colour-bubble overlay, opacity slider so the underlying
molecule stays visible.

## Ribbon thickness modulation (sausage)

**What.** The backbone ribbon with per-residue width driven by a
scalar array. Thickness encodes magnitude, existing ribbon colour
(by residue type or DSSP) encodes a second channel — two
channels for free.

**Why.** Per-residue scalars (DSSP H-bond energy, summed
bonded-energy, per-residue Mahalanobis, chemical-shift
perturbation) want a residue-scoped render. The ribbon is drawn
anyway; widening it for strong-signal regions is the cleanest
visual hierarchy from a distance and is the classic NMR-review
trick for chemical shift perturbation.

**How to implement.** `vtkProteinRibbonFilter` exposes per-CA
width via a scalar array on the ribbon polydata. Overlay class
`QtRibbonThicknessOverlay` reads a per-residue array (either an
H5 dataset directly or a derived per-residue summary from
per-atom fields), maps to a width range (min ≈ 0.3×, max ≈ 2.0×
default), writes to the ribbon source before re-running the
filter. Independent of the existing ribbon colour-by-SS mode so
the two channels stay separable.

## Protein-wide volumetric isosurface

**What.** A 3D scalar field over the protein bounding box,
sampled on a grid and contoured at one or more values. For the
total B-field this is the cumulative magnetic-environment
picture — the per-ring butterfly is the same primitive at a
narrower scale.

**Why.** Per-atom vector arrows answer "what does atom i see";
the volumetric isosurface answers "what does the magnetic
environment of this fold look like". Both are interesting;
reviewers with NMR backgrounds expect the second from textbooks
(Haigh-Mallion / Johnson-Bovey picture around a ring system).

**How to implement.** Protein-wide grid (e.g. 64³ voxels over
the bounding box, adjustable), per-voxel scalar from the same
evaluators that drive the per-ring butterfly
(`QtBiotSavartCalc::EvaluateBField`, summed over every ring).
`vtkContourFilter` at ±threshold, two actors (shielded /
deshielded, sky-blue / coral), translucency for depth cues.
Overlay class `QtProteinFieldOverlay`, same pattern as
`QtFieldGridOverlay` with a wider grid and summed sources.
Cost: ~100–300 ms per frame at 64³ × ~20 rings; cache per-ring
base grids and sum into the shared grid if this is tight.

## Projection / multi-field charts in a dedicated dock

**What.** Any view that does not project to a single atom: UMAP
of the AIMNet2 embedding, Ramachandran scatter, DSSP heatmap,
K=6 ring bar chart, per-group stacked bar chart for predictions,
system-energy line chart.

**How.** New `QDockWidget` per view, each with its own paint
surface. Projection dock uses a precomputed UMAP (one-off at
file-load, ~10 s). Ramachandran uses a `QChart` with a scatter
series. DSSP heatmap uses a `QGraphicsView` grid or a custom
`QQuickItem`. Keep them separable so an adviser can open one
dock without loading the others.

## Current scalar line chart

**What.** `QtAtomTimeSeriesDock` as it stands 2026-04-17. One
per-atom scalar vs frame_time, cursor at current frame. Desc
table for field selection. See `notes/TIME_SERIES_EXPANSION.md`
for the next scope-expansion step.

---

# Index — "what should this field look like?"

The index groups the ~200 fields by *primary* modality so a
reader planning the next overlay knows what surfaces need
building first. Only the modality marked Primary is listed; a
field with useful secondary modalities is covered once in its
own group section above.

## Superquadric tensor glyph at atom (new overlay)

Per-atom rank-2 tensors. All require the `QtTensorGlyphOverlay`
class described above.

- `/ring_current/bs_shielding`, `hm_shielding`, `rs_shielding`
- `/bond_aniso/mc_shielding` and its five sub-decompositions
  (backbone, sidechain, aromatic, CO_nearest, CN_nearest)
- `/quadrupole/pq_shielding`
- `/dispersion/disp_shielding`
- `/hbond/hbond_shielding`, `nearest_spherical`
- `/efg/coulomb_shielding`, `/efg/aimnet2_shielding`
- `/predictions/raw_T2`, `/predictions/mlp_T2` (as prediction
  colour)

## Haeberlen principal-axis glyph at atom (new overlay, same
  class with EFG mode)

Traceless symmetric EFGs. V_zz arrow + η disk.

- `/efg/coulomb_total`, `coulomb_backbone`, `coulomb_aromatic`
- `/efg/apbs_efg`
- `/efg/aimnet2_total`, `aimnet2_backbone`, `aimnet2_aromatic`
- `/water/efg`, `/water/efg_first`

## Vector arrow at atom (new overlay)

Per-atom 3-vectors.

- `/efg/E_total`, `E_backbone`, `E_sidechain`, `E_aromatic`,
  `E_solvent`, `apbs_efield`
- `/water/efield`, `efield_first`, `dipole_vector`,
  `surface_normal`
- `/bond_aniso/dir_nearest_CO` (line segment, not magnitude)
- `/hbond/nearest_dir` (line segment, not magnitude)
- `/sasa/normal`
- `/ring_current/total_B_field` (secondary — volumetric
  isosurface is primary)

## Atom colour bubble (new overlay)

Per-atom scalars where the protein-wide picture is the point
*and* the protein is small enough / the contrast stark enough
that a sea of coloured spheres reads cleanly. For the continuous-
surface version, see SAS surface colouring below.

- `/atoms/eneg_sum_1`, `eneg_sum_2`, `bfs_decay`,
  `partial_charge`
- `/ring_current/G_iso_var_8A`, `n_rings_{3A,5A,8A,12A}`,
  `mean_ring_dist`, `nearest_ring_atom`
- `/efg/E_magnitude`, `E_bond_proj`, `E_backbone_frac`
- `/bond_aniso/co_sum`, `cn_sum`, `sidechain_sum`, `aromatic_sum`,
  `co_nearest`, `nearest_CO_dist`, `nearest_CN_dist`
- `/hbond/nearest_dist`, `inv_d3`
- `/sasa/sasa`
- `/water/n_first`, `n_second`, `half_shell_asymmetry`,
  `dipole_cos`, `nearest_ion_dist`, `sasa_asymmetry`,
  `sasa_first_shell_n`
- `/charges/aimnet2_charge`, `eeq_charge`, `eeq_cn`
- `/bonded_energy/bond`, `angle`, `urey_bradley`, `proper_dih`,
  `improper_dih`, `cmap`

## SAS surface colouring (new overlay)

Per-atom scalars where the protein-wide spatial pattern is the
story and a continuous coloured surface beats a sea of spheres.

- `/projections/mahalanobis`
- `/projections/residual_water_efg`, `residual_water_first`,
  `residual_aimnet2_efg`, `residual_hydration`, `residual_charges`
- `/water/sasa_dipole_align`, `sasa_dipole_cohere`
- `/bonded_energy/total` (protein-wide strain map)

Derived (not raw H5, painted onto the surface):

- Ridge residual: summed measured T0 minus `/predictions/raw_T0`
- Per-residue chemical-shift perturbation (WT vs mutant)

## Ribbon thickness modulation (new overlay)

Per-residue scalars where the structural context is the story.

- `/dssp/hbond_energy`

Derived:

- Per-residue sum of `/bonded_energy/total`
- Per-residue mean of `/projections/mahalanobis`
- Per-residue chemical-shift perturbation

## Protein-wide volumetric isosurface (new overlay)

3D scalar fields over the bounding box.

- `/ring_current/total_B_field` — |B| isosurface contoured at
  ±0.1 mT, shielded (sky blue) / deshielded (coral)

## Projection / n-D chart (new dock)

Views that are not atom-local.

- `/aimnet2_embedding/aim` — UMAP/PCA projection dock
- `/projections/subspace_coords` — per-element 2D scatter dock
- `/per_ring/*` — K=6 bar chart dock, per-atom per-calculator
- `/ring_current/bs_T0_per_type`, `hm_T0_per_type`,
  `/quadrupole/pq_T0_per_type`, `/dispersion/disp_T0_per_type`
  — per-type bar charts per atom
- `/predictions/raw_{ring_current,efg,bond_aniso,quadrupole,
  dispersion,hbond}`, and norm_*, water/charges/sasa
  counterparts — stacked-bar additive-decomposition dock
- `/predictions/mlp_kernel_weights`, `/projections/mlp_gate_weights`
  — 55-kernel heat-strip dock
- `/dihedrals/phi`, `psi` — Ramachandran dock
- `/dihedrals/chi1..chi4` (with `chi*_{cos,sin}`) — circular
  histogram dock per residue
- `/dssp/ss8` — residue × time heatmap dock
- `/energy/*` (except `box`, `virial`, `pressure_tensor`) —
  system-energy dock
- `/energy/virial`, `/energy/pressure_tensor` — 3×3 matrix view
  with diagonal/off-diagonal split
- `/ring_geometry/data` — ring-scoped time series (scope
  expansion)

## Scalar line chart (existing `QtAtomTimeSeriesDock`)

Everything else suitable for a strip chart. Fields already in
the current Desc table are noted in
`notes/TIME_SERIES_EXPANSION.md`. Extending the Desc table is a
one-line-per-field operation where the `QtFrame` accessor exists,
and a ~10-line operation where it doesn't.

Primary modality for a line chart applies when the field is an
atom-local, per-frame scalar **and** the relevant story is its
variability over time (not the protein-wide pattern at a moment).
For scalars where the protein-wide snapshot is the story, the
primary is the colour bubble.

---

# Surfacing discipline

Three rules the reader should enforce when a new modality is
wired up.

1. **Field name on screen.** Inspector rows and chart titles
   must show the H5 dataset path (e.g. `/efg/apbs_efg`) alongside
   the human label. Advisers reading a paper will go look at the
   H5; the group path is the contract (feedback_ui_data).
2. **Units visible always.** Every axis label, every colour bar,
   every glyph-scale legend prints units. No bare numbers. ppm,
   Tesla, V/Å, V/Å², Å, rad, kJ/mol, kcal/mol (DSSP only).
3. **T2 is sacred.** A scalar summary of a tensor (T0, |T2|, ‖T‖)
   may appear *in addition to* a tensor glyph, never in place of
   it. A field marked in this glossary as tensor-glyph-primary
   that is shown only as a scalar is a regression
   (feedback_t2_sacred).

Sources informing best-case recommendations:

- [Kindlmann, G. (2004). Superquadric Tensor Glyphs. VisSym
  2004](https://cgl.ethz.ch/teaching/scivis_common/Literature/kindlmann04.pdf)
- [Schultz, T.; Kindlmann, G. (2010). Superquadric Glyphs for
  Symmetric Second-Order Tensors. IEEE TVCG](http://people.cs.uchicago.edu/~glk/sqd/schultzTVCG10SuperquadricTensorGlyphs.pdf)
- [Chemical Shift Tensor Conventions — Haeberlen / Mehring /
  Maryland (Eichele, CCP14)](http://www.ccp14.ac.uk/ccp/web-mirrors/klaus_eichele_software/klaus/nmr/conventions/csa/csa.html)
- [Further Conventions for NMR Shielding and Chemical Shifts
  (IUPAC Recommendations 2008)](https://bmrb.io/standards/iupac_2008.pdf)
- [Harris et al. TensorView: A software tool for displaying NMR
  tensors. Conc Magn Reson A (2018)](https://par.nsf.gov/servlets/purl/10099402)
- [Sun et al. (2020). Density functional theory-based electric
  field gradient database. Sci Data 7, 314.](https://www.nature.com/articles/s41597-020-00707-8)
