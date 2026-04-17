# What the Extraction Computes and How We Look At It

The extraction engine computes ~100 physical quantities per atom per
frame. This document organises them by what they measure and what
questions they answer when combined with experimental RefDB shifts
across 685 proteins.

---

## I. Ring current fields (4 independent models of the same physics)

### Biot-Savart (BS)
- **What it computes:** Magnetic field from a circular current loop
  at each aromatic ring. T0 is isotropic ring current shift. T2 is
  the directional structure.
- **Per-ring-type:** PHE, TYR, TRP_benzene, TRP_pyrrole,
  TRP_perimeter, HIS/HID/HIE. Eight separate T0 and T2 values.
- **Per-ring sparse:** Individual ring's BS contribution with
  distance, cylindrical coordinates (rho, z, theta), McConnell
  factor, all in ring_contributions.
- **Distance dependence:** r^-3 (verified at -3.04 in calibration).
- **What to look for in the ensemble:** Which rings dominate at each
  atom? Does the dominant ring change across conformations? Do
  stacked rings (multiple rings contributing) show cooperative
  effects?

### Haigh-Mallion (HM)
- **What it computes:** Surface-integral ring current model. Uses
  actual ring geometry (dihedral angles between ring bonds and
  atom vector), not circular approximation.
- **Why both BS and HM:** They agree to 1% at >5A (dipole limit)
  but differ in near-field. The DIVERGENCE between BS and HM at
  short range is itself a signal -- it measures how much the
  non-dipolar character of the ring matters for that atom.
- **What to look for:** Atoms where BS and HM disagree across the
  ensemble. The disagreement pattern encodes ring geometry
  sensitivity. HM should be better for HIE (asymmetric 5-ring).

### Pi-Quadrupole (PQ)
- **What it computes:** Electric quadrupole interaction from the
  pi-electron cloud. Different physics from ring current (electric
  not magnetic). Distance dependence r^-5 (faster falloff).
- **What to look for:** PQ dominates at very short range but dies
  fast. In the ensemble, PQ variance should be highest for atoms
  that move in and out of the near-field. PQ vs BS ratio across
  frames measures the short-range/long-range balance.

### Dispersion (Disp)
- **What it computes:** van der Waals contact tensor (r^-8 T2,
  r^-6 scalar). The scalar and tensor are DIFFERENT quantities
  (not T0/T2 decomposition of the same tensor).
- **What to look for:** Dispersion is the shortest-range kernel.
  It's effectively a contact indicator. In the ensemble, disp
  variance tracks ring-atom contact dynamics.

### Ring Susceptibility (RingSusc)
- **What it computes:** Diamagnetic susceptibility anisotropy of
  the aromatic ring. A different angular dependence from BS.
- **What to look for:** How it correlates with BS across the
  ensemble. If perfectly correlated, it's redundant. If the
  correlation breaks in specific geometries, the susceptibility
  model captures something BS misses.

### The three-view convergence/divergence signal

BS, HM, and EFG_aro all measure the same geometric relationship
(R^2 ~0.77 for H, cos ~0.93 with DFT target). In the far field
(>5A) they must agree -- all three approach the same dipole limit.
In the near field they can diverge because they encode the ring
geometry differently:

- BS assumes circular current loop -- fails for asymmetric rings
  (HIE calibration R^2 = 0.062)
- HM uses actual ring bond geometry -- handles asymmetry better
- EFG uses point charges at atom positions -- no current model at all

Across the ensemble, tracking WHERE and WHEN these three views
diverge is a primary investigation target. The divergence pattern
tells us:
- Which ring types break the circular-loop approximation
- Which atom positions are sensitive to model choice
- Which conformations push atoms into the near-field regime

This is not a model limitation to work around. It is a physics
finding about the validity range of the ring current models.

---

## II. Electrostatic fields (a third geometric view)

### The key finding: EFG is geometric, not charge-dependent

The calibration showed BS, HM, and EFG_aro all score R^2 ~0.77 for
H independently (0.772, 0.784, 0.766). They are three views of the
same ring geometry through different physics:

- BS: magnetic dipole field from ring current
- HM: surface integral of ring current
- EFG: electric field gradient from ring atom charges

The angular correlation with the DFT target is:
- EFG_aro (topology charges): cos = 0.9296
- MopacEFG_aro (MOPAC charges): cos = 0.9323

Difference: 0.003. The EFG tensor direction comes from the geometric
T-tensor (3 r_a r_b - r^2 delta_ab)/r^5, which depends only on atom
POSITIONS, not charges. Charge quality is irrelevant to the angular
structure. The forward selection artefact (MOPAC enters first, topology
looks useless) obscured this in the old analysis.

### Coulomb E-field and EFG
- **What it computes:** Electric field (E) and electric field
  gradient (EFG) at each atom from all other atoms' partial charges.
  Decomposed: backbone contribution vs aromatic contribution.
- **Aromatic EFG:** A geometric proxy for ring current angular
  structure. The ring atoms that carry charge are the same atoms
  that carry ring current. The T-tensor from those atoms captures
  the same spatial relationship.
- **Backbone EFG:** Changes with backbone conformation (phi/psi).
  This is NOT a ring current proxy -- it's independent physics
  from backbone charge rearrangement. Dominant for C shifts.
- **What to look for in the ensemble:**
  - EFG_aro should co-vary tightly with BS/HM (same geometry).
    Where it DOESN'T is where charge asymmetry matters.
  - EFG_bb variance across frames, stratified by element. C and N
    should show high EFG_bb variance (backbone conformation
    dominates their shifts). H should show lower EFG_bb variance.
  - The element-dependent variance pattern IS the physics.

---

## III. Bond anisotropy (McConnell)

### What it computes
Magnetic shielding from the anisotropy of nearby chemical bonds.
Each bond acts as a magnetic dipole. Decomposed into 5 categories:
- Peptide C=O (strongest, partial double bond)
- Peptide C-N (partial double bond)
- Backbone other
- Sidechain C=O
- Aromatic bonds

### What to look for
- Peptide bond categories should be conformationally STABLE (the
  peptide bond is rigid). If they vary, it's because the RELATIVE
  geometry changes, not the bond itself.
- Aromatic bond category should track ring motion (same structural
  modes as BS/HM).
- Cross-correlation between McConnell aromatic and BS: both respond
  to ring motion but through different physics (magnetic dipole of
  bonds vs ring current). The residual between them is the
  multi-body effect BS captures and McConnell doesn't.

---

## IV. H-bond geometry

### What it computes
Nearest H-bond distance, count within 3.5A, H-bond shielding tensor.

### What to look for
- H-bond formation/breaking across the ensemble. Binary signal:
  H-bond present or absent.
- HN shift residuals should correlate with H-bond occupancy.
  Literature: HN shifts move ~1 ppm per H-bond (Wagner et al. 1983).
- H-bond dynamics vs ring current dynamics: do they correlate
  (same structural motion) or are they independent?

---

## V. Backbone geometry (DSSP)

### What it computes
Phi, psi, SASA, secondary structure classification per frame.

### What to look for
- DSSP is the CONTEXT for every other kernel. Stratify all kernel
  analyses by secondary structure (helix/sheet/loop).
- Phi/psi directly modulate backbone shifts. The centring in Layer 2
  removes the mean per residue type but not the phi/psi effect
  within a residue type. This is the dominant confound for CA/N
  analysis.
- SASA: surface exposure changes which kernels are relevant. Buried
  atoms see more ring current. Surface atoms see more solvent.

---

## VI. Per-ring sparse data (the richest view)

### What it gives us
59 columns per (atom, ring) pair. For each atom near each ring:
- Distance, cylindrical coordinates (rho, z, theta)
- McConnell factor ((3cos^2(theta)-1)/r^3)
- Every ring kernel's individual contribution (BS, HM, PQ, Chi,
  Disp) as 9-component tensors
- Azimuthal angle (cos_phi, sin_phi) relative to ring vertex 0
- Gaussian density, exponential decay

### Why this is the richest view
Everything else is a SUM over rings. This is the individual terms.
Across the ensemble, we can track:
- How each ring-atom pair's geometry evolves
- Which pairs are stable vs dynamic
- Whether individual ring contributions correlate with each other
  (cooperative ring currents from stacked rings)
- The azimuthal dependence (phi) which the calibration never used

### What to look for
- Pairs where distance crosses a threshold (gating)
- Pairs where theta flips sign (atom crosses ring plane)
- Correlation between different rings' contributions at the same
  atom (multi-ring effects)
- Azimuthal angle dependence -- does cos_phi/sin_phi correlate with
  the shift residual beyond what distance explains?

---

## VII. Things we compute that the calibration never used

These are available in the extraction but were not part of the
55-kernel calibration:

1. **T1 components** of every kernel. Encodes geometric asymmetry
   of each kernel-atom interaction. 3 extra values per kernel.
   The calibration discarded T1. The ensemble analysis can look at
   T1 variance and T1-shift correlation.

2. **BS total B-field vector** (3 components). The actual magnetic
   field, not the shielding tensor derived from it. Magnitude and
   direction.

3. **BS ring_counts** at 4 distance thresholds (3/5/8/12A). A
   coarse-grained ring proximity measure.

4. **McConnell nearest CO/CN distance and direction**. Not just the
   tensor but the geometry of the dominant bond.

5. **Coulomb scalars**: E magnitude, bond projection, backbone
   fraction, aromatic E fraction. Decomposition of the electrostatic
   environment.

6. **H-bond scalars**: nearest distance, 1/r^3 factor, count.

7. **Azimuthal angle** in ring_contributions. The phi angle around
   the ring axis, relative to vertex 0. This breaks the cylindrical
   symmetry of the BS model and could reveal non-circular-loop
   effects.

Each of these is a potential finding: "T1 of BS correlates with
shift residual in a way that T2 alone doesn't capture" would be
a new observation about asymmetric ring current effects.

---

## How this connects to the RefDB ground truth

For each matched atom, the experimental shift residual (after
residue-type centering) is the TARGET. Every quantity above is a
MEASUREMENT. The investigation is: which measurements explain which
part of the target, for which elements, at which distances, in which
structural contexts?

The 28 calibration realities told us the instrument is accurate.
The ensemble analysis tells us what the instrument sees in nature.
