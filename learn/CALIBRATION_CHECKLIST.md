# Calibration Checklist: Every Assumed Constant

Every hardcoded value in the 8 calculators, where it comes from,
what might modulate it, and how to check. Not ML feature hunting --
specific physics questions with specific measurements.

Last updated: 2026-04-05.

---

## Distance Boundaries (when to stop computing)

### [ ] RING_CALC_CUTOFF = 15.0 Å
- **Source**: Literature convention (spherical)
- **Physics**: The iso-|B| surface from a magnetic dipole is NOT
  spherical. It follows a sextic surface: r^6 = C^6 × (3z² + r²).
  Polar-to-equatorial aspect ratio = 4^(1/6) = 1.260.
- **What modulates**: Nothing -- pure geometry of the dipole field.
  But the SHAPE of the cutoff matters: spherical cutoff at 15 Å
  either wastes computation on weak in-plane contributions or misses
  strong axial contributions at the same distance.
- **Sweep**: 3 shapes × 8 scale values = 48 configurations.
  Sphere, sextic, ellipsoid. Compare residuals.
- **MOPAC relevance**: None directly. This is pure field geometry.

### [ ] McConnell bond cutoff = 10.0 Å (via spatial index)
- **Source**: Convention. 1/r^3 decay makes 10 Å contribution ~0.1%
  of 3 Å.
- **Physics**: Backbone amides at 10 Å still feel C=O anisotropy.
  The question is whether including them helps or adds noise.
- **Sweep**: 5, 7, 10, 15, 20 Å. Plot residual vs cutoff.
- **MOPAC relevance**: Bond orders at different distances could show
  whether distant bonds have meaningfully different anisotropy.

### [ ] DISP_VERTEX_R_CUT = 5.0 Å, R_SWITCH = 4.3 Å
- **Source**: CHARMM convention (Brooks et al. 1983).
- **Physics**: 1/r^6 means 5 Å is genuinely negligible. Low priority.
- **Sweep**: Maybe. Dispersion is a small effect.

### [ ] HBOND_MAX_DIST = 50.0 Å
- **Source**: Effectively infinite (H-bond geometry from DSSP has its
  own criteria).
- **Physics**: Not a real cutoff. DSSP identifies H-bonds; this just
  prevents numerical issues.

---

## Filter Parameters (when to reject a kernel evaluation)

### [ ] DipolarNearField extent_factor = 0.5 × source_extent
- **Source**: Physical reasoning -- multipole expansion invalid inside
  source distribution.
- **Physics**: "Inside" depends on the electron distribution, not just
  the geometric extent. For a PHE ring with radius 1.4 Å, the
  exclusion zone is 1.4 Å from center. But π-electrons extend above
  and below -- the exclusion should be shaped, not spherical.
- **Sweep**: 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5.
  **Prime sweep target.** Expect a clear minimum.
- **MOPAC relevance**: Orbital populations (s-pop, p-pop) on ring
  atoms indicate electron distribution extent. If p-orbital population
  differs between ring types, the near-field extent should differ too.

### [ ] Sequential exclusion min_sep = 2 residues
- **Source**: Convention. Through-bond coupling not modelled by
  through-space kernel.
- **Physics**: Residues i±1 are connected by peptide bonds whose
  through-bond effect dominates the through-space dipolar kernel.
  But residues i±2 have real through-space contributions. Is 2 right?
- **Sweep**: 1, 2, 3, 4. Small sweep.
- **MOPAC relevance**: Bond orders across the i→i+1 peptide bond
  measure how strongly coupled adjacent residues are. If bond order
  drops sharply at i+2, the exclusion is justified.

### [ ] RingBondedExclusion (topology-based)
- **Source**: Physics. Atoms that are ring vertices or bonded to them
  experience through-bond effects, not through-space dipolar fields.
- **Physics**: Correct by construction. The TRP5+TRP6=TRP9 additivity
  test confirms this filter eliminates the surface integral excess.
- **Sweep**: Not needed. Binary: on or off. On is correct.

---

## Geometric Parameters (shape of the source)

### [ ] BS lobe_offset per ring type
- **Source**: Haigh & Mallion 1979 (via Johnson-Bovey model).
  PHE/TYR/TRP6 = 0.64, TRP5 = 0.52, HIS = 0.50, TRP9 = 0.60 Å.
- **Physics**: Distance of π-electron current loops above/below ring
  plane. Determines the B-field shape: larger offset → sharper lobes.
  T0 constrains intensity. T2 constrains lobe offset.
- **Directly learnable**: The linear fit (Level A) learns intensity.
  The T2 residual after fitting constrains lobe offset because
  changing d changes the angular pattern.
- **MOPAC relevance**: **Directly relevant.** p-orbital population on
  ring carbons (from mopac_scalars.npy) indicates how much electron
  density is above/below the ring plane. If TRP5 has different p-pop
  than PHE, the lobe offset should differ, and MOPAC measures this
  independently from the shielding fit.
- **Sweep**: Per ring type. 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 Å.

### [ ] McConnell dipole centre = bond geometric midpoint
- **Source**: Convention. The effective magnetic dipole is assumed at
  the bond midpoint.
- **Physics**: Electron density in C=O is asymmetric (shifted toward O
  due to electronegativity). The effective dipole centre should be
  shifted along the bond axis.
- **MOPAC relevance**: **Directly relevant.** Bond order from MOPAC
  gives a continuous measure of electron distribution asymmetry.
  A C=O bond with order 1.75 vs 1.60 has different electron density
  centre. Delta(bond_order) between WT and ALA traces redistribution.
  The midpoint_shift parameter should correlate with bond order
  asymmetry.
- **Coupled constraint**: midpoint_shift must improve T0, T1, AND T2
  simultaneously. If it helps T0 but ruins T2, the model rejects it.
- **Sweep**: shift = 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 (fraction toward
  more electronegative atom).

### [ ] HM quadrature thresholds: L1 = 2.0 Å, L2 = 1.0 Å
- **Source**: Numerical accuracy (adaptive quadrature).
- **Physics**: Not physics. Integration precision. Verify convergence
  but don't sweep.

---

## Electrostatic Assumptions

### [ ] Coulomb charges = ff14SB fixed per atom type
- **Source**: AMBER force field.
- **Physics**: ff14SB charges do NOT respond to conformation or
  local environment. Every GLY CA has the same charge regardless
  of what's nearby.
- **MOPAC relevance**: **Core contribution.** MOPAC gives QM-derived
  charges that respond to the electronic environment. Recalculating
  Coulomb EFG with MOPAC charges and comparing to ff14SB tells you
  whether charge polarisation matters for ring current physics.
  If delta(EFG_mopac) correlates with the DFT T2 residual better
  than delta(EFG_ff14sb), charge polarisation is load-bearing.

### [ ] Bond anisotropy Δχ = one value per bond category
- **Source**: Literature values per category (CO, CN, aromatic, etc.).
- **Physics**: Δχ should depend on the actual bond character, not
  just the category. A C=O with bond order 1.8 has different
  anisotropy than one with 1.6.
- **MOPAC relevance**: **Directly relevant.** Bond order × category
  gives a continuous Δχ estimate. If the McConnell T2 residual
  correlates with bond order perturbation near the mutation site,
  that's a finding: fixed Δχ per category misses environment-
  dependent anisotropy.

---

## Ring Proximity Shells (reporting, not physics)

### [x] Shell boundaries = 3.0, 5.0, 8.0, 12.0 Å
- **Source**: Arbitrary binning for ring proximity counting.
- **Physics**: Not load-bearing. Used for diagnostic grouping only.
- **No sweep needed.**

---

## Concrete Description Strategy

### Filter rejection logging

The KernelEvaluationFilter system already tracks per-filter rejection
counts. Extending it to write **exhaustive structured descriptions**
enables a second analytical path distinct from numerical analysis:

Per evaluation (accepted or rejected):
```
atom_index, atom_element, atom_role, residue_type, residue_index,
source_type (ring/bond), source_index, ring_type (if ring),
bond_category (if bond), distance, rho, z, cone_angle,
source_extent, filter_name, decision (accept/reject),
reject_reason (if rejected)
```

This enables:
- Pattern recognition: "every rejected atom near HIS is backbone N
  within 2 residues" → finding about sequential exclusion
- Audit trail: click an atom in the viewer, see its evaluation log
- Comparison: same atom, same ring, different cutoff shape → what
  changed in the acceptance set?

### MOPAC loading per conformation

MOPAC data belongs on the ProteinConformation (one per pose),
not the Protein (one per sequence). Architecture:

- MOPAC runs offline on each conformation's XYZ file
- Output files stored alongside pose data
- MopacResult (ConformationResult) reads at extraction time
- Loader reads MOPAC files when loading a trajectory pose set
- Each MDFrameConformation gets its own charges and bond orders
- ORCA validates on one frame; MOPAC tracks the trend across all 10
