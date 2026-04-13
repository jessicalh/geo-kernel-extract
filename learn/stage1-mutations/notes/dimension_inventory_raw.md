# Dimension Inventory

2026-04-12, revised 2026-04-13 with 720-protein final numbers
and proper physical-dimension vocabulary.

---

## The dimensions are physics, not tools

The extractor computes 121 T2 kernels.  These are not 121
independent dimensions.  They are multiple ways of accessing a
small number of physical dimensions — the distinct geometric
sources of shielding tensor perturbation at an atom.

The vocabulary: a **dimension** is a physical quantity with its
own angular symmetry.  A **kernel** is one way to compute a
projection onto that dimension.  Multiple kernels can access the
same dimension (BS and HM both access ring current).  The question
is: how many independent dimensions does the T2 shielding
perturbation actually have?

---

## The three primary dimensions

### Dimension 1: Ring current magnetic field

The circulating pi-electrons in an aromatic ring produce a magnetic
dipole field.  T2 angular pattern: (1 - 3cos^2 theta) dipolar.
Distance: r^-3.

Accessed by: Biot-Savart (wire-loop model), Haigh-Mallion (surface
integral model), ring susceptibility (point dipole at ring center).
All three compute the same physics with different approximations.
BS and HM agree to 0.1% in the far field (cosine 0.9999) and 0.2%
near field (cosine 0.9977).

Per-element contribution (fair-set R-squared, 720 proteins):
H = dominant (0.928 total, ring current carries most).
C = minor (0.075 from BS alone).
N, O = moderate (~0.08-0.13 from BS alone).

### Dimension 2: Electric field gradient

Partial charges on all atoms create an electric field gradient
tensor at each atom.  Traceless symmetric (pure T2) by Gauss's
law.  Distance: r^-3 (point charge sum).

Accessed by: Coulomb EFG with three charge sources —
- ff14SB (fixed, invariant to environment)
- MOPAC PM7 (Mulliken charges, respond to electronic environment)
- AIMNet2 (neural network charges, respond to geometry)

These are not three dimensions.  They are three instruments
measuring the same physical quantity with different accuracy.
ff14SB and MOPAC are nearly parallel (cosine 0.99 on raw aromatic
EFG).  AIMNet2 is orthogonal to both (cosine 0.34, near random
in 5D) — it captures a different projection of the charge-
polarization response.

Per-element contribution:
H: EFG_aro alone = 0.711 (ff14SB), 0.757 (MOPAC).
C: EFG_aro alone = 0.215 (ff14SB), 0.388 (MOPAC).
N: EFG_aro alone = 0.068 (ff14SB), 0.066 (MOPAC) — both weak.
O: EFG_aro alone = 0.160 (ff14SB), 0.150 (MOPAC).

### Dimension 3: Bond magnetic anisotropy

Each covalent bond has an anisotropic magnetic susceptibility that
creates a dipolar field at nearby atoms.  The McConnell tensor
M_ab/r^3 has T0+T1+T2 (asymmetric, non-traceless).  Distance:
r^-3.

Accessed by: McConnell (categorical Delta_chi per bond type),
MOPAC McConnell (Wiberg bond-order-weighted, continuous).  The
bond-order weighting accesses the charge-polarization dimension
through bond character rather than through point charges.

Per-element contribution:
C: MC_aromatic_total enters forward selection at rank 4 for
   geo-only.  Bond anisotropy matters for carbon because the
   removed aromatic bonds had the largest Delta_chi.
H, N, O: bond anisotropy contributes < 0.01 individually.

---

## The charge-polarization sub-dimension

The gap between static charges (ff14SB) and responsive charges
(MOPAC, AIMNet2) defines a sub-dimension: how the electronic
structure responds to the local environment.

For the mutation delta, this sub-dimension contributes:
- H: +0.007 (full minus geo-only fair R-squared).  Negligible.
- C: **+0.197**.  The dominant gap.
- N: +0.035.
- O: +0.016.

Carbon's large gap means the charge redistribution from removing
an aromatic ring changes the EFG angular pattern in ways that
geometry alone cannot predict.  This is accessed by MOPAC charges
in the mutation analysis.  For ensemble work (Stage 2), AIMNet2
charges (0.17s vs 10 min) access the same underlying dimension —
the charge response to conformational change — through a different
projection.

The mutation-specific projection (aromatic ring removal) favours
MOPAC.  The conformational projection (backbone flexibility)
may favour AIMNet2.  The dimension is the same; the projection
onto it depends on the perturbation.

---

## Secondary dimensions

### Pi-quadrupole

Electric quadrupole from the pi-electron cloud.  r^-5, pure T2.
Shares geometry with Dimension 2 (charge-based) but has different
angular symmetry (quadrupolar vs dipolar).  Contributes modestly
for N (PQ_total enters forward selection at rank 2 in geo-only).

### Dispersion

Van der Waals from ring vertices.  r^-6/r^-8.  Enters forward
selection at rank 2-5 for C, N, O in geo-only basis.  Dispersion
group R-squared for O = 0.272 (larger than ring current 0.231 or
EFG 0.165 in geo-only).

Dispersion was originally listed as a minor effect.  It is not.
The dispersion kernel captures short-range angular structure
through DispChi (dispersion scalar x ring susceptibility T2),
which accesses the ring-proximity dimension at close range where
the r^-6 scalar is large.

### Ring type specificity

Which aromatic ring type was removed (PHE, TYR, TRP, HIE).
Accessed as a categorical variable (mutation type identity).
Adds +0.027 for H, +0.121 for C, +0.220 for N, +0.195 for O.

This is Case (1995)'s ring-type-specific intensity factors
extended to T2.  It is a real physical variable: different ring
types have different current intensities, charge distributions,
and susceptibility anisotropies.

---

## The normalization dimension

Per-protein normalization strips magnitude (dominated by protein
size and ring count) to isolate angular structure.  The stripped
magnitude is real signal — kernel scale factors add +0.06 R-squared
uniformly.

The effect on the eigenspectrum:
- Raw: top 10 PCs = 99.3% variance.  One dominant source.
- Normalised: top 10 PCs = 24%.  Angular structure in the tail.

The prediction lives in the fine angular structure, not the
dominant modes.  For H, the top 10 variance PCs carry 0% of the
ridge weight.  The ridge sees the angular disagreement between
kernel families — the thing that distinguishes ring current from
EFG at each atom.  This is the angular fine structure that IS
the physics.

---

## Proven zeros (dimensions tested and rejected)

These are candidate dimensions that turned out to be projections
of the three primary dimensions — the kernels already encode them.

| Candidate | R-squared | Why zero |
|-----------|----------|---------|
| MOPAC valence | 0.000 | Sum of bond orders = bond count. Kernels know geometry. |
| Molecular dipole | 0.000 | Global property. Per-atom kernels don't need it. |
| Bond order scalar | 0.000 | Scalar projection of Dimension 3. T2 carries more. |
| H-bond shielding | 0.002 | Unchanged in mechanical mutants (backbone fixed). |
| Solvation EFG (DeltaAPBS) | 0.005 | 2% of direct aromatic EFG. Local, through-space effect. |
| SASA | < 0.001 | Solvent exposure unchanged by sidechain-only mutation. |
| DSSP SS8 | 0.001 | Backbone structure unchanged in mechanical mutants. |
| DSSP H-bond energy | < 0.001 | Same as H-bond — backbone fixed. |
| AIMNet2 charge (scalar) | < 0.001 | Scalar loses angular info. EFG tensor is the right form. |
| MOPAC charge (scalar) | < 0.001 | Same. |

Many of these are correctly zero BECAUSE of the mechanical mutant
design.  SASA, DSSP, H-bond don't change when only the sidechain
changes.  They are not physically irrelevant — they are irrelevant
to THIS perturbation.  For conformational variation (Stage 2) they
may matter.

---

## Dia/para decomposition (from DFT)

The DFT decomposes shielding into diamagnetic (sigma_d) and
paramagnetic (sigma_p) contributions.  For the mutation delta:

| Element | |delta_dia| | |delta_para| | |delta_total| |
|---------|-----------|-------------|--------------|
| H | 6.44 ppm | 6.49 ppm | 0.78 ppm |
| C | 7.42 | 7.67 | 1.54 |
| N | 6.96 | 7.27 | 1.58 |
| O | 6.55 | 7.15 | 2.13 |

The dia and para perturbations are individually ~7 ppm but nearly
cancel.  The geometric kernels predict the net (R-squared 0.70 for
H) but not the individual channels (R-squared < 0.05).

The para/dia variance ratio increases H→O (1.02→1.20), confirming
Saito et al. 2010 in the T2 delta channel: paramagnetic variability
grows for heavier atoms.

The cancellation means the mutation delta is a small angular
residual on top of a large self-cancelling perturbation.  This is
why per-protein normalization is necessary and why the per-protein
ceiling exists.

---

## Distance-dependent dimensionality (110-protein numbers, re-run on 720 pending)

Raw kernel space has 3 predictive dimensions at every distance
band.  R-squared varies:

| Element | 0-4 A | 4-8 A | 8-12 A | 12+ A |
|---------|-------|-------|--------|-------|
| H | 0.919 | 0.959 | 0.969 | 0.908 |
| C | 0.501 | 0.556 | 0.951 | 0.853 |
| N | 0.661 | 0.224 | 0.906 | 0.818 |
| O | NaN | 0.195 | 0.745 | 0.748 |

The same 3 dimensions work everywhere for hydrogen.  Carbon gets
better at distance (near-field hard, far-field easy).  Nitrogen
near-field has 8 predictive dimensions — the multi-mechanism
complexity concentrates close to the ring.

---

## Final numbers (720 proteins, 446,006 atoms)

| Element | n | Fair R-squared | Geo-only | Gap |
|---------|--------|--------------|---------|------|
| H | 230,135 | 0.928 | 0.921 | +0.007 |
| C | 133,488 | 0.562 | 0.365 | +0.197 |
| N | 39,954 | 0.380 | 0.346 | +0.035 |
| O | 42,429 | 0.382 | 0.366 | +0.016 |
| Weighted | 446,006 | **0.718** | 0.650 | |

Down from 0.818 on 110 proteins.  The physics is the same.  The
cross-protein generalisation is harder with more structural
diversity.  H held well (0.949→0.928).  Heavy atoms dropped:
the per-protein ceiling (0.81 within) hasn't changed — the gap
is cross-protein variation that local geometric kernels cannot
and should not capture.
