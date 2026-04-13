# Physics Group Grounding

2026-04-13.  For each kernel group: what physical interaction, what
equation, what determines the T2 direction, and what the 400-protein
data shows.  This is the viva-level account.

All cosine numbers are mean |cos| across matched atoms.  Random in
5D = 0.45.  All R² are ridge on raw (unnormalised) kernels.

---

## 1. Ring current (17 kernels, R² = 0.88 H / 0.13 C / 0.12 N / 0.20 O)

### The physics

Circulating pi-electrons in an aromatic ring produce a magnetic
dipole field.  The shielding tensor at a probe atom is:

    G_ab = -n_b * B_a * PPM_FACTOR

where n is the ring normal and B is the B-field from the current
loop (Biot-Savart) or surface integral (Haigh-Mallion).  This is
a rank-1 outer product: n ⊗ B.  It has T0 (the isotropic ring
current shift), T1 (from the asymmetric n/B coupling), and T2.

Pople 1956: concept.  Johnson & Bovey 1958: JB wire model.
Haigh & Mallion 1979: surface integral.  Case 1995: calibration.

### What determines the T2 direction

Two vectors: the ring normal n and the B-field direction at the
atom.  The B-field depends on the atom's position relative to the
ring geometry (center, vertices, lobe offset d).  In the far field,
B is purely along the ring normal (magnetic dipole limit) and the
T2 has the (1 - 3cos²θ) angular pattern.  In the near field, B
has transverse components that depend on the specific ring vertex
geometry.

### What the data shows

- BS and HM within-group cosine = 0.534.  They are the SAME physics
  with different numerics — the within-group spread comes from
  per-ring-type differences (PHE ring is different geometry from TRP).
- Ring current → target cosine: 0.70 for H, 0.60-0.64 for heavy
  atoms.  The T2 direction from ring current aligns with the DFT
  target more for H than for heavy atoms.
- Group R² = 0.88 for H (dominant), 0.13 for C (minor).
  Consistent with Saito 2010: proton shielding is diamagnetic-
  dominated, and ring current is the dominant through-space
  diamagnetic perturbation.
- BS-PHE and HM-PHE have cosine 0.70 with target, R² = 0.24.
  A single ring type's ring current explains 24% of H variance.

---

## 2. Coulomb EFG — ff14SB charges (3 kernels, R² = 0.72 H / 0.22 C / 0.07 N / 0.17 O)

### The physics

Point charges on all atoms create an electric field gradient at
each atom:

    V_ab = sum_j q_j * (3 d_a d_b / r^5 - delta_ab / r^3)

This is the dipolar kernel K_ab weighted by charge q_j and summed
over all source atoms.  The result is symmetric and traceless by
Gauss's law (div E = 0 outside charges).  Pure T2.

Buckingham 1960: sigma_T2 = gamma * V_T2.

### What determines the T2 direction

The angular distribution of charges around the atom.  If the
charge environment is asymmetric (e.g., near a carbonyl C=O with
q_C = +0.6, q_O = -0.6), the EFG tensor points along the charge
asymmetry axis.  The T2 direction encodes WHERE the dominant
charges are, angularly.

The backbone/aromatic decomposition: EFG_bb from backbone charges
(unchanged in mutation), EFG_aro from aromatic charges (removed
in mutation).  EFG_aro dominates because the aromatic charges ARE
the perturbation.

### What the data shows

- EFG → target cosine = 0.93 for H, 0.77 for O.  Tightest
  alignment of any group.  The charge distribution asymmetry
  around each atom points in almost exactly the direction of the
  DFT target T2.
- ff14SB ↔ MOPAC EFG cross-group cosine = 0.58-0.67.  Same
  equation, different charges.  The angular patterns are
  correlated but not identical — charge polarization rotates
  the EFG.
- EFG_aro alone: R² = 0.72 for H.  A SINGLE kernel (5 components)
  explains 72% of H T2 variance.  This is remarkable and means
  the charge distribution around the removed ring encodes most
  of the angular structure.
- For N: R² = 0.07.  EFG is weak.  Nitrogen's lone pair and
  paramagnetic-term dominance mean the electrostatic perturbation
  has a different effect on the tensor than for H.

---

## 3. Coulomb EFG — MOPAC charges (3 kernels, R² = 0.76 H / 0.39 C / 0.07 N / 0.16 O)

### The physics

Same equation as ff14SB EFG.  Different charges: PM7 Mulliken
charges from a semi-empirical QM calculation.  These respond to
the local electronic environment — a backbone nitrogen near a
charged sidechain has a different MOPAC charge than the same
nitrogen far from charged groups.

### What determines the T2 direction

Same as ff14SB, but the charge MAGNITUDES differ.  MOPAC charges
capture electronic redistribution from the aromatic ring.  This
rotates the EFG tensor slightly (ff14SB↔MOPAC cosine = 0.58-0.67)
and changes the magnitude.

### What the data shows

- MOPAC EFG R² = 0.39 for C vs ff14SB = 0.22.  75% improvement.
  The charge polarization response to the aromatic ring removal
  is specifically important for carbon.
- For H: MOPAC = 0.76 vs ff14SB = 0.72.  Modest improvement.
  H is already well-served by the angular pattern; magnitude
  refinement helps less.
- For N: 0.07 vs 0.07.  Neither charge source helps.  The EFG
  dimension is weak for nitrogen regardless of charge quality.
- MOPAC is the charge-polarization dimension.  The gap
  (MOPAC - ff14SB) = +0.17 for C, +0.04 for H, ≈0 for N/O.

---

## 4. Bond magnetic anisotropy — McConnell (7 kernels, R² = 0.12 H / 0.02 C / 0.03 N / 0.04 O)

### The physics

Each covalent bond has an anisotropic magnetic susceptibility
(Δχ: different parallel vs perpendicular to the bond axis).  This
creates a dipolar shielding field at nearby atoms.

    sigma_ab = (Δχ/3) * M_ab / r^3

where M is the full McConnell tensor:

    M_ab = 9 cos_θ d_a b_b - 3 b_a b_b - (3 d_a d_b - delta_ab)

Asymmetric (d_a b_b ≠ b_a d_b).  Non-traceless (Tr = 3cos²θ - 1).
T0 + T1 + T2 all non-zero.

McConnell 1957.  Full derivation in GEOMETRIC_KERNEL_CATALOGUE.md.

### What determines the T2 direction

Two vectors: the bond axis b and the atom-to-midpoint direction d.
Each bond creates a different T2 pattern because each bond points
in a different direction.  The total McConnell T2 is a SUM over
all bonds — the direction depends on which bonds dominate
(determined by Δχ and distance).

Per-category decomposition: backbone (unchanged in mutation),
aromatic (removed in mutation), CO/CN nearest (proximity effects).
MC_aromatic_total carries the mutation signal.

### What the data shows

- Bond aniso → target cosine = 0.40.  Near random (0.45).  The
  bond anisotropy angular pattern is genuinely independent of the
  DFT target direction for most atoms.
- Within-group cosine = 0.43.  The different bond categories (CO,
  CN, backbone, aromatic) produce moderately correlated T2 patterns
  — they share the 1/r^3 dipolar structure but point along
  different bond axes.
- Group R² is low everywhere (0.02-0.12).  In mechanical mutants,
  most bonds don't change.  The signal comes from MC_aromatic_total
  (the removed aromatic bonds) — enters forward selection at
  rank 3 for C (delta +0.015).
- Bond aniso ↔ ring current cosine = 0.40.  Independent.  Magnetic
  dipole from a bond midpoint vs magnetic dipole from a ring center
  — different source positions, different T2 directions.

---

## 5. Bond anisotropy — MOPAC weighted (6 kernels, R² = 0.11 H / 0.02 C / 0.03 N / 0.04 O)

### The physics

Same McConnell equation, but each bond contribution weighted by
its Wiberg bond order — a continuous quantum measure of electron
sharing.  A C=O with bond order 1.85 contributes 1.85× what a
C-O with order 1.0 contributes.

### What determines the T2 direction

Same geometry as unweighted McConnell, but the RELATIVE
contributions of bonds change.  Strong double bonds dominate over
weak single bonds.  This rotates the total T2 vector.

### What the data shows

- MC ↔ MopacMC cosine = 0.44-0.47.  Modestly correlated.  Bond-
  order weighting changes which bonds dominate, producing a
  measurably different angular pattern.
- For N: MopacMC_total enters forward selection at rank 4
  (delta +0.017).  The peptide bond's partial double-bond character
  (Wiberg 1.3-1.5) means bond-order weighting captures electronic
  coupling through the peptide bond that categorical Δχ misses.
  This is why MopacMC appears for N but not for H — nitrogen
  feels the peptide bond's electronic structure directly.

---

## 6. Pi-quadrupole (9 kernels, R² = 0.02 H / 0.01 C / 0.07 N / 0.03 O)

### The physics

The pi-electron cloud above and below an aromatic ring creates
a quadrupole electric field.  Derived from Stone's T-tensor
formalism:

    G_ab = 105 dn² d_a d_b / r^9 - 30 dn (n_a d_b + n_b d_a) / r^7
           - 15 d_a d_b / r^7 + 6 n_a n_b / r^5
           + delta_ab (3/r^5 - 15 dn²/r^7)

where dn = d · n (height above ring plane).  Symmetric, traceless.
Pure T2.  Leading decay: r^-5.

Stone 2013, Ch. 3.  Buckingham 1959 for the scalar E-field version.

### What determines the T2 direction

The ring normal n AND the height above the ring plane.  The angular
pattern has terms involving n_a n_b that the dipolar kernel doesn't
have — quadrupolar angular dependence is genuinely different from
dipolar.

### What the data shows

- PQ ↔ ring current cosine = 0.45.  Random.  Dipole (r^-3) and
  quadrupole (r^-5) produce orthogonal angular patterns.  This is
  the direct empirical confirmation of Stone 2013's mathematics:
  different multipolar orders are angularly independent.
- PQ ↔ EFG cosine = 0.42.  Also random.  Despite both being
  "electric" effects, the quadrupole from the ring and the EFG
  from point charges have different angular symmetry.
- For N: PQ_total enters forward selection at rank 2 (delta +0.043).
  This is the largest single delta after EFG_aro for nitrogen.
  The r^-5 decay means PQ is only significant for atoms close to
  rings — and nitrogen atoms near aromatic rings are in the
  regime where multiple multipolar orders (dipole + quadrupole)
  contribute comparable signal.
- Measured distance slope: -5.05 (theory -5).  Three-significant-
  figure agreement with the quadrupolar prediction.

---

## 7. Dispersion (8 kernels, R² = 0.40 H / 0.10 C / 0.02 N / 0.06 O)

### The physics

London van der Waals attraction between the atom and each ring
vertex:

    K_disp_ab = sum_vertices C6 * (3 d_a d_b / r^8 - delta_ab / r^6)

The tensor part goes as r^-8 (anisotropic).  The scalar part
as r^-6 (isotropic).  A CHARMM switching function tapers to zero
at 5Å.

London 1937.  Brooks et al. 1983 for the switching function.

The DispChi kernel is: disp_scalar × ring susceptibility T2.
This is a PRODUCT — the r^-6 scalar measures proximity and the
ring susceptibility T2 provides the angular pattern.  DispChi is
not independent angular information; it's ring susceptibility
GATED by dispersion proximity.

### What determines the T2 direction

For the raw dispersion tensor: the direction from atom to ring
vertices, weighted by C6/r^8.  Only significant at short range
(< 5Å) due to the switching function.

For DispChi: the ring susceptibility direction (ring normal
dependent), activated only when the atom is close enough for
the dispersion scalar to be non-zero.

### What the data shows

- Dispersion → target cosine = 0.83-0.85 for H.  HIGHER than
  any ring current kernel (0.70).  This is surprising: dispersion
  points in the right direction more than ring current does.
- But per-kernel R² = 0.13.  Dispersion only activates for atoms
  within 5Å of a ring vertex.  It's angularly correct but
  spatially sparse.
- Dispersion ↔ EFG cosine = 0.53.  The highest cross-group cosine
  in the matrix (tied with solvation↔EFG).  Dispersion and EFG
  share geometric information because the C6 coefficients are
  partially charge-dependent.
- Dispersion ↔ ring current cosine = 0.50.  Moderate.  DispChi
  uses the ring susceptibility T2 direction, which is ring-normal-
  dependent like ring current.
- Group R² = 0.40 for H.  Third after ring current (0.88) and
  MOPAC EFG (0.76).  Dispersion is NOT a minor effect for H — it
  captures 40% of variance with 8 kernels.
- For N: group R² = 0.015.  Dispersion is irrelevant.  The
  switching function limits it to atoms within 5Å, and at that
  range N is in the multi-mechanism regime where dispersion is
  one of many weak contributors.
- Forward selection: Disp enters rank 2-5 for C and O.  It acts
  as a proximity-gated angular indicator — when the atom is close
  enough to the ring for dispersion to activate, the DispChi
  angular pattern is informative.

---

## 8. Solvation (2 kernels, R² = 0.01 everywhere)

### The physics

APBS solves the Poisson-Boltzmann equation for the solvated E-field.
DeltaAPBS = APBS_EFG - Coulomb_EFG: the difference between the
solvated and vacuum electric field gradient.  Measures what the
solvent dielectric does to the EFG.

### What determines the T2 direction

The shape of the protein surface relative to the atom.  The solvent
dielectric screens the vacuum EFG, and the screening depends on
how much solvent is between the atom and each charge.

### What the data shows

- R² ≈ 0.01 everywhere.  2% of direct aromatic EFG.  The solvation
  perturbation to the EFG is negligible for mechanical mutants
  because the solvent boundary barely changes when a sidechain is
  replaced with alanine.
- Solvation ↔ EFG cosine = 0.53-0.55.  The highest cross-group
  cosine.  This is expected: DeltaAPBS is a CORRECTION to the
  Coulomb EFG, so it's angularly correlated with the thing it
  corrects.
- Correctly identified as a negative control: unchanged backbone
  means unchanged solvation.

---

## Cross-group independence structure

The group-to-group cosine matrix has a clear structure:

**Near-random (cos ≈ 0.40, independent):**
- Bond aniso ↔ everything else.  Bond axes point in locally
  random directions relative to ring geometry.
- MOPAC bond ↔ everything except unweighted bond.

**Moderately correlated (cos ≈ 0.45-0.55):**
- Ring current ↔ dispersion (0.50).  Both depend on ring proximity
  and ring normal.
- EFG ↔ dispersion (0.53).  Charge-C6 correlation.
- EFG ↔ solvation (0.53).  Solvation is a correction to EFG.
- Ring current ↔ quadrupole (0.45).  Different multipolar orders
  but same ring geometry.

**Correlated (cos ≈ 0.58-0.67):**
- ff14SB EFG ↔ MOPAC EFG (0.58-0.67).  Same equation, different
  charges.  Higher for N (0.67) where both charge sources produce
  similar EFGs because there are fewer aromatic charges near N
  atoms.

**Within-group:**
- Ring current internal (0.53).  Per-ring-type differences within
  the same physics.
- Quadrupole internal (0.56).  Per-ring-type differences in the
  quadrupole pattern.
- Dispersion internal (0.47-0.56).  Per-ring-type differences.
  Higher for N (0.56) — dispersion kernels for different ring
  types agree more for N because fewer N atoms are close enough
  for dispersion to activate.

This independence structure is ELEMENT-INDEPENDENT.  The geometry
of the kernel space doesn't change from H to O.  What changes is
how much each group contributes to the TARGET — and that's
determined by the element's electronic structure.
