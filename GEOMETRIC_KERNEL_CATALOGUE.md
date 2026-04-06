# Geometric Kernel Catalogue

What each classical calculator computes, why, and what angular structure
(T2) it produces. Parameter values are not specified here — they are
learnable weights on these geometric kernels. The kernels ARE the physics.

This document is the mathematical foundation for implementing the
eight classical calculators. It describes the geometry-to-tensor
mapping of each physical effect, independent of parameterisation.

---

## Unit Convention (enforced, from old project)

All geometric calculations in SI at the computation boundary, then
stored in mixed units matching the old project:

- Positions: Angstroms (input and storage)
- Current: nanoamperes (ring current intensity parameter)
- B-field from Biot-Savart: Tesla (SI, computed internally)
- Geometric kernel G: dimensionless (after PPM_FACTOR multiplication)
- Dipolar kernel K: Angstrom^-3 (1/r^3 with r in Angstroms)
- EFG: V/Angstrom^2 (from charge in e, distance in Angstroms)
- Shielding contribution: ppm (kernel × parameter)

Unit conversion for Biot-Savart:
```
mu_0/(4*pi) = 1e-7 T·m/A (exact in SI)
1 Angstrom = 1e-10 m
1 nanoampere = 1e-9 A
```

Positions converted to metres, current to amperes, B computed in
Tesla. Then: G_ab = n_b × B_a × PPM_FACTOR (PPM_FACTOR = 1e6).
G is dimensionless.

For dipolar kernels (McConnell, ring susceptibility, H-bond):
K computed directly in Angstrom^-3 (no SI conversion needed —
positions already in Angstroms).

---

## The Shared Dipolar Kernel

Four of eight calculators share one geometric core:

```
K_ab(d) = 3 d_a d_b / r^5 - delta_ab / r^3
```

where d = r_atom - r_source, r = |d|.

### Properties of K

- **Symmetric**: K_ab = K_ba
- **Traceless**: Tr(K) = 3|d|^2/r^5 - 3/r^3 = 0
- **Units**: Angstrom^-3 (when d is in Angstroms)
- **Rank**: full rank-2 tensor. Since Tr = 0 and K is symmetric,
  the SphericalTensor decomposition gives T0 = 0, T1 = 0,
  all information lives in T2 (5 components).
- **Angular pattern**: The (3cos^2 theta - 1) dependence on the
  angle between d and any reference direction.

### K is the gradient of the gradient of 1/r

```
K_ab = -del_a del_b (1/r) = 3 d_a d_b / r^5 - delta_ab / r^3
```

This is the same tensor that appears in:
- Magnetic dipole fields (McConnell, ring susceptibility)
- Electric field gradients (Coulomb EFG)
- Any point-source dipolar interaction

### The McConnell scalar is NOT the trace

The McConnell scalar factor:
```
f(d, b_hat) = (3 cos^2 theta - 1) / r^3
            = sum_ab K_ab * b_a * b_b
```

where theta is the angle between d and the bond axis b_hat.
This is the double contraction of K with the bond direction,
not the trace of K (which is zero).

The isotropic shielding from McConnell (what NMR measures as
a chemical shift perturbation) is:
```
sigma_iso = (Delta_chi / 3) * f = (Delta_chi / 3) * (3cos^2 theta - 1) / r^3
```

This is non-zero despite K being traceless, because it comes from
contracting K with the anisotropic susceptibility direction, not
from taking the trace.

### What K's T2 encodes physically

K_ab at an atom tells you: "if there is a magnetic dipole at the
source point, what angular pattern of shielding does it create here?"
The T2 components encode the orientation of this dipolar field pattern
relative to the lab frame.

Near a C=O bond midpoint, K's T2 shows the angular distribution
of bond anisotropy shielding. Near a ring center, K's T2 shows the
angular distribution of ring susceptibility shielding.

### Four calculators using K

| Calculator | Source point | Reference direction | Physical constant |
|------------|-------------|--------------------|--------------------|
| McConnell | bond midpoint | bond axis (for scalar f) | Delta_chi per bond category |
| RingSusceptibility | ring center | ring normal (for scalar f) | Delta_chi_ring per ring type |
| HBond | H-bond partner | D-H-A direction | eta per element |
| Coulomb EFG | each charged atom | (none — sum over all) | charges (from ff14SB/xTB) |

---

## Calculator 1 & 2: Ring Current (Biot-Savart and Haigh-Mallion)

These two calculators model the SAME physical effect (pi-electron
circulation in aromatic rings) with different mathematical
approximations.

### Biot-Savart: Johnson-Bovey Double-Loop Model

**Physical picture**: The pi-electrons above and below an aromatic ring
are modelled as two current loops at height ±d from the ring plane.
Each loop carries half the total current I/2.

**Geometric kernel** (per ring, at atom position r):

Step 1: Compute B-field from each wire segment (Biot-Savart law):
```
B_segment = (mu_0 / 4pi) * I * (dl x dA) / |dl x dA|^2
            * (dl . dA/|dA| - dl . dB/|dB|)
```
where:
- dl = b - a (wire segment vector)
- dA = r - a, dB = r - b (displacement from endpoints to field point)
- All in SI (metres, amperes → Tesla)

Step 2: Sum over all segments of both loops:
```
B_total = sum over segments of upper loop (z = +d, current I/2)
        + sum over segments of lower loop (z = -d, current I/2)
```

Step 3: Construct geometric kernel:
```
G_ab = -n_b * B_a * PPM_FACTOR
```
where n is the ring normal. The minus sign comes from the shielding
tensor definition: sigma_ab = -dB_a^sec / dB_{0,b}.

**Units**: B in Tesla, G dimensionless.

**Sign convention**: sigma_ab = I_type * G_ab, where I_type is the
ring current intensity (negative for diamagnetic ring currents).
With the minus sign in G, this gives sigma > 0 above the ring
(shielded) with literature I values. Verified analytically:
I=-12, atom 3A above PHE -> sigma = +1.40 ppm (Case 1995).

**T2 structure of Biot-Savart**:

G_ab = n_b * B_a is an outer product of two vectors. This is
rank-1 by construction (the matrix has rank 1, meaning only one
non-zero eigenvalue if B is parallel to n, or at most two if B has
components transverse to n).

Decomposition:
- T0 = Tr(G)/3 = n . B * PPM_FACTOR / 3. Non-zero when B has a
  component along the ring normal. This IS the isotropic ring current
  shift.
- T1 = antisymmetric part of G = (n x B) components / 2 * PPM_FACTOR.
  Non-zero because G is NOT symmetric (G_ab = n_b B_a ≠ n_a B_b
  in general).
- T2 = traceless symmetric part. Contains angular information about
  how the B-field projects onto the ring normal.

The T2 of Biot-Savart depends on:
- The direction of B relative to n (changes with position around ring)
- The lobe offset d (changes the B-field shape: larger d concentrates
  the field into sharper lobes above/below the ring)
- The ring geometry (center, vertices, radius)

The T2 does NOT depend on the intensity I (which scales magnitude
uniformly). So T2 constrains the GEOMETRIC parameters (d, ring
geometry) while T0 constrains the MAGNITUDE parameter (I).

### Haigh-Mallion: Surface Integral

**Physical picture**: The pi-electron cloud is modelled as a uniformly
magnetized disk. The shielding tensor comes from integrating the
dipolar kernel over the ring surface.

**Geometric kernel** (per ring, at atom position r):
```
H_ab(r) = integral_S [ 3 rho_a rho_b / rho^5 - delta_ab / rho^3 ] dS
```
where rho = r - r_s (field point minus surface point), and the
integral is over the ring polygon area.

Numerical implementation: fan triangulation from ring centroid,
7-point Gaussian quadrature per triangle (Stroud T2:5-1). Adaptive
subdivision when atom is within 2.0 A (level 1) or 1.0 A (level 2)
of the ring.

**Tensor structure of the raw integral H_ab**:

The integrand is the dipolar kernel K_ab, which is symmetric and
traceless. The integral preserves both properties (trace is linear):
Tr(H) = integral Tr(K) dS = integral 0 dS = 0. H_ab is symmetric
and traceless: pure T2, T0 = 0, T1 = 0. Verified at machine
precision on 1UBQ (max asymmetry 3.9e-17, max |Tr| 1.2e-14) and
A0A7C5FAR6 (max |Tr| 7.1e-15). Units: Angstrom^-1.

**Full shielding tensor construction** (Boyd & Skrynnikov JACS 2002):

The raw integral H is NOT the shielding tensor. The induced ring
current depends on the flux through the ring (n . B_0), so the
shielding tensor sigma_ab = -dB_a^sec / dB_{0,b} must carry n_b.
The magnetic field from the magnetised surface is V_a = sum_c H_ac n_c
(i.e. V = H . n). The full shielding kernel is:

```
G_ab = n_b * V_a = n_b * (H . n)_a
```

This is rank-1 (same structure as BS), with:
- T0 = n . V / 3 = n^T H n / 3 (the HM isotropic shift, non-zero
  despite H being traceless — same mechanism as McConnell scalar)
- T1 != 0 (from antisymmetric n_b V_a coupling)
- T2 != 0 (angular anisotropy)

Both H (pure T2 geometric kernel) and G (full shielding kernel with
all irreps) are stored. The model can use either or both.

**T2 structure comparison (BS vs HM)**:

BS and HM model the same physical effect with different mathematical
approximations (line integral vs surface integral). Both produce
rank-1 shielding tensors G = n (x) V, but compute V differently:
- BS: V = B (wire segment B-field in Tesla, then * PPM_FACTOR)
- HM: V = H . n (surface integral contracted with normal, in A^-1)

Whether these agree or disagree on T2, and how that depends on
distance and geometry, is an empirical question. Far-field T0 sign
agreement verified at 99.3% (878/884 pairs beyond 8A on 1UBQ).

### BS-HM relationship

Both BS and HM have intensity parameters per ring type.
Traditionally set equal (same I for both). The model can learn
different intensities for BS vs HM. Whether this reveals
distance-dependent accuracy differences between the two
approximations is something the data will show.

### Fused ring representation (TRP)

TRP produces three ring objects: TRP5 (pyrrole, I=-6.72), TRP6
(benzene, I=-12.48), TRP9 (indole perimeter, I=-19.2). Note
I(TRP9) = I(TRP5) + I(TRP6) exactly. The model receives all
three as independent per-type features.

Measured on 720 proteins (87,234 atom-TRP triplets):
- BS: T0(TRP5) + T0(TRP6) / T0(TRP9) = 1.000. Wire-segment
  superposition is exact — the JB model is additive for
  non-overlapping current loops.
- HM: same ratio = 1.000. Both models are perfectly additive
  after RingBondedExclusionFilter excludes ring vertices and their
  bonded neighbours from through-space evaluation.

**Correction (2026-04-02):** The previously reported HM ratio of 1.127
(13% surface integral excess) was an artifact. Atoms on the shared
CD2-CE2 edge were being evaluated by both TRP5 and TRP6 surface
integrals while sitting inside the source distribution. The
RingBondedExclusionFilter (topology check via bond graph) correctly
excludes these atoms, and the excess disappears.

The model sees TRP5+TRP6 as the component decomposition and TRP9
as the whole-system representation.

### Sign convention

G_ab = -n_b * B_a * PPM_FACTOR (BS) and G_ab = -n_b * V_a (HM).
The minus sign comes from the shielding tensor definition:
sigma_ab = -dB_a^sec / dB_{0,b}.

With this convention, sigma = I * G gives the correct physical sign
using literature ring current intensities (I < 0 for diamagnetic).
Verified analytically: I=-12, atom 3A above PHE centre gives
sigma = +1.40 ppm (shielded), matching Case (1995).
In-plane at 5A gives sigma = -0.16 ppm (deshielded).

---

## Calculator 3: McConnell Bond Anisotropy

**Physical picture**: Each covalent bond has a magnetic susceptibility
that is anisotropic — different parallel and perpendicular to the
bond axis. This creates a dipolar shielding field at nearby nuclei.

**The old formula is incomplete.** The commonly implemented version
stores only the symmetric traceless dipolar kernel K_ab. The full
shielding tensor from first principles (derived below in "RESOLVED:
McConnell Full Shielding Tensor Derivation") is:

```
sigma_ab = (Delta_chi / 3) * M_ab / r^3
```

where M_ab is the full McConnell tensor:

```
M_ab = 9 cos_theta * d_hat_a * b_hat_b
     - 3 b_hat_a * b_hat_b
     - (3 d_hat_a * d_hat_b - delta_ab)
```

with d = r_atom - r_midpoint, r = |d|, d_hat = d/r,
cos_theta = d_hat . b_hat, b_hat = bond direction unit vector.

**Three terms, three roles**:
- 9 cos_theta d_hat_a b_hat_b: asymmetric coupling → T1
- -3 b_hat_a b_hat_b: bond-axis projection → contributes to T0
- -(3 d_hat_a d_hat_b - delta_ab): the dipolar kernel K → T2

**Tensor structure** (per bond):
- T0 = (Delta_chi / 3)(3cos^2 theta - 1) / r^3 — the McConnell
  isotropic shift. Non-zero. Falls directly from the trace of M.
- T1 != 0 — from the asymmetric d_hat/b_hat coupling. This is
  physically real: magnetic interactions produce asymmetric shielding.
- T2 != 0 — angular anisotropy. Richer than just K because M
  includes b_hat-dependent terms.

**What to store** (per bond per atom):
- The full tensor M_ab / r^3 (asymmetric, non-traceless, 9 components)
- SphericalTensor decomposition of M/r^3 (T0 + T1 + T2)
- The symmetric traceless kernel K_ab (useful as a pure geometric
  feature, independent of bond direction)
- The McConnell scalar f = (3cos^2 theta - 1) / r^3 = Tr(M) / r^3
  (redundant with T0 of M, but explicit for feature clarity)

**Why McConnell should dominate T2** (physical reasoning):
1. Bonds are everywhere in a protein (hundreds, vs a few rings)
2. The 1/r^3 decay is steep, so nearby bonds contribute strongly
3. Each bond has a DIFFERENT direction, creating complex angular patterns
4. The anisotropy Delta_chi of C=O bonds is large

**Per-category decomposition**: Bond categories (PeptideCO, PeptideCN,
BackboneOther, SidechainCO, Aromatic, SidechainOther) have different
Delta_chi values. The features accumulate per category. This
decomposition tells the model WHICH BONDS contribute to shielding.

**The midpoint_shift parameter**: The effective dipole center need not
be at the geometric midpoint. For C=O bonds, the electron density is
asymmetric (shifted toward O). A midpoint_shift != 1.0 moves the
effective source along the bond axis, changing d and therefore changing
ALL of T0, T1, T2 simultaneously. The coupled constraint (Gemini's
point 3) means midpoint_shift must improve the fit in all three irreps
or the model rejects it as unphysical.

**Unit prefactor**: The full expression is:

```
sigma_ab = -(mu_0 / 4pi) * (Delta_chi / 3) * M_ab / r^3
```

In our unit system (distances in Angstroms, Delta_chi in cm^3/mol),
the prefactor absorbs unit conversions. The exact numerical prefactor
depends on how Delta_chi is expressed. For implementation: compute
M_ab / r^3 as the pure geometric kernel (Angstrom^-3), multiply by
(Delta_chi / 3) and the appropriate unit conversion to get ppm.
Pin down the exact conversion factor during implementation by
matching against the old code's output.

---

## Calculator 4: Coulomb / Buckingham EFG

**Physical picture**: Partial charges on all atoms create an electric
field E and electric field gradient (EFG) tensor V at each atom. The
Buckingham effect converts these to shielding perturbations.

**Geometric kernel** (at atom i, from all other atoms j):

Electric field (rank-1, Vec3):
```
E_a(i) = sum_j q_j * (r_j - r_i)_a / |r_j - r_i|^3
```

Electric field gradient (rank-2, Mat3, TRACELESS by Gauss's law):
```
V_ab(i) = sum_j q_j * [3 (r_j-r_i)_a (r_j-r_i)_b / |r_j-r_i|^5
                       - delta_ab / |r_j-r_i|^3]
```

Note: V_ab is the SAME dipolar kernel K_ab, weighted by charge q_j
and summed over all sources. Each individual term is traceless; the
sum is traceless. V is ALWAYS traceless for external charges (Laplace
equation: del^2 V = 0 outside charge sources).

**Shielding contribution** (Buckingham 1960):
```
sigma_iso = A_elem * E_z + B_elem * E_z^2      (T0, ppm)
sigma_T2 = gamma_elem * V_T2                    (T2, ppm)
```
where E_z is the E-field projected along the atom's primary bond axis.

**T2 structure**: V_ab is traceless symmetric — pure T2 (like K).
The T2 of the Coulomb EFG encodes the angular asymmetry of the
charge distribution around the atom. gamma converts this directly
to shielding T2.

**Relationship to McConnell**: McConnell T2 comes from bond magnetic
anisotropy. Coulomb T2 comes from the charge environment. These are
different physical sources that produce tensors of the same form
(dipolar kernel) at the same atoms. Whether they reinforce or oppose
depends on the specific geometry and charges. At a backbone amide H
near a C=O bond, both effects are present — the McConnell dipolar
tensor from the bond's magnetic anisotropy and the EFG from the C=O
partial charges. The actual relationship (reinforcing, opposing, or
orthogonal) at each atom is what the data will show.

**Decomposition**: E and V are decomposed by source:
- Backbone atoms: E_backbone, V_backbone
- Sidechain atoms: E_sidechain, V_sidechain
- Aromatic atoms: E_aromatic, V_aromatic
- Solvent: E_solvent = APBS_E - vacuum_E (what solvation adds)

**APBS relationship**: APBS solves the full Poisson-Boltzmann equation,
giving the SOLVATED E-field and EFG. Vacuum Coulomb gives the gas-phase
values. The difference is the solvation contribution. Both are features.
APBS is already computed in Layer 0; vacuum Coulomb from CoulombResult
provides the decomposition APBS cannot.

---

## Calculator 5: Pi-Quadrupole EFG

**Physical picture**: The pi-electron cloud above and below an aromatic
ring creates a quadrupole electric field. This is the ELECTRIC analogue
of ring susceptibility (which is the MAGNETIC analogue).

**Geometric kernel** (per ring, at atom position r):

The EFG from a point axial quadrupole Theta along ring normal n at the
ring center. Derived from Stone's T-tensor formalism (Theory of
Intermolecular Forces, OUP 2013, Ch. 3):

V_ab = -(1/3) Theta_cd T_abcd, where T_abcd = d^4(1/r)/dx_a dx_b dx_c dx_d.

For Theta_cd = Theta(3 n_c n_d - delta_cd)/2, the delta contraction
vanishes by Laplace (T_abcc = 0), giving:

```
G_ab = 105 dn^2 d_a d_b / r^9
     - 30 dn (n_a d_b + n_b d_a) / r^7
     - 15 d_a d_b / r^7
     + 6 n_a n_b / r^5
     + delta_ab (3/r^5 - 15 dn^2/r^7)
```
where d = r - ring_center, r = |d|, dn = d . n (height above ring plane).

The -Theta/2 prefactor is absorbed into the learnable parameter Q_type.

**Properties** (verified on 723 proteins, 1.03M atom-ring pairs):
- Symmetric: G_ab = G_ba (max asymmetry 2.2e-16)
- Traceless: Tr(G) = 0 (Laplace equation; max |Tr| = 5.1e-15)
- Pure T2: T0 = 0, T1 = 0
- Leading radial decay: 1/r^5

**Note on the original formula**: The earlier version of this catalogue
had a formula with terms going as 1/r^7, 1/r^5, 1/r^3. That formula
had Tr = 6(3cos^2 theta - 1)/r^3 != 0, violating Laplace's equation.
The correct formula (above) has terms going as 1/r^9, 1/r^7, 1/r^5,
and is traceless. The error was in the order of the derivatives (the
EFG from a quadrupole involves 4th derivatives of 1/r, not 2nd).

The scalar version (Buckingham A-term, E-field from quadrupole):
```
V_scalar = (3 cos^2 theta - 1) / r^4
```
where theta = angle(d, ring normal). This contributes to T0 via the
Buckingham mechanism sigma_iso = A * E_z. It is stored separately from
the EFG tensor because A and gamma are different learnable parameters.

**T2 structure**: Pure T2 (the tensor is traceless and symmetric).
The T2 pattern differs from the dipolar kernel (1/r^3) because of the
additional ring-normal-dependent terms. T2 independence verified:
PQ vs McConnell 0.38, PQ vs Coulomb 0.40, PQ vs RingChi 0.37,
PQ vs BiotSavart 0.57 (all well below 0.9).

**Per-ring-type pattern**: Five-membered rings (TRP5, HIE) produce
larger signals — atoms approach the smaller ring center more closely,
and the 1/r^5 decay amplifies this. TRP9 (fused perimeter) produces
the weakest signal per pair. This is consistent with the physics:
the quadrupole approximation is better for larger rings where the
field point is farther from the source.

The height_offset parameter shifts the effective quadrupole source
above the ring plane, changing the T2 angular pattern. This is
directly constrained by T2.

---

## Calculator 6: Ring Susceptibility Anisotropy

**Physical picture**: The aromatic ring as a whole has a bulk diamagnetic
susceptibility anisotropy (chi_parallel ≠ chi_perpendicular). Modelled
as a point magnetic dipole at the ring center. The "bond direction"
in the McConnell derivation becomes the RING NORMAL n_hat.

**Full shielding tensor** (same derivation as McConnell with b_hat → n_hat):
```
sigma_ab = prefactor * (Delta_chi_ring / 3) *
    [9 cos_theta * d_hat_a * n_b / r^3
   - 3 n_a n_b / r^3
   - (3 d_hat_a d_hat_b - delta_ab) / r^3]
```
where d = r_atom - ring_center, cos_theta = d_hat . n_hat.

**Tensor structure**: Same as McConnell — asymmetric, non-traceless:
- T0 = (Delta_chi_ring / 3)(3cos^2 theta - 1) / r^3 (non-zero)
- T1 != 0 (from asymmetric d_hat / n_hat coupling)
- T2 != 0 (angular anisotropy)

The only adjustable parameter is Delta_chi_ring (magnitude).
No geometric parameter changes the SHAPE of the angular pattern
(unlike McConnell's midpoint_shift or BS's lobe_offset). The
point-dipole geometry is fully determined by ring center and normal.

**Implication**: If the T2 residual for ring susceptibility is
systematically wrong for a ring type, the point-dipole approximation
itself fails. The fix would be a DISTRIBUTED susceptibility model
(susceptibility spread over ring atoms rather than concentrated at
center). That would be a significant physics finding.

---

## Calculator 7: London Dispersion

**Physical picture**: Van der Waals attraction between the atom and
each ring vertex, via the C6 dispersion coefficient.

**Geometric kernel** (per ring, summed over vertices):
```
K_disp_ab = sum_vertices C6 * (3 d_a d_b / r^8 - delta_ab / r^6)
disp_scalar = sum_vertices C6 / r^6
```
where d = r_atom - r_vertex, r = |d|.

**T2 structure**: Similar to the dipolar kernel but with 1/r^6 decay
(isotropic) and 1/r^8 (anisotropic). The anisotropic part is small
relative to the isotropic part. T2 is present but typically dominated
by McConnell and ring current T2.

Distance filter: R_MIN = 1.5 A to R_CUT = 5.0 A per vertex.

---

## Calculator 8: H-Bond Dipolar

**Physical picture**: The magnetic interaction with an H-bond partner
creates a dipolar shielding contribution. The "bond direction" in the
McConnell derivation becomes the D-H...A direction (or a suitable
unit vector from donor to acceptor).

**Full shielding tensor**: Same derivation as McConnell with the
coupling direction being the H-bond direction h_hat:
```
sigma_ab = prefactor * (eta / 3) *
    [9 cos_theta * d_hat_a * h_b / r^3
   - 3 h_a h_b / r^3
   - (3 d_hat_a d_hat_b - delta_ab) / r^3]
```

The Cornilescu & Bax (1999) isotropic form sigma_iso = cos^2(alpha)/r^3
(where alpha is the D-H...A angle) is the T0 component of this tensor,
with the angular factor being (3cos^2 theta - 1)/r^3 projected along h.

**Tensor structure**: Asymmetric, non-traceless (same as McConnell
and ring susceptibility). T0, T1, T2 all non-zero.

**T2 physics**: Backbone H-bonds (more linear, alpha near 180°)
produce a T2 pattern aligned with the backbone direction. Sidechain
H-bonds (variable geometry) produce more scattered T2 patterns. The
backbone/sidechain weight ratio (w_bb / w_sc) is weakly constrained
by T2 data because the two populations have different angular
distributions.

---

## Calculator 9: MOPAC Coulomb EFG

**Physical picture**: Same as Calculator 4 — partial charges create an
electric field E and electric field gradient tensor V at each atom via
the Coulomb interaction. The difference is the charge source.
CoulombResult uses ff14SB fixed charges (one value per atom type,
invariant across conformations). MopacCoulombResult uses PM7 Mulliken
charges from MOPAC — quantum-mechanically derived charges that respond
to the local electronic environment at each conformation. A backbone
nitrogen near a charged sidechain has a different MOPAC charge than the
same nitrogen far from charged groups. ff14SB cannot see this.

**Geometric kernel**: Identical to Calculator 4:

Electric field (rank-1, Vec3):
```
E_a(i) = ke * sum_{j!=i} q_mopac_j * (r_i - r_j)_a / |r_i - r_j|^3
```

Electric field gradient (rank-2, Mat3, TRACELESS by Gauss's law):
```
V_ab(i) = ke * sum_{j!=i} q_mopac_j * [3 (r_i-r_j)_a (r_i-r_j)_b / |r_i-r_j|^5
                                        - delta_ab / |r_i-r_j|^3]
```

The charges q_mopac_j are different from q_ff14sb_j. Because the charge
distribution differs, the angular pattern of the summed EFG differs —
which atoms dominate the sum at each field point changes. This is not a
re-parameterization of Coulomb. It is a different kernel — the field
that nature produces at this geometry given the actual electron
distribution, not the force field's approximation of it.

**Units**: V/A (E-field), V/A^2 (EFG). Same as CoulombResult.

**T2 structure**: V_ab is traceless symmetric — pure T2 (same as
Calculator 4). The T2 of the MOPAC EFG encodes the angular asymmetry
of the QM charge distribution around the atom. Where ff14SB assigns
uniform charges to all GLY CA atoms regardless of environment, MOPAC
captures the charge redistribution from nearby polar groups, aromatic
rings, and hydrogen bonds.

**T2 independence** (verified on A0A062V9G2, 732 atoms):
MopacCoulomb vs ff14SB Coulomb: |cos| = 0.85. Highly correlated (same
underlying geometry) but not redundant — 15% angular difference from
charge polarisation. The model can distinguish them.

**Decomposition**: Same as Calculator 4:
- Backbone atoms: E_backbone, V_backbone
- Sidechain atoms: E_sidechain, V_sidechain
- Aromatic atoms: E_aromatic, V_aromatic
No APBS solvent subtraction (APBS was solved with ff14SB charges; the
MOPAC-charge solvent field would require a separate PB solve).

**Shielding contribution** (Buckingham 1960, same form as Calculator 4):
```
sigma_iso = A_elem * E_z + B_elem * E_z^2      (T0, ppm)
sigma_T2 = gamma_mopac * V_T2                   (T2, ppm)
```

The model learns gamma_mopac — the Buckingham coupling coefficient for
the QM-charge EFG. The comparison between the ff14SB gamma and MOPAC
gamma weights is itself a finding about whether charge polarisation
matters for the angular structure of shielding.

**Relationship to Calculator 4**: CoulombResult and MopacCoulombResult
compute the same physics with different input data. CoulombResult is the
force field baseline — what you get from fixed charges. MopacCoulombResult
is the QM-derived version — what you get from charges that respond to
conformation. The delta between them isolates the effect of charge
polarisation on the angular EFG pattern. Near a mutation site where
the aromatic ring's removal redistributes electron density, the two
EFGs diverge. That divergence is the charge redistribution signal.

**Dependencies**: MopacResult (conformation electronic structure,
provides mopac_charge on ConformationAtom), SpatialIndexResult.

**Filter set**: SelfSourceFilter (field undefined at source itself).
No spatial cutoff — N^2 sum over all atoms, same as Calculator 4.

---

## Calculator 10: MOPAC Bond-Order-Weighted Anisotropy

**Physical picture**: Same as Calculator 3 — bond magnetic
susceptibility anisotropy creates a dipolar shielding field at nearby
nuclei. The kernel is the same full McConnell tensor M_ab / r^3. The
difference is that each bond's contribution is weighted by its MOPAC
Wiberg bond order — a continuous quantum-mechanical measure of how
much electron density the two bonded atoms share.

**Geometric kernel** (per bond, at atom position r):

```
total_ab = sum_bonds  bond_order_k * M_ab_k / r_k^3
```

where M_ab is the full McConnell tensor from Calculator 3:

```
M_ab = 9 cos_theta * d_hat_a * b_hat_b
     - 3 b_hat_a * b_hat_b
     - (3 d_hat_a * d_hat_b - delta_ab)
```

The bond order is not a parameter. It is a measured property of the
bond at this conformation — part of the geometric reality, like the
charge in Calculator 9. A C=O with bond order 1.85 contributes 1.85×
the kernel of a C=O with order 1.0. Bonds with negligible electron
sharing (order < 0.01) are skipped — they are electronically
insignificant.

**Units**: Angstrom^-3 (same as Calculator 3, before Dchi multiplication).

**What bond-order weighting changes physically**:

In unweighted McConnell, every C=O contributes equally regardless of
its actual electron distribution. In bond-order-weighted McConnell,
double bonds dominate over single bonds. The key consequence: the
RELATIVE contributions of individual bonds to the angular sum change.
A strong C=O near the field atom may now dominate where previously it
was one among many equal contributors. This changes the total T2
angular pattern — the direction of the summed tensor rotates because
different bonds carry different weight.

Near a mutation site where bond character changes (the aromatic ring's
influence on nearby C=O bonds), the weighted and unweighted sums
diverge. This divergence IS the electron redistribution signal in the
bond anisotropy channel.

**Tensor structure** (per bond): Same as Calculator 3 — asymmetric,
non-traceless, T0+T1+T2 all non-zero. The bond-order weighting is a
scalar multiplier on each bond's tensor, so it does not change the
irrep structure of individual terms. It changes the linear combination
of terms in the sum.

**T2 independence** (verified on A0A062V9G2, 732 atoms):
MopacMcConnell vs unweighted McConnell: |cos| = 0.55. Substantially
different angular patterns — well above random (0.36) but well below
redundant (1.0). Bond-order weighting does not merely rescale
McConnell; it reweights which bonds matter, producing genuinely
different angular structure.

**Per-category decomposition**: Same 5 categories as Calculator 3
(PeptideCO, PeptideCN, BackboneOther, SidechainCO, Aromatic). Each
category accumulates the bond-order-weighted kernel sum. The model
learns Delta_chi_mopac per category.

The comparison between unweighted and weighted McConnell at each bond
category answers: does a continuous description of bond character
improve the angular structure beyond categorical Delta_chi? If the
weighted kernel reduces the T2 residual, fixed Delta_chi per category
was the bottleneck — the actual anisotropy varies with electron
distribution, and MOPAC measures that variation.

**Dependencies**: MopacResult (conformation electronic structure,
provides TopologyBondOrder(bond_index) for each covalent bond),
SpatialIndexResult (bond midpoint search within cutoff),
GeometryResult (bond midpoints and directions).

**Filter set**: SelfSourceFilter + DipolarNearFieldFilter (source
extent = bond length). Same as Calculator 3.

**Cutoff**: 10.0 A from bond midpoint (same as Calculator 3).

---

## The T2 Hierarchy (what we expect, to be verified)

The old bilinear model's alpha coefficients suggested a ranking of
T2 contributions. These numbers came from one model on one dataset
and should be treated as HYPOTHESES, not facts. Our correctly-
implemented calculators with full tensors and DFT ground truth will
produce the actual ranking.

**Expected ranking** (from physical reasoning, not old model output):

| Calculator | T2 mechanism | Why it should matter | Key geometric parameter |
|------------|-------------|---------------------|------------------------|
| McConnell | Full dipolar tensor from bonds | Many bonds, steep 1/r^3 decay, large Dchi for C=O | midpoint_shift |
| Coulomb EFG | Charge environment asymmetry | Traceless by Gauss's law — pure T2 | gamma (EFG→T2 coupling) |
| Ring current (BS) | n ⊗ B outer product | Near aromatic rings only | lobe_offset d |
| Ring current (HM) | Surface integral | Same physics as BS, different approximation | (intensity) |
| Ring susceptibility | Full dipolar tensor from ring center | Same kernel as McConnell, fewer sources | (none — pure geometry) |
| Pi-quadrupole | Quadrupole EFG from ring | 1/r^4 decay, ring-normal dependence | height_offset |
| Dispersion | 1/r^8 anisotropy | Small anisotropic contribution | (small) |
| H-bond | Full dipolar tensor from partner | Fewer sources than McConnell | backbone/sidechain weight |
| MOPAC Coulomb | QM-charge EFG asymmetry | Charge polarisation changes angular EFG | gamma_mopac |
| MOPAC McConnell | Bond-order-weighted dipolar | Electron sharing reweights bond contributions | Delta_chi_mopac |

McConnell should dominate T2 because: (a) there are hundreds of
bonds vs a few rings, (b) 1/r^3 decay is steep so nearby bonds
contribute strongly, (c) each bond has a different direction creating
complex angular patterns, (d) C=O bond anisotropy is large.

Whether Coulomb EFG opposes McConnell, and by how much, is an
empirical question. It's physically plausible (charges on C and O
create an EFG that partially opposes the bond's magnetic anisotropy
at nearby atoms) but the magnitude depends on the actual charges
and geometry.

---

## Build Order (from shared kernels)

### Phase 1: The dipolar kernel (one implementation, four users)

Implement K_ab(d) = 3 d_a d_b / r^5 - delta_ab / r^3 once, tested
with hand-calculated values for known geometry. Verify:
- Traceless (T0 = 0)
- Symmetric (K_ab = K_ba)
- Correct T2 components for d = (1, 0, 0): K should be diagonal
  with K_11 = 2/r^3, K_22 = K_33 = -1/r^3
- McConnell scalar f = K contracted with b_hat × b_hat matches
  (3cos^2 theta - 1)/r^3

### Phase 2: McConnellResult

Uses the dipolar kernel from Phase 1. Per-bond, per-category.
The dominant T2 contributor — if this is wrong, everything downstream
is wrong. Verify against the old McConnell calculator output on the
test protein.

### Phase 3: CoulombResult

Same dipolar kernel, charge-weighted sum. Verify EFG is traceless.
Verify decomposition (backbone/sidechain/aromatic). Compare vacuum
EFG to APBS EFG: solvation screening should reduce magnitudes.

### Phase 4: RingSusceptibilityResult and HBondResult

Same dipolar kernel from ring center / H-bond partner. Simple
applications of the verified kernel.

### Phase 5: BiotSavartResult

The JB model. Different kernel (line integral, not dipolar). Port
from old code with exact unit conversions. Verify: proton 3A above
PHE ring → sigma > 0. Compare G tensor to old code output on the
test protein.

### Phase 6: HaighMallionResult

Surface integral. Different kernel (area integral of dipolar).
Verify against BS at large distance (should converge). Verify T2
differs from BS at close range (the opposing-T2 finding).

### Phase 7: PiQuadrupoleResult and DispersionResult

Smaller effects. Build last.

---

## Verification Strategy

For each calculator, two levels of verification:

### Level 1: Analytical test geometry

Choose a simple geometry where the answer is known analytically.
Example for dipolar kernel: atom at (3, 0, 0), source at origin.
K should be diag(2/(3^3), -1/(3^3), -1/(3^3)) = diag(0.0741, -0.0370, -0.0370).

### Level 2: Cross-check against old code

Run both old and new calculators on A0A7C5FAR6 WT. Compare Mat3
outputs element by element. They should match within floating point
tolerance. Differences indicate bugs.

### Level 3: Physics check against DFT

The test pair A0A7C5FAR6 has FOUR aromatic→ALA mutations:
TRP3, TYR10, HIE12, PHE31 — all aromatic residues removed.

WT: 543 atoms, 35 residues (MET LYS TRP ARG CYS ... HIE ... PHE ...)
ALA: 501 atoms, 35 residues (MET LYS ALA ARG CYS ... ALA ... ALA ...)

This means: the WT-ALA delta isolates the TOTAL ring current
contribution from all four aromatic rings simultaneously. It does
NOT isolate single ring types. For ring-type-specific analysis,
single-site mutants would be needed (from the larger 560-pair set).

The protonation state differs: WT has HIE (epsilon-protonated HIS),
ALA mutant has no histidine. The delta includes the protonation
change — this is physical reality, not an artefact.

Atom matching: backbone atoms match 1:1 across all 35 residues.
Sidechain atoms match at non-mutated residues (31 of 35).
At the 4 mutated residues, only the backbone atoms correspond.
Heavy-atom position matching (not name matching) is required.

---

## RESOLVED: McConnell Full Shielding Tensor Derivation

The dipolar kernel K_ab = 3 d_a d_b / r^5 - delta_ab / r^3 is
traceless. But McConnell produces measurable isotropic shifts.
The spec's formula sigma_ab = Delta_chi * K_ab gives T0 = 0.

This is because Delta_chi * K is the GEOMETRIC KERNEL, not the
full shielding tensor. The full tensor comes from the magnetic
dipole interaction, derived here from first principles.

### Derivation

A bond with axial susceptibility anisotropy Delta_chi along
direction b_hat in external field B_0 acquires induced moment
m = chi . B_0. The susceptibility tensor:

```
chi_cd = chi_bar delta_cd + (Delta_chi / 3)(3 b_c b_d - delta_cd)
```

The shielding tensor sigma_ab = -dB_a^sec / dB_0b gives:

```
sigma_ab = prefactor * Sum_c chi_cb * K_ca
```

where K_ca = (3 d_hat_c d_hat_a - delta_ca) / r^3.

The isotropic susceptibility (chi_bar) term gives chi_bar * K_ab,
which is traceless. No T0 from an isotropic magnetic sphere.

The anisotropic term (Delta_chi) expands to:

```
sigma_ab = prefactor * (Delta_chi / 3) *
    [9 cos_theta * d_hat_a * b_b / r^3
   - 3 b_a b_b / r^3
   - (3 d_hat_a d_hat_b - delta_ab) / r^3]
```

where cos_theta = d_hat . b_hat.

### Trace verification

```
Tr(sigma) = prefactor * (Delta_chi / 3) *
    [9 cos^2_theta / r^3 - 3 / r^3 - 0]
  = prefactor * Delta_chi * (3 cos^2_theta - 1) / r^3
```

T0 = Tr(sigma) / 3 = prefactor * (Delta_chi / 3) * (3 cos^2_theta - 1) / r^3.

This IS the McConnell formula. T0 is non-zero. Mystery solved.

### Key finding: the tensor is ASYMMETRIC

The term 9 cos_theta * d_hat_a * b_b has d_hat in position a and
b_hat in position b. Swapping a, b gives a different value.

The full McConnell shielding tensor is NOT symmetric. It has:
- T0 != 0 (the McConnell isotropic shift)
- T1 != 0 (from asymmetric coupling of d_hat and b_hat)
- T2 != 0 (the angular anisotropy)

This is physically correct: magnetic interactions CAN produce
asymmetric shielding tensors (the antisymmetric part relates to
paramagnetic contributions and exists in the DFT output too).

### What this means for implementation

1. The old code's K_ab (symmetric, traceless) is INCOMPLETE. It
   captures only part of the T2 and none of the T0 or T1.

2. The full formula per bond is:
   ```
   sigma_ab = (Delta_chi / 3) * [9 cos_theta * d_hat_a * b_hat_b
              - 3 b_hat_a * b_hat_b
              - (3 d_hat_a * d_hat_b - delta_ab)] / r^3
   ```
   (with appropriate unit/sign prefactor)

3. This gives mc_shielding_contribution with non-zero T0, T1, T2.
   McConnell WILL reduce the T0 residual. McConnell WILL contribute
   to T1. The T2 is richer than just the dipolar kernel.

4. The McConnell scalar f = (3 cos^2_theta - 1) / r^3 is correctly
   derived as T0 * 3 / Delta_chi. No separate formula needed — it
   falls out of the tensor trace.

5. Both the FULL tensor and the symmetric dipolar kernel K should
   be stored. K is still useful as a feature (the pure geometric
   part). The full tensor is what the residual subtraction needs.

### What this means for other dipolar calculators

Ring susceptibility and H-bond also use the dipolar interaction,
but with different source properties:

- Ring susceptibility: the "bond direction" is the ring normal.
  Same derivation applies — the full tensor is asymmetric with
  non-zero T0.

- H-bond: the coupling direction is the D-H...A vector.
  Same structure applies.

- Coulomb EFG: different physics (electric, not magnetic). The EFG
  tensor V_ab = Sum q_j * K_ab(j) IS symmetric and traceless
  (Gauss's law). No asymmetric term. Coulomb T0 comes from the
  Buckingham conversion (A * E_z + B * E_z^2), not from the EFG trace.

So: dipolar MAGNETIC calculators (McConnell, ring susceptibility,
H-bond) produce asymmetric non-traceless tensors. The ELECTRIC
calculator (Coulomb EFG) produces symmetric traceless tensors.
Different physics, different tensor structure.
