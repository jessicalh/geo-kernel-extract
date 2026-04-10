# Calculator Parameter API

How tuneable parameters enter each classical calculator. Each calculator
has a typed parameter struct with literature defaults. The calibration
pipeline (Python e3nn in learn/c_equivariant/) discovers optimal values
by training against DFT WT-ALA delta tensors across 723 proteins.
Calibrated values enter the C++ system as TOML configuration overriding
the literature defaults.

This document defines the parameters, their equations, and their
physical meaning for all 8 classical calculators (+ 2 MOPAC-derived).

Aligned with CONSTITUTION.md (2026-03-29 revision) and OBJECT_MODEL.md.

---

## Architecture Overview

Each calculator defines a typed parameter struct (e.g., BiotSavartParams).
Literature defaults are the starting point. TOML configuration can
override any parameter. The calculator does not know or care where the
values came from — it computes the geometric kernel with whatever
parameters it receives.

The calibration loop:
```
C++ geometric kernels → NPY features (via WriteFeatures)
→ Python e3nn calibration against DFT WT-ALA deltas
→ TOML parameter file → C++ (next extraction run)
```

The typed parameter interface is designed to support parameter delivery
from any source — currently TOML configuration from the external
calibration pipeline, potentially an internal e3nn model
(ParameterCorrectionResult) in future. This design discipline ensures
full tensor output and T2 completeness: every calculator MUST produce
the same typed output regardless of where its parameters came from.

All code paths produce identical output types: full Mat3 + SphericalTensor
per atom, per source. Changing parameters does not change the output
structure — only the tensor values.

---

## CalculatorParameters: The Type

Each calculator defines a typed parameter struct. This is NOT a generic
map. It is a typed, compile-time-checked struct specific to each
calculator.

```cpp
struct BiotSavartParams {
    // Per-calculator-specific fields with literature defaults
    array<double, 8> intensity;   // from ring type class
    array<double, 8> lobe_offset; // from ring type class
    double stacking_scale = 1.0;
};
```

Each calculator has a `DefaultParameters()` method returning literature
values with references. TOML configuration can override any field.

Parameters are GLOBAL (per ring type, per bond category, per element),
not per-atom. The calibration pipeline discovers that "PHE ring current
intensity should be -11.3 nA/T instead of -12.0 nA/T" — a single
number that applies to every PHE ring in every protein. Per-atom
variation is handled by the geometry, not by per-atom parameter tuning.

---

## Calculator 1: Biot-Savart Ring Current (BiotSavartResult)

### The equation

Johnson-Bovey double-loop model. For atom at position r relative to
ring with center c, normal n, radius R, and pi-electron lobe offset d:

```
B(r) = (mu_0 * I) / (4 * pi) * sum_loop integral[ dl x r' / |r'|^3 ]
```

where the two loops are at z = +d and z = -d above and below the ring
plane.

The geometric kernel (separated from intensity):

```
G_ab(r) = -(1/I) * n_b * B_a(r) * PPM_FACTOR
```

The shielding contribution from ring type t:

```
sigma_t = I_t * G_T0_t          (isotropic, scalar)
T2_t    = I_t * G_T2_t          (anisotropic, 5-component)
```

Total ring current shielding:

```
sigma_RC = sum_t sum_rings_of_type_t I_t * G(r, ring)
```

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| PHE benzene intensity | I_PHE | -12.0 | nA/T | Giessner-Prettre & Pullman 1969 |
| TYR phenol intensity | I_TYR | -11.28 | nA/T | Giessner-Prettre & Pullman 1969 |
| TRP benzene intensity | I_TRP6 | -12.48 | nA/T | Giessner-Prettre & Pullman 1969 |
| TRP pyrrole intensity | I_TRP5 | -6.72 | nA/T | Giessner-Prettre & Pullman 1969 |
| TRP perimeter intensity | I_TRP9 | -19.2 | nA/T | Estimated from TRP5 + TRP6 |
| HIS imidazole intensity | I_HIS | -5.16 | nA/T | Giessner-Prettre & Pullman 1969 |
| HID imidazole intensity | I_HID | -5.16 | nA/T | Same as HIS (sparse data) |
| HIE imidazole intensity | I_HIE | -5.16 | nA/T | Same as HIS (sparse data) |
| PHE lobe offset | d_PHE | 0.64 | Angstroms | Johnson & Bovey 1958 |
| TYR lobe offset | d_TYR | 0.64 | Angstroms | Johnson & Bovey 1958 |
| TRP6 lobe offset | d_TRP6 | 0.64 | Angstroms | Johnson & Bovey 1958 |
| TRP5 lobe offset | d_TRP5 | 0.52 | Angstroms | Johnson & Bovey 1958 |
| TRP9 lobe offset | d_TRP9 | 0.60 | Angstroms | Estimated |
| HIS lobe offset | d_HIS | 0.50 | Angstroms | Johnson & Bovey 1958 |
| HID lobe offset | d_HID | 0.50 | Angstroms | Johnson & Bovey 1958 |
| HIE lobe offset | d_HIE | 0.50 | Angstroms | Johnson & Bovey 1958 |

Total: 16 parameters (8 intensities + 8 lobe offsets).

Ring current is element-INDEPENDENT (same B-field regardless of
observing nucleus). These 16 parameters apply to ALL atoms.

### TOML-configurable parameters

```cpp
struct BiotSavartParams {
    bool has_corrections = false;

    // Corrected ring current intensities (nA/T)
    // Index by RingTypeIndex (0..7)
    array<double, 8> intensity;  // default: from ring type class

    // Corrected JB lobe offsets (Angstroms)
    array<double, 8> lobe_offset;  // default: from ring type class

    // Multiplicative stacking correction (dimensionless)
    // Applied when multiple rings are within 8A of each other.
    // Accounts for mutual ring polarisation (Haigh-Mallion 1979).
    // Default: 1.0 (no correction)
    double stacking_scale;
};
```

The calibration pipeline discovers these by fitting against DFT WT-ALA
delta tensors. The pipeline sees the geometric kernels G and learns
what intensity values best explain the data. Calibrated values override
the literature defaults via TOML configuration.

### How the calculator uses it

The calculator reads its parameter struct (from TOML or defaults) and
computes the kernel. It does not know where the values came from:

```cpp
auto params = LoadBiotSavartParams(toml_config);  // or DefaultParams()

for (auto& ring : conf.AllRings()) {
    int type = ring.TypeIndexAsInt();
    double I = params.intensity[type];
    double d = params.lobe_offset[type];

    for (auto& atom : conf.AtomsNear(ring.center, 15.0)) {
        Vec3 B = JohnsonBoveyBField(atom.position, ring, d);
        Mat3 G = GeometricKernel(B, ring.normal, I);
        // ... store G_tensor, G_spherical, B_field ...
    }
}
```

### Measuring calibration effect

Running the pipeline with default parameters and again with calibrated
parameters produces a per-atom delta:

```
delta_BS[atom] = BS_calibrated[atom].T0 - BS_default[atom].T0
```

If zero everywhere, the literature values are optimal. If consistently
nonzero for a specific ring type, the calibration pipeline has
discovered a systematic correction to that ring type's intensity.

### Can T2 corrections feed back into the calculation?

**Partially, through lobe offset.** The T2 output of the calibration pipeline
encodes the angular pattern of shielding around each ring. The JB lobe
offset d controls the angular shape of the B-field (higher d = more
concentrated lobes, sharper angular dependence). If the model's T2
prediction for atoms near a ring systematically shows sharper angular
variation than the default d produces, this constrains d.

The constraint is indirect: the calibration model outputs a corrected d (L=0 scalar),
not a T2-dependent d. But the calibration model learned d by fitting its T2 predictions
to the DFT T2, so the T2 physics IS encoded in the corrected d value.

**Direct T2 feedback is not computationally useful here.** The Biot-Savart
equation has a fixed functional form: B(r) from a current loop. The
angular pattern of B is determined entirely by the geometry (r, d, R)
and cannot be independently modulated by a T2 tensor without breaking
the physical model. The T2 output is a prediction to be compared with
the calculator's T2 output, not an input to the calculator.

---

## Calculator 2: Haigh-Mallion Surface Integral (HaighMallionResult)

### The equation

Haigh-Mallion (1973) surface integral over the ring polygon:

```
sigma_HM = I_HM * integral_S [ (1/|r-r'|^3) * (r-r') . dS ]
```

Implemented as 7-point Gaussian quadrature on fan triangulation of the
ring polygon. Produces a full rank-2 tensor (NOT rank-1 like Biot-Savart).

The key difference from Biot-Savart: HM uses a surface integral that
naturally produces rank-2 anisotropy. BS uses a line integral that
produces rank-1 (dipolar) anisotropy. The two make opposing T2
predictions at the same geometry (INTERFACE_AND_PHYSICS_NOTES finding).

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| PHE HM intensity | J_PHE | -12.0 | nA/T | Same as BS (Haigh & Mallion 1979) |
| TYR HM intensity | J_TYR | -11.28 | nA/T | Same as BS |
| TRP6 HM intensity | J_TRP6 | -12.48 | nA/T | Same as BS |
| TRP5 HM intensity | J_TRP5 | -6.72 | nA/T | Same as BS |
| TRP9 HM intensity | J_TRP9 | -19.2 | nA/T | Same as BS |
| HIS HM intensity | J_HIS | -5.16 | nA/T | Same as BS |
| HID HM intensity | J_HID | -5.16 | nA/T | Same as BS |
| HIE HM intensity | J_HIE | -5.16 | nA/T | Same as BS |
| Quadrature order | N_quad | 7 | - | Haigh & Mallion 1979 |

Total: 9 parameters (8 intensities + 1 quadrature order).

Note: HM intensities are traditionally set equal to BS intensities.
The calibration pipeline can learn DIFFERENT HM intensities, which would
mean the surface integral model benefits from different effective
current strengths than the line integral model. This is physically
reasonable: the two models approximate different aspects of the
actual current distribution.

### TOML-configurable parameters

```cpp
struct HaighMallionParams {
    bool has_corrections = false;

    // Corrected HM ring current intensities (nA/T)
    array<double, 8> intensity;  // default: same as BS defaults

    // BS/HM blending weight (dimensionless, 0.0 to 1.0)
    // At what distance does HM dominate BS?
    // Default: not used (both run independently)
    // The calibration pipeline can learn that for atoms within X Angstroms,
    // HM is more accurate, and beyond X, BS is more accurate.
    // This is expressed as a distance-decay crossover length.
    double crossover_distance;  // Angstroms, default: 0.0 (disabled)
};
```

### How the calculator uses it

```cpp
auto params = LoadHaighMallionParams(toml_config);  // or DefaultParams()

for (auto& ring : conf.AllRings()) {
    int type = ring.TypeIndexAsInt();
    double J = params.intensity[type];

    for (auto& atom : conf.AtomsNear(ring.center, 15.0)) {
        Mat3 hm_tensor = HaighMallionSurfaceIntegral(
            atom.position, ring.geometry, N_QUAD);
        hm_tensor *= J;
        // ... store hm_tensor, hm_spherical, hm_B_field ...
    }
}
```

### Can T2 corrections feed back?

**Yes, indirectly through intensity separation.** The key finding from
INTERFACE_AND_PHYSICS_NOTES is that BS and HM make opposing T2
predictions. The calibration pipeline, by learning separate I_BS and I_HM
intensities, effectively learns the optimal blend of the two models'
T2 patterns. If the calibration model sets I_HM higher than I_BS for a ring type,
it is saying "the surface integral's angular pattern is more accurate
than the line integral's for this ring type."

**Direct T2 feedback is not possible.** The HM surface integral has a
fixed functional form determined by geometry. You cannot inject a T2
correction into the integrand.

---

## Calculator 3: McConnell Bond Anisotropy (McConnellResult)

### The equation

McConnell (1957) magnetic susceptibility anisotropy of covalent bonds:

```
T_ab = Delta_chi * (3 * d_a * d_b / r^5 - delta_ab / r^3)
```

where:
- Delta_chi is the magnetic susceptibility anisotropy of the bond
- d = r_atom - r_bond_midpoint (displacement vector)
- r = |d| (distance to bond midpoint)

The McConnell scalar factor (derived FROM the tensor trace):

```
f_McConnell = (3 * cos^2(theta) - 1) / r^3
```

where theta is the angle between the displacement vector and the
bond axis.

Category-decomposed:

```
sigma_bond = sum_bonds [ Delta_chi_cat * f_McConnell(r, theta) ]
```

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| Peptide C=O anisotropy | Delta_chi_CO | -13.0e-6 | cm^3/mol | Flygare 1971 |
| Peptide C-N anisotropy | Delta_chi_CN | -5.5e-6 | cm^3/mol | Flygare 1971 |
| Sidechain C=O anisotropy | Delta_chi_SCO | -13.0e-6 | cm^3/mol | Same as backbone CO |
| Aromatic C-C anisotropy | Delta_chi_arom | -8.0e-6 | cm^3/mol | Estimated |
| Sidechain other anisotropy | Delta_chi_other | -2.0e-6 | cm^3/mol | Estimated |

Total: 5 parameters (one per BondCategory that has anisotropy).

Note: Delta_chi values from gas-phase measurements on small molecules.
Protein environment may modify effective anisotropy.

### TOML-configurable parameters

```cpp
struct McConnellParams {
    bool has_corrections = false;

    // Corrected bond magnetic susceptibility anisotropy (cm^3/mol)
    // Indexed by BondCategory:
    //   [0] PeptideCO, [1] PeptideCN, [2] BackboneOther,
    //   [3] SidechainCO, [4] Aromatic, [5] SidechainOther
    array<double, 6> delta_chi;

    // Effective bond length correction (dimensionless multiplier on r)
    // Accounts for the fact that the effective dipole center is not
    // exactly at the bond midpoint. Default: 1.0 (midpoint).
    // Values < 1.0 move the effective center toward atom A.
    // Values > 1.0 move it toward atom B.
    // Per-category, since peptide CO has an asymmetric electron density.
    array<double, 6> midpoint_shift;  // default: all 1.0
};
```

### How the calculator uses it

```cpp
auto params = LoadMcConnellParams(toml_config);  // or DefaultParams()

for (auto& bond : conf.AllBonds()) {
    int cat = static_cast<int>(bond.category);
    double dchi = params.delta_chi[cat];
    double shift = params.midpoint_shift[cat];

    Vec3 midpoint = bond.Midpoint(positions);
    if (shift != 1.0) {
        Vec3 dir = bond.Direction(positions);
        midpoint += (shift - 1.0) * 0.5 * bond.Length(positions) * dir;
    }

    for (auto& atom : conf.AtomsNear(midpoint, 10.0)) {
        Vec3 d = atom.position - midpoint;
        double r = d.norm();
        Mat3 T = dchi * (3.0 * d * d.transpose() / pow(r, 5)
                       - Mat3::Identity() / pow(r, 3));
        // ... store T, SphericalTensor, mcconnell_scalar ...
    }
}
```

### Can T2 corrections feed back?

**Yes, through midpoint_shift.** The McConnell equation's T2 pattern is
controlled by where the effective dipole center sits along the bond.
If the calibration pipeline's T2 prediction for atoms near a C=O bond shows
a systematic tilt relative to what the midpoint geometry predicts, the
model can encode this as a midpoint_shift parameter. This changes the
angular pattern of the (3cos^2 theta - 1) / r^3 term without changing
the equation's form.

**Also through Delta_chi itself.** Different Delta_chi values produce
different T2 magnitudes. The T2 output gives the calibration model a stronger
constraint on Delta_chi than T0 alone, because T2 depends on both
the magnitude AND the angular distribution of the anisotropy.

This is the calculator where T2 feedback is MOST useful. McConnell
was identified as the dominant T2 contributor (|alpha| = 3.61 in
the INTERFACE_AND_PHYSICS_NOTES analysis), meaning the calibration model has
the most T2 signal to learn from here.

---

## Calculator 4: Coulomb Electric Field Gradient (CoulombResult)

### The equation

Buckingham (1960) electric field effect on NMR shielding:

```
sigma_Buckingham = A_elem * E_z + B_elem * E_z^2
```

where E_z is the E-field component along the atom's primary bond axis
(for H atoms: the X-H bond direction).

For the full EFG tensor:

```
E_a = sum_j [ q_j * (r_j - r_i)_a / |r_j - r_i|^3 ]
V_ab = sum_j [ q_j * (3*(r_j-r_i)_a*(r_j-r_i)_b / |r_j-r_i|^5
                      - delta_ab / |r_j-r_i|^3) ]
```

T2 coupling:

```
T2_Buckingham = gamma_elem * EFG_T2
```

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| H linear Buckingham | A_H | -1.14 | ppm/(V/A) | Case 1995 (unit-converted) |
| C linear Buckingham | A_C | -2.28 | ppm/(V/A) | 2x A_H (polarisability ratio) |
| N linear Buckingham | A_N | -3.42 | ppm/(V/A) | 3x A_H |
| O linear Buckingham | A_O | -5.70 | ppm/(V/A) | 5x A_H |
| S linear Buckingham | A_S | -9.12 | ppm/(V/A) | 8x A_H |
| H quadratic Buckingham | B_H | -0.111 | ppm/(V/A)^2 | Buckingham 1960 |
| C quadratic Buckingham | B_C | -0.222 | ppm/(V/A)^2 | Estimated 2x B_H |
| N quadratic Buckingham | B_N | -0.333 | ppm/(V/A)^2 | Estimated 3x B_H |
| O quadratic Buckingham | B_O | -0.555 | ppm/(V/A)^2 | Estimated 5x B_H |
| S quadratic Buckingham | B_S | -0.888 | ppm/(V/A)^2 | Estimated 8x B_H |
| H EFG coupling | gamma_H | 0.5 | ppm/(V/A^2) | Estimated |
| C EFG coupling | gamma_C | 1.0 | ppm/(V/A^2) | Estimated |
| N EFG coupling | gamma_N | 1.5 | ppm/(V/A^2) | Estimated |
| O EFG coupling | gamma_O | 2.5 | ppm/(V/A^2) | Estimated |
| S EFG coupling | gamma_S | 4.0 | ppm/(V/A^2) | Estimated |

Total: 15 parameters (5 A + 5 B + 5 gamma, all element-dependent).

Note: Buckingham IS element-dependent (unlike ring current). The
atomic polarisability and paramagnetic susceptibility determine how
strongly the electron cloud responds to the external E-field. The
ratios above follow chemical reasoning but are approximate. The
parameter correction model refines these.

The Coulomb calculator itself (charge summation) has no tuneable
parameters. The charges come from ChargeAssignmentResult (ff14SB)
or MopacResult (PM7 Mulliken). What the calibration pipeline tunes is how the E-field/EFG
TRANSLATES into shielding perturbation.

### TOML-configurable parameters

```cpp
struct CoulombParams {
    bool has_corrections = false;

    // Element-dependent Buckingham linear coefficients (ppm/(V/A))
    // Index by Element enum: H=0, C=1, N=2, O=3, S=4
    array<double, 5> A;

    // Element-dependent Buckingham quadratic coefficients (ppm/(V/A)^2)
    array<double, 5> B;

    // Element-dependent EFG->T2 coupling (ppm/(V/A^2))
    array<double, 5> gamma;

    // Charge scaling factor (dimensionless)
    // If the calibration pipeline finds that ff14SB charges systematically
    // over/underestimate the effective E-field, this global scale
    // corrects for it. Default: 1.0
    double charge_scale;

    // Spatial cutoff (Angstroms). Default: 20.0
    // E ~ 1/r^2, EFG ~ 1/r^3. At 20A with 0.5e charge, E = 0.018 V/A.
    // Uses SpatialIndexResult::AtomsWithinRadius for the inner loop
    // instead of the N^2 all-pairs sum. Set to 999 for full sum.
    // TOML key: coulomb_efield_cutoff
    double efield_cutoff;  // default: 20.0
};
```

### How the calculator uses it

The CoulombResult itself computes E and EFG from raw charges and
geometry. It does NOT apply Buckingham coefficients. The Buckingham
translation happens at FeatureExtractionResult time, using the
CoulombResult's stored E-field and the CoulombParams:

```cpp
// At feature extraction time:
Vec3 E = conf.Coulomb().EFieldAt(atom_idx);
Mat3 EFG = conf.Coulomb().EFGAt(atom_idx);
Element elem = protein.AtomAt(atom_idx).element;
int e = static_cast<int>(elem);

double E_z = E.dot(bond_direction);

// Buckingham T0:
double sigma_buck = params.A[e] * E_z + params.B[e] * E_z * E_z;

// Buckingham T2:
SphericalTensor efg_sph = SphericalTensor::Decompose(EFG);
array<double,5> T2_buck;
for (int m = 0; m < 5; m++) T2_buck[m] = params.gamma[e] * efg_sph.T2[m];
```

This separation means: the CoulombResult stores physics (E-field, EFG),
the Buckingham parameters convert physics to shielding. The two are
independent. Different Buckingham parameters do not require recomputing
the E-field.

### Can T2 corrections feed back?

**Yes, directly through gamma.** The gamma parameter maps EFG (a rank-2
tensor from the physics) to T2 shielding (a rank-2 observable). The
parameter correction model's T2 prediction provides a direct constraint on gamma:

```
T2_predicted = gamma_elem * EFG_T2_computed
```

If the calibration model predicts a specific T2 pattern at an atom where the
computed EFG_T2 is known, gamma is determined by the ratio. This is
the most direct T2-to-parameter feedback in the entire system.

**Also through charge_scale.** If the calibration pipeline finds that the
ff14SB E-field is systematically too large (or too small), the
charge_scale corrects this. The T2 pattern is preserved (EFG geometry
does not change with scaling), but the magnitude changes.

---

## Calculator 5: Pi-Quadrupole EFG (PiQuadrupoleResult)

### The equation

Second derivative of the axial quadrupole potential from pi-electron
density above and below the ring:

```
V_quad = Q * (3*cos^2(theta) - 1) / r^4
```

Full EFG tensor:

```
V_ab = Q * [ 15*dn^2*d_a*d_b/r^7 - 3*dn*(n_a*d_b + n_b*d_a)/r^5
           - 3*d_a*d_b/r^5 + delta_ab*(3*dn^2/r^5 - 1/r^3) ]
```

where:
- Q is the ring quadrupole moment
- dn = d . n (displacement projected onto ring normal)
- d = r_atom - r_center
- r = |d|

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| PHE quadrupole moment | Q_PHE | -8.5 | Buckingham (D*A) | Battaglia et al. 1981 |
| TYR quadrupole moment | Q_TYR | -7.5 | D*A | Estimated (OH reduces) |
| TRP6 quadrupole moment | Q_TRP6 | -9.0 | D*A | Estimated |
| TRP5 quadrupole moment | Q_TRP5 | -4.0 | D*A | Estimated |
| TRP9 quadrupole moment | Q_TRP9 | -13.0 | D*A | Estimated |
| HIS quadrupole moment | Q_HIS | -3.5 | D*A | Estimated (2N reduces) |
| HID quadrupole moment | Q_HID | -3.5 | D*A | Same as HIS |
| HIE quadrupole moment | Q_HIE | -3.5 | D*A | Same as HIS |

Total: 8 parameters (one per ring type).

Note: experimental quadrupole moments for isolated amino acid rings
are sparse. Literature values are for benzene (Q = -8.48 D*A from
Battaglia et al.). Protein ring values are estimated by analogy.
These are prime candidates for parameter correction model correction.

### TOML-configurable parameters

```cpp
struct PiQuadrupoleParams {
    bool has_corrections = false;

    // Corrected quadrupole moments per ring type (D*A)
    array<double, 8> Q;  // default: from literature estimates above

    // Effective height of quadrupole above ring plane (Angstroms)
    // The axial quadrupole model assumes the charge distribution
    // is exactly in the ring plane. Real pi-electrons extend above
    // and below. This parameter shifts the effective source.
    // Default: 0.0 (in-plane approximation)
    array<double, 8> height_offset;
};
```

### Can T2 corrections feed back?

**Yes, through height_offset.** The T2 pattern of the quadrupole EFG
depends sensitively on the effective height of the quadrupole source.
An in-plane quadrupole (height=0) produces a different angular T2
pattern than an offset quadrupole (height=0.3A). The calibration pipeline's
T2 prediction constrains this parameter.

**Also through Q magnitude.** The T2 magnitude is proportional to Q.
But T0 is also proportional to Q, so Q is already constrained by T0.
The additional information from T2 is the ANGULAR pattern, which
constrains height_offset independently of Q.

---

## Calculator 6: Ring Susceptibility Anisotropy (RingSusceptibilityResult)

### The equation

McConnell (1957) applied to the whole ring's magnetic susceptibility
anisotropy, treated as a point dipole at the ring center:

```
T_ab = Delta_chi_ring * (3*d_a*d_b/r^5 - delta_ab/r^3)
```

McConnell scalar:

```
f = (3*cos^2(theta) - 1) / r^3
```

This is the same functional form as bond McConnell but applied from
the ring center rather than from each bond midpoint. It captures the
bulk diamagnetic anisotropy of the ring.

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| PHE ring chi anisotropy | Dchi_PHE | -94.0e-6 | cm^3/mol | Flygare & Benson 1971 |
| TYR ring chi anisotropy | Dchi_TYR | -85.0e-6 | cm^3/mol | Estimated |
| TRP6 ring chi anisotropy | Dchi_TRP6 | -97.0e-6 | cm^3/mol | Estimated |
| TRP5 ring chi anisotropy | Dchi_TRP5 | -55.0e-6 | cm^3/mol | Estimated |
| TRP9 ring chi anisotropy | Dchi_TRP9 | -152.0e-6 | cm^3/mol | Sum TRP5+TRP6 |
| HIS ring chi anisotropy | Dchi_HIS | -40.0e-6 | cm^3/mol | Estimated |
| HID ring chi anisotropy | Dchi_HID | -40.0e-6 | cm^3/mol | Same as HIS |
| HIE ring chi anisotropy | Dchi_HIE | -40.0e-6 | cm^3/mol | Same as HIS |

Total: 8 parameters (one per ring type).

Note: the benzene experimental value is -94.0e-6 cm^3/mol. Protein
ring values are scaled by ring size and electron count.

### TOML-configurable parameters

```cpp
struct RingSusceptParams {
    bool has_corrections = false;

    // Corrected ring susceptibility anisotropy per ring type (cm^3/mol)
    array<double, 8> delta_chi;  // default: from literature above
};
```

This calculator has the simplest parameter correction: one scalar per
ring type. No geometric parameters to tune (the geometry is fully
determined by the ring center and atom position).

### Can T2 corrections feed back?

**Only through Delta_chi magnitude.** The dipolar tensor from a point
source has a fixed angular pattern: (3cos^2 theta - 1)/r^3 is purely
geometric. Delta_chi scales the magnitude uniformly. There is no
parameter that changes the angular pattern.

However, if the calibration pipeline's T2 prediction systematically deviates
from the point-dipole angular pattern, this signals that the point-dipole
approximation is inadequate for that ring type. The correction would be
to extend the physical model (e.g., use a distributed dipole over the ring atoms
rather than a point at the center), not to tune a parameter.

---

## Calculator 7: London Dispersion (DispersionResult)

### The equation

Per-vertex van der Waals dispersion interaction:

```
T_disp_ab = sum_vertices [ C6 * (3*d_a*d_b/r^8 - delta_ab/r^6) ]
disp_scalar = sum_vertices [ C6 / r^6 ]
```

where:
- C6 is the dispersion coefficient (from VdW parameters)
- d = r_atom - r_vertex
- r = |d|
- Distance filter: R_MIN = 1.5A to R_CUT = 5.0A

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| C6 coefficient (C-C) | C6_CC | 1.75e-3 | kcal*A^6/mol | AMBER ff14SB |
| C6 coefficient (C-N) | C6_CN | 1.50e-3 | kcal*A^6/mol | AMBER ff14SB |
| C6 coefficient (C-H) | C6_CH | 0.90e-3 | kcal*A^6/mol | AMBER ff14SB |
| Minimum distance | R_MIN | 1.5 | Angstroms | Constitution |
| Maximum distance | R_CUT | 5.0 | Angstroms | Constitution |
| Dispersion-to-shielding coupling | alpha_disp | 0.1 | ppm*A^6 | Estimated |

Total: 6 parameters.

Note: C6 coefficients come from the force field VdW parameters.
The tuneable physics is alpha_disp: how much does a given 1/r^6
dispersion contact actually perturb the shielding? This is poorly
known experimentally.

### TOML-configurable parameters

```cpp
struct DispersionParams {
    bool has_corrections = false;

    // Dispersion-to-shielding coupling factor (ppm*A^6)
    // Element-dependent: dispersion affects different nuclei differently
    array<double, 5> alpha_disp;  // index by Element, default: 0.1

    // Corrected distance cutoffs (Angstroms)
    double r_min;  // default: 1.5
    double r_cut;  // default: 5.0
};
```

### Can T2 corrections feed back?

**Minimally.** Dispersion is primarily isotropic (1/r^6 dominates
over the anisotropic 3*d_a*d_b/r^8 - delta_ab/r^6 term). The T2
component exists but is typically small compared to other calculators.
The calibration pipeline's T2 prediction near aromatic rings is dominated by
ring current and McConnell T2, making it difficult to isolate the
dispersion T2 contribution.

The alpha_disp parameter is mainly constrained by T0, not T2.

---

## Calculator 8: H-Bond Dipolar (HBondResult)

### The equation

Dipolar tensor from H-bond geometry:

```
T_ab = eta * (3*d_a*d_b/r^5 - delta_ab/r^3)
```

Isotropic contribution:

```
sigma_hbond = eta * cos^2(alpha) / r^3
```

where:
- eta is the H-bond shielding coupling constant
- d = r_atom - r_partner (displacement to H-bond partner)
- r = |d|
- alpha is the D-H...A angle

### Tuneable parameters

| Parameter | Symbol | Default | Unit | Reference |
|-----------|--------|---------|------|-----------|
| H-bond coupling constant | eta | 2.0 | ppm*A^3 | Barfield & Karplus 1969 |
| Distance exponent | n | 3 | - | Classical dipolar (1/r^3) |
| Angular function shape | cos_power | 2 | - | cos^2(alpha) |
| Maximum H-bond distance | d_max | 3.5 | Angstroms | HBOND_COUNT_RADIUS |
| Backbone H-bond weight | w_bb | 1.0 | dimensionless | Default: no distinction |
| Sidechain H-bond weight | w_sc | 1.0 | dimensionless | Default: no distinction |

Total: 6 parameters.

Note: the classical form uses cos^2(alpha)/r^3 but the actual
distance and angular dependence may differ. The calibration pipeline can
learn effective exponents.

### TOML-configurable parameters

```cpp
struct HBondParams {
    bool has_corrections = false;

    // H-bond shielding coupling constant (ppm*A^3)
    // Element-dependent (N-H...O vs O-H...O have different couplings)
    array<double, 5> eta;  // index by Element of the observed atom

    // Effective distance exponent (default: 3.0, classical dipolar)
    // The calibration pipeline may learn that 1/r^2.8 or 1/r^3.2 fits
    // better than the strict 1/r^3 dipolar form.
    double distance_exponent;

    // Backbone vs sidechain H-bond weights (dimensionless)
    double w_backbone;   // default: 1.0
    double w_sidechain;  // default: 1.0
};
```

### Can T2 corrections feed back?

**Yes, through the angular function.** The T2 pattern of the H-bond
dipolar tensor is controlled by the angular distribution of shielding.
If the calibration pipeline predicts a T2 that deviates from cos^2(alpha)
angular dependence, this constrains the effective angular function.
However, changing the angular function (e.g., cos^1.8 instead of
cos^2) would require modifying the equation form, not just a parameter.

The more practical T2 feedback is through the backbone/sidechain
weight ratio (w_backbone / w_sidechain). Backbone and sidechain
H-bonds have different angular distributions (backbone H-bonds are
more linear), so their T2 patterns differ. The ratio between the
weights is constrained by the T2 data.

---

## Summary: Parameter Count and Classification

### Element-independent parameters (same for all nuclei)

| Calculator | Parameters | Count |
|------------|-----------|-------|
| Biot-Savart | I[8], d[8], stacking_scale | 17 |
| Haigh-Mallion | J[8], crossover_distance | 9 |
| McConnell | Delta_chi[6], midpoint_shift[6] | 12 |
| Pi-Quadrupole | Q[8], height_offset[8] | 16 |
| Ring Susceptibility | Delta_chi_ring[8] | 8 |
| Dispersion | r_min, r_cut | 2 |
| H-Bond | distance_exponent, w_bb, w_sc | 3 |
| **Subtotal** | | **67** |

### Element-dependent parameters (different per nucleus type)

| Calculator | Parameters | Count |
|------------|-----------|-------|
| Coulomb/Buckingham | A[5], B[5], gamma[5], charge_scale | 16 |
| Dispersion | alpha_disp[5] | 5 |
| H-Bond | eta[5] | 5 |
| **Subtotal** | | **26** |

### Grand total: 93 parameters

This is physically correct: element-independent parameters (ring
current, bond anisotropy) are NOT multiplied by element count.
Element-dependent parameters (Buckingham, dispersion coupling)
ARE per-element.

Compare with the naive approach: 93 vs the 95-parameter
element-dependent model from v1 (which put element dependence in
the wrong place by making ring current intensities element-dependent).
Similar parameter count, but physically correct decomposition.

---

## Parameter Delivery

Each calculator has a `DefaultParams()` method returning literature values.
TOML configuration can override any field. The calculator receives a
parameter struct and computes — it does not know whether the values are
literature defaults or calibrated overrides.

```cpp
class BiotSavartResult : public ConformationResult {
public:
    static BiotSavartParams DefaultParams() {
        BiotSavartParams p;
        for (int t = 0; t < 8; t++) {
            p.intensity[t] = RingTypeDefaults::Intensity(t);
            p.lobe_offset[t] = RingTypeDefaults::LobeOffset(t);
        }
        p.stacking_scale = 1.0;
        return p;
    }

    // params may come from DefaultParams() or TOML configuration
    void ComputeAllRings(ProteinConformation& conf,
                         const BiotSavartParams& params);
};
```

---

## T2 Tensor Feedback: Full Analysis

### The question

Can the calibration pipeline's T2 (L=2, anisotropy) output feed BACK into
the classical calculations, not just compare against them?

### The answer: it depends on the calculator's parameterisation

**Direct T2 feedback (T2 constrains a parameter):**

1. **McConnell**: T2 constrains Delta_chi AND midpoint_shift. The
   dipolar tensor's angular pattern is sensitive to where the effective
   dipole sits along the bond. The model's T2 prediction provides
   angular information that T0 alone cannot.

2. **Coulomb/Buckingham**: T2 constrains gamma (the EFG-to-T2 coupling).
   This is the most direct feedback: the model predicts T2, the
   calculator computes EFG_T2, and gamma = T2_predicted / EFG_T2.

3. **Pi-Quadrupole**: T2 constrains height_offset. The quadrupole
   EFG's angular pattern changes with the effective height of the
   charge distribution above the ring plane.

**Indirect T2 feedback (T2 constrains a parameter through fitting):**

4. **Biot-Savart**: T2 constrains lobe_offset d through the angular
   shape of the B-field. The calibration model learns d by fitting its T2 predictions
   to DFT T2. The constraint is mediated through the fitting process,
   not through a direct formula.

5. **Haigh-Mallion**: T2 constrains the BS/HM intensity ratio. Since
   BS and HM produce opposing T2 patterns, the calibration model's T2 prediction
   determines the optimal blend.

**No T2 feedback (parameter does not affect angular pattern):**

6. **Ring Susceptibility**: Delta_chi_ring scales magnitude uniformly.
   The angular pattern is pure (3cos^2 theta - 1)/r^3 geometry.

7. **Dispersion**: T2 component is too small to constrain parameters
   independently.

8. **H-Bond**: T2 constrains backbone/sidechain ratio weakly, but the
   signal is small compared to other calculators.

### Can T2 be injected as a TENSOR input to the calculation?

No, in the strict sense. Each calculator has a fixed functional form
(Biot-Savart law, Coulomb's law, dipolar tensor formula). These produce
tensors from geometry and parameters. Injecting an external T2 tensor
INTO the calculation would break the physical model: you would be
computing a Biot-Savart field that is not actually a Biot-Savart field.

The correct flow is:

```
C++ calculators → full tensor output (NPY)
                       |
DFT delta tensors → also exported (NPY)
                       |
                       v
Python calibration pipeline: compare calculator T2 vs DFT T2
                       |
                       v
Calibrated parameters (TOML) → C++ (next run)
```

The T2 residual (DFT delta T2 minus classical kernel T2) is computed
in the calibration pipeline. It tells the downstream e3nn model what
physics the classical calculators are NOT capturing. The calculators
produce it; they do not consume it.

### The exception: the multiplicative stacking scale

There is ONE place where a tensor-derived quantity feeds back into
a calculation: the stacking_scale parameter in BiotSavartParams.

When two rings are stacked (normal-normal dot product near 1.0,
center distance < 5A), their mutual magnetic susceptibility changes.
The calibration pipeline can learn a stacking_scale > 1.0 or < 1.0 that
modifies the effective ring current when multiple rings overlap.

This stacking_scale is a scalar (L=0) derived from the baseline
calibration model's analysis of ring-stacking T2 patterns. The calibration model notices
that stacked-ring T2 predictions are systematically different from
the sum of individual ring T2 predictions, and encodes this as
a multiplicative correction.

This is L=0 feedback (a scalar scale factor), but it was LEARNED from
T2 data. The T2 information flows through the training process into
L=0 parameters, not directly into the calculation.

---

## Calibration Pipeline (learn/c_equivariant/)

The calibration pipeline is a Python e3nn equivariant model. It reads
the NPY feature arrays that the C++ system writes (via WriteFeatures)
and trains against DFT WT-ALA delta tensors. See learn/JOURNAL.md for
training rounds and results.

The pipeline discovers optimal parameter values. After training, the
calibrated parameters are written to a TOML configuration file that the
C++ system reads on the next extraction run. No model inference happens
inside the C++ library.

### What the pipeline sees

```
Input: per-atom NPY arrays from C++ (46 arrays per conformation)
  - L=0 scalars: ring counts, distances, McConnell, etc.
  - L=2 tensors: per-calculator T2 components
  - DFT delta tensors (when available via OrcaShieldingResult)

Target: DFT WT-ALA delta T2
Measure: T2 R² across 723 proteins
```

### How parameters flow back

Calibrated values are extracted from the trained model as fixed numbers
and written to TOML. The C++ system reads them at startup. If no TOML
override exists for a parameter, the literature default is used.
