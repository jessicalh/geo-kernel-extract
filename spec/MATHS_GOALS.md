# Mathematical Goals: What This System Must Get Right

This document specifies what the NMR shielding tensor system must
compute correctly, why, and how we validate it. It references the
existing specification documents and identifies where the old system
went wrong.

This is the document for the 4-way critique: the human researcher,
the primary AI collaborator, an independent Opus review, and a Gemini
review. Each reviewer validates the maths from their perspective.

---

## What We Are Doing (and what we were doing wrong)

### The thesis hypothesis

Classical electromagnetic calculations (ring currents, bond anisotropy,
electric field gradients, etc.) can be refined by training equivariant
tensor product models against DFT shielding tensors. The decomposed
per-calculator results serve as physics-grounded features for
understanding shielding perturbations at ALL irrep levels (L=0, L=1,
L=2), not just isotropic shifts.

This is NOT an NMR chemical shift predictor. It is an exploration of
how classical physics maps onto quantum-mechanical shielding, where
it succeeds, and where it fails — diagnosed by the L=2 angular
residual.

### What the old system got wrong

The old system (biot-savart/) had TWO codebases that diverged:

1. **v1**: Good feature engineering (136 features), but old physics
   conventions, string-based atom identity, and bugs in the molecular
   graph BFS that produced wrong through-bond distances
   (UNIFICATION_BUGS.md Bug 2).

2. **v2**: Corrected physics, clean object model, but could not extract
   features because the port from v1 was incomplete and the remaining
   features used string-based patterns that broke when the object model
   changed (UNIFICATION_BUGS.md Bugs 1, 3).

3. **The histidine disaster** (Bug 1): v2 normalised HID/HIE/HIP to
   HIS, destroying protonation state and producing wrong ring type
   assignments for every histidine in every protein. All v2 ring current
   calculations on HIS-containing proteins were wrong.

4. **Three different T2 normalizations**: The design review
   (DESIGN_REVIEW_PRELIM.md) found three incompatible T2 normalization
   conventions in the codebase. The new system uses isometric
   normalization (verified in SphericalTensor::Decompose, 10 tests pass).

5. **Ring current intensities differ between modules**: irreps.py used
   Giessner-Prettre 1969 values. ring_relevance.py used fitted values.
   The C++ RingType table used DFT-calibrated values. The training data
   was produced with the Giessner-Prettre values. Changing them changes
   which atoms enter training.

The rewrite starts from first principles with one sign convention, one
normalization, one set of ring type constants, and typed objects that
cannot silently lose protonation state.

---

## The Three Mathematical Pillars

### Pillar 1: T0 parameter corrections (scalar tuning)

The ParameterCorrectionResult (e3nn model) learns corrected scalar
parameters for each classical calculator. These are L=0 outputs:
ring current intensities, bond anisotropies, Buckingham coefficients.
93 parameters total (CALCULATOR_PARAMETER_API.md).

**What must be correct:**
- Each parameter connects to a real physical quantity with known units
- Default values match published literature (with references)
- The parameter enters the calculator equation at the right place
- Changing the parameter changes the calculator output proportionally

**Validation:** For each of the 93 parameters:
- The equation in which it appears is stated with a literature reference
- The default value is stated with a literature reference
- A unit test verifies that doubling the parameter doubles the relevant
  output component (linearity test for linear parameters; correct
  functional dependence for nonlinear ones like lobe offset)

**References:**
- CALCULATOR_PARAMETER_API.md: all 93 parameters with equations
- CONSTITUTION.md: per-calculator minimum output specifications
- OBJECT_MODEL.md: ConformationAtom shielding_contribution fields

### Pillar 2: T2 angular residual (diagnostic and derivation)

After classical calculators run with corrected T0 parameters, the L=2
residual reveals where the classical angular physics breaks down. This
is the thesis's primary analytical result.

**What must be correct:**
- SphericalTensor decomposition is mathematically correct (isometric
  normalization, roundtrip verified, 10 tests)
- Each calculator produces full rank-2 tensor output (not just scalar)
- Shielding contributions are in ppm, comparable across calculators
- The residual subtraction is correct: per-atom, per-irrep, per-calculator

**T2 derivation goals:**
- McConnell midpoint_shift: T2 pattern constrains where the effective
  dipole center sits along the bond. McConnell is the dominant T2
  contributor (|alpha| = 3.61 from DESIGN_REVIEW_PRELIM.md analysis)
- Biot-Savart lobe_offset: T2 pattern constrains the pi-electron
  lobe height. Different d produces different angular dependence
- HM/BS crossover: at what distance does the surface integral (HM)
  angular pattern better match DFT than the line integral (BS)?
- Pi-quadrupole height_offset: T2 pattern constrains the effective
  height of the quadrupole source above the ring plane
- Coulomb gamma: direct T2-to-parameter channel. gamma maps EFG to
  T2 shielding. The model's T2 prediction constrains gamma.

**Per-calculator T2 explained fields:**
- Each calculator stores how much T2 residual it reduced
  (OBJECT_MODEL.md: per-calculator shielding_contribution)
- BS-HM T2 disagreement: Biot-Savart and Haigh-Mallion make opposing
  T2 predictions at the same geometry. The disagreement vector per
  atom is a diagnostic and a feature.

**Validation:**
- T2 residual map across 553 proteins, per ring type, per calculator
- T2 residual as a function of distance to nearest ring
- Does BS or HM better explain T2 near PHE rings? Near TRP?
- Where is T2 residual large far from rings? (missing physics)

**References:**
- CONSTITUTION.md: "The T2 residual: where classical angular physics
  breaks down" section
- CALCULATOR_PARAMETER_API.md: T2 feedback analysis per calculator
- OBJECT_MODEL.md: shielding_contribution fields, per-calculator T2

### Pillar 3: Classical calculation correctness

Each calculator must produce the correct answer from its equation
before any parameter correction is applied. The default-parameter
path must match published analytical results for known test geometries.

**The eight calculators and their equations:**

#### 1. Biot-Savart ring current
Johnson-Bovey double-loop model (Johnson & Bovey 1958).

```
B(r) = (mu_0 * I) / (4 * pi) * sum_loop integral[ dl x r' / |r'|^3 ]
G_ab(r) = -(1/I) * n_b * B_a(r) * PPM_FACTOR
sigma_t = I_t * G_T0_t
```

Sign convention: G_ab negative → sigma positive above ring (shielded).
Verification: proton 3A above PHE ring on normal axis has sigma > 0.

#### 2. Haigh-Mallion surface integral
Haigh & Mallion 1973/1979. 7-point Gaussian quadrature on fan
triangulation.

```
sigma_HM = I_HM * integral_S [ (1/|r-r'|^3) * (r-r') . dS ]
```

Key difference from BS: produces rank-2 tensor (not rank-1). BS and HM
make opposing T2 predictions at the same geometry.

#### 3. McConnell bond anisotropy
McConnell 1957. Dipolar kernel from bond midpoint.

```
T_ab = Delta_chi * (3 * d_a * d_b / r^5 - delta_ab / r^3)
f_McConnell = (3 * cos^2(theta) - 1) / r^3   [derived FROM tensor trace]
```

The scalar factor is DERIVED from the full tensor. Never computed
instead of the tensor.

Per-bond-category decomposition: PeptideCO, PeptideCN, BackboneOther,
SidechainCO, Aromatic, SidechainOther.

#### 4. Coulomb electric field gradient (Buckingham)
Buckingham 1960.

```
sigma_Buckingham = A_elem * E_z + B_elem * E_z^2     (T0)
T2_Buckingham = gamma_elem * EFG_T2                   (T2)
E_a = sum_j [ q_j * (r_j - r_i)_a / |r_j - r_i|^3 ]
V_ab = sum_j [ q_j * (3*(r_j-r_i)_a*(r_j-r_i)_b / |r_j-r_i|^5
                      - delta_ab / |r_j-r_i|^3) ]
```

Element-dependent coefficients. E-field from ChargeAssignmentResult
(ff14SB) or MopacResult (PM7 Mulliken). E-field decomposed: backbone, sidechain, aromatic,
solvent (APBS - vacuum).

#### 5. Pi-quadrupole EFG
Buckingham 1959. Axial quadrupole from pi-electron density.

```
V_quad = Q * (3*cos^2(theta) - 1) / r^4
V_ab = Q * [ 15*dn^2*d_a*d_b/r^7 - 3*dn*(n_a*d_b + n_b*d_a)/r^5
           - 3*d_a*d_b/r^5 + delta_ab*(3*dn^2/r^5 - 1/r^3) ]
```

#### 6. Ring susceptibility anisotropy
McConnell 1957 applied to whole-ring diamagnetic susceptibility.

```
T_ab = Delta_chi_ring * (3*d_a*d_b/r^5 - delta_ab/r^3)
```

Same dipolar form as bond McConnell but from ring center.

#### 7. London dispersion
Per-vertex van der Waals. London 1937.

```
T_disp_ab = sum_vertices [ C6 * (3*d_a*d_b/r^8 - delta_ab/r^6) ]
disp_scalar = sum_vertices [ C6 / r^6 ]
```

#### 8. H-bond dipolar
Dipolar tensor to H-bond partner. Barfield & Karplus 1969, Cornilescu
& Bax 1999.

```
T_ab = eta * (3*d_a*d_b/r^5 - delta_ab/r^3)
sigma_hbond = eta * cos^2(alpha) / r^3
```

**What must be correct for each calculator:**
- The equation matches the literature reference
- The sign convention matches (CONSTITUTION.md: NMR convention,
  shielding positive)
- Full rank-2 tensor output (Mat3 + SphericalTensor)
- Default parameters from published sources
- Unit test with known analytical geometry

**References:**
- CONSTITUTION.md: per-calculator minimum output, sign convention
- CALCULATOR_PARAMETER_API.md: equations, default values, references
- EXTRACTION_ORDER.md: what each result stores on ConformationAtom
- advisor_precis.tex: calculator table with references

---

## What the External Tools Must Provide

The classical calculators need correct inputs from the external tools.
If the inputs are wrong, the calculators are precisely wrong.

### APBS (Poisson-Boltzmann solvation)

**Provides:** Per-atom solvated E-field (Vec3) and EFG (Mat3 + SphericalTensor)
**Used by:** CoulombResult (Buckingham E-field effect), all solvent features
**Must be correct:**
- E-field units: V/Angstrom
- EFG is traceless (Laplacian of potential is zero in charge-free region)
- Solvation effect is real: E_solvated differs from E_vacuum by the
  dielectric screening factor (~20x near the surface)
- Grid resolution sufficient (current: 161^3, ~0.625A spacing)
- NaN/Inf sanitised

**Validation:**
- Compare APBS E-field at a known atom to the old system's output
- Verify tracelessness of EFG (trace < epsilon)
- Verify E_solvated != E_vacuum (solvation effect present)
- Spot-check against published APBS results for a standard protein

### MOPAC (PM7+MOZYME conformation electronic structure)

**Provides:** Per-atom Mulliken charges, s/p orbital populations,
per-bond Wiberg bond orders, heat of formation
**Used by:** MopacCoulombResult (QM-charge EFG kernel),
MopacMcConnellResult (bond-order-weighted anisotropy kernel),
ParameterCorrectionResult (input feature), feature extraction,
MutationDeltaResult (delta_mopac_charge)
**Must be correct:**
- Charges reasonable (|q| < 3.0 per atom)
- Wiberg bond orders: ~1 for single bonds, ~2 for double, ~1.5 for aromatic
- Orbital populations: s_pop > 0 for all atoms, p_pop > 0 for heavy atoms
- Heat of formation negative for proteins

**Validation:**
- Charge range check on 1UBQ (1231 atoms, ~80s)
- Bond order sanity check on known bonds
- Compare MOPAC charges to ff14SB charges — they should correlate but differ
  (MOPAC captures polarisation, ff14SB does not)
- TopologyBondOrder(bi) agrees with BondOrder(atom_a, atom_b) for all bonds

### DSSP (secondary structure)

**Provides:** Per-residue SS assignment, phi/psi angles, SASA, H-bond partners
**Used by:** Feature extraction (structural context), ParameterCorrectionResult
**Must be correct:**
- SS assignments match published 1UBQ secondary structure
- Phi/psi in radians, range [-pi, pi]
- SASA in Angstroms^2, non-negative
- H-bond partners are valid residue indices

**Validation:**
- Compare to RCSB PDB secondary structure annotation for 1UBQ
- Verify known helix (residues 23-34) and sheet regions

### ff14SB charges (tleap/prmtop)

**Provides:** Per-atom partial charges and VdW radii
**Used by:** CoulombResult (Coulomb E-field), APBS (PQR input)
**Must be correct:**
- Charges from the ff14SB force field, not made up
- Total charge is integer (net formal charge of protein)
- Backbone atoms have expected charge ranges
- Charges depend on protonation state (HID vs HIE have different charges)

**Validation:**
- Total charge check
- Known backbone N/CA/C/O charge ranges
- Compare to published ff14SB charge tables
- Verify protonation-dependent charges differ

---

## The Environment Floor (what classical calculators must beat)

The ParameterCorrectionResult is trained on APBS + MOPAC + DSSP + DFT.
Its prediction already incorporates environment tensor data at all
irrep levels:

- APBS EFG: L=2 angular structure of solvation environment
- MOPAC charges and bond orders: L=0 conformation electronic structure
- DFT shielding: L=0+L=1+L=2 full tensor baseline

A classical calculator earns its place by REDUCING the residual beyond
what these environment tools already capture. If APBS EFG T2 already
explains the angular structure, adding Coulomb EFG doesn't help.

**The test for each calculator:**
1. Residual before attaching = model prediction from environment only
2. Attach calculator with default parameters
3. Residual after attaching = prediction minus calculator contribution
4. If residual decreased: calculator adds physics
5. If residual unchanged: calculator is redundant with environment data
6. Measure at L=0 AND L=2 independently

This is CONSTITUTION.md "The bar for classical calculators" and
"Environment tools set the floor."

---

## Mutant Isolation: The Training Signal

The training signal is NOT "predict the DFT shielding." It is
"predict what changes when you mutate."

- WT-ALA: ring removal → isolates ring current contribution
- ASP→ASN: charge flip → isolates Buckingham contribution
- TYR→PHE: ring modification → isolates intensity from geometry
- LYS→MET: salt-bridge break → isolates solvation contribution

Each mutant type gives the model a controlled experiment. Classical
features extracted on BOTH WT and mutant structures, differenced, and
compared to the DFT delta.

**What must be correct:**
- Atom matching between WT and mutant is position-based, not name-based
- Delta computation: WT shielding minus mutant shielding, per atom
- Delta decomposition via SphericalTensor (same normalization as features)
- The right DFT level on both sides (r2SCAN for mutant pairs, same
  basis set, delta cancels basis set errors)

**References:**
- CONSTITUTION.md: "Mutant Comparison: First-Class Operation"
- CONSTITUTION.md: MutantProteinConformationComparison sealed ABC
- REWRITE_DECISIONS.md: mutant strategy (3 new categories planned)

---

## Summary: What We Validate and When

### Before any calculator (Layer 0):
- [x] SphericalTensor roundtrip (DONE, 10 tests)
- [x] PDB loading with correct residue/ring/bond types (DONE)
- [x] DSSP assignments match known 1UBQ structure (DONE, 6 tests)
- [x] APBS E-field is real PB solve, not vacuum (DONE, 3 tests)
- [ ] APBS EFG is traceless (verify)
- [x] ff14SB charges: total charge is integer on protonated structures
      (DONE — PrmtopChargeSource on ORCA prmtop gives -4.0 exact;
      GmxTprChargeSource on CHARMM .tpr gives -2.0 exact)
- [x] MOPAC charges computed on 1UBQ (DONE — 1231 atoms, charges
      reasonable, heat of formation negative, bond orders present)
- [x] MOPAC bond orders reasonable (DONE — covalent bonds > 0.1,
      topology bridge verified, sorted descending per atom)
- [x] Protonation detection matches hydrogen atoms (DONE, 11 tests —
      works on PDB-standard and CHARMM naming, PROPKA agreement tested)
- [x] ORCA shielding tensors load correctly (DONE — dia+para+total,
      543 atoms element-verified, T0 matches raw output)
- [x] ORCA mutant comparison (DONE — MutationDeltaResult with rich
      atom matching, delta tensors T0+T1+T2, ring proximity with
      cylindrical coords, distance decay curve. 465 pairs batch tested.)
- [x] MD trajectory loading (DONE — GromacsEnsembleLoader reads TPR
      via libgromacs, 10 MDFrameConformations per protein, full
      pipeline per frame, 283 tests pass)

### For each classical calculator (Layer 1+):
All 8 calculators pass items 1-6. Item 7 awaits ParameterCorrectionResult.
- [x] Equation matches literature reference
- [x] Sign convention matches (proton above PHE → sigma > 0)
- [x] Full rank-2 tensor output (T2 non-zero near relevant sources)
- [x] Default parameters from published sources (geometric kernels only)
- [x] Analytical test geometry produces correct value
- [x] Shielding contribution stored on ConformationAtom
- [ ] Residual decreases when calculator attaches (needs ParameterCorrectionResult)

### McConnell (first calculator, 2026-04-02):
- [x] Equation: full asymmetric tensor M from first-principles derivation
      (see GEOMETRIC_KERNEL_CATALOGUE.md). NOT the incomplete Δχ×K formula.
- [ ] Sign convention (not yet verified — needs parameterised output)
- [x] Full rank-2 tensor output: T0+T1+T2 all non-zero. T2 > T0 as expected.
- [ ] Default parameters (not yet applied — geometric kernel only)
- [x] Analytical test geometry: atom at (3,0,0), bond at origin along z.
      f = -1/27, K traceless, T0 = f verified at machine precision.
- [x] Shielding contribution stored: mc_shielding_contribution on ConformationAtom
- [ ] Residual decrease (needs ParameterCorrectionResult)
- [x] Batch validated: 465 proteins, 288,741 matched atoms, 0 failures
- [x] KernelFilterSet: SelfSourceFilter + DipolarNearFieldFilter(bond_length)

### CoulombResult (2026-04-02):
- [x] Equation: E = ke * sum q_j r/r³, EFG = ke * sum q_j K_ab. Units V/A.
- [x] Full rank-2 EFG tensor: traceless (Gauss's law) to 1.4e-14.
- [x] Decomposition: backbone + sidechain + aromatic = total (exact to 3.2e-14)
- [x] E-field magnitudes: grand mean 2.81 V/A (consistent with Case 1995)
- [x] All 465 proteins: integer total charge from prmtop, 0 failures
- [x] Solvent contribution: APBS (V/A) - vacuum (V/A), both same units
- [x] KernelFilterSet: SelfSourceFilter (i≠j)
- [x] T2 independence: |cos| = 0.47 vs McConnell, 0.40 vs RingChi
- [ ] Buckingham shielding contribution (T0 from E_z needs parameters)

### RingSusceptibilityResult (2026-04-02):
- [x] Equation: full McConnell tensor with b_hat → ring normal n_hat
- [x] T0 = f identity: verified at 3.1e-15 across 914,155 ring-atom pairs
- [x] Full rank-2 tensor: T0+T1+T2 non-zero. T2/T0 ratio ~2.4 as expected.
- [x] Per-ring-type validated: PHE/TYR/TRP6 consistent, TRP5/HIE 30% higher,
      TRP9 2x (fused ring). Physics-consistent size dependence.
- [x] DFT proximity: 6.9x near/far ratio at mutation sites
- [x] KernelFilterSet: DipolarNearFieldFilter(ring_diameter)
- [x] T2 independence: |cos| = 0.41 vs McConnell, 0.40 vs Coulomb
- [x] Batch validated: 465 proteins, 0 failures

### HBondResult (2026-04-02):
- [x] Equation: full McConnell tensor with b_hat → H-bond direction h_hat
- [x] H-bonds from DSSP Kabsch-Sander criterion, resolved to atom positions
- [x] Full rank-2 tensor: T0+T1+T2 non-zero. T2/T0 ratio ~2.4.
- [x] Per-SS validated: helix 6809, sheet 3491, coil 8168 donors
- [x] KernelFilterSet: SelfSourceFilter + DipolarNearFieldFilter(N...O distance)
      Effect: max |T2| 1908 → 0.78 A^-3 (near-field atoms excluded)
- [x] T2 independence: |cos| = 0.38 vs McConnell, 0.40 vs Coulomb, 0.38 vs RingChi
- [x] Batch validated: 466 proteins, 0 failures

### BiotSavartResult (2026-04-02):
- [x] Equation: Johnson-Bovey double-loop, B-field in SI, G = -n⊗B × PPM_FACTOR
- [x] Sign convention: I=-12, atom 3A above PHE → sigma = +1.40 ppm (shielded)
- [x] Full rank-1 tensor: G_ab = n_b * B_a, stored as full Mat3 + SphericalTensor
- [x] Per-ring attribution with cylindrical B-field components (B_n, B_rho, B_phi)
- [x] KernelFilterSet: DipolarNearFieldFilter + RingBondedExclusionFilter
- [x] T2 independence: |cos| = 0.37-0.38 vs McConnell/Coulomb/HBond
- [x] BS-HM T2 redundancy: cos = 0.999 (same angular structure)
- [x] Fused ring additivity: T0(TRP5)+T0(TRP6)/T0(TRP9) = 1.000
- [x] Batch validated: 720 proteins, 1.03M atom-ring pairs, 0 failures

### HaighMallionResult (2026-04-02):
- [x] Equation: surface integral with 7-point Gaussian quadrature (Stroud T2:5-1)
- [x] Adaptive subdivision at 2.0A and 1.0A for near-field atoms
- [x] TWO tensors: H (raw integral, symmetric traceless A⁻¹) + G (rank-1 kernel)
- [x] Boyd-Skrynnikov construction: G = -n⊗(H·n)
- [x] H traceless to 2.8e-16 across all ring types
- [x] KernelFilterSet: DipolarNearFieldFilter + RingBondedExclusionFilter
- [x] Fused ring: ratio 1.000 after RingBondedExclusionFilter (was 1.127)
- [x] Batch validated: 720 proteins, 0 failures

### PiQuadrupoleResult (2026-04-02):
- [x] Equation: Stone T-tensor EFG from point axial quadrupole (4th derivative of 1/r)
- [x] Traceless at 5e-15 across 1.03M atom-ring pairs (pure T2, T0=0)
- [x] 1/r^5 leading decay, finite-difference verified
- [x] KernelFilterSet: DipolarNearFieldFilter + RingBondedExclusionFilter
      (max |T2| reduced 11x: 7.39→0.66 A^-5 after audit)
- [x] T2 independence: |cos| = 0.38 vs McConnell, 0.40 vs Coulomb
- [x] DFT proximity: 29.0x near/far ratio at mutation sites
- [x] Batch validated: 723 proteins, 0 failures

### DispersionResult (2026-04-02):
- [x] Equation: London 1/r^6 kernel summed over ring vertices
- [x] CHARMM switching function (C1 continuous, R_switch=4.3A, R_cut=5.0A)
- [x] Through-bond vertex exclusion + DipolarNearFieldFilter
- [x] DFT proximity: 204x near/far ratio (strongest proximity signal)
- [x] T2 independence verified (small T2 contribution, dominated by others)
- [x] Batch validated: 723 proteins, 124K contact pairs, 0 failures

### For the ParameterCorrectionResult (not yet implemented):
- [ ] Environment features correctly extracted
- [ ] Model predicts full tensor (L=0+L=1+L=2)
- [ ] 93 parameter corrections are physically reasonable
- [ ] Residual tracking per irrep works
- [ ] T2 residual analysis produces interpretable patterns

### For the complete system:
- [ ] WT-ALA delta features match DFT delta within model error
- [ ] Sign accuracy ≥ 89% at 95% coverage (old best: 89%)
- [ ] T0 R² ≥ 0.875 (old bilinear best)
- [ ] T2 residual map reveals where classical physics breaks down
- [ ] Each calculator demonstrably reduces the residual
