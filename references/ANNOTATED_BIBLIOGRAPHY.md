# Annotated Bibliography

Compiled for the NMR shielding tensor calibration thesis. References
marked **[CORE]** are the 12-15 that directly underpin the work.
Remaining entries are the reading pool (~90 candidates).

---

## A. Classical NMR Shielding Theory

The physical models implemented as geometric kernels in this system.
These papers contain the equations in `spec/MATHS_GOALS.md`.

### Ring Current Models

**[CORE] A1. Johnson, C.E. & Bovey, F.A. (1958)** "Calculation of
nuclear magnetic resonance spectra of aromatic hydrocarbons." *J. Chem.
Phys.* 29, 1012--1014. DOI: 10.1063/1.1744645

The double-loop Biot-Savart model: two current loops at +/-d from
the aromatic plane carry pi-electron current, and the B-field at
a probe nucleus is computed via the Biot-Savart law. This is the
physical model inside `BiotSavartResult.cpp`. Simple, analytical
where possible, and the standard ring-current model for small
molecules. Its tensor structure (rank-1 B-field, rank-2 kernel via
outer product with ring normal) is well-defined.

**[CORE] A2. Haigh, C.W. & Mallion, R.B. (1979)** "Ring current
theories in nuclear magnetic resonance." *Prog. Nucl. Magn. Reson.
Spectrosc.* 13, 303--344.

Comprehensive review and definitive formulation of the surface
integral approach. Instead of wire loops, integrates the dipolar
kernel over the ring surface, producing a rank-2 tensor directly.
The 7-point Gaussian quadrature (Stroud T2:5-1) used in
`HaighMallionResult.cpp` comes from the numerical treatment here.
Key theoretical point: BS and HM make opposing T2 predictions at
the same geometry (MATHS_GOALS.md line 127), which is a diagnostic
the DFT calibration can resolve.

A3. Haigh, C.W. & Mallion, R.B. (1972) "New tables of 'ring current'
shielding in proton magnetic resonance." *Org. Magn. Reson.* 4,
203--228.

The original numerical tables. Superseded by the 1979 review but
historically important as the first practical implementation.

A4. Giessner-Prettre, C. & Pullman, B. (1969) "On the ring current
effect in nucleic acids." *C. R. Acad. Sci. Ser. D* 268, 1115--1117.

Ring current intensity factors for nucleic acid bases computed from
SCF molecular orbital theory. The ring current intensities in the
codebase (PHE: -12.0 nA, TYR: -11.28 nA, etc.) derive from this
lineage. `MATHS_GOALS.md` line 59 notes these as the baseline
values against which DFT-calibrated values are compared.

A5. Giessner-Prettre, C. & Pullman, B. (1987) "Intermolecular nuclear
shielding from the ring current in nucleic acids and aromatic amino
acids." *Q. Rev. Biophys.* 20, 113--172. DOI: 10.1017/S0033583500004169

Extended treatment with intensity factors for all standard aromatic
amino acids. The later values here refine the 1969 estimates.

A6. Pople, J.A. (1956) "Proton magnetic resonance of hydrocarbons."
*J. Chem. Phys.* 24, 1111. DOI: 10.1063/1.1742701

The original "ring current" concept: pi-electron circulation in
aromatic rings produces a magnetic dipole that shifts nearby proton
resonances. Pople's point-dipole approximation is the simplest
model (1/r^3 axial dependence). Too crude for quantitative work
but establishes the physical picture every subsequent model refines.

A7. Waugh, J.S. & Fessenden, R.W. (1957) "Nuclear resonance spectra
of hydrocarbons: the free electron model." *J. Am. Chem. Soc.* 79,
846--849.

Free-electron model for ring current intensities. Historical
precursor to the Giessner-Prettre SCF calculations. Sets the scale
of ring current effects (~1-3 ppm for protons near aromatic rings).

**A8. Boyd, J. & Skrynnikov, N.R. (2002)** "Calculation of 1H chemical
shift anisotropies in proteins." *J. Am. Chem. Soc.* 124, 1832--1833.
DOI: 10.1021/ja017638f

Full shielding tensor from ring magnetization: sigma_ab = -n_b(H.n)_a.
Cited in GEOMETRIC_KERNEL_CATALOGUE.md line 218. Critical for
connecting the Haigh-Mallion surface integral to a measurable CSA
tensor. This paper bridges the isotropic ring current shift
(everyone computes) to the full rank-2 tensor (which your T2
analysis requires).

### Bond Magnetic Anisotropy

**[CORE] A9. McConnell, H.M. (1957)** "Theory of nuclear magnetic
shielding in molecules. I. Long-range dipolar shielding of protons."
*J. Chem. Phys.* 27, 226--229. DOI: 10.1063/1.1743676

The original derivation of the bond magnetic anisotropy contribution
to NMR shielding. The dipolar kernel T_ab = Delta_chi * (3d_a*d_b/r^5
- delta_ab/r^3) implemented in `McConnellResult.cpp` is equation (7)
of this paper. Also applied to whole-ring susceptibility in
`RingSusceptibilityResult.cpp`. The McConnell factor
(3cos^2(theta)-1)/r^3 is the scalar trace, but the full tensor is
what the code computes.

A10. Flygare, W.H. (1974) "Magnetic interactions in molecules and an
analysis of molecular electronic charge distribution from magnetic
parameters." *Chem. Rev.* 74, 653--687. DOI: 10.1021/cr60292a003

Comprehensive treatment of molecular magnetic susceptibility tensors
and their relationship to electronic structure. Provides the
theoretical basis for bond susceptibility anisotropy parameters
(Delta_chi values) used in the McConnell kernel.

### Electric Field Effects

**[CORE] A11. Buckingham, A.D. (1960)** "Chemical shifts in the
nuclear magnetic resonance spectra of molecules containing polar
groups." *Can. J. Chem.* 38, 300--307. DOI: 10.1139/v60-040

Establishes that nearby electric fields perturb shielding linearly
(sigma = A*E_z) and that the electric field gradient produces a
rank-2 tensor contribution (sigma_ab = gamma*V_ab). This is the
theoretical foundation for `CoulombResult.cpp` and the Buckingham
coefficients in `CALCULATOR_PARAMETER_API.md`. The E-field/EFG
decomposition maps directly to T0/T2 decomposition in the code.

A12. Buckingham, A.D. (1959) "Solvent effects in nuclear magnetic
resonance." *J. Chem. Phys.* 30, 1580. DOI: 10.1063/1.1730242

Electric quadrupole contribution from pi-electrons. The
`PiQuadrupoleResult` calculator implements this model:
V_quad = Q*(3cos^2(theta)-1)/r^4. Dropped from the calibration
pipeline (r < 0.03 with DFT) but retained in the extraction code.

### Dispersion and Hydrogen Bonding

A13. London, F. (1937) "The general theory of molecular forces."
*Trans. Faraday Soc.* 33, 8--26. DOI: 10.1039/TF9373300008

The London dispersion force: C6/r^6 interaction between
fluctuating dipoles. The `DispersionResult` kernel tensorises
this as C6*(3d_a*d_b/r^8 - delta_ab/r^6) with a CHARMM switching
function from Brooks et al. (1983).

A14. Barfield, M. & Karplus, M. (1969) "Contact contributions to
nuclear magnetic resonance chemical shifts." *J. Am. Chem. Soc.* 91,
1--16.

Hydrogen bond contributions to chemical shifts via through-space
dipolar interactions. The `HBondResult.cpp` calculator uses the
same dipolar kernel form as McConnell but with the H-bond vector.
Cited in `MATHS_GOALS.md` line 231.

A15. Cornilescu, G., Hu, J.-S. & Bax, A. (1999) "Identification of
the hydrogen bonding network in a protein by scalar couplings."
*J. Am. Chem. Soc.* 121, 2949--2950. DOI: 10.1021/ja9902221

Experimental validation of H-bond geometry effects on NMR
observables. The cos^2(alpha)/r^3 angular dependence of H-bond
shielding used in the code is grounded in this work.

---

## B. DFT Shielding Calculations

The quantum-mechanical reference that the geometric kernels are
calibrated against. The 720 WT-ALA mutant pairs use ORCA DFT.

**[CORE] B1. de Dios, A.C., Pearson, J.G. & Oldfield, E. (1993)**
"Secondary and tertiary structural effects on protein NMR chemical
shifts: an ab initio approach." *Science* 260, 1491--1496. DOI:
10.1126/science.8502992

Landmark: first ab initio prediction of 1H, 13C, 15N, 19F shifts
in proteins. Shows that phi/psi torsions dominate 13C shielding,
hydrogen bonding dominates 1HN, and long-range electrostatic fields
dominate 19F. Establishes the physical decomposition of shielding
contributions that the geometric kernels recapitulate. Your
calibration is essentially measuring these same contributions but
with a controlled mutation experiment instead of perturbation theory.

B2. Oldfield, E. (1995) "Chemical shifts and three-dimensional
protein structures." *J. Biomol. NMR* 5, 217--225. DOI:
10.1007/BF00211749

Review of early ab initio protein shift calculations. Argues that
shifts encode structural information beyond NOEs. Complements
de Dios 1993 with broader perspective.

B3. Wolinski, K., Hinton, J.F. & Pulay, P. (1990) "Efficient
implementation of the gauge-including atomic orbital method for NMR
chemical shift calculations." *J. Am. Chem. Soc.* 112, 8251--8260.
DOI: 10.1021/ja00179a005

The GIAO (Gauge-Including Atomic Orbital) method that makes DFT
shielding calculations practical. Your ORCA reference calculations
use GIAO. Understanding the method is necessary for understanding
the systematic errors in the DFT reference.

B4. Cheeseman, J.R., Trucks, G.W., Keith, T.A. & Frisch, M.J.
(1996) "A comparison of models for calculating nuclear magnetic
resonance shielding tensors." *J. Chem. Phys.* 104, 5497--5509.
DOI: 10.1063/1.471789

Systematic benchmark of HF, MP2, and DFT methods for shielding
tensors. Establishes that DFT with GIAO is accurate enough for
protein-scale calculations. Sets expectations for the accuracy of
your DFT reference.

B5. He, X., Wang, B. & Merz, K.M. Jr. (2009) "Protein NMR chemical
shift calculations based on the automated fragmentation QM/MM
approach." *J. Phys. Chem. B* 113, 10380--10388. DOI:
10.1021/jp901992p

AF-QM/MM: each residue treated at QM level, rest as point charges.
Represents the computational cost (~hours per protein) that your
classical kernel approach aims to replace with seconds-scale
evaluation.

B6. Mulder, F.A.A. & Filatov, M. (2010) "NMR chemical shift data
and ab initio shielding calculations: emerging tools for protein
structure determination." *Chem. Soc. Rev.* 39, 578--590. DOI:
10.1039/B811366C

Review bridging ab initio shielding theory and protein NMR
applications. Good overview of the accuracy-cost tradeoff that
motivates your hybrid classical-kernel + DFT-calibration approach.

### Alternative gauge methods (GIAO's parallel tradition)

B7. Keith, T.A. & Bader, R.F.W. (1993) "Calculation of magnetic
response properties using a continuous set of gauge
transformations." *Chem. Phys. Lett.* 210, 223--231. DOI:
10.1016/0009-2614(93)89127-4

CSGT: position-dependent gauge transformation, alternative to
GIAO. ORCA's documentation recommends GIAO as the sole origin for
shielding calculations, but CSGT is competitive and occasionally
correlates better with experiment. Context for understanding what
r²SCAN-via-ORCA is actually computing. Paywalled.

B8. Kutzelnigg, W. (1980) "Theory of magnetic susceptibilities
and NMR chemical shifts in terms of localized quantities." *Isr.
J. Chem.* 19, 193--200.

Origin paper for IGLO (Individual Gauge for Localized Orbitals).
Library needed.

B9. Hansen, A.E. & Bouman, T.D. (1985) "Localized orbital/local
origin method for calculation and analysis of NMR chemical
shieldings and magnetic susceptibilities." *J. Chem. Phys.* 82,
5035--5047. DOI: 10.1063/1.448625

LORG: localized orbitals with per-orbital local origins. Third
alternative gauge method. Paywalled.

**B10. Helgaker, T., Jaszuński, M. & Ruud, K. (1999)** "Ab
initio methods for the calculation of NMR shielding and indirect
spin-spin coupling constants." *Chem. Rev.* 99, 293--352. DOI:
10.1021/cr960017t

Canonical review of magnetic property calculation methods: GIAO,
CSGT, IGLO, LORG, method-basis benchmarks. The reference for
understanding what r²SCAN-via-ORCA is actually computing.
Surfaced in `spec/PHYSICS_FOUNDATIONS.md` Session 0.2. Paywalled;
library needed.

B11. Furness, J.W., Kaplan, A.D., Ning, J., Perdew, J.P. & Sun,
J. (2020) "Accurate and numerically efficient r²SCAN
meta-generalized gradient approximation." *J. Phys. Chem. Lett.*
11, 8208--8215. DOI: 10.1021/acs.jpclett.0c02405

The meta-GGA functional behind the 260-DFT calibration.
Background for understanding shielding-specific biases of r²SCAN
vs. B3LYP/PBE0. Paywalled.

### Residue-specific DFT CSA surveys (Oldfield lineage, Thread 2)

**[CORE] B12. Havlin, R.H., Le, H., Laws, D.D., de Dios, A.C. &
Oldfield, E. (1997)** "An ab initio quantum chemical
investigation of carbon-13 NMR shielding tensors in glycine,
alanine, valine, isoleucine, serine, and threonine: Comparisons
between helical and sheet tensors, and the effects of χ1 on
shielding." *J. Am. Chem. Soc.* 119, 11951--11958. DOI:
10.1021/ja971796d

Foundational DFT survey of Cα tensor principal values across
secondary structures and χ1 rotamers for six residues.
Establishes that Cα ¹³C CSA elements are individually sensitive
to (φ, ψ, χ₁) — not just the isotropic shift. The physics our
per-atom-type stratified calibration implicitly depends on is
mapped here. Paywalled.

B13. Havlin, R.H., Laws, D.D., Bitter, H.-M.L., Sanders, L.K.,
Sun, H., Grimley, J.S., Wemmer, D.E., Pines, A. & Oldfield, E.
(2001) "Carbon-13 NMR shielding in the twenty common amino acids:
Comparisons with experimental results in proteins." *J. Am. Chem.
Soc.* 123, 10362--10369. DOI: 10.1021/ja011863a

Extends B12 to all 20 amino acids. Reports ~4-5 ppm isotropic
shift between β-sheet and α-helix residues. Direct quantitative
target for our per-stratum predictions. Paywalled.

### Magnetic response beyond shielding: magnetizabilities, NICS, aromaticity

Surfaced in `spec/PHYSICS_FOUNDATIONS.md` Session 0.2. The same
GIAO / second-order response-theory machinery that produces our
shielding ground truth also produces magnetizabilities,
nucleus-independent chemical shifts (NICS), and magnetic-exaltation
aromaticity indices. Same physics, different probe points. A
reviewer in the magnetizability / aromaticity-criterion community
would expect us to know this literature exists.

B14. Schleyer, P.v.R., Maerker, C., Dransfeld, A., Jiao, H. & van
Eikema Hommes, N.J.R. (1996) "Nucleus-Independent Chemical Shifts:
A Simple and Efficient Aromaticity Probe." *J. Am. Chem. Soc.* 118,
6317--6318. DOI: 10.1021/ja960582d

NICS origin paper: evaluate the shielding at a dummy nucleus placed
at a chosen probe point (classically the ring centroid, later
NICS(1) at +1 Å above the plane). Same GIAO machinery as our DFT
ground truth, different probe location. Our calibrated ring-current
kernels implicitly measure the same induced circulation NICS
indexes; the difference is that NICS evaluates at one point while
our kernels evaluate at every atom across the protein. Paywalled.

**B15. Gershoni-Poranne, R. & Stanger, A. (2015)** "Magnetic
criteria of aromaticity." *Chem. Soc. Rev.* 44, 6597--6615. DOI:
10.1039/c5cs00114e (verify)

Critical review of magnetic-based aromaticity indicators:
susceptibility exaltation, NICS variants (NICS(0), NICS(1)zz,
NICS-scan), ring-current density maps, anisotropy of induced
current density. The closest literature sibling to what our
calibrated BS / HM intensities capture per ring type. Cite when we
position the calibration as measuring ring-susceptibility
contributions. Paywalled.

B16. Kahn, K. et al. (2022) "Molecular Magnetizabilities Computed
Via Finite Fields: Assessing Alternatives to MP2 and Revisiting
Magnetic Exaltations in Aromatic and Antiaromatic Species." Open
access (PMC8903098). DOI and full author list to verify against PDF.

Modern benchmark of ab initio magnetizability methods; finds
κ-OOMP2 best among MP2 alternatives; revisits magnetic-exaltation
aromaticity. Companion to B15 on the methods side — where B15
reviews criteria, B16 assesses which ab initio method delivers them
reliably. PDF fetched during Session 0.2.

B17. Helgaker, T., Coriani, S., Jørgensen, P., Kristensen, K.,
Olsen, J. & Ruud, K. (2012) "Recent Advances in Wave Function-Based
Methods of Molecular-Property Calculations." *Chem. Rev.* 112,
543--631. DOI: 10.1021/cr2002239 (verify)

Extended companion to B10. Where B10 (1999) is the canonical
shielding-methods review, B17 covers the broader magnetic and
electric response-property landscape 13 years later: state of the
field, CCSD(T) vs DFT trade-offs, relativistic corrections, and
the finite-field / sum-over-states alternatives to GIAO. Paywalled;
library needed.

B18. Flygare, W.H. & Benson, R.C. (1971) "The molecular Zeeman
effect and magnetic susceptibilities of small organic molecules."
*Mol. Phys.* 20, 225--250. DOI to verify.

Pre-DFT experimental anchor: magnetizability anisotropies measured
directly via molecular Zeeman spectroscopy on small organic
molecules. Our classical per-atom magnetic-anisotropy kernels live
between the small-molecule scale B18 measures and the bulk-protein
scale Babaei et al. 2017 (N3) computes; B18 is the experimental
floor of that hierarchy. Paywalled; library needed.

### r²SCAN functional: benchmarks for our DFT reference

Extends B11 Furness 2020 (r²SCAN functional paper) and B17
Helgaker-Coriani 2012 (methods review) with recent r²SCAN-specific
NMR shielding benchmarks. Our 260-DFT calibration batch was run
with r²SCAN-via-ORCA; these entries ground the systematic-bias
question at the functional-specific level.

B19. (2022) "Optimization of the r²SCAN-3c Composite
Electronic-Structure Method for Use with Slater-Type Orbital
Basis Sets." *J. Phys. Chem. A*. DOI: 10.1021/acs.jpca.2c02951.
Open access (PMC9255700). (Author list to verify — Grimme group
likely.)

r²SCAN-3c composite method: r²SCAN + tailored mTZ2P STO basis +
D4 dispersion + gCP geometric counterpoise + SR-ZORA relativity.
Benchmarked on GMTKN55 with accuracy on par with or better than
M06-2X-D3(0)/TZP at a fraction of the cost. **Relevance:** a
plausible fast composite method for future expanded-DFT coverage
during the 2026-06 μs-harvester phase where per-pose DFT budget
is the binding constraint (see `spec/MICROSECOND_MD_HARVESTER_
2026-04-22.md`). Related-family member to the r²SCAN on our
current calibration set, not the same functional.

B20. (2025) "Accurate predictions of Chemical Shifts with the
rSCAN and r²SCAN mGGA exchange-correlation functionals."
*Faraday Discussions*. DOI: 10.1039/d4fd00142g. (Author list to
verify.)

Direct r²SCAN-for-shielding benchmark on halide and oxide
inorganic compounds: r²SCAN meta-GGA delivers significantly
improved NMR shielding accuracy over GGA at marginal
computational-cost increase. **Tethered to MATHS_GOALS pillar 3
(classical calculation correctness vs DFT reference):** evidence
that our choice of r²SCAN as calibration functional sits in a
validated performance regime for shielding, even though the
benchmark system differs from protein backbones. The quality of
our 260-DFT reference is downstream of this performance floor.

---

## C. Semi-Empirical Shift Predictors

The existing methods your calibrated kernels replace or improve upon.
These use the same physical decomposition but with fixed parameters.

**[CORE] C1. Case, D.A. (1995)** "Calibration of ring-current effects
in proteins and nucleic acids." *J. Biomol. NMR* 6, 341--346. DOI:
10.1007/BF00197633

The paper most directly comparable to your thesis. DFT shielding
calculations for methane near aromatic rings, compared against
Johnson-Bovey and Haigh-Mallion models. Establishes calibrated ring
current intensity factors. Your 720-mutant DFT calibration is a
modern, tensor-valued, environment-aware extension of exactly this
experiment. Case calibrates isotropic shifts against a single-molecule
DFT; you calibrate full T2 tensors against 720 protein-embedded DFT
calculations.

C2. Osapay, K. & Case, D.A. (1991) "A new analysis of proton
chemical shifts in proteins." *J. Am. Chem. Soc.* 113, 9436--9444.

Adds peptide magnetic anisotropy and electrostatic contributions to
ring-current models, improving R from 0.70 to 0.88 for 5678 protons
across 17 proteins. Demonstrates the additive model of shielding
contributions. Your kernel architecture is the tensorial version of
this additive decomposition.

C3. Xu, X.P. & Case, D.A. (2001) "Automated prediction of 15N,
13Calpha, 13Cbeta and 13C' chemical shifts in proteins using a
density functional database." *J. Biomol. NMR* 21, 321--333. DOI:
10.1023/A:1013324104681

SHIFTS 4.0: DFT-computed database of 1335 peptides for backbone
shift prediction. R = 0.90-0.99. The DFT-database hybrid approach:
your calibration pipeline is philosophically similar but uses
mutation deltas instead of a peptide fragment database.

**C4. Shen, Y. & Bax, A. (2010)** "SPARTA+: a modest improvement in
empirical NMR chemical shift prediction by means of an artificial
neural network." *J. Biomol. NMR* 48, 13--22. DOI:
10.1007/s10858-010-9433-9

Neural network trained on 580 proteins: phi/psi/chi1 angles,
H-bond geometry, ring currents (Haigh-Mallion), electric fields.
The feature decomposition is the same as your kernels. SPARTA+
is the most widely used empirical predictor and the baseline your
system must exceed. Reports per-nucleus RMSDs: 1HN 0.49 ppm,
13CA 0.94 ppm, 15N 2.45 ppm.

C5. Han, B., Liu, Y., Ginzinger, S.W. & Wishart, D.S. (2011)
"SHIFTX2: significantly improved protein chemical shift prediction."
*J. Biomol. NMR* 50, 43--57. DOI: 10.1007/s10858-011-9478-4

Hybrid predictor: sequence-based + structure-based models trained on
>190 proteins. 26% better correlation than SPARTA+. Uses chi2/chi3
angles, solvent accessibility, H-bond geometry. The current
state-of-the-art empirical baseline.

C6. Neal, S., Nip, A.M., Zhang, H. & Wishart, D.S. (2003)
"Rapid and accurate calculation of protein 1H, 13C and 15N
chemical shifts." *J. Biomol. NMR* 26, 215--240.

Original SHIFTX: explicit decomposition into ring current, electric
field, H-bond, and solvent terms. Architecturally analogous to
your kernel decomposition but with scalar output.

C7. Kohlhoff, K.J., Robustelli, P., Cavalli, A., Salvatella, X. &
Vendruscolo, M. (2009) "Fast and accurate predictions of protein
NMR chemical shifts from interatomic distances." *J. Am. Chem. Soc.*
131, 13894--13895. DOI: 10.1021/ja903772t

CamShift: shifts as polynomial functions of interatomic distances,
differentiable for structure refinement. Distance-based formulation
connects to how your near-field kernels operate within cutoff radii.

C8. Sahakyan, A.B. & Vendruscolo, M. (2013) "Analysis of the
contributions of ring current and electric field effects to the
chemical shifts of RNA bases." *J. Phys. Chem. B* 117, 1989--1998.
DOI: 10.1021/jp3057306

Already in your `references/` directory. Decomposes RNA shifts into
ring current and electric field contributions, validates the
additive model, and provides empirical tests of the Haigh-Mallion
formulation against DFT.

---

## D. Shielding Tensor Theory and Spherical Tensors

The mathematical formalism for the rank-2 tensor output.
"T2 is sacred" — these papers define what that means.

**[CORE] D1. Haeberlen, U. (1976)** *High Resolution NMR in Solids:
Selective Averaging.* Academic Press. ISBN: 978-0120255610.

Defines the Haeberlen convention for shielding tensor principal
components and the irreducible decomposition into isotropic (rank-0),
antisymmetric (rank-1), and symmetric traceless (rank-2) parts.
The T0/T1/T2 decomposition throughout the codebase follows this
formalism. Your thesis must declare which convention it follows;
Haeberlen is the standard choice.

D2. Mehring, M. (1983) *Principles of High Resolution NMR in
Solids.* 2nd ed. Springer. ISBN: 978-3642687563.

Complementary treatment of spherical tensor representations for
NMR interactions. Covers Mehring notation for CSA parameters.
IUPAC recommends both Haeberlen and Mehring conventions; read
both to make a deliberate choice.

**[CORE] D3. Facelli, J.C. (2011)** "Chemical shift tensors: Theory
and application to molecular structural problems." *Prog. Nucl.
Magn. Reson. Spectrosc.* 58, 176--201. DOI:
10.1016/j.pnmrs.2010.10.003

Single best review of shielding tensor theory: irreducible
spherical tensor decomposition, GIAO methods, applications across
1H/13C/15N/17O. Covers the connection between tensor symmetry and
molecular structure. Your thesis theory chapter can lean heavily on
this. Section 3 on the irreducible decomposition is directly
relevant to why T2 (the symmetric traceless part) carries
geometric information that T0 (the isotropic shift) does not.

D4. Anet, F.A.L. & O'Leary, D.J. (1991) "The shielding tensor.
Part I: Understanding its symmetry properties." *Concepts Magn.
Reson.* 3, 193--214. DOI: 10.1002/cmr.1820030403

Tutorial on shielding tensor symmetry. Key point: NMR spectra are
sensitive only to the symmetric part; the antisymmetric part
contributes to relaxation only. Relevant for understanding what
your equivariant model can learn from shift data vs. what requires
relaxation data.

D5. Zare, R.N. (1988) *Angular Momentum: Understanding Spatial
Aspects in Chemistry and Physics.* Wiley. ISBN: 978-0471858928.

Textbook covering spherical tensor operators, Clebsch-Gordan
coefficients, and the Wigner-Eckart theorem with chemical physics
applications. The mathematical language for your L=0/L=2
decomposition and the equivariant tensor products in `model.py`.

D6. Silver, B.L. (1976) *Irreducible Tensor Methods: An Introduction
for Chemists.* Academic Press. ISBN: 978-0126437607.

More accessible introduction to irreducible tensors than Zare.
Particularly good for the Wigner-Eckart theorem's implications
for selection rules — which tensor products can contribute to which
output ranks.

**D7. Grisafi, A., Wilkins, D.M., Csanyi, G. & Ceriotti, M. (2018)**
"Symmetry-Adapted Machine Learning for Tensorial Properties of
Atomistic Systems." *Phys. Rev. Lett.* 120, 036002. DOI:
10.1103/PhysRevLett.120.036002

Derives tensor kernels adapted to rotational symmetry as the
natural generalisation of SOAP for predicting tensorial (not
just scalar) properties. This is the key theoretical bridge
between scalar molecular descriptors and full shielding tensor
prediction. Directly relevant to your T2 output: shows how to
build regression models that respect the transformation properties
of rank-2 tensor targets. Proves that you need equivariant
features to predict equivariant outputs.

D8. Grisafi, A. & Ceriotti, M. (2019) "Incorporating long-range
physics in atomic-scale machine learning." *J. Chem. Phys.* 151,
204105. DOI: 10.1063/1.5128375

Extends symmetry-adapted kernels to long-range electrostatic
contributions using multipole-like features. Relevant because your
Coulomb/Buckingham kernels are the physics-prescribed version of
what this paper learns.

---

## E. Equivariant Neural Networks

The architecture family for the correction head and the upstream
message-passing model.

**[CORE] E1. Thomas, N., Smidt, T., Kearnes, S. et al. (2018)**
"Tensor field networks: Rotation- and translation-equivariant neural
networks for 3D point clouds." *arXiv:1802.08219*.

The foundational architecture for SE(3)-equivariant molecular
networks. Introduces locally equivariant layers that accept and
output scalars (L=0), vectors (L=1), and higher-order tensors
(L=2, ...) using Clebsch-Gordan tensor products as the nonlinearity.
Your shielding tensor output (L=0 isotropic + L=2 anisotropy) is
naturally expressed in this framework. The theoretical basis for
why your `EquivariantCorrectionHead` uses `o3.FullyConnectedTensorProduct`.

**[CORE] E2. Geiger, M. & Smidt, T. (2022)** "e3nn: Euclidean Neural
Networks." *arXiv:2207.09453*.

The software framework your model uses. Implements E(3)-equivariant
operations: tensor products, spherical harmonics, Clebsch-Gordan
decompositions. Cite as the implementation reference for the tensor
products in `model.py`. The irrep string notation ("40x0e + 9x2e")
comes from e3nn.

E3. Kondor, R. & Trivedi, S. (2018) "On the generalization of
equivariance and convolution in neural networks to the action of
compact groups." *ICML 2018*. arXiv:1802.03690.

Proves that convolutional structure is both sufficient and
necessary for equivariance under compact group actions. The
theoretical justification for why Clebsch-Gordan tensor products
are the right nonlinearity — there is no other choice that
preserves SO(3) equivariance.

**[CORE] E4. Batzner, S., Musaelian, A., Sun, L. et al. (2022)**
"E(3)-equivariant graph neural networks for data-efficient and
accurate interatomic potentials." *Nat. Commun.* 13, 2453. DOI:
10.1038/s41467-022-29939-5

NequIP: E(3)-equivariant convolutions on geometric tensor features
achieve state-of-the-art with up to 1000x less training data. The
data efficiency argument is critical for your 720-mutant set where
each DFT calculation costs hours. NequIP demonstrates that
equivariant architectures extract more from limited data because
they don't waste capacity learning rotational symmetry.

E5. Batatia, I., Kovacs, D.P., Simm, G.N.C. et al. (2022) "MACE:
Higher Order Equivariant Message Passing Neural Networks for Fast
and Accurate Force Fields." *NeurIPS 2022*. arXiv:2206.07697.

Extends message passing to four-body equivariant interactions via
Clebsch-Gordan tensor products. The connection to the atomic
cluster expansion (ACE) makes this relevant to understanding which
body-order terms your kernels capture: ring currents are
fundamentally many-body (atom + ring geometry), which is why a
simple pair potential can't reproduce them.

E6. Fuchs, F.B., Worrall, D.E., Fischer, V. & Welling, M. (2020)
"SE(3)-Transformers: 3D Roto-Translation Equivariant Attention
Networks." *NeurIPS 2020*. arXiv:2006.10503.

Equivariant self-attention: invariant attention weights with
equivariant value embeddings. The attention mechanism offers an
alternative to message passing for learning which kernels matter
at which distances — related to the "scalar MLP produces per-kernel
weights" architecture of your `KernelMixingHead`.

E7. Musaelian, A., Batzner, S., Johansson, A. et al. (2023)
"Learning local equivariant representations for large-scale
atomistic dynamics." *Nat. Commun.* 14, 579. DOI:
10.1038/s41467-023-36329-y

Allegro: strictly local equivariant architecture (no global message
passing), demonstrated on 100M atoms. Relevant because your
calibration assumes per-atom independence — Allegro validates that
locality is not a crippling restriction for equivariant models.

E8. Schutt, K.T., Unke, O.T. & Gastegger, M. (2021) "Equivariant
message passing for the prediction of tensorial properties and
molecular spectra." *ICML 2021*. arXiv:2102.03150.

PaiNN: message passing extended to equivariant vector features.
Predicts rank-2 tensor quantities (dipole moments, polarisability).
The most directly relevant architecture for your output: PaiNN
shows how equivariant message passing produces tensor predictions,
which is what your upstream model needs to do.

E9. Bronstein, M.M., Bruna, J., Cohen, T. & Velickovic, P. (2021)
"Geometric Deep Learning: Grids, Groups, Graphs, Geodesics, and
Gauges." *arXiv:2104.13478*.

Comprehensive survey unifying CNNs, GNNs, Transformers under a
common geometric framework. The "5G" taxonomy provides vocabulary
for situating your work: you operate on graphs with SO(3) group
symmetry, and your kernels are gauge-equivariant by construction.
Good citation for the thesis introduction.

### Recent protein-SE(3)-equivariant applications

Surfaced in `spec/PHYSICS_FOUNDATIONS.md` Session 0.6. E1-E9
establish the equivariant-ML framework; E10-E14 capture the more
recent wave of SE(3)-equivariant architectures applied to
protein-scale prediction (interactions, binding sites, ternary
complexes). Our Stage 3 learning framework will sit alongside this
community even when specialised to shielding rather than
interaction prediction.

**E10. Satorras, V.G., Hoogeboom, E. & Welling, M. (2021)** "E(n)
Equivariant Graph Neural Networks." *ICML 2021*. arXiv:2102.09844.

EGNN: the simpler equivariant GNN that drops high-order tensor
products in favour of scalar/vector message passing. The baseline
that later architectures (VN-EGNN in E12, the SE(3)-equivariant
protein applications in E11 / E13) compare against. Under-cited in
E1-E9 given how frequently it is the reference point.

E11. (2025) "SE(3)-Equivariant Geometric Deep Learning for
Structure-Aware Prediction of Compound-Protein Interactions."
arXiv:2504.04654. (EquiCPI — author list to verify against arXiv.)

SE(3)-equivariant GNN applied to compound-protein interaction
prediction at 2025 scale. Evidence that the equivariant approach
has consolidated from potential-energy surfaces (NequIP / MACE)
toward protein-scale structural tasks. Paired with E13 as the
2025-vintage SE(3)-equivariant protein applications surfaced in
Session 0.6.

E12. "E(3)- and SE(3)-Equivariant Graph Neural Networks with
Virtual Nodes Enhance Protein Binding Site Identification."
(OpenReview.) VN-EGNN — authors, year, and venue to verify.

Virtual-node extension of EGNN: learnable nodes at non-atom
positions that preserve equivariance. Directly relevant to us
because virtual nodes are the clean way to evaluate classical
kernels at non-atom probe points (ring centres, bond midpoints,
volumetric grids) inside an equivariant learning framework. Worth
a close read if Stage 3 approaches a GNN-for-shielding-tensor
architecture that wants a thesis-viz probe-point vocabulary.

E13. (2025) "SE(3)-equivariant ternary complex prediction towards
target protein degradation." *Nat. Commun.* DOI:
10.1038/s41467-025-61272-5. (Author list to verify against article.)

SE(3)-equivariant prediction of protein-ligand-ligase ternary
complexes. Pairs with E11 as 2025 SE(3)-equivariant protein-scale
applications surfaced in Session 0.6. Different biological target,
same architectural family.

E14. Krapp, L.F. et al. (2023) "E(3) equivariant graph neural
networks for robust and accurate protein-protein interaction site
prediction." *PLOS Comput. Biol.* (Author list and DOI to verify.)

E(3)-equivariant GNN for protein-protein interaction site
prediction. Completes the E10-E14 arc: simple EGNN (E10) →
virtual-node extension (E12) → protein-compound application (E11)
→ protein-ligand-ligase ternary (E13) → protein-protein (E14).

### Equivariant tensor-valued and time-series prediction (community 6)

Surfaced in `spec/PHYSICS_FOUNDATIONS.md` Session 0.5 and chased
in depth 2026-04-22. E15-E23 cover the community our Stage 3
problem sits inside: predicting tensor-valued observables
(T0+T1+T2 irreps, not scalars) with equivariance preserved across
architectures that may be temporal (autoregressive trajectory,
neural-operator over time), spectral (frequency-dependent), or
static per-frame with downstream accumulation. Our
shielding-tensor-per-MD-frame task is most directly a per-frame
prediction at high cardinality (~10⁹ atom-frames per protein
trajectory) where the architecture choice decides whether the
T2 residual analysis at the heart of the thesis is well-posed.
J2 Venetos 2023 (silicate shift tensors) remains the closest
NMR-specific prior art; E15-E17 are the survey anchors,
E18-E20 are recent tensor-message-passing architectures, and
E21-E23 are the trajectory / time-series precedents.

**E15. Kondor, R. (2025)** "The principles behind equivariant
neural networks for physics and chemistry." *PNAS* 122(41),
e2415656122. DOI: 10.1073/pnas.2415656122. Open access
(PMC12541325).

Single-author pedagogical review tying the tensor-as-array vs
tensor-as-physicist distinction to irreducible representations.
Spherical harmonics as the natural generalisation of Fourier
series from circle to sphere; Clebsch-Gordan decomposition as
what turns generic multi-arrays into physically meaningful
irreps. The reference to cite when the thesis theory chapter
needs a single pedagogical anchor for why we track T0 / T1 / T2
separately. (Note: `PHYSICS_FOUNDATIONS.md` 0.4 initially
attributed this paper to "Smidt et al." — it is single-author
Kondor; 0.4 corrected 2026-04-22.)

E16. Gerken, J.E., Aronsson, J., Carlsson, O., Linander, H.,
Ohlsson, F., Petersson, C. & Persson, D. (2023) "Geometric deep
learning and equivariant neural networks." *Artif. Intell. Rev.*
56, 14605--14662. DOI: 10.1007/s10462-023-10502-7. arXiv:2105.13926.

Comprehensive survey of group-equivariant and gauge-equivariant
neural networks. Develops gauge-equivariant CNNs on arbitrary
manifolds via principal bundles and associated vector bundles.
Companion to E9 Bronstein et al. 2021 on the formal mathematical
side — where E9 gives the 5G taxonomy, E16 gives the derivations.

E17. Fei, J. & Deng, Z. (2024) "Rotation invariance and
equivariance in 3D deep learning: a survey." *Artif. Intell. Rev.*
DOI: 10.1007/s10462-024-10741-2. (Volume / article number to
verify.)

Survey scoped specifically to rotation invariance / equivariance
for 3D deep learning — adjacent to E16 but restricted to the
3D-structure case we actually care about. Organises architectures
by whether inputs, hidden layers, or outputs carry the
equivariance. Our Stage 3 approach is the output-equivariant
category.

**E18. Wang, J., Wang, Y., Zhang, H., Yang, Z., Liang, Z., Shi,
J., Wang, H.-T., Xing, D. & Sun, J. (2024)** "E(n)-Equivariant
Cartesian tensor message passing interatomic potential." *Nat.
Commun.* 15, 7607. DOI: 10.1038/s41467-024-51886-6. Open access.

HotPP: E(n)-equivariant message passing with arbitrary-order
Cartesian tensor node embeddings and messages. Three equivariant
operations: linear combinations of same-order tensors,
contraction of two tensors, partial derivative with respect to
another Cartesian tensor. Predicts dipole moments, polarisability
tensors, and IR/Raman spectra. Architectural precedent for
carrying T2+ tensor messages through a message-passing network in
Cartesian form rather than spherical-harmonic irreps —
relevant when the thesis considers architectural alternatives to
e3nn-style irrep decomposition.

**E19. Hsu, T.-W., Fang, Z., Bansil, A. & Yan, Q. (2026)**
"Accurate prediction of tensorial spectra using equivariant graph
neural network." *Nat. Commun.* DOI: 10.1038/s41467-026-69159-9.
arXiv:2505.04862.

TSENN: maps crystal structures to photon-frequency-dependent
optical (dielectric) tensors. Explicitly encodes isotropic
sequential scalar components and anisotropic sequential tensor
components into spherical tensor representations for crystalline
symmetry-aware prediction. The "anisotropic-sequential-tensors
approach" referred to in `PHYSICS_FOUNDATIONS.md` 0.5 is TSENN;
identity confirmed 2026-04-22. Trained on 1,432 bulk
semiconductors; MAE 0.127 on dielectric tensor. Closest conceptual
neighbour to our per-atom anisotropic shielding tensor trajectories
even though the sequence dimension differs (frequency for TSENN,
time for us). Code: qmatyanlab/TSENN.

E20. (2024) "Tensor improve equivariant graph neural network for
molecular dynamics prediction." *Comput. Biol. Chem.*? (Authors,
journal, DOI to verify against ScienceDirect pii S1476927124000410.)
Colloquial: TEGNN.

Tensor-Improved Equivariant GNN for MD prediction. Introduces
equivariant locally-complete frames into scalar-only equivariant
GNNs, allowing projection of tensor information of a given order
onto the frame with equivariant information transfer preserved.
Evaluated on N-body simulations and MD17. Citation details
surfaced only from web-search snippets during depth pass 2026-04-22;
primary record verification pending.

**E21. Xu, M., Han, J., Lou, A., Kossaifi, J., Ramanathan, A.,
Azizzadenesheli, K., Leskovec, J., Ermon, S. & Anandkumar, A.
(2024)** "Equivariant Graph Neural Operator for Modeling 3D
Dynamics." *ICML 2024*. arXiv:2401.11037.

EGNO: explicitly models dynamics as trajectories (not next-step
prediction) by formulating dynamics as a function over time and
learning neural operators to approximate it. Equivariant temporal
convolutions parameterised in Fourier space, SE(3)-equivariance
preserved across time steps. Applied to particle simulations,
human motion capture, and molecular dynamics. One of the clearer
architectural precedents for a community-6 learner that ingests a
shielding-tensor trajectory and produces time-resolved /
autocorrelation-function outputs. Code: MinkaiXu/EGNO.

E22. Thiemann, F.L., Reschützegger, T., Esposito, M., Taddese, T.,
Olarte-Plata, J.D. & Martelli, F. (2025) "Force-Free Molecular
Dynamics Through Autoregressive Equivariant Networks."
arXiv:2503.23794.

TrajCast: autoregressive equivariant message-passing network
that directly updates atomic positions and velocities without
solving equations of motion. Forecast intervals up to 30× larger
than traditional MD time-steps. Distinct from EGNO (which models
dynamics as neural operators on function spaces) by going
step-by-step but force-free. Recent trajectory-prediction
architecture family. Code at IBM/trajcast under CC-BY-NC-ND.

E23. Gregory, W.G., Tonelli-Cueto, J., Marshall, N.F., Lee, A.S.
& Villar, S. (2026) "Tensor learning with orthogonal, Lorentz,
and symplectic symmetries." *ICLR 2026*. arXiv:2406.01552.

Equivariant ML under classical Lie groups (orthogonal, Lorentz,
symplectic) combined with the path-signatures approach for
temporal tensor data. The only community-6 entry here that
treats tensor time series specifically through rough-path
theory rather than time-convolutions or autoregressive steps.
Complementary to E21 / E22 architectures on the mathematical side.
Original v1 title (June 2024) was "Learning equivariant tensor
functions with applications to sparse vector recovery"; retitled
for ICLR 2026 acceptance (v2, Feb 2026) — 0.5 referred to it by
the v1 title.

### Methods neighbours for Stage 3 (Lie-group GPs, SE(3) diffusion, neural operators)

Surfaced in `spec/PHYSICS_FOUNDATIONS.md` Session 0.7 methods
queries 2026-04-22 pm. E24-E28 are architectural neighbours a
Stage 3 learner might borrow from — mathematical cousins of the
equivariant-GNN architectures in E1-E23 rather than closely-
aligned protein-shielding precedents. Forward-looking survey;
inclusion is "committee would expect us to know these exist,"
not "these solve our problem."

**E24. Azangulov, I., Borovitskiy, V., Terenin, A., Mostowsky,
P., Deisenroth, M.P. & Durrande, N. (2024)** "Stationary Kernels
and Gaussian Processes on Lie Groups and their Homogeneous Spaces
I: the Compact Case." *J. Mach. Learn. Res.* 25, paper 23-0115.
arXiv:2208.14960. (Full author list to verify.)

Gaussian processes on SO(n), SE(n), and related Lie groups via
heat / Matérn kernels adapted to the group structure. Software at
geometric-kernels.github.io and github.com/imbirik/LieStationaryKernels.
**Calculator hint / Stage 3 alternative:** GP regression on SE(3)
configuration space is an alternative to equivariant-GNN for
learning the kernel-feature → shielding-tensor map, especially
where uncertainty quantification matters. Not a new classical
calculator — a Stage 3 model architecture neighbour paired with
F6 Chirikjian / F7 Helgason on the formal-theory side.

E25. Yim, J., Trippe, B.L., De Bortoli, V., Mathieu, E., Doucet,
A., Barzilay, R. & Jaakkola, T. (2023) "SE(3) diffusion model
with application to protein backbone generation." arXiv:2302.02277.
(Full author list to verify.)

FrameDiff: SE(3)-equivariant diffusion over multiple residue
frames for backbone generation. Architectural cousin to E10-E14
protein equivariant GNNs, but generative rather than
discriminative. Less directly applicable to our shielding-
prediction task; cited for completeness of the 2023-2024
SE(3)-for-proteins landscape and because the μs-harvester phase
may eventually incorporate conformational sampling via diffusion.

E26. Yim, J. et al. (2024) "Diffusion models in protein structure
and docking." *WIREs Comput. Mol. Sci.* DOI: 10.1002/wcms.1711.
(Co-authors to verify.)

Recent WIREs review of diffusion models in protein structure and
docking. Entry point to the E25 family of work. Skippable for the
immediate shielding pipeline; flagged as the up-to-date survey to
consult if/when diffusion sampling becomes part of the harvester
phase.

**E27. Li, Z., Kovachki, N., Azizzadenesheli, K., Liu, B.,
Bhattacharya, K., Stuart, A. & Anandkumar, A. (2021)** "Fourier
Neural Operator for Parametric Partial Differential Equations."
*ICLR 2021*. arXiv:2010.08895.

FNO: learns mappings between function spaces via Fourier-space
integral kernels. **Directly grounds E21 EGNO's architectural
ancestry:** EGNO's equivariant temporal convolution in Fourier
space inherits FNO's core insight. **Calculator hint:** if a
future volumetric-shielding-field extension is ever wanted (e.g.,
predicting the BS / HM butterfly isosurfaces that `QtBiotSavartCalc`
/ `QtHaighMallionCalc` in h5-reader currently render from
closed-form math), FNO applied to Cartesian shielding-field grids
is the natural architectural starting point — function-to-function
learning rather than per-atom regression.

E28. Lu, L., Jin, P., Pang, G., Zhang, Z. & Karniadakis, G.E.
(2021) "Learning nonlinear operators via DeepONet based on the
universal approximation theorem of operators." *Nat. Mach. Intell.*
3, 218--229. DOI: 10.1038/s42256-021-00302-5.

DeepONet: the neural-operator archetype alongside E27 FNO.
Branch-trunk decomposition for function-to-function learning.
Complementary to FNO with different inductive biases (branch
handles arbitrary discretisations; FNO handles periodic structure
natively). Together E27-E28 bracket the neural-operator design
space a Stage 3 volumetric-field learner might borrow from.

---

## F. Geometric Descriptors and Molecular Kernels

The descriptor design tradition your geometric kernels extend.

**[CORE] F1. Bartok, A.P., Kondor, R. & Csanyi, G. (2013)** "On
representing chemical environments." *Phys. Rev. B* 87, 184115. DOI:
10.1103/PhysRevB.87.184115

SOAP: atom-density correlation kernels built from spherical harmonics
and radial bases provide complete, smooth, rotationally invariant
descriptors. The theoretical ancestor of your approach. SOAP encodes
the same spatial information your kernels encode, but in a learned
rather than physics-prescribed basis. Key distinction: SOAP is
general-purpose and invariant; your kernels are physics-specific
and equivariant (they transform as tensors). Your thesis should
articulate why physics-prescribed equivariant kernels are preferable
to general-purpose invariant descriptors for this problem.

F2. Behler, J. & Parrinello, M. (2007) "Generalized neural-network
representation of high-dimensional potential-energy surfaces." *Phys.
Rev. Lett.* 98, 146401. DOI: 10.1103/PhysRevLett.98.146401

Atom-centered symmetry functions for neural network potentials.
Establishes the paradigm of local atomic descriptors feeding into
neural networks. The radial/angular symmetry function decomposition
maps to your distance/angle-dependent kernel terms. Historically
important as the paper that made ML interatomic potentials practical.

F3. Ceriotti, M., Willatt, M.J. & Csanyi, G. (2019) "Machine
learning of atomic-scale properties based on physical principles."
*Handbook of Materials Modeling.* Springer. arXiv:1901.10971.

Review connecting SOAP to other density-based descriptors. Practical
guidance on how body-order truncation, angular resolution, and
radial basis completeness affect accuracy. Applicable to your kernel
design choices (cutoff radius, ring type resolution, bond category
decomposition).

F4. Drautz, R. (2019) "Atomic cluster expansion for accurate and
transferable interatomic potentials." *Phys. Rev. B* 99, 014104.
DOI: 10.1103/PhysRevB.99.014104

ACE: complete, systematically improvable basis for local atomic
environments. MACE (E5) extends ACE to equivariant message passing.
Understanding ACE body-order helps rationalise which physical
effects need which descriptor complexity.

F5. Willatt, M.J., Musil, F. & Ceriotti, M. (2019) "Atom-density
representations for machine learning." *J. Chem. Phys.* 150, 154110.
DOI: 10.1063/1.5090481

Unifies SOAP, ACE, and moment tensors as projections of atomic
density. Shows how invariant and equivariant representations are
related. Helpful for understanding why your L=2 kernel vectors are
equivariant density projections in disguise.

### SE(3) and Lie-group formalism (robotics crossover)

Surfaced in `spec/PHYSICS_FOUNDATIONS.md` Session 0.6. F1-F5 cover
descriptor design built on spherical harmonics and atom-density
projections; F6-F7 ground those choices in the broader Lie-group
formalism the robotics and drug-discovery communities both use.
Our classical kernels are implicit functions on SE(3) configuration
space; this literature is the rigorous vocabulary.

F6. Chirikjian, G.S. (2009, 2012) *Stochastic Models, Information
Theory, and Lie Groups, Volumes I and II.* Birkhäuser / Springer.
(Vol. I ISBN 978-0-8176-4802-2, 2009; Vol. II ISBN
978-0-8176-4943-2, 2012. Verify editions / subtitles.)

Robotics-origin reference text for SE(3) kernel mathematics:
Gaussian processes on Lie groups, the exponential and log maps,
invariant measures, adjoint representation, stochastic differential
equations on manifolds. The formal machinery that underlies
SE(3)-equivariant GNN design (E10-E14) and that our classical
kernels embed implicitly. Library needed.

F7. Helgason, S. (2000) *Groups and Geometric Analysis: Integral
Geometry, Invariant Differential Operators, and Spherical
Functions.* American Mathematical Society, Mathematical Surveys
and Monographs vol 83. ISBN 978-0-8218-2673-7. (Originally
Academic Press 1984; verify reprint details.)

Foundational harmonic analysis on Lie groups. Covers the theoretical
apparatus (spherical functions, invariant differential operators,
the Plancherel formula on semisimple groups) that grounds SOAP
(F1), tensor field networks (E1), and e3nn (E2) at the level of
completeness / equivariance proofs — the textbook to cite if a
committee member asks why equivariant features form a closed space
under the SO(3) action. Library needed.

---

## G. GNN Message Passing for Molecular Properties

Architecture family for the upstream NMR shift predictor.

G1. Gilmer, J., Schoenholz, S.S., Riley, P.F., Vinyals, O. &
Dahl, G.E. (2017) "Neural message passing for quantum chemistry."
*ICML 2017*. arXiv:1704.01212.

Unifies GNNs under the MPNN framework: message, update, readout.
The abstraction your upstream model extends with SO(3)-valued
messages. Essential background citation.

G2. Schutt, K.T., Kindermans, P.J., Sauceda, H.E. et al. (2017)
"SchNet: A continuous-filter convolutional neural network for
modeling quantum interactions." *NeurIPS 2017*. arXiv:1706.08566.

Continuous-filter convolutions for atoms: learns distance-dependent
interaction functions. The invariant predecessor to PaiNN.
Understanding SchNet clarifies what equivariance adds.

G3. Gasteiger, J., Gross, J. & Gunnemann, S. (2020) "Directional
Message Passing for Molecular Graphs." *ICLR 2020*.
arXiv:2003.03123.

DimeNet: messages embed bond angles via spherical Bessel functions.
Intermediate step between invariant SchNet and fully equivariant
PaiNN/NequIP. Shows that angular information matters for molecular
properties, which your thesis demonstrates more explicitly via
T2 residuals.

G4. Klicpera, J., Becker, F. & Gunnemann, S. (2021) "GemNet:
Universal directional graph neural networks for molecules." *NeurIPS
2021*. arXiv:2106.08903.

Extends DimeNet to include dihedral angles (four-body interactions).
Relevant to the body-order discussion: some shielding contributions
(ring current from a specific ring on a specific atom) are
intrinsically many-body.

---

## H. MD Simulation for NMR Prediction

Using molecular dynamics to improve shift prediction through
conformational averaging.

**[CORE] H1. Markwick, P.R.L., Cervantes, C.F., Abel, B.L. et al.
(2010)** "Enhanced conformational space sampling improves the
prediction of chemical shifts in proteins." *J. Am. Chem. Soc.* 132,
1220--1221. DOI: 10.1021/ja9093692

Accelerated MD + SHIFTX: 20% improvement in shift RMSD over
standard MD, 28% over static X-ray. Demonstrates that conformational
averaging is not optional for accurate shift prediction. Your
approach of extracting kernels from 10 MD poses per RefDB structure
is directly motivated by this finding.

H2. Li, D.W. & Bruschweiler, R. (2010) "Certification of molecular
dynamics trajectories with NMR chemical shifts." *J. Phys. Chem.
Lett.* 1, 246--248. DOI: 10.1021/jz9001345

Back-calculated shifts as a quantitative score for validating MD
trajectories. Ensemble averages systematically improve over
individual snapshots. Provides the validation framework for your
short MD runs: if MD-averaged kernel features predict shifts better
than static-structure features, the MD is contributing real
conformational information.

H3. Robustelli, P., Stafford, K.A. & Palmer, A.G. III (2012)
"Interpreting protein structural dynamics from NMR chemical shifts."
*J. Am. Chem. Soc.* 134, 6365--6374. DOI: 10.1021/ja300265w

MD-averaged CamShift vs. experiment quantifies conformational
dynamics populations. Relevant to understanding what your
MD-averaged features encode beyond static structure.

H4. Robustelli, P., Piana, S. & Shaw, D.E. (2018) "Developing a
molecular dynamics force field for both folded and disordered protein
states." *Proc. Natl. Acad. Sci. USA* 115, E4758--E4766. DOI:
10.1073/pnas.1800690115

Benchmarks six force fields against >9000 experimental observables
including NMR shifts. Develops a99SB-disp. Critical for
understanding force-field-dependent systematic errors: your MD uses
ff14SB (via GROMACS topology), so understanding its biases matters.

H5. Yi, Y., Negrut, A., Bhate, M.P. & McDermott, A.E. (2024)
"Predicted and Experimental NMR Chemical Shifts at Variable
Temperatures." *J. Am. Chem. Soc.* 146, 16012--16024. DOI:
10.1021/jacs.4c02092

QM/MM over 20 MD snapshots: R = 0.73 for 15N, 0.94 for 13Calpha.
Recent validation that averaging over MD frames outperforms single
conformations. Directly supports your 10-pose averaging approach.

H6. Vendruscolo, M., Paci, E., Dobson, C.M. & Karplus, M. (2001)
"Three key residues form a critical contact network in a protein
folding transition state." *Nature* 409, 641--645. DOI:
10.1038/35054591

Uses MD with chemical shift restraints for structural refinement.
Establishes the bidirectional relationship between MD and NMR data
that your pipeline exploits in one direction (MD -> features ->
shifts).

### MD convergence and sampling quality

**[CORE] H7. Grossfield, A. & Zuckerman, D.M. (2009)** "Quantifying
uncertainty and sampling quality in biomolecular simulations."
*Annu. Rep. Comput. Chem.* 5, 23--48. DOI:
10.1016/S1574-1400(09)00502-7

Standard reference on MD convergence: block averaging, effective
sample size, decorrelation times, and the distinction between "long
simulation" and "converged observable." Every `TrajectoryResult`
field in the rollup spec section 5 needs a diagnostic from this
playbook — at 50 ns per protein, "we ran a 50 ns trajectory" is
not a convergence claim. Open access (PMC2865156).

H8. Lyman, E. & Zuckerman, D.M. (2007) "On the structural
convergence of biomolecular simulations by determination of the
effective sample size." *J. Phys. Chem. B* 111, 12876--12882.
DOI: 10.1021/jp073061t

Effective sample size (ESS) formalism: counts statistically
independent configurations via structural-histogram variance.
Directly applicable to our 50 ns × 600-frame fleet trajectories;
runnable as a standalone audit over the analysis H5s. Paywalled.

H9. Zhang, X., Bhatt, D. & Zuckerman, D.M. (2010) "Automated
sampling assessment for molecular simulations using the effective
sample size." *J. Chem. Theory Comput.* 6, 3048--3057. DOI:
10.1021/ct1002384

Practical automated ESS computation. Matched tooling to H8.
(Note: my earlier draft of this entry mis-attributed it to Lyman
& Zuckerman — the 2010 JCTC paper is Zhang–Bhatt–Zuckerman.
Corrected 2026-04-22.)

**H10. Showalter, S.A. & Brüschweiler, R. (2007)** "Validation of
molecular dynamics simulations of biomolecules using NMR spin
relaxation as benchmarks: application to the AMBER99SB force
field." *J. Chem. Theory Comput.* 3, 961--975. DOI:
10.1021/ct7000045

The NMR-relaxation benchmark for protein MD. Establishes that
back-calculated S² from MD vs. measured S² from NMR is the
validation standard, at multi-hundred-ns trajectory lengths. Our
ff14SB / CHARMM36m fleet at 50 ns sits at or below this baseline
protocol length — S² correlations against experiment expected in
the R²≈0.7 regime per H13. Paywalled.

H11. Romo, T.D. & Grossfield, A. (2011) "Block covariance overlap
method and convergence in molecular dynamics simulation." *J.
Chem. Theory Comput.* 7, 2464--2472. DOI: 10.1021/ct2002754

Method for error bars on trajectory averages via block covariance
overlap. Applicable to every covariance field in
`TrajectoryResult` (position_covariance, residue_residue_
covariance, per-observable cov). Paywalled.

H12. Romo, T.D. & Grossfield, A. (2014) "Unknown unknowns: the
challenge of systematic and statistical error in molecular
dynamics simulations." *Biophys. J.* 106, 1553--1554. DOI:
10.1016/j.bpj.2014.03.007. Open access (PMC4008789).

Conceptual anchor for the thesis honesty section — distinguishes
bounded-by-protocol errors from errors visible only via
cross-protocol comparison. Framework for scoping what
`TrajectoryResult` fields can claim vs. what the raw signal
unambiguously provides.

H13. Hoffmann, F., Mulder, F.A.A. & Schäfer, L.V. (2018)
"Narrowing the gap between experimental and computational
determination of methyl group dynamics in proteins." *Phys. Chem.
Chem. Phys.* 20, 24577--24590. DOI: 10.1039/C8CP03915A

State-of-the-art MD-vs-NMR methyl S² comparison: R²≈0.7, per-
residue median absolute deviation ~0.06 (CH3 S²axis). Reports that
μs-scale trajectories are required for full methyl S² convergence
but that current-force-field benchmarks locate the practical
"optimal" single-trajectory regime at ~50 ns. Our 50 ns MD window
sits in that zone for fast methyl motions; the μs regime remains
aspirational.

H15. (2018) "Validating Molecular Dynamics Simulations against
Experimental Observables in Light of Underlying Conformational
Ensembles." *J. Phys. Chem. B* 122. DOI: 10.1021/acs.jpcb.8b02144.
(Author list to verify.)

Package-comparison study (AMBER, GROMACS, NAMD, CHARMM-adjacent)
showing that different MD engines reproduce experimental
observables for engrailed homeodomain and RNase H equally well
overall, but with subtle differences in the underlying
conformational distributions — experiment cannot always
distinguish ensembles with matching observable means. **Load-
bearing caveat for T2-residual interpretation (MATHS_GOALS pillar
2):** when the residual is small, it may mean the ensemble-average
σ(T2) is right; it may also mean different ensembles produce the
same observable. Strengthens the motivation for `TrajectoryResult`
covariance fields alongside means (rollup 13.3b explicitly).

### CSA-averaged-over-MD parameterisation (Thread 2 surfaced)

H14. Jordan, J.B., Rule, G.S. & Tjandra, N. (2007)
"Parameterization of peptide 13C carbonyl chemical shielding
anisotropy in molecular dynamics simulations." *ChemPhysChem* 8,
1489--1495. DOI: 10.1002/cphc.200700003

DFT on 103 N-methylacetamide conformers → parametric surface →
CSA predictions along MD trajectories. Directly parallel to our
approach but with smaller body-order and smaller scale: our
calibrated kernel-based σ over full-protein MD is the
protein-scale, body-order-rich version. Precedent to cite
explicitly. Paywalled.

---

## I. NMR Reference Databases

The experimental data sources.

**[CORE] I1. Zhang, H., Neal, S. & Wishart, D.S. (2003)** "RefDB:
A database of uniformly referenced protein chemical shifts." *J.
Biomol. NMR* 25, 173--195. DOI: 10.1023/A:1022836027055

RefDB corrects systematic referencing errors in BMRB: ~25% of 13C
and 27% of 15N entries had significant referencing problems. Your
~700 curated RefDB structures inherit this quality control. The
re-referencing procedure is important to document because it means
your experimental targets are consistently referenced, unlike raw
BMRB entries.

I2. Ulrich, E.L. et al. (2008) "BioMagResBank." *Nucleic Acids
Res.* 36, D402--D408. DOI: 10.1093/nar/gkm957

The upstream data source for RefDB. Provenance citation for your
experimental chemical shift labels.

I3. Wishart, D.S. & Nip, A.M. (1998) "Protein chemical shift
analysis: a practical guide." *Biochem. Cell Biol.* 76, 153--163.
DOI: 10.1139/o98-038

Practical guide to chemical shift data quality. Useful for
understanding the noise floor in your experimental targets.

---

## J. ML for NMR Chemical Shifts

Direct precedents: machine learning applied to NMR prediction.

**[CORE] J1. Yang, Z., Chakraborty, M. & White, A.D. (2021)**
"Predicting chemical shifts with graph neural networks." *Chem. Sci.*
12, 10802--10809. DOI: 10.1039/D1SC01895G

GNNs for NMR shifts: captures H-bonding, secondary structure,
works on arbitrary molecular graphs. Computes 1M shifts in ~5
seconds. The direct GNN-for-NMR-shifts prior work. Your approach
adds two things: (1) physics-informed tensor kernel features instead
of learned-from-scratch, and (2) full tensor output instead of
scalar isotropic shifts.

**[CORE] J2. Venetos, M.C., Wen, M. & Persson, K.A. (2023)**
"Machine Learning Full NMR Chemical Shift Tensors of Silicon Oxides
with Equivariant Graph Neural Networks." *J. Phys. Chem. A* 127,
2388--2398. DOI: 10.1021/acs.jpca.2c07530

The closest precedent to your project: equivariant GNN predicting
full 29Si chemical shift tensors (not just isotropic shifts) in
silicates, achieving 1.05 ppm MAE. The key differences: (1) they
work on inorganic silicates, you work on proteins; (2) they learn
from scratch, you provide physics-prescribed kernel features;
(3) your T2 angular residual diagnostic has no analogue in their
work.

J3. Guan, Y., Sowndarya, S.V.S., Gallegos, L.C., St. John, P.C.
& Paton, R.S. (2021) "Real-time prediction of 1H and 13C chemical
shifts with DFT accuracy using a 3D graph neural network." *Chem.
Sci.* 12, 12012--12026. DOI: 10.1039/D1SC03343C

3D GNN predicting small-molecule shifts to DFT accuracy (~0.2 ppm
for 1H) in milliseconds. Demonstrates that geometric features
carry enough information to reproduce DFT. Different domain
(small molecules vs proteins) but validates the concept.

J4. Unke, O.T., Chmiela, S., Sauceda, H.E. et al. (2021) "Machine
learning force fields." *Chem. Rev.* 121, 10142--10186. DOI:
10.1021/acs.chemrev.0c01111

Comprehensive review of ML potentials: symmetry functions, message
passing, equivariant architectures, training strategies. Though
focused on forces, the architecture and symmetry discussions apply
directly to your tensor prediction problem. Good reference for the
thesis background on why equivariance matters.

J5. Li, J. & Harrison, R.J. (2023) "NMR chemical shift prediction
using a message passing neural network." *J. Chem. Inf. Model.* 63,
2710--2718. DOI: 10.1021/acs.jcim.3c00165

MPNN for 1H/13C shifts with attention to local chemical environment.
Standard GNN approach without equivariance. Useful as an ablation
comparator: does equivariance help for shift prediction?

J6. (2024) "Accurate Prediction of NMR Chemical Shifts:
Integrating DFT Calculations with Three-Dimensional Graph Neural
Networks." *J. Chem. Theory Comput.* DOI: 10.1021/acs.jctc.4c00422.
(Author list to verify.)

DFT + 3D GNN hybrid for chemical shift prediction. **Closest-
prior-art on the Stage 3 architecture axis:** our kernel-feature +
ridge / e3nn calibration pipeline is in the same design family
(physics-informed features + GNN-style learner, trained against
DFT). J6 positions us against a 2024 benchmark; J2 Venetos 2023
(silicates, full tensor) remains closest on the "full tensor
output" axis. Together J2 and J6 bracket the Stage 3 design space.

---

## K. Force Fields and Molecular Simulation Infrastructure

Background citations for the simulation layer.

K1. Maier, J.A. et al. (2015) "ff14SB: improving the accuracy of
protein side chain and backbone parameters from ff99SB." *J. Chem.
Theory Comput.* 11, 3696--3713. DOI: 10.1021/acs.jctc.5b00255

The force field providing partial charges for your Coulomb kernel
(via `ChargeAssignmentResult`). Understanding ff14SB charge
parameterisation is necessary for understanding the quality of your
electrostatic features.

K2. Brooks, B.R. et al. (1983) "CHARMM: A program for macromolecular
energy, minimization, and dynamics calculations." *J. Comput. Chem.*
4, 187--217. DOI: 10.1002/jcc.540040211

Source of the switching function used in `DispersionResult.cpp`
(lines 42-48). The smooth taper from R_switch to R_cut prevents
discontinuities in the dispersion kernel.

K3. Stewart, J.J.P. (2013) "Optimization of parameters for
semiempirical methods VI: more modifications to the NDDO
approximations and re-optimization of parameters." *J. Mol. Model.*
19, 1--32. DOI: 10.1007/s00894-012-1667-x

MOPAC PM7: the semiempirical method providing Mulliken charges and
Wiberg bond orders for `MopacCoulombResult` and
`MopacMcConnellResult`. PM7+MOZYME scales linearly for large
systems.

K4. Abraham, M.J. et al. (2015) "GROMACS: High performance molecular
simulations through multi-level parallelism from laptops to
supercomputers." *SoftwareX* 1-2, 19--25. DOI:
10.1016/j.softx.2015.06.001

The MD engine. Your `GromacsEnsembleLoader` reads TPR topology and
trajectory data from GROMACS.

K5. Řezáč, J. & Hobza, P. (2016) "Semiempirical Quantum
Mechanical Methods for Noncovalent Interactions for Chemical and
Biochemical Applications." *Chem. Rev.* 116, 5038--5071. DOI:
10.1021/acs.chemrev.5b00584. (Author list / volume / page range
to verify.)

Comprehensive Chem. Rev. on semiempirical methods (PM7, OMx,
DFTB, GFN-xTB) for noncovalent interactions in biomolecular
systems. **Context for K3 Stewart 2013 PM7:** grounds the
expected quality of MOPAC-derived Mulliken charges and Wiberg
bond orders that feed MopacCoulombResult and MopacMcConnellResult
(Calculators 9 and 10 in GEOMETRIC_KERNEL_CATALOGUE.md). Scopes
what PM7 captures (protonation-dependent charge redistribution,
noncovalent geometry) vs what it misses (transition-metal
coordination, dispersion beyond D4).

---

## L. Protein NMR and Structure

Broader context for why NMR shielding prediction matters.

L1. Williamson, M.P. (1990) "Secondary-structure dependent chemical
shifts in proteins." *Biopolymers* 29, 1423--1431. DOI:
10.1002/bip.360291009

Establishes the relationship between secondary structure and
chemical shifts. Shows that alpha-helices and beta-sheets have
characteristic shift patterns. Background for why your backbone
features (which you deliberately exclude from the calibration
model because "backbone doesn't change" in mechanical mutants)
would matter for the upstream predictor.

L2. Wishart, D.S. & Sykes, B.D. (1994) "Chemical shifts as a tool
for structure determination." *Methods Enzymol.* 239, 363--392. DOI:
10.1016/S0076-6879(94)39014-2

Review establishing chemical shifts as structural restraints.
The motivation for your entire project: if shifts encode structure,
and geometric kernels encode the physics of that encoding, then
calibrated kernels provide better structure-shift relationships.

L3. Cavalli, A., Salvatella, X., Dobson, C.M. & Vendruscolo, M.
(2007) "Protein structure determination from NMR chemical shifts."
*Proc. Natl. Acad. Sci. USA* 104, 9615--9620. DOI:
10.1073/pnas.0610313104

CS-ROSETTA predecessor: determines protein structure using only
chemical shifts. Demonstrates that shifts contain enough
information for full structure determination, motivating the
development of better shift predictors.

L4. Shen, Y. et al. (2008) "Consistent blind protein structure
generation from NMR chemical shift data." *Proc. Natl. Acad. Sci.
USA* 105, 4685--4690. DOI: 10.1073/pnas.0800256105

CS-ROSETTA: automated protein structure determination from shifts.
The downstream application that benefits from more accurate shift
prediction. Better kernels -> better shifts -> better structures.

L5. Robustelli, P., Kohlhoff, K., Cavalli, A. & Vendruscolo, M.
(2010) "Using NMR chemical shifts as structural restraints in
molecular dynamics simulations of proteins." *Structure* 18,
923--933. DOI: 10.1016/j.str.2010.04.016

Uses back-calculated shifts as MD restraints. The bidirectional
MD-NMR relationship: shifts guide MD, MD improves shifts. Your
pipeline exploits this in the forward direction.

### CSA tensors measured in proteins (Thread 2 surfaced)

Per-residue CSA tensor principal components are the experimental
bench for our predicted σ T2. These references establish what is
measured, at what resolution, and in which test proteins.

**[CORE] L6. Wylie, B.J., Sperling, L.J., Nieuwkoop, A.J.,
Franks, W.T., Oldfield, E. & Rienstra, C.M. (2011)** "Ultrahigh
resolution protein structures using NMR chemical shift tensors."
*Proc. Natl. Acad. Sci. USA* 108, 16974--16979. DOI:
10.1073/pnas.1103728108. Open access (PMC3193204).

Per-residue measured CSA tensors used as direct structural
restraints, resolving structure at ultra-high resolution.
Demonstrates that tensor principal values carry geometric
information beyond isotropic shifts — **the same thesis claim
our T2 work makes, independently tested by the solid-state
community**. The closest conceptual precedent for why T2 is the
thesis argument.

L7. Yao, L., Grishaev, A., Cornilescu, G. & Bax, A. (2010)
"Site-specific backbone amide 15N chemical shift anisotropy
tensors in a small protein from liquid crystal and
cross-correlated relaxation measurements." *J. Am. Chem. Soc.*
132, 4295--4309. DOI: 10.1021/ja910186u

Per-residue amide ¹⁵N CSA tensor components for GB3 from two
orthogonal methods (liquid-crystal alignment + CCR). Finds
α-helix ¹⁵N CSA ≈ −173 ± 7 ppm, β-sheet ≈ −162 ± 6 ppm. Direct
bench for rollup W3 (NH DD/CSA CCR): the prediction depends on
¹⁵N CSA magnitude, so our per-residue CSA must reproduce this
~10 ppm secondary-structure split for the CCR rate to be
trustworthy. Paywalled.

L8. Cornilescu, G. & Bax, A. (2000) "Measurement of proton,
nitrogen, and carbonyl chemical shielding anisotropies in a
protein dissolved in a dilute liquid crystalline phase." *J. Am.
Chem. Soc.* 122, 10143--10154. DOI: 10.1021/ja0016194

Simultaneous 1H, 15N, 13C' CSAs measured in the same protein via
liquid-crystal alignment. Multi-nucleus per-residue benchmark.
Paywalled.

L9. Wylie, B.J., Franks, W.T. & Rienstra, C.M. (2006)
"Determinations of 15N chemical shift anisotropy magnitudes in a
uniformly 15N,13C-labeled microcrystalline protein by 3D
magic-angle spinning nuclear magnetic resonance spectroscopy."
*J. Phys. Chem. B* 110, 10926--10936. DOI: 10.1021/jp060507h

Site-resolved 15N CSA magnitudes for GB1 via 3D MAS NMR.
Companion to L6 and key reference for the carbonyl CSA finding
that **δ₂₂ of ¹³C' depends significantly on CO···HN H-bond
length** — direct test target for our HBond kernel + McConnell
kernel combined contribution to carbonyl C' σ. Paywalled.

L10. Sitkoff, D. & Case, D.A. (1998) "Theories of chemical shift
anisotropies in proteins and nucleic acids." *Prog. Nucl. Magn.
Reson. Spectrosc.* 32, 165--190. DOI: 10.1016/S0079-6565(98)00013-2

Case-group review of ab initio CSA in biomolecules: GIAO
methodology, conformation dependence, H-bond effects, solvent
contributions. Paired with B1 de Dios-Pearson-Oldfield 1993 for
the ab initio-CSA lineage. **Tethered to T2-residual diagnostic
(MATHS_GOALS pillar 2):** this is the ab initio landscape our
per-atom σ(T0+T1+T2) ridge-calibrated prediction sits inside; the
L=2 component is what the review treats as "the tensor" and what
our angular-residual exposes. Paywalled.

L11. Brender, J.R., Taylor, D.M. & Ramamoorthy, A. (2001)
"Orientation of Amide-Nitrogen-15 Chemical Shift Tensors in
Peptides: A Quantum Chemical Study." *J. Am. Chem. Soc.* 123,
914--922. DOI: 10.1021/ja001980q

DFT study establishing the orientation of the ¹⁵N CSA principal
axis system relative to the N-H bond vector in model peptides.
Companion to L7 Yao-Bax 2010 (CSA magnitudes) on the orientation
side. **Tethered to rollup W3 CCR derivation:** the N-H dipolar /
¹⁵N CSA CCR rate back-calculation requires both tensor magnitudes
AND tensor orientations relative to the bond — L7 gives
magnitudes, L11 gives orientations, our pipeline produces both
via per-frame σ tensor + per-frame N-H bond vector. Paywalled.

**[CORE] L12. Loth, K., Pelupessy, P. & Bodenhausen, G. (2005)**
"Chemical Shift Anisotropy Tensors of Carbonyl, Nitrogen, and
Amide Proton Nuclei in Proteins through Cross-Correlated
Relaxation in NMR Spectroscopy." *J. Am. Chem. Soc.* 127,
6062--6068. DOI: 10.1021/ja042863o

Canonical solution-state amide CSA dataset: principal components
and orientations of C', ¹⁵N, and HN CSA tensors at 64 amide bonds
in human ubiquitin, from 14 complementary auto- and cross-
correlated relaxation rates. The closest-scale external validation
bench for our per-atom σ prediction on an individual test protein:
ubiquitin is in the calibration set, L12 reports per-residue
tensors, and the observable is back-calculable from per-frame
σ(T0+T1+T2) + per-frame bond-vectors. **Directly tethered to 13.2
R1/R2/R1ρ row and rollup W3 (`nh_dipolar_csa_ccr_rate`).**
Paywalled.

---

## M. Protein Dynamics & NMR Relaxation

Motion-visible observables and the theory underpinning them. This
section is the theoretical backing for the `TrajectoryResult`
object in `spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md`
section 5, and for honest scoping of what the 50 ns MD window can
and cannot observe. Added 2026-04-22 during Session 0 literature
pass for `spec/PHYSICS_FOUNDATIONS.md`.

### Foundational relaxation theory

**[CORE] M1. Lipari, G. & Szabo, A. (1982)** "Model-free approach
to the interpretation of nuclear magnetic resonance relaxation in
macromolecules. 1. Theory and range of validity; 2. Analysis of
experimental results." *J. Am. Chem. Soc.* 104, 4546--4559 and
4559--4570. DOI: 10.1021/ja00381a009 and 10.1021/ja00381a010

The S² / τe model-free framework. Defines the generalized order
parameter S² as the long-time limit of the bond-vector orientation
autocorrelation function, and τe as the effective correlation time
of the motion decaying to that limit. Our
`bond_vector_autocorrelation` and derived `order_parameter_S2` /
`lipari_szabo_tau_e` in `TrajectoryResult` (rollup section 5) are
the MD-computed version of these quantities. Validity hinges on
timescale separation between internal motion and overall tumbling.
**Library needed** (two-part JACS 1982).

M2. Palmer, A.G. III (2004) "NMR characterization of the dynamics
of biomacromolecules." *Chem. Rev.* 104, 3623--3640. DOI:
10.1021/cr030413t

Review covering spin-relaxation-based protein dynamics from ps to
ms. Canonical citation for the timescale taxonomy that should
scope every `TrajectoryResult` field. Paywalled.

M3. Palmer, A.G. III (2014) "Chemical exchange in biomacromolecules:
past, present, and future." *J. Magn. Reson.* 241, 3--17. DOI:
10.1016/j.jmr.2014.01.008

CPMG / CEST / R1ρ for μs-ms chemical exchange dynamics. Defines
the exchange-regime observables our 50 ns window cannot reach by
rate, though populations at slow interconversion may be partially
accessible. Paywalled.

**M4. Henzler-Wildman, K. & Kern, D. (2007)** "Dynamic
personalities of proteins." *Nature* 450, 964--972. DOI:
10.1038/nature06522

Canonical review of the "function = dynamics" framing. Concept
anchor for rollup section 13.1 opening — what MD brings that
static QM cannot. **Library needed** (paywalled *Nature*).

**[CORE] M5. Frauenfelder, H., Sligar, S.G. & Wolynes, P.G.
(1991)** "The energy landscapes and motions of proteins." *Science*
254, 1598--1603. DOI: 10.1126/science.1749933

Conformational substates view: proteins as ensembles of nearly
isoenergetic basins with hierarchical motions between them.
Conceptual anchor for our `chi_density_grid`, `phi_psi_density`,
and substate-population interpretations in `TrajectoryResult`.
**Library needed** (paywalled *Science*).

### Order parameter as entropy meter — the Wand programme

M6. Wand, A.J. (2001) "Dynamic activation of protein function: a
view emerging from NMR spectroscopy." *Nat. Struct. Biol.* 8,
926--931. DOI: 10.1038/nsb1101-926

Establishes measured S² as a probe of protein dynamics and
function. Paywalled.

**M7. Kasinath, V. & Wand, A.J. (2013)** "Microscopic insights
into the NMR relaxation-based protein conformational entropy
meter." *J. Am. Chem. Soc.* 135, 15092--15100. DOI:
10.1021/ja405200u

The calmodulin calibration that converts methyl S² into a ΔS
entropy estimate. Validation opportunity for our predicted S²
from 50 ns MD against published calmodulin + target-peptide S²
datasets. Paywalled.

**M8. Frederick, K.K., Marlow, M.S., Valentine, K.G. & Wand, A.J.
(2007)** "Conformational entropy in molecular recognition by
proteins." *Nature* 448, 325--329. DOI: 10.1038/nature05959

The paper that made conformational entropy measurable via NMR.
**Library needed** (paywalled *Nature*). Related open-access
companion at PMC4156320.

M9. Sharp, K.A., O'Brien, E., Kasinath, V. & Wand, A.J. (2015)
"On the relationship between NMR-derived amide order parameters
and protein backbone entropy changes." *Proteins* 83, 922--930.
DOI: 10.1002/prot.24789. Open access (PMC4400257).

Develops a calibration curve for backbone ΔS vs. Δ⟨O²NH⟩ via MD
ensembles. Methodology for connecting our predicted S² fields to
a backbone-entropy story, companion to the methyl-entropy meter
(M7).

### Cross-correlated relaxation

**[CORE] M10. Tugarinov, V. & Kay, L.E. (2003)** "Cross-correlated
relaxation enhanced 1H-13C NMR spectroscopy of methyl groups in
very high molecular weight proteins and protein complexes." *J.
Am. Chem. Soc.* 125, 13868--13878. DOI: 10.1021/ja030153x

Methyl TROSY exploits intra-methyl DD / DD cross-correlation.
Back-calculating CCR rates from an MD trajectory requires per-frame
σ tensor coupled to dipolar bond vectors, correctly time-
correlated. **Our pipeline produces calibrated per-frame
σ(T0+T1+T2) at ~4000-atom scale — directly testable against
published methyl CCR rates.** Paywalled.

M11. Ferrage, F., Piserchio, A., Cowburn, D. & Ghose, R. (2008)
"Direct measurement of dipole-dipole/CSA cross-correlated
relaxation by a constant-time experiment." *J. Magn. Reson.* 192,
302--313. DOI: 10.1016/j.jmr.2008.03.011. Open access
(PMC2542487).

Experimental methodology for DD/CSA CCR. Back-calculation target
dataset for our pipeline.

M12. Fushman, D., Tjandra, N. & Cowburn, D. (1998) "Direct
measurement of 15N chemical shift anisotropy in solution." *J.
Am. Chem. Soc.* 120, 10947--10952. DOI: 10.1021/ja981961e

Amide N-15 CSA measured in solution, residue by residue. Test
case for whether our classical-kernel predictions reproduce the
observed residue-to-residue CSA variation. Paywalled.

### RDCs and ensemble dynamics

**M13. Tolman, J.R., Flanagan, J.M., Kennedy, M.A. & Prestegard,
J.H. (1995)** "Nuclear magnetic dipole interactions in
field-oriented proteins: information for structure determination
in solution." *Proc. Natl. Acad. Sci. USA* 92, 9279--9283. DOI:
10.1073/pnas.92.20.9279. Open access (PMC40968).

Seminal RDC paper. Establishes that alignment-media orientation
information carries structure-and-dynamics content. Paywalled
(old PNAS).

M14. Meiler, J., Prompers, J.J., Peti, W., Griesinger, C. &
Brüschweiler, R. (2001) "Model-free approach to the dynamic
interpretation of residual dipolar couplings in globular proteins."
*J. Am. Chem. Soc.* 123, 6098--6107. DOI: 10.1021/ja010002z

RDC-based model-free. Connects bond-vector orientation tensor to
measured RDCs. Direct test target for our proposed
`bond_orientation_tensor (n_bonds, 3, 3)` field (rollup 13.3a).
Paywalled.

M15. Lakomek, N.A., Walter, K.F., Farès, C., Lange, O.F.,
de Groot, B.L., Grubmüller, H., Brüschweiler, R., Munk, A.,
Becker, S., Meiler, J. & Griesinger, C. (2008) "Self-consistent
residual dipolar coupling based model-free analysis for the robust
determination of nanosecond to microsecond protein dynamics." *J.
Biomol. NMR* 41, 139--155. DOI: 10.1007/s10858-008-9244-4. Open
access (PMC2480484).

The SCRM (self-consistent RDC-based model-free) method — removes
the influence of structural noise in RDC model-free analysis.
Method applicable to our MD-derived orientation tensors; direct
test target for `bond_orientation_tensor`.

M16. Prestegard, J.H., Bougault, C.M. & Kishore, A.I. (2004)
"Residual dipolar couplings in structure determination of
biomolecules." *Chem. Rev.* 104, 3519--3540. DOI: 10.1021/cr030419i

Comprehensive RDC review. Paywalled.

### Ring flips — slow dynamics our window cannot observe

M17. Akke, M. & Weininger, U. (2023) "NMR studies of aromatic
ring flips to probe conformational fluctuations in proteins."
*J. Phys. Chem. B* 127, 591--599. DOI: 10.1021/acs.jpcb.2c07258.
Open access (PMC9884080).

Current review of ring-flip NMR. Slow flips k ≈ 10¹-10² s⁻¹ are
measurable via CEST / CPMG; fast flips ≫10³ s⁻¹ are inferred from
chemical-shift equivalence but remain rate-hidden. **Our 50 ns MD
window sees neither regime's kinetics.** Honesty bound on
`ring_normal_autocorrelation` claims in rollup section 5: captures
sub-ns ring-normal breathing only.

M18. Skrynnikov, N.R. et al. (2013) "Slow aromatic ring flips
detected despite near-degenerate NMR frequencies of the exchanging
nuclei." *J. Phys. Chem. B* 117, 9886--9900. DOI: 10.1021/jp4058065

CEST-based detection of slow ring flips at near-degenerate shifts.
Paywalled.

### Essential dynamics and collective motion

**[CORE] M19. Amadei, A., Linssen, A.B.M. & Berendsen, H.J.C.
(1993)** "Essential dynamics of proteins." *Proteins* 17,
412--425. DOI: 10.1002/prot.340170408

PCA of MD as essential dynamics. Most of the trajectory variance
lives in ~20 modes. Backing for our `pca_eigenmodes` /
`pca_eigenvalues` fields in `TrajectoryResult`. **Library needed**
(pre-2000 journal).

M20. David, C.C. & Jacobs, D.J. (2014) "Principal component
analysis: a method for determining the essential dynamics of
proteins." *Methods Mol. Biol.* 1084, 193--226. DOI:
10.1007/978-1-62703-658-0_11. Open access (PMC4676806).

Practical tutorial with scripts. Implementation reference for the
PCA reduction.

### Memory kernel / non-Markovian

M21. Kou, S.C. & Xie, X.S. (2004) "Generalized Langevin equation
with fractional Gaussian noise: subdiffusion within a single
protein molecule." *Phys. Rev. Lett.* 93, 180603. DOI:
10.1103/PhysRevLett.93.180603

Evidence for non-exponential (power-law) memory in single-molecule
protein dynamics. Justifies shipping `bond_vector_autocorrelation`
as a function of lags rather than reducing to a single
Lipari-Szabo τe. Paywalled.

M22. Zwanzig, R. (2001) *Nonequilibrium Statistical Mechanics.*
Oxford University Press. ISBN: 978-0195140187.

Projection-operator / generalized Langevin treatment. Canonical
textbook for memory-kernel formalism. **Library needed** (book).

### Dynamic allostery

**M23. Cooper, A. & Dryden, D.T.F. (1984)** "Allostery without
conformational change: a plausible model." *Eur. Biophys. J.* 11,
103--109. DOI: 10.1007/BF00276625

Origin paper: allostery through entropy, not structural change.
The conceptual backing for why `residue_residue_covariance`
remains load-bearing even when mean positions are static.
**Library needed** (pre-2000 journal).

### Paramagnetic NMR — parallel angular-kernel literature

M24. Bertini, I., Luchinat, C. & Parigi, G. (2002) "Magnetic
susceptibility in paramagnetic NMR." *Prog. Nucl. Magn. Reson.
Spectrosc.* 40, 249--273. DOI: 10.1016/S0079-6565(02)00002-X

Canonical reference for magnetic-susceptibility-tensor and
pseudocontact-shift (PCS) formalism. PCS is structurally
identical to our `RingSusceptibilityResult` angular kernel ((3cos²θ
- 1)/r³ axial term plus rhombic term), with a different physical
source (paramagnetic centre vs. aromatic ring susceptibility).
Parallel validation opportunity via proteins with lanthanide tags
at measured positions. Paywalled.

M25. Nitsche, C. & Otting, G. (2021) "Pseudocontact shifts in
biomolecular NMR spectroscopy." *Chem. Rev.* 121, 8150--8172.
DOI: 10.1021/acs.chemrev.1c00796

Modern PCS review. Paywalled.

M26. Orton, H.W., Huber, T. & Otting, G. (2020) "Paramagpy:
software for fitting magnetic susceptibility tensors using
paramagnetic effects measured in NMR spectra." *Magn. Reson.* 1,
1--12. DOI: 10.5194/mr-1-1-2020.

Tool for fitting χ tensor from experimental PCS / PRE / CCR data —
the inverse problem analogue to our DFT-to-kernel fit. Open access.

### Stat-mech foundations (TCF, linear response, memory kernels — community 4)

Surfaced in `spec/PHYSICS_FOUNDATIONS.md` Session 0.3 and flagged
in 0.8 as a genuinely under-covered community at the foundational
level: section M already carried Lipari-Szabo (M1) and the
memory-kernel applications (M21 Kou-Xie, M22 Zwanzig book), but
not the primary stat-mech papers those sit on top of. Added
2026-04-22 during community-4 depth pass. These four entries
ground the `TrajectoryResult` autocorrelation fields (rollup
section 5) in the stat-mech theory that makes them physically
meaningful quantities rather than summary statistics.

**M27. Kubo, R. (1957)** "Statistical-mechanical theory of
irreversible processes. I. General theory and simple applications
to magnetic and conduction problems." *J. Phys. Soc. Jpn.* 12,
570--586. DOI: 10.1143/JPSJ.12.570

Canonical paper for linear-response theory and the
fluctuation-dissipation theorem. Derives how equilibrium
time-correlation functions determine macroscopic response
coefficients. The stat-mech foundation under every
`TrajectoryResult` autocorrelation field — our bond-vector /
kernel / σ-tensor autocorrelations are FDT-adjacent quantities
by construction. Paywalled; library needed.

**M28. Berne, B.J. & Harp, G.D. (1970)** "On the Calculation of
Time Correlation Functions." *Adv. Chem. Phys.* 17, 63--227. DOI:
10.1002/9780470143636.ch3

Canonical review chapter establishing the TCF formalism in MD.
Systematic treatment of correlation-function computation from
trajectory data: estimators, error bars, decorrelation times,
sampling-density requirements. Directly relevant to how we
compute `bond_vector_autocorrelation` and related
`TrajectoryResult` fields from finite 50 ns trajectories. A
companion short paper "Time-Correlation Functions, Memory
Functions, and Molecular Dynamics" *Phys. Rev. Lett.* / *Phys.
Rev. A* 2, 975 (1970) (DOI: 10.1103/PhysRevA.2.975) treats the
memory-kernel side. Library needed.

**M29. Zwanzig, R. (1961)** "Memory Effects in Irreversible
Thermodynamics." *Phys. Rev.* 124, 983--992. DOI:
10.1103/PhysRev.124.983

The primary memory-function paper underlying M22 (Zwanzig 2001
textbook). Derives the non-Markovian generalised kinetic equation
where response follows cause by a time-convolution with a
memory kernel that is itself a TCF in the rates of change of
phase functions. Rationalises why single-τ Lipari-Szabo (M1) is
a limit case and why our `bond_vector_autocorrelation` field
ships as a function of lag rather than as a collapsed single
time. Complements M21 Kou-Xie 2004 (experimental power-law memory
evidence) and M22 Zwanzig 2001 book (textbook formalism).
Paywalled; library needed.

M30. Ditler, J. & Luber, S. (2022) "Vibrational spectroscopy by
means of first-principles molecular dynamics simulations." *WIREs
Comput. Mol. Sci.* 12, e1605. DOI: 10.1002/wcms.1605

Recent review of the compute-observable-per-frame +
autocorrelate-to-get-spectrum workflow for IR, VCD, Raman, ROA,
SFG, and nonlinear spectroscopies from first-principles MD. The
methodological sibling of our pipeline one step across from NMR:
the dipole / polarisability / transition-moment tensor computed
per frame and autocorrelated becomes a vibrational spectrum,
where the shielding tensor computed per frame and autocorrelated
becomes CCR rates and related observables (rollup W3
derivations). Architectural precedent for how the
spectroscopy-from-MD community handles convergence, finite-time
artefacts, and tensor-valued autocorrelations. Open access via
ZORA (UZH institutional repository).

**M31. Green, M.S. (1954)** "Markoff Random Processes and the
Statistical Mechanics of Time-Dependent Phenomena. II.
Irreversible Processes in Fluids." *J. Chem. Phys.* 22,
398--413. DOI: 10.1063/1.1740082

The "Green" in "Green-Kubo." Derives transport coefficients
(viscosity, diffusion, heat conductivity) as time integrals of
equilibrium autocorrelation functions. Companion to M27 Kubo 1957.
**Tethered to rollup W3 and 13.2 R1/R2/R1ρ row:** our σ-tensor
autocorrelation → CCR rate is formally a Green-Kubo calculation of
an NMR-relaxation-adjacent coefficient; spectral densities J(ω) at
Larmor frequencies are Fourier transforms of the σ autocorrelation.
**Calculator hint:** a `GreenKuboSpectralDensityResult` (per-atom
J(ω) at named Larmor frequencies, computed from the per-atom σ(t)
autocorrelation as the Fourier-transform pair) would be close to
free given the autocorrelation is already the planned
`TrajectoryResult` field, and would feed R1 / R2 / R1ρ predictions
directly per 13.2. Paywalled; library needed.

---

## N. Magnetic Biophysics — Bulk Protein χ Anisotropy and PCS Validation

The per-atom magnetic-anisotropy kernels (McConnell,
RingSusceptibility, HBond) have a macroscopic sibling in the
MRI / biophysics tradition: bulk protein χ anisotropy computed
from structural subunits. Added 2026-04-22 during Session 0
thread-3 literature pass — surfaced as genuinely under-covered by
the bibliography cross-check in `spec/PHYSICS_FOUNDATIONS.md`
section 0.8.

### Diamagnetic anisotropy at bulk protein scale

**[CORE] N1. Worcester, D.L. (1978)** "Structural origins of
diamagnetic anisotropy in proteins." *Proc. Natl. Acad. Sci. USA*
75, 5475--5477. DOI: 10.1073/pnas.75.11.5475. Open access
(PMC392987; PMID 281695).

Foundational structural analysis of per-atom / per-bond
contributions to bulk protein χ anisotropy. α-helix large
anisotropy traced to axial peptide-bond alignment; aromatic-ring
contributions identified. Same physics as our McConnell and
ring-current kernels summed macroscopically. (Title was
previously mis-attributed to Pauling in earlier `PHYSICS_
FOUNDATIONS.md` 0.1 and in one memory note; corrected
2026-04-22.)

**[CORE] N2. Pauling, L. (1979)** "Diamagnetic anisotropy of the
peptide group." *Proc. Natl. Acad. Sci. USA* 76, 2293--2294.
DOI: 10.1073/pnas.76.5.2293. Open access (PMC383585).

Simple theory of planar peptide-group anisotropy. Gives the
classical molar value Δχ_peptide = −5.36 × 10⁻⁶ cm·g⁻¹·sec⁻¹
emu — a direct quantitative benchmark for our calibrated
McConnell peptide-CO / PeptideCN intensities under an appropriate
unit conversion.

N3. Babaei, M., Jones, I.C., Dayal, K. & Mauter, M.S. (2017)
"Computing the diamagnetic susceptibility and diamagnetic
anisotropy of membrane proteins from structural subunits." *J.
Chem. Theory Comput.* 13(6), 2945--2953. DOI:
10.1021/acs.jctc.6b01251.

Modern hierarchical method: per-amino-acid → secondary →
tertiary χ tensors summed in a unit coordinate system to give
bulk protein χ. **Directly parallel to our per-atom kernel
strategy, summed for MRI applications.** Reports that
aromatic-residue content and β-barrel topology dominate volumetric
χ anisotropy. Paywalled. **Authorship corrected 2026-04-22** —
earlier drafts of `PHYSICS_FOUNDATIONS.md` 0.1 and the conversation
carried a stale "Palmer et al." attribution; the correct first
author is Babaei (group at Carnegie Mellon). Note: Arthur Palmer
(M2, M3) is a separate NMR-relaxation theorist — not this group.

N4. Li, W., Liu, C., Duong, T.Q., van Zijl, P.C.M. & Li, X.
(2017) "Susceptibility tensor imaging (STI) of the brain." *NMR
in Biomedicine* 30(4), e3540. DOI: 10.1002/nbm.3540. Open access
(PMC5083244).

Review of STI methodology. Key claim for our purposes: molecular
magnetic susceptibility is the tensor sum of per-bond
contributions; in myofibers (unlike brain white matter) peptide-
bond χ anisotropy dominates the macroscopic signal. Bridges our
McConnell peptide kernel to an entirely different measurement
tradition. Authorship verified 2026-04-22.

### Paramagnetic PCS validation datasets (parallel angular kernel)

Cross-reference: M24 Bertini-Luchinat-Parigi 2002, M25 Nitsche-
Otting 2021, M26 Orton-Huber-Otting Paramagpy 2020 — formalism
and tooling.

A cluster of published datasets gives explicit lanthanide-tag χ
tensor parameters on standard test proteins (ubiquitin, GB1). They
enable forward-prediction of PCSs at backbone positions from our
RingSusceptibility-style angular kernel, and comparison to measured
— an external validation independent of our DFT oracle.
Representative entries:

N5. Graham, B. et al. (2011) "DOTA-amide lanthanide tag for
reliable generation of pseudocontact shifts in protein NMR
spectra" on ubiquitin A28C. *Bioconjugate Chem.* (PubMed 21877751;
DOI 10.1021/bc200353c). Paywalled.

N6. Joss, D. et al. (2019) "A sterically overcrowded,
isopropyl-substituted, lanthanide-chelating tag for protein PCS
NMR spectroscopy" benchmarked on ubiquitin S57C and hCA II S166C.
*Chem. Eur. J.* (PubMed 31199526; DOI 10.1002/chem.201901692).
Paywalled.

N7. Tharayil, S.M., Mahawaththa, M.C., Loh, C.-T., Adekoya, I.
& Otting, G. (2021) "Phosphoserine for the generation of
lanthanide-binding sites on proteins for paramagnetic nuclear
magnetic resonance spectroscopy." *Magn. Reson.* 2, 1--13. Open
access (PMC10539748).

Phosphoserine-based Ln-binding sites on ubiquitin and GB1; full
χ tensor parameters published per mutant. Direct external
benchmark for PCS forward-prediction from our RingSusceptibility
angular kernel given a known external χ tensor. Authorship
verified 2026-04-22.

(N5 and N6 authorship remain tentative — hand-verify when those
entries are elevated to specific validation targets.)

---

## Summary: Core 12 References

These define the equations, the data, the architecture, and the
closest prior art. A thesis committee member who reads only these
12 papers would understand what you built and why.

| # | Ref | Why core |
|---|-----|----------|
| 1 | A2 Haigh-Mallion 1979 | Ring current surface integral — your HM kernel |
| 2 | A9 McConnell 1957 | Bond anisotropy dipolar kernel — your MC kernel |
| 3 | A11 Buckingham 1960 | E-field/EFG shielding — your Coulomb kernel |
| 4 | B1 de Dios/Oldfield 1993 | Ab initio proof that shifts decompose into your kernel terms |
| 5 | C1 Case 1995 | The experiment yours extends from scalar to tensor |
| 6 | D1 Haeberlen 1976 | T0/T1/T2 decomposition convention |
| 7 | D3 Facelli 2011 | Shielding tensor review — thesis theory chapter |
| 8 | E1 Thomas 2018 | Tensor field networks — equivariant architecture basis |
| 9 | E2 e3nn 2022 | The framework your model uses |
| 10 | E4 NequIP 2022 | Data efficiency of equivariant models |
| 11 | F1 Bartok/SOAP 2013 | Geometric descriptor theory your kernels extend |
| 12 | H1 Markwick 2010 | MD averaging improves shift prediction |
| 13 | I1 RefDB 2003 | Your experimental data source |
| 14 | J2 Venetos 2023 | Closest prior art: equivariant GNN for full shift tensors |

Add C4 (SPARTA+) or C5 (SHIFTX2) as the empirical baseline
comparator and D7 (Grisafi tensor kernels) if the committee has
a methods focus. That gives 14-16 core citations.
