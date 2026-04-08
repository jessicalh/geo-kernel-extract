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
