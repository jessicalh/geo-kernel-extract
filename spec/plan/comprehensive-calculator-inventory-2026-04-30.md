# Comprehensive Calculator Inventory — nmr-extract (present and planned)

**Date:** 2026-04-30
**Status:** Authored for adversarial-agent vetting. Lists every calculator that exists in `src/` plus every calculator named as planned across the running notes, the PLANNED_CALCULATORS docs, the literature pass, and the cache findings. Each entry: plain English sentence + physics sentence + references.

**Vetting target:** an adversarial research agent should pull each cited reference, confirm the physics statement is correct (sign, units, derivation), and flag any discrepancy. References include local-cache filenames where applicable; otherwise the standard citation suffices.

**Reading conventions:**
- **PRESENT** — calculator exists in `src/`. File name listed.
- **PLANNED** — calculator is in PLANNED_CALCULATORS docs, running_plan_notes, or literature pass; not yet implemented.
- **ORPHAN** — schema slot in H5 with no producer in `src/`. User direction (2026-04-30): "we probably want those" — treat as PLANNED.
- **EXISTING / IMPROVEMENT** — improvement to an existing calculator from literature; named here as a calculator-pass item.

---

## Section 1 — Topology and setup utilities (PRESENT)

These are not "physics calculators" in the geometric-kernel sense; they prepare typed substrate the kernels need. Listed for completeness.

### 1.1 `GeometryResult` (PRESENT, `src/GeometryResult.h`)

**Plain English:** Computes ring centers / ring normals (SVD-based), bond lengths, bond midpoints, bond directions, and global geometric properties for the conformation. First result attached to any conformation; required by everything downstream.

**Physics:** Pure Euclidean geometry — ring center is the centroid of ring atom positions, ring normal is the smallest singular vector of the centered ring atom matrix, bond midpoint is the mean of bonded atom positions. No electronic structure, no field theory.

**References:** None required (linear algebra).

### 1.2 `SpatialIndexResult` (PRESENT, `src/SpatialIndexResult.h`)

**Plain English:** Builds nanoflann KD-trees over atom positions, ring centers, and bond midpoints; provides radius search and k-nearest queries; populates 15 Å neighbour lists per atom.

**Physics:** Spatial indexing is a computational utility, not physics. The 15 Å cutoff is the project constitution's ring-current cutoff radius (`spec/CONSTITUTION.md`).

**References:** nanoflann library (Blanco, https://github.com/jlblancoc/nanoflann), descended from Muja & Lowe's FLANN.

### 1.3 `EnrichmentResult` (PRESENT, `src/EnrichmentResult.h`)

**Plain English:** Assigns per-atom AtomRole (backbone/sidechain/aromatic/etc.) and categorical booleans using the typed bond graph and backbone index cache.

**Physics:** None. Pure topology classification using the typed substrate (no string comparisons after PDB load).

**References:** None.

### 1.4 `MolecularGraphResult` (PRESENT, `src/MolecularGraphResult.h`)

**Plain English:** BFS traversal of the bond graph to compute through-bond features per atom: graph distance to nearest ring / N / O, electronegativity sums at depth 1 and 2, π-bond count within 3 hops, conjugation flag, BFS decay weight.

**Physics:** Through-bond electronic communication via π-conjugation is real; the BFS decay (exp(-d / 4.0 bonds)) is a heuristic for through-bond effect attenuation, not a derived constant.

**References:** Pauling 1932 (electronegativity scale, used by `eneg_sum_*`); Mulliken 1934 (alternative scale).

### 1.5 `BuildResult` (PRESENT, `src/BuildResult.h`)

**Plain English:** Result-class wrapper around the build/load pipeline outputs (Protein, ProteinConformation, charges, topology). Plumbing for the typed object model.

**Physics:** None.

**References:** None.

### 1.6 `ProtonationDetectionResult` (PRESENT, `src/ProtonationDetectionResult.h`)

**Plain English:** Reports the prepared protonation state of titratable residues (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/ARN) by checking the presence of titratable hydrogens in the typed substrate.

**Physics:** Protonation states map onto distinct chemical species at the residue level: HIS protonation alters the imidazole ring electron density, CYX (disulfide) vs CYS (thiol) alters sulfur σ⁻ partial charge, ASP↔ASH vs GLU↔GLH alters carboxylate vs carboxylic-acid character, LYN vs LYS alters amine vs ammonium character.

**References:** Markley 1998 IUPAC nomenclature for NMR (atom-name conventions).

### 1.7 `DemoResult` (PRESENT, `src/DemoResult.h`)

**Plain English:** Trivial test calculator (per-atom distance to nearest ring center). Exists as a worked example of the ConformationResult pattern.

**Physics:** None — just `argmin(d(atom, ring_center_k))` over rings.

**References:** None.

---

## Section 2 — Charge sources (PRESENT)

### 2.1 `ChargeAssignmentResult` (PRESENT, `src/ChargeAssignmentResult.h`)

**Plain English:** Projects the prepared `ForceFieldChargeTable` (ff14SB / CHARMM36m / AMBER prmtop) onto each `ConformationAtom`, providing per-atom partial charge (in elementary charges) and PB radius.

**Physics:** Force-field partial charges are derived empirically (RESP fits to HF/6-31G(d) electrostatic potentials for AMBER ff14SB; CGenFF for CHARMM36m). They are a lookup, not a per-conformation calculation; CoulombResult and ApbsFieldResult are the consumers that compute fields from them.

**References:** Maier et al. 2015 *J Chem Theory Comput* 11, 3696 (ff14SB); Vanommeslaeghe et al. 2010 *J Comput Chem* 31, 671 (CHARMM36); Bayly et al. 1993 *J Phys Chem* 97, 10269 (RESP).

### 2.2 `EeqResult` (PRESENT, `src/EeqResult.h`)

**Plain English:** Computes geometry-dependent partial charges via extended electronegativity equilibration (D4 parameters), solving a single linear system per frame for charges that minimise `E(q) = Σ χᵢqᵢ + ½Σ ηᵢqᵢ² + ½Σ qᵢqⱼγ(Rᵢⱼ)` subject to charge neutrality.

**Physics:** EEQ is a charge-equalization model: each atom has electronegativity χ and chemical hardness η; charges flow until chemical potential is equalised, with Coulomb interaction via the Ohno-Klopman kernel `γ(R) = 1/√(R² + 1/(ηᵢ·ηⱼ))`. CN-dependent shifts κ make χ context-aware. Output is a fast (no QM, no SCF) per-conformation charge set.

**References:** Caldeweyher et al. 2019 *J Chem Phys* 150, 154122 (EEQ + D4); Klopman 1964 *J Am Chem Soc* 86, 4550 (Ohno-Klopman kernel); Ohno 1964 *Theor Chim Acta* 2, 219.

### 2.3 `MopacResult` (PRESENT, `src/MopacResult.h`)

**Plain English:** Runs PM7 + MOZYME semiempirical QM via libmopac in-process, providing per-atom Mulliken charges, orbital populations (s, p), Wiberg bond orders (continuous, conformation-dependent), heat of formation, and dipole moment.

**Physics:** PM7 is a parameterised neglect-of-diatomic-differential-overlap (NDDO) Hamiltonian fit to thousands of experimental enthalpies and geometries. MOZYME enables linear-scaling SCF for large biomolecules. The Mulliken charges and Wiberg bond orders are conformation-dependent QM observables (not force-field constants), capturing the real electronic redistribution under the local geometry.

**References:** Stewart 2013 *J Mol Model* 19, 1 (PM7); Stewart 1996 *Int J Quantum Chem* 58, 133 (MOZYME); Wiberg 1968 *Tetrahedron* 24, 1083 (Wiberg bond order); Mulliken 1955 *J Chem Phys* 23, 1833 (Mulliken population).

### 2.4 `AIMNet2Result` (PRESENT, `src/AIMNet2Result.h`)

**Plain English:** Runs AIMNet2 neural-network potential via libtorch on a CUDA GPU, producing per-atom Hirshfeld charges, an atomic environment embedding (the "aim" vector), and a Coulomb EFG tensor decomposed by source (backbone, aromatic) using the same dipolar kernel as CoulombResult.

**Physics:** AIMNet2 is a message-passing neural network trained on millions of conformations across the organic-chemistry chemical space (~25M conformations / ~2M molecules in the published model) at ωB97M-D3/def2-TZVPP level. The Hirshfeld charges are QM-derived per-conformation; the embedding is a learned atom-in-molecule descriptor. EFG comes from these charges via the same Stone T-tensor formalism as CoulombResult.

**References:** Anstine et al. 2024 *J Chem Inf Model* 64, 4382 (AIMNet2); Hirshfeld 1977 *Theor Chim Acta* 44, 129 (Hirshfeld partitioning); Behler & Parrinello 2007 *Phys Rev Lett* 98, 146401 (NN potential foundation).

---

## Section 3 — External tool wrappers (PRESENT)

### 3.1 `ApbsFieldResult` (PRESENT, `src/ApbsFieldResult.h`)

**Plain English:** Runs APBS to solve the linearised Poisson-Boltzmann equation in the protein + solvent system (charges from `ChargeAssignmentResult`, PB radii from same), then extracts per-atom electric field (V/Å) and electric field gradient (V/Å²) by central-difference interpolation on the potential grid; both are traceless-projected.

**Physics:** Linearised PB is `∇·(ε(r)∇φ(r)) - ε_s κ²(r) φ(r) = -ρ_protein(r)`, where ε(r) is the position-dependent dielectric (protein interior ~2-4, solvent ~78), κ is the Debye-Hückel screening parameter, and ρ_protein is the protein charge distribution. The solution φ(r) is the potential the ions and water see; E = -∇φ; EFG = -∇⊗∇φ. This is the solvent-mediated electric field — the physics that vacuum Coulomb misses.

**References:** Baker et al. 2001 *PNAS* 98, 10037 (APBS); Honig & Nicholls 1995 *Science* 268, 1144 (PB for biology); Buckingham 1960 *Can J Chem* 38, 300 (E-field shielding coupling).

### 3.2 `OrcaShieldingResult` (PRESENT, `src/OrcaShieldingResult.h`)

**Plain English:** Loads per-atom DFT shielding tensors (diamagnetic, paramagnetic, total) from an external ORCA NMR calculation output file and stores them as Mat3 + SphericalTensor on each `ConformationAtom`.

**Physics:** DFT shielding is `σ_αβ = ∂²E/∂B_α∂μ_β` evaluated at the SCF-converged electronic structure with GIAO (gauge-including atomic orbitals) for translational invariance. Diamagnetic part is the closed-form expectation value over the ground state; paramagnetic part involves the response to the magnetic perturbation via CPSCF. This is the ground-truth target the geometric kernels are calibrated against.

**References:** Neese 2012 *WIREs Comput Mol Sci* 2, 73 (ORCA); Ditchfield 1974 *Mol Phys* 27, 789 (GIAO origin); Helgaker et al. 1999 *Chem Rev* 99, 293 (NMR DFT review); Schreckenbach & Ziegler 1995 *J Phys Chem* 99, 606 (DFT shielding implementation).

---

## Section 4 — Geometric kernels (the load-bearing physics) (PRESENT)

### 4.1 `BiotSavartResult` (PRESENT, `src/BiotSavartResult.h`)

**Plain English:** For each aromatic ring within 15 Å of each atom, computes the ring-current B-field via the Johnson-Bovey double-loop wire-segment model and projects it onto the ring normal to produce a rank-2 shielding tensor `G_ab = -n_b · B_a · PPM_FACTOR`; full T0+T1+T2 decomposition is preserved per Boyd-Skrynnikov 2002 eq 3.

**Physics:** Aromatic π-electrons circulate in response to an external field B₀ perpendicular to the ring plane, producing a secondary B-field at nearby nuclei via classical Biot-Savart. The Johnson-Bovey model represents the current as two loops at ±d above/below the ring plane (each carrying I/2) with Biot-Savart line integral over wire segments. The rank-1 outer product `n ⊗ B` (with appropriate sign and factor) is the full ring-current shielding tensor; T0 = (n·B)·PPM_FACTOR/3 (the classic isotropic shift); T2 carries the angular anisotropy.

**References:** Johnson & Bovey 1958 *J Chem Phys* 29, 1012 (double-loop model); Boyd & Skrynnikov 2002 *JACS* 124, 1832 (full tensor formula, σ_xz extension; `references-text/boyd-skrynnikov-2002-ring-current-chemical-shielding-anisotropy-text-1.txt`); Pople 1956 *J Chem Phys* 24, 1111 (point-dipole precursor); Case 1995 *J Biomol NMR* 6, 341 (calibration vs DFT).

### 4.2 `HaighMallionResult` (PRESENT, `src/HaighMallionResult.h`)

**Plain English:** Same physics target as BiotSavart but uses surface integration: for each ring, integrates the symmetric traceless dipolar kernel `(3 ρ_a ρ_b / ρ⁵ - δ_ab / ρ³)` over the ring's surface via fan triangulation + Gaussian quadrature with adaptive subdivision, then constructs the shielding tensor as `G_ab = n_b · (H · n)_a`.

**Physics:** Same diamagnetic ring-current physics as BiotSavart but solved as a surface integral over the ring's enclosed area rather than a line integral over the wire. The surface formulation is more accurate at short distances where the wire singularity matters; both should agree at large distance. Whether they agree on the T2 anisotropy across the dynamic range is the open empirical question per the calculator docstring.

**References:** Haigh & Mallion 1979 *Prog NMR Spectrosc* 13, 303 (HM model); Haigh & Mallion 1972 *Mol Phys* 22, 955 (early formulation); Boyd & Skrynnikov 2002 *JACS* 124, 1832 (cross-check via tensor structure).

### 4.3 `RingSusceptibilityResult` (PRESENT, `src/RingSusceptibilityResult.h`)

**Plain English:** Treats each aromatic ring as a single point magnetic dipole at its center (with dipole moment along the ring normal) and computes the dipolar shielding tensor at each atom; same kernel structure as McConnellResult but with bond direction replaced by ring normal.

**Physics:** The ring's bulk diamagnetic susceptibility anisotropy Δχ produces a magnetic moment μ = Δχ · B₀ / (4π) along the ring normal under an external field. The induced field at a nearby atom is the standard magnetic dipole field; the shielding contribution is `M_ab/r³` with `M_ab = 9 cosθ d̂_a n_b - 3 n_a n_b - (3 d̂_a d̂_b - δ_ab)` per the geometric kernel catalogue. T0+T1+T2 all non-zero; asymmetric.

**References:** Pople-Schneider-Bernstein 1959 *High-Resolution NMR* (point-dipole ring current); Pascal 1910 *Ann Chim Phys* (diamagnetic susceptibility); Lazzeretti 2000 *Prog NMR Spectrosc* 36, 1 (ring current review).

### 4.4 `McConnellResult` (PRESENT, `src/McConnellResult.h`)

**Plain English:** For each bond within 10 Å of each atom, computes the full McConnell shielding tensor (asymmetric, non-traceless, T0+T1+T2 all non-zero) and the symmetric traceless dipolar kernel; tracks per-category totals (PeptideCO, PeptideCN, BackboneOther, SidechainCO, Aromatic) and nearest CO/CN.

**Physics:** A bond with magnetic susceptibility anisotropy Δχ_∥ - Δχ_⊥ acts as a point magnetic dipole at the bond midpoint, oriented along the bond axis. Induced field at nearby nuclei follows the dipolar law; shielding is `M_ab/r³` with M derived in the project's GEOMETRIC_KERNEL_CATALOGUE.md. The McConnell scalar `f = (3 cos²θ - 1)/r³` is the double contraction of K with the bond direction, NOT the tensor trace — the trace comes from the full M.

**References:** McConnell 1957 *J Chem Phys* 27, 226 (original); Pople 1956 *J Chem Phys* 24, 1111 (dipole approximation context); ApSimon & Beierbeck 1971 *J Chem Soc Perkin Trans 2* 1655 (parameterisation).

### 4.5 `MopacMcConnellResult` (PRESENT, `src/MopacMcConnellResult.h`)

**Plain English:** Same McConnell kernel as McConnellResult but each bond's contribution is weighted by its MOPAC Wiberg bond order; the model learns Δχ per category for this weighted-tensor variant.

**Physics:** A double bond (Wiberg ~1.8) has different magnetic susceptibility anisotropy than a single bond (Wiberg ~1.0); peptide CO with partial double-bond character (Wiberg 1.4-1.5) sits in between. Weighting the McConnell sum by per-bond Wiberg order modulates the angular pattern by the actual electron-density-derived bond character rather than a fixed force-field assumption.

**References:** McConnell 1957 *J Chem Phys* 27, 226; Wiberg 1968 *Tetrahedron* 24, 1083 (Wiberg index); Stewart 2013 (MOPAC PM7).

### 4.6 `CoulombResult` (PRESENT, `src/CoulombResult.h`) — **gated off in trajectory mode by `skip_coulomb=true`**

**Plain English:** For each atom, sums the vacuum Coulomb electric field E and electric field gradient tensor V from all other partial charges; decomposes by source (backbone, sidechain, aromatic) and computes solvent contribution as `apbs_efield − vacuum_E_total`. Disabled in `PerFrameExtractionSet` because APBS supersedes at N > 1000 atoms.

**Physics:** Standard Coulomb electrostatics: `E_a(i) = k_e Σ q_j · r_a / r³`, `V_ab(i) = k_e Σ q_j · (3 r_a r_b / r⁵ − δ_ab / r³)`. Each EFG term is traceless (Gauss's law); the sum is traceless after floating-point projection. In a vacuum (no dielectric screening) the field is unmediated; in solution, ions and waters reorient to screen it (which is why APBS-supersedes is the per-frame default).

**References:** Buckingham 1960 *Can J Chem* 38, 300 (`references-text/buckingham-1960-chemical-shifts-polar-groups-text-{1-3}.txt`; σ^(1)·E formula and coefficient); Stone 2013 *The Theory of Intermolecular Forces* OUP, ch. 3 (T-tensor formalism); Case 1995 *J Biomol NMR* 6, 341 (-3.0 ppm/au for nucleic acids).

### 4.7 `MopacCoulombResult` (PRESENT, `src/MopacCoulombResult.h`) — **trajectory mode only when MOPAC ran**

**Plain English:** Same Coulomb kernel as CoulombResult but reads `mopac_charge` (PM7 Mulliken, conformation-dependent) instead of fixed force-field `partial_charge`; the T2 angular pattern differs because MOPAC charges respond to the local electronic environment.

**Physics:** The Coulomb kernel is unchanged; only the charge source changes from force-field-fixed to QM-redistributed-per-conformation. A polarisable atom in a hydrogen-bonded environment carries different charge than the same atom in vacuo; the resulting field pattern at neighbours encodes this electronic response.

**References:** Same as CoulombResult plus Stewart 2013 (PM7 charges).

### 4.8 `PiQuadrupoleResult` (PRESENT, `src/PiQuadrupoleResult.h`)

**Plain English:** For each aromatic ring within range of each atom, treats the π-electron cloud as a point axial quadrupole at the ring center (along the ring normal) and computes the resulting EFG tensor at the atom; pure T2 (symmetric, traceless).

**Physics:** The π-electron density above/below the ring plane is a quadrupolar distribution, parameterised as Θ_cd = Θ(3 n_c n_d − δ_cd)/2. The potential is `Φ = Θ(3 cos²θ − 1)/(2r³)`, and the EFG `V_ab = -∂²Φ/∂x_a∂x_b` follows via Stone's rank-4 T-tensor with the trace contraction vanishing by Laplace (∇²(1/r) = 0 outside sources). Leading radial decay is 1/r⁵. The Buckingham γ-coefficient couples this EFG to ³C/¹⁵N/¹⁷O shielding via `Δσ ∝ γ · V`.

**References:** Stone 2013 *The Theory of Intermolecular Forces* OUP, ch. 3 (T-tensor); Buckingham 1960 *Can J Chem* 38, 300 (γ coefficient for EFG-shielding coupling); Sahakyan & Vendruscolo 2013 *J Phys Chem B* 117, 1989 (`references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-{1-4}.txt`; π-electron-quadrupole framing for RNA ¹⁷O).

### 4.9 `DispersionResult` (PRESENT, `src/DispersionResult.h`)

**Plain English:** For each aromatic ring within range of each atom, sums anisotropic dispersion tensors from each ring vertex atom using `K_ab = S(r) · (3 d_a d_b / r⁸ − δ_ab / r⁶)` with a CHARMM-style smooth switching function from 4.3 Å to 5.0 Å; ring-vertex through-bond exclusion via the bond graph.

**Physics:** London dispersion (instantaneous-dipole / induced-dipole) attractive interaction goes as 1/r⁶ in energy; its anisotropic contribution to shielding goes via the polarisability tensor coupling. The 1/r⁶ scalar captures the isotropic attractive part; the 1/r⁸ rank-2 part captures the angular dependence from the discrete vertex positions. The CHARMM C¹ switching function ensures continuity of forces across MD frames (no discontinuity at cutoff).

**References:** London 1937 *Trans Faraday Soc* 33, 8 (dispersion theory); Brooks et al. 1983 *J Comput Chem* 4, 187 (CHARMM switching function); Grimme et al. 2010 *J Chem Phys* 132, 154104 (modern dispersion corrections).

### 4.10 `HBondResult` (PRESENT, `src/HBondResult.h`)

**Plain English:** For each atom involved in a backbone H-bond identified by DSSP, computes the dipolar shielding tensor from the donor→acceptor partner using the same kernel structure as McConnell with bond direction replaced by H-bond vector.

**Physics:** Backbone N–H···O=C hydrogen bonds carry a partial covalent character; the donor-acceptor magnetic anisotropy contributes to nearby nucleus shielding via the dipolar kernel. H-bond identification uses the Kabsch-Sander electrostatic energy criterion (E = q·q · [1/r_NO + 1/r_CH − 1/r_OH − 1/r_CN] · f, with q = 0.42, 0.20 e and 332 kcal/mol·Å conversion).

**References:** McConnell 1957 *J Chem Phys* 27, 226 (kernel form); Kabsch & Sander 1983 *Biopolymers* 22, 2577 (DSSP H-bond definition); Joosten et al. 2011 *Nucleic Acids Res* 39, D411 (libdssp implementation).

---

## Section 5 — Solvent and environment (PRESENT)

### 5.1 `WaterFieldResult` (PRESENT, `src/WaterFieldResult.h`)

**Plain English:** For each protein atom, sums the explicit electric field and EFG from TIP3P water charges (O = -0.834e, H = +0.417e) within cutoff; tracks first-shell (≤ 3.5 Å from O) and second-shell (3.5–5.5 Å) separately.

**Physics:** Water's H-bond network produces a structured field that no continuum dielectric model captures — orientation fluctuations, cavities, bridging waters, structural waters all matter at the atomic scale. The explicit-water field is the per-frame ground truth that APBS approximates with its dielectric continuum.

**References:** Jorgensen et al. 1983 *J Chem Phys* 79, 926 (TIP3P parameters); Buckingham 1960 *Can J Chem* 38, 300 (E-field shielding coupling, applied to explicit water); Sahakyan & Vendruscolo 2013 (¹⁷O sensitivity to E-field up to 80 ppm).

### 5.2 `HydrationShellResult` (PRESENT, `src/HydrationShellResult.h`)

**Plain English:** Per-atom hydration geometry: half-shell asymmetry (UCSB method — fraction of first-shell waters on solvent-exposed side via protein-COM reference), mean water dipole orientation, nearest ion distance and charge.

**Physics:** Local dielectric environment is asymmetric for surface atoms — water cavities differ from hydrophobic pockets. The half-shell asymmetry encodes whether an atom faces water or protein interior; water dipole orientation reflects local H-bond ordering. Both quantities encode dielectric structure invisible to geometry-only kernels.

**References:** Hummer-Garde-Pratt 1996 *J Phys Chem* 100, 1206 (explicit-water statistics framework); 2026 *Chem Rev* — Akram et al. (recommendation of 2 closest waters per polar atom for QM/MM).

### 5.3 `HydrationGeometryResult` (PRESENT, `src/HydrationGeometryResult.h`)

**Plain English:** Replaces HydrationShellResult's protein-COM reference with the atom's actual SASA-derived surface normal, giving cleaner first-shell water dipole alignment, dipole coherence, and half-shell asymmetry per polar atom.

**Physics:** The "solvent-exposed" direction for a buried atom is degenerate; for a surface atom, the SASA outward normal (mean direction of non-occluded Fibonacci probe points) is the correct local reference frame. Water dipole alignment with this normal is a direct measure of local interfacial structuring.

**References:** Shrake & Rupley 1973 *J Mol Biol* 79, 351 (SASA Fibonacci sampling); same hydration framework refs as 5.2.

---

## Section 6 — Per-residue / per-frame structural (PRESENT)

### 6.1 `DsspResult` (PRESENT, `src/DsspResult.h`)

**Plain English:** Per-residue secondary structure assignment (H/G/I/E/B/T/S/C), φ/ψ dihedrals, SASA, and backbone H-bond donor/acceptor partners via libdssp.

**Physics:** DSSP's secondary-structure characters reflect H-bond pattern + backbone geometry: H = 4-helix, G = 3₁₀ helix, I = π-helix, E = β-strand, B = isolated bridge, T = turn, S = bend, C = coil. The H-bond identification uses Kabsch-Sander electrostatic energy < -0.5 kcal/mol.

**References:** Kabsch & Sander 1983 *Biopolymers* 22, 2577 (DSSP); Joosten et al. 2011 *Nucleic Acids Res* 39, D411 (libdssp).

### 6.2 `SasaResult` (PRESENT, `src/SasaResult.h`)

**Plain English:** Per-atom solvent-accessible surface area via Shrake-Rupley with ~92 Fibonacci-lattice probe points and r_probe = 1.4 Å; emits per-atom SASA (Å²) and outward surface normal (mean direction of non-occluded probe points).

**Physics:** SASA is the sphere-of-radius (r_vdW + r_probe) area not occluded by other atoms. The 1.4 Å probe is the canonical water-radius. Outward normal is purely geometric — direction toward the average non-occluded probe point.

**References:** Shrake & Rupley 1973 *J Mol Biol* 79, 351 (Shrake-Rupley); Lee & Richards 1971 *J Mol Biol* 55, 379 (SASA concept).

---

## Section 7 — Energy and dynamics (PRESENT)

### 7.1 `BondedEnergyResult` (PRESENT, `src/BondedEnergyResult.h`)

**Plain English:** Per-atom decomposition of bonded energy terms (bond stretch, angle bend, Urey-Bradley, proper/improper dihedral, CMAP) computed from CHARMM36m / ff14SB parameters and current positions; energy split evenly among participating atoms.

**Physics:** Force-field bonded terms are fitted harmonic / cosine functions:
- Bond: `½ k_b (r − r₀)²`
- Angle: `½ k_θ (θ − θ₀)²`
- Urey-Bradley: `½ k_UB (r_13 − r_UB,0)²`
- Proper dihedral: `Σ_n V_n (1 + cos(nφ − γ_n))`
- Improper: `½ k_ψ (ψ − ψ₀)²`
- CMAP: 2D φ/ψ correction grid (CHARMM36)

**References:** Maier et al. 2015 *J Chem Theory Comput* 11, 3696 (ff14SB); Best et al. 2012 *J Chem Theory Comput* 8, 3257 (CHARMM36m + CMAP); Mackerell et al. 2004 *J Comput Chem* 25, 1400 (CMAP).

### 7.2 `GromacsEnergyResult` (PRESENT, `src/GromacsEnergyResult.h`)

**Plain English:** Reads per-frame whole-system energy components from the GROMACS .edr file (PME short-range + reciprocal, LJ, total potential, temperature, pressure) at the conformation's time and stores as scalars.

**Physics:** GROMACS computes the full simulation energy per step. PME (Particle-Mesh Ewald) splits long-range Coulomb into short-range direct + long-range reciprocal sums for tractability. LJ short-range is the standard Lennard-Jones 12-6 potential. These are aggregate diagnostics, not per-atom.

**References:** Darden, York & Pedersen 1993 *J Chem Phys* 98, 10089 (PME); Lennard-Jones 1924 *Proc R Soc A* 106, 463 (LJ); Abraham et al. 2015 *SoftwareX* 1, 19 (GROMACS modern).

---

## Section 8 — Comparison (PRESENT)

### 8.1 `MutationDeltaResult` (PRESENT, `src/MutationDeltaResult.h`)

**Plain English:** Stores per-atom delta tensors between the WT conformation and a paired mutant conformation: shielding deltas (DFT), APBS field deltas, MOPAC charge deltas, DSSP/graph deltas, plus per-removed-ring cylindrical coordinates (z, ρ, θ, McConnell factor, exponential decay) for the mutation site.

**Physics:** A WT-to-Ala mutation removes the sidechain past Cβ; for aromatic residues (PHE/TYR/TRP/HIS) the removed ring's contributions to neighbour shielding can be computed in isolation. The cylindrical decomposition relative to each removed ring isolates the geometric kernel each ring contributed before mutation; comparison against actual DFT mutation deltas is the calibration target for Stage 1.

**References:** None specific to "mutation deltas as calibration probe" — the methodology is project-internal (`learn/CLAUDE.md`, learn/EXTRACTOR_DESIGN.md). Adjacent: Markwick et al. 2010 *JACS* 132, 1220 (enhanced conformational sampling of IκBα via accelerated MD; chemical shifts as ensemble averages).

---

## Section 9 — Trajectory aggregators (PRESENT)

### 9.1 `BsWelfordTrajectoryResult` (PRESENT, `src/BsWelfordTrajectoryResult.h`)

**Plain English:** Running mean / variance / min / max of BiotSavart shielding T0 (isotropic) and |T2| (magnitude) per atom across all trajectory frames, plus frame-to-frame T0 delta statistics, via Welford's online algorithm.

**Physics:** Welford's algorithm gives numerically stable single-pass running mean/variance: `M_n = M_{n-1} + (x_n - M_{n-1})/n`, `M2_n = M2_{n-1} + (x_n - M_{n-1})(x_n - M_n)`. T0/|T2| moments characterise the per-atom shielding distribution over the trajectory; delta statistics characterise frame-to-frame fluctuation.

**References:** Welford 1962 *Technometrics* 4, 419; Knuth 1997 *The Art of Computer Programming* vol 2.

### 9.2 `BondLengthStatsTrajectoryResult` (PRESENT, `src/BondLengthStatsTrajectoryResult.h`)

**Plain English:** Per-bond mean/std/min/max bond length and frame-to-frame length-delta statistics across the trajectory.

**Physics:** Bond-length distribution width reflects the harmonic potential `½ k_b (r − r₀)²` thermal fluctuation `<Δr²> = k_BT/k_b`. Frame-to-frame deltas at 40 ps cadence sample the integrated bond-vibration motion (vibrations are femtosecond, so deltas at 40 ps reflect already-averaged quasi-equilibrium fluctuation).

**References:** Same as BondedEnergyResult.

### 9.3 `BsAnomalousAtomMarkerTrajectoryResult` (PRESENT, `src/BsAnomalousAtomMarkerTrajectoryResult.h`)

**Plain English:** Per-atom per-frame outlier marker against the BiotSavart T0 distribution: emits `BsAnomalyHighT0` records when z-score > +threshold and `BsAnomalyLowT0` when z < -threshold; pushes to the per-atom event bag.

**Physics:** Anomalies in the per-frame ring-current shielding (relative to the trajectory's running mean and variance) flag specific frames where the atom transiently entered an unusual geometric configuration relative to nearby aromatic rings — e.g., approaching a stacking event or a transient close contact. The z-score threshold is configurable.

**References:** Standard outlier-detection methodology; not a physics paper. The framing as "events" relates to extreme-value statistics for trajectory analysis.

### 9.4 `BsShieldingTimeSeriesTrajectoryResult` (PRESENT, `src/BsShieldingTimeSeriesTrajectoryResult.h`)

**Plain English:** Dense per-atom per-frame buffer of the BiotSavart shielding contribution as a SphericalTensor; emitted at end of trajectory as `(N, T, 9)` float64 with e3nn-compatible irrep layout (`0e+1o+2e` for magnetic kernels).

**Physics:** Storing the full per-frame T0+T1+T2 trajectory enables downstream temporal analyses (autocorrelation, spectral density, S²) without re-running the kernel calculator. Parity is `0e+1o+2e` because B is an axial pseudovector in e3nn convention; Coulomb / APBS EFG (which are symmetric traceless) use `0e+1e+2e`.

**References:** Geiger & Smidt 2022 *arXiv* 2207.09453 (e3nn library); Kondor 2025 (irreducible decompositions of tensorial physical observables).

### 9.5 `BsT0AutocorrelationTrajectoryResult` (PRESENT, `src/BsT0AutocorrelationTrajectoryResult.h`)

**Plain English:** Per-atom normalised autocorrelation function `ρ(k) = C(k)/C(0)` of the BiotSavart T0 shielding contribution at lags k = 0, 1, …, N_LAGS-1, using the biased estimator that guarantees |ρ(k)| ≤ 1.

**Physics:** Autocorrelation `C(k) = (1/N) Σ_t (x_t − ⟨x⟩)(x_{t+k} − ⟨x⟩)` is the standard time-domain measure of trajectory persistence. Biased estimator (divide by N rather than N-k) guarantees positive-definiteness; the unbiased one drifts outside [-1, 1] when trajectory mean differs from sub-window mean.

**References:** Berne & Harp 1970 *Adv Chem Phys* 17, 63 (TCF in MD); Allen & Tildesley 2017 *Computer Simulation of Liquids* OUP, ch. 8 (autocorrelation estimators).

### 9.6 `ChiRotamerSelectionTrajectoryResult` (PRESENT, `src/ChiRotamerSelectionTrajectoryResult.h`)

**Plain English:** Per-residue chi-angle rotamer transition detector: classifies each frame's χ_k into rotamer bins (g+, t, g−) and pushes a SelectionRecord onto the run-scope bag whenever the bin changes from the previous frame.

**Physics:** Sidechain χ angles populate ~120°-spaced rotamer wells separated by ~3 kcal/mol barriers; transitions happen on 100 ps – 10 ns timescales for common residues. Detecting transitions per residue per frame enables event-based analyses (e.g., rotamer-population reweighting, dwell-time distributions).

**References:** Lovell et al. 2000 *Proteins* 40, 389 (rotamer libraries); Ponder & Richards 1987 *J Mol Biol* 193, 775 (rotamer concept).

### 9.7 `PositionsTimeSeriesTrajectoryResult` (PRESENT, `src/PositionsTimeSeriesTrajectoryResult.h`)

**Plain English:** Dense per-atom per-frame position buffer; emitted as `(N, T, 3)` float64. Worked example for Vec3-typed time-series TR pattern.

**Physics:** Just stored coordinates. No new physics; enables downstream time-series analyses without re-reading XTC.

**References:** None.

---

## Section 10 — Calculator-pass IMPROVEMENTS to existing kernels (literature-driven, planned)

These are not new kernels; they are corrections / refinements to existing calculators motivated by specific literature findings.

### 10.1 Ring-normal stability fix (IMPROVEMENT to BiotSavartResult + HaighMallionResult)

**Plain English:** Replace the single 3-atom ring normal with the average of two independently computed normals (e.g., atoms 1-2-3 and atoms 4-5-6 of a six-membered ring), which is more robust to MD out-of-plane fluctuations.

**Physics:** A 6-membered aromatic ring sampled from MD has thermal out-of-plane atomic displacements of ~0.05–0.1 Å. The 3-atom-derived normal can swing ~5° on these displacements, propagating directly into ring-current tensor angular error. Two-normal averaging cancels first-order out-of-plane noise.

**References:** Sahakyan & Vendruscolo 2013 *J Phys Chem B* 117, 1989 — explicit recommendation (`references-text/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-text-{1-4}.txt`).

### 10.2 Single-loop default audit (IMPROVEMENT to BiotSavartResult)

**Plain English:** Audit the loop separation parameter (default S in the Johnson-Bovey two-loop model) against the Agarwal 1977 [10]-paracyclophane benchmark; single-loop (S=0) fits experiment with r=0.9928 vs two-loop (S=1.28 Å) with r=0.8883.

**Physics:** Johnson-Bovey originally proposed two loops at ±0.64 Å above/below the ring plane to model the spatial extent of the π density. Agarwal's [10]-paracyclophane geometry probes the ring at known positions; their experimental NMR shifts test which loop separation gives the right field. Single-loop gives a much higher fit quality for this benchmark.

**References:** Agarwal et al. 1977 *Can J Chem* 55, 2575 (`references-text/agarwal-1977-ring-currents-local-anisotropy-paracyclophane-text-{1-3}.txt`); Johnson & Bovey 1958 *J Chem Phys* 29, 1012.

### 10.3 H-bond geometry: angle θ over distance d (IMPROVEMENT or NEW per-atom feature)

**Plain English:** When characterising H-bonds for kernel input, use the donor-H-acceptor angle θ as the primary geometric feature; the distance d is much weaker as a CS predictor.

**Physics:** In MD trajectories, hydrogen-bond *angle* θ couples strongly to backbone ψ, which itself drives ¹⁵N chemical shifts; the full ψ excursion produces ~25 ppm shifts in flexible residues, ~5 ppm per 5° of ψ change in less-flexible residues. H-bond *distance* d fluctuations correlate only weakly with shifts. The dipolar shielding contribution from an H-bond partner scales as `(3 cos²θ - 1)/r³`; angular dependence dominates over radial when r is roughly fixed.

**References:** Yi-McDermott 2024 *J Phys Chem Lett* 15, 2270 (`references-text/yi-mcdermott-2024-temperature-shifts-conformational-dynamics-text-{1-3}.txt`).

### 10.4 Buckingham σ^(1)·E coefficient verification (validation, not new kernel)

**Plain English:** Validate that our existing CoulombResult / WaterFieldResult / ApbsFieldResult E-fields produce the expected isotropic shielding contribution per Buckingham's analytical coefficient: Δσ ≈ −2 × 10⁻¹² E_z − 10⁻¹⁸ E². At E ~ 7 Å from a unit charge, the linear term is ~0.2 ppm — about 20× the quadratic term.

**Physics:** Buckingham 1960 derived `σ = σ_0 + σ_b + σ_a + σ_w + σ_E` for solvent shielding; the σ_E term is the linear E-field response with coefficient −3.0 ppm/au (nucleic acids, Case 1995) or −3.4 ppm/au (proteins). At known geometries, the field-derived shielding should match this analytical coefficient.

**References:** Buckingham 1960 *Can J Chem* 38, 300 (`references-text/buckingham-1960-chemical-shifts-polar-groups-text-{1-3}.txt`); Case 1995 *J Biomol NMR* 6, 341.

### 10.5 Per-element CSA validation tests (test additions, not new kernel)

**Plain English:** Add unit-test cases that assert the rank-2 tensor magnitude reaches literature values at known geometries: 16.6 ppm CSA for an N–H···π hydrogen bond from ring-current alone; up to 80 ppm for ¹⁷O via electric field. Catches T2 normalisation or sign bugs.

**Physics:** Boyd-Skrynnikov 2002 published the closed-form σ_xz formula extending Johnson-Bovey σ_zz to the full rank-2 tensor; at G42 of fibronectin type-2, ring-current-dominant N–H···π bond gives σ_rc CSA contribution of 16.6 ppm. Sahakyan-Vendruscolo 2013 reports up to 80 ppm ¹⁷O sensitivity via E-field for RNA bases.

**References:** Boyd & Skrynnikov 2002 *JACS* 124, 1832; Sahakyan & Vendruscolo 2013 *J Phys Chem B* 117, 1989.

### 10.6 EF dominance audit for ¹³C / ¹⁵N / ¹⁷O (audit, not new kernel)

**Plain English:** Verify that our calibration weighting reflects the literature finding that electric-field effects dominate over ring-current for heavy nuclei: for RNA bases, ¹³C ring-current R = 0.257 vs EF R = 0.702.

**Physics:** Buckingham's σ^(1)·E coefficient |γ_C| ~ 1 ppm/au is much larger than the ring-current-induced shift for typical geometries; for ¹³C, ¹⁵N, ¹⁷O the EF sensitivity dominates, while for ¹H ring-current dominates. This is reflected in our existing `learn/docs/element_physics_2026-04-10.md` per-element R² table.

**References:** Sahakyan & Vendruscolo 2013 *J Phys Chem B* 117, 1989; learn/docs/element_physics_2026-04-10.md (project-internal validation).

---

## Section 11 — NEW geometric kernels (planned, literature-driven)

### 11.1 `Ch3ShiftMethylKernelResult` (PLANNED)

**Plain English:** Methyl-specific kernel decomposing CH₃ proton shielding into ring-current rotamer-averaged + dihedral correction + ring contribution + magnetic anisotropy + electric field + phenomenological distance terms; mirrors Sahakyan 2011's CH3Shift formula.

**Physics:** `δ = δ_rc^rot + Δδ_dih + Δδ_ring + Δδ_ma + Δδ_EF + Δδ_dist`. Each term is a polynomial function of interatomic distances within a 6.5 Å sphere (active region) plus a 1.8 Å neutral region around the methyl proper. Differentiable so usable as restraints. RMSD reaches 0.133–0.198 ppm on Ala/Thr/Val/Leu/Ile methyls.

**References:** Sahakyan, Vranken, Cavalli & Vendruscolo 2011 *J Biomol NMR* 50, 331 (`references-text/sahakyan-2011-methyl-chemical-shifts-proteins-text-{1-6}.txt`).

### 11.2 `LarsenProcs15CorrectionsResult` (PLANNED, multi-term)

**Plain English:** Implements the full ProCS15 long-range correction set as named additions to the tripeptide-baseline shielding: Δσ_BB^(i±1) (next-nearest backbone), Δσ_HB (H-bond), Δσ_HαB (Hα bonding), Δσ_sc (sidechain effects).

**Physics:** Larsen 2015 ProCS15 trains 2.35M OPBE/6-31G(d,p)//PM6 calculations on tripeptides + H-bonding model systems, decomposing per-atom shielding into σ_BB (the tripeptide DFT lookup) plus four named correction terms. Each correction encodes a specific local-to-medium-range structural effect; together they reach RMSD 2.2 ppm ¹³C, 0.7 ppm ¹H, 4.8 ppm ¹⁵N on protein structures.

**References:** Larsen, Bratholm, Christensen, Channir & Jensen 2015 *PeerJ* 3:e1344 (`references-text/larsen-2015-procs15-dft-chemical-shift-predictor-text-{1-7}.txt`).

### 11.3 `PelloniDifferentialBiotSavartResult` (PLANNED, exploratory)

**Plain English:** Per-volume-element formulation of the Biot-Savart ring-current shielding: re-formulates the wire integral as a per-volume-element shielding density Σ_αβ(r), enabling near-vs-far decomposition (near segment dominates and deshields for benzene ¹H; far segment dominates and shields for ¹³C).

**Physics:** The Biot-Savart line integral can be re-cast as an integral over volume elements of the current density J(r), giving a per-element shielding contribution. Integrating in spatial regions reveals which parts of the ring current produce shielding vs deshielding for a specific probe nucleus geometry. Pelloni 2004 reports a near/far decomposition for benzene; the π contribution at the near segment is approximately +19 ppm for σ_zz^C, with the corresponding deshielding contribution at σ_zz^H of comparable magnitude (sign-flipped per the local-vs-far geometry).

**References:** Pelloni, Ligabue & Lazzeretti 2004 *Org Lett* 6, 4451; Lazzeretti 2004 *Phys Chem Chem Phys* 6, 217 (current density formalism).

### 11.4 `MoynaRingCurrentComparisonResult` (PLANNED, comparative)

**Plain English:** Ensemble comparison of three ring-current methods (HM, JB, Pople-PD-new) on per-frame data; reports per-method R² against DFT, RMS gap, and reparameterised constants.

**Physics:** Moyna 1998 directly compared HM (r=0.854), JB (r=0.846), Pople point-dipole (PD-new, r=0.846) on 11406 protein ¹H shifts, with ~5% RMS gap between methods. Reparameterisation of the PD-new constant against HM-fitted gives B = 21.82 ppm (vs Perkins-Dwek's 27.41).

**References:** Moyna, Zauhar, Williams, Nachman & Scott 1998 *J Chem Inf Comput Sci* 38, 702.

### 11.5 `ChargeDifferentialKernelResult` (PLANNED, generic)

**Plain English:** Per-atom Δkernel(frame) = kernel(frame) − kernel(reference) for each kernel-output T2 tensor; reference is configurable (minimum-energy frame from MD, OF3-predicted "exact-to-bottle" structure, trajectory-mean kernel, or first post-equilibration frame).

**Physics:** Static-vs-dynamic decomposition: the reference state's kernel value is the "static" prediction; deviations across the trajectory are the "dynamic" contribution. Each kernel improvement (calculator-pass items 10.1–10.6) reduces residual magnitude in the differential, providing a direct quantitative "show some effect" measure.

**References:** Robustelli, Stafford & Palmer III 2012 *JACS* 134, 6365 — interpretation pivot for MD-vs-X-ray shielding deltas (`references-text/robustelli-stafford-palmer-2012-protein-dynamics-shifts-text-{1-4}.txt`); Li & Brüschweiler 2012 *J Biomol NMR* 54, 257 (static a^MD vs averaged a^Xray distinction).

### 11.6 `BackboneNHUnitVectorResult` (PLANNED, trivial)

**Plain English:** Per-residue per-frame backbone N–H unit vector in the protein-aligned frame (after backbone-superimposed rotation removal); feeds Lipari-Szabo S² and direct R1/R2/NOE back-calculation.

**Physics:** N–H bond vector autocorrelation `C(τ) = ⟨P₂(μ̂(0)·μ̂(τ))⟩` (with P₂ the second-order Legendre) is the foundational time series for NMR relaxation observables. Removing overall protein tumbling via backbone alignment gives the internal-motion-only autocorrelation that Lipari-Szabo S² fits.

**References:** Lipari & Szabo 1982 *JACS* 104, 4546 (model-free); Lesovoy & Orekhov 2025 *Int J Mol Sci* 26, 8917.

### 11.7 `HydrationShellPolarAtomResult` (PLANNED, refinement of HydrationShellResult)

**Plain English:** Per-polar-atom 2-closest-water selector and tracker: stores indices, distances, donor-H-acceptor angles, and exchange events for the 2 closest waters within ~3.5 Å of each protein polar atom across the trajectory.

**Physics:** The 2026 Chemical Reviews recommendation (Akram et al.) is that explicit-solvent CS calculations need only the two closest water molecules combined with PCM to recover most of the solvent shift. Tracking these waters per polar atom per frame gives the explicit-water analog of APBS, with exchange events flagged.

**References:** Akram et al. 2026 *Chem Rev* (recommendation); fragment-QM/MM literature (3.5 Å cutoff convention).

---

## Section 12 — Trajectory aggregators (PLANNED)

### 12.1 `BlockAveragedConvergenceResult` (PLANNED, gating diagnostic)

**Plain English:** Per-atom per-tensor-component block-averaged standard error of the mean as a function of block length; directly answers "is the arithmetic-mean shielding tensor over our trajectory converged?"

**Physics:** Grossfield-Zuckerman 2009 block averaging: divide the trajectory into B contiguous blocks; compute per-block means; the SEM at block size T is `σ(block_means)/√B`. As block size grows past the integrated autocorrelation time, SEM stops decreasing, giving the apparent independent-sample count `N_eff = (T_total / SEM(T_max)²) · σ²_full`.

**References:** Grossfield & Zuckerman 2009 *Annu Rep Comput Chem* 5, 23.

### 12.2 `SigmaLipariSzaboResult` (PLANNED)

**Plain English:** Per-atom Lipari-Szabo fit `C(t) = S² + (1−S²)·exp(−t/τ_e)` of the shielding-tensor (or N–H bond vector) autocorrelation; emits (S², τ_e) per atom per irrep component.

**Physics:** Lipari-Szabo "model-free" approach assumes internal motion is fast and uncorrelated with overall tumbling; S² is the squared order parameter (long-time plateau of C(t)) and τ_e is the effective internal-motion timescale. For backbone N–H, S² values 0.85-0.9 indicate rigid; 0.4-0.6 indicate flexible loops.

**References:** Lipari & Szabo 1982 *JACS* 104, 4546 (model-free); Clore et al. 1990 *JACS* 112, 4989 (extended LS-2 for bimodal motions); Kasinath & Wand 2013 *Chem Rev* 113, 9173 (S² as entropy meter).

### 12.3 `SigmaEssentialDynamicsResult` (PLANNED, exploratory)

**Plain English:** PCA / essential dynamics decomposition of per-atom shielding-tensor trajectory: stack centered (T, N·9) tensor matrix, SVD across frames, extract top-k modes (spatial pattern + time-course) ordered by variance explained.

**Physics:** Essential dynamics extracts collective modes from time-indexed protein data. Applied to coordinates, it recovers backbone collective motions; applied to shielding tensors, it recovers collective shielding-fluctuation modes. Cumulative variance curve diagnoses whether the trajectory is low-rank (a few modes carry most variance) or high-dimensional.

**References:** Amadei, Linssen & Berendsen 1993 *Proteins* 17, 412.

### 12.4 `SigmaTuckerDecompositionResult` (PLANNED, exploratory)

**Plain English:** Tucker decomposition (higher-order SVD) of the (T, N, K, 9) shielding tensor trajectory (T = frames, N = atoms, K = kernels, 9 = irrep components); produces a core tensor plus mode matrices along each axis.

**Physics:** Tucker decomposition extends SVD to higher-order tensors; the core tensor encodes joint factor interactions, mode matrices encode per-axis patterns. Applied to multi-way time series (frames × atoms × kernels × tensor-components), it can reveal patterns that 2D PCA averages out — e.g., "mode 3 is ring-current-dominated aromatic shielding fluctuating with backbone φ/ψ flips on timescale τ_3."

**References:** Kolda & Bader 2009 *SIAM Rev* 51, 455; Tucker 1966 *Psychometrika* 31, 279 (original).

### 12.5 `GreenKuboSpectralDensityResult` (PLANNED)

**Plain English:** Per-atom spectral density `J(ω)` at named Larmor frequencies (typical 600 / 700 / 850 MHz) computed as the Fourier transform of the per-atom σ(t) autocorrelation; direct input to R1 / R2 / R1ρ relaxation predictions.

**Physics:** By Wiener-Khinchin, the Fourier transform of an autocorrelation function is the spectral density. For NMR relaxation, R1 ∝ J(ω₀) + ··· at the nuclear Larmor frequency; J(ω) at multiple frequencies sets the full relaxation tensor. The shielding-tensor autocorrelation gives the CSA contribution to relaxation specifically.

**References:** Kubo 1957 *J Phys Soc Jpn* 12, 570 (linear response); Berne & Harp 1970 *Adv Chem Phys* 17, 63 (TCF in MD); Tugarinov & Kay 2003 *J Am Chem Soc* 125, 13868 (CCR rates).

### 12.6 `KernelTimeSeriesMomentsResult` (PLANNED, simple)

**Plain English:** Per-atom per-kernel mean, variance, skewness, kurtosis, and selected percentiles of each kernel time-series across the trajectory.

**Physics:** Higher-order moments characterise the shape of the kernel-value distribution per atom: mean is the static prediction, variance is fluctuation magnitude, skewness flags asymmetric distributions (e.g., bimodal stack/unstack ring states), kurtosis flags heavy-tailed events (rare large excursions).

**References:** Standard distributional moment definitions; not a physics paper.

### 12.7 `DihedralAutocorrelationResult` (PLANNED)

**Plain English:** Per-residue per-dihedral (φ, ψ, χ_k) autocorrelation function across the trajectory; supports rotamer-transition-rate estimates and direct comparison with NMR-derived dihedral order parameters.

**Physics:** Dihedral autocorrelation `C_φ(τ) = ⟨cos(φ(0)−φ(τ))⟩` gives the dihedral relaxation timescale per residue. For sidechain χ angles, decay timescales of 100 ps – 10 ns are typical, reflecting rotamer-well dwell times. For backbone φ/ψ in folded regions, decay is mostly within picosecond librations.

**References:** Brüschweiler 1998 *Curr Opin Struct Biol* 8, 168 (NMR dihedral order parameters); Lai & Brooks III 2024 *J Phys Chem B* 128, 10813 (S² / dihedral coupling).

### 12.8 `DihedralDistributionKDEResult` (PLANNED, Gadanecz-style)

**Plain English:** Per-residue kernel density estimate of the (φ, ψ) joint distribution across the trajectory; emit as (N_residue × N_grid_φ × N_grid_ψ) array.

**Physics:** Per-residue (φ, ψ) populates discrete Ramachandran basins; the trajectory KDE captures the basin populations + transitions. Gadanecz 2026 EDMD uses the negative log of the KDE as a Boltzmann-inversion potential for restraint-driven structure refinement.

**References:** Gadanecz, Fazekas, Menyhárd & Perczel 2026 *J Chem Inf Model* 66, 2844 (EDMD); Ramachandran et al. 1963 *J Mol Biol* 7, 95 (Ramachandran plot original).

### 12.9 `HBondGeometryDistributionResult` (PLANNED)

**Plain English:** Per-amide H-bond donor-acceptor distance d, donor-H-acceptor angle θ, and partner-residue identity time series across the trajectory; emit per-amide histograms + transition events.

**Physics:** Backbone H-bond fluctuations (geometry + partner identity) are the primary source of backbone amide ¹⁵N shielding variance per Yi-McDermott 2024. Per-amide distribution of (d, θ) over the trajectory plus partner-residue tracking gives the data to compare against published amide ¹⁵N relaxation observables.

**References:** Yi-McDermott 2024 *J Phys Chem Lett* 15, 2270.

---

## Section 13 — External-tool wrappers (PLANNED)

### 13.1 `LegolasShiftPredictorResult` (PLANNED)

**Plain English:** Per-frame backbone chemical-shift prediction via LEGOLAS (TorchANI-based DNN trained on DFT shifts); outputs HA, HN, CA, CB, C′, N scalar shifts per residue per frame as an extra feature channel for the model.

**Physics:** LEGOLAS is a TorchANI message-passing neural network trained on a large set of DFT-computed protein backbone shifts. It learns to predict CS from atomic environments; integrating its output as a per-frame feature gives the GNN a "DFT-distilled prior" channel alongside the physics-decomposed kernel features.

**References:** Darrows, Kodituwakku, Xue, Pickering, Terrel & Roitberg 2025 *J Chem Theory Comput* 21, 4266; Smith et al. 2017 *Chem Sci* 8, 3192 (TorchANI).

### 13.2 `UCBShift2Result` (PLANNED)

**Plain English:** Per-frame backbone + sidechain CS prediction via UCBShift 2.0 (transfer-learning ML predictor); outputs scalar shifts as another comparison channel.

**Physics:** UCBShift uses a transfer-prediction approach combining sequence- and structure-based ML; UCBShift 2.0 extends to side-chain predictions. As with LEGOLAS, integration is as a feature channel rather than a comparison tool.

**References:** Ptaszek et al. 2024 *JACS* 146, 31733.

---

## Section 14 — Schema-without-producer ORPHANS (treat as PLANNED per user direction 2026-04-30)

### 14.1 `RingExponentialSumResult` (ORPHAN — schema slot; user wants implemented)

**Plain English:** Per-atom exponential-distance-weighted sum of ring contributions: `G_iso = Σ_rings exp(-r_ring/scale) · g_iso(ring, atom)`, similar for the T2 tensor and an 8 Å variance variant.

**Physics:** Distance-weighted ring-influence aggregation: each ring contributes with weight exp(-r/scale), bounding the contribution at long range and giving a smooth proxy for "effective ring proximity." T2 version weights the rank-2 tensor; the 8 Å variance variant captures local-cluster-level variability. Schema slots `G_iso_exp_sum`, `G_T2_exp_sum`, `G_iso_var_8A` were declared in `ConformationAtom` but no calculator currently writes them.

**References:** Reference TBD — the form is consistent with general distance-weighted aggregation; the specific scale parameter and weight choice need literature anchor.

### 14.2 `RingDistanceStatsResult` (ORPHAN — schema slot)

**Plain English:** Per-atom mean ring-center distance and nearest-ring-atom distance across all rings within the cutoff; simple geometric scalars.

**Physics:** Two scalars: `mean_ring_distance = ⟨r_ring⟩` over rings within cutoff; `nearest_ring_atom_distance = min_atoms_in_rings d`. No new physics; geometric summary statistics. Schema slots declared but no producer.

**References:** None required (geometry).

---

## Section 15 — References by year (cross-reference table)

```
1957  McConnell — bond magnetic anisotropy formalism                         §4.4
1958  Johnson & Bovey — ring current double-loop wire model                   §4.1
1960  Buckingham — σ^(1)·E linear electric-field shielding                    §4.6, 10.4
1962  Welford — running variance algorithm                                     §9.1
1963  Ramachandran — φ/ψ dihedral plot                                         §12.8
1964  Klopman / Ohno — Ohno-Klopman γ kernel                                   §2.2
1966  Tucker — multi-way decomposition                                         §12.4
1968  Wiberg — bond order index                                                §4.5
1970  Berne & Harp — TCF in MD                                                 §9.5, 12.5
1971  Lee & Richards — SASA concept                                            §6.2
1972  Haigh & Mallion — surface-integral ring current                          §4.2
1973  Shrake & Rupley — Fibonacci-lattice SASA                                 §6.2, 5.3
1977  Agarwal — [10]-paracyclophane benchmark                                  §10.2
1979  Haigh & Mallion — full HM model                                          §4.2
1982  Lipari & Szabo — model-free relaxation                                   §11.6, 12.2
1983  Kabsch & Sander — DSSP H-bond definition                                 §4.10, 6.1
1983  Brooks et al. — CHARMM switching function                                §4.9
1983  Jorgensen et al. — TIP3P water                                           §5.1
1990  Clore et al. — extended Lipari-Szabo (LS-2)                              §12.2
1993  Amadei et al. — essential dynamics                                       §12.3
1993  Bayly et al. — RESP charges                                              §2.1
1993  Darden, York, Pedersen — PME                                             §7.2
1995  Case — ring current calibration vs DFT                                   §4.1, 10.4
1995  Schreckenbach & Ziegler — DFT shielding implementation                   §3.2
1996  Hummer-Garde-Pratt — explicit-water statistics                           §5.2
1998  Moyna et al. — HM/JB/PD ring-current comparison                         §11.4
1998  Markley — IUPAC NMR atom names                                           §1.6
1999  Helgaker et al. — DFT NMR review                                         §3.2
2001  Baker et al. — APBS                                                      §3.1
2002  Boyd & Skrynnikov — full ring-current tensor (σ_xz formula)              §4.1, 10.5
2003  Tugarinov & Kay — CCR rates                                              §12.5
2004  Pelloni, Ligabue, Lazzeretti — differential Biot-Savart                 §11.3
2007  Behler & Parrinello — NN potential foundation                            §2.4
2019  Caldeweyher et al. — D4 / EEQ                                            §2.2
2009  Grossfield & Zuckerman — block averaging                                 §12.1
2009  Kolda & Bader — tensor decompositions review                            §12.4
2010  Vanommeslaeghe et al. — CHARMM CGenFF                                    §2.1
2010  Markwick et al. — AMD on IκBα                                            §8.1
2011  Sahakyan et al. — CH3Shift                                              §11.1
2011  Joosten et al. — libdssp                                                 §4.10, 6.1
2011  Han et al. — SHIFTX2 (¹⁵N ring-current isotropic <0.6% bound)            §18.4
2012  Robustelli, Stafford, Palmer — MD-averaged CS interpretation             §11.5
2012  Li & Brüschweiler — PPM predictor / static a^MD                          §11.5
2012  Best et al. — CHARMM36m                                                  §7.1
2012  Neese — ORCA                                                             §3.2
2013  Sahakyan & Vendruscolo — RNA EFG dominance + ring-normal stability      §10.1, 10.5, 10.6
2013  Stewart — PM7 / MOZYME                                                   §2.3, 4.5
2013  Stone — Theory of Intermolecular Forces (T-tensor)                       §4.6, 4.8
2013  Kasinath & Wand — S² entropy meter                                       §12.2
2015  Maier et al. — ff14SB                                                    §2.1, 7.1
2015  Larsen et al. — ProCS15                                                  §11.2
2015  Abraham et al. — GROMACS modern                                          §7.2
2017  Smith et al. — TorchANI                                                  §13.1
2017  Allen & Tildesley — Computer Simulation of Liquids (autocorrelation)     §9.5
2022  Geiger & Smidt — e3nn library                                            §9.4
2024  Anstine et al. — AIMNet2                                                 §2.4
2024  Ptaszek et al. — UCBShift 2.0                                            §13.2
2024  Lai & Brooks III — S² convergence for relaxation                         §12.7
2024  Yi-McDermott — H-bond angle θ over distance                             §10.3, 12.9
2025  Darrows et al. — LEGOLAS                                                §13.1
2025  Lesovoy & Orekhov — N-H autocorrelation R1/R2/NOE                        §11.6
2025  Kondor — irreducible decompositions                                      §9.4
2026  Gadanecz et al. — EDMD                                                  §12.8
2026  Akram et al. — Computational NMR review (2-water rec)                    §11.7
```

---

## Section 16 — What this document is NOT

- Not a sequencing plan. The order in which calculators land is a separate Chunk 1 question.
- Not a test plan. Per-calculator validation tests (Boyd-Skrynnikov geometric benchmarks, Agarwal paracyclophane, Buckingham σ^(1) coefficient checks) need their own doc.
- Not a calibration architecture. The "static a^MD vs averaged a^Xray" discipline (Li-Brüschweiler 2012) is a methodology decision that lives in `learn/CLAUDE.md` and the ridge-regression code.
- Not a complete pipeline. The Chunk 2 OF3 + tripeptide assembler is documented elsewhere; this doc lists calculators inside nmr-extract, not the model layer or the structure-prep layer.

## Section 17 — Reviewer pass: clean, truthful description for every calculator

**The goal.** A clean, truthful description of the purpose and functionality of every calculator listed — both PRESENT (already in `src/`) and PLANNED (we want to add). Truthful in two senses:

- **Implementation truth.** For PRESENT items, does the description match what the code actually does? For PLANNED items, does the proposed implementation shape make sense given the dependencies and the existing architecture?
- **Concept truth.** Is the physics / math statement correct on its own terms? Are the references cited correctly? Does the named mechanism actually do what the description says?

This applies equally to PRESENT and PLANNED items. Every calculator stays in the inventory. The reviewer's job is to confirm that what we say about each is right — not to recommend dropping, consolidating, or triaging any of them.

**For each item, please check:**

1. **Purpose (plain English).** Does the description accurately say what this calculator does (PRESENT) or what we want it to do (PLANNED)? Propose a more accurate phrasing if the current one is vague or wrong about purpose.

2. **Physics / concept truth.** Is the physics one-liner mathematically and physically correct?
   - Equations: signs, factors, units, tensor indices
   - Mechanism: the named mechanism (ring current, EFG, dipolar bond anisotropy, dispersion, π-quadrupole, etc.) — accurately described?
   - Sign conventions: BiotSavart's `-n_b · B_a · PPM_FACTOR`; traceless projection of Coulomb EFG; McConnell's full M_ab vs symmetric-traceless K_ab — distinctions stated correctly?
   - For PLANNED items: is the proposed derivation consistent with the cited literature?

3. **Implementation truth.** For PRESENT items: does the named source file (`src/<Name>Result.h`) exist; does its actual content / docstring agree with our description here? For PLANNED items: does the proposed implementation shape (inputs, kernel form, dependencies, output type) make sense as a calculator we could build?

4. **References.** For each citation:
   - Year, journal volume, page correct?
   - Authors right?
   - Does the cited paper actually contain the claim attributed to it?
   - For local-cache references, is the filename accurate?

5. **Anything missing.** If a calculator exists in `src/` (under any name, not just `*Result.h`) and contributes physics that isn't here, flag for inclusion. If a literature finding from the cache surfaced a named calculator that isn't in §10-14, flag for inclusion. Adding entries is in scope.

**On Section 18 (atom-type dimensional contribution).** The user is aware that the Stage 1 analysis backing §18 has known limitations — the granularity gap (type-level not per-AMBER-name) is one example, but there are others, and the Stage 1 plumbing is something the user already plans to clean up. **Note any limitations briefly, but do not spend significant attention there.** The §18 numbers are reproductions from existing analysis; they should be taken as-is for this pass. What matters for the reviewer is whether the §18.4 narrative accurately reflects what the §18.2 / §18.3 tables show — not whether the underlying analysis is the strongest possible.

**Out of scope:**

- Recommending we drop items
- Recommending we consolidate items
- Calling items "redundant" or "duplicate" — multiple ways of accessing the same physical dimension is exactly what §18.1 documents and is part of the design
- Gating on "is this useful for the thesis" or "is this worth the effort"
- Proposing a different inventory shape
- Deep critique of the Stage 1 analysis underlying §18

The shape is fixed: every calculator has a place. The reviewer's question is whether the description of each one is **true and clear enough**.

---

## Section 18 — Atom-type dimensional contribution

Source: `learn/stage1-mutations/notes/atom_type_stratification.md`, `learn/stage1-mutations/notes/dimension_inventory_raw.md`, `learn/stage1-mutations/notes/master_chart.md`, `learn/docs/element_physics_2026-04-10.md`. All numbers are 720 proteins × 446K atoms; ridge regression with per-protein kernel normalisation and mutation-type categorical scalars (the "fair" R²).

### 18.1 The actual physical dimensions accessed

Per `dimension_inventory_raw.md`, the 121 T2 kernels in the extractor are not 121 independent dimensions. They are multiple ways of accessing **three primary physical dimensions**, plus secondary dimensions from specific physics:

```
DIMENSION                      PHYSICS                         CALCULATORS THAT ACCESS IT

1. Ring current B-field        π-electron circulation →        BiotSavartResult
                                magnetic dipole field at        HaighMallionResult
                                nearby nuclei. Distance r⁻³.    RingSusceptibilityResult
                                T2 angular pattern              (3 calculators, same physics,
                                (1 − 3cos²θ) dipolar.            different approximations)

2. Electric field gradient      Partial-charge sum at every     CoulombResult (ff14SB charges)
                                atom; pure T2 (traceless        MopacCoulombResult (PM7 charges)
                                symmetric by Gauss's law).      AIMNet2Result EFG path
                                Distance r⁻³.                    ApbsFieldResult (Poisson-
                                                                 Boltzmann, solvent-aware)

3. Bond magnetic anisotropy     Each covalent bond carries     McConnellResult
                                Δχ susceptibility anisotropy;  MopacMcConnellResult
                                dipolar field at nearby atoms;  HBondResult
                                full T0+T1+T2 tensor.            (per-bond / per-H-bond)
                                Distance r⁻³.

SECONDARY DIMENSIONS

π-electron quadrupole           T2 from quadrupolar π-cloud,   PiQuadrupoleResult
                                pure T2 with r⁻⁵ leading
                                radial decay.

Dispersion                      Van der Waals from ring         DispersionResult
                                vertices; r⁻⁶ scalar +
                                anisotropic r⁻⁸ tensor.

Charge-polarisation sub-dim     The gap between fixed (ff14SB)  Manifests as MopacCoulomb
                                and responsive (MOPAC, AIMNet2) − Coulomb difference
                                charges.

Ring-type-specificity           PHE / TYR / TRP / HIE differ    Mutation-type categorical
                                in current intensity, charge    (a scalar input to the
                                distribution, susceptibility    ridge regression, not a
                                anisotropy.                      kernel)

Solvent                         Solvent dielectric, ions,       ApbsFieldResult
                                explicit waters.                 WaterFieldResult
                                                                 HydrationShellResult
                                                                 HydrationGeometryResult
```

### 18.2 Per-element R² by kernel group (720 proteins, fair R²)

```
ELEMENT  n           ring_current   efg     bond_aniso   BS    HM    EFG_aro  ALL
H        230,135     0.909          0.770   0.139        0.772 0.784 0.766    0.935
C        133,488     0.156          0.431   0.041        0.097 0.103 0.390    0.546
N        39,954      0.325          0.100   0.070        0.140 0.133 0.079    0.434
O        42,429      0.317          0.196   0.055        0.169 0.169 0.155    0.387
all      446,006     0.213          0.318   0.027        0.171 0.174 0.302    0.385
```

(EFG_aro = single MOPAC EFG aromatic kernel, the standout single feature.)

### 18.3 Per-atom-type R² (the actual structure: backbone vs sidechain stratification)

```
ATOM TYPE      n          RAW R²    +SCALES   +MUT      FAIR R²

H              230,135    0.921     —         —         0.928

C              133,488    0.514     0.529     0.512     0.562
   CA          29,944     0.541     0.577     0.597     0.627
   C=O         29,944     0.361     0.411     0.430     0.463
   CB          27,429     0.578     0.597     0.604     0.647
   C side      46,171     0.660     0.690     0.694     0.729  ← strong

N              39,954     0.210     0.292     0.345     0.380
   N bb        29,944     0.212     0.301     0.351     0.387  ← weak
   N side      10,010     0.588     0.762     0.870     0.887  ← second-best after H

O              42,429     0.231     0.303     0.358     0.382
   O bb        29,944     0.253     0.338     0.395     0.422
   O side      12,485     0.241     0.373     0.543     0.566
```

### 18.4 What each calculator family contributes per atom type — narrative

**HYDROGEN (R² = 0.928, 17 ring-current kernels carry most of it)**

- Ring current (BS, HM, RingSusc): **dominant** mechanism. R² ≈ 0.78 from ring-current alone.
- EFG (Coulomb, ApbsField, MopacCoulomb, AIMNet2 EFG): strong **angular alignment** with target (cos→target ≈ 0.93, the highest of any group); EFG_aro single kernel reaches R² ≈ 0.71-0.76.
- Dispersion: meaningful (R² ≈ 0.39); highest cos→target ≈ 0.85, but spatially sparse (5 Å cutoff).
- Bond anisotropy (McConnell): **minor for H** because hydrogen has only one short bond; doesn't feel bond anisotropy directly.
- Para-quadrupole: small contribution.
- HBond: ~0 (mechanical mutants don't change backbone H-bonds).

**CARBON (R² = 0.562 pooled; CA = 0.627, CB = 0.647, C-side = 0.729, C=O = 0.463 — paramagnetic-dominated)**

- EFG (especially MOPAC-EFG): **dominant**, R² ≈ 0.43. The charge redistribution from ring removal changes the EFG angular pattern; MOPAC charges capture this, ff14SB doesn't. The MOPAC−ff14SB gap (+0.155) is the **charge-polarisation dimension** specifically for carbon.
- Ring current: **minor** (R² ≈ 0.16). Carbon shielding is paramagnetic-term-dominated; ring-current magnetic field affects diamagnetic which is secondary for C.
- Dispersion: enters forward selection at rank 4-5; short-range angular structure captures ring proximity for C.
- Bond anisotropy: McC_aromatic_total enters forward selection at rank 3 — the only bond kernel that matters (backbone unchanged in mechanical mutants).
- **C=O (carbonyl) drags carbon down** (0.463 vs 0.729 for sidechain C): sp² C with C=O π* excitation is paramagnetic-dominated in a way the geometric kernels capture poorly.
- **Sidechain C is strong**: closer to removed aromatic ring, more diverse angular geometry.

**NITROGEN (R² = 0.380 pooled; N bb = 0.387 (weak), N side = 0.887 (second-best after H))**

- Ring current: contributes (R² ≈ 0.33); collective from 17 kernels, no single one dominates. Isotropic effect on ¹⁵N is negligible (<0.6% per Han 2011) but T2 anisotropic effect is measurable.
- EFG: **weak for N** (R² ≈ 0.10). ff14SB and MOPAC give similar weak R². N lone pair + paramagnetic dominance mean the electrostatic perturbation has different effect than for H or C.
- π-Quadrupole: **most important for N** of any element. PQ_total enters forward selection at rank 2 (+0.043). r⁻⁵ quadrupole from ring's π-cloud provides angular info at distances where N atoms sit. PQ ⊥ BS (cos = 0.45, near-random in 5D) — genuinely independent dimension.
- Dispersion: normalisation reveals dispersion for N (4.5× increase after norm); enters 5 of last 10 forward-selection steps.
- MOPAC bond order: enters rank 4. Peptide-bond partial-double character (Wiberg 1.3-1.5) couples N electronically to the mutation site. Only element where bond-order weighting matters.
- **N side is second-best after H** because sidechain N (HID/HIE/ARG/LYS terminal) is in direct contact with the aromatic ring and ring identity matters strongly here (mutation type +Mut delta = +0.150 for N side).
- **N bb is hard** because backbone geometry is fixed in mechanical mutants — backbone N sees the perturbation only through long-range fields.

**OXYGEN (R² = 0.382 pooled; O bb = 0.422, O side = 0.566)**

- Ring current: contributes (R² ≈ 0.32); near-field RBF kernels (close to ring) help.
- EFG: moderate (R² ≈ 0.20).
- Dispersion: **largest contributor for O** in geo-only basis (R² ≈ 0.23, larger than ring current 0.20 or EFG 0.16 in geo-only — see master_chart). DispChi (dispersion scalar × RingSusc T2) accesses ring-proximity dimension where r⁻⁶ scalar is large.
- Mutation type sensitivity for O side: **+0.209**, the largest of any atom type. SER/THR/TYR sidechain O atoms are direct H-bond partners and ring-orientation reporters.

### 18.5 Per-element forward selection top 5 (which kernels enter the ridge first)

```
HYDROGEN (R² = 0.935, 4 of top 5 are ring-current)
1. MopacEFG_aro      +0.766    ← single kernel does most of the work
2. RbfBsNear_ring0   +0.097    ← ring-current near-field
3. RbfBsMid_ring0    +0.012
4. RbfBsMid_ring1    +0.010
5. BS_TYR            +0.006

CARBON (R² = 0.546)
1. MopacEFG_aro      +0.390    ← MOPAC-derived EFG dominates
2. EFG_aro           +0.037    ← ff14SB EFG adds angular detail
3. MC_aromatic_total +0.017    ← removed aromatic bond anisotropy
4. BS_HIE            +0.010    ← ring-type-specific
5. HM_TRP_benzene    +0.007

NITROGEN (R² = 0.434, multi-mechanism — no single dominant)
1. EFG_aro           +0.080    ← weak start; nothing as dominant as for H or C
2. PQ_total          +0.061    ← π-quadrupole — distinctive for N
3. BS_TRP_benzene    +0.026
4. MopacMC_total     +0.015    ← MOPAC bond-order weighting
5. RingSusc_total    +0.015

OXYGEN (R² = 0.387)
1. EFG_aro           +0.171    ← electrostatic dominant
2. RbfBsNear_ring0   +0.051    ← near-field ring current
3. MopacEFG_aro      +0.020
4. AngBsAxial_ring0  +0.017
5. HM_HIE            +0.009
```

### 18.6 Distance dependence of dimensionality (110-protein, re-run on 720 pending)

Three predictive dimensions at every distance band; R² varies:

```
ELEMENT   0-4 Å   4-8 Å   8-12 Å   12+ Å
H         0.919   0.959   0.969    0.908
C         0.501   0.556   0.951    0.853
N         0.661   0.224   0.906    0.818
O         NaN     0.195   0.745    0.748
```

Same 3 dimensions work everywhere for hydrogen. Carbon gets better at distance (near-field hard, far-field easy). Nitrogen near-field has 8 predictive dimensions — multi-mechanism complexity concentrates close to the ring.

### 18.7 Granularity gap — type level, not atom-name level

The stratification above is at AMBER-atom-type level:

```
Carbon:    CA  /  C=O  /  CB  /  C_side  (4 buckets; "C_side" = CG/CD/CE/CZ/CH all pooled)
Nitrogen:  N_bb  /  N_side                (2 buckets; "N_side" = ND/NE/NZ/NH all pooled)
Oxygen:    O_bb  /  O_side                (2 buckets; "O_side" = OG/OH/OD/OE all pooled)
Hydrogen:  one group
```

A finer per-AMBER-name stratification (e.g., separate buckets for PHE's CG / CD1 / CD2 / CE1 / CE2 / CZ; LYS's CG / CD / CE distinct; HIS's CD2 / CE1 / NE2 / ND1 distinct) **does not currently exist** in the analysis output. The existing script `learn/src/actual_physics/atom_type_calibration.py` does the type-bucketing shown above and stops there.

This is a known gap. The finer granularity would distinguish, for example:
- HIS NE2 vs ND1 (the two ring nitrogens — chemically inequivalent in HID vs HIE)
- TRP's CG / CD1 / CD2 / NE1 / CE2 / CE3 / CZ2 / CZ3 / CH2 (each at a distinct ring position)
- LYS CB / CG / CD / CE (gradient of distance from charged Nζ)
- ARG CB / CG / CD (gradient of distance from guanidinium)

Each of these would have qualitatively different relationships to ring-current, EFG, and bond-anisotropy contributions — pooling them under "C_side" or "N_side" mixes signals.

Per the user (2026-04-30): finer per-atom-name results would be ideal but don't currently exist. The bucketing in §18.3 is what the existing data supports.

### 18.8 Provenance for §18

All numbers in §18.2 - §18.6 are reproductions from the named `learn/` documents. Adversarial verification of these numbers means re-running the cited scripts (`learn/src/actual_physics/element_physics.py`, `learn/src/actual_physics/atom_type_calibration.py`) against the canonical 720-protein extraction — not in scope for the inventory-correctness vetting pass.

What IS in scope: the **per-calculator narrative** in §18.4 — does the description of how each calculator family contributes to each atom type accurately reflect what the data tables show? For example, "MOPAC bond order is the only element where bond-order weighting matters; this is for nitrogen" — is that a correct interpretation of the table?
