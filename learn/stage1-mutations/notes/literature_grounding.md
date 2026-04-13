# Literature Grounding for Stage 1 Claims

2026-04-12.  Assessment of what's established, what's extended, and
what's novel in the argument chain.

---

## Well-established (cite, don't derive)

### Different mechanisms produce different shielding contributions

Standard NMR theory.  Each mechanism has its own functional form:

- Ring current: n x B from wire loops or surface integral.
  Pople 1956 (concept), Johnson & Bovey 1958 (JB model),
  Haigh & Mallion 1979 (surface integral), Case 1995 (calibration).
- Electric field: sigma = A*E_z + gamma*V_ab.
  Buckingham 1960.
- Bond anisotropy: dipolar kernel from bond susceptibility.
  McConnell 1957.
- Different multipolar orders have different angular patterns.
  Stone 2013, Ch. 3 (T-tensor formalism).  Mathematical fact of
  the spherical harmonic expansion.

### Element-dependent mechanisms

Saito, Ando & Ramamoorthy 2010 (Prog NMR Spectrosc 57:181):
"for a 1H nucleus the relative chemical shift is more significantly
determined by contributions from sigma_d and sigma' as compared to
sigma_p.  For atoms with 2p electrons such as 13C, 15N and 17O...
the relative chemical shift is predominantly contributed by sigma_p."

This is the fundamental reason H responds to different physics than
C/N/O.  Diamagnetic term (sigma_d) is sensitive to through-space
fields (ring current, EFG).  Paramagnetic term (sigma_p) is sensitive
to local electronic structure (bond character, orbital symmetry).

Open access: PMC2905606.

### Tensor orientation varies site-specifically

Hall & Fushman 2006 (JACS 128:7855):
- 15N CSA in protein GB3 varies from -111 to -241 ppm
- Tensor orientation (beta_z angle) varies 7.5 to 27.6 degrees
- Phe52 and Trp43 are explicit outliers (largest |CSA|)
- They measured the variation but did NOT explain it mechanistically

Open access: PMC2519110.

Key quote for us: they observed aromatic-residue-dependent 15N CSA
outliers but could not attribute them to a specific mechanism.  Our
kernel decomposition provides the mechanistic explanation — at least
for the aromatic contribution.

### Ring current affects the full CSA tensor, not just isotropic shift

Boyd & Skrynnikov 2002 (JACS 124:1832):
Derived full shielding tensor from ring magnetization:
sigma_ab = -n_b * (H . n)_a.  CSA contributions up to 16.6 ppm
for backbone HN near aromatic rings.  This is the key paper
connecting ring currents to the TENSOR (not just the scalar shift).

We already have this PDF in references/.

### Ring current vs electric field: separately calculable, both matter

Sahakyan & Vendruscolo 2013 (JPC-B 117:1989):
Decomposed RNA base shifts into ring current and EFG contributions.
Found EFG dominates for non-hydrogen atoms.  But: this is for
ISOTROPIC shifts only.

Case 1995 (J Biomol NMR 6:341):
Calibrated ring current intensity factors against DFT.  Found
fitting ring current + electrostatic simultaneously gives the best
result.  Also isotropic only.

Xu & Case 2002 (Biopolymers 65:408):
DFT decomposition of 15N and 13C shifts into hydrogen bonding
(up to 8 ppm), backbone conformation (up to 13 ppm), and
neighborhood effects (up to 22 ppm).  Established the additive
model.  Isotropic only.

---

## Extended to T2 (cite the isotropic version + note the extension)

### Per-element mechanism dominance in the T2 channel

Literature: ring current dominates H isotropic shifts (Boyd 2002),
EFG dominates heavy atom isotropic shifts (Sahakyan 2013).

Our extension: this holds in the T2 channel for H (R²=0.909 from
ring current) and C (R²=0.431 from EFG), but breaks down for N
(multi-mechanism, no single family dominates) and O (dispersion
and EFG both contribute).

The N result connects to Hall & Fushman 2006: they observed 15N
CSA outliers near aromatics but couldn't explain them.  Our forward
selection shows N needs 5 kernel families (EFG, PQ, BS, MC,
RingSusc) each contributing 0.015-0.080 — consistent with the
paramagnetic-term dominance that makes 15N sensitive to multiple
perturbation channels (Saito 2010).

### Calibrated coefficients as physical constants

Case 1995 calibrated isotropic ring current intensity factors on
single-molecule DFT.  We calibrate the full T2 weight vector on
protein-embedded DFT across 110 proteins.  The extension is:
scalar → tensor, single molecule → 69K atoms, ring current only
→ all mechanism families.

### BS-HM agreement in the T2 channel

Case 1995 found no significant difference between JB and HM for
isotropic shifts.  We extend this to T2: magnitude ratio 1.009,
cosine similarity 0.9977-0.9999.

---

## Appears to be novel (needs careful framing)

### Tensor DIRECTION distinguishes physical mechanisms

The claim: |cos(BS, EFG)| = 0.684 across 58,193 atom-ring pairs.
The two kernel families produce T2 tensors pointing in different
directions at the same atom because they arise from different
multipolar orders.

Literature status: The MATHS is established (different multipole
orders → different angular dependence, Stone 2013).  The ISOTROPIC
decomposition exists (Sahakyan 2013).  Boyd 2002 derived the ring
current TENSOR.  But: nobody has systematically measured the angular
separation between kernel families in the T2 channel across a large
protein dataset.  This measurement is new.

How to frame: "The angular independence of different kernel families
in the T2 channel follows from the distinct multipolar symmetries
of the underlying interactions (Stone 2013).  We quantify this
independence empirically: mean |cos(BS, EFG)| = 0.684 across
58,193 atom-ring pairs..."

### Per-element T2 decomposition

Nobody has stratified the ANISOTROPIC shielding perturbation by
element and physical mechanism family.  The isotropic version exists
(Sahakyan 2013 for RNA, Case 1995 for protons).  The tensor version
at this scale is new.

### 3 effective T2 dimensions in the kernel space

PCA on geometric kernel T2 vectors showing 3 dominant blocks of 5
eigenvalues (one block per physical source geometry) is new.

### Ridge weight orthogonal to dominant variance PCs

The finding that prediction lives in the angular fine structure
(tail eigenvalues) rather than the dominant modes is a structural
finding about the kernel space geometry.  No precedent found in
NMR literature.  May have analogues in chemometrics / partial least
squares literature (PLS vs PCA) but we haven't looked.

### Dispersion as a significant T2 contributor for N and O

Dispersion is treated as minor in the isotropic shift literature.
Our T2 forward selection shows it enters at rank 2-5 for C, N, O.
This may be a genuine T2 finding — the r^-6 scalar is small but
the tensor (r^-8 anisotropic part, combined with ring susceptibility
via DispChi) captures short-range angular structure that isotropic
analysis misses.

---

## Papers we need but don't have

### Have as PDFs (in references/)

- Buckingham 1960 (A11) — buckingham-1960-chemical-shifts-polar-groups.pdf
- Case 1995 (C1) — case-1995-ring-current-calibration.pdf
- Boyd & Skrynnikov 2002 (A8) — calculations-of-the-contribution-of...pdf
- Sahakyan & Vendruscolo 2013 (C8) — sahakyan-vendruscolo-2013-...pdf

### Open access (fetch from PMC)

- Saito et al. 2010 — PMC2905606
- Facelli 2011 — PMC3058154
- Hall & Fushman 2006 — PMC2519110

### Need from library

- Pople 1956, J. Chem. Phys. 24, 1111.  DOI: 10.1063/1.1742701
- Johnson & Bovey 1958, J. Chem. Phys. 29, 1012.  DOI: 10.1063/1.1744645
- McConnell 1957, J. Chem. Phys. 27, 226.  DOI: 10.1063/1.1743676
- Haigh & Mallion 1979, Prog. NMR Spectrosc. 13, 303.  (no DOI found)
- Stone 2013, Theory of Intermolecular Forces, 2nd ed., OUP.
  ISBN: 978-0199672394.  BOOK — need Ch. 3 on T-tensors.
- Yao et al. 2010, JACS 132, 10866.  DOI: 10.1021/ja104711v
  (site-specific 15N CSA, ~166 ppm)
- Han et al. 2011, J. Biomol. NMR 50, 43.  DOI: 10.1007/s10858-011-9478-4
  (SHIFTX2, negligible ring current on 15N)
- Xu & Case 2002, Biopolymers 65, 408.  DOI: 10.1002/bip.10276
  (DFT decomposition of 15N/13C)
- Poon et al. 2004, JPC-B 108, 16577.  DOI: 10.1021/jp0370670
  (15N CSA and backbone angles)
- London 1937, Trans. Faraday Soc. 33, 8.  DOI: 10.1039/TF9373300008
  (dispersion force, if we feature dispersion prominently)
- Haeberlen 1976, High Resolution NMR in Solids, Academic Press.
  ISBN: 978-0120255610.  BOOK — for T0/T1/T2 convention declaration.

### Needed for the nonlinear signal finding

- Ramsey 1950, Phys. Rev. 78, 699.  DOI: 10.1103/PhysRev.78.699
  The paramagnetic shielding term: sigma_p involves 1/DeltaE energy
  denominators and orbital angular momentum matrix elements.  The
  product terms (geometry x geometry / DeltaE) produce nonlinear
  kernel interactions at paramagnetic-dominated atoms (N, O).
  This grounds the RF finding: nitrogen's +0.169 nonlinear signal
  comes from interaction effects in the paramagnetic term that
  linear ridge cannot capture.  Need at write time.

### Possible additions (found during search)

- Plasser 2021, Eur. J. Org. Chem. 2021.
  DOI: 10.1002/ejoc.202100352.
  "Visualisation of Chemical Shielding Tensors (VIST) to Elucidate
  Aromaticity."  Couldn't access but title suggests tensor direction
  used to interpret mechanism.  Worth checking.

- Tang & Case 2011, J. Biomol. NMR.  "Calculation of Chemical Shift
  Anisotropy in Proteins."  PMC3196061.  Hydrogen bonding dominates
  CSA variation, but no mechanism decomposition.  Secondary reference.
