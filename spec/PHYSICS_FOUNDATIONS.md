# Physics Foundations — NMR Shielding at the Protein Scale

**Status:** IN PROGRESS. Session 0 (landscape) started 2026-04-22. The
three drafting sessions (see `spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md`
section 13.11) have not begun. Do not treat the contents of this
document as drafted yet; most of it is placeholders or landscape
pointers.

**Purpose:** replace "the assistant's training plus a handful of
references" with a properly-referenced physics theory chapter for
every mechanism the extractor touches. Each section does four things:
origin and history, formal derivation with primary citations, our
classical approximation and where it breaks down, and time-resolved
vs snapshot distinction.

**Companion documents:**
- `GEOMETRIC_KERNEL_CATALOGUE.md` — what we compute. This document is
  the *why*.
- `MATHS_GOALS.md` — validation plan. This document justifies the
  plan.
- `CALCULATOR_PARAMETER_API.md` — 93 parameters with equations. This
  document puts those references in their theoretical context.
- `references/ANNOTATED_BIBLIOGRAPHY.md` — paper-level index. Every
  citation below has a PDF in `references/` and a bibliography entry.

---

## 0 — Session 0 landscape (in progress, 2026-04-22 pm)

Method per user research philosophy: "ask the biggest questions you
are addressing and see what people have to say, even if it doesn't
have immediate application to what you are doing. That way we have
checked out the haystack for needles. Then we can come back and build
an entire little mini-castle of needles and pins, and know we are
not leaving out things that people will see and note."

Queries run and their shape below. Session 0 goal: sweep broadly, no
NMR qualifier, surface literature the focused searches would have
missed. Output is this landscape list plus papers fetched to
`references/` and `ANNOTATED_BIBLIOGRAPHY.md` updates.

### 0.1 Query: diamagnetic susceptibility of proteins — biophysics

Surfaced a coherent literature area that NMR shift prediction people
largely don't cite. The MRI / magnetic-biophysics community has been
computing bulk susceptibility anisotropy of proteins for decades.

Key finds:
- **Worcester, D.L. (1978)** "Structural origins of diamagnetic
  anisotropy in proteins." *Proc. Natl. Acad. Sci. USA* 75,
  5475--5477. DOI: 10.1073/pnas.75.11.5475 (PMID 281695;
  PMC392987). Foundational structural work linking peptide-bond
  orientation and aromatic-ring per-atom anisotropy to bulk protein
  anisotropy — especially the α-helix axial-peptide alignment. Same
  physics as our McConnell peptide CO / PeptideCN bond kernel and
  our ring-current kernels, summed to the bulk scale. **Correction
  2026-04-22:** earlier draft of this section mis-attributed this
  paper to Pauling and conflated it with Pauling's separate 1979
  piece — see next bullet.
- **Pauling, L. (1979)** "Diamagnetic anisotropy of the peptide
  group." *Proc. Natl. Acad. Sci. USA* 76, 2293--2294. DOI:
  10.1073/pnas.76.5.2293 (PMC383585). Simple theory of planar
  peptide-group magnetic anisotropy, giving the classical molar
  value Δχ_peptide = −5.36 × 10⁻⁶ cm·g⁻¹·sec⁻¹ emu. Direct
  quantitative target: our calibrated McConnell peptide-CO /
  PeptideCN intensities should map to this per-bond anisotropy
  under an appropriate unit conversion.
- **Babaei, M., Jones, I.C., Dayal, K. & Mauter, M.S. (2017)**
  "Computing the Diamagnetic Susceptibility and Diamagnetic
  Anisotropy of Membrane Proteins from Structural Subunits." *J.
  Chem. Theory Comput.* 13(6), 2945--2953. DOI:
  10.1021/acs.jctc.6b01251. Modern computational realisation:
  hierarchical subunit summation of per-amino-acid χ tensors into
  full-protein bulk χ. Directly parallel to our per-atom kernels
  but summed for MRI contrast prediction. **Authorship corrected
  2026-04-22** — earlier drafts here and in conversation carried a
  stale "Palmer et al." attribution; the correct first author is
  Babaei. (Arthur Palmer, the NMR-relaxation theorist cited in
  bib M2 / M3, is a different researcher.)
- **Liu (2010) and follow-up work** — Susceptibility Tensor Imaging
  (STI) in MRI. "Magnetic Susceptibility Anisotropy of Human Brain
  in vivo and its Molecular Underpinnings" (PMC3254777). Establishes
  that macroscopic brain MRI susceptibility anisotropy originates
  from oriented molecular-scale magnetic anisotropy of peptide bonds
  and lipid bilayers. Parallel to our argument at molecular scale.
- **Bras et al. (2014)** "The Diamagnetic Susceptibility of the
  Tubulin Dimer" — protein-scale susceptibility calculation from
  structural subunits.

**Implication for our work:** the bulk-scale biophysics literature
uses the same per-atom / per-bond magnetic-anisotropy physics we do.
Our WT-ALA delta isolates the aromatic contribution at the atomic
scale; the MRI community sees the same physics manifested
macroscopically. This is territory a reviewer in the broader field
would expect us to know exists.

### 0.2 Query: classical electromagnetic calculation of magnetic properties at molecular scale

The textbook classical EM doesn't directly apply — at molecular scale
quantum mechanical effects dominate. But the "hybrid classical kernel
+ quantum-calibrated parameters" approach we use is a recognised
strategy.

Key finds:
- **Kahn et al. (2022)** "Molecular Magnetizabilities Computed Via
  Finite Fields: Assessing Alternatives to MP2 and Revisiting Magnetic
  Exaltations in Aromatic and Antiaromatic Species"
  (PMC8903098). Modern benchmark of methods for computing molecular
  magnetizabilities; discusses magnetic exaltation as an aromaticity
  criterion directly tied to ring currents; finds κ-OOMP2 is the best
  MP2-based method.
- **Helgaker, Jaszuński & Ruud (1999)** "Ab Initio Methods for the
  Calculation of NMR Shielding and Indirect Spin–Spin Coupling
  Constants." *Chem. Rev.* 99, 293. **Canonical review of magnetic
  property calculations.** Needs to be fetched and read in Session 1.
- **Helgaker, Coriani et al. (2012)** "Recent Advances in Wave
  Function-Based Methods of Molecular-Property Calculations."
  *Chem. Rev.* 112, 543. Extended review covering current state of
  the field.
- **Gershoni-Poranne & Stanger (2015)** "Magnetic criteria of
  aromaticity." *Chem. Soc. Rev.* 44, 6597. Critical review of
  magnetizability exaltation and related aromaticity measures.
- **London (1937)** "Théorie quantique des courants interatomiques
  dans les combinaisons aromatiques." *J. Phys. Radium* 8, 397 — the
  origin paper for GIAO (London orbitals) and for classical ring
  current theory. Already in our bibliography target; this surfaced
  its dual role in modern quantum magnetic property calculations.
- **Flygare & Benson (1971)** "The molecular Zeeman effect and
  magnetic susceptibilities of small organic molecules." *Mol.
  Phys.* 20, 225. Historical magnetizability calculations.

**Implication:** the magnetizability literature is a parallel
tradition to chemical-shift prediction, using the same underlying
physics (second-order energy response to magnetic field), but with
different foci (aromaticity criteria, method benchmarking). Relevant
to part 2 (ring currents) and part 5 (higher multipoles) of the
physics foundations document.

### 0.3 Query: time-correlation functions in molecular dynamics simulation

Surfaced the stat-mech foundation underneath Lipari-Szabo. The MD
community's canonical framework for extracting observable dynamics
from trajectories is broader than the NMR-relaxation corner of it.

Key finds:
- **Berne & Harp (1970)** "Time-Correlation Functions, Memory
  Functions, and Molecular Dynamics." *Phys. Rev. A* 2, 975. Canonical
  paper; underlies everything in the field.
- **Kubo (1957) and follow-up** — linear response theory and the
  fluctuation-dissipation theorem connect autocorrelation functions
  to macroscopic response coefficients. This is the stat-mech bridge
  between our per-frame signals and their macroscopic meaning.
- **Zwanzig memory function formalism** — generalises correlation
  functions beyond exponential decay. Relevant when proteins show
  non-exponential relaxation (they do).
- **Ditler & Luber (2022)** "Vibrational spectroscopy by means of
  first-principles molecular dynamics simulations." *WIREs Comput.
  Mol. Sci.* 12, e1605. Directly parallel to our approach: compute an
  observable per frame (dipole / IR spectrum) and autocorrelate to
  get a physical quantity. Recent review worth fetching.
- **GROMACS correlation function documentation** — practical
  implementation, error bars, cutoffs.

**Implication:** our proposed bond-vector autocorrelations and
kernel-time-series treatment aren't an NMR-specific invention;
they're a standard MD-community technique whose theoretical
underpinning is stat-mech linear-response / fluctuation-dissipation.
The physics-foundations chapter on time-resolved dynamics (part 8)
should ground in Berne-Harp / Zwanzig before jumping to Lipari-Szabo.

### 0.4 Query: tensor formalism in geometric deep learning — "tensor-as-array" vs "tensor-as-physicist"

User-raised 2026-04-22 pm: we use "tensor" ambiguously in our own
work — sometimes a Cartesian Mat3, sometimes a SphericalTensor with
irreducible T0+T1+T2, sometimes an e3nn-style indexed multi-array.
The equivariant-ML community has formalised this distinction and we
should cite their framework rather than carrying it informally.

Key finds:
- **Geiger, Smidt et al. (2022)** "e3nn: Euclidean Neural Networks."
  arXiv:2207.09453. We already use this library in
  `learn/c_equivariant/`; the paper itself is the canonical
  specification of the tensor formalism: irreducible representations
  of O(3), `TensorProduct` via Clebsch-Gordan coefficients,
  equivariance by construction.
- **Kondor, R. (2025)** "The principles behind equivariant neural
  networks for physics and chemistry." *PNAS* 122(41), e2415656122.
  DOI: 10.1073/pnas.2415656122. Open access (PMC12541325). The
  pedagogical review of why irreducible representations matter when
  tensor values have physical meaning. **Makes the tensor-as-array
  vs tensor-as-physicist distinction explicit** and resolves it via
  irreps. Canonical review. (Attribution corrected 2026-04-22 pm
  from "Smidt et al." during community 6 depth pass — the paper is
  single-author Kondor.)
- **Gerken et al. (2023)** "Geometric deep learning and equivariant
  neural networks." *Artif. Intell. Rev.* 56, 14605. Comprehensive
  survey.
- **Batzner et al. (2022)** "E(3)-equivariant graph neural networks
  for data-efficient and accurate interatomic potentials" (NequIP).
  *Nat. Commun.* 13, 2453. Applied directly to molecular property
  prediction.
- **Thomas et al. (2018)** "Tensor Field Networks." arXiv:1802.08219.
  The tensor-field-network paper that precedes e3nn.
- **Cormorant, SEGNN, EGNN** — the family of equivariant
  architectures; all build on the same irreducible-representation
  machinery.

**Implication:** when we write PHYSICS_FOUNDATIONS part 1 (shielding
fundamentals) and the T2/T1/T0 decomposition discussion, the formal
framework to cite is the equivariant-ML irreducible-representation
machinery, not just the physics textbook tensor decomposition.
Field manifests (rollup section 13.8) can carry explicit irrep
metadata (`0e`, `1o`, `2e` in e3nn notation) so downstream ML
operates on typed physical tensors not arbitrary arrays. Probably
merits a dedicated appendix in `PHYSICS_FOUNDATIONS.md` on "tensor
formalism and its two meanings in our work."

### 0.5 Query: equivariant tensor time-series learning — anisotropic data in training over time

User-raised 2026-04-22 pm: this is in-progress research, not a
settled field; explicitly flagged as "hard stuff and likely in
progress." Surfaced enough to cite and to frame our Stage 3 GNN /
time-series work against.

Key finds:
- **arXiv:2406.01552 (2024)** "Learning equivariant tensor
  functions with applications to sparse vector recovery."
  Universally expressive equivariant architectures exploiting
  diagonal actions of orthogonal, Lorentz, and symplectic groups;
  applications include time series.
- **arXiv:2108.09541** "Rotation equivariant operators for machine
  learning on scalar and vector fields." Tensor-field convolutions,
  linear equivariant operators.
- **Nat. Commun. (2026)** "Accurate prediction of tensorial spectra
  using equivariant graph neural network"
  (DOI 10.1038/s41467-026-69159-9). The *tensorial* part matters —
  predicting tensor-valued observables from an equivariant GNN.
  Directly applicable to Stage 3 if the learning system goes
  GNN-for-shielding-tensor.
- **Springer AIR survey (2024)** "Rotation invariance and
  equivariance in 3D deep learning: a survey." Generic survey of
  the field.
- **Chen-Cai-OSU/awesome-equivariant-network** GitHub collection.
  Curated paper list; entry point for the field.
- **Anisotropic-sequential-tensors approach** (noted in one of the
  search results): encodes isotropic sequential scalar components
  and anisotropic sequential tensor components into spherical
  tensor representations for symmetry-aware prediction under
  crystalline symmetry constraints. Nearest conceptual neighbour to
  our use case — trajectories of anisotropic shielding tensors per
  atom.

**Implication:** Stage 3 GNN work (or equivalent time-series
learner) can look to this literature for architectures that handle
tensor-valued sequences. None of it is yet settled, so we bring our
own use-case constraints and cite the closest existing work as
direction, not authority. This is the literature that feeds the
13.11 Session 3 drafting work for parts 8 (time-resolved) and 10
(connections to code).

### 0.6 Query: geometric kernels in robotics and protein structure — SE(3), Lie groups

User intuition 2026-04-22 pm: "I bet robots have geometric kernels.
I bet other people have sought geometric models for proteins which
extend beyond my graduate structural biology course. Those
abstractions might be gold to us (or not but we should check)."
Confirmed gold.

Key finds:
- **EquiCPI (2025)** "SE(3)-Equivariant Geometric Deep Learning for
  Structure-Aware Prediction of Compound-Protein Interactions."
  arXiv:2504.04654.
- **Nat. Commun. (2025)** "SE(3)-equivariant ternary complex
  prediction towards target protein degradation" (DOI
  10.1038/s41467-025-61272-5). Equivariant prediction of molecular
  complexes.
- **VN-EGNN** "E(3)- and SE(3)-Equivariant Graph Neural Networks
  with Virtual Nodes Enhance Protein Binding Site Identification"
  (OpenReview). Relevant because *virtual nodes* in the equivariant
  framework are how you place kernel-probe points at arbitrary
  spatial positions rather than only at atom centres — directly
  applicable if we want to evaluate our classical kernels at
  non-atom positions (volumetric grids, ring centres, bond
  midpoints).
- **Satorras et al. (2022)** EGNN "E(n) Equivariant Graph Neural
  Networks." Foundational simpler equivariant GNN.
- **PLOS Comp. Biol. (2023)** "E(3) equivariant graph neural
  networks for robust and accurate protein-protein interaction site
  prediction." Binding-site class of problems.
- **Chirikjian** "Stochastic Models, Information Theory, and Lie
  Groups" — robotics-origin reference for SE(3) kernel mathematics
  (Gaussian processes on Lie groups, invariant measures, exponential
  map).
- **Helgason** "Groups and Geometric Analysis" — foundational
  harmonic analysis on Lie groups.

**Implication:** there is a rich, mature vocabulary —
SE(3) / tangent spaces / Lie algebra so(3) / exponential map /
adjoint representation / invariant measures / Gaussian processes on
Lie groups — that protein computational work and robotics both use.
Our classical kernels are implicit functions on SE(3)-configuration
space; parameterising them explicitly on SE(3) (rather than as
ad-hoc functions of distance + angle) generalises cleanly and has
been actively developed in both robotic manipulation and drug
discovery communities. Stage 3 GNN work should borrow this formalism
explicitly rather than re-deriving it in ad-hoc ring-current form.

### 0.7 Queries still to run (pending)

**Physics haystacks:**
- Paramagnetism in biomolecules — metal centres, radicals,
  paramagnetic relaxation enhancement and pseudocontact shifts, even
  though our calibration set is diamagnetic.
- GIAO and alternatives / gauge-origin problem in magnetic property
  calculation — CSGT, IGLO, LORG, LRESC, recent comparisons. Where
  r2SCAN sits in the method landscape.
- MD-derived molecular property calculation — general literature on
  compute-on-frames-average across chemistry. Error bars,
  ergodicity checks, frame-selection criteria.
- Semiempirical methods for magnetic properties — what PM7 / MOZYME
  are known to capture vs miss beyond charges.
- Ensemble averaging of computed molecular observables — stat-mech
  view across chemistry.
- Chemical shielding anisotropy as a field — angular (T2) components
  in their own right, beyond isotropic prediction.
- Classical predictors of magnetic shielding as a research area —
  historical and current, any field framing.

**Methods / formalism haystacks (user-raised 2026-04-22 pm, second
batch, three covered above as 0.4-0.6):**
- Gaussian processes on Lie groups / SE(3) kernels — robotics
  manipulation literature to read explicitly.
- Diffusion models for protein structure generation (AlphaFold 3,
  Chroma, FrameDiff, DiffDock) — not directly our domain but the
  SE(3)-equivariant diffusion formalism has cross-applicability.
- Neural operator methods for PDE solution on molecular domains —
  physics-informed neural operators, DeepONet, Fourier Neural
  Operator. If we ever wanted to predict full field maps rather
  than per-atom values.

### 0.8 Shape of the haystack — corrected after bibliography cross-check

**Correction, 2026-04-22 pm during pre-clear cleanup.** During the
first and second passes of Session 0 I made claims about "thin
paper coverage" in our bibliography for several communities without
cross-referencing against `references/ANNOTATED_BIBLIOGRAPHY.md`.
That bibliography is substantial — 879 lines, 14 core references
across 12 sections (A through L) — and already strongly covers
several of the communities I surfaced. Corrected analysis below.

**Physics-side research communities:**

1. **NMR chemical shift prediction** — Case, Mulder, Sitkoff,
   Oldfield, de Dios, Hansen, Fushman. **Bibliography coverage
   strong** (sections A, B, C, J).
2. **Magnetic biophysics / MRI susceptibility imaging** — Pauling
   lineage, Liu / STI, Palmer et al. Different framing, same
   per-atom physics summed to the macroscopic scale.
   **Bibliography coverage genuinely thin.** Surfaced by query 0.1.
   **Needs new bibliography section** (call it section M:
   Magnetic Biophysics / MRI Susceptibility).
3. **Quantum-chemistry magnetic-property-calculation** — Helgaker,
   Ruud, Coriani, Gauss, Schleyer, Gershoni-Poranne. GIAO/CSGT/IGLO
   methods, magnetizabilities, NICS, aromaticity via magnetic
   response. **Bibliography coverage thin** — section B covers DFT
   shielding calculations but not the broader magnetic-property
   literature (magnetizabilities, NICS, aromaticity-via-magnetic-
   response). Surfaced by query 0.2. **Needs bibliography extension
   in section B or new section.**
4. **Stat-mech / MD time-correlation-function theory** — Berne,
   Zwanzig, Kubo, Palmer (Arthur — NMR relaxation theorist),
   Ditler-Luber. **Bibliography coverage thin at the stat-mech
   foundations level** — section H covers MD-for-NMR but not the
   broader TCF / linear-response theory (Berne-Harp 1970, Kubo 1957,
   Zwanzig memory-function formalism). Surfaced by query 0.3. **Needs
   bibliography extension.**

**Methods / formalism-side research communities:**

5. **Equivariant ML / irreducible-representation tensor formalism**
   — Geiger, Smidt, Batzner, Thomas. **Bibliography coverage
   strong** — section E has [CORE] E1 Thomas 2018 Tensor Field
   Networks, E2 Geiger-Smidt 2022 e3nn, E4 Batzner 2022 NequIP, plus
   other entries. Query 0.4 genuinely added: Kondor 2025 PNAS
   (canonical pedagogical review, now E15 in the bib; originally
   recorded in 0.4 as "Smidt et al." — attribution corrected
   2026-04-22 during community 6 depth pass), Gerken et al. 2023
   survey (now E16). The broader framework was already well
   covered; my "thin" claim in 0.4 was wrong.
6. **Equivariant tensor time-series learning** — **Bibliography
   coverage partial** — section J has [CORE] J2 Venetos 2023
   "equivariant GNN for full shift tensors" which is the closest
   prior art to our Stage 3 tensor-prediction problem. Query 0.5
   genuinely added: Nat. Commun. 2026 tensorial spectra,
   arXiv 2406.01552. Field is young so the truly new papers are
   in-progress; the conceptual framework is already bibliographed.
   My "thin because field is young" claim holds only for the most
   recent work.
7. **SE(3) / Lie-group geometric deep learning and robotics
   crossover** — **Bibliography coverage thin on the robotics
   origin and recent protein-structure applications.** Section F
   covers geometric descriptors (SOAP etc.) but not the SE(3)
   Lie-group formalism explicitly; section G (GNN message passing)
   may touch it but the cross-community robotics side (Chirikjian,
   Helgason) is not represented. Recent protein-SE(3)-equivariant
   papers (EquiCPI 2025, SE(3) ternary complex 2025, VN-EGNN) are
   too recent to be there. **Needs bibliography extension** —
   probably in section F or a new "SE(3) and geometric protein
   learning" section.

**Corrected summary:** of seven research communities identified,
**three are genuinely under-covered in the bibliography and need
new sections or substantial additions** (community 2 magnetic
biophysics, community 3 QM magnetic-property-calculation at the
magnetizabilities / aromaticity level, community 7 SE(3) robotics /
recent protein-SE(3)-equivariant). **Two are partially
under-covered** (community 4 stat-mech TCF foundations beyond
MD-for-NMR, community 6 recent equivariant tensor time-series
papers). **Two are well covered in the bibliography already**
(community 1 NMR shift prediction, community 5 equivariant ML core).

**Implication for Session 0 completion:** the bibliography addenda
are smaller than the first-pass claim suggested. The three
drafting sessions (13.11) should cite liberally against existing
bibliography entries and add the genuinely new papers as focused
extensions to appropriate sections, not as a wholesale rebuild.

### 0.9 Integration — synergies, validation benches, and caveats

Session 0 turned into three interleaved literature threads, walked
against concrete `TrajectoryResult` and kernel-catalogue contents.
Captured here as the handoff map for Sessions 1–3 drafting. Added
2026-04-22 pm.

**Thread 1 — Limits of reductive MD + what only exists in motion.**
Two questions, answered separately.

*(A) Reductive-MD convergence.* The 50 ns window adopted by the user
is at the practical optimum for fast methyl motions per current-
force-field benchmarks (Hoffmann 2018 bib H13); well short of the μs
regime that resolves slow ring flips and chemical exchange. Effective
sample size diagnostics (Lyman-Zuckerman 2007 H8; Zhang-Bhatt-
Zuckerman 2010 H9; Romo-Grossfield 2011 block covariance H11) are
runnable against the existing 10-protein analysis H5s.

*(B) Observables visible only in motion.* Cross-correlated relaxation
(Tugarinov-Kay 2003 M10; Ferrage et al. 2008 M11), RDCs in aligned
media (Tolman et al. 1995 M13; Lakomek et al. 2008 SCRM M15), measured
conformational entropy (Kasinath-Wand 2013 entropy meter M7; Sharp-
O'Brien-Kasinath-Wand 2015 backbone ΔS M9), dynamic allostery via
directional covariance (Amadei-Linssen-Berendsen 1993 M19; Cooper-
Dryden 1984 M23), memory-kernel non-exponentiality (Kou-Xie 2004 M21,
Zwanzig 2001 M22).

**Thread 2 — CSA as a field.** The σ-tensor principal components our
kernels predict per atom have their own experimental and ab-initio
tradition measuring exactly those values.

*Ab initio:* Havlin-Oldfield 1997 / 2001 (B12 / B13) established that
Cα CSA principal elements are individually sensitive to (φ, ψ, χ₁),
not just isotropic shift.

*Measured:* Wylie et al. 2011 PNAS (L6) used per-residue measured CSA
tensors as ultra-high-resolution structural restraints — the same
thesis claim our T2 work makes, tested by the solid-state NMR
community. Yao-Grishaev-Cornilescu-Bax 2010 (L7) reports α-helix ¹⁵N
CSA ≈ −173 ± 7 ppm vs β-sheet ≈ −162 ± 6 ppm — quantitative per-
stratum benchmark. Wylie-Franks-Rienstra 2006 (L9) finds δ₂₂ of ¹³C'
dominated by CO···HN H-bond length — direct test for our HBond kernel
contribution to carbonyl σ. Jordan-Rule-Tjandra 2007 (H14) is the
smaller-scale, lower-body-order precedent for our calibrated-kernel-
over-MD CSA pipeline.

**Thread 3 — Magnetic biophysics bridge.** Our per-atom magnetic-
anisotropy kernels have a macroscopic sibling: bulk protein χ
anisotropy computed in the MRI / biophysics tradition.

*Historical anchors:* Worcester 1978 (N1) traces α-helix anisotropy
to axial peptide-bond alignment. Pauling 1979 (N2) gives the classical
molar value Δχ_peptide = −5.36 × 10⁻⁶ cm·g⁻¹·sec⁻¹ emu.

*Modern realisation:* Babaei et al. 2017 (N3) sums per-subunit χ
hierarchically to bulk protein χ; aromatic content + β-barrel topology
dominate volumetric anisotropy. Li et al. 2017 (N4) STI review shows
the same per-bond-χ summation at tissue scale (peptide-bond-dominated
in myofibers).

### External validation benches surfaced

| Bench | Requires | Primary reference |
|---|---|---|
| Bulk protein χ anisotropy from per-atom kernel sum | Orientation-aware summation | Babaei et al. 2017 (N3) |
| Per-peptide Δχ vs calibrated McConnell CO/CN | McConnell intensity × unit conversion | Pauling 1979 (N2) |
| Per-residue CSA tensor principal components | σ prediction per backbone atom | Wylie 2011 (L6); Yao-Bax 2010 (L7) |
| α-helix vs β-sheet ¹⁵N CSA split (~11 ppm) | Per-stratum σ prediction | Yao-Bax 2010 (L7) |
| δ₂₂(¹³C') vs CO···HN distance | HBond kernel contribution to carbonyl σ | Wylie 2006 (L9) |
| NH DD/CSA CCR rate (rollup W3) | σ tensor + NH vector per frame | Tugarinov-Kay 2003 (M10) |
| RDC residuals after S² | `bond_orientation_tensor` (rollup W1) | Lakomek 2008 (M15) |
| Methyl entropy meter | `order_parameter_S2` per methyl | Kasinath-Wand 2013 (M7) |
| Backbone ΔS from amide S² | `order_parameter_S2` per NH | Sharp et al. 2015 (M9) |
| PCS from external Ln-tag | Angular kernel forward-prediction | Paramagpy (M26) + ubiquitin/GB1 datasets (N5–N7) |

### Honesty caveats — what a 50 ns window cannot observe

| Observable | Required window | Honest scoping |
|---|---|---|
| Slow aromatic ring flips (k ≈ 10¹–10² s⁻¹) | ms+ | `ring_normal_autocorrelation` captures sub-ns breathing only; flip kinetics absent (Akke-Weininger 2023 M17). |
| Fast aromatic ring flips (k ≫ 10³ s⁻¹) | sub-ms | Rate not directly measurable in 50 ns even when fast. |
| Chemical exchange (CPMG / CEST / R1ρ) | μs–ms | Rates inaccessible; slow-interconversion populations partly accessible (Palmer 2014 M3). |
| μs-timescale loop / domain motions | μs+ | Backbone S² for those sites not converged at 50 ns (H8, H13). |
| Overall tumbling τ_c | τ_c itself often > 50 ns | Underdetermined from trajectory; literature τ_c per protein used for W3 derivations. |

### Where each thread stitches into sections 1–10

- **Section 1 (Shielding fundamentals):** Havlin-Oldfield 1997 / 2001
  tensor-principal-components framing (B12, B13); Wylie 2011 PNAS
  demonstration that T2 carries structural information (L6).
- **Section 2 (Ring current physics):** Pauling 1979 per-peptide Δχ
  anchor (N2); Babaei et al. 2017 hierarchical-summation bridge from our
  per-atom kernels to macroscopic χ (N3).
- **Section 3 (Bond anisotropy / McConnell):** calibrated per-bond Δχ
  vs Pauling 1979 −5.36 × 10⁻⁶ emu target (N2).
- **Section 6 (H-bond):** Wylie 2006 δ₂₂(¹³C') vs CO···HN length
  (L9) as the direct HBond-kernel diagnostic.
- **Section 8 (Time-resolved NMR / MD payload):** thread 1's core
  territory. Berne-Harp / Kubo / Zwanzig (M22) stat-mech foundations;
  Lipari-Szabo 1982 (M1); Palmer reviews (M2 / M3); Henzler-Wildman
  & Kern 2007 (M4); Frauenfelder 1991 (M5); ring flips (M17); dynamic
  allostery (M23, M19). The `TrajectoryResult` fields (rollup section
  5, amended 2026-04-22) are the computational realisation.
- **Section 9 (Proxies vs first principles):** our McConnell /
  RingSusceptibility kernels as per-atom classical proxies for what
  Helgaker-Jaszuński-Ruud 1999 (B10) magnetic-property-calculation
  methods would yield; Babaei et al. 2017 (N3) is the intermediate scale.
- **Section 10 (Connections to code and H5):** each `TrajectoryResult`
  field maps to its thread-1 literature anchor; each kernel maps to
  its thread-2/3 validation bench.

### Session 0 deliverables as of 2026-04-22 pm

- **Bibliography additions:** B7–B13 (GIAO alternatives + r²SCAN +
  Oldfield CSA surveys), H7–H14 (MD convergence, sampling quality,
  CSA-on-MD parameterisation), L6–L9 (measured CSA tensors), M1–M26
  (new section: Protein Dynamics & NMR Relaxation), N1–N7 (new
  section: Magnetic Biophysics). Every added entry carries a DOI
  except Kutzelnigg 1980 (B8, old *Isr. J. Chem.* — library-only)
  and Zwanzig 2001 (M22, book with ISBN).
- **Rollup amendments** (`spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md`):
  `residue_residue_covariance` retyped scalar → (R, R, 3, 3) primary
  with Frobenius-norm scalar as derived companion;
  `bond_orientation_tensor` typing confirmed against Lakomek 2008
  SCRM; new per-amide cross-correlated-relaxation block
  (`nh_dipolar_csa_angle_autocorr`, `nh_dipolar_csa_ccr_rate`,
  `nh_dipolar_csa_principal_axis`) with three named caveats on
  τ_c convergence, T2-magnitude threshold, CSA-magnitude dependence.
- **Error correction:** Session 0.1 attribution mistake fixed —
  "Structural origins of diamagnetic anisotropy in proteins" is
  Worcester 1978 (PMC392987, PMID 281695), not Pauling; Pauling 1979
  is "Diamagnetic anisotropy of the peptide group" (PMC383585).
  Memory note updated in tandem.

### Three thesis-grade synergies

1. **NH DD/CSA cross-correlated relaxation back-calculation (rollup
   W3).** Per-frame σ(T0 + T1 + T2) + per-frame bond vectors +
   trajectory-length autocorrelation, at ~4000-atom scale — nothing
   else in the literature produces this from a calibrated classical
   model. External bench: published Γ_DD,CSA for GB3 / ubiquitin.
2. **Per-residue CSA tensor validation against measured components.**
   Predict σ₁₁ / σ₂₂ / σ₃₃ for backbone C' and amide N in ubiquitin
   / GB3 / GB1; compare to Wylie 2011 (L6), Yao-Bax 2010 (L7),
   Wylie 2006 (L9). The exact thesis claim that T2 carries structural
   information, externally validated by the solid-state NMR community.
3. **Per-atom-to-bulk χ anisotropy reconstruction.** Orient-aware sum
   of our McConnell + RingSusceptibility + HBond kernels over a
   protein, compared to Babaei et al. 2017 hierarchical-subunit bulk χ. A
   first-of-its-kind bridge claim: our calibrated per-atom kernels
   recover the macroscopic diamagnetic anisotropy independently
   computed by the MRI community.

Handoff: Sessions 1–3 drafting can begin against this scaffold. The
bibliography addenda are smaller than the first-pass claim suggested,
but three strong thesis-grade synergies give the thesis novelty a
clear shape. The literature pass does not replace the drafting — it
prepares it.

---

## 1 — Shielding fundamentals

*(Session 1 target. Placeholder.)*

σ tensor definition and irreducible decomposition (T0+T1+T2), sign
conventions, chemical shift vs shielding, time-averaging and
motional narrowing, DFT shielding via GIAO, r2SCAN functional, the
WT-ALA delta as aromatic isolator. Key references to write against:
Ramsey, Facelli, Ditchfield, Furness et al. 2020, plus Helgaker-
Jaszuński-Ruud 1999 surfaced in Session 0.

---

## 2 — Ring current physics

*(Session 1 target. Placeholder.)*

Aromatic π-electron origin (Pauling, London), quantum treatment
(Pople 1956), Johnson-Bovey classical loop, Haigh-Mallion surface
integral, Case calibration. Modern reviews (Mulder, Lampert). Plus
the aromaticity-via-magnetizability-exaltation lineage surfaced in
Session 0 (Gershoni-Poranne & Stanger 2015, Kahn et al. 2022).

Connection to MRI susceptibility biophysics (Session 0 finding):
what we compute per-atom is what the bulk-susceptibility-anisotropy
community sums up.

---

## 3 — Bond anisotropy (McConnell)

*(Session 2 target. Placeholder.)*

McConnell 1957, ApSimon-Craig 1967, Bothner-By. Full asymmetric
tensor derivation. Peptide-bond diamagnetic anisotropy links to
Pauling 1979 and to the bulk-protein-susceptibility literature.

---

## 4 — Electric field effects (Buckingham)

*(Session 2 target. Placeholder.)*

Buckingham 1960, Raynes 1962, Baker 2001 (APBS), Stewart 2013 (PM7).

---

## 5 — Aromatic quadrupole

*(Session 2 target. Placeholder.)*

Stone T-tensor, Pyykkö for quadrupolar-nucleus interpretation.

---

## 6 — H-bond shielding

*(Session 2 target. Placeholder.)*

Cornilescu-Bax 1999.

---

## 7 — Dispersion

*(Session 2 target. Placeholder.)*

London 1937 (re-used here), Brooks 1983.

---

## 8 — Time-resolved NMR and the MD-derived payload

*(Session 3 target. Heaviest fresh literature work. Placeholder.)*

Start from stat-mech foundations (Berne-Harp 1970, Kubo linear
response, Zwanzig memory functions — all Session-0 surfaces) before
descending to NMR-specific Lipari-Szabo 1982 and Palmer relaxation
reviews. Then RDCs (Tjandra-Bax), and the MD-community general
approach to ensemble-averaged observables (Ditler-Luber 2022 analog).

This chapter is the theoretical backing for the TrajectoryResult
fields in the rollup spec section 5.

---

## 9 — Proxies vs first principles

*(Session 3 target. Placeholder.)*

Where we use classical proxies (Coulomb as E-field proxy, point-
dipole approximations in ring susceptibility, fixed ff14SB/CHARMM
charges as electronic-structure proxy), what rigorous first-
principles would require, and the explicit scope decision: our
classical contributors must be understood rigorously and honestly;
first-principles QM beyond the classical kernels is outside scope
(dft-ex1's territory — read this project before writing the section).

---

## 10 — Connections to code and H5

*(Session 3 target. Placeholder.)*

Each calculator → its physics section. Each TrajectoryResult field →
its time-resolved-physics section. Each signal → physical
interpretation anchored to a reference.

---

## Session 0 progress log

- **2026-04-22 pm — first pass.** Three broad physics queries run
  (0.1 diamagnetic-protein biophysics, 0.2 classical EM molecular
  magnetic properties, 0.3 MD time-correlation functions). One paper
  fetched and summarised (Kahn et al. 2022, magnetizabilities —
  PMC8903098). Brain-susceptibility paper incidentally fetched and
  noted (PMC3254777). Landscape populated at 0.1–0.3 plus
  shape-of-haystack section.
- **2026-04-22 pm — second pass, user-raised methods/formalism
  queries.** Three further queries run (0.4 tensor formalism in
  geometric deep learning / tensor-as-array vs tensor-as-physicist;
  0.5 equivariant tensor time-series learning / anisotropic data in
  training; 0.6 SE(3) / Lie-group kernels in robotics and protein
  structure). All three returned substantial literature: query 0.4
  returned the canonical equivariant-ML framework we already depend
  on via e3nn but had not cited; 0.5 returned active in-progress
  work (the user predicted this); 0.6 returned the mature SE(3)
  robotics + recent protein-equivariant-GNN literature (the user
  called this correctly). Landscape updated to eight subsections;
  shape-of-haystack now identifies seven distinct research
  communities, six under-covered in current bibliography.
- **2026-04-22 pm — bibliography cross-check during pre-clear
  cleanup.** Read `references/ANNOTATED_BIBLIOGRAPHY.md` structure
  (879 lines, 14 core references, sections A-L). Discovered that
  several "thin coverage" claims in 0.4 and 0.6 were wrong —
  equivariant ML (E) and geometric descriptors (F) are strong
  sections of the existing bib. Corrected shape-of-haystack in 0.8.
  Actual under-covered communities are 3 of 7 genuinely thin
  (magnetic biophysics / QM magnetic properties / SE(3) robotics),
  2 of 7 partially thin (stat-mech TCF foundations, recent
  equivariant tensor time-series), 2 of 7 already well covered.
- **Remaining for Session 0.** Pending queries listed in 0.7 —
  seven physics haystacks plus three additional methods/formalism
  haystacks (Lie-group Gaussian processes, SE(3) diffusion, neural
  operators). Bibliography addenda to write, focused on the three
  genuinely under-covered communities (M: magnetic biophysics;
  B-extension or new section: QM magnetic-property-calculation at
  the magnetizabilities / NICS / aromaticity level; F-extension or
  new section: SE(3) robotics and recent protein-SE(3)-equivariant
  work). A checkpoint-and-resume in a fresh session is the right
  next move.
- **2026-04-22 pm — 0.7 pending queries close-out.** All seven
  physics haystacks + three methods haystacks from 0.7 addressed
  with tethered bibliography additions. Entries landed across
  four sections, every new entry carrying explicit back-references
  to the 13.2 handles table, T2-residual diagnostic (MATHS_GOALS
  pillar 2), and/or specific `TrajectoryResult` / rollup W-series
  field per the 13.4 posture directive.

  **Added this pass:**
  - Section L (CSA as a field): **L10-L12** — Sitkoff-Case 1998
    ab initio CSA review, Brender-Taylor-Ramamoorthy 2001 ¹⁵N
    CSA orientations, Loth-Pelupessy-Bodenhausen 2005 ubiquitin
    amide CSA from CCR (CORE, closest-scale external validation
    on our calibration-set protein).
  - Section H (MD convergence): **H15** — 2018 JPC B package-
    comparison "Validating MD vs Observables in Light of
    Ensembles" — load-bearing caveat for T2-residual interpretation.
  - Section M (Stat-mech foundations, community 4): **M31** —
    Green 1954 JCP 22, 398 — the "Green" in "Green-Kubo," pairs
    with M27 Kubo 1957. Embeds the **first explicit calculator
    hint**: a `GreenKuboSpectralDensityResult` producing J(ω) at
    Larmor frequencies as the Fourier dual of the σ-tensor
    autocorrelation TrajectoryResult field.
  - Section B (r²SCAN benchmarks): **B19-B20** — Grimme r²SCAN-3c
    composite (relevant to μs-harvester DFT budget) + 2025
    Faraday Discussions r²SCAN-for-shielding benchmark (grounds
    our calibration functional's performance regime).
  - Section J (ML for NMR): **J6** — 2024 JCTC DFT + 3D GNN
    hybrid, closest-prior-art on the Stage 3 architecture axis.
  - Section K (force fields / infra): **K5** — Řezáč-Hobza 2016
    Chem. Rev. on semiempirical methods, contexts K3 PM7 for
    charge quality feeding MopacCoulomb / MopacMcConnell.
  - Section E (Methods neighbours for Stage 3, 0.7 methods
    haystacks): **E24-E28** — Azangulov-Borovitskiy et al. 2024
    stationary kernels on Lie groups (Lie-group GPs; calculator
    hint for GP-on-SE(3) alternative to equivariant-GNN);
    Yim-Trippe et al. 2023 FrameDiff + Yim 2024 WIREs review
    (SE(3) diffusion); Li et al. 2021 FNO + Lu et al. 2021
    DeepONet (neural operators; calculator hint for FNO-
    volumetric-shielding-field extension of `QtBiotSavartCalc`
    closed-form math).

  **Calculator hints surfaced this pass (for user review):**
  1. **`GreenKuboSpectralDensityResult`** — J(ω) per atom at named
     Larmor frequencies, computed as the Fourier-transform dual of
     the planned σ(t) autocorrelation in TrajectoryResult. Feeds
     R1/R2/R1ρ predictions directly per 13.2. Close to free if
     the autocorrelation is already computed.
  2. **`PseudocontactShiftResult`** — same K_ab dipolar kernel as
     `RingSusceptibilityResult`, but with an external χ tensor
     (from paramagpy fit) and lanthanide-tag position as input.
     Directly validates against N5-N7 PCS datasets. Follows the
     `MopacCoulombResult` / `MopacMcConnellResult` "same kernel,
     different source data" precedent already in the
     GEOMETRIC_KERNEL_CATALOGUE.
  3. **Per-secondary-structure CSA stratification** — not a new
     calculator, a diagnostic analysis pass: our predicted σ
     tensor stratified by DSSP SS should reproduce the Yao-Bax
     2010 (L7) α-helix vs β-sheet ¹⁵N CSA split (~11 ppm). Thesis
     finding in its own right if reproduced; T2-residual flag if
     not.
  4. **Lie-group GP regression on SE(3)** (E24) — Stage 3 model
     architecture alternative to equivariant-GNN, with natural
     uncertainty quantification.
  5. **FNO for volumetric shielding field** (E27) — Stage 3
     alternative to per-atom regression if the viewer-rendered
     BS / HM butterfly isosurfaces are ever to be learned rather
     than computed closed-form.

  **Still open:** none of the 0.7 queries remain unaddressed.
  Citation polish pending on ~15 entries marked "to verify"
  (grep `to verify` in ANNOTATED_BIBLIOGRAPHY.md). Session 0
  landscape is effectively complete; Sessions 1-3 drafting can
  begin against this scaffold.
- **2026-04-22 pm — community 4 depth pass.** Chased stat-mech /
  MD time-correlation-function theory in depth per user request.
  Bibliography gains section M subsection "Stat-mech foundations
  (TCF, linear response, memory kernels — community 4)" with
  entries **M27-M30**: M27 Kubo 1957 (linear response + FDT,
  *J. Phys. Soc. Jpn.* 12, 570); M28 Berne-Harp 1970 (TCF
  review chapter, *Adv. Chem. Phys.* 17, 63-227 — canonical
  citation; the short companion at *Phys. Rev. A* 2, 975 also
  referenced inline in the entry); M29 Zwanzig 1961 primary
  memory-function paper (*Phys. Rev.* 124, 983 — the primary
  paper complementing the M22 Zwanzig 2001 book); M30
  Ditler-Luber 2022 (WIREs Comput. Mol. Sci. review of
  compute-observable-per-frame-then-autocorrelate workflow for
  vibrational spectroscopies — architectural sibling one step
  across from NMR). All four citations DOI-verified by web search
  during the pass. Community 4 is now closed; no remaining
  genuinely under-covered communities from the 0.8 list.
  Remaining Session 0 moves: the 0.7 pending queries
  (paramagnetism in biomolecules, GIAO alternatives survey,
  MD-derived property calculation, semiempirical magnetic,
  ensemble averaging, CSA as a field, classical predictors as
  a research area, Lie-group GPs, SE(3) diffusion, neural
  operators).
- **2026-04-22 pm — community 6 depth pass.** Chased "equivariant
  tensor time-series learning" in depth per user request.
  Bibliography gains section E subsection **E15-E23** —
  "Equivariant tensor-valued and time-series prediction
  (community 6)". Breakdown: **surveys** (E15 Kondor 2025 PNAS
  pedagogical — the misattribution flagged in 0.4 corrected
  during this pass; E16 Gerken 2023 Artif. Intell. Rev.; E17
  Fei & Deng 2024 3D equivariance survey); **tensor-message-
  passing architectures** (E18 HotPP — Wang et al. 2024
  Nat. Commun. Cartesian tensor messages; E19 TSENN — Hsu et al.
  2026 Nat. Commun. tensorial spectra, confirmed as the
  "anisotropic-sequential-tensors approach" referred to in 0.5;
  E20 TEGNN 2024 with citation details still to verify);
  **trajectory / time-series architectures** (E21 EGNO — Xu et al.
  ICML 2024 neural operators over time, Fourier-space equivariant
  temporal conv; E22 TrajCast — Thiemann et al. 2025 autoregressive
  equivariant MD; E23 Gregory-Villar et al. 2026 ICLR — tensor
  time series via path signatures, retitled from the 0.5 arXiv
  entry). Several entries carry "to verify" flags for
  volume/author-list polish; each flag is inline in the entry.
  Remaining Session 0 moves: the 0.7 pending queries (paramagnetism
  in biomolecules, GIAO alternatives survey, MD-derived property
  calculation, semiempirical magnetic, ensemble averaging, CSA as
  a field, classical predictors as a research area, Lie-group GPs,
  SE(3) diffusion, neural operators).
- **2026-04-22 pm — bibliography round-out (Communities 3 and 7).**
  Three genuinely under-covered communities flagged in 0.8 now have
  starting-place bibliography entries: community 2 (magnetic
  biophysics) was already covered by section N in the earlier pass;
  community 3 (QM magnetic-property-calculation at the
  magnetizabilities / NICS / aromaticity level) now has section B
  extension **B14-B18** — new subsection "Magnetic response beyond
  shielding: magnetizabilities, NICS, aromaticity" (Schleyer 1996
  NICS, Gershoni-Poranne & Stanger 2015, Kahn 2022, Helgaker-Coriani
  2012, Flygare-Benson 1971); community 7 (SE(3) robotics and
  recent protein-SE(3)-equivariant) now has section E extension
  **E10-E14** — new subsection "Recent protein-SE(3)-equivariant
  applications" (Satorras EGNN 2021, EquiCPI 2025, VN-EGNN, SE(3)
  ternary complex 2025, Krapp PPI 2023) plus section F extension
  **F6-F7** — new subsection "SE(3) and Lie-group formalism
  (robotics crossover)" (Chirikjian two-volume, Helgason).
  Citation polish pending: several entries carry "DOI to verify"
  / "author list to verify" flags where they were added from
  Session 0.6 mentions rather than from fetched PDFs. Organisation
  work (community 6 equivariant tensor time-series — partial — and
  the 0.7 pending queries) is the next Session 0 move.
