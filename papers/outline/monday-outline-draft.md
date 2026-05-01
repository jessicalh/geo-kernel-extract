# Literature-Review Outline — Draft for Monday submission

**Author:** Jessica Hansberry
**Date:** 2026-04-27

Five topical sections. Each opens with a brief orientation paragraph, lists the relevant works from the accompanying annotated bibliography under a "Key references" heading, and sketches the argument flow in bullets. Placeholders marked `[...]` are points for the author to fill in or extend with additional material.

---

## 1. NMR Physics

*Framing paragraph.* Nuclear magnetic resonance reports on local electronic environment through several distinct observables — chemical shifts, scalar couplings, dipolar couplings, relaxation rates — each with its own sensitivity to structure and dynamics. For computational work that integrates molecular-dynamics simulation with NMR measurement, the central physical fact is that different observables average over different timescales. Chemical-shift averaging occurs over microseconds to milliseconds through exchange coalescence and scales with the magnetic field, while spin relaxation reports on picosecond-to-nanosecond motions. Conflating these regimes is a common pitfall; distinguishing them is a prerequisite for any argument about what MD trajectories can and cannot predict about NMR measurements.

**Key references**
- Bryant 1983 — "The NMR Time Scale." Timescale framework for each observable.
- Facelli 2011 — bridges from observable to theory (also under §3).

**Argument flow**
- NMR observables as indirect probes of atomic-level structure.
- The distinct averaging windows of shift, coupling, and relaxation (Bryant 1983).
- Implications for MD-simulation-based NMR prediction: what the sampled trajectory is averaging and at what effective bandwidth.
- [Author note: one additional sentence on the motivation for working with shielding tensors rather than scalar shifts goes here, framed as "because the tensor carries structural information the scalar trace averages over."]

---

## 2. NMR and Machine Learning

*Framing paragraph.* Machine-learning applications in NMR have grown rapidly since the 2010s, encompassing automated spectrum analysis (peak picking, chemical-shift assignment, structure determination), structure-to-property prediction for proteins, small molecules, and solids, and spectrum reconstruction and quality enhancement. For work at the intersection of geometric prediction and tensor-valued observables, the key theoretical development is the formal framework of equivariant neural networks, which constrains the functional form of any model whose inputs and outputs are required to transform predictably under rotation, translation, or permutation.

**Key references**
- Kondor 2025 — equivariant-neural-network principles for physics and chemistry.
- [Author note: if the program's rubric has room, a recent review such as Klukowski, Riek & Güntert (2025) *Prog. NMR Spectrosc.* 148–149, 101575 gives the application-side landscape (SHIFTX, SPARTA+, UCBShift, ShiftML, and related). That paper is on disk if needed.]

**Argument flow**
- Why chemical-shift prediction is a mature applied-ML problem for proteins but data-scarce compared to structure prediction.
- Evolution from empirical regression predictors (PROSHIFT, SPARTA+, SHIFTX lineage) through DFT-database approaches (PROCS15, ShiftML) to modern graph and equivariant networks.
- Kondor 2025 as the design-space reference: Peter-Weyl decomposition into irreducible representations, Clebsch-Gordan product as the canonical polynomial equivariant nonlinearity.
- Why tensor-output protein shift prediction remains open territory, even given the full ML landscape above.

---

## 3. Shielding Physics

*Framing paragraph.* The nuclear magnetic shielding tensor σ_αβ is a second-derivative response of the electronic energy to the nuclear magnetic moment and the external magnetic field. Its value at a given nucleus depends on the full electronic structure of the molecule plus any external electric and magnetic perturbations. Classical decompositions into contributions from ring currents, local bond anisotropies, through-space electric-field effects, and solvent reaction fields provide the physical vocabulary that quantitative shift-prediction work inherits. The tensor nature of the observable matters: isotropic-trace agreement can hide tensor-component disagreement, so tensor-level comparison is the more rigorous validation.

**Key references**
- Buckingham 1960 — linear electric-field contribution to shielding in polar-group environments.
- Facelli 2011 — modern review of shielding-tensor theory and computation.

**Argument flow**
- Shielding-tensor definition and irreducible-representation decomposition (T0 = isotropic trace; T1 = antisymmetric part; T2 = traceless symmetric part, containing the anisotropy).
- Classical decomposition of the shielding contribution by physical mechanism (ring current, electric field, bond anisotropy, paramagnetic, diamagnetic).
- Buckingham 1960's linear-in-E term as the archetype for classical electric-field modelling and the numerical anchor for its coefficient.
- Facelli 2011's argument that tensor-component comparison is diagnostically more rigorous than isotropic-trace comparison, because of fortuitous error cancellation in the trace.
- Expected accuracies for DFT calculation of shielding (¹H 2–3%, ¹³C 1.4–1.9%, ¹⁵N approximately 10%, ¹⁷O approximately 14%).

---

## 4. Geometric Kernels

*Framing paragraph.* Geometric kernels are closed-form functions of atomic geometry that estimate specific physical contributions to nuclear magnetic shielding without requiring per-evaluation quantum-mechanical calculation. Ring-current models in particular (Pople point-dipole, Waugh-Fessenden-Johnson-Bovey loop-current, Haigh-Mallion Hückel surface-integral) have been in use since the 1950s and can be combined with Buckingham-style electric-field terms to capture the dominant polar-group contributions. Calibration of these kernels against density-functional-theory shieldings places them on an accuracy footing sufficient for biomolecular shift prediction, and in some regimes the classical kernels are competitive with, or comparable to, empirical database predictors.

**Key references**
- Case 1995 — DFT calibration of ring-current and electric-field kernels for proteins and nucleic acids.
- Sahakyan & Vendruscolo 2013 — extension and comparison across three ring-current models in RNA bases.
- Gershoni-Poranne & Stanger 2015 — magnetic aromaticity criteria, including the NICS_zz tensor-component argument that parallels this thesis's T2 focus.

**Argument flow**
- Historical lineage: Pople 1956 → Johnson-Bovey 1958 → Haigh-Mallion 1980.
- Case 1995's DFT-calibrated intensity factors as the standard reference values for biomolecular work; the equal empirical performance of Johnson-Bovey and Haigh-Mallion for scalar shifts.
- Sahakyan-Vendruscolo 2013's RNA extension; the finding that electric-field contributions dominate ring-current contributions for heavy nuclei.
- Gershoni-Poranne 2015's NICS_zz argument as an independent precedent, in the aromaticity subfield, for tensor-component analysis being more informative than isotropic-trace analysis.
- [Author note: the novel claim the thesis argues into is that these classical kernels, when evaluated per-frame on MD trajectories and calibrated against DFT shielding tensors, correlate with the DFT signal at the tensor-component level — not only on the isotropic trace. One or two lines setting that up here.]

---

## 5. Physics Extractions from MD Trajectories

*Framing paragraph.* Predicting NMR observables from molecular simulation requires treating the computed observable as a dynamical quantity rather than a single-structure property. For chemical shifts, this means averaging computed per-frame shieldings over conformational substates sampled by the trajectory; for tensor observables, it means preserving the tensor structure through the averaging. Published work demonstrates that microsecond-scale classical molecular dynamics, together with picosecond-scale first-principles molecular dynamics to capture fast bond vibrations, measurably improves the agreement of DFT-computed chemical shifts with experimental solid-state measurements. Extending this approach from small peptides in the solid state to proteins in explicit solvent is the operational scope of physics-extraction work in this thesis.

**Key references**
- de Gortari et al. 2010 — DFT + MD averaging for MLF tripeptide in the solid state; direct methodological ancestor for MD-averaged shift prediction.
- Sahakyan & Vendruscolo 2013 — multi-kernel calibration in structural contexts comparable to trajectory analysis (also under §4).

**Argument flow**
- The operational scope: per-frame extraction of computed observables over an MD trajectory, with the trajectory sampling conformational substates below the NMR coalescence limit (links back to §1, Bryant 1983).
- de Gortari 2010 as the quantitative demonstration: approximately 1 ppm improvement in ¹³C RMSD from microsecond-scale averaging on a tripeptide, against first-principles and classical-MD backgrounds.
- The averaging framework is conservative: arithmetic mean over Boltzmann-distributed frames (when classical MD samples are converged), with block-averaging diagnostics establishing the independent-sample count.
- The claim shape of this thesis's physics-extraction work is that per-frame classical-kernel outputs correlate with DFT-computed shieldings across several dimensions — mean, fluctuation amplitude, principal-mode structure — not that they match pointwise. Correlation at the r² level is the validation currency; the kernel extractors are demonstrating that they saw the signal, not that they reproduce quantum-mechanical values.
- [Author note: paragraph on where this thesis extends beyond the published ancestors (DFT-grounded shielding on real proteins in explicit solvent rather than tripeptide in solid state) fits here.]

---

## Scope and Limitations

[Author note: brief paragraph acknowledging what the literature does not yet cover and where this thesis sits. Suggested content: the gap between existing ML-based chemical-shift predictors (which predict isotropic shifts as scalars, from static structures) and DFT-grounded tensor calculations on MD-averaged protein geometries; the data-scarcity constraint in NMR that precludes AlphaFold-scale ML approaches; the novelty of correlating classical-kernel extractions with DFT-computed shieldings on a per-frame basis at protein scale.]

---

## Writing Conventions for This Review

Throughout this review the comparison between classical kernel predictions and DFT-computed shieldings is characterised as **correlation** rather than **match**. The evaluation currency is r² (or Pearson correlation) sufficient to demonstrate that the kernel extractors have captured real physical signal, rather than pointwise numerical agreement on shielding values. This reflects the thesis's actual epistemic claim: that a deep classical extraction from atomic geometry reproduces the signal-level content of quantum-mechanical shieldings, not their absolute accuracy.
