# Entry options menu — annotated bibliography

For each of the 8 papers below: a single Sentence 1 (system + method, mostly fixed by what the paper is) and menus for Sentence 2 (result) and Sentence 3 (authors' interpretation), plus short relevance-tag options. Pick, recombine, rewrite. Sentence forms follow the rubric example: long compound sentences in plain English, no quoted equations, no cross-paper synthesis.

---

## 1. Markwick et al. 2010

**Citation.** Markwick, P.R.L., Cervantes, C.F., Abel, B.L., Komives, E.A., Blackledge, M. & McCammon, J.A. (2010) Enhanced conformational space sampling improves the prediction of chemical shifts in proteins. *Journal of the American Chemical Society* 132(4), 1220–1221.

**Source label.** Primary source.

**Header (one-line topic) — options.**
- (A) Evidence that averaging predicted chemical shifts over a sampled molecular dynamics ensemble agrees with experiment better than predicting from a single static structure.
- (B) Evidence that conformational sampling improves chemical-shift prediction in a real protein, with the gain concentrated in regions that are independently shown to be dynamic.

**Sentence 1 (system + method).**
This study examines IκBα, an ankyrin-repeat protein with regions of slow conformational exchange, by computing backbone chemical shifts for each frame of an accelerated molecular dynamics ensemble using the empirical SHIFTX predictor and averaging across the trajectory.

**Sentence 2 (result) — options.**
- (A) Chemical-shift root-mean-square deviation against experiment fell from 8.19 ppm for the X-ray structure to 5.92 ppm for the optimal accelerated trajectory, with the largest single-channel improvement on backbone nitrogen, and the spatial pattern of improvement matched independently measured regions of millisecond exchange.
- (B) The crystal structure gave a chemical-shift error around 8 ppm; plain molecular dynamics improved this slightly; accelerated dynamics — which biases dihedral-angle sampling to escape low-energy basins faster — brought the error down to about 6 ppm, with the largest improvement on backbone nitrogen and on residues already known to be dynamic from spin-relaxation measurements.
- (C) The improvement in shift agreement was largest in regions of the protein where independent relaxation measurements showed slow conformational exchange, and was negligible in a comparison protein (ubiquitin) which lacks such regions.

**Sentence 3 (interpretation) — options.**
- (A) It is proposed that conformational sampling, rather than improvements in the chemical-shift predictor itself, accounts for the gain in agreement, since the empirical predictor and the protein were held fixed across the comparison.
- (B) The authors note that no single sampled frame gave as good a chemical-shift prediction as the trajectory average, indicating that the gain comes from the averaging rather than from finding a better representative structure.
- (C) It is concluded that chemical shifts in proteins should be treated as ensemble averages rather than single-structure properties when comparing to experiment.

**Relevance tag — options.**
- (A) Relevant to molecular-dynamics-based prediction of NMR observables in proteins.
- (B) Bears on whether conformational sampling is needed for protein chemical-shift prediction.

---

## 2. Yates and Bartók 2024

**Citation.** Yates, J.R. and Bartók, A.P. (2024) Accurate predictions of chemical shifts with the rSCAN and r²SCAN mGGA exchange–correlation functionals. *Faraday Discussions* 255, 192–202.

**Source label.** Primary source.

**Header — options.**
- (A) Evidence that the rSCAN and r²SCAN meta-GGA functionals predict NMR chemical shifts substantially more accurately than the older PBE functional, at modest extra computational cost.
- (B) Direct benchmark of the rSCAN and r²SCAN functionals for NMR shielding prediction against experimental chemical shifts of inorganic compounds.

**Sentence 1 (system + method).**
This study benchmarks two recent meta-GGA exchange–correlation functionals — rSCAN and r²SCAN — by computing nuclear magnetic shieldings for a set of inorganic halide and oxide compounds within the GIPAW pseudopotential framework, and comparing the predicted chemical shifts to experimentally measured values.

**Sentence 2 (result) — options.**
- (A) The rSCAN and r²SCAN functionals reduce root-mean-square errors against experiment to roughly half the values obtained with the older PBE generalised-gradient functional, at modest additional cost (about twice the time of PBE for a small unit cell, dropping to 1.25 times for a larger one).
- (B) The authors report that r²SCAN gives results comparable in accuracy to rSCAN but is more numerically stable to changes in the integration grid, making it the more practical choice of the two for routine use.
- (C) The correlation between predicted and experimentally measured chemical shifts approaches the theoretical maximum, with substantially smaller deviations than were achievable with the previous generation of functionals.

**Sentence 3 (interpretation) — options.**
- (A) It is concluded that the rSCAN and r²SCAN functionals offer a clear improvement in chemical-shift accuracy over the widely-used PBE functional at modest additional cost, supporting their adoption for routine NMR shielding calculations.
- (B) The authors argue that the improvement carries through to spectrum assignment in experimental investigations, where the smaller deviation from observed shifts allows more reliable identification of resonances in complex compounds.
- (C) It is shown that the meta-GGA functionals close much of the accuracy gap between generalised-gradient methods and more expensive hybrid functionals for NMR shielding, without the cost of exact-exchange evaluation.

**Relevance tag — options.**
- (A) Relevant to the choice of exchange–correlation functional in DFT chemical-shift calculations.
- (B) Bears on which density functional should be used for NMR shielding calculations targeting experimental accuracy.

---

## 3. Boyd and Skrynnikov 2002

**Citation.** Boyd, J. and Skrynnikov, N.R. (2002) Calculations of the contribution of ring currents to the chemical shielding anisotropy. *Journal of the American Chemical Society* 124(9), 1832–1833.

**Source label.** Primary source.

**Header — options.**
- (A) Evidence that the Johnson–Bovey ring-current model can be extended from the isotropic chemical shift to the full rank-2 shielding tensor, with the additional tensor component validated against full DFT calculation.
- (B) Extension of the classical Johnson–Bovey ring-current formula from a scalar isotropic shift to the complete rank-2 nuclear magnetic shielding tensor.

**Sentence 1 (system + method).**
This study derives a closed-form expression for the previously-missing tensor component of the ring-current shielding contribution, completing the rank-2 form of the Johnson–Bovey loop-current model, and tests the result by comparing the analytical prediction against full DFT shielding calculations on an N-methyl acetamide probe placed near a benzene ring.

**Sentence 2 (result) — options.**
- (A) The combined analytical expression — local shielding from the probe placed near a non-aromatic analogue, plus the new ring-current tensor formula — closely tracks the full DFT calculation across rotation of the probe, with isotropic-shift agreement to roughly 1 ppm and matching tensor anisotropy.
- (B) Application to the fibronectin type-2 module backbone amide proton of residue G42, which forms a hydrogen bond into the face of a phenylalanine ring, gives a chemical-shielding-anisotropy contribution of about 17 ppm from the ring current alone, in agreement with full-fragment DFT.
- (C) For backbone carbonyl carbon and amide nitrogen sites, the ring-current contribution to the shielding anisotropy is reported as typically less than 1 ppm, indicating that the ring current is a major contributor to anisotropy only at proton sites in close contact with aromatic rings.

**Sentence 3 (interpretation) — options.**
- (A) It is shown that the Johnson–Bovey ring-current model can be extended to predict the full rank-2 shielding tensor, and that the analytical formula reproduces full DFT to within roughly 1 ppm at proton sites near aromatic rings.
- (B) The authors argue that ring-current contributions to chemical shielding anisotropy are negligible for backbone heavy nuclei but can dominate for protons in N–H⋯π contacts, with implications for the interpretation of relaxation experiments in proteins.
- (C) It is proposed that local-environment shielding and ring-current contributions can be separated by comparing model fragments with and without the aromatic ring, and combined additively to recover the full tensor.

**Relevance tag — options.**
- (A) Relevant to extending classical ring-current kernels from the isotropic shift to the full shielding tensor.
- (B) Bears on whether ring-current models can predict tensor components rather than only the isotropic shift.

---

## 4. Case 1995

**Citation.** Case, D.A. (1995) Calibration of ring-current effects in proteins and nucleic acids. *Journal of Biomolecular NMR* 6, 341–346.

**Source label.** Primary source.

**Header — options.**
- (A) Calibrates two empirical ring-current models and a Buckingham-style electric-field term against direct DFT shielding calculations on a model system relevant to proteins and nucleic acids.
- (B) Direct DFT calibration of the empirical ring-current and electric-field formulas used in protein and nucleic acid chemical-shift prediction.

**Sentence 1 (system + method).**
This study computes the chemical shielding of a methane probe placed at hundreds of positions in and around ten aromatic ring systems found in proteins and nucleic acids — benzene, indole, imidazole, the nucleobases, and others — using density-functional theory, and fits the resulting shieldings to two empirical ring-current models (Johnson-Bovey loop currents and Haigh-Mallion's surface-integral expression) together with a separate Buckingham-style electric-field term.

**Sentence 2 (result) — options.**
- (A) Both ring-current models reproduce the calculated shieldings with correlation coefficients of 0.94 to 0.99 and root-mean-square differences below 0.1 ppm; the two are statistically indistinguishable for predicting the isotropic chemical shift.
- (B) The fit gives updated ring-current intensity factors that are substantially larger than the previously used Giessner-Prettre and Pullman values for histidine and tryptophan, and an electric-field coefficient consistent with Buckingham's earlier theoretical estimate.
- (C) The fits show that the Johnson-Bovey and Haigh-Mallion models agree to within 0.1 ppm root-mean-square difference for the isotropic shift, despite using different physical pictures of the ring current.

**Sentence 3 (interpretation) — options.**
- (A) It is proposed that the updated ring-current intensity factors should replace the older Giessner-Prettre and Pullman values used in protein chemical-shift prediction software.
- (B) The author argues that the choice between Johnson-Bovey and Haigh-Mallion is a matter of convenience rather than accuracy for predicting isotropic shifts in proteins, since both fit the calculated DFT values equally well.
- (C) It is shown that an electric-field correction with a coefficient near three parts per million per atomic unit accounts for polarisation effects neglected by the ring-current term alone.

**Relevance tag — options.**
- (A) Relevant to the calibration of geometric chemical-shift kernels against quantum-mechanical reference calculations.
- (B) Bears directly on the empirical formulas used in protein and nucleic acid chemical-shift prediction software.

---

## 5. Kellner et al. 2025

**Citation.** Kellner, M., Holmes, J.B., Rodriguez-Madrid, R., Viscosi, F., Zhang, Y., Emsley, L. and Ceriotti, M. (2025) A deep learning model for chemical shieldings in molecular organic solids including anisotropy. *Journal of Physical Chemistry Letters*, 8714–8722.

**Source label.** Primary source.

**Header — options.**
- (A) Evidence that a machine-learning model can predict the full anisotropic chemical shielding tensor of atoms in molecular organic solids to an accuracy approaching that of the underlying DFT reference calculations.
- (B) Demonstration of a machine-learning model that predicts the complete shielding tensor — not just the isotropic chemical shift — for atoms in organic solids.

**Sentence 1 (system + method).**
This study introduces ShiftML3.0, a deep-learning model that predicts nuclear magnetic chemical shieldings for atoms in molecular organic solids, including the full anisotropic shielding tensor — extending the previous ShiftML and ShiftML2 models which were trained to predict only the isotropic chemical shift.

**Sentence 2 (result) — options.**
- (A) The model achieves root-mean-square errors against experiment of approximately 0.5, 2.4, and 7.2 parts per million for proton, carbon-13, and nitrogen-15 nuclei respectively — close to the errors of the underlying density-functional reference calculations themselves (0.5, 2.3, 5.8 ppm).
- (B) By predicting the full shielding tensor rather than the isotropic value alone, the model recovers tensor anisotropy and orientation with comparable accuracy, allowing simulation of solid-state NMR observables that depend on tensor structure rather than only on the trace.
- (C) The authors report that the model's accuracy approaches that of the density-functional reference itself across all three nuclei tested, indicating that further improvement is now limited by the reference calculations rather than the machine-learning architecture.

**Sentence 3 (interpretation) — options.**
- (A) It is argued that the gap between machine-learning-predicted and density-functional-computed shieldings has now closed to within the accuracy of the underlying functional itself, and that further progress requires improving the reference calculations rather than the machine-learning method.
- (B) The authors argue that tensor-output prediction — rather than isotropic-shift prediction alone — is necessary for simulating the observables of solid-state NMR experiments that resolve anisotropy, and that machine learning can now provide this at low computational cost.
- (C) It is proposed that the same architecture can be extended to other anisotropic NMR observables, such as electric field gradients, by making the equivariant network output the appropriate irreducible-representation components.

**Relevance tag — options.**
- (A) Relevant to the prediction of full shielding tensors, rather than isotropic shifts alone, by machine-learning methods.
- (B) Bears on whether machine-learning architectures can recover the tensor structure of nuclear magnetic shielding to chemical accuracy.

---

## 6. Han et al. 2011 (SHIFTX2)

**Citation.** Han, B., Liu, Y., Ginzinger, S.W. and Wishart, D.S. (2011) SHIFTX2: significantly improved protein chemical shift prediction. *Journal of Biomolecular NMR* 50(1), 43–57.

**Source label.** Secondary source (bioinformatics / software methods). [Reasoning: SHIFTX2 integrates existing methods (SHIFTX, SHIFTY, BLAST, REPTree via WEKA) with a hand-curated training set; the work is software/methods integration rather than physical-chemistry derivation or measurement, so secondary in a physical-chemistry-discipline sense even though it is original peer-reviewed research in the bioinformatics sense.]

**Header — options.**
- (A) Evidence that a hybrid sequence-and-structure machine-learning predictor reaches substantially higher accuracy than its predecessors for protein backbone chemical shifts.
- (B) Description and benchmark of SHIFTX2, a hybrid machine-learning predictor of protein chemical shifts.

**Sentence 1 (system + method).**
This study introduces SHIFTX2, a chemical-shift prediction program for proteins that combines a sequence-based component (which aligns the query against a database of proteins with known shifts) with a structure-based machine-learning component (boosted regression trees over features including backbone and side-chain torsion angles, solvent accessibility, hydrogen-bond geometry, pH, and temperature).

**Sentence 2 (result) — options.**
- (A) On a benchmark of high-resolution X-ray structures with assigned chemical shifts, SHIFTX2 reaches correlation coefficients between predicted and observed backbone shifts of roughly 0.97 to 0.999 across the six backbone nuclei, with root-mean-square errors ranging from about 0.12 ppm for protons up to about 1.1 ppm for nitrogen.
- (B) The authors report that SHIFTX2 is substantially more accurate than the next-best protein chemical-shift predictor — up to 26% better in correlation coefficient and with root-mean-square errors several times smaller — and is also faster and covers a wider variety of shift types.
- (C) The hybrid combination of sequence-based and structure-based sub-modules outperforms either component on its own across all six backbone shift channels.

**Sentence 3 (interpretation) — options.**
- (A) It is argued that the improvement comes from three sources together: a larger and higher-quality training set of around 200 proteins, additional structural features such as hydrogen-bond geometry and solvent accessibility, and the combination of sequence-based and structure-based methods rather than reliance on either alone.
- (B) The authors propose that this accuracy now makes chemical-shift prediction practically useful for protein structure determination, refinement, and validation — applications that had previously been limited by predictor accuracy.
- (C) It is concluded that hybrid predictors of this form set the practical accuracy bar that subsequent protein chemical-shift methods must clear.

**Relevance tag — options.**
- (A) Relevant to the empirical / machine-learning predictor lineage that physics-kernel approaches contrast against.
- (B) Bears on what level of accuracy is currently achievable in protein chemical-shift prediction by methods other than physics-based extraction.

---

## 7. Sahakyan & Vendruscolo 2013

**Citation.** Sahakyan, A.B. & Vendruscolo, M. (2013) Analysis of the contributions of ring current and electric field effects to the chemical shifts of RNA bases. *The Journal of Physical Chemistry B* 117, 1989–1998.

**Source label.** Primary source.

**Header — options.**
- (A) Compares three ring-current models and a separate electric-field term against DFT shielding calculations on RNA bases.
- (B) Calibrates ring-current and electric-field contributions to chemical shifts for the four nuclei observed in RNA bases.

**Sentence 1 (system + method).**
This study builds a database of inter-base geometries from RNA crystal structures, computes chemical shielding for each arrangement using density-functional theory, and fits the calculated shifts to three ring-current models — Pople's point-dipole, Johnson-Bovey's loop currents, and Haigh-Mallion's surface integral — together with a Buckingham-style electric-field term, separately for proton, carbon-13, nitrogen-15, and oxygen-17 nuclei.

**Sentence 2 (result) — options.**
- (A) For protons, the combined ring-current and electric-field model reaches a correlation of 0.96 against the calculated shifts; for the heavier nuclei, the electric-field contribution dominates the ring-current contribution, with correlations from electric field alone exceeding those from ring current alone for carbon and nitrogen.
- (B) The three ring-current models give similar correlations for protons but diverge in their geometric factors, and the authors derive empirical conversion equations between the simpler Pople formula and the more accurate Haigh-Mallion expression.
- (C) The effect of an adjacent ring on the chemical shift of a nucleus was found to reach approximately 5, 35, and 27 parts per million for proton, carbon-13, and nitrogen-15 respectively at hydrogen-bond distances, with electric-field effects on oxygen-17 reaching 80 ppm.

**Sentence 3 (interpretation) — options.**
- (A) It is shown that ring-current contributions dominate the chemical shifts of light nuclei, while electric-field contributions dominate for heavier nuclei — a hierarchy that should be respected when calibrating empirical predictors.
- (B) The authors argue that the Pople point-dipole model, when corrected by their derived conversion factors, can recover Haigh-Mallion-level accuracy at lower computational cost.
- (C) It is proposed that the standard practice of computing a ring's normal vector from three ring atoms is unstable to out-of-plane atomic fluctuations, and that averaging two independently computed normals gives a more reliable estimate.

**Relevance tag — options.**
- (A) Relevant to per-nucleus calibration of ring-current and electric-field kernels in biomolecular chemical-shift prediction.
- (B) Bears on which physical contributions dominate the chemical shift of each kind of nucleus.

---

## 8. Klukowski, Riek and Güntert 2025

**Citation.** Klukowski, P., Riek, R. and Güntert, P. (2025) Machine learning in NMR spectroscopy. *Progress in Nuclear Magnetic Resonance Spectroscopy* 148–149, 101575.

**Source label.** Secondary source. (Review.)

**Header — options.**
- (A) Comprehensive review of machine-learning applications across NMR spectroscopy, covering experimental-data analysis, structure-based property prediction, and spectrum reconstruction.
- (B) Recent review of how machine learning has been applied across NMR — from peak picking and assignment to chemical-shift prediction and data reconstruction.

**Sentence 1 (system + method).**
This review surveys machine-learning applications across nuclear magnetic resonance spectroscopy, covering automated analysis of experimental spectra, prediction of chemical shifts and other observables from molecular structure, and reconstruction or denoising of measured data, with an audience of NMR practitioners new to machine learning.

**Sentence 2 (result) — options.**
- (A) The review compiles a chronological map of methods from 2000 through 2025 and a tabular comparison of structure-based chemical-shift predictors, finding that regression-based methods such as SHIFTX2 and UCBShift remain competitive with deep-learning methods on protein backbones, while graph neural networks operating directly on molecular geometry offer the more promising direction for generalisation beyond hand-crafted features.
- (B) Across more than two decades of methods, the review documents a persistent data-scarcity problem in NMR — multi-dimensional time-domain spectra are deposited for fewer than 1.5% of database entries — and notes that this constraint, rather than algorithmic limitations, has slowed the adoption of large-scale deep-learning approaches.
- (C) The review reports that for protein chemical shifts, deep-learning methods have not yet displaced regression-based predictors trained on hand-crafted features, but that solid-state methods such as ShiftML have been extended to molecular-dynamics snapshots — establishing a precedent for trajectory-based shift prediction.

**Sentence 3 (interpretation) — options.**
- (A) It is argued that data scarcity in NMR — far below the scale that has driven progress in image and language tasks — makes novel-angle contributions, rather than scaled-up training, the realistic route forward for the field.
- (B) The authors argue that machine learning's role in NMR is most useful as an integration engine for the many geometric and spectroscopic probes a single molecule provides, rather than as a replacement for physics-based predictors.
- (C) It is concluded that learning directly on molecular-graph geometry — rather than on hand-crafted feature sets — gives the best route to predictors that generalise beyond the proteins seen at training time.

**Relevance tag — options.**
- (A) Relevant to placing classical chemical-shift prediction in the current machine-learning landscape.
- (B) Bears on what role machine learning currently plays in the prediction of NMR observables.

---

## Notes on assembly

- **Word budget per entry.** Citation does not count. Body should be ~120–150 words: header (~15–25 words) + Sentence 1 (~30–45) + Sentence 2 (~30–45) + Sentence 3 (~25–35) + relevance tag (~10–15). The example sits at 119 words; aim near the floor.
- **Source label.** Bare line, just `Primary source.` or `Secondary source.` Full stop. The header sentence follows on the same line per the example.
- **Tense.** Past tense for what was done; present tense for what is proposed or argued.
- **No equations, no acronyms unexplained at first use, no proper-noun terms-of-art that aren't.**
- **Five primary required.** Final eight (alphabetical): Boyd-Skrynnikov (P), Case (P), Han (S — bioinformatics methods), Kellner (P), Klukowski (S, advisor review), Markwick (P, advisor), Sahakyan-V (P), Yates-Bartók (P) → 6 P / 2 S. Compliant.
