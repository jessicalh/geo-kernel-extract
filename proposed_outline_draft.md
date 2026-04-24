# proposed_outline_draft

**What this is:** the working lit-review / thesis-supporting-paper outline. A long-lived document, not the Monday program-rubric submission. Organises our reading and provides the scaffolding for the eventual thesis-supporting paper. Grows as we read more material; "good ideas worth adding" from reference sessions get folded here.

**What this is NOT:** Monday's program-rubric submission. That submission will be constrained to whatever format the program requires (8-item annotated bib already drafted; lit-review outline bounded by rubric). Monday's submission draws from this document but will be shorter and more structured to rubric. See final section for that distinction.

**Legend:**
- ✓ — PDF on disk, summary in `references-meta/`
- ⧖ — acquisition pending, see `references/PENDING_ACQUISITIONS.md`

**Last meaningful update:** 2026-04-24.

---

## Positioning

**Landscape reality.** Protein chemical-shift prediction is, from a working perspective, solved: SHIFTX and its descendants give useful numbers for structural biology work, and the community treats the problem as handled. Papers nevertheless keep appearing, for three reasons: (a) AlphaFold made ML-for-structural-biology high-status and many groups want to replicate that story in adjacent sub-fields; (b) students (including this one) want ML experience alongside bio; (c) it is a tidy boxed area of inquiry. This landscape has an uncomfortable structural feature: big-data ML, the approach that drove AlphaFold's success, is not available here — even using every entry in the BMRB there are not enough NMR records to train an AlphaFold-scale or LLM-scale system. What counts as a contribution is therefore not "can we match SHIFTX" (matched, by many groups) nor "can we outscale them with data" (cannot), but **a genuinely novel angle**.

**Our novel angle, three-fold (the thesis claim).**
1. Deep classical extraction of per-atom shielding information from scalar through T2, quickly, on frames — with DFT-validated signal and mathematically-validated dimensional structure.
2. Relating those extractions back to physics — what each kernel measures, what T2 components represent, how the Stage 1 weights partition across kernel families.
3. Pushing the resulting time-series data at ML, alongside a traditional GNN for the same problem. The ML backing is ESM-2-style treatment of OpenFold3 embeddings, *not* ESM-2 itself (advisor rejected ESM-2 on two grounds: airless sequence-only modelling without geometry, and failure to use mmseqs2 for test/train partitioning in the original paper).

**Two working advantages that let us occupy that angle.**

*Small / exploratory.* This is a 1-year MSc project, so the bar for publication-grade completeness is set to accept calculated-risk directions that produce partial but carefully-reported results. That licenses methods the mainstream NMR literature treats as too rigorous or borrowed-from-elsewhere — memory-kernel decomposition, Markov State Models, tensor factorisation, state-space / dynamic-factor modelling, coupled-mode time-series analysis. Our downside is bounded; the upside is a thesis chapter with novel methodology.

*Large / tooled-up.* Very few people come at this combination of data and methodology with our specific kit: a C++ library that produces per-frame kernels cheaply, an already-replicated ESM-2 pipeline being ported to OpenFold3 embeddings on DGX Spark hardware, modern equivariant / tensor-factorisation tools that have stabilised in the last two to three years, and the engineering discipline to drop new kernels in quickly and re-run. Data generation keeps pace with analysis; methods whose empirical payoff has historically been out of reach for small groups are within reach for us.

**Calibration on "low bar."** The advisor's framing is sincere in the sense that exploratory first-draft work is welcomed. It is not sincere in the sense of "anything goes." When someone working ten steps away from Rosalind Franklin's old lab sets a low bar, the implicit standard is higher than it sounds; the licence is to take calculated risks, not to ship sub-interesting work. "We only need a result" is never the right framing for a scope decision.

**What this means for the outline.** Stage 2 is framed not as "can we reproduce the ancestor papers' averaging result" (we must, and will, as a conservative baseline) but as "what does per-frame shielding-tensor time-series analysis show, given tools the ancestors did not have." Stage 3 is framed not as "can we train a GNN that predicts shifts" (SHIFTX already works) but as "does our deep physics-grounded extraction provide signal an ML model can use to recover interpretable structure at a scale where big-data ML cannot." The averaging regime of de Gortari 2010 and Markwick 2010 is the baseline our method must not regress from; the exploratory arm is where the contribution lives.

---

## Section 1 — Introduction

The NMR chemical shift as a structural probe with unmatched sensitivity to local chemical environment. The gap between measurable per-atom shifts and the atomistic-conformational reality they sample. The thesis addresses that gap across three coordinated stages.

**Cross-paper consistency worth highlighting (thesis-motivating observation):** backbone ¹⁵N is the hardest element-channel in our Stage 1 ridge regression (R² = 0.387, pooled; all other element channels sit higher) AND it is the largest beneficiary of MD averaging in Markwick 2010 (2.89 → 1.84 ppm RMSD, the biggest gap across five nuclei). Same nucleus, same story from two independent angles: backbone ¹⁵N carries dynamics, which is why it is hard for static predictors and why ensemble averaging helps most there. This is the empirical thread that connects Stages 1 and 2.

---

## Section 2 — NMR physics as backdrop

Papers anchoring the section:
- **Ramsey 1950** *Phys. Rev.* 78, 699 (⧖, BL library) — paramagnetic shielding framework that made the chemical shift a calculable quantity.
- **Buckingham 1960** *Can. J. Chem.* 38, 300 (✓) — power-series decomposition of shielding in the local electric field; solvent and polar-group classification still in use.
- **Bryant 1983** *J. Chem. Educ.* 60, 933 (✓) — timescale framework: chemical-shift averaging operates at μs–ms (coalescence), relaxation at ps–ns; the two regimes report different physics and must not be conflated.
- **Facelli 2011** *Prog. NMR Spectrosc.* 58, 176 (✓) — modern tensor-theory anchor; the "fortuitous cancellation of errors in individual tensor components can lead to artificially high agreement in the isotropic chemical shifts" quote directly supports rank-2 analysis over isotropic trace.

Arc of the section: shielding is a tensor defined by the second derivative of the electronic energy with respect to nuclear magnetic moment and external field; it admits a power-series expansion in local field (Buckingham); the observable we measure is an ensemble- and time-average over conformational substates out to the μs–ms coalescence limit (Bryant); tensor-component analysis is diagnostically stronger than isotropic-trace comparison (Facelli).

---

## Section 3 — Stage 1 landscape: static kernels from static structures

**Ancestors on disk (summaries available):**
- **Case 1995** *J. Biomol. NMR* 6, 341 (✓) — DFT-calibrated ring-current intensity factors and electric-field coefficient via methane probes. Direct thesis ancestor. 421 shifts across ~10 ring types; Johnson-Bovey and Haigh-Mallion statistically indistinguishable for scalar shifts.
- **Sahakyan-Vendruscolo 2013** *J. Phys. Chem. B* 117, 1989 (✓) — RNA-base extension of the Case apparatus. Establishes heavy-nucleus dominance of electric-field effects over ring currents (¹³C, ¹⁵N, ¹⁷O).
- **Moyna 1998** *J. Chem. Inf. Comput. Sci.* 38, 702 (✓) — systematic comparison of HM, JB, and point-dipole models on 2992 protein ¹H shifts; all three cap at ~0.22 ppm RMS, the classical-kernel accuracy floor.
- **Pelloni 2004** *Org. Lett.* 6, 4451 (✓) — differential-Biot-Savart interpretation of ring-current signs; nearest-element dominance as a sanity-check intuition for kernel outputs.
- **Agarwal 1977** *Can. J. Chem.* 55, 2575 (✓) — [10]-paracyclophane as an out-of-plane-proton validation system; partitions observed shifts into ~65% ring current / ~35% local bond anisotropy; demonstrates the JB two-loop separation parameter is numerically unnecessary when local anisotropy is handled separately.

**Foundational papers pending acquisition:**
- Pople 1956 *J. Chem. Phys.* 24, 1111 (⧖)
- Johnson-Bovey 1958 *J. Chem. Phys.* 29, 1012 (⧖)
- Haigh-Mallion 1980 *Prog. NMR Spectrosc.* 13, 303 (⧖)
- Lazzeretti 2000 *Prog. NMR Spectrosc.* 36, 1 (⧖) — the pedagogical compilation of ring-current theory

**Other pending context:** Ösapay-Case 1991/1994 (⧖, the experimental protein-shift dataset used by downstream work), Williamson-Asakura 1991/1992/1995 (⧖, Cα calibration), SHIFTX/SHIFTS/SPARTA predictors (⧖).

**Where the thesis sits:**
- 720 proteins / 446K atoms / 55 simultaneous kernels — the largest calibration of this kind.
- Rank-2 shielding tensors with T0/T1/T2 decomposition, not scalar shifts.
- Per-element AND per-atom-type stratification. Key finding (2026-04-15): backbone N is the hardest atom-type (R² = 0.387), but sidechain N is the second-best (R² = 0.887). Element-pooled analysis had hidden this.
- Overall R² = 0.818 (pooled).

**Gaps the thesis addresses:**
- Tensor components (T2) carry information the isotropic trace averages away — Facelli's argument, Gershoni-Poranne's NICS_zz analog, our kernels' output discipline.
- Explicit-solvent handling — the gap Case 1995 closes his paper on.
- Atom-type-resolved analysis — largely absent from the predictor-ecosystem literature.

---

## Section 4 — Stage 2 landscape: per-frame shielding time series from MD

**Reframe from the morning version:** the object of Stage 2 is the per-frame shielding-tensor **time series** — for each protein, a (T, N, 9) array of shielding-tensor components across T frames, N atoms, and 9 Cartesian tensor components (decomposable into T0/T1/T2 irreducible representations). With the 55 kernels along one more axis, the analysis-ready object is (T, N, K, 9). Averaging over the T axis collapses this to the observable the NMR community has historically treated as the target; it is one valid diagnostic but not the end of the analysis. The ancestor literature works almost exclusively in the averaged regime, because until recently the data volume, tooling, and computational routes to per-frame QM made the full time series impractical. Two of those three constraints have loosened; the third (QM per frame at protein scale) is exactly the gap our library closes. This reframes de Gortari 2010 and Markwick 2010 from "precedents for what we're doing" to "conservative baselines we must not regress from, carrying only part of the available signal."

### 4.1 Averaging-regime ancestors (conservative baseline)

- **de Gortari 2010** *JACS* 132, 5993 (✓) — MLF tripeptide in solid state. DFT-grounded (GIPAW/PBE via CASTEP). 5 μs classical MD rotameric averaging improves ¹³C RMSD from 4.2 ppm (static) to 3.0 ppm. Demonstrates the MD-averaging principle with DFT-quality shift calculation, but on a 22-carbon peptide in a crystal.
- **Markwick 2010** *JACS Commun.* 132, 1220 (✓) — IκBα ankyrin repeat protein in solution. Empirical SHIFTX predictor (not DFT). AMD with dihedral boost; cumulative RMSD 8.19 → 5.92 ppm. Demonstrates the MD-averaging principle on real proteins with per-residue spatial cross-validation against R1/R2 ms-exchange regions (Cervantes 2009).

Complementary scope of the two ancestors:

| | de Gortari 2010 | Markwick 2010 |
|---|---|---|
| Shift predictor | DFT (GIPAW/PBE) | SHIFTX (empirical) |
| System | tripeptide (22 ¹³C) | protein (~140 residues) |
| Phase | solid state (P2₁2₁2₁) | solution |
| MD type | classical + CP-MD | AMD (enhanced sampling) |
| MD length | 5 μs classical | 20 × 10 M AMD steps |
| Solvent | crystal | explicit TIP3P |
| Result | 4.2 → 3.0 ppm | 8.19 → 5.92 ppm |
| Time axis | collapsed to average | collapsed to average |

Their results are our conservative-baseline territory: if our MD-averaged DFT-quality shielding tensors do not recover the 0.5–1 ppm gain these papers report over static structures, something is wrong upstream.

### 4.2 Kinds of means — because the kind matters

In the averaging regime, the literature uses several distinct weightings, each appropriate for different sampling conditions:

- **Simple arithmetic mean** over Boltzmann-sampled frames — what plain MD delivers if sampling is converged. Our 25 ns AMBER ff14SB with 600 frames is in this regime if convergence checks pass.
- **Free-energy-weighted mean** (Markwick) — required when enhanced sampling distorts the Boltzmann distribution; would apply if we move to AMD or metadynamics.
- **Rotameric / conformational-substate mean** (de Gortari) — identify substates by dihedral-angle clustering, average within each, weight by substate population. Two-level decomposition.
- **Lipari-Szabo separable mean** — fast ps-ns internal motion characterised by an order parameter S², slow motion factored out. Relaxation-centric historical formalism.
- **Block-averaged mean with convergence diagnostics** (Grossfield-Zuckerman 2009 H7 [CORE], ⧖) — standard MD discipline for error estimation. Particularly relevant at 2 GB per-protein time-series volumes where whether the data is 2 GB of independent signal or 2 GB of correlated noise is a real question.

### 4.3 Time-series primitives from the NMR-MD tradition

One level up from the mean:

- **Lipari & Szabo 1982** *JACS* 104, 4546 (⧖, M1 [CORE]) — autocorrelation function C(t) of relaxation-active motions, Fourier-transformed to spectral density J(ω) for relaxation rate prediction. The canonical NMR-MD time-series object.
- **Amadei-Linssen-Berendsen 1993** (⧖, M19 [CORE]) — Essential Dynamics: PCA on atomic coordinate trajectory to extract collective motion modes. Standard template for decomposing a T × N trajectory into mode structure; generalises to any time-indexed signal, not just coordinates.
- **Tugarinov & Kay 2003** (⧖, M10 [CORE]) — cross-correlated relaxation; using multiple NMR observables jointly to constrain dynamics.
- **Kasinath & Wand 2013** (⧖, M7) — order parameter as an entropy meter; biological interpretation of S² useful for thesis introduction framing.
- **Showalter & Brüschweiler 2007** (⧖, H10) — validation of MD force fields against NMR relaxation data; checks whether the dynamics our ff14SB is producing matches what NMR sees.
- **Berne & Harp 1970** (⧖, M28 [CORE]) — foundational time-correlation-function treatment in MD; the statistical-mechanics grounding.

### 4.4 Methods mature elsewhere, underused for shielding time series

These are the directions licensed by the "small + tooled-up" positioning. Each is mature in some adjacent field but (as far as we can tell) has not been applied to per-frame protein shielding-tensor series at protein scale.

- **Tensor decomposition (Tucker, CP, hierarchical)** on the (T, N, K, 9) object — the NYC-taxi-tensor analogy already in your Stage 2 stats memory. Tensorly, rTensor, or PyTorch-based implementations are stable. Modes of the decomposition should be interpretable along kernel, residue-type, or secondary-structure lines if the physics is coherent.
- **Markov State Models** (PyEMMA, deeptime) — discretise the trajectory into kinetic states and learn a transition matrix; state-resolved shielding becomes the observable. MSM literature mature; per-frame computed-property integration is less developed.
- **Dynamic factor models / state-space models** — Kalman-style latent-variable time-series modelling. Common in econometrics, climate, fMRI; rare in MD for non-coordinate observables.
- **Wavelet / time-frequency decomposition** — localise shielding fluctuations in time and frequency; relevant if there is reason to suspect non-stationarity across the 25 ns.
- **Memory-kernel / Mori-Zwanzig projection** — non-Markovian decomposition of the shielding time series into systematic, memory, and noise components. Rigorous; less work has been done applying it to computed-property series than to thermodynamic variables.
- **Generalised-Langevin / fluctuation-dissipation treatments** — connect the shielding fluctuations to an effective temperature or response function.

The common theme is that the data object (T, N, K, 9) is structurally the same kind of thing these fields already analyse. The content — per-frame shielding tensor — is novel; the methodology is not.

### 4.5 Where the thesis sits

- Data object: (T, N, K, 9) per protein. T ≈ 600 frames, N = hundreds of atoms, K = 55 kernels, 9 = shielding-tensor Cartesian components (or equivalent T0/T1/T2 decomposition).
- Solution-phase proteins with cluster DFT at r²SCAN in ORCA on clipped environments (the Stage-2 DFT reference).
- Plain 25 ns AMBER ff14SB MD with explicit TIP3P + APBS — no enhanced sampling, so the Boltzmann approximation is intact and reweighting is not required.
- 10 calibration proteins × 260 DFTs complete (2026-04-18); extending to the full 685-protein fleet as MDs land.
- Conservative-baseline analysis: arithmetic-mean and substate-mean shielding tensors vs. the DFT reference, cross-validated spatially against BMRB R1/R2 where available.
- Exploratory analysis: tensor decomposition of (T, N, K, 9), MSM-style kinetic decomposition, autocorrelation / spectral density per-atom-per-component, PCA/ED on shielding trajectories by analogy to coordinate ED.

### 4.6 Validation strategy

**Conservative (baseline):**
- Per-residue averaging-improvement regions should spatially coincide with R1/R2 ms-exchange regions (Markwick template, adapted from Cervantes 2009 on IκBα).
- Is a diagnostic, not a training target.
- Cross-stage corroboration of the ¹⁵N-hard finding from Stage 1 (backbone N hardest, sidechain N second-best).
- Audit of the 10 calibration proteins for dynamical diversity — Markwick's ubiquitin-vs-IκBα contrast shows improvement scales with heterogeneity; a uniformly rigid calibration set would mute Stage 2 signal.

**Exploratory (extension):**
- Tensor-decomposition modes must have chemical interpretation — kernel-aligned, residue-type-aligned, or secondary-structure-aligned. If the top modes of a Tucker decomposition of (T, N, K, 9) don't admit physical reading, the factorisation is capturing statistical structure without physics and must be reported that way.
- Autocorrelation decay times per atom should cluster by residue type or secondary-structure context in interpretable ways, and (where possible) correlate with NMR-derived motional timescales in the ns-ps window.
- MSM state assignments should correlate with structural features (rotamer states, hydrogen-bond patterns, secondary-structure fluctuations) rather than appearing as statistical artefacts.

A null result on the exploratory analysis is still a thesis chapter — reporting that the MD-time-series signal does not carry structure beyond the mean, carefully framed, is itself a methodological finding at our data volume.

---

## Section 5 — Existing ML systems for NMR shift prediction and structural biology, and their basis

This section is a prior-art survey of the ML systems a thesis in this area must position against. It is deliberately separate from the stage-by-stage landscape because these systems cut across stages: classical shift predictors are Stage 1 competitors, structure-prediction embedding models are Stage 3 infrastructure, and full-tensor equivariant systems are the closest architectural precedent for what we might build.

**The uncomfortable truth of the landscape (from Positioning, repeated here for the lit-review reader):** SHIFTX and its descendants solve protein chemical-shift prediction well enough for working structural biology. Papers keep appearing because the field is boxed, tidy, and adjacent to AlphaFold's high-status ML-for-structural-biology domain. But big-data ML cannot scale the way AlphaFold did, because the BMRB even in aggregate is too small for AlphaFold-scale training. The field therefore selects for novel angles.

### 5.1 Empirical shift predictors (the solved-in-practice tradition)

Basis: regression on a database of measured protein shifts, using sequence and/or static-structure features. Output is usually isotropic shift per nucleus.

- **SHIFTX2** (Han-Wishart 2011) and predecessors SHIFTX, SHIFTY, RefDB (Wishart lineage). Hybrid empirical + machine-learning regression on BMRB-derived training sets. The field's de facto working tool. (⧖)
- **SPARTA+** (Shen-Bax 2010) — neural network with structural features, trained on refined structure–shift pairs. Bax-group lineage. (⧖)
- **CamShift** (Kohlhoff-Vendruscolo 2009) — empirical predictor inspired by classical ring-current and electric-field terms, integrated into ensemble-based structure determination. Vendruscolo lineage; part of what became CHESHIRE. (⧖)
- **PROSHIFT** (Meiler 2003) — neural network; historical but cited. (⧖)
- **UCBShift** (Han lab, Berkeley, 2020s) — deep-learning predictor using multiple structural features. (⧖)
- **CSpred** — deep-learning chemical-shift predictor; part of the 2020s generation. (⧖)

What these give us: a floor. If our MD-averaged Stage 1 + Stage 2 residuals against experiment are much worse than SHIFTX2's, something upstream is broken. If they are comparable, the ML systems become the like-for-like comparator for our classical + DFT-calibrated approach.

### 5.2 DFT-tensor-based approaches

Basis: pre-computed DFT shielding tensors on structural fragments or a fragment database, with interpolation or fragment-matching as the prediction mechanism.

- **PROCS15** (Kjaergaard and co-workers) — database of DFT-calculated shielding tensors on protein backbone and sidechain fragments. Jessica pulled the full PROCS15 tensor set into the project database during the reconnaissance phase. This was scaffolding for the original "ML on poses + PROCS15" direction, which the pivot to r²SCAN + per-frame extraction superseded. PROCS15 may re-enter as an independent DFT-tensor cross-check for Stage 1 and Stage 2 results. (⧖ at paper level; database is locally accessible)
- **CHESHIRE** (Cavalli-Salvatella-Dobson-Vendruscolo 2007) — chemical-shift-driven protein structure determination. Combines empirical predictors with ensemble MD. Vendruscolo-group lineage. (⧖)
- **ShiftML** (Paruzzo-Hofstetter-Emsley and co-workers) — kernel-based ML for solid-state NMR shielding tensors. Different regime (solid state) but methodologically instructive because they do output tensors; worth reading for how they validate tensor predictions. (⧖)

### 5.3 Structure-prediction embedding models (Stage 3 infrastructure)

Basis: transformer or attention architectures trained on massive protein sequence (ESM-2) or sequence + structure (AlphaFold, OpenFold, RoseTTAFold) data. Outputs that matter for us are the learned embeddings, which can be used as features for downstream tasks like shift prediction.

- **AlphaFold2 / AlphaFold3** (Jumper 2021 / Abramson 2024) — the field-defining demonstration that attention-based models trained on massive data predict protein structure. Infrastructure for structure-dependent downstream tasks; establishes the high-status ML narrative around which our field orbits. (⧖)
- **OpenFold / OpenFold3** (Ahdritz et al.) — open reimplementations of AlphaFold with trainable weights and modifiable architectures. **The OpenFold3 embeddings are our current Stage 3 ML backing**, being redone on DGX Spark hardware to submittable standard using the ESM-2-replication pipeline Jessica built over the break. (⧖)
- **ESM-2** (Lin-Akin-Verkuil et al. 2023) — protein language model, sequence-only. Jessica replicated it; advisor rejected it for our work on two grounds: airless sequence-only modelling without geometry, and failure in the original paper to use mmseqs2 to partition test and training data (risk of data leakage). In our landscape as a contrast cite, not as a backing model. (⧖)
- **RoseTTAFold** (Baek et al. 2021) — contemporary to AlphaFold2, similar embedding character. Contextual only. (⧖)

### 5.4 Equivariant ML for NMR observables

Basis: neural networks with built-in SO(3) or E(3) equivariance, so that tensor-valued outputs transform consistently with input geometry. Closest architectural precedent to what a Stage 3 full-tensor model would look like.

- **Venetos et al. 2023** — equivariant GNN predicting full NMR shielding tensors on silicon oxides. The direct architectural precedent; differs by system class (inorganic oxides, not proteins). (⧖)
- **Yang-Chakraborty-White 2021** — earlier GNN for scalar protein shifts; scalar-only so not a tensor precedent but relevant as a protein-GNN reference. (⧖)
- **Kondor 2025** (✓) — the mathematical framework surveyed in Section 6.
- **Kondor-Trivedi 2018** (✓, supplementary) — the generalisation theorem underpinning the framework.
- Equivariant-NN family (TFN, e3nn, NequIP, MACE) — see Section 6 for detail. Most relevant here as "the toolbox Venetos 2023 drew from, which a Stage 3 model would also draw from."

### 5.5 Recent lit reviews and surveys on ML and protein structure

Needed for the outline to position our work in the current (2023–2025) ML-for-structural-biology conversation. We do not have these on disk yet; acquisition priority raised because the Stage 3 and Positioning arguments need to cite recent review-level framing.

- **Jumper et al. 2021** *Nature* — AlphaFold2 paper. (⧖)
- **Abramson et al. 2024** *Nature* — AlphaFold3. (⧖)
- **Bronstein et al. 2021** *Geometric Deep Learning* — the "5G" taxonomy; foundational review. (⧖)
- **Gerken et al. 2023** *Artif. Intell. Rev.* — geometric deep learning and equivariant neural networks survey (referenced by Kondor 2025's E16). (⧖)
- **A current review of ML for chemical shifts** — needs search; the Venetos 2023 paper likely references the relevant one. (⧖)
- **A current review of ML for protein structure** — the AlphaFold family papers + any Nat. Rev. Mol. Cell Biol. or similar review. (⧖)

### 5.6 Where the thesis sits in this existing-systems landscape

- Against empirical predictors (SHIFTX2 et al.): we do not compete on the isotropic-shift accuracy metric at which they are well-tuned. We compete on tensor-level structure and on per-frame time-series signal that their static-input architectures cannot take as input.
- Against DFT-tensor-database methods (PROCS15, CHESHIRE): we replace the database with a per-frame DFT calibration computed fresh on our MD-extracted geometries, and we retain the full time series rather than averaging before prediction.
- Against structure-prediction embeddings (AlphaFold / OpenFold3): we *use* these, not compete with them. OpenFold3 embeddings are our Stage 3 ML input; the prediction target is our DFT-calibrated shielding tensors.
- Against equivariant GNNs for NMR (Venetos 2023): closest architectural competitor. Their system is inorganic; ours is protein. Their output is the NMR tensor directly; ours may be either the tensor or the residual from our kernel decomposition. The question Stage 3 addresses is whether starting from kernel residuals rather than raw geometry lets a smaller model recover the same signal, interpretably.

---

## Section 6 — Stage 3 landscape: model evaluation against our data

**Ancestors on disk:**
- **Kondor 2025** *PNAS* 122, e2415656122 (✓) — principles of equivariant neural networks for physics and chemistry. Peter-Weyl decomposition into irreps, Clebsch-Gordan product as the canonical polynomial equivariant nonlinearity. Design-space reference for any SO(3)-equivariant tensor-output model.
- **Kondor-Trivedi 2018** *ICML* (✓, supplementary) — the generalisation theorem underpinning the 2025 review.

**Existing ML infrastructure we build on (detail in Section 5):**
- OpenFold3 embeddings, applied via the ESM-2-replication pipeline on DGX Spark hardware. Produces per-residue feature vectors that feed a regression / GNN head targeting our DFT-calibrated shielding tensors and/or the kernel-residual signal from Stage 1.
- Traditional GNN for the same problem as a comparator — standard GNN architecture, same input features, same target. Reports on both regimes let the thesis say which design gets more out of the data.

**Pending — the equivariant-NN toolbox (Section 5.4 detail, also used here):**
- Thomas-Smidt 2018 (⧖, TFN)
- Geiger-Smidt 2022 (⧖, e3nn)
- Batzner 2022 (⧖, NequIP)
- Batatia 2022 (⧖, MACE)
- Cohen-Welling 2016 (⧖, group-equivariant CNN origin)

**Where the thesis sits:**
- Two models, reported together: an OpenFold3-embedding-based regression / GNN head targeting the shielding tensor (or kernel residuals), and a traditional GNN as a like-for-like comparator. Both take the Stage 2 time-series-aware data and output predictions; both are evaluated on the same held-out protein set.
- Whether residuals from the Stage 1 kernel ridge are noise-like (ridge suffices end of story) or still carry angular signal (Stage 3 learning is justified) is an open empirical question that Stage 2 and 3 together answer.
- If the learning side succeeds, the contribution is a full-tensor protein shielding prediction at our scale from a novel input (kernel-decomposed time-series + embeddings) — unexplored territory.

---

## Section 7 — The thesis contribution

The three-fold claim (the "schtick"):

1. **Deep, fast, validated classical extraction.** Per-atom shielding information extracted from scalar through T2, on frames, via a C++ kernel library producing hundreds of thousands of evaluations cheaply. The extraction is validated on two fronts: signal-level against DFT shielding tensors (Stage 1 calibration, R² = 0.818 and ongoing Stage 2 per-frame DFTs), and mathematically against known dimensional structure (T0/T1/T2 irreducible representations, verified by construction). This is the conservative technical contribution and is already in credible first-draft state.

2. **Physics-grounded interpretation of the extractions.** Each kernel is a concrete physical mechanism (ring currents, electric-field effects, bond anisotropy, hydration, etc.) whose contribution to shielding is known from 60 years of NMR theory. Stage 1's per-element and per-atom-type stratification shows which mechanisms dominate where. The interpretability that follows — backbone N hard because of dynamics, sidechain N second-best because of hydration and EFG, etc. — is the physics side of the contribution and distinguishes this work from black-box ML predictors whose weights do not correspond to physical mechanisms.

3. **Pushing the time-series data at ML in two modes.** The (T, N, K, 9) Stage 2 time-series data is used as input to (a) an OpenFold3-embedding-based learning pipeline (Jessica's ESM-2 replication style, ported to OpenFold3, running on DGX Spark) and (b) a traditional GNN for the same problem, reported as a comparator. Both models target the shielding tensor and/or the kernel-ridge residual from Stage 1. Reporting both modes lets the thesis say what kind of ML architecture gets more out of DFT-grounded physics-kernel-decomposed time-series data than what kind gets from raw geometry.

Additional supporting claims:
- Rank-2 shielding tensor output preserved end-to-end through kernels, calibration, MD trajectory, and any downstream ML model.
- Explicit-solvent DFT-grounded shielding at protein scale, in solution — the cell Case 1995 flagged as open, de Gortari 2010 did not reach (solid state), and Markwick 2010 did not reach (empirical predictor).
- Calibration against ORCA r²SCAN cluster DFT on clipped environments, scaling to proteins where GIPAW would not.
- Evaluation against classical ring-current / electric-field kernels and (Section 5's landscape of) existing ML systems on common protein data.

We are past the fallback point. The three-fold claim is the deliverable; the conservative baseline work sits underneath it as evidence, not as a back-up plan.

**Claim-shape convention.** Throughout the thesis, the kernel-versus-DFT comparison is phrased as *correlation* (r² good enough to reflect that our kernel extractors saw a signal), not *match* (pointwise numerical agreement). The validation currency is correlation on several dimensions — mean shielding, S²_σ amplitude, principal-mode structure, tensor-factor decomposition — not prediction accuracy against experiment. "The kernels saw the signal" is the honest version of the Stage 2 claim.

---

## Growing edges / TODO

Things to thicken as material lands:

**Foundational and ancestor-lineage:**
- Acquire the foundational ring-current set (Pople, JB, HM 1980 and earlier) — BL library trip.
- Acquire Ramsey 1950 for Section 2 depth.
- Acquire Ösapay-Case 1991/1994 as the experimental-dataset anchor cited across the empirical-predictor lineage.
- Acquire Vila-Arnautova-Scheraga fragment-DFT protein shift papers (Facelli 2011 refs 201–203) as a direct methodological parallel to our cluster-DFT approach.

**Stage 2 time-series methodology (priority raised after 4.3 / 4.4 reframe):**
- Acquire Lipari-Szabo 1982 (M1 [CORE]) for the canonical NMR time-series formalism.
- Acquire Grossfield-Zuckerman 2009 (H7 [CORE]) for MD convergence and block-averaging discipline — directly relevant to whether 2 GB of time-series data is 2 GB of signal.
- Acquire Amadei-Linssen-Berendsen 1993 (M19 [CORE]) as the Essential Dynamics / PCA-on-trajectory template we intend to generalise to shielding.
- Acquire Tugarinov-Kay 2003 (M10 [CORE]) for cross-correlated relaxation framing.
- Acquire Kasinath-Wand 2013 (M7) for order-parameter-as-entropy-meter framing.
- Acquire Berne-Harp 1970 (M28 [CORE]) for statistical-mechanical grounding of time-correlation functions in MD.
- Acquire Showalter-Brüschweiler 2007 (H10) for validating MD force fields against NMR dynamics data.
- Consider for the exploratory analysis tools: a canonical Markov State Model paper (Prinz et al. or Chodera's reviews), a tensor-decomposition reference (Kolda-Bader 2009 SIAM Review), and a state-space-model textbook chapter if we commit to that direction.

**Stage 3 learning side:**
- Acquire the equivariant-NN family (arXiv-accessible without BL trip, mostly — Thomas-Smidt 2018, Geiger-Smidt 2022, Batzner 2022, Batatia 2022, Cohen-Welling 2016).
- Acquire Venetos 2023 as the closest prior architecture for full-tensor NMR prediction.

**Structural decisions:**
- Decide whether Section 2 should absorb a subsection on magnetic-aromaticity criteria (Gershoni-Poranne 2015) or keep that as contextual depth in Section 3.
- Decide which of Section 4.4's exploratory methods get space in the thesis body versus appendix — likely tensor decomposition in body, MSM / memory-kernel / state-space in appendix unless one produces a striking result.
- Fold any new Stage-1 or Stage-2 findings that emerge from ongoing analysis.

---

## Relationship to Monday's program-rubric submission

Monday's submission is a separate, program-bounded document:
- An 8-item annotated bibliography (already drafted; the 8 are locked: Bryant, Facelli, Buckingham, Case, Gershoni-Poranne, de Gortari, Sahakyan-Vendruscolo, Kondor 2025).
- A lit-review outline bounded by the program's rubric (specifics not yet in hand).

All 8 Monday-submission papers appear in this proposed outline. The framing and supporting context in this document will inform Monday's writing but won't fit it. Monday's submission is a one-off deliverable; this document outlives it and supports the thesis-chapter writing that follows.
