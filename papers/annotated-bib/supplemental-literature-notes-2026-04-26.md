# Supplemental Literature Notes

**Author:** Jessica Hansberry
**Date:** 2026-04-26

Working notes from a wide search around protein NMR chemical-shift prediction, ML predictors, DFT/shielding validation, and MD/NMR dynamics. This is a supplement to the Monday annotated bibliography, not a final bibliography and not an evidence-card file.

Coverage discipline: the literature-review job is to show competent command of the existing ground before making any claim about what this project adds. Candidate novelty language should stay quiet: "we did this specific thing; this may not have been done in quite this form" is the outer edge unless the evidence base says more.

Project-specific Stage 2 framing: Stage 2 is not a separate attempt to prove a full relaxation/NMR-dynamics theory. It is close in spirit to Stage 1: extract tensor-valued classical/geometric quantities from structures or frames, compare their dimensions and relationships against rSCAN shielding calculations across atom classes, and report what relationships appear. Correlations with experimental relaxation/order-parameter data would be valuable if they appear, but they are a bonus discussion point, not the central result being promised.

Status labels:

- **On disk:** PDF already present in `references/` and should have or receive the normal `references-text/`, `references-images/`, and `references-meta/` companions.
- **Open, not ingested:** full text appears openly available, but the PDF is not yet in the local corpus.
- **Pending acquisition:** likely needs Birkbeck/British Library/institutional access or a manual fetch.
- **Metadata only:** found by bibliographic record/abstract; do not cite for claim-bearing prose until acquired and processed.

Evidence-card discipline for later thesis argument: nothing in this file is citable in that mode unless the full text is on disk and a card has been made from it.

---

## Modern Protein Chemical-Shift Predictors

### Ptaszek, A.L., Li, J., Konrat, R., Platzer, G. & Head-Gordon, T. (2024) "UCBShift 2.0: Bridging the Gap from Backbone to Side Chain Protein Chemical Shift Prediction for Protein Structures." *Journal of the American Chemical Society* 146(46), 31733-31745. DOI: 10.1021/jacs.4c10474.

**Status:** Open, not ingested. Sources found: ACS landing page, PubMed, PMC author manuscript.

Extends UCBShift from backbone atoms to side-chain chemical shifts. The method keeps the two-module UCBShift architecture: a transfer module using sequence and structural alignments against experimental reference shifts, plus a machine-learning module using curated structural and chemistry features. This is a strong current comparator because it is not presented as sequence-only or end-to-end black-box learning; it is a modern feature-rich predictor with explicit curation and transfer logic. It also gives a recent benchmark against SHIFTX2 for side-chain atoms.

Quiet take: useful for describing the current practical state of protein shift prediction and for making sure our work is not framed as "we built another shift predictor." The relevant contrast is narrower: their target is accurate whole-protein shift prediction from structures; our interest is whether explicit tensor/geometric feature extraction from structures and trajectories exposes signal that can feed learning and be checked against quantum calculations.

Source URLs:
- https://pubs.acs.org/doi/abs/10.1021/jacs.4c10474
- https://pubmed.ncbi.nlm.nih.gov/39531038/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11784523/

### Li, J., Bennett, K.C., Liu, Y., Martin, M.V. & Head-Gordon, T. (2020) "Accurate prediction of chemical shifts for aqueous protein structure on 'Real World' data." *Chemical Science* 11, 3180-3191. DOI: 10.1039/C9SC06561J.

**Status:** In `references/incoming/` as `li-2020-ucbshift-real-world-data.pdf`; not yet ingested.

Introduces UCBShift, combining sequence/structure transfer with random-forest prediction from extracted features. The paper is especially useful because it foregrounds "real-world" benchmark behavior, chemical-shift data curation, and the effect of outlier filtering. It explicitly discusses RefDB, SPARTA+, SHIFTX2, feature extraction, sequence homology transfer, structural homology transfer, and the ambiguity introduced by arbitrary post-prediction filtering.

Quiet take: this belongs high in the supplemental bibliography because it is adjacent to the data and benchmarking problem rather than merely another model paper. It is a good source for discussing why BMRB/RefDB-based learning requires care, and why a predictor can look different under heavily curated versus less filtered test conditions.

Source URL:
- https://pubs.rsc.org/en/content/articlehtml/2020/sc/c9sc06561j

### Zhu, H., Hu, L., Yang, Y. & Chen, Z. (2025) "A novel approach to protein chemical shift prediction from sequences using a protein language model." *Digital Discovery* 4, 331-337. DOI: 10.1039/D4DD00367E.

**Status:** In `references/incoming/` as `zhu-2025-plm-cs-sequence-chemical-shifts.pdf`; not yet ingested.

Introduces PLM-CS, a sequence-only predictor using frozen ESM2-650M embeddings followed by a transformer predictor. The training set is RefDB-derived, with one model per backbone atom type. The paper reports that PLM-CS is less accurate than SHIFTX2 on the SHIFTX test set, but more stable on a custom solution-NMR structure test set. The ablation study reports that ESM embeddings outperform a one-hot amino-acid baseline for most backbone atoms, while Cbeta remains strongly amino-acid-type dependent.

Quiet take: useful as a modern contrast case. It shows that sequence embeddings contain shift-relevant information, but the authors' own comparison still leaves room for explicit structural and expert-selected features. The solution-NMR versus crystal/AlphaFold comparison may be useful when discussing structure source and ensemble mismatch.

Source URL:
- https://pubs.rsc.org/en/content/articlehtml/2025/dd/d4dd00367e

### Yang, Z., Chakraborty, M. & White, A.D. (2021) "Predicting chemical shifts with graph neural networks." *Chemical Science* 12, 10802-10809. DOI: 10.1039/D1SC01895G.

**Status:** In `references/incoming/` as `yang-white-2021-gnn-chemical-shifts.pdf`; not yet ingested.

Message-passing GNN for chemical shifts, trained with RefDB, SHIFTX, and a small organic-molecule dataset. The input graph uses atom identities plus covalent and non-bonded edges with distances. The paper emphasizes differentiability with respect to structure, speed, non-bonded effects, hydrogen-bond downfield shifts, and the possibility of use inside simulation or structure-inference workflows. It contrasts accurate but feature/homology-dependent methods such as SHIFTX2 with differentiable but narrower methods such as CamShift.

Quiet take: important near-neighbor for the eventual GNN stage. It does not do the same thing as the extractor, but it sets a precedent for graph-based shift prediction that can use distances and non-bonded contacts without hand-engineered protein features. It also gives language for why differentiability and simulation compatibility matter.

Source URL:
- https://pubs.rsc.org/en/content/articlehtml/2021/sc/d1sc01895g

### Han, B., Liu, Y., Ginzinger, S.W. & Wishart, D.S. (2011) "SHIFTX2: significantly improved protein chemical shift prediction." *Journal of Biomolecular NMR* 50, 43-57. DOI: 10.1007/s10858-011-9478-4.

**Status:** In `references/incoming/` as `han-2011-shiftx2-protein-chemical-shift-prediction.pdf`; not yet ingested.

Hybrid predictor combining sequence-based transfer with structure-based prediction, using a large training set and engineered features including torsions, solvent accessibility, hydrogen bonding, pH, and temperature. Still a practical baseline and a named comparator for UCBShift, PLM-CS, and many later papers.

Quiet take: this is a necessary acquisition/ingestion item because many newer papers define themselves against it. The relevant use is not to spend pages retelling SHIFTX2, but to know exactly what feature families and evaluation choices the field already treats as standard.

Source URL:
- https://link.springer.com/article/10.1007/s10858-011-9478-4

### Shen, Y. & Bax, A. (2010) "SPARTA+: a modest improvement in empirical NMR chemical shift prediction by means of an artificial neural network." *Journal of Biomolecular NMR* 48, 13-22. DOI: 10.1007/s10858-010-9433-9.

**Status:** Open HTML available via PMC; PDF download blocked from this environment. Manual fetch recommended.

Neural-network predictor trained on a carefully pruned database of high-resolution X-ray structures and backbone/Cbeta chemical shifts. Includes structural descriptors for backbone/side-chain conformation, hydrogen bonding, electric fields, and ring-current effects.

Quiet take: useful because it is an older ML method that already included classical physical descriptors. It helps avoid overstating novelty around "physics features plus ML"; the narrower question is what our tensor-preserving extractor and MD/DFT validation add.

Source URLs:
- https://pubmed.ncbi.nlm.nih.gov/20628786/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC2935510/

### Kohlhoff, K.J., Robustelli, P., Cavalli, A., Salvatella, X. & Vendruscolo, M. (2009) "Fast and accurate predictions of protein NMR chemical shifts from interatomic distances." *Journal of the American Chemical Society* 131, 13894-13895. DOI: 10.1021/ja903772t.

**Status:** Pending acquisition or open SI only; ACS article may need institutional access.

CamShift expresses chemical shifts through polynomial functions of interatomic distances. Because the functions are fast and differentiable, the method became useful in refinement and simulation contexts.

Quiet take: necessary background for the "distance-only differentiable predictor" line. It is also a useful foil for explicit geometric/tensor features: CamShift is intentionally compact and differentiable, while our extraction is broader and more diagnostic.

Source URL:
- https://pubs.acs.org/doi/10.1021/ja903772t

### Larsen, A.S., Bratholm, L.A., Christensen, A.S., Channir, M. & Jensen, J.H. (2015) "ProCS15: a DFT-based chemical shift predictor for backbone and Cbeta atoms in proteins." *PeerJ* 3, e1344. DOI: 10.7717/peerj.1344.

**Status:** In `references/incoming/` as `larsen-2015-procs15-dft-chemical-shift-predictor.pdf`; not yet ingested.

DFT-based protein chemical-shift predictor for backbone and Cbeta atoms. It is based on roughly 2.35 million OPBE/6-31G(d,p)//PM6 calculations on tripeptides and hydrogen-bonding model systems, and reports that the predictor computes isotropic shielding values for a protein structure in less than a second. The paper also states that structures and DFT calculations, including full NMR shielding tensors, were made available.

Quiet take: a must-have for this project because it is closer to the DFT-feature-extraction part than pure BMRB-learning papers. The reported availability of full shielding tensors is particularly relevant for tensor discussion and comparison practices.

Source URL:
- https://peerj.com/articles/1344/

### Bratholm, L.A., Christensen, A.S., Hamelryck, T., Jensen, J.H. & co-workers (2017) "Protein structure refinement using a quantum mechanics-based chemical shielding predictor." *Chemical Science* 8, 2061-2072. DOI: 10.1039/C6SC04344E.

**Status:** In `references/incoming/` as `bratholm-2017-procs15-qm-chemical-shielding-refinement.pdf`; not yet ingested.

Uses the ProCS15 quantum-mechanics-based shielding predictor in protein structure refinement. This is adjacent to ProCS15 itself but more directly concerned with what chemical shielding information can do to structures when used operationally.

Quiet take: useful bridge paper between DFT-derived shift prediction and structure/trajectory work. It is likely worth processing with ProCS15 rather than later as a separate literature island.

Source URL:
- https://pubs.rsc.org/en/content/articlehtml/2017/sc/c6sc04344e

---

## Databases, Re-Referencing, and Validation

### Zhang, H., Neal, S. & Wishart, D.S. (2003) "RefDB: A database of uniformly referenced protein chemical shifts." *Journal of Biomolecular NMR* 25, 173-195. DOI: 10.1023/A:1022836027055.

**Status:** Open or accessible, not ingested.

RefDB is a secondary database of reference-corrected protein chemical shifts derived from BMRB. The correction workflow used SHIFTX to predict expected shifts from available structures, then compared those predictions to BMRB assignments to identify referencing offsets and possible errors.

Quiet take: foundational data-source paper. It should be acquired before writing any detailed claim about BMRB/RefDB training data, re-referencing, or the limits of using RefDB-derived labels.

Source URLs:
- https://refdb.wishartlab.com/
- https://bmrb.io/published/refDB/
- https://link.springer.com/article/10.1023/A:1022836027055

### Rieping, W. & Vranken, W.F. (2010) "Validation of archived chemical shifts through atomic coordinates." *Proteins* 78(11), 2482-2489. DOI: 10.1002/prot.22756.

**Status:** In `references/incoming/` as `paruzzo-2018-shiftml-molecular-solids.pdf`; not yet ingested.

Presents VASCO, a method for validating and correcting archived chemical shifts using atomic coordinates, solvent-accessible surface area, and secondary structure. The paper is useful because it gives a coordinate-aware alternative to purely distributional checks and explicitly discusses per-atom Z scores and error estimates for correction values.

Quiet take: valuable for data-quality discussion. It also connects to one of the recurring issues in this project: chemical-shift disagreement can reflect assignment/reference problems, structural mismatch, dynamics, or model error; a useful workflow keeps those possibilities separate.

Source URL:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC2970900/

### Dashti, H., Tonelli, M., Lee, W., Westler, W.M., Cornilescu, G., Ulrich, E.L. & Markley, J.L. (2016) "Probabilistic validation of protein NMR chemical shift assignments." *Journal of Biomolecular NMR* 64, 17-25. DOI: 10.1007/s10858-015-0007-8.

**Status:** Open HTML available via PMC/ACS; direct PDF download blocked from this environment. Manual fetch recommended.

Introduces ARECA, a probabilistic method for validating chemical-shift assignments using NOE data. It is complementary to methods based only on chemical-shift statistics or 3D structures, and explicitly warns that local physical effects can make valid assignments look like statistical outliers.

Quiet take: useful guardrail for any future discussion of "outliers." The source supports a cautious distinction between assignment error, referencing error, unusual local environment, and predictor failure.

Source URLs:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4744101/
- https://pubmed.ncbi.nlm.nih.gov/26724815/

### Ginzinger, S.W., Gerick, F., Coles, M. & Heun, V. (2007) "CheckShift: automatic correction of inconsistent chemical shift referencing." *Journal of Biomolecular NMR* 39, 223-227. DOI: 10.1007/s10858-007-9191-5.

**Status:** Metadata only.

Reference-correction method that does not require an available structure. It is likely useful as a short methods/context citation around preprocessing chemical-shift data.

Quiet take: lower priority than RefDB, VASCO, and ARECA, but worth acquiring if data-cleaning discussion grows.

Source URL:
- https://pubmed.ncbi.nlm.nih.gov/17899394/

---

## DFT, Shielding Tensors, and ML Against Quantum Labels

### Paruzzo, F.M., Hofstetter, A., Musil, F., De, S., Ceriotti, M. & Emsley, L. (2018) "Chemical shifts in molecular solids by machine learning." *Nature Communications* 9, 4501. DOI: 10.1038/s41467-018-06972-x.

**Status:** Open HTML available via PMC; PDF download blocked from this environment. Manual fetch recommended.

ShiftML uses local-environment descriptors and Gaussian process regression to predict GIPAW DFT chemical shifts for molecular solids. Because experimental chemical-shift databases for molecular solids are too small, the model is trained on DFT-calculated shifts. The paper explicitly notes that using computed shifts avoids ambiguity from dynamics/distributions in the experimental label.

Quiet take: important method analogue for "train or validate against computed quantum labels." Not a protein paper, but methodologically relevant to using rSCAN/DFT as a physics check rather than treating experimental shifts as the only truth source.

Source URL:
- https://www.nature.com/articles/s41467-018-06972-x

### Han, C., Zhang, D., Xia, S. & Zhang, Y. (2024) "Accurate Prediction of NMR Chemical Shifts: Integrating DFT Calculations with Three-Dimensional Graph Neural Networks." *Journal of Chemical Theory and Computation* 20(12), 5250-5258. DOI: 10.1021/acs.jctc.4c00422.

**Status:** Open HTML available via PMC; PDF download blocked from this environment. Manual fetch recommended.

Introduces CSTShift for small molecules. The model combines 3D graph neural networks with DFT-calculated shielding tensor descriptors. It is not a protein predictor, but it is directly relevant to the idea that shielding tensor descriptors can be useful ML inputs rather than merely intermediate calculations.

Quiet take: methodologically close to the "DFT-derived tensor information as ML signal" line. Useful as an adjacent-field citation if we later discuss why tensor features might be informative to a GNN.

Source URL:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11209944/

### Jonas, E., Kuhn, S. & Schlorer, N. (2022) "Prediction of chemical shift in NMR: A review." *Magnetic Resonance in Chemistry* 60(11), 1021-1031. DOI: 10.1002/mrc.5234.

**Status:** Metadata only; likely pending acquisition.

Review of data-driven chemical-shift prediction methods, databases, strengths, and pitfalls. Appears broader than proteins, but useful for positioning data-driven NMR prediction without over-relying on protein-specific reviews.

Quiet take: likely useful in the literature review as background, but not as central as the primary protein papers.

Source URL:
- https://pubmed.ncbi.nlm.nih.gov/34787335/

### Klukowski, P., Riek, R. & Guntert, P. (2025) "Machine learning in NMR spectroscopy." *Progress in Nuclear Magnetic Resonance Spectroscopy* 148-149, 101575. DOI: 10.1016/j.pnmrs.2025.101575.

**Status:** On disk.

Recent broad review of ML in NMR, including chemical-shift assignment, structure determination, prediction, denoising, signal detection, and repositories. Local PDF and companion summary already exist as `klukowski-2025-machine-learning-nmr-spectroscopy`.

Quiet take: useful for a high-level "what ML is doing in NMR now" citation. It should not be made to carry detailed claims about the specific predictor lineage where primary papers are available.

Local files:
- `references/klukowski-2025-machine-learning-nmr-spectroscopy.pdf`
- `references-meta/klukowski-2025-machine-learning-nmr-spectroscopy-summary.txt`

---

## Dynamics, Ensembles, and NMR/MD Time Series

### Markwick, P.R.L., Cervantes, C.F., Abel, B.L., Komives, E.A., Blackledge, M. & McCammon, J.A. (2010) "Enhanced conformational space sampling improves the prediction of chemical shifts in proteins." *Journal of the American Chemical Society* 132, 1220-1221. DOI: 10.1021/ja9093692.

**Status:** On disk.

Combines accelerated MD ensembles with SHIFTX prediction for IkappaBalpha and reports improved agreement with experimental chemical shifts at an optimal acceleration level. The paper is compact but useful because it ties chemical-shift prediction to conformational sampling and dynamics.

Quiet take: already one of the most directly useful local sources for Stage 2. It supports looking at frames/ensembles rather than a single static structure, without requiring us to make broad claims about all dynamic shift prediction.

Local files:
- `references/markwick-2010-enhanced-conformational-sampling-chemical-shifts.pdf`
- `references-meta/markwick-2010-enhanced-conformational-sampling-chemical-shifts-summary.txt`

### Lipari, G. & Szabo, A. (1982) "Model-free approach to the interpretation of nuclear magnetic resonance relaxation in macromolecules. 1. Theory and range of validity." *Journal of the American Chemical Society* 104, 4546-4559. DOI: 10.1021/ja00381a009.

### Lipari, G. & Szabo, A. (1982) "Model-free approach to the interpretation of nuclear magnetic resonance relaxation in macromolecules. 2. Analysis of experimental results." *Journal of the American Chemical Society* 104, 4559-4570. DOI: 10.1021/ja00381a010.

**Status:** Pending acquisition.

Canonical model-free relaxation papers defining generalized order parameters and effective correlation times for macromolecular internal motion. Necessary if the thesis uses S2/order-parameter language in a substantive way.

Quiet take: acquire before writing relaxation theory. These are foundations, not optional review citations.

Source URLs:
- https://pubs.acs.org/doi/10.1021/ja00381a009
- https://pubs.acs.org/doi/10.1021/ja00381a010

### Showalter, S.A. & Bruschweiler, R. (2007) "Validation of Molecular Dynamics Simulations of Biomolecules Using NMR Spin Relaxation as Benchmarks: Application to the AMBER99SB Force Field." *Journal of Chemical Theory and Computation* 3, 961-975. DOI: 10.1021/ct7000045.

**Status:** Pending acquisition.

Uses NMR spin relaxation as a benchmark for MD simulations. Relevant to any methods section that compares MD-derived motional quantities against relaxation observables.

Quiet take: this is a methods anchor for Stage 2/3 validation rather than a chemical-shift predictor paper.

Source URL:
- https://pubs.acs.org/doi/10.1021/ct7000045

### Maragakis, P. et al. (2008) "Microsecond molecular dynamics simulation shows effect of slow loop dynamics on backbone amide order parameters of proteins." *Journal of Physical Chemistry B* 112, 6155-6158. DOI: 10.1021/jp077018h.

**Status:** Open HTML available via PMC; PDF download blocked from this environment. Manual fetch recommended if this becomes more than context.

Microsecond MD study connecting slow loop dynamics with backbone amide order parameters. Useful as a scale reference for what longer simulations can reveal in NMR-relevant dynamics.

Quiet take: useful if comparing what a 25 ns trajectory can and cannot reasonably represent. Do not overuse unless the Stage 2 discussion actually needs that comparison.

Source URL:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC2805408/

### Linke, M., Kofinger, J. & Hummer, G. (2018) "Rotational Diffusion Depends on Box Size in Molecular Dynamics Simulations." *Journal of Physical Chemistry Letters* 9, 2874-2878. DOI: 10.1021/acs.jpclett.8b01090.

**Status:** Pending acquisition.

Shows that rotational diffusion from periodic MD has finite-size effects and gives a hydrodynamic correction. Relevant if relaxation/backcalculation uses rotational diffusion, tumbling time, or box-size-sensitive quantities.

Quiet take: only needed if rotational diffusion enters the analysis. If it does, this is a practical methods citation, not background ornament.

Source URL:
- https://pubs.acs.org/doi/10.1021/acs.jpclett.8b01090

### Torchia, D.A. (2015) "NMR Studies of Dynamic Biomolecular Conformational Ensembles." *Progress in Nuclear Magnetic Resonance Spectroscopy* 84-85, 14-32. DOI: 10.1016/j.pnmrs.2014.11.001.

**Status:** Open HTML available via PMC; PDF not yet fetched.

Review of NMR observables and dynamic biomolecular ensembles, including relaxation, residual interactions, chemical shifts, and MD-supported ensemble interpretation.

Quiet take: useful for a readable bridge between NMR dynamics observables and computational ensemble thinking. It should sit behind the primary relaxation and MD-validation papers.

Source URL:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4325279/

---

## Ring Current, Electric Fields, and Biomolecular Physical Features

### Sahakyan, A.B. & Vendruscolo, M. (2013) "Analysis of the Contributions of Ring Current and Electric Field Effects to the Chemical Shifts of RNA Bases." *Journal of Physical Chemistry B* 117, 1989-1998. DOI: 10.1021/jp3057306.

**Status:** On disk.

Already in the local corpus. Useful because it compares Pople and Haigh-Mallion ring-current models, unifies sign conventions, and studies the coupling of ring-current and electric-field factors in RNA base shifts. It is not a protein paper, but it is directly relevant to the physical decomposition of chemical-shift effects in biomolecules.

Quiet take: useful for the physical-feature argument, especially the point that ring-current and electric-field geometric factors can be coupled and nucleus-dependent.

Local files:
- `references/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions.pdf`
- `references-meta/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-summary.txt`

### Sahakyan, A.B., Vranken, W.F., Cavalli, A. & Vendruscolo, M. (2011) "Structure-based prediction of methyl chemical shifts in proteins." *Journal of Biomolecular NMR* 50, 331-346. DOI to verify.

**Status:** Metadata only.

CH3Shift method for methyl chemical shifts. It includes ring current, magnetic anisotropy, electric field, rotameric type, and dihedral angle effects with polynomial functions of interatomic distances.

Quiet take: likely useful because it sits between protein side-chain shifts, classical physical terms, and distance functions. Acquire if methyl/side-chain discussion grows.

Source URL:
- https://pubmed.ncbi.nlm.nih.gov/21748266/

### Williamson, M.P. & Asakura, T. (1993) "Empirical comparisons of models for chemical-shift calculation in proteins." *Journal of Magnetic Resonance, Series B* 101, 63-71. DOI: 10.1006/jmrb.1993.1008.

**Status:** Pending acquisition.

Compares Johnson-Bovey and Haigh-Mallion ring-current models, magnetic anisotropy models, and electric-field models using protein CalphaH shifts. This looks like a useful historical bridge between physical-feature models and later empirical predictors.

Quiet take: likely worth acquiring because it directly compares the same families of terms that appear in the extractor.

Source URL:
- https://www.sciencedirect.com/science/article/pii/S1064186683710083

### Case, D.A. (2013) "Chemical shifts in biomolecules." *Current Opinion in Structural Biology* 23, 172-176. DOI: 10.1016/j.sbi.2013.01.007.

**Status:** Pending acquisition.

Short review linking chemical shifts, empirical predictors, quantum calculations, ring-current/magnetic susceptibility terms, electric-field effects, and conformational dynamics.

Quiet take: useful as a compact review for committee-readable framing. It should not replace primary sources for equations.

Source URL:
- https://www.sciencedirect.com/science/article/abs/pii/S0959440X13000237

---

## Structure Prediction and AI/NMR Context

### Klukowski, P. et al. (2023) "Time-optimized protein NMR assignment with an integrative deep learning approach using AlphaFold and chemical shift prediction." *Science Advances* 9, eadi9323. DOI to verify.

**Status:** Open HTML available via PMC; Science direct PDF returned 403 from this environment. Manual fetch recommended if needed.

Integrates ARTINA with AlphaFold and UCBShift to reduce experimental assignment burden. Uses structure and predicted shift distributions as optional inputs to assignment workflows.

Quiet take: useful context for current AI-assisted NMR workflows. It is probably not central to the extractor argument, but it helps place chemical-shift prediction in the practical NMR pipeline rather than only as a benchmark problem.

Source URL:
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10664993/

### Klukowski, P. et al. (2022) "Rapid protein assignments and structures from raw NMR spectra with the deep learning technique ARTINA." *Nature Communications* 13, 6151. DOI: 10.1038/s41467-022-33879-5.

**Status:** In `references/incoming/` as `klukowski-2022-artina-raw-nmr-deep-learning.pdf`; not yet ingested.

ARTINA uses deep learning for automated protein NMR assignment and structure determination from spectra. Its graph-based chemical shift refinement step uses BMRB fragment statistics and graph structure.

Quiet take: adjacent rather than central. Useful if describing the broader modern NMR/AI ecosystem.

Source URL:
- https://www.nature.com/articles/s41467-022-33879-5

### Fowler, N.J. & Williamson, M.P. (2022) "The accuracy of protein structures in solution determined by AlphaFold and NMR." *Structure* 30, 925-933. DOI to verify.

**Status:** Metadata only.

Appears in the PLM-CS reference list. Potentially relevant to the solution-NMR structure versus AlphaFold/crystal-structure mismatch discussion.

Quiet take: acquire only if the thesis needs to talk about AlphaFold structures as inputs to shift predictors or as comparison structures for NMR-derived ensembles.

Source trail:
- Cited by Zhu et al. 2025 PLM-CS.

---

## Acquisition / Ingestion Priorities from This Pass

Downloaded into `references/incoming/` on 2026-04-26 and verified with `file`/`pdfinfo`:

- `li-2020-ucbshift-real-world-data.pdf`
- `zhu-2025-plm-cs-sequence-chemical-shifts.pdf`
- `yang-white-2021-gnn-chemical-shifts.pdf`
- `han-2011-shiftx2-protein-chemical-shift-prediction.pdf`
- `larsen-2015-procs15-dft-chemical-shift-predictor.pdf`
- `bratholm-2017-procs15-qm-chemical-shielding-refinement.pdf`
- `paruzzo-2018-shiftml-molecular-solids.pdf`
- `klukowski-2022-artina-raw-nmr-deep-learning.pdf`

Manual fetch list from this pass:

- **UCBShift 2.0** — ACS PDF returned 403; PMC page blocked by reCAPTCHA from this environment. Fetch via Birkbeck/ACS or browser-accessible PMC author manuscript.
- **SPARTA+** — PMC has HTML and a PDF link, but direct PDF returned challenge/HTML; Springer endpoint returned article HTML, not PDF.
- **CSTShift** — PMC/ACS full text available, but direct PDF from PMC returned challenge/HTML. Fetch manually if the tensor-descriptor GNN angle is pursued.
- **VASCO** — PMC full text available, Wiley PDF returned 403 and PMC direct PDF returned challenge/HTML.
- **ARECA** — PMC full text available, Springer endpoint returned article HTML and PMC direct PDF returned challenge/HTML.
- **Maragakis 2008** — PMC full text available, direct PDF returned challenge/HTML.
- **Klukowski 2023 ARTINA+AlphaFold** — Science Advances PDF returned 403; PMC full text should be browser-accessible.
- **RefDB and CheckShift** — Springer endpoints returned article HTML, not PDF.
- **Lipari-Szabo 1982 I/II, Showalter-Bruschweiler 2007, Linke-Kofinger-Hummer 2018, Williamson-Asakura 1993, Case 2013** — likely need institutional access.

1. **UCBShift 2.0**: current state-of-field comparator for side-chain and whole-protein prediction.
2. **UCBShift 2020**: real-world data and filtering discussion; directly useful for BMRB/RefDB caution.
3. **ProCS15**: DFT-based protein chemical-shift predictor; possible tensor data availability.
4. **Yang/White GNN 2021**: differentiable graph predictor near the planned GNN stage.
5. **PLM-CS 2025**: sequence-only ESM contrast; useful but not central.
6. **SHIFTX2 and SPARTA+**: baseline predictor papers; needed because every modern comparison points at them.
7. **RefDB, VASCO, ARECA**: data-quality and validation foundations.
8. **Lipari-Szabo 1982 I/II and Showalter-Bruschweiler 2007**: acquire before writing relaxation/order-parameter sections.
9. **ShiftML and CSTShift**: adjacent-field support for ML against quantum labels and tensor descriptors.
10. **Williamson-Asakura 1993 and Case 2013 review**: physical-feature lineage and compact review support.

---

## Notes for Later Processing

- When PDFs are acquired, place them in `references/incoming/` and ingest with `scripts/references/ingest_pdf.sh <pdf_path> <basename>`.
- Suggested basenames should follow the existing pattern, e.g. `ptaszek-li-2024-ucbshift-2-side-chain-prediction`, `li-2020-ucbshift-real-world-data`, `yang-white-2021-gnn-chemical-shifts`, `zhu-2025-plm-cs-sequence-chemical-shifts`.
- Do not migrate any item from this note into claim-bearing thesis prose until its full text is in `references/`, companion files exist, and the relevant claims have been checked against the local text/PDF.
