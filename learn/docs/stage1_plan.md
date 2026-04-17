# Stage 1: Mutation Calibration — Thesis Plan

## What Stage 1 is

The thesis chapter answering: **what does an always-T2 kernel feature
extraction tell us about the physics of shielding?**  Aromatic mutant
DFTs provide concrete statistical assessment of same.

The contribution is not new physics.  It is showing that established
shielding physics — ring current, EFG, bond anisotropy — expressed as
rank-2 tensor kernels, survives quantitative statistical assessment
against DFT at protein scale (69K atoms, 110→723 proteins).  The
isotropic decomposition exists in the literature (Case 1995, Sahakyan
2013).  The T2-channel decomposition per element at this scale is new.

## Three stages of the project (context)

1. **Stage 1 (here):** Calibrate geometric kernels against DFT using
   723 WT-ALA mechanical mutants.  Per-element ridge.  Thesis chapter.
2. **Stage 2 (in progress):** Ensemble poses from GROMACS + DFT.
   Extend calibration to conformational variation.
3. **Stage 3 (future, not learn's job):** Full extraction from 685
   RefDB proteins → GNN NMR prediction model.

## Current state

- Calibration settled: R²=0.818 on 110 proteins (AzimuthalExtraction)
- 28 physics realities verified analytically + tensor-directly
- 7 R figures, partial LaTeX draft
- 3 newer analyses done (geometry-only basis, dimensionality,
  accumulation basis) — JSON output exists, no LaTeX/figures yet
- Working notes in `stage1-mutations/notes/` from 2026-04-12 session

## What happens next

### 1. Final extraction (723 proteins)

Re-run `extract.py --run FinalMutantExtraction --resume` over all 723
mutant pairs with current extractor (including non-water calculators).
Same scripts, final numbers.  The 110-protein results are structurally
correct — the physics won't change, the numbers will tighten.

### 2. Formalise the analysis

Re-run all `src/actual_physics/` scripts on the 723-protein data.
Add one new analysis: project PCA dimensions onto physics groups to
verify the "3 dimensions = 3 physical mechanisms" claim.  Write proper
R figures for dimensionality, geometry-only comparison, coefficient
structure.

### 3. Write the chapter

Working directory: `stage1-mutations/`.  Structure follows the
argument chain (see working notes).

---

## The chapter structure

The chapter is a characterisation of the instrument, not a result.
It answers: what does the tool see at each element, how many
independent angular dimensions does it resolve, and why?

### What the tool is

10 physical mechanisms expressed as rank-2 tensor kernels.  Each
mechanism has a known multipolar order (cited: Pople, Johnson-Bovey,
Buckingham, McConnell, Stone).  Different multipolar orders produce
different angular symmetries (Stone 2013 Ch. 3).  The kernels
confirm the predictions at 3 significant figures (BS -3.04, PQ -5.05).

### What the tool sees, per element

The element-dependent dimensionality IS the physics story:

- **H (20 predictive dimensions, R²=0.928):** The tool sees hydrogen
  richly.  Ring current and EFG both contribute, with 20 independent
  angular dimensions all carrying small slices of the total.  The
  prediction lives in the fine angular disagreements between kernel
  families.

- **C (6 dimensions, R²=0.562):** EFG-dominated.  Charge polarization
  contributes +0.197 (the largest gap, accessed via MOPAC).  Fewer
  angular dimensions because one mechanism dominates.

- **N (3 blurred dimensions, R²=0.380):** No single mechanism dominates.
  Five families each contribute 0.015-0.080.  The 3 dimensions are
  blurred mixtures.  The tool sees nitrogen in the corner of its eye.

- **O (12 dimensions, R²=0.382):** Dispersion-driven near-field (8 dims
  at 0-4Å).  Paramagnetic dominance (para/dia ratio 1.20).  Complex
  angular landscape.

### What the tool sees, per atom type (2026-04-15)

The element-pooled numbers hide large variation between backbone
and sidechain atoms.  Stratifying by AMBER atom name:

- **C=O (R²=0.463):** The hardest carbon.  Paramagnetic sp2 with
  C=O π* excitation.  Drags the carbon average down from 0.63-0.73
  for other carbon types.
- **C side (R²=0.729):** Sidechain carbons (CG/CD/CE/CZ).  Best
  predicted carbon type — closest to the removed ring.
- **N bb (R²=0.387):** Backbone amide nitrogen.  This is what made
  "nitrogen" look hard.
- **N side (R²=0.887):** Sidechain nitrogen (HIS/ARG/LYS).  Second-
  best atom type after hydrogen.  Extremely sensitive to mutation
  type (+0.150 from ring identity alone).
- **O side (R²=0.566):** Sidechain oxygen.  Most sensitive to ring
  identity of any atom type (+0.209 from mutation type).

The atom-type decomposition must be presented alongside the element
numbers.  See `atom_type_stratification.md` for full tables.

### Why the numbers differ

The electronic structure of each element determines which geometric
perturbations it responds to.  Diamagnetic-dominated elements (H)
respond to through-space magnetic fields (ring current) which the
wire-loop/surface-integral models capture in many angular dimensions.
Paramagnetic-dominated elements (C, N, O) respond to local electronic
structure (bond character, orbital symmetry) which the geometric
kernels capture partially (EFG, McConnell) or barely (N's blurred
multi-mechanism response).  Cited: Saito et al. 2010 for the
dia/para distinction; confirmed in the T2 delta channel by our
dia/para decomposition (para/dia variance ratio H=1.02 → O=1.20).

### What the tool cannot see

Per-protein ceiling: R²=0.81 within a protein, 0.35 across proteins.
The gap is the global electrostatic environment (protein shape, bulk
dielectric) that local geometric kernels do not capture.  This
motivates ensemble conformational sampling (Stage 2).

The charge-polarization dimension (+0.197 for C) is currently
accessible only via MOPAC (10 min/protein).  For ensemble work,
Drude-FF MD or EEQ charges provide affordable paths.

## The argument chain

Every step is either CITED (established physics) or SHOWN (our data).

1. Each mechanism has a known multipolar order.  **Cited**: Pople 1956,
   Johnson & Bovey 1958, Buckingham 1960, McConnell 1957, Stone 2013.
2. Different multipolar orders produce different angular symmetries in
   rank-2 tensors.  **Cited**: Stone 2013 Ch. 3.
3. The kernels confirm the multipolar predictions.  **Shown**: distance
   slopes at 3 significant figures (BS -3.04, PQ -5.05).
4. The kernel families are angularly independent.  **Shown**:
   |cos(BS, EFG)| = 0.684 across 58K pairs.
5. The DFT target aligns with the correct mechanism per element.
   **Shown + Cited**: H=ring current (Boyd 2002), C=EFG (Sahakyan 2013).
6. The calibrated coefficients are physical constants.  **Shown**:
   per-element ridge, weighted R²=0.718 (720 proteins, 446K atoms).
7. Nonlinear signal follows paramagnetic ordering.  **Shown + Cited**:
   N=+0.169, C=+0.128, O=+0.013, H=+0.002. Ramsey 1950, Saito 2010.

## Status (2026-04-13)

**Analysis: complete.**  All scripts run on 720 proteins.  All JSONs
current.  verify_numbers.py passes.  Completeness checks done (LPOCV,
bootstrap, RF nonlinear).  13 working notes in stage1-mutations/notes/.
master_chart.md has the full picture.

**Writing: not started.**  The notes are working documents, not thesis
prose.  Writing dive planned for ~2 months from now.  Will need:
Ramsey 1950 for the nonlinear argument, per-ring-type detail in the
dimension tree, Mathematica vetting of ridge algebra, R figures from
stage1_figures.R.

**This work is WRAPPED.**  stage1-mutations/ is frozen for the thesis.
New analysis (ensemble, MD, Stage 2) goes in new directories.  The
numbers here must continue to pass verify_numbers.py.

## Lessons for Stage 2

The following findings from Stage 1 directly inform ensemble analysis:

1. **Per-element everything.**  Never pool.  The physics is
   element-dependent.  The dimensionality, the dominant groups, the
   nonlinear structure, the overfitting behaviour — all differ by
   element.

2. **Normalisation is physics.**  Per-protein normalisation separates
   magnitude from angular structure.  For ensemble data, per-frame
   normalisation will play the same role: strip frame-to-frame
   magnitude variation to reveal conformational angular signal.

3. **Carbon needs charge polarisation.**  +0.197 gap.  Gating
   recovers +0.128 (nonlinear from geometry alone).  The remaining
   ~0.07 needs EEQ or Drude charges.  AIMNet2 EFG is orthogonal
   (cos 0.34) — wrong projection for mutation delta, may be
   different for conformational variation.

4. **Nitrogen is nonlinear.**  Per-element gated model for N will
   find signal that ridge misses.  5 mechanism families interact
   through the paramagnetic term.

5. **Dispersion is real for O after normalisation.**  0.058→0.234.
   DispChi (dispersion × ring susceptibility) carries angular
   structure at short range.  For ensemble analysis, this dimension
   should be preserved.

6. **The cosine independence structure is fixed.**  It doesn't
   depend on the protein, the element, or the normalisation.
   The geometry of the kernel space is intrinsic.  Only the
   projection onto the target changes.

## Literature status

Full assessment in `stage1-mutations/notes/literature_grounding.md`.
Summary:

**Established** (cite, don't derive): mechanism physics, element-
dependent diamagnetic/paramagnetic dominance (Saito 2010), site-
specific tensor variability (Hall & Fushman 2006).

**Extended to T2** (cite isotropic + note extension): per-element
mechanism dominance, calibrated coefficients, BS-HM agreement.

**Novel** (frame carefully): T2 angular separation between kernel
families, per-element T2 decomposition at scale, eigenspectrum
structure, ridge weight orthogonal to variance PCs, dispersion as
significant T2 contributor for N/O.

### Papers we have (PDFs in references/)

- Buckingham 1960, Case 1995, Boyd & Skrynnikov 2002,
  Sahakyan & Vendruscolo 2013

### Open access (download from PMC in browser)

- Saito et al. 2010 — PMC2905606
- Facelli 2011 — PMC3058154
- Hall & Fushman 2006 — PMC2519110

### Need from library

- Pople 1956.  DOI: 10.1063/1.1742701
- Johnson & Bovey 1958.  DOI: 10.1063/1.1744645
- McConnell 1957.  DOI: 10.1063/1.1743676
- Haigh & Mallion 1979.  Prog. NMR Spectrosc. 13, 303.
- Stone 2013.  Theory of Intermolecular Forces, 2nd ed.  OUP.
  ISBN 978-0199672394.  Need Ch. 3.
- Yao et al. 2010.  DOI: 10.1021/ja104711v
- Han et al. 2011.  DOI: 10.1007/s10858-011-9478-4
- Xu & Case 2002.  DOI: 10.1002/bip.10276
- Haeberlen 1976.  High Resolution NMR in Solids.  Academic Press.
  ISBN 978-0120255610.

---

## Rigour requirements

### If we talk about it, we cite it

Every assertion about physics must be either:
- A citation to established literature, OR
- An empirical result from our data (with the analysis script named)

No handwaving.  No "it is well known that."  The specific paper, the
specific equation, the specific data table.

### Per-element, always

Never pool elements.  Every result is reported per element.  If
something only works for hydrogen, say "works for hydrogen."  The
pooled view hides the physics (demonstrated: pooled R²=0.385,
weighted per-element=0.718, gap=0.333).

### Distinguish established, extended, and novel

The chapter must be clear about which claims are:
- **Established**: citing known results
- **Extended**: taking known isotropic results to the T2 channel
- **Novel**: new findings from our analysis

### Honest about limitations

- Per-protein ceiling: 0.81 within, 0.35 across
- Near-field accuracy below 4 Angstroms
- HIE is the worst ring type (self-fit R²=0.062)
- 110 proteins until final extraction completes
- The "3 dimensions correspond to 3 mechanisms" claim needs the
  PC-to-physics-group projection (not yet done)

### Reproducible

Every number traces to a named script in `src/actual_physics/` and
a config in `src/calibration.toml`.  The reproduction commands are
in `src/actual_physics/OVERVIEW.md`.

---

## Key findings from 2026-04-12 analysis session

These are in `stage1-mutations/notes/` (3 files).  Headlines:

1. **Dispersion is bigger than expected.** Enters forward selection at
   rank 2-5 for C, N, O.  Group R² for O = 0.272 (higher than ring
   current 0.231 or EFG 0.165 in geometry-only basis).

2. **Per-protein normalization reveals angular structure.** Raw kernels:
   top 10 PCs = 99.3% variance, predictive plateau k=3.  Normalized:
   top 10 PCs = 24%, plateau k=46 for H.  Scale factors carry +0.06
   R² uniformly.

3. **ff14SB vs MOPAC EFG is a proper dimension.** Cosine similarity
   0.50 for H (nearly orthogonal), 0.61 for C, 0.90 for N.  MOPAC
   wins by angular direction, not by adding independent signal after
   geometry.

4. **Ridge weight lives orthogonal to dominant variance PCs.** For H,
   top 10 variance PCs carry 0% of ridge weight.  Prediction is in
   the angular fine structure between kernel families, not in the
   shared geometric information.

5. **Proven zeros.** MOPAC valence, bond order, molecular dipole each
   add +0.000 for every element.  HBond = 0.002.  DeltaAPBS = 0.005.

6. **Dia/para cancellation** (2026-04-13).  Parsed raw orca dia/para
   from .out files, validated to 0.003 ppm against C++.  Mutation
   perturbs dia and para by ~7 ppm each; they cancel to ~1-2 ppm
   total.  Kernels see the net (R²=0.70-0.73 for H) but not the
   channels individually (R²<0.05).  Para/dia variance ratio
   1.02→1.06→1.07→1.20 for H→C→N→O confirms Saito 2010 in T2.

7. **AIMNet2 EFG is orthogonal to Coulomb EFG** (2026-04-13).  Cosine
   ~0.34 (near random in 5D) despite R²=0.53 for H.  ff14SB↔MOPAC
   EFG cosine = 0.99 (nearly identical on raw aromatic).  AIMNet2
   sees different physics but adds nothing beyond MOPAC for prediction.

8. **Spherical tensor audit** (2026-04-13).  Full audit of isometric
   normalization convention across C++, Python SDK, and analysis code.
   No inconsistencies found.  Convention is T2[m=-2..+2] with sqrt(2)
   and sqrt(3/2) factors matching Types.cpp.

---

## File map

### Working notes (from 2026-04-12 session)
- `stage1-mutations/notes/dimension_inventory_raw.md`
- `stage1-mutations/notes/argument_chain.md`
- `stage1-mutations/notes/literature_grounding.md`

### Analysis scripts (run these on 723 data)
- `src/actual_physics/OVERVIEW.md` — describes all 7 scripts
- `src/actual_physics/element_physics.py`
- `src/actual_physics/twenty_realities.py`
- `src/actual_physics/tensor_realities.py`
- `src/actual_physics/export_for_r.py`
- `src/actual_physics/clean_calibration.py`
- `src/actual_physics/physics_calibration.py`
- `src/actual_physics/per_element_calibration.py`
- `src/actual_physics/geometry_only_basis.py`
- `src/actual_physics/dimensionality_test.py`
- `src/actual_physics/accumulation_basis.py`

### Existing output (on 110 proteins — replace with 723)
- `src/output/actual_physics/` — JSON results from newer analyses
- `src/output/secondary/` — CSV/JSON from element_physics, realities

### Existing thesis material (update with 723 numbers)
- `docs/realities_latex.tex` — partial LaTeX (28 realities)
- `R/twenty_eight_realities.R` — 7 publication figures

### Full reference
- `docs/twenty_eight_realities_2026-04-10.md` — validation document
- `docs/element_physics_2026-04-10.md` — per-element decomposition
- `docs/calibrated_weights_2026-04-10.md` — weight vectors
- `EXPERIMENTS.md` — experiment log
- `references/ANNOTATED_BIBLIOGRAPHY.md` — 80+ papers catalogued
