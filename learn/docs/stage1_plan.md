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
   per-element ridge, R²=0.818.

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
