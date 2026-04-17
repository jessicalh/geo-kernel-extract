# learn/docs/ — Document Index

**Read `stage1_plan.md` first.**

## Plan and procedure

- **stage1_plan.md** — Chapter plan, argument chain, rigour
  requirements, literature status, key findings, file map.
- **extraction_procedure.md** — How to set up and run a named
  extraction over the 723 mutant pairs.

## Analysis results (on 110 proteins — update to 723)

- **twenty_eight_realities_2026-04-10.md** — Full validation: 28
  physics realities verified analytically and tensor-directly.
- **element_physics_2026-04-10.md** — Per-element kernel decomposition
  with literature grounding.
- **calibrated_weights_2026-04-10.md** — Per-element ridge weight
  vectors and their physical interpretation.

## Thesis draft

- **realities_latex.tex** — Partial LaTeX section (28 realities).
  Needs 723-protein final numbers.

## Working notes

Session-by-session working notes live in
`stage1-mutations/notes/`.  Current contents:

- **dimension_inventory_raw.md** — Complete inventory of dimensions
  examined, including surprises (dispersion, normalization, ff14SB
  vs MOPAC, ridge weight in tail PCs, proven zeros).
- **argument_chain.md** — The central thesis claim and its 6-step
  grounding (each step cited or shown).
- **literature_grounding.md** — What's established, extended to T2,
  or novel.  Papers we have, need from PMC, need from library.
- **reproducible_analysis_plan.md** — Four-phase plan to turn the
  iterative discovery into defensible methodology: full-space
  characterisation, group LASSO, nonlinear check, stability.
- **dia_para_decomposition.md** — Validated WT-ALA delta of DFT
  dia/para components.  Dia and para cancel (~7+7→1 ppm).  Kernels
  see the net, not the channels.  Para/dia ratio confirms Saito 2010.
- **new_arrays_assessment.md** — All new Stage1Results arrays tested.
  Every new scalar adds zero.  AIMNet2 EFG orthogonal to Coulomb.
  Charges diverge for N/O.  SASA/DSSP/charges for Stage 2 not Stage 1.
- **dimensionality_honest.md** — Element-specific predictive dims:
  H=20, C=6, N=3, O=12.  Not "3 dimensions."  The element-dependent
  complexity IS the physics story.
- **physics_group_grounding.md** — For each of 8 kernel groups: the
  physical interaction, equation, what determines T2 direction, and
  what 400-protein data shows.  Viva-level.  Group-to-group cosine
  matrix interpretation.  Independence structure is element-invariant;
  contributions are element-dependent.
- **normalisation_physics.md** — Cosine matrix invariant to
  normalisation (angular structure is intrinsic).  R² changes:
  H loses (-0.06), N/O gain (+0.05/+0.06).  Dispersion explodes
  for O (0.058→0.234).  Normalisation separates "how far" from
  "which direction."
- **dimension_tree.md** — Full tree from T2 perturbation down to
  each contributing kernel, with cosines at each branch, citations,
  equations.  The tree is element-invariant; the contributions are
  element-dependent.
- **physics_analysis_bridge.md** — How physics groups connect to
  analysis dimensions.  Per-element: which groups contribute how
  many dims.  Bridging table.  Why H=20 (ring current diversity),
  C=6 (EFG dominated), N=3 (everything blurred), O=12 (dispersion
  after normalisation).
- **completeness_and_nonlinear.md** — LPOCV, bootstrap, nonlinear.
  N=+0.169, C=+0.128, O=+0.013, H=+0.002.  Nonlinear signal follows
  paramagnetic ordering (Saito 2010 / Ramsey 1950).
- **atom_type_stratification.md** — Per-atom-type ridge (CA, C=O, CB,
  C side, N bb, N side, O bb, O side).  Sidechain N = 0.887,
  carbonyl C = 0.463.  "Nitrogen is hard" was wrong.
- **master_chart.md** — The complete picture.  Per-element tree with
  all physics groups, R², cosines, dimensionality, nonlinear signal,
  dia/para, proven zeros, limitations.  Citation map at end.

## Historical

Session notes, next-session prompts, and non-Stage-1 docs are in
`docs/bones/`.
