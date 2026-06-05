# Codex — modern rank-2 tensor (T2 / NMR shielding) visualisation: survey + recommendation

## (Standard qt6-cpp + VTK-docs preamble is prepended above. PLUS explicit BROAD web-search license for this task — see below. READ-ONLY research+author: no code, no git. Cite every claim with a URL.)

## You have broad web-search license for this task
Search the open web freely — papers, tool docs, library docs/examples, conference talks — for
how rank-2 symmetric tensors (and spherical/irrep signals) are visualised, classic through
ML-era. Cite every claim with a URL; quote sources where you can. "I don't know / couldn't
find it" beats a confident guess.

## Why
The reader must visualise the **rank-2 shielding tensor (T2)** per atom — T2 is the thesis
argument and may not be collapsed to a scalar. The default plan was a classic eigenvalue-scaled
**ellipsoid glyph**, but the lead questions whether that's the best/most modern choice: "the
rise of machine learning has led to new paper visualisations, some remarkably effective."
She's right to ask — ellipsoids are perceptually weak for planar/linear anisotropy, and
equivariant-ML work renders exactly these objects (spherical tensors / irreps) very
effectively. **Crucially, this project's Stage-3 model is e3nn (equivariant), so an
e3nn-native spherical-harmonic rendering would be both modern AND thematically coherent with
the model.** Survey the field and recommend, before we build.

## Survey (answer each, cited)
1. **Classic tensor glyphs.** Ellipsoid (eigenvector axes × eigenvalue magnitudes) and its
   well-documented failure modes (ambiguity when two eigenvalues are close — linear vs planar
   reads the same). **Superquadric tensor glyphs** (Kindlmann 2004, "Superquadric Tensor
   Glyphs") as the perceptual fix; Reynolds/Haber glyphs. Cite the superquadric advantages.
2. **NMR-specific tensor surfaces.** The chemical-shielding-tensor "**ovaloid**"/shielding
   surface (σ as a function of orientation, σ11/σ22/σ33), and how sign (shielding vs
   deshielding) is shown. Tools: TensorView / TensorView for MATLAB, MagresView, EFGShield,
   Simpson/SOMA, VMD tensor plugins. How solid-state-NMR / NMR-crystallography papers draw
   shielding & EFG tensors. Cite.
3. **ML-era / equivariant / spherical-signal renderings.** How do equivariant-network papers
   visualise spherical tensors / irreps? The **spherical-harmonic "signal on S²"** rendering:
   the tensor plotted as a real function on the sphere, surface deformed by magnitude and
   coloured by sign. Look at: **e3nn** docs/examples (`e3nn.io`, the `SphericalTensor` / signal
   plotting, `s2grid`), Tensor Field Networks (Thomas 2018), spherical CNNs (Cohen 2018),
   Tensor Field / equivariant viz in recent papers. Note e3nn's own visualisation conventions
   and whether they're directly reusable. Cite.
4. **Effectiveness + practicality.** For each candidate, judge how well it conveys: (a)
   anisotropy magnitude, (b) orientation, (c) **sign** (shielding/deshielding — colour, not
   shape), (d) at **scale** (a scene of many atoms, not one hero glyph — overdraw, clutter,
   GPU cost). Note any perceptual studies.
5. **VTK-renderability.** Which are buildable in this reader's VTK pipeline at per-atom scale:
   superquadrics via `vtkSuperquadricSource` / `vtkTensorGlyph`; spherical-harmonic deformed
   spheres via a parametric/`vtkProgrammableSource`-warped sphere coloured by sign; ovaloids.
   Cite the relevant VTK classes (you have VTK-docs license) and flag cost.

## Output — `notes/T2_VISUALISATION_PRIORART_2026-06-05.md`
Per-method summary (cited), then a **RECOMMENDATION**: the 1-2 best options for a per-atom,
at-scale, VTK-rendered T2 visualisation in the reader, weighing perceptual effectiveness, the
e3nn-coherence angle, and build cost. State clearly what each costs to build and what it needs
from the data (full Mat3? eigendecomposition? the T2 spherical components we already have?).
No code, no git. This is a recommendation for the lead to choose from.
