# Planned Calculators — Session 0 Ideas (2026-04-22 pm)

**Status:** Note, not spec. Five calculator / diagnostic / model-
architecture ideas surfaced during the Session 0 literature pass
for `spec/PHYSICS_FOUNDATIONS.md`. Captured here so the ideas are
preserved through subsequent sessions without becoming crib-sheet
gospel. Re-read against the current pipeline and bibliography
before implementing any of them. None are committed design.

**Provenance:** all five came out of the 2026-04-22 pm bibliography
round-outs for communities 2, 3, 4, 6, 7 plus the 0.7 queries. Each
entry names the bibliography origins so a future session can trace
the reasoning back to source.

**Posture per rollup 13.4:** these are shielding-from-maths
extensions of the existing kernel catalogue and adjacent analysis,
*not* NMR-surface features. Ideas 1 and 2 are genuinely new
calculators that fit the existing `ConformationResult` /
`TrajectoryResult` architecture. Idea 3 is a zero-new-calculator
diagnostic pass over existing results. Ideas 4 and 5 are Stage 3
model-architecture neighbours.

**Companion docs:**
- `references/ANNOTATED_BIBLIOGRAPHY.md` — bibliography entries
  cited by each idea below.
- `spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md` section 13.2
  (non-shielding NMR observables → pipeline handles table) and
  13.4 (posture directive).
- `spec/PHYSICS_FOUNDATIONS.md` section 0.9 (external validation
  benches table) and progress log.
- `GEOMETRIC_KERNEL_CATALOGUE.md` — existing kernels these
  calculators reuse or extend.
- `OBJECT_MODEL.md` — `ConformationResult` / `TrajectoryResult`
  patterns.

---

## 1. `GreenKuboSpectralDensityResult` (new, trajectory-level)

### What
Per-atom spectral density `J(ω)` at a named set of Larmor
frequencies (600 / 800 / 950 MHz typical), computed as the
Fourier-transform pair of the per-atom σ(t) autocorrelation
carried on `TrajectoryResult`.

### Why
`J(ω)` at Larmor frequencies is the direct input to R1 / R2 / R1ρ
relaxation predictions per rollup 13.2 NMR-observable table. The
time-domain handle is already the planned
`TrajectoryResult.sigma_tensor_autocorrelation`; `J(ω)` is the
frequency-domain dual at no additional physics cost — the Fourier
transform of an autocorrelation is the spectral density by
Wiener–Khinchin.

### Type
`TrajectoryResult` field (per-atom, per-irrep, discretised over a
named frequency set) or a post-trajectory `ConformationResult`
subclass. Placement belongs with the trajectory-level work in
section 5 of the rollup, not yet begun.

### Dependencies
- `TrajectoryResult.sigma_tensor_autocorrelation` (planned in
  rollup section 5)
- A TOML-configurable list of Larmor frequencies

### Origin
- M27 Kubo 1957 — linear response / FDT
- **M31 Green 1954** — the "Green" in Green-Kubo
- M28 Berne-Harp 1970 — TCF-in-MD formalism
- M10 Tugarinov-Kay 2003, M11 Ferrage 2008 — CCR-rate observables
- rollup W3 (`nh_dipolar_csa_ccr_rate`) — the downstream consumer

### Cost estimate (rough)
If σ(t) is already an array of `SphericalTensor` per atom per
frame, `J(ω)` is one FFT per atom per irrep component per Larmor
frequency. Negligible vs the autocorrelation itself.

### Risk
Convergence. `J(ω)` at low ω requires long time windows; at 50 ns
our window is too short for slow tumbling (τ_c often > 50 ns).
Honesty caveat: usable only for fast regimes. Document what's
trustworthy and what requires the μs-harvester phase.

---

## 2. `PseudocontactShiftResult` (new, conformation-level)

### What
Pseudocontact shift (PCS) per atom given an external χ tensor
(from paramagpy fit) anchored at a lanthanide-tag position.
Reuses the K_ab dipolar kernel already implemented for
`RingSusceptibilityResult`, with source = lanthanide position and
anisotropy direction = eigenvectors of χ.

### Why
Directly produces a classical-physics prediction of the PCS
observable for proteins tagged with lanthanide chelators at known
positions — e.g. ubiquitin S57C, GB1 Q32C. Paramagpy fits the χ
tensor from experimental PCS (M26); feeding that tensor back into
our pipeline closes the loop and gives a per-atom classical
prediction with an external ground truth. First-class external
validation handle per rollup 13.4 posture.

### Type
New `ConformationResult` class, patterned after
`RingSusceptibilityResult`. Full rank-2 tensor output (Mat3 +
`SphericalTensor`) per the minimum representation contract.

### Dependencies
- Existing K_ab dipolar kernel machinery
- External inputs per conformation:
  - `lanthanide_position: Vec3` (Angstroms)
  - `chi_tensor: Mat3` (from paramagpy fit, units per paramagpy)
- `SpatialIndexResult`
- `GeometryResult` (atom positions)

### Origin
- **M24 Bertini-Luchinat-Parigi 2002** — PCS formalism (identical
  kernel to RingSusceptibility, different source)
- M25 Nitsche-Otting 2021 — modern PCS review
- M26 Orton-Huber-Otting 2020 — paramagpy, the χ-tensor fitter
- N5 (Aime-Barge 2018) — ubiquitin S57C with known Δχ
- N6 (Joss-Häussinger 2019) — ubiquitin + hCA II benchmarks
- N7 Tharayil et al. 2021 — phosphoserine lanthanide sites on
  ubiquitin and GB1 with published χ tensors

### Precedent
`MopacCoulombResult` and `MopacMcConnellResult` follow the same
"same kernel, different source data" pattern. PseudocontactShift
extends that pattern from QM-derived charges and bond orders to
paramagnetic magnetic susceptibility.

### Risk
Motional averaging of the tag. Flexible lanthanide tags partially
average the PCS, introducing an effective distance dependence
different from the rigid-tag prediction. L12 Loth 2005 and M24
Bertini 2002 both discuss. Our prediction is the rigid-tag limit;
the comparison to experiment reveals tag flexibility.

### Status
Low implementation risk: kernel exists, interface is analogous to
an existing calculator, validation benches have published χ
tensors.

---

## 3. Per-secondary-structure CSA stratification (diagnostic pass)

### What
Post-prediction analysis: stratify the predicted σ tensor per-atom
by DSSP 8-class secondary structure, reporting per-stratum mean
and variance of principal components (σ₁₁, σ₂₂, σ₃₃) and of the
irrep decomposition (T0, T2).

### Why
L7 Yao-Bax 2010 reports ¹⁵N CSA magnitudes of −173 ± 7 ppm in
α-helix vs −162 ± 6 ppm in β-sheet — a ~11 ppm systematic split.
If our prediction reproduces this split without being tuned for
secondary structure, that is a thesis finding. If the residual is
large specifically at the α/β transition, the failure points at a
specific T2 angular breakdown — backbone kernel interaction with
H-bonding topology differs between α and β.

### Type
**Not a new calculator.** Python-side analysis over existing
results (DsspResult + OrcaShieldingResult + per-kernel shielding
contributions). Lives with the learn/ calibration pipeline.

### Dependencies
- `DsspResult` (exists)
- `OrcaShieldingResult` per-atom (exists)
- Per-kernel `shielding_contribution` fields (exist)

### Origin
- **L7 Yao-Bax 2010** — the α/β ¹⁵N CSA split target
- L6 Wylie 2011 — full-tensor per-residue CSA
- L12 Loth-Pelupessy-Bodenhausen 2005 — ubiquitin amide CSA
  per-residue
- B12-B13 Havlin-Oldfield surveys — ab initio precedent for
  per-stratum expectation
- MATHS_GOALS pillar 2 (T2-residual diagnostic)

### Status
Implementable whenever someone opens the learn/ calibration
notebook with the existing fleet calibration. Trivial code.
Thesis-finding impact depends on what the numbers say.

---

## 4. Lie-group GP regression on SE(3) (Stage 3 model neighbour)

### What
Alternative to equivariant-GNN for learning the kernel-feature →
shielding-tensor map: Gaussian process regression with a heat /
Matérn kernel on SO(3) or SE(3) as the input space. Provides
natural per-point posterior variance.

### Why
For a thesis that measures T2 angular residuals, per-residual
uncertainty quantification is worth considering. GP output comes
with posterior mean AND variance by construction; e3nn / ridge
currently produce point estimates only. A model neighbour, not
a replacement.

### Type
**Stage 3 model architecture alternative**, not a classical
calculator. Lives with the learn/c_equivariant/ (or a sibling
learn/c_gp/) pipeline.

### Dependencies
- Existing kernel features from the extraction pipeline
- `geometric-kernels` Python library or equivalent
- Per-atom input points embedded into SE(3) / SO(3) configuration
  space (ring frames, bond frames)

### Origin
- **E24 Azangulov-Borovitskiy et al. 2024** — stationary kernels
  on Lie groups
- F6 Chirikjian — SE(3) Gaussian processes in robotics
- F7 Helgason — harmonic analysis on Lie groups foundation

### Risk
GP scaling is O(N³) in training points; 10⁹ atom-frames over the
full trajectory fleet is infeasible for a full GP. Sparse GP or
inducing-point variants (Boots-group sparse GP on matrix Lie
groups) required. May be practical only on per-atom-type
subsets.

### Status
Exploratory. Not on the near-term Stage 3 roadmap. Flagged for
consideration when uncertainty quantification becomes important
to a specific thesis claim.

---

## 5. FNO for volumetric shielding field (Stage 3 / UI neighbour)

### What
A Fourier Neural Operator learning the map from protein structure
(Cartesian input) to a volumetric shielding field (Cartesian
output), trained against classical-kernel field evaluations on a
grid.

### Why
The h5-reader viewer renders BS / HM butterfly isosurfaces via
closed-form `QtBiotSavartCalc` / `QtHaighMallionCalc`. For large
proteins, volumetric evaluation at high resolution is expensive.
An FNO trained once from classical evaluations could render faster
at per-use time while preserving the calibrated physics. This is
a viewer / h5-reader performance concern, **not** a thesis-claim
concern.

### Type
UI / viewer-adjacent tool, not a pipeline `ConformationResult`.

### Dependencies
- Existing classical kernel evaluators (`QtBiotSavartCalc`,
  `QtHaighMallionCalc`) as training-data generators
- A neural-operator library (e.g., `neuraloperator/neuraloperator`)
- Cartesian grid sampling of the volumetric field

### Origin
- **E21 EGNO 2024** — equivariant temporal conv in Fourier space
  (architectural ancestor)
- E27 Li FNO 2021 — canonical FNO paper
- E28 Lu DeepONet 2021 — architectural cousin

### Risk
Training data is cheap to generate but the quality bar is set by
the closed-form math itself — an FNO that matches BS / HM to
visual-render precision may not add value over the closed form.
Worth only if the closed-form evaluation is a performance
bottleneck in practice.

### Status
Speculative. Not on any roadmap.

---

## How to use this document

- Before implementing any of the above, re-read against the
  current state of `OBJECT_MODEL.md`, `spec/CONSTITUTION.md`, and
  `GEOMETRIC_KERNEL_CATALOGUE.md`. Ideas here are frozen
  2026-04-22; the pipeline is not.
- Ideas 1 and 2 are genuinely new calculators that fit the
  existing architecture. Lowest implementation risk.
- Idea 3 is a zero-new-calculator diagnostic pass that could run
  against existing fleet calibration data.
- Ideas 4 and 5 are Stage 3 model-architecture neighbours.
  Evaluate alongside the equivariant-GNN Stage 3 design when that
  work begins.
- None of these are committed design. The note exists so the ideas
  are not lost between sessions; the user decides whether any
  graduate to a full spec pass.

## Amendments

*(None yet. Append amendments below, never rewrite above.)*
