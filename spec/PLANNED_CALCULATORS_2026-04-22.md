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

### Amendment 2026-05-08 — `PlanarGeometryResult` (new, conformation-level)

#### What
Per-frame `ConformationResult` companion to the `LegacyAmber` substrate's
typed planar-group / ring-position fields. Per-conformation values for:

- ω (peptide-bond planarity dihedral) per peptide bond.
- Ring-flip state per Ring.
- sp2 pyramidalization per planar-group atom.
- Ring-pucker phase (Pro pyrrolidine and other saturated rings).

The substrate carries the typed *classification*
(`PlanarGroupKind`, `PlanarStereo`, `RingPosition` — landed in
the 2026-05-05 → 2026-05-08 topology slice). This calculator carries
the actual *deviation from canonical* in each frame.

#### Why
The substrate names the chemistry ("this atom is in a peptide-amide
planar group, expected E-stereo"); only the conformation knows whether
the actual ω deviates from 180° in this frame, or whether a ring is
flipped. Stratification analyses (per-SS CSA, mutation-shift mechanism,
ProCS15-class rebuilds) need both pieces. ProCS15 is the cautionary
tale: fixing ω = 180° in the substrate scan loses 5–14% of backbone
shielding signal. Capturing per-frame deviation makes the calculator–
substrate pair complete and is the load-bearing reason the substrate
slice locked in `PlanarStereo` / `PlanarGroupKind` enums.

#### Type
`ConformationResult` subclass. Per-frame field on `ConformationAtom`
for sp2 pyramidalization; per-bond / per-residue companion store for ω;
per-`Ring` field for flip + pucker.

#### Dependencies
- `GeometryResult` (positions, derived planes / dihedrals).
- `LegacyAmberTopology::HasAtomSemantic()` (substrate must be present;
  `PlanarGroupKind` / `PlanarStereo` / `RingPosition` drive which atoms
  / bonds / rings this calculator visits).

#### NPY emission (planned)
Per-frame planar-geometry NPYs paired with the per-protein
`atoms_category_info.npy`:

- `omega_actual.npy`, `omega_deviation.npy` — per peptide bond
- Ring-flip state and pucker phase per ring
- Per-atom sp2 pyramidalization

Schema details preserved in the topology-substrate implementation plan
(retired to `spec/plan/bones/` 2026-05-08; that doc is the prose
source for the design intent above).

#### Origin
Topology slice (Bundle B/C, 2026-05-05 → 2026-05-08). The substrate
side landed (`LegacyAmberSemanticTables`, `atoms_category_info.npy`,
typed identity matchup); the per-frame conformation companion did NOT.
Captured here so the design intent survives the bones retirement.

#### Cost estimate (rough)
One pass over residues for ω (per-peptide-bond dihedral, four atoms);
one pass over rings for flip / pucker (Cremer–Pople or simpler proxy);
one pass over planar-group atoms for sp2 pyramidalization (out-of-plane
distance from the three substituents' plane). All conformation-only,
no expensive electronic-structure work. Should fit comfortably in the
existing per-frame `OperationRunner::Run` budget.

#### Risk
Low. The math is well-trodden (ω is a standard Ramachandran companion
angle; Cremer–Pople pucker has been in MD analysis for decades). The
substrate's typed driver (PlanarGroupKind / RingPosition) means the
calculator does not need to recompute "is this atom in a peptide
amide?" — the substrate already says so.

#### Status
**PENDING.** Substrate landed in commits `ba647cd` (atom_semantic
store), `6ec9bff` (Bundle C ring substrate), `8accdb6`
(CategoryInfoProjection NPY). Per-frame `PlanarGeometryResult` not yet
implemented. Earlier audit docs (now in bones/) sometimes referred to
this as "Phase 1 companion landed" — the substrate side did, the
per-frame Result did not. This entry is the canonical pending pointer.

---

### Amendment 2026-05-08(b) — `AIMNet2PolarisabilityResult` (new, conformation-level)

#### What
Per-atom charge sensitivity computed by a **single autograd backward
pass** through the TorchScript AIMNet2 model:
`d(q_i)/d(r_i)` as a per-atom 3-vector (gradient of atom i's predicted
charge with respect to its own position), reduced to a scalar magnitude
or kept as a vector — both useful for stratification.

This is the autograd path documented in `src/AIMNet2Result.cpp` lines
465–490 ("Charge sensitivity: autograd path (future)") that has been
PENDING since the AIMNet2 calculator landed.

#### Why
Polarisability is a real physical dimension in shielding prediction —
Stage 1 dimension inventory shows charge polarisation contributes
**+0.197 R²** for carbon, the dominant element-specific gap between
geometry-only and charge-aware models. The deleted random-bulk-
perturbation approach was correctly retired (random splats are non-
comparable across frames and lie about solvation when applied to other
conformations). The autograd path replaces it with a deterministic,
per-conformation, physically meaningful quantity at a fraction of the
cost.

The substrate landed in 2026-05-08 (`atoms_category_info.npy`) gives
codex per-atom typed identity columns, so polarisability-by-atom-type /
by-ring-position / by-protonation-state stratification becomes a single
boolean-mask operation for the thesis chapter.

#### Type
`ConformationResult` subclass. Per-atom field on `ConformationAtom` for
`aimnet2_polarisability_vector` (Vec3) and / or
`aimnet2_polarisability_scalar` (norm).

#### Dependencies
- `AIMNet2Result` (uses the same loaded model; Compute either runs the
  forward pass with `requires_grad=True` on coordinates, OR caches the
  required intermediate state at AIMNet2Result::Compute time and runs
  backward separately).
- `Session::Aimnet2Model()` (model lifetime).

#### NPY emission (planned)
- `aimnet2_polarisability.npy` (N, 3) — per-atom gradient vector
- `aimnet2_polarisability_scalar.npy` (N,) — norm of the per-atom
  gradient vector
  
Both indexed by atom row, parallel to existing `aimnet2_charges.npy`.

#### Pre-flight blocker — `.jpt` requires_grad support

The TorchScript model export at
`/shared/2026Thesis/nmr-shielding/data/models/aimnet2_wb97m_0.jpt`
must propagate gradients through the coordinate input tensor. This is
**unverified**. The calculator should NOT be written until a 10-line
Python check confirms:

```python
import torch
m = torch.jit.load("data/models/aimnet2_wb97m_0.jpt")
coords = torch.randn(1, N, 3, requires_grad=True)
# … run forward …
charges.sum().backward()
assert coords.grad is not None and coords.grad.norm() > 0
```

If the check passes, the C++ calculator is ~50 lines (forward without
`NoGradGuard`, single backward, copy gradients off CPU). If the check
fails, the AIMNet2 weights need to be re-exported from the original
PyTorch with grad-tracking enabled — that's a model-asset task
upstream of the calculator.

#### Cost estimate (rough)

Per-frame:
- Forward pass on 4000-atom protein: ~1.4 s (linear extrapolation from
  the polarisability-roadmap's 0.17 s @ 479 atoms).
- Backward pass: 2–3× forward typical for graph-NN-with-attention →
  ~3–5 s.
- Total: ~5–6 s per frame.

Aggregate for plausible thesis-chapter scopes:

| Target | Frames | Hours at 5 s/frame |
|---|---|---|
| 1 representative protein × 25 ns (1250 frames) | 1.25K | ~1.7 h |
| 10 calibration proteins × 25 ns | 12.5K | ~17 h (overnight) |
| 685-protein fleet × DFT pose set (26 frames each) | 17.8K | ~25 h |
| Full 685-protein fleet × every frame | 857K | ~50 days (not affordable) |

The 2-page thesis chapter budget is the 1-protein or 10-protein scope —
overnight on one GPU.

#### Origin
- `spec/POLARISABILITY_ROADMAP_2026-04-13.md` is the broader
  five-approaches landscape (items 1, 2, 4 landed: SasaResult surface
  normal, HydrationGeometryResult, EeqResult; items 3, 5 remain).
- `src/AIMNet2Result.cpp:465-490` carries the autograd-path TODO with
  the design described verbatim.
- Original (now-deleted) random-bulk-perturbation approach was retired
  on 2026-04-12 (perturbation method deleted) and the SDK slot cleaned
  up in commit `7390a6d` (2026-05-05). Earlier estimates of the
  polarisability cost were from per-atom × 3-direction × 2-sided
  perturbation (~6N forward passes per frame, blew up by ~4 orders of
  magnitude); the autograd path is **one** backward pass per frame
  total — orders of magnitude cheaper, hence the recovery of this
  calculator as a tractable thesis-chapter target.

#### Risk
Medium. The `.jpt` autograd-support unknown is the only structural
blocker. If it works, the calculator is straightforward; if it
doesn't, the project needs an upstream model regeneration. Either
way, the FIRST step is the 10-line Python check, NOT writing C++.

#### Status
**PENDING.** Pre-flight `.jpt` requires_grad check is the gating
investigation. Calculator implementation, NPY emission, SDK wrapper
all downstream of that check.

---

### Amendment 2026-05-08(c) — Illustrative test peptides for workbook chapters

When the calculator-walkthrough chapters land in the thesis, two test
peptides serve as illustrative substrates for the Mathematica
workbooks that visualise each calculator's per-atom output. Pinning
the choice here so the calculator work picks the same substrates and
the workbook chapters all share a comparable geometric reference.

#### Primary: Trp-cage (TC5b)

Sequence `NLYIQWLKDGGPSSGRPPPS`, 20 residues. Real folded mini-protein;
PDB `1L2Y` (NMR ensemble, Neidigh / Andersen 2002) and `2JOF` (high-
resolution variant). Covers Trp (fused indole rings), Tyr (aromatic
+ hydroxyl), Pro (saturated 5-ring + no backbone H — four prolines,
unusually rich), Gly (HA2/HA3 inversion case), Asp / Lys / Arg
(charged sidechains), Asn / Gln / Ser (polar amides + hydroxyl),
Leu / Ile / Val (β-branched + methyl super-groups).

Used as the **primary illustrative substrate**: folded geometry,
heavily NMR-validated, established benchmark in MD/folding/shielding
literature, citation-friendly. Workbook figures from a real folded
peptide carry more weight than an extended-conformation synthetic.

Coverage gap: A, C, E, F, H, M, T, V are absent (8 of 20 canonical
residues). The chemistry features driving those (methyl, disulfide,
carboxylate, simple aromatic 6-ring, imidazole variants, sulfur,
β-hydroxyl, branched β) are partly compensated by Trp's indole +
Tyr's aromatic-OH + Asp/Asn carboxyl/amide. But for full demonstration
of e.g. PHE-vs-TYR ring-current contrasts, or HID/HIE/HIP variant
chemistry, Trp-cage isn't enough.

#### Companion: synthetic all-canonical peptide

Designed 22-residue peptide for chemistry features Trp-cage doesn't
reach. Suggested sequence: `ACDEFGHIKLMNPQRSTVWYCV` — one of each
canonical residue, plus a second Cys (so the CYX disulfide variant
has a partner) and a second Val (so the methyl super-group QG
aggregator has two methyl-pair instances).

Linear extended conformation built via `tleap` from the sequence; no
fold required because the workbook chapters need the chemistry, not
the structure. The PDB and topology files commit to
`tests/data/illustrative_peptides/` (or wherever the calculator work
chooses) so all workbooks load the same reference geometry.

Generation deferred until the calculator-walkthrough work begins.
~5-line tleap script; the geometry is stable across rebuilds because
tleap deterministically applies ff14SB ideal angles + the standard
extended-chain φ/ψ.

#### Variant peptides

HID / HIE / HIP, ASH, GLH, LYN, ARN, TYM, CYX-as-disulfide-pair are
mutually exclusive within a single residue — a HIS is one variant at
a time. Don't try to cram all variants into one peptide. The natural
shape is **multiple peptide instances**: swap one residue in either
Trp-cage or the synthetic peptide, re-run the calculator, compare
outputs side-by-side. The categorical record from
`atoms_category_info.npy` makes "this is the HID instance, this is
the HIE instance" trivially identifiable in the join.

---

### Amendment 2026-05-09(a) — `CSAPrincipalAxisResult` (new, conformation-level)

#### What
Per-atom diagonalisation of the σ tensor. Stores principal components
σ₁₁, σ₂₂, σ₃₃ (sorted by convention) and the orthonormal principal-
axis-system (PAS) orientations as a Mat3 of column eigenvectors per
atom. Source σ tensor is either `OrcaShieldingResult` (where DFT is
present) or the calibrated kernel sum (everywhere else); the calculator
is source-agnostic — it consumes whatever Mat3 σ tensor is attached.

#### Why
The full Mat3 + SphericalTensor representation is preserved end-to-end
per the T2-is-sacred contract, but downstream comparisons to the
literature (L7 Yao-Bax 2010 per-stratum ¹⁵N CSA magnitudes, L12
Loth-Pelupessy-Bodenhausen 2005 per-bond ubiquitin tensors,
L6 Wylie 2011 full-tensor per-residue CSA) are quoted in PAS terms.
Storing the diagonalisation alongside the lab-frame Mat3 makes the
literature comparison a direct table join rather than an on-the-fly
re-diagonalisation per consumer. Rollup W3 CCR derivation also wants
σ-in-PAS as input.

#### Type
`ConformationResult` subclass. Per-atom field on `ConformationAtom`
for `csa_principal_components` (Vec3, sorted) and
`csa_principal_axes` (Mat3, columns are eigenvectors aligned to the
sorted principal components).

#### Dependencies
- A populated σ tensor on `ConformationAtom` (either
  `OrcaShieldingResult` or the calibrated kernel sum, when a
  runtime calibration path lands).
- `LegacyAmberTopology` substrate present so per-stratum analyses
  (per-element, per-`BackboneRole`, per-`RingPositionLabel`) can join
  the PAS output against typed identity columns from
  `atoms_category_info.npy`.

#### NPY emission (planned)
- `csa_principal_components.npy` (N, 3), float64 — sorted σ₁₁ ≤
  σ₂₂ ≤ σ₃₃ per atom (sign convention per CONSTITUTION).
- `csa_principal_axes.npy` (N, 3, 3), float64 — column-eigenvector
  Mat3 per atom; column k aligns with `csa_principal_components[k]`.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.1` — the section
  this amendment migrates.
- L6 Wylie 2011 (full-tensor per-residue CSA), L7 Yao-Bax 2010 (α/β
  ¹⁵N CSA split), L11 Brender-Taylor-Ramamoorthy 2001 (amide
  orientation conventions), L12 Loth-Pelupessy-Bodenhausen 2005
  (ubiquitin per-bond principal components) — all in
  `references/ANNOTATED_BIBLIOGRAPHY.md`.

#### Cost estimate (rough)
One symmetric-3×3 eigendecomposition per atom per frame. Negligible
relative to any kernel evaluation; ~µs per atom on commodity CPU.
Per-frame cost for a 4000-atom protein is dominated by memory traffic,
not arithmetic.

#### Risk
Low. Numerical convention (sort order, eigenvector sign-fixing across
near-degenerate components) needs to be locked in CONSTITUTION when the
slice lands; otherwise downstream consumers will see arbitrary sign
flips on near-axially-symmetric tensors. The literature quoting
convention (Haeberlen vs Mehring vs Herzfeld-Berger) must be picked
once and documented in the H5 sign-convention metadata
(C.6 in the new I/O schema doc).

#### Status
PENDING. Migrated from DESIDERATA §A.1 as part of the 2026-05-09
doc-cleanup pass. Slice gated on (a) the calibrated kernel-sum σ
tensor being a stable per-frame field and (b) sign-convention pick
in CONSTITUTION.

---

### Amendment 2026-05-09(b) — `AmideTensorGeometryResult` (new, conformation-level)

#### What
Per-amide bond, expresses the per-atom σ tensor (carried on N and
HN) in the local N-H frame following the Brender-Taylor-Ramamoorthy
2001 convention (L11). One Mat3 per amide bond representing σ in the
local frame, plus the rotation Mat3 that takes lab-frame σ to the
local frame. The amide bond inventory is driven by typed substrate:
peptide-bond-amide planar groups identified by
`AtomSemanticTable::planar_group_kind == PlanarGroupKind::PeptideAmide`
(or whatever the typed enum spelling is at slice time), not by atom
name matching.

#### Why
13.2's R1/R2/R1ρ handle wants σ-in-bond-frame, not σ-in-lab-frame.
CCR back-calculation (A10) consumes σ-in-amide-frame as input. Storing
this once at calculator time means downstream relaxation predictions
join against a stable typed field rather than re-rotating the lab-frame
tensor per consumer. The L11 convention is the documented standard
for amide-CSA work; concretising it as a typed field aligns the
pipeline with the convention used by experimentally-anchored
literature comparisons (L12 Loth ubiquitin amide CSA per residue).

#### Type
`ConformationResult` subclass. Per-amide-bond record (one entry per
backbone N-H pair that is part of a `PeptideAmide` planar group):
- `amide_sigma_local_frame` (Mat3) — σ tensor expressed in the local
  N-H frame
- `amide_lab_to_local` (Mat3) — rotation taking the lab-frame Mat3
  on the N atom to the local frame

Storage layout follows existing per-bond / per-residue companion
patterns (one parallel store keyed by residue index, with a "no
amide" sentinel for the N-terminal residue and prolines).

#### Dependencies
- `GeometryResult` (positions of N, HN, and the carbonyl C of the
  preceding residue, defining the amide plane).
- Per-atom σ Mat3 on N and HN (kernel-sum calibrated or
  `OrcaShieldingResult`).
- `LegacyAmberTopology::SemanticAt(ai)` for the `PlanarGroupKind`
  driver; `BackboneRole` for finding the N / HN / preceding-C atom
  triple via typed lookup, not string matching.

#### NPY emission (planned)
- `amide_sigma_local_frame.npy` (R, 3, 3), float64 — per residue
  index, with NaN-filled entries where the residue has no peptide-amide
  group (N-terminus, prolines).
- `amide_lab_to_local.npy` (R, 3, 3), float64 — paired rotation Mat3
  per residue index, NaN-filled where absent.

R = number of residues. Pairs cleanly with `atoms_category_info.npy`'s
per-residue row index for downstream join.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.1` — the section
  this amendment migrates (A.1 captures both A.1 CSAPrincipalAxis and
  A.2 AmideTensorGeometry).
- L11 Brender-Taylor-Ramamoorthy 2001 — the orientation convention
  the calculator implements.
- L12 Loth-Pelupessy-Bodenhausen 2005 — ubiquitin amide CSA
  per-residue dataset for validation against this calculator's output.
- A.10 `CCRRateResult` (Amendment 2026-05-09(h)) — direct downstream
  consumer.

#### Cost estimate (rough)
Per amide bond per frame: one 3×3 rotation construction (cross product
+ Gram-Schmidt) + one Mat3 conjugation (R σ Rᵀ). Few µs per bond.
Total per-frame cost for a 200-residue protein is in the millisecond
range; negligible.

#### Risk
Low. Convention pick must be locked once (Brender-Taylor-Ramamoorthy
local-frame axis definitions: x along N-H, z normal to the amide
plane) and documented in the H5 sign-convention metadata. The
"no amide here" sentinel pattern (N-terminus, prolines) needs the
same NaN-fill discipline as other absent-by-chemistry fields on the
substrate side.

#### Status
PENDING. Migrated from DESIDERATA §A.1 (A.2 sub-item) as part of the
2026-05-09 doc-cleanup pass. Slice gated on (a) the calibrated σ
tensor field landing per atom, (b) Amendment 2026-05-09(a)
CSAPrincipalAxis convention pick, (c) typed `PlanarGroupKind` driver
exposing the peptide-amide inventory cleanly.

---

### Amendment 2026-05-09(c) — `BulkSusceptibilityAccumulator` (new, conformation-level)

#### What
Protein-level (single-record-per-conformation) magnetic susceptibility
χ tensor accumulated from per-atom contributions. Sums McConnell
local-anisotropy contributions, ring-current susceptibility from
`RingSusceptibilityResult`, and H-bond susceptibility from `HBondResult`
to produce a 3×3 χ tensor for the whole protein, plus its
diagonalisation (principal Δχ) and SphericalTensor decomposition. Babaei
et al. 2017 (N3) hierarchical-subunit approach is the prose reference
for the accumulation strategy.

#### Why
Methodologically novel output: the project's DFT-calibrated per-atom
kernels rolled up to a bulk observable measured by a different
literature (peptide Δχ via Pauling 1979 N2; bulk-protein anisotropy
via STI MRI per Li 2017 N4). No other group produces this specific
bridge from calibrated classical kernels to MRI-scale susceptibility.
First-class thesis-claimable output per rollup section 13 synergy 3.
Validates the per-atom kernel set against an integrated independent
measurement, strengthening the "kernels carry the signal" claim
beyond per-atom DFT residual comparison.

#### Type
`ConformationResult` subclass storing **a single record per
conformation**, not a per-atom field. Fields:
- `bulk_chi_tensor` (Mat3) — accumulated lab-frame χ tensor
- `bulk_chi_principal` (Vec3) — sorted principal Δχ values
- `bulk_chi_axes` (Mat3) — column-eigenvector PAS
- `bulk_chi_irrep` (SphericalTensor) — T0/T1/T2 decomposition

#### Dependencies
- `McConnellResult` (per-atom local anisotropy contributions).
- `RingSusceptibilityResult` (per-ring current susceptibility).
- `HBondResult` (per-H-bond susceptibility contribution).
- `LegacyAmberTopology` substrate so the accumulator can mask
  contributions by `RingPositionLabel`, `PolarHKind`, etc. when
  the methods chapter wants per-stratum decomposition of the bulk
  observable.

#### NPY emission (planned)
- `bulk_chi_tensor.npy` (3, 3), float64 — per-conformation Mat3.
- `bulk_chi_principal.npy` (3,), float64 — sorted principal Δχ.
- `bulk_chi_axes.npy` (3, 3), float64 — PAS Mat3.
- `bulk_chi_irrep.npy` (9,) or shape-matching SphericalTensor layout
  (per the established irrep convention C.5).

One-record-per-conformation NPYs; the trajectory roll-up handles
per-frame stacking via the standard trajectory layout.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.2` — the section
  this amendment migrates.
- N2 Pauling 1979 — peptide Δχ external bench.
- N3 Babaei et al. 2017 — hierarchical subunit accumulation
  approach.
- N4 Li 2017 — STI / MRI bulk-protein anisotropy bench.
- B18 — methodological context for χ accumulation.

#### Cost estimate (rough)
Linear in atoms + rings + H-bonds; the contributing per-atom kernels
are already evaluated. Accumulation is a single sum-and-diagonalise
pass per frame. Negligible vs the contributing kernel evaluations.

#### Risk
Medium. The accumulation rule (which per-atom kernels contribute to
χ, with what sign/scaling, and how H-bond pairs vs ring currents
double-count avoidance is handled) is a methods-paper-grade
decision. Methodological-novelty risk: if the accumulator says one
thing and the experimental Δχ literature says another, the residual
itself becomes a thesis discussion item — that's the positive case;
the discipline is to not retrofit the accumulation rule to chase the
external bench. Also: per-frame χ assumes zero-velocity-correlation
between contributions; in practice some pair-correlation correction
may be needed (literature has variants).

#### Status
PENDING. Migrated from DESIDERATA §A.2 as part of the 2026-05-09
doc-cleanup pass. Slice gated on (a) all three contributing
calculators (McConnell, RingSusceptibility, HBond) being per-frame
stable and producing comparable-units outputs, and (b) accumulation-
rule pick locked once with a literature-cite trail.

---

### Amendment 2026-05-09(d) — `NICSProbeEvaluator` (new, conformation-level)

#### What
Evaluate the existing K_ab kernel family (BS, HM, McConnell, Coulomb,
HBond) at **arbitrary non-atom probe positions**, not just at the
existing `ConformationAtom` positions. Probes are typed input: one
record per probe carries position (Vec3), source-atom-set list (which
atoms contribute), and a probe-kind tag. Default probe inventory
covers NICS(0) at ring centroids and NICS(1)zz at ±1 Å along the ring
normal (Schleyer 1996 B14); user-supplied probes (E2 GeigerSmidt-class
virtual nodes, viewer-volumetric-grid points) attach via the same
contract.

#### Why
**Cross-cutting unification.** The h5-reader currently re-implements
volumetric kernel evaluation in
`QtBiotSavartCalc` / `QtHaighMallionCalc`, parallel to the library's
per-atom evaluation. NICSProbeEvaluator makes the library the single
source of truth for kernel evaluation at any point in space; the
viewer becomes a probe-position consumer rather than a parallel
re-implementation. NICS observables are a literature-standard
aromaticity diagnostic (Schleyer 1996); having them as first-class
output enables direct comparison with QM-NICS literature for ring-
current calibration. Pairs with the FNO volumetric calculator hint
(this document section 5) as training-data generator.

#### Type
`ConformationResult` subclass with per-probe records (a list, not a
per-atom field on `ConformationAtom`). Each record stores the probe
metadata (position, kind, source-set) plus the per-kernel evaluation
output (Mat3 + SphericalTensor for the irrep-decomposable kernels).

#### Dependencies
- All currently-implemented K_ab kernel evaluators (BS, HM, McConnell,
  Coulomb, HBond), used as functions over `(probe_position,
  source_atom_list)`.
- `RingTopology` from `LegacyAmberTopology` — drives the default
  ring-centroid + ring-normal probe generation via typed
  `ring.centroid_atoms` / `ring.plane_normal` (no atom-name matching).
- `SpatialIndexResult` for probe-source distance queries.

#### NPY emission (planned)
- `probe_positions.npy` (P, 3) — per-probe Cartesian position.
- `probe_kind.npy` (P,) — typed probe-kind enum (e.g.
  `RingCentroidNicsZero`, `RingNormalPlusOneAng`,
  `RingNormalMinusOneAng`, `UserProbe`).
- `probe_kernel_<name>.npy` (P, 3, 3) — per-probe Mat3 per kernel.
- `probe_kernel_<name>_irrep.npy` (P, 9) — irrep decomposition.

P = number of probes. Schema integrates with C.5 irrep tagging and
C.6 units / sign metadata in the new I/O schema doc.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.3` — the section
  this amendment migrates.
- B14 Schleyer 1996 — NICS observables.
- B15 — magnetic-aromaticity context.
- E12 — virtual-node philosophy in equivariant networks.
- Calculator hint #5 in this document (FNO volumetric) — direct
  training-data consumer.

#### Cost estimate (rough)
One kernel evaluation per probe per kernel per frame; cost scales
like (P × N) for each K_ab kernel where N is source-atom count,
identical scaling to per-atom evaluation but with P = O(rings) for
default NICS coverage. For a 4000-atom protein with ~40 rings → 120
probes (centroid + ±1 Å), cost is ~3% of the per-atom evaluation
time per kernel. Negligible.

#### Risk
Low to medium. Architectural risk: the probe-evaluation API design
must be re-usable across all current K_ab kernels without per-kernel
duplication; this is the B.1 "Unified `KernelSource` hierarchy"
DESIDERATA item being implicitly required. If the unified-source
slice isn't in place, each kernel evaluator needs a parallel "evaluate
at this Vec3" overload — workable but architectural debt.

#### Status
PENDING. Migrated from DESIDERATA §A.3 as part of the 2026-05-09
doc-cleanup pass. Slice gated on (a) probe-evaluation API design,
ideally in concert with the B.1 unified-source hierarchy (in
DESIDERATA Section B / migration target
`comprehensive-calculator-inventory-2026-04-30.md` Section 10),
(b) probe-kind enum lock-in.

---

### Amendment 2026-05-09(e) — `ParamagneticRelaxationEnhancementResult` (PRE) (new, conformation-level)

#### What
Per-atom paramagnetic relaxation enhancement: R1 / R2 contribution
from a paramagnetic centre (typically a lanthanide tag) via 1/r⁶
geometry — Solomon-Bloembergen-Morgan (SBM) formalism. Different
observable from the PCS in Section 2 of this document
(`PseudocontactShiftResult`); shared input surface (lanthanide-tag
position, electronic relaxation time τ_s). Pairs with PCS — both
calculators consume the same external-input contract and validate
against the same Tharayil et al. 2021 (N7) ubiquitin / GB1 datasets.

#### Why
PCS and PRE together provide two independent paramagnetic observables
from one tag position. PRE is dominated by 1/r⁶ near the tag; PCS is
1/r³ with anisotropy direction from χ tensor. The pair gives both a
distance (PRE) and an oriented-distance (PCS) constraint, and
literature benches publish both observables side-by-side per tagged
mutant. Closing both in the pipeline gives a stronger external
validation handle than either alone, per rollup 13.4 first-class
external bench posture.

#### Type
`ConformationResult` subclass. Per-atom field on `ConformationAtom`
for `pre_r1` (scalar, Hz contribution) and `pre_r2` (scalar, Hz
contribution). Optionally `pre_distance_to_tag` (scalar, Å) as a
geometric companion field for sanity plots.

#### Dependencies
- Per-conformation external input: `lanthanide_position` (Vec3,
  Angstroms) and electronic relaxation time `tau_s` (seconds, per-tag
  literature value).
- `GeometryResult` for atom positions.
- `SpatialIndexResult` for tag-to-atom distance lookup at scale.
- The same input surface as A5 PCS (Section 2 of this document) — see
  `PseudocontactShiftResult` for the input contract template.

#### NPY emission (planned)
- `pre_r1.npy` (N,), float64 — per-atom R1 enhancement (Hz).
- `pre_r2.npy` (N,), float64 — per-atom R2 enhancement (Hz).
- `pre_distance_to_tag.npy` (N,), float64 — per-atom distance (Å)
  to the tag position. Kept as separate NPY for diagnostic plots
  (the 1/r⁶ rolloff inspection).

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.4` — the section
  this amendment migrates (A.4 captures both A.5 PCS — covered in
  Section 2 of this document — and A.6 PRE).
- M24 Bertini-Luchinat-Parigi 2002 — paramagnetic-NMR formalism
  including SBM derivation.
- M25 Nitsche-Otting 2021 — modern PCS/PRE review.
- N5 Aime-Barge 2018, N6 Joss-Häussinger 2019, N7 Tharayil et al.
  2021 — published ubiquitin/GB1 benchmark datasets with paired
  PCS+PRE observables.
- Pairs with Section 2 of this document (`PseudocontactShiftResult`,
  Amendment 2026-04-22 corpus) — same input surface, complementary
  observable.

#### Cost estimate (rough)
Per atom per frame: one distance computation + scalar SBM evaluation.
Cheaper than any K_ab kernel sum (no source-atom loop; only one
source = the tag). Negligible.

#### Risk
Low to medium. **Tag motion**: flexible tags partially average the
PRE just like the PCS; SBM with a static tag position is the rigid-
limit prediction. Honesty caveat for the methods text. **Order-of-
magnitude sanity**: PRE near-zone enhancements can be huge (kHz at
sub-5 Å), and the literature reports both R1 and R2 contributions —
double-check sign and unit conventions during slice. Tag-electronic-
relaxation-time τ_s is per-tag and a literature-lookup, not a
property of the protein.

#### Status
PENDING. Migrated from DESIDERATA §A.4 (A.6 sub-item) as part of the
2026-05-09 doc-cleanup pass. Pairs with PCS (this document §2). Slice
gated on the same lanthanide-tag external-input contract that PCS
needs — when one lands, the other should follow in the same slice.

---

### Amendment 2026-05-09(f) — `MemoryKernelExtractionResult` (new, trajectory-level)

#### What
Zwanzig projection-operator memory kernel extracted from per-atom
σ(t) autocorrelation. Complements `GreenKuboSpectralDensityResult`
(this document §1) on the non-Markovian side: where J(ω) is the
frequency-domain Wiener-Khinchin dual of the autocorrelation, the
memory kernel is the time-domain non-Markovian-correction handle
extracted via the Volterra equation from the same autocorrelation.

#### Type
`TrajectoryResult` field. Per-atom, per-irrep, on a discretised
time-lag grid (same grid as the autocorrelation). Consumed by
analysis scripts; not a per-frame field on `ConformationAtom`.

#### Dependencies
- `TrajectoryResult.sigma_tensor_autocorrelation` (the planned
  trajectory-scope autocorrelation field; same input as J(ω)).
- The C.8 trajectory-schema formalisation (lag-grid layout,
  per-irrep ordering) — must be settled before this lands.

#### NPY emission (planned)
- `memory_kernel.npy` (N, irrep, lag), float64 — per-atom
  per-irrep kernel on the lag grid.

Schema lag layout follows the same convention picked in C.8 for
all autocorrelation-derived trajectory NPYs.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.5` — the section
  this amendment migrates (A.5 also covers A.7
  `GreenKuboSpectralDensityResult`, already captured in this
  document §1).
- M29 Zwanzig 1961 — primary projection-operator paper.
- M22 Zwanzig 2001 — textbook treatment.
- M21 Kou-Xie — power-law-memory evidence in biomolecular dynamics.

#### Cost estimate (rough)
Per atom per irrep: one numerical Volterra-equation solve over the
lag grid given the autocorrelation. Cheap (linear in lag count) once
the autocorrelation is computed. Storage is the same per-atom,
per-irrep, per-lag block as the autocorrelation itself.

#### Risk
Medium. **Convergence honesty**: the memory kernel from a finite-
window autocorrelation is intrinsically noisy at long lags; document
trustworthy-lag-range explicitly per the same discipline as the
J(ω) low-frequency caveat in this document §1. **Numerical stability
of Volterra inversion**: ill-conditioned at low signal-to-noise lags
— pick a regularised solver and document the regularisation choice
in CONSTITUTION when the slice lands.

#### Status
PENDING. Migrated from DESIDERATA §A.5 (A.8 sub-item) as part of the
2026-05-09 doc-cleanup pass. Slice gated on (a) the autocorrelation
field landing as a stable `TrajectoryResult` field, (b) the C.8
trajectory-schema formalisation lag-grid convention pick.

---

### Amendment 2026-05-09(g) — `ErgodicityMetricResult` (new, trajectory-level, diagnostic)

#### What
Per-atom Thirumalai ergodic measure plus block-averaged standard-
deviation per `TrajectoryResult` field. Formalises the H7-H13
convergence discipline as first-class output. One entry per
trajectory (not per frame); attaches a single companion record per
analyzed trajectory carrying per-atom + per-field convergence
diagnostics.

#### Why
**Honesty-caveat machinery for every other trajectory-level field.**
Rollup 13.3 captures the project posture: "shipped by named default,
not rule-coverage" — convergence diagnostics ride alongside every
trajectory output by default rather than as an opt-in analysis
script. Without this, J(ω), CCR rates, memory kernels, and any other
trajectory-derived observables look more confident than they are.
Thirumalai's measure plus block-SD is the literature-standard pair:
the measure quantifies how far from ergodic the sampling is (relative
to the ergodic limit), block-SD quantifies the practical uncertainty
on whatever derived observable was computed from that sampling.

#### Type
`TrajectoryResult` companion (one record per trajectory). Per-atom
ergodic measure plus per-field block-SD. Diagnostic, not consumed by
downstream calculators.

#### Dependencies
- A complete trajectory of per-atom σ tensors (the same time series
  the autocorrelation is computed from).
- A list of "fields to measure" — every per-atom trajectory-derivable
  scalar / tensor field. Implementation-time decision: is the field
  list discovered from the `TrajectoryResult` schema directly, or
  declared by name in TOML?

#### NPY emission (planned)
- `ergodic_measure.npy` (N,), float64 — per-atom Thirumalai measure
  for the σ tensor time series.
- `block_sd_<field>.npy` (N,), float64 — per-atom block-SD for each
  named trajectory-derived field. One NPY per consumed field.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.5` — the section
  this amendment migrates (A.5 also covers A.7 J(ω) and A.8 memory
  kernel; A.9 ergodicity is the diagnostic-style sibling).
- H7-H13 — convergence-discipline thread in the bibliography.
- Rollup 13.3 — "shipped by named default" posture directive.

#### Cost estimate (rough)
Per atom: one Thirumalai-measure pass over the time series (linear
in trajectory length). Per (atom, field): one block-SD pass over a
block-size grid. Total cost is linear in (atoms × frames), comparable
to the autocorrelation cost itself.

#### Risk
Low. The metrics are well-understood; the project risk is **scope
creep** — the temptation to inflate this into a full convergence-
analysis subsystem rather than the disciplined "honesty caveat
emitted by named default" output. The discipline is: ship it, plot
it, don't make it the centre of the chapter.

#### Status
PENDING. Migrated from DESIDERATA §A.5 (A.9 sub-item) as part of the
2026-05-09 doc-cleanup pass. Slice gated on (a) the trajectory-scope
σ time-series field being stable, (b) the field-list-to-measure
contract being settled (TOML-driven vs schema-discovered).

---

### Amendment 2026-05-09(h) — `CCRRateResult` (new, trajectory-level)

#### What
Auto- and cross-correlated relaxation rates (DD/DD and DD/CSA) per
amide-bond per atom. Computed from the per-atom σ tensor (in the
amide local frame, A.2 — Amendment 2026-05-09(b) above), the per-
frame bond vector, and a Larmor-frequency input. Consumes
`GreenKuboSpectralDensityResult` (this document §1) as intermediate.
Directly produces the experimental observable in L12 Loth 2005's
ubiquitin 64-bond dataset (rollup W3 handle).

#### Why
**Thesis-grade external validation: the project's top-end NMR
observable.** CCR rates are a literature-quoted observable (L12 Loth
2005 publishes 64 amide bonds with experimental DD/CSA cross-rates
on ubiquitin); back-calculating them from our pipeline closes the
loop kernel → calibrated σ → bond-frame projection → autocorrelation
→ J(ω) → CCR. Each upstream is independently testable; the CCR rate
is the integrated thesis-claim observable. Consuming A.2 amide-frame
σ + A.7 J(ω) + per-frame bond vector means CCR is the natural
endpoint of the trajectory-scope chain.

#### Type
`TrajectoryResult` field. Per-amide-bond record (one per
peptide-amide planar group, driven by typed substrate per A.2).
Stores DD/DD auto-rate, DD/CSA cross-rate, plus diagnostic
intermediates (which J(ω) values were consumed at which Larmor
frequency).

#### Dependencies
- A.2 `AmideTensorGeometryResult` (Amendment 2026-05-09(b)) for σ
  in the local N-H frame.
- A.7 `GreenKuboSpectralDensityResult` (this document §1) for J(ω)
  at the requested Larmor frequencies.
- Per-frame bond-vector time series (already implicit in
  `GeometryResult` over the trajectory; possibly a derived
  `TrajectoryResult` field).
- TOML-configured Larmor-frequency list (shared with A.7).

#### NPY emission (planned)
- `ccr_dd_dd.npy` (R,), float64 — per residue DD/DD auto-rate (Hz).
- `ccr_dd_csa.npy` (R,), float64 — per residue DD/CSA cross-rate
  (Hz).
- `ccr_diagnostic_jw_consumed.npy` (R, F), float64 — per residue,
  per Larmor frequency, the J(ω) value consumed (for traceability).

R = residues; F = Larmor-frequency count. NaN-fill for residues
with no peptide-amide group.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.6` — the section
  this amendment migrates.
- M10 Tugarinov-Kay 2003, M11 Ferrage 2008, M12, M13 — CCR-rate
  formalism.
- L12 Loth-Pelupessy-Bodenhausen 2005 — ubiquitin 64-bond
  experimental dataset.
- Rollup W3 — `nh_dipolar_csa_ccr_rate` handle.
- Pairs with Amendment 2026-05-09(i) `BenchmarkBackCalculationResult`
  — Loth-2005-ubiquitin-CCR is one of the named benches.

#### Cost estimate (rough)
Per amide bond per frame: a few scalar combinations of (σ-in-local,
J(ω), bond-vector). Negligible. The cost is in the upstream
autocorrelation → J(ω); CCR itself adds a thin combination layer.

#### Risk
Medium. **Convention-dense**: every literature group has a slightly
different CCR rate definition (sign of the cross-rate, prefactor on
the dipolar coupling, choice of Larmor frequency for cross-J(ω)).
Pin one convention (recommend M11 Ferrage 2008) and document in
CONSTITUTION + sign-convention metadata when the slice lands.
**Validation honesty**: matching L12 Loth's 64 bonds is the headline
thesis claim — discipline-of-evidence applies, "correlate not match"
is the framing (per project memory).

#### Status
PENDING. Migrated from DESIDERATA §A.6 as part of the 2026-05-09
doc-cleanup pass. Slice gated on (a) Amendment 2026-05-09(b) amide
local-frame σ landing, (b) GreenKuboSpectralDensityResult (§1)
landing, (c) per-frame bond-vector trajectory field, (d) CCR-rate
convention pick.

---

### Amendment 2026-05-09(i) — `BenchmarkBackCalculationResult` (new, per-bench conformation-level)

#### What
Per **named validation bench**, compute the project pipeline's
prediction and store it alongside the published experimental value
for direct residual comparison. One attach per bench; queryable by
bench name. Initial roster:

- `loth_2005_ubiquitin_ccr` — DD/CSA cross-rates per amide bond
  (consumes Amendment 2026-05-09(h) `CCRRateResult`).
- `yao_bax_2010_15n_csa_split` — α/β ¹⁵N CSA stratification
  (consumes Amendment 2026-05-09(a) CSA principal components +
  DSSP).
- `babaei_2017_bulk_chi` — bulk Δχ per protein (consumes Amendment
  2026-05-09(c) `BulkSusceptibilityAccumulator`).
- `tharayil_2021_pcs` — PCS per atom per tagged mutant (consumes
  PCS, this document §2).
- `wylie_2006_d22_co_hn` — δ22 carbonyl + amide-H shift bench
  (consumes the calibrated kernel-sum σ).

Turns the `PHYSICS_FOUNDATIONS.md` §0.9 validation-benches table
into concrete, consumed, tested pipeline output.

#### Why
The validation-benches table already exists as project intent; the
gap is that nothing in the pipeline materialises (prediction,
experimental, residual) per bench in a queryable form. Per-bench
H5 slices (C.7 in the new I/O schema doc — see below) consume this
result; per-bench thesis plots become a single H5 group open per
chapter. Discipline of evidence: "this is what we predicted; this
is what the experiment says; this is the residual" — one record
per atom per bench, no further bookkeeping.

#### Type
`ConformationResult` subclass; one **per-bench attach** per
conformation. Pattern: bench name + per-atom-or-per-bond record
storing `prediction`, `experimental`, `residual`, plus per-bench
metadata (units, sign convention, citation). Multiple attaches
allowed (one per bench in the consumed roster).

#### Dependencies
- Per-bench upstream calculator landing (see roster above; each
  bench cites its consumer Amendment).
- `ExperimentalReferenceLoader` (DESIDERATA C.2; planned for
  consolidation into the new `spec/I_O_AND_SCHEMA_2026-05-09.md`
  doc) — provides the typed experimental-reference store keyed by
  bench name.
- `LegacyAmberTopology` substrate so per-bench atom matching uses
  typed `(residue_index, AtomMechanicalIdentity)` joins (the
  project discipline since the 2026-05-06–08 slices), not string
  comparison against legacy bench-published atom names.

#### NPY emission (planned)
Per-bench NPY group named `bench_<name>_*`:
- `bench_<name>_prediction.npy`, `bench_<name>_experimental.npy`,
  `bench_<name>_residual.npy` — same per-atom-or-per-bond record
  shape per bench.
- `bench_<name>_meta.json` — citation, units, sign convention,
  matched-atom count, unmatched-atom list.

The C.7 "Validation-bench H5 slices" item (DESIDERATA §C.2,
migration target = new I/O schema doc) consumes these NPYs into
H5 groups `/benches/<name>/`.

#### Origin
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md §A.7` — the section
  this amendment migrates.
- `PHYSICS_FOUNDATIONS.md §0.9` — the validation-benches table
  this amendment makes concrete.
- Each bench cites its own literature: L12 Loth 2005, L7 Yao-Bax
  2010, N3 Babaei 2017, N7 Tharayil 2021, L6 Wylie 2006/2011.
- Pairs with C.2 `ExperimentalReferenceLoader` and C.7
  `validation-bench H5 slices`, both consolidated into the new
  `spec/I_O_AND_SCHEMA_2026-05-09.md` doc planned in the
  2026-05-09 cleanup pass.

#### Cost estimate (rough)
Per bench per conformation: one filter pass over the per-atom field
matching the bench's atom-set; one residual computation. Negligible.
Cost is dominated by the upstream calculators (CCR, BulkSus, PCS,
PRE, CSA) whose outputs this consumes.

#### Risk
Medium. **Atom-identity matching**: legacy benches publish atom names
in conventions that pre-date our typed substrate (BMRB / RefDB
remapping is the existing case study, see project memory
`feedback_huxley_data_discipline` and `project_nmr_forensics`).
Per-bench loaders must remap to typed
`(residue_index, AtomMechanicalIdentity)` carefully; unmatched-atom
lists must surface, not be silently dropped (Huxley discipline).
**Convention drift**: each bench has its own units/sign/Larmor-
frequency assumptions; per-bench metadata must lock these explicitly
so the residual is computed with matched conventions on both sides.

#### Status
PENDING. Migrated from DESIDERATA §A.7 as part of the 2026-05-09
doc-cleanup pass. Pairs with C.2 `ExperimentalReferenceLoader` and
C.7 bench H5 slices in the new
`spec/I_O_AND_SCHEMA_2026-05-09.md` doc landing in the same cleanup
pass. Slice gated on (a) per-bench upstream calculators landing, (b)
`ExperimentalReferenceLoader` typed external-data contract landing,
(c) per-bench atom-identity matching loader (one per bench).
