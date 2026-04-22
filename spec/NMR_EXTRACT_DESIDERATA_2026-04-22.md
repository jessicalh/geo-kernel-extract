# nmr-extract Desiderata — Session 0 consolidated pass (2026-04-22 pm)

**Status:** Note, not spec. Consolidates every calculator,
variation, input/output surface, diagnostic, and architectural
idea surfaced during the 2026-04-22 Session 0 literature pass for
`spec/PHYSICS_FOUNDATIONS.md`. Captured here *before* any
organisation so implicit triage doesn't discard options that the
user or a future session might want to build.

**Posture (user, 2026-04-22 pm, close-to-verbatim):** "In many
respects the time we spend building calculators is 'free' compared
to the time required to run stuff. We just have to organise it."
The cost model for this project is run-time on protein fleets and
DFT queues, not developer-time on library extensions. A useful-to-
anyone calculator that costs a few days to write and a few μs per
atom-frame to evaluate is a net win whenever its output sees a
single downstream consumer. This document does not pre-triage.
Organisation and sequencing come next, against a complete list —
this one.

**Cost-model detail — MOPAC (user, 2026-04-22 pm):** MOPAC
(PM7 + MOZYME) runs ~45 s per ~889-atom protein per conformation —
expensive vs frame cadence, cheap vs DFT. **Not workable per frame
across a 25 ns trajectory** (~1200 frames); workable only on sparse
selected frames (the DFT pose set, checkpoints, user-selected
snapshots). Implication for every item below: trajectory-level
autocorrelation / covariance / CCR-rate / spectral-density fields
(A7-A10, plus the existing rollup section 5 fields) must be
driven by the per-frame-dense ff14SB + APBS classical kernel
stack, not by MOPAC. MOPAC-derived kernels (`MopacCoulombResult`,
`MopacMcConnellResult`) provide sparse "ground-truth-adjacent"
signals on the 260 DFT pose set and on future μs-harvester
selected poses. **The calibration hope** is that those sparse
MOPAC signals and the dense ff14SB classical signals agree on
the pose set — if they do, the per-frame ff14SB kernels carry
enough of the MOPAC signal that the dense trajectory outputs
are trustworthy; if they diverge, a sparse-to-dense MOPAC
correction at the trajectory level becomes an additional design
question. See diagnostic D5 below.

**Companion docs:**
- `spec/PLANNED_CALCULATORS_2026-04-22.md` — five ideas already
  captured in more detail there, marked `[P]` below when
  referenced.
- `references/ANNOTATED_BIBLIOGRAPHY.md` — every tether below
  points into specific entries there.
- `spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md` — 13.2 handles
  table, 13.4 posture directive, 13.5 kernel-kernel covariance
  signal, section 5 `TrajectoryResult` scaffold.
- `spec/PHYSICS_FOUNDATIONS.md` — 0.9 validation benches table,
  Session 0 progress log.
- `GEOMETRIC_KERNEL_CATALOGUE.md` — existing kernels this doc
  extends or reuses.
- `OBJECT_MODEL.md`, `spec/CONSTITUTION.md` — architectural
  invariants every new calculator must respect.

**Existing library baseline (verified 2026-04-22):** 33
`*Result.h` classes in `src/`. `ChargeSource.h` already
establishes the typed-source-abstraction pattern. `GeometryChoice`
records runtime decisions. `KernelEvaluationFilter` holds
per-calculator filter sets. No `TrajectoryResult` or
autocorrelation / covariance / time-series machinery yet — all
trajectory-level work below is genuinely new territory. NICS /
virtual-node / probe-point concepts likewise absent. Architectural
extensions below build on the established patterns, not on
greenfield.

---

## A. New calculators / `ConformationResult` additions

### A.1 From CSA-as-a-field (L6-L12)

**A1. `CSAPrincipalAxisResult`** (new, conformation-level)
Diagonalise the per-atom σ tensor (from `OrcaShieldingResult` or
calibrated kernel sum), store principal components (σ₁₁, σ₂₂,
σ₃₃) and principal-axis-system orientations per atom. Mat3 is
already present; this just stores the diagonalisation product
alongside for direct comparison to L7 Yao-Bax 2010 per-stratum
magnitudes and L12 Loth-Pelupessy-Bodenhausen 2005 ubiquitin
per-bond principal components. Rollup W3 CCR derivation needs this.

**A2. `AmideTensorGeometryResult`** (new, conformation-level)
Per-amide bond: σ tensor expressed in the local N-H frame
(Brender-Taylor-Ramamoorthy 2001 convention, L11). Concretises
the L11 orientation into a typed tensor field instead of scattered
angle metadata, and closes the 13.2 R1/R2/R1ρ handle: CCR
back-calculation wants σ-in-bond-frame, not σ-in-lab-frame.

### A.2 From magnetic biophysics (N1-N3, B18)

**A3. `BulkSusceptibilityAccumulator`** (new, conformation-level,
accumulator-style)
Protein-level χ anisotropy from per-atom McConnell +
RingSusceptibility + HBond sums (Babaei et al. 2017 N3
hierarchical subunit approach). External benches: per-peptide Δχ
vs Pauling 1979 N2 (−5.36 × 10⁻⁶ emu) and bulk protein anisotropy
vs STI MRI (N4 Li 2017). **Methodologically novel output** — our
DFT-calibrated per-atom kernels → a bulk observable a different
literature measures. No other group produces this bridge from
calibrated classical kernels to MRI-scale χ. Directly thesis-
claimable (rollup section 13 synergy 3).

### A.3 From magnetic aromaticity / NICS (B14-B15)

**A4. `NICSProbeEvaluator`** (new, conformation-level, or
generalized as `VirtualNodeKernel`)
Evaluate BS / HM / McConnell / Coulomb kernels at arbitrary
non-atom probe positions. Enables NICS(0) at ring centroids and
NICS(1)zz at ±1 Å offsets (Schleyer 1996 B14). Natural companion
to E12 VN-EGNN virtual-node philosophy. **Cross-cutting benefit:**
unifies the volumetric shielding-field computation the h5-reader
viewer currently does via closed-form `QtBiotSavartCalc` /
`QtHaighMallionCalc` — the library becomes the single source of
truth for kernel evaluation at any point in space, not just at
atoms. Pairs with calculator hint #5 (FNO volumetric) as training-
data generator.

### A.4 From paramagnetic NMR (M24-M26, N5-N7)

**A5. `PseudocontactShiftResult`** [P] — see
`PLANNED_CALCULATORS_2026-04-22.md` §2.

**A6. `ParamagneticRelaxationEnhancementResult`** (PRE)
(new, conformation-level)
R1 / R2 enhancement from a paramagnetic centre via 1/r⁶
geometry — Solomon-Bloembergen-Morgan formalism. Different
observable from PCS (A5), shared geometry (distance to
paramagnetic centre). Pairs with A5; both validate against
Tharayil et al. 2021 N7 datasets.

### A.5 From stat-mech foundations (M27-M31)

**A7. `GreenKuboSpectralDensityResult`** [P] — see
`PLANNED_CALCULATORS_2026-04-22.md` §1.

**A8. `MemoryKernelExtractionResult`** (new, trajectory-level)
Zwanzig projection-operator memory kernel extracted from the
per-atom σ(t) autocorrelation. Complements A7 J(ω) on the
non-Markovian side. M29 Zwanzig 1961 primary, M22 Zwanzig 2001
textbook, M21 Kou-Xie power-law-memory evidence all tethered
directly.

**A9. `ErgodicityMetricResult`** (new, trajectory-level,
diagnostic-style)
Thirumalai ergodic measure per atom + block-averaged SD per
`TrajectoryResult` field. Formalises H7-H13 convergence
discipline as first-class output, per rollup 13.3's "shipped by
named default, not rule-coverage." One entry per trajectory;
provides honesty-caveat machinery for every other trajectory-
level field.

### A.6 From CCR and solution-state amide NMR (M10-M13, L12)

**A10. `CCRRateResult`** (new, trajectory-level)
Compute auto- and cross-correlated relaxation rates (DD/DD and
DD/CSA) from per-atom σ tensor × per-frame bond vector × Larmor
frequency. Consumes A7 J(ω) as intermediate. Directly produces the
observable compared to L12 Loth 2005 ubiquitin 64-bond dataset
(W3 handle). Thesis-grade validation: our pipeline's top-end
NMR observable.

### A.7 From cross-cutting validation (J-section, H-section,
L-section together)

**A11. `BenchmarkBackCalculationResult`** (new, per-bench
conformation-level)
Per named validation bench, compute our pipeline's prediction and
store alongside the published experimental value for direct
residual comparison. One attach per bench; queryable by bench
name. Initial roster: Loth-2005-ubiquitin-CCR, Yao-Bax-2010-15N-
CSA-split, Babaei-2017-bulk-χ, Tharayil-2021-PCS, Wylie-2006-
δ22-CO-HN. Turns the PHYSICS_FOUNDATIONS 0.9 validation-benches
table into concrete, consumed, tested pipeline output.

### A.8 From Stage 3 forward-looking (already captured)

**A12. Lie-group GP on SE(3)** [P] — Stage 3 Python, out of
nmr-extract scope, recorded for completeness.

**A13. FNO volumetric field** [P] — Stage 3 Python / viewer-
adjacent, out of nmr-extract scope, recorded for completeness.

---

## B. Variations on existing calculators

**B1. Unified `KernelSource` hierarchy** (architectural)
Formalise the "same K_ab, different source" pattern already
exploited by `MopacCoulombResult` and `MopacMcConnellResult`.
Parallels the existing `ChargeSource.h` hierarchy. Typed sources:
`RingSusceptibilitySource`, `BondSource`, `ChargedAtomSource`,
`LanthanideTagSource`, `NICSProbeSource`, `VirtualNodeSource`. A5
PCS, A6 PRE, A3 BulkSusceptibility, A4 NICS all slot in without
kernel-code duplication. Enables future sources (e.g., electronic
spin centres, unusual residues) without architectural churn.

**B2. Smooth cutoffs on all K_ab kernels**
`DispersionResult` already uses R_switch → R_cut tapering (per K2
Brooks CHARMM). Apply to McConnell, RingSusceptibility, HBond,
Coulomb. Plausibly reduces T2 residual artifacts at cutoff
boundaries and unifies cutoff discipline across the family.
MATHS_GOALS pillar 2 connection.

**B3. Multipole-expanded Coulomb**
Truncated dipole / quadrupole expansion of the Coulomb EFG for
far-field contributors; complements the near-field point-charge
sum. Reduces O(N²) cost at large systems and may regularise T2 at
far range. B17 Helgaker-Coriani 2012 covers the relevant
methodology.

**B4. Distributed ring current**
Replace the uniform-plane Johnson-Bovey circulation with per-
vertex current density. Addresses protonation-state-dependent
ring asymmetry (HID vs HIE) and heteroatom-localised current
(HisImidazole) that the uniform model blurs. M24 Bertini 2002
discusses anisotropic χ in heterocycles.

**B5. H-bond cooperativity**
Network-level correction over `HBondResult` output.
Cornilescu-Bax 1999 notes cooperativity matters — α-helix amide
H-bonds act differently from isolated H-bonds. Currently
HBondResult is pairwise. Post-pass using DSSP network topology.

**B6. Close MopacCoulomb / MopacMcConnell shielding_contribution
gap**
OBJECT_MODEL marks `mopac_coulomb_shielding_contribution` and
`mopac_mc_shielding_contribution` as TBD. Implementation-pass
rather than conceptually new calculator, but completes the
shielding-contribution coverage for Calculators 9 and 10 from
GEOMETRIC_KERNEL_CATALOGUE. **MOPAC cost context:** these fields
populate only on sparse selected frames (the 260 DFT pose set,
checkpoints), not per-frame across trajectories — see the
cost-model posture note above. Their primary use is the signal-
reconciliation diagnostic D5 against the per-frame-dense ff14SB
classical kernels, not as per-frame trajectory output.

**B7. C8 dispersion term**
Higher-order dispersion in `DispersionResult`; current is C6 only.
Literature in F4 Drautz ACE, D4 dispersion corrections.

---

## C. New input / output surfaces

### C.1 Input

**C1. ORCA Wiberg / Mayer bond orders from existing DFT runs**
`%output Print[P_Mayer] 1` in ORCA `.out` files. Wherever ORCA is
already run for shielding (the 260 calibration set, future μs-
harvester DFT), skip redundant MOPAC runs for bond orders. Faster
calibration path; MOPAC stays needed for the non-DFT pipeline.

**C2. `ExperimentalReferenceLoader`** (typed)
Read BMRB / RefDB shifts, L6-L12 CSA tables, N5-N7 PCS tables,
M10-M12 CCR datasets into a typed `ExperimentalReference` object.
Attached per ProteinBuildContext or as per-atom lookup. Enables
per-atom residual computation inside C++ rather than Python-only
post-pass. Pairs with A11 benchmark back-calculation.

**C3. Lanthanide-tag position + χ-tensor input**
For A5 `PseudocontactShiftResult`. N7 Tharayil et al. 2021
publishes per-mutant χ tensors for ubiquitin + GB1 — direct
consumable.

**C4. ORCA-computed magnetizability extraction**
ORCA can compute molecular magnetizability directly. Ingest for
comparison to A3 `BulkSusceptibilityAccumulator` output as a
same-DFT cross-check at the protein scale.

### C.2 Output

**C5. Irrep-tagged H5 metadata** (e3nn convention)
`0e`, `1o`, `2e` labels per tensor field (E2 Geiger-Smidt).
Downstream learners consume H5 as typed irrep tensors without
manual decoding. Zero-cost schema change; big convenience.

**C6. Machine-readable units + sign-convention metadata**
Per-H5-field unit and sign-convention attributes. Currently
documented in source comments only. JSON or H5 attributes;
eliminates "is this ppm or ppm × something" class of bug across
downstream Python / viewer consumers.

**C7. Validation-bench H5 slices**
Pre-computed groups `/benches/loth_2005/`, `/benches/yao_bax_2010/`,
etc., each carrying per-atom prediction + experimental + residual.
Pairs with A11. Bench-scoped slicing lets thesis-plot scripts open
only what they need.

**C8. `TrajectoryResult` serialization formalization**
Once-and-done schema per rollup 13.3 before implementation begins.
Named-convention-vs-schema decisions: layout of autocorrelation
lags, covariance tensor ordering, per-field block-SD fields (A9),
irrep tagging (C5). Touches every subsequent trajectory
calculator.

---

## D. Diagnostic analyses (not new calculators)

**D1. Per-secondary-structure CSA stratification** [P] — see
`PLANNED_CALCULATORS_2026-04-22.md` §3.

**D2. Per-calculator T2-residual map over the 260 DFT set**
First-class pre-computed output rather than notebook artifact.
MATHS_GOALS pillar 2 made physical: per-atom T2 residual
(DFT delta minus classical kernel prediction) binned by atom
type, secondary structure, distance-to-ring, etc. Thesis-grade
diagnostic.

**D3. Per-kernel-pair T2 correlation map**
`|cos|` between T2 of each calculator pair, per atom, over the
260 DFT set. Verifies or overturns the independence claims
already written in GEOMETRIC_KERNEL_CATALOGUE (e.g., PQ vs
McConnell ~0.38, MopacCoulomb vs Coulomb ~0.85). Auto-computed
from extraction; becomes a living thesis figure.

**D4. Ergodicity + block-SD per TrajectoryResult field** —
covered by A9.

**D5. MOPAC-vs-ff14SB kernel-signal reconciliation on the 260 DFT
poses**
Does the MOPAC-derived T2 signal (`MopacCoulombResult` +
`MopacMcConnellResult`) carry information beyond the ff14SB-driven
classical Coulomb / McConnell kernels on the *same* poses, and if
so how much? Concrete diagnostic: per-atom |cos| between (MopacCoulomb
T2 minus Coulomb T2) and (DFT delta T2 minus Coulomb T2), and the
analogous McConnell pair. Pairs with D3 (per-kernel-pair T2
correlation map) on the reconciliation side. **Calibration-hope
test (user 2026-04-22 pm):** if the sparse MOPAC signals and the
dense ff14SB signals agree on the pose set, per-frame ff14SB
kernels are sufficient for trajectory-level fields and the μs
harvester does not need per-frame MOPAC. If they diverge in a
structured way, the divergence itself is a finding (electron
polarisation carries trajectory-relevant information), and a
sparse-to-dense MOPAC correction becomes a design question to
raise. Either outcome is thesis-useful.

---

## E. Architectural / cross-cutting

**E1. `TrajectoryResult` as first-class object**
Already decided in rollup section 5. Flagging because every
trajectory-level item (A7-A10, C8, D4, the TrajectoryResult
fields per rollup 13.3) depends on this being designed once and
correctly.

**E2. Event menu hookable extractors** (rollup 13.7)
Ring flip detection (M17 Akke-Weininger), rotamer transitions
(Ramachandran / χ₁ bin-crossings), H-bond break/form events.
Threshold-parameterised menu under `/derived_events/`. Hookable in
both C++ and Python so a threshold sweep is one invocation, not
a regeneration pass.

**E3. Probe-point evaluation API**
Unified surface supporting A4 NICS, C7 volumetric training-data
generation, viewer butterfly rendering, user-defined lattice
sampling. Consolidates kernel evaluation into one non-atom-
aware code path.

**E4. `KernelSource` hierarchy formalisation** — covered by B1.

**E5. Calibration-ready vs raw-geometric separation schema-wide**
Per-calculator `shielding_contribution` is present for 8/10
classical calculators; B6 closes the MOPAC gap. Formalise the
schema: every calculator outputs (raw geometric kernel) and
(parameter-weighted shielding contribution in ppm), both stored,
both exported. No leakage of parameters into raw kernel fields.
CONSTITUTION §"Calculator Shielding Contribution Contract"
already states this; formalisation is enforcement.

**E6. Per-bench pipeline mode**
`nmr_extract --bench loth-2005` runs only what that bench needs
(subset of calculators + the required Experimental loader +
BenchmarkBackCalculationResult for that bench). Small CLI change,
big workflow effect for validation work. Pairs with C7.

---

## F. Scope boundaries — explicitly NOT in nmr-extract

- ML training pipeline (`learn/`)
- Visualization (`h5-reader/`, `ui/`)
- Equivariant-NN architecture or inference (`learn/c_equivariant/`)
- GP regression on SE(3) / Lie-group models (hypothetical
  `learn/c_gp/`)
- SE(3) diffusion sampling
- FNO training (Stage 3 Python territory)
- Python SDK structural changes (`python/` is a parallel project
  with its own spec)
- Wholesale BMRB / RefDB ingestion as a database (python/-side
  unless per-atom experimental values attach to typed objects per
  C2)

### Explicit exclusions from this pass (user, 2026-04-22 pm)

- **GFN-xTB as a `ChargeSource`** — considered and ruled out;
  problematic in our context per user guidance 2026-04-22. Do not
  re-propose without explicit user revisit.

---

## Thesis-angle highlights

Items most directly tied to a claimable thesis finding or
external validation bench:

| Idea | Thesis handle | Bibliography tether |
|---|---|---|
| A3 BulkSusceptibilityAccumulator | Calibrated per-atom kernels → MRI-scale bulk χ observable; independent validation | Pauling 1979 (N2), Babaei 2017 (N3), Li 2017 STI (N4) |
| A10 CCRRateResult + L12 bench | Per-frame σ × NH vector → ubiquitin CCR rates per 64 bonds | Loth 2005 (L12), Ferrage 2008 (M11), Tugarinov-Kay 2003 (M10) |
| D1 per-SS CSA stratification | α/β 15N CSA 11 ppm split reproduced or T2-residual exposed | Yao-Bax 2010 (L7), Wylie 2011 (L6) |
| D3 per-kernel-pair T2 correlation map | GEOMETRIC_KERNEL_CATALOGUE independence claims auto-validated | 260 DFT calibration set |
| A4 NICSProbeEvaluator | Connects our pipeline to aromaticity-criterion literature | Schleyer 1996 (B14), Gershoni-Poranne 2015 (B15) |
| A11 BenchmarkBackCalculationResult | Turns 0.9 validation-benches table into pipeline output | PHYSICS_FOUNDATIONS §0.9 |

---

## How to use this document

- This note is a consolidated starting place, not a prioritised
  roadmap. The user's next step per their 2026-04-22 framing is to
  organise the work against use cases and the thesis narrative —
  sequence, risk, dependency, overlap with existing spec — without
  pre-triage discarding anything from the list.
- `[P]` items (A5, A7, A12, A13, D1) are captured in more detail
  in `PLANNED_CALCULATORS_2026-04-22.md`; this doc references
  rather than repeats.
- Each item names its bibliography origin. If the physics or
  validation story behind an item is unclear, the bib entry is the
  starting reading, not this doc.
- Every item should be re-evaluated against the current state of
  `OBJECT_MODEL.md`, `spec/CONSTITUTION.md`, and
  `GEOMETRIC_KERNEL_CATALOGUE.md` before implementation. Ideas
  frozen 2026-04-22; the pipeline is not.
- Nothing here is committed design. The user decides what
  graduates to a full spec pass and in what order.

## Amendments

*(None yet. Append below, never rewrite above.)*
