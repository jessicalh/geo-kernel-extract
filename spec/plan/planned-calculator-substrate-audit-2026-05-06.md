# Planned-calculator substrate audit — 2026-05-06

## Purpose

Phase 1 of the substrate runtime integration is about to start. The
substrate (`src/SemanticEnums.h:781-809` — `AtomSemanticTable` — plus
the lookup-key `AtomMechanicalIdentity` at lines 848-863) has been
designed and committed in a generated table at
`src/generated/LegacyAmberSemanticTables.{h,cpp}`. What is **not** yet
in place is runtime exposure: `LegacyAmberTopology` (`src/LegacyAmberTopology.h`)
currently holds covalent topology + invariant FF-numerical fields, but
no per-atom semantic record. Phase 1 will add that exposure.

The risk this audit addresses: Phase 1 narrows the accessor design to
what *current* calculators consume, locking out *planned* calculators
that would benefit from richer substrate exposure. Per
`feedback_planned_calculators_stay_planned`, the planned-calculator
list is kept on principle. This audit characterises which planned
calculators have substrate touchpoints, which substrate fields are
"symbolic richness only", whether any substrate gaps exist, and which
shape of accessor best serves the planned set.

## Scope conventions

- **Planned calculators** in this audit: the union of
  - `spec/PLANNED_CALCULATORS_2026-04-22.md` (5 ideas)
  - `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md` (~40 ideas in Sections A-F)
  - `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` (8 ideas)
  - `spec/plan/comprehensive-calculator-inventory-2026-04-30.md` Sections 10-14 (improvements + new + orphans)
- **Substrate fields** in this audit: the 14 typed fields of
  `AtomSemanticTable` enumerated below (with their `SemanticEnums.h` line numbers).
- **Existing calculators**: a thin reference column in §2 only.
  Existing-calculator state is not re-audited here; per the task brief,
  use `spec/plan/mechanical-identity-model-and-audit-2026-05-05.md`
  for that.
- **No new calculator proposals.** The list is the input.
- **No code changes.** This is a research artifact.

## The 14 substrate fields, indexed for cross-reference

`AtomSemanticTable` (`src/SemanticEnums.h:781-809`) holds these fields.
Per `topology-substrate-implementation-plan-2026-05-05.md` lines 92-117,
six are mechanical-identity (used as the `LookupBy` key) and the rest
are chemistry-substrate. I list all 14 in lookup order, with the
accessor-name shape Phase 1 will likely expose:

| # | Field | Type / values | Lines |
|---|---|---|---|
| F1 | `element` | `Element` (H/C/N/O/S/Unknown) | `SemanticEnums.h:786`; via `Types.h` |
| F2 | `locant` | `Locant` (None/Alpha/.../Eta) | `SemanticEnums.h:119-128, 787` |
| F3 | `branch` | `BranchAddress {outer, inner}` | `SemanticEnums.h:152-162, 788` |
| F4 | `di_index` | `DiastereotopicIndex` (None/Position2/Position3) | `SemanticEnums.h:178-182, 789` |
| F5 | `backbone_role` | `BackboneRole` (None/N/CA/C/O/AmideH/AlphaH) | `SemanticEnums.h:89-99, 790` |
| F6 | `prochiral` | `ProchiralStereo` (NotProchiral/ProR/ProS/Unassigned) | `SemanticEnums.h:206-213, 794` |
| F7 | `planar_group` | `PlanarGroupKind` (10 values) | `SemanticEnums.h:231-298, 795` |
| F8 | `planar_stereo` | `PlanarStereo` (NA/E/Z/Unspecified) | `SemanticEnums.h:319-324, 796` |
| F9 | `pseudoatom` | `PseudoatomMembership {kind,locant,branch,super}` | `SemanticEnums.h:354-383, 797` |
| F10 | `polar_h` | `PolarHKind` (13 values incl. AmineNH, AromaticOH, etc.) | `SemanticEnums.h:406-479, 798` |
| F11 | `ring_position` | `RingPosition {primary, secondary}` of `RingMembership` | `SemanticEnums.h:496-632, 799` |
| F12 | `aromatic` | `bool` | `SemanticEnums.h:800` |
| F13 | `formal_charge` | `int8_t` | `SemanticEnums.h:801` |
| F14 | `is_exchangeable` | `bool` (derived from `polar_h != NotPolar`) | `SemanticEnums.h:802` |

Two further fields appear in the table but are only documented as
runtime-relevant: `equivalence_class : uint8_t` (line 803, RDKit
canonical-rank-within-residue), and the explicit note (lines 804-808)
that **Hybridisation** (`Types.h::Hybridisation`) and **`BondOrderMask`**
(`SemanticEnums.h:669-676`) are intentionally **not** placed in the
runtime record yet — they live in the generation log and are surfaced
"when a calculator requires" them. I treat both as substrate-adjacent
in §3 below: a planned calculator that would consume them is a Phase 1
substrate-shape question.

`RingPosition` is itself a struct of two `RingMembership` records
(`SemanticEnums.h:600-607`); each `RingMembership` carries
`{ring : RingSystemKind, position : RingPositionLabel, ring_size,
aromatic, planar, n_heteroatoms}`. When this audit talks about
"`ring_position` consumers", that means consumers of any subfield of
`RingMembership` reachable through `primary` or `secondary`.

---

## §1 Calculator → substrate-field map

For each planned calculator, the substrate fields it would consume +
why. Calculators are grouped by the doc they are originally proposed in;
duplicates across docs are noted.

### 1.1 PLANNED_CALCULATORS_2026-04-22.md (5 ideas)

#### P1.1 GreenKuboSpectralDensityResult (trajectory-level)
**Substrate-irrelevant.** This is a Fourier transform of σ(t)
autocorrelation per atom; Wiener–Khinchin gives the spectral density
J(ω) at named Larmor frequencies. The physics commitment is on the
autocorrelation substrate (a `TrajectoryResult` field), not on
chemistry typing. No field consumed; no field needed.

#### P1.2 PseudocontactShiftResult (conformation-level)
**Light substrate dependence.** PCS = Δχ-tensor at lanthanide tag
contracted with a 1/r³ dipolar-kernel. The kernel reuses
`RingSusceptibilityResult`'s K_ab (per the planned doc's §2 "Type"
clause). Substrate touchpoints:
- `F1 element` — for diagnostic stratification (PCS sensitivity by
  element) but **not** for the kernel itself; the kernel sees positions
  + Δχ tensor.
- `F5 backbone_role` — useful for thesis-grade per-atom-class PCS
  reporting (vs. published per-residue PCS tables, e.g. Tharayil 2021
  N7), not for the calculation.

The lanthanide tag's χ-tensor and position are *external* to the
substrate (loaded per conformation, per `PLANNED_CALCULATORS_2026-04-22.md`
lines 110-114). They do not extend the substrate.

#### P1.3 Per-secondary-structure CSA stratification (diagnostic)
**Substrate-critical.** The diagnostic explicitly stratifies by atom
class within secondary structure. From
`PLANNED_CALCULATORS_2026-04-22.md` lines 158-162: "L7 Yao-Bax 2010
reports ¹⁵N CSA magnitudes of −173 ± 7 ppm in α-helix vs −162 ± 6 ppm
in β-sheet". The "¹⁵N" is element + role: backbone N specifically.
Fields needed:
- `F1 element` — to filter ¹⁵N (or ¹³C, ¹H per nucleus).
- `F5 backbone_role` — to isolate `BackboneRole::Nitrogen` vs.
  side-chain nitrogens (AmineNH, GuanidiniumNH, etc.).
- `F11 ring_position` — for ring-atom strata (Phe Cδ vs. Cε vs. Cζ
  per the L6 Wylie 2011 full-tensor data).
- `F10 polar_h` — for proton-donor vs. methyl-H vs. aromatic-H strata
  in proton stratification cuts.
- `F7 planar_group` — for amide vs. carboxylate vs. guanidinium
  carbon strata (¹³C C=O versus C-zeta of Arg, etc.).
- `F12 aromatic` — boolean fallback when only "aromatic vs. not"
  matters, e.g. ring-current vs. non-ring contributions.

Joined to `DsspResult.secondary_structure` (per atom's residue), these
fields make `(element, role-class, secondary structure)` strata
machine-computable. This is the substrate's clearest direct payoff.

#### P1.4 Lie-group GP regression on SE(3) (Stage 3 model neighbour)
**Substrate-irrelevant for the kernel itself**; uses kernel features
already extracted. The GP's input space is SE(3) embeddings of ring/bond
frames (`PLANNED_CALCULATORS_2026-04-22.md` lines 213-216). Substrate
fields might enter as **stratification labels** for separate GPs per
atom-type subset (per the planned doc's "May be practical only on
per-atom-type subsets" risk note, line 226); fields F1+F5+F11+F12 then
become the strata, identical to P1.3.

This calculator lives in `learn/`, not the C++ extraction library
(per the planned doc's "Type" line 205). Substrate exposure to
Python is via the H5 emission boundary, which Session E handles.

#### P1.5 FNO for volumetric shielding field (Stage 3 / UI neighbour)
**Substrate-irrelevant.** Trains on classical kernel evaluations on a
3D grid. Inputs are Cartesian; output is volumetric. No per-atom
chemistry needed (per `PLANNED_CALCULATORS_2026-04-22.md` lines 248-261).
Lives in viewer/`learn/`, not the C++ library.

### 1.2 NMR_EXTRACT_DESIDERATA_2026-04-22.md (Sections A-F, ~40 ideas)

#### A.1 — CSA-as-a-field (L6-L12)

**A1. CSAPrincipalAxisResult** — Diagonalisation of σ tensor.
Substrate-irrelevant for the diagonalisation; substrate-driven for the
**reporting** (it would naturally stratify by F5 + F11 + F7 the same
way P1.3 does). Per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:74-82`.

**A2. AmideTensorGeometryResult** — σ tensor in local N-H frame
(Brender-Taylor-Ramamoorthy convention, L11). Substrate touchpoints:
- `F5 backbone_role == AmideHydrogen` to gate which atoms get the
  per-amide tensor stored. The amide N-H frame is a backbone primitive.
- `F10 polar_h == BackboneAmide` is the equivalent gate from the H side.
Per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:83-89`.

#### A.2 — Magnetic biophysics (N1-N3, B18)

**A3. BulkSusceptibilityAccumulator** — Σ over per-atom McConnell +
RingSusceptibility + HBond contributions. Substrate-light:
- `F12 aromatic` — for verification that the per-ring summation
  covers exactly the aromatic atoms.
- `F11 ring_position` — same purpose.
- `F1 element` — element-stratified Δχ contributions.
Most of the work is in the per-atom kernels already; this calculator
sums them. Per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:91-103`.

#### A.3 — Magnetic aromaticity / NICS (B14-B15)

**A4. NICSProbeEvaluator** (or generalised `VirtualNodeKernel`) —
Evaluate kernels at non-atom probe positions (NICS(0), NICS(1)zz).
**Substrate-light, ring-coupled.** The probe position is geometric
(ring centroid + offset), but the relevant chemistry classification is:
- `F11 ring_position.primary.ring` — to restrict the probe set to a
  per-ring centroid for `NICSProbeEvaluator` (Phe vs. Tyr vs. Trp 5/6
  ring vs. His imidazole — distinct centroids).
- `F12 aromatic` — to gate the calculation's applicability.

Per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:107-117`. The probe-kernel
evaluation itself doesn't need substrate; the **selection of probe
positions** does.

#### A.4 — Paramagnetic NMR (M24-M26, N5-N7)

**A5. PseudocontactShiftResult** — same as P1.2 above. Substrate
fields F1, F5 for stratification only.

**A6. ParamagneticRelaxationEnhancementResult (PRE)** — 1/r⁶
geometry from paramagnetic centre (Solomon-Bloembergen-Morgan).
Substrate-light:
- `F1 element` — proton vs. heteronucleus PRE rate.
- `F10 polar_h` — possibly for amide-H specific PRE reporting (Amide
  H exchange interacts with PRE).
Per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:124-130`.

#### A.5 — Stat-mech foundations (M27-M31)

**A7. GreenKuboSpectralDensityResult** — same as P1.1.
Substrate-irrelevant.

**A8. MemoryKernelExtractionResult** — Mori-Zwanzig kernel from σ(t)
autocorrelation. Substrate-irrelevant. Per
`NMR_EXTRACT_DESIDERATA_2026-04-22.md:137-142`.

**A9. ErgodicityMetricResult** — Thirumalai ergodic measure +
block-averaged SD. Substrate-irrelevant; trajectory-level diagnostic.
Per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:144-150`.

#### A.6 — CCR / amide NMR (M10-M13, L12)

**A10. CCRRateResult** — Auto- and cross-correlated relaxation rates
(DD/DD and DD/CSA) per atom. Substrate touchpoints:
- `F5 backbone_role == AlphaCarbon` — for backbone CSA-DD CCR per
  the Tugarinov-Kay formalism.
- `F10 polar_h == BackboneAmide` — for amide H-N CCR (the L12
  Loth 2005 ubiquitin 64-bond data).
- `F1 element` — for nucleus-specific CCR (¹⁵N, ¹³Cα, etc.).
Per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:154-162`.

#### A.7 — Cross-cutting validation

**A11. BenchmarkBackCalculationResult** — per-bench prediction +
experimental + residual. Substrate touchpoints **per bench**:
- Yao-Bax 2010 ¹⁵N CSA bench: needs F5 + F1 (backbone N).
- Loth 2005 ubiquitin CCR bench: F5 + F10 (amide H + backbone N).
- Tharayil 2021 PCS bench: F1 + F5 (element-stratified).
- Wylie 2006 δ22 CO/HN: F5 (backbone C and amide H).
- Babaei 2017 bulk-χ: F12 + F11 (aromatic atoms).
A11 is essentially a per-bench-customised stratification harness; its
substrate consumption is the union of A1-A10's. Per
`NMR_EXTRACT_DESIDERATA_2026-04-22.md:166-174`.

#### A.8 — Stage 3 forward (already noted)

A12 (Lie-group GP), A13 (FNO) — same as P1.4, P1.5.
Substrate-irrelevant for the kernels.

#### B — Variations on existing calculators

**B1. Unified `KernelSource` hierarchy** — architectural; no per-atom
fields. Per lines 188-196.

**B2. Smooth cutoffs on K_ab kernels** — modifies kernel, not chemistry.
Per lines 198-204.

**B3. Multipole-expanded Coulomb** — no substrate (geometry only).
Per lines 206-210.

**B4. Distributed ring current** — substrate-coupled:
- `F11 ring_position.primary.position` — distinguishes vertex types
  on the ring (per-vertex current density), e.g. Trp pyrrole NE1 vs.
  CD1, His variant atoms.
- `F11 ring_position.primary.ring` — for ring-type-specific
  per-vertex models (HID vs HIE asymmetry per the planned-doc text).
- `F11 RingMembership.n_heteroatoms` — to know whether a vertex is
  C, N, or other.
This is a meaningful substrate consumer. Per lines 212-217.

**B5. H-bond cooperativity** — network correction over `HBondResult`:
- `F10 polar_h == BackboneAmide` — to identify H-bond donors typed.
- `F10 polar_h == SidechainPrimaryAmide` for Asn/Gln side chains.
DSSP topology already provides the network; the substrate makes the
typing typed instead of string-keyed. Per lines 219-223.

**B6. Close MopacCoulomb / MopacMcConnell shielding_contribution gap** —
implementation gap, no chemistry. Per lines 225-236.

**B7. C8 dispersion term** — no substrate (geometry + element). Per
lines 238-241.

#### C — New input/output surfaces

**C1. ORCA Wiberg/Mayer bond orders** — input to MopacResult-class
calculators; substrate-irrelevant. Per lines 248-252.

**C2. ExperimentalReferenceLoader** — typed external bench loader.
**Substrate-coupled** because the per-atom lookup needs an identity
key. Either F2+F3+F4+F5 (mechanical identity) or BMRB/IUPAC name
projection (Session E). Per lines 254-259.

**C3. Lanthanide-tag position + χ-tensor input** — external. No
substrate. Per lines 261-263.

**C4. ORCA-computed magnetizability extraction** — external.
Substrate-irrelevant. Per lines 265-268.

**C5. Irrep-tagged H5 metadata** — output schema; no substrate. Per
lines 273-277.

**C6. Machine-readable units + sign-convention metadata** — schema;
no substrate. Per lines 279-282.

**C7. Validation-bench H5 slices** — output of A11; reuses A11's
substrate consumption. Per lines 284-288.

**C8. TrajectoryResult serialisation formalisation** — schema;
no substrate. Per lines 290-295.

#### D — Diagnostic analyses

**D1. Per-secondary-structure CSA stratification** — same as P1.3
above. Substrate-critical (F1, F5, F11, F10, F7, F12).

**D2. Per-calculator T2-residual map over 260 DFT set** — needs the
same stratification dimensions as D1: F1 + F5 + F11 + F7 to bin
residuals by atom-type. Per lines 304-309.

**D3. Per-kernel-pair T2 correlation map** — substrate-light
(stratification only); F1 + F5 for headline element/role bins. Per
lines 311-315.

**D4. Ergodicity + block-SD per TrajectoryResult field** — covered by
A9; substrate-irrelevant.

**D5. MOPAC-vs-ff14SB kernel-signal reconciliation** — substrate-light;
might use F1 + F5 to publish per-element-and-role correlation tables.
Per lines 321-337.

#### E — Architectural

**E1-E6** — architectural; none are per-atom calculators that consume
substrate fields directly. E2 (event-menu hookable extractors) is the
only nontrivial one: ring-flip and rotamer-transition events
**publish** events keyed by atom (per `NMR_EXTRACT_DESIDERATA_2026-04-22.md:349-354`).
The atom keys could be typed-identity rather than name strings; the
event detection itself is geometric.

#### F — Out of scope

Section F (line 381 onwards) lists ML/visualisation/Python items
explicitly outside `nmr-extract`. None affect substrate design.

### 1.3 PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md (8 ideas)

#### TS1. BlockAveragedConvergenceResult
**Substrate-irrelevant.** Block-averaged SEM per tensor component.
Stratification by F1 + F5 is useful for the *report* (per-element
convergence) but the calculation is not chemistry-typed. Per
`PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md:17-40`.

#### TS2. SigmaLipariSzaboResult
**Substrate-light.** Lipari-Szabo S² fit per atom per irrep. The
fit itself is substrate-irrelevant. Stratification F5 (backbone CA / N
S²) + F1 + F11 is the comparison-publication path. Per
`PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md:44-68`.

#### TS3. SigmaEssentialDynamicsResult
**Substrate-irrelevant** for the SVD. Stratification of mode-loadings
per element/role is the report; not the calculation. Per lines 72-95.

#### TS4. SigmaTuckerDecompositionResult
**Substrate-irrelevant** for the decomposition. Tucker mode loadings
along the kernel-family axis are pre-existing kernel labels; along the
atom axis they are atoms (any stratification is post-hoc, downstream).
Per lines 99-124.

#### TS5. CrossCorrelatedRelaxationResult
Same as A10 above. F5 + F10 + F1 for amide H-N CCR specifically.
Per lines 128-152.

#### TS6. SigmaMSMResult
**Substrate-irrelevant** for the MSM construction. Stratification of
slow modes by element/role is post-hoc reporting. Per lines 156-180.

#### TS7. SigmaMemoryKernelResult
Same as A8. Substrate-irrelevant. Per lines 184-207.

#### TS8. Aromatic-H geometry sanity-check
**Substrate-light.** Probe geometries for Case 1995 / Agarwal 1977
benchmarks are aromatic-H-context. Substrate-coupled in two ways:
- `F12 aromatic` + `F1 element == H` to identify "aromatic H" atoms
  in real proteins for validation against the synthetic probes.
- `F11 ring_position.primary.position` to compare e.g. Phe Hδ vs.
  Hε vs. Hζ shifts at known geometries.
Per lines 211-235.

### 1.4 comprehensive-calculator-inventory-2026-04-30.md (Sections 10-14)

#### Section 10 — Existing-calculator improvements

**10.1 Ring-normal stability fix** (BS + HM) — geometric. Substrate-irrelevant.
**10.2 Single-loop default audit** — calibration; substrate-irrelevant.
**10.3 H-bond geometry (θ over d)** — kernel improvement on existing
H-bond surface; could be substrate-aware via F10 (BackboneAmide vs.
SidechainPrimaryAmide → different θ statistics) but the calculation
is geometry. Per `comprehensive-calculator-inventory-2026-04-30.md:384-390`.
**10.4-10.6** — validation passes; no substrate.

#### Section 11 — New geometric kernels

**11.1 Ch3ShiftMethylKernelResult** — methyl-specific. **Substrate-critical:**
- `F9 pseudoatom` — `PseudoatomKind::M` identifies methyl-group Hs
  unambiguously. Sahakyan 2011's CH3Shift formula sums over methyl
  protons; the substrate's pseudoatom field IS the methyl-membership
  mask.
- `F1 element` — H + the parent C.
- `F2 locant` — to know which methyl per residue (Leu Mδ1/Mδ2,
  Val Mγ1/Mγ2 etc.).
- `F3 branch.outer` — to disambiguate Leu Mδ1 vs Mδ2 (1 vs 2).
Per lines 420-426.

**11.2 LarsenProcs15CorrectionsResult** — ProCS15 long-range corrections.
**Substrate-coupled:**
- `F5 backbone_role` — Δσ_BB^(i±1) is backbone-only.
- `F10 polar_h == BackboneAmide` — Δσ_HB H-bond donor identification.
- `F5 backbone_role == AlphaHydrogen` — Δσ_HαB.
- (sidechain) F2 + F11 + F12 — Δσ_sc identification of sidechain
  atom classes.
Per lines 428-434.

**11.3 PelloniDifferentialBiotSavartResult** — re-formulates Biot-Savart
as per-volume-element. **Ring-substrate-coupled:**
- `F11 ring_position.primary.ring` — selects the ring system to integrate.
- `F11 ring_position.primary.position` — tells the integrator
  near-vs-far per vertex (the Pelloni decomposition Pelloni 2004
  performs).
Per lines 436-442.

**11.4 MoynaRingCurrentComparisonResult** — comparison of HM/JB/PD on
ring atoms. F12 + F11 to filter eligible atoms; otherwise substrate-light.
Per lines 444-450.

**11.5 ChargeDifferentialKernelResult** — Δkernel = kernel(frame) -
kernel(reference). Substrate-irrelevant per atom (it's geometric);
stratification F1 + F5 is the report. Per lines 452-458.

**11.6 BackboneNHUnitVectorResult** — per-residue per-frame N-H unit
vector. **Substrate-critical:**
- `F5 backbone_role == Nitrogen` and `BackboneRole::AmideHydrogen` —
  identifies the N and H per residue.
- `F10 polar_h == BackboneAmide` — equivalent gate from H side.
- (negation) Pro is excluded since Pro has no `AmideHydrogen`; the
  substrate's per-atom encoding gives this exclusion for free.
Per lines 460-466.

**11.7 HydrationShellPolarAtomResult** — 2-closest-water selector per
polar atom. **Substrate-critical:**
- `F10 polar_h != NotPolar` — gates which atoms get water-tracked.
  This is the canonical polar-atom mask.
- `F14 is_exchangeable` — derived; same answer.
- `F7 planar_group == Carboxylate` — for stratifying carboxylate
  oxygens.
Per lines 468-474.

#### Section 12 — Trajectory aggregators

**12.1 BlockAveragedConvergenceResult** — same as TS1. Substrate-light.
**12.2 SigmaLipariSzaboResult** — same as TS2.
**12.3 SigmaEssentialDynamicsResult** — same as TS3.
**12.4 SigmaTuckerDecompositionResult** — same as TS4.
**12.5 GreenKuboSpectralDensityResult** — same as P1.1 / A7.
**12.6 KernelTimeSeriesMomentsResult** — moments of kernel time-series
per atom. Substrate-irrelevant for moments; F1 + F5 for stratified report.
Per lines 520-526.
**12.7 DihedralAutocorrelationResult** — backbone φ/ψ + side-chain χ_k.
**Substrate-coupled:**
- `F5 backbone_role` — to identify φ/ψ atom triples.
- `F2 locant` + `F3 branch` — to identify χ_k atom-quadruples per
  side chain (existing AminoAcidType chi-tables become projection
  layers).
Per lines 528-534.
**12.8 DihedralDistributionKDEResult** — same identification needs as
12.7. Substrate-coupled F5 + F2 + F3.
Per lines 536-542.
**12.9 HBondGeometryDistributionResult** — same as B5 above. F10 +
DSSP topology. Per lines 544-550.

#### Section 13 — External-tool wrappers

**13.1 LegolasShiftPredictorResult** — DNN; uses Element + position.
Substrate-irrelevant for the network forward pass; F1 + F5 for the
emit/compare layer.
Per lines 556-562.

**13.2 UCBShift2Result** — same shape as 13.1. Substrate-irrelevant.
Per lines 564-570.

#### Section 14 — Schema-without-producer ORPHANS

**14.1 RingExponentialSumResult** — Σ_rings exp(-r/scale)·g(ring,atom).
**Ring-substrate-coupled:**
- `F11 ring_position` — to identify ring sets (the rings being summed
  are detected via substrate-driven ring construction per Q3 of the
  pivot synthesis).
- F12 fallback for the binary "is in ring" gate.
Per lines 576-582.

**14.2 RingDistanceStatsResult** — mean ring-center distance + nearest
ring-atom distance. Same substrate touchpoints as 14.1. Per lines 584-590.

### 1.5 Summary table — per-calculator substrate-touchpoint

This is the consolidated `calculator → fields` index, sorted by
substrate intensity. **Substrate-critical** = the calculator's identity
gate is one or more substrate fields. **Substrate-coupled** =
stratification or selection but the kernel is field-free.
**Substrate-irrelevant** = no substrate touch.

| Calculator | Doc | Touchpoint class | Fields |
|---|---|---|---|
| 11.1 Ch3ShiftMethylKernel | inv §11.1 | **Critical** | F9, F1, F2, F3 |
| 11.6 BackboneNHUnitVector | inv §11.6 | **Critical** | F5 (∨ F10) |
| 11.7 HydrationShellPolarAtom | inv §11.7 | **Critical** | F10 (∨ F14) |
| P1.3 / D1 Per-SS CSA stratification | planned §3, des D1 | **Critical** | F1, F5, F7, F10, F11, F12 |
| A2 AmideTensorGeometry | des §A.1 | **Critical** | F5 (∨ F10) |
| 11.2 LarsenProcs15Corrections | inv §11.2 | **Critical** | F5, F10, F2, F11, F12 |
| 11.3 PelloniDifferentialBiotSavart | inv §11.3 | **Critical** | F11 (ring + position) |
| B4 Distributed ring current | des §B | **Critical** | F11 (per-vertex), F12 |
| A4 NICSProbeEvaluator | des §A.3 | **Coupled** | F11 (ring), F12 |
| A10 / TS5 CCRRate | des §A.6, ts §5 | **Coupled** | F5, F10, F1 |
| A11 BenchmarkBackCalculation | des §A.7 | **Coupled** | union of A1-A10 |
| B5 H-bond cooperativity | des §B | **Coupled** | F10 |
| 14.1 RingExponentialSum | inv §14.1 | **Coupled** | F11, F12 |
| 14.2 RingDistanceStats | inv §14.2 | **Coupled** | F11, F12 |
| 12.7 DihedralAutocorrelation | inv §12.7 | **Coupled** | F5, F2, F3 |
| 12.8 DihedralDistributionKDE | inv §12.8 | **Coupled** | F5, F2, F3 |
| TS8 Aromatic-H sanity-check | ts §8 | **Coupled** | F12, F1, F11 |
| C2 ExperimentalReferenceLoader | des §C.1 | **Coupled** | F2, F3, F4, F5 (lookup-key) |
| 11.4 MoynaRingCurrentComparison | inv §11.4 | **Light** | F12, F11 |
| A1 CSAPrincipalAxis | des §A.1 | **Light** (report only) | F1, F5, F11 |
| A3 BulkSusceptibility | des §A.2 | **Light** (verify) | F12, F11, F1 |
| A5/P1.2 PseudocontactShift | des §A.4, planned §2 | **Light** (report) | F1, F5 |
| A6 ParamagneticRelaxationEnhancement | des §A.4 | **Light** (report) | F1, F10 |
| TS1 BlockAveragedConvergence | ts §1 | **Light** (report) | F1, F5 |
| TS2 SigmaLipariSzabo | ts §2 | **Light** (report) | F5, F1, F11 |
| 12.5 GreenKuboSpectralDensity (P1.1, A7) | planned §1, ts §1 | **Irrelevant** | — |
| TS3 SigmaEssentialDynamics | ts §3 | **Irrelevant** | — |
| TS4 SigmaTuckerDecomposition | ts §4 | **Irrelevant** | — |
| TS6 SigmaMSM | ts §6 | **Irrelevant** | — |
| A8/TS7 MemoryKernelExtraction | des §A.5, ts §7 | **Irrelevant** | — |
| A9 ErgodicityMetric | des §A.5 | **Irrelevant** | — |
| 10.1-10.6 Improvements | inv §10 | **Irrelevant** | — |
| 11.5 ChargeDifferentialKernel | inv §11.5 | **Irrelevant** (per-atom) | — |
| 12.6 KernelTimeSeriesMoments | inv §12.6 | **Irrelevant** | — |
| 13.1-13.2 External-tool wrappers | inv §13 | **Irrelevant** (NN inputs) | — |
| P1.4 Lie-group GP | planned §4 | **Irrelevant** (Stage 3) | — |
| P1.5 / A13 FNO volumetric | planned §5 | **Irrelevant** (UI/Stage 3) | — |
| B1, B2, B3, B6, B7 | des §B | **Irrelevant** | — |
| C1, C3-C8 | des §C | **Irrelevant** | — |
| D2-D5 | des §D | **Light** (report) | F1, F5 |
| E1-E6 | des §E | **Irrelevant** (architectural) | — |

**Counts:** 8 substrate-critical, 9 substrate-coupled, 8 substrate-light,
~18 substrate-irrelevant. The substrate is most useful for new
**chemistry-typed kernels** (methyl, amide, ring-vertex, polar-atom)
and for **stratification harnesses** (per-SS CSA, benchmark
back-calculation). Time-series + statistical-decomposition calculators
mostly bypass substrate entirely.

---

## §2 Field → consumers map (reverse)

For each substrate field, the planned calculators that would consume
it, plus a brief existing-calculator presence note (drawn from the
chemistry-question audits 2026-05-05 and the comprehensive inventory).

### F1 `element` (`Element`)

**Planned consumers:** P1.3/D1, A1, A3, A6, A10, A11, B5, 11.1, 11.7,
12.7, 12.8, TS1-TS3, TS5, TS8, 13.1-13.2 (most as stratification, some
as kernel selection — proton vs heavy-nucleus paths).

**Existing consumers:** Element is the most-consumed property on
`Atom` already. Every kernel branches on `Element::H` vs heavy. F1 is
already the most-used field in the project; substrate exposure is
redundant with `Atom::element`. **The substrate's `F1` is structurally
identical to the existing `Atom::element`.** Its substrate role is as
part of the `AtomMechanicalIdentity` lookup key, not as a calculator
surface. Calculators continue to read `atom.element` directly.

### F2 `locant` (`Locant`)

**Planned consumers:** 11.1 (Ch3Shift methyl identification), 11.2
(Larsen sidechain corrections), 12.7 / 12.8 (chi-angle atom selection),
C2 (lookup key).

**Existing consumers:** none yet. Methyl, χ-angle, and side-chain
atom enumeration are currently string-keyed via
`AminoAcidType::rings[].atom_names` / `chi_angles[].atoms` (`const
char*` literals).

### F3 `branch` (`BranchAddress`)

**Planned consumers:** 11.1 (Leu Mδ1/Mδ2, Val Mγ1/Mγ2 disambiguation),
12.7 / 12.8 (χ-angle disambiguation), C2 (lookup key).

**Existing consumers:** none yet. Branch is a distinguishing axis only
exercised by methyl-pair and Arg side-chain Hs; existing calculators
don't differentiate at this level.

### F4 `di_index` (`DiastereotopicIndex`)

**Planned consumers:** C2 (lookup key); none of the planned calculators
in §1 explicitly stratify by Position2 vs. Position3.

**Existing consumers:** none.

**Note:** This field is **lookup-key-essential** but **calculator-light**.
A calculator that wanted prochiral asymmetry (e.g. a methylene shielding
asymmetry diagnostic) would consume it; none of the current planned
calculators commit to that. The `ProchiralStereo` field (F6) carries
the chemistry; F4 carries the IUPAC numeric label.

### F5 `backbone_role` (`BackboneRole`)

**Planned consumers:** P1.3/D1, A1, A2, A10/TS5, A11, B5, 11.2, 11.6,
12.7, 12.8, TS1, TS2, D2, D3, D5 — the **most-consumed substrate field**
across the planned set. Backbone-vs-sidechain stratification is in
nearly every backbone-NMR-observable kernel.

**Existing consumers:** Backbone caches on `Residue` (`res.CA`,
`res.N`, etc., per OBJECT_MODEL.md and chemistry-question-2). Currently
populated by string equality `name == "CA"` etc. Phase 1 will repopulate
these from `BackboneRole` typed equality (per the pivot synthesis §13.6).
Once that lands, every existing backbone-aware calculator implicitly
consumes F5 through the cache.

### F6 `prochiral` (`ProchiralStereo`)

**Planned consumers:** none of the planned calculators in §1 explicitly
stratify by ProR vs ProS. A future "methylene asymmetry diagnostic"
(not currently in the planned set) would consume it.

**Existing consumers:** none.

**Symbolic-richness-only candidate.** F6 is a clean CIP encoding with
high chemistry value (Markley 1998 Figure 1 worked example, Gly HA2 = ProR,
HA3 = ProS); it deserves the symbolic richness defence on chemistry-
literature grounds — every methylene-asymmetry calculation in the
literature uses CIP labels — even though no planned calculator names
it explicitly.

### F7 `planar_group` (`PlanarGroupKind`)

**Planned consumers:** P1.3/D1, A11 (per-bench amide / carboxylate
strata), 11.7 (carboxylate oxygen subclass for water tracking).

**Existing consumers:** none. Currently amide / carboxylate / aromatic
identification is via residue-name string match (e.g. ASP/GLU oxygens)
in calculator code — though most calculators don't actually distinguish
at this level today.

**Higher-yield-than-it-looks.** The 10-value taxonomy
(PeptideAmide / SidechainAmide / Guanidinium / Imidazole / Aromatic6Ring
/ Aromatic5Ring / Carboxylate / AromaticHydroxyl / AromaticOxide) maps
1-to-1 with chemistry literature on planar-group shielding; A11
benchmark back-calculation will consume it implicitly any time a bench
references an amide-N CSA or a carboxylate ¹³C shift.

### F8 `planar_stereo` (`PlanarStereo`)

**Planned consumers:** none in the planned calculator set explicitly
consume E/Z. C2 (`ExperimentalReferenceLoader`) needs it indirectly:
BMRB shift tables sometimes label HD21 vs HD22 separately, and the
BMRB convention is encoded by F8 (per
`topology-encoding-dependencies-2026-05-05.md` §C.2).

**Existing consumers:** none.

**Symbolic-richness with literature standing.** Same defence as F6:
Markley §2.1.4 + BMRB atom_nom.tbl establish E/Z labelling; the
encoding is correct (per the §C.2 corrections). No planned calculator
names it; a future Asn/Gln HD21/HD22 differential-shift diagnostic
would. **Symbolic-richness-only candidate.**

### F9 `pseudoatom` (`PseudoatomMembership`)

**Planned consumers:** 11.1 (Ch3Shift — methyl Pseudoatom::M is the
canonical methyl mask). 11.4 (Moyna ring comparison — pseudoatom R
for the ring grouping). C2 (lookup-key when an experiment reports a
pseudoatom). The 11.1 use is **the** strong consumer.

**Existing consumers:** none.

### F10 `polar_h` (`PolarHKind`)

**Planned consumers:** P1.3/D1, A2, A6, A10/TS5, A11, B5, 11.2, 11.6
(via amide-H gate), 11.7 (polar-atom mask), TS5, 12.9 (HBondGeometry).
13-value taxonomy (BackboneAmide / SidechainPrimaryAmide / IndoleNH /
AmmoniumNH / GuanidiniumNH / ImidazoleNH / CarboxylOH /
HydroxylOH_Aliphatic / HydroxylOH_Aromatic / ThiolSH / AmineNH /
OtherPolarH) is **second-most-consumed** after F5.

**Existing consumers:** none yet — current H-bond and polar-atom
identification is via element + bond-graph + sometimes residue-name
checks. Phase 1 makes the typed gate trivial.

**SHIFTX2 grounding** (Han et al., per `SemanticEnums.h:399`): "trains
separate models for HD21/HD22/HE21/HE22/HH11/HH12/HH21/HH22, confirming
the chemistry distinction matters at the shift level". The substrate
PolarHKind taxonomy maps directly to SHIFTX2's separate-model groups.
Anyone running `BenchmarkBackCalculationResult` against a SHIFTX2-style
output stratifies on this field.

### F11 `ring_position` (`RingPosition`)

**Planned consumers:** P1.3/D1, A4, A11, B4, 11.3, 11.4, 14.1, 14.2,
TS8 — the **third-most-consumed substrate field**. Every ring-aware
kernel and every per-vertex ring-current variant.

**Existing consumers:** Ring objects (5 ring-current calculators —
BiotSavart, HaighMallion, RingSusceptibility, Dispersion, +1) consume
`Ring::atom_indices` — a typed surface (chemistry-question-3 §0).
The current populator (`Protein::DetectAromaticRings`) is
string-keyed and Phase 1 (per the pivot synthesis §13.4) replaces it
with `ConstructRingsFromSubstrate`, reading F11 to group atoms by
`(residue index, RingSystemKind)`. After Phase 1, every existing ring
calculator implicitly consumes F11 via the populator.

### F12 `aromatic`

**Planned consumers:** P1.3/D1, A3, A4, B4, 11.4, 14.1, 14.2, TS8 —
broad usage as a coarse aromatic-vs-not gate. Lower information than
F11 (which carries ring-system + position) but cheap and intuitive.

**Existing consumers:** Ring objects gate aromaticity through
`AminoAcidType::is_aromatic` + `Ring::Intensity()`; F12 is the per-atom
restatement.

### F13 `formal_charge`

**Planned consumers:** A11 (per-bench charged-residue strata), per
the C.1 dependency-doc note ("substrate is the topology label, not the
live electron distribution") any calculator using formal_charge needs
to know it's Lewis-localised. The **CoulombResult / MopacCoulombResult**
kernels currently use AMBER-table charges; formal_charge is the
substrate's complementary integer-Lewis encoding.

**Existing consumers:** none of the geometric kernels consume integer
formal_charge directly (charges come from `ChargeAssignmentResult`
floats). F13 is mostly substrate-side bookkeeping. A future
`ChargeDeltaKernel` (no current planned-calculator entry) that wanted
"Lewis charge - AMBER charge" delta diagnostics would consume it.

### F14 `is_exchangeable` (derived from F10)

**Planned consumers:** 11.7 (HydrationShellPolarAtom) explicitly. All
F10 consumers can read F14 as the binary fallback.

**Existing consumers:** none.

### Field-utilisation summary

| Field | Planned-critical consumers | Planned-coupled+light | Existing consumers (after Phase 1) |
|---|---|---|---|
| F1 element | many (stratification) | many | EVERY KERNEL (already on Atom::element) |
| F2 locant | 11.1 | 11.2, 12.7-8, C2 | Phase 2 (replaces AtomRole-based string match) |
| F3 branch | 11.1 | 12.7-8, C2 | none today |
| F4 di_index | (none) | C2 | (none) |
| F5 backbone_role | 11.6, A2, P1.3, 11.2 | many | backbone caches (after Phase 1) |
| F6 prochiral | (none) | (none) | (none) |
| F7 planar_group | P1.3, 11.7 | A11 | (none) |
| F8 planar_stereo | (none) | C2 indirect | (none) |
| F9 pseudoatom | 11.1, 11.4 | (none) | (none) |
| F10 polar_h | 11.7, P1.3, A2, B5, 11.2 | many | Phase 2 (replaces H-bond name match) |
| F11 ring_position | P1.3, B4, 11.3 | A4, A11, 11.4, 14.1-2 | Phase 1 ring populator + 5 ring kernels |
| F12 aromatic | P1.3, B4 | A3, A4, 11.4, 14.1-2, TS8 | Ring populator (binary fallback) |
| F13 formal_charge | (none) | A11 indirect | (none) |
| F14 is_exchangeable | 11.7 | F10 fallback users | (none) |

**Fields no current or near-term planned calculator consumes (the
"symbolic-richness-only" set):** F4 `di_index`, F6 `prochiral`,
F8 `planar_stereo`, F13 `formal_charge`. All four have strong
chemistry-literature standing per the citations baked into
`SemanticEnums.h` (Markley 1998 Fig 1 for F4, F6; CCD
`pdbx_stereo_config` for F8; per `topology-encoding-dependencies-2026-05-05.md`
§C.1 Lewis-localised convention for F13). Per
`feedback_planned_calculators_stay_planned`, these stay; their defence
is symbolic richness + chemistry-literature standing.

---

## §3 Substrate gaps

Are there planned calculators that need information **not** in the
current 14-field substrate?

### Gap G1 — Hybridisation
**Planned consumer:** 11.3 (PelloniDifferentialBiotSavart) implicitly
needs sp2/sp3 classification per ring atom for the per-volume-element
integration. 11.1 (Ch3Shift) needs sp3 confirmation for methyl C.
Some Δσ_sc terms in 11.2 (Larsen) need it for sidechain C/N typing.

**Substrate status:** **Already encoded but not in the runtime record.**
Per `SemanticEnums.h:804-808`: "Hybridisation (from Types.h) and
BondOrderMask are not included in the runtime record at the moment;
they are available via the generation log if needed for analysis. The
runtime is intentionally minimal until a calculator requires a
specific field on the substrate."

**Proposal:** When 11.1 / 11.2 / 11.3 land, promote `Hybridisation`
into `AtomSemanticTable` as F15. The generator already produces it
(via RDKit `Atom::getHybridization` — `SemanticSource::RDKit_Hybridisation`,
`SemanticEnums.h:697`). One field, one accessor, no new enum.

### Gap G2 — Bond order to neighbour (BondOrderMask)
**Planned consumer:** 11.2 (Larsen Δσ_sc) for distinguishing aromatic
vs single-bond sidechain neighbours; B4 (distributed ring current)
for vertex-pair bond-order weighting; 14.1 (RingExponentialSum)
for distance-decay weighting where bond-order matters at near contact.

**Substrate status:** **Already designed but not in the runtime record.**
Same status as G1 — `BondOrderMask` (`SemanticEnums.h:669-676`) is
defined but excluded from `AtomSemanticTable` at lines 804-808.

**Proposal:** Same as G1 — promote when a consuming calculator lands.
Note that bond-order information **is** accessible via
`CovalentTopology::Bonds()` (current `LegacyAmberTopology::Bonds()`
accessor returns the bond list; `Bond::category` is typed), so the
need is more for atom-local convenience than for missing chemistry.
A calculator can already reach bond orders via `Bonds().BondsAt(atom_index)`.

### Gap G3 — Magnetic susceptibility tensor (paramagnetic χ)
**Planned consumer:** P1.2 / A5 (PseudocontactShiftResult) — Δχ tensor
anchored at lanthanide-tag position; A6 (PRE) — paramagnetic centre.

**Substrate status:** **Conformation-side, not substrate-side.** Per
the planned doc text (`PLANNED_CALCULATORS_2026-04-22.md` lines 110-114),
the χ tensor is loaded **per conformation**: "lanthanide_position: Vec3,
chi_tensor: Mat3 (from paramagpy fit)". The substrate is per-residue/-atom
chemistry; tag position is per-conformation. **Not a substrate gap.**

The matching SDK input surface is C3 (`Lanthanide-tag position +
χ-tensor input`). It belongs as a `ConformationResult` external-input
loader, not on `LegacyAmberTopology`.

### Gap G4 — Per-atom anisotropic thermal factor
**Planned consumer:** None of the planned calculators in §1 explicitly
consume B-factors. A "B-factor-aware diagnostic" was raised in the
task brief as an example; the comprehensive inventory + planned docs
do not contain such a calculator entry today.

**Substrate status:** Not on the substrate; not a per-atom-class
chemistry field. B-factors are crystal-structure-side data carried by
`Atom` (existing `b_factor` field per OBJECT_MODEL.md `## Atom`
section). **Not a substrate gap; existing `Atom` field.**

### Gap G5 — Per-atom solvent-accessibility class
**Planned consumer:** 11.7 (HydrationShellPolarAtom) needs polar-atom
classification (covered by F10). A future "solvent-accessibility-aware
NH chemical-shift correction" (e.g. a Wishart-style protocol) is **not**
in the planned set today.

**Substrate status:** SasaResult (existing) and HydrationShellResult
(existing) are conformation-side. Per-atom-static "is this atom
typically buried" is structure-prediction territory. **Not a substrate
gap; conformation-side or external.**

### Gap G6 — Per-residue protonation history (time-dependent pH)
**Planned consumer:** "Time-dependent protonation under constant-pH MD"
flagged as Phase 2 work in `topology-and-identity-pivot-synthesis-2026-05-05.md`
§14. Substrate is invariant by design.

**Substrate status:** Real gap, but **Phase 2 / out-of-scope for Phase 1**.
A `TitrationStateTrajectoryResult` companion would carry per-frame
protonation; the substrate stays invariant.

### Gap G7 — His tautomer "confidence"
**Planned consumer:** Diagnostic — "the 'confidence' of HIS → HIE
default isn't surfaced" per
`topology-and-identity-pivot-synthesis-2026-05-05.md` §14.

**Substrate status:** Information lives in `SemanticProvenance.confidence`
(per `SemanticEnums.h:716-723`) — **already in the substrate** at
provenance level. Per the implementation plan, provenance is **not
in the runtime record** (`topology-substrate-implementation-plan-2026-05-05.md`
lines 775-779). A planned calculator wanting to surface it needs the
provenance pipeline lifted from the generation log into runtime, or
needs the provenance recorded at H5-emission boundary.

**Proposal:** Same as G1/G2 — promote when consuming work lands.

### Gap G8 — Per-frame planar-group geometry (ω, ring-flip, ring-pucker)
**Planned consumer:** B4 distributed ring current (HID vs HIE
asymmetry needs per-frame planarity), 12.7/12.8 (chi-angle work).

**Substrate status:** **Already designed and explicitly committed as
the conformation-side companion**:
`topology-substrate-implementation-plan-2026-05-05.md` §"Architecture
in one paragraph" line 28: "Conformation-dependent geometry (peptide
ω, ring-flip state, sp2 pyramidalization, ring-pucker phase) lands as
a per-frame `PlanarGeometryResult` ConformationResult companion".
**Not a substrate gap; PlanarGeometryResult lands paired with substrate
in C+D (the same Phase 1 sequencing).**

### Gap G9 — Cross-residue topology (peptide-bond i to i+1)
**Planned consumer:** 11.2 LarsenProcs15Corrections explicitly: "Δσ_BB^(i±1)".
12.7-8 dihedral work. A2 amide tensor geometry needs the C(i-1) →
N(i) → CA(i) → C(i) frame.

**Substrate status:** Not in `AtomSemanticTable` (per-atom record). But
**already accessible** via `Residue` cross-residue navigation in
`Protein` (existing residue-list ordering + per-residue backbone
caches). **Not a substrate gap; Residue-side and orthogonal to
AtomSemanticTable.**

### Gap G10 — `equivalence_class` (already in `AtomSemanticTable.equivalence_class`)
The `equivalence_class : uint8_t` field at `SemanticEnums.h:803` is
populated from RDKit `CanonicalRankAtoms` (per `SemanticSource::
RDKit_CanonicalRank`, line 696). Not yet exposed as F15 in this audit
(I treated the 14 listed in the implementation plan).

**Planned consumer:** `BenchmarkBackCalculationResult` per-bench could
use it for "atoms NMR-equivalent within residue" stratification (e.g.
Tyr Hδ1+Hδ2 symmetry-related). 11.4 Moyna ring comparison potentially.

**Substrate status:** **Present in the substrate; needs accessor.** Per
the implementation plan line 117 it is part of the runtime record; the
14-field count in this audit excluded it because the original task
brief framed the 14 as the symbolic surface. A 15th-field accessor is
straightforward.

### Gap summary

| Gap | Field needed | Substrate status | Proposed extension |
|---|---|---|---|
| G1 | Hybridisation | designed, not runtime | Promote to F15 when 11.1/11.2/11.3 land |
| G2 | BondOrderMask | designed, not runtime | Promote when needed; bonds graph already accessible |
| G3 | χ-tensor at tag | per-conformation | NOT a substrate gap (PCS conformation input) |
| G4 | B-factor | on Atom | NOT a substrate gap |
| G5 | solvent class | conformation-side | NOT a substrate gap |
| G6 | titration-state(t) | Phase 2 | NOT Phase 1 |
| G7 | provenance.confidence | designed, not runtime | Promote at H5 emission |
| G8 | per-frame planarity | PlanarGeometryResult | LANDED as Phase 1 companion |
| G9 | cross-residue topo | on Residue | NOT a substrate gap |
| G10 | equivalence_class | already in record | Add accessor in Phase 1 |

**Real Phase 1 gaps:** none that block Phase 1's accessor design.
The substrate is well-shaped; the four "symbolic-richness-only" fields
are defensible per chemistry-literature; the deferred fields
(Hybridisation, BondOrderMask) have specified promotion paths when
consumers materialise.

---

## §4 Accessor-design recommendation

Given §1+§2+§3, what should `LegacyAmberTopology` expose for Phase 1?
The choices framed in the task brief:

- **(a)** Per-field accessors only (`PolarHKindAt(atom_idx)` etc.)
- **(b)** Whole-row accessor only (`SemanticAt(atom_idx) -> const AtomSemanticTable&`)
- **(c)** Both: whole-row primary + per-field shortcuts.

### Analysis against planned calculators

The planned-calculator set divides cleanly:

**Calculators that consume one or two fields naturally use per-field
accessors.** Examples:
- 11.6 BackboneNHUnitVector — needs `BackboneRoleAt(ai)` only.
- 11.7 HydrationShellPolarAtom — needs `PolarHKindAt(ai)` only.
- 11.1 Ch3Shift — needs `PseudoatomAt(ai)` (struct read) + `LocantAt(ai)`
  + `BranchAt(ai)`.

For these, a per-field accessor produces calculator code that reads
plainly:

```cpp
if (tp.LegacyAmber().PolarHKindAt(ai) != PolarHKind::NotPolar) { ... }
```

This is the form `topology-substrate-implementation-plan-2026-05-05.md`
line 84 already pictures: "when a calculator needs to know 'is this
atom a backbone amide H?', it reads
`tp.LegacyAmber().PolarHKindAt(atom_idx) == PolarHKind::BackboneAmide`."

**Calculators that consume multiple fields per atom on a hot path
benefit from whole-row.** Examples:
- P1.3/D1 per-SS CSA stratification — F1, F5, F7, F10, F11, F12 in
  one stratum-keying pass per atom. Per-field would be 6 method calls;
  whole-row is one reference.
- 11.2 LarsenProcs15Corrections — 5+ fields per atom in the long-range
  correction loop.
- A11 BenchmarkBackCalculation per bench — variable but typically multi-field.

For these, a whole-row reference is simpler and faster:

```cpp
const AtomSemanticTable& sem = tp.LegacyAmber().SemanticAt(ai);
auto stratum = std::tuple{a.element, sem.backbone_role,
                          sem.ring_position.primary.ring,
                          sem.polar_h, sem.planar_group, sem.aromatic};
```

**No planned calculator wants to receive an `AtomSemanticTable`
reference as a struct argument** that isn't bound to the calculator's
own per-atom loop. The struct is a record, not a conceptual unit
calculators pass around.

### Recommendation: option (c), whole-row primary + per-field shortcuts

Specifically:

1. **Whole-row primary** — `const AtomSemanticTable& SemanticAt(size_t atom_idx) const`.
   This is the canonical accessor. It is the right shape because:
   - The substrate IS a structured record per
     `feedback_no_attach_lifecycle_for_invariant_data` ("plain const
     fields populated at construction"). Returning the record by const
     reference is the pattern.
   - Multi-field calculators get one cache-line read.
   - Single-field calculators pay nothing extra (the struct is
     ~16 bytes; the field they want is one offset into it).
   - The accessor name is unambiguous on the call site
     (`SemanticAt(ai)` reads as "fetch this atom's substrate row").

2. **Per-field shortcuts** for the high-frequency typed gates that
   the calculator code reads as predicates:
   - `BackboneRoleAt(size_t ai) const -> BackboneRole`
   - `PolarHKindAt(size_t ai) const -> PolarHKind`
   - `PlanarGroupAt(size_t ai) const -> PlanarGroupKind`
   - `RingPositionAt(size_t ai) const -> const RingPosition&` (struct, returned by ref)
   - `IsAromaticAt(size_t ai) const -> bool`
   - `IsExchangeableAt(size_t ai) const -> bool`

   These are PATTERNS.md "objects answer questions about themselves"
   shortcuts; they read as predicates in calculator code:

   ```cpp
   if (tp.LegacyAmber().IsExchangeableAt(ai)) { /* polar-atom path */ }
   ```

   Per-field shortcuts for low-utilisation fields (F4 `di_index`, F6
   `prochiral`, F8 `planar_stereo`, F13 `formal_charge`) — if there
   is no near-term planned consumer, **do not introduce shortcut
   accessors yet**. They can be added when a consumer lands; whole-row
   `SemanticAt(ai).di_index` works in the meantime. Don't proliferate
   shortcuts ahead of consumers.

3. **Identity access — separate accessor**:
   - `const AtomMechanicalIdentity& IdentityAt(size_t ai) const`
   - This is for graph-navigation / equality testing per
     `mechanical-identity-model-and-audit-2026-05-05.md` §1.5. It
     stays distinct from substrate access; identity is the locator,
     substrate is the answer.

### Why not option (a) only (per-field)

- Multi-field calculators get N method calls per atom per loop
  iteration. P1.3 stratification at 6 fields × ~10000 atoms × hundreds
  of strata is real volume.
- The substrate **is** a record; insisting on field-accessor-only
  hides that fact and forces the proliferation of accessor names per
  field even when no calculator needs them.
- F11 `ring_position` is a struct of structs; `RingPositionAt` *must*
  return a reference anyway (or a tuple of 4+ fields, which is the
  same thing badly). The whole-row form already exists for this field.

### Why not option (b) only (whole-row)

- Calculator code reads better with named predicates (`IsAromaticAt`,
  `IsExchangeableAt`).
- The PATTERNS.md "objects answer questions about themselves" rule
  is best honoured by typed self-questions, not by a generic
  whole-row fetch.

### Concrete header sketch (advisory only — Phase 1 implementation owns the precise signature)

```cpp
class LegacyAmberTopology final : public ProteinTopology {
public:
    // ... existing accessors ...

    // ── Substrate access ────────────────────────────────────────────
    // Whole-row primary (per AtomSemanticTable in SemanticEnums.h:781).
    const AtomSemanticTable& SemanticAt(size_t atom_idx) const;

    // Identity (lookup-key + equality-test surface; per the audit).
    const AtomMechanicalIdentity& IdentityAt(size_t atom_idx) const;

    // Per-field shortcuts for high-frequency typed gates.
    BackboneRole       BackboneRoleAt(size_t ai) const;
    PolarHKind         PolarHKindAt(size_t ai) const;
    PlanarGroupKind    PlanarGroupAt(size_t ai) const;
    const RingPosition& RingPositionAt(size_t ai) const;
    bool               IsAromaticAt(size_t ai) const;
    bool               IsExchangeableAt(size_t ai) const;

private:
    std::vector<AtomSemanticTable>      semantic_;   // size = atom_count
    std::vector<AtomMechanicalIdentity> identity_;   // size = atom_count
    // ... existing private fields ...
};
```

The internal storage as `std::vector<AtomSemanticTable>` is implied
by Phase 1 §13.6 in the pivot synthesis (`SetAtomSemantic(std::move(...))`
at `topology-and-identity-pivot-synthesis-2026-05-05.md:896`).

---

## §5 Planned-calculator scope check — substrate-blocked vs. ready

A planned calculator is **substrate-blocked** if it cannot run until
substrate exposure lands at runtime. **Ready** = it can run today
(modulo non-substrate dependencies). This section identifies which
calculators move from blocked to ready when Phase 1 lands.

### Blocked → ready (Phase 1 unblocks)

These planned calculators **cannot start** until `SemanticAt` /
field accessors exist. Phase 1 directly unblocks them:

| Calculator | Blocked because |
|---|---|
| 11.6 BackboneNHUnitVector | Needs `BackboneRoleAt` (∨ `PolarHKindAt`) for the gate |
| 11.7 HydrationShellPolarAtom | Needs `PolarHKindAt` for polar-atom mask |
| 11.1 Ch3ShiftMethylKernel | Needs `PseudoatomAt` + `LocantAt` for methyl identification |
| A2 AmideTensorGeometry | Needs `BackboneRoleAt` for amide-pair selection |
| B4 Distributed ring current | Needs `RingPositionAt` per-vertex labels |
| 11.3 PelloniDifferentialBiotSavart | Needs `RingPositionAt` for near/far decomposition |
| P1.3/D1 Per-SS CSA stratification | Needs F1+F5+F7+F10+F11+F12 in stratification harness |
| A1 CSAPrincipalAxis (report) | Needs F5+F11 for stratified report (calc itself is ready) |
| 11.2 LarsenProcs15Corrections | Needs F5, F10, F2, F11 for term gating |
| 11.4 MoynaRingCurrentComparison | Needs F12, F11 (substrate-light but uses) |
| TS8 Aromatic-H sanity-check | Needs F12+F11+F1 for protein-side validation set |
| A11 BenchmarkBackCalculation | Per-bench substrate consumption |
| 14.1 RingExponentialSum (orphan) | Needs F11/F12 for ring identification |
| 14.2 RingDistanceStats (orphan) | Needs F11/F12 |
| C2 ExperimentalReferenceLoader | Needs F2+F3+F4+F5 (lookup-key) for shift-table joining |

That is **15 of the ~50 planned calculators** that Phase 1 directly
unblocks. The Phase 1 substrate exposure is the gate.

### Ready today (substrate-irrelevant or kernel-only)

The substrate-irrelevant calculators in §1.5's "Substrate-Irrelevant"
class can ship before Phase 1, modulo their own non-substrate
dependencies:

- TS1 BlockAveragedConvergence, TS3 SigmaEssentialDynamics, TS4
  SigmaTuckerDecomposition, TS6 SigmaMSM, TS7 SigmaMemoryKernel
- A7/P1.1 GreenKuboSpectralDensity, A8 MemoryKernelExtraction,
  A9 ErgodicityMetric
- 10.1-10.6 (existing-calculator improvements)
- 11.5 ChargeDifferentialKernel, 12.6 KernelTimeSeriesMoments
- 13.1 LegolasShiftPredictor, 13.2 UCBShift2 (NN-input wrappers)
- D2-D5 diagnostics (substrate-light report only)
- B1, B2, B3, B6, B7 architectural / improvement
- E1-E6 architectural

These do **not** wait for Phase 1, but most ALSO do not need it.

### Substrate-coupled with ready paths

A few in §1.5's "Substrate-coupled" class have an existing-calculator
fallback that can run today:

- A4 NICSProbeEvaluator — could ship using existing `Ring` objects
  (typed) for ring centroids; F11 is a refinement, not a gate.
- A10/TS5 CCRRate — existing backbone caches give F5 substitute;
  F10 substitute via element + bond graph. Can ship before Phase 1.
- B5 H-bond cooperativity — existing DSSP topology is the network
  primary; F10 is a typing improvement, not a gate.
- 12.7 DihedralAutocorrelation, 12.8 DihedralDistributionKDE —
  existing AminoAcidType chi-tables are the fallback; F2+F3+F5
  refinement.
- 12.9 HBondGeometryDistribution — same as B5.

These can start with the existing surface and migrate to substrate
later.

### Substrate-blocked summary

**Phase 1 unblocks ~15 planned calculators.** That is the largest
single-event shift in calculator readiness in the planned set. Phase 1
is the right place to land it; the dependency graph is well-shaped.

---

## §6 Phase 1 + Session E context

Per `topology-substrate-implementation-plan-2026-05-05.md` §Sequencing:

- **Session B**: generator + tables + provenance — LANDED.
- **Session C+D**: substrate fields on `LegacyAmberTopology` +
  PlanarGeometryResult per-frame companion — Phase 1 work.
- **Session E**: IUPAC + BMRB + cross-system projections — projection
  layer.

The user's distinction (per task brief "Phase 1 + Session E context"):
calculators consume **typed substrate**; projections are
**output-side**. Do any planned calculators interact with naming
projections directly?

### Inspection of planned calculators against Session E projections

Session E produces (per the implementation-plan §Sequencing):
- `IupacName(atom_idx)`, `BmrbName(atom_idx)`, `AmberName(atom_idx)`
- H5 emission `/atoms/legacy_amber/iupac_name`, `bmrb_name`

**Calculators that consume names directly (none should, but check):**

- **C2 ExperimentalReferenceLoader** — must join experimental shift
  tables (BMRB, RefDB) keyed by atom names. **This is the one
  legitimate name-projection consumer.** Per the audit's
  "Crystal Projection Rule" (per
  `mechanical-identity-model-and-audit-2026-05-05.md` §1.7 and
  `amber-implementation-plan-2026-04-29.md`): names are
  pure functions on typed substrate. C2's join is one such projection
  at the wire boundary; the calculator body sees typed-identity-keyed
  data. Session E exposes the projection function the loader needs.

- **A11 BenchmarkBackCalculationResult** — references published
  experiments by atom name (Loth 2005 by ubiquitin amide name, Wylie
  2006 by atom name, etc.). Same as C2: the wire join uses
  `BmrbName(ai)` from Session E; the calculator body operates on
  typed identity.

- **All other planned calculators** — operate on typed substrate
  fields, not on names. No interaction with Session E.

**Calculators that emit names directly:** none of the planned
calculators in §1 are themselves name emitters. Naming projection is
output-side per the user's framing; calculators sit upstream of it.

### Conclusion

Phase 1 unblocks ~15 calculators directly via substrate exposure.
Session E unblocks **two more** calculators (C2 and A11 to a
bench-by-bench degree) by providing the wire-format projection layer
they need for external-data joins. The two pieces of work are
complementary, not competing: most planned calculators want Phase 1
(substrate fields); few of them want Session E (projection names).

The substrate's calculator-side payoff is firmly in Phase 1's
court. Session E's payoff is in **input loaders** (C2) and in
**output emission** (the H5 columns the SDK reads). Phase 1's
accessor design does not need to anticipate Session E.

---

## Closing notes

- **The 14-field substrate is well-shaped for the planned calculator
  set.** No structural gap blocks Phase 1.
- **Three deferred fields** (Hybridisation, BondOrderMask,
  provenance.confidence) have specified promotion paths when
  consumers materialise; Phase 1 should not pre-emptively wire them.
- **Four "symbolic-richness-only" fields** (F4 di_index, F6
  prochiral, F8 planar_stereo, F13 formal_charge) have chemistry-
  literature defences; per `feedback_planned_calculators_stay_planned`
  they stay.
- **Whole-row primary + per-field shortcuts** is the right accessor
  shape per the planned-calculator distribution: ~40% of substrate-
  consumers want one or two predicates; ~30% want multi-field
  stratification at high volume; both are best served by the dual
  surface.
- **15 planned calculators are Phase-1-substrate-blocked.** That is
  the work Phase 1 unblocks.
- **Session E (naming projections) interacts with two planned
  calculators (C2, A11)** at the wire-boundary join, not at
  calculator surface. Session E and Phase 1 are independent.

## Provenance

- `src/SemanticEnums.h` (the 14 fields + identity tuple)
- `src/generated/LegacyAmberSemanticTables.h` (runtime API)
- `src/LegacyAmberTopology.h` (current accessor surface)
- `OBJECT_MODEL.md` (Atom + Residue + ConformationResult)
- `PATTERNS.md` §"naming boundary" + §"Objects answer questions"
- `spec/PLANNED_CALCULATORS_2026-04-22.md` (5 ideas)
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md` (~40 ideas)
- `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` (8 ideas)
- `spec/plan/comprehensive-calculator-inventory-2026-04-30.md` (Sections 10-14)
- `spec/plan/topology-substrate-implementation-plan-2026-05-05.md`
- `spec/plan/topology-encoding-dependencies-2026-05-05.md` (§C, §H)
- `spec/plan/topology-and-identity-pivot-synthesis-2026-05-05.md`
- `spec/plan/mechanical-identity-model-and-audit-2026-05-05.md`
- `spec/plan/chemistry-question-3-ring-topology-2026-05-05.md`
