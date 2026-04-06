# Object Model Critique — Code vs Document

**Date**: 2026-04-03
**Purpose**: Audit of OBJECT_MODEL.md (1914 lines) against actual C++
code. Same categories as Constitution critique: TRUE, STALE,
ASPIRATIONAL, WRONG.

The object model is what agents code against. If it says "property X
exists on class Y" and it doesn't, agents waste time looking for it
or create shadow structures. If it describes convenience wrappers that
don't exist, agents try to call them and fail. Every claim must be
verifiable.

---

## Core Types (lines 15-48)

**TRUE.** Vec3, Mat3, SphericalTensor, FieldValue, CalculatorId all
match Types.h exactly.

---

## KernelEvaluationFilter (lines 50-84)

**TRUE.** KernelEvaluationContext, KernelEvaluationFilter ABC,
KernelFilterSet all match KernelEvaluationFilter.h exactly.

**MISSING from doc**: `source_ring_index` field on
KernelEvaluationContext (used by RingBondedExclusionFilter). Added
after the doc was written.

**MISSING from doc**: `RingBondedExclusionFilter` — the fourth
concrete filter, added during the calculator audit (2026-04-02). It's
a topological check (ring vertices + bonded neighbours) distinct from
the distance-based DipolarNearFieldFilter. Used by BS, HM, PQ, RS,
Disp. This is a significant omission — it's the filter that fixed the
HM fused ring artifact.

**MISSING from doc**: `conf` pointer on KernelEvaluationContext
(optional, for future physics-aware filters).

**Filter sets per calculator table**: STALE. Missing:
- BS: DipolarNearField + RingBondedExclusion
- HM: DipolarNearField + RingBondedExclusion
- PQ: DipolarNearField + RingBondedExclusion
- RS: DipolarNearField + RingBondedExclusion
- Disp: DipolarNearField + through-bond vertex exclusion (inline)
- HBond: SelfSource + DipolarNearField

---

## AtomRole (lines ~160-215)

**TRUE.** All 17 roles match Types.h AtomRole enum and
EnrichmentResult assignment logic.

---

## BondCategory (lines ~215-250)

**TRUE.** All 8 categories match.

---

## Bond (lines ~250-280)

**TRUE for topology.** Bond.h has order, category, atom_a, atom_b.

**STALE**: Doc says Bond has "length, direction, midpoint." These are
conformation-dependent and live on ProteinConformation (bond_lengths,
bond_directions, bond_midpoints vectors set by GeometryResult), not
on Bond.

---

## Residue (lines ~282-326)

**MOSTLY TRUE.** Residue.h has type, protonation_variant_index, seq,
chain_id, insertion_code, backbone indices (N, CA, C, O, H, HA, CB),
chi atom indices.

**STALE details**:
- `ResidueProtonation` value type described in doc does not exist as a
  separate struct. Protonation variant is a single int
  (protonation_variant_index) on Residue.
- Doc mentions `formal_charge` on Residue — actual: it's on the
  ProtonationState decision, not on Residue itself.

---

## ChargeSource hierarchy (lines ~328-360)

**TRUE.** ForceField enum, ToolContext enum, ChargeSource abstract
class — all match ChargeSource.h.

**MISSING from doc**: `PreloadedChargeSource` (wraps pre-extracted
charges from TPR via libgromacs). Added with fleet loader.

---

## Protein (lines ~363-420)

**MOSTLY TRUE.** Heap-only, non-copyable, owns conformations via
factory methods.

**STALE API names**: Doc lists `CrystalConformation()` — actual is
`CrystalConf()`. Doc lists `MDFrames()` — actual is
`MDFrameAt(i)`/`MDFrameCount()`. Doc lists `Predictions()` — actual is
`PredictionAt(i)`/`PredictionCount()`.

**MISSING from doc**: `Conformation()` — the primary conformation
accessor that all 198 test calls use. Returns ProteinConformation&
regardless of subtype.

### FinalizeConstruction (line 384)

**TRUE for the three-layer pipeline:**
1. CacheResidueBackboneIndices (symbolic)
2. DetectAromaticRings (symbolic)
3. CovalentTopology::Resolve (geometric)

**MISSING from doc**: CovalentTopology as a named concept. The doc
describes bond detection via OpenBabel but doesn't mention
CovalentTopology.cpp/h as the geometry→topology boundary object.

---

## Atom (lines ~423-460)

**TRUE.** Atom.h matches: element, pdb_atom_name, residue_index,
bond_indices, parent_atom_index. Factory via Atom::Create().

---

## ConformationAtom (lines ~462-770)

**MOSTLY TRUE.** ConformationAtom.h is the wide table with all fields.

**STALE fields in doc that don't exist in code:**
- `hybridisation` is listed as on ConformationAtom — it IS there,
  confirmed.
- Fields for `ParameterCorrectionResult` exist as PLACEHOLDERS
  (correction_predicted_delta, correction_residual_T0/T1/T2,
  correction_sigma_T0/T2) but the result class that would fill them
  does not exist.
- Fields for `PredictionResult` exist as PLACEHOLDERS (predicted_T0,
  predicted_T2, confidence, tier) but the result class does not exist.

**MISSING from doc (fields that exist in code but not in doc):**
- `per_type_pq_scalar_sum[8]`, `per_type_pq_T2_sum[8][5]` — per-type
  PiQuadrupole accumulation
- `per_type_disp_scalar_sum[8]`, `per_type_disp_T2_sum[8][5]` —
  per-type Dispersion accumulation
- `bs_shielding_contribution`, `hm_shielding_contribution`,
  `mc_shielding_contribution`, `coulomb_shielding_contribution`,
  `hbond_shielding_contribution`, `piquad_shielding_contribution`,
  `ringchi_shielding_contribution`, `disp_shielding_contribution` —
  per-calculator shielding contribution SphericalTensors
  (these are the fields WriteFeatures outputs to NPY)

The shielding_contribution fields are significant — they're the primary
output of each calculator and the basis for the T2 independence
analysis. The doc's ConformationAtom section was written before all 8
calculators were implemented.

---

## Ring (lines ~797-1015)

**TRUE.** Ring class hierarchy matches Ring.h exactly. All 8 types,
virtual properties, accumulated fields.

**MISSING from doc**: `LiteratureIntensity()` virtual method — returns
the published value separately from `Intensity()` (which may be
corrected in future).

---

## ProteinConformation (lines ~1048-1300)

### Hierarchy (line 1048)

**STALE.** Same issue as Constitution: doc shows intermediate classes
(ExperimentalConformation, ComputedConformation, NMRConformation,
MinimisedConformation) that don't exist. Actual flat hierarchy:
ProteinConformation → Crystal, Prediction, MDFrame, Derived.

### Core geometry (line 1073)

**STALE**: Lists `protonation` as a ProteinConformation property.
Protonation is per-Protein (determines atom list), not per-conformation.

### Named result accessors (line 1113)

**ENTIRELY STALE.** The 18 named convenience wrappers listed
(Dssp(), ApbsField(), BiotSavart(), etc.) DO NOT EXIST. All access
is via `Result<T>()` template. This was the deliberate design choice
— PATTERNS.md says the template pattern is the ONLY template in the
system, and named wrappers were never implemented.

**This is the most misleading section in the document.** An agent
reading this will try to call `conformation.BiotSavart()` and get a
compile error.

### Query patterns (line 1196)

**PARTIALLY IMPLEMENTED** — same assessment as Constitution critique.
KD-tree queries work. Combined role+spatial queries require manual
intersection.

### Copy policies (line 1265)

**STALE.** GeometryOnly/Full copy policies not implemented. Protein
is non-copyable. Same as Constitution.

---

## ConformationResult types (lines 1315-1914)

### Implemented and accurate:
- GeometryResult, DsspResult, ChargeAssignmentResult, EnrichmentResult,
  SpatialIndexResult, ApbsFieldResult, MolecularGraphResult,
  BiotSavartResult, HaighMallionResult, PiQuadrupoleResult,
  RingSusceptibilityResult, DispersionResult, McConnellResult,
  CoulombResult, HBondResult, OrcaShieldingResult, XtbChargeResult
- All 17 result types exist as .h/.cpp pairs with correct dependencies.

### Missing from code:
- **MutationDeltaResult** — EXISTS in code but MISSING from doc. It's
  a ConformationResult attached to the WT conformation after comparing
  WT and mutant. Has atom matching, delta computation, ring proximity
  analysis. ~470 lines of code with its own Compute() method.
- **ProtonationDetectionResult** — EXISTS in code but MISSING from doc.
  Detects protonation state from hydrogen atoms present.

### Aspirational (described but not implemented):
- **ParameterCorrectionResult** (lines 1856-1908) — 52 lines of
  detailed API spec describing e3nn integration, 93 parameters, 8
  calculator-specific correction structs, residual tracking,
  SubtractCalculatorContribution. NONE of this exists in code.
  ConformationAtom has placeholder fields but no result class.
- **FeatureExtractionResult** (line 1910) — does not exist. Features
  are output via distributed WriteFeatures() per ConformationResult.
- **PredictionResult** (line 1913) — does not exist.

---

## Enforcement: Framework Stores (line 1657)

**SAME as Constitution critique.** The typed store interface
(StoreRingContribution) with automatic logging does not exist.
Calculators write directly to ConformationAtom fields.

---

## Summary: What Needs Changing

### Accurate and critical (DO NOT TOUCH):
- Core Types (SphericalTensor, FieldValue, CalculatorId)
- KernelEvaluationFilter framework (update to add RingBondedExclusion)
- AtomRole and BondCategory taxonomies
- Atom identity model
- Ring class hierarchy
- All 17 implemented ConformationResult descriptions
- ConformationAtom field inventory (update to include shielding
  contribution fields and per-type accumulation arrays)

### Must fix:
- **Named result accessors table**: DELETE or rewrite as "All access
  is via Result<T>(). No convenience wrappers exist." This is the
  single most misleading section.
- **ProteinConformation hierarchy**: Replace with actual flat hierarchy.
- **Protonation on conformation**: Move to Protein.
- **Bond geometry on Bond**: Clarify that length/direction/midpoint
  are per-conformation (on ProteinConformation), not on Bond.
- **Copy policies**: Mark as not implemented, describe rebuild pattern.
- **Framework stores section**: Rewrite to match actual mechanism.
- **ParameterCorrectionResult**: Mark clearly as ASPIRATIONAL spec for
  future ML integration. Do not remove — it's the design for the next
  phase — but label honestly.
- **FeatureExtractionResult, PredictionResult**: Same treatment.

### Must add:
- RingBondedExclusionFilter + source_ring_index on context
- CovalentTopology as a named concept
- PreloadedChargeSource
- MutationDeltaResult and ProtonationDetectionResult
- Per-calculator shielding_contribution fields on ConformationAtom
- Per-type accumulation arrays (PQ, Disp) on ConformationAtom
- Protein::Conformation() accessor
- WriteFeatures pattern (distributed, per-ConformationResult)
- GromacsEnsembleLoader as a loading path
