# Meta Documentation Review

**DO NOT READ THIS IF YOU ARE A WORKING AGENT.** This directory
contains audit artifacts from a documentation review session
(2026-04-03). It is for the human author and her review partner.
If you are implementing code, building features, or fixing bugs,
you do not need anything in here. Read spec/INDEX.md instead.

---

## What this is

A documentation coherency review comparing every spec document
against the actual C++ code. The goal: make the spec documents true
again after intense implementation work outran the docs.

## Files

### Critiques (code-vs-doc audits)
- **CONSTITUTION_CRITIQUE.md** — section-by-section audit of
  CONSTITUTION.md. Each claim marked TRUE, STALE, ASPIRATIONAL,
  or WRONG with code evidence.
- **OBJECT_MODEL_CRITIQUE.md** — same treatment for OBJECT_MODEL.md.
  The named accessor table (18 convenience wrappers that don't exist)
  is the worst offender.
- **PATTERNS_CRITIQUE.md** — same treatment for PATTERNS.md. Core
  lessons all true. 7 missing patterns identified.

### Proposed edits (replacement text for Constitution)
- **CONSTITUTION_PROPOSED_EDITS.md** — first draft with reasoning
  notes explaining WHY each edit was made. Contains forward-looking
  language and TBDs. The working document.
- **CONSTITUTION_FINAL_EDITS.md** — clean replacement text in the
  constitution's own voice. True statements only. No TBDs. This is
  what gets applied to CONSTITUTION.md.

## Decisions made in this review

1. **Constitution speaks only truth.** No "planned" markers, no TBDs.
   When something gets built, it goes in. Until then, it stays out.
   Aspirational content lives in CALCULATOR_PARAMETER_API.md and
   MATHS_GOALS.md.

2. **Three foundational docs need collaborative editing.** Constitution,
   Object Model, Patterns. The critiques are the checklists. The final
   edits file has the replacement text for Constitution.

3. **Supporting docs updated directly.** INDEX.md, LAYER0_PLAN.md,
   MATHS_GOALS.md, APPLIED_MATHS_FIXES.md were updated in the review
   session because they are primarily consumed by agents and needed
   to be factually correct immediately.

4. **Geometric Kernel Catalogue untouched.** The maths is the maths.

5. **Gemini review triaged.** 5 of 15 files contain real validation
   data. The rest is flattery, duplication, or stale bug reports.
   Gold files in /mnt/extfast/gemini-nmr/:
   - baseline_consistency_analysis.md
   - kernel_independence_report.md
   - fused_ring_additivity_validation.md
   - mathematical_validation_report.md
   - numerical_stability_audit.md

## What's left to do

### Category 1: Mechanical fixes (can be done by any context with
   the critique files as checklists)
- Fix accessor names in Constitution and Object Model
- Fix test counts, source line counts
- Remove named accessor table from Object Model (or replace with
  "all access via Result<T>()")
- Add missing code descriptions to Object Model: RingBondedExclusion-
  Filter, CovalentTopology, MutationDeltaResult, ProtonationDetection-
  Result, shielding_contribution fields, PreloadedChargeSource,
  GromacsEnsembleLoader, WriteFeatures pattern

### Category 2: Judgement edits (replacement text ready)
- Apply CONSTITUTION_FINAL_EDITS.md to CONSTITUTION.md
- Apply equivalent edits to OBJECT_MODEL.md (same issues: hierarchy,
  protonation, copy semantics, named accessors, planned results)
- Update PATTERNS.md with missing patterns from PATTERNS_CRITIQUE.md

### Category 3: Code work (separate sessions)
- Formal charges through single pipeline path
- xTB/APBS batch validation
- Protonation builder
- ParameterCorrectionResult (e3nn integration)
