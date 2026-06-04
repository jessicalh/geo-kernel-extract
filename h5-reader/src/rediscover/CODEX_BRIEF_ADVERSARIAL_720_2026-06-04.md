# Codex brief — thoughtful review of the 720 static-ingest design (indicative maths + nothing real left on the floor)

> **Historical session brief — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Status: **the lead fired this.** Branch `h5-reader-pysr-spike` — NEVER touch git (lead owns ALL git). DOCS ONLY.

Be a **thoughtful scientist**, not a rubber stamp. You respect our engineering rules (the foreclosure, the
discipline) — of course — but that is your *secondary* concern. Your MAIN worry is **the maths**, and whether
this design risks **leaving real relationships on the floor**: real signal that is on the table *even with the
small given sample* (720 proteins, a sparse between-axis, indicative-not-proof) that a careful design would
surface and a careless one would quietly miss.

**Know this: we WANT the design exactly as it is now described** — the `StaticPoseConformation` cohort ingest,
the lean per-(protein, atom) C++-emitted substrate, the model-is-spine shape. This is **not** a redesign
review and **not** a "should we do something else" review. We are committed to this shape. The question is
narrow and important:

> **Is what is described OK maths *for what it is* — INDICATIVE (statistical position, not proof) — and free of
> basic maths errors? And does it leave any real relationship on the floor that the small sample could actually
> show us?**

## CONTEXT (read for intent + the constraints you're respecting)
- `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md` (original brief), `SPEC_720_STATIC_INGEST_2026-06-04.md` (the
  design under review), `REVIEW_opus_SPEC_720_2026-06-04.md` (prior review — go deeper, don't repeat).
- The reporting standard: **statistical position + determinability, indicative not proof** (`NOW.md`,
  memory `feedback_law_as_statistical_position`); no-simplification-bias (present the real complexity).
- Grounding for the maths: `ExtractionSupport.{h,cpp}` (T0/T2 target build + dia/para), `SphericalBasis.*`
  (frozen `get_C`), `PerAtomSubstrate.{h,cpp}` (the reductions/aggregation/exclusions), the field catalog
  (`../io/QtFieldCatalog*`) + `../../../python/nmr_extract/_catalog.py` (the full producer feature surface,
  ~430 trajectory columns), `GUIDANCE.md`, `PARTITION_FILTER_DESIGN.md`, `NODE_STORE_CONTRACT_2026-06-02.md`.

## THE REVIEW
1. **Indicative maths soundness.** For its actual purpose — INDICATIVE between-axis statistical position on a
   small, sparse-between sample — is the design's maths sound and *honest about being indicative* (not dressed
   up as proof)? Look at: the per-(protein, atom) aggregation; the C++ reductions; the **total-T2 target +
   dia/para split**; the T0/T2 basis (frozen `get_C`); the between-axis accumulator hooks (centering / null /
   strata / support-flagging); the small-N statistical shape. Is it valid *as indicative work*?
2. **Basic maths errors.** Hunt for any fundamental mistake — a wrong reduction, a basis/normalisation slip, a
   centering or null that would bias the between-axis estimate, a per-(protein, atom) aggregation that
   double-counts or drops, a DFT-target/frame handling that corrupts T2, a self/same-residue exclusion that
   severs a real relationship. Cite specifics with file/section refs.
3. **Real relationships left on the floor — your main worry.** Even with the small sample, what real,
   *capturable* relationships could this design MISS? Think hard about: over-aggregation that washes out a real
   effect; a stratum or grouping that hides a relationship the data could show; a producer feature *family*
   that carries a real relationship but the lean menu omits, so we never get to SEE it; a reduction that
   collapses a distinction (per-type, per-mechanism, orientation, dia/para, T1) that actually matters. NAME
   what would be left on the floor, and the **minimal, rule-respecting** addition that would surface it (a C++
   emitted scalar/sidecar or named reduction — never a Python model). *(This subsumes feature-completeness:
   enumerate the real static-relevant feature families vs the trajectory ~430 and flag any real-signal family
   the lean menu drops — but frame it as "would we miss a relationship," not as a disk-budget question.)*
4. **Rules respected (brief, secondary).** Confirm the foreclosure / no-discovery / extractor-SACRED /
   model-is-spine hold. One paragraph; this is not the focus.

## OUTPUT + DISCIPLINE
Write `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md`. DOCS ONLY: write exactly that one file; read the
codebase/substrate/catalog freely; change NO code, run NO 720/build/fit, run NO git. End by printing a summary:
**(a)** is the described maths OK-for-indicative + free of basic errors (yes/no + any error named); **(b)** what
real relationships, if any, would be left on the floor + the minimal rule-respecting fix for each; **(c)** the
one-paragraph rules check.
