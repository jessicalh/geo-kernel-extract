# Review items to assess — 2026-05-03

Triage doc, not a punch list.

Source: critical review by a general-purpose agent against
`spec/INDEX.md`, `PATTERNS.md`, `OBJECT_MODEL.md`, the recent
`spec/plan/` docs, and the topology / readback / frame-pull code as
landed at `HEAD = 776ca75`. Captured here so the substantive findings
don't get lost while we work through them. Several items require
conversation about what we are really doing — the agent flagged what
it could not resolve without conversational history.

OpenAI's parallel review may add items or refine the framing of
these. This doc is the staging surface for both.

The OBJECT_MODEL.md trim pass and the `spec/plan/` consolidation pass
are both queued; many items below will likely fold into one of those.
Items are grouped by the kind of work they prompt: decisions, deeper
questions, or prose-level cleanup.

---

## A. Decisions needed

The agent could not resolve these without conversational history. They
require us to say what we are actually doing.

### A1. `CalculatorContract.h` — adopt or remove?

> Resolved 2026-05-03 — see "Triage outcomes" section at end.

`src/CalculatorContract.h:1-43` exists. It is declared as the explicit
topology + charge contract by `spec/plan/legacy-amber-implementation-
brief-2026-04-29.md:25-32` and `openai-5.5-strong-architecture-
layout.md:328-335`. No calculator uses it: `grep CalculatorContract`
returns only the file itself.

This is a shape-without-rule. A future calculator author reads the
brief, sees the contract declared, looks for an example, finds none.
Two coherent answers exist:

- **"The contract is the new norm."** Then we need the first adopter,
  plus an exemplar pointer in INDEX.
- **"The contract is a future migration; for now bind concretely."**
  Then the file should probably go away until a real consumer
  arrives, or be marked clearly as forward-looking with no current
  consumers.

Decision needed.

### A2. `OBJECT_MODEL.md` "Pending full schema" / IUPAC-vocabulary block

`OBJECT_MODEL.md:2486-2498` lists pending H5 emission for typed
LegacyAmberTopology semantic fields: `Locant`, `BranchIndex`,
`DiastereotopicIndex`, `ProchiralStereo`, `PlanarStereo`,
`PseudoatomClass`, `PolarHKind`, `RingPosition`.

These names overlap heavily with the IUPAC-attempt vocabulary that
was reverted 2026-04-27. The block could be:

- The current authoritative pending list (LegacyAmberTopology will
  acquire these typed semantic fields).
- A leftover from the 2026-04-26 IUPAC episode that survived the
  revert and should be deleted during the OBJECT_MODEL.md trim.

Decision needed before the trim pass.

### A3. "Deferred work pointers" at `OBJECT_MODEL.md:2503-2521`

Five named items: LegacyAmberTopology + calculator-contract template,
full TrajectoryResult catalog, `TopologyTrajectoryResult`, byte-parity
validation pass, NVRTC rpath fix.

The agent flagged uncertainty about whether this list is current or
has been overtaken by the 2-protein dense validation strategy and
other shifts in the calibration plan since 2026-04-24. The NVRTC
item is resolved per `reference_nvrtc_rpath_fix` memory entry but
the doc doesn't reflect it.

Refresh needed during the trim pass.

---

## B. Items that look like prose fixes but expose a deeper question

### B1. Entry-point doc topology gap

The 2026-04-29 / 2026-05-02 entities — `LegacyAmberTopology`,
`LegacyAmberInvariants`, `GromacsToAmberReadbackBlock`,
`GromacsFramePullResult` — are well-explained in their headers and in
`spec/plan/`, but absent from `doc/ARCHITECTURE.md`, `PATTERNS.md`,
and `OBJECT_MODEL.md` as load-bearing entities. INDEX Tier 1 doesn't
lead to them. Rules currently live in code header comments + memory
entries; a future author following INDEX won't find them in the
canonical living docs.

**Surface fix.** Short `OBJECT_MODEL.md` section, "Invariant FF data
and trajectory-source capture," with one-paragraph contracts for the
three new entities + the readback-as-compiler-trace rule.

**Underlying question.** Are these entities settled enough to seat
in the canonical doc? The 2026-05-02 `LegacyAmberInvariants` shape
went through two iterations within one session (wrapper-with-
Attach/Has → plain const fields). The plain-fields shape is current
truth; before canonising, are we confident we're finished iterating,
or is one more shape change plausible?

### B2. Disulfide authority populator rule fragility

`LegacyAmberInvariants.has_disulfide_authority` (added 2026-05-03)
correctly disambiguates the three semantic states, and
`FullSystemReader.cpp:925` sets it correctly. The rule "every
populator of `disulfide_pairs` MUST also set
`has_disulfide_authority` whenever upstream chemistry exists" lives
only in the `LegacyAmberTopology.h` header comment.

Today there is one populator. When the AMBER-PRMTOP-direct path
lands (implied by `amber-implementation-plan-2026-04-29.md` Step 5,
`AmberPrmtopChargeSource`), it becomes a second populator. If the
second populator forgets, geometric SG-SG inference silently re-takes
authority on PRMTOP-direct loads.

**Surface fix.** Add the populator-must-set rule to the
`LegacyAmberInvariants` struct comment as an invariant that future
populators must respect.

**Underlying question.** The same shape will recur across other
`LegacyAmberInvariants` fields (`disulfide_pairs` is the first;
`mass`, `ff_atom_type_index`, etc. all carry "has_authority" potential
when multiple loaders exist). Do we want a per-field `has_X` flag
pattern, or a different shape for source-of-truth tracking?

### B3. PATTERNS §17 duplication-vs-chaining scope ambiguity

The "duplication preferred over chaining" rule appears in PATTERNS
§17 and the 2026-04-24 trajectory-scope addition. It's a
trajectory-scope rule. PATTERNS doesn't say so contrastively;
conformation-scope `ConformationResult` dependencies are routinely
chained (`OBJECT_MODEL.md:1377-1404` lists ~20 chained dependencies).

**Surface fix.** One-paragraph note in the addition's preamble or in
pattern 17: "trajectory-scope-specific; conformation-scope
dependencies remain chained per pattern 4."

**Underlying question.** As TrajectoryResults that depend on
ConformationResults proliferate (most of the PerFrameExtractionSet
catalog), will the boundary stay clean, or are there cases where the
right answer is genuinely "duplicate the conformation-scope work as
trajectory-scope state" rather than "chain through the
ConformationResult"? Current exemplars all chain.

---

## C. Prose-level cleanup

Small, mostly drift. Likely to fold into the `OBJECT_MODEL.md` trim
pass and / or the `spec/plan/` consolidation pass.

- **C1.** PATTERNS §13 "three shapes" vs the 2026-04-24 addition's
  "two shapes." `OBJECT_MODEL.md:1509` still says "Three coexisting
  field shapes (see PATTERNS.md §13)"; the §13 reference now points
  at a definition that says two. Add "(superseded by 2026-04-24
  addition; see end of this file)" to the OBJECT_MODEL.md headline.

- **C2.** `md.tpr` vs `production.tpr` drift. `OBJECT_MODEL.md:1469`
  and `OBJECT_MODEL.md:2585` (in the addition) both say `md.tpr`;
  current code uses `production.tpr` (`TrajectoryProtein.cpp:42`).
  Sibling to the queued `nmr_extract.cpp` fix.

- **C3.** `GromacsFramePullResult` rule not in `OBJECT_MODEL.md`.
  Header explains the contract; the canonical doc doesn't carry the
  "per-frame trajectory-source data goes here, not as direct
  `ProteinConformation` fields" rule. Folds into B1.

- **C4.** `FullSystemReader.cpp:716` — one-line comment at the
  `canonical_three = entry.canonical_three;` line: "string dies
  here; never lands on `Residue` or any persistent object." Anchors
  the readback-as-compiler-trace rule at the actual leak point a
  future author would be working in.

- **C5.** `AtomEvent::emitter` vs `kind`. `OBJECT_MODEL.md`'s
  2026-04-24 addition explains the distinction; the only code
  example is the two-kinds case (`BsAnomalousAtomMarker`).
  Single-kind emitters need a sentence: "single-kind emitters set
  `emitter == kind`."

- **C6.** `INDEX.md` "Where to start, by task" matrix is missing an
  "AMBER trajectory loader / readback feature" row pointing at
  `spec/plan/gromacs-to-amber-readback-block-design-2026-05-02.md`,
  `tests/test_amber_streaming.cpp`, `tests/test_amber_trajectory.cpp`.

- **C7.** PATTERNS.md "String traversals for atom or structure
  identity" anti-pattern section already has the "PDB LOADING
  BOUNDARY" pattern documented elsewhere; would benefit from a
  back-reference naming the legitimate construction-boundary
  consumers (`ResolveProtonationStates`, `DetectAromaticRings`,
  `CacheResidueBackboneIndices`) so a future author isn't tempted
  to extend the exception. "Consolidate" rather than "add."

---

## Connection to queued work

- **OBJECT_MODEL.md trim pass** (target ~30K tokens): items A2, A3,
  B1, C1, C2, C3, C5 will all be touched. Decisions on A2 and A3
  should land before the trim runs.
- **`spec/plan/` consolidation pass**: items A1 and A2 should be
  resolved before consolidation, since their answers affect what
  the consolidated plan says. C6 folds into the INDEX update that
  the consolidation will produce.
- **INDEX.md updates queued for 2026-05-02 work**: C6 folds in. The
  task matrix should also be checked for any other gaps a future
  author would hit.

---

## What this doc is not

A punch list. Several items — particularly A1, A2, A3 and the
underlying questions in B1 and B2 — need conversation about what we
are actually doing. The right next step on each of those is to
discuss the question, not act on the surface fix.

When OpenAI's parallel review lands, items it surfaces should be
appended in the same A / B / C structure (or a new section if the
shape doesn't fit), not interleaved with Claude's. Keeping the
provenance visible helps with triangulation.

---

## Triage outcomes — 2026-05-03 conversation

Decisions made and new cleanup goals surfaced while discussing the
agent's findings. Sequenced after A/B/C so the original
agent-surfaced framing stays legible.

### T1. A1 resolved — `CalculatorContract` as typedef-only convention

CHARMM is dead-letter — not "quarantined for some legacy invocation"
but "no live tool path produces a non-AMBER substrate." The current
prep+run pipeline funnels GROMACS-prepared trajectories
(AMBER ff14SB, with CHARMM-style decisions corrected back to AMBER
form) into a single canonical substrate type. That is the
unification win behind the "AMBER-only" decision.

With one topology type in the live tree, the typed enforcement value
of `CalculatorContract`'s `static_assert` collapses — it's checking
a tautology. Decision: `CalculatorContract` becomes a pure
typedef-as-documentation. Each calculator declares
`using Contract = CalculatorContract<LegacyAmberTopology, ChargeTableT>`
at its class header, where `ChargeTableT` is `ForceFieldChargeTable`
for charge-consuming calculators (Coulomb, ApbsField) or
`NoChargeTable` for purely geometric calculators (BiotSavart,
HaighMallion, McConnell, RingSusceptibility, PiQuadrupole,
Dispersion). The discipline is the substrate-declaration; the type
system does no enforcement work.

Implied cleanup work:

- Shrink `src/CalculatorContract.h` to the bare struct. Remove the
  `static_assert` and the `RequiredTopology<>` /
  `RequiredChargeTable<>` accessor templates — they protect against
  a multi-topology world that won't exist.
- Add `using Contract = ...` typedef to each existing calculator
  (~22 files). One-line addition per calculator.
- Document the typedef-as-discipline rule in PATTERNS.md so a future
  calculator author knows to declare it.
- Add an exemplar pointer in INDEX naming the first calculator to
  adopt the typedef as the template for new work.

Consequence over time: `ConformationAtom::partial_charge` and
`pb_radius` projections become deprecation-eligible as Contract
adoption grows. Calculators with `ChargeTable = ForceFieldChargeTable`
read directly from the table (which lives on `Protein`); the
ConformationAtom projection is for non-Contract code only. Once
every calculator declares a Contract, the projections are dead by
construction. Worth being explicit about: ConformationAtom
semantically splits into "computed fields from CRs" (live forever)
and "projections from invariant substrate" (deprecation-eligible).

Adjacent cleanup unlocked by CHARMM-dead framing:

- `src/xtc_reader.h` — superseded by libgromacs-direct TRR reading.
- `src/GromacsEnsembleLoader.{h,cpp}` — older CHARMM-XTC ensemble
  loader.
- `src/ChargeSource.cpp::GmxTprChargeSource::LoadCharges` subprocess
  path — replaced by libgromacs-direct.

Currently framed as "quarantined-legacy, do not invoke from new
paths." With CHARMM dead-letter, these are dead branches from a
dead-letter scope — stronger candidates for actual deletion than
"quarantined" suggests.

One small wrinkle to remember: a `<LegacyAmberTopology,
ForceFieldChargeTable>`-declaring calculator running against a
`--pdb` path (no readback, no FF data) is not a typed mismatch —
it's the same type with sparse fill (FF-numerical fields empty per
`LegacyAmberTopology.h:114-122`). The empty-vector signal is
conventional, not type-encoded. Calculators gate on
`!Mass().empty()` etc. when they need FF data; that gate stays
conventional. Contract doesn't change it.

### T2. Project-wide test inventory

Post-CHARMM-retirement, the test landscape splits unevenly:

- **Live AMBER tests against the new prep+run path.**
  `tests/test_amber_streaming.cpp`, `tests/test_amber_trajectory.cpp`.
  Hit `tests/data/fleet_amber/` (1P9J_5801 + 1Z9B_6577, gitignored,
  ~6.2 GB).
- **CHARMM-era trajectories and the tests written for them.**
  `tests/data/fleet/` is the CHARMM/XTC fixture; the 5 currently-
  failing FleetLoader tests, `SmokeTest.NoDft`,
  `JobSpecE2E.FleetLibraryDirect`. Per CLAUDE.md (2026-04-29 block),
  these failures are "documented and inert" pending the AMBER
  fixture replacement — they're testing dead-letter scope.
- **`tests/data/fleet_test_large/`** (~295 MB, gitignored) — older
  MD trajectory fixtures of indeterminate current relevance.
- **`tests/data/sdk_geo_only/`** — SDK reader fixture; consumed by
  Python tests; status post-AMBER-cutover unclear.
- **Smoke test framework.** `tests/golden/smoke/`, `tests/golden/
  blessed/`. The 28 NPY binary-diff failures (deferred 2026-04-24)
  are against blessed baselines that may or may not still represent
  intended output. Bless-deferred decision is documented in memory
  (`project_smoke_test_bless_deferred_20260424`).

A bunch of trajectory-style tests claim to test the system but are
testing dead-letter scope, against fixtures that may or may not be
worth keeping. Unclear at a project level which tests are live,
which are inert, which could be re-derived against AMBER fixtures,
and which are testing things we no longer do.

Need: a project-wide test inventory that, for every test file in
`tests/`, answers:

- What does it actually test?
- Against what fixture (location, format, force-field family)?
- At what scope (smoke / unit / integration / acceptance)?
- Is the fixture live (AMBER, current prep path), dead-letter
  (CHARMM-era), or fixture-free?
- If dead-letter: delete, or re-derive against an AMBER fixture?

Output should make it clear which tests are testing live scope and
which are testing dead-letter scope. That gates several downstream
decisions:

- Which fixtures to delete from disk (CHARMM-era ones are large).
- Which tests to delete or re-derive.
- Whether the smoke-test framework needs a new AMBER-based smoke
  set (and a fresh bless, replacing the 28 deferred NPY diffs).
- Whether SDK reader tests against `tests/data/sdk_geo_only/`
  remain meaningful.

This unblocks the "small fixes / improved tests / docs cruft" work
the current session opened with — knowing what tests we actually
have and run is the precondition for the rest. It also closes the
loop on the 5 inert FleetLoader failures: rather than living in
documentation as "inert pending fixture replacement," each test gets
a real disposition (delete / re-derive / repurpose).

The inventory itself is the artifact. Producing it is a separate
piece of work; this entry just captures the goal so it does not get
lost.
