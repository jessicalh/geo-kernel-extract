# Cold-Read Review — 2026-04-24

A fresh-reader evaluation of whether the docs in `spec/INDEX.md` give a new
contributor what they need to (a) add a `ConformationResult` subclass and
(b) add a `TrajectoryResult` subclass. Observation only.

---

## 1. Reading path taken

In order opened:

1. `spec/INDEX.md` — recommended reading order, the three tiers — **clear**.
2. `doc/ARCHITECTURE.md` (Tier 1 item 1) — system map: Protein/conformation
   model, provenance paths, calculator pipeline, filters + GeometryChoice,
   PATTERNS correspondence table, file manifest. **Clear.** The PATTERNS
   correspondence table (`doc/ARCHITECTURE.md:468-484`) is particularly
   useful for bridging architecture prose to actual header lines.
3. `PATTERNS.md` (Tier 1 item 2) — patterns 1-12 for conformation scope,
   13-18 twice for trajectory scope (earlier section + "2026-04-24
   library-reference" addition at `PATTERNS.md:1059`).
   **Clear but noisy — both sections are live and that was flagged to me
   up front, otherwise I would have wasted time deciding which to trust.**
4. `OBJECT_MODEL.md` (Tier 2 item 3) — read the Core Types, up through
   `ConformationResult (ABC)` at `:1304`, the earlier "Trajectory-scope
   entities" section at `:1377-~1683`, and the 2026-04-24 addition at
   `:2469-2993`. Jumped past the historical "copy-and-modify" sections as
   instructed at the top of the file. **Clear.**
5. `spec/CONSTITUTION.md` §§ "Per-calculator minimum output" (`:682-743`)
   and the anti-simplification rule. **Clear.**
6. `src/ConformationResult.h` — 47 lines, the full ABC. **Clear.**
7. `src/ProteinConformation.h:60-102` — `AttachResult`, `Result<T>`,
   `HasResult<T>`. **Clear.**
8. `src/ProteinConformation.cpp:22-59` — verify AttachResult actually does
   the singleton + dependency check the docs promise. **Clear.**
9. `src/McConnellResult.h` + `.cpp` — conformation-scope exemplar. Covers
   `Dependencies()`, static `Compute()` factory, per-atom field writes,
   KernelFilterSet assembly, GeometryChoiceBuilder recording, tracelessness
   projection, `WriteFeatures` NPY emission. **Clear.**
10. `src/OperationRunner.cpp:30-165` — how `Compute` results are attached in
    order. **Clear.**
11. `src/TrajectoryResult.h` — the trajectory ABC (102 lines, heavily
    commented). **Clear.**
12. `src/TrajectoryProtein.h` — full header including `AttachResult`,
    `DispatchCompute`, `FinalizeAllResults`, `AdoptDenseBuffer<T>`,
    `GetDenseBuffer<T>`, templated `Result<T>`. **Clear.**
13. `src/TrajectoryAtom.h` — short (72 lines), only BS Welford fields +
    events bag. **Clear.**
14. `src/Trajectory.h` — 8-phase `Run` docstring + `TrajectoryEnv`. **Clear.**
15. `src/BsWelfordTrajectoryResult.h` + `.cpp` — AV exemplar; full
    cross-result-read marker discipline visible. **Clear.**
16. `src/BsT0AutocorrelationTrajectoryResult.h` + `.cpp` — FO exemplar with
    `DenseBuffer<double>` transfer, biased ACF physics commitment. **Clear.**
17. `src/PositionsTimeSeriesTrajectoryResult.h` + `.cpp` — FO with no
    ConformationResult dependency; simplest FO exemplar. **Clear.**
18. `src/RunConfiguration.h` + `.cpp:120-157` — factory registration,
    `AddTrajectoryResultFactory`, `RequireConformationResult`. **Clear.**
19. `CMakeLists.txt:195-242` — source-list registration (single-target pattern
    from PATTERNS.md confirmed). **Clear.**
20. `tests/` listing — noticed **no `test_*_trajectory_result.cpp`** anywhere.
    **Confusing** — tests exist for every conformation-scope calculator
    (e.g. `test_mcconnell_result.cpp`) but nothing for BsWelford, BsT0Acf,
    BondLengthStats, Positions, ChiRotamer, BsShieldingTimeSeries, or
    BsAnomalousAtomMarker. No doc tells a new reader what the test
    contract for a `TrajectoryResult` is supposed to look like.

I did not read: Tier 3 physics docs (MATHS_GOALS, GEOMETRIC_KERNEL_CATALOGUE),
EXTRACTION_ORDER.md beyond the graph section, CALCULATOR_PARAMETER_API.md,
TEST_FRAMEWORK.md, `spec/pending_include_trajectory_scope_2026-04-22.md`
(correctly flagged as working notes), ENSEMBLE_MODEL, USE_CASES, any of the
session-dated state records, or ui/.

---

## 2. Could I add a new `ConformationResult`? **Yes.**

The ABC surface is small (`src/ConformationResult.h:23-45`), four overrides
max: `Name()`, `Dependencies()`, optional `WriteFeatures()`, and a static
`Compute()` factory. `McConnellResult.{h,cpp}` demonstrates everything the
docs promise:

- **ABC interface.** Documented in `OBJECT_MODEL.md:1304-1374`, verified in
  the header. No hidden virtual method.
- **Static `Compute()` factory.** `McConnellResult.cpp:104-306` — takes
  `ProteinConformation&`, returns `unique_ptr`. PATTERNS.md §5.
- **`Dependencies()` returns `vector<type_index>`.**
  `McConnellResult.cpp:18-23`. `ProteinConformation::AttachResult`
  (`:22-59`) actually enforces both singleton and dependency checks with a
  diagnostic message that lists what IS attached — promise made in docs,
  verified in code.
- **Field ownership on `ConformationAtom`.** Documented as a field table
  per result in `OBJECT_MODEL.md:658-788`. Each field has a `Source result`
  column; McConnell writes `mcconnell_co_sum`, `mc_shielding_contribution`,
  `T2_{CO,CN,backbone,sidechain,aromatic}_*`. PATTERNS.md §§3, 7, 11 spell
  out the one-writer-per-field singleton rule.
- **Geometric output + shielding contribution contract.** PATTERNS.md §6 +
  `CONSTITUTION.md:682-743` (per-calculator minimum output) + the T2
  completeness subsection (PATTERNS §11 and the dedicated "T2 Completeness"
  block at `PATTERNS.md:520-561`). Both Mat3 AND SphericalTensor stored;
  every tensor yields `shielding_contribution` in ppm alongside its
  natural-units form.
- **`WriteFeatures` serialisation.** `McConnellResult.cpp:362-400` + docs
  at `INDEX.md:56-69` (plus the Python SDK `_catalog.py` registration
  requirement).
- **Numerical stability.** `PATTERNS.md:898-1011` enumerates
  MIN_DISTANCE, KernelFilterSet-before-compute, tracelessness projection,
  near-field stability per kernel, unit chains in comments, NaN/Inf
  absent-not-faked. Filter-set construction is visible in
  `McConnellResult.cpp:121-128`.
- **GeometryChoice recording.** Required by `INDEX.md:42-54` and
  `PATTERNS.md:1012-1020`; the pattern is visible at
  `McConnellResult.cpp:172-182, 291-296`.
- **Test contract.** `tests/test_mcconnell_result.cpp` exists as a precedent
  but no doc enumerates what a conformation-scope test MUST assert.
  `PATTERNS.md "Lessons Learned §21"` (Steps 1-7 calculator-validation
  process) comes closest to spelling this out and is the best thing I
  read.

**Gaps I would still want to ask about:**
- Where exactly the *new field list* gets added to `ConformationAtom.h` —
  no explicit instruction, but `OBJECT_MODEL.md:604-609` hints.
- Whether `WriteFeatures` has a naming convention for NPY filenames (the
  existing ones are ad-hoc: `mc_shielding.npy`, `coulomb_efg.npy`, etc.).
  The docs say to register with Python `_catalog.py` — agreed — but give no
  example of that registration.
- Whether the new result needs to appear in `OperationRunner::Run`
  explicitly; reading `OperationRunner.cpp` answers yes, the runner is an
  explicit sequence rather than a topological sort of the dependency
  graph. `EXTRACTION_ORDER.md:30-40` prose says "pipeline resolves the
  order" but the code is a flat sequence with skip flags. Mild contradiction
  noted below.

---

## 3. Could I add a new `TrajectoryResult`? **Yes, with caveats.**

The ABC is 102 lines of heavily-commented header; six virtuals, only three
pure (`TrajectoryResult.h:27-99`). AV + FO exemplars cover all shapes I
plausibly need:

- **ABC interface.** `Name`, `Dependencies`, `Compute(conf, tp, traj,
  frame_idx, time_ps)`, `Finalize(tp, traj)`, optional `WriteFeatures` and
  `WriteH5Group`. The header comment names these by role and spells out
  the intended reads/writes for Compute explicitly (`:42-61`).
- **`Create(tp)` factory.** Documented in the 2026-04-24 addition
  (`OBJECT_MODEL.md:2625-2628`) + PATTERNS.md §14/15. The rationale for
  `Create` instead of `Compute` is that allocation is sized from
  `tp.AtomCount()`, `tp.BondCount()`, `tp.RingCount()`, and seed-precedes-
  attach means the factory sees a finalized Protein
  (`PATTERNS.md:1140-1155`). `BsWelfordTrajectoryResult::Create` at `.cpp:27-34`
  shows the clean shape.
- **`Compute` signature + contract.** Visible at every exemplar; the
  doc-stated contract (conf for this-frame inputs, tp for running outputs,
  traj for run-scope env/selections) is demonstrated concretely at
  `BsWelfordTrajectoryResult.cpp:44-67` with an explicit comment block
  naming the READS and WRITES.
- **`Finalize`.** `BsWelford.cpp:122-138` (AV — convert m2 to std).
  `BsT0Autocorrelation.cpp:90-168` (FO — allocate `DenseBuffer<double>`,
  walk per-atom history, compute biased ACF, transfer via
  `tp.AdoptDenseBuffer<double>`). `PositionsTimeSeries.cpp:57-82` (FO —
  flatten per-atom growing buffers into `DenseBuffer<Vec3>`, transfer).
- **AV vs FO lifecycle shapes.** Explicitly named in
  `OBJECT_MODEL.md:1524-1526` + `:2645-2653`, the 2026-04-24 addition
  reinforces the same. Exemplars map cleanly: AV writes
  `TrajectoryAtom` fields in place; FO appends to an internal buffer and
  materialises at Finalize.
- **Internal state vs TrajectoryAtom fields.** PATTERNS.md §13 is loud
  about this: accumulator implementation (Welford structs, history
  buffers, rolling windows) stays INSIDE the TR; finalized per-atom output
  lands on `TrajectoryAtom`. Exemplars comply; the pattern is visible.
- **`DenseBuffer<T>` ownership transfer at Finalize.**
  `src/TrajectoryProtein.h:131-145` — `AdoptDenseBuffer<T>(buffer, owner)`,
  `GetDenseBuffer<T>(owner)` keyed by owning TR's `type_index`.
  `OBJECT_MODEL.md:1659-1683` (earlier) and `:2869-2915` (addition) both
  enumerate the permissible T: `double`, `Vec3`, `SphericalTensor`, `Mat3`
  reserved.
- **H5 group emission.** Per-TR `WriteH5Group(tp, file)`; explicit
  named-component flattening required (`.x()/.y()/.z()`, `.T0/.T1[k]/.T2[k]`,
  `m(i,j)`) — no `reinterpret_cast`. Comment at
  `PositionsTimeSeries.cpp:94-104` is the canonical statement of this
  discipline.
- **Dependency declaration against TR or ConformationResult types.**
  `BsWelfordTrajectoryResult.cpp:22-24` returns
  `typeid(BiotSavartResult)` — a ConformationResult type — and the docs
  say Phase 4 validates this against `config.RequiredConformationResultTypes()`.
  `BsAnomalousAtomMarker` depends on another TR. Both shapes are documented
  and exemplified.
- **Factory registration in `RunConfiguration`.**
  `RunConfiguration.cpp:131-154` shows the six current factories registered
  into `PerFrameExtractionSet` as lambdas. `AddTrajectoryResultFactory`,
  `RequireConformationResult` are the API.

**Gaps I would still want to ask about:**
- **No TR tests exist in `tests/`.** Nothing in the docs tells me what a
  trajectory-result test should assert: frame-0 semantics, Finalize
  idempotency, AV field validity, H5 round-trip, NPY round-trip, etc.
  Conformation-scope has 47 `test_*.cpp` files; trajectory-scope has zero.
  `spec/TEST_FRAMEWORK.md` is referenced but I did not open it; it is not
  clear whether smoke tests cover trajectory output.
- **Cross-result read marker discipline** (PATTERNS.md §17) is documented
  and exemplified well (BsAnomalousAtomMarker), but the actual protocol of
  "comment in writer header, comment in reader header, marker at each read
  line" is easy to miss. A checklist would help.
- **Events bag H5 emission is not yet in the tree.**
  `OBJECT_MODEL.md:2582-2583` + `:2667` both say "per-atom events bag is
  not currently H5-emitted." A new scan-mode TR that wants its events to
  survive to the output file has no documented path. I can guess from
  the selection bag pattern but it is not stated.
- **No documented file where a new TR's factory gets registered.** I
  inferred it from reading `RunConfiguration.cpp`. A one-line "to add a
  new TR to production, append a factory in `PerFrameExtractionSet`" would
  close the loop.
- `ScanForDftPointSet` is flagged for removal in two places; OK.

---

## 4. The 2026-04-24 addition-form restatement: useful or confusing?

**Useful, but the coexistence is a hidden tax.**

Both sections are live in both `OBJECT_MODEL.md` (earlier
`:1377-~1683`, addition `:2469-2993`) and `PATTERNS.md` (earlier §§13-18
at `:197-316`, addition `:1059-1266`). The addition header explicitly
says "trust this as current truth; the originals are retained for design
history." Because the task brief warned me this was coming, I reached
for the addition when looking up concrete code-facing answers (method
signatures, exact field names, H5 schema) and reached for the earlier
sections when I wanted design motivation (why RecordBag exists at two
scopes, why attach-order-is-dispatch-order matters).

**When they disagree, the addition wins by stated policy.** I found no
outright factual contradiction between them on the trajectory-scope
shapes. Both say the same thing about:
- TrajectoryProtein.Seed precedes TR attach
- AV vs FO lifecycles
- DenseBuffer element types observed (double, Vec3, SphericalTensor, Mat3 reserved)
- RecordBag at two scopes
- Handler as pure reader, Run as orchestrator

Minor drift I noticed:
- PATTERNS.md §13 (earlier) lists THREE coexisting shapes on TrajectoryAtom
  (typed accumulator fields, typed struct vectors for known-shape
  per-source data, event bag). The addition at
  `PATTERNS.md:1088-1107` and `OBJECT_MODEL.md:2574-2606` lists only TWO
  (event bag + finalized rollup fields). The third shape ("typed struct
  vectors like `RingNeighbourhoodTrajectoryStats`") is described in the
  earlier section as aspirational with no current fields; the addition
  simply drops it. A cold reader noticing the "three shapes" claim and
  not finding a third concrete example in `TrajectoryAtom.h` would
  conclude the earlier section is describing a future state. With the
  brief's flag, I knew to treat the addition as current.

**Did reading both waste time?** Slightly. The addition is the
library-reference form I actually needed to add a TR. The earlier section
has design motivation (esp. §17's "duplicate rather than chain" rule +
worked case) that the addition compresses. If only one could remain, I
would want the addition plus §17's rationale section, not the earlier
structural prose.

**Would I have known which to trust without the task brief?** The header
on both additions states it explicitly ("Where this section disagrees
with ... trust this section as current truth"). So yes, the docs
self-disclose — but a first-time reader not primed by the brief will
likely read top-down and burn the context on the earlier section before
reaching the addition. Reordering or inlining the addition into the
primary section would save that cost.

---

## 5. Claims in the docs that contradict the code

Observed divergences:

1. **`doc/ARCHITECTURE.md:487` claims "54 headers in src/".** Actual count
   in `src/` is higher (88 per `INDEX.md:131`, which says "88 headers, 70
   .cpp files"). Not load-bearing, just stale.

2. **`EXTRACTION_ORDER.md:30-40` ("The pipeline resolves extraction order
   from the dependency graph") vs `src/OperationRunner.cpp:77-165`.** The
   code is a flat hand-ordered sequence of `TimedAttach` calls with skip
   flags, not a graph-resolving scheduler. The `Dependencies()` vector is
   checked at `ProteinConformation::AttachResult` for correctness (fail if
   missing) but does not influence ordering; ordering is hardcoded in
   `OperationRunner::Run`. A new contributor following
   EXTRACTION_ORDER.md literally would expect to just declare
   `Dependencies()` and be slotted in; in fact they must also modify
   `OperationRunner.cpp`.

3. **`OBJECT_MODEL.md:384-397` describes factory methods
   `protein.AddCrystalConformation(...)`, `protein.CrystalConformation()`
   (exactly one or throws), `protein.NMRConformations()`, etc.** The
   surrounding note at `:28-36` says "copy-and-modify pattern: superseded",
   which flags some of this as design history, but the factory methods
   themselves are not explicitly disclaimed. A reader trying to find
   `AddCrystalConformation` and finding instead that `BuildFromPdb`
   creates crystal conformations via a builder-function path (per
   `ARCHITECTURE.md:118-138`) will be momentarily confused.

4. **`OBJECT_MODEL.md:204-225` (AtomRole assignment table).** The table
   lists role assignment criteria; some rows mention residue index caches
   (`Residue.H`, `Residue.HA`) that are real, and others use phrasing like
   "PDB atom name" via hints in the last column. I did not verify these
   against `EnrichmentResult.cpp` and flag this only because a reader
   primed by PATTERNS.md anti-string-dispatch rules would wonder whether
   the "Residue.H index (not PRO)" column is really index-based or
   name-based. Not a contradiction I can confirm, just an ambiguity.

5. **`OBJECT_MODEL.md:1596-1601` (earlier) calls `ScanForDftPointSet`
   "slated for removal".** The 2026-04-24 addition at `:2798` and the
   PATTERNS.md addition repeat this. Three places, same statement. Not a
   contradiction — consistent. Flagging because the INDEX brief also
   mentions it, and the redundancy might be worth consolidating.

6. **Trajectory tests: documentation mentions `TEST_FRAMEWORK.md` and
   smoke_tests for output validation.** No test files in `tests/` target
   TrajectoryResults by name. The test surface is promised
   (`PATTERNS.md §14` mentions "testable in isolation"; §17 is explicit:
   "independent, testable in isolation") but not delivered in this
   directory. I did not open `tests/regression/` or smoke tests — they
   might cover this at the NPY round-trip level.

---

## 6. Areas where the docs clearly helped

- **`doc/ARCHITECTURE.md` section 5 "Relationship to PATTERNS.md"
  (`:461-484`).** The correspondence table mapping each pattern to a
  specific header/line was the fastest way to feel oriented. I used it
  to decide which headers to read first.
- **`PATTERNS.md §§1-12` for conformation scope.** Small, numbered,
  concrete. Each pattern has a code snippet. §11 ("Both representations
  always present") + the "T2 Completeness" block (`:520-561`) were the
  single most clarifying explanation of why every tensor has companion
  SphericalTensor and why scalar-only output is rejected.
- **`TrajectoryResult.h` header docstring (`:3-12`).** The invariant stated
  at the top ("accumulator state lives inside the TR subclass, not on
  TrajectoryAtom") made the AV/FO distinction click immediately when I
  then read BsWelford.
- **The "CROSS-RESULT READ" marker pattern in `BsWelfordTrajectoryResult.h`
  (`:7-17`) paired with `BsWelfordTrajectoryResult.cpp:44-67`.** Shows
  exactly how to declare, cross-read, and comment. If I needed to build
  a cross-reading TR, this is enough to work from.
- **`OBJECT_MODEL.md` Tier-2 field tables (`:658-788`).** The
  `Source result` column on every per-atom field is a massive aid to
  recovering "who writes this". PATTERNS.md §7 promised this; the tables
  deliver.
- **`CONSTITUTION.md §682-743` minimum per-calculator output contract.**
  Hard constraint, listed calculator-by-calculator, consistent with the
  anti-simplification rule at `:744-758`. Clear and terminal: if in
  doubt, this is what every calculator owes.
- **The 2026-04-24 addition in `OBJECT_MODEL.md` §§ TrajectoryProtein,
  TrajectoryAtom, TrajectoryResult, DenseBuffer.** As library-reference
  prose it is the right shape for the task. The H5 layout block at
  `:2953-2981` is exactly what a new TR author needs.

---

## 7. Other observations

- **`INDEX.md` is good but long.** The "When implementing a calculator"
  block at `:42-54` is a de-facto checklist I appreciated. A parallel
  "When implementing a TrajectoryResult" block would be useful and is
  conspicuously absent.
- **The `spec/INDEX.md:97-107` "All Documents" listing is a minefield of
  session-dated, partially-historical specs.** ENSEMBLE_MODEL.md flagged
  as "partially historical"; pending_include_trajectory_scope noted
  as "not authoritative"; pending_decisions_20260423.md referenced from
  `RunConfiguration.cpp:165-169` but not listed in INDEX.md at all.
  The presence of SESSION-dated landing-state + refactor-gaps docs
  indicates a recent major refactor still settling; a cold reader
  arriving today could waste significant time on history.
- **No documentation of the `spec/pending_decisions_20260423.md` file.**
  It is referenced from code comments (`RunConfiguration.cpp:165-169`,
  and the FullFatFrameExtraction note in OBJECT_MODEL.md `:2800`) but
  never listed in INDEX. If I were adding a TR that needed MOPAC I would
  need this doc; I cannot find it in the index.
- **Pleasant surprise: the `BsT0AutocorrelationTrajectoryResult.h` header
  (`:4-40`) is itself a worked-example template**, written to be cloned.
  The "clone this file and swap the scalar source" instruction is
  exactly what a new author wants to read. If every exemplar header had
  that paragraph, the extensibility story would tell itself.
- **The INDEX says "Single TPR parse" (`:148`) and describes production
  canonical stride=2.** Useful framing that lands the architecture in
  real fleet numbers without digging into ENSEMBLE_MODEL.
- **The single-target CMake discipline** (one library, `nmr_shielding`,
  all `.cpp` files listed explicitly) is stated in PATTERNS.md `:592-595`
  and verified in `CMakeLists.txt:195-242`. This is the sort of rule
  easy to miss and important to honour.
- **The "no utility namespace", "no template metaprogramming beyond
  Result<T>()", "no exception hierarchies", "no adapter/wrapper classes"
  anti-patterns** in PATTERNS.md are numerous and consistent. They
  reinforce a clear architectural tone. The prose about "the name is the
  thing" at `PATTERNS.md:9-16` is a distinct voice that sets
  expectations — one of the better project-intro opening paragraphs I
  have read.

---

## Overall impression

A reader could reasonably add either calculator type on the strength of
this documentation and the exemplar code, **provided they are flagged to
read the 2026-04-24 addition sections for trajectory-scope work and
provided they accept that tests for the trajectory scope do not yet
exist**. The gap between the conformation-scope docs (very mature, battle-
tested, redundant in a useful way) and the trajectory-scope docs (mature
but bifurcated, with design-history sections still live and no test
surface) is the single most visible artifact of being in the middle of
a recent refactor. The docs are honest about this — every divergence is
flagged at the top of its section — but a first-time reader would still
save meaningful context-window if the earlier sections were either
archived or compressed to a short "design history" appendix, and the
trajectory-scope docs promoted to the same first-class status as the
conformation-scope ones. The headline instruction "INDEX.md first" is
the right one; what follows rewards trust.
