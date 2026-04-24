# Trajectory Framework Sketch

A clean-room implementation of the trajectory-scope framework described
by `spec/WIP_OBJECT_MODEL.md`. The sandbox is designed to be diffed
against `src/` so the reader can see which parts of `src/` are
pattern-essential (both places have them) versus residue from earlier
sessions (only in `src/`, not here).

Read the spec first. Do not read `src/` to calibrate against this
sketch — that would defeat the exercise.

## What is implemented

| File | Purpose (spec §) |
|------|------------------|
| `Stubs.h` | Minimal forward declarations / tiny stubs for `Protein`, `ProteinConformation`, `ConformationAtom`, `SphericalTensor`, `Vec3`, `Mat3`, `RunOptions`, and forward declarations for `BiotSavartResult`, `ConformationResult`, `AIMNet2Model`, `OperationRunner`, etc. Real bodies live in `src/`. |
| `DenseBuffer.h` | Template dense buffer + `DenseBufferBase` polymorphic base. Per §4 and Appendix C. |
| `FrameSelectionRecord.h` | The record shape `SelectionEmittingTrajectoryResult` emits. Per §5 and Appendix D. |
| `TrajectoryAtom.h` | Per-atom trajectory-scope data store. Private ctor, `friend TrajectoryProtein`. Only fields written by `BsWelfordTrajectoryResult` declared. Per §3. |
| `TrajectoryResult.h` | Base class with four virtuals (`Name`, `Dependencies`, `Compute`, `Finalize`) plus optional `WriteFeatures` / `WriteH5Group`. Per §4. |
| `SelectionEmittingTrajectoryResult.h` | Opt-in mixin interface. Per Appendix D. |
| `TrajectoryProtein.{h,cpp}` | Wraps Protein; holds `atoms_`; attached TrajectoryResults (singleton-per-class, attach-order preserved); adopted dense buffers. Per §3. |
| `RunConfiguration.{h,cpp}` | Typed first-class run shape. Three named static factories. Per §6. |
| `RunContext.h` | Class (not struct), constructor takes config + output dir; optional extras; no mutable cache. Per §6. Header-only because the spec's body is fully inline. |
| `Trajectory.{h,cpp}` | Process entity with 5-phase `Run()`. Per §5. Also performs Appendix D's selection-record collection in Phase 5. |
| `BsWelfordTrajectoryResult.{h,cpp}` | The one concrete worked example. Compute and Finalize bodies are the spec's §4 pseudocode verbatim. |

Fifteen files total. Five headers pair with four .cpps; the rest are
header-only where the spec puts the bodies inline.

## What is intentionally omitted

Things the spec explicitly says do not belong to the framework, or
does not demand:

1. **No `AllWelfords()` or any central field enumeration.** §3
   anti-pattern 3 is explicit: "Serialization is by each
   TrajectoryResult's WriteH5Group — no central enumeration, no
   parallel list to maintain."
2. **No `AccumulateFrame` on TrajectoryProtein.** §3 anti-pattern 2
   forbids the monolithic dispatcher. The per-frame loop lives in
   `Trajectory::Run` Phase 4 and dispatches polymorphically through
   `TrajectoryResult::Compute`.
3. **No per-field pointer list, no `std::vector<std::unique_ptr<Accumulator>>` on TrajectoryAtom.** §3 anti-pattern 1 names that
   "cleaner-C++" direction as the thing to avoid.
4. **No `TrajectoryBond` first-class store on TrajectoryProtein.** §3
   "Per-bond trajectory state: option (b) decided 2026-04-22" settles
   this: per-bond state lives internal to bond-scope TrajectoryResults,
   not as a parallel typed entity. Appendix B is preserved as the
   rejected option (a); the sandbox does not implement it.
5. **No `TrajectoryRing`.** §3's bond decision note: "No
   `TrajectoryRing` at this pass. ... Deferred."
6. **No `NmrAtomIdentity`, no `IupacAtomPosition` enum, no typed-slicing
   orthogonal indices on Protein.** §2 flags that subsystem as
   "PROPOSAL PENDING USER REVIEW" with "Do not implement §2 / Appendix
   A / Appendix E until user sign-off". Appendix H Stage A also
   "Deferred" to post-rollup. The sandbox respects that status.
7. **No NPY or H5 writers.** §7 commits to a writer shape (top-level
   `WriteTrajectoryH5` orchestrator that iterates attached Results)
   but the HDF5 emission is not the framework — it is writer concern.
   `Trajectory::WriteH5` is an empty body; `TrajectoryResult::WriteH5Group`
   has a default empty body. The framework commits to where H5 lives,
   not how it serialises.
8. **No concrete TrajectoryResult subclasses beyond `BsWelfordTrajectoryResult`.** The task says one worked example. Appendix F
   lists ~30 subclasses as an implementation checklist; that is
   catalog content, not framework. Adding them is the per-batch work
   of Stage 1b in Appendix H.
9. **No `GromacsFrameHandler` implementation.** §9 keeps it named
   format-specific ("stays Gromacs-named"). The sandbox forward-declares
   a minimal interface matching the five methods `Trajectory::Run`
   calls (`Open`, `Next`, `LastConformation`, `LastTimePs`,
   `LastXtcIndex`) and leaves them as no-op stubs in `Trajectory.cpp`.
10. **No `GromacsEnsembleLoader`.** §5 "EDR preload timing" names it
    as the thing that builds `TrajectoryProtein` from a TPR; the
    construction entry point is the single-arg `TrajectoryProtein`
    ctor here. A real loader wraps that.
11. **No `MopacVsFf14SbReconciliationTrajectoryResult` or any MOPAC-family
    class.** These are `FullFatFrameExtraction`-specific entries in
    Appendix F. Out of scope per point 8.
12. **No `EnrichmentResult`, no OpenBabel plumbing.** Appendix E
    explicitly decides no migration: "OpenBabel and EnrichmentResult
    stay unchanged." That is a decision about the *existing* code, not
    an addition to the framework.
13. **No tests, no CMakeLists.txt, no build.** Task said: "Do NOT commit
    to git. Do NOT build. Do NOT write tests."
14. **No catch blocks in framework code.** PATTERNS.md §9 "Return codes,
    not exceptions" for calculator/result/pipeline code. `Trajectory::Run`
    does `throw std::runtime_error` because the spec §5 concrete body
    does; see "Spec contradictions noted" below.

## Choices the spec did not fully specify

Each is the simplest thing that honours the rest of the design.

### 1. `std::span` not available in C++17

The spec snippets use `std::span<T>` in `DenseBuffer::AtomSlice`
(Appendix C) and `SelectionEmittingTrajectoryResult::SelectionRecords()`
(Appendix D). PATTERNS.md §"C++ Rules" says C++17 only.

- `DenseBuffer::AtomSlice` is not used by any framework-scope call
  site the spec shows (the two worked §4 Compute bodies use `At()`
  only). Omitted. `RawData()` + `AtomCount()` + `StridePerAtom()`
  cover the caller needs.
- `SelectionEmittingTrajectoryResult::SelectionRecords()` returns
  `const std::vector<FrameSelectionRecord>&` instead of
  `std::span<const FrameSelectionRecord>`. Same semantics for the
  single caller (the Phase 5 collection loop in `Trajectory::Run`
  per Appendix D).

Alternative considered: add `std::span`-shim. Declined — not in the
spec, adds complexity the framework does not require.

### 2. Attach-order storage

§3 code block names `std::vector<TrajectoryResult*> results_attach_order_`
as a member of TrajectoryProtein alongside the `unordered_map`. §5
"Ordering guarantees" makes this load-bearing. The sketch uses exactly
the pair (map for lookup, vector for order) the spec names. No
softening.

Alternative considered: keep only the map and sort at iteration.
Declined — §5 explicitly requires "attach order equals dispatch
order", which means "the order the caller attached them in Phase 1",
not a sort key.

### 3. `Compute` in `BsWelfordTrajectoryResult::MeanAtAtom` convenience methods

§4 shows `double MeanAtAtom(size_t atom_idx) const;` with no
TrajectoryProtein argument. The field it answers from lives on
TrajectoryAtom, which is owned by TrajectoryProtein. Spec doesn't show
a back-pointer on the Result. The sketch implements the method with
the spec signature and returns 0.0; noted in source that "a real
implementation would pair this with a back-pointer design the spec
does not commit to."

Alternatives considered:
- Add a TrajectoryProtein argument (changes the signature).
- Hold a back-pointer on the Result (the spec does not show one).
- Drop the method (the spec example shows it).

None is comfortable. The simplest-and-faithful choice is to keep the
signature and note it.

### 4. `AttachResult` dependency check

§4's attach code block shows a per-dep validation loop at attach time.
§12 item 4 RESOLVED: "stays in §5 Phase 2. Reason: attach doesn't need
to know about run configuration; Phase 2 is the one place with both
attached results and run config in scope." This sketch follows the
resolution: `AttachResult` does only the singleton check; all dep
validation (TrajectoryResult-typed and ConformationResult-typed) is
Phase 2's responsibility in `Trajectory::Run`. No attach-time dep
iteration — the spec resolution supersedes its own earlier code-block
example.

Alternative considered: keep the attach-time partial dep check for
already-attached TrajectoryResult types. Declined — either the check
belongs at attach (spec §12 resolution says no) or it doesn't.
Duplicating halfway is worse than full Phase 2.

### 5. `Compute` in `TrajectoryResult` base: virtual, not static

§4's prose explains this explicitly: "ConformationResult's Compute is
a static factory on the subclass because conformation-scope compute
is one-shot... For multi-phase streaming compute, the result is
constructed empty at attach and called polymorphically per frame."
Base has virtual Compute. `Create()` factories on subclasses stay
static. The sketch follows both.

### 6. `RunConfiguration`'s required set is empty in the sandbox

The spec's §6 factory bodies populate `required_conf_result_types_`
with things like `typeid(BiotSavartResult)`. `typeid` requires a
complete type, and `Stubs.h` forward-declares `BiotSavartResult`
without a body. The sketch leaves `required_conf_result_types_` empty
in all three factories; the framework shape still compiles and
Phase 2 dep validation still fires for declared TrajectoryResult
deps. Real code with real `BiotSavartResult` included populates this
set.

### 7. No `OperationLog::Log` calls in `AttachResult`

§4's concrete `AttachResult` body calls
`OperationLog::Log(OperationLog::Level::Warning, LogResultAttach,
"AttachResult", "rejected " + name + ": already attached");` on reject
and similarly on success. `OperationLog` is pre-existing infrastructure
(`src/OperationLog.h`), not framework surface. The sketch returns false
silently on reject and true silently on success. A real implementation
forwards-declares OperationLog alongside the other library forward-decls
in Stubs.h and calls through.

### 8. `Trajectory::WriteH5` and `TrajectoryResult::WriteH5Group` bodies are empty

Per §7 and spec §5, the framework commits to *where* serialisation
lives (on the owning entity) and *how the call flows* (top-level
writer iterates `ResultsInAttachOrder`). The actual HDF5 writes
depend on the HighFive library and on dataset-layout decisions
(Appendix C) that are writer implementation, not framework design.
The sketch leaves the bodies empty and documents the responsibility.

## Spec contradictions noted

### Exceptions vs return codes

PATTERNS.md §9 "Return codes, not exceptions" is strong: "Do not
create exception class hierarchies. Do not use std::expected. Do not
use try/catch in calculator or result code." `TrajectoryResult` /
`TrajectoryProtein` / `RunConfiguration` / `RunContext` subclasses are
all result/pipeline code under this rule.

WIP §5's `Trajectory::Run` concrete body throws `std::logic_error` and
`std::runtime_error`. WIP §3's `TrajectoryProtein::AttachResult` code
block in §4 returns `bool` but uses `OperationLog::Log` for
diagnostics — no throw there.

Choice: `TrajectoryProtein::AttachResult` returns `bool` (matches §3).
`Trajectory::Run` throws (matches §5's concrete code). The contradiction
is between two parts of the trajectory-scope spec itself; the sketch
sides with each section's own code, leaving the cross-section tension
visible rather than resolving it one way.

### `typeid` on forward-declared types

Dependencies() return types like `typeid(BiotSavartResult)`, and
§6 factories populate required sets with `typeid(BiotSavartResult)`
etc. These types are ConformationResult subclasses living in
`src/`. `typeid` needs a complete type. The sketch gets around this
by including only `typeid(BiotSavartResult)` (the one Dependencies()
this sandbox calls) — since the Dependencies return is compiled on
instantiation and BiotSavartResult is a forward declaration in
Stubs.h, `typeid` of it compiles as long as the translation unit is
pointed at the forward declaration (which works for a sketch, would
fail in a real link step).

This is a sketch-only convenience. In production code every cited
type must be a complete type at the point of `typeid` use.

### Zero `using namespace` in headers — kept

PATTERNS.md §"Headers". Followed.

### `auto` usage

PATTERNS.md: "Use for iterator types, make_unique results, and
range-for loop variables. Do NOT use for function return types."
Followed — `auto` appears only in range-for loops and
`std::make_unique` contexts in Trajectory.cpp and TrajectoryProtein.cpp.

### `ProteinConformation::AtomAt(i)` in Compute

§4 example reads `conf.AtomAt(i).bs_shielding_contribution.T0`. Stubs.h
provides a `ProteinConformation::AtomAt(i)` that returns a reference
to a sentinel `ConformationAtom`. The sketch compiles; the real type
returns a real atom. No framework difference.

## What to diff the sandbox against

When comparing `analysis-speculative/trajectory_framework_sketch/` to
`src/Trajectory*.{h,cpp}` / `src/RunConfiguration.{h,cpp}` /
`src/RunContext.{h,cpp}` / `src/BsWelfordTrajectoryResult.{h,cpp}` /
`src/DenseBuffer.h` / `src/TrajectoryMoments.h`:

1. **Things that appear in both** are probably spec-essential. The
   five-phase Run, the attach-order preservation, the singleton-per-class
   discipline, the Compute/Finalize virtuals, the private ctor on
   TrajectoryAtom with friend TrajectoryProtein, the DenseBuffer
   template, the SelectionEmittingTrajectoryResult mixin with its
   dynamic_cast collection loop.
2. **Things only in `src/` but not here** are candidates for "was
   this pattern-essential, or drift?" Typical suspects: centralised
   field enumeration methods (anti-pattern 3), `AccumulateFrame`-shaped
   dispatch (anti-pattern 2), accumulator-state on per-atom structs
   (anti-pattern 1), `TrajectoryBond` first-class store
   (rejected per §3 decision), `NmrAtomIdentity` typed slicing
   (Stage A deferred), MOPAC-family TrajectoryResult classes
   (not in this sketch but catalogued in Appendix F).
3. **Things only here but not in `src/`** are gaps in `src/` that
   the spec demands — if any. None expected, because the sketch is
   the floor of what the spec requires and `src/` should exceed
   that.

The spec is the authority. If `src/` diverges from the spec, the
diagnosis happens against the spec's text, not against this sketch
alone.

## One reminder

The worked example `BsWelfordTrajectoryResult` implements just the
pattern of a TrajectoryResult: a pair of factory + Compute + Finalize
bodies driving two of TrajectoryAtom's fields. The ~30 entries in
Appendix F follow the same pattern (AV Welford template) or the
§4 autocorrelation template (FO + dense buffer). This sketch does not
multiply them — that is Stage 1b grind.
