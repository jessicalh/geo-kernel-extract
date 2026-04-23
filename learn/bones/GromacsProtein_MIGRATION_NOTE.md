# Migration: GromacsProtein family → TrajectoryProtein family

**Moved:** 2026-04-23 (trajectory-scope refactor, WIP §3 Stage 1 landing)

**Files in this directory:**
- `GromacsProtein.h`, `GromacsProtein.cpp` — old adapter
- `GromacsProteinAtom.h` — old per-atom Welford bag
- `GromacsRunContext.h`, `GromacsRunContext.cpp` — dissolved
- `GromacsFinalResult.h`, `GromacsFinalResult.cpp` — dissolved (CSV output)

**Why moved:** per `spec/WIP_OBJECT_MODEL.md` §3, these classes embodied four
anti-patterns:

1. Accumulator state (`Welford`, `DeltaTracker`, `TransitionCounter`) baked
   into the per-atom struct `GromacsProteinAtom`. Now replaced by
   `TrajectoryAtom` (typed summary fields only) with accumulator state
   living inside `TrajectoryResult` subclasses.
2. Monolithic per-frame dispatcher `GromacsProtein::AccumulateFrame` that
   knew every subsystem. Replaced by
   `TrajectoryProtein::DispatchCompute` which iterates polymorphic
   `TrajectoryResult::Compute` calls.
3. Hand-enumerated column list `GromacsProteinAtom::AllWelfords()` for
   serialization. Replaced by per-result `TrajectoryResult::WriteH5Group`
   (each result writes its own group).
4. Parallel `Gromacs*` typed hierarchy as a substitute for canonical
   trajectory-scope types. Replaced by `TrajectoryProtein` +
   `TrajectoryAtom` as peers of `ProteinConformation` +
   `ConformationAtom` at trajectory scope.

Additionally, `GromacsRunContext` dissolved: bonded parameters moved to
`TrajectoryProtein` (topology-scope, set once), preloaded EDR frames moved
to `Trajectory` (process-scope, bound to handler for per-frame lookup),
cursor state became private to `GromacsFrameHandler` (retained as the
XTC/TPR reader — format-specific role stays named after the format).

**Replacement map:**

| Old | New |
|---|---|
| `GromacsProtein` | `TrajectoryProtein` |
| `GromacsProteinAtom` (with Welford fields) | `TrajectoryAtom` (typed output fields only) |
| `GromacsProteinAtom::AllWelfords()` | per-`TrajectoryResult` `WriteH5Group` |
| `GromacsProtein::AccumulateFrame` | `TrajectoryProtein::DispatchCompute` |
| `GromacsRunContext::bonded_params` | `TrajectoryProtein::BondedParams()` |
| `GromacsRunContext::edr_frames` | `Trajectory::EnergyAtTime()` (preloaded) |
| `GromacsRunContext::AdvanceFrame` / `current_energy` | handler queries `Trajectory::EnergyAtTime(time_ps)` per frame |
| `GromacsFinalResult` (CSV writer) | dissolved; H5 via `TrajectoryProtein::WriteH5` + `Trajectory::WriteH5` |

The `GromacsFrameHandler` and `FullSystemReader` classes stay named
after the Gromacs format (they *are* format-specific readers) but now
operate on `TrajectoryProtein` instead of `GromacsProtein`.

**Byte-parity with pre-refactor output:** not preserved in the landing
commit. Only `BsWelfordTrajectoryResult` is implemented as a concrete
TrajectoryResult — the other ~40 Welford fields are no longer
accumulated. Restoring byte-parity requires migrating each remaining
Welford into its own TrajectoryResult subclass; follow-up sessions
work through `spec/WIP_OBJECT_MODEL.md` Appendix F catalog.

**See also:** `spec/WIP_OBJECT_MODEL.md` §3 (anti-patterns), §4
(TrajectoryResult pattern), §5 (Trajectory::Run), §9 (rename map),
Appendix F (TrajectoryResult catalog).
