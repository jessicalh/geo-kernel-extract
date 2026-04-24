# Pending decisions — 2026-04-23

Items surfaced during the trajectory-scope doc wrongness pass that need
user choice or defer-to-later implementation. Each entry: what, current
state in the tree, user's current direction (if any), priority.

Companion to `spec/doc_wrongness_20260423.md`. The wrongness file lists
contradictions for in-place fixing; this file lists items where fixing
requires a call rather than just removing the stale claim.

---

## 1. AllWelfords revival

- **What**: revive `PATTERNS.md §24`'s `AllWelfords` enumeration on
  `TrajectoryAtom` (or equivalent). Single source of truth for per-atom
  Welford field naming + pointer-to-member; the writing TR filters by
  its own `type_index` when serialising.
- **Scope**: Welford-specific, as the original §24 had it. Not a general
  "named fields" abstraction.
- **Current state**: pattern is documented in `PATTERNS.md §24` as live,
  but the class it attached to (`GromacsProteinAtom`) moved to
  `learn/bones/` during the trajectory refactor. No equivalent on the
  current `TrajectoryAtom`. Each `*WelfordTrajectoryResult::WriteH5Group`
  spells out its dataset names by hand — the drift risk §24 was
  solving.
- **Direction**: not yet chosen.
- **Priority**: low. No blocking consumer today; becomes worth it as
  Hm / McConnell / Sasa Welford TRs land.

## 2. ChiRotamerSelectionTrajectoryResult — is this the right scan-trajectory shape?

- **What**: class lives at
  `src/ChiRotamerSelectionTrajectoryResult.{h,cpp}` but is not attached
  by any `RunConfiguration` factory. `ScanForDftPointSet()` attaches
  only `BsWelford`. Currently orphaned.
- **Current state**: class is a scan-mode selection emitter; pushes
  per-residue chi rotamer transitions into
  `traj.MutableSelections()`. Exists as code but is not exercised in
  production or in any attached configuration.
- **Direction**: unsure — may not be the right selection-trajectory
  shape for the scan-for-DFT use case at all.
- **Priority**: sit. Scan-for-DFT bit is ≥2 weeks out soonest (μs-MD
  harvester timeline).

## 3. FullFatFrameExtraction MOPAC ConformationResult dependencies

- **What**: `RunConfiguration::FullFatFrameExtraction()` sets
  `skip_mopac = false` but does not add `MopacResult`,
  `MopacCoulombResult`, or `MopacMcConnellResult` to
  `required_conf_result_types_`. An attached MOPAC-family
  `TrajectoryResult` declaring those as `Dependencies()` would fail
  Phase 4 validation today.
- **Current state**: deferred in the code. Comment notes the set
  "just lists more than required" is harmless — true — but the
  converse ("just lists fewer than required") is what would bite when
  a MOPAC-family TR is added.
- **Direction**: include MOPAC deps. Won't run often (MOPAC every
  frame is ~15 h per protein, so `FullFatFrameExtraction` is for
  selected-frame MOPAC only, which makes the cost acceptable).
- **Priority**: medium. Wire before any MOPAC-family
  `TrajectoryResult` lands.

---

END — extend with other pending items as the wrongness-fix pass identifies them.
