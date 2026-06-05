# Codex — adversarial critique of the Stage-1 registry plan

## (Preamble prepended. VERIFY FROM THE CODE BEFORE TRUSTING THE PLAN, EVER — cite file:line. Read-only: NO code, NO git, NO build.)

You are the adversarial check before the lead commits to an architecture plan an opus agent
wrote for introducing a typed visualisation registry into this Qt6/VTK reader. The author is a
collaborator; your job is to find where the plan is WRONG, INFEASIBLE, or UNDER-SPECIFIED —
before we build from it. Verify every load-bearing claim against the actual code; do not take
the plan's word.

## Read
- `notes/STAGE1_REGISTRY_PLAN_2026-06-05.md` — THE plan to critique.
- `notes/VISUALISATION_AUDIT_AND_MUSTHAVES_2026-06-05.md` — the audit it builds on.
- The cited code: `src/model/DisplayModeCapability.h`, `DashboardSignal.{h,cpp}`,
  `DashboardPanelModel.cpp`, `DashboardDisplayController.{h,cpp}`, `SignalDisplayDialog.cpp`,
  `DashboardSelectionController.cpp`, `TrajectorySignalCatalog.cpp`, `AbstractStripPanel.h`,
  `SceneRevealOverlay.h`, `TrajectoryFieldAvailability.h`, `CMakeLists.txt`.

## Verify the plan's load-bearing claims FROM THE CODE (confirm or refute, file:line)
1. `DisplayModeCapabilityFor` is consulted at ~6 sites and can be WRAPPED (body re-pointed) with
   all consumers unchanged. List the actual call sites; confirm wrapping preserves behaviour.
2. The three flags (`hasVisibleSurface`/`buildsPanelWidget`/`emitsPanelRef`) are each read by a
   DIFFERENT consumer; the registry must preserve all three. Confirm which consumer reads which.
3. `static.tensor` is "registered-but-inert" (`emitsPanelRef=true`, builds nothing) — a naive
   registered⇒offerable rule wrongly exposes it. Confirm.
4. The removal cascade (`DashboardSelectionController`) is bidirectional + re-entrancy-guarded
   (`signalsBeingRemoved_`) and refs are keyed by `displayModeId` STRING; Stage 1 keeps the
   string as the ref key so the cascade + strip-history persistence are byte-for-byte unchanged.
   Confirm this invariant actually holds under the plan's changes.
5. The empty-reality-check uses `TrajectoryFieldAvailability` canonical+storage-path lookup +
   `Absent`/`AllMissing` states + an `alternatesTried` list. Confirm the availability table
   actually exposes what's needed, and that the UDP / `smokeSummary` log path exists.

## Critique (find what's wrong or missing)
- Are the 8 migration steps ACTUALLY independently buildable/green, or does any have a hidden
  ordering dependency / a step that won't compile alone?
- The startup-validation rule: will it FALSE-POSITIVE (refuse to start) on descriptors whose
  data is legitimately absent for THIS run? The empty case must be a runtime LOG, not a
  startup FAIL — confirm the plan draws that line correctly and won't brick a normal load.
- Registry singleton + `VisualizationDefinition` lifetime/ownership/thread: static-init order,
  who owns the definitions, when registered, `ASSERT_THREAD` correctness. Sound?
- What landmine did the plan MISS beyond its 7? (lifetime, static init, the dialog substring
  matching, panel-ref keying, MOC/CMake additions, reveal-overlay coupling.)
- Anything that would make a codex IMPLEMENTING this stumble or sprawl.

## Output
Write `notes/STAGE1_REGISTRY_PLAN_CRITIQUE_2026-06-05.md`: findings (file:line, severity
CRITICAL/MAJOR/MINOR, the fix), the missed landmines, and a one-line verdict —
PLAN-IS-SOUND / FIX-PLAN-FIRST (list must-fixes) / RETHINK. Read-only. No code. No git.
