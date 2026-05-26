# Planning Artifacts

Small directory. Contains:

- `cleanup-config-cmake-audit-2026-05-25.md` — working inventory for the
  general cleanup, runtime/config cleanup, and CMake/install cleanup pass.
- `proper-program-phase-2026-05-25.md` — working tracker for the
  pre-Doxygen "proper program" doc/comment cleanup phase.
- `bones/` — retired design records, working notes, session handoffs.
  **Archaeology only — do not consult to drive new work.** Decisions
  themselves are in `master`; bones is prose history.

That's it. The project is code-complete: the planned-calculator docs
(inventory, substrate audit, the PLANNED_CALCULATORS /
POLARISABILITY_ROADMAP snapshots) and the landed design notes
(welford-data-shape, test-suite-realignment, adversarial-review-prompt)
were retired to `bones/` once the work landed. Nothing is exempt.

## Discipline

New planning docs are a high-friction move. Prefer:

- Inlining into `OBJECT_MODEL.md` / `PATTERNS.md` / `CONSTITUTION.md`
  (or the matching subproject's `CLAUDE.md`).
- Writing a memory entry.
- For in-flight refactors: keep working notes in conversation, bones
  once the work lands.

A new standalone planning doc lives here only on explicit user request,
and gets bones'd on landing. The pattern of leaving session-handoff
docs around becomes its own debt within weeks.
