# Planning Artifacts

Small directory. Contains:

- `comprehensive-calculator-inventory-2026-04-30.md` — full
  planned-calculator inventory. This is what we're building.
- `planned-calculator-substrate-audit-2026-05-06.md` — substrate ↔
  planned-calculator mapping.
- `bones/` — retired design records, working notes, session handoffs.
  **Archaeology only — do not consult to drive new work.** Decisions
  themselves are in `master`; bones is prose history.

That's it. The 2026-05-13 audit-and-bones pass moved all other
session-time planning artifacts to `bones/`. Earlier planned-calc
snapshots live at the spec root: `spec/PLANNED_CALCULATORS_2026-04-22.md`,
`spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`,
`spec/POLARISABILITY_ROADMAP_2026-04-13.md`.

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
