# Document Index

Quick pointer for fresh readers. **Read CLAUDE.md first** — it's the
durable architectural record; this index is just a map to the rest.

## Foundationals (read when relevant)

- `README.md` — what the project is.
- `CLAUDE.md` — Claude-session reading order + project state.
- `doc/ARCHITECTURE.md` — system map with line-cited code pointers.
- `PATTERNS.md` — load-bearing patterns; anti-patterns; C++ rules.
- `OBJECT_MODEL.md` — typed object model; class/property/unit catalogue.
- `GEOMETRIC_KERNEL_CATALOGUE.md` — the 55 kernels.

## Spec foundationals

- `spec/CONSTITUTION.md` — supreme constraints.
- `spec/MATHS_GOALS.md` — what we're solving mathematically.
- `spec/PHYSICS_FOUNDATIONS.md` — physics underpinning each kernel.

## Active planning (small)

- `spec/plan/README.md` — guide to what's live in `spec/plan/`.
- `spec/plan/comprehensive-calculator-inventory-2026-04-30.md` — the
  full planned-calculator list. This is what we're building.
- `spec/plan/planned-calculator-substrate-audit-2026-05-06.md` —
  substrate ↔ planned-calculator mapping.
- `spec/PLANNED_CALCULATORS_2026-04-22.md` +
  `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` +
  `spec/POLARISABILITY_ROADMAP_2026-04-13.md` — earlier planned-calc
  snapshots; kept per user direction (planned-calc docs are exempt
  from retirement).

## Tracking

- `TENTATIVE_OUTSTANDING_ISSUES.md` — known issues, open decisions,
  TODOs, gotchas. Subsumes the former `KNOWN_BUGS.md`,
  `FIX_TESTS.md`, and `pending_decisions_20260423.md` (all bones'd).

## Subprojects

See **CLAUDE.md** §"Subprojects" — `ui/`, `h5-reader/`, `python/`,
`learn/`, `fileformat/`. Each that has its own `CLAUDE.md` owns
authoritative rules for work inside it.

## Archaeology

- `spec/plan/bones/` — retired design docs, working notes, session
  handoffs. **Do not consult to drive new work.** The decisions
  themselves are in `master`; bones is prose history only.
- `spec/meta-docs-review/` — 2026-04-03 documentation-audit artifacts.

## Discipline

The doc tree got chokingly large by mid-May 2026 (47 active spec/ +
spec/plan/ docs). The 2026-05-13 audit-and-bones pass reduced this
to ~10 active foundational docs + the planned-calc set. Future
sessions: prefer inlining into the foundational docs above OR
writing a memory entry. New standalone `spec/*` or `doc/*` files
are a high-friction move that creates the same problem again.
