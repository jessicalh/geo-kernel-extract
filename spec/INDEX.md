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
- `spec/APPLIED_MATHEMATICS.md` — inventory of the numerical methods in
  `src/`, grouped by Numerical Recipes category, with `file:line` and which
  library (Eigen / nanoflann / libtorch / external) does the heavy lifting.

## Active planning (small)

- `spec/plan/README.md` — guide to what's live in `spec/plan/`.
- `spec/plan/docker-producer-handoff-2026-05-30.md` - Docker producer
  appliance handoff: PG18 build recipe, disk policy, runtime checks, and
  future maintainer notes.

The extractor is code-complete *as the existing producer*; a **forward
build** (kernel redesign, equivariant model, revamped stats) is live and
ahead of the rest of the tree. The older planned-calculator docs described
work now built and are retired history, not implementation authority.

## Forward build (live — read before kernel/model/stats work)

The current design effort lives under `doc/emerging/` and supersedes the
old "three stages" with **three Parts** (law study / tensor predictor /
shift predictor) and **The Three + the cage** kernel set. Reconcile by
reading code, not the older docs above.

- `doc/emerging/CONTROLLING_SPEC.md` — the spine: the three Parts, The
  Three, pinned-not-learned, the scope of the forward build.
- `doc/emerging/kernel_design/CONTINUITY.md` — running session state.
- `doc/emerging/DEFERRED_LEDGER.md` — every parked / contingent / follow-up
  item.
- `doc/emerging/kernel_design/` — per-kernel specs (ring / McConnell /
  charge-EFG) + the rhombic spec.

## Tracking

- `TENTATIVE_OUTSTANDING_ISSUES.md` — known issues, open decisions,
  TODOs, gotchas. Subsumes the former `KNOWN_BUGS.md`,
  `FIX_TESTS.md`, and `pending_decisions_20260423.md` (all bones'd).

## Subprojects

See **CLAUDE.md** §"Subprojects" — `ui/`, `h5-reader/`, `python/`,
`learn/`, `fileformat/`. Each that has its own `CLAUDE.md` owns
authoritative rules for work inside it.

## Archaeology

- Retired design docs, working notes, and session handoffs are prose
  history only. **Do not consult them to drive new work.** The decisions
  themselves are in `master`.
- `spec/meta-docs-review/` — 2026-04-03 documentation-audit artifacts.

## Discipline

The doc tree got chokingly large by mid-May 2026 (47 active spec/ +
spec/plan/ docs). The 2026-05-13 audit-and-bones pass reduced this
to ~10 active foundational docs + the planned-calc set. Future
sessions: prefer inlining into the foundational docs above OR
writing a memory entry. New standalone `spec/*` or `doc/*` files
are a high-friction move that creates the same problem again.
