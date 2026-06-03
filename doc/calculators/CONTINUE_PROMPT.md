# Catch-up prompt — calculator documentation effort

Paste the block below at the start of a fresh session to catch up on the docs. Full
orientation is in `HANDOFF.md`, the briefs, and the memory store; this is the directed
catch-up.

---

You are picking up a documentation effort for the NMR shielding extractor: short,
high-quality LaTeX/PDF technical notes under `doc/calculators/`, one per physics
calculator. As of **2026-06-03 the two big series are done and committed** (`2500bfa`,
branch `h5-reader-pysr-spike`, **not pushed**):

- 12 per-calculator **exposition docs** (`<slug>.tex/.pdf`) — the physics.
- 12 per-calculator **numerical-methods addenda** (`<slug>-numerics.tex/.pdf`) — the fancy
  maths, code-grounded. **HaighMallion is the locked exemplar.**

Plus the `use_guide` and the object-model `architecture` note.

## Read first, in order
1. `doc/calculators/HANDOFF.md` — full state, the four-pass pipeline that built the addenda,
   what's next, mechanics.
2. Memory (start the session at `/shared/2026Thesis/nmr-shielding/` for the root keyspace):
   `project_numerics_addenda`, `project_calculator_exposition_docs`,
   `feedback_technical_writing_standard`, `feedback_plain_spoken_but_charming`,
   `feedback_token_economy_codex_codes`.
3. The briefs (the standards): `NUMERICS_BRIEF.md`, `EDITOR_BRIEF.md`, `CHECK_BRIEF.md`,
   `GLOSSARY_REFINE_BRIEF.md`, exposition `BRIEF.md`.
4. The locked exemplar: `doc/calculators/haighmallion-numerics.tex` + `.pdf`.

## The standard (full in the briefs)
Numerical-Recipes register but **transparency not advice**; **tell how the numbers move**;
**every descriptor earned**; **state exactly what the code proves** (fitted != used; <= vs
<; "rank at most one"; shielding != shift). Glossary first beats run **picture +
grounding** — grounded metaphor, then the precise description in parens (the exemplar's
bar-magnet). Geometry gets a **TikZ figure** (Mermaid is for flowcharts). The code is the
final authority. Format: `\usepackage{nmrdoc}` + the pastel breakable `codebox`; hand-fit
code lines; keep it short.

## Open work
- **Jessica's human pass** over the 12 addenda — light hand edits. This is hers, not an
  agent's. Per-pass diff snapshots are in `.pre-check/` / `.pre-glossrefine/` /
  `.pre-opus-edit/` (gitignored) if you need to compare a passage across passes.
- **Physics backgrounders** — the next series, **scope TBD as of 2026-06-03**. Do not start
  without scoping with Jessica first.
- **Object-model tone pass** — `architecture.tex`; Jessica leads by hand, **sober**
  register (not the numerics' warmer voice).

## Mechanics
Codex authors/checks at `gpt-5.5 xhigh`: `codex exec --cd <repo>
-c model_reasoning_effort="xhigh" --dangerously-bypass-approvals-and-sandbox "$(cat prompt)"
</dev/null > log 2>&1`, background, kill by **precise PID** (never a broad `pkill`). Use
`--cd` and absolute paths — the Bash cwd persists between calls and a stray `cd` makes codex
silently no-op. Editorial / glossary taste is done by **opus subagents** (the Agent tool,
`model: opus`) with the relevant brief; codex does the source-grind. Token economy: the
orchestrator saves tokens for briefs, judgment, and reading reports.
