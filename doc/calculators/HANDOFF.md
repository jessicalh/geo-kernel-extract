# Handoff — calculator documentation effort (2026-06-03)

Where this stands and how to restart. Durable standards live in the briefs and the memory
store; this is the orientation.

## What exists (committed)

`doc/calculators/`, all on shared `nmrdoc.sty` (Helvetica 11pt, thin margins, pastel
page-breakable `codebox`):

- **12 per-calculator exposition docs** (`<slug>.tex/.pdf`) — the physics of each
  calculator: QM origin -> T0/T1/T2 -> method -> "geometric kernel, not shielding". Done,
  source-verified.
- **`use_guide`** — the `nmr_extract` operator guide. Done. (Flags real `CLAUDE.md`
  5-mode-spec vs code mismatches: single `--stride` not independent m/n; `.trr` required;
  `--mopac` asymmetry; `--pH` parsed-but-unused.)
- **`architecture`** — object-model / runtime note. Approved version committed, plus a
  sober grace/expand WIP draft. **The object-model tone pass is still open** — Jessica
  leads it by hand, **sober** register (not the numerics' warmer voice). The first draft
  overshot chipper and was pulled back.
- **12 numerical-methods addenda** (`<slug>-numerics.tex/.pdf`) — **DONE and committed**
  (`2500bfa`, branch `h5-reader-pysr-spike`, **not pushed**). One per calculator: the fancy
  maths given a Numerical-Recipes-grade, code-grounded, code-review treatment.
  **HaighMallion is the locked exemplar** — read it first; the other eleven were calibrated
  to it.
- **Briefs (the standards):** exposition `BRIEF.md` + `REVIEW_BRIEF.md`;
  `USE_GUIDE_BRIEF.md`; `ARCHITECTURE_BRIEF.md`; numerics `NUMERICS_BRIEF.md` (author
  standard, quotes the exemplar); `EDITOR_BRIEF.md` (opus editorial pass); `CHECK_BRIEF.md`
  (the four-lens final codex check); `GLOSSARY_REFINE_BRIEF.md` (metaphor-then-grounding).

## How the 12 addenda were produced (the pipeline)

The early plan was HIGH-TOUCH per-doc (HaighMallion took ~8 hand rounds). That was retired
for the rest in favour of *powering through* with a four-pass pipeline:

1. **Author** — codex `gpt-5.5 xhigh` from `NUMERICS_BRIEF.md` + the exemplar. (The ring
   family got a normal-effort author then an xhigh improve pass; the other eight a single
   xhigh author.)
2. **Opus editor** — one opus subagent per doc (`EDITOR_BRIEF.md`): condescension /
   AI-speak / ungrounded complexity / unexplained terms / tone / effectiveness. Caught two
   real factual errors.
3. **Codex check** — one codex `xhigh` per doc (`CHECK_BRIEF.md`), full source license:
   (1) facts, (2) science-writing tone (cut loose "hard"; keep terms of art),
   (3) read-as-the-reader, (4) "would the sentence stand if the topic were beets or car
   parts". Many factual fixes.
4. **Glossary-refine** — opus subagents (`GLOSSARY_REFINE_BRIEF.md`): pass 3 had stripped
   the picture-first glossary metaphors; the resolution is **metaphor + precise grounding
   in parentheses** (the exemplar's small-bar-magnet pattern), with the overwrought
   metaphors trimmed.

Per-pass pristine snapshots are in `doc/calculators/.pre-check/`, `.pre-glossrefine/`,
`.pre-opus-edit/` (gitignored) for diffing any entry both ways.

## What's next

1. **Jessica's human pass** over the 12 addenda — light hand edits. The automated passes
   have taken them as far as they sensibly go; the floor is hers.
2. **Physics backgrounders** — a new series, **scope TBD as of 2026-06-03**. The next major
   work after the addenda settle. Do not start without scoping with Jessica.
3. **Object-model tone pass** (above) — still open; Jessica leads, sober.

## The standard, in one breath (full in the briefs)

Numerical-Recipes *register* but **transparency not advice**; **tell how the numbers
move** (no handwave, no name-drop-and-defer); **every descriptor earned**; **state exactly
what the code proves** (fitted != used; <= vs <; "rank at most one"; shielding != shift).
Glossary first beats run **picture + grounding** — a grounded metaphor, then the precise
description right after (often in parens); the exemplar's bar-magnet entry is the model.
Geometry gets a **TikZ figure** (Mermaid is for flowcharts). Recast rep-theory shorthand to
plain ("the nine components of a 3x3 tensor split as 1+3+5"); name a contrast's unnamed
counterpart in the prose.

## Mechanics

- Compile: `pdflatex <slug>.tex` (twice if it has a `\tableofcontents`).
- **Codex:** `codex exec --cd <repo> -c model_reasoning_effort="xhigh"
  --dangerously-bypass-approvals-and-sandbox "$(cat prompt)" </dev/null > log 2>&1`, run in
  the background; kill by **precise PID** (never a broad `pkill` — concurrent codex sessions
  run). **Use `--cd` and absolute paths**: the Bash tool's cwd persists between calls, so a
  stray earlier `cd` makes relative-path redirects resolve wrong and codex silently no-ops
  (no log file, exit 0).
- **Opus passes** (editorial / glossary taste) are spawned as opus subagents (the Agent
  tool, `model: opus`) with the relevant brief — taste is opus's job; codex does the
  source-grind.
- The gitignored `prompts/` directory is the retrospective data bank.

## Memory pointers

`project_numerics_addenda` (updated: pipeline + status + codex mechanics),
`project_calculator_exposition_docs`, `feedback_technical_writing_standard` (sober-
architecture note + the glossary metaphor-then-grounding lesson),
`feedback_plain_spoken_but_charming`, `feedback_token_economy_codex_codes`.
