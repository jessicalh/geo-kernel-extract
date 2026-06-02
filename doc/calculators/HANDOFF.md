# Handoff — calculator documentation effort (2026-06-02)

Where this stands and how to restart. The durable standards live in the briefs and in
the memory store; this is the orientation.

## What exists (committed)

`doc/calculators/` holds the documentation set, all on the shared `nmrdoc.sty` style
(Helvetica 11pt, thin margins, pastel page-breakable `codebox`):

- **12 per-calculator exposition docs** (`<slug>.tex/.pdf`) — the physics of each
  calculator: QM origin -> T0/T1/T2 -> method -> "geometric kernel, not shielding".
  Done, each authored and source-verified.
- **`use_guide`** — the `nmr_extract` operator guide. Done. It flags real
  `CLAUDE.md` 5-mode-spec vs code mismatches (single `--stride` not independent m/n;
  `.trr` required, not `.xtc` rejected; `--mopac` asymmetry; `--pH` parsed but unused) —
  worth reconciling spec or code someday.
- **`architecture`** — the object-model / runtime tech note. The approved version is
  committed. A **sober grace/expand WIP draft** is also committed (see below): it adds
  the *why* in three lean spots without removing anything.
- **`haighmallion-numerics`** — the FIRST numerics addendum and the **locked
  calibration target / exemplar** for the numerics series. Read it first; it is what
  "good" looks like for that series.
- **Briefs (the standards):** `BRIEF.md` (exposition author), `REVIEW_BRIEF.md`
  (verification pass), `USE_GUIDE_BRIEF.md`, `ARCHITECTURE_BRIEF.md`, and
  `NUMERICS_BRIEF.md` (distilled, quotes the exemplar).

## What's next

1. **Three more numerics addenda — the ring-current family: McConnell,
   RingSusceptibility, BiotSavart.** They share the most with HaighMallion (the dipolar
   kernel, the tensor decomposition, ring geometry), so the calibration transfers
   hardest.
   - **HIGH-TOUCH.** author -> review is NOT enough. Each needs a hands-on revision
     cycle: de-overwrought pass, composed-not-breezy register, three-beat glossary, a
     TikZ figure where the method is geometric, then an evaluative review-and-fix.
     HaighMallion took ~8 rounds.
   - Heavy codex spend -> fresh 5-hour window, a few at a time, not all at once.
   - Author/review prompt scaffolding is staged in `doc/calculators/.codex-*`
     (gitignored). An author reads `NUMERICS_BRIEF.md` + the exemplar.
2. **Object-model doc — Jessica leads the tone.** The committed approved version is
   excellent; the WIP draft adds the *why* (the identity/geometry/history separation;
   the buffer discipline; named operations). She sets the register by writing passages
   herself; match whatever she sets. **Default SOBER for this doc** — not the numerics
   addenda's warmer voice. (The first draft overshot into chipper and was pulled back.)

## The standard, in one breath (full version in `NUMERICS_BRIEF.md`)

Numerical-Recipes *register* but **transparency, not advice** (describe the code's
choices; never advise). **Tell how the numbers move** — no handwave, no
name-drop-and-defer. **Every descriptor earned by the topic** — the fizz in the drink,
nothing chipper or AI-flavored. **State exactly what the code proves** (fitted != used;
<= vs <; "rank at most one"; shielding != shift). Glossary entries run
**picture / example / formal**. A geometric method gets a **TikZ figure** (Mermaid is
for flowcharts, not geometry).

## Mechanics

- Compile: `pdflatex <slug>.tex` (twice if it has a `\tableofcontents`).
- Mermaid: `mmdc -i d.mmd -o d.pdf -p .puppeteerrc.json --pdfFit`.
- Codex: `codex exec --dangerously-bypass-approvals-and-sandbox "$(cat prompt.txt)"
  </dev/null > log 2>&1`, run in the background; kill by **precise PID**, never a broad
  `pkill` (other codex sessions run concurrently).
- The gitignored `prompts/` directory is the retrospective data bank — one entry per
  prompt: purpose / what's good / what we tuned / reusable lesson.

## Memory pointers

`project_numerics_addenda`, `project_calculator_exposition_docs`,
`feedback_technical_writing_standard` (includes the sober-architecture-tone note),
`feedback_plain_spoken_but_charming`.
