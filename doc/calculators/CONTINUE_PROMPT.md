# Continuation prompt — calculator documentation effort

Paste the block below at the start of a fresh session to continue the work.
Orientation and detail are in `HANDOFF.md` and the briefs; this is the directed prompt.

---

You are continuing a documentation effort for the NMR shielding extractor: short,
high-quality LaTeX/PDF technical notes under `doc/calculators/`, one per physics
calculator, plus an object-model architecture note. The physics exposition docs are
done. The active work is the **numerical-methods addenda** (one per calculator) and a
tone pass on the **object-model** note.

## Read first, in order
1. `doc/calculators/HANDOFF.md` — current state, what's next, mechanics.
2. `doc/calculators/NUMERICS_BRIEF.md` — the standard for the addenda (quotes the exemplar).
3. `doc/calculators/haighmallion-numerics.tex` + `.pdf` — the **locked calibration
   target**; this is what "good" looks like. Read it before writing.
4. Memory (start the session at `/shared/2026Thesis/nmr-shielding/` for the root
   keyspace): `project_numerics_addenda`, `project_calculator_exposition_docs`,
   `feedback_technical_writing_standard`, `feedback_plain_spoken_but_charming`,
   `feedback_token_economy_codex_codes`.
5. Per addendum: the calculator's source (the final authority) and
   `spec/APPLIED_MATHEMATICS.md`.

## The rules (full version in `NUMERICS_BRIEF.md`)
- Numerical-Recipes **register**, but **transparency, not advice** — explain how *this
  code's* numbers move; describe its choices as fact; never advise ("you should",
  "prefer", "for better stability").
- **Tell how the numbers move** — no handwave, no name-drop-and-defer.
- **Every descriptor earned by the topic.** The fizz is in the drink. Nothing chipper,
  existential, or AI-flavored. "We" not "you". Ground variables; don't announce them.
- **State exactly what the code proves**, never a notch stronger: fitted != used, <= vs
  <, "rank at most one", shielding != shift.
- **Glossary = three beats**: illustrative picture / a real example / a quotable formal
  definition.
- **Geometry gets a TikZ figure** (Mermaid is for flowcharts, not geometry).
- **The code is the final authority.** Format: `\usepackage{nmrdoc}` (shared style +
  pastel breakable `codebox`); load `siunitx`/`bm`/`tikz` as needed; hand-fit code
  lines; keep it short.

## The task and how to work it
- **Three more numerics — the ring-current family: McConnell, RingSusceptibility,
  BiotSavart** (most shared with HaighMallion). One or two at a time.
- **HIGH-TOUCH:** codex author -> codex review is NOT enough. Each needs a hands-on
  revision cycle — author, then a de-overwrought / register / glossary / figure
  grace-pass by hand, then an evaluative codex review (truth / accuracy / writing),
  then fixes. HaighMallion took ~8 rounds. Stay in the loop with Jessica on register
  and taste.
- **Object-model note:** the committed version is approved; a sober grace/expand draft
  exists. Jessica sets that tone by hand — match her register, keep it **sober** (not
  the numerics' warmer voice).
- **Token economy:** codex does the authoring/verifying grind; the orchestrator saves
  tokens for briefs, judgment, and reading reports. Launch codex in the background; kill
  by **precise PID**, never a broad `pkill` (other codex sessions run concurrently).
  Heavy spend -> fresh 5-hour window, a few docs at a time.

When you have picked a calculator, stage its author prompt (the existing
`doc/calculators/.codex-*` scaffolding is the model), launch codex, and run the cycle.
