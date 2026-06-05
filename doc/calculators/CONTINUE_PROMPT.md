# Continuation prompt — the thesis/book chapters

Paste the block below at the start of a fresh session. Full orientation is in `HANDOFF.md`,
the briefs, and the memory store; this is the directed prompt.

---

You are continuing a thesis/book documentation effort for the NMR shielding extractor, under
Jessica's imprint. As of 2026-06-04 the **three core chapters are done** — each truth-vetted,
polished, committed, with a `-backup.tex` beside it:

- **`physics-architecture.tex`** (`f785fe4`) — the unified view: the calculators as classical
  reductions of one quantum object (the Ramsey current response); the recurring kernel
  `D_ab = ∂_a∂_b(1/ρ)`.
- **`categorical-grounding.tex`** (`6ca47d6`) — the categorical view: by atom type and
  geometry; why each measure exists; the determinant literature; which calculator reaches for
  which nucleus.
- **`method-and-proxy.tex`** (`5c4bbf6`) — the intellectual shape: how each method PROXIES the
  physics, literal→metaphorical; the geometric kernel as the question posed to nature,
  calibration as the answer.

Also done/committed: 12 per-calculator exposition docs, 12 numerics addenda, the verified
reference markers (`backgrounder-refs/`), the briefs, and `OUTLINE.md` (rough notes only, NOT
a decided plan).

## Read first
1. `doc/calculators/HANDOFF.md` — full state + mechanics.
2. The three chapters above (read the PDFs).
3. Memory (start at repo root for the keyspace): `project_thesis_book`,
   `project_physics_backgrounders`, `feedback_draft_as_scaffold_her_words` (the OUTLINE
   boundary), `feedback_technical_writing_standard`, `feedback_single_focus_passes`,
   `feedback_correlate_not_match`, `feedback_run_the_algorithm_get_a_cookie`,
   `project_origin_true_tale`.
4. The briefs: `NUMERICS_BRIEF`, `BACKGROUNDERS_SEED`, `CHECK_BRIEF`, `EDITOR_BRIEF`,
   `GLOSSARY_REFINE_BRIEF`.

## The rules (load-bearing)
- **Voice:** Carl Sagan willing to use equations — every equation term-by-term, every variable
  grounded; the wonder is in the physics being true and connected, not in adjectives;
  respect-for-the-reader is the humanism. Audience: Jessica — fluent in maths and code; the
  PHYSICS is the gap, never the maths.
- kernels-not-shielding (geometric kernels until calibration); correlate, not pointwise match;
  state exactly what the code proves.
- **Every chapter goes through a two-pass truth-vet — one opus (physics + register), one codex
  (code-faithfulness + references) — BEFORE Jessica reads it ("before I fall in love"). Then,
  on her word, a light surgical tightening pass: backup first, edit in the doc's own voice,
  recompile, commit.** References: Crossref (DOI) + Europe PMC (abstracts); reasonable, not
  draconian; never guess a DOI; Semantic Scholar rate-limits (429).
- **THE OUTLINE IS DISCUSS-FIRST, NOT DRAFT-AND-PRESENT.** Chapters are scaffolds you may draft
  for her to cut against; the outline of the whole argument (intro spine, poster panels, lit
  review) is hers to architect, built TOGETHER, one-sentence telegraphic. Do not hand her a
  finished structure. (A session got ahead of this on 2026-06-04 and she pulled it back hard.)

## The framing (modest, by design)
Working title: **"Investigation of Classical Frame Calculations in an NMR Learning Model"** —
investigation-framed (the advisor's "a model that shows a result doesn't have to be good").
The arc is **static → dynamic → utility**: context → kernels-and-physics (the three chapters;
geometry / category / math illustrate strongly) → mutants → ridge → static contributors
(Stage 1) → 1P9J → in-frame, with a real equivariant model + PySR + the unifying equation
(Stage 2) → the actual model: does tossing these in the mix help, basically (Stage 3) →
architecture pictures. **Stage 3 is Stage 2 with a model attached** (+ maybe another learning
layer); nascent by EOD 2026-06-05 out of the rediscovery work. The quiet heart: PySR
rediscovering something shaped like the `D_ab` unifying kernel — the equation falling out of
the data, not imposed.

## Open work
- **The outline** (poster + thesis + telegraphic lit-review) — build TOGETHER, not unilateral.
- **30-second fix:** `physics-architecture.tex` line ~909 has the same kjær typo
  `134, 144104` → `044514` (Crossref-confirmed). Backup + correct + recompile + commit.
- **Stage 3** nascent by EOD 2026-06-05 (the in-frame equivariant + PySR work) → the
  "does it help" panel.
- **Per-calculator backgrounders** — designed and HELD (`BACKGROUNDERS_SEED.md`); a third
  per-calculator series (~2pp reaction-paper READMEs grounding the physics); awaits the framing
  the Stage-2/3 results may give.

## Mechanics
Codex at `gpt-5.5 xhigh`: `codex exec --cd <repo> -c model_reasoning_effort="xhigh"
--dangerously-bypass-approvals-and-sandbox "$(cat prompt)" </dev/null > log 2>&1`, background,
kill by precise PID, always `--cd` + absolute paths (Bash cwd persists between calls).
Editorial/vet taste = opus subagents (Agent tool, `model: opus`). Per-pass diff snapshots are
gitignored; `prompts/` is the retrospective data bank.
