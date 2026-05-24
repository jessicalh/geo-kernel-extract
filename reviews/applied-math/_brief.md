# Adversarial algorithm review — brief

You are reviewing one or more source files in the `nmr_shielding` C++ library
(geometric NMR chemical-shielding kernels on protein structures). The named
file(s) are inlined below, line-numbered. Review the inlined text; you are not
expected to read the repo.

## Standing assumptions — read first
- The code is believed **CODE-COMPLETE** and in production use.
- We **live with some shortcomings**. Changing an algorithm and re-validating
  it in our test + extraction pipeline is **very expensive**.
- So: **do NOT propose algorithm rewrites, redesigns, new abstractions, new
  helper layers, or library swaps.** The mathematical approach is fixed. This
  is not "you could instead use X."
- Prefer **"I don't know — check X"** over confident speculation about physics,
  conventions, or external behaviour. Do not invent bugs.

## The central question — spend most of your effort here

**Does this code tell a coherent story to a human, or is it only intelligible
to an AI or after long study?**

Read each file the way a chemist or physicist who knows the science — but not
this codebase — would read it, top to bottom, once. Then answer honestly:

- Can they follow the **through-line** — setup → the physical/mathematical
  steps in a sensible order → the result — without jumping around, holding a
  dozen cryptic one-letter symbols in their head, or reverse-engineering
  intent from the math?
- Or does it read as **symbol soup / premature compression / cleverness over
  clarity** — correct, but legible only to a machine or after an hour of study?

Quote the specific places where the story breaks, and for each say briefly
what it would take to make that passage read coherently (a rename, a named
intermediate, a regrouping, a signpost comment) — **without changing the
computation.**

## Review axes, in priority order

**1. Coherent story / readability (primary).** Does the file read as a
narrative? Flag: passages you had to decode, ordering that hides the logic,
dense expressions that fuse several meaningful steps, magic numbers without a
named meaning, control flow that obscures the math. For the worst offenders,
give a one-line "as written, a reader must …" so we feel the cost.

**2. Naming carries meaning.** Can a domain expert tell from each name what it
holds — the **quantity, its units, its frame**? Flag single-letter or
AI-abstraction names where a meaningful name exists; suggest one (1–3 words).
Names a physicist would recognize beat short names. Do not rename for style
where the current name is already clear.

**3. Visible math structure (grouping).** Are the stages grouped and sequenced
so the shape of the computation is obvious (e.g. build tensor → symmetrize →
de-trace → decompose, as distinct visible steps)? Flag tangled blocks, one-
liners that should be split into named steps, and intermediate quantities that
deserve a name.

**4. Function / method naming.** Does each name say what it computes and what
it returns (including units / which quantity)? Flag vague or misleading names;
suggest a clearer one.

**5. Comments as signposts.** We want **simple, direct, grounded, non-AI**
comments: a **2–4 word** label on each non-obvious math block (e.g.
`// traceless projection`, `// fan-triangulate ring`, `// reflection guard`,
`// near-field cutoff`). Flag: verbose/explanatory AI-style comments, comments
that restate the code, history/process prose, stale or wrong comments, and
non-obvious blocks with no label. Give the terse replacement.

**6. Correctness (secondary — still report it).** If you see an actual bug
(sign, convention, frame/transpose direction, index/bounds, units, missing
degeneracy guard, accumulation hazard), report it with `file:line`. Do not
hunt for correctness at the expense of axes 1–5, and do not manufacture
uncertainty into a "bug."

## Output format
Lead with a short verdict (2–4 sentences) on the central question for each
file: does it tell a coherent story, and if not, why. Then findings grouped by
the six axes above (readability first, correctness last). Each finding:
`file:line — issue — suggestion`. If a file is clean on an axis, say so in one
line. **Do not modify any files — review only.**
