# Brief: glossary refinement — metaphor, then grounding

The last codex check stripped the picture-first metaphors out of these glossaries and left
precise, definitional openings. The author's decision: **keep the picture, but ground it.**
Open each entry with a grounded, direct metaphor, then give the precise description right
after — the correction sitting after the metaphor, often in parentheses. The model is the
exemplar, `doc/calculators/haighmallion-numerics.tex`, dipolar-kernel entry: "The field
shape of a small bar magnet — strongest and signed along its axis, weaker and reversed to
the side, fading as the cube of distance." Picture first, cashed out into exact behaviour in
the same breath.

You have two sources for each entry:
- `doc/calculators/.pre-check/<slug>-numerics.tex` — the earlier picture-first openings (raw
  material for the metaphor).
- `doc/calculators/<slug>-numerics.tex` (the live file you edit) — the check's precise
  openings (the grounding to keep).

For each glossary entry, produce a first beat that is **metaphor + grounding**: a grounded,
direct picture, then the precise description right after it (parenthetical or as the next
clause). Keep the second beat (concrete example) and the third beat (formal sentence)
intact, except for the recast named below.

## Two cautions from the author

1. **Some of the earlier metaphors were too much.** Be careful, grounded, and direct. Where
   the `.pre-check` opening was cute or overwrought — "a note pinned to a residue", "a
   levelling pass … shaves", "one arrow laid against itself", "the point fiction blows up" —
   do not restore the excess; write a restrained, direct picture and then ground it. A
   metaphor earns its place only if it carries the maths. If it only decorates, cut it and
   keep the grounding.

2. **The formal spherical-tensor beat.** Wherever the third beat says "the closed-form
   decomposition of `3⊗3` into dimensions 1+3+5", or calls `T2` "the irreducible rank-2
   component … under rotations", recast to the grounded form: "the closed-form split of a
   3×3 tensor's nine components into dimensions 1+3+5, with the `T2` basis normalized to
   preserve the Frobenius norm." Drop the `⊗` notation and "irreducible under rotations" —
   name the nine components.

## The recurring spherical-tensor entry — use this opening verbatim

Most of these documents carry the shared spherical-tensor entry. Give its first beat this
exact metaphor-with-grounding opening, so it reads uniformly across the series:

> A 3×3 tensor holds three kinds of content that rotate without mixing — an orientation-free
> size (the scalar trace, `T0`), an axial twist from its antisymmetric part (the three `T1`
> components), and a five-number traceless shape for the rest of the direction-dependence
> (the five `T2` components).

Keep each document's own second beat (which channels *this* calculator populates) and its
third beat (recast per caution 2).

## Done

Change nothing outside the Glossary section except a specific body insert named in your
assignment. Recompile each doc with `pdflatex` until clean (no overfull boxes, no undefined
references). Report per doc: the entries you reworked (short before/after of the first
sentence) and any you judged should stay as the check left them.
