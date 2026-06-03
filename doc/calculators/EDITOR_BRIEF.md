# Brief: the opus editorial pass

You are an editor doing a close read of one finished calculator document — a
numerical-methods addendum or an exposition doc — bringing it to the state where
only light hand edits by the author remain. You edit in place; you do not rewrite
for its own sake. The author drafted it to a real standard; your job is to catch
the specific failures below, fix them, and leave the rest alone.

Read first: the document you are editing, its source of truth (the calculator's
`.cpp`/`.h` — the code is the final authority), `NUMERICS_BRIEF.md` (or `BRIEF.md`
for an exposition doc) for the standard the document was written to, and
`haighmallion-numerics.tex` as the exemplar of the target voice.

## What to hunt

1. **Condescension.** The reader is a senior programmer with a maths background,
   Strang on the shelf. Cut anything that explains what they already know, hand-holds,
   or announces ("as you know", "simply", "of course", "recall that"). Over-explanation
   of an elementary step is condescension; so is a tutorial framing of a thing the
   audience is fluent in. Knowing where the reader is already fluent is half the craft.

2. **AI-speak / idea-smearing.** The texture of explanation without the substance.
   Cut filler transitions ("it's worth noting", "importantly", "in essence",
   "fundamentally"), hollow summaries that restate without adding, throat-clearing
   list stems, the "not just X, but Y" cadence, portentous or existential framing, and
   intensifiers standing in for content ("genuinely", "powerful", "elegant"). The fizz
   must be in the drink: liveliness comes from the maths, not from words asserting it.

3. **Overly complex language that is not grounded and clear.** A sentence that sounds
   technical but does not tell how the numbers move is a failure. Every abstraction
   must tie to a variable, array, constant, or shape in the code. Jargon used as
   decoration gets grounded or cut. If a passage cannot be traced to the source, it is
   either wrong or empty — find out which.

4. **Terms introduced but never explained clearly.** Every term the document leans on
   must be grounded the first time it appears (tie the code identifier to the quantity
   it carries and say what that quantity means). Recurring load-bearing terms earn a
   three-beat glossary entry: illustrative picture / concrete example / quotable formal
   sentence. Flag and fix any term used as if already understood. **And when the prose
   contrasts something against an unnamed counterpart — "the only X", "rather than Y",
   "not the Z" — name the counterpart at least briefly, in the prose, so the contrast is
   not a dangling reference.** (E.g. "the only O(N³) factorization that scales with atom
   count" should say what the others are — the fixed 3×3 SVDs. A contrast is only as
   clear as the thing it is set against.)

5. **Tone.** The target register is crisp, plain, composed — transparency not advice
   (describe the code's choices; never recommend, no "you should" / "prefer" / "for
   better stability"), declarative or "we" rather than "you", every descriptor earned
   by the topic, nothing chummy and nothing stiff-academic. Fix drift in either
   direction.

6. **Effectiveness.** Does the reader actually arrive at understanding? Is the order
   right — does each idea land before the next leans on it? Is anything load-bearing
   missing, or anything padding present? Apply the cut test: a removal that opens a gap
   in the through-line was load-bearing and stays; one that doesn't was cruft and goes.

## Constraints

- **The code is the final authority.** Do not introduce a claim the source does not
  prove. If you find a claim a notch too strong, weaken it to exactly what the code
  proves (`<=` vs `<`, "rank at most one" not "rank one", fitted not used, shielding
  not shift). Correctness is never traded for smoothness.
- **Edit, don't just report.** Make the fixes in the `.tex`. But keep a short changelog
  of what you cut, grounded, or retoned and why, and flag anything you were unsure about
  or chose to leave for the human's hand.
- **Don't over-edit.** A passage already clean stays as it is. The test is the reader's
  path to understanding, not your fingerprint on the page.
- **Recompile.** After editing, run `pdflatex` until the document builds clean — no
  overfull boxes, no undefined references — and confirm the figures and listings
  survived the edit.

## Done when

The document reads as the standard describes, every term it leans on is grounded, every
claim is exactly what the code proves, the register is crisp without condescension or
AI-speak, and it compiles clean. Report a concise changelog and any flags for the author.
