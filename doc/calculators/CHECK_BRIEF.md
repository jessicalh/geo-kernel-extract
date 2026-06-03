# Brief: the final codex check

This is the last automated pass before the author's own hand edit. Four questions, in
order. The first is a correctness gate; the rest ask whether the prose earns its place as
science writing for the reader we have in mind. **Fix what is unambiguous; flag the
judgment calls for the author's hand** — and when you flag, be specific: quote the
sentence, name the problem.

You have **full license to go check the source.** The calculator's `.cpp`/`.h` (and the
constants in `data/calculator_params.toml`, `src/PhysicalConstants.h`, the shared
decomposition in `src/Types.h`) are the final authority. Read whatever you need to confirm
a claim.

## 1. Is everything factually true?

Go verify. Every constant, shape, inequality (`<=` vs `<`), sign, axis range, array
dimension, and "the code does X" claim must match the source as it stands. Where a claim is
a notch stronger than the code proves, weaken it to exactly what the code proves. Where it
is simply wrong, fix it. A documentation claim the code does not back is the worst defect
here — hunt it.

## 2. Is the tone right for science writing — even in a README?

These are extended READMEs with science content: not a textbook, not a blog post, but what
a careful science writer would put in a serious README. The failure to catch is the loose
evaluative word that smuggles in an unstated metric or metaphor. The canonical example is
**"hard."** A science writer does not call something "hard" and leave it — they either
(a) name a metric of the hardness and explain it briefly (in parentheses), or (b) make the
metaphor explicit. Either way it is a lot to say, so the discipline is to be **careful and
economical** — often the right move is to cut the loose word, not to expand it.

Terms of art stay and are good: a "hard cutoff", a "stiff" system, a "well-conditioned"
matrix are precise to the field. The target is the word used *colloquially where the field
would demand precision* — "hard", "easy", "natural", "simply", "obviously", "powerful",
"nicely", "robust" used as vibe rather than as a defined property. Flag each; fix the clear
ones by cutting or grounding; leave the ones that need a real metric explained for the
author, with a note on which metric would ground them.

## 3. Read it as the reader, not as yourself.

Set aside your own training. Take the reader the document is written for: a year of
calculus-based physics with Strang on the shelf, plus the instincts of a senior programmer
with a maths background — and nothing more. Read each passage as that person. Does it make
sense to them, first time through, without the knowledge you happen to have that they do
not? Where a step assumes something that reader was not given, name it. Where a term is
leaned on before it is grounded for that reader, flag it.

## 4. Would the sentence stand if the topic were beets, or car parts?

Grant the terms of art — assume the jargon is fine. Now ask whether the sentence *itself* is
well constructed underneath the vocabulary. Mentally replace "spherical tensor" with "beet"
and "Cholesky factorization" with "carburetor": is the sentence still grammatical, clear,
and well-built — subject and verb close and in agreement, each modifier attached to the
noun it means, no garden path, no clause stacked on clause until the point is lost, parallel
things in parallel form? Jargon hides bad construction; strip it in your head and the
scaffolding must stand on its own. Fix the broken scaffolding. The topic is allowed to be
science, but the sentence has to be a good sentence about anything.

## What to do

- **Fix the unambiguous:** factual errors; overclaims (weaken to exactly what the code
  proves); clearly-loose evaluative words that have an obvious precise replacement or a
  clean cut; broken sentence construction.
- **Flag the judgment calls:** a loose word that needs a real metric explained (say which
  metric); a reader-comprehension gap that needs the author's decision on how much to add;
  a sentence whose fix would change the meaning. Quote the sentence; name the problem.
- Keep every change true to the code and in the document's crisp register. Introduce no new
  claims. Do not pad — economy is part of the standard.
- Recompile with `pdflatex` until clean (no overfull boxes, no undefined references), with
  figures and listings intact.
- **Report**, grouped by the four questions: the fixes you made, and the flags for the
  author, each with the quoted sentence. Do not commit.
