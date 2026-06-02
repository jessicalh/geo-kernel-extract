# Brief: the verification pass — a restrained copy editor and fact-checker

You're reviewing one finished document. Someone else wrote it carefully, to the
standard in `doc/calculators/BRIEF.md`, and did good work. Your job is not to
rewrite it. It is to read it the way a trusted copy editor and fact-checker reads a
colleague's chapter: closely, with the source open beside you, fixing what is wrong
or off-standard with the lightest possible touch and flagging what you cannot settle
yourself. Restraint is the whole posture — a sentence that is correct and clear
stays, even if you would have written it differently.

Take your time. You have source access and there is no rush.

## Read first

- `doc/calculators/BRIEF.md` — the rules the document was written to. You are not
  re-deriving them; you are holding the document to them.
- the document under review (the `.tex` you are given).
- the calculator's source — the final authority. Where the document and the code
  disagree, the code is right and the document is wrong; correct it.

## Fact-check — this is the heart of the pass

The author read the source once; you read it again, independently, and cross-check.
Go claim by claim:

- Every equation: does it match what the code computes — the terms, the signs, the
  normalization?
- Every number: each cutoff, threshold, normalization factor, and parameter value
  stated in the text — find it in the source or `data/calculator_params.toml` and
  confirm it, digit for digit.
- Every description of behavior: filters, accumulations, order of operations, output
  field names and shapes — confirm each against the code.
- Units on every quantity.

Where a claim is right, leave it. Where it is wrong, correct it to match the code and
note the source location. Where you cannot find it in the source, do not assert it —
flag it as unverified. "The source doesn't show this" is a real and useful finding.

## Copy-edit — restrained, against the standard

Hold the prose to `BRIEF.md`: real-plain grounding rather than the fake-plain
register; an analogy only where it maps exactly and survives the
exactness/clarity/truth test (flag a forced or loosely-fitting one); no superlatives
or absolutes; no idea-smearing; no unexplained leap for the Strang-level reader; no
inward project jargon carried across. Fix clear lapses with minimal edits. Where it
is a judgment call, flag it rather than impose your taste.

## Restraint — what not to do

Do not paraphrase for flavor. Do not add flourish, transitions, or summary the author
chose to omit. Do not homogenize the voice. Do not expand scope. You are removing
error and friction, not leaving your fingerprints.

## After you edit

Recompile with `pdflatex`. Confirm the document still fits five pages with no overfull
boxes and no undefined references. If an edit pushed it over length or broke the
build, resolve that before you finish.

## Report

Close with a short, plain account: what you verified against source, what you
corrected (each with the source location), and your open flags — claims you could not
verify, and judgment calls you left for the author or for Jessica. No superlatives in
the report either; just what you found.
