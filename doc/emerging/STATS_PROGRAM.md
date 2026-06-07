# Stats Program — run under the guardrail

*Draft scaffold, 2026-06-06, for Jessica to shape and place. This is the discipline the
go-forward stats path runs under, and the program it runs. It exists so a future Claude/codex
instance does not walk in and start throwing good work out and narrowing left and right.*

---

## The guardrail (the governing discipline)

The one line, embedded, the north star of this whole path:

> **Say the sensible things at the right confidence, and leave the uncertain ones alone.**

Everything below is that line, unfolded.

- **The paper is not the code.** We did two years of programmer work in two months, and genuinely
  did it. That means some of what we built will have been silly — a botched T2 somewhere, a kernel
  resting on a bet that did not come in. That is fine, and it is not a debt. We do not talk about it
  in the paper. The paper is honestly based on the results, for **the part that turned out to pay off
  in truth.** Nothing in the repository obligates us to account for everything in it — only to be
  right about what we present.

- **Judge against *sensible*, not *proven*.** Sensible = decent applied maths, "see if this makes
  real sense" lab work, at the level of this degree. Not Nobel prizes. Far more clears the bar of
  *sensible* than clears the bar of *proven*, and the things that clear it are kept and presented
  carefully. Judging against *proven* is what turns the tool into a monster: against that bar almost
  everything fails, so it throws almost everything out. The bar is *sensible*.

- **You may leave a maybe as a maybe.** Footnote it, or do not mention it. You are not obligated to
  adjudicate it. The compulsion to resolve every uncertain thing — and resolving-down means deleting
  — is the engine of the narrowing. Resist it. Do not throw out the sensible in order to chase
  certainty about the uncertain.

- **A kernel that contributed once is presumptively real signal.** Before calling anything null, vet
  the pipeline — frame, units, self-bonded exclusions, target match — for a bug. Re-verify the
  current state before narrowing. Never compound on a single retraction. Never act on stale info (if
  a memory or note names a number, confirm it still holds before you build on it or narrow from it).

- **Do not load the dice — the analysis hopes nothing.** Jessica, who carries the accounting
  burden, may privately prefer that a speculative kernel fail: a failure is quietly dropped, a
  success must be explained. That is a human feeling and it is allowed — but it **never enters the
  analysis**, and the tool especially must hope for *nothing*, because a tool that hopes for an
  outcome will find it. Every kernel gets the same honest sweep and the same R². A speculative
  *success* does warrant scrutiny for a confound — but with the identical rigour owed to a *failure*
  (vet the pipeline for a bug before calling it null). Symmetric honesty: neither cheer nor jeer,
  neither cheaply keep nor cheaply kill.

---

## The level of scrutiny (the actual bar)

Calibrate to the real thing, not an imagined tribunal. Birkbeck does not even keep MSc theses — they
are deleted after a year. The people who judge them are serious — some are FRS — and the degree
stamps *University of London biologist* on your forehead, so it is real: it must be **heads-up**
(competent, correct, nothing a knowledgeable reader would wince at) and a good **demonstration of
mastery of the subject** — it is a Master's. But it is not expected to be Claude Shannon's thesis, or
anything near it, or even in the same field as one. The bar is **high on competence, correctness and
mastery; low — by design — on novelty and significance.** Right and heads-up; modest and unoriginal
on purpose.

## What we are ignorant of (both of us)

We are ignorant — the human and the tool both. The tool has access to the whole human corpus but
focuses on one thing at a time, which is its own kind of ignorance: narrow focus wearing the
confidence of broad access. Hold that; it is where the tool's wrong-certainty comes from.

A current, honest example: **we have not fully answered what our geometric kernel process encodes.**
We got closer, limned the algorithm, and established that it is a real and describable thing for
*some* kernels — the clean l=2 ones (EFG). For the *others* — the three-term McConnell, ring as the
far-field shadow of a Biot–Savart computation — it would be genuinely hard to say what we are really
doing. We do not know. The `GEOMETRIC_KERNEL_ACCOUNT` is therefore **complete for the solid kernels
and open for the funky ones**, and that is fine to say in exactly those words.

And, held open and *not* adjudicated: those same not-fully-understood kernels are the ones that
perform poorly in e3nn, and that might not be a coincidence — a kernel resting on assumptions that do
not cleanly hold could both resist accounting and fail to predict, for a single underlying reason.
This is precisely a *maybe*. We note it; we do not resolve it; and — per the dice rule above — we do
not let it tempt us into hoping the funky kernels fail so the puzzle goes away. They get the same
honest sweep as everything else, and whatever the sweep says, it says.

---

## What we know (the firm structure)

We know only a few things. We build on exactly those.

**1. The big model is a deliverable in its own right.** We want the best predictor we can get — all
the goodies, AIMNet2 and the works — per element and per IUPAC position, sidechains included. We do
not care that the neighbour's lost dog is in the stew: **if fido tastes good — if the model predicts
— the model is a deliverable**, and its ingredients need not be individually justified. Stand it up
first as an *interpolator* (within the data we hold), then test it as a *model* against a few fresh
extractions and short DFTs. "The model kinda gives us thirty solid percent" is a real win, said
plainly and not oversold.

**2. Kernels are not alike — classify them, up front.** They rest on different assumptions. Some are
so solid there is no question they go in (EFG, the established physics). Some are speculative bets —
the funky ones, resting on assumptions that may not hold (e.g. the three-term McConnell, which
carries T0+T1+T2 and is richer than a clean l=2; ring, whose clean form is only a far-field shadow
of a Biot–Savart computation). The sweep must say which kernel is which, by its assumptions, before
it reports a number.

**3. Ablate the funky ones off the top — then ablate everything.** The speculative kernels are
factored in and ablated *off the top*, prominently and first, so we can see whether they pay (and
hope they do not). In the end we need the **full ablation set** — every feature in and out — but the
funky ones lead, because their result is the one that decides scrutiny-or-drop.

**4. Every result comes back stratified, on the page, every time.** Per element alone, *and* per
IUPAC type including sidechains — a **full chart, here, every single time, even with small n.** Never
a pooled single number. Stage 1 told us **type is everything**: the pooled "nitrogen is hard" was an
artifact — backbone N is genuinely hard (R² ≈ 0.387), sidechain N is nearly the best (R² ≈ 0.887).
A pooled number hides exactly the structure that is the finding. We live with small n; this is a
little exploratory thesis, and small n honestly reported beats a pooled number that lies by
averaging. The chart is not a summary of the result — the chart *is* the result.

**5. The sweep procedure.** Start with **complete PySR** for anything we might see. Then **try
equivariance** where it makes sense. **Always get an R².** Run it all before judging any of it; the
sweep *generates* the landscape, and only afterward do we *select* the sensible and footnote or drop
the rest. The sweep does not pre-narrow ("skip that, it won't work") — pre-narrowing is the monster
deciding before the data does.

---

## How the two halves fit

The sweep (§5) generates; the guardrail selects. The sweep is exhaustive on purpose — it runs every
type, every ablation, always an R², and pre-judges nothing. The guardrail is calibrated on purpose —
afterward it keeps the sensible, footnotes or omits the maybe, and never destroys a contributor to
chase certainty. Neither half narrows. The big model (§1) is the deliverable that does not need its
stew explained; the stratified charts (§4) are the honest, type-resolved record that does. Between
the two, we say the sensible things at the right confidence, and leave the uncertain ones alone.
