# Brief: the numerical-methods addendum for one calculator

You're writing a short addendum to a calculator's exposition document --- a separate
document, same visual format --- that takes the *fancy maths* in the calculator's code
and gives it a code-review treatment, grounded in the source: the factorization, the
quadrature, the closed-form decomposition, the finite-difference stencil, the
constrained solve --- whatever in this calculator would earn a section in a
numerical-methods book.

The worked exemplar of this standard is `doc/calculators/haighmallion-numerics.tex`. It
took many rounds to reach the register below; read it before you start and treat it as
what "good" looks like for this series. The quotations in this brief are from it.

Take your time. The bar is high and the document is short; both are the point.

## Who reads this

The exposition document's reader (a year of calculus-based physics, Strang on the
shelf) with, added, the instincts of a senior programmer who has a maths background and
an applied brain. They know SVD and Cholesky and quadrature exist; they do not need the
textbook derivation. What they want is how *this code's* numbers actually move through
the technique. Write at the clarity of K&R, Jon Bentley, and Steve McConnell.

## What this addendum is for --- a technique inventory

This is an inventory of the numerical techniques this calculator runs, not a spotlight
on the one most impressive. Cover them --- every genuinely numerical technique the
calculator uses, including the constructions it *builds up and then consumes* (the
geometry it assembles, the triangulation or grid it forms, the neighbour structure it
queries), not only the headline method. Skip the truly pedestrian --- a dot product, a
norm, a running sum need no treatment --- but do not skip a real technique because it is
supporting rather than central. Convey relative importance and keep it clear as you go:
open with a sentence or two mapping the techniques the document covers and their weight,
then take them in a sensible order. `spec/APPLIED_MATHEMATICS.md` inventories the
numerical methods the system uses; cross it with this calculator's source.

Leaving out is the craft. Cut the pedestrian, the redundant, the line a reader could
reconstruct --- but what remains must still be a coherent image of how the calculator
works, complete enough to enable full understanding. The test of a cut is whether the
picture stands whole without it: a removal that opens a gap in the through-line was
load-bearing, and stays. Brevity that costs the coherent whole is not brevity, it is
damage.

## Numerical Recipes is the touchstone --- for register, not for stance

Numerical Recipes (Press, Teukolsky, Vetterling, Flannery) is the touchstone for the
*register*: the implementer's-eye explanation, with real code, that a working programmer
could follow to see how a method actually works. Match that clarity and that
concreteness.

The stance, though, is different. Numerical Recipes gives advice --- it recommends
methods, tells you what to prefer and what to avoid. This addendum gives none. The job
is to make the choices in *this* code transparent: which method it uses, how it works on
these numbers, what a guard prevents, when an approximation stops holding, how the cost
scales --- every one stated as a fact about the code, not as guidance to a reader writing
their own. "These are hard omissions --- the skipped terms leave the sum, and nothing is
reweighted to make up for them." "The adaptive part is a heuristic, not an error
estimate." Describe the choice and its consequence; do not recommend, do not rank
alternatives, do not address the reader's future code. No "you should," no "prefer," no
"for better stability."

## The register: composed, not breezy

Direct and plain, but composed --- the K&R / Bentley / McConnell clarity nudged a touch
formal, never breezy.

- **"We," not "you."** Take the reader through the problem as a shared look at it, not
  as instructions: "We have five to nine atoms lying close to a common plane, and we
  want the normal of that plane."
- **Variables already have a programmer's attention.** Ground a code identifier the
  first time it appears --- tie it to the symbol it carries --- but do not announce it.
  "`sep` is \(\bm\rho = \bm r - \bm s\), the vector from the surface point to the field
  point" earns its keep; "Here `cfg` is the `CalculatorConfig` accessor" does not,
  because a programmer reading `cfg("…")` already knows what it is.
- **Every descriptor is earned by the topic; the fizz must be in the drink.** The
  liveliness comes from the content --- the least-variance normal, the spike at the
  pole, the double-count behind a \(\sqrt2\) --- not from words telling the reader the
  content is lively. Cut the intensifiers and the performance: "genuinely," "cleverly,"
  "simply," "just," "the one honest thing," "the whole craft," "the picture is simple,"
  "worth a second look," "the payoff is concrete." Nothing existential, nothing
  AI-flavored.
- **Keep the metaphors that carry the maths; cut the ones that only decorate.** Earned:
  "the \(\sqrt2\) on each off-diagonal pays for double counting: \(S_{xy}\) sits in the
  matrix twice, as \(S_{xy}\) and \(S_{yx}\), so it owes \(2S_{xy}^2\) to the Frobenius
  sum"; "the field shape of a point magnet." The test of a cut is the coherent image:
  does the picture still stand whole without it?

## The standard --- tell how the numbers move

This is the whole job, so it gets the most words.

- **For each technique: one or two plainspoken paragraphs and a clearly explained code
  block.** The paragraphs are a code review, not a derivation --- what the technique
  does to the numbers here, why it is the right tool for this calculator's shape of
  problem, and where it would break. The code block is the illuminating lines from the
  source (distil if the real code is cluttered, and say that you have), with each step's
  effect on the numbers stated.
- **Concrete and grounded.** Name the actual shapes and sizes from this code --- the EEQ
  system is \((N{+}1)\times(N{+}1)\) for \(N\) atoms; the triangle quadrature is a
  7-point rule, at most two subdivisions, so at most 112 points. When you write a
  symbol, say which variable or array in the code holds it. The relationship between the
  technique and the notation must be concrete and traceable, never gestured at.
- **No "this is an X, see reference Y."** Do not name the algorithm and defer to a book.
  The reader can find the reference; your job is to make *this code's* use of it legible
  --- how the input becomes the output, step by step. Naming-and-deferring is the
  failure mode this addendum exists to avoid.
- **Tell how the numbers move; do not handwave the formula.** "We solve \(Ax=b\) by
  Cholesky" is a handwave. The treatment is: \(A\) is symmetric positive-definite
  because [reason from the code]; we factor \(A = LL^{\mathsf T}\) with \(L\)
  lower-triangular; forward- then back-substitution recovers \(x\); here is the call and
  what it returns; for this calculator \(A\) is [shape], so [cost / conditioning note].
  That is the level.
- **Ground every variable, in the maths and in the code alike.** When a name first
  appears, whether \(\rho\) or `sep`, tie the code identifier to the quantity it carries
  and say what that quantity means. A reader should never meet a variable --- in a
  formula or in a listing --- that the prose has not named.
- **Do not list rules and let the meaning be assumed.** A block of tensor components, a
  set of quadrature points, a row of normalization factors --- for each, say what they
  are, what they mean in plain language, and what they mean numerically. The five
  \(T_2\) factors are the worked example: they carry the orientation-dependent content
  of the tensor, and the \(\sqrt2\) on each off-diagonal is there because a symmetric
  matrix counts every off-diagonal twice in its Frobenius norm.
- **State exactly what the code proves, never a notch more.** The boundary cases are
  where overclaim hides. Distinguish what the code *fits* from what it *uses*: the SVD
  fits the ring *normal* (the axis the result is later projected onto), while the
  integration *surface* is the actual centre-to-vertex fan --- "piecewise planar through
  the buckled positions, never flattened onto the fitted plane." Honour the boundaries
  the code holds: "rank at most one," not "rank one"; "not greater than the ring
  radius," not "under it"; "beyond degree five it is not guaranteed," not "only
  approximates degree six"; the isotropic *shielding*, not the *shift*. A claim one
  notch stronger than the code is a bug.

## The shared spherical-tensor decomposition

Every calculator decomposes a \(3\times3\) tensor into T0/T1/T2 through the same closed
form (`SphericalTensor::Decompose`, `src/Types.h`). It is itself a clean numerical
example --- closed-form, with the isometric (real-spherical-harmonic) normalization
chosen so the five-vector's squared length equals the Frobenius norm of the symmetric
traceless part. Give it a compact treatment (the closed form, the normalization
constants, which components this calculator populates); you need not re-derive its full
provenance in every addendum. Spend the words on the technique distinctive to this
calculator.

## The glossary --- three beats

Collect the recurring terms in a short glossary at the end. No row of rules is
self-explanatory; give each term three beats:

1. an illustrative first sentence that makes the reader *see* the thing;
2. a concrete example that earns the middle slot --- a real instance, and the best ones
   deepen the entry;
3. a formal sentence rigorous enough to quote.

The dipolar-kernel entry is the shape to match:

> **Dipolar kernel.** The field shape of a small bar magnet --- strongest and signed
> along its axis, weaker and reversed to the side, fading as the cube of distance. Along
> the axis it is \(+2/\rho^3\), broadside \(-1/\rho^3\), and zero at the magic angle
> (\(54.7^\circ\)) where the two balance. Formally, it is the Hessian of \(1/\rho\), the
> symmetric traceless tensor \((3\rho_a\rho_b - \delta_{ab}\rho^2)/\rho^5\).

Keep it to the terms the document leans on.

## Figures --- geometry deserves a picture

When the method is geometric --- a fan, a mesh, a finite-difference stencil, a sampling
pattern --- reach for a small figure; it lands in a glance what prose labours at.
Mermaid is the wrong instrument (it draws flowcharts and sequence diagrams, not
geometry); use TikZ. Keep it to what a picture conveys better than words --- a dataflow
box-chart of a linear pipeline does not earn its place. The figure is part of the
verification surface: it must match the code, and the review will check it against the
source.

## The code is the final authority

The calculator's source is the truth. Describe what the code computes --- the actual
stencil, the actual factorization call, the actual quadrature constants --- not the
idealized algorithm. Where the code departs from the textbook version (a guard, a clamp,
a fixed iteration count, a non-obvious ordering), that departure is exactly what is worth
explaining.

## Read these

- the calculator's source (the authority).
- `spec/APPLIED_MATHEMATICS.md` --- the inventory of numerical methods used and where.
- `src/Types.h` (SphericalTensor) for the shared decomposition.
- the calculator's exposition document (`doc/calculators/<slug>.tex`) --- so this
  addendum complements it and does not repeat the physics.
- `doc/calculators/haighmallion-numerics.tex` --- the worked exemplar of this standard.

## Format

- `\usepackage{nmrdoc}` gives the shared format and the `codebox` environment for code:
  a pale-yellow box that breaks gracefully across a page. Put each listing in a
  `codebox`; load `siunitx`, `bm`, and `tikz` (for a figure) as the document needs.
- Show the illuminating lines, not the whole function, and label the distillation: "The
  listings below are distilled from the source --- loops and bookkeeping trimmed --- but
  the constants, guards, and arithmetic are the implementation's."
- Hand-fit the lines so nothing wraps mid-expression; the box does not auto-wrap. Distil
  a too-long line (a short local alias, said) rather than letting it break.
- Keep the document short --- brevity is part of the spec. Write to
  `doc/calculators/<slug>-numerics.tex`, title via `\NmrTitle{...}`, compile clean.

## Done when

- Each technique has a grounded, how-the-numbers-move treatment with a code block a
  senior programmer with maths could follow; no name-drop-and-defer, no handwave.
- The register is composed, not breezy: "we" not "you", every descriptor earned, no
  intensifiers, no advice --- transparency only. Every claim is exactly what the code
  proves, no notch stronger.
- Lists carry plain and numerical meaning; the glossary's entries run picture / example
  / formal; a geometric method gets its figure.
- The shapes and notation are the code's actual ones; it is short; it compiles in the
  shared format with `codebox` listings.
