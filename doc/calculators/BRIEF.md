# Brief: an expository document for one physics calculator

You're writing one self-contained document that explains a single calculator in
this NMR chemical-shielding pipeline: what it computes, the physics it descends
from, the rank-2 tensor it produces, and the applied math it actually runs. One
calculator per document. Yours is: **{CALCULATOR}** (source:
`src/{CALCULATOR}Result.{h,cpp}`).

Take your time and do a really good job. There's one of these per calculator and
they'll sit in the thesis as the reference a reader reaches for when they want to
know what a kernel means and where it comes from. There's no deadline on your end
and no reason to economize effort — a careful, re-read, re-edited document is the
whole point of the assignment.

## Who reads this

Someone who finished a year of calculus-based physics and math through linear
algebra and vector calculus, and keeps Gilbert Strang on the shelf — both
*Differential Equations and Linear Algebra* and *Linear Algebra and Learning from
Data*. They remember the material the way you remember a course you took two years
ago: the names are familiar, the book is within reach, the fluency is gone. So you
can throw them into real concepts — a change of basis, an eigendecomposition, a
Green's function, an irreducible representation, a gradient flow. The depth is not
the problem. What is not allowed is using a concept without explaining it, in clear
intuitive terms, as you go. When you reach for a tool, name it, say in a clause what
it does and why it's the right tool here, then use it. You are not re-teaching
linear algebra and you are not avoiding it — you give the intuitive handle as each
piece comes up, for a reader who can pick the rest up from Strang.

## The two ways this writing fails — guard against both

1. **Speaking down.** Assuming the reader is a novice, padding with reassurance and
   over-explanation, dropping to a simplified high-school picture. It reads as
   condescending and it buries the content. That register works for some material;
   it does not work here.
2. **The jargon party.** Stacking technical terms into prose that scans like a poem
   — rhythmic, dense, and opaque to the person who would actually use it. A
   practitioner reads it and cannot extract what you did.

The target sits between these: keep the technical language, ground every term of it.

## The standard for grounding

Do not eschew technical language — ground it. For every equation, say in plain,
*real* terms what each symbol means and what any specialized operator does. Plain
and real, not the fake-plain register where a sentence has the shape of an
explanation and conveys nothing. "The Laplacian, a fundamental operator, captures
the field's essential behavior" explains nothing. Real grounding says what the
thing measures: "∇² at a point sums how far the field there sits below or above
the average of its immediate neighbors — positive in a valley, negative on a peak."
If a symbol has units, give them. If a constant comes from somewhere, say where.
If an operator is specialized (a spherical-tensor decomposition, a switching taper,
a quadrature rule), say what it takes in and what it returns before you use it.

## Analogy as a precise model — reached for with confidence, verified by reflection

Turing's *The Chemical Basis of Morphogenesis* is the model here. It carries a
pure-mathematics treatment of reaction and diffusion, and alongside it passages
where an analogy is built that is a precise model of that mathematics — the picture
and the equations are the same structure said two ways, each part of the picture
corresponding to a term in the math. That is what good teaching can do: not
decoration laid over the math, but a second exact description of it the reader can
hold onto.

This is a tool to reach for, not a quota to fill. Most terms need only the plain
grounding described above. Deploy an analogy when you have one that maps exactly and
you can see that it clarifies — when there is confidence and clarity. A forced or
merely evocative analogy is worse than none: it leaves the reader holding a picture
that does not correspond to the thing. And do not reach for a stock teaching analogy
because it comes to hand — the worn logic puzzle, the river-crossing riddle, the
off-the-shelf metaphor borrowed from another field. Familiarity is not exactness; an
analogy earns its place only when it is drawn from the structure of this physics
itself.

The method this demands runs against the grain of how an AI writer works, and you
have to do it on purpose. You cannot write these explanations because they sound
good — sounding good is not evidence of anything. Build the analogy, then reflect it
back to yourself: re-derive it, set it beside the actual model, and test it on three
counts. Exactness — does each element of the analogy correspond to a real element of
the physics or math, with the relations among them preserved? Clarity — does it make
the thing easier to hold, or only easier to nod at? Truth — is it true, or only
evocative? Keep the analogies that survive all three; discard the rest, including the
ones you liked.

When the mapping is exact and you have it in hand, develop the analogy rather than
name it once and move on. The physics in these calculators often hands you one. The
angular factor (3cos²θ − 1)/r³ that recurs here is the field pattern of a magnetic
dipole — it is proportional to the second Legendre polynomial P₂(cosθ), zero at the
magic angle near 54.7°, and it is the same rank-2 object the T₂ part of the
decomposition carries. A document that writes "dipolar" once and moves on leaves that
intuition on the table; one that builds the picture — the lobes, the node, why a
source sitting at the magic angle contributes nothing — gives the reader the
geometric reason the formula has the shape it does. Reach for it where the
correspondence is that tight, run it through the exactness / clarity / truth check,
then use it. Where no exact analogy is at hand, the plain grounding stands on its
own; do not manufacture one to fill the space.

## The code is the final authority

The source — `src/{CALCULATOR}Result.{h,cpp}` and what it calls — is the final
authority, above every spec document, above the literature, above your own memory
of how this physics is usually written. The document describes what *this code*
does. Where the code and a spec doc disagree, the code wins. Where the code and the
textbook formula disagree, the code wins, and you note the difference rather than
quietly writing the textbook version.

## Read these to understand — but do not borrow their voice

Read the material below to understand the physics, the method, and the intent.
These are our internal documents: they are inward-facing and dense, written for us,
not for your reader. **Do not carry their vocabulary or their framing across.** Your
job is to bridge from this material to the reader cleanly — translate each in-house
term into grounded language — without dropping to a simplified high-school
explanation and without dissolving into jargon. When a project-specific term is
genuinely needed (e.g. "kernel"), introduce it in grounded terms the first time you
use it, then use it. Do not import in-house shorthand at all (calcset/`.LGS`, "the
ruler", "Environment Floor", "Pillar 1/2/3", "T2 angular residual"): say the
underlying thing in plain terms instead. Methodology terms like the wild-type-minus-
alanine difference, or the mutant pose pair, get grounded the first time, not parroted.

- `src/{CALCULATOR}Result.{h,cpp}` — the ground truth for the math. Read it first.
- `src/Types.h` ~lines 256–320 — `SphericalTensor`: the exact T0/T1/T2 definition,
  the isometric (real-spherical-harmonic) normalization that preserves the L2 norm,
  and the 9-double layout. This struct is the single owner of that decomposition.
  Describe T0/T1/T2 as the code defines them.
- `spec/PHYSICS_FOUNDATIONS.md` — §1 shielding fundamentals, plus the section for
  this mechanism: ring current §2, bond anisotropy / McConnell §3, Buckingham
  electric-field §4, aromatic quadrupole §5, H-bond §6, dispersion §7.
- `spec/APPLIED_MATHEMATICS.md` — §2 (closed-form spherical-tensor decomposition)
  and the entry for any other numerical operation this calculator uses.
- `spec/MATHS_GOALS.md` — the T0 corrections / T2 angular residual / classical-
  correctness material, the environment baseline the calculators must beat, and what
  the external tools must provide. (Read for intent; do not import its section names.)
- `data/calculator_params.toml` — the real parameter values. Cite the actual
  numbers with units. Do not invent constants or quote remembered ones.
- `spec/CONSTITUTION.md` — the project's claim-discipline conventions.

## What the document covers, in this order

1. **What it computes, and why it bears on shielding** (~½ page). The physical
   effect this calculator approximates, and how that effect changes the magnetic
   shielding a nucleus feels. Name the observable.
2. **Where shielding comes from** (~1 page). The minimal quantum mechanics: a
   nucleus in field B₀ feels B₀(1−σ); the electrons respond to B₀ and produce a
   secondary field at the nucleus; σ is the proportionality, and it is a 3×3 tensor
   because the secondary field's size and direction depend on how the molecule sits
   in B₀. Then connect that general picture to *this* calculator's specific
   mechanism.
3. **The tensor and its T0/T1/T2 split** (~1 page). Why σ is rank-2; what the 3×3
   matrix holds; the unique decomposition into T0 (isotropic, trace/3 — the part
   that survives molecular tumbling and sets the measured isotropic shift), T1
   (antisymmetric, 3 components, usually unobservable in standard NMR), and T2
   (symmetric traceless, 5 components — the orientation-dependent anisotropy). Use
   the code's normalization. State which of T0/T1/T2 this calculator populates (some
   emit pure T2).
4. **The method this calculator runs** (~1.5–2 pages). The actual model and its
   equations, taken from the source. Ground every term. Name the parameters and
   cite their real values from the TOML. Identify the specialized operations. State
   the approximations the model makes and where they stop holding.
5. **External inputs, and what the output is** (~½ page). Name the external data
   this calculator consumes — MOPAC charges, AIMNet2 charges/embedding, ff14SB
   charges, APBS field, DSSP — say the role each plays, and that its internal method
   is out of scope (one-line pointer, no explanation). Then state plainly: the
   output is a geometric quantity, not a shielding. Calibration against the DFT
   reference is what turns it into a shielding contribution; before that the number
   is geometric.

## Content rules specific to this project — these are correctness, not style

Honor each as a concept; do not import the in-house label for it.

- **A geometric quantity, not a shielding.** The calculator computes something
  derived from geometry. Never write that it "predicts the chemical shift." It
  produces a number that *correlates with* a DFT shielding contribution once
  calibrated.
- **Correlate, do not match.** Frame agreement with DFT as correlation strong
  enough to show the geometry carries the signal — never as the calculator
  reproducing DFT pointwise.
- **The anisotropy is the argument.** The rank-2 part (the 5 traceless-symmetric
  components) is the thesis. Do not collapse the tensor to a scalar, and do not
  imply the isotropic part is the whole story.
- **No simplification bias.** Do not average away element- or atom-type-dependent
  structure to get a tidy sentence. If an effect differs across atom types, that
  difference is content.
- **Do not assert physics from model fit.** R² and calibration quality are
  diagnostics that the descriptor set carries signal; they are not physical
  conclusions about the molecule.

## Method: write, then re-read and edit

Write the draft. Then read the whole thing again and edit specifically for:

- **Idea-smearing** — places where the prose blurs two ideas together or gestures
  at a concept instead of stating it. Sharpen to one idea per sentence, each stated
  plainly.
- **Clarity gaps** — any step the year-of-physics reader could not follow from what
  comes before it. Close the gap; do not paper over it with confident phrasing.
- **Superlatives and absolutes — remove them, all of them.** No "fundamental,"
  "powerful," "crucial," "essential," "the key," "the only," "canonical,"
  "simply," "just," "obviously," "clearly," "seamless," "robust," "elegant." If a
  claim is true, the plain statement of it is stronger than the inflated one.

## LaTeX and PDF

- One file. `article`, 11pt, `amsmath`/`amssymb`, `siunitx` for units, `bm` for
  vectors/tensors, `geometry` for sane margins. No exotic packages.
- Hard limit: five pages. Compile with `pdflatex`; check the page count and that
  there are no overfull boxes and no undefined references. Iterate until it
  compiles clean and fits.
- Number equations; define notation once; units through `siunitx`.
- Emit `{calculator}.tex` and the compiled `{calculator}.pdf`.

## Done when

- The arc runs QM origin → tensor → T0/T1/T2 → this calculator's method → external
  inputs → geometric-quantity-not-shielding.
- Every symbol in every equation is grounded in real terms, with units, and
  parameter values match `data/calculator_params.toml`.
- It bridges from our internal docs to the reader without borrowing in-house
  vocabulary, without high-school dilution, and without dissolving into jargon.
- The edit pass ran: no idea-smearing, no clarity gaps, no superlatives or absolutes.
- ≤ 5 pages, compiles clean. Source is the final authority throughout.
