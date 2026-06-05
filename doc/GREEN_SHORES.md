# Green Shores
### The solid ground — the science and the math, described as they are

Before we push off into harder water — the outline-and-progress document,
telegraphic and tabular, that has to say *what we are doing and how and why* — this
is a last clear sight of the coast we already know. It is the enjoyable country: the
science and the math set down as they actually are. Nothing here is a plan or a
promise. It is description, finished and in good order, and we will see it again,
because everything the thesis argues sails from this harbour.

## The one idea

Every magnetic and electrostatic calculator in the library is a classical shadow of
a single quantum object — the Ramsey current response, the way a molecule's electrons
circulate when a field is switched on. Reduce that object to a point and a direction
and one kernel keeps falling out:

```
D_ab = ∂_a ∂_b (1/ρ) = (3 ρ̂_a ρ̂_b − δ_ab) / ρ³
```

the second derivative of a Coulomb potential. Its eigenvalues are `+2/ρ³` along the
separation and `−1/ρ³` across it (twice). Contract it onto an axis and the angular
factor `(3cos²θ − 1)/ρ³` appears — the McConnell scalar, the dipolar coupling, the
electric-field gradient, all the same shape wearing different clothes. The rank-2
tensor it produces splits, as any 3×3 does, into **1 + 3 + 5**: the trace **T0**, the
antisymmetric **T1**, and the five-component **T2** that carries the angular story the
thesis is about. That decomposition is preserved end to end. We never collapse it to a
number.

## What we have described

**The three core chapters** — the synthesis, each truth-vetted and committed, each
with a `-backup.tex`:

- `physics-architecture.tex` (`f785fe4`) — the unified view: every calculator as a
  reduction of the one current response, organised around the kernel `D_ab`.
- `categorical-grounding.tex` (`6ca47d6`) — the categorical view: by atom type and
  geometry, why each measure exists, which tool reaches for which nucleus, and the
  determinant literature behind each.
- `method-and-proxy.tex` (`5c4bbf6`) — the intellectual shape: how each method
  *proxies* the physics, from the literal to the metaphorical; the geometric kernel as
  the question put to nature, the calibration as nature's answer.

**The twelve exposition docs** — one per calculator: QM origin → T0/T1/T2 → method →
"this is a geometric kernel, not a shielding." Source-verified.

**The twelve numerics addenda** (`2500bfa`) — the same twelve given a
Numerical-Recipes-grade, code-grounded reading of the actual arithmetic. HaighMallion
is the locked exemplar; the other eleven were calibrated to it.

**The operator guide and the architecture note** — how the extractor is run, and how
the `Protein` / `ProteinConformation` / `ConformationResult` object model holds
together.

A descriptive census of the twelve — what each reaches for, and its one physical idea.
This table *describes*; it does not yet *weigh*:

| Calculator | Reaches for | The one physical idea |
|---|---|---|
| McConnell | ring-current shifts | point magnetic dipole at the ring centre; `(3cos²θ−1)/ρ³` |
| RingSusceptibility | aromatic rings | anisotropic susceptibility tensor of the ring loop |
| BiotSavart | aromatic rings | the loop current's field, by direct integration |
| HaighMallion | aromatic rings | geometric ring-current factor summed over the ring bonds |
| APBSField | charged & polar sites | Poisson–Boltzmann reaction field of the solvated protein |
| WaterField | polar sites | electrostatic field of the explicit water |
| EEQ | all charges | electronegativity-equalisation partial charges |
| PiQuadrupole | aromatic / π systems | the electric quadrupole of the π cloud |
| HBond | H-bonded H/N/O | geometric hydrogen-bond field |
| LarsenHBond | amide protons | Larsen's calibrated H-bond shift (in ppm — the one exception) |
| Dispersion | all atoms | the London/dispersion contribution |
| HydrationShell | surface atoms | the geometry of the solvent shell around each atom |

## Where the shore ends

This coast is **description**. It says what each tool computes and why the physics
says so. It does not say whether the tools, taken together, predict a shielding better
than nothing — that is the model-utility question, and it is genuinely open. Stage 1
has already shown the math carries signal: the per-element, per-atom-type ridge reaches
R² = 0.818 on the 110-protein fair set (69K atoms) and holds at 0.718 across the full
720 proteins / 446K atoms, the drop being cross-protein generalisation rather than a
physics failure, over a 126-kernel T2 layout; and Stage 2's 1P9J campaign is in frame.
But the *argument* that ties captured signal to utility — and the
terse progress document that tracks what we are doing, how, and why — is the sea we are
about to sail. We named the coast so we can find it again when the fog comes in.

---

*The shore is the science and math as they are; the sea is the outline-and-progress
document, built together. This file is the waypoint between them.*
