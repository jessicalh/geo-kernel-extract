# Predicting NMR chemical shielding from protein geometry — methods & progress
### A short précis (draft, 2026-06-04)

> **Draft status (trued 2026-06-04).** Stage 1 numbers have been refreshed from
> the `learn/` material, the Stage 3 section has been rewritten around the
> two-model e3nn pipeline, and the 1P9J interpolation v0 result is built:
> T2 component R2 ≈ 0.48, |T2| r ≈ 0.82, T0 modulation R2 ≈ 0.45. Keep the
> 1P9J-within caveat: no transferability, between-axis, BMRB validation, or
> process-model claim from this figure. Markdown now for fast iteration → LaTeX/PDF for the desk once content locks.
> This is a draft-to-best scaffold; the final words are the author's.

## The question

Can the **NMR chemical-shielding tensor** of a nucleus in a protein be accounted for, from
**protein geometry alone**, by a set of **classical physics kernels** — ring current, electric-field
gradient from charges, bond/group anisotropy, hydrogen-bonding — calibrated against quantum (DFT)
shielding? The aim is **physics explanation**, not a black-box predictor: we want to know *which*
geometric mechanisms carry the signal, and *how much*.

The project takes the advisor's own framing seriously — *a model that shows a result need not be a
good model* — and pivoted early from "machine learning on docked poses + an empirical predictor"
to **geometric kernels calibrated against r²SCAN DFT**. Three stages fell out of the coupling between
*extracting* the full rank-2 shielding tensor and *validating* it:

1. **Stage 1 — mutations** (settled): is the kernel set rich enough to carry the signal across many proteins?
2. **Stage 2 — the laws** (in progress): what *equation* sits behind the kernels, and does it recover DFT?
3. **Stage 3 — the learning model** (forward arc): turn the explanation into a transferable *predictor*.

A constant throughout: we keep the **full rank-2 tensor** (the `T0/T1/T2` spherical decomposition), never
collapsing it to a scalar. The angular structure — the `T2` part — *is* the argument; averaging it away
would discard the physics we are trying to explain.

## Stage 1 — mutations (settled)

*[Numbers recovered from the `learn/` archive 2026-06-04 and cross-checked against current work. The
0.818-vs-0.718 reconciliation below corrects a conflation in our running notes (incl. CLAUDE.md) — worth
confirming before it goes in the thesis.]*

The dataset is **723 BMRB single-chain proteins** (720 complete), each prepared as a **WT pose and an
all-aromatics→ALA mechanical mutant** (backbone fixed); the target is the **WT−ALA static DFT
shielding-delta** per matched atom, with the mutated residue itself excluded — so the target is the
*effect on the surrounding protein*. On this, **per-element, per-atom-type ridge regression** over the
**126-kernel T2 layout** ("55 kernels" was a core-subset label) reaches **R² = 0.818 on the 110-protein
fair set** (69K atoms) and **holds at 0.718 on the full 720 proteins / 446K atoms** — the drop is
cross-protein generalisation, not a physics failure.

The headline lesson runs *against* simplification. "Nitrogen is hard" is an **element-pooling artifact**:
**backbone N is genuinely hard (R² = 0.39)** while **sidechain N is the second-easiest atom type (0.89)**,
behind only H (0.93) — and *every* element is easier in the sidechain than the backbone (e.g. C=O at 0.46,
the hardest carbon, vs sidechain C at 0.73). The **effective-dimension count differs by element — H ≈ 20,
C ≈ 6, N ≈ 3, O ≈ 12** — and that spread *is* the story, not noise to average. R² here is a **diagnostic**
that the kernels carry the signal, never the thing we optimize.

Crucially, Stage 1 already measured the gap that defines Stage 2: a per-protein ceiling of **~0.81 within
a protein but only ~0.35 across proteins**. That within/between split is not an artifact of one protein —
it is the central quantitative motivation for everything that follows.

## Stage 2 — the larger equation, and what we are fitting

### One classical object behind the calculators

The through-space part of the shielding tensor is a **local sum over sources** of a single dipolar kernel:

$$\sigma_{ab}(i)\;\approx\;\sum_{\text{sources }s} I_s \, D_{ab}(\mathbf{r}_{is}),
\qquad
D_{ab}(\boldsymbol{\rho}) \;=\; \partial_a\partial_b\frac{1}{\rho} \;=\; \frac{3\,\hat\rho_a\hat\rho_b - \delta_{ab}}{\rho^{3}} .$$

`D_ab` is symmetric, **traceless**, with eigenvalues `(+2,−1,−1)/ρ³` and a **magic-angle node at
cos²θ = 1/3** — the familiar `(3cos²θ − 1)/r³` angular law. The claim that organizes Stage 2 is that the
**classical calculators are projections of this one object** onto different source axes and source lists:

| calculator | what it sums | the projection |
|---|---|---|
| charge / EFG | point charges `q_j` | `Σ_j q_j · D_ab(r_ij)` |
| McConnell (bond anisotropy) | bond susceptibility `Δχ` about the bond axis `ĥ` | `Σ Δχ (3ĥĥ − δ)/r³` |
| H-bond (geometric) | donor→acceptor dipolar axis | same dipolar form |
| ring current | aromatic ring loop (Biot–Savart / Johnson–Bovey) | a *related* rank-≤1 object; Pople `(3cos²θ−1)/r³` is its far-field limit |

So the "larger equation" is **shared angular form, per-mechanism radial intensity** — and the
intensity `I_s` (the coefficient that turns a geometric shadow into ppm) **is the deliverable**: the
*calibration* is a fitted law plus its coefficient, not a good R².

### How we grade — statistical position, not literature-match

We do not ask "does the coefficient equal the textbook value." We ask **where our desired (literature)
form falls in the distribution of fits the data admits** — its *statistical position* — and **whether we
can determine what drives a recovery**. A result is reportable when it is (1) fairly indicative and not an
artifact, and (2) robust in context — with probability framing and honest caveats, never a binary
pass/fail. *Correlate, do not match.*

### Progress on 1P9J (the within instrument)

On a single protein with a deep trajectory (1P9J, ~15 ns, 751 DFT frames), the **within-geometry axis**
carries the statistics. There:

- **charge q/r³** and **ring current** both **recover** (the clean "cookies");
- the **unified `D_ab`-sum combine recovers the through-space tensor** where the individual shadows cannot
  isolate — and it is a **real combine, not charge in disguise** (it is carried by the MOPAC-field and
  McConnell projections, with charge contributing ≈0). If it holds, this is the **deepest** result: the
  *one object* reconstructs what no single shadow does alone.

**Honest framing.** One protein gives a *within* instrument; its *between-atom / between-protein*
(transferability) axis is thin by construction. So 1P9J grades the within axis; the between axis is
deferred — by design — to Stage 2's transferability pilot.

### Anticipated Stage 2 results

*[To be sharpened from `SPEC_720_STATIC_INGEST` + `SPEC_MATHS_FIX_AND_REVIEW` as they land.]*

The **720-WT statics pilot** (the same r²SCAN DFT we already have — *no new quantum chemistry*) is the
between-axis instrument: many proteins, many rings, absolute shieldings. We expect it to (a) earn
**cross-protein validation** of the charge law (its strongest, static, axis), (b) **fatten ring current's
thin between-axis**, and (c) let us state the per-mechanism coefficients **in physical units** with
population-level confidence — turning today's within-axis cookies into a defensible **equations table**.

## Stage 3 — the learning model (the forward arc)

The forward step flips the ethos from **explanation to prediction** — here **R² *is* the metric**. It is
**two GPU models**, not one, and not a staged crawl toward one.

**Model 1 — the shielding emulator.** A per-atom **equivariant GNN (e3nn)**, applied **locally, per frame**,
predicting the full shielding **tensor** (T0 + T2; the antisymmetric T1 only if we want it — shielding is
non-symmetric) from local geometry plus the per-atom feature pile (kernel shadows, AIMNet2, MOPAC, APBS,
dihedrals/motion, exposure). It trains **only on DFT** — 1P9J's many frames *plus* the 720 static poses — so
the physics stays clean; and once trained it predicts shielding **where there is no DFT.** Equivariant-GNN is
the established route to full shielding tensors (an e3nn model on ²⁹Si tensors cut error ~53% over an
invariant baseline — the empirical case for equivariance on exactly this problem). The **readout irreps are
the crux**, and the **binding constraint is the DFT corpus's size and diversity, not the architecture** —
which is exactly why the 720 *is* the corpus.

**Model 2 — the shift predictor.** A second GPU model on the **no-DFT fleet** (~656 proteins) that eats the
**full per-frame firehose** — Model 1's predicted shieldings *plus* the feature outputs across the trajectory
— and is **trained to predict the measured shift** against **BMRB**. The leap from a per-frame quantum tensor
to a solution-averaged observable (the ensemble/reference mapping) is **learned, not hand-averaged.** Because
the fleet has no DFT, Model 1's prediction fills the gap — and we **ablate Model 1 in and out to test whether
the distilled DFT earns its place.**

The two-model split is the elegant part: the **clean DFT trains the physics (Model 1)** while the **noisy
measured shift trains the observable (Model 2)** — so the real, messy experiment never corrupts the physics
model. The arc is **distill → transfer → test**, and **1P9J is the calibration case** because it carries
*both* grounds — DFT to check Model 1 directly, BMRB to check Model 2 — before the same machinery rolls out
over the 656 where only BMRB exists. **Un-phased:** two real models, **make Model 1 first** — it de-risks the
effort and scopes Model 2 exactly.

**The running prototype.** Model 1's first instance already exists — a per-frame equivariant predictor that
recovers 1P9J's held-out shielding tensor (T2 R² ≈ 0.48, |T2| r ≈ 0.82): a software tool that *does*
equivariant prediction on a protein. *[Figure: the capability graph.]*

## The arc, in one line

**signal → equation calibrations (each law with its coefficient + confidence) → ensemble → equivariant
transferability pilot.** Stage 1 showed the signal is there; Stage 2 is writing the equations and earning
their confidence; Stage 3 turns them into a transferable predictor.
