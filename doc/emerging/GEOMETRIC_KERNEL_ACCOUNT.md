# The Geometric Kernel, Literally — What the Compression Is and How It Is Calculated

*A working account, written in full, every term defined. This is the substance we
need in order to **account for** the geometric kernels in the thesis — not to
**claim** anything about them, but to describe, truly and literally, what they are
and how they are computed. An MSc in methods for structural biology does not owe a
mechanistic explanation of why a result holds; it owes a responsible, non-hand-waving
account of what each piece is and what it does. This document is that account for the
geometric kernels. It is written without digest, on purpose: nothing here is
telegraphed, and every symbol is defined where it first appears.*

---

## 0. What this document is, and is not

It is an **account**: a literal statement of what the geometric kernels are, as
mathematical objects, and what operation each one performs on a structural frame.
It is **not a claim** that they explain shielding, nor that the rank-2 signal is
physically load-bearing for any particular atom. Those are empirical questions held
elsewhere; here we only say, exactly, what the thing is. In the current forward
build these kernels are kept and emitted, not used as Step-1 inputs.

An honest note on origin, so the account does not have to invent a principled one:
the geometric kernels are in this project because they were offered early and kept on
a hunch — *"keep that, who knows."* That contingency does not weaken anything below,
because the account rests entirely on **what the kernels are**, which does not care
how they came to be here. A true description of an object owes nothing to the story
of its arrival.

Two phrases recur and are worth fixing now. **"Claim"** = an assertion we would have
to defend as true of our data. **"Account"** = a true, situated description of what a
component is and does, which stands on its own correctness rather than on a contest
with a skeptic. The whole document is in the second register.

---

## 1. Definitions and notation

We define everything used below, in the order it is needed.

**Atoms and positions.** We study one *target atom*, written **i**, and the atoms
around it, the *source atoms*, each written **j**. Each atom has a position in
three-dimensional space: **r_i** is the position vector of the target, **r_j** the
position vector of a source. Positions are taken from a single structural *frame* — a
single conformation, e.g. one snapshot from a molecular-dynamics (MD) trajectory.

**Separation.** The *separation vector* from source j to target i is

  **ρ** = r_i − r_j .

Its *length* (the distance between the two atoms) is the non-negative number

  **ρ** = |ρ| = √(ρ·ρ) ,

and its *direction* is the *unit vector*

  **ρ̂** = ρ / ρ , so that |ρ̂| = 1 .

We write the Cartesian components of these vectors with a lower index taken from
{x, y, z}: **ρ_a** is the a-th component of ρ, and **ρ̂_a** is the a-th component of
ρ̂. The symbol **a** (or **b**) is therefore an index ranging over the three spatial
axes.

**Neighbourhood and cutoff.** Not all atoms in the structure are summed over; only
those within a chosen distance of the target. We write **N(i)** for this
*neighbourhood* — the set of source atoms j retained for target i — and the distance
bound that defines it is the *cutoff*. The cutoff is a chosen number; that choice is
itself part of the account (it sets what counts as "context"), but it is not part of
the kernel's mathematical form.

**Source weight.** Each source j enters the calculation through a single number,
its *weight*, written **w_j**. What this number physically *is* differs by calculator:
it may be an electric charge, an induced magnetic moment, a ring-current intensity
factor. The point of the definition is that the geometry of the neighbourhood enters
through ρ_ij and the chemistry of each source enters *only* through the scalar w_j —
nothing else about source j survives into the kernel.

**Kronecker delta.** **δ_ab** is the *Kronecker delta*: the number 1 when a = b and 0
when a ≠ b. As a 3×3 array it is the identity matrix. It is what "subtract the
isotropic part" looks like in index notation.

**Partial derivative.** **∂_a** means ∂/∂x_a — the rate of change of a quantity with
respect to the a-th spatial coordinate of the target. Applying it twice, ∂_a∂_b,
gives a *second* derivative, one of the nine ways to differentiate twice with respect
to the three coordinates; these nine numbers form a 3×3 array indexed by (a, b).

**Multipole order (degree).** The integer **l** = 0, 1, 2, 3, … labels how an angular
quantity behaves under rotation. l = 0 is a *scalar* (unchanged by rotation), l = 1 is
a *vector* (three components, rotating as an arrow), l = 2 is a *rank-2 symmetric
traceless tensor* (five independent components, the quadrupolar pattern), and so on.
The number of independent components at order l is **2l + 1**: one for l = 0, three
for l = 1, five for l = 2. The *real spherical harmonics* **Y_l^m(ρ̂)** are the
standard set of angular functions of a direction at each order l (the index m runs
over the 2l + 1 components). They are the canonical "angular shapes" out of which any
direction-dependent quantity can be built.

**The magic angle.** The l = 2 angular shape (3cos²θ − 1) is zero when
3cos²θ − 1 = 0, i.e. cos²θ = 1/3, i.e. **θ ≈ 54.7°**. This is the *magic angle*; here
θ is the angle between the separation direction ρ̂ and a chosen reference axis (for a
ring current, the ring's normal). The l = 2 pattern is positive along the axis,
negative in the perpendicular plane, and crosses zero at the magic angle.

**Tensor vocabulary.** A *rank-2 tensor* is a 3×3 array of numbers M_ab. Its *trace*
is the sum of its diagonal, Σ_a M_aa. It is *traceless* if that sum is zero. It is
*symmetric* if M_ab = M_ba. A symmetric traceless 3×3 array has exactly **five**
independent numbers, which is why it carries the l = 2 (five-component) content.

**Spherical-tensor decomposition (T0 / T1 / T2).** Any 3×3 array can be split, with no
loss, into three pieces that do not mix under rotation: the *isotropic* part **T0**
(the trace; one number; l = 0), the *antisymmetric* part **T1** (three numbers; l = 1),
and the *symmetric traceless* part **T2** (five numbers; l = 2). When we say a kernel
"preserves T2," we mean its five-component l = 2 piece is kept rather than averaged
away into the single number T0. T2 *is* the anisotropy.

**Chemical shielding.** The quantity we are ultimately about is the *nuclear magnetic
shielding* σ at an atom — a rank-2 tensor describing how the atom's electrons screen
it from an applied magnetic field. By Ramsey's theory it splits into a *diamagnetic*
part **σ_dia** (a relatively simple ground-state electron-density term) and a
*paramagnetic* part **σ_para** (a harder term involving a sum over electronic
excited states, large for atoms with low-lying excitations and valence p-electron
angular momentum). The *isotropic* shielding is (1/3) of the trace of σ; the
*anisotropy* is what is left, the T2 part. This vocabulary matters in §6.

---

## 2. The literal object: a multipole field of the weighted neighbourhood

Start from the *Newtonian / Coulomb potential* of a unit source, the function 1/ρ of
the separation length. Every geometric kernel in this project is a spatial derivative
of this one function, summed over the weighted neighbourhood.

**First derivative.** Differentiating 1/ρ once with respect to the target coordinate
x_a:

  ∂_a (1/ρ) = − ρ_a / ρ³ = − ρ̂_a / ρ² .

This is (minus) the *field* of the source: it points along the separation and falls
off as 1/ρ². It is an l = 1 (vector) quantity.

**Second derivative.** Differentiating once more, with respect to x_b:

  ∂_b ∂_a (1/ρ) = ∂_b ( − ρ_a / ρ³ )
        = − δ_ab / ρ³ + 3 ρ_a ρ_b / ρ⁵
        = ( 3 ρ̂_a ρ̂_b − δ_ab ) / ρ³ .

Call this 3×3 array the **dipolar / field-gradient kernel**:

  **D_ab(ρ) = ( 3 ρ̂_a ρ̂_b − δ_ab ) / ρ³ .**

Three facts about D_ab, each checkable by hand:

1. It is **symmetric** (D_ab = D_ba) and **traceless**: its trace is
   (3 ρ̂·ρ̂ − 3)/ρ³ = (3 − 3)/ρ³ = 0. It is therefore a pure l = 2 object — five
   independent numbers, no isotropic part. (The tracelessness is the statement that
   1/ρ satisfies Laplace's equation away from the origin: ∇²(1/ρ) = 0.)
2. Its **angular shape**, read along any single axis, is (3cos²θ − 1)/ρ³ — the magic-
   angle pattern of §1. This is the form that appears in the ring-current, electric-
   field-gradient, and dipolar literatures alike; here it is simply the second
   derivative of 1/ρ.
3. Its **eigenvalues** are +2/ρ³ along ρ̂ and −1/ρ³ in the two perpendicular
   directions — the characteristic (+2, −1, −1) signature of a traceless rank-2
   tensor with one distinguished axis.

**The general kernel.** Each calculator forms, at the target i, a sum over its
neighbourhood of such a derivative, weighted by the per-source scalar:

  **K_i = Σ_{j ∈ N(i)} w_j · ∂^(l)( 1/ρ_ij ) ,**

where ∂^(l) denotes the l-th derivative (l = 1 gives the field; l = 2 gives D_ab). For
the rank-2 calculators this is, explicitly,

  **K_i = Σ_{j ∈ N(i)} w_j · ( 3 ρ̂_ij ρ̂_ij − δ ) / ρ_ij³ .**

That is the literal object: **the l-th multipole field at the target atom, produced by
the weighted source distribution of its neighbourhood.** For weights equal to charges,
the l = 2 version is exactly the *electric field gradient* (EFG) at i — a standard,
named physical quantity.

---

## 3. The literal compression: neighbourhood → low-order multipole moment

Now read K_i not as a field but as a *featurizer* — a map from "the neighbourhood" to
"a fixed, small set of numbers."

**The input** to the map is the neighbourhood as raw data: a set of source positions
relative to the target, {ρ_ij : j ∈ N(i)}, each tagged with its scalar weight w_j.
This is, in general, a large and variable amount of information — a labelled point
cloud whose size depends on how many atoms are nearby.

**The output** is the order-l multipole moment of that weighted cloud — a *fixed* set
of 2l + 1 numbers:

  - l = 0 → 1 number (the isotropic moment, T0);
  - l = 1 → 3 numbers (a vector, T1);
  - l = 2 → 5 numbers (the symmetric traceless tensor, T2).

**What the compression keeps.** It keeps the angular moments of the neighbourhood up
to order l, each computed under a fixed radial discount of 1/ρ^(l+1) (1/ρ² for the
field, 1/ρ³ for the gradient). For the rank-2 kernels this is precisely "how
anisotropic is the arrangement of matter around this atom, and along which axes," with
nearer sources counting far more than far ones.

**What the compression discards** — and this list is the honest core of the account:

1. **All angular structure above order l.** The map projects onto l ≤ 2 and is blind
   to l = 4, l = 6, … Two neighbourhoods with the *same* rank-2 moment but different
   fine angular structure map to the *same* K_i. The kernel cannot tell them apart.
2. **The radial detail.** Because the radial weighting is a fixed 1/ρ^(l+1) and then
   summed, one strong near source and a cloud of weak far ones can produce the same
   moment. The map cannot recover "how many, how far."
3. **Source identity beyond w_j.** Whatever a source atom is — its element, its
   chemistry — enters only through the single scalar weight. Two chemically different
   sources with the same weight at the same place are identical to the kernel.

So the one-line literal answer to *"what is encoded?"* is: **the low-order (l ≤ 2)
multipole moment of the weighted local environment, under a fixed inverse-power radial
kernel — and nothing else.** It is a *lossy, fixed-basis, low-angular-order* summary
of structural context.

---

## 4. Per calculation: which order each kernel is

"Geometric kernel" is not a single l. The family spans orders, unified only by being
derivatives of 1/ρ over the weighted neighbourhood. By physical type:

- **Electric field gradient (EFG), dipolar through-space, ring current
  (McConnell / Pople point-dipole):** order **l = 2**. The angular form is
  (3cos²θ − 1)/ρ³. The weight w_j is a charge (EFG), an induced magnetic moment
  (dipolar), or a ring-current intensity factor (ring). The ring-current case is
  D_ab *contracted onto the ring normal* — i.e. the rank-2 kernel evaluated with θ
  measured from the ring's normal axis — which is why its scalar reduction is the
  familiar (3cos²θ − 1)/ρ³ of Johnson–Bovey and Haigh–Mallion.
- **Buckingham electric-field effect:** order **l = 1**. This kernel uses the *field*
  itself, ~cosθ/ρ² (the first derivative of 1/ρ), not its gradient. It is a genuinely
  different multipole order, and it is worth keeping the distinction visible: not every
  "geometric kernel" we carry is rank-2.

**Verification flag (a pursuit item, not a claim).** The *mathematical shape* above is
standard physics and is what we are accounting for. The *exact* per-calculator weight
w_j, the precise radial power, and the precise contraction used in our code must be
read off `src/` before any specific per-calculator line becomes thesis prose. The
have-it-to-cite-it rule applies to our own source as much as to the literature: this
section states the form of the operation, not the implementation's constants. Reading
`src/` to pin each calculator to its exact emitted expression is explicit next-step
work.

---

## 5. The shared object: geometric moment with physical names

Here is the resolution of a tension we circled for a long time — *is the rank-2 signal
a real physical mechanism, or is it a context-encoding featurization in the manner of a
multiple-sequence alignment (MSA)?* The literal maths keeps those readings tied to the
same emitted object.

The multipole moment K_i can be read in two compatible ways:

- a named low-order field in standard physics vocabulary (e.g. the electric field
  gradient) when the scalar weights match that source model; and
- a **lossy, fixed-basis geometric descriptor** — the order-l projection of the
  weighted neighbourhood point cloud, a compression of structural context; this is the
  "featurization / MSA-like" reading.

These are not separable source tensors. The object emitted by the code is the
weighted-neighbourhood moment; chemistry enters through the scalar weight w_j. That is
the account to state.

---

## 6. Why T2 helped some atoms and not others — the Stage-1 account

The Stage-1 ridge result found, empirically, that the T2 (rank-2) content
improved the prediction of some atoms' shielding and not others'. The literal account
of *what* that means and *why*:

T2 is the l = 2 multipole field at i. By §3 it can carry exactly the part of the
shielding that is governed by that l = 2 through-space field — and nothing the field
gradient cannot see. So, reading off standard NMR (Ramsey) theory for *which* atoms
have shielding of that kind:

- **Where T2 helped — protons.** A ¹H nucleus has almost no local paramagnetic term
  (σ_para is small; hydrogen has no valence p-electron angular momentum to mix
  low-lying excited states). A proton's shielding is therefore *made of* through-space
  contributions — exactly the neighbourhood field gradient the kernel computes. The
  classical kernel tradition is, historically, a proton technology for precisely this
  reason. The compression and the physics line up, and T2 carries real signal.
- **Where T2 did not help — backbone C / N.** A backbone ¹³C or ¹⁵N has a large local
  paramagnetic term: its shielding (200–300 ppm in range) is dominated by *local*
  electronic structure — hybridisation, bond angles, the excited-state sum — which the
  neighbourhood multipole moment does not encode at all. Whatever real through-space
  l = 2 contribution exists is a thin slice under that local quantum bulk. Or the
  relevant neighbourhood structure lives above l = 2, where the projection is blind.
  Either way, the kernel cannot reach the dominant term, so T2 adds little.

This is a hedged "this or that" of exactly the kind an account is allowed — *this* =
l = 2-through-space-governed; *that* = local-quantum-governed (or above-l = 2) — and it
names what each is and how to tell them apart. The test is built into the maths: the
multipole moment can only ever carry the former, so an atom it fails to predict is
telling us its shielding lives where l ≤ 2 geometry does not.

*(Standard-theory flag: the σ_dia / σ_para split and the smallness of proton σ_para
are established NMR theory, not a claim of ours. We borrow them to read the result; we
do not assert them as findings.)*

---

## 7. Why the model "recovers" a T2 equation at all — the circularity

This is the section that keeps us honest about the symbolic-regression and equivariant
results, and it is the consequence the maths makes unavoidable.

When a model — PySR distilling a closed form, or an equivariant network — appears to
"recover" the (3cos²θ − 1)/ρ³ form, we must ask what the recovery actually demonstrates.
If the gradient field was supplied as a feature, a model that reports a
gradient-field-shaped relationship is, to a large degree, **handing us back the form of
its own input.**

State it as the separation it forces:

- **Recovery of the form is largely tautological.** If the input feature already has
  the shape (3cos²θ − 1)/ρ³, then a method that reproduces that shape (especially when
  scored against the producer's own clean kernel) is confirming that we computed what we
  computed. This is why our own findings file calls the "equation fell out" result a
  *pipeline check* and warns *do not headline it*. The literal reason is here: the form
  was supplied, so its return is expected, not discovered.
- **Only the match to the DFT target is non-circular.** The single part that is a real
  statement about physics is whether K_i (with some coefficient) correlates with the
  *independently computed DFT shielding*. That correlation is not handed to us by our
  own input; it is the empirical, faint, atom-dependent question of §6, and it is the
  only part an ablation can adjudicate.

So the account of "why do we get a T2 equation match out of the model" is conditional:
if a gradient field was put in, form-recovery measures the pipeline; target-match
measures the physics. Keeping those two apart keeps the account honest.

---

## 8. What this licenses, and what to pursue

**What we may say (account, not claim).**

- The geometric kernels are low-order (l ≤ 2) multipole moments of the weighted local
  neighbourhood — cheap, fixed-basis, equivariant compressions of structural context,
  computable per frame in seconds to hours. This is a methods statement, not a boast.
- Each is a weighted-neighbourhood moment with standard physical names in some source
  models (§5).
- T2 carries the part of shielding that is l = 2-through-space-governed, which is why it
  helps protons and the exposed and not the local-quantum-dominated backbone (§6).

**What remains an empirical question, held elsewhere (ours to trace, not asserted
here).**

- *Relevance*: if the kernel features are tried in a downstream predictor, the honest
  test is ablation and target-match to DFT — not form-recovery (§7).
- *The per-atom check*: if a model uses these kernels, compare its successes with the
  l = 2-through-space-governed account in §6.

**Verification to-dos (the pursuit).**

- Read `src/` to pin each calculator to its exact emitted expression — weight w_j,
  radial power, contraction — replacing the physical-shape statements of §4 with the
  implementation's literal forms.
- Confirm how T0 / T1 / T2 are actually decomposed and emitted in the output, so the
  five-number / three-number / one-number counts of §3 match the real array layout.
- When the equivariant run lands, check the per-atom result against §6's prediction and
  record which way it fell.

**The claim / account split, restated.** Account for the kernels as emitted geometric
features. Do not make them the Step-1 foundation.

---

*Pursue from here. Nothing above is digested; if a future reader needs the short
version, it is the title of §3 — the kernel is the low-order multipole moment of the
weighted neighbourhood — but the point of this file is to keep the long version, with
its definitions and its discards intact.*
