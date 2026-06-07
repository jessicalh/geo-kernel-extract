# Kernel design: the aromatic π-system electric-quadrupole effect on shielding

Scope: one input feature for the e3nn equivariant shielding predictor — the
through-space shielding a nucleus feels from the **electric (electrostatic)
quadrupole moment of a nearby aromatic π system**. This is the *electrostatic*
face of an aromatic ring (its permanent charge distribution), **not** the
magnetic ring current and **not** a point-charge field. The job is to state the
defensible, standard way to compute and featurize it; web-grounded first, our
code read last and described neutrally.

Terms, defined as they appear:

- **Molecular electric quadrupole moment** Θ: the leading nonzero electric
  multipole of a symmetric, charge-neutral, dipole-free molecule like benzene.
  It is a **rank-2, symmetric, traceless** Cartesian tensor (five independent
  components); for an aromatic ring it is axially symmetric about the ring
  normal, so one scalar Θ_zz (the "Θ" below) describes it. Benzene's is large
  and negative, Θ_zz(C₆H₆) ≈ −29 × 10⁻⁴⁰ C·m² — negative charge above and below
  the ring face, compensating positive charge in the plane
  ([Battaglia et al., *Mol. Phys.* 1980 / quadrupole moments of benzene & C₆F₆](https://pubs.rsc.org/en/content/articlelanding/1980/f2/f29807600648)).
- **EFG** (electric-field gradient): V_αβ = ∂²φ/∂r_α∂r_β at the nucleus — a
  rank-2 symmetric traceless tensor (l=2) for any field sourced by charge
  *external* to the nucleus (Laplace).
- **Buckingham field effect**: the standard mechanism by which an external
  electric field **E** and field-gradient **∇E** perturb a nucleus's shielding,
  σ = σ⁰ + (response)·E + (response)·∇E + … (Buckingham 1960). This is the same
  mechanism the charge/EFG kernel uses; here the *source* of E and ∇E is the
  ring's quadrupole rather than point charges.
- **l / irrep**: SO(3)/O(3) irreducible representation; `0e` scalar, `1o`/`1e`
  vector, `2e` the 5-component quadrupolar tensor. e3nn consumes these.

---

## 1. What the effect is, and how the literature computes it

An aromatic ring has no net charge and no net dipole, but it has a substantial
**permanent electric quadrupole moment** — the famous "negative face, positive
edge" charge distribution that drives cation–π and π–π (T-shaped, parallel-
displaced) recognition
([Cation–π interaction — Wikipedia](https://en.wikipedia.org/wiki/Cation%E2%80%93%CF%80_interaction);
[Dougherty, *Cation–π in chemistry & biology*, Science 1996](https://www.science.org/doi/10.1126/science.271.5246.163)).
That static charge distribution produces, at a nearby nucleus, an **electric
field and an electric-field gradient**, and those perturb the nucleus's shielding
through the Buckingham field-effect mechanism — exactly the same physical channel
as a point charge, but driven by the ring's *quadrupole* rather than by monopoles.

### 1a. The field of a point quadrupole (the clean part)

For a charge distribution whose leading moment is the quadrupole Θ, the
multipole expansion gives a potential that falls as **1/R³**, with angular shape
P₂(cosθ) = (3cos²θ−1)/2 for the axial case; the **electric field** (its gradient)
falls as **1/R⁴**, and the **field-gradient (EFG)** as **1/R⁵**
([Electric quadrupole — potential ∝ 1/R³, field ∝ 1/R⁴, gradient ∝ 1/R⁵](https://en.wikipedia.org/wiki/Quadrupole);
[Multipole expansion of the electric field — Citizendium](https://en.citizendium.org/wiki/Multipole_expansion_of_electric_field)).
The axial quadrupole potential is

```
φ(R,θ)  =  Θ · (3cos²θ − 1) / (4 R³)          (θ from the ring normal)
```

This is standard electrostatics, parameter-free given the geometry and Θ. The
EFG it produces is obtained by twice-differentiating φ — cleanly done with
Stone's interaction (T-) tensors, V_αβ = −(Θ/2)·T_αβγδ n_γ n_δ
([Stone, *The Theory of Intermolecular Forces*, OUP — T-tensor machinery](https://global.oup.com/academic/product/the-theory-of-intermolecular-forces-9780199672394)).
Because Θ is traceless and the field is source-free at the nucleus, the
resulting EFG is **symmetric and traceless — a pure l=2 object**, the same
mathematical class as the point-charge EFG (§ charge_efg.md).

### 1b. The mechanism on shielding (Buckingham), and that it is a *named*,
separate contribution in the shift literature

The way an external electrostatic field changes shielding is Buckingham's field
expansion: the **linear field E** (l=1) couples to the dipole shielding
polarizability, and the **field-gradient ∇E** (l=2) couples to a corresponding
response tensor
([Buckingham 1960, *Can. J. Chem.*](https://cdnsciencepub.com/doi/10.1139/v60-040);
[de Dios, Pearson & Oldfield, *Science* 1993 — charge-field perturbation for protein shifts](https://www.science.org/doi/10.1126/science.8502992)).
For an aromatic ring, the *source* of that E and ∇E is the ring's own charge
distribution, whose leading term is the quadrupole. The aromatic-shift
literature treats this **electric/electrostatic contribution as distinct from the
magnetic ring current**:

- Abraham's "proton chemical shifts in NMR" series decomposes aromatic ¹H shifts
  into **separate** ring-current (magnetic anisotropy), π-electron-density, and
  **electric-field** terms — the electric-field term arising from the ring's
  charge distribution / C–H dipoles, modelled via the Buckingham field effect,
  not from the ring current
  ([Abraham et al., *J. Chem. Soc. Perkin Trans. 2* 2000, Part 14](https://pubs.rsc.org/en/content/articlelanding/2000/p2/a907830d)).
- RNA/protein shift analyses likewise partition the total into **ring-current
  and intramolecular-electrostatic (electric-field)** contributions and find they
  differ in magnitude and even dominate for different nuclei
  ([Case et al., *J. Phys. Chem. B* 2012 — ring-current vs electric-field contributions to RNA-base shifts](https://pubs.acs.org/doi/10.1021/jp3057306)).

So the π-quadrupole effect is a real, literature-recognised, **electrostatic**
contribution, computed by (i) representing the ring as a quadrupole source,
(ii) evaluating its E / ∇E at the nucleus, (iii) feeding the Buckingham response.
What is *not* standard is treating the ring strictly as a *single point*
quadrupole — most explicit protein/RNA work uses **distributed partial atomic
charges** on the ring atoms and sums the Coulomb field/EFG, of which the point
quadrupole is the leading far-field multipole (see §1d).

### 1c. How this differs from — and what it adds beyond — ring current and charge/EFG

This is the crux of the kernel, so it is worth stating sharply. There are three
*physically distinct* through-space effects an aromatic ring exerts:

| Effect | Physical origin | Field that acts | Mechanism | Distance |
|---|---|---|---|---|
| **Ring current** (ring.md, mcconnell.md) | *Induced* circulation of π electrons in B₀ | *Magnetic* secondary field | magnetic, σ from induced currents | dipole ~1/R³ |
| **Charge / EFG** (charge_efg.md) | *Permanent* point partial charges (monopoles) | *Electric* E (1/R²) and EFG (1/R³) | Buckingham field effect | monopole-sourced |
| **π-quadrupole** (this) | *Permanent* quadrupole of the neutral π charge cloud | *Electric* E (1/R⁴) and EFG (1/R⁵) | Buckingham field effect | quadrupole-sourced |

- **vs ring current:** different *physics entirely* — magnetic (induced,
  field-dependent currents) vs electrostatic (permanent charge distribution). A
  ring at zero applied field still has its quadrupole; it has no ring current.
  They share the ring geometry and the (3cos²θ−1) angular skeleton (both are
  axial l=2 about the normal), which makes them look alike on a plot, but they
  are not the same term and do not double-count *each other*. The honest caveat:
  any *empirical* ring-current calibration fit to DFT/experiment may have already
  absorbed part of the electrostatic effect into its fitted intensity, so the two
  are separable in principle but can be entangled in fitted coefficients.

- **vs charge/EFG:** *same mechanism, different multipole order of the same
  charge distribution.* If the charge/EFG kernel already sums a Coulomb field
  from **partial atomic charges placed on the aromatic ring atoms**, then the
  ring's quadrupole is **already contained** in that sum — the point quadrupole
  is just the leading far-field multipole of those same ring charges. In that
  case a separate point-quadrupole term is a **redundant re-expansion of charge
  the EFG kernel already carries: a double-count.** The π-quadrupole term earns
  independent signal *only* if (a) the charge/EFG kernel **excludes** aromatic-
  ring atoms from its charge sum (so the ring's electrostatics are delegated to
  this kernel), **or** (b) the ring is represented by a *better* quadrupole than
  the force-field partial charges give — e.g. a literature/QM molecular Θ that
  captures the true π-cloud anisotropy the point charges approximate poorly at
  short range. This is the single most important design decision for the kernel
  and must be made explicit (see §3, §4).

**The one-line answer to "what does it add":** beyond ring current it adds the
ring's *permanent electrostatic* field-effect (a different physical channel);
beyond charge/EFG it adds **nothing unless the ring atoms are removed from the
point-charge sum or the quadrupole is sourced more accurately than the partial
charges** — otherwise it is the same charge counted twice.

### 1d. Point quadrupole vs distributed charges (the near-field ceiling)

A single point quadrupole at the ring centre is a multipole expansion truncated
at l=2. Like every multipole truncation, it is a *far-field* statement: accurate
when R ≫ ring size, degrading as the nucleus approaches the ring. For benzene
(radius ≈ 1.4 Å) many protein ring–nucleus contacts sit at 3–5 Å — the regime
where the point-multipole picture is coarse and the honest source model is the
**distributed ring charges** (partial charges on the six/five atoms, summed
Coulomb field and EFG), of which the point quadrupole is only the leading term.
This is the same near-field caveat the ring-current and McConnell kernels carry,
and the same resolution: prefer a distributed source at contact range, keep the
point quadrupole as the clean far-field statement.

---

## 2. Featurizing for e3nn — which irreps

The effect enters shielding through **two multipole orders of the field**, and
each is a clean irrep:

- the quadrupole's **electric field E** at the nucleus is a vector → **`1o`**
  (l=1, odd parity), magnitude ∝ Θ/R⁴;
- the quadrupole's **EFG** at the nucleus is symmetric-traceless (Laplace) →
  **`2e`** (l=2, even parity), the 5-component tensor, magnitude ∝ Θ/R⁵.

This mirrors the charge/EFG design exactly (`1o ⊕ 2e`), and for the same reason:
the Buckingham expansion couples to **both** the field and its gradient, the two
orders carry independent, non-redundant geometric information (different distance
weighting, different l), and an equivariant network forms the rotational
invariants (E·b̂, E:V, …) itself by tensor product rather than us pre-projecting
into a local frame. Emit raw in the molecular frame, geometry unscaled.

Whether to emit **both** `1o` and `2e` or **only `2e`** is a real choice. The
EFG `2e` is the cleaner, more-distinctive object (pure l=2 by Laplace, the part a
tensor model exists to use); the `1o` field is also legitimate but at this
multipole order falls fast (1/R⁴) and overlaps more with what the charge/EFG
`1o` already carries. The defensible default is to emit the **`2e` EFG as the
primary feature** and the **`1o` field as a parallel channel** the fit can weigh
or drop — consistent with the charge/EFG sibling, which feeds both orders and
lets the model decide. The clean single-irrep statement of *this* kernel is the
`2e`.

This is the field-standard recipe for NMR tensors into equivariant GNNs:
decompose to spherical-tensor irreps and let Clebsch–Gordan tensor products
combine them, never nine raw Cartesian numbers and never one collapsed scalar
([Venetos et al., equivariant GNN for full ²⁹Si shift tensors, *JPCA* 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/);
[Ben Mahmoud et al., solid-state NMR parameters via spherical-tensor decomposition, arXiv:2412.15063](https://arxiv.org/abs/2412.15063);
[e3nn irreps guide](https://docs.e3nn.org/en/stable/guide/irreps.html)).

---

## 3. The defensible, sane recommendation

**Compute, per target nucleus, the electrostatic field of each nearby aromatic
ring's quadrupole and emit its EFG as a `2e` irrep (primary), with the field as a
parallel `1o` channel — but first settle the double-count boundary with the
charge/EFG kernel, because that boundary is what makes this term carry
independent signal or not.**

Concretely:

1. **Decide the charge-partition first (load-bearing).** Pick one and state it:
   - *Delegated:* the charge/EFG kernel **excludes aromatic-ring atoms** from its
     point-charge sum, and this kernel carries all aromatic-ring electrostatics
     via the quadrupole. Clean separation, no double-count.
   - *Upgrade:* the charge/EFG kernel keeps ring atoms (as partial charges), and
     this kernel sources Θ from a **molecular/QM quadrupole** (literature benzene
     Θ scaled per ring type, or a per-ring QM moment) that is a *better* source
     than the partial charges — justified only if it demonstrably captures π-cloud
     anisotropy the partial charges miss at short range.
   - *Redundant (avoid):* both kernels sum the same ring partial charges — this is
     a double-count and the π-quadrupole term should then be zeroed, exactly as an
     aromatic Δχ term is zeroed when a ring-current kernel is present (mcconnell.md
     §1.6).

2. **Source model — distributed charges preferred at contact range; point
   quadrupole as the clean far-field form.** At 3–5 Å the point quadrupole is
   coarse; the honest production source is the summed Coulomb field/EFG of the
   ring's atomic partial charges (the "Upgrade" path's QM Θ is the alternative
   when a better moment exists). Keep the point-quadrupole closed form as the
   documented far-field anchor and sanity check.

3. **Output object per nucleus.** Sum the EFG contributions of all aromatic rings
   within a generous cutoff (KD-tree neighbourhood; the field is shorter-ranged
   than the dipole, 1/R⁵, so the cutoff can be tighter than the ring-current one),
   keep symmetric-traceless, decompose:
   - `2e`: the 5-component l=2 EFG tensor (primary feature).
   - `1o`: the quadrupole's electric field at the nucleus (parallel channel).
   Emit **per aromatic ring type** (Phe, Tyr, Trp 5-/6-membered, His) so the model
   can carry a different effective Θ per chemistry, and emit the **geometry
   unscaled** (Å⁻⁵ for the EFG, Å⁻⁴ for the field) so the quadrupole magnitude Θ
   (and its Buckingham response coefficient) is a **fitted/calibrated coefficient**
   or a parallel literature-scaled channel — not a hard-coded constant. This is the
   house pattern: emit the clean geometry, fit the uncertain scale.

4. **Fast linalg where it buys defensibility, not novelty.** KD-tree ring
   neighbourhoods; closed-form Stone T-tensor (or distributed Coulomb sum) per
   ring; explicit traceless projection after accumulation to kill floating-point
   trace drift; one fixed, documented Cartesian→real-l=2 map. All cheap and exact.

**Why this and not the alternatives:**
- *A scalar (3cos²θ−1)/R⁴ alone:* the orientation-skeleton of the effect, but it
  is one `0e`-like number handed to a model whose purpose is orientation. Keep it
  only as a sanity/baseline scalar derived from the tensor, never as the primary
  feature.
- *Point quadrupole as production source at contact range:* it is the leading
  multipole, not the field — defensible as the far-field anchor, coarse at 3–5 Å.
  Use the distributed charges (or a justified QM Θ) for production, the point form
  for the check.
- *Re-deriving the ring's electrostatics from DFT per nucleus:* that is the ground
  truth and belongs as the **target/anchor**, not as a hand-built input — and the
  fixed-DFT campaign wrote shielding only, so no new QM field/EFG is available
  anyway. Literature/QM Θ and partial charges are the established surrogates.

---

## 4. Where the physics is clean, and where it is not

**Clean:**
- **The target object's form is unambiguous:** the quadrupole-sourced EFG is a
  symmetric, traceless l=2 tensor (Laplace), decomposing exactly to `2e`; the
  field is `1o`. No judgement call on the irreps. The (3cos²θ−1) angular skeleton
  and the 1/R³ / 1/R⁴ / 1/R⁵ scaling of potential / field / EFG are textbook
  electrostatics.
- **The mechanism is established and named:** Buckingham field effect; the
  aromatic-shift literature (Abraham; Case et al.) treats the ring's electric
  contribution as a real term distinct from the ring current.
- **It is genuinely a different physical channel from the ring current** —
  electrostatic vs magnetic — so it is not redundant *with the ring-current
  kernel* by construction.

**Not clean (honest uncertainties, stated as footnotes not bugs):**
- **Double-counting against charge/EFG is the central risk.** If the charge/EFG
  kernel already sums the ring's partial atomic charges, the quadrupole is the
  same charge re-expanded — redundant. The term carries independent signal only
  under an explicit charge-partition decision (§3.1). This is a *bookkeeping*
  correctness point, not a physics ambiguity, but it is the thing most likely to
  make the kernel null or misleading if left implicit.
- **Quadrupole magnitude Θ and the Buckingham response.** The geometry is exact;
  the *scale* (the ring's effective Θ per type, and the shielding-response
  coefficient) is the soft part — sourced from literature benzene Θ, QM moments,
  or fitted, and treated as a calibrated coefficient with a stated band, not a
  constant.
- **Point vs distributed source at contact range.** The point quadrupole is a
  far-field truncation; at the 3–5 Å contacts that matter for proteins it is
  coarse, and a distributed-charge (or QM-moment) source is the honest production
  choice. Same ceiling the ring-current and McConnell kernels carry.
- **Entanglement with empirical ring-current calibrations.** A ring-current
  intensity fitted to data may have silently absorbed some of the electrostatic
  effect; the two are separable in principle but can co-vary in fitted
  coefficients. Worth disclosing, not fixable here.

If forced to the single most defensible package: **the ring-quadrupole EFG as a
`2e` feature, emitted per ring type with geometry unscaled and Θ as a fitted/
literature-scaled coefficient, under an explicit decision that the charge/EFG
kernel does not also carry the same ring charges** — with the `1o` field as a
parallel channel and the point-quadrupole closed form as the far-field anchor.

---

## 5. Neutral note: current code vs this recommendation

Descriptive only, no verdict. The library (`src/PiQuadrupoleResult.{h,cpp}`)
computes, per (aromatic ring, atom) pair within a KD-tree cutoff, the EFG of a
**point axial quadrupole at the ring centre along the ring normal**, via Stone's
T-tensor: G_αβ = T_αβγδ n_γ n_δ — a symmetric, traceless, pure-l=2 tensor with
1/R⁵ leading decay (it documents T0=T1=0 and the Laplace tracelessness). The
−Θ/2 prefactor is deliberately **not** applied; the geometry is emitted unscaled
in Å⁻⁵ and the per-ring-type quadrupole strength is left as a learnable/fitted
coefficient. It also stores the scalar (3cos²θ−1)/R⁴, labelled the **Buckingham
A-term** (the quadrupole's *electric field* feeding an isotropic/T0 shift), as a
separate quantity from the EFG→T2 coupling. It accumulates the total EFG per atom
and decomposes to a `SphericalTensor`, and emits `pq_shielding` (N×9 packed
T0/T1/T2), `pq_per_type_T0` (N×8, the per-ring-type A-term scalar), and
`pq_per_type_T2` (N×40, per-ring-type T2). Filters include a dipolar/near-field
guard (with a comment noting the quadrupole's larger convergence radius) and a
ring-bonded topological exclusion; the h5-reader mirrors this as
`QtPiQuadrupoleGroup`. Relative to §3: the current path already emits the clean
**`2e` EFG** as a per-type tensor with **unscaled geometry and a fitted scale**,
already carries the **field-order scalar separately** (the A-term, the `1o`-family
content in collapsed form), and already keeps it as its own kernel — closely
matching the recommendation's shape. The descriptive gaps are that (i) the
field-order content is emitted as a *scalar* (3cos²θ−1)/R⁴ rather than the full
`1o` vector field, (ii) the source is the **point** quadrupole rather than a
distributed-charge field at contact range (the code itself flags the larger
near-field convergence radius), and (iii) **the charge-partition decision against
the charge/EFG kernel — whether those ring atoms are also summed there — is not
settled in this file**, which §3.1/§4 name as the load-bearing choice for whether
the term carries independent signal. This paragraph is a description, not a
verdict.

---

## Sources

- [Battaglia / quadrupole moments of benzene & hexafluorobenzene — *J. Chem. Soc. Faraday Trans. 2* 1980](https://pubs.rsc.org/en/content/articlelanding/1980/f2/f29807600648)
- [Dougherty, *Cation–π interactions in chemistry and biology* — Science 1996](https://www.science.org/doi/10.1126/science.271.5246.163)
- [Cation–π interaction — Wikipedia (benzene quadrupole, negative-face charge distribution)](https://en.wikipedia.org/wiki/Cation%E2%80%93%CF%80_interaction)
- [Electric quadrupole — potential 1/R³, field 1/R⁴, gradient 1/R⁵; traceless rank-2 tensor — Wikipedia](https://en.wikipedia.org/wiki/Quadrupole)
- [Multipole expansion of the electric field — Citizendium](https://en.citizendium.org/wiki/Multipole_expansion_of_electric_field)
- [Buckingham 1960, chemical shifts of molecules with polar groups (field expansion) — *Can. J. Chem.*](https://cdnsciencepub.com/doi/10.1139/v60-040)
- [de Dios, Pearson & Oldfield, charge-field perturbation for protein shifts — Science 1993](https://www.science.org/doi/10.1126/science.8502992)
- [Abraham et al., proton chemical shifts: ring current vs π-electron vs electric-field terms — *J. Chem. Soc. Perkin Trans. 2* 2000, Part 14](https://pubs.rsc.org/en/content/articlelanding/2000/p2/a907830d)
- [Case et al., ring-current vs electric-field contributions to RNA-base shifts — *J. Phys. Chem. B* 2012](https://pubs.acs.org/doi/10.1021/jp3057306)
- [Stone, *The Theory of Intermolecular Forces* (T-tensor machinery for multipole fields) — OUP](https://global.oup.com/academic/product/the-theory-of-intermolecular-forces-9780199672394)
- [Venetos et al., equivariant GNN for full ²⁹Si NMR shift tensors — *JPCA* 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/)
- [Ben Mahmoud et al., GNN prediction of solid-state NMR parameters via spherical-tensor decomposition — arXiv:2412.15063](https://arxiv.org/abs/2412.15063)
- [e3nn irreps guide](https://docs.e3nn.org/en/stable/guide/irreps.html)
</content>
</invoke>
