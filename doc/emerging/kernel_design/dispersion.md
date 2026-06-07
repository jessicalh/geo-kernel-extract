# Kernel design: the van der Waals / dispersion contribution to shielding

Scope: one input feature for the e3nn equivariant shielding predictor — the
effect of London (van der Waals) dispersion on a nucleus's magnetic shielding.
The brief is "what *should* we compute, the defensible standard way, in
equivariant form", web-grounded first, code-read second. This is a small,
subtle effect; the document says so honestly and does not inflate it.

Register: **modest and honest.** Of the through-space terms this is the one the
literature treats as marginal, hard to separate from steric/electrostatic
neighbours, and at leading order a near-contact *isotropic deshielding*. The
right output here may be a single `0e` scalar with a stated ceiling — and an
honest "it may not earn its place." That is a complete, defensible answer, not a
failure to design.

Terms are defined on first use. Nothing below is novel.

- **Shielding tensor** σ: the rank-2 object DFT computes per nucleus; the
  dispersion effect is *one additive contribution* to it. "Dispersion's
  contribution to σ" means that additive piece, not a standalone observable.
- **Dispersion / London force**: the attractive interaction between two
  non-bonded groups arising from *correlated quantum fluctuations* of their
  electron clouds — instantaneous-dipole/induced-dipole, energy ∝ −C₆/R⁶.
- **Mean-square fluctuating field** ⟨E²⟩: the time-averaged square of the
  fluctuating electric field a perturber's electron cloud imposes at the
  target nucleus. This is the quantity dispersion-shift theory is built on.
- **Irrep / l**: an irreducible representation of O(3), labelled by l (l=0
  scalar `0e`, l=1 vector, l=2 the 5-component quadrupole-like object) and
  parity. e3nn consumes features written in these.

---

## 1. What the effect is, and how the literature computes it

### 1.1 The physical picture

A non-bonded group sitting a few ångström from a nucleus does two things to
that nucleus's shielding that are *not* captured by the static electrostatic
field (the charge/EFG kernel) or by magnetic anisotropy (McConnell/ring):

1. **A fluctuating-field (true dispersion) effect.** The perturber's electron
   cloud fluctuates; those fluctuations impose a *fluctuating* electric field
   at the target nucleus. The nucleus's shielding responds quadratically to
   that field, and the time-average of the square is non-zero even though the
   mean field may be small. This is the genuine London-dispersion contribution
   to shielding.
2. **A short-range "steric"/overlap (Pauli-repulsion) effect.** At true contact
   the electron clouds overlap and distort; this deshields the compressed bond.
   In the experimental literature this is lumped together with (1) under "the
   van der Waals shift" and is hard to separate from it
   ([steric compression deshielding, Org. Lett. / PMC8154237](https://pmc.ncbi.nlm.nih.gov/articles/PMC8154237/)).

Both are real; both are *small* relative to ring-current, electrostatic, and
H-bond terms; and the field's honest position is that they are the residual,
not a leading term.

### 1.2 The standard theory: the Buckingham B-term and σ_W = −B⟨E²⟩

The canonical framing is **Marshall–Pople / Raynes–Buckingham–Bernstein (RBB)**,
within the same Buckingham (1960) electric-field expansion of the shielding that
the electrostatic kernel uses
([Buckingham 1960, Can. J. Chem.](https://cdnsciencepub.com/doi/10.1139/v60-040)).
Expand σ in the electric field **E** at the nucleus:

```
σ = σ⁰  +  A·E_z  +  B·⟨E²⟩  +  …
```

- the **A·E_z** linear term is the electrostatic/LEFE leg — *that is the
  charge-EFG kernel's territory*, not this one;
- the **B·⟨E²⟩** quadratic term is where dispersion lives. Dispersion supplies a
  *fluctuating* field whose mean is ~0 but whose **mean square ⟨E²⟩ is not**, so
  the van der Waals shift is

```
σ_W  =  B · ⟨E²⟩          (commonly written σ_W = −B⟨E²⟩ with B > 0)
```

B is the nucleus's quadratic shielding-polarizability coefficient (a *response
property*, tabulated/fitted per nucleus type, larger for ¹H than for heavy
nuclei in a given bond). The sign convention makes σ_W a **deshielding**
(positive chemical shift) for the usual B; this matches the experimental "vdW
deshielding on close contact" rule of thumb
([gas-phase intermolecular NMR, RBB binary-collision model, arXiv:1107.5261](https://arxiv.org/pdf/1107.5261);
[factors affecting ¹H shift — vdW/steric deshielding](https://conductscience.com/1h-nuclear-magnetic-resonance-nmr-chemical-shifts)).

### 1.3 The geometric content: ⟨E²⟩ as an R⁻⁶ near-neighbour sum

The dispersion ⟨E²⟩ at the nucleus, summed over surrounding perturbers j, has
the **same R⁻⁶ distance law as the dispersion energy itself** — it is built from
the same correlated polarizability product. Homer & Mohammadi's *generalized
polyatomic London dispersion theorem* gives exactly this: a net-attractive
polyatomic sum used to compute ¹H NMR gas-to-solution vdW shifts with good
agreement to experiment
([Homer & Mohammadi 1988, J. Chem. Soc. Faraday Trans. 1, 84, 2959](https://pubs.rsc.org/en/content/articlelanding/1988/f1/f19888402959)).
The defensible *geometric* descriptor is therefore

```
⟨E²⟩(target)  ≈  Σ_j  w_j / R_j⁶          (per-perturber weight w_j ∝ polarizability/ionization terms)
```

i.e. a **scalar, isotropic, short-ranged R⁻⁶ sum over near neighbours** — the
same shape as a Lennard-Jones attractive sum. The R⁻⁶ falloff makes this a
*contact* term: it is dominated by the single closest few non-bonded atoms and
is negligible beyond ~5–6 Å. The leading object is `0e`.

### 1.4 Is there an anisotropic part?

Yes, but it is a sub-leading correction, not the main effect. The gas-phase
*virial-shielding* literature extends the bare R⁻⁶ dispersion term with R⁻³ and
R⁻¹⁰ companions and notes that the **shielding (and the magnetizability) acquire
anisotropic pieces** beyond the isotropic R⁻⁶ core
([weak intermolecular interactions in gas-phase NMR, arXiv:1107.5261](https://arxiv.org/pdf/1107.5261)).
The honest reading: the **dominant, well-established dispersion contribution to
shielding is isotropic** (the ⟨E²⟩·B scalar); an anisotropic l=2 piece exists at
the level of higher virial terms and orientation-dependent polarizability, but
it is smaller, less standardised, and entangled with the steric-overlap term. So
the clean statement is "primarily `0e`, with a real but weak `2e` tail."

### 1.5 How DFT (our ground truth) holds dispersion — and the fixed-DFT bind

This is the load-bearing honesty for our setup. **In a QM shielding calculation
there is no separate "dispersion shielding" channel.** Dispersion enters in two
non-separable ways:

- it shapes the **geometry** (dispersion-corrected functionals — D3/D2/TS —
  pull non-bonded contacts to their right distances, and the shielding then
  reflects that geometry)
  ([dispersion-corrected DFT for NMR crystallography, JPCA 2020](https://pubs.acs.org/doi/10.1021/acs.jpca.0c06372);
  [DSD-PBEP86 dispersion-corrected shifts beat MP2](https://par.nsf.gov/servlets/purl/10165326));
- it is folded into the **total electron response** that produces σ. The total σ
  the DFT writes already *contains* whatever dispersion does, inseparably from
  everything else.

Our DFT campaign emitted **total shielding only** and is not re-runnable. So we
**cannot** extract a per-nucleus "DFT dispersion contribution" to train against
or to validate a dispersion kernel pointwise. A dispersion feature can only ever
be a **geometric surrogate** computed from coordinates + tabulated constants,
offered to e3nn as one more input channel; the model decides whether it carries
residual signal once ring/electrostatic/H-bond terms are in. This is the
strongest reason the term is "may not earn its place" rather than "must be in."

### 1.6 Double-counting: the honest separation problem

The R⁻⁶ near-contact scalar is **correlated with several terms already in the
feature set**, and the literature is explicit that vdW/steric deshielding is
"difficult to distinguish from" magnetic anisotropy and electrostatics
([PMC8154237](https://pmc.ncbi.nlm.nih.gov/articles/PMC8154237/)):

- a close aromatic contact contributes to **ring-current** *and* to this R⁻⁶ sum;
- a close polar group contributes to the **charge/EFG** field *and* to it;
- the steric-overlap part overlaps the **H-bond** kernel at H-bond distances.

This is not a reason to omit it, but it *is* the reason to (a) emit it as its
own labelled channel so the fit can decide its marginal value against the
others, and (b) not claim a clean physical "dispersion shielding" readout from
it. Empirical predictors (SHIFTX2, SPARTA, ProCS15) carry proximity/steric
descriptors of exactly this flavour as minor correction terms, not headline
physics
([SHIFTX2](https://pmc.ncbi.nlm.nih.gov/articles/PMC3085061/);
[ProCS15](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4662583/)).

---

## 2. Turning it into an equivariant e3nn feature

e3nn consumes features as irreps (`Nx0e + Mx1o + Kx2e`), each transforming by
the Wigner-D matrices under rotation; physical scalars/tensors can be fed as
steerable node features and the network forms invariants by tensor product
([e3nn irreps guide](https://docs.e3nn.org/en/stable/guide/irreps.html);
[SEGNN, Brandstetter et al., ICLR 2022](https://arxiv.org/abs/2110.02905)).

The dispersion effect's irrep content follows directly from §1:

- **Primary: `0e`.** The ⟨E²⟩·B scalar is a rotational invariant by construction
  (⟨E²⟩ is a scalar; the leading vdW shift is isotropic deshielding). This is the
  honest, single-irrep statement of the effect. A scalar fed to an equivariant
  network is perfectly legitimate — it simply enters as an l=0 node feature and
  modulates the message passing; it does not need to be a tensor to be useful.
- **Optional, weak: `2e`.** If one wants the anisotropic tail of §1.4, build the
  symmetric-traceless R⁻⁶/R⁻⁸ dispersion tensor (the dipolar-shaped
  `3 d̂⊗d̂ − I` kernel weighted by C₆/R⁶, summed over neighbours), keep its
  symmetric part, and emit the five l=2 components. State plainly that this is
  the sub-leading correction, more model-risk than the scalar, and may be noise.
- **No `1o`.** There is no defensible vector (l=1) dispersion shielding object —
  the effect is even-parity and, at leading order, scalar. Do not synthesise one.

Parity: a shielding contribution is even under inversion, so `0e` / `2e` are the
correct parities; an outer-product-of-displacements construction must have its
parity declared correctly to e3nn or the equivariance guarantee is silently
wrong (a correctness point, not a physics ambiguity).

---

## 3. The defensible, sane recommendation

**Compute, per target nucleus, the dispersion ⟨E²⟩ near-contact scalar as a `0e`
feature, geometry-only and unscaled, over *all* non-bonded near neighbours — and
treat it as a candidate the fit may reject.** Concretely:

1. **Geometry — an isotropic R⁻⁶ near-neighbour sum, not a ring-only sum.**
   The physics (§1.3) is a sum over *every* non-bonded perturber, weighted by a
   per-element dispersion coefficient:

   ```
   D₀(target) = Σ_j  C₆(target, type_j) · S(R_j) / R_j⁶      → emit as 0e (Å⁻⁶, unscaled)
   ```

   over all non-bonded atoms j within a KD-tree cutoff (R⁻⁶ makes ~5–6 Å ample),
   with the standard near-bond/self exclusions and a smooth switching taper S(R)
   so the truncated sum does not jump as atoms cross the cutoff. Emit the
   geometry **unscaled** (Å⁻⁶); the C₆/B prefactor is the calibration layer
   (point 3). This is the LJ-attractive-shaped descriptor every empirical
   predictor uses for proximity, computed honestly.

2. **Stratify the scalar by perturber element/chemistry, as parallel `0e`
   channels** (e.g. C, N, O, S, aromatic-C separately), rather than one pooled
   number. Different elements have different C₆/polarizability, and keeping them
   separate (a) lets the fit weight by chemistry and (b) lets the aromatic-C
   channel be zeroed if it double-counts the ring kernel (§1.6). More channels,
   not more assumptions.

3. **The scale (C₆ / the Buckingham B) is a fitted/calibrated coefficient with a
   stated source, not a hard constant.** Per-element C₆ are tabulated (Grimme D3
   dispersion coefficients are the obvious, defensible, widely-used source; or
   atomic-polarizability/ionization estimates per the London theorem of §1.3),
   and the nuclear B-coefficient varies by nucleus type and is only loosely
   known. Emit unscaled geometry and let the per-channel C₆·B product be a
   calibratable coefficient — exactly as ring-current intensity and McConnell Δχ
   are treated in the sibling kernels. The honest state is "the *shape* (R⁻⁶,
   isotropic, contact-dominated) is solid; the *magnitude* is a fitted scale."

4. **Optional anisotropic `2e` channel — only if it earns it.** If wanted, also
   emit the symmetric-traceless dispersion tensor's five l=2 components per
   element-channel (§2). Label it the sub-leading tail; expect it to be weak and
   be ready for the fit to drop it. Do not make it the primary object.

5. **Fast linalg where it buys defensibility, not novelty.** KD-tree
   neighbourhood; a per-pair scalar (and optional 3×3) accumulation; a fixed
   Cartesian→real-l=2 projection done once. All cheap and exact; nothing here is
   a new method.

**Why this shape and not the alternatives:**
- *A scalar is the right primary, not a cop-out.* The effect is isotropic at
  leading order; handing e3nn a `0e` ⟨E²⟩ descriptor is the faithful form. A
  tensor would over-promise structure the physics doesn't robustly supply.
- *All non-bonded neighbours, not aromatic rings only.* Dispersion is a
  general close-contact interaction (every non-bonded atom disperses); confining
  it to aromatic-ring vertices models only the aromatic slice of a general
  effect and conflates it with the ring-current geometry. The web physics is an
  all-neighbour R⁻⁶ sum.
- *Geometric surrogate, not a DFT-validated readout.* Because the DFT never
  emitted a separable dispersion shielding (§1.5), this feature is a coordinate
  descriptor offered to the model, validated only by whether it adds marginal
  fit — never matched pointwise to a "true" dispersion shielding.

---

## 4. Where the physics is clean, and where it is not

**Clean (high confidence):**
- **The geometric shape.** "Isotropic, contact-dominated, R⁻⁶ sum over near
  non-bonded neighbours" is the well-established form of the dispersion ⟨E²⟩ and
  maps cleanly to a `0e` feature. The geometry, given coordinates, is exact.
- **The irrep.** Leading order is a scalar (`0e`), even parity — unambiguous.
  No frame, no projection, no judgement call.

**Not clean (state as footnotes, not bugs):**
- **Magnitude / scale.** The C₆ coefficients and especially the nuclear
  B-coefficient are tabulated/fitted with real uncertainty; the right design
  makes the scale a calibrated per-channel coefficient, not a constant. The
  shape is solid; the size is fitted.
- **Separability from neighbouring terms.** vdW/steric deshielding is
  *explicitly* hard to separate from magnetic anisotropy, electrostatics, and
  H-bonding (§1.6). Its marginal value over those terms is an empirical question
  for the fit, and the honest prior is "small, possibly redundant."
- **The anisotropic tail.** A `2e` part exists (higher virial terms) but is
  weak, less standardised, and entangled with steric overlap. Defensible to
  offer as a parallel channel; not defensible to lean on.
- **No DFT anchor.** Unlike the cleaner kernels, there is no separable QM
  dispersion shielding in our fixed campaign to check against (§1.5). The
  feature is a surrogate, full stop.

**Single most defensible package:** a per-element-stratified, geometry-only
`0e` R⁻⁶ near-contact scalar (⟨E²⟩-shaped), with the C₆·B scale left as a
calibratable per-channel coefficient, and the aromatic channel zeroable to
avoid ring double-counting. The optional `2e` tail is a legitimate but weak
add-on. **And honestly: this term may not earn its place once ring,
electrostatic, and H-bond kernels are in — that is an acceptable outcome for a
small effect, and the design above is built so the fit can say so cleanly.**

---

## 5. Neutral note: current code vs this recommendation

Descriptive only, no verdict. The current library (`src/DispersionResult.{h,cpp}`)
computes a dispersion kernel **summed over aromatic-ring vertices** within range
of each atom: per vertex, with unit C₆ = 1, it forms the traceless tensor
`K_ab = S(r)·(3 d_a d_b / r⁸ − δ_ab / r⁶)` and the isotropic scalar `S(r)/r⁶`
(both Å⁻⁶), where S(r) is a CHARMM-form smooth switching taper between a
configurable onset and cutoff. It applies a ring near-field filter and a
through-bond vertex exclusion, accumulates per ring-neighbourhood and per ring
type, decomposes to a SphericalTensor, and emits `disp_shielding` (N,9),
`disp_per_type_T0` (N,8), and `disp_per_type_T2` (N,40); the C₆ scale is left
as "learnable per ring type" and the geometry is emitted unscaled. Relative to
the recommendation, the descriptive differences are: (a) the present sum is over
**aromatic-ring vertices** whereas the web physics is an **all-non-bonded-neighbour**
R⁻⁶ sum (the ring-only form models the aromatic slice of a general close-contact
effect); (b) the present object leads with the **full traceless tensor + T2**
whereas the recommendation leads with the **isotropic `0e` ⟨E²⟩ scalar** and
treats the `2e` tail as the weak, optional add-on; and (c) the present
stratification is **per ring type** whereas the recommendation stratifies by
**perturber element/chemistry** with an aromatic channel zeroable against the
ring kernel. The switching-taper, KD-tree neighbourhood, traceless decomposition,
and unscaled-geometry/learnable-scale discipline are already in place. This
paragraph is description, not a verdict.

---

## Sources

- Buckingham, *Chemical shifts in the NMR spectra of molecules containing polar
  groups*, Can. J. Chem. 38, 300 (1960) — the field expansion (A·E + B·E²).
  [cdnsciencepub 10.1139/v60-040](https://cdnsciencepub.com/doi/10.1139/v60-040)
- *Weak intermolecular interactions in gas-phase NMR* — RBB binary-collision
  model, σ_W from ⟨E²⟩, R⁻³/R⁻⁶/R⁻¹⁰ virial-shielding terms, anisotropic pieces.
  [arXiv:1107.5261](https://arxiv.org/pdf/1107.5261)
- Homer & Mohammadi, *NMR van der Waals chemical shifts from a generalized
  polyatomic London dispersion theorem*, J. Chem. Soc. Faraday Trans. 1, 84,
  2959 (1988) — net-attractive polyatomic R⁻⁶ sum for ¹H vdW shifts.
  [RSC F19888402959](https://pubs.rsc.org/en/content/articlelanding/1988/f1/f19888402959)
- *The effect of dispersion interaction on nuclear magnetic shielding* — pair
  shift ∝ van der Waals energy × unperturbed shielding (closure approximation).
  [ScienceDirect 002223646990081X](https://www.sciencedirect.com/science/article/abs/pii/002223646990081X)
- Steric-compression / vdW deshielding, hard to separate from anisotropy —
  [Org. Lett., PMC8154237](https://pmc.ncbi.nlm.nih.gov/articles/PMC8154237/);
  [¹H shift factors incl. vdW/steric](https://conductscience.com/1h-nuclear-magnetic-resonance-nmr-chemical-shifts)
- Dispersion in DFT shielding (not separable; geometry + total response):
  [dispersion-corrected DFT for NMR crystallography, JPCA 2020](https://pubs.acs.org/doi/10.1021/acs.jpca.0c06372);
  [DSD-PBEP86 dispersion-corrected shifts](https://par.nsf.gov/servlets/purl/10165326)
- Empirical predictors carrying proximity/steric terms —
  [SHIFTX2](https://pmc.ncbi.nlm.nih.gov/articles/PMC3085061/);
  [ProCS15](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4662583/)
- e3nn irreps; SEGNN (physical tensors as steerable features) —
  [e3nn docs](https://docs.e3nn.org/en/stable/guide/irreps.html);
  [Brandstetter et al., ICLR 2022](https://arxiv.org/abs/2110.02905)
