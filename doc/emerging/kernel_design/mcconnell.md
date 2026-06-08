Historical -- not current truth; see doc/emerging/kernel_design/mcconnell_spec.md and doc/emerging/CONTROLLING_SPEC.md.

# McConnell neighbour magnetic-anisotropy — kernel design

What this is: a design note for **one input feature** to an e3nn equivariant
shielding predictor — the through-space shielding a nucleus feels from the
anisotropic magnetic susceptibility (Δχ) of a neighbouring bond or group. The
brief is "what *should* we compute, the defensible standard way", not a verdict
on what we currently compute. It is web-grounded first, code-read second.

> **GROUNDING UPDATE — 2026-06-06.** An independent code + web + corpus grounding pass
> (**[`mcconnell_grounding_agent1.md`](mcconnell_grounding_agent1.md)**) confirmed the
> replace-`M`-with-the-clean-propagator direction and **corrected a parity error in this note**:
> the McConnell tensor is **even** — emit `0e ⊕ 1e ⊕ 2e`, **never `1o`** (`1o ⊗ 1o → 0e ⊕ 1e ⊕ 2e`
> because parity multiplies). The PCS scalar is a genuine `0e` coupling, **not** the trace of the
> `2e`. Decision: we **emit the full irrep set always** (0e ⊕ 1e ⊕ 2e), not a pure-`2e` surface —
> all the layers for all the things. MOPAC Wiberg bond order enters as a **calibrated source-strength**
> term (two basis channels: fixed source + bond-order-weighted), not as "Δχ ∝ Wiberg" physics. The
> real, cited, vetted **McConnell spec** built from that grounding supersedes the recommendation
> below; this note is kept as the reasoning trail.

Terms are defined as they appear. "Δχ" is the anisotropy of the magnetic
susceptibility tensor of the neighbouring group (units below). "Source" = the
neighbouring bond/group whose induced moment does the shielding. "Target" = the
nucleus being shielded. "T0/T1/T2" = the rank-0 (scalar), rank-1 (vector), and
rank-2 (symmetric traceless) irreducible pieces of a Cartesian 3×3 tensor;
equivalently the e3nn irreps `0e`, `1o`/`1e`, `2e`.

---

## 1. The effect, and how the literature says to compute it

### 1.1 The physical picture

Put a group with anisotropic magnetic susceptibility near a nucleus. The
external field B₀ induces a magnetic moment **m** in that group's electrons.
Because the group's susceptibility is a *tensor* (it responds more strongly
along some molecular axes than others), the induced moment is not parallel to
B₀. That misaligned induced moment is a magnetic dipole; its stray field at the
target nucleus adds to or subtracts from B₀, shifting the resonance. This is the
"neighbour anisotropy" or "magnetic anisotropy" contribution, and for aromatic
groups it is the same physics as the ring-current effect — those are the groups
with the largest susceptibility anisotropy
([Case, *Chemical shifts in biomolecules*, PMC3877577](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/);
[RSC Adv. 2014, C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b)).

The standard treatment is the **point-dipole approximation**: collapse the
group's induced currents to a point dipole at the group centre, then evaluate
its dipole field at the target. H. M. McConnell wrote this down in 1957
([J. Chem. Phys. 27, 226](https://doi.org/10.1063/1.1743676)).

### 1.2 The general (full-tensor) McConnell form

The honest, general statement — the one modern biomolecular-shift reviews use —
is a **tensor** equation, not the scalar magic-angle formula. The contribution
to the *shielding tensor* of the target from a remote group is the contraction
of the dipole propagator with the group's susceptibility tensor **χ**. Writing
**r** for the source-centre→target vector, r = |**r**|, **I** the identity, and
⊗ the outer product
([Case, PMC3877577](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)):

```
σ(neighbour)  =  (1 / 4π) · [ χ / r³  −  3 (r ⊗ r) χ / r⁵ ]
```

i.e. the dipole-field propagator **D**(**r**) = (3 **r̂**⊗**r̂** − **I**)/r³
acting on **χ**. The key facts that make this a *clean* statement:

- The **isotropic part of χ produces no net isotropic shift** — it cancels on
  orientational averaging. So the observable depends on the **anisotropy** Δχ,
  not on χ itself ([Case, PMC3877577](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)).
  This is *why* the feature is a Δχ feature.
- The propagator **D**(**r**) is exactly a rank-2, symmetric, traceless object —
  pure spherical-harmonic ℓ=2 in **r̂**. The angular content of the *clean*
  effect is therefore ℓ=2 and nothing else.

### 1.3 Where (3cos²θ−1)/r³ comes from, and what it drops

The familiar textbook form

```
σ ∝ Δχ · (1 − 3cos²θ) / (3 r³)        (chemical-shift sign; θ = angle from the group axis)
```

is **not** a different physics. It is the full-tensor form above under two
collapses, both of which throw information away:

1. **Axial symmetry of the source.** Assume the group's χ tensor is axially
   symmetric, so a single scalar Δχ = χ∥ − χ⊥ describes it (true-ish for C≡C,
   roughly for C=C/C=O along their local axes; *false* for a fully rhombic
   group). This reduces the five components of χ to one number.
2. **Isotropic (powder/tumbling) averaging.** Average over molecular
   orientations relative to B₀. This projects the rank-2 angular dependence onto
   the second Legendre polynomial P₂(cosθ) = (3cos²θ−1)/2 — the *only* survivor
   of the tensor under isotropic averaging.

The "magic angle" θ = 54.7° where the effect vanishes is just the zero of
P₂(cosθ); the same zero that kills dipolar couplings and CSA under magic-angle
spinning ([Wikipedia, *Magic angle*](https://en.wikipedia.org/wiki/Magic_angle)).
So (3cos²θ−1)/r³ is the **doubly-collapsed, orientation-averaged scalar** of a
genuinely rank-2 tensor effect. It is correct as far as it goes; it is *lossy*
by construction. For an axially-symmetric source it loses nothing beyond the
averaging; for a rhombic source it also loses the in-plane (rhombic) asymmetry.

### 1.4 Non-axial / rhombic sources, and the cleanest modern statement

When the source is not axial you need the full χ. The cleanest modern statement
of exactly this rank-2 through-space effect comes from the **pseudocontact-shift
(PCS)** literature — the paramagnetic analogue is *the same dipole-on-anisotropic-χ
maths*, just with a paramagnetic χ. There the field-standard form is written
directly in **second-rank spherical harmonics**
([Parigi, Benda et al., *Pseudocontact shifts and paramagnetic susceptibility…*,
arXiv:1804.09055](https://arxiv.org/abs/1804.09055);
[Benda et al., *Solution of a Puzzle*, PMC7584370](https://pmc.ncbi.nlm.nih.gov/articles/PMC7584370/)):

```
σ_point  =  (1 / 4π r³) · Σ_{m=-2..2}  χ_m  Y₂^m(r̂)
```

where χ_m (m = −2…+2) are the **five irreducible spherical components** of the
susceptibility tensor, and Y₂^m(**r̂**) are the second-rank spherical harmonics
of the source→target direction. Equivalently the five χ_m collapse to:

- **axiality** Δχ_ax = χ_zz − ½(χ_xx + χ_yy),
- **rhombicity** Δχ_rh = χ_xx − χ_yy,
- and **three Euler angles** orienting the χ principal-axis frame.

This is the rigorous, complete form (Kurland–McGarvey lineage,
[J. Magn. Reson. 2, 286 (1970)]). The axial-only formula is its Δχ_rh = 0,
orientation-averaged special case. **This spherical-harmonic statement is the
bridge to e3nn**, because Y₂^m *is* the ℓ=2 irrep e3nn consumes (§3).

### 1.5 The honest ceiling: the point-dipole model is itself an approximation

Two independent strands of the literature say the McConnell point-dipole picture
is a *useful semi-empirical model, not ground truth*, and it is worth carrying
this as a stated limit rather than pretending the feature is exact:

- **Spatial-extent failure.** Collapsing a delocalised current loop to a point
  is only valid when r ≫ group size. Near the group (the regime that matters for
  protein contacts at 3–5 Å) the point-dipole field is wrong; treating the
  source as a probability/current *distribution* gives materially different
  results ([*Pseudocontact shifts from mobile spin labels*, arXiv:1607.08869](https://arxiv.org/abs/1607.08869)).
  This is the classical motivation for ApSimon-style double-summation tables and
  for the GIAO "through-space NMR shielding" (TSNMR/ICSS) maps that replaced the
  single point dipole with computed iso-shielding surfaces
  ([Klod & Kleinpeter, *J. Mol. Struct. THEOCHEM*; Modgraph alkene/carbonyl notes](https://www.sciencedirect.com/science/article/abs/pii/S1093326302001201)).
- **Decomposition is not clean.** GIAO/NICS analyses (Lazzeretti and others)
  found that the conventional "π-electron anisotropy cone" for C=C and aromatic
  rings does not survive a careful orbital decomposition — much of the in-plane
  deshielding / out-of-plane shielding is actually a **σ** contribution, not the
  π-anisotropy the cone model attributes it to
  ([*Is the conventional interpretation… correct?*, Chem. Eur. J. 2011,
  PMID 22135110](https://pubmed.ncbi.nlm.nih.gov/22135110/)).

Neither strand says "don't use McConnell". They say: it is a **calibrated
geometric descriptor of a real field**, exact only in the far field and for an
axial source, and the Δχ values are fitted effective parameters, not
first-principles constants.

### 1.6 Where the Δχ values come from

Δχ for organic groups is tabulated, not derived per-structure. The canonical
sources are diamagnetic-anisotropy compilations and NMR-fit values: Pople's
point-dipole estimates, ApSimon, Schneider, Zürcher, Williamson–Asakura (the
protein-shift lineage), and Abraham's re-fits. Units are usually molar cgs
(10⁻⁶ cm³ mol⁻¹) or per-molecule (10⁻³⁰ cm³ ≈ 10⁻⁶ Å³). Representative anchor:
Δχ(C=C) ≈ −12.1 × 10⁻⁶ cm³ mol⁻¹ (perpendicular term)
([Modgraph, *Chemical shifts of alkenes*](http://www.modgraph.co.uk/Downloads/Chemical%20shifts%20of%20alkenes.pdf)).
The values **disagree across sources by factors of ~2–5** for the same group
(e.g. carbonyl C=O anisotropy spans ~+2 to +24 × 10⁻⁶ cm³ mol⁻¹ across
Williamson / Abraham / ApSimon / Schneider) — this scatter is itself the honest
statement of how well-determined the parameter is. Aromatic Δχ is the
macroscopic face of the ring current; if a ring-current kernel is already in the
feature set, an aromatic Δχ term **double-counts** the same physics and should
be zeroed.

---

## 2. The defensible, sane recommendation

### 2.1 What to compute — the physical quantity

Compute the **full-tensor McConnell contribution** per (source group, target
atom) pair and accumulate it, then expose its irreducible pieces. Concretely,
for each neighbouring group within a KD-tree cutoff:

```
σ_pair(source→target)  =  Δχ_tensor(source)  contracted with  D(r) = (3 r̂⊗r̂ − I)/r³
```

Do **not** start from (3cos²θ−1)/r³. Start from the dipole propagator **D**(**r**)
— which is exactly ℓ=2 in the geometry — and the source susceptibility tensor
Δχ(source). The scalar magic-angle form is then recoverable as the
axial-isotropic contraction, but you keep the rank-2 information you would
otherwise have thrown away. This is well within our linalg budget (a 3×3
contraction per neighbour pair over a KD-tree neighbourhood).

**Source susceptibility:** prefer the **axially-symmetric Δχ per bond category**
as the default (one Δχ per category along the local bond/group axis), because
that is what is tabulated and traceable. Build the source tensor as

```
Δχ_tensor = Δχ · ( ŝ ⊗ ŝ − I/3 )
```

where **ŝ** is the group's local symmetry axis (bond direction / carbonyl axis /
ring normal). This is honest: it is the axial model, but it enters the *tensor*
propagator, so the **rank-2 angular structure is preserved even though the
source model is axial**. Carry rhombic Δχ_rh only for groups where a tabulated
rhombic value exists (rare; most organic-group tables are axial). Where no
rhombic value exists, set Δχ_rh = 0 and say so — "I don't know the rhombicity"
is the correct state, not a fabricated number.

**Δχ values:** pull from the tabulated semi-empirical sets (Williamson–Asakura
for the peptide C=O / C−N lineage is the most defensible protein-specific
choice; Abraham for neutral sidechain amides). Treat them as **fitted effective
parameters with a stated source and a stated uncertainty band**, not constants.
The right design move is to make Δχ a *named per-category parameter the model can
calibrate*, exactly as ring-current intensities are treated — the geometry is
emitted unscaled (Å⁻³) and the Δχ scaling is a calibration layer, so the
parameter scatter in §1.6 becomes a fitted coefficient rather than a hard-coded
guess.

**Does the full tensor beat the point-dipole scalar?** Yes, for *this* use, and
the reason is specific: we are feeding an **equivariant** model. The scalar
(3cos²θ−1)/r³ is the orientation-*averaged* projection — it has thrown away
exactly the angular degrees of freedom an equivariant network exists to use. The
full propagator costs essentially nothing more and hands the network the honest
ℓ=2 field. The point-dipole *spatial* approximation (point vs distributed
source) is a separate, deeper limitation that the full tensor does **not** fix —
that one is a stated ceiling (§1.5), not something to engineer away here.

### 2.2 What form / which irreps for e3nn

A Cartesian rank-2 tensor decomposes as 3 ⊗ 3 = **1 ⊕ 3 ⊕ 5** =
`0e ⊕ 1(o/e) ⊕ 2e`
([e3nn docs](https://docs.e3nn.org/); [Tensor Field Networks, arXiv:1802.08219](https://arxiv.org/abs/1802.08219)):

- **0e** — the trace (rank-0 scalar).
- **1e** — the antisymmetric part (rank-1, a pseudovector). The rank-2
  *shielding/susceptibility* tensor is even, so its irreps are `0e ⊕ 1e ⊕ 2e`; a tensor built as an
  outer product of polar vectors is **also even** — `1o ⊗ 1o → 0e ⊕ 1e ⊕ 2e` because parity
  multiplies (odd × odd = even). So there is **no `1o`** here — emit even irreps. (Corrected
  2026-06-06; the old "outer product → `1o`" line was wrong — see the banner +
  [`mcconnell_grounding_agent1.md`](mcconnell_grounding_agent1.md).)
- **2e** — the symmetric traceless part (rank-2). **This is the physically
  clean piece** of the McConnell effect: the dipole propagator **D**(**r**) is
  exactly symmetric-traceless, so the honest, single-ℓ form of the feature is
  the **2e** spherical-tensor (the five χ_m·Y₂^m components of §1.4).

**Recommended emit (most defensible):**

1. **Primary feature — the ℓ=2 spherical tensor (`2e`, 5 numbers).** Form the
   symmetric-traceless McConnell tensor (the propagator-on-Δχ result, summed
   over neighbours per category), decompose to its five real spherical-harmonic
   components, emit as `2e`. This is the clean, honest, single-irrep statement of
   the effect. Emit it **per source category** (peptide C=O, peptide C−N,
   sidechain C=O, aromatic-if-not-double-counted) so the model can carry
   different Δχ per chemistry, and emit the unscaled Å⁻³ geometry so Δχ stays a
   calibratable coefficient.
2. **Secondary feature — the scalar (`0e`, 1 number) only if you also emit a
   non-clean (trace-bearing) variant.** A *pure* propagator-on-Δχ tensor is
   traceless, so its `0e` is ~0 and carries nothing; don't emit a fake scalar. If
   a particular construction does carry a trace, emit it as `0e` and label it as
   the non-dipole residue, not as "the McConnell scalar".
3. **Vector (`1`) — emit only if the construction is genuinely asymmetric**, and
   then *only with its parity stated*. For the clean symmetric propagator this is
   zero. Do not synthesise a `1o` from `cosθ·d̂` cross terms and call it physics
   unless you can name what asymmetry it represents.

This mirrors how serious people feed NMR rank-2 tensors to equivariant GNNs:
decompose the Cartesian tensor into its irreducible spherical pieces and learn /
feed each by irrep — exactly the recipe in
[*GNN predictions of solid-state NMR parameters from spherical tensor
decomposition*, arXiv:2412.15063](https://arxiv.org/abs/2412.15063). The
field-standard is "spherical-tensor decomposition, one irrep at a time", not
"nine raw Cartesian numbers" and not "one collapsed scalar".

### 2.3 Why not just (3cos²θ−1)/r³

Because it is the orientation-averaged scalar of a tensor we are handing to a
model whose entire purpose is to use orientation. It is defensible as a *sanity
scalar* or a baseline; it is the wrong *primary* feature for e3nn. Use it as the
`0e` baseline if you want a scalar, derived from the same tensor, clearly labeled
as the averaged projection.

---

## 3. Where the physics is clean and where it is not

**Clean:**
- The **rank-2 (ℓ=2 / 2e) angular structure** is exact and parameter-free given
  geometry. The dipole propagator (3 **r̂**⊗**r̂** − **I**)/r³ is rigorous
  classical electromagnetism. This part we can stand behind fully.
- The **spherical-harmonic form** σ ∝ Σ χ_m Y₂^m(**r̂**)/r³ is the field-standard
  rigorous statement (PCS literature) and maps one-to-one onto an e3nn `2e`
  feature.

**Not clean (state these as footnotes, not bugs):**
- **Δχ magnitude.** Tabulated values disagree by ~2–5× across sources (§1.6).
  The defensible response is to make Δχ a *calibrated per-category coefficient*
  with a reported band, not a fixed constant. The geometry is exact; the scale
  is fitted.
- **Point vs distributed source.** The point-dipole collapse is only valid in
  the far field; at protein contact distances (3–5 Å) it is in its weakest
  regime. We do not fix this; we state it. The honest options are (a) keep the
  point dipole and disclose the near-field caveat, or (b) move to a
  distributed/GIAO TSNMR treatment — which is a *different, heavier method*, off
  the table for a "standard, non-novel" feature.
- **Rhombicity.** Most organic-group Δχ tables are axial. Carrying Δχ_rh = 0 is
  honest absence, not an approximation we chose for convenience. If a rhombic
  value isn't tabulated, the rhombic irrep components are genuinely unknown.
- **Parity of the irreps (corrected 2026-06-06).** The shielding/susceptibility tensor is even, AND
  the geometric propagator built from displacement outer products is **also even** —
  `1o ⊗ 1o → 0e ⊕ 1e ⊕ 2e` (parity multiplies). So emit **`0e`, `1e`, `2e` — never `1o`** for these
  pieces; declaring `1o` would silently break equivariance. (The earlier claim that an
  outer-product-of-displacements tensor is odd was wrong — see
  [`mcconnell_grounding_agent1.md`](mcconnell_grounding_agent1.md), which cites e3nn's own parity
  rules, Geiger–Smidt, and Ben Mahmoud 2024.) Correctness point: get it right and state it.
- **Aromatic double-counting.** If a ring-current kernel is present, aromatic Δχ
  is the *same physics* and must be zeroed in the McConnell term, else the model
  sees the aromatic field twice.

**The honest one-line summary:** the geometry/angular form is exact and maps
cleanly to a single e3nn irrep (`2e`); the *magnitude* (Δχ) and the *point-source
and axial assumptions* are the soft parts, and the right design treats Δχ as a
calibrated coefficient and discloses the source-model limits rather than hiding
them.

---

## 4. Neutral note: current compute vs this recommendation

Descriptive only, no verdict. The current library
(`src/McConnellResult.{h,cpp}`) computes, per neighbouring bond, a **three-term
tensor**

```
M = 9 cosθ (d̂ ⊗ b̂)  −  3 (b̂ ⊗ b̂)  −  (3 d̂ ⊗ d̂ − I),   stored as M/r³,
```

where **d̂** is the bond-midpoint→atom unit vector and **b̂** the bond direction.
By construction the first term (9 cosθ d̂⊗b̂) is generally asymmetric and so
carries a rank-1 (T1) part; the second (−3 b̂⊗b̂) is symmetric and carries T0+T2;
the third (−(3 d̂⊗d̂ − I)) is the symmetric-traceless dipolar propagator K,
i.e. the clean T2. The code keeps K separately, computes the scalar
f = (3cos²θ−1)/r³ as the double-contraction of K with the bond direction, and
emits a 9-component full-M decomposition (`mc_shielding`), per-category
symmetrised-traceless T2 blocks (`mc_category_T2`), and the scalar category sums
(`mc_scalars`). Susceptibility scaling is **not** applied by the producer; the
geometry is emitted unscaled in Å⁻³, and a separate rediscover-side note
(`MCCONNELL_DCHI_LITERATURE.md`) holds candidate per-category Δχ values and the
molar→ppm prefactor. The recommendation in §2 differs principally in *which
tensor* is treated as the physical object: it takes the **clean symmetric-traceless
propagator contracted with an (axial) source Δχ tensor as the primary,
single-irrep (`2e`) feature**, derives the scalar as its orientation-averaged
`0e` projection, and treats Δχ as a calibrated per-category coefficient — whereas
the current path emits the full three-term M (carrying T0+T1+T2 together) as the
primary tensor alongside the separate K and scalar. The ~55% reconstruction gap
our earlier fitting hit is consistent with the primary emitted object being a
mixed-rank tensor rather than the clean ℓ=2 piece, but that is an observation to
test, not a conclusion to assert here.

---

## Sources

- McConnell, *J. Chem. Phys.* 27, 226 (1957) — original point-dipole model.
  [DOI 10.1063/1.1743676](https://doi.org/10.1063/1.1743676)
- Case, *Chemical shifts in biomolecules* — full-tensor McConnell form,
  isotropic-χ cancellation. [PMC3877577](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)
- Parigi, Benda et al., *Pseudocontact shifts and paramagnetic susceptibility…*
  — second-rank spherical-harmonic form, five χ_m, axiality/rhombicity.
  [arXiv:1804.09055](https://arxiv.org/abs/1804.09055)
- Benda et al., *Solution of a Puzzle* — QM confirms the semi-empirical PCS
  tensor form. [PMC7584370](https://pmc.ncbi.nlm.nih.gov/articles/PMC7584370/)
- *Pseudocontact shifts from mobile spin labels* — point vs distributed source.
  [arXiv:1607.08869](https://arxiv.org/abs/1607.08869)
- *Is the conventional interpretation of the anisotropic effects… correct?*,
  Chem. Eur. J. 2011 — GIAO/NICS critique of the π-anisotropy cone.
  [PMID 22135110](https://pubmed.ncbi.nlm.nih.gov/22135110/)
- Klod & Kleinpeter / Modgraph — TSNMR/ICSS through-space shielding maps,
  Δχ tabulations. [THEOCHEM](https://www.sciencedirect.com/science/article/abs/pii/S1093326302001201),
  [Modgraph alkenes](http://www.modgraph.co.uk/Downloads/Chemical%20shifts%20of%20alkenes.pdf)
- RSC Adv. 2014 — neighbour-anisotropy review, axial Δχ definition, magic angle.
  [C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b)
- e3nn docs; Tensor Field Networks (arXiv:1802.08219) — rank-2 → 0⊕1⊕2 irreps.
  [e3nn](https://docs.e3nn.org/), [arXiv:1802.08219](https://arxiv.org/abs/1802.08219)
- Ben Mahmoud et al., *GNN predictions of solid-state NMR parameters from
  spherical tensor decomposition* — the field-standard "decompose to irreps"
  recipe for NMR tensors. [arXiv:2412.15063](https://arxiv.org/abs/2412.15063)
