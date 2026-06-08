Historical -- not current truth; see doc/emerging/kernel_design/ring_spec.md and doc/emerging/CONTROLLING_SPEC.md.

# Kernel design: the aromatic ring-current contribution to shielding

Scope: one input feature for the e3nn equivariant shielding predictor — the
through-space magnetic effect of an aromatic ring on a nearby nucleus. This
document answers "what *should* we compute, and in what equivariant form,"
grounded in the open literature first. It is a design, not a verdict on our
current code; our code is summarised neutrally at the end.

Terms used throughout:
- **Shielding tensor** σ: the rank-2 object relating the secondary field a
  nucleus feels to the applied field, σ_ab = −∂B_a^sec/∂B_{0,b}. It is what
  DFT computes per nucleus; the ring current is *one additive contribution* to
  it. By convention here "the ring's contribution to σ" means that additive
  piece, not a standalone observable.
- **Isotropic part**: the trace, σ_iso = (1/3)Tr σ — the one number that
  shows up in a liquid-state chemical shift.
- **Anisotropy / CSA**: everything in σ beyond the trace; it survives in solids,
  in residual CSA, and is the orientation-dependent part a tensor model needs.
- **Irrep / l**: an irreducible representation of the rotation group SO(3),
  labelled by integer l (l=0 scalar, l=1 vector, l=2 the 5-component
  quadrupole-like object). e3nn consumes features written in these.

---

## 1. What the ring-current effect is, and how the field computes it

When an aromatic ring sits in the applied field B₀, the delocalised π electrons
circulate, and that circulation produces a secondary magnetic field in the
space around the ring. A nucleus in that space feels the secondary field added
to B₀, which shifts its resonance. This is a *through-space* effect: it depends
only on where the nucleus sits relative to the ring, not on any bond to it.
([Wikipedia: aromatic ring current](https://en.wikipedia.org/wiki/Aromatic_ring_current))

The literature offers a ladder of models, cheapest/coarsest to most complete.

### 1a. Point magnetic dipole (Pauling 1936; London; McConnell 1957)

The ring is collapsed to a single induced magnetic dipole at the ring centre,
with an *anisotropic* magnetic susceptibility — large perpendicular to the ring
plane (χ⊥, along the ring normal), small in-plane (χ∥). The dipole field of an
axially-symmetric source gives the through-space shielding at a nucleus a
distance R from the centre, at polar angle θ measured from the ring normal:

  Δσ_iso(R, θ) = (Δχ / 3) · (1 − 3cos²θ) / R³ ,   Δχ = χ∥ − χ⊥

This is McConnell's neighbouring-anisotropy equation specialised to a ring; it
is the origin of the famous (1 − 3cos²θ)/R³ angular shape and of the magic angle
θ ≈ 54.7° where the isotropic effect changes sign.
([RSC Adv. 2014, C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b);
[Chem LibreTexts: diamagnetic anisotropy](https://chem.libretexts.org/Bookshelves/Organic_Chemistry/Map%3A_Organic_Chemistry_(Bruice)/14%3A_NMR_Spectroscopy/14.08%3A_Diamagnetic_Anisotropy))

Crucial point for us: the (1 − 3cos²θ)/R³ scalar is *only the isotropic trace*
of the dipole's contribution. The point-dipole model already carries a full
rank-2 tensor — the dipolar field of a susceptibility tensor is itself a tensor
field — but the single number above is what you get *after* taking the trace
and assuming axial symmetry. The "(3cos²θ−1)/r³ is only the far-field shadow"
framing in our internal note is the same observation: it is the leading
*isotropic, far-field* term, not the whole object.

### 1b. Johnson–Bovey double current loop (1958)

The ring is modelled as two coaxial circular current loops, one above and one
below the ring plane (representing the π density lobes), each carrying half the
ring current. The secondary field at a point is obtained in closed form via
complete elliptic integrals. This is a genuine *distributed-source* Biot–Savart
field, valid much closer to the ring than the point dipole. The classic
Johnson–Bovey output is tabulated as the isotropic shift vs (ρ, z) in ring-plane
cylindrical coordinates.
([ResearchGate: Johnson-Bovey applied to aromatic amino acids](https://www.researchgate.net/publication/22232702_Application_of_ring-current_theory_based_on_the_Johnson-Bovey_equation_to_the_aromatic_amino_acids))

### 1c. Haigh–Mallion (1971, 1972)

A quantum-mechanical-derived geometric model: the shielding is a sum over ring
bonds of a signed factor S_ij built from the *areas* of triangles formed by the
nucleus's projection and adjacent ring atoms, weighted by ring-current
intensity. It is the workhorse for protein/nucleic-acid ring-current shift
prediction and is well calibrated per ring type.
([Haigh & Mallion 1972, Org. Magn. Reson.](https://onlinelibrary.wiley.com/doi/10.1002/mrc.1270040203);
calibration for biomolecules: [Case, J. Biomol. NMR 1995](https://link.springer.com/article/10.1007/BF00197633))

**Key shared limitation, stated plainly in the literature:** the Waugh–Fessenden,
Johnson–Bovey, and Haigh–Mallion analytical expressions all quantify *only the
isotropic component* of the shielding tensor. They were built to predict a
liquid-state shift, which is a trace. They do not, as published, give the
anisotropy.
([Boyd & Skrynnikov, JACS 2002, 124, 1832](https://pubmed.ncbi.nlm.nih.gov/11866578/))

### 1d. Full shielding tensor (Boyd & Skrynnikov 2002)

This is the reference for the *tensor* we actually want. Boyd & Skrynnikov
extend the classical machinery with an additional analytical expression so that
the ring current's contribution to the **full shielding tensor** — hence its
anisotropy (CSA) — can be computed, and validate it against DFT (benzene ring
present vs a non-aromatic analogue plus their ring-current term; "excellent
agreement"). They apply it to backbone amide CSA, including an N–H···aromatic
case. This is the established, peer-reviewed way to get the *tensor*, not just
the shift, from a classical ring-current treatment.
([Boyd & Skrynnikov, JACS 2002](https://pubmed.ncbi.nlm.nih.gov/11866578/))

### 1e. Full induced-current / NICS / ipsocentric (the DFT ground truth)

The most complete picture computes the actual induced current density j(r) in
the molecule (ipsocentric/CTOCD-DZ partitions it cleanly into π and σ
contributions) and integrates Biot–Savart over it to get the secondary field,
i.e. NICS as a spatial field. This is the gold standard and is what our DFT
campaign effectively contains. Two facts from this literature matter for the
design:
- Close to the ring the field is *not* a single dipole. NICS maxima sit above
  the atoms and migrate toward the centre as you rise off the plane; there is no
  one-to-one map from a NICS value back to a current distribution. The
  distributed structure is real at short range.
  ([Stanger, ChemPhysChem 2023](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cphc.202300080))
- NICS mixes π ring current with σ-framework and local contributions; "ring
  current" in the strict sense is the π part.
  ([gqcg NICS notes](https://gqcg.github.io/services/nics/))

### 1f. How far the dipole is trustworthy (the multipole question)

Treating the ring as a point dipole is a multipole expansion truncated at the
leading term. For a current loop of radius a the dipole form is accurate to
~2% only for R ≳ 5a, ~10% error around R ≈ 3a, and fails badly for R ≲ a.
Because the ring is axially symmetric, the magnetic *quadrupole* term vanishes;
the first correction is the octupole, falling as 1/R⁵. ([dipole-approximation
range, standard electrodynamics result, e.g.
[Physics LibreTexts §5.4](https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Essential_Graduate_Physics_-_Classical_Electrodynamics_(Likharev)/05%3A_Magnetism/5.04%3A_Magnetic_Dipole_Moment_and_Magnetic_Dipole_Media)])

For a benzene ring a ≈ 1.4 Å, so R ≳ 5a means R ≳ 7 Å for clean dipolar
behaviour. Many protein ring–nucleus contacts of interest sit at 3–5 Å —
squarely in the regime where the point dipole is a *coarse* approximation and a
distributed-source field (Johnson–Bovey / Haigh–Mallion / Boyd–Skrynnikov) is
the honest choice.

---

## 2. Turning a through-space field into an equivariant e3nn feature

e3nn consumes features written as irreps: lists like `Nx0e + Mx1o + Kx2e`,
each block a definite-l, definite-parity object that transforms by the Wigner-D
matrices under rotation. ([e3nn irreps guide](https://docs.e3nn.org/en/stable/guide/irreps.html))
Physical (not just geometric) tensor quantities can be fed as steerable
node/edge features — that is precisely the SEGNN result: covariant vectors and
tensors as inputs improve E(3)-equivariant message passing.
([Brandstetter et al., "Geometric and Physical Quantities improve E(3)
Equivariant Message Passing", ICLR 2022](https://arxiv.org/abs/2110.02905))

The decomposition we need is standard. A **symmetric** rank-2 tensor (which the
shielding-tensor *contribution* we care about is, once we keep the symmetric
part — only the symmetric part of σ affects NMR observables) splits as:

  symmetric 3×3  →  l=0 (trace, the isotropic part, `0e`)
                 ⊕  l=2 (5-component traceless-symmetric, the anisotropy, `2e`)

i.e. **`1x0e + 1x2e`**. A general (non-symmetric) rank-2 tensor would also carry
an l=1 (`1o`) antisymmetric block, but the physically meaningful shielding
contribution is symmetric, so l=1 is not part of the target object. This is
exactly the decomposition equivariant-shielding-tensor ML uses.
([Venetos et al., "ML Full NMR Chemical Shift Tensors of Silicon Oxides with
Equivariant GNNs," J. Phys. Chem. A 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/);
[GeqShift carbohydrate shifts, RSC Adv. 2024](https://pubs.rsc.org/en/content/articlehtml/2024/ra/d4ra03428g))
Parity: under inversion the shielding tensor is even, so even-parity irreps
(`0e`, `2e`) are correct.

So the ring-current contribution is **not** a clean single-l multipole. The
honest statement is: it is dominantly an l=0 + l=2 object (the isotropic shift
plus the dipolar/CSA anisotropy), with the l=2 part the physically interesting
one the tensor model exists to capture, and with sub-leading distributed-source
structure at short range that no rank-2 form captures exactly.

---

## 3. The defensible, sane recommendation

**Compute the ring-current contribution to the shielding tensor as a per-nucleus
symmetric rank-2 tensor using a distributed-source classical model, decompose it
into `1x0e + 1x2e`, and feed both blocks to e3nn as a steerable input feature.**
Concretely:

1. **Field model — distributed source, not point dipole.** Use a
   Johnson–Bovey double-loop (or equivalently a numerically integrated
   Haigh–Mallion / Biot–Savart over the ring) to get the secondary field, and
   build the *tensor* contribution via the Boyd–Skrynnikov (2002) full-tensor
   expression rather than only the isotropic shift. This is the established way
   to get the anisotropy at the ranges (3–5 Å) where protein nuclei actually
   sit, where the point dipole is only ~10% accurate. Cost is not the
   constraint for us; the distributed model is the defensible default and the
   point dipole becomes a documented far-field check, not the production form.

2. **Per ring type, use literature-calibrated ring-current intensities** (Phe,
   Tyr, Trp 5- and 6-membered, His; nucleic-acid bases if relevant) from the
   biomolecular calibrations (Case 1995 lineage). Keep the intensity as a named,
   externalised per-ring-type scalar — it is the one empirical knob, and it is
   exactly the kind of literature-scaled coefficient the project already treats
   as a calibration rather than a fit.

3. **Output object per nucleus.** Sum the tensor contributions of all rings
   within a generous cutoff (the field is long-ranged; ~15 Å with the dipolar
   near-field handled by the distributed model is sane), keep the symmetric
   part, and decompose:
   - `0e`: σ_iso, the isotropic ring-current shift (the classical Haigh–Mallion
     number, recovered as the trace — a built-in sanity anchor against decades
     of literature).
   - `2e`: the 5-component traceless-symmetric anisotropy (the CSA contribution
     Boyd–Skrynnikov add).
   Feed the concatenation `1x0e + 1x2e` as a steerable node feature on the
   nucleus. Because it is already in irreps and already equivariant, e3nn ingests
   it directly via tensor products — no per-atom local frame is imposed
   (consistent with the project's equivariance discipline; the raw geometry +
   tensor go in and equivariance handles rotation).

4. **Fast linalg where it buys defensibility, not novelty.** A KD-tree gives
   each nucleus its ring neighbourhood; the per-ring field/tensor is closed-form
   or a small fixed quadrature; the symmetric-traceless → l=2 projection is a
   fixed 5×9 (or Cartesian→spherical) linear map. All of this is cheap and
   exact. The spatial/tensor substrate we already have is exactly what this
   needs; nothing here is a new method.

**Why this and not the alternatives:**
- *Point dipole only* (the (1−3cos²θ)/R³ scalar): too coarse at protein contact
  ranges, and it discards the anisotropy — it would hand e3nn a single `0e`
  number when the whole reason to use a tensor model is the `2e` part.
- *NICS / full induced-current as the input feature*: this is the right *ground
  truth*, but recomputing induced current density per nucleus as an *input*
  re-does the DFT we are trying to predict; it belongs as the *target/anchor*,
  not as a hand-built classical input feature. Boyd–Skrynnikov is the
  established classical surrogate of exactly that field.
- *Higher multipoles (octupole, 1/R⁵)*: available in principle but
  over-engineering — at the ranges where the dipole fails, switching to the
  distributed-source field is the standard fix and is cleaner than a multipole
  series. Mention it only as the reason the point dipole alone is insufficient.

---

## 4. Where the physics is clean, and where it is not

**Clean:**
- The *form of the target object* is unambiguous: a symmetric rank-2 tensor,
  decomposing exactly into `1x0e + 1x2e` with even parity. No judgement call.
- The *far-field limit*: at R ≳ 7 Å the point-dipole (1−3cos²θ)/R³ tensor is a
  genuine, well-understood leading term, and the isotropic trace matches the
  classical ring-current shift literature. This is a solid sanity anchor.
- The *isotropic part at all ranges* is the best-established quantity in the
  field — Haigh–Mallion/Johnson–Bovey have been calibrated for decades.

**Not clean (honest uncertainties):**
- **The anisotropy is the weakest-pinned part.** It is exactly what the
  classical models historically omitted; only Boyd–Skrynnikov (2002) supplies a
  validated analytical tensor form, and it is one paper's construction validated
  on benzene + amide cases, not a broadly stress-tested standard. The `2e` block
  is defensible but carries more model risk than the `0e` block.
- **Short range (3–5 Å) is genuinely a regime, not a point.** Here the field is
  distributed and not exactly any rank-2 form; the dipole is ~10% off, the
  distributed models are better but still classical surrogates of the true
  induced current, and π/σ separation (which the classical "ring current" assumes
  cleanly) is itself only clean in the ipsocentric DFT picture.
- **Intensity calibration is per-ring-type empirical**, drawn from biomolecular
  fits; it is a known and reasonable input, but it is a calibration, not a
  first-principles constant.

If forced to name the single most defensible package: **a distributed-source
field (Johnson–Bovey / Haigh–Mallion) with the Boyd–Skrynnikov full-tensor
extension, emitted as `1x0e + 1x2e`.** If one wanted the absolutely-safe minimal
claim, the `0e` isotropic part alone is rock-solid; the value of the tensor
model lives in the `2e` part, which is the defensible-but-less-certain piece —
and that is the honest trade.

---

## 5. Neutral note: current code vs this recommendation

We currently compute three parallel ring objects (`src/`), all decomposed to
spherical tensors per nucleus. `BiotSavartResult` builds a Johnson–Bovey
double-loop B-field and a rank-1 kernel G_ab = −n_b B_a, decomposed to T0/T1/T2,
and its header cites Boyd & Skrynnikov 2002 eq. 3. `HaighMallionResult` builds a
surface-integral dipolar kernel contracted with the ring normal into a rank-1
G_ab = −V_a n_b, likewise T0/T1/T2 and citing the same paper. `RingSusceptibilityResult`
evaluates a closed-form McConnell-style tensor M_ab/r³ (with the three terms
9cosθ d̂⊗n̂, −3 n̂⊑n̂, −(3d̂⊗d̂−I)) plus the scalar (3cos²θ−1)/r³, also decomposed
to spherical tensors. So our existing surface keeps the full tensor (T0/T1/T2),
already cites the Boyd–Skrynnikov tensor reference, uses KD-tree ring
neighbourhoods and near-field/topological filters, and carries both a distributed
(Johnson–Bovey, Haigh–Mallion) and a point-dipole-family (susceptibility) path.
Relative to the recommendation above, the descriptive gap is in *which object is
treated as the production tensor* and *whether the l=0/l=2 split is the explicit
emitted feature*: the recommendation names the distributed-source field with the
Boyd–Skrynnikov full-tensor construction as the production form and the symmetric
`1x0e + 1x2e` decomposition as the e3nn feature, where our current code emits
several rank-1/full-tensor variants in parallel without one being designated the
canonical equivariant input. This is a description, not a verdict.

---

### Sources
- [Aromatic ring current — Wikipedia](https://en.wikipedia.org/wiki/Aromatic_ring_current)
- [McConnell point-dipole / neighbour anisotropy — RSC Adv. 2014, C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b)
- [Diamagnetic anisotropy — Chem LibreTexts 14.8](https://chem.libretexts.org/Bookshelves/Organic_Chemistry/Map%3A_Organic_Chemistry_(Bruice)/14%3A_NMR_Spectroscopy/14.08%3A_Diamagnetic_Anisotropy)
- [Johnson–Bovey applied to aromatic amino acids — ResearchGate](https://www.researchgate.net/publication/22232702_Application_of_ring-current_theory_based_on_the_Johnson-Bovey_equation_to_the_aromatic_amino_acids)
- [Haigh & Mallion, new ring-current tables — Org. Magn. Reson. 1972](https://onlinelibrary.wiley.com/doi/10.1002/mrc.1270040203)
- [Case, calibration of ring-current effects in proteins and nucleic acids — J. Biomol. NMR 1995](https://link.springer.com/article/10.1007/BF00197633)
- [Boyd & Skrynnikov, ring-current contribution to chemical shielding anisotropy (full tensor) — JACS 2002, 124, 1832](https://pubmed.ncbi.nlm.nih.gov/11866578/)
- [Stanger, NICS at small distances from the molecular plane — ChemPhysChem 2023](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cphc.202300080)
- [NICS / induced-current background — gqcg](https://gqcg.github.io/services/nics/)
- [Dipole-approximation range / magnetic multipoles — Physics LibreTexts 5.4](https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Essential_Graduate_Physics_-_Classical_Electrodynamics_(Likharev)/05%3A_Magnetism/5.04%3A_Magnetic_Dipole_Moment_and_Magnetic_Dipole_Media)
- [e3nn irreps guide](https://docs.e3nn.org/en/stable/guide/irreps.html)
- [Brandstetter et al., Geometric and Physical Quantities improve E(3) Equivariant Message Passing (SEGNN) — ICLR 2022](https://arxiv.org/abs/2110.02905)
- [Venetos et al., ML full NMR shift tensors with equivariant GNNs — J. Phys. Chem. A 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/)
- [GeqShift, E(3)-equivariant carbohydrate NMR shift prediction — RSC Adv. 2024](https://pubs.rsc.org/en/content/articlehtml/2024/ra/d4ra03428g)
