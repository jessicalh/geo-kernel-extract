# Kernel design: the hydrogen-bond contribution to shielding

Scope: kept/emitted caged H-bond and Larsen artifacts for the effect a hydrogen
bond has on the magnetic shielding of the nuclei it involves, principally the
donor proton (amide HN, and Hα as a weak donor). They are not Step-1 inputs. This
document answers what the emitted descriptors represent, grounded in the open
literature first; our code is summarised neutrally at the end.

This kernel is the odd sibling. The other through-space kernels (ring current,
McConnell, charge/EFG) are *fields* whose primary object is a rank-2 tensor with
a physically interesting `2e` part. The dominant, best-established H-bond term is
different: it is an **empirical scalar shielding contribution in ppm** — a `0e`
object — read off a calibrated table of H-bond geometry. There *is* a real
anisotropic piece (the H-bond changes the proton CSA), but it is secondary, less
cleanly separable, and partly the same physics the charge/EFG kernel already
carries. The design below keeps the scalar as the spine and treats the
anisotropy and the geometry descriptors as honest, clearly-labelled additions.

Terms used throughout:
- **Donor / acceptor**: in an X–H···Y hydrogen bond, X–H is the donor (here
  almost always N–H, sometimes Cα–H) and Y is the acceptor (here a carbonyl,
  hydroxyl, or carboxylate **oxygen**).
- **H-bond geometry coordinates**: the H···O distance r (sometimes the heavy-atom
  N···O distance), the angle at the acceptor (∠H···O=C), and a dihedral about
  the acceptor (∠H···O=C–N etc.). These are the three knobs the empirical tables
  are built on.
- **Shielding tensor** σ, **isotropic part** σ_iso = (1/3)Tr σ, **CSA /
  anisotropy** = everything beyond the trace. As elsewhere: "the H-bond's
  contribution to σ" means the *additive piece* attributable to the H-bond, not
  a standalone observable.
- **Irrep / l**: an irreducible representation of O(3), labelled by l (l=0
  scalar `0e`; l=1 vector `1o`/`1e`; l=2 the 5-component quadrupole-like object
  `2e`). If consumed by an equivariant model, features are written in these.

---

## 1. What the effect is, and how the literature computes it

### 1a. The dominant term is an empirical scalar in ppm

When an amide N–H donates a hydrogen bond to an acceptor oxygen, the donor proton
**deshields** — its isotropic shift moves downfield, and the size of the move
tracks the H-bond geometry: shorter H···O distance and more linear geometry give
a larger deshielding. The mechanism, as the field now reads it, is a depletion of
electron density around the proton (a Pauli-principle / charge-transfer effect)
rather than a long-range field acting on a polarizable nucleus
([Jensen et al., *J. Phys. Chem. Lett.* 2018 — 1H-shielding of H-bond donors
reflects the Pauli principle](https://pubs.acs.org/doi/10.1021/acs.jpclett.8b01502)).
The relationship is monotone but **not globally linear**: locally linear over a
narrow distance window, curved (closer to exponential) over the full range
([same](https://pubs.acs.org/doi/10.1021/acs.jpclett.8b01502)).

The protein-prediction line turns this into a tabulated, additive **scalar**
correction. The canonical machinery is **Barfield (2002)**, who computed amide-1H
shifts on formamide/NMA H-bond model dimers by DFT and fitted the result to the
H-bond distance and angle
([Barfield, *J. Am. Chem. Soc.* 2002, 124, 4158](https://pubs.acs.org/doi/10.1021/ja903772t)
is the CamShift line; the original Barfield fit underlies the QM predictors
below). **Larsen et al. (ProCS15, 2013/2015)** — the reference our calculator is
named for — make this explicit: the H-bond contribution is a lookup over
geometry,

```
Δσ_HB(r_HO, θ, ρ)      (amide-H donor)
Δσ_HαB(r_HαO, θ, ρ)    (Cα-H donor)
```

with an **exponential dependence on the H···O bond length**, an angular
dependence at the acceptor, a dihedral, and **a separate table per acceptor
oxygen type** (carbonyl, hydroxyl, carboxylate). Donor protons with no protein
H-bond partner get a fixed **Δσ_w = 2.07 ppm** water-exposure offset
([Christensen/Larsen et al., *PLOS ONE* 2013 — amide-proton shifts from QM](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877219/);
ProCS15: [Larsen et al., *PeerJ* 2015, DOI 10.7717/peerj.1344]).
Every serious empirical protein-shift predictor carries an H-bond term of this
shape — SHIFTX/SHIFTX2 has explicit "ring current, electric field and hydrogen
bonding" terms; SPARTA corrects HN/Hα secondary shifts using H-bond length;
CamShift folds H-bond distance into polynomial features (dHB², dHB⁻¹, dHB⁻², dHB⁻³)
([SHIFTX2, *J. Biomol. NMR* 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3085061/);
[Kohlhoff et al., CamShift, *JACS* 2009](https://pubs.acs.org/doi/10.1021/ja903772t)).

The unambiguous statement: **the established, well-calibrated H-bond effect on
shielding is an isotropic scalar Δσ(geometry) in ppm.** It is `0e`. It is also,
already, a *calibration* in the project's sense — a fitted law with a coefficient
— not an unscaled geometric kernel; the ProCS15 grid output *is* the ppm answer.

### 1b. There is a genuine anisotropic part, oriented by the N–H bond

The scalar is the trace; the H-bond also changes the *shape* of the proton
shielding tensor, and that change is real and measured. Solution and solid-state
NMR plus DFT on NMA···water find the amide-¹H CSA span Δσ ≈ 9–10 ppm, and it
**grows with H-bond strength**: ≈ 8.2 ppm (α-helix) vs ≈ 10.3 ppm (β-sheet),
because shorter H-bonds deshield more
([Yao, Grishaev, Bax et al., *JACS* 2010 — impact of H-bonding on amide ¹H CSA](https://pmc.ncbi.nlm.nih.gov/articles/PMC2915638/)).
Crucially, the change is **anisotropic in a frame fixed to the N–H bond**: the two
tensor components roughly *perpendicular* to N–H (σ_ZZ off the peptide plane,
σ_YY in-plane) are the ones sensitive to H-bond length and angle, while the
component roughly *parallel* to N–H (σ_XX) is nearly insensitive. The fitted
exponential distance-sensitivities differ by component (≈ 3.0, 2.8, 1.1 Å⁻¹ for
ZZ, YY, XX) — the perpendicular directions respond ~3× more strongly
([same JACS 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2915638/)).

So the H-bond contribution is **not purely `0e`.** It has a real `2e` part whose
orientation is set by the **N–H bond axis** (and secondarily the peptide-plane
normal and the H-bond direction) — *not* a freely tumbling tensor. The mechanism
the same authors give is electronic: the H-bond polarizes electron density along
N–H, perturbing the induced field perpendicular to N–H more than parallel — an
*electronic-polarization* picture, not cleanly an external-field-gradient one.

### 1c. The electric-field-of-the-H-bond view — and the overlap to flag

There is a second, older way to frame the same proton response: as a **linear
electric-field effect (LEFE)** in the Buckingham expansion. A field component
along the X–H bond draws electron density off the proton and deshields it,
Δσ ≈ −A·E_∥ − B·E² ([Buckingham 1960, *Can. J. Chem.*](https://cdnsciencepub.com/doi/10.1139/v60-040);
[de Dios/Oldfield charge-field program]). The H-bond acceptor's charge produces
exactly such a field along N–H, so part of the H-bond deshielding *can* be
written as the electric-field effect of the acceptor.

**This is the overlap to flag, and it is partial, not total:**
- The **charge/EFG kernel already computes the through-space electric field E
  (`1o`) and EFG (`2e`)** from surrounding point charges at every nucleus
  (see `charge_efg.md`). The component of the H-bond effect that *is* the Coulomb
  field of the acceptor's charges is **inside that kernel already**. Re-adding a
  field-derived H-bond term on top of the charge/EFG field would double-count it.
- But the empirical/Barfield/Larsen H-bond term is **not** just the point-charge
  field. It is fitted to *DFT on the bonded complex*, so it captures the
  short-range **charge-transfer / Pauli / covalent** part of the interaction —
  electron density actually moving off the proton — which a sum over external
  point charges does **not** represent (the proton sits *inside* the polarized
  density, not in an external field of it). The literature is explicit that the
  donor-proton deshielding is dominated by this local depletion, with acceptor
  lone-pair/bond contributions to the proton site negligible
  ([Jensen et al. 2018](https://pubs.acs.org/doi/10.1021/acs.jpclett.8b01502)).

So: the H-bond term and the charge/EFG term **share their long-range
electrostatic component and differ in their short-range covalent component.**
They are neither independent nor redundant. The honest design does not add them
as if independent and does not assume one subsumes the other; it emits both and
lets the fit weigh the shared vs unique content (and discloses the overlap).

---

## 2. Emitted equivariant descriptors

Equivariant models consume features as irreps (`Nx0e + Mx1o + Kx2e`), each a
definite-l, definite-parity object transforming by the Wigner-D matrices
([e3nn irreps guide](https://docs.e3nn.org/en/stable/guide/irreps.html)).
Physical scalars, vectors and tensors can all be represented as steerable
node/edge features ([SEGNN, Brandstetter et al., ICLR 2022](https://arxiv.org/abs/2110.02905)).

The H-bond effect decomposes into three honestly-different pieces:

1. **The scalar shielding Δσ (`0e`).** The Larsen/Barfield ppm value is rotation-
   invariant by construction — it is already a number, the trace of the
   contribution. It is a clean `0e` node feature on the affected nucleus. This is
   the spine of the feature and the only part with decades of calibration behind it.

2. **The geometric H-bond descriptors.** The quantities the law is *built on* are
   natural equivariant descriptors; if consumed, they let a model learn its own
   geometry response rather than only using the pre-baked ppm:
   - H···O distance r → `0e` (a scalar; or several `0e` powers r⁻¹, r⁻², r⁻³, the
     CamShift-style basis).
   - The H-bond **direction** (donor-H → acceptor-O unit vector) → a `1o` edge/node
     vector, raw in the molecular frame. The network forms the angle invariants
     (with the N–H axis, the C=O axis) by tensor product — do **not** pre-reduce
     to ∠ and feed a cosine, which imposes a frame and discards equivariant content.
   - The N–H bond direction is already available to the network as edge geometry;
     no special handling needed beyond emitting the H-bond direction alongside it.

3. **The anisotropic part (`2e`), as a parallel channel, clearly labelled
   uncertain.** The H-bond's change to the proton CSA is a symmetric rank-2 object
   oriented by the N–H frame (§1b). Where we can build it — i.e. where the
   producer carries a *tensor* readout of the H-bond contribution rather than only
   the scalar — keep its symmetric part and decompose `1x0e + 1x2e`, emitting the
   `2e` block. Parity is even (a shielding tensor is inversion-even), so `2e` is
   correct. This is defensible but is the **soft** part of the feature: see §4.

The dominant, defensible object is therefore **`1x0e` (the scalar) plus the
geometric descriptors (`0e` distance powers + `1o` H-bond direction)**, with the
`2e` anisotropy as an honest parallel channel rather than the headline. This is
the genuine physical departure from the ring/McConnell/EFG siblings: those are
`2e`-led because the physics *is* a through-space tensor field; this one is
`0e`-led because the established physics is an empirical scalar law, with the
tensor part real but secondary.

---

## 3. The defensible, sane recommendation

**Emit the H-bond contribution as a `0e` scalar plus geometric descriptors, with
the anisotropy and the literature ppm as labelled parallel channels:**

1. **Primary: the geometric descriptors, raw and equivariant.** Per donor proton
   (amide HN; Hα as a weak second donor), for its identified acceptor(s) within a
   sane cutoff:
   - **H-bond direction** → `1o` (donor-H → acceptor-O unit vector, molecular
     frame, unprojected).
   - **H···O distance** → `0e`, emitted unscaled (Å) and, optionally, as the
     CamShift power basis (r⁻¹, r⁻², r⁻³) so a future model has easy access to the
     forms the empirical laws use.
   - **Acceptor type** as a categorical/`0e` channel (carbonyl / hydroxyl /
     carboxylate), since the empirical tables are per-acceptor-type.
   This is the law's own input, emitted geometry-unscaled, leaving a future model
   able to learn the H-bond response and form its own angle invariants.

2. **Carry the literature scalar Δσ(geometry) as a `0e` channel / cross-check.**
   The Barfield/Larsen ppm value (per acceptor type, including the 2.07 ppm
   water-exposure offset for unbonded amide protons) is a calibrated `0e` feature.
   It is the established law made explicit; emit it as a scalar channel *and*
   keep it as the cross-check that the geometric descriptors reproduce the known
   law. Because it is already calibrated ppm, it is a calibration, not an unscaled
   kernel — label it as such.

3. **Anisotropy as a parallel `2e` channel, only where a tensor readout exists,
   labelled uncertain.** If the producer carries a tensor form of the H-bond
   contribution (oriented by the N–H frame), emit `1x0e + 1x2e` so the `2e` CSA
   change is available. Treat it as a hypothesis channel future ablation can use
   or drop, not as a settled object — its separability from the rest of the proton
   CSA is the soft part (§4).

4. **Disclose and do not double-count the charge/EFG overlap.** The long-range
   electrostatic part of the H-bond effect is already in the charge/EFG kernel's
   `1o`/`2e` field at the proton. Do **not** also add a *field-derived* H-bond
   term — the H-bond channel's distinct content is the **short-range
   covalent/Pauli** part the point-charge field misses, carried by the
   geometry descriptors and the Barfield/Larsen scalar. Emit both the charge/EFG
   field and the H-bond geometry/scalar; let future ablation resolve the shared
   component; say in the feature's documentation that they overlap.

**Why this and not the alternatives:**
- *A `2e`-led tensor kernel (treating H-bond like McConnell/dipolar)*: the
  dominant established physics is a scalar empirical law, not a through-space
  dipolar tensor field. Leading with a dipolar-form `2e` would assert a tensor
  shape the literature does not support as the primary effect for the donor proton.
- *Only the pre-baked ppm scalar*: fine as a calibration channel, but using only
  the cooked number throws away the geometry a future model could use, and bakes
  one paper's fit in as ground truth.
- *Re-deriving the H-bond effect from the electric field*: that part is the
  charge/EFG kernel's job and double-counts; the H-bond channel earns its place by
  the short-range covalent part, which the field does not carry.

Fast linalg where it buys defensibility: a KD-tree finds donor→acceptor-O pairs
within the cutoff; the direction vector and distance are trivial; the acceptor-type
class is a typed-model query; the Cartesian→`2e` projection (where a tensor exists)
is a fixed linear map. Nothing here is novel.

---

## 4. Where the physics is clean, and where it is not

**Clean:**
- **The scalar is real and well-calibrated.** Δσ(geometry) in ppm for the donor
  proton is among the best-established empirical shift terms in protein NMR
  (Barfield → SHIFTX/SPARTA/CamShift → ProCS15). As a `0e` feature it is solid.
- **The geometry is exact.** H···O distance, the acceptor angle, the H-bond
  direction vector are unambiguous given coordinates; the `1o` direction and `0e`
  distance descriptors are clean and equivariant.
- **The acceptor-type stratification is principled**, not ad hoc: the empirical
  tables genuinely differ by carbonyl / hydroxyl / carboxylate.

**Not clean (honest uncertainties):**
- **The anisotropy is the soft part.** That the H-bond changes the proton CSA is
  measured and oriented by the N–H frame (the `2e` exists and is not isotropic).
  But *separating* the H-bond's contribution to the proton tensor from the rest of
  that tensor is a model decomposition, resting on DFT difference calculations on
  model complexes, not a parameter-free split. The `2e` channel is defensible as a
  hypothesis, weaker than the `0e` scalar.
- **The overlap with charge/EFG is partial and not crisply partitioned.** "How
  much of the H-bond effect is the acceptor's electric field (already in the
  charge/EFG kernel) versus short-range covalent/Pauli depletion (unique to this
  channel)" has no clean numeric answer at this scale — it depends on the charge
  set and the H-bond. The defensible move is to emit both and disclose the
  overlap, not to claim a partition.
- **The empirical law is a fit, with model-dependent calibration.** ProCS15 is
  OPBE/6-31G(d,p) on capped model systems with continuum solvation, and its ppm is
  *that* method's answer — useful as a calibrated cross-check, not first-principles
  truth, and explicitly **not** the project's r²SCAN target (do not fold ProCS15
  ppm into the DFT-target definition).
- **The water-exposure term (2.07 ppm) is a coarse stand-in** for a real
  solvent-averaged effect our fixed-geometry, single-structure features cannot
  resolve; it is a sane default for unbonded amide protons, not a measurement.

If forced to name the single most defensible package: **the geometric H-bond
descriptors (`1o` direction + `0e` distance powers + acceptor-type) with the
Barfield/Larsen scalar as a `0e` calibration channel**, the `2e` CSA change as a
labelled parallel hypothesis where a tensor readout exists, and the charge/EFG
overlap disclosed so the field part is not double-counted. The `0e` scalar is the
rock-solid core; the `2e` is the defensible-but-soft extension — and that is the
honest trade, the mirror image of the ring/McConnell kernels where the `2e` was
the headline and the `0e` the anchor.

---

## 5. Neutral note: current code vs this recommendation

Descriptive only, no verdict. Two parallel H-bond objects exist in the library,
and a rediscover grounding note already separates them.

`src/LarsenHBondShieldingResult.{h,cpp}` implements the **ProCS15 / Barfield
scalar law**: it enumerates geometric donor→acceptor-O pairs over a spatial index
(amide HN and any Hα as donors; backbone-carbonyl / sidechain-carbonyl-as-backbone
/ hydroxyl / carboxylate acceptor classes), computes Larsen's `(r_HO, θ_HOC,
ρ_HOCX)` geometry, looks up the per-class DFT grids, and accumulates per-atom
contributions following Larsen Table 2 (which target atom gets which of the
1°/2° HB and HαB terms). It applies the 2.07 ppm water term to amide protons with
no geometric partner. The grids return **ppm**, so this calculator is the
project's documented "ppm exception" — calibrated shielding directly, not an
unscaled geometric kernel. Although the named law is isotropic, the parsed grids
carry full tensors, so the producer stores and emits each per-class contribution
as a 9-component `T0,T1,T2` `SphericalTensor` in ppm (rotated into the protein
lab frame), plus the scalar water term and a pair count.

`src/HBondResult.{h,cpp}` is a **different** object: a McConnell/dipolar-form
geometric kernel `M_ab/r³` with the bond direction replaced by the donor→acceptor
H-bond direction, over DSSP/Kabsch-Sander-resolved backbone H-bonds, decomposed to
T0/T1/T2 (it can carry all three). This is a through-space dipolar-tensor treatment
of the H-bond geometry, not the empirical scalar law.

The rediscover grounding note (`h5-reader/src/rediscover/LARSEN_GROUNDING.md`)
already records the key design point this document reaches independently: the named
Larsen law is an **isotropic ppm scalar** over H-bond geometry, the library grid
output is *already* the calibrated answer (so consuming it directly is circular for
law recovery), and a recovery substrate must emit the **geometric H-bond features
spine-side** (Larsen's `H···O`, `H–O–C` angle, `H–O–C–third` dihedral, plus any
generic `D–A`/`D–H–A` descriptors). It flags that the older rediscover stub treated
`larsen_hbond` as a generic `T0+T2` dipolar source-sum, which is stale relative to
the scalar producer.

Relative to §3, the descriptive gaps are: (a) *which* object is the canonical
feature — the recommendation leads with the **geometric descriptors + the
calibrated `0e` scalar**, where the library carries both a scalar-law calculator
(Larsen, ppm, tensor-stored) and a separate dipolar-kernel calculator (HBond),
without one designated as the current model input; (b) the **`2e` anisotropy** is
present in the stored Larsen tensor but its status as a labelled parallel
hypothesis channel (vs the scalar spine) is a framing the current emit does not
make explicit; and
(c) the **charge/EFG overlap** — the recommendation's disclose-and-don't-double-count
point — is not currently annotated at the artifact boundary. This paragraph is a
description, not a verdict.

---

### Sources
- [Jensen et al., 1H-shielding of H-bond donors reflects the Pauli principle — *J. Phys. Chem. Lett.* 2018](https://pubs.acs.org/doi/10.1021/acs.jpclett.8b01502)
- [Christensen/Larsen et al., amide-proton chemical shifts from quantum mechanics — *PLOS ONE* 2013 / PMC3877219](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877219/)
- Larsen et al., ProCS15 — *PeerJ* 2015, DOI 10.7717/peerj.1344 (the calculator's named reference; see `references/larsen-2015-procs15-dft-chemical-shift-predictor.pdf`)
- Barfield, interresidue couplings and donor-1H shifts in H-bond regions — *J. Am. Chem. Soc.* 2002, 124, 4158 (original DFT fit underlying the QM predictors)
- [Kohlhoff et al., CamShift — interatomic-distance shift prediction, H-bond polynomial terms — *JACS* 2009](https://pubs.acs.org/doi/10.1021/ja903772t)
- [Han et al., SHIFTX2 — additive ring-current / electric-field / H-bond terms — *J. Biomol. NMR* 2011 / PMC3085061](https://pmc.ncbi.nlm.nih.gov/articles/PMC3085061/)
- [Yao, Grishaev, Bax et al., impact of H-bonding on amide ¹H CSA (anisotropy, N–H frame) — *JACS* 2010 / PMC2915638](https://pmc.ncbi.nlm.nih.gov/articles/PMC2915638/)
- [Buckingham, chemical shift of molecules with polar groups (linear electric-field effect) — *Can. J. Chem.* 1960](https://cdnsciencepub.com/doi/10.1139/v60-040)
- [e3nn irreps guide](https://docs.e3nn.org/en/stable/guide/irreps.html)
- [Brandstetter et al., Geometric and Physical Quantities improve E(3) Equivariant Message Passing (SEGNN) — ICLR 2022](https://arxiv.org/abs/2110.02905)
- Sibling design notes: `charge_efg.md` (the overlapping electric-field/EFG kernel), `mcconnell.md`, `ring.md`; rediscover grounding `h5-reader/src/rediscover/LARSEN_GROUNDING.md`
