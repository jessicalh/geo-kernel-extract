# McConnell Neighbour Anisotropy — Feature Spec

*Shareable spec for one of the three featured kernels, and the
template for the other two. Built on `mcconnell_structured_grounding.md` (which carries the full
cites and the opus-adversarial pass) and `mcconnell_grounding_agent1.md`; decisions are settled in
`CONTROLLING_SPEC.md`. Figures are marked for generation.*

## 1. The effect, and why it's lovely

A protein nucleus is never alone. A few ångström away sits a peptide carbonyl, a bond, a small cloud
of electrons. When the spectrometer's field **B₀** washes over the molecule, those neighbouring
electrons stir — they circulate, and the little current becomes a tiny magnet. Here is the part that
makes it *visible*: the neighbour's electrons do not respond equally in every direction. They are
easier to push along some molecular axes than others — they are **anisotropic** — so the tiny magnet
they become does not line up with B₀. It points slightly askew. And a tilted magnet, even a faint
one, has a stray field; that field reaches across the gap to our nucleus and nudges its resonance up
or down. McConnell wrote this down in 1957. The quiet miracle is that we can compute the nudge from
almost nothing — just *where the atoms are*.

The clean statement. Write **r** for the vector from the neighbour to our nucleus, **r̂** its
direction, *r* its length. The contribution to the nucleus's **shielding tensor** is the neighbour's
susceptibility **Q** seen through the **dipole propagator**:

```
    D(r) = (3 r̂ r̂ᵀ − I) / r³          σ = D(r) · Q
```

`D(r)` is the field of a point dipole, and it is *exactly* a rank-2, symmetric, traceless object —
pure ℓ=2 in the geometry, and nothing else. The entire angular content of the clean effect is sitting
in plain sight in the geometry. [Case 2013; Suturina & Kuprov 2016]

A 3×3 tensor decomposes, in the language an equivariant network speaks, into three irreducible pieces:

```
    σ  →  0e  ⊕  1e  ⊕  2e          (1 + 3 + 5 numbers)
```

- the **2e** is the clean angular heart — the five components of the ℓ=2 field;
- the **0e** is a genuine scalar coupling, the pseudocontact/McConnell scalar `n·Q·n / r³` — *not* the
  trace of the 2e (a pure 2e is traceless, so that would be zero);
- the **1e** is the antisymmetric remainder: because `D` and `Q` need not commute, their product
  carries a small pseudovector — the antisymmetric part of the shielding response, real if
  conventionally hard to observe.

All three are **even** under inversion — `0e/1e/2e`, **never `1o`** (a shielding tensor is even; two
polar-vector factors give `1o ⊗ 1o → 0e ⊕ 1e ⊕ 2e`, because parity multiplies). This is not pedantry:
declare the wrong parity and an equivariant model's guarantee is silently broken. [Ben Mahmoud 2024;
Geiger & Smidt 2022]

## 2. What we compute

We replace the calculator's current three-term tensor `M` — which jumbled rank-0, 1, and 2 together,
was asymmetric and non-traceless, and is *not* a defensible object to hand an equivariant model — with
the clean propagator-on-susceptibility above. Per (target nucleus *i*, neighbour source *s*):

```
    A_is = D(r_is) · Q̂_s
```

where `Q̂_s` is the **unit-strength shape** of the source's susceptibility (axial by default):

```
    Q̂_s = (û_s û_sᵀ − I/3)        ( + a rhombic term ONLY where a real value exists )
```

`û_s` is the source's local axis (bond direction / carbonyl axis). We **emit the full even irrep set
`0e ⊕ 1e ⊕ 2e`**, per source category, in **two channels**:

- a **fixed-source** channel: `Σ A_is`;
- a **bond-order** channel: `Σ (bond_order_s) · A_is`, weighting each source by its measured **MOPAC
  Wiberg bond order** (a rotationally-invariant, QM-derived measure of electron sharing).

A model or fit reads the fixed and bond-order channels per category, with susceptibility scales pinned
from literature/physics rather than learned as free coefficients. **Δχ is not fitted from this sample:**
the geometry is emitted as a unit shape (Å⁻³) and scaled by defended values. Aromatic
sources are kept as a category but set to **zero** while the ring-current kernels are active — that
anisotropy is the *same physics* the ring kernels already carry, and we do not count it twice.

The MOPAC bond order is an honest claim: a measured source-strength modulator beside the fixed channel,
**not** a law that susceptibility is linear in bond order — we found no citation that proves such a
law, so we do not assert one. [OpenMOPAC]

## 3. The structural biology that earns its place

A protein is built from exactly the groups this effect cares about: peptide bonds, backbone and
sidechain carbonyls, sidechain unsaturation, aromatic rings. Each responds anisotropically to the
field, and the literature has measured that response at protein scale — peptide planes carry a real
diamagnetic anisotropy [Worcester 1978; Pauling 1979], and a whole protein's magnetic-susceptibility
tensor can be reconstructed from peptide-bond and aromatic-sidechain contributions [Tjandra 1996].
Case's biomolecular-shift work makes the neighbour-anisotropy term a standard ingredient of any
honest shielding decomposition [Case 2013]. The ring's share is handled by
Biot–Savart and Haigh–Mallion; **this kernel is the non-aromatic bond/group anisotropy** — the rest of
the protein's through-space magnetic story.

## 4. What is problematic — said before the examiner says it

None of these is a reason to drop the feature. They are the limits we name first.

1. **The susceptibility magnitudes are not settled constants.** Tabulated Δχ disagree by factors of
   two to five across sources for the same group; our held sources already show it (Pauling rejects an earlier
   ester-based value as unreliable and gives a corrected peptide value differing in *sign* and size). The defensible response is
   the one we took: emit a **calibrated source-shape basis**, not a fixed Δχ table, and let DFT fix the
   scale. *(We do not yet hold primary-verified Δχ numbers — the spec defends shapes, not constants. An
   honest gap, named.)*
2. **Rhombicity is usually unknown.** Most organic-group tables are axial. We carry a rhombic term only
   where a real value exists, and mark its absence in metadata — "I don't know the rhombicity" is the
   correct state, not a fabricated number. [Suturina & Kuprov 2016]
3. **The source is not really a point.** The point-dipole collapse is a far-field approximation —
   Suturina shows distributed-source corrections exceeding 1 ppm out to **~10 Å**, well past the 3–5 Å
   of a close protein contact. We disclose this; we do not engineer it away. [Suturina & Kuprov 2016]
4. **The clean cone is not clean at the orbital level.** GIAO/NICS decompositions show the textbook
   π-anisotropy "cone" does not survive a clean orbital decomposition — the shielding is set by σ/π
   and tensor-component balances, not by a single π-current label. [Baranac-Stojanović 2012]

## 5. A picture showing why — the danger zone

> **[FIGURE 1 — the propagator butterfly.]** The ℓ=2 angular pattern `(3cos²θ−1)/r³` around a source:
> the two-lobe shielding/deshielding "butterfly", with the magic-angle cone (54.7°) where it vanishes.
> *Caption: the clean angular heart of the effect — and the surface the model actually sees.*

> **[FIGURE 2 — the danger zone.]** McConnell's predicted shielding above a C=C vs the GIAO-SCF truth
> as a proton approaches: the two agree in the far field and **split, even in sign, below ~3 Å.**
> *Caption: where the point dipole lies. Built from Martin & Brown 2000.*

Below about 3 Å, the point-dipole picture can get the sign **wrong**. Martin & Brown put the McConnell
shielding-cone next to GIAO-SCF and experiment for a proton close over a C=C: McConnell predicts
shielding where GIAO and experiment show strong *de*shielding — because the near field also carries
electric-field, orbital, and dispersion terms the magnetic model omits. Protein contacts live in this
regime. The honest build response is not to hide it behind a cleaner tensor — we **emit the
point-dipole basis, record the near-field distances and exclusions, calibrate against DFT, and state
plainly that below ~3 Å this is a semi-empirical far-field descriptor.** The irrep correction makes the
feature *equivariant and defensible*; it does **not** make the point-dipole *exact*.

## 6. The fields produced, and what the model does with them

Per source category `c` ∈ {`peptide_co`, `peptide_cn`, `backbone_other`, `sidechain_co`,
`sidechain_other`, `disulfide`, `aromatic_zeroed` (zeroed while ring active)}, two channels (`fixed`, `bo`),
each carrying the full even irrep set:

```
    mc_<c>_0e_<ch>   (N,1)   PCS scalar coupling
    mc_<c>_1e_<ch>   (N,3)   antisymmetric pseudovector
    mc_<c>_2e_<ch>   (N,5)   ℓ=2 symmetric-traceless (real spherical)
```

stored as packed `(N,9)` arrays per (category, channel), units Å⁻³, parity declared **even**. The
aromatic-zeroed rows are *kept* (as zeros, while ring is active) so the no-double-count decision stays
auditable in the feature chart. [The full row list and the `SphericalTensor` component order are in
`mcconnell_structured_grounding.md` → "Fields Produced".]

These arrays are emitted for diagnostics and ablations. They are not Step-1 inputs. Step 2 and Step 3
may use them only if ablation earns them.

---

### Build notes (the adversarial flags, folded)

- The `2e` real-spherical normalization constants must be **verified against the actual
  `SphericalTensor::Decompose`**, not trusted from any transcription — the rotation-equivariance test is
  the backstop.
- The `1e` is the genuine **antisymmetric part** of the shielding response (real, conventionally
  near-unobservable); we carry it per full-irrep-always and say what it is rather than dress it.
- The form leans on the **pseudocontact-shift** literature (Suturina & Kuprov) for the clean spherical
  math — valid because PCS is the *same* dipole-on-anisotropic-χ maths; Case is the diamagnetic
  biomolecular anchor.
- The seven-category scheme is an **engineering/topology choice**, not literature-derived — stated.
- This spec defends the **form**; whether the clean tensor recovers *more signal* than the old `M` is
  a later empirical question, measured later, not claimed here.
- The full build-time numeric checks — rotation-equivariance, parity-declaration, PCS-scalar
  consistency (`T0(D·Q) = trace/3` vs `n·Q·n/r³`), aromatic-zero, MOPAC-channel separation, and the
  near-field audit — are enumerated in `mcconnell_structured_grounding.md` → "Implementation Checks".
  Build against that list, not this prose.

### Cites
Case 2013 (PMC3877577); Suturina & Kuprov 2016; Ben Mahmoud 2024; Geiger & Smidt 2022; OpenMOPAC
BONDS/ALLBONDS; Worcester 1978; Pauling 1979; Tjandra 1996; Martin & Brown 2000; Baranac-Stojanović
2012. *(Full chunk/DOI links in `mcconnell_structured_grounding.md`. Buckingham 1960 — the
electric-field σ_E term — belongs to the charge/EFG kernel, not here, and is cited there.)*
