# Positional / equivariant encoding — review + actionables (2026-06-07)

**Provenance:** a ChatGPT 5.5 Pro summary (a textbook tour of geometric encoding for equivariant models)
→ an Opus-max **blind** code review of this project's geometric encoding (no docs, no code edits) → this
synthesis. The striking thing first: **blind, with no docs, the review independently reached most of the
conclusions this project bled for** — parity (`1e`-not-`1o`), the don't-double-feed circularity,
kernels-as-ℓ≤2-multipoles, the radial-as-the-fittable-law. That is the same two-path comfort the field
survey and the lit-sweep gave: an uninvested path landing where we landed. A meta-instance of the very
principle the project rests on.

## The ChatGPT summary

Faithful textbook tour — relative positions kill translation, `|r|` invariant / `r̂` equivariant, the
e3nn/TFN convolution stated correctly, the `rrᵀ − |r|²/3 ↔ ℓ=2` Cartesian↔spherical bridge right, the
three-level taxonomy (invariant → scalar-vector EGNN/PaiNN → full steerable TFN/NequIP/MACE/Allegro) the
right map. **Its one fatal blind spot is parity.** For a magnetic-response tensor that is the whole
ballgame: following the generic summary walks a builder straight into the `1o` bug (feed the antisymmetric
part as polar `1o`, silently break reflection equivariance). Our "never `1o`" rule is therefore **a
theorem, not a convention** — and the reality-check found that exact ghost living in our own SDK/H5 labels.

## What the blind review got right (and corroborates us)

- **Parity as a derived result.** Shielding is `B_induced = −σ·B_external`, an axial→axial response, so `σ`
  is even: `0e ⊕ 1e ⊕ 2e`; the antisymmetric (rank-1) part is an axial vector → `1e`, never `1o`.
- **Don't double-feed the ℓ≤2 kernels into an ℓ≤2 first conv** — a steerable conv's first layer computes
  `Y^(ℓ≤2)(r̂)`, the same object the kernels are; feed them in and the network trivially relearns them and
  you lose attribution. (Our circularity / dual-nature point, reached blind.)
- **Kernels are the low-order multipole moments of the weighted neighbourhood.**
- **Radial-as-fittable-law; hyperparameters-as-physics-claims** (below).

## New insights the review adds

- **`T1` is intrinsically a ≥2-body geometric quantity — the best single idea.** One edge gives `r̂⊗r̂`,
  symmetric, no antisymmetric part. The `1e` needs *two non-collinear directions* (their cross-product is
  the pseudovector). That is exactly where *our* McConnell `1e` comes from: `D(r)·Q̂` is two symmetric
  tensors that don't commute, and the non-commutation is the **`r̂ × û` cross-term** (edge direction ×
  source-bond axis). The review reverse-engineered the provenance of our `1e`. **Extension (dragon
  logic):** the `1e` is the *least* observable piece of all — antisymmetric shielding is essentially
  unmeasurable even in the solid state. Full-irrep correctly *emits* it (the DFT tensor carries it), but
  if the `2e` averages out for the solution shift, the `1e` is gone twice over. Kept for completeness +
  DFT-anchored evidence, not because anyone measures it.
- **The mirror-reflection test** — reflect the structure, confirm the pseudo-components flip with the
  right sign. The cheap runtime guard for the silent-`1o`-no-runtime-error class. The unit test "T2 is
  sacred" actually deserves.
- **TensorNet ≈ our object model.** Its hidden state is a rank-2 Cartesian tensor decomposed into
  scalar / antisymmetric-vector / symmetric-traceless = our `T0/T1/T2` natively, **no Clebsch–Gordan /
  Wigner, few message-passing steps**. An isomorphism to our `Mat3 ⊕ SphericalTensor` duality. Given the
  fido reframe, a perfectly good *lighter* predictor to try. Caveat (the review's own): "the natural first
  thing to try given the object model," not "it'll beat e3nn." (Higher-rank escape hatch if ever needed:
  Zaverkin et al.'s irreducible-Cartesian MPNN, NeurIPS 2024.)
- **Radial seeded from the physical power laws** — the network-level version of "a fittable law is the
  only actual calibration." Seed the radial basis with `1/r³` (ring/dipolar/EFG) as fixed channels + a
  learned residual, so the learned part reads as a *correction to the textbook*, not a black box. The
  right home for the found-fit-space.
- **Hyperparameters as physics claims** — single-layer / lmax=2 / fixed-radial ≈ the classical kernels;
  each knob beyond is a physical hypothesis (more layers = beyond-first-shell range; higher tensor-product
  body-order = cooperative many-body shielding, ring stacking, H-bond networks). The architecture
  hyperparameters ARE the answer to "what does the model add over the kernels" — the thesis question — so
  they belong in the write-up as claims, not in a config as defaults. The architecture-level version of
  Part 2's internal ablation.

## The one thing it could not see — the dragon

The review's spine is "the thesis claim is the T2 angular residual, so carry ℓ=2 throughout." That is
**exactly right for Parts 1–2** — the DFT target IS the full tensor, the `2e` survives at the static
anchors, the architecture effort pays there. But the **solution-shift observable is the `0e`**; the `2e`
is orthogonal to it and averages out under tumbling (the dragon). So the dragon refines every architecture
recommendation **per-part**: full steerable / carry-ℓ=2 for the DFT-anchored *evidence* (Parts 1–2);
**not** for the solution-shift *deployment* (Part 3 = `0e`, fido, likely lower-altitude) — which is
precisely our fido reframe (equivariance earned-only-where-it-helps). Blind, the review would build the
full tower everywhere; we now know that's the right tool for the evidence and over-engineering for the
deployment. Review and dragon are complementary: it couldn't see the dragon without our docs; we couldn't
have its `T1`-corollary without its fresh eyes.

**Refinement to its don't-double-feed point:** *don't* double-feed for **Part 1**, where attribution is
the whole point (raw-geometry-only vs kernels-with-conv-ablated = the maths-methods experiment). But for
**Part 2 fido**, double-feeding is *fine* — predicting, not attributing. The confound is real only where
we make the "what did the model add" claim, which is Part 1.

## Actionables

1. **Add the mirror-reflection parity test** — the runtime guard for the convention-mismatch class.
2. **Seed the radial from `1/rⁿ` + learned residual** — calibration-as-correction; the found-fit-space's home.
3. **Write the architecture hyperparameters (lmax / depth / body-order) as physics claims**, in the Part-2
   methodology chapter.
4. **Keep the double-feed split honest:** clean (no double-feed) for Part-1 attribution; fido (double-feed
   fine) for Part-2 prediction.
5. **Put TensorNet on the "first architectures to try" list** — lighter, ℓ≤2-Cartesian, isomorphic to our
   `Mat3 ⊕ SphericalTensor` object model; an interpretability gift. (e3nn stays the reference; this is the
   lower-friction candidate.)
6. **Carry the per-part altitude:** full steerable for the DFT-anchored tensor work (Parts 1–2); not
   mandated for the solution-shift deployment (Part 3, fido).

## The shape of it

A genuinely good review whose best parts are real gifts (the `T1`-≥2-body corollary, the mirror test, the
radial seeding, TensorNet-as-object-model), and whose only gap is *exactly* the dragon. **That is the best
shape a blind review can return in: it agrees with everything defensible and misses only the thing we
ourselves discovered.** It confirms we're standing on recognised ground while *we* hold the one piece of
local knowledge it couldn't supply.

*More passes and thought are warranted — but this, in exact detail, reads as a map off the limb: standard
symmetry machinery + standard NMR physics, combined honestly, with the dragon as the one piece of local
knowledge no blind reviewer could supply.*
