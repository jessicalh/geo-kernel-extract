# Ring current — shareable spec

*The aromatic ring-current effect on nuclear shielding, by two independent paths: Biot–Savart
(Johnson–Bovey) and Haigh–Mallion.*

**Provenance.** Third of the Three to reach spec (McConnell was the template; EFG the second). Folds the
first-stage grounding (`bs_hm_grounding_agent1.md`), its structured grounding
(`bs_hm_structured_grounding.md`), and the 2026-06-07 codex confirm pass (which served as the ring
reality-check). **Register: confirm-and-refine, and the one that wants the most care** — the producer
already computes both kernels; the surgery is naming, parity metadata, a sign pass, and switching on
dropped channels, *not* new physics. **Citation note:** the Johnson–Bovey (1958) and Haigh–Mallion
(1971/1980) primaries are frayed (not held on disk); we cite *through* the held bridges — **Case 1995,
Moyna 1998, Perkins 1977, Agarwal 1977, Boyd–Skrynnikov 2002** — and the vet pass verifies the held-file
basis, the way it did for EFG.

---

## 1. The wonder, and the clean maths

Put an aromatic ring in a magnetic field and its delocalised π electrons circulate — a **ring current** —
which raises its own secondary magnetic field. A nucleus sitting *above or below* the ring sees that
secondary field oppose the applied one (shielded, upfield); a nucleus *in the ring plane* sees it add
(deshielded, downfield). This is the famous **ring-current "butterfly"** — one of the largest and most
recognisable through-space effects in protein NMR, and the reason aromatic side chains leave such a
distinctive fingerprint on nearby proton shifts.

Two classical models compute it, by genuinely different mathematics — and that is the point:

- **Biot–Savart / Johnson–Bovey.** Model the ring current as physical current loops (a double loop at
  `±lobe_offset` above and below the ring plane), integrate the **Biot–Savart** field of the wire
  segments at the target, and read off the shielding. The isotropic Johnson–Bovey shift is the classic
  result; **Boyd–Skrynnikov (2002, eq. 3)** is the tensor lift that assembles the *full rank-2* shielding
  contribution from it (not a separate calculator — the formula that turns the scalar into the tensor).
- **Haigh–Mallion.** A geometric model: a signed-area / surface integral over the ring, weighted by the
  dipolar `(3cos²θ−1)/r³` kernel. An entirely independent derivation of the *same* physical quantity.

Their agreement *by different math* is a **two-path validation** (`feedback_two_path_validation`) — but an
honest one. Where Johnson–Bovey and Haigh–Mallion concur on **sign, lobe structure, scalar shift, and
principal-tensor orientation**, we have a sanity check no single calculator gives. The check is
**strongest near the ring**, where the wire-integral and the surface-integral genuinely diverge; **far
from the ring both approach the same dipolar limit, so agreement there is partly expected, not deep
validation**. The value is the near-ring sign/topology/orientation agreement, reported as **signed
residuals**, not absolute convergence. (Held bridges: Case 1995 calibrates both against DFT, no
significant scalar-fit difference for protein aromatics; Moyna 1998 finds the two within ~5% RMS on
protein proton shifts.)

**Parity.** The ring-current shielding is a *magnetic* response tensor, so it is **even**: `0e ⊕ 1e ⊕ 2e`,
**never `1o`**. Here McConnell's all-even rule **does** apply (unlike the EFG, whose E-field is a genuine
polar `1o`) — the ring contribution is a shielding tensor, not an electrostatic input. The antisymmetric
part is an axial vector → `1e`, never `1o`.

---

## 2. What we compute

- **Biot–Savart** (`BiotSavartResult`): the double-loop Johnson–Bovey field in SI at the target, then
  `G_ab = −n_b · B_a · PPM_FACTOR` (ring normal ⊗ B-field, shielding-sign convention), decomposed to
  `T0/T1/T2`. Computed at **unit current** (`I = 1 nA`) — the geometric kernel is independent of
  intensity.
- **Haigh–Mallion** (`HaighMallionResult`): a surface integral of the dipolar kernel over a
  fan-triangulated ring area, contracted with the ring normal, then the same `G_ab = −n_b V_a` shielding
  tensor. Emitted in geometric units (`Å⁻¹`), unscaled.

Both are **source-based** kernels following the shared rule-1 pattern (see the surgery guide): **emit the
unit kernel; the physical scale rides separately.** The literature **ring-current intensity** is the
scale (`ring.Intensity()`, per ring type, cited, with a verified worked example: `I = −12 nA` for PHE,
3 Å above the ring → `+1.40 ppm` shielded; `σ = I·G`, `I < 0` diamagnetic). The code records the
intensity in provenance but **computes with `I = 1` and does not bake it into the emit** — the
calibration **learns** the coefficient (the same mirror as McConnell; the literature intensity is the
**direction-4 physical hypothesis**, not a pinned scale). The rediscover layer *does* carry a
literature-scaled `jb_T*` via `LiteratureIntensity()`, separate from the unit-current producer emit.

**The Haigh–Mallion naming caveat (carry it honestly).** Our `HaighMallionResult` is a surface-integral
*tensor* variant, **not** provably the literal published *scalar* HM signed-area bond-sum
`Σ S_ij (1/r_i³ + 1/r_j³)` (Case eq. 5 / Moyna). Same family, not proven the same formula. The in-scope
move: benchmark our tensor's **scalar part `T0`** against the published bond-sum in populated geometries —
*correlate, don't match* (up to scale/sign). If it correlates → document it as **"Haigh–Mallion extended
to the full tensor"** (the published HM is scalar-only; our `T2` is the principled extension, the name
earned). If not → rename **"HM-style surface kernel."** Keep both `hm_H` (raw symmetric-traceless surface
kernel) and `hm_G` (the shielding tensor).

---

## 3. The structural biology

Aromatic side chains — **Phe, Tyr, Trp, His** — are the sources, and their ring currents are the textbook
explanation for the dramatic upfield (and occasional downfield) proton shifts seen near aromatic rings in
folded proteins. A methyl group parked over a Phe face can shift several tenths of a ppm; a backbone
amide swung into an aromatic plane deshields. Ring-current corrections are a standard term in every
empirical protein-shift predictor, precisely because the effect is large, long-range, and unmistakably
geometric. That geometry is what makes it a clean equivariant feature: the same `G` that estimates the
secondary field also encodes the directional arrangement of the aromatic neighbourhood.

---

## 4. What's clean, and what's not

**Clean (high confidence):**
- **The geometric kernel.** Given the ring geometry and a current, the Johnson–Bovey field and the
  Haigh–Mallion surface integral are well-defined and equivariant. Rotate the protein, both rotate
  correctly.
- **The two-path topology.** Outside the pathological near-ring zone, BS and HM agree on sign, lobe
  structure, scalar `T0`, and principal `T2` orientation — the agreement that buys the sanity check.

**Less clean (state the honest range):**
- **The intensity constants.** The current defaults read like older Giessner-Prettre / Pullman-style
  factors scaled to Phe `= −12 nA/T` (`CalculatorConfig.cpp`). Case 1995 reports larger, model-specific
  fitted factors, especially for His and the Trp five-membered ring, and **separates JB from HM
  calibration**. Treat the in-code constants as **hypotheses**, not best-practice values — the calibration
  learns the real coefficient.
- **The five-membered-ring lobe offsets.** Perkins 1977 / Case 1995 put the Johnson–Bovey loop height at
  ~0.64 Å; the code uses 0.50–0.52 Å for His/Trp-5. Without a stronger cite, these are hypotheses too.
- **The `2e` near the ring.** Legitimate BS↔HM divergence occurs directly above/below the π cloud, near
  vertices, in-plane, and in hetero/fused rings, where the point-loop and surface models genuinely differ
  (Perkins 1977 preferred JB above/below aromatic rings; Agarwal 1977 showed local anisotropy can be ~35%
  of the observed shift). The `2e` is defensible but less certain than the `0e` there.
- **The HM-formula equivalence.** Until the `T0`-vs-bond-sum benchmark is run, "Haigh–Mallion" is a claim
  about family, not identity (see §2).

---

## 5. The danger-zone picture

Two pictures, and the second is the one that has bitten before:

- **The near-ring pathology.** Right above the π cloud, at the vertices, and in the ring plane, the
  point-loop and surface integrals diverge and the far-field intuition breaks. This is *expected*
  physics, not a bug — the move is to **stratify** the agreement diagnostics by `(ρ, z, θ, ring_type)`
  and *report* where BS and HM part ways, rather than average over it.
- **The sign-and-units gauntlet across consumers.** Ring-current has **historical sign issues**, and the
  live risk is not in the kernel (the producer's `G = −n·B·PPM_FACTOR`, `I < 0` diamagnetic convention is
  documented and tested, and HM now matches BS) — it is in the **four places one quantity passes
  through**: the static catalog (`ppm_T_per_nA`), the BS H5 time-series (currently mislabelled `ppm`), the
  h5-reader viz (which multiplies `LiteratureIntensity()`), and the learn-side fit (which learns the
  coefficient). A unit-current kernel, a literature-scaled prediction, and a fitted shift are **three
  different numbers**; conflating them is how a sign or a factor goes silently wrong. The discipline:
  every emit carries its true unit and scale state, and **sign-convention verification is the first
  surgical step**, not an afterthought.

---

## 6. Fields, irreps, and the model

**The irreps:** full even `0e ⊕ 1e ⊕ 2e` from *both* calculators (the `2e` is the Boyd–Skrynnikov tensor
lift). Never `1o`. If a downstream observable wants symmetric shielding only, expose a secondary `0e+2e`
view — but do not replace the full emit.

**The dragon, for ring specifically.** Ring's `0e` is the **classic, genuinely-measured ring-current
shift** — the one place among the Three where the isotropic part is itself a real solution observable
(this is where ring earns its seat in Part 3). The `2e` is the DFT-angular evidence (Parts 1–2). T2 stays
sacred as that evidence, never sold as the solution-shift predictor; the solution-relevant ring signal is
the `0e`.

**The two-path, as a deliverable.** Today BS↔HM agreement exists only as a *test*
(`tests/test_batch_biot_savart_haigh_mallion.cpp` — per-ring + distance-binned convergence, the `T2`
cosine, explicitly "report, don't assert"). It is **not emitted**, and it currently compares **absolute**
values. If the thesis is to *report* the two-path validation, that comparison must become an emit carrying
**signed** residuals and explicit **sign agreement** (alongside scaled `T0` R², full-9 cosine, `T2`
cosine, norm ratios, binned by geometry) — and is most informative **stratified near the ring**, where the
agreement is non-trivial. A small, named to-build.

**The surgery (confirm-and-refine — sign first, then labels and channels):**
1. **Sign-convention verification — the first step.** Re-derive `σ = I·G` against Boyd–Skrynnikov and the
   PHE worked example before anything moves; trace the sign end-to-end through every consumer that applies
   intensity. The producer is consistent today; the job is to keep it so through the relabelling.
2. **Parity metadata.** BS/HM (and ring-chi) H5 time-series headers say `0e+1o+2e` — **wrong; should be
   `0e+1e+2e`** (the antisymmetric part is `1e`, even). Fix the writers.
3. **Naming / units.** `bs_shielding(.npy)` and the `"ppm"` H5 unit attr read as *final shielding* when
   the value is a **unit-current kernel** — rename/relabel to `ppm_T_per_nA` (and the `nA` → `nA·T⁻¹`
   source-label). HM emits geometric `Å⁻¹`, not ppm — label it so.
4. **Switch on the dropped channels, separably.** Per-type `T1` (BS and HM emit per-type `T0`+`T2`, drop
   `T1`); per-ring `B_field` / `B_cylindrical`; HM `hm_B_field`. Emit per-type **full 9**, separably;
   re-extract is cheap.
5. **The HM `T0` benchmark** (§2) — settles the name.
6. **The two-path emit** (above) — if we want to report the validation.
7. **Ring-axis contract — confirmed aromatic-only** (vet, 2026-06-07). `Protein::RingCount()` is
   aromatic-only (`Protein.cpp:287`) and the sidecar manifest's `ring_contributions` references the
   aromatic-ring axis (`TopologySidecar.cpp:573`); no fixture check needed. (The saturated-Pro index 8 is
   a ring-*type* enum value, not an extra emitted ring row.)

**Rule 1, restated:** unit kernel emitted; the literature intensity is the **direction-4 physical
hypothesis** to check the model's discovery against; the calibration learns the coefficient. Calibrated
variants should be named by set (raw unit-current / `Case95-JB` / `Case95-HM`), never silently mixed.

**Cross-kernel:** ring (magnetic) and charge/EFG (electric) **do not** overlap at the mechanism level —
they coexist; their `2e` parts merely *correlate* as geometric descriptors near aromatics, which the law
study **reports** (collinearity) rather than zeroes. The real magnetic overlap is **ring ↔ McConnell**,
handled by McConnell zeroing its aromatic category (not here).

**What this spec deliberately does NOT do:** pick final intensity constants (the calibration learns them);
assert HM-formula identity before the benchmark; claim the tensor predicts the solution shift (the dragon
— the `0e` is the measured ring-current shift, the `2e` is DFT-angular evidence); or treat the near-ring
divergence as a bug (it's stratified and reported).
