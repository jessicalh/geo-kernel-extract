# Outline — Introduction and Poster

Draft scaffold (yours to cut against). It combines the three core chapters with Stage 1,
the ridge calibration model, and the Stage-3 ablations.

## The three core chapters (the conceptual spine)

- **Physics Architecture** — the unified view: the calculators as classical reductions of
  one quantum object (the Ramsey current response); the recurring kernel
  `D_ab = ∂_a∂_b(1/ρ)`.
- **Categorical Grounding** — the categorical view: by atom type and geometry; why each
  measure exists; the determinant literature; which calculator reaches for which nucleus.
- **Method and Proxy** — the intellectual shape: what each method proxies; the
  literal→metaphorical spectrum; the geometric kernel as the question posed to nature,
  calibration as the answer.

---

## Introduction — outline

1. **The question, and why it resists a single answer.**
   - NMR chemical shielding is an exquisite local-structure readout; predicting it from
     geometry is the goal.
   - It is a *categorical, multi-determinant* problem — not one derivation from quantum
     chemistry. Different nuclei have different dominant determinants.
   - The state of the art is the geometry→shift machine (SHIFTX, SPARTA, ProCS): accurate,
     largely empirical, black-box.

2. **The stance: a mechanistic, interpretable decomposition.**
   - Not a black box — a suite of physically-motivated geometric kernels, each a *named
     classical reduction* of a real shielding effect (ring current, bond/group anisotropy,
     electric field / field-gradient, hydrogen bond, dispersion, multipole, hydration).
   - The contribution is not R²; it is *what carries the signal, and why*.
   - The three chapters establish the foundation:
     - *Architecture* — what they are: shadows of one quantum response; the recurring kernel.
     - *Categorical* — why they exist and where they sit: the per-nucleus, per-geometry
       determinants the field actually reads.
     - *Method & Proxy* — what kind of move each is: what its method proxies,
       literal→metaphorical; the kernel as the question, calibration as the listening.

3. **The empirical foundation — do the kernels carry signal? (Stage 1.)**
   - Per-element, per-atom-type ridge on 720 proteins / 446K atoms, 55 kernels,
     R² = 0.818 (settled).
   - The "nitrogen is hard" resolution: an element-pooling artifact — backbone N is hard
     (R² ≈ 0.39), sidechain N is among the best (R² ≈ 0.89). The per-atom-type complexity
     *is* the story (H ≈ 20 effective dims, C ≈ 6, N ≈ 3, O ≈ 12).
   - The learning model is ridge — a fittable calibration, not an MLP (tested and rejected).
     Calibration is what turns geometric kernels into shielding.

4. **The dynamical dimension. (Stage 2.)**
   - The observable is an ensemble/time average; MD (GROMACS, ff14SB) is the sampler; the
     calculators read each frame.
   - The DFT reference (r²SCAN on the MD frames) factors the MD-geometry model *out* of the
     kernel-vs-DFT calibration; the "MD is a model" worry re-enters only versus experiment.
   - 1P9J as the over-examined single system; the network of physical laws in one trajectory.

5. **The evaluation — what matters. (Stage 3, ablations.)**
   - Ablations: remove kernels / components to test which carry *independent* signal versus
     redundant geometry.
   - The model evaluation: where the signal lives, per atom type. *(Results pending — panel
     placeholder until Stage 3 runs.)*

6. **The contribution.**
   - A mechanistic, interpretable, honestly-proxied decomposition of protein shielding: the
     unity (architecture), the determinants (categorical), the intellectual shape of each
     tool (method & proxy) — grounded by the Stage-1 calibration and tested by the Stage-3
     ablations.
   - Not a predictor chasing R²; a way of asking the universe what carries the shielding
     signal, and being exact about how each tool poses its question.

---

## Poster — outline (panels)

- **Title + one-line thesis.** "A mechanistic, interpretable decomposition of protein NMR
  shielding into physically-grounded geometric kernels — what each is, why it exists, and
  what it proxies."
- **Panel 1 — The hidden object.** The shielding tensor / Ramsey current response; the
  calculators as shadows of one object. *(Architecture figure: one kernel → three families →
  12 tools.)*
- **Panel 2 — The categorical reality.** A nucleus → determinants → calculators table
  (Cα/Cβ ← φ/ψ + secondary structure; HN ← H-bond + ring + field; ¹⁵N multi-determinant;
  aromatic-near ← ring current). The field reads by atom and geometry.
- **Panel 3 — Method & proxy.** The literal→metaphorical spectrum with the 12 tools placed
  (BiotSavart literal … HBond borrowed). The kernel is the question; calibration is the
  answer.
- **Panel 4 — It carries signal (Stage 1).** R² = 0.818; the per-atom-type bar (backbone N
  hard, sidechain N best); ridge, not MLP.
- **Panel 5 — What matters (Stage 3 ablations).** Which kernels carry independent signal, by
  atom type. *(Placeholder — results pending.)*
- **Take-home.** Geometry asks the question; calibration listens. Mechanistic and honest, not
  a black box.

---

## Notes / open decisions

- The Stage-3 ablation panels (intro §5, poster Panel 5) are structure-only until the
  ablations run; they inherit the Stage-2 trajectory results.
- The triad of chapters is the methods/conceptual core; Stage 1 supplies the "it works"
  evidence and the model; Stage 3 supplies "what matters."
- Voice throughout: the same Sagan-with-equations / mechanistic-not-predictive register as
  the chapters.
