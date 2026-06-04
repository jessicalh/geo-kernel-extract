# Thesis overview — all three stages (advisor précis + prep)

Stashed 2026-06-04 to carry further. The **précis is the desk doc**; the rest is prep behind it.

## Files
- **`PRECIS_ADVISOR_2026-06-04.md`** — the overview/précis (the desk-facing doc). Body complete across all
  three stages. **What's left:** slot the figure (placeholder marked in Stage 3), a pass in your voice
  (it's a draft-to-best scaffold), and a LaTeX/PDF render for the desk.
- `SPEC_LEARNING_MODEL_2026-06-04.md` — the deeper internal learning-model spec (technical companion the
  précis summarizes). *Note: its "v0→v5 staged" framing predates the un-phasing — the précis Stage-3 is the
  current, un-phased two-model version; the spec should be de-phased to match when you revise.*
- `LEARN_STAGE1_UNDERSTANDING_2026-06-04.md` — the recovered Stage-1 numbers (R² 0.818 on 110 / 0.718 at
  720-scale; 126 kernels; the nitrogen-pooling-artifact; the 0.81-within/0.35-across ceiling).
- `NOTES_INTERP_LIVE_PREDICTION_2026-06-04.md` — the interpolator + the **tabled** live-reader-prediction
  design (server vs libtorch vs the in-reader trap; "a server is its own complexity"; the killer-app idea).
- `interp_1p9j_advisor_graph.{png,pdf}` — the **capability graph** (the Stage-3 figure): equivariant
  shielding-tensor prediction on 1P9J, held-out T2 R² ≈ 0.48.

## The shape it lands on (Stage 3)
Two GPU models, **un-phased**, distill → transfer → test:
- **Model 1 — shielding emulator:** per-atom **e3nn** equivariant GNN, per-frame/local, predicts the full
  shielding tensor (T0+T2; T1 only if wanted — shielding is non-symmetric). **Trained only on DFT** (1P9J's
  many frames + the 720) so the physics stays clean; predicts shielding where there's no DFT. *Make this one
  first — it's the prototype.* e3nn (not MACE — the direct-precedent SOTA for shielding tensors); the readout
  irreps are the crux; the **DFT corpus's size/diversity (the 720) is the binding constraint**, not the architecture.
- **Model 2 — shift predictor:** runs on the **no-DFT fleet (~656)**, eats the per-frame **firehose** (Model 1's
  shieldings + features), **trained to predict the measured shift** vs BMRB; ablate Model 1 in/out to test
  whether the distilled DFT earns its place. The clean DFT trains the physics; the noisy shift trains the
  observable — the experiment never corrupts the physics model.
- **1P9J is the calibration case** (has both DFT and BMRB); then the 656.

## Grounding (lit scan, 2026-06-04)
Full-tensor shielding SOTA is equivariant GNNs — e3nn (²⁹Si tensors, ~53% error cut vs invariant) and PET
(ShiftML3). MACE has reached NMR only as a scalar descriptor; MACE-class nets emit per-atom rank-2 tensors
for polarizability/Born charges, so MACE-on-shielding would be a defensible novelty — but **e3nn is the
direct precedent**, which is why we use it. (Sources in the scan; pull arXiv:2412.15063 + the MACE-MDP PDF
for verbatim numbers.)
