# Directory Set: Biot-Savart NMR Shielding Tensor Prediction

Four directories constitute the complete project state as of 2026-03-31:

1. **working-reference/** — Frozen snapshot of a known-good state. Do not
   modify. Reference for reproducing published results.

2. **biot-savart/** — Current state with all contradictions exposed.
   v1/v2 divergence documented. Bugs catalogued. Goals NOT updated.
   This is the honest broken state we are rewriting from.

3. **rewrite-spec/** — Design documents for the clean-room rewrite.
   Multiple perspectives (Opus, Gemini, xHigh). Adversarial reviews.
   Architecture decisions. Feature specifications. No code here.

4. **mutant-build/** — New mutant DFT calculations.
   Charge-flip (ASP→ASN), ring modification (TYR→PHE),
   salt-bridge breakers (LYS→MET). Structure prep and ORCA inputs.
