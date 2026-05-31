# Next-session kickoff — rediscover, through the Python side

Continuing the `rediscover` work on branch `h5-reader-pysr-spike`.

READ FIRST, in order: `h5-reader/src/rediscover/STATE.md` (current state +
the Codex must-fix list + build/run commands + the 1P9J calcset path), then
`GUIDANCE.md` and `DESIGN.md` in the same dir. Memory:
`project_rediscover_state`, `feedback_build_inmemory_export_dont_relitigate`.
Do NOT re-open the design or re-litigate the fitter choice.

Status: the extractor builds, runs on 1P9J, and emits the four CSVs; Codex
confirmed the scalar physics is correct. Build/run ONLY works in the lead
session (subagents' Bash is sandboxed from the toolchain):
`cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests`
with the sandbox override. Restore agent build-agency if you want agents to
compile.

Decision recorded: aromatic-H stratum = **C–H only** (HA/H4/H5) for now —
with the honest caveat that the stacked narrowing choices may limit
statistical depth (report effective N).

## Arc for this session

1. **Apply the value-affecting Codex fixes, re-run, re-emit.** From the
   STATE.md list, the ones that change the numbers a fitter sees:
   (a) exclude the self/bonded ring from the ring-current aggregate — emit
   `ring_index` + `is_self_or_bonded`, split `sum_dipolar_all` vs
   `sum_dipolar_producer_valid`; (b) aromatic-H frame anchor → typed CG/CD2
   per `substrate_conventions` (emit `frame_anchor_atom_index`); (c) flip the
   SVD ring normal to the canonical traversal; (d) McConnell: emit the
   bond-axis vector in the local frame + endpoint indices, and record the
   cutoff (use 10 Å to match the producer). Re-run on 1P9J. Have Codex
   critique the fixes.

2. **The Python side — is there anything here?** Run fitters on the CSVs to
   see whether signal drops out. Start cheap, per stratum, on the aggregated
   tables: ridge and scalar SR (PySR) of `σ_iso` (T0, basis-free) against the
   summed `(3cos²θ−1)/r³` + per-ring-type / per-bond-category sums — does the
   ring-current / McConnell form and a literature-plausible coefficient
   appear? Claim shape is **correlate, not match**. The fitter is OPEN: these
   scalar looks first; set up toward the equivariant / sum-pooling path for
   the T2 story — but only after the T2 Cartesian-frame caveat in STATE.md is
   resolved (T0 is safe regardless of frame).

3. **Report honestly, depth in view.** The substrate is deliberately narrowed
   (single protein, DFT subset, C–H only, cutoffs). Report effective N per
   stratum alongside any fit, and read a weak/null result as a possible
   statistical-depth limit, not necessarily absent physics. Disappointments
   are findings. No hyperbole.

## Discipline
Build/verify in the lead session (or Codex); reuse the reader; additive
edits only; GUI untouched; truthful docs. The full Codex must-fix list,
the T2-frame resolution, and the 1P9J path are in `STATE.md`.
