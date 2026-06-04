# Session checkpoint — 2026-06-04 (state stash; context-decay insurance)

A single comprehensive checkpoint of this session's arc + state so a fresh session (or post-context-decay)
resumes losslessly. It POINTS at the detail docs and SYNTHESISES the current state + open threads. Branch
`h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR; **the LEAD owns ALL git.** Mechanism: codex grinds,
opus reviews, lead vets + owns git; the maths discussion is lead+Claude.

## The session arc
Stashed the analysis/spec work as **specs** (codex briefs → opus reviews → sharpen), built the **interp
capability demo**, stashed the **thesis overview**, built the **720 ingest in a worktree**, then pivoted to
the **reader (LGS + UI)** and evaluated its state. We are now on docs/state before the interactive UI sessions.

## Track state

### Specs — DONE, stashed (`h5-reader/src/rediscover/`)
Four specs, codex-drafted + opus-reviewed + sharpened: `SPEC_720_STATIC_INGEST`, `SPEC_MATHS_FIX_AND_REVIEW`,
`SPEC_LEARNING_MODEL`, `SPEC_INTERP_1P9J_SCOPING` (+ `REVIEW_opus_SPEC_*` and the codex adversarial 720 review).
**Parked maths verdict re-runs** (the corrected LOAO re-fit + the e3nn clean-protocol re-run) remain parked —
see `MATHS_AUDIT_CHECKLIST` + `SPEC_MATHS_FIX_AND_REVIEW`. The maths *discussion* (VALIDITY-vs-ATTRIBUTION) is
lead+Claude, deferred.

### 720 ingest — HALF DONE (built, not run, not merged)
- **BUILT** in the worktree `/shared/2026Thesis/nmr-720` (branch `wt-720-build`). Full static-pose ingest:
  `StaticPoseConformation`, strict exact-path NPY loader, `--static-cohort … --memory-strategy …` CLI, the
  **all-atom** emitter with **BMRB/IUPAC + element/role/stratum labels on every row**, the six sidecar
  families + **T2-orientation handling**, the manifest, both memory backends. Compiles
  (`build/linux-rwdi/h5reader_extract`). **1P9J oracle-parity PASSES.**
- **PENDING:** the actual 720 run (curated data not yet available — on the spinner). The merge
  (`wt-720-build` → `h5-reader-pysr-spike`) is the lead's, after review.
- **MERGE-REVIEW FLAGS:** (1) a `nanoflann.hpp` **symlink** codex used as a build-dep workaround — violates
  no-symlinks; needs a real copy / CMake include fix. (2) two honest open-spec calls codex made: the cohort
  manifest key shape, and thin-support threshold = 10 (recorded in the manifest).
- **Dual objective:** the between-axis **stats** (Stage 2) AND the **Model-1 training corpus** (Stage 3) —
  1P9J's substrate + the 720 traversal = the corpus. Docs: `SPEC_720_STATIC_INGEST`, `CODEX_BRIEF_720_BUILD`,
  `REVIEW_codex_adversarial_SPEC_720`.

### Interpolator — DONE demo; live-prediction TABLED (wanted back, not today)
- **BUILT:** source-e3nn v0, equivariant shielding-tensor prediction on 1P9J. Capability graph (held-out
  **T2 R² ≈ 0.48**, |T2| r ≈ 0.82). In `/tmp/rediscover-runs/2026-06-04-interp-1p9j-e3nn/` + standalone
  `analysis/interp_1p9j_graph.py`.
- **TABLED (lead wants it back, NOT today):** the **live-reader-prediction** (the chewer killer-app — the C++
  reader asks a GPU model live). Design notes: `NOTES_INTERP_LIVE_PREDICTION` (server vs libtorch vs
  in-reader; "a server is its own complexity"; the C++↔GPU boundary + format contract are the open design).
- **OPEN:** the `analysis/equiv_t2_e3nn.py` build-diff (+441/−49 from the interp build churned the exemplar;
  the re-cut graph went to a standalone script so NO further churn — but verify the build's exemplar changes
  didn't disturb the clean-protocol exemplar the parked maths e3nn re-runs depend on; consider a derived script).

### Thesis overview — STASHED for the lead (`doc/thesis-overview/`)
`PRECIS_ADVISOR` (the desk doc: Stage 1/2/3) + the learning-model spec + the Stage-1 understanding + the
interp notes + the capability graph + a README. Body complete; **lead finalises** (figure-slot, her voice,
LaTeX). **Stage 3 = two-model e3nn pipeline:** Model 1 (e3nn shielding emulator, DFT-trained on 1P9J + 720)
→ Model 2 (shift predictor, BMRB-trained on the no-DFT ~656 fleet, eats Model 1's firehose, ablate Model 1
in/out). **e3nn, not MACE** (lit scan: e3nn is the direct-precedent SOTA for shielding tensors; MACE-on-shielding
is unprecedented — a defensible novelty if ever chosen). The **0.818 (110-protein) / 0.718 (720-scale)**
correction is propagated to CLAUDE.md + memory. *[Worth a memory update to `project_stage3_equivariant_gnn`:
the two-model split + e3nn-not-MACE — held pending the lead's MEMORY.md curation.]*

### Reader / UI — the NEXT phase (several interactive sessions)
- The reader is a **real Qt/VTK app**, not a prototype: typed `.LGS` load, the model spine, the VTK scene,
  playback, picking, ordered selection, a modern camera composer, the transform layer, overlays + docks + REST.
  Qt discipline mostly honoured. Evals: `notes/STABILISATION_FEATURE_EVAL_2026-06-04.md` (stabilisation is real
  — `TransformedConformation` + `CameraComposer`; plane lock tested) + `notes/UI_STATE_OVERVIEW_2026-06-04.md`
  (full doc/code state, the infelicities, the install gap).
- **MAIN QUESTION (in flight):** did the long CLI/rediscover work (shared model/io) break the reader? Codex
  build-verify `b0gpfdrxr` running; the lead runs the runtime (she knows the app well).
- **Three goals (interactive sessions):**
  1. **Open 1P9J + work + not-suck** — the infelicities (the overview's list ≈ the lead's): debug-flavoured
     toolbar; "Instrument" user-facing label; no selection Clear; per-frame overlay log noise; the default
     `dssp_chi` dashboard signal; **`QtFrame::eeqCharge` returns placeholder 0** (could mislead — don't ship
     that if EEQ shows). In concert.
  2. **The RADIUS view (headline, wanted badly):** pick atom + radius → show only that atom + the residues
     within R. NOT built; clear narrow path (a `displayMolecule_` + `AtomFilter::WithinRadiusOf` + source↔display
     index map + isolate-toggle/radius-spinbox/clear UI + optional REST). Builds on the real stabilisation
     infra. Risk: picker/selection over the filtered molecule — mitigate by picking over source positions + the map.
  3. **Linux install** — mostly easy (lead is prepping a test machine): `install(TARGETS)` + `$ORIGIN` rpath +
     **AppImage** bundling Qt/VTK/HDF5; none exists yet. (The qt6-cpp skill is Windows-targeted; codex has a
     Linux-tailored version. Win/Mac use Pro Qt licences; Linux is open-source Qt.)
- **In flight:** codex doc-decruft (`bcfazv9yw` → `notes/DOC_DECRUFT_TRUE_LIST`) + codex reader-build-verify (`b0gpfdrxr`).

## Open threads / next steps
**ORDER (lead, 2026-06-04): doc-fix → integrate the 720 (merge `wt-720-build`) → the interactive UI sessions.**
Integrating the 720 first lands the shared model/io changes into the main tree before UI work, so the merge
can't disrupt it later, and the reader is verified against the final state.

- **720:** lead merges `wt-720-build` after the symlink fix + review; the run waits for the curated data.
- **Interp:** the equiv-diff cleanup; the live-prediction (tabled, wanted back) — a focused chewer session.
- **Maths:** the parked verdict re-runs, on resume.
- **Overview:** lead finalises the précis.
- **Reader:** act on the doc-decruft list; confirm did-we-break-it; then the interactive UI sessions
  (infelicities + radius view + AppImage).

## Disciplines (load-bearing)
Branch never merge/switch; **lead owns ALL git**; codex grinds + opus reviews + lead vets; plan-vet risky
loops; **no symlinks**; model-is-spine / no Python second model / lean disk <15 G / extractor SACRED; the
maths discussion is lead+Claude; e3nn not hand-rolled; T2 sacred; for reader code the **qt6-cpp skill** +
CENSUS/ACONNECT/ASSERT_THREAD + UDP-log-first + one-persistent-QTimer.
