# Codex brief — BUILD the 720 static-pose ingest (worktree; lead green-lit)

> **Historical session brief — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Status: **lead GREEN-LIT.** Runs in the **`wt-720-build` worktree** (its own `build/` dir). Branch never
merge/switch. **LEAD owns ALL git** — you build + gate in the worktree, you do NOT commit; the lead merges
`wt-720-build` into `h5-reader-pysr-spike` when green. Build the whole ingest — **no phasing.**

## THE DUAL OBJECTIVE (state it in the manifest; it drives the design)
The 720 ingest delivers TWO things from ONE substrate:
1. **Stage-2 between-axis statistics** — the cross-protein stats instrument (the original purpose).
2. **The Stage-3 Model-1 training corpus** — the SAME per-(protein, atom) substrate (full feature pile +
   DFT shielding target + geometry) IS the training data for **Model 1**, the per-frame equivariant
   *shielding* predictor. **1P9J's substrate (already emitted) + the 720 traversal = the training corpus.**
   Model 1 trains on all the DFT shielding we have (1P9J's many frames + the 720 static poses).
The adversarial fold already requires the full feature families + the T2-orientation handling — exactly what
the equivariant Model 1 consumes. So the stats requirements and the corpus requirements ALIGN: one
substrate, both jobs, no conflict. Record both objectives in `static_manifest.json`.

## NON-NEGOTIABLES (lead gate-check — hard, not preferences)
1. **ALL ATOMS — no toy, no backbone-only, no stratum subset.** Every atom of every protein is a row
   (per-(protein, atom), all atoms) — the full corpus, not a demo slice. The interp v0's ring/backbone-only
   scoping does NOT apply here.
2. **Both labelings on every atom row, across the whole corpus — NEVER sliced at the end.** Each atom row
   carries BOTH (a) its **BMRB/IUPAC** identity (for the experimental-shift / Model-2 alignment) AND (b) its
   **atom category** (element / atom-role / stratum). Downstream training must be able to stratify by EITHER
   or BOTH **throughout** — per-category and per-BMRB/IUPAC across the whole thing, never pooled-then-sliced
   post-hoc (the Stage-1 per-type lesson; Build-2's H anti-prediction was the bill for pool-then-slice).
3. **No second protein model in Python — not even an aggregate** (LAW 2; already in the spec). Python only
   fits the lean emitted substrate.
4. **The corpus feeds e3nn Model 1** — the per-frame equivariant *shielding-tensor* predictor (**e3nn**, the
   direct-precedent SOTA for shielding tensors; NOT MACE). So the substrate must be in the form e3nn consumes
   — positions + per-atom features + the T2-orientation handling + the DFT target (already required by the fold).

If any of these cannot be satisfied from the spec, STOP and flag it — do NOT improvise a toy.

## READ FIRST (by absolute path from this worktree — the docs live in the main tree)
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/SPEC_720_STATIC_INGEST_2026-06-04.md` — the
  build plan (adversarial-folded; the requirements are in it).
- `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md` + `REVIEW_opus_SPEC_720_2026-06-04.md` (same dir) — the vetted requirements.
- The grounding code the spec cites: `RunData`, `../model/SingleConformation`, `../io/QtProteinLoader::LoadPose`,
  `../io/QtTopologySidecar`, `../io/QtNpyReader`, `../io/DftShieldingLoader`, `ResidentIndexes`, `Verbs`, `PerAtomSubstrate`.

## BUILD (whole, per the spec)
- `StaticPoseConformation` + generalize the trajectory-only `RunLoader` (drop the `Trajectory`-kind gate +
  the `h5()` requirement; `static_frame_index = 0`).
- The **strict expected-field loader** — resolve each required NPY by exact path, log-and-stop; NEVER the
  `FrameNpyLoader` glob/enumerate.
- The `StaticCohortAccumulator` with **both** backends (all-resident + streaming); default UNSET — the
  memory strategy stays open (maths-dependent).
- The emit: the lean per-(protein, atom) substrate + the grouped C++ sidecars the fold requires (ring
  valid/self split; per-type/per-category mechanism; electrostatic source slabs; hbond/solvation; separable
  AIMNet2 embedding) + the **T2-orientation handling** (local-frame / invariant T2; NO Python tensor
  rotation) + the DFT target (raw ORCA total/dia/para). The manifests + the disk preflight (< 15 G).

## GATES (build green in the worktree)
- Compiles in the worktree's own `build/<preset>` dir.
- **1P9J oracle-parity (ring + McConnell/MOPAC-Mc, v1) PASSES** — build + test against 1P9J (we have it).
  The curated 720 data is **pending** (on the spinner / awaiting the producer MOPAC re-run), so the actual
  720 traversal is a later run; this build must be **ready and oracle-clean** for when the data lands.
- The FORECLOSURE holds (no terabyte of Python, no second protein/spatial model in Python — not even an
  aggregate; spatial work in C++). Extractor SACRED, model-is-spine, no file discovery, frozen `get_C`.

## OUTPUT
The built + 1P9J-oracle-parity-passing ingest; a build summary (what's built, the oracle result, what's
tested vs pending-on-720-data, the dual objective recorded in the manifest); and any spec gaps you hit.
Leave everything uncommitted in the worktree.
