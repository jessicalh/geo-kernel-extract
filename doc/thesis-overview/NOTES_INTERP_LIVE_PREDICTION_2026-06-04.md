# Interpolator + live-reader-prediction — notes & decision (2026-06-04)

**Status:** the interpolator is BUILT (capability graph done). The **live-reader-prediction** idea is
**TABLED for tonight** — it could be the killer app, but it deserves a focused session, not a tired one.
This doc stashes the thinking so we resume without re-deriving. Tonight's reader work is **LGS + UI**.

Related: `BUILD_INTERP_RESULT_2026-06-04.md`, `SPEC_INTERP_1P9J_SCOPING_2026-06-04.md`,
[[project_stage3_equivariant_gnn]], [[feedback_model_is_spine]].

## What's built (the interpolator)
- **Source-e3nn v0:** equivariant shielding-tensor prediction on 1P9J — trains on a within-protein frame
  subset, predicts held-out frames. Clean held-out protocol (blocked/purged, train-only centering, frozen
  `get_C`, 5-component T2 intact). Bounded re-emit (141k source rows, not the 558k all-atom).
- **Result (held-out):** T2 component R² **+0.478**, |T2| r **+0.819**, T0 modulation R² **+0.454**.
- **Advisor graph:** capability-first 2×2, "Equivariant shielding-tensor prediction — 1P9J", from the
  standalone `analysis/interp_1p9j_graph.py`. Graph + predictions in `/tmp/rediscover-runs/2026-06-04-interp-1p9j-e3nn/`.

## Why it matters — three framings (lead)
1. **A test of building an equivariant predictor of any kind** — built, runs, recovers. Capability proven.
2. **The advisor demo** — the data doing work: predicts within the protein from the protein's own data.
3. **The commercial per-protein bootstrap** — DFT a *subset* of a protein's frames + our tools → predict the
   rest of *that* protein at ~50% (the held-out T2 R² ≈ 0.48 *is* this number). Not good; better than nothing; real.

## The killer-app idea — live prediction FROM the reader
The win is NOT a precomputed lookup. It is the C++ reader holding **live state** and **asking a learned model
on a GPU at runtime**, getting back a non-deterministic prediction. That is the engineering win *and* the
first real exercise of the **chewer boundary** the whole Stage-3 vision rests on.

## The C++ ↔ GPU-model query — design space (the notes)
**Hard constraint:** the model CANNOT live in the reader. The reader must stay lean + cross-platform
(Win/Mac/Linux); a GPU torch model = CUDA + torch libs, and **Macs have no CUDA**. So inference is a
**separate process regardless** — on a Mac the reader just talks to the GPU box over the wire (same shape).

Options:
- **A) Python prediction server + reader as a thin C++ socket client.** Reader stays PURE C++ (opens a socket
  — already does, for UDP logging — and asks). Python is isolated in a tiny server (the model's native
  runtime). **Caveat (lead): a "hacked server" is its OWN ongoing complexity** — process lifecycle, protocol
  drift, error handling, versioning — invisible when optimizing for fast+simple, real later. Not free.
- **B) libtorch + TorchScript in a separate C++ inference process.** No Python, but a heavy CUDA dependency
  (plus the known broken `fs::remove_all` gotcha — [[feedback_libtorch_broken_filesystem]]), e3nn ops are
  finicky to TorchScript, and it is the **same architecture** (separate process + socket) for marginal
  benefit (deletes a ~50-line Python file).
- **C) THE TRAP — libtorch IN the reader.** Breaks cross-platform (CUDA / no Mac) and bloats the lean reader. Don't.

**Decision axis:** "no Python *in the reader*" (either server works; the reader stays clean) vs "no Python
*anywhere in the path*" (→ the libtorch C++ server). Layered on top: a server is its own complexity, not free.

**Ready if/when we resume:** the per-(atom, frame) predictions already exist (41 atoms × 500 frames,
train/held-out flagged). **Open:** the C++↔model boundary choice (server vs libtorch) + the request/response
**format contract** (current-state-in → prediction-out; JSON prototype vs compact binary).

## Decision (2026-06-04)
**TABLED for tonight.** Could be the killer app (reader asks a live GPU model — the chewer / Stage-3 seed).
Revisit in a focused chewer/reader-prediction session. Tonight: the reader story rides the **static
capability graph**; the engineering build of the live query waits.
