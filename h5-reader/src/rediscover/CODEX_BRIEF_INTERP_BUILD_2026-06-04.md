# Codex brief — BUILD the 1P9J interpolation tool + the advisor graph (the lead GREEN-LIT this)

Status: **GREEN-LIT by the lead (decision #3, 2026-06-04). Build it.** Branch `h5-reader-pysr-spike` —
NEVER merge/switch/rebase/PR; the lead owns ALL git (you build + gate, you do NOT commit).

You own the grind; the lead vets + judges. This is the advisor payoff: a **graph on the desk that shows we
are HEADING toward the Stage-3 maths** — even a simple, honest one is the win. It IS (mostly) the equivariant
model, extending the existing e3nn exemplar — **never a hand-rolled second model**. Realistic labelling:
within-axis, geometry-sampler, **direction not destination**.

## READ FIRST
- `SPEC_INTERP_1P9J_SCOPING_2026-06-04.md` — the build plan (it doubles as this build's spec).
- `REVIEW_opus_SPEC_INTERP_2026-06-04.md` — **APPLY its 3 sharpenings** (below).
- The e3nn exemplar to EXTEND (do NOT re-author): `analysis/equiv_t2_e3nn.py`, `analysis/e3nn_protocol.py`
  (the clean blocked/purged protocol — USE it), `analysis/change_of_basis.py` (frozen get_C), `analysis/FINDINGS.md`.
- Top of `STATE.md` (STAGE 2 block).

## THE THREE OPUS SHARPENINGS — bake them in
1. **The per-source CSVs are GONE from disk → this is a gated PRODUCER RE-EMISSION, not "packaging."**
   Regenerate the 1P9J per-source CSVs (`ring_current_sources.csv`, `*_target_local_T2.npy`) by re-running
   the rediscover per-source emit — **extractor UNMODIFIED (SACRED)**, **bounded to the dominance-clean
   ring/backbone exemplar subsets ONLY** (NOT all 558,360 atom-frame rows — that is the forbidden terabyte
   path). Byte-parity/additive sub-gate where an oracle exists; print the row bound + a `rows×cols×8B` size
   estimate BEFORE writing; lean disk.
2. **Anchor the gate to the CLEAN number THIS build produces — NOT the leaky-suspect 0.466/0.757.** That
   headline is pre-protocol-fix; tag any pre-fix number "leaky-suspect" everywhere. The recovery the graph
   reports is whatever the CLEAN protocol (`e3nn_protocol.py`: blocked/purged split, train-only centering)
   gives. The build produces its own honest number; the ≥0.25 "worth showing" gate is anchored to that.
3. **Panel C anti-flatter is MANDATORY.** Baseline-restored σ_iso flatters (the per-atom chemical baseline
   swamps the modulation). Put the **modulation-R² on-panel**; demote restored-σ_iso to a labelled inset.

## BUILD STEPS
1. **Re-emit** the bounded per-source 1P9J exemplars (sharpening #1).
2. **Train the source-e3nn v0** — extend `equiv_t2_e3nn.py`, clean protocol, held-out (blocked/purged),
   **5-component T2** (frozen `get_C`, NEVER scalar-collapse). Report held-out T2 recovery + |T2| corr.
3. **The 2×2 advisor graph** (the deliverable): (A) held-out predicted-vs-DFT T2 scatter; (B) |T2| recovery;
   (C) σ_iso/T0 companion WITH the anti-flatter fix (#3); (D) the honest **"direction not destination"**
   caption. Use the project chart tooling (R / matplotlib); emit PNG + PDF.

## GATES (commit only green; the lead gates)
- Byte-parity/additive on the re-emit; held-out blocked/purged; **disk: `df` before any write, abort if
  free < 20 G, total rediscover < 15 G, deletes by EXPLICIT FULL PATH only, never `/shared`**; SETI (report
  the honest number, no overclaim). What proves it works = the graph + the held-out recovery the clean
  protocol yields.
- torch/e3nn env: use the existing `analysis/venv`; set `LD_LIBRARY_PATH` for torch cu130
  (the nvidia cu13 lib) before running — see the project's CUDA note.

## DISCIPLINE (load-bearing)
- **Extractor SACRED** (re-run unmodified). **Model is e3nn — extend the exemplar, NEVER hand-roll a second
  model.** **T2 sacred** (never scalar-collapse). **NO protein/spatial model in Python — not even a
  secondary aggregate; NO terabyte; Python consumes the emitted substrate only; spatial work stays in C++.**
  Frozen `get_C`. Never open `trajectory.h5` in Python. **NO git** (leave everything uncommitted — the lead
  owns git). Branch never merge/switch.
- Realistic labelling mandatory: within-axis, geometry-sampler, correlate-not-match, direction-not-destination.

## OUTPUT
The graph (PNG+PDF) + a short `BUILD_INTERP_RESULT_2026-06-04.md`: the honest held-out number, the caveats,
where the graph is, and what each panel shows. End by printing a 4-6 line summary (the number, GO-worthy or
not, the graph path, anything that didn't ground).
