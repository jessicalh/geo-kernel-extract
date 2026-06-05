# Codex mission — restore + declutter + smoke-test the `learn/` run-making tool

Working root: `/shared/2026Thesis/nmr-shielding` (branch `h5-reader-pysr-spike`). You own the grind; the
lead vets + owns ALL git. Honest stakes: the frozen 2026-05-13/14 720 NPYs have drifted too far from the
current schema to adapt cheaply (no `.LGS`, 9-col EFG, AIMNet2 naming, FF14SB gaps). The lead's call —
**re-run the 720 fresh** instead of building an old-schema adapter, so the rediscover static-ingest eats
current-schema producer output the way it was built to. Everything's ready to fire that run **except the
run-driver itself**: it lives on `/mnt/expansion`, not the live tree. Your job: get it home, decrufted, and
**proven on 1–2 proteins** so the lead can fire a clean overnight 720 with confidence. You do NOT fire the
full run — that's her green light.

## What I already verified (don't re-derive — build on it)
- **Inputs present, at scale.** `calibration/{ID}/` holds the WT/ALA pairs: `{ID}_WT.xyz/.prmtop/_WT_nmr.out`
  + ALA equivalents (sampled `calibration/A0A062V9G2/`, dated Mar 19). ~723–724 protein dirs live under
  `calibration/` alongside `features/`.
- **Producer is current-schema.** The EFG schema rev `b580ebe` (2026-05-18, "all 5 calculators drop
  structurally-zero T0+T1") is an ancestor of `h5-reader-pysr-spike` → `build/nmr_extract` at HEAD emits
  current 5-col EFG + the current field set the rediscover consumer expects.
- **Disk:** 232 G free on `/` (88% used). A full 720 output tree is ~10–30 G — fine, not near exhaustion.
  Confirm the actual per-run size and don't silently approach the wall.
- **The gap:** `learn/extract.py` + `learn/docs/extraction_procedure.md` are **absent** from the live tree.
  They exist at `/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/`. That copy is ~3 weeks stale.
- **MOPAC** on/off is the `--mopac` CLI flag (the toml only carries mopac *parameters*). The 2026-05-13/14
  run wrote MOPAC arrays for all 720, so the tool fired it then — confirm it still does and make it certain.

## The tool (what it does)
`learn/extract.py` is the batch run-driver. Per protein it constructs:
```
build/nmr_extract --mutant \
  --wt  calibration/{ID}/{ID}_WT \
  --ala calibration/{ID}/{ID}_ALA \
  --output calibration/features/{run}/{ID} \
  --config data/calculator_params.toml
```
over the 723 WT/ALA pairs, with `--stage1-audit` (requires WT+ALA ORCA `*_nmr.out`, validates Stage-1 +
topology-sidecar outputs after each protein). Invocation form: `python3 learn/extract.py --run <name>
--resume --stage1-audit` from the repo root. It writes run logs (`extract_log.jsonl`,
`extract_background.log`, `extract_background.pid`) — clarify exactly how it backgrounds/daemonizes.

## Tasks
1. **Read first.** `/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/CLAUDE.md` (the learn/
   subproject rules) + `extract.py` + `docs/extraction_procedure.md`.
2. **Restore** the run-driver + its *direct* dependencies (helper modules, config, the procedure doc,
   `learn/CLAUDE.md`) to the live `/shared/2026Thesis/nmr-shielding/learn/`. `learn/` is gitignored — restoring
   files is NOT a git op; do not commit. **Leave the analysis workups** (`stage1-topology-workup-*`, derived
   tables, compendia) on `/mnt/expansion` — scope to the run-MAKING tool, not the analysis.
3. **Declutter, don't rewrite.** Fix stale absolute paths, dead imports, broken refs so it runs from the live
   repo root. No feature additions, no re-architecture, no new Python model. Report exactly what you changed
   and why.
4. **Build the producer.** `cmake --build build --target nmr_extract -j$(nproc)` at HEAD. **Do NOT modify
   `src/`** — the extractor is SACRED; build it, never touch it. If it doesn't build, STOP and report (don't
   "fix" src/).
5. **Confirm MOPAC ON.** Inspect the args `extract.py` actually builds; verify `--mopac` is passed for the
   mutant jobs (or that `--mutant` defaults it on). The 720s MUST get MOPAC (lead policy). Make it explicit
   and certain; if ambiguous, surface it — don't guess.
6. **Smoke-test on 1–2 proteins ONLY** (e.g. `--run SmokeTest_20260604` over a 1–2 ID subset, or the resume
   mechanism on a fresh tiny run). Confirm end-to-end: job succeeds, MOPAC arrays written, and the output is
   **current schema** — 5-col EFG, the current field set. Cross-check one protein's emitted field
   names/shapes against what the worktree static-ingest consumer expects (`h5-reader/src/rediscover/` +
   `RunData::staticPoseExpectedFields()` / `QtFieldCatalog.gen.h`). Keep it lean — 1–2 proteins, no large
   data.

## Hard boundaries
- **Do NOT fire the full 723 run.** Smoke-test 1–2 only. The full overnight fire is the lead's green light.
- **Do NOT modify `src/`** (SACRED extractor) — build only.
- **Do NOT touch git** — the lead owns ALL git. `learn/` is gitignored so restoring files is fine, but no
  `git add/commit/branch/switch/merge/rebase/reset`. Branch stays `h5-reader-pysr-spike`.
- **No sprawl** — clean/restore the tool, don't grow it; no second Python model, no terabyte anything.
- **Locate before absent** — if an expected input/file is missing, log the exact path and STOP; no glob,
  no try-and-fail.

## Report (conclude with this)
1. **Restored:** what files you brought to live `learn/`, and what you deliberately left on `/mnt/expansion`.
2. **Cleaned:** every change you made to make it run, with the reason.
3. **MOPAC:** confirmed-on evidence (the exact arg the tool passes).
4. **Smoke test:** the run command, the result, proof the output is current-schema (field names + a couple of
   shapes, EFG width), MOPAC arrays present.
5. **The fire command:** the EXACT `python3 learn/extract.py …` line for the full 720, how it backgrounds,
   where it writes (`calibration/features/{run}/`), estimated runtime (~16 h) + output size/disk.
6. **Risks / open questions** for the lead before she fires.
