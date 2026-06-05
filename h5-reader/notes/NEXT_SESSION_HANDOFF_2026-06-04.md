# Next-session handoff (2026-06-04) → integrate the 720, make the combined model, do the UI

> **UPDATE 2026-06-04 (evening) — partly superseded; see `notes/720_RUN_STATE_AND_PLAN_2026-06-04.md`.**
> The version-surgery / adapt-the-frozen-NPYs premise below is **MOOT**: we re-ran the 720 **fresh**
> (current schema). Run `Stage1BMRB_20260604` is cooking (task `b41z04q8e`, ~8–12 h); first 12 validated at
> ORCA-print tolerance. The C++ ingest needs its CRITICAL fixes (see the worktree review) before it makes a
> trustworthy corpus. The reader-UI section (Track B) stands as written.

Paste to a fresh session. Branch `h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR to master. **The LEAD
owns ALL git** (you edit/build; the lead commits, except a commit she explicitly delegates). Mechanism: codex
grinds + opus reviews + lead vets; **the maths discussion is lead+Claude, not codex**; plan-vet risky loops.

## STEP 0 — ORIENT (read, in order)
1. `h5-reader/notes/FRESH_LOOK_2026-06-04.md` — a codex fresh-eyes readiness review (must-fix-first + the
   version-surgery specifics). **Read it first — it's the go/no-go.**
2. `h5-reader/notes/SESSION_CHECKPOINT_2026-06-04.md` — the full prior-session state.
3. `h5-reader/src/rediscover/LEARN_720_RECIPE_AND_DATA_2026-06-04.md` — the frozen 720 data (usable, skip the re-run) + the recipe.
4. Reader: `notes/UI_STATE_OVERVIEW_2026-06-04.md` + `notes/STABILISATION_FEATURE_EVAL_2026-06-04.md`.
5. Arc: `src/rediscover/SPEC_720_STATIC_INGEST` + `SPEC_LEARNING_MODEL` + `doc/thesis-overview/PRECIS_ADVISOR` (the two-model arch).

## THE PLAN — two tracks, can interleave (A = codex/science, B = interactive UI)

### A. The combined 1P9J+720 Model-1 (the empirical Stage-3 test — "or we have no idea what we have")
1. **Integrate the 720.** Merge `wt-720-build` (the built static-pose ingest) → `h5-reader-pysr-spike`.
   Prep first: merge-readiness/conflict check (tracked io/model/rediscover changes vs the main tree);
   a POST-merge reader build re-verify (pre-merge baseline PASSED, `b0gpfdrxr`). The nanoflann symlink is
   MOOT (untracked, not merged; main tree has the real header). Open-spec calls recorded: cohort-manifest
   shape, thin-support threshold 10. **The lead merges.**
2. **Run the ingest on the FROZEN data — NO re-run.** Frozen 720 NPYs at
   `/shared/2026Thesis/nmr-shielding/calibration/features/Stage1BMRB_20260513_topology` — 720 complete,
   MOPAC + **absolute** ORCA DFT shielding, dated 2026-05-13/14. **VERSION-SURGERY (the lead's "version
   pass"):** the frozen data predates current EFG/schema/static-ingest conventions → a small **old-schema
   exporter/adapter** converts it to the current all-atom absolute-shielding substrate, recording old-schema
   provenance (see FRESH_LOOK for the exact field/schema deltas). Re-run (`python3 learn/extract.py --run
   <new> --resume --stage1-audit`, MOPAC ON, from the repo root) only for the exact current schema or a real
   calculator bug. **Target = the RAW absolute DFT shielding** (the ingest's raw dia/para/total), **NOT** the
   derived `target_matrix.npy` (that is the Stage-1 mutation-delta).
3. **Make the combined model — honestly.** Union the 1P9J substrate (the interp v0's data) + the 720
   substrate → train the **e3nn** shielding emulator (Model 1) → report held-out **shielding** recovery:
   within-1P9J AND **across the 720** (the transferability we've had no instrument for). Extend the existing
   `analysis/equiv_t2_e3nn.py` (the interp v0 IS Model 1's first instance). **e3nn, NOT MACE** (lit scan: e3nn
   is the direct-precedent SOTA for shielding tensors; MACE-on-shielding is unprecedented). Honest first-pass
   report (the stale-data caveat). "A model that shows a result need not be good." Should be short if lucky.
   - Two-model arch (context): Model 1 = shielding emulator (DFT-trained, this step). Model 2 = shift
     predictor (BMRB-trained on the no-DFT ~656 fleet, eats Model 1's firehose, ablate) — LATER.

### B. The reader UI (desk-ready for the advisor)
- **Confirm runtime:** the reader build PASSED (no CLI breakage); the lead should RUN it (load 1P9J) to confirm.
- **Open 1P9J + not-suck — infelicities, in concert** (the overview list ≈ the lead's): debug-flavoured
  toolbar; "Instrument" user-facing label; no selection Clear; per-frame overlay log noise; the default
  `dssp_chi` dashboard signal; **`QtFrame::eeqCharge` returns placeholder 0** (don't ship if EEQ shows). The lead drives.
- **The RADIUS view (headline, wanted badly):** pick atom + radius → show only that atom + the residues
  within R. Clear narrow path (UI_STATE_OVERVIEW §2): a `displayMolecule_` + `AtomFilter::WithinRadiusOf(atom,
  R, residue-expanded)` + source↔display index map + isolate-toggle/radius-spinbox/clear UI + optional REST.
  Builds on the stabilisation infra (`TransformedConformation` + `CameraComposer`). Risk: picking over the
  filtered molecule — pick over source positions + the map.
- **Linux install:** AppImage (`install(TARGETS)` + `$ORIGIN` rpath + bundle Qt/VTK/HDF5), tested on a clean
  Ubuntu (the lead is prepping a machine). Mostly-easy; Linux = open-source Qt (Pro licences are Win/Mac).
- Reader code: the **`qt6-cpp` skill** + CENSUS/ACONNECT/ASSERT_THREAD + UDP-log-first + one-persistent-QTimer.

## PARKED / LATER
- Maths verdict re-runs (LOAO re-fit + e3nn clean re-run) — parked; maths is lead+Claude.
- Interp `equiv_t2_e3nn.py` build-diff (+441/−49) — verify it didn't disturb the clean-protocol exemplar; consider a derived script.
- Model 2 (shift predictor) — after Model 1. Live-reader-prediction / the chewer killer-app — tabled
  (`NOTES_INTERP_LIVE_PREDICTION`), the lead wants it back.
- The advisor précis — the lead finalises (`doc/thesis-overview/`).

## DISCIPLINES
Branch never merge/switch; lead owns ALL git; codex grinds + opus reviews + lead vets; no symlinks;
model-is-spine / no Python second model / lean disk <15 G / extractor SACRED; e3nn not hand-rolled; T2 sacred.

## GIT STATE (as of `dfc8f51`)
Decruft + session docs committed (`dfc8f51`). **Uncommitted in the working tree** (lead's to handle): the
interp `analysis/equiv_t2_e3nn.py` + `analysis/interp_1p9j_graph.py` (the interp churn/graph); `doc/calculators/*`
(concurrent effort); `h5-reader/CLAUDE.md` + `.gitignore` + `CMakeLists.txt` (pre-existing); the new
`LEARN_720_RECIPE_AND_DATA` + `FRESH_LOOK` reports. The `wt-720-build` worktree (the built 720 ingest) awaits merge.
