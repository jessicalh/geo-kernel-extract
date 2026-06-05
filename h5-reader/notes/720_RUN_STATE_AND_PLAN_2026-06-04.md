# 720 run — state + the C / Python / model plan (2026-06-04, evening)

**Authoritative current state.** Supersedes the version-surgery premise in `NEXT_SESSION_HANDOFF_2026-06-04.md`
(we re-ran fresh, not adapt-the-frozen-NPYs). Branch `h5-reader-pysr-spike`; lead owns ALL git; codex grinds +
opus reviews + lead vets; the maths/model discussion is lead+Claude.

## Where the 720 work IS

**Decision (2026-06-04):** re-run the 720 fresh rather than adapt the stale 2026-05-13/14 frozen NPYs. The
frozen data had drifted too far from current schema (no `.LGS`, 9-col EFG, AIMNet2 naming, FF14SB gaps); a
fresh run emits current-schema output the rediscover ingest consumes drop-in — the entire version-surgery is
gone.

**The run — COOKING.** `Stage1BMRB_20260604`, background task `b41z04q8e`:
- `python3 -u learn/extract.py --run Stage1BMRB_20260604 --resume --stage1-audit --config data/calculator_params.toml`
  (from repo root) → `calibration/features/Stage1BMRB_20260604/{ID}/`; logs `extract_background.log` +
  `extract_log.jsonl`.
- 720 complete pairs (3 missing-NMR excluded: A0A075FQU3, A0A7C4ZM98, A0A7J2L4W1). MOPAC ON (`--mutant`
  default-on; `--no-mopac` is the off-switch; `--mopac` is trajectory-only). AIMNet2 ON (extract.py self-sets
  `LD_LIBRARY_PATH` for the nvrtc/CUDA libs; model `data/models/aimnet2_wb97m_0.jpt`; RTX 5090 free).
- Pace ~30–45 s/protein → **~8–12 h** (not the 16 h ceiling); first ~50 land ~30 min in. Per-protein timeout
  600 s (large proteins fail boundedly, never hang).

**Validation — first 12 clean.** Validated against `h5-reader/src/io/QtFieldCatalog.gen.h`: no missing required
fields, correct shapes, EFG width-5, AIMNet2 present, non-degenerate, and **`total == dia+para` componentwise
within ORCA print tolerance**. The DFT target is consistent. Tolerance note: `total−(dia+para)` caps at
**1.633e-3** (identical across proteins = ORCA's 3-decimal print rounding, NOT corruption) on magnitudes
~700–1170 → the correct check tolerance is ~2–5e-3, never exact equality.

**Run status + known issue (2026-06-05 ~00:30).** ~198/720 logged; **154 OK, 44 FAILED.** Every failure is
the identical `[Errno 13] Permission denied: '…/build/nmr_extract'` — an **operational collision**, not a
calculator / pipeline / data bug: the producer binary was being **rebuilt** concurrently, and any protein that
exec'd `build/nmr_extract` during a relink window got `EACCES`. The original Stage-1 720 was 100%; this is
self-inflicted. **Zero data loss** — a failed protein wrote no `pos.npy`, so `--resume` re-runs exactly those
while skipping the good ones. **Plan (lead's call 2026-06-05): let it keep rolling**; finish the 720 issues
later with a final `--resume` pass against a stable binary. Robust alternative if rebuilds continue: snapshot
`build/nmr_extract` → a stable path, repoint the run, restart `--resume` (decouples the run from rebuilds).

**The run-tool.** `learn/extract.py` (+ docs + `learn/CLAUDE.md`) restored to the live tree from
`/mnt/expansion/...release-cleanup-20260528/learn` (it had been spinner-only), decrufted (current AIMNet2 field
names, complete-pair filtering, foreground + nohup wrapper). It self-configures the AIMNet2 `LD_LIBRARY_PATH`.
`learn/` is gitignored.

## How the rest works — C / Python / model, working off the first 50

**C++ ingest** (`/shared/2026Thesis/nmr-720`, branch `wt-720-build`, **UNCOMMITTED**): the `--static-cohort`
static-pose ingest that loads producer NPYs into the typed model and emits the per-(protein,atom) e3nn
substrate. Built, compiles, trajectory oracle passes. **The adversarial review
(`h5-reader/src/rediscover/REVIEW_worktree720_AND_MERGE_PLAN_2026-06-04.md`) says NOT merge-safe as-is:**
- CRITICAL — DFT target only trace-validated → add componentwise + element-order check **at ORCA-print
  tolerance (~2–5e-3)** (the first-12 validation pins this number).
- CRITICAL — strict topology soft-fails bonds/rings → enforce in `LoadPoseStrict`.
- HIGH — the oracle never exercises the static path → add a real static-vs-trajectory oracle.
- MED — FF14SB not loaded for statics (schema still emits the fields); 720-cardinality not enforced.
- LOW — nanoflann symlink → real vendoring.
- Bones right: T2 not scalar-collapsed, raw-ORCA-not-delta, all-atom labels, no Python model, disk gate <15 G.
- **Split:** mechanical fixes (componentwise/element-order, strict topology, nanoflann) = **codex**;
  static-oracle design + FF14SB call = **lead+Claude**. PENDING the lead's green.

**Python (e3nn first-pass):** extend `analysis/equiv_t2_e3nn.py` (the 1P9J interp v0) to the combined 1P9J +
first-50 substrate. Consumes the **C++-emitted substrate only** — never a Python scratch-read of NPYs; no
spatial work / no protein model in Python (model-is-spine).

**Model contract (lead+Claude):**
- Target = **raw absolute σ** (the DFT total), NOT delta, NOT de-meaned.
- **Discover within/between, don't force it** — no per-atom de-meaning, no architected split. Within (held-out
  1P9J frames) and between (held-out across proteins) are lenses measured *off* the trained model, not
  training constraints.
- e3nn (not MACE); T2-valued (σ_iso/T0 alongside as the scalar). Equivariant → no imposed local frame.
- The **first-50** gives the first across-protein transferability signal — the advisor down-payment graph. "A
  model that shows a result need not be good" — direction, not destination. "Not wired in" = a preliminary
  50-protein graph, distinct from the eventual full-720 model.

**Two-model arc (context):** Model 1 = this e3nn shielding emulator (DFT-trained) → Model 2 = shift predictor
(BMRB, no-DFT fleet, eats Model 1's firehose, ablate) — later.

## Consolidation (after ingest fixes) — lead's git
`REVIEW_worktree720_AND_MERGE_PLAN_2026-06-04.md` Part 2: an 8-step safe sequence to bring all scattered work
to one tree, nothing lost. The 720 build is uncommitted on `wt-720-build` (base `04fcb3f`, an ancestor of
`dfc8f51`) → commit-then-merge, conflict-free (code vs docs disjoint). Open lead decisions: Welford
agent-branches (keep `10e7dd1` / the uncommitted Eeq/HBond/Sasa edits?), dirty target docs, 720-cardinality
mode, worktree pruning.

## Pending / next
- [ ] Ingest CRITICAL fixes (codex mechanical + lead+Claude on oracle/FF14SB) — **pending green**.
- [ ] Fixed ingest on the first-50 → e3nn first-pass transferability graph (advisor down-payment).
- [ ] Model contract discussion (lead+Claude).
- [ ] Safe-merge to one tree (lead's git, after fixes).
- Parked: maths verdict re-runs; Model 2 (shift predictor); live-reader-prediction.
- **Active now: UI work** (radius view + infelicities + AppImage) while the run chews.
