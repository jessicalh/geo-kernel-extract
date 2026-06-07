# Fresh-Look Prep — a grumpy-but-helpful cold-read assessment (2026-06-05)

**Who wrote this:** a pretend brand-new Claude session, opened cold tomorrow morning, no
memory of today, told to "continue serious work." I read the docs in newcomer order, then
checked the load-bearing claims against the code. Verdict up front: **the code is in good
shape and the technical doc descriptions are accurate; the problem is almost entirely
*staleness of status* and *doc sprawl* — the docs describe yesterday's git/process reality,
and there is no single "start here tomorrow" pointer.** Every gripe below has a concrete fix.

Read-only pass. I changed no code and ran no git mutations.

---

## TL;DR — the 3 things that would burn a fresh session, and the minimal fix

1. **The spine doc lies about git state.** `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md`
   says Stage 1 is "**BUILT + opus SHIP, UNCOMMITTED**" and "**Uncommitted — lead owns
   git**" (lines 26, 32) and the stabilisation floor is `33c2f50` (line 9). **Reality:
   `a6193dd` is committed AND pushed** — it contains the entire Stage-1 registry
   (`VisualizationRegistry`, all 7 concrete definitions, `DisplayModeCapability.cpp`) *and*
   the 67-line update to this very spine doc. So a newcomer reads "uncommitted, go commit
   it," runs `git status`, finds a clean tree, and wastes 20 minutes reconciling.
   **Fix:** flip the Stage-1 line to "COMMITTED `a6193dd` (pushed)"; bump the floor from
   `33c2f50` to `a6193dd`.

2. **There is no single resume pointer.** `h5-reader/notes/` has **40 files all dated
   2026-06-05.** The CLAUDE.md sends you to `notes/POLISH_BACKLOG.md` and `notes/SCOPE.md`
   (neither is the live spine), and the project CLAUDE.md doesn't name today's spine doc at
   all. The actual live doc is `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md`, but you only
   learn that by guessing or by `ls -t`. **Fix:** one line at the top of
   `h5-reader/CLAUDE.md`: "Live reader state + next step: `notes/BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md`
   (read its CURRENT STATE block at the very top; ignore the SUPERSEDED mid-sections)."

3. **The spine doc is a 22 KB stratigraphic dig** where the *bottom* is newest. Five
   "CURRENT/SESSION END/DUCK-WALK/BUILT" layers stack up, and only the reader who already
   knows the punchline can tell that the top block (line 3) and the very last block
   (line 228, "COMMITTED `33c2f50`") win, while the middle ("UNCOMMITTED temporal
   stabilisation," "duck walks at window=32," lines 156-193) is dead history. The
   `[SUPERSEDED]` banner at line 90 is good practice — but it's the *only* one, and three
   other stale blocks have no such banner. **Fix:** see the truthfulness table below; banner
   or delete the stale blocks.

---

## 1. Orientation — could a newcomer find where things are and what state we're in?

**Mostly yes, eventually, but with avoidable detours.**

- The project `CLAUDE.md` is excellent and current on the big picture (5-mode spec, three
  stages, subproject map). No complaints there.
- `h5-reader/CLAUDE.md` is solid on *rules* (read-only, model-is-spine, Qt discipline) but
  **points at stale entry docs**: it cites `notes/POLISH_BACKLOG.md`, `notes/SCOPE.md`, and
  the project CLAUDE.md cites `notes/POLISH_BACKLOG.md` as "next sessions." None of those is
  the live work doc. A newcomer following the documented breadcrumb does **not** arrive at
  `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md`.
- **First thing that confuses:** the sheer count of same-date docs with no manifest. There's
  a "Doc map" inside the spine doc (line 53) but it's incomplete (lists ~12 of the 40 notes;
  omits e.g. `DESK_READY_PUNCHLIST`, `UI_SURGERY`, `NEXT_SESSION_UI_HANDOFF`,
  `SELECTED_METRICS_DESIGN`, all the `CODEX_BRIEF_*`) and you have to already be reading the
  spine doc to find it.
- **Where you'd waste time:** trying to figure out whether Stage 1 needs committing (it
  doesn't), and whether the "temporal stabilisation / duck-walk" saga is live work (it's
  done and shipped). Both are answered correctly *somewhere* in the spine doc, just not
  where you look first.

**The good news a newcomer would also discover:** once you find the spine doc's top block
and its final "BUILT + VERIFIED" block, the technical narrative is honest and matches the
code. The `[SUPERSEDED]` marker on the windowed-stabilisation section (line 90) is exactly
the right move and should be the template for the rest.

---

## 2. TRUTHFULNESS FLAGS (the important section) — `file:line`, doc claim vs reality

Paths are relative to `h5-reader/notes/` unless noted. I verified each against the code/disk.

| # | Where | Doc claims | Reality (verified) | Fix |
|---|---|---|---|---|
| **T1** | `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md:26` and `:30-32` | Stage 1 registry "BUILT + opus SHIP, **UNCOMMITTED**"; "**Uncommitted — lead owns git.** Before/with the commit: the 2 opus test-asks + the empty-check refinements." | **Committed + pushed as `a6193dd`.** The commit ships `VisualizationRegistry.{h,cpp}`, `VisualizationDefinition.{h,cpp}`, `DisplayModeCapability.cpp`, and all 7 concrete defs. Working tree is clean (only untracked `doc/searchresults1/`). | Change to "COMMITTED `a6193dd` (pushed)." Note which of the "2 opus test-asks + empty-check refinements" actually made it in vs deferred — the commit landed without an explicit confirmation here (see T9). |
| **T2** | `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md:9-10` | "Floor for new work is now `33c2f50`, not `141e7b6`." | Floor is now **`a6193dd`** (two commits past `33c2f50`: `c0cce8c` rediscover, then `a6193dd` registry). | Bump the stated floor to `a6193dd`. |
| **T3** | `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md:156-193` (SESSION END STATE block) | "**UNCOMMITTED in the working tree:** the temporal stabilisation (`TransformedConformation` precompute + quaternion-smoothed (R,T))… NOT committed." + "Live now: reader on `:1`… window=32 (duck-walking). `141e7b6` floor." | Fully **superseded and committed**. The windowed (R,T) smoothing was the *bug*; the shipped fix is iterative-mean + centroid-pin in `33c2f50`. The "duck-walking window=32" live state is two commits stale. | This whole block needs a `[SUPERSEDED — see top + the COMMITTED 33c2f50 block]` banner like line 90 has, or deletion. As written it reads as current and points at the wrong floor + a broken live instance. |
| **T4** | `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md:126-139` (Build plan) | "Floor before this work: `141e7b6`." Steps 1-4 framed as the upcoming plan; step 1 "Temporal stabilisation… Tune the smoothing window on `:1`." | Step 1 (stabilisation) is **done and shipped** (`33c2f50`). Floor is `a6193dd`. Steps 2-4 (named two-mode switch, depth-aware pick, the big cut) are genuinely still open — but the framing buries that. | Mark step 1 done; restate the real remaining work (2-4) under the top CURRENT STATE so the open items aren't hidden under a stale floor line. |
| **T5** | `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md:61,139,193` | Three separate "floor = `141e7b6`" assertions across the doc. | `141e7b6` is 4 commits back. | Single source of truth for the floor at the top only; strip the repeated stale ones. |
| **T6** | `DFT_LGS_REGISTRATION_SCOPE_2026-06-05.md` — header `:5` says "**SCOPE / DESIGN only. No code, no git, no build in this pass.**" but the spine doc `:41-51` says the fix is **done**. | The note reads as scope-only; the spine doc says "SCOPED + FIXED," "660→751 DFT frames `[0..1500]`," "lean LGS stashed → `.LGS.back`." | **The fix landed on disk** (verified): the calcset dir now holds exactly one `.LGS` (`…dense-mopac-live-orca.LGS`) + a `…with-dft.LGS.back`; the dense LGS has **752 `frame_index` lines** (751 frames) with the last at `f001500`. So the spine doc is right, but the SCOPE note's own status header is now stale/misleading if read standalone. | Add a one-line "**STATUS UPDATE (later 2026-06-05): executed — 751 frames registered `[0..1500]`, Gap B resolved (one `.LGS` + `.back`); see spine doc.**" to the top of the SCOPE note. |
| **T7** | `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md:50-51` | "The empty-check still false-flags `orca_dft` until the Stage-2 sidecar fix." | **True and correctly stale-flagged.** Verified `collectExpectedButEmpty()` (`src/app/DashboardDisplayController.cpp:2438-2482`) walks `allDescriptorList()` and consults only `recordForDescriptor`/`recordForStoragePath` — it is NOT run-mode aware and does NOT consult `dft.frames[]`. So `orca_dft` (a sidecar, absent from the H5 availability table) and mutant-only NPYs on a trajectory will both still false-flag. | No fix needed — this one is honest. Keep it. (Good example of the discipline the stale blocks lack.) |
| **T8** | `STAGE1_REGISTRY_PLAN_REVISED_2026-06-05.md` — header `:1-5` "Plan only. No code changed." | It's a plan. | The plan was **executed** and matches the code closely (registry, out-of-line gate, 6 static rows, `strip.` prefix arm, `visibleOfferable`/`trackedButHidden`/`unresolvedStaticModes`, `recordForStoragePath` accessor, `StripComponent` enum, `TensorGlyphVisualization` inert via `tensorGlyphGestureEnabled=false`). All present in `src/model/`. | Add a one-line "**IMPLEMENTED in `a6193dd`** — see `src/model/VisualizationRegistry.cpp`" so a newcomer doesn't re-plan a done thing. |
| **T9** | `BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md:30-32` and `:257` | "Before/with the commit: the 2 opus test-asks (`visibleOfferable`/`trackedButHidden`/`collectExpectedButEmpty` coverage) + the empty-check refinements." The codex deviations note (`:257`) lists what landed but doesn't confirm the test-asks. | The methods (`visibleOfferable`, `trackedButHidden`, `collectExpectedButEmpty`) all **exist** in code. Whether dedicated *tests* for the tri-state coverage were added — and whether the "empty-check refinements" (run-mode / sidecar awareness) were folded in — is **NOT** confirmed by the committed-deviations note, and T7 shows the sidecar refinement is explicitly still pending. So at least part of the "before the commit" checklist was **deferred past the commit**, not done before it. | State plainly: the registry shipped in `a6193dd`; the empty-check run-mode/sidecar refinements are **Stage 2, not yet done**. Confirm separately whether the tri-state coverage tests exist (I did not find them by name; worth a `ctest -N \| grep -i visualiz` before claiming). |
| **T10** | `learn/CLAUDE.md:65-70, 101` | "Weighted R² = **0.718** (720 proteins)" in one place but "**Do not try to improve R². 0.818 is the answer.**" two lines later, and `:74` "Was 0.818 on 110 proteins." | Both numbers are real (0.718 full-720, 0.818 110-protein fair set — the project CLAUDE.md states this clearly). But `learn/CLAUDE.md:101` baldly says "0.818 is the answer" with no qualifier, contradicting its own headline 0.718. A newcomer doesn't know which is "the thesis number." | One-word fix: `:101` → "**0.718 (full 720) is the thesis number; 0.818 was the 110-protein fair set.** Do not try to improve it." Align with project CLAUDE.md. |
| **T11** | `learn/CLAUDE.md:14, 134-136` | "calibration/{ID}/ (**723** proteins)"; "AzimuthalExtraction/ — current extraction (**112/723 done**)"; "GatedCalibration/ — previous full extraction (723 proteins)." | The project CLAUDE.md and memory say the settled number is **720** proteins / 446K atoms (effective fleet count 676 elsewhere). `learn/CLAUDE.md` still says 723 and references a half-finished "AzimuthalExtraction 112/723" that reads as live but is almost certainly an abandoned scratch run. | Reconcile to 720 (or explain the 723→720 drop once). Mark AzimuthalExtraction dead or delete the line; a newcomer will otherwise think there's a 112/723 extraction to resume. |

**Net truthfulness read:** the *code-level* technical claims (DFT bridge mechanics,
stabilisation math, registry shape, empty-check behaviour) are **accurate** — I checked them
and they hold. The lies are all in the *status/process layer*: git committed-vs-uncommitted,
which floor, which session-end live state, and a couple of drifted protein counts in
`learn/`.

---

## 3. Up-to-date — what status lags reality

- **Stage 1 reader registry: built + committed + pushed (`a6193dd`).** Docs say
  uncommitted. (T1, T2.)
- **Stabilisation: shipped (`33c2f50`), default window 0, centroid invariant proven.** The
  mid-doc "uncommitted / duck-walking / floor 141e7b6" blocks lag by 4 commits. (T3-T5.)
- **DFT registration: executed.** 660→751 frames `[0..1500]` now in the dense `.LGS`; Gap B
  (two `.LGS` in one dir) resolved to one `.LGS` + one `.LGS.back`. Verified on disk. The
  SCOPE note still reads "design only." (T6.)
- **The overnight/720 run:** the 720-WT static work lives in a **separate worktree at
  `/shared/2026Thesis/nmr-720` on branch `wt-720-build`** — and that branch is at `04fcb3f`,
  i.e. **well behind** the reader branch's stabilisation+registry work. No reader doc tells a
  newcomer this worktree exists or where it is; you only find it via `git worktree list`. The
  rediscover docs (`src/rediscover/STATE.md` top synthesis, trued 2026-06-04) are the real
  home for that thread and are reasonably current, but the *pointer* from the reader side is
  missing. (See gap G3 below.)
- **DFT campaign itself: 100% (751/751).** Confirmed: `/shared/2026Thesis/1p9j-orcas/jobs`
  has 751 job dirs, last `…_f001500_t15000.0`. Spine doc and SCOPE note agree.

---

## 4. Gaps a newcomer NEEDS but isn't told

- **G1 — The LGS data model + `lgs_write.py`.** The DFT bridge hinges on `manifest.dft.frames[]`
  (the `frame_index → meta.json` map) and the registrar `tools/lgs_write.py`. Neither is in
  the reader CLAUDE.md or any "start here" doc; you have to read the 417-line SCOPE note to
  learn that `lgs_write.py --force --dft-jobs-dir …` is how the `.LGS` is regenerated, that it
  filters to `orca_exit_code==0`, and that it's idempotent. **A newcomer asked to "re-register
  the latest DFT frames" has no breadcrumb to the tool.** Fix: name `tools/lgs_write.py` and
  the `dft.frames[]` contract in the reader CLAUDE.md's "DFT" sentence.

- **G2 — The empty-check's two false-positive modes.** Verified in code (T7): a registered
  descriptor with **no availability record** is assumed present (`allowsDescriptor` returns
  `true` on a null record, `TrajectoryFieldAvailability.h:133-138`), so the check only fires on
  a *positive* empty state (`Absent/NoFramePayload/AllMissing`). The two things it currently
  *wrongly* flags: (a) `orca_dft:*` — a frame-local **sidecar** the H5 availability table never
  saw; (b) **mutant-only NPYs** (`npy:orca_*`, `wt/mut/delta_shielding`) on a trajectory run,
  which are honestly absent there. Both await the Stage-2 run-mode + sidecar awareness fix. A
  newcomer debugging "why is the reader screaming `viz_expected_but_empty` about `orca_dft`?"
  needs this stated as known/expected, not a bug to chase. Fix: a 3-line "Known false-positives
  (Stage-2 TODO)" note next to the empty-check description.

- **G3 — Where the 720 lives + that it's a separate, behind worktree.** See §3. Fix: one line
  in the spine doc or reader CLAUDE.md: "720-WT static work = worktree `/shared/2026Thesis/nmr-720`
  (`wt-720-build`, currently behind the reader branch); state in `src/rediscover/STATE.md`."

- **G4 — The mutant-eval vs trajectory-DFT shielding distinction.** The spine doc states it
  (line 49-51) but it's load-bearing and easy to miss: `npy:orca_*` / `npy:{wt,mut,delta}_shielding`
  are **MUTANT-eval-only** (mode-4) and honestly absent on a trajectory; `orca_dft:*` is the
  **trajectory** DFT (the H5/sidecar path, now registered). Mixing these up will send a newcomer
  hunting for trajectory shielding in the wrong arrays. Fix: keep it, but lift it into the
  reader CLAUDE.md "Scope" or a one-paragraph DFT subsection so it isn't buried in a session log.

- **G5 — The "tomorrow plan" itself.** The next concrete reader steps (Stage 2 Atom Colour +
  empty-check false-positive fixes; Stage 3 Tensor glyph ovaloid; named two-mode switch;
  depth-aware pick; the big cut of ~14 dead `static.*` modes) are real and scattered across the
  spine doc top block, the (stale-framed) build plan, and the VISUALISATION_AUDIT "DECISIONS
  OWED." There's no single ordered "do this next" list that isn't tangled with dead history.
  Fix: a 5-line "NEXT (ordered)" block at the very top of the spine doc, under CURRENT STATE.

---

## 5. The grumpy verdict — what would most slow me down, and the minimal fix

1. **"Go commit Stage 1" — but it's already pushed.** *(T1/T2.)* Fix: two-line status flip +
   floor bump to `a6193dd`. Highest-value single edit in this whole report.

2. **Four floor citations, three of them stale (`141e7b6`, `33c2f50`).** I can't trust *any*
   "floor" line without `git log`, which defeats the doc. Fix: one floor, at the top, = `a6193dd`.

3. **The spine doc's middle is a graveyard with one tombstone.** Only the windowed-stabilisation
   section is marked `[SUPERSEDED]`; the SESSION END STATE + build-plan + DUCK-WALK blocks read
   as live and point at a broken `:1` instance and the wrong floor. Fix: banner or delete blocks
   T3/T4; this doc wants a hard "everything below the line is history" rule.

4. **No front door.** `ls notes/` returns 40 same-date files and the CLAUDE.md breadcrumbs point
   elsewhere. I'd spend my first 15 minutes just deciding what to read. Fix: one pointer line in
   `h5-reader/CLAUDE.md` to the spine doc + "read only its top CURRENT STATE block."

5. **`learn/` number drift (723 vs 720, 0.818 vs 0.718, a phantom 112/723 extraction).** Minor
   for the reader thread, but if I wander into calibration I'll mis-cite the thesis number or try
   to resume a dead run. Fix: T10/T11 one-liners.

**Bottom line, said plainly:** nobody lied about the *physics* or the *code* — those are honest
and verified. What's stale is the bookkeeping: the docs were written across a long day and the
last two commits (`33c2f50` stabilisation, `a6193dd` registry) outran the prose. Flip ~6 status
lines, banner ~3 dead blocks, add ~4 pointer lines, and a fresh session lands clean instead of
spending its first hour reconciling the docs against `git log`.
