# UI surgery — 2026-06-05 (running log; document-as-we-go)

The reader is the advisor-facing vetting surface; the lead's read is that the 3D + display have "real real
issues" and that UI was "added whenever it felt good" — accreted, not designed. This is the running record of
the coherence surgery: decisions, diagnoses, the cut. Companion to `POLISH_BACKLOG.md` (broad debt),
`UI_STATE_OVERVIEW_2026-06-04.md` (survey), and the pending `UI_COHERENCE_REVIEW_2026-06-05.md` (codex map).

## Organizing principle (lead, 2026-06-05)
**No buttons — or fields — for what we don't need.** Strip the UI to what's actually used. This is the lens
for the whole cut: the tumbling transforms, the empty data fields, the accreted toolbar all go.

## Decisions (made)
- **Transforms → collapse to ONE switch: backbone-fit ↔ all-atom-fit. Identity + Center-COM removed.** Both
  fits earn their keep (all-atom = full-structure alignment; backbone = clean sidechain-motion viewing), so
  one switch stays — the minimal *needed* switch, not cruft. The tumbling modes (Identity, Center-COM) come
  out. **Default on open: backbone-fit** (rock-stable backbone, clean dynamics); flip to all-atom on the
  switch. No tumbling, ever. (Lead 2026-06-05: "Fit reference is the only default we need… we do not need
  tumbling here"; "we don't need buttons for what we don't need"; "we need both so there is one switch left.")
- **Clean startup — DONE.** Inspector/Selection/Strip docks no longer open on launch (`ReaderMainWindow`:
  `setVisible(false)` ×3 + `kSettingsVersion` 1→2; window geometry still persists). See `POLISH_BACKLOG.md`.

## Diagnoses (findings)
- **Kabsch is textbook-correct + safe** (`FitTargetMath.h::ComputeSubsetTransform`). Centroids → centered
  covariance `H = P·Qᵀ` → SVD → `R = V·diag(1,1,sign det)·Uᵀ` → `T = cr − R·cc`; rank + `det(R)=1` guards. One
  R, one T applied uniformly to every atom → **internal flexing fully preserved; it cannot hold atoms in
  place** (a rigid transform never changes relative positions). All-atom Fit-reference never goes
  rank-degenerate, so the degenerate-freeze fallback never fires. → Fit-reference is safe as the permanent mode.
- **What each transform did** (for the record): Identity = raw (drift + tumble); Center-COM = de-drift, still
  tumbles (R = identity, T = −COM); Fit-reference = Kabsch-align all atoms to a reference frame (de-drift +
  de-tumble); Fit-backbone = Kabsch over the backbone subset. Coherent as a set, but unlabeled + inert on
  static poses + the display-vs-source split below. We keep only Fit-reference.
- **Clipping bug — ROOT-CAUSED.** `MoleculeScene` line 136 turns VTK auto-clipping OFF (correctly:
  `vtkActor::GetBounds()` is stale, pinned at frame-0 per `feedback_vtk_bounds_cache`). Clipping is reset by
  hand (`ResetCameraClippingRange(current-frame bounds + 5 Å)`) only on **frame change** (`setFrame`:428) and
  **plane-lock** (:499). The **zoom/dolly** path (`CameraInputFilter::handleWheel`:141 + drag-dolly:122) only
  `requestRender(CameraInput)` — **no clipping refresh** — so the near-clip stays at the old (farther) camera
  distance and slices the molecule as you push in. That is exactly "the box cuts off when you zoom in and the
  protein disappears." **Fix:** re-derive near/far in the **render path** (cache per-frame bounds; reset
  clipping every render from the live camera) so no gesture can skip it.
- **Display-vs-source coherence (the central issue).** The transform moves the **displayed** positions
  (`TransformedConformation::atomPosition` = `R·raw + T`) but NOT the calculator/snapshot/H5 data
  (`loadSnapshot` forwards source-frame snapshots unchanged, by design). So any overlay / picker / marker that
  reads SOURCE positions desyncs from the displayed (aligned) molecule. With Fit-reference now ALWAYS on, this
  goes from "only when toggled" to "always load-bearing." **Codex is mapping exactly which consumers read
  source vs the transformed seam** (`bfd0llv22` → `UI_COHERENCE_REVIEW_2026-06-05.md`).

## The cut (surgery plan)
1. **Fit-reference permanent; remove the transform menu** (`ReaderMainWindow` transform actions out;
   `TransformedConformation` default → `FitReference`, reference frame 0).
2. **Unify the display pipeline on the aligned positions** — overlays + picker + measurement/reveal markers
   all read the transformed seam, so always-on Fit-reference is coherent. Scope from codex's map.
3. **Clipping render-path fix** — re-derive clipping every render.
4. **Strip empty fields** — detect empty/sentinel/all-zero fields on load; show nothing with no data
   (lead's rule). Mechanism from codex's map.
5. **Pane toggles** — a clear, obvious way to bring the docks back (post the hide-on-startup change).
6. **Cut accreted cruft** — the Cruft Register (`QtSelectionOverlay` dormant, `/proc/self/statm` in render
   code, B-field `dynamic_cast`, render-drop PROBE, default `dssp_chi`, "Instrument" label, …).

## Codex coherence map — RECEIVED 2026-06-05 (`UI_COHERENCE_REVIEW_2026-06-05.md`)
- **Relief on the scary part:** the display pipeline is ALREADY unified on the `atomPosition` seam — scene,
  overlays, picker, REST `/positions` all read the wrapper polymorphically. So Fit-permanent is coherent out
  of the box; the only source-vs-display divergence is snapshot/calculator *data* staying in source frame (by
  design) → needs a LABEL, not a pipeline rebuild. (Step 2 of the cut shrinks to labeling.)
- **Clipping:** confirms the root-cause + fix exactly — a scene-owned `syncCameraClippingRange()` from cached
  per-frame bounds, called after EVERY camera write; the plane-lock special-case collapses into it.
- **Empty fields → the mechanism:** a shared **loaded-run availability layer** built on load (Absent /
  NoFramePayload / AllMissing / AllZeroStructural / AllZeroObserved / Available), consumed by the inspector
  (adders return `bool` → no dash-only groups) AND the metric catalog (offer only available channels) + gate
  the default `dssp_chi`. One shared layer, not per-panel filters. This is "find the empty stuff on load."
- **Panes:** the toggle buttons EXIST (relabeled `toggleViewAction`s) but are buried at the dense toolbar's
  far end with no View menu → add **View → Panels** + clearer grouping (same `QAction`s).
- **RESOLVED (lead 2026-06-05):** need BOTH fits → one switch backbone-fit ↔ all-atom-fit, default backbone;
  tumbling (Identity, Center-COM) out.

## Pending
- **Codex implements the full cut in one pass** (review-as-spec + the lead's calls); Claude reviews the diff;
  rebuild + relaunch; lead tests.
- Lead's pace call (2026-06-05): **cut everything together** — no piecemeal clipping fix.

## Log
- 2026-06-05 ~00:xx — created. Transforms decided (Fit-reference permanent, menu out); Kabsch verified;
  clipping root-caused; display-vs-source coherence identified as the central surgery. Codex coherence review
  in flight.
- 2026-06-05 — **codex executed the full cut** (`bkl0l0uni`, build green): clipping fix, availability layer,
  View→Panels, transform one-switch (backbone default + "All-atom fit" toggle), cruft cut. (The `equiv_t2_e3nn.py`
  churn in the diff is the *prior interp work*, not this surgery.)
- 2026-06-05 — **Claude review of the availability layer:** conservative — hides only Absent / NoFramePayload /
  AllMissing; **all-zero-observed is SHOWN**, so it can't hide real data. Flags: (a) because all-zero-observed
  shows + `structuralZeroDescriptor()` is a stub (→false), it may NOT hide the `eeqCharge=0` / always-zero
  fields the lead wanted gone — a judgment call (is "all zeros" empty?); (b) `classifyFramePayloads` scans
  EVERY sampled snapshot at load (eager, not lazy) → possible slow open on many-frame runs.
- 2026-06-05 — **Opus adversarial review in flight** (`a69…`). Next: synthesize agent + Claude findings → fix
  anything real → rebuild → relaunch the vetted surgery for the lead to test.
- 2026-06-05 — **Opus adversarial verdict: high quality, no CRITICALs.** Confirmed correct: clipping
  (render-path sync covers every camera write by construction), transform switch (exhaustive, no dangling
  modes, mode not persisted), View→Panels (shared QActions), cruft cuts; availability over-hide logic
  fail-safe. One HIGH (**H1**): eager all-sampled-snapshot scan at open → **FIXED** (`TrajectoryFieldAvailability.h`,
  cap at 4 representative rows — presence is uniform per-frame). Non-blocking follow-ups: **M1** dead
  `npy:mc_scalars` inspector gate id (fail-safe); **M2** transform-action-sync emitted before the sync slot
  connects (works, incidental); **L1** `QtSelectionOverlay` unwired but still built by the top-level
  `CMakeLists` (outside codex's edit boundary).
- 2026-06-05 — **Relaunched the vetted + H1-fixed surgery** (`b7fi9td52`): opens clean, 846 atoms, playing,
  interactive — H1 fix held (no open stall). Lead testing. (A too-broad `pkill` self-matched the prior
  relaunch command and killed it; clean relaunch fixed it.)
- 2026-06-05 — **lead live-test feedback, round 1** (4 UI-polish issues): (1) one toolbar holds everything →
  split into **Playback + Tools** rows; (2) **Panels toggle shows nothing** — docks open at ~0 width
  (`minimumSize(0,0)` override + discarded saved width); (3) **"Free" appears twice** (the `Free` radio + a
  redundant `cameraStateLabel`); (4) the default `npy:dssp_chi` metric seed is always pre-selected → **no
  clean unselected state**. → **codex on all four** (`b21gcvvvi`, brief `CODEX_BRIEF_UI_TOOLBAR_FIXES_2026-06-05.md`),
  with exact line-level diagnoses. Claude reviews → rebuild → relaunch.
- 2026-06-05 — **round-1 fixes landed** (`b21gcvvvi`, build-verified, diff clean): toolbar split (Playback +
  Tools rows); panels resize-to-360 + raise on show (minimumSize(0,0) override removed); `cameraStateLabel_`
  removed (Free radio kept); `npy:dssp_chi` seed removed + metric dialog opens unselected
  (`SignalDisplayDialog` `clearSelection`/`clearCurrentIndex` instead of `selectRow(0)`). Relaunched
  (`bfs1vhs5f`) — up clean, 846 atoms, lead re-testing.
- 2026-06-05 — **lead: panels STILL don't appear** ("I select a field, I add it, and alas no panel") — the
  targeted dock-show fix did not resolve the Add→panel flow. → fired a **read-only comprehensive UI-paths
  review** (`bpgh47r1x`, brief `CODEX_BRIEF_UI_PATHS_REVIEW_2026-06-05.md`) mapping the metric→panel→dock
  visibility path + the others end-to-end; we fix from the map, not another guess. Leading hypothesis: the
  Strip dock starts hidden and nothing shows/raises it on Add, so the panel renders invisibly.
- 2026-06-05 — **UI-paths map back** (`bpgh47r1x` → `UI_PATHS_REVIEW_2026-06-05.md`). Root cause CONFIRMED at
  code level: add-metric builds the panel but it renders into the **hidden Strip dock, never shown on Add** →
  invisible. Other funk: some Add modes (`static.tensor/table/atomColor/...`) produce no visible surface;
  Focus/Newman camera buttons don't apply until the next frame (plane-lock does); pick updates the inspector
  but doesn't reveal it. Working paths (camera gestures, transform, overlays, pick→marker) are sound.
- 2026-06-05 — **Lead frame for the eval:** *"things mostly work; it is just too funky to hand to someone
  because too many things don't"* → the bar is **HANDABLE, not perfect**; fix the handability blockers, not
  depth. → fired **opus check of the map** (`ac6…`, handability-lensed) per the lead's call. Sequence after it
  clears: **checkpoint commit (lead's git)** → apply de-funk fixes → rebuild → relaunch.
- 2026-06-05 — **opus check verdict: GO** (`ac6…`) — every diagnosis confirmed against the code; refinement:
  the dead-end modes are MORE pervasive than the map said (`static.table` is offered + enabled on nearly every
  descriptor and renders nothing), so **"disable show-nothing modes" is CO-EQUAL with the dock fix**, not
  optional. Risk low + contained (UI/dialog/catalog only).
- 2026-06-05 — **checkpoint committed `7b3012d`** (vetted surgery + fixes + session docs; restore point before
  the next fixes). NB git-process: even with permission, ASK/show scope before committing — lead correction.
- 2026-06-05 — **lead revision + go:** do NOT auto-pop the Strip dock on Add — make the dock-SELECT (Panels
  toggle) reliably work instead (root cause per opus: `DashboardStripDock` has no `setMinimumWidth` → collapses
  to ~0). → codex on 3 fixes (`bgcjgs1y7`, brief `CODEX_BRIEF_UI_HANDABILITY_FIXES_2026-06-05.md`):
  (1) dock-select reliable — min-width + queued resize/raise, no pop; (2) disable dead-end mode checkboxes;
  (3) Focus/Newman immediate. Defer pick→inspector. **Gate relaunch on a runtime confirm the Strip dock opens
  at ~360 wide.**
- 2026-06-05 — **handability fixes landed + relaunched** (`bgcjgs1y7`): dock min-width 260 + queued resize/raise
  (no auto-pop), dead-mode checkboxes gated by `modeHasVisibleSurface`, Focus/Newman immediate.
- 2026-06-05 — **lead: panels STILL greyed/empty** ("selected another field, still no panels enabled… the
  panel selector is constantly greyed out"). Reframe: the **selected-metrics MODEL is murky** — UI greys
  disconnected from reality (silent failure). Lead's target model: ONE selected-metrics list (source of truth);
  add→grow, close→shrink; **panel up-by-default unless the user disabled it**; selector REFLECTS the list. PLUS
  **ultrathink the relation to DATA AVAILABILITY** (`TrajectoryFieldAvailability` as the gate — untangle
  "display mode has no renderer" from "data unavailable"; the dead-mode-greying conflated them) **and the
  model-is-spine UI DATA-FLOW** (H5→typed model→views; selection = model state, views derive). → codex
  **design review** (READ-ONLY, qt-skill-loaded, portable Windows/macOS first-class) `bshk32hbp` →
  `SELECTED_METRICS_DESIGN_2026-06-05.md`. Then **opus-check → build the unified model**. Git: ask first.
  Checkpoint floor `7b3012d`.
- 2026-06-05 — **opus check: GO, staged** (`a2fb…`); diagnosis CORRECT, STEERABLE not rebuild. **Greyed selector
  = a ROUND-1 SELF-INFLICTED regression** (the `clearSelection`/`clearCurrentIndex` "clean-unselected" line we
  added is what greys all mode boxes) → fix: auto-select the PICKER candidate while keeping the ADDED-metrics
  list empty (both asks hold). **Latent bug:** `isPanelMode` (controller) vs `isPanelDisplayMode` (panel-model)
  ALREADY disagree on `static.tensor` (intentional — scene overlay); the capability-table unify must PRESERVE
  that (`isPanelRef` vs `hasStripWidget`), not flatten.
- 2026-06-05 — **DECISION (lead): full coherence is the rescue** (band-aids won't hold; the thrash came from
  treating view-state as model-truth). **NEXT SESSION → `NEXT_SESSION_UI_HANDOFF_2026-06-05.md`** — carries the
  deep-dive mandate, the instrument-REST/logging-first directive, the opus-ordered 5-step build
  (selector-reflection[LOW] → capability-table[LOW] → availability-through-path[MED] → dock-visibility[MED] →
  named controller[HIGH, LAST]), the 2 product calls (dock-visibility rule; prune-vs-keep stale data), and
  state. Floor `7b3012d`; ask before git.
