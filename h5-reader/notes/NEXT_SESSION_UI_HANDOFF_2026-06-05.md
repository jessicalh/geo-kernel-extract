# NEXT SESSION — h5-reader full-coherence build (deep-dive-first, instrument-first)

Paste to a fresh session started at `/shared/2026Thesis/nmr-shielding`. Branch `h5-reader-pysr-spike` — never
merge/switch/rebase. **The LEAD owns ALL git; ASK before ANY state-changing git even when it's fine ("no off
we go"); NEVER destructive git (a session once deleted a day's work).** Mechanism: codex grinds + opus reviews
+ lead vets. **Portable Qt6 project (Windows/macOS/Linux) — the Windows/macOS source is first-class.**

## STEP 0 — DEEP INVESTIGATORY DIVE (mandatory; the lead's explicit ask — do NOT skim a summary)
You must KNOW what last session knew — not because this doc told you, but because you READ the documents AND
CHECKED THE CODE and re-derived it yourself. Read IN FULL, in order, then verify each load-bearing claim
against the cited source BEFORE touching anything:
1. `h5-reader/notes/UI_SURGERY_2026-06-05.md` — the running surgery log (the spine; read the whole Log).
2. `h5-reader/notes/SELECTED_METRICS_DESIGN_2026-06-05.md` — the model design (the rescue plan).
3. `h5-reader/notes/UI_PATHS_REVIEW_2026-06-05.md` + `UI_COHERENCE_REVIEW_2026-06-05.md` — the path/coherence maps.
Then OPEN THE CODE and confirm with your own eyes:
- the **two-model split** — `DashboardSignalModel` (selected signals) vs `DashboardPanelModel` (placement refs),
  synced by the anonymous `DashboardSignalPanelCoordinator` in `ReaderMainWindow.cpp`;
- the **greyed selector** — `SignalDisplayDialog::onCandidateSelectionChanged` disables all mode boxes when no
  candidate row is current, after `refreshCatalog`'s `clearSelection`/`clearCurrentIndex` (← the round-1
  "clean unselected" change; same line, opposite asks);
- the **availability-vs-renderer conflation** — `modeHasVisibleSurface` (renderer) vs
  `TrajectoryFieldAvailability` (data), collapsed into one grey;
- the **`static.tensor` predicate disagreement** — `isPanelMode` (controller) vs `isPanelDisplayMode`
  (panel-model) already differ (intentional).
If your reading and the code disagree, the code wins — re-derive. Do NOT build on the summary alone.

## STEP 1 — INSTRUMENT FIRST (REST + logging) so you verify autonomously
The lead's unlock: with REST + logging hooked up properly you can take the whole build all the way WITHOUT
four human visual round-trips. The reader already has a loopback `RestServer` + UDP/stderr structured logs.
Extend them so YOU drive and observe the metric/panel UI headlessly:
- REST to: add a metric (descriptor + anchor + modes), remove a metric, LIST the selected-metrics, query dock
  visibility + width, query the rendered panel/strip state, set/clear dock visibility, toggle a display mode.
- Logging to: emit the selected-metrics list + dock visibility + panel-render outcome on every relevant
  transition, so a log tail proves "add X → list=[X] → dock visible@~360 → panel painted."
- Then each STEP-2 fix is verified by you (launch on `:1`, drive via REST, assert the log), not by asking the
  lead to look. Reader build/launch loop is in `POLISH_BACKLOG.md` (top dated block).

## STEP 2 — BUILD FULL COHERENCE (the lead's decision; opus-verified SAFE ORDER — each step ships)
Do NOT follow the design's own Fix-Plan order (it front-loads the riskiest controller work). Use THIS:
1. **[LOW] Selector reflection** — auto-select the PICKER candidate row when one exists (re-enables the mode
   boxes) WHILE keeping the ADDED-metrics list empty (the lead's clean start — a reconciliation, both asks
   hold). + split the disabled tooltip: "no renderer yet" vs "data unavailable". `SignalDisplayDialog`,
   ~20–30 lines, dialog-local. The lead's immediate pain — ship + verify FIRST.
2. **[LOW] Shared mode-capability table** — unify `modeHasVisibleSurface` + the two panel-mode predicates into
   one table. ⚠️ `isPanelMode` (controller) and `isPanelDisplayMode` (panel-model) ALREADY disagree on
   `static.tensor` (intentional — tensor is a scene overlay, deferred). PRESERVE the asymmetry (encode
   `isPanelRef` vs `hasStripWidget` as data); do NOT flatten to one bool or rendering silently breaks.
3. **[MED] Availability through the whole path** — add an availability role to `DashboardSignalModel`; validate
   adds against catalog/availability; carry availability on existing rows so a metric whose descriptor later
   returns null shows "unavailable" instead of silently vanishing (the current `if (!descriptor) continue;`
   drop in `DashboardDisplayController`).
4. **[MED] Dock-visibility derivation** — `dockVisible = selectedCount>0 && !userDisabled`. PRODUCT CALL owed
   (below) — get the lead's answer first.
5. **[HIGH] Named `DashboardSelectionController` / single source of truth** — promote `DashboardSignalModel`
   to the ONE selected-metrics authority; demote `DashboardPanelModel` to placement; move the anonymous
   coordinator's cleanup loop into the named controller (`addMetric`/`removeMetric`/`setMetricMode`/
   `removePanel`/`clearAll`); route the dialog's two-step writes through it. LAST, deliberately — the
   dog's-breakfast risk if rushed. Preserve the `static.tensor` asymmetry + the rebuild-race ordering
   (`DashboardDisplayController` ~:1229).
Each step: codex grinds → opus/you review → verify via REST+logs → relaunch only as needed. Floor `7b3012d`.
Ask the lead before any commit.

## PRODUCT CALLS — RESOLVED by the lead (2026-06-05)
- **Dock visibility — DECIDED:** on metric-add, if the dock is hidden, REVEAL it. "None? Now one? You get a
  dock." Any add-while-no-dock shows the dock (not just the first metric). That single add event is the ONLY
  trigger needed — no complex derived-visibility rule, no `userDisabled` suppression. (Supersedes the earlier
  "no auto-pop": auto-popup on the add-while-hidden event is exactly what's wanted.)
- **Run / `.LGS` load — DECIDED:** all-new. A run load resets the selected-metrics list + dock state (clean
  slate), so there is no mid-run "stale metric whose data vanished" prune problem — the selected list is
  per-run; reload is the reset point. Step-3 availability work still matters for what's OFFERABLE at add time.

## STATE
- **Checkpoint `7b3012d`** = vetted surgery (clipping, transform one-switch, availability layer, View→Panels,
  cruft) + handability fixes (dock min-width 260, dead-mode greying, Focus/Newman immediate) + docs. The
  restore FLOOR. Working tree on top of it: uncommitted handability fixes already in `7b3012d`; the
  `equiv_t2_e3nn.py` / `rediscover/` changes are SEPARATE (interp + 720 work), not the reader.
- **720 run** `Stage1BMRB_20260604` (check `calibration/features/Stage1BMRB_20260604/extract_log.jsonl`):
  fresh current-schema producer NPYs; ~44 early failures were the `nmr_extract` rebuild collision
  (`Permission denied`), recoverable via `--resume` against a stable `build/nmr_extract`. Combined-model /
  first-pass work is downstream + separate from the reader. See `720_RUN_STATE_AND_PLAN_2026-06-04.md`.

## DISCIPLINES
Branch never merge/switch; lead owns ALL git, ASK before any state-changing git, never destructive; codex
grinds + opus reviews + lead vets; **qt6-cpp skill before any Qt code** (CENSUS_REGISTER/ACONNECT/ASSERT_THREAD
/one-persistent-QTimer/UDP-log-first); **model-is-spine** (selection = model state, views DERIVE; no parallel
state machine on widgets); portable Qt (Windows/macOS first-class); T2 sacred; never open `trajectory.h5` in
Python.
