# Codex — STEP 2.1: selector reflection (the greyed mode-boxes) + a picker probe

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`). **Portable Qt6
(Windows/macOS/Linux) — cross-platform source is first-class.** Scope: `h5-reader/src` ONLY. Never link the
`nmr_shielding` library, never write H5.

**Verify from the code before changing anything, ever.** Cites are `file:line` against the current working
tree; if one is stale, follow the code and tell me.

## Why this matters
The reader "mostly works but is too funky to hand off," and the worst of it: the Metrics dialog's display-mode
checkboxes are **greyed out disconnected from reality** — you pick a descriptor and the boxes stay dead. We
traced the root cause and it's a small, self-inflicted one. You're fixing it, and adding a tiny probe so the
fix is verifiable headlessly (no human at the screen). This is the first *visible* coherence fix; get it
clean and the user gets a selector that reflects what they picked.

## The root cause (verified in code — confirm, then fix)
The candidate (picker) side greys ALL mode boxes whenever no descriptor row is *current*, and two code paths
deliberately leave no row current:
- `SignalDisplayDialog::refreshCatalog()` clears selection + current index after loading descriptors
  (`src/app/SignalDisplayDialog.cpp:965-968`), then calls `onCandidateSelectionChanged()`.
- `SignalDisplayDialog::onAnchorSelectionChanged()` — when the candidate proxy has rows but the current index
  is invalid/unmapped, it ALSO clears selection + current index (`src/app/SignalDisplayDialog.cpp:1048-1055`),
  then calls `onCandidateSelectionChanged()`.
- `onCandidateSelectionChanged()` (`src/app/SignalDisplayDialog.cpp:1060`): with no current row, `record` is
  null → every candidate mode box gets `modeId=""` → `supported=false` → **disabled** with the single tooltip
  "This display mode does not have an implemented visible renderer" (`:1069-1076`). The auto-check fallback is
  guarded by `if (!checkedOne && record)` so it never fires either.

(These two clear sites are the prior session's "clean unselected" intent — they over-reached and greyed the
picker. The ADDED-metrics / active list staying empty is correct and must be preserved; only the PICKER side
is wrong.)

## PART A — the fix (dialog-local, ~25-40 lines)
1. **Auto-select the picker candidate row 0** at BOTH clear sites: when the candidate proxy has rows and there
   is no valid current candidate row, set the candidate view's current index to **proxy row 0** (a real,
   mapped row) instead of `clearSelection()/clearCurrentIndex()`. After this, `onCandidateSelectionChanged()`
   sees a non-null `record` and the mode boxes reflect that descriptor. **Do NOT auto-select an active/added
   row** — the active list (`activeView`) stays empty on a clean start (that's the lead's clean-unselected
   ask; this is the reconciliation — picker has a live default row, added list stays empty).
   - Refactor the "select row 0 if rows exist and none current" into one small private helper
     (e.g. `ensureCandidateRowSelected()`) and call it from both `refreshCatalog()` and
     `onAnchorSelectionChanged()` so they can't drift.
2. **Split the disabled tooltip** in `onCandidateSelectionChanged()` so the greyed reason is truthful per
   axis (the three predicates disagree — STEP 1's `/dashboard/state` already exposes them). For each disabled
   box, set the reason:
   - `modeId` empty (descriptor doesn't advertise this kind) → "This descriptor does not offer that display
     mode."
   - `modeId` present but `!DashboardModeHasVisibleSurface(modeId)` → "This display mode has no implemented
     renderer yet."
   - (If you can cleanly reach loaded-run availability for the descriptor here: data unavailable → "The data
     for this descriptor is not available in this run." If not cleanly reachable, skip availability for now —
     do NOT fabricate; leave a `// TODO(availability-step)`.)
   Use the now-global `DashboardModeHasVisibleSurface(...)` (STEP 1 lifted it out of the anonymous namespace
   in `SignalDisplayDialog.cpp`). Keep the enabled-box tooltip as-is.

No other behavior change. Do not touch the active-side gating logic, the panel/controller predicates, or the
two-model coordinator — those are later steps.

## PART B — the picker probe (so I verify headlessly)
The dialog is modeless and owned by `ReaderMainWindow` (`SignalDisplayDialog`). Add a small REST surface so a
headless caller can open the picker and read its live selector state. Verify the real open path first:
`ReaderMainWindow::onOpenSignalDisplays()` (the Metrics action) and its focus gate (the Metrics action is
enabled only when `AtomSelection` has focus — read `ReaderMainWindow.cpp` around the Metrics action wiring).
- **POST `/dashboard/picker/open`** `{"atom": int (optional)}` — if `atom` given, pick it first (reuse the
  existing pick path, e.g. the same call `/selection/pick` uses) so selection has focus and the gate opens;
  then open the dialog through the SAME path the Metrics action uses (show/raise + `refreshCatalog()`). Return
  the picker state (same shape as GET below). If the focus gate blocks opening, return that reason — don't
  force around it.
- **GET `/dashboard/picker`** — read the dialog's live widget state and return:
  ```json
  { "open": bool, "candidate_row": int|null, "anchor": {...},
    "modes": [ { "kind": "...", "mode_id": "...", "enabled": bool, "checked": bool, "reason": "..." } ] }
  ```
  Read the dialog's `candidateModes` boxes (`enabled()`, `isChecked()`, the `modeId` property, the tooltip/
  reason). Add a minimal `const` accessor on `ReaderMainWindow`/`SignalDisplayDialog` if needed; keep the
  dialog the owner of its state (RestServer just reads it).
- Also emit one `h5reader.dashboard` (`cDash`) log line at the tail of `onCandidateSelectionChanged()`:
  `event=picker_selector_state current_row=<r|none> enabled=[kind,...] disabled=[kind:reason,...]`
  so the fix is also provable from a log tail.

## Discipline
`CENSUS_REGISTER` in any new QObject ctor; `ACONNECT` (not raw `connect`); `ASSERT_THREAD(this)` on
thread-sensitive methods; no new `QTimer`; full error handling at the JSON boundary; portable; match existing
idioms. The probe is additive observability — no behavior change beyond PART A.

## Output
Implement, build green (`cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target
h5reader -j$(nproc)`). **Do NOT launch, do NOT git.** Append a short section to
`notes/INSTRUMENT_DASHBOARD_2026-06-05.md` (the new `/dashboard/picker*` routes + the selector-state log line,
with example `curl`/tail), and note anything in the code that contradicted this brief (code wins).
