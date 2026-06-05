# Codex — STEP 2.5b: dock reveal-on-add + run-load reset (the lead's resolved rules)

## FIRST — read the qt6-cpp skill before touching any code (full filesystem access — open them)
- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md` ← signal/slot wiring, app startup/shutdown, dock/window state

Honor it: `ACONNECT` (never raw `connect`), `ASSERT_THREAD(this)` on UI-thread methods, **no new `QTimer`**
(the existing `revealDockQueued` singleShot is fine — reuse it), portable Qt (Windows/macOS/Linux first-class).

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`; floor `8007e2e` +
uncommitted 2.2/2.3/2.5a). Scope `h5-reader/src` ONLY. **Verify from the code before changing anything.**

## Why this matters
2.5a made `DashboardSelectionController` the single owner of selection and gave it a `selectedCountChanged(int)`
signal + `clearAllMetrics()` — currently unconsumed *on purpose*. This step consumes them to implement the
two dock-lifecycle rules the lead decided, completing the comprehensive fix. These are HER exact calls —
implement them literally, do not add cleverness:

> **Dock reveal:** "on metric-add, if the dock is hidden, REVEAL it. None? Now one? You get a dock. Any
> add-while-no-dock shows the dock (not just the first metric). That single add event is the ONLY trigger —
> **no complex derived-visibility rule, no `userDisabled` suppression.**"

> **Run / .LGS load:** "all-new — a run load resets the selected-metrics list + dock state (clean slate). The
> selected list is per-run; reload is the reset point."

## Read first and confirm
- `DashboardSelectionController` (2.5a) — its `selectedCountChanged(int)` / `selectionChanged()` emission
  points and `clearAllMetrics()`.
- `ReaderMainWindow` — the dock members (`dashboardStripDock_`), `revealDockQueued(QDockWidget*)`, the
  controller instantiation, and the load/wiring sequence (where the protein/conformation/models get wired
  after a run loads).
- **Where the selected list gets populated on load today** — there is a metric that appears
  *non-deterministically* across launches (`h5:dssp8_transition`), and it is NOT a hard-coded `addSignal` in
  `ReaderMainWindow` (already checked — none there). Find the actual source (likely a `QSettings`
  `restoreState`/dashboard-state restore, or `.LGS`/session persistence). Confirm it before neutralizing it.

## Deliverable
1. **Dock reveal-on-add** (in `ReaderMainWindow`, app-level — the dock is the window's): connect the
   controller's add signal so that **when a metric is added and `dashboardStripDock_` is hidden, call
   `revealDockQueued(dashboardStripDock_)`**. Trigger on *add specifically* (a count increase), not on
   removes — either by tracking the previous count from `selectedCountChanged(int)` (reveal when
   `count > prevCount` and the dock is hidden), or by a dedicated `metricAdded` signal if you find that
   cleaner. **No `userDisabled` flag, no derived `count>0 && !disabled` rule** — the lead explicitly rejected
   that. Just: add happened + dock hidden → reveal. Log `event=dock_reveal_on_add count=N`.
2. **Run-load reset to clean slate**: on a run/.LGS load, the selected-metrics list starts **empty** and the
   dashboard strip dock starts **hidden**. Implement by neutralizing whatever currently restores/seeds the
   selection on load (preferred: do not restore dashboard *selection* from persisted state — keep restoring
   window geometry/dock layout if that's separate), or by calling `controller->clearAllMetrics()` immediately
   after the load+wire completes if the seed can't be cleanly suppressed at the source. The dock stays hidden
   until the first add (rule 1 then reveals it). Log `event=selection_reset_on_load`. Confirm the
   non-deterministic `h5:dssp8_transition` seed is gone — every launch opens with `selected: []`.
3. Keep `revealDockQueued`'s existing queued resize/raise behavior (the all-tabified/all-hidden Qt corner the
   prior session hit). Don't regress dock width (min 260 / resize-to-360).
4. **Make the REAL Add path testable headlessly** (this is critical — our prior verification drove
   `controller->addMetric` directly via REST and never exercised the dialog's actual Add button, so a bug in
   the dialog handler would be invisible to the harness). Add **`POST /dashboard/picker/add`** that invokes
   the dialog's genuine `SignalDisplayDialog::onAddSelected()` slot — the exact code the user's "Add Signal"
   button runs — using whatever is currently selected in the dialog (call `/dashboard/picker/open` first to
   set up the anchor/descriptor). Return the resulting picker state + selected count + the dock's
   visible/width. This lets the harness prove the **user-equivalent** flow: open picker → (row auto-selected,
   mode checked) → **click Add via this route** → metric added AND dock auto-revealed — with NO manual
   `/dashboard/dock` call. If `onAddSelected` can't be invoked cleanly (e.g. it reads widget state that REST
   can't populate), say so and tell me exactly what the Add button depends on.

## Preserve
- 2.5a behavior (add/remove/mode/close-panel through the controller) — untouched.
- `/dashboard/dock` REST + the toolbar Panels toggles still work (a user can still manually hide/show).
- Window geometry/dock-layout persistence (if separate from selection) is fine to keep — only the *selection*
  resets on load.

## Output
Implement, build green (`cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-rwdi --target
h5reader -j$(nproc)`). **Do NOT launch, do NOT git.** Append a `notes/INSTRUMENT_DASHBOARD_2026-06-05.md`
section (the reveal-on-add hook + the reset-on-load mechanism + where the old seed came from), and note
anything that contradicted this brief (code wins — tell me what you found about the seed source).
