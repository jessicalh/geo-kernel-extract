# Codex brief — named two-mode stabilisation switch (2026-06-06)

> **PREAMBLE (prepended verbatim from `notes/CODEX_PREAMBLE_QT_VTK.md`):**

## FIRST — load the qt6-cpp skill (before any project code), and you have VTK-docs license
**Before touching any project code, READ these in full** (you have filesystem access — open them):
- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/3d-vtk.md`  ← VTK 9.5 + Qt: camera, pipeline, interactor
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md`  ← signal/slot, ownership, startup/shutdown
- `/home/jessica/.claude/skills/qt6-cpp/references/model-view.md`  ← custom models, roles, delegates

This is a **portable Qt6/VTK** reader (Windows / macOS / Linux all first-class — treat the cross-platform source
as real). **The skill's SKILL.md is written Windows-first (MSVC/`C:\` paths, `Enter-VsDevShell`, minidumps) — that
framing is fine to read past; the DISCIPLINE is what transfers and is mandatory here on Linux too.** Honour the
idioms: `CENSUS_REGISTER(this)` in every QObject ctor, `ACONNECT` (never raw `connect`), `ASSERT_THREAD(this)` on
thread-sensitive methods, **one** app-owned `QTimer` max, full error handling at every external-library seam,
UDP-log-first when debugging.

**You have EXPLICIT LICENSE to consult VTK documentation ONLINE.** For any camera / renderer / interactor-style /
widget question, look up the authoritative VTK class references (vtk.org/doc/nightly) and cite the doc URL.

**Verify from the code before changing anything, ever** — cite `file:line`. Root
`/shared/2026Thesis/nmr-shielding/h5-reader`, scope `h5-reader/src` ONLY; never read/link the `nmr_shielding`
library, never write H5, never trigger extraction.

---

## The task — a single, self-naming stabilisation-mode button ("seeing the button")

Hi codex — small, contained, high-taste UI change. The reader's toolbar has a single checkable `QAction` labelled
**"All-atom fit"** that toggles the display *stabilisation* between two modes. Two problems the lead wants fixed:
(1) the button names the **target action**, not the **mode you're currently in**, so you can't tell which mode is
active; (2) "All-atom fit" is jargon. The fix: **one toolbar button whose TEXT names the *current* stabilisation
mode**, with a swap affordance and a tooltip explaining both. This is a **presentation change over already-correct
machinery** — do NOT change the transform math, the `TransformedConformation` model, or the mode semantics.

### The two modes (already wired — keep semantics exactly)
- `TransformedConformation::Mode::FitSubset` (backbone subset) = **"Locked backbone"** — the startup default.
- `TransformedConformation::Mode::FitReference` (all atoms) = **"Kabsch with give"**.
- enum: `src/model/TransformedConformation.h:64-67`.

### Verify these anchors before editing (cite `file:line` in your writeup)
- Action creation: `src/app/ReaderMainWindow.cpp:914-920` — `transformFitAction_`, currently checkable, label
  "All-atom fit", `toggled` → `onTransformFitToggled`.
- The `transformChanged` sync lambda: `src/app/ReaderMainWindow.cpp:139-147` — currently
  `transformFitAction_->setChecked(mode==FitReference)` on every `transformChanged`, plus
  `scene_->refreshCurrentFrame()`.
- Startup mode set + fallback: `src/app/ReaderMainWindow.cpp:124-138`.
- The toggle handler: `src/app/ReaderMainWindow.cpp:1162-1183` — `onTransformFitToggled(bool)`.
- `setMode` emits `transformChanged`: `src/model/TransformedConformation.cpp:193,212`; signal decl `.h:117`.
  **REST also changes mode** via `POST /transform` (`src/app/RestServer.cpp:~808`, `{"kind":"all_atom_fit"|"backbone_fit"}`)
  and that path ALSO emits `transformChanged`. **Therefore the button label MUST be driven from the
  `transformChanged` lambda, not only from the click handler** — else a REST-driven switch leaves the label stale.

### The change — `ReaderMainWindow.{h,cpp}` ONLY
1. **Make `transformFitAction_` non-checkable**: remove `setCheckable(true)` / `setChecked(false)` (`:915-916`).
   Connect `QAction::triggered` (not `toggled`) via `ACONNECT` to a new slot.
2. **New click slot** (replace `onTransformFitToggled(bool)` with e.g. `onTransformFitClicked()`):
   `ASSERT_THREAD(this)`; guard `!transformed_||!loaded_||!loaded_->protein`. Read `transformed_->mode()`:
   - if `FitReference` → switch to `FitSubset` using `TransformedConformation::BackboneSubset(*loaded_->protein)`,
     **keeping the existing guard**: if `subset.size() < 3`, `qCWarning(cWindow)` and stay `FitReference`
     (`setMode(FitReference,0)`), return.
   - else (`FitSubset`) → `setMode(FitReference, 0)`.
   No `setChecked` calls anywhere.
3. **Add a private helper** `void updateFitModeLabel()`: reads `transformed_->mode()` and sets the action text +
   tooltip:
   - `FitSubset`   → text `"Mode: Locked backbone  ⇄"`
   - `FitReference` → text `"Mode: Kabsch with give  ⇄"`
   - tooltip (names both modes):
     `"Stabilisation mode — click to switch.\nLocked backbone: Kabsch fit of the backbone (industry standard) — removes global tumbling; the backbone holds still while sidechains/loops move.\nKabsch with give: all-atom fit — removes tumbling but lets real internal motion show."`
   - Glyph: use `⇄` (U+21C4). Only if you have a concrete reason it won't render in the default Qt toolbar font on
     a target platform, use `↔` (U+2194) and say why.
4. **Drive the label from `transformChanged`**: in the `:139-147` lambda, replace the `setChecked(...)` block with a
   call to `updateFitModeLabel()`. Keep `scene_->refreshCurrentFrame()`. (No `QSignalBlocker` needed — the action
   is no longer checkable and `setText` doesn't re-enter the handler.)
5. **Initial label sync**: startup `setMode` (`:129`/`:133`) emits `transformChanged` BEFORE this lambda is
   connected, so call `updateFitModeLabel()` ONCE explicitly right after the `ACONNECT` (after `:147`) to set the
   opening label. In the startup fallback block (`:134-137`), drop the `setChecked(true)` (label now handled by the
   helper); keep the `qCWarning` + `setMode(FitReference,0)`.
6. **Header** (`ReaderMainWindow.h`): rename the slot decl (`onTransformFitToggled` → `onTransformFitClicked`,
   no bool param) and add the `updateFitModeLabel()` private method decl. Update any other references.
7. **Do NOT touch the toolbar's global `toolButtonStyle`** — the action already renders its text; other toolbar
   buttons must be unaffected.

### Discipline
`ACONNECT` (never raw `connect`) for the `triggered` connection; `ASSERT_THREAD(this)` at the top of the new slot
(the old handler had it); keep `CENSUS_REGISTER`/`Q_LOGGING_CATEGORY` patterns intact. No `QTimer`. No model / H5 /
library changes. T2 / kernels untouched (not in scope here).

### Build + verify
- From `h5-reader/`: build the existing preset incrementally — `cmake --build build/linux-rwdi`
  (if not configured, `cmake --preset linux-rwdi` first). 
- Tests: `ctest --test-dir build/linux-rwdi --output-on-failure` — the model/registry tests must stay green.
- Confirm a clean build. You do NOT need to launch the GUI — opus verifies visually on a headless `:99` instance.

### Writeup back to me
The `file:line` anchors you verified, the diff summary, build result (green or the errors), ctest result, and any
judgment calls (especially the glyph). A positive account of what you did. **If anything in this brief contradicts
the code you find, STOP and report — do not guess.**
