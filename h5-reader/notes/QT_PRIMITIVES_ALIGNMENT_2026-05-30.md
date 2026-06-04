# Qt primitives alignment — chewer + viewer refactor scope

Date: 2026-05-30 (end-of-session capture)
Status: planning doc; captures the canonical-Qt alignment items
surfaced during the chewer design sessions PLUS the pending viewer
animation refactor scope. Companion to
`h5-reader/spec/reader_as_platform_2026-05-29.md`,
`h5-reader/spec/substrate_calculator_framework_2026-05-30.md`, and
the existing `h5-reader/notes/ROBUSTNESS_BACKLOG_2026-05-30.md`.

June 4 status: partly superseded. Some primitive/refactor intent landed or
shifted with the viewport/camera/dashboard work; `QSettings` persistence is now
present. Treat the old "ALL committed work" language below as design intent,
not current task state.

| Area | 2026-06-04 status |
|---|---|
| Viewport/camera refactor | partly landed; see `VIEWPORT_OBSERVATIONS_2026-05-30.md` and UI state |
| QSettings persistence | done |
| Chewer/future batch primitives | verify against current specs before acting |
| Session work-order text | historical only |

## Frame

The chewer redesign converged on a canonical-Qt pattern
(`QtConcurrent::run + QPromise + QFutureWatcher`) after noticing we
were designing custom plumbing for a problem Qt had already solved.
That recognition prompted a sweep of the rest of the project for
other "we built custom code for X" patterns where Qt has a canonical
answer. This doc captures the resulting inventory: 10 items, ordered
by cost.

It also captures a separate but related decision: **we will be
redoing the protein animation display** to fix two long-standing bugs
(the failure-to-paint-some-frames issue and the inability to remain
centred on an atom+bond pair properly during scrub/playback). That
work is in the same neighborhood as the Conformation /
DftShieldingStore frame-source providers, which means several
"medium-cost Qt-alignment" items become natural to land alongside the
viewer animation refactor rather than as separate passes.

## The 10 items (ALL committed work)

All ten are planned work, not a "do these now, maybe later for the
rest" list. Categorization below is purely sequencing — when each
lands in the implementation flow — not a priority filter. Per the
2026-05-30 session: "These are all worth pursuing. We are looking at
substantial changes anyway."

### Group A — land in the spec docs in this pass (mechanical)

**1. Provenance manifest format: YAML → JSON**

The substrate's provenance sidecar was specced as `manifest.yaml`.
YAML needs an external dependency (yaml-cpp). The existing
`StructuredLogger` already emits JSON; `QJsonDocument` is Qt's
canonical JSON primitive. Switching the manifest to `manifest.json`
eliminates a dependency, matches the existing UDP log format,
consumers ingest the manifest with the same JSON tooling.

**Where to land**: `substrate_conventions_2026-05-30.md` "Provenance
metadata sidecar" section + the file naming example.

**2. Chewer run IDs: `QUuid::createUuid()`**

The REST API plans `POST /chewer/runs` returning `{"run_id": "uuid"}`.
Qt has `QUuid::createUuid()` — canonical UUID generation.
Specification: each chewer run gets a `QUuid` as its run_id; REST
serializes via `QUuid::toString(QUuid::WithoutBraces)`.

**Where to land**: `reader_as_platform_2026-05-29.md` REST API
section.

**3. Protein load as `QFuture<LoadedProtein>`**

The chewer's load phase should explicitly use `QFuture<LoadedProtein>`
returned from `QtConcurrent::run(loaderPool, loadProtein, path)`.
Then `loaded.waitForFinished()` inside `RunChewer` (the chewer
already runs on a worker thread; blocking on the load future is
fine), OR streamed via `QPromise<LoadEvent>` for per-frame progress
during the load. Matches the chewer's own QtConcurrent pattern; one
async-load idiom across the codebase.

**Where to land**: `substrate_calculator_framework_2026-05-30.md`
loader-calculator section; `reader_as_platform_2026-05-29.md` load
section.

**4. `QCommandLineParser` for `--batch` CLI args**

`main_reader.cpp` today has positional + `--rest` handling, probably
ad-hoc. Adding `--batch`, `--script`, `--protein`, `--output-dir`,
`--log-file` etc. is canonical-Qt via `QCommandLineParser`. Provides
`--help` and `--version` for free, proper argument validation,
consistent error messages.

**Verify first**: check whether main_reader.cpp already uses
`QCommandLineParser`. If yes, just add the new options. If no, the
switch is a small cleanup worth doing alongside adding the chewer
flags.

**Where to land**: `reader_as_platform_2026-05-29.md` CLI mode
section.

### Group B — paired with the viewer animation refactor (committed work)

These two items are committed; the viewer animation refactor (see
"Viewer animation refactor scope" below) is the natural moment to
land them because both touch the frame-source residency machinery
that the animation refactor is already opening up. Not optional, not
"defer indefinitely" — they ride the same touch.

**5. `Conformation::requestSnapshot` / `snapshotReady` → `QFuture<QtConformationSnapshot>`**

The viewer's current async-snapshot-load API is custom (`snapshot()`
returns null-or-resident; `requestSnapshot()` triggers load;
`snapshotReady(frame)` signal fires when done). Canonical Qt is
`QFuture<shared_ptr<const QtConformationSnapshot>>` returned from a
`loadSnapshotAsync(frame)` method. Consumers attach
`QFutureWatcher` and react to `finished`.

**Wins**:
- Explicit lifetime: the QFuture is the handle; watcher destruction
  cleans up cleanly
- Composes with cancellation: `future.cancel()` aborts a pending
  snapshot load that's already been superseded (solves the "user
  scrubs past a frame before its load completes" wasted-work problem)
- One async-load idiom across the codebase: same as the chewer

**Cost**: meaningful viewer refactor. Multiple subscribers (inspector
dock, dashboard controller) currently react to `snapshotReady`;
they'd switch to QFutureWatcher.

**Pair with viewer animation refactor.**

**6. `DftShieldingStore::requestFrame` / `frameReady` → `QFuture<DftShieldingFrame>`**

Identical pattern to item 5, same canonical replacement, same
refactor cost. Pair with the viewer animation refactor for the same
reason.

### Group C — chewer-foundation work (committed; lands with the chewer)

**7. `ErrorBus` — audit + confirm or replace**

The reader has `src/diagnostics/ErrorBus.cpp` for cross-thread error
reporting. The chewer will emit errors via `ChewerEvent::Error` AND
through the existing diagnostics layer; we need to know whether
ErrorBus is the canonical wrapper over Qt's standard error-handling
(`qWarning`/`qCritical` + `qInstallMessageHandler` + cross-thread
signal delivery) or whether it's custom plumbing reinventing what Qt
already does.

Action: 15-minute audit before chewer code lands. If thin wrapper,
confirm in `substrate_calculator_framework_2026-05-30.md` that chewer
errors route through it. If custom plumbing, either cleanup-in-place
or document the divergence so chewer doesn't accidentally double up
on its own error channel.

**8. `QFileSystemWatcher` for transform-script live-reload**

Wins development workflow: edit a transform script, save, the
chewer's next iteration picks up the change without restart.
`QFileSystemWatcher` on the configured transform-script directory
fires `fileChanged(QString)` → chewer marks the script's module for
reload via `py::module_::reload`.

Committed scope; lands as part of the chewer's transform-loader work
(not deferred indefinitely). For the dogfood case this is "load
once, run, exit"; for the iterate-on-PySR-features case it's a
quality-of-life win that costs ~30 lines of code.

### Group D — already in adjacent specs (committed; not re-scoped here)

**9. `QSettings` for any persistent state**

Already in `ROBUSTNESS_BACKLOG_2026-05-30.md` item 6 with the
`QMainWindow::saveState/restoreState` + version-tagged-blob pattern
as the canonical approach. When settings persistence lands, it lands
that way. No re-spec needed here.

**10. `QStandardPaths` for default output locations**

Already in `substrate_conventions_2026-05-30.md` as the canonical
primitive for any standard-location lookup. When the chewer needs a
default output directory, it uses
`QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation)`.
No re-spec needed here.

## Viewer animation refactor scope (new — not yet specced)

Separate from the chewer work but in the same architectural
neighborhood. Two long-standing bugs to fix:

- **Failure to paint some frames during scrub or playback.** Already
  noted in code comments and the recent robustness audit (item 6,
  `MoleculeScene::setFrame` diagnostic probe around line 501 — the
  `mapper_->Modified()` workaround for intermittent atom-render
  drop). The probe is "for now" code that hasn't been resolved into
  a settled fix. The animation refactor is the natural moment to
  reproduce the bug deliberately, identify the actual cause, fix it
  cleanly, and remove the probe.

- **Inability to remain centred on an atom+bond properly during
  scrub/playback.** The user picks an atom (or atom+bond pair) and
  wants the camera to track it as frames advance. Currently
  insufficient or absent. Requires a per-frame camera-target update
  hook in `MoleculeScene` driven by the atom selection + the bond
  vector if the user picked a bond. Likely needs interpolation
  between frames (smooth tracking, not jerk-per-frame) and
  consideration of orientation (track-the-bond-axis vs track-just-
  the-atom).

**The pairing with items 5 + 6**: both bugs live in `MoleculeScene`
and the frame-advance signal chain. The async-snapshot-load refactor
(item 5) is in the same chain. Doing items 5 + 6 alongside the
animation refactor means we only touch those subscriber chains once.

**Worth a dedicated design doc** when we get to it
(`spec/viewer_animation_refactor_2026-XX-XX.md` likely). Out of scope
for this session.

## Suggested work order

Dependency-flow sequencing, not priority filtering. All ten items
are committed.

1. **Group A landed in this pass** — items 1, 2, 3, 4 (manifest
   JSON, QUuid run IDs, QFuture<LoadedProtein>, QCommandLineParser).
   Mechanical doc edits already applied to the spec layer.

2. **Item 7 (ErrorBus audit)** — 15 minutes before chewer code
   lands. Outcome: confirm-canonical or fix-as-part-of-chewer-foundation.

3. **Chewer implementation starts** — Group A conventions in place;
   item 8 (`QFileSystemWatcher` for script live-reload) lands as part
   of the chewer's transform-loader work; items 5 + 6 stay paired
   with the viewer animation refactor (the chewer's worker thread
   doesn't touch the viewer's snapshot/DFT residency, so the
   pairing is clean).

4. **Viewer animation refactor** — when prioritized: fixes the
   paint-some-frames bug + atom+bond centring, AND lands items 5 + 6
   (Conformation → `QFuture<QtConformationSnapshot>`, DftShieldingStore
   → `QFuture<DftShieldingFrame>`) as part of the same architectural
   touch. Item 6 in `ROBUSTNESS_BACKLOG_2026-05-30.md` (the
   MoleculeScene rendering probe) gets resolved as part of this work
   too — same neighborhood.

5. **Item 9 (`QSettings`) lands when settings persistence is added**
   per the robustness backlog's existing schedule (judgment-call
   item there; the canonical pattern is already named).

6. **Item 10 (`QStandardPaths`) lands the first time the chewer
   needs a default output location** — pattern already documented in
   conventions doc.

## What this doesn't change

- The chewer's QtConcurrent + QPromise + QFutureWatcher design from
  the platform + framework docs.
- The substrate's eager-everything memory model.
- The DRY discipline (reuse h5-reader's existing typed structures).
- The conventions doc (modulo item 1: manifest format).
- The reader-untouchable analog (chewer additions in `src/analysis/`).
- The foundational refactor (commit `518855e`) — already landed
  Changes A + B.

## Cross-references

- `h5-reader/spec/INDEX.md` — entry point for the chewer spec layer
- `h5-reader/spec/reader_as_platform_2026-05-29.md` — platform
  architecture; items 2, 3, 4 land here
- `h5-reader/spec/substrate_calculator_framework_2026-05-30.md` —
  calculator framework; item 3 lands here
- `h5-reader/spec/substrate_conventions_2026-05-30.md` — conventions
  doc; item 1 lands here
- `h5-reader/notes/ROBUSTNESS_BACKLOG_2026-05-30.md` — robustness
  fix list; item 6 there (MoleculeScene rendering probe) pairs with
  the viewer animation refactor
- `qt6-cpp` skill at `/home/jessica/.claude/skills/qt6-cpp/` —
  canonical patterns for every item above; `references/architecture.md`
  + `references/3d-vtk.md` particularly relevant for the viewer
  animation refactor
