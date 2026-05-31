# Reader-as-platform architecture (design intent)

Date: 2026-05-29
Status: settled architecture across the 2026-05-29 conversation.
Companion to `stage2_pysr_campaign_2026-05-29.md` — this doc is the
**platform architecture**; that doc is the Stage 2 campaign which uses
the platform. Both written so an agent can read either and have full
context.

Headline lens added 2026-05-31 after the lens emerged in conversation
and was not load-bearing in the original draft; see memory entry
`feedback_model_is_spine`.

## What this is, in one paragraph

The chewer is **scripting against the live C++ protein model**. The
reader already hosts that model — typed atoms, residues, rings, bonds,
kernels — in memory, populated from H5 via existing loaders. The
chewer's pybind11 surface exposes that model to Python; transforms are
short scripts that query it in its native NMR-chemistry vocabulary
(`atom.is_aromatic_HE()`, `ring.member_atoms`,
`frame.kernel('BS', atom).T0`, `residue.phi(t)`) and emit derived rows.
Output format (Parquet, Postgres, Arrow IPC, static snapshots) is
downstream of the model — chosen for distribution and consumption, not
for query semantics. Spatial structure lives in the running model; the
durable output is whatever subset of "answers from the model" needs to
be persisted for downstream tools or methods-paper reviewers.

Everything below — the QtConcurrent execution model, the eager-loaded
substrate features list, the calculator framework, the CLI/REST/GUI
triggers — is implementation of this lens. Read it through this frame.

## Goal

Move a lot of data through the h5-reader app for sophisticated analysis.

**Two applications to maintain: extract (producer) and reader (consumer
+ analysis platform).** The reader binary handles two modes through
the same machinery:

- **Interactive (GUI)**: today's behavior — load a protein, scrub
  frames, inspect signals, do thesis-vetting work
- **Batch (CLI)**: new capability — load a protein, run registered
  Python transform scripts in a frame-by-frame walk, emit per-transform
  Parquet outputs

Both modes use the **same C++ chewer**, triggered by either a toolbar
button (GUI) or a CLI flag (`--batch`).

## Benchmark

This is **not VLC**. The workload is compute-heavy, not bandwidth-heavy.
Data fits in RAM after load (~270 MB for 1P9J). Bytes flow once at
load; multiple consumers then read from the same immutable buffer at
their natural pace. Indexing (KD-trees, ring centroids) and computation
(kernel queries, finite differences, transform features) dominate cost
— and **they're allowed to lag** the UI.

The right reference is scientific viewers with heavy compute over
shared immutable data: **ParaView, MeshLab, Avogadro2, Chimera**.
Not real-time streaming media players.

This matters because it sets the simplicity bar: we do NOT need
producer-consumer pipelines with inter-stage backpressure, lock-free
MPMC queues, frame budgets, or per-stage thread pools. We need shared
immutable data + multiple independent readers + Qt's standard
threading and signaling primitives.

> **Simple but not too simple.** Don't add architecture the workload
> doesn't require. Don't strip out machinery a long-running Qt process
> legitimately needs.

---

## Architectural primitives

### Shared immutable data

```cpp
struct LoadedProtein {
    QtProtein protein;
    QtTopology topology;
    // typed buffers from QtTrajectoryH5: positions, kernel TRs, etc.
    // immutable after construction
};

std::shared_ptr<const LoadedProtein> loaded;
```

- `LoadedProtein` is read-only after construction
- Every consumer holds `std::shared_ptr<const LoadedProtein>`
- **No locks anywhere**, because nothing mutates after load
- UI consumer reads positions + topology + per-frame kernel values for
  rendering
- Chewer consumer reads the same buffers for transform iteration
- Independent reads at independent paces; no coordination

### One chewer, three triggers

There is one `RunChewer` function in C++ (see "Execution model" below
for the full QtConcurrent + QPromise + QFutureWatcher pattern). Three
ways to fire it: CLI batch, REST endpoint, GUI toolbar button. The
CLI is the natural entry point because it's simplest (pure
`QCoreApplication`, no widgets); REST and GUI reuse the same function
and the same `QFuture`-based lifecycle.

**The deliverable is dogfood, not slicing.** The first goal is "we ran
a real Stage 2 transform against 1P9J and got real PySR equations
out, end-to-end" — not "pure chewer machinery shipped as a milestone."
CLI batch may be the natural first entry point, but if the toolbar
button is trivial to add at the same time, add it. If REST is easy to
wire because the server is already there, wire it. Don't gate one on
the others. Real use surfaces real issues; ceremonial milestones
don't.

**GUI mode**: toolbar button. User has a protein loaded, selects a set
of transform scripts from a dialog, clicks Run. The button handler
calls `QtConcurrent::run(chewerPool_, RunChewer, config)` and binds a
`QFutureWatcher<ChewerEvent>` to the returned `QFuture`; the watcher
delivers progress + events to the GUI via queued signals.

**CLI mode (the first mode to build)**: `--batch` flag, **pure
headless**. No `QApplication`, no widgets, no VTK initialization. Just
`QCoreApplication` + a dedicated `QThreadPool` + the chewer function.

```bash
h5reader --batch \
  --script transforms/bs_aromatic_he.py \
  --protein /path/to/1p9j.h5 \
  --output-dir runs/1p9j_2026-06-XX/ \
  [--log-file path/to/log.jsonl]    # optional; default is UDP-only
```

CLI parsing uses `QCommandLineParser` (canonical Qt) — provides
`--help` and `--version` for free, proper argument validation,
consistent error messages. The dispatcher main parses the mode
(`--batch` vs interactive) BEFORE constructing `QApplication`/VTK, so
batch mode never initializes GUI state.

The CLI entry point (`main_batch.cpp`) constructs `QCoreApplication`,
sets up a `QFutureWatcher<ChewerEvent>` for progress logging
(UDP primary; stderr only for phase-boundary events, NOT the per-frame
hot-path emits — see the framework doc's "Logging — UDP primary"
section), calls `QtConcurrent::run(pool, RunChewer, config)`,
blocks on `QCoreApplication::exec()` until the watcher's `finished`
signal arrives (or until SIGINT triggers `chewerFuture.cancel()`
through the existing `ShutdownSignals` infrastructure). Progress to
UDP log primary; stderr for phase-boundary events only (zero-write
hot path per the framework doc's logging spec); no progress dialog
(none would render — no GUI).

This mode is the smallest self-contained slice that proves the whole
architecture: substrate + QtConcurrent worker + pybind11 + Python
execution + 3-level cancel. Develop and test this first, in isolation
from any GUI coordination concerns. Once it works end-to-end against
the 1P9J fixture, the same `RunChewer` function is reused for REST
and GUI without architectural changes.

**REST mode**: `POST /chewer/runs` to a running `h5reader` instance
(GUI or headless). Returns a `run_id`; subsequent `GET /chewer/runs/{id}`,
`POST /chewer/runs/{id}/cancel`, `GET /chewer/runs/{id}/outputs`
drive and observe the run remotely. Same `RunChewer` function under
the hood; the REST handler constructs the config from the request
body and launches via `QtConcurrent::run(pool, RunChewer, config)`,
storing the returned `QFuture` indexed by run_id for subsequent
cancel / status calls.

Same function, same QFuture-based lifecycle, same event stream. Three
ways to fire it.

### REST API (new endpoints on existing `RestServer`)

The current `RestServer` is described in the h5-reader notes as "a
test surface, not a public analysis API." The Stage 2 platform
promotes it to a real analysis API by adding chewer-control endpoints:

```
POST /chewer/runs                              # start a run
  body: {"script": "transforms/bs_aromatic_he.py",
         "protein": "/path/to/1p9j.h5",
         "output_dir": "runs/foo/"}
  returns: {"run_id": "<QUuid>", "state": "starting"}

# run_id is QUuid::createUuid().toString(QUuid::WithoutBraces) — canonical
# Qt UUID generation; consumers treat as opaque string. The chewer's
# QFuture and QFutureWatcher are stored in a map keyed by run_id for
# subsequent /chewer/runs/{id}/* dispatch.

GET  /chewer/runs                              # list runs
GET  /chewer/runs/{id}                         # state + progress
GET  /chewer/runs/{id}/log                     # recent UDP log lines for this run
POST /chewer/runs/{id}/cancel                  # triggers 3-level cancel
GET  /chewer/runs/{id}/outputs                 # paths + row counts
GET  /chewer/runs/{id}/transform/{name}/state  # Python-registered diagnostic endpoint
```

Python-side endpoint registration is exposed via the substrate:

```python
def main(traj, output_dir):
    traj.register_endpoint("bs_aromatic_he/sample_count",
                           lambda: {"rows": len(self.rows)})
    # ... iteration emits rows; the registered endpoint returns the
    # current count whenever queried over REST ...
```

The substrate's `register_endpoint(path, callable)` forwards to the
running `QtHttpServer` (which the REST server already wraps). Python
doesn't see Qt directly; just registers a callable that returns
JSON-serializable data.

### Testing strategy (Python heavy)

The platform refactor entails writing substantial Python:

1. **Transform scripts** during the walk — emit rows, poll
   `traj.cancel_requested()`, optionally register diagnostic endpoints
2. **SDK scripts** post-walk — consume Parquets, fit models (PySR /
   sklearn / MACE), write results
3. **REST integration tests** (pytest with `requests`/`httpx`) — drive
   end-to-end runs against a launched `h5reader` process: start a
   run, poll until done, fetch outputs, assert shapes

Testing infrastructure:

- Existing `tests/rest/` pytest scaffolding extends to cover new
  chewer endpoints
- New `tests/transforms/` for unit-testing transforms against fixture
  trajectories (no REST involved; instantiate the substrate directly
  via Python bindings)
- New `tests/sdk/` for unit-testing SDK scripts against fixture
  Parquets (no h5reader involved; pure Python + Arrow)
- Integration tests that exercise the full pipeline via REST: launch
  h5reader, POST a chewer run, poll, fetch outputs, run an SDK script
  on the outputs, assert results
- CI: the existing `xvfb-run` infrastructure handles the launched
  h5reader for REST tests; the transform/SDK unit tests are
  display-free

### Execution model — QtConcurrent + QPromise + QFutureWatcher

The chewer uses Qt's canonical concurrency primitives. NO custom
`QObject + QThread` subclass, NO custom cancel atomic, NO bespoke
signal/slot wiring. Reframed 2026-05-30 after noticing we were
designing custom plumbing for a problem Qt already solved.

The chewer is a function executed via `QtConcurrent::run` on a
dedicated `QThreadPool`. It receives a `QPromise<ChewerEvent>&` that
carries cancel state + a result stream back to the owner. The owner
holds a `QFutureWatcher<ChewerEvent>` that delivers events to the GUI
thread automatically.

```cpp
namespace h5reader::analysis {

// Typed event stream from chewer worker back to owner.
struct ChewerEvent {
    enum Kind {
        ProteinLoaded, CalculatorStarted, CalculatorFinished,
        TransformStarted, TransformFinished, Error
    };
    Kind     kind;
    QString  name;           // calculator name, transform name, protein id
    qint64   elapsed_ms = 0;
    QString  detail;
};

struct ChewerConfig {
    QStringList               proteinPaths;
    QList<TransformSpec>      transforms;
    QString                   outputDir;
    QPointer<QThreadPool>     parallelPool;   // for per-calculator parallelism
};

// The chewer's main function. Runs on a thread from a dedicated
// QThreadPool. promise carries cancel state via promise.isCanceled()
// and streams typed events via promise.addResult(ChewerEvent{...}).
// Numeric progress via promise.setProgressValue(framesDone).
void RunChewer(QPromise<ChewerEvent>& promise, ChewerConfig config);

}  // namespace h5reader::analysis
```

Owner side (GUI or batch main):

```cpp
// Dedicated pool — NOT QThreadPool::globalInstance() (which serves
// all of Qt's internal users; a long-running chewer would block its
// slots).
chewerPool_ = new QThreadPool(this);
chewerPool_->setMaxThreadCount(1);  // single long-running run

// Launch. QFuture is the handle for cancel + waitForFinished.
chewerFuture_ = QtConcurrent::run(chewerPool_, RunChewer, config);

// Watcher delivers events + numeric progress to GUI thread.
chewerWatcher_ = std::make_unique<QFutureWatcher<ChewerEvent>>(this);
CENSUS_REGISTER(chewerWatcher_.get());
ACONNECT(chewerWatcher_.get(), &QFutureWatcher<ChewerEvent>::resultReadyAt,
         this, &ReaderMainWindow::onChewerEvent);
ACONNECT(chewerWatcher_.get(), &QFutureWatcher<ChewerEvent>::finished,
         this, &ReaderMainWindow::onChewerFinished);
ACONNECT(chewerWatcher_.get(), &QFutureWatcher<ChewerEvent>::progressValueChanged,
         progressDialog_, &QProgressDialog::setValue);
chewerWatcher_->setFuture(chewerFuture_);
```

**What this gets us for free** (each line replaces custom plumbing
we'd otherwise have written):

| Custom plumbing we'd have built | What `QtConcurrent + QPromise + QFutureWatcher` gives us |
|---|---|
| `Chewer : public QObject` subclass | Free function `RunChewer` — no class |
| `QThread` subclass + lifecycle | `QtConcurrent::run` handles thread dispatch |
| Custom `cancelFlag_` atomic + `requestCancel()` direct method | `QFuture::cancel()` + `QPromise::isCanceled()` polling |
| Custom `progress(int,int)` signal | `QPromise::setProgressValue` → `QFutureWatcher::progressValueChanged` |
| Custom per-event signals (`started`, `proteinFinished`, etc.) | `QPromise::addResult(ChewerEvent)` → `QFutureWatcher::resultReadyAt` — one discriminated-union event type carries them all |
| Extending ACONNECT to accept `Qt::QueuedConnection` | QFutureWatcher's signals are documented as emitted on the thread that called `setFuture` (the GUI thread); ACONNECT with AutoConnection just works correctly |
| Custom shutdown/wait wiring for `QThread::quit() + wait()` | `chewerFuture_.cancel(); chewerFuture_.waitForFinished(timeoutMs);` |

**Why this is the right model**: every Qt application doing
"long-running work with progress + cancel + result" uses this exact
pattern. We aren't extending Qt; we're using its standard concurrency
primitives the way the qt6-cpp skill documents. Anyone reading the
chewer code recognizes the QFuture/QFutureWatcher pattern
immediately.

### Cancellation — three layers via QPromise + Python interrupt

The 3-level escalation stays the same; only the cooperative layer
moves from a custom atomic to QPromise's built-in cancel:

**Level 1 — cooperative (immediate)**: GUI calls `chewerFuture_.cancel()`.
`QFuture::cancel()` sets the promise's cancel state atomically. The
chewer worker polls `promise.isCanceled()` at frame boundaries; on
detection, finishes current frame, commits output sinks via
`QSaveFile::commit()` so partial outputs are atomic, returns from
`RunChewer`.

**Level 2 — escalate to KeyboardInterrupt (after ~5 sec if Level 1 didn't take)**:
GUI thread (via `QTimer::singleShot` from the cancel handler) calls
`PyErr_SetInterruptEx(SIGINT)`. Thread-safe per Python C API, no GIL
needed. Raises `KeyboardInterrupt` at the next opcode boundary in the
running Python interpreter. The transform script's try/except
KeyboardInterrupt handler does cleanup, returns. C++ side sees the
exception via pybind11; chewer continues to graceful shutdown path.

**Level 3 — log + abandon (after another ~10 sec)**: emit error,
`chewerFuture_.waitForFinished(timeoutMs)` may return false; log
warning; never `QThread::terminate()` (documented as dangerous).
Python interpreter cleaned up at process exit.

The GUI thread orchestrates the escalation directly (it has the
`chewerFuture_` handle); no need for queued slots to deliver cancel
to the worker.

### Tight inner loop (no Qt machinery per item)

Inside `RunChewer`, the per-(atom, frame) work is plain C++. No
signal emits per atom, no event posts per frame. The promise's
result-stream and progress-value are coarse-grained:

```cpp
void RunChewer(QPromise<ChewerEvent>& promise, ChewerConfig config) {
    py::scoped_interpreter interp;  // one per chewer run; lives on the worker thread

    for (auto& proteinPath : config.proteinPaths) {
        if (promise.isCanceled()) return;

        auto loaded = loadProtein(proteinPath);
        promise.addResult({ChewerEvent::ProteinLoaded, proteinPath, 0});
        promise.setProgressRange(0, loaded->frameCount());

        // Substrate calculators run in dependency order (see calculator framework doc).
        for (auto& calc : calculators) {
            if (promise.isCanceled()) return;
            promise.addResult({ChewerEvent::CalculatorStarted, calc.name(), 0});
            const auto t0 = std::chrono::steady_clock::now();
            calc.onLoad(*loaded);
            calc.onCompute(*loaded, promise);   // promise threaded through for cancel + progress
            const auto elapsed = elapsedMsSince(t0);
            promise.addResult({ChewerEvent::CalculatorFinished, calc.name(), elapsed});
        }

        // Run Python transforms (each pulls GIL via py::gil_scoped_acquire when calling into Python)
        for (auto& tspec : config.transforms) {
            if (promise.isCanceled()) return;
            promise.addResult({ChewerEvent::TransformStarted, tspec.name, 0});
            const auto t0 = std::chrono::steady_clock::now();
            runOnePythonTransform(*loaded, tspec, promise);
            promise.addResult({ChewerEvent::TransformFinished, tspec.name, elapsedMsSince(t0)});
        }
    }
}
```

Calculators check `promise.isCanceled()` at frame boundaries and call
`promise.setProgressValue(framesDone)` for numeric progress. Their
`onCompute` signature takes `QPromise<ChewerEvent>&` directly so they
can both check cancel and emit per-calculator events if desired (see
the calculator framework doc for the full base class shape).

The inner loop is straight iteration; the only Qt machinery per item
is the atomic-read in `promise.isCanceled()` (cheap, lock-free).
`QPromise::setProgressValue` and `QPromise::addResult` are coarse —
called per ~100 frames and per calculator/protein boundary
respectively, never per atom.

### Coarse coordination via QPromise results

Events flow from the chewer worker to the GUI thread via the typed
`ChewerEvent` discriminated union, streamed via
`promise.addResult(ChewerEvent{...})`. QFutureWatcher emits
`resultReadyAt(int index)` on the GUI thread for each new event; the
GUI handler indexes into the watcher's future to retrieve the
ChewerEvent and dispatch on its `kind`.

Numeric progress for the progress dialog uses the parallel channel
`promise.setProgressValue(framesDone)` → QFutureWatcher's
`progressValueChanged(int)`. Built-in to QFutureWatcher; no custom
signal needed.

Event budget per 30-minute walk: ~hundreds of `ChewerEvent` entries
+ a few thousand progress value updates. Both well within
QFutureWatcher's delivery capacity; both delivered to the GUI's
event loop without blocking the worker.

### File access: pragmatic mix

Use Qt's APIs **where natural for our code**:

- `QFile` for reading our own files
- `QSaveFile` for atomic writes (Parquet output, JSON results, status
  logs) — write to temp, rename on `commit()`, never leave half-files
  on crash
- `QDir` / `QFileInfo` for path operations
- `QString` for paths throughout, converted to `std::string` only at
  third-party boundaries
- `QStandardPaths` for any standard-location lookups

**Third-party libs (HighFive, Arrow, pybind11) using their own I/O is
fine, not cruft to convert.** They take `std::string` paths and do
their own internal POSIX file handling. The pragmatic rule: where we
do our own file ops, use Qt; where a third-party lib already handles
I/O, let it.

Specifically `QSaveFile` for output: `Chewer` writes results via
`QSaveFile` (typed-buffer → Arrow → `QSaveFile`, or whatever the layer
shape is). On crash or cancel, the temp file is just abandoned; the
target file is never partial.

The audit point (per the user, 2026-05-29): if h5-reader currently
end-runs Qt's file APIs for OUR-CODE file access (existing
`std::ifstream`/`fopen` usage), that's a deviation that should change
during the platform refactor. Third-party uses are exempt.

### UI freedom during chewer runs

While the chewer runs, the user can still **scrub the loaded protein
in the 3D view**. Because the data is shared-immutable, this works for
free:

- UI reads positions for frame M, renders
- Chewer reads positions for frame N (where N is somewhere in its
  iteration), runs transforms
- They read from the same buffers; no locks, no coordination
- The chewer's progress bar updates from coarse queued signals
- The UI's playback controller and the chewer's iteration counter are
  independent state

The 30-minute walk is fine. UI stays responsive throughout.

### Indexing/compute lag and natural level

Per the user (2026-05-29): "indexing and calculating can lag, find
their natural level."

- The chewer is greedy: it wants to chew through every (atom, frame)
  as fast as it can
- The UI is interactive: it wants to render the current frame as fast
  as the user needs
- These run at **independent rates** on their own threads
- No artificial throttling, no frame budgets, no backpressure machinery
- The compute finds its natural level given the work; the UI finds its
  natural level given user input
- They're decoupled by the shared-immutable-data contract

The Qt event system is for **coordination at sane rates**, not for
per-item flow. The chewer has its own work pump (plain C++ iteration);
the UI has Qt's event loop. Both have queues, served by appropriate
threads. The event system never has a 20-second slot to choke on
because no slot does 20 seconds of work — they're all coarse UI
updates.

---

## One walk, many transforms

The walk is the expensive thing — load + iterate + indexing. One walk
through 40K atoms × 20K frames cannot be repeated for each downstream
model. So:

**Register N transforms in one walk**; each emits to its own Parquet
sink. The single iteration feeds N outputs.

```bash
h5reader --batch --protein 1p9j.h5 \
  --transforms transforms/bs_aromatic_he.py \
               transforms/mc_backbone_hn.py \
               transforms/apbs_efg.py \
               transforms/larsen_hbond.py \
               transforms/ridge_baseline_features.py \
               transforms/mace_descriptor_features.py \
  --output-dir runs/1p9j_2026-06-XX/
```

One protein load. One traversal. Six Parquets. Each consumed
independently by its own downstream SDK script (PySR, ridge, MACE
training, etc.).

**The killer facility**: walking is rare and expensive; modeling is
frequent and cheap. The platform pays the walk once; the Parquets are
the durable contract for every subsequent model.

---

## Python is a data transform during the walk

Python scripts loaded by the chewer are **pure data transformations**:
given typed access to the substrate, emit rows. No model fitting, no
PySR, no sklearn, no MACE during the walk. The substrate is rich (the
"erector set" — KD-tree queries, topology queries, finite differences,
all available); the transform just outputs derived rows.

After the walk, a separate Python SDK (independent processes) consumes
the per-transform Parquets to fit models. This separation eliminates
GIL contention with PySR's Julia subprocess, PyTorch's thread pool,
etc., during the walk — modeling is post-walk, on its own time.

See `stage2_pysr_campaign_2026-05-29.md` for the transform schema and
SDK pattern.

---

## What stays put

- **Strip UI** (`AbstractStripPanel`, `StripStackWidget`, the
  dashboard signal infrastructure that just landed) — leave alone
- **Producer** (`src/` extractor) — not modified per
  `feedback_extractor_untouchable`
- **Existing typed model**: `QtProtein`, `QtTopology`, `QtFrame`,
  `QtTrajectoryH5`, `TrajectoryConformation`, the `QtRing` hierarchy
  — already the substrate the chewer needs
- **DftShielding cross-tool typed-JSON pattern** — model for cross-
  tool data exchange
- **CENSUS/ACONNECT/ASSERT_THREAD discipline** — chewer participates

## What changes

- **CMake refactor**: extract `h5reader-core` static library from the
  current monolithic target. Model, IO, diagnostics, and new analysis
  primitives live in the library. One `h5reader` binary still — the
  GUI binary links `h5reader-core + Qt Widgets + VTK + pybind11`.
- **New `src/analysis/` directory**: substrate primitives
  (`QtSpatialIndex` wrapping vendored nanoflann, `QtWindowDerivatives`,
  `QtTopologyQueries`), chewer machinery (`Chewer`, `QtTransformLoader`,
  `QtParquetSink`), pybind11 bindings (`bindings.cpp`), batch entry
  (`main_batch.cpp`).
- **New `extern/nanoflann/nanoflann.hpp`** (vendored single header,
  BSD license, ~5K LOC). Used for per-frame KD-trees in
  `QtSpatialIndex`. NOT VTK's spatial locators — VTK's
  rendering-pipeline overhead (vtkSmartPointer refcounting, vtkIdList
  allocations per query, vtkPolyData wrappers) is wasted weight for
  tight analytical loops at 750-20000 frames × per-frame index
  builds. nanoflann is the standard tool for low-dimensional spatial
  indexing in MD/computational-physics codebases (freud uses it under
  the hood). VTK stays for what it's actually good at: rendering,
  scene interaction, picking the molecular structure in the 3D view.
- **`main_batch.cpp` alongside `main_reader.cpp`**: the `--batch` flag
  dispatches to the batch main; otherwise the GUI main runs.
- **Optional pybind11 + Python embed**: gated `if(pybind11_FOUND)`;
  the chewer is built only when present (per the dft-ex1 pattern).
  Adviser distribution without Python is byte-identical to today.
- **Load on a worker thread**: `QtProteinLoader::Load()` runs on a
  worker; UI shows progress; protein appears as `shared_ptr<const
  LoadedProtein>` when ready. (Audit point: today this is likely
  synchronous on the main thread.)
- **File API audit + cleanup**: convert OUR-code file access to Qt
  APIs where natural (per the pragmatic rule above); leave third-party
  I/O alone.

## What's deferred

- **Movie capture per protein** (future capability; VTK already linked
  so no architecture change required)
- **Live in-app PySR overlays** (the "interactive payoff" pattern;
  works fine after the platform is solid)
- **Live experiment dock with cross-pollination** (different from a
  toolbar button + progress dialog; the cross-pollination case is for
  later)
- **Expert system + narrative engine** (Stage 3+)

---

# Substrate eager-precomputation set (Phase 1 + 2 substrate work)

The user's "prechewed everything" framing means the C++ substrate at
load time materializes a rich working surface so Python transforms
become small declarative scripts rather than feature-engineering
exercises. Compute time at load is NOT a constraint (user-confirmed:
"we don't actually care if it sits and chews for a day"); the
substrate optimizes for transform simplicity and PySR rediscovery odds,
not for load speed.

Settled per the 2026-05-30 substrate-tiers audit (Codex + Claude
general-purpose, with `feedback_appropriateness_lens` framing). The
list below is the merged tier set after agent feedback; see the
companion `substrate_conventions_2026-05-30.md` for the convention
calls that every feature here depends on.

## Phase 1 gate — conventions document precedes all substrate code

**Both audit agents independently demanded this.** Before any of the
features below is implemented, `substrate_conventions_2026-05-30.md`
must capture every convention call: spherical harmonic ordering and
phase, dipolar form sign, local frame handedness per atom class with
edge-case policies, charge source naming, residual baseline naming,
default-cutoff policy (none — required parameter everywhere),
multipole conventions, Cremer-Pople convention. Skipping this gate
means three high-impact features (SH, local frames, multipoles)
silently disagree across transforms.

The conventions doc is in place at
`analysis-speculative/substrate_conventions_2026-05-30.md` and serves
as the authoritative reference for every substrate accessor's
contract.

## Substrate features the chewer materializes

These sit alongside the baseline (positions, kernel time series,
per-frame KD-trees, ring centroids, H-bond candidates, derivatives,
top-K neighbors, ~10 standard derived series, DFT shielding access).

### Identity surface (Phase 2)

- **Atom-class chemistry predicates** — `atom.is_aromatic_HE()`,
  `is_backbone_HN()`, `is_carbonyl_O()`, `is_amide_N()`,
  `is_charged_N()`, `is_methyl_H()`, etc. Predicate dispatch only
  (no string keys). Already half-there in `QtAtom` typed identity;
  promote to documented public predicate surface.
- **Symmetry equivalence classes** — `traj.symmetry_class(atom)` and
  `traj.equivalent_atoms(atom)`. Computed from explicit topology
  rules (NOT geometry). Existing `QtAtom::equivalenceClass int8_t`
  promoted to typed enum. **No `symmetry_averaged_kernel` method**
  — averaging is a transform decision per audit feedback.
- **Prochiral markers** — `atom.prochirality_label` typed enum
  (`None`, `ProR`, `ProS`). Existing `QtAtom::diIndex` promoted to
  documented public surface. Required for methylene fits per
  `feedback_residual_as_ml_feature`.
- **Source identity columns** for every nearest/top-K feature —
  neighbor atom id, residue id, atom class, ring class. Plus a
  per-frame "did nearest identity change this frame" flag so
  transforms can detect neighbor-discontinuity artifacts that
  mislead PySR.

### Local frames per atom class (Phase 2)

- **`traj.local_frame(atom, frame_t) -> LocalFrame`** returning typed
  `LocalFrame { Vec3 z; Vec3 x; Vec3 y; bool is_valid; FrameVariant variant; }`.
  Right-handed throughout.
- Atom classes supported: HN, HE/HD (aromatic), HA, Cα, C=O. Edge-case
  policies for each (N-terminus, proline, pre-proline cis peptide, Gly
  prochiral pair, TRP bridgehead with dual ring membership, HID/HIE
  tautomer anchoring on chemistry-typed atom not protonated N) all
  documented in `substrate_conventions_2026-05-30.md`.
- **`traj.neighbor_in_local_frame(atom, target, frame_t) -> Vec3`** —
  the canonical anisotropic-feature primitive. Replaces transforms
  rebuilding residue geometry inline.

### Per-residue context joined to atoms (Phase 2)

- **`atom.residue_phi(frame_t)`, `residue_psi(frame_t)`,
  `residue_chi(frame_t, i)`** — radians. Convenience pair
  `residue_phi_sincos(frame_t) -> (sin, cos)` to avoid ±π
  discontinuity in PySR fits.
- **`atom.residue_dssp(frame_t) -> DsspCode`** — typed enum. Per
  convention: stratification only, NOT a PySR feature.
- **`atom.residue_is_h_bond_donor(frame_t)`,
  `is_h_bond_acceptor(frame_t)`** — derived from geometric H-bond
  candidates, NOT from DSSP.

### Geometric primitives + tensor operators (Phase 3)

Most of this is method-based (not pre-named columns) so the convention
is enforced in one C++ implementation per the audit recommendation:

- **`traj.dipolar_kernel(atom, source_pt, ref_dir, frame_t) -> float`**
  — returns `(3cos²θ − 1) / r³` per the conventions doc sign choice.
  The canonical anisotropic-scalar primitive.
- **`traj.dipolar_tensor(atom, source_pt, frame_t) -> Mat3`** —
  returns the full `K = (3 r̂r̂ᵀ − I) / r³` tensor in lab frame.
  Codex's catch; required for T2-form rediscovery per `feedback_t2_sacred`.
  Stage 2's thesis claim lives in T2; scalar Pople is necessary but
  not sufficient.
- **`traj.dipolar_tensor_T2(atom, source_pt, frame_t) -> SphericalTensor`**
  — the T2[5] decomposition matching `SphericalDecomposition.cpp`
  basis ordering.
- **`traj.dipolar_tensor_local(atom, source_pt, frame_t, frame_class) -> Mat3`**
  — same tensor projected into the atom's local frame for the
  specified atom class.
- **`traj.spherical_harmonics(atom, target, frame_t, l_max=2) -> numpy.ndarray`**
  — returns Y_l^m for the (atom → target) direction. SciPy real-form
  convention per the conventions doc. Skip Y₀⁰ (constant). Match
  `SphericalDecomposition.cpp` basis ordering for cross-consistency
  with T2 residual fits.
- **`traj.ring_susceptibility_projection(atom, ring, frame_t) -> Mat3`**
  — Claude's catch; useful for HM/BS distinction experiments.

### Integrated/aggregated kernels (Phase 3, Claude's catch)

- **`traj.integrated_dipolar(atom, target_predicate, cutoff_Å, frame_t) -> float`**
  — `Σ_{neighbors satisfying predicate, within cutoff} (3cos²θ − 1) / r³`.
  The form the integrated McConnell/EFG kernel takes in literature.
  PySR struggles to discover neighbor-enumeration; this gives it the
  integrated form as a single column. **High-value addition.**
- **`traj.integrated_dipolar_tensor(atom, target_predicate, cutoff_Å, frame_t) -> Mat3`**
  — tensor version (sum of per-pair K matrices).

### Charge multipoles (Phase 3, with charge source explicit)

Three separate methods per charge source (FF14SB / MOPAC / AIMNet2),
no fallback:

- **`traj.charge_within(atom, frame_t, cutoff_Å, charge_source) -> float`**
  — net charge in the local sphere (target atom excluded, optionally
  own-residue excluded)
- **`traj.charge_dipole_within(atom, frame_t, cutoff_Å, charge_source) -> Vec3`**
  — local charge cloud dipole moment
- **`traj.charge_quadrupole_within(atom, frame_t, cutoff_Å, charge_source) -> Mat3`**
  — traceless quadrupole tensor per the conventions doc
- **`traj.field_at(atom, frame_t, charge_source) -> Vec3`** — electric
  field at the atom from the charge cloud (Buckingham first-order
  feature)
- **`traj.field_gradient_at(atom, frame_t, charge_source) -> Mat3`**
  — EFG tensor (Buckingham second-order feature)

### DFT comparison surface (Phase 4)

- **`traj.dft_shielding(atom, frame_t) -> SphericalTensor or None`**
  — total shielding; T2-preserving
- **`traj.dft_shielding_components(atom, frame_t) -> {total, dia, para}`**
  — typed dict with all three components
- **`traj.dft_availability(atom) -> numpy.ndarray[bool]`** — per-frame
  mask
- **`traj.literature_residual(atom, frame_t, kernel_name, literature_constant_id) -> SphericalTensor`**
  — subtracts the literature equation's prediction
- **`traj.stage1_ridge_residual(atom, frame_t, kernel_name) -> SphericalTensor`**
  — subtracts what Stage 1's per-stratum ridge predicts

NEVER expose a method named `kernel_residual` — kernels aren't
shielding per `feedback_kernel_not_shielding`; "residual" requires
explicit baseline naming.

### Time-domain statistics (Phase 4)

Pre-computed per (atom, kernel) and per (atom, derived geometric series):
- mean, std (always)
- autocorrelation time per the conventions doc estimator (always)
- skew, kurtosis (always, but with a `n_samples_for_estimate` column
  so SDK scripts can downweight at 750-frame scale)
- peak frequency + spectral entropy (opt-in, marked by transform)

Sentinel-aware Welford per `feedback_conditional_welford_for_sentinels`
for series with NaN observations (e.g. distance-to-nearest when no
target in range).

### Rolling/lagged feature banks (Phase 4, Codex's catch)

- **`traj.kernel_lagged(atom, kernel_name, lag_frames) -> numpy.ndarray`**
- **`traj.kernel_rolling_mean(atom, kernel_name, window_frames) -> numpy.ndarray`**
- **`traj.kernel_rolling_std(atom, kernel_name, window_frames) -> numpy.ndarray`**
- Same for standard derived geometric series

Equation-forms section already names lagged/dynamics-averaged forms as
campaign targets; precomputing the named operators saves transforms
from rolling their own.

### Contact persistence summaries (Phase 4)

- **`traj.fraction_frames_within(atom, target_predicate, cutoff_Å) -> float`**
  — cutoff required, no default
- **`traj.mean_distance_to_nearest(atom, target_predicate) -> float`**
- **`traj.distance_to_nearest_percentile(atom, target_predicate, q) -> float`**
  — q ∈ [0, 1]
- **`traj.has_any_neighbor_ever(atom, target_predicate, cutoff_Å) -> bool`**
  — filter for unfittable atoms (Claude's catch)

### H-bond geometry in Larsen terms (Phase 4, Codex's catch)

- **`traj.h_bond_partners(atom, frame_t) -> list[HBondInfo]`** where
  `HBondInfo { Atom partner; HBondRole role; float distance_HA;
  float angle_DHA_rad; LarsenClass class; bool water_eligible; }`
- **`traj.h_bond_count(atom, frame_t) -> int`**

### Aligned RMSF (Phase 4)

- **`traj.aligned_rmsf(atom) -> float`** — scalar per atom after
  documented structural alignment to a reference frame
- **`traj.local_frame_rmsf(atom) -> Vec3`** — RMSF components in
  the atom's local frame (decomposes into "along bond" vs
  "perpendicular" mobility)

**Dropped from initial proposal** (both audit agents flagged as
cargo-cult; the 2026-05-30 canonical-rediscovery re-evaluation
confirmed drops because no canonical NMR shielding equation
references atom-path Frenet-Serret quantities):
- `curvature(atom)` — lab-frame Frenet-Serret is timestep- and
  thermostat-sensitive; not interpretable as shielding mechanism
- `torsion(atom)` — same
- `arc_length(atom)` — monotone in time; competes with `frame_index`
  as a trivial covariate
- Higher-power inverse distances (r⁻², r⁻⁴, r⁻⁵, r⁻⁷) — no canonical
  shielding equation uses these powers; offering them risks PySR
  finding spurious non-integer exponents (crank-discovery failure
  mode for a methods paper)
- DSSP as a PySR feature column — per-stratum rediscovery (run Pople
  separately on helix-HN vs sheet-HN, compare coefficients) achieves
  the Wishart-style finding without giving PySR a categorical
  shortcut. DSSP remains as a stratification axis only.

If a specific result later motivates any of these with explicit
canonical-literature anchor + fixture showing timestep-stability,
reconsider then.

### Per-pair distance tensor materialization at 1P9J scale (Phase 3, unconditional addition 2026-05-30)

At 1P9J's specific scale (~850 atoms × 1500 frames), the full
atom-pair distance tensor is `850² × 1500 × 8 bytes ≈ 8.7 GB` — fine
on the 128 GB RAM minimum. Substrate materializes it eagerly at load
and exposes:

- **`traj.pair_distance_tensor() -> numpy.ndarray[N_atom, N_atom, N_frame]`**
  — zero-copy view; symmetric along the first two axes
- **`traj.pair_distance(i, j, frame_t) -> float`** — O(1) lookup
  convenience
- **`traj.pair_distance_matrix(frame_t) -> numpy.ndarray[N_atom, N_atom]`**
  — per-frame slice

This is not a new feature exposure (doesn't give PySR new columns
PySR didn't already have); it's a storage/access optimization that
makes per-pair queries O(1) instead of per-call KD-tree query. At
larger protein scales (Stage 3 / fleet), the same API switches to
memory-mapped storage — that's a separate engineering call when
fleet scaling actually happens.

### Bond strain features (Phase 3, scope-gated → unconditional 2026-05-30)

Per the user direction "science lotto tickets" framing 2026-05-30:
the campaign explicitly adds **vibrational-correction rediscovery**
(Sundholm-Gauss-Schäfer 1996; Jameson's NMR perspective book) to its
kernel scope. Substrate exposes bond-strain features supporting this
rediscovery:

- **`traj.bond_strain(atom, bond_idx, frame_t) -> float`** — returns
  `(length_t − length_mean) / length_mean`. Data already in
  `BondLengthStatsTrajectoryResult` (`length_mean`, `length_std`,
  `length_min`, `length_max`, `length_delta_mean`, `length_delta_std`
  computed at frame 0); ~50 lines of accessor wrapping the existing
  buffers.
- **`traj.bond_strain_along_bond(atom, frame_t) -> Vec3`** — sum of
  strain vectors from all of atom's bonds, projected to per-atom
  Vec3. Direct PySR feature for atoms where bond-dipole modulation
  couples to shielding (HA via Buckingham first-order; HN via
  amide-bond π conjugation in McConnell).

Memory: 850 atoms × ~4 bonds × 1500 frames × 8 bytes ≈ 40 MB.
Trivial.

### Neighbor residue context (Phase 3, scope-gated → unconditional 2026-05-30)

Per the user direction "science lotto tickets" framing 2026-05-30:
the campaign explicitly adds **Wishart-Sykes-Richards 1991 chemical
shift index rediscovery** to its kernel scope. Substrate exposes
neighbor-residue context supporting this rediscovery:

- **`atom.neighbor_residue_phi(frame_t, offset) -> float`** for offset
  ∈ {−2, −1, +1, +2}
- **`atom.neighbor_residue_psi(frame_t, offset) -> float`** for same
- **`atom.neighbor_residue_chi(frame_t, offset, chi_index) -> float`**
  for sidechain context

Plus sin/cos convenience pairs to avoid ±π discontinuity per the
conventions doc.

Memory: 850 atoms × 1500 frames × ~12 scalars (4 offsets × 3
angles) × 8 bytes ≈ 120 MB. Small.

**Boundary handling**: at chain termini, neighbor accessors return
NaN for missing offsets. Transforms handle via the standard
sentinel-aware Welford pattern.

## What changes (revised)

This expands the prior "What changes" section. The substrate's new
work items, dependency-ordered:

1. **`substrate_conventions_2026-05-30.md`** committed first (the
   Phase 1 gate, already done)
2. **CMake refactor**: extract `h5reader-core` static library; one
   `h5reader` binary still
3. **Identity exposures** in `src/analysis/identity/`: atom-class
   predicates, symmetry classes, prochiral markers, source identity
   columns. Mostly wrappers over existing `QtAtom` / `QtTopology`
   typed identity; promote to documented public surface.
4. **Per-residue context** accessors in `src/analysis/per_residue/`:
   φ, ψ, χᵢ, sin/cos forms, DSSP, H-bond donor/acceptor flags
5. **`LocalFrame` struct + per-class `local_frame()` methods** in
   `src/analysis/local_frame/`. Edge-case handling per conventions
   doc. Per-class fixture tests for HN, HE-aromatic, HA-glycine, Cα,
   C=O handedness — all in `tests/local_frame_tests.cpp`.
6. **Geometric methods + tensor operators** in `src/analysis/dipolar/`:
   `dipolar_kernel`, `dipolar_tensor`, `dipolar_tensor_T2`,
   `dipolar_tensor_local`, `spherical_harmonics`, integrated
   variants, ring susceptibility projection. Y_l^m fixture against
   known direction `Y₂⁰(ẑ) = (1/4)√(5/π)`.
7. **DFT access + residual surfaces** in `src/analysis/dft/`:
   `dft_shielding`, `dft_shielding_components`, `dft_availability`,
   `literature_residual`, `stage1_ridge_residual`. T2-preserving
   throughout per `feedback_t2_sacred`.
8. **Cremer-Pople ring pucker** in `src/analysis/ring_pucker/`: port
   the corrected formula from `PlanarGeometryResult.cpp` with the
   regression test per `feedback_adversarial_review_physics`.
9. **H-bond geometry** in `src/analysis/h_bond/`: typed `HBondInfo`,
   Larsen classifications, water eligibility per the geometric
   H-bond candidate machinery.
10. **Charge multipoles** in `src/analysis/charge_multipoles/`: three
    charge-source-explicit method families.
11. **Contact persistence + time-domain stats** in
    `src/analysis/stats/`: aggregation helpers; sentinel-aware
    Welford reused from existing implementation.
12. **Rolling/lagged feature banks** in `src/analysis/temporal/`.
13. **Aligned RMSF** in `src/analysis/rmsf/`. Document the alignment
    reference (probably first frame, or median-frame, or
    minimum-energy frame; pick one per conventions).
14. **pybind11 bindings** in `src/analysis/bindings/`: typed Python
    classes exposing all of the above.
15. **Provenance sidecar** writer in `src/analysis/manifest/`:
    `manifest.json` per walk capturing every convention call (JSON via `QJsonDocument`, not YAML).

The original "What changes" items (CMake refactor, main_batch.cpp,
pybind11 + Python embed, async loader, file API cleanup, reader-
untouchable analog) remain unchanged.

---

## Reader-untouchable analog

Per `feedback_extractor_untouchable` (producer not modified during
viewer/reader feature work), the platform refactor has an analogous
discipline:

**Stage 2 platform additions live in `src/analysis/` and new top-level
files (`main_batch.cpp`, `bindings.cpp`). The GUI surface (`src/app/`
widgets, scenes, overlays, REST endpoints, playback controller) is
NOT modified unless the change is explicitly an interactive-payoff
feature.**

The CMake split makes this enforceable: `h5reader-core` test builds
without anything in `src/app/` proves no leakage. The GUI tests
(`h5reader_rest_smoke`) catch regressions in the GUI surface.

---

## The audit deliverable

Agents are tasked with: read existing h5-reader code; compare to this
design; identify gaps per dimension. The dimensions:

1. **I/O async-ness** — is load synchronous on the main thread? are
   per-frame NPYs loaded synchronously?
2. **Thread affinity** — what's the affinity of every QObject?
   anywhere `moveToThread` is called? any `Qt::QueuedConnection` ?
3. **Frame-advance fanout** — how many consumers connect to
   `frameChanged`? what's the cost?
4. **VTK pipeline efficiency** — one actor or N actors per molecule?
   incremental or rebuild per frame? bounds cache implications?
5. **Per-frame allocation profile** — heap allocations / QObject
   constructions per frame change
6. **Memory layout** — SoA vs AoS; cache-friendliness for the chewer's
   linear iteration vs the UI's random access
7. **Backpressure** — does fast scrubbing pile up render events?
8. **Lifecycle and shutdown** — clean shutdown? race conditions? clean
   thread join?
9. **Diagnostics** — per-frame cost visibility? queue depth metrics?
10. **File access** — Qt APIs for our code? third-party usage
    appropriate?
11. **Scalability projection** — 5000 atoms × 20K frames, what breaks
    first?

For each dimension: current state in named files/classes, target
shape per this doc, gap assessment (small/medium/large), specific
changes. Then synthesis: overall distance, highest-leverage changes,
risk to current functionality, sequence that works (not a phase
ladder; actual dependency ordering).

The agent prompt is composed separately; this doc is the
specification it audits against.

---

## Cross-references

- Companion: `stage2_pysr_campaign_2026-05-29.md` — the Stage 2
  campaign that consumes this platform
- `feedback_extractor_untouchable` — analogous reader-touchable
  discipline
- `feedback_qt_discipline` — CENSUS/ACONNECT/ASSERT_THREAD,
  UDP-log-first
- `feedback_qt_citizen` — portable Qt, no platform libs, respect
  Qt/VTK idioms
- `feedback_vtk_bounds_cache` — known VTK pipeline confusion;
  symptom of the architecture being audited
- `feedback_udp_logging` — UDP for per-evaluation diagnostics; stderr
  for startup/summary
- `feedback_surface_complex_data` — observability principle for the
  interactive features when they earn their weight
- `project_h5reader_thesis_vetting_tool` — h5-reader role as trust
  boundary; adviser-facing reliability requirement
- dft-ex1 precedent at `/shared/dft-ex1/` — pybind11 embed conditional
  pattern (`SCF_ENABLE_PYTHON`), sibling typed JSON exchange model
- `project_h5_north_star` — Qt-native trajectory protein architecture
  reference
