# Substrate calculator framework

Date: 2026-05-30
Status: settled design. Companion to
`reader_as_platform_2026-05-29.md` (the platform architecture that
uses this framework), `substrate_conventions_2026-05-30.md` (the
convention calls each calculator must respect), and
`stage2_pysr_campaign_2026-05-29.md` (the campaign whose transforms
consume the calculators' outputs).

## Why this exists

Substrate calculators differ from nmr_extract calculators in three
load-bearing ways. The framework exists so we don't reinvent those
three differences in every calculator implementation.

**nmr_extract calculator pattern** (one-shot, batch):
- Given a conformation, compute the result, write to H5, done.
- No cancellation — runs to completion or fails.
- No query surface — output is an H5 group consumers read later.
- No multi-stage lifecycle — load + compute is one operation.

**Substrate calculator pattern** (multi-phase, in-process):
- **Load phase**: input data from ANY of the four input sources
  becomes available (see "What loads from where" below). The
  calculator initializes.
- **Compute phase**: eager precomputation across all frames.
  **Cancellable** — the chewer's atomic flag is checked at frame
  boundaries.
- **Query phase**: Python transforms query the precomputed data via
  the calculator's pybind11 bindings. Reads are lock-free (data is
  immutable post-compute).
- **Cleanup**: destroyed when the protein is dropped.

## DRY discipline — reuse h5-reader's existing typed structures and loaders

**Hard rule (settled 2026-05-30 after concrete code-reading pass)**:
substrate calculators NEVER roll their own structural view of the
data. Every typed data structure they consume comes from h5-reader's
existing model layer. Every loader call goes through h5-reader's
existing static loader functions. The chewer adds bulk-parallel
orchestration around the existing loaders, NOT replacement parsers
or replacement typed structures.

The viewer (the existing `Conformation` / `TrajectoryConformation` /
`DftShieldingStore` machinery) is **NOT modified**. Its design
discipline — "one resident frame at a time; persistent history lives
in strip ChannelBuffers, not here" (per `Conformation.h:13-14` and
`DftShieldingStore.h:18-20`) — stays intact. The chewer goes AROUND
this residency policy, not THROUGH it. The chewer owns its own
`std::vector<shared_ptr<const QtConformationSnapshot>>` storage,
populated via parallel calls to the existing `FrameNpyLoader::LoadSnapshotDir`
static method.

### What is reused as-is from h5-reader

| h5-reader component | Reused for |
|---|---|
| `h5reader::io::FrameNpyLoader::LoadSnapshotDir(dir, protein, frameIdx, timePs)` | Per-frame NPY snapshot loading; chewer calls in parallel loop over all frames. ALREADY STATIC; directly callable. |
| `h5reader::io::QtTrajectoryH5(path)` | trajectory.h5 reader. Constructor takes a path; eager bulk load by design. Chewer instantiates one per protein. |
| `h5reader::io::QtTopologySidecar` | Topology sidecar loader. Direct construction. |
| `h5reader::model::QtConformationSnapshot` | Per-frame snapshot data type. Pure-data class (NOT QObject), `shared_ptr<const>` semantics, designed for worker-thread fill (Conformation.h:12-14: "PLAIN immutable-after-load data behind shared_ptr<const> — no QObject, no thread affinity, so a future worker can fill one and hand it to the GUI"). The chewer IS that future worker. |
| `h5reader::model::QtProtein`, `QtTopology`, `QtAtom`, `QtRing`, `QtResidue` | Typed identity model; every typed query reuses these. |
| `h5reader::model::DftShieldingFrame` (in `DftShielding.h`) | Pure-data DFT frame type; chewer's vector holds `shared_ptr<const DftShieldingFrame>`. |
| `h5reader::io::OrcaShieldingParser` | DFT .out file parser; chewer uses for every DFT frame. |
| `h5reader::io::QtProteinLoader::Load(h5_path)` | Top-level orchestrator returning `QtLoadResult { QtProtein, Conformation }`. Chewer uses for the protein; discards the bundled Conformation (it builds its own bulk storage). |

### What is NOT reused (and why)

| h5-reader component | Why the chewer doesn't reuse |
|---|---|
| `h5reader::model::Conformation` | One-resident-at-a-time residency policy explicitly forbids accumulating snapshots (Conformation.h:13-14). The chewer needs all frames resident; goes around. |
| `h5reader::model::TrajectoryConformation` | Same residency reason; also bundles snapshot lifecycle that the chewer doesn't need (chewer has all snapshots already). |
| `h5reader::model::DftShieldingStore` | Same one-resident-at-a-time discipline (DftShieldingStore.h:18-20). Chewer builds its own DFT frame vector. |

The new code in the substrate is small: a `LoadedProtein` carrier
struct + ~30-line loader calculators that call the existing static
loaders in parallel + the new derived calculators (KD-trees, local
frames, dipolar tensors, bond strain, neighbor context, etc.) that
the viewer doesn't have. Everything else reuses h5-reader's typed
data layer directly.

## h5-reader changes required BEFORE chewer work lands

Two small refactors in h5-reader to enable DRY reuse. Both preserve
the viewer's behavior unchanged.

### Change A: extract DFT load+validate to a static helper

`DftShieldingStore::loadAndValidate(originalIndex)` is currently
**private**. Extract to a static helper in a sibling header so the
chewer can call it directly without going through `DftShieldingStore`
(which has the one-resident discipline).

Proposed shape:

```cpp
// h5-reader/src/io/DftShieldingLoader.h (new)
namespace h5reader::io {
class DftShieldingLoader {
public:
    // Parse + validate one DFT job. Returns null and logs at the seam
    // on any failure (missing meta.json, parser hole, atom-count
    // mismatch, dia+para != total).
    static std::shared_ptr<const h5reader::model::DftShieldingFrame>
    LoadAndValidate(std::size_t originalIndex,
                    const QString& jobsDir,
                    const h5reader::model::QtProtein* protein);
};
}
```

`DftShieldingStore::loadAndValidate` becomes a thin wrapper that
calls `DftShieldingLoader::LoadAndValidate`. Viewer behavior
unchanged; chewer can now bulk-load all DFT frames in parallel by
calling the static helper.

Effort: ~30 lines moved + ~10 lines wrapper. Half-hour.

### Change B: extract sampled-rows + original-frame-index logic

`TrajectoryConformation` computes two things in its constructor that
the chewer also needs:
- `sampledRows_` — which H5 rows have per-frame NPY snapshots (scans
  the perFrameNpysDir).
- `originalFrameIndex(frame)` — H5 row → original XTC index mapping
  used as the key for NPY dir names and DFT job dir names.

Extract both into static helpers callable independently:

```cpp
// h5-reader/src/io/TrajectoryFrameMap.h (new)
namespace h5reader::io {
class TrajectoryFrameMap {
public:
    // Scan perFrameNpysDir; return sorted H5 rows that have a
    // frame_NNNNNN/ subdirectory.
    static std::vector<std::size_t>
    ScanSampledRows(const QString& perFrameNpysDir, const QtTrajectoryH5& h5);

    // H5 row -> original XTC frame index (from frame_indices map).
    static std::size_t
    OriginalIndex(std::size_t row, const QtTrajectoryH5& h5);
};
}
```

`TrajectoryConformation`'s constructor calls these instead of doing
the computation inline. Viewer behavior unchanged; chewer reuses the
same helpers.

Effort: ~30 lines moved. Half-hour.

### Change C (deferred — use the existing path with Conformation discarded)

`QtProteinLoader::Load(h5_path)` returns `QtLoadResult { QtProtein,
Conformation }`. The chewer wants the protein + paths only, not the
Conformation. Decision: use the existing Load, take the protein,
discard the Conformation. Minor one-time cost (one TrajectoryConformation
constructed then dropped) vs the cost of a new sibling API
(QtProteinLoader::LoadForBulk). Defer adding the sibling API unless
profiling shows the discard cost matters.

### Order

Changes A and B land BEFORE any chewer code. Both are mechanical,
both preserve viewer behavior, both have zero risk to the existing
codebase. They go in one commit ("Extract DFT load+validate and
trajectory-frame-map helpers for chewer reuse") with the existing
viewer tests confirming no behavior change.

## What loads from where (input source map)

The substrate has FOUR input sources, not just the trajectory H5.
Different calculators consume different sources; the framework
exposes them all on the `LoadedProtein` carrier so calculators just
ask for what they need.

| Source | What it carries | Loaded by |
|--------|----------------|-----------|
| `trajectory.h5` | positions per frame, the 13 kernel TR groups (BS, HM, MC, APBS EFG, water field, hydration shell, MOPAC variants, Larsen H-bond, dispersion, etc.), Welford rollups, the 5 new TRs (iRED, KernelDynamics, ReorientationalDynamics, DihedralAutocorrelation, KernelCoherence), DSSP/dihedral/ring-pucker per-frame, BondLengthStatsTrajectoryResult | `PositionsLoader`, `KernelTRLoader` (selectively per `REQUIRED_KERNELS`), `BondLengthStatsLoader`, etc. |
| **Topology sidecar** | atoms (identity, element, residue), bonds, rings (with member_atoms), residues, **static FF14SB charges from prmtop**, ring chemistry classification | `TopologyLoader` |
| **Per-frame NPYs** (`per_frame_npys/`) | per-atom snapshots with the heavier per-frame data: **MOPAC charges per frame**, **AIMNet2 charges + 256-dim embeddings per frame**, **APBS electric field per atom per frame**, anything else the producer emits as per-frame full-fidelity snapshots | `PerFrameNpyLoader` |
| **DFT job directory** (`dft/jobs/<job>/`) | Per-job `meta.json` orchestration file (`files.out_primary` points at the successful ORCA output) + the ORCA `.out` file itself carrying the shielding tensors. Per-frame paths come from the calcset's `.LGS` `dft.frames[]` array (`frame_index → meta_json`); see `spec/CALCSET_MANIFEST.md`. Parsed via the existing `OrcaShieldingParser` (NOT a JSON parser for the shielding data — meta.json is orchestration only). 750 DFT-stride frames at every other NPY frame. | `DftShieldingLoader` (Change A — wraps the existing `OrcaShieldingParser`) |

**Critical difference from h5-reader**: the per-frame NPY snapshots
are LAZY in h5-reader (loaded one at a time on user scrub via
`Conformation::requestSnapshot`). In the chewer's "prechewed
everything" framing, they're **eagerly loaded** at the start of the
compute phase. For 1P9J at ~850 atoms × 1500 NPY frames × ~few MB
per snapshot, that's ~few GB resident — comfortable on the 128 GB RAM
minimum.

`LoadedProtein` is therefore a carrier struct holding the immutable
results of bulk-parallel calls to h5-reader's existing static loaders:

```cpp
struct LoadedProtein {
    std::unique_ptr<h5reader::model::QtProtein>                                  protein;
        // from QtProteinLoader::Load (taking the protein, discarding the bundled Conformation)
    std::unique_ptr<h5reader::io::QtTrajectoryH5>                                h5;
        // direct construction with h5_path (extracted from the QtLoadResult)
    std::vector<std::shared_ptr<const h5reader::model::QtConformationSnapshot>>  snapshots;
        // populated via FrameNpyLoader::LoadSnapshotDir in parallel loop;
        // indexed by H5 row; null for rows with no NPY snapshot on disk
    std::vector<std::shared_ptr<const h5reader::model::DftShieldingFrame>>       dft_frames;
        // populated via DftShieldingLoader::LoadAndValidate (Change A) in parallel loop;
        // indexed by ORIGINAL XTC frame index; null for absent
    std::vector<std::size_t> sampled_rows;     // from TrajectoryFrameMap::ScanSampledRows (Change B)
    std::vector<std::size_t> original_indices; // H5 row -> original XTC index, per Change B
    // Immutable after load. shared_ptr<const LoadedProtein> handed to all consumers.
};
```

**Every type referenced above is an existing h5-reader type.** The
chewer constructs none of its own data structures for the loaded
data. The DRY discipline applies all the way down.

Each calculator's `onLoad()` and `onCompute()` receives this carrier
and reads whichever fields it needs. Calculators that need positions
read `loaded.h5->position(atomIdx, t)` (the existing QtTrajectoryH5
API); calculators that need per-frame AIMNet2 charges read
`loaded.snapshots[t]->column(io::FieldKind::Aimnet2Charge)` via
QtConformationSnapshot's existing typed-column accessor; calculators
that compute DFT residuals read DFT shielding from `loaded.dft_frames[origIdx]`.

**Loader calculators are thin wrappers around h5-reader's existing
static loaders.** Each is typically ~30-50 lines: parallel reads via
`QtConcurrent::map` over a frame range, calling the existing static
loader per frame, storing in the `LoadedProtein`'s vector. The example
for `PerFrameNpyLoader`:

```cpp
// src/analysis/calculators/per_frame_npy_loader.cpp (the WHOLE compute fn)
void PerFrameNpyLoader::onCompute(LoadedProtein& loaded,
                                  QPromise<ChewerEvent>& promise) {
    log_phase_start("compute");
    const auto t0 = std::chrono::steady_clock::now();

    const std::size_t n_frames = loaded.h5->frameCount();
    loaded.snapshots.assign(n_frames, nullptr);

    // sampled_rows tells us which H5 rows have a per-frame NPY directory
    // on disk (via the existing Change-B helper). Skip the rest.
    //
    // Canonical Qt parallelism: QtConcurrent::map over a dedicated
    // QThreadPool (passed in via the chewer's ChewerConfig.parallelPool).
    // Each iteration is a per-frame snapshot load via the existing
    // FrameNpyLoader::LoadSnapshotDir static helper. The map iterates
    // in parallel; we check promise.isCanceled() per iteration and skip
    // the read when set.
    std::atomic<std::size_t> progress_counter{0};
    QFuture<void> mapFuture = QtConcurrent::map(
        *loaded.parallelPool,
        loaded.sampled_rows,
        [&](std::size_t row) {
            if (promise.isCanceled()) return;
            const std::size_t orig = loaded.original_indices[row];
            const QString npyDir = npyDirForOriginalIndex(perFrameNpysDir_, orig);
            // Direct call to the existing static loader — no rewriting of
            // the NPY parser or the typed snapshot structure.
            loaded.snapshots[row] = h5reader::io::FrameNpyLoader::LoadSnapshotDir(
                npyDir, loaded.protein.get(), orig,
                loaded.h5->timePsAtFrame(row));
            const std::size_t done = ++progress_counter;
            if (done % 100 == 0)
                log_progress("frame=" + std::to_string(done) +
                             " of " + std::to_string(loaded.sampled_rows.size()));
        });
    mapFuture.waitForFinished();

    const auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    log_phase_done("compute", ms);
}
```

`QtConcurrent::map` is the canonical Qt primitive for parallel
collection processing. It respects the pool we pass in, integrates
with the rest of Qt's concurrency story (the inner `QFuture<void>`
could also be cancelled independently if needed), and avoids the
`std::execution::par_unseq` problems with side-effectful Qt I/O
inside the loop body that an earlier draft flagged.

That's the entire compute phase for a loader calculator. The actual
NPY parsing, the FieldKind resolution, the typed `QtConformationSnapshot`
construction, the per-column data layout — ALL of that is in
`FrameNpyLoader::LoadSnapshotDir`, which already exists in h5-reader
and serves the viewer. The chewer adds the parallel-loop orchestration
and stores the results in the chewer's own vector. Same pattern for
`DftShieldingLoader` (parallel calls to `DftShieldingLoader::LoadAndValidate`
per Change A), `KernelTRLoader` (extracts requested groups from the
already-loaded `QtTrajectoryH5`), etc.

**Downstream calculators depend on the loaders** by name and consume
the existing typed structures. E.g. `ChargeMultipoles_FF14SB` reads
charges from `loaded.protein->ff14sbCharges()` (existing API on QtProtein);
`ChargeMultipoles_MOPAC` reads from
`loaded.snapshots[row]->column(io::FieldKind::MopacCharge)` (existing
QtConformationSnapshot accessor); `ChargeMultipoles_AIMNet2` reads
from the AIMNet2 FieldKind column on the same snapshots; `DftResiduals`
reads DFT shielding from `loaded.dft_frames[origIdx]` (existing
`DftShieldingFrame` type) and kernel values from
`loaded.h5->bsShielding(...)` (existing QtTrajectoryH5 accessor).

**No new typed structures for any of the loaded data.** The chewer's
typed surface to Python via pybind11 is built on h5-reader's existing
typed model layer.

The chewer has ~20+ such calculators (per-frame KD-tree, ring
centroids, H-bond candidates, local frames per atom class, dipolar
tensors, integrated dipolar sums, ring susceptibility projection,
charge multipoles in three charge sources × three moments, DFT
residual surfaces, per-residue context, contact persistence,
time-domain stats, rolling/lagged feature banks, aligned RMSF, bond
strain, neighbor residue context, symmetry classes, prochiral
markers, etc.). Without a framework, each calculator reinvents
lifecycle, cancel discipline, progress reporting, and pybind11
exposure.

## The interface

```cpp
class SubstrateCalculator {
public:
    // Quick — record sizes, allocate output buffers, no heavy work yet.
    virtual void onLoad(const LoadedProtein& protein) = 0;

    // Do the work. MUST check promise.isCanceled() at frame boundaries
    // (or finer if a single frame's work is long). MAY call
    // promise.setProgressValue(framesDone) for numeric progress and
    // promise.addResult(ChewerEvent{...}) for typed events to the GUI.
    // MAY emit UDP progress via the inherited log_* helpers for
    // human-readable debug output alongside the QPromise events.
    //
    // Canonical Qt pattern: QPromise carries cancel state + result
    // stream; QFuture::cancel() (called from any thread by the owner)
    // sets the promise's cancel state atomically; the calculator polls
    // it via promise.isCanceled().
    virtual void onCompute(const LoadedProtein& protein,
                           QPromise<ChewerEvent>& promise) = 0;

    // Calculator identity + dependency graph (for the scheduler)
    virtual std::string name() const = 0;
    virtual std::vector<std::string> depends_on() const = 0;

    virtual ~SubstrateCalculator() = default;

protected:
    // UDP progress helpers (use the existing StructuredLogger)
    void log_phase_start(const std::string& phase) const;
    void log_progress(const std::string& message) const;
    void log_phase_done(const std::string& phase, double elapsed_ms) const;
    void log_cancel() const;
    void log_error(const std::string& message) const;
};
```

The base class is **deliberately small**: three required virtual
methods + five logging helpers. No std::function callbacks. No
template metaprogramming. No registration macros. The calculator's
pybind11 bindings are NOT a virtual method on the base — each
calculator's query surface is its own; the binding code lives in the
same file as the calculator.

## Concrete example: RingCentroids

The whole calculator fits in one file. Reader knows everything about
it by reading this one file.

```cpp
// src/analysis/calculators/ring_centroids.h
#pragma once
#include "substrate_calculator.h"
#include "../../model/Vec3.h"

namespace h5reader::analysis {

class RingCentroids : public SubstrateCalculator {
public:
    std::string name() const override { return "RingCentroids"; }
    std::vector<std::string> depends_on() const override {
        return {"PositionsLoader", "TopologyLoader"};
    }

    void onLoad(const LoadedProtein& protein) override;
    void onCompute(const LoadedProtein& protein,
                   QPromise<ChewerEvent>& promise) override;

    // Query surface (Python-visible; declared here, bound below)
    Vec3 centroid_at(size_t ring_id, size_t frame_t) const;
    const std::vector<Vec3>& centroid_series(size_t ring_id) const;

private:
    size_t n_rings_ = 0;
    size_t n_frames_ = 0;
    std::vector<std::vector<Vec3>> centroids_;  // [ring_id][frame_t]
};

}  // namespace h5reader::analysis
```

```cpp
// src/analysis/calculators/ring_centroids.cpp
#include "ring_centroids.h"
#include "../substrate_calculator_bindings.h"  // helper for pybind11 inline
#include <chrono>

namespace h5reader::analysis {

void RingCentroids::onLoad(const LoadedProtein& loaded) {
    n_rings_ = loaded.protein->ringCount();         // existing QtProtein API
    n_frames_ = loaded.h5->frameCount();             // existing QtTrajectoryH5 API
    centroids_.assign(n_rings_, std::vector<Vec3>(n_frames_));
    log_phase_start("load");
    log_phase_done("load", 0.0);
}

void RingCentroids::onCompute(const LoadedProtein& loaded,
                              QPromise<ChewerEvent>& promise) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    log_phase_start("compute");

    for (size_t t = 0; t < n_frames_; ++t) {
        // Canonical Qt cancel check: QPromise::isCanceled() is an
        // atomic read internally. QFuture::cancel() from any thread
        // sets it; we poll at frame boundaries.
        if (promise.isCanceled()) {
            log_cancel();
            return;
        }
        for (size_t r = 0; r < n_rings_; ++r) {
            const auto& ring = loaded.protein->ring(r);  // existing QtProtein API
            Vec3 c{0, 0, 0};
            for (auto atomIdx : ring.member_atoms) {
                // existing QtTrajectoryH5 position accessor — atom-major buffer
                c = c + loaded.h5->position(atomIdx, t);
            }
            centroids_[r][t] = c / static_cast<double>(ring.member_atoms.size());
        }
        // Numeric progress for the watcher's progressValueChanged
        // signal (drives the GUI's QProgressDialog).
        promise.setProgressValue(static_cast<int>(t));
        if (t % 100 == 0) {
            log_progress("frame=" + std::to_string(t) + " of " + std::to_string(n_frames_));
        }
    }

    const auto t1 = clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    log_phase_done("compute", ms);
}

Vec3 RingCentroids::centroid_at(size_t ring_id, size_t frame_t) const {
    return centroids_.at(ring_id).at(frame_t);
}

const std::vector<Vec3>& RingCentroids::centroid_series(size_t ring_id) const {
    return centroids_.at(ring_id);
}

}  // namespace h5reader::analysis

// pybind11 bindings — inline, same file
namespace py = pybind11;
void bind_ring_centroids(py::module_& m) {
    py::class_<h5reader::analysis::RingCentroids>(m, "RingCentroids")
        .def("centroid_at", &h5reader::analysis::RingCentroids::centroid_at)
        .def("centroid_series", &h5reader::analysis::RingCentroids::centroid_series,
             py::return_value_policy::reference_internal);
}
```

```cpp
// tests/calculators/ring_centroids_tests.cpp
#include "../../src/analysis/calculators/ring_centroids.h"
#include <QtTest>

class RingCentroidsTests : public QObject {
    Q_OBJECT
private slots:
    void emptyRingTrivialCase();
    void singleAromaticRingPHE();
    void cancellationStopsAtFrameBoundary();
};
// ... tests follow ...
```

That's the full lifecycle for one calculator in three files — header,
cpp+bindings, tests. ~150 lines total. Reader knows the input, the
output, the dependencies, the cancel behavior, the progress log
shape, and the pybind11 surface by reading these three files. Nothing
hidden in a central registry.

## The scheduler

The framework's scheduler runs inside `RunChewer` (see the platform
doc's "Execution model" section). It does three things only:

1. Topologically sorts calculators by `depends_on()` declarations.
   Errors loudly on cycles or unknown dependencies.
2. Calls `onLoad()` for every calculator in dep order.
3. Calls `onCompute()` for every calculator in dep order, passing the
   `QPromise<ChewerEvent>&` for cancel checks + progress reporting.
   May parallelize independent calculators via `QtConcurrent::map`
   over a sub-pool when their `depends_on()` sets don't overlap.

```cpp
class SubstrateCalculatorScheduler {
public:
    void register_calculator(std::unique_ptr<SubstrateCalculator> calc);

    // Both phases take the QPromise that QtConcurrent::run handed to
    // the surrounding RunChewer function. promise carries cancel +
    // progress + event-stream channels back to the owner.
    void run_load_phase(const LoadedProtein& protein,
                        QPromise<ChewerEvent>& promise);
    void run_compute_phase(const LoadedProtein& protein,
                           QPromise<ChewerEvent>& promise);

    template <typename T>
    T& get(const std::string& name) const;  // for pybind11 binders
private:
    std::vector<std::unique_ptr<SubstrateCalculator>> calculators_;
    std::vector<size_t> compute_order_;  // computed at register-time
};
```

~250 lines for the entire scheduler + base class + UDP logging
helpers. The scheduler is **boring on purpose** — no clever dispatch,
no graph optimization, no async magic. Cancel propagation is just
"pass the QPromise through to each calculator's onCompute." Read it
once, trust it.

## Cancel propagation

The cancel path is the canonical Qt pattern: the chewer's owner calls
`chewerFuture.cancel()`; QPromise's cancel state is set atomically;
the calculator polls `promise.isCanceled()` inside its loop. No
custom CancelToken class, no separate atomic flag, no signal/slot for
cancel.

```cpp
// Inside any calculator's onCompute:
for (size_t t = 0; t < n_frames_; ++t) {
    if (promise.isCanceled()) {           // canonical Qt cancel check
        log_cancel();
        return;
    }
    // ... do work ...
}
```

Two distinct cancel layers exist, both already in the 3-level cancel
design (per the platform doc):

- **C++ compute-phase cancel**: `promise.isCanceled()` polled at
  frame boundaries inside `onCompute()`. Cooperative. The owner's
  `chewerFuture.cancel()` (called from any thread — typically the
  GUI thread via the Cancel button or `aboutToQuit` handler) sets
  the promise's cancel state atomically. The calculator sees it on
  its next check; bails cleanly, logs `cancel`, returns. No Python
  involvement — the 3-level escalation doesn't apply at this layer
  because there's no Python interpreter to interrupt.
- **Python transform-phase cancel** (after compute is done): the
  full 3-level pattern (cooperative `promise.isCanceled()` check
  before each transform iteration + `PyErr_SetInterruptEx(SIGINT)`
  if the Python script doesn't notice + log+abandon if even that
  doesn't take) as documented in the platform doc. C++ substrate is
  already finished computing; this layer protects the Python
  iteration only.

Both layers route through the same QPromise the calculator's
`onCompute` received. No alternate channel; no class duplication.

The chewer's `requestCancel()` slot routes to the right layer
depending on which phase is active. Inside compute phase: just set
the C++ atomic. Inside transform phase: do the 3-level escalation.

## Logging — UDP primary, sync sinks optional (zero-write hot path)

Per user direction (2026-05-30): logging is **UDP-primary, with
synchronous sinks (stderr, file) explicitly optional and OFF by
default for the chewer's hot path**. The point is zero-write logging
— the chewer must not pay a synchronous disk/stderr write per progress
emit. UDP is fire-and-forget (`QUdpSocket::writeDatagram`); it doesn't
block on consumer readiness; it's the only logging primitive cheap
enough to be in the chewer's per-100-frames inner-loop emit path.

**Existing reader's StructuredLogger behavior**: today emits to BOTH
stderr (synchronous `fprintf` + `fflush`) AND UDP on every call. The
synchronous stderr write is fine for the GUI app's startup/shutdown
event rate; it would throttle the chewer's emit rate. Required change
to the existing logger: a per-emit-call sink-suppression knob (e.g.
`StructuredLogger::Emit(..., SinkPolicy { include_stderr=false,
include_file=false })` with the default for chewer-path calls being
"UDP only").

Format JSON (matches the existing logger's output) — each datagram is
a self-contained JSON object with structured fields the receiver
parses directly:

```json
{"ts":"2026-05-30T...","severity":"INFO","category":"calculator",
 "thread":"chewer-worker","calculator":"RingCentroids",
 "phase":"compute","frame":372,"total":750}
```

Examples of what the chewer emits:

```json
{"calculator":"RingCentroids","phase":"load"}
{"calculator":"RingCentroids","phase":"compute"}
{"calculator":"RingCentroids","frame":100,"total":750}
{"calculator":"RingCentroids","frame":200,"total":750}
{"calculator":"RingCentroids","frame":300,"total":750}
...
{"calculator":"RingCentroids","phase":"done","elapsed_ms":4231}
```

Receivers parse JSON directly. The existing `udp_listen.py` helper
pretty-prints these for humans; downstream tooling can ingest the
JSON stream without re-parsing a "human-readable" line that was just
the JSON's `message` field. (Earlier drafts of this doc proposed
`key=value` text — corrected 2026-05-30 after verifying the actual
logger format is JSON.)

Why not numeric percentages: a 750-frame compute phase progressing
through 49.6% conveys less than `frame=372 of 750 (RingCentroids)`,
which names the calculator AND the position. The latter is also
filterable. Numeric percent bars are a UI affordance (the
QFutureWatcher's `progressValueChanged` signal drives the progress
dialog from `QPromise::setProgressValue`) — the structured log is
for diagnostics, not UI.

**Sink policy decided per call site**:
- **Chewer's per-100-frames `log_progress` calls**: UDP only, no
  sync sinks. Hot path; can fire thousands of times per run.
- **Chewer's `log_phase_start` / `log_phase_done` / `log_cancel` /
  `log_error` calls**: UDP + stderr (low frequency; once per calculator;
  the human running --batch wants to see them).
- **Reader's existing startup/shutdown / error / lifecycle logging**:
  UDP + stderr (the GUI app's existing behavior; unchanged).
- **File sink**: never default-on. Available as an `--log-file
  PATH` CLI option for batch runs where the user explicitly wants
  the JSON stream persisted to disk for later analysis. Writes via
  `QSaveFile` so the file isn't corrupted on crash.

**Concurrency note**: the existing logger serializes via mutex; the
chewer's parallel calculator threads emitting concurrently throttle
on the mutex. For the per-100-frames rate (~tens of emits per second
total across all parallel calculators), this is well within the
mutex's capacity. If we ever need much higher emit rates, the
canonical Qt answer is `qInstallMessageHandler` + a lock-free SPSC
queue draining to UDP on a dedicated logger thread — defer until
profiling demands it.

The base class's `log_*` helpers do the UDP write with the right sink
policy; each calculator emits whatever its specific progress shape
is. Helpers are inherited from `SubstrateCalculator` so calculators
don't need to know about the logger directly.

## Readability discipline

This is the load-bearing section. The framework's whole point is
making the 20+ calculators easy to read, modify, and maintain. The
following anti-patterns are forbidden by convention:

**Anti-patterns to avoid**:

- **Deep inheritance hierarchies.** `SubstrateCalculator` is the only
  base. No intermediate "AbstractRingCalculator" or
  "TensorCalculatorBase" — duplication beats indirection at this
  scale.
- **Template metaprogramming for dependency tracking.** Dependencies
  are declared as a `std::vector<std::string>` returned from
  `depends_on()`. If a calculator gets a typo, the scheduler errors
  loudly at registration time.
- **std::function callbacks in the lifecycle.** The three virtual
  methods are concrete; tracing control flow means reading the
  scheduler's `run_compute_phase()` — which is ~30 lines you can read
  in one sitting.
- **Registration macros that hide the dependency graph.**
  `register_calculator()` is a plain function call. The chewer's
  initialization explicitly creates and registers each calculator —
  the registration block is the canonical place to see the
  substrate's dependency graph.
- **Cross-cutting concerns hidden in base classes.** The base class
  is small and concrete; if you need to know what a calculator does
  during compute, you read the calculator, not the base.
- **Per-binding files (`bindings.cpp` central registry).** pybind11
  bindings live in the same .cpp as the calculator. The chewer's
  top-level `bind_all_calculators(py::module_& m)` function calls
  each calculator's `bind_*(m)` function — one line per calculator,
  no hidden registration.
- **String-keyed runtime config**. If a calculator has parameters
  (e.g., a window size for a smoothing helper), they're constructor
  arguments to typed C++ types, not loaded from a config file at
  runtime. Constants in the conventions doc are constants in the
  code.

**Patterns to prefer**:

- **One calculator per file.** Header + cpp+bindings + tests in a
  parallel test file. Reader knows everything about the calculator
  from those three files.
- **Explicit dependency declarations.** A static list of calculator
  names; the scheduler validates them. Typos surface at registration,
  not at runtime.
- **Linear data flow inside each phase.** `onCompute()` reads
  immutable inputs, writes the calculator's own output buffers, exits.
  No callbacks, no events emitted to other calculators.
- **Concrete progress logging.** Each calculator emits whatever
  progress shape makes sense for its loop structure (frame-based,
  ring-based, atom-pair-based, etc.). Helpers from the base class do
  the UDP write.
- **Tests in `tests/calculators/<name>_tests.cpp` parallel to the
  calculator.** Per-calculator test target; runs in the standard
  ctest suite. CI catches regressions per calculator.
- **The framework code itself stays small.** ~300 lines total for
  scheduler + base + cancel + logging. Adding the 21st calculator
  doesn't grow the framework.

## File organization

```
h5-reader/
├── src/analysis/
│   ├── substrate_calculator.h          # ~50 lines: base class
│   ├── substrate_calculator.cpp        # ~80 lines: log helpers
│   ├── substrate_calculator_scheduler.h
│   ├── substrate_calculator_scheduler.cpp   # ~150 lines: topo sort + run
│   ├── cancel_token.h
│   ├── bind_calculators.cpp            # ~30 lines: calls each bind_* fn
│   └── calculators/                    # one file per calculator
│       │   # Loaders (run first; pull from the four input sources)
│       ├── topology_loader.{h,cpp}             # topology sidecar
│       ├── positions_loader.{h,cpp}            # trajectory.h5
│       ├── kernel_tr_loader.{h,cpp}            # trajectory.h5 (selective per REQUIRED_KERNELS)
│       ├── per_frame_npy_loader.{h,cpp}        # per_frame_npys/ (MOPAC + AIMNet2 + APBS + embeddings)
│       ├── dft_shielding_loader.{h,cpp}        # DFT JSON sidecar
│       ├── bond_length_stats_loader.{h,cpp}    # trajectory.h5 BondLengthStats
│       │   # Geometric primitives (depend on positions + topology)
│       ├── per_frame_kd_tree.{h,cpp}
│       ├── ring_centroids.{h,cpp}
│       ├── h_bond_candidates.{h,cpp}
│       ├── local_frames.{h,cpp}                # all atom classes in one file
│       ├── dipolar_tensors.{h,cpp}
│       ├── integrated_dipolar_sums.{h,cpp}
│       ├── ring_susceptibility_projection.{h,cpp}
│       │   # Charge / field calculators (each per charge source)
│       ├── charge_multipoles_ff14sb.{h,cpp}    # consumes topology
│       ├── charge_multipoles_mopac.{h,cpp}     # consumes per-frame NPYs
│       ├── charge_multipoles_aimnet2.{h,cpp}   # consumes per-frame NPYs
│       │   # DFT comparison
│       ├── dft_residuals.{h,cpp}               # consumes DFT JSON + kernel computations
│       │   # Per-residue context + symmetry
│       ├── per_residue_context.{h,cpp}
│       ├── neighbor_residue_context.{h,cpp}
│       ├── symmetry_classes.{h,cpp}
│       ├── prochiral_markers.{h,cpp}
│       │   # Trajectory-level derived quantities
│       ├── contact_persistence.{h,cpp}
│       ├── time_domain_stats.{h,cpp}
│       ├── rolling_lagged_banks.{h,cpp}
│       ├── aligned_rmsf.{h,cpp}
│       ├── bond_strain.{h,cpp}
│       └── per_atom_derived_series.{h,cpp}     # the 10+ canonical series
└── tests/calculators/                  # parallel test tree
    ├── positions_loader_tests.cpp
    ├── topology_loader_tests.cpp
    └── ... (one per calculator)
```

~20-25 calculators, each in a self-contained file. The whole substrate
is `src/analysis/calculators/` + ~5 framework files + parallel tests.

## Provenance integration

Each calculator's per-run metadata feeds the chewer's `manifest.json`
provenance sidecar:

```yaml
calculators:
  - name: RingCentroids
    version: 1
    elapsed_ms: 4231
    cancelled: false
    frames_processed: 750
    output_size_bytes: 18000
  - name: LocalFrames
    version: 1
    elapsed_ms: 8420
    cancelled: false
    frames_processed: 750
    output_size_bytes: 102000
  # ... per calculator
```

The base class's logging helpers record these stats; the scheduler
aggregates them into the manifest at the end of compute phase. Methods
paper supplementary materials can cite which calculators ran with
what timing.

## What's NOT in this framework

- **Producer-style one-shot calculators.** The existing nmr_extract
  pattern stays as-is in `src/`; not modified per
  `feedback_extractor_untouchable`.
- **Python-side transform infrastructure.** Transforms are scripts in
  Python that USE calculator outputs; they don't subclass
  `SubstrateCalculator`. The chewer's transform-running machinery is
  separate.
- **Dynamic calculator registration from config.** The chewer's
  initialization code explicitly creates each calculator. No plugins,
  no DLL loading, no config-driven enablement. To add a calculator,
  add a file + add one line to the chewer's init.
- **Per-calculator GUI integration.** Calculators are headless; if a
  calculator's output should be displayable, that's a separate
  h5-reader UI feature using the same data via the same pybind11
  bindings.
- **Async/coroutine-based calculator implementations.** The scheduler
  may run independent calculators in parallel via `QtConcurrent::map`,
  but each calculator's `onCompute()` is a synchronous function. No
  C++ coroutines, no continuation chains.

## Cross-references

- `reader_as_platform_2026-05-29.md` — the platform architecture
  that owns the chewer that owns the scheduler that owns the
  calculators
- `substrate_conventions_2026-05-30.md` — every convention call each
  calculator must respect (SH basis ordering, dipolar form,
  local-frame definitions, charge sources, default-cutoff policy,
  etc.)
- `stage2_pysr_campaign_2026-05-29.md` — the campaign whose
  transforms query the calculators' outputs
- `feedback_qt_discipline` — CENSUS_REGISTER / ACONNECT /
  ASSERT_THREAD discipline that applies to QObject calculators (most
  calculators are NOT QObjects; only those with Qt-signal surfaces
  need to be)
- `feedback_udp_logging` — UDP port 9997 is the primary debug
  channel; calculator progress logs ride the same socket
- `feedback_extractor_untouchable` — nmr_extract producer pattern is
  preserved; this framework is for the consumer side only
- `feedback_methods_accumulate` — methodologies coexist; calculators
  for different methods (e.g., Larsen H-bond geometric vs DSSP
  H-bond detection) run side-by-side, never replace
