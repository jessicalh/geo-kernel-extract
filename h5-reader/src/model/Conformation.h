// Conformation — abstract base: a viewable sequence of >=1 frames of one
// protein, AND the working buffer for per-frame full-fidelity detail.
//
// Two concrete subclasses ("subclass, don't hack" — feedback_vtk_subclassing):
//   TrajectoryConformation — N frames backed by trajectory.h5 (dense, eager).
//   SingleConformation     — one pose, no H5 (a --orca / --mutant / --pdb run).
//
// WHY this is a QObject (the manager role lives here — decided 2026-05-26, see
// memory project_h5reader_buffering_decision_20260526): the conformation owns
// the snapshot buffer and emits the readiness signal. The buffered piles
// (QtConformationSnapshot) stay PLAIN immutable-after-load data behind
// shared_ptr<const> — no QObject, no thread affinity, so a future prefetch
// worker can fill one and hand it to the GUI lock-free. "pin" a frame = keep
// its shared_ptr (the LRU cannot evict what a holder keeps alive); "observe" =
// drop it.
//
// WHY it stays a thin FACADE (codex review, 2026-05-26): the DENSE per-frame
// data (positions, time-series) lives in the subclass
// (TrajectoryConformation / QtTrajectoryH5). When windowed-H5 streaming lands
// for the microsecond trajectories, its residency policy stays THERE — never
// here. This base manages only the sparse snapshot buffer + the frame-access
// seam, so it never becomes a god-manager of two residency policies.
//
// Cross-object holders of a Conformation* MUST use QPointer<Conformation>
// (crash-diagnosis.md QPointer discipline) — it is a QObject with a lifetime
// shorter than the app.

#pragma once

#include "Types.h"  // Vec3

#include <QObject>

#include <cstddef>
#include <deque>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace h5reader::model {

class QtProtein;
class QtConformationSnapshot;
class TrajectoryConformation;

class Conformation : public QObject {
    Q_OBJECT

public:
    ~Conformation() override;

    Conformation(const Conformation&) = delete;
    Conformation& operator=(const Conformation&) = delete;

    // ----- Identity / topology spine (shared, loaded once) -----
    const QtProtein* protein() const { return protein_; }
    std::size_t ringCount() const;

    // ----- Frame-access seam (both subclasses satisfy) -----
    virtual std::size_t frameCount() const = 0;
    virtual double timePicoseconds(std::size_t frame) const = 0;

    // Per-atom position for `frame`. The ONE position seam both run shapes
    // share — rendering (MoleculeScene) and geometry overlays read this:
    // trajectory → the resident H5; single pose → the snapshot's Pos column.
    virtual Vec3 atomPosition(std::size_t frame, std::size_t atomIdx) const = 0;

    // Typed downcast without RTTI ("objects answer questions about
    // themselves"): non-null only for the H5-backed trajectory. Consumers that
    // need the dense per-frame time-series (the ribbon's per-frame DSSP, the
    // time-series dock, QtFrame) gate on this; single-pose returns nullptr and
    // those features simply do not engage.
    virtual const TrajectoryConformation* asTrajectory() const { return nullptr; }

    // ----- Full-fidelity detail (the buffered pile) -----
    // Cached-or-null; NEVER blocks. null == not resident: call
    // requestSnapshot() and react to snapshotReady(). Holding the returned
    // shared_ptr pins the snapshot against LRU eviction.
    std::shared_ptr<const QtConformationSnapshot> snapshot(std::size_t frame) const;

    // Ensure `frame`'s snapshot is loaded, then emit snapshotReady(frame).
    // v1: SYNCHRONOUS — loads + caches + emits before returning. The committed
    // prefetch increment moves the load onto a worker behind this SAME
    // signature; consumers already cope with snapshot()==null and react to
    // snapshotReady, so they do not change. Idempotent for a resident frame.
    void requestSnapshot(std::size_t frame);

signals:
    // Emitted when snapshot(frame) has become non-null (resident or just
    // loaded). A failed load logs at the loader seam and emits nothing.
    void snapshotReady(std::size_t frame);

protected:
    // cacheCapacity = max resident snapshots in the LRU. SingleConformation
    // passes 1 (its one pose, effectively pinned).
    Conformation(const QtProtein* protein, std::size_t cacheCapacity);

    // Where a frame's snapshot comes from — the ONE virtual that varies
    // behaviour (no factories, no handles, no pluggable store: honours the
    // no-pluggable-interfaces rule). Returns null on load failure (logged at
    // the seam). Called OUTSIDE the cache lock (it does file IO).
    virtual std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t frame) = 0;

    const QtProtein* protein_ = nullptr;

private:
    // Frame-keyed LRU of immutable snapshots. The mutex is present from day one
    // even though v1 loads synchronously on the GUI thread: the shared_ptr
    // control block is atomic, but the map/deque bookkeeping is not, so the
    // moment prefetch goes async this is already the guarded container (one
    // fewer 2am use-after-free — crash-diagnosis discipline).
    mutable std::mutex cacheMutex_;
    std::unordered_map<std::size_t, std::shared_ptr<const QtConformationSnapshot>> cache_;
    std::deque<std::size_t> lru_;  // front = least-recently-used
    std::size_t cacheCapacity_;

    void touch_(std::size_t frame);  // move to MRU; caller holds cacheMutex_
};

}  // namespace h5reader::model
