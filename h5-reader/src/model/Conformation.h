// A viewable sequence of one or more conformations of one protein.
// TrajectoryConformation reads trajectory.h5; SingleConformation holds one
// pose. Dense trajectory data stays in those subclasses. This base owns only
// the current per-frame calculator snapshot and its readiness signal.
//
// Snapshot access is synchronous and confined to this QObject's thread.
// Long-lived QObject owners should hold Conformation through QPointer.

#pragma once

#include "Types.h"  // Vec3

#include <QObject>

#include <cstddef>
#include <memory>
#include <optional>

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

    // Identity and topology shared by every frame.
    const QtProtein* protein() const { return protein_; }
    std::size_t ringCount() const;

    // Frame access implemented by both concrete subclasses.
    virtual std::size_t frameCount() const = 0;
    virtual double timePicoseconds(std::size_t frame) const = 0;

    // Original trajectory index shared by H5 rows, frame NPY directories, and
    // DFT jobs. A single pose uses the identity mapping.
    virtual std::size_t originalFrameIndex(std::size_t frame) const { return frame; }

    // Resolve an original trajectory frame identity to the exact row accepted
    // by atomPosition(), timePicoseconds(), and requestSnapshot(). The base
    // identity implementation covers a single pose; sampled trajectories
    // override this with their H5 frame map.
    virtual std::optional<std::size_t>
    frameRowForOriginalIndex(std::size_t originalFrame) const {
        if (originalFrame >= frameCount() || originalFrameIndex(originalFrame) != originalFrame)
            return std::nullopt;
        return originalFrame;
    }

    // Per-atom position for `frame`. The ONE position seam both run shapes
    // share — rendering (MoleculeScene) and geometry overlays read this:
    // trajectory → the resident H5; single pose → the snapshot's Pos column.
    virtual Vec3 atomPosition(std::size_t frame, std::size_t atomIdx) const = 0;

    // Non-null only for an H5-backed trajectory. Consumers of dense time-series
    // data use this to distinguish a trajectory from a single pose.
    virtual const TrajectoryConformation* asTrajectory() const { return nullptr; }

    // Current resident snapshot, or null when another frame is resident.
    std::shared_ptr<const QtConformationSnapshot> snapshot(std::size_t frame) const;

    // Load `frame` synchronously when needed, then emit snapshotReady(frame).
    // Re-requesting the resident frame is idempotent and still emits the signal.
    void requestSnapshot(std::size_t frame);

signals:
    // Emitted when snapshot(frame) has become non-null (resident or just
    // loaded). A failed load logs at the loader seam and emits nothing.
    void snapshotReady(std::size_t frame);

protected:
    Conformation(const QtProtein* protein);

    // Load one frame from the concrete backing store. Failures are logged at
    // the loader boundary and represented by a null pointer.
    virtual std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t frame) = 0;

    const QtProtein* protein_ = nullptr;

private:
    std::optional<std::size_t> residentSnapshotFrame_;
    std::shared_ptr<const QtConformationSnapshot> residentSnapshot_;
};

}  // namespace h5reader::model
