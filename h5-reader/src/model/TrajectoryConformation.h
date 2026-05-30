// TrajectoryConformation — the H5-backed conformation: N frames of dense,
// eager-loaded per-TR time-series (positions + ~56 calculator groups), with
// sparse per-frame snapshot detail loaded as current-frame source data through
// the base Conformation seam.
//
// This is the former QtConformation, lifted under the Conformation base. It
// keeps the dense layer it always owned (QtTrajectoryH5 + the QtFrame view);
// the base adds the snapshot buffer. Windowed-H5 streaming, when it lands,
// belongs HERE (or in QtTrajectoryH5) — never in the base facade.

#pragma once

#include "Conformation.h"

#include "../io/QtTrajectoryH5.h"
#include "../io/TrajectoryFrameMap.h"
#include "QtFrame.h"

#include <QString>
#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

namespace h5reader::model {

class QtProtein;
class QtConformationSnapshot;

class TrajectoryConformation final : public Conformation {
    Q_OBJECT

public:
    // perFrameNpysDir: absolute path to the run's frame snapshot NPY directory,
    // or empty if the --trajectory run did not emit per-frame NPYs (then
    // snapshots simply never become resident -- dense H5 still drives playback).
    TrajectoryConformation(const QtProtein* protein,
                           std::unique_ptr<h5reader::io::QtTrajectoryH5> h5,
                           QString perFrameNpysDir);
    ~TrajectoryConformation() override;

    // ----- Conformation seam -----
    std::size_t frameCount() const override { return h5_->frameCount(); }
    double timePicoseconds(std::size_t frame) const override;
    // H5 row -> original (XTC) frame index via the frame_indices map; identity
    // fallback if the row is out of range. This is the key the DFT job dirs and
    // per-frame npys are named by.
    std::size_t originalFrameIndex(std::size_t frame) const override {
        return h5reader::io::TrajectoryFrameMap::OriginalIndex(frame, *h5_);
    }
    Vec3 atomPosition(std::size_t frame, std::size_t atomIdx) const override;
    const TrajectoryConformation* asTrajectory() const override { return this; }

    // ----- Dense H5 layer (trajectory-only) -----
    const h5reader::io::QtTrajectoryH5* h5() const { return h5_.get(); }
    QtFrame frame(std::size_t t) const;  // thin dense per-atom view over the H5

    // Nearest H5 frame-row that has a per-frame snapshot on disk — lets the UI
    // point at where full per-atom detail IS sampled when the cursor sits
    // between emit-stride frames. nullopt if the run emitted no frame snapshots.
    std::optional<std::size_t> nearestSampledFrame(std::size_t frame) const;
    const std::vector<std::size_t>& sampledFrameRows() const { return sampledRows_; }

protected:
    std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t frame) override;

private:
    std::unique_ptr<h5reader::io::QtTrajectoryH5> h5_;
    QString perFrameNpysDir_;
    std::vector<std::size_t> sampledRows_;  // sorted H5 rows that have a frame snapshot dir
};

}  // namespace h5reader::model
