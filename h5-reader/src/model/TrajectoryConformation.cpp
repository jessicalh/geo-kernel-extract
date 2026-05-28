// TrajectoryConformation — implementation.

#include "TrajectoryConformation.h"

#include "../io/FrameNpyLoader.h"
#include "QtConformationSnapshot.h"
#include "QtProtein.h"

#include <QDir>
#include <QFileInfo>

#include <algorithm>
#include <cstdint>
#include <unordered_map>

namespace h5reader::model {

TrajectoryConformation::TrajectoryConformation(const QtProtein* protein,
                                               std::unique_ptr<h5reader::io::QtTrajectoryH5> h5,
                                               QString perFrameNpysDir)
    : Conformation(protein),
      h5_(std::move(h5)),
      perFrameNpysDir_(std::move(perFrameNpysDir)) {
    setObjectName(QStringLiteral("TrajectoryConformation"));

    // Index which H5 rows have a per-frame snapshot dir, so the UI can point at
    // the nearest sampled frame when parked between emit-stride frames. Scan the
    // frame_NNNNNN dirs once (O(emitted), not O(frames)) and map each
    // original index back to its H5 row via frame_indices.
    if (!perFrameNpysDir_.isEmpty()) {
        const auto& fidx = h5_->frameIndices();
        std::unordered_map<std::uint64_t, std::size_t> origToRow;
        origToRow.reserve(fidx.size());
        for (std::size_t row = 0; row < fidx.size(); ++row)
            origToRow.emplace(fidx[row], row);
        const QStringList dirs = QDir(perFrameNpysDir_)
                                     .entryList(QStringList{QStringLiteral("frame_*")},
                                                QDir::Dirs | QDir::NoDotAndDotDot);
        for (const QString& d : dirs) {
            bool ok = false;
            const std::uint64_t orig = d.mid(6).toULongLong(&ok);  // skip the "frame_" prefix
            if (!ok)
                continue;
            const auto it = origToRow.find(orig);
            if (it != origToRow.end())
                sampledRows_.push_back(it->second);
        }
        std::sort(sampledRows_.begin(), sampledRows_.end());
    }
}

TrajectoryConformation::~TrajectoryConformation() = default;

double TrajectoryConformation::timePicoseconds(std::size_t frame) const {
    const auto& ft = h5_->frameTimes();
    return frame < ft.size() ? ft[frame] : 0.0;
}

Vec3 TrajectoryConformation::atomPosition(std::size_t frame, std::size_t atomIdx) const {
    const auto* pos = h5_->positions();
    return pos ? pos->at(atomIdx, frame) : Vec3::Zero();
}

QtFrame TrajectoryConformation::frame(std::size_t t) const {
    return QtFrame(this, t);
}

std::optional<std::size_t> TrajectoryConformation::nearestSampledFrame(std::size_t frame) const {
    if (sampledRows_.empty())
        return std::nullopt;
    const auto hi = std::lower_bound(sampledRows_.begin(), sampledRows_.end(), frame);
    if (hi == sampledRows_.begin())
        return *hi;
    if (hi == sampledRows_.end())
        return sampledRows_.back();
    const auto lo = std::prev(hi);
    return (frame - *lo <= *hi - frame) ? *lo : *hi;
}

std::shared_ptr<const QtConformationSnapshot> TrajectoryConformation::loadSnapshot(std::size_t frame) {
    if (perFrameNpysDir_.isEmpty())
        return nullptr;  // run emitted no per-frame NPY snapshots; detail unavailable

    // The frame dir is keyed by the ORIGINAL (XTC) frame index, zero-padded to
    // six digits (frame_NNNNNN), matching FrameNpyEmitter's layout. The H5 row
    // `frame` maps to that original index via frame_indices.
    const auto& idx = h5_->frameIndices();
    const std::size_t orig = frame < idx.size() ? static_cast<std::size_t>(idx[frame]) : frame;
    const QString dir = QStringLiteral("%1/frame_%2").arg(perFrameNpysDir_).arg(orig, 6, 10, QLatin1Char('0'));

    // Per-frame NPYs may be emitted at a stride. Check the documented frame
    // directory name and return null silently when absent, so scrubbing past
    // unsampled frames does not spam the loader's
    // "directory does not exist" report. A present-but-malformed dir still
    // reports through the FrameNpyLoader seam.
    if (!QFileInfo::exists(dir))
        return nullptr;

    return h5reader::io::FrameNpyLoader::LoadSnapshotDir(dir, protein_, orig, timePicoseconds(frame));
}

}  // namespace h5reader::model
