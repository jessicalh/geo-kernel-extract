#pragma once

#include "QtTrajectoryH5.h"

#include <QString>

#include <cstddef>
#include <vector>

namespace h5reader::io {

class TrajectoryFrameMap {
public:
    // Scan perFrameNpysDir for frame_NNNNNN/ subdirectories; return sorted H5
    // rows that have a snapshot dir on disk. Empty perFrameNpysDir returns
    // empty vector (--trajectory runs may not emit per-frame NPYs).
    static std::vector<std::size_t>
    ScanSampledRows(const QString& perFrameNpysDir, const QtTrajectoryH5& h5);

    // H5 row -> original XTC frame index via the frame_indices map; identity
    // fallback if the row is out of range.
    static std::size_t
    OriginalIndex(std::size_t row, const QtTrajectoryH5& h5);
};

}  // namespace h5reader::io
