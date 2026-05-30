#include "TrajectoryFrameMap.h"

#include <QDir>
#include <QStringList>

#include <algorithm>
#include <cstdint>
#include <unordered_map>

namespace h5reader::io {

std::vector<std::size_t>
TrajectoryFrameMap::ScanSampledRows(const QString& perFrameNpysDir, const QtTrajectoryH5& h5) {
    std::vector<std::size_t> sampledRows;
    if (perFrameNpysDir.isEmpty())
        return sampledRows;

    const auto& fidx = h5.frameIndices();
    std::unordered_map<std::uint64_t, std::size_t> origToRow;
    origToRow.reserve(fidx.size());
    for (std::size_t row = 0; row < fidx.size(); ++row)
        origToRow.emplace(fidx[row], row);

    const QStringList dirs = QDir(perFrameNpysDir)
                                 .entryList(QStringList{QStringLiteral("frame_*")},
                                            QDir::Dirs | QDir::NoDotAndDotDot);
    for (const QString& d : dirs) {
        bool ok = false;
        const std::uint64_t orig = d.mid(6).toULongLong(&ok);  // skip the "frame_" prefix
        if (!ok)
            continue;
        const auto it = origToRow.find(orig);
        if (it != origToRow.end())
            sampledRows.push_back(it->second);
    }
    std::sort(sampledRows.begin(), sampledRows.end());
    return sampledRows;
}

std::size_t TrajectoryFrameMap::OriginalIndex(std::size_t row, const QtTrajectoryH5& h5) {
    const auto& idx = h5.frameIndices();
    return row < idx.size() ? static_cast<std::size_t>(idx[row]) : row;
}

}  // namespace h5reader::io
