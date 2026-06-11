#pragma once

#include "AnalysisBody.h"

#include <QString>

#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

struct SpineProbeDatasetResult {
    bool passed = false;
    QString dataset_label;
    QString manifest_path;
    std::size_t field_count = 0;
    std::size_t failed_fields = 0;
    std::vector<QString> failed_stems;
};

SpineProbeDatasetResult RunSpineReachabilityProbe(const Body& body,
                                                  const QString& datasetLabel,
                                                  const QString& outDir,
                                                  QString* err_out = nullptr);

}  // namespace h5reader::rediscover
