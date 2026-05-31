// OutputManifest — per-run manifest for the rediscover substrate outputs.

#pragma once

#include "ExtractionSupport.h"
#include "RecordSink.h"

#include <QString>
#include <QStringList>

#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

struct OutputEntry {
    QString relationship;
    QString relationshipKind;
    QString sourcesCsv;
    QString aggregatedCsv;
    QStringList sidecarNpys;
    std::size_t cases = 0;
    std::size_t sourceRows = 0;
    std::size_t aggregatedRows = 0;
};

bool WriteOutputManifest(const QString& outDir, const std::vector<OutputEntry>& outputs,
                         const DftFrameAlignment& alignment, int rc, QString* err_out = nullptr);

}  // namespace h5reader::rediscover
