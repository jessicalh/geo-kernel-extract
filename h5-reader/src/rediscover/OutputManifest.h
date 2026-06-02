// OutputManifest — per-run manifest for the rediscover substrate outputs.

#pragma once

#include "ExtractionSupport.h"
#include "RecordSink.h"

#include <QMap>
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

    // Provenance (#51 §4). The row-alignment contract (which lets a consumer
    // load a sidecar NPY without the CSV) + per-feature support counts (nonzero /
    // present-row counts per emitted feature family) so a downstream fitter can
    // report effective-N WITHOUT re-deriving it (a Python recompute is forbidden).
    // Empty ⇒ not provided for this relationship.
    QString rowAlignmentContract;
    QString normalization;                 // e.g. "raw" — un-normalized lab frame
    QMap<QString, std::size_t> featureSupport;
};

bool WriteOutputManifest(const QString& outDir, const std::vector<OutputEntry>& outputs,
                         const DftFrameAlignment& alignment, int rc, QString* err_out = nullptr);

}  // namespace h5reader::rediscover
