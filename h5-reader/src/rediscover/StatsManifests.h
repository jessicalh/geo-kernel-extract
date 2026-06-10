#pragma once

#include "ConditioningSpec.h"
#include "RowDesign.h"
#include "RowDesignEmitter.h"
#include "RunData.h"

#include <QString>

namespace h5reader::rediscover {

struct RowDesignManifestContext {
    QString datasetId;
    QString rootPath;
    bool fixture = false;
};

bool WriteRowDesignManifests(const QString& outDir,
                             const std::vector<RowColumnSpec>& schema,
                             const RowDesignStats& stats,
                             const RunData& run,
                             const ConditioningSpec& spec,
                             const RowDesignManifestContext& context,
                             QString* err_out = nullptr);

}  // namespace h5reader::rediscover
