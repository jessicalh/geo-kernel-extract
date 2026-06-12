#pragma once

#include "ConditioningSpec.h"

#include <QString>

#include <cstddef>

namespace h5reader::rediscover {

struct ConsolidatedEmitOptions {
    QString root720;
    QString run1p9j;
    QString outDir;
    ConditioningSpec conditioningSpec = ConditioningSpec::Default();
};

struct ConsolidatedEmitStats {
    std::size_t rows = 0;
    std::size_t scopedFieldCount = 0;
};

bool RunConsolidatedEmit(const ConsolidatedEmitOptions& options,
                         ConsolidatedEmitStats* stats,
                         QString* err_out);

}  // namespace h5reader::rediscover
