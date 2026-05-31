// AnalysisBody — the immutable resident body captured by rediscover closures.

#pragma once

#include "Catalog.h"
#include "ResidentIndexes.h"
#include "RunData.h"

namespace h5reader::rediscover {

struct Body {
    const RunData& run;
    const ResidentIndexes& idx;
    const Catalog& catalog;
};

}  // namespace h5reader::rediscover
