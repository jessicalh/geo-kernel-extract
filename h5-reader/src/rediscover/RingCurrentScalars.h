#pragma once

#include "LiteratureConstants.h"

#include "../model/Types.h"

#include <cmath>
#include <cstddef>
#include <limits>

namespace h5reader::rediscover {

static_assert(static_cast<int>(model::RingTypeIndex::PheBenzene) == 0, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::ProPyrrolidine)
              == model::kAromaticRingTypeCount,
              "saturated proline ring must stay outside aromatic per-type BS slots");

inline double RingForwardContributionPpm(const model::SphericalTensor& literatureScaledKernel) {
    return literatureScaledKernel.T0;
}

inline double RingPerTypeT0Ppm(const double* unitCurrentT0ByType, std::size_t count) {
    if (!unitCurrentT0ByType || count < static_cast<std::size_t>(model::kAromaticRingTypeCount))
        return std::numeric_limits<double>::quiet_NaN();

    double out = 0.0;
    for (int t = 0; t < model::kAromaticRingTypeCount; ++t) {
        const double unitT0 = unitCurrentT0ByType[static_cast<std::size_t>(t)];
        if (!std::isfinite(unitT0)) return std::numeric_limits<double>::quiet_NaN();
        out += unitT0 * RingIntensity(static_cast<model::RingTypeIndex>(t)).value;
    }
    return out;
}

}  // namespace h5reader::rediscover
