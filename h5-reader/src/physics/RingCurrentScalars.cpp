#include "RingCurrentScalars.h"

#include <cmath>
#include <limits>

namespace h5reader::physics {

double RingPerTypeT0Ppm(const double* unitCurrentT0ByType, std::size_t count) {
    static_assert(static_cast<int>(model::RingTypeIndex::PheBenzene) == 0,
                  "ring intensity slot mismatch");
    static_assert(static_cast<int>(model::RingTypeIndex::ProPyrrolidine)
                      == model::kAromaticRingTypeCount,
                  "saturated proline ring must stay outside aromatic slots");

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

model::SphericalTensor ScaleSphericalTensor(const model::SphericalTensor& t, double scale) {
    model::SphericalTensor out;
    out.T0 = scale * t.T0;
    for (int i = 0; i < 3; ++i)
        out.T1[static_cast<std::size_t>(i)] = scale * t.T1[static_cast<std::size_t>(i)];
    for (int i = 0; i < 5; ++i)
        out.T2[static_cast<std::size_t>(i)] = scale * t.T2[static_cast<std::size_t>(i)];
    return out;
}

}  // namespace h5reader::physics
