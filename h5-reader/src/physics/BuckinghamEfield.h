#pragma once

#include "LocalFrameBasis.h"

#include "../model/Types.h"

#include <string_view>

namespace h5reader::physics {

struct BuckinghamFieldProjection {
    bool present = false;
    model::Vec3 efield_local = model::Vec3::Zero();
    double e_parallel = 0.0;
    double e_magnitude = 0.0;
    double e_magnitude_squared = 0.0;
};

BuckinghamFieldProjection ProjectBuckinghamEfield(const LocalFrame& frame,
                                                  const model::Vec3& efield_lab);
double BuckinghamShiftPpm(double e_parallel, double a, double b = 0.0);
double BuckinghamShiftPpm(model::Element element,
                          std::string_view residueType,
                          std::string_view atomName,
                          std::string_view frameKind,
                          double e_parallel);

}  // namespace h5reader::physics
