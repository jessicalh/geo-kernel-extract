#pragma once

#include "LiteratureAccessors.h"

#include "../model/Types.h"

#include <cstddef>

namespace h5reader::physics {

inline double RingForwardContributionPpm(const model::SphericalTensor& literatureScaledKernel) {
    return literatureScaledKernel.T0;
}

double RingPerTypeT0Ppm(const double* unitCurrentT0ByType, std::size_t count);
model::SphericalTensor ScaleSphericalTensor(const model::SphericalTensor& t, double scale);

}  // namespace h5reader::physics
