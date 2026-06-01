// RingCurrentKernel — source-level Johnson-Bovey kernel emission for rediscover.
//
// The analysis consumer must not reconstruct the J-B/Biot-Savart field in
// Python. This helper keeps the source-level double-loop physics in the C++
// spine and returns the tensor in the rediscovery/library spherical basis.

#pragma once

#include "AnalysisBody.h"
#include "LocalFrameBasis.h"

#include "../model/QtRing.h"
#include "../model/Types.h"

#include <cstddef>

namespace h5reader::rediscover {

model::SphericalTensor JohnsonBoveySourceUnitKernelLocal(const Body& body,
                                                         const LocalFrame& frame,
                                                         const model::Vec3& atomPos,
                                                         std::size_t ringIdx,
                                                         const model::QtRing& ring,
                                                         std::size_t h5Row);

model::SphericalTensor ScaleSphericalTensor(const model::SphericalTensor& t, double scale);

}  // namespace h5reader::rediscover
