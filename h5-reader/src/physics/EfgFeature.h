#pragma once

#include "LocalFrameBasis.h"

#include "../model/Types.h"

#include <array>

namespace h5reader::physics {

std::array<double, 5> DecomposeEfgLibraryT2(const model::Mat3& efg);
model::Mat3 ReconstructLibraryT2(const std::array<double, 5>& t2);
std::array<double, 5> RotateLibraryT2ToLocal(const LocalFrame& frame,
                                             const std::array<double, 5>& labT2);
double T2Magnitude(const std::array<double, 5>& t2);
bool FiniteT2(const std::array<double, 5>& t2);

}  // namespace h5reader::physics
