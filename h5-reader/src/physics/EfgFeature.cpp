#include "EfgFeature.h"

#include "SphericalBasis.h"

#include <cmath>

namespace h5reader::physics {

std::array<double, 5> DecomposeEfgLibraryT2(const model::Mat3& efg) {
    return DecomposeLibrary(efg).T2;
}

model::Mat3 ReconstructLibraryT2(const std::array<double, 5>& t2) {
    return ReconstructLibraryT2Matrix(t2);
}

std::array<double, 5> RotateLibraryT2ToLocal(const LocalFrame& frame,
                                             const std::array<double, 5>& labT2) {
    if (!frame.is_valid) return {};
    const model::Mat3 local = frame.TensorToLocal(ReconstructLibraryT2(labT2));
    return DecomposeEfgLibraryT2(local);
}

double T2Magnitude(const std::array<double, 5>& t2) {
    double s = 0.0;
    for (double v : t2) s += v * v;
    return std::sqrt(s);
}

bool FiniteT2(const std::array<double, 5>& t2) {
    for (double v : t2)
        if (!std::isfinite(v)) return false;
    return true;
}

}  // namespace h5reader::physics
