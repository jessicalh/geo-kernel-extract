#include "BuckinghamEfield.h"

#include "LiteratureAccessors.h"

#include <cmath>
#include <limits>

namespace h5reader::physics {

BuckinghamFieldProjection ProjectBuckinghamEfield(const LocalFrame& frame,
                                                  const model::Vec3& efield_lab) {
    BuckinghamFieldProjection out;
    if (!frame.is_valid || !std::isfinite(efield_lab.x()) || !std::isfinite(efield_lab.y())
        || !std::isfinite(efield_lab.z()))
        return out;

    out.present = true;
    out.efield_local = frame.ToLocal(efield_lab);
    out.e_parallel = out.efield_local.z();
    out.e_magnitude_squared = efield_lab.squaredNorm();
    out.e_magnitude = std::sqrt(out.e_magnitude_squared);
    return out;
}

double BuckinghamShiftPpm(double e_parallel, double a, double b) {
    if (!std::isfinite(e_parallel) || !std::isfinite(a) || !std::isfinite(b))
        return std::numeric_limits<double>::quiet_NaN();
    return -a * e_parallel - b * e_parallel * e_parallel;
}

double BuckinghamShiftPpm(model::Element element,
                          std::string_view residueType,
                          std::string_view atomName,
                          std::string_view frameKind,
                          double e_parallel) {
    const double a = BuckinghamA(element, residueType, atomName, frameKind).value;
    const double b = BuckinghamB(element, residueType, atomName, frameKind).value;
    return BuckinghamShiftPpm(e_parallel, a, b);
}

}  // namespace h5reader::physics
