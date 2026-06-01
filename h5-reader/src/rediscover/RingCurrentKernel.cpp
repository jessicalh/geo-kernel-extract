#include "RingCurrentKernel.h"

#include "SphericalBasis.h"

#include "../calculators/QtBiotSavartCalc.h"
#include "../calculators/QtPhysicalConstants.h"
#include "../model/ConformationGeometry.h"
#include "../model/Conformation.h"

#include <vector>

namespace h5reader::rediscover {

model::SphericalTensor JohnsonBoveySourceUnitKernelLocal(const Body& body,
                                                         const LocalFrame& frame,
                                                         const model::Vec3& atomPos,
                                                         std::size_t ringIdx,
                                                         const model::QtRing& ring,
                                                         std::size_t h5Row) {
    const model::Conformation& conf = *body.run.conformation;
    const model::RingGeometry& geo = body.idx.ringGeometry.at(ringIdx, h5Row);
    const std::vector<model::Vec3> vertices = model::RingVertices(conf, ringIdx, h5Row);
    const model::Vec3 bUnitCurrent =
        calculators::EvaluateBField(atomPos, geo, vertices, ring.JohnsonBoveyLobeOffset(), 1.0);

    model::Mat3 sigma = model::Mat3::Zero();
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            sigma(a, b) = -geo.normal(b) * bUnitCurrent(a) * calculators::PPM_FACTOR;
        }
    }

    const model::Mat3 local = frame.is_valid ? frame.TensorToLocal(sigma) : sigma;
    return DecomposeLibrary(local);
}

model::SphericalTensor ScaleSphericalTensor(const model::SphericalTensor& t, double scale) {
    model::SphericalTensor out;
    out.T0 = scale * t.T0;
    for (int i = 0; i < 3; ++i) out.T1[static_cast<std::size_t>(i)] = scale * t.T1[static_cast<std::size_t>(i)];
    for (int i = 0; i < 5; ++i) out.T2[static_cast<std::size_t>(i)] = scale * t.T2[static_cast<std::size_t>(i)];
    return out;
}

}  // namespace h5reader::rediscover
