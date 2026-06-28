#include "McConnellLiteratureKernel.h"

#include "LiteratureAccessors.h"
#include "SphericalBasis.h"
#include "constants/LiteratureConstants.h"

#include <cmath>

namespace h5reader::physics {
namespace {

bool finiteVec(const model::Vec3& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

}  // namespace

bool McConnellLiteratureCategory(model::BondCategory category) {
    switch (category) {
    case model::BondCategory::PeptideCO:
    case model::BondCategory::PeptideCN:
    case model::BondCategory::SidechainCO:
    case model::BondCategory::Aromatic:
        return true;
    default:
        return false;
    }
}

double McConnellDeltaChiQ(model::BondCategory category) {
    return McConnellDeltaChi(category).value;
}

double McConnellMolarPrefactor() {
    return nmr::constants::kMcConnellMolarPrefactor.value;
}

double McConnellForwardContributionPpm(const model::SphericalTensor& literatureScaledKernel) {
    return literatureScaledKernel.T0;
}

model::SphericalTensor McConnellSourceLiteratureKernelLocal(const McConnellSourceGeometry& source,
                                                            bool* present) {
    if (present) *present = false;
    model::SphericalTensor out;
    if (!(source.r > 1e-9) || !finiteVec(source.disp_local) || !finiteVec(source.bond_axis_local))
        return out;

    if (!McConnellLiteratureCategory(source.category)) return out;

    const double axisNorm = source.bond_axis_local.norm();
    if (!(axisNorm > 1e-9)) return out;

    const model::Vec3 dHat = -source.disp_local / source.r;
    const model::Vec3 bHat = source.bond_axis_local / axisNorm;
    const double cosTheta = dHat.dot(bHat);
    const double r3 = source.r * source.r * source.r;
    if (!(r3 > 0.0) || !std::isfinite(cosTheta)) return out;

    const model::Mat3 mCode =
        (9.0 * cosTheta * dHat * bHat.transpose()
         - 3.0 * bHat * bHat.transpose()
         - (3.0 * dHat * dHat.transpose() - model::Mat3::Identity()))
        / r3;

    const double q = McConnellDeltaChiQ(source.category);
    const double scale = -McConnellMolarPrefactor() * q / 3.0;
    const model::Mat3 sigma = scale * mCode;

    if (present) *present = true;
    return DecomposeLibrary(sigma);
}

}  // namespace h5reader::physics
