#include "McConnellLiteratureKernel.h"

#include "LiteratureConstants.h"
#include "SphericalBasis.h"

#include <cmath>

namespace h5reader::rediscover {

namespace {
bool finiteVec(const Vec3& v) {
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
    return kMcConnellMolarPrefactor.value;
}

double McConnellForwardContributionPpm(const model::SphericalTensor& literatureScaledKernel) {
    return literatureScaledKernel.T0;
}

model::SphericalTensor McConnellSourceLiteratureKernelLocal(const SourceSlot& source,
                                                            bool* present) {
    if (present) *present = false;
    model::SphericalTensor out;
    if (source.kind != SourceKind::Bond || !(source.r > 1e-9)
        || !finiteVec(source.disp_local) || !finiteVec(source.bond_axis_local)) {
        return out;
    }

    const auto category = static_cast<model::BondCategory>(source.bond_category);
    if (!McConnellLiteratureCategory(category)) return out;

    const double axisNorm = source.bond_axis_local.norm();
    if (!(axisNorm > 1e-9)) return out;

    const Vec3 dHat = -source.disp_local / source.r;  // bond midpoint -> target atom
    const Vec3 bHat = source.bond_axis_local / axisNorm;
    const double cosTheta = dHat.dot(bHat);
    const double r3 = source.r * source.r * source.r;
    if (!(r3 > 0.0) || !std::isfinite(cosTheta)) return out;

    const Mat3 mCode =
        (9.0 * cosTheta * dHat * bHat.transpose()
         - 3.0 * bHat * bHat.transpose()
         - (3.0 * dHat * dHat.transpose() - Mat3::Identity()))
        / r3;

    // The project-canonical M tensor has T0=f=(3cos^2-1)/r^3. The literature
    // shielding scalar is sigma=-prefactor*q*f/3, so scale M by -prefactor*q/3.
    // DecomposeLibrary keeps that signed isotropic scalar in T0 while deriving
    // T2 from the traceless symmetric part.
    const double q = McConnellDeltaChiQ(category);
    const double scale = -McConnellMolarPrefactor() * q / 3.0;
    const Mat3 sigma = scale * mCode;

    if (present) *present = true;
    return DecomposeLibrary(sigma);
}

}  // namespace h5reader::rediscover
