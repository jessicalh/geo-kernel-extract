#include "McConnellLiteratureKernel.h"

#include "SphericalBasis.h"

#include <cmath>

namespace h5reader::rediscover {

namespace {

constexpr double kAvogadro = 6.02214076e23;

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
    switch (category) {
    case model::BondCategory::PeptideCO:
        return 2.41;
    case model::BondCategory::PeptideCN:
        return -5.42;
    case model::BondCategory::SidechainCO:
        return 2.41;
    case model::BondCategory::Aromatic:
        return 0.0;  // RING carries the aromatic pi current; avoid double-counting.
    default:
        return 0.0;
    }
}

double McConnellMolarPrefactor() {
    return 1.0e24 / kAvogadro;
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
    // shielding scalar is sigma=-prefactor*q*f/3, so scale M by -prefactor*q/3
    // before taking the traceless PCS tensor. The trace removal leaves T2
    // unchanged and makes mc_lit_T0 an explicit near-zero audit channel.
    const double q = McConnellDeltaChiQ(category);
    const double scale = -McConnellMolarPrefactor() * q / 3.0;
    Mat3 sigma = scale * mCode;
    sigma -= (sigma.trace() / 3.0) * Mat3::Identity();

    if (present) *present = true;
    return DecomposeLibrary(sigma);
}

}  // namespace h5reader::rediscover
