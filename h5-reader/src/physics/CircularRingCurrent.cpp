#include "CircularRingCurrent.h"

#include "SphericalBasis.h"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>

namespace h5reader::physics {

namespace {

constexpr double kElementaryChargeC = 1.602176634e-19;
constexpr double kElectronMassKg = 9.1093837015e-31;
constexpr double kMu0TeslaMetrePerAmpere = 1.25663706212e-6;
constexpr double kAngstromToMetre = 1.0e-10;
constexpr double kNanoampereToAmpere = 1.0e-9;
constexpr double kPpm = 1.0e6;
constexpr double kPi = 3.141592653589793238462643383279502884;

constexpr double kBenzeneCurrentNanoamperePerTesla =
    -6.0 * kElementaryChargeC * kElementaryChargeC / (4.0 * kPi * kElectronMassKg) / kNanoampereToAmpere;

struct CylindricalField {
    double radial = 0.0;
    double axial = 0.0;
    bool valid = false;
};

CylindricalField singleLoop(double radiusA, double currentNanoamperePerTesla, double radialA, double axialA) {
    CylindricalField out;
    if (!(radiusA > 0.0) || !std::isfinite(radiusA) || !std::isfinite(currentNanoamperePerTesla) || !std::isfinite(radialA)
        || !std::isfinite(axialA)) {
        return out;
    }

    const double radius = radiusA * kAngstromToMetre;
    const double rho = std::abs(radialA) * kAngstromToMetre;
    const double z = axialA * kAngstromToMetre;
    const double current = currentNanoamperePerTesla * kNanoampereToAmpere;

    if (rho <= 1.0e-14 * radius) {
        out.axial = kMu0TeslaMetrePerAmpere * current * radius * radius / (2.0 * std::pow(radius * radius + z * z, 1.5));
        out.valid = std::isfinite(out.axial);
        return out;
    }

    const double alphaSq = (radius + rho) * (radius + rho) + z * z;
    const double betaSq = (radius - rho) * (radius - rho) + z * z;
    if (betaSq <= std::pow(1.0e-16 * radius, 2.0))
        return out;

    const double parameter = std::clamp(4.0 * radius * rho / alphaSq, 0.0, 1.0);
    const double modulus = std::sqrt(parameter);
    const double K = std::comp_ellint_1(modulus);
    const double E = std::comp_ellint_2(modulus);
    const double rootAlpha = std::sqrt(alphaSq);
    const double common = kMu0TeslaMetrePerAmpere * current / (2.0 * kPi * rootAlpha);

    out.radial = common * z / (rho * betaSq) * (-K * betaSq + E * (radius * radius + rho * rho + z * z));
    out.axial = common * (K + E * (radius * radius - rho * rho - z * z) / betaSq);
    out.valid = std::isfinite(out.radial) && std::isfinite(out.axial);
    return out;
}

}  // namespace

std::optional<CircularRingPlane> FitCircularRingPlane(const std::vector<model::Vec3>& vertices) {
    if (vertices.size() < 3)
        return std::nullopt;

    CircularRingPlane plane;
    for (const model::Vec3& vertex : vertices) {
        if (!vertex.allFinite())
            return std::nullopt;
        plane.center += vertex;
    }
    plane.center /= static_cast<double>(vertices.size());

    Eigen::MatrixXd centered(static_cast<Eigen::Index>(vertices.size()), 3);
    for (std::size_t i = 0; i < vertices.size(); ++i)
        centered.row(static_cast<Eigen::Index>(i)) = vertices[i] - plane.center;

    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinV);
    if (svd.matrixV().cols() < 3 || svd.singularValues().size() < 2 || !(svd.singularValues()(1) > 1.0e-12)) {
        return std::nullopt;
    }
    plane.normal = svd.matrixV().col(2);
    const double normalNorm = plane.normal.norm();
    if (!(normalNorm > 1.0e-12) || !std::isfinite(normalNorm))
        return std::nullopt;
    plane.normal /= normalNorm;

    double squaredOffsets = 0.0;
    for (const model::Vec3& vertex : vertices) {
        const double offset = (vertex - plane.center).dot(plane.normal);
        squaredOffsets += offset * offset;
    }
    plane.planeRmsA = std::sqrt(squaredOffsets / static_cast<double>(vertices.size()));
    return plane;
}

std::optional<CircularRingParameters> CandidateACircularParameters(model::RingTypeIndex type, int protonationVariant) {
    // Frozen circular Johnson--Bovey candidate used by the Reader provenance
    // figures. The radii, +/-0.64 A loop spacing, and Christensen-family
    // relative currents are the same declared table used by the independently
    // replayable Candidate-A analysis; they are not fitted in the Reader.
    double radiusA = 0.0;
    double relativeCurrent = 0.0;
    switch (type) {
    case model::RingTypeIndex::PheBenzene:
        radiusA = 1.39;
        relativeCurrent = 1.13;
        break;
    case model::RingTypeIndex::TyrPhenol:
        radiusA = 1.39;
        relativeCurrent = 0.91;
        break;
    case model::RingTypeIndex::TrpBenzene:
        radiusA = 1.39;
        relativeCurrent = 1.18;
        break;
    case model::RingTypeIndex::TrpPyrrole:
        radiusA = 1.182;
        relativeCurrent = 1.06;
        break;
    case model::RingTypeIndex::HisImidazole:
        if (protonationVariant != 2)
            return std::nullopt;
        radiusA = 1.182;
        relativeCurrent = 1.27;
        break;
    case model::RingTypeIndex::HidImidazole:
    case model::RingTypeIndex::HieImidazole:
        radiusA = 1.182;
        relativeCurrent = 1.25;
        break;
    case model::RingTypeIndex::TrpPerimeter:
    case model::RingTypeIndex::ProPyrrolidine:
        return std::nullopt;
    }

    return CircularRingParameters{
        radiusA,
        0.64,
        kBenzeneCurrentNanoamperePerTesla * relativeCurrent,
    };
}

std::optional<model::Vec3> EvaluateCircularDoubleLoopField(const model::Vec3& pointA,
                                                           const CircularRingPlane& plane,
                                                           const CircularRingParameters& parameters) {
    const model::Vec3& normalLike = plane.normal;
    const double normalNorm = normalLike.norm();
    if (!(normalNorm > 1.0e-12) || !std::isfinite(normalNorm) || !(parameters.radiusA > 0.0)
        || !std::isfinite(parameters.lobeOffsetA) || !std::isfinite(parameters.currentNanoamperePerTesla)) {
        return std::nullopt;
    }

    const model::Vec3 normal = normalLike / normalNorm;
    const model::Vec3 displacement = pointA - plane.center;
    const double zA = displacement.dot(normal);
    const model::Vec3 inPlane = displacement - zA * normal;
    const double rhoA = inPlane.norm();
    model::Vec3 rhoHat = model::Vec3::Zero();
    if (rhoA > 1.0e-14 * parameters.radiusA)
        rhoHat = inPlane / rhoA;

    const double halfCurrent = 0.5 * parameters.currentNanoamperePerTesla;
    const CylindricalField plus = singleLoop(parameters.radiusA, halfCurrent, rhoA, zA - parameters.lobeOffsetA);
    const CylindricalField minus = singleLoop(parameters.radiusA, halfCurrent, rhoA, zA + parameters.lobeOffsetA);
    if (!plus.valid || !minus.valid)
        return std::nullopt;

    return (plus.radial + minus.radial) * rhoHat + (plus.axial + minus.axial) * normal;
}

std::optional<model::SphericalTensor>
EvaluateCircularShielding(const model::Vec3& pointA, const CircularRingPlane& plane, const CircularRingParameters& parameters) {
    const model::Vec3& normalLike = plane.normal;
    const double normalNorm = normalLike.norm();
    if (!(normalNorm > 1.0e-12) || !std::isfinite(normalNorm))
        return std::nullopt;
    const model::Vec3 normal = normalLike / normalNorm;
    const CircularRingPlane normalizedPlane{plane.center, normal, plane.planeRmsA};
    const auto field = EvaluateCircularDoubleLoopField(pointA, normalizedPlane, parameters);
    if (!field)
        return std::nullopt;

    const model::Mat3 shielding = -kPpm * (*field) * normal.transpose();
    return DecomposeLibrary(shielding);
}

}  // namespace h5reader::physics
