#include "BiotSavartRingCurrent.h"

#include "QtPhysicalConstants.h"

#include "../physics/SphericalBasis.h"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace h5reader::calculators {

namespace {

using model::Mat3;
using model::RingGeometry;
using model::Vec3;

// Finite straight-wire field, in SI units:
// B = (mu0/4pi) I (dl x dA)/|dl x dA|^2
//     (dl.dA/|dA| - dl.dB/|dB|).
Vec3 wireSegmentField(const Vec3& startMetres, const Vec3& endMetres,
                      double currentAmperes, const Vec3& pointMetres) {
    const Vec3 segment = endMetres - startMetres;
    const Vec3 fromStart = pointMetres - startMetres;
    const Vec3 fromEnd = pointMetres - endMetres;

    const double startDistance = fromStart.norm();
    const double endDistance = fromEnd.norm();
    if (startDistance < BS_WIRE_ENDPOINT_GUARD || endDistance < BS_WIRE_ENDPOINT_GUARD) {
        return Vec3::Zero();
    }

    const Vec3 cross = segment.cross(fromStart);
    const double crossSquared = cross.squaredNorm();
    if (crossSquared < BS_WIRE_AXIS_GUARD)
        return Vec3::Zero();

    const double scale = BIOT_SAVART_PREFACTOR * currentAmperes / crossSquared;
    const double endpointProjection =
        segment.dot(fromStart) / startDistance
        - segment.dot(fromEnd) / endDistance;
    return scale * endpointProjection * cross;
}

double pointToSegmentDistance(const Vec3& point, const Vec3& start, const Vec3& end) {
    const Vec3 segment = end - start;
    const double lengthSquared = segment.squaredNorm();
    if (lengthSquared <= std::numeric_limits<double>::epsilon())
        return (point - start).norm();
    const double fraction =
        std::clamp((point - start).dot(segment) / lengthSquared, 0.0, 1.0);
    return (point - (start + fraction * segment)).norm();
}

double minimumWireDistance(const std::vector<Vec3>& vertices,
                           const Vec3& normal,
                           double lobeOffsetA,
                           const Vec3& pointA) {
    const Vec3 offset = normal * lobeOffsetA;
    double minimum = std::numeric_limits<double>::infinity();
    for (std::size_t edge = 0; edge < vertices.size(); ++edge) {
        const std::size_t next = (edge + 1) % vertices.size();
        minimum = std::min(
            minimum,
            pointToSegmentDistance(pointA, vertices[edge] + offset,
                                   vertices[next] + offset));
        minimum = std::min(
            minimum,
            pointToSegmentDistance(pointA, vertices[edge] - offset,
                                   vertices[next] - offset));
    }
    return minimum;
}

RingGeometry fitProducerRingGeometry(const std::vector<Vec3>& vertices) {
    RingGeometry geometry;
    for (const Vec3& vertex : vertices)
        geometry.center += vertex;
    geometry.center /= static_cast<double>(vertices.size());

    Eigen::MatrixXd centered(static_cast<Eigen::Index>(vertices.size()), 3);
    for (std::size_t i = 0; i < vertices.size(); ++i)
        centered.row(static_cast<Eigen::Index>(i)) =
            (vertices[i] - geometry.center).transpose();

    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeFullV);
    geometry.normal = svd.matrixV().col(2);
    const Vec3 firstEdge = vertices[1] - vertices[0];
    const Vec3 secondEdge = vertices[2] - vertices[0];
    if (geometry.normal.dot(firstEdge.cross(secondEdge)) < 0.0)
        geometry.normal = -geometry.normal;

    for (const Vec3& vertex : vertices)
        geometry.radius += (vertex - geometry.center).norm();
    geometry.radius /= static_cast<double>(vertices.size());
    return geometry;
}

std::vector<Vec3> offsetVertices(const std::vector<Vec3>& vertices, const Vec3& offset) {
    std::vector<Vec3> result;
    result.reserve(vertices.size());
    for (const Vec3& vertex : vertices)
        result.push_back(vertex + offset);
    return result;
}

}  // namespace

const char* BiotSavartSampleStatusName(BiotSavartSampleStatus status) {
    switch (status) {
    case BiotSavartSampleStatus::Evaluated:
        return "evaluated";
    case BiotSavartSampleStatus::OutsideSpatialCutoff:
        return "outside_spatial_cutoff";
    case BiotSavartSampleStatus::WireSingularity:
        return "wire_singularity";
    case BiotSavartSampleStatus::InvalidInput:
        return "invalid_input";
    }
    return "invalid_input";
}

std::optional<BiotSavartRingCurrent>
BiotSavartRingCurrent::Build(const std::vector<Vec3>& orderedVerticesA, double lobeOffsetA) {
    if (orderedVerticesA.size() < 3 || !std::isfinite(lobeOffsetA))
        return std::nullopt;
    for (const Vec3& vertex : orderedVerticesA) {
        if (!vertex.allFinite())
            return std::nullopt;
    }

    RingGeometry geometry = fitProducerRingGeometry(orderedVerticesA);
    if (!geometry.center.allFinite() || !geometry.normal.allFinite()
        || !std::isfinite(geometry.radius)) {
        return std::nullopt;
    }
    return BiotSavartRingCurrent(orderedVerticesA, geometry, lobeOffsetA);
}

BiotSavartRingCurrent::BiotSavartRingCurrent(std::vector<Vec3> verticesA,
                                             RingGeometry geometry,
                                             double lobeOffsetA)
    : verticesA_(std::move(verticesA))
    , geometry_(std::move(geometry))
    , lobeOffsetA_(lobeOffsetA) {}

BiotSavartRingSample BiotSavartRingCurrent::evaluate(const Vec3& pointA,
                                                     double ringCurrentIntensityNanoamperesPerTesla,
                                                     bool includeSegments) const {
    BiotSavartRingSample sample;
    sample.pointA = pointA;
    sample.ringCurrentIntensityNanoamperesPerTesla =
        ringCurrentIntensityNanoamperesPerTesla;
    if (!pointA.allFinite()
        || !std::isfinite(ringCurrentIntensityNanoamperesPerTesla))
        return sample;

    sample.distanceToCenterA = (pointA - geometry_.center).norm();
    if (sample.distanceToCenterA > RING_CURRENT_CUTOFF) {
        sample.status = BiotSavartSampleStatus::OutsideSpatialCutoff;
        return sample;
    }

    sample.minimumWireDistanceA =
        minimumWireDistance(verticesA_, geometry_.normal, lobeOffsetA_, pointA);
    const double wireDistanceGuardA = BS_WIRE_ENDPOINT_GUARD / ANGSTROMS_TO_METRES;
    if (sample.minimumWireDistanceA <= wireDistanceGuardA) {
        sample.status = BiotSavartSampleStatus::WireSingularity;
        return sample;
    }

    const Vec3 offsetA = geometry_.normal * lobeOffsetA_;
    // Evaluate the induced current at B0 = 1 T. The resulting field is the
    // secondary-field response in T/T.
    const double halfCurrentAmperes =
        0.5 * ringCurrentIntensityNanoamperesPerTesla
        * NANOAMPERES_TO_AMPERES;
    const Vec3 pointMetres = pointA * ANGSTROMS_TO_METRES;
    if (includeSegments)
        sample.segments.reserve(verticesA_.size() * 2);

    for (std::size_t edge = 0; edge < verticesA_.size(); ++edge) {
        const std::size_t next = (edge + 1) % verticesA_.size();
        for (const int loopSide : {1, -1}) {
            const Vec3 loopOffsetA = static_cast<double>(loopSide) * offsetA;
            const Vec3 startA = verticesA_[edge] + loopOffsetA;
            const Vec3 endA = verticesA_[next] + loopOffsetA;
            const Vec3 field = wireSegmentField(
                startA * ANGSTROMS_TO_METRES,
                endA * ANGSTROMS_TO_METRES,
                halfCurrentAmperes,
                pointMetres);
            sample.inducedFieldTeslaPerTesla += field;
            if (includeSegments) {
                sample.segments.push_back(
                    BiotSavartSegmentContribution{edge, loopSide, startA, endA, field});
            }
        }
    }

    // With B evaluated per unit applied field, sigma_ab = -B_a n_b 10^6.
    sample.shieldingTensorPpm =
        -sample.inducedFieldTeslaPerTesla
        * geometry_.normal.transpose()
        * PPM_FACTOR;
    sample.spherical = physics::DecomposeLibrary(sample.shieldingTensorPpm);
    sample.status = BiotSavartSampleStatus::Evaluated;
    return sample;
}

std::vector<Vec3> BiotSavartRingCurrent::upperLoopVerticesA() const {
    return offsetVertices(verticesA_, geometry_.normal * lobeOffsetA_);
}

std::vector<Vec3> BiotSavartRingCurrent::lowerLoopVerticesA() const {
    return offsetVertices(verticesA_, -geometry_.normal * lobeOffsetA_);
}

}  // namespace h5reader::calculators
