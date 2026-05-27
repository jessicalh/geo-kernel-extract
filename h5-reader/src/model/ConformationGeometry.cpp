// ConformationGeometry — implementation. Lifted from QtFrame's former
// private FitRingGeometry + ringVertices so single pose and trajectory
// share one interpretation of ring geometry.

#include "ConformationGeometry.h"

#include "Conformation.h"
#include "QtProtein.h"
#include "AtomSelection.h"  // GeometryKind enumerators (None/Distance/Angle/Dihedral)

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace h5reader::model {

std::vector<Vec3> RingVertices(const Conformation& conf, std::size_t ringIdx, std::size_t frame) {
    std::vector<Vec3> verts;
    const QtProtein* p = conf.protein();
    if (!p || ringIdx >= p->ringCount())
        return verts;
    const QtRing& ring = p->ring(ringIdx);
    verts.reserve(ring.atomIndices.size());
    for (int32_t a : ring.atomIndices)
        verts.push_back(conf.atomPosition(frame, static_cast<std::size_t>(a)));
    return verts;
}

RingGeometry FitRingGeometry(const std::vector<Vec3>& verts) {
    RingGeometry g;
    if (verts.empty())
        return g;

    // Center = centroid.
    g.center = Vec3::Zero();
    for (const auto& v : verts)
        g.center += v;
    g.center /= static_cast<double>(verts.size());

    // Normal = smallest singular vector of the centered point matrix.
    if (verts.size() >= 3) {
        Eigen::Matrix<double, Eigen::Dynamic, 3> M(static_cast<Eigen::Index>(verts.size()), 3);
        for (std::size_t i = 0; i < verts.size(); ++i)
            M.row(static_cast<Eigen::Index>(i)) = (verts[i] - g.center).transpose();
        Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, 3>> svd(M, Eigen::ComputeFullV);
        g.normal = svd.matrixV().col(2);
        if (g.normal.norm() > 1e-12)
            g.normal.normalize();
    }

    // Radius = mean distance from center.
    double r_sum = 0.0;
    for (const auto& v : verts)
        r_sum += (v - g.center).norm();
    g.radius = r_sum / static_cast<double>(verts.size());
    return g;
}

RingGeometry RingGeometryAt(const Conformation& conf, std::size_t ringIdx, std::size_t frame) {
    const QtProtein* p = conf.protein();
    if (!p || ringIdx >= p->ringCount())
        return {};
    return FitRingGeometry(RingVertices(conf, ringIdx, frame));
}

// --- Geometry of an ordered atom selection -------------------------------

namespace {
constexpr double kRadToDeg = 57.295779513082320876798;  // 180 / pi
// Below this (Å) a difference vector has no reliable direction, so the angle
// / dihedral it would define is undefined rather than noisy.
constexpr double kMinVecNorm = 1e-9;
}  // namespace

double Distance(const Vec3& a, const Vec3& b) {
    return (b - a).norm();
}

double AngleDegrees(const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3   u  = a - b;  // vertex at the middle atom b
    const Vec3   v  = c - b;
    const double nu = u.norm();
    const double nv = v.norm();
    if (nu < kMinVecNorm || nv < kMinVecNorm)
        return std::numeric_limits<double>::quiet_NaN();
    // Clamp guards floating-point drift just outside acos's [-1, 1] domain.
    const double cosang = std::clamp(u.dot(v) / (nu * nv), -1.0, 1.0);
    return std::acos(cosang) * kRadToDeg;  // [0, 180]
}

double DihedralDegrees(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    // Signed dihedral about the b-c axis (Blondel-Karplus atan2 convention).
    const Vec3   b1  = b - a;
    const Vec3   b2  = c - b;
    const Vec3   b3  = d - c;
    const Vec3   n1  = b1.cross(b2);  // normal of plane (a, b, c)
    const Vec3   n2  = b2.cross(b3);  // normal of plane (b, c, d)
    const double b2n = b2.norm();
    if (n1.norm() < kMinVecNorm || n2.norm() < kMinVecNorm || b2n < kMinVecNorm)
        return std::numeric_limits<double>::quiet_NaN();
    // m completes a right-handed frame {n1, m, b2_hat}; (x = n1.n2, y = m.n2)
    // makes atan2 return the signed angle in (-180, 180].
    const Vec3   b2hat = b2 / b2n;
    const Vec3   m     = n1.cross(b2hat);
    const double x     = n1.dot(n2);
    const double y     = m.dot(n2);
    return std::atan2(y, x) * kRadToDeg;  // (-180, 180]
}

GeometryMeasurement Measure(const Conformation& conf, std::size_t frame,
                            const std::vector<std::size_t>& atoms) {
    GeometryMeasurement r;  // {None, 0.0, false}
    const QtProtein*    p = conf.protein();
    if (!p || frame >= conf.frameCount())
        return r;
    for (std::size_t idx : atoms)
        if (idx >= p->atomCount())
            return r;  // a stale/out-of-range index makes the whole tuple invalid

    const auto P = [&](std::size_t i) { return conf.atomPosition(frame, atoms[i]); };

    switch (atoms.size()) {
    case 2:
        r.kind  = GeometryKind::Distance;
        r.value = Distance(P(0), P(1));
        break;
    case 3:
        r.kind  = GeometryKind::Angle;
        r.value = AngleDegrees(P(0), P(1), P(2));
        break;
    case 4:
        r.kind  = GeometryKind::Dihedral;
        r.value = DihedralDegrees(P(0), P(1), P(2), P(3));
        break;
    default:
        return r;  // 0, 1, or >4 atoms define no measurement
    }
    r.valid = std::isfinite(r.value);
    return r;
}

}  // namespace h5reader::model
