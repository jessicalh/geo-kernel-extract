// FitTargetMath — pure-function math for the camera-composer's per-frame
// fit anchors. No QObject, no VTK, no I/O — by design, so the math is
// unit-testable from `tests/scene_math_tests` with Qt6::Test + Eigen
// only. Modelled on PlaneFrameMath.h.
//
// One free function per CameraMode kind that needs more than trivial
// position lookup; each returns std::nullopt for degenerate inputs
// (zero-length axis, collinear atoms, underdetermined Kabsch). Callers
// (CameraComposer) decide how to handle degeneracy — typically fall back
// to CameraMode::Free for that frame and log a warning.
//
// Kabsch implementation lifted from
// TransformedConformation::KabschFit (TransformedConformation.cpp:228-276).
// Keeping the math in a header lets the camera composer reuse the same
// algorithm without depending on the decorator type; the decorator's
// static method now delegates to ComputeSubsetTransform.

#pragma once

#include "../model/QtAtom.h"
#include "PlaneFrameMath.h"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::math {

struct AtomAnchor {
    model::Vec3 focal = model::Vec3::Zero();   // F: where the camera looks
    // sight axis for Bond/Dihedral; normal for Plane; std::nullopt for Atom
    std::optional<model::Vec3> axis;
    // full orthonormal frame for Plane; std::nullopt otherwise
    std::optional<PlaneFrame> frame;
    // For Dihedral: the (a - b) direction projected perpendicular to axis,
    // useful as the natural ViewUp for a Newman projection. std::nullopt
    // for Atom / Bond / Plane.
    std::optional<model::Vec3> viewUp;
};

// Kabsch fit output: rotation + translation such that
// R * current[i] + T approximates reference[i].
struct Transform3D {
    model::Mat3 R = model::Mat3::Identity();
    model::Vec3 T = model::Vec3::Zero();
};

// Atom — F = positions[0]; axis = std::nullopt; frame = std::nullopt.
// Trivial, but typed-function shape so the composer dispatches uniformly.
inline std::optional<AtomAnchor> ComputeAtomAnchor(
    const std::array<model::Vec3, 1>& positions) {
    AtomAnchor a;
    a.focal = positions[0];
    return a;
}

// Bond — F = midpoint; axis = (b - a).normalized(); frame = std::nullopt.
// Returns std::nullopt for zero-length bonds (a == b).
inline std::optional<AtomAnchor> ComputeBondAnchor(
    const std::array<model::Vec3, 2>& positions) {
    constexpr double kAxisMin = 1e-9;
    const model::Vec3 axis = positions[1] - positions[0];
    if (axis.norm() < kAxisMin)
        return std::nullopt;
    AtomAnchor a;
    a.focal = 0.5 * (positions[0] + positions[1]);
    a.axis  = axis.normalized();
    return a;
}

// Dihedral — F = midpoint(b, c); axis = (c - b).normalized().
// viewUp = (a - b) projected perpendicular to axis, normalised — the
// natural Newman-projection orientation. Returns std::nullopt for a
// zero-length central bond. If (a - b) is parallel to axis, viewUp
// falls back to the projection of (d - c) perpendicular to axis; if
// THAT is also parallel, viewUp is std::nullopt and the composer keeps
// the gesture's current ViewUp.
inline std::optional<AtomAnchor> ComputeDihedralAnchor(
    const std::array<model::Vec3, 4>& positions) {
    constexpr double kAxisMin = 1e-9;
    const model::Vec3& a = positions[0];
    const model::Vec3& b = positions[1];
    const model::Vec3& c = positions[2];
    const model::Vec3& d = positions[3];

    const model::Vec3 rawAxis = c - b;
    if (rawAxis.norm() < kAxisMin)
        return std::nullopt;
    const model::Vec3 axis = rawAxis.normalized();

    AtomAnchor out;
    out.focal = 0.5 * (b + c);
    out.axis  = axis;

    // ViewUp from (a - b), projected perpendicular to axis. If degenerate,
    // try (d - c). Leave nullopt if both legs are parallel to the axis.
    auto project = [&](const model::Vec3& v) -> std::optional<model::Vec3> {
        const model::Vec3 perp = v - v.dot(axis) * axis;
        if (perp.norm() < kAxisMin)
            return std::nullopt;
        return perp.normalized();
    };
    if (auto up = project(a - b))
        out.viewUp = up;
    else if (auto up = project(d - c))
        out.viewUp = up;
    return out;
}

// Plane — thin wrapper around computePlaneFrame from PlaneFrameMath.h.
// The wrapping lets the composer route every kind through one helper
// family; the underlying frame still feeds the existing sign-continuity
// guard (chooseContinuousNormal).
inline std::optional<AtomAnchor> ComputePlaneAnchor(
    const std::array<model::Vec3, 3>& positions) {
    auto frame = computePlaneFrame(positions);
    if (!frame)
        return std::nullopt;
    AtomAnchor a;
    a.focal = frame->origin;
    a.axis  = frame->z;   // plane normal
    a.frame = *frame;
    return a;
}

// Subset Kabsch — minimises sum of squared distances between
// (R * current[i] + T) and reference[i]. Lifted verbatim from
// TransformedConformation::KabschFit; both inputs must be the SAME
// length and the SAME atoms in the SAME order. Returns nullopt when
// the system is underdetermined (n < 3). The composer uses this to
// stabilise the camera against a subset of atoms (e.g. the backbone),
// matching what TransformedConformation::Mode::FitSubset does for
// positions but applied to the camera state instead.
inline std::optional<Transform3D> ComputeSubsetTransform(
    const std::vector<model::Vec3>& current,
    const std::vector<model::Vec3>& reference) {
    const std::size_t n = current.size();
    if (n != reference.size() || n < 3)
        return std::nullopt;

    model::Vec3 cc = model::Vec3::Zero();
    model::Vec3 cr = model::Vec3::Zero();
    for (std::size_t i = 0; i < n; ++i) {
        cc += current[i];
        cr += reference[i];
    }
    cc /= static_cast<double>(n);
    cr /= static_cast<double>(n);

    Eigen::MatrixXd P(3, n);
    Eigen::MatrixXd Q(3, n);
    for (std::size_t i = 0; i < n; ++i) {
        P.col(static_cast<Eigen::Index>(i)) = current[i] - cc;
        Q.col(static_cast<Eigen::Index>(i)) = reference[i] - cr;
    }

    Eigen::Matrix3d H = P * Q.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::Matrix3d& U = svd.matrixU();
    const Eigen::Matrix3d& V = svd.matrixV();

    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    const double d = (V * U.transpose()).determinant();
    D(2, 2) = (d < 0.0) ? -1.0 : 1.0;

    Transform3D out;
    out.R = V * D * U.transpose();
    out.T = cr - out.R * cc;
    return out;
}

// View-up orthogonalisation — gives back a ViewUp made strictly
// perpendicular to a view direction. Mirrors vtkCamera::OrthogonalizeViewUp
// (one Gram-Schmidt step). Returns std::nullopt if the input ViewUp is
// parallel to the view direction (caller picks a fallback). Used by the
// composer when assembling the final absolute camera write.
inline std::optional<model::Vec3> OrthogonalizeViewUp(
    const model::Vec3& viewDir, const model::Vec3& candidateUp) {
    constexpr double kPerpMin = 1e-9;
    const model::Vec3 perp = candidateUp - candidateUp.dot(viewDir) * viewDir;
    if (perp.norm() < kPerpMin)
        return std::nullopt;
    return perp.normalized();
}

}  // namespace h5reader::math
