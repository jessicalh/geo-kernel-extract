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
// (R * current[i] + T) and reference[i]. Lifted verbatim from the original
// TransformedConformation::KabschFit; both inputs must be the SAME
// length and the SAME atoms in the SAME order. Returns nullopt when
// the system is underdetermined (n < 3) or rank-deficient (collinear
// or coplanar subset; see "Rank degeneracy" below). The composer uses
// this to stabilise the camera against a subset of atoms (e.g. the
// backbone), matching what TransformedConformation::Mode::FitSubset
// does for positions but applied to the camera state instead.
// TransformedConformation::KabschFit delegates here (Codex finding #6)
// so both code paths share one degeneracy policy.
//
// Rank degeneracy (Codex finding #4):
//
// Kabsch is underdetermined when the input subset lies on a line
// (sigma[1] ~ 0 and sigma[2] ~ 0) or in a plane (sigma[2] ~ 0). In
// those cases Eigen's SVD returns SOME orthonormal null-space basis,
// but the choice is numerically arbitrary — the rotation about the
// degenerate axis can flip sign frame-to-frame even when the input
// moved only slightly. The visible failure is sudden 180° camera flips
// from a tiny molecular twitch.
//
// After SVD we inspect the singular values. If sigma[1] < kRankRelTol *
// sigma[0] OR sigma[2] < kRankRelTol * sigma[0], the fit is
// rank-deficient and the function returns std::nullopt. The relative
// threshold (default 1e-3) is scale-invariant — what matters is the
// shape of the subset, not its absolute size in Ångströms.
// Documentation of the chosen threshold lives at the constant.
//
// Determinant guard: even after the det-sign correction, accumulated
// floating-point error can leave det(R) slightly off from +1. We assert
// |det(R) - 1| < kDetTol; if the guard fires, we treat as rank-deficient
// and return std::nullopt rather than ship a bad rotation downstream
// (callers — CameraComposer::writeSubset, TransformedConformation —
// choose freeze-on-degenerate as their policy).
inline std::optional<Transform3D> ComputeSubsetTransform(
    const std::vector<model::Vec3>& current,
    const std::vector<model::Vec3>& reference) {
    // Relative singular-value threshold for the rank check. 1e-3 means
    // "the smallest covariance singular value must be at least 0.1 %
    // of the largest"; well-conditioned subsets (backbone of a normal
    // protein) sit at sigma[2]/sigma[0] ~ O(0.1-1.0). Anywhere near the
    // threshold is geometrically pathological (planar or near-linear
    // subset). 1e-3 chosen by default per the implementation prompt:
    // permissive enough to never reject the backbone, strict enough to
    // catch the documented failure modes (4-atom linear subset, 3-atom
    // collinear subset, planar 5-residue ring). Adjust upward if rank-
    // deficient false negatives are observed in the harness; downward
    // if backbone fits start freezing in practice.
    constexpr double kRankRelTol = 1e-3;
    // Determinant tolerance: SVD-derived rotations from well-conditioned
    // covariance matrices land within 1e-12; 1e-6 is many orders of
    // magnitude above noise and a clear "something went wrong" signal.
    constexpr double kDetTol = 1e-6;

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
    const Eigen::Vector3d sigma = svd.singularValues();

    // Rank-degenerate check: singular values are returned in descending
    // order by Eigen's JacobiSVD; check sigma[1] and sigma[2] against
    // sigma[0]. A zero sigma[0] means the centred coordinates were all
    // zero (every atom at the centroid) — definitionally rank-deficient.
    if (sigma(0) < kRankRelTol)
        return std::nullopt;
    if (sigma(1) < kRankRelTol * sigma(0))
        return std::nullopt;
    if (sigma(2) < kRankRelTol * sigma(0))
        return std::nullopt;

    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    const double d = (V * U.transpose()).determinant();
    D(2, 2) = (d < 0.0) ? -1.0 : 1.0;

    Transform3D out;
    out.R = V * D * U.transpose();
    out.T = cr - out.R * cc;

    // Belt-and-suspenders: the det-sign correction above should always
    // give det(R) = +1, but accumulated floating-point error during SVD
    // can leave it slightly off. Above kDetTol indicates SVD trouble
    // we'd rather not pipe to the camera.
    if (std::abs(out.R.determinant() - 1.0) > kDetTol)
        return std::nullopt;
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

// safeViewUp — guaranteed-non-degenerate ViewUp picker.
//
// Returns a unit vector strictly perpendicular to `sight`, derived from
// the caller's `preferred` direction projected perpendicular to sight.
// If `preferred` is parallel or anti-parallel to `sight` (cross product
// magnitude below tolerance), falls back deterministically by trying
// world Z, then world X, then world Y until one of them is non-degenerate
// against the sight axis. Never returns a vector parallel to sight; never
// returns std::nullopt; the world-axis fallback chain guarantees a hit
// because no single direction can be parallel to all three world axes
// simultaneously.
//
// Deterministic fallback policy (so a future debugger can reproduce):
//   1. project `preferred` perpendicular to sight; if norm above tol, use it
//   2. else try world Z = (0, 0, 1)
//   3. else try world X = (1, 0, 0)
//   4. else use world Y = (0, 1, 0) (always succeeds — sight can be
//      parallel to at most one of Z/X/Y, never all three)
//
// Fixes Codex finding #3 (dihedral sight-parallel ViewUp). Used by every
// write*() path that sets view-up; replaces ad-hoc fallback-to-(0,1,0)
// branches that left the camera with a degenerate up at sight directions
// aligned with world Y.
inline model::Vec3 safeViewUp(const model::Vec3& sight,
                                const model::Vec3& preferred) {
    constexpr double kPerpMin = 1e-9;

    // Normalise the sight axis; treat zero-sight as a programmer error
    // and pick world Y arbitrarily so we never crash on a bad caller.
    model::Vec3 sightUnit = sight;
    if (sightUnit.norm() < kPerpMin)
        return model::Vec3(0.0, 1.0, 0.0);
    sightUnit.normalize();

    // Try the preferred direction first.
    if (auto orthog = OrthogonalizeViewUp(sightUnit, preferred))
        return *orthog;

    // Deterministic fallback chain: world Z, world X, world Y in that
    // order. The first axis whose cross product with sight is above
    // tolerance wins. World Y is the last-resort because the caller's
    // ad-hoc fallback was (0,1,0) and we want different behaviour when
    // sight aligns with Y (the original failure mode).
    const std::array<model::Vec3, 3> fallbacks{{
        model::Vec3(0.0, 0.0, 1.0),
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(0.0, 1.0, 0.0),
    }};
    for (const auto& axis : fallbacks) {
        if (auto orthog = OrthogonalizeViewUp(sightUnit, axis))
            return *orthog;
    }
    // Mathematically unreachable: sight cannot be parallel to all three
    // world axes. Return world Y as a final safety so the function never
    // ill-conditions a downstream camera write.
    return model::Vec3(0.0, 1.0, 0.0);
}

}  // namespace h5reader::math
