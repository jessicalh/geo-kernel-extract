// TensorGlyphMath — pure-function eigendecomposition for the
// bond-orientation Mat3 glyph rendered by SceneRevealOverlay.
//
// Decomposes a symmetric 3x3 tensor (the producer's
// bond_orientation_tensor ⟨u⊗u⟩, sampled per bond vector per
// trajectory) into:
//
//   - three principal axes (orthonormal Vec3, sorted by descending
//     eigenvalue magnitude),
//   - three radii (positive double, = sqrt(max(eigenvalue, eps))).
//
// The radii feed a vtkSphereSource scaled+rotated by a vtkTransform
// to draw the ellipsoid attached to the highlighted bond. Sorting by
// descending eigenvalue puts the dominant axis first; the glyph's
// "elongation" axis matches the bond vector for highly ordered
// (S² → 1) systems.
//
// Eigen3 is the only dependency; no VTK, no Qt — by design so this
// can be unit-tested headlessly in scene_math_tests.
//
// L-3a (2026-05-29).

#pragma once

#include "../model/QtAtom.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <array>
#include <cmath>

namespace h5reader::math {

struct TensorEllipsoid {
    // Principal axes — orthonormal, sorted by descending eigenvalue.
    std::array<model::Vec3, 3> axes{};
    // Radii along the principal axes; same sort order as `axes`.
    // Always positive — clamped to `eps_radius` at the floor so a
    // near-zero eigenvalue still produces a visible disc rather than a
    // degenerate line/dot. The clamp magnitude is tiny relative to a
    // typical unit-vector tensor whose eigenvalues sum to 1.
    std::array<double, 3> radii{};
};

// Decompose a row-major flat 3x3 tensor into principal axes + radii.
// `flat` indices: [0..2]=row0, [3..5]=row1, [6..8]=row2. The matrix is
// assumed symmetric; non-symmetric input is symmetrised as
// 0.5 * (M + M^T) before solving so a numeric perturbation in the
// producer doesn't poison the result.
//
// NOTE on coordinate frame: the producer's bond_orientation_tensor is
// computed in the BODY frame — that is, the trajectory is rigidly
// aligned to its frame-0 reference each frame, and ⟨u⊗u⟩ accumulates
// in that aligned frame. So the eigenvectors returned here are
// expressed in the body frame's basis, NOT the lab frame at the
// current playback frame. Use composeEllipsoidTransform below to
// render correctly in the current scene — passing the body-frame
// eigenvectors directly into a VTK world transform produces a glyph
// whose orientation is wrong once the molecule rotates relative to
// frame 0 (Codex NOW-4, 2026-05-29).
inline TensorEllipsoid decomposeSymmetric3x3(const std::array<double, 9>& flat,
                                             double eps_radius = 1e-6) {
    Eigen::Matrix3d M;
    M << flat[0], flat[1], flat[2],
         flat[3], flat[4], flat[5],
         flat[6], flat[7], flat[8];
    // Symmetrise — the bond_orientation_tensor IS symmetric by
    // construction (it's ⟨u⊗u⟩), but a producer float-to-double or
    // packing artefact could leave it ~ 1e-15 asymmetric. The
    // SelfAdjointEigenSolver assumes symmetry; feeding asymmetric
    // input is undefined.
    M = 0.5 * (M + M.transpose());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(M);
    Eigen::Vector3d eigenvalues = solver.eigenvalues();
    Eigen::Matrix3d eigenvectors = solver.eigenvectors();

    // SelfAdjointEigenSolver returns eigenvalues in ascending order;
    // we want descending so the largest axis is axes[0] (matches the
    // glyph rendering convention — dominant axis = first).
    std::array<int, 3> order = {2, 1, 0};

    TensorEllipsoid out;
    for (int i = 0; i < 3; ++i) {
        const int idx = order[i];
        const Eigen::Vector3d v = eigenvectors.col(idx);
        out.axes[i] = model::Vec3(v.x(), v.y(), v.z());
        const double lambda = eigenvalues(idx);
        // sqrt(eigenvalue) gives ellipsoid radii. Clamp at eps_radius
        // so a numerically-zero eigenvalue still produces a visible
        // (but small) sliver. Negative eigenvalues should not occur
        // for a PSD orientation tensor; if they do (numeric noise),
        // we floor them too.
        out.radii[i] = std::sqrt(std::max(lambda, eps_radius * eps_radius));
    }
    return out;
}

// Build a row-major 4×4 transform matrix that places an ellipsoid in
// the current scene: principal axis aligned with currentBondDir,
// secondary axes rotated consistently, eigenvalue-derived radii
// (sqrt eigenvalue) scaled by `scale`, translated to `midpoint`.
//
// This handles the body→lab coord-frame issue (Codex NOW-4): the body
// frame's z is approximately the reference bond direction; aligning
// the largest eigenvector with the CURRENT bond direction
// approximates the body→lab rotation without needing the full Kabsch
// rotation matrix per frame. The in-plane rotation around the bond
// axis is unconstrained by this approximation but the ellipsoid's
// secondary radii still scale correctly.
//
// Edge cases:
// - bondDir near-parallel to eigenvectors[0]: identity rotation.
// - bondDir near-antiparallel: 180° rotation around any perpendicular.
// - bondDir of zero length: caller must skip (returns identity here
//   for defensiveness but the ellipsoid would render at midpoint with
//   body-frame orientation — caller should test bondDir.norm() > 0
//   before invoking).
inline std::array<double, 16>
composeEllipsoidTransform(const TensorEllipsoid& eig,
                          const model::Vec3& currentBondDir,
                          const model::Vec3& midpoint,
                          double scale) {
    // Rotation R that maps eig.axes[0] (body-frame primary) to
    // currentBondDir. Using Eigen::Quaterniond::FromTwoVectors for
    // numerical stability (handles the antiparallel case correctly).
    const Eigen::Vector3d body0(eig.axes[0].x(), eig.axes[0].y(), eig.axes[0].z());
    const Eigen::Vector3d bondDirE(currentBondDir.x(), currentBondDir.y(), currentBondDir.z());
    const double bondNorm = bondDirE.norm();
    Eigen::Matrix3d R;
    if (bondNorm < 1e-9) {
        R.setIdentity();
    } else {
        const Eigen::Vector3d bondUnit = bondDirE / bondNorm;
        R = Eigen::Quaterniond::FromTwoVectors(body0, bondUnit).toRotationMatrix();
    }

    // Lab-frame eigenvectors after applying R.
    const Eigen::Vector3d body1(eig.axes[1].x(), eig.axes[1].y(), eig.axes[1].z());
    const Eigen::Vector3d body2(eig.axes[2].x(), eig.axes[2].y(), eig.axes[2].z());
    const Eigen::Vector3d lab0 = R * body0;  // == bondUnit (verification)
    const Eigen::Vector3d lab1 = R * body1;
    const Eigen::Vector3d lab2 = R * body2;

    const double r0 = eig.radii[0] * scale;
    const double r1 = eig.radii[1] * scale;
    const double r2 = eig.radii[2] * scale;

    std::array<double, 16> M{};
    // Row-major layout (matches vtkMatrix4x4::SetElement(row, col, v)).
    // Columns 0..2 are the scaled lab-frame eigenvectors; column 3 is
    // the translation; row 3 is the homogeneous (0 0 0 1).
    M[0]  = lab0.x() * r0; M[1]  = lab1.x() * r1; M[2]  = lab2.x() * r2; M[3]  = midpoint.x();
    M[4]  = lab0.y() * r0; M[5]  = lab1.y() * r1; M[6]  = lab2.y() * r2; M[7]  = midpoint.y();
    M[8]  = lab0.z() * r0; M[9]  = lab1.z() * r1; M[10] = lab2.z() * r2; M[11] = midpoint.z();
    M[12] = 0.0;           M[13] = 0.0;           M[14] = 0.0;           M[15] = 1.0;
    return M;
}

}  // namespace h5reader::math
