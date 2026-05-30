// h5reader_scene_math_tests — QtTest binary for the plane-frame math
// extracted from MoleculeScene.
//
// Covers:
//   computePlaneFrame on canonical 3-atom triangles → orthonormal basis
//   computePlaneFrame degenerate inputs (collinear, coincident, zero edge)
//   chooseContinuousNormal flips sign when dot(prev, current) < 0
//   chooseContinuousNormal leaves the input alone otherwise
//
// These are pure functions over h5reader::model::Vec3; no VTK, no Qt
// widget code. The binary builds against Qt6::Test + Eigen3 only.

#include "app/FitTargetMath.h"
#include "app/PlaneFrameMath.h"
#include "app/TensorGlyphMath.h"

#include <QObject>
#include <QtTest>

#include <array>
#include <cmath>
#include <vector>

using namespace h5reader;

// QFETCH/addColumn for model::Vec3 needs registration as a metatype.
// Local-only: header-side Q_DECLARE_METATYPE would force every consumer
// of QtAtom.h to know about it.
Q_DECLARE_METATYPE(h5reader::model::Vec3)

namespace {

bool nearly(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}

bool isUnitVector(const model::Vec3& v, double tol = 1e-9) {
    return nearly(v.norm(), 1.0, tol);
}

bool isPerpendicular(const model::Vec3& a, const model::Vec3& b, double tol = 1e-9) {
    return std::abs(a.dot(b)) <= tol;
}

}  // namespace

class SceneMathTests : public QObject {
    Q_OBJECT

private slots:
    // computePlaneFrame: canonical XY-plane triangle, normal points +z.
    void testComputePlaneFrame_canonicalXyTriangle();

    // computePlaneFrame: tilted triangle in arbitrary orientation —
    // basis is orthonormal, origin is the centroid.
    void testComputePlaneFrame_tiltedTriangle();

    // computePlaneFrame: x is along (b - a), y completes the right-handed
    // frame in the plane.
    void testComputePlaneFrame_xAxisIsAlongFirstEdge();

    // computePlaneFrame: degenerate inputs return nullopt.
    void testComputePlaneFrame_collinearReturnsNullopt();
    void testComputePlaneFrame_coincidentReturnsNullopt();
    void testComputePlaneFrame_zeroFirstEdgeReturnsNullopt();

    // chooseContinuousNormal: data-driven over the four interesting cases.
    void testChooseContinuousNormal_data();
    void testChooseContinuousNormal();

    // chooseContinuousNormal: a near-orthogonal current vector
    // (dot ≈ 0) is treated as same-side and returned unchanged.
    void testChooseContinuousNormal_nearOrthogonalLeavesAlone();

    // L-3a tensor-glyph eigendecomposition. Validates against a
    // diagonal tensor (axes = identity, radii = sqrt(diag)) and an
    // off-diagonal case with known eigenvalues.
    void testTensorGlyph_diagonalAxisAligned();
    void testTensorGlyph_offDiagonalKnownEigenvalues();
    // NOW-4 (2026-05-29): composeEllipsoidTransform body→lab rotation.
    // Body primary eigenvector along +z; current bond along +x →
    // primary axis ends up along +x in the world transform.
    void testTensorGlyph_composeAlignsPrimaryWithBondDir();
    void testTensorGlyph_composeAntiparallelHandled();

    // FitTargetMath — per-mode fit anchors. Pure functions over Vec3;
    // tested without a renderer or loaded conformation. Each function
    // is one variant of CameraComposer::write's dispatch table.
    void testFitTarget_atomAnchorReturnsPosition();
    void testFitTarget_bondAnchorMidpointAndAxis();
    void testFitTarget_bondAnchorZeroLengthReturnsNullopt();
    void testFitTarget_dihedralAnchorTransPlanar();
    void testFitTarget_dihedralAnchorZeroCentralBondReturnsNullopt();
    void testFitTarget_planeAnchorWrapsComputePlaneFrame();
    void testFitTarget_subsetKabschIdentity();
    void testFitTarget_subsetKabschKnownRotation();
    void testFitTarget_subsetKabschNlessThan3Nullopt();
    void testFitTarget_orthogonalizeViewUp();

    // writeSubset's rotation half — verify that R^T applied to the
    // captured camera-relative vector recovers the camera's pose in
    // the rotated frame. This is the math the CameraComposer's
    // writeSubset uses to keep the camera following the molecule's
    // rotation; tested as a pure invariant on ComputeSubsetTransform's
    // output so we don't drag in vtkRenderer or QtProtein.
    void testFitTarget_subsetCameraFollowsRotation();
    void testFitTarget_subsetCameraIdentityIsIdentity();

    // Codex finding #4 — rank-degenerate Kabsch must freeze
    // (return nullopt) rather than ship a numerically-arbitrary rotation
    // from SVD's null-space.
    void testFitTarget_subsetKabschCollinearReturnsNullopt();
    void testFitTarget_subsetKabschCoplanarReturnsNullopt();
    void testFitTarget_subsetKabschWellConditionedAccepted();

    // Codex finding #3 — safeViewUp guarantees non-degenerate output.
    // Trivial cases pass through; pathological (parallel) inputs fall
    // back deterministically.
    void testFitTarget_safeViewUpPassThroughNonDegenerate();
    void testFitTarget_safeViewUpDegenerateFallsBackToWorldZ();
    void testFitTarget_safeViewUpSightAlongZFallsBackToWorldX();
    void testFitTarget_safeViewUpDeterministic();

    // Codex finding #1 — dihedral sight-axis sign continuity uses an
    // EXPLICIT stored reference, not implicit feedback through the live
    // camera. Verify that an input axis crossing through perpendicular
    // to the stored reference produces a continuous post-flip output
    // sequence (no 180° flip on the boundary frame).
    void testFitTarget_dihedralContinuityNoFlipAtBoundary();
};

// ---- computePlaneFrame: canonical XY triangle ---------------------------

void SceneMathTests::testComputePlaneFrame_canonicalXyTriangle() {
    const std::array<model::Vec3, 3> positions{{
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(2.0, 0.0, 0.0),
        model::Vec3(0.0, 2.0, 0.0),
    }};

    const auto frame = math::computePlaneFrame(positions);
    QVERIFY(frame.has_value());

    // Origin at the centroid.
    QVERIFY(nearly(frame->origin.x(), 2.0 / 3.0));
    QVERIFY(nearly(frame->origin.y(), 2.0 / 3.0));
    QVERIFY(nearly(frame->origin.z(), 0.0));

    // x is along (b - a) normalised — here (1, 0, 0).
    QVERIFY(nearly(frame->x.x(), 1.0));
    QVERIFY(nearly(frame->x.y(), 0.0));
    QVERIFY(nearly(frame->x.z(), 0.0));

    // z is the cross product normalised — here (0, 0, 1).
    QVERIFY(nearly(frame->z.x(), 0.0));
    QVERIFY(nearly(frame->z.y(), 0.0));
    QVERIFY(nearly(frame->z.z(), 1.0));

    // y = z × x — here (0, 1, 0).
    QVERIFY(nearly(frame->y.x(), 0.0));
    QVERIFY(nearly(frame->y.y(), 1.0));
    QVERIFY(nearly(frame->y.z(), 0.0));
}

// ---- computePlaneFrame: tilted triangle ---------------------------------

void SceneMathTests::testComputePlaneFrame_tiltedTriangle() {
    const std::array<model::Vec3, 3> positions{{
        model::Vec3(1.0,  2.0, 3.0),
        model::Vec3(3.0,  2.0, 4.0),
        model::Vec3(1.5,  5.0, 3.5),
    }};

    const auto frame = math::computePlaneFrame(positions);
    QVERIFY(frame.has_value());

    QVERIFY(isUnitVector(frame->x));
    QVERIFY(isUnitVector(frame->y));
    QVERIFY(isUnitVector(frame->z));
    QVERIFY(isPerpendicular(frame->x, frame->y));
    QVERIFY(isPerpendicular(frame->y, frame->z));
    QVERIFY(isPerpendicular(frame->z, frame->x));

    // Right-handed: x × y = z.
    const model::Vec3 cross = frame->x.cross(frame->y);
    QVERIFY(nearly(cross.x(), frame->z.x(), 1e-12));
    QVERIFY(nearly(cross.y(), frame->z.y(), 1e-12));
    QVERIFY(nearly(cross.z(), frame->z.z(), 1e-12));

    // Origin is the centroid.
    QVERIFY(nearly(frame->origin.x(), (1.0 + 3.0 + 1.5) / 3.0));
    QVERIFY(nearly(frame->origin.y(), (2.0 + 2.0 + 5.0) / 3.0));
    QVERIFY(nearly(frame->origin.z(), (3.0 + 4.0 + 3.5) / 3.0));
}

// ---- computePlaneFrame: x axis is (b - a) normalised --------------------

void SceneMathTests::testComputePlaneFrame_xAxisIsAlongFirstEdge() {
    const std::array<model::Vec3, 3> positions{{
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(3.0, 4.0, 0.0),     // |b-a| = 5
        model::Vec3(0.0, 0.0, 1.0),
    }};

    const auto frame = math::computePlaneFrame(positions);
    QVERIFY(frame.has_value());
    QVERIFY(nearly(frame->x.x(), 0.6));
    QVERIFY(nearly(frame->x.y(), 0.8));
    QVERIFY(nearly(frame->x.z(), 0.0));
}

// ---- computePlaneFrame: degenerate inputs ------------------------------

void SceneMathTests::testComputePlaneFrame_collinearReturnsNullopt() {
    const std::array<model::Vec3, 3> positions{{
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(2.0, 0.0, 0.0),
    }};
    QVERIFY(!math::computePlaneFrame(positions).has_value());
}

void SceneMathTests::testComputePlaneFrame_coincidentReturnsNullopt() {
    const std::array<model::Vec3, 3> positions{{
        model::Vec3(1.0, 2.0, 3.0),
        model::Vec3(1.0, 2.0, 3.0),     // coincides with [0]
        model::Vec3(5.0, 6.0, 7.0),
    }};
    QVERIFY(!math::computePlaneFrame(positions).has_value());
}

void SceneMathTests::testComputePlaneFrame_zeroFirstEdgeReturnsNullopt() {
    // a == b → (b-a) is zero, x cannot be defined.
    const std::array<model::Vec3, 3> positions{{
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(1.0, 0.0, 0.0),
    }};
    QVERIFY(!math::computePlaneFrame(positions).has_value());
}

// ---- chooseContinuousNormal -------------------------------------------

void SceneMathTests::testChooseContinuousNormal_data() {
    QTest::addColumn<model::Vec3>("previous");
    QTest::addColumn<model::Vec3>("current");
    QTest::addColumn<model::Vec3>("expected");

    QTest::newRow("identical")    << model::Vec3(0, 0, 1) << model::Vec3(0, 0, 1) << model::Vec3(0, 0, 1);
    QTest::newRow("opposite")     << model::Vec3(0, 0, 1) << model::Vec3(0, 0, -1) << model::Vec3(0, 0, 1);
    QTest::newRow("near-aligned") << model::Vec3(0, 0, 1) << model::Vec3(0.1, 0.1, 0.99)
                                  << model::Vec3(0.1, 0.1, 0.99);
    QTest::newRow("near-flipped") << model::Vec3(0, 0, 1) << model::Vec3(0.1, 0.1, -0.99)
                                  << model::Vec3(-0.1, -0.1, 0.99);
    QTest::newRow("tilted-axis")  << model::Vec3(1.0, 0.5, 0.2) << model::Vec3(-0.9, -0.4, -0.1)
                                  << model::Vec3(0.9, 0.4, 0.1);
}

void SceneMathTests::testChooseContinuousNormal() {
    QFETCH(model::Vec3, previous);
    QFETCH(model::Vec3, current);
    QFETCH(model::Vec3, expected);
    const model::Vec3 result = math::chooseContinuousNormal(previous, current);
    QVERIFY(nearly(result.x(), expected.x(), 1e-12));
    QVERIFY(nearly(result.y(), expected.y(), 1e-12));
    QVERIFY(nearly(result.z(), expected.z(), 1e-12));
}

void SceneMathTests::testChooseContinuousNormal_nearOrthogonalLeavesAlone() {
    // dot ≈ 0 — by convention "not strictly negative", so the current
    // vector is returned unchanged. This is the contract: only a
    // *negative* dot triggers the flip.
    const model::Vec3 previous(0, 0, 1);
    const model::Vec3 current(1, 0, 1e-15);
    const model::Vec3 result = math::chooseContinuousNormal(previous, current);
    QVERIFY(nearly(result.x(), current.x(), 1e-12));
    QVERIFY(nearly(result.y(), current.y(), 1e-12));
    QVERIFY(nearly(result.z(), current.z(), 1e-12));
}

// ---- TensorGlyphMath: diagonal tensor ----------------------------------

void SceneMathTests::testTensorGlyph_diagonalAxisAligned() {
    // Diagonal {0.7, 0.2, 0.1} — eigenvectors are the standard basis;
    // eigenvalues sorted descending = {0.7, 0.2, 0.1}; radii =
    // sqrt(eigenvalues).
    const std::array<double, 9> diag = {
        0.7, 0.0, 0.0,
        0.0, 0.2, 0.0,
        0.0, 0.0, 0.1,
    };
    const auto e = math::decomposeSymmetric3x3(diag);

    QVERIFY(nearly(e.radii[0], std::sqrt(0.7), 1e-9));
    QVERIFY(nearly(e.radii[1], std::sqrt(0.2), 1e-9));
    QVERIFY(nearly(e.radii[2], std::sqrt(0.1), 1e-9));

    // First axis aligned with +x or -x (eigenvector sign is arbitrary).
    QVERIFY(nearly(std::abs(e.axes[0].x()), 1.0, 1e-9));
    QVERIFY(nearly(std::abs(e.axes[1].y()), 1.0, 1e-9));
    QVERIFY(nearly(std::abs(e.axes[2].z()), 1.0, 1e-9));

    // Orthonormal frame.
    QVERIFY(isUnitVector(e.axes[0]));
    QVERIFY(isUnitVector(e.axes[1]));
    QVERIFY(isUnitVector(e.axes[2]));
    QVERIFY(isPerpendicular(e.axes[0], e.axes[1]));
    QVERIFY(isPerpendicular(e.axes[1], e.axes[2]));
    QVERIFY(isPerpendicular(e.axes[2], e.axes[0]));
}

// ---- TensorGlyphMath: off-diagonal with known eigenvalues --------------

void SceneMathTests::testTensorGlyph_offDiagonalKnownEigenvalues() {
    // Symmetric 2x2 block {{4, 1}, {1, 3}} embedded in z=2 → analytic
    // eigenvalues are 2 (for z block) and (7 ± √5) / 2 = {4.618, 2.382}.
    const std::array<double, 9> M = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 0.0,
        0.0, 0.0, 2.0,
    };
    const auto e = math::decomposeSymmetric3x3(M);

    // 2x2 xy block {{4,1},{1,3}} has eigenvalues (7 ± √5) / 2; the z
    // block contributes a third eigenvalue of 2. Sorted descending:
    // 4.618 (xy) > 2.382 (xy) > 2.000 (z) — so the z eigenvector
    // lands at axes[2], NOT axes[1].
    const double lambda_xy_max = (7.0 + std::sqrt(5.0)) / 2.0;  // ≈ 4.618
    const double lambda_xy_min = (7.0 - std::sqrt(5.0)) / 2.0;  // ≈ 2.382
    const double lambda_z      = 2.0;

    QVERIFY(nearly(e.radii[0], std::sqrt(lambda_xy_max), 1e-9));
    QVERIFY(nearly(e.radii[1], std::sqrt(lambda_xy_min), 1e-9));
    QVERIFY(nearly(e.radii[2], std::sqrt(lambda_z),      1e-9));

    QVERIFY(isUnitVector(e.axes[0]));
    QVERIFY(isUnitVector(e.axes[1]));
    QVERIFY(isUnitVector(e.axes[2]));
    QVERIFY(isPerpendicular(e.axes[0], e.axes[1]));
    QVERIFY(isPerpendicular(e.axes[1], e.axes[2]));
    QVERIFY(isPerpendicular(e.axes[2], e.axes[0]));

    // The z-block eigenvector is the smallest (axes[2]); xy-block
    // eigenvectors (axes[0] and axes[1]) lie in the xy plane.
    QVERIFY(nearly(std::abs(e.axes[2].z()), 1.0, 1e-9));
    QVERIFY(nearly(e.axes[0].z(), 0.0, 1e-9));
    QVERIFY(nearly(e.axes[1].z(), 0.0, 1e-9));
}

void SceneMathTests::testTensorGlyph_composeAlignsPrimaryWithBondDir() {
    // Body-frame eigenvectors: primary along +z, secondaries +y, +x.
    // Radii: 2, 1, 0.5. Scale 1.0. Midpoint (5, 5, 5). Bond along +x.
    math::TensorEllipsoid eig;
    eig.axes = {model::Vec3(0, 0, 1), model::Vec3(0, 1, 0), model::Vec3(1, 0, 0)};
    eig.radii = {2.0, 1.0, 0.5};

    const model::Vec3 bondDir(1, 0, 0);
    const model::Vec3 mid(5.0, 5.0, 5.0);
    const auto M = math::composeEllipsoidTransform(eig, bondDir, mid, 1.0);

    // Column 0 (primary axis scaled by radii[0]=2) should now point
    // along +x in the world transform, magnitude 2.
    QVERIFY(nearly(M[0 * 4 + 0], 2.0, 1e-9));
    QVERIFY(nearly(M[1 * 4 + 0], 0.0, 1e-9));
    QVERIFY(nearly(M[2 * 4 + 0], 0.0, 1e-9));

    // Translation column = midpoint.
    QVERIFY(nearly(M[0 * 4 + 3], 5.0, 1e-12));
    QVERIFY(nearly(M[1 * 4 + 3], 5.0, 1e-12));
    QVERIFY(nearly(M[2 * 4 + 3], 5.0, 1e-12));
    QVERIFY(nearly(M[3 * 4 + 3], 1.0, 1e-12));

    // Homogeneous row.
    QVERIFY(nearly(M[3 * 4 + 0], 0.0, 1e-12));
    QVERIFY(nearly(M[3 * 4 + 1], 0.0, 1e-12));
    QVERIFY(nearly(M[3 * 4 + 2], 0.0, 1e-12));

    // Lengths of secondary columns = radii[1], radii[2].
    const double s1 = std::sqrt(M[0 * 4 + 1] * M[0 * 4 + 1] +
                                M[1 * 4 + 1] * M[1 * 4 + 1] +
                                M[2 * 4 + 1] * M[2 * 4 + 1]);
    const double s2 = std::sqrt(M[0 * 4 + 2] * M[0 * 4 + 2] +
                                M[1 * 4 + 2] * M[1 * 4 + 2] +
                                M[2 * 4 + 2] * M[2 * 4 + 2]);
    QVERIFY(nearly(s1, 1.0, 1e-9));
    QVERIFY(nearly(s2, 0.5, 1e-9));
}

void SceneMathTests::testTensorGlyph_composeAntiparallelHandled() {
    // Body primary along +z, bond along -z (antiparallel) — historically
    // the FromTwoVectors degenerate case. Should still produce a valid
    // rotation matrix (180° around any axis perpendicular to +z).
    math::TensorEllipsoid eig;
    eig.axes = {model::Vec3(0, 0, 1), model::Vec3(0, 1, 0), model::Vec3(1, 0, 0)};
    eig.radii = {1.0, 1.0, 1.0};

    const model::Vec3 bondDir(0, 0, -1);
    const model::Vec3 mid(0, 0, 0);
    const auto M = math::composeEllipsoidTransform(eig, bondDir, mid, 1.0);

    // Column 0 (primary, scaled by radii[0]=1) should point along -z.
    QVERIFY(nearly(M[0 * 4 + 0], 0.0, 1e-9));
    QVERIFY(nearly(M[1 * 4 + 0], 0.0, 1e-9));
    QVERIFY(nearly(M[2 * 4 + 0], -1.0, 1e-9));

    // Translation column = origin.
    QVERIFY(nearly(M[0 * 4 + 3], 0.0));
    QVERIFY(nearly(M[1 * 4 + 3], 0.0));
    QVERIFY(nearly(M[2 * 4 + 3], 0.0));
}

// ---- FitTargetMath: per-mode anchor functions --------------------------

void SceneMathTests::testFitTarget_atomAnchorReturnsPosition() {
    const std::array<model::Vec3, 1> positions{{model::Vec3(1.0, 2.0, 3.0)}};
    const auto anchor = math::ComputeAtomAnchor(positions);
    QVERIFY(anchor.has_value());
    QVERIFY(nearly(anchor->focal.x(), 1.0));
    QVERIFY(nearly(anchor->focal.y(), 2.0));
    QVERIFY(nearly(anchor->focal.z(), 3.0));
    QVERIFY(!anchor->axis.has_value());
    QVERIFY(!anchor->frame.has_value());
}

void SceneMathTests::testFitTarget_bondAnchorMidpointAndAxis() {
    // Bond along +x of length 2, atoms at (1,1,1) and (3,1,1).
    // Expected: focal = (2,1,1); axis = (1,0,0).
    const std::array<model::Vec3, 2> positions{{
        model::Vec3(1.0, 1.0, 1.0),
        model::Vec3(3.0, 1.0, 1.0),
    }};
    const auto anchor = math::ComputeBondAnchor(positions);
    QVERIFY(anchor.has_value());
    QVERIFY(nearly(anchor->focal.x(), 2.0));
    QVERIFY(nearly(anchor->focal.y(), 1.0));
    QVERIFY(nearly(anchor->focal.z(), 1.0));
    QVERIFY(anchor->axis.has_value());
    QVERIFY(nearly(anchor->axis->x(), 1.0));
    QVERIFY(nearly(anchor->axis->y(), 0.0));
    QVERIFY(nearly(anchor->axis->z(), 0.0));
}

void SceneMathTests::testFitTarget_bondAnchorZeroLengthReturnsNullopt() {
    const std::array<model::Vec3, 2> positions{{
        model::Vec3(1.0, 2.0, 3.0),
        model::Vec3(1.0, 2.0, 3.0),
    }};
    QVERIFY(!math::ComputeBondAnchor(positions).has_value());
}

void SceneMathTests::testFitTarget_dihedralAnchorTransPlanar() {
    // Canonical trans dihedral on the XY plane. (a) and (d) on opposite
    // sides of the central (b->c) axis along +z. Axis = (c-b)/||(c-b)||
    // = +z. viewUp = (a-b) projected perp to +z = (a-b) restricted to
    // XY plane, then normalised. With a = (-1, 0, 0), b = (0, 0, 0),
    // c = (0, 0, 1), d = (1, 0, 1), viewUp = -x (because (a-b)=(-1,0,0)).
    const std::array<model::Vec3, 4> positions{{
        model::Vec3(-1.0, 0.0, 0.0),
        model::Vec3( 0.0, 0.0, 0.0),
        model::Vec3( 0.0, 0.0, 1.0),
        model::Vec3( 1.0, 0.0, 1.0),
    }};
    const auto anchor = math::ComputeDihedralAnchor(positions);
    QVERIFY(anchor.has_value());
    // Focal = midpoint of (b, c) = (0, 0, 0.5).
    QVERIFY(nearly(anchor->focal.x(), 0.0));
    QVERIFY(nearly(anchor->focal.y(), 0.0));
    QVERIFY(nearly(anchor->focal.z(), 0.5));
    // Axis = +z.
    QVERIFY(anchor->axis.has_value());
    QVERIFY(nearly(anchor->axis->x(), 0.0));
    QVERIFY(nearly(anchor->axis->y(), 0.0));
    QVERIFY(nearly(anchor->axis->z(), 1.0));
    // viewUp = -x (normalised).
    QVERIFY(anchor->viewUp.has_value());
    QVERIFY(nearly(anchor->viewUp->x(), -1.0));
    QVERIFY(nearly(anchor->viewUp->y(), 0.0));
    QVERIFY(nearly(anchor->viewUp->z(), 0.0));
}

void SceneMathTests::testFitTarget_dihedralAnchorZeroCentralBondReturnsNullopt() {
    // b == c → central axis is zero-length.
    const std::array<model::Vec3, 4> positions{{
        model::Vec3(1.0, 2.0, 3.0),
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(4.0, 5.0, 6.0),
    }};
    QVERIFY(!math::ComputeDihedralAnchor(positions).has_value());
}

void SceneMathTests::testFitTarget_planeAnchorWrapsComputePlaneFrame() {
    const std::array<model::Vec3, 3> positions{{
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(2.0, 0.0, 0.0),
        model::Vec3(0.0, 2.0, 0.0),
    }};
    const auto anchor = math::ComputePlaneAnchor(positions);
    QVERIFY(anchor.has_value());
    QVERIFY(anchor->axis.has_value());
    QVERIFY(anchor->frame.has_value());
    // Focal = centroid (2/3, 2/3, 0).
    QVERIFY(nearly(anchor->focal.x(), 2.0 / 3.0));
    QVERIFY(nearly(anchor->focal.y(), 2.0 / 3.0));
    QVERIFY(nearly(anchor->focal.z(), 0.0));
    // Axis = +z (plane normal).
    QVERIFY(nearly(anchor->axis->z(), 1.0));
}

void SceneMathTests::testFitTarget_subsetKabschIdentity() {
    // current == reference → R = I, T = 0.
    std::vector<model::Vec3> ref = {
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(0.0, 1.0, 0.0),
        model::Vec3(0.0, 0.0, 1.0),
    };
    std::vector<model::Vec3> cur = ref;
    const auto t = math::ComputeSubsetTransform(cur, ref);
    QVERIFY(t.has_value());
    // R == I, T == 0.
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            QVERIFY(nearly(t->R(i, j), (i == j) ? 1.0 : 0.0, 1e-9));
        }
        QVERIFY(nearly(t->T(i), 0.0, 1e-9));
    }
}

void SceneMathTests::testFitTarget_subsetKabschKnownRotation() {
    // current = reference rotated +90° about z (about origin). Kabsch
    // recovers R = rotation by -90° about z so that R * current ≈
    // reference.
    //
    // Reference is a non-coplanar tetrahedron-ish shape so the
    // Codex-finding-#4 rank check (sigma[2] > tol) accepts it; a flat
    // square in the XY plane would correctly be rejected as
    // rank-deficient. The rotation math here is identical to the planar
    // case — the +z atom just keeps sigma[2] well-defined.
    std::vector<model::Vec3> ref = {
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(0.0, 1.0, 0.0),
        model::Vec3(-1.0, 0.0, 0.0),
        model::Vec3(0.0, -1.0, 0.0),
        model::Vec3(0.0, 0.0, 1.0),   // out of XY plane → full rank
    };
    // Rotate ref by +90° about z to get current.
    std::vector<model::Vec3> cur;
    cur.reserve(ref.size());
    for (const auto& p : ref) {
        cur.emplace_back(-p.y(), p.x(), p.z());
    }
    const auto t = math::ComputeSubsetTransform(cur, ref);
    QVERIFY(t.has_value());
    // R should rotate +90° about z (cur->ref means cur->ref = +90)
    // Actually R * cur ≈ ref, where cur is ref rotated +90 about z, so
    // R must rotate -90 about z, i.e. R = [[0, 1, 0], [-1, 0, 0],
    // [0, 0, 1]].
    QVERIFY(nearly(t->R(0, 0), 0.0, 1e-9));
    QVERIFY(nearly(t->R(0, 1), 1.0, 1e-9));
    QVERIFY(nearly(t->R(1, 0), -1.0, 1e-9));
    QVERIFY(nearly(t->R(1, 1), 0.0, 1e-9));
    QVERIFY(nearly(t->R(2, 2), 1.0, 1e-9));
}

void SceneMathTests::testFitTarget_subsetKabschNlessThan3Nullopt() {
    std::vector<model::Vec3> ref = {model::Vec3(0, 0, 0), model::Vec3(1, 0, 0)};
    std::vector<model::Vec3> cur = ref;
    QVERIFY(!math::ComputeSubsetTransform(cur, ref).has_value());
}

void SceneMathTests::testFitTarget_orthogonalizeViewUp() {
    // viewDir along +z; viewUp candidate (0.3, 0.4, 0.5) — orthogonal
    // version drops the +z component, normalises (0.3, 0.4, 0).
    const model::Vec3 viewDir(0.0, 0.0, 1.0);
    const model::Vec3 candidate(0.3, 0.4, 0.5);
    const auto orthog = math::OrthogonalizeViewUp(viewDir, candidate);
    QVERIFY(orthog.has_value());
    QVERIFY(nearly(orthog->z(), 0.0, 1e-9));
    QVERIFY(nearly(orthog->norm(), 1.0, 1e-9));
    // Parallel input returns nullopt.
    const model::Vec3 parallel(0.0, 0.0, 5.0);
    QVERIFY(!math::OrthogonalizeViewUp(viewDir, parallel).has_value());
}

// ---- writeSubset rotation half --------------------------------------------
//
// The Kabsch fit returns R, T such that R * current[i] + T approximates
// reference[i]. The CameraComposer's writeSubset uses R^T as the
// "body-to-current" rotation so that the camera, captured in the
// reference frame, follows the molecule into its current orientation.
//
// Two invariants tested here:
//   1. Identity input → R^T is identity, camera-relative vector
//      unchanged.
//   2. Rotated input → R^T applied to a captured camera-relative
//      vector recovers a vector that, added to the current centroid,
//      gives the same camera position relative to the rotated
//      subset that the captured camera had relative to the reference
//      subset.

void SceneMathTests::testFitTarget_subsetCameraIdentityIsIdentity() {
    // Reference and current both equal — Kabsch returns identity.
    // Camera captured at (0, 0, 10) looking at centroid (0, 0, 0),
    // up = (0, 1, 0). Rotated by R^T = I, the camera-relative vector
    // is unchanged.
    //
    // Non-coplanar tetrahedron-style ref so Codex finding #4's rank
    // check accepts it; the rotation math is identical to a planar set
    // but the +z atom keeps sigma[2] well-conditioned.
    std::vector<model::Vec3> ref = {
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(0.0, 1.0, 0.0),
        model::Vec3(-1.0, 0.0, 0.0),
        model::Vec3(0.0, 0.0, 1.0),
    };
    std::vector<model::Vec3> cur = ref;
    const auto t = math::ComputeSubsetTransform(cur, ref);
    QVERIFY(t.has_value());

    const model::Mat3 Rinv = t->R.transpose();
    const model::Vec3 capturedCamRel(0.0, 0.0, 10.0);   // pos - centroid
    const model::Vec3 capturedUp(0.0, 1.0, 0.0);

    const model::Vec3 newCamRel = Rinv * capturedCamRel;
    const model::Vec3 newUp     = Rinv * capturedUp;

    QVERIFY(nearly(newCamRel.x(), 0.0, 1e-9));
    QVERIFY(nearly(newCamRel.y(), 0.0, 1e-9));
    QVERIFY(nearly(newCamRel.z(), 10.0, 1e-9));
    QVERIFY(nearly(newUp.x(), 0.0, 1e-9));
    QVERIFY(nearly(newUp.y(), 1.0, 1e-9));
    QVERIFY(nearly(newUp.z(), 0.0, 1e-9));
}

void SceneMathTests::testFitTarget_subsetCameraFollowsRotation() {
    // Reference subset on the XY plane, centroid at origin. Apply a
    // known +90° rotation about z to every atom to get the current
    // subset. The Kabsch fit returns R such that R * current + T = ref;
    // that R is a -90° rotation about z (it rotates current back to
    // ref). R^T is then a +90° rotation about z — which is exactly
    // the rotation the molecule underwent.
    //
    // Captured camera at (5, 0, 0) looking at centroid (0, 0, 0) with
    // up = (0, 0, 1). The captured camera-relative vector is (5, 0, 0)
    // (camera was 5 units along +x from the centroid).
    //
    // Apply R^T (+90° about z): (5, 0, 0) -> (0, 5, 0). Adding to the
    // current centroid (still 0,0,0 because the rotation preserves it
    // for an origin-centred subset) gives the new camera position
    // (0, 5, 0) — which is exactly where the camera "should" be to
    // see the rotated subset from the same relative angle.
    // Non-coplanar ref (octahedron — square plus apices at +z and -z)
    // so Codex finding #4's rank check accepts it. The +90°-about-z
    // rotation preserves the centroid AT THE ORIGIN (symmetric apices
    // cancel) and the math is identical to a planar square; the apices
    // just keep sigma[2] non-zero.
    std::vector<model::Vec3> ref = {
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(0.0, 1.0, 0.0),
        model::Vec3(-1.0, 0.0, 0.0),
        model::Vec3(0.0, -1.0, 0.0),
        model::Vec3(0.0, 0.0, 1.0),
        model::Vec3(0.0, 0.0, -1.0),
    };
    // Rotate ref +90° about z to make current.
    std::vector<model::Vec3> cur;
    cur.reserve(ref.size());
    for (const auto& p : ref) {
        cur.emplace_back(-p.y(), p.x(), p.z());
    }
    const auto t = math::ComputeSubsetTransform(cur, ref);
    QVERIFY(t.has_value());

    const model::Mat3 Rinv = t->R.transpose();
    const model::Vec3 capturedCamRel(5.0, 0.0, 0.0);    // along +x
    const model::Vec3 capturedUp(0.0, 0.0, 1.0);

    const model::Vec3 newCamRel = Rinv * capturedCamRel;
    const model::Vec3 newUp     = Rinv * capturedUp;

    // R^T applied to (5, 0, 0) — for a +90° rotation about z — gives
    // (0, 5, 0): the captured-camera direction rotated by the same
    // angle the subset rotated.
    QVERIFY(nearly(newCamRel.x(), 0.0, 1e-9));
    QVERIFY(nearly(newCamRel.y(), 5.0, 1e-9));
    QVERIFY(nearly(newCamRel.z(), 0.0, 1e-9));

    // Up vector unchanged: rotation was about z, up was along z.
    QVERIFY(nearly(newUp.x(), 0.0, 1e-9));
    QVERIFY(nearly(newUp.y(), 0.0, 1e-9));
    QVERIFY(nearly(newUp.z(), 1.0, 1e-9));

    // The new sight direction (focal - position) should be in the same
    // body-relative direction (always pointing toward the subset
    // centroid from the rotated camera position).
    model::Vec3 newCentroid = model::Vec3::Zero();
    for (const auto& p : cur) newCentroid += p;
    newCentroid /= static_cast<double>(cur.size());
    QVERIFY(nearly(newCentroid.x(), 0.0, 1e-9));
    QVERIFY(nearly(newCentroid.y(), 0.0, 1e-9));
    QVERIFY(nearly(newCentroid.z(), 0.0, 1e-9));

    const model::Vec3 newPos    = newCentroid + newCamRel;
    const model::Vec3 newSight  = (newCentroid - newPos).normalized();
    // Captured sight was (centroid - pos).normalized() = (-5,0,0)/5 = -x.
    // R^T applied to -x is -y (the rotated equivalent).
    QVERIFY(nearly(newSight.x(), 0.0, 1e-9));
    QVERIFY(nearly(newSight.y(), -1.0, 1e-9));
    QVERIFY(nearly(newSight.z(), 0.0, 1e-9));
}

// ---- Codex finding #4: rank-degenerate Kabsch -------------------------

void SceneMathTests::testFitTarget_subsetKabschCollinearReturnsNullopt() {
    // 4 atoms strictly on the x-axis — rank 1 (only one non-zero
    // singular value). Rotation around the x-axis is underdetermined;
    // the function must freeze (return nullopt) rather than ship a
    // numerically-arbitrary rotation from SVD's null space.
    std::vector<model::Vec3> ref = {
        model::Vec3(0.0, 0.0, 0.0),
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(2.0, 0.0, 0.0),
        model::Vec3(3.0, 0.0, 0.0),
    };
    // Translate AND rotate slightly: even rotated, the line is rank-1,
    // so the fit is still underdetermined about the line's axis.
    std::vector<model::Vec3> cur;
    cur.reserve(ref.size());
    for (const auto& p : ref) {
        // Translate by (10, 5, 2), rotate around z by a tiny angle so
        // it's not strictly the same axis but still collinear.
        const double theta = 0.0;  // pure translation; result is still 1D
        const double x = std::cos(theta) * p.x() - std::sin(theta) * p.y();
        const double y = std::sin(theta) * p.x() + std::cos(theta) * p.y();
        cur.emplace_back(x + 10.0, y + 5.0, p.z() + 2.0);
    }
    QVERIFY(!math::ComputeSubsetTransform(cur, ref).has_value());
}

void SceneMathTests::testFitTarget_subsetKabschCoplanarReturnsNullopt() {
    // 4 atoms in the XY plane — rank 2. Rotation about the plane
    // normal is underdetermined; freeze.
    std::vector<model::Vec3> ref = {
        model::Vec3( 1.0,  0.0, 0.0),
        model::Vec3( 0.0,  1.0, 0.0),
        model::Vec3(-1.0,  0.0, 0.0),
        model::Vec3( 0.0, -1.0, 0.0),
    };
    std::vector<model::Vec3> cur = ref;  // identity input is fine, but...
    // ...rank is still 2 because all atoms have z=0. Identity input
    // would otherwise produce R=I, T=0; the rank check catches it
    // before the SVD does.
    QVERIFY(!math::ComputeSubsetTransform(cur, ref).has_value());
}

void SceneMathTests::testFitTarget_subsetKabschWellConditionedAccepted() {
    // 4 atoms forming a non-degenerate tetrahedron — rank 3, all
    // singular values well above threshold. Fit accepted.
    std::vector<model::Vec3> ref = {
        model::Vec3(1.0, 0.0, 0.0),
        model::Vec3(0.0, 1.0, 0.0),
        model::Vec3(0.0, 0.0, 1.0),
        model::Vec3(-1.0, -1.0, -1.0),
    };
    std::vector<model::Vec3> cur = ref;
    const auto t = math::ComputeSubsetTransform(cur, ref);
    QVERIFY(t.has_value());
    // Identity input -> R = I, det = +1, R^T*R = I.
    QVERIFY(nearly(t->R.determinant(), 1.0, 1e-9));
    const model::Mat3 RtR = t->R.transpose() * t->R;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            QVERIFY(nearly(RtR(i, j), (i == j) ? 1.0 : 0.0, 1e-9));
}

// ---- Codex finding #3: safeViewUp deterministic perpendicular fallback

void SceneMathTests::testFitTarget_safeViewUpPassThroughNonDegenerate() {
    // sight along +z; preferred along +y. Already perpendicular; the
    // result should equal preferred (normalised).
    const model::Vec3 sight(0.0, 0.0, 1.0);
    const model::Vec3 preferred(0.0, 1.0, 0.0);
    const model::Vec3 up = math::safeViewUp(sight, preferred);
    QVERIFY(nearly(up.x(), 0.0, 1e-9));
    QVERIFY(nearly(up.y(), 1.0, 1e-9));
    QVERIFY(nearly(up.z(), 0.0, 1e-9));
    QVERIFY(isUnitVector(up));
    QVERIFY(isPerpendicular(up, sight));
}

void SceneMathTests::testFitTarget_safeViewUpDegenerateFallsBackToWorldZ() {
    // sight along +y; preferred ALSO along +y (parallel). Old code
    // would have returned (0, 1, 0) which is parallel to sight (bad).
    // safeViewUp's fallback chain tries world Z first → (0, 0, 1).
    const model::Vec3 sight(0.0, 1.0, 0.0);
    const model::Vec3 preferred(0.0, 1.0, 0.0);
    const model::Vec3 up = math::safeViewUp(sight, preferred);
    QVERIFY(isUnitVector(up));
    QVERIFY(isPerpendicular(up, sight, 1e-9));
    // Specifically world Z because it's first in the fallback chain
    // and is perpendicular to sight=+y.
    QVERIFY(nearly(up.x(), 0.0, 1e-9));
    QVERIFY(nearly(up.y(), 0.0, 1e-9));
    QVERIFY(nearly(up.z(), 1.0, 1e-9));
}

void SceneMathTests::testFitTarget_safeViewUpSightAlongZFallsBackToWorldX() {
    // sight along +z; preferred ALSO along +z. World Z fallback would
    // also be parallel — skip to world X = (1, 0, 0).
    const model::Vec3 sight(0.0, 0.0, 1.0);
    const model::Vec3 preferred(0.0, 0.0, 1.0);
    const model::Vec3 up = math::safeViewUp(sight, preferred);
    QVERIFY(isUnitVector(up));
    QVERIFY(isPerpendicular(up, sight, 1e-9));
    // World X is the second fallback.
    QVERIFY(nearly(up.x(), 1.0, 1e-9));
    QVERIFY(nearly(up.y(), 0.0, 1e-9));
    QVERIFY(nearly(up.z(), 0.0, 1e-9));
}

void SceneMathTests::testFitTarget_safeViewUpDeterministic() {
    // Same inputs -> same output, every call. The deterministic
    // contract is what lets future debuggers reproduce a camera state
    // when chasing a regression.
    const model::Vec3 sight(0.0, 1.0, 0.0);
    const model::Vec3 preferred(0.0, 1.0, 0.0);
    const model::Vec3 up1 = math::safeViewUp(sight, preferred);
    const model::Vec3 up2 = math::safeViewUp(sight, preferred);
    QVERIFY(nearly(up1.x(), up2.x(), 0.0));
    QVERIFY(nearly(up1.y(), up2.y(), 0.0));
    QVERIFY(nearly(up1.z(), up2.z(), 0.0));
}

// ---- Codex finding #1: dihedral sight-axis sign continuity ------------

void SceneMathTests::testFitTarget_dihedralContinuityNoFlipAtBoundary() {
    // Simulate the dihedralLastDirection_ guard math (lifted from
    // CameraComposer::writeDihedral). The guard is: if the new axis
    // dots negative against the stored reference, flip its sign; then
    // update the reference. Without the guard, a 180° flip happens
    // when the natural cross-product crosses through the perpendicular
    // boundary.
    //
    // Scenario: a sequence of axes that drift through a configuration
    // where the natural axis flips sign. The guard should produce a
    // continuous post-flip sequence (each adjacent dot positive).
    std::optional<model::Vec3> lastDir;
    std::vector<model::Vec3> rawAxisSequence = {
        model::Vec3( 0.0,  0.0,  1.0),       // first frame
        model::Vec3( 0.1,  0.0,  0.99).normalized(),
        model::Vec3( 0.5,  0.0,  0.87).normalized(),  // near perpendicular
        model::Vec3(-0.5,  0.0, -0.87).normalized(),  // natural flip here
        model::Vec3(-0.1,  0.0, -0.99).normalized(),
    };
    std::vector<model::Vec3> postGuard;
    for (const auto& axis : rawAxisSequence) {
        model::Vec3 sightDir = axis;
        if (lastDir && sightDir.dot(*lastDir) < 0.0) {
            sightDir = -sightDir;
        }
        lastDir = sightDir;
        postGuard.push_back(sightDir);
    }
    // Adjacent frames must all dot positive — no 180° flips.
    for (std::size_t i = 1; i < postGuard.size(); ++i) {
        QVERIFY(postGuard[i].dot(postGuard[i - 1]) > 0.0);
    }
    // The "natural flip" frame (index 3) had its raw axis pointing
    // OPPOSITE the prior frame's stored direction; the guard must have
    // flipped its sign to maintain continuity.
    QVERIFY(rawAxisSequence[3].dot(postGuard[3]) < 0.0);
}

QTEST_GUILESS_MAIN(SceneMathTests)

#include "scene_math_tests.moc"
