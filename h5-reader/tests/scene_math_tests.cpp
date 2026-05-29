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

#include "app/PlaneFrameMath.h"
#include "app/TensorGlyphMath.h"

#include <QObject>
#include <QtTest>

#include <array>
#include <cmath>

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

QTEST_GUILESS_MAIN(SceneMathTests)

#include "scene_math_tests.moc"
