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

QTEST_GUILESS_MAIN(SceneMathTests)

#include "scene_math_tests.moc"
