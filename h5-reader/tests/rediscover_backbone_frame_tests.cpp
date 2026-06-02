// rediscover_backbone_frame_tests — pin the broad-backbone local frames
// (BuildBackboneNFrame / Ca / CarbonylC / CarbonylO / Ha) against KNOWN
// synthetic geometry. Pure Vec3 math: no calcset, no H5, no model — runs
// standalone under ctest, the SphericalBasis √6 fixture as precedent.
//
// The backbone frames are the genuinely-new pure-geometric piece of
// broad-backbone, and they are the exact bug class that produced the
// aromatic-H positional-anchor time bomb (wrong axis / handedness / wrong
// anchor). So the convention is pinned by ASSERTIONS, not by hoping a broad
// run looks right.
//
// For each class: axes unit-length + mutually orthogonal; right-handed
// (x × y ≈ z); z along the DEFINED convention for that class; the x azimuth
// points toward the typed anchor direction (its in-plane projection); and
// frame_valid == true. Plus a degenerate case (coincident/collinear anchors)
// → frame_valid == false with NO NaN poisoning.

#include "rediscover/BroadBackboneSink.h"
#include "rediscover/LocalFrameBasis.h"
#include "rediscover/SphericalBasis.h"

#include <QtTest>

#include <cmath>
#include <limits>
#include <vector>

using h5reader::model::Mat3;
using h5reader::model::SphericalTensor;
using h5reader::model::Vec3;
using h5reader::rediscover::BroadAggregate;
using h5reader::rediscover::BuildBackboneCaFrame;
using h5reader::rediscover::BuildBackboneCarbonylCFrame;
using h5reader::rediscover::BuildBackboneCarbonylOFrame;
using h5reader::rediscover::BuildBackboneHaFrame;
using h5reader::rediscover::BuildBackboneNFrame;
using h5reader::rediscover::DecomposeLibrary;
using h5reader::rediscover::FrameVariant;
using h5reader::rediscover::LocalFrame;
using h5reader::rediscover::ReduceBroadBackboneSources;
using h5reader::rediscover::SourceKind;
using h5reader::rediscover::SourceSlot;

namespace {

constexpr double kTol = 1e-9;

// Axes orthonormal + right-handed (y = z × x, x = y × z) to tolerance.
void assertOrthonormalRightHanded(const LocalFrame& f) {
    QVERIFY(std::abs(f.x.norm() - 1.0) < kTol);
    QVERIFY(std::abs(f.y.norm() - 1.0) < kTol);
    QVERIFY(std::abs(f.z.norm() - 1.0) < kTol);
    QVERIFY(std::abs(f.x.dot(f.y)) < kTol);
    QVERIFY(std::abs(f.x.dot(f.z)) < kTol);
    QVERIFY(std::abs(f.y.dot(f.z)) < kTol);
    // Right-handed: x × y ≈ z.
    const Vec3 cross = f.x.cross(f.y);
    QVERIFY((cross - f.z).norm() < kTol);
}

// No coordinate of any axis is NaN/Inf (the fail-loud-not-silent guard).
void assertNoNaN(const LocalFrame& f) {
    for (const Vec3* v : {&f.x, &f.y, &f.z})
        for (int i = 0; i < 3; ++i) QVERIFY(std::isfinite((*v)[i]));
}

// The in-plane projection of `dirLab` (the lab anchor direction) onto the
// frame's x–y plane points along +x (positive x-component, ~zero y-component).
void assertXAzimuthToward(const LocalFrame& f, const Vec3& dirLab) {
    const double xComp = dirLab.dot(f.x);
    const double yComp = dirLab.dot(f.y);
    QVERIFY(xComp > 0.0);                  // x points toward the anchor side
    QVERIFY(std::abs(yComp) < kTol);       // anchor has no y component (in x–z plane)
}

void assertVecNear(const Vec3& actual, const Vec3& expected, double tol = kTol) {
    QVERIFY((actual - expected).norm() < tol);
}

void assertMatNear(const Mat3& actual, const Mat3& expected, double tol = kTol) {
    QVERIFY((actual - expected).norm() < tol);
}

Mat3 fixedRotation() {
    const Vec3 axis = Vec3(1.0, 2.0, 3.0).normalized();
    const double angle = 0.7;
    const Mat3 k = (Mat3() << 0.0, -axis.z(), axis.y(),
                                axis.z(), 0.0, -axis.x(),
                                -axis.y(), axis.x(), 0.0).finished();
    return Mat3::Identity() + std::sin(angle) * k
           + (1.0 - std::cos(angle)) * (k * k);
}

SourceSlot source(SourceKind kind, double dipolar) {
    SourceSlot s;
    s.kind = kind;
    s.dipolar = dipolar;
    return s;
}

SourceSlot charge(double q, const Vec3& disp) {
    SourceSlot s;
    s.kind = SourceKind::Charge;
    s.source_q_e = q;
    s.disp_local = disp;
    s.r = disp.norm();
    return s;
}

}  // namespace

class RediscoverBackboneFrameTests : public QObject {
    Q_OBJECT

private slots:
    void nFrameInterior();
    void nFrameNTerminus();
    void caFrame();
    void carbonylCFrame();
    void carbonylOFrame();
    void haFrame();
    void degenerateCoincident();
    void degenerateCollinearCa();
    void rotationEquivarianceInvariance();
    void broadReducerSumsSources();
};

// ── N frame: z = unit(CA − N); x in-plane toward C_prev ─────────────────────
void RediscoverBackboneFrameTests::nFrameInterior() {
    const Vec3 nPos(0, 0, 0);
    const Vec3 caPos(2, 0, 0);     // CA along +x → z must be +x
    const Vec3 cPrev(0, 3, 0);     // C_prev along +y from N → x azimuth toward +y

    const LocalFrame f = BuildBackboneNFrame(nPos, caPos, cPrev, /*c_prev_valid=*/true);
    QVERIFY(f.is_valid);
    QCOMPARE(static_cast<int>(f.variant), static_cast<int>(FrameVariant::BackboneN));
    assertNoNaN(f);
    assertOrthonormalRightHanded(f);
    // z along N→CA (+x).
    QVERIFY((f.z - Vec3(1, 0, 0)).norm() < kTol);
    // x azimuth toward the N→C_prev direction (+y), in-plane projected.
    assertXAzimuthToward(f, cPrev - nPos);
}

void RediscoverBackboneFrameTests::nFrameNTerminus() {
    const Vec3 nPos(0, 0, 0);
    const Vec3 caPos(0, 0, 2);     // CA along +z → z must be +z
    const Vec3 cOwn(0, 4, 0);      // C_own along +y → x azimuth toward +y

    const LocalFrame f = BuildBackboneNFrame(nPos, caPos, cOwn, /*c_prev_valid=*/false);
    QVERIFY(f.is_valid);
    QCOMPARE(static_cast<int>(f.variant), static_cast<int>(FrameVariant::BackboneN_NTerminus));
    assertNoNaN(f);
    assertOrthonormalRightHanded(f);
    QVERIFY((f.z - Vec3(0, 0, 1)).norm() < kTol);
    assertXAzimuthToward(f, cOwn - nPos);
}

// ── Cα frame: z = bisector(Cα→N, Cα→C); x in-plane toward N ──────────────────
void RediscoverBackboneFrameTests::caFrame() {
    const Vec3 caPos(0, 0, 0);
    const Vec3 nPos(1, 0, 0);      // Cα→N = +x
    const Vec3 cPos(0, 1, 0);      // Cα→C = +y → bisector = (+x+y)/√2

    const LocalFrame f = BuildBackboneCaFrame(caPos, nPos, cPos);
    QVERIFY(f.is_valid);
    QCOMPARE(static_cast<int>(f.variant), static_cast<int>(FrameVariant::BackboneCA));
    assertNoNaN(f);
    assertOrthonormalRightHanded(f);
    // z = unit bisector of (Cα→N=+x) and (Cα→C=+y) = (1,1,0)/√2.
    const Vec3 zExpect = Vec3(1, 1, 0).normalized();
    QVERIFY((f.z - zExpect).norm() < kTol);
    // x azimuth toward Cα→N (+x), in-plane projected (its component on z removed).
    assertXAzimuthToward(f, nPos - caPos);
}

// ── carbonyl C frame: z = unit(O − C); x in-plane toward CA ──────────────────
void RediscoverBackboneFrameTests::carbonylCFrame() {
    const Vec3 cPos(0, 0, 0);
    const Vec3 oPos(0, 0, 3);      // C→O = +z → z must be +z (McConnell reference)
    const Vec3 caPos(5, 0, 0);     // C→CA = +x → x azimuth toward +x

    const LocalFrame f = BuildBackboneCarbonylCFrame(cPos, oPos, caPos);
    QVERIFY(f.is_valid);
    QCOMPARE(static_cast<int>(f.variant), static_cast<int>(FrameVariant::BackboneCarbonylC));
    assertNoNaN(f);
    assertOrthonormalRightHanded(f);
    QVERIFY((f.z - Vec3(0, 0, 1)).norm() < kTol);
    assertXAzimuthToward(f, caPos - cPos);
}

// ── carbonyl O frame: z = unit(C − O); x in-plane toward CA ──────────────────
void RediscoverBackboneFrameTests::carbonylOFrame() {
    const Vec3 oPos(0, 0, 0);
    const Vec3 cPos(0, 0, 2);      // O→C = +z → z must be +z
    const Vec3 caPos(0, 4, 2);     // CA−C = +y → x azimuth toward +y

    const LocalFrame f = BuildBackboneCarbonylOFrame(oPos, cPos, caPos);
    QVERIFY(f.is_valid);
    QCOMPARE(static_cast<int>(f.variant), static_cast<int>(FrameVariant::BackboneCarbonylO));
    assertNoNaN(f);
    assertOrthonormalRightHanded(f);
    QVERIFY((f.z - Vec3(0, 0, 1)).norm() < kTol);
    // x azimuth is toward (CA − C), the in-plane reference, projected.
    assertXAzimuthToward(f, caPos - cPos);
}

// ── HA frame: z = unit(HA − Cα); x in-plane toward N ─────────────────────────
void RediscoverBackboneFrameTests::haFrame() {
    const Vec3 caPos(0, 0, 0);
    const Vec3 haPos(0, 0, 1);     // Cα→HA = +z → z must be +z
    const Vec3 nPos(2, 0, 0);      // Cα→N = +x → x azimuth toward +x

    const LocalFrame f = BuildBackboneHaFrame(haPos, caPos, nPos);
    QVERIFY(f.is_valid);
    QCOMPARE(static_cast<int>(f.variant), static_cast<int>(FrameVariant::BackboneHA));
    assertNoNaN(f);
    assertOrthonormalRightHanded(f);
    QVERIFY((f.z - Vec3(0, 0, 1)).norm() < kTol);
    assertXAzimuthToward(f, nPos - caPos);
}

// ── Degenerate: coincident z-defining atoms → invalid, no NaN ────────────────
void RediscoverBackboneFrameTests::degenerateCoincident() {
    const Vec3 p(1, 2, 3);
    // N frame with N == CA (zero N→CA bond).
    const LocalFrame fN = BuildBackboneNFrame(p, p, Vec3(0, 5, 0), true);
    QVERIFY(!fN.is_valid);
    QCOMPARE(static_cast<int>(fN.variant), static_cast<int>(FrameVariant::Invalid));
    assertNoNaN(fN);
    // HA frame with HA == Cα (zero Cα→HA bond).
    const LocalFrame fHa = BuildBackboneHaFrame(p, p, Vec3(2, 0, 0));
    QVERIFY(!fHa.is_valid);
    assertNoNaN(fHa);
    // Carbonyl C frame with C == O (zero carbonyl bond).
    const LocalFrame fC = BuildBackboneCarbonylCFrame(p, p, Vec3(5, 0, 0));
    QVERIFY(!fC.is_valid);
    assertNoNaN(fC);
    // N frame with the in-plane reference parallel to z (C_prev colinear on N→CA):
    // ref = (C_prev − N) ∥ z ⇒ no in-plane axis ⇒ invalid, no NaN.
    const LocalFrame fNpar = BuildBackboneNFrame(Vec3(0, 0, 0), Vec3(2, 0, 0),
                                                 Vec3(5, 0, 0), true);
    QVERIFY(!fNpar.is_valid);
    assertNoNaN(fNpar);
}

// ── Degenerate: Cα bisector collinear (N, Cα, C anti-parallel) → invalid ─────
void RediscoverBackboneFrameTests::degenerateCollinearCa() {
    const Vec3 caPos(0, 0, 0);
    const Vec3 nPos(1, 0, 0);      // Cα→N = +x
    const Vec3 cPos(-1, 0, 0);     // Cα→C = −x → bisector = 0 (degenerate)
    const LocalFrame f = BuildBackboneCaFrame(caPos, nPos, cPos);
    QVERIFY(!f.is_valid);
    QCOMPARE(static_cast<int>(f.variant), static_cast<int>(FrameVariant::Invalid));
    assertNoNaN(f);
}

void RediscoverBackboneFrameTests::rotationEquivarianceInvariance() {
    const Vec3 caPos(1.0, -2.0, 0.5);
    const Vec3 nPos(2.5, -1.0, 0.2);
    const Vec3 cPos(0.2, 1.2, 1.8);
    const LocalFrame f = BuildBackboneCaFrame(caPos, nPos, cPos);
    QVERIFY(f.is_valid);

    const Mat3 r = fixedRotation();
    const LocalFrame rf = BuildBackboneCaFrame(r * caPos, r * nPos, r * cPos);
    QVERIFY(rf.is_valid);

    assertVecNear(rf.x, r * f.x);
    assertVecNear(rf.y, r * f.y);
    assertVecNear(rf.z, r * f.z);

    const Vec3 dispLab(3.4, -2.2, 5.0);
    assertVecNear(rf.ToLocal(r * dispLab), f.ToLocal(dispLab));

    Mat3 sigma;
    sigma << 11.0, 1.2, -0.4,
             -0.7, 8.0, 2.1,
             0.3, -1.5, 6.5;
    const Mat3 rotatedSigma = r * sigma * r.transpose();
    assertMatNear(rf.TensorToLocal(rotatedSigma), f.TensorToLocal(sigma));

    const SphericalTensor lab = DecomposeLibrary(sigma);
    const SphericalTensor rotatedLab = DecomposeLibrary(rotatedSigma);
    QVERIFY(std::abs(rotatedLab.T0 - lab.T0) < kTol);
    QVERIFY(std::abs(rotatedLab.T2Magnitude() - lab.T2Magnitude()) < kTol);
}

void RediscoverBackboneFrameTests::broadReducerSumsSources() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<SourceSlot> sources;
    sources.push_back(source(SourceKind::Ring, 2.5));
    sources.push_back(source(SourceKind::Ring, -0.25));
    sources.push_back(source(SourceKind::Ring, nan));  // filtered before run; ignored if present
    sources.push_back(source(SourceKind::Bond, 4.0));
    sources.push_back(source(SourceKind::Bond, -1.25));
    sources.back().mc_source_is_self_or_bonded = true;
    sources.push_back(charge(2.0, Vec3(1.0, 2.0, 2.0)));
    sources.push_back(charge(-0.5, Vec3(-2.0, 0.0, 1.0)));

    const BroadAggregate agg =
        ReduceBroadBackboneSources(sources, 8.0, 9.0, 10.0, QStringLiteral("ff14sb"), 0.5);

    QCOMPARE(agg.ring_n, 2);
    QVERIFY(std::abs(agg.ring_sum_dipolar - 2.25) < kTol);
    QCOMPARE(agg.bond_n, 2);
    QVERIFY(std::abs(agg.bond_sum_dipolar - 2.75) < kTol);
    QCOMPARE(agg.bond_n_valid, 1);
    QVERIFY(std::abs(agg.bond_sum_dipolar_valid - 4.0) < kTol);
    QCOMPARE(agg.charge_n, 2);
    QCOMPARE(agg.charge_source, QStringLiteral("ff14sb"));
    QVERIFY(std::abs(agg.ring_cutoff_A - 8.0) < kTol);
    QVERIFY(std::abs(agg.bond_cutoff_A - 9.0) < kTol);
    QVERIFY(std::abs(agg.charge_cutoff_A - 10.0) < kTol);
    QVERIFY(std::abs(agg.mc_near_field_ratio - 0.5) < kTol);

    Vec3 expectedField = Vec3::Zero();
    Vec3 expectedMu = Vec3::Zero();
    for (const SourceSlot& s : sources) {
        if (s.kind != SourceKind::Charge) continue;
        const double r3 = s.r * s.r * s.r;
        expectedField += (-s.source_q_e / r3) * s.disp_local;
        expectedMu += s.source_q_e * s.disp_local;
    }
    assertVecNear(agg.charge_field_local, expectedField);
    assertVecNear(agg.charge_mu_local, expectedMu);
    QVERIFY(std::abs(agg.charge_field_z - expectedField.z()) < kTol);
    QVERIFY(std::abs(agg.charge_field_mag - expectedField.norm()) < kTol);
}

QTEST_GUILESS_MAIN(RediscoverBackboneFrameTests)

#include "rediscover_backbone_frame_tests.moc"
