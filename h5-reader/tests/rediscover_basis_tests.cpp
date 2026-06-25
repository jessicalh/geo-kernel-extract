// rediscover_basis_tests — pin SphericalBasis::DecomposeLibrary against
// analytic values, matching the library's T2 component order
// ([xy, yz, zz, xz, xx−yy], isometric real-SH normalization).
//
// The load-bearing fixture (GUIDANCE.md / DESIGN.md): the dipolar tensor
// 3ẑẑᵀ − I = diag(−1,−1,2) has trace 0 ⇒ T0 = 0; Szz = 2 ⇒
// T2[2] = √(3/2)·2 = √6 ≈ 2.449 (NOT √(3/2), which would be (3ẑẑᵀ−I)/2);
// all other T2 components 0; T1 = 0.

#include "rediscover/SphericalBasis.h"
#include "rediscover/RunData.h"

#include <QtTest>

#include <cmath>
#include <memory>
#include <stdexcept>

using h5reader::model::Mat3;
using h5reader::model::SphericalTensor;
using h5reader::model::Vec3;
using h5reader::rediscover::DftFrameSet;
using h5reader::rediscover::DecomposeLibrary;
using h5reader::rediscover::ReconstructLibraryT2Matrix;

class RediscoverBasisTests : public QObject {
    Q_OBJECT

private slots:
    void dipolarZZ();
    void offDiagonalComponents();
    void isometricNormPreserved();
    void traceIsotropic();
    void reconstructT2RoundTrip();
    void molecularFrameRoundTripAuditMath();
    void caseHunterSelectionGuardRejectsDftRead();
};

void RediscoverBasisTests::dipolarZZ() {
    Mat3 m = Mat3::Zero();
    m(0, 0) = -1.0;
    m(1, 1) = -1.0;
    m(2, 2) = 2.0;  // 3 zz^T - I

    const SphericalTensor st = DecomposeLibrary(m);

    QVERIFY(std::abs(st.T0) < 1e-12);
    for (double v : st.T1) QVERIFY(std::abs(v) < 1e-12);

    // T2[2] (m=0, zz) = sqrt(3/2) * Szz = sqrt(3/2)*2 = sqrt(6).
    QVERIFY(std::abs(st.T2[2] - std::sqrt(6.0)) < 1e-9);
    // Every other T2 component is zero (purely axial about z).
    QVERIFY(std::abs(st.T2[0]) < 1e-12);
    QVERIFY(std::abs(st.T2[1]) < 1e-12);
    QVERIFY(std::abs(st.T2[3]) < 1e-12);
    QVERIFY(std::abs(st.T2[4]) < 1e-12);
}

void RediscoverBasisTests::offDiagonalComponents() {
    // Symmetric off-diagonal probe: each component lands in the named slot
    // with the sqrt(2) isometric factor, library order [xy, yz, zz, xz, xx-yy].
    Mat3 m = Mat3::Zero();
    m(0, 1) = m(1, 0) = 1.0;  // Sxy = 1
    m(1, 2) = m(2, 1) = 2.0;  // Syz = 2
    m(0, 2) = m(2, 0) = 3.0;  // Sxz = 3

    const SphericalTensor st = DecomposeLibrary(m);
    const double s2 = std::sqrt(2.0);
    QVERIFY(std::abs(st.T2[0] - s2 * 1.0) < 1e-9);  // xy
    QVERIFY(std::abs(st.T2[1] - s2 * 2.0) < 1e-9);  // yz
    QVERIFY(std::abs(st.T2[3] - s2 * 3.0) < 1e-9);  // xz
    QVERIFY(std::abs(st.T0) < 1e-12);
}

void RediscoverBasisTests::isometricNormPreserved() {
    // Frobenius norm of the traceless-symmetric part equals |T2|.
    Mat3 m;
    m << 1.0, 0.4, -0.7,
         0.4, -2.0, 0.9,
        -0.7, 0.9, 1.0;  // symmetric; trace 0 -> already traceless

    const SphericalTensor st = DecomposeLibrary(m);

    // Build the traceless-symmetric part explicitly.
    const double tr = m.trace() / 3.0;
    Mat3 s = 0.5 * (m + m.transpose());
    s(0, 0) -= tr; s(1, 1) -= tr; s(2, 2) -= tr;
    double frob2 = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) frob2 += s(i, j) * s(i, j);

    double t2sq = 0.0;
    for (double v : st.T2) t2sq += v * v;
    QVERIFY(std::abs(t2sq - frob2) < 1e-9);
}

void RediscoverBasisTests::traceIsotropic() {
    Mat3 m = Mat3::Identity() * 30.0;  // pure isotropic
    const SphericalTensor st = DecomposeLibrary(m);
    QVERIFY(std::abs(st.T0 - 30.0) < 1e-12);
    QVERIFY(st.T2Magnitude() < 1e-12);
}

void RediscoverBasisTests::reconstructT2RoundTrip() {
    const std::array<double, 5> t2 = {0.25, -1.5, std::sqrt(6.0), 0.75, -0.4};
    const Mat3 m = ReconstructLibraryT2Matrix(12.0, t2);
    const SphericalTensor st = DecomposeLibrary(m);

    QVERIFY(std::abs(st.T0 - 12.0) < 1e-9);
    for (double v : st.T1) QVERIFY(std::abs(v) < 1e-12);
    for (std::size_t i = 0; i < t2.size(); ++i)
        QVERIFY(std::abs(st.T2[i] - t2[i]) < 1e-9);

    const Mat3 dipolar = ReconstructLibraryT2Matrix({0.0, 0.0, std::sqrt(6.0), 0.0, 0.0});
    QVERIFY(std::abs(dipolar(0, 0) + 1.0) < 1e-9);
    QVERIFY(std::abs(dipolar(1, 1) + 1.0) < 1e-9);
    QVERIFY(std::abs(dipolar(2, 2) - 2.0) < 1e-9);
}

void RediscoverBasisTests::molecularFrameRoundTripAuditMath() {
    Vec3 x(1.0, 2.0, 0.5);
    x.normalize();
    Vec3 zSeed(-0.25, 0.4, 1.0);
    Vec3 z = zSeed - zSeed.dot(x) * x;
    z.normalize();
    Vec3 y = z.cross(x);
    y.normalize();
    z = x.cross(y);
    z.normalize();

    Mat3 axes;
    axes.col(0) = x;
    axes.col(1) = y;
    axes.col(2) = z;
    QVERIFY((axes.transpose() * axes - Mat3::Identity()).cwiseAbs().maxCoeff() < 1e-13);
    QVERIFY(std::abs(axes.determinant() - 1.0) < 1e-13);

    Mat3 sigma;
    sigma << 101.0, 1.2, -0.7,
             0.4, 94.0, 2.1,
            -0.2, 1.4, 88.0;
    const Mat3 sigmaMol = axes.transpose() * sigma * axes;
    const Mat3 roundTrip = axes * sigmaMol * axes.transpose();
    QVERIFY((roundTrip - sigma).cwiseAbs().maxCoeff() < 1e-13);
}

void RediscoverBasisTests::caseHunterSelectionGuardRejectsDftRead() {
    DftFrameSet dft;
    auto frame = std::make_shared<h5reader::model::DftShieldingFrame>();
    frame->atoms.resize(1);
    dft.Insert(7, frame);
    (void)dft.frameCount();

    bool threw = false;
    try {
        const DftFrameSet::SelectionReadGuard guard(dft);
        (void)dft.AtomShielding(0, 7);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    QVERIFY(threw);

    QVERIFY(dft.AtomShielding(0, 7) != nullptr);
}

QTEST_GUILESS_MAIN(RediscoverBasisTests)

#include "rediscover_basis_tests.moc"
