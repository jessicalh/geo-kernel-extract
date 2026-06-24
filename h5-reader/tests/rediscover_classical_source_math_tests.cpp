#include "rediscover/ClassicalSourceMath.h"
#include "rediscover/LiteratureConstants.h"
#include "rediscover/McConnellLiteratureKernel.h"
#include "rediscover/RingCurrentScalars.h"
#include "rediscover/TensorConventionGuard.h"

#include "calculators/QtBiotSavartCalc.h"

#include <QtTest/QtTest>

#include <QByteArray>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace h5reader::rediscover;
namespace calculators = h5reader::calculators;
namespace model = h5reader::model;

namespace {

double finiteOlsSlope(const std::vector<double>& x, const std::vector<double>& y) {
    double sx = 0.0;
    double sy = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < x.size() && i < y.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
        sx += x[i];
        sy += y[i];
        ++n;
    }
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);
    double sxx = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size() && i < y.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
        const double dx = x[i] - mx;
        sxx += dx * dx;
        sxy += dx * (y[i] - my);
    }
    return sxx > 0.0 ? sxy / sxx : std::numeric_limits<double>::quiet_NaN();
}

SourceSlot peptideCoSource(double r, double cosTheta) {
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const model::Vec3 dHat(sinTheta, 0.0, cosTheta);
    SourceSlot s;
    s.kind = SourceKind::Bond;
    s.disp_local = -r * dHat;
    s.r = r;
    s.cos_theta = cosTheta;
    s.dipolar = (3.0 * cosTheta * cosTheta - 1.0) / (r * r * r);
    s.bond_category = static_cast<int>(model::BondCategory::PeptideCO);
    s.bond_axis_local = model::Vec3::UnitZ();
    return s;
}

}  // namespace

class RediscoverClassicalSourceMathTests : public QObject {
    Q_OBJECT

private slots:
    void scaleFactorIsSdRatioNotSlope() {
        const std::vector<double> qm = {-1.0, -1.0, 1.0, 1.0};
        const double q = std::sqrt(0.75);
        const std::vector<double> cl = {
            0.1 * (-0.5 - q),
            0.1 * (-0.5 + q),
            0.1 * (0.5 - q),
            0.1 * (0.5 + q),
        };

        const double scale = SdRatioScaleFactor(qm, cl);
        const double slope = finiteOlsSlope(cl, qm);
        QVERIFY(std::abs(scale - 10.0) < 1e-12);
        QVERIFY(std::abs(slope - 5.0) < 1e-12);
        QVERIFY(std::abs(scale - slope) > 1.0);
    }

    void ringContributionUsesSignedIsotropicKernel() {
        model::SphericalTensor fixed;
        fixed.T0 = -2.5;
        fixed.T2 = {3.0, 4.0, 0.0, 0.0, 0.0};

        const double t2Magnitude = std::sqrt(3.0 * 3.0 + 4.0 * 4.0);
        QCOMPARE(RingForwardContributionPpm(fixed), -2.5);
        QVERIFY(RingForwardContributionPpm(fixed) != t2Magnitude);
    }

    void mcconnellContributionUsesSignedIsotropicKernel() {
        model::SphericalTensor fixed;
        fixed.T0 = -1.25;
        fixed.T2 = {3.0, 4.0, 0.0, 0.0, 0.0};

        const double t2Magnitude = std::sqrt(3.0 * 3.0 + 4.0 * 4.0);
        QCOMPARE(McConnellForwardContributionPpm(fixed), -1.25);
        QVERIFY(McConnellForwardContributionPpm(fixed) != t2Magnitude);
    }

    void sharedClassicalSigmaUsesSignedMopacFieldAndAllTerms() {
        const ClassicalSigmaResult out = ComputeClassicalSigma({
            {100.0, true},
            {-2.0, true},
            {3.0, true},
            {0.5, true},
            {4.0, true},
            {-1.5, true},
            {2.25, true},
        });
        QVERIFY(out.buckingham_linear.present);
        QVERIFY(out.buckingham_quadratic.present);
        QVERIFY(out.buckingham.present);
        QCOMPARE(out.buckingham_linear.value, 6.0);
        QCOMPARE(out.buckingham_quadratic.value, -2.0);
        QCOMPARE(out.buckingham.value, 4.0);
        QCOMPARE(out.sigma_cl.value, 108.75);
    }

    void sharedClassicalSigmaOmitsAbsentBuckinghamAxis() {
        const ClassicalSigmaResult out = ComputeClassicalSigma({
            {100.0, true},
            {std::numeric_limits<double>::quiet_NaN(), false},
            {3.0, true},
            {0.5, true},
            {4.0, true},
            {-1.5, true},
            {2.25, true},
        });
        QVERIFY(!out.buckingham.present);
        QCOMPARE(out.sigma_cl.value, 104.75);
    }

    void mcconnellProducerT0ScalingUsesCategoryDeltaChi() {
        const double co = McConnellProducerT0ToPpm(model::BondCategory::PeptideCO, 0.125);
        const double cn = McConnellProducerT0ToPpm(model::BondCategory::PeptideCN, 0.125);
        const double aromatic = McConnellProducerT0ToPpm(model::BondCategory::Aromatic, 0.125);
        QVERIFY(std::abs(co - (-kMcConnellMolarPrefactor.value
                               * McConnellDeltaChi(model::BondCategory::PeptideCO).value
                               * 0.125)) < 1e-15);
        QVERIFY(co < 0.0);
        QVERIFY(cn > 0.0);
        QCOMPARE(aromatic, 0.0);
    }

    void mcconnellMagicAngleSignFollowsSignedT0() {
        constexpr double r = 3.0;
        const double prefactor = McConnellMolarPrefactor();
        const double deltaChi = McConnellDeltaChi(model::BondCategory::PeptideCO).value;

        bool facePresent = false;
        const model::SphericalTensor face =
            McConnellSourceLiteratureKernelLocal(peptideCoSource(r, 1.0), &facePresent);
        QVERIFY(facePresent);
        const double faceExpected = -prefactor * deltaChi * (3.0 * 1.0 * 1.0 - 1.0)
                                    / (3.0 * r * r * r);
        QVERIFY(std::abs(face.T0 - faceExpected) < 1e-14);
        QVERIFY(face.T0 < 0.0);

        bool inPlanePresent = false;
        const model::SphericalTensor inPlane =
            McConnellSourceLiteratureKernelLocal(peptideCoSource(r, 0.0), &inPlanePresent);
        QVERIFY(inPlanePresent);
        const double inPlaneExpected = -prefactor * deltaChi * (3.0 * 0.0 * 0.0 - 1.0)
                                       / (3.0 * r * r * r);
        QVERIFY(std::abs(inPlane.T0 - inPlaneExpected) < 1e-14);
        QVERIFY(inPlane.T0 > 0.0);

        const double magicCos = 1.0 / std::sqrt(3.0);
        bool magicPresent = false;
        const model::SphericalTensor magic =
            McConnellSourceLiteratureKernelLocal(peptideCoSource(r, magicCos), &magicPresent);
        QVERIFY(magicPresent);
        QVERIFY(std::abs(magic.T0) < 1e-14);
    }

    void ringCurrentSignConventionShieldsFaceDeshieldsInPlane() {
        model::RingGeometry geo;
        geo.center = model::Vec3::Zero();
        geo.normal = model::Vec3::UnitZ();
        geo.radius = 1.4;

        std::vector<model::Vec3> vertices;
        vertices.reserve(6);
        constexpr double kPi = 3.14159265358979323846264338327950288;
        for (int i = 0; i < 6; ++i) {
            const double theta = (2.0 * kPi * static_cast<double>(i)) / 6.0;
            vertices.emplace_back(geo.radius * std::cos(theta),
                                  geo.radius * std::sin(theta),
                                  0.0);
        }

        const double intensity = RingIntensity(model::RingTypeIndex::PheBenzene).value;
        const double lobeOffset = JohnsonBoveyLobeOffset(model::RingTypeIndex::PheBenzene).value;
        const model::SphericalTensor face =
            calculators::EvaluateShielding(model::Vec3(0.0, 0.0, 3.0), geo, vertices,
                                           lobeOffset, intensity);
        const model::SphericalTensor inPlane =
            calculators::EvaluateShielding(model::Vec3(3.0, 0.0, 0.0), geo, vertices,
                                           lobeOffset, intensity);

        QVERIFY2(face.T0 > 0.0,
                 qPrintable(QStringLiteral("face T0 should be shielding-positive, got %1")
                                .arg(face.T0, 0, 'g', 17)));
        QVERIFY2(inPlane.T0 < 0.0,
                 qPrintable(QStringLiteral("in-plane T0 should be deshielding-negative, got %1")
                                .arg(inPlane.T0, 0, 'g', 17)));
    }

    void ringPerTypeWeightingUsesLiteratureIntensities() {
        std::array<double, model::kAromaticRingTypeCount> perType{};
        perType[static_cast<std::size_t>(model::RingTypeIndex::PheBenzene)] = 0.25;
        const double phe = RingPerTypeT0Ppm(perType.data(), perType.size());

        perType.fill(0.0);
        perType[static_cast<std::size_t>(model::RingTypeIndex::HisImidazole)] = 0.25;
        const double his = RingPerTypeT0Ppm(perType.data(), perType.size());

        perType.fill(0.0);
        perType[static_cast<std::size_t>(model::RingTypeIndex::TrpPyrrole)] = 0.25;
        const double trpPyrrole = RingPerTypeT0Ppm(perType.data(), perType.size());

        const double hisRatio = std::abs(his / phe);
        const double trpPyrroleRatio = std::abs(trpPyrrole / phe);
        QVERIFY(std::abs(hisRatio - (5.16 / 12.0)) < 1e-12);
        QVERIFY(std::abs(trpPyrroleRatio - (6.72 / 12.0)) < 1e-12);
    }

    void larsenTermParticipatesInClassicalFold() {
        const double withoutLarsen =
            FoldClassicalSigma(100.0, {{1.0, true}, {2.0, true}, {7.0, false}});
        const double withLarsen =
            FoldClassicalSigma(100.0, {{1.0, true}, {2.0, true}, {7.0, true}});
        QVERIFY(std::abs(withLarsen - 110.0) < 1e-12);
        QVERIFY(std::abs((withLarsen - withoutLarsen) - 7.0) < 1e-12);
    }

    void pasConventionRequiresPrincipalDescendingOverride() {
        qunsetenv("H5READER_PAS_SHAPE_CONVENTION_OVERRIDE");
        QString err;
        QVERIFY(AssertPasShapeConventionEnv(&err));
        QVERIFY(err.isEmpty());

        qputenv("H5READER_PAS_SHAPE_CONVENTION_OVERRIDE",
                QByteArrayLiteral("signed_sorted"));
        QVERIFY(!AssertPasShapeConventionEnv(&err));
        QVERIFY(err.contains(QStringLiteral("unsupported PAS shape convention")));

        qputenv("H5READER_PAS_SHAPE_CONVENTION_OVERRIDE",
                QByteArrayLiteral("haeberlen_distance_from_isotropic_v1"));
        err.clear();
        QVERIFY(!AssertPasShapeConventionEnv(&err));
        QVERIFY(err.contains(QStringLiteral("unsupported PAS shape convention")));

        qputenv("H5READER_PAS_SHAPE_CONVENTION_OVERRIDE",
                QByteArrayLiteral("principal_shielding_descending_v1"));
        err.clear();
        QVERIFY(AssertPasShapeConventionEnv(&err));
        QVERIFY(err.isEmpty());
        qunsetenv("H5READER_PAS_SHAPE_CONVENTION_OVERRIDE");
    }
};

QTEST_MAIN(RediscoverClassicalSourceMathTests)
#include "rediscover_classical_source_math_tests.moc"
