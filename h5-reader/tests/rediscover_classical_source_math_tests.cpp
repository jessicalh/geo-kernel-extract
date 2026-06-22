#include "rediscover/ClassicalSourceMath.h"
#include "rediscover/TensorConventionGuard.h"

#include <QtTest/QtTest>

#include <QByteArray>

#include <cmath>
#include <limits>
#include <vector>

using namespace h5reader::rediscover;

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

    void larsenTermParticipatesInClassicalFold() {
        const double withoutLarsen =
            FoldClassicalSigma(100.0, {{1.0, true}, {2.0, true}, {7.0, false}});
        const double withLarsen =
            FoldClassicalSigma(100.0, {{1.0, true}, {2.0, true}, {7.0, true}});
        QVERIFY(std::abs(withLarsen - 110.0) < 1e-12);
        QVERIFY(std::abs((withLarsen - withoutLarsen) - 7.0) < 1e-12);
    }

    void pasConventionRefusesSignedSortedOverride() {
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
        QVERIFY(AssertPasShapeConventionEnv(&err));
        QVERIFY(err.isEmpty());
        qunsetenv("H5READER_PAS_SHAPE_CONVENTION_OVERRIDE");
    }
};

QTEST_MAIN(RediscoverClassicalSourceMathTests)
#include "rediscover_classical_source_math_tests.moc"
