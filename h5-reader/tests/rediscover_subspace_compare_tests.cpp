#include "rediscover/SubspaceCompare.h"

#include <QtTest/QtTest>

#include <cmath>
#include <numeric>
#include <vector>

using namespace h5reader::rediscover;

class RediscoverSubspaceCompareTests : public QObject {
    Q_OBJECT

private slots:
    void identicalTwoDimensionalSubspacesHaveUnitCanonicalCorrelations() {
        std::vector<double> x(96);
        std::vector<double> y(96);
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double t = static_cast<double>(i) / 9.0;
            x[i] = std::sin(t);
            y[i] = std::cos(0.7 * t);
        }
        std::vector<double> bx(x.size());
        std::vector<double> by(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            bx[i] = 2.0 * x[i] + 0.25 * y[i];
            by[i] = -0.5 * x[i] + y[i];
        }
        std::vector<std::size_t> rows(x.size());
        std::iota(rows.begin(), rows.end(), 0);

        const SubspaceCompareResult r =
            CompareSubspaces({QStringLiteral("a"),
                              {{QStringLiteral("x"), x}, {QStringLiteral("y"), y}}},
                             {QStringLiteral("b"),
                              {{QStringLiteral("bx"), bx}, {QStringLiteral("by"), by}}},
                             rows);
        QVERIFY(r.computed);
        QCOMPARE(r.provenance, QStringLiteral("svd_subspace_compare_v1"));
        QCOMPARE(r.active_dim_a, 2);
        QCOMPARE(r.active_dim_b, 2);
        QVERIFY(r.basis_dim_a >= 1);
        QVERIFY(r.basis_dim_b >= 1);
        QVERIFY(r.max_canonical_corr > 0.999);
        QVERIFY(r.mean_canonical_corr > 0.999);
    }

    void constantColumnsAreDroppedBeforeBasisSelection() {
        std::vector<double> x(32);
        std::vector<double> c(32, 7.0);
        for (std::size_t i = 0; i < x.size(); ++i)
            x[i] = static_cast<double>(i);
        std::vector<std::size_t> rows(x.size());
        std::iota(rows.begin(), rows.end(), 0);

        const SubspaceCompareResult r =
            CompareSubspaces({QStringLiteral("a"),
                              {{QStringLiteral("variable"), x}, {QStringLiteral("constant"), c}}},
                             {QStringLiteral("b"),
                              {{QStringLiteral("same_variable"), x}}},
                             rows);
        QVERIFY(r.computed);
        QCOMPARE(r.active_dim_a, 1);
        QCOMPARE(r.active_dim_b, 1);
        QVERIFY(r.dropped_channels_a.contains(QStringLiteral("constant")));
        QCOMPARE(r.basis_dim_a, 1);
        QCOMPARE(r.basis_dim_b, 1);
        QVERIFY(r.max_canonical_corr > 0.999);
    }

    void scalarFamiliesStillUseSvdProvenance() {
        std::vector<double> x(48);
        std::vector<double> y(48);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = static_cast<double>(i % 11) - 5.0;
            y[i] = -3.0 * x[i] + 2.0;
        }
        std::vector<std::size_t> rows(x.size());
        std::iota(rows.begin(), rows.end(), 0);

        const SubspaceCompareResult r =
            CompareSubspaces({QStringLiteral("a"), {{QStringLiteral("x"), x}}},
                             {QStringLiteral("b"), {{QStringLiteral("y"), y}}},
                             rows);
        QVERIFY(r.computed);
        QCOMPARE(r.provenance, QStringLiteral("svd_subspace_compare_v1"));
        QCOMPARE(r.basis_dim_a, 1);
        QCOMPARE(r.basis_dim_b, 1);
        QVERIFY(r.max_canonical_corr > 0.999);
        QVERIFY(r.min_angle_deg < 0.01);
    }
};

QTEST_MAIN(RediscoverSubspaceCompareTests)
#include "rediscover_subspace_compare_tests.moc"
