#include "rediscover/SubspaceCompare.h"

#include <QtTest/QtTest>

#include <QFile>
#include <QJsonObject>
#include <QString>
#include <QTextStream>

#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace h5reader::rediscover;

namespace {

double pearson(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    double sx = 0.0;
    double sy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sx += x[i];
        sy += y[i];
    }
    const double mx = sx / static_cast<double>(x.size());
    const double my = sy / static_cast<double>(y.size());
    double sxx = 0.0;
    double syy = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    const double denom = std::sqrt(sxx * syy);
    return denom > 0.0 ? sxy / denom : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace

class RediscoverNoRelabelTests : public QObject {
    Q_OBJECT

private slots:
    void canonicalCorrelationJsonComesFromSubspaceComputation() {
        std::vector<double> x(96);
        std::vector<double> y(96);
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double t = static_cast<double>(i) / 7.0;
            x[i] = std::sin(t);
            y[i] = std::cos(t);
        }
        std::vector<double> bx(x.size());
        std::vector<double> by(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            bx[i] = x[i] + y[i];
            by[i] = x[i] - y[i];
        }
        std::vector<std::size_t> rows(x.size());
        std::iota(rows.begin(), rows.end(), 0);

        const double scalarR = pearson(x, bx);
        QVERIFY(std::isfinite(scalarR));
        QVERIFY(std::abs(scalarR) < 0.90);

        const SubspaceCompareResult r =
            CompareSubspaces({QStringLiteral("a"),
                              {{QStringLiteral("x"), x}, {QStringLiteral("y"), y}}},
                             {QStringLiteral("b"),
                              {{QStringLiteral("bx"), bx}, {QStringLiteral("by"), by}}},
                             rows);
        const QJsonObject json = SubspaceCompareJson(r);
        const double maxCanonical = json.value(QStringLiteral("max_canonical_corr")).toDouble();
        QVERIFY(maxCanonical > 0.999);
        QVERIFY(std::abs(maxCanonical - std::abs(scalarR)) > 0.10);
        QCOMPARE(json.value(QStringLiteral("provenance")).toString(),
                 QStringLiteral("svd_subspace_compare_v1"));
    }

    void analysisAtomDoesNotMapPearsonToCanonicalFields() {
        QFile f(QStringLiteral(H5READER_SOURCE_DIR "/src/rediscover/AnalysisAtom.cpp"));
        QVERIFY2(f.open(QIODevice::ReadOnly | QIODevice::Text), qPrintable(f.errorString()));
        const QString source = QString::fromUtf8(f.readAll());

        QVERIFY2(!source.contains(QStringLiteral("basis_dim_a\"), 1")),
                 "basis_dim_a must not be hard-coded from scalar Pearson paths");
        QVERIFY2(!source.contains(QStringLiteral("basis_dim_b\"), 1")),
                 "basis_dim_b must not be hard-coded from scalar Pearson paths");
        QVERIFY2(!source.contains(QStringLiteral("explained_fraction_a\"), 1.0")),
                 "explained_fraction_a must come from SVD energy only");
        QVERIFY2(!source.contains(QStringLiteral("explained_fraction_b\"), 1.0")),
                 "explained_fraction_b must come from SVD energy only");
        QVERIFY2(!source.contains(QStringLiteral("max_canonical_corr\"), std::abs(r)")),
                 "scalar Pearson r must not be relabeled max_canonical_corr");
        QVERIFY2(!source.contains(QStringLiteral("mean_canonical_corr\"), std::abs(r)")),
                 "scalar Pearson r must not be relabeled mean_canonical_corr");
        QVERIFY2(!source.contains(QStringLiteral("contraction_canonical_correlations")),
                 "BS/HM contraction Pearson audits must not be named canonical correlations");

        QVERIFY(source.contains(QStringLiteral("field_efg_pearson_r_unthresholded")));
        QVERIFY(source.contains(QStringLiteral("r_bs_hm_contraction")));
        QVERIFY(source.contains(QStringLiteral("bs_hm_divergence_to_ring_contribution_ratio")));
        QVERIFY(source.contains(QStringLiteral("independent_forms_checked_distinct_kernel_forms_shared_geometry_projection_decomposition")));
        QVERIFY(source.contains(QStringLiteral("subspaceCompareObject")));
    }

    void exposurePassSurfacesUseBoundedGrains() {
        QFile atom(QStringLiteral(H5READER_SOURCE_DIR "/src/rediscover/AnalysisAtom.cpp"));
        QVERIFY2(atom.open(QIODevice::ReadOnly | QIODevice::Text), qPrintable(atom.errorString()));
        const QString source = QString::fromUtf8(atom.readAll());

        QVERIFY(source.contains(QStringLiteral("per_type_tensor_reveals.csv")));
        QVERIFY(source.contains(QStringLiteral("ring_well_target_eta2.csv")));
        QVERIFY(source.contains(QStringLiteral("serial_recurrence_summary.csv")));
        QVERIFY(source.contains(QStringLiteral("per_atom_per_type_record")));
        QVERIFY(source.contains(QStringLiteral("per_atom_ring_well_target")));
        QVERIFY(source.contains(QStringLiteral("no_per_partner_rows;no_frame_index")));
        QVERIFY(source.contains(QStringLiteral("pyramidalization_oop_A")));
        QVERIFY(source.contains(QStringLiteral("E_parallel_XH")));
        QVERIFY(source.contains(QStringLiteral("pas_delta11")));
    }

    void cohortStaticRelationshipRelabelDetectorIsLegacyAware() {
        QFile cohort(QStringLiteral(H5READER_SOURCE_DIR "/src/rediscover/CohortContextAccumulator.cpp"));
        QVERIFY2(cohort.open(QIODevice::ReadOnly | QIODevice::Text), qPrintable(cohort.errorString()));
        const QString source = QString::fromUtf8(cohort.readAll());

        QVERIFY(source.contains(QStringLiteral("cohort_static_source_relationships.csv")));
        QVERIFY(source.contains(QStringLiteral("legacyEquivalentColumnForChannel")));
        QVERIFY(source.contains(QStringLiteral("channel_vs_sigma[%1]")));
        QVERIFY(source.contains(QStringLiteral("full_tensor_r2")));
        QVERIFY(source.contains(QStringLiteral("magnitude_only_r2")));
        QVERIFY(source.contains(QStringLiteral("direction_only_r2")));
        QVERIFY(!source.contains(QStringLiteral("blocked pending predecessor chi1")));
        QVERIFY(source.contains(QStringLiteral("full_run_reference=1148/1494=0.768")));
    }
};

QTEST_MAIN(RediscoverNoRelabelTests)
#include "rediscover_no_relabel_tests.moc"
