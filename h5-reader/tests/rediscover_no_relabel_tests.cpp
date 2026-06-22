#include <QtTest/QtTest>

#include <QFile>
#include <QString>
#include <QTextStream>

class RediscoverNoRelabelTests : public QObject {
    Q_OBJECT

private slots:
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
        QVERIFY(source.contains(QStringLiteral("subspaceCompareObject")));
    }
};

QTEST_MAIN(RediscoverNoRelabelTests)
#include "rediscover_no_relabel_tests.moc"
