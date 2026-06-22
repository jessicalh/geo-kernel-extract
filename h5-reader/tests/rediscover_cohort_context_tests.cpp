#include "rediscover/CohortContextAccumulator.h"
#include "rediscover/DeltaRunData.h"
#include "rediscover/ScopedProducerCatalog.h"

#include <QtTest/QtTest>

#include <cmath>
#include <utility>
#include <vector>

using namespace h5reader::rediscover;

namespace {

StaticNpyArray tensorArray(const QString& stem, std::size_t rows) {
    StaticNpyArray a;
    a.stem = stem;
    a.path = stem + QStringLiteral(".npy");
    a.rows = rows;
    a.cols = 9;
    a.values.assign(rows * 9, 0.0);
    for (std::size_t r = 0; r < rows; ++r) {
        a.values[r * 9 + 0] = 1.0 + static_cast<double>(r);
        a.values[r * 9 + 4] = 2.0 + static_cast<double>(r);
        a.values[r * 9 + 8] = 3.0 + static_cast<double>(r);
    }
    return a;
}

void put(RunData& run, h5reader::io::FieldKind kind, StaticNpyArray a) {
    run.producerArrays[static_cast<int>(kind)] = std::move(a);
}

}  // namespace

class RediscoverCohortContextTests : public QObject {
    Q_OBJECT

private slots:
    void keyBuilderUsesFullIupacIdentity();
    void supportCredentialDoesNotAuthorizeInsufficientCouplings();
    void accumulatorCountsProteinResidency();
    void permutationNullUsesProteinLabelShuffle();
    void helixDipoleHasExpectedSign();
    void deltaLoaderIsNarrowAndRefusesMissingArrays();
};

void RediscoverCohortContextTests::keyBuilderUsesFullIupacIdentity() {
    Axis2ContextKeyFields od1;
    od1.element = QStringLiteral("O");
    od1.residue_type = QStringLiteral("ASP");
    od1.atom_name = QStringLiteral("OD1");
    od1.hyb = QStringLiteral("sp2");
    od1.contact_class = QStringLiteral("neg_charged");
    od1.dihedral_region = QStringLiteral("beta");
    od1.SS = QStringLiteral("sheet");

    Axis2ContextKeyFields od2 = od1;
    od2.atom_name = QStringLiteral("OD2");

    const Axis2ContextKey k1 = BuildAxis2ContextKey(od1);
    const Axis2ContextKey k2 = BuildAxis2ContextKey(od2);
    QVERIFY(k1.canonical != k2.canonical);
    QVERIFY(k1.canonical.contains(QStringLiteral("atom_name=OD1")));
    QVERIFY(k1.canonical.contains(QStringLiteral("residue_type=ASP")));
    QVERIFY(k1.canonical.contains(QStringLiteral("dihedral_region=beta")));
    QVERIFY(k1.canonical.contains(QStringLiteral("SS=sheet")));
    QVERIFY(!k1.canonical.contains(QStringLiteral("graph_stratum")));
    QVERIFY(!k1.canonical.contains(QStringLiteral("n_frames")));

    od1.SS.clear();
    const Axis2ContextKey missing = BuildAxis2ContextKey(od1);
    QVERIFY(missing.canonical.contains(QStringLiteral("SS=unknown")));
}

void RediscoverCohortContextTests::supportCredentialDoesNotAuthorizeInsufficientCouplings() {
    const SupportThresholds t{3, 12};
    const SupportCredential insufficient = CredentialSupport(2, true, t);
    QCOMPARE(insufficient.support_name, QStringLiteral("insufficient"));
    QVERIFY(!insufficient.may_emit_coupling);
    QVERIFY(!insufficient.may_emit_full_subspace);

    const SupportCredential reduced = CredentialSupport(6, true, t);
    QCOMPARE(reduced.support_name, QStringLiteral("reduced"));
    QVERIFY(reduced.may_emit_coupling);
    QVERIFY(!reduced.may_emit_full_subspace);

    const SupportCredential full = CredentialSupport(12, true, t);
    QCOMPARE(full.support_name, QStringLiteral("full"));
    QVERIFY(full.may_emit_coupling);
    QVERIFY(full.may_emit_full_subspace);
}

void RediscoverCohortContextTests::accumulatorCountsProteinResidency() {
    Axis2ContextKeyFields f;
    f.element = QStringLiteral("N");
    f.residue_type = QStringLiteral("GLY");
    f.atom_name = QStringLiteral("N");
    f.hyb = QStringLiteral("sp2");
    f.contact_class = QStringLiteral("polar");
    f.dihedral_region = QStringLiteral("alphaR");
    f.SS = QStringLiteral("helix");

    CohortContextAccumulator acc;
    for (int i = 0; i < 3; ++i) {
        CohortSample s;
        s.key = BuildAxis2ContextKey(f);
        s.protein_id = QStringLiteral("p%1").arg(i);
        s.sigma_iso = 100.0 + i;
        acc.push(s);
    }
    QCOMPARE(acc.cellCount(), std::size_t(1));
    QCOMPARE(acc.sampleCount(), std::size_t(3));
    QCOMPARE(acc.cells().begin()->second.proteins.size(), qsizetype(3));
}

void RediscoverCohortContextTests::permutationNullUsesProteinLabelShuffle() {
    const std::vector<double> x = {1, 2, 3, 4, 5, 6};
    const std::vector<double> y = {2, 4, 6, 8, 10, 12};
    QCOMPARE(LinearSlope(x, y), 2.0);
    const PermutationNull p = ProteinLabelPermutationNull(x, y, 31, 123u);
    QCOMPARE(p.permutation_K, 31);
    QVERIFY(p.perm_p > 0.0);
    QVERIFY(p.perm_p <= 1.0);
    QVERIFY(std::isfinite(p.obs_slope_z));
}

void RediscoverCohortContextTests::helixDipoleHasExpectedSign() {
    HelixDipoleInput input;
    input.ca_z_A = {-3.0, -1.5, 0.0, 1.5};
    input.target_z_A = 5.0;
    QVERIFY(ComputeHelixDipoleField(input) > 0.0);
    input.target_z_A = -5.0;
    QVERIFY(ComputeHelixDipoleField(input) < 0.0);
}

void RediscoverCohortContextTests::deltaLoaderIsNarrowAndRefusesMissingArrays() {
    QVERIFY(!IsScopedProducerField(h5reader::io::FieldSpecFor(h5reader::io::FieldKind::DeltaShielding)));
    QVERIFY(DeltaFieldNames().contains(QStringLiteral("delta_shielding")));

    RunData missing;
    QString err;
    QVERIFY(!LoadDeltaRunData(missing, &err));
    QVERIFY(err.contains(QStringLiteral("wt_shielding_diamagnetic")));

    RunData valid;
    put(valid, h5reader::io::FieldKind::WtShieldingDiamagnetic,
        tensorArray(QStringLiteral("wt_shielding_diamagnetic"), 3));
    put(valid, h5reader::io::FieldKind::WtShieldingParamagnetic,
        tensorArray(QStringLiteral("wt_shielding_paramagnetic"), 3));
    put(valid, h5reader::io::FieldKind::MutShieldingDiamagnetic,
        tensorArray(QStringLiteral("mut_shielding_diamagnetic"), 3));
    put(valid, h5reader::io::FieldKind::MutShieldingParamagnetic,
        tensorArray(QStringLiteral("mut_shielding_paramagnetic"), 3));
    put(valid, h5reader::io::FieldKind::DeltaShielding,
        tensorArray(QStringLiteral("delta_shielding"), 3));
    const auto delta = LoadDeltaRunData(valid, &err);
    QVERIFY(delta.has_value());
    QCOMPARE(delta->matched_axis_count, std::size_t(3));
    QCOMPARE(delta->wt_n, std::size_t(3));
    QCOMPARE(delta->ala_n, std::size_t(3));
}

QTEST_GUILESS_MAIN(RediscoverCohortContextTests)

#include "rediscover_cohort_context_tests.moc"
