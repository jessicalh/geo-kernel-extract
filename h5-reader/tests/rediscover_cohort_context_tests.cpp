#include "rediscover/CohortContextAccumulator.h"
#include "rediscover/DeltaRunData.h"
#include "rediscover/ScopedProducerCatalog.h"

#include <QtTest/QtTest>

#include <QFile>

#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#ifdef __linux__
#include <unistd.h>
#endif

using namespace h5reader::rediscover;
namespace model = h5reader::model;

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

std::size_t residentSetBytes() {
#ifdef __linux__
    QFile statm(QStringLiteral("/proc/self/statm"));
    if (!statm.open(QIODevice::ReadOnly | QIODevice::Text)) return 0;
    const QList<QByteArray> fields = statm.readAll().trimmed().split(' ');
    if (fields.size() < 2) return 0;
    bool ok = false;
    const qulonglong residentPages = fields[1].toULongLong(&ok);
    if (!ok) return 0;
    const long pageSize = sysconf(_SC_PAGESIZE);
    if (pageSize <= 0) return 0;
    return static_cast<std::size_t>(residentPages) * static_cast<std::size_t>(pageSize);
#else
    return 0;
#endif
}

}  // namespace

class RediscoverCohortContextTests : public QObject {
    Q_OBJECT

private slots:
    void keyBuilderUsesFullIupacIdentity();
    void supportCredentialDoesNotAuthorizeInsufficientCouplings();
    void accumulatorCountsProteinResidency();
    void accumulatorDoesNotRetainSamples();
    void accumulatorResidentSetIsFlatInSampleCount();
    void permutationNullUsesProteinLabelShuffle();
    void permutationNullLeavesSentinelBlankWhenNotComputable();
    void helixDipoleHasExpectedSign();
    void deltaLoaderIsNarrowAndRefusesMissingArrays();
    void foldAdapterProjectsMolecularComponents();
    void classicalAgreementUsesPerProteinPoints();
    void syntheticFoldedDeltaRidgeMatchesHandSlope();
    void syntheticDistantNonzeroRidgeIsCharacterizedNotGated();
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

void RediscoverCohortContextTests::accumulatorDoesNotRetainSamples() {
    Axis2ContextKeyFields f;
    f.element = QStringLiteral("C");
    f.residue_type = QStringLiteral("ALA");
    f.atom_name = QStringLiteral("CA");
    f.hyb = QStringLiteral("sp3");
    f.contact_class = QStringLiteral("hydrophobic");
    f.dihedral_region = QStringLiteral("AlphaR");
    f.SS = QStringLiteral("helix");

    CohortContextAccumulator acc;
    for (int i = 0; i < 1000; ++i) {
        CohortSample s;
        s.key = BuildAxis2ContextKey(f);
        s.protein_id = QStringLiteral("p");
        s.sigma_iso = static_cast<double>(i);
        s.channels.insert(QStringLiteral("apbs_E_mag"), static_cast<double>(i % 7));
        acc.push(s);
    }
    const CohortCellTruth& cell = acc.cells().begin()->second;
    QCOMPARE(acc.sampleCount(), std::size_t(1000));
    QCOMPARE(cell.sample_count, std::size_t(1000));
    QCOMPARE(cell.retainedSampleCount(), std::size_t(0));
    QVERIFY(cell.retainedAccumulatorValueCount() < cell.sample_count);
}

void RediscoverCohortContextTests::accumulatorResidentSetIsFlatInSampleCount() {
#ifndef __linux__
    QSKIP("RSS flatness check uses /proc/self/statm on Linux");
#else
    Axis2ContextKeyFields f;
    f.element = QStringLiteral("C");
    f.residue_type = QStringLiteral("ALA");
    f.atom_name = QStringLiteral("CA");
    f.hyb = QStringLiteral("sp3");
    f.contact_class = QStringLiteral("hydrophobic");
    f.dihedral_region = QStringLiteral("AlphaR");
    f.SS = QStringLiteral("helix");
    const Axis2ContextKey key = BuildAxis2ContextKey(f);

    CohortContextAccumulator acc;
    auto pushSample = [&](int i) {
        CohortSample s;
        s.key = key;
        s.protein_id = QStringLiteral("p");
        s.sigma_iso = static_cast<double>(i % 97);
        s.channels.insert(QStringLiteral("apbs_E_mag"), static_cast<double>(i % 13));
        s.channels.insert(QStringLiteral("ring_bs_iso"), static_cast<double>(i % 17));
        acc.push(s);
    };

    for (int i = 0; i < 10000; ++i) pushSample(i);
    const std::size_t baselineBytes = residentSetBytes();
    QVERIFY(baselineBytes > 0);

    for (int i = 10000; i < 160000; ++i) pushSample(i);
    const std::size_t afterBytes = residentSetBytes();
    QVERIFY(afterBytes > 0);

    const CohortCellTruth& cell = acc.cells().begin()->second;
    QCOMPARE(cell.sample_count, std::size_t(160000));
    QCOMPARE(cell.retainedSampleCount(), std::size_t(0));
    QVERIFY(cell.retainedAccumulatorValueCount() < 2048);

    const std::size_t rssDelta = afterBytes > baselineBytes ? afterBytes - baselineBytes : 0;
    QVERIFY2(rssDelta < 16u * 1024u * 1024u,
             qPrintable(QStringLiteral("RSS grew by %1 bytes after 150000 extra same-cell samples")
                            .arg(rssDelta)));
#endif
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

void RediscoverCohortContextTests::permutationNullLeavesSentinelBlankWhenNotComputable() {
    const PermutationNull p = ProteinLabelPermutationNull({1.0, 2.0}, {3.0, 4.0}, 31, 123u);
    QCOMPARE(p.permutation_K, 31);
    QVERIFY(std::isnan(p.perm_p));
    QVERIFY(std::isnan(p.null_slope_mean));
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

void RediscoverCohortContextTests::foldAdapterProjectsMolecularComponents() {
    model::Mat3 raw = model::Mat3::Zero();
    raw(0, 0) = 1.0;
    raw(1, 1) = 2.0;
    raw(2, 2) = 3.0;

    model::Mat3 axes = model::Mat3::Zero();
    axes.col(0) = model::Vec3::UnitY();
    axes.col(1) = model::Vec3::UnitX();
    axes.col(2) = model::Vec3::UnitZ();

    const Axis2FoldedTensor folded = FoldAxis2TensorChannels(raw, std::optional<model::Mat3>(axes));
    QVERIFY(folded.molecular_frame_projected);
    QCOMPARE(folded.sigma_iso, 2.0);
    QCOMPARE(folded.mol_components[0], 2.0);
    QCOMPARE(folded.mol_components[1], 1.0);
    QVERIFY(folded.projection.contains(QStringLiteral("R^T*T*R")));
}

void RediscoverCohortContextTests::classicalAgreementUsesPerProteinPoints() {
    Axis2ContextKeyFields f;
    f.element = QStringLiteral("H");
    f.residue_type = QStringLiteral("GLY");
    f.atom_name = QStringLiteral("H");
    f.hyb = QStringLiteral("s");
    f.contact_class = QStringLiteral("polar");
    f.dihedral_region = QStringLiteral("AlphaR");
    f.SS = QStringLiteral("helix");

    CohortContextAccumulator acc;
    for (int i = 0; i < 4; ++i) {
        CohortSample s;
        s.key = BuildAxis2ContextKey(f);
        s.protein_id = QStringLiteral("p%1").arg(i);
        s.sigma_iso = 10.0 + 2.0 * i;
        s.channels.insert(QStringLiteral("apbs_E_mag"), static_cast<double>(i));
        s.channels.insert(QStringLiteral("ring_bs_iso"), 10.0 + static_cast<double>(i));
        s.channels.insert(QStringLiteral("mc_lit_iso"), static_cast<double>(i));
        acc.push(s);
    }

    const ClassicalAgreementStats stats =
        ComputeClassicalAgreementForCell(acc.cells().begin()->second, 0.0);
    QVERIFY(std::isfinite(stats.r));
    QVERIFY(stats.r > 0.999);
    QCOMPARE(stats.slope, 1.0);
    QCOMPARE(stats.rmsd, 0.0);
    QVERIFY(std::isfinite(stats.residual_sd));
}

void RediscoverCohortContextTests::syntheticFoldedDeltaRidgeMatchesHandSlope() {
    model::Mat3 axes = model::Mat3::Zero();
    axes.col(0) = model::Vec3::UnitY();
    axes.col(1) = model::Vec3::UnitX();
    axes.col(2) = model::Vec3::UnitZ();

    std::vector<double> deltaChannel;
    std::vector<double> deltaSigma;
    for (int i = 1; i <= 4; ++i) {
        model::Mat3 wt = model::Mat3::Zero();
        model::Mat3 ala = model::Mat3::Zero();
        wt(0, 0) = static_cast<double>(i);
        wt(1, 1) = 2.0 * static_cast<double>(i);
        wt(2, 2) = 3.0 * static_cast<double>(i);
        ala(0, 0) = 0.0;
        ala(1, 1) = 0.0;
        ala(2, 2) = 0.0;
        const Axis2FoldedTensor wtFold = FoldAxis2TensorChannels(wt, std::optional<model::Mat3>(axes));
        const Axis2FoldedTensor alaFold = FoldAxis2TensorChannels(ala, std::optional<model::Mat3>(axes));
        deltaChannel.push_back(wtFold.mol_components[0] - alaFold.mol_components[0]);
        deltaSigma.push_back(2.0 * (wtFold.sigma_iso - alaFold.sigma_iso));
    }
    QCOMPARE(LinearSlope(deltaChannel, deltaSigma), 2.0);
}

void RediscoverCohortContextTests::syntheticDistantNonzeroRidgeIsCharacterizedNotGated() {
    const DistantRidgeCharacterization flagged =
        CharacterizeDistantNonzeroRidge(3,
                                        QStringLiteral("distant_from_all_sites"),
                                        0.25,
                                        QStringLiteral("delta_mol_xx_projected_refolded"));
    QVERIFY(flagged.flagged);
    QCOMPARE(flagged.distant_zero_check,
             QStringLiteral("flagged_nonzero_distant_from_all_sites"));
    QCOMPARE(flagged.characterization, QStringLiteral("characterized_not_gated"));
    QCOMPARE(flagged.nonzero_channel, QStringLiteral("delta_mol_xx_projected_refolded"));

    const DistantRidgeCharacterization near =
        CharacterizeDistantNonzeroRidge(3,
                                        QStringLiteral("near_mutation_site"),
                                        0.25,
                                        QStringLiteral("delta_mol_xx_projected_refolded"));
    QVERIFY(!near.flagged);
    QCOMPARE(near.distant_zero_check, QStringLiteral("ok"));
    QVERIFY(near.characterization.isEmpty());
}

QTEST_GUILESS_MAIN(RediscoverCohortContextTests)

#include "rediscover_cohort_context_tests.moc"
