// h5reader_model_tests — QtTest binary for h5-reader model-layer logic.
//
// One QObject test class, private slots discovered by QTEST_GUILESS_MAIN.
// Coverage: anchor variant equality + Ring-widening axis match; QtTopology
// ring-axis normalisation; DashboardSignalModel / DashboardPanelModel
// signal emission contracts (via QSignalSpy); TrajectorySignalCatalog
// descriptor auto-promotion + canSample contract.
//
// Data-driven tests use _data() slots + QFETCH for the ring-axis and
// anchor-equality tables. Emission tests use QSignalSpy.
//
// Build: linked by CMakeLists.txt under BUILD_TESTING AND Qt6Test_FOUND.
// Run:   ctest -R h5reader_model_tests --output-on-failure
//        ./h5reader_model_tests -functions
//        ./h5reader_model_tests testRingAxisNormalization

#include "model/DashboardPanelModel.h"
#include "model/DashboardSignal.h"
#include "model/DashboardSignalModel.h"
#include "model/QtRing.h"
#include "model/QtTopology.h"
#include "model/SignalDictionary.h"
#include "model/TrajectorySignalCatalog.h"

#include <QObject>
#include <QSignalSpy>
#include <QString>
#include <QUuid>
#include <QVector>
#include <QtTest>

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

using namespace h5reader::model;

// QFETCH columns and QSignalSpy argument extraction go through QVariant,
// which requires QMetaType registration for these non-Qt types. Local-only
// (header-side Q_DECLARE_METATYPE would force every TU to know about them).
Q_DECLARE_METATYPE(h5reader::model::SignalAnchor)
Q_DECLARE_METATYPE(h5reader::model::SignalBinding)
Q_DECLARE_METATYPE(h5reader::model::DashboardDisplayRef)

namespace {

DashboardSignal makeSignal(const QString& label, std::size_t atom = 0) {
    DashboardSignal signal;
    signal.label = label;
    signal.binding.sourceKind = SignalSourceKind::DenseH5Trajectory;
    signal.binding.descriptorId = QStringLiteral("h5:positions");
    signal.binding.conceptKey = QStringLiteral("positions");
    signal.binding.displayModeId = QStringLiteral("strip.vector.component");
    signal.binding.anchor = AtomAnchor{atom};
    signal.displayModeIds = {QStringLiteral("strip.vector.component")};
    signal.nativeAxis = SignalAxis::Atom;
    signal.requiredAnchor = SignalAxis::Atom;
    signal.valueShape = SignalValueShape::Vector3;
    return signal;
}

// One aromatic + one saturated ring so absoluteRingIndex has both subkinds
// to translate. The atom-count is 0 because these tests only exercise the
// ring-axis arithmetic; per-atom reverse-index caches aren't read here.
QtTopology makeRingTopology() {
    std::vector<std::unique_ptr<QtRing>> rings;
    auto aromatic = CreateQtRing(RingTypeIndex::PheBenzene);
    aromatic->ringId = 0;
    aromatic->nativeAxisIndex = 0;
    aromatic->ringKind = RingKind::Aromatic;
    rings.push_back(std::move(aromatic));

    auto saturated = CreateQtRing(RingTypeIndex::ProPyrrolidine);
    saturated->ringId = 1;
    saturated->nativeAxisIndex = 0;
    saturated->ringKind = RingKind::Saturated;
    rings.push_back(std::move(saturated));

    return QtTopology(0, {}, std::move(rings), {}, 1, 1);
}

DashboardDisplayRef makeRef(const QUuid& signalId,
                            const char* mode,
                            const char* channel) {
    return DashboardDisplayRef{signalId,
                               QString::fromLatin1(mode),
                               QString::fromLatin1(channel)};
}

}  // namespace

class DashboardModelTests : public QObject {
    Q_OBJECT

private slots:
    // ---- anchor variant + axis-matching ---------------------------------

    void testAnchorEquality_data();
    void testAnchorEquality();
    void testAnchorMatchesAxis_data();
    void testAnchorMatchesAxis();
    // Direct test of the axis-pair widening helper that both the
    // controller's AnchorMatchesAxis and the picker dialog's
    // anchorAxisCanSatisfy share. Lockstep contract; if the two ever
    // diverge again the dialog-side picker gap returns.
    void testAxisCanSatisfy_data();
    void testAxisCanSatisfy();

    // ---- ring axis normalisation ----------------------------------------

    void testRingAxisNormalization_data();
    void testRingAxisNormalization();

    // ---- DashboardSignalModel emission contracts ------------------------

    void testSignalModel_addEmitsSignalAdded();
    void testSignalModel_removeEmitsSignalRemoved();
    void testSignalModel_clearEmitsAllRemovedIds();
    void testSignalModel_findSignalById();

    // ---- DashboardPanelModel emission contracts -------------------------

    void testPanelModel_addDisplayRefEmitsChanges();
    void testPanelModel_removeDisplayRefForSignalCleansAllPanels();
    void testPanelModel_clearPreservesOnePanelAndEmitsRefs();
    void testPanelModel_signalReferenceCountTracksCorrectly();

    // ---- TrajectorySignalCatalog auto-promotion + canSample -------------

    void testCatalog_descriptorAutoPromotionForKnownPaths();
    void testCatalog_canSampleRequiresStripMode();
    void testCatalog_canSampleFalseForStaticMode();

    // ---- Per-TR catalog presence (one row per new TR landing) ----------

    void testCatalog_iredS2DescriptorPresent();
    void testCatalog_kernelDynamicsDescriptorsPresent();
    void testCatalog_reorientDescriptorsPresent();
    void testCatalog_dihedralAutocorrDescriptorsPresent();
    void testCatalog_kernelCoherenceDescriptorPresent();

    // ---- Phase H: lockstep regression ---------------------------------
    // Iterate every Valid temporal DenseH5Trajectory descriptor; assert
    // canSample(synthetic binding) returns true. Catches the "added a
    // descriptor but forgot to register its storagePath in kDensePaths
    // AND/OR forgot to add the denseH5Plan branch" failure mode.
    void testCatalog_allValidTemporalDenseH5DescriptorsAreSampleable();
};

// ---- anchor variant + axis-matching -------------------------------------

void DashboardModelTests::testAnchorEquality_data() {
    QTest::addColumn<SignalBinding>("a");
    QTest::addColumn<SignalBinding>("b");
    QTest::addColumn<bool>("expected");

    auto withAnchor = [](const SignalAnchor& anchor) {
        SignalBinding b;
        b.descriptorId = QStringLiteral("test");
        b.anchor = anchor;
        return b;
    };

    QTest::newRow("none-none")       << withAnchor(NoneAnchor{}) << withAnchor(NoneAnchor{}) << true;
    QTest::newRow("atom-equal")      << withAnchor(AtomAnchor{3}) << withAnchor(AtomAnchor{3}) << true;
    QTest::newRow("atom-different")  << withAnchor(AtomAnchor{3}) << withAnchor(AtomAnchor{4}) << false;
    QTest::newRow("atom-vs-residue") << withAnchor(AtomAnchor{3}) << withAnchor(ResidueAnchor{3}) << false;
    QTest::newRow("residue-equal")   << withAnchor(ResidueAnchor{7}) << withAnchor(ResidueAnchor{7}) << true;
    QTest::newRow("tuple-equal")     << withAnchor(AtomTupleAnchor{{1, 2, 3}}) << withAnchor(AtomTupleAnchor{{1, 2, 3}}) << true;
    QTest::newRow("tuple-order")     << withAnchor(AtomTupleAnchor{{1, 2, 3}}) << withAnchor(AtomTupleAnchor{{3, 2, 1}}) << false;
    QTest::newRow("bond-equal")      << withAnchor(BondAnchor{12}) << withAnchor(BondAnchor{12}) << true;
    QTest::newRow("bondvec-equal")    << withAnchor(BondVectorAnchor{7, 1}) << withAnchor(BondVectorAnchor{7, 1}) << true;
    QTest::newRow("bondvec-diff-res") << withAnchor(BondVectorAnchor{7, 1}) << withAnchor(BondVectorAnchor{8, 1}) << false;
    QTest::newRow("bondvec-diff-kind")<< withAnchor(BondVectorAnchor{7, 1}) << withAnchor(BondVectorAnchor{7, 2}) << false;
    QTest::newRow("bondvec-vs-bond")  << withAnchor(BondVectorAnchor{7, 1}) << withAnchor(BondAnchor{7}) << false;
    QTest::newRow("ring-equal")      << withAnchor(RingAnchor{2}) << withAnchor(RingAnchor{2}) << true;
    QTest::newRow("aromatic-equal")  << withAnchor(AromaticRingAnchor{0}) << withAnchor(AromaticRingAnchor{0}) << true;
    QTest::newRow("saturated-equal") << withAnchor(SaturatedRingAnchor{0}) << withAnchor(SaturatedRingAnchor{0}) << true;
    QTest::newRow("aromatic-vs-sat") << withAnchor(AromaticRingAnchor{0}) << withAnchor(SaturatedRingAnchor{0}) << false;
    QTest::newRow("ring-vs-aromatic")<< withAnchor(RingAnchor{0}) << withAnchor(AromaticRingAnchor{0}) << false;
    QTest::newRow("ringpair-equal")  << withAnchor(RingContributionPairAnchor{5}) << withAnchor(RingContributionPairAnchor{5}) << true;
    QTest::newRow("membership-equal")<< withAnchor(RingMembershipAnchor{8}) << withAnchor(RingMembershipAnchor{8}) << true;
    QTest::newRow("mutation-equal")  << withAnchor(MutationMatchPairAnchor{1}) << withAnchor(MutationMatchPairAnchor{1}) << true;
    QTest::newRow("protein-equal")   << withAnchor(ProteinAnchor{}) << withAnchor(ProteinAnchor{}) << true;
    QTest::newRow("system-equal")    << withAnchor(SystemAnchor{}) << withAnchor(SystemAnchor{}) << true;
    QTest::newRow("event-equal")     << withAnchor(EventAnchor{}) << withAnchor(EventAnchor{}) << true;
}

void DashboardModelTests::testAnchorEquality() {
    QFETCH(SignalBinding, a);
    QFETCH(SignalBinding, b);
    QFETCH(bool, expected);
    QCOMPARE(a == b, expected);
}

void DashboardModelTests::testAnchorMatchesAxis_data() {
    QTest::addColumn<SignalAnchor>("anchor");
    QTest::addColumn<int>("axisInt");
    QTest::addColumn<bool>("expected");

    // The Ring-widening rule: Aromatic and Saturated subkinds satisfy a
    // Ring-required descriptor. The reverse does NOT widen.
    QTest::newRow("none-axis-passes-all")    << SignalAnchor{NoneAnchor{}} << int(SignalAxis::None) << true;
    QTest::newRow("atom-axis-match")         << SignalAnchor{AtomAnchor{0}} << int(SignalAxis::Atom) << true;
    QTest::newRow("residue-axis-match")      << SignalAnchor{ResidueAnchor{0}} << int(SignalAxis::Residue) << true;
    QTest::newRow("aromatic-satisfies-ring") << SignalAnchor{AromaticRingAnchor{0}} << int(SignalAxis::Ring) << true;
    QTest::newRow("saturated-satisfies-ring")<< SignalAnchor{SaturatedRingAnchor{0}} << int(SignalAxis::Ring) << true;
    QTest::newRow("ring-does-not-satisfy-aromatic") << SignalAnchor{RingAnchor{0}} << int(SignalAxis::AromaticRing) << false;
    QTest::newRow("atom-mismatch-residue")   << SignalAnchor{AtomAnchor{0}} << int(SignalAxis::Residue) << false;
    // The BondVector-widening rule: a Residue anchor satisfies a
    // BondVector-required descriptor (picking residue Lys17 surfaces its
    // N-H / Cα-Hα / C=O candidates). The reverse does NOT widen.
    QTest::newRow("bondvec-axis-match")       << SignalAnchor{BondVectorAnchor{0, 1}} << int(SignalAxis::BondVector) << true;
    QTest::newRow("residue-satisfies-bondvec")<< SignalAnchor{ResidueAnchor{0}} << int(SignalAxis::BondVector) << true;
    QTest::newRow("bondvec-no-residue")       << SignalAnchor{BondVectorAnchor{0, 1}} << int(SignalAxis::Residue) << false;
    QTest::newRow("atom-no-bondvec")          << SignalAnchor{AtomAnchor{0}} << int(SignalAxis::BondVector) << false;
}

void DashboardModelTests::testAnchorMatchesAxis() {
    QFETCH(SignalAnchor, anchor);
    QFETCH(int, axisInt);
    QFETCH(bool, expected);
    const auto axis = static_cast<SignalAxis>(axisInt);
    QCOMPARE(AnchorMatchesAxis(anchor, axis), expected);
}

void DashboardModelTests::testAxisCanSatisfy_data() {
    QTest::addColumn<int>("selectedInt");
    QTest::addColumn<int>("requiredInt");
    QTest::addColumn<bool>("expected");

    // None on the required side is always permissive (filter-off case).
    QTest::newRow("none-required-accepts-any-atom") << int(SignalAxis::Atom) << int(SignalAxis::None) << true;

    // Identity matches.
    QTest::newRow("atom-atom") << int(SignalAxis::Atom) << int(SignalAxis::Atom) << true;
    QTest::newRow("bondvec-bondvec") << int(SignalAxis::BondVector) << int(SignalAxis::BondVector) << true;

    // Ring widening (existing).
    QTest::newRow("aromatic-satisfies-ring") << int(SignalAxis::AromaticRing) << int(SignalAxis::Ring) << true;
    QTest::newRow("saturated-satisfies-ring") << int(SignalAxis::SaturatedRing) << int(SignalAxis::Ring) << true;

    // BondVector widening (the dialog-side gap that Codex NOW-3 caught).
    QTest::newRow("residue-satisfies-bondvec") << int(SignalAxis::Residue) << int(SignalAxis::BondVector) << true;

    // Negative cases — neither widening flips the other direction.
    QTest::newRow("ring-no-aromatic") << int(SignalAxis::Ring) << int(SignalAxis::AromaticRing) << false;
    QTest::newRow("bondvec-no-residue") << int(SignalAxis::BondVector) << int(SignalAxis::Residue) << false;
    QTest::newRow("atom-no-bondvec") << int(SignalAxis::Atom) << int(SignalAxis::BondVector) << false;
}

void DashboardModelTests::testAxisCanSatisfy() {
    QFETCH(int, selectedInt);
    QFETCH(int, requiredInt);
    QFETCH(bool, expected);
    QCOMPARE(AxisCanSatisfy(static_cast<SignalAxis>(selectedInt),
                            static_cast<SignalAxis>(requiredInt)),
             expected);
}

// ---- ring axis normalisation --------------------------------------------

void DashboardModelTests::testRingAxisNormalization_data() {
    QTest::addColumn<int>("axisInt");
    QTest::addColumn<int>("subIndex");
    QTest::addColumn<int>("expected");  // -1 = nullopt

    // Layout: aromatic at row 0, saturated at row 1.
    QTest::newRow("ring-absolute-0")    << int(QtRingAxis::Ring) << 0 << 0;
    QTest::newRow("ring-absolute-1")    << int(QtRingAxis::Ring) << 1 << 1;
    QTest::newRow("ring-overflow")      << int(QtRingAxis::Ring) << 2 << -1;
    QTest::newRow("aromatic-0")         << int(QtRingAxis::AromaticRing) << 0 << 0;
    QTest::newRow("aromatic-overflow")  << int(QtRingAxis::AromaticRing) << 1 << -1;
    QTest::newRow("saturated-0-offset") << int(QtRingAxis::SaturatedRing) << 0 << 1;
    QTest::newRow("saturated-overflow") << int(QtRingAxis::SaturatedRing) << 1 << -1;
}

void DashboardModelTests::testRingAxisNormalization() {
    QFETCH(int, axisInt);
    QFETCH(int, subIndex);
    QFETCH(int, expected);
    const QtTopology topology = makeRingTopology();
    const auto axis = static_cast<QtRingAxis>(axisInt);
    const auto result = topology.absoluteRingIndex(axis, static_cast<std::size_t>(subIndex));
    if (expected < 0) {
        QVERIFY(!result.has_value());
    } else {
        QVERIFY(result.has_value());
        QCOMPARE(*result, static_cast<std::size_t>(expected));
    }
}

// ---- DashboardSignalModel emission contracts ----------------------------

void DashboardModelTests::testSignalModel_addEmitsSignalAdded() {
    DashboardSignalModel model;
    QSignalSpy spy(&model, &DashboardSignalModel::signalAdded);
    QVERIFY(spy.isValid());

    const QUuid id = model.addSignal(makeSignal(QStringLiteral("first")));
    QVERIFY(!id.isNull());
    QCOMPARE(spy.count(), 1);
    QCOMPARE(spy.first().at(0).toUuid(), id);
    QCOMPARE(model.rowCount(), 1);
}

void DashboardModelTests::testSignalModel_removeEmitsSignalRemoved() {
    DashboardSignalModel model;
    const QUuid id = model.addSignal(makeSignal(QStringLiteral("doomed")));

    QSignalSpy spy(&model, &DashboardSignalModel::signalRemoved);
    QVERIFY(spy.isValid());
    QVERIFY(model.removeSignal(id));
    QCOMPARE(spy.count(), 1);
    QCOMPARE(spy.first().at(0).toUuid(), id);
    QCOMPARE(model.rowCount(), 0);

    // Removing a non-existent id is a no-op (returns false, no emit).
    QVERIFY(!model.removeSignal(id));
    QCOMPARE(spy.count(), 1);
}

void DashboardModelTests::testSignalModel_clearEmitsAllRemovedIds() {
    DashboardSignalModel model;
    const QUuid first  = model.addSignal(makeSignal(QStringLiteral("first")));
    const QUuid second = model.addSignal(makeSignal(QStringLiteral("second"), /*atom=*/1));
    const QUuid third  = model.addSignal(makeSignal(QStringLiteral("third"),  /*atom=*/2));

    QSignalSpy spy(&model, &DashboardSignalModel::signalRemoved);
    QVERIFY(spy.isValid());
    model.clear();

    QCOMPARE(model.rowCount(), 0);
    QCOMPARE(spy.count(), 3);

    QVector<QUuid> emitted;
    for (const QList<QVariant>& call : spy)
        emitted.push_back(call.at(0).toUuid());
    QVERIFY(emitted.contains(first));
    QVERIFY(emitted.contains(second));
    QVERIFY(emitted.contains(third));

    // clear() on empty model is a no-op.
    QSignalSpy spyEmpty(&model, &DashboardSignalModel::signalRemoved);
    model.clear();
    QCOMPARE(spyEmpty.count(), 0);
}

void DashboardModelTests::testSignalModel_findSignalById() {
    DashboardSignalModel model;
    const QUuid id = model.addSignal(makeSignal(QStringLiteral("findable")));
    const DashboardSignal* found = model.signalById(id);
    QVERIFY(found != nullptr);
    QCOMPARE(found->id, id);
    QCOMPARE(found->label, QStringLiteral("findable"));

    QVERIFY(model.signalById(QUuid::createUuid()) == nullptr);
}

// ---- DashboardPanelModel emission contracts -----------------------------

void DashboardModelTests::testPanelModel_addDisplayRefEmitsChanges() {
    DashboardPanelModel panels;
    const QUuid panelId = panels.activePanelId();
    QVERIFY(!panelId.isNull());

    const QUuid signalId = QUuid::createUuid();
    const DashboardDisplayRef ref = makeRef(signalId, "strip.scalar", "value");

    QSignalSpy refsChangedSpy(&panels, &DashboardPanelModel::displayRefsChanged);
    QVERIFY(refsChangedSpy.isValid());
    QVERIFY(panels.addDisplayRef(panelId, ref));
    QCOMPARE(refsChangedSpy.count(), 1);
    QCOMPARE(panels.signalReferenceCount(signalId), 1);

    // Adding the same ref again is rejected, no emission.
    QVERIFY(!panels.addDisplayRef(panelId, ref));
    QCOMPARE(refsChangedSpy.count(), 1);
}

void DashboardModelTests::testPanelModel_removeDisplayRefForSignalCleansAllPanels() {
    DashboardPanelModel panels;
    const QUuid first  = panels.activePanelId();
    const QUuid second = panels.addPanel(QStringLiteral("second"));

    const QUuid signalA = QUuid::createUuid();
    const QUuid signalB = QUuid::createUuid();
    QVERIFY(panels.addDisplayRef(first,  makeRef(signalA, "strip.scalar", "value")));
    QVERIFY(panels.addDisplayRef(second, makeRef(signalA, "strip.rollup", "mean")));
    QVERIFY(panels.addDisplayRef(second, makeRef(signalB, "strip.scalar", "value")));

    QSignalSpy removedSpy(&panels, &DashboardPanelModel::displayRefRemoved);
    QVERIFY(removedSpy.isValid());

    const int removed = panels.removeDisplayRefsForSignal(signalA);
    QCOMPARE(removed, 2);
    QCOMPARE(panels.signalReferenceCount(signalA), 0);
    QCOMPARE(panels.signalReferenceCount(signalB), 1);
    QCOMPARE(removedSpy.count(), 2);
}

void DashboardModelTests::testPanelModel_clearPreservesOnePanelAndEmitsRefs() {
    DashboardPanelModel panels;
    const QUuid firstPanel  = panels.activePanelId();
    const QUuid secondPanel = panels.addPanel(QStringLiteral("second"));
    const QUuid signalId    = QUuid::createUuid();

    const DashboardDisplayRef refA = makeRef(signalId, "strip.scalar", "value");
    const DashboardDisplayRef refB = makeRef(signalId, "strip.rollup", "mean");
    QVERIFY(panels.addDisplayRef(firstPanel,  refA));
    QVERIFY(panels.addDisplayRef(secondPanel, refB));

    QSignalSpy removedSpy(&panels, &DashboardPanelModel::displayRefRemoved);
    QVERIFY(removedSpy.isValid());

    panels.clear();

    QCOMPARE(panels.rowCount(), 1);
    QVERIFY(!panels.activePanelId().isNull());
    QCOMPARE(panels.signalReferenceCount(signalId), 0);
    QCOMPARE(removedSpy.count(), 2);

    QVector<DashboardDisplayRef> emitted;
    for (const QList<QVariant>& call : removedSpy)
        emitted.push_back(call.at(1).value<DashboardDisplayRef>());
    QVERIFY(emitted.contains(refA));
    QVERIFY(emitted.contains(refB));
}

void DashboardModelTests::testPanelModel_signalReferenceCountTracksCorrectly() {
    DashboardPanelModel panels;
    const QUuid panelId = panels.activePanelId();
    const QUuid signalId = QUuid::createUuid();

    QCOMPARE(panels.signalReferenceCount(signalId), 0);
    QVERIFY(panels.addDisplayRef(panelId, makeRef(signalId, "strip.scalar", "value")));
    QCOMPARE(panels.signalReferenceCount(signalId), 1);
    QVERIFY(panels.addDisplayRef(panelId, makeRef(signalId, "strip.rollup", "mean")));
    QCOMPARE(panels.signalReferenceCount(signalId), 2);

    QVERIFY(panels.removeDisplayRef(panelId, makeRef(signalId, "strip.scalar", "value")));
    QCOMPARE(panels.signalReferenceCount(signalId), 1);
}

// ---- TrajectorySignalCatalog auto-promotion + canSample -----------------

void DashboardModelTests::testCatalog_descriptorAutoPromotionForKnownPaths() {
    const TrajectorySignalCatalog catalog;

    // h5:positions is in the hasImplementedTemporalSampler allowlist —
    // makeDescriptor flips it from default Gap/Pending to Valid/None.
    const SignalDescriptor* positions = catalog.findDescriptor(QStringLiteral("h5:positions"));
    QVERIFY(positions != nullptr);
    QCOMPARE(positions->samplingStatus, SampleStatus::Valid);
    QCOMPARE(positions->samplingGapReason, GapReason::None);

    // topology:atoms is non-temporal — flips to NotAvailable/NotApplicable
    // so the picker shows it honestly even though it can't be sampled as a strip.
    const SignalDescriptor* atoms = catalog.findDescriptor(QStringLiteral("topology:atoms"));
    QVERIFY(atoms != nullptr);
    QCOMPARE(atoms->samplingStatus, SampleStatus::NotAvailable);
    QCOMPARE(atoms->samplingGapReason, GapReason::NotApplicable);
}

void DashboardModelTests::testCatalog_canSampleRequiresStripMode() {
    const TrajectorySignalCatalog catalog;

    DisplaySignalBinding stripBinding;
    stripBinding.sourceKind = SignalSourceKind::DenseH5Trajectory;
    stripBinding.descriptorId = QStringLiteral("h5:positions");
    stripBinding.conceptKey = QStringLiteral("positions");
    stripBinding.displayModeId = QStringLiteral("strip.vector.component");
    stripBinding.anchor = AtomAnchor{0};

    QVERIFY(catalog.canSample(stripBinding));
}

void DashboardModelTests::testCatalog_iredS2DescriptorPresent() {
    const TrajectorySignalCatalog catalog;

    const SignalDescriptor* ired = catalog.findDescriptor(QStringLiteral("h5:ired_s2"));
    QVERIFY(ired != nullptr);
    QCOMPARE(ired->sourceKind, SignalSourceKind::DenseH5Trajectory);
    QCOMPARE(ired->nativeAxis, SignalAxis::BondVector);
    QCOMPARE(ired->requiredAnchor, SignalAxis::BondVector);
    QCOMPARE(ired->valueShape, SignalValueShape::Scalar);
    QCOMPARE(ired->storagePath, QStringLiteral("/trajectory/ired_order_parameters"));
    QCOMPARE(ired->samplingStatus, SampleStatus::Valid);
    QCOMPARE(ired->samplingGapReason, GapReason::None);
    QVERIFY(ired->staticModes.contains(QStringLiteral("static.bar.sequence")));
}

void DashboardModelTests::testCatalog_kernelDynamicsDescriptorsPresent() {
    const TrajectorySignalCatalog catalog;

    // ACF curve descriptor
    const SignalDescriptor* acf = catalog.findDescriptor(QStringLiteral("h5:kernel_dynamics_acf"));
    QVERIFY(acf != nullptr);
    QCOMPARE(acf->valueShape, SignalValueShape::CurveOverLag);
    QCOMPARE(acf->storagePath, QStringLiteral("/trajectory/kernel_dynamics"));
    QVERIFY(acf->staticModes.contains(QStringLiteral("static.curve.lag.animated")));
    QCOMPARE(acf->channels.size(), 13);
    QCOMPARE(acf->channels.first().id, QStringLiteral("bs_T0"));

    // PSD curve descriptor
    const SignalDescriptor* psd = catalog.findDescriptor(QStringLiteral("h5:kernel_dynamics_psd"));
    QVERIFY(psd != nullptr);
    QCOMPARE(psd->valueShape, SignalValueShape::Spectrum);
    QVERIFY(psd->staticModes.contains(QStringLiteral("static.spectrum.power")));

    // Three scalar reductions (per-class block on atom axis)
    for (const char* id : {"h5:kernel_dynamics_decay_time",
                           "h5:kernel_dynamics_peak_freq",
                           "h5:kernel_dynamics_spectral_centroid"}) {
        const SignalDescriptor* scalar = catalog.findDescriptor(QString::fromLatin1(id));
        QVERIFY2(scalar != nullptr, id);
        QCOMPARE(scalar->valueShape, SignalValueShape::PerClassBlock);
        QCOMPARE(scalar->channels.size(), 13);
    }
}

void DashboardModelTests::testCatalog_reorientDescriptorsPresent() {
    const TrajectorySignalCatalog catalog;
    // Five scalar descriptors all on BondVector axis.
    for (const char* id : {"h5:reorient_s2", "h5:reorient_tau_e",
                           "h5:reorient_r1", "h5:reorient_r2", "h5:reorient_noe"}) {
        const SignalDescriptor* d = catalog.findDescriptor(QString::fromLatin1(id));
        QVERIFY2(d != nullptr, id);
        QCOMPARE(d->nativeAxis, SignalAxis::BondVector);
        QCOMPARE(d->valueShape, SignalValueShape::Scalar);
        QCOMPARE(d->storagePath, QStringLiteral("/trajectory/reorientational_dynamics"));
        QVERIFY(d->staticModes.contains(QStringLiteral("static.bar.sequence")));
    }
    // Two TCF curve descriptors.
    for (const char* id : {"h5:reorient_acf_internal", "h5:reorient_acf_lab"}) {
        const SignalDescriptor* d = catalog.findDescriptor(QString::fromLatin1(id));
        QVERIFY2(d != nullptr, id);
        QCOMPARE(d->valueShape, SignalValueShape::CurveOverLag);
        QVERIFY(d->staticModes.contains(QStringLiteral("static.curve.lag.animated")));
    }
    // L-3a: Mat3 orientation-tensor descriptor (ellipsoid glyph in
    // the 3-D scene). Carries static.tensor as its primary mode.
    const SignalDescriptor* tensor = catalog.findDescriptor(QStringLiteral("h5:reorient_orientation_tensor"));
    QVERIFY(tensor != nullptr);
    QCOMPARE(tensor->valueShape, SignalValueShape::Mat3PerRow);
    QVERIFY(tensor->staticModes.contains(QStringLiteral("static.tensor")));
    // L-3b: FixedFreqBlock J(ω) descriptor (5 KTB Larmor combinations,
    // NH only). Carries static.fixed_freq for the dedicated panel.
    const SignalDescriptor* j = catalog.findDescriptor(QStringLiteral("h5:reorient_spectral_density"));
    QVERIFY(j != nullptr);
    QCOMPARE(j->valueShape, SignalValueShape::FixedFreqBlock);
    QVERIFY(j->staticModes.contains(QStringLiteral("static.fixed_freq")));
}

void DashboardModelTests::testCatalog_dihedralAutocorrDescriptorsPresent() {
    const TrajectorySignalCatalog catalog;
    for (const char* id : {"h5:dihedral_phi_corr_time", "h5:dihedral_psi_corr_time"}) {
        const SignalDescriptor* d = catalog.findDescriptor(QString::fromLatin1(id));
        QVERIFY2(d != nullptr, id);
        QCOMPARE(d->nativeAxis, SignalAxis::Residue);
        QCOMPARE(d->valueShape, SignalValueShape::Scalar);
        QCOMPARE(d->storagePath, QStringLiteral("/trajectory/dihedral_autocorrelation"));
    }
    for (const char* id : {"h5:dihedral_phi_acf", "h5:dihedral_psi_acf"}) {
        const SignalDescriptor* d = catalog.findDescriptor(QString::fromLatin1(id));
        QVERIFY2(d != nullptr, id);
        QCOMPARE(d->valueShape, SignalValueShape::CurveOverLag);
    }
    // L-2a chi composite descriptors — PerClassBlock scalar + CurveOverLag
    // with 4 chi channels (chi0..chi3). Per-channel dispatch lives in the
    // controller's panel builders.
    const SignalDescriptor* chiScalar = catalog.findDescriptor(QStringLiteral("h5:dihedral_chi_corr_time"));
    QVERIFY(chiScalar != nullptr);
    QCOMPARE(chiScalar->valueShape, SignalValueShape::PerClassBlock);
    QCOMPARE(chiScalar->storagePath, QStringLiteral("/trajectory/dihedral_autocorrelation"));
    QCOMPARE(chiScalar->channels.size(), 4);
    QCOMPARE(chiScalar->channels.first().id, QStringLiteral("chi0"));
    QVERIFY(chiScalar->staticModes.contains(QStringLiteral("static.bar.sequence")));

    const SignalDescriptor* chiAcf = catalog.findDescriptor(QStringLiteral("h5:dihedral_chi_acf"));
    QVERIFY(chiAcf != nullptr);
    QCOMPARE(chiAcf->valueShape, SignalValueShape::CurveOverLag);
    QCOMPARE(chiAcf->channels.size(), 4);
    QVERIFY(chiAcf->staticModes.contains(QStringLiteral("static.curve.lag.animated")));
}

void DashboardModelTests::testCatalog_kernelCoherenceDescriptorPresent() {
    const TrajectorySignalCatalog catalog;
    const SignalDescriptor* d = catalog.findDescriptor(QStringLiteral("h5:kernel_coherence"));
    QVERIFY(d != nullptr);
    QCOMPARE(d->nativeAxis, SignalAxis::Atom);
    QCOMPARE(d->valueShape, SignalValueShape::Matrix);
    QCOMPARE(d->storagePath, QStringLiteral("/trajectory/kernel_coherence"));
    QVERIFY(d->staticModes.contains(QStringLiteral("static.chord.coupling")));
    QCOMPARE(d->channels.size(), 13);
}

void DashboardModelTests::testCatalog_allValidTemporalDenseH5DescriptorsAreSampleable() {
    const TrajectorySignalCatalog catalog;
    int checked = 0;
    int failed = 0;
    QStringList failureMessages;
    for (const SignalDescriptor& d : catalog.descriptorList()) {
        if (d.sourceKind != SignalSourceKind::DenseH5Trajectory) continue;
        if (d.samplingStatus != SampleStatus::Valid) continue;
        if (d.temporalModes.isEmpty()) continue;

        DisplaySignalBinding binding;
        binding.sourceKind = d.sourceKind;
        binding.descriptorId = d.id;
        binding.conceptKey = d.conceptKey;
        binding.displayModeId = d.temporalModes.first();
        // Anchor matching the descriptor's required axis. Synthetic row 0.
        switch (d.requiredAnchor) {
        case SignalAxis::Atom:    binding.anchor = AtomAnchor{0}; break;
        case SignalAxis::Residue: binding.anchor = ResidueAnchor{0}; break;
        case SignalAxis::Bond:    binding.anchor = BondAnchor{0}; break;
        case SignalAxis::Ring:    binding.anchor = RingAnchor{0}; break;
        case SignalAxis::AromaticRing: binding.anchor = AromaticRingAnchor{0}; break;
        case SignalAxis::SaturatedRing: binding.anchor = SaturatedRingAnchor{0}; break;
        case SignalAxis::RingMembership: binding.anchor = RingMembershipAnchor{0}; break;
        case SignalAxis::RingContributionPair: binding.anchor = RingContributionPairAnchor{0}; break;
        case SignalAxis::MutationMatchPair: binding.anchor = MutationMatchPairAnchor{0}; break;
        case SignalAxis::Protein: binding.anchor = ProteinAnchor{}; break;
        case SignalAxis::System:  binding.anchor = SystemAnchor{}; break;
        case SignalAxis::Event:   binding.anchor = EventAnchor{}; break;
        case SignalAxis::AtomTuple: binding.anchor = AtomTupleAnchor{{0}}; break;
        case SignalAxis::BondVector: binding.anchor = BondVectorAnchor{0, 1}; break;
        case SignalAxis::None:    binding.anchor = NoneAnchor{}; break;
        }
        ++checked;
        if (!catalog.canSample(binding)) {
            ++failed;
            failureMessages << QStringLiteral("  %1 (path=%2, mode=%3)")
                                   .arg(d.id, d.storagePath, binding.displayModeId);
        }
    }
    if (failed > 0) {
        const QString details = failureMessages.join(QStringLiteral("\n"));
        QFAIL(qPrintable(QStringLiteral(
            "Lockstep: %1 of %2 Valid temporal DenseH5Trajectory descriptors are not "
            "sampleable. Likely missing kDensePaths entry or denseH5Plan branch.\n%3")
                                .arg(failed).arg(checked).arg(details)));
    }
    QVERIFY(checked > 0);  // Sanity: we actually iterated some descriptors.
}

void DashboardModelTests::testCatalog_canSampleFalseForStaticMode() {
    const TrajectorySignalCatalog catalog;

    DisplaySignalBinding staticBinding;
    staticBinding.sourceKind = SignalSourceKind::Topology;
    staticBinding.descriptorId = QStringLiteral("topology:atoms");
    staticBinding.conceptKey = QStringLiteral("topology.atoms");
    staticBinding.displayModeId = QStringLiteral("static.table");
    staticBinding.anchor = AtomAnchor{0};

    QVERIFY(!catalog.canSample(staticBinding));
}

QTEST_GUILESS_MAIN(DashboardModelTests)

#include "dashboard_model_tests.moc"
