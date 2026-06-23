// display_policy_tests -- headless QtTest for the PURE per-field presentation
// policy (model/DisplayPolicy.h). A2-step-1: the displayability rule that makes
// the structural topology tables + the 256-d embedding honestly non-displayable
// (they were offered dead/nonsense modes by the mechanical-by-shape assignment),
// while keeping the real per-bond bond_length scalar displayable. Links Qt6::Test
// only (header-only, like the other policy/taxonomy tests).

#include "model/DisplayPolicy.h"

#include <QObject>
#include <QtTest>

using namespace h5reader::model;

namespace {
SignalDescriptor field(const char* family, SignalValueShape shape) {
    SignalDescriptor d;
    d.family = QString::fromLatin1(family);
    d.valueShape = shape;
    return d;
}
}  // namespace

class DisplayPolicyTests : public QObject {
    Q_OBJECT
private slots:
    void embeddingNotDisplayable();
    void topologyTablesNotDisplayable();
    void topologyBondLengthStaysDisplayable();
    void physicsFieldsDisplayable();
    void rollupSummariesNotDisplayable();
};

void DisplayPolicyTests::embeddingNotDisplayable() {
    // The 256-d AIMNet2 embedding is an ML feature, not a plottable curve.
    QVERIFY(!IsDashboardDisplayable(field("aimnet2", SignalValueShape::Embedding)));
}

void DisplayPolicyTests::topologyTablesNotDisplayable() {
    // The structural tables (atoms/residues/bonds/rings/ring_membership) are
    // Category-shape topology -- the molecule, not a dashboard series.
    QVERIFY(!IsDashboardDisplayable(field("topology", SignalValueShape::Category)));
}

void DisplayPolicyTests::topologyBondLengthStaysDisplayable() {
    // topology:bond_length is a real per-bond Scalar -- it MUST stay displayable
    // (this is why the rule gates on Category, not the whole topology family).
    QVERIFY(IsDashboardDisplayable(field("topology", SignalValueShape::Scalar)));
}

void DisplayPolicyTests::physicsFieldsDisplayable() {
    QVERIFY(IsDashboardDisplayable(field("biot_savart", SignalValueShape::SphericalTensor)));
    QVERIFY(IsDashboardDisplayable(field("sasa", SignalValueShape::Scalar)));
    QVERIFY(IsDashboardDisplayable(field("coulomb", SignalValueShape::EfgT2)));
    QVERIFY(IsDashboardDisplayable(field("aimnet2", SignalValueShape::Vector3)));
}

void DisplayPolicyTests::rollupSummariesNotDisplayable() {
    // Rollup-moment summaries (welford mean/var/count, *.stats, autocorrelation)
    // are whole-trajectory statistics -- in a temporal strip they draw a flat
    // line, so they are de-stripped pending a static mean +/- std readout and are
    // not dashboard signals for now. RollupMoments is exclusively this family.
    QVERIFY(!IsDashboardDisplayable(field("biot_savart", SignalValueShape::RollupMoments)));
    QVERIFY(!IsDashboardDisplayable(field("sasa", SignalValueShape::RollupMoments)));
}

QTEST_APPLESS_MAIN(DisplayPolicyTests)
#include "display_policy_tests.moc"
