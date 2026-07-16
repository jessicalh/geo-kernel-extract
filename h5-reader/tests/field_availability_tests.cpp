// field_availability_tests -- headless QtTest for the PURE classification
// policy of TrajectoryFieldAvailability (the dashboard "is this field live?"
// gate). Build() itself is coupled to the dense-H5 / per-frame-NPY I/O it
// probes, so it is integration-tested via the REST /catalog harness; this unit
// test pins the pure pieces the leg-#1 fix turns on:
//   * classifyDft       -- live ORCA from the DFT job count (the old code routed
//                          ORCA through the NPY probe and always read Absent);
//   * classifyTopology  -- the startup-loaded spine, honest about EMPTY tables
//                          (rings on a ring-less peptide -> AllMissing, not the
//                          old silent Available);
//   * isVisibleState    -- the visible/hidden contract the consumers rely on;
//   * the conservative missing-record default (Absent), which closes the
//     "unhandled source kind silently reads Available" leak.
//
// Links Qt6::Core + Qt6::Test only (mirrors h5reader_model_tests' minimal set):
// this TU never ODR-uses Build or the dense-H5 / frame-NPY probes, so it does
// not drag in the H5/conformation link web.

#include "model/TrajectoryFieldAvailability.h"

#include <QObject>
#include <QtTest>

using namespace h5reader::model;
using State = h5reader::model::TrajectoryFieldAvailabilityState;
using Extent = h5reader::model::TrajectoryFieldAvailability::TopologyExtent;

namespace {
State dftState(std::size_t jobs) {
    return TrajectoryFieldAvailability::classifyDft(jobs).state;
}
State topoState(const char* stem, const Extent& t) {
    return TrajectoryFieldAvailability::classifyTopology(QString::fromLatin1(stem), t).state;
}
}  // namespace

class FieldAvailabilityTests : public QObject {
    Q_OBJECT

private slots:
    void dftLiveWhenJobsExist();
    void dftAbsentWhenNoCampaign();
    void experimentalShieldingMlAvailabilityIsExplicit();
    void topologyLiveWhenTableNonEmpty();
    void topologyAllMissingWhenTableEmpty();
    void topologyBondLengthTracksBonds();
    void topologyUnknownStemTracksSpine();
    void visibleStateContract();
    void missingRecordDefaultsAbsent();
};

void FieldAvailabilityTests::dftLiveWhenJobsExist() {
    const auto r = TrajectoryFieldAvailability::classifyDft(751);
    QVERIFY(r.state == State::Available);
    QCOMPARE(r.finiteSamples, qsizetype(751));
    QVERIFY(TrajectoryFieldAvailability::isVisibleState(r.state));
}

void FieldAvailabilityTests::dftAbsentWhenNoCampaign() {
    // No DFT campaign (tier A / extractor-only) -> honestly Absent, not faked.
    QVERIFY(dftState(0) == State::Absent);
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(dftState(0)));
}

void FieldAvailabilityTests::experimentalShieldingMlAvailabilityIsExplicit() {
    const auto off = TrajectoryFieldAvailability::classifyExperimentalShieldingMl(false);
    QVERIFY(off.state == State::Absent);
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(off.state));

    const auto on = TrajectoryFieldAvailability::classifyExperimentalShieldingMl(true);
    QVERIFY(on.state == State::NoFramePayload);
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(on.state));
}

void FieldAvailabilityTests::topologyLiveWhenTableNonEmpty() {
    Extent t;
    t.atoms = 304; t.bonds = 310; t.residues = 20; t.rings = 8; t.ringMemberships = 46;
    QVERIFY(topoState("atoms", t) == State::Available);
    QVERIFY(topoState("residues", t) == State::Available);
    QVERIFY(topoState("bonds", t) == State::Available);
    QVERIFY(topoState("rings", t) == State::Available);
    QVERIFY(topoState("ring_membership", t) == State::Available);
    QCOMPARE(TrajectoryFieldAvailability::classifyTopology(QStringLiteral("rings"), t).finiteSamples,
             qsizetype(8));
}

void FieldAvailabilityTests::topologyAllMissingWhenTableEmpty() {
    // A ring-less peptide: rings / ring_membership genuinely empty -> honest
    // AllMissing (hidden), NOT the Available the old silent default would give.
    Extent t;
    t.atoms = 304; t.bonds = 310; t.residues = 20; t.rings = 0; t.ringMemberships = 0;
    QVERIFY(topoState("rings", t) == State::AllMissing);
    QVERIFY(topoState("ring_membership", t) == State::AllMissing);
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(topoState("rings", t)));
    QVERIFY(topoState("atoms", t) == State::Available);  // spine still live
}

void FieldAvailabilityTests::topologyBondLengthTracksBonds() {
    // topology:bond_length is derived from the bond table.
    Extent t;
    t.bonds = 310;
    QVERIFY(topoState("bond_length", t) == State::Available);
    t.bonds = 0;
    QVERIFY(topoState("bond_length", t) == State::AllMissing);
}

void FieldAvailabilityTests::topologyUnknownStemTracksSpine() {
    // A future/unknown topology stem is live iff the spine itself loaded (atoms).
    Extent t;
    t.atoms = 304;
    QVERIFY(topoState("future_table", t) == State::Available);
    t.atoms = 0;
    QVERIFY(topoState("future_table", t) == State::AllMissing);
}

void FieldAvailabilityTests::visibleStateContract() {
    QVERIFY(TrajectoryFieldAvailability::isVisibleState(State::Available));
    QVERIFY(TrajectoryFieldAvailability::isVisibleState(State::AllZeroObserved));
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(State::Absent));
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(State::AllMissing));
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(State::NoFramePayload));
    QVERIFY(!TrajectoryFieldAvailability::isVisibleState(State::AllZeroStructural));
}

void FieldAvailabilityTests::missingRecordDefaultsAbsent() {
    // An empty gate (no records) now reports Absent for any id -- the
    // conservative default that closes the "unhandled kind reads Available" leak.
    const TrajectoryFieldAvailability empty;
    QVERIFY(empty.stateForDescriptor(QStringLiteral("anything")) == State::Absent);
    QVERIFY(empty.recordForDescriptor(QStringLiteral("anything")) == nullptr);
}

QTEST_APPLESS_MAIN(FieldAvailabilityTests)
#include "field_availability_tests.moc"
