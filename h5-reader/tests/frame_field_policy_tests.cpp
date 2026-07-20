#include "io/FrameFieldPolicy.h"

#include <QObject>
#include <QtTest>

#include <array>

using namespace h5reader::io;

namespace {

QString stemFor(FieldKind kind) {
    const std::string_view stem = FieldSpecFor(kind).stem;
    return QString::fromLatin1(stem.data(), static_cast<qsizetype>(stem.size()));
}

bool isExpectedMopacModelInput(FieldKind kind) {
    constexpr std::array<FieldKind, 6> inputs = {
        FieldKind::MOPACChargesFullPrecision,
        FieldKind::MOPACBondValenciesFullPrecision,
        FieldKind::MOPACAtomSPopulation,
        FieldKind::MOPACAtomPPopulation,
        FieldKind::MOPACAtomDPopulation,
        FieldKind::MOPACLewisBondCount,
    };
    for (const FieldKind input : inputs) {
        if (input == kind)
            return true;
    }
    return false;
}

}  // namespace

class FrameFieldPolicyTests : public QObject {
    Q_OBJECT

private slots:
    void calculatorFieldsLoad();
    void topologyAndSpatialStorageStayOut();
    void mopacDirectLoadsExactlyModelInputs();
};

void FrameFieldPolicyTests::calculatorFieldsLoad() {
    for (const FieldSpec& spec : kFieldCatalog) {
        if (spec.kind == FieldKind::AtomsCategoryInfo || spec.group == FieldGroup::Topology
            || spec.group == FieldGroup::SpatialIndex || spec.group == FieldGroup::MOPACDirect) {
            continue;
        }
        QVERIFY2(ShouldLoadFrameField(spec.kind), qPrintable(stemFor(spec.kind)));
    }
}

void FrameFieldPolicyTests::topologyAndSpatialStorageStayOut() {
    QVERIFY(!ShouldLoadFrameField(FieldKind::Count));
    QVERIFY(!ShouldLoadFrameField(FieldKind::AtomsCategoryInfo));
    QVERIFY(!ShouldLoadFrameField(FieldKind::SpatialNeighbors));
    QVERIFY(!ShouldLoadFrameField(FieldKind::Residues));
    QVERIFY(!ShouldLoadFrameField(FieldKind::Bonds));
    QVERIFY(!ShouldLoadFrameField(FieldKind::Rings));
    QVERIFY(!ShouldLoadFrameField(FieldKind::RingMembership));
}

void FrameFieldPolicyTests::mopacDirectLoadsExactlyModelInputs() {
    int loaded = 0;
    int excluded = 0;
    for (const FieldSpec& spec : kFieldCatalog) {
        if (spec.group != FieldGroup::MOPACDirect)
            continue;
        const bool expected = isExpectedMopacModelInput(spec.kind);
        QCOMPARE(ShouldLoadFrameField(spec.kind), expected);
        QCOMPARE(IsMopacModelInput(spec.kind), expected);
        expected ? ++loaded : ++excluded;
    }
    QCOMPARE(loaded, 6);
    QCOMPARE(excluded, 33);
}

QTEST_APPLESS_MAIN(FrameFieldPolicyTests)
#include "frame_field_policy_tests.moc"
