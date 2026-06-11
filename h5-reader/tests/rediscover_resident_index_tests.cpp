#include "rediscover/RunQuery.h"

#include "rediscover/Catalog.h"
#include "rediscover/ResidentIndexes.h"
#include "rediscover/RunData.h"

#include <QtTest>

using h5reader::rediscover::AcceptsFrameSelectors;
using h5reader::rediscover::ApplyAtomPrefilters;
using h5reader::rediscover::BackboneAngleSignPolicy;
using h5reader::rediscover::Body;
using h5reader::rediscover::Catalog;
using h5reader::rediscover::ChemicalCategorySelector;
using h5reader::rediscover::ChooseBackboneSignPolicy;
using h5reader::rediscover::DecodeSs8;
using h5reader::rediscover::DihedralBinSelector;
using h5reader::rediscover::DihedralKind;
using h5reader::rediscover::DihedralRangeSelector;
using h5reader::rediscover::DihedralState;
using h5reader::rediscover::FixedDihedralBin;
using h5reader::rediscover::IupacNameSelector;
using h5reader::rediscover::ResidentIndexes;
using h5reader::rediscover::RunData;
using h5reader::rediscover::SecondaryStructure3;
using h5reader::rediscover::SecondaryStructure8;
using h5reader::rediscover::SecondaryStructureSelector;
using h5reader::rediscover::SecondaryStructureState;
using h5reader::rediscover::Selector;

class RediscoverResidentIndexTests : public QObject {
    Q_OBJECT

private slots:
    void ss8DecodingCarriesExplicitUnknown();
    void fixedDihedralBinsAreDeterministic();
    void backboneSignPolicyIsChosenFromParityErrors();
    void residentIndexSelectorsUseAtomAndFramePhases();
};

void RediscoverResidentIndexTests::ss8DecodingCarriesExplicitUnknown() {
    std::array<double, 8> zero = {};
    const auto unknown = DecodeSs8(zero);
    QVERIFY(!unknown.observed);
    QCOMPARE(unknown.ss8, SecondaryStructure8::Unknown);
    QCOMPARE(unknown.ss3, SecondaryStructure3::Unknown);

    std::array<double, 8> helix = {};
    helix[0] = 1.0;
    const auto h = DecodeSs8(helix);
    QVERIFY(h.observed);
    QCOMPARE(h.ss8, SecondaryStructure8::H);
    QCOMPARE(h.ss3, SecondaryStructure3::Helix);

    std::array<double, 8> coil = {};
    coil[7] = 1.0;
    const auto c = DecodeSs8(coil);
    QVERIFY(c.observed);
    QCOMPARE(c.ss8, SecondaryStructure8::C);
    QCOMPARE(c.ss3, SecondaryStructure3::Coil);
}

void RediscoverResidentIndexTests::fixedDihedralBinsAreDeterministic() {
    constexpr double pi = 3.14159265358979323846264338327950288;
    QCOMPARE(FixedDihedralBin(-2.0 * pi / 3.0), int8_t(0));
    QCOMPARE(FixedDihedralBin(0.0), int8_t(1));
    QCOMPARE(FixedDihedralBin(2.0 * pi / 3.0), int8_t(2));
}

void RediscoverResidentIndexTests::backboneSignPolicyIsChosenFromParityErrors() {
    QCOMPARE(ChooseBackboneSignPolicy(10.0, 0.01, 12),
             BackboneAngleSignPolicy::NpyIsNegatedIupac);
    QCOMPARE(ChooseBackboneSignPolicy(0.01, 10.0, 12),
             BackboneAngleSignPolicy::NpyIsIupac);
    QCOMPARE(ChooseBackboneSignPolicy(0.5, 0.5, 12),
             BackboneAngleSignPolicy::Unresolved);
    QCOMPARE(ChooseBackboneSignPolicy(0.0, 0.0, 0),
             BackboneAngleSignPolicy::Unresolved);
}

void RediscoverResidentIndexTests::residentIndexSelectorsUseAtomAndFramePhases() {
    RunData run;
    ResidentIndexes indexes;
    indexes.iupacNames.reset(2);
    indexes.iupacNames.add(0, QStringLiteral("CA"));
    indexes.iupacNames.add(1, QStringLiteral("CB"));
    indexes.chemicalCategories.reset(2);
    indexes.chemicalCategories.add(0, QStringLiteral("backbone"));
    indexes.chemicalCategories.add(1, QStringLiteral("sidechain"));
    indexes.secondaryStructure.reset(2, 2);
    indexes.secondaryStructure.set(0, 1,
                                   SecondaryStructureState{SecondaryStructure8::H,
                                                           SecondaryStructure3::Helix,
                                                           true});
    indexes.dihedrals.reset(2, 2);
    indexes.dihedrals.set(DihedralKind::Phi, 0, 1,
                          DihedralState{0.0, FixedDihedralBin(0.0), true});

    Catalog catalog(run);
    Body body{run, indexes, catalog};

    const std::vector<Selector> selectors = {
        IupacNameSelector(QStringLiteral("CA")),
        ChemicalCategorySelector(QStringLiteral("backbone")),
        SecondaryStructureSelector(SecondaryStructure3::Helix),
        DihedralBinSelector(DihedralKind::Phi, 1),
        DihedralRangeSelector(DihedralKind::Phi, -0.1, 0.1),
    };

    const std::vector<std::size_t> scope = {0, 1};
    const std::vector<std::size_t> atoms = ApplyAtomPrefilters(body, scope, selectors);
    QCOMPARE(atoms.size(), std::size_t(1));
    QCOMPARE(atoms[0], std::size_t(0));
    QVERIFY(!AcceptsFrameSelectors(body, 0, 0, selectors));
    QVERIFY(AcceptsFrameSelectors(body, 0, 1, selectors));
}

QTEST_GUILESS_MAIN(RediscoverResidentIndexTests)

#include "rediscover_resident_index_tests.moc"
