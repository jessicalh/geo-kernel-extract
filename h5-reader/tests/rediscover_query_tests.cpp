#include "rediscover/RunQuery.h"

#include "rediscover/Catalog.h"
#include "rediscover/LiteratureConstants.h"
#include "rediscover/ResidentIndexes.h"
#include "rediscover/RunData.h"

#include <QtTest>

using h5reader::rediscover::AcceptsFrameSelectors;
using h5reader::rediscover::ApplyAtomPrefilters;
using h5reader::rediscover::Body;
using h5reader::rediscover::Catalog;
using h5reader::rediscover::FieldRef;
using h5reader::rediscover::GatherField;
using h5reader::rediscover::GatherRank;
using h5reader::rediscover::ResidentIndexes;
using h5reader::rediscover::RunData;
using h5reader::rediscover::Selector;
using h5reader::rediscover::StaticNpyArray;

class RediscoverQueryTests : public QObject {
    Q_OBJECT

private slots:
    void fieldKindGatherUsesProducerArrayComponents();
    void catalogScalesRingShieldingFields();
    void selectorEvaluationIsTwoPhase();
};

void RediscoverQueryTests::fieldKindGatherUsesProducerArrayComponents() {
    RunData run;
    StaticNpyArray efield;
    efield.stem = QStringLiteral("apbs_E");
    efield.rows = 2;
    efield.cols = 3;
    efield.values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    run.producerArrays[static_cast<int>(h5reader::io::FieldKind::APBSE)] = efield;

    ResidentIndexes indexes;
    Catalog catalog(run);
    Body body{run, indexes, catalog};

    FieldRef ref = FieldRef::Producer(h5reader::io::FieldKind::APBSE);
    QCOMPARE(ref.rank, GatherRank::Vec3);
    const auto gathered = GatherField(body, ref, 1, 0);
    QVERIFY2(gathered.present, qPrintable(gathered.absence_reason));
    QCOMPARE(gathered.values.size(), std::size_t(3));
    QCOMPARE(gathered.values[0], 4.0);
    QCOMPARE(gathered.values[1], 5.0);
    QCOMPARE(gathered.values[2], 6.0);
}

void RediscoverQueryTests::catalogScalesRingShieldingFields() {
    RunData run;

    StaticNpyArray bs;
    bs.stem = QStringLiteral("bs_shielding");
    bs.rows = 1;
    bs.cols = 9;
    bs.values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    run.producerArrays[static_cast<int>(h5reader::io::FieldKind::BSShielding)] = bs;

    StaticNpyArray ringContrib;
    ringContrib.stem = QStringLiteral("ring_contributions");
    ringContrib.rows = 1;
    ringContrib.cols = 40;
    ringContrib.values.resize(40);
    for (std::size_t i = 0; i < ringContrib.values.size(); ++i)
        ringContrib.values[i] = static_cast<double>(i);
    run.producerArrays[static_cast<int>(h5reader::io::FieldKind::RingContributions)] = ringContrib;

    ResidentIndexes indexes;
    Catalog catalog(run);
    Body body{run, indexes, catalog};
    const double factor = h5reader::rediscover::RingCurrentPpmFactor();

    const auto tensor = catalog.valueTensor(body, h5reader::rediscover::ArrayId::KernelBs, 0, 0);
    QVERIFY(std::abs(tensor.T0 - 1.0 * factor) < 1e-12);
    QVERIFY(std::abs(tensor.T1[2] - 4.0 * factor) < 1e-12);
    QVERIFY(std::abs(tensor.T2[4] - 9.0 * factor) < 1e-12);

    const auto t2 = catalog.valueT2(body, h5reader::rediscover::ArrayId::KernelBs, 0, 0);
    QVERIFY(std::abs(t2[0] - 5.0 * factor) < 1e-12);
    QVERIFY(std::abs(t2[4] - 9.0 * factor) < 1e-12);

    const auto bsComponent = catalog.value(body, h5reader::io::FieldKind::BSShielding, 0, 0, 8);
    QVERIFY(bsComponent.has_value());
    QVERIFY(std::abs(*bsComponent - 9.0 * factor) < 1e-12);

    const auto rc8 = catalog.value(body, h5reader::io::FieldKind::RingContributions, 0, 0, 8);
    const auto rc9 = catalog.value(body, h5reader::io::FieldKind::RingContributions, 0, 0, 9);
    const auto rc17 = catalog.value(body, h5reader::io::FieldKind::RingContributions, 0, 0, 17);
    const auto rc18 = catalog.value(body, h5reader::io::FieldKind::RingContributions, 0, 0, 18);
    const auto rc27 = catalog.value(body, h5reader::io::FieldKind::RingContributions, 0, 0, 27);
    const auto rc35 = catalog.value(body, h5reader::io::FieldKind::RingContributions, 0, 0, 35);
    const auto rc36 = catalog.value(body, h5reader::io::FieldKind::RingContributions, 0, 0, 36);
    QVERIFY(rc8.has_value() && rc9.has_value() && rc17.has_value() && rc18.has_value()
            && rc27.has_value() && rc35.has_value() && rc36.has_value());
    QVERIFY(std::abs(*rc8 - 8.0) < 1e-12);
    QVERIFY(std::abs(*rc9 - 9.0 * factor) < 1e-12);
    QVERIFY(std::abs(*rc17 - 17.0 * factor) < 1e-12);
    QVERIFY(std::abs(*rc18 - 18.0) < 1e-12);
    QVERIFY(std::abs(*rc27 - 27.0 * factor) < 1e-12);
    QVERIFY(std::abs(*rc35 - 35.0 * factor) < 1e-12);
    QVERIFY(std::abs(*rc36 - 36.0) < 1e-12);
}

void RediscoverQueryTests::selectorEvaluationIsTwoPhase() {
    RunData run;
    ResidentIndexes indexes;
    Catalog catalog(run);
    Body body{run, indexes, catalog};

    const Selector evenAtom = Selector::Atom(
        QStringLiteral("even_atom"),
        [](const Body&, std::size_t atom) { return atom % 2 == 0; },
        [](const Body&, std::size_t atom, std::size_t) {
            return QStringLiteral("atom_%1").arg(atom);
        });
    const Selector frameOne = Selector::Frame(
        QStringLiteral("frame_one"),
        [](const Body&, std::size_t, std::size_t frame) { return frame == 1; },
        [](const Body&, std::size_t, std::size_t frame) {
            return QStringLiteral("frame_%1").arg(frame);
        });

    const std::vector<std::size_t> scope = {0, 1, 2, 3};
    const std::vector<Selector> selectors = {evenAtom, frameOne};
    const std::vector<std::size_t> atoms = ApplyAtomPrefilters(body, scope, selectors);
    QCOMPARE(atoms.size(), std::size_t(2));
    QCOMPARE(atoms[0], std::size_t(0));
    QCOMPARE(atoms[1], std::size_t(2));
    QVERIFY(!AcceptsFrameSelectors(body, 0, 0, selectors));
    QVERIFY(AcceptsFrameSelectors(body, 0, 1, selectors));
}

QTEST_GUILESS_MAIN(RediscoverQueryTests)

#include "rediscover_query_tests.moc"
