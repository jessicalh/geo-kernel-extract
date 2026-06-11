#include "rediscover/RunQuery.h"

#include "rediscover/Catalog.h"
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
