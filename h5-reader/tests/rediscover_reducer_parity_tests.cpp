#include "rediscover/ReducerParity.h"

#include <QtTest>

using h5reader::rediscover::CompareReducerAggregateSnapshots;
using h5reader::rediscover::ReducerAggregateSnapshot;
using h5reader::rediscover::ReducerParityOptions;

class RediscoverReducerParityTests : public QObject {
    Q_OBJECT

private slots:
    void aggregateComparisonAcceptsParity();
    void aggregateComparisonReportsDrift();
};

void RediscoverReducerParityTests::aggregateComparisonAcceptsParity() {
    ReducerAggregateSnapshot a;
    a.present = true;
    a.sum_all = 1.0;
    a.sum_valid = 0.5;
    a.n_valid = 2;
    a.per_type = {0.25, 0.25};

    QStringList mismatches;
    CompareReducerAggregateSnapshots(a, a, QStringLiteral("same"),
                                     ReducerParityOptions{}, &mismatches);
    QVERIFY(mismatches.isEmpty());
}

void RediscoverReducerParityTests::aggregateComparisonReportsDrift() {
    ReducerAggregateSnapshot a;
    a.present = true;
    a.sum_all = 1.0;
    a.sum_valid = 0.5;
    a.n_valid = 2;
    a.per_type = {0.25, 0.25};

    ReducerAggregateSnapshot b = a;
    b.sum_valid = 0.75;
    b.per_type[1] = 0.5;

    QStringList mismatches;
    CompareReducerAggregateSnapshots(a, b, QStringLiteral("drift"),
                                     ReducerParityOptions{}, &mismatches);
    QVERIFY(mismatches.join(QLatin1Char('\n')).contains(QStringLiteral("sum_valid mismatch")));
    QVERIFY(mismatches.join(QLatin1Char('\n')).contains(QStringLiteral("per_type[1] mismatch")));
}

QTEST_GUILESS_MAIN(RediscoverReducerParityTests)

#include "rediscover_reducer_parity_tests.moc"
