#include "../src/model/Conformation.h"
#include "../src/model/QtConformationSnapshot.h"
#include "../src/model/TransformedConformation.h"

#include <QSemaphore>
#include <QSignalSpy>
#include <QTest>
#include <QThread>
#include <atomic>
#include <stdexcept>

using namespace h5reader::model;

struct ReadControl {
    QSemaphore entered;
    QSemaphore release;
    std::atomic<int> reads{0};
    bool block = false;
    bool fail = false;
    bool values = false;
};

class TestConformation final : public Conformation {
public:
    explicit TestConformation(std::shared_ptr<ReadControl> control)
        : Conformation(nullptr), control_(std::move(control)) {}
    std::size_t frameCount() const override { return 100; }
    double timePicoseconds(std::size_t frame) const override { return double(frame); }
    Vec3 atomPosition(std::size_t, std::size_t) const override { return Vec3::Zero(); }
protected:
    SnapshotReader snapshotReader(std::size_t frame) const override {
        return [control = control_, frame] () -> std::shared_ptr<const QtConformationSnapshot> {
            ++control->reads;
            control->entered.release();
            if (control->block)
                control->release.acquire();
            if (control->fail)
                throw std::runtime_error("test read failure");
            if (control->values) {
                auto snapshot = std::make_shared<QtConformationSnapshot>(nullptr, 15 * frame, 150.0 * frame);
                auto& column = snapshot->mutableColumn(h5reader::io::FieldKind::BSTotalB);
                column.present = true;
                column.rows = 1;
                column.cols = 3;
                column.data = {double(frame), double(frame + 1), double(frame + 2)};
                return snapshot;
            }
            return {};
        };
    }
private:
    std::shared_ptr<ReadControl> control_;
};

class ConformationAsyncTests : public QObject {
    Q_OBJECT
private slots:
    void queuedFramesRetainTheirOwnDataAndOriginalIndices() {
        auto control = std::make_shared<ReadControl>();
        control->values = true;
        TestConformation source(control);
        TransformedConformation view(&source);
        std::vector<std::size_t> completed;
        connect(&view, &Conformation::snapshotReady, this, [&](std::size_t frame) {
            const auto snapshot = view.snapshot(frame);
            QVERIFY(snapshot);
            QCOMPARE(snapshot->frameIndex(), 15 * frame);
            QCOMPARE(snapshot->timePs(), 150.0 * frame);
            QCOMPARE(snapshot->column(h5reader::io::FieldKind::BSTotalB).data,
                     (std::vector<double>{double(frame), double(frame + 1), double(frame + 2)}));
            completed.push_back(frame);
        });
        std::vector<std::size_t> requested;
        for (std::size_t frame = 0; frame < 100; ++frame) {
            // Out-of-order scrubbing, including requests through both wrappers.
            const auto row = (frame * 37) % 100;
            requested.push_back(row);
            view.requestSnapshotAsync(row);
            source.requestSnapshotAsync(row);
        }
        QTRY_COMPARE_WITH_TIMEOUT(completed.size(), requested.size(), 15000);
        QCOMPARE(completed, requested);
        QCOMPARE(control->reads.load(), 100);
        QVERIFY(!view.isBusy());
    }

    void missingFramesCompleteAndRequestsAreDeduplicated() {
        auto control = std::make_shared<ReadControl>();
        TestConformation source(control);
        TransformedConformation view(&source);
        QSignalSpy ready(&view, &Conformation::snapshotReady);
        QSignalSpy idle(&source, &Conformation::becameIdle);
        source.requestSnapshotAsync(0);
        view.requestSnapshotAsync(0);
        view.requestSnapshotAsync(1);
        view.requestSnapshotAsync(2);
        QTRY_COMPARE_WITH_TIMEOUT(ready.count(), 3, 5000);
        QCOMPARE(control->reads.load(), 3);
        QCOMPARE(ready[0][0].value<std::size_t>(), std::size_t(0));
        QCOMPARE(ready[1][0].value<std::size_t>(), std::size_t(1));
        QCOMPARE(ready[2][0].value<std::size_t>(), std::size_t(2));
        QVERIFY(!view.snapshot(2));
        QVERIFY(!view.isBusy());
        QCOMPARE(idle.count(), 1);
        view.requestSnapshotAsync(2);
        QCOMPARE(ready.count(), 4);
        QCOMPARE(control->reads.load(), 3);
    }

    void cancellationKeepsEventLoopFreeAndDiscardsResult() {
        auto control = std::make_shared<ReadControl>();
        control->block = true;
        TestConformation source(control);
        QSignalSpy ready(&source, &Conformation::snapshotReady);
        QSignalSpy idle(&source, &Conformation::becameIdle);
        source.requestSnapshotAsync(0);
        source.requestSnapshotAsync(1);
        QVERIFY(control->entered.tryAcquire(1, 5000));
        source.cancelPending();
        QVERIFY(source.isBusy());
        bool guiCallback = false;
        QMetaObject::invokeMethod(&source, [&] {
            guiCallback = true;
            control->release.release();
        }, Qt::QueuedConnection);
        QTRY_COMPARE_WITH_TIMEOUT(idle.count(), 1, 5000);
        QVERIFY(guiCallback);
        QCOMPARE(ready.count(), 0);
        QCOMPARE(control->reads.load(), 1);
        control->block = false;
        source.requestSnapshotAsync(2);
        QTRY_COMPARE_WITH_TIMEOUT(ready.count(), 1, 5000);
        QCOMPARE(control->reads.load(), 2);
    }

    void readExceptionCompletesWithoutStrandingWaiters() {
        auto control = std::make_shared<ReadControl>();
        control->fail = true;
        TestConformation source(control);
        QSignalSpy ready(&source, &Conformation::snapshotReady);
        source.requestSnapshotAsync(1);
        QTRY_COMPARE_WITH_TIMEOUT(ready.count(), 1, 5000);
        QVERIFY(!source.snapshot(1));
        QVERIFY(!source.isBusy());
    }
};

QTEST_GUILESS_MAIN(ConformationAsyncTests)
#include "conformation_async_tests.moc"
