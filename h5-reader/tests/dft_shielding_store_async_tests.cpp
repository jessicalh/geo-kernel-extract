// Standalone lifecycle tests: link DftShieldingStore.cpp, ObjectCensus.cpp,
// their required moc output, Qt Core/Test, and Eigen. This file supplies the
// loader double and thread-guard category; do not link the real loader/logger.
// Worker gates and event loops are signal-driven, with no timing assumptions.

#include "model/DftShieldingStore.h"
#include "model/QtProtein.h"
#include "io/DftShieldingLoader.h"
#include "diagnostics/ThreadGuard.h"

#include <QEventLoop>
#include <QScopeGuard>
#include <QSemaphore>
#include <QtTest>

#include <functional>
#include <stdexcept>

using h5reader::model::DftPart;
using h5reader::model::DftScalar;
using h5reader::model::DftShieldingFrame;
using h5reader::model::DftShieldingStore;
using h5reader::model::QtProtein;
using FramePtr = std::shared_ptr<const DftShieldingFrame>;

namespace {
std::function<FramePtr(const QString&, const QtProtein*)> load;

FramePtr makeFrame(double value) {
    auto frame = std::make_shared<DftShieldingFrame>();
    frame->valid = true;
    frame->atoms.resize(1);
    frame->atoms[0].total.T0 = value;
    return frame;
}

std::vector<h5reader::io::DftFrame> jobs() {
    std::vector<h5reader::io::DftFrame> frames;
    for (int index = 1; index <= 4; ++index) {
        h5reader::io::DftFrame frame;
        frame.frame_index = index;
        frame.meta_json_abspath = QString::number(index);
        frames.push_back(frame);
    }
    return frames;
}

void drain(DftShieldingStore& store) {
    QEventLoop loop;
    QObject::connect(&store, &DftShieldingStore::becameIdle, &loop, &QEventLoop::quit);
    while (store.isBusy())
        loop.exec();
}
}  // namespace

namespace h5reader::io {
FramePtr DftShieldingLoader::LoadAndValidate(const QString& path, const QtProtein* protein) {
    return load(path, protein);
}
}  // namespace h5reader::io

namespace h5reader::diagnostics {
Q_LOGGING_CATEGORY(cThreadGuard, "h5reader.thread.test")
}  // namespace h5reader::diagnostics

class DftShieldingStoreAsyncTests : public QObject {
    Q_OBJECT
private slots:
    void init();
    void deduplicatesAllFramesAndKeepsOneResident();
    void announcesAndCachesGapsAndExceptions();
    void cancellation_data();
    void cancellation();
    void cancellationAfterWorkerExitBeforeDelivery();
    void reentrantRequestsKeepPublishedFrameStable();
    void cancelFromFrameReady();
    void idleCanRestartOrDestroyStore();
    void destroyFromFrameReady_data();
    void destroyFromFrameReady();
    void destructorFallbackWaitsWithoutPublishing();
    void synchronousBatchRemainsAvailable();
};

void DftShieldingStoreAsyncTests::init() {
    load = [](const QString& path, const QtProtein*) { return makeFrame(path.toDouble()); };
}

void DftShieldingStoreAsyncTests::deduplicatesAllFramesAndKeepsOneResident() {
    QtProtein protein;
    QStringList parsed;
    bool workerInputsCorrect = true;
    load = [&](const QString& path, const QtProtein* input) {
        workerInputsCorrect = workerInputsCorrect
            && input == &protein && QThread::currentThread() != QCoreApplication::instance()->thread();
        parsed.push_back(path);
        return makeFrame(path.toDouble());
    };
    DftShieldingStore store(&protein, jobs());
    std::vector<std::size_t> ready;
    QSignalSpy idle(&store, &DftShieldingStore::becameIdle);
    connect(&store, &DftShieldingStore::frameReady, this, [&](std::size_t index) {
        QVERIFY(QThread::currentThread() == QCoreApplication::instance()->thread());
        QVERIFY(store.isBusy());
        for (auto* worker : store.findChildren<QThread*>())
            QVERIFY(worker->isFinished());
        ready.push_back(index);
        QCOMPARE(store.sample(index, 0, DftPart::Total, DftScalar::IsotropicT0).value(),
                 double(index));
        for (std::size_t other = 1; other <= 4; ++other)
            QCOMPARE(store.frame(other) != nullptr, other == index);
    });
    for (const std::size_t index : {1, 2, 1, 3, 2, 4, 4})
        store.requestFrameAsync(index);
    QVERIFY(store.isBusy());
    QVERIFY(ready.empty());
    drain(store);
    QCOMPARE(ready, (std::vector<std::size_t>{1, 2, 3, 4}));
    QCOMPARE(parsed, (QStringList{"1", "2", "3", "4"}));
    QVERIFY(workerInputsCorrect);
    QCOMPARE(idle.count(), 1);
    store.requestFrameAsync(4);
    QCOMPARE(ready.back(), std::size_t{4});
    QCOMPARE(ready.size(), std::size_t{5});
    QCOMPARE(parsed.size(), 4);
    QVERIFY(!store.isBusy());
}

void DftShieldingStoreAsyncTests::announcesAndCachesGapsAndExceptions() {
    QStringList parsed;
    load = [&](const QString& path, const QtProtein*) -> FramePtr {
        parsed.push_back(path);
        if (path == "1") return nullptr;
        if (path == "2") throw std::runtime_error("injected parser failure");
        if (path == "3") throw 42;
        return makeFrame(4);
    };
    DftShieldingStore store(nullptr, jobs());
    std::vector<std::size_t> ready;
    connect(&store, &DftShieldingStore::frameReady, this, [&](std::size_t index) {
        ready.push_back(index);
        QCOMPARE(store.frame(index) != nullptr, index == 4);
    });
    for (const std::size_t index : {1, 2, 3, 99, 4})
        store.requestFrameAsync(index);
    drain(store);
    QCOMPARE(ready, (std::vector<std::size_t>{1, 2, 3, 99, 4}));
    for (const std::size_t index : {1, 2, 3, 99}) {
        QCOMPARE(store.hasFailedFrame(index), index != 99);
        store.requestFrameAsync(index);
        QCOMPARE(ready.back(), index);
        QVERIFY(!store.isBusy());
    }
    QCOMPARE(ready.size(), std::size_t{9});
    QCOMPARE(parsed, (QStringList{"1", "2", "3", "4"}));
}

void DftShieldingStoreAsyncTests::cancellation_data() {
    QTest::addColumn<int>("outcome");
    QTest::newRow("success") << 0;
    QTest::newRow("absent") << 1;
    QTest::newRow("exception") << 2;
}

void DftShieldingStoreAsyncTests::cancellation() {
    QFETCH(int, outcome);
    QtProtein protein;
    DftShieldingStore store(&protein, jobs());
    std::vector<std::size_t> ready;
    connect(&store, &DftShieldingStore::frameReady, this,
            [&](std::size_t index) { ready.push_back(index); });
    store.requestFrame(3);
    QSignalSpy idle(&store, &DftShieldingStore::becameIdle);
    QSemaphore release;
    bool guiRanDuringLoad = false;
    load = [&](const QString&, const QtProtein* input) -> FramePtr {
        QMetaObject::invokeMethod(this, [&] {
            const auto unblock = qScopeGuard([&] { release.release(); });
            guiRanDuringLoad = true;
            QVERIFY(store.isBusy());
            QCOMPARE(idle.count(), 0);
            store.cancelPending();
            store.cancelPending();
            for (const std::size_t index : {1, 2, 3, 99})
                store.requestFrameAsync(index);
            store.requestFrame(4);  // Batch requests are also ignored while cancelling.
            QVERIFY(store.isBusy());
            QVERIFY(store.frame(3));
            QCOMPARE(ready, (std::vector<std::size_t>{3}));
        }, Qt::QueuedConnection);
        release.acquire();
        // The owner still retains the protein while cancellation drains.
        if (input->atomCount() != 0)
            throw std::runtime_error("unexpected topology");
        if (outcome == 1) return nullptr;
        if (outcome == 2) throw std::runtime_error("cancelled parser failure");
        return makeFrame(1);
    };
    store.requestFrameAsync(1);
    store.requestFrameAsync(2);
    drain(store);
    QVERIFY(guiRanDuringLoad);
    QCOMPARE(idle.count(), 1);
    QCOMPARE(ready, (std::vector<std::size_t>{3}));
    QVERIFY(store.frame(3));
    QVERIFY(!store.hasFailedFrame(1));
    QVERIFY(!store.isBusy());
    store.cancelPending();
    QCOMPARE(idle.count(), 1);
    init();
    store.requestFrameAsync(1);
    store.requestFrameAsync(2);
    drain(store);
    QCOMPARE(ready, (std::vector<std::size_t>{3, 1, 2}));
}

void DftShieldingStoreAsyncTests::cancellationAfterWorkerExitBeforeDelivery() {
    DftShieldingStore store(nullptr, jobs());
    QSignalSpy ready(&store, &DftShieldingStore::frameReady);
    QSignalSpy idle(&store, &DftShieldingStore::becameIdle);
    store.requestFrameAsync(1);
    auto* worker = store.findChild<QThread*>();
    QVERIFY(worker);
    // Test-only wait: force finished to be queued but not yet delivered.
    QVERIFY(worker->wait());
    QVERIFY(store.isBusy());
    QCOMPARE(ready.count(), 0);
    store.cancelPending();
    store.requestFrameAsync(2);
    drain(store);
    QCOMPARE(ready.count(), 0);
    QCOMPARE(idle.count(), 1);
    QVERIFY(!store.frame(1));
    store.requestFrameAsync(2);
    drain(store);
    QCOMPARE(ready.count(), 1);
    QVERIFY(store.frame(2));
}

void DftShieldingStoreAsyncTests::reentrantRequestsKeepPublishedFrameStable() {
    DftShieldingStore store(nullptr, jobs());
    std::vector<std::size_t> ready;
    connect(&store, &DftShieldingStore::frameReady, this, [&](std::size_t index) {
        ready.push_back(index);
        if (index == 1) {
            for (const std::size_t requested : {1, 3, 3, 99})
                store.requestFrameAsync(requested);
        }
    });
    connect(&store, &DftShieldingStore::frameReady, this, [&](std::size_t index) {
        QCOMPARE(store.frame(index) != nullptr, index != 99);
        QVERIFY(store.isBusy());
    });
    QSignalSpy idle(&store, &DftShieldingStore::becameIdle);
    store.requestFrameAsync(1);
    store.requestFrameAsync(2);
    drain(store);
    QCOMPARE(ready, (std::vector<std::size_t>{1, 2, 3, 99}));
    QCOMPARE(idle.count(), 1);
}

void DftShieldingStoreAsyncTests::cancelFromFrameReady() {
    DftShieldingStore store(nullptr, jobs());
    QSignalSpy ready(&store, &DftShieldingStore::frameReady);
    connect(&store, &DftShieldingStore::frameReady, this, [&](std::size_t) {
        store.cancelPending();
        store.requestFrameAsync(3);
        QVERIFY(store.isBusy());
    });
    store.requestFrameAsync(1);
    store.requestFrameAsync(2);
    drain(store);
    QCOMPARE(ready.count(), 1);
    QVERIFY(store.frame(1));
    QVERIFY(!store.isBusy());
}

void DftShieldingStoreAsyncTests::idleCanRestartOrDestroyStore() {
    auto store = std::make_unique<DftShieldingStore>(nullptr, jobs());
    int idleCount = 0;
    QEventLoop loop;
    connect(store.get(), &DftShieldingStore::becameIdle, this, [&] {
        QVERIFY(!store->isBusy());
        if (++idleCount == 1)
            store->requestFrameAsync(2);
        else {
            store.reset();
            loop.quit();
        }
    });
    store->requestFrameAsync(1);
    loop.exec();
    QCOMPARE(idleCount, 2);
    QVERIFY(!store);
}

void DftShieldingStoreAsyncTests::destructorFallbackWaitsWithoutPublishing() {
    QtProtein protein;
    bool parsed = false;
    load = [&](const QString&, const QtProtein* input) {
        parsed = input == &protein && input->atomCount() == 0;
        return makeFrame(1);
    };
    int readyCount = 0;
    {
        DftShieldingStore store(&protein, jobs());
        connect(&store, &DftShieldingStore::frameReady, this,
                [&](std::size_t) { ++readyCount; });
        store.requestFrameAsync(1);
        store.requestFrameAsync(2);
    }
    QVERIFY(parsed);
    QCoreApplication::sendPostedEvents();
    QCOMPARE(readyCount, 0);
}

void DftShieldingStoreAsyncTests::destroyFromFrameReady_data() {
    QTest::addColumn<int>("index");
    QTest::newRow("worker") << 1;
    QTest::newRow("cached-gap") << 99;
}

void DftShieldingStoreAsyncTests::destroyFromFrameReady() {
    QFETCH(int, index);
    auto store = std::make_unique<DftShieldingStore>(nullptr, jobs());
    if (index == 99)
        store->requestFrame(99);
    QEventLoop loop;
    connect(store.get(), &DftShieldingStore::frameReady, this, [&](std::size_t) {
        store.reset();
        loop.quit();
    });
    store->requestFrameAsync(static_cast<std::size_t>(index));
    if (store)
        loop.exec();
    QVERIFY(!store);
}

void DftShieldingStoreAsyncTests::synchronousBatchRemainsAvailable() {
    int parses = 0;
    bool onGuiThread = false;
    load = [&](const QString&, const QtProtein*) {
        ++parses;
        onGuiThread = QThread::currentThread() == QCoreApplication::instance()->thread();
        return makeFrame(1);
    };
    DftShieldingStore store(nullptr, jobs());
    QSignalSpy ready(&store, &DftShieldingStore::frameReady);
    QSignalSpy idle(&store, &DftShieldingStore::becameIdle);
    for (const std::size_t index : {1, 1, 99, 99})
        store.requestFrame(index);
    QCOMPARE(parses, 1);
    QCOMPARE(ready.count(), 4);
    QCOMPARE(idle.count(), 0);
    QVERIFY(onGuiThread);
    QVERIFY(!store.isBusy());
    QVERIFY(!store.frame(1));
}

QTEST_GUILESS_MAIN(DftShieldingStoreAsyncTests)
#include "dft_shielding_store_async_tests.moc"
