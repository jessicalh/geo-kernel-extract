#include "diagnostics/StructuredLogger.h"

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEvent>
#include <QHostAddress>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLoggingCategory>
#include <QThread>
#include <QUdpSocket>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

using h5reader::diagnostics::StructuredLogger;

namespace {
constexpr const char* WorkerMessage = "structured logger worker delivery";
constexpr const char* AfterApplicationMessage = "structured logger after application destruction";
constexpr const char* LateMessage = "structured logger static destructor delivery";
std::atomic<int> previousCalls{0};
std::atomic<int> afterApplicationCalls{0};
std::atomic<int> lateCalls{0};
bool checkLateShutdown = false;

[[noreturn]] void fail(const char* reason) {
    std::fprintf(stderr, "structured logger lifetime regression: %s\n", reason);
    std::fflush(stderr);
    // Fail without running the shutdown sentinel recursively.
    std::_Exit(EXIT_FAILURE);
}

void require(bool condition, const char* reason) {
    if (!condition)
        fail(reason);
}

void previousHandler(QtMsgType, const QMessageLogContext&, const QString& message) {
    ++previousCalls;
    if (message == QLatin1String(AfterApplicationMessage))
        ++afterApplicationCalls;
    if (message == QLatin1String(LateMessage))
        ++lateCalls;
}

void verifyCategoryLookup() {
    using namespace h5reader::diagnostics;
    require(LogCategoryMaskFor("h5reader.scene") == (kCatRender | kCatFrame),
            "known category lookup changed or its storage expired");
    require(LogCategoryMaskFor("qt.core.plugin.loader") == 0,
            "unknown categories must remain ungated");
    require(LogCategoryMaskFor(nullptr) == 0, "null category must remain ungated");
}

// Constructed before main, so its destructor runs after function-local static
// state first initialized by Install/Emit. This models Qt's late plugin logs.
struct LateShutdownSentinel {
    ~LateShutdownSentinel() {
        if (!checkLateShutdown)
            return;
        verifyCategoryLookup();
        // Inspect without accidentally masking a missing restoration.
        const auto handler = qInstallMessageHandler(&previousHandler);
        require(handler == &previousHandler, "previous handler was not restored before static teardown");
        QMessageLogger(nullptr, 0, nullptr, "qt.core.plugin.loader").warning("%s", LateMessage);
        require(lateCalls.load() == 1, "late message did not reach previous handler exactly once");
        std::fprintf(stderr, "late structured logger check passed\n");
        std::fflush(stderr);
    }
} lateShutdownSentinel;

bool receiveWorkerMessage(QUdpSocket& receiver) {
    QElapsedTimer deadline;
    deadline.start();
    while (deadline.elapsed() < 3000) {
        // Worker Emit enqueues the UDP write onto the logger's application thread.
        QCoreApplication::sendPostedEvents(StructuredLogger::Instance(), QEvent::MetaCall);
        QCoreApplication::processEvents();
        if (!receiver.hasPendingDatagrams())
            receiver.waitForReadyRead(20);
        while (receiver.hasPendingDatagrams()) {
            QByteArray bytes;
            bytes.resize(static_cast<qsizetype>(receiver.pendingDatagramSize()));
            require(receiver.readDatagram(bytes.data(), bytes.size()) == bytes.size(),
                    "could not read logger UDP datagram");
            const auto record = QJsonDocument::fromJson(bytes).object();
            if (record.value(QStringLiteral("message")).toString() == QLatin1String(WorkerMessage))
                return record.value(QStringLiteral("severity")).toString() == QStringLiteral("info");
        }
    }
    return false;
}
} // namespace

int main(int argc, char** argv) {
    const bool directExit = argc == 2 && std::strcmp(argv[1], "--exit") == 0;
    require(argc == 1 || directExit, "only optional argument --exit is supported");
    qInstallMessageHandler(&previousHandler);
    verifyCategoryLookup();
    {
        QCoreApplication app(argc, argv);
        QThread::currentThread()->setObjectName(QStringLiteral("logger-test-main"));
        QUdpSocket receiver;
        require(receiver.bind(QHostAddress(QHostAddress::LocalHost), quint16(0)),
                "could not bind loopback UDP test receiver");
        const QByteArray destination = QByteArray("127.0.0.1:") + QByteArray::number(receiver.localPort());
        qputenv("H5READER_LOG_UDP", destination);
        StructuredLogger::Install();
        auto* logger = StructuredLogger::Instance();
        require(logger != nullptr, "Install did not create logger");
        StructuredLogger::Install();
        require(StructuredLogger::Instance() == logger, "Install is not idempotent");
        const int previousBeforeWorker = previousCalls.load();
        std::unique_ptr<QThread> worker(QThread::create([] {
            QThread::currentThread()->setObjectName(QStringLiteral("logger-test-worker"));
            qInfo("%s", WorkerMessage);
        }));
        worker->start();
        require(worker->wait(3000), "worker logger did not finish");
        require(receiveWorkerMessage(receiver), "worker log was not delivered through queued UDP write");
        require(previousCalls.load() == previousBeforeWorker,
                "installed logger unexpectedly bypassed to previous handler");
        worker.reset();
        checkLateShutdown = true;
        if (directExit) {
            // Mirrors QCommandLineParser help/version paths: stack QApp/logger
            // teardown is bypassed, so process-exit cleanup must unhook safely.
            std::exit(EXIT_SUCCESS);
        }
    }
    require(StructuredLogger::Instance() == nullptr,
            "logger singleton remains nonnull after application destruction");
    verifyCategoryLookup();
    qInfo("%s", AfterApplicationMessage);
    qWarning("%s", AfterApplicationMessage);
    require(afterApplicationCalls.load() == 2,
            "post-application info/warning did not reach previous handler");
    return EXIT_SUCCESS;
}
