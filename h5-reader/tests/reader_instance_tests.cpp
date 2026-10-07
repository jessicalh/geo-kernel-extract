#include "app/ReaderInstance.h"

#include <QDir>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLocalSocket>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QThread>
#include <QtTest>

using h5reader::app::ReaderInstance;

namespace {
class Secondary final : public QThread {
public:
    Secondary(const QString& root, const QString& path) : root_(root), path_(path) {}
    ~Secondary() override { wait(); }
    ReaderInstance::Result result = ReaderInstance::Result::Failed;
    QString error;

private:
    void run() override {
        ReaderInstance instance(nullptr, root_);
        result = instance.start(path_);
        error = instance.errorString();
    }
    QString root_;
    QString path_;
};
}

class ReaderInstanceTests final : public QObject {
    Q_OBJECT
private slots:
    void ownershipShutdownAndRestart() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        {
            ReaderInstance primary(nullptr, directory.path());
            QCOMPARE(primary.start({}), ReaderInstance::Result::Primary);
            QSignalSpy opened(&primary, &ReaderInstance::openRequested);
            ReaderInstance rest(nullptr, directory.path());
            QCOMPARE(rest.start("run.LGS", false), ReaderInstance::Result::Failed);
            QVERIFY(rest.errorString().contains("--rest"));
            QCOMPARE(opened.size(), 0);

            primary.shutdown();
            primary.shutdown();
            ReaderInstance duringCleanup(nullptr, directory.path());
            QCOMPARE(duringCleanup.start({}, false), ReaderInstance::Result::Failed);
            QVERIFY(QFileInfo::exists(directory.filePath("reader-instance.lock")));
        }
        QVERIFY(!QFileInfo::exists(directory.filePath("reader-instance.lock")));
        ReaderInstance restarted(nullptr, directory.path());
        QCOMPARE(restarted.start({}), ReaderInstance::Result::Primary);
    }

    void forwardsAndAcknowledges_data() {
        QTest::addColumn<QString>("path");
        QTest::newRow("activate") << QString();
        QTest::newRow("file") << QStringLiteral("protein ") + QChar(0x03b1) + QStringLiteral("/run.LGS");
    }

    void forwardsAndAcknowledges() {
        QFETCH(QString, path);
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        ReaderInstance primary(nullptr, directory.path());
        QCOMPARE(primary.start({}), ReaderInstance::Result::Primary);
        QSignalSpy opened(&primary, &ReaderInstance::openRequested);
        Secondary secondary(directory.path(), path);
        secondary.start();
        QTRY_VERIFY_WITH_TIMEOUT(secondary.isFinished(), 10000);
        QVERIFY(secondary.wait());
        QVERIFY2(secondary.result == ReaderInstance::Result::Forwarded, qPrintable(secondary.error));
        QCOMPARE(opened.size(), 1);
        QCOMPARE(opened.first().first().toString(), path.isEmpty() ? QString() : QFileInfo(path).absoluteFilePath());
    }

    void fragmentedRequest() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        ReaderInstance primary(nullptr, directory.path());
        QCOMPARE(primary.start({}), ReaderInstance::Result::Primary);
        QSignalSpy opened(&primary, &ReaderInstance::openRequested);
        QLocalSocket client;
        client.connectToServer(directory.filePath("reader-instance"));
        QTRY_COMPARE_WITH_TIMEOUT(client.state(), QLocalSocket::ConnectedState, 5000);
        auto* server = primary.findChild<QLocalServer*>();
        QVERIFY(server);
        QTRY_COMPARE_WITH_TIMEOUT(server->findChildren<QLocalSocket*>().size(), 1, 5000);
        auto* accepted = server->findChild<QLocalSocket*>();

        const QString path = directory.filePath("a quoted \"name\".LGS");
        const QByteArray request = QJsonDocument(QJsonObject{{"path", path}}).toJson(QJsonDocument::Compact) + '\n';
        const auto split = request.size() / 2;
        QCOMPARE(client.write(request.first(split)), split);
        QTRY_COMPARE_WITH_TIMEOUT(accepted->bytesAvailable(), split, 5000);
        QCOMPARE(opened.size(), 0);
        QCOMPARE(client.bytesAvailable(), qint64(0));

        QCOMPARE(client.write(request.sliced(split)), request.size() - split);
        QTRY_COMPARE_WITH_TIMEOUT(opened.size(), 1, 5000);
        QTRY_COMPARE_WITH_TIMEOUT(client.bytesAvailable(), qint64(1), 5000);
        QCOMPARE(client.readAll(), QByteArray("1"));
        QCOMPARE(opened.first().first().toString(), path);
    }

    void abandonedRequestDoesNotOpenLater() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        ReaderInstance primary(nullptr, directory.path());
        QCOMPARE(primary.start({}), ReaderInstance::Result::Primary);
        QSignalSpy opened(&primary, &ReaderInstance::openRequested);
        QLocalSocket client;
        client.connectToServer(directory.filePath("reader-instance"));
        QTRY_COMPARE_WITH_TIMEOUT(client.state(), QLocalSocket::ConnectedState, 5000);
        auto* server = primary.findChild<QLocalServer*>();
        QVERIFY(server);
        QTRY_COMPARE_WITH_TIMEOUT(server->findChildren<QLocalSocket*>().size(), 1, 5000);

        const QByteArray abandoned = "{\"path\":\"abandoned.LGS\"}\n";
        QCOMPARE(client.write(abandoned), abandoned.size());
        // Keep the primary event loop unserviced until the sender has disconnected.
        while (client.bytesToWrite() > 0)
            QVERIFY(client.waitForBytesWritten(5000));
        client.abort();

        const QString path = directory.filePath("requested.LGS");
        Secondary secondary(directory.path(), path);
        secondary.start();
        QTRY_VERIFY_WITH_TIMEOUT(secondary.isFinished(), 10000);
        QVERIFY(secondary.wait());
        QVERIFY2(secondary.result == ReaderInstance::Result::Forwarded, qPrintable(secondary.error));
        QCOMPARE(opened.size(), 1);
        QCOMPARE(opened.first().first().toString(), path);
    }
};

QTEST_GUILESS_MAIN(ReaderInstanceTests)
#include "reader_instance_tests.moc"
