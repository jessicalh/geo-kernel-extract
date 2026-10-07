#include <sciencefiles/BundleCache.h>
#include <sciencefiles/Download.h>

#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QSignalSpy>
#include <QSslCertificate>
#include <QSslSocket>
#include <QTemporaryDir>
#include <QtTest>

#include <chrono>
#include <memory>

using sciencefiles::Bundle;
using sciencefiles::BundleCache;
using sciencefiles::Download;
using namespace std::chrono_literals;

namespace {
constexpr qsizetype PayloadBytes = 2 * 1024 * 1024;
constexpr qint64 PrefixBytes = 512 * 1024;
const QByteArray ETag = "\"fixture-v1\"";
const QByteArray LastModified = "Mon, 02 Feb 2026 12:00:00 GMT";

QByteArray payload(int version = 1) {
    QByteArray bytes(PayloadBytes, Qt::Uninitialized);
    for (qsizetype i = 0; i < bytes.size(); ++i)
        bytes[i] = static_cast<char>((i * 17 + i / 251 + (version - 1) * 73) % 251);
    return bytes;
}

QByteArray contents(const QString &path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        QTest::qFail(qPrintable(path + ": " + file.errorString()), __FILE__, __LINE__);
        return {};
    }
    const auto bytes = file.readAll();
    if (file.error() != QFileDevice::NoError)
        QTest::qFail(qPrintable(file.errorString()), __FILE__, __LINE__);
    return bytes;
}

qint64 received(const QSignalSpy &progress) {
    return progress.isEmpty() ? 0 : progress.last().at(0).toLongLong();
}

Download::Result result(const QSignalSpy &finished) {
    return finished.last().at(0).value<Download::Result>();
}

BundleCache::Result cacheResult(const QSignalSpy &finished) {
    return finished.last().at(0).value<BundleCache::Result>();
}

QStringList entries(const QString &root) {
    return QDir(root).entryList(QDir::AllEntries | QDir::Hidden | QDir::NoDotAndDotDot, QDir::Name);
}
} // namespace

class ResumableDownloadTests final : public QObject {
    Q_OBJECT
    QTemporaryDir directory_;
    QUrl base_;
    QSslConfiguration tls_;
    QString protocol_;
    QString id_;
    QString root_;
    QString destination_;

    QUrl url(const QString &mode) const {
        return base_.resolved(QUrl("files/" + id_ + "/" + mode));
    }

    Bundle bundle(const QString &mode) const {
        return {"sample-v1", url(mode), "run.LGS", PayloadBytes + 2};
    }

    void readRecords(QJsonArray &records) {
        QNetworkAccessManager network;
        QNetworkRequest request(base_.resolved(QUrl("records/" + id_)));
        request.setSslConfiguration(tls_);
        request.setTransferTimeout(5s);
        auto *reply = network.get(request);
        QTRY_VERIFY_WITH_TIMEOUT(reply->isFinished(), 10000);
        QCOMPARE(reply->error(), QNetworkReply::NoError);
        QCOMPARE(reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt(), 200);
        const auto document = QJsonDocument::fromJson(reply->readAll());
        QVERIFY(document.isArray());
        records = document.array();
    }

    void verifyRequest(int index, qint64 offset, int status, const QByteArray &validator = ETag) {
        QJsonArray records;
        readRecords(records);
        QVERIFY(records.size() > index);
        const auto record = records[index].toObject();
        QCOMPARE(record["protocol"].toString(), protocol_);
        QCOMPARE(record["status"].toInt(), status);
        QCOMPARE(record["range"].toString(), offset > 0 ? QStringLiteral("bytes=%1-").arg(offset) : QString());
        QCOMPARE(record["if_range"].toString(), offset > 0 ? QString::fromLatin1(validator) : QString());
    }

    void verifyRequestCount(int count) {
        QJsonArray records;
        readRecords(records);
        QCOMPARE(records.size(), count);
    }

    void releasePrefix() {
        QNetworkAccessManager network;
        QNetworkRequest request(base_.resolved(QUrl("release/" + id_)));
        request.setSslConfiguration(tls_);
        request.setTransferTimeout(5s);
        auto *reply = network.post(request, QByteArray());
        QTRY_VERIFY_WITH_TIMEOUT(reply->isFinished(), 10000);
        QCOMPARE(reply->error(), QNetworkReply::NoError);
        QCOMPARE(reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt(), 204);
    }

    void verifyPartial(const QString &destination, const QByteArray &expected) {
        QVERIFY(!QFileInfo::exists(destination));
        const auto prefix = contents(destination + ".part");
        QVERIFY(!prefix.isEmpty());
        QVERIFY(prefix.size() < expected.size());
        QCOMPARE(prefix, expected.first(prefix.size()));
        QVERIFY(QJsonDocument::fromJson(contents(destination + ".json")).isObject());
    }

    void verifySaved(const QByteArray &expected) {
        QCOMPARE(contents(destination_), expected);
        QVERIFY(!QFileInfo::exists(destination_ + ".part"));
        QVERIFY(!QFileInfo::exists(destination_ + ".json"));
        QCOMPARE(entries(root_), QStringList{"download.bin"});
    }

    void stopDownload(Download &download, const QSignalSpy &finished, int expectedCount) {
        QSignalSpy stopped(&download, &Download::shutdownFinished);
        download.shutdown();
        download.shutdown();
        QTRY_COMPARE_WITH_TIMEOUT(stopped.count(), 1, 10000);
        QCOMPARE(finished.count(), expectedCount);
        QVERIFY(download.isShutdown());
        QVERIFY(!download.active());
        QVERIFY(!download.start(url("hold"), destination_, true));
    }

    void seedCancelledPrefix(Download &download, QSignalSpy &finished, const QString &mode) {
        QSignalSpy progress(&download, &Download::progress);
        QVERIFY(download.start(url(mode), destination_, true));
        QTRY_VERIFY_WITH_TIMEOUT(received(progress) >= PrefixBytes, 10000);
        download.cancel();
        download.cancel();
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QCOMPARE(result(finished), Download::Result::Cancelled);
        QVERIFY(!download.active());
        verifyPartial(destination_, payload());
    }

    void verifyPublishedBundle(BundleCache &cache, const Bundle &item) {
        QCOMPARE(cache.cachedPath(item), QDir(root_).filePath(item.key));
        QCOMPARE(contents(QDir(root_).filePath(item.key + "/run.LGS")), QByteArray("{}"));
        QCOMPARE(contents(QDir(root_).filePath(item.key + "/data.bin")), payload());
        QCOMPARE(entries(root_), QStringList{item.key});
    }

private slots:
    void initTestCase() {
        QVERIFY(directory_.isValid());
        base_ = QUrl(qEnvironmentVariable("H5READER_DOWNLOAD_URL"));
        QVERIFY2(base_.scheme() == "https" && !base_.host().isEmpty(),
                 "Run through tests/download's Go harness or the h5reader_resumable_download_tests CTest entry.");
        const auto certificates = QSslCertificate::fromData(contents(qEnvironmentVariable("H5READER_DOWNLOAD_CA")));
        QVERIFY(!certificates.isEmpty());
        tls_ = QSslConfiguration::defaultConfiguration();
        tls_.addCaCertificates(certificates);
        QVERIFY(QSslSocket::supportsSsl());
        protocol_ = qEnvironmentVariable("H5READER_DOWNLOAD_HTTP2") == "true" ? "HTTP/2.0" : "HTTP/1.1";
    }

    void init() {
        id_ = QString::fromLatin1(QTest::currentTestFunction()) + "-" + QString::fromLatin1(QTest::currentDataTag());
        root_ = directory_.filePath(id_);
        QVERIFY(QDir().mkpath(root_));
        destination_ = QDir(root_).filePath("download.bin");
    }

    void slowProgressKeepsMainThreadResponsive() {
        Download download(nullptr, tls_, 500ms);
        QSignalSpy finished(&download, &Download::finished);
        QSignalSpy progress(&download, &Download::progress);
        int servicedEvents = 0;
        bool onMainThread = true;
        connect(&download, &Download::progress, &download, [&](qint64 count, qint64 total) {
            onMainThread = onMainThread && QThread::currentThread() == thread();
            if (count > 0 && count < total) {
                QMetaObject::invokeMethod(&download, [&] {
                    if (download.active())
                        ++servicedEvents;
                }, Qt::QueuedConnection);
            }
        });
        QElapsedTimer elapsed;
        elapsed.start();
        QVERIFY(download.start(url("slow"), destination_, true));
        QVERIFY(!download.start(url("hold"), QDir(root_).filePath("duplicate.bin"), true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        QVERIFY(elapsed.elapsed() > 500);
        QVERIFY(onMainThread);
        QVERIFY(servicedEvents > 2);
        QVERIFY(progress.count() > 2);
        qint64 previous = 0;
        for (const auto &sample : progress) {
            QVERIFY(sample[0].toLongLong() >= previous);
            QCOMPARE(sample[1].toLongLong(), qint64(PayloadBytes));
            previous = sample[0].toLongLong();
        }
        QCOMPARE(previous, qint64(PayloadBytes));
        verifySaved(payload());
        stopDownload(download, finished, 1);
        verifyRequestCount(1);
        verifyRequest(0, 0, 200);
    }

    void disconnectedPrefixResumesExactly() {
        // Bound failure even if the TLS backend delays reporting the peer closure.
        Download download(nullptr, tls_, 1s);
        QSignalSpy finished(&download, &Download::finished);
        QSignalSpy progress(&download, &Download::progress);
        QVERIFY(download.start(url("disconnect"), destination_, true));
        QTRY_VERIFY_WITH_TIMEOUT(received(progress) >= PrefixBytes, 10000);
        releasePrefix();
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QCOMPARE(result(finished), Download::Result::NetworkError);
        QVERIFY(!finished.first().at(1).toString().isEmpty());
        QVERIFY(!download.active());
        verifyPartial(destination_, payload());
        const qint64 offset = QFileInfo(destination_ + ".part").size();
        QCOMPARE(offset, PrefixBytes);
        QVERIFY(download.start(url("disconnect"), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 2, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload());
        stopDownload(download, finished, 2);
        verifyRequestCount(2);
        verifyRequest(0, 0, 200);
        verifyRequest(1, offset, 206);
    }

    void cancelAndRetry_data() {
        QTest::addColumn<QString>("mode");
        QTest::addColumn<QByteArray>("validator");
        QTest::newRow("etag") << QStringLiteral("hold") << ETag;
        QTest::newRow("last-modified") << QStringLiteral("date-validator") << LastModified;
    }

    void cancelAndRetry() {
        QFETCH(QString, mode);
        QFETCH(QByteArray, validator);
        Download download(nullptr, tls_);
        QSignalSpy finished(&download, &Download::finished);
        seedCancelledPrefix(download, finished, mode);
        const qint64 offset = QFileInfo(destination_ + ".part").size();
        QVERIFY(download.start(url(mode), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 2, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload());
        stopDownload(download, finished, 2);
        verifyRequestCount(2);
        verifyRequest(1, offset, 206, validator);
        if (mode == "date-validator") {
            QJsonArray records;
            readRecords(records);
            QCOMPARE(records.size(), 2);
            QVERIFY(records[0].toObject()["etag"].toString().isEmpty());
            QCOMPARE(records[0].toObject()["last_modified"].toString(), QString::fromLatin1(LastModified));
        }
    }

    void shutdownAndNewInstanceResume() {
        {
            Download download(nullptr, tls_);
            QSignalSpy finished(&download, &Download::finished);
            QSignalSpy stopped(&download, &Download::shutdownFinished);
            QSignalSpy progress(&download, &Download::progress);
            QStringList events;
            connect(&download, &Download::finished, this, [&] { events << "finished"; });
            connect(&download, &Download::shutdownFinished, this, [&] { events << "stopped"; });
            QVERIFY(download.start(url("hold"), destination_, true));
            QTRY_VERIFY_WITH_TIMEOUT(received(progress) >= PrefixBytes, 10000);
            download.shutdown();
            download.shutdown();
            QTRY_COMPARE_WITH_TIMEOUT(stopped.count(), 1, 10000);
            QCOMPARE(finished.count(), 1);
            QCOMPARE(result(finished), Download::Result::Cancelled);
            QCOMPARE(events, (QStringList{"finished", "stopped"}));
            QVERIFY(download.isShutdown());
            QVERIFY(!download.active());
            QVERIFY(!download.start(url("hold"), destination_, true));
        }
        verifyPartial(destination_, payload());
        const qint64 offset = QFileInfo(destination_ + ".part").size();
        Download restarted(nullptr, tls_);
        QSignalSpy finished(&restarted, &Download::finished);
        QVERIFY(restarted.start(url("hold"), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload());
        stopDownload(restarted, finished, 1);
        verifyRequestCount(2);
        verifyRequest(1, offset, 206);
    }

    void fullResponseRestartsWithoutConcatenation_data() {
        QTest::addColumn<QString>("mode");
        QTest::addColumn<int>("version");
        QTest::newRow("changed-etag") << QStringLiteral("changed") << 2;
        QTest::newRow("ignored-range") << QStringLiteral("ignored-range") << 1;
    }

    void fullResponseRestartsWithoutConcatenation() {
        QFETCH(QString, mode);
        QFETCH(int, version);
        Download download(nullptr, tls_);
        QSignalSpy finished(&download, &Download::finished);
        seedCancelledPrefix(download, finished, mode);
        const qint64 offset = QFileInfo(destination_ + ".part").size();
        QVERIFY(download.start(url(mode), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 2, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload(version));
        stopDownload(download, finished, 2);
        verifyRequestCount(2);
        verifyRequest(1, offset, 200);
    }

    void rejectedResumePreservesValidatedPrefix_data() {
        QTest::addColumn<QString>("mode");
        QTest::addColumn<int>("status");
        QTest::newRow("wrong-offset") << QStringLiteral("wrong-offset") << 206;
        QTest::newRow("wrong-total") << QStringLiteral("wrong-total") << 206;
        QTest::newRow("wrong-etag") << QStringLiteral("wrong-etag") << 206;
        QTest::newRow("transient-503") << QStringLiteral("unavailable") << 503;
    }

    void rejectedResumePreservesValidatedPrefix() {
        QFETCH(QString, mode);
        QFETCH(int, status);
        Download download(nullptr, tls_);
        QSignalSpy finished(&download, &Download::finished);
        seedCancelledPrefix(download, finished, mode);
        const auto prefix = contents(destination_ + ".part");
        QVERIFY(download.start(url(mode), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 2, 10000);
        QCOMPARE(result(finished), Download::Result::HttpError);
        QVERIFY(!finished.last().at(1).toString().isEmpty());
        verifyPartial(destination_, payload());
        QCOMPARE(contents(destination_ + ".part"), prefix);
        QVERIFY(download.start(url(mode), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 3, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload());
        stopDownload(download, finished, 3);
        verifyRequestCount(3);
        verifyRequest(1, prefix.size(), status);
        verifyRequest(2, prefix.size(), 206);
    }

    void range416DiscardsPrefixAndAllowsFreshRetry() {
        Download download(nullptr, tls_);
        QSignalSpy finished(&download, &Download::finished);
        seedCancelledPrefix(download, finished, "unsatisfiable");
        const qint64 offset = QFileInfo(destination_ + ".part").size();
        QVERIFY(download.start(url("unsatisfiable"), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 2, 10000);
        QCOMPARE(result(finished), Download::Result::HttpError);
        QVERIFY(entries(root_).isEmpty());
        QVERIFY(download.start(url("unsatisfiable"), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 3, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload());
        stopDownload(download, finished, 3);
        verifyRequestCount(3);
        verifyRequest(1, offset, 416);
        verifyRequest(2, 0, 200);
    }

    void inactivityTimeoutPreservesPrefixForRetry() {
        Download download(nullptr, tls_, 500ms);
        QSignalSpy finished(&download, &Download::finished);
        QSignalSpy progress(&download, &Download::progress);
        QVERIFY(download.start(url("stall"), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QCOMPARE(result(finished), Download::Result::NetworkError);
        QVERIFY(received(progress) >= PrefixBytes);
        verifyPartial(destination_, payload());
        const qint64 offset = QFileInfo(destination_ + ".part").size();
        QVERIFY(download.start(url("stall"), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 2, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload());
        stopDownload(download, finished, 2);
        verifyRequestCount(2);
        verifyRequest(1, offset, 206);
    }

    void existingDestinationIsNeverReplaced() {
        QFile existing(destination_);
        QVERIFY(existing.open(QIODevice::WriteOnly));
        QCOMPARE(existing.write("previous complete file"), qint64(22));
        existing.close();
        Download download(nullptr, tls_);
        QSignalSpy finished(&download, &Download::finished);
        QVERIFY(download.start(url("hold"), destination_, true));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QCOMPARE(result(finished), Download::Result::FileError);
        QVERIFY(!finished.first().at(1).toString().isEmpty());
        verifySaved("previous complete file");
        stopDownload(download, finished, 1);
        verifyRequestCount(0);
    }

    void resumeOptOutDoesNotKeepSidecars() {
        Download download(nullptr, tls_);
        QSignalSpy finished(&download, &Download::finished);
        QSignalSpy progress(&download, &Download::progress);
        // Paced data gives Qt's coalesced downloadProgress signal time to fire.
        QVERIFY(download.start(url("slow"), destination_));
        QTRY_VERIFY_WITH_TIMEOUT(received(progress) > 0, 10000);
        QVERIFY(download.active());
        download.cancel();
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QCOMPARE(result(finished), Download::Result::Cancelled);
        QVERIFY(entries(root_).isEmpty());
        QVERIFY(download.start(url("slow"), destination_));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 2, 10000);
        QCOMPARE(result(finished), Download::Result::Saved);
        verifySaved(payload());
        stopDownload(download, finished, 2);
        verifyRequestCount(2);
        verifyRequest(1, 0, 200);
    }

    void bundleTransferResumesAndPublishesOnce_data() {
        QTest::addColumn<QString>("interruption");
        QTest::newRow("cancel-retry") << QStringLiteral("cancel");
        QTest::newRow("disconnect-retry") << QStringLiteral("disconnect");
        QTest::newRow("shutdown-new-instance") << QStringLiteral("shutdown");
    }

    void bundleTransferResumesAndPublishesOnce() {
        QFETCH(QString, interruption);
        const bool shutdown = interruption == "shutdown";
        const auto item = bundle("bundle-hold");
        auto cache = std::make_unique<BundleCache>(root_, nullptr, tls_);
        QSignalSpy finished(cache.get(), &BundleCache::finished);
        QSignalSpy progress(cache.get(), &BundleCache::progress);
        QSignalSpy stopped(cache.get(), &BundleCache::shutdownFinished);
        QVERIFY(cache->fetch(item));
        QTRY_VERIFY_WITH_TIMEOUT(received(progress) >= PrefixBytes, 10000);
        QVERIFY(cache->cachedPath(item).isEmpty());
        if (shutdown)
            cache->shutdown();
        else if (interruption == "disconnect")
            releasePrefix();
        else
            cache->cancel();
        // BundleCache uses the production 60-second inactivity timeout, not Download's test override.
        const int finishTimeout = interruption == "disconnect" ? 75000 : 10000;
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, finishTimeout);
        QCOMPARE(cacheResult(finished), interruption == "disconnect" ? BundleCache::Result::TransferError
                                                                    : BundleCache::Result::Cancelled);
        QVERIFY(finished.first().at(2).toString().isEmpty());
        QVERIFY(cache->cachedPath(item).isEmpty());
        QCOMPARE(entries(root_), (QStringList{".download.archive.json", ".download.archive.part"}));
        const QString partial = QDir(root_).filePath(".download.archive.part");
        const qint64 offset = QFileInfo(partial).size();
        QVERIFY(offset > 0);
        QVERIFY(QJsonDocument::fromJson(contents(QDir(root_).filePath(".download.archive.json"))).isObject());
        if (shutdown) {
            QTRY_COMPARE_WITH_TIMEOUT(stopped.count(), 1, 10000);
            QVERIFY(cache->isShutdown());
            cache = std::make_unique<BundleCache>(root_, nullptr, tls_);
        }
        QSignalSpy retried(cache.get(), &BundleCache::finished);
        QVERIFY(cache->fetch(item));
        QTRY_COMPARE_WITH_TIMEOUT(retried.count(), 1, 10000);
        QCOMPARE(cacheResult(retried), BundleCache::Result::Ready);
        QCOMPARE(retried.first().at(2).toString(), QDir(root_).filePath(item.key));
        verifyPublishedBundle(*cache, item);
        QSignalSpy finalStop(cache.get(), &BundleCache::shutdownFinished);
        cache->shutdown();
        QTRY_COMPARE_WITH_TIMEOUT(finalStop.count(), 1, 10000);
        QCOMPARE(retried.count(), 1);
        QCOMPARE(finished.count(), shutdown ? 1 : 2);
        verifyRequestCount(2);
        verifyRequest(1, offset, 206);
    }

    void bundleExtractionCancellationNeverPublishes_data() {
        QTest::addColumn<bool>("shutdown");
        QTest::newRow("cancel-extraction") << false;
        QTest::newRow("shutdown-extraction") << true;
    }

    void bundleExtractionCancellationNeverPublishes() {
        QFETCH(bool, shutdown);
        const auto item = bundle("bundle-complete");
        auto cache = std::make_unique<BundleCache>(root_, nullptr, tls_);
        QSignalSpy finished(cache.get(), &BundleCache::finished);
        QSignalSpy stopped(cache.get(), &BundleCache::shutdownFinished);
        bool interrupted = false;
        const auto connection = connect(cache.get(), &BundleCache::stateChanged, this, [&](BundleCache::State state) {
            if (state == BundleCache::State::Extracting && !interrupted) {
                interrupted = true;
                if (shutdown)
                    cache->shutdown();
                else
                    cache->cancel();
            }
        });
        QVERIFY(cache->fetch(item));
        QTRY_COMPARE_WITH_TIMEOUT(finished.count(), 1, 10000);
        QVERIFY(interrupted);
        QCOMPARE(cacheResult(finished), BundleCache::Result::Cancelled);
        QVERIFY(finished.first().at(2).toString().isEmpty());
        QVERIFY(cache->cachedPath(item).isEmpty());
        QVERIFY(QDir(root_).entryList(QDir::Dirs | QDir::Hidden | QDir::NoDotAndDotDot).isEmpty());
        disconnect(connection);
        if (shutdown) {
            QTRY_COMPARE_WITH_TIMEOUT(stopped.count(), 1, 10000);
            QVERIFY(cache->isShutdown());
            cache = std::make_unique<BundleCache>(root_, nullptr, tls_);
        }
        QSignalSpy retried(cache.get(), &BundleCache::finished);
        QVERIFY(cache->fetch(item));
        QTRY_COMPARE_WITH_TIMEOUT(retried.count(), 1, 10000);
        QCOMPARE(cacheResult(retried), BundleCache::Result::Ready);
        verifyPublishedBundle(*cache, item);
        QSignalSpy finalStop(cache.get(), &BundleCache::shutdownFinished);
        cache->shutdown();
        QTRY_COMPARE_WITH_TIMEOUT(finalStop.count(), 1, 10000);
        QCOMPARE(retried.count(), 1);
        QCOMPARE(finished.count(), shutdown ? 1 : 2);
        verifyRequest(0, 0, 200);
    }
};

QTEST_GUILESS_MAIN(ResumableDownloadTests)
#include "resumable_download_tests.moc"
