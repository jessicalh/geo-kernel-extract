#include <sciencefiles/Download.h>

#include <QFile>
#include <QDir>
#include <QFileInfo>
#include <QScopeGuard>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QtTest>

using sciencefiles::Download;

class LocalFileDownloadTests final : public QObject {
    Q_OBJECT

    static bool writeFile(const QString &path, const QByteArray &data) {
        QFile file(path);
        return file.open(QIODevice::WriteOnly) && file.write(data) == data.size() && file.flush();
    }

    static QByteArray readFile(const QString &path) {
        QFile file(path);
        return file.open(QIODevice::ReadOnly) ? file.readAll() : QByteArray();
    }

private slots:
    void copiesReadOnlySourceWithByteProgress() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString source = directory.filePath("source.zip");
        const QString destination = directory.filePath("copy.zip");
        QByteArray bytes(1024 * 1024 + 17, 'a');
        bytes[bytes.size() - 1] = 'z';
        QVERIFY(writeFile(source, bytes));

        const auto originalPermissions = QFile::permissions(source);
        [[maybe_unused]] const auto restorePermissions =
            qScopeGuard([&] { QFile::setPermissions(source, originalPermissions); });
        QVERIFY(QFile::setPermissions(source, QFileDevice::ReadOwner | QFileDevice::ReadUser |
                                              QFileDevice::ReadGroup | QFileDevice::ReadOther));
        QVERIFY(!(QFile::permissions(source) & QFileDevice::WriteUser));

        Download download;
        QSignalSpy progress(&download, &Download::progress);
        QSignalSpy finished(&download, &Download::finished);
        QVERIFY(download.start(QUrl::fromLocalFile(source), destination));
        QTRY_COMPARE_WITH_TIMEOUT(finished.size(), 1, 10000);
        QCOMPARE(finished.first().first().value<Download::Result>(), Download::Result::Saved);
        QCOMPARE(readFile(destination), bytes);
        QCOMPARE(readFile(source), bytes);
        QVERIFY(!(QFile::permissions(source) & QFileDevice::WriteUser));

        QVERIFY(!progress.isEmpty());
        qint64 previous = -1;
        for (const auto &sample : progress) {
            const qint64 received = sample.at(0).toLongLong();
            QCOMPARE(sample.at(1).toLongLong(), static_cast<qint64>(bytes.size()));
            QVERIFY(received > previous);
            previous = received;
        }
        QCOMPARE(previous, static_cast<qint64>(bytes.size()));
    }

    void copiesEmptySource() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString source = directory.filePath("empty.zip");
        const QString destination = directory.filePath("copy.zip");
        QVERIFY(writeFile(source, {}));

        Download download;
        QSignalSpy progress(&download, &Download::progress);
        QSignalSpy finished(&download, &Download::finished);
        QVERIFY(download.start(QUrl::fromLocalFile(source), destination));
        QTRY_COMPARE_WITH_TIMEOUT(finished.size(), 1, 10000);
        QCOMPARE(finished.first().first().value<Download::Result>(), Download::Result::Saved);
        QVERIFY(QFileInfo::exists(destination));
        QCOMPARE(QFileInfo(destination).size(), qint64(0));
        QCOMPARE(progress.size(), 1);
        QCOMPARE(progress.first().at(0).toLongLong(), qint64(0));
        QCOMPARE(progress.first().at(1).toLongLong(), qint64(0));
    }

    void cancellationPreservesSourceAndDestination() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString source = directory.filePath("source.zip");
        const QString destination = directory.filePath("copy.zip");
        const QByteArray bytes(8 * 1024 * 1024, 'x');
        QVERIFY(writeFile(source, bytes));
        QVERIFY(writeFile(destination, "previous"));

        Download download;
        QSignalSpy finished(&download, &Download::finished);
        QVERIFY(download.start(QUrl::fromLocalFile(source), destination));
        download.cancel();
        QTRY_COMPARE_WITH_TIMEOUT(finished.size(), 1, 10000);
        QCOMPARE(finished.first().first().value<Download::Result>(), Download::Result::Cancelled);
        QCOMPARE(readFile(source), bytes);
        QCOMPARE(readFile(destination), QByteArray("previous"));
    }

    void rejectsSourceAsDestination() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString source = directory.filePath("source.zip");
        QVERIFY(writeFile(source, "original"));

        Download download;
        QSignalSpy finished(&download, &Download::finished);
        QVERIFY(download.start(QUrl::fromLocalFile(source), source));
        QTRY_COMPARE_WITH_TIMEOUT(finished.size(), 1, 10000);
        QCOMPARE(finished.first().first().value<Download::Result>(), Download::Result::InvalidRequest);
        QCOMPARE(readFile(source), QByteArray("original"));
    }

    void cancelDuringCopyThenRetryLeavesOnlyTheCompleteCopy() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString source = directory.filePath("source.tar.xz");
        const QString destination = directory.filePath("copy.tar.xz");
        const QByteArray bytes(64 * 1024 * 1024, 'x');
        QVERIFY(writeFile(source, bytes));
        Download download;
        QSignalSpy finished(&download, &Download::finished);
        const auto cancel = connect(&download, &Download::progress, &download,
                                    [&download](qint64 count, qint64 total) {
            if (count > 0 && count < total)
                download.cancel();
        });
        QVERIFY(download.start(QUrl::fromLocalFile(source), destination));
        QTRY_COMPARE_WITH_TIMEOUT(finished.size(), 1, 10000);
        QCOMPARE(finished.first().first().value<Download::Result>(), Download::Result::Cancelled);
        QVERIFY(!QFileInfo::exists(destination));
        QCOMPARE(readFile(source), bytes);
        disconnect(cancel);
        QVERIFY(download.start(QUrl::fromLocalFile(source), destination));
        QTRY_COMPARE_WITH_TIMEOUT(finished.size(), 2, 10000);
        QCOMPARE(finished.last().first().value<Download::Result>(), Download::Result::Saved);
        QCOMPARE(readFile(destination), bytes);
        QCOMPARE(readFile(source), bytes);
        QCOMPARE(QDir(directory.path()).entryList(QDir::Files | QDir::Hidden),
                 (QStringList{"copy.tar.xz", "source.tar.xz"}));
    }
};

QTEST_GUILESS_MAIN(LocalFileDownloadTests)
#include "local_file_download_tests.moc"
