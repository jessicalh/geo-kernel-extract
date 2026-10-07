#include "app/ReaderCollection.h"
#include "app/ReaderCollectionDialog.h"
#include "app/ReaderMainWindow.h"
#include "app/QtPlaybackController.h"
#include "app/SceneVideoExporter.h"
#include "model/AtomSelection.h"
#include "model/TransformedConformation.h"

#include <QDateTime>
#include <QDialog>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QHostAddress>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QProgressBar>
#include <QSignalSpy>
#include <QSlider>
#include <QSettings>
#include <QStandardPaths>
#include <QSurfaceFormat>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QVTKOpenGLNativeWidget.h>
#include <QtTest>

using h5reader::app::ParseReaderCollection;
using h5reader::app::ReaderCollection;
using h5reader::app::ReaderCollectionDialog;
using h5reader::app::ReadReaderDocument;

class TrajectoryLibraryTests final : public QObject {
    Q_OBJECT
    QTemporaryDir directory_;
    QString catalog_;
    QString installed_;
    QString buffer_;
    QJsonObject catalogRoot_;
    ReaderCollection collection_;

    static void writeFile(const QString& path, const QByteArray& data) {
        QVERIFY(QDir().mkpath(QFileInfo(path).absolutePath()));
        QFile file(path);
        QVERIFY(file.open(QIODevice::WriteOnly));
        QCOMPARE(file.write(data), data.size());
        file.close();
        QCOMPARE(file.error(), QFileDevice::NoError);
    }

    QString bufferedRun(const QString& key = QStringLiteral("downloaded-v1")) const {
        return QDir(buffer_).filePath(key + "/run.LGS");
    }

    QString archivePath(const QString& key) const {
        return QFileInfo(catalog_).dir().filePath("packages/" + key + ".tar.xz");
    }

    void verifySourceFiles() const {
        for (const auto& entry : collection_.entries) {
            QFile archive(entry.archiveUrl.toLocalFile());
            QVERIFY(archive.open(QIODevice::ReadOnly));
            QCOMPARE(archive.readAll(), QByteArray("source-") + entry.key.toUtf8());
        }
        QFile installed(QDir(installed_).filePath("included-v1/run.LGS"));
        QVERIFY(installed.open(QIODevice::ReadOnly));
        QCOMPARE(installed.readAll(), QByteArray("{}"));
        QFile catalog(catalog_);
        QVERIFY(catalog.open(QIODevice::ReadOnly));
        QCOMPARE(QJsonDocument::fromJson(catalog.readAll()).object(), catalogRoot_);
    }

private slots:
    void initTestCase() {
        QCoreApplication::setOrganizationName("h5reader-tests");
        QCoreApplication::setApplicationName("trajectory-library");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory_.path());
        QStandardPaths::setTestModeEnabled(true);
    }

    void init() {
        QVERIFY(directory_.isValid());
        QSettings().clear();
        const QString test = QString::fromLatin1(QTest::currentTestFunction())
            + "-" + QString::fromLatin1(QTest::currentDataTag());
        catalog_ = directory_.filePath(test + "/Reader.lgs");
        installed_ = directory_.filePath(test + "/installed");
        buffer_ = directory_.filePath(test + "/buffer");
        QJsonArray entries{
            QJsonObject{{"key", "included-v1"}, {"title", "Trp-cage"},
                        {"group", "Small proteins"}, {"description", "Folding trajectory"},
                        {"organism", "Synthetic construct"}, {"pdb", "1L2Y"},
                        {"keywords", QJsonArray{"miniprotein", "folding"}}},
            QJsonObject{{"key", "downloaded-v1"}, {"title", "Ubiquitin"},
                        {"group", "MD trajectories"}, {"description", "100-frame production trajectory"},
                        {"organism", "Homo sapiens"}, {"pdb", "1UBQ"},
                        {"keywords", QJsonArray{"relaxation", "backbone dynamics"}}}};
        for (int i = 0; i < entries.size(); ++i) {
            auto entry = entries[i].toObject();
            const QString key = entry["key"].toString();
            const QByteArray archive = QByteArray("source-") + key.toUtf8();
            writeFile(archivePath(key), archive);
            entry["archive"] = "packages/" + key + ".tar.xz";
            entry["entry_point"] = "run.LGS";
            entry["frames"] = 100;
            entry["archive_bytes"] = static_cast<qint64>(archive.size());
            entry["expanded_bytes"] = 2;
            entries[i] = entry;
        }
        catalogRoot_ = {{"schema_version", 1}, {"kind", "collection"},
                        {"title", "Reader test collection"}, {"entries", entries}};
        writeFile(catalog_, QJsonDocument(catalogRoot_).toJson());
        writeFile(QDir(installed_).filePath("included-v1/run.LGS"), "{}");
        writeFile(bufferedRun(), "{}");
        QString error;
        const auto collection = ParseReaderCollection(catalogRoot_, QUrl::fromLocalFile(catalog_), &error);
        QVERIFY2(collection.has_value(), qPrintable(error));
        collection_ = *collection;
    }

    void installedExampleIsOfflineAndProtected() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QString error;
        QVERIFY2(dialog.openTrajectory("included-v1", &error), qPrintable(error));
        QCOMPARE(opened.count(), 1);
        QCOMPARE(opened.first().first().toString(), QDir(installed_).filePath("included-v1/run.LGS"));
        QCOMPARE(opened.first().at(1).toString(), QStringLiteral("included-v1"));
        QVERIFY(!dialog.clearTrajectory("included-v1", &error));
        QVERIFY(error.contains("installation"));
        QVERIFY(!dialog.isBusy());
        dialog.runFailed("included-v1", "The installed run could not be loaded.");
        QTRY_VERIFY(!dialog.isBusy());
        verifySourceFiles();
    }

    void bufferedExampleUsesTheCatalogKeyAndOpenSignal() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QString error;
        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        QCOMPARE(opened.count(), 0);
        QVERIFY(dialog.isBusy());
        QTRY_COMPARE(opened.count(), 1);
        QCOMPARE(opened.first().first().toString(), bufferedRun());
        QCOMPARE(opened.first().at(1).toString(), QStringLiteral("downloaded-v1"));
        QVERIFY(dialog.state()["error"].toString().isEmpty());
        QVERIFY(!dialog.isBusy());
        const auto entry = dialog.state()["entries"].toArray()[1].toObject();
        QCOMPARE(entry["key"].toString(), QStringLiteral("downloaded-v1"));
        QVERIFY(entry["downloaded"].toBool());
        QCOMPARE(entry["directory"].toString(), QFileInfo(bufferedRun()).absolutePath());
        QCOMPARE(dialog.findChild<QLabel*>("collectionStatus")->text(), QStringLiteral("Ready."));
        verifySourceFiles();
    }

    void emptyDialogHasNoBuiltInCatalog() {
        ReaderCollectionDialog dialog({}, nullptr, buffer_, installed_);
        const auto state = dialog.state();
        QVERIFY(state["entries"].toArray().isEmpty());
        QVERIFY(state["source"].toString().isEmpty());
        QVERIFY(state["error"].toString().isEmpty());
        QVERIFY(!state["busy"].toBool());
        auto* table = dialog.findChild<QTableWidget*>("collectionTable");
        auto* open = dialog.findChild<QPushButton*>("openCollectionRun");
        QVERIFY(table);
        QVERIFY(open);
        QCOMPARE(table->rowCount(), 0);
        QVERIFY(!open->isEnabled());
        QString error;
        QVERIFY(!dialog.openTrajectory("included-v1", &error));
        QVERIFY(!error.isEmpty());
    }

    void trajectoriesGetAccessorCreatesAnEmptyLazyDialog() {
        h5reader::app::ReaderMainWindow window;
        QVERIFY(!window.findChild<ReaderCollectionDialog*>());
        // GET /api/trajectories calls this accessor without requesting a catalog.
        auto* dialog = window.trajectoryLibrary();
        QVERIFY(dialog);
        const auto state = dialog->state();
        QVERIFY(state["entries"].isArray());
        QVERIFY(state["entries"].toArray().isEmpty());
        QVERIFY(state["source"].toString().isEmpty());
        QVERIFY(!state["busy"].toBool());
        QVERIFY(!state["visible"].toBool());
        QCOMPARE(window.findChild<ReaderCollectionDialog*>(), dialog);
        QCOMPARE(window.trajectoryLibrary(), dialog);
        QCOMPARE(window.findChildren<ReaderCollectionDialog*>().size(), 1);
        QVERIFY(!dialog->isBusy());
        QVERIFY(!dialog->isVisible());
    }

    void emptyWindowOpensWebsiteOnlyWhenChosen() {
        h5reader::app::ReaderMainWindow window;
        auto* local = window.findChild<QPushButton*>("OpenLocalRunButton");
        auto* website = window.findChild<QPushButton*>("OpenWebsiteButton");
        QVERIFY(local);
        QVERIFY(website);
        QVERIFY(!window.findChild<ReaderCollectionDialog*>());

        window.show();
        website->click();
        auto* library = window.findChild<ReaderCollectionDialog*>();
        QVERIFY(library);
        QVERIFY(library->isVisible());
        QCOMPARE(library->state()["phase"].toString(), QStringLiteral("catalog"));
        QSignalSpy stopped(library, &ReaderCollectionDialog::shutdownFinished);
        window.close();
        QTRY_COMPARE(stopped.count(), 1);
        QTRY_VERIFY(!window.isVisible());
        QVERIFY(library->isShutdown());
        QVERIFY(!library->isBusy());
        QVERIFY(library->state()["entries"].toArray().isEmpty());
        QVERIFY(library->state()["error"].toString().isEmpty());
    }

    void incompleteBufferCanBeCleared() {
        QVERIFY(QFile::remove(bufferedRun()));
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QString error;
        QVERIFY2(dialog.clearTrajectory("downloaded-v1", &error), qPrintable(error));
        QVERIFY(dialog.isClearing());
        QTRY_VERIFY(!dialog.isBusy());
        QVERIFY(!QFileInfo::exists(QFileInfo(bufferedRun()).absolutePath()));
        verifySourceFiles();
    }

    void loadedDataCannotBeCleared_data() {
        QTest::addColumn<QString>("relativeRun");
        QTest::newRow("directory") << QStringLiteral("downloaded-v1");
        QTest::newRow("manifest") << QStringLiteral("downloaded-v1/run.LGS");
        QTest::newRow("nested-directory") << QStringLiteral("downloaded-v1/data");
    }

    void loadedDataCannotBeCleared() {
        QFETCH(QString, relativeRun);
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        dialog.setCurrentRun(QDir(buffer_).filePath(relativeRun));
        QString error;
        QVERIFY(!dialog.clearTrajectory("downloaded-v1", &error));
        QVERIFY(error.contains("different trajectory"));
        QVERIFY(!dialog.isBusy());
        QVERIFY(QFileInfo::exists(bufferedRun()));
        dialog.setCurrentRun(QDir(installed_).filePath("included-v1"));
        QVERIFY2(dialog.clearTrajectory("downloaded-v1", &error), qPrintable(error));
        QTRY_VERIFY(!dialog.isBusy());
        QVERIFY(!QFileInfo::exists(bufferedRun()));
        verifySourceFiles();
    }

    void cancellingPendingOpenDoesNotLoadTheTrajectory() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QString error;
        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        QCOMPARE(opened.count(), 0);
        dialog.cancel();
        const auto* progress = dialog.findChild<QProgressBar*>();
        QVERIFY(progress);
        QCOMPARE(progress->minimum(), 0);
        QCOMPARE(progress->maximum(), 0);
        const auto* cancel = dialog.findChild<QPushButton*>("cancelCollectionTransfer");
        QVERIFY(cancel);
        QVERIFY(!cancel->isEnabled());
        QTRY_VERIFY(!dialog.isBusy());
        QCOMPARE(opened.count(), 0);
        QVERIFY(QFileInfo::exists(bufferedRun()));
        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        QTRY_COMPARE(opened.count(), 1);
        QCOMPARE(opened.first().at(1).toString(), QStringLiteral("downloaded-v1"));
        verifySourceFiles();
    }

    void readerCannotLoadAnEntryPendingClear() {
        h5reader::app::ReaderMainWindow window;
        auto* library = window.trajectoryLibrary();
        auto collection = collection_;
        collection.entries = {collection_.entries[1]};
        const QString key = "pending-clear-" + QFileInfo(directory_.path()).fileName();
        collection.entries[0].key = key;
        QString error;
        QVERIFY2(library->setCollection(collection, &error), qPrintable(error));
        const QString path = QDir(library->state()["cache_root"].toString()).filePath(key + "/run.LGS");
        writeFile(path, "{}");
        QVERIFY2(library->clearTrajectory(key, &error), qPrintable(error));
        QVERIFY(library->isClearing());
        QVERIFY(QFileInfo::exists(path));
        QVERIFY(!window.loadRunPath(path));
        QVERIFY(window.lastLoadError().contains("finish clearing"));
        QTRY_VERIFY(!library->state()["busy"].toBool());
        QVERIFY(!QFileInfo::exists(path));
        QVERIFY(!window.loadRunPath(path));
        QVERIFY(!window.lastLoadError().contains("finish clearing"));
        library->shutdown();
        QTRY_VERIFY(library->isShutdown());
    }

    void readerCannotReplaceRunDuringActiveExport() {
        const QString fixture = qEnvironmentVariable("H5READER_REST_FIXTURE");
        if (fixture.isEmpty())
            QSKIP("Set H5READER_REST_FIXTURE to exercise the active export guard.");
        h5reader::app::ReaderMainWindow window;
        QVERIFY2(window.loadRunPath(fixture), qPrintable(window.lastLoadError()));
        window.show();
        QVERIFY(QTest::qWaitForWindowExposed(&window));
        QVERIFY(window.startRestServer(QHostAddress::LocalHost, 0) != 0);
        auto* exporter = window.findChild<h5reader::app::SceneVideoExporter*>();
        QVERIFY(exporter);
        auto* conformation = window.transformedConformation();
        QSignalSpy finished(exporter, &h5reader::app::SceneVideoExporter::finished);
        h5reader::app::SceneVideoExportRequest request;
        request.outputPath = QFileInfo(catalog_).dir().filePath("active-export.mp4");
        QString error;
        QVERIFY2(exporter->start(request, &error), qPrintable(error));
        QVERIFY(exporter->isActive());

        const bool opened = window.openDocumentPath(fixture);
        const QString loadError = window.lastLoadError();
        const bool unchanged = window.transformedConformation() == conformation;
        const bool stillActive = exporter->isActive();
        exporter->requestStop(false);
        QTRY_VERIFY_WITH_TIMEOUT(!exporter->isActive(), 10000);
        QCOMPARE(finished.count(), 1);

        QVERIFY(!opened);
        QCOMPARE(loadError, QStringLiteral("Finish the active export before opening another document."));
        QVERIFY(unchanged);
        QVERIFY(stillActive);
    }

    void runReloadClosesGoToAtom() {
        const QString fixture = qEnvironmentVariable("H5READER_REST_FIXTURE");
        if (fixture.isEmpty())
            QSKIP("Set H5READER_REST_FIXTURE to exercise a real run reload.");
        h5reader::app::ReaderMainWindow window;
        QVERIFY2(window.loadRunPath(fixture), qPrintable(window.lastLoadError()));
        bool dialogFound = false;
        bool reloaded = false;
        bool rejected = false;
        QMetaObject::invokeMethod(&window, [&] {
            for (auto* dialog : window.findChildren<QDialog*>()) {
                if (dialog->windowTitle() != QStringLiteral("Go to atom"))
                    continue;
                dialogFound = true;
                QSignalSpy finished(dialog, &QDialog::finished);
                reloaded = window.loadRunPath(fixture);
                rejected = finished.count() == 1
                    && finished.first().first().toInt() == QDialog::Rejected;
                if (!rejected)
                    dialog->reject();
                break;
            }
        }, Qt::QueuedConnection);
        QVERIFY(QMetaObject::invokeMethod(&window, "onGoToAtomTriggered",
                                          Qt::DirectConnection));
        QVERIFY(dialogFound);
        QVERIFY2(reloaded, qPrintable(window.lastLoadError()));
        QVERIFY(rejected);
    }

    void inspectorLoadsDetailOnlyAfterNavigationSettles() {
        const QString fixture = qEnvironmentVariable("H5READER_REST_FIXTURE");
        if (fixture.isEmpty())
            QSKIP("Set H5READER_REST_FIXTURE to exercise real frame navigation.");
        h5reader::app::ReaderMainWindow window;
        QVERIFY2(window.loadRunPath(fixture), qPrintable(window.lastLoadError()));
        auto* selection = window.findChild<h5reader::model::AtomSelection*>();
        auto* playback = window.findChild<h5reader::app::QtPlaybackController*>();
        auto* slider = window.findChild<QSlider*>();
        auto* conformation = window.transformedConformation();
        QVERIFY(selection);
        QVERIFY(playback);
        QVERIFY(slider);
        QVERIFY(conformation);
        QVERIFY(playback->frameCount() > 4);
        selection->applyPick(0, Qt::NoModifier);
        QVERIFY(conformation->snapshot(0));
        QSignalSpy snapshots(conformation, &h5reader::model::Conformation::snapshotReady);
        QVERIFY(snapshots.isValid());

        playback->setFrame(2);
        QCOMPARE(snapshots.count(), 1);
        QVERIFY(conformation->snapshot(2));
        const auto tree = window.inspectorTreeJson();
        QVERIFY(tree[0].toObject()["value"].toString().startsWith("frame 3 /"));
        for (const auto& child : tree[0].toObject()["children"].toArray())
            QVERIFY(child.toObject()["field"].toString() != "Per-frame detail");

        snapshots.clear();
        playback->playForward();
        playback->setFrame(3);
        QCOMPARE(snapshots.count(), 0);
        QVERIFY(!conformation->snapshot(3));
        playback->pause();
        QCOMPARE(snapshots.count(), 1);
        QVERIFY(conformation->snapshot(3));

        snapshots.clear();
        slider->setSliderDown(true);
        playback->setFrame(4);
        QCOMPARE(snapshots.count(), 0);
        QVERIFY(!conformation->snapshot(4));
        slider->setSliderDown(false);
        QCOMPARE(snapshots.count(), 1);
        QVERIFY(conformation->snapshot(4));
    }

    void shutdownDoesNotOpenAQueuedCacheHit() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QSignalSpy stopped(&dialog, &ReaderCollectionDialog::shutdownFinished);
        QString error;
        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        QCOMPARE(opened.count(), 0);
        dialog.shutdown();
        QTRY_COMPARE(stopped.count(), 1);
        QVERIFY(dialog.isShutdown());
        QVERIFY(!dialog.isBusy());
        QCOMPARE(opened.count(), 0);
        QVERIFY(QFileInfo::exists(bufferedRun()));
        verifySourceFiles();
    }

    void searchMatchesMetadata_data() {
        QTest::addColumn<QString>("query");
        QTest::addColumn<int>("row");
        QTest::newRow("title") << QStringLiteral("TRP-CAGE") << 0;
        QTest::newRow("key") << QStringLiteral("downloaded-v1") << 1;
        QTest::newRow("group") << QStringLiteral("small proteins") << 0;
        QTest::newRow("description") << QStringLiteral("production") << 1;
        QTest::newRow("keyword") << QStringLiteral("RELAXATION") << 1;
        QTest::newRow("pdb") << QStringLiteral("1ubq") << 1;
        QTest::newRow("organism") << QStringLiteral("homo sapiens") << 1;
        QTest::newRow("multiword-across-fields") << QStringLiteral("  1UBQ\tHOMO  dynamics \nrelaxation ") << 1;
        QTest::newRow("all-words-required") << QStringLiteral("1UBQ folding") << -1;
        QTest::newRow("no-match") << QStringLiteral("not in this catalog") << -1;
    }

    void searchMatchesMetadata() {
        QFETCH(QString, query);
        QFETCH(int, row);
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        auto* search = dialog.findChild<QLineEdit*>("collectionSearch");
        auto* table = dialog.findChild<QTableWidget*>("collectionTable");
        auto* open = dialog.findChild<QPushButton*>("openCollectionRun");
        QVERIFY(search);
        QVERIFY(table);
        QVERIFY(open);
        QCOMPARE(search->placeholderText(), QStringLiteral("Search BMRB, protein, organism, or PDB"));
        QCOMPARE(table->rowCount(), 2);
        QCOMPARE(table->columnCount(), 6);
        QCOMPARE(table->item(1, 2)->text(), QStringLiteral("Homo sapiens"));
        QCOMPARE(table->item(1, 3)->text(), QStringLiteral("1UBQ"));
        search->setText(query);
        for (int i = 0; i < table->rowCount(); ++i)
            QCOMPARE(table->isRowHidden(i), i != row);
        QCOMPARE(open->isEnabled(), row >= 0);
        if (row >= 0)
            QCOMPARE(table->currentRow(), row);
        else {
            QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
            open->click();
            QCOMPARE(opened.count(), 0);
        }
        search->clear();
        QVERIFY(!table->isRowHidden(0));
        QVERIFY(!table->isRowHidden(1));
        QVERIFY(open->isEnabled());
    }

    void sameCatalogParsesFromDiskAndHttps() {
        QString error;
        const auto document = ReadReaderDocument(catalog_, &error);
        QVERIFY2(document.has_value(), qPrintable(error));
        QVERIFY(document->collection);
        const QUrl source(QStringLiteral("https://example.test/reader/Reader.lgs"));
        const auto remote = ParseReaderCollection(catalogRoot_, source, &error);
        QVERIFY2(remote.has_value(), qPrintable(error));
        QCOMPARE(document->collection->source, QUrl::fromLocalFile(catalog_));
        QCOMPARE(remote->source, source);
        QCOMPARE(remote->title, document->collection->title);
        QCOMPARE(remote->entries.size(), document->collection->entries.size());
        for (int i = 0; i < remote->entries.size(); ++i) {
            const auto& localEntry = document->collection->entries[i];
            const auto& remoteEntry = remote->entries[i];
            QCOMPARE(localEntry.key, remoteEntry.key);
            QCOMPARE(localEntry.title, remoteEntry.title);
            QCOMPARE(localEntry.group, remoteEntry.group);
            QCOMPARE(localEntry.description, remoteEntry.description);
            QCOMPARE(localEntry.organism, remoteEntry.organism);
            QCOMPARE(localEntry.pdb, remoteEntry.pdb);
            QCOMPARE(localEntry.keywords, remoteEntry.keywords);
            QCOMPARE(localEntry.entryPoint, remoteEntry.entryPoint);
            QCOMPARE(localEntry.frames, remoteEntry.frames);
            QCOMPARE(localEntry.archiveBytes, remoteEntry.archiveBytes);
            QCOMPARE(localEntry.expandedBytes, remoteEntry.expandedBytes);
            QCOMPARE(localEntry.archiveUrl, QUrl::fromLocalFile(archivePath(localEntry.key)));
            QCOMPARE(remoteEntry.archiveUrl,
                     QUrl("https://example.test/reader/packages/" + remoteEntry.key + ".tar.xz"));
        }
        QCOMPARE(remote->entries[1].organism, QStringLiteral("Homo sapiens"));
        QCOMPARE(remote->entries[1].pdb, QStringLiteral("1UBQ"));
        QCOMPARE(remote->entries[1].keywords, (QStringList{"relaxation", "backbone dynamics"}));

        ReaderCollectionDialog dialog({}, nullptr, buffer_, installed_);
        QVERIFY2(dialog.loadCatalog(QUrl::fromLocalFile(catalog_), &error), qPrintable(error));
        const auto localEntries = dialog.state()["entries"].toArray();
        QVERIFY2(dialog.setCollection(*remote, &error), qPrintable(error));
        QCOMPARE(dialog.state()["source"].toString(), source.toString());
        const auto remoteEntries = dialog.state()["entries"].toArray();
        QCOMPARE(remoteEntries.size(), localEntries.size());
        for (int i = 0; i < localEntries.size(); ++i) {
            auto localEntry = localEntries[i].toObject();
            auto remoteEntry = remoteEntries[i].toObject();
            localEntry.remove("url");
            remoteEntry.remove("url");
            QCOMPARE(localEntry, remoteEntry);
        }
        verifySourceFiles();
    }

    void unsafeArchivePathsAreRejected_data() {
        QTest::addColumn<QString>("archive");
        QTest::newRow("parent") << QStringLiteral("../outside.tar.xz");
        QTest::newRow("nested-parent") << QStringLiteral("packages/../../outside.tar.xz");
        QTest::newRow("backslash-parent") << QStringLiteral("..\\outside.tar.xz");
        QTest::newRow("absolute-posix") << QStringLiteral("/outside.tar.xz");
        QTest::newRow("absolute-windows") << QStringLiteral("C:/outside.tar.xz");
        QTest::newRow("drive-relative") << QStringLiteral("C:outside.tar.xz");
        QTest::newRow("unc") << QStringLiteral("\\\\server\\share\\outside.tar.xz");
        QTest::newRow("network-path") << QStringLiteral("//example.test/outside.tar.xz");
        QTest::newRow("file-url") << QStringLiteral("file:///C:/outside.tar.xz");
        QTest::newRow("https-url") << QStringLiteral("https://example.test/outside.tar.xz");
        QTest::newRow("http-url") << QStringLiteral("http://example.test/outside.tar.xz");
        QTest::newRow("query") << QStringLiteral("packages/run.tar.xz?other.tar.xz");
        QTest::newRow("fragment") << QStringLiteral("packages/run.tar.xz#other.tar.xz");
    }

    void unsafeArchivePathsAreRejected() {
        QFETCH(QString, archive);
        auto root = catalogRoot_;
        auto entries = root["entries"].toArray();
        auto entry = entries[0].toObject();
        entry["archive"] = archive;
        entries[0] = entry;
        root["entries"] = entries;
        for (const auto& source : {QUrl::fromLocalFile(catalog_),
                                   QUrl(QStringLiteral("https://example.test/reader/Reader.lgs"))}) {
            QString error;
            QVERIFY(!ParseReaderCollection(root, source, &error));
            QVERIFY2(!error.isEmpty(), qPrintable(source.toString()));
        }
        verifySourceFiles();
    }

    void unsafeCatalogKeysAreRejected_data() {
        QTest::addColumn<QString>("key");
        QTest::newRow("empty") << QString();
        QTest::newRow("current") << QStringLiteral(".");
        QTest::newRow("parent") << QStringLiteral("..");
        QTest::newRow("traversal") << QStringLiteral("../other");
        QTest::newRow("backslash") << QStringLiteral("..\\other");
        QTest::newRow("drive") << QStringLiteral("C:other");
        QTest::newRow("duplicate") << QStringLiteral("downloaded-v1");
        QTest::newRow("duplicate-case") << QStringLiteral("DOWNLOADED-V1");
    }

    void unsafeCatalogKeysAreRejected() {
        QFETCH(QString, key);
        auto root = catalogRoot_;
        auto entries = root["entries"].toArray();
        auto entry = entries[0].toObject();
        entry["key"] = key;
        entries[0] = entry;
        root["entries"] = entries;
        QString error;
        QVERIFY(!ParseReaderCollection(root, QUrl::fromLocalFile(catalog_), &error));
        QVERIFY(!error.isEmpty());
    }

    void runDocumentIsNotACollection() {
        const QString run = QFileInfo(catalog_).dir().filePath("source/run.LGS");
        writeFile(run, R"({"schema_version":1,"kind":"trajectory"})");
        QString error;
        const auto document = ReadReaderDocument(run, &error);
        QVERIFY2(document.has_value(), qPrintable(error));
        QVERIFY(!document->collection);
        ReaderCollectionDialog dialog({}, nullptr, buffer_, installed_);
        QVERIFY(!dialog.loadCatalog(QUrl::fromLocalFile(run), &error));
        QVERIFY(error.contains("not a collection"));
        QVERIFY(dialog.state()["entries"].toArray().isEmpty());
        QVERIFY(QFileInfo::exists(run));
    }

    void collectionOpensBufferedRunsAndFiltersThePicker() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        auto* search = dialog.findChild<QLineEdit*>("collectionSearch");
        auto* table = dialog.findChild<QTableWidget*>("collectionTable");
        auto* open = dialog.findChild<QPushButton*>("openCollectionRun");
        QVERIFY(search);
        QVERIFY(table);
        QVERIFY(open);
        QSignalSpy selected(&dialog, &ReaderCollectionDialog::openRequested);
        QCOMPARE(table->rowCount(), 2);
        search->setText("Homo sapiens");
        QVERIFY(table->isRowHidden(0));
        QVERIFY(!table->isRowHidden(1));
        QVERIFY(open->isEnabled());
        open->click();
        QTRY_COMPARE(selected.size(), 1);
        QCOMPARE(selected.first().first().toString(), bufferedRun());
        QCOMPARE(selected.first().at(1).toString(), QStringLiteral("downloaded-v1"));
        dialog.runOpened("downloaded-v1");
        dialog.setCurrentRun(QFileInfo(bufferedRun()).absolutePath());
        search->setText("no matches");
        QVERIFY(!open->isEnabled());
        search->setText("Trp-cage");
        QVERIFY(open->isEnabled());
        open->click();
        QTRY_COMPARE(selected.size(), 2);
        QCOMPARE(selected.last().first().toString(), QDir(installed_).filePath("included-v1/run.LGS"));
        QCOMPARE(selected.last().at(1).toString(), QStringLiteral("included-v1"));
        QVERIFY(QFileInfo::exists(bufferedRun()));
        dialog.setCurrentRun(QDir(installed_).filePath("included-v1"));
        dialog.runOpened("included-v1");
        QTRY_VERIFY(!dialog.isBusy());
        QVERIFY(!QFileInfo::exists(bufferedRun()));
        verifySourceFiles();
    }

    void collectionReportsMissingPackage() {
        const QString root = directory_.filePath("missing-collection-entry");
        const QString index = QDir(root).filePath("Reader Library.lgs");
        const QString key = "missing-package-" + QFileInfo(directory_.path()).fileName();
        const QJsonArray entries{
            QJsonObject{{"key", key}, {"title", "BMRB 68"},
                        {"group", "MD trajectories"}, {"frames", 100},
                        {"archive", "packages/bmr68-reader-v1.tar.xz"},
                        {"entry_point", "run.LGS"}, {"archive_bytes", 12}, {"expanded_bytes", 2}}};
        writeFile(index, QJsonDocument(QJsonObject{{"schema_version", 1}, {"kind", "collection"},
                                                   {"title", "Reader datasets"}, {"entries", entries}}).toJson());

        h5reader::app::ReaderMainWindow window;
        window.show();
        QVERIFY2(window.openDocumentPath(index), qPrintable(window.lastLoadError()));
        auto* dialog = window.findChild<h5reader::app::ReaderCollectionDialog*>();
        QVERIFY(dialog);
        QVERIFY(dialog->isVisible());
        auto* open = dialog->findChild<QPushButton*>("openCollectionRun");
        QVERIFY(open);
        open->click();
        QVERIFY(dialog->isVisible());
        auto* status = dialog->findChild<QLabel*>("collectionStatus");
        QVERIFY(status);
        QVERIFY(status->text().contains("missing or incomplete"));
    }

    void collectionCancelAndFailedLoadNeverTouchTheArchive() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy selected(&dialog, &ReaderCollectionDialog::openRequested);
        QString error;
        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        dialog.cancelPendingOpen();
        QTRY_VERIFY(!dialog.isBusy());
        QCOMPARE(selected.size(), 0);
        QVERIFY(QFileInfo::exists(bufferedRun()));
        verifySourceFiles();

        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        QTRY_COMPARE(selected.size(), 1);
        dialog.runFailed("downloaded-v1", "The run manifest is invalid.");
        QTRY_VERIFY(!dialog.isBusy());
        QCOMPARE(dialog.state()["error"].toString(), QStringLiteral("The run manifest is invalid."));
        QVERIFY(!QFileInfo::exists(bufferedRun()));
        verifySourceFiles();
        const auto* status = dialog.findChild<QLabel*>("collectionStatus");
        QVERIFY(status);
        QCOMPARE(status->text(), QStringLiteral("The run manifest is invalid."));
    }

    void failedPackageCanBeRetriedWithoutOpeningOrDuplicatingIt() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QString error;
        QVERIFY(dialog.clearTrajectory("downloaded-v1", &error));
        QTRY_VERIFY(!dialog.isBusy());
        auto* table = dialog.findChild<QTableWidget*>("collectionTable");
        auto* open = dialog.findChild<QPushButton*>("openCollectionRun");
        QVERIFY(table && open);
        table->selectRow(1);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        for (int attempt = 0; attempt < 2; ++attempt) {
            open->click();
            QVERIFY(dialog.isBusy());
            QVERIFY(!open->isEnabled());
            QTRY_VERIFY(!dialog.isBusy());
            QVERIFY(!dialog.state()["error"].toString().isEmpty());
            QCOMPARE(open->text(), QStringLiteral("Retry"));
            QVERIFY(open->isEnabled());
            QCOMPARE(opened.count(), 0);
            QVERIFY(QDir(buffer_).entryList(QDir::AllEntries | QDir::Hidden | QDir::NoDotAndDotDot).isEmpty());
            verifySourceFiles();
        }
    }

    void failedCatalogOffersRetry() {
        ReaderCollectionDialog dialog({}, nullptr, buffer_, installed_);
        auto* open = dialog.findChild<QPushButton*>("openCollectionRun");
        QVERIFY(open);
        QVERIFY(dialog.loadCatalog(QUrl("https://127.0.0.1:1/Reader.lgs")));
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 10000);
        QVERIFY(!dialog.state()["error"].toString().isEmpty());
        QCOMPARE(open->text(), QStringLiteral("Retry collection"));
        QVERIFY(open->isEnabled());
        open->click();
        QVERIFY(dialog.isBusy());
        QVERIFY(!open->isEnabled());
        dialog.cancel();
        QTRY_VERIFY(!dialog.isBusy());
        QVERIFY(!open->isEnabled());
        QCOMPARE(dialog.state()["entries"].toArray().size(), 0);
        verifySourceFiles();
    }

    void failedCatalogLeavesExistingEntriesUsable() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        auto* open = dialog.findChild<QPushButton*>("openCollectionRun");
        QVERIFY(open);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QVERIFY(dialog.loadCatalog(QUrl("https://127.0.0.1:1/Reader.lgs")));
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 10000);
        QVERIFY(!dialog.state()["error"].toString().isEmpty());
        QCOMPARE(open->text(), QStringLiteral("Open"));
        QVERIFY(open->isEnabled());
        open->click();
        QCOMPARE(opened.count(), 1);
        QVERIFY(!dialog.isBusy());
        verifySourceFiles();
    }

    void collectionReplacementRetainsPreviousBuffer_data() {
        QTest::addColumn<bool>("https");
        QTest::newRow("local-to-local") << false;
        QTest::newRow("local-to-https") << true;
    }

    void collectionReplacementRetainsPreviousBuffer() {
        QFETCH(bool, https);
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QString error;
        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        QTRY_COMPARE(opened.count(), 1);
        dialog.runOpened("downloaded-v1");
        dialog.setCurrentRun(QFileInfo(bufferedRun()).absolutePath());

        auto root = catalogRoot_;
        auto entry = root["entries"].toArray()[1].toObject();
        entry["key"] = "replacement-v2";
        entry["archive"] = "packages/replacement-v2.tar.xz";
        root["title"] = "Replacement collection";
        root["entries"] = QJsonArray{entry};
        const QString replacementCatalog = QFileInfo(catalog_).dir().filePath("replacement/Reader.lgs");
        const QUrl source = https ? QUrl(QStringLiteral("https://example.test/reader/Reader.lgs"))
                                  : QUrl::fromLocalFile(replacementCatalog);
        if (https) {
            const auto replacement = ParseReaderCollection(root, source, &error);
            QVERIFY2(replacement.has_value(), qPrintable(error));
            QVERIFY2(dialog.setCollection(*replacement, &error), qPrintable(error));
        } else {
            writeFile(replacementCatalog, QJsonDocument(root).toJson());
            QVERIFY2(dialog.loadCatalog(source, &error), qPrintable(error));
        }
        QCOMPARE(dialog.state()["source"].toString(), source.toString());
        QCOMPARE(dialog.state()["entries"].toArray().size(), 1);
        QCOMPARE(dialog.state()["entries"].toArray()[0].toObject()["key"].toString(),
                 QStringLiteral("replacement-v2"));
        QVERIFY(QFileInfo::exists(bufferedRun()));
        QVERIFY(!dialog.isBusy());

        // Ready is only a load request. A rejected replacement must keep the current run.
        const QString replacementRun = bufferedRun("replacement-v2");
        writeFile(replacementRun, "{}");
        QVERIFY2(dialog.openTrajectory("replacement-v2", &error), qPrintable(error));
        QTRY_COMPARE(opened.count(), 2);
        QCOMPARE(opened.last().first().toString(), replacementRun);
        QCOMPARE(opened.last().at(1).toString(), QStringLiteral("replacement-v2"));
        QVERIFY(QFileInfo::exists(bufferedRun()));
        dialog.runFailed("replacement-v2", "Replacement run could not be loaded.");
        QTRY_VERIFY(!dialog.isBusy());
        QVERIFY(!QFileInfo::exists(replacementRun));
        QVERIFY(QFileInfo::exists(bufferedRun()));
        verifySourceFiles();

        writeFile(replacementRun, "{}");
        QVERIFY2(dialog.openTrajectory("replacement-v2", &error), qPrintable(error));
        QTRY_COMPARE(opened.count(), 3);
        QVERIFY(QFileInfo::exists(bufferedRun()));
        dialog.setCurrentRun(QFileInfo(replacementRun).absolutePath());
        dialog.runOpened("replacement-v2");
        QVERIFY(dialog.isClearing());
        QTRY_VERIFY(!dialog.isBusy());
        QVERIFY(!QFileInfo::exists(bufferedRun()));
        QVERIFY(QFileInfo::exists(replacementRun));
        QVERIFY(!dialog.clearTrajectory("replacement-v2", &error));
        QVERIFY(error.contains("different trajectory"));
        verifySourceFiles();
    }

    void mainWindowReusesDialogAcrossCollections() {
        h5reader::app::ReaderMainWindow window;
        QVERIFY2(window.openDocumentPath(catalog_), qPrintable(window.lastLoadError()));
        auto* dialog = window.trajectoryLibrary();
        QCOMPARE(dialog->state()["source"].toString(), QUrl::fromLocalFile(catalog_).toString());
        auto replacement = catalogRoot_;
        replacement["title"] = "Another collection";
        replacement["entries"] = QJsonArray{catalogRoot_["entries"].toArray()[1]};
        const QString path = QFileInfo(catalog_).dir().filePath("another/Reader.lgs");
        writeFile(path, QJsonDocument(replacement).toJson());
        QVERIFY2(window.openDocumentPath(path), qPrintable(window.lastLoadError()));
        QCOMPARE(window.trajectoryLibrary(), dialog);
        QCOMPARE(window.findChildren<ReaderCollectionDialog*>().size(), 1);
        QCOMPARE(dialog->windowTitle(), QStringLiteral("Another collection"));
        QCOMPARE(dialog->state()["source"].toString(), QUrl::fromLocalFile(path).toString());
        QCOMPARE(dialog->state()["entries"].toArray().size(), 1);
        QVERIFY(dialog->isVisible());
    }

    void busyDialogRejectsCollectionReplacement() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QString error;
        const auto originalEntries = dialog.state()["entries"].toArray();
        QVERIFY2(dialog.openTrajectory("downloaded-v1", &error), qPrintable(error));
        QVERIFY(!dialog.setCollection({}, &error));
        QVERIFY(!error.isEmpty());
        QVERIFY(!dialog.loadCatalog(QUrl::fromLocalFile(catalog_), &error));
        QCOMPARE(dialog.state()["entries"].toArray(), originalEntries);
        dialog.cancel();
        QTRY_VERIFY(!dialog.isBusy());
        QCOMPARE(opened.count(), 0);
        QVERIFY(QFileInfo::exists(bufferedRun()));
    }

    void closingDuringHttpsCatalogRequest_data() {
        QTest::addColumn<bool>("shutdown");
        QTest::newRow("close-dialog") << false;
        QTest::newRow("shutdown") << true;
    }

    void closingDuringHttpsCatalogRequest() {
        QFETCH(bool, shutdown);
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy opened(&dialog, &ReaderCollectionDialog::openRequested);
        QSignalSpy stopped(&dialog, &ReaderCollectionDialog::shutdownFinished);
        const auto originalEntries = dialog.state()["entries"].toArray();
        dialog.show();
        QString error;
        QVERIFY2(dialog.loadCatalog(QUrl(QStringLiteral("https://semantic.construction/files/Reader.lgs")),
                                    &error), qPrintable(error));
        QCOMPARE(dialog.state()["phase"].toString(), QStringLiteral("catalog"));
        QVERIFY(dialog.isBusy());

        // Cancel before pumping events: this exercises QNetworkReply lifetime, not server availability.
        if (shutdown)
            dialog.shutdown();
        else
            dialog.close();
        QTRY_VERIFY(!dialog.isBusy());
        if (shutdown) {
            QTRY_COMPARE(stopped.count(), 1);
            QVERIFY(dialog.isShutdown());
        } else {
            QVERIFY(!dialog.isVisible());
            QVERIFY(!dialog.isShutdown());
            QCOMPARE(stopped.count(), 0);
        }
        QCOMPARE(opened.count(), 0);
        QCOMPARE(dialog.state()["entries"].toArray(), originalEntries);
        QCOMPARE(dialog.state()["source"].toString(), collection_.source.toString());
        QVERIFY(dialog.state()["error"].toString().isEmpty());
        const auto* status = dialog.findChild<QLabel*>("collectionStatus");
        QVERIFY(status);
        QCOMPARE(status->text(), QStringLiteral("Cancelled."));
        QVERIFY(QFileInfo::exists(bufferedRun()));
        verifySourceFiles();
        if (!shutdown) {
            QVERIFY2(dialog.loadCatalog(QUrl::fromLocalFile(catalog_), &error), qPrintable(error));
            dialog.shutdown();
            QTRY_COMPARE(stopped.count(), 1);
        }
    }

    void realPackageLoadsThroughTheReader() {
        const QString index = qEnvironmentVariable("H5READER_COLLECTION_FIXTURE");
        if (index.isEmpty())
            QSKIP("Set H5READER_COLLECTION_FIXTURE to check a published package.");
        QString error;
        const auto document = h5reader::app::ReadReaderDocument(index, &error);
        QVERIFY2(document.has_value(), qPrintable(error));
        QVERIFY(document->collection);
        QCOMPARE(document->collection->entries.size(), 1);
        QVERIFY(document->collection->entries[0].archiveUrl.isLocalFile());
        const QFileInfo source(document->collection->entries[0].archiveUrl.toLocalFile());
        QVERIFY(source.isFile());
        const qint64 sourceBytes = source.size();
        const QDateTime sourceModified = source.lastModified();
        QTemporaryDir buffer;
        QVERIFY(buffer.isValid());
        h5reader::app::ReaderMainWindow window;
        ReaderCollectionDialog dialog(*document->collection, nullptr, buffer.path(),
                                      directory_.filePath("empty-installation"));
        bool finished = false;
        bool loaded = false;
        QObject::connect(&dialog, &h5reader::app::ReaderCollectionDialog::openRequested,
                         &window, [&](const QString& lgs, const QString& key) {
                             loaded = window.loadRunPath(lgs, false);
                             if (loaded)
                                 dialog.runOpened(key);
                             else
                                 dialog.runFailed(key, window.lastLoadError());
                             finished = true;
                         });
        auto* open = dialog.findChild<QPushButton*>("openCollectionRun");
        QVERIFY(open);
        open->click();
        QVERIFY(dialog.isBusy());
        QTRY_VERIFY_WITH_TIMEOUT(finished || !dialog.isBusy(), 600000);
        QVERIFY2(finished, qPrintable(dialog.findChild<QLabel*>("collectionStatus")->text()));
        QVERIFY2(loaded, qPrintable(window.lastLoadError()));
        QCOMPARE(dialog.findChild<QLabel*>("collectionStatus")->text(), QStringLiteral("Ready."));
        QVERIFY(window.uiStateJson().value("loaded").toBool());
        QCOMPARE(QFileInfo(source.absoluteFilePath()).size(), sourceBytes);
        QCOMPARE(QFileInfo(source.absoluteFilePath()).lastModified(), sourceModified);
    }

    void shutdownRejectsNewWork() {
        ReaderCollectionDialog dialog(collection_, nullptr, buffer_, installed_);
        QSignalSpy stopped(&dialog, &ReaderCollectionDialog::shutdownFinished);
        dialog.shutdown();
        QString error;
        QVERIFY(!dialog.openTrajectory("included-v1", &error));
        QVERIFY(!dialog.clearTrajectory("downloaded-v1", &error));
        QVERIFY(!dialog.setCollection({}, &error));
        QVERIFY(!dialog.loadCatalog(QUrl::fromLocalFile(catalog_), &error));
        QCOMPARE(stopped.count(), 0);
        QTRY_COMPARE(stopped.count(), 1);
        QVERIFY(dialog.isShutdown());
        dialog.shutdown();
        QCOMPARE(stopped.count(), 1);
        verifySourceFiles();
    }

    void invalidCatalogIsVisible_data() {
        QTest::addColumn<QByteArray>("json");
        QTest::newRow("malformed-json") << QByteArray("not JSON");
        QTest::newRow("old-array-catalog") << QByteArray("[]");
        const QJsonObject entry{{"key", "sample-v1"}, {"title", "Sample"}, {"frames", 100},
                                {"archive", "packages/sample-v1.tar.xz"}, {"entry_point", "run.LGS"},
                                {"archive_bytes", 100}, {"expanded_bytes", 200}};
        QJsonObject root{{"kind", "collection"}, {"schema_version", 2},
                         {"title", "Reader collection"}, {"entries", QJsonArray{entry}}};
        QTest::newRow("wrong-version") << QJsonDocument(root).toJson();
        root["schema_version"] = 1;
        root.remove("title");
        QTest::newRow("missing-title") << QJsonDocument(root).toJson();
        root["title"] = "Reader collection";
        root["entries"] = QJsonArray();
        QTest::newRow("empty-entries") << QJsonDocument(root).toJson();
        root["entries"] = QJsonArray{QJsonObject()};
        QTest::newRow("invalid-entry") << QJsonDocument(root).toJson();
    }

    void invalidCatalogIsVisible() {
        QFETCH(QByteArray, json);
        const QString invalid = QFileInfo(catalog_).dir().filePath("invalid.lgs");
        writeFile(invalid, json);
        ReaderCollectionDialog dialog({}, nullptr, buffer_, installed_);
        QString error;
        QVERIFY(!dialog.loadCatalog(QUrl::fromLocalFile(invalid), &error));
        QVERIFY(!error.isEmpty());
        QCOMPARE(dialog.state()["error"].toString(), error);
        QVERIFY(dialog.state()["entries"].toArray().isEmpty());
        QVERIFY(!dialog.isBusy());
        const auto* status = dialog.findChild<QLabel*>("collectionStatus");
        QVERIFY(status);
        QCOMPARE(status->text(), error);
        QVERIFY(!dialog.openTrajectory("included-v1", &error));

        QVERIFY2(dialog.loadCatalog(QUrl::fromLocalFile(catalog_), &error), qPrintable(error));
        QVERIFY(dialog.state()["error"].toString().isEmpty());
        const auto entries = dialog.state()["entries"].toArray();
        QVERIFY(!dialog.loadCatalog(QUrl::fromLocalFile(invalid), &error));
        QCOMPARE(dialog.state()["entries"].toArray(), entries);
        QCOMPARE(dialog.state()["source"].toString(), collection_.source.toString());
        QCOMPARE(status->text(), error);
        QVERIFY(QFileInfo::exists(bufferedRun()));
        verifySourceFiles();
    }
};

int main(int argc, char** argv) {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication application(argc, argv);
    application.setQuitOnLastWindowClosed(false);
    TrajectoryLibraryTests tests;
    return QTest::qExec(&tests, argc, argv);
}
#include "trajectory_library_tests.moc"
