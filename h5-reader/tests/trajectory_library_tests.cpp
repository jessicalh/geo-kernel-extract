#include "app/TrajectoryLibraryDialog.h"
#include "app/ReaderMainWindow.h"

#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLineEdit>
#include <QPushButton>
#include <QProgressBar>
#include <QSignalSpy>
#include <QSettings>
#include <QSurfaceFormat>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QVTKOpenGLNativeWidget.h>
#include <QtTest>

using h5reader::app::TrajectoryLibraryDialog;

class TrajectoryLibraryTests final : public QObject {
    Q_OBJECT
    QTemporaryDir directory_;
    QString catalog_;
    QString installed_;
    QString cache_;

    static void writeFile(const QString& path, const QByteArray& data) {
        QVERIFY(QDir().mkpath(QFileInfo(path).absolutePath()));
        QFile file(path);
        QVERIFY(file.open(QIODevice::WriteOnly));
        QCOMPARE(file.write(data), data.size());
        file.close();
        QCOMPARE(file.error(), QFileDevice::NoError);
    }

private slots:
    void initTestCase() {
        QCoreApplication::setOrganizationName("h5reader-tests");
        QCoreApplication::setApplicationName("trajectory-library");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory_.path());
    }

    void init() {
        QVERIFY(directory_.isValid());
        const QString test = QString::fromLatin1(QTest::currentTestFunction());
        catalog_ = directory_.filePath(test + "/trajectories.json");
        installed_ = directory_.filePath(test + "/installed");
        cache_ = directory_.filePath(test + "/cache");
        QJsonArray entries;
        for (const QString key : {QStringLiteral("included"), QStringLiteral("downloaded")}) {
            entries.append(QJsonObject{{"key", key},
                                       {"title", key},
                                       {"description", "Test trajectory"},
                                       {"url", "https://semantic.construction/files/sample.zip"},
                                       {"frames", 100},
                                       {"entry_point", "run.LGS"},
                                       {"archive_bytes", 100},
                                       {"expanded_bytes", 2}});
        }
        writeFile(catalog_, QJsonDocument(entries).toJson());
        writeFile(QDir(installed_).filePath("included/run.LGS"), "{}");
        writeFile(QDir(cache_).filePath("downloaded/run.LGS"), "{}");
    }

    void installedExampleIsOfflineAndProtected() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY(dialog.openTrajectory("included", &error));
        QCOMPARE(opened.count(), 1);
        QCOMPARE(opened.first().first().toString(), QDir(installed_).filePath("included/run.LGS"));
        QVERIFY(!dialog.clearTrajectory("included", &error));
        QVERIFY(error.contains("installation"));
        QVERIFY(QFileInfo::exists(QDir(installed_).filePath("included/run.LGS")));
        QVERIFY(!dialog.state()["busy"].toBool());
    }

    void cachedExampleUsesTheSameOpenSignal() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY(dialog.openTrajectory("downloaded", &error));
        QTRY_COMPARE(opened.count(), 1);
        QCOMPARE(opened.first().first().toString(), QDir(cache_).filePath("downloaded/run.LGS"));
        QVERIFY(dialog.state()["error"].toString().isEmpty());
    }

    void loadedDataCannotBeCleared() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        auto* table = dialog.findChild<QTableWidget*>("trajectoryCatalog");
        auto* clear = dialog.findChild<QPushButton*>("clearDownload");
        QVERIFY(table);
        QVERIFY(clear);
        QVERIFY(!clear->isEnabled());
        table->selectRow(1);
        QVERIFY(clear->isEnabled());
        dialog.setCurrentRun(QDir(cache_).filePath("downloaded"));
        QVERIFY(!clear->isEnabled());
        QString error;
        QVERIFY(!dialog.clearTrajectory("downloaded", &error));
        QVERIFY(error.contains("different trajectory"));
        dialog.setCurrentRun(QDir(installed_).filePath("included"));
        QVERIFY(clear->isEnabled());
        QVERIFY(dialog.clearTrajectory("downloaded", &error));
        QTRY_VERIFY(!dialog.state()["busy"].toBool());
        QVERIFY(!clear->isEnabled());
        QVERIFY(!QFileInfo::exists(QDir(cache_).filePath("downloaded")));
        QVERIFY(QFileInfo::exists(QDir(installed_).filePath("included/run.LGS")));
    }

    void cancellingPendingOpenDoesNotLoadTheTrajectory() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY(dialog.openTrajectory("downloaded", &error));
        dialog.cancel();
        const auto* progress = dialog.findChild<QProgressBar*>();
        QVERIFY(progress);
        QCOMPARE(progress->minimum(), 0);
        QCOMPARE(progress->maximum(), 0);
        const auto* cancel = dialog.findChild<QPushButton*>("cancelDownload");
        QVERIFY(cancel);
        QVERIFY(!cancel->isEnabled());
        dialog.shutdown();
        QTRY_VERIFY(dialog.isShutdown());
        QCOMPARE(opened.count(), 0);
        QVERIFY(QFileInfo::exists(QDir(cache_).filePath("downloaded/run.LGS")));
    }

    void readerCannotLoadAnEntryPendingClear() {
        QSettings().setValue("viewer/library/cacheRoot", cache_);
        h5reader::app::ReaderMainWindow window;
        auto* library = window.trajectoryLibrary();
        const auto entries = library->state()["entries"].toArray();
        QString path;
        QString key;
        for (const auto& value : entries) {
            const auto entry = value.toObject();
            if (!entry["included"].toBool()) {
                key = entry["key"].toString();
                path = QDir(cache_).filePath(key + "/" + entry["entry_point"].toString());
                break;
            }
        }
        QVERIFY(!path.isEmpty());
        writeFile(path, "{}");
        QString error;
        QVERIFY(library->clearTrajectory(key, &error));
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

    void shutdownDoesNotOpenAQueuedCacheHit() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY(dialog.openTrajectory("downloaded", &error));
        dialog.shutdown();
        QTRY_VERIFY(dialog.isShutdown());
        QCOMPARE(opened.count(), 0);
    }

    void searchFiltersRows() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        auto* search = dialog.findChild<QLineEdit*>("trajectorySearch");
        auto* table = dialog.findChild<QTableWidget*>("trajectoryCatalog");
        QVERIFY(search);
        QVERIFY(table);
        QCOMPARE(table->rowCount(), 2);
        search->setText("INCLUDED");
        QVERIFY(!table->isRowHidden(0));
        QVERIFY(table->isRowHidden(1));
        search->clear();
        QVERIFY(!table->isRowHidden(1));
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        auto* open = dialog.findChild<QPushButton*>("openTrajectory");
        QVERIFY(open);
        search->setText("downloaded");
        QVERIFY(!open->isEnabled());
        open->click();
        QCOMPARE(opened.count(), 0);
        table->selectRow(1);
        QVERIFY(open->isEnabled());
        search->setText("no matches");
        QVERIFY(!open->isEnabled());
    }

    void choosingCurrentCacheDoesNotNestDirectories() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        auto* button = dialog.findChild<QPushButton*>("cacheLocation");
        QVERIFY(button);
        button->click();
        auto* chooser = dialog.findChild<QFileDialog*>();
        QVERIFY(chooser);
        QVERIFY(QMetaObject::invokeMethod(chooser, "fileSelected", Q_ARG(QString, cache_)));
        QCOMPARE(dialog.state()["cache_root"].toString(), QDir(cache_).absolutePath());
        QVERIFY(!dialog.state()["shutdown"].toBool());
        chooser->reject();
    }

    void shutdownRejectsNewWork() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        QSignalSpy stopped(&dialog, &TrajectoryLibraryDialog::shutdownFinished);
        dialog.shutdown();
        QString error;
        QVERIFY(!dialog.openTrajectory("included", &error));
        QVERIFY(!dialog.clearTrajectory("downloaded", &error));
        QCOMPARE(stopped.count(), 0);
        QTRY_COMPARE(stopped.count(), 1);
        QVERIFY(dialog.isShutdown());
        dialog.shutdown();
        QCOMPARE(stopped.count(), 1);
    }

    void invalidCatalogIsVisible() {
        writeFile(catalog_, "not JSON");
        TrajectoryLibraryDialog dialog(nullptr, catalog_, cache_, installed_);
        QVERIFY(!dialog.state()["error"].toString().isEmpty());
        QVERIFY(dialog.state()["entries"].toArray().isEmpty());
        QString error;
        QVERIFY(!dialog.openTrajectory("included", &error));
        QVERIFY(!error.isEmpty());
    }
};

int main(int argc, char** argv) {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication application(argc, argv);
    TrajectoryLibraryTests tests;
    return QTest::qExec(&tests, argc, argv);
}
#include "trajectory_library_tests.moc"
