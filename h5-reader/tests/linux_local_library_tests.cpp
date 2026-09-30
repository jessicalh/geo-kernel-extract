#include "app/TrajectoryLibraryDialog.h"

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLineEdit>
#include <QPushButton>
#include <QSettings>
#include <QSignalSpy>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QtTest>

using h5reader::app::TrajectoryLibraryDialog;

class LinuxLocalLibraryTests final : public QObject {
    Q_OBJECT
    QTemporaryDir temporary_;
    QString sourceRoot_;
    QString copyRoot_;
    QString catalog_;
    QString sourceRelative_;
    QByteArray previousSource_;
    QByteArray previousCopies_;
    bool hadSource_ = false;
    bool hadCopies_ = false;

    static void writeFile(const QString& path, const QByteArray& contents) {
        QVERIFY(QDir().mkpath(QFileInfo(path).absolutePath()));
        QFile file(path);
        QVERIFY2(file.open(QIODevice::WriteOnly), qPrintable(file.errorString()));
        QCOMPARE(file.write(contents), contents.size());
    }

    QJsonObject entry(const QString& source = {}) const {
        return {{"key", "local-test-protein"},
                {"member", "test-protein"},
                {"dataset_id", "test-trajectory"},
                {"protein_id", "test-protein"},
                {"title", "Local test protein"},
                {"description", "A trajectory available without a network connection"},
                {"group", "MD"},
                {"frames", 1},
                {"dft_frames", 0},
                {"source_lgs", source.isEmpty() ? sourceRelative_ : source}};
    }

    void writeCatalog(const QJsonArray& entries) {
        writeFile(catalog_, QJsonDocument(QJsonObject{
            {"schema_version", 1}, {"kind", "h5reader-local-library"}, {"datasets", entries},
        }).toJson());
    }

    static int visibleRows(const QTableWidget* table) {
        int count = 0;
        for (int row = 0; row < table->rowCount(); ++row)
            if (!table->isRowHidden(row))
                ++count;
        return count;
    }

private slots:
    void initTestCase() {
        QVERIFY(temporary_.isValid());
        QCoreApplication::setOrganizationName("h5reader-tests");
        QCoreApplication::setApplicationName("linux-local-library");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, temporary_.path());
        hadSource_ = qEnvironmentVariableIsSet("H5READER_LOCAL_LIBRARY_ROOT");
        hadCopies_ = qEnvironmentVariableIsSet("H5READER_LOCAL_COPY_ROOT");
        previousSource_ = qgetenv("H5READER_LOCAL_LIBRARY_ROOT");
        previousCopies_ = qgetenv("H5READER_LOCAL_COPY_ROOT");
    }

    void init() {
        const QString test = QString::fromLatin1(QTest::currentTestFunction());
        sourceRoot_ = temporary_.filePath(test + "/source-drive");
        copyRoot_ = temporary_.filePath(test + "/local-copies");
        catalog_ = temporary_.filePath(test + "/linux-local-library.json");
        sourceRelative_ = QStringLiteral("01_Reader/datasets/MD_173/test-protein/source.lgs");
        QVERIFY(QDir().mkpath(sourceRoot_));
        qputenv("H5READER_LOCAL_LIBRARY_ROOT", sourceRoot_.toUtf8());
        qputenv("H5READER_LOCAL_COPY_ROOT", copyRoot_.toUtf8());
        writeCatalog({entry()});
        writeFile(QDir(sourceRoot_).filePath(sourceRelative_), "original source LGS bytes\n");
    }

    void cleanupTestCase() {
        if (hadSource_)
            qputenv("H5READER_LOCAL_LIBRARY_ROOT", previousSource_);
        else
            qunsetenv("H5READER_LOCAL_LIBRARY_ROOT");
        if (hadCopies_)
            qputenv("H5READER_LOCAL_COPY_ROOT", previousCopies_);
        else
            qunsetenv("H5READER_LOCAL_COPY_ROOT");
    }

    void fullCatalogHas176SelectableDatasetsAndSearch() {
        const QString publishedCatalog = QDir(QStringLiteral(H5READER_SOURCE_DIR))
                                             .filePath("resources/linux-local-library.json");
        TrajectoryLibraryDialog dialog(nullptr, publishedCatalog);
        const auto state = dialog.state();
        QVERIFY(state["local_library"].toBool());
        QVERIFY(state["open_in_place"].toBool());
        QVERIFY(state["cache_root"].toString().isEmpty());
        QVERIFY2(state["error"].toString().isEmpty(), qPrintable(state["error"].toString()));
        const auto rows = state["entries"].toArray();
        QCOMPARE(rows.size(), 176);
        int dftFrames = 0;
        int frames = 0;
        for (const auto& value : rows) {
            const auto row = value.toObject();
            QVERIFY(!row["source_lgs"].toString().isEmpty());
            frames += row["frames"].toInt();
            dftFrames += row["dft_frames"].toInt();
        }
        QCOMPARE(frames, 24553);
        QCOMPARE(dftFrames, 7253);
        auto* table = dialog.findChild<QTableWidget*>("trajectoryCatalog");
        auto* search = dialog.findChild<QLineEdit*>("trajectorySearch");
        QVERIFY(table);
        QVERIFY(search);
        QCOMPARE(table->rowCount(), 176);
        QCOMPARE(visibleRows(table), 176);
        search->setText("bmr10005");
        QCOMPARE(visibleRows(table), 1);
        search->setText("trpcage");
        QCOMPARE(visibleRows(table), 1);
        search->setText("chignolin");
        QCOMPARE(visibleRows(table), 1);
        search->setText("no-such-protein-in-the-library");
        QCOMPARE(visibleRows(table), 0);
        search->clear();
        QCOMPARE(visibleRows(table), 176);
    }

    void directSourceOpenNeedsNoDownloadOrWritableSource() {
        QVERIFY(QFile::setPermissions(QDir(sourceRoot_).filePath(sourceRelative_),
                                     QFileDevice::ReadOwner | QFileDevice::ReadGroup | QFileDevice::ReadOther));
        TrajectoryLibraryDialog dialog(nullptr, catalog_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY2(dialog.openTrajectory("local-test-protein", &error), qPrintable(error));
        QTRY_COMPARE(opened.size(), 1);
        QCOMPARE(opened.first().first().toString(), QDir(sourceRoot_).filePath(sourceRelative_));
        QVERIFY(!dialog.state()["busy"].toBool());
        QVERIFY(dialog.state()["cache_root"].toString().isEmpty());
        QVERIFY(dialog.state()["open_in_place"].toBool());
        QFile original(QDir(sourceRoot_).filePath(sourceRelative_));
        QVERIFY(original.open(QIODevice::ReadOnly));
        QCOMPARE(original.readAll(), QByteArray("original source LGS bytes\n"));
        QVERIFY(!QFileInfo::exists(QDir(copyRoot_).filePath("local-test-protein/run.LGS")));
    }

    void missingSourceReportsErrorWithoutStartingDownload() {
        QVERIFY(QFile::remove(QDir(sourceRoot_).filePath(sourceRelative_)));
        TrajectoryLibraryDialog dialog(nullptr, catalog_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY(!dialog.openTrajectory("local-test-protein", &error));
        QVERIFY(!error.isEmpty());
        QVERIFY(!dialog.state()["busy"].toBool());
        QCOMPARE(opened.size(), 0);
    }

    void sourceIsProtectedFromClear() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_);
        QString error;
        QVERIFY(!dialog.clearTrajectory("local-test-protein", &error));
        QVERIFY(!error.isEmpty());
        QFile original(QDir(sourceRoot_).filePath(sourceRelative_));
        QVERIFY(original.open(QIODevice::ReadOnly));
        QCOMPARE(original.readAll(), QByteArray("original source LGS bytes\n"));
    }

    void pathEscapesAreRejected_data() {
        QTest::addColumn<QString>("path");
        QTest::newRow("parent") << QStringLiteral("../outside.lgs");
        QTest::newRow("nested-parent") << QStringLiteral("01_Reader/../../outside.lgs");
        QTest::newRow("absolute") << QStringLiteral("/tmp/outside.lgs");
        QTest::newRow("URL") << QStringLiteral("file:///tmp/outside.lgs");
    }

    void pathEscapesAreRejected() {
        QFETCH(QString, path);
        writeCatalog({entry(path)});
        TrajectoryLibraryDialog dialog(nullptr, catalog_);
        QVERIFY(!dialog.state()["error"].toString().isEmpty());
        QVERIFY(dialog.state()["entries"].toArray().isEmpty());
        QString error;
        QVERIFY(!dialog.openTrajectory("local-test-protein", &error));
    }

    void symlinkOutsideSourceCannotBeOpened() {
        const QString source = QDir(sourceRoot_).filePath(sourceRelative_);
        const QString unrelated = temporary_.filePath("outside-source.lgs");
        writeFile(unrelated, "unrelated file\n");
        QVERIFY(QFile::remove(source));
        QVERIFY(QFile::link(unrelated, source));
        TrajectoryLibraryDialog dialog(nullptr, catalog_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY(!dialog.openTrajectory("local-test-protein", &error));
        QCOMPARE(opened.size(), 0);
        QVERIFY(!error.isEmpty());
    }

    void legacyCachedRunNeverShadowsExistingSource_data() {
        QTest::addColumn<bool>("withReceipt");
        QTest::newRow("no-receipt") << false;
        QTest::newRow("completed-receipt") << true;
    }

    void legacyCachedRunNeverShadowsExistingSource() {
        QFETCH(bool, withReceipt);
        const QString cached = QDir(copyRoot_).filePath("local-test-protein/run.LGS");
        writeFile(cached, "legacy cached run");
        if (withReceipt) {
            writeFile(QDir(copyRoot_).filePath("local-test-protein/local-copy.json"),
                      QJsonDocument(QJsonObject{{"schema_version", 1}, {"completed", true},
                          {"key", "local-test-protein"}, {"type", "h5reader-local-copy"}}).toJson());
        }
        TrajectoryLibraryDialog dialog(nullptr, catalog_);
        QSignalSpy opened(&dialog, &TrajectoryLibraryDialog::openRequested);
        QString error;
        QVERIFY2(dialog.openTrajectory("local-test-protein", &error), qPrintable(error));
        QTRY_COMPARE(opened.size(), 1);
        QCOMPARE(opened.first().first().toString(), QDir(sourceRoot_).filePath(sourceRelative_));
        const auto state = dialog.state();
        QVERIFY(state["local_library"].toBool());
        QVERIFY(state["open_in_place"].toBool());
        QVERIFY(state["cache_root"].toString().isEmpty());
        const auto row = state["entries"].toArray().first().toObject();
        QVERIFY(!row["downloaded"].toBool());
        QVERIFY(!row.contains("copied"));
        QCOMPARE(row["directory"].toString(), sourceRoot_);
        QCOMPARE(row["entry_point"].toString(), sourceRelative_);
        QVERIFY(!dialog.clearTrajectory("local-test-protein", &error));
        QVERIFY(QFileInfo::exists(cached));
        QVERIFY(QFileInfo::exists(QDir(sourceRoot_).filePath(sourceRelative_)));
    }

    void localChooserOnlyOffersOpen() {
        TrajectoryLibraryDialog dialog(nullptr, catalog_);
        QVERIFY(!dialog.findChild<QPushButton*>("copyTrajectory"));
        auto* open = dialog.findChild<QPushButton*>("openTrajectory");
        QVERIFY(open);
        QCOMPARE(open->text(), QStringLiteral("Open"));
        for (const char* name : {"clearDownload", "cancelDownload", "cacheLocation"}) {
            auto* control = dialog.findChild<QPushButton*>(name);
            QVERIFY(control);
            QVERIFY2(control->isHidden(), name);
        }
    }

};

int main(int argc, char** argv) {
    QApplication application(argc, argv);
    LinuxLocalLibraryTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "linux_local_library_tests.moc"
