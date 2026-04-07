// Viewer entry point.
//
// Follows the Qt skill startup sequence:
//   1. QApplication (must be first)
//   2. Library logging (OperationLog — the ONE logging system)
//   3. MainWindow
//   4. show()
//   5. Deferred auto-load via QTimer::singleShot (event loop must be running)
//   6. app.exec()

#include <QApplication>
#include <QSurfaceFormat>
#include <QMetaType>
#include <QTimer>
#include <QDir>
#include <QFileInfo>
#include <QCommandLineParser>
#include <QMessageBox>
#include <QVTKOpenGLNativeWidget.h>

#include "MainWindow.h"
#include "ComputeWorker.h"
#include "RestServer.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"

// udp_log shim: legacy viewer code calls this. Forwards to OperationLog.
#include <cstdarg>
extern "C" void udp_log(const char* fmt, ...) {
    char buf[2048];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    nmr::OperationLog::Info(nmr::LogViewer, "Viewer", buf);
}

int main(int argc, char* argv[]) {
    // 1. QApplication — must be first
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc, argv);
    app.setApplicationName("NMR Shielding Tensor Viewer");
    app.setApplicationVersion("2.0");

    // 2. Library environment + logging
    nmr::RuntimeEnvironment::Load();
    nmr::OperationLog::ConfigureUdp("127.0.0.1", 9998);
    nmr::OperationLog::SetChannelMask(nmr::LogAll);
    nmr::OperationLog::LogSessionStart();

    // 3. Parse command line (Qt-compliant)
    QCommandLineParser parser;
    parser.setApplicationDescription("NMR shielding tensor viewer");
    parser.addHelpOption();
    parser.addOption({{"d", "dir"}, "Initial directory for File > Open", "path"});
    parser.addOption({{"p", "pdb"}, "Load a PDB file on startup", "file"});
    parser.addOption({"protein", "Load a protein directory on startup", "dir"});
    parser.addOption({"rest-port", "REST API port (default 9147, 0 to disable)", "port", "9147"});
    parser.addPositionalArgument("path", "PDB file or protein directory to load");
    parser.process(app);

    QString initialDir = parser.value("dir");
    QString loadPdb = parser.value("pdb");
    QString loadProteinDir = parser.value("protein");

    // Bare positional argument: auto-detect file vs directory
    if (loadPdb.isEmpty() && loadProteinDir.isEmpty() && !parser.positionalArguments().isEmpty()) {
        QString path = parser.positionalArguments().first();
        QFileInfo fi(path);
        if (fi.isDir())
            loadProteinDir = path;
        else if (fi.isFile())
            loadPdb = path;
    }

    // 4. Register metatypes for cross-thread signals
    qRegisterMetaType<ViewerLoadRequest>("ViewerLoadRequest");
    qRegisterMetaType<ComputeResult>("ComputeResult");

    // 5. Create and show window
    MainWindow window(initialDir);
    QObject::connect(&app, &QCoreApplication::aboutToQuit, &window, &MainWindow::shutdown);
    window.show();

    // 5b. REST API server for programmatic control
    quint16 restPort = parser.value("rest-port").toUShort();
    RestServer* restServer = nullptr;
    if (restPort > 0) {
        restServer = new RestServer(&window, restPort, &app);
        if (restServer->actualPort() > 0) {
            window.setWindowTitle(QString("Protein Tensor Viewer — REST :%1")
                .arg(restServer->actualPort()));
        }
    }

    // 6. Validate paths before deferring load
    if (!loadPdb.isEmpty() && !QFileInfo::exists(loadPdb)) {
        QMessageBox::critical(&window, "File Not Found",
            QString("PDB file does not exist:\n%1").arg(loadPdb));
        return EXIT_FAILURE;
    }
    if (!loadProteinDir.isEmpty() && !QDir(loadProteinDir).exists()) {
        QMessageBox::critical(&window, "Directory Not Found",
            QString("Protein directory does not exist:\n%1").arg(loadProteinDir));
        return EXIT_FAILURE;
    }

    // Deferred auto-load — event loop must be running first
    if (!loadProteinDir.isEmpty()) {
        nmr::OperationLog::Info(nmr::LogViewer, "main",
            "auto-loading protein dir: " + loadProteinDir.toStdString());
        QTimer::singleShot(0, &window, [&window, loadProteinDir]() {
            window.loadProteinDir(loadProteinDir.toStdString());
        });
    } else if (!loadPdb.isEmpty()) {
        nmr::OperationLog::Info(nmr::LogViewer, "main",
            "auto-loading PDB: " + loadPdb.toStdString());
        QTimer::singleShot(0, &window, [&window, loadPdb]() {
            window.loadPdb(loadPdb.toStdString());
        });
    }

    nmr::OperationLog::Info(nmr::LogViewer, "main", "entering event loop");
    return app.exec();
}
