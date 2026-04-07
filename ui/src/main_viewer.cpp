// Viewer entry point.
//
// Command line parsed via JobSpec (shared with nmr_extract).
// Qt-specific options (--rest-port) handled separately.
//
// Startup sequence:
//   1. QApplication (must be first)
//   2. Library logging (OperationLog — the ONE logging system)
//   3. Parse command line (JobSpec + Qt extras)
//   4. MainWindow
//   5. show()
//   6. Deferred auto-load via QTimer::singleShot (event loop must be running)
//   7. app.exec()

#include <QApplication>
#include <QSurfaceFormat>
#include <QMetaType>
#include <QTimer>
#include <QDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QVTKOpenGLNativeWidget.h>

#include "MainWindow.h"
#include "ComputeWorker.h"
#include "RestServer.h"
#include "JobSpec.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"
#include "CalculatorConfig.h"

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


// ============================================================================
// Extract --rest-port from argv before JobSpec parsing.
// JobSpec ignores unknown flags, but we need to pull this out for Qt.
// ============================================================================

static quint16 ExtractRestPort(int argc, char* argv[]) {
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--rest-port") == 0)
            return static_cast<quint16>(std::atoi(argv[i + 1]));
    }
    return 9147;  // default
}


int main(int argc, char* argv[]) {
    // 1. QApplication — must be first
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc, argv);
    app.setApplicationName("Protein Tensor Viewer");
    app.setApplicationVersion("2.0");

    // 2. Library environment + logging
    nmr::RuntimeEnvironment::Load();
    nmr::OperationLog::ConfigureUdp("127.0.0.1", 9998);
    nmr::OperationLog::SetChannelMask(nmr::LogAll);
    nmr::OperationLog::LogSessionStart();

    // 3. Parse command line via JobSpec
    //    --output is optional for the viewer (no output = just visualise).
    //    --rest-port is viewer-specific, extracted separately.
    auto spec = nmr::ParseJobSpec(argc, argv);
    quint16 restPort = ExtractRestPort(argc, argv);

    // Handle --help: print both JobSpec usage and viewer extras
    if (spec.mode == nmr::JobMode::None && spec.error.empty()) {
        nmr::PrintJobSpecUsage(argv[0]);
        fprintf(stderr,
            "\nViewer options:\n"
            "  --rest-port PORT   REST API port (default 9147, 0 to disable)\n"
            "\n"
            "  The viewer accepts all the same modes as nmr_extract.\n"
            "  --output is optional: omit it to just visualise.\n");
        return 0;
    }

    // Parse error
    if (!spec.Ok()) {
        fprintf(stderr, "ERROR: %s\n\n", spec.error.c_str());
        nmr::PrintJobSpecUsage(argv[0]);
        return 1;
    }

    // Validate file existence
    if (!nmr::ValidateJobSpec(spec)) {
        QMessageBox::critical(nullptr, "Validation Error",
            QString::fromStdString(spec.error));
        return 1;
    }

    // Print warnings
    for (const auto& w : spec.warnings)
        fprintf(stderr, "WARNING: %s\n", w.c_str());

    // Load calculator config if specified
    if (!spec.config_path.empty())
        nmr::CalculatorConfig::Load(spec.config_path);

    // 4. Register metatypes for cross-thread signals
    qRegisterMetaType<nmr::JobSpec>("nmr::JobSpec");
    qRegisterMetaType<ComputeResult>("ComputeResult");

    // 5. Create and show window
    MainWindow window;
    QObject::connect(&app, &QCoreApplication::aboutToQuit, &window, &MainWindow::shutdown);
    window.show();

    // 5b. REST API server
    RestServer* restServer = nullptr;
    if (restPort > 0) {
        restServer = new RestServer(&window, restPort, &app);
        if (restServer->actualPort() > 0) {
            window.setWindowTitle(QString("Protein Tensor Viewer — REST :%1")
                .arg(restServer->actualPort()));
        }
    }

    // 6. Deferred auto-load from JobSpec — all modes supported
    if (spec.mode != nmr::JobMode::None) {
        QTimer::singleShot(0, &window, [&window, spec]() {
            nmr::OperationLog::Info(nmr::LogViewer, "main",
                "auto-loading from JobSpec");
            window.loadFromJobSpec(spec);
        });
    }

    nmr::OperationLog::Info(nmr::LogViewer, "main", "entering event loop");
    return app.exec();
}
