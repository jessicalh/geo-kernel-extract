// h5-reader — entry point.
//
// Step 3: open an H5, decode the typed model, construct the main
// window with the VTK viewport, show it. Playback controller animates
// the protein through the trajectory. Still no overlays, no atom
// inspector, no REST server — those ride on top of this substrate.
//
// Startup order (matches the existing viewer; see
// feedback_qt_discipline):
//   1. CrashHandler before QApplication (dump works even if QApp ctor
//      crashes).
//   2. QSurfaceFormat::setDefaultFormat before QApplication.
//   3. QApplication.
//   4. StructuredLogger (installs qInstallMessageHandler).
//   5. Warm ErrorBus and ObjectCensus singletons.
//   6. Parse CLI.
//   7. Load H5 via QtProteinLoader (synchronous on the GUI thread for
//      now; a worker isn't justified yet at typical H5 sizes).
//   8. Construct ReaderMainWindow, connect aboutToQuit → shutdown.
//   9. Deferred show() via QTimer::singleShot(0, ...) so the event loop
//      is running before any first render.
//  10. app.exec().

#include "app/ReaderMainWindow.h"
#include "diagnostics/CrashHandler.h"
#include "diagnostics/ErrorBus.h"
#include "diagnostics/ObjectCensus.h"
#include "diagnostics/ShutdownSignals.h"
#include "diagnostics/StructuredLogger.h"
#include "io/QtProteinLoader.h"

#include <QApplication>
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QFile>
#include <QFileInfo>
#include <QLoggingCategory>
#include <QSurfaceFormat>
#include <QThread>
#include <QTimer>


#include <QVTKOpenGLNativeWidget.h>

#include <vtkSMPTools.h>

#ifdef _WIN32
#  include <QPalette>
#  include <QStyleFactory>
#  include <QStyleHints>
#endif

Q_LOGGING_CATEGORY(cLifecycle, "h5reader.lifecycle")

int main(int argc, char* argv[]) {
    // 1. Crash handler before anything else.
    h5reader::diagnostics::CrashHandler::Install();

    // 2a. Force the desktop OpenGL driver before QApplication. Qt 6's
    //     auto-detect can fall back to ANGLE (software D3D-on-OpenGL)
    //     on Windows with some AMD driver combinations — and the
    //     symptom is exactly what we saw: VTK rendering CPU-bound,
    //     no GPU acceleration, "not notably faster" no matter what
    //     compile flags we add. Strix Halo's Radeon 8060S driver is
    //     modern; force the native path. Linux/macOS ignore this.
    //     See references/windows-gotchas.md §8.
    QCoreApplication::setAttribute(Qt::AA_UseDesktopOpenGL);

    // 2b. Surface format before QApplication. Must come AFTER
    //     setAttribute so the format request uses the chosen backend.
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    // 3. QApplication.
    QApplication app(argc, argv);
    app.setApplicationName(QStringLiteral("h5 reader"));
    app.setApplicationVersion(QStringLiteral(H5READER_VERSION));
    app.setOrganizationName(QStringLiteral("Beardsley Lab"));
    QThread::currentThread()->setObjectName(QStringLiteral("gui"));

    // 4. Structured logger. Install BEFORE any qInfo() we want routed
    //    through the stderr+UDP handler — earlier calls go through Qt's
    //    default handler which on Windows GUI builds may drop them.
    h5reader::diagnostics::StructuredLogger::Install();

#ifdef _WIN32
    // Windows 11's "modern" QStyle renders QStyle::SP_Media* icons in
    // the OS text colour. When the toolbar inherits the OS dark-mode
    // background but Qt's palette is still its default light Fusion
    // palette (or vice versa), the VCR icons collapse to black-on-black
    // — only the highlight on hover makes them visible. Fix: install
    // Fusion AND a matching palette so widgets and icons share one
    // colour scheme. We pick palette to match the OS preference.
    QApplication::setStyle(QStyleFactory::create(QStringLiteral("Fusion")));

    const bool osDark =
        QGuiApplication::styleHints()->colorScheme() == Qt::ColorScheme::Dark;
    QPalette pal;
    if (osDark) {
        // "Dark Fusion" — the standard Qt dark palette in the wild.
        pal.setColor(QPalette::Window,          QColor(53, 53, 53));
        pal.setColor(QPalette::WindowText,      Qt::white);
        pal.setColor(QPalette::Base,            QColor(35, 35, 35));
        pal.setColor(QPalette::AlternateBase,   QColor(53, 53, 53));
        pal.setColor(QPalette::ToolTipBase,     Qt::white);
        pal.setColor(QPalette::ToolTipText,     Qt::white);
        pal.setColor(QPalette::Text,            Qt::white);
        pal.setColor(QPalette::Button,          QColor(53, 53, 53));
        pal.setColor(QPalette::ButtonText,      Qt::white);
        pal.setColor(QPalette::BrightText,      Qt::red);
        pal.setColor(QPalette::Link,            QColor(42, 130, 218));
        pal.setColor(QPalette::Highlight,       QColor(42, 130, 218));
        pal.setColor(QPalette::HighlightedText, Qt::black);
        pal.setColor(QPalette::Disabled, QPalette::Text,       QColor(127, 127, 127));
        pal.setColor(QPalette::Disabled, QPalette::ButtonText, QColor(127, 127, 127));
    } else {
        // Fusion's default light palette is already fine; nothing to do.
        pal = QApplication::palette();
    }
    QApplication::setPalette(pal);
    qInfo().noquote() << "UI: Fusion style +"
                      << (osDark ? "dark" : "light")
                      << "palette installed";
#endif

    // VTK SMP backend: STDThread parallelises vtkStreamTracer + the
    // tube filter the butterfly overlay relies on, plus other VTK
    // filters used by the field-grid overlay. Default is "Sequential";
    // VTK ships with STDThread compiled in (see vtkSMP.h). Setting it
    // here before any filter executes lets every per-frame update use
    // all cores. SetBackend is silent on failure (unrecognised string
    // is a no-op), so read GetBackend() back and log it.
    vtkSMPTools::SetBackend("STDThread");
    qInfo().noquote() << "VTK SMP backend:" << vtkSMPTools::GetBackend()
                      << "(max threads"
                      << vtkSMPTools::GetEstimatedNumberOfThreads() << ")";

    // 5. Warm singletons. Install the Unix-signal bridge so Ctrl-C
    //    routes through QCoreApplication::quit(), which fires
    //    aboutToQuit → window->shutdown() (VTK finalise before GL
    //    context teardown). Without this, SIGINT hard-kills the
    //    process and VTK destructors touch dead GL state.
    (void)h5reader::diagnostics::ErrorBus::Instance();
    (void)h5reader::diagnostics::ObjectCensus::Instance();
    h5reader::diagnostics::InstallShutdownSignalHandlers();

    qCInfo(cLifecycle).noquote() << "h5reader" << H5READER_VERSION << "starting" << "| Qt" << QT_VERSION_STR
                                 << "| thread=" << QThread::currentThread()->objectName();

    // 6. Parse CLI.
    QCommandLineParser cli;
    cli.setApplicationDescription(QStringLiteral("Qt/VTK trajectory reader for nmr-extract analysis H5 files."));
    cli.addHelpOption();
    cli.addVersionOption();
    cli.addPositionalArgument(QStringLiteral("run_path"),
                              QStringLiteral("A run directory (trajectory or single-pose) "
                                             "or a trajectory.h5 file."),
                              QStringLiteral("<run_path>"));
    const QCommandLineOption restOption(
        QStringLiteral("rest"),
        QStringLiteral("Start embedded HTTP test surface on 127.0.0.1:<port>; "
                       "port 0 = kernel-pick (printed as H5READER_REST_PORT=NNNNN on stderr). "
                       "Window stays open until the process is signalled. "
                       "Replaces the retired --dashboard-path-smoke and --camera-plane-lock-smoke runners; "
                       "see h5-reader/tests/rest/ for the pytest suite that drives this surface."),
        QStringLiteral("port"));
    cli.addOption(restOption);
    cli.process(app);

    const bool runRest = cli.isSet(restOption);
    bool restPortOk = false;
    const quint16 restPort = static_cast<quint16>(
        runRest ? cli.value(restOption).toUInt(&restPortOk) : 0u);
    if (runRest && !restPortOk) {
        qCCritical(cLifecycle).noquote()
            << "--rest <port> must be a non-negative integer (0 = kernel-pick)";
        return 1;
    }

    const QStringList args = cli.positionalArguments();
    if (args.isEmpty()) {
        qCCritical(cLifecycle).noquote()
            << "No run path given. Usage: h5reader <run-dir | trajectory.h5>";
        return 1;
    }
    const QString& runPath = args.first();
    if (!QFileInfo::exists(runPath)) {
        qCCritical(cLifecycle).noquote() << "Run path not found:" << runPath;
        return 2;
    }

    // 7. Load — sniff trajectory vs single-pose by documented convention.
    qCInfo(cLifecycle).noquote() << "loading" << runPath;
    auto loaded = h5reader::io::QtProteinLoader::LoadRunPath(runPath);
    if (!loaded.ok) {
        qCCritical(cLifecycle).noquote() << "Load failed:" << loaded.error;
        return 3;
    }
    if (loaded.decodeWarnings > 0) {
        qCWarning(cLifecycle).noquote() << "Decode completed with" << loaded.decodeWarnings << "warnings";
    }

    // 8. Construct window. aboutToQuit → window->shutdown() runs the
    //    VTK-finalise-before-GL-context-destruction sequence.
    auto* window = new h5reader::app::ReaderMainWindow(std::move(loaded));
    QObject::connect(&app, &QCoreApplication::aboutToQuit, window, &h5reader::app::ReaderMainWindow::shutdown);

    // 9. Deferred show — event loop must be running before first render.
    if (runRest) {
        QTimer::singleShot(0, window, [window, restPort]() {
            window->show();
            qCInfo(cLifecycle).noquote() << "window shown for REST surface on port" << restPort;
            const quint16 bound = window->startRestServer(restPort);
            if (bound == 0) {
                qCCritical(cLifecycle).noquote()
                    << "REST server failed to bind; exiting";
                QCoreApplication::exit(6);
            }
        });
    } else {
        QTimer::singleShot(0, window, [window]() {
            window->show();
            qCInfo(cLifecycle).noquote() << "window shown";
        });
    }

    // 10. Event loop.
    qCInfo(cLifecycle).noquote() << "entering event loop";
    const int rc = app.exec();
    qCInfo(cLifecycle).noquote() << "event loop exited with rc=" << rc;

    // The window is Qt-managed; aboutToQuit fired its shutdown(). We
    // still delete it explicitly for cleanliness — ObjectCensus will
    // empty out on next dump.
    delete window;
    return rc;
}
