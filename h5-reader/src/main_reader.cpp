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
#include "diagnostics/StructuredLogger.h"
#include "diagnostics/ShutdownSignals.h"
#include "io/QtProteinLoader.h"

#include <QApplication>
#include <QCommandLineParser>
#include <QFile>
#include <QLoggingCategory>
#include <QSurfaceFormat>
#include <QThread>
#include <QTimer>

#include <QVTKOpenGLNativeWidget.h>

Q_LOGGING_CATEGORY(cLifecycle, "h5reader.lifecycle")

int main(int argc, char* argv[]) {
    // 1. Crash handler before anything else.
    h5reader::diagnostics::CrashHandler::Install();

    // 2. Surface format before QApplication.
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    // 3. QApplication.
    QApplication app(argc, argv);
    app.setApplicationName(QStringLiteral("h5 reader"));
    app.setApplicationVersion(QStringLiteral(H5READER_VERSION));
    app.setOrganizationName(QStringLiteral("Beardsley Lab"));
    QThread::currentThread()->setObjectName(QStringLiteral("gui"));

    // 4. Structured logger.
    h5reader::diagnostics::StructuredLogger::Install();

    // 5. Warm singletons. Install the Unix-signal bridge so Ctrl-C
    //    routes through QCoreApplication::quit(), which fires
    //    aboutToQuit → window->shutdown() (VTK finalise before GL
    //    context teardown). Without this, SIGINT hard-kills the
    //    process and VTK destructors touch dead GL state.
    (void)h5reader::diagnostics::ErrorBus::Instance();
    (void)h5reader::diagnostics::ObjectCensus::Instance();
    h5reader::diagnostics::InstallShutdownSignalHandlers();

    qCInfo(cLifecycle).noquote()
        << "h5reader" << H5READER_VERSION << "starting"
        << "| Qt" << QT_VERSION_STR
        << "| thread=" << QThread::currentThread()->objectName();

    // 6. Parse CLI.
    QCommandLineParser cli;
    cli.setApplicationDescription(
        QStringLiteral("Qt/VTK trajectory reader for nmr-extract analysis H5 files."));
    cli.addHelpOption();
    cli.addVersionOption();
    cli.addPositionalArgument(QStringLiteral("h5_path"),
        QStringLiteral("Path to an analysis H5 file."),
        QStringLiteral("<h5_path>"));
    cli.process(app);

    const QStringList args = cli.positionalArguments();
    if (args.isEmpty()) {
        qCCritical(cLifecycle).noquote()
            << "No H5 path given. Usage: h5reader <path/to/*_analysis.h5>";
        return 1;
    }
    const QString h5Path = args.first();
    if (!QFile::exists(h5Path)) {
        qCCritical(cLifecycle).noquote() << "H5 file not found:" << h5Path;
        return 2;
    }

    // 7. Load.
    qCInfo(cLifecycle).noquote() << "loading" << h5Path;
    auto loaded = h5reader::io::QtProteinLoader::Load(h5Path);
    if (!loaded.ok) {
        qCCritical(cLifecycle).noquote() << "Load failed:" << loaded.error;
        return 3;
    }
    if (loaded.decodeErrors > 0) {
        qCWarning(cLifecycle).noquote()
            << "Decode completed with" << loaded.decodeErrors << "warnings";
    }

    // 8. Construct window. aboutToQuit → window->shutdown() runs the
    //    VTK-finalise-before-GL-context-destruction sequence.
    auto* window = new h5reader::app::ReaderMainWindow(std::move(loaded));
    QObject::connect(&app, &QCoreApplication::aboutToQuit,
                     window, &h5reader::app::ReaderMainWindow::shutdown);

    // 9. Deferred show — event loop must be running before first render.
    QTimer::singleShot(0, window, [window]() {
        window->show();
        qCInfo(cLifecycle).noquote() << "window shown";
    });

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
