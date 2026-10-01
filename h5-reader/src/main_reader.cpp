// h5-reader entry point. Startup order matters for crash capture and the VTK
// OpenGL context; the window itself may start empty or with one calcset.

#include "app/ReaderMainWindow.h"
#include "diagnostics/CrashHandler.h"
#include "diagnostics/ErrorBus.h"
#include "diagnostics/ObjectCensus.h"
#include "diagnostics/ShutdownSignals.h"
#include "diagnostics/StructuredLogger.h"

#include <QApplication>
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QFileInfo>
#include <QHostAddress>
#include <QLoggingCategory>
#include <QSurfaceFormat>
#include <QThread>


#include <QVTKOpenGLNativeWidget.h>

#include <vtkSMPTools.h>

#include <limits>

#ifdef _WIN32
#  include <QPalette>
#  include <QStyle>
#  include <QStyleFactory>
#  include <QStyleHints>
#endif

Q_LOGGING_CATEGORY(cLifecycle, "h5reader.lifecycle")

int main(int argc, char* argv[]) {
    // Install crash capture before QApplication construction.
    h5reader::diagnostics::CrashHandler::Install();

    // VTK requires desktop OpenGL; select it before QApplication.
    QCoreApplication::setAttribute(Qt::AA_UseDesktopOpenGL);

    // Request VTK's format after selecting the OpenGL backend.
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    QApplication app(argc, argv);
    app.setApplicationName(QStringLiteral("h5 reader"));
    app.setApplicationVersion(QStringLiteral(H5READER_VERSION));
    app.setOrganizationName(QStringLiteral("Beardsley Lab"));
    QThread::currentThread()->setObjectName(QStringLiteral("gui"));

    // Install structured logging before emitting application diagnostics.
    h5reader::diagnostics::StructuredLogger::Install();

#ifdef _WIN32
    // Keep enabled, checked, and disabled controls distinct on Windows.
    QApplication::setStyle(QStyleFactory::create(QStringLiteral("Fusion")));
    QGuiApplication::styleHints()->setColorScheme(Qt::ColorScheme::Light);

    QPalette pal;
    pal.setColor(QPalette::Window,          QColor(239, 239, 239));
    pal.setColor(QPalette::WindowText,      Qt::black);
    pal.setColor(QPalette::Base,            Qt::white);
    pal.setColor(QPalette::AlternateBase,   QColor(247, 247, 247));
    pal.setColor(QPalette::Text,            Qt::black);
    pal.setColor(QPalette::Button,          QColor(239, 239, 239));
    pal.setColor(QPalette::ButtonText,      Qt::black);
    pal.setColor(QPalette::BrightText,      Qt::red);
    pal.setColor(QPalette::Highlight,       QColor(48, 140, 198));
    pal.setColor(QPalette::HighlightedText, Qt::white);
    pal.setColor(QPalette::ToolTipBase,     QColor(255, 255, 225));
    pal.setColor(QPalette::ToolTipText,     Qt::black);
    pal.setColor(QPalette::Disabled, QPalette::WindowText, QColor(160, 160, 160));
    pal.setColor(QPalette::Disabled, QPalette::Text,       QColor(160, 160, 160));
    pal.setColor(QPalette::Disabled, QPalette::ButtonText, QColor(160, 160, 160));
    QApplication::setPalette(pal);
    qInfo().noquote() << "UI: Fusion + forced light palette (traditional) installed";
#endif

    // Enable VTK's threaded filters and report the backend actually selected.
    vtkSMPTools::SetBackend("STDThread");
    qInfo().noquote() << "VTK SMP backend:" << vtkSMPTools::GetBackend()
                      << "(max threads"
                      << vtkSMPTools::GetEstimatedNumberOfThreads() << ")";

    (void)h5reader::diagnostics::ErrorBus::Instance();
    (void)h5reader::diagnostics::ObjectCensus::Instance();

    qCInfo(cLifecycle).noquote() << "h5reader" << H5READER_VERSION << "starting" << "| Qt" << QT_VERSION_STR
                                 << "| thread=" << QThread::currentThread()->objectName();

    QCommandLineParser cli;
    cli.setApplicationDescription(QStringLiteral("Qt/VTK trajectory reader for nmr-extract analysis H5 files."));
    cli.addHelpOption();
    cli.addVersionOption();
    cli.addPositionalArgument(QStringLiteral("run_path"),
                              QStringLiteral("A calcset directory or .LGS calcset manifest."),
                              QStringLiteral("<run_path>"));
    const QCommandLineOption restOption(
        QStringLiteral("rest"),
        QStringLiteral("Start the embedded HTTP surface on <port>; "
                       "port 0 = kernel-pick (printed as H5READER_REST_PORT=NNNNN on stderr). "
                       "The window remains available until normal application shutdown."),
        QStringLiteral("port"));
    const QCommandLineOption restAddressOption(
        QStringLiteral("rest-address"),
        QStringLiteral("Address for the embedded HTTP surface. The default is 127.0.0.1; "
                       "use 0.0.0.0 for all IPv4 interfaces on a trusted network."),
        QStringLiteral("address"),
        QStringLiteral("127.0.0.1"));
    cli.addOption(restOption);
    cli.addOption(restAddressOption);
    cli.process(app);

    const bool runRest = cli.isSet(restOption);
    bool restPortOk = false;
    const uint restPortValue = runRest ? cli.value(restOption).toUInt(&restPortOk) : 0u;
    if (runRest && (!restPortOk || restPortValue > std::numeric_limits<quint16>::max())) {
        qCCritical(cLifecycle).noquote()
            << "--rest <port> must be an integer from 0 through 65535 (0 = kernel-pick)";
        return 1;
    }
    if (!runRest && cli.isSet(restAddressOption)) {
        qCCritical(cLifecycle).noquote() << "--rest-address requires --rest <port>";
        return 1;
    }
    QHostAddress restAddress;
    if (runRest && !restAddress.setAddress(cli.value(restAddressOption))) {
        qCCritical(cLifecycle).noquote()
            << "--rest-address must be a literal IPv4 or IPv6 address";
        return 1;
    }
    const quint16 restPort = static_cast<quint16>(restPortValue);

    const QStringList args = cli.positionalArguments();
    QString runPath;
    if (!args.isEmpty()) {
        runPath = args.first();
    } else if (runRest) {
        qCCritical(cLifecycle).noquote()
            << "No run path given. Usage: h5reader <calcset-dir | calcset.lgs>";
        return 1;
    }
    if (!runPath.isEmpty() && !QFileInfo::exists(runPath)) {
        qCCritical(cLifecycle).noquote() << "Run path not found:" << runPath;
        return 2;
    }

    // Finalize VTK while its OpenGL context is still valid.
    auto* window = new h5reader::app::ReaderMainWindow();
    QObject::connect(&app, &QCoreApplication::aboutToQuit, window, &h5reader::app::ReaderMainWindow::shutdown);
    h5reader::diagnostics::InstallShutdownSignalHandlers();

    if (!runPath.isEmpty() && !window->loadRunPath(runPath)) {
        qCCritical(cLifecycle).noquote() << "Load failed:" << window->lastLoadError();
        delete window;
        return 3;
    }

    // Show through the event queue so the event loop is live before first render.
    if (runRest) {
        QMetaObject::invokeMethod(window, [window, restAddress, restPort]() {
            window->show();
            qCInfo(cLifecycle).noquote() << "window shown for REST surface"
                                         << "| address=" << restAddress.toString()
                                         << "| port=" << restPort;
            const quint16 bound = window->startRestServer(restAddress, restPort);
            if (bound == 0) {
                qCCritical(cLifecycle).noquote()
                    << "REST server failed to bind; exiting";
                QCoreApplication::exit(6);
            }
        }, Qt::QueuedConnection);
    } else {
        QMetaObject::invokeMethod(window, [window]() {
            window->show();
            qCInfo(cLifecycle).noquote() << "window shown";
        }, Qt::QueuedConnection);
    }

    qCInfo(cLifecycle).noquote() << "entering event loop";
    const int rc = app.exec();
    qCInfo(cLifecycle).noquote() << "event loop exited with rc=" << rc;

    delete window;
    return rc;
}
