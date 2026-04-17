#include "ErrorBus.h"

#include <QCoreApplication>
#include <QLoggingCategory>
#include <QMetaType>
#include <QThread>

namespace h5reader::diagnostics {

Q_LOGGING_CATEGORY(cErrorBus, "h5reader.errors")

namespace {
ErrorBus* g_instance = nullptr;
}  // namespace

ErrorBus* ErrorBus::Instance() {
    if (!g_instance) {
        g_instance = new ErrorBus();
        if (auto* app = QCoreApplication::instance()) {
            g_instance->setParent(app);
        }
    }
    return g_instance;
}

ErrorBus::ErrorBus() {
    qRegisterMetaType<Severity>("h5reader::diagnostics::Severity");
}

void ErrorBus::Report(Severity severity,
                      const QString& source,
                      const QString& message,
                      const QString& values) {
    // Log first so the UDP stream sees every error, even if no UI is
    // connected yet.
    switch (severity) {
        case Severity::Info:
            qCInfo(cErrorBus).noquote()
                << "[" << source << "]" << message
                << (values.isEmpty() ? QString() : QStringLiteral(" | ") + values);
            break;
        case Severity::Warning:
            qCWarning(cErrorBus).noquote()
                << "[" << source << "]" << message
                << (values.isEmpty() ? QString() : QStringLiteral(" | ") + values);
            break;
        case Severity::Error:
        case Severity::Fatal:
            qCCritical(cErrorBus).noquote()
                << "[" << source << "]" << message
                << (values.isEmpty() ? QString() : QStringLiteral(" | ") + values);
            break;
    }

    // Marshal the signal to the instance's thread. QMetaObject::invokeMethod
    // with QueuedConnection is correct regardless of the caller's thread.
    auto* inst = Instance();
    QMetaObject::invokeMethod(
        inst,
        [inst, severity, source, message, values]() {
            emit inst->errorReported(severity, source, message, values);
        },
        Qt::QueuedConnection);
}

}  // namespace h5reader::diagnostics
