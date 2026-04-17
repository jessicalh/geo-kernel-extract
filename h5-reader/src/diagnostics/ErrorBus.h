// ErrorBus — central sink for user-visible errors and degraded-state reports.
//
// Non-UI components (H5 loader, kernel evaluators, REST server, worker
// threads) call ErrorBus::Report(...) from any thread. The UI connects
// slots to errorReported() to populate the status bar, a modal dialog,
// or a persistent error panel.
//
// Every Report() also emits the same information through the
// StructuredLogger, so log consumers see errors in the UDP stream without
// the UI having to mirror them.

#pragma once

#include <QObject>
#include <QString>

namespace h5reader::diagnostics {

enum class Severity {
    Info,       // informational — status bar only
    Warning,    // degraded state — note in a panel, not a dialog
    Error,      // user action failed — dialog + panel
    Fatal,      // unrecoverable — dialog, then exit
};

class ErrorBus final : public QObject {
    Q_OBJECT
public:
    static ErrorBus* Instance();

    // Thread-safe — internally posts the signal via QueuedConnection when
    // called from a non-GUI thread. source and values may be empty.
    static void Report(Severity severity,
                       const QString& source,
                       const QString& message,
                       const QString& values = QString());

signals:
    void errorReported(h5reader::diagnostics::Severity severity,
                       QString source,
                       QString message,
                       QString values);

private:
    ErrorBus();
};

}  // namespace h5reader::diagnostics

Q_DECLARE_METATYPE(h5reader::diagnostics::Severity)
