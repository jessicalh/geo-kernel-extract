// StructuredLogger — per-category logging with stderr + UDP output.
//
// Installed once at startup via Install(). Takes over Qt's global message
// handler (qInstallMessageHandler) so every qDebug/qInfo/qWarning/qCritical
// call is enriched with timestamp, thread name, category, file/line, then
// written to stderr (human-readable single line) AND to a UDP datagram
// (structured JSON, one datagram per message).
//
// The UDP stream is the primary debugging channel. udp_listen.py on port
// 9997 tails it during development. When the reader misbehaves, the log
// stream is consulted BEFORE code changes. See feedback_qt_discipline.
//
// Target host/port default to 127.0.0.1:9997. Override via environment
// variable H5READER_LOG_UDP="host:port".

#pragma once

#include <QObject>
#include <QHostAddress>
#include <QMutex>
#include <QString>
#include <QUdpSocket>

namespace h5reader::diagnostics {

class StructuredLogger final : public QObject {
    Q_OBJECT
public:
    // Call once from main() AFTER QCoreApplication exists. Subsequent calls
    // are no-ops.
    static void Install();

    // Singleton pointer. Null before Install().
    static StructuredLogger* Instance();

    // Emit one message. Called from the Qt global message handler and
    // directly from the diagnostics macros. Thread-safe — a mutex guards
    // the UDP socket.
    void Emit(QtMsgType type,
              const char* category,
              const QString& message,
              const char* file,
              int line,
              const char* function);

private:
    explicit StructuredLogger(const QHostAddress& host, quint16 port);
    ~StructuredLogger() override = default;

    QUdpSocket   udp_;
    QHostAddress host_;
    quint16      port_;
    QMutex       lock_;
};

}  // namespace h5reader::diagnostics
