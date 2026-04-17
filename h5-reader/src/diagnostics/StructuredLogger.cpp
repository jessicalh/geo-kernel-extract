#include "StructuredLogger.h"

#include <QCoreApplication>
#include <QDateTime>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QMutexLocker>
#include <QString>
#include <QStringList>
#include <QThread>

#include <cstdio>
#include <cstdlib>

namespace h5reader::diagnostics {

namespace {
StructuredLogger* g_instance = nullptr;

// Qt's global message handler. Routes every qDebug/qInfo/qWarning/
// qCritical/qFatal call through the logger so we get consistent
// stderr + UDP output.
void QtMessageHandler(QtMsgType type,
                      const QMessageLogContext& ctx,
                      const QString& msg) {
    if (auto* inst = StructuredLogger::Instance()) {
        inst->Emit(type, ctx.category, msg, ctx.file, ctx.line, ctx.function);
    } else {
        // Before Install() — keep noise on stderr at least.
        std::fprintf(stderr, "%s\n", qUtf8Printable(msg));
    }
}

const char* SeverityTag(QtMsgType type) {
    switch (type) {
        case QtDebugMsg:    return "debug";
        case QtInfoMsg:     return "info";
        case QtWarningMsg:  return "warning";
        case QtCriticalMsg: return "error";
        case QtFatalMsg:    return "fatal";
    }
    return "?";
}
}  // namespace

void StructuredLogger::Install() {
    if (g_instance) return;

    QHostAddress host(QStringLiteral("127.0.0.1"));
    quint16      port = 9997;

    if (const char* env = std::getenv("H5READER_LOG_UDP")) {
        QString s = QString::fromUtf8(env);
        const int colon = s.lastIndexOf(':');
        if (colon > 0 && colon + 1 < s.size()) {
            QHostAddress h;
            const QString hostStr = s.left(colon);
            if (h.setAddress(hostStr)) host = h;
            bool ok = false;
            const quint16 p = s.mid(colon + 1).toUShort(&ok);
            if (ok) port = p;
        }
    }

    g_instance = new StructuredLogger(host, port);
    if (auto* app = QCoreApplication::instance()) {
        g_instance->setParent(app);   // lifetime follows QApp
    }

    qInstallMessageHandler(&QtMessageHandler);
    qInfo().noquote() << "StructuredLogger installed — UDP target"
                      << host.toString() << ":" << port;
}

StructuredLogger* StructuredLogger::Instance() { return g_instance; }

StructuredLogger::StructuredLogger(const QHostAddress& host, quint16 port)
    : host_(host), port_(port) {
    // No bind — writeDatagram() picks an ephemeral outbound port.
}

void StructuredLogger::Emit(QtMsgType type,
                            const char* category,
                            const QString& message,
                            const char* file,
                            int line,
                            const char* function) {
    const char* severity  = SeverityTag(type);
    const QString threadName =
        QThread::currentThread()->objectName().isEmpty()
            ? QStringLiteral("unnamed")
            : QThread::currentThread()->objectName();
    const QString cat = QString::fromUtf8(category ? category : "default");

    // Human-readable line to stderr.
    std::fprintf(stderr, "[%s] %s/%s: %s\n",
                 severity,
                 qUtf8Printable(threadName),
                 qUtf8Printable(cat),
                 qUtf8Printable(message));
    std::fflush(stderr);

    // Structured JSON datagram.
    QJsonObject obj;
    obj["ts"]       = QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs);
    obj["severity"] = QString::fromLatin1(severity);
    obj["category"] = cat;
    obj["thread"]   = threadName;
    obj["message"]  = message;
    if (file)     { obj["file"] = QString::fromUtf8(file); obj["line"] = line; }
    if (function)   obj["function"] = QString::fromUtf8(function);

    const QByteArray json = QJsonDocument(obj).toJson(QJsonDocument::Compact);
    {
        QMutexLocker lk(&lock_);
        const qint64 sent = udp_.writeDatagram(json, host_, port_);
        if (sent != json.size()) {
            std::fprintf(stderr,
                         "StructuredLogger: UDP send failed — sent=%lld "
                         "expected=%lld error=\"%s\"\n",
                         static_cast<long long>(sent),
                         static_cast<long long>(json.size()),
                         qUtf8Printable(udp_.errorString()));
        }
    }

    if (type == QtFatalMsg) std::abort();
}

}  // namespace h5reader::diagnostics
