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

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>

namespace h5reader::diagnostics {

namespace {
StructuredLogger* g_instance = nullptr;
std::atomic<std::uint32_t> g_categoryMask{kDefaultLogMask};

// Category-name to bit table. The keys match the strings passed to
// Q_LOGGING_CATEGORY in each module. Categories without a row here are
// emitted unconditionally — silent drops are worse than verbose logs.
//
// New categories that should be mask-gated: add a row here, declare the
// bit constant in StructuredLogger.h, and surface the symbolic name in
// the symbolic-names table below.
const std::unordered_map<std::string, std::uint32_t>& CategoryBitTable() {
    static const std::unordered_map<std::string, std::uint32_t> table = {
        // Render scheduler + EndEvent observer
        {"h5reader.scene",           kCatRender | kCatFrame},
        // Camera composer
        {"h5reader.composer",        kCatCamera},
        // Camera input filter (mouse / wheel)
        {"h5reader.caminput",        kCatCamera},
        // Atom picker
        {"h5reader.picker",          kCatPicker},
        // Transform decorator
        {"h5reader.transform",       kCatTransform},
        // REST surface (lifecycle + control)
        {"h5reader.rest",            kCatHealth},
        // Lifecycle (startup / shutdown)
        {"h5reader.lifecycle",       kCatHealth},
        // Main window
        {"h5reader.window",          kCatHealth},
        // Diagnostics infrastructure
        {"h5reader.threadguard",     kCatHealth},
        {"h5reader.connections",     kCatHealth},
        // Overlay categories
        {"h5reader.overlay.measurement", kCatOverlay | kCatPicker},
    };
    return table;
}

// Symbolic names for REST surface. Order matches the bit positions.
struct SymbolEntry { const char* name; std::uint32_t bit; };
const SymbolEntry kSymbolTable[] = {
    {"RENDER",    kCatRender},
    {"FRAME",     kCatFrame},
    {"CAMERA",    kCatCamera},
    {"SCRUB",     kCatScrub},
    {"OVERLAY",   kCatOverlay},
    {"PICKER",    kCatPicker},
    {"TRANSFORM", kCatTransform},
    {"HEALTH",    kCatHealth},
};

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
                      << host.toString() << ":" << port
                      << "| category mask = 0x" << QString::number(g_categoryMask.load(), 16);
}

StructuredLogger* StructuredLogger::Instance() { return g_instance; }

void StructuredLogger::SetCategoryMask(std::uint32_t mask) {
    g_categoryMask.store(mask);
}

std::uint32_t StructuredLogger::CategoryMask() {
    return g_categoryMask.load();
}

std::uint32_t StructuredLogger::MaskFromSymbolicNames(const QStringList& names) {
    std::uint32_t out = 0;
    for (const QString& n : names) {
        const QString upper = n.trimmed().toUpper();
        for (const auto& e : kSymbolTable) {
            if (upper == QLatin1String(e.name)) {
                out |= e.bit;
                break;
            }
        }
    }
    return out;
}

QStringList StructuredLogger::SymbolicNamesFromMask(std::uint32_t mask) {
    QStringList out;
    for (const auto& e : kSymbolTable) {
        if (mask & e.bit)
            out.append(QString::fromLatin1(e.name));
    }
    return out;
}

std::uint32_t LogCategoryMaskFor(const char* categoryName) {
    if (!categoryName) return 0;
    const auto& table = CategoryBitTable();
    auto it = table.find(std::string(categoryName));
    if (it == table.end()) return 0;
    return it->second;
}

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
    // Severity gate: warning/critical/fatal bypass the category mask.
    const bool alwaysEmit = (type == QtWarningMsg)
                             || (type == QtCriticalMsg)
                             || (type == QtFatalMsg);

    if (!alwaysEmit) {
        const std::uint32_t catBits = LogCategoryMaskFor(category);
        if (catBits != 0) {
            const std::uint32_t mask = g_categoryMask.load();
            if ((catBits & mask) == 0)
                return;   // gated out
        }
        // Categories without a bit mapping fall through (emit
        // unconditionally — silent drops are worse than verbose logs).
    }

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
