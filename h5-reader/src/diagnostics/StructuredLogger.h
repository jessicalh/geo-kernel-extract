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
//
// Category bitmask (per spec/viewport_pipeline_2026-05-30.md §8 decision 3):
// each Q_LOGGING_CATEGORY name maps to a single bit; the logger consults
// the mask at emit time and drops messages whose category bit is clear.
// ERROR / WARNING / CRITICAL severities bypass the gate (always emit).
// The mask is runtime-settable via SetCategoryMask() and the REST surface
// at /log/mask. Default mask is 0xFE — everything except RENDER (the
// volume one; turn on when debugging the render scheduler).
//
// Category bit table:
//   RENDER     0x01  every render trigger + EndEvent observer line; HIGH
//   FRAME      0x02  setFrame entry / fan-out / per-frame snapshot; MEDIUM
//   CAMERA     0x04  camera mode changes + per-frame absolute writes; MEDIUM
//   SCRUB      0x08  slider drag events; bursty during interaction
//   OVERLAY    0x10  visibility toggles + overlay refreshCurrentFrame; LOW-MED
//   PICKER     0x20  pick + selection events; LOW
//   TRANSFORM  0x40  TransformedConformation mode change + Kabsch recompute
//   HEALTH     0x80  startup / shutdown / REST lifecycle / snapshot success
//
// Q_LOGGING_CATEGORY names are mapped to bits in Categories.cpp (the
// LogCategoryMaskFor() helper). Categories without a mapping are
// emitted unconditionally (the conservative default — silent drops are
// worse than verbose logs).

#pragma once

#include <QHostAddress>
#include <QMutex>
#include <QObject>
#include <QString>
#include <QStringList>
#include <QUdpSocket>

#include <atomic>
#include <cstdint>

namespace h5reader::diagnostics {

// Bit positions; powers of two so consumers can OR them.
constexpr std::uint32_t kCatRender    = 0x01;
constexpr std::uint32_t kCatFrame     = 0x02;
constexpr std::uint32_t kCatCamera    = 0x04;
constexpr std::uint32_t kCatScrub     = 0x08;
constexpr std::uint32_t kCatOverlay   = 0x10;
constexpr std::uint32_t kCatPicker    = 0x20;
constexpr std::uint32_t kCatTransform = 0x40;
constexpr std::uint32_t kCatHealth    = 0x80;

// Default mask: everything except the high-volume RENDER category.
constexpr std::uint32_t kDefaultLogMask = 0xFE;

class StructuredLogger final : public QObject {
    Q_OBJECT
public:
    // Call once from main() AFTER QCoreApplication exists. Subsequent calls
    // are no-ops.
    static void Install();

    // Singleton pointer. Null before Install().
    static StructuredLogger* Instance();

    // Runtime category-mask setter / reader. Atomic so REST handlers
    // can flip it from any thread without disturbing in-flight Emit
    // calls. Default value is kDefaultLogMask = 0xFE.
    static void SetCategoryMask(std::uint32_t mask);
    static std::uint32_t CategoryMask();

    // Resolve a symbolic name list (["RENDER", "FRAME", ...]) into the
    // bitmask. Unknown names are skipped. Empty list returns 0.
    static std::uint32_t MaskFromSymbolicNames(const QStringList& names);

    // Inverse: bitmask -> symbolic name list. Useful for REST GET.
    static QStringList SymbolicNamesFromMask(std::uint32_t mask);

    // Emit one message. Called from the Qt global message handler and
    // directly from the diagnostics macros. Thread-safe — a mutex guards
    // the UDP socket. Filters by category bitmask except for warning/
    // critical/fatal severities (those always emit).
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

// Resolve a category name (the Q_LOGGING_CATEGORY's name argument, e.g.
// "h5reader.scene") to its mask bit. Returns 0 if the category isn't
// mapped — in which case the logger emits unconditionally (silent drops
// are worse than verbose logs).
std::uint32_t LogCategoryMaskFor(const char* categoryName);

}  // namespace h5reader::diagnostics
