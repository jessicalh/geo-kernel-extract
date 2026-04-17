// ConnectionAuditor — audited signal/slot connections.
//
// ACONNECT(sender, signal, receiver, slot) wraps QObject::connect and logs
// the connection to the structured logger. The UDP stream then records the
// full wiring of the application at runtime — when a signal "doesn't fire",
// the log tells you whether the connection was ever made.
//
// Cleanup is Qt's — the connection auto-disconnects when either endpoint
// is destroyed. The auditor does NOT own lifetime; it only records the
// successful or failed connection at setup time.

#pragma once

#include <QLoggingCategory>
#include <QMetaObject>
#include <QObject>
#include <QString>

namespace h5reader::diagnostics {

Q_DECLARE_LOGGING_CATEGORY(cConnections)

namespace ConnectionAuditorImpl {

// The macro below calls this with a lambda that invokes QObject::connect
// with the original sender/signal/receiver/slot. This indirection lets
// the macro accept any signal/slot signature (including lambdas and
// member-pointer-taken method references) without needing the auditor
// to be a template on those types.
template <class Fn>
QMetaObject::Connection LogAndConnect(QObject* sender,
                                      const char* senderName,
                                      const char* signalText,
                                      QObject* receiver,
                                      const char* receiverName,
                                      const char* slotText,
                                      const char* file,
                                      int line,
                                      Fn&& makeConnection) {
    const auto conn = makeConnection();
    if (conn) {
        qCInfo(cConnections).nospace().noquote()
            << "[aconnect] "
            << (senderName ? senderName : "?")
            << "::" << (signalText ? signalText : "?")
            << "  ->  "
            << (receiverName ? receiverName : "?")
            << "::" << (slotText ? slotText : "?")
            << "  @ " << (file ? file : "?") << ":" << line;
    } else {
        qCCritical(cConnections).nospace().noquote()
            << "[aconnect FAILED] "
            << (senderName ? senderName : "?")
            << "::" << (signalText ? signalText : "?")
            << "  ->  "
            << (receiverName ? receiverName : "?")
            << "::" << (slotText ? slotText : "?")
            << "  @ " << (file ? file : "?") << ":" << line;
    }
    return conn;
}

inline const char* ClassName(const QObject* o) {
    return o ? o->metaObject()->className() : "nullptr";
}

}  // namespace ConnectionAuditorImpl

}  // namespace h5reader::diagnostics

// ACONNECT(sender, signal, receiver, slot [, type])
// Logs the connection to the structured logger. Returns the
// QMetaObject::Connection so callers may disconnect if needed.
#define ACONNECT(sender, signal, receiver, slot) \
    ::h5reader::diagnostics::ConnectionAuditorImpl::LogAndConnect( \
        (sender), \
        ::h5reader::diagnostics::ConnectionAuditorImpl::ClassName(sender), \
        #signal, \
        (receiver), \
        ::h5reader::diagnostics::ConnectionAuditorImpl::ClassName(receiver), \
        #slot, \
        __FILE__, __LINE__, \
        [&]() { return QObject::connect((sender), (signal), (receiver), (slot)); })
