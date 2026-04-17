// ObjectCensus — global registry of live QObjects, for crash-time dumps.
//
// Every QObject constructor should register itself with CENSUS_REGISTER(this).
// On crash, the CrashHandler dumps the set to the companion .txt file so a
// reviewer sees which classes were alive at the crash point.
//
// Objects are automatically removed from the registry when they emit
// QObject::destroyed, so the set reflects live objects up to Qt's normal
// lifetime guarantees.
//
// The signal-handler Dump path does NOT take the mutex (signals can fire
// on any thread at any point) — it reads addresses best-effort. Classnames
// via vtable access would risk a nested crash and are omitted.

#pragma once

#include <QObject>
#include <QSet>
#include <QMutex>

namespace h5reader::diagnostics {

class ObjectCensus final {
public:
    static ObjectCensus* Instance();

    // Register a live object. Safe to call from any thread.
    static void Register(QObject* obj);

    // Dump live object addresses to a file descriptor. Async-signal-safe
    // best-effort — does NOT take the mutex. Called from CrashHandler.
    static void Dump(int fd);

    // Count live objects. Takes the mutex; GUI-safe.
    static int LiveCount();

private:
    ObjectCensus() = default;
    ~ObjectCensus() = default;

    QMutex          lock_;
    QSet<QObject*>  live_;
};

}  // namespace h5reader::diagnostics

// Usage: CENSUS_REGISTER(this); in every QObject constructor.
#define CENSUS_REGISTER(obj) \
    ::h5reader::diagnostics::ObjectCensus::Register(obj)
