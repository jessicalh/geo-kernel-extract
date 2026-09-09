// ObjectCensus — global registry of live QObjects, for crash-time dumps.
//
// Registered QObject addresses are appended to the companion crash report.
//
// Objects are automatically removed from the registry when they emit
// QObject::destroyed, so the set reflects live objects up to Qt's normal
// lifetime guarantees.
//
// Dump deliberately avoids the mutex because a crash may interrupt a registry
// update. The resulting address list is best effort; class-name lookup is
// omitted because virtual dispatch is unsafe after memory corruption.

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

    // Dump live object addresses to a file descriptor without locking.
    static void Dump(int fd);

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
