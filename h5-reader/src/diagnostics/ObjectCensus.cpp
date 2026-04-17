#include "ObjectCensus.h"

#include <QMutexLocker>

#include <cstdio>
#include <cstring>

#ifndef _WIN32
#  include <unistd.h>
#endif

namespace h5reader::diagnostics {

namespace {
ObjectCensus* g_instance = nullptr;

// Signal-handler-safe write of a C string. The result is deliberately
// discarded — we have no useful recovery path during crash dumping.
void SafeWrite(int fd, const char* s) {
    if (!s) return;
#ifndef _WIN32
    const ssize_t rc = ::write(fd, s, std::strlen(s));
    (void)rc;
#else
    std::fputs(s, stderr);  // TODO: Win32 WriteFile when CrashHandler lands
#endif
}
}  // namespace

ObjectCensus* ObjectCensus::Instance() {
    if (!g_instance) g_instance = new ObjectCensus();
    return g_instance;
}

void ObjectCensus::Register(QObject* obj) {
    if (!obj) return;

    auto* inst = Instance();
    {
        QMutexLocker lk(&inst->lock_);
        inst->live_.insert(obj);
    }

    // Auto-remove on destruction. The lambda captures inst by value
    // because ObjectCensus lives for the process lifetime.
    QObject::connect(obj, &QObject::destroyed,
                     [inst](QObject* dying) {
                         QMutexLocker lk(&inst->lock_);
                         inst->live_.remove(dying);
                     });
}

void ObjectCensus::Dump(int fd) {
    // ASYNC-SIGNAL-SAFE: no mutex, no QString, no heap allocation.
    // We iterate the set without locking — on a crash the state may be
    // inconsistent, but reading addresses is cheap and addresses can be
    // matched later against memory maps or the live symbol table.
    auto* inst = g_instance;
    if (!inst) {
        SafeWrite(fd, "  (no census instance)\n");
        return;
    }

    int count = 0;
    for (const auto& obj : inst->live_) {
        if (!obj) continue;
        char buf[64];
        const int n = std::snprintf(buf, sizeof(buf), "  %p\n", (void*)obj);
        if (n > 0) {
#ifndef _WIN32
            const ssize_t rc = ::write(fd, buf, static_cast<size_t>(n));
            (void)rc;
#else
            std::fputs(buf, stderr);
#endif
        }
        ++count;
    }

    char tail[64];
    const int n = std::snprintf(tail, sizeof(tail),
                                 "  (%d live objects)\n", count);
    if (n > 0) {
#ifndef _WIN32
        const ssize_t rc = ::write(fd, tail, static_cast<size_t>(n));
        (void)rc;
#else
        std::fputs(tail, stderr);
#endif
    }
}

int ObjectCensus::LiveCount() {
    auto* inst = Instance();
    QMutexLocker lk(&inst->lock_);
    return inst->live_.size();
}

}  // namespace h5reader::diagnostics
