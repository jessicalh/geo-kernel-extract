#include "ObjectCensus.h"

#include <QMutexLocker>

#include <cstdio>
#include <cstring>

#ifndef _WIN32
#  include <unistd.h>
#else
#  include <io.h>
#endif

namespace h5reader::diagnostics {

namespace {
ObjectCensus* g_instance = nullptr;

// Minimal crash-path write. Its result is discarded because recovery is not
// possible while producing the final report.
void SafeWrite(int fd, const char* s) {
    if (!s) return;
#ifndef _WIN32
    const ssize_t rc = ::write(fd, s, std::strlen(s));
    (void)rc;
#else
    const int rc = ::_write(fd, s, static_cast<unsigned int>(std::strlen(s)));
    (void)rc;
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
    // Do not take the mutex on the crash path. The state may be
    // inconsistent, but reading addresses is cheap and addresses can be
    // matched later against memory maps or the live symbol table.
    const auto* inst = g_instance;
    if (!inst) {
        SafeWrite(fd, "  (no census instance)\n");
        return;
    }

    int count = 0;
    for (const auto& obj : inst->live_) {
        if (!obj) continue;
        char buf[64];
        const int n = std::snprintf(buf, sizeof(buf), "  %p\n",
                                    static_cast<const void*>(obj));
        if (n > 0) {
#ifndef _WIN32
            const ssize_t rc = ::write(fd, buf, static_cast<size_t>(n));
            (void)rc;
#else
            const int rc = ::_write(fd, buf, static_cast<unsigned int>(n));
            (void)rc;
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
        const int rc = ::_write(fd, tail, static_cast<unsigned int>(n));
        (void)rc;
#endif
    }
}

}  // namespace h5reader::diagnostics
