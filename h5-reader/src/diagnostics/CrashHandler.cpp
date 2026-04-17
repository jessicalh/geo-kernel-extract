#include "CrashHandler.h"
#include "ObjectCensus.h"

#include <QByteArray>
#include <QDir>
#include <QStandardPaths>
#include <QString>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef __linux__
#  include <execinfo.h>
#  include <fcntl.h>
#  include <signal.h>
#  include <sys/types.h>
#  include <time.h>
#  include <unistd.h>
#  ifndef PATH_MAX
#    define PATH_MAX 4096
#  endif
#elif defined(__APPLE__)
#  include <execinfo.h>
#  include <fcntl.h>
#  include <signal.h>
#  include <sys/types.h>
#  include <sys/syslimits.h>
#  include <time.h>
#  include <unistd.h>
#endif

namespace h5reader::diagnostics {

namespace {
// Populated at Install() time. QString allocations are NOT signal-safe;
// we must pre-resolve the dump directory into a plain C buffer.
#if defined(__linux__) || defined(__APPLE__)
char g_dumpDir[PATH_MAX] = {0};
#else
char g_dumpDir[4096] = {0};
#endif

std::atomic<bool> g_installed{false};

#if defined(__linux__) || defined(__APPLE__)
// Async-signal-safe write of a C string. Return value is intentionally
// ignored — no useful recovery path in the crash handler.
void SafeWrite(int fd, const char* s) {
    if (!s) return;
    const ssize_t rc = ::write(fd, s, std::strlen(s));
    (void)rc;
}

void SafeWriteInt(int fd, long v) {
    char buf[32];
    const int n = std::snprintf(buf, sizeof(buf), "%ld", v);
    if (n > 0) {
        const ssize_t rc = ::write(fd, buf, static_cast<size_t>(n));
        (void)rc;
    }
}

const char* SignalName(int signum) {
    switch (signum) {
        case SIGSEGV: return "SIGSEGV";
        case SIGBUS:  return "SIGBUS";
        case SIGABRT: return "SIGABRT";
        case SIGILL:  return "SIGILL";
        case SIGFPE:  return "SIGFPE";
        default:      return "UNKNOWN";
    }
}

void PosixHandler(int signum, siginfo_t* info, void* /*ucontext*/) {
    // Compose filename: {g_dumpDir}/crash_{pid}_{time}.txt
    char path[PATH_MAX];
    const long now = static_cast<long>(::time(nullptr));
    const long pid = static_cast<long>(::getpid());

    int n = std::snprintf(path, sizeof(path),
                          "%s/crash_%ld_%ld.txt",
                          g_dumpDir[0] ? g_dumpDir : "/tmp",
                          pid, now);
    if (n <= 0) { path[0] = '\0'; }

    int fd = ::open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) fd = STDERR_FILENO;

    SafeWrite(fd, "h5reader crash: signal ");
    SafeWriteInt(fd, signum);
    SafeWrite(fd, " (");
    SafeWrite(fd, SignalName(signum));
    SafeWrite(fd, ") at addr ");
    char addr[32];
    const int an = std::snprintf(addr, sizeof(addr), "%p",
                                  info ? info->si_addr : nullptr);
    if (an > 0) {
        const ssize_t rc = ::write(fd, addr, static_cast<size_t>(an));
        (void)rc;
    }
    SafeWrite(fd, "\npid=");
    SafeWriteInt(fd, pid);
    SafeWrite(fd, "\n\nBacktrace:\n");

    constexpr int MAX_FRAMES = 64;
    void* frames[MAX_FRAMES];
    const int nFrames = ::backtrace(frames, MAX_FRAMES);
    ::backtrace_symbols_fd(frames, nFrames, fd);

    SafeWrite(fd, "\nObjectCensus snapshot (best effort):\n");
    ObjectCensus::Dump(fd);

    SafeWrite(fd, "\nEnd of crash report.\n");

    if (fd != STDERR_FILENO) {
        ::close(fd);
        SafeWrite(STDERR_FILENO, "h5reader: crash dump written to ");
        SafeWrite(STDERR_FILENO, path);
        SafeWrite(STDERR_FILENO, "\n");
    }

    // Re-raise with default disposition so the OS/debugger sees the real
    // signal and a core dump is produced if enabled.
    ::signal(signum, SIG_DFL);
    ::raise(signum);
}
#endif  // __linux__ || __APPLE__
}  // namespace

void CrashHandler::SetDumpDirectory(const QString& path) {
    QDir().mkpath(path);   // create if missing
    const QByteArray ba = path.toUtf8();
    const size_t cap = sizeof(g_dumpDir);
    if (static_cast<size_t>(ba.size()) + 1 > cap) {
        std::fprintf(stderr,
                     "CrashHandler: dump path too long (%zu bytes, max %zu)\n",
                     static_cast<size_t>(ba.size()), cap - 1);
        return;
    }
    std::memcpy(g_dumpDir, ba.constData(), ba.size());
    g_dumpDir[ba.size()] = '\0';
}

void CrashHandler::Install() {
    if (g_installed.exchange(true)) return;

    if (g_dumpDir[0] == '\0') {
        const QString defaultDir =
            QStandardPaths::writableLocation(QStandardPaths::GenericDataLocation)
            + "/h5reader/crashes";
        SetDumpDirectory(defaultDir);
    }

#if defined(__linux__) || defined(__APPLE__)
    struct sigaction sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sa_sigaction = &PosixHandler;
    sa.sa_flags = SA_SIGINFO | SA_RESETHAND;
    sigemptyset(&sa.sa_mask);

    sigaction(SIGSEGV, &sa, nullptr);
    sigaction(SIGBUS,  &sa, nullptr);
    sigaction(SIGILL,  &sa, nullptr);
    sigaction(SIGFPE,  &sa, nullptr);
    sigaction(SIGABRT, &sa, nullptr);
#elif defined(_WIN32)
    // Windows: SetUnhandledExceptionFilter + MiniDumpWriteDump is the
    // correct path. Not yet implemented for commit 1 — Linux is the
    // active build target. Tracked in notes/.
    std::fprintf(stderr,
                 "CrashHandler: Windows crash capture not yet implemented\n");
#endif
}

}  // namespace h5reader::diagnostics
