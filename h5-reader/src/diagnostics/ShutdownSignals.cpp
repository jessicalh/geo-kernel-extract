#include "ShutdownSignals.h"

#include <QCoreApplication>
#include <QLoggingCategory>

#include <atomic>
#include <cstring>

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#else
#  include <QSocketNotifier>
#  include <signal.h>
#  include <sys/socket.h>
#  include <sys/types.h>
#  include <unistd.h>
#endif

namespace h5reader::diagnostics {

namespace {
Q_LOGGING_CATEGORY(cSignals, "h5reader.signals")

std::atomic<bool> g_installed{false};

#ifndef _WIN32
// ---------------------------------------------------------------------------
// POSIX — self-pipe + QSocketNotifier.
// ---------------------------------------------------------------------------
int g_sigPipe[2] = {-1, -1};

// Async-signal-safe. write(2) is on the allowed-functions list; no
// other library calls happen here.
void PosixHandler(int sig) {
    const char byte = static_cast<char>(sig & 0xFF);
    const ssize_t rc = ::write(g_sigPipe[1], &byte, 1);
    (void)rc;
}

void InstallPosix() {
    if (::socketpair(AF_UNIX, SOCK_STREAM, 0, g_sigPipe) != 0) {
        qCCritical(cSignals).noquote()
            << "socketpair() failed — terminal signals will not trigger clean quit";
        return;
    }

    auto* notifier = new QSocketNotifier(g_sigPipe[0], QSocketNotifier::Read,
                                         QCoreApplication::instance());
    QObject::connect(notifier, &QSocketNotifier::activated,
                     QCoreApplication::instance(),
                     [](QSocketDescriptor /*sd*/, QSocketNotifier::Type /*t*/) {
        char byte = 0;
        const ssize_t rc = ::read(g_sigPipe[0], &byte, 1);
        (void)rc;
        qCInfo(cSignals).noquote()
            << "received signal" << static_cast<int>(byte) << "— quitting";
        if (auto* app = QCoreApplication::instance()) app->quit();
    });

    struct sigaction sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sa_handler = &PosixHandler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;   // don't interrupt I/O; don't reset on use

    sigaction(SIGINT,  &sa, nullptr);
    sigaction(SIGTERM, &sa, nullptr);

    qCInfo(cSignals).noquote()
        << "POSIX: SIGINT/SIGTERM → QCoreApplication::quit bridge installed";
}

#else
// ---------------------------------------------------------------------------
// Windows — SetConsoleCtrlHandler.
//
// The handler runs on a separate thread the OS creates when the
// console event fires. QMetaObject::invokeMethod with
// Qt::QueuedConnection is documented thread-safe and will post the
// quit() call onto the GUI thread's event loop. Returning TRUE from
// the handler tells Windows we've taken responsibility; the OS gives
// our process a few seconds to complete before hard-terminating.
// ---------------------------------------------------------------------------
BOOL WINAPI Win32Handler(DWORD dwCtrlType) {
    const char* name = "?";
    switch (dwCtrlType) {
        case CTRL_C_EVENT:        name = "CTRL_C";        break;
        case CTRL_BREAK_EVENT:    name = "CTRL_BREAK";    break;
        case CTRL_CLOSE_EVENT:    name = "CTRL_CLOSE";    break;
        case CTRL_LOGOFF_EVENT:   name = "CTRL_LOGOFF";   break;
        case CTRL_SHUTDOWN_EVENT: name = "CTRL_SHUTDOWN"; break;
        default: return FALSE;   // not ours; let the default handler run
    }

    // qCInfo is not guaranteed async-safe from a non-Qt thread the OS
    // spawned, but it uses Qt's logging machinery which is mutex-
    // protected internally. Safe in practice; if a future stress test
    // shows otherwise, we drop the log line and rely on the GUI-side
    // lifecycle logs after quit posts.
    qCInfo(cSignals).noquote()
        << "Windows console event" << name << "— posting quit";

    if (auto* app = QCoreApplication::instance()) {
        QMetaObject::invokeMethod(app, "quit", Qt::QueuedConnection);
    }

    // TRUE = handled. Windows waits a short time for the process to
    // exit before hard-terminating it, which is enough for our
    // aboutToQuit → shutdown() path.
    return TRUE;
}

void InstallWin32() {
    if (!SetConsoleCtrlHandler(&Win32Handler, TRUE)) {
        qCCritical(cSignals).noquote()
            << "SetConsoleCtrlHandler() failed — terminal signals will not "
               "trigger clean quit. err="
            << static_cast<unsigned long>(GetLastError());
        return;
    }
    qCInfo(cSignals).noquote()
        << "Windows: SetConsoleCtrlHandler → QCoreApplication::quit bridge installed";
}
#endif  // _WIN32

}  // namespace

void InstallShutdownSignalHandlers() {
    if (g_installed.exchange(true)) return;
#ifdef _WIN32
    InstallWin32();
#else
    InstallPosix();
#endif
}

}  // namespace h5reader::diagnostics
