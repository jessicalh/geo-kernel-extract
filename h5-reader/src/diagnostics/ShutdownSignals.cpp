#include "ShutdownSignals.h"

#include <QApplication>
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

void CloseApplicationWindows() {
    QApplication::closeAllWindows();
}

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
            << "received signal" << static_cast<int>(byte)
            << "- requesting graceful close";
        CloseApplicationWindows();
    });

    struct sigaction sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sa_handler = &PosixHandler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;   // don't interrupt I/O; don't reset on use

    sigaction(SIGINT,  &sa, nullptr);
    sigaction(SIGTERM, &sa, nullptr);

    qCInfo(cSignals).noquote()
        << "POSIX: SIGINT/SIGTERM graceful-close bridge installed";
}

#else
// ---------------------------------------------------------------------------
// Windows — SetConsoleCtrlHandler.
//
// Ctrl-C and Ctrl-Break arrive on a thread created by Windows.
// QMetaObject::invokeMethod with Qt::QueuedConnection posts the close
// request onto the GUI thread. Console-window close, logoff, and system
// shutdown do not provide a reliable lifetime for queued GUI work, so this
// bridge does not claim to handle them.
// ---------------------------------------------------------------------------
BOOL WINAPI Win32Handler(DWORD dwCtrlType) {
    if (dwCtrlType != CTRL_C_EVENT && dwCtrlType != CTRL_BREAK_EVENT)
        return FALSE;

    if (auto* app = QCoreApplication::instance()) {
        QMetaObject::invokeMethod(app, [dwCtrlType]() {
            qCInfo(cSignals).noquote()
                << "Windows console event"
                << (dwCtrlType == CTRL_C_EVENT ? "CTRL_C" : "CTRL_BREAK")
                << "- requesting graceful close";
            CloseApplicationWindows();
        }, Qt::QueuedConnection);
        return TRUE;
    }
    return FALSE;
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
        << "Windows: SetConsoleCtrlHandler graceful-close bridge installed";
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
