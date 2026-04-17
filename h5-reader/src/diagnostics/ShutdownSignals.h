// ShutdownSignals — bridge terminal / OS shutdown signals to a clean
// Qt quit, cross-platform.
//
// Qt's GUI apps do not install terminal-signal handlers by default; on
// POSIX, Ctrl-C kills the process before aboutToQuit fires, which
// means the VTK-finalise-before-GL-context-destruction sequence
// (ReaderMainWindow::shutdown) never runs. On Windows the same thing
// happens when the host console issues CTRL_C_EVENT or CTRL_CLOSE_EVENT.
//
// Platform strategies:
//   POSIX  — self-pipe + QSocketNotifier. Async-signal-safe handler
//            writes a byte; the notifier wakes on the GUI thread and
//            calls QCoreApplication::quit().
//   Win32  — SetConsoleCtrlHandler posts quit via
//            QMetaObject::invokeMethod(..., Qt::QueuedConnection), which
//            is documented thread-safe. The handler returns TRUE so
//            the OS waits for our quit path to complete.
//
// Call InstallShutdownSignalHandlers() once, after QApplication
// construction and before event-loop entry. Idempotent.

#pragma once

namespace h5reader::diagnostics {

// Handles SIGINT + SIGTERM on POSIX; Ctrl-C / Ctrl-Break / console
// close / log-off / system-shutdown on Windows. All paths cause
// QCoreApplication::quit() to be posted on the event loop.
void InstallShutdownSignalHandlers();

}  // namespace h5reader::diagnostics
