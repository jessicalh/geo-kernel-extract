// CrashHandler — capture unhandled signals / exceptions into a companion
// .txt file with a stack trace and an object-census snapshot.
//
// Linux today: sigaction() on SIGSEGV/SIGBUS/SIGABRT/SIGILL/SIGFPE +
// backtrace()/backtrace_symbols_fd(). macOS when built: same sigaction
// path. Windows when built: SetUnhandledExceptionFilter +
// MiniDumpWriteDump (Dbghelp). One API regardless of platform:
// CrashHandler::Install(). Platform specifics are hidden.
//
// The handler writes to {dumpDir}/crash_{pid}_{timestamp}.txt and also
// prints the path to stderr so the operator knows where to look.
// Default dumpDir is the user's writable generic-data location under
// "h5reader/crashes"; override with SetDumpDirectory() before Install().

#pragma once

#include <QString>

namespace h5reader::diagnostics {

class CrashHandler {
public:
    // Install once, as early in main() as possible — ideally before
    // QApplication construction. Subsequent calls are no-ops.
    static void Install();

    // Override where dump files are written. Takes effect only if called
    // before Install(). Directory is created if it does not exist.
    static void SetDumpDirectory(const QString& path);

private:
    CrashHandler() = delete;
};

}  // namespace h5reader::diagnostics
