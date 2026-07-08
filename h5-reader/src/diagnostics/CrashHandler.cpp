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

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <Windows.h>
#  include <DbgHelp.h>
#  include <fcntl.h>
#  include <io.h>
#  include <share.h>
#  include <sys/stat.h>
#  ifndef EXCEPTION_FATAL_APP_EXIT
#    define EXCEPTION_FATAL_APP_EXIT 0x40000015
#  endif
#endif

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

#ifdef _WIN32
std::atomic<bool> g_windowsCrashWritten{false};

bool Utf8ToWide(const char* src, wchar_t* dst, size_t dstCount) {
    if (!src || !dst || dstCount == 0)
        return false;
    const int rc = MultiByteToWideChar(CP_UTF8, 0, src, -1,
                                       dst, static_cast<int>(dstCount));
    if (rc <= 0) {
        dst[0] = L'\0';
        return false;
    }
    return true;
}

void WinWrite(int fd, const char* s) {
    if (!s)
        return;
    const int rc = ::_write(fd, s, static_cast<unsigned int>(std::strlen(s)));
    (void)rc;
}

void WinWriteInt(int fd, unsigned long v) {
    char buf[32];
    const int n = std::snprintf(buf, sizeof(buf), "%lu", v);
    if (n > 0) {
        const int rc = ::_write(fd, buf, static_cast<unsigned int>(n));
        (void)rc;
    }
}

const char* WindowsExceptionName(DWORD code) {
    switch (code) {
        case EXCEPTION_ACCESS_VIOLATION: return "EXCEPTION_ACCESS_VIOLATION";
        case EXCEPTION_ARRAY_BOUNDS_EXCEEDED: return "EXCEPTION_ARRAY_BOUNDS_EXCEEDED";
        case EXCEPTION_DATATYPE_MISALIGNMENT: return "EXCEPTION_DATATYPE_MISALIGNMENT";
        case EXCEPTION_FLT_DENORMAL_OPERAND: return "EXCEPTION_FLT_DENORMAL_OPERAND";
        case EXCEPTION_FLT_DIVIDE_BY_ZERO: return "EXCEPTION_FLT_DIVIDE_BY_ZERO";
        case EXCEPTION_FLT_INEXACT_RESULT: return "EXCEPTION_FLT_INEXACT_RESULT";
        case EXCEPTION_FLT_INVALID_OPERATION: return "EXCEPTION_FLT_INVALID_OPERATION";
        case EXCEPTION_FLT_OVERFLOW: return "EXCEPTION_FLT_OVERFLOW";
        case EXCEPTION_FLT_STACK_CHECK: return "EXCEPTION_FLT_STACK_CHECK";
        case EXCEPTION_FLT_UNDERFLOW: return "EXCEPTION_FLT_UNDERFLOW";
        case EXCEPTION_ILLEGAL_INSTRUCTION: return "EXCEPTION_ILLEGAL_INSTRUCTION";
        case EXCEPTION_IN_PAGE_ERROR: return "EXCEPTION_IN_PAGE_ERROR";
        case EXCEPTION_INT_DIVIDE_BY_ZERO: return "EXCEPTION_INT_DIVIDE_BY_ZERO";
        case EXCEPTION_INT_OVERFLOW: return "EXCEPTION_INT_OVERFLOW";
        case EXCEPTION_INVALID_DISPOSITION: return "EXCEPTION_INVALID_DISPOSITION";
        case EXCEPTION_NONCONTINUABLE_EXCEPTION: return "EXCEPTION_NONCONTINUABLE_EXCEPTION";
        case EXCEPTION_PRIV_INSTRUCTION: return "EXCEPTION_PRIV_INSTRUCTION";
        case EXCEPTION_STACK_OVERFLOW: return "EXCEPTION_STACK_OVERFLOW";
        case EXCEPTION_FATAL_APP_EXIT: return "EXCEPTION_FATAL_APP_EXIT";
        default: return "UNKNOWN";
    }
}

bool IsLikelyFatalWindowsException(DWORD code) {
    switch (code) {
        case EXCEPTION_ACCESS_VIOLATION:
        case EXCEPTION_ARRAY_BOUNDS_EXCEEDED:
        case EXCEPTION_DATATYPE_MISALIGNMENT:
        case EXCEPTION_FLT_DENORMAL_OPERAND:
        case EXCEPTION_FLT_DIVIDE_BY_ZERO:
        case EXCEPTION_FLT_INEXACT_RESULT:
        case EXCEPTION_FLT_INVALID_OPERATION:
        case EXCEPTION_FLT_OVERFLOW:
        case EXCEPTION_FLT_STACK_CHECK:
        case EXCEPTION_FLT_UNDERFLOW:
        case EXCEPTION_ILLEGAL_INSTRUCTION:
        case EXCEPTION_IN_PAGE_ERROR:
        case EXCEPTION_INT_DIVIDE_BY_ZERO:
        case EXCEPTION_INT_OVERFLOW:
        case EXCEPTION_INVALID_DISPOSITION:
        case EXCEPTION_NONCONTINUABLE_EXCEPTION:
        case EXCEPTION_PRIV_INSTRUCTION:
        case EXCEPTION_STACK_OVERFLOW:
        case EXCEPTION_FATAL_APP_EXIT:
            return true;
        default:
            return false;
    }
}

void BuildWindowsCrashPaths(char* txtPath, size_t txtCap,
                            char* dmpPath, size_t dmpCap) {
    SYSTEMTIME st;
    GetLocalTime(&st);
    const DWORD pid = GetCurrentProcessId();
    const char* dir = g_dumpDir[0] ? g_dumpDir : ".";
    std::snprintf(txtPath, txtCap,
                  "%s/crash_%lu_%04u%02u%02u_%02u%02u%02u.txt",
                  dir, static_cast<unsigned long>(pid),
                  static_cast<unsigned>(st.wYear),
                  static_cast<unsigned>(st.wMonth),
                  static_cast<unsigned>(st.wDay),
                  static_cast<unsigned>(st.wHour),
                  static_cast<unsigned>(st.wMinute),
                  static_cast<unsigned>(st.wSecond));
    std::snprintf(dmpPath, dmpCap,
                  "%s/crash_%lu_%04u%02u%02u_%02u%02u%02u.dmp",
                  dir, static_cast<unsigned long>(pid),
                  static_cast<unsigned>(st.wYear),
                  static_cast<unsigned>(st.wMonth),
                  static_cast<unsigned>(st.wDay),
                  static_cast<unsigned>(st.wHour),
                  static_cast<unsigned>(st.wMinute),
                  static_cast<unsigned>(st.wSecond));
}

bool WriteWindowsMiniDump(const wchar_t* dumpPath, EXCEPTION_POINTERS* ep) {
    HANDLE file = CreateFileW(dumpPath, GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS,
                              FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE)
        return false;

    MINIDUMP_EXCEPTION_INFORMATION exceptionInfo;
    exceptionInfo.ThreadId = GetCurrentThreadId();
    exceptionInfo.ExceptionPointers = ep;
    exceptionInfo.ClientPointers = FALSE;

    const MINIDUMP_TYPE dumpType = static_cast<MINIDUMP_TYPE>(
        MiniDumpWithFullMemory
        | MiniDumpWithHandleData
        | MiniDumpWithThreadInfo
        | MiniDumpWithProcessThreadData
        | MiniDumpWithFullMemoryInfo
        | MiniDumpWithTokenInformation);

    const BOOL ok = MiniDumpWriteDump(GetCurrentProcess(),
                                      GetCurrentProcessId(),
                                      file,
                                      dumpType,
                                      ep ? &exceptionInfo : nullptr,
                                      nullptr,
                                      nullptr);
    CloseHandle(file);
    return ok == TRUE;
}

void WriteWindowsCrashReport(EXCEPTION_POINTERS* ep, const char* origin) {
    if (g_windowsCrashWritten.exchange(true))
        return;

    char txtPath[4096];
    char dmpPath[4096];
    BuildWindowsCrashPaths(txtPath, sizeof(txtPath), dmpPath, sizeof(dmpPath));

    wchar_t txtWide[4096];
    wchar_t dmpWide[4096];
    const bool txtOk = Utf8ToWide(txtPath, txtWide, sizeof(txtWide) / sizeof(txtWide[0]));
    const bool dmpOk = Utf8ToWide(dmpPath, dmpWide, sizeof(dmpWide) / sizeof(dmpWide[0]));

    int fd = -1;
    if (txtOk) {
        const errno_t openErr = ::_wsopen_s(&fd, txtWide,
                                            _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY,
                                            _SH_DENYWR,
                                            _S_IREAD | _S_IWRITE);
        if (openErr != 0)
            fd = -1;
    }
    if (fd < 0)
        fd = 2;

    const DWORD code = (ep && ep->ExceptionRecord)
        ? ep->ExceptionRecord->ExceptionCode
        : 0;
    const void* address = (ep && ep->ExceptionRecord)
        ? ep->ExceptionRecord->ExceptionAddress
        : nullptr;

    WinWrite(fd, "h5reader crash: Windows exception\n");
    WinWrite(fd, "origin=");
    WinWrite(fd, origin ? origin : "unknown");
    WinWrite(fd, "\npid=");
    WinWriteInt(fd, static_cast<unsigned long>(GetCurrentProcessId()));
    WinWrite(fd, "\ntid=");
    WinWriteInt(fd, static_cast<unsigned long>(GetCurrentThreadId()));
    WinWrite(fd, "\nexception_code=0x");
    char codeBuf[32];
    const int cn = std::snprintf(codeBuf, sizeof(codeBuf), "%08lx",
                                 static_cast<unsigned long>(code));
    if (cn > 0)
        ::_write(fd, codeBuf, static_cast<unsigned int>(cn));
    WinWrite(fd, " (");
    WinWrite(fd, WindowsExceptionName(code));
    WinWrite(fd, ")\nexception_address=");
    char addrBuf[32];
    const int an = std::snprintf(addrBuf, sizeof(addrBuf), "%p", address);
    if (an > 0)
        ::_write(fd, addrBuf, static_cast<unsigned int>(an));

    bool dumpWritten = false;
    if (dmpOk)
        dumpWritten = WriteWindowsMiniDump(dmpWide, ep);
    WinWrite(fd, "\nminidump=");
    WinWrite(fd, dumpWritten ? dmpPath : "(failed)");
    WinWrite(fd, "\ntext_report=");
    WinWrite(fd, txtPath);
    WinWrite(fd, "\n\nObjectCensus snapshot (best effort):\n");
    ObjectCensus::Dump(fd);
    WinWrite(fd, "\nEnd of crash report.\n");

    if (fd != 2)
        ::_close(fd);

    std::fprintf(stderr,
                 "h5reader: crash report written to %s; minidump %s\n",
                 txtPath,
                 dumpWritten ? dmpPath : "failed");
    std::fflush(stderr);
}

LONG WINAPI VectoredCrashHandler(EXCEPTION_POINTERS* ep) {
    const DWORD code = (ep && ep->ExceptionRecord)
        ? ep->ExceptionRecord->ExceptionCode
        : 0;
    if (IsLikelyFatalWindowsException(code))
        WriteWindowsCrashReport(ep, "vectored exception");
    return EXCEPTION_CONTINUE_SEARCH;
}

LONG WINAPI UnhandledCrashHandler(EXCEPTION_POINTERS* ep) {
    WriteWindowsCrashReport(ep, "unhandled exception");
    return EXCEPTION_EXECUTE_HANDLER;
}
#endif  // _WIN32

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
    // Use both hooks: the vectored handler catches fatal first-chance
    // exceptions before CRT filters can interfere, while the unhandled
    // filter covers exceptions that reach process termination.
    AddVectoredExceptionHandler(1, &VectoredCrashHandler);
    SetUnhandledExceptionFilter(&UnhandledCrashHandler);
#endif
}

}  // namespace h5reader::diagnostics
