#include "OperationLog.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <nlohmann/json.hpp>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

namespace nmr {

// Static members
std::string OperationLog::host_;
int OperationLog::port_ = 0;
int OperationLog::socket_ = -1;
bool OperationLog::udpConfigured_ = false;
std::ofstream OperationLog::fileStream_;
bool OperationLog::fileConfigured_ = false;
uint32_t OperationLog::channelMask_ = LogAll;


void OperationLog::ConfigureUdp(const std::string& host, int port) {
    host_ = host;
    port_ = port;
    socket_ = ::socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);

    if (socket_ >= 0) {
        udpConfigured_ = true;
        Log(Level::Info, "OperationLog::ConfigureUdp",
            "UDP logging to " + host + ":" + std::to_string(port));
    } else {
        Log(Level::Error, "OperationLog::ConfigureUdp",
            "Failed to create UDP socket, falling back to stderr");
    }
}


bool OperationLog::IsUdpConfigured() { return udpConfigured_; }
bool OperationLog::IsFileConfigured() { return fileConfigured_; }


void OperationLog::ConfigureFile(const std::string& path) {
    if (fileConfigured_) CloseFile();
    fileStream_.open(path, std::ios::out | std::ios::trunc);
    if (fileStream_.is_open()) {
        fileConfigured_ = true;
        Log(Level::Info, "OperationLog::ConfigureFile",
            "file logging to " + path);
    } else {
        Log(Level::Error, "OperationLog::ConfigureFile",
            "failed to open " + path);
    }
}


void OperationLog::CloseFile() {
    if (fileConfigured_) {
        fileStream_.close();
        fileConfigured_ = false;
    }
}


static std::string CurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;

    std::tm tm_buf;
    localtime_r(&time, &tm_buf);

    std::ostringstream oss;
    oss << std::put_time(&tm_buf, "%Y-%m-%dT%H:%M:%S")
        << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
}


static const char* LevelString(OperationLog::Level level) {
    switch (level) {
        case OperationLog::Level::Debug:   return "debug";
        case OperationLog::Level::Info:    return "info";
        case OperationLog::Level::Warning: return "warn";
        case OperationLog::Level::Error:   return "ERROR";
    }
    return "unknown";
}


static std::string BuildJson(OperationLog::Level level,
                              const std::string& operation,
                              const std::string& detail) {
    // ordered_json keeps ts/level/op/detail order; compact dump() matches
    // the historic line format and escapes all control chars correctly
    // (the old hand loop missed everything below \n).
    nlohmann::ordered_json j;
    j["ts"]     = CurrentTimestamp();
    j["level"]  = LevelString(level);
    j["op"]     = operation;
    j["detail"] = detail;
    // error_handler_t::replace: logging must never throw. detail/operation
    // carry non-guaranteed-UTF-8 bytes (file paths, PQerrorMessage, e.what())
    // and the default strict dump() throws type_error.316 on a bad byte —
    // an exception out of a logging call, often on an error path. replace
    // substitutes U+FFFD instead. The old hand escaper copied raw bytes and
    // could not throw; this preserves that infallibility.
    return j.dump(-1, ' ', false, nlohmann::ordered_json::error_handler_t::replace);
}


// Last resort when the normal emit path throws (e.g. bad_alloc building the
// line): emit a schema-shaped marker through one ::write, so the stream
// records that an entry was lost rather than the program terminating or
// losing the event silently. This path uses no C++ heap (no std::string /
// stream / container) and throws no C++ exception out of this noexcept
// fallback: chrono::now + localtime_r + strftime/snprintf all write into
// fixed stack buffers (the normal path allocates only because it uses
// ostringstream). libc time/locale formatting may touch its own internal
// state, but it will not raise a C++ exception here. The line carries a real
// "ts" so it matches the BuildJson schema downstream consumers parse.
static void EmitDroppedMarker() noexcept {
    char ts[32] = "unknown";
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  now.time_since_epoch()).count() % 1000;
    std::tm tm_buf;
    if (localtime_r(&t, &tm_buf)) {
        char stamp[20];
        if (std::strftime(stamp, sizeof(stamp), "%Y-%m-%dT%H:%M:%S", &tm_buf) > 0)
            std::snprintf(ts, sizeof(ts), "%s.%03d", stamp, static_cast<int>(ms));
    }
    char line[200];
    int n = std::snprintf(line, sizeof(line),
        "{\"ts\":\"%s\",\"level\":\"ERROR\",\"op\":\"OperationLog\","
        "\"detail\":\"log entry dropped: exception during emit\"}\n", ts);
    // n<0 can't happen with this fixed format, but clamp defensively; consume
    // write()'s warn_unused_result (nothing safe to do if even this fails).
    std::size_t len = (n > 0 && static_cast<std::size_t>(n) < sizeof(line))
                          ? static_cast<std::size_t>(n) : sizeof(line) - 1;
    if (::write(STDERR_FILENO, line, len) < 0) { /* give up */ }
}

void OperationLog::Log(Level level, const std::string& operation,
                        const std::string& detail) noexcept {
    // Unthrowable by contract (see header). The only realistic throw source
    // is BuildJson allocating (ostringstream / nlohmann); the sinks are
    // syscall / default-stream writes that do not throw. Contain everything:
    // a logging call must never unwind out of a destructor or an in-flight
    // error path.
    try {
        std::string json = BuildJson(level, operation, detail);
        if (fileConfigured_)
            SendFile(json);
        if (udpConfigured_)
            SendUdp(json);
        // Warn/Error always go to stderr — cannot be silently lost if listener is down.
        // Info/Debug go to stderr only as fallback when UDP is not configured.
        if (!udpConfigured_ || level == Level::Warning || level == Level::Error)
            SendStderr(json);
    } catch (...) {
        EmitDroppedMarker();
    }
}

void OperationLog::LogScopeEnd(const std::string& operation,
                               long long elapsed_ms) noexcept {
    // The "unthrowable log entry form" for Scope::~Scope: the destructor is
    // implicitly noexcept, and building the op/detail strings at the call site
    // would let a bad_alloc escape the destructor → std::terminate. So the
    // destructor passes only its already-allocated operation_ (by ref) and a
    // scalar; ALL string assembly happens here, inside the contained boundary.
    try {
        // Channelled (LogAll), matching the BEGIN line emitted by Scope's ctor
        // via Info(): both honour channelMask_, so BEGIN/END stay paired —
        // when Info is masked off, neither emits.
        Log(Level::Info, LogAll, operation + " [END]",
            "elapsed=" + std::to_string(elapsed_ms) + "ms");
    } catch (...) {
        EmitDroppedMarker();
    }
}

void OperationLog::Log(Level level, uint32_t channel,
                        const std::string& operation,
                        const std::string& detail) noexcept {
    // Warnings and errors ALWAYS emit regardless of channel mask.
    if (level == Level::Info || level == Level::Debug) {
        if ((channelMask_ & channel) == 0) return;
    }
    Log(level, operation, detail);
}

void OperationLog::SetChannelMask(uint32_t mask) { channelMask_ = mask; }
uint32_t OperationLog::GetChannelMask() { return channelMask_; }

void OperationLog::LoadChannelConfig(const std::string& tomlPath) {
    std::string path = tomlPath;
    if (path.empty()) {
        const char* env_path = std::getenv("NMR_TOOLS_TOML");
        if (env_path && *env_path) {
            path = env_path;
        } else {
            const char* home = std::getenv("HOME");
            if (home) path = std::string(home) + "/.nmr_tools.toml";
        }
    }
    if (path.empty()) return;

    std::ifstream in(path);
    if (!in.is_open()) return;

    std::string udp_host;
    int udp_port = 0;
    std::string log_file;

    std::string line;
    bool inLogging = false;
    while (std::getline(in, line)) {
        if (line.find("[logging]") != std::string::npos) { inLogging = true; continue; }
        if (line.find('[') != std::string::npos && line.find("[logging]") == std::string::npos)
            inLogging = false;
        if (!inLogging) continue;

        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        while (!key.empty() && (key.back() == ' ' || key.back() == '\t')) key.pop_back();
        while (!val.empty() && (val.front() == ' ' || val.front() == '\t')) val = val.substr(1);
        while (!val.empty() && (val.back() == ' ' || val.back() == '\t' || val.back() == '"')) val.pop_back();
        while (!val.empty() && val.front() == '"') val = val.substr(1);

        if (key == "channels") {
            // Fail loud on every way the value can be wrong: unparseable,
            // trailing junk after the number, or wider than the 32-bit mask
            // (stoul returns unsigned long; the cast would silently truncate).
            try {
                std::size_t pos = 0;
                unsigned long parsed = std::stoul(val, &pos, 0);
                if (pos != val.size()) {
                    Warn("OperationLog::LoadChannelConfig",
                         "trailing junk after [logging].channels value, keeping default: \"" + val + "\"");
                } else if (parsed > 0xFFFFFFFFUL) {
                    Warn("OperationLog::LoadChannelConfig",
                         "[logging].channels exceeds the 32-bit mask, keeping default: \"" + val + "\"");
                } else {
                    channelMask_ = static_cast<uint32_t>(parsed);
                }
            } catch (...) {
                Warn("OperationLog::LoadChannelConfig",
                     "unparseable [logging].channels value, keeping default: \"" + val + "\"");
            }
        } else if (key == "udp_host") {
            udp_host = val;
        } else if (key == "udp_port") {
            // Same discipline: unparseable, trailing junk, or out of the valid
            // 1-65535 port range (htons would otherwise silently truncate).
            try {
                std::size_t pos = 0;
                int parsed = std::stoi(val, &pos);
                if (pos != val.size()) {
                    Warn("OperationLog::LoadChannelConfig",
                         "trailing junk after [logging].udp_port value, UDP disabled: \"" + val + "\"");
                } else if (parsed < 1 || parsed > 65535) {
                    Warn("OperationLog::LoadChannelConfig",
                         "[logging].udp_port out of range 1-65535, UDP disabled: \"" + val + "\"");
                } else {
                    udp_port = parsed;
                }
            } catch (...) {
                Warn("OperationLog::LoadChannelConfig",
                     "unparseable [logging].udp_port value, UDP disabled: \"" + val + "\"");
            }
        } else if (key == "file") {
            log_file = val;
        }
    }

    if (const char* env_host = std::getenv("NMR_LOG_UDP_HOST"); env_host && *env_host)
        udp_host = env_host;
    if (const char* env_port = std::getenv("NMR_LOG_UDP_PORT"); env_port && *env_port) {
        try {
            int parsed = std::stoi(env_port);
            if (parsed >= 1 && parsed <= 65535) udp_port = parsed;
        } catch (...) {
            Warn("OperationLog::LoadChannelConfig",
                 "unparseable NMR_LOG_UDP_PORT value, UDP TOML/default retained");
        }
    }
    if (const char* env_file = std::getenv("NMR_LOG_FILE"); env_file && *env_file)
        log_file = env_file;

    // Configure UDP if both host and port are specified.
    if (!udp_host.empty() && udp_port > 0)
        ConfigureUdp(udp_host, udp_port);

    // Configure file logging if path is specified.
    if (!log_file.empty())
        ConfigureFile(log_file);
}

void OperationLog::LogSessionStart() {
    std::string active;
    auto bit = [&](uint32_t b, const char* name) {
        if (channelMask_ & b) {
            if (!active.empty()) active += "|";
            active += name;
        }
    };
    bit(LogDiag0, "DIAG0"); bit(LogDiag1, "DIAG1");
    bit(LogBondClassify, "bond"); bit(LogNaming, "naming");
    bit(LogRingDetect, "ring"); bit(LogCalcDispatch, "dispatch");
    bit(LogCalcBiotSavart, "BS"); bit(LogCalcHaighMal, "HM");
    bit(LogCalcMcConnell, "MC"); bit(LogCalcOther, "calc");
    bit(LogFileIO, "file"); bit(LogDatabase, "db");
    bit(LogAtomMapping, "mapping"); bit(LogDSSP, "dssp");
    bit(LogProtonation, "prot"); bit(LogConformation, "conf");
    bit(LogResultAttach, "attach"); bit(LogAPBS, "apbs");
    bit(LogToolExec, "tool"); bit(LogMopac, "mopac");
    bit(LogCharges, "charges"); bit(LogViewer, "viewer");

    std::ostringstream hex;
    hex << std::hex << std::uppercase << channelMask_;
    Log(Level::Info, "OperationLog::SessionStart",
        "channels=0x" + hex.str() + " (" + active + ")");
}


void OperationLog::SendUdp(const std::string& json) {
    struct sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(static_cast<uint16_t>(port_));
    inet_pton(AF_INET, host_.c_str(), &addr.sin_addr);
    sendto(socket_, json.c_str(), json.size(), 0,
           reinterpret_cast<struct sockaddr*>(&addr), sizeof(addr));
}


void OperationLog::SendFile(const std::string& json) {
    fileStream_ << json << "\n";
    // No flush per message — hot-path writes must not block.
    // OS flushes on buffer fill and on CloseFile().
}


void OperationLog::SendStderr(const std::string& json) {
    std::cerr << json << "\n";
}


// Scoped logger
OperationLog::Scope::Scope(const std::string& operation, const std::string& detail)
    : operation_(operation)
    , start_(std::chrono::steady_clock::now()) {
    OperationLog::Info(operation + " [BEGIN]", detail);
}

OperationLog::Scope::~Scope() {
    auto elapsed = std::chrono::steady_clock::now() - start_;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    // Pass operation_ by ref + the scalar ms only — no string assembly in this
    // (implicitly noexcept) destructor body. LogScopeEnd builds and emits the
    // line inside its own contained boundary so nothing can escape here.
    OperationLog::LogScopeEnd(operation_, static_cast<long long>(ms));
}

}  // namespace nmr
