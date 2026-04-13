#include "OperationLog.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
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
    std::string escaped;
    for (char c : detail) {
        if (c == '"') escaped += "\\\"";
        else if (c == '\\') escaped += "\\\\";
        else if (c == '\n') escaped += "\\n";
        else escaped += c;
    }

    return "{\"ts\":\"" + CurrentTimestamp() +
           "\",\"level\":\"" + LevelString(level) +
           "\",\"op\":\"" + operation +
           "\",\"detail\":\"" + escaped + "\"}";
}


void OperationLog::Log(Level level, const std::string& operation,
                        const std::string& detail) {
    std::string json = BuildJson(level, operation, detail);
    if (fileConfigured_)
        SendFile(json);
    if (udpConfigured_)
        SendUdp(json);
    // Warn/Error always go to stderr — cannot be silently lost if listener is down.
    // Info/Debug go to stderr only as fallback when UDP is not configured.
    if (!udpConfigured_ || level == Level::Warning || level == Level::Error)
        SendStderr(json);
}

void OperationLog::Log(Level level, uint32_t channel,
                        const std::string& operation,
                        const std::string& detail) {
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
        const char* home = std::getenv("HOME");
        if (home) path = std::string(home) + "/.nmr_tools.toml";
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
            try {
                channelMask_ = static_cast<uint32_t>(std::stoul(val, nullptr, 0));
            } catch (...) {}
        } else if (key == "udp_host") {
            udp_host = val;
        } else if (key == "udp_port") {
            try { udp_port = std::stoi(val); } catch (...) {}
        } else if (key == "file") {
            log_file = val;
        }
    }

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
    OperationLog::Info(operation_ + " [END]",
                        "elapsed=" + std::to_string(ms) + "ms");
}

}  // namespace nmr
