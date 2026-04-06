#pragma once
//
// OperationLog: structured UDP logging for every meaningful operation.
//
// This is not optional debugging. This is the system's external memory.
// When an agent builds layers of complexity and can no longer hold the
// full state in context, the operation log shows what actually happened.
//
// Every operation that changes state, computes a result, loads data,
// calls an external tool, or encounters an error emits a structured
// JSON message over UDP. A listener captures the stream.
//
// The log is ALWAYS ON. The cost of a UDP send is negligible compared
// to any calculation this library performs. The cost of NOT logging
// is hours of debugging silent failures.
//
// Message format (JSON over UDP or stderr):
// {
//   "ts": "2026-04-01T18:34:56.789",
//   "level": "info",
//   "op": "DsspResult::Compute",
//   "detail": "protein=1UBQ, residues=76, ss_assigned=76"
// }
//

#include <string>
#include <chrono>
#include <cstdint>

namespace nmr {

// Log channel bitmask. Set in TOML [logging].channels.
// Warnings and errors ALWAYS log regardless of channel mask.
// Info-level logs only fire if their channel bit is set.
enum LogChannel : uint32_t {
    LogDiag0          = 1 << 0,   // Reserved: diagnostic slot 0
    LogDiag1          = 1 << 1,   // Reserved: diagnostic slot 1
    LogBondClassify   = 1 << 2,   // Bond classification (per-bond, chatty)
    LogNaming         = 1 << 3,   // NamingRegistry translations
    LogRingDetect     = 1 << 4,   // Ring detection and geometry
    LogCalcDispatch   = 1 << 5,   // Calculator dispatch
    LogCalcBiotSavart = 1 << 6,   // BiotSavart detail
    LogCalcHaighMal   = 1 << 7,   // HaighMallion detail
    LogCalcMcConnell  = 1 << 8,   // McConnell detail
    LogCalcOther      = 1 << 9,   // PQ, RSA, LD, HB, Coulomb
    LogFileIO         = 1 << 10,  // PDB/XYZ/XTC file reading
    LogDatabase       = 1 << 11,  // Database queries
    LogAtomMapping    = 1 << 12,  // Atom correspondence / WT-ALA matching
    LogDSSP           = 1 << 13,  // DSSP computation
    LogProtonation    = 1 << 14,  // Protonation detection
    LogConformation   = 1 << 15,  // Conformation add/modify
    LogResultAttach   = 1 << 16,  // ConformationResult attachment
    LogAPBS           = 1 << 17,  // APBS electrostatics
    LogToolExec       = 1 << 18,  // External tool execution
    LogMopac          = 1 << 19,  // MOPAC computation
    LogCharges        = 1 << 20,  // Charge assignment
    LogViewer         = 1 << 21,  // Viewer events
    LogAll            = 0xFFFFFFFF
};

class OperationLog {
public:
    enum class Level { Debug, Info, Warning, Error };

    // Configure UDP destination. Call once at startup.
    // If not called, logs go to stderr.
    static void ConfigureUdp(const std::string& host, int port);

    // Set the channel mask. Only Info on enabled channels is emitted.
    // Warn and Error ALWAYS emit regardless of mask.
    static void SetChannelMask(uint32_t mask);
    static uint32_t GetChannelMask();

    // Load channel mask from TOML [logging].channels.
    static void LoadChannelConfig(const std::string& tomlPath = "");

    // Log session start: records active channels.
    static void LogSessionStart();

    // Log with channel (Info checks mask, Warn/Error always emit).
    static void Log(Level level, uint32_t channel,
                    const std::string& operation,
                    const std::string& detail);

    // Log without channel (always emits — for Warn/Error).
    static void Log(Level level,
                    const std::string& operation,
                    const std::string& detail);

    // Convenience
    static void Info(const std::string& op, const std::string& detail) {
        Log(Level::Info, LogAll, op, detail);
    }
    static void Info(uint32_t ch, const std::string& op, const std::string& detail) {
        Log(Level::Info, ch, op, detail);
    }
    static void Warn(const std::string& op, const std::string& detail) {
        Log(Level::Warning, op, detail);
    }
    static void Error(const std::string& op, const std::string& detail) {
        Log(Level::Error, op, detail);
    }

    // Scoped logger: logs BEGIN on construction, END with elapsed time on destruction.
    class Scope {
    public:
        Scope(const std::string& operation, const std::string& detail = "");
        ~Scope();
    private:
        std::string operation_;
        std::chrono::steady_clock::time_point start_;
    };

    static bool IsUdpConfigured();

private:
    static void SendUdp(const std::string& json);
    static void SendStderr(const std::string& json);

    static std::string host_;
    static int port_;
    static int socket_;
    static bool udpConfigured_;
    static uint32_t channelMask_;
};

}  // namespace nmr
