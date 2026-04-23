#pragma once
//
// errors.h — shared error codes for the library.
//
// Return convention: 0 = success, non-zero = error. Each numeric range
// is owned by a subsystem so a code uniquely identifies the site. The
// caller's job is to check != 0 and propagate; the emitting site
// writes a diagnostic to OperationLog::Error at the moment of failure.
// Status-as-int lets helpers compose without allocating, and lets us
// chain ops with a simple non-zero-propagates macro.
//
// No exception throwing in library code (PATTERNS.md §9). Callers of
// the trajectory framework are expected to check statuses and stop.
//
// If we later need richer status (bitmask of warnings, severity bits,
// category bits), the top 16 bits of a 32-bit status are reserved for
// that use; the bottom 16 bits always carry the ErrorCode integer.
// For now only ErrorCode values are returned.
//

#include <cstdint>

namespace nmr {

// Status is a plain int. 0 = ok. Non-zero = one of the ErrorCode values
// below (bottom 16 bits) plus any future status bits (top 16 bits).
using Status = std::int32_t;

// Error code integers. Grouped by subsystem so new codes slot in
// without colliding. Reserve 0x0000 - 0x00FF for the framework
// (Session, Trajectory, TrajectoryProtein, handler). Reserve
// 0x0100+ for subsystems and later categories.
enum ErrorCode : std::int32_t {
    kOk = 0,

    // Framework: Session / configuration / startup
    kSessionTomlReadFailed       = 0x0001,
    kSessionAimnet2LoadFailed    = 0x0002,
    kSessionCalculatorConfigBad  = 0x0003,
    kSessionValidationFailed     = 0x0004,

    // Framework: TrajectoryProtein construction
    kTprReadFailed               = 0x0010,
    kProteinBuildFailed          = 0x0011,
    kChargesExtractFailed        = 0x0012,
    kBondedParamsExtractFailed   = 0x0013,

    // Framework: Trajectory / handler
    kXtcOpenFailed               = 0x0020,
    kEdrLoadFailed               = 0x0021,
    kFrameReadFailed             = 0x0022,
    kPbcFixFailed                = 0x0023,
    kFrameExtractFailed          = 0x0024,
    kProteinAtomCountMismatch    = 0x0025,

    // Framework: per-frame calculator pipeline
    kCalculatorPipelineFailed    = 0x0030,

    // Framework: attach / dependency
    kAttachRejectedSingleton     = 0x0040,
    kAttachDependencyUnmet       = 0x0041,
    kConfigRequiresAimnet2       = 0x0042,

    // Framework: I/O
    kH5WriteFailed               = 0x0050,
    kNpyWriteFailed              = 0x0051,

    // Reserved for future expansion: 0x0060 - 0x00FF
};

}  // namespace nmr
