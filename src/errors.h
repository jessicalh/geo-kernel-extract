#pragma once
//
// Error codes for the library. 0 = success, non-zero = an ErrorCode
// below, grouped into per-subsystem ranges; the value identifies the
// failure category. The caller checks != 0 and propagates.
//

#include <cstdint>

namespace nmr {

// Status: 0 = ok, else one of the ErrorCode values below.
using Status = std::int32_t;

// Error code integers, grouped by subsystem.
enum ErrorCode : std::int32_t {
    kOk = 0,

    kSessionTomlReadFailed       = 0x0001,
    kSessionAimnet2LoadFailed    = 0x0002,
    kSessionCalculatorConfigBad  = 0x0003,
    kSessionValidationFailed     = 0x0004,
    kSessionTripeptideDbLoadFailed = 0x0005,
    kSessionLarsenHBondGridLoadFailed = 0x0006,

    kTprReadFailed               = 0x0010,
    kProteinBuildFailed          = 0x0011,
    kChargesExtractFailed        = 0x0012,
    kBondedParamsExtractFailed   = 0x0013,

    kXtcOpenFailed               = 0x0020,
    kEdrLoadFailed               = 0x0021,
    kFrameReadFailed             = 0x0022,
    kPbcFixFailed                = 0x0023,
    kFrameExtractFailed          = 0x0024,
    kProteinAtomCountMismatch    = 0x0025,

    kCalculatorPipelineFailed    = 0x0030,

    kAttachRejectedSingleton     = 0x0040,
    kAttachDependencyUnmet       = 0x0041,
    kConfigRequiresAimnet2       = 0x0042,

    kH5WriteFailed               = 0x0050,
    kNpyWriteFailed              = 0x0051,
};

}  // namespace nmr
