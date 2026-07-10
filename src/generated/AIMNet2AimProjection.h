#pragma once
//
// Committed AIMNet2 aim-projection basis data.
//
// This header exposes only the immutable table contract. The generated/
// translation unit contains the literal float32 coefficients and their
// generation provenance; all runtime projection arithmetic belongs to the
// owning AIMNet2Result calculator.

#include "ConformationAtom.h"

#include <cstddef>
#include <string_view>

namespace nmr {

inline constexpr std::size_t kAimnet2AimProjectionElementSlots = 5;
inline constexpr std::size_t kAimnet2AimProjectionDims =
    AIMNET2_AIM_PROJECTION_DIMS;
inline constexpr std::string_view kAimnet2AimProjectionBasisId =
    "splitmix64_0xA17E20260708_achlioptas_32x256_element_HCNOS";

static_assert(kAimnet2AimProjectionDims == 32,
              "the committed AIMNet2 projection basis has 32 rows");
static_assert(AIMNET2_AIM_DIMS == 256,
              "the committed AIMNet2 projection basis has 256 columns");

enum class Aimnet2AimProjectionElementSlot : std::size_t {
    H = 0,
    C = 1,
    N = 2,
    O = 3,
    S = 4,
};

extern const float kAimnet2AimProjectionBasis
    [kAimnet2AimProjectionElementSlots]
    [kAimnet2AimProjectionDims]
    [AIMNET2_AIM_DIMS];

}  // namespace nmr
