#include "physics/SpecialFunctions.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

int main() {
    struct Reference {
        double modulus;
        double first;
        double second;
    };
    // Independent references generated with mpmath 1.3.0 at 100 decimal
    // digits: ellipk(mpf(k)**2), ellipe(mpf(k)**2), where k is the exact
    // binary64 input. mpmath takes the parameter; our API takes the modulus.
    const std::array<Reference, 11> references{{
        {0.0, 1.5707963267948966, 1.5707963267948966},
        {1e-10, 1.5707963267948966, 1.5707963267948966},
        {0.01, 1.5708355989121523, 1.5707570561503852},
        {0.1, 1.5747455615173560, 1.5668619420216683},
        {0.5, 1.6857503548125961, 1.4674622093394272},
        {0.70710678118654757, 1.8540746773013719, 1.3506438810476755},
        {0.9, 2.2805491384227703, 1.1716970527816142},
        {0.99, 3.3566005233611920, 1.0284758090288040},
        {0.9999, 5.6451482168297478, 1.0005145000837812},
        {0.999999999999, 14.855242389793775, 1.0000000000143550},
        {0.99999999999999989, 19.408121055678471, 1.0000000000000020},
    }};
    int failures = 0;
    const auto check = [&](double actual, double expected, const char* kind, double k) {
        // Standard libraries may first round 1-k*k near the singularity;
        // allow that conditioning loss while retaining a stringent reference.
        if (!std::isfinite(actual) || std::abs(actual - expected) > 2e-13 * std::max(1.0, std::abs(expected))) {
            std::cerr << kind << " mismatch at modulus " << k << ": " << actual << " vs " << expected << '\n';
            ++failures;
        }
    };
    for (const auto& reference : references) {
        check(h5reader::physics::detail::CompleteEllipticFirstKind(reference.modulus),
              reference.first, "K", reference.modulus);
        check(h5reader::physics::detail::CompleteEllipticSecondKind(reference.modulus),
              reference.second, "E", reference.modulus);
    }
#if defined(H5READER_USE_BOOST_MATH)
    // The fallback must allow the caller to reject a rounded singularity
    // without an exception escaping the circular-field calculation.
    if (!std::isinf(h5reader::physics::detail::CompleteEllipticFirstKind(1.0))) {
        std::cerr << "K(1) must report the singularity as infinity\n";
        ++failures;
    }
    check(h5reader::physics::detail::CompleteEllipticSecondKind(1.0), 1.0, "E", 1.0);
#endif
    if (failures != 0)
        return 1;
    std::cout << "Complete elliptic integrals agree with independent high-precision references\n";
    return 0;
}
