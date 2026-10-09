#pragma once

#include <cmath>

#if defined(H5READER_USE_BOOST_MATH)
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#endif

namespace h5reader::physics::detail {

#if defined(H5READER_USE_BOOST_MATH)
// A point arbitrarily close to the wire can round its modulus to one.
// Return a nonfinite value there so the field's existing isfinite check can
// reject it, instead of letting Boost's default overflow exception escape.
using EllipticPolicy = boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::ignore_error>,
    boost::math::policies::overflow_error<boost::math::policies::ignore_error>>;
#endif

// Both implementations take the modulus k, not the parameter m = k*k.
inline double CompleteEllipticFirstKind(double modulus) {
#if defined(H5READER_USE_BOOST_MATH)
    return boost::math::ellint_1(modulus, EllipticPolicy{});
#else
    return std::comp_ellint_1(modulus);
#endif
}

inline double CompleteEllipticSecondKind(double modulus) {
#if defined(H5READER_USE_BOOST_MATH)
    return boost::math::ellint_2(modulus, EllipticPolicy{});
#else
    return std::comp_ellint_2(modulus);
#endif
}

}  // namespace h5reader::physics::detail
