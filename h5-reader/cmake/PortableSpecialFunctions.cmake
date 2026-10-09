include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

# libc++ does not provide all of the C++17 mathematical special functions.
# Detect the actual library capability so other platforms keep their existing
# implementation without acquiring another dependency.
check_cxx_source_compiles([=[
    #include <cmath>
    int main(int argc, char**) {
        const double k = argc > 1 ? 0.5 : 0.75;
        return std::isfinite(std::comp_ellint_1(k) + std::comp_ellint_2(k)) ? 0 : 1;
    }
]=] H5READER_HAS_STD_COMPLETE_ELLIPTIC_INTEGRALS)

add_library(h5reader_special_functions INTERFACE)
target_compile_features(h5reader_special_functions INTERFACE cxx_std_17)

if(H5READER_HAS_STD_COMPLETE_ELLIPTIC_INTEGRALS)
    message(STATUS "Complete elliptic integrals: C++ standard library")
else()
    find_path(H5READER_BOOST_MATH_INCLUDE_DIR
        NAMES boost/math/special_functions/ellint_1.hpp
        DOC "Boost.Math include directory (standalone headers are sufficient)")
    if(NOT H5READER_BOOST_MATH_INCLUDE_DIR)
        message(FATAL_ERROR
            "This C++ standard library lacks complete elliptic integrals. "
            "Provide Boost.Math standalone headers via H5READER_BOOST_MATH_INCLUDE_DIR.")
    endif()

    cmake_push_check_state(RESET)
    set(CMAKE_REQUIRED_INCLUDES "${H5READER_BOOST_MATH_INCLUDE_DIR}")
    set(CMAKE_REQUIRED_DEFINITIONS -DBOOST_MATH_STANDALONE)
    check_cxx_source_compiles([=[
        #include <boost/math/special_functions/ellint_1.hpp>
        #include <boost/math/special_functions/ellint_2.hpp>
        #include <boost/math/policies/policy.hpp>
        int main(int argc, char**) {
            using namespace boost::math::policies;
            using Policy = policy<domain_error<ignore_error>, overflow_error<ignore_error>>;
            const double k = argc > 1 ? 0.5 : 0.75;
            return boost::math::ellint_1(k, Policy{}) > boost::math::ellint_2(k, Policy{}) ? 0 : 1;
        }
    ]=] H5READER_HAS_USABLE_BOOST_MATH)
    cmake_pop_check_state()
    if(NOT H5READER_HAS_USABLE_BOOST_MATH)
        message(FATAL_ERROR
            "Boost.Math at H5READER_BOOST_MATH_INCLUDE_DIR does not support the "
            "required standalone elliptic integrals. See the CMake configure log.")
    endif()

    target_include_directories(h5reader_special_functions SYSTEM INTERFACE
        "${H5READER_BOOST_MATH_INCLUDE_DIR}")
    target_compile_definitions(h5reader_special_functions INTERFACE
        H5READER_USE_BOOST_MATH=1 BOOST_MATH_STANDALONE)
    message(STATUS "Complete elliptic integrals: Boost.Math standalone headers")
endif()
