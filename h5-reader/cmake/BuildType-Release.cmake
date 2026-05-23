# BuildType-Release.cmake — LTO + size/perf flags for Release builds.
#
# RWDI is the default dev preset on every platform and gets CMake's
# native -O2 -g -DNDEBUG flags without a module. Release is for
# distribution binaries: full optimization, link-time optimization
# when supported, stripped symbols (CMake default for Release).

include_guard(GLOBAL)

include(CheckIPOSupported)

function(h5reader_apply_build_type_settings target)
    check_ipo_supported(RESULT _ipo_supported OUTPUT _ipo_error LANGUAGES CXX)
    if(_ipo_supported)
        set_target_properties(${target} PROPERTIES
            INTERPROCEDURAL_OPTIMIZATION TRUE)
        message(STATUS "h5reader: Release build; LTO/IPO enabled")
    else()
        message(STATUS "h5reader: Release build; LTO/IPO unavailable (${_ipo_error})")
    endif()
endfunction()
