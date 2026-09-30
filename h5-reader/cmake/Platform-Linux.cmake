# Platform-Linux.cmake — Linux-specific build settings for h5reader.
#
# Dependency stack:
#   Qt 6.10 OSS (validated with 6.10.2; includes FFmpeg multimedia plugin)
#   VTK 9.5 (built from source against the selected Qt prefix)
#   HDF5 1.10 (apt: libhdf5-dev)
#   Eigen 3.4 (apt: libeigen3-dev)
#
# Set H5READER_QT_DIR and H5READER_VTK_DIR for the selected toolchain.
#
# Exposes a single function consumed by CMakeLists.txt after the
# target is created:
#   h5reader_apply_platform_target_settings(<target>)

include_guard(GLOBAL)

function(h5reader_apply_platform_target_settings target)
    # HDF5 include-order fix. HighFive picks H5Dvlen_reclaim (<1.12)
    # vs H5Treclaim (>=1.12) from the HDF5 header it sees first.
    # Transitive includes can surface a different HDF5 version. Prefer the
    # headers selected by find_package, including explicit custom prefixes.
    if(HDF5_INCLUDE_DIRS)
        target_include_directories(${target}
            BEFORE PRIVATE ${HDF5_INCLUDE_DIRS})
    endif()

    target_compile_options(${target} PRIVATE
        -Wall -Wextra -Wpedantic -Wno-unused-parameter)

    # Crash-handler runtime needs dlsym / dladdr.
    target_link_libraries(${target} PRIVATE ${CMAKE_DL_LIBS})
endfunction()
