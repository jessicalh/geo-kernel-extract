# Platform-Linux.cmake — Linux-specific build settings for h5reader.
#
# Pinned dependency stack (per notes/BUILD_LAYOUT_PLAN_2026-05-23.md):
#   Qt 6.4 (apt: qt6-base-dev qt6-charts-dev)
#   VTK 9.5 (built from source, typically at $ENV{HOME}/VTK/)
#   HDF5 1.10 (apt: libhdf5-dev)
#   Eigen 3.4 (apt: libeigen3-dev)
#
# Same presets apply to DGX Spark (Linux aarch64) — Ubuntu 24 LTS with
# the same apt package names.
#
# Exposes a single function consumed by CMakeLists.txt after the
# target is created:
#   h5reader_apply_platform_target_settings(<target>)

include_guard(GLOBAL)

function(h5reader_apply_platform_target_settings target)
    # HDF5 include-order fix. HighFive picks H5Dvlen_reclaim (<1.12)
    # vs H5Treclaim (>=1.12) from the HDF5 header it sees first.
    # System lib is 1.10; transitive include paths from micromamba /
    # OpenBabel / etc. can surface 1.14 headers. Force the system
    # 1.10 path first on this target only.
    if(EXISTS /usr/include/hdf5/serial)
        target_include_directories(${target}
            BEFORE PRIVATE /usr/include/hdf5/serial)
    endif()

    target_compile_options(${target} PRIVATE
        -Wall -Wextra -Wpedantic -Wno-unused-parameter)

    # Crash-handler runtime needs dlsym / dladdr.
    target_link_libraries(${target} PRIVATE ${CMAKE_DL_LIBS})
endfunction()
