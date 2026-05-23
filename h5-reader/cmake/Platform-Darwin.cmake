# Platform-Darwin.cmake — macOS-specific build settings for h5reader.
#
# Target stack (per notes/BUILD_LAYOUT_PLAN_2026-05-23.md):
#   Qt Pro latest (installer at $ENV{HOME}/Qt/<ver>/macos)
#   VTK latest, built from source (typically $ENV{HOME}/VTK/)
#   HDF5 1.14 (brew: hdf5)
#   Eigen 3.4 (brew: eigen)
#
# Stub at migration step 1. Populated in step 4 when the Mac
# environment is first synced. Defines the same function surface
# as Platform-Linux so CMakeLists.txt is platform-agnostic.

include_guard(GLOBAL)

function(h5reader_apply_platform_target_settings target)
    target_compile_options(${target} PRIVATE
        -Wall -Wextra -Wpedantic -Wno-unused-parameter)

    # Crash-handler runtime needs dlsym / dladdr.
    target_link_libraries(${target} PRIVATE ${CMAKE_DL_LIBS})

    # No HDF5 BEFORE-include workaround: brew HDF5 1.14 matches what
    # HighFive expects (H5Treclaim path).
endfunction()
