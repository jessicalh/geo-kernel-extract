# Platform-Windows.cmake — Windows-specific build settings for h5reader.
#
# Target stack (per notes/BUILD_LAYOUT_PLAN_2026-05-23.md):
#   Qt Pro latest (installer at C:\Qt\<ver>\msvc2022_64)
#   VTK latest, built from source
#   HDF5 1.14 (vcpkg: hdf5)
#   Eigen 3.4 (vcpkg: eigen3)
#
# Stub at migration step 1. Populated in step 5 when the Windows
# environment is first synced. Defines the same function surface as
# Platform-Linux so CMakeLists.txt is platform-agnostic.

include_guard(GLOBAL)

function(h5reader_apply_platform_target_settings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    endif()

    # Crash-handler minidump support.
    target_link_libraries(${target} PRIVATE Dbghelp)

    # No HDF5 BEFORE-include workaround: vcpkg HDF5 1.14 matches what
    # HighFive expects (H5Treclaim path).
endfunction()
