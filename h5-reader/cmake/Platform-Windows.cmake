# Platform-Windows.cmake — Windows-specific build settings for h5reader.
#
# Target stack (per notes/BUILD_LAYOUT_PLAN_2026-05-23.md):
#   Qt Pro 6.10.x   (installer → C:\Qt\<ver>\msvc2022_64)
#   VTK 9.5+        (built from source, CMake INSTALL prefix → C:\Projects\VTK)
#   HDF5 1.14       (vcpkg: hdf5)
#   Eigen 3.4       (vcpkg: eigen3)
#   CMake + Ninja   (Qt installer bundles both at C:\Qt\Tools\)
#
# Defaults below match the standard developer environment described
# in the qt6-cpp skill notes. Override per-machine via:
#   cmake --preset win-rwdi `
#         -DH5READER_QT_DIR="C:\Qt\6.10.2\msvc2022_64" `
#         -DH5READER_VTK_DIR="C:\Projects\VTK"
#
# vcpkg integration (HDF5 + Eigen3 via find_package): pass the
# toolchain file at first configure; it's cached after that:
#   cmake --preset win-rwdi `
#         -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake"
# Or set $env:VCPKG_ROOT in your PowerShell profile and rely on the
# auto-detection below.

include_guard(GLOBAL)

# Qt Pro install — prefer the newest 6.10.x we can find. Falls back
# to 6.9.x / 6.8.x for older installer caches.
if(NOT H5READER_QT_DIR)
    foreach(_qt_ver IN ITEMS
            6.10.3 6.10.2 6.10.1 6.10.0
            6.9.3  6.9.2  6.9.1  6.9.0
            6.8.0)
        if(EXISTS "C:/Qt/${_qt_ver}/msvc2022_64")
            set(H5READER_QT_DIR "C:/Qt/${_qt_ver}/msvc2022_64"
                CACHE PATH "Qt installation root")
            message(STATUS
                "Platform-Windows: auto-detected Qt at ${H5READER_QT_DIR}")
            break()
        endif()
    endforeach()
endif()

# VTK install dir (after a from-source build's INSTALL step).
if(NOT H5READER_VTK_DIR)
    foreach(_vtk_dir IN ITEMS
            "C:/Projects/VTK"
            "C:/VTK"
            "C:/Program Files/VTK")
        if(EXISTS "${_vtk_dir}")
            set(H5READER_VTK_DIR "${_vtk_dir}"
                CACHE PATH "VTK install root")
            message(STATUS
                "Platform-Windows: auto-detected VTK at ${H5READER_VTK_DIR}")
            break()
        endif()
    endforeach()
endif()

# vcpkg toolchain auto-detect IF $env:VCPKG_ROOT is set AND
# CMAKE_TOOLCHAIN_FILE was not already supplied (cmake processes
# the toolchain file before this module loads on subsequent
# configures, so the cache var will already exist).
if(NOT CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{VCPKG_ROOT})
    set(_vcpkg_toolchain
        "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
    if(EXISTS "${_vcpkg_toolchain}")
        message(STATUS
            "Platform-Windows: VCPKG_ROOT detected but CMAKE_TOOLCHAIN_FILE "
            "was not set at configure time. Re-run with "
            "-DCMAKE_TOOLCHAIN_FILE=${_vcpkg_toolchain} to pick up "
            "vcpkg-provided HDF5 + Eigen3.")
    endif()
endif()

function(h5reader_apply_platform_target_settings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    endif()

    # Crash-handler minidump support.
    target_link_libraries(${target} PRIVATE Dbghelp)
endfunction()
