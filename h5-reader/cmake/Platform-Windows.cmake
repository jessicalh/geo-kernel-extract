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
    file(GLOB _qt_candidates LIST_DIRECTORIES true "C:/Qt/6.*/msvc2022_64")
    list(SORT _qt_candidates COMPARE NATURAL ORDER DESCENDING)
    foreach(_qt_dir IN LISTS _qt_candidates)
        if(EXISTS "${_qt_dir}")
            set(H5READER_QT_DIR "${_qt_dir}"
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
    # qt_add_executable() correctly makes release-style Windows builds GUI
    # applications. Keep Debug as a console subsystem so structured logging is
    # visible when the Windows machine is doing first-pass bring-up.
    set_target_properties(${target} PROPERTIES
        WIN32_EXECUTABLE $<NOT:$<CONFIG:Debug>>)

    if(MSVC)
        # MSVC RWDI defaults to /O2 /Ob1 /Zi with no /GL — so cl never
        # sees across TU boundaries, and the BS/HM kernel call chain
        # (RecomputeRingScalars → EvaluateShielding → JohnsonBoveyField
        # → WireSegmentField, four TUs deep) pays a function-call
        # frame per kernel evaluation. With the field grid overlay
        # turned on that's 128k×4 frames per playback frame, which we
        # measured at ~2500 ms on Strix Halo. Whole-program optimisation
        # / link-time codegen lets cl inline the chain into the inner
        # loop and use the Eigen Vec3 specialisations end-to-end.
        #
        # BuildType-Release.cmake already enables IPO for Release. Use
        # the config-specific target property here so this works for both
        # single-config Ninja and Visual Studio / Ninja Multi-Config.
        include(CheckIPOSupported)
        check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_err LANGUAGES CXX)
        if(_ipo_ok)
            set_target_properties(${target} PROPERTIES
                INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
            message(STATUS
                "h5reader: Windows RelWithDebInfo; LTO/IPO enabled (/GL + /LTCG)")
        else()
            message(STATUS
                "h5reader: Windows RelWithDebInfo; LTO/IPO unavailable (${_ipo_err})")
        endif()
    endif()

    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)

        # MSVC's <cmath> does not expose M_PI etc. without this define.
        # The overlay code (QtBFieldStreamOverlay, QtRingPolygonOverlay)
        # uses M_PI; setting it on the command line guarantees it is
        # visible before any TU includes <cmath>. Harmless on Linux/macOS;
        # we only set it here so the Linux/macOS modules stay untouched.
        target_compile_definitions(${target} PRIVATE _USE_MATH_DEFINES)

        # Some Windows headers leak min/max macros that collide with
        # std::min / std::max and Qt's containers. NOMINMAX disables them.
        target_compile_definitions(${target} PRIVATE NOMINMAX)

        # AVX2 baseline for MSVC. Eigen's vectorised paths and the
        # Biot-Savart inner loops (QtBiotSavartCalc, QtHaighMallionCalc,
        # QtBFieldStreamOverlay grid eval) all benefit. Strix Halo /
        # Zen 4 / Zen 5 / any Haswell+ Intel supports AVX2 natively;
        # the adviser-class Win11 machines this binary targets all
        # post-date Haswell (2013). On pre-Haswell hosts the binary
        # will refuse to start with an "illegal instruction" — same
        # behaviour as Linux when /march=haswell+ is set.
        #
        # Skipping /arch:AVX512 deliberately: Zen 5 supports it but
        # not all reviewer machines will, and the marginal win over
        # AVX2 for our workload is small. Revisit if profiling shows
        # the BS grid eval is still the bottleneck.
        target_compile_options(${target} PRIVATE /arch:AVX2)
    endif()

    # Crash-handler minidump support.
    target_link_libraries(${target} PRIVATE Dbghelp)
endfunction()
