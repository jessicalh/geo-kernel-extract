# Platform-Darwin.cmake — macOS-specific build settings for h5reader.
#
# Target stack (per notes/BUILD_LAYOUT_PLAN_2026-05-23.md):
#   Qt Pro 6.10.x (installer at $ENV{HOME}/Qt/<ver>/macos)
#   VTK 9.5+ built from source (typically $ENV{HOME}/VTK/)
#   HDF5 1.14 (brew: hdf5)
#   Eigen 3.4 (brew: eigen)
#
# Exact macOS install prefixes can still be supplied at configure time:
#   cmake --preset mac-rwdi \
#         -DH5READER_QT_DIR="$HOME/Qt/6.10.2/macos" \
#         -DH5READER_VTK_DIR="$HOME/VTK"
#
# Defines the same function surface as Platform-Linux/Windows so
# CMakeLists.txt is platform-agnostic.

include_guard(GLOBAL)

if(NOT H5READER_QT_DIR)
    file(GLOB _qt_candidates LIST_DIRECTORIES true "$ENV{HOME}/Qt/6.*/macos")
    list(SORT _qt_candidates COMPARE NATURAL ORDER DESCENDING)
    foreach(_qt_dir IN LISTS _qt_candidates)
        if(EXISTS "${_qt_dir}")
            set(H5READER_QT_DIR "${_qt_dir}"
                CACHE PATH "Qt installation root")
            message(STATUS
                "Platform-Darwin: auto-detected Qt at ${H5READER_QT_DIR}")
            break()
        endif()
    endforeach()
endif()

if(NOT H5READER_VTK_DIR)
    foreach(_vtk_dir IN ITEMS
            "$ENV{HOME}/VTK"
            "$ENV{HOME}/Projects/VTK"
            "/opt/vtk"
            "/opt/homebrew/opt/vtk")
        if(EXISTS "${_vtk_dir}")
            set(H5READER_VTK_DIR "${_vtk_dir}"
                CACHE PATH "VTK install root")
            message(STATUS
                "Platform-Darwin: auto-detected VTK at ${H5READER_VTK_DIR}")
            break()
        endif()
    endforeach()
endif()

if(NOT H5READER_HDF5_DIR)
    foreach(_hdf5_dir IN ITEMS
            "/opt/homebrew/opt/hdf5"
            "/usr/local/opt/hdf5")
        if(EXISTS "${_hdf5_dir}")
            set(H5READER_HDF5_DIR "${_hdf5_dir}"
                CACHE PATH "HDF5 installation root")
            message(STATUS
                "Platform-Darwin: auto-detected HDF5 at ${H5READER_HDF5_DIR}")
            break()
        endif()
    endforeach()
endif()

if(NOT H5READER_EIGEN_DIR)
    foreach(_eigen_dir IN ITEMS
            "/opt/homebrew/share/eigen3/cmake"
            "/usr/local/share/eigen3/cmake")
        if(EXISTS "${_eigen_dir}")
            set(H5READER_EIGEN_DIR "${_eigen_dir}"
                CACHE PATH "Eigen3 CMake package directory")
            message(STATUS
                "Platform-Darwin: auto-detected Eigen at ${H5READER_EIGEN_DIR}")
            break()
        endif()
    endforeach()
endif()

function(h5reader_apply_platform_target_settings target)
    target_compile_options(${target} PRIVATE
        -Wall -Wextra -Wpedantic -Wno-unused-parameter)

    # Crash-handler runtime needs dlsym / dladdr.
    target_link_libraries(${target} PRIVATE ${CMAKE_DL_LIBS})

    # No HDF5 BEFORE-include workaround: brew HDF5 1.14 matches what
    # HighFive expects (H5Treclaim path).
endfunction()
