# DiscoverDeps.cmake — single block of dependency discovery for
# h5reader, shared across all three platforms.
#
# Consumes per-environment hint cache variables before find_package()
# so a user can point the build at a non-default Qt / VTK / HDF5 /
# Eigen install without editing the CMakeLists or platform module:
#
#   -DH5READER_QT_DIR=$HOME/Qt/6.10.2/macos
#   -DH5READER_VTK_DIR=$HOME/VTK
#   -DH5READER_HDF5_DIR=/opt/homebrew/Cellar/hdf5/1.14.3
#   -DH5READER_EIGEN_DIR=/opt/homebrew/share/eigen3/cmake
#
# Platform-${CMAKE_SYSTEM_NAME}.cmake may set defaults for these
# before this module is included.

include_guard(GLOBAL)

if(H5READER_QT_DIR)
    list(PREPEND CMAKE_PREFIX_PATH "${H5READER_QT_DIR}")
endif()
if(H5READER_VTK_DIR)
    list(PREPEND CMAKE_PREFIX_PATH "${H5READER_VTK_DIR}")
endif()
if(H5READER_HDF5_DIR)
    set(HDF5_ROOT "${H5READER_HDF5_DIR}")
endif()
if(H5READER_EIGEN_DIR)
    set(Eigen3_DIR "${H5READER_EIGEN_DIR}")
endif()

find_package(Qt6 6.4 REQUIRED COMPONENTS
    Widgets
    OpenGLWidgets
    Network
    Charts
)

find_package(VTK 9 REQUIRED COMPONENTS
    IOChemistry
    DomainsChemistry
    DomainsChemistryOpenGL2
    FiltersSources
    FiltersCore
    FiltersGeneral
    RenderingOpenGL2
    RenderingCore
    RenderingAnnotation
    GUISupportQt
    InteractionStyle
    ImagingHybrid
    FiltersFlowPaths
    CommonColor
)

find_package(HDF5 REQUIRED COMPONENTS C)

# Eigen: apt/brew ship 3.4.x, vcpkg has moved to 5.0.x. The API we use
# (Vector3d, Matrix3d, dynamic Matrix, JacobiSVD with ComputeFullV,
# Eigen::Index) is stable from 3.4 through the 5.x line, so any of them
# is fine. We deliberately do NOT use a find_package version range
# ("3.4...<6"): CMake ranges require the package's ConfigVersion script
# to opt into range handling, and Eigen 3.4.0's shipped
# Eigen3ConfigVersion.cmake predates that — it rejects the range outright
# (verified against apt /usr/share/eigen3 on the build host). A bare
# SameMajorVersion request ("3.4") conversely rejects vcpkg's 5.0.x. So:
# find unversioned, then assert the floor ourselves. Identical on all
# three platforms.
find_package(Eigen3 REQUIRED NO_MODULE)
if(Eigen3_VERSION VERSION_LESS 3.4)
    message(FATAL_ERROR
        "Eigen >= 3.4 required (found ${Eigen3_VERSION} at ${Eigen3_DIR}).")
endif()
message(STATUS "h5reader: Eigen ${Eigen3_VERSION}")

# HighFive — header-only, vendored in the parent tree.
set(HIGHFIVE_INCLUDE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../extern/HighFive/include"
    CACHE PATH "Path to HighFive include directory")
if(NOT EXISTS "${HIGHFIVE_INCLUDE_DIR}/highfive/H5File.hpp")
    message(FATAL_ERROR
        "HighFive headers not found at ${HIGHFIVE_INCLUDE_DIR}. "
        "Clone https://github.com/BlueBrain/HighFive into extern/HighFive "
        "or pass -DHIGHFIVE_INCLUDE_DIR=...")
endif()
