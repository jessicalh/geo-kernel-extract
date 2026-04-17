# h5-reader — standalone Qt6/VTK reader for analysis H5 files.
#
# Builds independently of the nmr_shielding library. Consumes the frozen
# fileformat/analysis_file.{h,cpp} via source include (same pattern as
# the existing nmr-viewer). HighFive is header-only and pointed at via
# HIGHFIVE_INCLUDE_DIR. HDF5 C library is a system install.
#
# Runs on Linux, macOS, and Windows. See CMakePresets.json.

cmake_minimum_required(VERSION 3.21)
project(h5reader LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type when not configuring via presets.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    message(STATUS "Build type defaulted to RelWithDebInfo")
endif()

# Qt6 auto-generated code — must be set BEFORE add_executable, otherwise
# Q_OBJECT classes silently fail to generate moc and the link fails with
# unresolved vtable symbols. Hard-won lesson from the existing viewer.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Emit compile_commands.json so clangd and other language servers resolve
# Qt/VTK include paths automatically. Symlinked to project root by the
# launcher script for IDE pickup.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================================
# Dependencies
# ============================================================================

find_package(Qt6 6.4 REQUIRED COMPONENTS
    Widgets
    OpenGLWidgets
    Network
    Charts
)

# VTK component list mirrors the existing viewer. Add components as new
# overlays demand them; don't broaden to "find_package(VTK REQUIRED)" —
# module autoinit needs the exact set we link against.
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

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

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

# Frozen H5 (de)serialiser — compiled directly in, not linked as a library.
set(FILEFORMAT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../fileformat")
set(FILEFORMAT_SRC "${FILEFORMAT_DIR}/analysis_file.cpp")
if(NOT EXISTS "${FILEFORMAT_SRC}")
    message(FATAL_ERROR
        "fileformat/analysis_file.cpp not found at ${FILEFORMAT_SRC}. "
        "h5-reader must be built within the nmr-shielding tree so "
        "fileformat/ is available as a sibling directory.")
endif()

# ============================================================================
# Target
# ============================================================================

add_executable(h5reader
    # Diagnostics — installed in main() before any Qt construction.
    src/diagnostics/StructuredLogger.h
    src/diagnostics/StructuredLogger.cpp
    src/diagnostics/ErrorBus.h
    src/diagnostics/ErrorBus.cpp
    src/diagnostics/ObjectCensus.h
    src/diagnostics/ObjectCensus.cpp
    src/diagnostics/CrashHandler.h
    src/diagnostics/CrashHandler.cpp
    src/diagnostics/ConnectionAuditor.h
    src/diagnostics/ThreadGuard.h
    src/diagnostics/Categories.cpp
    src/diagnostics/ShutdownSignals.h
    src/diagnostics/ShutdownSignals.cpp

    # Typed object model — QtProtein / QtConformation / QtFrame.
    src/model/Types.h
    src/model/QtAtom.h
    src/model/QtResidue.h
    src/model/QtBond.h
    src/model/QtRing.h
    src/model/QtRing.cpp
    src/model/QtProtein.h
    src/model/QtProtein.cpp
    src/model/QtConformation.h
    src/model/QtConformation.cpp
    src/model/QtFrame.h
    src/model/QtFrame.cpp

    # H5 → typed-model boundary (the ONE place strings / ordinals decode).
    src/io/QtNamingRegistry.h
    src/io/QtNamingRegistry.cpp
    src/io/QtProteinLoader.h
    src/io/QtProteinLoader.cpp

    # App layer — window, playback, scene, overlays.
    src/app/QtPlaybackController.h
    src/app/QtPlaybackController.cpp
    src/app/MoleculeScene.h
    src/app/MoleculeScene.cpp
    src/app/ReaderMainWindow.h
    src/app/ReaderMainWindow.cpp
    src/app/QtRingPolygonOverlay.h
    src/app/QtRingPolygonOverlay.cpp
    src/app/QtBackboneRibbonOverlay.h
    src/app/QtBackboneRibbonOverlay.cpp
    src/app/QtFieldGridOverlay.h
    src/app/QtFieldGridOverlay.cpp
    src/app/QtBFieldStreamOverlay.h
    src/app/QtBFieldStreamOverlay.cpp
    src/app/Colormap.h

    # Picker + selection + inspector dock (step 6).
    src/app/QtAtomPicker.h
    src/app/QtAtomPicker.cpp
    src/app/QtSelectionOverlay.h
    src/app/QtSelectionOverlay.cpp
    src/app/QtAtomInspectorDock.h
    src/app/QtAtomInspectorDock.cpp
    src/app/QtAtomTimeSeriesDock.h
    src/app/QtAtomTimeSeriesDock.cpp

    # Closed-form kernel evaluators — free functions, thread-safe.
    src/calculators/QtPhysicalConstants.h
    src/calculators/QtBiotSavartCalc.h
    src/calculators/QtBiotSavartCalc.cpp
    src/calculators/QtHaighMallionCalc.h
    src/calculators/QtHaighMallionCalc.cpp

    # Frozen H5 serialiser (source include, not a library)
    ${FILEFORMAT_SRC}

    # Entry point
    src/main_reader.cpp
)

target_include_directories(h5reader PRIVATE
    src
    ${FILEFORMAT_DIR}
    ${HIGHFIVE_INCLUDE_DIR}
)

# Hard-won include-order fix from the existing viewer. HighFive selects
# between H5Dvlen_reclaim (<1.12) and H5Treclaim (>=1.12) based on the
# HDF5 header version it sees. On systems where OpenBabel pulls in 1.14
# headers via a micromamba path while the link library is system 1.10,
# HighFive emits H5Treclaim and the link fails. Prepending the system
# HDF5 include path forces the old API to match the link library.
# See feedback_hdf5_include_order.md.
if(UNIX AND NOT APPLE AND EXISTS /usr/include/hdf5/serial)
    target_include_directories(h5reader BEFORE PRIVATE /usr/include/hdf5/serial)
endif()

target_link_libraries(h5reader PRIVATE
    Qt6::Widgets
    Qt6::OpenGLWidgets
    Qt6::Network
    Qt6::Charts
    Eigen3::Eigen
    ${VTK_LIBRARIES}
    ${HDF5_C_LIBRARIES}
)

# Without this, VTK object factories are not registered at runtime and
# rendering silently fails. qt6-cpp skill: "#1 cause of VTK compiles but
# renders nothing." DO NOT REMOVE.
vtk_module_autoinit(TARGETS h5reader MODULES ${VTK_LIBRARIES})

# Platform-conditional crash-capture linkage.
#   Linux:  libdl for dladdr inside backtrace_symbols_fd.
#   macOS:  dl symbols live in libSystem (linked by default); CMAKE_DL_LIBS
#           resolves to empty so the line becomes a no-op.
#   Windows: Dbghelp for future MiniDumpWriteDump support.
target_link_libraries(h5reader PRIVATE ${CMAKE_DL_LIBS})
if(WIN32)
    target_link_libraries(h5reader PRIVATE Dbghelp)
endif()

target_compile_definitions(h5reader PRIVATE
    H5READER_VERSION="0.4.0"
)

# Suppress one cosmetic warning in the frozen fileformat/analysis_file.cpp:
# ReadFloat() is defined but unused in the current H5 schema (the aim
# embedding is handled specially). The library is read-only from our side;
# we do not edit fileformat/. When that file is revisited upstream, remove
# this suppression.
if(NOT MSVC)
    set_source_files_properties(${FILEFORMAT_SRC}
        PROPERTIES COMPILE_OPTIONS "-Wno-unused-function")
endif()

# Compiler hygiene
if(MSVC)
    target_compile_options(h5reader PRIVATE /W4 /permissive-)
else()
    target_compile_options(h5reader PRIVATE
        -Wall -Wextra -Wpedantic -Wno-unused-parameter)
endif()

message(STATUS "h5reader: Qt ${Qt6_VERSION}, VTK ${VTK_VERSION}, HDF5 ${HDF5_VERSION}")
