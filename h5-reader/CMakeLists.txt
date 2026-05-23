# h5-reader — standalone Qt6/VTK reader for the per-TR-emits-its-own-group
# trajectory.h5 + 5-NPY sidecar format produced by `nmr_extract --trajectory`
# post-2026-05-13.
#
# Builds independently of the nmr_shielding library.
#
# HighFive is header-only and pointed at via HIGHFIVE_INCLUDE_DIR
# (vendored at ../extern/HighFive/include). HDF5 C library is a system
# install. Eigen3 + Qt6 are system-found.
#
# Runs on Linux, macOS, and Windows. See CMakePresets.json.

cmake_minimum_required(VERSION 3.21)
project(h5reader LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    message(STATUS "Build type defaulted to RelWithDebInfo")
endif()

# Qt6 auto-gen — must be set BEFORE add_executable, otherwise Q_OBJECT
# classes silently fail to generate moc and the link fails with
# unresolved vtable symbols. Hard-won lesson from the viewer.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

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

    # Typed object model — QtProtein / QtConformation / QtFrame
    # mirroring the nmr_shielding library shape, populated from the
    # sidecar substrate + per-TR H5 buffers.
    src/model/Types.h
    src/model/QtSemanticEnums.h
    src/model/QtChainAddress.h
    src/model/QtAminoAcidType.h
    src/model/QtAminoAcidType.cpp
    src/model/QtAtom.h
    src/model/QtAtomNames.h
    src/model/QtResidue.h
    src/model/QtResidueNames.h
    src/model/QtResidueNames.cpp
    src/model/QtBond.h
    src/model/QtRing.h
    src/model/QtRing.cpp
    src/model/QtRingMembership.h
    src/model/QtTopology.h
    src/model/QtTopology.cpp
    src/model/QtProtein.h
    src/model/QtProtein.cpp

    # Per-TR buffer family (eager-loaded from trajectory.h5)
    src/model/QtTimeSeriesBuffers.h
    src/model/QtPerResidueBuffers.h
    src/model/QtWelfordBuffers.h
    src/model/QtSpecialBuffers.h
    src/model/QtSelectionBag.h
    src/model/QtSelectionBag.cpp

    src/model/QtConformation.h
    src/model/QtConformation.cpp
    src/model/QtFrame.h
    src/model/QtFrame.cpp
    src/model/QtFrameAtomView.h

    # H5 + sidecar boundary
    src/io/QtNpyRecords.h
    src/io/QtNpyReader.h
    src/io/QtNpyReader.cpp
    src/io/QtManifest.h
    src/io/QtManifest.cpp
    src/io/QtEnumVocab.h
    src/io/QtEnumVocab.cpp
    src/io/QtTopologySidecar.h
    src/io/QtTopologySidecar.cpp
    src/io/QtTrajectoryH5.h
    src/io/QtTrajectoryH5.cpp
    # QtNamingRegistry.{h,cpp} retired 2026-05-23 — vestigial CHARMM↔IUPAC
    # translation; the new sidecar format carries typed AmberAtomName /
    # IupacAtomName / BmrbAtomName per atom in atoms_category_info.npy.
    # No consumers. Files remain on disk for now; deletion in a future session.
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

    src/app/QtAtomPicker.h
    src/app/QtAtomPicker.cpp
    src/app/QtSelectionOverlay.h
    src/app/QtSelectionOverlay.cpp
    src/app/QtAtomInspectorDock.h
    src/app/QtAtomInspectorDock.cpp
    src/app/QtAtomTimeSeriesDock.h
    src/app/QtAtomTimeSeriesDock.cpp

    # Closed-form kernel evaluators
    src/calculators/QtPhysicalConstants.h
    src/calculators/QtBiotSavartCalc.h
    src/calculators/QtBiotSavartCalc.cpp
    src/calculators/QtHaighMallionCalc.h
    src/calculators/QtHaighMallionCalc.cpp

    # Entry point
    src/main_reader.cpp
)

target_include_directories(h5reader PRIVATE
    src
    ${HIGHFIVE_INCLUDE_DIR}
)

# Hard-won include-order fix from the existing viewer. HighFive picks
# H5Dvlen_reclaim vs H5Treclaim from the HDF5 header version it sees;
# system 1.10 vs micromamba 1.14 mismatch causes link failures. Prepend
# the system 1.10 include path.
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
target_link_libraries(h5reader PRIVATE ${CMAKE_DL_LIBS})
if(WIN32)
    target_link_libraries(h5reader PRIVATE Dbghelp)
endif()

target_compile_definitions(h5reader PRIVATE
    H5READER_VERSION="0.5.0-rewrite"
)

# Compiler hygiene
if(MSVC)
    target_compile_options(h5reader PRIVATE /W4 /permissive-)
else()
    target_compile_options(h5reader PRIVATE
        -Wall -Wextra -Wpedantic -Wno-unused-parameter)
endif()

message(STATUS "h5reader: Qt ${Qt6_VERSION}, VTK ${VTK_VERSION}, HDF5 ${HDF5_VERSION}")
