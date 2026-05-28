# h5-reader — standalone Qt6/VTK reader for the per-TR-emits-its-own-group
# trajectory.h5 + 5-NPY sidecar format produced by `nmr_extract --trajectory`
# post-2026-05-13.
#
# Builds independently of the nmr_shielding library on Linux, macOS, and
# Windows. Per-platform dependency paths and target flags live in
# cmake/Platform-${CMAKE_SYSTEM_NAME}.cmake; find_package() calls in
# cmake/DiscoverDeps.cmake. See notes/BUILD_LAYOUT_PLAN_2026-05-23.md.

cmake_minimum_required(VERSION 3.21)
project(h5reader LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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
# Platform module + dependency discovery + per-build-type flags
# ============================================================================

include(Platform-${CMAKE_SYSTEM_NAME})
include(DiscoverDeps)

# BuildType-Debug / BuildType-Release exist; RelWithDebInfo + MinSizeRel
# use CMake's defaults (no module). OPTIONAL so missing modules don't
# error — the function-defined check below decides whether to apply.
include(BuildType-${CMAKE_BUILD_TYPE} OPTIONAL)

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

    # Typed object model — QtProtein / Conformation hierarchy / QtFrame
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

    src/model/Conformation.h
    src/model/Conformation.cpp
    src/model/TrajectoryConformation.h
    src/model/TrajectoryConformation.cpp
    src/model/SingleConformation.h
    src/model/SingleConformation.cpp
    src/model/ConformationGeometry.h
    src/model/ConformationGeometry.cpp
    src/model/SpectralAnalysis.h
    src/model/SpectralAnalysis.cpp
    src/model/SignalDictionary.h
    src/model/SignalDictionary.cpp
    src/model/DashboardSignal.h
    src/model/DashboardSignal.cpp
    src/model/DashboardSignalModel.h
    src/model/DashboardSignalModel.cpp
    src/model/TrajectorySignalCatalog.h
    src/model/TrajectorySignalCatalog.cpp
    src/model/SignalTimeSeries.h
    src/model/SignalTimeSeries.cpp
    src/model/StripCalculation.h
    src/model/GeometryStrips.h
    src/model/GeometryStrips.cpp
    src/model/DftSigmaStrips.h
    src/model/DftSigmaStrips.cpp
    src/model/StripChartChannel.h
    src/model/DftShielding.h
    src/model/SphericalDecomposition.h
    src/model/SphericalDecomposition.cpp
    src/model/DftShieldingStore.h
    src/model/DftShieldingStore.cpp
    src/model/AtomSelection.h
    src/model/AtomSelection.cpp
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
    src/io/FrameNpyLoader.h
    src/io/FrameNpyLoader.cpp
    src/io/OrcaShieldingParser.h
    src/io/OrcaShieldingParser.cpp

    # App layer — window, playback, scene, overlays.
    src/app/QtPlaybackController.h
    src/app/QtPlaybackController.cpp
    src/app/TimeViewportController.h
    src/app/TimeViewportController.cpp
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
    src/app/MeasurementOverlay.h
    src/app/MeasurementOverlay.cpp
    src/app/SceneRevealOverlay.h
    src/app/SceneRevealOverlay.cpp
    src/app/SelectionDock.h
    src/app/SelectionDock.cpp
    src/app/NearbySignalModel.h
    src/app/NearbySignalModel.cpp
    src/app/SignalDisplayDialog.h
    src/app/SignalDisplayDialog.cpp
    src/app/DashboardDisplayController.h
    src/app/DashboardDisplayController.cpp
    src/app/DashboardStripDock.h
    src/app/DashboardStripDock.cpp
    src/app/SignalPickerDock.h
    src/app/SignalPickerDock.cpp
    src/app/QtAtomInspectorDock.h
    src/app/QtAtomInspectorDock.cpp
    src/app/QtAtomTimeSeriesDock.h
    src/app/QtAtomTimeSeriesDock.cpp
    src/app/StripStackWidget.h
    src/app/StripStackWidget.cpp
    src/app/StripChartDock.h
    src/app/StripChartDock.cpp

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

# Platform-specific tweaks: warning flags, crash-handler libs, HDF5
# include-order fix on Linux, etc. Defined in
# cmake/Platform-${CMAKE_SYSTEM_NAME}.cmake.
h5reader_apply_platform_target_settings(h5reader)

# Build-type tweaks: sanitizers in Debug, LTO in Release, nothing
# in RelWithDebInfo / MinSizeRel. Function only exists when the
# matching BuildType-*.cmake module was found.
if(COMMAND h5reader_apply_build_type_settings)
    h5reader_apply_build_type_settings(h5reader)
endif()

target_compile_definitions(h5reader PRIVATE
    H5READER_VERSION="0.5.0-rewrite"
)

message(STATUS "h5reader: Qt ${Qt6_VERSION}, VTK ${VTK_VERSION}, HDF5 ${HDF5_VERSION}")
