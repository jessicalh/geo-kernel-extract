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

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================================
# Platform module + dependency discovery + per-build-type flags
# ============================================================================

include(Platform-${CMAKE_SYSTEM_NAME})
include(DiscoverDeps)
qt_standard_project_setup()

# BuildType-Debug / BuildType-Release exist; RelWithDebInfo + MinSizeRel
# use CMake's defaults (no module). OPTIONAL so missing modules don't
# error — the function-defined check below decides whether to apply.
include(BuildType-${CMAKE_BUILD_TYPE} OPTIONAL)

# ============================================================================
# Target
# ============================================================================

qt_add_executable(h5reader
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
    src/model/DashboardSignal.h
    src/model/DashboardSignal.cpp
    src/model/DashboardSignalModel.h
    src/model/DashboardSignalModel.cpp
    src/model/DashboardPanelModel.h
    src/model/DashboardPanelModel.cpp
    src/model/TrajectorySignalCatalog.h
    src/model/TrajectorySignalCatalog.cpp
    src/model/SignalTimeSeries.h
    src/model/SignalTimeSeries.cpp
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
    src/io/DftShieldingLoader.h
    src/io/DftShieldingLoader.cpp
    src/io/TrajectoryFrameMap.h
    src/io/TrajectoryFrameMap.cpp

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
    src/app/QtAtomInspectorDock.h
    src/app/QtAtomInspectorDock.cpp
    src/app/AbstractStripPanel.h
    src/app/AbstractStripPanel.cpp
    src/app/SequenceBarPanel.h
    src/app/SequenceBarPanel.cpp
    src/app/PowerSpectrumPanel.h
    src/app/PowerSpectrumPanel.cpp
    src/app/LagDecayPanel.h
    src/app/LagDecayPanel.cpp
    src/app/ChordCouplingPanel.h
    src/app/ChordCouplingPanel.cpp
    src/app/FixedFreqPanel.h
    src/app/FixedFreqPanel.cpp
    src/app/TensorGlyphMath.h
    src/app/StripStackWidget.h
    src/app/StripStackWidget.cpp
    src/app/RestServer.h
    src/app/RestServer.cpp

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
    Qt6::HttpServer
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


include(CTest)
if(BUILD_TESTING AND Qt6Test_FOUND)
    qt_add_executable(h5reader_model_tests
        tests/dashboard_model_tests.cpp
        src/model/DashboardSignal.h
        src/model/DashboardSignal.cpp
        src/model/DashboardSignalModel.h
        src/model/DashboardSignalModel.cpp
        src/model/DashboardPanelModel.h
        src/model/DashboardPanelModel.cpp
        src/model/QtRing.h
        src/model/QtRing.cpp
        src/model/QtTopology.h
        src/model/QtTopology.cpp
        src/model/SignalDictionary.h
        src/model/TrajectorySignalCatalog.h
        src/model/TrajectorySignalCatalog.cpp
    )
    target_include_directories(h5reader_model_tests PRIVATE src)
    target_link_libraries(h5reader_model_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_model_tests COMMAND h5reader_model_tests)

    # Pure-function plane geometry tests (PlaneFrameMath.h). Header-only
    # math; no VTK / Qt widgets / I/O — links Qt6::Test + Eigen3 only.
    qt_add_executable(h5reader_scene_math_tests
        tests/scene_math_tests.cpp
        src/app/PlaneFrameMath.h
    )
    target_include_directories(h5reader_scene_math_tests PRIVATE src)
    target_link_libraries(h5reader_scene_math_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_scene_math_tests COMMAND h5reader_scene_math_tests)

    qt_add_executable(h5reader_io_tests
        tests/qttrajectoryh5_tests.cpp
        src/io/QtTrajectoryH5.h
        src/io/QtTrajectoryH5.cpp
        src/model/QtSelectionBag.h
        src/model/QtSelectionBag.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_io_tests PRIVATE
        src
        ${HIGHFIVE_INCLUDE_DIR}
    )
    target_link_libraries(h5reader_io_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
        ${HDF5_C_LIBRARIES}
    )
    h5reader_apply_platform_target_settings(h5reader_io_tests)
    add_test(NAME h5reader_io_tests COMMAND h5reader_io_tests)

    qt_add_executable(h5reader_app_tests
        tests/dashboard_controller_tests.cpp
        src/app/DashboardDisplayController.h
        src/app/DashboardDisplayController.cpp
        src/app/AbstractStripPanel.h
        src/app/AbstractStripPanel.cpp
        src/app/SequenceBarPanel.h
        src/app/SequenceBarPanel.cpp
        src/app/PowerSpectrumPanel.h
        src/app/PowerSpectrumPanel.cpp
        src/app/LagDecayPanel.h
        src/app/LagDecayPanel.cpp
        src/app/ChordCouplingPanel.h
        src/app/ChordCouplingPanel.cpp
        src/app/FixedFreqPanel.h
        src/app/FixedFreqPanel.cpp
        src/app/SceneRevealOverlay.h
        src/app/SceneRevealOverlay.cpp
        src/model/Conformation.h
        src/model/Conformation.cpp
        src/model/ConformationGeometry.h
        src/model/ConformationGeometry.cpp
        src/model/DashboardSignal.h
        src/model/DashboardSignal.cpp
        src/model/DashboardSignalModel.h
        src/model/DashboardSignalModel.cpp
        src/model/DashboardPanelModel.h
        src/model/DashboardPanelModel.cpp
        src/model/AtomSelection.h
        src/model/AtomSelection.cpp
        src/model/TrajectorySignalCatalog.h
        src/model/TrajectorySignalCatalog.cpp
        src/model/SignalTimeSeries.h
        src/model/SignalTimeSeries.cpp
        src/model/DftShielding.h
        src/model/DftShieldingStore.h
        src/model/DftShieldingStore.cpp
        src/model/SphericalDecomposition.h
        src/model/SphericalDecomposition.cpp
        src/io/DftShieldingLoader.h
        src/io/DftShieldingLoader.cpp
        src/model/QtProtein.h
        src/model/QtProtein.cpp
        src/model/QtRing.h
        src/model/QtRing.cpp
        src/model/QtTopology.h
        src/model/QtTopology.cpp
        src/model/QtAminoAcidType.h
        src/model/QtAminoAcidType.cpp
        src/model/QtResidueNames.h
        src/model/QtResidueNames.cpp
        src/io/OrcaShieldingParser.h
        src/io/OrcaShieldingParser.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
        src/diagnostics/ObjectCensus.h
        src/diagnostics/ObjectCensus.cpp
        src/diagnostics/Categories.cpp
    )
    target_include_directories(h5reader_app_tests PRIVATE
        src
        ${HIGHFIVE_INCLUDE_DIR}
    )
    target_link_libraries(h5reader_app_tests PRIVATE
        Qt6::Core
        Qt6::Widgets
        Qt6::Test
        Eigen3::Eigen
        ${VTK_LIBRARIES}
    )
    vtk_module_autoinit(TARGETS h5reader_app_tests MODULES ${VTK_LIBRARIES})
    h5reader_apply_platform_target_settings(h5reader_app_tests)
    add_test(NAME h5reader_app_tests COMMAND h5reader_app_tests)

    qt_add_executable(h5reader_dft_shielding_loader_tests
        tests/dft_shielding_loader_tests.cpp
        src/io/DftShieldingLoader.h
        src/io/DftShieldingLoader.cpp
        src/io/OrcaShieldingParser.h
        src/io/OrcaShieldingParser.cpp
        src/model/DftShielding.h
        src/model/SphericalDecomposition.h
        src/model/SphericalDecomposition.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_dft_shielding_loader_tests PRIVATE src)
    target_link_libraries(h5reader_dft_shielding_loader_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    h5reader_apply_platform_target_settings(h5reader_dft_shielding_loader_tests)
    add_test(NAME h5reader_dft_shielding_loader_tests COMMAND h5reader_dft_shielding_loader_tests)

    qt_add_executable(h5reader_trajectory_frame_map_tests
        tests/trajectory_frame_map_tests.cpp
        src/io/TrajectoryFrameMap.h
        src/io/TrajectoryFrameMap.cpp
        src/io/QtTrajectoryH5.h
        src/io/QtTrajectoryH5.cpp
        src/model/QtSelectionBag.h
        src/model/QtSelectionBag.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_trajectory_frame_map_tests PRIVATE
        src
        ${HIGHFIVE_INCLUDE_DIR}
    )
    target_link_libraries(h5reader_trajectory_frame_map_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
        ${HDF5_C_LIBRARIES}
    )
    h5reader_apply_platform_target_settings(h5reader_trajectory_frame_map_tests)
    add_test(NAME h5reader_trajectory_frame_map_tests COMMAND h5reader_trajectory_frame_map_tests)
endif()

# REST end-to-end smoke: pytest drives the built h5reader binary via its
# --rest <port> HTTP surface against a fixture trajectory. Replaces the
# bespoke --dashboard-path-smoke / --camera-plane-lock-smoke CLI runners.
# Requires Python3 + pytest + httpx + Pillow on PATH; H5READER_REST_FIXTURE
# env var (or the default below) must point at a valid trajectory directory.
if(BUILD_TESTING AND Python3_Interpreter_FOUND)
    set(H5READER_REST_FIXTURE_DEFAULT
        "/shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft"
        CACHE PATH "Default H5READER_REST_FIXTURE for the REST smoke test")
    add_test(NAME h5reader_rest_smoke
             COMMAND ${Python3_EXECUTABLE} -m pytest
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/rest -v --tb=short)
    # We deliberately do NOT set QT_QPA_PLATFORM=offscreen here — VTK's
    # QOpenGLWidget needs a real GL FBO, which the offscreen platform
    # does not provide. The pytest fixture detects missing DISPLAY and
    # wraps the binary in `xvfb-run -s "-screen 0 1280x720x24 +extension GLX"`
    # with LIBGL_ALWAYS_SOFTWARE=1 (Mesa software path) for true headless
    # runs. When DISPLAY is set (dev machine, ctest inherited env), the
    # binary runs against the real display.
    set_tests_properties(h5reader_rest_smoke PROPERTIES
        ENVIRONMENT
            "H5READER_BINARY=$<TARGET_FILE:h5reader>;H5READER_REST_FIXTURE=${H5READER_REST_FIXTURE_DEFAULT}"
        LABELS "rest;e2e"
        TIMEOUT 300
    )
endif()

message(STATUS "h5reader: Qt ${Qt6_VERSION}, VTK ${VTK_VERSION}, HDF5 ${HDF5_VERSION}")
