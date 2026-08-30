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
set(H5READER_SHARED_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../shared")

if(WIN32)
    option(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML
        "Require and package the Experimental Shielding ML model/runtime" ON)
    option(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML_ROCM
        "Require and package the native ROCm Experimental Shielding ML runtime" ON)
endif()

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

qt_add_library(h5reader_core STATIC
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
    src/model/TransformedConformation.h
    src/model/TransformedConformation.cpp
    src/model/ConformationGeometry.h
    src/model/ConformationGeometry.cpp
    src/model/RingNullCollar.h
    src/model/RingNullCollar.cpp
    src/model/RingCurrentFaceCollar.h
    src/model/RingCurrentFaceCollar.cpp
    src/model/SpectralAnalysis.h
    src/model/SpectralAnalysis.cpp
    src/model/SignalDictionary.h
    src/model/DashboardSignal.h
    src/model/DashboardSignal.cpp
    src/model/DashboardSignalModel.h
    src/model/DashboardSignalModel.cpp
    src/model/DashboardPanelModel.h
    src/model/DashboardPanelModel.cpp
    src/model/DisplayModeCapability.cpp
    src/model/VisualizationDefinition.h
    src/model/VisualizationDefinition.cpp
    src/model/VisualizationRegistry.h
    src/model/VisualizationRegistry.cpp
    src/model/StripVisualization.h
    src/model/StripVisualization.cpp
    src/model/SequenceBarVisualization.h
    src/model/SequenceBarVisualization.cpp
    src/model/LagCurveVisualization.h
    src/model/LagCurveVisualization.cpp
    src/model/ChordCouplingVisualization.h
    src/model/ChordCouplingVisualization.cpp
    src/model/FixedFrequencyVisualization.h
    src/model/FixedFrequencyVisualization.cpp
    src/model/PowerSpectrumVisualization.h
    src/model/PowerSpectrumVisualization.cpp
    src/model/TensorGlyphVisualization.h
    src/model/TensorGlyphVisualization.cpp
    src/model/NewmanVisualization.h
    src/model/NewmanVisualization.cpp
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
    src/model/ExperimentalShieldingMlStore.h
    src/model/ExperimentalShieldingMlStore.cpp
    src/model/AtomSelection.h
    src/model/AtomSelection.cpp
    src/model/QtFrame.h
    src/model/QtFrame.cpp

    # H5 + sidecar boundary
    src/io/QtNpyRecords.h
    src/io/QtNpyReader.h
    src/io/QtNpyReader.cpp
    src/io/QtManifest.h
    src/io/QtManifest.cpp
    src/io/CalcsetManifest.h
    src/io/CalcsetManifest.cpp
    src/io/QtEnumVocab.h
    src/io/QtEnumVocab.cpp
    src/io/QtTopologySidecar.h
    src/io/QtTopologySidecar.cpp
    src/io/QtTrajectoryH5.h
    src/io/QtTrajectoryH5.cpp
    src/io/QtProteinLoader.h
    src/io/QtProteinLoader.cpp
    src/io/FrameFieldPolicy.h
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
    src/app/CameraMode.h
    src/app/OrientationPolicy.h
    src/app/FitTargetMath.h
    src/app/CameraComposer.h
    src/app/CameraComposer.cpp
    src/app/CameraAnchorHelper.h
    src/app/CameraAnchorHelper.cpp
    src/app/QuietTrackballStyle.h
    src/app/QuietTrackballStyle.cpp
    src/app/CameraInputFilter.h
    src/app/CameraInputFilter.cpp
    src/app/MoleculeScene.h
    src/app/MoleculeScene.cpp
    src/app/NewmanProjection.h
    src/app/NewmanProjection.cpp
    src/app/NewmanPaint.h
    src/app/NewmanPaint.cpp
    src/app/MeasurementsDock.h
    src/app/MeasurementsDock.cpp
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
    src/app/OccupancyShellsMath.h
    src/app/QtAtomTrajectoryOverlay.h
    src/app/QtAtomTrajectoryOverlay.cpp
    src/app/AtomTrackOverlay.h
    src/app/AtomTrackOverlay.cpp
    src/app/HeroshotButterflyOverlay.h
    src/app/HeroshotButterflyOverlay.cpp
    src/app/HeroshotTensorPairOverlay.h
    src/app/HeroshotTensorPairOverlay.cpp
    src/app/Colormap.h

    src/app/QtAtomPicker.h
    src/app/QtAtomPicker.cpp
    src/app/MeasurementOverlay.h
    src/app/MeasurementOverlay.cpp
    src/app/SceneRevealOverlay.h
    src/app/SceneRevealOverlay.cpp
    src/app/TensorGlyphActor.h
    src/app/TensorGlyphActor.cpp
    src/app/CsaTensorOverlay.h
    src/app/CsaTensorOverlay.cpp
    src/app/TensorGhostTrail.h
    src/app/TensorGhostTrail.cpp
    src/app/AngleCollarActor.h
    src/app/AngleCollarActor.cpp
    src/app/NearbySignalModel.h
    src/app/NearbySignalModel.cpp
    src/app/SignalDisplayDialog.h
    src/app/SignalDisplayDialog.cpp
    src/app/DashboardDisplayController.h
    src/app/DashboardDisplayController.cpp
    src/app/DashboardSelectionController.h
    src/app/DashboardSelectionController.cpp
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
    src/app/NewmanPanel.h
    src/app/NewmanPanel.cpp
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

    # Physics — pure shared math reused by the viewer to derive validated
    # metrics on loaded npy (no emit path).
    src/physics/LiteratureAccessors.h
    src/physics/SphericalBasis.h
    src/physics/SphericalBasis.cpp
    src/physics/LocalFrameBasis.h
    src/physics/LocalFrameBasis.cpp
    src/physics/ClassicalSourceMath.h
    src/physics/ClassicalSourceMath.cpp
    src/physics/CircularRingCurrent.h
    src/physics/CircularRingCurrent.cpp
    src/physics/RingCurrentScalars.h
    src/physics/RingCurrentScalars.cpp
    src/physics/McConnellLiteratureKernel.h
    src/physics/McConnellLiteratureKernel.cpp
    src/physics/BuckinghamEfield.h
    src/physics/BuckinghamEfield.cpp
    src/physics/EfgFeature.h
    src/physics/EfgFeature.cpp

)

target_include_directories(h5reader_core PUBLIC
    src
    src/physics
    ${H5READER_SHARED_INCLUDE_DIR}
    ${HIGHFIVE_INCLUDE_DIR}
)

target_link_libraries(h5reader_core PUBLIC
    Qt6::Widgets
    Qt6::OpenGLWidgets
    Qt6::Network
    Qt6::Charts
    Qt6::HttpServer
    Eigen3::Eigen
    ${VTK_LIBRARIES}
    ${HDF5_C_LIBRARIES}
)

# Platform-specific tweaks: warning flags, crash-handler libs, HDF5
# include-order fix on Linux, etc. Defined in
# cmake/Platform-${CMAKE_SYSTEM_NAME}.cmake.
h5reader_apply_platform_target_settings(h5reader_core)

# Build-type tweaks: sanitizers in Debug, LTO in Release, nothing
# in RelWithDebInfo / MinSizeRel. Function only exists when the
# matching BuildType-*.cmake module was found.
if(COMMAND h5reader_apply_build_type_settings)
    h5reader_apply_build_type_settings(h5reader_core)
endif()

# Version macro is consumed by RestServer.cpp and main_reader.cpp; needs
# to be visible to every TU in the core lib. PUBLIC so executables that
# link the lib (and tests, eventually) also see it.
target_compile_definitions(h5reader_core PUBLIC
    H5READER_VERSION="0.5.0-rewrite"
)

# ============================================================================
# GUI reader executable — links h5reader_core; entry is main_reader.cpp.
# A second executable (or a dispatcher main.cpp picking between
# main_reader / main_batch) lands when the chewer's batch entry point
# arrives; both will link the same static core.
# ============================================================================

qt_add_executable(h5reader
    src/main_reader.cpp
)

target_link_libraries(h5reader PRIVATE h5reader_core)

# vtk_module_autoinit applies to the EXECUTABLE that links the static
# lib, not to the lib itself — the autoinit translation unit would be
# stripped from the static archive otherwise. qt6-cpp skill: "#1 cause
# of VTK compiles but renders nothing." DO NOT REMOVE.
vtk_module_autoinit(TARGETS h5reader MODULES ${VTK_LIBRARIES})

h5reader_apply_platform_target_settings(h5reader)

# Build-type tweaks (sanitizer flags etc.) must also apply to the
# executable so its link command picks up the sanitizer runtime that
# the lib's compiled .o files reference.
if(COMMAND h5reader_apply_build_type_settings)
    h5reader_apply_build_type_settings(h5reader)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(h5reader_physics_tests
        src/physics/PhysicsTests.cpp
        src/physics/SphericalBasis.cpp
        src/physics/LocalFrameBasis.cpp
        src/physics/ClassicalSourceMath.cpp
        src/physics/CircularRingCurrent.cpp
        src/physics/RingCurrentScalars.cpp
        src/physics/McConnellLiteratureKernel.cpp
        src/physics/BuckinghamEfield.cpp
        src/physics/EfgFeature.cpp
    )
    set_target_properties(h5reader_physics_tests PROPERTIES
        AUTOMOC OFF
        AUTOUIC OFF
        AUTORCC OFF
    )
    target_include_directories(h5reader_physics_tests PRIVATE
        src
        src/physics
        ${H5READER_SHARED_INCLUDE_DIR}
    )
    target_link_libraries(h5reader_physics_tests PRIVATE
        Eigen3::Eigen
    )
    target_compile_definitions(h5reader_physics_tests PRIVATE
        H5READER_REPO_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
    )
    add_test(NAME h5reader_physics_tests COMMAND h5reader_physics_tests)
endif()

if(BUILD_TESTING AND Qt6Test_FOUND)
    qt_add_executable(h5reader_model_tests
        tests/dashboard_model_tests.cpp
        src/model/DashboardSignal.h
        src/model/DashboardSignal.cpp
        src/model/DashboardSignalModel.h
        src/model/DashboardSignalModel.cpp
        src/model/DashboardPanelModel.h
        src/model/DashboardPanelModel.cpp
        src/model/DisplayModeCapability.cpp
        src/model/VisualizationDefinition.h
        src/model/VisualizationDefinition.cpp
        src/model/VisualizationRegistry.h
        src/model/VisualizationRegistry.cpp
        src/model/StripVisualization.h
        src/model/StripVisualization.cpp
        src/model/SequenceBarVisualization.h
        src/model/SequenceBarVisualization.cpp
        src/model/LagCurveVisualization.h
        src/model/LagCurveVisualization.cpp
        src/model/ChordCouplingVisualization.h
        src/model/ChordCouplingVisualization.cpp
        src/model/FixedFrequencyVisualization.h
        src/model/FixedFrequencyVisualization.cpp
        src/model/PowerSpectrumVisualization.h
        src/model/PowerSpectrumVisualization.cpp
        src/model/TensorGlyphVisualization.h
        src/model/TensorGlyphVisualization.cpp
        src/model/NewmanVisualization.h
        src/model/NewmanVisualization.cpp
        src/model/QtRing.h
        src/model/QtRing.cpp
        src/model/QtTopology.h
        src/model/QtTopology.cpp
        src/model/SignalDictionary.h
        src/model/TrajectorySignalCatalog.h
        src/model/TrajectorySignalCatalog.cpp
        src/diagnostics/ObjectCensus.h
        src/diagnostics/ObjectCensus.cpp
        src/diagnostics/Categories.cpp
    )
    target_include_directories(h5reader_model_tests PRIVATE
        src
        src/physics
        ${H5READER_SHARED_INCLUDE_DIR}
    )
    target_link_libraries(h5reader_model_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_model_tests COMMAND h5reader_model_tests)

    # Pure-function plane geometry tests (PlaneFrameMath.h +
    # FitTargetMath.h). Header-only math; no VTK / Qt widgets / I/O —
    # links Qt6::Test + Eigen3 only.
    qt_add_executable(h5reader_scene_math_tests
        tests/scene_math_tests.cpp
        src/app/PlaneFrameMath.h
        src/app/FitTargetMath.h
    )
    target_include_directories(h5reader_scene_math_tests PRIVATE src)
    target_link_libraries(h5reader_scene_math_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_scene_math_tests COMMAND h5reader_scene_math_tests)

    # Pure-function CSA / molecular-frame math (CsaShape.h + MolecularFrame.h) --
    # the viewer-side tensor decomposition/projection grabbed from AnalysisAtom,
    # for the standalone "look at a bare extractor run" path. Header-only math;
    # links Qt6::Test + Eigen3 only.
    qt_add_executable(h5reader_csa_math_tests
        tests/csa_math_tests.cpp
        src/model/CsaShape.h
        src/model/MolecularFrame.h
    )
    target_include_directories(h5reader_csa_math_tests PRIVATE src)
    target_link_libraries(h5reader_csa_math_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_csa_math_tests COMMAND h5reader_csa_math_tests)

    # Availability-gate classification policy (TrajectoryFieldAvailability) --
    # the pure classifyDft / classifyTopology / visible-state / conservative-
    # default logic of leg #1 (make the catalog tell the truth about what is
    # live). Build() itself is coupled to H5/NPY I/O and is REST-integration-
    # tested; this pins the header-only policy. Links Qt6::Test only.
    qt_add_executable(h5reader_field_availability_tests
        tests/field_availability_tests.cpp
        src/model/TrajectoryFieldAvailability.h
    )
    target_include_directories(h5reader_field_availability_tests PRIVATE
        src
        src/physics
        ${H5READER_SHARED_INCLUDE_DIR}
    )
    target_link_libraries(h5reader_field_availability_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_field_availability_tests
             COMMAND h5reader_field_availability_tests)

    # Metric taxonomy (mechanism -> concept -> form) -- the pure classification
    # that organizes the flat catalog into a navigable, hypothesis-first tree
    # (the metrics-clarity de-sprawl). Header-only; links Qt6::Test only.
    qt_add_executable(h5reader_metric_taxonomy_tests
        tests/metric_taxonomy_tests.cpp
        src/model/MetricTaxonomy.h
    )
    target_include_directories(h5reader_metric_taxonomy_tests PRIVATE src)
    target_link_libraries(h5reader_metric_taxonomy_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_metric_taxonomy_tests
             COMMAND h5reader_metric_taxonomy_tests)

    # Display policy (per-field presentation limits) -- A2: the data-driven
    # displayability rule (topology tables + the 256-d embedding non-displayable,
    # so the audit splits intentional non-display from real dead options). Pure.
    qt_add_executable(h5reader_display_policy_tests
        tests/display_policy_tests.cpp
        src/model/DisplayPolicy.h
    )
    target_include_directories(h5reader_display_policy_tests PRIVATE src)
    target_link_libraries(h5reader_display_policy_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_display_policy_tests
             COMMAND h5reader_display_policy_tests)

    qt_add_executable(h5reader_frame_field_policy_tests
        tests/frame_field_policy_tests.cpp
        src/io/FrameFieldPolicy.h
        src/io/QtFieldCatalog.gen.h
    )
    target_include_directories(h5reader_frame_field_policy_tests PRIVATE src)
    target_link_libraries(h5reader_frame_field_policy_tests PRIVATE
        Qt6::Core
        Qt6::Test
    )
    add_test(NAME h5reader_frame_field_policy_tests
             COMMAND h5reader_frame_field_policy_tests)

    # Pure-function density / HDR math used by the selected-atom trajectory
    # envelope. Header-only math over model::Vec3/Mat3; links Qt6::Test +
    # Eigen3 only. The viewer contours the resulting field downstream, where
    # VTK is the right tool.
    qt_add_executable(h5reader_occupancy_shells_math_tests
        tests/occupancy_shells_math_tests.cpp
        src/app/OccupancyShellsMath.h
    )
    target_include_directories(h5reader_occupancy_shells_math_tests PRIVATE src)
    target_link_libraries(h5reader_occupancy_shells_math_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_occupancy_shells_math_tests
             COMMAND h5reader_occupancy_shells_math_tests)

    # CalcsetManifest tests — `.LGS` JSON loader + path validation
    # (replaces the retired ReaderInputManifest TOML loader). Pure Qt /
    # Qt::Test; no VTK.
    qt_add_executable(h5reader_calcset_manifest_tests
        tests/calcset_manifest_tests.cpp
        src/io/CalcsetManifest.h
        src/io/CalcsetManifest.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_calcset_manifest_tests PRIVATE
        src
    )
    target_link_libraries(h5reader_calcset_manifest_tests PRIVATE
        Qt6::Core
        Qt6::Test
    )
    add_test(NAME h5reader_calcset_manifest_tests
             COMMAND h5reader_calcset_manifest_tests)

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

    qt_add_executable(h5reader_npy_reader_tests
        tests/npy_reader_tests.cpp
        src/io/QtNpyReader.h
        src/io/QtNpyReader.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_npy_reader_tests PRIVATE src)
    target_link_libraries(h5reader_npy_reader_tests PRIVATE
        Qt6::Core
        Qt6::Test
    )
    h5reader_apply_platform_target_settings(h5reader_npy_reader_tests)
    add_test(NAME h5reader_npy_reader_tests COMMAND h5reader_npy_reader_tests)

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
        src/app/NewmanPanel.h
        src/app/NewmanPanel.cpp
        src/app/NewmanPaint.h
        src/app/NewmanPaint.cpp
        src/app/NewmanProjection.h
        src/app/NewmanProjection.cpp
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
        src/model/DisplayModeCapability.cpp
        src/model/VisualizationDefinition.h
        src/model/VisualizationDefinition.cpp
        src/model/VisualizationRegistry.h
        src/model/VisualizationRegistry.cpp
        src/model/StripVisualization.h
        src/model/StripVisualization.cpp
        src/model/SequenceBarVisualization.h
        src/model/SequenceBarVisualization.cpp
        src/model/LagCurveVisualization.h
        src/model/LagCurveVisualization.cpp
        src/model/ChordCouplingVisualization.h
        src/model/ChordCouplingVisualization.cpp
        src/model/FixedFrequencyVisualization.h
        src/model/FixedFrequencyVisualization.cpp
        src/model/PowerSpectrumVisualization.h
        src/model/PowerSpectrumVisualization.cpp
        src/model/TensorGlyphVisualization.h
        src/model/TensorGlyphVisualization.cpp
        src/model/NewmanVisualization.h
        src/model/NewmanVisualization.cpp
        src/model/AtomSelection.h
        src/model/AtomSelection.cpp
        src/model/TrajectorySignalCatalog.h
        src/model/TrajectorySignalCatalog.cpp
        src/model/SignalTimeSeries.h
        src/model/SignalTimeSeries.cpp
        src/model/DftShielding.h
        src/model/DftShieldingStore.h
        src/model/DftShieldingStore.cpp
        src/model/ExperimentalShieldingMlStore.h
        src/model/ExperimentalShieldingMlStore.cpp
        src/model/SphericalDecomposition.h
        src/model/SphericalDecomposition.cpp
        src/io/DftShieldingLoader.h
        src/io/DftShieldingLoader.cpp
        src/io/CalcsetManifest.h
        src/io/CalcsetManifest.cpp
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
        src/physics
        ${H5READER_SHARED_INCLUDE_DIR}
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
    target_include_directories(h5reader_dft_shielding_loader_tests PRIVATE
        src
        src/physics
        ${H5READER_SHARED_INCLUDE_DIR}
    )
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

    qt_add_executable(h5reader_ring_current_case_tests
        tests/ring_current_case_tests.cpp
    )
    target_link_libraries(h5reader_ring_current_case_tests PRIVATE
        h5reader_core
        Qt6::Test
    )
    vtk_module_autoinit(TARGETS h5reader_ring_current_case_tests MODULES ${VTK_LIBRARIES})
    h5reader_apply_platform_target_settings(h5reader_ring_current_case_tests)
    add_test(NAME h5reader_ring_current_case_tests
             COMMAND h5reader_ring_current_case_tests)
    set_tests_properties(h5reader_ring_current_case_tests PROPERTIES TIMEOUT 900)
endif()

# REST end-to-end smoke: pytest drives the built h5reader binary via its
# --rest <port> HTTP surface against a fixture trajectory. Replaces the
# bespoke --dashboard-path-smoke / --camera-plane-lock-smoke CLI runners.
# Requires Python3 + pytest + httpx + Pillow on PATH; H5READER_REST_FIXTURE
# env var (or the default below) must point at a valid trajectory directory.
if(BUILD_TESTING AND Python3_Interpreter_FOUND)
    add_test(NAME h5reader_field_catalog_generated_check
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/gen_field_catalog.py
                     --check)
    set_tests_properties(h5reader_field_catalog_generated_check PROPERTIES
        LABELS "catalog;contract"
        TIMEOUT 60
    )

    if(WIN32)
        set(H5READER_REST_FIXTURE_DEFAULT
            "C:/projects/reader-data/1p9j-calibration-with-dft/1p9j-july20260717.LGS"
            CACHE FILEPATH "Default H5READER_REST_FIXTURE for the REST smoke test")
        set(H5READER_F006_STATIC_FIXTURE_DEFAULT
            "C:/projects/reader-data/july-20260717-staging/full720-f006-A0A822IKI2.LGS"
            CACHE FILEPATH "Complete July full720 single-pose fixture for F006 inference tests")
    else()
        set(H5READER_REST_FIXTURE_DEFAULT
            "/shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft"
            CACHE PATH "Default H5READER_REST_FIXTURE for the REST smoke test")
    endif()
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

    if(WIN32)
        add_test(NAME h5reader_rest_installed_ml_stale_dev_env
                 COMMAND ${Python3_EXECUTABLE} -m pytest
                         ${CMAKE_CURRENT_SOURCE_DIR}/tests/rest/test_experimental_shielding_ml_runtime.py
                         -v --tb=short)
        set_tests_properties(h5reader_rest_installed_ml_stale_dev_env PROPERTIES
            ENVIRONMENT
                "H5READER_BINARY=$<TARGET_FILE:h5reader>;H5READER_REST_FIXTURE=${H5READER_F006_STATIC_FIXTURE_DEFAULT};H5READER_EXPECT_ML_SUCCESS=1;H5READER_EXPECT_STALE_DEV_ML_FALLBACK=1;H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL=C:/h5reader-missing-experimental-ml/model.ts;H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST=C:/h5reader-missing-experimental-ml/manifest.json;H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER=C:/h5reader-missing-experimental-ml/infer.exe"
            LABELS "rest;e2e;ml"
            TIMEOUT 120
        )

        add_test(NAME h5reader_rest_ml_cpu_fallback
                 COMMAND ${Python3_EXECUTABLE} -m pytest
                         ${CMAKE_CURRENT_SOURCE_DIR}/tests/rest/test_experimental_shielding_ml_runtime.py
                         -v --tb=short)
        set_tests_properties(h5reader_rest_ml_cpu_fallback PROPERTIES
            ENVIRONMENT
                "H5READER_BINARY=$<TARGET_FILE:h5reader>;H5READER_REST_FIXTURE=${H5READER_F006_STATIC_FIXTURE_DEFAULT};H5READER_EXPECT_ML_SUCCESS=1;HIP_VISIBLE_DEVICES=-1;H5READER_EXPECT_ML_CPU_FALLBACK=1"
            LABELS "rest;e2e;ml"
            TIMEOUT 120
        )

        if(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML_ROCM)
            add_test(NAME h5reader_rest_ml_rocm
                     COMMAND ${Python3_EXECUTABLE} -m pytest
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests/rest/test_experimental_shielding_ml_runtime.py
                             -v --tb=short)
            set_tests_properties(h5reader_rest_ml_rocm PROPERTIES
                ENVIRONMENT
                    "H5READER_BINARY=$<TARGET_FILE:h5reader>;H5READER_REST_FIXTURE=${H5READER_F006_STATIC_FIXTURE_DEFAULT};H5READER_EXPECT_ML_SUCCESS=1;H5READER_EXPERIMENTAL_SHIELDING_ML_DEVICE=rocm;H5READER_EXPECT_ML_ACTIVE_DEVICE=rocm"
                LABELS "rest;e2e;ml;rocm"
                TIMEOUT 120
            )
        endif()
    endif()

    set(H5READER_CHIGNOLIN_FIXTURE_DEFAULT
        "C:/projects/reader-data/chignolin-oneoff/chignolin-july20260717.LGS"
        CACHE FILEPATH "Optional chignolin/CLN025 dense-DFT .LGS for REST fixture smoke")
    add_test(NAME h5reader_rest_chignolin_fixture
             COMMAND ${Python3_EXECUTABLE} -m pytest
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/rest/test_chignolin_fixture.py -v --tb=short)
    set_tests_properties(h5reader_rest_chignolin_fixture PROPERTIES
        ENVIRONMENT
            "H5READER_BINARY=$<TARGET_FILE:h5reader>;H5READER_REST_FIXTURE=${H5READER_CHIGNOLIN_FIXTURE_DEFAULT}"
            LABELS "rest;e2e;chignolin"
            TIMEOUT 180
    )

    set(H5READER_TRPCAGE_FIXTURE_DEFAULT
        "C:/projects/reader-data/trpcage-tc5b-1l2y/trpcage-july20260717.LGS"
        CACHE FILEPATH "Trp-cage dense-DFT .LGS for REST fixture smoke")
    add_test(NAME h5reader_rest_trpcage_fixture
             COMMAND ${Python3_EXECUTABLE} -m pytest
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/rest/test_trpcage_fixture.py -v --tb=short)
    set_tests_properties(h5reader_rest_trpcage_fixture PROPERTIES
        ENVIRONMENT
            "H5READER_BINARY=$<TARGET_FILE:h5reader>;H5READER_REST_FIXTURE=${H5READER_TRPCAGE_FIXTURE_DEFAULT}"
        LABELS "rest;e2e;trpcage"
        TIMEOUT 180
    )
endif()

# ============================================================================
# Install / package
# ============================================================================

install(TARGETS h5reader
    RUNTIME DESTINATION .
    COMPONENT Runtime
)

if(WIN32)
    set(H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR
        "ml/experimental_shielding_ml"
        CACHE STRING "Install subdirectory for the Experimental Shielding ML runtime")

    if(NOT H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT
            AND DEFINED ENV{H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT})
        file(TO_CMAKE_PATH "$ENV{H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}"
             H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT)
    endif()
    set(H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}"
        CACHE PATH
        "Directory containing the Experimental Shielding ML bundle")
    if(H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT
            AND NOT H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT)
        set(H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/rocm")
    endif()
    set(H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT}"
        CACHE PATH "Directory containing the native ROCm inference runtime")

    if(H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT)
        set(H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/model.ts"
            CACHE FILEPATH "TorchScript model for Experimental Shielding ML" FORCE)
        set(H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/manifest.json"
            CACHE FILEPATH "Provenance manifest for Experimental Shielding ML" FORCE)
        set(H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/infer.exe"
            CACHE FILEPATH "CPU inference helper for Experimental Shielding ML" FORCE)
        foreach(_torch_candidate IN ITEMS
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}"
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/torch/lib"
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/libtorch/lib")
            if(EXISTS "${_torch_candidate}/c10.dll")
                set(H5READER_LIBTORCH_CPU_BIN_DIR "${_torch_candidate}"
                    CACHE PATH "LibTorch CPU DLL directory for Experimental Shielding ML" FORCE)
                break()
            endif()
        endforeach()
    endif()

    set(H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL}"
        CACHE FILEPATH "TorchScript model for Experimental Shielding ML")
    set(H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST}"
        CACHE FILEPATH "Provenance manifest for Experimental Shielding ML")
    set(H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER}"
        CACHE FILEPATH "CPU inference helper for Experimental Shielding ML")
    set(H5READER_LIBTORCH_CPU_BIN_DIR
        "${H5READER_LIBTORCH_CPU_BIN_DIR}"
        CACHE PATH "LibTorch CPU DLL directory for Experimental Shielding ML")

    if(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML)
        foreach(_ml_asset IN ITEMS
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL}"
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST}"
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER}")
            if(NOT EXISTS "${_ml_asset}")
                message(FATAL_ERROR
                    "Experimental Shielding ML packaging is required, but this asset is missing: ${_ml_asset}\n"
                    "Set H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT to a portable bundle directory, "
                    "or set H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL, "
                    "H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST, "
                    "H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER, and "
                    "H5READER_LIBTORCH_CPU_BIN_DIR explicitly.")
            endif()
        endforeach()
        if(NOT EXISTS "${H5READER_LIBTORCH_CPU_BIN_DIR}")
            message(FATAL_ERROR
                "Experimental Shielding ML packaging is required, but LibTorch CPU DLL dir is missing: ${H5READER_LIBTORCH_CPU_BIN_DIR}\n"
                "Place the LibTorch CPU DLLs in H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT, "
                "H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT/torch/lib, or set "
                "H5READER_LIBTORCH_CPU_BIN_DIR explicitly.")
        endif()

        if(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML_ROCM)
            foreach(_rocm_asset IN ITEMS
                    infer.exe
                    c10.dll
                    c10_hip.dll
                    torch_cpu.dll
                    torch_hip.dll
                    caffe2_nvrtc.dll
                    aotriton_v2.dll
                    amdhip64_7.dll
                    amd_comgr0702.dll
                    hiprtc0702.dll
                    MIOpen.dll
                    rocblas.dll
                    libhipblaslt.dll
                    rocblas/library/TensileLibrary_lazy_gfx1151.dat
                    hipblaslt/library)
                if(NOT EXISTS
                        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT}/${_rocm_asset}")
                    message(FATAL_ERROR
                        "ROCm Experimental Shielding ML packaging is required, but this asset is missing: "
                        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT}/${_rocm_asset}")
                endif()
            endforeach()
        endif()

        set(_h5reader_libtorch_cpu_dlls)
        foreach(_torch_dll_name IN ITEMS
                c10.dll
                torch.dll
                torch_cpu.dll
                torch_global_deps.dll
                libiomp5md.dll
                libiompstubs5md.dll
                uv.dll)
            set(_torch_dll "${H5READER_LIBTORCH_CPU_BIN_DIR}/${_torch_dll_name}")
            if(NOT EXISTS "${_torch_dll}")
                message(FATAL_ERROR
                    "Experimental Shielding ML packaging is required, but LibTorch DLL is missing: ${_torch_dll}")
            endif()
            list(APPEND _h5reader_libtorch_cpu_dlls "${_torch_dll}")
        endforeach()

        set(_h5reader_experimental_shielding_ml_runtime_files
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL}"
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST}"
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER}"
            ${_h5reader_libtorch_cpu_dlls}
        )

        add_custom_target(h5reader_deploy_experimental_shielding_ml
            COMMAND "${CMAKE_COMMAND}" -E make_directory
                    "$<TARGET_FILE_DIR:h5reader>/${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}"
            COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                    ${_h5reader_experimental_shielding_ml_runtime_files}
                    "$<TARGET_FILE_DIR:h5reader>/${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}"
            DEPENDS ${_h5reader_experimental_shielding_ml_runtime_files}
            COMMENT "Deploying Experimental Shielding ML runtime next to h5reader"
            VERBATIM
        )
        add_dependencies(h5reader h5reader_deploy_experimental_shielding_ml)

        if(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML_ROCM)
            file(GLOB_RECURSE _h5reader_rocm_runtime_assets
                CONFIGURE_DEPENDS
                LIST_DIRECTORIES false
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT}/*")
            set(_h5reader_rocm_runtime_stamp
                "${CMAKE_CURRENT_BINARY_DIR}/h5reader_experimental_shielding_ml_rocm.stamp")
            add_custom_command(
                OUTPUT "${_h5reader_rocm_runtime_stamp}"
                COMMAND "${CMAKE_COMMAND}" -E make_directory
                        "$<TARGET_FILE_DIR:h5reader>/${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}/rocm"
                COMMAND "${CMAKE_COMMAND}" -E copy_directory
                        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT}"
                        "$<TARGET_FILE_DIR:h5reader>/${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}/rocm"
                COMMAND "${CMAKE_COMMAND}" -E touch "${_h5reader_rocm_runtime_stamp}"
                DEPENDS
                    ${_h5reader_rocm_runtime_assets}
                COMMENT "Deploying native ROCm Experimental Shielding ML runtime"
                VERBATIM)
            add_custom_target(h5reader_deploy_experimental_shielding_ml_rocm
                DEPENDS "${_h5reader_rocm_runtime_stamp}")
            add_dependencies(h5reader h5reader_deploy_experimental_shielding_ml_rocm)
        endif()

        install(FILES
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL}"
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST}"
            ${_h5reader_libtorch_cpu_dlls}
            DESTINATION "${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}"
            COMPONENT Runtime
        )
        install(PROGRAMS
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER}"
            DESTINATION "${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}"
            COMPONENT Runtime
        )
        if(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML_ROCM)
            install(DIRECTORY
                "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_ROOT}/"
                DESTINATION "${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}/rocm"
                COMPONENT Runtime)
        endif()
    endif()

    if(NOT H5READER_VCPKG_BIN_DIR)
        if(DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET
                AND EXISTS "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin")
            set(H5READER_VCPKG_BIN_DIR
                "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin"
                CACHE PATH "vcpkg runtime DLL directory")
        elseif(EXISTS "C:/Projects/vcpkg/installed/x64-windows/bin")
            set(H5READER_VCPKG_BIN_DIR
                "C:/Projects/vcpkg/installed/x64-windows/bin"
                CACHE PATH "vcpkg runtime DLL directory")
        endif()
    endif()

    if(NOT H5READER_MSVC_RUNTIME_DIR)
        set(_h5reader_msvc_redist_dirs)
        file(GLOB _h5reader_msvc_redist_roots
            LIST_DIRECTORIES true
            "C:/Program Files/Microsoft Visual Studio/*/*/VC/Redist/MSVC")
        foreach(_h5reader_msvc_redist_root IN LISTS _h5reader_msvc_redist_roots)
            if(EXISTS "${_h5reader_msvc_redist_root}")
                file(GLOB _h5reader_msvc_redist_candidates
                    LIST_DIRECTORIES true
                    "${_h5reader_msvc_redist_root}/*/x64/Microsoft.VC*.CRT")
                list(APPEND _h5reader_msvc_redist_dirs
                    ${_h5reader_msvc_redist_candidates})
            endif()
        endforeach()
        if(_h5reader_msvc_redist_dirs)
            list(SORT _h5reader_msvc_redist_dirs ORDER DESCENDING)
            list(GET _h5reader_msvc_redist_dirs 0 _h5reader_msvc_redist_dir)
            set(H5READER_MSVC_RUNTIME_DIR
                "${_h5reader_msvc_redist_dir}"
                CACHE PATH "MSVC runtime DLL directory")
        endif()
    endif()

    if(NOT H5READER_QT_DIR OR NOT EXISTS "${H5READER_QT_DIR}/bin/windeployqt.exe")
        message(FATAL_ERROR
            "Windows packaging requires H5READER_QT_DIR with bin/windeployqt.exe; "
            "set -DH5READER_QT_DIR=C:/Qt/6.10.2/msvc2022_64")
    endif()
    if(NOT H5READER_VTK_DIR OR NOT EXISTS "${H5READER_VTK_DIR}/bin")
        message(FATAL_ERROR
            "Windows packaging requires H5READER_VTK_DIR with a bin directory; "
            "set -DH5READER_VTK_DIR=C:/Projects/VTK")
    endif()
    if(NOT H5READER_VCPKG_BIN_DIR OR NOT EXISTS "${H5READER_VCPKG_BIN_DIR}")
        message(FATAL_ERROR
            "Windows packaging requires H5READER_VCPKG_BIN_DIR; set it to the vcpkg x64-windows bin directory")
    endif()
    if(NOT H5READER_MSVC_RUNTIME_DIR OR NOT EXISTS "${H5READER_MSVC_RUNTIME_DIR}")
        message(FATAL_ERROR
            "Windows packaging requires H5READER_MSVC_RUNTIME_DIR; set it to the MSVC x64 Microsoft.VC*.CRT directory")
    endif()

    # CTest is part of the supported Windows build, so tests must locate the
    # same Qt/VTK/vcpkg runtime DLLs without relying on a caller-maintained PATH.
    get_property(_h5reader_windows_tests DIRECTORY PROPERTY TESTS)
    foreach(_h5reader_test IN LISTS _h5reader_windows_tests)
        set_property(TEST "${_h5reader_test}" APPEND PROPERTY
            ENVIRONMENT_MODIFICATION
                "PATH=path_list_prepend:${H5READER_QT_DIR}/bin"
                "PATH=path_list_prepend:${H5READER_VTK_DIR}/bin"
                "PATH=path_list_prepend:${H5READER_VCPKG_BIN_DIR}")
    endforeach()

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/InstallWindowsRuntimeDeps.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/InstallWindowsRuntimeDeps.cmake"
        @ONLY
    )
    install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/InstallWindowsRuntimeDeps.cmake"
            COMPONENT Runtime)
endif()

set(CPACK_PACKAGE_NAME "h5reader")
set(CPACK_PACKAGE_VENDOR "h5reader")
set(CPACK_PACKAGE_VERSION "0.5.0")
if(WIN32)
    # Keep NSIS staging paths below legacy MAX_PATH: ROCm kernel filenames are long.
    set(CPACK_PACKAGE_FILE_NAME "h5reader-${CPACK_PACKAGE_VERSION}-win64")
else()
    set(CPACK_PACKAGE_FILE_NAME
        "h5reader-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif()
set(CPACK_PACKAGE_INSTALL_DIRECTORY "h5reader")
set(CPACK_MONOLITHIC_INSTALL ON)
set(CPACK_PACKAGE_EXECUTABLES "h5reader" "h5reader")

if(WIN32)
    set(CPACK_CREATE_DESKTOP_LINKS "h5reader")
    set(CPACK_NSIS_EXECUTABLES_DIRECTORY ".")
    find_program(H5READER_MAKENSIS_EXECUTABLE
        NAMES makensis
        HINTS
            "$ENV{ProgramFiles}/NSIS"
            "$ENV{ProgramFiles\(x86\)}/NSIS"
    )
    if(H5READER_MAKENSIS_EXECUTABLE)
        set(CPACK_GENERATOR "NSIS;ZIP")
        set(CPACK_NSIS_DISPLAY_NAME "h5reader")
        set(CPACK_NSIS_PACKAGE_NAME "h5reader")
        set(CPACK_NSIS_MODIFY_PATH OFF)
    else()
        set(CPACK_GENERATOR "ZIP")
    endif()
else()
    set(CPACK_GENERATOR "TGZ")
endif()
include(CPack)

message(STATUS "h5reader: Qt ${Qt6_VERSION}, VTK ${VTK_VERSION}, HDF5 ${HDF5_VERSION}")
