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
    src/io/CalcsetManifest.h
    src/io/CalcsetManifest.cpp
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

)

target_include_directories(h5reader_core PUBLIC
    src
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


# ============================================================================
# Rediscover — headless feature extractor (experimental branch
# h5-reader-pysr-spike). QCoreApplication only: NO widgets, NO VTK. Links
# h5reader_core + Eigen + the vendored nanoflann header. Loads one 1P9J
# calcset into an all-frames-resident RunData and writes the ring-current
# and McConnell feature/target CSVs. Not wired into the GUI; one-shot study
# substrate. See src/rediscover/{GUIDANCE,DESIGN}.md.
# ============================================================================

qt_add_executable(h5reader_extract
    src/rediscover/SphericalBasis.h
    src/rediscover/SphericalBasis.cpp
    src/rediscover/LocalFrameBasis.h
    src/rediscover/LocalFrameBasis.cpp
    src/rediscover/RediscoverTypes.h
    src/rediscover/TemporalIndex.h
    src/rediscover/TemporalIndex.cpp
    src/rediscover/TypedAtomIndex.h
    src/rediscover/TypedAtomIndex.cpp
    src/rediscover/SpatialIndexSet.h
    src/rediscover/SpatialIndexSet.cpp
    src/rediscover/RingGeometryCache.h
    src/rediscover/RingGeometryCache.cpp
    src/rediscover/ResidentIndexes.h
    src/rediscover/ResidentIndexes.cpp
    src/rediscover/Catalog.h
    src/rediscover/Catalog.cpp
    src/rediscover/AnalysisBody.h
    src/rediscover/ChargeStore.h
    src/rediscover/ChargeStore.cpp
    src/rediscover/OutputManifest.h
    src/rediscover/OutputManifest.cpp
    src/rediscover/RunData.h
    src/rediscover/RunData.cpp
    src/rediscover/RecordSink.h
    src/rediscover/RecordSink.cpp
    src/rediscover/ExtractionSupport.h
    src/rediscover/ExtractionSupport.cpp
    src/rediscover/EfgFeatureSink.h
    src/rediscover/EfgFeatureSink.cpp
    src/rediscover/EfgFeature.h
    src/rediscover/EfgFeature.cpp
    src/rediscover/Aimnet2FeatureSink.h
    src/rediscover/Aimnet2FeatureSink.cpp
    src/rediscover/Aimnet2Feature.h
    src/rediscover/Aimnet2Feature.cpp
    src/rediscover/BuckinghamEfieldSink.h
    src/rediscover/BuckinghamEfieldSink.cpp
    src/rediscover/BuckinghamEfield.h
    src/rediscover/BuckinghamEfield.cpp
    src/rediscover/RediscoveryExtraction.h
    src/rediscover/FrameSpatialIndex.h
    src/rediscover/FrameSpatialIndex.cpp
    src/rediscover/RingCurrentKernel.h
    src/rediscover/RingCurrentKernel.cpp
    src/rediscover/McConnellLiteratureKernel.h
    src/rediscover/McConnellLiteratureKernel.cpp
    src/rediscover/RingCurrentNeighborhood.h
    src/rediscover/RingCurrentNeighborhood.cpp
    src/rediscover/McConnellNeighborhood.h
    src/rediscover/McConnellNeighborhood.cpp
    src/rediscover/ChargeDipoleNeighborhood.h
    src/rediscover/ChargeDipoleNeighborhood.cpp
    # Functional API (SURFACE_DESIGN Layers 1-3): primitive verbs, curried
    # relationship combinators, the pure iterated-closure engine, and the two
    # composed relationships (ring_current / mcconnell) the driver runs.
    src/rediscover/Verbs.h
    src/rediscover/Verbs.cpp
    src/rediscover/Relationship.h
    src/rediscover/Relationship.cpp
    src/rediscover/RelationshipEngine.h
    src/rediscover/RelationshipEngine.cpp
    src/rediscover/ComposedRelationships.h
    src/rediscover/ComposedRelationships.cpp
    # broad_backbone — the first composed relationship over EVERY backbone atom
    # × all mechanisms (rings + aniso bonds + charge field), with its own
    # heterogeneous two-kind carrier (the target-repeat fix). Backbone frames
    # added in LocalFrameBasis.{h,cpp} above.
    src/rediscover/BroadBackboneSink.h
    src/rediscover/BroadBackboneSink.cpp
    src/rediscover/BroadBackbone.h
    src/rediscover/BroadBackbone.cpp
    src/rediscover/AllAtomEquivariantSink.h
    src/rediscover/AllAtomEquivariantSink.cpp
    src/rediscover/AllAtomEquivariant.h
    src/rediscover/AllAtomEquivariant.cpp
    src/rediscover/main_extract.cpp
)
target_include_directories(h5reader_extract PRIVATE
    src
    ${CMAKE_CURRENT_SOURCE_DIR}/extern
    ${HIGHFIVE_INCLUDE_DIR}
)
target_link_libraries(h5reader_extract PRIVATE
    h5reader_core
    Eigen3::Eigen
)
h5reader_apply_platform_target_settings(h5reader_extract)
if(COMMAND h5reader_apply_build_type_settings)
    h5reader_apply_build_type_settings(h5reader_extract)
endif()

include(CTest)
if(BUILD_TESTING AND Qt6Test_FOUND)
    # SphericalBasis library-order decomposition fixture (the foundation).
    qt_add_executable(h5reader_rediscover_tests
        tests/rediscover_basis_tests.cpp
        src/rediscover/SphericalBasis.h
        src/rediscover/SphericalBasis.cpp
    )
    target_include_directories(h5reader_rediscover_tests PRIVATE src)
    target_link_libraries(h5reader_rediscover_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_rediscover_tests COMMAND h5reader_rediscover_tests)

    # Broad-backbone local-frame fixtures (broad-backbone spike, 2026-06-01).
    # Pure synthetic Vec3 geometry over LocalFrameBasis — no calcset / H5 /
    # model. Pins the N / Cα / C=O / HA frame conventions + handedness +
    # degenerate-fail-loud by assertion (the aromatic-H positional-anchor bug
    # class). Links only LocalFrameBasis + Eigen + Qt::Test, like the basis test.
    qt_add_executable(h5reader_rediscover_backbone_frame_tests
        tests/rediscover_backbone_frame_tests.cpp
        src/rediscover/LocalFrameBasis.h
        src/rediscover/LocalFrameBasis.cpp
        src/rediscover/SphericalBasis.h
        src/rediscover/SphericalBasis.cpp
    )
    target_include_directories(h5reader_rediscover_backbone_frame_tests PRIVATE src)
    target_link_libraries(h5reader_rediscover_backbone_frame_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_rediscover_backbone_frame_tests
             COMMAND h5reader_rediscover_backbone_frame_tests)

    # APBS EFG per_atom_feature fixtures: the EFG T2 basis equals the DFT
    # target basis (DecomposeLibrary), rotation preserves |T2|, and the sink
    # writes finite row-aligned (rows,5) feature/target sidecars.
    qt_add_executable(h5reader_rediscover_efg_tests
        tests/rediscover_efg_tests.cpp
        src/rediscover/EfgFeatureSink.h
        src/rediscover/EfgFeatureSink.cpp
        src/rediscover/SphericalBasis.h
        src/rediscover/SphericalBasis.cpp
        src/io/QtNpyReader.h
        src/io/QtNpyReader.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_rediscover_efg_tests PRIVATE src)
    target_link_libraries(h5reader_rediscover_efg_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_rediscover_efg_tests
             COMMAND h5reader_rediscover_efg_tests)

    # APBS E-field per_atom_feature fixtures: the spine-side local-frame
    # projection is E_local.z, |E| is invariant, target T2 keeps the √6
    # library fixture, and T1 is emitted as an unverified audit payload only.
    qt_add_executable(h5reader_rediscover_buckingham_tests
        tests/rediscover_buckingham_tests.cpp
        src/rediscover/BuckinghamEfieldSink.h
        src/rediscover/BuckinghamEfieldSink.cpp
        src/rediscover/LocalFrameBasis.h
        src/rediscover/LocalFrameBasis.cpp
        src/rediscover/SphericalBasis.h
        src/rediscover/SphericalBasis.cpp
        src/io/QtNpyReader.h
        src/io/QtNpyReader.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_rediscover_buckingham_tests PRIVATE src)
    target_link_libraries(h5reader_rediscover_buckingham_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_rediscover_buckingham_tests
             COMMAND h5reader_rediscover_buckingham_tests)

    qt_add_executable(h5reader_rediscover_all_atom_equivariant_tests
        tests/rediscover_all_atom_equivariant_tests.cpp
        src/rediscover/AllAtomEquivariantSink.h
        src/rediscover/AllAtomEquivariantSink.cpp
        src/io/QtNpyReader.h
        src/io/QtNpyReader.cpp
        src/diagnostics/ErrorBus.h
        src/diagnostics/ErrorBus.cpp
    )
    target_include_directories(h5reader_rediscover_all_atom_equivariant_tests PRIVATE src)
    target_link_libraries(h5reader_rediscover_all_atom_equivariant_tests PRIVATE
        Qt6::Core
        Qt6::Test
        Eigen3::Eigen
    )
    add_test(NAME h5reader_rediscover_all_atom_equivariant_tests
             COMMAND h5reader_rediscover_all_atom_equivariant_tests)

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
