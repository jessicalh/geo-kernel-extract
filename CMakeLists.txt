cmake_minimum_required(VERSION 3.20)
project(nmr-shielding LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type: optimised with debug symbols.
# Override with: cmake -DCMAKE_BUILD_TYPE=Release ..
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    message(STATUS "Build type: RelWithDebInfo (default)")
endif()

# ============================================================================
# Dependencies
# ============================================================================

find_package(Eigen3 REQUIRED)
find_package(dssp QUIET)  # brings cifpp::cifpp with mm module
if(NOT dssp_FOUND)
    # Manual fallback for machines without dssp cmake config (Ubuntu 22.04).
    # Requires libdssp.a + dssp.hpp installed to /usr/local, cifpp from apt.
    find_library(DSSP_LIB dssp PATHS /usr/local/lib REQUIRED)
    find_path(DSSP_INCLUDE dssp.hpp PATHS /usr/local/include REQUIRED)
    find_package(cifpp REQUIRED)
    add_library(dssp::dssp STATIC IMPORTED)
    set_target_properties(dssp::dssp PROPERTIES
        IMPORTED_LOCATION "${DSSP_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${DSSP_INCLUDE}"
        INTERFACE_LINK_LIBRARIES "cifpp::cifpp"
    )
    message(STATUS "dssp: manual fallback — ${DSSP_LIB}")
endif()

# OpenBabel3 for bond perception and hybridisation
# Installed in miniforge3. The cmake config is there but does not
# export include dirs properly, so we set them explicitly.
set(OPENBABEL_INCLUDE "/usr/include/openbabel3"
    CACHE PATH "OpenBabel3 include directory")
set(OPENBABEL_LIB "/usr/lib/x86_64-linux-gnu/libopenbabel.so"
    CACHE FILEPATH "OpenBabel3 library")
if(NOT EXISTS ${OPENBABEL_LIB})
    message(FATAL_ERROR "OpenBabel3 not found at ${OPENBABEL_LIB}")
endif()
message(STATUS "OpenBabel3: ${OPENBABEL_LIB}")

# GROMACS 2026 for TPR reading (fleet loader)
# Linked against the installed libgromacs.so; internal headers from
# the source tree because GROMACS does not install them as public API.
set(GROMACS_LIB "/home/jessica/gromacs/build/gromacs-2026/lib/libgromacs_mpi.so"
    CACHE FILEPATH "GROMACS library")
set(GROMACS_SRC "/home/jessica/gromacs/src/gromacs-2026.0/src"
    CACHE PATH "GROMACS source tree (for internal headers)")
set(GROMACS_BUILD "/home/jessica/gromacs/build/gromacs-2026"
    CACHE PATH "GROMACS build tree (for generated headers)")
if(NOT EXISTS ${GROMACS_LIB})
    message(WARNING "libgromacs not found — fleet loader will not link")
endif()
message(STATUS "GROMACS: ${GROMACS_LIB}")

# Richardson lab reduce (hydrogen placement for protonation)
# Pre-built as static library from /home/jessica/builds/reduce-src
set(REDUCE_SRC "/home/jessica/builds/reduce-src"
    CACHE PATH "reduce source tree")
set(REDUCE_LIB "/home/jessica/builds/reduce-src/build/reduce_src/libreducelib.a"
    CACHE FILEPATH "reduce static library")
set(REDUCE_LIBPDB "/home/jessica/builds/reduce-src/build/libpdb/libpdb++.a"
    CACHE FILEPATH "reduce libpdb++ static library")
set(REDUCE_TOOLCLASSES "/home/jessica/builds/reduce-src/build/toolclasses/libtoolclasses.a"
    CACHE FILEPATH "reduce toolclasses static library")
set(REDUCE_HET_DICT "${REDUCE_SRC}/reduce_wwPDB_het_dict.txt"
    CACHE FILEPATH "reduce het dictionary")
if(NOT EXISTS ${REDUCE_LIB})
    message(WARNING "reducelib not found — protonation will not build")
endif()
message(STATUS "reduce: ${REDUCE_LIB}")

# MOPAC (PM7+MOZYME semiempirical) — linked library, replaces xTB subprocess
set(MOPAC_INCLUDE "/home/jessica/micromamba/envs/mm/include"
    CACHE PATH "MOPAC include directory (mopac.h)")
set(MOPAC_LIB "/home/jessica/micromamba/envs/mm/lib/libmopac.so"
    CACHE FILEPATH "MOPAC shared library")
if(NOT EXISTS ${MOPAC_LIB})
    message(WARNING "libmopac not found at ${MOPAC_LIB} — MopacResult will return nullptr")
endif()
message(STATUS "MOPAC: ${MOPAC_LIB}")

# GTest for unit tests
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# ============================================================================
# Header-only external libraries
# ============================================================================

add_library(extern_headers INTERFACE)
target_include_directories(extern_headers INTERFACE ${CMAKE_SOURCE_DIR}/extern)

# HighFive — header-only HDF5 C++ interface (extern/HighFive)
# Link HDF5 directly to avoid cmake pulling in curl/ssl transitive deps
add_library(highfive INTERFACE)
target_include_directories(highfive INTERFACE
    ${CMAKE_SOURCE_DIR}/extern/HighFive/include
    /usr/include/hdf5/serial)
target_link_libraries(highfive INTERFACE
    /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
    /usr/lib/x86_64-linux-gnu/libcurl.so
    /usr/lib/x86_64-linux-gnu/libz.so)

# ============================================================================
# xdrfile — XTC trajectory I/O (ported from fes-sampler, GROMACS XDR format)
# ============================================================================

add_library(xdrfile STATIC
    ${CMAKE_SOURCE_DIR}/extern/xdrfile/xdrfile.c
    ${CMAKE_SOURCE_DIR}/extern/xdrfile/xdrfile_xtc.c
)
target_include_directories(xdrfile PUBLIC ${CMAKE_SOURCE_DIR}/extern/xdrfile)

# ============================================================================
# libtorch (CUDA) — AIMNet2 neural network charge calculator
# ============================================================================

set(CMAKE_PREFIX_PATH
    "/home/jessica/.local/lib/python3.12/site-packages/torch/share/cmake"
    ${CMAKE_PREFIX_PATH})
find_package(Torch REQUIRED)
message(STATUS "Torch: ${TORCH_LIBRARIES}")
# libtorch headers use some constructs that trigger these warnings
# Suppress warnings from libtorch headers
set(TORCH_CXX_FLAGS -Wno-unused-parameter -Wno-sign-compare -Wno-deprecated-declarations)

# ============================================================================
# APBS C bridge — compiled as C, no C++, includes FETK headers
# ============================================================================

add_library(apbs_bridge STATIC
    ${CMAKE_SOURCE_DIR}/src/apbs_bridge.c
)
target_include_directories(apbs_bridge PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    /usr/include/apbs
    /usr/include/fetk
    /usr/include/maloc
)
# FETK libs (mc, punc) are separate on Ubuntu 24.04 but bundled
# inside libapbs on 22.04.  Link them only if found.
find_library(FETK_MC mc)
find_library(FETK_PUNC punc)
set(APBS_LIBS apbs_routines apbs_mg apbs_generic apbs_pmgc maloc)
if(FETK_MC)
    list(APPEND APBS_LIBS ${FETK_MC})
endif()
if(FETK_PUNC)
    list(APPEND APBS_LIBS ${FETK_PUNC})
endif()
target_link_libraries(apbs_bridge PRIVATE ${APBS_LIBS})

# ============================================================================
# Main library
# ============================================================================

add_library(nmr_shielding
    src/Types.cpp
    src/AminoAcidType.cpp
    src/Atom.cpp
    src/Ring.cpp
    src/Protein.cpp
    src/CovalentTopology.cpp
    src/ProteinConformation.cpp
    src/PdbFileReader.cpp
    src/GeometryResult.cpp
    src/DemoResult.cpp
    src/NamingRegistry.cpp
    src/IupacAtomIdentity.cpp
    src/DsspResult.cpp
    src/ChargeAssignmentResult.cpp
    src/ApbsFieldResult.cpp
    src/ProtonationDetectionResult.cpp
    src/ProtonationState.cpp
    src/PropkaProtonator.cpp
    src/KamlProtonator.cpp
    src/ChargeSource.cpp
    src/OrcaRunLoader.cpp
    src/GromacsEnsembleLoader.cpp
    src/SasaResult.cpp
    src/TrajectoryResult.cpp
    src/TrajectoryProtein.cpp
    src/Trajectory.cpp
    src/RunConfiguration.cpp
    src/Session.cpp
    src/BsWelfordTrajectoryResult.cpp
    src/PositionsTimeSeriesTrajectoryResult.cpp
    src/ChiRotamerSelectionTrajectoryResult.cpp
    src/GromacsFrameHandler.cpp
    src/GromacsEnergyResult.cpp
    src/FullSystemReader.cpp
    src/WaterFieldResult.cpp
    src/HydrationShellResult.cpp
    src/HydrationGeometryResult.cpp
    src/EeqResult.cpp
    src/OrcaShieldingResult.cpp
    src/MutationDeltaResult.cpp
    src/McConnellResult.cpp
    src/CoulombResult.cpp
    src/RingSusceptibilityResult.cpp
    src/HBondResult.cpp
    src/BiotSavartResult.cpp
    src/HaighMallionResult.cpp
    src/PiQuadrupoleResult.cpp
    src/DispersionResult.cpp
    src/KernelEvaluationFilter.cpp
    src/MopacResult.cpp
    src/MopacCoulombResult.cpp
    src/MopacMcConnellResult.cpp
    src/SpatialIndexResult.cpp
    src/EnrichmentResult.cpp
    src/MolecularGraphResult.cpp
    src/RuntimeEnvironment.cpp
    src/OperationLog.cpp
    src/OperationRunner.cpp
    src/ConformationResult.cpp
    src/GeometryChoice.cpp
    src/ReduceProtonation.cpp
    src/CalculatorConfig.cpp
    src/JobSpec.cpp
    src/AIMNet2Result.cpp
    src/BondedEnergyResult.cpp
)

target_include_directories(nmr_shielding PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/extern
    ${OPENBABEL_INCLUDE}
    ${MOPAC_INCLUDE}
    ${REDUCE_SRC}/reduce_src
    ${REDUCE_SRC}/libpdb
    ${REDUCE_SRC}/toolclasses
    ${GROMACS_SRC}
    ${GROMACS_SRC}/gromacs/fileio/include
    ${GROMACS_SRC}/gromacs/topology/include
    ${GROMACS_SRC}/gromacs/mdtypes
    ${GROMACS_SRC}/gromacs/utility/include
    ${GROMACS_SRC}/gromacs/math/include
    ${GROMACS_BUILD}/api/legacy/include
    ${GROMACS_SRC}/../api/legacy/include
    ${TORCH_INCLUDE_DIRS}
)

target_compile_options(nmr_shielding PRIVATE -Wall -Wextra -Wpedantic -march=native -flto ${TORCH_CXX_FLAGS})
target_compile_definitions(nmr_shielding PUBLIC
    NMR_DATA_DIR="${CMAKE_SOURCE_DIR}/data"
    HAS_OPENBABEL=1
    BOOLPREDEFINED=1
    CHARFUNCMACROS=1
    BRACKETOPERPARMS=1
    LISTFRIENDFIX=1
    INCTEMPLATEDEFNS=1
    LEFT_JUSTIFY_NUC_RES_OK=1
    AROMATICS_ACCEPT_HBONDS=1
    HET_DICTIONARY="${REDUCE_HET_DICT}"
)

target_link_libraries(nmr_shielding PUBLIC
    Eigen3::Eigen
    dssp::dssp
    cifpp::cifpp
    ${OPENBABEL_LIB}
    ${GROMACS_LIB}
    ${REDUCE_LIB}
    ${REDUCE_LIBPDB}
    ${REDUCE_TOOLCLASSES}
    ${MOPAC_LIB}
    apbs_bridge
    extern_headers
    highfive
    xdrfile
    ${TORCH_LIBRARIES}
)

# ============================================================================
# CLI executable
# ============================================================================

add_executable(nmr_extract src/nmr_extract.cpp)
target_link_libraries(nmr_extract PRIVATE nmr_shielding)


# ============================================================================
# Tests — tiered by cost and use case
#
#   unit_tests        <5s    No file I/O, pure logic
#   structure_tests   ~5min  Static PDB path (1UBQ, P84477), no MOPAC
#   trajectory_tests  ~5min  GROMACS fleet loader, trajectory path
#   mopac_tests       ~1hr   MOPAC semiempirical (test_mopac, mutation, pipeline)
#   batch_tests       ~2hr   Multi-protein sweep across consolidated/
#
# Separate executables (unchanged):
#   smoke_tests, fleet_smoke_tests, fes_fleet_smoke_tests,
#   job_spec_tests, gromacs_streaming_tests, water_field_tests
# ============================================================================

enable_testing()
include(GoogleTest)

# Helper: define a test executable with standard boilerplate
macro(nmr_test_target NAME)
    add_executable(${NAME}
        tests/test_main.cpp
        tests/TestEnvironment.cpp
        ${ARGN}
    )
    target_compile_definitions(${NAME} PRIVATE
        NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    )
    target_link_libraries(${NAME} PRIVATE nmr_shielding GTest::gtest)
    gtest_add_tests(TARGET ${NAME})
endmacro()

# ── Unit tests (<5s) — no protein loading, pure logic ──────────────────
nmr_test_target(unit_tests
    tests/test_spherical_tensor.cpp
    tests/test_amino_acid.cpp
    tests/test_ring_hierarchy.cpp
    tests/test_object_model.cpp
    tests/test_naming_registry.cpp
    tests/test_iupac_atom_identity.cpp
    tests/test_runtime_environment.cpp
    tests/test_calculator_config.cpp
)

# ── Structure tests (~5 min) — static PDB path, no MOPAC ──────────────
nmr_test_target(structure_tests
    tests/test_pdb_loading.cpp
    tests/test_geometry_result.cpp
    tests/test_demo_result.cpp
    tests/test_two_conformations.cpp
    tests/test_traversal_dump.cpp
    tests/test_dssp_result.cpp
    tests/test_protonation_detection.cpp
    tests/test_atom_flat.cpp
    tests/test_foundation_results.cpp
    tests/test_apbs_field_result.cpp
    tests/test_apbs_wired.cpp
    tests/test_apbs_ff14sb.cpp
    tests/test_protonation_pipeline.cpp
    tests/test_mcconnell_result.cpp
    tests/test_coulomb_result.cpp
    tests/test_ring_susceptibility_result.cpp
    tests/test_hbond_result.cpp
    tests/test_biot_savart_result.cpp
    tests/test_haigh_mallion_result.cpp
    tests/test_pi_quadrupole_result.cpp
    tests/test_dispersion_result.cpp
    tests/test_calculation_runner.cpp
    tests/test_pipeline_and_sample.cpp
    tests/test_write_features.cpp
)

# ── Trajectory tests (~5 min) — GROMACS fleet loader path ─────────────
nmr_test_target(trajectory_tests
    tests/test_fleet_loader.cpp
)

# ── MOPAC tests (~1 hr) — semiempirical charge calculation ────────────
nmr_test_target(mopac_tests
    tests/test_mopac_result.cpp
    tests/test_full_pipeline.cpp
    tests/test_mutation_delta.cpp
)

# ── Batch tests (~2 hr) — multi-protein sweep across consolidated ─────
nmr_test_target(batch_tests
    tests/test_batch_mcconnell.cpp
    tests/test_batch_coulomb_ringchi.cpp
    tests/test_batch_biot_savart_haigh_mallion.cpp
    tests/test_batch_piquad_disp.cpp
)

# ============================================================================
# Smoke tests — end-to-end pipeline validation (~5 min)
# ============================================================================

add_executable(smoke_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_smoke.cpp
)
target_compile_definitions(smoke_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(smoke_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_add_tests(TARGET smoke_tests)

# ============================================================================
# Fleet smoke test — separate exe (~40 min, 10 GROMACS poses x full pipeline)
# ============================================================================

add_executable(fleet_smoke_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_smoke_fleet.cpp
)
target_compile_definitions(fleet_smoke_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(fleet_smoke_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_add_tests(TARGET fleet_smoke_tests)

# ============================================================================
# JobSpec tests — separate exe (~10 min, not the 2.5-hour full suite)
# ============================================================================

add_executable(job_spec_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_job_spec.cpp
)
target_compile_definitions(job_spec_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(job_spec_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_add_tests(TARGET job_spec_tests)

# ============================================================================
# Gromacs streaming test — trajectory streaming proof of concept
# ============================================================================

add_executable(gromacs_streaming_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_gromacs_streaming.cpp
)
target_compile_definitions(gromacs_streaming_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(gromacs_streaming_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_add_tests(TARGET gromacs_streaming_tests)

# ============================================================================
# Water field test — explicit solvent calculators on full-system trajectory
# ============================================================================

add_executable(water_field_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_water_field.cpp
)
target_compile_definitions(water_field_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(water_field_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_add_tests(TARGET water_field_tests)

# ============================================================================
# FES fleet smoke test — fes-sampler naming, 2 proteins x 6 poses (~15 min)
# ============================================================================

add_executable(fes_fleet_smoke_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_smoke_fes_fleet.cpp
)
target_compile_definitions(fes_fleet_smoke_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(fes_fleet_smoke_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_add_tests(TARGET fes_fleet_smoke_tests)

# ============================================================================
# Example driver
# ============================================================================

add_executable(nmr_pass0_demo tests/pass0_demo.cpp)
target_link_libraries(nmr_pass0_demo PRIVATE nmr_shielding)

# ============================================================================
# Qt/VTK viewer (optional — only built if Qt6 and VTK 9 are found)
# ============================================================================

find_package(Qt6 COMPONENTS Widgets OpenGLWidgets Network QUIET)
find_package(VTK 9 COMPONENTS
    IOChemistry DomainsChemistry DomainsChemistryOpenGL2
    FiltersSources FiltersCore FiltersGeneral
    RenderingOpenGL2 RenderingCore RenderingAnnotation
    GUISupportQt InteractionStyle ImagingHybrid
    FiltersFlowPaths CommonColor
    QUIET)

if(Qt6_FOUND AND VTK_FOUND)
    set(CMAKE_AUTOMOC ON)
    add_executable(nmr-viewer
        ui/src/main_viewer.cpp
        ui/src/MainWindow.cpp
        ui/src/ComputeWorker.cpp
        ui/src/RestServer.cpp
        ui/src/BackboneRibbonOverlay.cpp
        ui/src/ButterflyOverlay.cpp
        ui/src/TensorGlyph.cpp
        ui/src/EllipsoidGlyph.cpp
        ui/src/FieldOverlay.cpp
        ui/src/FieldGridOverlay.cpp
        ui/src/IsosurfaceOverlay.cpp
        ui/src/RingCurrentOverlay.cpp
        ui/src/PeptideBondOverlay.cpp
        # Standalone analysis H5 reader (read-only consumer of the viewer).
        # Library stays Qt-free; Qt binding happens in ComputeWorker/MainWindow.
        fileformat/analysis_file.cpp
    )
    target_include_directories(nmr-viewer PRIVATE
        ${CMAKE_SOURCE_DIR}/ui/src
        ${CMAKE_SOURCE_DIR}/fileformat
    )
    # HighFive reads the HDF5 header version and picks between H5Dvlen_reclaim
    # (<1.12) and H5Treclaim (>=1.12).  OpenBabel's transitive micromamba
    # include path contributes HDF5 1.14 headers; the link library is the
    # system libhdf5 1.10 which lacks H5Treclaim.  Prepending system HDF5
    # here keeps HighFive emitting the old API matching the link library.
    target_include_directories(nmr-viewer BEFORE PRIVATE
        /usr/include/hdf5/serial
    )
    target_compile_definitions(nmr-viewer PRIVATE
        NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    )
    target_link_libraries(nmr-viewer PRIVATE
        nmr_shielding
        highfive
        Qt6::Widgets Qt6::OpenGLWidgets Qt6::Network
        ${VTK_LIBRARIES}
    )
    vtk_module_autoinit(TARGETS nmr-viewer MODULES ${VTK_LIBRARIES})
    message(STATUS "Qt6 + VTK found — building nmr-viewer")
else()
    message(STATUS "Qt6 or VTK not found — skipping nmr-viewer")
endif()
