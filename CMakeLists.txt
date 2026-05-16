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
# RUNPATH: prepend /usr/lib/x86_64-linux-gnu so the system libpq +
# libssl + libcrypto trio resolves consistently. Other linked libs
# (libmopac, libreducelib, libgromacs_mpi, libtorch, libnvrtc, ...)
# bring micromamba and torch lib dirs along as RUNPATH entries
# automatically; without this prepend the loader picks micromamba's
# libpq → micromamba's libssl (OpenSSL 3.2-3.6 symbol requirements)
# but system libcrypto (3.0 symbols only), and the load fails. We
# don't USE OpenSSL — libpq's NEEDED entries pull libssl/libcrypto
# at process startup regardless of runtime config — but the two
# halves of the SSL pair must come from the same install.
list(PREPEND CMAKE_BUILD_RPATH "/usr/lib/x86_64-linux-gnu")
list(PREPEND CMAKE_INSTALL_RPATH "/usr/lib/x86_64-linux-gnu")
set(CMAKE_BUILD_RPATH_USE_ORIGIN OFF)

# ============================================================================
# Dependencies
# ============================================================================

find_package(Eigen3 REQUIRED)
find_package(dssp QUIET)  # brings cifpp::cifpp with mm module

# PostgreSQL libpq — required for TripeptideDftTable (ProCS15 tripeptide
# DFT data lookup against the local tensorcs15 replica). The DSN is
# read from ~/.nmr_tools.toml [databases].tensorcs15. No DSN configured
# means TripeptideBackboneShieldingResult skips at Compute time;
# linking is unconditional so the calculator code is always present.
# PostgreSQL libpq — explicit system path. We do NOT use find_package
# here because find_package picks up micromamba's libpq, which pulls
# in a newer libssl (OpenSSL 3.2-3.6 symbols) that doesn't match the
# system's libcrypto (3.0 symbols). Loading the test binary fails. We
# use system libpq instead, which pulls system libssl + libcrypto
# consistently. We don't USE SSL at all (DSN is a Unix socket) but
# libpq's NEEDED entries pull libssl/libcrypto at process startup
# regardless of runtime config — so the versions must match.
set(PostgreSQL_LIBRARY "/usr/lib/x86_64-linux-gnu/libpq.so"
    CACHE FILEPATH "System libpq (avoids micromamba OpenSSL mismatch)")
set(PostgreSQL_INCLUDE_DIR "/usr/include/postgresql"
    CACHE PATH "System libpq headers")
if(NOT EXISTS ${PostgreSQL_LIBRARY})
    message(FATAL_ERROR "System libpq not found at ${PostgreSQL_LIBRARY}; "
                        "install libpq-dev or set PostgreSQL_LIBRARY explicitly")
endif()
if(NOT EXISTS ${PostgreSQL_INCLUDE_DIR}/libpq-fe.h)
    message(FATAL_ERROR "System libpq headers not found at "
                        "${PostgreSQL_INCLUDE_DIR}; install libpq-dev")
endif()
add_library(PostgreSQL::PostgreSQL UNKNOWN IMPORTED)
set_target_properties(PostgreSQL::PostgreSQL PROPERTIES
    IMPORTED_LOCATION "${PostgreSQL_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PostgreSQL_INCLUDE_DIR}")
message(STATUS "PostgreSQL libpq: ${PostgreSQL_LIBRARY} (explicit system path)")
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
    src/RingTopology.cpp
    src/Protein.cpp
    src/LegacyAmberTopology.cpp
    src/CovalentTopology.cpp
    src/ForceFieldChargeTable.cpp
    src/generated/LegacyAmberSemanticTables.cpp
    src/ProteinConformation.cpp
    src/PdbFileReader.cpp
    src/GeometryResult.cpp
    src/DemoResult.cpp
    src/NamingRegistry.cpp
    src/DsspResult.cpp
    src/ChargeAssignmentResult.cpp
    src/ApbsFieldResult.cpp
    src/ProtonationDetectionResult.cpp
    src/ProtonationState.cpp
    src/PropkaProtonator.cpp
    src/KamlProtonator.cpp
    src/ChargeSource.cpp
    src/AmberChargeResolver.cpp
    src/AmberLeapInput.cpp
    src/AmberPreparedChargeSource.cpp
    src/OrcaRunLoader.cpp
    src/SasaResult.cpp
    src/TrajectoryResult.cpp
    src/TrajectoryProtein.cpp
    src/Trajectory.cpp
    src/RunConfiguration.cpp
    src/Session.cpp
    src/BsWelfordTrajectoryResult.cpp
    src/BsShieldingTimeSeriesTrajectoryResult.cpp
    src/HmShieldingTimeSeriesTrajectoryResult.cpp
    src/TripeptideBackboneShieldingTimeSeriesTrajectoryResult.cpp
    src/TripeptideBackboneResidualVecTimeSeriesTrajectoryResult.cpp
    src/TripeptideNeighborShieldingTimeSeriesTrajectoryResult.cpp
    src/TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult.cpp
    src/TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult.cpp
    src/TripeptideBackboneMethodTagTimeSeriesTrajectoryResult.cpp
    src/LarsenHBondWaterTermTimeSeriesTrajectoryResult.cpp
    src/LarsenHBondCountTimeSeriesTrajectoryResult.cpp
    src/LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult.cpp
    src/LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult.cpp
    src/LarsenHBond1pHaBShieldingTimeSeriesTrajectoryResult.cpp
    src/LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult.cpp
    src/BsAnomalousAtomMarkerTrajectoryResult.cpp
    src/BsT0AutocorrelationTrajectoryResult.cpp
    src/BondLengthStatsTrajectoryResult.cpp
    src/PositionsTimeSeriesTrajectoryResult.cpp
    src/ChiRotamerSelectionTrajectoryResult.cpp
    src/GromacsFrameHandler.cpp
    src/GromacsEnergyResult.cpp
    src/GromacsFramePullResult.cpp
    src/GromacsToAmberReadbackBlock.cpp
    src/FramePdbEmitter.cpp
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
    src/CategoryInfoProjection.cpp
    src/TopologySidecar.cpp
    src/OperationLog.cpp
    src/OperationRunner.cpp
    src/ConformationResult.cpp
    src/GeometryChoice.cpp
    src/ReduceProtonation.cpp
    src/CalculatorConfig.cpp
    src/JobSpec.cpp
    src/AIMNet2Result.cpp
    src/AIMNet2PolarisabilityResult.cpp
    src/PlanarGeometryResult.cpp
    src/BondedEnergyResult.cpp
    src/TripeptideDftTable.cpp
    src/TripeptidePoseAssembler.cpp
    src/TripeptideBackboneShieldingResult.cpp
    src/TripeptideNeighborShieldingResult.cpp
    src/LarsenResidue.cpp
    src/LarsenHBondGrid.cpp
    src/LarsenHBondShieldingResult.cpp
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
    PostgreSQL::PostgreSQL
    apbs_bridge
    extern_headers
    highfive
    xdrfile
    ${TORCH_LIBRARIES}
)

# Match HDF5 1.10 headers to the 1.10 .so for every TU in the library —
# OpenBabel/MOPAC's micromamba include path contributes 1.14 headers
# which causes HighFive to emit H5Treclaim and bakes mismatched
# H5_VERS_* into compilation units; the link library is system
# libhdf5 1.10. Without this prepend, the library's WriteH5 paths
# (TrajectoryProtein, *TimeSeriesTrajectoryResult, etc.) compile against
# 1.14 and crash at runtime when HDF5 init detects the mismatch.
# Same hazard pattern as the gromacs_streaming_tests / nmr-viewer
# fixes below.
target_include_directories(nmr_shielding BEFORE PUBLIC
    /usr/include/hdf5/serial
)

# ============================================================================
# CLI executable
# ============================================================================

add_executable(nmr_extract src/nmr_extract.cpp)
target_link_libraries(nmr_extract PRIVATE nmr_shielding)
# Same HDF5-1.10-headers prepend for nmr_extract.cpp itself (the
# RunTrajectory WriteH5 site is in this TU).
target_include_directories(nmr_extract BEFORE PRIVATE
    /usr/include/hdf5/serial
)


# ============================================================================
# Tests — tiered by cost and use case
#
#   unit_tests        <5s    No file I/O, pure logic
#   structure_tests   ~5min  Static PDB path (1UBQ, P84477), no MOPAC
#   trajectory_tests  ~5min  AMBER trajectory loader (fleet_amber/)
#   mopac_tests       ~1hr   MOPAC semiempirical (test_mopac, mutation, pipeline)
#   batch_tests       ~2hr   Multi-protein sweep across consolidated/
#
# Separate executables (unchanged):
#   smoke_tests, job_spec_tests, gromacs_streaming_tests
#
# Retired 2026-05-04 (Phase 1 CHARMM-dead-letter cleanup, sources in
# tests/bones/): fleet_smoke_tests, fes_fleet_smoke_tests, water_field_tests.
# water_field_tests will be re-derived against fleet_amber/ in Phase 2.
# ============================================================================

enable_testing()
include(GoogleTest)

# Resolve PyTorch's bundled cu13 nvrtc-builtins directory at configure
# time. AIMNet2 / libtorch CUDA JIT dlopens libnvrtc-builtins.so.13.0
# from this directory and the dlopen chain does not consult the test
# binary's RPATH (only LD_LIBRARY_PATH). Every test target must run
# with this on LD_LIBRARY_PATH; this is enforced via the gtest_discover_tests
# PROPERTIES ENVIRONMENT below, applied uniformly so manual invocation
# (./build/<binary>) and ctest invocation share the same env source.
# See reference_nvrtc_rpath_fix and scripts/run_with_cuda_env.sh
# (the same path resolution, used for manual runs of the binaries).
execute_process(
    COMMAND python3 -c "import nvidia.cu13; print(list(nvidia.cu13.__path__)[0])"
    OUTPUT_VARIABLE NMR_NVIDIA_CU13_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _python_cu13_status
)
if(NOT _python_cu13_status EQUAL 0 OR NOT NMR_NVIDIA_CU13_DIR)
    message(FATAL_ERROR
        "Could not locate the bundled PyTorch cu13 directory via "
        "'python3 -c \"import nvidia.cu13\"'. Tests that exercise "
        "AIMNet2 will fail at runtime without LD_LIBRARY_PATH set "
        "to this directory's lib/. Install via "
        "'pip install --user nvidia-cu13' or equivalent. See "
        "memory entry reference_nvrtc_rpath_fix for the full picture.")
endif()
set(NMR_NVRTC_LIB_DIR "${NMR_NVIDIA_CU13_DIR}/lib")
message(STATUS "Tests: bundled cu13 lib dir = ${NMR_NVRTC_LIB_DIR}")

# Test environment string. Two prepends to LD_LIBRARY_PATH, in order:
#
#  1. /usr/lib/x86_64-linux-gnu so the system libpq + libssl + libcrypto
#     trio resolves consistently. The binary's RUNPATH includes
#     micromamba lib (for libmopac, libreducelib, etc.), and without
#     this prepend the loader picks micromamba's libpq → micromamba's
#     libssl (OpenSSL 3.2-3.6 symbol requirements) but system libcrypto
#     (3.0 symbols only) — load fails. We don't USE OpenSSL (DSN is a
#     Unix socket) but libpq's NEEDED entries pull libssl/libcrypto at
#     process startup regardless. Prepending /usr/lib forces both halves
#     of the SSL pair to the same install.
#
#  2. cu13 lib so libnvrtc.so.13's dlopen of libnvrtc-builtins.so.13.0
#     succeeds at AIMNet2 inference time.
#
# $ENV{LD_LIBRARY_PATH} captures the configure-time shell's existing
# path (orca, system cuda); we don't displace it.
set(NMR_TEST_ENVIRONMENT
    "LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${NMR_NVRTC_LIB_DIR}:$ENV{LD_LIBRARY_PATH}"
)

# Helper: define a test executable with standard boilerplate. Every
# test target gets the AIMNet2-compatible environment uniformly. No
# escape hatch — manual `./build/<binary>` invocation goes through
# scripts/run_with_cuda_env.sh; ctest invocation uses this property.
macro(nmr_test_target NAME)
    add_executable(${NAME}
        tests/test_main.cpp
        tests/TestEnvironment.cpp
        ${ARGN}
    )
    target_compile_definitions(${NAME} PRIVATE
        NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    )
    target_include_directories(${NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/tests
    )
    target_link_libraries(${NAME} PRIVATE nmr_shielding GTest::gtest)
    gtest_discover_tests(${NAME} PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")
endmacro()

# ── Unit tests (<5s) — no protein loading, pure logic ──────────────────
nmr_test_target(unit_tests
    tests/test_spherical_tensor.cpp
    tests/test_amino_acid.cpp
    tests/test_ring_hierarchy.cpp
    tests/test_object_model.cpp
    tests/test_naming_registry.cpp
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
    tests/test_amber_charge_resolver.cpp
    tests/test_amber_leap_input.cpp
    tests/test_amber_prepared_charge_source.cpp
    tests/test_apbs_field_result.cpp
    tests/test_apbs_wired.cpp
    tests/test_apbs_ff14sb.cpp
    tests/test_aimnet2_polarisability.cpp
    tests/test_planar_geometry_result.cpp
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
    tests/test_tripeptide_backbone_shielding.cpp
    tests/test_tripeptide_neighbor_shielding.cpp
    tests/test_larsen_residue_perception.cpp
    tests/test_larsen_residue_against_source_log.cpp
    tests/test_larsen_residue_ser_sidechain.cpp
    tests/test_larsen_residue_wl_ambiguity.cpp
    tests/test_larsen_hbond_grid.cpp
    tests/test_atom_semantic_predicates.cpp
    tests/test_larsen_hbond_shielding.cpp
)

# ── Trajectory tests (~5 min) — AMBER trajectory loader (fleet_amber/) ─
nmr_test_target(trajectory_tests
    tests/test_amber_trajectory.cpp
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
    tests/BlessCompare.cpp
    tests/test_smoke.cpp
)
target_compile_definitions(smoke_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(smoke_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_discover_tests(smoke_tests PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# String barrier discipline test (lightweight; greps source tree).
# Asserts that gemmi / RDKit / cifpp use stays bounded per the
# substrate-side string barrier in
# spec/plan/topology-substrate-implementation-plan-2026-05-05.md.
# ============================================================================

add_executable(string_barrier_tests
    tests/test_string_barrier.cpp
)
target_compile_definitions(string_barrier_tests PRIVATE
    NMR_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
target_link_libraries(string_barrier_tests PRIVATE
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(string_barrier_tests PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")

# ============================================================================
# Topology semantic substrate integration tests
#
# Loads representative AMBER-prepared protein PDBs (1ubq + fleet_amber 1Z9B,
# 1P9J), runs Protein::FinalizeConstruction, and audits the per-atom
# AtomSemanticTable substrate populated on LegacyAmberTopology. Spot-checks
# NTERM/CTERM cap composition, GLY HA2 prochiral, HIS variant chemistry,
# Pro res.H == NONE, and the Bundle B post-protonation re-canonicalisation
# of LYS-labelled-LYN-chemistry residues (1Z9B). Per
# spec/plan/topology-encoding-dependencies-2026-05-05.md §H.5.
# ============================================================================

add_executable(topology_semantic_integration_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/topology/test_legacy_amber_semantic_integration.cpp
)
target_compile_definitions(topology_semantic_integration_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(topology_semantic_integration_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/tests
)
target_link_libraries(topology_semantic_integration_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_discover_tests(topology_semantic_integration_tests
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")

# ============================================================================
# CategoryInfoProjection tests
#
# Loads 1UBQ via BuildFromProtonatedPdb, configures the projection from
# RuntimeEnvironment::BmrbAtomNom(), and asserts atom_nom.tbl-driven
# IUPAC/BMRB queries (Gly Hα inversion, methylene stereo, residue
# variant 3-letter mapping). NPY file lands with a structured-dtype
# header naming all expected fields.
# ============================================================================

add_executable(category_info_projection_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_category_info_projection.cpp
)
target_compile_definitions(category_info_projection_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(category_info_projection_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/tests
)
target_link_libraries(category_info_projection_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_discover_tests(category_info_projection_tests
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")

# ============================================================================
# TopologySidecar tests
#
# Verifies that bonds.npy, rings.npy, ring_membership.npy, and
# extraction_manifest.json land alongside atoms_category_info.npy with
# axis sizes matching the protein's typed topology surface, bond
# endpoints within the atom axis, and idempotent byte-identical
# re-emission of the topology files. See src/TopologySidecar.{h,cpp}.
# ============================================================================

add_executable(topology_sidecar_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_topology_sidecar.cpp
)
target_compile_definitions(topology_sidecar_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(topology_sidecar_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/tests
)
target_link_libraries(topology_sidecar_tests PRIVATE
    nmr_shielding
    GTest::gtest
)
gtest_discover_tests(topology_sidecar_tests
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")

# fleet_smoke_tests retired 2026-05-04 to tests/bones/test_smoke_fleet.cpp
# (CHARMM/XTC fleet/1A6J_5789 fixture).

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
gtest_discover_tests(job_spec_tests PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")

# ============================================================================
# Gromacs streaming test — trajectory streaming proof of concept
# ============================================================================

add_executable(gromacs_streaming_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_amber_streaming.cpp
    tests/test_frame_pdb_emitter.cpp
)
target_compile_definitions(gromacs_streaming_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
# Match HDF5 1.10 headers to the 1.10 .so — the same hazard documented
# on nmr-viewer below: OpenBabel's micromamba include path contributes
# HDF5 1.14 headers which emit H5Treclaim, but the link library is
# system libhdf5 1.10. The per-TR H5 round-trip test reads via HighFive
# directly from the test TU so this header ordering matters here.
target_include_directories(gromacs_streaming_tests BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(gromacs_streaming_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(gromacs_streaming_tests PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")

# ============================================================================
# HmShieldingTimeSeries tests
#
# Discipline trio + synthetic unit + integration for the FO-pattern
# dense-buffer TR that captures HaighMallionResult::hm_shielding_contribution
# into the /trajectory/hm_shielding_time_series/ H5 group.
# ============================================================================

add_executable(hm_shielding_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_hm_shielding_time_series.cpp
)
target_compile_definitions(hm_shielding_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(hm_shielding_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(hm_shielding_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(hm_shielding_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# TripeptideBackboneShieldingTimeSeries tests
#
# Discipline trio + synthetic unit + integration sanity for the FO-pattern
# dense-buffer TR that captures Larsen 2015 σ_BB^i into the
# /trajectory/tripeptide_bb_shielding_time_series/ H5 group.
# ============================================================================

add_executable(tripeptide_bb_shielding_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_tripeptide_bb_shielding_time_series.cpp
)
target_compile_definitions(tripeptide_bb_shielding_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(tripeptide_bb_shielding_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(tripeptide_bb_shielding_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(tripeptide_bb_shielding_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# TripeptideBackboneResidualVecTimeSeries tests
#
# Vec3-buffer TR; captures tripeptide_bb_residual_vec (post-Kabsch
# positional residual, polar 3-vector in protein lab frame) into
# /trajectory/tripeptide_bb_residual_vec_time_series/. Discipline trio
# + synthetic + log-overages integration.
# ============================================================================

add_executable(tripeptide_bb_residual_vec_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_tripeptide_bb_residual_vec_time_series.cpp
)
target_compile_definitions(tripeptide_bb_residual_vec_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(tripeptide_bb_residual_vec_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(tripeptide_bb_residual_vec_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(tripeptide_bb_residual_vec_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# TripeptideNeighborShieldingTimeSeries tests
#
# Same shape as the BB TR target above, against the Δσ_BB^{i±1} time
# series.
# ============================================================================

add_executable(tripeptide_neighbor_shielding_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_tripeptide_neighbor_shielding_time_series.cpp
)
target_compile_definitions(tripeptide_neighbor_shielding_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(tripeptide_neighbor_shielding_time_series_tests BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(tripeptide_neighbor_shielding_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(tripeptide_neighbor_shielding_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# TripeptideNeighborResidualVecPrev / Next TimeSeries tests
#
# Vec3-buffer TRs against the per-direction (i-1, i+1) neighbor residual
# fields. Source field is NaN-initialised per frame by the calculator;
# integration tests are NaN-tolerant and log magnitude overages.
# ============================================================================

add_executable(tripeptide_neighbor_residual_vec_prev_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_tripeptide_neighbor_residual_vec_prev_time_series.cpp
)
target_compile_definitions(tripeptide_neighbor_residual_vec_prev_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(tripeptide_neighbor_residual_vec_prev_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(tripeptide_neighbor_residual_vec_prev_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(tripeptide_neighbor_residual_vec_prev_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


add_executable(tripeptide_neighbor_residual_vec_next_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_tripeptide_neighbor_residual_vec_next_time_series.cpp
)
target_compile_definitions(tripeptide_neighbor_residual_vec_next_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(tripeptide_neighbor_residual_vec_next_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(tripeptide_neighbor_residual_vec_next_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(tripeptide_neighbor_residual_vec_next_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# TripeptideBackboneMethodTagTimeSeries tests
#
# int-buffer (uint8) (N, T) flat emission TR; captures the categorical
# tripeptide_bb_method_tag (0=no_match, 1=opbe, 2=pbe_ser). Establishes
# the int-buffer pattern for the bundle's other discrete-categorical TRs.
# ============================================================================

add_executable(tripeptide_bb_method_tag_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_tripeptide_bb_method_tag_time_series.cpp
)
target_compile_definitions(tripeptide_bb_method_tag_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(tripeptide_bb_method_tag_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(tripeptide_bb_method_tag_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(tripeptide_bb_method_tag_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# LarsenHBondWaterTermTimeSeries tests
#
# Scalar-double (N, T) flat TR; captures larsen_hbond_water_term (Larsen
# Δσ_w 2.07 ppm isotropic on HN atoms with no H-bond pairs).
# Establishes the scalar-buffer pattern for Larsen scalar/count TRs.
# ============================================================================

add_executable(larsen_hbond_water_term_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_larsen_hbond_water_term_time_series.cpp
)
target_compile_definitions(larsen_hbond_water_term_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(larsen_hbond_water_term_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(larsen_hbond_water_term_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(larsen_hbond_water_term_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# LarsenHBondCountTimeSeries tests
#
# int-buffer (N, T) flat TR; captures larsen_hbond_n_pairs (per-atom
# H-bond pair count summed across 4 Larsen classes).
# ============================================================================

add_executable(larsen_hbond_count_time_series_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_larsen_hbond_count_time_series.cpp
)
target_compile_definitions(larsen_hbond_count_time_series_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_include_directories(larsen_hbond_count_time_series_tests
    BEFORE PRIVATE
    /usr/include/hdf5/serial
)
target_link_libraries(larsen_hbond_count_time_series_tests PRIVATE
    nmr_shielding
    GTest::gtest
    highfive
)
gtest_discover_tests(larsen_hbond_count_time_series_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")


# ============================================================================
# Larsen per-class H-bond shielding TimeSeries tests (4 TRs)
#
# Four SphericalTensor TRs cloned from the BB shielding pattern, each
# against one of the four Larsen contribution classes
# (1°HB, 2°HB, 1°HαB, 2°HαB). Same (N, T, 9) flat-double emission.
# ============================================================================

foreach(_larsen_class IN ITEMS 1pHB 2pHB 1pHaB 2pHaB)
    set(_target larsen_hbond_${_larsen_class}_shielding_time_series_tests)
    add_executable(${_target}
        tests/test_main.cpp
        tests/TestEnvironment.cpp
        tests/test_larsen_hbond_${_larsen_class}_shielding_time_series.cpp
    )
    target_compile_definitions(${_target} PRIVATE
        NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    )
    target_include_directories(${_target} BEFORE PRIVATE
        /usr/include/hdf5/serial
    )
    target_link_libraries(${_target} PRIVATE
        nmr_shielding
        GTest::gtest
        highfive
    )
    gtest_discover_tests(${_target}
        DISCOVERY_MODE PRE_TEST
        PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")
endforeach()

# water_field_tests retired 2026-05-04 to tests/bones/test_water_field.cpp
# (CHARMM/XTC fleet_test_fullsys fixture + xtc_reader.h + GromacsEnsembleLoader).
# To be re-derived against fleet_amber/ in Phase 2.

# fes_fleet_smoke_tests retired 2026-05-04 to tests/bones/test_smoke_fes_fleet.cpp
# (fes-sampler retired; fixture in tests/bones/data/fes_fleet/).

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


# ============================================================================
# Build-time generator for LegacyAmberTopology semantic tables (opt-in).
#
# Links cifpp + RDKit. Reads CCD, runs perception, emits the typed-enum
# table at src/generated/LegacyAmberSemanticTables.cpp. The runtime
# library does NOT link RDKit -- the string barrier is enforced by
# the linker. See tools/topology/README.md and
# spec/plan/topology-substrate-implementation-plan-2026-05-05.md.
# ============================================================================

option(NMR_BUILD_TABLE_GENERATOR
    "Build the LegacyAmberTopology semantic tables generator (requires RDKit)"
    OFF)
if(NMR_BUILD_TABLE_GENERATOR)
    message(STATUS "Building topology table generator (NMR_BUILD_TABLE_GENERATOR=ON)")
    add_subdirectory(tools/topology)
endif()
