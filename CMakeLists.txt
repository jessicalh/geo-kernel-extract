cmake_minimum_required(VERSION 3.20)
project(nmr-shielding LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Always emit compile_commands.json so clang-tidy / IWYU / scan-build
# / clangd can find translation-unit-specific flags. Lives in the
# build directory; scripts/lint/* + clangd consume it from there.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
    message(FATAL_ERROR "libgromacs not found at ${GROMACS_LIB} — required "
                        "(the library links it unconditionally for trajectory "
                        "reading). Set -DGROMACS_LIB=... (and GROMACS_SRC / "
                        "GROMACS_BUILD), e.g. in CMakeUserPresets.json.")
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
# reduce links three archives unconditionally (reducelib + libpdb++ +
# toolclasses, see target_link_libraries below) — all required.
foreach(_reduce_archive REDUCE_LIB REDUCE_LIBPDB REDUCE_TOOLCLASSES)
    if(NOT EXISTS ${${_reduce_archive}})
        message(FATAL_ERROR "reduce archive ${_reduce_archive}=${${_reduce_archive}} "
                            "not found — required (linked unconditionally for "
                            "protonation). Set -D${_reduce_archive}=... (and "
                            "REDUCE_SRC), e.g. in CMakeUserPresets.json.")
    endif()
endforeach()
message(STATUS "reduce: ${REDUCE_LIB}")

# MOPAC (PM7+MOZYME semiempirical) — linked library, replaces xTB subprocess
set(MOPAC_INCLUDE "/home/jessica/micromamba/envs/mm/include"
    CACHE PATH "MOPAC include directory (mopac.h)")
set(MOPAC_LIB "/home/jessica/micromamba/envs/mm/lib/libmopac.so"
    CACHE FILEPATH "MOPAC shared library")
if(NOT EXISTS ${MOPAC_LIB})
    message(FATAL_ERROR "libmopac not found at ${MOPAC_LIB} — required "
                        "(the library links it unconditionally; the --mopac "
                        "toggle gates whether MOPAC RUNS, not whether it links). "
                        "Set -DMOPAC_LIB=... (and MOPAC_INCLUDE), e.g. in "
                        "CMakeUserPresets.json.")
endif()
message(STATUS "MOPAC: ${MOPAC_LIB}")

# GTest for unit tests. Prefer a system/installed GTest; otherwise fetch
# it from GitHub. The fetch needs network at first configure, so it is
# announced (not silent) and gated by NMR_FETCH_GTEST: default ON keeps
# the current behaviour, OFF requires a system GTest (for offline /
# restricted-network configures).
option(NMR_FETCH_GTEST "Fetch googletest from GitHub when no system GTest is found" ON)
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    if(NOT NMR_FETCH_GTEST)
        message(FATAL_ERROR "GTest not found and NMR_FETCH_GTEST=OFF — install "
                            "a system GTest (e.g. libgtest-dev) or re-enable "
                            "-DNMR_FETCH_GTEST=ON to fetch it from GitHub.")
    endif()
    message(STATUS "GTest not found — fetching googletest v1.14.0 from GitHub "
                   "(needs network; set -DNMR_FETCH_GTEST=OFF to require a "
                   "system GTest instead).")
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
    src/Cli/Parse.cpp
    src/Cli/Validate.cpp
    src/Cli/PrintUsage.cpp
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
    src/HmWelfordTrajectoryResult.cpp
    src/McConnellWelfordTrajectoryResult.cpp
    src/EeqWelfordTrajectoryResult.cpp
    src/SasaWelfordTrajectoryResult.cpp
    src/HBondCountWelfordTrajectoryResult.cpp
    src/BsShieldingTimeSeriesTrajectoryResult.cpp
    src/HmShieldingTimeSeriesTrajectoryResult.cpp
    src/McConnellShieldingTimeSeriesTrajectoryResult.cpp
    src/PiQuadrupoleShieldingTimeSeriesTrajectoryResult.cpp
    src/RingSusceptibilityShieldingTimeSeriesTrajectoryResult.cpp
    src/DispersionShieldingTimeSeriesTrajectoryResult.cpp
    src/HBondShieldingTimeSeriesTrajectoryResult.cpp
    src/SasaTimeSeriesTrajectoryResult.cpp
    src/AIMNet2ChargeTimeSeriesTrajectoryResult.cpp
    src/AIMNet2EmbeddingTimeSeriesTrajectoryResult.cpp
    src/AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult.cpp
    src/AIMNet2ChargeResponseGradientWelfordTrajectoryResult.cpp
    src/ApbsEfieldTimeSeriesTrajectoryResult.cpp
    src/ApbsEfgTimeSeriesTrajectoryResult.cpp
    src/MopacChargeWelfordTrajectoryResult.cpp
    src/MopacBondOrderWelfordTrajectoryResult.cpp
    src/MopacCoulombShieldingTimeSeriesTrajectoryResult.cpp
    src/MopacMcConnellShieldingTimeSeriesTrajectoryResult.cpp
    src/MopacVsFf14SbReconciliationTrajectoryResult.cpp
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
    src/GromacsEnergyTimeSeriesTrajectoryResult.cpp
    src/BondedEnergyTimeSeriesTrajectoryResult.cpp
    src/GromacsFramePullResult.cpp
    src/GromacsToAmberReadbackBlock.cpp
    src/FrameNpyEmitter.cpp
    src/FramePdbEmitter.cpp
    src/FullSystemReader.cpp
    src/WaterFieldResult.cpp
    src/WaterFieldTimeSeriesTrajectoryResult.cpp
    src/WaterFieldWelfordTrajectoryResult.cpp
    src/HydrationShellResult.cpp
    src/HydrationGeometryResult.cpp
    src/HydrationGeometryTimeSeriesTrajectoryResult.cpp
    src/HydrationGeometryWelfordTrajectoryResult.cpp
    src/HydrationShellTimeSeriesTrajectoryResult.cpp
    src/HydrationShellWelfordTrajectoryResult.cpp
    src/DihedralTimeSeriesTrajectoryResult.cpp
    src/DihedralBinTransitionTrajectoryResult.cpp
    src/Dssp8TimeSeriesTrajectoryResult.cpp
    src/Dssp8TransitionTrajectoryResult.cpp
    src/RingPuckerTimeSeriesTrajectoryResult.cpp
    src/JCouplingTimeSeriesTrajectoryResult.cpp
    src/RingNeighbourhoodTrajectoryStats.cpp
    src/RmsdTrackingTrajectoryResult.cpp
    src/RmsdSpikeSelectionTrajectoryResult.cpp
    src/DftPoseCoordinatorTrajectoryResult.cpp
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
    src/AIMNet2Result.cpp
    src/AIMNet2ChargeResponseGradientResult.cpp
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
# Tests
#
# Seven binaries, organised by scope:
#
#   unit_tests                 pure logic, no protein loading (<5s)
#   conformation_scope_tests   per-Conf-Result + per-calc + substrate + cli
#   trajectory_scope_tests     every TR + trajectory loader + per-TR discipline
#   mopac_tests                MOPAC ConfResult + MOPAC TR family
#                              (separate because fixtures are slow / optional)
#   batch_tests                multi-protein sweeps (~2 hr; separate to keep
#                              day-to-day ctest fast)
#   smoke_tests                end-to-end bless-baseline + smoke
#   string_barrier_tests       discipline grep on the source tree
#                              (no nmr_shielding link; hand-rolled below)
#
# One macro `nmr_test()` defines the first six. Discovery in PRE_TEST mode
# (AIMNet2/torch needs LD_LIBRARY_PATH set at test run time, before binary
# enumeration). Only `fast` is a sanctioned LABEL; hardware / dependency
# gating is NOT a test classifier (see memory entry
# `feedback_aimnet2_required_no_weasel`).
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
#     this prepend the loader picks micromamba's libpq -> micromamba's
#     libssl (OpenSSL 3.2-3.6 symbol requirements) but system libcrypto
#     (3.0 symbols only) -- load fails. We don't USE OpenSSL (DSN is a
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

# ----------------------------------------------------------------------------
# nmr_test(<NAME>
#     SOURCES <file.cpp> ...    test sources beyond the universal two
#     [LABELS fast]             only "fast" is sanctioned (absence = full)
#     [LINK   <lib> ...]        extras beyond nmr_shielding + GTest::gtest
# )
#
# Every nmr_test binary:
#   - builds from tests/test_main.cpp + tests/TestEnvironment.cpp + SOURCES;
#   - gets NMR_TEST_DATA_DIR + tests/ on the include path;
#   - links nmr_shielding + GTest::gtest + LINK;
#   - is discovered via gtest_discover_tests in DISCOVERY_MODE PRE_TEST with
#     ENVIRONMENT = NMR_TEST_ENVIRONMENT (the libssl + cu13 LD_LIBRARY_PATH
#     that AIMNet2's torch runtime needs at test invocation).
#
# When LINK contains "highfive", /usr/include/hdf5/serial is prepended to
# the include path so HighFive picks up the 1.10 ABI matching the system
# libhdf5 link (OpenBabel's micromamba contributes 1.14 headers that
# would emit symbols the link library doesn't have).
# ----------------------------------------------------------------------------
macro(nmr_test NAME)
    set(_multi SOURCES LABELS LINK)
    cmake_parse_arguments(_NT "" "" "${_multi}" ${ARGN})

    add_executable(${NAME}
        tests/test_main.cpp
        tests/TestEnvironment.cpp
        ${_NT_SOURCES}
    )
    target_compile_definitions(${NAME} PRIVATE
        NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    )

    list(FIND _NT_LINK "highfive" _hf_idx)
    if(NOT _hf_idx EQUAL -1)
        target_include_directories(${NAME} BEFORE PRIVATE /usr/include/hdf5/serial)
    endif()
    target_include_directories(${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/tests)

    target_link_libraries(${NAME} PRIVATE
        nmr_shielding
        GTest::gtest
        ${_NT_LINK}
    )

    set(_props ENVIRONMENT "${NMR_TEST_ENVIRONMENT}")
    if(_NT_LABELS)
        list(APPEND _props LABELS ${_NT_LABELS})
    endif()
    gtest_discover_tests(${NAME}
        DISCOVERY_MODE PRE_TEST
        PROPERTIES ${_props}
    )
endmacro()


# -- unit_tests -------------------------------------------------------------
# Pure logic; no protein loading, no torch init. <5s. LABELS fast.
nmr_test(unit_tests
    LABELS fast
    SOURCES
        tests/test_amino_acid.cpp
        tests/test_calculator_config.cpp
        tests/test_naming_registry.cpp
        tests/test_object_model.cpp
        tests/test_ring_hierarchy.cpp
        tests/test_runtime_environment.cpp
        tests/test_spherical_tensor.cpp
)


# -- conformation_scope_tests ----------------------------------------------
# Every per-Conf-Result and per-calculator unit test; substrate / topology /
# CategoryInfo projection; cli parse. Static PDB path, no MOPAC, no MD
# trajectory. Loads torch via nmr_shielding (AIMNet2 path).
nmr_test(conformation_scope_tests
    SOURCES
        tests/test_aimnet2_charge_response_gradient.cpp
        tests/test_amber_charge_resolver.cpp
        tests/test_amber_leap_input.cpp
        tests/test_amber_prepared_charge_source.cpp
        tests/test_apbs_ff14sb.cpp
        tests/test_apbs_field_result.cpp
        tests/test_apbs_wired.cpp
        tests/test_atom_flat.cpp
        tests/test_atom_semantic_predicates.cpp
        tests/test_biot_savart_result.cpp
        tests/test_calculation_runner.cpp
        tests/test_category_info_projection.cpp
        tests/test_cli_parse.cpp
        tests/test_coulomb_result.cpp
        tests/test_demo_result.cpp
        tests/test_dispersion_result.cpp
        tests/test_dssp_result.cpp
        tests/test_foundation_results.cpp
        tests/test_geometry_result.cpp
        tests/test_haigh_mallion_result.cpp
        tests/test_hbond_result.cpp
        tests/test_larsen_hbond_grid.cpp
        tests/test_larsen_hbond_shielding.cpp
        tests/test_larsen_residue_against_source_log.cpp
        tests/test_larsen_residue_perception.cpp
        tests/test_larsen_residue_ser_sidechain.cpp
        tests/test_larsen_residue_wl_ambiguity.cpp
        tests/test_mcconnell_result.cpp
        tests/test_pdb_loading.cpp
        tests/test_pi_quadrupole_result.cpp
        tests/test_pipeline_and_sample.cpp
        tests/test_planar_geometry_result.cpp
        tests/test_protonation_detection.cpp
        tests/test_ring_susceptibility_result.cpp
        tests/test_topology_sidecar.cpp
        tests/test_traversal_dump.cpp
        tests/test_tripeptide_backbone_shielding.cpp
        tests/test_tripeptide_neighbor_shielding.cpp
        tests/test_two_conformations.cpp
        tests/test_write_features.cpp
        tests/topology/test_legacy_amber_semantic_integration.cpp
)


# -- trajectory_scope_tests ------------------------------------------------
# Every TrajectoryResult test + trajectory loader (AMBER + Gromacs streaming)
# + FramePdbEmitter. Per-TR discipline trio + synthetic + integration tests
# all live here. Links highfive for H5RoundTrip checks. MOPAC-dependent TRs
# live in mopac_tests; batch sweeps live in batch_tests.
nmr_test(trajectory_scope_tests
    LINK highfive
    SOURCES
        tests/TestConfig.cpp
        tests/test_aimnet2_charge_time_series.cpp
        tests/test_aimnet2_embedding_charge_response_gradient_time_series.cpp
        tests/test_amber_streaming.cpp
        tests/test_amber_trajectory.cpp
        tests/test_apbs_efg_time_series.cpp
        tests/test_apbs_efield_time_series.cpp
        tests/test_bonded_energy_time_series.cpp
        tests/test_bs_anomalous_atom_marker.cpp
        tests/test_bs_welford.cpp
        tests/test_dihedral_bin_transition.cpp
        tests/test_dihedral_time_series.cpp
        tests/test_disp_shielding_time_series.cpp
        tests/test_dssp8_time_series.cpp
        tests/test_dssp8_transition.cpp
        tests/test_eeq_welford.cpp
        tests/test_frame_pdb_emitter.cpp
        tests/test_gromacs_energy_time_series.cpp
        tests/test_hbond_count_welford.cpp
        tests/test_hbond_shielding_time_series.cpp
        tests/test_hm_shielding_time_series.cpp
        tests/test_hm_welford.cpp
        tests/test_hydration_geometry_time_series.cpp
        tests/test_hydration_geometry_welford.cpp
        tests/test_hydration_shell_time_series.cpp
        tests/test_hydration_shell_welford.cpp
        tests/test_j_coupling_time_series.cpp
        tests/test_larsen_hbond_1pHB_shielding_time_series.cpp
        tests/test_larsen_hbond_1pHaB_shielding_time_series.cpp
        tests/test_larsen_hbond_2pHB_shielding_time_series.cpp
        tests/test_larsen_hbond_2pHaB_shielding_time_series.cpp
        tests/test_larsen_hbond_count_time_series.cpp
        tests/test_larsen_hbond_water_term_time_series.cpp
        tests/test_mc_shielding_time_series.cpp
        tests/test_mc_welford.cpp
        tests/test_piquad_shielding_time_series.cpp
        tests/test_ring_neighbourhood_trajectory_stats.cpp
        tests/test_ring_pucker_time_series.cpp
        tests/test_ringchi_shielding_time_series.cpp
        tests/test_rmsd_spike_dft_coord.cpp
        tests/test_rmsd_tracking.cpp
        tests/test_sasa_time_series.cpp
        tests/test_sasa_welford.cpp
        tests/test_tripeptide_bb_method_tag_time_series.cpp
        tests/test_tripeptide_bb_residual_vec_time_series.cpp
        tests/test_tripeptide_bb_shielding_time_series.cpp
        tests/test_tripeptide_neighbor_residual_vec_next_time_series.cpp
        tests/test_tripeptide_neighbor_residual_vec_prev_time_series.cpp
        tests/test_tripeptide_neighbor_shielding_time_series.cpp
        tests/test_water_field_time_series.cpp
        tests/test_water_field_welford.cpp
)


# -- mopac_tests -----------------------------------------------------------
# MOPAC ConfResult (PM7+MOZYME charges, ~10 min/run) + MOPAC-dependent TRs
# (charge welford, bond order welford, coulomb / mcconnell shielding TS,
# ff14SB reconciliation) + full_pipeline + mutation_delta. Kept isolated
# so day-to-day ctest can skip the slow (~10 min/run) MOPAC fixtures;
# libmopac itself is a required link dep (checked fail-loud at configure).
nmr_test(mopac_tests
    LINK highfive
    SOURCES
        tests/test_full_pipeline.cpp
        tests/test_mopac_bond_order_welford.cpp
        tests/test_mopac_charge_welford.cpp
        tests/test_mopac_coulomb_shielding_time_series.cpp
        tests/test_mopac_mc_shielding_time_series.cpp
        tests/test_mopac_result.cpp
        tests/test_mopac_vs_ff14sb_reconciliation.cpp
        tests/test_mutation_delta.cpp
)


# -- batch_tests -----------------------------------------------------------
# Multi-protein sweeps across consolidated/ (~2 hr wall). Isolated so
# day-to-day ctest never touches them; run explicitly via
# `ctest -R batch_tests`.
nmr_test(batch_tests
    SOURCES
        tests/test_batch_biot_savart_haigh_mallion.cpp
        tests/test_batch_coulomb_ringchi.cpp
        tests/test_batch_mcconnell.cpp
        tests/test_batch_piquad_disp.cpp
)


# -- smoke_tests -----------------------------------------------------------
# End-to-end pipeline validation against bless baselines
# (tests/golden/blessed/). Includes BlessCompare comparator.
nmr_test(smoke_tests
    SOURCES
        tests/BlessCompare.cpp
        tests/test_smoke.cpp
)


# ----------------------------------------------------------------------------
# string_barrier_tests
#
# Standalone discipline test: greps the source tree for cifpp / RDKit /
# gemmi use outside the substrate-generator boundary. Does NOT link
# nmr_shielding (linking would force libtorch + every external dep for
# a grep that runs in <1 s) and uses NMR_SOURCE_DIR instead of
# NMR_TEST_DATA_DIR (its fixture is the tree itself). Hand-rolled by
# design -- the nmr_test() macro above is for tests that exercise the
# library.
# ----------------------------------------------------------------------------
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
gtest_discover_tests(string_barrier_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES ENVIRONMENT "${NMR_TEST_ENVIRONMENT}"
)


# ============================================================================
# Example driver -- builds against the library; not a test.
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

option(NMR_BUILD_VIEWER "Build nmr-viewer (Qt6+VTK)" OFF)

# ui/ still depends on the legacy JobSpec API and has not yet been migrated
# to nmr::cli. NMR_BUILD_VIEWER defaults OFF so the top-level build stays
# green; flip it ON once ui/main_viewer.cpp + MainWindow + ComputeWorker +
# RestServer have been ported.
if(NMR_BUILD_VIEWER AND Qt6_FOUND AND VTK_FOUND)
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
        # Typed read boundary for `nmr_extract --trajectory` H5 output.
        # See ui/src/TrajectoryH5.h.
        ui/src/TrajectoryH5.cpp
    )
    target_include_directories(nmr-viewer PRIVATE
        ${CMAKE_SOURCE_DIR}/ui/src
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
    message(STATUS "Building nmr-viewer (NMR_BUILD_VIEWER=ON, Qt6+VTK found)")
elseif(NMR_BUILD_VIEWER)
    message(STATUS "NMR_BUILD_VIEWER=ON but Qt6/VTK not found — skipping nmr-viewer")
else()
    message(STATUS "Skipping nmr-viewer (NMR_BUILD_VIEWER=OFF; ui/ pending nmr::cli migration)")
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
