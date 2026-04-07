cmake_minimum_required(VERSION 3.20)
project(nmr-shielding LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ============================================================================
# Dependencies
# ============================================================================

find_package(Eigen3 REQUIRED)
find_package(dssp REQUIRED)  # brings cifpp::cifpp with mm module

# OpenBabel3 for bond perception and hybridisation
# Installed in miniforge3. The cmake config is there but does not
# export include dirs properly, so we set them explicitly.
set(OPENBABEL_INCLUDE "/home/jessicalh/miniforge3/include/openbabel3"
    CACHE PATH "OpenBabel3 include directory")
set(OPENBABEL_LIB "/home/jessicalh/miniforge3/lib/libopenbabel.so"
    CACHE FILEPATH "OpenBabel3 library")
if(NOT EXISTS ${OPENBABEL_LIB})
    message(FATAL_ERROR "OpenBabel3 not found at ${OPENBABEL_LIB}")
endif()
message(STATUS "OpenBabel3: ${OPENBABEL_LIB}")

# GROMACS 2026 for TPR reading (fleet loader)
# Linked against the installed libgromacs.so; internal headers from
# the source tree because GROMACS does not install them as public API.
set(GROMACS_LIB "/home/jessicalh/gromacs/lib/libgromacs.so"
    CACHE FILEPATH "GROMACS library")
set(GROMACS_SRC "/home/jessicalh/builds/gromacs-2026.0/src"
    CACHE PATH "GROMACS source tree (for internal headers)")
set(GROMACS_BUILD "/home/jessicalh/builds/gromacs-2026.0/build"
    CACHE PATH "GROMACS build tree (for generated headers)")
if(NOT EXISTS ${GROMACS_LIB})
    message(WARNING "libgromacs not found — fleet loader will not link")
endif()
message(STATUS "GROMACS: ${GROMACS_LIB}")

# Richardson lab reduce (hydrogen placement for protonation)
# Pre-built as static library from /home/jessicalh/builds/reduce-src
set(REDUCE_SRC "/home/jessicalh/builds/reduce-src"
    CACHE PATH "reduce source tree")
set(REDUCE_LIB "/home/jessicalh/builds/reduce-src/build/reduce_src/libreducelib.a"
    CACHE FILEPATH "reduce static library")
set(REDUCE_LIBPDB "/home/jessicalh/builds/reduce-src/build/libpdb/libpdb++.a"
    CACHE FILEPATH "reduce libpdb++ static library")
set(REDUCE_TOOLCLASSES "/home/jessicalh/builds/reduce-src/build/toolclasses/libtoolclasses.a"
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
target_link_libraries(apbs_bridge PRIVATE
    apbs_routines apbs_mg apbs_generic apbs_pmgc mc punc maloc
)

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
)

target_compile_options(nmr_shielding PRIVATE -Wall -Wextra -Wpedantic)
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
)

# ============================================================================
# CLI executable
# ============================================================================

add_executable(nmr_extract src/nmr_extract.cpp)
target_link_libraries(nmr_extract PRIVATE nmr_shielding)

# ============================================================================
# Unit tests
# ============================================================================

enable_testing()

add_executable(nmr_tests
    tests/test_main.cpp
    tests/TestEnvironment.cpp
    tests/test_spherical_tensor.cpp
    tests/test_amino_acid.cpp
    tests/test_ring_hierarchy.cpp
    tests/test_object_model.cpp
    tests/test_pdb_loading.cpp
    tests/test_geometry_result.cpp
    tests/test_demo_result.cpp
    tests/test_two_conformations.cpp
    tests/test_traversal_dump.cpp
    tests/test_naming_registry.cpp
    tests/test_iupac_atom_identity.cpp
    tests/test_dssp_result.cpp
    tests/test_apbs_field_result.cpp
    tests/test_apbs_wired.cpp
    tests/test_protonation_detection.cpp
    tests/test_mopac_result.cpp
    tests/test_atom_flat.cpp
    tests/test_runtime_environment.cpp
    tests/test_foundation_results.cpp
    tests/test_apbs_ff14sb.cpp
    tests/test_full_pipeline.cpp
    tests/test_protonation_pipeline.cpp
    tests/test_mutation_delta.cpp
    tests/test_mcconnell_result.cpp
    tests/test_batch_mcconnell.cpp
    tests/test_coulomb_result.cpp
    tests/test_ring_susceptibility_result.cpp
    tests/test_batch_coulomb_ringchi.cpp
    tests/test_hbond_result.cpp
    tests/test_biot_savart_result.cpp
    tests/test_haigh_mallion_result.cpp
    tests/test_batch_biot_savart_haigh_mallion.cpp
    tests/test_pi_quadrupole_result.cpp
    tests/test_dispersion_result.cpp
    tests/test_batch_piquad_disp.cpp
    tests/test_pipeline_and_sample.cpp
    tests/test_calculation_runner.cpp
    tests/test_write_features.cpp
    tests/test_fleet_loader.cpp
    tests/test_calculator_config.cpp
)

target_compile_definitions(nmr_tests PRIVATE
    NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)
target_link_libraries(nmr_tests PRIVATE
    nmr_shielding
    GTest::gtest
)

include(GoogleTest)
gtest_add_tests(TARGET nmr_tests)

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
    )
    target_include_directories(nmr-viewer PRIVATE
        ${CMAKE_SOURCE_DIR}/ui/src
    )
    target_compile_definitions(nmr-viewer PRIVATE
        NMR_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    )
    target_link_libraries(nmr-viewer PRIVATE
        nmr_shielding
        Qt6::Widgets Qt6::OpenGLWidgets Qt6::Network
        ${VTK_LIBRARIES}
    )
    vtk_module_autoinit(TARGETS nmr-viewer MODULES ${VTK_LIBRARIES})
    message(STATUS "Qt6 + VTK found — building nmr-viewer")
else()
    message(STATUS "Qt6 or VTK not found — skipping nmr-viewer")
endif()
