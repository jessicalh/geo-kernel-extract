# tools/topology/CMakeLists.txt
#
# Build-time generator for the LegacyAmberTopology semantic tables.
# Linked against cifpp (CCD reader) and RDKit (CIP labeler, canonical
# rank, aromatic + hybridisation perception). Reads
# data/ccd/components.cif, runs reconciliation per the precedence
# table specified in spec/plan/topology-substrate-implementation-plan-2026-05-05.md,
# and emits src/generated/LegacyAmberSemanticTables.cpp plus a
# structured process log for inspection.
#
# This generator is OPT-IN: it is built only when
# -DNMR_BUILD_TABLE_GENERATOR=ON is passed at cmake configure time.
# The runtime library does not link RDKit.
#
# Run:
#   cmake -B build-gen -DNMR_BUILD_TABLE_GENERATOR=ON
#   cmake --build build-gen --target build_semantic_tables
#   ./build-gen/tools/topology/build_semantic_tables \
#       --ccd data/ccd/components.cif \
#       --output src/generated/LegacyAmberSemanticTables.cpp \
#       --log    src/generated/LegacyAmberSemanticTables.log.txt
#

# RDKit installation. Cache variables are the primary override surface;
# environment variables seed the cache at first configure. RDKit can come
# from a system package, vcpkg, conda/micromamba, or a source install.
if(NOT RDKIT_ROOT AND NOT "$ENV{RDKIT_ROOT}" STREQUAL "")
    set(RDKIT_ROOT "$ENV{RDKIT_ROOT}" CACHE PATH "RDKit installation root" FORCE)
endif()
if(NOT RDKIT_ROOT AND NOT "$ENV{NMR_RDKIT_ROOT}" STREQUAL "")
    set(RDKIT_ROOT "$ENV{NMR_RDKIT_ROOT}" CACHE PATH "RDKit installation root" FORCE)
endif()
if(NOT RDKIT_INCLUDE AND NOT "$ENV{NMR_RDKIT_INCLUDE}" STREQUAL "")
    set(RDKIT_INCLUDE "$ENV{NMR_RDKIT_INCLUDE}" CACHE PATH "RDKit include directory containing GraphMol/GraphMol.h" FORCE)
endif()
if(NOT RDKIT_LIB_DIR AND NOT "$ENV{NMR_RDKIT_LIB_DIR}" STREQUAL "")
    set(RDKIT_LIB_DIR "$ENV{NMR_RDKIT_LIB_DIR}" CACHE PATH "RDKit library directory" FORCE)
endif()
set(RDKIT_ROOT "${RDKIT_ROOT}"
    CACHE PATH "RDKit installation root")
set(RDKIT_INCLUDE "${RDKIT_INCLUDE}"
    CACHE PATH "RDKit include directory containing GraphMol/GraphMol.h")
set(RDKIT_LIB_DIR "${RDKIT_LIB_DIR}"
    CACHE PATH "RDKit library directory")

if(NOT RDKIT_INCLUDE)
    find_path(RDKIT_INCLUDE
        NAMES GraphMol/GraphMol.h
        HINTS "${RDKIT_ROOT}"
        PATH_SUFFIXES include/rdkit rdkit
        DOC "RDKit include directory containing GraphMol/GraphMol.h")
endif()
if(NOT RDKIT_LIB_DIR AND RDKIT_ROOT)
    foreach(_rdkit_lib_candidate "${RDKIT_ROOT}/lib" "${RDKIT_ROOT}/lib64")
        if(EXISTS "${_rdkit_lib_candidate}")
            set(RDKIT_LIB_DIR "${_rdkit_lib_candidate}"
                CACHE PATH "RDKit library directory" FORCE)
            break()
        endif()
    endforeach()
endif()
if(NOT RDKIT_LIB_DIR AND NOT "$ENV{CONDA_PREFIX}" STREQUAL ""
   AND EXISTS "$ENV{CONDA_PREFIX}/lib")
    set(RDKIT_LIB_DIR "$ENV{CONDA_PREFIX}/lib"
        CACHE PATH "RDKit library directory" FORCE)
endif()
if(NOT EXISTS "${RDKIT_INCLUDE}/GraphMol/GraphMol.h")
    message(FATAL_ERROR
        "RDKit headers not found at ${RDKIT_INCLUDE}. "
        "Set RDKIT_INCLUDE, RDKIT_ROOT, NMR_RDKIT_INCLUDE, or NMR_RDKIT_ROOT.")
endif()
message(STATUS "RDKit headers: ${RDKIT_INCLUDE}")
message(STATUS "RDKit libs:    ${RDKIT_LIB_DIR}")

# RDKit is split into many shared libraries; we use the subset needed
# for chemistry perception + CIP labelling.
set(_rdkit_components
    GraphMol
    FileParsers
    SmilesParse
    MolStandardize
    CIPLabeler
    DataStructs
    RDGeneral
    RDGeometryLib)
set(RDKIT_LIBRARIES)
foreach(_component ${_rdkit_components})
    find_library(_rdkit_lib_${_component}
        NAMES RDKit${_component}
        HINTS "${RDKIT_LIB_DIR}"
        DOC "RDKit ${_component} library")
    if(NOT _rdkit_lib_${_component})
        message(FATAL_ERROR "RDKit ${_component} library not found. "
                            "Set RDKIT_LIB_DIR or NMR_RDKIT_LIB_DIR.")
    endif()
    list(APPEND RDKIT_LIBRARIES "${_rdkit_lib_${_component}}")
endforeach()

add_executable(build_semantic_tables
    build_semantic_tables.cpp
)

target_compile_features(build_semantic_tables PRIVATE cxx_std_17)
target_compile_options(build_semantic_tables PRIVATE -Wall -Wextra)

target_include_directories(build_semantic_tables PRIVATE
    ${RDKIT_INCLUDE}
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(build_semantic_tables PRIVATE
    ${RDKIT_LIBRARIES}
    cifpp::cifpp
    Eigen3::Eigen          # Needed by src/Types.h (Vec3, Mat3) -- the
                           # generator includes Types.h to share the
                           # Element + Hybridisation enums with the
                           # runtime. Eigen does not cross any string
                           # barrier; only the typed-enum subset of
                           # Types.h is consumed.
)

# Set rpath so the binary can find RDKit shared libraries at runtime
# without LD_LIBRARY_PATH gymnastics.
set_target_properties(build_semantic_tables PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${RDKIT_LIB_DIR}"
)
