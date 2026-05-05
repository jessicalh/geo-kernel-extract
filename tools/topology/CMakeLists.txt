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

# RDKit -- micromamba install path used as default; can be overridden
# by passing -DRDKIT_ROOT=/path at configure time.
set(RDKIT_ROOT "/home/jessica/micromamba/envs/mm"
    CACHE PATH "RDKit installation root (must contain include/rdkit + lib)")
set(RDKIT_INCLUDE "${RDKIT_ROOT}/include/rdkit")
set(RDKIT_LIB_DIR "${RDKIT_ROOT}/lib")
if(NOT EXISTS "${RDKIT_INCLUDE}/GraphMol/GraphMol.h")
    message(FATAL_ERROR
        "RDKit headers not found at ${RDKIT_INCLUDE}. "
        "Set -DRDKIT_ROOT=/path/to/rdkit-root that contains include/rdkit + lib.")
endif()
message(STATUS "RDKit headers: ${RDKIT_INCLUDE}")
message(STATUS "RDKit libs:    ${RDKIT_LIB_DIR}")

# RDKit is split into many shared libraries; we use the subset needed
# for chemistry perception + CIP labelling.
set(RDKIT_LIBRARIES
    ${RDKIT_LIB_DIR}/libRDKitGraphMol.so
    ${RDKIT_LIB_DIR}/libRDKitFileParsers.so
    ${RDKIT_LIB_DIR}/libRDKitSmilesParse.so
    ${RDKIT_LIB_DIR}/libRDKitMolStandardize.so
    ${RDKIT_LIB_DIR}/libRDKitCIPLabeler.so
    ${RDKIT_LIB_DIR}/libRDKitDataStructs.so
    ${RDKIT_LIB_DIR}/libRDKitRDGeneral.so
    ${RDKIT_LIB_DIR}/libRDKitRDGeometryLib.so
)
foreach(_lib ${RDKIT_LIBRARIES})
    if(NOT EXISTS ${_lib})
        message(FATAL_ERROR "RDKit library not found: ${_lib}")
    endif()
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
