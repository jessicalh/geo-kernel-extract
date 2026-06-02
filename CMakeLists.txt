cmake_minimum_required(VERSION 3.20)
project(nmr-shielding LANGUAGES C CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Always emit compile_commands.json so clang-tidy / IWYU / scan-build
# / clangd can find translation-unit-specific flags. Lives in the
# build directory; scripts/lint/* + clangd consume it from there.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(NMR_PROFILE "local" CACHE STRING
    "Configuration profile: local or producer-full")
set_property(CACHE NMR_PROFILE PROPERTY STRINGS local producer-full)
option(NMR_REQUIRE_PORTABLE_DEPS
    "Fail configure when producer dependencies resolve through machine-local paths"
    OFF)
set(NMR_PORTABLE_DEP_ALLOW_REGEXES "" CACHE STRING
    "Semicolon-separated regexes for additional dependency paths allowed by NMR_REQUIRE_PORTABLE_DEPS")
if(NMR_PROFILE STREQUAL "producer-full")
    set(NMR_REQUIRE_PORTABLE_DEPS ON CACHE BOOL
        "Fail configure when producer dependencies resolve through machine-local paths" FORCE)
    set(NMR_USE_INSTALL_DATA_DIR ON CACHE BOOL
        "Compile NMR_DATA_DIR to the install data directory instead of the source tree" FORCE)
    set(NMR_ENABLE_NATIVE_ARCH OFF CACHE BOOL
        "Compile nmr_shielding with -march=native on compilers that support it" FORCE)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "/opt/nmr-shielding" CACHE PATH
            "Install path prefix" FORCE)
    endif()
elseif(NOT NMR_PROFILE STREQUAL "local")
    message(FATAL_ERROR
        "Unknown NMR_PROFILE='${NMR_PROFILE}'. Expected 'local' or "
        "'producer-full'.")
endif()

set(NMR_INSTALL_DATA_DIR "${CMAKE_INSTALL_DATAROOTDIR}/nmr-shielding/data"
    CACHE STRING "Relative install directory for nmr-shielding runtime data")
option(NMR_USE_INSTALL_DATA_DIR
    "Compile NMR_DATA_DIR to the install data directory instead of the source tree"
    OFF)
if(NMR_USE_INSTALL_DATA_DIR)
    if(IS_ABSOLUTE "${NMR_INSTALL_DATA_DIR}")
        set(NMR_RUNTIME_DATA_DIR "${NMR_INSTALL_DATA_DIR}")
    else()
        set(NMR_RUNTIME_DATA_DIR "${CMAKE_INSTALL_PREFIX}/${NMR_INSTALL_DATA_DIR}")
    endif()
else()
    set(NMR_RUNTIME_DATA_DIR "${CMAKE_SOURCE_DIR}/data")
endif()

# Default build type: optimised with debug symbols.
# Override with: cmake -DCMAKE_BUILD_TYPE=Release ..
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    message(STATUS "Build type: RelWithDebInfo (default)")
endif()

# ============================================================================
# Machine-local dependency knobs
#
# Public override surface:
#   1. CMake cache variables (-D..., CMakeUserPresets.json)
#   2. Matching NMR_* environment variables at first configure
#   3. Standard system discovery / compatibility probes
#
# The committed build files must not encode a username or a private absolute
# path. Keep per-machine choices in the cache, environment, or ignored
# CMakeUserPresets.json.
# ============================================================================

function(nmr_cache_from_env cache_var env_var cache_type help_text)
    if((NOT DEFINED ${cache_var}
        OR "${${cache_var}}" STREQUAL ""
        OR "${${cache_var}}" MATCHES "-NOTFOUND$")
       AND NOT "$ENV{${env_var}}" STREQUAL "")
        set(${cache_var} "$ENV{${env_var}}" CACHE ${cache_type} "${help_text}" FORCE)
    endif()
endfunction()

function(nmr_cache_first_existing cache_var cache_type help_text)
    if(DEFINED ${cache_var}
       AND NOT "${${cache_var}}" STREQUAL ""
       AND NOT "${${cache_var}}" MATCHES "-NOTFOUND$")
        return()
    endif()
    foreach(_candidate ${ARGN})
        if(NOT "${_candidate}" STREQUAL "" AND EXISTS "${_candidate}")
            set(${cache_var} "${_candidate}" CACHE ${cache_type} "${help_text}" FORCE)
            return()
        endif()
    endforeach()
endfunction()

function(nmr_require_existing cache_var description)
    if(NOT DEFINED ${cache_var}
       OR "${${cache_var}}" STREQUAL ""
       OR "${${cache_var}}" MATCHES "-NOTFOUND$"
       OR NOT EXISTS "${${cache_var}}")
        message(FATAL_ERROR
            "${description} not found at '${${cache_var}}'. Set ${cache_var} "
            "or the documented matching NMR_* environment variable.")
    endif()
endfunction()

function(nmr_require_contains cache_var relative_path description)
    if(NOT DEFINED ${cache_var}
       OR "${${cache_var}}" STREQUAL ""
       OR "${${cache_var}}" MATCHES "-NOTFOUND$"
       OR NOT EXISTS "${${cache_var}}/${relative_path}")
        message(FATAL_ERROR
            "${description} not found under '${${cache_var}}' "
            "(expected ${relative_path}). Set ${cache_var} or the "
            "documented matching NMR_* environment variable.")
    endif()
endfunction()

set(NMR_MACHINE_LOCAL_PATH_REGEX "^(/home/|/shared/|/mnt/|/opt/orca/|/usr/local/cuda/)")
set(_nmr_portability_errors)

function(nmr_require_portable_paths label)
    if(NOT NMR_REQUIRE_PORTABLE_DEPS)
        return()
    endif()

    foreach(_path IN LISTS ARGN)
        if("${_path}" STREQUAL "" OR "${_path}" MATCHES "-NOTFOUND$")
            continue()
        endif()
        if(NOT IS_ABSOLUTE "${_path}")
            continue()
        endif()

        set(_allowed OFF)
        foreach(_allow_regex IN LISTS NMR_PORTABLE_DEP_ALLOW_REGEXES)
            if(NOT "${_allow_regex}" STREQUAL "" AND "${_path}" MATCHES "${_allow_regex}")
                set(_allowed ON)
            endif()
        endforeach()

        if(NOT _allowed AND "${_path}" MATCHES "${NMR_MACHINE_LOCAL_PATH_REGEX}")
            list(APPEND _nmr_portability_errors "${label}: ${_path}")
        endif()
    endforeach()
    set(_nmr_portability_errors "${_nmr_portability_errors}" PARENT_SCOPE)
endfunction()

function(nmr_fail_on_portability_errors)
    if(NOT NMR_REQUIRE_PORTABLE_DEPS)
        return()
    endif()

    if(NOT NMR_USE_INSTALL_DATA_DIR)
        message(FATAL_ERROR
            "NMR_REQUIRE_PORTABLE_DEPS=ON requires "
            "NMR_USE_INSTALL_DATA_DIR=ON so installed builds do not compile "
            "NMR_DATA_DIR to the source tree.")
    endif()

    if(NMR_ENABLE_NATIVE_ARCH)
        message(FATAL_ERROR
            "NMR_REQUIRE_PORTABLE_DEPS=ON requires "
            "NMR_ENABLE_NATIVE_ARCH=OFF so binaries are not tied to this CPU.")
    endif()

    if(_nmr_portability_errors)
        list(REMOVE_DUPLICATES _nmr_portability_errors)
        string(REPLACE ";" "\n  " _nmr_portability_report "${_nmr_portability_errors}")
        message(FATAL_ERROR
            "NMR_REQUIRE_PORTABLE_DEPS=ON but machine-local dependency paths "
            "were resolved:\n  ${_nmr_portability_report}\n"
            "Set these dependencies to a deliberate producer prefix or use "
            "NMR_PROFILE=local for machine-local development.")
    endif()
endfunction()

nmr_cache_from_env(NMR_SYSTEM_LIBRARY_DIR NMR_SYSTEM_LIBRARY_DIR PATH
    "System library directory to prepend to RUNPATH / LD_LIBRARY_PATH for libpq/OpenSSL consistency")
if(NOT NMR_SYSTEM_LIBRARY_DIR)
    if(CMAKE_LIBRARY_ARCHITECTURE AND EXISTS "/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}")
        set(NMR_SYSTEM_LIBRARY_DIR "/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
            CACHE PATH "System library directory to prepend to RUNPATH / LD_LIBRARY_PATH for libpq/OpenSSL consistency" FORCE)
    elseif(EXISTS "/usr/lib/x86_64-linux-gnu")
        set(NMR_SYSTEM_LIBRARY_DIR "/usr/lib/x86_64-linux-gnu"
            CACHE PATH "System library directory to prepend to RUNPATH / LD_LIBRARY_PATH for libpq/OpenSSL consistency" FORCE)
    endif()
endif()
option(NMR_INSTALL_RPATH_USE_ORIGIN
    "Add $ORIGIN-relative install RUNPATH for bundled producer dependencies"
    ON)
option(NMR_INSTALL_RPATH_USE_LINK_PATH
    "Allow CMake to append external link directories to installed RUNPATH"
    OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH "${NMR_INSTALL_RPATH_USE_LINK_PATH}")

# ============================================================================
# RUNPATH: prepend the selected system library directory so the system libpq +
# libssl + libcrypto trio resolves consistently. Other linked libs
# (libmopac, libreducelib, libgromacs_mpi, libtorch, libnvrtc, ...)
# bring their own dependency lib dirs along as RUNPATH entries
# automatically; without this prepend the loader can pick libpq/libssl
# from one dependency environment but libcrypto from the host system,
# and the load fails. We
# don't USE OpenSSL — libpq's NEEDED entries pull libssl/libcrypto
# at process startup regardless of runtime config — but the two
# halves of the SSL pair must come from the same install.
if(NMR_SYSTEM_LIBRARY_DIR)
    list(PREPEND CMAKE_BUILD_RPATH "${NMR_SYSTEM_LIBRARY_DIR}")
    list(APPEND CMAKE_INSTALL_RPATH "${NMR_SYSTEM_LIBRARY_DIR}")
endif()
if(NMR_INSTALL_RPATH_USE_ORIGIN)
    list(PREPEND CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set(CMAKE_BUILD_RPATH_USE_ORIGIN OFF)

# ============================================================================
# Dependencies
# ============================================================================

find_package(Eigen3 REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

nmr_cache_from_env(NMR_DSSP_ROOT NMR_DSSP_ROOT PATH
    "DSSP installation root; optional hint before find_package(dssp)")
if(NMR_DSSP_ROOT)
    list(PREPEND CMAKE_PREFIX_PATH "${NMR_DSSP_ROOT}")
endif()
find_package(dssp QUIET)  # brings cifpp::cifpp with mm module

# PostgreSQL libpq — required for TripeptideDftTable (ProCS15 tripeptide
# DFT data lookup against the local tensorcs15 replica). The DSN is
# read from ~/.nmr_tools.toml [databases].tensorcs15. No DSN configured
# means TripeptideBackboneShieldingResult skips at Compute time;
# linking is unconditional so the calculator code is always present.
# PostgreSQL libpq — explicitly resolved so libpq, libssl, and libcrypto
# come from one installation. We do NOT use find_package here because
# dependency environments can provide a libpq that pulls a different
# OpenSSL pair than the host linker resolves by default.
nmr_cache_from_env(PostgreSQL_LIBRARY NMR_POSTGRESQL_LIBRARY FILEPATH
    "libpq shared library")
nmr_cache_from_env(PostgreSQL_INCLUDE_DIR NMR_POSTGRESQL_INCLUDE_DIR PATH
    "libpq include directory containing libpq-fe.h")
find_library(PostgreSQL_LIBRARY
    NAMES pq
    PATHS "${NMR_SYSTEM_LIBRARY_DIR}"
    DOC "libpq shared library")
find_path(PostgreSQL_INCLUDE_DIR
    NAMES libpq-fe.h
    PATH_SUFFIXES postgresql
    DOC "libpq include directory containing libpq-fe.h")
if(NOT PostgreSQL_LIBRARY OR NOT EXISTS "${PostgreSQL_LIBRARY}")
    message(FATAL_ERROR "libpq not found. Set PostgreSQL_LIBRARY or "
                        "NMR_POSTGRESQL_LIBRARY to the desired libpq.so.")
endif()
if(NOT PostgreSQL_INCLUDE_DIR OR NOT EXISTS "${PostgreSQL_INCLUDE_DIR}/libpq-fe.h")
    message(FATAL_ERROR "libpq headers not found. Set PostgreSQL_INCLUDE_DIR "
                        "or NMR_POSTGRESQL_INCLUDE_DIR to the directory "
                        "containing libpq-fe.h.")
endif()
add_library(PostgreSQL::PostgreSQL UNKNOWN IMPORTED)
set_target_properties(PostgreSQL::PostgreSQL PROPERTIES
    IMPORTED_LOCATION "${PostgreSQL_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PostgreSQL_INCLUDE_DIR}")
message(STATUS "PostgreSQL libpq: ${PostgreSQL_LIBRARY}")
if(NOT dssp_FOUND)
    # Manual fallback for machines without dssp cmake config (Ubuntu 22.04).
    # Requires libdssp + dssp.hpp, cifpp from apt or another CMAKE_PREFIX_PATH.
    nmr_cache_from_env(DSSP_LIB NMR_DSSP_LIB FILEPATH
        "DSSP library for manual fallback")
    nmr_cache_from_env(DSSP_INCLUDE NMR_DSSP_INCLUDE PATH
        "DSSP include directory containing dssp.hpp")
    find_library(DSSP_LIB
        NAMES dssp
        HINTS "${NMR_DSSP_ROOT}"
        PATH_SUFFIXES lib lib64
        DOC "DSSP library for manual fallback")
    find_path(DSSP_INCLUDE
        NAMES dssp.hpp
        HINTS "${NMR_DSSP_ROOT}"
        PATH_SUFFIXES include
        DOC "DSSP include directory containing dssp.hpp")
    if(NOT DSSP_LIB OR NOT DSSP_INCLUDE)
        message(FATAL_ERROR "dssp not found via find_package or manual "
                            "fallback. Set NMR_DSSP_ROOT, DSSP_LIB, and/or "
                            "DSSP_INCLUDE.")
    endif()
    find_package(cifpp REQUIRED)
    add_library(dssp::dssp STATIC IMPORTED)
    set_target_properties(dssp::dssp PROPERTIES
        IMPORTED_LOCATION "${DSSP_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${DSSP_INCLUDE}"
        INTERFACE_LINK_LIBRARIES "cifpp::cifpp"
    )
    message(STATUS "dssp: manual fallback — ${DSSP_LIB}")
endif()

# OpenBabel3 for bond perception and hybridisation. The cmake config is
# not reliable across package managers, so resolve the include + library
# explicitly but allow any install root.
nmr_cache_from_env(NMR_OPENBABEL_ROOT NMR_OPENBABEL_ROOT PATH
    "OpenBabel installation root")
nmr_cache_from_env(OPENBABEL_INCLUDE NMR_OPENBABEL_INCLUDE PATH
    "OpenBabel3 include directory containing openbabel/mol.h")
nmr_cache_from_env(OPENBABEL_LIB NMR_OPENBABEL_LIB FILEPATH
    "OpenBabel3 shared library")
find_path(OPENBABEL_INCLUDE
    NAMES openbabel/mol.h
    HINTS "${NMR_OPENBABEL_ROOT}"
    PATH_SUFFIXES include/openbabel3 openbabel3
    DOC "OpenBabel3 include directory containing openbabel/mol.h")
find_library(OPENBABEL_LIB
    NAMES openbabel
    HINTS "${NMR_OPENBABEL_ROOT}"
    PATH_SUFFIXES lib lib64
    DOC "OpenBabel3 shared library")
if(NOT OPENBABEL_INCLUDE OR NOT EXISTS "${OPENBABEL_INCLUDE}/openbabel/mol.h")
    message(FATAL_ERROR "OpenBabel3 headers not found. Set OPENBABEL_INCLUDE "
                        "or NMR_OPENBABEL_INCLUDE.")
endif()
if(NOT OPENBABEL_LIB OR NOT EXISTS "${OPENBABEL_LIB}")
    message(FATAL_ERROR "OpenBabel3 library not found. Set OPENBABEL_LIB "
                        "or NMR_OPENBABEL_LIB.")
endif()
message(STATUS "OpenBabel3: ${OPENBABEL_LIB}")

# GROMACS 2026 for TPR reading (fleet loader)
# Linked against the installed libgromacs.so; internal headers from
# the source tree because GROMACS does not install them as public API.
nmr_cache_from_env(NMR_GROMACS_ROOT NMR_GROMACS_ROOT PATH
    "GROMACS install root; optional hint for libgromacs")
nmr_cache_from_env(GROMACS_LIB NMR_GROMACS_LIB FILEPATH
    "GROMACS library")
nmr_cache_from_env(GROMACS_SRC NMR_GROMACS_SRC PATH
    "GROMACS source tree (for internal headers)")
nmr_cache_from_env(GROMACS_BUILD NMR_GROMACS_BUILD PATH
    "GROMACS generated header root")
find_library(GROMACS_LIB
    NAMES gromacs_mpi gromacs
    HINTS "${NMR_GROMACS_ROOT}"
    PATH_SUFFIXES lib lib64
    DOC "GROMACS library")
nmr_cache_first_existing(GROMACS_LIB FILEPATH "GROMACS library"
    "$ENV{HOME}/gromacs/build/gromacs-2026/lib/libgromacs_mpi.so"
    "$ENV{HOME}/build/gmx2026/gromacs-2026-build/lib/libgromacs_mpi.so"
    "/opt/gromacs/2026.0/lib/libgromacs_mpi.so")
nmr_cache_first_existing(GROMACS_SRC PATH "GROMACS source tree (for internal headers)"
    "$ENV{HOME}/gromacs/src/gromacs-2026.0/src"
    "$ENV{HOME}/build/gmx2026/gromacs-2026.0/src")
nmr_cache_first_existing(GROMACS_BUILD PATH "GROMACS generated header root"
    "$ENV{HOME}/gromacs/build/gromacs-2026"
    "$ENV{HOME}/build/gmx2026/gromacs-2026-build")
if(NOT GROMACS_LIB OR NOT EXISTS "${GROMACS_LIB}")
    message(FATAL_ERROR "libgromacs not found at ${GROMACS_LIB} — required "
                        "(the library links it unconditionally for trajectory "
                        "reading). Set GROMACS_LIB or NMR_GROMACS_LIB.")
endif()
if(NOT GROMACS_SRC OR NOT EXISTS "${GROMACS_SRC}/gromacs")
    message(FATAL_ERROR "GROMACS source headers not found at ${GROMACS_SRC}. "
                        "Set GROMACS_SRC or NMR_GROMACS_SRC.")
endif()
if(NOT GROMACS_BUILD OR NOT EXISTS "${GROMACS_BUILD}/api/legacy/include/gromacs/libgromacs_export.h")
    message(FATAL_ERROR "GROMACS generated headers not found at ${GROMACS_BUILD}. "
                        "Set GROMACS_BUILD or NMR_GROMACS_BUILD.")
endif()
message(STATUS "GROMACS: ${GROMACS_LIB}")

# Richardson lab reduce (hydrogen placement for protonation)
nmr_cache_from_env(REDUCE_SRC NMR_REDUCE_SRC PATH
    "reduce source tree")
nmr_cache_from_env(REDUCE_LIB NMR_REDUCE_LIB FILEPATH
    "reduce static library")
nmr_cache_from_env(REDUCE_LIBPDB NMR_REDUCE_LIBPDB FILEPATH
    "reduce libpdb++ static library")
nmr_cache_from_env(REDUCE_TOOLCLASSES NMR_REDUCE_TOOLCLASSES FILEPATH
    "reduce toolclasses static library")
nmr_cache_first_existing(REDUCE_SRC PATH "reduce source tree"
    "$ENV{HOME}/builds/reduce-src")
nmr_cache_first_existing(REDUCE_LIB FILEPATH "reduce static library"
    "${REDUCE_SRC}/build/reduce_src/libreducelib.a")
nmr_cache_first_existing(REDUCE_LIBPDB FILEPATH "reduce libpdb++ static library"
    "${REDUCE_SRC}/build/libpdb/libpdb++.a")
nmr_cache_first_existing(REDUCE_TOOLCLASSES FILEPATH "reduce toolclasses static library"
    "${REDUCE_SRC}/build/toolclasses/libtoolclasses.a")
nmr_cache_from_env(REDUCE_HET_DICT NMR_REDUCE_HET_DICT FILEPATH
    "reduce het dictionary")
set(REDUCE_HET_DICT "${REDUCE_SRC}/reduce_wwPDB_het_dict.txt"
    CACHE FILEPATH "reduce het dictionary")
# reduce links three archives unconditionally (reducelib + libpdb++ +
# toolclasses, see target_link_libraries below) — all required.
foreach(_reduce_archive REDUCE_LIB REDUCE_LIBPDB REDUCE_TOOLCLASSES)
    if(NOT EXISTS "${${_reduce_archive}}")
        message(FATAL_ERROR "reduce archive ${_reduce_archive}=${${_reduce_archive}} "
                            "not found — required (linked unconditionally for "
                            "protonation). Set ${_reduce_archive} or the "
                            "matching NMR_REDUCE_* environment variable.")
    endif()
endforeach()
if(NOT EXISTS "${REDUCE_HET_DICT}")
    message(FATAL_ERROR "reduce het dictionary not found at ${REDUCE_HET_DICT}. "
                        "Set REDUCE_HET_DICT or NMR_REDUCE_HET_DICT.")
endif()
message(STATUS "reduce: ${REDUCE_LIB}")

# MOPAC (PM7+MOZYME semiempirical) — linked library, replaces xTB subprocess
nmr_cache_from_env(NMR_MOPAC_ROOT NMR_MOPAC_ROOT PATH
    "MOPAC installation root")
if(NOT NMR_MOPAC_ROOT AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
    set(NMR_MOPAC_ROOT "$ENV{CONDA_PREFIX}" CACHE PATH "MOPAC installation root" FORCE)
elseif(NOT NMR_MOPAC_ROOT AND NOT "$ENV{MAMBA_ROOT_PREFIX}" STREQUAL ""
       AND EXISTS "$ENV{MAMBA_ROOT_PREFIX}/envs/mm")
    set(NMR_MOPAC_ROOT "$ENV{MAMBA_ROOT_PREFIX}/envs/mm" CACHE PATH "MOPAC installation root" FORCE)
endif()
nmr_cache_from_env(MOPAC_INCLUDE NMR_MOPAC_INCLUDE PATH
    "MOPAC include directory (mopac.h)")
nmr_cache_from_env(MOPAC_LIB NMR_MOPAC_LIB FILEPATH
    "MOPAC shared library")
find_path(MOPAC_INCLUDE
    NAMES mopac.h
    HINTS "${NMR_MOPAC_ROOT}"
    PATH_SUFFIXES include
    DOC "MOPAC include directory (mopac.h)")
find_library(MOPAC_LIB
    NAMES mopac
    HINTS "${NMR_MOPAC_ROOT}"
    PATH_SUFFIXES lib lib64
    DOC "MOPAC shared library")
if(NOT MOPAC_INCLUDE OR NOT EXISTS "${MOPAC_INCLUDE}/mopac.h")
    message(FATAL_ERROR "MOPAC headers not found. Set MOPAC_INCLUDE or "
                        "NMR_MOPAC_INCLUDE.")
endif()
if(NOT MOPAC_LIB OR NOT EXISTS "${MOPAC_LIB}")
    message(FATAL_ERROR "libmopac not found at ${MOPAC_LIB} — required "
                        "(the library links it unconditionally; the --mopac "
                        "toggle gates whether MOPAC RUNS, not whether it links). "
                        "Set MOPAC_LIB or NMR_MOPAC_LIB.")
endif()
message(STATUS "MOPAC: ${MOPAC_LIB}")

# ============================================================================
# Header-only external libraries
# ============================================================================

add_library(extern_headers INTERFACE)
target_include_directories(extern_headers INTERFACE ${CMAKE_SOURCE_DIR}/extern)

# HighFive — header-only HDF5 C++ interface (extern/HighFive)
# Link HDF5 directly so the exact HDF5/curl/zlib installation is visible
# in the cache and can be pinned per machine.
nmr_cache_from_env(NMR_HDF5_INCLUDE_DIR NMR_HDF5_INCLUDE_DIR PATH
    "HDF5 C include directory containing hdf5.h")
nmr_cache_from_env(NMR_HDF5_C_LIBRARY NMR_HDF5_C_LIBRARY FILEPATH
    "HDF5 C library")
nmr_cache_from_env(NMR_CURL_LIBRARY NMR_CURL_LIBRARY FILEPATH
    "libcurl library linked with HDF5")
nmr_cache_from_env(NMR_Z_LIBRARY NMR_Z_LIBRARY FILEPATH
    "zlib library linked with HDF5")
find_path(NMR_HDF5_INCLUDE_DIR
    NAMES hdf5.h
    PATHS "/usr/include/hdf5/serial"
    PATH_SUFFIXES hdf5/serial
    DOC "HDF5 C include directory containing hdf5.h")
find_library(NMR_HDF5_C_LIBRARY
    NAMES hdf5
    PATHS "${NMR_SYSTEM_LIBRARY_DIR}/hdf5/serial"
    PATH_SUFFIXES hdf5/serial
    DOC "HDF5 C library")
find_library(NMR_CURL_LIBRARY
    NAMES curl
    PATHS "${NMR_SYSTEM_LIBRARY_DIR}"
    DOC "libcurl library linked with HDF5")
find_library(NMR_Z_LIBRARY
    NAMES z
    PATHS "${NMR_SYSTEM_LIBRARY_DIR}"
    DOC "zlib library linked with HDF5")
nmr_require_contains(NMR_HDF5_INCLUDE_DIR hdf5.h "HDF5 C headers")
nmr_require_existing(NMR_HDF5_C_LIBRARY "HDF5 C library")
nmr_require_existing(NMR_CURL_LIBRARY "libcurl library")
nmr_require_existing(NMR_Z_LIBRARY "zlib library")
add_library(highfive INTERFACE)
target_include_directories(highfive INTERFACE
    ${CMAKE_SOURCE_DIR}/extern/HighFive/include
    ${NMR_HDF5_INCLUDE_DIR})
target_link_libraries(highfive INTERFACE
    ${NMR_HDF5_C_LIBRARY}
    ${NMR_CURL_LIBRARY}
    ${NMR_Z_LIBRARY})

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

nmr_cache_from_env(NMR_TORCH_CMAKE_PREFIX_PATH NMR_TORCH_CMAKE_PREFIX_PATH PATH
    "Directory containing TorchConfig.cmake, or a prefix with share/cmake/Torch")
if(NMR_TORCH_CMAKE_PREFIX_PATH)
    list(PREPEND CMAKE_PREFIX_PATH "${NMR_TORCH_CMAKE_PREFIX_PATH}")
elseif(NOT Torch_DIR)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import pathlib, torch; print(pathlib.Path(torch.__file__).parent / 'share' / 'cmake')"
        OUTPUT_VARIABLE _nmr_torch_cmake_from_python
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _nmr_torch_python_status
    )
    if(_nmr_torch_python_status EQUAL 0 AND EXISTS "${_nmr_torch_cmake_from_python}/Torch/TorchConfig.cmake")
        set(NMR_TORCH_CMAKE_PREFIX_PATH "${_nmr_torch_cmake_from_python}"
            CACHE PATH "Directory containing TorchConfig.cmake, or a prefix with share/cmake/Torch" FORCE)
        list(PREPEND CMAKE_PREFIX_PATH "${NMR_TORCH_CMAKE_PREFIX_PATH}")
        message(STATUS "Torch hint from ${Python3_EXECUTABLE}: ${NMR_TORCH_CMAKE_PREFIX_PATH}")
    endif()
endif()
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
nmr_cache_from_env(NMR_APBS_INCLUDE_DIR NMR_APBS_INCLUDE_DIR PATH
    "APBS include root containing apbs/apbs.h")
nmr_cache_from_env(NMR_FETK_INCLUDE_DIR NMR_FETK_INCLUDE_DIR PATH
    "FETK include root containing fetk/")
nmr_cache_from_env(NMR_MALOC_INCLUDE_DIR NMR_MALOC_INCLUDE_DIR PATH
    "MALOC include root containing maloc/")
find_path(NMR_APBS_INCLUDE_DIR
    NAMES apbs/apbs.h
    DOC "APBS include root containing apbs/apbs.h")
find_path(NMR_FETK_INCLUDE_DIR
    NAMES fetk/punc/punc.h
    DOC "FETK include root containing fetk/")
find_path(NMR_MALOC_INCLUDE_DIR
    NAMES maloc/maloc.h
    DOC "MALOC include root containing maloc/")
nmr_require_contains(NMR_APBS_INCLUDE_DIR apbs/apbs.h "APBS headers")
nmr_require_contains(NMR_MALOC_INCLUDE_DIR maloc/maloc.h "MALOC headers")

set(_nmr_apbs_include_dirs
    ${CMAKE_SOURCE_DIR}/src
    ${NMR_APBS_INCLUDE_DIR}
    ${NMR_APBS_INCLUDE_DIR}/apbs
    ${NMR_MALOC_INCLUDE_DIR}
    ${NMR_MALOC_INCLUDE_DIR}/maloc
)
if(NMR_FETK_INCLUDE_DIR AND NOT "${NMR_FETK_INCLUDE_DIR}" MATCHES "-NOTFOUND$")
    nmr_require_contains(NMR_FETK_INCLUDE_DIR fetk/punc/punc.h "FETK headers")
    list(APPEND _nmr_apbs_include_dirs
        ${NMR_FETK_INCLUDE_DIR}
        ${NMR_FETK_INCLUDE_DIR}/fetk)
endif()
target_include_directories(apbs_bridge PRIVATE ${_nmr_apbs_include_dirs})
# FETK libs (mc, punc) are separate on Ubuntu 24.04 but bundled
# inside libapbs on 22.04.  Link them only if found.
find_library(FETK_MC mc)
find_library(FETK_PUNC punc)
nmr_cache_from_env(NMR_APBS_ROUTINES_LIBRARY NMR_APBS_ROUTINES_LIBRARY FILEPATH
    "APBS routines library")
nmr_cache_from_env(NMR_APBS_MG_LIBRARY NMR_APBS_MG_LIBRARY FILEPATH
    "APBS multigrid library")
nmr_cache_from_env(NMR_APBS_GENERIC_LIBRARY NMR_APBS_GENERIC_LIBRARY FILEPATH
    "APBS generic library")
nmr_cache_from_env(NMR_APBS_PMGC_LIBRARY NMR_APBS_PMGC_LIBRARY FILEPATH
    "APBS PMGC library")
nmr_cache_from_env(NMR_MALOC_LIBRARY NMR_MALOC_LIBRARY FILEPATH
    "MALOC library")
find_library(NMR_APBS_ROUTINES_LIBRARY NAMES apbs_routines DOC "APBS routines library")
find_library(NMR_APBS_MG_LIBRARY NAMES apbs_mg DOC "APBS multigrid library")
find_library(NMR_APBS_GENERIC_LIBRARY NAMES apbs_generic DOC "APBS generic library")
find_library(NMR_APBS_PMGC_LIBRARY NAMES apbs_pmgc DOC "APBS PMGC library")
find_library(NMR_MALOC_LIBRARY NAMES maloc DOC "MALOC library")
nmr_require_existing(NMR_APBS_ROUTINES_LIBRARY "APBS routines library")
nmr_require_existing(NMR_APBS_MG_LIBRARY "APBS multigrid library")
nmr_require_existing(NMR_APBS_GENERIC_LIBRARY "APBS generic library")
nmr_require_existing(NMR_APBS_PMGC_LIBRARY "APBS PMGC library")
nmr_require_existing(NMR_MALOC_LIBRARY "MALOC library")
set(APBS_LIBS
    ${NMR_APBS_ROUTINES_LIBRARY}
    ${NMR_APBS_MG_LIBRARY}
    ${NMR_APBS_GENERIC_LIBRARY}
    ${NMR_APBS_PMGC_LIBRARY}
    ${NMR_MALOC_LIBRARY})
if(FETK_MC)
    nmr_require_existing(FETK_MC "FETK mc library")
    list(APPEND APBS_LIBS ${FETK_MC})
endif()
if(FETK_PUNC)
    nmr_require_existing(FETK_PUNC "FETK punc library")
    list(APPEND APBS_LIBS ${FETK_PUNC})
endif()
target_link_libraries(apbs_bridge PRIVATE ${APBS_LIBS})

nmr_require_portable_paths("PostgreSQL"
    "${PostgreSQL_LIBRARY}" "${PostgreSQL_INCLUDE_DIR}")
nmr_require_portable_paths("DSSP"
    "${NMR_DSSP_ROOT}" "${DSSP_LIB}" "${DSSP_INCLUDE}")
nmr_require_portable_paths("OpenBabel"
    "${NMR_OPENBABEL_ROOT}" "${OPENBABEL_INCLUDE}" "${OPENBABEL_LIB}")
nmr_require_portable_paths("GROMACS"
    "${NMR_GROMACS_ROOT}" "${GROMACS_LIB}" "${GROMACS_SRC}" "${GROMACS_BUILD}")
nmr_require_portable_paths("reduce"
    "${REDUCE_SRC}" "${REDUCE_LIB}" "${REDUCE_LIBPDB}"
    "${REDUCE_TOOLCLASSES}" "${REDUCE_HET_DICT}")
nmr_require_portable_paths("MOPAC"
    "${NMR_MOPAC_ROOT}" "${MOPAC_INCLUDE}" "${MOPAC_LIB}")
nmr_require_portable_paths("HDF5"
    "${NMR_HDF5_INCLUDE_DIR}" "${NMR_HDF5_C_LIBRARY}"
    "${NMR_CURL_LIBRARY}" "${NMR_Z_LIBRARY}")
nmr_require_portable_paths("Torch"
    "${NMR_TORCH_CMAKE_PREFIX_PATH}" "${Torch_DIR}" "${Caffe2_DIR}"
    ${TORCH_INCLUDE_DIRS} ${TORCH_LIBRARIES})
nmr_require_portable_paths("CUDA toolkit"
    "${CMAKE_CUDA_COMPILER}" "${CUDA_TOOLKIT_ROOT_DIR}"
    "${CUDA_NVCC_EXECUTABLE}" "${CUDA_NVRTC_LIB}"
    "${CUDA_CUDART}" "${CUDA_CUDART_LIBRARY}"
    "${CUDA_nvrtc_LIBRARY}" "${CUDA_cudart_LIBRARY}")
nmr_require_portable_paths("APBS/MALOC/FETK"
    "${NMR_APBS_INCLUDE_DIR}" "${NMR_FETK_INCLUDE_DIR}"
    "${NMR_MALOC_INCLUDE_DIR}" ${APBS_LIBS})

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
    src/KernelDynamicsTrajectoryResult.cpp
    src/KernelCoherenceTrajectoryResult.cpp
    src/IRedOrderParameterTrajectoryResult.cpp
    src/ReorientationalDynamicsTrajectoryResult.cpp
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
    src/DihedralAutocorrelationTrajectoryResult.cpp
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

option(NMR_ENABLE_NATIVE_ARCH
    "Compile nmr_shielding with -march=native on compilers that support it"
    ON)
option(NMR_ENABLE_LTO
    "Compile nmr_shielding with -flto on compilers that support it"
    ON)
set(NMR_EXTRA_CXX_FLAGS "" CACHE STRING
    "Additional C++ compile flags for nmr_shielding, semicolon-separated")
nmr_fail_on_portability_errors()
set(_nmr_compile_options -Wall -Wextra -Wpedantic ${TORCH_CXX_FLAGS})
if(NMR_ENABLE_NATIVE_ARCH)
    list(APPEND _nmr_compile_options -march=native)
endif()
if(NMR_ENABLE_LTO)
    list(APPEND _nmr_compile_options -flto)
endif()
if(NMR_EXTRA_CXX_FLAGS)
    list(APPEND _nmr_compile_options ${NMR_EXTRA_CXX_FLAGS})
endif()
target_compile_options(nmr_shielding PRIVATE ${_nmr_compile_options})
target_compile_definitions(nmr_shielding PUBLIC
    NMR_DATA_DIR="${NMR_RUNTIME_DATA_DIR}"
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

# Match HDF5 headers to the selected HDF5 library for every TU in the
# library. Dependency include paths can contribute headers from a different
# HDF5 install, which causes HighFive to emit symbols the selected link
# library may not have and bakes mismatched H5_VERS_* into compilation units.
# Without this prepend, WriteH5 paths can crash at runtime when HDF5 init
# detects the mismatch. Same hazard pattern as the trajectory tests.
target_include_directories(nmr_shielding BEFORE PUBLIC
    ${NMR_HDF5_INCLUDE_DIR}
)

# ============================================================================
# CLI executable
# ============================================================================

add_executable(nmr_extract src/nmr_extract.cpp)
target_link_libraries(nmr_extract PRIVATE nmr_shielding)
# Same selected-HDF5-headers prepend for nmr_extract.cpp itself (the
# RunTrajectory WriteH5 site is in this TU).
target_include_directories(nmr_extract BEFORE PRIVATE
    ${NMR_HDF5_INCLUDE_DIR}
)

add_custom_target(nmr_shielding_doctor
    COMMAND bash "${CMAKE_SOURCE_DIR}/scripts/nmr-shielding-doctor.sh" "$<TARGET_FILE:nmr_extract>"
    DEPENDS nmr_extract
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Run lightweight nmr-shielding producer installability checks")

install(TARGETS nmr_extract
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(PROGRAMS
        scripts/nmr-shielding-doctor.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    RENAME nmr-shielding-doctor)
install(PROGRAMS
        scripts/tensorcs15-check.py
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    RENAME nmr-tensorcs15-check)
install(PROGRAMS
        scripts/run_with_cuda_env.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    RENAME nmr-run-with-cuda-env)
install(FILES
        data/bmrb_atom_nom.tbl
        data/calculator_params.toml
        data/ff14sb_params.dat
    DESTINATION "${NMR_INSTALL_DATA_DIR}")
install(FILES
        data/models/aimnet2_wb97m_0.jpt
    DESTINATION "${NMR_INSTALL_DATA_DIR}/models")
install(DIRECTORY
        data/larsen_hbond_grids
    DESTINATION "${NMR_INSTALL_DATA_DIR}")
install(FILES
        config/nmr_tools.toml.example
    DESTINATION "${CMAKE_INSTALL_SYSCONFDIR}/nmr-shielding")
install(FILES
        config/tensorcs15.manifest.toml
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/nmr-shielding")


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
# enumeration). Labels are suite/tier selectors, not skip policy: tests still
# decide fixture availability at runtime and skip with a reason when optional
# local data is absent.
# ============================================================================

enable_testing()

# GTest for unit tests. Prefer a system/installed GTest; otherwise fetch
# it from GitHub. The fetch needs network at first configure, so it is
# announced (not silent) and gated by NMR_FETCH_GTEST: default ON keeps
# the current behaviour, OFF requires a system GTest (for offline /
# restricted-network configures). This lives after the producer dependency
# closure checks so NMR_PROFILE=producer-full fails on non-portable producer
# paths before considering test-only downloads.
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
nmr_cache_from_env(NMR_NVIDIA_CU13_DIR NMR_NVIDIA_CU13_DIR PATH
    "Python nvidia.cu13 package directory")
nmr_cache_from_env(NMR_NVRTC_LIB_DIR NMR_NVRTC_LIB_DIR PATH
    "Directory containing libnvrtc-builtins.so.13.0")
if(NOT NMR_NVRTC_LIB_DIR AND NMR_NVIDIA_CU13_DIR)
    set(NMR_NVRTC_LIB_DIR "${NMR_NVIDIA_CU13_DIR}/lib"
        CACHE PATH "Directory containing libnvrtc-builtins.so.13.0" FORCE)
endif()
if(NOT NMR_NVRTC_LIB_DIR)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import nvidia.cu13; print(list(nvidia.cu13.__path__)[0])"
        OUTPUT_VARIABLE NMR_NVIDIA_CU13_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _python_cu13_status
    )
    if(_python_cu13_status EQUAL 0 AND NMR_NVIDIA_CU13_DIR)
        set(NMR_NVIDIA_CU13_DIR "${NMR_NVIDIA_CU13_DIR}"
            CACHE PATH "Python nvidia.cu13 package directory" FORCE)
        set(NMR_NVRTC_LIB_DIR "${NMR_NVIDIA_CU13_DIR}/lib"
            CACHE PATH "Directory containing libnvrtc-builtins.so.13.0" FORCE)
    endif()
endif()
if(NOT NMR_NVRTC_LIB_DIR OR NOT EXISTS "${NMR_NVRTC_LIB_DIR}")
    message(FATAL_ERROR
        "Could not locate the bundled PyTorch cu13 directory via "
        "'${Python3_EXECUTABLE} -c \"import nvidia.cu13\"'. Tests that "
        "exercise AIMNet2 will fail at runtime without LD_LIBRARY_PATH set "
        "to this directory's lib/. Set NMR_NVRTC_LIB_DIR directly, set "
        "NMR_NVIDIA_CU13_DIR, or install nvidia-cu13 for the configured "
        "Python interpreter. See "
        "memory entry reference_nvrtc_rpath_fix for the full picture.")
endif()
nmr_require_portable_paths("cu13 nvrtc builtins"
    "${NMR_NVIDIA_CU13_DIR}" "${NMR_NVRTC_LIB_DIR}")
nmr_fail_on_portability_errors()
message(STATUS "Tests: bundled cu13 lib dir = ${NMR_NVRTC_LIB_DIR}")

# Test environment string. Two prepends to LD_LIBRARY_PATH, in order:
#
#  1. NMR_SYSTEM_LIBRARY_DIR so the selected libpq + libssl + libcrypto
#     trio resolves consistently. The binary's RUNPATH can include
#     dependency lib dirs (for libmopac, libreducelib, etc.), and without
#     this prepend the loader can pick libpq/libssl from one environment
#     but libcrypto from the host system -- load fails. We don't USE OpenSSL (DSN is a
#     Unix socket) but libpq's NEEDED entries pull libssl/libcrypto at
#     process startup regardless. Prepending the selected system library
#     directory forces both halves of the SSL pair to the same install.
#
#  2. cu13 lib so libnvrtc.so.13's dlopen of libnvrtc-builtins.so.13.0
#     succeeds at AIMNet2 inference time.
#
# $ENV{LD_LIBRARY_PATH} captures the configure-time shell's existing
# path (orca, system cuda); we don't displace it.
set(_nmr_test_ld_library_path_entries)
if(NMR_SYSTEM_LIBRARY_DIR)
    list(APPEND _nmr_test_ld_library_path_entries "${NMR_SYSTEM_LIBRARY_DIR}")
endif()
list(APPEND _nmr_test_ld_library_path_entries "${NMR_NVRTC_LIB_DIR}")
string(REPLACE ";" ":" _nmr_test_ld_library_path "${_nmr_test_ld_library_path_entries}")
if(NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
    set(_nmr_test_ld_library_path "${_nmr_test_ld_library_path}:$ENV{LD_LIBRARY_PATH}")
endif()
set(NMR_TEST_ENVIRONMENT "LD_LIBRARY_PATH=${_nmr_test_ld_library_path}")

# ----------------------------------------------------------------------------
# nmr_test(<NAME>
#     SOURCES <file.cpp> ...    test sources beyond the universal two
#     [LABELS <label>]          one primary CTest suite/tier label
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
# When LINK contains "highfive", NMR_HDF5_INCLUDE_DIR is prepended to
# the include path so HighFive picks up headers matching the selected
# libhdf5 link.
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
        target_include_directories(${NAME} BEFORE PRIVATE ${NMR_HDF5_INCLUDE_DIR})
    endif()
    target_include_directories(${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/tests)

    target_link_libraries(${NAME} PRIVATE
        nmr_shielding
        GTest::gtest
        ${_NT_LINK}
    )

    list(LENGTH _NT_LABELS _nt_label_count)
    if(_nt_label_count GREATER 1)
        message(FATAL_ERROR
            "nmr_test(${NAME}) accepts one primary CTest label. "
            "Multiple labels do not round-trip reliably through "
            "gtest_discover_tests(... DISCOVERY_MODE PRE_TEST).")
    endif()

    if(_NT_LABELS)
        gtest_discover_tests(${NAME}
            DISCOVERY_MODE PRE_TEST
            PROPERTIES
                ENVIRONMENT "${NMR_TEST_ENVIRONMENT}"
                LABELS "${_NT_LABELS}"
        )
    else()
        gtest_discover_tests(${NAME}
            DISCOVERY_MODE PRE_TEST
            PROPERTIES
                ENVIRONMENT "${NMR_TEST_ENVIRONMENT}"
        )
    endif()
endmacro()


# -- unit_tests -------------------------------------------------------------
# Pure logic; no protein loading, no torch init. <5s.
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
    LABELS conformation
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
    LABELS trajectory
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
        tests/test_dynamics_observables.cpp
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
    LABELS mopac
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
# day-to-day ctest never touches them; run explicitly via `ctest -L batch`.
nmr_test(batch_tests
    LABELS batch
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
    LABELS smoke
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
    PROPERTIES
        ENVIRONMENT "${NMR_TEST_ENVIRONMENT}"
        LABELS fast
)


# ============================================================================
# Example driver -- builds against the library; not a test.
# ============================================================================

add_executable(nmr_pass0_demo tests/pass0_demo.cpp)
target_link_libraries(nmr_pass0_demo PRIVATE nmr_shielding)


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
