# Native Linux runtime bundle. Include after install(TARGETS h5reader ...).
# Qt owns plugin discovery/deployment; the launcher owns the process environment.
include_guard(GLOBAL)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    return()
endif()
if(CMAKE_CROSSCOMPILING)
    message(FATAL_ERROR "Linux runtime deployment requires a native Linux build")
endif()
if(Qt6_VERSION VERSION_LESS 6.10)
    message(FATAL_ERROR
        "Linux runtime deployment requires Qt 6.10 or newer "
        "(the deployment plugin selection API was added in Qt 6.10).")
endif()
if(NOT TARGET h5reader OR NOT COMMAND qt6_generate_deploy_script)
    message(FATAL_ERROR "DeployLinux requires the h5reader target and Qt deployment APIs")
endif()

set_target_properties(h5reader PROPERTIES
    INSTALL_RPATH "$ORIGIN/lib"
    INSTALL_RPATH_USE_LINK_PATH FALSE
    BUILD_RPATH_USE_ORIGIN TRUE)

option(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML
    "Require and package the Experimental Shielding ML model/runtime" ON)
if(NOT H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT
        AND DEFINED ENV{H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT})
    file(TO_CMAKE_PATH "$ENV{H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}"
        H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT)
endif()
set(H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT
    "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}" CACHE PATH
    "Portable Linux ML bundle containing infer, model.ts, manifest.json, and lib/")
if(DEFINED H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR
        AND NOT H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR STREQUAL "ml/experimental_shielding_ml")
    message(FATAL_ERROR
        "Linux runtime discovery requires H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR=ml/experimental_shielding_ml")
endif()
set(H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR "ml/experimental_shielding_ml")

set(_h5reader_linux_ml_enabled FALSE)
set(_h5reader_linux_ml_libraries "")
if(H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT)
    foreach(_asset IN ITEMS infer model.ts manifest.json
            lib/libc10.so lib/libtorch.so lib/libtorch_cpu.so lib/libtorch_global_deps.so)
        set(_asset_path "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/${_asset}")
        if(NOT EXISTS "${_asset_path}" OR IS_DIRECTORY "${_asset_path}")
            message(FATAL_ERROR "The Linux Experimental Shielding ML bundle is missing ${_asset_path}")
        endif()
        file(SIZE "${_asset_path}" _asset_size)
        if(_asset_size EQUAL 0)
            message(FATAL_ERROR "The Linux Experimental Shielding ML asset is empty: ${_asset_path}")
        endif()
    endforeach()
    file(READ "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/infer"
        _helper_magic LIMIT 4 HEX)
    if(NOT _helper_magic STREQUAL "7f454c46")
        message(FATAL_ERROR "The Linux ML helper must be a native ELF executable: infer")
    endif()
    find_program(_h5reader_test_program NAMES test REQUIRED)
    execute_process(COMMAND "${_h5reader_test_program}" -x
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/infer"
        RESULT_VARIABLE _helper_not_executable)
    if(_helper_not_executable)
        message(FATAL_ERROR "The Linux ML helper is not executable; run chmod +x on infer")
    endif()
    file(READ "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/manifest.json" _manifest)
    string(JSON _manifest_type ERROR_VARIABLE _manifest_error TYPE "${_manifest}")
    if(_manifest_error OR NOT _manifest_type STREQUAL "OBJECT")
        message(FATAL_ERROR "The ML manifest must contain a JSON object: ${_manifest_error}")
    endif()
    unset(_manifest)

    # Copy the entire dedicated bundle, including all Torch/OpenMP shared objects
    # and license notices. A stamp avoids recopying hundreds of MB on each build.
    file(GLOB_RECURSE _h5reader_linux_ml_assets CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/*")
    file(GLOB _h5reader_linux_ml_libraries CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE
        "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/lib/*.so*")
    set(_h5reader_linux_ml_stamp
        "${CMAKE_CURRENT_BINARY_DIR}/h5reader-experimental-shielding-ml-$<CONFIG>.stamp")
    add_custom_command(OUTPUT "${_h5reader_linux_ml_stamp}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory
            "$<TARGET_FILE_DIR:h5reader>/${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}"
        COMMAND "${CMAKE_COMMAND}" -E copy_directory
            "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}"
            "$<TARGET_FILE_DIR:h5reader>/${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}"
        COMMAND "${CMAKE_COMMAND}" -E touch "${_h5reader_linux_ml_stamp}"
        DEPENDS ${_h5reader_linux_ml_assets}
        COMMENT "Staging Linux Experimental Shielding ML runtime"
        VERBATIM)
    add_custom_target(h5reader_deploy_experimental_shielding_ml
        DEPENDS "${_h5reader_linux_ml_stamp}")
    add_dependencies(h5reader h5reader_deploy_experimental_shielding_ml)
    install(DIRECTORY "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/"
        DESTINATION "${H5READER_EXPERIMENTAL_SHIELDING_ML_INSTALL_DIR}"
        COMPONENT Runtime USE_SOURCE_PERMISSIONS)
    set(_h5reader_linux_ml_enabled TRUE)
elseif(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML)
    message(FATAL_ERROR
        "Linux packages require the Experimental Shielding ML runtime. "
        "Set H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT to a portable Linux bundle, "
        "or explicitly set H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML=OFF for development.")
endif()

# The OpenSSL backend and Qt's FFmpeg stubs can load these libraries dynamically,
# so include them even when no executable DT_NEEDED entry mentions them.
find_package(OpenSSL 3 REQUIRED)
find_program(H5READER_PATCHELF_EXECUTABLE NAMES patchelf REQUIRED
    DOC "patchelf used to remove development paths from Linux runtime bundles")
configure_file("${CMAKE_CURRENT_LIST_DIR}/InstallLinuxRuntime.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/InstallLinuxRuntime.in.cmake" @ONLY)
file(READ "${CMAKE_CURRENT_BINARY_DIR}/InstallLinuxRuntime.in.cmake" _h5reader_deploy_content)
# The versioned function preserves the script's semicolon lists and variables;
# the versionless macro expands ARGV before forwarding it on Qt 6.10.
qt6_generate_deploy_script(TARGET h5reader OUTPUT_SCRIPT _h5reader_linux_deploy_script
    CONTENT "${_h5reader_deploy_content}")
install(SCRIPT "${_h5reader_linux_deploy_script}" COMPONENT Runtime)

configure_file("${CMAKE_CURRENT_LIST_DIR}/LaunchLinux.sh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/launch-h5reader" @ONLY)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/launch-h5reader"
    DESTINATION . COMPONENT Runtime)

# Upstream notices are separate from the application's own licensing decision.
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../extern/HighFive/LICENSE"
    "${CMAKE_CURRENT_SOURCE_DIR}/../extern/HighFive/AUTHORS.txt"
    DESTINATION licenses/HighFive COMPONENT Runtime OPTIONAL)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/extern/nanoflann/COPYING"
    DESTINATION licenses/nanoflann COMPONENT Runtime OPTIONAL)
get_filename_component(_h5reader_vtk_prefix "${VTK_DIR}/../../.." ABSOLUTE)
if(IS_DIRECTORY "${_h5reader_vtk_prefix}/share/licenses/VTK")
    install(DIRECTORY "${_h5reader_vtk_prefix}/share/licenses/VTK/"
        DESTINATION licenses/VTK COMPONENT Runtime)
endif()
set(H5READER_THIRD_PARTY_NOTICES_DIR "" CACHE PATH
    "Prepared upstream license texts, notices and provenance for the Linux bundle")
if(H5READER_THIRD_PARTY_NOTICES_DIR)
    if(NOT IS_DIRECTORY "${H5READER_THIRD_PARTY_NOTICES_DIR}")
        message(FATAL_ERROR "H5READER_THIRD_PARTY_NOTICES_DIR must name an existing directory")
    endif()
    install(DIRECTORY "${H5READER_THIRD_PARTY_NOTICES_DIR}/"
        DESTINATION licenses/upstream COMPONENT Runtime)
endif()
