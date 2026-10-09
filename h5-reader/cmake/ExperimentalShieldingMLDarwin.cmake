# Stage the native CPU helper and unchanged model in standard macOS bundle
# locations. The input keeps the portable infer/model.ts/manifest.json/lib layout.
include_guard(GLOBAL)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    return()
endif()

option(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML
    "Require and package the Experimental Shielding ML model/runtime" ON)
if(NOT H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT
        AND DEFINED ENV{H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT})
    file(TO_CMAKE_PATH "$ENV{H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}"
        H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT)
endif()
set(H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT
    "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}" CACHE PATH
    "Native macOS CPU ML bundle containing infer, model.ts, manifest.json, and lib/")
set(_h5reader_darwin_ml_enabled FALSE)
if(NOT H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT)
    if(H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML)
        message(FATAL_ERROR
            "macOS packages require the Experimental Shielding ML runtime. "
            "Set H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT to a native macOS bundle, "
            "or explicitly set H5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML=OFF for development.")
    endif()
    return()
endif()

foreach(_asset IN ITEMS infer model.ts manifest.json
        lib/libc10.dylib lib/libtorch.dylib lib/libtorch_cpu.dylib)
    set(_asset_path "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/${_asset}")
    if(NOT EXISTS "${_asset_path}" OR IS_DIRECTORY "${_asset_path}")
        message(FATAL_ERROR "The macOS Experimental Shielding ML bundle is missing ${_asset_path}")
    endif()
    file(SIZE "${_asset_path}" _asset_size)
    if(_asset_size EQUAL 0)
        message(FATAL_ERROR "The macOS Experimental Shielding ML asset is empty: ${_asset_path}")
    endif()
endforeach()
find_program(H5READER_LIPO_EXECUTABLE NAMES lipo REQUIRED)
find_program(H5READER_OTOOL_EXECUTABLE NAMES otool REQUIRED)
find_program(H5READER_INSTALL_NAME_TOOL_EXECUTABLE NAMES install_name_tool REQUIRED)
find_program(H5READER_CODESIGN_EXECUTABLE NAMES codesign REQUIRED)
set(_ml_architectures "${CMAKE_OSX_ARCHITECTURES}")
if(NOT _ml_architectures)
    set(_ml_architectures "${CMAKE_SYSTEM_PROCESSOR}")
endif()
file(GLOB _ml_libraries CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE
    "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/lib/*.dylib")
foreach(_binary IN ITEMS "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/infer" ${_ml_libraries})
    execute_process(COMMAND "${H5READER_LIPO_EXECUTABLE}" -verify_arch
        ${_ml_architectures} "${_binary}"
        RESULT_VARIABLE _arch_error ERROR_VARIABLE _arch_message)
    if(_arch_error)
        message(FATAL_ERROR "ML binary does not support ${_ml_architectures}: ${_binary}: ${_arch_message}")
    endif()
endforeach()
execute_process(COMMAND test -x "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/infer"
    RESULT_VARIABLE _helper_not_executable)
if(_helper_not_executable)
    message(FATAL_ERROR "The macOS ML helper is not executable: infer")
endif()
file(READ "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/manifest.json" _ml_manifest)
string(JSON _manifest_type ERROR_VARIABLE _manifest_error TYPE "${_ml_manifest}")
if(_manifest_error OR NOT _manifest_type STREQUAL "OBJECT")
    message(FATAL_ERROR "The ML manifest must contain a JSON object: ${_manifest_error}")
endif()

set(_ml_stage_script [=[
set(bundle "$<TARGET_BUNDLE_DIR:h5reader>/Contents")
set(root "@H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT@")
file(MAKE_DIRECTORY "${bundle}/Helpers" "${bundle}/Frameworks"
    "${bundle}/Resources/ml/experimental_shielding_ml")
file(COPY "${root}/model.ts" "${root}/manifest.json"
    DESTINATION "${bundle}/Resources/ml/experimental_shielding_ml")
file(COPY "${root}/infer" DESTINATION "${bundle}/Helpers")
file(GLOB libraries LIST_DIRECTORIES FALSE "${root}/lib/*.dylib")
file(COPY ${libraries} DESTINATION "${bundle}/Frameworks")
if(IS_DIRECTORY "${root}/licenses")
    file(COPY "${root}/licenses/" DESTINATION "${bundle}/Resources/licenses")
endif()
set(helper "${bundle}/Helpers/infer")
execute_process(COMMAND "@H5READER_OTOOL_EXECUTABLE@" -l "${helper}"
    OUTPUT_VARIABLE load_commands COMMAND_ERROR_IS_FATAL ANY)
if(NOT load_commands MATCHES "path @loader_path/../Frameworks ")
    execute_process(COMMAND "@H5READER_INSTALL_NAME_TOOL_EXECUTABLE@"
        -add_rpath "@loader_path/../Frameworks" "${helper}"
        COMMAND_ERROR_IS_FATAL ANY)
endif()
# Changing Mach-O load commands invalidates its existing signature. Sign the
# staged copy only; the input toolchain and runtime are never modified.
execute_process(COMMAND "@H5READER_CODESIGN_EXECUTABLE@" --force --sign - "${helper}"
    COMMAND_ERROR_IS_FATAL ANY)
]=])
string(CONFIGURE "${_ml_stage_script}" _ml_stage_script @ONLY)
set(_ml_stage_script_path "${CMAKE_CURRENT_BINARY_DIR}/StageDarwinMl-$<CONFIG>.cmake")
file(GENERATE OUTPUT "${_ml_stage_script_path}" CONTENT "${_ml_stage_script}")
file(GLOB_RECURSE _ml_assets CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE
    "${H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT}/*")
set(_ml_stamp "${CMAKE_CURRENT_BINARY_DIR}/h5reader-darwin-ml-$<CONFIG>.stamp")
add_custom_command(OUTPUT "${_ml_stamp}"
    COMMAND "${CMAKE_COMMAND}" -P "${_ml_stage_script_path}"
    COMMAND "${CMAKE_COMMAND}" -E touch "${_ml_stamp}"
    DEPENDS ${_ml_assets} "${_ml_stage_script_path}"
    COMMENT "Staging native macOS Experimental Shielding ML runtime"
    VERBATIM)
add_custom_target(h5reader_deploy_experimental_shielding_ml DEPENDS "${_ml_stamp}")
add_dependencies(h5reader h5reader_deploy_experimental_shielding_ml)
set(_h5reader_darwin_ml_enabled TRUE)
