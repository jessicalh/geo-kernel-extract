# Keep macOS packaging independent from Linux ELF deployment. Qt's deployment
# API invokes the selected Qt SDK's macdeployqt, including its code signing.
include_guard(GLOBAL)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    return()
endif()
if(CMAKE_CROSSCOMPILING)
    message(FATAL_ERROR "macOS runtime deployment requires a native macOS build")
endif()
get_target_property(_reader_is_bundle h5reader MACOSX_BUNDLE)
if(NOT _reader_is_bundle)
    message(FATAL_ERROR "macOS runtime deployment requires an h5reader app bundle")
endif()

include(ExperimentalShieldingMLDarwin)
set_target_properties(h5reader PROPERTIES
    INSTALL_RPATH "@executable_path/../Frameworks"
    INSTALL_RPATH_USE_LINK_PATH FALSE)

# Optional, independently built Qt 6.12.0 Cocoa ownership backport. Deploy into
# a clean staging prefix; NO_OVERWRITE keeps macdeployqt from replacing it with
# the SDK plugin. Other platforms and ordinary SDK deployment are unchanged.
set(H5READER_MACOS_COCOA_PLUGIN "" CACHE FILEPATH
    "Qt 6.12.0 Cocoa plugin built by tools/macos/qt-cocoa with QTBUG-149612 fixed")
set(_darwin_cocoa_copy "")
set(_darwin_cocoa_no_overwrite "")
if(H5READER_MACOS_COCOA_PLUGIN)
    if(NOT Qt6_VERSION VERSION_EQUAL "6.12.0")
        message(FATAL_ERROR "The Cocoa ownership backport must match Qt 6.12.0 exactly")
    endif()
    if(NOT EXISTS "${H5READER_MACOS_COCOA_PLUGIN}")
        message(FATAL_ERROR "Build H5READER_MACOS_COCOA_PLUGIN before deploying Reader")
    endif()
    get_filename_component(H5READER_MACOS_COCOA_PLUGIN "${H5READER_MACOS_COCOA_PLUGIN}" ABSOLUTE)
    set(_darwin_cocoa_copy
        "file(MAKE_DIRECTORY \"\${_bundle}/Contents/PlugIns/platforms\")\nfile(COPY \"${H5READER_MACOS_COCOA_PLUGIN}\" DESTINATION \"\${_bundle}/Contents/PlugIns/platforms\")\n")
    set(_darwin_cocoa_no_overwrite "    NO_OVERWRITE\n")
endif()

# CMake removes development RPATHs from the installed executable. Discover the
# native dylib closure from the build binaries while their original RPATHs are
# available; macdeployqt's -libpath does not resolve these third-party @rpaths.
set(_darwin_helper_argument "")
set(_darwin_source_helper_argument "")
if(_h5reader_darwin_ml_enabled)
    set(_darwin_helper_argument
        "    ADDITIONAL_EXECUTABLES \"\${QT_DEPLOY_PREFIX}/$<TARGET_BUNDLE_DIR_NAME:h5reader>/Contents/Helpers/infer\"\n")
    set(_darwin_source_helper_argument
        "\"$<TARGET_BUNDLE_DIR:h5reader>/Contents/Helpers/infer\"")
endif()
qt6_generate_deploy_script(TARGET h5reader OUTPUT_SCRIPT _darwin_deploy_script
    CONTENT "
set(_bundle \"\${QT_DEPLOY_PREFIX}/$<TARGET_BUNDLE_DIR_NAME:h5reader>\")
file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES \"$<TARGET_FILE:h5reader>\" ${_darwin_source_helper_argument}
    RESOLVED_DEPENDENCIES_VAR _native_resolved
    UNRESOLVED_DEPENDENCIES_VAR _native_unresolved
    PRE_EXCLUDE_REGEXES \"^/System/Library/\" \"^/usr/lib/\"
    POST_EXCLUDE_REGEXES \"^/System/Library/\" \"^/usr/lib/\")
if(_native_unresolved)
    message(FATAL_ERROR \"Unresolved build runtime dependencies: \${_native_unresolved}\")
endif()
foreach(_dependency IN LISTS _native_resolved)
    if(_dependency MATCHES \"[.]dylib$\")
        file(COPY \"\${_dependency}\" DESTINATION \"\${_bundle}/Contents/Frameworks\"
            FOLLOW_SYMLINK_CHAIN)
    endif()
endforeach()
${_darwin_cocoa_copy}
qt6_deploy_runtime_dependencies(
    EXECUTABLE \"$<TARGET_BUNDLE_DIR_NAME:h5reader>\"
${_darwin_helper_argument}${_darwin_cocoa_no_overwrite}    DEPLOY_TOOL_OPTIONS -no-strip)

set(_main_executable \"\${_bundle}/Contents/MacOS/$<TARGET_FILE_NAME:h5reader>\")
set(_executables \"\${_main_executable}\")
if(EXISTS \"\${_bundle}/Contents/Helpers/infer\")
    list(APPEND _executables \"\${_bundle}/Contents/Helpers/infer\")
endif()
file(GLOB_RECURSE _plugins LIST_DIRECTORIES FALSE
    \"\${_bundle}/Contents/PlugIns/*.dylib\")
file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES \${_executables}
    MODULES \${_plugins}
    BUNDLE_EXECUTABLE \"\${_main_executable}\"
    RESOLVED_DEPENDENCIES_VAR _resolved
    UNRESOLVED_DEPENDENCIES_VAR _unresolved
    PRE_EXCLUDE_REGEXES \"^/System/Library/\" \"^/usr/lib/\"
    POST_EXCLUDE_REGEXES \"^/System/Library/\" \"^/usr/lib/\")
if(_unresolved)
    message(FATAL_ERROR \"Unresolved macOS runtime dependencies: \${_unresolved}\")
endif()
file(REAL_PATH \"\${_bundle}\" _bundle_real)
foreach(_dependency IN LISTS _resolved)
    file(REAL_PATH \"\${_dependency}\" _dependency_real)
    string(FIND \"\${_dependency_real}\" \"\${_bundle_real}/\" _in_bundle)
    if(NOT _in_bundle EQUAL 0)
        message(FATAL_ERROR \"macOS runtime dependency remains outside the bundle: \${_dependency}\")
    endif()
endforeach()
message(STATUS \"macOS runtime bundle validated at \${_bundle}\")
")
set(_darwin_notices "$<TARGET_BUNDLE_DIR_NAME:h5reader>/Contents/Resources/licenses")
if(H5READER_MACOS_COCOA_PLUGIN)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/tools/macos/qt-cocoa/QTBUG-149612.patch"
        DESTINATION "${_darwin_notices}/Qt-Cocoa-backport" COMPONENT Runtime)
endif()
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../extern/HighFive/LICENSE"
    "${CMAKE_CURRENT_SOURCE_DIR}/../extern/HighFive/AUTHORS.txt"
    DESTINATION "${_darwin_notices}/HighFive" COMPONENT Runtime OPTIONAL)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/extern/nanoflann/COPYING"
    DESTINATION "${_darwin_notices}/nanoflann" COMPONENT Runtime OPTIONAL)
set(H5READER_THIRD_PARTY_NOTICES_DIR "" CACHE PATH
    "Prepared upstream license texts, notices and provenance for the macOS bundle")
if(H5READER_THIRD_PARTY_NOTICES_DIR)
    if(NOT IS_DIRECTORY "${H5READER_THIRD_PARTY_NOTICES_DIR}")
        message(FATAL_ERROR "H5READER_THIRD_PARTY_NOTICES_DIR must name an existing directory")
    endif()
    install(DIRECTORY "${H5READER_THIRD_PARTY_NOTICES_DIR}/"
        DESTINATION "${_darwin_notices}/upstream" COMPONENT Runtime)
endif()
# Deploy and sign after every resource has been installed into the bundle.
install(SCRIPT "${_darwin_deploy_script}" COMPONENT Runtime)
