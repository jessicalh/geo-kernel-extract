# Build the Qt 6.12.0 Cocoa plugin with the proposed QTBUG-149612 ownership fix.
# Sources and build products stay here; the selected Qt SDK is read-only input.
cmake_minimum_required(VERSION 3.25)
project(H5ReaderCocoaFix VERSION 6.12.0 LANGUAGES C CXX OBJCXX)

if(NOT APPLE OR CMAKE_CROSSCOMPILING)
    message(FATAL_ERROR "This standalone dependency build requires native macOS")
endif()
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install" CACHE PATH "Install prefix" FORCE)
endif()
set(QT_BUILD_INTERNALS_NO_FORCE_SET_INSTALL_PREFIX ON)
find_package(Qt6 6.12.0 EXACT REQUIRED COMPONENTS Core Gui BuildInternals)
include(QtBuild)
find_package(Qt6CorePrivate 6.12.0 EXACT REQUIRED)
find_package(Qt6GuiPrivate 6.12.0 EXACT REQUIRED)
set(QT_GENERATE_SBOM OFF)
set(QT_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(INSTALL_PLUGINSDIR plugins)

set(H5READER_QTBASE_SOURCE_DIR "" CACHE PATH "Unmodified Qt 6.12.0 qtbase sources")
set(_cocoa_source "${H5READER_QTBASE_SOURCE_DIR}/src/plugins/platforms/cocoa")
if(NOT EXISTS "${_cocoa_source}/qcocoaaccessibilityelement.mm")
    message(FATAL_ERROR "Set H5READER_QTBASE_SOURCE_DIR to the installer-provided Qt 6.12.0 qtbase sources")
endif()
file(SHA256 "${_cocoa_source}/qcocoaaccessibilityelement.mm" _source_sha256)
if(NOT _source_sha256 STREQUAL "458a782ea0df93ad9c3c452f7e1b50dd1923fffc41ae8b0fdc15f045e0566ea7")
    message(FATAL_ERROR "Unexpected Cocoa source revision; review the backport for this Qt build")
endif()
get_filename_component(_sdk_prefix "${Qt6_DIR}/../../.." REALPATH)
file(REAL_PATH "${CMAKE_INSTALL_PREFIX}" _install_prefix)
string(FIND "${_install_prefix}/" "${_sdk_prefix}/" _installs_into_sdk)
file(REAL_PATH "${CMAKE_CURRENT_BINARY_DIR}" _binary_dir)
file(REAL_PATH "${H5READER_QTBASE_SOURCE_DIR}" _qt_source_dir)
string(FIND "${_binary_dir}/" "${_qt_source_dir}/" _builds_into_sdk_source)
string(FIND "${_binary_dir}/" "${_sdk_prefix}/" _builds_into_sdk)
string(FIND "${_install_prefix}/" "${_qt_source_dir}/" _installs_into_sdk_source)
if(_installs_into_sdk EQUAL 0 OR _builds_into_sdk_source EQUAL 0
        OR _builds_into_sdk EQUAL 0 OR _installs_into_sdk_source EQUAL 0)
    message(FATAL_ERROR "Use build and install directories outside the Qt SDK")
endif()

set(_copy "${CMAKE_CURRENT_BINARY_DIR}/qt-src")
file(COPY "${_cocoa_source}" DESTINATION "${_copy}/src/plugins/platforms")
file(COPY "${H5READER_QTBASE_SOURCE_DIR}/tests/auto/other/qaccessibilitymac"
    DESTINATION "${_copy}/tests/auto/other")
find_program(_patch_program patch REQUIRED)
execute_process(COMMAND "${_patch_program}" --batch -p1 -i "${CMAKE_CURRENT_SOURCE_DIR}/QTBUG-149612.patch"
    WORKING_DIRECTORY "${_copy}" RESULT_VARIABLE _patch_result
    OUTPUT_VARIABLE _patch_output ERROR_VARIABLE _patch_error)
if(NOT _patch_result EQUAL 0)
    message(FATAL_ERROR "Qt ownership backport did not apply: ${_patch_output}${_patch_error}")
endif()

# Reuse Qt's source list, feature conditions, resources, and plugin build rules.
# Rename only the CMake target to avoid the installed plugin's imported target.
set(_cocoa_cmake "${_copy}/src/plugins/platforms/cocoa/CMakeLists.txt")
file(READ "${_cocoa_cmake}" _cocoa_cmake_text)
string(REPLACE QCocoaIntegrationPlugin H5ReaderCocoaIntegrationPlugin _cocoa_cmake_text "${_cocoa_cmake_text}")
string(REPLACE "qt_internal_add_plugin(H5ReaderCocoaIntegrationPlugin"
    "qt_internal_add_plugin(H5ReaderCocoaIntegrationPlugin\n    SKIP_INSTALL"
    _cocoa_cmake_text "${_cocoa_cmake_text}")
file(WRITE "${_cocoa_cmake}" "${_cocoa_cmake_text}")
foreach(_framework Foundation AppKit Carbon CoreServices CoreVideo IOKit IOSurface Metal QuartzCore UniformTypeIdentifiers ScreenCaptureKit)
    find_library(FW${_framework} ${_framework} REQUIRED)
endforeach()
add_subdirectory("${_copy}/src/plugins/platforms/cocoa" cocoa)
set_target_properties(H5ReaderCocoaIntegrationPlugin PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "@loader_path/../../Frameworks")
install(FILES "$<TARGET_FILE:H5ReaderCocoaIntegrationPlugin>" DESTINATION plugins/platforms)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/QTBUG-149612.patch" DESTINATION share/h5reader-qt-cocoa)

option(H5READER_COCOA_SANITIZERS "Instrument this plugin and its native regression with ASan/UBSan" OFF)
if(H5READER_COCOA_SANITIZERS)
    target_compile_options(H5ReaderCocoaIntegrationPlugin PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(H5ReaderCocoaIntegrationPlugin PRIVATE -fsanitize=address,undefined)
endif()

option(H5READER_COCOA_REGRESSION_TESTS "Build the upstream native Mac accessibility regression" ON)
if(H5READER_COCOA_REGRESSION_TESTS)
    enable_testing()
    find_package(Qt6 6.12.0 EXACT REQUIRED COMPONENTS Widgets Test)
    qt_standard_project_setup()
    qt_add_executable(tst_qaccessibilitymac "${_copy}/tests/auto/other/qaccessibilitymac/tst_qaccessibilitymac.mm")
    target_link_libraries(tst_qaccessibilitymac PRIVATE Qt6::Widgets Qt6::Test
        "-framework AppKit" "-framework ApplicationServices" "-framework Security")
    if(H5READER_COCOA_SANITIZERS)
        target_compile_options(tst_qaccessibilitymac PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(tst_qaccessibilitymac PRIVATE -fsanitize=address,undefined)
    endif()
    add_test(NAME qt_macos_accessibility COMMAND tst_qaccessibilitymac)
    set_tests_properties(qt_macos_accessibility PROPERTIES
        ENVIRONMENT "QT_PLUGIN_PATH=${CMAKE_CURRENT_BINARY_DIR}/plugins" TIMEOUT 55)
endif()
