cmake_minimum_required(VERSION 3.24)
project(ScienceFilesQt LANGUAGES CXX)

find_package(Qt6 6.8 REQUIRED COMPONENTS Core Network)
find_package(LibArchive REQUIRED)

add_library(science_files_qt STATIC
    include/sciencefiles/Download.h
    include/sciencefiles/BundleCache.h
    src/Download.cpp
    src/BundleCache.cpp
    src/ArchiveTask.cpp src/ArchiveTask.h
    src/BundlePaths.cpp src/BundlePaths.h)
add_library(ScienceFiles::Qt ALIAS science_files_qt)
set_target_properties(science_files_qt PROPERTIES AUTOMOC ON)
target_compile_features(science_files_qt PUBLIC cxx_std_17 PRIVATE cxx_std_20)
target_include_directories(science_files_qt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(science_files_qt PUBLIC Qt6::Core Qt6::Network
    PRIVATE LibArchive::LibArchive)
if(MSVC)
    target_compile_options(science_files_qt PRIVATE /W4 /WX /permissive- /utf-8)
else()
    target_compile_options(science_files_qt PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
