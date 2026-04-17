cmake_minimum_required(VERSION 3.16)
project(analysis_fileformat CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# HighFive — header-only, lives in the main project's extern/
set(HIGHFIVE_DIR "${CMAKE_SOURCE_DIR}/../extern/HighFive/include")

# HDF5 C library (system)
find_package(HDF5 REQUIRED COMPONENTS C)

# Library
add_library(analysis_file analysis_file.cpp)
target_include_directories(analysis_file PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${HIGHFIVE_DIR}
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(analysis_file PUBLIC ${HDF5_C_LIBRARIES})

# Roundtrip test
add_executable(roundtrip_test roundtrip_test.cpp)
target_link_libraries(roundtrip_test PRIVATE analysis_file)

# Roundtrip writer (generates persistent roundtripped H5s)
add_executable(roundtrip_write roundtrip_write.cpp)
target_link_libraries(roundtrip_write PRIVATE analysis_file)
