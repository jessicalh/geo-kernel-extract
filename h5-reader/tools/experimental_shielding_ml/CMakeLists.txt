cmake_minimum_required(VERSION 3.24)
project(h5reader_experimental_shielding_infer LANGUAGES CXX)

find_package(Torch REQUIRED)

if(MSVC AND TARGET torch_hip)
    set_property(TARGET torch_hip PROPERTY INTERFACE_COMPILE_OPTIONS
        "$<$<COMPILE_LANGUAGE:CXX>:/permissive->"
        /EHsc
        /bigobj
        -D__HIP_PLATFORM_AMD__=1
        -DCUDA_HAS_FP16=1
        -DUSE_ROCM
        -D__HIP_NO_HALF_OPERATORS__=1
        -D__HIP_NO_HALF_CONVERSIONS__=1
        -DTORCH_HIP_VERSION=702
        -DCAFFE2_USE_MIOPEN
        -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_HIP
        -DHIPBLAS_V2
        -DHIP_ENABLE_WARP_SYNC_BUILTINS)

    set(_exception_source
        "${TORCH_INSTALL_PREFIX}/include/c10/util/Exception.h")
    if(NOT EXISTS "${_exception_source}")
        message(FATAL_ERROR
            "ROCm Exception header not found: ${_exception_source}")
    endif()
    file(READ "${_exception_source}" _exception_text)
    set(_value_error_declaration
        "class C10_API ValueError : public Error {\n  using Error::Error;\n};")
    string(FIND "${_exception_text}" "${_value_error_declaration}"
        _value_error_offset)
    if(_value_error_offset EQUAL -1)
        message(FATAL_ERROR
            "ROCm ValueError compatibility point changed; review the upstream header")
    endif()
    set(_value_error_compat
        "class ValueError : public Error {\n public:\n  ValueError(SourceLocation source_location, std::string msg)\n      : Error(source_location, std::move(msg)) {}\n  using Error::Error;\n};")
    string(REPLACE "${_value_error_declaration}" "${_value_error_compat}"
        _exception_text "${_exception_text}")
    set(_torch_compat_include
        "${CMAKE_CURRENT_BINARY_DIR}/torch-compat/include")
    file(MAKE_DIRECTORY "${_torch_compat_include}/c10/util")
    file(WRITE
        "${_torch_compat_include}/c10/util/Exception.h"
        "${_exception_text}")
endif()

add_executable(infer infer.cpp)
target_compile_features(infer PRIVATE cxx_std_17)
target_link_libraries(infer PRIVATE "${TORCH_LIBRARIES}")

if(_torch_compat_include)
    target_include_directories(infer BEFORE PRIVATE "${_torch_compat_include}")
endif()

if(MSVC)
    target_compile_options(infer PRIVATE /EHsc /bigobj /wd4267)
endif()
