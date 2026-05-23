# BuildType-Debug.cmake — sanitizer + assertion flags for Debug builds.
#
# Linux + macOS: AddressSanitizer + UndefinedBehaviorSanitizer, glibc
# assertions on, frame pointers preserved. Compile + link flags must
# match or the linker drops the sanitizer runtime.
#
# Windows MSVC: AddressSanitizer is supported on MSVC 17.0+ but has
# real runtime caveats (only /MD works, conflicts with iterator
# debugging, requires the ASan DLL on PATH). Opt-in via
# -DH5READER_MSVC_ASAN=ON because most Debug usage on Windows is just
# stepping through with the debugger, not sanitizer-driven.

include_guard(GLOBAL)

function(h5reader_apply_build_type_settings target)
    if(MSVC)
        if(H5READER_MSVC_ASAN)
            target_compile_options(${target} PRIVATE /fsanitize=address)
            target_compile_definitions(${target} PRIVATE
                _DISABLE_VECTOR_ANNOTATION=1
                _DISABLE_STRING_ANNOTATION=1)
            message(STATUS
                "h5reader: MSVC AddressSanitizer enabled. "
                "Ensure clang_rt.asan_dynamic-x86_64.dll is on PATH.")
        else()
            message(STATUS
                "h5reader: Debug build; MSVC ASan off by default "
                "(pass -DH5READER_MSVC_ASAN=ON to enable).")
        endif()
    else()
        # GCC / Clang / AppleClang. -O1 gives readable stacks while
        # keeping sanitizer warnings actionable; -O0 generates such
        # bloated output that UBSan trips on stack-protector noise.
        target_compile_options(${target} PRIVATE
            -O1
            -g3
            -fno-omit-frame-pointer
            -fsanitize=address
            -fsanitize=undefined)
        target_link_options(${target} PRIVATE
            -fsanitize=address
            -fsanitize=undefined)
        target_compile_definitions(${target} PRIVATE
            _GLIBCXX_ASSERTIONS=1)
        message(STATUS
            "h5reader: Debug build; ASan + UBSan enabled "
            "(LD_PRELOAD not needed for static-linked sanitizers).")
    endif()
endfunction()
