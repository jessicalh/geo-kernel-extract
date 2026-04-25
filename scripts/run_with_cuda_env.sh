#!/bin/bash
#
# run_with_cuda_env.sh — universal test / binary launcher.
#
# Resolves PyTorch's bundled cu13 lib directory and prepends it to
# LD_LIBRARY_PATH before exec'ing the wrapped command. Required for
# any binary that touches libtorch CUDA JIT — AIMNet2 in this
# project. Without this, libnvrtc-builtins.so.13.0 fails to load and
# the call site throws an nvrtc compile error.
#
# Why this is a wrapper, not CMAKE_BUILD_RPATH:
#   When libtorch_cuda dlopens libnvrtc.so.13, libnvrtc.so.13 in
#   turn dlopens libnvrtc-builtins.so.13.0 from a search that does
#   NOT consult the calling executable's RPATH — only the loading
#   library's own RPATH or LD_LIBRARY_PATH. The bundled libnvrtc.so
#   .13 has empty RPATH, so LD_LIBRARY_PATH is the actual path.
#
# Why a wrapper, not a setenv() in C++ startup:
#   Setting environment variables from inside a process does not
#   reliably affect the dynamic loader's library search path; the
#   loader caches resolution at process start. Setting the env
#   *before* exec'ing the binary is the path that actually works.
#
# Why a discovery script and not hardcoded paths:
#   PyTorch is pip-installed under the user's site-packages; its
#   exact path varies by user and machine. Hardcoding it would
#   break deployment to other batcave-class machines. Resolving via
#   Python at launch is one location-sensitive lookup.
#
# Usage:
#   scripts/run_with_cuda_env.sh ./build/some_test [args...]
#   scripts/run_with_cuda_env.sh ./build/nmr_extract --trajectory ...
#
# CMake registers this as TEST_LAUNCHER for every gtest_discover_tests
# target so `ctest` invokes every test through it. Direct
# `./build/<binary>` invocation is NOT the documented path.

set -e

# Resolve the bundled cu13 directory via Python. Returns the package
# directory; the libraries we need are in the lib/ subdirectory.
CU13_DIR=$(python3 -c "import nvidia.cu13; print(list(nvidia.cu13.__path__)[0])" 2>/dev/null)

if [ -n "$CU13_DIR" ] && [ -d "${CU13_DIR}/lib" ]; then
    export LD_LIBRARY_PATH="${CU13_DIR}/lib:${LD_LIBRARY_PATH}"
else
    echo "run_with_cuda_env: nvidia.cu13 not found via Python; AIMNet2 / libtorch CUDA JIT will fail." >&2
    echo "run_with_cuda_env: install via 'pip install --user nvidia-cu13' or equivalent." >&2
fi

exec "$@"
