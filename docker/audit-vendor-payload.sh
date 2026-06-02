#!/usr/bin/env bash
# Validate the curated Docker vendor payload before building deps-vendor.Dockerfile.

set -u

vendor_dir=${1:-docker/vendor}
deps_dir=${vendor_dir%/}/deps
tensor_dir=${vendor_dir%/}/tensorcs15
failures=0
warnings=0

ok() { printf 'OK    %s\n' "$*"; }
warn() { warnings=$((warnings + 1)); printf 'WARN  %s\n' "$*" >&2; }
fail() { failures=$((failures + 1)); printf 'FAIL  %s\n' "$*" >&2; }

check_file() {
    local path=$1 label=${2:-$1}
    if [ -f "$path" ]; then ok "$label"; else fail "$label missing: $path"; fi
}

check_exec() {
    local path=$1 label=${2:-$1}
    if [ -x "$path" ]; then ok "$label"; else fail "$label missing/not executable: $path"; fi
}

check_dir() {
    local path=$1 label=${2:-$1}
    if [ -d "$path" ]; then ok "$label"; else fail "$label missing: $path"; fi
}

check_glob() {
    local pattern=$1 label=$2
    local matches=()
    shopt -s nullglob
    matches=( $pattern )
    shopt -u nullglob
    if [ "${#matches[@]}" -gt 0 ]; then ok "$label (${#matches[@]})"; else fail "$label missing: $pattern"; fi
}

check_dir "$deps_dir" "vendor deps dir"

symlink_report=$(find "$vendor_dir" -type l -print 2>/dev/null | head -20)
if [ -n "$symlink_report" ]; then
    fail "vendor payload contains symlinks/links; copy with docker/collect-local-vendor-deps.sh --copy to dereference them"
    printf '%s\n' "$symlink_report" >&2
else
    ok "vendor payload has no symlinks"
fi

check_exec "$deps_dir/gromacs/bin/gmx_mpi" "GROMACS gmx_mpi executable"
check_dir "$deps_dir/gromacs/share/gromacs/top" "GROMACS topology/share data"
check_file "$deps_dir/gromacs/lib/libgromacs_mpi.so" "GROMACS libgromacs_mpi"
check_file "$deps_dir/gromacs/lib/libmuparser.so.2" "GROMACS muparser runtime library"
check_file "$deps_dir/gromacs-src/CMakeLists.txt" "GROMACS full source root"
check_dir "$deps_dir/gromacs-src/src/gromacs" "GROMACS source internal headers"
check_file "$deps_dir/gromacs-src/api/legacy/include/gromacs/utility/basedefinitions.h" "GROMACS legacy API headers"
check_file "$deps_dir/gromacs-build/api/legacy/include/gromacs/libgromacs_export.h" "GROMACS generated libgromacs export header"
check_file "$deps_dir/gromacs-build/api/legacy/include/gromacs/version.h" "GROMACS generated version header"
intermediate_report=$(find "$deps_dir/gromacs-build" \( -name CMakeCache.txt -o -name CMakeFiles -o -name Makefile -o -name '*.o' -o -name '*.a' \) -print -quit 2>/dev/null || true)
if [ -n "$intermediate_report" ]; then
    fail "GROMACS generated-header payload contains compiler/build intermediates: $intermediate_report"
else
    ok "GROMACS generated-header payload excludes compiler/build intermediates"
fi
check_exec "$deps_dir/orca/orca" "ORCA executable"
check_dir "$deps_dir/orca/lib" "ORCA library directory"
check_file "$deps_dir/reduce-src/build/reduce_src/libreducelib.a" "reduce reducelib archive"
check_file "$deps_dir/reduce-src/build/libpdb/libpdb++.a" "reduce libpdb archive"
check_file "$deps_dir/reduce-src/build/toolclasses/libtoolclasses.a" "reduce toolclasses archive"
check_file "$deps_dir/reduce-src/reduce_wwPDB_het_dict.txt" "reduce het dictionary"
check_file "$deps_dir/chem-env/include/mopac.h" "MOPAC header"
check_glob "$deps_dir/chem-env/lib/libmopac.so*" "MOPAC shared library"
check_exec "$deps_dir/chem-env/bin/mopac" "MOPAC executable"
check_exec "$deps_dir/chem-env/bin/tleap" "AmberTools tleap"
check_file "$deps_dir/torch/share/cmake/Torch/TorchConfig.cmake" "Torch CMake package"
check_glob "$deps_dir/torch/lib/libtorch*.so*" "Torch shared libraries"
check_file "$deps_dir/nvidia/cu13/lib/libnvrtc-builtins.so.13.0" "CUDA 13.0 nvrtc builtins"
check_glob "$deps_dir/nvidia/cu13/lib/libcudart.so*" "CUDA runtime library"
check_glob "$deps_dir/nvidia/cudnn/lib/libcudnn*.so*" "cuDNN libraries"
check_glob "$deps_dir/nvidia/cusparselt/lib/libcusparseLt*.so*" "cuSPARSELt libraries"
check_glob "$deps_dir/nvidia/nccl/lib/libnccl*.so*" "NCCL libraries"
check_glob "$deps_dir/nvidia/nvshmem/lib/libnvshmem*.so*" "NVSHMEM libraries"

check_file "$tensor_dir/tensorcs15.dump" "tensorcs15 PostgreSQL dump"

if [ -f "$vendor_dir/manifest.toml" ]; then
    ok "vendor manifest: $vendor_dir/manifest.toml"
else
    warn "vendor manifest absent: $vendor_dir/manifest.toml"
fi

printf 'SUMMARY failures=%d warnings=%d\n' "$failures" "$warnings"
[ "$failures" -eq 0 ]
