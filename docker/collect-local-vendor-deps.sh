#!/usr/bin/env bash
# Build a local-to-vendor copy plan for the producer Docker dependency payload.
# Default mode is dry-run: it prints what would be copied and reports sizes.
# Use --copy only after reviewing the plan and disk space.

set -euo pipefail

cache=build/CMakeCache.txt
vendor_dir=docker/vendor
copy=0

usage() {
    cat <<'USAGE'
usage: docker/collect-local-vendor-deps.sh [--copy] [--cache build/CMakeCache.txt] [--vendor-dir docker/vendor]

Default is dry-run. --copy dereferences links while populating docker/vendor/deps from the
current local CMake/Torch environment.
USAGE
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --copy) copy=1 ;;
        --cache) cache=$2; shift ;;
        --vendor-dir) vendor_dir=$2; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

if [ ! -f "$cache" ]; then
    echo "missing CMake cache: $cache" >&2
    exit 2
fi

cache_get() {
    local key=$1
    sed -n "s#^${key}:[^=]*=##p" "$cache" | tail -n 1
}

py_path() {
    python3 - "$1" <<'PY'
import pathlib
import sys
kind = sys.argv[1]
if kind == "torch":
    import torch
    print(pathlib.Path(torch.__file__).resolve().parent)
elif kind == "nvidia":
    import nvidia
    print(pathlib.Path(list(nvidia.__path__)[0]).resolve())
elif kind == "torch_version":
    import torch
    print(torch.__version__)
elif kind == "torch_cuda":
    import torch
    print(torch.version.cuda)
else:
    raise SystemExit(2)
PY
}

add_plan() {
    local src=$1 dst=$2 label=$3
    printf '%s\t%s\t%s\n' "$src" "$dst" "$label" >>"$plan_file"
}

copy_path() {
    local src=$1 dst=$2 label=$3
    if [ ! -e "$src" ]; then
        echo "MISSING $label: $src" >&2
        return 1
    fi
    if [ "$copy" -eq 1 ]; then
        if [ "${dst%/}" != "$dst" ]; then
            mkdir -p "$dst"
            rsync -aL "$src" "$dst"
        else
            mkdir -p "$(dirname "$dst")"
            if [ -f "$src" ]; then
                cp -L "$src" "$dst"
            else
                rsync -aL "$src" "$dst"
            fi
        fi
    else
        local size_src=$src
        if [ -f "$src" ]; then
            size_src=$(realpath -e "$src")
        fi
        if command -v du >/dev/null 2>&1; then
            du -shL "$size_src" 2>/dev/null | awk -v label="$label" -v dst="$dst" '{print "DRY   " label ": " $1 " -> " dst}'
        else
            echo "DRY   $label -> $dst"
        fi
    fi
}

vendor_dir=${vendor_dir%/}
deps_dir=$vendor_dir/deps
plan_file=$(mktemp)
trap 'rm -f "$plan_file"' EXIT

GROMACS_LIB=$(cache_get GROMACS_LIB)
GROMACS_SRC=$(cache_get GROMACS_SRC)
GROMACS_BUILD=$(cache_get GROMACS_BUILD)
GROMACS_SOURCE_ROOT=$(dirname "$GROMACS_SRC")
GROMACS_RUNTIME_PREFIX=${NMR_GROMACS_RUNTIME_PREFIX:-/tmp/nmr-gromacs-install}
if [ ! -x "$GROMACS_RUNTIME_PREFIX/bin/gmx_mpi" ]; then
    GROMACS_RUNTIME_PREFIX=$(dirname "$(dirname "$GROMACS_LIB")")
fi
ORCA_ROOT=${NMR_ORCA_ROOT:-/opt/orca}
REDUCE_SRC=$(cache_get REDUCE_SRC)
MOPAC_INCLUDE=$(cache_get MOPAC_INCLUDE)
MOPAC_LIB=$(cache_get MOPAC_LIB)

torch_root=$(py_path torch)
nvidia_root=$(py_path nvidia)
torch_version=$(py_path torch_version)
torch_cuda=$(py_path torch_cuda)
torch_cuda_re=${torch_cuda//./\.}
gromacs_version=$("$GROMACS_RUNTIME_PREFIX/bin/gmx_mpi" --version 2>/dev/null || true)
if [ -z "$gromacs_version" ]; then
    echo "unable to run GROMACS runtime probe: $GROMACS_RUNTIME_PREFIX/bin/gmx_mpi --version" >&2
    exit 2
fi
if ! printf '%s\n' "$gromacs_version" | grep -Eq "CUDA runtime:[[:space:]]*${torch_cuda_re}($|[^0-9])"; then
    echo "GROMACS CUDA runtime does not match Torch CUDA ${torch_cuda}; set NMR_GROMACS_RUNTIME_PREFIX to the CUDA-matched GROMACS install prefix" >&2
    printf '%s\n' "$gromacs_version" | sed -n '/CUDA compiler:/p;/CUDA runtime:/p;/CUDA targets:/p' >&2
    exit 2
fi
if ! printf '%s\n' "$gromacs_version" | grep -Eq "CUDA compiler:.*${torch_cuda_re}"; then
    echo "GROMACS CUDA compiler does not match Torch CUDA ${torch_cuda}; set NMR_GROMACS_RUNTIME_PREFIX to the CUDA-matched GROMACS install prefix" >&2
    printf '%s\n' "$gromacs_version" | sed -n '/CUDA compiler:/p;/CUDA runtime:/p;/CUDA targets:/p' >&2
    exit 2
fi

echo "Torch: ${torch_version} (CUDA ${torch_cuda})"
echo "Vendor dir: $vendor_dir"
if [ "$copy" -eq 0 ]; then
    echo "Mode: dry-run (pass --copy to write payloads)"
else
    echo "Mode: copy"
fi

add_plan "$GROMACS_RUNTIME_PREFIX/" "$deps_dir/gromacs/" "GROMACS installed runtime prefix"
add_plan "$GROMACS_SOURCE_ROOT/" "$deps_dir/gromacs-src/" "GROMACS full source root"
add_plan "$GROMACS_BUILD/api/legacy/include/" "$deps_dir/gromacs-build/api/legacy/include/" "GROMACS generated legacy API headers"
add_plan "$REDUCE_SRC/" "$deps_dir/reduce-src/" "reduce source/build tree"
add_plan "$ORCA_ROOT/" "$deps_dir/orca/" "ORCA installed runtime prefix"

chem_env_root=$(dirname "$(dirname "$MOPAC_LIB")")
if [ ! -x "$chem_env_root/bin/mopac" ]; then
    echo "MISSING MOPAC executable in chem-env: $chem_env_root/bin/mopac" >&2
fi
if [ ! -x "$chem_env_root/bin/tleap" ]; then
    echo "MISSING AmberTools tleap in chem-env: $chem_env_root/bin/tleap" >&2
fi
add_plan "$chem_env_root/" "$deps_dir/chem-env/" "chem-env prefix (MOPAC + AmberTools)"

add_plan "$torch_root/" "$deps_dir/torch/" "Torch package root"
for name in cu13 cudnn cusparselt nccl nvshmem; do
    add_plan "$nvidia_root/$name/" "$deps_dir/nvidia/$name/" "NVIDIA Python package: $name"
done

printf '\nCopy plan:\n'
while IFS=$'\t' read -r src dst label; do
    copy_path "$src" "$dst" "$label"
done <"$plan_file"

mkdir -p "$vendor_dir"
manifest=$vendor_dir/manifest.toml
if [ "$copy" -eq 1 ]; then
    cat >"$manifest" <<EOF
# Generated by docker/collect-local-vendor-deps.sh
producer_torch = "${torch_version}"
producer_torch_cuda = "${torch_cuda}"
nvrtc_builtins = "${deps_dir}/nvidia/cu13/lib/libnvrtc-builtins.so.13.0"
cmake_cache = "${cache}"

[gromacs]
runtime_prefix = "${deps_dir}/gromacs"
gmx = "${deps_dir}/gromacs/bin/gmx_mpi"
topology_data = "${deps_dir}/gromacs/share/gromacs/top"
lib = "${deps_dir}/gromacs/lib/libgromacs_mpi.so"
muparser = "${deps_dir}/gromacs/lib/libmuparser.so.2"
source_root = "${deps_dir}/gromacs-src"
source_headers = "${deps_dir}/gromacs-src/src/gromacs"
legacy_api_headers = "${deps_dir}/gromacs-src/api/legacy/include/gromacs"
generated_headers = "${deps_dir}/gromacs-build/api/legacy/include/gromacs"

[orca]
prefix = "${deps_dir}/orca"
executable = "${deps_dir}/orca/orca"
lib_dir = "${deps_dir}/orca/lib"
EOF
    echo "WROTE $manifest"
else
    echo "DRY   manifest would be written to $manifest"
fi
