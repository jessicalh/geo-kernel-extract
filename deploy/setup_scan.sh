#!/bin/bash
#
# setup_scan.sh — Build nmr_extract on a scan machine (Ubuntu 22.04).
#
# Prerequisites already present:
#   - Eigen3 (apt)
#   - GROMACS 2026 library + source/build trees
#   - libmopac + headers
#   - cmake, g++, make, git
#
# Ubuntu 22.04 differences from 24.04:
#   - No libdssp-dev (must build dssp+cifpp from source)
#   - No libfetk-dev (use libmaloc-dev for APBS link)
#
# Usage:
#   ssh scan1
#   bash ~/nmr-shielding/deploy/setup_scan.sh
#
# After running, nmr_extract is at ~/nmr-shielding/build/nmr_extract
#

set -euo pipefail

NMR_PROJECT_ROOT=${NMR_PROJECT_ROOT:-$HOME/nmr-shielding}
NMR_BUILD_PARENT=${NMR_BUILD_PARENT:-$HOME/builds}
NMR_BUILD_DIR=${NMR_BUILD_DIR:-$NMR_PROJECT_ROOT/build}
NMR_DSSP_SRC=${NMR_DSSP_SRC:-$NMR_BUILD_PARENT/dssp-src}
NMR_REDUCE_SRC=${NMR_REDUCE_SRC:-$NMR_BUILD_PARENT/reduce-src}
NMR_DSSP_REPO=${NMR_DSSP_REPO:-https://github.com/PDB-REDO/dssp.git}
NMR_REDUCE_REPO=${NMR_REDUCE_REPO:-https://github.com/rlabduke/reduce.git}
NMR_SCAN_INSTALL_APT=${NMR_SCAN_INSTALL_APT:-1}
NMR_SCAN_DISABLE_KITWARE_REPO=${NMR_SCAN_DISABLE_KITWARE_REPO:-$NMR_SCAN_INSTALL_APT}
# Source builds remain the scan-machine default. Set these to 0 only on
# Docker/offline/pre-provisioned machines where CMake should use existing deps.
NMR_SCAN_BUILD_DSSP=${NMR_SCAN_BUILD_DSSP:-1}
NMR_SCAN_BUILD_REDUCE=${NMR_SCAN_BUILD_REDUCE:-1}
NMR_SCAN_INSTALL_PREFIX=${NMR_SCAN_INSTALL_PREFIX:-/usr/local}
NMR_SCAN_APT_PACKAGES=${NMR_SCAN_APT_PACKAGES:-"libopenbabel-dev libcifpp-dev libapbs-dev libmaloc-dev libboost-dev"}
NMR_GROMACS_ROOT=${NMR_GROMACS_ROOT:-/opt/gromacs/2026.0}
NMR_GROMACS_LIB=${NMR_GROMACS_LIB:-$NMR_GROMACS_ROOT/lib/libgromacs_mpi.so}
NMR_GROMACS_SRC=${NMR_GROMACS_SRC:-$HOME/build/gmx2026/gromacs-2026.0/src}
NMR_GROMACS_BUILD=${NMR_GROMACS_BUILD:-$HOME/build/gmx2026/gromacs-2026-build}
NMR_DSSP_ROOT=${NMR_DSSP_ROOT:-$NMR_SCAN_INSTALL_PREFIX}
NMR_MOPAC_ROOT=${NMR_MOPAC_ROOT:-${CONDA_PREFIX:-}}
NMR_MOPAC_INCLUDE=${NMR_MOPAC_INCLUDE:-${NMR_MOPAC_ROOT:+$NMR_MOPAC_ROOT/include}}
NMR_MOPAC_LIB=${NMR_MOPAC_LIB:-${NMR_MOPAC_ROOT:+$NMR_MOPAC_ROOT/lib/libmopac.so}}
NMR_OPENBABEL_LIB=${NMR_OPENBABEL_LIB:-}
NMR_SCAN_EXTRA_TARGETS=${NMR_SCAN_EXTRA_TARGETS:-}

# Disable broken kitware repo if present
if [ "$NMR_SCAN_DISABLE_KITWARE_REPO" = "1" ] && [ -f /etc/apt/sources.list.d/kitware.list ]; then
    sudo mv /etc/apt/sources.list.d/kitware.list /etc/apt/sources.list.d/kitware.list.disabled
fi

if [ "$NMR_SCAN_INSTALL_APT" = "1" ]; then
    echo "=== Installing apt packages ==="
    export DEBIAN_FRONTEND=noninteractive
    sudo apt-get update -qq
    sudo apt-get install -y -qq $NMR_SCAN_APT_PACKAGES
else
    echo "=== Skipping apt packages (NMR_SCAN_INSTALL_APT=$NMR_SCAN_INSTALL_APT) ==="
fi

# ============================================================================
# Build dssp from source (no libdssp-dev on 22.04)
# ============================================================================
echo "=== Building dssp from source ==="
mkdir -p "$NMR_BUILD_PARENT"
if [ "$NMR_SCAN_BUILD_DSSP" != "1" ]; then
    echo "dssp source build skipped (NMR_SCAN_BUILD_DSSP=$NMR_SCAN_BUILD_DSSP)"
elif [ ! -f "$NMR_SCAN_INSTALL_PREFIX/lib/libdssp.a" ] && [ ! -f "$NMR_SCAN_INSTALL_PREFIX/lib/libdssp.so" ]; then
    cd "$NMR_BUILD_PARENT"
    if [ ! -d "$NMR_DSSP_SRC" ]; then
        git clone "$NMR_DSSP_REPO" "$NMR_DSSP_SRC"
    fi
    cd "$NMR_DSSP_SRC"
    mkdir -p build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$NMR_SCAN_INSTALL_PREFIX" 2>&1 | tail -5
    make -j$(nproc) 2>&1 | tail -5
    sudo make install 2>&1 | tail -3
    echo "dssp installed"
else
    echo "dssp already installed"
fi

# ============================================================================
# Build reduce from source
# ============================================================================
echo "=== Building reduce from source ==="
if [ "$NMR_SCAN_BUILD_REDUCE" != "1" ]; then
    echo "reduce source build skipped (NMR_SCAN_BUILD_REDUCE=$NMR_SCAN_BUILD_REDUCE)"
elif [ ! -f "$NMR_REDUCE_SRC/build/reduce_src/libreducelib.a" ]; then
    cd "$NMR_BUILD_PARENT"
    if [ ! -d "$NMR_REDUCE_SRC" ]; then
        git clone "$NMR_REDUCE_REPO" "$NMR_REDUCE_SRC"
    fi
    cd "$NMR_REDUCE_SRC"
    mkdir -p build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
    make -j$(nproc) reducelib 2>&1 | tail -5
    echo "reduce built: $(ls -la reduce_src/libreducelib.a)"
else
    echo "reduce already built"
fi

# ============================================================================
# Configure and build nmr_extract
# ============================================================================
echo "=== Configuring CMake ==="
cd "$NMR_PROJECT_ROOT"
mkdir -p "$NMR_BUILD_DIR"

REDUCE_LIB="$NMR_REDUCE_SRC/build/reduce_src/libreducelib.a"
REDUCE_LIBPDB="$NMR_REDUCE_SRC/build/libpdb/libpdb++.a"
REDUCE_TOOLCLASSES="$NMR_REDUCE_SRC/build/toolclasses/libtoolclasses.a"

# OpenBabel on Ubuntu 22.04
if [ -n "$NMR_OPENBABEL_LIB" ]; then
    OB_LIB=$NMR_OPENBABEL_LIB
elif [ -f /usr/lib/x86_64-linux-gnu/libopenbabel.so ]; then
    OB_LIB="/usr/lib/x86_64-linux-gnu/libopenbabel.so"
elif [ -f /usr/lib/x86_64-linux-gnu/libopenbabel.so.7 ]; then
    OB_LIB="/usr/lib/x86_64-linux-gnu/libopenbabel.so.7"
else
    echo "ERROR: libopenbabel not found"
    exit 1
fi

cmake_args=(
    -DCMAKE_BUILD_TYPE=Release
    -DGROMACS_LIB="$NMR_GROMACS_LIB"
    -DGROMACS_SRC="$NMR_GROMACS_SRC"
    -DGROMACS_BUILD="$NMR_GROMACS_BUILD"
    -DNMR_DSSP_ROOT="$NMR_DSSP_ROOT"
    -DOPENBABEL_LIB="$OB_LIB"
    -DREDUCE_SRC="$NMR_REDUCE_SRC"
    -DREDUCE_LIB="$REDUCE_LIB"
    -DREDUCE_LIBPDB="$REDUCE_LIBPDB"
    -DREDUCE_TOOLCLASSES="$REDUCE_TOOLCLASSES"
)
[ -n "${NMR_MOPAC_INCLUDE:-}" ] && cmake_args+=(-DMOPAC_INCLUDE="$NMR_MOPAC_INCLUDE")
[ -n "${NMR_MOPAC_LIB:-}" ] && cmake_args+=(-DMOPAC_LIB="$NMR_MOPAC_LIB")

cmake -S "$NMR_PROJECT_ROOT" -B "$NMR_BUILD_DIR" "${cmake_args[@]}" \
    2>&1 | tail -20

build_targets=(nmr_extract)
for target in $NMR_SCAN_EXTRA_TARGETS; do
    build_targets+=("$target")
done

echo "=== Building ${build_targets[*]} ==="
cmake --build "$NMR_BUILD_DIR" --target "${build_targets[@]}" -j"$(nproc)" 2>&1 | tail -10

echo ""
echo "=== Done ==="
echo "Binary: $NMR_BUILD_DIR/nmr_extract"
echo ""
echo "Test with:"
echo "  NMR_NVRTC_LIB_DIR=/path/to/nvidia/cu13/lib \\"
echo "    scripts/run_with_cuda_env.sh $NMR_BUILD_DIR/<test-target>"
