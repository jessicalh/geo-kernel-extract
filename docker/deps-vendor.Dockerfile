# syntax=docker/dockerfile:1.7

# Build the dependency-closed base image consumed by producer-full.Dockerfile.
# The apt layer supplies ordinary OS/runtime/build packages. The non-ordinary
# producer dependencies are expected under docker/vendor/deps and are copied as
# curated vendored artifacts.
ARG NMR_DEPS_BASE_IMAGE=nvidia/cuda:13.0.2-cudnn-devel-ubuntu24.04
ARG NMR_DEPS_PREFIX=/opt/nmr-shielding/deps

FROM ${NMR_DEPS_BASE_IMAGE}
ARG NMR_DEPS_PREFIX
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    gnupg \
 && install -d /usr/share/postgresql-common/pgdg \
 && curl -fsSL https://www.postgresql.org/media/keys/ACCC4CF8.asc -o /usr/share/postgresql-common/pgdg/apt.postgresql.org.asc \
 && printf '%s\n' 'deb [signed-by=/usr/share/postgresql-common/pgdg/apt.postgresql.org.asc] https://apt.postgresql.org/pub/repos/apt noble-pgdg main' > /etc/apt/sources.list.d/pgdg.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    cmake \
    coreutils \
    findutils \
    g++ \
    gcc \
    gfortran \
    libcifpp-dev \
    libapbs-dev \
    libfetk-dev \
    libblas-dev \
    libboost-dev \
    libcurl4-openssl-dev \
    dssp \
    libeigen3-dev \
    libgtest-dev \
    libhdf5-dev \
    liblapack-dev \
    libmaloc-dev \
    libopenbabel-dev \
    libopenmpi-dev \
    libpq-dev \
    ninja-build \
    pkg-config \
    postgresql-18 \
    postgresql-client-18 \
    procps \
    python3 \
    python3-numpy \
    python3-psycopg2 \
    tini \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

RUN if [ ! -f /usr/include/apbs/apbscfg.h ]; then \
        printf '%s\n' \
            '/* Generated for APBS 3.4.1 on Ubuntu 24.04 container builds. */' \
            '#ifndef APBSCFG_H' \
            '#define APBSCFG_H' \
            '#define PACKAGE_STRING "APBS 3.4.1"' \
            '#define APBS_VERSION "3.4.1"' \
            '#define HAVE_MG_SOLVER 1' \
            '#define HAVE_UNISTD_H 1' \
            '#define HAVE_STDLIB_H 1' \
            '#define HAVE_STRING_H 1' \
            '#endif /* APBSCFG_H */' \
            > /usr/include/apbs/apbscfg.h; \
    fi

COPY docker/vendor/deps/ ${NMR_DEPS_PREFIX}/
COPY docker/vendor/tensorcs15/ /opt/nmr-shielding/vendor/tensorcs15/

RUN first_link="$(find "${NMR_DEPS_PREFIX}" /opt/nmr-shielding/vendor/tensorcs15 -type l -print -quit)" \
 && test -z "$first_link" || (echo "Vendored payload contains a symlink/link: $first_link" >&2; exit 1)

RUN test -x "${NMR_DEPS_PREFIX}/gromacs/bin/gmx_mpi" \
 && test -d "${NMR_DEPS_PREFIX}/gromacs/share/gromacs/top" \
 && test -f "${NMR_DEPS_PREFIX}/gromacs/lib/libgromacs_mpi.so" \
 && test -f "${NMR_DEPS_PREFIX}/gromacs/lib/libmuparser.so.2" \
 && test -f "${NMR_DEPS_PREFIX}/gromacs-src/CMakeLists.txt" \
 && test -d "${NMR_DEPS_PREFIX}/gromacs-src/src/gromacs" \
 && test -f "${NMR_DEPS_PREFIX}/gromacs-src/api/legacy/include/gromacs/utility/basedefinitions.h" \
 && test -f "${NMR_DEPS_PREFIX}/gromacs-build/api/legacy/include/gromacs/libgromacs_export.h" && test -f "${NMR_DEPS_PREFIX}/gromacs-build/api/legacy/include/gromacs/version.h" && test -z "$(find "${NMR_DEPS_PREFIX}/gromacs-build" \( -name CMakeCache.txt -o -name CMakeFiles -o -name Makefile -o -name '*.o' -o -name '*.a' \) -print -quit)" \
 && test -x "${NMR_DEPS_PREFIX}/orca/orca" \
 && test -f "${NMR_DEPS_PREFIX}/reduce-src/build/reduce_src/libreducelib.a" \
 && test -f "${NMR_DEPS_PREFIX}/reduce-src/build/libpdb/libpdb++.a" \
 && test -f "${NMR_DEPS_PREFIX}/reduce-src/build/toolclasses/libtoolclasses.a" \
 && test -f "${NMR_DEPS_PREFIX}/reduce-src/reduce_wwPDB_het_dict.txt" \
 && test -f "${NMR_DEPS_PREFIX}/chem-env/include/mopac.h" \
 && test -f "${NMR_DEPS_PREFIX}/chem-env/lib/libmopac.so" \
 && test -x "${NMR_DEPS_PREFIX}/chem-env/bin/mopac" \
 && test -x "${NMR_DEPS_PREFIX}/chem-env/bin/tleap" \
 && test -f "${NMR_DEPS_PREFIX}/torch/share/cmake/Torch/TorchConfig.cmake" \
 && test -f "${NMR_DEPS_PREFIX}/nvidia/cu13/lib/libnvrtc-builtins.so.13.0" \
 && test -d "${NMR_DEPS_PREFIX}/nvidia/cudnn/lib" \
 && test -d "${NMR_DEPS_PREFIX}/nvidia/cusparselt/lib" \
 && test -d "${NMR_DEPS_PREFIX}/nvidia/nccl/lib" \
 && test -d "${NMR_DEPS_PREFIX}/nvidia/nvshmem/lib"

RUN printf '%s\n' \
    "${NMR_DEPS_PREFIX}/gromacs/lib" \
    "${NMR_DEPS_PREFIX}/torch/lib" \
    "${NMR_DEPS_PREFIX}/nvidia/cu13/lib" \
    "${NMR_DEPS_PREFIX}/nvidia/cudnn/lib" \
    "${NMR_DEPS_PREFIX}/nvidia/cusparselt/lib" \
    "${NMR_DEPS_PREFIX}/nvidia/nccl/lib" \
    "${NMR_DEPS_PREFIX}/nvidia/nvshmem/lib" \
    > /etc/ld.so.conf.d/nmr-shielding-deps.conf \
 && ldconfig

ENV NMR_DEPS_PREFIX=${NMR_DEPS_PREFIX} \
    NMR_NVRTC_LIB_DIR=${NMR_DEPS_PREFIX}/nvidia/cu13/lib \
    NMR_SYSTEM_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial \
    NMR_DEPS_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial:${NMR_DEPS_PREFIX}/gromacs/lib:${NMR_DEPS_PREFIX}/torch/lib:${NMR_DEPS_PREFIX}/nvidia/cu13/lib:${NMR_DEPS_PREFIX}/nvidia/cudnn/lib:${NMR_DEPS_PREFIX}/nvidia/cusparselt/lib:${NMR_DEPS_PREFIX}/nvidia/nccl/lib:${NMR_DEPS_PREFIX}/nvidia/nvshmem/lib:${NMR_DEPS_PREFIX}/chem-env/lib
