# syntax=docker/dockerfile:1.7

# This Dockerfile packages nmr_extract from a dependency-closed producer base.
# The dependency base may itself be built from curated vendored artifacts, but
# this source image must not accidentally copy host-local dependency paths.
ARG NMR_BASE_IMAGE=nmr-shielding-producer-deps:local
ARG NMR_PREFIX=/opt/nmr-shielding
ARG NMR_DEPS_PREFIX=/opt/nmr-shielding/deps

FROM ${NMR_BASE_IMAGE} AS build
ARG NMR_PREFIX
ARG NMR_DEPS_PREFIX
ARG NMR_CMAKE_GENERATOR=Ninja
ARG NMR_BUILD_TYPE=Release
ARG NMR_ENABLE_LTO=OFF

WORKDIR /src/nmr-shielding
COPY . .

RUN if [ -f /etc/ld.so.conf.d/nmr-shielding-deps.conf ]; then \
        sed -i "\#^${NMR_DEPS_PREFIX}/orca/lib\$#d" /etc/ld.so.conf.d/nmr-shielding-deps.conf; \
    fi \
 && ldconfig

RUN test -d "${NMR_DEPS_PREFIX}" || \
    (echo "Missing producer dependency prefix: ${NMR_DEPS_PREFIX}" >&2; exit 1)

RUN cmake -S . -B /tmp/nmr-shielding-build -G "${NMR_CMAKE_GENERATOR}" \
    -DCMAKE_BUILD_TYPE="${NMR_BUILD_TYPE}" \
    -DCMAKE_INSTALL_PREFIX="${NMR_PREFIX}" \
    -DCMAKE_INSTALL_SYSCONFDIR=/etc \
    -DNMR_PROFILE=producer-full \
    -DNMR_PORTABLE_DEP_ALLOW_REGEXES="^/usr/local/cuda/" \
    -DNMR_FETCH_GTEST=OFF \
    -DNMR_ENABLE_LTO="${NMR_ENABLE_LTO}" \
    -DNMR_SYSTEM_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu \
    -DNMR_GROMACS_ROOT="${NMR_DEPS_PREFIX}/gromacs" \
    -DGROMACS_LIB="${NMR_DEPS_PREFIX}/gromacs/lib/libgromacs_mpi.so" \
    -DGROMACS_SRC="${NMR_DEPS_PREFIX}/gromacs-src/src" \
    -DGROMACS_BUILD="${NMR_DEPS_PREFIX}/gromacs-build" \
    -DREDUCE_SRC="${NMR_DEPS_PREFIX}/reduce-src" \
    -DREDUCE_LIB="${NMR_DEPS_PREFIX}/reduce-src/build/reduce_src/libreducelib.a" \
    -DREDUCE_LIBPDB="${NMR_DEPS_PREFIX}/reduce-src/build/libpdb/libpdb++.a" \
    -DREDUCE_TOOLCLASSES="${NMR_DEPS_PREFIX}/reduce-src/build/toolclasses/libtoolclasses.a" \
    -DREDUCE_HET_DICT="${NMR_DEPS_PREFIX}/reduce-src/reduce_wwPDB_het_dict.txt" \
    -DNMR_MOPAC_ROOT="${NMR_DEPS_PREFIX}/chem-env" \
    -DNMR_TORCH_CMAKE_PREFIX_PATH="${NMR_DEPS_PREFIX}/torch/share/cmake" \
    -DNMR_NVRTC_LIB_DIR="${NMR_DEPS_PREFIX}/nvidia/cu13/lib"

RUN cmake --build /tmp/nmr-shielding-build --target nmr_extract -j"$(nproc)"
RUN cmake --install /tmp/nmr-shielding-build
RUN sed -i '1s@^#!.*python3$@#!/usr/bin/python3@' "${NMR_PREFIX}/bin/nmr-tensorcs15-check"

FROM ${NMR_BASE_IMAGE} AS runtime
ARG NMR_PREFIX
ARG NMR_DEPS_PREFIX

COPY --from=build ${NMR_PREFIX}/bin/ ${NMR_PREFIX}/bin/
COPY --from=build ${NMR_PREFIX}/share/ ${NMR_PREFIX}/share/
COPY --from=build /etc/nmr-shielding/ /etc/nmr-shielding/
RUN if [ -f /etc/ld.so.conf.d/nmr-shielding-deps.conf ]; then \
        sed -i "\#^${NMR_DEPS_PREFIX}/orca/lib\$#d" /etc/ld.so.conf.d/nmr-shielding-deps.conf; \
    fi \
 && ldconfig
COPY docker/nmr-container-entrypoint.sh /usr/local/bin/nmr-container-entrypoint
RUN chmod +x /usr/local/bin/nmr-container-entrypoint

ENV PATH="${NMR_PREFIX}/bin:${NMR_DEPS_PREFIX}/gromacs/bin:${NMR_DEPS_PREFIX}/orca:${NMR_DEPS_PREFIX}/chem-env/bin:${PATH}" \
    NMR_PREFIX="${NMR_PREFIX}" \
    NMR_DEPS_PREFIX="${NMR_DEPS_PREFIX}" \
    NMR_TOOLS_TOML=/etc/nmr-shielding/nmr_tools.toml \
    NMR_NVRTC_LIB_DIR="${NMR_DEPS_PREFIX}/nvidia/cu13/lib" \
    NMR_GMX="${NMR_DEPS_PREFIX}/gromacs/bin/gmx_mpi" \
    NMR_ORCA="${NMR_DEPS_PREFIX}/orca/orca" \
    NMR_SYSTEM_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial \
    NMR_DEPS_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial:${NMR_DEPS_PREFIX}/gromacs/lib:${NMR_DEPS_PREFIX}/torch/lib:${NMR_DEPS_PREFIX}/nvidia/cu13/lib:${NMR_DEPS_PREFIX}/nvidia/cudnn/lib:${NMR_DEPS_PREFIX}/nvidia/cusparselt/lib:${NMR_DEPS_PREFIX}/nvidia/nccl/lib:${NMR_DEPS_PREFIX}/nvidia/nvshmem/lib:${NMR_DEPS_PREFIX}/chem-env/lib"

RUN if [ ! -f /etc/nmr-shielding/nmr_tools.toml ]; then \
        cp /etc/nmr-shielding/nmr_tools.toml.example /etc/nmr-shielding/nmr_tools.toml; \
    fi

ENTRYPOINT ["tini", "--", "nmr-container-entrypoint"]
CMD ["--help"]
