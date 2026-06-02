# Vendored Producer Dependencies

This directory is intentionally small in git. Before building
`docker/deps-vendor.Dockerfile`, populate it with the curated dependency payloads
for the `nmr_extract` producer appliance.

Required layout:

```text
docker/vendor/
  deps/
    gromacs/                 installed GROMACS prefix with bin/gmx_mpi, share/gromacs/top, libs
    gromacs-src/             matching full GROMACS source root, kept for rebuilds
    gromacs-build/api/legacy/include/ matching generated headers only
    reduce-src/              reduce tree with built static archives and het dict
    orca/                    ORCA prefix with orca executable and lib/
    chem-env/                single chemistry env prefix with bin/mopac and bin/tleap
    torch/                   libtorch 2.11.0+cu130 prefix with TorchConfig.cmake
    nvidia/cu13/lib/         CUDA 13.0 user-space libraries and nvrtc builtins
    nvidia/cudnn/lib/        cuDNN library family required by Torch cu130
    nvidia/cusparselt/lib/   cuSPARSELt library family required by Torch cu130
    nvidia/nccl/lib/         NCCL library family required by Torch cu130
    nvidia/nvshmem/lib/      NVSHMEM library family required by Torch cu130
  tensorcs15/
    tensorcs15.dump          PostgreSQL dump restored at container startup
```

The dependency base image copies `deps/` to `/opt/nmr-shielding/deps` and
`tensorcs15/` to `/opt/nmr-shielding/vendor/tensorcs15`. The default container
startup path restores and validates this dump, so it is required for a
self-contained delivery image. Keep the payload curated and versioned outside
the source repository until we decide how to store large binary artifacts.

Keep source trees that are needed to rebuild the image on another CUDA host, but
do not vendor compiler intermediates. The only retained GROMACS build outputs
are the generated legacy API headers under `gromacs-build/api/legacy/include/`.

GROMACS and ORCA are required runtime tools for the producer image. They are not
optional extras: the container startup doctor checks that `gmx_mpi` and `orca`
resolve to executable files inside the vendored dependency prefix.

Use `docker/collect-local-vendor-deps.sh` for a dry-run copy plan from the
current local CMake cache. The collector probes `gmx_mpi --version` and rejects a GROMACS CUDA runtime/compiler that does not match the Torch CUDA line; set `NMR_GROMACS_RUNTIME_PREFIX` to the CUDA-matched GROMACS install prefix when collecting.
Use `--copy` only after reviewing the plan and disk
space. Use `docker/audit-vendor-payload.sh` to validate the populated payload.
The vendored payload must not contain symlinks; the collector dereferences
source links during copy and the auditor rejects any remaining links.

## AIMNet2 CUDA Contract

AIMNet2 consumes CUDA through Torch/libtorch. The current producer Torch is
`2.11.0+cu130`, reports `torch.version.cuda == 13.0`, and requires
`libnvrtc-builtins.so.13.0`. Vendor the CUDA 13.0/Torch cu130 user-space
closure even when the host driver reports a newer maximum CUDA version. The
NVIDIA container runtime supplies the host driver at run time.
