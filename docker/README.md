# Producer Docker

This directory contains the Docker packaging path for `nmr_extract` only. It
intentionally excludes `h5-reader` and all desktop UI dependencies.

The target is a self-contained producer appliance: instantiate the container,
let the entrypoint prepare the local runtime services, then run `nmr_extract`.

## Images

Build the dependency base from curated vendored payloads:

```bash
docker build \
  -f docker/deps-vendor.Dockerfile \
  -t nmr-shielding-producer-deps:local .
```

For disk-conservative rebuilds, trade time for space by pruning BuildKit cache before
large builds and removing stale runtime images before replacing them:

```bash
docker builder prune -af
docker image rm nmr-shielding:producer-full 2>/dev/null || true
```

Then build the producer image:

```bash
docker build \
  -f docker/producer-full.Dockerfile \
  --build-arg NMR_BASE_IMAGE=nmr-shielding-producer-deps:local \
  -t nmr-shielding:producer-full .
```

Prune BuildKit cache again after a successful producer rebuild if disk headroom
matters more than incremental rebuild speed:

```bash
docker builder prune -af
```

`producer-full.Dockerfile` configures with `NMR_PROFILE=producer-full`, installs
into `/opt/nmr-shielding`, and uses the same portability and doctor checks as
the non-container install contract.

## Vendored Payloads

Populate `docker/vendor/` before building the dependency base. See
`docker/vendor/README.md` for the required layout. The important point is that
all non-ordinary producer dependencies are deliberate payloads, not accidental
copies from a workstation path.

Start with a dry-run plan from the current local build cache:

```bash
# Set this to the CUDA-matched GROMACS install prefix.
NMR_GROMACS_RUNTIME_PREFIX=docker/vendor/deps/gromacs \
  docker/collect-local-vendor-deps.sh
```

After reviewing disk use and paths, populate the payload explicitly. The copy step dereferences source links; the vendored payload must not contain symlinks.


```bash
# Point this at the CUDA-matched GROMACS install prefix; the collector rejects a CUDA mismatch with Torch.
NMR_GROMACS_RUNTIME_PREFIX=docker/vendor/deps/gromacs \
  docker/collect-local-vendor-deps.sh --copy
```

Then validate it before building the base image:

```bash
docker/audit-vendor-payload.sh
```

Export the frozen tensorcs15 database into the vendor payload before the final
delivery build. The export script reads the DSN from `NMR_TENSORCS15_DSN` or
`~/.nmr_tools.toml` and does not print it:

```bash
docker/export-tensorcs15-dump.sh
```

Expected installed dependency prefix inside the image:

```text
/opt/nmr-shielding/deps/
  gromacs/                 # runnable GROMACS prefix: bin/gmx_mpi, share/gromacs/top, libs
  gromacs-src/             # full GROMACS source root retained for rebuilds
  gromacs-build/api/legacy/include/  # generated headers only; no compiler intermediates
  reduce-src/
  orca/                    # runnable ORCA prefix: orca plus lib/
  chem-env/               # single prefix with MOPAC + AmberTools/tleap
  torch/                  # libtorch 2.11.0+cu130
  nvidia/cu13/lib/        # CUDA 13.0 + libnvrtc-builtins.so.13.0
  nvidia/cudnn/lib/
  nvidia/cusparselt/lib/
  nvidia/nccl/lib/
  nvidia/nvshmem/lib/
```


## Source Data

The producer image installs runtime source data from `data/`, including
`data/models/aimnet2_wb97m_0.jpt` and `data/larsen_hbond_grids/`. The
producer-specific Docker ignore file keeps the large vendored dependency payload
out of the source build context while still including these runtime data files.

## CUDA/AIMNet2

AIMNet2 remains the primary CUDA consumer for nmr_extract. It runs through Torch/libtorch, and the current
producer Torch build is `2.11.0+cu130` with `torch.version.cuda == 13.0`. Pin
the container user-space CUDA/Torch closure to CUDA 13.0 (`cu130`) and run the
container with the NVIDIA container runtime so the host driver is injected. A
host reporting CUDA 13.1 is acceptable; it does not change the user-space ABI
that AIMNet/Torch should load.

GROMACS is a separate CUDA consumer. For final reproducibility, prefer a GROMACS
runtime built inside this Docker toolchain so its embedded device code, RUNPATH,
and CUDA provenance match the CUDA user-space carried by the image. The fallback
is a tool-specific GROMACS CUDA closure behind a wrapper; do not solve this by
putting a second CUDA or ORCA library tree on the global loader path.

## Runtime

The final image entrypoint is `nmr-container-entrypoint`. It:

- creates `/etc/nmr-shielding/nmr_tools.toml` if needed;
- puts the vendored `gmx_mpi` and `orca` executables on `PATH`;
- starts a local PostgreSQL instance when `NMR_CONTAINER_ENABLE_POSTGRES=1`;
- restores `/opt/nmr-shielding/vendor/tensorcs15/tensorcs15.dump` once when it
  is present;
- runs `nmr-shielding-doctor --strict` by default;
- runs `nmr-tensorcs15-check --quiet` by default;
- dispatches to `nmr_extract` through `nmr-run-with-cuda-env`.

Useful knobs:

```bash
NMR_CONTAINER_ENABLE_POSTGRES=0      # use an external tensorcs15 DSN
NMR_CONTAINER_RUN_DOCTOR=0           # skip startup doctor check
NMR_CONTAINER_RUN_TENSORCS15_CHECK=0 # skip DB manifest validation
NMR_CONTAINER_CHECK_ORCA_SWAP=1     # warn if ORCA is present but no swap is visible
NMR_CONTAINER_REQUIRE_SWAP=1        # fail startup if no swap is visible
NMR_CONTAINER_ORCA_SWAP_MIN_GIB=64  # optional visible-swap floor for ORCA runs
NMR_TMPDIR=/scratch/nmr_shielding   # mount this to fast local scratch for tool temp files
NMR_TENSORCS15_DSN=...               # override database connection
NMR_TENSORCS15_DUMP=/path/in/image   # override restore dump path
```

## ORCA Swap And Scratch

ORCA should use host-managed swap. Do not rely on the container creating or
enabling swap internally; `swapon` needs elevated container privileges and is a
host capacity decision. On 64 GB hosts, create and enable a host swapfile sized
for the ORCA workload, then let Docker expose that swap to the container. If a
container memory limit is used, set Docker swap policy deliberately, for example
`--memory=64g --memory-swap=-1` for unlimited host swap within the cgroup.

Use a bind-mounted scratch directory for ORCA-heavy temporary files and point
`NMR_TMPDIR` at it. The entrypoint and doctor read `/proc/swaps`; by default
they warn when ORCA is installed but no swap is visible. Set
`NMR_CONTAINER_REQUIRE_SWAP=1` and optionally
`NMR_CONTAINER_ORCA_SWAP_MIN_GIB=64` to make this a startup gate.

Host-side swapfile example:

```bash
sudo fallocate -l 128G /swapfile-orca
sudo chmod 600 /swapfile-orca
sudo mkswap /swapfile-orca
sudo swapon /swapfile-orca
```

Container run shape with explicit swap/scratch policy:

```bash
docker run --rm --gpus all \
  --memory=64g --memory-swap=-1 \
  -v /scratch/nmr_orca:/scratch/nmr_orca \
  -e NMR_TMPDIR=/scratch/nmr_orca \
  -e NMR_CONTAINER_REQUIRE_SWAP=1 \
  -e NMR_CONTAINER_ORCA_SWAP_MIN_GIB=64 \
  nmr-shielding:producer-full \
  <nmr_extract mode arguments>
```


Example run shape:

```bash
docker run --rm --gpus all \
  -v "$PWD/input:/input:ro" \
  -v "$PWD/output:/output" \
  nmr-shielding:producer-full \
  <nmr_extract mode arguments>
```

## Delivery Acceptance

The delivery image is not considered self-contained until default startup can
initialize PostgreSQL, restore `tensorcs15.dump`, pass the tensorcs15 manifest
check, see the NVIDIA runtime, see ORCA swap, and run GROMACS/ORCA from the
container PATH.

Run this on the build machine, the DGX Spark, and the advisor machine:

```bash
NMR_DELIVERY_SWAP_MIN_GIB=64 docker/delivery-smoke.sh
```

For a real extraction smoke, mount the repository fixture or a project PDB and
enable the optional PDB path:

```bash
NMR_DELIVERY_RUN_PDB_SMOKE=1 \
NMR_DELIVERY_PDB_INPUT=tests/data/1ubq_protonated.pdb \
NMR_DELIVERY_SWAP_MIN_GIB=64 \
docker/delivery-smoke.sh
```
