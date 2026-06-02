# Producer Docker Handoff - 2026-05-30

This is the handoff note for the single-container `nmr_extract` producer
appliance. It complements `docker/README.md` and `doc/INSTALLABILITY_PLAN.md`;
those remain the procedural and installability references.

## Scope

The Docker appliance is for `nmr_extract` only. It deliberately excludes
`h5-reader`, desktop UI dependencies, and reader/chewer development tooling.

The release shape is one image that can:

- start a local PostgreSQL instance in the container;
- restore the frozen tensorcs15 dump when needed;
- run `nmr-shielding-doctor --strict`;
- run `nmr-tensorcs15-check --quiet`;
- dispatch the requested `nmr_extract` job.

Compose may be useful for development, but the adviser-facing release contract
is the single-container appliance unless that requirement changes.

## Current Known-Good State

Built locally on 2026-05-30:

- dependency image: `nmr-shielding-producer-deps:local`
- producer image: `nmr-shielding:producer-full`
- PostgreSQL packages: `postgresql-18`, `postgresql-client-18`, `libpq5`,
  `libpq-dev`
- observed PostgreSQL version: `18.4-1.pgdg24.04+1`
- tensorcs15 dump embedded at:
  `/opt/nmr-shielding/vendor/tensorcs15/tensorcs15.dump`
- dump size observed in image: about 11 GB

Important: an older local image existed with Ubuntu default PostgreSQL 16. The
PG18 Dockerfile is intentional. Always verify the built image rather than
trusting the tag:

```bash
docker run --rm --entrypoint /bin/bash nmr-shielding:producer-full -lc \
  '/usr/lib/postgresql/18/bin/postgres --version && dpkg -l postgresql-18 postgresql-client-18 libpq5 libpq-dev'
```

## Build Recipe

Validate the vendored payload first:

```bash
docker/audit-vendor-payload.sh
```

Build the dependency base:

```bash
docker build \
  -f docker/deps-vendor.Dockerfile \
  -t nmr-shielding-producer-deps:local .
```

After the dependency build completes, prune completed BuildKit cache if disk
headroom matters. This keeps source and `docker/vendor/` intact:

```bash
docker builder prune -f
```

Build the producer image:

```bash
docker build \
  -f docker/producer-full.Dockerfile \
  --build-arg NMR_BASE_IMAGE=nmr-shielding-producer-deps:local \
  -t nmr-shielding:producer-full .
```

If disk is tight before the producer rebuild, remove only stale producer
runtime images that are about to be replaced:

```bash
docker image rm nmr-shielding:producer-full
```

Do not remove `nmr-shielding-producer-deps:local` before the producer rebuild,
and do not delete `docker/vendor/`.

## Disk Policy

This build uses large transient layers. The dependency image build can leave
about 100 GB of completed BuildKit cache; prune it right after the dependency
image is successfully tagged.

Good trims:

- `docker builder prune -f` after a completed large build;
- `docker image rm nmr-shielding:producer-full` before replacing a stale
  producer image;
- checking `docker system df` and `df -h /var/lib/docker` before and during
  large image export/unpack phases.

Bad trims:

- deleting source files;
- deleting `docker/vendor/`;
- deleting the dependency image before the producer image is rebuilt;
- pruning active BuildKit cache during a live build unless the filesystem is
  clearly about to run out of space.

During the 2026-05-30 PG18 dependency rebuild, free space dropped near 70 GB
during image unpack. Removing the stale PG16-era producer image recovered
headroom without touching source or vendor payloads.

## PostgreSQL Runtime Contract

`docker/nmr-container-entrypoint.sh` owns the in-container PostgreSQL lifecycle.

Defaults:

- `NMR_CONTAINER_ENABLE_POSTGRES=1`
- `PGDATA=/var/lib/postgresql/data`
- `PGDATABASE=tensorcs15`
- `PGUSER=nmr_shielding`
- `NMR_TENSORCS15_DSN="host=/var/run/postgresql dbname=tensorcs15 user=nmr_shielding"`
- `NMR_TENSORCS15_DUMP=/opt/nmr-shielding/vendor/tensorcs15/tensorcs15.dump`
- restore marker: `${PGDATA}/.tensorcs15-restored`

The entrypoint initializes PGDATA with trust auth for local container use,
listens only on the Unix socket directory, creates the role/database if needed,
restores the dump once, runs checks, then dispatches the job.

The entrypoint must remain a supervisor when it starts PostgreSQL. Do not
replace that path with a plain `exec`, or PostgreSQL will bypass the cleanup
trap. The current script waits for the child command and stops PostgreSQL on
normal exit or signal-driven termination.

For full restore testing, use a mounted PGDATA volume or an explicitly managed
scratch location. Avoid doing the 11 GB dump restore repeatedly into anonymous
Docker writable layers.

## Smoke Checks

Fast static image check:

```bash
docker run --rm --entrypoint /bin/bash nmr-shielding:producer-full -lc \
  'set -e; /usr/lib/postgresql/18/bin/postgres --version; test -f /opt/nmr-shielding/vendor/tensorcs15/tensorcs15.dump; grep -q terminate_child /usr/local/bin/nmr-container-entrypoint'
```

Bounded PostgreSQL startup check without restoring tensorcs15:

```bash
docker run --rm \
  --tmpfs /var/lib/postgresql/data:rw,size=1g \
  -e NMR_CONTAINER_RUN_DOCTOR=0 \
  -e NMR_CONTAINER_RUN_TENSORCS15_CHECK=0 \
  -e NMR_CONTAINER_CHECK_ORCA_SWAP=0 \
  -e NMR_TENSORCS15_DUMP=/missing/tensorcs15.dump \
  nmr-shielding:producer-full \
  psql -h /var/run/postgresql -U nmr_shielding -d tensorcs15 -tAc \
  "select current_database() || ' ' || current_user || ' ' || current_setting('server_version');"
```

Fast entrypoint/`nmr_extract` check without PostgreSQL:

```bash
docker run --rm \
  -e NMR_CONTAINER_ENABLE_POSTGRES=0 \
  -e NMR_CONTAINER_RUN_TENSORCS15_CHECK=0 \
  -e NMR_CONTAINER_CHECK_ORCA_SWAP=0 \
  nmr-shielding:producer-full --help >/dev/null
```

Full delivery acceptance is still:

```bash
NMR_DELIVERY_SWAP_MIN_GIB=64 docker/delivery-smoke.sh
```

That full check is intentionally heavier: it should initialize PostgreSQL,
restore tensorcs15, pass the manifest check, see the NVIDIA runtime, see ORCA
swap, and run GROMACS/ORCA.

## Future Docker Manager Notes

- Keep PostgreSQL pinned through PGDG package names (`postgresql-18`,
  `postgresql-client-18`). Do not silently fall back to Ubuntu's default
  `postgresql` metapackage.
- Keep DSNs out of logs. The export and check scripts should report that a DSN
  exists, not print credentials.
- Keep the vendored dependency payload deliberate. `docker/vendor/` is not
  build trash.
- Keep source-image context hygiene. The dependency Dockerfile consumes
  `docker/vendor/`; the producer Dockerfile consumes the dependency image and
  should not recopy vendored binary payloads into `/src`.
- Keep ORCA swap host-managed. The container should check/warn/fail based on
  visible swap policy, not create swap internally.
- Keep GROMACS CUDA, Torch CUDA, and user-space CUDA provenance aligned.
  Torch is currently `2.11.0+cu130`; host CUDA reporting a newer minor version
  does not change the image's CUDA 13.0 user-space ABI.
- Treat `nmr-shielding-doctor --strict` and `nmr-tensorcs15-check --quiet` as
  release gates, but keep expensive runtime probes opt-in.
- Avoid moving this into `h5-reader` or reader specs. The boundary is the
  producer appliance for `nmr_extract`.
