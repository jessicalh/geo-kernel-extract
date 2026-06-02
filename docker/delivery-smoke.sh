#!/usr/bin/env bash
# Acceptance smoke test for the Docker producer delivery image.

set -euo pipefail

IMAGE=${NMR_DELIVERY_IMAGE:-nmr-shielding:producer-full}
SWAP_MIN_GIB=${NMR_DELIVERY_SWAP_MIN_GIB:-64}
GPU_ARGS=${NMR_DOCKER_GPU_ARGS:---gpus all}
OUT_ROOT=${NMR_DELIVERY_OUT_ROOT:-"$(pwd)/artifacts/docker-delivery-smoke"}
RUN_PDB_SMOKE=${NMR_DELIVERY_RUN_PDB_SMOKE:-0}
PDB_SMOKE_INPUT=${NMR_DELIVERY_PDB_INPUT:-"$(pwd)/tests/data/1ubq_protonated.pdb"}

log() {
    printf 'delivery-smoke: %s\n' "$*" >&2
}

fail() {
    log "$*"
    exit 2
}

bool_enabled() {
    case "${1:-}" in
        1|true|TRUE|yes|YES|on|ON) return 0 ;;
        *) return 1 ;;
    esac
}

docker_run() {
    # shellcheck disable=SC2086
    docker run --rm $GPU_ARGS "$@"
}

docker image inspect "$IMAGE" >/dev/null || fail "image not found: $IMAGE"

case "$SWAP_MIN_GIB" in
    ""|*[!0-9]*) fail "invalid NMR_DELIVERY_SWAP_MIN_GIB: $SWAP_MIN_GIB" ;;
esac

mkdir -p "$OUT_ROOT"

log "1/6 host GPU visibility through NVIDIA container runtime"
docker_run \
    -e NMR_CONTAINER_ENABLE_POSTGRES=0 \
    -e NMR_CONTAINER_RUN_DOCTOR=0 \
    -e NMR_CONTAINER_RUN_TENSORCS15_CHECK=0 \
    "$IMAGE" nvidia-smi

log "2/6 fast strict doctor plus nmr_extract help"
docker_run \
    -e NMR_CONTAINER_ENABLE_POSTGRES=0 \
    -e NMR_CONTAINER_RUN_TENSORCS15_CHECK=0 \
    -e NMR_CONTAINER_REQUIRE_SWAP=1 \
    -e NMR_CONTAINER_ORCA_SWAP_MIN_GIB="$SWAP_MIN_GIB" \
    "$IMAGE" --help >/dev/null

log "3/6 GROMACS is runnable"
docker_run \
    -e NMR_CONTAINER_ENABLE_POSTGRES=0 \
    -e NMR_CONTAINER_RUN_DOCTOR=0 \
    -e NMR_CONTAINER_RUN_TENSORCS15_CHECK=0 \
    -e NMR_CONTAINER_REQUIRE_SWAP=1 \
    -e NMR_CONTAINER_ORCA_SWAP_MIN_GIB="$SWAP_MIN_GIB" \
    "$IMAGE" gmx_mpi --version >/dev/null

log "4/6 ORCA is runnable from PATH"
set +e
orca_output=$(docker_run \
    -e NMR_CONTAINER_ENABLE_POSTGRES=0 \
    -e NMR_CONTAINER_RUN_DOCTOR=0 \
    -e NMR_CONTAINER_RUN_TENSORCS15_CHECK=0 \
    -e NMR_CONTAINER_REQUIRE_SWAP=1 \
    -e NMR_CONTAINER_ORCA_SWAP_MIN_GIB="$SWAP_MIN_GIB" \
    "$IMAGE" orca 2>&1)
orca_status=$?
set -e
if [ "$orca_status" -ne 0 ] && ! printf '%s\n' "$orca_output" | grep -qi 'input'; then
    printf '%s\n' "$orca_output" >&2
    fail "ORCA smoke failed unexpectedly"
fi

log "5/6 default self-contained startup, tensorcs15 restore, and manifest check"
docker_run \
    -e NMR_CONTAINER_REQUIRE_SWAP=1 \
    -e NMR_CONTAINER_ORCA_SWAP_MIN_GIB="$SWAP_MIN_GIB" \
    "$IMAGE" --help >/dev/null

log "6/6 optional PDB extraction smoke"
if bool_enabled "$RUN_PDB_SMOKE"; then
    [ -f "$PDB_SMOKE_INPUT" ] || fail "PDB smoke input not found: $PDB_SMOKE_INPUT"
    smoke_dir="${OUT_ROOT}/$(date +%Y%m%d-%H%M%S)"
    mkdir -p "$smoke_dir"
    docker_run \
        -v "$PDB_SMOKE_INPUT:/input/input.pdb:ro" \
        -v "$smoke_dir:/output" \
        -e NMR_CONTAINER_REQUIRE_SWAP=1 \
        -e NMR_CONTAINER_ORCA_SWAP_MIN_GIB="$SWAP_MIN_GIB" \
        "$IMAGE" \
        --protonated-pdb /input/input.pdb \
        --no-mopac \
        --aimnet2 /opt/nmr-shielding/share/nmr-shielding/data/models/aimnet2_wb97m_0.jpt \
        --output /output
    find "$smoke_dir" -maxdepth 1 -name '*.npy' -type f -print -quit | grep -q . \
        || fail "PDB smoke produced no NPY files in $smoke_dir"
else
    log "skipped; set NMR_DELIVERY_RUN_PDB_SMOKE=1 to run it"
fi

log "PASS $IMAGE"
