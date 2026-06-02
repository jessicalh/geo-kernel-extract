#!/usr/bin/env bash
# Container entrypoint for the nmr_extract producer appliance.
# Starts the local tensorcs15 PostgreSQL service when requested, validates the
# install, then dispatches to nmr_extract via the CUDA library wrapper.

set -euo pipefail

NMR_PREFIX=${NMR_PREFIX:-/opt/nmr-shielding}
NMR_DEPS_PREFIX=${NMR_DEPS_PREFIX:-/opt/nmr-shielding/deps}
NMR_POSTGRES=${NMR_CONTAINER_ENABLE_POSTGRES:-1}
NMR_RUN_DOCTOR=${NMR_CONTAINER_RUN_DOCTOR:-1}
NMR_RUN_TENSORCS15_CHECK=${NMR_CONTAINER_RUN_TENSORCS15_CHECK:-1}
NMR_CHECK_ORCA_SWAP=${NMR_CONTAINER_CHECK_ORCA_SWAP:-1}
NMR_REQUIRE_SWAP=${NMR_CONTAINER_REQUIRE_SWAP:-0}
NMR_ORCA_SWAP_MIN_GIB=${NMR_CONTAINER_ORCA_SWAP_MIN_GIB:-0}
PGDATA=${PGDATA:-/var/lib/postgresql/data}
PGDATABASE=${PGDATABASE:-tensorcs15}
PGUSER=${PGUSER:-nmr_shielding}
NMR_TENSORCS15_DUMP=${NMR_TENSORCS15_DUMP:-/opt/nmr-shielding/vendor/tensorcs15/tensorcs15.dump}
NMR_TENSORCS15_RESTORE_MARKER=${NMR_TENSORCS15_RESTORE_MARKER:-${PGDATA}/.tensorcs15-restored}

export NMR_PREFIX NMR_DEPS_PREFIX
export NMR_SYSTEM_LIBRARY_PATH=${NMR_SYSTEM_LIBRARY_PATH:-/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial}
export NMR_DEPS_LIBRARY_PATH=${NMR_DEPS_LIBRARY_PATH:-${NMR_SYSTEM_LIBRARY_PATH}:${NMR_DEPS_PREFIX}/gromacs/lib:${NMR_DEPS_PREFIX}/torch/lib:${NMR_DEPS_PREFIX}/nvidia/cu13/lib:${NMR_DEPS_PREFIX}/nvidia/cudnn/lib:${NMR_DEPS_PREFIX}/nvidia/cusparselt/lib:${NMR_DEPS_PREFIX}/nvidia/nccl/lib:${NMR_DEPS_PREFIX}/nvidia/nvshmem/lib:${NMR_DEPS_PREFIX}/chem-env/lib}
export PATH="${NMR_DEPS_PREFIX}/gromacs/bin:${NMR_DEPS_PREFIX}/orca:${NMR_DEPS_PREFIX}/chem-env/bin:${PATH}"
export LD_LIBRARY_PATH="${NMR_DEPS_LIBRARY_PATH}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export NMR_TOOLS_TOML=${NMR_TOOLS_TOML:-/etc/nmr-shielding/nmr_tools.toml}
export NMR_NVRTC_LIB_DIR=${NMR_NVRTC_LIB_DIR:-${NMR_DEPS_PREFIX}/nvidia/cu13/lib}
export NMR_AIMNET2_MODEL=${NMR_AIMNET2_MODEL:-${NMR_PREFIX}/share/nmr-shielding/data/models/aimnet2_wb97m_0.jpt}
export NMR_MOPAC=${NMR_MOPAC:-${NMR_DEPS_PREFIX}/chem-env/bin/mopac}
export NMR_TLEAP=${NMR_TLEAP:-${NMR_DEPS_PREFIX}/chem-env/bin/tleap}
export NMR_GMX=${NMR_GMX:-${NMR_DEPS_PREFIX}/gromacs/bin/gmx_mpi}
export NMR_ORCA=${NMR_ORCA:-${NMR_DEPS_PREFIX}/orca/orca}
export NMR_FF14SB_PARAMS=${NMR_FF14SB_PARAMS:-${NMR_PREFIX}/share/nmr-shielding/data/ff14sb_params.dat}
export NMR_BMRB_ATOM_NOM=${NMR_BMRB_ATOM_NOM:-${NMR_PREFIX}/share/nmr-shielding/data/bmrb_atom_nom.tbl}
export NMR_TMPDIR=${NMR_TMPDIR:-/tmp/nmr_shielding}
export NMR_TENSORCS15_DSN=${NMR_TENSORCS15_DSN:-host=/var/run/postgresql dbname=${PGDATABASE} user=${PGUSER}}

postgres_started=0
child_pid=0

log() {
    printf 'nmr-container: %s\n' "$*" >&2
}

bool_enabled() {
    case "${1:-}" in
        1|true|TRUE|yes|YES|on|ON) return 0 ;;
        *) return 1 ;;
    esac
}

swap_total_kib() {
    if [ -r /proc/swaps ]; then
        awk 'NR > 1 { total += $3 } END { printf "%d\n", total + 0 }' /proc/swaps
    else
        printf '0\n'
    fi
}

ensure_orca_swap() {
    bool_enabled "$NMR_CHECK_ORCA_SWAP" || return 0
    [ -x "$NMR_ORCA" ] || return 0

    case "$NMR_ORCA_SWAP_MIN_GIB" in
        ""|*[!0-9]*)
            log "invalid NMR_CONTAINER_ORCA_SWAP_MIN_GIB: $NMR_ORCA_SWAP_MIN_GIB"
            exit 2
            ;;
    esac

    local total_kib min_kib
    total_kib=$(swap_total_kib)
    min_kib=$((NMR_ORCA_SWAP_MIN_GIB * 1024 * 1024))

    if [ "$total_kib" -le 0 ]; then
        if bool_enabled "$NMR_REQUIRE_SWAP"; then
            log "ORCA swap required but no swap is visible in /proc/swaps"
            exit 2
        fi
        log "warning: ORCA is available but no swap is visible in /proc/swaps"
        return 0
    fi

    if [ "$min_kib" -gt 0 ] && [ "$total_kib" -lt "$min_kib" ]; then
        if bool_enabled "$NMR_REQUIRE_SWAP"; then
            log "ORCA swap requirement not met: ${total_kib} KiB visible, need ${min_kib} KiB"
            exit 2
        fi
        log "warning: ORCA swap below requested floor: ${total_kib} KiB visible, requested ${min_kib} KiB"
    fi
}

postgres_bin_dir() {
    if command -v postgres >/dev/null 2>&1; then
        dirname "$(command -v postgres)"
        return 0
    fi
    local candidate
    for candidate in /usr/lib/postgresql/*/bin; do
        if [ -x "${candidate}/postgres" ]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done
    return 1
}

run_as_postgres() {
    if [ "$(id -u)" -eq 0 ]; then
        su -s /bin/sh postgres -c "$*"
    else
        sh -c "$*"
    fi
}

restore_tensorcs15_dump() {
    local pg_bin=$1
    local restore_log
    restore_log=$(mktemp /tmp/nmr_tensorcs15_restore.XXXXXX.log)

    set +e
    run_as_postgres "'${pg_bin}/pg_restore' --no-owner --role='${PGUSER}' -d '${PGDATABASE}' '${NMR_TENSORCS15_DUMP}'" \
        2> >(tee "$restore_log" >&2)
    local restore_status=$?
    set -e

    if [ "$restore_status" -eq 0 ]; then
        rm -f "$restore_log"
        return 0
    fi

    if grep -q 'setval: value 0 is out of bounds for sequence "raw_dft_calculations_calc_id_seq"' "$restore_log" \
        && grep -q 'pg_restore: warning: errors ignored on restore: 1' "$restore_log"; then
        log "warning: accepting tensorcs15 restore despite raw_dft_calculations_calc_id_seq setval(0) error"
        rm -f "$restore_log"
        return 0
    fi

    if grep -q 'input file appears to be a text format dump' "$restore_log"; then
        rm -f "$restore_log"
        run_as_postgres "'${pg_bin}/psql' -v ON_ERROR_STOP=1 -d '${PGDATABASE}' -f '${NMR_TENSORCS15_DUMP}'"
        return 0
    fi

    rm -f "$restore_log"
    return "$restore_status"
}

stop_postgres() {
    if [ "$postgres_started" -eq 1 ]; then
        local pg_bin
        pg_bin=$(postgres_bin_dir || true)
        if [ -n "$pg_bin" ]; then
            log "stopping PostgreSQL"
            run_as_postgres "'${pg_bin}/pg_ctl' -D '${PGDATA}' -m fast -w stop" || true
        fi
    fi
}

terminate_child() {
    local signal=$1
    if [ "$child_pid" -gt 0 ] && kill -0 "$child_pid" 2>/dev/null; then
        kill "-${child_pid}" 2>/dev/null || kill "$child_pid" 2>/dev/null || true
    fi
}

trap stop_postgres EXIT
trap 'terminate_child TERM' TERM
trap 'terminate_child INT' INT

ensure_runtime_config() {
    mkdir -p /etc/nmr-shielding "$NMR_TMPDIR"
    if [ ! -f "$NMR_TOOLS_TOML" ]; then
        if [ -f /etc/nmr-shielding/nmr_tools.toml.example ]; then
            cp /etc/nmr-shielding/nmr_tools.toml.example "$NMR_TOOLS_TOML"
        else
            cat >"$NMR_TOOLS_TOML" <<EOF
mopac = "${NMR_DEPS_PREFIX}/chem-env/bin/mopac"
tleap = "${NMR_DEPS_PREFIX}/chem-env/bin/tleap"
gmx = "${NMR_DEPS_PREFIX}/gromacs/bin/gmx_mpi"
orca = "${NMR_DEPS_PREFIX}/orca/orca"
ff14sb_params = "${NMR_FF14SB_PARAMS}"
tmpdir = "${NMR_TMPDIR}"
bmrb_atom_nom = "${NMR_BMRB_ATOM_NOM}"

[databases]
tensorcs15 = "${NMR_TENSORCS15_DSN}"
EOF
        fi
    fi
}

ensure_postgres() {
    bool_enabled "$NMR_POSTGRES" || return 0

    local pg_bin
    pg_bin=$(postgres_bin_dir) || {
        log "PostgreSQL requested but postgres binary was not found"
        exit 2
    }

    mkdir -p "$PGDATA" /var/run/postgresql
    if [ "$(id -u)" -eq 0 ]; then
        chown -R postgres:postgres "$PGDATA" /var/run/postgresql
    fi

    if [ ! -s "${PGDATA}/PG_VERSION" ]; then
        log "initializing PostgreSQL data directory at ${PGDATA}"
        run_as_postgres "'${pg_bin}/initdb' -D '${PGDATA}' --auth-local=trust --auth-host=trust"
        printf "listen_addresses = ''\nunix_socket_directories = '/var/run/postgresql'\n" >>"${PGDATA}/postgresql.conf"
        if [ "$(id -u)" -eq 0 ]; then
            chown postgres:postgres "${PGDATA}/postgresql.conf"
        fi
    fi

    log "starting PostgreSQL"
    run_as_postgres "'${pg_bin}/pg_ctl' -D '${PGDATA}' -w start"
    postgres_started=1

    run_as_postgres "'${pg_bin}/psql' -v ON_ERROR_STOP=1 -d postgres -tc \"SELECT 1 FROM pg_roles WHERE rolname='${PGUSER}'\" | grep -q 1 || '${pg_bin}/createuser' '${PGUSER}'"
    run_as_postgres "'${pg_bin}/psql' -v ON_ERROR_STOP=1 -d postgres -tc \"SELECT 1 FROM pg_database WHERE datname='${PGDATABASE}'\" | grep -q 1 || '${pg_bin}/createdb' -O '${PGUSER}' '${PGDATABASE}'"

    if [ -f "$NMR_TENSORCS15_DUMP" ] && [ ! -f "$NMR_TENSORCS15_RESTORE_MARKER" ]; then
        log "restoring tensorcs15 from ${NMR_TENSORCS15_DUMP}"
        restore_tensorcs15_dump "$pg_bin"
        touch "$NMR_TENSORCS15_RESTORE_MARKER"
        if [ "$(id -u)" -eq 0 ]; then
            chown postgres:postgres "$NMR_TENSORCS15_RESTORE_MARKER"
        fi
    elif [ ! -f "$NMR_TENSORCS15_DUMP" ]; then
        log "tensorcs15 dump not found at ${NMR_TENSORCS15_DUMP}; expecting mounted/restored database or optional checks disabled"
    fi
}

run_checks() {
    bool_enabled "$NMR_RUN_DOCTOR" && nmr-shielding-doctor --strict "${NMR_PREFIX}/bin/nmr_extract"
    if bool_enabled "$NMR_RUN_TENSORCS15_CHECK"; then
        nmr-tensorcs15-check --manifest "${NMR_PREFIX}/share/nmr-shielding/tensorcs15.manifest.toml" --quiet
    fi
}

normalize_args() {
    if [ "$#" -eq 0 ]; then
        set -- --help
    fi
    case "$1" in
        -* ) set -- nmr_extract "$@" ;;
    esac
    printf '%s\0' "$@"
}

ensure_runtime_config
ensure_orca_swap
ensure_postgres
run_checks

mapfile -d '' argv < <(normalize_args "$@")
command=("${argv[@]}")
if [ "${argv[0]}" = "nmr_extract" ]; then
    command=(nmr-run-with-cuda-env "${argv[@]}")
fi

if [ "$postgres_started" -eq 0 ]; then
    exec "${command[@]}"
fi

set +e
"${command[@]}" &
child_pid=$!
wait "$child_pid"
status=$?
child_pid=0
set -e
exit "$status"
