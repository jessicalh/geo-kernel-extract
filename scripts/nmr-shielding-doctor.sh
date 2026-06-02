#!/usr/bin/env bash
#
# Lightweight producer installability preflight.
# Default mode is non-destructive: no builds, no tests, no MOPAC/AIMNet runs,
# no database connection attempts.

set -u

deep=0
strict=${NMR_DOCTOR_STRICT:-0}
extract_arg=""

usage() {
    cat <<'USAGE'
usage: nmr-shielding-doctor [--deep] [--strict] [path/to/nmr_extract]

Default checks inspect configured paths, dynamic links, and obvious runtime
closure issues. --strict promotes machine-local dependency leaks to failures.
--deep enables opt-in runtime probes that may connect to local services.
USAGE
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --deep) deep=1 ;;
        --strict) strict=1 ;;
        -h|--help) usage; exit 0 ;;
        *) extract_arg=$1 ;;
    esac
    shift
done
case "$strict" in
    1|true|TRUE|yes|YES|on|ON) strict=1 ;;
    *) strict=0 ;;
esac

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(CDPATH= cd -- "${SCRIPT_DIR}/.." && pwd)
LOCAL_PATH_RE='^/(home|shared|mnt)/|^/opt/orca/|^/usr/local/cuda/'
LOCAL_RUNPATH_RE='/home/|/shared/|/mnt/|/opt/orca/|/usr/local/cuda/'

failures=0
warnings=0

info() { printf 'INFO  %s\n' "$*"; }
ok() { printf 'OK    %s\n' "$*"; }
warn() { warnings=$((warnings + 1)); printf 'WARN  %s\n' "$*" >&2; }
fail() { failures=$((failures + 1)); printf 'FAIL  %s\n' "$*" >&2; }
leak() {
    if [ "$strict" = "1" ]; then
        fail "$@"
    else
        warn "$@"
    fi
}

have_cmd() { command -v "$1" >/dev/null 2>&1; }
is_enabled() {
    case "${1:-}" in
        1|true|TRUE|yes|YES|on|ON) return 0 ;;
        *) return 1 ;;
    esac
}


check_file() {
    local label=$1 path=$2
    if [ -n "$path" ] && [ -f "$path" ]; then
        ok "$label: $path"
    else
        fail "$label missing${path:+: $path}"
    fi
}

check_exec_file() {
    local label=$1 path=$2
    if [ -n "$path" ] && [ -x "$path" ]; then
        ok "$label: $path"
    else
        fail "$label missing/not executable${path:+: $path}"
    fi
}

check_optional_file() {
    local label=$1 path=$2
    if [ -z "$path" ]; then
        warn "$label not configured"
    elif [ -f "$path" ]; then
        ok "$label: $path"
    else
        fail "$label configured but missing: $path"
    fi
}

check_optional_dir() {
    local label=$1 path=$2
    if [ -z "$path" ]; then
        warn "$label not configured"
    elif [ -d "$path" ]; then
        ok "$label: $path"
    else
        fail "$label configured but missing: $path"
    fi
}

check_temp_dir() {
    local label=$1 path=$2 parent
    if [ -z "$path" ]; then
        warn "$label not configured"
    elif [ -d "$path" ]; then
        if [ -w "$path" ]; then
            ok "$label: $path"
        else
            fail "$label is not writable: $path"
        fi
    else
        parent=$(dirname -- "$path")
        if [ -d "$parent" ] && [ -w "$parent" ]; then
            ok "$label can be created by runtime: $path"
        else
            fail "$label cannot be created: $path"
        fi
    fi
}

swap_total_kib() {
    if [ -r /proc/swaps ]; then
        awk 'NR > 1 { total += $3 } END { printf "%d\n", total + 0 }' /proc/swaps
    else
        printf '0\n'
    fi
}

check_orca_swap() {
    local orca_path=$1
    [ -n "$orca_path" ] && [ -x "$orca_path" ] || return 0

    local min_gib total_kib min_kib total_gib
    min_gib=${NMR_ORCA_SWAP_MIN_GIB:-${NMR_CONTAINER_ORCA_SWAP_MIN_GIB:-0}}
    case "$min_gib" in
        ""|*[!0-9]*)
            fail "invalid ORCA swap floor: $min_gib"
            return 0
            ;;
    esac

    total_kib=$(swap_total_kib)
    min_kib=$((min_gib * 1024 * 1024))
    total_gib=$((total_kib / 1024 / 1024))

    if [ "$total_kib" -le 0 ]; then
        if is_enabled "${NMR_CONTAINER_REQUIRE_SWAP:-${NMR_DOCTOR_REQUIRE_SWAP:-0}}"; then
            fail "ORCA swap required but no swap is visible in /proc/swaps"
        else
            warn "ORCA executable is present but no swap is visible in /proc/swaps"
        fi
    elif [ "$min_kib" -gt 0 ] && [ "$total_kib" -lt "$min_kib" ]; then
        if is_enabled "${NMR_CONTAINER_REQUIRE_SWAP:-${NMR_DOCTOR_REQUIRE_SWAP:-0}}"; then
            fail "ORCA swap below requested floor: ${total_gib} GiB visible, ${min_gib} GiB requested"
        else
            warn "ORCA swap below requested floor: ${total_gib} GiB visible, ${min_gib} GiB requested"
        fi
    else
        ok "ORCA swap visible: ${total_gib} GiB"
    fi
}

first_existing_file() {
    for path in "$@"; do
        if [ -n "${path:-}" ] && [ -f "$path" ]; then
            printf '%s\n' "$path"
            return 0
        fi
    done
    return 1
}

toml_get() {
    local file=$1 dotted=$2
    [ -f "$file" ] || return 1
    have_cmd python3 || return 1
    python3 - "$file" "$dotted" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
dotted = sys.argv[2].split(".")
try:
    import tomllib
except ModuleNotFoundError:
    sys.exit(1)

try:
    data = tomllib.loads(path.read_text())
    value = data
    for key in dotted:
        value = value[key]
except Exception:
    sys.exit(1)

if isinstance(value, (str, int, float)):
    print(value)
else:
    sys.exit(1)
PY
}

first_nonempty() {
    for value in "$@"; do
        if [ -n "${value:-}" ]; then
            printf '%s\n' "$value"
            return 0
        fi
    done
    return 1
}

resolve_cmd_path() {
    local configured=$1 name=$2
    if [ -n "$configured" ]; then
        printf '%s\n' "$configured"
    elif have_cmd "$name"; then
        command -v "$name"
    fi
}

resolve_relative_to() {
    local base_file=$1 value=$2
    [ -n "$value" ] || return 1
    case "$value" in
        /*) printf '%s\n' "$value" ;;
        *) printf '%s/%s\n' "$(dirname -- "$base_file")" "$value" ;;
    esac
}

extract="${extract_arg:-${NMR_EXTRACT:-}}"
if [ -z "$extract" ]; then
    if have_cmd nmr_extract; then
        extract=$(command -v nmr_extract)
    elif [ -x "${ROOT_DIR}/build/nmr_extract" ]; then
        extract="${ROOT_DIR}/build/nmr_extract"
    fi
fi

if [ -n "$extract" ] && [ -x "$extract" ]; then
    ok "nmr_extract executable: $extract"
else
    fail "nmr_extract executable not found; pass a path or set NMR_EXTRACT"
fi

if [ -n "$extract" ] && [ -x "$extract" ] && have_cmd ldd; then
    missing=$(ldd "$extract" 2>/dev/null | awk '/not found/ {print}')
    if [ -n "$missing" ]; then
        fail "dynamic libraries missing:"
        printf '%s\n' "$missing" >&2
    else
        ok "dynamic linker reports no missing libraries"
    fi

    suspicious=$(ldd "$extract" 2>/dev/null |
        awk '/=>/ {print $3}' |
        grep -E "$LOCAL_PATH_RE" || true)
    if [ -n "$suspicious" ]; then
        leak "dynamic libraries resolve through machine-local paths:"
        printf '%s\n' "$suspicious" >&2
    fi
else
    warn "ldd unavailable or nmr_extract missing; skipped dynamic link check"
fi

if [ -n "$extract" ] && [ -x "$extract" ] && have_cmd readelf; then
    runpath=$(readelf -d "$extract" 2>/dev/null | grep 'RUNPATH' || true)
    if printf '%s\n' "$runpath" | grep -E "$LOCAL_RUNPATH_RE" >/dev/null; then
        leak "RUNPATH contains machine-local paths: $runpath"
    elif [ -n "$runpath" ]; then
        ok "RUNPATH present without obvious local path leak"
    else
        warn "no RUNPATH reported"
    fi
else
    warn "readelf unavailable or nmr_extract missing; skipped RUNPATH check"
fi

tools_toml="${NMR_TOOLS_TOML:-${HOME:-}/.nmr_tools.toml}"
if [ -f "$tools_toml" ]; then
    ok "runtime TOML: $tools_toml"
else
    warn "runtime TOML not found: $tools_toml"
fi

toml_mopac=$(toml_get "$tools_toml" mopac 2>/dev/null || true)
toml_tleap=$(toml_get "$tools_toml" tleap 2>/dev/null || true)
toml_gmx=$(toml_get "$tools_toml" gmx 2>/dev/null || true)
toml_orca=$(toml_get "$tools_toml" orca 2>/dev/null || true)
toml_ff14sb=$(toml_get "$tools_toml" ff14sb_params 2>/dev/null || true)
toml_tmpdir=$(toml_get "$tools_toml" tmpdir 2>/dev/null || true)
toml_bmrb=$(toml_get "$tools_toml" bmrb_atom_nom 2>/dev/null || true)
toml_larsen=$(toml_get "$tools_toml" larsen_hbond_grids 2>/dev/null || true)
toml_dsn=$(toml_get "$tools_toml" databases.tensorcs15 2>/dev/null || true)

mopac_path=$(resolve_cmd_path "$(first_nonempty "${NMR_MOPAC:-}" "$toml_mopac" || true)" mopac)
tleap_configured=$(first_nonempty "${NMR_TLEAP:-}" "$toml_tleap" || true)
if [ -z "$tleap_configured" ] && [ -n "${AMBERHOME:-}" ]; then
    tleap_configured="${AMBERHOME%/}/bin/tleap"
fi
tleap_path=$(resolve_cmd_path "$tleap_configured" tleap)
gmx_path=$(resolve_cmd_path "$(first_nonempty "${NMR_GMX:-}" "${NMR_GROMACS:-}" "$toml_gmx" || true)" gmx_mpi)
orca_path=$(resolve_cmd_path "$(first_nonempty "${NMR_ORCA:-}" "$toml_orca" || true)" orca)

check_optional_file "MOPAC binary" "$mopac_path"
check_optional_file "AmberTools tleap" "$tleap_path"
check_exec_file "GROMACS gmx_mpi" "$gmx_path"
check_exec_file "ORCA executable" "$orca_path"
check_orca_swap "$orca_path"

data_dir="${NMR_DATA_DIR:-}"
if [ -z "$data_dir" ]; then
    if [ -f "${ROOT_DIR}/data/ff14sb_params.dat" ]; then
        data_dir="${ROOT_DIR}/data"
    else
        data_dir="${ROOT_DIR}/share/nmr-shielding/data"
    fi
fi

ff14sb_path=$(first_nonempty "${NMR_FF14SB_PARAMS:-}" "$toml_ff14sb" "${data_dir}/ff14sb_params.dat" || true)
check_file "ff14SB params" "$ff14sb_path"

tmpdir=$(first_nonempty "${NMR_TMPDIR:-}" "$toml_tmpdir" "/tmp/nmr_shielding" || true)
check_temp_dir "temp directory" "$tmpdir"

bmrb_path=$(first_nonempty "${NMR_BMRB_ATOM_NOM:-}" "$toml_bmrb" || true)
check_optional_file "BMRB atom-name table" "$bmrb_path"

larsen_dir=$(first_nonempty "${NMR_LARSEN_HBOND_GRIDS:-}" "$toml_larsen" || true)
check_optional_dir "Larsen H-bond grids" "$larsen_dir"

dsn_present=0
if [ -n "${NMR_TENSORCS15_DSN:-}" ] || [ -n "$toml_dsn" ]; then
    dsn_present=1
fi
if [ "$dsn_present" -eq 1 ]; then
    ok "tensorcs15 DSN configured (value not printed)"
else
    warn "tensorcs15 DSN not configured"
fi

tensorcs15_manifest=$(first_existing_file \
    "${NMR_TENSORCS15_MANIFEST:-}" \
    "${ROOT_DIR}/config/tensorcs15.manifest.toml" \
    "${ROOT_DIR}/share/nmr-shielding/tensorcs15.manifest.toml" \
    "/etc/nmr-shielding/tensorcs15.manifest.toml" || true)
check_file "tensorcs15 manifest" "$tensorcs15_manifest"

calculator_config="${NMR_CALCULATOR_CONFIG:-}"
if [ -z "$calculator_config" ]; then
    if [ -f "${data_dir}/calculator_params.toml" ]; then
        calculator_config="${data_dir}/calculator_params.toml"
    elif [ -f "${ROOT_DIR}/data/calculator_params.toml" ]; then
        calculator_config="${ROOT_DIR}/data/calculator_params.toml"
    fi
fi

if [ -n "$calculator_config" ]; then
    check_file "calculator params" "$calculator_config"
else
    fail "calculator params not found"
fi

model_path="${NMR_AIMNET2_MODEL:-}"
if [ -z "$model_path" ] && [ -n "$calculator_config" ]; then
    model_value=$(toml_get "$calculator_config" aimnet2_model_path 2>/dev/null || true)
    model_path=$(resolve_relative_to "$calculator_config" "$model_value" || true)
fi
check_optional_file "AIMNet2 model" "$model_path"

if [ -n "${NMR_NVRTC_LIB_DIR:-}" ]; then
    nvrtc_dir=$NMR_NVRTC_LIB_DIR
elif [ -n "${NMR_NVIDIA_CU13_DIR:-}" ]; then
    nvrtc_dir="${NMR_NVIDIA_CU13_DIR%/}/lib"
elif have_cmd python3; then
    cu13_dir=$(python3 - <<'PY' 2>/dev/null || true
try:
    import nvidia.cu13
    print(list(nvidia.cu13.__path__)[0])
except Exception:
    pass
PY
)
    nvrtc_dir=${cu13_dir:+${cu13_dir}/lib}
else
    nvrtc_dir=""
fi

if [ -n "${nvrtc_dir:-}" ] && [ -d "$nvrtc_dir" ]; then
    if ls "$nvrtc_dir"/libnvrtc-builtins.so* >/dev/null 2>&1; then
        ok "cu13/nvrtc builtins directory: $nvrtc_dir"
    else
        fail "cu13 directory lacks libnvrtc-builtins: $nvrtc_dir"
    fi
else
    fail "cu13/nvrtc builtins directory not found"
fi

if [ "$deep" -eq 1 ]; then
    tensorcs15_checker=$(first_existing_file \
        "${ROOT_DIR}/scripts/tensorcs15-check.py" \
        "${SCRIPT_DIR}/nmr-tensorcs15-check" \
        "${SCRIPT_DIR}/tensorcs15-check.py" || true)
    if [ -z "$tensorcs15_checker" ] && have_cmd nmr-tensorcs15-check; then
        tensorcs15_checker=$(command -v nmr-tensorcs15-check)
    fi

    if [ "$dsn_present" -eq 1 ] && [ -n "$tensorcs15_checker" ] && have_cmd python3; then
        if python3 "$tensorcs15_checker" \
                --manifest "$tensorcs15_manifest" --quiet; then
            ok "tensorcs15 manifest check passed"
        else
            fail "tensorcs15 manifest check failed"
        fi
    else
        warn "tensorcs15 manifest check skipped; DSN, checker, or python3 not available"
    fi

    warn "--deep currently does not launch AIMNet2, MOPAC, or tleap; future opt-in probes go here"
fi

if [ "$failures" -gt 0 ]; then
    printf 'SUMMARY failed=%d warnings=%d\n' "$failures" "$warnings" >&2
    exit 1
fi

printf 'SUMMARY failed=0 warnings=%d\n' "$warnings"
exit 0
