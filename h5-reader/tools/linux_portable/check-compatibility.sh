#!/bin/sh
# Small, local-only compatibility check. Never opens trajectory datasets.
set -eu
unset LD_LIBRARY_PATH LD_PRELOAD
ulimit -c 0 2>/dev/null || true
umask 077
package_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
workspace=${H5READER_WORKSPACE:-${XDG_DATA_HOME:-$HOME/.local/share}/h5reader-portable}
report=
backend=auto
show_report=0
usage() {
    printf '%s\n' 'Usage: /bin/sh check-compatibility.sh [--package-dir DIR] [--workspace DIR] [--backend auto|singularity|apptainer] [--report FILE] [--show-report]'
}
fail_early() { printf 'Reader compatibility check: %s\n' "$*" >&2; exit 1; }
while [ "$#" -gt 0 ]; do
    case "$1" in
        --package-dir|--workspace|--report|--backend)
            [ "$#" -ge 2 ] || fail_early "Missing value for $1"
            case "$1" in
                --package-dir) package_dir=$2 ;;
                --workspace) workspace=$2 ;;
                --report) report=$2 ;;
                --backend) backend=$2 ;;
            esac
            shift 2 ;;
        --show-report) show_report=1; shift ;;
        --help|-h) usage; exit 0 ;;
        *) usage >&2; fail_early "Unknown option: $1" ;;
    esac
done
case "$backend" in auto|apptainer|singularity) ;; *) fail_early 'Unknown container backend.' ;; esac
package_dir=$(CDPATH= cd -- "$package_dir" && pwd -P) || fail_early 'Package directory is unavailable.'
source_root=$(realpath -m -- "$package_dir/../..")
case "$workspace" in /*) ;; *) fail_early 'Workspace must be an absolute local path.' ;; esac
workspace=$(realpath -m -- "$workspace")
case "$workspace/" in "$source_root/"*) fail_early 'Workspace must be outside the source drive.' ;; esac
case "$source_root/" in "$workspace/"*) fail_early 'Workspace must not contain the source drive.' ;; esac
case "$workspace$package_dir" in *:*|*,*) fail_early 'Paths cannot contain a colon or comma.' ;; esac
state_dir=$(realpath -m -- "$workspace/state")
temporary_dir=$(realpath -m -- "$workspace/tmp")
case "$state_dir/" in "$workspace/"*) ;; *) fail_early 'Workspace state directory points outside the workspace.' ;; esac
case "$temporary_dir/" in "$workspace/"*) ;; *) fail_early 'Workspace temporary directory points outside the workspace.' ;; esac
mkdir -p -- "$state_dir" "$temporary_dir" || fail_early 'Cannot write to the local workspace.'
if [ -z "$report" ]; then
    report="$state_dir/compatibility-$(date +%Y%m%d-%H%M%S)-$$.txt"
fi
case "$report" in /*) ;; *) fail_early 'Report path must be absolute.' ;; esac
report=$(realpath -m -- "$report")
case "$report/" in "$source_root/"*) fail_early 'The report must be saved outside the source drive.' ;; esac
check_dir=$(mktemp -d "$temporary_dir/compatibility.XXXXXXXX") || fail_early 'Cannot create the small local check directory.'
cleanup() { rm -rf -- "$check_dir"; }
trap cleanup EXIT
trap 'exit 130' HUP INT TERM
exec 3>&1
exec > "$report" 2>&1
printf '%s\n' 'H5 Reader compatibility check' '============================'
printf 'Date: %s\n' "$(date -Iseconds)"
printf 'Package: %s\nLocal workspace: %s\n' "$package_dir" "$workspace"
printf '%s\n\n' 'Only synthetic test files are used. No trajectory files are opened or copied. This report stays on this computer.'
printf '%s\n' 'Host facts (host glibc is recorded, not used as the Reader requirement):'
uname -smr
if [ -r /etc/os-release ]; then
    sed -n '/^PRETTY_NAME=/p; /^ID=/p; /^VERSION_ID=/p' /etc/os-release
fi
getconf GNU_LIBC_VERSION 2>/dev/null || true
if command -v getenforce >/dev/null 2>&1; then printf 'SELinux: '; getenforce || true; fi
printf 'Desktop session: %s\nDISPLAY: %s\n' "${XDG_SESSION_TYPE:-unknown}" "${DISPLAY:-not set}"
xauthority_file=${XAUTHORITY:-$HOME/.Xauthority}
if [ -r "$xauthority_file" ]; then
    printf 'X11 authorization file: readable (%s); cookie contents are not recorded.\n' "$xauthority_file"
else
    printf 'X11 authorization file: unavailable (%s); the actual Qt launch will test access.\n' "$xauthority_file"
fi
printf '\n'
finish() {
    result=$1
    summary=$2
    printf '\n%s\n' "$summary"
    printf '%s\n' 'This checks container startup, read-only binding, local writes and Qt/X11 startup. Full rendering, inference, and your actual drive still need normal Reader acceptance.'
    printf '\nReport: %s\n' "$report"
    printf '%s\nReport: %s\n' "$summary" "$report" >&3
    if [ "$show_report" = 1 ] && [ -n "${DISPLAY:-}" ] && command -v zenity >/dev/null 2>&1; then
        env GSK_RENDERER=cairo LIBGL_ALWAYS_SOFTWARE=1 zenity --text-info \
            --title='Reader compatibility report' --width=820 --height=620 \
            --filename="$report" >/dev/null 2>&1 || true
    elif [ "$show_report" = 1 ] && [ -n "${DISPLAY:-}" ] && command -v kdialog >/dev/null 2>&1; then
        env LIBGL_ALWAYS_SOFTWARE=1 kdialog --title='Reader compatibility report' \
            --textbox "$report" 820 620 >/dev/null 2>&1 || true
    elif [ "$show_report" = 1 ] && [ -n "${DISPLAY:-}" ] && command -v xmessage >/dev/null 2>&1; then
        xmessage -title 'Reader compatibility report' -geometry 820x620 -file "$report" >/dev/null 2>&1 || true
    elif [ "$show_report" = 1 ] && [ -n "${DISPLAY:-}" ] && command -v xdg-open >/dev/null 2>&1; then
        xdg-open "$report" >/dev/null 2>&1 || true
    elif [ "$show_report" = 1 ]; then
        cat "$report" >&3
        if [ -t 0 ]; then
            printf '\nPress Enter to close this report. ' >&3
            read -r ignored || true
        fi
    fi
    exit "$result"
}
[ "$(uname -s)" = Linux ] && [ "$(uname -m)" = x86_64 ] || \
    finish 1 'FAIL: this package requires x86_64 Linux.'
image="$package_dir/payload/reader.sif"
[ -r "$image" ] || finish 1 'FAIL: payload/reader.sif is missing or unreadable.'
runtime=
if [ "$backend" = singularity ]; then
    runtime=$(command -v singularity || true)
elif [ "$backend" = apptainer ]; then
    runtime=$(command -v apptainer || true)
elif command -v apptainer >/dev/null 2>&1; then
    runtime=$(command -v apptainer)
elif command -v singularity >/dev/null 2>&1; then
    runtime=$(command -v singularity)
fi
if [ -z "$runtime" ]; then
    printf 'Requested container backend: %s\n' "$backend"
    printf '%s\n' 'Apptainer/Singularity was not found on PATH.' \
        'Your site may provide it through environment modules. Load the site-provided Apptainer or Singularity module in a terminal, then rerun this check there.' \
        'No shell profile was sourced and no module or system package was installed.'
    [ -z "${MODULESHOME:-}" ] || printf 'MODULESHOME: %s\n' "$MODULESHOME"
    [ -z "${LOADEDMODULES:-}" ] || printf 'Loaded modules: %s\n' "$LOADEDMODULES"
    finish 1 'NEEDS SITE SETUP: Apptainer/Singularity is not available in this session. The fallback runtime was not installed.'
fi
printf 'Container command: %s\n' "$runtime"
version_text=$("$runtime" --version 2>&1) || finish 1 'FAIL: the installed container command cannot report its version.'
printf '%s\n' "$version_text"
version=$(printf '%s\n' "$version_text" | sed -n 's/^[^0-9]*\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p')
major=${version%%.*}
minor=${version#*.}
case "$version" in
    ''|*[!0-9.]*) finish 1 'FAIL: could not identify the container version.' ;;
esac
case "$version_text" in
    *[Aa]pptainer*) [ "$major" -ge 1 ] || finish 1 'FAIL: Apptainer 1.0 or newer is required.' ;;
    *) [ "$major" -gt 3 ] || { [ "$major" -eq 3 ] && [ "$minor" -ge 7 ]; } || \
        finish 1 'FAIL: Singularity 3.7 or newer is required by the portable launcher.' ;;
esac
if ! "$runtime" exec --help > "$check_dir/exec-help.txt" 2>&1; then
    cat "$check_dir/exec-help.txt"
    finish 1 'FAIL: container exec is unavailable.'
fi
for flag in --cleanenv --containall --no-home --no-mount --pwd --workdir --bind --env; do
    if ! grep -F -q -- "$flag" "$check_dir/exec-help.txt"; then
        finish 1 "FAIL: the installed container command does not advertise required option $flag."
    fi
done
printf '%s\n' 'PASS: container version and required command options.'
# --no-nv / --no-rocm are hidden options in supported older releases; the real
# execution below verifies them instead of requiring their presence in help.
unset APPTAINER_BIND APPTAINER_BINDPATH APPTAINER_MOUNT
unset SINGULARITY_BIND SINGULARITY_BINDPATH SINGULARITY_MOUNT
unset APPTAINER_NV APPTAINER_NVCCLI APPTAINER_ROCM
unset SINGULARITY_NV SINGULARITY_NVCCLI SINGULARITY_ROCM
synthetic_source="$check_dir/source"
synthetic_work="$check_dir/work"
mkdir "$synthetic_source" "$synthetic_work"
printf '%s\n' reader-compatibility-marker > "$synthetic_source/marker"
for directory in home config cache data state run tmp output; do
    mkdir "$synthetic_work/$directory"
done
mkdir -p "$synthetic_work/tmp/apptainer/tmp/.X11-unix"
tr -d '-' < /proc/sys/kernel/random/uuid > "$synthetic_work/state/machine-id"
if [ -r "$xauthority_file" ]; then cp -- "$xauthority_file" "$synthetic_work/run/Xauthority"; fi
export APPTAINER_CACHEDIR="$synthetic_work/cache/apptainer" APPTAINER_TMPDIR="$synthetic_work/tmp"
export SINGULARITY_CACHEDIR="$synthetic_work/cache/apptainer" SINGULARITY_TMPDIR="$synthetic_work/tmp"
run_container() {
    set -- "$image" "$@"
    if [ -n "${DISPLAY:-}" ] && [ -d /tmp/.X11-unix ]; then
        set -- --bind /tmp/.X11-unix:/tmp/.X11-unix:ro \
            --env "DISPLAY=$DISPLAY" --env XAUTHORITY=/workspace/run/Xauthority "$@"
    fi
    # Every bind is a small synthetic directory or the display socket. No
    # scientific-data path is passed to the runtime.
    timeout 90s "$runtime" exec --cleanenv --containall --no-home --no-nv --no-rocm \
        --no-mount hostfs,cwd,home,sys --pwd /workspace/output \
        --workdir "$synthetic_work/tmp/apptainer" \
        --bind "$synthetic_source:/provenance:ro" --bind "$synthetic_work:/workspace:rw" \
        --bind "$synthetic_work/state/machine-id:/etc/machine-id:ro" "$@"
}
printf '\n%s\n' 'Checking the packaged userspace and synthetic mounts (90-second limit):'
if run_container /bin/sh -ec '
    /bin/true
    printf "Packaged userspace: "; getconf GNU_LIBC_VERSION
    test "$(cat /provenance/marker)" = reader-compatibility-marker
    if (printf forbidden >> /provenance/marker) 2>/dev/null; then
        echo "FAIL: synthetic source mount is writable."; exit 1
    fi
    printf local-write-ok > /workspace/output/probe
    test "$(cat /workspace/output/probe)" = local-write-ok
    rm /workspace/output/probe
    echo "PASS: image executes; synthetic source is read-only; local workspace is writable."
'; then
    :
else
    status=$?
    printf 'Container exit status: %s\n' "$status"
    finish 1 'FAIL: the image could not pass startup and mount checks. The runtime error above identifies the failure; ask the site administrator about its container policy.'
fi
if [ -z "${DISPLAY:-}" ]; then
    finish 2 'PARTIAL PASS: container and mounts work. Run this check from the advisor desktop to test Qt/X11; DISPLAY is not set in this session.'
fi
printf '\n%s\n' 'Checking the actual Reader binary and Qt/X11 (version-only launch, no trajectory opens):'
if run_container /opt/h5reader/guest-start.sh --version; then
    finish 0 'PASS: packaged userspace, read-only source binding, local workspace and Reader Qt/X11 startup work on this computer.'
else
    status=$?
    printf 'Reader startup exit status: %s\n' "$status"
    finish 1 'FAIL: the container works, but Reader Qt/X11 startup failed. See the exact loader/display error above; check the desktop session and X11 authorization.'
fi
