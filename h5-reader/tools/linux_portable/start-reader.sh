#!/bin/sh
# Offline launcher. It reads the source drive and extracts only onto local storage.
set -eu
unset LD_LIBRARY_PATH LD_PRELOAD
fail() {
    printf '%s\n' "H5 Reader: $*" >&2
    if [ "${H5READER_PORTABLE_NO_DIALOG:-0}" != 1 ] && \
            command -v zenity >/dev/null 2>&1 && [ -n "${DISPLAY:-}" ]; then
        env GSK_RENDERER=cairo LIBGL_ALWAYS_SOFTWARE=1 zenity --error --title='H5 Reader' --text="$*" 2>/dev/null || true
    fi
    exit 1
}
bundle_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
source_root=${H5READER_SOURCE_ROOT:-}
workspace=${H5READER_WORKSPACE:-${XDG_DATA_HOME:-$HOME/.local/share}/h5reader-portable}
backend=${H5READER_PORTABLE_BACKEND:-auto}
noninteractive_setup=${H5READER_PORTABLE_NONINTERACTIVE_SETUP:-0}
while [ "$#" -gt 0 ]; do
    case "$1" in
        --source-root) [ "$#" -ge 2 ] || fail 'Missing source root'; source_root=$2; shift 2 ;;
        --workspace) [ "$#" -ge 2 ] || fail 'Missing workspace'; workspace=$2; shift 2 ;;
        --backend) [ "$#" -ge 2 ] || fail 'Missing backend'; backend=$2; shift 2 ;;
        --non-interactive-setup) noninteractive_setup=1; shift ;;
        --) shift; break ;;
        *) break ;;
    esac
done
[ "$(uname -m)" = x86_64 ] || fail 'This Reader package requires x86_64 Linux.'
[ -n "${DISPLAY:-}" ] || fail 'Start Reader from an X11 or XWayland desktop session.'
if [ -z "$source_root" ]; then
    source_root=$(CDPATH= cd -- "$bundle_dir/../.." && pwd -P)
fi
[ -d "$source_root" ] || fail 'The source drive is not available.'
source_root=$(CDPATH= cd -- "$source_root" && pwd -P)
case "$workspace" in /*) ;; *) fail 'The workspace must be an absolute path.' ;; esac
case "$source_root$workspace" in *:*|*,*) fail 'Source and workspace paths cannot contain a colon or comma.' ;; esac
# Resolve existing ancestors before creating anything: a workspace link into the
# source drive must not turn a read-only browsing session into archive writes.
workspace=$(realpath -m -- "$workspace")
case "$workspace/" in "$source_root/"*) fail 'Choose a writable workspace outside the source drive.' ;; esac
case "$source_root/" in "$workspace/"*) fail 'Choose a workspace that does not contain the source drive.' ;; esac
mkdir -p -- "$workspace" || fail 'Cannot create the local workspace.'
workspace=$(CDPATH= cd -- "$workspace" && pwd -P)
case "$workspace/" in "$source_root/"*) fail 'The workspace resolves inside the source drive.' ;; esac
case "$source_root/" in "$workspace/"*) fail 'The workspace contains the source drive.' ;; esac
for directory in home config cache data state run tmp output runtime; do
    resolved=$(realpath -m -- "$workspace/$directory")
    case "$resolved/" in "$workspace/"*) ;; *) fail 'A workspace directory points outside the workspace.' ;; esac
    mkdir -p -- "$workspace/$directory" || fail 'Cannot prepare the local workspace.'
done
chmod 700 "$workspace/run"
if [ ! -s "$workspace/state/machine-id" ]; then
    tr -d '-' < /proc/sys/kernel/random/uuid > "$workspace/state/machine-id"
fi
payload="$bundle_dir/payload"
container_runtime=
case "$backend" in auto|apptainer|singularity|proot) ;; *) fail 'Unknown runtime backend.' ;; esac
if [ "$backend" = singularity ]; then
    container_runtime=$(command -v singularity || true)
elif [ "$backend" = apptainer ]; then
    container_runtime=$(command -v apptainer || true)
elif [ "$backend" = auto ]; then
    if command -v apptainer >/dev/null 2>&1; then container_runtime=$(command -v apptainer)
    elif command -v singularity >/dev/null 2>&1; then container_runtime=$(command -v singularity)
    fi
fi
if [ -n "$container_runtime" ]; then
    [ -f "$payload/reader.sif" ] || fail 'The Apptainer image is missing from this package.'
    xauthority_file=${XAUTHORITY:-$HOME/.Xauthority}
    if [ -r "$xauthority_file" ]; then
        cp -- "$xauthority_file" "$workspace/run/Xauthority"
    fi
    if [ "${H5READER_PORTABLE_KEEP_STDIO:-0}" != 1 ]; then
        exec >> "$workspace/state/reader.log" 2>&1
    fi
    # Suppress inherited container options, including GPU injection and extra
    # host binds. The only scientific-data mount is explicitly read-only.
    unset APPTAINER_BIND APPTAINER_BINDPATH APPTAINER_MOUNT
    unset SINGULARITY_BIND SINGULARITY_BINDPATH SINGULARITY_MOUNT
    unset APPTAINER_NV APPTAINER_NVCCLI APPTAINER_ROCM
    unset SINGULARITY_NV SINGULARITY_NVCCLI SINGULARITY_ROCM
    export APPTAINER_CACHEDIR="$workspace/cache/apptainer" APPTAINER_TMPDIR="$workspace/tmp"
    export SINGULARITY_CACHEDIR="$workspace/cache/apptainer" SINGULARITY_TMPDIR="$workspace/tmp"
    # --containall replaces /tmp with this work directory. Create the socket
    # bind target there before Apptainer applies the nested X11 bind.
    mkdir -p "$workspace/tmp/apptainer/tmp/.X11-unix"
    printf '%s\n' 'H5READER_PORTABLE_BACKEND=apptainer' >&2
    if "$container_runtime" exec --cleanenv --containall --no-home \
        --no-mount hostfs,bind-paths,cwd,home,sys --pwd /workspace/output \
        --workdir "$workspace/tmp/apptainer" \
        --bind "$source_root:/provenance:ro" --bind "$workspace:/workspace:rw" \
        --bind "$workspace/state/machine-id:/etc/machine-id:ro" \
        --bind /tmp/.X11-unix:/tmp/.X11-unix:ro \
        --env "DISPLAY=$DISPLAY" --env XAUTHORITY=/workspace/run/Xauthority \
        "$payload/reader.sif" /opt/h5reader/guest-start.sh "$@"; then
        exit 0
    fi
    fail "Apptainer could not run Reader. Details are in $workspace/state/reader.log. The optional PRoot fallback is available with --backend proot."
elif [ "$backend" = apptainer ] || [ "$backend" = singularity ]; then
    fail 'Apptainer/Singularity is not available on this computer.'
fi
[ -r "$payload/SHA256SUMS" ] || fail 'The portable runtime is incomplete.'
runtime_id=$(sed -n 's/  reader-userspace.tar.gz$//p' "$payload/SHA256SUMS")
case "$runtime_id" in ''|*[!0-9a-f]*) fail 'The runtime checksum manifest is invalid.' ;; esac
[ "${#runtime_id}" -eq 64 ] || fail 'The runtime checksum manifest is invalid.'
runtime="$workspace/runtime/$runtime_id"
if [ ! -f "$runtime/.complete" ]; then
    setup_text="Reader is on your external drive. Without Apptainer/Singularity, this computer needs a compatible application runtime on local storage.

One-time setup: unpack about 1.3 GB of application software to:
$workspace

Later launches reuse it. Your trajectories stay on the external drive and are opened directly. No trajectories will be copied."
    progress_kind=none
    if [ "$noninteractive_setup" != 1 ]; then
        if command -v zenity >/dev/null 2>&1; then
            env GSK_RENDERER=cairo LIBGL_ALWAYS_SOFTWARE=1 zenity --question \
                --title='Set up Reader on this computer — once' --width=560 \
                --ok-label='Set up Reader' --cancel-label=Cancel --text="$setup_text" || exit 130
            progress_kind=zenity
        elif command -v kdialog >/dev/null 2>&1; then
            env LIBGL_ALWAYS_SOFTWARE=1 kdialog --title='Set up Reader on this computer — once' \
                --yes-label='Set up Reader' --no-label=Cancel --yesno "$setup_text" || exit 130
            progress_kind=kdialog
        elif command -v xmessage >/dev/null 2>&1; then
            xmessage -center -buttons 'Set up Reader:0,Cancel:1' "$setup_text" || exit 130
            progress_kind=xmessage
        elif [ -t 0 ]; then
            printf '%s\nSet up Reader? [y/N] ' "$setup_text" >&2
            read -r answer
            case "$answer" in y|Y|yes|YES) ;; *) exit 130 ;; esac
        else
            fail 'First-time fallback setup requires a confirmation dialog. Run this launcher from a terminal, or use a desktop with zenity, kdialog or xmessage.'
        fi
    fi
    temporary="$workspace/runtime/.extract-$$"
    mkdir "$temporary" || fail 'Cannot create runtime extraction directory.'
    worker_pid= progress_pid=
    cleanup_setup() {
        [ -z "$worker_pid" ] || kill "$worker_pid" 2>/dev/null || true
        [ -z "$progress_pid" ] || kill "$progress_pid" 2>/dev/null || true
        [ -z "$worker_pid" ] || wait "$worker_pid" 2>/dev/null || true
        [ -z "$progress_pid" ] || wait "$progress_pid" 2>/dev/null || true
        rm -rf -- "$temporary"
    }
    trap cleanup_setup EXIT
    trap 'exit 130' HUP INT TERM
    case "$progress_kind" in
        zenity)
            mkfifo "$temporary/progress"
            env GSK_RENDERER=cairo LIBGL_ALWAYS_SOFTWARE=1 zenity --progress --pulsate \
                --auto-close --title='Setting up Reader application' --width=520 \
                --text='Checking and unpacking the application. Trajectories stay on the external drive.' \
                < "$temporary/progress" &
            progress_pid=$!
            exec 3> "$temporary/progress"
            ;;
        kdialog)
            env LIBGL_ALWAYS_SOFTWARE=1 kdialog --title='Setting up Reader application' --ok-label=Cancel \
                --msgbox 'Checking and unpacking the application. Trajectories stay on the external drive. Close this message to cancel.' &
            progress_pid=$! ;;
        xmessage)
            xmessage -center -buttons Cancel:1 'Setting up Reader application. Checking and unpacking application files; trajectories stay on the external drive.' &
            progress_pid=$! ;;
    esac
    finish_worker() {
        while kill -0 "$worker_pid" 2>/dev/null; do
            if [ -n "$progress_pid" ] && ! kill -0 "$progress_pid" 2>/dev/null; then
                kill "$worker_pid" 2>/dev/null || true
                wait "$worker_pid" 2>/dev/null || true
                worker_pid=
                exit 130
            fi
            sleep 0.2
        done
        if wait "$worker_pid"; then worker_pid=; else worker_pid=; return 1; fi
    }
    # Detect blocked ptrace / noexec before reading and unpacking the large
    # application archive. The test executes only the host's /bin/true.
    proot_expected=$(sed -n 's/  proot-x86_64$//p' "$payload/SHA256SUMS")
    proot_actual=$(sha256sum "$payload/proot-x86_64")
    [ "${proot_actual%% *}" = "$proot_expected" ] || fail 'The PRoot checksum failed.'
    cp "$payload/proot-x86_64" "$temporary/proot"
    chmod 755 "$temporary/proot"
    "$temporary/proot" --version >/dev/null 2>&1 || fail 'The workspace must permit execution of Linux programs.'
    export PROOT_TMP_DIR="$workspace/tmp"
    "$temporary/proot" -r / -w / /bin/true >/dev/null 2>&1 || \
        fail 'This computer blocks the fallback runtime (ptrace). No application archive was unpacked.'
    (cd "$payload" && exec sha256sum --check --status SHA256SUMS) &
    worker_pid=$!
    finish_worker || fail 'The portable runtime checksum failed.'
    mkdir "$temporary/rootfs"
    tar --extract --gzip --file "$payload/reader-userspace.tar.gz" \
        --directory "$temporary/rootfs" --no-same-owner &
    worker_pid=$!
    finish_worker || fail 'Cannot unpack the portable runtime.'
    if [ "$progress_kind" = zenity ]; then
        printf '100\n' >&3
        exec 3>&-
        wait "$progress_pid" 2>/dev/null || true
    elif [ -n "$progress_pid" ]; then
        kill "$progress_pid" 2>/dev/null || true
        wait "$progress_pid" 2>/dev/null || true
    fi
    progress_pid=
    rm -f "$temporary/progress"
    : > "$temporary/.complete"
    if [ -e "$runtime" ]; then
        [ -f "$runtime/.complete" ] || fail 'A previous runtime preparation is incomplete.'
        rm -rf -- "$temporary"
    else
        mv -- "$temporary" "$runtime"
    fi
    trap - EXIT HUP INT TERM
fi
if [ "${H5READER_PORTABLE_KEEP_STDIO:-0}" != 1 ]; then
    exec >> "$workspace/state/reader.log" 2>&1
fi
printf '%s\n' 'H5READER_PORTABLE_BACKEND=proot' >&2
unset LD_LIBRARY_PATH LD_PRELOAD PROOT_LOADER PROOT_LOADER_32
export PROOT_TMP_DIR="$workspace/tmp"
# Deliberately do not bind /dev, /sys, the user's home, or the host library tree.
# Mesa's bundled CPU renderer needs no GPU devices. PRoot does not provide a
# security boundary. Reader opens scientific data for reading in this workflow.
xauthority_file=${XAUTHORITY:-$HOME/.Xauthority}
if [ -r "$xauthority_file" ]; then
    cp -- "$xauthority_file" "$workspace/run/Xauthority"
    export XAUTHORITY=/workspace/run/Xauthority
fi
# Probe ptrace before opening the GUI; managed RHEL/SELinux policy can restrict it.
"$runtime/proot" -r "$runtime/rootfs" /bin/true >/dev/null 2>&1 || \
    fail 'This Linux system blocks the portable runtime (ptrace). Its administrator may need to allow PRoot.'
exec "$runtime/proot" -r "$runtime/rootfs" -w /workspace/output \
    -b "$source_root:/provenance" -b "$workspace:/workspace" \
    -b "$workspace/state/machine-id:/etc/machine-id" \
    -b /proc:/proc -b /dev/null:/dev/null -b /dev/zero:/dev/zero \
    -b /dev/urandom:/dev/urandom -b /dev/random:/dev/random \
    -b "$workspace/tmp:/tmp" -b "$workspace/run:/dev/shm" \
    -b /tmp/.X11-unix:/tmp/.X11-unix /opt/h5reader/guest-start.sh "$@"
