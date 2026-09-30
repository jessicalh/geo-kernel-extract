#!/bin/sh
# Copy only the carried container engine to a local executable path.
set -eu
bundle=$1
workspace=$2
mode=${3:-0}
fail() { printf '%s\n' "$*" >&2; exit 1; }
case "$workspace" in *[[:space:]]*) fail 'Bundled Apptainer needs a local workspace path without whitespace. Choose such a workspace or use compatibility mode.' ;; esac
identity=$(cat "$bundle/runtime/apptainer.identity")
case "$identity" in ''|*[!0-9a-f]*) fail 'Invalid bundled engine identity.' ;; esac
[ "${#identity}" -eq 64 ] || fail 'Invalid bundled engine identity.'
engine="$workspace/data/apptainer-$identity"
[ ! -L "$engine" ] || fail 'The bundled engine cache must not be a symbolic link.'
if [ -f "$engine/.complete" ]; then printf '%s\n' "$engine/bin/apptainer"; exit 0; fi
[ ! -e "$engine" ] || fail "A previous bundled-engine cache is incomplete: $engine"
[ "$mode" != check ] || fail 'The bundled engine is available. Open Start Reader - Apptainer.desktop once to prepare its local application files, then rerun this check.'
message="Prepare Apptainer on this computer — once

Copy about 163 MB of container-engine software to:
$workspace/data

The local engine can run even when the external drive blocks program execution.

The Reader application image and all trajectories stay on your drive. No trajectory files are copied. Later launches reuse this small local engine."
if [ "$mode" != 1 ]; then
    if command -v zenity >/dev/null 2>&1; then
        env GSK_RENDERER=cairo LIBGL_ALWAYS_SOFTWARE=1 zenity --question --title='Prepare Apptainer — once' --width=560 --ok-label='Set up Apptainer' --cancel-label=Cancel --text="$message" || exit 130
    elif command -v kdialog >/dev/null 2>&1; then
        env LIBGL_ALWAYS_SOFTWARE=1 kdialog --title='Prepare Apptainer — once' --yes-label='Set up Apptainer' --no-label=Cancel --yesno "$message" || exit 130
    elif command -v xmessage >/dev/null 2>&1; then
        xmessage -center -buttons 'Set up Apptainer:0,Cancel:1' "$message" || exit 130
    elif [ -t 0 ]; then
        printf '%s\nSet up Apptainer? [y/N] ' "$message" >&2
        read -r answer
        case "$answer" in y|Y|yes|YES) ;; *) exit 130 ;; esac
    else
        fail 'Bundled Apptainer setup needs a confirmation dialog or interactive terminal.'
    fi
fi
mkdir -p "$workspace/data"
temporary=$(mktemp -d "$workspace/data/.apptainer-XXXXXXXX")
trap 'rm -rf -- "$temporary"' EXIT
trap 'exit 130' HUP INT TERM
cp -R --preserve=mode,timestamps -- "$bundle/runtime/apptainer/." "$temporary/"
"$temporary/bin/apptainer" --version >/dev/null || fail 'This carried engine needs an EL8-compatible host (glibc 2.28+) and executable local storage. Try the site Singularity launcher or compatibility mode.'
: > "$temporary/.complete"
if ! mv -T --no-clobber -- "$temporary" "$engine"; then
    [ ! -L "$engine" ] && [ -f "$engine/.complete" ] || fail 'Could not publish the local engine cache.'
fi
if [ -d "$temporary" ]; then
    [ ! -L "$engine" ] && [ -f "$engine/.complete" ] || fail 'The bundled engine cache changed during setup; retry with a fresh local workspace.'
    rm -rf -- "$temporary"
fi
trap - EXIT HUP INT TERM
printf '%s\n' "$engine/bin/apptainer"
