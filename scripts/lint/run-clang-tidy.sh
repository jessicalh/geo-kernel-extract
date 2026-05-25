#!/usr/bin/env bash
# scripts/lint/run-clang-tidy.sh
#
# Full-tree clang-tidy sweep (or single-TU for incremental work).
# Emits JSON to artifacts/quality/clang-tidy-<sha>.json so downstream
# baseline-diff scripts can compute per-file finding counts.
#
# Usage:
#   scripts/lint/run-clang-tidy.sh                  # full sweep
#   scripts/lint/run-clang-tidy.sh src/Foo.cpp      # single TU
#
# Requires:
#   build/compile_commands.json (from CMAKE_EXPORT_COMPILE_COMMANDS=ON)
#   .clang-tidy at repo root
#   run-clang-tidy (LLVM-supplied parallel runner; from clang-tools)
#
# Parallel runner uses all cores by default. Override with -j arg
# passed through verbatim, e.g. `scripts/lint/run-clang-tidy.sh -j 8`.

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

if [[ ! -f build/compile_commands.json ]]; then
    echo "ERROR: build/compile_commands.json missing. Run cmake first." >&2
    exit 1
fi

SHA="$(git rev-parse --short HEAD)"
ARTIFACTS_DIR="artifacts/quality"
# Reporter, not ratchet: one overwritable report (SHA recorded inside the
# JSON, not in the filename) so runs don't accrete per-SHA dumps.
OUT_RAW="${ARTIFACTS_DIR}/clang-tidy-latest.txt"
OUT_JSON="${ARTIFACTS_DIR}/clang-tidy-latest.json"
mkdir -p "${ARTIFACTS_DIR}"

# Default parallelism: leave 8 cores of headroom on big boxes, all
# cores otherwise. Caller can override with -j on the command line.
TOTAL_CORES="$(nproc)"
DEFAULT_J=$(( TOTAL_CORES > 16 ? TOTAL_CORES - 8 : TOTAL_CORES ))

# Separate user args into:
#   PASSTHRU_ARGS — flags consumed by run-clang-tidy itself (e.g. -j)
#   TARGETS       — positional file/regex args
PASSTHRU_ARGS=()
TARGETS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|-j*|-extra-arg=*|-config-file=*|-clang-tidy-binary=*|-export-fixes=*|-fix|-format)
            PASSTHRU_ARGS+=("$1"); shift ;;
        *)
            TARGETS+=("$1"); shift ;;
    esac
done

# If no -j was passed in, supply our default.
if ! printf '%s\n' "${PASSTHRU_ARGS[@]}" | grep -qE '^-j'; then
    PASSTHRU_ARGS+=("-j" "${DEFAULT_J}")
fi

# Build file regex: when single-TU, regex matches anything ending in
# the given path (compile_commands.json uses absolute paths); when no
# args, regex matches our in-tree TUs only (src/tests/fileformat).
if [[ ${#TARGETS[@]} -gt 0 ]]; then
    # Strip ./ prefix if present, then build alternation matching the
    # absolute path's suffix.
    ESCAPED=()
    for t in "${TARGETS[@]}"; do
        t="${t#./}"
        ESCAPED+=("$(printf '%s' "$t" | sed 's,/,\\/,g; s,\.,\\.,g')")
    done
    FILE_REGEX="(${ESCAPED[0]}$"
    for ((i=1; i<${#ESCAPED[@]}; i++)); do
        FILE_REGEX+="|${ESCAPED[i]}$"
    done
    FILE_REGEX+=")"
else
    FILE_REGEX='^.*/(src|tests|fileformat)/[^/]+\.(cpp|cc)$'
fi

echo "clang-tidy parallel sweep at SHA ${SHA} (${DEFAULT_J}-way default)"

# run-clang-tidy emits per-TU output to stdout interleaved. -quiet
# suppresses the per-file progress chatter.
run-clang-tidy \
    -p build \
    -quiet \
    "${PASSTHRU_ARGS[@]}" \
    "${FILE_REGEX}" > "${OUT_RAW}" 2>&1 || true

# Convert the raw clang-tidy output to per-file per-check JSON.
python3 - <<EOF > "${OUT_JSON}"
import json, re, sys, os
root = os.getcwd()
warn_re = re.compile(
    r'^(?P<path>[^:]+):(?P<line>\d+):(?P<col>\d+): '
    r'(?P<level>warning|error): '
    r'(?P<msg>.*?) \[(?P<check>[a-zA-Z0-9_.,-]+)\]\s*$'
)
findings = []
with open("${OUT_RAW}") as f:
    for line in f:
        m = warn_re.match(line)
        if not m:
            continue
        path = m.group("path")
        if path.startswith(root + "/"):
            path = path[len(root)+1:]
        # Skip transitive header findings — only count in-tree user code.
        if not path.startswith(("src/", "tests/", "fileformat/")):
            continue
        findings.append({
            "file": path,
            "line": int(m.group("line")),
            "col":  int(m.group("col")),
            "level": m.group("level"),
            "check": m.group("check"),
            "msg":  m.group("msg"),
        })
# Deduplicate: same (file, line, col, check) reported by multiple
# parallel workers when one TU header-includes another.
seen = set()
deduped = []
for f in findings:
    key = (f["file"], f["line"], f["col"], f["check"])
    if key in seen:
        continue
    seen.add(key)
    deduped.append(f)
print(json.dumps({"sha": "${SHA}", "findings": deduped}, indent=2))
EOF

n=$(python3 -c "import json; print(len(json.load(open('${OUT_JSON}'))['findings']))")
echo "clang-tidy: ${n} findings (raw=${OUT_RAW}, json=${OUT_JSON})"
exit 0
