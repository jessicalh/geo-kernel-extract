#!/usr/bin/env bash
# scripts/lint/run-iwyu.sh
#
# include-what-you-use sweep — flags includes that bring in more than
# what's used, and missing includes that work only via transitive
# pull-in. IWYU is opinionated and often wrong on Eigen / HighFive
# (template-heavy headers); use as a HINT generator, not a hard gate.
#
# Usage:
#   scripts/lint/run-iwyu.sh                  # full sweep
#   scripts/lint/run-iwyu.sh src/Foo.cpp ...
#
# Requires:
#   include-what-you-use + iwyu_tool installed
#   build/compile_commands.json

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

if [[ ! -f build/compile_commands.json ]]; then
    echo "ERROR: build/compile_commands.json missing. Run cmake first." >&2
    exit 1
fi

SHA="$(git rev-parse --short HEAD)"
ARTIFACTS_DIR="artifacts/quality"
OUT_TXT="${ARTIFACTS_DIR}/iwyu-latest.txt"   # single overwritable report
OUT_JSON="${ARTIFACTS_DIR}/iwyu-latest.json"
mkdir -p "${ARTIFACTS_DIR}"

# iwyu_tool consumes compile_commands.json and dispatches IWYU per TU.
# --output-format=clang gives one-finding-per-line format easy to parse.
# Per-TU we want to lint src/ tests/ fileformat/ only, not _deps,
# extern/, or build-generated TUs.

if [[ $# -gt 0 ]]; then
    iwyu_tool -p build "$@" > "${OUT_TXT}" 2>&1 || true
else
    # Build a target list from the compile DB.
    mapfile -t TARGETS < <(python3 -c '
import json, os
root = os.getcwd()
with open("build/compile_commands.json") as f:
    db = json.load(f)
for entry in db:
    src = entry["file"]
    if "/_deps/" in src or "/build/" in src or "/extern/" in src:
        continue
    rel = os.path.relpath(src, root)
    if rel.startswith(("src/", "tests/", "fileformat/")):
        print(rel)
' | sort -u)
    iwyu_tool -p build "${TARGETS[@]}" > "${OUT_TXT}" 2>&1 || true
fi

# Convert IWYU output to a JSON summary of per-file "should add" /
# "should remove" counts.
python3 - <<EOF > "${OUT_JSON}"
import json, os, re
root = os.getcwd()
findings = []
current_file = None
mode = None  # "add" or "remove" within a per-file block
with open("${OUT_TXT}") as f:
    for line in f:
        m = re.match(r'^(.+?) should add these lines:\s*$', line)
        if m:
            current_file = m.group(1)
            if current_file.startswith(root + "/"):
                current_file = current_file[len(root)+1:]
            mode = "add"
            continue
        m = re.match(r'^(.+?) should remove these lines:\s*$', line)
        if m:
            current_file = m.group(1)
            if current_file.startswith(root + "/"):
                current_file = current_file[len(root)+1:]
            mode = "remove"
            continue
        if current_file and mode and line.strip() and not line.startswith("---") and not line.startswith("The full include-list"):
            findings.append({
                "file": current_file,
                "check": "iwyu-" + mode,
                "msg": line.strip(),
            })
print(json.dumps({"sha": "${SHA}", "findings": findings}, indent=2))
EOF

n=$(python3 -c "import json; print(len(json.load(open('${OUT_JSON}'))['findings']))")
echo "IWYU: ${n} include-hygiene hints (txt=${OUT_TXT}, json=${OUT_JSON})"
exit 0
