#!/usr/bin/env bash
# scripts/lint/run-cppcheck.sh
#
# Full-tree cppcheck sweep. Complementary to clang-tidy: cppcheck's
# flow analysis catches a different bug class (esp. memory leaks,
# uninit reads, integer overflow paths). Emits XML to
# artifacts/quality/cppcheck-<sha>.xml.
#
# Usage:
#   scripts/lint/run-cppcheck.sh           # full sweep
#   scripts/lint/run-cppcheck.sh src/Foo.cpp ...
#
# Requires:
#   cppcheck installed
#   build/compile_commands.json

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

if [[ ! -f build/compile_commands.json ]]; then
    echo "ERROR: build/compile_commands.json missing. Run cmake first." >&2
    exit 1
fi

SHA="$(git rev-parse --short HEAD)"
ARTIFACTS_DIR="artifacts/quality"
OUT_XML="${ARTIFACTS_DIR}/cppcheck-${SHA}.xml"
OUT_JSON="${ARTIFACTS_DIR}/cppcheck-${SHA}.json"
mkdir -p "${ARTIFACTS_DIR}"

# Standard cppcheck check set:
#   --enable=all: warnings/style/perf/portability/info/missingInclude
#   --inconclusive: allow lower-confidence findings (we triage)
#   --suppress=missingIncludeSystem: don't pester about std headers
#   --suppress=unusedFunction: many fn defined for ABI we don't always see
#   --suppress=unmatchedSuppression: harmless meta-suppress
#
# Pin to C++17 explicitly to match CMakeLists; without it cppcheck
# sometimes defaults too high or low.
#
# --project consumes compile_commands.json directly so include paths
# + defines come from the real build.

if [[ $# -gt 0 ]]; then
    CPPCHECK_ARGS=(--enable=all --inconclusive --std=c++17
                   --suppress=missingIncludeSystem
                   --suppress=unusedFunction
                   --suppress=unmatchedSuppression
                   --xml --xml-version=2
                   -I src -I fileformat
                   "$@")
else
    CPPCHECK_ARGS=(--enable=all --inconclusive --std=c++17
                   --suppress=missingIncludeSystem
                   --suppress=unusedFunction
                   --suppress=unmatchedSuppression
                   --xml --xml-version=2
                   --project=build/compile_commands.json
                   "-i${PWD}/build"
                   "-i${PWD}/extern")
fi

echo "cppcheck sweep at SHA ${SHA}..."
# cppcheck XML goes to stderr; reroute to file.
cppcheck "${CPPCHECK_ARGS[@]}" 2> "${OUT_XML}" || true

# Convert XML to a JSON list of findings (file/line/severity/id/msg).
python3 - <<EOF > "${OUT_JSON}"
import json, os, xml.etree.ElementTree as ET
root = os.getcwd()
tree = ET.parse("${OUT_XML}")
findings = []
for err in tree.iter("error"):
    sev = err.get("severity", "")
    cid = err.get("id", "")
    msg = err.get("msg", "")
    for loc in err.findall("location"):
        path = loc.get("file", "")
        if path.startswith(root + "/"):
            path = path[len(root)+1:]
        # Filter out paths from build/ or external (paranoia in case
        # cppcheck didn't honor the -i filter).
        if path.startswith(("build/", "extern/", "/usr/")):
            continue
        findings.append({
            "file": path,
            "line": int(loc.get("line", "0") or 0),
            "col":  int(loc.get("column", "0") or 0),
            "level": sev,
            "check": cid,
            "msg":  msg,
        })
print(json.dumps({"sha": "${SHA}", "findings": findings}, indent=2))
EOF

n=$(python3 -c "import json; print(len(json.load(open('${OUT_JSON}'))['findings']))")
echo "cppcheck: ${n} findings (xml=${OUT_XML}, json=${OUT_JSON})"
exit 0
