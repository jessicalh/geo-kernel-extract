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
# Reporter, not ratchet: one overwritable report (SHA recorded inside the
# JSON, not in the filename) so runs don't accrete per-SHA dumps.
OUT_XML="${ARTIFACTS_DIR}/cppcheck-latest.xml"
OUT_JSON="${ARTIFACTS_DIR}/cppcheck-latest.json"
mkdir -p "${ARTIFACTS_DIR}"

# W3 discipline (mirrors .clang-tidy 2026-05-24): bugs + UB + dead code
# only; style / performance / portability preferences out of scope.
#
#   --enable=warning           : real-bug-class findings cppcheck is
#                                confident about (knownConditionTrueFalse,
#                                negativeContainerIndex, syntaxError,
#                                uninit reads, memleak, resourceLeak, ...)
#   --enable=style             : needed because cppcheck's dead-code
#                                checks (unreadVariable, redundantAssignment,
#                                unreachableCode, unusedStructMember,
#                                unusedScopedObject, unusedLabel) sit in
#                                the "style" severity bucket. Noisy style
#                                preferences from this bucket are
#                                suppressed individually below.
#   --enable=unusedFunction    : cross-TU dead-function detector. The
#                                only cppcheck check that crosses TU
#                                boundaries; clang-tidy cannot do this.
#                                Previously suppressed by mistake.
#   (no --inconclusive)        : W3 = high-confidence only. The
#                                inconclusive flag was producing noise
#                                we triaged repeatedly without acting on.
#   --suppress=uninitMemberVar : codebase initializes POD fields at use
#                                site, not declaration. Same rationale as
#                                .clang-tidy's disable of
#                                cppcoreguidelines-pro-type-member-init.
#   --suppress=...             : pure-preference style IDs that fired
#                                in the previous baseline without
#                                indicating a bug.
#
# Pin to C++17 to match CMakeLists. --project consumes compile_commands.json
# directly so include paths + defines come from the real build.

CPPCHECK_BASE=(
    --enable=warning,style,unusedFunction
    --std=c++17
    --suppress=missingIncludeSystem
    --suppress=unmatchedSuppression
    # Init-at-use-site is codebase convention.
    --suppress=uninitMemberVar
    # Style preferences below — not bugs, not dead code.
    --suppress=useStlAlgorithm
    --suppress=noExplicitConstructor
    --suppress=constVariable
    --suppress=constVariableReference
    --suppress=constVariablePointer
    --suppress=constParameterReference
    --suppress=constParameterPointer
    --suppress=constParameter
    --suppress=functionStatic
    --suppress=functionConst
    --suppress=funcArgNamesDifferent
    --suppress=passedByValue
    --suppress=uselessCallsSubstr
    --suppress=shadowFunction
    --suppress=variableScope
    --suppress=cstyleCast
    --suppress=postfixOperator
    --xml --xml-version=2
)

if [[ $# -gt 0 ]]; then
    CPPCHECK_ARGS=("${CPPCHECK_BASE[@]}" -I src -I fileformat "$@")
else
    CPPCHECK_ARGS=("${CPPCHECK_BASE[@]}"
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
        # Only keep findings that live inside the project tree —
        # cppcheck walks the full include closure (reduce/, gromacs/,
        # libtorch/, openbabel/, miniforge3/...) and reports
        # unusedFunction on every vendored header. Those aren't our
        # debt to manage. Mirrors clang-tidy's HeaderFilterRegex.
        if not path.startswith(("src/", "tests/", "fileformat/")):
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
