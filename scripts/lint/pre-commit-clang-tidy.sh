#!/usr/bin/env bash
# scripts/lint/pre-commit-clang-tidy.sh
#
# Pre-commit hook: run clang-tidy on each changed C++ TU, count its
# warnings, and compare to the per-file count in
# artifacts/quality/baseline.json. Fail (= block the commit) on any
# file whose count exceeds its baseline anchor.
#
# Called by pre-commit with the list of staged C++ files (after the
# files: regex filter in .pre-commit-config.yaml).
#
# Side-effect-free: reads baseline.json, runs clang-tidy, no writes.
#
# Counting rule: we count only warnings whose SITE is in the changed
# file itself. Transitive-header warnings (from src/Foo.h included by
# src/Bar.cpp) are attributed to Foo.h in baseline.json, not to Bar.cpp,
# so the per-file ratchet stays apples-to-apples.

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

if [[ $# -eq 0 ]]; then
    exit 0
fi

if [[ ! -f build/compile_commands.json ]]; then
    echo "ERROR: build/compile_commands.json missing. Run cmake first." >&2
    exit 1
fi

if [[ ! -f artifacts/quality/baseline.json ]]; then
    echo "ERROR: artifacts/quality/baseline.json missing." >&2
    echo "Bootstrap with:" >&2
    echo "  scripts/lint/run-all-cpp.sh --skip-scan-build" >&2
    echo "  scripts/lint/baseline.py --update" >&2
    exit 1
fi

ROOT="$(pwd)"

# Run clang-tidy per file. For each TU, count only warnings whose
# warning-site file matches the TU path itself (not transitively-
# included headers — those are attributed to the header in baseline.json).
declare -A CURRENT
for tu in "$@"; do
    # `|| true` on the grep: empty match (TU has 0 in-TU warnings — the
    # happy case after a successful clean-up pass) returns exit 1, which
    # the `set -o pipefail` above would otherwise propagate and kill the
    # script before we ever reach the comparison.
    n=$( { clang-tidy -p build --quiet "$tu" 2>/dev/null \
            | grep -E "^${ROOT}/${tu}:[0-9]+:[0-9]+: warning:" || true; } \
        | wc -l)
    CURRENT["$tu"]="$n"
done

# Compare to baseline.json.
fail=0
while IFS=$'\t' read -r tu current_n; do
    base_n=$(python3 -c "
import json
b = json.load(open('artifacts/quality/baseline.json'))
print(b['tools'].get('clang-tidy', {}).get('${tu}', {}).get('total', 0))
")
    if (( current_n > base_n )); then
        echo "  REGRESSION: ${tu}  baseline=${base_n}  current=${current_n}  (+$((current_n - base_n)))" >&2
        fail=1
    fi
done < <(for tu in "${!CURRENT[@]}"; do
    printf "%s\t%s\n" "$tu" "${CURRENT[$tu]}"
done)

if (( fail )); then
    echo "" >&2
    echo "Pre-commit clang-tidy ratchet failed — fix the new finding(s)," >&2
    echo "or re-anchor with deliberate review:" >&2
    echo "  scripts/lint/run-clang-tidy.sh && scripts/lint/baseline.py --update" >&2
    exit 1
fi

echo "clang-tidy ratchet OK ($# file(s) checked)."
exit 0
