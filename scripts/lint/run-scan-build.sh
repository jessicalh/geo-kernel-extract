#!/usr/bin/env bash
# scripts/lint/run-scan-build.sh
#
# Clang static analyzer (scan-build) sweep. Runs the analyzer as part
# of a fresh CMake build, intercepting compiler invocations. Catches
# path-sensitive bugs (use-after-move, null deref, unreachable code,
# integer overflow, etc.) that clang-tidy's syntactic checks miss.
#
# This is the slowest tool — full rebuild required, single-threaded
# analyzer on each TU. Plan ~30 min on a clean repo.
#
# Output is HTML at artifacts/quality/scan-build-<sha>/index.html plus
# a JSON summary derived from the HTML.

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

SHA="$(git rev-parse --short HEAD)"
ARTIFACTS_DIR="artifacts/quality"
# single overwritable report + reused build tree (no per-SHA accretion)
OUT_DIR="${ARTIFACTS_DIR}/scan-build-latest"
OUT_JSON="${ARTIFACTS_DIR}/scan-build-latest.json"
BUILD_DIR="build-scan-latest"
mkdir -p "${ARTIFACTS_DIR}"
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

SCAN_BUILD="${SCAN_BUILD:-scan-build-18}"

echo "scan-build sweep at SHA ${SHA} (this takes a while...)"

# scan-build wraps cmake AND make so analyzer state is captured for
# every TU compiled. --use-cc/cxx forces the wrapper compiler that
# emits the analyzer plists.
${SCAN_BUILD} \
    -o "${OUT_DIR}" \
    --use-cc=clang \
    --use-c++=clang++ \
    --status-bugs \
    --keep-going \
    cmake -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Debug \
                            -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
                            . >> "${ARTIFACTS_DIR}/scan-build-cmake-latest.log" 2>&1 || true

${SCAN_BUILD} \
    -o "${OUT_DIR}" \
    --use-cc=clang \
    --use-c++=clang++ \
    --status-bugs \
    --keep-going \
    cmake --build "${BUILD_DIR}" -- -j"$(nproc)" 2>> "${ARTIFACTS_DIR}/scan-build-build-latest.log" || true

# scan-build creates a dated subdir under OUT_DIR. Flatten by counting.
n_html=$(find "${OUT_DIR}" -name "*.html" -type f 2>/dev/null | grep -cE "report-[a-z0-9]+\.html" || echo 0)
echo "scan-build: ${n_html} bug reports in HTML"

# Generate a minimal JSON summary by scraping the HTML index.
python3 - <<EOF > "${OUT_JSON}"
import json, os, re, glob
out_dir = "${OUT_DIR}"
findings = []
# Each dated subdir holds an index.html with a table of bugs.
for idx in glob.glob(os.path.join(out_dir, "*/index.html")):
    with open(idx) as f:
        html = f.read()
    # Rough scrape: each bug row has class "bt", file in <td>foo.cpp</td>.
    rows = re.findall(
        r'class="bt"[^>]*>(?P<category>[^<]*)</td>.*?'
        r'<td[^>]*>(?P<file>[^<]+\.(?:cpp|h|cc|hpp))</td>.*?'
        r'<td[^>]*>(?P<line>\d+)</td>',
        html, re.S)
    for cat, path, line in rows:
        findings.append({
            "file":  path.strip(),
            "line":  int(line),
            "level": "warning",
            "check": "scan-build:" + cat.strip(),
            "msg":   cat.strip(),
        })
print(json.dumps({"sha": "${SHA}", "findings": findings}, indent=2))
EOF

# Clean up: keep HTML for spot-investigation, remove the throwaway
# build tree (we rebuild for each scan-build run).
rm -rf "${BUILD_DIR}"

n_json=$(python3 -c "import json; print(len(json.load(open('${OUT_JSON}'))['findings']))")
echo "scan-build: ${n_json} findings (html=${OUT_DIR}, json=${OUT_JSON})"
exit 0
