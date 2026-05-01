#!/usr/bin/env bash
# qwen_summary_batch.sh — run qwen_summary.py over every PDF in references/
# that lacks a canonical -summary.txt. Idempotent: skips papers whose
# -summary-qwen-test.txt already exists. Continues on errors.
#
# Usage:
#     scripts/references/qwen_summary_batch.sh
#
# Logs to /tmp/qwen-batch-<timestamp>.log (regenerable, kept out of the repo).
# SI files (basename ending in -SI) are skipped per WORKFLOW.md.

set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

ts="$(date +%Y%m%d-%H%M%S)"
LOG="/tmp/qwen-batch-${ts}.log"

PDFS_TMP="$(mktemp)"
SUMS_TMP="$(mktemp)"
trap 'rm -f "$PDFS_TMP" "$SUMS_TMP"' EXIT

ls references/*.pdf 2>/dev/null | xargs -n1 basename | sed 's/\.pdf$//' | sort > "$PDFS_TMP"
ls references-meta/*-summary.txt 2>/dev/null | xargs -n1 basename | sed 's/-summary\.txt$//' | sort > "$SUMS_TMP"

backlog=$(comm -23 "$PDFS_TMP" "$SUMS_TMP" | grep -v -- '-SI$' || true)
n_total=$(printf '%s\n' "$backlog" | grep -c . || true)

echo "$(date) batch start: $n_total papers in backlog" | tee -a "$LOG"
echo "log: $LOG"

ok=0; skip=0; fail=0

while IFS= read -r basename; do
    [[ -z "$basename" ]] && continue
    test_file="references-meta/${basename}-summary-qwen-test.txt"
    if [[ -s "$test_file" ]]; then
        echo "SKIP $basename (test summary exists)" | tee -a "$LOG"
        skip=$((skip + 1))
        continue
    fi
    echo "=== $basename ===" | tee -a "$LOG"
    if python3 scripts/references/qwen_summary.py "$basename" >>"$LOG" 2>&1; then
        ok=$((ok + 1))
    else
        rc=$?
        echo "FAIL $basename (exit $rc)" | tee -a "$LOG"
        fail=$((fail + 1))
    fi
done <<< "$backlog"

echo "$(date) batch done: ok=$ok skip=$skip fail=$fail of $n_total" | tee -a "$LOG"
