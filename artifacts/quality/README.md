# artifacts/quality/ — static-analysis baseline + per-sha reports

This directory holds the C++ static-analysis ratchet anchor and per-sha
tool reports.

## Files

### Committed (the contract)

- `baseline.json` — per-file finding counts for each tool, anchored at
  a specific SHA. The pre-commit hook compares against this; new
  findings in a changed file beyond its baseline entry block the
  commit. Re-anchor with `scripts/lint/baseline.py --update` only
  after deliberate review (cleanup work, intentional accepted churn).
- `README.md` — this file.

### Gitignored (regenerable)

- `clang-tidy-<sha>.{txt,json}` — raw clang-tidy output + parsed JSON.
- `cppcheck-<sha>.{xml,json}` — cppcheck XML report + parsed JSON.
- `iwyu-<sha>.{txt,json}` — include-what-you-use output + parsed JSON.
- `scan-build-<sha>/` — scan-build HTML bug-report tree, plus
  `scan-build-<sha>.json` summary. Large; clean up old SHAs to save
  space.
- `scan-build-cmake-<sha>.log`, `scan-build-build-<sha>.log` — logs
  from the scan-build run's cmake + build phases.

## Workflow

```
# Run a full C++ sweep (clang-tidy + cppcheck + IWYU + scan-build):
scripts/lint/run-all-cpp.sh

# Faster: skip scan-build (which rebuilds + runs the analyzer per TU):
scripts/lint/run-all-cpp.sh --skip-scan-build

# Single TU:
scripts/lint/run-clang-tidy.sh src/Foo.cpp
scripts/lint/run-cppcheck.sh   src/Foo.cpp
scripts/lint/run-iwyu.sh       src/Foo.cpp

# Re-anchor baseline.json from the just-emitted reports:
scripts/lint/baseline.py --update

# Verify (used by pre-commit):
scripts/lint/baseline.py --check                    # all files
scripts/lint/baseline.py --check src/Foo.cpp ...    # changed files only
```

## Ratchet philosophy

The baseline is the floor, not the target. The pre-commit gate
ENFORCES "do not regress" — a changed file may not finish a commit
with more findings than its baseline entry. The TARGET is zero, and
cleanup work whittles the baseline down opportunistically:

1. Pick a file with high baseline count.
2. Fix some/all of its findings.
3. Re-run the relevant tool: `scripts/lint/run-clang-tidy.sh src/X.cpp`
4. Re-anchor: `scripts/lint/baseline.py --update`
5. Commit the fixes + the baseline.json change in the same commit so
   the new floor moves down atomically.

Re-anchoring UPWARD (deliberately accepting more findings) is allowed
but should be rare — only when a refactor surfaces previously-hidden
patterns OR when a check is reclassified after review. Document in
the commit message why the baseline rose.

## Tools active in the ratchet

| Tool                   | Strength                                          | Notes                                                              |
|------------------------|---------------------------------------------------|--------------------------------------------------------------------|
| clang-tidy 18          | Syntactic + semantic lint, many check families    | `.clang-tidy` at repo root; ~100 families enabled                  |
| cppcheck 2.13          | Flow analysis (memleak, uninit, overflow paths)   | `--enable=all --inconclusive`; consumes `compile_commands.json`    |
| include-what-you-use   | Header hygiene (missing / extra includes)         | Opinionated on template-heavy headers; treat as HINT not gate      |
| clang scan-build 18    | Path-sensitive bug analysis (use-after-move, etc) | Slowest; full rebuild per run; HTML + JSON summary                 |

## Future expansion

- Python (Phase 2): ruff + mypy + bandit → `artifacts/quality/ruff-*.json` etc.
- R (Phase 2): lintr
- Universal (Phase 3): shellcheck, codespell, markdownlint
