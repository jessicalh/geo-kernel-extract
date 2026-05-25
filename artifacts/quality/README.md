# artifacts/quality/ — static-analysis reports

C++ static-analysis is a **reporting tool, not a ratchet** — see
`CLEANUP_WORKFLOW.md`. No commit-gate, no `baseline.json` anchor (that
machinery was removed 2026-05-25). Everything here is regenerable and
gitignored except this README.

## Files (all gitignored, one overwritable report per tool)

- `clang-tidy-latest.{txt,json}` — raw clang-tidy output + parsed JSON.
- `cppcheck-latest.{xml,json}` — cppcheck XML report + parsed JSON.
- `iwyu-latest.{txt,json}` — include-what-you-use output + parsed JSON.
- `scan-build-latest/` — scan-build HTML bug-report tree + JSON summary.

The committed SHA the report was generated at is recorded *inside* the
JSON (`"sha"` field), not in the filename — so reruns overwrite rather
than accrete.

## Run

```
scripts/lint/run-clang-tidy.sh [src/Foo.cpp]
scripts/lint/run-cppcheck.sh   [src/Foo.cpp]
scripts/lint/run-iwyu.sh       [src/Foo.cpp]
scripts/lint/run-scan-build.sh
scripts/lint/run-all-cpp.sh    [--skip-scan-build]
```

Findings are triaged by a human and fixed in reviewed chunks; `-fix` is
not used.
