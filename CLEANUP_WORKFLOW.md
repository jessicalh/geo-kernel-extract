# C++ static analysis — reporters, not a ratchet

Lint runs against this tree as a **reporting tool**: you run it, you read
the findings, you fix what's worth fixing. There is **no commit-gate, no
pre-commit hook, and no baseline ratchet** (that machinery was removed
2026-05-25 — it welded lint into version control and made messes). The
build itself carries `-Wall -Wextra -Wpedantic`; these tools go further.

## Run

```
scripts/lint/run-clang-tidy.sh            # full sweep (parallel)
scripts/lint/run-clang-tidy.sh src/X.cpp  # single TU
scripts/lint/run-cppcheck.sh              # full sweep
scripts/lint/run-iwyu.sh                  # include-what-you-use (hint)
scripts/lint/run-scan-build.sh            # deep analyzer (slow; full rebuild)
scripts/lint/run-all-cpp.sh               # all of the above
```

Each writes a single, overwritable report to `artifacts/quality/`
(`<tool>-latest.{txt,json,xml}`, gitignored). The committed SHA is
recorded *inside* the JSON, not in the filename, so runs don't accrete.

## Check set

`.clang-tidy` is a **W3-discipline allowlist** (MSVC terms): bug detectors
only — bugprone / clang-analyzer / concurrency + a couple of
cppcoreguidelines/performance checks. Style, modernize, naming, readability,
unused-includes, magic-numbers are explicitly out of scope. cppcheck runs
at `warning,style,unusedFunction` (note: `unusedFunction` is false-positive
-prone across TU / subproject / virtual boundaries — treat with skepticism).

## Discipline

Findings are triaged by a human, fixed in reviewed chunks, never auto-applied
(`-fix` is not used) and never block a commit.
