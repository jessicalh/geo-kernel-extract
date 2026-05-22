# Cleanup-phase workflow — C++ static analysis ratchet

Bootstrapped 2026-05-22 (post-baseline-closeout, sha `76e7068`).
Live record of how lint / static-analysis runs against this tree
and how the no-regression ratchet works.

## Tools at a glance

| Tool             | Trigger              | Output                                         |
|------------------|----------------------|------------------------------------------------|
| clang-format 18  | `git clang-format`, pre-commit on staged hunks | format-on-write fixes in place         |
| clang-tidy 18    | manual sweep + pre-commit ratchet | `artifacts/quality/clang-tidy-<sha>.json` |
| cppcheck 2.13    | manual sweep         | `artifacts/quality/cppcheck-<sha>.xml/.json`   |
| IWYU 0.21        | manual sweep (hint, not gate) | `artifacts/quality/iwyu-<sha>.txt/.json` |
| scan-build 18    | manual sweep (slow; full rebuild) | `artifacts/quality/scan-build-<sha>/` (HTML+JSON) |

## Daily commit flow

```
# When you stage C++ changes:
git add src/MyChange.cpp src/MyChange.h
git commit -m "..."     # pre-commit hook fires automatically:
                        #   - trailing whitespace stripped
                        #   - clang-format on staged hunks
                        #   - clang-tidy on changed files, count
                        #     compared to artifacts/quality/baseline.json;
                        #     blocks if any file count exceeded its anchor

# If blocked: fix the new findings, re-stage, re-commit.
# If the new findings are deliberate (refactor surfaces a pattern):
scripts/lint/baseline.py --update      # re-anchors the ratchet at NEW state
git add artifacts/quality/baseline.json
git commit ...
```

## Cleanup work — whittling down the baseline

Pick a target file with high finding count, fix some/all, re-anchor:

```
# 1. See top-offender files in the current baseline:
python3 -c "
import json
b = json.load(open('artifacts/quality/baseline.json'))
ct = b['tools']['clang-tidy']
for path, d in sorted(ct.items(), key=lambda kv: -kv[1]['total'])[:20]:
    print(f'  {d[\"total\"]:4d}  {path}')
"

# 2. See WHICH checks fire on a specific file:
clang-tidy -p build src/Foo.cpp 2>/dev/null | \
    grep -oP '\[[a-z]+-[a-z0-9-]+\]' | sort | uniq -c | sort -rn

# 3. Fix some findings. clang-tidy can auto-fix a subset:
clang-tidy -p build --fix --fix-errors src/Foo.cpp     # CAREFUL; review diff

# 4. Re-run + re-anchor:
scripts/lint/run-clang-tidy.sh src/Foo.cpp
scripts/lint/baseline.py --update

# 5. Commit fixes + baseline.json change together so the floor moves
#    down atomically.
```

## Full sweep (run periodically; baseline is anchored from a full sweep)

```
scripts/lint/run-all-cpp.sh --skip-scan-build   # ~25 min on 268 TUs
scripts/lint/baseline.py --update               # re-anchors baseline

# Full sweep WITH scan-build (slowest tool, adds ~30 min):
scripts/lint/run-all-cpp.sh
```

## Adding scan-build to the baseline (one-off)

Scan-build wasn't in the original 2026-05-22 baseline because its
~30 min rebuild cost is high for a bootstrap. When ready:

```
scripts/lint/run-scan-build.sh
scripts/lint/baseline.py --update          # picks up the new tool
git add artifacts/quality/baseline.json
git commit -m "Extend ratchet baseline to include scan-build"
```

## Bypass (don't, except for the cases below)

```
git commit --no-verify       # skips the pre-commit gate
```

Legitimate uses:
- The baseline gate is wrong and you're committing the fix (e.g.
  the gate script itself).
- An emergency hotfix where the lint issue is unrelated and time is
  short — file a follow-up to fix the regression.

If you find yourself reaching for `--no-verify` more than once or
twice in a session, that's a smell: the rules are wrong, the
baseline is wrong, or the workflow is fighting you. Stop and
re-anchor instead.

## Phase 2/3 (not yet bootstrapped)

- **Phase 2 (Python)**: ruff (format + lint with most rule families
  enabled) + mypy (strict on `nmr_extract` SDK, looser on
  `learn/` scripts) + bandit (security). Adds
  `artifacts/quality/ruff-*.json`, `mypy-*.json`, `bandit-*.json`.

- **Phase 3 (Universal)**: shellcheck on `scripts/*.sh`, codespell on
  `**/*.{cpp,h,py,md}`, markdownlint on `spec/`, `doc/`, root
  `*.md`. Adds `artifacts/quality/{shellcheck,codespell,markdown}-*.json`.

## Drifting away from the ratchet

If the project grows out of this ratchet (different compiler version,
new check families released, etc.):

```
# Compare what's enabled now vs what would be enabled with new clang-tidy:
clang-tidy --list-checks -p build > /tmp/now.txt
# ... upgrade ...
clang-tidy --list-checks -p build > /tmp/new.txt
diff /tmp/now.txt /tmp/new.txt
```

Then update `.clang-tidy` deliberately and re-anchor.

## Related files

- `.clang-format`, `.clang-tidy`, `.pre-commit-config.yaml`,
  `CMakePresets.json` — at repo root.
- `scripts/lint/*` — runners + baseline tool + pre-commit shim.
- `artifacts/quality/` — reports + baseline.json + README.
- `artifacts/quality/baseline.json` — committed; the ratchet anchor.
