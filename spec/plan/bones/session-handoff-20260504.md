# Session handoff — 2026-05-04

**Pickup brief.** This session landed commit `b9636f0` (Phase 1
cleanup + bless v2 tolerance + AIMNet2 wired into smoke), then ran
a read-only review agent against it. The agent surfaced real findings;
the user reviewed, made decisions, and elected to wind down before
acting on them rather than burn more context. This handoff is the
TODO list for the next session. Reading order at the bottom.

This session also turned up the AIMNet2 weasel-monkey pattern (see
`feedback_aimnet2_required_no_weasel`): every prior builder that
"helpfully" tries to remove AIMNet2 because of the GPU/setup step
gets caught mid-flight by the user. Vigilance.

---

## Decisions made this session (carry into next)

### D1. AIMNet2 is required in ALL production paths, not just mutant + trajectory

Strengthens the existing `project_aimnet2_contract_20260426` memory.
The contract's earlier scope ("Other modes accept --aimnet2 but do
not require it; contract is at the production-pipeline level
(mutant + trajectory)") is superseded. Every mode that exercises
the calculator pipeline (`--pdb`, `--protonated-pdb`, `--orca`,
`--mutant`, `--trajectory`) requires AIMNet2 to load. Setup is a
one-time stable step on every machine; no excuse to weasel out.

### D2. Production baseline for tests = MOPAC-excluded only

Per `feedback_no_convenience_exclusions`. Every test that exercises
the calculator pipeline starts from the production baseline:
`PerFrameExtractionSet` for trajectory scope, the
calculator-stack-minus-MOPAC for conformation scope. The only
legitimate exclusion is MOPAC (cost ~92s on 1UBQ-class). APBS +
AIMNet2 + everything else stays on by default.

**Refinement (user 2026-05-04):** the one-shot test set
(`--protonated-pdb` + bare `--pdb` flow) needs ONE MOPAC test in
the whole set as a canary, so MOPAC silently breaking gets caught.
Existing `mopac_tests` executable likely covers this; verify.

### D3. Trajectory directory layout = Round-3 Option B fleet_amber

Canonical assumption for trajectory-related tests: layout looks
like `tests/data/fleet_amber/{1P9J_5801,1Z9B_6577}/prep_run_*/
batcave_local_15ns_optB_*/production` with `production.{tpr,trr,
xtc,edr}`. The location is in `testpaths.toml` keys
`fleet_amber_<id>_subpath`. **Non-applicable runs may exist next to
the canonical ones** (the `_backup_round1_*` and `_backup_round2_*`
dirs); the latest are the right ones. We can revisit the directory
convention later, but adapt to ONE shape at a time, not chase
multiple.

### D4. `aimnet2_charge_sensitivity` cruft removal — green-lit

Goes now. The catalog `LEGACY` entry in `_catalog.py:183` and the
test assertion in `test_load.py:483-487` (`test_aimnet2_specs_registered`)
that pins it both come out. **Aspirational replacement** (some
ensemble-Welford or autograd derivative of `aimnet2_charges` that
captures genuine sensitivity) is a separate thinking task; goal
isn't a like-for-like swap, it's understanding what the field was
trying to express. Handwavy until later.

### D5. Wind-down before more rm/edits

User flagged context-load this session ("we've backed up to the
spinner"). No rm/edits beyond saving state. Cleanup work — including
the agent's ranked TODOs — moves to the next session.

---

## TODO list for next session, sequenced

### P0 — small quick wins (≤30 min each)

- **EDR `md.edr` → `production.edr` fix.** `src/nmr_extract.cpp:239`
  hardcodes `md.edr` instead of using `spec.traj_edr` (which JobSpec
  set to `production.edr` and validated). Per D3, layout assumption
  is `production.edr`; fix is one-line: pass `spec.traj_edr`. Verify
  with a smoke run on the 1P9J or 1Z9B fixture.

- **`aimnet2_charge_sensitivity` removal** (D4):
  - Delete the `ArraySpec("aimnet2_charge_sensitivity", ...)` line
    from `python/nmr_extract/_catalog.py:183` (and the comment block
    above it).
  - Delete or revise `test_aimnet2_specs_registered` in
    `python/tests/test_load.py:483-487` to drop the legacy stem.
  - Drop the `AIMNet2ChargeSensitivity` import from
    `python/nmr_extract/_catalog.py:38` if no other consumer.
  - Re-run pytest — should still be clean.
  - Memory note: the aspirational replacement is a separate item,
    not a blocker.

### P1 — AIMNet2 always-required (D1) — production CLI path LANDED 2026-05-04 evening

**Landed in commit on 2026-05-04 evening (after the wind-down was
nominally over):**

- `src/JobSpec.cpp` `ValidateJobSpec`: every non-None mode now
  requires `aimnet2_model_path` to resolve (via `--aimnet2` CLI or
  `CalculatorConfig::GetString("aimnet2_model_path")` TOML fallback)
  AND for the file to exist. Hard-fail with a clear diagnostic
  pointing at the contract memory.
- `src/nmr_extract.cpp` main(): `CalculatorConfig::Load` now runs
  before `ValidateJobSpec` so the TOML fallback resolves. The
  AIMNet2 load block simplified to unconditional load (validation
  already verified the path).
- `tests/test_main.cpp`: now also loads `CalculatorConfig` from the
  project's `data/calculator_params.toml` so `JobSpecE2E.*` tests
  see the same TOML defaults the production binary does.
- `data/calculator_params.toml:83`: already had
  `aimnet2_model_path = "/shared/2026Thesis/nmr-shielding/data/models/aimnet2_wb97m_0.jpt"`
  — TOML default works on batcave; install env will need the same
  key (or pass `--aimnet2`).
- Memory: `feedback_aimnet2_required_no_weasel` (added in this
  session) captures the vigilance rule + the all-paths scope
  generalisation. The earlier `project_aimnet2_contract_20260426`
  is now partly superseded on scope (was "mutant + trajectory only";
  new contract is "every non-None mode") but the rest of that memory
  stays valid.

**Verification:** all 21 `SmokeTest|JobSpec*` tests pass after the
change. `JobSpecE2E.MutantEndToEnd`, `OrcaEndToEnd`, `PdbEndToEnd`,
`ProtonatedPdbEndToEnd` all run the production CLI path with the new
validation in place.

**Still TODO (P2 below) — per-test-site AIMNet2 wire-in:** the 18+
test sites that build `RunOptions` directly (bypassing JobSpec) still
don't load AIMNet2 into their RunOptions. Validation passes, but the
calculator pipeline inside those test bodies runs without AIMNet2.
That's the per-test-site issue and is the next piece of work. See
P2.

### P2 — Production baseline operationalisation (D2)

The 18+ test sites the agent flagged that build pipelines without
AIMNet2 — list verbatim:

- `tests/test_calculation_runner.cpp:38, 72, 124, 129` (3 sites)
- `tests/test_pipeline_and_sample.cpp:69, 104, 184, 249, 291, 326,
  393, 411` (8 sites)
- `tests/test_write_features.cpp:38` (1 site)
- `tests/test_job_spec.cpp:256, 285, 312, 346, 352` (5 sites; the
  five `JobSpecE2E` end-to-end tests)
- `tests/test_amber_streaming.cpp:138, 168, 225, 263, 303` (5 sites)

Plan:
- **Define the rule visibly.** A short note in
  `spec/TEST_FRAMEWORK.md` (or wherever the test tier doc lives)
  saying "production baseline = full calculator stack except MOPAC;
  one MOPAC canary per executable is fine; tests must opt-in to
  exclusions, not the other way around."
- **Build a shared helper** (e.g., `tests/ProductionTestSession.{h,
  cpp}` or extension to `TestEnvironment`) that returns a Session
  with AIMNet2 loaded, and a `MakeProductionRunOptions(session)`
  helper that returns RunOptions with the canonical baseline. Test
  sites use the helper instead of `RunOptions{}` literal.
- **Walk the 18+ sites** and convert. `JobSpecE2E.MutantEndToEnd`
  is the most contract-relevant — it's exactly what
  `project_aimnet2_contract_20260426` named.
- **`tests/test_amber_streaming.cpp`**: the trajectory tests
  currently `skip_mopac/coulomb/apbs/dssp = true`. Per D2 only MOPAC
  may be skipped; APBS + DSSP + AIMNet2 all stay on. RunConfiguration
  at the per-test sites needs revisiting.
- **Add ONE MOPAC canary** (D2 refinement). Likely
  `JobSpecE2E.PdbEndToEnd` or `ProtonatedPdbEndToEnd` — pick the
  fastest small-protein test in the one-shot set, leave MOPAC on.
  Document the canary intent in a comment.

### P3 — Test surface cleanup (Phase-1 leftovers)

- **`tests/golden/blessed/fleet/` retire.** 151 historical files
  tracked under `tests/golden/blessed/` despite the gitignore (it
  doesn't retroactively untrack). Use `git rm --cached -r tests/
  golden/blessed/fleet/` to drop from index. Commit separately so
  the diff is readable. The gitignore is already in place; this
  just realigns the index.
- **`tests/regression/run_regression.sh` retire.** Broken (Use
  Case B at `:46` calls `--orca "$ORCA_DIR"` without `--root NAME`
  which current JobSpec rejects) AND orphan (zero callers). Use
  Case D already in bones. Cleanest move: `git rm
  tests/regression/run_regression.sh` plus the `output_fleet/`
  artifact dir.
- **`tests/data/sdk_geo_only/1Q8K/` frames 002-006 retire.** 32
  NPYs each, no AIMNet2, no test consumes them. Stale orphans. Per
  D5 not now; in next session check whether any future test wants
  cross-frame validation; if not, delete.
- **Decide on `tests/data/external/1OKH_4587_protonated.pdb`.**
  `gmx_protonated`. FF intent flagged TBD in
  `spec/plan/test_inventory_2026-05-03.md:22` since 2026-05-03;
  consumed by `test_protonation_detection.cpp` and
  `test_foundation_results.cpp`. If CHARMM-era, retire to bones
  per the dead-letter argument; if AMBER, document the provenance
  and call it live.

### P4 — Polish

- **Tighten `min_npy_files` thresholds** in `tests/test_smoke.cpp`.
  Currently `nodft` requires ≥40 (actual 57); `withdft` requires
  ≥45 (actual 60). Set to ~95% of expected (55 / 58). Catches
  silent partial-output regressions.
- **Stop blessing `log.jsonl`** (or compare it with sensible
  timestamp normalisation). Currently dead weight in the bless
  dirs; future readers may think log content is contractually
  pinned when it isn't.
- **Add the "warn-once if bless_policy.toml missing" log line** to
  `tests/BlessCompare.cpp::LoadTable` per the comment in
  `BlessCompare.h:69` that promised it.
- **Hoist 6 PROPKA tests' shared skip condition** to suite-level.
  Six separate test cases skipping on the same `propka3` absence
  is duplication; one suite-level skip would cover it. Same shape
  for the 3 KaML tests.
- **Tombstone `project_iupac_topology_landed_20260426` memory.**
  The artifacts it claims (AtomTopology / IupacAtomName /
  AtomReference) were reverted 2026-04-27. Memory still claims
  they landed. Either remove the memory or add a stale-by-revert
  header.

### P5 — Aspirational / thinking tasks (not for next session)

- **What replaces `aimnet2_charge_sensitivity`?** D4 gave the
  removal green-light but flagged the question. The original was a
  perturbation method. Ensemble Welford on `aimnet2_charge` over
  MD frames is one direction; autograd-derivative on selected
  frames (already wired per the catalog comment) is another. What
  do we want this to MEAN scientifically before we wire something
  back in?
- **`MutationDeltaResult` dia/para completeness debt** per
  `project_dft_compare_calculator_completeness`. Memory says we
  should add WT total + mutant total + dia + para tensors plus three
  deltas. Untouched. Thesis-relevant correctness, not just a test
  concern.
- **Three-TOML-parser consolidation.** Low priority — they all work
  today — but they've already drifted. Lifting one shared helper
  (or vendoring toml++) would prevent the next surprise.
- **F6 systemic audit shape.** The contract-vs-code drift table in
  the agent's report (`spec/plan/test_review_20260504.md` §4) is
  one finding. The pattern likely recurs as the project accumulates
  more contract memories. Worth a periodic sweep tied to
  doc-cleanup, or a "verify before write" rule for new contract
  memories.

---

## Reading order for next session

1. **This doc** (`spec/plan/session-handoff-20260504.md`) —
   start here.
2. **`spec/plan/test_review_20260504.md`** — the agent's full
   review report. The agent worked file:line-cited; the report is
   the source of truth for what's where.
3. **`tests/TEST_HEALTH.md`** — current test status snapshot
   (commit `b9636f0`). F1-F6 findings + how-to-refresh.
4. **`tests/golden/blessed/BLESS_NOTES.md`** — bless v2 history.
5. **`spec/plan/session-handoff-20260503.md`** — yesterday's plan,
   which `b9636f0` executed.
6. **Memories (auto-loaded):** `project_aimnet2_contract_20260426`,
   `feedback_aimnet2_required_no_weasel` (new), the
   `feedback_no_convenience_exclusions`, plus the standard reading
   order in `spec/INDEX.md`.

---

## What this session left UNDONE

Nothing in production state is broken that wasn't broken before.
The audit surfaced gaps; the user chose to wind down before
acting. The `b9636f0` checkin landed cleanly, tests are green
(347/347 ctest, 64/68 pytest, all skips with reasons), bless v2
+ AIMNet2-in-smoke are durable. The next session has a clean
starting point and a sequenced TODO.

The only thing genuinely deferred-with-impact is the JobSpec
mutant `--aimnet2` gate (P1). Until it lands, a
`nmr_extract --mutant ...` run without `--aimnet2` produces
contract-violating output silently. Smoke catches it; the CLI
binary doesn't. Land this early next session.
