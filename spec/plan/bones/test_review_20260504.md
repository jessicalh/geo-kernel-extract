# Test Infrastructure + Contract Review — 2026-05-04

Read-only review by general-purpose agent against commit `b9636f0`
("Phase 1 cleanup + bless v2 tolerance + AIMNet2 wired into smoke").
Brief in `spec/plan/session-handoff-20260504.md`. ~11 minutes /
128 tool uses.

User responses to the findings (recorded inline in the conversation
that produced this artifact) are summarised at the end of
`session-handoff-20260504.md`. The report itself is reproduced
verbatim below.

---

## Overall verdict

The Phase-1 cleanup landed coherently; bless v2 + AIMNet2-in-smoke
close the two named gaps cleanly. But the contract-without-enforcement
audit (F6) is bigger than the AIMNet2 instance: at least three other
test paths run pipelines without AIMNet2 in violation of the same
memory and the no-convenience-exclusions rule, and `--mutant`/`--orca`
don't actually require `--aimnet2 MODEL` despite the contract memory
saying they do. The bless v2 contract itself is solid; gaps are
around what's NOT comparing (log.jsonl, structurally-grown fixtures),
parallel TOML parsers that have already drifted, and an orphan
blessed-fleet directory that survived Phase 1.

---

## 1 — Bless v2 design review

**Contract correctness.** The shape (Identical → WithinTolerance →
Drifted; ZeroOutput / DtypeOrShapeMismatch / ReadFailed as orthogonal
failure modes) is the right enumeration. The `BlessVerdict` enum at
`tests/BlessCompare.h:44-51` cleanly separates "drift inside spec"
from "drift outside spec" from "structural" from "I/O" failure.
Diagnostic strings carry magnitudes (`tests/BlessCompare.cpp:319-336`)
which is the actual win over `BINARY DIFF: file (X bytes vs X bytes)`.

**ZeroOutput sanity at the right granularity.** The check at
`tests/BlessCompare.cpp:285-296` is asymmetric on purpose: only flag
if blessed had data and run lost it, not the other way around.
That's what makes `[arrays.dssp_chi] min_nonzero_fraction = 0.0`
legitimate — sparsity is per-array; the policy file is an honest
enumeration of "this output is sparse-by-construction." Per-array
policies in `bless_policy.toml` are documented inline with reasons,
which is exactly the maintainability shape.

**Default 0.05 defensibility.** For protein-scale (~hundreds to
thousands of atoms), 5% nonzero fraction means at least one in
twenty atoms must carry signal in any array supposed to be dense.
Looking at the override list: every per-array exception is a
sparse-by-construction case (DSSP one-hot, aromatic-only EFG slots,
mutation-delta arrays nonzero only at the mutated residue).
Default 0.05 is defensible for the existing dense calculator
outputs (BS, HM, McConnell, Coulomb shielding, AIMNet2 EFG full).
One concern: 0.05 is fine on 1UBQ-class (1231 atoms), but on a
100-atom test fixture (e.g. ALA capped tripeptide), 0.05 means 5
atoms — variance may be high enough that one of the override-needing
arrays trips. Worth thinking about if/when dense-validation moves
beyond the two large fixtures.

**Things missing or that will bite.**

- **`log.jsonl` is blessed but not compared.**
  `tests/golden/blessed/{nodft,withdft}/` contains `log.jsonl`
  (149 lines for nodft, similar for withdft), but `BinaryCompare`
  at `tests/test_smoke.cpp:316-389` only iterates `NpyFiles()`
  (extension filter at `:102`). The blessed log file is dead weight.
  Not a correctness problem, but if a future reader thinks log
  content is part of the contract, that's misleading. Either
  compare it (with sensible normalization for timestamps) or stop
  blessing it.
- **No "missing in run" failure mode.** At
  `tests/test_smoke.cpp:322-325`, missing-in-run fires
  `EXPECT_TRUE(run_files.count(f))` so it does fail the test.
  But "new file in run, not in blessed" only logs a `NEW file` line
  and doesn't fail (`:326-330`). For an extraction-pipeline
  regression that adds an unintended NPY, that's a silent
  forward-only addition. May be intentional (forward-compat); worth
  thinking about.
- **`FilesByteIdentical` short-circuit
  (`tests/BlessCompare.cpp:202-218`) reads both files twice for
  the byte-identical case** — once for tellg-comparison, once via
  memcmp loop. Tiny inefficiency.
- **The `default` policy section's per-array inheritance only works
  if `[default]` precedes any `[arrays.*]`**
  (`tests/BlessCompare.cpp:434` comment). This is a parser
  convention, not a parser invariant. A future re-ordering of the
  TOML would silently zero per-array policies; no test catches this.
- **`min_nonzero_fraction` default applied silently when no policy
  file.** If `bless_policy.toml` is absent (or path mistyped),
  `LoadTable` returns the empty-default table at
  `tests/BlessCompare.cpp:399-401` — every array gets the default
  policy and no warning. The "warn once per process" intent in the
  header (`tests/BlessCompare.h:69`) isn't actually implemented.
- **NPY reader supports four dtypes (`<f8`, `<f4`, `<i4`, `<i8`).**
  All current `NpyWriter` outputs hit those. If a calculator adds
  `<u4` (e.g. unsigned counts) the bless framework silently calls
  `ReadFailed`. Worth flagging — `BondLengthStats` `n_frames` is
  `uint32_t` per memory entry pattern, and trajectory H5 outputs
  can drift. Not a today-bug; a watch.
- **Parallel implementations now exist for the 0.05 default.**
  Hardcoded in two places — `tests/BlessCompare.h:41`
  (`BlessPolicy` struct default) and
  `tests/golden/blessed/bless_policy.toml:32`
  (`min_nonzero_fraction = 0.05`). If one is changed without the
  other, the `[default]` section in TOML overrides the C++ default
  for callers passing the toml path; callers passing empty path get
  the C++ default. Not load-bearing today (the toml path is always
  passed in `test_smoke.cpp:333-336`) but worth knowing.

**Beyond smoke.** The framework will work fine for any test that
needs to compare an NPY directory against a blessed one, but the
per-array policy is keyed on filename stem with no scoping to what
test produced it. If `tests/golden/blessed/nodft/coulomb_shielding.npy`
and a future
`tests/golden/blessed/some_other_test/coulomb_shielding.npy` need
different tolerances, the policy can't express that —
`arrays.coulomb_shielding` would apply globally. May be premature
concern but the design is single-namespace.

---

## 2 — AIMNet2 contract enforcement coverage

**Where AIMNet2 IS enforced (good):**

- Smoke `SmokeTest::SetUpTestSuite` (`tests/test_smoke.cpp:158-168`)
  hard-fails if the model is missing. Both `NoDft` and `WithDft`
  exercise the AIMNet2 path (`tests/test_smoke.cpp:422, 459`),
  confirmed in the most recent log files at
  `tests/golden/smoke/2026-05-04_175506/nodft/log.jsonl` and
  `2026-05-04_175613/withdft/log.jsonl` showing `AIMNet2Result wrote
  5 arrays`.
- Trajectory mode via `RunConfiguration::PerFrameExtractionSet`
  (`src/RunConfiguration.cpp:104, 122`) — `SetRequiresAimnet2(true)`
  plus `RequireConformationResult(typeid(AIMNet2Result))`, validated
  at `src/Trajectory.cpp:145-150` with `kConfigRequiresAimnet2`. So
  a `--trajectory` run with a Session that didn't `LoadAimnet2Model`
  hard-fails at Phase 4.
- Python loader at `python/nmr_extract/_protein.py:299-302` raises
  `FileNotFoundError` on any missing required NPY. With
  `_catalog.py:172-180` setting all five AIMNet2 arrays
  `required=True`, baseline/mutant fixtures lacking them properly
  skip via the fixture handlers' `try/except` at
  `tests/test_load.py:55-58, 64-68`.

**Where AIMNet2 is NOT enforced but the contract memory says it
should be:**

- **`--mutant` mode does not require `--aimnet2 MODEL` at JobSpec
  validation time.** `src/JobSpec.cpp:325-328` for
  `JobMode::Mutant` only validates the WT and ALA ORCA file pairs.
  The AIMNet2 path check at `:292-298` is "if specified must
  exist" — never required. The contract memory
  `project_aimnet2_contract_20260426` explicitly says: *"`--aimnet2
  MODEL` is now required for `--mutant` runs. Validate-time hard-fail
  with diagnostic. Mirrors the `_nmr.out` gate added earlier in the
  same session."* That gate was either never landed or was reverted
  with the IUPAC topology rollback. The TEST_HEALTH F6 names
  AIMNet2 specifically as the one that got partially fixed; this is
  the unfinished half.
- **`src/nmr_extract.cpp:339-347` silently skips AIMNet2 if neither
  CLI `--aimnet2` nor TOML `aimnet2_model_path` is set.** No
  warning, no error, no log. A user running `nmr_extract --mutant
  ...` without specifying the model produces an output that's
  missing the contract-required arrays, the C++ side never
  complains, and the Python SDK loader at
  `python/nmr_extract/_protein.py:299` raises `FileNotFoundError`
  on read-back. The error is far from where the actual mistake was.
- **`--orca` mode same shape.** Same JobSpec fall-through at
  `:322-323`; same silent skip in `nmr_extract.cpp` at `:339`. The
  contract memory is more lax on `--orca` ("Other modes accept
  `--aimnet2` but do not require it"), but the smoke `WithDft` test
  enforces it through the test fixture, suggesting the contract
  intent has tightened beyond the memory's scope.

**Pipeline-running tests that don't load AIMNet2 (per
`feedback_no_convenience_exclusions`):**

- **`tests/test_calculation_runner.cpp` lines 38, 72, 124, 129** —
  three tests build `RunOptions` without `aimnet2_model`. Pipelines
  run without AIMNet2 silently.
- **`tests/test_pipeline_and_sample.cpp` lines 69, 104, 184, 249,
  291, 326, 393, 411** — eight call sites; some pass empty
  `RunOptions{}` directly (e.g. `:184`).
- **`tests/test_write_features.cpp` line 38** — single test, builds
  `RunOptions` without AIMNet2.
- **`tests/test_job_spec.cpp` lines 256, 285, 312, 346, 352** — five
  `JobSpecE2E` end-to-end tests. `MutantEndToEnd` (`:327-365`) is
  the most relevant — it's exactly the use case the contract names,
  and the test runs the pipeline without AIMNet2.
- **`tests/test_amber_streaming.cpp` lines 138, 168, 225, 263, 303**
  — five `nmr::Session session;` constructions, none load AIMNet2.
  The `RunConfiguration` is custom-built per-test with
  `skip_mopac/coulomb/apbs/dssp = true`. None call
  `SetRequiresAimnet2(true)`. This is a deliberate convenience
  exclusion in violation of the no-convenience-exclusions feedback
  memory.

The two convenience-exclusion memories together
(`feedback_no_convenience_exclusions` at 9 days old,
`project_aimnet2_contract_20260426` at 8 days old) make the rule
concrete: AIMNet2 is baseline, only MOPAC is genuinely prohibitive,
tests must start from the production baseline. The `b9636f0`
checkin closed the smoke gap but left at least 18+ test sites
running pipelines without AIMNet2.

---

## 3 — Test fixture discipline

**Stale fixtures relative to the AIMNet2 `required=True` contract:**

- **`baseline_features/P84477/`** (57 files) — flagged in TEST_HEALTH
  F5 as deferred. No AIMNet2 NPYs (confirmed by listing). Consumed
  by `test_load.py::TestMopac.test_*` which now skip cleanly via
  fixture handler. Status correctly documented.
- **`/tmp/sdk_mutant_test`** — flagged in F5; not on disk; pytest
  skips. Status correctly documented.
- **`tests/data/sdk_geo_only/1Q8K/1Q8K_10023_02..06_*` (frames
  002-006)** — only 32 NPYs each, no AIMNet2 NPYs. Only frame 001
  was regenerated for the contract. The `test_load.py::geo` fixture
  loads only frame 001 (`GEO_ONLY` constant at `:36`). So the other
  five frames are accumulated cruft that no test currently consumes
  — they can be deleted, OR the geo fixture should load all six and
  validate cross-frame consistency. Currently they sit on disk as
  stale orphans.
- **`tests/data/orca/A0A7C5FAR6_*`** — used by JobSpec E2E and many
  calculator tests. PDBs, prmtops, XYZ, NMR `.out` from before
  AIMNet2 contract; never extracted via the pipeline (raw fixture,
  not blessed output) so AIMNet2 doesn't apply. OK as-is.
- **`tests/data/external/1OKH_4587_protonated.pdb`** — flagged "FF
  unclear" in `spec/plan/test_inventory_2026-05-03.md` line 22.
  Consumed by `test_protonation_detection.cpp` and
  `test_foundation_results.cpp`. Provenance is marked TBD; nobody
  re-verified.
- **`tests/data/external/1UBQ.pdb`** — raw crystal PDB, no FF
  context needed. OK.

**Other fixture observations:**

- **`tests/golden/blessed/fleet/frame_001/`** survived Phase-1
  cleanup as a 53-file orphan. The CHARMM `tests/data/fleet/` was
  retired in `b9636f0`, but `tests/golden/blessed/fleet/` remains.
  It's git-tracked (151 files in `tests/golden/blessed/` total per
  `git ls-files`) despite the `.gitignore` line
  `tests/golden/blessed/`. The gitignore doesn't retroactively
  untrack. No live test consumes this dir; it's pure cruft.
- **`tests/data/fleet_amber/_backup_round1_*/`** and
  **`_backup_round2_*/`** — gitignored backup MD trees from prior
  round runs. Not a fixture problem per se but they take up disk
  and aren't named in any spec doc.
- **`tests/regression/baseline_orca/` and `output_fleet/`** — the
  regression script at `tests/regression/run_regression.sh` has its
  own baseline + a leftover `output_fleet/` dir from before the
  fleet retirement. The script is not invoked by CMake/CTest (no
  `add_test` reference), so it's effectively orphan. Use Case D was
  excerpted to bones; Use Case B (lines 38-85) appears to call
  `nmr_extract --orca "$ORCA_DIR"` (`:46`) without `--root NAME`,
  which the current JobSpec rejects (`src/JobSpec.cpp:96-100`
  requires `--root`). The regression script is broken as well as
  orphan.

---

## 4 — Contract-without-enforcement audit (F6)

Walking the memory entries dated 2026-04-15 onward with code-side
claims:

**`project_aimnet2_contract_20260426`** — three layers claimed:
- Catalog `required=True`: **enforced**
  (`python/nmr_extract/_catalog.py:172-180`).
- SDK loader `Protein.aimnet2: AIMNet2Group`: **partially**.
  `python/nmr_extract/_protein.py:404-414` still does
  `aimnet2 = None; if "aimnet2_charges" in available: aimnet2 =
  AIMNet2Group(...)`. The `None` default means a Protein constructed
  from a non-required path retains the Optional shape. The hard-fail
  comes from the catalog `required=True` raising before construction,
  not from the dataclass annotation. Functionally it's enforced;
  structurally the field still allows None.
- JobSpec mutant mode `--aimnet2 MODEL` required: **NOT enforced**.
  `src/JobSpec.cpp:325-328` for `JobMode::Mutant` doesn't check.
  This is the contract-not-landed instance most closely paralleling
  the catalog gap that F1 just fixed.

**`project_iupac_topology_landed_20260426`** — four artifacts
claimed (commits `fc76f47 / 145d2cc / 6448272`):
- AtomTopology + IupacAtomName + AtomReference: **REVERTED**.
  Confirmed: no `AtomReference.h`, no `AtomTopology.h`, no
  `IupacAtomName.h` in `src/`. The 2026-04-27 revert dropped all of
  it. This memory entry is stale-by-revert; CLAUDE.md says don't
  even investigate. The "three debts" listed in this memory
  (pro-R/pro-S verification, variant overrides, DFT-compare
  completeness) are also moot at the implementation level because
  the substrate they targeted no longer exists.

**`project_dft_compare_calculator_completeness`** — claimed
`MutationDeltaResult` carries only `delta_shielding`; should add
WT total + mutant total + dia + para for both:
- The file `src/MutationDeltaResult.{h,cpp}` exists. The
  completeness claim hasn't been audited in this review (would
  require reading those files); the claim is "must add", not "is
  enforced." The mutant validation gate
  (`project_mutant_validation_gate_20260426`) — two FULL 723-pair
  runs comparing outputs — was never run. So this is an unfinished
  thread, not a contract violation per se.

**`project_charmm_retired_amber_only_2026-05-02`** —
quarantined-legacy code listed:
- `src/xtc_reader.h`: **gone from src/, in bones/** (verified
  `tests/bones/test_water_field.cpp` includes it).
- `src/GromacsEnsembleLoader.{h,cpp}`: **gone from src/, in bones/**.
- `src/ChargeSource.cpp::GmxTprChargeSource::LoadCharges` subprocess
  path: **gone, comment at `src/ChargeSource.cpp:262-263` confirms
  retirement to bones**. Compliant.

**`project_catch_all_frame_pull_cr`** — claimed
`GromacsFramePullResult` lives at `src/GromacsFramePullResult.{h,cpp}`
with `velocities`/`box_matrix`/per-bond-lengths/etc., attached when
`RunOptions::velocities` or `box_matrix` non-null:
- `src/GromacsFramePullResult.{h,cpp}` exists.
  `src/OperationRunner.cpp:197-200` attaches conditionally. The
  "all derived geometry" claim (per-bond/per-angle/per-dihedral)
  wasn't audited at the field level in this review.

**`feedback_capture_at_the_boundary`** +
**`feedback_no_attach_lifecycle_for_invariant_data`** +
**`feedback_readback_block_is_a_compiler_trace`** — three companion
architectural rules:
- `LegacyAmberTopology` plain fields, no Attach/Has/Optional:
  confirmed by `grep -l "MutableLegacyAmber\|AttachAmberFFData\|
  HasAmberFFData\|AmberFFData"` returning empty (the May 2
  wrapper-revert landed).
- `GromacsToAmberReadbackBlock` exists at
  `src/GromacsToAmberReadbackBlock.{h,cpp}`. Compliant.

**`feedback_no_convenience_exclusions`** — claimed only MOPAC may
be excluded from test pipelines; APBS + AIMNet2 are baseline:
- **NOT enforced**. The 18+ test sites listed under question 2
  violate this rule. The rule is an explicit feedback memory dated
  9 days ago. The rule applies "When building a test
  `RunConfiguration` or `RunOptions`", and dozens of tests do so
  without AIMNet2.

**`feedback_test_invocation_via_ctest`** +
**`reference_nvrtc_rpath_fix`** — claimed: tests must run via ctest
(CMake attaches `LD_LIBRARY_PATH` per-test), never
`./build/<binary>`. The pattern is enforced via CMake
`gtest_discover_tests ... PROPERTIES ENVIRONMENT`. Not auditable
from source alone, but TEST_HEALTH §How-to-refresh uses the correct
`ctest --test-dir build` invocation.

**`feedback_no_pdbfixer_mid_pipeline`** — no concrete code claim;
rule for what NOT to do. Auditable only by absence: `grep -r
pdbfixer src/` — not done in this review; flagged for any future
audit.

**`feedback_no_regex_structural_biology`** — rule, not concrete
claim. Auditable by spot-checking parser sites; out of scope here.

**`project_findorcafiles_bug`** (28 days old) — claimed 9 sites
doing alphabetical glob:
- `src/nmr_extract.cpp:43 FindOrcaFiles()` — **dissolved**. The
  current `JobSpec` uses `ExpandOrcaRoot()` at
  `src/JobSpec.cpp:38-44` with `--root NAME` deterministically
  expanding to `{root}.xyz`, `{root}.prmtop`, `{root}_nmr.out`.
  Compliant.
- The 6 batch-test `FindNmrOutput()` sites: not audited in this
  review. Worth a follow-up grep.

**`project_no_symlinks_rule`** — rule; auditable via `find
/shared/2026Thesis/nmr-shielding -type l`. Not done.

**`project_aimnet2_contract_20260426`** sub-claim about regenerated
`tests/data/sdk_geo_only/1Q8K/.../`:
- Frame 001 regenerated (44 files including AIMNet2). Frames
  002-006 are 32 files each, NOT regenerated. The contract memory
  said "regenerated with AIMNet2 included (49 NPYs +
  atom_name.csv)" — current frame 001 has 43 npy + atom_name.csv =
  44 entries. Either the count drifted (some npy were removed since
  2026-04-26) or the memory's count was off. Either way, frames
  002-006 are stale.

**Summary of contract-vs-code drift findings:**

| Memory entry | Code-side claim | Status |
|---|---|---|
| `project_aimnet2_contract_20260426` | catalog `required=True` | **enforced** |
| `project_aimnet2_contract_20260426` | SDK `Protein.aimnet2` non-Optional | **partially** (still `None` default; hard-fail comes from catalog) |
| `project_aimnet2_contract_20260426` | JobSpec mutant requires `--aimnet2` | **NOT enforced** |
| `project_iupac_topology_landed_20260426` | 3 commits with new types | **REVERTED** — entry is stale by revert |
| `project_dft_compare_calculator_completeness` | MutationDeltaResult must add 6 tensors | **not audited at field level** |
| `project_charmm_retired_amber_only_2026-05-02` | three quarantined-legacy modules retired | **compliant** |
| `feedback_capture_at_the_boundary` + `feedback_no_attach_lifecycle...` | LegacyAmberTopology plain fields | **compliant** |
| `feedback_readback_block_is_a_compiler_trace` | GromacsToAmberReadbackBlock exists | **compliant** |
| `feedback_no_convenience_exclusions` | MOPAC-only exclusions in tests | **NOT enforced** (18+ sites) |
| `project_findorcafiles_bug` Pri0 | FindOrcaFiles must die | **compliant** |
| `project_findorcafiles_bug` Med | 6 batch-test glob sites | **not audited** |
| `project_aimnet2_contract_20260426` sub | sdk_geo_only/1Q8K regenerated | **frame 001 only; frames 002-006 stale** |

The dominant pattern is "memory says X is enforced, code mostly
enforces X but one path was missed." For AIMNet2 specifically the
JobSpec gate is the live miss.

---

## 5 — Anything else weird, surprising, or worth tightening

**Three parallel TOML parsers, all hand-rolled, all subtly
different:**
- `src/CalculatorConfig.cpp:137-174` — supports comment stripping,
  `=` split, double + string values; no sections, no nested tables.
  Will silently mis-parse a `[default]` header (treats `[` as no
  `=`, skips the line).
- `tests/TestEnvironment.cpp:50-90` — same shape, double-quote
  strip on values, no sections.
- `tests/BlessCompare.cpp:389-454` — supports sections (`[default]`
  and `[arrays.<stem>]`), inheritance from default to per-array,
  scientific notation. Most complete.

The three implementations have already drifted: BlessCompare has
section support; the other two would silently mishandle the same
TOML file. If a future test reads `bless_policy.toml` via
`TestEnvironment` it'd parse `[default]` as a no-op line and the
`rtol` line as an entry called `rtol`. If `testpaths.toml` ever
gets a section header, `TestEnvironment` would silently skip the
entry above it. Worth one shared minimal TOML parser; toml++ is
header-only and small.

**Smoke test SetUpTestSuite skip behavior.**
`tests/test_smoke.cpp:158-168` uses
`ASSERT_FALSE(model_path.empty())` and `ASSERT_TRUE(fs::exists(...))`
— `ASSERT_*` in `SetUpTestSuite` aborts the entire suite if either
fails. That's correct policy (the contract says AIMNet2 is
required), but it means a machine without the model file gets a
hard fail rather than a clean skip. May be intentional; just note
that "missing AIMNet2 model on this machine" no longer surfaces as
a per-test SKIP — it crashes the suite. Compare to the GTEST_SKIP
pattern at `:404, 443`.

**AmberStreaming tests as a baseline mismatch.** Per
`feedback_no_convenience_exclusions` and
`project_aimnet2_contract_20260426`, the trajectory baseline IS
AIMNet2-required (`PerFrameExtractionSet` sets it). But
`tests/test_amber_streaming.cpp` builds custom `RunConfiguration`
objects from scratch with `skip_mopac = skip_coulomb = skip_apbs =
skip_dssp = true` and never calls `SetRequiresAimnet2(true)`. The
intent is "test the trajectory mechanics, not the calculators" —
but per the no-convenience-exclusions rule, the only legitimate
exclusion is MOPAC. APBS + AIMNet2 should run. The four discipline
tests (BondLengthStatsFrame0Semantics, FinalizeIdempotency,
H5RoundTrip, end-to-end) all skip-stack the calculators.

**JobSpecE2E tests.** Five tests run real builders + OperationRunner
on real fixtures (`PdbEndToEnd`, `ProtonatedPdbEndToEnd`,
`OrcaEndToEnd`, `MutantEndToEnd`). None load AIMNet2.
`MutantEndToEnd` is precisely the case the contract memory targets.
Same baseline-mismatch pattern.

**`aimnet2_charge_sensitivity` legacy entry.** TEST_HEALTH F4 marks
this `INFORMATIONAL`. The `python/tests/test_load.py:483-487` test
`test_aimnet2_specs_registered` explicitly asserts the legacy stem
is in CATALOG. So removing it from the catalog would break a test,
even though no current code produces the file. The legacy is now
contractually preserved by a test assertion. Not necessarily wrong,
but worth noting the test enforces the legacy.

**Pre-existing BlessNotes vs reality.**
`tests/golden/blessed/BLESS_NOTES.md:42-47` says "nodft/ went from
53 → 58 NPYs" — actual is 57 files (58 minus log.jsonl when
filtered). `withdft/` "from 56 → 61" — actual is 61 files (60 NPYs
+ log.jsonl). Counts in BLESS_NOTES are off by the log.jsonl
bookkeeping; minor.

**Smoke test `min_npy_files` thresholds.**
`tests/test_smoke.cpp:433` (`RunSmoke("nodft", conf, opts, 14, 40,
blessed)`) requires ≥40 NPYs; actual blessed is 57 NPYs.
`tests/test_smoke.cpp:478` (`RunSmoke("withdft", ..., 15, 45, ...)`)
requires ≥45; actual is 60. The thresholds are liberal cushions
but they're far below current production counts; a calculator
silently dropping 5-10 NPYs would still pass the count check and
need to be caught by per-file structural check.

**`tests/golden/blessed/` .gitignore + 151 tracked files.** The
`.gitignore` line `tests/golden/blessed/` was added at some point
but doesn't retroactively untrack the 151 historical files
(including the entire `fleet/frame_001/` orphan). Updating the
gitignore doesn't clean up — `git rm --cached` would. The current
state is hybrid: some blessed dirs gitignored, some tracked, no
clear intent visible from inspection.

**`tests/regression/run_regression.sh` is broken AND orphan.** Not
invoked by CMake/CTest. Use Case B at `:46` calls `--orca
"$ORCA_DIR"` without `--root NAME` (current JobSpec rejects). Use
Case D excerpted to bones. The script's `output_fleet/` artifact
dir from a prior fleet run is still there. Either fix-and-wire-in
or retire entirely.

**`gmx_protonated` / `1OKH_4587` provenance.**
`spec/plan/test_inventory_2026-05-03.md:22` flagged "FF unclear;
needs check" and "TBD." Eight days later, status is still TBD. Two
tests consume it (`test_protonation_detection.cpp`,
`test_foundation_results.cpp`). If it's CHARMM-era like the rest of
the fleet retirement, the same dead-letter argument applies.

**Skip count consistency.** TEST_HEALTH §C++ ctest detail lists 11
skips; the breakdown (`AmberPreparedChargeStep5NegativeTest.LoadCharges*`
× 2, `PropkaTest.*` × 6, `KamlTest.*` × 3) sums to 11. The Propka
skip count of 6 is concerning — that's 6 separate test cases all
skipped on the same condition. If one Propka assertion suite covers
what's needed, six test cases asserting essentially the same skip
is duplication. Worth checking whether they share a fixture-class
condition that should hoist to suite-level.

**`Trajectory.cpp:239` uses `md.edr` but JobSpec sets
`production.edr`.** `src/nmr_extract.cpp:239` constructs `Trajectory
traj(spec.traj_xtc, spec.traj_tpr, fs::path(spec.traj_dir) /
"md.edr")` — hardcoded "md.edr". But `src/JobSpec.cpp:161` sets
`spec.traj_edr = dir + "/production.edr"`. So the EDR file path in
spec is set but never used; the EDR path passed to Trajectory is
recomputed inline as `md.edr`. Either intentional fallback or a bug.
The spec field is accessed by `ValidateJobSpec` at `:336-337` for
existence-check — so we validate `production.edr` exists, then
construct the Trajectory pointing at `md.edr` (which may not exist,
may be a different file).

---

## What I would do next (suggestions, ranked)

1. **Audit the JobSpec mutant `--aimnet2` gate.** Either land it
   (matches the contract memory) or revise the memory and the smoke
   wire-in to match the actual rule. The current state — smoke
   enforces it, JobSpec lets it slide silently — is the worst of
   both worlds because it lets users produce contract-violating
   output and then be told by the SDK loader after the fact.

2. **Decide what "production baseline" means for non-smoke tests,
   then enforce it.** Per `feedback_no_convenience_exclusions` plus
   the AIMNet2 contract, the rule is "MOPAC-only exclusions." The
   18+ sites currently violating that rule either need updating
   (AIMNet2 model loaded via process Session, like smoke does) or
   the rule needs revising. Pick one.

3. **Settle the EDR path inconsistency in `nmr_extract.cpp:239` vs
   `JobSpec.cpp:161`.** Five-minute investigation; either there's a
   rationale for re-deriving the path OR there's a bug where
   validation says "production.edr exists" and then we open
   "md.edr."

4. **Retire `tests/golden/blessed/fleet/` and the
   `tests/regression/output_fleet/` cruft.** The Phase-1 cleanup
   intentionally retired the consumers; the orphan artifacts can
   follow. While at it, make a call on
   `tests/regression/run_regression.sh` (fix and wire in vs delete).

5. **Audit `gmx_protonated` / `1OKH_4587_protonated.pdb`
   provenance.** Two tests consume it; status unclear since
   2026-05-03. Either confirm AMBER provenance, mark it for
   retirement, or document the FF intent.

6. **Decide what to do about frames 002-006 of `sdk_geo_only/1Q8K`.**
   They're stale relative to the AIMNet2 contract and aren't
   consumed by any test. Either delete or extend the SDK fixture to
   load and validate cross-frame.

7. **Consolidate the three TOML parsers.** Lowest priority because
   they "work" today, but the drift is real and a single shared
   minimal parser (or toml++) would prevent the next surprise.

8. **Tighten the `smoke_tests` `min_npy_files` thresholds.** Setting
   them to ~95% of expected (e.g., 55 for nodft, 58 for withdft)
   catches partial-output regressions early. Current 40/45 are far
   enough below 57/60 that real regressions can pass.

9. **Address the F6 audit pattern systemically.** The four-line
   audit table in §4 above is the shape of one finding. The same
   pattern (memory entry asserts code-side enforcement that wasn't
   fully landed) likely recurs as the project accumulates contract
   memories. Worth a "verify before write" rule for new contract
   memories, or a periodic sweep tied to doc-cleanup.

10. **Skim `MutationDeltaResult` for the dia/para completeness debt.**
    `project_dft_compare_calculator_completeness` flagged 6 tensors
    + 3 deltas to add. Untouched in the audit. This is a thesis-
    relevant correctness gap, not just a test concern.
