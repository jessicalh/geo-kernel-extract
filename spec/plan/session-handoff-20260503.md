# Session handoff — 2026-05-03

**Pickup brief.** Push of 11 commits to origin landed (cleanup + AMBER
readback + gap fix). Body of the session was triage of the topology
work via a critical agent + conversation. Outcome: a concrete cleanup
plan to retire dead-letter scope (CHARMM-era tests + fixtures + loader
code, fes-sampler entirely), revive the AMBER counterparts of what the
dead-letter tests covered (loader assertions on the two cases, solvent
calculator validation), pick up `CalculatorContract` adoption as a
typedef-only convention along the way, and re-bless the smoke tests.
Then forward to IUPAC projection output and additional topology
internals once the substrate is solid.

This handoff is the work plan for the next pass. Reading order at the
bottom.

---

## State at end of session

**HEAD:** `776ca75` on `master`. `origin/master` matches (push +
filter-repo + gap-fix + README-restore landed earlier in this session;
11 commits pushed total).

**Working tree:** three plan docs written this session, all untracked:

- `spec/plan/review_items_to_assess.md` — triage doc from agent review +
  this session's conversation (A/B/C structure + Triage Outcomes
  section at end)
- `spec/plan/test_inventory_2026-05-03.md` — per-cpp disposition of all
  51 test files + fixture taxonomy
- `spec/plan/session-handoff-20260503.md` — this file

Plus `spec/plan/session-handoff-20260502-evening.md` still untracked
from yesterday. `Testing/` is build artifact (separate cleanup).

**Backups present:**

- `/mnt/expansion/nmr-shielding-backup-pre-push-20260503/` — full
  validated 97 GB copy of nmr-shielding before the push. Defensive,
  can be deleted whenever.
- `data/ccd/components.cif` + `.gz` are gitignored on disk (~465 MB);
  re-fetchable from the URL in `data/ccd/README.md`.

---

## What was decided this session

### D1. `CalculatorContract` = typedef-as-documentation

CHARMM is dead-letter. The current prep+run pipeline funnels everything
into AMBER ff14SB through GROMACS readback. With one topology kind in
the live tree, the typed enforcement value of `CalculatorContract`'s
`static_assert` collapses (it's checking a tautology). Decision: the
contract becomes a pure typedef each calculator declares at its class
header. The static_assert and accessor templates come out of
`src/CalculatorContract.h`. The discipline is a convention in
`PATTERNS.md`, not a type-system gate.

The typedef gets added to each calculator at the natural moment its
test is touched in Phase 2 / 3 below. No separate "Contract pass"
across all calculators in a single sitting. (Captured as T1 in
`review_items_to_assess.md`.)

### D2. fes-sampler is retired

We don't do fes anymore. Test, fixture, executable, memory entry, and
any stragglers come out.

### D3. `fleet_amber/` trajectories are full-system

Verified mid-session: `production.trr` is 1.8 GB with positions +
velocities, `production.xtc` is 546 MB with `compressed-x-grps =
System` (`production.mdp:13`). `GromacsFrameHandler.cpp:65,178-179`
reads the full system and splits via `SysReader().ExtractFrame(...)`
into protein slice + solvent slice every frame. `traj.env_.solvent`
carries the per-frame solvent data. The water/hydration/GROMACS-energy
calculators can consume it directly on the AMBER path. **No new
fixture preparation needed** — the rewrite of `test_water_field.cpp`
goes against existing `fleet_amber/`.

### D4. The two test cases as scientific anchors

`fleet_amber/1P9J_5801` and `fleet_amber/1Z9B_6577` are deliberately
distinct in chemistry:

- **1P9J_5801**: 3 disulfides (CYX), HISH→HIP residues. The case the
  readback-block HISH unblock landed against. Has published
  S²/T1/T2 NMR relaxation data.
- **1Z9B_6577**: no CYS at all. The case the
  `has_disulfide_authority` gap fix (2026-05-03) was written against —
  authority-says-zero must demote any geometric Disulfide tag. Also
  has published relaxation data.

A loader/fleet test on these two carries the chemistry-authority
chain end-to-end: 1P9J validates the populated path, 1Z9B validates
the empty-authority path. Different from `test_amber_trajectory.cpp`,
which validates the topology-build itself.

### D5. 28 deferred NPY bless is not OK

The 2026-04-24 deferral has accumulated past its purpose. AMBER
readback work has settled. The smoke tests will be re-blessed in
Phase 4 below — either against current code with documented drift, or
the bless-per-commit framing gets retired in favor of milestone
blessing.

---

## The plan, in order

### Phase 1 — Move dead-letter to bones/

Goal: the live test/fixture/source surface contains only what actually
runs against the AMBER substrate or against fixture-free unit logic.

**Test sources to retire:**

- `tests/test_fleet_loader.cpp` (4 of 5 inert FleetLoader failures)
- `tests/test_smoke_fleet.cpp` (the 5th)
- `tests/test_smoke_fes_fleet.cpp` (fes retired)
- `tests/test_water_field.cpp` (current form; reborn against
  `fleet_amber/` in Phase 2)
- The fleet-using subset of `tests/test_job_spec.cpp` (split or
  rewrite; the rest stays)

**Fixtures to retire:**

- `tests/data/fleet/` (CHARMM/XTC: 1A6J_5789 + 1AEP_4814)
- `tests/data/fes_fleet/` (1HD6_4820 + 1I8X_4351)
- `tests/data/fleet_test_fullsys/` (1ZR7_6721 CHARMM full-system)
- `tests/data/fleet_test_large/` (PLUMED ensemble; orphan, no test
  references)

**Quarantined-legacy code to retire:**

- `src/xtc_reader.h`
- `src/GromacsEnsembleLoader.{h,cpp}`
- The `GmxTprChargeSource::LoadCharges` subprocess path in
  `src/ChargeSource.cpp` (the typed loader path stays; only the
  subprocess call goes)

**CMakeLists.txt adjustments:**

- Drop the `fes_fleet_smoke_tests`, `fleet_smoke_tests`, and
  `water_field_tests` executable blocks (water_field comes back in
  Phase 2 against AMBER).
- `trajectory_tests` keeps only `test_amber_trajectory.cpp` (and
  whatever Phase 2 produces); rename if cleaner.
- Sweep `target_link_libraries` and `target_include_directories` for
  references to anything moved.

**Smoke run history:**

- Prune pre-2026-04-30 dated dirs in `tests/golden/smoke/` (28 of
  them; they pre-date current scope and the bless rework in Phase 4
  will produce fresh ones).

**Memory:**

- Retire `project_fes_sampler_state.md` from MEMORY.md and
  the memory file itself.

**Sweep before declaring Phase 1 done:**

```
grep -rn 'fes_fleet\|fes_sampler\|FES_' src/ tests/ python/ scripts/
grep -rn 'xtc_reader\|GromacsEnsembleLoader\|GmxTprChargeSource' src/ tests/ python/
grep -rn 'fleet_test_fullsys\|fleet_test_large\|fleet_data' src/ tests/ python/
```

Catch any remaining live consumers before deleting.

**Open question:** does "bones" mean `learn/bones/` (existing,
historical session notes) or `tests/bones/` (existing, has
`test-parameters-initial/` etc.) or a new `bones/` at repo root? Pick
one before moving.

**Open question:** move-and-keep-history (`git mv`) versus delete-and-
tombstone? Move preserves blame; delete is cleaner. Probably move,
since `bones/` connotes "kept for reference, not active."

### Phase 2 — Re-derive against AMBER fixtures

**New loader/fleet test on the two AMBER cases:**

Either a new `tests/test_amber_fleet_loader.cpp` or extension of
`tests/test_amber_trajectory.cpp`. Decide based on file size + topical
fit. Asserts on both 1P9J_5801 and 1Z9B_6577:

- TRR frame count vs expected (production length × cadence;
  per `testpaths.toml` Round-3 Option B 15 ns at `nstxout = 5000`).
- EDR row count + time alignment with TRR.
- Velocity presence: TRR has velocities, NaN-sentinel survives
  read-back, XTC stripping the velocity slot would fail (regression
  for the 2026-05-02 NaN-sentinel detection fix in
  `GromacsFrameHandler.cpp`).
- HISH residues in 1P9J resolve to HIP variants via readback (the
  readback-block unblock case).
- Disulfide bond count: 3 in 1P9J with the authority list applied;
  0 in 1Z9B with `has_disulfide_authority = true` and the empty list
  demoting any geometric tag (the gap-fix regression case).
- `PerFrameExtractionSet` runs end-to-end on a small frame stride
  (smoke-style, not the full 17-minute production sweep).

**Solvent calculator validation on AMBER:**

New `tests/test_water_field.cpp` (replaces the dead-letter version).
Uses `Trajectory::Run` + `GromacsFrameHandler` on `fleet_amber/
1P9J_5801` (or both fixtures). Drops `xtc_reader.h` and
`GromacsEnsembleLoader` includes entirely. Exercises:

- `WaterFieldResult` — sanity-check water field at backbone amide H
  positions, water count near protein scales sensibly.
- `HydrationShellResult` — first-shell water count > 0,
  consistent with protein surface area.
- `GromacsEnergyResult` — per-atom energy decomposition sums
  consistent with EDR aggregate.

This is the "we need to know AMBER works for solvent" verification.
The three calculators get their `Contract` typedef during the rewrite
(D1 / Phase 3).

**JobSpec test split:**

`tests/test_job_spec.cpp` mixes live-AMBER and fleet-using assertions.
The fleet-using subset (currently `JobSpecE2E.FleetLibraryDirect`) goes
to bones with the rest of the fleet test material, or gets re-derived
against `fleet_amber/`. The non-fleet assertions stay live.

### Phase 3 — Calculator `Contract` pass

Driven by Phase 2's natural touchpoints, then completed:

- The three solvent calculators get `using Contract = ...` during
  Phase 2's water-field rewrite.
- Other calculators get the typedef as the smoke test and structure
  tests are revisited in Phase 4.
- `src/CalculatorContract.h` shrinks to the bare struct: drop
  `static_assert`, drop `RequiredTopology<>` and
  `RequiredChargeTable<>` accessor templates. Comment becomes
  "convention typedef, type system does no enforcement work."
- `PATTERNS.md` gets a new pattern (or extension of existing): every
  calculator declares a `Contract` typedef. The typedef is the
  documentation of substrate dependence; reviewers check that a
  calculator's accesses match its declared contract.
- `INDEX.md` "When implementing a calculator" section gets an
  exemplar pointer naming the first calculator to adopt the typedef
  as the template for new work.

The 22 existing calculators get the typedef one-by-one, paired with
whatever test or feature work happens to touch them. No "all 22 in
one sitting" pass. (D1 captured this — the typedef rides along, it
doesn't drive its own phase.)

### Phase 4 — Re-bless smoke

The smoke test framework's blessed baselines are at
`tests/golden/blessed/withdft/`; 28 NPY binary-diff failures are
deferred per `project_smoke_test_bless_deferred_20260424` memory.

Steps:

- Diagnose: which 28 NPYs differ, by what magnitude, against which
  current-code behavior? Spot-check a representative few to
  understand whether drift is expected (e.g. the AMBER readback
  changing chemistry decisions, the disulfide authority changing
  which bonds get tagged) or unexpected.
- Re-bless against current code on `1ubq_protonated.pdb` +
  `consolidated/`. Document the drift in a methods-text-suitable
  comment alongside the bless.
- Decide on the bless rhythm going forward: per-commit binary-diff
  bless, or milestone bless? The deferral has effectively made the
  decision (~10 days without per-commit bless and the workflow
  hasn't broken). Codify whatever ends up working.

If the drift turns out to be unexplainable (a calculator producing
genuinely wrong output, not just different output), that's the actual
gate — bless waits until the bug is fixed. The AMBER readback +
disulfide changes are the most likely sources of expected drift; an
investigation pass should be able to attribute most diffs to those.

### Phase 5 — Forward (out of scope for this plan)

After Phases 1-4, the live surface is clean and the AMBER substrate
is verified. Then:

- **IUPAC projection output.** The Crystal Projection Rule
  (`amber-implementation-plan-2026-04-29.md:1383-1424`): names are
  pure functions over the typed substrate, never cached strings.
  IUPAC and BMRB projections become emission-side functions on
  `LegacyAmberTopology`. H5 schema gets the projected name columns
  as separate datasets derived from typed enums.

- **`LegacyAmberTopology` semantic fields.** Decision needed first
  (review_items_to_assess.md A2): is the OBJECT_MODEL.md "Pending
  full schema" block (Locant / BranchIndex / DiastereotopicIndex /
  ProchiralStereo / PlanarStereo / PseudoatomClass / PolarHKind /
  RingPosition) the live forward intent, or IUPAC scar tissue from
  the reverted 2026-04-26 attempt? Phase 5 either implements those
  fields or removes them from the doc.

- **Additional topology internals.** Construction-time work that
  enriches `LegacyAmberTopology` with whatever the projection layer
  needs but doesn't yet have. Driven by what calculators end up
  reaching for after Phase 3 — the substrate grows in response to
  actual demand, not speculation.

---

## Open questions for next session

1. **bones location.** `learn/bones/`, `tests/bones/`, or a new
   `bones/` at repo root? Pick one before Phase 1 starts.
2. **Move vs delete.** `git mv` to bones (preserves blame) versus
   delete-with-tombstone (cleaner)?
3. **Loader test file shape.** New `test_amber_fleet_loader.cpp` or
   extend `test_amber_trajectory.cpp`? Look at file sizes + topical
   fit.
4. **Bless investigation depth.** Spot-check a few drifted NPYs and
   re-bless wholesale, or full audit of all 28?
5. **A2 decision** (review_items_to_assess.md): IUPAC-vocabulary
   "Pending full schema" block — current intent, or scar tissue?
   Gates Phase 5.
6. **OpenAI's parallel review** — when its findings arrive, append to
   `review_items_to_assess.md` per the doc's own instructions.

---

## Reading order for next session

1. **This doc** (`spec/plan/session-handoff-20260503.md`) — start here.
2. `spec/plan/review_items_to_assess.md` — triage findings + Triage
   Outcomes section (D1, T2 are captured there; D2-D5 are this
   session's).
3. `spec/plan/test_inventory_2026-05-03.md` — per-cpp test
   disposition.
4. `spec/plan/session-handoff-20260502-evening.md` — yesterday's
   handoff. AMBER readback context, what shipped 2026-05-02.

After those, the standard reading order from `spec/INDEX.md` Tier 1
+ Tier 2 covers the substrate as it stands. The trio
`LegacyAmberTopology.h`, `GromacsToAmberReadbackBlock.h`,
`GromacsFramePullResult.h` is the architectural artifact produced
this week.

---

## Functionality to recover

From dead-letter, will be reborn in Phase 2:

- Fleet loader assertions (frame counts, EDR alignment, velocity
  presence, chemistry authority)
- Water field, hydration shell, GROMACS energy calculator validation
- Per-frame extraction smoke at the loader level

Per the 2026-04-13 OBJECT_MODEL.md note, `WaterFieldResult`,
`HydrationShellResult`, `GromacsEnergyResult`, plus `SasaResult`,
`AIMNet2Result`, the extended `DsspResult`, and `BondedEnergyResult`
were already documented as live; the dead-letter test infrastructure
was the only thing exercising the explicit-solvent ones. Restoring
the test path restores the validation guarantee.

The AMBER substrate gives us this for free — we built it; the work is
just rewiring the tests to use it.
