# TENTATIVE_OUTSTANDING_ISSUES.md

**Status:** tentative. The items below were surfaced by a four-agent
harvest pass over ~30 docs being retired to `spec/plan/bones/`
(2026-05-13). Each was nominally verified against source — agents
made 60-90 tool calls each, citing `file:line` evidence — but **the
verification was not re-audited**. Treat every entry as "agent says
this is true; spot-check before acting."

**Framing:** most items here are **payable debt**, not load-bearing
flags. The bugs section (2 items, OI-001 / OI-002) is worth
attention; the rest are routine items that get paid down as work
happens. Don't read every item every session — read them when
starting work that intersects.

Subsumes the former `spec/plan/bones/KNOWN_BUGS.md`, `spec/plan/bones/FIX_TESTS.md`, and
`spec/plan/bones/pending_decisions_20260423.md`. When confirming or paying down
an item, move it to **Resolved** with the resolving commit. Append
new items at the bottom of the matching section.

Entry conventions:
- **OI-NNN** stable ID.
- **Verified (per agent):** `file:line` or grep evidence per the
  harvest pass.
- **Source:** original doc where the item surfaced (most are bones'd).
- **Action:** one-line if obvious; otherwise the item sits.

---

## Bugs (broken behaviour, reproducible)

### OI-001 — ORCA shielding-load failure is silent-Ok

- **Verified:** `src/OperationRunner.cpp:276-284` — failure logs
  `OperationLog::Error` but does NOT set `out.error` and does NOT
  return; the `RunResult` returns Ok-shaped with no shielding result.
  A typo'd `--orca` path silently produces a WT-only output that
  looks like success.
- **Source:** KNOWN_BUGS.md "Regression 1 part 2" (originally cited
  `:218-223`; numbering shifted but the bug stands).
- **Action:** Set `out.error = "ORCA shielding load failed for " +
  opts.orca_nmr_path; return out;` at the failure site.

### OI-002 — APBS global C state does not reinitialise between calls

- **Verified:** `src/apbs_bridge.c` — per-call ctors and dtors are paired
  correctly, but FETK library-level globals (Vmem accounting, Vnm
  logging) are not reset between `apbs_solve()` invocations. Latent.
  Not exercised today because `ui/` / `h5-reader/` do not run the
  calculator pipeline; a session-reuse viewer would hit it.
- **Source:** `spec/plan/bones/DEPENDENCIES.md:62`.
- **Action:** Document as known limitation if any consumer ever runs
  APBS twice in one process. Real fix needs FETK-state-reset upstream.

---

## Open decisions (need user judgment)

### OI-010 — `CalculatorContract.h` adoption vs removal

- **Verified:** `src/CalculatorContract.h` (42 lines) exists; `grep -rn
  CalculatorContract src/ tests/` returns only the file itself. Zero
  calculator adopters across all 22+ result classes.
- **Source:** `spec/plan/openai-5.5-strong-architecture-layout.md:328-338,
  752-790`; also flagged in `review_items_to_assess.md` (T1 resolved as
  "typedef-as-documentation," not yet executed).
- **Action:** Either land typedef-as-discipline on each calculator
  OR remove the unused header.

### OI-011 — `ChiRotamerSelectionTrajectoryResult` is orphan

- **Verified:** `src/ChiRotamerSelectionTrajectoryResult.{h,cpp}` exist;
  `src/RunConfiguration.cpp:71-77` only attaches `BsWelford` to
  `ScanForDftPointSet`. Comment in the .cpp explicitly: "(ChiRotamerSelection
  et al.) are a pending-decision item."
- **Source:** `spec/plan/bones/pending_decisions_20260423.md` item 2 +
  `spec/plan/bones/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md` G4.
- **Action:** Decide if this is the right scan-for-DFT selection shape
  or design something else. Scan-for-DFT is ≥2 weeks out; can sit.

### OI-012 — `FullFatFrameExtraction` missing MOPAC ConformationResult dependencies

- **Verified:** `src/RunConfiguration.cpp:173-180` — inherits
  `PerFrameExtractionSet` (no MOPAC deps), only flips
  `skip_mopac = false`. Comment acknowledges the gap.
- **Source:** `spec/plan/bones/pending_decisions_20260423.md` item 3 (direction
  chosen: include MOPAC deps; not yet landed).
- **Action:** Wire `MopacResult`, `MopacCoulombResult`,
  `MopacMcConnellResult` into `required_conf_result_types_` BEFORE
  any MOPAC-family TrajectoryResult lands.

### OI-013 — `AllWelfords` revival on `TrajectoryAtom`

- **Verified:** `grep -n AllWelfords src/*.{h,cpp}` returns zero hits.
  `PATTERNS.md:1107-1116` still presents the pattern as live. Each
  `*WelfordTrajectoryResult::WriteH5Group` hand-spells dataset names —
  exactly the drift the pattern was meant to prevent.
- **Source:** `spec/plan/bones/pending_decisions_20260423.md` item 1.
- **Action:** Either revive AllWelfords on `TrajectoryAtom` or
  downgrade the §24 pattern in `PATTERNS.md` to "not revived in
  trajectory scope." Low priority; no blocking consumer yet.

### OI-014 — `--trajectory --analysis` is stubbed; AnalysisWriter dissolution incomplete

- **Verified:** `src/nmr_extract.cpp:309-315` — `RunAnalysis()` is
  hard-stubbed: *"disabled pending the dissolution of AnalysisWriter
  into per-Result H5 emitters."* `learn/bones/` holds the retired
  AnalysisWriter; per-result `WriteH5Group` is the target surface.
- **Source:** `spec/plan/bones/TRAJECTORY_REFACTOR_GAPS_2026-04-23.md` G6.
- **Action:** Implement `WriteH5Group` per-result emitters across
  attached TRs; remove the stub.

### OI-015 — Narrow `RunConfiguration` in trajectory tests slides suite away from production

- **Verified:** `tests/test_amber_streaming.cpp:84-86, 129-131` and
  `tests/test_frame_pdb_emitter.cpp:71-72` set `skip_apbs = true`,
  `skip_coulomb = true`, `skip_dssp = true`. None of these tests
  exercise the full `PerFrameExtractionSet`.
- **Source:** `spec/plan/bones/FIX_TESTS.md` §1 (§2 NVRTC fix landed 2026-04-25
  and unblocked the dread of moving these to full pipeline).
- **Action:** ~700-line refactor: build shared `ProductionTestSession`
  helper, convert sites, leave one MOPAC canary. Document baseline
  rule in `TEST_FRAMEWORK.md`.

### OI-016 — H5 schema lacks irrep / units / sign-convention metadata

- **Resolved:** commits `f2781da` + `dc50917` (2026-05-13 topology
  sidecar landing). `ArraySpec` gained `native_axis`, `irreps`,
  `units`, `sign_convention`, `tensor_rank`, `parity`, `mechanism`,
  populated for all ~108 CATALOG entries. R-side regex-mechanism
  refactor (OI-120) is the downstream consumer.

### OI-017 — Prochiral methylene H disambiguation in `LarsenResidue`

- **Verified:** `src/LarsenResidue.cpp:382-445` — methyl-H pseudoatom
  collapse implemented; no CIP perception for pro-R/pro-S yet.
  Spatial-tiebreak is current behavior.
- **Source:** `spec/plan/bones/larsen-residue-design-2026-05-11.md:486-524`.
- **Action:** Open only if calibration shows systematic pro-R/pro-S
  bias on methylene Hα/HB/HG/HD/HE pairs. Trigger is
  calibration-residual asymmetry.

---

## TODOs (work known to be needed)

### OI-030 — Layer-1 Welford `TrajectoryResult` clones (~9 classes)

- **Verified:** `ls src/*Welford*TrajectoryResult*` returns ONLY
  `BsWelfordTrajectoryResult`. HmWelford, McWelford, CoulombWelford,
  HBondWelford, PiQuadWelford, RingSuscWelford, DispersionWelford,
  TotalGWelford, AimnetPredictedWelford — none landed.
- **Source:** `spec/plan/bones/TRAJECTORY_RESULT_PLAN_2026-04-24.md:107-152` and
  `spec/plan/bones/pending_include_trajectory_scope_2026-04-22.md:3553-3604` (~23
  of 30 cataloged TR classes never landed; superseded in priority by
  the 2026-05 Larsen tripeptide work).
- **Action:** Add as Layer-1 fan-out when trajectory-scope work resumes.

### OI-031 — `TripeptideShieldingTimeSeriesTrajectoryResult` queued

- **Verified:** `src/Trajectory.cpp:170` comment notes
  `TripeptideShieldingTimeSeriesTrajectoryResult` is queued but not
  implemented; per-frame tensors land on `ConformationAtom` only,
  no H5/NPY surface mirroring `BsShieldingTimeSeriesTrajectoryResult`.
- **Source:** `spec/plan/bones/session_2026-05-10_part2_tripeptide_calculators.md:208-214`.
- **Action:** Mirror the `Bs*TimeSeries...` shape for tripeptide BB +
  Neighbor tensors when trajectory-scope tripeptide output is needed.

### OI-032 — `AmideTensorGeometry` calculator (A.2) blocked on Loth 2005 PDF

- **Verified:** `grep AmideTensorGeometry src/` returns nothing.
  Recommended frame in design: z along N–H, y normal to amide plane,
  x = y × z. Loth 2005 not in `references/`.
- **Source:** `spec/plan/bones/session_2026-05-10_calculator_queue_and_data_backend.md:118-128`.
- **Action:** Acquire Loth 2005, verify convention, then implement.

### OI-033 — Canonicalization JSON log writer not landed

- **Verified:** `grep "canonicalization_record\|CanonicalizationRecord"
  src/ProteinBuildContext.h src/FullSystemReader.cpp` returns empty.
  `GromacsToAmberReadbackBlock` covers the protein-side merge but
  does not emit the JSON canonicalization log.
- **Source:** `spec/plan/bones/naming-canonical-vocabulary-2026-04-30.md:648-703`
  steps NC1–NC7.
- **Action:** Emit
  `prep_run_*/.../gromacs_to_amber_readback_block.json` per the design.

### OI-034 — Multiple test sites build pipelines without AIMNet2

- **Verified:** `grep "Aimnet2Model()\|aimnet2_model"
  tests/test_calculation_runner.cpp tests/test_pipeline_and_sample.cpp
  tests/test_write_features.cpp` returns empty. Full tree census not
  yet done.
- **Source:** `spec/plan/bones/session-handoff-20260504.md:136-170`.
- **Action:** Census the tests/ tree, build a shared
  `ProductionTestSession` helper, convert the sites, leave one MOPAC
  canary in the one-shot set.

### OI-035 — `tests/golden/blessed/fleet/` still tracked despite gitignore

- **Verified:** `ls tests/golden/blessed/fleet/` shows `frame_001`;
  `git ls-files tests/golden/blessed/` returns 53+ tracked files
  including the fleet subtree. `.gitignore` is in place but does not
  retroactively untrack.
- **Source:** `spec/plan/bones/session-handoff-20260504.md:173-180` (P3).
- **Action:** `git rm --cached -r tests/golden/blessed/fleet/`,
  commit separately for readable diff.

### OI-036 — `tests/regression/run_regression.sh` orphan

- **Verified:** Script exists; `tests/regression/output_fleet/` has
  frame artifacts; `grep add_test` returns no references; Use Case B
  invocation at `:46` is also broken (calls `--orca "$ORCA_DIR"`
  without `--root NAME`, current `JobSpec` rejects).
- **Source:** `spec/plan/bones/session-handoff-20260504.md:179-184` (P3).
- **Action:** `git rm tests/regression/run_regression.sh` and the
  `output_fleet/` artifact dir.

### OI-037 — Tighten `min_npy_files` thresholds in `tests/test_smoke.cpp`

- **Verified:** `tests/test_smoke.cpp:433`
  `RunSmoke("nodft", conf, opts, 14, 40, ...)` and `:478`
  `RunSmoke("withdft", conf, opts, 15, 45, ...)`. Actual outputs are
  57 and 60; thresholds at 40 and 45 catch only catastrophic regression.
- **Source:** `spec/plan/bones/session-handoff-20260504.md:200-203` (P4).
- **Action:** Set to ~95% of expected (55 / 58) to catch partial-output
  regressions.

### OI-038 — `tests/BlessCompare.cpp::LoadTable` missing warn-once for missing `bless_policy.toml`

- **Verified:** `tests/BlessCompare.cpp:389` `LoadTable` returns empty
  default; `grep warn tests/BlessCompare.cpp` returns nothing. Comment
  at `BlessCompare.h:69` promised the warn line.
- **Source:** `spec/plan/bones/session-handoff-20260504.md:209-212`.
- **Action:** Add the warn-once log line.

### OI-039 — `tests/external/1OKH_4587_protonated.pdb` FF provenance unresolved

- **Verified:** `tests/test_protonation_detection.cpp:48,51` and
  `tests/test_foundation_results.cpp:564,569` still consume
  `GmxProtonated()`; `tests/testpaths.toml:11` maps it. CHARMM-era
  vs AMBER-era origin undocumented.
- **Source:** `spec/plan/bones/test_inventory_2026-05-03.md:22, 210-213`.
- **Action:** Verify provenance and document, or retire if CHARMM-era.

### OI-040 — SDK test fixture missing for water / hydration / gromacs_energy NPYs

- **Verified:** `find tests/data/ -name "water*" -o -name "hydr*"
  -o -name "gromacs_energy*"` returns zero. Catalog entries exist in
  `_catalog.py` but no fixture covers them.
- **Source:** `spec/OUTSTANDING_GROMACS_PATH.md:56-62`.
- **Action:** Run a small trajectory extraction, commit the
  water/hydration/gromacs_energy NPYs as a fixture.

### OI-041 — `HydrationGeometryResult` and `EeqResult` lack trajectory-level accumulation

- **Verified:** Both exist as `ConformationResult` classes; no
  `Accumulate`/trajectory-handler wiring. Integration target is now
  `TrajectoryResult` (not the retired `GromacsProtein.AccumulateFrame`).
- **Source:** `spec/plan/bones/TRAJECTORY_EXTRACTION.md:97-98, 127-131`.
- **Action:** Build per-frame `*TrajectoryResult` wrappers when
  trajectory accumulation is needed.

### OI-042 — Ring-normal stability audit pending (Sahakyan-Vendruscolo 2013)

- **Verified:** `src/BiotSavartResult.cpp` and `src/HaighMallionResult.cpp`
  use single-normal ring computation. Two-normal averaging (their fix
  for out-of-plane fluctuations in MD) is not implemented.
- **Source:** `spec/plan/bones/md-rerun-685-discussion-priors-2026-04-30.md:226-235`.
- **Action:** Pre-MDP-lock audit on a 2-protein subset; histogram Δθ
  between independent ring normals.

### OI-044 — `BsShieldingTimeSeriesTrajectoryResult::Finalize` latent UB

- **Verified:** `src/BsShieldingTimeSeriesTrajectoryResult.cpp:62-88`
  uses the same swap-then-leave-n_frames pattern that produced the
  Finalize-idempotency UB in `TripeptideBackbone/NeighborShieldingTimeSeries`
  (fixed at `bb657f5`). Bs's pattern reads `src[f]` on empty
  `per_atom_shielding_[i]` on any second `Finalize` call.
- **Source:** discovered during pilot cleanup; same code shape, fixed
  in the Tripeptide TRs via bounds-check, not yet propagated to Bs.
- **Action:** apply the same bounds-check pattern when Bs is next
  touched. No test currently exercises double-Finalize on Bs, so the
  bug is latent.

### OI-045 — Tripeptide / Larsen TR bundle (RESOLVED 2026-05-15)

- **Status:** **RESOLVED.** All 10 queued TRs landed across 6 commits
  on 2026-05-15:
  - `5c7488e` BB residual_vec (Vec3, establishes pattern)
  - `0279a97` Neighbor residual_vec prev + next (Vec3 ×2; NaN-tolerant)
  - `ccfd60b` BB method_tag (uint8, establishes int-buffer pattern)
  - `32a96aa` Larsen water_term (double, establishes scalar pattern)
  - `3b27497` Larsen count (int)
  - `48d5788` Larsen 1pHB/2pHB/1pHaB/2pHaB shielding (4 SphericalTensor TRs)
- **Tests:** 50 new tests across 7 test targets (5 tests per TR ×
  10 TRs), all green on 1P9J fleet_amber. Pattern 15 discipline trio
  (synthetic + frame-0 + idempotency + H5 round-trip) + integration
  with log-overages per `feedback_log_overages_dont_assert`.
- **Docs:** OBJECT_MODEL.md catalog updated (both `Known types` table
  at line 1666 and `PerFrameExtractionSet` table at line 2524).
  Memory entry `project_tripeptide_tr_pilot_landed` updated.
- **Remaining:** Adversarial codex + opus review on the bundle
  (task #133).

### OI-043 — `tests/data/fleet_amber/_backup_round{1,2}_*` disk hygiene

- **Verified:** Both `_backup_round1_15ns_2.5nm_padding_20260501/` and
  `_backup_round2_15ns_bumped_padding_spec_cadence_20260501/` still
  on disk under `tests/data/fleet_amber/`.
- **Source:** `spec/plan/bones/test_inventory_2026-05-03.md:225-227`.
- **Action:** Remove (these are the rejected pre-Option-B MDP attempts;
  current canonical is the `optB` subdirectory).

---

## Gotchas (won't break code but bite if violated)

### OI-050 — Running `./build/<binary>` directly bypasses AIMNet2 / NVRTC env

- **Verified:** `scripts/run_with_cuda_env.sh` exists; `CMakeLists.txt:374,406`
  applies it via `set_tests_properties(... ENVIRONMENT ...)` for ctest.
  Direct binary invocation silently fails AIMNet2 path at first inference.
- **Source:** `spec/plan/bones/TEST_FRAMEWORK.md:53-57`, `feedback_test_invocation_via_ctest`
  memory.
- **Action:** Use `ctest -R` or `scripts/run_with_cuda_env.sh
  ./build/<binary>`. Discipline only.

### OI-052 — Disulfide-authority populator rule fragility

- **Verified:** `src/FullSystemReader.cpp:925` is currently the only
  populator of `LegacyAmberInvariants::disulfide_pairs`. The
  AMBER-PRMTOP-direct load path is not yet implemented; when it is,
  it must also populate. If it forgets, geometric SG-SG inference
  silently re-takes authority.
- **Source:** `review_items_to_assess.md:111-130`.
- **Action:** Add a code-side comment in `LegacyAmberTopology.h` near
  the disulfide_pairs field flagging populator discipline; consider
  a runtime assert when any AMBER load path completes.

### OI-053 — Larsen H-bond grid frame convention (parser-side)

- **Verified:** `src/LarsenHBondGrid.{h,cpp}` enforce the
  donor-frame canonicalization; parser at
  `scripts/larsen_hbond_grid_parse/parse_larsen_hbond_grids.py`
  documents `ρ` is sign-flipped relative to standard IUPAC.
  Get this wrong once and every contribution has silent angular garbage.
- **Source:** `spec/plan/bones/larsen-hbond-shielding-design-2026-05-11.md:298-319`.
- **Action:** Discipline only — verified covered by the canonicalization
  helpers; future grid extensions must use the same convention.

### OI-054 — Event-extractor duplication risk vs `ChiRotamerSelectionTrajectoryResult`

- **Verified:** `src/ChiRotamerSelectionTrajectoryResult.h:42` exists
  (orphan per OI-011). `RingFlipEventTrajectoryResult`,
  `RotamerTransitionEventTrajectoryResult`, `HBondEventTrajectoryResult`
  do NOT exist in src/.
- **Source:** `spec/plan/bones/DIAGNOSTICS_AND_WORKFLOWS_2026-05-09.md:224-229`.
- **Action:** Anyone writing a rotamer-transition event extractor
  must audit `ChiRotamerSelectionTrajectoryResult` for overlap first.

---

## Doc/code drift (where docs lie about code)

### OI-070 — `OBJECT_MODEL.md` "Compute vs Attach" sloppy paraphrase

- **Verified:** `OBJECT_MODEL.md:1462` says computation happens in
  the static `Compute()` factory. `:1472` (two lines later)
  colloquialises this as "the attach method is where computation
  happens." `src/ProteinConformation.cpp:22-56` confirms AttachResult
  does singleton check + dep check + store only. Not a real
  contradiction; sloppy paraphrase that could mislead a careful reader.
- **Source:** `spec/plan/bones/doc_wrongness_20260423.md`.
- **Action:** Tighten `OBJECT_MODEL.md:1472` to clarify Compute is
  the factory, AttachResult merely stores. One-line edit.

### OI-071 — `OBJECT_MODEL.md` "md.tpr" vs current "production.tpr"

- **Verified:** `OBJECT_MODEL.md:1567` and `:2683` both say `md.tpr`;
  `src/TrajectoryProtein.cpp:42` uses `production.tpr`. Path
  convention was changed during the AMBER readback work (CHARMM-era
  used `md.*`, AMBER Option B uses `production.*`).
- **Source:** `review_items_to_assess.md:163-168`.
- **Action:** Two-line edit to `OBJECT_MODEL.md`.

---

## Known limitations (accepted scope; calibration absorbs)

### OI-090 — Larsen H-bond reference subtraction uses r-max-edge proxy

- **Verified:** Per `project_larsen_acceptor_reference_subtraction`
  memory and `spec/plan/bones/larsen-hbond-shielding-design-2026-05-11.md:741-755`.
  Bias 10-30× on N/CA/C acceptor tensors; calibration absorbs.
  True monomer DFT not in published Larsen archive.
- **Source:** Documented in the design doc + memory entry.
- **Revisit trigger:** Calibration residuals stratifying by
  acceptor class would point back here.

### OI-091 — Larsen sidechain primary amide acceptor approximation

- **Verified:** ASN OD1 / GLN OE1 use the NMA acceptor grid as a
  documented approximation. `src/SemanticEnums.h:934-940` carries
  the `IsSidechainAmideOxygen()` predicate; no separate sidechain-amide
  grid in Larsen's published set.
- **Source:** `spec/plan/bones/larsen-hbond-shielding-design-2026-05-11.md:813-816`.
- **Revisit trigger:** Direct sidechain-amide DFT scan if it ever
  becomes available; or if calibration residuals stratify the way.

---

## Resolved (kept for provenance)

- **OI-016** — `ArraySpec` irrep / units / sign-convention metadata.
  Resolved by topology sidecar landing 2026-05-13 (`f2781da` +
  `dc50917`). Still listed under "Open decisions" above with the
  resolved note in place.

---

## Pointers

- Architecture record: `spec/plan/openai-5.5-strong-architecture-layout.md`
- Operational reference: `project_proteintopology_architecture` memory entry
- Larsen calculator landing: `project_larsen_hbond_calculator` + `project_larsen_residue_model` memory entries
- 1UBQ-Stage-2 narrowing: `project_1p9j_study_system` memory entry
- Fleet status (685-stop, OF3 recovery): CLAUDE.md Stage 2 section + `project_three_stages` memory
