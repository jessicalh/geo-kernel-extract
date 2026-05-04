# Session handoff — 2026-05-02 evening

**Pickup brief.** Read this first; it points at what's done, what's pending, and the one open decision blocking the push.

---

## State at end of session

**Local HEAD:** `e343b5e` "AMBER readback: HISH unblock, disulfide authority, plain-fields LegacyAmberTopology, frame catch-all CR"

**origin/master:** `f8729e3` (post-revert clean; unchanged from previous session).

**Local commits ahead of origin (8):** `ce8a2e6`, `4c22528`, `86d3e12`, `1ea0dd0`, `a8c94d9`, `88d813e`, `590a754`, `e343b5e`. All landed clean locally; tests pass.

**Working tree:** uncommitted changes implementing a small follow-up gap fix (see below). Not yet committed because the push is blocked.

```text
Uncommitted (gap fix, ready to commit):
  M src/LegacyAmberTopology.h    has_disulfide_authority field added
  M src/FullSystemReader.cpp     sets has_disulfide_authority = (readback != nullptr)
  M src/Protein.cpp              FinalizeConstruction gates on bool, not vector emptiness
  M src/CovalentTopology.{h,cpp} OverrideDisulfides(empty) now means authority says zero
  M tests/test_object_model.cpp  regression tests for authority-zero vs no-authority
```

---

## What `e343b5e` accomplished

The block-as-architecture for AMBER trajectory loads. Concretely:

- **HISH unblock** — `GromacsToAmberReadbackBlock` reads pdb2gmx's chemistry decisions from `topol.top` rtp comments + the TPR bonded list, applies typed facts to the existing object model, dies. `--trajectory` AMBER loads now resolve `HISH/HISD/HISE/CYS-with-no-HG` to canonical chemistry without expanding `NamingRegistry`. The 3 previously HISH-blocked tests pass.
- **Disulfide authority** — `CovalentTopology::OverrideDisulfides` consumes the readback's SG-SG pairs as authority. Geometric SG-SG inference at `CovalentTopology.cpp:134` stops being source of truth on the trajectory path. For 1P9J: 3 authoritative pairs, geometric agreement, no overrides/demotions/force-adds — but the source-of-truth chain is now principled.
- **Plain-fields LegacyAmberTopology** — `AmberFFData` wrapper struct + Attach/Has/Optional + `MutableLegacyAmber` removed. Replaced by plain const fields populated through a `LegacyAmberInvariants` value-pack at construction.
- **Frame catch-all CR** — `GromacsFramePullResult` holds per-frame velocities + box (positions stay constitutive on conformation; per-frame source-specific data goes through a CR so PDB and GROMACS load paths share the same `ProteinConformation` shape).
- **Test reshape** — `test_gromacs_streaming.cpp` (CHARMM/XTC) retired; `test_amber_streaming.cpp` (AMBER/TRR via 1P9J fixture) replaces it. 5 tests, all pass.
- **Path convention** — `TrajectoryProtein::BuildFromTrajectory` uses `production.tpr` (was `md.tpr`); reads `topol.top` from `<production_dir>/../topol.top`.
- **Real bug fixes** (per OpenAI review) — variant-mapper indices for ASH/GLH/CYX/CYM/LYN, TRR velocity-presence detection (NaN sentinel, was unconditionally true), `LYSN` variant lookup, `<cstdio>` include.

**Test results at HEAD:**
- 68/68 wide regression sweep passes (AmberCharge*, AmberFlatTable*, AmberPreparation*, AmberLeapInput*, AmberPreparedCharge*, ChargeFF14SB*, ApbsFF14SB*, PrmtopChargeTest, OrcaRunTest, ChargesReturned, AmberTrajectoryFixture).
- 5/5 AmberStreaming tests pass (~17 minutes through the full PerFrameExtractionSet on 1P9J's TRR trajectory).
- 6/6 AmberTrajectoryFixture (was 3/6 before HISH unblock).

---

## The gap fix (uncommitted, ready)

OpenAI's 2026-05-02 evening review flagged: empty `LegacyAmberInvariants::disulfide_pairs` is currently overloaded — it means BOTH "no upstream authority" (PDB load) AND "authority says zero pairs" (TPR with no disulfides). This means the override can't demote a false geometric S-S when GROMACS explicitly says zero, and it can't distinguish extraction failure from a real zero.

**Fix:** add `bool has_disulfide_authority` to `LegacyAmberInvariants`. `FullSystemReader::BuildProtein` sets it `true` when `readback != nullptr` regardless of pair count. `Protein::FinalizeConstruction` gates `OverrideDisulfides` on the bool, not on `!disulfide_pairs.empty()`. `CovalentTopology::OverrideDisulfides(empty)` no longer returns early; once the caller has established authority, an empty pair list means "authority says zero disulfides", so any geometry-derived `BondCategory::Disulfide` tags are demoted.

Regression tests added in `tests/test_object_model.cpp`:

- `ObjectModel.DisulfideAuthorityZeroDemotesGeometry`: synthetic close SG-SG geometry + authority-present empty list -> zero disulfide bonds.
- `ObjectModel.NoDisulfideAuthorityKeepsGeometry`: same geometry + no authority -> one geometric disulfide bond.

Builds clean. Focused tests pass:

```text
cmake --build build --target unit_tests trajectory_tests -j2
./build/unit_tests --gtest_filter=ObjectModel.DisulfideAuthorityZeroDemotesGeometry:ObjectModel.NoDisulfideAuthorityKeepsGeometry
./build/trajectory_tests --gtest_filter=AmberTrajectoryFixtureTest.*
```

**Test footing still pending for 2026-05-03.** The gap fix now has a
synthetic regression and the AMBER trajectory fixture suite still passes,
but this should get a firmer test footing before the topology/readback
work is considered settled:

- finalized 1P9J topology should assert exactly 3 `BondCategory::Disulfide`
  bonds, matching 6 CYX residues and 3 TPR authority pairs;
- finalized 1Z9B topology should assert authority-present/zero disulfides
  and zero final `BondCategory::Disulfide` bonds;
- add or keep direct disagreement regressions for all authority outcomes:
  force-add missing geometric bond, override wrongly categorized geometric
  bond, and demote false geometric disulfide;
- add a construction-boundary consistency guard or test for CYX endpoint
  count vs authoritative pair count, so partial extraction cannot silently
  look like a legitimate zero.

Held uncommitted because the push is blocked and we don't want to compound.

---

## The push blocker

Pushing to `origin/master` (`https://github.com/jessicalh/geo-kernel-extract.git`) fails:

```
remote: warning: File data/ccd/components.cif.gz is 88.39 MB ...
remote: error: File data/ccd/components.cif is 375.69 MB; this exceeds GitHub's file size limit of 100.00 MB
remote: error: GH001: Large files detected.
! [remote rejected] master -> master (pre-receive hook declined)
```

**Cause:** commit `ce8a2e6` ("intermediate checkin: AMBER slice + PHASE 0 partial + day's planning + naming spec", 2026-05-01) committed `data/ccd/components.cif` (376 MB) and `data/ccd/components.cif.gz` (89 MB). Both exceed GitHub's 100 MB hard limit. The push is rejected.

`git filter-repo` is installed at `/home/jessica/.local/bin/git-filter-repo`.

**User direction (this session, end):** "this data we can keep locally that does not need to be in github." → remove from history, keep locally. "surgical work is not ideal" → use the cleanest non-tedious path.

**Cleanest path next session:**

```
# 1. Commit the gap fix as its own atomic commit (small, no CCD touch).
git add src/LegacyAmberTopology.h src/FullSystemReader.cpp src/Protein.cpp \
        src/CovalentTopology.h src/CovalentTopology.cpp tests/test_object_model.cpp
git commit -m "$(cat <<EOF
LegacyAmberInvariants: has_disulfide_authority flag

Distinguishes "no upstream authority" (PDB load, default) from
"authority says zero pairs" (TPR with no disulfides). Prior shape
overloaded empty disulfide_pairs to mean both, so a false geometric
S-S could not be demoted when GROMACS explicitly recorded zero.

  src/LegacyAmberTopology.h
    bool has_disulfide_authority = false; on LegacyAmberInvariants,
    with comment explaining the three semantic states.

  src/FullSystemReader.cpp
    BuildProtein sets has_disulfide_authority = (readback != nullptr).

  src/Protein.cpp
  src/CovalentTopology.{h,cpp}
  tests/test_object_model.cpp
    FinalizeConstruction gates OverrideDisulfides on the bool, not
    on !disulfide_pairs.empty(). Empty + authority-true demotes any
    geometric Disulfide tag.

Doesn't bite our current fixtures (1P9J: 3 real disulfides; 1Z9B:
no CYS at all). Provenance fix; sets up correct behaviour for
no-disulfide AMBER inputs going forward.
EOF
)"

# 2. Add data/ccd/ to .gitignore, commit.
echo "data/ccd/" >> .gitignore
git add .gitignore
git commit -m "gitignore: data/ccd/ — PDB CCD reference data, kept local"

# 3. Filter-repo to drop the CCD files from history.
git filter-repo --path data/ccd --invert-paths --force

# 4. The local data/ccd/ files stay on disk (filter-repo doesn't
#    touch the working tree). The .gitignore prevents re-adding.

# 5. Re-add origin (filter-repo strips remotes by default).
git remote add origin https://github.com/jessicalh/geo-kernel-extract.git

# 6. Force-push.
git push --force-with-lease origin master
```

**Force-push warning:** `master` will be rewritten. Anyone else with a clone of this remote will need to reset their local. There are no other collaborators on this repo per project memory; the only consumers are the user's own machines.

**Alternative if you want to defer the rewrite:** push to a side branch on a personal remote with LFS enabled (option 4 from this session's discussion). Lower risk; doesn't preserve to `origin/master` until later.

---

## What is decided vs open

**Decided this session (in `e343b5e` + the uncommitted gap fix):**

- The `GromacsToAmberReadbackBlock` exists, parses `topol.top`, applies typed facts, dies. Compiler trace, not live state.
- `LegacyAmberTopology` holds invariant FF-numerical data + disulfide pairing as plain fields populated at construction via the value-pack.
- `GromacsFramePullResult` is the catch-all for per-frame GROMACS data. Not on ProteinConformation as direct fields.
- `production.tpr/.trr/.edr` is the canonical naming. `md.*` is legacy.
- Disulfide pairing comes from TPR bonded list + readback CYX, not geometric inference, on the trajectory path.
- The verification triad (NPY identity + UDP log + GeometryChoice) stays untouched — it IS the project's evidence layer.

**Open / next session:**

1. **Push blocker resolution.** Per user: filter-repo + force-push as above. Once unblocked, commit the gap fix and push.
2. **OpenAI #2 unfinished.** `nmr_extract.cpp` still passes `md.edr` to `Trajectory`. JobSpec was updated in commit `590a754` to use `production.*`; `nmr_extract.cpp` was not. Small, follow-up.
3. **Cap pseudo-residues** (NHE/ACE/NME) — `AminoAcid` enum extension if needed for future fixtures. Not blocking current work.
4. **Orphan `ConformationAtom` schema fields** — user direction recorded in `running_plan_notes.md`: "we probably want those." Implement producers, don't drop schema. Not blocking.
5. **Duplicated PRMTOP parsing** between prepared AMBER and existing charge paths — local cleanup. Not blocking.

**NOT touched (and shouldn't be without explicit per-change agreement):** `OBJECT_MODEL.md`, `PATTERNS.md`. The session learned the hard way: foundational doc edits without per-change agreement are the slop pattern.

---

## Memory entries that document the architecture

Loaded automatically next session. The four most relevant for this work:

- `feedback_capture_at_the_boundary` — source-available data is not speculative; capture all of it at every read boundary
- `feedback_no_attach_lifecycle_for_invariant_data` — invariant data on typed objects as plain fields populated at construction; no Attach/Has/Optional/Mutable accessors
- `feedback_readback_block_is_a_compiler_trace` — `GromacsToAmberReadbackBlock` is import-time-only; calculators see only compiled facts
- `project_charmm_retired_amber_only_2026-05-02` — project state: AMBER-only via TRR + libgromacs-direct; quarantined-legacy code listed

---

## Files listing

```
Committed at HEAD `e343b5e` (28 files, 1924 inserts, 1636 deletes):
  CLAUDE.md                                                          modified (additive 2026-05-02 update block)
  CMakeLists.txt                                                     modified (test target rename)
  spec/plan/gromacs-to-amber-readback-block-design-2026-05-02.md     new
  src/AminoAcidType.{cpp,h}                                          modified
  src/CovalentTopology.{cpp,h}                                       modified (DisulfidePair + OverrideDisulfides)
  src/FullSystemReader.{cpp,h}                                       modified (readback + disulfide scan)
  src/GromacsFrameHandler.cpp                                        modified (NaN sentinel)
  src/GromacsFramePullResult.{h,cpp}                                 new
  src/GromacsToAmberReadbackBlock.{h,cpp}                            new
  src/LegacyAmberTopology.{cpp,h}                                    modified (plain fields)
  src/OperationRunner.{cpp,h}                                        modified (catch-all CR attach)
  src/Protein.{cpp,h}                                                modified (FinalizeConstruction value-pack)
  src/ProteinConformation.h                                          modified (vel/box removed)
  src/RunConfiguration.cpp                                           modified
  src/Trajectory.{cpp,h}                                             modified (env_ vel/box)
  src/TrajectoryProtein.cpp                                          modified (production.* + readback parse)
  tests/test_amber_streaming.cpp                                     new (replaces test_gromacs_streaming.cpp)
  tests/test_amber_trajectory.cpp                                    modified (parses topol.top + passes block)
  tests/test_gromacs_streaming.cpp                                   DELETED (CHARMM/XTC retired)

Uncommitted (gap fix):
  src/LegacyAmberTopology.h
  src/FullSystemReader.cpp
  src/Protein.cpp

Memory entries (in ~/.claude/.../memory/, not in repo):
  feedback_capture_at_the_boundary.md                                new
  feedback_no_attach_lifecycle_for_invariant_data.md                 new
  feedback_readback_block_is_a_compiler_trace.md                     new (renamed from feedback_canonicalization_record_is_a_compiler_trace.md)
  project_charmm_retired_amber_only_2026-05-02.md                    new
  MEMORY.md                                                          updated (4 new index lines)
```

---

## Don't forget

- **Don't try to push without addressing the CCD blocker first.** Same error every time.
- **Don't commit the gap fix until the rewrite path is decided** — the gap fix would just become a 9th commit ahead of origin, unhelpful until the path forward is settled. (Alternatively: commit it BEFORE filter-repo, since filter-repo doesn't touch source-only commits.)
- **Don't touch OBJECT_MODEL.md or PATTERNS.md** without explicit per-change agreement. The session learned this the hard way.
- The `Testing/` directory in the build tree is a CTest temp directory, not source. Already not staged.

---

## Suggested first 30 seconds of next session

1. Read this file.
2. `git status` to confirm the uncommitted gap fix is intact.
3. `git log --oneline origin/master..HEAD` to confirm the 8 ahead.
4. Decide: filter-repo + force-push, or another option from this doc.
5. Then proceed.
