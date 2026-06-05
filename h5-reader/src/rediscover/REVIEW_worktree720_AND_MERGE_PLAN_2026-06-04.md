# Review: worktree-720 static ingest + safe consolidation plan

Scope honored: reviewed `/shared/2026Thesis/nmr-720` only as code under review, did not edit that tree or `src/`, used git read-only, did not investigate the IUPAC-topology episode, and re-ran the oracle-parity gate.

## Part 1 - Findings

### CRITICAL - DFT dia/para/total identity is not enforced componentwise

`RunData.cpp:150-163` loads the canonical static target through `DftShieldingLoader::LoadAndValidate`, then `PerAtomSubstrate.cpp:4314-4329` emits raw total/dia/para and decomposed T0/T2 as the training target. The loader only checks `total == dia + para` through `a.total.T0` at `DftShieldingLoader.cpp:95-104`; it never checks the raw 3x3 components that are later exported. It also does not compare parsed DFT element order to topology element order, despite `DftShielding.h:20-24` carrying the element for that cross-check.

Why it matters: a tensor component/sign/order corruption can pass when the trace is right. That would poison T2 targets silently while leaving T0 and atom count apparently valid.

Fix: in `DftShieldingLoader::LoadAndValidate`, for every atom, compare all nine `total_raw(i,j)` components against `dia_raw(i,j) + para_raw(i,j)` within the ORCA print tolerance; also compare `a.element` to `protein->atom(i).element` and fail static target load on mismatch.

### CRITICAL - The strict static topology path still soft-fails bonds/rings/ring membership

`QtProteinLoader::LoadPoseStrict` calls the existing sidecar loader at `QtProteinLoader.cpp:229-247`. That sidecar still treats `bonds.npy` as soft-fail at `QtTopologySidecar.cpp:476-493`, `rings.npy` as soft-fail at `QtTopologySidecar.cpp:495-521`, and missing `ring_membership.npy` as warning-only at `QtTopologySidecar.cpp:533-540`, then returns `ok` at `QtTopologySidecar.cpp:580-581`. The strict NPY loader then uses topology-derived expected row counts; when ring topology is missing, `FrameNpyLoader.cpp:104-113` cannot enforce the ring-axis contract.

Why it matters: the 720 static substrate depends on ring/bond source geometry and ring membership. A pose can be accepted with empty rings/bonds, producing rows with missing or zero source structure instead of failing loudly.

Fix: add a strict sidecar mode for static rediscover that requires bonds, rings, and ring membership to exist, decode, and match manifest/topology counts. Use that mode from `LoadPoseStrict`; warning-only behavior can remain for the viewer path.

### HIGH - The rerun oracle passes, but it does not test the static ingest

I built/checked `h5-reader/build/linux-rwdi/h5reader_extract` (`ninja: no work to do`) and reran:

```text
python3 h5-reader/src/rediscover/analysis/oracle_parity.py \
  --bin h5-reader/build/linux-rwdi/h5reader_extract \
  --run /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft/1p9j-calibration-with-dft.LGS \
  --work /tmp/rediscover-parity-wt720-review --case all --mc-cutoff 8.0
```

Result: PASS. The gate reported identical ring/Mc CSVs and T2 NPYs: 141000 ring-source rows, 20500 ring aggregate rows, 812205 McConnell source rows, and 26000 McConnell aggregate rows.

But `oracle_parity.py:36-39` and `oracle_parity.py:107-121` only run `h5reader_extract --run ... --engine composed|procedural`. The new static path is behind `--static-cohort` at `main_extract.cpp:271-312`, so this gate never exercises `RunStaticCohortSubstrateEmit`, `RunLoader::LoadStaticPoseExact`, `FrameNpyLoader::LoadSnapshotExactFields`, or `StaticWriter`. The static manifest even records `oracle_parity.v1_ring_mc_mopac_mc_gate = "not_run_by_static_emit"` at `PerAtomSubstrate.cpp:4641-4643`, while also marking `stage2_between_axis_statistics_ready` and `stage3_model1_e3nn_training_corpus_ready` true at `PerAtomSubstrate.cpp:4591-4592`.

Why it matters: this oracle would still pass if the strict static path emitted wrong rows. It proves composed-vs-procedural trajectory parity, not 1P9J static-vs-trajectory equivalence.

Fix: add and run the contract oracle: one 1P9J trajectory frame and an equivalent static-pose producer output, joined on `(protein_id, atom_index)` vs `(protein_id, selected_h5_row, atom_index)`, with row-order-independent comparisons and explicit max deltas per gated family. Do not mark the static manifest training-ready until that gate passes.

### MED - Static poses never load FF14SB charges, but the static schema emits FF14SB fields

`RunData.cpp:178-187` loads `topol.top` FF14SB charges only for trajectory manifests. `CalcsetManifest.h:101-105` gives `SinglePose` no topology-top path, and `QtTopologySidecar` does not populate `QtAtom::partialCharge` (`QtTopologySidecar.cpp:401-425`, `QtAtom.h:65-66`). `Catalog.cpp:161-196` therefore marks `Ff14sbCharge` available only when topology atoms already have partial charges. The static rows/sidecars still emit FF14SB charge and field columns at `PerAtomSubstrate.cpp:4013-4024` and `PerAtomSubstrate.cpp:4296-4298`.

Why it matters: every static row will report FF14SB absent unless some upstream, non-modeled path injects charges. The electrostatic slab is then incomplete relative to the review-folded spec, and this can be mistaken for true missing physics rather than a loader omission.

Fix: add a lead-vetted static FF14SB charge source: either extend the single-pose manifest with an exact `topology_top` path and reuse `LoadFf14sbChargesFromTopol`, or require/validate an exact static FF14SB charge NPY. Record the source in `static_manifest.json`.

### MED - Production 720 cardinality is not enforced

`readStaticCohortManifest` requires non-empty `proteins[]` and unique IDs at `PerAtomSubstrate.cpp:3795-3816`, but it does not enforce the production contract of exactly 720 curated proteins. `RunStaticCohortSubstrateEmit` accepts that arbitrary count into `stats.protein_count` at `PerAtomSubstrate.cpp:4932-4938`.

Why it matters: a truncated cohort can produce a formally successful static substrate and even mark the manifest ready.

Fix: add an explicit run mode or manifest flag: production mode requires exactly 720 unique protein IDs; oracle/dev mode may allow one or a small controlled set and must label itself non-production.

### LOW - `nanoflann.hpp` in `wt-720-build` is a symlink

Confirmed: `/shared/2026Thesis/nmr-720/h5-reader/extern/nanoflann/nanoflann.hpp -> ../../../extern/nanoflann.hpp`. That violates the no-symlinks rule. The consolidation tree already has a real untracked vendored copy at `/shared/2026Thesis/nmr-shielding/h5-reader/extern/nanoflann/{COPYING,nanoflann.hpp}`.

Fix: vendor real files, including license, under `h5-reader/extern/nanoflann/`; do not commit the symlink.

## Sound paths checked

- Exact static field loading is directionally right: `FrameNpyLoader.cpp:60-61` resolves `<FieldSpec.stem>.npy` directly, and `FrameNpyLoader.cpp:241-247` fails missing/malformed required fields without directory enumeration.
- Canonical DFT source is raw ORCA via `.LGS` DFT meta, not mutation deltas: `RunData.cpp:138-167` requires a loadable static DFT frame 0, and `PerAtomSubstrate.cpp:4255-4258` records `raw_orca_out`.
- T2 is not scalar-collapsed in the emitter: raw target T2, local-frame T2, frame axes, and invariant norms are written at `PerAtomSubstrate.cpp:4317-4340`; local T2 is gated through `targetLocalT2` at `PerAtomSubstrate.cpp:4044-4048`.
- Static rows carry all-atom identity, BMRB/IUPAC labels, element, role, and stratum at `PerAtomSubstrate.cpp:4385-4406` and `PerAtomSubstrate.cpp:4248-4310`.
- Memory strategies remain C++-side: `all_resident` keeps `StaticResidentBundle`s in memory at `PerAtomSubstrate.cpp:4958-4964`; `streaming` reloads per protein and emits through the same C++ writer at `PerAtomSubstrate.cpp:5110-5115`. No Python second model is introduced.
- The disk gate exists and fails over 15 GiB at `PerAtomSubstrate.cpp:4972-4978`.

Verdict: not safe to merge as-is for production 720 training corpus. Fix the CRITICAL DFT validation and strict topology enforcement first, add a meaningful static-vs-trajectory oracle gate, and replace the nanoflann symlink with real vendored files before merge.

## Part 2 - Inventory

Read-only commands run: `git worktree list`, `git branch -vv`, `git status --short --branch`, `git diff --stat`, ancestry checks, and non-mutating log/show summaries.

### Worktrees and branches

- `/shared/2026Thesis/nmr-shielding` on `h5-reader-pysr-spike` at `dfc8f51`, ahead of origin by 49. Dirty: `.gitignore`, `CMakeLists.txt`, calculator docs, `h5-reader/CLAUDE.md`, `h5-reader/src/rediscover/analysis/equiv_t2_e3nn.py`, plus many untracked docs/scripts, the brief, this report, ORCA pilot files, `tools/apbs_radii_ab_workaround.cpp`, and real untracked `h5-reader/extern/nanoflann/{COPYING,nanoflann.hpp}`. Classification: real target work; preserve before merges.
- `/shared/2026Thesis/nmr-720` on `wt-720-build` at `04fcb3f`, which is already an ancestor of `h5-reader-pysr-spike`. Dirty: the 12 static ingest files plus untracked symlinked `h5-reader/extern/`. Classification: real uncommitted static ingest; not merge-safe until findings above are fixed.
- `worktree-agent-a85186b89732f6e97` at `10e7dd1` is clean and NOT an ancestor of `h5-reader-pysr-spike`. Tip commit: "Commit C: tensor Welford TRs - dxdt + per-dataset units + frame_index_range"; touches `src/BsWelfordTrajectoryResult.*`, `src/HmWelfordTrajectoryResult.*`, and `src/McConnellWelfordTrajectoryResult.*`. Classification: real committed work to consolidate if lead still wants Welford C.
- `worktree-agent-adb024158ea589895` branch tip `d22d04a` is already an ancestor of target, but the worktree has uncommitted Welford edits in `src/EeqWelfordTrajectoryResult.*`, `src/HBondCountWelfordTrajectoryResult.*`, and `src/SasaWelfordTrajectoryResult.*` (408 insertions, 225 deletions). Classification: real uncommitted work; lead decision required.
- `worktree-agent-a2320bf1880a690d5` at `9a1eeb1` and `worktree-agent-a6e4e859818582d26` at `ad17e07` are clean and already ancestors of target. Classification: consolidated; no merge needed.
- `worktree-agent-a0a9b04b`, `worktree-agent-a8d1b7b8`, `worktree-agent-abbddb5f`, and `worktree-agent-ae011756` all point to `d036bea` (viewer MOPAC work), already an ancestor of target; several paths are prunable/missing. Classification: duplicate/superseded worktrees; lead may prune later.
- `worktree-agent-ab2fe753eab00021f` at `ee40983`, `master` at `65c9c47`, and `backup-before-filter-repo` at `cf7796f` are already ancestors of target. Classification: already consolidated or archival; no merge needed.
- Detached `/tmp/bitdiff/precleanup` and `/tmp/h5review` are prunable/missing. Classification: prunable after lead review; no merge work.

## Proposed safe consolidation sequence

This is a plan only; do not execute git from this review.

1. In `/shared/2026Thesis/nmr-shielding`, preserve the current dirty target state first: commit it on `h5-reader-pysr-spike` or stash it with a named stash. Include this report and the real nanoflann vendoring if the lead wants them retained. Expected conflicts: none; this is just protecting existing work.
2. In `/shared/2026Thesis/nmr-720`, fix the blockers before committing: componentwise DFT identity + element alignment, strict topology sidecar mode, static-vs-trajectory oracle, FF14SB static charge source or explicit absence decision, production 720 cardinality mode, and replace/omit the nanoflann symlink. Re-run the existing ring/Mc oracle and the new static oracle. Expected conflicts: none inside the worktree; this is local cleanup before consolidation.
3. Commit the cleaned `wt-720-build` work on `wt-720-build`. Do not commit the symlink. Expected merge conflict into target after step 1: low, because `04fcb3f` is an ancestor of `dfc8f51` and the static ingest touches `h5-reader/src/io`, `h5-reader/src/model/QtProtein.h`, and `h5-reader/src/rediscover` code while `dfc8f51` was docs-focused.
4. Merge `wt-720-build` into `h5-reader-pysr-spike` from a clean target tree, or cherry-pick the one cleaned static-ingest commit if the lead wants tighter history. No force, no rebase required. Verify build plus both oracle gates after merge.
5. Consolidate `worktree-agent-a85186b89732f6e97`: merge or cherry-pick `10e7dd1` onto `h5-reader-pysr-spike`. Expected conflicts: low to medium; it touches Welford result files and may interact with target CMake/test evolution, but its branch parent prep is already in target.
6. Decide on `worktree-agent-adb024158ea589895` uncommitted Welford edits. If kept, commit them in that worktree first, then merge/cherry-pick onto target after `10e7dd1` so all six Welford families align. Expected conflicts: low with `10e7dd1` because it touches Eeq/HBond/Sasa while `10e7dd1` touches Bs/Hm/McConnell, but shared patterns/tests may need review.
7. For already-contained agent branches (`d036bea`, `9a1eeb1`, `ad17e07`, `d22d04a`, `ee40983`) and ancestor branches (`master`, `backup-before-filter-repo`), do not merge anything. After the lead confirms no local-only files are needed, she may prune stale worktrees/branches. No deletion should be automatic.
8. After all chosen work lands, run a final clean-tree gate set: configure/build `h5reader_extract`, existing ring/Mc `oracle_parity.py`, the new static oracle, and targeted Welford tests if Welford branches are merged.

Open decisions for lead:

- Whether `worktree-agent-a85186...` Welford C is still desired.
- Whether to keep and finish the uncommitted `agent-adb024...` Eeq/HBond/Sasa Welford edits.
- Whether current dirty target docs/scripts/ORCA pilot files are to be committed, stashed, or split.
- Whether production static mode must hard-require exactly 720 immediately, or use an explicit `oracle/dev` exception.
- When to prune duplicate/prunable worktrees; do not delete them as part of the merge.
