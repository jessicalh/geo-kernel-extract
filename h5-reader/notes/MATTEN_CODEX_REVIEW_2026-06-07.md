# Codex review: MatTen predictor plan, 2026-06-07

Scope: read-only review of `h5-reader/notes/MATTEN_PREDICTOR_PLAN_2026-06-07.md`.
I checked the plan's code assertions against the repository, the 1P9J ORCA tree,
and the 720-WT worktree at a shallow level. I did not edit source or the plan.

## Executive read

The plan is solid on the reader-side mechanics: `target_T1` is computed upstream
and then dropped by the all-atom equivariant sink; there is no lab-frame per-atom
`positions` sidecar; the all-atom equivariant target is cleanly lab-frame and
not a per-atom local-frame target; and both proposed emits are additive against
data already resident in the reader.

The main correction is scope. I found strong local evidence for the 1P9J
full-tensor DFT campaign, and code support for static 720-style ORCA tensors when
the `_nmr.out` files exist. I did not find local evidence that a populated
720-WT full-tensor DFT shielding corpus currently exists in the checked 720
worktree. Until that is audited, v1 should be framed as a 1P9J single-protein
conformational-generalization pilot, not a transferable leave-proteins-out
predictor.

The second correction is guardrail strength. The "same frame/order by
construction" claim is mostly supported, but the reader DFT loader does not
currently compare ORCA elements against topology elements. That check exists in
the static producer-side ORCA path, not in the reader trajectory DFT path. I
would keep the element check as a required fail-loud guard before any 720 or
re-prepped data is trusted.

## 1. Reader-side emits

### `target_T1` is genuinely absent from the all-atom equivariant sink

Confirmed. The all-atom target CSV header includes `dft_sigma_iso`, the 3x3 raw
tensor, and `dft_total_T2_0..4`, but no `dft_total_T1_*`
(`h5-reader/src/rediscover/AllAtomEquivariantSink.cpp:93-100`). The sidecar list
contains `target_T2`, `target_sigma_iso`, and `target_raw`, but no `target_T1`
(`h5-reader/src/rediscover/AllAtomEquivariantSink.cpp:146-160`). Commit writes
those same sidecars and no T1 sidecar
(`h5-reader/src/rediscover/AllAtomEquivariantSink.cpp:318-356`).

T1 is not missing upstream. `BuildTarget` copies the raw ORCA matrices and runs
`DecomposeLibrary` for total/dia/para
(`h5-reader/src/rediscover/ExtractionSupport.cpp:45-60`). `DecomposeLibrary`
computes `T1` as the Levi-Civita dual of the antisymmetric part
(`h5-reader/src/rediscover/SphericalBasis.cpp:13-17`). So the plan's "T1 dropped
silently" statement is correct for consumers that use the named irrep sidecars.

One implementation detail to preserve: the current sink uses a split predicate.
`dftPresent` is `row.target.present && finiteT2(...)`
(`h5-reader/src/rediscover/AllAtomEquivariantSink.cpp:183-187`), but sigma iso
and raw are appended under only `row.target.present`
(`h5-reader/src/rediscover/AllAtomEquivariantSink.cpp:210-215`,
`h5-reader/src/rediscover/AllAtomEquivariantSink.cpp:246-250`). The plan is
right to unify the presence gate, but I would define it as `target.present &&
finite(total_raw) && finite(T0/T1/T2)` rather than only `finite(total_raw)`.

### There is no lab-frame per-atom `positions` array emitted by this sink

Confirmed. The target record carries identity fields and `DftTarget`, but no
position field (`h5-reader/src/rediscover/AllAtomEquivariantSink.h:25-37`). The
source record carries only relational geometry such as `disp`, `r`,
`orientation_a`, and `orientation_b`
(`h5-reader/src/rediscover/AllAtomEquivariantSink.h:68-129`). The target sidecar
vectors are target T2, sigma iso, raw, APBS/AIMNet/MOPAC payloads, and no
positions vector (`h5-reader/src/rediscover/AllAtomEquivariantSink.h:162-174`).
The source CSV header likewise has `disp_x/y/z` and orientations, not node
positions (`h5-reader/src/rediscover/AllAtomEquivariantSink.cpp:120-136`).

The data is already resident. The catalog declares `ArrayId::Positions` as dense
H5 data when the conformation exists
(`h5-reader/src/rediscover/Catalog.cpp:76-82`), `verbs::pos` returns it
(`h5-reader/src/rediscover/Verbs.h:65-67`,
`h5-reader/src/rediscover/Verbs.cpp:17-26`), and trajectory positions are
required and shape-validated on reader load
(`h5-reader/src/io/QtTrajectoryH5.cpp:461-486`,
`h5-reader/src/io/QtTrajectoryH5.cpp:2312-2314`,
`h5-reader/src/io/QtTrajectoryH5.cpp:88-111`). So emitting row-aligned target
positions is an additive sink change, not a new reader architecture.

Adapter caveat: a `row_id` is one atom/frame target row, not one molecule. The
traversal is DFT rows outer, atoms inner
(`h5-reader/src/rediscover/RelationshipEngine.h:79-86`), and the all-atom
stratum is all atoms (`h5-reader/src/rediscover/AllAtomEquivariant.cpp:220-223`).
A `(rows,3)` position sidecar is usable, but the Python adapter must group rows
by `h5_row` or `original_index`, sort by `atom_index`, and build one whole
molecule per frame. I would state that explicitly in the plan and assert it in
the adapter tests.

### The local-frame target is confined to the broad-backbone sink

Confirmed. The all-atom equivariant sink says it is lab-frame and emits no
per-atom local frame (`h5-reader/src/rediscover/AllAtomEquivariantSink.h:1-6`).
The all-atom traversal returns a default invalid local frame
(`h5-reader/src/rediscover/AllAtomEquivariant.cpp:226-228`) and logs
`local_frame=none` (`h5-reader/src/rediscover/AllAtomEquivariant.cpp:687-690`).
`BuildTarget` only computes `total_local` when the frame is valid
(`h5-reader/src/rediscover/ExtractionSupport.cpp:62-66`), so the all-atom target
does not rotate DFT tensors into per-atom local frames.

The local T2 sidecar is in `BroadBackboneSink`: it declares
`*_aggregated_target_local_T2.npy`
(`h5-reader/src/rediscover/BroadBackboneSink.cpp:191-204`), computes it from
`rec.target.total_local` only when `local_frame_valid`
(`h5-reader/src/rediscover/BroadBackboneSink.cpp:317-325`), and writes it at
commit (`h5-reader/src/rediscover/BroadBackboneSink.cpp:353-358`). The plan's
quarantine claim is correct.

### Are the emits structurally harder than claimed?

No. The position emit and T1 emit are additive sink-side payloads. The only
design choice worth making up front is file shape:

- `positions.npy` as `(rows,3)` is simplest and matches target row order, but
  requires explicit frame grouping in Python.
- A frame-level `(dft_frames, atoms, 3)` sidecar would be more directly MatTen
  shaped, but would be a slightly different sink contract.

Given the existing CSV already carries `row_id`, `atom_index`, `h5_row`, and
`original_index`, row-aligned `(rows,3)` is fine if the adapter enforces the
grouping invariant.

## 2. Frame and atom-order consistency

The plan's construction story is mostly supported.

Producer-side trajectory mode configures per-frame PDB emission during traversal:
PDBs go to `output/pdbs`, NPYs to `output/npys`, and cadence is governed by
`--stride` (`src/nmr_extract.cpp:249-264`). The PDB filename encodes frame and
time (`src/FramePdbEmitter.cpp:102-108`). `EmitOnePdb` writes the current
conformation, iterating `ai = 0..natoms-1` and writing `conf.AtomAt(ai)` in
topology order (`src/FramePdbEmitter.cpp:176-220`,
`src/FramePdbEmitter.cpp:227-249`). `OnFrame` calls `EmitOnePdb` only for the
frames dispatched by the trajectory runner (`src/FramePdbEmitter.cpp:300-310`).

Reader-side DFT loading is also keyed rather than discovered by glob. `.LGS`
`dft.frames[]` entries carry `frame_index` and resolved `meta_json`
(`h5-reader/src/io/CalcsetManifest.h:39-51`,
`h5-reader/src/io/CalcsetManifest.h:110-120`), and manifest parsing requires a
`dft.frames` array with integer `frame_index` plus `meta_json`
(`h5-reader/src/io/CalcsetManifest.cpp:459-525`). `DftShieldingLoader` then reads
`files.out_primary` from that meta JSON and resolves it relative to the job dir,
with no globbing (`h5-reader/src/io/DftShieldingLoader.cpp:45-76`). `RunData`
loads each manifest DFT frame keyed by `frame_index`
(`h5-reader/src/rediscover/RunData.cpp:117-131`) and builds a frame map only
when the trajectory frame-index basis is `trr_frame_index`
(`h5-reader/src/rediscover/RunData.cpp:46-67`).

The ORCA parser captures the first Cartesian coordinate block in Angstroms
(`h5-reader/src/io/OrcaShieldingParser.cpp:46-70`), parses chemical shielding
nucleus index and element (`h5-reader/src/io/OrcaShieldingParser.cpp:77-100`),
stores total/dia/para raw tensors (`h5-reader/src/io/OrcaShieldingParser.cpp:99-110`),
and attaches `orca_coord` by nucleus index / coordinate order
(`h5-reader/src/io/OrcaShieldingParser.cpp:119-122`). The model type documents
that this index is intended to be ORCA nucleus index == emitted-PDB atom order ==
topology order (`h5-reader/src/model/DftShielding.h:47-49`).

The Kabsch check exists, but it is diagnostic only. It compares ORCA input coords
against H5 positions for each DFT row, computes optimal rotation/RMSD, and records
mean/max angle/RMSD (`h5-reader/src/rediscover/ExtractionSupport.cpp:101-162`).
It does not reject on a threshold.

What I would change in the plan: do not downgrade the element check to "only
belt-and-suspenders" yet. The reader loader validates atom count, parser holes,
and total == dia + para (`h5-reader/src/io/DftShieldingLoader.cpp:78-105`), but
it does not compare `DftAtomShielding.element` to the topology atom. The
producer/static ORCA path does have that hard element check
(`/shared/2026Thesis/nmr-720/src/OrcaShieldingResult.cpp:181-204`), which is
probably why the assumption feels safe. The reader DFT path should get the same
cheap guard before it is used to justify new data.

Kabsch enforcement can remain diagnostic for the existing 1P9J campaign if the
actual audit stats are near zero, but for future 720-WT or re-prepped trajectory
data I would require a recorded threshold in the substrate build log or coverage
manifest.

## 3. Full-tensor DFT target coverage

### 1P9J

The 1P9J campaign exists locally at `/shared/2026Thesis/1p9j-orcas/`. Shallow
filename counts under `jobs/` found:

- 751 directories with both a `_meta.json` and a `_nmr.out`.
- 751 `.pdb` files.
- 751 `.xyz` files.
- 761 raw `*_nmr.out` files total, which implies some extra/retry outputs and
  reinforces why the reader's `meta.json -> files.out_primary` rule matters.

This matches the project note that Stage 2 narrowed to 1P9J and consolidated a
751-frame r2SCAN/def2-SVP ORCA campaign to `/shared/2026Thesis/1p9j-orcas/`
(`/shared/2026Thesis/nmr-720/CLAUDE.md:75-79`). I did not parse all 751 outputs
or verify every alignment statistic in this review, so the next step should be a
machine-readable "parsed frames / atoms / elements / alignment" audit, not more
manual inspection.

### 720-WT

The 720 code supports full-tensor static ORCA shielding when the files exist.
The Python catalog has optional `orca_total`, `orca_diamagnetic`, and
`orca_paramagnetic` 9-component arrays with shielding irreps
(`/shared/2026Thesis/nmr-720/python/nmr_extract/_catalog.py:430-436`).
`OrcaShieldingResult` parses `_nmr.out`, requires atom count agreement, verifies
ORCA element order against topology, stores total/dia/para tensors, and writes
`orca_total.npy`, `orca_diamagnetic.npy`, and `orca_paramagnetic.npy`
(`/shared/2026Thesis/nmr-720/src/OrcaShieldingResult.cpp:152-179`,
`/shared/2026Thesis/nmr-720/src/OrcaShieldingResult.cpp:181-204`,
`/shared/2026Thesis/nmr-720/src/OrcaShieldingResult.cpp:216-252`). The mutation
delta path preserves total/dia/para WT, mutant, and delta shielding tensors
(`/shared/2026Thesis/nmr-720/src/MutationDeltaResult.cpp:458-481`,
`/shared/2026Thesis/nmr-720/src/MutationDeltaResult.cpp:681-704`).

That is capability, not coverage. The checked 720 worktree is on branch
`wt-720-build` and is dirty in reader/extraction files. I found no `.LGS` files
within `maxdepth 6`, no protein subdirectories under
`/shared/2026Thesis/nmr-720/calibration`, and zero `*_WT_nmr.out` or
`*_ALA_nmr.out` files under that calibration tree. The default consolidated
source path used by `calibration/populate.py` is
`/shared/2026Thesis/consolidated`, and that directory is absent in this
environment (`/shared/2026Thesis/nmr-720/calibration/populate.py:35-38`).

More importantly, `populate.py` explicitly says `_nmr.out` is optional for a
"complete" pair: required files are only WT/ALA `.xyz` and `.prmtop`
(`/shared/2026Thesis/nmr-720/calibration/populate.py:39-47`). The README frames
the 720 stage as calibration against DFT WT-ALA deltas, not necessarily absolute
WT full-tensor targets for all proteins
(`/shared/2026Thesis/nmr-720/README.md:11-17`). `CLAUDE.md` also says the
685-protein fleet was stopped and residual-fleet DFT was deferred
(`/shared/2026Thesis/nmr-720/CLAUDE.md:80-83`).

Conclusion: I would not let the plan claim a transferable 720-WT full-tensor DFT
target set based on the current local evidence. The build order should put a
coverage audit first:

1. Enumerate candidate 720 proteins.
2. Count complete WT `_nmr.out` roots.
3. Parse every target with the same loader/convention.
4. Record per-protein, per-element, per-atom counts and parse failures.
5. Only then decide whether the split is leave-proteins-out or 1P9J-only.

Until that exists, the honest v1 is "single-protein conformational
generalization on 1P9J." A later Trp-cage DFT run, if spent, should be treated
as an external sanity check, not hyperparameter-tuning data.

## 4. Equivariance and parity

The plan's target parity is physically right: shielding is a rank-2 response from
axial magnetic field to axial induced field, so the tensor is parity-even, with
irreducible pieces `0e + 1e + 2e`. The code agrees at the target-decomposition
level: T1 is explicitly called an antisymmetric pseudovector
(`h5-reader/src/rediscover/SphericalBasis.cpp:13-17`). The active Python catalog
also documents SphericalTensor irreps as `0e + 1e + 2e`
(`python/nmr_extract/_catalog.py:77-80`) and sets `_SHIELD_IRREPS = "0e + 1e + 2e"`
(`python/nmr_extract/_catalog.py:130-133`). The all-atom source code has a
useful parity contract warning that only `disp_*` is a genuine polar `1o`, while
ring normals and bond axes must not be consumed as polar vectors
(`h5-reader/src/rediscover/AllAtomEquivariant.h:7-23`).

There is a migration hazard. Older docs/comments still show or imply `1o` for
shielding T1 (`python/doc/index.rst:17-19`,
`h5-reader/src/model/QtTimeSeriesBuffers.h:66-70`,
`h5-reader/src/rediscover/CODEX_BRIEF_PIECE3B_2026-06-03.md:57-63`). Some vector
metadata also mixes `irreps="1e"` with `parity="odd"` for fields in the catalog.
The MatTen adapter should not infer target parity from stale docs or generic
Vec3 metadata. It should declare the shielding output irreps explicitly and have
rotation plus reflection tests.

Reflection test detail: for inversion, an axial T1 pseudovector should not flip.
For a mirror/improper rotation, it transforms as `det(R) * R * v`; for example
reflecting x -> -x maps an axial vector with `diag(1,-1,-1)`. The test should
encode that exact rule, not just "pseudo-components flip."

The lmax statement needs one refinement. Output `l <= 2` is correct because a
3x3 tensor decomposes into `0/1/2`. But hidden/model `lmax = 2` is not guaranteed
by physics alone. The response tensor is rank 2; the angular complexity of a
many-atom protein environment is not necessarily exhausted by hidden `lmax=2`.
Using `lmax=2` as the default is defensible and cheaper, but I would phrase
higher intermediate irreps as a sensitivity/compute tradeoff, not as physically
irrelevant.

## 5. Scientific requirements

What is ready:

- Geometry-only v1 is the right thesis claim. Feeding the hand-built kernels into
  MatTen would blur the result; using ridge-on-kernels as a baseline is fine
  because it is a comparator, not model input.
- The split logic is directionally right: temporal/block holdout for 1P9J,
  leave-proteins-out only if a real multi-protein full-tensor DFT corpus exists.
- Per-irrep and per-element scaling are not over-engineering. Without them, T0
  will dominate and the T2/T1 story becomes hard to defend.
- Cutoff sensitivity is scientifically necessary because NMR shielding has both
  short-range electronic and longer-range ring-current contributions.
- Reporting seed/fold variance is necessary for a small, autocorrelated DFT set.

What I would add before calling the first number defensible:

1. A target coverage manifest: parseable frames, atom counts, element counts,
   NaN/gap counts, and alignment stats. This should be produced before training.
2. A hard split artifact committed beside the substrate. For 1P9J, use temporal
   blocks with gaps; random per-row splits are invalid.
3. A same-trajectory autocorrelation baseline, such as element/atom-name mean
   plus "nearby frame" or block-aware baselines, so the MatTen number is not just
   exploiting adjacent-frame similarity.
4. Metrics by irrep and element: T0, T1 vector/norm, T2 component/norm, raw 3x3
   Cartesian MAE/RMSE/R2, and full-tensor reconstruction error after inverse
   transform.
5. Round-trip tests: raw 3x3 -> `0e+1e+2e` -> raw 3x3, rotation equivariance,
   reflection parity, and row-grouped molecule reconstruction from the emitted
   substrate.
6. One compact cutoff/depth grid, not a large architecture sweep. The DFT set is
   the limiting resource; the first result should be robust, not exhaustively
   optimized.

## Recommended plan edits

1. Keep the two reader emits in scope. They are additive and well justified.
2. Add one sentence to the positions emit: "Rows must be regrouped by `h5_row` or
   `original_index` and sorted by `atom_index` to create one MatTen molecule per
   DFT frame."
3. Strengthen the reader guardrails: add topology-vs-ORCA element validation and
   record Kabsch thresholds in the substrate audit before future campaigns.
4. Move the 720-WT transferability claim behind an explicit coverage audit. The
   current local evidence supports 1P9J; it does not prove a 720-WT full-tensor
   DFT target set exists.
5. Rephrase hidden `lmax=2` as a default/sensitivity choice, while keeping output
   irreps fixed at `0e + 1e + 2e`.
6. Make the reflection test precise for axial vectors: `v' = det(R) R v`.

Bottom line: build the emits, but do not spend Jessica's next DFT/modeling budget
on a transferable 720 story until target coverage is pinned. The ready, honest
v1 is a 1P9J full-tensor conformational-generalization pilot with strong
baselines and explicit equivariance/parity tests.
