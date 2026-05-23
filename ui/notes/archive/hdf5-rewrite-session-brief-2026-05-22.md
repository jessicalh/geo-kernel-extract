# nmr-viewer HDF5 reader rewrite — session brief

**Written:** 2026-05-22, after commit `353aad4` (viewer re-alignment + audit).
**Status:** Cut-and-paste prompt for the *next* session. The user will copy
the body below into a fresh Claude conversation. Pre-Phase-1 of the work.
**Author note:** Grounded in actual file probes done in the originating
session — schema dump done via h5py against a real production fixture,
producer side read at file:line, sibling subproject (`h5-reader/`) state
confirmed by reading its includes. A first draft of this brief was written
from grep counts and was structurally wrong about the format situation;
this revision corrects that.

---

You're picking up the nmr-viewer at `/shared/2026Thesis/nmr-shielding/`.
Previous session (commit `353aad4`) re-aligned the viewer with the current
library, wired multicast UDP logging, added inspector sections for every
per-atom calculator the live pipeline produces, and stood up
`ui/tests/test_inspector.py` via two new REST endpoints (`atom_dump`,
`list_atoms`). That work is complete and validated.

**THE JOB FOR THIS SESSION:** the viewer's `--analysis-h5 PATH` reader is
broken at the format level, not just the field level. The old AnalysisWriter
output that `fileformat/analysis_file.h` was built to read is no longer
produced by anything in the library. The library switched to a
per-TR-emits-its-own-group architecture some time after the topology
substrate work landed. Both Qt projects (this viewer and the `h5-reader/`
subproject) are reading a dead format.

The user describes this as "basically a rewrite." It is.

---

## Ground truth — verified by reading source + h5py-probing a real file

### THE OLD FORMAT (dead)

- `fileformat/analysis_file.h` declared `AnalysisFile` with 22 nested structs
  (`ring_current`, `efg`, `bond_aniso`, `hbond`, `sasa`, `water`, `charges`,
  `aimnet2_embedding`, `bonded_energy`, `dihedrals`, `dssp`, …).
- `AnalysisWriter` in the main library wrote that shape.
- `--trajectory --analysis` was the producer.
- That entire mode is now **STUBBED**. See `src/nmr_extract.cpp:356-366`:
  > `"ERROR: --trajectory --analysis is disabled pending the dissolution`
  > `of AnalysisWriter into per-Result H5 emitters."`
- `fileformat/analysis_file.h` still compiles and can still parse old files.
  The viewer's `ComputeWorker.cpp:260-263` calls `AnalysisFile::ReadH5(spec
  .analysis_h5_path)`. **No current pipeline feeds it a file in that shape.**
- `h5-reader/src/io/QtProteinLoader.h` still includes `analysis_file.h` and
  holds `shared_ptr<const AnalysisFile>`. Same dead reader; same stale
  subproject. Only commit on `h5-reader/src/` is `41b5123` (the original
  landing).

### THE NEW FORMAT (live)

**Producer:** `nmr_extract --trajectory` (no `--analysis`) at
`src/nmr_extract.cpp:340-344`:

```cpp
HighFive::File file(h5_path, HighFive::File::Truncate);
traj.WriteH5(file);    // → src/Trajectory.cpp:376
tp.WriteH5(file);      // → src/TrajectoryProtein.cpp:247
                       //   + each attached TR's WriteH5Group
```

**Verified by h5py on a real production file.** Schema:

```
/atoms/                       (from TrajectoryProtein::WriteH5)
  element            [N]      int32 (Element enum ordinal)
  pdb_atom_name      [N]      object (string)
  residue_index      [N]      uint64
  attrs: protein_id, n_atoms, finalized

/trajectory/                  (from Trajectory::WriteH5)
  attrs: xtc_path, tpr_path, edr_path, configuration

/trajectory/frames/
  time_ps            [T]      float64
  original_index     [T]      uint64
  attrs: n_frames

/trajectory/selections/<mangled-typeid>/  (zero or more groups)
  frame_idx          [R]      uint64
  time_ps            [R]      float64
  reason             [R]      object (string)
  metadata_json      [R]      object (one JSON dict per record)
  attrs: n_records, metadata_encoding

/trajectory/<tr_name>/        (ONE GROUP PER ATTACHED TR — 49 possible)
  Per-TR shape, NOT central-doc'd. Read each TR's WriteH5Group.
```

**THREE representative shapes verified against the May-8 1P9J fixture:**

```
/trajectory/positions/        (PositionsTimeSeriesTrajectoryResult)
  xyz                [N,T,3]  float64
  frame_indices      [T]      uint64
  frame_times        [T]      float64
  attrs: result_name, n_atoms, n_frames, finalized

/trajectory/bs_shielding_time_series/  (BsShieldingTimeSeriesTrajectoryResult)
  xyz                [N,T,9]  float64
  frame_indices, frame_times
  attrs: result_name, n_atoms, n_frames, finalized, irrep_layout,
         normalization, parity, units
  irrep_layout = 'T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2'

/trajectory/bs_welford/       (BsWelfordTrajectoryResult)
  n_frames_per_atom  [N]      uint64
  t0_mean, t0_m2, t0_std, t0_min, t0_max,
      t0_min_frame, t0_max_frame      [N] float64 (or uint64 for *_frame)
  t0_delta_mean, t0_delta_std         [N] float64
  t2mag_{mean,m2,std,min,max,...}     [N] float64
  attrs: result_name, n_frames, finalized, ddof, mean_dt_ps,
         frame_index_range, irrep_layout_t1, irrep_layout_t2, units
  units = 'ppm_T_per_nA'
```

Additional shapes you'll find:

- `bond_length_stats`   (per-bond rollup, NOT per-atom)
- `bs_t0_autocorrelation` (`rho[N, L]`, `lag_times_ps[L]`)
- `dihedral_*`, `dssp8_*`, `ring_pucker`, `jcoupling`, `rmsd_*`
- `aimnet2_charge_*`, `aimnet2_embedding_*`, `aimnet2_charge_response_*`
- `apbs_efield_*`, `apbs_efg_*`
- `water_field_*`, `hydration_shell_*`, `hydration_geometry_*`
- `tripeptide_bb_*`, `tripeptide_neighbor_*`, `larsen_hbond_*`
- `mopac_*`
- `bonded_energy_*`, `gromacs_energy_*`

**TOPOLOGY is NOT in trajectory.h5.** It's emitted by `TopologySidecar`
(`src/TopologySidecar.h`) as NPY+JSON files in a SEPARATE `output_dir`:

```
residues.npy (structured)
bonds.npy    (structured)
rings.npy
ring_membership.npy
extraction_manifest.json
```

And critically: the trajectory mode in `nmr_extract.cpp:330-345` does **NOT**
call `TopologySidecar::WriteFeatures`. Only the per-conformation modes
(Pdb/Orca/Mutant at lines 79/132/187/255) do.

→ For an `--analysis-h5` viewer load, **there is no guaranteed topology
sidecar to read**. Topology comes from the live PDB build.

---

## Concrete fixtures to point at

**A REAL production trajectory.h5 to test against:**

```
/shared/2026Thesis/nmr-shielding/molprobity_runs/1P9J_5801/extract/trajectory.h5
```

- 60 MB, 846 atoms, 751 frames (15 ns @ 20 ps stride)
- May 8 vintage, so partial TR set (only 6 groups under `/trajectory/`):
  `bond_length_stats`, `bs_shielding_time_series`, `bs_t0_autocorrelation`,
  `bs_welford`, `frames`, `positions`
- Source xtc/tpr/edr in the file's `/trajectory` attributes
- **Best starting fixture** because it's REAL output, not synthetic

**To produce a CURRENT-schema trajectory.h5** (with all currently-attached
TRs), launch `nmr_extract` against the bundled inputs:

```
tests/data/fleet_amber/1P9J_5801/prep_run_*/batcave_local_15ns_*/
  production.{tpr,trr,edr}
tests/data/fleet_amber/1Z9B_6577/prep_run_*/batcave_local_15ns_*/
  production.{tpr,trr,edr}
```

Recipe (verify exact flags via `nmr_extract --help`; the trajectory mode is
described at `src/nmr_extract.cpp:330-345`):

```bash
scripts/run_with_cuda_env.sh build/nmr_extract \
    --trajectory tests/data/fleet_amber/1P9J_5801/prep_run_.../batcave_local_15ns_.../  \
    --output /tmp/1P9J_test \
    --aimnet2 data/models/aimnet2_wb97m_0.jpt
# → produces /tmp/1P9J_test/trajectory.h5 with the full current TR set
```

Inspect either file at the shell with:

```bash
python3 -c "import h5py; f=h5py.File('PATH','r'); print(list(f.keys())); print(list(f['trajectory'].keys()))"
```

A separate fixture exists for the per-TR roundtrip test pattern at
`tests/test_amber_streaming.cpp:310-333` — writes a TR's H5 group to `/tmp`
and reads it back. Good shape to mirror in the new viewer reader's unit tests.

---

## Investigative reading order (one pass, no code yet)

### Producer side (where the canonical schema lives)

1. `src/Trajectory.cpp:376-460` (`Trajectory::WriteH5`)
2. `src/TrajectoryProtein.cpp:247-300` (`TrajectoryProtein::WriteH5`)
3. `src/nmr_extract.cpp:330-345` (the only caller of those two)
4. **ONE TR per shape category** — enough to see the pattern:
   - Welford rollup: `src/BsWelfordTrajectoryResult.cpp::WriteH5Group`
   - Tensor TS: `src/BsShieldingTimeSeriesTrajectoryResult.cpp`
   - Scalar TS: `src/SasaTimeSeriesTrajectoryResult.cpp`
   - Per-bond stats: `src/BondLengthStatsTrajectoryResult.cpp`
   - Autocorrelation: `src/BsT0AutocorrelationTrajectoryResult.cpp`
   - Selection bag: `src/ChiRotamerSelectionTrajectoryResult.cpp`
     (writes via `traj.MutableSelections()`; read surface lives in
     `Trajectory::WriteH5` selections loop, NOT in the TR's own `WriteH5Group`)
5. `src/TopologySidecar.{h,cpp}` (the separate NPY+JSON path)

### Consumer side (what's currently there, all dead/stale)

6. `ui/src/ComputeWorker.cpp:248-345` (current Phase 2b H5 path)
7. `ui/src/ComputeWorker.h::AnalysisBinding` (current binding shape)
8. `ui/src/MainWindow.cpp::populateTimeSeries` (current "Time Series (H5)"
   tab; search `ts_` or `Time Series (H5)`)
9. `fileformat/analysis_file.h` (the dead reader; 312 lines, 21 sections)
10. `h5-reader/src/io/QtProteinLoader.h` (sibling subproject, same dead reader)
11. `h5-reader/src/model/QtFrame.h` (the QtFrame slab-accessor pattern —
    not a model to follow because it consumes AnalysisFile, but the
    one-frame-at-a-time accessor style may translate)

### Existing validation infra (extend, don't dismantle)

12. `ui/src/RestServer.cpp::cmdAtomDump` (atom_dump endpoint; the
    atom-shaped JSON consumers expect to see H5 time-series data
    arrive in some form)
13. `ui/tests/test_inspector.py` (currently 30 PASS / 12 NOTE / 0 WARN
    on 1UBQ; will need H5 sections added when reader lands)

---

## Scope decisions to surface before writing code

### 1. Single-frame vs. all-frames vs. rollups-only

The viewer is single-conformation per UI_ROADMAP. Three reasonable subsets
for what the new reader actually surfaces:

- **(a) Frame 0 only** (parity with the old behaviour) — read `/atoms`, then
  read frame-0 slab from each `/trajectory/<ts>/xyz` array and the entire
  `/trajectory/<welford>/` N-shaped arrays. Don't load the T axis of
  time-series.
- **(b) Rollups + frame 0** — load Welford rollups in full; load frame 0
  of time-series. Inspector shows "frame 0 value AND ensemble mean ± std"
  per kernel.
- **(c) Full trajectory load** — load every `(N, T, 9)` array. Inspector
  gains a frame slider. Crosses into trajectory-viewer scope that
  UI_ROADMAP explicitly punts.

**(b) is the most useful per-thesis** ("does my live calculator at this
pose agree with the recorded ensemble mean?") and stays under the single-
conf scope umbrella. Confirm with user.

### 2. AnalysisBinding contract — keep, retrofit, or replace

Today `AnalysisBinding` holds `shared_ptr<const AnalysisFile>` + a `libToH5`
identity map. The `libToH5` piece is still meaningful (the new `/atoms`
group has its own atom axis; could disagree with the protein's). The
`AnalysisFile` pointer is dead. Either:

- replace with a typed `TrajectoryFile` holding `HighFive::File` + N + T
  + per-group lazy accessors
- or scrap the binding type and have `ComputeWorker` stash the
  `HighFive::File` and pull data inline

The first preserves the "identity check at bind, lazy reads later"
discipline. The second is more direct. Confirm.

### 3. Shared reader library with h5-reader, or independent

`h5-reader` uses the same dead format and faces the same rewrite. Two options:

- Build a shared reader library (probably under `fileformat/` since that's
  where the I/O layer was) that both viewer and h5-reader depend on. Per
  `ui/CLAUDE.md` "fileformat is frozen during feature sessions" — **explicitly
  ASK the user whether this is schema-change time.**
- Each subproject implements its own reader. More code, more drift, but
  stays in each subproject's scope rules.

This is a coordination decision, not a technical one.

### 4. TopologySidecar — read it when present, or always ignore?

The viewer builds topology from its live PDB load. The sidecar NPYs would
only be useful for cross-check (does the recorded bond list match my live
one?). Trajectory mode doesn't even write them today. Probably **IGNORE for
now**; revisit if needed.

### 5. atom_dump REST endpoint surface for H5 data

Either:

- grow `atom_dump` with an optional `"h5"` block (single endpoint, conditional
  payload)
- add a separate `atom_h5_dump` endpoint
- leave REST as it is; H5 access stays GUI-only

Test surface considerations: `test_inspector.py` currently has `send_cmd`
one-socket-per-call (fine up to ~5000 atoms per the docstring); per-atom H5
dump adds another full pass over atoms. The N-shaped Welford datasets are
read-once-then-slab; the `(N,T,9)` time-series are larger but you only need
a single frame slab. Cost analysis worth doing before deciding.

### 6. The dead AnalysisFile reader — leave, tear out, or quarantine

The reader still compiles and the standalone `fileformat/` build still
produces `libanalysis_file.a` + `roundtrip_test`. It serves exactly ZERO
live consumers (writer is stubbed, viewer is being rewritten, h5-reader is
stale). Three options:

- **tear out:** remove `fileformat/` from CMakeLists, remove the include
  from ComputeWorker. Cleanest; removes 312 lines of dead-format guarantee.
- **leave:** lets future archaeology read old fleet H5s. Costs only library
  size.
- **quarantine:** move to `fileformat/legacy/`, mark deprecated, update
  README.

This is a user call. Ask before doing.

---

## Deliverable shape

### Phase 1 — investigation + plan (this session)

- A written plan with concrete `file:line` refs, same shape as the audit at
  commit `353aad4` (CRITICAL / MAJOR / MINOR buckets, opinionated, single-
  highest-priority pick).
- Answers to all six scope decisions above, gathered by reading the
  codebase and asking the user.
- A proposed schema for: the new reader's public API, what the Inspector
  tab shows, whether atom_dump grows, whether anything in `fileformat/` moves.
- **The user reviews the plan. No code yet.**

### Phase 2 — execution (only after sign-off)

- Small focused commits. Build clean after each.
- `test_inspector.py` grows H5 coverage when the reader lands; every check
  follows `feedback_log_overages_dont_assert` (report, don't assert; only
  structural invariants ERROR).
- Audit-agent pass before calling it done (mirror the Plan-agent pass we
  did at `353aad4`).

---

## What works, don't dismantle

- `nmr::Session` is built once in `main_viewer.cpp` and held by MainWindow
  + ComputeWorker. AIMNet2 model + tripeptide DFT libpq connection + Larsen
  H-bond grids live there.
- Multicast UDP logging to `239.255.0.1:9998` — both viewer Log dock and
  `udp_listen.py` co-receive. `~/.nmr_tools.toml` carries the multicast
  destination.
- REST endpoints `atom_dump` and `list_atoms`. Inspector / atom_dump /
  test_inspector.py are a triple per `ui/CLAUDE.md`.
- ComputeWorker / workerThread lifecycle uses Qt's canonical
  `QThread::finished → deleteLater` pattern.
- Launch always through `scripts/run_with_cuda_env.sh` (cu13
  `LD_LIBRARY_PATH` dependency).
- Library skip flags: `skip_mopac=true`, `skip_coulomb=true`. Tripeptide
  DFT and Larsen H-bond DO run via Session-loaded backends. Don't change
  these without asking.

---

## Before you speak — sanity checks against drift

- There are 49 `*TrajectoryResult.h` headers in `src/` as of `353aad4`.
  Confirm by `ls src/*TimeSeriesTrajectoryResult.h src/*WelfordTrajectoryResult.h | wc -l`.
  Significantly fewer or more means something shifted; surface.
- The real fixture `trajectory.h5` has 6 TR groups (May 8 vintage, partial
  TR set). A fresh `nmr_extract` run on the same xtc/tpr/edr should produce
  more groups today. If it doesn't, something is not attaching that you'd
  expect to be.
- `nmr_extract --trajectory --analysis` STILL errors out with the stub
  message. Confirm.
- `fileformat/analysis_file.cpp::ReadH5` still compiles. Confirm by
  `cd build && cmake --build . --target analysis_file`.

You are starting fresh. The user trusts you to read carefully and propose
before acting. Begin with Phase 1.
