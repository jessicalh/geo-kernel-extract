# HDF5 rewrite — Phase 1 plan

**Status:** Pre-execution plan. Reviewed scope decisions locked
(2026-05-23 AskUserQuestion). NO CODE YET.

**Locks (from the AskUserQuestion):**

1. **Inspector scope** — Rollups + frame 0. The "live calc vs ensemble"
   diagnostic. Welford `mean ± std` per metric, plus frame-0 slab for
   each `(N,T,9)` / `(N,T,3)` / `(N,T)` time series.
2. **Binding fate** — Replace `shared_ptr<const AnalysisFile>` inside
   `AnalysisBinding` with a typed `TrajectoryH5` object. Identity-check
   role (`libToH5`, name-mismatch logging) preserved.
3. **Sharing** — Independent in `ui/`. No `fileformat/` modifications.
   `h5-reader/` adopts in its own future session.
4. **REST surface** — Grow `cmdAtomDump` with optional `"h5"` block,
   present iff the binding `Valid()`. `test_inspector.py` extends in
   lockstep.

**Pre-locked (not asked):**

- **`TopologySidecar`** — Ignore. Viewer builds topology from live PDB.
- **`fileformat/analysis_file.{h,cpp}`** — Out of scope this session
  (`ui/CLAUDE.md`: fileformat is frozen during feature sessions). We
  drop the include from `ui/` only. Library-window decides leave /
  tear-out / quarantine.

---

## The single highest-priority pick

**Build `ui/src/TrajectoryH5.h` + `.cpp` — the typed read boundary for
the new per-TR-emits-its-own-group format. Everything else in this
plan cascades from it.**

The dead `AnalysisFile` was structurally correct (typed structs, one
place that knew the schema). The new `TrajectoryH5` is the same shape
against the new schema. It is not a wrapper/adapter — it IS the
typed boundary at the H5 ingress, symmetric to `PdbFileReader` at
structure ingress.

---

## Concrete fixture for development + validation

`/shared/2026Thesis/nmr-shielding/molprobity_runs/1P9J_5801/extract/trajectory.h5`

- 60 MB, 846 atoms, 751 frames (15 ns @ 20 ps stride)
- Sparse TR set (May 8 vintage): only 6 groups present
  (`bond_length_stats`, `bs_shielding_time_series`, `bs_t0_autocorrelation`,
  `bs_welford`, `frames`, `positions`)
- Real production output — first fixture to point at.

For full-current-schema fixture, a fresh `nmr_extract --trajectory`
on `tests/data/fleet_amber/1P9J_5801/` would produce a file with the
complete current TR catalog. That run is part of Phase 2 setup, not
Phase 1.

---

## Proposed public API: `TrajectoryH5`

```cpp
// ui/src/TrajectoryH5.h
//
// Typed boundary at the trajectory.h5 read ingress. Eager-loads
// the bounded-cost slices we need (Welford N-shaped + frame-0 slab
// of each TS) at construction; closes the HighFive::File handle.
// All accessors return std::optional<> where the underlying group
// is absent from the file (sparse TR sets are normal).
//
// Throws on construction if the file is missing the structural
// minimum (/atoms, /trajectory/frames). Otherwise sparse-tolerant.

namespace nmr_ui {

struct WelfordPair { double mean = 0, std = 0; };

struct BsWelfordRow {
    WelfordPair t0;
    WelfordPair t2magnitude;
    // T1/T2 per-component pairs available if needed; inspector reads
    // the two headline numbers above.
};
struct HmWelfordRow      { /* identical shape to Bs */ };
struct McWelfordRow      { /* identical shape; units Å⁻³ */ };
struct EeqWelfordRow     { WelfordPair charge; };
struct SasaWelfordRow    { WelfordPair sasa; };
struct HBondCountWelfordRow { WelfordPair count; WelfordPair occupancy_fraction; };

class TrajectoryH5 {
public:
    explicit TrajectoryH5(const std::string& path);    // throws on structural failure

    // Top-level metadata
    size_t  AtomCount()  const { return n_atoms_; }
    size_t  FrameCount() const { return n_frames_; }
    double  FrameTimePs(size_t t) const { return frame_times_[t]; }

    // /atoms identity (for cross-check against the library Protein)
    int          ElementAt(size_t i)        const { return atom_element_[i]; }
    const std::string& AtomNameAt(size_t i) const { return atom_name_[i];    }
    size_t       ResidueIndexAt(size_t i)   const { return residue_index_[i]; }

    // /trajectory attrs
    const std::string& XtcPath() const { return xtc_path_; }
    const std::string& TprPath() const { return tpr_path_; }
    const std::string& EdrPath() const { return edr_path_; }

    // Per-Welford rollups at the atom axis (frame-0-independent).
    // Each returns nullopt iff the group is not in this file.
    std::optional<BsWelfordRow>          BsWelford(size_t atom_idx) const;
    std::optional<HmWelfordRow>          HmWelford(size_t atom_idx) const;
    std::optional<McWelfordRow>          McWelford(size_t atom_idx) const;
    std::optional<EeqWelfordRow>         EeqWelford(size_t atom_idx) const;
    std::optional<SasaWelfordRow>        SasaWelford(size_t atom_idx) const;
    std::optional<HBondCountWelfordRow>  HBondCountWelford(size_t atom_idx) const;

    // Frame-0 slabs of each time series.
    std::optional<SphericalTensor> BsShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> HmShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> McShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> PiQuadShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> RingChiShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> DispShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> HBondShieldingFrame0(size_t atom_idx) const;
    std::optional<double>          SasaFrame0(size_t atom_idx) const;
    std::optional<double>          Aimnet2ChargeFrame0(size_t atom_idx) const;
    std::optional<Vec3>            ApbsEfieldFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> TripeptideBbShieldingFrame0(size_t atom_idx) const;
    std::optional<Vec3>            TripeptideBbResidualVecFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> TripeptideNeighborShieldingFrame0(size_t atom_idx) const;
    std::optional<Vec3>            TripeptideNeighborResidualVecPrevFrame0(size_t atom_idx) const;
    std::optional<Vec3>            TripeptideNeighborResidualVecNextFrame0(size_t atom_idx) const;
    std::optional<uint8_t>         TripeptideBbMethodTagFrame0(size_t atom_idx) const;
    std::optional<double>          LarsenHBondWaterTermFrame0(size_t atom_idx) const;
    std::optional<int32_t>         LarsenHBondCountFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> LarsenHBond1pHBShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> LarsenHBond2pHBShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> LarsenHBond1pHaBShieldingFrame0(size_t atom_idx) const;
    std::optional<SphericalTensor> LarsenHBond2pHaBShieldingFrame0(size_t atom_idx) const;
    std::optional<Vec3>            PositionsFrame0(size_t atom_idx) const;

    // Sparse-set introspection: what groups did this file actually have?
    bool HasGroup(const std::string& name) const;
    const std::vector<std::string>& GroupsPresent() const { return groups_present_; }

private:
    // Eager-loaded data. All bounded: at 1000 atoms, 24 TS groups,
    // 6 Welford groups → total ~5 MB.
    size_t                              n_atoms_ = 0;
    size_t                              n_frames_ = 0;
    std::vector<double>                 frame_times_;
    std::vector<int>                    atom_element_;
    std::vector<std::string>            atom_name_;
    std::vector<size_t>                 residue_index_;
    std::string                         xtc_path_, tpr_path_, edr_path_;
    std::vector<std::string>            groups_present_;

    // Per-Welford eager buffers, parallel to atom axis (N-sized vectors).
    std::optional<std::vector<BsWelfordRow>>         bs_welford_;
    // ... one per Welford group ...

    // Per-TS frame-0 slabs (atom-axis, N-sized vectors).
    std::optional<std::vector<SphericalTensor>>      bs_shielding_f0_;
    // ... one per TS group ...
};

}  // namespace nmr_ui
```

**Key shape choices:**

- **Eager load at construction.** Bounded cost; closes the HighFive
  handle immediately. No lifetime / threading complications later.
- **`std::optional<>` per group.** The May-8 fixture has 6 of 24 TS
  groups present. Sparse-tolerance is the default, not an error case.
- **One typed accessor per TR group.** No generic `Get(name)` adapter
  surface. Each TR's schema is known at compile time on the consumer
  side. Adding a new TR = adding one accessor + one private member.
  Same shape as `ConformationAtom` field growth.
- **No `HighFive::File` member.** The class is a snapshot, not a
  handle. The handle lives only during the constructor.

---

## CRITICAL bucket — must land for the rewrite to work

| # | Item | Files | Notes |
|---|------|-------|-------|
| C1 | New `ui/src/TrajectoryH5.{h,cpp}` reader | NEW | The typed boundary at H5 ingress. ~24 accessor methods, one per TR group attached today. |
| C2 | `AnalysisBinding` swap — replace `shared_ptr<const AnalysisFile> h5` with `shared_ptr<const TrajectoryH5> h5`; drop `#include "analysis_file.h"` | `ui/src/ComputeWorker.h:25,71` | Identity-check role stays (`libToH5`, `nameMismatches`). |
| C3 | `ComputeWorker.cpp` Phase 2b — rewrite to `TrajectoryH5(path)` (throw → log + skip), atom-count + per-index element check against `protein.AtomAt(i).element` | `ui/src/ComputeWorker.cpp:240-347` | Same identity-check shape; new accessor names. Drop `<stdexcept>` boundary comment about `AnalysisFile`. |
| C4 | `MainWindow::populateTimeSeries` — rewrite from `AnalysisFile`'s 22 nested structs to `TrajectoryH5`'s typed accessors | `ui/src/MainWindow.cpp:2009-2178` | Inspector reads `mean ± std` row + `frame 0` row per metric. Each `optional<>` decides whether a section appears. |
| C5 | Build wiring — add `ui/src/TrajectoryH5.cpp` to viewer target; remove `analysis_file` from viewer link list (verify standalone `fileformat/` target untouched) | `CMakeLists.txt` (viewer target only — bottom of file) | Per `ui/CLAUDE.md`: viewer target only, no library-target edits. |

---

## MAJOR bucket — clear value, not blocking

| # | Item | Files | Notes |
|---|------|-------|-------|
| M1 | `cmdAtomDump` grows optional `"h5"` block when `analysisBinding_.Valid()`. Mirrors the inspector's section shape (welford rollups + frame-0 slabs). | `ui/src/RestServer.cpp::cmdAtomDump` | Per `ui/CLAUDE.md` "Inspector / atom_dump / test_inspector are a triple." Section gated on `Valid()`. |
| M2 | `test_inspector.py` grows H5 coverage section. Per `feedback_log_overages_dont_assert`: report, never assert. Structural invariants only ERROR. | `ui/tests/test_inspector.py` | Drives the M1 endpoint. Verifies sparse-set tolerance (the May-8 fixture exercises this). |
| M3 | `MainWindow.h` — drop `#include "analysis_file.h"`; forward-declare `TrajectoryH5` or include it. | `ui/src/MainWindow.h:10` | Cosmetic but it's the obvious dead-include cleanup once C2/C3 land. |
| M4 | Inspector display polish — show `frame 0` line + `mean ± std` line per metric, side by side. Two `QTreeWidgetItem` children per metric. | `MainWindow.cpp::populateTimeSeries` | This is the visual realisation of "rollups + frame 0". |
| M5 | UDP log on H5 read: log groups present + groups absent (sparse-set diagnostic). One line at bind. | `ComputeWorker.cpp` Phase 2b | Mirrors the existing nameMismatches logging. |

---

## MINOR bucket — polish / good hygiene

| # | Item | Files | Notes |
|---|------|-------|-------|
| m1 | Maybe rename `AnalysisBinding` → `TrajectoryBinding` (the H5 it holds is no longer "analysis"). | `ui/src/ComputeWorker.h`, all references | Cosmetic. Defer if churn cost exceeds value. |
| m2 | `tsSphItem` / `tsVec3Item` / `tsScalarItem` helpers — these are still useful for the new inspector rows. Audit, possibly tweak signatures (no flat vector, just typed values). | `MainWindow.cpp` | Most of the diff at C4 is replacing the data source, not the rendering helpers. |
| m3 | `--analysis-h5 PATH` CLI flag rename to `--trajectory-h5 PATH`? Probably no — keep the old flag for muscle memory, mention the format change in `--help`. | `src/JobSpec.h:77` (read-only — comment only) | Documentation-side; library is not modified. |
| m4 | Drop `AnalysisFile::AIM_DIMS` ref at `MainWindow.cpp:2157`; replace with a constant from `TrajectoryH5` if the embedding TS is still emitted, else delete the section. | `MainWindow.cpp` | Captured at C4 anyway. |

---

## Inspector tab content (the visual realisation)

Each section appears in the QTreeWidget IFF `TrajectoryH5::<group>()`
returns a value at this atom. Two rows per metric: `frame 0` and
`mean ± std`.

```
H5
  protein_id        1P9J
  frame             0 of 751
  frame time (ps)   0.000
  n_atoms           846   n_frames  751
  library atom idx  42    H5 atom idx  42
  atom name (H5)    CA    element (Z)  6
  groups present    bs_welford, bs_shielding_time_series, ...
  atom-name mismatches (total)  0

Ring Current (BS)
  T0 frame 0        +1.42 ppm
  T0 mean ± std     +1.51 ± 0.08 ppm
  |T2| frame 0      0.046
  |T2| mean ± std   0.052 ± 0.004

Ring Current (HM)        (only if hm_welford or hm_shielding_time_series present)
  ...

McConnell                (only if mc_welford or mc_shielding_time_series present)
  ...

SASA                     (only if sasa_welford or sasa_time_series present)
  frame 0           23.4 Å²
  mean ± std        24.1 ± 0.6 Å²

H-bond count             (only if hbond_count_welford present)
  ...

AIMNet2 charge           (only if aimnet2_charge_time_series present)
  ...

APBS efield              (only if apbs_efield_time_series present)
  ...

Tripeptide / Larsen      (only if those TS groups present)
  ...

Positions                (only if positions present)
  frame 0  (X, Y, Z)
```

Each section is its own `if (auto x = h5.BsWelford(i)) {...}` block —
no central registry, no string dispatch. Adding a new section after
a new TR lands is one `if` block.

---

## `atom_dump` REST growth

`RestServer.cpp::cmdAtomDump` already has 14 typed sections. The
`"h5"` block becomes section 15, conditional on
`mw_->analysisBinding_.Valid()`:

```json
{
  "atom_index": 42,
  "identity": {...},
  "amber_substrate": {...},
  ...,
  "h5": {
    "groups_present": ["bs_welford", "bs_shielding_time_series", ...],
    "frame_times": { "n_frames": 751, "first": 0.0, "last": 15000.0 },
    "bs_welford": { "t0_mean": 1.51, "t0_std": 0.08, ... },
    "bs_shielding_frame_0": { "T0": 1.42, "T2": [..., ..., ..., ..., ...] },
    ...
  }
}
```

Sparse-set tolerant: absent groups absent from JSON, not present
with null. `test_inspector.py` checks structural invariants
(`groups_present` is a list of strings, every key in `h5.*` matches
a group in `groups_present`), reports any out-of-range numbers as
NOTE, asserts nothing beyond structural shape.

---

## Phase 2 execution order

1. **C1** — `TrajectoryH5.h` + `.cpp` with the BsWelford + BsShielding
   + Sasa + frame metadata accessors only. Build clean. Write a small
   smoke test (one fixture, one atom_idx, expected mean / frame-0 values).
2. **C2 + C3 + C5** — Swap the binding + ComputeWorker Phase 2b. Build
   clean. Verify the May-8 1P9J fixture loads, identity-check passes,
   nameMismatches still logs, UDP log shows groups present.
3. **C4 + M4** — Inspector tab rewrite. Launch viewer, double-click
   atoms, verify display.
4. **C1 expansion** — add the remaining TR accessors (HM, McConnell,
   PiQuad, RingChi, Disp, HBond, AIMNet2 charge, APBS efield,
   Tripeptide, Larsen). One commit per logical TR family.
5. **M1 + M2** — `atom_dump` H5 block + `test_inspector.py` extension.
   Targeting 30+ NEW pass cases on the May-8 fixture; 0 ERROR.
6. **M3 + M5** — Dead-include cleanup + sparse-set UDP logging.
7. **Audit-agent pass** — Mirror the `353aad4` audit shape. Identify
   any MAJORs to act on before close-out.

---

## What works, don't dismantle (carryover from brief)

- `nmr::Session` lifecycle (process resources)
- Multicast UDP logging
- REST endpoints `atom_dump` and `list_atoms`
- ComputeWorker / workerThread Qt-canonical lifecycle
- `scripts/run_with_cuda_env.sh` launch discipline
- Library skip flags (`skip_mopac=true`, `skip_coulomb=true`)

---

## Open questions / decisions deferred to Phase 2

None blocking. The four locked decisions and two pre-locks are
sufficient to write the code. Surface any new design surface that
emerges during Phase 2 to the user before committing.

---

## Out of scope (explicit)

- Modifying `fileformat/analysis_file.{h,cpp}` (frozen during feature
  sessions per `ui/CLAUDE.md`).
- Modifying `src/` (library off-limits per `ui/CLAUDE.md`).
- Modifying `h5-reader/` (separate subproject; future session).
- Trajectory frame slider in the inspector (`UI_ROADMAP` punts this).
- `TopologySidecar` consumption (library builds topology from PDB).
- New extraction runs (viewer is strictly read-only w.r.t. H5).
