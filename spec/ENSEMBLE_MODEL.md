# Ensemble Extraction Model

**Trajectory-scope object model landed.** The
`TrajectoryProtein` + `TrajectoryAtom` + `TrajectoryResult` +
`Trajectory` shape is documented in `OBJECT_MODEL.md`
("Trajectory-scope entities" section) and `PATTERNS.md` §§13-18;
the old `GromacsProtein` / `GromacsProteinAtom` / `AnalysisWriter`
classes are in `learn/bones/`. Sections of THIS document that
describe the old object-model shape (`GromacsProtein`,
`GromacsProteinAtom`, two-pass scan/extract ownership) are
historical. Non-object-model content here — EDR handling, streaming
mechanics, convergence discussion, test coverage — remains relevant.
For design working-notes see
`spec/pending_include_trajectory_scope_2026-04-22.md`. When the
WIP folds in, this document will be reconciled.

**Status:** IMPLEMENTED and TESTED. Two-pass streaming architecture
works end-to-end. Three tests pass on 1ZR7_6721 full-system data
(479 protein atoms, 3525 water molecules, 25 ions). .h5 master file
writing is IMPLEMENTED (2026-04-13, HighFive v2.6.2).
WriteH5 produces positions (T,N,3), 48-column Welford rollup,
per-bond stats, frame times. SDK `load_trajectory()` reads it.

Design revised 2026-04-11/12 after critical review of original
EnsembleConformation/observer design. See memory:
`project_ensemble_design_decisions.md`.

---

## Architecture: GromacsProtein pattern

No abstract EnsembleConformation. No observer/accumulator ABC.
Four concrete classes married to GROMACS:

### GromacsProtein

Wraps a real Protein. NOT a feature extractor.

Two build paths:

- **Build(FleetPaths)**: pre-extracted PDB poses. Protein is fully
  constructed (bonds, rings, conformations) immediately. InitAtoms()
  called at the end.

- **BuildFromTrajectory(dir_path)**: dir_path must contain md.tpr,
  md.xtc, md.edr. FullSystemReader parses the TPR once — atom ranges,
  bonded parameters, and protein construction all from one parse.
  EDR preloaded for O(1) per-frame lookup. Does NOT finalize the
  protein -- FinalizeProtein() must be called after the first XTC
  frame provides real coordinates for bond detection.

**Deferred finalization pattern (trajectory path):**

1. `BuildFromTrajectory` reads TPR: atom count known, topology known,
   charges loaded. Protein has atoms but NO bonds, NO rings, NO
   conformations.
2. `GromacsFrameHandler::Open` reads first XTC frame, PBC-fixes
   protein coordinates via MoleculeWholer, splits via FullSystemReader.
3. `gp.FinalizeProtein(positions, time_ps)`:
   - `protein_->FinalizeConstruction(positions)` -- bond detection
   - `protein_->AddMDFrame(positions, ...)` -- conformation 0 (permanent)
   - `InitAtoms()` -- creates GromacsProteinAtom Welford accumulators

After finalization, conformation 0 lives in the Protein's
conformations_ vector. It is permanent and feeds the .h5 master file.

Streaming frames are **free-standing ProteinConformations** created
by GromacsFrameHandler via the public constructor with a
`const Protein*` back-pointer. They are NOT added to the Protein's
conformations_ vector. They point at the Protein for topology lookups,
are computed on by all calculators, then freed.

This works because:
- ProteinConformation constructor is public, no side effects
- No calculator accesses conformations through the Protein's vector
- OperationRunner::Run takes ProteinConformation& directly
- Verified: zero hidden coupling in all calculators and loaders

Holds: the real Protein, ChargeSource from TPR, FullSystemReader
(topology + frame splitting), per-atom GromacsProteinAtom accumulators.

**AccumulateFrame(conf, frame_idx)**: called by GromacsFrameHandler
after each frame's calculators have run. Reads ConformationAtom fields
(positions, SASA, water counts, E-field magnitudes, hydration geometry)
into per-atom Welford accumulators. Also tracks dry/exposed frame counts.

**WriteCatalog(path)**: writes atom_catalog.csv with per-atom Welford
summaries (mean, std, min/max with frame indices) for all tracked
quantities.

**SelectFrames(max_frames)**: selects frames for full extraction from
accumulated stats. Returns sorted, deduplicated frame indices.

### FullSystemReader

Reads the full-system GROMACS TPR topology to identify atom ranges
and splits full-system XTC frames into protein + solvent.

**ReadTopology(tpr_path)**: parses the TPR to build a SystemTopology:
- `protein_start`, `protein_count` -- atom range in full-system frame
- `water_O_start`, `water_count` -- water molecules (3 atoms each)
- `ion_start`, `ion_count` -- counter-ions
- `total_atoms`
- Per-water charges (`water_O_charge`, `water_H_charge`)
- Per-ion charges and atomic numbers

**ExtractFrame(full_frame_xyz, protein_positions, solvent)**: given
a full-system XTC frame (all atoms, in nm), extracts:
- `protein_positions` -- protein atom coords converted to Angstroms
- `solvent` -- SolventEnvironment with waters and ions

Owned by GromacsProtein. Borrowed by GromacsFrameHandler via
`gp.sys_reader()` for each frame.

### SolventEnvironment

Per-frame explicit solvent data, built by FullSystemReader::ExtractFrame.

- **WaterMolecule**: O/H1/H2 positions (Angstroms) + O_charge + H_charge
  (elementary charges, from TPR, typically TIP3P: O=-0.834, H=+0.417).
  Has `Dipole()` method returning the molecular dipole vector.

- **Ion**: position (Angstroms) + charge + atomic_number
  (e.g. Na+=11/+1, Cl-=17/-1).

- **SolventEnvironment**: vectors of WaterMolecule + Ion, plus
  `water_O_positions` for spatial indexing. `Empty()`/`WaterCount()`/
  `IonCount()` accessors.

Passed to calculators via `RunOptions.solvent` (pointer, null = no
solvent calculators). WaterFieldResult and HydrationShellResult check
`opts.solvent && !opts.solvent->Empty()` before computing.

### GromacsFrameHandler

Owns the XTC stream and MoleculeWholer (PBC fix). Borrows GromacsProtein
for accumulation, sys_reader, and protein access.

**Open(xtc_path, tpr_path)**: opens XTC for streaming, creates
MoleculeWholer from TPR, verifies protein atom count agreement between
MoleculeWholer and FullSystemReader. Reads first frame to finalize
protein:

1. Read first XTC frame
2. PBC fix on protein portion via `wholer_->make_whole()`
3. Write fixed coords back into full frame
4. Split via `gp.sys_reader().ExtractFrame()` -> protein_pos + solvent
5. `gp.FinalizeProtein(protein_pos, time_ps)`
6. Run all calculators on conformation 0 via `OperationRunner::Run`
7. `gp.AccumulateFrame(conf0, 0)` -- frame 0 into Welford accumulators

**Next(opts, output_dir, accumulate)**: processes the next frame:

1. Read full-system XTC frame (streaming, one at a time)
2. PBC fix on protein coordinates via MoleculeWholer
3. Put fixed coords back, split via FullSystemReader
4. Create free-standing `ProteinConformation(&protein, positions)`
   held as `last_conf_` member (unique_ptr)
5. Run calculators via `OperationRunner::Run` with `opts.solvent = &solvent`
6. If `accumulate=true`: `gp.AccumulateFrame(conf, frame_idx)`
7. If `output_dir` non-empty: write NPY to `output_dir/frame_NNNN/`
8. Conformation persists in `last_conf_` until next Next() or Skip()

**conformation()**: returns const ref to the most recently processed
conformation with all calculator results attached. Valid after Next()
returns true, invalidated by next Next(), Skip(), or Reopen().
This is how AnalysisWriter harvests per-frame data — it calls
`handler.conformation()` between Next() calls.

Returns false at EOF.

**Skip()**: skip one frame without processing. Returns false at EOF.

**Reopen()**: reopen XTC from the start for pass 2. Skips first frame
(already processed during Open). Resets frame_idx to 0.

**Two-pass scan/extract pattern:**

- **Pass 1 (scan)**: `Open()` + `Next(scan_opts)` in a loop. Lightweight
  calculators only (skip MOPAC, APBS, Coulomb, DSSP). Welford
  accumulators collect per-atom statistics. No NPY written.

- **Between passes**: `SelectFrames()` picks winner frame indices from
  accumulated stats. `Reopen()` resets XTC to start.

- **Pass 2 (extract)**: loop through frames. For selected frames:
  `Next(extract_opts, output_dir, /*accumulate=*/false)` runs full
  calculators and writes NPY. For others: `Skip()`. The
  `accumulate=false` flag prevents double-counting frames already
  accumulated during pass 1.

### GromacsFinalResult

Runs after all frames. Currently:
- Writes `atom_catalog.csv` from GromacsProtein's per-atom Welford
  accumulators (RMSF, water_n_first, water_emag, SASA, half_shell,
  dipole_cos, nearest_ion -- each with mean/std/min/max/frame indices,
  plus dry/exposed frame counts).
- Logs what was processed (protein_id, frame count, atom count,
  conformation count).

IMPLEMENTED (2026-04-13): `{protein_id}.h5` master file written via
GromacsProtein::WriteH5 (HighFive v2.6.2). Contains positions (T,N,3),
48-column Welford rollup, per-bond stats, frame times. SDK
`load_trajectory()` reads it. See "Master file" section below.

NOT yet implemented:
- Winner frame selection/movement from temp to output
- Temp directory cleanup

### GromacsProteinAtom

Per-atom trajectory-level Welford accumulators. Private constructor --
only GromacsProtein can create these via InitAtoms(). Atom identity
goes through the Protein back-pointer, NOT duplicated here.

48 named Welford accumulators in `AllWelfords()` (single source of
truth for CSV columns, H5 rollup, and SDK column names). Adding a
Welford = add the field + one line in AllWelfords(); nothing else
changes. Categories (all using Welford online accumulator with
mean/variance/min/max/frame indices):

- Ring current (6): `bs_T0`, `bs_T2mag`, `hm_T0`, `hm_T2mag`,
  `rs_T0`, `rs_T2mag`
- McConnell (4): `mc_T0`, `mc_T2mag`, `mc_aromatic`, `mc_backbone`
- PiQuad + Dispersion (4): `pq_T0`, `pq_T2mag`, `disp_T0`, `disp_T2mag`
- H-bond (2): `hbond_inv_d3`, `hbond_count`
- APBS (3): `apbs_emag`, `apbs_efg_T0`, `apbs_efg_T2mag`
- AIMNet2 (3): `aimnet2_charge`, `aimnet2_efg_T0`, `aimnet2_efg_T2mag`
- SASA (4): `sasa`, `sasa_normal_x`, `sasa_normal_y`, `sasa_normal_z`
- Water (7): `water_n_first`, `water_n_second`, `water_emag`,
  `water_emag_first`, `water_efield_x`, `water_efield_y`, `water_efield_z`
- Hydration (3): `half_shell`, `dipole_cos`, `nearest_ion_dist`
- DSSP (7): `phi_cos`, `psi_cos`, `dssp_hbond_energy`,
  `chi1_cos`, `chi2_cos`, `chi3_cos`, `chi4_cos`
- Bond geometry (1): `mean_bond_angle_cos`
- Deltas (4): `bs_T0_delta`, `aimnet2_charge_delta`, `sasa_delta`,
  `water_n_first_delta` (frame-to-frame DeltaTracker inner Welfords)

Additional per-atom accumulators NOT in AllWelfords():

- Positional dynamics: `position_x`/`y`/`z` for RMSF computation
  (separate because WriteH5 writes positions as (T,N,3) not as rollup)
- Frame counters: `n_frames_dry` (water_n_first == 0),
  `n_frames_exposed` (water_n_first >= 4)
- TransitionCounters: `chi_transitions[4]`, `ss8_transitions`
- `n_bonded_neighbors` (constant across frames, from bond graph)

`RMSF()` returns sqrt(var_x + var_y + var_z).

---

## Three lifetimes (unchanged)

| Object | Lifetime | Holds |
|---|---|---|
| Protein / Atom | Eternal | Identity, topology, bonds, rings |
| ProteinConformation / ConformationAtom | One frame | Positions, computed fields |
| GromacsProtein accumulation state | Across frames | Welford means, histograms, frame indices |

The third lifetime is NOT a separate object model (no EnsembleAtom,
no EnsembleConformation). It's internal state on GromacsProtein,
accessed by GromacsFrameHandler methods.

---

## Conformation 0

0 is 0. First GROMACS frame. No Boltzmann minimum selection.

In the trajectory path, conformation 0 is created during
`FinalizeProtein()` after the first XTC frame provides PBC-fixed
positions. Bond detection runs on these positions. Lives in the
Protein's conformations_ vector via AddMDFrame. Always valid,
always in memory. EnrichmentResult from conformation 0 feeds the
.h5 master file (atom_role, hybridisation, etc.).

For `--pdb` or `--mutant` modes: everything works exactly as before.
The existing path is untouched.

---

## Master file: {protein_id}.h5

HDF5 file containing trajectory-level data. IMPLEMENTED 2026-04-13
via GromacsProtein::WriteH5 (HighFive v2.6.2). The current file
contains positions, Welford rollup, and bond stats. The topology
section below is the original plan for conformation-0 data — not
yet written alongside the trajectory rollup.

```
{protein_id}.h5
├── topology/          <- from Protein object
│   ├── element          (N,) int32
│   ├── residue_type     (N,) int32
│   ├── atom_to_residue  (N,) int32
│   ├── bonds            (B, 2) int32
│   ├── bond_type        (B,) int32
│   ├── ring_atoms       (P, 2) int32
│   └── ring_type        (R,) int32
├── classification/    <- from EnrichmentResult on conformation 0
│   ├── atom_role        (N,) int32
│   ├── hybridisation    (N,) int32
│   ├── graph_dist_ring  (N,) int32
│   └── is_conjugated    (N,) int8
├── graph/             <- from SpatialIndexResult on conformation 0
│   ├── radius_edges     (E, 2) int32
│   └── radius_distances (E,) float64
└── attrs: n_atoms, n_residues, n_rings, pdb_id
```

The SDK (nmr_extract) reads both .h5 and NPY files, presents one
typed object. One SDK, one load() call.

---

## Solvent calculators

WaterFieldResult and HydrationShellResult are BUILT and produce
NPY files per frame when run via OperationRunner:

**WaterFieldResult** (5 NPY files):
- `water_efield.npy` -- (N, 3) per-atom E-field vector from all waters
- `water_efg.npy` -- (N, 9) per-atom EFG tensor from all waters
- `water_efg_first.npy` -- (N, 9) EFG from first-shell waters only
- `water_efield_first.npy` -- (N, 3) E-field from first-shell only
- `water_shell_counts.npy` -- (N, 2) first/second shell water counts

**HydrationShellResult** (1 NPY file):
- `hydration_shell.npy` -- (N, 4) half-shell asymmetry, dipole cos,
  nearest ion distance, and shell density

SolventEnvironment reaches these calculators via
`RunOptions.solvent` pointer, set by GromacsFrameHandler in
ProcessFrame. OperationRunner checks
`opts.solvent && !opts.solvent->Empty()` before dispatching.

**Not yet integrated into the accumulator pattern:**

- **HydrationGeometryResult** -- BUILT. Produces per-atom hydration
  geometry features (coordination number, tetrahedral order parameter,
  angular distribution). Needs wiring into AccumulateFrame +
  corresponding Welfords on GromacsProteinAtom.

- **EeqResult** -- BUILT. Electronegativity-equilibration charges.
  Needs wiring into AccumulateFrame + corresponding Welfords.

These are the next integration items. See
`spec/OUTSTANDING_GROMACS_PATH.md` for the full task list.

---

## Charge sensitivity

Per-atom charge sensitivity is NOT a calculator feature. It is NOT
on Atom (topology) and NOT on ConformationAtom.

The perturbation approach (10 random bulk displacements, Welford
variance) was removed. It produced conformation-specific values
via random splats that were non-comparable across frames and lied
about solvation when applied to other conformations.

Two quantities replace it, both computed by GromacsFrameHandler:

1. **Ensemble charge variance** -- Welford on aimnet2_charge across
   all frames. Free (charges already computed per frame). Output:
   `ensemble_charges_var.npy`.

2. **Autograd charge sensitivity** -- d(charges)/d(positions) via
   backward pass through the .jpt model. Deterministic, exact,
   per-conformation. Validated: the model supports autograd on
   coords. Cost: N backward passes per frame (~5s for 5000 atoms).
   Run on selected frames only (conformation 0, survivors).

---

## What does NOT change

- Protein / Atom: untouched
- ConformationAtom: untouched (one frame's data)
- ConformationResult: untouched (computes from one frame)
- OperationRunner::Run: untouched (sequences one frame)
- WriteFeatures per calculator: untouched
- Python SDK: loads whatever NPY files are present
- Existing `--pdb`, `--mutant`, `--orca` modes: untouched

---

## Superseded designs (historical)

The following concepts from the original 2026-04-10 design were
evaluated and rejected during the 2026-04-11/12 review:

- **EnsembleConformation / EnsembleAtom** -- replaced by accumulation
  state on GromacsProtein. No separate object model needed.
- **EnsembleResult / observer ABC** -- replaced by direct methods on
  GromacsFrameHandler. No virtual dispatch.
- **TrajectoryObserver / EvaluationAccumulator** -- same.
- **ConformationList with invalidation** -- replaced by free-standing
  conformations that die naturally after processing.
- **GraphExportResult** -- bond graph goes in .h5 from Protein topology.
- **charge_sensitivity on Atom** -- moved to evaluator/handler output.
  Perturbation approach deleted.
- **PDB-based protein building for trajectory** -- BuildProteinFromTpr
  replaced PDB builder after atom count crash (commit 0b829b7).

---

## Test data

1ZR7_6721 (479 protein atoms, 3525 water molecules, 25 ions) at
`tests/data/fleet_test_fullsys/1ZR7_6721/` with 5 walkers, XTC
trajectories, run_params (prod_fullsys.tpr).

Three tests in `test_gromacs_streaming.cpp`:

1. **FleetBuildAndCompute** -- fleet-path Build, run calculators
   on conformation 0.

2. **TrajectoryBuildAndScan** -- BuildFromTrajectory + Open +
   scan 10 frames with lightweight calculators. Verifies deferred
   finalization (bonds=0 before Open, bonds>0 after). Checks
   accumulator frame counts match. Writes atom_catalog.csv.

3. **TrajectoryExtractSelectedFrames** -- full two-pass pattern.
   Scan 20 frames, Reopen, extract frames 3 and 7 with full
   calculators and NPY output. Verifies frame_0003/ and frame_0007/
   directories contain NPY files.

1Q8K_10023 (4876 atoms, 26 rings) at `tests/data/fleet_test_large/`
for fleet-path testing with 5 walkers and initial-samples PDB poses.

The fleet tree is at `/shared/fleet/clean_fleet/results/`.
Do not modify it.
