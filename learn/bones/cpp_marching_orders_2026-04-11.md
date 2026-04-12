# C++ Marching Orders: Ensemble Extraction

**PARTIALLY SUPERSEDED (2026-04-12).** This document was the starting
point for the ensemble design. A critical review session (2026-04-11/12)
revised several major decisions:

- **charge_sensitivity**: moved OFF Atom, perturbation DELETED.
  Replaced by ensemble charge variance + autograd on selected frames.
- **GraphExportResult**: REMOVED. Bond graph goes in {protein_id}.h5.
- **TrajectoryObserver / EvaluationAccumulator ABC**: REPLACED by
  GromacsProtein + GromacsFrameHandler + GromacsFinalResult (concrete,
  no ABC).
- **ConformationList / invalidation**: REMOVED. Free-standing
  conformations via public ProteinConformation constructor.
- **Write-all-to-disk strategy**: REPLACED by accumulate-in-memory,
  re-load winners from XTC at the end.

**Current truth:** `spec/ENSEMBLE_MODEL.md` (revised) and
`learn/docs/next_session_2026-04-12.md`.

**What is still valid below:** Change 1 (AIMNet2Result — DONE),
Change 2 (DsspResult extension — DONE, SasaResult — DONE),
TOML/GeometryChoice discipline, calculator parameter specs.
Change 3 (streaming architecture) is superseded.

**Date:** 2026-04-11
**Context:** Session 1 of ensemble accumulation design. These are the
C++ changes needed to produce the GNN training package from fleet MD
trajectories (685 proteins × ~500 subsampled frames).

---

## INVIOLABLE: Configuration and GeometryChoice discipline

**Every cutoff, radius, probe size, neighbour count, or distance
threshold in every new calculator MUST:**

1. **Be a named constant in `calculator_params.toml`**, not a
   hardcoded literal.  No magic numbers.  If it has a unit, the
   TOML key documents the unit.

2. **Have a `GeometryChoice` entry** recorded on the
   `ProteinConformation` after use, identical to the pattern used
   by BiotSavartResult, McConnellResult, CoulombResult, and every
   other existing calculator.  The GeometryChoice records what
   parameter was used and why.

This applies to:
- AIMNet2 short-range cutoff (read from .jpt model, but the
  value is recorded as a GeometryChoice)
- AIMNet2 long-range cutoff (DSF, default 15Å, TOML-configurable)
- AIMNet2 max_nb and max_nb_lr (TOML-configurable)
- SASA probe radius (default 1.4Å, TOML)
- SASA vdW radii table source (Bondi vs ff14SB, TOML)
- SASA n_points (default 92, TOML)
- SASA neighbour cutoff (derived from max_vdw + probe, but the
  max_vdw is from the table which is in TOML)
- Radius graph cutoff (default 10Å, TOML)
- Ring contribution cutoff (already in TOML, unchanged)
- Ensemble frame subsampling stride (TOML)

If a future session agent adds a calculator with a hardcoded
distance, that is a bug.  The TOML is the single source of truth.
The GeometryChoice is the audit trail.

---

## Architecture: calculators first, loop later

**DO NOT IMPLEMENT A STREAMING LOOP OR ENSEMBLE ACCUMULATORS.**

Each new calculator (AIMNet2Result, SasaResult, ChiAngleResult) is a
ConformationResult that runs on any single ProteinConformation, just
like CoulombResult or BiotSavartResult.  Same pattern: Dependencies,
Compute, WriteFeatures.  The ensemble loop and Welford accumulators
are a SEPARATE FUTURE PROJECT that will be designed after all
calculators are working and validated individually.

The streaming loop pseudocode below is CONTEXT for why these
calculators are being written.  It is NOT a spec to implement now.

```
// FUTURE — not current scope
// for each XTC frame:
//     conf = make_conformation(protein, frame_coords)
//     OperationRunner::Run(conf, geometry_only_opts)
//     AIMNet2Result::Compute(conf)
//     SasaResult::Compute(conf)
//     ChiAngleResult::Compute(conf)
//     observers...
```

---

## Change 1: AIMNet2Result

**What:** Neural network charge calculator using libtorch + TorchScript.
Produces per-atom Hirshfeld charges, 256-dim aim embedding, and
AIMNet2 Coulomb EFG tensor.  A ConformationResult like any other.

**Design decisions (settled 2026-04-11 session 2):**

- **Married to AIMNet2.**  No NeuralChargeProvider interface, no
  factory, no pluggable abstraction.  Direct AIMNet2-specific code.
  Named fields: `aimnet2_*` everywhere.  If a better model arrives,
  write a new calculator — a plugin swap hides the differences.
- **CUDA mandatory.**  libtorch CUDA build, no CPU fallback.  RTX 5090
  is the deployment target.
- **Model loaded once** as a static resource (like the ff14SB parameter
  file).  Not passed per-Compute.
- **charge_sensitivity on Atom** (Protein level, permanent).  Computed
  on first AIMNet2Result::Compute via one `if` check.  10 bulk
  perturbations, ~2s.  Not a framework, not a lazy-init abstraction.
  An intrinsic property of the atom's chemical identity, like
  hybridisation — computed from geometry but stored as topology.
- **Geometry-only path.**  Runs at 0.17s/frame alongside BS/HM/MC/PQ
  (~1.3s total).  NOT the MOPAC path.  No MOPAC dependency.
- **PBC infrastructure.**  Port xdrfile, xtc_reader.h, pbc_whole.h
  from fes-sampler VERBATIM.  The MoleculeWholer PBC fix was
  validated on 685 proteins.  Do not redesign.

**Why:** Recovers the carbon EFG dimension lost without MOPAC.
ff14SB EFG direction cosine with MOPAC = 0.709 for C (45° wrong).
AIMNet2 = 0.971 (14° off). Timing: 0.17s/frame on RTX 5090 for
4876 atoms. Eliminates MOPAC from ensemble pipeline.

**Model file:** `aimnet2_wb97m_0.jpt` (TorchScript archive, ~50MB).
Downloaded from github.com/zubatyuk/aimnet-model-zoo.

**Dependencies:** libtorch (CUDA), linked via CMake `find_package(Torch)`.
SpatialIndexResult (nanoflann KD-tree), EnrichmentResult (backbone/
aromatic classification for EFG decomposition).

**Input tensors (built in C++):**

| Tensor | Shape | Type | Source |
|--------|-------|------|--------|
| coord | (N, 3) | float32 | ProteinConformation positions |
| numbers | (N,) | int32 | Protein atom elements |
| charge | (1,) | float32 | 0.0 (neutral protein) |
| mol_idx | (N,) | int32 | all zeros (one molecule) |
| nbmat | (N+1, max_nb) | int32 | nanoflann radius search at model cutoff (~5Å) |
| nbmat_lr | (N+1, max_nb_lr) | int32 | nanoflann radius search at 15Å (long-range Coulomb) |
| cutoff_lr | (1,) | float32 | 15.0 (use DSF Coulomb, not default inf) |

**Cutoff:** Read from .jpt model via `model->attr("cutoff").toDouble()`
(returns 5.0Å — confirmed working).  Also `cutoff_lr` available from
model (defaults to inf; we override to 15.0 for DSF).  Both in TOML
and GeometryChoice per the rules.

**Neighbour list construction:**

The model stores its own cutoff (~5Å).  For each atom i, find all
atoms j within cutoff using the existing nanoflann KD-tree.  Build
a dense padded half-neighbour matrix:

```cpp
// For each atom i
auto results = kdtree.radiusSearch(pos[i], cutoff_sq);
for (auto& [j, dist] : results) {
    if (j > i) {  // half list only
        nbmat[i][count_i++] = j;
        nbmat[j][count_j++] = i;
    }
}
// Pad remaining slots with N (sentinel)
// Last row (index N) is all-sentinel padding row
```

The model expects this padded (N+1, max_nb) int32 matrix.

**Output:** Per-atom charges from `result["charges"]`, aim embedding
from internal representation.  Then compute Coulomb EFG using existing
CoulombResult patterns but with AIMNet2 charges instead of ff14SB.

**Binary validation gate (non-negotiable):**

The Python AIMNet2 reference implementation is installed at
`/tmp/dxtb_test/bin/python3` with AIMNet2 from `/tmp/aimnet2_repo/`.
(Reinstall if the venv is gone: `pip install torch torch-cluster
numba requests ase` then `pip install -e /tmp/aimnet2_repo` — or
clone fresh from `github.com/isayevlab/AIMNet2`.)

The model file auto-downloads on first use to
`aimnet2calc/assets/aimnet2/aimnet2_wb97m_0.jpt`.  This is the
same .jpt file the C++ libtorch integration loads.

For EVERY test protein in the test trajectory set:

1. Run Python AIMNet2 on the same PDB/coordinates, extract:
   - charges (N,) float32
   - aim embedding (N, 256) float32
2. Run C++ AIMNet2Result on the same coordinates, extract the
   same arrays from the NPY output.
3. Compare element-by-element.  Gate: max |diff| < 1e-5 for
   float32 quantities.  Both implementations load the same .jpt
   model and build the same neighbour list — if they disagree,
   the C++ neighbour list or tensor formatting is wrong.
4. Compute Coulomb EFG from both charge sets using identical
   positions.  EFG arrays must be binary-identical (both use
   float64 math on identical float32 charges cast to float64).

Write `tests/validate_aimnet2.py`: takes a protein directory,
runs both paths, reports pass/fail per array with max absolute
difference.  This becomes part of the smoke test suite.

This is the standard we use for every C++ calculator.  The
compiler is the adversary.  Binary comparison is the gate.

**Output files:**
- `aimnet2_charges.npy` — (N,) float64, per-atom Hirshfeld charges
- `aimnet2_aim.npy` — (N, 256) float64, learned electronic structure embedding
  (AIMNet2's internal per-atom representation — encodes hybridisation,
  polarisability, conjugation, charge transfer.  191/256 dims active for C.
  Geometry-dependent, changes per frame.  Zero extra cost — already in
  the forward pass.  The forward bet.)
- `aimnet2_efg.npy` — (N, 9) float64, full SphericalTensor from AIMNet2 charges
- `aimnet2_efg_aromatic.npy` — (N, 9) float64, aromatic decomposition
- `aimnet2_efg_backbone.npy` — (N, 9) float64, backbone decomposition

**charge_sensitivity (on Atom, computed once):**
- `aimnet2_charge_sensitivity.npy` — (N,) float64, per-atom
  intrinsic polarisability proxy.

  Method: 10 bulk perturbations (all atoms shifted by random
  0.1Å displacements simultaneously), compute AIMNet2 charges
  for each, take per-atom variance.  This works because the
  charge response is overwhelmingly local: tested at 0.046 mean
  |dq| within 3Å vs 0.00006 at 6-10Å (750× falloff).  The
  off-diagonal contamination from neighbours moving is negligible.

  Cost: 10 forward passes × 0.17s = ~2s for 4876 atoms.
  Validated on 14 calibration proteins, 2845 carbon atoms.

  OBJECT MODEL: Computed on first AIMNet2Result::Compute call.
  Checks if atom.aimnet2_charge_sensitivity is populated; if not,
  runs the 10 perturbations and stores the result on Atom (Protein
  level, permanent).  One if statement.  Available to all subsequent
  conformations without recomputation.

  WHY THIS MATTERS: as a kernel interaction term (sensitivity ×
  kernel T2), it adds +0.08 R² for carbon (0.361 → 0.441).
  That's the third-largest single improvement we've found after
  kernel scales (+0.06) and mutation type (+0.12).  The sensitivity
  lets the model weight kernels differently for sp2 vs sp3 carbons
  — high-polarisability carbons (carbonyl, aromatic) respond
  differently to the same geometric field than low-polarisability
  ones (methyl, Cα).  The scalar value alone adds zero; it's the
  INTERACTION with the kernels that carries signal.  112× variation
  across carbon atoms.

  BONUS: if the .jpt model supports autograd (test early),
  d(charges)/d(positions) via one backward pass gives the exact
  per-atom charge sensitivity tensor without finite differencing.
  Worth 5 minutes of testing.  If it works, replaces the 10
  perturbations with one backward pass (~0.3s).

---

## Change 2: Small calculators — 2 new Results, 2 extensions

No libtorch, no GPU.  Positions + nanoflann + lookup tables.

Minimise new ConformationResults.  Fold features into existing
writers where the existing Result already iterates the right
objects and has the right dependencies.

### 2a. DsspResult extension → residue environment writer

DsspResult already iterates residues, knows backbone atoms (N, CA,
C, O), tracks previous residue for H-bonds, and has the Protein
reference.  Extend it to be the single residue environment writer.

**Add to existing DsspResult::WriteFeatures:**

| Feature | Shape per residue | Notes |
|---------|------------------|-------|
| Full 8-class SS one-hot | (8,) | H/G/I/E/B/T/S/C.  Currently 2 cols. |
| DSSP H-bond energies | (4,) | NH→O_1, NH→O_2, O→NH_1, O→NH_2 |
| Chi angles (chi1-4) | (4,) cos + (4,) sin + (4,) exists = (12,) | From Residue::chi atom indices.  NaN for missing. |
| HSE | (3,) | up_count, down_count, angle.  CA→CB direction, nanoflann heavy atom count in each hemisphere.  GLY = NaN. |
| S2 order parameter | (1,) | Zhang & Brüschweiler 2002 contact model.  Heavy atom density around backbone N-H.  Leave in/out as model experiment. |
| Random coil shifts | (6,) | Wishart 1995 lookup: 20 AA × 6 nuclei.  Proline correction.  No geometry — pure residue_type lookup. |

HSE and S2 need SpatialIndexResult (nanoflann).  Add as a
dependency of DsspResult if not already present.

Output: expand `dssp_backbone.npy` or write a companion
`dssp_residue_environment.npy` with the new columns.

### 2b. EnrichmentResult extension

Already computed, never written.  Add WriteFeatures:

| File | Shape | Notes |
|------|-------|-------|
| `atom_role.npy` | (N,) int32 | BackboneN, SidechainC, AromaticH, MethylH, etc. |
| `hybridisation.npy` | (N,) int32 | sp, sp2, sp3 |
| `graph_dist_ring.npy` | (N,) int32 | Bond hops to nearest ring atom |
| `is_conjugated.npy` | (N,) int8 | Boolean |
| `bfs_to_nearest_ring_atom.npy` | (N,) int32 | BFS distance |

### 2c. SasaResult (NEW, ~100 lines)

Genuinely new — per-ATOM Shrake-Rupley is a different algorithm at
a different granularity than DSSP's per-residue SASA.

**Algorithm:** For each atom, distribute ~92 points on a sphere of
radius (r_vdW + r_probe).  For each point, check if any neighbour
atom's vdW sphere occludes it (nanoflann query).
SASA = fraction_exposed × 4π(r_vdW + r_probe)².

**Parameters (all TOML):**
- r_probe = 1.4 Å (water)
- vdW radii: standard table (Bondi or ff14SB)
- n_points = 92 (Fibonacci lattice on sphere)

**Dependencies:** SpatialIndexResult.

**Output:** `atom_sasa.npy` — (N,) float64, Å² per atom.

### 2d. GraphExportResult (NEW, ~50 lines)

Graph structure output for the GNN.  Not a shielding calculator —
just topology and spatial structure.

**Bond graph (from topology, written once):**
- `bonds.npy` — (B, 2) int32, atom index pairs
- `bond_categories.npy` — (B,) int32, bond type enum

**Radius graph (from positions via nanoflann):**
- `radius_edges.npy` — (E, 2) int32, atom pairs within cutoff
- `radius_distances.npy` — (E,) float64

**Parameters (TOML):** radius graph cutoff (default 10Å).
**Dependencies:** SpatialIndexResult.

---

**Summary: 2 new Results + 2 extensions = all UCBSHIFT features covered.**

| Result | Status | Features added |
|--------|--------|---------------|
| DsspResult | Extend | 8-class SS, H-bond energies, chi1-4, HSE, S2, random coil |
| EnrichmentResult | Extend | atom_role, hybridisation, graph topology scalars |
| SasaResult | New | Per-atom Shrake-Rupley |
| GraphExportResult | New | bonds.npy, radius_edges.npy |

---

## Change 3: TrajectoryRunner + TrajectoryObserver + EvaluationAccumulators

**Status:** Architecture for the ensemble streaming loop.  Implement
AFTER all calculators are working.  The GROMACS trajectory loading
(XTC iteration + PBC fix) is a prerequisite and comes first.

---

### The object model

A **ProteinConformationObserver** is a third-lifetime object.
It is not the Protein (topology, eternal).  It is not a
ProteinConformation (one frame, born and dies).  It watches
conformations pass by and builds an understanding of the set.

```
Protein (eternal)
    ↓ owns topology, atoms, bonds, rings
ProteinConformation (per frame)
    ↓ owns positions, all ConformationResults
    ↓ AIMNet2Result, BiotSavartResult, DsspResult, etc.
    ↓ each result attaches computed fields to ConformationAtoms
ProteinConformationObserver (across frames)
    ↓ reads from live conformations before they are discarded
    ↓ holds its own accumulated state
    ↓ may hold whole conformations it decided to keep
```

**The conformation is a real object.**  It carries everything we
computed about the protein at one moment in time.  The observer
reads from ConformationResults on the live conformation — the
full tensors, the charges, the EFG, the DSSP, everything.  It
does not get a pre-digested summary.  It harvests what it needs
from the real thing.

**The observer may hold whole conformations.**  The diversity
filter is an observer that claims conformations instead of
letting them die.  Those survivors are first-class
ProteinConformations with all their ConformationResults intact.
They can be interrogated later for analyses we haven't thought
of yet — like running charge sensitivity perturbations on
survivors, or computing a new feature after the main loop.

**The observer holds whatever it needs.**  Welford accumulators
for statistics.  A buffer of whole conformations for survivors.
Histograms for DSSP.  Whatever serves the accumulation.  Its
internal state is its own business.

### The streaming loop

```cpp
Protein protein = BuildFromGromacs(tpr);
AIMNet2Model aimnet2 = AIMNet2Model::Load(jpt_path);
RunOptions opts = { .aimnet2_model = &aimnet2, ... };

// TrajectoryRunner owns the loop and all I/O.
// TrajectoryObserver holds evaluation results.
// Accumulators hold running statistics.  Nobody holds conformations.

TrajectoryRunner runner(protein, opts, temp_dir, output_dir);

// Observer with accumulators — evaluation only, no I/O
TrajectoryObserver observer(protein);
observer.AddAccumulator<KernelMomentsAccumulator>(protein);
observer.AddAccumulator<AimEmbeddingAccumulator>(protein);
observer.AddAccumulator<ScalarMomentsAccumulator>(protein);
observer.AddAccumulator<DsspHistogramAccumulator>(protein);
observer.AddAccumulator<RingGeometryAccumulator>(protein);
observer.AddAccumulator<FieldMomentsAccumulator>(protein);
observer.AddAccumulator<PhysicsSurvivorAccumulator>(protein, calibration);
observer.AddAccumulator<DiversityAccumulator>(protein);

runner.Run(observer, xtc_paths, colvar_paths, temperature);

// Inside runner.Run():
//
//   weights = ReadColvarWeights(colvar_paths, temperature)
//   XtcReader reader(xtc_path)
//   MoleculeWholer wholer(protein)
//
//   for each frame:
//       positions = wholer.MakeWhole(reader.Positions())
//       conf = MakeConformation(protein, positions)
//       OperationRunner::Run(conf, opts)          // all calculators
//       WriteAllFeatures(conf, temp_dir/frame_i)  // full NPY to temp
//       observer.Observe(conf, frame_i, weights[frame_i])  // accumulators read live conf
//       // conf dies here — memory freed
//
//   observer.Finalize()
//   exemplars = observer.CollectAllExemplars()     // tagged frame indices
//   for each exemplar: promote temp/frame_{idx} → output/survivor_{tag}
//   observer.WriteAllResults(output_dir)           // Welford stats etc
//   delete temp                                    // non-winners gone
```

### Object model

```
TrajectoryRunner
    role: owns the loop and ALL filesystem I/O
    streams: XTC frames via XtcReader + MoleculeWholer
    per frame: builds conformation, calls OperationRunner::Run,
               writes full per-frame results to temp storage,
               hands live conformation to TrajectoryObserver
    finalize: calls TrajectoryObserver::Finalize(),
              gets exemplar indices back,
              promotes winners from temp to output (tagged),
              writes observer's accumulated statistics to output,
              deletes non-winner temp data
    The observer never touches the filesystem.  The runner does all I/O.

TrajectoryObserver
    role: holds evaluation results — per-conformation and overall
    owns: vector<EvaluationAccumulator*>
    per frame: passes live conformation to each accumulator's Observe(),
               holds per-conformation evaluation results
    finalize: calls each accumulator's Finalize(),
              then each accumulator's ListExemplars()
              returns exemplar indices + tags to the runner
    Does NOT hold conformations in memory.
    Does NOT touch the filesystem.

EvaluationAccumulator (ABC)
    role: one axis of evaluation across the trajectory
    Observe(): reads from the live conformation, compares to
               its own running state from all previous frames,
               updates statistics (Welford, argmax, histograms)
    Finalize(): computes final statistics from accumulated state
    ListExemplars(): returns tagged frame indices — this
                     accumulator's survivors, the frames that
                     best exemplify its dimension of the physics
    WriteResults(): writes accumulated statistics to output
    Holds: running statistics + frame indices.  NEVER holds
           conformations in memory.  If an accumulator starts
           holding conformations, that is a bug.

Concrete accumulators:
    KernelMomentsAccumulator    — Welford on kernel T2 tensors
    AimEmbeddingAccumulator     — Welford on (N, 256) aim vectors
    ScalarMomentsAccumulator    — Welford on charges, SASA; circular for angles
    DsspHistogramAccumulator    — per-residue SS class counts
    RingGeometryAccumulator     — Welford on per-ring center/normal/radius
    FieldMomentsAccumulator     — Welford on E-field and EFG tensors
    PhysicsSurvivorAccumulator  — holds calibration coefficients,
                                  tracks argmax per element per kernel
                                  dimension, exemplars are the frames
                                  that define each physics extreme
    DiversityAccumulator        — Mahalanobis from running mean,
                                  exemplars are geometrically diverse
                                  fallback picks
```

### Interfaces

```cpp
class EvaluationAccumulator {
public:
    virtual ~EvaluationAccumulator() = default;

    // Called for each frame.  The conformation is LIVE — all
    // ConformationResults attached, all ConformationAtom fields
    // populated.  Read whatever you need.
    virtual void Observe(const ProteinConformation& conf,
                         size_t frame_index,
                         double boltzmann_weight) = 0;

    // Called once after all frames.
    virtual void Finalize() = 0;

    // This accumulator's survivors: frame indices + tags.
    // Tags describe WHY each frame was selected.
    struct Exemplar {
        size_t frame_index;
        std::string tag;  // e.g. "H_bs_phe_max", "N_pq_total_peak"
    };
    virtual std::vector<Exemplar> ListExemplars() const = 0;

    // Write accumulated statistics (not exemplars — those are
    // full per-frame outputs already on disk).
    virtual int WriteResults(const std::string& output_dir) const = 0;
};


class TrajectoryObserver {
    Protein& protein_;
    std::vector<std::unique_ptr<EvaluationAccumulator>> accumulators_;

public:
    // Per frame: each accumulator reads from the live conformation.
    // No I/O here — the runner handles disk writes.
    void Observe(const ProteinConformation& conf,
                 size_t frame_index,
                 double boltzmann_weight) {
        for (auto& acc : accumulators_)
            acc->Observe(conf, frame_index, boltzmann_weight);
    }

    // After all frames: finalize each accumulator.
    void Finalize() {
        for (auto& acc : accumulators_)
            acc->Finalize();
    }

    // Collect all exemplars from all accumulators.
    // Returns tagged frame indices — the runner handles promotion.
    std::vector<EvaluationAccumulator::Exemplar> CollectAllExemplars() const {
        std::vector<EvaluationAccumulator::Exemplar> all;
        for (auto& acc : accumulators_)
            for (auto& ex : acc->ListExemplars())
                all.push_back(ex);
        return all;
    }

    // Write all accumulated statistics.  Called by the runner.
    void WriteAllResults(const std::string& output_dir) const {
        for (auto& acc : accumulators_)
            acc->WriteResults(output_dir);
    }
};
```

### What each observer accumulates

**KernelMomentsObserver:**
Reads each calculator's shielding SphericalTensor (full 9 components)
from ConformationAtom.  Welford mean + M2, both Boltzmann-weighted
and unweighted.  Also accumulates AIMNet2 aim embedding (N, 256).

**RingGeometryObserver:**
Reads ring center/normal/radius from GeometryResult on the
conformation.  Welford moments per ring.

**ScalarMomentsObserver:**
Reads per-atom scalars: AIMNet2 charges, SASA, backbone dihedrals.
Welford for non-angular quantities, circular Welford (accumulate
sin/cos) for phi/psi/chi.

**DsspHistogramObserver:**
Reads DSSP secondary structure classification per residue.
Accumulates a count histogram per residue (H/G/I/E/B/T/S/C).
Finalize divides by frame count → fractions.

**PhysicsInformedSurvivorObserver:**
Holds the calibration ridge coefficients and per-element kernel
importance rankings.  Tracks running argmax per element per
kernel dimension — which frame index maximised what.

**Does NOT hold conformations in memory.**  Every frame's full
NPY output is written to disk as it passes through (~3ms per
frame write).  The observer holds indices and running statistics
only.

During streaming:
1. Each frame's full results are written to `frame_{index}/`
2. Observer projects kernel values into per-element prediction
   space using calibration coefficients
3. Updates running argmax: "frame 247 is the new BS_PHE peak
   for hydrogen"
4. Updates Welford accumulators
5. Conformation object dies (memory freed)

At Finalize:
1. Observer examines its running argmax table and diversity scores
2. Identifies winner frames per element per dimension
3. Winners' complete results are ALREADY ON DISK
4. Creates survivor entries (symlink or copy from frame_{index}/
   to survivor_H_bs_phe_max/)
5. Optionally: go back to winner frames and run additional analysis
   (e.g., charge sensitivity perturbations on just those frames)

**All choices are preserved.**  If the DFT validation later reveals
a dimension nobody anticipated, go back to the per-frame output
on disk and find the winners for that dimension.  No re-extraction.
The disk is the survivor buffer.

The number of survivors is driven by the physics:
- H needs ~3 (ring current is 2 kernels, plus near/far field)
- C needs the AIMNet2 EFG extremes
- N needs ~15 (multi-mechanism, each kernel's peak)
- Plus a few for raw geometric diversity as fallback

Disk: write all frames per protein as temporary working space.
After Finalize, keep winners + Welford output (~200MB per protein),
delete the per-frame data, start the next protein.  Temporary
disk usage is large per protein but doesn't accumulate.

### The observer as experiment designer

The observer sees the entire trajectory.  Its output supports
multiple competing training strategies from ONE extraction run:

**Strategy 1: Ensemble statistics.**  Welford means + variances
as node features on one graph.  The main design from session 1.

**Strategy 2: UCSB-comparable baseline.**  Boltzmann minimum
frame only, per-residue features (phi/psi/chi + ring current T0
+ E-field + H-bond + SASA).  Exactly what SPARTA+/ShiftX2
would compute.  The "beat them on their own terms" control.

**Strategy 3: Physics-informative multi-pose.**  The tagged
survivors, each as a separate training example with full features.
Same experimental target for all conformations of the same protein.
GNN sees real conformations, not statistics.

**Strategy 4: Refined geometry (no kernels).**  Dihedral ensemble
statistics + per-ring distance/theta/rho/z distributions.  No T2
tensors, no calibrated kernels.  Isolates whether our geometric
description beats UCSB even without the physics engine.

The comparison between strategies IS thesis content.  One
extraction run, four experimental designs, four results chapters.

### Visualisation output

The observer writes per-frame per-dimension time series:
per-element ridge projection magnitude, kernel dominance scores,
ring geometry evolution.  The VTK viewer (nmr-viewer) renders
these as animations of the protein's physics evolving through
its trajectory.  Thesis figures from the same extraction.

**FieldMomentsObserver:**
Reads per-atom E-field vector and EFG tensor (both ff14SB Coulomb
and AIMNet2).  Welford moments.

### What this writes

```
output_dir/
  # Welford statistics
  ensemble_kernel_mean_boltz.npy    (N, K, 9)
  ensemble_kernel_mean_unwtd.npy    (N, K, 9)
  ensemble_kernel_var_boltz.npy     (N, K, 9)
  ensemble_kernel_var_unwtd.npy     (N, K, 9)
  ensemble_aim_mean.npy             (N, 256)
  ensemble_aim_var.npy              (N, 256)
  ensemble_ring_geometry.npy        (R, 14) mean+var center/normal/radius
  ensemble_charges_mean.npy         (N,)
  ensemble_charges_var.npy          (N,)   ← this IS the dynamics charge sensitivity
  ensemble_sasa_mean.npy            (N,)
  ensemble_sasa_var.npy             (N,)
  ensemble_dssp_histogram.npy       (N_res, 8) SS class fractions
  ensemble_field_moments.npy        (N, ...) E-field and EFG moments
  ensemble_meta.npy                 n_frames, n_eff, weight_entropy

  # Survivors: full single-frame output per survivor
  survivor_00/
    pos.npy, bs_shielding.npy, aimnet2_charges.npy, ... (all NPYs)
  survivor_01/
    ...
  survivor_15/
    ...

  # Static (from Protein topology + conformation zero)
  element.npy, residue_type.npy, atom_role.npy, hybridisation.npy
  aimnet2_charge_sensitivity.npy
  bonds.npy, bond_categories.npy
  radius_edges.npy, radius_distances.npy
  kernel_scales.npy
```

### GROMACS trajectory loading (prerequisite)

The streaming loop needs:

1. **XTC reader** — ported from fes-sampler (`xtc_reader.h` +
   `xdrfile/`).  VERBATIM.  Reads one frame at a time.

2. **PBC fix** — ported from fes-sampler (`pbc_whole.h`,
   MoleculeWholer).  VERBATIM.  The frankenmolecule fix was
   validated on 685 proteins.  Without it, atoms from different
   periodic images get mixed and you get 72Å RMSDs.

3. **COLVAR reader** — parse Boltzmann weights from walker
   COLVAR files (6-column text: time, rmsd, rg, bias, rbias, rct).
   Weight = exp(rbias / kT), normalised via log-sum-exp.
   Temperature from plumed.dat.

4. **Frame subsampling** — configurable stride (default every ~10th
   frame from 5005 → ~500 frames).  TOML parameter.

These are infrastructure, not calculators.  They come before the
observer framework because the observers need frames to observe.

---

## What we are retiring

**MOPAC is retired from the ensemble pipeline.**  AIMNet2 replaces
it for charges and charge-dependent EFG at 0.17s/frame vs 10+ min.
MOPAC remains available for single-structure calibration extraction
(the existing `--mutant` mode) but is not called in the streaming loop.

**APBS is retired from the ensemble pipeline.**  DeltaAPBS_EFG
contributed R²=0.005.  Not worth the 15-30s/frame cost.

The geometry-only calculators (BS, HM, MC, PQ, RingSusc, Disp,
HBond, Coulomb) remain unchanged.  They run at ~1.3s/frame and
provide the core kernel features.

---

## Calibration warm start

The per-element ridge coefficients from the 723-protein calibration
(R²=0.818) are shipped as a metadata file alongside the extraction.
The GNN uses these to pre-weight the atom-ring edge features
(Option A: edge feature = calibrated_coefficient × kernel T2).

Files (computed once, not per protein):
- `calibration/ridge_coefficients.npy` — (4, 55, 5, 5) per-element
- `calibration/kernel_importance.npy` — (4, 55) rankings
- `calibration/ring_intensities.npy` — (8,) I per ring type

---

## Priority order

1. **AIMNet2Result** — DONE.  libtorch CUDA, binary validated.

2. **GROMACS trajectory infrastructure** — XTC reader + PBC fix
   ported from fes-sampler VERBATIM, COLVAR weight reader.  Next
   because PBC fix is critical and the session has context.

3. **Small calculators (Change 2):**
   - DsspResult extension (8-class SS, H-bond energies, chi1-4,
     HSE, S2, random coil shifts)
   - EnrichmentResult extension (atom_role, hybridisation, graph
     topology scalars)
   - SasaResult (new, per-atom Shrake-Rupley)
   - GraphExportResult (new, bonds + radius edges)

4. **TrajectoryRunner + TrajectoryObserver + EvaluationAccumulators
   (Change 3).**  The streaming architecture.  Implement AFTER all
   calculators work so there's something to observe.

5. **Integration test** — full streaming loop on fleet test protein
   (1Q8K_10023, ~500 frames).

6. **Fleet extraction** — 685 proteins, ~500 frames each.

---

## Issues from review (2026-04-11)

1. **aim embedding accumulator:** BoltzmannKernelMoments or a new
   TensorMoments observer must accumulate the (N, 256) aim embedding.
   Output: `ensemble_aim_mean.npy` (N, 256), `ensemble_aim_var.npy` (N, 256).

2. **Kernel scale factors:** Per-protein normalization scale factors
   (what normalize_kernels strips) must be written alongside ensemble
   arrays.  Static per-protein: `kernel_scales.npy` (K,).

3. **atom_role and hybridisation export:** Static per-atom arrays from
   EnrichmentResult.  Currently computed but not written.  Add:
   `atom_role.npy` (N,) int32, `hybridisation.npy` (N,) int32.

4. **Read cutoff from .jpt model:** Use `model.attr("cutoff")` via
   libtorch C++ API.  Do not hardcode 5.0.

5. **nbmat deduplication:** Build in one pass.  Iterate atoms 0..N-1,
   for each nanoflann result (i,j) with j>i, push j into row i and
   i into row j.  Each pair appears exactly once in the search results
   because nanoflann is symmetric.

6. **Circular Welford for angles:** Chi angles, phi, psi all wrap at
   ±π.  Use circular statistics: accumulate sin and cos separately,
   mean = atan2(mean_sin, mean_cos), variance from resultant length.

7. **Radius graph from representative frame:** Write once from the
   Boltzmann minimum frame, not per-frame.  The graph topology is
   fixed; node features change.

8. **DispChi in ring_contributions:** disp_scalar × chi.T2 is the
   validated dispersion feature (from LESSONS_FROM_CALIBRATION).
   Already in per-ring output.  Confirm it's in the accumulator.

9. **Graph topology scalars:** graph_dist_ring, is_conjugated,
   bfs_to_nearest_ring_atom — static per-atom, already computed in
   C++ enrichment.  Need WriteFeatures.  Add alongside atom_role.

---

## Speculative: per-atom polarisability from AIMNet2

AIMNet2 is a differentiable model.  Per-atom energy is available.
The polarisability tensor α = -d²E/dF² (second derivative of energy
with respect to electric field).  If computed via autograd, this
gives per-atom polarisability anisotropy at no extra frame cost.

This would be a genuinely novel feature — the Buckingham γ coefficient
that converts EFG to shielding, estimated per-atom instead of
per-element.  Not yet validated.  The long bet.
