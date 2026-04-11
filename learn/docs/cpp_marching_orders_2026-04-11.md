# C++ Marching Orders: Ensemble Extraction

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

## Change 2: Per-atom SASA (~100 lines)

**What:** Shrake-Rupley solvent-accessible surface area per atom.

**Why:** UCBShift/SPARTA+ use SASA as a key feature.  DSSP gives
per-residue SASA (broadcast to all atoms in residue — same value).
Per-atom SASA distinguishes buried backbone atoms from exposed ones
within the same residue.

**Algorithm:** For each atom, distribute ~92 points on a sphere of
radius (r_vdW + r_probe).  For each point, check if any neighbour
atom's vdW sphere occludes it (neighbour query via nanoflann).
SASA = fraction_exposed × 4π(r_vdW + r_probe)².

**Parameters:**
- r_probe = 1.4 Å (water)
- vdW radii: standard table (Bondi or ff14SB)
- n_points = 92 (Fibonacci lattice on sphere)
- Neighbour cutoff = max_vdw + max_vdw + 2 * r_probe ≈ 6 Å

**Dependencies:** SpatialIndexResult (for nanoflann KD-tree).

**Output:** `atom_sasa.npy` — (N,) float64, Å² per atom.

---

## Change 3: Chi angles (~30 lines)

**What:** Sidechain dihedral angles (chi1-chi4) per residue.

**Why:** Backbone phi/psi are already in dssp_backbone.npy.  Sidechain
chi angles are missing.  UCBShift uses chi1/chi2.  For the ensemble,
chi angle distributions (rotamer populations) are physically meaningful.

**Implementation:** The infrastructure exists — `Residue::chi[4]` holds
4-atom index tuples, resolved at PDB load time from `AminoAcidType::
chi_angles` definitions.  Compute the dihedral from the 4 atom
positions (atan2 of cross products).  Write per-residue.

**Output:** `chi_angles.npy` — (N_residues, 4) float64, radians.
NaN for chi angles that don't exist (e.g., GLY has 0, ALA has 0,
VAL has 1, PHE has 2).

---

## Change 4: Bond graph export (~20 lines)

**What:** Covalent bond pairs from topology.

**Why:** The GNN needs the bond graph as an edge type.  We have bonds
internally (McConnell iterates them) but never write them.

**Output:**
- `bonds.npy` — (B, 2) int32, atom index pairs
- `bond_categories.npy` — (B,) int32, bond type enum (PeptideCO,
  PeptideCN, BackboneOther, SidechainCO, Aromatic, SidechainOther)

---

## Change 5: Radius graph export (~30 lines)

**What:** All atom pairs within a configurable cutoff.

**Why:** Standard GNN input.  The GNN needs (i, j) pairs with
distances.  Currently the GNN would have to recompute this from
positions, which is redundant — we already have nanoflann.

**Parameters:** cutoff from config (default 10 Å).

**Output:**
- `radius_edges.npy` — (E, 2) int32, atom index pairs
- `radius_distances.npy` — (E,) float64

---

## Change 6: Ensemble observer framework — FUTURE, DO NOT IMPLEMENT NOW

**Status:** Design reference only.  Implement AFTER all calculators
(AIMNet2, SASA, chi angles, graph exports) are working and validated.
The streaming loop and Welford accumulators are a separate project.

**What:** The streaming accumulation from ENSEMBLE_MODEL.md, with
revisions based on today's analysis.

**Observers:**

1. **BoltzmannKernelMoments** — Per-atom, per-kernel, full 9-component
   tensor Welford (mean + M2).  Both Boltzmann-weighted and unweighted.
   Includes all geometry-only calculators AND AIMNet2 EFG.

2. **RingGeometryStats** — Per-ring center/normal/radius Welford moments.

3. **ScalarMoments** — Per-atom Welford for: AIMNet2 charges, SASA,
   chi angles, backbone dihedrals.

4. **DsspEnsembleStats** — Per-residue secondary structure histogram
   (fraction of frames in each DSSP class).

5. **FrameDiversityFilter** — Keeps 10-20 diverse survivor
   conformations (Mahalanobis distance from running kernel mean).
   Survivors get full WriteFeatures — their per-atom features become
   additional node channels for the GNN.

6. **CoulombFieldMoments** — Per-atom E-field vector and EFG tensor
   Welford moments (ff14SB and AIMNet2 separately).

**Output per ensemble:**
- `ensemble_kernel_mean_boltz.npy` — (N, K, 9) Boltzmann mean
- `ensemble_kernel_mean_unwtd.npy` — (N, K, 9) unweighted mean
- `ensemble_kernel_var_boltz.npy` — (N, K, 9) Boltzmann variance
- `ensemble_kernel_var_unwtd.npy` — (N, K, 9) unweighted variance
- `ensemble_ring_geometry.npy` — (R, 14) mean + var of center/normal/radius
- `ensemble_dssp_histogram.npy` — (N_res, 8) SS class fractions
- `ensemble_n_frames.npy` — (1,) total frames processed
- `ensemble_n_eff.npy` — (1,) effective sample size from weights
- `ensemble_weight_entropy.npy` — (1,) entropy of Boltzmann weights
- `survivor_*/` — subdirectories with full single-frame features for
  each survivor conformation

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

**Current scope — calculators only:**

1. AIMNet2Result — libtorch CUDA + model load + neighbour lists +
   ConformationResult + EFG + charge_sensitivity + binary validation.
   Includes porting xdrfile/xtc_reader.h/pbc_whole.h from fes-sampler.
2. Chi angles + per-atom SASA + bond/radius graph export + atom_role
   + hybridisation + graph topology scalars (small calculators)

**Future scope — NOT current work:**

3. Ensemble observer framework (Welford accumulators — needs all
   calculators working first so it has something to observe)
4. XTC streaming loop integration
5. Full ensemble run on fleet test protein (1Q8K_10023)
6. Full fleet extraction (685 proteins, ~500 subsampled frames each)

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
