# Remaining Work: GromacsProtein Ensemble Extraction

Cross-conformational features and accumulation that need to be
implemented in GromacsFrameHandler / GromacsFinalResult.

Extracted from cpp_marching_orders_2026-04-11.md (now in bones/).

---

## Accumulation (GromacsFrameHandler)

Per-frame, the handler reads from the live conformation and updates
Welford accumulators. These are direct methods, not abstract observers.

| Quantity | Source on ConformationAtom | Accumulation | Output |
|----------|--------------------------|--------------|--------|
| AIMNet2 charges | aimnet2_charge | Welford mean+var | ensemble_charges_{mean,var}.npy |
| AIMNet2 aim embedding | aimnet2_aim[256] | Welford mean+var | ensemble_aim_{mean,var}.npy |
| Kernel T2 tensors | ring_neighbours[].G_spherical etc | Welford mean+var per kernel | ensemble_kernel_{mean,var}.npy |
| SASA | atom_sasa (from SasaResult) | Welford mean+var | ensemble_sasa_{mean,var}.npy |
| DSSP SS class | secondary_structure char | Count histogram per residue | ensemble_dssp_histogram.npy |
| Ring geometry | ring_geometries[] | Welford on center/normal/radius | ensemble_ring_geometry.npy |
| E-field / EFG | coulomb_E_total, coulomb_EFG_* | Welford moments | ensemble_field_moments.npy |
| Backbone angles | phi, psi (from DsspResult) | Circular Welford (sin/cos) | ensemble_angles.npy |
| Chi angles | dssp_chi.npy values | Circular Welford | included in ensemble_angles |

### Circular Welford for angles
Phi, psi, chi wrap at +/-pi. Accumulate sin and cos separately.
Mean = atan2(mean_sin, mean_cos). Variance from resultant length.

### Charge sensitivity from ensemble
Per-atom charge variance across all frames IS the charge sensitivity.
No extra computation — it's ensemble_charges_var.npy.

---

## Frame selection (GromacsFrameHandler)

Ridge-informed: holds calibration coefficients from the 723-protein
calibration (R²=0.818). For each frame, projects kernel values into
per-element prediction space. Tracks running argmax per element per
kernel dimension. Winners = frames that define physics extremes.

| Selection criterion | What it tracks | N survivors |
|---------------------|---------------|-------------|
| Per-element kernel peaks | argmax of calibrated projection | ~15-20 |
| Geometric diversity | Mahalanobis from running mean | ~5-10 |
| Boltzmann minimum | lowest free energy | 1 |

Winner frames are re-loaded from XTC at the end and get full
WriteAllFeatures.

---

## Finalization (GromacsFinalResult)

1. All accumulators finalize (variance = M2 / (count - 1))
2. Write ensemble NPY files from accumulated state
3. Collect winner frame indices from selection criteria
4. Re-load winner frames from XTC (read_xtc_frame at specific time)
5. PBC fix, create conformation, run all calculators
6. WriteAllFeatures for each winner → output/survivor_{tag}/
7. Write {protein_id}.h5 from Protein topology + conformation 0
8. Log summary

---

## .h5 master file (GromacsFinalResult)

Requires HighFive (header-only C++ HDF5 wrapper) in extern/.

```
{protein_id}.h5
├── topology/          ← from Protein
│   ├── element, residue_type, atom_to_residue
│   ├── bonds, bond_type
│   └── ring_atoms, ring_type
├── classification/    ← from EnrichmentResult on conformation 0
│   ├── atom_role, hybridisation
│   ├── graph_dist_ring, is_conjugated
└── graph/             ← from SpatialIndexResult on conformation 0
    ├── radius_edges, radius_distances
```

SDK (nmr_extract) reads both .h5 and NPY, presents one typed object.

---

## Fleet cost reality

At ~45s/frame (dominated by CoulombResult at 25s), 500 frames per
protein = ~6.25 hrs per protein. 685 proteins on 4 machines = ~45 days.
This is tight for a thesis timeline.

Autograd polarisability adds ~5s/frame. Total ~50s/frame. Scales linearly.

**CoulombResult is the bottleneck.** The current implementation is O(N²)
all-pairs within a 20A cutoff. For 4876 atoms that's ~2.9M pairs.
Options to investigate:
- Reduce Coulomb cutoff (currently 20A — does R² degrade at 15A?)
- Spatial decomposition (cell list instead of nanoflann radius search)
- Skip CoulombResult for streaming frames, run only on survivors
- Ask GROMACS for per-atom Coulomb field directly (it already computes this)

Selection of a representative subset of the 685 proteins may be needed.
~100 proteins at 500 frames each = ~260 hrs = ~3 days on 4 machines.
This is feasible. The subset must cover the element/ring-type diversity
of the full fleet.

---

## Boltzmann weights

Read from walker COLVAR files (6-column text: time, rmsd, rg, bias,
rbias, rct). Weight = exp(rbias / kT). Temperature from plumed.dat.
Normalised via log-sum-exp.

---

## TOML parameters for ensemble

| Parameter | Default | Notes |
|-----------|---------|-------|
| ensemble_frame_stride | 10 | Subsample XTC (5005 → ~500 frames) |
| ensemble_max_survivors | 30 | Max winner frames to keep |

---

## Autograd charge sensitivity (optional, on selected frames)

The .jpt model supports backward() on charges (validated).
d(charges)/d(positions) via N backward passes. ~5s for 5000 atoms.
Run on conformation 0 and survivors only.

TOML: aimnet2_sensitivity_autograd = false (default off).
