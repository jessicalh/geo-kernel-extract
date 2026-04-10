# Ensemble Object Model

**Status:** Planned. Design established 2026-04-10.

---

## Three lifetimes

The system has three distinct object lifetimes for per-atom data:

| Object | Lifetime | Holds | Example |
|---|---|---|---|
| Protein / Atom | Eternal | Identity, topology, residue, bonds, rings | element, pdb_atom_name, residue_index |
| ProteinConformation / ConformationAtom | One frame | Positions, instantaneous computed fields | bs_shielding_contribution, coulomb_E_total |
| **EnsembleConformation / EnsembleAtom** | **Across all frames** | **Accumulated statistics** | **Boltzmann-weighted kernel T2 mean, variance** |

Protein owns what doesn't change. ConformationAtom owns what one
geometry produces. EnsembleConformation owns what 5005 geometries
collectively reveal.

**Atom index is the join key across all three.**

---

## EnsembleConformation

Parallel to ProteinConformation. Same atom count, same indexing.
Attaches to the Protein the same way a ProteinConformation does.

Where a ProteinConformation holds positions and computed fields
from one frame, an EnsembleConformation holds the reduced
statistics from all frames: means, variances, counts, extrema,
correlation matrices — whatever the EnsembleResults computed.

It does not hold positions (there is no single position for an
ensemble). It does not duplicate identity (that's on Protein).

---

## EnsembleResult

Parallel to ConformationResult. Same patterns:

- `Dependencies()` — declares what other EnsembleResults it reads
- `WriteFeatures()` — writes NPY arrays
- Attaches to EnsembleConformation
- Chaining works: EnsembleResult B can depend on EnsembleResult A

The difference: a ConformationResult is populated by a single
`Compute(conf)` call. An EnsembleResult is populated by many
`Observe(conf, frame_i, weight)` calls followed by one
`Finalize()`.

### Concrete EnsembleResults (planned)

| Result | Accumulates | Output arrays |
|---|---|---|
| BoltzmannKernelMoments | Per-atom Welford mean + M2 for each kernel T2, weighted by Boltzmann probability | (N, K, 5) mean, (N, K, 5) variance |
| DsspEnsembleStats | Per-residue backbone angle circular statistics | (R, 5) mean phi/psi/sasa + circular variance |
| RingGeometryStats | Per-ring center/normal distributions | (rings, 6) mean center + normal, (rings, 6) variance |
| CoulombFieldMoments | Per-atom E-field and EFG statistics | (N, 3) mean E, (N, 9) mean EFG, variances |
| FrameDiversityScores | Per-frame Mahalanobis distance from running kernel mean | (F,) score per frame |

---

## The streaming loop

```
protein = BuildFromGromacs(tpr)
ensemble = EnsembleConformation(protein)  // same atom count

for each XTC frame:
    conf = make_conformation(protein, frame_coords)
    OperationRunner::Run(conf, opts)      // geometric stack

    for each observer on ensemble:
        observer.Observe(conf, frame_i, boltzmann_weight)

    // conf discarded here (unless a filter claims it)

for each observer on ensemble:
    observer.Finalize()

ensemble.WriteAllFeatures(output_dir)     // same pattern as ProteinConformation
```

The Protein is built once from the TPR. Each XTC frame produces
a temporary ProteinConformation that is fully computed (all
geometric calculators attached), observed by the ensemble, then
discarded. The EnsembleConformation survives with the accumulated
statistics.

---

## Filters: observers that keep survivors

A filter is an observer that sometimes holds onto a
ProteinConformation instead of letting it be discarded:

```
void KernelDiversityFilter::Observe(conf, frame_i, weight) {
    score = score_diversity(conf)
    if score > worst_survivor:
        evict worst
        adopt conf into buffer
```

The filter's buffer holds 10-20 full ProteinConformations with
all ConformationResults attached. These survivors are first-class
citizens:

- A second-pass calculator can attach new results to them
  (e.g. run MOPAC on just the 20 winners)
- An EnsembleResult can read across the survivors to compute
  cross-conformation statistics on the diverse subset
- WriteAllFeatures works on each survivor individually

The buffer size is bounded: 10-20 conformations at ~464 MB
worst case (4876 atoms) = 5-10 GB. Fits in 64 GB.

---

## Memory model

| Component | Memory | Notes |
|---|---|---|
| Protein (topology) | ~50 MB | Atoms, bonds, rings, residues. Built once. |
| One ProteinConformation | ~464 MB | Largest protein (4876 atoms). Transient per frame. |
| EnsembleConformation accumulators | ~20 MB | Welford mean + M2, independent of frame count. |
| Filter buffer (20 survivors) | ~5-10 GB | Full conformations with results. |
| COLVAR weights | ~1 MB | 5005 doubles. |

Total steady state during streaming: Protein + one transient
conformation + accumulators + filter buffer. Well within 64 GB.

---

## Boltzmann weights

Read directly from walker COLVAR files (6-column text).
Temperature from plumed.dat. Weight = exp(rbias / kT),
normalised via log-sum-exp.

No fes-sampler dependency. The COLVAR format is the contract.
fes-sampler's dihedral selection remains a standalone tool;
its 6 picks are the control for the ensemble experiment.

---

## Output and downstream

EnsembleConformation participates in WriteAllFeatures the
same way ProteinConformation does. Produces NPY arrays.
The Python SDK loads them. Atom index joins to identity
arrays (element.npy, residue_type.npy).

Downstream consumers:

- **R ridge regression**: ensemble-averaged kernels + RefDB
  experimental shifts, per element per atom type
- **learn/src analysis tools**: same iter_proteins pattern,
  reading ensemble features instead of single-pose features
- **GNN graph builder** (nmr-training Piece 3): ensemble
  features as additional node attributes

---

## Test data

The largest fleet protein (1Q8K_10023, 4876 atoms, 26 rings) is
copied to `tests/data/fleet_test_large/` with all 5 walkers,
COLVAR files, XTC trajectories, initial-samples (6 fes-sampler
poses), and run_params. This is the benchmark protein for memory
and timing measurements.

Empirical measurements (2026-04-10, batcave, Ryzen 9950X 128 GB):
- Per-conformation RSS: ~464 MB (geometry-only, no MOPAC/APBS)
- Per-frame wall time: ~1.3s (with 20A Coulomb cutoff)
- Filter buffer of 20 survivors: ~9 GB worst case

The fleet tree is at `/shared/fleet/clean_fleet/results/`.
Do not modify it. The test copy is the safe playground.

---

## What does NOT change

- Protein / Atom: untouched
- ConformationAtom: untouched (still holds one frame's data)
- ConformationResult: untouched (still computes from one frame)
- OperationRunner::Run: untouched (still sequences one frame)
- WriteFeatures per calculator: untouched
- Python SDK: loads whatever NPY files are present; new arrays
  get new ArraySpec entries in _catalog.py

The wall between Protein and ConformationAtom stays exactly
where it is. EnsembleConformation is the third object on the
other side of a new wall: across-frame statistics that belong
to neither identity nor any single geometry.
