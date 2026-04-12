# Ensemble Extraction Model

**Status:** IMPLEMENTED (partial). GromacsProtein + GromacsFrameHandler +
GromacsFinalResult are built and tested. Streaming test passes on
1Q8K_10023 (10 XTC frames, 14 results per frame, 36 NPY arrays).
Accumulation logic and .h5 writer are next.

Design revised 2026-04-11/12 after critical review of original
EnsembleConformation/observer design. See memory:
`project_ensemble_design_decisions.md`.

---

## Architecture: GromacsProtein pattern

No abstract EnsembleConformation. No observer/accumulator ABC.
Three concrete classes married to GROMACS:

### GromacsProtein

Wraps a real Protein. Conformation 0 (GROMACS frame 0) lives in
the Protein's conformations_ vector via AddMDFrame. It is permanent
and feeds the .h5 master file.

Streaming frames are **free-standing ProteinConformations** created
via the public constructor with a `const Protein*` back-pointer.
They are NOT added to the Protein's conformations_ vector. They
point at the Protein for topology lookups, are computed on by all
calculators, then freed.

This works because:
- ProteinConformation constructor is public, no side effects
- No calculator accesses conformations through the Protein's vector
- OperationRunner::Run takes ProteinConformation& directly
- Verified: zero hidden coupling in all calculators and loaders

Holds: the real Protein, ChargeSource from TPR, frame manifest
(paths or indices), accumulation state.

### GromacsFrameHandler

Opens XTC trajectory, applies PBC fix (MoleculeWholer from
fes-sampler, VERBATIM), converts nm→A, creates free-standing
MDFrameConformations. For each frame:

1. Read XTC frame, PBC fix, build positions
2. Create free-standing MDFrameConformation(&protein, positions, ...)
3. OperationRunner::Run — all calculators attach
4. Extract what accumulation needs from the live conformation
5. Conformation dies — memory freed

No disk writes during streaming (unless debugging). The handler
accumulates in memory. Winners are re-loaded from XTC at the end.

### GromacsFinalResult

Runs after all frames. Computes final statistics from accumulated
state. Writes:
- `{protein_id}.h5` — topology + EnrichmentResult from conformation 0
- Ensemble NPY files (Welford means/variances, DSSP histograms, etc.)
- Winner frames: re-load from XTC, run calculators, WriteAllFeatures

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
Lives in the Protein's conformations_ vector via AddMDFrame.
Always valid, always in memory. EnrichmentResult from conformation 0
feeds the .h5 master file (atom_role, hybridisation, etc.).

For `--pdb` or `--mutant` modes: everything works exactly as before.
The existing path is untouched.

---

## Master file: {protein_id}.h5

HDF5 file containing protein-level data that no calculator writes
as NPY. Written once by GromacsFinalResult from conformation 0 +
Protein topology.

```
{protein_id}.h5
├── topology/          ← from Protein object
│   ├── element          (N,) int32
│   ├── residue_type     (N,) int32
│   ├── atom_to_residue  (N,) int32
│   ├── bonds            (B, 2) int32
│   ├── bond_type        (B,) int32
│   ├── ring_atoms       (P, 2) int32
│   └── ring_type        (R,) int32
├── classification/    ← from EnrichmentResult on conformation 0
│   ├── atom_role        (N,) int32
│   ├── hybridisation    (N,) int32
│   ├── graph_dist_ring  (N,) int32
│   └── is_conjugated    (N,) int8
├── graph/             ← from SpatialIndexResult on conformation 0
│   ├── radius_edges     (E, 2) int32
│   └── radius_distances (E,) float64
└── attrs: n_atoms, n_residues, n_rings, pdb_id
```

The SDK (nmr_extract) reads both .h5 and NPY files, presents one
typed object. One SDK, one load() call.

---

## Charge sensitivity

Per-atom charge sensitivity is NOT a calculator feature. It is NOT
on Atom (topology) and NOT on ConformationAtom.

The perturbation approach (10 random bulk displacements, Welford
variance) was removed. It produced conformation-specific values
via random splats that were non-comparable across frames and lied
about solvation when applied to other conformations.

Two quantities replace it, both computed by GromacsFrameHandler:

1. **Ensemble charge variance** — Welford on aimnet2_charge across
   all frames. Free (charges already computed per frame). Output:
   `ensemble_charges_var.npy`.

2. **Autograd charge sensitivity** — d(charges)/d(positions) via
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

- **EnsembleConformation / EnsembleAtom** — replaced by accumulation
  state on GromacsProtein. No separate object model needed.
- **EnsembleResult / observer ABC** — replaced by direct methods on
  GromacsFrameHandler. No virtual dispatch.
- **TrajectoryObserver / EvaluationAccumulator** — same.
- **ConformationList with invalidation** — replaced by free-standing
  conformations that die naturally after processing.
- **GraphExportResult** — bond graph goes in .h5 from Protein topology.
- **charge_sensitivity on Atom** — moved to evaluator/handler output.
  Perturbation approach deleted.

---

## Test data

1Q8K_10023 (4876 atoms, 26 rings) at `tests/data/fleet_test_large/`
with 5 walkers, XTC trajectories, initial-samples (6 PDB poses),
run_params. Streaming test processes 10 XTC frames with all 14
calculators (including SasaResult) at ~45s/frame.

The fleet tree is at `/shared/fleet/clean_fleet/results/`.
Do not modify it.
