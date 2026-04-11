# Next Session Prompt (after 2026-04-11 design sessions)

## Reading order — DO THIS FIRST

### 1. Understand the system (before touching AIMNet2)

Read the spec index for orientation:
- `spec/INDEX.md` — reading order, document tiers, what to read when

Then read an existing calculator to learn the ConformationResult pattern:
- `src/CoulombResult.h` + `src/CoulombResult.cpp` — the closest
  analog to AIMNet2Result.  Understand Dependencies(), Compute(),
  WriteFeatures().  See how it queries SpatialIndexResult for
  neighbour lists, stores results on ConformationAtom, decomposes
  EFG by source (backbone/aromatic).
- `src/ConformationResult.h` — base class interface
- `src/ConformationAtom.h` — where per-frame computed fields live.
  All `aimnet2_*` fields will be added here.
- `src/SpatialIndexResult.h` + `.cpp` — nanoflann KD-tree.
  AIMNet2 neighbour lists are built from this.
- `src/EnrichmentResult.h` + `.cpp` — backbone/aromatic/sidechain
  classification.  AIMNet2 EFG decomposition depends on this.

Understand the object model walls:
- `spec/OBJECT_MODEL.md` — every class, property, type, unit.
  Code against this.
- **Protein = identity/topology ONLY.  Never put geometry on it.**
  Geometry lives on ProteinConformation.  The wall is absolute.
- charge_sensitivity is the ONE exception: an intrinsic property
  of the atom (like hybridisation) that is computed from geometry
  but stored permanently on Atom.  See design rationale below.

Understand the parameter discipline:
- `data/calculator_params.toml` — every cutoff, radius, threshold.
  AIMNet2 params go here.
- `spec/GEOMETRY_CHOICE_BRIEF.md` — GeometryChoice recording spec.
  Every parameter used during Compute() gets a GeometryChoice entry.

### 2. Understand the AIMNet2 spec

- `learn/docs/cpp_marching_orders_2026-04-11.md` — the full spec.
  Read Change 1 (AIMNet2Result) thoroughly.  Changes 2-5 (SASA,
  chi, graphs) are separate work.  Change 6 (ensemble observers)
  is FUTURE — do not implement.

### 3. Understand the PBC fix to port

- `/shared/2026Thesis/fes-sampler/src/pbc_whole.h` — MoleculeWholer.
  Port VERBATIM.  Do not redesign.
- `/shared/2026Thesis/fes-sampler/src/xtc_reader.h` — XTC frame I/O.
  Port VERBATIM.
- `/shared/2026Thesis/fes-sampler/xdrfile/` — bundled C library.
  Copy as-is.
- `/shared/2026Thesis/fes-sampler/CLAUDE.md` — context on the fix
  (frankenmolecule PBC failures on 685 proteins).

---

## What to build: AIMNet2Result calculator

**AIMNet2Result is a ConformationResult.  It runs on any single
ProteinConformation.  Same pattern as CoulombResult.  Nothing more.**

### DO NOT implement

- Streaming loops, frame iteration, ensemble processing
- Welford accumulators, EnsembleConformation, EnsembleResult
- "Conformation zero" as a concept
- Abstract NeuralChargeProvider interfaces, factories, or pluggable
  abstractions of any kind
- Model-passing ceremonies — the model is loaded once as a static
  resource

### Design decisions (settled, non-negotiable)

- **Married to AIMNet2.**  No abstractions.  Direct `aimnet2_*`
  named fields everywhere.  If a better model arrives, write a
  new calculator.  A plugin swap hides the differences.
- **CUDA mandatory.**  `find_package(Torch REQUIRED)`.  No CPU
  fallback.  RTX 5090 is the target.
- **charge_sensitivity on Atom.**  First Compute checks if Atom
  has it.  If not, runs 10 bulk perturbations (~2s), stores on
  Atom permanently.  One if statement.  Not a framework.
- **PBC from fes-sampler VERBATIM.**  Copy the fix.  Do not
  redesign.  The multi-chain handling, topology trimming, and
  do_pbc_mtop call were validated on 685 production proteins.
- **Binary validation gate.**  C++ must match Python AIMNet2 to
  1e-5 for charges, binary-identical for EFG.

### Implementation order

1. Port fes-sampler PBC infrastructure into nmr-shielding
   (xdrfile as static lib, xtc_reader.h, pbc_whole.h)
2. libtorch CUDA CMake integration + model load + toy forward pass
3. Neighbour list construction from nanoflann (nbmat + nbmat_lr,
   half-list, padded, sentinel=N)
4. AIMNet2Result::Compute — build input tensors, forward pass,
   extract charges + aim embedding, compute Coulomb EFG.
   Store all on ConformationAtom.  Record GeometryChoices.
5. charge_sensitivity lazy-init on Atom
6. WriteFeatures (6 NPY files)
7. Binary validation against Python reference
   (tests/validate_aimnet2.py)

### Key references

| What | Where |
|------|-------|
| Closest analog calculator | `src/CoulombResult.h/.cpp` |
| Base class | `src/ConformationResult.h/.cpp` |
| Per-frame atom fields | `src/ConformationAtom.h` |
| Nanoflann KD-tree | `src/SpatialIndexResult.h/.cpp` |
| Backbone/aromatic classification | `src/EnrichmentResult.h/.cpp` |
| Object model | `spec/OBJECT_MODEL.md` |
| Parameter patterns | `data/calculator_params.toml` |
| GeometryChoice spec | `spec/GEOMETRY_CHOICE_BRIEF.md` |
| NPY writer | `src/NpyWriter.h` |
| Operation runner | `src/OperationRunner.h/.cpp` |
| PBC fix | `/shared/2026Thesis/fes-sampler/src/pbc_whole.h` |
| XTC reader | `/shared/2026Thesis/fes-sampler/src/xtc_reader.h` |
| xdrfile library | `/shared/2026Thesis/fes-sampler/xdrfile/` |
| CMake build | `CMakeLists.txt` |
| AIMNet2 full spec | `learn/docs/cpp_marching_orders_2026-04-11.md` |
| AIMNet2 Python ref | `/tmp/aimnet2_repo/` (github.com/isayevlab/AIMNet2) |
| AIMNet2 model | `aimnet2_wb97m_0.jpt` (auto-downloads on first use) |
| Model cutoff API | `model->attr("cutoff").toDouble()` = 5.0 |
| Test protein | `tests/data/fleet_test_large/` (1Q8K_10023, 4876 atoms) |

### AIMNet2 Python environment

- Venv: `/tmp/dxtb_test/` (may need reinstall if ephemeral)
- Install: `pip install torch torch-cluster numba requests ase`
  then `pip install -e /tmp/aimnet2_repo`
- Or clone fresh: `github.com/isayevlab/AIMNet2`

---

## Context (why, not what to implement)

### What happened on 2026-04-11

Session 1 established mathematical foundation for ensemble NMR
prediction (6 findings — see ensemble_session_2026-04-11.md).
Key result: AIMNet2 recovers carbon EFG dimension (cos 0.71→0.97
vs MOPAC) at 0.17s/frame.

Session 2 refined the implementation plan:
- Separated calculator from ensemble/accumulator scope
- Married to AIMNet2 (no abstractions)
- charge_sensitivity on Atom (topology-like)
- PBC from fes-sampler verbatim
- Eliminated streaming loop from current scope

### Future work (NOT current scope)

- Ensemble observer framework (spec/ENSEMBLE_MODEL.md — reference only)
- Streaming XTC loop
- Welford accumulators
- Full fleet extraction (685 proteins)
