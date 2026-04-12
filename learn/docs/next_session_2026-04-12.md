# Next Session Prompt (after 2026-04-11/12 design + implementation sessions)

## What is DONE

### AIMNet2Result — COMPLETE
Implemented, binary-validated, wired into OperationRunner. Runs at
0.17s/frame on RTX 5090. Produces charges, aim embedding (256-dim),
Coulomb EFG (total/backbone/aromatic). 5 NPY files.

### charge_sensitivity — RESOLVED
- Removed from Atom (was violating the topology wall)
- Removed from ConformationAtom (not a per-frame calculator output)
- Perturbation approach DELETED (random splats non-comparable across
  conformations, lies about solvation)
- Replaced by two quantities in GromacsFrameHandler:
  1. Ensemble charge variance (Welford on aimnet2_charge, free)
  2. Autograd d(charges)/d(positions) on selected frames (validated:
     .jpt supports backward(), N backward passes per frame)
- Field removed from Atom.h. Perturbation code removed from
  AIMNet2Result.cpp. TOML params removed.

### PBC + XTC infrastructure — PORTED
xdrfile (static lib), xtc_reader.h, pbc_whole.h all from fes-sampler.
Built and linked. Used by GromacsFrameHandler.

### Streaming trajectory processing — WORKING
GromacsProtein + GromacsFrameHandler + GromacsFinalResult.
- Free-standing conformations (public constructor, const Protein*)
- NOT in the Protein's conformations_ vector
- All 14 calculators run per frame (including SasaResult)
- 36 NPY arrays written per frame
- Conformation freed after each frame
- Tested: 10 XTC frames of 1Q8K_10023 (4876 atoms), PASSED

### SasaResult — NEW CALCULATOR
Per-atom Shrake-Rupley SASA. Bondi radii, Fibonacci sphere (92 pts),
TOML params (sasa_probe_radius, sasa_n_points). GeometryChoice recorded.
Wired into OperationRunner. Writes atom_sasa.npy.

### DsspResult extension — NEW OUTPUT
Was: 1 file (dssp_backbone.npy, 5 cols). Now: 4 files:
- dssp_backbone.npy (N, 5) — unchanged
- dssp_ss8.npy (N, 8) — full 8-class SS one-hot
- dssp_hbond_energy.npy (N, 4) — H-bond energies
- dssp_chi.npy (N, 12) — chi1-4 cos/sin/exists

### libtorch filesystem clash — IDENTIFIED
libtorch ships its own std::filesystem that shadows the system one.
fs::remove_all resolves to a broken stub. Use POSIX system() calls
for filesystem ops in code linked against libtorch.

---

## What to do NEXT

### 1. SDK updates (Python) — DONE (2026-04-12)
New arrays registered in _catalog.py: atom_sasa, dssp_ss8,
dssp_hbond_energy, dssp_chi. aimnet2_charge_sensitivity kept as
LEGACY entry for backward compat with old extractions.
Still needed: update python/API.md with the new arrays, run
`python -m pytest python/tests/` to verify.

### 2. HighFive integration
Add HighFive (header-only C++ HDF5 wrapper) to extern/.
Write the .h5 master file in GromacsFinalResult::Finalize.

### 3. Accumulation in GromacsFrameHandler
Add Welford accumulators for per-atom quantities:
- aimnet2_charge → ensemble_charges_{mean,var}.npy
- Kernel T2 tensors → ensemble_kernel_{mean,var}.npy
- aim embedding → ensemble_aim_{mean,var}.npy
- SASA → ensemble_sasa_{mean,var}.npy
- DSSP counts → ensemble_dssp_histogram.npy

Ridge-informed frame selection using calibration coefficients.

### 4. Constitution + Object Model docs
Update to reflect:
- charge_sensitivity no longer on Atom
- SasaResult as new calculator
- DsspResult extended output
- GromacsProtein pattern (not EnsembleConformation)

### 5. Full fleet extraction
685 proteins x ~500 frames on 4 machines.

---

## Key design decisions (settled 2026-04-11/12)

### Free-standing conformations
Streaming frames are created via the public ProteinConformation
constructor with a const Protein* back-pointer. They are never
added to the Protein's conformations_ vector. This was verified:
zero hidden coupling in all calculators and loaders. The Protein
stays immutable. No ConformationList, no invalidation, no records.

### GromacsProtein pattern (not evaluator ABC)
No abstract observer/evaluator interface. Three concrete classes
married to GROMACS. If another trajectory format appears, write
its own cousin, not a framework.

### {protein_id}.h5 master file
HDF5 containing topology + EnrichmentResult classifications from
conformation 0. Bond graph from Protein, not a ConformationResult.
SDK reads both .h5 and NPY, presents one typed object.

### Conformation 0
0 is 0. First GROMACS frame. No selection. Always valid, always in
memory, feeds the .h5. For --pdb/--mutant modes: untouched.

### No disk during streaming
Don't write every frame to disk. Accumulate what the handler needs
from the live conformation, then free it. Winners are re-loaded
from XTC at the end.

---

## Key references

| What | Where |
|------|-------|
| Ensemble design decisions | memory: project_ensemble_design_decisions.md |
| Ensemble model spec | spec/ENSEMBLE_MODEL.md (updated) |
| Streaming test | tests/test_gromacs_streaming.cpp |
| GromacsProtein | src/GromacsProtein.h/.cpp |
| GromacsFrameHandler | src/GromacsFrameHandler.h/.cpp |
| GromacsFinalResult | src/GromacsFinalResult.h/.cpp |
| SasaResult | src/SasaResult.h/.cpp |
| AIMNet2Result | src/AIMNet2Result.h/.cpp |
| DsspResult (extended) | src/DsspResult.cpp (WriteFeatures) |
| XTC reader | src/xtc_reader.h |
| PBC fixer | src/pbc_whole.h |
| Calculator params | data/calculator_params.toml |
| SDK catalog | python/nmr_extract/_catalog.py |
| Marching orders (historical) | learn/docs/cpp_marching_orders_2026-04-11.md |
| Autograd test | learn/src/actual_physics/test_aimnet2_autograd.py |
| Test protein | tests/data/fleet_test_large/ (1Q8K_10023) |
