# Trajectory-Driven Extraction: Full-System GROMACS Path

**Status: implemented (2026-04-13).  All solvent calculators built and
wired.  See OperationRunner.cpp for dispatch order.**

## Context

The fleet has 685 proteins × full MD trajectories (CHARMM36m, PME,
explicit TIP3P water, ~150K atoms per system).  The pre-extracted
6-pose PDB path (fleet mode) works and must continue to work for
paper tests.  But the production run for the thesis points directly
at the trajectory with all atoms — protein, water, ions.

We have 30 days and 4 machines.  Frame selection is not settled.
The strategy is: extract everything we can per frame, write NPY,
let the calibration pipeline (per-element ridge) find the
dimensions.  Cost per frame may increase; that trades against
being more selective about which frames to extract.
(IMPLEMENTED 2026-04-12: GromacsProtein::SelectFrames() uses
Welford variance stats from the scan pass to select up to N
representative frames. See src/GromacsProtein.h.)

## What the full-system trajectory gives us (per protein atom)

### Already built (geometry-only, protein atoms)

8 classical calculators (ring current, bond anisotropy, EFG,
quadrupole, dispersion, H-bond, ring susceptibility), SASA,
DSSP (8-class SS, chi1-4, H-bond energies).
55 kernels.  R²=0.818 on mutants.  3 identified T2 dimensions.

### Already built (charges from TPR, no water needed)

- APBS solvated E-field + EFG (continuum Poisson-Boltzmann, 4s)
- Coulomb vacuum E-field + EFG (ff charges, 25s — retired from
  canonical path, available via --no-coulomb inversion)
- AIMNet2 charges + EFG (CUDA, geometry-only neural net)

### GromacsEnergyResult

- Per-frame Coulomb-SR, Coulomb-recip,
  LJ-SR, potential, temperature, pressure, volume from .edr.
  Aggregate (not per-atom) — useful for frame characterisation
  and selection.

### Built: explicit solvent calculators (full-system .xtc)

All built and wired into OperationRunner.  Dispatch is conditional
on `opts.solvent && !opts.solvent->Empty()`.

#### 1. WaterFieldResult (DONE)

Per-atom Coulomb E-field + EFG from TIP3P water charges within
cutoff.  KernelFilterSet: MinDistanceFilter.  5 NPY arrays.
TOML: water_efield_cutoff, water_first_shell_cutoff,
water_second_shell_cutoff, singularity_guard_distance.

#### 2. HydrationShellResult (DONE)

Per-atom hydration shell geometry using protein-COM reference.
Half-shell asymmetry, dipole orientation, nearest ion.  1 NPY.
TOML: water_first_shell_cutoff (shared), hydration_ion_cutoff.

#### 3. HydrationGeometryResult (DONE, 2026-04-13)

Per-atom water polarisation using SASA-derived surface normal as
reference frame (replaces COM direction).  Writes
water_polarization.npy (N, 10): net dipole vector(3), surface
normal(3), half-shell asymmetry, dipole alignment, coherence,
first-shell count.  Depends on SasaResult.  No KernelFilterSet.

#### 4. EeqResult (DONE, 2026-04-13)

D4 extended electronegativity equilibration (Caldeweyher 2019).
Geometry-dependent charges from N×N Cholesky solve.  Writes
eeq_charges.npy + eeq_cn.npy.  Pure C++/Eigen.  D4 parameters
in PhysicalConstants.h.  TOML: eeq_total_charge, eeq_cn_steepness,
eeq_cn_cutoff, eeq_charge_clamp.

#### IonFieldResult (not built — subsumed)

Ion E-field contribution is handled by nearest_ion_distance and
nearest_ion_charge in HydrationShellResult.  A separate IonFieldResult
was not needed — ions are too few (~10 per system) for a standalone
Coulomb sum to add meaningful signal beyond nearest-ion distance.

### From trajectory dynamics (multi-frame, accumulated)

Implemented via GromacsProteinAtom Welford accumulators (48 tracked
quantities).  See spec/ENSEMBLE_MODEL.md for the full inventory.

Includes: per-atom RMSF, kernel T0+|T2| variance (6 ring current,
4 McConnell, 2 PiQuad, 2 dispersion), AIMNet2 charge variance,
APBS field variance, water E-field mean+variance, SASA, hydration,
DSSP phi/psi/chi/SS, bond angles, 4 frame-to-frame delta trackers.

**Not yet integrated:** HydrationGeometryResult (water polarisation)
and EeqResult (EEQ charges).  See spec/OUTSTANDING_GROMACS_PATH.md.

## Data requirement: full-system .xtc

**RESOLVED.** Fleet re-run uses `compressed-x-grps = System` in
the .mdp.  All atoms (protein + water + ions) saved to .xtc.
Restart from equilibrated .gro with converged PLUMED biasing.

## Implementation status

### Step 1: Full-system trajectory reader — DONE

FullSystemReader reads all atoms from .xtc.  TPR topology
identifies protein, water O/H, and ion atoms.  Protein positions
go to existing calculators; water/ion positions packaged as
SolventEnvironment and passed via RunOptions.solvent.

### Step 2: Solvent calculators — DONE

WaterFieldResult, HydrationShellResult, HydrationGeometryResult,
EeqResult as ConformationResult subclasses.  Wired into
OperationRunner.  Hard-fail policy: if solvent data is provided,
solvent calculators MUST succeed.

### Step 3: NPY output + SDK — DONE

74 arrays registered in _catalog.py.  WaterPolarizationGroup and
EeqGroup in the Python SDK.  Column layouts in python/API.md.

### Step 4: GromacsProtein accumulation — NEXT

Wire HydrationGeometryResult and EeqResult into GromacsFrameHandler
(Welford accumulators), H5 master file, WriteCatalog CSV, and SDK
trajectory loader.  See spec/OUTSTANDING_GROMACS_PATH.md.

### Step 5: Calibration

Run per-element ridge with expanded feature set on 723 mutations.
Compare R² with and without solvent + EEQ features.

## Timing estimates

At 4876 protein atoms and ~50,000 water atoms:
- WaterFieldResult with 15Å cutoff: ~500-2000 water neighbors per
  protein atom.  O(N_protein × k_water) ≈ 4876 × 1000 ≈ 5M pair
  evaluations.  With spatial indexing, estimated 5-15s per frame.
- HydrationShellResult: same neighbor search, lighter computation.
  Estimated 2-5s per frame.
- IonFieldResult: ~10 ions, trivial.

Total per frame with all calculators: ~30-40s (vs 21s now).
At 685 proteins × N frames per protein, 4 machines, 30 days:
if N=100 frames: 685 × 100 × 40s / 4 = ~19 days.  Feasible.
if N=500 frames: ~95 days.  Need frame selection.

## Relation to existing paths

The pre-extracted PDB fleet path (--fleet --tpr --poses) is
unchanged.  It reads protein-only PDBs, runs the existing
calculators, writes NPY.  This path is for paper tests and
backward compatibility.

The trajectory-driven path is new.  It reads the full-system
trajectory, runs all calculators (existing + solvent), writes
an expanded NPY set.  Frame selection: IMPLEMENTED 2026-04-12.
