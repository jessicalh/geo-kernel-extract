# Trajectory-Driven Extraction: Full-System GROMACS Path

**Status: design (2026-04-12).  GromacsEnergyResult (.edr reader) built
and wired.  Full-system calculators below are next.**

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

### Built this session

- GromacsEnergyResult: per-frame Coulomb-SR, Coulomb-recip,
  LJ-SR, potential, temperature, pressure, volume from .edr.
  Aggregate (not per-atom) — useful for frame characterisation
  and selection.

### To build: explicit solvent calculators (need full-system .xtc)

These require reading water and ion positions from the trajectory,
not just protein atoms.  The current XtcReader reads protein-only
.xtc (water stripped by harvest).  The trajectory-driven path reads
the original full-system .xtc.

#### 1. WaterFieldResult — per-atom explicit-solvent E-field and EFG

For each protein atom i, sum over water atoms j within cutoff R:

    E_i = Σ_j  q_j · r_ij / |r_ij|³
    V_ij = q_j · (3 r_ij r_ij^T / |r_ij|⁵  −  I / |r_ij|³)

Same Coulomb kernel as CoulombResult but with water charges
(TIP3P: O = -0.834e, H = +0.417e) instead of protein charges.
Cutoff ~10-15Å.  Spatial indexing over water atoms.

Output: per-atom E-field (Vec3), EFG (Mat3 + SphericalTensor),
decomposed into first-shell (< 3.5Å) and outer-shell contributions.

This is what APBS approximates with a smooth dielectric.  The
explicit field includes water orientation fluctuations, cavities,
bridging water, structural water — none of which the continuum
model captures.

#### 2. HydrationShellResult — per-atom water packing geometry

For each protein atom i, characterise the water environment:

- n_water_first_shell: count of water O within 3.5Å
- n_water_second_shell: count of water O within 3.5-5.5Å
- half_shell_asymmetry: fraction of first-shell waters on the
  solvent-exposed side vs buried side (UCSB half-shell method)
- mean_water_dipole_cos: average cos(angle) between water dipole
  and atom-water vector (orientation order parameter)
- nearest_ion_distance: distance to closest Na+/Cl- (if present)
- nearest_ion_charge: charge of that ion (+1 or -1)

Output: per-atom scalars (6+ columns).

This captures burial/exposure at atomic resolution — not a smooth
surface (SASA) but the actual water packing.  An atom in a water
cavity vs in a hydrophobic pocket vs at an interface have different
first-shell counts and dipole orientations.

#### 3. IonFieldResult — per-atom ion contribution

Same kernel as WaterFieldResult but over ion atoms only.
Ions carry integer charges; a nearby Na+ creates a ~10× stronger
field than a water molecule at the same distance.

Output: per-atom E-field from ions, nearest-ion distance and charge.

### From trajectory dynamics (multi-frame, accumulated)

These need the accumulator / Welford pattern across frames:

- Per-atom RMSF (positional fluctuation)
- Per-atom kernel variance (how much each kernel changes across frames)
- Per-atom water residence time (how long first-shell waters persist)
- Per-atom charge variance (AIMNet2 charge across conformations)

## Data requirement: full-system .xtc

The original fleet runs used `compressed-x-grps = Protein` in the
.mdp, so only protein coordinates were saved to .xtc. Water and
ion positions were discarded during simulation. The .gro files
have one full-system frame (final snapshot).

**Fix**: re-run production with `compressed-x-grps = System`.
The simulations are equilibrated; this is a restart, not a new
run. The .gro provides the restart state. PLUMED biasing is
converged. Same physics, just saving all atoms this time.

The crashed fleet machine needs recovery first. The surviving
machines can restart immediately with the .mdp change.

## Implementation path

### Step 1: Full-system trajectory reader

Extend or replace XtcReader to read all atoms from .xtc.  TPR
provides topology: which atoms are protein, which are water O/H,
which are ions.  Atom selection from TPR atom types.

The reader hands protein positions to the existing calculators
AND water/ion positions to the new solvent calculators.

### Step 2: Solvent calculators

WaterFieldResult, HydrationShellResult, IonFieldResult as
ConformationResult subclasses.  Each takes the water/ion
positions as input (passed via RunOptions or a SolventEnvironment
object).

### Step 3: NPY output + SDK

Register all new arrays in _catalog.py and wire into load().
Column layout documented in API.md.

### Step 4: Calibration

Run the existing learn pipeline (per-element ridge) with the
expanded feature set.  Compare R² with and without solvent
features.  If new dimensions appear in T2 analysis, that is a
thesis finding.

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
