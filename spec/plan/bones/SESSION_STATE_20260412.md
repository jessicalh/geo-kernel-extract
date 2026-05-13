# Session State: 2026-04-12 morning

## What was built this session

### Code (all compiles, library + executables build clean)

1. **--no-coulomb flag**: JobSpec.h/cpp, OperationRunner.h/cpp, nmr_extract.cpp.
   Coulomb is 25s at 4876 atoms. APBS is 4s. APBS is canonical.
   Verified: 21s/pose with --no-coulomb vs 46s with Coulomb.

2. **GromacsEnergyResult**: reads .edr files, extracts per-frame
   Coulomb-SR, Coulomb-recip, LJ-SR, potential, temperature, pressure,
   volume. Tested on 1Q8K_10023 walker_0 — Coul_total=-2.38M kJ/mol.
   Wired into OperationRunner (fires when edr_path is set).

3. **SolventEnvironment + FullSystemReader**: data struct for water
   molecules (O/H positions + TIP3P charges) and ions. Reader parses
   TPR topology to identify protein/water/ion atom ranges, splits
   full-system XTC frames.

4. **WaterFieldResult**: per-atom E-field and EFG from explicit water
   charges within 15A cutoff. First-shell (3.5A) decomposition.
   Shell counts. Same Coulomb kernel as CoulombResult applied to
   water charges. Writes 5 NPY files.

5. **HydrationShellResult**: per-atom half-shell asymmetry (UCSB
   method), water dipole orientation order parameter, nearest ion
   distance and charge. Writes 1 NPY file (N, 4).

6. **ConformationAtom fields**: 14 new fields for water E-field/EFG,
   shell counts, hydration geometry, ion proximity. All per-conformation
   per-atom. Nothing on Protein.

7. **SDK fixes from earlier**: _protein.py wired dssp_ss8, dssp_hbond_energy,
   dssp_chi, sasa. _catalog.py corrected required flags. API.md updated.
   Tests pass (62/62).

### Docs

- spec/TIMING_4876_ATOMS.md — measured per-calculator times
- spec/TRAJECTORY_EXTRACTION.md — full design for trajectory path
- spec/USE_CASES.md — updated canonical extraction mode
- OBJECT_MODEL.md, PATTERNS.md — "update pending" notes
- spec/EXTRACTION_SDK.md — ring_contributions 57→59 fixed
- spec/TEST_FRAMEWORK.md — SasaResult + AIMNet2Result added to pipeline order
- doc/ARCHITECTURE.md — line number fix

## Late session additions

- **SASA on ConformationAtom**: atom_sasa field, SasaResult writes
  to it during Compute. Private vector killed. One path.
- **1ZR7_6721 full-system test data**: 479 protein atoms, 3525 water
  molecules, 25 ions, ~1001 frames per walker, 5 walkers.
  First explicit solvent extraction PASSED — 7s/frame, all 17
  calculators including WaterFieldResult + HydrationShellResult.
- **Literature search**: per-atom explicit water features for NMR
  prediction is novel. Nobody has done this.
- **Harvest bug found**: compressed-x-grps = Protein in prod.mdp
  stripped water from .xtc. Fleet re-run needed with System.
  Runbook at fleet-manager-tasks/RERUN_FULL_SYSTEM.md.
- **Disk budget**: 22 MB/frame at 4876 atoms (ring_contributions
  dominates). 5000 frames × large protein = 116 GB. Scout mode
  needed for frame selection before full extraction.

## Next session: --scout mode in nmr_extract

The scout is a run mode, not a separate tool. Same JobSpec, same
OperationRunner. It:
1. Reads full-system XTC + TPR
2. Streams all frames with minimal pipeline (skip APBS/Coulomb/MOPAC)
3. Per atom: Welford accumulators for water shell count, water
   E-field magnitude, SASA
4. Writes one CSV catalog: per-atom statistics with min/max frame
   indices for targeted frame selection
5. User inspects catalog in R, picks 50 frames for DFT
6. Then runs full extraction on those 50 frames only

Frame selection strategy: 40 by water diversity (atoms with biggest
solvation changes, frames at their extremes) + 10 by Boltzmann
weight. COLVAR has RMSD/Rg/rbias for free. .edr has Coulomb/LJ
energies for free.

## KNOWN BUG: --trajectory atom count mismatch

BuildFromTrajectory uses BuildFromProtonatedPdb to construct the
Protein, which drops atoms the TPR has (479 in TPR, 461 from PDB
builder). ScanFrame creates a ProteinConformation with 479 positions
from the XTC but the Protein has 461 atoms. Segfault.

**Fix**: BuildFromTrajectory should use GromacsEnsembleLoader's
TPR-based protein construction (which already handles CHARMM36m
naming, all atom types). The TPR is authoritative, not the PDB.
This is the same builder the fleet path uses. The reference PDB
should not be the protein construction source for the trajectory
path.

## What is NOT done

### Solvent calculators need full-system .xtc

The harvest script stripped water from .xtc files. The .gro has all
154,944 atoms (water + ions + protein), and the .tpr has the full
topology. But the .xtc is protein-only (4876 atoms).

**Fix**: on a fleet machine, re-extract with water:
```
gmx trjconv -f /path/to/original/md.xtc -s md.tpr -o full_system.xtc -pbc whole
```
Select "System" (not "Protein") when prompted.

Or fix the harvest script to not strip water.

### Test cleanup (tasks 3, 4, 5 from original list)

- Terminology: "streaming" → "trajectory-based" not yet done in code/docs
- Test restructuring for canonical fast path not yet done
- gromacs_streaming_tests updated to use APBS + GromacsEnergy + skip_coulomb

### SDK registration for new NPY files

New arrays from WaterFieldResult (water_efield, water_efg, water_efg_first,
water_efield_first, water_shell_counts) and HydrationShellResult
(hydration_shell) and GromacsEnergyResult (gromacs_energy) need entries
in _catalog.py and _protein.py. Deferred until we have test data to
validate against.

## Key findings

- APBS is 6x faster than Coulomb at 4876 atoms (4s vs 25s)
- APBS + all calculators = 21s/pose (no Coulomb, no MOPAC)
- .edr has real PME electrostatic energies (-2.38M kJ/mol)
- .gro has full system (154,944 atoms) — water IS in the tree
- .xtc is protein-only — harvest bug stripped water
- learn/actual_physics R²=0.818, 3 identified T2 dimensions
- Solvent features (water E-field, hydration shell) could carry
  new dimensions not captured by geometry-only kernels
