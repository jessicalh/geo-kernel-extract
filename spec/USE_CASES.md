# Use Cases

Written down verbatim from discussion 2026-04-03. Agreed between
Jessica and Claude before any spec or implementation work begins.

---

## Use Case A: UI loads a bare PDB

I am in front of the UI. I pick a PDB. It is not protonated. We
protonate it, and assign formal charges, and run all the things
except DFT load. Then we safely expose all the calculation results
to the UI caller and optionally write it all to a directory.

**Input:** one unprotonated PDB, protonation method (or default)
**Output:** all calculation results exposed to UI caller, optional write to dir
**No DFT. One conformation.**

1. Load PDB -> Protein (heavy atoms, no H, variant_index all -1)
2. Protonate: predict pKa (PROPKA or specified method), add H atoms
   -> new Protein with correct atoms and variant_index
3. Assign charges from ff14SB (variant_index now set -> correct
   variant-specific charges)
4. Net charge computed from charge sum (should be integer)
5. Pipeline: foundation -> charges -> xTB -> APBS -> 8 calculators
6. Results exposed to UI caller
7. Optionally write to directory

---

## Use Case B: UI loads a single DFT result

I am in front of the UI and I pick a single DFT. We load it with
formal charges from DFT and run all the things including DFT
extract and expose the calculation results to the UI caller and
potentially write it all to a dir.

**Input:** ORCA run (XYZ + prmtop + NMR output), already protonated by tleap
**Output:** all calculation results + DFT tensors exposed to UI caller, optional write

1. Load ORCA -> Protein with H atoms, variant_index from AMBER labels
2. Charges from PrmtopChargeSource (authoritative)
3. Net charge from prmtop charge sum
4. Pipeline: foundation -> charges -> xTB -> APBS -> 8 calculators
5. Load OrcaShieldingResult (DFT tensors)
6. Results exposed to UI caller
7. Optionally write to directory

---

## Use Case C: CLI mutant pair for training data

I am creating a model for physics OR a model for helping calculators
refine themselves. I load a DFT mutant pair from the command line
and do all the things including DFT and save all the things to a
directory.

**DFTs are required for training models for internal use -- both
mutants.** For the other use cases they are optional but should be
used if present, always.

**Input:** WT + ALA ORCA runs (both with prmtop, XYZ, NMR output)
**Output:** everything written to directory

1. Load WT -> Protein, charges from prmtop, net charge
2. Pipeline on WT conformation
3. Load OrcaShieldingResult on WT
4. Load ALA -> separate Protein, charges from prmtop, net charge
5. Pipeline on ALA conformation
6. Load OrcaShieldingResult on ALA
7. MutationDeltaResult::Compute(wt_conf, ala_conf) -> attaches to WT
8. Write everything to directory

---

## Use Case D: Feature extraction for ML (static structures)

I am providing geometry kernel output to a high-level e3nn NMR
prediction (or other physical prediction model) and am acting as a
feature extractor. I have an optional DFT, a set of optional PDBs
which may or may not be protonated.

If I have a DFT I include it in the path. If I have PDBs that are
not protonated I protonate them with the method specified on the
command line. Then I write it all to disk.

**Input:** mix of static structure sources, specified on command line
**Output:** everything written to disk

For each input item:
- If DFT: use case B path (load ORCA, pipeline, DFT extract), write
- If PDB needing protonation: use case A path (protonate with
  CLI-specified method, pipeline), write
- If PDB already protonated: load, assign charges, pipeline, write

---

## Use Case E: Trajectory extraction (GROMACS ensemble)

I have a full-system GROMACS trajectory — protein + explicit water
+ ions — and I want to extract geometric kernels across all frames.

**Input:** protein directory containing md.tpr, md.xtc, md.edr (all required)
**Output:** per-frame NPY arrays, per-atom trajectory catalog (CSV),
H5 master file with Welford rollup statistics

    nmr_extract --trajectory /path/to/protein_dir \
                --no-mopac --no-coulomb --output /path/to/output

Two-pass architecture:
1. **Pass 1 (streaming):** GromacsFrameHandler reads all XTC frames,
   applies PBC fix, runs lightweight calculators per frame (classical
   kernels, AIMNet2, APBS, DSSP, SASA, WaterField, HydrationShell,
   HydrationGeometry, EEQ). Accumulators on GromacsProteinAtom track
   Welford statistics (mean + std) across frames. Frame-to-frame
   DeltaTrackers monitor charge fluctuation, SASA change, water
   exchange rate, SS transitions, chi rotamer flips.
2. **Pass 2 (selected frames):** extract selected representative
   frames with full calculators, write NPY per frame.

GromacsProtein writes the H5 master file (per-frame positions,
48-column Welford rollup, bond length statistics, frame times) and
CSV catalog. The Python SDK reads both via `load_trajectory()`.

This is CLI-only. The viewer returns immediately if trajectory mode
is dispatched — it is a batch operation, not interactive

---

## Use Case F: Analysis trajectory (exhaustive per-frame H5)

I have a GROMACS trajectory and I want to see EVERYTHING — every
calculator's output at every sampled frame, per-ring decompositions,
full SphericalTensors, enrichment flags, topology — in one
self-contained H5 file for R analysis, time series model training,
and GNN message design.

**Input:** protein directory containing md.tpr, md.xtc, md.edr (all required),
AIMNet2 model (required)
**Output:** `{protein_id}_analysis.h5` — exhaustive per-frame data,
PDB snapshots at ~1ns intervals (ORCA inputs)

    nmr_extract --trajectory --analysis /path/to/protein_dir \
                --aimnet2 model.jpt --output /path/to/output

Single-pass architecture with stride 2 (~625 of 1250 frames):
all calculators run every sampled frame (APBS, AIMNet2 mandatory,
MOPAC skipped). AnalysisWriter buffers per-frame data in memory,
writes the complete H5 at end. ~1.7 GB per protein at 4000 atoms.
PDB snapshots written at ~1ns intervals for DFT (ORCA) input.

The H5 is self-contained: typed topology (bonds, rings, enrichment),
per-frame physics groups (ring_current, efg, bond_aniso, quadrupole,
dispersion, hbond, sasa, water, charges), per-ring K=6 decomposition,
AIMNet2 256-dim embedding, per-frame ring geometry, positions.

See **spec/ANALYSIS_TRAJECTORY_2026-04-14.md** for the full design.

---

## Common observations

Effectively we are either a library that is run on something in a
UI when it is loaded (and the bar goes up while we work), or
something run on a command line over and over on DFT results,
things we might want to try new protonation on, or trajectories.

Our only inputs for each of these cases are in tests now so it
should not be hard to picture the full set.

We are nearly done. We should not screw up our object model -- it
should become cleaner in small ways. But in terms of paths through
what we do, that is what we do.

If we create a model that calculators can draw on, it will be its
own thing. This system can generate results for it but does not and
should not say how that might look.

### Two executables, one parser

Both `nmr_extract` and `nmr-viewer` parse command lines via
`JobSpec` (src/JobSpec.h). The viewer accepts all the same modes
as nmr_extract plus `--rest-port`. The difference:

- **nmr_extract:** headless, writes NPY and exits. MOPAC and
  Coulomb on by default. All 5 modes. The production tool.
- **nmr-viewer:** Qt/VTK GUI. MOPAC and Coulomb hardcoded off
  (interactive, not batch). `--output` optional — omit to just
  visualise. Trajectory mode returns immediately (batch only).
  REST API on port 9147 for programmatic control.

---

## Canonical extraction mode (--no-mopac --no-coulomb)

APBS (solvated Poisson-Boltzmann) is the preferred electrostatic
calculator. It is faster than vacuum Coulomb at N > 1000 atoms
(4s vs 25s at 4800 atoms — fixed 161^3 grid vs O(N*k) pairwise)
and produces solvated fields. Vacuum Coulomb is retained for
comparison tests only.

**What runs:** Geometry, SpatialIndex, Enrichment, DSSP,
BiotSavart, HaighMallion, McConnell, RingSusceptibility,
PiQuadrupole, Dispersion, APBS (solvated E-field + EFG),
HBond, SASA.  AIMNet2 if model provided.

**What is skipped:** MopacResult, MopacCoulombResult,
MopacMcConnellResult, CoulombResult. The Python SDK handles
absent MOPAC/Coulomb arrays — those groups are `None` on the
Protein dataclass.

**Typical invocation (trajectory):**

    nmr_extract --trajectory /path/to/protein_dir \
                --no-mopac --no-coulomb --output /path/to/output

**Typical invocation (single protein):**

    nmr_extract --pdb protein.pdb --no-mopac --no-coulomb --output /path/to/output

All skip flags are available for all modes (PDB, ORCA, mutant, trajectory)
and in the CLI. The viewer hardcodes MOPAC and Coulomb off; use `--no-apbs`
to skip APBS. See spec/TIMING_4876_ATOMS.md for measured per-calculator
times on 4876-atom protein.

---

## Corollaries

1. At least one protonation option and charge assignment mechanism
   must be working. PROPKA is an external Python binary call.

2. xTB and APBS must be working first-class citizens, not optional
   tools that silently fail.

3. There must be a CLI and UI "main" for the use cases above.

4. That main should not accumulate pass-to-pass cruft. The
   orthogonality we have now must be fully preserved.

5. DRY but KISS. No simplification that splits things up too much
   between UI and CLI.

6. This is a single milestone.
