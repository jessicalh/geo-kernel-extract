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

## Use Case D: Feature extraction for ML

I am providing geometry kernel output to a high-level e3nn NMR
prediction (or other physical prediction model) and am acting as a
feature extractor. I have an optional DFT, a set of optional PDBs
which may or may not be protonated, and an optional set of extracted
poses each with probability (see tests).

If I have a DFT I include it in the path. If I have PDBs that are
not protonated I protonate them with the method specified on the
command line. If I have a bunch of trajectory poses I run the
calculators on all of them. Then I write it all to disk.

**Input:** mix of sources, specified on command line
**Output:** everything written to disk

For each input item:
- If DFT: use case B path (load ORCA, pipeline, DFT extract), write
- If PDB needing protonation: use case A path (protonate with
  CLI-specified method, pipeline), write
- If PDB already protonated: load, assign charges, pipeline, write
- If trajectory (GROMACS): load protein once, run pipeline on each
  frame, write

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

---

## Geometry-only mode (--no-mopac --no-apbs)

For large ensembles (e.g. 5000 MD frames at 4000 atoms),
MOPAC (~10 min/frame) and APBS (~15-30 s/frame) dominate wall
time. The `--no-mopac` and `--no-apbs` flags skip these external
solvers while retaining all purely geometric calculators.

**What runs:** Geometry, SpatialIndex, Enrichment, DSSP,
BiotSavart, HaighMallion, McConnell, RingSusceptibility,
PiQuadrupole, Dispersion, Coulomb (topology charges), HBond.
8 calculators, sub-second per frame.

**What is skipped:** MopacResult, MopacCoulombResult,
MopacMcConnellResult, ApbsFieldResult. The Python SDK handles
absent MOPAC/APBS arrays — those groups are `None` on the
Protein dataclass.

**Typical invocation (fleet):**

    nmr_extract --fleet --tpr topology.tpr --poses /path/to/poses \
                --no-mopac --no-apbs --output /path/to/output

Both flags are available for all modes (PDB, ORCA, mutant, fleet)
and in both the CLI and UI.

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
