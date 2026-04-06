# Builder and Pipeline Unification — Tentative Spec

**Status: TENTATIVE — for discussion before implementation.**
**References: spec/USE_CASES.md (agreed use cases)**

---

## 1. What this milestone delivers

A single path from any input source to fully computed results:

```
Build(source) → BuildResult{protein, charges, net_charge}
    ↓
OperationRunner::Run(conf, opts) → all results on conformation
    ↓
expose to caller / write to disk
```

All four use cases (A-D from USE_CASES.md) use this path.
CLI and UI call the same library functions.
xTB and APBS are first-class citizens, not optional.
Protonation via linked reducelib (Richardson lab C++ library).
At least one protonation path works end-to-end for bare PDBs.

---

## 2. BuildResult — the common output of all builders

Every loading path produces the same thing:

```cpp
struct BuildResult {
    std::unique_ptr<Protein> protein;       // fully constructed, topology resolved
    std::unique_ptr<ChargeSource> charges;  // matched to this protein's atoms
    int net_charge = 0;                     // formal charge (integer)
    bool ok = false;
    std::string error;
};
```

The Protein inside BuildResult has:
- All atoms present (including H if protonated)
- variant_index set on titratable residues
- Symbolic topology resolved (backbone cache, rings)
- Covalent topology resolved (bonds, H parents)
- At least one conformation created

The ChargeSource matches the protein's atoms. The net_charge is
the integer formal charge determined by the protonation state.

---

## 3. What is CONSTANT across all builders

Every builder does these steps, in this order:

```
1. Read source → raw atoms with source-specific names
2. Translate names to canonical (NamingRegistry, declared ToolContext)
3. Add atoms and residues to Protein (canonical names)
4. Ensure protonation is determined (H atoms present, variant_index set)
5. Symbolic topology: CacheResidueBackboneIndices + DetectAromaticRings
6. Geometric topology: CovalentTopology::Resolve from positions
7. Create conformation(s) on Protein
8. Produce ChargeSource and compute net_charge
```

Steps 5-6 are FinalizeConstruction (unchanged).
Step 8 is new — loaders currently don't all do this.

---

## 4. What VARIES per builder

| Step | PDB builder | ORCA builder | GROMACS builder |
|------|-------------|--------------|-----------------|
| 1. Parser | cifpp | prmtop + XYZ | libgromacs TPR + pose PDBs |
| 2. ToolContext | Standard | Amber | Charmm |
| 4. Protonation | OpenBabel AddHydrogens(pH) | Already done by tleap | Already done by pdb2gmx |
| 7. Conformations | AddConformation (1) | AddPrediction (1) | AddMDFrame (N) |
| 8. ChargeSource | ParamFileChargeSource (ff14SB) | PrmtopChargeSource | PreloadedChargeSource |

---

## 5. The three builders (refactored from existing loaders)

### 5a. PDB builder (from PdbFileReader)

Current: `LoadProtein(path) → LoadResult{protein, ok, error}`
New: `BuildFromPdb(path, pH) → BuildResult`

New step between reading and FinalizeConstruction:
- If pH specified (non-NaN): call OpenBabel AddHydrogens(pH) on
  the heavy-atom structure
- Extract new atom positions (heavy + H) from OBMol
- H atoms get canonical names (OpenBabel uses standard PDB naming)
- Protonation variant detected from which H atoms are present
  (existing code in DetectAromaticRings and CacheResidueBackboneIndices)

After FinalizeConstruction:
- Create ParamFileChargeSource (ff14SB param file path from ToolPaths)
- Compute net_charge: either from ProtonationState::NetChargeForProtein()
  or from rounding the charge sum

If pH is NaN (no protonation requested): load as-is. Charges will
have variant_index=-1 defaults. Caller accepts reduced physics
fidelity.

**Protonation via reducelib (linked C++ library):**
reduce (Richardson lab, Duke University) is a C++ library for adding
hydrogens to PDB files. It handles HIS tautomer assignment, flip
optimization, and produces PDB-standard H atom names (H, HA, HD1,
HE2, etc.) that our existing protonation detection code recognises
without modification.

Source: https://github.com/rlabduke/reduce (C++11, CMake, MIT-like)
Already installed as binary at /home/jessicalh/amber24/bin/reduce.
Source cloned to /home/jessicalh/builds/reduce-src.
libreducelib.a builds clean on ARM64 GCC 13.

Integration: link reducelib as a static library. One thin wrapper
in our codebase:

```cpp
// The ONLY place reduce internals are touched.
// PDB string in → protonated PDB string out. One boundary crossing.
std::string ProtonateWithReduce(const std::string& pdb_content);
```

Internally: set reduce globals for BUILD mode, call inputModels(),
scanAndGroupRecords(), reduceList(), optimize(),
outputRecords_all_string(). The het dictionary
(reduce_wwPDB_het_dict.txt) ships with the reduce source.

After ProtonateWithReduce returns, our existing PdbFileReader parses
the protonated PDB string. DetectAromaticRings reads H atom names
(HD1/HE2) to determine HIS tautomer → variant_index. No new string
processing in our codebase. The naming boundary is unchanged.

PROPKA remains available as an alternative protonation decision tool
(external Python binary). The builder could take a ProtonationMethod
enum: Reduce, PROPKA, AlreadyProtonated.

### 5b. ORCA builder (from OrcaRunLoader)

Current: `LoadOrcaRun(files) → OrcaLoadResult{protein, ok, error}`
New: `BuildFromOrca(files) → BuildResult`

Changes:
- Create PrmtopChargeSource from files.prmtop_path
- Compute net_charge by summing prmtop charges (round to integer)
- Return both in BuildResult

The ORCA NMR path (for DFT loading) is NOT in BuildResult — it's
caller metadata. The caller passes it to the pipeline or loads
OrcaShieldingResult separately.

### 5c. GROMACS builder (from GromacsEnsembleLoader)

Current: `LoadFleetEnsemble(paths) → FleetLoadResult{protein, charges, net_charge, ok, error}`
New: `BuildFromGromacs(paths) → BuildResult`

Changes:
- Wrap existing charges vector in PreloadedChargeSource
- net_charge already computed — copy to BuildResult

This loader is closest to the target already. FleetLoadResult
already has charges and net_charge.

---

## 6. Formal charges — explicit, not incidental

The formal charge per residue is determined by protonation state:

| Residue | Default state | Formal charge | Protonated variant | Formal charge |
|---------|--------------|---------------|-------------------|---------------|
| LYS | protonated (NH3+) | +1 | LYN (neutral) | 0 |
| ARG | protonated (guanidinium+) | +1 | ARN (neutral) | 0 |
| HIS | HIE (neutral) | 0 | HIP (doubly protonated) | +1 |
| ASP | deprotonated (COO-) | -1 | ASH (protonated) | 0 |
| GLU | deprotonated (COO-) | -1 | GLH (protonated) | 0 |
| CYS | thiol (SH) | 0 | CYX (disulfide) | 0 |
| N-term | protonated (NH3+) | +1 | | |
| C-term | deprotonated (COO-) | -1 | | |

Net formal charge = sum over all residues. This is an integer by
definition. The force field partial charges sum to this integer.

The builder computes net_charge from the protonation state. xTB
receives it as --chrg. APBS receives per-atom partial charges that
sum to it. These are not independent — they are the same physical
quantity expressed differently.

ProtonationState::NetChargeForProtein() already exists and computes
this correctly. For paths where ProtonationState is not explicitly
constructed (ORCA, GROMACS), net_charge is computed by rounding the
partial charge sum. Both methods must agree.

---

## 7. OperationRunner — the single home for all ordered sequences

### Design principle

```cpp
// OperationRunner: the ONE home for all ordered sequences.
//
// THIS IS ONLY FOR ORDER. The conformation is the buffer where
// results accumulate. The runner does not hold state, does not
// cache intermediate results, does not decide what to compute.
// It sequences operations. If a step's prerequisites are on the
// conformation, the step runs. If not, it is skipped. Nothing
// runs backwards. Nothing runs twice. These are the ONLY
// variations — the presence or absence of upstream results.
//
// Every use case (UI, CLI, batch, training) calls one of
// these methods. The conformation arrives with its protein
// identity determined. The runner fills it with computed results
// in dependency order.
//
// All ordered sequences live here. Future agents add new
// sequences here, not in new files.
```

OperationRunner replaces BOTH Pipeline and CalculationRunner.
They were two halves of the same thing with duplicated order.

### RunOptions

```cpp
struct RunOptions {
    // Charges — REQUIRED for real physics.
    // Null = geometry-only mode (no Coulomb, no APBS, no xTB).
    const ChargeSource* charge_source = nullptr;
    int net_charge = 0;

    // xTB: runs when charges are available and atom count <= max.
    // Skip with logged warning above this threshold.
    size_t xtb_max_atoms = 1000;

    // DSSP: skip for structures without backbone.
    bool skip_dssp = false;

    // DFT: load ORCA shielding tensors after calculators.
    // Empty = no DFT.
    std::string orca_nmr_path;
};
```

### Sequences

**Run (the standard sequence):**

Every use case calls this for each conformation.

```
Tier 0: GeometryResult, SpatialIndexResult, EnrichmentResult
        DsspResult (unless skip_dssp)
        ChargeAssignmentResult (from charge_source)

Tier 0.5 (when charges available):
        xTB (if atom count <= xtb_max_atoms)
        APBS (failure is a hard error — charges were provided)

Tier 1: All 8 classical calculators
        Coulomb (requires charges)
        HBond (requires DSSP)

Tier 2 (optional):
        OrcaShieldingResult (if orca_nmr_path provided)
```

xTB and APBS are not optional flags. They run when charges exist.
xTB has a size gate. APBS failure stops with an error.

**RunMutantComparison:**

Calls Run on WT conformation, then Run on ALA conformation,
then attaches MutationDeltaResult to WT. This is use case C.
Same order within each run.

**RunEnsemble:**

Calls Run on each conformation of a protein (MD frames).
This is use case D with trajectories.

### What gets deleted

CalculationRunner.h/.cpp and Pipeline.h/.cpp are replaced by
OperationRunner.h/.cpp. All batch tests, fleet tests, and
full pipeline tests call OperationRunner.

The four batch test files' manual AttachResult chains are deleted.
They call OperationRunner::Run() instead.

---

## 8. Entry points — CLI and UI main

### Shared library API (both CLI and UI call these)

```cpp
// Build from any source — returns fully constructed protein + charges
BuildResult BuildFromPdb(const std::string& pdb_path, double pH = 7.0);
BuildResult BuildFromOrca(const OrcaRunFiles& files);
BuildResult BuildFromGromacs(const FleetPaths& paths);

// Run all calculators on one conformation
std::vector<std::string> RunClassicalCalculators(
    ProteinConformation& conf, const PipelineOptions& opts);

// Write all results to disk
int WriteAllFeatures(const ProteinConformation& conf,
                     const std::string& output_dir);
```

### CLI main (new executable: nmr_extract)

```
nmr_extract --pdb 1UBQ.pdb --pH 7.0 --output results/1UBQ/
nmr_extract --orca wt_dir/ --output results/wt/
nmr_extract --mutant wt_dir/ ala_dir/ --output results/delta/
nmr_extract --fleet paths.txt --output results/fleet/
```

Thin wrapper: parse args → build → pipeline → write.

### UI integration

The UI calls the same library functions. BuildResult gives it
a Protein with results on conformations. The UI reads results
via conf.Result<T>() and conf.HasResult<T>(). No special UI API
needed beyond what the library already provides.

---

## 9. What changes in existing code

### Files modified (surgery):
- src/PdbFileReader.cpp — add OpenBabel protonation, return BuildResult
- src/PdbFileReader.h — new signature
- src/OrcaRunLoader.cpp — add charge source + net_charge, return BuildResult
- src/OrcaRunLoader.h — new signature
- src/GromacsEnsembleLoader.cpp — wrap charges in ChargeSource, return BuildResult
- src/GromacsEnsembleLoader.h — new signature
- src/Pipeline.h — refined PipelineOptions, mandatory xTB/APBS
- src/Pipeline.cpp — add xTB/APBS/DFT to pipeline sequence
- Every test that calls LoadProtein/LoadOrcaRun/LoadFleetEnsemble
  (mechanical: .protein → .protein, add .charges/.net_charge)

### Files deleted:
- src/CalculationRunner.h
- src/CalculationRunner.cpp
- src/Pipeline.h
- src/Pipeline.cpp

### Files created:
- src/BuildResult.h — the common return type
- src/OperationRunner.h — the single home for all ordered sequences
- src/OperationRunner.cpp
- src/ReduceProtonation.h — thin wrapper around reducelib
- src/ReduceProtonation.cpp
- src/nmr_extract.cpp — CLI main (new executable)

### Files NOT changed:
- All calculator implementations (unchanged)
- ConformationResult framework (unchanged)
- ConformationAtom fields (unchanged)
- ProteinConformation (unchanged)
- Protein (unchanged — FinalizeConstruction stays)
- All ring types, bond types, atom roles (unchanged)

---

## 10. Verification strategy

### Before any changes:
- Capture current test pass count (283)
- Baseline NPY arrays exist (baseline_features/P84477/)
- Baseline batch outputs exist (baseline_batch_*.txt)

### After each step:
- All 283 tests pass
- NPY binary diff against baselines (physics unchanged)
- Batch test results match baselines

### After all changes:
- New tests for:
  - BuildFromPdb with protonation (H atoms added, variant_index set)
  - BuildFromPdb charges + net_charge correct
  - xTB runs on ORCA-sized proteins through Pipeline
  - APBS runs through Pipeline
  - CLI main produces correct output for each use case
- Batch tests run through Pipeline (not CalculationRunner)
- xTB and APBS assertions in batch tests

---

## 11. Implementation order (tentative)

1. Define BuildResult in src/BuildResult.h
2. Refactor GromacsEnsembleLoader → BuildFromGromacs (closest to target)
3. Refactor OrcaRunLoader → BuildFromOrca (add charges + net_charge)
4. Refactor PdbFileReader → BuildFromPdb (add OpenBabel protonation)
5. Refine PipelineOptions (mandatory xTB/APBS, net_charge, orca path)
6. Migrate batch tests from CalculationRunner to Pipeline
7. Delete CalculationRunner
8. Create CLI main (nmr_extract)
9. Final verification pass

Each step has passing tests before proceeding to the next.

---

## Resolved questions

1. **Protonation tool**: reduce (Richardson lab C++ library), linked
   as reducelib. Produces PDB-standard H naming. Tested on ARM64.
   No OpenBabel for protonation (its AddHydrogens names all H atoms
   just "H" — unusable). PROPKA available as alternative.

2. **Runner architecture**: OperationRunner is the ONE home for all
   ordered sequences. Not split into Pipeline + CalculationRunner.
   Multiple methods (Run, RunMutantComparison, RunEnsemble) but
   one file, one class, same discipline.

3. **Builder as free functions**: BuildFromPdb, BuildFromOrca,
   BuildFromGromacs. No ABC needed — callers always know their
   source type.

## Open questions

1. **APBS hard failure**: if APBS fails on a very large protein
   (grid memory), should the runner abort or continue with a
   warning? Proposed: hard error when charges were provided.

2. **xTB threshold**: 1000 atoms is a guess. Need empirical test
   to find the actual GFN2-xTB limit on ARM64 with 128GB RAM.

3. **reduce globals**: reducelib uses C-style globals for config.
   Our wrapper must set them before each call. Thread-unsafe but
   we are single-threaded per protein. Acceptable?

4. **reduce het dictionary**: reduce_wwPDB_het_dict.txt ships with
   the source (~500KB). Ship as data file in our repo, or reference
   from the installed amber24 path?
