# Loading Implementation Plan

**Date**: 2026-04-03
**Status**: CovalentTopology extracted and verified. This document
describes the remaining work: conformation types, loaders, fleet loading.

## What is done

### CovalentTopology (commit 9497c6a)

The geometry→topology boundary is a type. `CovalentTopology::Resolve()`
takes atoms, rings, residues, and positions. Returns bonds, per-atom
bond indices, and H parent assignments. Protein delegates bond access
through it. All 32 NPY feature files byte-identical to pre-refactor
baseline. 273 tests pass.

### Feature output (commit abf9ab9)

`WriteFeatures()` virtual method on each ConformationResult subclass.
All 8 calculators + ORCA write full SphericalTensor data to NPY files.
`ConformationResult::WriteAllFeatures()` traverses the conformation's
accumulated results and calls each one. This is the system's output
mechanism for the e3nn model.

### Baselines

Stashed at `baseline_stash/`:
- `baseline_features/P84477/`: 32 NPY arrays (351 atoms, all calculators)
- `baseline_features_log.txt`: operation logs with filter rejection counts
- `baseline_batch_*.txt`: 4 batch tests (720-723 proteins each)

After any change: `cmp baseline_stash/baseline_features/P84477/*.npy baseline_features/P84477/*.npy` must show 0 differences.

## What is NOT done

### 1. Conformation types (trivial)

Add to ProteinConformation.h:

```cpp
class PredictionConformation : public ProteinConformation {
public:
    PredictionConformation(const Protein* protein,
                           std::vector<Vec3> positions,
                           std::string method,
                           double confidence = std::nan(""));
    const std::string& PredictionMethod() const;
    double Confidence() const;
private:
    std::string method_;
    double confidence_;
};

class MDFrameConformation : public ProteinConformation {
public:
    MDFrameConformation(const Protein* protein,
                        std::vector<Vec3> positions,
                        int walker,
                        double time_ps,
                        double weight,
                        double rmsd_nm,
                        double rg_nm);
    int Walker() const;
    double TimePicoseconds() const;
    double BoltzmannWeight() const;
    double RmsdNanometres() const;
    double RadiusOfGyrationNm() const;
private:
    int walker_;
    double time_ps_;
    double weight_;
    double rmsd_nm_;
    double rg_nm_;
};
```

Add factory methods on Protein:

```cpp
PredictionConformation& AddPrediction(
    std::vector<Vec3> positions,
    std::string method,
    double confidence = std::nan(""));

MDFrameConformation& AddMDFrame(
    std::vector<Vec3> positions,
    int walker, double time_ps, double weight,
    double rmsd_nm, double rg_nm);

// Typed access
size_t MDFrameCount() const;
MDFrameConformation& MDFrameAt(size_t i);
size_t PredictionCount() const;
PredictionConformation& PredictionAt(size_t i);
```

Internally: `std::vector<size_t> md_frame_indices_`, same pattern
as `crystal_index_`.

No calculator sees the subtype. All calculators take
`ProteinConformation&`.

### 2. Fix OrcaRunLoader (small)

Change `AddCrystalConformation` to `AddPrediction` for AlphaFold
structures. These are computational predictions, not crystal data.

```cpp
// Before (lies about the type):
protein->AddCrystalConformation(positions, NaN, NaN, 0.0, stem);

// After (honest):
protein->AddPrediction(positions, "AlphaFold+tleap");
```

Tests that call `CrystalConf()` on OrcaRun-loaded proteins need
updating to use `PredictionAt(0)` or `ConformationAt(0)`.

### 3. Fleet loader (the main work)

New file: `src/GromacsEnsembleLoader.h/.cpp`

Input: one protein directory from `sampled_poses/`:
- `ensemble.json`: 10 poses with walker, time_ps, weight, rmsd, rg
- `pose_001.pdb` through `pose_010.pdb`: CHARMM36m-named PDBs
- TPR path for charges (from parent fleet directory)

Pipeline:
1. Read `pose_001.pdb` with ToolContext::Charmm declared
2. NamingRegistry translates residue names (HSD→HIS, HSE→HIS, etc.)
3. NamingRegistry translates atom names (HN→H, HB1→HB2, etc.)
4. Build atoms + residues with canonical names
5. Set protonation_variant_index from source residue names
   (HSD → variant 0 = delta, HSE → variant 1 = epsilon, etc.)
6. Call FinalizeConstruction (layers 2-3)
7. For each of the 10 poses:
   a. Read positions from pose PDB
   b. Verify atom count matches
   c. `protein->AddMDFrame(positions, walker, time_ps, weight, rmsd, rg)`
8. Return protein with 10 MDFrameConformations

#### CHARMM naming details

The NamingRegistry already has the residue and atom translations.
The fleet PDBs use CHARMM36m naming from pdb2gmx:

| CHARMM | Canonical | Meaning |
|--------|-----------|---------|
| HSD | HIS | delta-protonated |
| HSE | HIS | epsilon-protonated |
| HSP | HIS | doubly protonated |
| CYS2 | CYS | disulfide |
| ASPP | ASP | protonated |
| GLUP | GLU | protonated |
| HN | H | backbone amide |
| HB1 | HB2 | beta methylene |
| HB2 | HB3 | beta methylene |

The ToolContext is DECLARED (Charmm), not inferred. If an atom name
doesn't translate, that's an error.

#### Positions only from subsequent PDBs

For poses 2-10, we only need positions. The topology (atoms, residues,
rings, bonds) is from pose 1. Read the PDB, extract coordinates in
atom order. Atom count must match. If it doesn't, the PDB is corrupt.

Do NOT re-parse residues or re-detect rings for subsequent poses.
The topology is invariant.

#### Charges from TPR

```cpp
GmxTprChargeSource charges(gmx_binary, tpr_path, ForceField::CHARMM36m);
```

The TPR is the same for all poses (same simulation). Charges are
identical across frames.

### 4. CalculationRunner for fleet

Not a new Case — just a loop:

```cpp
GmxTprChargeSource charges(gmx_binary, tpr_path);
for (size_t i = 0; i < protein->MDFrameCount(); ++i) {
    auto& frame = protein->MDFrameAt(i);
    CalculationRunner::RunSingleDft(frame, charges, "");
    // Or write features immediately:
    ConformationResult::WriteAllFeatures(frame, output_dir + "/frame_" + ...);
}
```

Each frame gets its own independent calculation results. The
conformation IS the accumulator. The Boltzmann weight on the
MDFrameConformation is metadata for downstream ensemble averaging.

### 5. ProtonationDetectionResult (optional cleanup)

Currently a ConformationResult that writes Protein-level data through
const_cast. With the fleet loader setting protonation_variant_index
from CHARMM residue names at load time (step 5 above), this result
is redundant for fleet data.

For OrcaRunLoader, variant_index is already set from prmtop labels.
For LoadProtein (bare PDB), the fallback path in DetectAromaticRings
reads H atom names — this works but should eventually move to an
explicit step in FinalizeConstruction.

Elimination is desirable but not blocking. The const_cast is harmless
when there's only one conformation per protein (the current case for
ORCA data). It would be a problem with multiple conformations of the
same protein, but that doesn't happen with ORCA (each protein gets
its own Protein object).

## Testing protocol

After each change:

1. `make -j$(nproc)` — must compile
2. `./nmr_tests --gtest_filter='-Batch*'` — all non-batch tests pass
3. `./nmr_tests --gtest_filter='WriteFeatures*'` — generates NPY files
4. `cmp baseline_stash/.../*.npy baseline_features/P84477/*.npy` — byte-identical
5. Filter counts from log match baseline

For fleet loader specifically:
- Load one protein from `sampled_poses/`
- Verify atom count matches reference PDB
- Run pipeline on all 10 frames
- WriteFeatures for each frame
- Verify positions differ between frames (they're different poses)
- Verify topology is shared (same bond count, same ring count)

## Files to read

Before touching code, read:
- `spec/INDEX.md` — reading order
- `spec/PROTONATION_DESIGN_HISTORY.md` — why protonation matters
- `src/CovalentTopology.h` — the new topology type
- `src/OrcaRunLoader.cpp` — the existing loader pattern
- `src/NamingRegistry.cpp` — CHARMM translations
- `src/CalculationRunner.cpp` — the runner cases
- `GEOMETRIC_KERNEL_CATALOGUE.md` — what the calculators compute

## Hard constraints

- **No string inference.** ToolContext is declared by the caller.
- **No calculator changes.** Calculators take ProteinConformation&.
- **Protonation is not optional.** Every calculable Protein has H atoms.
- **Baselines must match.** Existing ORCA pipeline output unchanged.
- **NamingRegistry is the gatekeeper.** Unknown names are errors.
