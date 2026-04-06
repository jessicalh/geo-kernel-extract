# Session Work: Replace xTB with libmopac C++ Integration

## Goal

Remove XtbChargeResult and replace it with a MopacResult that calls
libmopac (OpenMOPAC C API) in-process via `mozyme_scf()`. This is the
only code deliverable for this session. Do not write the two downstream
calculators that will consume MopacResult — but you must understand
them deeply enough that the MopacResult you build exposes everything
they will need.

This is work we would like to do together, carefully. The specification
for this project exists as living documents — OBJECT_MODEL.md,
GEOMETRIC_KERNEL_CATALOGUE.md, EXTRACTION_ORDER.md, and others. There
is no separate spec for this session. Read the living documents, observe
where they are not perfectly up to date, and note the tension between
documentation and code as you go. We will need to describe our work
at the end, and the documents will be updated to reflect what we build
if we do well.

---

## Context: Why This Matters

The NMR shielding prediction system has 8 geometric kernel calculators
that compute the spatial shape of each shielding mechanism as rank-2
tensors. Two PLANNED new calculators (not built this session) will need
per-conformation QM-derived charges and bond orders from MOPAC. xTB was
the original provider but segfaults on >450 atoms (7% success rate on
our 725-protein dataset). MOPAC PM7+MOZYME works on everything at ~45s
per protein.

Currently MOPAC runs offline as a Python subprocess writing .mop/.out
files. This session integrates it as a first-class ConformationResult
so it participates in the dependency graph, its output is available
in-memory to future downstream calculators, and the offline Python
extraction step is eliminated.

The conda-installed mopac (23.2.2) at `/home/jessica/micromamba/envs/mm/`
already provides `lib/libmopac.so.2` and `include/mopac.h`. No build
from source is needed.

---

## MANDATORY READING (do all of this before writing any code)

You must deeply understand the architecture before touching it. Read
these files IN ORDER. Do not skim. The object model, topology, and
tensor decomposition are tightly interlocked. A change that seems
local can break invariants silently.

### Architecture & Physics (read fully)

1. `OBJECT_MODEL.md` — The Protein/ProteinConformation/ConformationResult
   hierarchy. Understand: Protein owns topology (atoms, bonds, residues,
   rings) but has NO positions. ProteinConformation owns positions and
   accumulates ConformationResult objects. ConformationAtom stores
   per-atom computed properties. This separation is the entire design.

2. `GEOMETRIC_KERNEL_CATALOGUE.md` — All 8 calculators: what physical
   effect each models, what tensor rank, what parameters. Understand
   the dipolar kernel K_ab = (3 d_a d_b - delta_ab) / r^3 that four
   calculators share. Understand T0/T1/T2 decomposition.

3. `EXTRACTION_ORDER.md` — How the dependency graph drives extraction.
   ConformationResults declare Dependencies(). The pipeline resolves
   order. A new result slots in by declaring what it needs.

4. `PATTERNS.md` — Code conventions: naming, error handling, logging,
   file I/O patterns. Follow these exactly.

5. `APPLIED_MATHS_FIXES.md` — Tensor decomposition details, sign
   conventions. The SphericalTensor T0+T1+T2 decomposition is used
   everywhere.

6. `CALCULATOR_PARAMETER_API.md` — How calculator parameters are
   declared, stored, and swept.

### Calibration & MOPAC Physics (read fully)

7. `learn/CALIBRATION_CHECKLIST.md` — Every assumed constant in the
   calculators and what MOPAC data could modulate it. Sections on
   Coulomb charges (ff14SB vs MOPAC), bond anisotropy (fixed delta-chi
   vs bond-order-dependent), lobe offset (p-orbital populations),
   midpoint shift (bond order asymmetry). THIS IS THE SCIENCE CASE
   for the two future kernels. You must understand this so MopacResult
   exposes everything those kernels will need.

8. `learn/ARGUMENT.md` — The thesis argument: DFT delta as a meter,
   training learns physics not AI, residual is the boundary of
   knowledge.

9. `learn/EXTRACTOR_DESIGN.md` — How the Python-side extraction
   orchestrates C++ binary + external tools.

### Topology & Bonds (read fully — critical for MOPAC integration)

10. `src/CovalentTopology.h` and `src/CovalentTopology.cpp` — How the
    codebase defines bonds. This is GRAPH topology from the PDB/force
    field: integer bonds, categorical (single/double/aromatic/peptide).
    MOPAC's Wiberg bond orders are CONTINUOUS (1.0, 1.5, 1.8) and
    represent electron density sharing, not graph connectivity. These
    are complementary, not competing. MopacResult must NOT modify or
    replace CovalentTopology. It provides additional continuous bond
    data that downstream calculators can optionally use.

11. `src/Bond.h` — The Bond struct: atom indices, category, ring
    membership. Again, this is topological. MOPAC bond orders annotate
    these bonds, they don't replace them.

12. `src/MolecularGraphResult.h` and `.cpp` — BFS graph traversal
    using CovalentTopology. Depends on ChargeAssignmentResult. Used
    for sequential exclusion logic (min_sep filter). MopacResult
    must not interfere with this.

### The Existing Calculators That Will Consume MOPAC Data (read to understand, not to build)

13. `src/CoulombResult.h` and `src/CoulombResult.cpp` — The existing
    Coulomb EFG calculator uses ff14SB fixed charges. A future
    MopacCoulombResult will do the same physics with MOPAC QM charges.
    MopacResult must expose per-atom charges in a form that makes this
    straightforward.

14. `src/McConnellResult.h` and `src/McConnellResult.cpp` — The existing
    bond anisotropy calculator uses fixed delta-chi per bond category.
    A future MopacBondAnisotropyResult will modulate delta-chi by MOPAC
    bond order. MopacResult must expose bond orders indexed by the SAME
    atom pairs that CovalentTopology uses, so the future calculator can
    look up "what is the MOPAC bond order for this Bond?".

15. `src/ChargeAssignmentResult.h` — How ff14SB charges currently get
    onto ConformationAtom. MopacResult provides a parallel charge source,
    not a replacement.

### The Code Being Replaced

16. `src/XtbChargeResult.h` and `src/XtbChargeResult.cpp` — Read
    completely. Understand the interface: per-atom charges, bond orders,
    HOMO-LUMO gap, WriteFeatures(). MopacResult provides the same data
    categories via libmopac in-process instead of xTB subprocess. The
    struct layout (MopacAtomResult, MopacBondOrder) should parallel
    XtbAtomResult/XtbBondOrder but will be richer (add orbital
    populations, heat of formation).

17. `tests/test_xtb_charge_result.cpp` — The test pattern. New
    MopacResult needs equivalent tests.

### Existing MOPAC Python Code (understand what it produces)

18. `learn/mopac_extract.py` — The offline MOPAC extraction script.
    Understand the MOPAC keywords (PM7 MOZYME 1SCF CHARGE=N BONDS
    MULLIK LET GEO-OK THREADS=8), the output parsing, what gets
    saved to .npy files. The C++ MopacResult replaces this for
    extraction, but the keyword choices encode physics decisions that
    must be preserved in the API call parameters.

    The .npy output files are:
    - `mopac_charges.npy` — (N,) Mulliken charges
    - `mopac_scalars.npy` — (N, 3) [charge, s_pop, p_pop]
    - `mopac_bond_orders.npy` — (B, 3) [atom_i, atom_j, bond_order]
    - `mopac_global.npy` — (4,) [heat_of_formation, homo_lumo_gap, 0, 0]

    MopacResult::WriteFeatures() must produce identical files.

### The libmopac C API

19. `/home/jessica/micromamba/envs/mm/include/mopac.h` — The actual
    API header. Key structs: mopac_system (input), mopac_properties
    (output), mozyme_state (MOZYME-specific). Key function:
    `mozyme_scf(mopac_system*, mozyme_state*, mopac_properties*)`.
    The API is diskless: pass coordinates + charge, get back charges
    + bond orders + heat of formation. No .mop file needed.

    NOTE: HOMO-LUMO gap is NOT exposed in the API properties struct.
    Eigenvalues are computed internally but not returned. For gap,
    write 0.0 to mopac_global.npy and add a TODO comment. It is not
    needed by either planned kernel.

---

## What To Build (this session)

### MopacResult (ConformationResult)

Replaces XtbChargeResult. Lives in `src/MopacResult.h` and
`src/MopacResult.cpp`.

**Dependencies**: none (runs from atom positions + net charge, same
as xTB did).

**Compute()** must:
- Build a `mopac_system` struct from ProteinConformation atom positions
  and elements. Set `model = 0` (PM7). Set charge from Protein metadata.
- Call `mozyme_scf()` for linear-scaling SCF (required for proteins).
- Extract from `mopac_properties`:
  - `charge[natom]` → per-atom Mulliken charges
  - `bond_index/bond_atom/bond_order` (CSC sparse) → Wiberg bond orders
  - `heat` → heat of formation (kcal/mol)
- Store results accessibly for future downstream ConformationResults:
  - Per-atom charges queryable by atom index
  - Bond orders queryable by atom pair (so a future calculator can ask
    "what is the MOPAC bond order for the bond between atoms i and j?"
    using the same atom indices as CovalentTopology)
  - Global scalars (heat of formation)
- Clean up with `destroy_mopac_properties()` / `destroy_mozyme_state()`.
- Handle errors: check `properties.error_msg`. Log via OperationLog.

**WriteFeatures()** must write the same .npy files that mopac_extract.py
currently produces so the Python pipeline can load them without changes.

**CRITICAL**: MOPAC's atom ordering may differ from the
ProteinConformation atom ordering if elements/coords are reordered.
Verify atom correspondence carefully. The API takes atoms in the order
you provide them — use ProteinConformation order and the correspondence
is 1:1. Do NOT sort or reorder atoms.

**Query interface** — design for the future calculators:
```cpp
double ChargeAt(size_t atom_index) const;
double BondOrder(size_t atom_a, size_t atom_b) const;  // returns 0.0 if no MOPAC bond
double HeatOfFormation() const;
const vector<MopacBondOrder>& BondOrders() const;  // full list
```

The BondOrder(a, b) lookup must be efficient (the future
MopacCoulombResult will call it for every atom pair within 15 Å, and
the future MopacBondAnisotropyResult will call it for every bond in
the CovalentTopology).

### CMakeLists.txt Changes

- Remove `src/XtbChargeResult.cpp` and `tests/test_xtb_charge_result.cpp`
- Add `src/MopacResult.cpp` and `tests/test_mopac_result.cpp`
- Link libmopac:
  ```cmake
  target_include_directories(... PRIVATE /home/jessica/micromamba/envs/mm/include)
  target_link_libraries(... /home/jessica/micromamba/envs/mm/lib/libmopac.so)
  ```
- The Fortran runtime (libgfortran) and OpenBLAS are pulled in
  transitively by libmopac.so.

### Remove xTB References

- Delete `src/XtbChargeResult.h`, `src/XtbChargeResult.cpp`
- Delete or repurpose `tests/test_xtb_charge_result.cpp`
- Remove any RuntimeEnvironment::Xtb() path configuration
- Update EXTRACTION_ORDER.md dependency graph
- Grep for "xtb" / "xTB" / "Xtb" across the entire project and
  clean up references (source, docs, CMake, Python)

---

## What NOT To Do

- Do NOT build MopacCoulombResult or MopacBondAnisotropyResult.
  Those are future work requiring physics discussion. But MopacResult
  must expose everything they will need.
- Do NOT modify CovalentTopology, Bond, or MolecularGraphResult.
  MOPAC bond orders are additional data, not a replacement for graph
  topology.
- Do NOT change the Protein class. Protein is sequence-level, not
  conformation-level. MOPAC results belong on ProteinConformation.
- Do NOT change the existing CoulombResult or McConnellResult.
- Do NOT add MOPAC as a build-time requirement (it's a runtime
  shared library). If libmopac.so is absent, MopacResult::Compute()
  should return nullptr with a log message, not a compile error.
- Do NOT change the Python training pipeline (learn/*.py). The new
  C++ result writes .npy files that the existing loader reads.
- Do NOT add features, abstractions, or "improvements" beyond what
  is specified here.

---

## Verification

1. Build compiles and links against libmopac.so
2. MopacResult::Compute() produces charges and bond orders on a
   small test protein (use 1CRN or similar, <500 atoms)
3. Charges agree with mopac_extract.py output (same protein, same
   coordinates, same net charge) to 1e-4
4. Bond orders agree to 1e-4
5. WriteFeatures() produces .npy files identical to mopac_extract.py
   output (byte-level comparison after accounting for float tolerance)
6. BondOrder(a, b) returns the correct value for bonds that exist
   in CovalentTopology
7. `grep -ri xtb src/ tests/ CMakeLists.txt` returns nothing
8. All existing tests still pass (no regressions)
9. The extraction pipeline (`learn/extract.py`) works end-to-end
   with the new MopacResult replacing XtbChargeResult

---

## Why You Must Understand The Future Kernels

You are not building MopacCoulombResult or MopacBondAnisotropyResult.
But you must understand them as if you were, because they define the
requirements for MopacResult's query interface.

**MopacCoulombResult** will:
- Depend on MopacResult
- Iterate over all atoms within 15 Å of each target atom (via
  SpatialIndexResult)
- For each neighbour, get its MOPAC charge and compute the dipolar
  EFG tensor: V_ab = q_j (3 d_a d_b - delta_ab) / r^3
- This means MopacResult.ChargeAt(index) must be O(1)

**MopacBondAnisotropyResult** will:
- Depend on MopacResult
- Iterate over all bonds in CovalentTopology within 10 Å of each
  target atom (via SpatialIndexResult, same pattern as McConnellResult)
- For each bond, get its MOPAC bond order and compute delta-chi as
  a function of (category, bond_order)
- This means MopacResult.BondOrder(atom_a, atom_b) must be O(1) and
  must use the SAME atom indices as CovalentTopology::Bond

Additionally, both future kernels may use:
- Orbital populations (s_pop, p_pop) for lobe offset calibration
- Heat of formation for conformation energy weighting
- Delta-charges (WT minus ALA) for mutation effect analysis

Expose all of these. Store all MOPAC output that the API returns.
The cost of storing unused data is zero; the cost of re-running
MOPAC because you didn't store something is 45 seconds per protein.

---

## Key File Paths

```
Project root:       /shared/2026Thesis/nmr-shielding/
Source:              /shared/2026Thesis/nmr-shielding/src/
Tests:              /shared/2026Thesis/nmr-shielding/tests/
Build:              /shared/2026Thesis/nmr-shielding/build/
Python extraction:  /shared/2026Thesis/nmr-shielding/learn/
Docs:               /shared/2026Thesis/nmr-shielding/spec/
DFT data:           /shared/2026Thesis/consolidated/
Features output:    /shared/2026Thesis/nmr-shielding/learn/features/FirstExtraction/
libmopac.so:        /home/jessica/micromamba/envs/mm/lib/libmopac.so.2
mopac.h:            /home/jessica/micromamba/envs/mm/include/mopac.h
```
