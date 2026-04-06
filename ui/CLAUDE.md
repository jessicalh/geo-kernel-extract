# UI Directory Rules

**STATUS (2026-04-05): Recovery in progress.** A previous session
destroyed two days of working UI source by copying stale gotham files
over it. The code has been partially reconstructed. Some files may
still reflect the old broken implementation. Files known to be stale
are quarantined in `warning_in_progress_ui_deletion_recovery/`.

**DO NOT modify any files outside of ui/**

This directory contains the Qt/VTK viewer. It links against the
nmr_shielding library as a consumer. The library is not your code.
You read its headers, you call its API, you do not change it.

---

## CRITICAL: Architecture fix required

The current viewer uses a `ViewerResults` buffer — ComputeWorker
copies every atom, bond, ring, and tensor into plain structs, then
the GUI thread reads from those structs. This is an adapter layer
that the CONSTITUTION explicitly forbids ("No ViewerProtein adapter
wrapping Protein").

**The fix:** ComputeWorker should hold the `shared_ptr<Protein>`
(from BuildResult). After OperationRunner::Run completes, the
Protein and its ProteinConformation are fully const. The worker
signals completion, and the MainWindow reads directly from the
library's const objects: `protein.AtomAt(i)`, `conf.AtomAt(i)`,
`protein.BondAt(i)`, `conf.ring_geometries[i]`, etc.

This is a one-way process — load, compute, then read. Nothing
flows backward. Everything is const after Run. There is no
threading hazard. The buffer adds complexity and loses information
(e.g., bond atom indices were converted to positions and then
O(N²) matched back — now fixed but the pattern is wrong).

**Do this FIRST before adding features.** Building on the buffer
makes the mess worse.

---

## What you own

- `ui/src/` — all viewer source files (Qt/VTK C++)
- `ui/CLAUDE.md` — this file
- `ui/udp_listen.py` — UDP log listener
- `ui/launch_viewer.sh` — launch script with REST examples

## What you don't own

- `src/` — the library. Read-only.
- `tests/` — library tests. Read-only.
- `spec/` — library specs. Read for context, don't modify.
- `CMakeLists.txt` — you may modify the viewer target at the
  bottom, but do not modify library targets or source lists.

---

## Key library API (CURRENT — not the old Pipeline/LoadProtein)

```cpp
#include "BuildResult.h"       // BuildFromPdb, BuildFromOrca
#include "PdbFileReader.h"     // BuildFromPdb(path, pH=7.0)
#include "OperationRunner.h"   // OperationRunner::Run(conf, opts)
#include "RuntimeEnvironment.h" // RuntimeEnvironment::Load() — call first

// === Load ===
// RuntimeEnvironment::Load() must be called once at startup.
auto buildResult = BuildFromPdb(path);          // protonates with reduce, assigns ff14SB charges
// buildResult.protein  — unique_ptr<Protein>, fully constructed
// buildResult.charges  — unique_ptr<ChargeSource>
// buildResult.net_charge — integer formal charge

auto& protein = *buildResult.protein;
auto& conf = protein.Conformation();            // the primary conformation

// === Run all calculators ===
RunOptions opts;
if (buildResult.charges) {
    opts.charge_source = buildResult.charges.get();
    opts.net_charge = buildResult.net_charge;
}
auto runResult = OperationRunner::Run(conf, opts);
// runResult.attached — list of what was computed
// runResult.Ok() — success check

// === Everything is now const. Read freely. ===

// Identity (from Protein — does not change with geometry)
for (size_t i = 0; i < protein.AtomCount(); i++) {
    const auto& id = protein.AtomAt(i);     // Atom: element, pdb_atom_name, residue_index
    const auto& atom = conf.AtomAt(i);      // ConformationAtom: position + all computed fields

    Vec3 pos = atom.Position();
    AtomRole role = atom.role;
    SphericalTensor bs = atom.bs_shielding_contribution;
    SphericalTensor mc = atom.mc_shielding_contribution;
    // ... all 8 calculator contributions as SphericalTensor
    Vec3 B = atom.total_B_field;
    Vec3 E = atom.coulomb_E_total;
}

// Bonds (from Protein topology — has atom INDICES)
for (size_t i = 0; i < protein.BondCount(); i++) {
    const auto& bond = protein.BondAt(i);
    // bond.atom_index_a, bond.atom_index_b — direct indices into atom arrays
    // bond.category — BondCategory enum (PeptideCO, Aromatic, etc.)
    // Positions: conf.AtomAt(bond.atom_index_a).Position()
}

// Rings (from Protein + conformation geometry)
for (size_t i = 0; i < protein.RingCount(); i++) {
    const auto& ring = protein.RingAt(i);
    const auto& geo = conf.ring_geometries[i];
    // ring.Intensity(), ring.TypeName(), ring.VertexAtomIndices()
    // geo.center, geo.normal, geo.radius
}

// Grid sampling (for isosurfaces, streamlines)
if (conf.HasResult<BiotSavartResult>()) {
    const auto& bs = conf.Result<BiotSavartResult>();
    auto st = bs.SampleShieldingAt(point);   // SphericalTensor
    auto B  = bs.SampleBFieldAt(point);      // Vec3
}
```

---

## External tools

The library manages all external tool calls. The viewer never
calls them directly. OperationRunner handles them via RunOptions.

| Tool | When | RunOptions field |
|------|------|-----------------|
| reduce | During BuildFromPdb (protonation) | automatic |
| OpenBabel | During BuildFromPdb (bond detection) | automatic |
| mkdssp | DsspResult (secondary structure) | skip_dssp=false (default) |
| xtb | XtbChargeResult (semiempirical) | xtb_max_atoms=1000 (size gate) |
| APBS | ApbsFieldResult (solvated E-field) | automatic when charges present |

All are installed. If a tool fails, OperationRunner logs the error
and continues — calculators that don't need it still run.

---

## Fast feedback loop

```bash
# Launch with UDP logging
bash ui/launch_viewer.sh                    # default: 1ubq
bash ui/launch_viewer.sh path/to/file.pdb   # specific PDB
bash ui/launch_viewer.sh path/to/dir/       # protein directory

# REST commands (port 9147)
echo '{"cmd":"status"}' | nc -q1 localhost 9147
echo '{"cmd":"screenshot","path":"/tmp/shot.png"}' | nc -q1 localhost 9147
echo '{"cmd":"reset_view"}' | nc -q1 localhost 9147
echo '{"cmd":"load_pdb","path":"/path/to/file.pdb"}' | nc -q1 localhost 9147
echo '{"cmd":"load_protein_dir","path":"/path/to/dir"}' | nc -q1 localhost 9147
echo '{"cmd":"set_overlay","mode":"classical"}' | nc -q1 localhost 9147
echo '{"cmd":"set_calculators","bs":true,"hm":false}' | nc -q1 localhost 9147
echo '{"cmd":"orbit","azimuth":30}' | nc -q1 localhost 9147
```

Use this loop. It is seconds, not minutes. Do not guess what the
viewer is doing — ask it.

---

## Current state (2026-04-05, session 4)

**ViewerResults adapter removed.** ComputeWorker holds
shared_ptr<Protein>. MainWindow reads library objects directly.
Overlays take (Protein&, ProteinConformation&). No intermediate
structs. Net -167 lines.

**Atom inspector added.** Double-click an atom → QTreeWidget dock
shows full object model: identity, charges, all 8 calculator
SphericalTensors, ring neighbours with G tensors and cylindrical
coords, bond neighbours with McConnell, vector fields, DSSP, ORCA
DFT. Yellow selection sphere highlights the picked atom.

**File menu stripped.** Load from command line only. Crash on
reload was happening because shared_ptr cleanup races with VTK
actor lifetime. Removing the reload path eliminates it.

## Rename: "Protein Tensor Viewer"

The window title says "NMR Shielding Tensor Viewer" but the viewer
mostly shows properties of geometric tensors — McConnell dipolar
kernels, ring current G tensors, Coulomb EFG, etc. It doesn't get
as far as NMR spectra. Rename to "Protein Tensor Viewer" in
MainWindow constructor (setWindowTitle) and main_viewer.cpp
(setApplicationName).

## Library gap: enum→string for AtomRole, BondCategory

Types.h has SymbolForElement, ThreeLetterCodeForAminoAcid,
AtomicNumberForElement. It does NOT have string conversion for
AtomRole or BondCategory enums. The viewer has local copies
(MainWindow.cpp, marked with GAP comment). These should move
to Types.h as NameForAtomRole() and NameForBondCategory() in
the same pattern as the existing enum→string functions.

## Next steps (in order)

1. **Fix crash on reload** — either fix shared_ptr/VTK cleanup
   or keep command-line-only loading
2. **Test loadProteinDir** with consolidated/ ORCA data
3. **Verify overlay rendering** — toggle modes via REST, screenshot
4. **Fun graphics** — tensor glyphs, isosurfaces, butterfly fields
