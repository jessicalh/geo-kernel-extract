# UI Directory Rules

**DO NOT modify any files outside of ui/**

This directory contains the Qt/VTK viewer. It links against the
nmr_shielding library as a consumer. The library is not your code.
You read its headers, you call its API, you do not change it.

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

| Tool | In viewer | In CLI (nmr_extract) |
|------|-----------|----------------------|
| reduce | automatic (BuildFromPdb) | automatic (BuildFromPdb) |
| OpenBabel | automatic (BuildFromPdb) | automatic (BuildFromPdb) |
| mkdssp | runs (skip_dssp=false) | runs (skip_dssp=false) |
| MOPAC | **never** — hardcoded off in ComputeWorker | runs unless `--no-mopac` |
| Coulomb | **never** — hardcoded off in ComputeWorker | runs unless `--no-coulomb` |
| APBS | runs unless `--no-apbs` | runs unless `--no-apbs` |
| AIMNet2 | auto from `AIMNET2_MODEL` env var | auto from `AIMNET2_MODEL` env var |

MOPAC (~10 min/protein) and Coulomb (vacuum charges, 25s) are both
suppressed in the viewer. APBS (solvated Poisson-Boltzmann, 4s) is
the canonical electrostatics path. AIMNet2 is loaded once and cached
across reloads within a session.

All are installed. If a tool fails, OperationRunner logs the error
and continues — calculators that don't need it still run.

---

## Fast feedback loop

```bash
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
echo '{"cmd":"export_features","path":"/tmp/features"}' | nc -q1 localhost 9147
echo '{"cmd":"orbit","azimuth":30}' | nc -q1 localhost 9147
```

Use this loop. It is seconds, not minutes. Do not guess what the
viewer is doing — ask it.

---

## Current state (2026-04-12)

**Stable.** Four interactive JobSpec modes: PDB, ProtonatedPDB,
ORCA, Mutant. Command-line parsing shared with nmr_extract via
JobSpec (src/JobSpec.h). Trajectory mode is CLI-only (two-pass
batch over full-system XTC) — ComputeWorker returns immediately
if Trajectory mode is dispatched. Fleet mode removed 2026-04-12;
use `nmr_extract --trajectory` for GROMACS ensemble extraction.

**Calculator profile (viewer).** MOPAC and Coulomb are hardcoded
off in ComputeWorker — they are batch/calibration-path tools, not
interactive. APBS (solvated Poisson-Boltzmann) is the canonical
electrostatics path and runs by default. AIMNet2 is auto-detected
from the `AIMNET2_MODEL` environment variable, loaded once, and
cached for the session. Pass `--no-apbs` to skip APBS.

**Log tab.** MainWindow binds UDP port 9998 directly and shows the
library log stream in real time. Do NOT run `udp_listen.py`
alongside the viewer — both compete for port 9998, and Linux
unicast UDP delivers each datagram to exactly one socket. See
ui/launch_viewer.sh. `udp_listen.py` is for batch/CLI sessions.

**Feature export.** Available via REST:
`{"cmd":"export_features","path":"/tmp/out"}`.

**ViewerResults adapter removed.** ComputeWorker holds
shared_ptr<Protein>. MainWindow reads library objects directly.
Overlays take (Protein&, ProteinConformation&). No intermediate
structs.

**Atom inspector.** Double-click an atom → QTreeWidget dock
shows full object model: identity, charges, all 8 calculator
SphericalTensors, ring neighbours with G tensors and cylindrical
coords, bond neighbours with McConnell, vector fields, DSSP, ORCA
DFT. Yellow selection sphere highlights the picked atom.

**Load from CLI only.** File menu stripped. Crash on reload was
happening because shared_ptr cleanup races with VTK actor lifetime.

## Library gap: enum→string for AtomRole, BondCategory

Types.h has SymbolForElement, ThreeLetterCodeForAminoAcid,
AtomicNumberForElement. It does NOT have string conversion for
AtomRole or BondCategory enums. The viewer has local copies
(MainWindow.cpp, marked with GAP comment). These should move
to Types.h as NameForAtomRole() and NameForBondCategory() in
the same pattern as the existing enum→string functions.

## Next steps

See UI_ROADMAP.md for the forward-looking visualization plan.
Remaining known issues:

1. **Crash on reload** — command-line-only loading avoids it
2. **Tensor glyph display broken** — eigenvector rendering needs fix
3. **Sidebar reorganization** — replace overlay modes with
   per-calculator toggles (see roadmap)
