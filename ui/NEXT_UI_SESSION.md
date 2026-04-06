# Message in a Bottle: Next UI Session

**Written:** 2026-04-05, end of session 4.
**Context:** You are waking up 4+ sessions from now. MOPAC is integrated,
2 new calculators exist, and CalculatedRegion objects are on the
conformation. You need to display all of this. Here is everything
you need to know.

## What you MUST read before touching code

This is not optional. The previous sessions that skimped on reading
destroyed two days of working code. You need 30+ minutes of reading
before you write a single line.

**In this order:**

1. **spec/INDEX.md** — then follow its reading order
2. **spec/CONSTITUTION.md** — sign conventions, inviolable rules,
   NO ADAPTERS (you will be tempted; don't)
3. **OBJECT_MODEL.md** — the concrete types. This is what you code against.
4. **PATTERNS.md** — what gets code rejected. Read the anti-patterns.
5. **ui/CLAUDE.md** — current state, what works, what's broken
6. **This file** — you're reading it

**Then the source headers, in this order:**

7. `src/BuildResult.h` — what BuildFromPdb returns
8. `src/OperationRunner.h` — Run() and RunOptions
9. `src/Protein.h` — AtomAt, BondAt, RingAt, Conformation()
10. `src/ProteinConformation.h` — AtomAt, HasResult, Result<T>, ring_geometries
11. `src/ConformationAtom.h` — every computed field. This IS the data.
12. `src/Bond.h` — atom_index_a/b, category, order, IsPeptideBond()
13. `src/Ring.h` — type hierarchy, Intensity(), TypeName(), virtual interface

**Then the viewer code:**

14. `ui/src/ComputeWorker.h` — shared_ptr<Protein>, ComputeResult
15. `ui/src/ComputeWorker.cpp` — Build + Run + grid sampling
16. `ui/src/MainWindow.h` — protein_, pickAtom, populateAtomInfo
17. `ui/src/MainWindow.cpp` — the whole thing, especially populateAtomInfo

**Then whatever MOPAC added** (you'll need to check what's new):

18. Whatever new ConformationResult holds MOPAC bond orders
19. Whatever new fields landed on ConformationAtom or Bond
20. The CalculatedRegion type and how it attaches to conformations
21. The 2 new calculator result types

## The architecture you're inheriting

The viewer holds `shared_ptr<Protein>`. After OperationRunner::Run
completes, everything is const. The viewer reads the library objects
directly — `protein.AtomAt(i)`, `conf.AtomAt(i)`, `protein.BondAt(i)`.
There is NO adapter layer. No ViewerResults. No intermediate structs
that copy library data.

When you add display for MOPAC bond orders: you read `protein.BondAt(i)`
and whatever new field is on Bond or ConformationAtom. You do NOT
create a ViewerMopacResult struct and copy into it.

When you add display for CalculatedRegion: you read it from the
conformation directly. You do NOT create a ViewerRegionResult.

The atom inspector (double-click → QTreeWidget) already shows the
full object model for a picked atom. You ADD sections to
`populateAtomInfo()` for the new data. That's it.

## The atom inspector pattern

`populateAtomInfo(size_t atomIndex)` in MainWindow.cpp is the template.
It reads identity from `protein.AtomAt(idx)`, computed data from
`conf.AtomAt(idx)`, cross-references ring neighbours by index back
to `protein.RingAt(rn.ring_index)`, etc. Every new data type follows
the same pattern: read from the library, populate the tree.

For CalculatedRegion: you'll want a separate tree section or maybe
a separate dock panel, since regions aren't per-atom — they're
spatial regions of the conformation. You might need a region list
view where clicking a region highlights the atoms within it and
shows the region's metadata (how it was applied, what it did).

## Command line

```
build-ui/nmr-viewer [options] [path]

Options:
  -p, --pdb <file>       Load a PDB file
  --protein <dir>        Load a protein directory (WT+ALA comparison)
  --rest-port <port>     REST API port (default 9147, 0 to disable)
  -d, --dir <path>       Initial directory (unused since file menu removed)

Positional:
  path                   Auto-detect: file → --pdb, directory → --protein
```

Launch script: `bash ui/launch_viewer.sh [path]`
(starts UDP log listener, sets DISPLAY=:1)

## REST API (port 9147)

```
{"cmd":"status"}
{"cmd":"screenshot","path":"/tmp/s.png"}
{"cmd":"reset_view"}
{"cmd":"load_pdb","path":"..."}
{"cmd":"load_protein_dir","path":"..."}
{"cmd":"set_overlay","mode":"classical"}   # none/heuristic/classical/dft_delta/residual
{"cmd":"set_calculators","bs":true,"hm":false}
{"cmd":"orbit","azimuth":30}
{"cmd":"show_rings","visible":true}
{"cmd":"show_bonds","visible":true}
```

Use this loop after every change. It is seconds. Do not guess
what the viewer is doing — ask it.

## How MOPAC data will fit the display

**Bond orders:** MOPAC produces Wiberg bond order indices for every
covalent bond. These are continuous values (0.0–3.0 typically).
Bond.h already has `BondOrder order` but that's a discrete enum.
MOPAC bond orders are a float — they'll likely be a new field on
Bond or on ConformationAtom's bond_neighbours.

Display options:
- **Bond tube thickness** proportional to bond order (single thin,
  double thick, aromatic intermediate). The PeptideBondOverlay
  already builds batched tube actors — extend with per-bond radius.
- **Bond coloring** by bond order on a colormap (cool→warm).
- **In the atom inspector:** add MOPAC bond order to each bond
  neighbour's tree entry.

**MOPAC charges:** These are a third charge set (alongside ff14SB
partial charges and xTB GFN2 charges). Display as a column in the
atom inspector's Charges section. Consider a charge comparison
overlay mode.

**CalculatedRegion:** These are spatial regions where a specific
calculation was applied. They have metadata about what was computed
and what the result was. Display options:
- **Region list panel** (separate dock or tab) showing all regions
  with their descriptions
- **Click a region** → highlight constituent atoms, show metadata
- **3D overlay** — convex hull or bounding box wireframe for each
  region, color-coded by calculation type

**2 new calculators:** Same pattern as the existing 8. Each writes
SphericalTensor fields to ConformationAtom. Add to:
- `checkedCalcT0()` / `checkedCalcST()` — include in the sum
- `updateOverlay()` — they participate in heuristic tier
- `populateAtomInfo()` — new tree section
- Physics checkboxes in the sidebar — add 2 more

## Known issues

- **Crash on reload:** Loading a second protein crashes. The
  shared_ptr<Protein> cleanup races with VTK actor lifetime.
  File menu was stripped to avoid this. Fix requires careful
  ordering of VTK actor removal before protein_ reset.
- **Right-side controls:** The overlay/glyph/isosurface controls
  will likely be redesigned. They work but are from the recovery
  sessions and reflect the old architecture's thinking.
- **DFT comparison mode:** loadProteinDir sets comparisonMode but
  the comparison computation was removed with ViewerResults. Needs
  re-adding — hold a second shared_ptr for the ORCA WT protein,
  do position matching in MainWindow.
- **Types.h gap:** No NameForAtomRole() or NameForBondCategory()
  in the library. Viewer has local copies marked with GAP comment.

## What NOT to do

- Do not create adapter/wrapper/bridge classes
- Do not copy library data into flat viewer structs
- Do not modify src/ — that is the library, read-only
- Do not trust ui/warning_in_progress_ui_deletion_recovery/
- Do not write physics summaries — read the living documents
- Do not copy files between directories without asking
- Do not skip the 30 minutes of reading
