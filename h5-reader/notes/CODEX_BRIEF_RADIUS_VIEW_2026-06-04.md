# Codex mission — the radius / local-isolation view (the reader's headline feature)

Working root: `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`; do not switch). You
own the grind; the lead vets + owns ALL git; I (opus) review your diff before it ships. Honest stakes: this is
the advisor-facing "look here, ignore the rest" interaction — the most-wanted missing reader capability. The
reader is **live and good** right now (1P9J loads, 846 atoms, picking works); this adds local isolation
cleanly without disturbing it.

## The goal
Pick an atom + a radius → the scene shows ONLY that atom and the residues within R of it (residue-expanded);
everything else hidden. A toggle turns it on/off, a radius spinbox sets R, a clear restores the full molecule.

## Read first
- `notes/UI_STATE_OVERVIEW_2026-06-04.md` §2 "Radius / Local-Isolation View" — the exact 8-step path + file
  refs. **This is your spec; follow it.**
- `src/app/MoleculeScene.{h,cpp}` — the scene owns `molecule_` (source `vtkMolecule`), `mapper_`
  (`vtkOpenGLMoleculeMapper`, attached via `SetInputData`), `actor_`. `Build()` builds `molecule_` once;
  `setFrame(t)` → `PushAtomPositions(t,bounds)` updates atom positions in place; overlays get `setFrame`.
- `src/app/ReaderMainWindow.cpp` — toolbar + selection dock wiring (where the controls go); `AtomSelection`
  (the focus atom); `QtAtomPicker` (double-click pick).
- `src/model/QtProtein.h`, `Conformation.h`, `TransformedConformation.{h,cpp}` — residue membership + the
  `atomPosition` seam (read DISPLAYED positions via the transformed conformation).
- The `qt6-cpp` skill (`references/3d-vtk.md`, `references/architecture.md`).

## THE TRAP (the overview's warning — design around it, do not fight it)
`vtkOpenGLMoleculeMapper` has NO clean per-atom visibility — hiding atoms orphans their bonds. So DO NOT
toggle atom visibility on the source molecule. Instead:
- Add a **second `displayMolecule_`** holding ONLY the included atoms + bonds where BOTH endpoints are included.
- Maintain index maps `sourceAtom -> displayAtom (or npos)` and `displayAtom -> sourceAtom`.
- Filter active → point `mapper_` at `displayMolecule_`; filter `All` → point it back at `molecule_`.

## The build (the overview's 8 steps)
1. `AtomFilter` state on `MoleculeScene`: `All` and `WithinRadiusOf(focusAtom, radius, residueExpanded=true)`.
2. `displayMolecule_` + the two index maps.
3. On filter change: compute membership from the CURRENT transformed conformation's positions — include every
   residue with ANY atom within `radius` of the focus atom, plus the focus atom's own residue unconditionally.
4. Rebuild `displayMolecule_` (`AppendAtom` for included; `AppendBond` only when both endpoints included).
   Point the mapper at it.
5. On each `setFrame`: update only the included display-atom positions. If the focus atom moves in playback,
   recompute membership + rebuild ONLY when the included SET changes (not every frame).
6. Keep selection/measurement/reveal overlays **source-index based** — they must keep drawing highlighted
   atoms even when the display molecule is filtered (overlay lookup stays over SOURCE positions).
7. UI near the selection: an **isolate toggle**, a **radius double-spinbox** (default ~6–8 Å), a **clear**.
   Enable the toggle only when a focus atom exists. Prefer icons/tooltips over long toolbar words.
8. REST (optional): `POST /local-view { atom, radius, residue_expanded }` + `POST /local-view/clear`, mirroring
   an existing `RestServer` handler.

## The other trap — picking
The picker ray-casts over atom positions. With a filtered display molecule, picking MUST stay over **source
positions** (+ the source↔display map), NOT display-molecule atom ids — else selection breaks. The overview
says the picker already casts over source positions; keep it so, and verify selection/inspector still resolve
the right SOURCE atom while isolated.

## qt6 discipline (non-negotiable)
- `CENSUS_REGISTER(this)` in any new QObject ctor; `ACONNECT(...)` for every connection; `ASSERT_THREAD(this)`
  atop any method mutating VTK (all GUI-thread).
- NO new persistent `QTimer` (the one legitimate timer is `QtPlaybackController`).
- All renders via `MoleculeScene::requestRender(...)` — never call `Render()` from new code.
- UDP-log filter changes (membership size, radius, focus atom) at info level.

## Boundaries
- **`h5-reader/src` only.** Do NOT touch the extractor (`/shared/2026Thesis/nmr-shielding/src`), the model's
  H5/calculator loading, or the rediscover code. Scene/window/UI feature — no model rewrite, no H5 changes.
- **No git** — the lead owns ALL git. Edit + build; do not commit.
- Build to verify: `cmake --build build/linux-rwdi --target h5reader -j$(nproc)` must succeed. **Do NOT launch
  the GUI** — the reader is currently running for the lead; I relaunch it.
- Model-is-spine: read positions/residues from the typed model (transformed conformation + `QtProtein`), never
  recompute geometry elsewhere.

## Report
Files changed + why; the membership algorithm; how the index maps thread through picking/overlays/`setFrame`;
the UI controls added; build result; how to test (pick → isolate → radius → clear); any seams you're unsure of.
