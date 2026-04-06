# Session prompt: UI recovery — crash at 10% during OperationRunner::Run

**This is a UI-only session. Do not modify anything in src/ or tests/.**

## Read first

1. Memory files (especially `project_ui_session2_state.md`)
2. `spec/SESSION_20260405_RECOVERY.md` — what was destroyed, what survives
3. `ui/CLAUDE.md` — you own ui/src/ only, library is read-only
4. `ui/LIBRARY_INTERFACE.md` — how the viewer calls the library

## Current state (commit 3f6fe62)

The UI reconstruction is done and builds clean. The viewer:
- Launches with Qt window, controls, REST server
- Loads 1ubq_protonated.pdb via BuildFromPdb (1237 bonds detected)
- Shows progress dialog, gets to ~10% (OperationRunner::Run phase)
- **CRASHES** during OperationRunner::Run

The crash is almost certainly xTB. OperationRunner::Run calls
XtbChargeResult::Compute after charges are assigned, and xTB has
known issues (segfaults >~450 atoms post-SCF on batcave). 1UBQ has
1231 atoms — well above the gate. But the xtb_max_atoms gate in
RunOptions defaults to 1000, so it SHOULD be skipped. If it's not
being skipped, or if the gate check has a bug, that's the crash.

## What to do

1. **Get UDP logging working first.** A Python UDP listener exists
   at `udp_listen.py`. Run it before the viewer. The library logs
   every calculator entry/exit to UDP port 9998. The last message
   before the crash tells you exactly which calculator died.

2. **Diagnose the crash.** Read `src/OperationRunner.cpp` to see
   the sequence. The viewer's ComputeWorker sets `opts.xtb_max_atoms`
   to the default (1000). If 1UBQ's 1231 atoms > 1000, xTB should
   be skipped with a log message. If you see "xTB skipped" in the
   UDP log, the crash is elsewhere. If you DON'T see it, the gate
   isn't working.

3. **Fix in ui/src/ComputeWorker.cpp only.** If xTB is the problem,
   the fix is to set `opts.xtb_max_atoms = 0` to disable it entirely
   for the viewer (it's expensive and not needed for visualization).
   Or set a sensible limit. Do NOT modify src/.

4. **Once it loads without crashing:** Take a REST screenshot
   (`{"cmd":"screenshot","path":"/tmp/nmr-viewer-shots/recovery.png"}`).
   The vtkMolecule is built from the object model in onComputeFinished.
   You should see the protein rendered in ball-and-stick.

5. **Then: overlay visibility.** The overlays compute but don't
   render visibly (same bug as the original binary). The field grid
   overlay creates VTK contour actors with real data but they never
   appear. Investigate with best VTK/Qt practices. This is the hard
   creative work.

6. **Visual features to add (one at a time):**
   - CSA principal axes as arrows at each atom (eigenvalues/eigenvectors
     already computed in ViewerAtomResult)
   - Conformal isosurfaces from the field grid data
   - Per-calculator tensor glyphs (SH surface from T2)
   - E-field and B-field vector arrows
   Each feature should be clean VTK, tested visually before moving on.

## Build command

```bash
cmake build-ui -DREDUCE_HET_DICT=/home/jessica/builds/reduce-src/reduce_wwPDB_het_dict.txt
cmake --build build-ui --target nmr-viewer -j$(nproc)
```

The CMakeLists.txt has gotham paths (jessicalh). The correct batcave
paths are cached in build-ui/CMakeCache.txt from the -D overrides
already applied. Don't delete the cache.

## DO NOT

- Modify src/, tests/, spec/ (except session notes)
- Copy files between directories
- Delete build/ or build-forensic/ (forensic artifacts from destroyed source)
- Build into build/ (protect the .o and .d files from the 10:40 working binary)
- Add features without understanding the physics (read GEOMETRIC_KERNEL_CATALOGUE.md)
