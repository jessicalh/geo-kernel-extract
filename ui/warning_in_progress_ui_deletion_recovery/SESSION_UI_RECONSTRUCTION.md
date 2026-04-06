# UI Reconstruction: Viewer Include Fixes Lost to Overwrites

## What Happened

On 2026-04-05, a session copied ui/src/ files from a stale directory
(nmr-shielding-clean) over the working versions in nmr-shielding.
The working versions had been adapted across two previous sessions to
build and run on batcave, wiring the viewer to work off the protein
object model (Protein, ProteinConformation, ConformationResult) instead
of raw PDB files. That work was never committed. It is now lost.

The viewer binary from the successful build still exists:
`build/nmr-viewer` (built 2026-04-05 10:40)

## What We Know

### The build dependency files survive

The .d files in `build/CMakeFiles/nmr-viewer.dir/ui/src/` record
every header that was actually included during the successful 10:40
compile. These are the ground truth for what the includes looked like.

For example, `ComputeWorker.cpp.o.d` shows that ComputeWorker.h
included `src/Types.h` (flat), NOT `Core/Types.h` (subdirectory).

### The CMake include paths

From `build/CMakeFiles/nmr-viewer.dir/flags.make`:
```
-I /shared/2026Thesis/nmr-shielding/ui/src
-I /shared/2026Thesis/nmr-shielding/src
-I /shared/2026Thesis/nmr-shielding/extern
```

So `#include "Types.h"` resolves to `src/Types.h`. The broken files
have `#include "Core/Types.h"` which does not resolve because there
is no `Core/` subdirectory.

### The pattern of fixes

The original UI code (from gotham, committed as `600ed7e`) used a
subdirectory include layout:
- `Core/Types.h` → should be `Types.h` (resolves via `-I src/`)
- `Viewer/MainWindow.h` → should be `MainWindow.h` (resolves via `-I ui/src/`)

The two previous sessions flattened these includes AND adapted the
code to work off the protein object model (ProteinConformation,
ConformationResult, ConformationAtom) instead of raw PDB data.

### What the broken files currently have

The files in ui/src/ right now are the ORIGINAL gotham versions with
`Core/` and `Viewer/` include prefixes. They contain the original
logic but not the batcave adaptations that wired them to the protein
model.

## Reconstruction Approach

1. Read every .d file in `build/CMakeFiles/nmr-viewer.dir/ui/src/`
   to determine the exact set of headers each source file included.

2. For each UI source file, compare the .d record against the current
   `#include` directives. Fix every include that doesn't match the
   .d record.

3. The include fixes alone will make it compile. But the code that
   wired the viewer to the protein model (using ProteinConformation
   instead of PDB) is harder to reconstruct. The .d files show WHICH
   headers were included (e.g., ProteinConformation.h, ConformationAtom.h,
   BiotSavartResult.h) which tells you what types the code used, but
   not the logic.

4. The existing binary at `build/nmr-viewer` can be run to see what
   the viewer actually did. The UI layout, menu items, and behaviour
   are observable even if the source is lost.

5. The git HEAD versions of the files (identical to what's on disk
   now) have the STRUCTURE of the UI code. The adaptations were
   modifications to that structure. Comparing what the code does now
   (loads PDB) against what it should do (work off Protein/
   ProteinConformation) plus knowing which result types were included
   (from .d files) constrains the reconstruction.

## Files affected

All files in `ui/src/` that have `Core/` or `Viewer/` include
prefixes. Check with:
```bash
grep -rn "Core/\|Viewer/" ui/src/
```

## DO NOT

- Do not touch any files outside ui/src/
- Do not modify the library (src/), tests, or learn/
- Do not "improve" the UI beyond restoring what was there
- Do not delete the .d files in build/ — they are the only record
