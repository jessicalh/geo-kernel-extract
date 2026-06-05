# Codex preamble — qt6-cpp skill + VTK-docs license (PREPEND TO EVERY h5-reader CODEX BRIEF)

## FIRST — load the qt6-cpp skill (before any project code), and you have VTK-docs license
**Before touching any project code, READ these in full** (you have filesystem access — open them):
- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/3d-vtk.md`  ← VTK 9.5 + Qt: camera, pipeline, interactor
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md`  ← signal/slot, ownership, startup/shutdown
- `/home/jessica/.claude/skills/qt6-cpp/references/model-view.md`  ← custom models, roles, delegates

This is a **portable Qt6/VTK** reader (Windows / macOS / Linux all first-class — treat the cross-platform source
as real). Honour the idioms: `CENSUS_REGISTER(this)` in every QObject ctor, `ACONNECT` (never raw `connect`),
`ASSERT_THREAD(this)` on thread-sensitive methods, **one** app-owned `QTimer` max, full error handling at every
external-library seam, UDP-log-first when debugging.

**You have EXPLICIT LICENSE to consult VTK documentation ONLINE** (you have network access — use it). For any
camera / renderer / clipping-range / interactor-style / `QVTKOpenGLNativeWidget` / `vtkWindowToImageFilter` /
pipeline question, look up the authoritative VTK class references (vtk.org/doc/nightly — `vtkCamera`,
`vtkRenderer`, `vtkInteractorStyle`, `vtkRenderWindow`, etc.) and Qt-VTK integration guides, so any judgment
about how VTK *actually behaves* is grounded, not guessed. **Cite the doc URL** for anything that rests on it.

**Verify from the code before changing anything, ever** — cite `file:line`. Root
`/shared/2026Thesis/nmr-shielding/h5-reader`, scope `h5-reader/src` ONLY; never read/link the `nmr_shielding`
library, never write H5, never trigger extraction.

---
(The task-specific brief follows below.)
