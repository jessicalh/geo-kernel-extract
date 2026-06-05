# Codex — SUPER-THOROUGH inventory of display modes, options & mechanisms (READ-ONLY)

## FIRST — load the qt6-cpp skill (before any project code), THEN you have VTK-docs license
**Before inspecting any project code, READ these in full** (you have filesystem access — open them):
- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/3d-vtk.md`  ← VTK 9.5 + Qt: camera, pipeline, interactor —
  **CENTRAL to this inventory; do not skip it.**
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md`
- `/home/jessica/.claude/skills/qt6-cpp/references/model-view.md`
The camera / clipping / render / transform mechanisms are VTK↔Qt integration; the `3d-vtk` reference is the lens.

**You have EXPLICIT LICENSE to consult VTK documentation ONLINE** (you have network access — use it). For the
camera / renderer / clipping-range / interactor-style / `QVTKOpenGLNativeWidget` / `vtkWindowToImageFilter`
parts, look up the authoritative VTK class references (vtk.org/doc/nightly, the `vtkCamera` / `vtkRenderer` /
`vtkInteractorStyle` / `vtkRenderWindow` pages) and Qt-VTK integration guides, so your judgment of whether a
mechanism is SOUND is grounded in how VTK actually behaves — not guessed. **Cite the doc URL** for any finding
that rests on it.

This is a **portable Qt6/VTK** reader (Windows/macOS/Linux first-class). Root
`/shared/2026Thesis/nmr-shielding/h5-reader`, scope `h5-reader/src`. **READ-ONLY: inventory + judgment only.
Do NOT edit, build, launch, or git.** **Verify every claim from the code** — cite `file:line`.

## Why
We're stripping the reader to a **basic, consistent screen**: pick one of **two** stabilisation modes, do
normal zoom/pan/maneuver and **keep that view**, and **select an atom** reliably. Everything else — the
multiple camera select-lock modes, residue-only views, half-working display modes — is a candidate for
removal. Before we cut, we must **understand exactly what exists, whether it does anything, and whether it's
sound.** The lead will weigh your judgments, not rubber-stamp them — so be concrete and honest, flag
"I'm not sure" where the code is ambiguous, and separate "works" from "looks like it should work."

## The two stabilisation modes we are keeping (frame your transform/camera findings around these)
1. **Locked backbone** — the industry-standard backbone least-squares superposition: each frame's backbone is
   fit to a reference so the backbone holds visually still while sidechains/loops move.
2. **Kabsch-stabilised "with give"** — removes global tumbling/rotation but lets all real internal motion show
   (the backbone still flexes). The lead calls the misbehaving variant the "**Kabsch gyroscope**" — it may
   leave a residual spin/precession. Determine, from the code, what each existing transform mode actually does
   and which of these two (if either) it corresponds to.

## Inventory — cover EVERY item in these five groups. For each: (a) what it's meant to do, (b) does it
## actually DO anything (functional / dead / half-wired), (c) is it SOUND or broken/janky, (d) `file:line`,
## (e) keep / cut / fix for the basic screen.

1. **Display modes** (metric display): every `strip.*` and `static.*` mode id. Cross-reference
   `model::DisplayModeCapabilityFor` (`src/model/DisplayModeCapability.h`) against the ACTUAL renderers
   (`DashboardDisplayController` rebuild paths, the `*Panel` classes, `StripStackWidget`). For each mode: does
   it produce a visible surface, render nothing, or half-render? Which descriptors offer it? Note the
   `static.tensor` deferred-glyph and any mode that's offered-but-dead.
2. **Toolbar / overlays** (`ReaderMainWindow` buildToolbar/buildUi): Play/step/slider/fps, the transform
   toggle, the camera-mode group (Focus / Newman / Plane lock / Free), All-atom-fit, Instrument (hidden?),
   Metrics, Panels, and the overlays Ribbon / Rings / Butterfly / B-field. Each: what it does, enabled-when,
   functional?, sound?
3. **Stabilisation / transform** (`TransformedConformation`, `FitTargetMath`, the transform actions): the
   available modes (Identity / CenterCom / FitReference / FitSubset and the UI's backbone-vs-all-atom), the
   Kabsch math, the reference-frame choice, and **how the user switches modes + whether the current mode is
   visible anywhere** (the lead says the locked-backbone ↔ Kabsch switch is "nowhere clear"). Investigate the
   "gyroscope": does the all-atom fit leave residual rotation a backbone fit wouldn't?
4. **Camera** (`CameraComposer`, `CameraInputFilter`, the modes Free / Atom / Bond / Dihedral / Plane /
   Subset): what each mode does, which need a selection precondition, whether free zoom/pan/orbit **persists**
   (holds the spot) or gets reset/overridden (e.g. by frame ticks, by Reveal clicks, by mode switches). The
   target end-state is: Free camera that holds its view + the two stabilisation modes; the select-lock camera
   modes are likely cuttable — say which are load-bearing vs removable.
5. **Atom selection** (`QtAtomPicker`, `AtomSelection`, inspector): the pick → select → focus → inspector
   flow. Does single-atom pick work reliably (hi-DPI, Y-flip, transformed positions)? Does selection survive
   frame changes and transform switches? Any inconsistency.

## Output
Write `notes/DISPLAY_MODE_INVENTORY_2026-06-05.md`: one section per group, a table or tight list per item with
the five fields above, then a **"Basic-screen recommendation"** summary — what to KEEP (the two stabilisation
modes, Free camera that holds, atom select), what to CUT (name them), what to FIX. Mark confidence; flag
anything you couldn't determine from the code. **READ-ONLY — understanding, not changes.**
