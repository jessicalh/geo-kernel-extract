# Codex mission — UI / 3D coherence review of the h5-reader reader (READ-ONLY)

Working root: `/shared/2026Thesis/nmr-shielding/h5-reader` (branch `h5-reader-pysr-spike`). **READ-ONLY review —
no edits, no build, no git.** You own the grind; the lead + I (opus) judge. Honest stakes: the reader is the
advisor-facing vetting surface, and the lead's read is that the 3D/display has "real real issues" and that UI
features were "added whenever it felt good" — accreted ad-hoc rather than designed. She wants a code-level
coherence review (to the degree possible without seeing the running UI) that separates designed intent from
accretion and drives the surgery. Claude is reviewing the transforms + the camera-clipping bug directly in
parallel; you cover the breadth.

## The lead's specific concerns — answer EACH at the code level, with `file:line`

1. **Display transforms are incoherent.** `TransformedConformation` offers Identity / Center COM / Fit
   reference / Fit backbone. For EACH option: what does it ACTUALLY do to the displayed coordinates (the
   math)? Do the four make sense as a coherent set, or do some overlap / under-define? Is the transform
   applied CONSISTENTLY — i.e. do the calculator/snapshot/H5 values stay in the source frame while only
   *display* positions move, or do they silently diverge? Flag redundant / ill-defined options.

2. **Camera clipping cutoff.** "The box cuts off when you zoom in and the protein disappears." Trace the
   clipping-range logic: `MoleculeScene::setFrame` → `PushAtomPositions` computes bounds and sets the
   renderer clipping range (see the `feedback_vtk_bounds_cache` comment). Why would zooming IN clip the
   molecule away? Is the near/far clip computed from DATA bounds rather than camera-to-data distance, not
   re-derived on camera dolly, or fighting `CameraInputFilter`/`CameraComposer`'s zoom writes? Pin the
   mechanism and the fix (likely: derive clip from camera distance, or `ResetCameraClippingRange` on dolly).

3. **Empty fields are shown.** The UI surfaces many fields with no data — `QtFrame::eeqCharge` returns a
   placeholder `0`; many `T0` read zero (traceless-by-construction, or out-of-cutoff); some fields may be
   all-zero/absent for a given protein. Review how fields are CHOSEN for display (inspector dock, time-series
   dock, dashboard signal catalog) — is there ANY emptiness filter today? The lead's rule: **do not show
   options/fields with nothing in them, even if we must DETECT the empty ones on load.** Recommend a concrete
   load-time emptiness/sentinel-detection mechanism and exactly where it plugs in.

4. **Pane toggles.** "Why is there no button to toggle the panes?" The dock toggles are `toggleViewAction`s
   appended to the toolbar (`ReaderMainWindow.cpp` ~:417-430, kPanels: Inspector/Selection/Strip). Are they
   actually present, visible, labelled, and functional? After the docks-now-start-hidden change, is there an
   OBVIOUS way to bring a pane back? Assess and recommend (labelled toolbar group / View menu / icons).

5. **Accreted-ad-hoc cruft.** The lead: features were "added whenever it felt good by [a] predecessor."
   Identify UI features/code that read as ad-hoc, incoherent, dead, or duplicated. The overview's Cruft
   Register is a start (`QtSelectionOverlay` dormant; `/proc/self/statm` in render code; B-field
   `dynamic_cast`; render-drop PROBE; default `npy:dssp_chi` signal; "Instrument" label). Add what you find.
   Each: keep / fix / cut, with reasoning.

## Read
- `notes/UI_STATE_OVERVIEW_2026-06-04.md` + `notes/POLISH_BACKLOG.md` (prior survey + infelicity backlog).
- `src/model/TransformedConformation.{h,cpp}`; `src/app/MoleculeScene.{h,cpp}`; `src/app/CameraComposer.*`,
  `src/app/CameraInputFilter.*`, `src/app/CameraMode.h`; `src/app/ReaderMainWindow.cpp`;
  `src/app/QtAtomInspectorDock.*` + the time-series/dashboard field tables; `src/model/QtFrame.*`;
  `src/model/SignalDictionary.*` / the dashboard signal catalog.
- The `qt6-cpp` skill (`references/3d-vtk.md`, `references/architecture.md`) as the discipline lens.

## Output
Write `notes/UI_COHERENCE_REVIEW_2026-06-05.md`: per concern — what the code actually does (cited), a
coherent-vs-accreted verdict, the bug/issue, and a concrete fix. Then a top-line **coherence verdict** on the
3D/display and a **prioritised surgery list** (what to cut into first). Specific, honest, cited — no hand-waving.

## Boundaries
READ-ONLY: no edits, no git, no build. `h5-reader/src` scope. Do not touch the extractor or the running
reader (the lead may have it open). If a behavior can't be settled from code alone, say so and say what a
runtime check would resolve it.
