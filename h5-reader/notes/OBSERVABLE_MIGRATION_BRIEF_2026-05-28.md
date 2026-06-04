# Observable migration brief: first safe slice

> **Historical — not current truth (trued 2026-06-04).** Design history only;
> current dashboard state lives in `UI_STATE_OVERVIEW_2026-06-04.md`.

Date: 2026-05-28

This brief exists to keep the h5-reader dashboard work from drifting into a
second UI system before the measurement/observable seam is proven.

## Required context

Anyone implementing this slice must read:

- `/home/jessica/.codex/skills/qt6-cpp/SKILL.md`
- `h5-reader/CLAUDE.md`
- `h5-reader/src/model/SignalDictionary.h`
- `h5-reader/src/model/StripCalculation.h`
- `h5-reader/src/model/DftSigmaStrips.{h,cpp}`
- `h5-reader/src/app/StripChartDock.{h,cpp}`
- `h5-reader/src/app/StripStackWidget.{h,cpp}`

Qt/VTK rules from the skill still apply: keep VTK mutations in
`MoleculeScene`, use Qt services and CMake conventions, and keep Windows/macOS
as active project targets rather than treating Linux as a different app.

## Product rule

The user-facing model is:

```text
AtomSelection + category -> observable candidates -> dashboard strip
```

Even residue, ring, bond-category, electrostatic, DFT, and AIMNet2 observables
begin from an atom selection or atom-selection-derived context. The app may
look around topology and data stores to discover what can be offered, but the
selection path stays atom-centered.

## Architecture target

Long term:

```text
ObservableCatalog/Provider
  -> ObservableBinding/SignalBinding
  -> FrameSourceProvider
  -> StripCalculation
  -> dashboard panel display
```

Current code names still use `SignalDictionary` and `SignalBinding`; do not do
a broad rename in this slice. `StripCalculation` is the QObject-owned,
non-widget calculation object with a binding, a source contract, a persistent
strip buffer, and Qt signals when sampled data changes.

## Frame-source rule

Every source of per-frame calculator data should be orthogonal:

```text
frame source provider -> frame-local source data -> observers/strips sample -> source data released
```

The persistent memory is the strip buffer, not the source reader. When a frame
is visited, each active strip asks its bound source for the display value it
needs and appends that scalar/vector-reduced value to its own buffer. Once the
active observers have sampled the frame, the loaded source frame can be dropped;
the original source stays on disk.

Do not solve this with LRU source-frame caches. That was an earlier framing and
is now intentionally rejected. A source may keep tiny identity/provenance maps
or negative "known absent" records, but parsed/calculator frame payloads should
be current-frame data. If the user wants a trail, the strip calculation stores
the trail.

This applies equally to:

- per-frame NPY directories through `FrameNpyLoader` / typed group views
- ORCA DFT frames through `DftShieldingStore`
- dense H5 trajectory groups as they are migrated away from eager all-data
  ownership
- future electrostatics, AIMNet2, H-bond, SASA, water, ring, and residue
  observables

DFT must not be a privileged dashboard path. It is just another frame source.
Likewise, the NPY loader must not become a special UI path. The orthogonal unit
is "load enough source data for this frame, let strips extract their display
values, release source data."

Known current mismatch: `QtTrajectoryH5` still eagerly loads many dense H5
groups into resident buffers. That is useful for the existing application and
for the adviser-safe executable, but the observable architecture should not
copy that shape. New strips should be written as observer-owned buffers fed by
frame-source providers.

## Next Required Milestone: Reader Inventory

Do not let geometry, DFT, or the current visible dashboard controls dominate
the design. They are only the proof cases. The complexity that matters is the
full set of possible reader inputs.

Before broadening the strip picker, collect an exhaustive inventory from the
authoritative reader surfaces:

- SDK `_catalog.py`: every NPY stem the extractor can emit, grouped by
  mechanism, native axis, wrapper, shape, units, required/optional status, and
  feature flag.
- Generated C++ `QtFieldCatalog.gen.h`: confirm it matches the SDK catalog and
  flag drift.
- `FrameNpyLoader` plus all `Qt*Group` snapshot views: what per-frame NPY data
  the app can already load and decode semantically.
- `QtTrajectoryH5` plus `Qt*TimeSeries` / special buffer classes: every dense H5
  group and channel the app can already read.
- DFT/ORCA readers and any other app-local readers: treat each as another
  orthogonal frame source, not as special dashboard code.

The output should be a table of candidate source signals, not UI. The key
questions are: source provider, native axis, anchor needed from selection,
display scalar/vector choices, gap semantics, units, source-attached mask, and
whether the value is frame-local, rollup/statistical, or system-wide. That list
is what should drive `StripCalculation` subclasses and the picker, not the few
strips we happened to build first.

## This slice

Migrate the existing geometry measurement chart path so current 2/3/4 atom
measurements travel through the same strip conventions as DFT.

Scope:

- Add concrete geometry strip classes over the existing `model::Measure`
  implementation.
- Preserve current distance/angle/dihedral behavior and labels.
- Keep the existing `StripChartDock` surface looking essentially the same.
- Keep the existing shared `TimeViewportController`.
- Keep reveal behavior working from strip bindings.
- Build the experimental target.

Definition of done:

- Selecting 2 atoms still charts distance.
- Selecting 3 atoms still charts angle.
- Selecting 4 atoms still charts dihedral.
- The geometry time-domain chart is backed by a concrete strip class, not
  bespoke sampling logic embedded in `StripChartDock`.
- Existing DFT strips still work.
- Build passes:

```bash
cmake --build h5-reader/build/codex-timeviewport --target h5reader -j 4
```

## Default non-goals

Avoid adding these in the first slice unless the implementer can state a direct
reason they are needed to prove the observable seam:

- dashboard tabs
- tab naming
- tab deletion
- add-strip dialog
- panel inventory model
- persistence
- new graph button behavior
- new residue-selection behavior
- new dihedral-selection UI
- new electrostatics, AIMNet2, ring, H-bond, or SASA observables
- new assumptions about dashboard/tab layout

This is not meant to cramp useful engineering judgment. A small UI addition is
acceptable if it serves the migration directly. The rule is that layout and
panel choices should be deliberate, not accidental fallout from moving geometry.

The phrase "do not improve dihedrals" means "do not let any current special
case become an excuse to design the future UI around it." Geometry is only the
first migration proof.

## Anti-fossilization rule

Do not make the future dialog or tabs implicit in today's code. In particular,
avoid APIs like:

```text
atomSignals(atom)
residueSignals(residue)
geometrySignals(tuple)
```

Those bake today's categories into tomorrow's UI. Future discovery should be
context-based:

```text
catalog.candidates(ObservableRequestContext)
```

But this slice should only move geometry toward that shape. It should not build
the full catalog or a panel system yet.

## Agent work rules

Agents should receive a narrow file ownership scope. They are not alone in the
codebase and must not revert edits by others.

Suggested split:

- Model worker: implement geometry strip classes only.
- Integration worker/main agent: wire existing `StripChartDock` geometry path
  through those classes only after reviewing the model worker output.
- Verification: build the experimental target and manually relaunch only if the
  user wants visual inspection.
