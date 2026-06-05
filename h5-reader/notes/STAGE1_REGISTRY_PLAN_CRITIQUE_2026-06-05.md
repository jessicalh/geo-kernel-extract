# Stage-1 Registry Plan Critique - 2026-06-05

Read-only adversarial review of
`notes/STAGE1_REGISTRY_PLAN_2026-06-05.md` against the requested `h5-reader/src`
code. No source code was changed, no git command was run, no build was triggered,
and the producer/library outside `h5-reader` was not read.

## Verdict

FIX-PLAN-FIRST: keep the wrap-then-peel direction, but fix the migration order,
the capability wrapper shape, the availability API, and the dialog/context
contract before implementation.

Must-fixes:

1. Move or soften startup validation so step 6 is actually green while hollow
   catalog strings still exist.
2. Make `DisplayModeCapabilityFor` an out-of-line wrapper or otherwise avoid a
   header include cycle/recursion.
3. Preserve the current `strip.` prefix rule during the wrap phase.
4. Add a real storage-path/alternate-path availability API before promising
   `alternatesTried`.
5. Specify dialog runtime context for scene/Dft/trajectory-dependent
   `offerable()` checks.
6. Preserve reorient composite-panel dispatch and active-panel ref filtering in
   the panel-definition routing step.
7. Update all duplicated CMake target source lists, not only the core flat list.

## Verified Ground Truth

### Capability Call Sites

`DisplayModeCapabilityFor` is inline and string keyed
(`src/model/DisplayModeCapability.h:25`). The direct source call sites are not
exactly the six listed by the plan:

- `src/app/SignalDisplayDialog.cpp:227` -
  `DashboardModeHasVisibleSurface`, an exported dialog helper, reads
  `hasVisibleSurface`.
- `src/app/SignalDisplayDialog.cpp:1195` - candidate-offer gate reads
  `hasVisibleSurface`.
- `src/app/SignalDisplayDialog.cpp:1407` - active-signal toggle gate reads
  `hasVisibleSurface`.
- `src/model/DashboardPanelModel.cpp:86` - `IsPanelDisplayMode` reads
  `emitsPanelRef`.
- `src/model/DashboardPanelModel.cpp:101` - display-ref emission reads
  `emitsPanelRef`.
- `src/app/DashboardDisplayController.cpp:70` -
  `IsRenderableDashboardPanelMode` reads `buildsPanelWidget`.
- `src/app/DashboardDisplayController.cpp:1606` - controller panel-build gate
  reads `buildsPanelWidget`.
- `src/model/DashboardSignalModel.cpp:411` - `ModeRenderabilityFor` reads all
  three flags and copies them to JSON-facing state
  (`src/model/DashboardSignalModel.cpp:412-417`).
- REST is indirect through `ModeRenderabilityFor`, not a direct capability
  call (`src/app/RestServer.cpp:411-419`).

Wrapping can preserve behavior, but only if all of these consumers keep seeing
the same triplet for every current string, including arbitrary `strip.*`.

### The Three Flags Are Independently Consumed

Confirmed, with one correction: the three flags are not only read by three
different consumers; `DashboardSignalModel::ModeRenderabilityFor` reads all
three for inventory surfaces (`src/model/DashboardSignalModel.cpp:411-417`),
and REST exposes those values (`src/app/RestServer.cpp:411-419`).

- `hasVisibleSurface`: dialog helper/gates
  (`src/app/SignalDisplayDialog.cpp:227`,
  `src/app/SignalDisplayDialog.cpp:1195-1199`,
  `src/app/SignalDisplayDialog.cpp:1407-1410`) plus JSON renderability
  (`src/model/DashboardSignalModel.cpp:414`).
- `buildsPanelWidget`: controller helper/gate
  (`src/app/DashboardDisplayController.cpp:69-70`,
  `src/app/DashboardDisplayController.cpp:1605-1607`) plus JSON renderability
  (`src/model/DashboardSignalModel.cpp:415`).
- `emitsPanelRef`: panel-model helper/ref emission
  (`src/model/DashboardPanelModel.cpp:85-86`,
  `src/model/DashboardPanelModel.cpp:101-103`) plus JSON renderability
  (`src/model/DashboardSignalModel.cpp:416`).

### `static.tensor` Is Registered But Inert

Confirmed. The capability row is `{false, false, true}`:
`src/model/DisplayModeCapability.h:37-40`. The controller explicitly says the
tensor-glyph trigger is intentionally omitted and only the API/pointer remains:
`src/app/DashboardDisplayController.cpp:1675-1693` and
`src/app/DashboardDisplayController.cpp:1748-1751`. The catalog still advertises
`static.tensor` in helpers/descriptors, for example
`src/model/TrajectorySignalCatalog.cpp:94-100`,
`src/model/TrajectorySignalCatalog.cpp:111-117`, and
`src/model/TrajectorySignalCatalog.cpp:815-840`.

A registered-equals-offerable rule would expose this incorrectly. The registry
needs a visible-offer gate separate from definition existence.

### Ref Identity And Removal Cascade

Confirmed. `DashboardDisplayRef` identity is exactly
`signalId + displayModeId + channelId` in equality
(`src/model/DashboardPanelModel.h:21-24`) and in `stableKey()`
(`src/model/DashboardPanelModel.cpp:80-82`). Panel refs use the literal
`"panel"` channel sentinel (`src/model/DashboardPanelModel.cpp:98-103`).

The cascade is bidirectional and guarded:

- `onSignalRemoved` inserts into `signalsBeingRemoved_`, removes refs, then
  removes the guard (`src/app/DashboardSelectionController.cpp:210-218`).
- `onDisplayRefRemoved` skips guarded signals and removes a signal only when
  its refcount reaches zero (`src/app/DashboardSelectionController.cpp:220-228`).
- Mode-specific ref removal compares the same `displayModeId` string
  (`src/model/DashboardPanelModel.cpp:408-425`).

Strip history also keys on the string `displayModeId`:
`src/app/DashboardDisplayController.cpp:55-58`,
`src/app/DashboardDisplayController.cpp:1716-1743`. Stage 1 can preserve this
invariant, but only if enum-ification never changes the stored mode string or
the produced channel id.

### Empty Reality Check Support

Partially confirmed, but the plan overstates the current API.

What exists:

- Availability states include `Absent`, `NoFramePayload`, `AllMissing`,
  `AllZeroStructural`, `AllZeroObserved`, and `Available`
  (`src/model/TrajectoryFieldAvailability.h:24-31`).
- Missing records default to available/sampleable
  (`src/model/TrajectoryFieldAvailability.h:123-139`).
- `Absent`, `NoFramePayload`, `AllMissing`, and `AllZeroStructural` are not
  visible/sampleable; `AllZeroObserved` is visible data
  (`src/model/TrajectoryFieldAvailability.h:142-158`).
- Records are keyed internally by descriptor and storage path during build
  (`src/model/TrajectoryFieldAvailability.h:91-112`).
- `recordForDescriptor` is public (`src/model/TrajectoryFieldAvailability.h:118-121`).
- `byStorage_` exists but is private with no accessor
  (`src/model/TrajectoryFieldAvailability.h:402-403`).

What does not exist:

- No public `recordForStoragePath` or equivalent.
- No `alternatesTried` field or alternate-path result object. A source search
  found no `alternate`/`alternates` API in the relevant model/app code.

The UDP/structured logging path does exist: `StructuredLogger` documents stderr
plus UDP output (`src/diagnostics/StructuredLogger.h:1-14`), installs a global
Qt message handler (`src/diagnostics/StructuredLogger.cpp:72-83`), defaults to
UDP port 9997 (`src/diagnostics/StructuredLogger.cpp:98-125`), and the dashboard
category is declared/defined at `src/diagnostics/DashboardLogging.h:11-14` and
`src/diagnostics/Categories.cpp:17-20`. `DashboardSmokeSummary` already carries
series sparseness fields (`src/app/DashboardDisplayController.h:35-90`) and is
filled in `smokeSummary()` (`src/app/DashboardDisplayController.cpp:1339-1375`).

## Findings

### CRITICAL: Step 6 Cannot Be Independently Green With Debug Startup Validation

The plan puts startup validation in step 6 and removes hollow strings in step 7.
At step 6, the catalog still contains unresolved hollow modes such as
`strip.spectrum`, `static.scalar`, `static.table`, and `static.atomColor`
(`src/model/TrajectorySignalCatalog.cpp:84-100`,
`src/model/TrajectorySignalCatalog.cpp:111-132`,
`src/model/TrajectorySignalCatalog.cpp:135-145`). Unknown static modes currently
return no capability (`src/model/DisplayModeCapability.h:43-47`).

The plan says validation asserts in debug for any unresolved mode and release
logs/continues. That means step 6 is not "independently buildable + green" if a
debug run or startup-covered test exercises the validation before step 7.

Fix: either move fatal validation after hollow-string removal, make step-6
validation log-only until cleanup lands, or register explicit tombstone/cut
definitions for hollow modes so validation can distinguish "known intentionally
cut" from "unregistered drift".

### CRITICAL: The Inline Capability Wrapper Has A Likely Include Cycle

`DisplayModeCapabilityFor` is currently inline in
`src/model/DisplayModeCapability.h:25-48`. The plan says to repoint that body to
`VisualizationRegistry::instance().capabilityForMode(mode)` while the new
definition interface also returns `DisplayModeCapability`. If
`VisualizationRegistry.h` includes `VisualizationDefinition.h`, and
`VisualizationDefinition.h` includes `DisplayModeCapability.h`, then making
`DisplayModeCapability.h` include the registry to inline-call it creates a
cycle. If the registry implementation calls back into `DisplayModeCapabilityFor`
for golden behavior, it can also recurse.

Fix: turn `DisplayModeCapabilityFor` into a declaration in the header and move
the body to a new `DisplayModeCapability.cpp` that includes the registry. Add
that file to every CMake target that compiles the model/app sources. Make the
golden equivalence test compare against an explicit expected table, not the
post-flip wrapper.

### MAJOR: The Current `strip.` Prefix Rule Must Be Preserved Explicitly

Today every mode beginning with `strip.` returns visible strip capability,
without enumeration (`src/model/DisplayModeCapability.h:25-27`). The plan says
the prefix rule becomes registered definitions' `capability()` +
`legacyModeIds()` output. A finite `legacyModeIds()` list cannot reproduce the
open-ended prefix rule.

Fix: during the wrap phase, `capabilityForMode()` must keep a prefix owner for
`strip.*` or an equivalent predicate. Exact legacy string lists are fine for
static modes, but not for the strip family until the peel phase deliberately
narrows the accepted ids.

### MAJOR: `offerable()` Is Under-Specified For Hidden But Real Modes

The plan defines `offerable(ctx, descriptor)` as `supports() && isAvailable()`,
but both `static.spectrum.power` and `static.tensor` are registered/tracked while
not picker-visible:

- `static.spectrum.power` is `{false, true, true}`
  (`src/model/DisplayModeCapability.h:31-33`).
- `static.tensor` is `{false, false, true}`
  (`src/model/DisplayModeCapability.h:37-40`).

If `offerable()` ignores `capability().hasVisibleSurface`, hidden modes become
visible in the dialog. If it includes `hasVisibleSurface`, the name
`offerable()` should say so explicitly.

Fix: define separate queries, for example `supporting()`, `visibleOfferable()`,
and `trackedButHidden()`, or make `offerable()` explicitly include
`capability().hasVisibleSurface`.

### MAJOR: Availability Alternates Are Promised But Not Exposed

The plan's `ExpectedButEmpty` records include `canonicalName` and
`alternatesTried`, and it says the table keys by descriptor id and storage path.
The descriptor API exists, but storage lookup is private
(`src/model/TrajectoryFieldAvailability.h:402-403`) and there is no alternate
path API.

Fix: add an availability lookup result type before implementing the
reality-check, for example:

- canonical descriptor id
- canonical storage path
- public `recordForStoragePath(path)`
- list of alternate paths checked
- final state and finite/nonzero counts

Without that API, implementers will either log only descriptor-level state or
invent ad hoc alternate probing in the controller.

### MAJOR: Expected-But-Empty Must Walk The Unfiltered Catalog

`TrajectorySignalCatalog::setFieldAvailability()` filters `descriptors_` to only
availability-allowed descriptors (`src/model/TrajectorySignalCatalog.cpp:1176-1188`).
Unavailable descriptors remain only in `allDescriptors_`, exposed by
`allDescriptorList()` (`src/model/TrajectorySignalCatalog.h:31-35`). If
`collectExpectedButEmpty()` walks `descriptorList()`/`descriptors()`, the very
descriptors it is meant to report can already be filtered out.

Fix: specify that startup structural validation and runtime expected-empty
collection use `allDescriptorList()`, while picker-visible rows continue to use
the filtered descriptor list.

### MAJOR: Dialog Context Cannot Currently Answer The Planned Scene/Dft Gates

`VisualizationContext` includes `hasSceneOverlay`, `hasTrajectory`, and
`hasDftStore`. The dialog is intentionally renderer-agnostic
(`src/app/SignalDisplayDialog.h:1-5`) and has setters for catalog, signal model,
panel model, selection controller, protein/conformation, and selection only
(`src/app/SignalDisplayDialog.h:38-43`). `ReaderMainWindow` passes
`SceneRevealOverlay` to `DashboardStripDock`, not to the dialog
(`src/app/ReaderMainWindow.cpp:416-423`;
`src/app/DashboardStripDock.cpp:215-218`). DFT store similarly routes to the
dock/controller, not the dialog in the cited interface.

Fix: either keep Stage-1 dialog offerability independent of scene/Dft runtime
pointers, or add an explicit app-owned visualization runtime context and feed it
to both dialog and controller. For `static.tensor`, do not rely on
`hasSceneOverlay` alone; Stage 1 also lacks the user gesture the controller says
is required (`src/app/DashboardDisplayController.cpp:1675-1685`).

### MAJOR: Panel Definition Dispatch Must Preserve The Reorient Composite Path

The plan says step 5 can route panel-build dispatch through `surface()==Panel`
definitions and keep pixels identical. The current dispatch has extra behavior
outside the simple `path && mode` ladder:

- It pre-scans active reorient scalar signals for a composite sequence bar,
  requiring `static.bar.sequence` and the active panel ref
  (`src/app/DashboardDisplayController.cpp:1540-1573`).
- It suppresses individual reorient sequence bars for absorbed signals
  (`src/app/DashboardDisplayController.cpp:1634-1643`).
- It builds the composite panel after the per-signal loop
  (`src/app/DashboardDisplayController.cpp:1710-1712`).
- The active-panel filter compares the literal
  `DashboardDisplayRef{signal.id, mode, "panel"}`
  (`src/app/DashboardDisplayController.cpp:1613-1619`).

Fix: define whether composite grouping lives outside the registry as a
coordinator, or whether `SequenceBarVisualization` owns both individual and
composite reorient dispatch. Do not reduce step 5 to a per-descriptor builder
lookup.

### MAJOR: CMake/Test Target Updates Are Incomplete

The plan mentions adding new files to the main flat source list
(`CMakeLists.txt:43-230`). But test targets compile duplicated source lists:

- `h5reader_model_tests` lists model sources directly
  (`CMakeLists.txt:505-520`).
- `h5reader_app_tests` lists controller/panel/model/app sources directly
  (`CMakeLists.txt:587-647`).

If `DashboardPanelModel.cpp`, `DashboardSignalModel.cpp`, or
`DashboardDisplayController.cpp` start depending on out-of-line registry or
definition code, those test targets will fail unless updated too.

Fix: include the new model registry/definition sources, the new
`DisplayModeCapability.cpp` if added, and any shared panel-definition sources in
all affected test executables, or refactor tests to link `h5reader_core` instead
of duplicating the core source list.

### MAJOR: Singleton Lifetime And Thread Contract Are Not Precise Enough

The plan proposes a process-wide singleton holding
`std::vector<std::unique_ptr<VisualizationDefinition>>`. That can work only if
the registry is immutable after construction and registration does not depend on
global/static object constructors. Current QObject-heavy code uses
`CENSUS_REGISTER(this)` in constructors and `ASSERT_THREAD(this)` in
thread-sensitive methods, for example `DashboardDisplayController`
(`src/app/DashboardDisplayController.cpp:1191-1199`) and
`SignalDisplayDialog` (`src/app/SignalDisplayDialog.cpp:705-710`). The registry
is not a QObject, so those tools do not apply automatically.

Fix: specify a function-local-static, immutable, non-QObject registry whose
constructor builds all definitions directly; no external registration from
static initializers. All methods should be const after construction. Any query
that depends on runtime Qt objects must receive them in `VisualizationContext`
and must not store them.

### MINOR: The Plan's "Exactly Six Sites" Claim Should Be Corrected

This is not architecture-breaking, but it is an implementation tripwire. The
plan omits `DashboardModeHasVisibleSurface` at
`src/app/SignalDisplayDialog.cpp:227`, treats REST as a direct capability site
when it is indirect through `ModeRenderabilityFor`, and does not call out the
exported helper definitions at `src/model/DashboardPanelModel.cpp:85-86` and
`src/app/DashboardDisplayController.cpp:69-70` as potential API surfaces.

Fix: update the plan's call-site inventory to the direct/indirect list above.

### MINOR: `VisualizationDefinition`'s Dependency Claim Is Too Narrow

The plan says `VisualizationDefinition` depends only on `SignalDescriptor`, but
the proposed interface also mentions `TrajectoryFieldAvailability` and returns
`DisplayModeCapability`. `SignalDescriptor` lives in `DashboardSignal.h`
(`src/model/DashboardSignal.h:137-167`), while `DisplayModeCapability` lives in
`DisplayModeCapability.h:10-14` and availability lives in
`TrajectoryFieldAvailability.h`. The include/forward-declare plan matters
because this type sits between model and app.

Fix: write the header dependency plan explicitly: forward-declare
`TrajectoryFieldAvailability`, include `DashboardSignal.h` for
`SignalDescriptor`, include `DisplayModeCapability.h` only where the complete
return type is needed, and avoid app-layer includes from model headers.

## Missed Landmines

- Dialog dynamic controls: `allModeKinds()` builds fixed checkbox/filter rows
  at construction (`src/app/SignalDisplayDialog.cpp:360-377`,
  `src/app/SignalDisplayDialog.cpp:801-810`,
  `src/app/SignalDisplayDialog.cpp:843-852`,
  `src/app/SignalDisplayDialog.cpp:874-876`). Replacing this with registry
  definitions affects UI construction, proxy filtering, picker JSON, checkbox
  properties, and active-mode fallback, not just the two visible gates.
- Dialog substring matching is embedded in the proxy filter, not only the
  helper functions (`src/app/SignalDisplayDialog.cpp:620-641`).
- `onAddSelected()` still calls `catalog->canBind(binding)` before adding a
  selected mode (`src/app/SignalDisplayDialog.cpp:1268-1281`). If the dialog
  starts proposing typed definitions whose legacy mode id is no longer in
  `descriptor.temporalModes/staticModes`, `canBind()` will reject them
  (`src/model/TrajectorySignalCatalog.cpp:1265-1279`).
- Removing hollow strings changes descriptor summaries/tooltips and REST mode
  inventories, because raw mode lists are displayed directly
  (`src/app/SignalDisplayDialog.cpp:437-465`,
  `src/app/RestServer.cpp:402-420`). That is intended later, but it is a visible
  behavior change and should not be described as pure dead-string removal.
- `DashboardSignalModel::addSignal()` defaults an empty mode list to
  `AllDisplayModes(descriptor)` (`src/model/DashboardSignalModel.cpp:512-516`).
  If any caller passes an empty list after the registry rewrite, hollow or hidden
  catalog strings can still be enabled until catalog cleanup is complete.
- `static.tensor` currently can survive as a panel ref if added
  programmatically because `emitsPanelRef` is true
  (`src/model/DisplayModeCapability.h:37-40`,
  `src/model/DashboardPanelModel.cpp:101-103`). The registry must preserve that
  tracking behavior without making it user-offerable.
- The build uses C++17 (`CMakeLists.txt:15-17`). New enum headers using
  `std::uint8_t`, `std::unique_ptr`, and `std::vector` need explicit standard
  includes; do not rely on transitive Qt includes.
- The structured log exists, but adding an expected-empty vector to
  `DashboardSmokeSummary` requires REST/report serialization updates anywhere
  the summary is exposed. The plan names REST generally but should list the
  concrete summary serialization site before implementation.

## Startup Validation Boundary

The plan draws the right conceptual line: missing runtime data should be a
runtime log, not a startup failure. The proposed validation checks catalog
authoring strings against registered definitions, not file availability, so it
should not false-positive just because a descriptor's data is absent in this run.

However, the code has two traps:

- If validation walks filtered `descriptors()` after
  `setFieldAvailability()`, absent descriptors are already removed
  (`src/model/TrajectorySignalCatalog.cpp:1176-1188`), so validation and
  expected-empty logging can both miss cases.
- If validation walks `allDescriptorList()`, it will also see currently hollow
  mode strings until step 7 removes or tombstones them. With the planned
  `Q_ASSERT`, that bricks debug startup during step 6.

Fix the list choice and severity before implementation.

## One-Line Verdict

FIX-PLAN-FIRST: the architecture direction is usable, but the current plan is
not yet a safe implementation brief because steps 3, 5, and 6 have hidden
compile/order/runtime dependencies and the empty-reality-check API does not
exist as described.
