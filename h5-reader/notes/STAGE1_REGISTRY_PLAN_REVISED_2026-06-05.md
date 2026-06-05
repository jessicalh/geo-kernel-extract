# Stage-1 Typed Visualisation Registry — REVISED Implementation Plan (2026-06-05)

Plan only. No code changed, no git run, no build triggered. Scope is `h5-reader/src`
only. The repo-root `nmr_shielding` library was never read or referenced. No H5 write
path is touched; no extraction is triggered. T2 is sacred throughout.

This revises [STAGE1_REGISTRY_PLAN_2026-06-05.md](STAGE1_REGISTRY_PLAN_2026-06-05.md)
to resolve every must-fix in
[STAGE1_REGISTRY_PLAN_CRITIQUE_2026-06-05.md](STAGE1_REGISTRY_PLAN_CRITIQUE_2026-06-05.md)
(verdict: FIX-PLAN-FIRST). The wrap-then-peel direction stands; the migration order,
the capability-wrapper TU layout, the availability API, the offer tri-state, the
dialog/context contract, the composite-panel dispatch, the singleton/thread contract,
and the CMake/test source lists are corrected and grounded in the actual code.

---

## CHANGES VS CRITIQUE

| # | Must-fix (critique) | Resolution in this plan | Code that makes the fix correct |
|---|---|---|---|
| 1 | **CRITICAL** — step-6 `Q_ASSERT` startup validation bricks debug startup while hollow strings still exist | **Reorder + downgrade.** Validation is introduced **WARN-only** (never `Q_ASSERT`/`qFatal`) at the *same step it lands*, walking the **unfiltered** catalog. It is hardened to a debug `Q_ASSERT` only in the **final** step, **after** hollow strings are peeled. Every step is independently green; a normal/debug load never bricks. The empty/hollow reality is always a runtime LOG. | Unknown static modes already return `{}` from `DisplayModeCapabilityFor` ([DisplayModeCapability.h:47](../src/model/DisplayModeCapability.h)); hollow strings live in catalog helpers ([TrajectorySignalCatalog.cpp:90,96,98,99,114-116](../src/model/TrajectorySignalCatalog.cpp)); `Q_ASSERT` against them mid-migration would abort a debug run. |
| 2 | **CRITICAL** — inline capability wrapper → include cycle/recursion | **`DisplayModeCapabilityFor` becomes a *declaration*; body moves to new `DisplayModeCapability.cpp`** that includes the registry. The struct + the `strip.` predicate stay header-side. `VisualizationDefinition.h` includes `DisplayModeCapability.h` for the return type; nothing in `DisplayModeCapability.h` includes the registry. The registry impl builds capabilities from an **explicit table**, never by calling back into `DisplayModeCapabilityFor`. No cycle, no recursion. | `DisplayModeCapabilityFor` is `inline` in the header today ([DisplayModeCapability.h:25-48](../src/model/DisplayModeCapability.h)); `DashboardSignal.h` (home of `SignalDescriptor`) does **not** include `DisplayModeCapability.h` ([DashboardSignal.h:1-13](../src/model/DashboardSignal.h)), so the only cycle risk is making the header call the registry. |
| 3 | **MAJOR** — the open-ended `strip.` prefix rule must be preserved | **`capabilityForMode()` keeps an explicit `strip.*` prefix arm** (first check, before the table), exactly as today. `legacyModeIds()` enumerates only the **finite static** ids; the strip family is matched by predicate, not by list, through all of Stage 1. | Every `strip.`-prefixed mode returns visible-strip caps by prefix, un-enumerated ([DisplayModeCapability.h:26-27](../src/model/DisplayModeCapability.h)). A finite list cannot reproduce this. |
| 4 | **MAJOR** — `offerable()` must not expose inert-but-real `static.tensor` / `static.spectrum.power` | **Three explicit queries** replace one fuzzy `offerable()`: `supporting()` (structural), `visibleOfferable()` (`supports && isAvailable && capability().hasVisibleSurface`), `trackedButHidden()` (`emitsPanelRef && !hasVisibleSurface`). The dialog calls `visibleOfferable()`. Inert-registered ≠ offerable. | `static.tensor` = `{false,false,true}` and `static.spectrum.power` = `{false,true,true}` ([DisplayModeCapability.h:31-40](../src/model/DisplayModeCapability.h)); both are picker-hidden yet ref-tracked. |
| 5 | **MAJOR** — the availability "alternates" API does not exist | **Drop `alternatesTried`; design on the real API.** Add a small public accessor `recordForStoragePath(path)` to `TrajectoryFieldAvailability` (the one minimal real addition; `byStorage_` already exists, just private). The reality-check uses `recordForDescriptor` then `recordForStoragePath` as the locate step. No invented alternate-probe list. | `byStorage_` is populated but private with no accessor ([TrajectoryFieldAvailability.h:402-403](../src/model/TrajectoryFieldAvailability.h)); `recordForDescriptor` is public ([:118-121](../src/model/TrajectoryFieldAvailability.h)); states are `Absent/NoFramePayload/AllMissing/AllZeroStructural/AllZeroObserved/Available` ([:24-31](../src/model/TrajectoryFieldAvailability.h)); a missing record means "assume present" ([:128-140](../src/model/TrajectoryFieldAvailability.h)). No `alternate*` symbol exists anywhere. |
| 6 | **MAJOR** — expected-but-empty must walk the UNFILTERED catalog | **Both startup validation and `collectExpectedButEmpty()` iterate `catalog->allDescriptorList()`**, never `descriptorList()`/`descriptors()`. Picker-visible rows keep using the filtered list. | `setFieldAvailability → rebuildVisibleDescriptors` filters `descriptors_` to availability-allowed only; the absent ones survive only in `allDescriptors_`, exposed by `allDescriptorList()` ([TrajectorySignalCatalog.cpp:1182-1189](../src/model/TrajectorySignalCatalog.cpp); [TrajectorySignalCatalog.h:33](../src/model/TrajectorySignalCatalog.h)). |
| 7 | **MAJOR** — dialog can't answer scene/DFT/trajectory gates | **An app-owned `VisualizationRuntime` is assembled in `ReaderMainWindow` and handed to BOTH the dialog and the controller.** Stage-1 dialog offerability uses only `hasTrajectory` (derivable) + availability; `hasSceneOverlay`/`hasDftStore` are carried for the controller's `Scene`-surface gates, not the dialog. `static.tensor` stays non-offerable via a Stage-3 gesture flag that is `false` in Stage 1 (not via `hasSceneOverlay` alone). | The dialog is "intentionally renderer-agnostic" with no scene/DFT setters ([SignalDisplayDialog.h:1-5,38-43](../src/app/SignalDisplayDialog.h)); scene/DFT are wired to the dock/controller, not the dialog ([ReaderMainWindow.cpp:245-263](../src/app/ReaderMainWindow.cpp)). The tensor trigger requires a not-yet-built user gesture ([DashboardDisplayController.cpp:1675-1693](../src/app/DashboardDisplayController.cpp)). |
| 8 | **MAJOR** — preserve the reorient composite-panel dispatch | **The composite grouping stays a controller-level coordinator, OUTSIDE the per-definition lookup.** `SequenceBarVisualization` owns *individual* reorient/iRED/dihedral sequence-bar dispatch; the **pre-scan group selection, the absorbed-signal suppression, and the post-loop composite build remain in `rebuild()`** and are explicitly out of the "route step-5 through definitions" rule. The `{signal.id, mode, "panel"}` active-panel filter is preserved verbatim. | The composite is selected by a pre-scan ([DashboardDisplayController.cpp:1535-1574](../src/app/DashboardDisplayController.cpp)), individual reorient bars are suppressed for absorbed signals ([:1640](../src/app/DashboardDisplayController.cpp)), the composite is built post-loop ([:1710-1713](../src/app/DashboardDisplayController.cpp)), and the active-panel ref is the literal `{signal.id, mode, "panel"}` ([:1613-1619](../src/app/DashboardDisplayController.cpp)). |
| 9 | **MAJOR** — CMake/test source lists incomplete | **Every new TU is added to the core lib AND to each test target that compiles the touched sources** (`h5reader_model_tests`, `h5reader_app_tests`). New `DisplayModeCapability.cpp` + registry/definition TUs are listed for both. Full touch-list in §5. | Test targets duplicate source lists: `h5reader_model_tests` lists `DashboardPanelModel.cpp`/`DashboardSignalModel.cpp`/`TrajectorySignalCatalog.cpp` ([CMakeLists.txt:505-520](../CMakeLists.txt)); `h5reader_app_tests` lists the controller + panel/model/catalog ([:587-647](../CMakeLists.txt)). Out-of-line registry code breaks them unless added. |
| 10 | **MAJOR** — singleton lifetime + thread contract imprecise | **Function-local-static, non-QObject, immutable-after-construction.** The ctor builds every definition directly in code (no static-initializer registration, no plugin list). All query methods `const`. No Qt runtime object is stored; runtime pointers arrive only inside `VisualizationContext`. Thread contract: build-once before first use on the GUI thread; reads are `const` and re-entrant. Documented; the non-QObject registry can't use `CENSUS_REGISTER`/`ASSERT_THREAD`, so the contract is stated as a comment + enforced by immutability. | Current QObject code uses `CENSUS_REGISTER`/`ASSERT_THREAD` ([DashboardDisplayController.cpp:1191-1199](../src/app/DashboardDisplayController.cpp); [SignalDisplayDialog.cpp:705-710](../src/app/SignalDisplayDialog.cpp)); the registry is intentionally not a QObject (value/strategy), so those macros don't apply — immutability is the substitute guarantee. C++17 ([CMakeLists.txt:15](../CMakeLists.txt)) requires explicit `<cstdint>/<memory>/<vector>` includes. |
| 11 | **MINOR** — "exactly six sites" + over-narrow dependency claim | **Corrected.** The capability gate is consulted at **8 direct call statements across 4 files** plus REST indirectly through `ModeRenderabilityFor`. The dependency note is rewritten: `VisualizationDefinition.h` includes `DashboardSignal.h` (for `SignalDescriptor`) **and** `DisplayModeCapability.h` (for the `capability()` return), and **forward-declares** `TrajectoryFieldAvailability` + `VisualizationContext`. | Direct sites: controller [70](../src/app/DashboardDisplayController.cpp) & [1606](../src/app/DashboardDisplayController.cpp); dialog [227](../src/app/SignalDisplayDialog.cpp), [1195](../src/app/SignalDisplayDialog.cpp), [1407](../src/app/SignalDisplayDialog.cpp); panel-model [86](../src/model/DashboardPanelModel.cpp) & [101](../src/model/DashboardPanelModel.cpp); signal-model [411](../src/model/DashboardSignalModel.cpp). REST is via `ModeRenderabilityFor` ([RestServer.cpp](../src/app/RestServer.cpp)). |

### Additional landmines from the critique, dispositioned

- **Dialog is more than two gates.** `allModeKinds()` ([SignalDisplayDialog.cpp:360-377](../src/app/SignalDisplayDialog.cpp)), the proxy `filterAcceptsRow` substring arms ([:620-642](../src/app/SignalDisplayDialog.cpp)), the candidate/active checkbox sets, `modeSummary`/tooltip, and the `modeForKind` fallback are ALL part of the step-6 dialog rewrite. The plan treats step 6 as a full dialog presentation rewrite, not a two-line gate swap (§5 step 6, R3).
- **`onAddSelected → canBind`.** The dialog must keep proposing **legacy mode ids that the descriptor actually lists**; `canBind` rejects any mode not in `temporalModes/staticModes` ([SignalDisplayDialog.cpp:1268-1281](../src/app/SignalDisplayDialog.cpp) → [TrajectorySignalCatalog.cpp:1265-1280](../src/model/TrajectorySignalCatalog.cpp)). Stage 1 keeps the legacy id as the binding key, so `canBind` still passes (§5 step 6 rule).
- **`addSignal` empty-list default.** `addSignal` defaults an empty mode list to `AllDisplayModes(descriptor)` ([DashboardSignalModel.cpp:512-516](../src/model/DashboardSignalModel.cpp)); until catalog cleanup (final step) a caller passing `{}` can still enable hollow strings. The honest-picker guarantee is therefore only complete **after** the peel step, not at step 6 — stated in §5.
- **Removing hollow strings is a visible change**, not pure dead-string removal: descriptor summaries/tooltips and REST inventories list raw mode ids ([SignalDisplayDialog.cpp:437-465](../src/app/SignalDisplayDialog.cpp); [RestServer.cpp:402-420](../src/app/RestServer.cpp)). Called out in §5 final step.
- **`DashboardSmokeSummary` serialization.** Adding an expected-empty vector requires updating wherever the summary is serialized to REST/report. The concrete site is `smokeSummary()` ([DashboardDisplayController.cpp:1339-1375](../src/app/DashboardDisplayController.cpp)) + the REST emitter; named in §4.

---

## 1. Types

New header/impl pair `src/model/VisualizationDefinition.{h,cpp}`, namespace
`h5reader::model`. **Not a `QObject`** — value/strategy objects, queried synchronously,
consistent with the existing free-function gate. The registry singleton (§2) is the one
object with lifetime.

### 1.0 Header dependency plan (resolves #2 + #11)

`VisualizationDefinition.h`:

```cpp
#pragma once
#include "DashboardSignal.h"          // SignalDescriptor, ChannelDescriptor (complete)
#include "DisplayModeCapability.h"    // DisplayModeCapability (complete; needed as return type)
#include <QString>
#include <QStringList>
#include <QVector>
#include <cstdint>                    // std::uint8_t  (C++17: no transitive guarantee)

namespace h5reader::model {
class TrajectoryFieldAvailability;    // forward-declare ONLY — no include
struct VisualizationContext;          // defined below, before the class
// ... enums + VisualizationContext + VisualizationDefinition ...
}
```

`DisplayModeCapability.h` stays a **leaf** header: it includes only `<QString>` /
`<array>` and now **declares** (not defines) `DisplayModeCapabilityFor`. It includes
nothing from the registry. This is the structural guarantee against the cycle.

### 1.1 `VisualizationType` and `DisplaySurface`

```cpp
enum class VisualizationType : std::uint8_t {
    TemporalStrip,   // generic ChannelBuffer strip (scalar/vector/tensor-derived)
    TensorGlyph,     // Stage 3 — rank-2 scene reveal (registered inert in Stage 1)
    AtomColor,       // Stage 2 — per-atom/residue colour map
    SequenceBar,     // SequenceBarPanel
    LagCurve,        // LagDecayPanel
    ChordCoupling,   // ChordCouplingPanel
    FixedFrequency,  // FixedFreqPanel
    PowerSpectrum    // PowerSpectrumPanel (registered, capability hasVisibleSurface=false)
};

enum class DisplaySurface : std::uint8_t {
    Strip,   // appended into StripStackWidget as a TemporalStripPanel-backed track
    Panel,   // an owned AbstractStripPanel pushed via setOwnedPanels()
    Scene    // a VTK overlay reveal in MoleculeScene / SceneRevealOverlay
};
```

`PowerSpectrum` is included because a real panel exists
([PowerSpectrumPanel.h:24](../src/app/PowerSpectrumPanel.h)) even though it stays
picker-hidden ([DisplayModeCapability.h:33](../src/model/DisplayModeCapability.h)).
Both enums get a `ToString(...)` overload added beside the existing `ToString` family in
`DashboardSignal.cpp` (the family begins at [DashboardSignal.cpp:24](../src/model/DashboardSignal.cpp)).

### 1.2 `VisualizationContext` and `VisualizationDefinition`

```cpp
// Runtime facts a definition needs to answer isAvailable() WITHOUT touching H5
// or the library. Assembled once by ReaderMainWindow (the app-owned
// VisualizationRuntime, §2.4) and passed by const-ref to dialog + controller.
struct VisualizationContext {
    const TrajectoryFieldAvailability* availability = nullptr; // present/absent/empty
    bool hasTrajectory = false;   // animated vs single pose
    bool hasDftStore   = false;   // ORCA frame source live (controller Scene gates)
    bool hasSceneOverlay = false; // SceneRevealOverlay wired (controller Scene gates)
    bool tensorGlyphGestureEnabled = false; // Stage-3 gesture; FALSE in Stage 1
};

class VisualizationDefinition {
public:
    virtual ~VisualizationDefinition() = default;

    virtual VisualizationType type() const = 0;
    virtual QString label() const = 0;
    virtual DisplaySurface surface() const = 0;

    // Structural: does this CLASS handle this descriptor's shape/axis at all?
    // Pure type/shape logic, no run data. supports()==false => not offered, silent.
    virtual bool supports(const SignalDescriptor& descriptor) const = 0;

    // Runtime: GIVEN supports()==true, is the backing data present in THIS run?
    // Uses ctx.availability (LOCATE-BEFORE-ABSENT, §4).
    virtual bool isAvailable(const VisualizationContext& ctx,
                             const SignalDescriptor& descriptor) const = 0;

    // Wrap-phase projection onto the legacy capability triplet consumed by the
    // dialog/panel-model/controller/signal-model. Deleted in the peel phase.
    virtual DisplayModeCapability capability() const = 0;

    // The finite STATIC mode-ids this definition owns during the wrap phase
    // (strip family is matched by prefix, NOT listed — see §2.2). Used by
    // startup validation (§2.3). Deleted in the peel phase.
    virtual QStringList legacyModeIds() const = 0;
};
```

`supports()` vs `isAvailable()` is the load-bearing split for the reality-check (§4):
`supports==false` is the structural cut (silent), `supports==true && isAvailable==false`
is expected-but-empty (logged). `capability()` and `legacyModeIds()` are **scaffolding**
that exist only while strings are still the ref key; they are deleted in the peel phase.

### 1.3 `StripVisualization` carries component policy as ENUMS (T2 sacred)

```cpp
// Replaces the strip.tensor.* / strip.vector.* string family.
enum class StripComponent : std::uint8_t {
    Auto,                                          // descriptor's natural scalar
    VectorX, VectorY, VectorZ, VectorMagnitude,    // replaces strip.vector.*
    TensorT0, TensorT1, TensorT2, TensorComponent  // replaces strip.tensor.*  (T2 SACRED)
};

class StripVisualization final : public VisualizationDefinition {
public:
    DisplaySurface surface() const override { return DisplaySurface::Strip; }
    VisualizationType type() const override { return VisualizationType::TemporalStrip; }
    // ... capability() => {true,false,false}; legacyModeIds() => {} (prefix-owned) ...
    QVector<StripComponent> componentsFor(const SignalDescriptor&) const;
};

struct StripSeriesRequest {        // typed successor to the displayModeId passed
    QString descriptorId;          // into sampleNpyValue() / sampleTensorValue() etc.
    ChannelDescriptor channel;
    StripComponent component = StripComponent::Auto;
};
```

`StripComponent` is the audit's demanded enum ("tensor/vector component policy as enums
instead of `strip.tensor.T2` text"). The four sample functions
([DashboardDisplayController.cpp:242-332](../src/app/DashboardDisplayController.cpp)) gain
a `StripComponent` overload; the existing string overload stays through Stage 1, reached
via a **1:1 `StripComponent ↔ mode-id` map** so output is byte-identical and the strip
**history key keeps using the same mode-id string** ([:55-58,1716-1743](../src/app/DashboardDisplayController.cpp)).
**`TensorT2` is first-class everywhere**; `tensor.T2Magnitude()` ([:303](../src/app/DashboardDisplayController.cpp))
is the magnitude *readout* of the preserved rank-2 object, never a collapse (T2 sacred).

### 1.4 Concrete definitions registered in Stage 1

| Definition (class) | `type()` | `surface()` | `supports()` keys on | `capability()` | Backed by (existing) |
|---|---|---|---|---|---|
| `StripVisualization` | `TemporalStrip` | `Strip` | any descriptor with a sampleable scalar/vector/tensor shape | `{true,false,false}` (prefix-owned) | `buildGenericTracks` ([:2368](../src/app/DashboardDisplayController.cpp)) |
| `SequenceBarVisualization` | `SequenceBar` | `Panel` | `storagePath ∈ {ired_order_parameters, reorientational_dynamics, dihedral_autocorrelation}` | `{true,true,true}` | `SequenceBarPanel` builders ([:1621-1655](../src/app/DashboardDisplayController.cpp)) |
| `LagCurveVisualization` | `LagCurve` | `Panel` | `storagePath ∈ {kernel_dynamics, reorientational_dynamics, dihedral_autocorrelation}` | `{true,true,true}` | `LagDecayPanel` builders ([:1630-1661](../src/app/DashboardDisplayController.cpp)) |
| `ChordCouplingVisualization` | `ChordCoupling` | `Panel` | `storagePath == kernel_coherence` | `{true,true,true}` | `ChordCouplingPanel` ([:1662](../src/app/DashboardDisplayController.cpp)) |
| `FixedFrequencyVisualization` | `FixedFrequency` | `Panel` | `storagePath == reorientational_dynamics` (J(ω)) | `{true,true,true}` | `FixedFreqPanel` ([:1666](../src/app/DashboardDisplayController.cpp)) |
| `PowerSpectrumVisualization` | `PowerSpectrum` | `Panel` | `storagePath == kernel_dynamics` | `{false,true,true}` (stays picker-hidden) | `PowerSpectrumPanel` ([:1626](../src/app/DashboardDisplayController.cpp)) |

The advanced-dynamics definitions are registered because real panels exist and the lead
has not ruled on them (audit DECISIONS OWED, item 2). Registration ≠ a keep decision;
deleting a definition removes it from the picker structurally.

**Stage 2 slot:** `AtomColorVisualization` (`type==AtomColor`, `surface==Scene`,
`supports()` = per-atom/residue scalar or `TensorT2` magnitude). Registered in Stage 2;
replaces hollow `static.atomColor`. **Do not design its colour/legend/mapper here.**

**Stage 3 slot:** `TensorGlyphVisualization` (`type==TensorGlyph`, `surface==Scene`,
`supports()` = rank-2 tensor descriptors). `capability()` = `{false,false,true}` to mirror
today's `static.tensor` row exactly (registered-but-inert; ref-tracked, never offered).
Its `isAvailable()` additionally requires `ctx.hasSceneOverlay && ctx.tensorGlyphGestureEnabled`,
the latter **false in Stage 1**. It drives the existing-but-unwired
`SceneRevealOverlay::revealTensor` ([SceneRevealOverlay.h:62](../src/app/SceneRevealOverlay.h));
the omitted trigger block ([DashboardDisplayController.cpp:1675-1693](../src/app/DashboardDisplayController.cpp))
is where its dispatch lands in Stage 3. **Do not design glyph policy here.** Stage 1 only
makes its inert status *typed and explicit*.

---

## 2. Registry + gate

### 2.1 `VisualizationRegistry`

`src/model/VisualizationRegistry.{h,cpp}`, namespace `h5reader::model`. A
**function-local-static, non-QObject, immutable-after-construction** registry holding
`std::vector<std::unique_ptr<VisualizationDefinition>>`, built once on first use.

```cpp
class VisualizationRegistry {
public:
    static const VisualizationRegistry& instance();   // function-local static; const

    // Structural: definitions whose supports()==true. (a)-cut input.
    QVector<const VisualizationDefinition*> supporting(const SignalDescriptor&) const;

    // PICKER gate: supports() && isAvailable(ctx) && capability().hasVisibleSurface.
    // This is what the dialog calls. Inert/hidden definitions never appear here.
    QVector<const VisualizationDefinition*> visibleOfferable(const VisualizationContext&,
                                                             const SignalDescriptor&) const;

    // Ref-tracked-but-not-picker-visible (emitsPanelRef && !hasVisibleSurface):
    // static.tensor + static.spectrum.power. Preserves their tracking without
    // making them user-offerable. (#4)
    QVector<const VisualizationDefinition*> trackedButHidden(const SignalDescriptor&) const;

    // Wrap-phase bridge: the capability a legacy mode-id projects. Built from an
    // EXPLICIT table + the strip.* prefix arm — never by calling back into
    // DisplayModeCapabilityFor (no recursion, #2). (#3 prefix preserved here.)
    DisplayModeCapability capabilityForMode(const QString& modeId) const;

    // Type/surface for a legacy mode-id (peel-phase successor to the dialog's
    // canonicalModeId/modeMatchesKind round-trip). nullptr for strip.* (prefix).
    const VisualizationDefinition* definitionForMode(const QString& modeId) const;

    // Startup self-check (§2.3): returns unresolved STATIC mode-ids. WARN-only
    // at the call site until the peel step. Walks allDescriptorList(). (#1,#6)
    QStringList unresolvedStaticModes(const TrajectorySignalCatalog&) const;

private:
    VisualizationRegistry();                       // builds every definition in code
    std::vector<std::unique_ptr<VisualizationDefinition>> defs_;
};
```

**Lifetime/thread contract (#10):** `instance()` returns a `const&` to a
function-local static, constructed on first use on the GUI thread (first reached during
`ReaderMainWindow` wiring, §2.4). The ctor builds all definitions directly — no static
initializers, no external registration, no plugin list. Every method is `const`; the
vector is never mutated after construction, so concurrent `const` reads are safe without
a lock. No Qt runtime object is stored; runtime pointers reach a definition only through
the `VisualizationContext` argument and are never retained. The registry is deliberately
not a `QObject`, so `CENSUS_REGISTER`/`ASSERT_THREAD` do not apply — immutability is the
substitute invariant, stated in a header comment.

### 2.2 Replacing `DisplayModeCapabilityFor` as the gate (cycle-free, #2 + #3)

`DisplayModeCapabilityFor` is consulted at these **8 direct call statements** (grep-verified,
correcting the plan's "six sites"):

1. `IsRenderableDashboardPanelMode` (`buildsPanelWidget`) — [DashboardDisplayController.cpp:70](../src/app/DashboardDisplayController.cpp)
2. controller panel-build gate (`buildsPanelWidget`) — [DashboardDisplayController.cpp:1606](../src/app/DashboardDisplayController.cpp)
3. `DashboardModeHasVisibleSurface` (`hasVisibleSurface`) — [SignalDisplayDialog.cpp:227](../src/app/SignalDisplayDialog.cpp)
4. dialog candidate-offer gate (`hasVisibleSurface`) — [SignalDisplayDialog.cpp:1195](../src/app/SignalDisplayDialog.cpp)
5. dialog active-toggle gate (`hasVisibleSurface`) — [SignalDisplayDialog.cpp:1407](../src/app/SignalDisplayDialog.cpp)
6. `IsPanelDisplayMode` (`emitsPanelRef`) — [DashboardPanelModel.cpp:86](../src/model/DashboardPanelModel.cpp)
7. display-ref emission (`emitsPanelRef`) — [DashboardPanelModel.cpp:101](../src/model/DashboardPanelModel.cpp)
8. `ModeRenderabilityFor` (all three flags → JSON) — [DashboardSignalModel.cpp:411](../src/model/DashboardSignalModel.cpp)

REST is **indirect** via `ModeRenderabilityFor` ([RestServer.cpp:402-420](../src/app/RestServer.cpp)).

**The TU move (the cycle fix):**

- `DisplayModeCapability.h` keeps `struct DisplayModeCapability` and the `detail` row
  type, and **changes** `inline DisplayModeCapability DisplayModeCapabilityFor(const QString&)`
  from an inline definition to a **plain declaration**.
- New `src/model/DisplayModeCapability.cpp` defines it as
  `return VisualizationRegistry::instance().capabilityForMode(mode);` and `#include`s
  `VisualizationRegistry.h`. Only this `.cpp` sees the registry; the header does not.
- `VisualizationRegistry::capabilityForMode` first checks the `strip.` prefix
  (returns `{true,false,false}`), then the explicit finite static table — **identical to
  today's two-stage logic** ([DisplayModeCapability.h:26-47](../src/model/DisplayModeCapability.h)).
  The strip prefix rule is preserved verbatim (#3); the static table mirrors the six rows.
- All 8 call sites keep calling `DisplayModeCapabilityFor` unchanged; only its body
  moved. The golden equivalence test compares `capabilityForMode(mode)` against an
  **explicit expected table** authored from the current header rows, **not** against the
  post-flip wrapper (so it can't tautologically pass).

### 2.3 Startup validation rule (WARN-only until the peel; #1)

New `VisualizationRegistry::unresolvedStaticModes(catalog)`: for every descriptor in
**`catalog.allDescriptorList()`** (#6), for every mode in `AllDisplayModes(descriptor)`
([DashboardSignal.cpp:285](../src/model/DashboardSignal.cpp)), skip `strip.*` (prefix-owned),
and for each remaining static mode require `definitionForMode(mode) != nullptr`. Collect
and return the unresolved set.

**Where it runs:** in `ReaderMainWindow`'s wiring path, right after the catalog is built
and availability set ([ReaderMainWindow.cpp:245-249](../src/app/ReaderMainWindow.cpp)) and
before the dialog/controller consume it. (The plan's `main_reader.cpp:204` was wrong —
that line only constructs `ReaderMainWindow`; the catalog is built inside it.)

**Severity, staged so every step is green (#1):**
- Steps where it is introduced (and through the dialog rewrite): **WARN-only.** Log each
  unresolved static mode at `qCWarning(diagnostics::cDash)` with
  `event=viz_unresolved_static_mode mode=…`. Never `Q_ASSERT`, never `qFatal`. A debug or
  release run continues; the registry already won't offer those modes.
- **Final (peel) step only:** after the hollow strings are removed from the catalog, the
  call becomes a debug `Q_ASSERT(unresolved.isEmpty())` (release stays log-and-continue).
  At that point the only way the assert can fire is genuine new drift — the intended
  signal — because all kept modes have definitions and the hollow ones are gone.

This is an authoring check (catalog strings vs registered definitions), **distinct** from
the runtime data reality-check (§4). It never false-positives on a run that simply lacks a
field's data, because it does not consult availability.

### 2.4 The app-owned `VisualizationRuntime` (resolves #7)

`ReaderMainWindow` already holds the protein, conformation, catalog, availability, DFT
store, and scene overlay. Add a small member that assembles a `VisualizationContext`
once these are wired ([ReaderMainWindow.cpp:245-263](../src/app/ReaderMainWindow.cpp)):

```cpp
model::VisualizationContext ctx;
ctx.availability      = fieldAvailability.get();
ctx.hasTrajectory     = conformation->asTrajectory() != nullptr;   // Conformation.h:80
ctx.hasDftStore       = (dftStore_ != nullptr);
ctx.hasSceneOverlay   = (sceneReveal_ != nullptr);
ctx.tensorGlyphGestureEnabled = false;                              // Stage 1: never
```

`ctx.hasTrajectory` is derivable from `Conformation::asTrajectory()`
([Conformation.h:80](../src/model/Conformation.h)); `frameCount()` is also available
([:59](../src/model/Conformation.h)). The same `ctx` is handed to the dialog (new setter
`SignalDisplayDialog::setVisualizationContext(const VisualizationContext&)`) and to the
controller. **Stage-1 dialog `visibleOfferable()` only consults `availability` +
`hasTrajectory`**; `hasSceneOverlay`/`hasDftStore`/`tensorGlyphGestureEnabled` matter only
to the controller's `Scene`-surface definitions. The dialog therefore stays
renderer-agnostic (no `SceneRevealOverlay*`/`DftShieldingStore*` ever enters it). The
context is passed by value/const-ref and not retained beyond the call, honouring #10.

---

## 3. Stage-1 scope (and only Stage 1)

**In scope:**
- The `VisualizationDefinition` hierarchy, `VisualizationType`, `DisplaySurface`,
  `StripComponent`, `VisualizationContext`, `VisualizationRegistry` (§1, §2).
- New `DisplayModeCapability.cpp` (out-of-line gate body, #2).
- The app-owned `VisualizationRuntime`/context assembly in `ReaderMainWindow` (#7).
- Migrate **strips** into `StripVisualization` with enum component policy, 1:1-mapped to
  the existing mode-id strings (history-key-stable).
- Register the existing panel definitions; route their build dispatch through the typed
  surface **except** the reorient composite coordinator, which stays in `rebuild()` (#8).
- Make hollow modes structurally un-offerable via the registry (no definition ⇒ not in
  `visibleOfferable()`): `static.table`, `static.scalar`, `static.efg`, `static.vector`,
  `static.vectorGlyph`, `static.category`, `static.per-class`, `static.rollup`,
  `static.event`, `static.system`, `static.embedding`, `static.geometry`,
  `static.topology`, and `strip.spectrum`-as-a-spectrum.
- The startup validation rule (§2.3, WARN-then-assert) and the runtime reality-check (§4).

**Out of scope (named, not designed):** Stage 2 Atom Colour; Stage 3 Tensor Glyph
(ovaloid/superquadric/ellipsoid) and its trigger gesture. Stage 1 only makes the tensor
case's inert status typed.

Stage 1 changes **no rendering**: every pixel that draws today draws the same after —
only the gate and the dispatch keys move.

---

## 4. The empty reality-check (on the REAL availability API)

The gate distinguishes two failures and treats them differently:

- **(a) No registered renderer** — `registry.supporting(descriptor)` is empty for a mode
  the catalog named. Structural cut: not offered. **Silent at runtime**; caught **once**,
  WARN-only, at startup (§2.3).
- **(b) Registered, `supports()==true`, data absent/empty at runtime** —
  `isAvailable(ctx, descriptor)==false`. The producer/data-integrity signal the reader
  exists to surface. **Logged.**

### Where it lives

`DashboardDisplayController::collectExpectedButEmpty()` — the controller is where the
catalog, live availability, and the per-descriptor decision already meet, and where
`smokeSummary()` ([DashboardDisplayController.cpp:1339-1375](../src/app/DashboardDisplayController.cpp))
already lives.

- It walks **`catalog_->allDescriptorList()`** (#6) — the UNFILTERED catalog — so the very
  descriptors that `rebuildVisibleDescriptors()` filtered out
  ([TrajectorySignalCatalog.cpp:1182-1189](../src/model/TrajectorySignalCatalog.cpp)) are
  exactly the ones it can still report.
- For each registered definition with `supports(descriptor)==true`, evaluate
  `isAvailable(ctx, descriptor)`. Each `supports && !isAvailable` becomes one
  `ExpectedButEmpty{ descriptorId, storagePath, visualizationType, canonicalState,
  storagePathState }` record (no `alternatesTried`, #5).
- Records are emitted to the **UDP/structured log** under
  `event=viz_expected_but_empty …` via `qCWarning(diagnostics::cDash)` (the project's
  primary debug channel; the category is declared at
  [DashboardLogging.h](../src/diagnostics/DashboardLogging.h)), and folded into
  `DashboardSmokeSummary` as a new `QVector<ExpectedButEmpty>` so the REST inventory picks
  them up. **The summary serialization site that must be updated is `smokeSummary()` +
  its REST emitter** — named here per the critique.
- It runs once per context change (catalog/availability set), **not** per frame — same
  cadence as `rebuild()`'s structural work, never on a `setFrame()` tick (one-QTimer
  discipline: `frameChanged` subscribers don't re-scan integrity).

### LOCATE-BEFORE-ABSENT on the real API (#5)

`isAvailable()` must not declare "empty" from a single canonical-name miss. The real API:

1. Look up the descriptor's record first: `availability->recordForDescriptor(descriptor.id)`
   (public, [TrajectoryFieldAvailability.h:118-121](../src/model/TrajectoryFieldAvailability.h)).
2. **No record ⇒ "assume present"** — `allowsDescriptor`/`canSampleDescriptor` return
   `true` on a null record ([:128-140](../src/model/TrajectoryFieldAvailability.h)). A bare
   miss is **not** an empty signal; do not log it.
3. If there is a record, also locate via the storage path:
   `availability->recordForStoragePath(descriptor.storagePath)` — **the one minimal real
   addition.** `byStorage_` is already populated during `Build`
   ([:110-111](../src/model/TrajectoryFieldAvailability.h)); it is currently private with
   no accessor ([:402-403](../src/model/TrajectoryFieldAvailability.h)). Add:

   ```cpp
   const TrajectoryFieldAvailabilityRecord* recordForStoragePath(const QString& path) const {
       const auto it = byStorage_.constFind(path);
       return it == byStorage_.constEnd() ? nullptr : &it.value();
   }
   ```

   This is additive, header-only, touches no producer, writes no H5. It is the "locate"
   half so a per-component descriptor can be cross-checked against its packed storage path.
4. Only a record whose `state ∈ {Absent, NoFramePayload, AllMissing}` — a *positive*
   "scanned and found empty" — qualifies as expected-but-empty. `AllZeroObserved` is data,
   not absence; `AllZeroStructural` is a separate structural-zero note, logged distinctly.
   These map directly to the existing states + `isVisibleState`
   ([:142-158](../src/model/TrajectoryFieldAvailability.h)).
5. The `ExpectedButEmpty` record carries both the descriptor-id state and the
   storage-path state, so the log line proves the locate step happened. "Absent under the
   canonical name" is never reported as "not in the data" without the storage-path lookup.

This keeps the check honest both ways: no crying wolf on descriptors a run simply omits,
no swallowing a genuinely-empty field the producer was expected to fill.

---

## 5. Migration approach — incremental, each step independently buildable + green

Invariant through Stage 1: **the `displayModeId` string stays the panel-ref key and the
on-the-wire identity** (removal cascade + persistence untouched, §6). The registry rides
alongside the strings as the authority; the strings peel later.

### File-by-file touch list (CMake-complete, #9)

**New files** (add to `h5reader_core` flat list at
[CMakeLists.txt:102-122 region](../CMakeLists.txt), beside the other `model/` entries):

- `src/model/VisualizationDefinition.h` / `.cpp` — types (§1).
- `src/model/VisualizationRegistry.h` / `.cpp` — registry + gate (§2).
- `src/model/DisplayModeCapability.cpp` — out-of-line gate body (#2). *(The header
  already exists and is already implicitly compiled via includers; the new `.cpp` must be
  added to every target that needs the symbol — see below.)*
- `src/model/StripVisualization.{h,cpp}`, `SequenceBarVisualization.{h,cpp}`,
  `LagCurveVisualization.{h,cpp}`, `ChordCouplingVisualization.{h,cpp}`,
  `FixedFrequencyVisualization.{h,cpp}`, `PowerSpectrumVisualization.{h,cpp}` — concrete
  definitions (may collapse to one `.cpp` if the lead prefers fewer files).

**Test targets** (each duplicates a source list — must be updated, #9):

- `h5reader_model_tests` ([CMakeLists.txt:505-520](../CMakeLists.txt)) compiles
  `DashboardPanelModel.cpp` / `DashboardSignalModel.cpp` / `TrajectorySignalCatalog.cpp`,
  all of which will reference the out-of-line `DisplayModeCapabilityFor`. **Add:**
  `DisplayModeCapability.cpp`, `VisualizationRegistry.{h,cpp}`,
  `VisualizationDefinition.{h,cpp}`, and every concrete-definition `.{h,cpp}`.
- `h5reader_app_tests` ([CMakeLists.txt:587-647](../CMakeLists.txt)) compiles the
  controller + panel/signal/catalog sources. **Add the same** new TUs.
- *(Simpler alternative the lead may prefer: have `h5reader_model_tests` /
  `h5reader_app_tests` link `h5reader_core` instead of re-listing core sources. That is a
  larger CMake change; the conservative per-target addition above is the default.)*

**Modified files:**

- `src/model/DisplayModeCapability.h` — `DisplayModeCapabilityFor` becomes a declaration;
  struct + `detail` row stay (#2).
- `src/model/TrajectoryFieldAvailability.h` — add `recordForStoragePath(path)` accessor
  (#5). Header-only, additive.
- `src/model/DashboardSignal.cpp` — add `ToString(VisualizationType)` /
  `ToString(DisplaySurface)` beside the existing family ([:24](../src/model/DashboardSignal.cpp)).
- `src/app/DashboardDisplayController.{h,cpp}` — (i) panel-build dispatch
  ([:1605-1673](../src/app/DashboardDisplayController.cpp)) asks the registry which
  `surface()==Panel` definition matches **except** the reorient composite coordinator,
  which stays as-is (#8); (ii) the four sample functions
  ([:242-332](../src/app/DashboardDisplayController.cpp)) gain `StripComponent` overloads;
  (iii) `channelsForMode`/`modeWantsChannel`/`canonicalModeChannel`
  ([:75-90](../src/app/DashboardDisplayController.cpp)) re-expressed via
  `StripVisualization::componentsFor`; (iv) add `collectExpectedButEmpty()` + the
  `ExpectedButEmpty` summary field (§4); (v) accept the `VisualizationContext` (#7).
- `src/model/DashboardPanelModel.cpp` — `canonicalModeChannel`/`modeWantsChannel`
  ([:21-36](../src/model/DashboardPanelModel.cpp)) and the `emitsPanelRef` branch
  ([:101](../src/model/DashboardPanelModel.cpp)) re-expressed via the registry, keeping the
  `{signalId, mode, channelId}` / `"panel"` ref shape **exactly**
  ([:102](../src/model/DashboardPanelModel.cpp)).
- `src/app/SignalDisplayDialog.{h,cpp}` — the big one: replace the anonymous
  `DisplayModeKind` enum and `canonicalModeId`/`modeMatchesKind`/`allModeKinds`/`modeForKind`/
  `modeSummary` and the **proxy `filterAcceptsRow` substring arms**
  ([:56-222,360-377,620-642](../src/app/SignalDisplayDialog.cpp)) with `VisualizationType`
  + registry queries; the offer gates ([:1195-1208,1399-1416](../src/app/SignalDisplayDialog.cpp))
  call `registry.visibleOfferable(ctx, descriptor)`; add
  `setVisualizationContext(...)` (#7). The dialog still proposes the descriptor's
  **legacy mode id** so `onAddSelected → canBind` keeps passing
  ([:1268-1281](../src/app/SignalDisplayDialog.cpp) → [TrajectorySignalCatalog.cpp:1265-1280](../src/model/TrajectorySignalCatalog.cpp)).
- `src/model/TrajectorySignalCatalog.cpp` — remove the hollow mode strings from the helper
  lists ([:90,96,98,99,114-116, …](../src/model/TrajectorySignalCatalog.cpp)) — done
  **last**, after the dialog no longer reads them.
- `src/app/ReaderMainWindow.cpp` — assemble the `VisualizationContext` and call
  `registry.unresolvedStaticModes(catalog)` in the wiring path
  ([:245-263](../src/app/ReaderMainWindow.cpp)); pass the context to dialog + controller.
- `src/model/DashboardSignalModel.cpp` ([:411](../src/model/DashboardSignalModel.cpp)) and
  `src/app/RestServer.cpp` — JSON renderability now sourced from the registry-backed
  `DisplayModeCapabilityFor` (mechanical; `ModeRenderabilityFor` body unchanged because
  the function it calls moved out-of-line). REST summary emitter extended for the new
  `ExpectedButEmpty` vector (§4).

### Ordered build steps (each independently buildable + green)

1. **Add types, no wiring.** Land `VisualizationDefinition.{h,cpp}`, the enums,
   `VisualizationContext`, `StripComponent`, and the `ToString` overloads. Add
   `recordForStoragePath` to availability. Nothing consumes them. Builds; behaviour
   identical.
2. **Add the registry + concrete definitions, mirroring today's table.** Each
   `capability()` + `legacyModeIds()` reproduces the exact static rows; the strip family
   is the prefix arm in `capabilityForMode`. Add the golden test:
   `capabilityForMode(mode) == <explicit expected table>` for the full catalog mode-id set
   (compared against a hand-authored table, not the wrapper — #2). Update the two test
   targets' source lists (#9). Still nothing wired; builds green.
3. **Move `DisplayModeCapabilityFor` out-of-line into `DisplayModeCapability.cpp`,
   delegating to the registry** (#2). All 8 consumers now read registry-backed answers
   through the unchanged function name. The step-2 golden test guarantees byte-identical
   behaviour. Add `DisplayModeCapability.cpp` to `h5reader_core` and both test targets
   (#9). Builds; UI unchanged.
4. **Enum-ify strip component policy** in controller + panel-model, behind the 1:1
   `StripComponent ↔ mode-id` map. String sample overloads still exist and are reached via
   the map; history key still uses the mode-id string. A regression test scrubs, rebuilds,
   asserts strip history length preserved (R4). Builds; strips render identically.
5. **Route `surface()==Panel` per-signal dispatch through definitions** instead of the
   `path && mode` ladder ([:1621-1671](../src/app/DashboardDisplayController.cpp)) — **but
   leave the reorient composite coordinator intact** (#8): the pre-scan group selection
   ([:1535-1574](../src/app/DashboardDisplayController.cpp)), absorbed-signal suppression
   ([:1640](../src/app/DashboardDisplayController.cpp)), and post-loop composite build
   ([:1710-1713](../src/app/DashboardDisplayController.cpp)) stay in `rebuild()`. The
   `{signal.id, mode, "panel"}` active-panel filter is unchanged. Same builders, selected
   by definition. Builds; panels render identically.
6. **Rewrite the dialog onto `VisualizationType` + `registry.visibleOfferable()`.** This is
   the full dialog presentation rewrite (enum, `modeForKind`, `allModeKinds`, the proxy
   filter substring arms, `modeSummary`, both gates). Wire `setVisualizationContext`.
   Introduce the startup validation call **WARN-only** (#1). The runtime reality-check
   (§4) lands here. Hollow modes stop being offered (no definition ⇒ not in
   `visibleOfferable()`). The dialog still proposes legacy ids (canBind-safe). Builds; the
   picker shows only backed visualisations. *(Note: because `addSignal` defaults an empty
   mode list to `AllDisplayModes`, the honest-picker guarantee is not fully closed until
   step 7 removes the hollow strings.)*
7. **Peel the hollow strings from the catalog helpers**
   ([:90,96,98,99,114-116, …](../src/model/TrajectorySignalCatalog.cpp)) and **harden
   startup validation to a debug `Q_ASSERT`** (release log-and-continue) now that no kept
   mode is unresolved and the hollow ones are gone (#1). This also changes descriptor
   summaries/tooltips + REST mode inventories (a visible change, not pure dead-string
   removal). Re-run validation: clean. Builds green.
8. **(Stage 1.5, optional within Stage 1)** begin peeling `legacyModeIds()`/`capability()`
   scaffolding once no consumer needs the string round-trip. Boundary into the typed-ref
   future (Stage 2+). Stop here for Stage 1 if the lead prefers.

Steps 1–3 are pure infrastructure (no UI change). Steps 4–5 are behaviour-preserving
dispatch swaps. Step 6 is the visible change (honest picker). Step 7 is cleanup + the
assert hardening. Each compiles and passes the existing suite on its own, and no step
(including a debug startup) bricks (#1).

---

## 6. Risks and how the plan de-risks each

**R1 — `DashboardSelectionController` removal cascade.** Bidirectional, re-entrancy-guarded:
`onSignalRemoved → removeDisplayRefsForSignal` (guarded by `signalsBeingRemoved_`),
`onDisplayRefRemoved → removeSignal` at refcount zero
([DashboardSelectionController.cpp:210-228](../src/app/DashboardSelectionController.cpp)).
*De-risk:* Stage 1 does not touch ref identity. Refs stay `{signalId, displayModeId,
channelId}` strings ([DashboardPanelModel.h:21-24](../src/model/DashboardPanelModel.h);
`stableKey()` [DashboardPanelModel.cpp:80-82](../src/model/DashboardPanelModel.cpp)); the
registry only changes which capability flag a mode reports. The cascade is in the no-touch
set.

**R2 — The panel-ref path.** The active-panel filter compares a literal
`{signal.id, mode, "panel"}` ([DashboardDisplayController.cpp:1613-1619](../src/app/DashboardDisplayController.cpp)),
and `emitsPanelRef` decides emission ([DashboardPanelModel.cpp:101](../src/model/DashboardPanelModel.cpp)).
*De-risk:* `capability().emitsPanelRef` is preserved per definition (panel defs true;
`TensorGlyphVisualization` true-but-inert; `PowerSpectrumVisualization` true). Step 5 swaps
builder *selection* but keeps the `"panel"` sentinel, so the comparison is unchanged. The
step-2 golden test pins this. `trackedButHidden()` (#4) preserves `static.tensor` /
`static.spectrum.power` tracking without offering them.

**R3 — The dialog mode-kind heuristics.** Substring matching
(`lower.contains("spectrum")`, [SignalDisplayDialog.cpp:188-196](../src/app/SignalDisplayDialog.cpp))
lives in BOTH the helpers and the proxy `filterAcceptsRow`
([:620-642](../src/app/SignalDisplayDialog.cpp)); `allModeKinds()` builds fixed rows at
construction ([:360-377](../src/app/SignalDisplayDialog.cpp)). *De-risk:* sequence the full
dialog rewrite **last** (step 6), after the registry is the gate (step 3) and dispatch is
typed (steps 4–5), so it is a presentation swap over a verified-stable backend.
`visibleOfferable()` replaces the `hasVisibleSurface && descriptorAvailable` compound
([:1199](../src/app/SignalDisplayDialog.cpp)) with one typed call; the `modeForKind`
fallback disappears because a definition either supports a descriptor or it doesn't.

**R4 — Strip history persistence.** `rebuild()` migrates buffer history by
`stripHistoryKey(signalId, channelId, displayModeId)`
([DashboardDisplayController.cpp:55-58,1716-1743](../src/app/DashboardDisplayController.cpp)).
*De-risk:* `StripComponent` enum-ification keeps a 1:1 map to the existing mode-id strings
(step 4); the history key keeps using those strings through Stage 1. A regression test
scrubs/rebuilds/asserts history length preserved, guarding step 4.

**R5 — `static.tensor` registered-but-inert.** Risk the typed registry "helpfully" offers
it once `TensorGlyphVisualization` exists. *De-risk:* it is reachable only via
`trackedButHidden()` (never `visibleOfferable()`), its `capability()` is `{false,false,true}`
(mirroring today), and its `isAvailable()` requires `hasSceneOverlay &&
tensorGlyphGestureEnabled` with the gesture **false in Stage 1**. Matches the deliberate
omission ([:1675-1693](../src/app/DashboardDisplayController.cpp)) but now typed/explicit.
Stage 3 flips the gesture; Stage 1 must not.

**R6 — Two copies of channel logic drifting.** `canonicalModeChannel`/`modeWantsChannel`
exist in both the controller ([:75-90](../src/app/DashboardDisplayController.cpp)) and the
panel model ([DashboardPanelModel.cpp:21-36](../src/model/DashboardPanelModel.cpp)).
*De-risk:* both re-expressed through the single `StripVisualization::componentsFor` in step
4 — a strict de-duplication, the one place the migration reduces code.

**R7 — Startup validation false-positive bricking a run.** *De-risk (the corrected core
fix):* validation is **WARN-only** from introduction (step 6) through the dialog rewrite,
and is hardened to a debug `Q_ASSERT` only in step 7 **after** the hollow strings are
removed. No step — debug or release — aborts. It walks `allDescriptorList()` (#6) so it
sees what it must, and consults no availability data so it never fires on an absent-field
run.

**R8 — `addSignal` empty-mode default.** `addSignal` defaults `{}` to
`AllDisplayModes(descriptor)` ([DashboardSignalModel.cpp:512-516](../src/model/DashboardSignalModel.cpp)).
*De-risk:* documented that the honest-picker guarantee closes at step 7 (catalog peel),
not step 6; until then a caller passing `{}` could re-enable a hollow string. No Stage-1
caller is changed to pass `{}`.

---

## Summary of guarantees

- **T2 sacred:** `StripComponent::TensorT2` is first-class; no path collapses the rank-2
  tensor to a scalar. Magnitude readouts remain readouts of the preserved object.
- **Producer untouched:** zero changes outside `h5-reader/src`; no H5 write; no extraction.
  The reality-check *reports* producer gaps to the log; it does not act on them.
- **Model-is-spine:** the registry is a thin authority over the existing typed model
  (`SignalDescriptor`, `AbstractStripPanel`, `SceneRevealOverlay`); it owns offerability,
  not data.
- **Cycle-free + green:** the gate body moves to a `.cpp` (no header→registry include);
  every step (including debug startup) compiles and passes the suite; startup validation is
  WARN-only until the hollow strings are gone; the only visible change is the honest picker.
