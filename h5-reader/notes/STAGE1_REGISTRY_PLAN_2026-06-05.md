# Stage-1 Typed Visualisation Registry — Implementation Plan (2026-06-05)

Plan only. No code was changed, no git run, no build triggered. Scope is
`h5-reader/src` only. The `nmr_shielding` library at the repo-root `src/` was
never read or referenced. No H5 write path is touched and no extraction is
triggered. This plan turns §3 of
[VISUALISATION_AUDIT_AND_MUSTHAVES_2026-06-05.md](VISUALISATION_AUDIT_AND_MUSTHAVES_2026-06-05.md)
(the `VisualizationDefinition` sketch) into a concrete, file-level Stage-1 plan.

Companion map: [DISPLAY_MODE_INVENTORY_2026-06-05.md](DISPLAY_MODE_INVENTORY_2026-06-05.md).

## Why this is the right cut, in one paragraph

The "what is a visualisation" decision today is a `QString` mode id (`strip.tensor.T2`,
`static.atomColor`, …) that flows from catalog helpers → dialog mode-kind heuristics →
panel refs → controller branch dispatch, gated by **one** string-keyed table,
`DisplayModeCapabilityFor` ([src/model/DisplayModeCapability.h:25](../src/model/DisplayModeCapability.h:25)).
That table is referenced from exactly six sites (verified by grep, listed in §2.2), which is
small enough to wrap cleanly. The audit's claim that a typed registry is "above" the existing
panel/overlay classes is correct: `AbstractStripPanel` and its subclasses are real and stay
([src/app/AbstractStripPanel.h:145](../src/app/AbstractStripPanel.h:145)); what is missing is a
typed object that owns *offerability*. The registry is that object, and it becomes the gate.

## What the code makes harder than the audit implies (flagged up front)

1. **The gate is consulted in three different shapes, not one.** The same
   `DisplayModeCapabilityFor` returns three independent bool flags
   (`hasVisibleSurface` / `buildsPanelWidget` / `emitsPanelRef`,
   [src/model/DisplayModeCapability.h:10-14](../src/model/DisplayModeCapability.h:10)), and each
   consumer reads a *different* flag: the dialog reads `hasVisibleSurface`
   ([src/app/SignalDisplayDialog.cpp:1196](../src/app/SignalDisplayDialog.cpp:1196),
   [:1407](../src/app/SignalDisplayDialog.cpp:1407)); the panel model reads `emitsPanelRef`
   ([src/model/DashboardPanelModel.cpp:86](../src/model/DashboardPanelModel.cpp:86),
   [:101](../src/model/DashboardPanelModel.cpp:101)); the controller reads `buildsPanelWidget`
   ([src/app/DashboardDisplayController.cpp:70](../src/app/DashboardDisplayController.cpp:70),
   [:1606](../src/app/DashboardDisplayController.cpp:1606)). The registry must preserve all three
   answers or panels/refs silently stop tracking. This is why Stage 1 keeps the capability struct
   as the registry's *output shape* (see §2.1) rather than replacing it wholesale.

2. **`static.tensor` deliberately carries `emitsPanelRef=true` but builds nothing**
   ([src/model/DisplayModeCapability.h:40](../src/model/DisplayModeCapability.h:40);
   omitted trigger at [src/app/DashboardDisplayController.cpp:1675-1693](../src/app/DashboardDisplayController.cpp:1675)).
   It is "registered but inert" so panel-ref *ownership* survives the deferred glyph gesture. A naive
   "registered ⇒ offerable" rule would wrongly expose it. The registry needs a tri-state per surface
   (offerable / owns-a-ref-but-not-offerable / structurally-absent), not a single bool — §1.1.

3. **Component policy is encoded as `displayModeId` string-equality in four sample functions**
   (`sampleNpyValue` [:258-270](../src/app/DashboardDisplayController.cpp:258), `sampleTensorValue`
   [:296-307](../src/app/DashboardDisplayController.cpp:296), `sampleT2Value`
   [:314](../src/app/DashboardDisplayController.cpp:314), `sampleVecValue`
   [:329](../src/app/DashboardDisplayController.cpp:329)) **and** in `canonicalModeChannel` /
   `modeWantsChannel`, which exist in **two** copies — controller
   ([src/app/DashboardDisplayController.cpp:75-90](../src/app/DashboardDisplayController.cpp:75)) and
   panel model ([src/model/DashboardPanelModel.cpp:21-36](../src/model/DashboardPanelModel.cpp:21)).
   The enum-ification (§1.3) has to land in all of them or the channel a strip plots silently
   disagrees with the channel the panel-ref tracks.

4. **`DisplayModeKind` is an anonymous-namespace enum private to the dialog**
   ([src/app/SignalDisplayDialog.cpp:56-70](../src/app/SignalDisplayDialog.cpp:56)) with its own
   string round-trip (`canonicalModeId` [:156](../src/app/SignalDisplayDialog.cpp:156),
   `modeMatchesKind` [:180](../src/app/SignalDisplayDialog.cpp:180)). It is *almost* a
   `VisualizationType` already, but it conflates UI-kind with mode-string and uses substring matching
   (`lower.contains("spectrum")`, [:188](../src/app/SignalDisplayDialog.cpp:188)). The registry's
   `VisualizationType` replaces it, but the dialog rewrite is the largest single touch (§5 step 6).

5. **The removal cascade is bidirectional and re-entrancy-guarded.** `onSignalRemoved` drops panel
   refs ([src/app/DashboardSelectionController.cpp:210-218](../src/app/DashboardSelectionController.cpp:210))
   while `onDisplayRefRemoved` removes the signal when its refcount hits zero
   ([:220-228](../src/app/DashboardSelectionController.cpp:220)), guarded by `signalsBeingRemoved_`
   ([:213](../src/app/DashboardSelectionController.cpp:213)). Refs are keyed by
   `(signalId, displayModeId, channelId)` ([DashboardDisplayRef::stableKey](../src/model/DashboardPanelModel.cpp:80)).
   **Stage 1 must keep `displayModeId` as the ref key string** so this cascade is byte-for-byte
   unchanged; the registry rides *alongside* the strings, it does not replace the ref identity yet
   (§5, incremental rule).

The net: this is a *wrap-then-peel* migration, not a rewrite. The registry becomes the single
authority that *answers* the capability question; the strings stay as opaque keys internally through
Stage 1 and are peeled later.

---

## 1. Types

New header/impl pair `src/model/VisualizationDefinition.{h,cpp}` in namespace
`h5reader::model` (it is queried by both `model/` and `app/`, and it depends only on
`SignalDescriptor` which already lives in `model/DashboardSignal.h`). It is **not** a `QObject`
(value/strategy objects, queried synchronously; no signals, no census — consistent with the existing
free-function gate). The registry singleton (§2) is the one object with lifetime.

### 1.1 `VisualizationType` and `DisplaySurface`

```cpp
enum class VisualizationType : std::uint8_t {
    TemporalStrip,     // generic ChannelBuffer strip (scalar/vector/tensor-derived)
    TensorGlyph,       // Stage 3 — rank-2 ellipsoid/superquadric scene reveal
    AtomColor,         // Stage 2 — per-atom/residue colour map
    SequenceBar,       // SequenceBarPanel
    LagCurve,          // LagDecayPanel
    ChordCoupling,     // ChordCouplingPanel
    FixedFrequency,    // FixedFreqPanel
    PowerSpectrum      // PowerSpectrumPanel (currently picker-hidden)
};

enum class DisplaySurface : std::uint8_t {
    Strip,   // appended into StripStackWidget as a TemporalStripPanel-backed track
    Panel,   // an owned AbstractStripPanel pushed via setOwnedPanels()
    Scene    // a VTK overlay reveal in MoleculeScene / SceneRevealOverlay
};
```

`VisualizationType` is the typed successor to the dialog's `DisplayModeKind`
([src/app/SignalDisplayDialog.cpp:56](../src/app/SignalDisplayDialog.cpp:56)); the audit's enum names
in §3 are adopted verbatim, with `PowerSpectrum` added because a real panel exists for it
([src/app/PowerSpectrumPanel.h:24](../src/app/PowerSpectrumPanel.h:24)) even though it stays
picker-hidden ([src/model/DisplayModeCapability.h:33](../src/model/DisplayModeCapability.h:33)).
`DisplaySurface` is the `surface()` split the audit asks for. Both get a `ToString(...)` overload
beside the existing ones in `DashboardSignal.cpp`
([src/model/DashboardSignal.cpp:24](../src/model/DashboardSignal.cpp:24)) for logging.

### 1.2 `VisualizationContext` and `VisualizationDefinition`

```cpp
// Runtime facts a definition needs to answer isAvailable() WITHOUT
// touching H5 or the library. Assembled by the controller/dialog from
// what they already hold.
struct VisualizationContext {
    const TrajectoryFieldAvailability* availability = nullptr; // canonical present/absent/empty
    bool hasTrajectory = false;                                // animated vs single pose
    bool hasDftStore   = false;                                // ORCA frame source live
    bool hasSceneOverlay = false;                              // SceneRevealOverlay wired (Scene types)
};

class VisualizationDefinition {
public:
    virtual ~VisualizationDefinition() = default;

    virtual VisualizationType type() const = 0;
    virtual QString label() const = 0;                  // dialog checkbox / filter label
    virtual DisplaySurface surface() const = 0;

    // Structural: does this definition's CLASS handle this descriptor's
    // shape/axis at all? Pure type/shape logic, no run data. This is the
    // (a)-cut in the reality-check: false => not offered, silently.
    virtual bool supports(const SignalDescriptor& descriptor) const = 0;

    // Runtime: GIVEN supports()==true, is the backing data present in
    // THIS run? Uses ctx.availability (LOCATE-BEFORE-ABSENT, §4). This is
    // the offer gate's runtime half and the input to the (b)-reality-check.
    virtual bool isAvailable(const VisualizationContext& ctx,
                             const SignalDescriptor& descriptor) const = 0;

    // The capability triplet this definition projects onto the legacy
    // gate consumers (dialog/panel-model/controller) during the wrap
    // phase. Default { offerable, false, false } for a plain strip;
    // panel/scene types override. Lets DisplayModeCapabilityFor() be
    // re-expressed as "ask the registry" without changing consumers.
    virtual DisplayModeCapability capability() const = 0;

    // The mode-id strings this definition OWNS during the incremental
    // phase (e.g. TemporalStrip tensor owns strip.tensor.T0/T1/T2/...).
    // Used by startup validation (§2.3) to prove every catalog string
    // resolves to exactly one definition. Removed in the peel phase.
    virtual QStringList legacyModeIds() const = 0;
};
```

`supports()` vs `isAvailable()` is the load-bearing split the lead's reality-check needs (§4):
`supports==false` is the *structural* cut (silent), `supports==true && isAvailable==false` is the
*expected-but-empty* case (logged). `capability()` and `legacyModeIds()` are **scaffolding** that
exist only while strings are still the ref key; they are deleted in the peel phase (§5, Stage 1.5+).

### 1.3 `StripVisualization` carries component policy as ENUMS

The strip definition's whole job is to replace `strip.tensor.T2`-style strings with a typed channel
selector. Two enums, both in `VisualizationDefinition.h`:

```cpp
// Replaces the strip.tensor.* / strip.vector.* string family.
enum class StripComponent : std::uint8_t {
    Auto,          // descriptor's natural scalar (replaces strip.scalar / first-value)
    VectorX, VectorY, VectorZ, VectorMagnitude,         // replaces strip.vector.*
    TensorT0, TensorT1, TensorT2, TensorComponent       // replaces strip.tensor.*  (T2 SACRED)
};

class StripVisualization final : public VisualizationDefinition {
public:
    DisplaySurface surface() const override { return DisplaySurface::Strip; }
    VisualizationType type() const override { return VisualizationType::TemporalStrip; }

    // The component policy this strip applies, chosen by descriptor.valueShape.
    // Returns the typed set the descriptor's channels expand into — the
    // typed successor to channelsForMode()/modeWantsChannel().
    QVector<StripComponent> componentsFor(const SignalDescriptor&) const;
};

// Typed sample request handed to the ChannelBuffer pipeline — the typed
// successor to the displayModeId string passed into sampleNpyValue() etc.
struct StripSeriesRequest {
    QString descriptorId;
    ChannelDescriptor channel;
    StripComponent component = StripComponent::Auto;
};
```

`StripComponent` is the enum the audit demands ("tensor/vector component policy as enums instead of
`strip.tensor.T2` text", audit §3). The four sample functions
([src/app/DashboardDisplayController.cpp:258-332](../src/app/DashboardDisplayController.cpp:258)) gain
an overload keyed on `StripComponent` instead of `displayModeId`; the existing string overload stays
during Stage 1 (called via a 1:1 `StripComponent ↔ mode-id` map, §5 step 4) and is deleted in the
peel. **T2 stays a first-class component (`TensorT2`) everywhere** — it is never folded into a scalar;
`tensor.T2Magnitude()` ([:303](../src/app/DashboardDisplayController.cpp:303)) is the magnitude
*readout* of the preserved rank-2 object, not a collapse, and the enum keeps that distinction
explicit (T2 sacred).

### 1.4 Concrete definitions registered in Stage 1

| Definition (class) | `type()` | `surface()` | `supports()` keys on | Backed by (existing) |
|---|---|---|---|---|
| `StripVisualization` | `TemporalStrip` | `Strip` | any descriptor with a sampleable scalar/vector/tensor shape | `buildGenericTracks` ([:2368](../src/app/DashboardDisplayController.cpp:2368)) |
| `SequenceBarVisualization` | `SequenceBar` | `Panel` | `storagePath ∈ {ired_order_parameters, reorientational_dynamics, dihedral_autocorrelation}` | `SequenceBarPanel` builders ([:1621-1655](../src/app/DashboardDisplayController.cpp:1621)) |
| `LagCurveVisualization` | `LagCurve` | `Panel` | `storagePath ∈ {kernel_dynamics, reorientational_dynamics, dihedral_autocorrelation}` | `LagDecayPanel` builders ([:1630-1661](../src/app/DashboardDisplayController.cpp:1630)) |
| `ChordCouplingVisualization` | `ChordCoupling` | `Panel` | `storagePath == kernel_coherence` | `ChordCouplingPanel` ([:1662](../src/app/DashboardDisplayController.cpp:1662)) |
| `FixedFrequencyVisualization` | `FixedFrequency` | `Panel` | `storagePath == reorientational_dynamics` (J(ω)) | `FixedFreqPanel` ([:1666](../src/app/DashboardDisplayController.cpp:1666)) |
| `PowerSpectrumVisualization` | `PowerSpectrum` | `Panel` | `storagePath == kernel_dynamics` | `PowerSpectrumPanel` ([:1626](../src/app/DashboardDisplayController.cpp:1626)); `capability()` keeps `hasVisibleSurface=false` (stays picker-hidden) |

The advanced-dynamics definitions are registered because real panels exist and the lead has not yet
ruled on them (audit "DECISIONS OWED", item 2). Registration ≠ a decision to keep them; it just makes
their existing behaviour go through the typed path. If the lead later cuts them, deleting the
definition removes them from the picker structurally (the whole point of the gate).

**Stage 2 slot:** `AtomColorVisualization : VisualizationDefinition` with `type()==AtomColor`,
`surface()==Scene`, `supports()` = per-atom/residue scalar or `TensorT2` magnitude descriptors.
Registered in Stage 2; replaces hollow `static.atomColor`
([src/model/TrajectorySignalCatalog.cpp:94](../src/model/TrajectorySignalCatalog.cpp:94)). **Do not
design its colour/legend/mapper surface here.**

**Stage 3 slot:** `TensorGlyphVisualization : VisualizationDefinition` with `type()==TensorGlyph`,
`surface()==Scene`, `supports()` = rank-2 tensor descriptors (`SphericalTensor` / `Mat3PerRow` /
`EfgT2`). Its `isAvailable()` additionally requires `ctx.hasSceneOverlay`. It drives the existing but
unwired `SceneRevealOverlay::revealTensor(tail, head, tensor[9], frame)`
([src/app/SceneRevealOverlay.h:62](../src/app/SceneRevealOverlay.h:62)); the omitted trigger block
([src/app/DashboardDisplayController.cpp:1675-1693](../src/app/DashboardDisplayController.cpp:1675))
is where its dispatch lands. **Do not design ovaloid/superquadric/ellipsoid policy here.** Stage 1
only ensures `static.tensor`'s current "registered-but-inert" status is expressed *typed* (a
registered `TensorGlyphVisualization` whose `isAvailable()` is gated on a not-yet-built gesture), so
it is honestly not-offered rather than hollow-offered.

---

## 2. Registry + gate

### 2.1 `VisualizationRegistry`

`src/model/VisualizationRegistry.{h,cpp}`, namespace `h5reader::model`. A process-wide singleton built
once at startup, holding `std::vector<std::unique_ptr<VisualizationDefinition>>`. It is the new home
of the offerability answer.

```cpp
class VisualizationRegistry {
public:
    static VisualizationRegistry& instance();          // built on first use

    // Offer query: definitions whose supports()==true for this descriptor.
    QVector<const VisualizationDefinition*> supporting(const SignalDescriptor&) const;

    // Offer + run query: supports() && isAvailable(ctx). The dialog's gate.
    QVector<const VisualizationDefinition*> offerable(const VisualizationContext&,
                                                      const SignalDescriptor&) const;

    // Wrap-phase bridge: the capability a given legacy mode-id projects.
    // Re-expresses DisplayModeCapabilityFor() as a registry lookup so the
    // six legacy consumers keep working unchanged (§2.2).
    DisplayModeCapability capabilityForMode(const QString& modeId) const;

    // Type/surface for a legacy mode-id (replaces the dialog's
    // canonicalModeId/modeMatchesKind round-trip in the peel phase).
    const VisualizationDefinition* definitionForMode(const QString& modeId) const;

    // Startup self-check (§2.3). Returns the unresolved mode-id list.
    QStringList validateAgainstCatalog(const TrajectorySignalCatalog&) const;
};
```

### 2.2 Replacing `DisplayModeCapabilityFor` as the gate

`DisplayModeCapabilityFor` ([src/model/DisplayModeCapability.h:25](../src/model/DisplayModeCapability.h:25))
is consulted at exactly these six sites (grep-verified):

1. dialog candidate-offer gate — [src/app/SignalDisplayDialog.cpp:1195](../src/app/SignalDisplayDialog.cpp:1195)
2. dialog active-signal toggle gate — [src/app/SignalDisplayDialog.cpp:1407](../src/app/SignalDisplayDialog.cpp:1407)
3. panel-ref routing — [src/model/DashboardPanelModel.cpp:86](../src/model/DashboardPanelModel.cpp:86) + [:101](../src/model/DashboardPanelModel.cpp:101)
4. controller panel-build gate — [src/app/DashboardDisplayController.cpp:70](../src/app/DashboardDisplayController.cpp:70) + [:1606](../src/app/DashboardDisplayController.cpp:1606)
5. signal-model JSON renderability — [src/model/DashboardSignalModel.cpp:411](../src/model/DashboardSignalModel.cpp:411)
6. REST JSON inventory — [src/app/RestServer.cpp:417](../src/app/RestServer.cpp:417) (via `ModeRenderabilityFor`)

**Wrap, don't rip:** the body of `DisplayModeCapabilityFor` is re-pointed to
`VisualizationRegistry::instance().capabilityForMode(mode)`. The function signature and all six call
sites stay identical, so the gate moves into the registry with a one-function change and zero consumer
churn. The `strip.` prefix rule and the six static rows currently inlined in the header
([:26-47](../src/model/DisplayModeCapability.h:26)) become the registered definitions'
`capability()` + `legacyModeIds()` outputs — same answers, typed source. This is the moment the
registry *becomes* the gate: after it, nothing can report a capability that no definition backs.

### 2.3 Startup validation rule

New `VisualizationRegistry::validateAgainstCatalog(catalog)`: for every descriptor, for every mode-id
in `AllDisplayModes(descriptor)` ([src/model/DashboardSignal.cpp:285](../src/model/DashboardSignal.cpp:285)),
assert the mode resolves to exactly one registered definition via `definitionForMode`. Any unresolved
mode-id is collected and returned.

**Where it runs:** `ReaderMainWindow` constructs the catalog and wires the dashboard
([src/main_reader.cpp:204](../src/main_reader.cpp:204)). The validation call goes in that wiring path,
right after the catalog is built and before the dialog/controller are handed it. There is **no**
existing startup self-check (grep for `qFatal`/`validateCatalog`/`Q_ASSERT` in startup found none), so
this is a new, additive hook.

**Severity:** loud. In debug builds, a non-empty unresolved list is a developer error — log each at
`qCCritical(diagnostics::cDash)` and `Q_ASSERT(unresolved.isEmpty())`. In release, log critical and
continue with those modes structurally un-offerable (the registry already won't offer them). This is a
*build/authoring* check — it catches a catalog helper advertising `static.table`
([src/model/TrajectorySignalCatalog.cpp:98](../src/model/TrajectorySignalCatalog.cpp:98)) with no
definition behind it, exactly the `static.table` / `static.atomColor` / `strip.spectrum` drift the
audit names. It is distinct from the *runtime data* reality-check (§4).

---

## 3. Stage-1 scope (and only Stage 1)

**In scope:**
- The `VisualizationDefinition` hierarchy, `VisualizationType`, `DisplaySurface`, `StripComponent`,
  `VisualizationContext`, `VisualizationRegistry` (§1, §2).
- Migrate **strips** into `StripVisualization` with enum component policy (§1.3).
- Register the existing panel definitions (SequenceBar / LagCurve / ChordCoupling / FixedFrequency /
  PowerSpectrum) so their build dispatch goes through the typed surface (§1.4) — behaviour-preserving.
- Make hollow modes structurally un-offerable: `static.table`, `static.scalar`, `static.efg`,
  `static.vector`, `static.vectorGlyph`, `static.category`, `static.per-class`, `static.rollup`,
  `static.event`, `static.system`, `static.embedding`, `static.geometry`, `static.topology`, and
  `strip.spectrum` as-a-spectrum (audit §5). They have no definition ⇒ `offerable()` returns nothing
  ⇒ the dialog can't show them. Their authoring strings are removed from the catalog helpers in the
  same step that removes the dialog's reliance on them (§5 step 6–7).
- The startup validation rule (§2.3) and the runtime reality-check (§4).

**Out of scope (named, not designed):**
- **Stage 2 — Atom Colour.** `AtomColorVisualization` slots in at §1.4's "Stage 2 slot": a new
  `surface()==Scene` definition + a new VTK colour/mapper/legend surface in `MoleculeScene`. Not now.
- **Stage 3 — Tensor glyph (ovaloid / superquadric / ellipsoid).** `TensorGlyphVisualization` slots
  in at §1.4's "Stage 3 slot" and the omitted-trigger block
  ([src/app/DashboardDisplayController.cpp:1675](../src/app/DashboardDisplayController.cpp:1675)),
  driving `SceneRevealOverlay::revealTensor` ([src/app/SceneRevealOverlay.h:62](../src/app/SceneRevealOverlay.h:62)).
  Stage 1 only makes its inert status *typed*. Not now.

Stage 1 changes **no rendering**: every pixel that draws today draws the same after, because the
panel/strip/overlay classes are untouched — only the *gate and the dispatch keys* move.

---

## 4. The empty reality-check (lead's hard requirement)

The gate must distinguish two failures and treat them differently:

- **(a) No registered renderer** — `registry.supporting(descriptor)` is empty for a mode the catalog
  named. This is a *structural* cut: the visualisation is simply not offered. **Silent is correct**
  (it is caught loudly *once*, at startup, by §2.3 — not per-descriptor at runtime).
- **(b) Registered and `supports()==true`, but the data is absent/empty at runtime** —
  `isAvailable(ctx, descriptor)==false`. This is **not** silently hidden. It is the producer/data-
  integrity signal the reader exists to surface ("we expected data for X here and it was empty"). It
  is **logged**.

### Where the check lives

The runtime reality-check belongs in `DashboardDisplayController` because that is where the catalog,
the live availability, and the per-descriptor offer decision already meet, and where the existing
`smokeSummary` data-integrity harness already lives
([src/app/DashboardDisplayController.h:35-90](../src/app/DashboardDisplayController.h:35),
`smokeSummary()` consumed by the dock and REST). Concretely:

- A new method `DashboardDisplayController::collectExpectedButEmpty()` walks the catalog's descriptors,
  and for each registered definition with `supports(descriptor)==true` evaluates
  `isAvailable(ctx, descriptor)`. Each `supports && !isAvailable` becomes one `ExpectedButEmpty{
  descriptorId, storagePath, visualizationType, canonicalName, alternatesTried, availabilityState }`
  record.
- The records are emitted **to the UDP/structured log** under a distinct event tag
  (`event=viz_expected_but_empty …`) via `qCWarning(diagnostics::cDash)`, the project's primary debug
  channel, and folded into the existing `DashboardSmokeSummary` as a new vector so the REST inventory
  and any report surface pick them up alongside the current sparseness stats (the summary already
  carries per-series gap accounting — this is the descriptor-level analogue).
- It runs once per context change (catalog/availability set), **not** per frame — same cadence as
  `rebuild()`'s structural work, never on a `setFrame()` tick (the one-QTimer discipline:
  `frameChanged` subscribers don't re-scan integrity).

### LOCATE-BEFORE-ABSENT (honoured)

`isAvailable()` must not declare "empty" from a single canonical-name miss. The availability table
already keys by descriptor id *and* by storage path
([src/model/TrajectoryFieldAvailability.h:118-140](../src/model/TrajectoryFieldAvailability.h:118),
`recordForDescriptor` / `byStorage_`), and distinguishes `Absent` / `NoFramePayload` / `AllMissing` /
`AllZeroStructural` / `AllZeroObserved` ([:142-158](../src/model/TrajectoryFieldAvailability.h:142)).
The rule:

1. Look up the descriptor's **canonical** storage path first (`recordForDescriptor`).
2. If no record, that is "not catalogued for this run" — **the absence of a record means "assume
   present"** today ([:130-138](../src/model/TrajectoryFieldAvailability.h:130) returns `true`), so a
   bare miss is *not* an empty signal; do not log it as one.
3. Only `state ∈ {Absent, NoFramePayload, AllMissing}` — a *positive* statement that the source was
   scanned and found empty — qualifies as "expected-but-empty". `AllZeroObserved` is data, not
   absence; `AllZeroStructural` is a separate (structural-zero) note, logged distinctly.
4. The `ExpectedButEmpty` record carries the canonical name **and** the `alternatesTried` list (e.g.
   for tensor descriptors that may live under a per-component vs packed path) so the log line proves
   the locate step happened. "Absent under the canonical name" is never reported as "not in the data"
   without the alternate having been checked.

This keeps the check honest in both directions: it does not cry wolf on descriptors a given run simply
doesn't include, and it does not swallow a genuinely-empty field the producer was expected to fill.

---

## 5. Migration approach — incremental, each step buildable and green

The migration is **wrap-then-peel**, and the invariant through Stage 1 is: **the `displayModeId`
string stays the panel-ref key and the on-the-wire identity** (so the removal cascade and persistence
are untouched, §6). The registry rides alongside the strings as the *authority*, then the strings are
peeled in a later stage.

### File-by-file touch list

New files (additive; register in the flat source list at
[CMakeLists.txt:100-220](../CMakeLists.txt:100), beside the other `model/` and `app/` entries):
- `src/model/VisualizationDefinition.h` / `.cpp` — types (§1).
- `src/model/VisualizationRegistry.h` / `.cpp` — registry + gate (§2).
- `src/model/StripVisualization.{h,cpp}`, `SequenceBarVisualization.{h,cpp}`,
  `LagCurveVisualization.{h,cpp}`, `ChordCouplingVisualization.{h,cpp}`,
  `FixedFrequencyVisualization.{h,cpp}`, `PowerSpectrumVisualization.{h,cpp}` — the concrete
  definitions (may share one `.cpp` if the lead prefers fewer files; one-per-type matches the existing
  one-panel-per-file convention, [CMakeLists.txt:206-216](../CMakeLists.txt:206)).

Modified files:
- `src/model/DisplayModeCapability.h` — body of `DisplayModeCapabilityFor` re-points to the registry
  (§2.2). Signature unchanged.
- `src/model/DashboardSignal.cpp` — add `ToString(VisualizationType)` / `ToString(DisplaySurface)`
  beside the existing `ToString` family ([:24](../src/model/DashboardSignal.cpp:24)).
- `src/app/DashboardDisplayController.{h,cpp}` — (i) panel-build dispatch
  ([:1605-1673](../src/app/DashboardDisplayController.cpp:1605)) asks the registry which
  `surface()==Panel` definition matches and calls *its* builder, instead of the
  `path == … && mode == …` ladder; (ii) the four sample functions
  ([:258-332](../src/app/DashboardDisplayController.cpp:258)) gain `StripComponent` overloads;
  (iii) `channelsForMode` / `modeWantsChannel` / `canonicalModeChannel`
  ([:75-90](../src/app/DashboardDisplayController.cpp:75), [:2416](../src/app/DashboardDisplayController.cpp:2416))
  re-expressed via `StripVisualization::componentsFor`; (iv) add `collectExpectedButEmpty()` + the
  `ExpectedButEmpty` summary field (§4).
- `src/model/DashboardPanelModel.cpp` — `canonicalModeChannel` / `modeWantsChannel`
  ([:21-36](../src/model/DashboardPanelModel.cpp:21)) and the `emitsPanelRef` branch
  ([:101](../src/model/DashboardPanelModel.cpp:101)) re-expressed via the registry, keeping the
  `(signalId, mode, channelId)` ref shape ([:102](../src/model/DashboardPanelModel.cpp:102)) **exactly**.
- `src/app/SignalDisplayDialog.cpp` — the big one: replace the anonymous `DisplayModeKind` enum and
  its `canonicalModeId` / `modeMatchesKind` / `allModeKinds` / `modeForKind`
  ([:56-222](../src/app/SignalDisplayDialog.cpp:56)) with `VisualizationType` + registry queries; the
  offer gates ([:1195-1208](../src/app/SignalDisplayDialog.cpp:1195),
  [:1399-1416](../src/app/SignalDisplayDialog.cpp:1399)) call `registry.offerable(ctx, descriptor)`;
  the candidate/active checkbox sets ([:805-810](../src/app/SignalDisplayDialog.cpp:805),
  [:847-852](../src/app/SignalDisplayDialog.cpp:847)) and the mode filter
  ([:874-876](../src/app/SignalDisplayDialog.cpp:874)) are built from registered definitions.
- `src/model/TrajectorySignalCatalog.cpp` — remove the hollow mode strings from the helper lists
  ([:84-206](../src/model/TrajectorySignalCatalog.cpp:84)) (done **last**, after the dialog no longer
  reads them, so each step stays green).
- `src/main_reader.cpp` / `src/app/ReaderMainWindow.cpp` — call
  `registry.validateAgainstCatalog(catalog)` in the catalog-wiring path
  ([src/main_reader.cpp:204](../src/main_reader.cpp:204)).
- `src/app/RestServer.cpp` ([:417](../src/app/RestServer.cpp:417)) and
  `src/model/DashboardSignalModel.cpp` ([:411](../src/model/DashboardSignalModel.cpp:411)) — JSON
  renderability flags now sourced from the registry-backed capability (mechanical; the
  `ModeRenderabilityFor` body changes, the JSON shape does not).

### Ordered build steps (each independently buildable + green)

1. **Add types, no wiring.** Land `VisualizationDefinition.h/.cpp`, the enums, `VisualizationContext`,
   and the `ToString` overloads. Nothing consumes them yet. Builds; behaviour identical.
2. **Add the registry + concrete definitions, populated to mirror today's table.** Each definition's
   `capability()` + `legacyModeIds()` reproduces the exact rows/prefix in
   [DisplayModeCapability.h:26-47](../src/model/DisplayModeCapability.h:26). Add unit coverage:
   `capabilityForMode(mode) == DisplayModeCapabilityFor(mode)` for the full mode-id set in the catalog
   (a golden equivalence test). Still nothing wired; builds green.
3. **Flip `DisplayModeCapabilityFor` to delegate to the registry** (§2.2). All six consumers now read
   registry-backed answers through the unchanged function. The step-2 equivalence test guarantees
   byte-identical behaviour. Builds; UI unchanged.
4. **Enum-ify strip component policy** inside the controller + panel-model, behind a 1:1
   `StripComponent ↔ mode-id` map so output is identical. The string sample overloads still exist and
   are called through the map. Builds; strips render identically.
5. **Route panel-build dispatch through `surface()==Panel` definitions** instead of the
   `path && mode` ladder ([:1621-1671](../src/app/DashboardDisplayController.cpp:1621)). Same builders,
   selected by definition. Builds; panels render identically.
6. **Rewrite the dialog onto `VisualizationType` + `registry.offerable()`.** This is where hollow
   modes stop being offered (no definition ⇒ not in the checkbox set). The startup validation (§2.3)
   and runtime reality-check (§4) land here. Builds; the picker now shows only backed visualisations.
7. **Peel the strings from the catalog helpers** ([:84-206](../src/model/TrajectorySignalCatalog.cpp:84)):
   delete `static.table` / `static.atomColor` / `strip.spectrum` / the dead `static.*` family. Because
   step 6 already stopped the dialog reading them, this is dead-string removal. Re-run startup
   validation: clean. Builds green.
8. **(Stage 1.5, optional within Stage 1) begin the peel of `legacyModeIds()`/`capability()`
   scaffolding** once no consumer needs the string round-trip. This is the boundary into the typed-ref
   future (Stage 2+ territory); stop here for Stage 1 if the lead prefers.

Steps 1–3 are pure infrastructure (no UI change). Steps 4–5 are behaviour-preserving dispatch swaps.
Step 6 is the visible change (honest picker). Step 7 is cleanup. Each compiles and passes the existing
suite on its own.

---

## 6. Risks and how the plan de-risks each

**R1 — The `DashboardSelectionController` removal cascade.** Bidirectional and re-entrancy-guarded:
`onSignalRemoved` → `removeDisplayRefsForSignal` (guarded by `signalsBeingRemoved_`,
[src/app/DashboardSelectionController.cpp:210-218](../src/app/DashboardSelectionController.cpp:210)),
`onDisplayRefRemoved` → `removeSignal` when refcount hits zero
([:220-228](../src/app/DashboardSelectionController.cpp:220)). *De-risk:* Stage 1 does **not** touch
the ref identity. Refs stay `(signalId, displayModeId, channelId)` strings
([DashboardDisplayRef](../src/model/DashboardPanelModel.cpp:80)); `DisplayRefsForSignal`
([:89](../src/model/DashboardPanelModel.cpp:89)) keeps emitting the same refs; the registry only
changes *which capability flag* a mode reports. The cascade code is in the no-touch set. Typed refs
are explicitly deferred past Stage 1 (step 8 boundary).

**R2 — The panel-ref path.** The active-panel filter compares a literal
`DashboardDisplayRef{signal.id, mode, "panel"}` in `rebuild()`
([src/app/DashboardDisplayController.cpp:1561-1565](../src/app/DashboardDisplayController.cpp:1561),
[:1616-1619](../src/app/DashboardDisplayController.cpp:1616)), and `emitsPanelRef` decides ref
emission ([src/model/DashboardPanelModel.cpp:101](../src/model/DashboardPanelModel.cpp:101)).
*De-risk:* `capability().emitsPanelRef` is preserved per definition (the panel definitions set it
true, `static.tensor`'s `TensorGlyphVisualization` keeps it true-but-inert, §1, flag #2). Step 5 swaps
the *builder selection* but keeps the `"panel"` sentinel ref shape, so the filter comparison is
unchanged. The step-2 equivalence test pins this.

**R3 — The dialog mode-kind heuristics.** Substring matching (`lower.contains("spectrum")`,
[src/app/SignalDisplayDialog.cpp:188-196](../src/app/SignalDisplayDialog.cpp:188)) and the default-
checked walk that deliberately orders Table last
([:360-377](../src/app/SignalDisplayDialog.cpp:360)) are brittle and the largest rewrite. *De-risk:*
sequence it **last** (step 6), after the registry is already the gate (step 3) and dispatch is typed
(steps 4–5), so the dialog rewrite is a presentation swap over a verified-stable backend. The
`offerable()` query replaces the `hasVisibleSurface && descriptorAvailable` compound
([:1199](../src/app/SignalDisplayDialog.cpp:1199)) with one typed call, removing the substring
matching entirely. The `modeForKind` fallback chain ([:209-218](../src/app/SignalDisplayDialog.cpp:209))
disappears because a definition either supports a descriptor or it doesn't.

**R4 — Strip history persistence.** `rebuild()` migrates per-series buffer history across rebuilds by a
string key `stripHistoryKey(signalId, channelId, displayModeId)`
([src/app/DashboardDisplayController.cpp:55-58](../src/app/DashboardDisplayController.cpp:55),
[:1716-1743](../src/app/DashboardDisplayController.cpp:1716)), and `SignalChannelKey::stableId()` feeds
the channel id ([:2394](../src/app/DashboardDisplayController.cpp:2394)). If the displayModeId or
channel id string changes under a series, its accumulated history is dropped. *De-risk:* the
`StripComponent` enum-ification keeps a 1:1 map to the existing mode-id strings (step 4), and the
history key continues to use those strings through Stage 1, so keys are stable across the migration.
The peel (step 8+) that would change the key is out of Stage-1 scope. A regression test that scrubs,
rebuilds, and asserts history length is preserved guards step 4.

**R5 — `static.tensor` is registered-but-inert (flag #2).** Risk that the typed registry "helpfully"
offers it once a `TensorGlyphVisualization` exists. *De-risk:* its `isAvailable()` is gated on
`ctx.hasSceneOverlay` **and** a Stage-3 gesture that does not exist yet, so `offerable()` returns false
and the dialog does not show it — matching today's deliberate omission
([:1675-1693](../src/app/DashboardDisplayController.cpp:1675)) but now *typed and explicit* instead of
a comment. Stage 3 flips the gesture flag; Stage 1 must not.

**R6 — Two copies of channel logic drifting.** `canonicalModeChannel`/`modeWantsChannel` exist in both
the controller ([:75-90](../src/app/DashboardDisplayController.cpp:75)) and the panel model
([src/model/DashboardPanelModel.cpp:21-36](../src/model/DashboardPanelModel.cpp:21)). *De-risk:* both
are re-expressed through the *single* `StripVisualization::componentsFor` in step 4, collapsing the
duplication — a strict improvement, and the one place the migration reduces rather than preserves code.

**R7 — Startup validation false-positive halting a real run.** A descriptor legitimately advertising a
mode whose definition isn't registered yet (mid-migration) would assert. *De-risk:* validation is
introduced in step 6, *after* every kept mode already has a definition (steps 2–5) and *before* the
hollow strings are removed (step 7), so the only unresolved modes at that point are precisely the
hollow ones being cut — the assert firing is the *intended* signal, and step 7 clears it. Release
builds log-and-continue rather than abort (§2.3).

---

## Summary of guarantees

- **T2 sacred:** `StripComponent::TensorT2` is a first-class channel; no path collapses the rank-2
  tensor to a scalar. Magnitude readouts remain readouts of the preserved object.
- **Producer untouched:** zero changes outside `h5-reader/src`; no H5 write; no extraction. The
  reality-check *reports* producer gaps to the log, it does not act on them.
- **Model-is-spine:** the registry is a thin authority over the existing typed model
  (`SignalDescriptor`, `AbstractStripPanel`, `SceneRevealOverlay`); it owns offerability, not data.
- **Incremental + green:** eight ordered steps, each compiling and passing the suite; the only visible
  change is at step 6 (an honest picker), backed by a verified-stable typed gate.
