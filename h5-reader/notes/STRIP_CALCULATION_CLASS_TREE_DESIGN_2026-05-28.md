# Strip Calculation Class Tree Design

Date: 2026-05-28

Pass: 2 of 3.

Inputs:

- `READER_INPUT_DISPLAY_INVENTORY_2026-05-28.md`
- `OBSERVABLE_MIGRATION_BRIEF_2026-05-28.md`
- current `StripCalculation`, `GeometryStrips`, `DftSigmaStrips`
- reader surfaces in `FrameNpyLoader`, `QtTrajectoryH5`, `DftShieldingStore`,
  and the typed `Qt*Group` / buffer classes

This is a design target, not the implementation patch.

## Design Goal

Build one consumer-side hierarchy for every trajectory, snapshot, or derived
signal in the data surface we already have. This is not a speculative provider
or plugin framework.

NPY snapshots, dense H5 trajectory products, ORCA/DFT frame data, topology, and
derived geometry are different source representations of one molecular
trajectory/snapshot workspace. They are not separate UI architectures. They
must all describe their available data through the same vocabulary:

- source residency
- native axis
- required anchor
- value shape
- physical units and display-unit conversion
- gap/source-attached semantics
- temporal strip modes
- static/current-frame display modes
- overlap/concept key when multiple sources represent the same physical
  measurement

The end user should see "what can I plot or display from this atom/residue/etc."
The app should know which fixed trajectory/snapshot source supplies it.

The simplifying constraint for pass 3 is explicit: the project data is always a
trajectory plus frame snapshots, with additional calculated values derived from
those sources. We should not add layers whose only justification is that a
different kind of application data might appear later.

## Scale-First Rule

Do not let the first rendered subset become the architecture. The signal source
surface is the whole known h5-reader data set from the inventory, not the six
signals that are easiest to render first.

The first implementation should separate two things:

- catalog coverage: descriptors for the known source families, axes, units,
  value shapes, and supported display modes should be broad from the start
- sampling/rendering coverage: actual frame sampling and display calculations
  can be enabled in staged slices

Unimplemented sampling is allowed to return an explicit `NotAvailable` or
`Pending` status. Missing catalog entries are not acceptable as a long-term
placeholder because they make the UI, sizing, and user workflow conform to the
small demo subset instead of the real data set.

The dialog and active signal model should be designed for hundreds of available
signals and many active signals:

- searchable/filterable Qt model/view tables, not hand-sized rows for the demo
- grouping by source family, axis, residue/atom context, value shape, and
  display mode
- stable row heights and bounded controls independent of how many signals are
  currently selected
- a strip area that scrolls/virtualizes active displays instead of assuming a
  small fixed stack
- display-mode toggles that work the same way for scalar, tensor, category,
  vector, and static overlay candidates

Implementation milestones must not make geometry/DFT/residue-dihedral the
controller shape. Early samplers can exist as reference bones behind the
calculation interface, but `DashboardDisplayController` should admit active
strip display modes generically and leave source-family knowledge to typed
recorders.

In practice:

- display controllers may route by display family: `strip.*`, `static.*`,
  `table.*`, overlay modes
- display controllers must not route by source family: geometry, DFT, H5, NPY,
  topology, electrostatics, H-bond, DSSP, etc.
- an active signal with a `strip.*` mode should create generic strip tracks
  from descriptor/channel metadata even before a sampler exists
- missing sampling is a pending/gap state on the track, not exclusion from the
  dashboard
- concrete first cases are reference implementations only if they sit behind a
  generic `StripCalculation`/recorder contract

This prevents the next implementation context from seeing the visible app and
concluding that dashboard work means adding one more geometry or shielding
branch.

## Plain Time-Series Path

Do not turn this into a broad recorder/sampler framework unless the code proves
it needs one. The dashboard strip is a retained time series for one
reader-exposed thing. The plain path is:

```text
DashboardSignal + displayMode + channel
  -> SignalBuffer
  -> small frame sample closure
  -> append value-or-gap on each frame
```

`SignalBuffer` is useful because it keeps durable frame-indexed values plus
sample/gap metadata. A virtual sampler catalog or source-family class tree is
not required for the first integration. Source-specific reading can enter as
small closures or helper functions behind this table.

The generic strip path still needs this information for every descriptor, not
just the first implemented source:

- descriptor id and concept key
- source kind and residency
- required anchor axis plus a resolver from the current focus/selection
- display mode id, with channel expansion rules for scalar, count, vector,
  tensor, per-class, embedding, rollup, and event shapes
- channel id, label, value shape, units, and display-unit conversion
- frame sample function returning value-or-gap with a gap reason
- source-attached/frame-available masks where applicable
- reveal target metadata: atom, residue, atom tuple, ring, or none
- transform eligibility: FFT/spectrum, rolling stats, histogram, range ops

A descriptor can be visible in the dashboard before its sample function exists.
That is a pending recorder, not a different UI path.

## Core Correction

The current concrete DFT and geometry strips sample their source directly. That
worked for the proof slice but should not scale.

The scalable path is:

```text
active calculations declare signal/source requirements
FrameSamplingCoordinator visits frame N
catalog prepares frame-local snapshot/DFT views if needed
calculations sample and append display values
frame-local source payloads are released
strip/static display state persists in calculation-owned buffers
```

Pass-3 scope is synchronous and single-threaded for sampling. The coordinator
may still batch work, but it should not introduce request generations,
cancellation tokens, or async machinery unless a source actually moves to a
worker thread. That complexity is justified only at a queued-signal boundary.

This prevents two actual problems in the current project:

- every strip independently loading/parsing the same source frame
- source readers becoming hidden long-lived caches of calculator payloads

Startup-loaded H5 data is still sampled through the same coordinator. It simply
has no per-frame source payload to release.

## Naming

Use "signal" for data identity and "calculation" for consumer-owned display
state.

- `SignalDescriptor`: static description of an available signal/channel.
- `SignalBinding`: descriptor plus concrete anchor selection.
- `DashboardSignal`: one user-selected signal instance. It owns the binding,
  display-mode choices, label, lifetime, and deletion state.
- `TrajectorySignalCatalog`: facade over the fixed trajectory/snapshot sources.
- `StripCalculation`: persistent temporal consumer. Owns strip buffers.
- `StaticDisplayCalculation`: persistent current-frame/static consumer. Owns
  overlay/table state, not source payload.

Avoid naming the calculation base after today's widgets. A calculation may feed
a strip, table, overlay, or event timeline, but it still belongs to the
trajectory/snapshot dashboard.

## Active Signals and Display Trees

The app shows a mutable set of `DashboardSignal` instances. A signal is the
thing the user adds or removes. It can be displayed as a strip, an atom/residue
overlay, a table, or several of those at the same time.

Keep the split clean:

1. `TrajectorySignalCatalog`: the data tree. It knows fixed trajectory/snapshot
   sources, descriptors, axes, units, gaps, reducers, and which display modes a
   signal supports. It owns no dashboard history and no VTK actors.
2. `DashboardSignalModel`: the active signal set. It owns user-selected
   bindings and their lifetime. It does not care whether a signal is currently
   rendered as strip, overlay, table, or all three.
3. Display calculation trees: peer renderer families for the active signals.
   `StripCalculation` owns temporal buffers and chart render data.
   `SpatialDisplayCalculation` / `StaticDisplayCalculation` owns current-frame
   overlay, glyph, color-map, table, and reveal/highlight state.

The strip chart tree and atom overlay tree are not children of each other.
Display instances are siblings keyed by the active signal id and display mode.
There is no special binding semantics between strip and overlay instances.
Turning an overlay or strip on/off creates or removes that display calculation;
removing the `DashboardSignal` removes all display calculations for that signal.

The UI panel/tab tree is separate layout state. It may contain strip widgets,
static tables, and overlay toggles, but it is not the signal hierarchy.

This matters for tensor/time-series data. A single descriptor such as an atom
shielding tensor time series can validly produce:

- a strip of `T0`, `T1` magnitude, `T2` component, or FFT power over frames
- a current-frame atom color map from the same reducer
- a current-frame tensor glyph/ellipsoid/principal-axis overlay
- a current-frame table using the same channel descriptors and units

The strip and overlay share the same active signal identity. Each display type
samples or renders through the catalog/coordinator path appropriate to its mode;
one display type should not read another display type's internal buffers.

## UI Replacement Boundary

The current dock/tab system is a prototype surface to learn from, not the
target architecture. The replacement UI should be:

- one always-available atom/frame inspector
- one integrated dashboard strip area fed by active signals
- one `SignalDisplayDialog` for adding/removing signals and choosing display
  modes
- scene overlays controlled by active signal display modes or explicit global
  scene-layer toggles

`SignalPickerDock`, `StripChartDock`, and `QtAtomTimeSeriesDock` should not
become permanent tab layers. Their useful pieces should be moved behind the
active signal architecture:

- `NearbySignalModel` becomes candidate-building logic for
  `SignalDisplayDialog`.
- `StripStackWidget`, `ChannelBuffer`, and current strip calculations become
  the first strip renderer/calculation components.
- current geometry, DFT sigma, residue dihedral, and FFT behavior become
  catalog-backed active signal/display modes.

The dialog mutates only `DashboardSignalModel`. It does not create charts,
overlays, or tables directly. A display controller listens to the active signal
model and creates/removes peer display calculations keyed by:

```text
DashboardSignal.id + displayModeId
```

This gives the user one list of active signals whose displays can come and go
without preserving the old tab model as a hidden permanent layer.

## Qt Ownership and Thread Rules

The catalog, coordinator, panel models, and calculations are QObjects where they
need signals, lifetime tracking, or UI ownership. Reader/source helpers do not
need to be QObjects unless they participate in Qt ownership or signals.

Rules:

- Parent calculations/models to the owning window, dock, panel model, or
  coordinator.
- Cross-object QObject references should be `QPointer`, not bare pointers.
- Signal connections must use a context object so teardown disconnects safely.
- Keep current pass-3 sampling on the GUI thread unless there is a measured
  need to move parsing work.
- If profiling forces a source read onto a worker thread, marshal immutable
  source data or final sampled values back to the GUI thread before touching Qt
  widgets or VTK scene state.
- Queued Qt signals should use Qt-friendly registered types. Prefer `quint64`
  for frame ids in signals over raw `std::size_t`.

This is not extra architecture; it is the Qt lifetime discipline needed to make
the unified hierarchy safe.

## Fixed Sources and Residency

```cpp
enum class SignalSourceKind {
    DenseH5Trajectory,
    FrameNpySnapshot,
    OrcaDftFrame,
    Topology,
    DerivedGeometry,
    SelectionEvents,
};

enum class SourceResidency {
    StartupLoaded,  // H5 time series, H5 rollups, topology, manifest
    FrameLoaded,    // per-frame NPY snapshot, ORCA/DFT frame
    Derived,        // geometry, derived ring counts, transforms
};
```

These are fixed source kinds for h5-reader, not extension points.

### Startup-Loaded Sources

Examples:

- `DenseH5TrajectorySource`
- `TopologySource`
- `H5RollupSource`
- `SelectionEventSource`

Data is already resident after run open. Sampling a frame should be cheap and
must still honor source-attached masks and gaps.

### Frame-Loaded Sources

Examples:

- `FrameNpySnapshotSource`
- `DftOrcaFrameSource`

Frame-loaded sources expose one current frame at a time. They may retain tiny
indexes, identity maps, or negative "known absent" records, but not parsed
calculator payloads as a trail.

Required source contract:

```cpp
class FrameNpySnapshotSource {
public:
    bool hasFrameSource(const FrameSampleContext& ctx) const;
    bool prepareFrame(const FrameSampleContext& ctx);
    bool frameResident(const FrameSampleContext& ctx) const;
    void releaseFrame(const FrameSampleContext& ctx);
    SignalSample sample(const DisplaySignalBinding& binding,
                        const FrameSampleContext& ctx) const;
};
```

`releaseFrame` is a semantic hook even if an early source no-ops. It keeps the
contract clear.

The DFT source follows the same shape. Do not introduce an abstract
`FrameLoadedSignalSource` base unless the two implementations first prove they
are duplicating meaningful code. If profiling forces one of these reads off the
GUI thread, add a small async adapter at that source boundary. Do not put
async/cancel complexity into the pass-3 base.

### Derived Sources

Examples:

- `GeometrySignalSource`
- `RingNeighbourDerivedSource`
- transforms such as FFT/spectrum

Derived sources read positions, topology, frame snapshots, H5 time series, or
calculation-owned buffers. They should declare those dependencies rather than
secretly reaching through the UI.

## Axes and Anchors

The old `SignalAnchorKind` is too small. Extend or replace it with:

```cpp
enum class SignalAxis {
    None,
    Atom,
    Residue,
    AtomTuple,
    Bond,
    Ring,
    AromaticRing,
    SaturatedRing,
    RingContributionPair,
    MutationMatchPair,
    Protein,
    System,
    Event,
};
```

Bindings need a typed anchor payload:

```cpp
struct AtomAnchor { std::size_t atom = 0; };
struct ResidueAnchor { std::size_t residue = 0; };
struct AtomTupleAnchor { std::vector<std::size_t> atoms; };
struct BondAnchor { std::size_t bond = 0; };
struct RingAnchor { std::size_t ring = 0; };
struct AromaticRingAnchor { std::size_t ring = 0; };
struct SaturatedRingAnchor { std::size_t ring = 0; };
struct RingContributionPairAnchor { std::size_t pair = 0; };
struct RingMembershipAnchor { std::size_t membership = 0; };
struct MutationMatchPairAnchor { std::size_t pair = 0; };
struct SystemAnchor {};
struct EventAnchor {};

using SignalAnchor = std::variant<AtomAnchor,
                                  ResidueAnchor,
                                  AtomTupleAnchor,
                                  BondAnchor,
                                  RingAnchor,
                                  AromaticRingAnchor,
                                  SaturatedRingAnchor,
                                  RingContributionPairAnchor,
                                  RingMembershipAnchor,
                                  MutationMatchPairAnchor,
                                  SystemAnchor,
                                  EventAnchor>;
```

The picker may derive anchors from a selected atom, but the binding should name
the actual axis. "Residue mode" is then just a residue-axis binding discovered
from the focused atom's residue.

## Value Shapes

```cpp
enum class SignalValueShape {
    Scalar,
    Count,
    Category,
    Vector3,
    EfgT2,
    SphericalTensor,
    TensorComponents,
    PerClassBlock,
    Embedding,
    RollupMoments,
    EventRecord,
};
```

Sampling returns a value plus status, never an ambiguous zero:

```cpp
enum class SampleStatus {
    Valid,
    Gap,              // source absent at this frame or conditional mask off
    NotAvailable,     // source/import set does not contain this signal
    Invalid,          // malformed or non-finite source value
};

enum class GapReason {
    None,
    SourceAbsent,
    FrameSourceAbsent,
    SourceMaskOff,
    AnchorUnavailable,
    NotApplicable,
    NaNSentinel,
    MalformedSource,
    Pending,
};

struct SignalSample {
    SampleStatus status = SampleStatus::Gap;
    GapReason gapReason = GapReason::None;
    SignalValueShape shape = SignalValueShape::Scalar;
    QVariant value;        // implementation can later replace with std::variant
    QString diagnostic;
};
```

For implementation, prefer typed C++ value structs over large `QVariant`
sprawl. The design point is that every sample carries status explicitly.

## Descriptors

`SignalDescriptor` is the stable entry in the picker/catalog.
Implementation should define `UnitSpec` and `ChannelDescriptor` before this
struct; they are shown in the Units section below for readability.

```cpp
struct SignalDescriptor {
    QString id;             // stable source-local signal id
    QString conceptKey;     // overlap key, e.g. "bs_shielding"
    SignalSourceKind sourceKind;
    QString importSet;      // SDK_NPY, TrajectoryH5, ORCA, Derived
    QString label;
    UnitSpec sourceUnits;
    UnitSpec defaultDisplayUnits;

    SourceResidency residency;
    SignalAxis nativeAxis;
    SignalAxis requiredAnchor;
    SignalValueShape valueShape;

    bool temporal = false;
    bool staticDisplay = false;
    bool sourceAttachedMask = false;
    bool frameLocalPayload = false;
    bool finiteScalarRequired = true;

    QStringList temporalModes;  // "scalar", "tensor.T0", "vector.magnitude"
    QStringList staticModes;    // "tensor", "vectorGlyph", "atomColor"
    QVector<ChannelDescriptor> channels;
};
```

Overlap is explicit but not forced. `conceptKey` says two descriptors represent
the same conceptual measurement or a closely related one. It does not require
the same source, sampling path, or display default.

Examples:

- `npy:bs_shielding` and `h5:bs_shielding_time_series` share concept
  `bs_shielding`.
- `h5:bs_welford` shares concept `bs_shielding.stats`, not `bs_shielding`
  directly.
- `orca_dft:total` and `npy:orca_total` share concept `orca_total`.
- `geometry:dihedral(atom_tuple)` has concept `geometry.dihedral`.

The migrated binding must include source and reducer identity, not only the
concept:

```cpp
struct DisplaySignalBinding {
    SignalSourceKind sourceKind;
    QString descriptorId;
    QString conceptKey;
    QString reducerId;
    QString displayModeId;
    SignalAnchor anchor;
    bool followsFocus = false;
};

struct DashboardSignal {
    QUuid id;
    DisplaySignalBinding binding;
    QString label;
    QStringList displayModeIds;  // e.g. "strip.tensor.T0", "static.tensor.glyph"
    bool enabled = true;
};
```

Keep an adapter from the current `SignalBinding`/`SignalKey` while migrating.
Do not rely on concept key alone when H5, NPY, and ORCA can all offer related
signals.

## Units

Units are part of the signal contract, not just axis-label text.

The SDK catalog already carries unit strings (`ppm`, `V/A`, `V/A^2`, `A`,
`radians`, `degrees`, `kJ/mol`, etc.). H5 buffers also carry display attrs for
many groups. Some imported records are mixed-unit composites, so units must be
available per descriptor, per channel, and per reducer output.

```cpp
enum class UnitDimension {
    Dimensionless,
    Length,
    Angle,
    MagneticShielding,
    ElectricField,
    ElectricFieldGradient,
    Charge,
    Energy,
    Temperature,
    Pressure,
    Volume,
    Time,
    Frequency,
    Power,
    Count,
    Tag,
};

struct UnitSpec {
    UnitDimension dimension = UnitDimension::Dimensionless;
    QString sourceSymbol;     // exact source/import symbol, e.g. "Angstrom^-3"
    QString displaySymbol;    // UI symbol, e.g. "A^-3" or "deg"
    double scaleToDisplay = 1.0;
    double offsetToDisplay = 0.0;
    bool convertible = true;
};

struct ChannelDescriptor {
    QString id;          // e.g. "temperature", "pressure_tensor.xx", "T0"
    QString label;
    SignalValueShape valueShape = SignalValueShape::Scalar;
    UnitSpec sourceUnits;
    UnitSpec defaultDisplayUnits;
};
```

Rules:

- Keep the source unit string from the reader/catalog for provenance.
- Store the display unit separately so UI labels and conversions are stable.
- Reducers must declare output units. Example: vector magnitude preserves the
  vector unit; FFT converts the x-axis to frequency and y-axis to power of the
  input unit; angle reducers should explicitly choose degrees/radians.
- Mixed-unit products such as Gromacs energy and ring-neighbourhood geometry
  must expose channel-level units. One descriptor-level unit is not enough.
- Static displays and temporal strips must both read units from the same
  descriptor/reducer path.
- Do not infer physical meaning from display text. Use `UnitDimension` and
  reducer metadata.
- Some legacy catalog strings use ASCII `A`/`Angstrom` while UI may display the
  Angstrom symbol. Keep source strings ASCII-safe in persisted IDs; display can
  use the current UI convention.
- For category/tag/count data, units are still explicit: count, tag, class, or
  dimensionless.

Open implementation choice: use a small local `UnitSpec` table first rather
than pulling in a general units library. This keeps Windows/macOS/Linux builds
simple and matches the Qt skill rule to avoid unnecessary platform-specific
dependencies.

## Trajectory Signal Catalog

Use one facade for consumer code:

```cpp
class TrajectorySignalCatalog : public QObject {
    Q_OBJECT
public:
    QVector<SignalDescriptor> descriptors() const;
    bool canSample(const DisplaySignalBinding& binding) const;

    bool prepareFrameSources(const FrameSampleContext& ctx,
                             const QVector<DisplaySignalBinding>& bindings);
    void releaseFrameSources(const FrameSampleContext& ctx);

    SignalSample sampleFrame(const DisplaySignalBinding& binding,
                             const FrameSampleContext& ctx) const;
    SignalSample readStatic(const DisplaySignalBinding& binding,
                            const FrameSampleContext& ctx) const;
    QVector<SignalSample> eventsInRange(const DisplaySignalBinding& binding,
                                        std::size_t firstFrame,
                                        std::size_t lastFrame) const;
};
```

The catalog owns or references the fixed source helpers:

- `FrameNpySnapshotSource`: frame-loaded; wraps
  `Conformation::requestSnapshot`, `QtConformationSnapshot`, and typed
  `Qt*Group` views.
- `DenseH5TrajectorySource`: startup-loaded; wraps `QtTrajectoryH5` buffers and
  source-attached masks.
- `DftOrcaFrameSource`: frame-loaded; wraps `DftShieldingStore`.
- `DerivedGeometrySource`: derived; wraps `ConformationGeometry::Measure`.
- `TopologySource`: startup-loaded/static; wraps `QtProtein`.
- `SelectionEventSource`: startup-loaded/range; wraps `QtSelectionBag` and H5
  transition/selections products.

Do not make one source class per field. The catalog should have one source
helper per actual reader/import surface, with many descriptors, then use
value-shape helpers to reduce duplication.

`readStatic` and `eventsInRange` are catalog operations because the dashboard
needs them, not because every source supports every mode. The catalog routes to
the correct fixed source and returns explicit `NotAvailable`/gap status when a
binding asks for an unsupported read shape. This avoids a plugin-style
`SignalProvider` hierarchy with broad virtual methods on every source.

Rollups, transitions, autocorrelation, and selections should not be forced
through one-frame sampling when their natural shape is static/statistical or
range/event based.

## Frame Sampling Coordinator

Responsibilities:

- own the list of active display calculations produced from
  `DashboardSignalModel`
- ask each calculation for source requirements for target frame
- prepare frame-loaded snapshot/DFT sources once per frame context
- call calculations to sample/update
- release frame-loaded payloads after all active consumers finish
- emit one UI update signal after the batch, not one storm per field

Frame identity must not be a naked `size_t`. Sources need both the H5 row and
the original trajectory frame key because current sources disagree:

- `Conformation::requestSnapshot()` currently takes the H5 row.
- DFT job directories and per-frame NPY directory names are keyed by original
  trajectory frame.
- H5 time-series buffers sample by row.

Use one context object:

```cpp
struct FrameSampleContext {
    std::size_t frameRow = 0;
    std::size_t originalFrame = 0;
    double trajectoryTimePs = 0.0;
    std::optional<double> sourceTimePs;
};
```

Coordinator sketch:

```cpp
class FrameSamplingCoordinator : public QObject {
    Q_OBJECT
public:
    void setCatalog(QPointer<TrajectorySignalCatalog> catalog);
    void setCalculations(QVector<QPointer<DisplayCalculation>> calculations);
    void visitFrame(const FrameSampleContext& ctx);

signals:
    void frameSampled(std::size_t frameRow);
};
```

The current synchronous implementation runs on the GUI thread. Final UI and
VTK mutations stay on the GUI thread per the Qt/VTK skill. If DFT/NPY parsing
later moves to workers, Qt queued-signal ordering means stale completions become
a real event-model issue; add request generations at that source adapter, not
before.

## Calculation Base Classes

### DisplayCalculation

```cpp
class DisplayCalculation : public QObject {
    Q_OBJECT
public:
    virtual QUuid signalId() const = 0;
    virtual DisplaySignalBinding binding() const = 0;
    virtual QVector<SignalSourceKind> requiredSourceKinds() const = 0;
    virtual void sampleFrame(const FrameSampleContext& ctx) = 0;

signals:
    void changed();
};
```

`StripCalculation` and `StaticDisplayCalculation` derive from this.

### StripCalculation

Already exists as the proof base. It should move toward:

```cpp
class StripCalculation : public DisplayCalculation {
public:
    virtual StripSpec spec() const = 0;
    virtual const ChannelBuffer& primaryBuffer() const = 0;
    virtual StripRenderData renderData() const = 0;
};
```

`extendToFrame` should eventually become coordinator-driven `sampleFrame`.
Direct source access inside concrete strips should be considered transitional.

### Temporal Strip Bases

Keep the class tree buffer-oriented. Do not create one C++ class for every
physics family or every tensor component.

Core temporal calculation classes:

- `SampledScalarStrip`: one scalar channel from any source/reducer.
- `SampledMultiChannelStrip`: several coordinated scalar channels sharing an
  anchor and source, such as per-class blocks or tensor components.
- `CategoryTimelineCalculation`: category/tag/event band over frames.
- `EventTimelineCalculation`: sparse event/rug strip and jump markers.
- `SpectrumTransformCalculation`: FFT/power spectrum over another sampled
  scalar calculation.

Tensor/vector/category/rollup differences should mostly live in
`ValueReducer` metadata and channel descriptors. Dedicated subclasses are
appropriate only when the sampling/display behavior is genuinely different.

### Static Display Bases

Static/current-frame display should use the same descriptors and catalog:

- `StaticTableCalculation`: scalar/tensor/per-class/system tables.
- `OverlayCalculation`: base for geometry/topology/vector/tensor overlays.
- `ColorMapCalculation`: atom/residue/bond/ring coloring by scalar/category.
- `EventListCalculation`: event lists and jump actions.
- `EmbeddingProjectionCalculation`: embedding table/projection/nearest display.

These own display state and emit changes. They do not hold source-frame payloads
after sampling.

## Value Reducers

A reducer maps a descriptor value shape to a display channel.

Examples:

- spherical tensor -> T0
- spherical tensor -> T1 magnitude
- spherical tensor -> T2 magnitude
- vector3 -> x/y/z
- vector3 -> magnitude
- EFG T2 -> component or magnitude
- per-class block -> selected class scalar
- category -> ordinal, occupancy, event, or color band
- embedding -> selected dimension or projection coordinate

Reducers should be explicit picker choices, not hidden string parsing.

```cpp
struct ValueReducer {
    QString id;        // "tensor.T0", "vector.norm", "class.tyrosine"
    QString label;
    SignalValueShape inputShape;
    SignalValueShape outputShape;  // usually Scalar for strips
    ChannelDescriptor outputChannel;
};
```

## Canonical Display Mode IDs

Use stable mode IDs so the picker and factories do not branch on labels.

Temporal examples:

- `strip.scalar`
- `strip.multiChannel`
- `strip.category.band`
- `strip.event.rug`
- `strip.vector.component`
- `strip.vector.magnitude`
- `strip.tensor.T0`
- `strip.tensor.T1.component`
- `strip.tensor.T1.magnitude`
- `strip.tensor.T2.component`
- `strip.tensor.T2.magnitude`
- `strip.efg.T2.component`
- `strip.efg.T2.magnitude`
- `strip.system.scalar`
- `strip.system.tensor.component`
- `strip.rollup.statistic`
- `strip.spectrum.power`

Static examples:

- `static.scalar.table`
- `static.scalar.colorMap`
- `static.category.label`
- `static.category.band`
- `static.vector.glyph`
- `static.tensor.table`
- `static.tensor.glyph`
- `static.geometryOverlay`
- `static.topologyOverlay`
- `static.perClass.table`
- `static.system.table`
- `static.rollup.table`
- `static.rollup.histogram`
- `static.event.list`
- `static.embedding.table`
- `static.embedding.projection`
- `static.embedding.nearest`

## Catalog and Picker

Use Qt model/view for the signal/display dialog.

- `TrajectorySignalCatalog`: aggregates descriptors from the fixed
  trajectory/snapshot sources.
- `DashboardSignalModel`: `QAbstractListModel` for the active signal set.
- `SignalCandidateModel`: `QAbstractTableModel` or `QAbstractListModel` for the
  candidate pane.
- `SignalCandidateFilterModel`: filter by current selection-derived context,
  source availability, axis, value shape, and display mode.
- `ActiveSignalDisplayModel`: optional proxy/table model for editing display
  modes per active signal if `DashboardSignalModel` needs to stay minimal.

Discovery flow:

```text
focused atom / selected atom tuple
  -> SignalRequestContext
  -> catalog.candidates(context)
  -> SignalDisplayDialog shows candidates and current active signals
  -> user chooses descriptor + reducer + display mode(s), or edits existing row
  -> DashboardSignal added to DashboardSignalModel
  -> display factories create strip/overlay/table calculations for those modes
```

Context must remain atom-centered for user workflow, but bindings must store
the actual axis.

## Factories

Factories take a `DashboardSignal`, descriptor, source kind, display mode, and
reducer and return display calculations.

The factory should be data-driven for common scalar/vector/tensor cases. Only
special cases should get dedicated subclasses:

- selected atom-tuple geometry
- event timelines
- embedding projections
- rich per-class table/stack display
- domain-specific static overlays

Most "100s of types" should not mean 100s of C++ classes. It should mean many
descriptors and reducers feeding a smaller class tree.

## Mapping From Inventory to Class Families

| Inventory family | Primary value shape | Temporal base | Static base |
| --- | --- | --- | --- |
| positions | vector3 | vector component/magnitude or derived geometry | vector/table/topology |
| topology | topology/category | event/category if time-varying; otherwise none | topology/label/overlay |
| shielding tensors | spherical tensor | tensor T0/T1/T2/magnitude | tensor table/glyph/color |
| EFG fields | EFG T2 | T2 component/magnitude | EFG table/glyph/color |
| electrostatic fields | vector3 + EFG + scalars | vector/scalar/tensor strips | vector glyph/field table |
| ring/bond per-type blocks | per-class block | per-class scalar strip | per-class table/bar |
| scalar physical values | scalar | scalar strip | scalar inspector/color |
| categories/tags | category | category/event strip | categorical label/band |
| residue angles | scalar angle/category | scalar/category strip | angle table/Rama/rotamer |
| system runtime | scalar/tensor system | system scalar/tensor strip | system dashboard table |
| rollups | rollup moments | rollup statistic strip | stats table/min-max markers |
| events/selections | event record | event timeline/rug | event list/jump action |
| embeddings | embedding | explicit dimension/projection only | projection/table/nearest |
| mutation deltas | tensor/scalar comparison | only if sequenced; otherwise static | comparison table/heatmap |

## Migration Steps

1. Add descriptor/value-shape enums and data structs beside `SignalDictionary`;
   keep compatibility with current `SignalKey`.
2. Add `TrajectorySignalCatalog` as a facade over the existing trajectory H5,
   snapshot, DFT, topology, and geometry readers. Do not create a plugin-style
   provider base first.
3. Populate the catalog broadly from the inventory: SDK/NPY catalog families,
   dense H5 product families, ORCA/DFT, topology, selections/events, and derived
   geometry. It is acceptable for many descriptors to report `Pending` sampling
   at first; it is not acceptable for the UI model to know only the first demo
   signals.
4. Add `DashboardSignalModel` as the active signal set with stable ids and
   display-mode toggles. This is the dashboard state; `AtomSelection` remains
   only focus/anchor input.
5. Implement `DenseH5TrajectorySource`, `FrameNpySnapshotSource`,
   `DftOrcaFrameSource`, and `DerivedGeometrySource` helpers with a small
   descriptor subset first.
6. Add `SignalDisplayDialog` using catalog candidates plus active signal rows.
   It replaces `SignalPickerDock` as the user path for choosing signals and
   display modes. The dialog must be model/view and scale to the complete
   descriptor set even while only some descriptors are renderable.
7. Add `FrameSamplingCoordinator` and move current geometry/DFT strips to sample
   through it.
8. Replace `StripChartDock` ownership logic with an integrated strip panel fed
   by active signal strip calculations. Keep `StripStackWidget` as the first
   renderer if it continues to fit.
9. Add generic scalar/vector/tensor strip bases.
10. Replace bespoke residue-dihedral tracks with descriptor-backed
   `ScalarStripCalculation` instances.
11. Move current static overlays/tables behind active signal display modes where
   they are signal-specific. Leave true global scene layers as global toggles.
12. Close sampling/rendering gaps family by family until every catalog
   descriptor with a meaningful display mode has an implementation.

## Non-Goals for Pass 3

- Do not rewrite `QtTrajectoryH5`.
- Do not remove existing docks.
- Do not build the complete picker UI before the descriptor/catalog path is
  proven.
- Do not implement every descriptor at once.
- Do not build source-frame caches or LRU schemes to make sampling convenient.
- Do not design for non-trajectory application data.
- Do not keep the current tab layout as a hidden permanent compatibility layer.

## Open Risks

- Dense H5 currently owns large startup-loaded payloads. The class tree should
  not make that worse, but it also should not pretend H5 is frame-loaded today.
- `SignalBinding` needs a richer anchor without creating a fragile variant mess.
- `QVariant` is convenient for Qt models but dangerous in hot sampling paths.
  Use typed reducers internally.
- Static overlays must mutate VTK/Qt scene objects only on the GUI thread.
- Source descriptors must be cheap to build and deterministic across Windows,
  macOS, and Linux.
