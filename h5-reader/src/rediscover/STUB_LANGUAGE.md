# Rediscover — STUB language (Phase 1 waterfall draft, 2026-05-31)

Status: **historical stub draft**; the current C++ surface is partly built in `RelationshipEngine`, `ComposedRelationships`, `ResidentIndexes`, `PerAtomSubstrate`, and `CaseHunter`. Positive 1P9J leave-atoms-out/between oracle claims below are superseded by the June 4 true-LOAO retraction.
Branch `h5-reader-pysr-spike`. This is the
complete-today skeleton of the multi-scenario surface specified in
`SURFACE_DESIGN.md`: the catalog/ArraySpec, the primitive verbs, the
combinator types, the `Relationship` bundle, the engine loop, the JSON
scenario-spec + manifest, a worked example for every one of the nine
relationship items, and a surfaced-issues list.

Discipline marker: **there is only today, no deferral.** Every item below
has a full stub + a worked example as a bundle of curried closures. Nothing
is "exploratory / TBD / later." Where an item does not fit the language, it
is written up as an ISSUE in the final section — a finding to fix, not a
punt.

This is a **stub**: signatures, declarations, schemas; the later built bodies live in the C++ files named above. It reads with
`SURFACE_DESIGN.md` (the design), `DESIGN.md` / `STATE.md` / `FINDINGS.md` /
`PLAN.md` (the one-off + the oracle numbers), and
`spec/substrate_conventions_2026-05-30.md` (the conventions). It does NOT
re-open the decided design; it stubs it.

What was verified by reading vs. inferred is called out per section. All the
one-off code (`RediscoverTypes.h`, `RunData.h`, `LocalFrameBasis.h`,
`ExtractionSupport.h`, `RecordSink.h`, `FrameSpatialIndex.h`,
`SphericalBasis.h`, `RingCurrentNeighborhood.cpp`,
`McConnellNeighborhood.cpp`, `main_extract.cpp`), the H5 boundary
(`io/QtTrajectoryH5.h`, `model/QtTimeSeriesBuffers.h`,
`model/QtSpecialBuffers.h`), the typed model (`model/QtAtom.h`,
`model/QtRing.h`, `model/Types.h`, `model/ConformationGeometry.h`), and the
diagnostics (`diagnostics/StructuredLogger.h`) were read in full for this
draft. The stub reuses their real names.

The shape is preserved exactly: **iterated closures / currying** is the
heart; relationships are **named bundles defined in code**, not a plugin
kit; **reader owns H5** (the body loads once, resident, immutable; Python
consumes only the emitted CSV/manifest); the interface is a **batch file**
(CLI in, `.log` + atomic CSV + `manifest.json` out).

---

## Part (a) — the stub C++ language

Namespace `h5reader::rediscover` throughout (the one-off's namespace). Plain
data classes (no `QObject`) except where a Qt primitive is the right tool
(`QSaveFile`, `QJsonDocument`, `QCommandLineParser`, `StructuredLogger`).

### a.0 — the resident immutable body

`RunData` (from `RunData.h`, read in full) is the resident body, loaded once
by `RunLoader::Load` and read-only thereafter. The language wraps it: the
catalog, the indexes, and the verbs all capture `const RunData&`. No
re-reads, no Python H5 path.

```cpp
// Already exists (RunData.h) — the body the whole language captures.
struct RunData {
    std::unique_ptr<model::QtProtein>     protein;        // identity + topology spine
    std::unique_ptr<model::Conformation>  conformation;   // owns the H5 trajectory
    io::CalcsetManifest                   manifest;
    DftFrameSet                           dft;            // sparse (atom,orig)->tensor
    FrameMap                              frameMap;       // H5 row -> orig -> dftRows()
    const model::TrajectoryConformation*  trajectory() const;  // typed downcast
    const io::QtTrajectoryH5*             h5() const;
};
```

The stub adds one resident sidecar — the day-1 index set (a.2) — built at
load, held alongside `RunData`. To keep `RunData` immutable and additive, the
indexes are a separate const carrier the engine receives next to the body
(NOT a mutable field bolted onto `RunData`):

```cpp
// NEW. Resident day-1 indexes, built once after RunLoader::Load, held const.
// Bundled with the body for the engine; never mutated during traversal.
struct ResidentIndexes;   // declared in a.2

struct Body {                          // the const pair the verbs/engine see
    const RunData&         run;
    const ResidentIndexes& idx;
    const Catalog&         catalog;    // a.1
};
```

> Naming note: `SURFACE_DESIGN.md` says "the body" and `at(atom,frame)`
> verbs take `Body`. `RunData` is the spine inside it. Stub keeps both: `Body`
> is the const triple {run, idx, catalog}; `RunData` is unchanged.

### a.1 — the data catalog (one ArraySpec per array; uniform `value()`)

`SURFACE_DESIGN.md` Layer 0. One entry per H5 dataset / NPY, so the access
language covers all of them uniformly. Mirrors
`python/nmr_extract/_catalog.py`'s one-`ArraySpec`-per-NPY discipline, on the
C++ body. The accessor closures bridge each `ArrayId` to its concrete typed
buffer in `QtTrajectoryH5` (verified: all source buffers below are accessors
on `io/QtTrajectoryH5.h`).

```cpp
// Typed id per dataset/NPY. One enumerator per resident array. Adding a
// datum is adding an enumerator + a registry line, NOT new grammar.
enum class ArrayId : int {
    Positions = 0,            // QtPositionsTimeSeries        (N,T,3)
    KernelBs,                 // QtShieldingTimeSeries bs     (N,T,9)  ring-current kernel
    KernelMc,                 // QtShieldingTimeSeries mc     (N,T,9)  McConnell kernel
    RingNeighbourhood,        // QtRingNeighbourhoodTimeSeries(N,T,11,4) + membership(N,11)
    ApbsEfg,                  // QtT2TimeSeries apbs_efg       (N,T,5)  EFG T2
    ApbsEfield,               // QtVec3TimeSeries apbs_efield  (N,T,3)  E-field l=1
    Aimnet2Charge,            // QtScalarTimeSeries            (N,T)    per-atom charge
    Aimnet2ChargeRespScalar,  // QtAimnet2ChargeResponseGradientTimeSeries.scalar (N,T)
    Aimnet2ChargeRespVector,  // QtAimnet2ChargeResponseGradientTimeSeries.vec    (N,T,3)
    Aimnet2Embedding,         // QtEmbeddingTimeSeries         (N,T,256) float32
    Ff14sbCharge,             // *** NOT in H5 — from prmtop/topol.top. See ISSUE C1.
    MopacCharge,              // *** NOT present in this calcset. See ISSUE C2.
    DftTotalRaw,              // sparse (atom,orig)->Mat3 total (DftFrameSet)
    DftDiaRaw,                // sparse                        diagnostic
    DftParaRaw,               // sparse                        diagnostic
};

// What axes an array is addressed on. Drives which value() overload is legal.
struct AxisSpec {
    bool atom  = false;
    bool frame = false;   // dense time-series; false for static (prmtop) / Welford
    bool slot  = false;   // ring_nbhd slot axis
    bool comp  = false;   // component axis (3 for Vec3, 5 for T2, 9 for tensor, 256 emb)
    int  comp_count = 0;  // 0 = scalar
};

enum class ArrayRank : int { Scalar, Vec3, T2_5, Tensor9, Embedding256, RingNbhd4 };
enum class ArrayDType : int { F64, F32, I32 };

// Where a sparse (DFT) or absent (mopac) array lives, so value() can answer
// "present?" honestly per (atom,frame) rather than faking a zero.
enum class ArrayResidence : int { DenseH5, StaticPrmtop, SparseDftByOriginal, Absent };

struct ArraySpec {
    ArrayId        id;
    QString        name;        // catalog key, == column-name stem in CSV
    ArrayRank      rank;
    AxisSpec       axes;
    ArrayDType     dtype;
    ArrayResidence residence;
    QString        unit;        // "" for dimensionless; e.g. "e", "e/Angstrom^2"
    // Bridge closure: captures Body, resolves (atom,frame[,slot][,comp]) to the
    // concrete typed buffer accessor. Set at registry build. NOT a body itself —
    // the stub declares the slot; the registry (a.1, no body here) wires it.
    std::function<double(const Body&, std::size_t atom, std::size_t frame,
                         int slot, int comp)> read_scalar;   // for Scalar/comp reads
    // present? — honest "no observation" answer (sparse DFT, absent mopac).
    std::function<bool(const Body&, std::size_t atom, std::size_t frame)> present;
};

// The catalog: every resident array, one ArraySpec each. Built day 1.
class Catalog {
public:
    const ArraySpec& spec(ArrayId id) const;          // typed lookup
    bool             has(ArrayId id) const;            // false for Absent (mopac)
    // ── the ONE uniform addressing verb (SURFACE_DESIGN Layer 0) ──
    double value(const Body&, ArrayId, std::size_t atom, std::size_t frame,
                 int slot = -1, int comp = -1) const;  // scalar | one component
    Vec3   valueVec3(const Body&, ArrayId, std::size_t atom, std::size_t frame) const;
    std::array<double,5> valueT2(const Body&, ArrayId, std::size_t atom, std::size_t frame) const;
    Mat3   valueTensor(const Body&, ArrayId, std::size_t atom, std::size_t frame) const;
    const float* valueEmbedding(const Body&, ArrayId, std::size_t atom,
                                std::size_t frame, std::size_t& n_dims_out) const;
    bool   present(const Body&, ArrayId, std::size_t atom, std::size_t frame) const;
private:
    std::vector<ArraySpec> specs_;  // indexed by ArrayId ordinal
};
```

### a.2 — day-1 indexes

`SURFACE_DESIGN.md` Layer 0 + "Additional indexing tools." Built at load, not
lazily — scoping must be O(log n) from frame 0. The one-off built the
McConnell KD-tree lazily per frame (`FrameSpatialIndex`, verified); the
language generalizes to per-cloud, per-frame, all built up front.

```cpp
// The four source-cloud kinds the verbs scope over (SURFACE_DESIGN Layer 0).
enum class CloudKind : int { Atoms = 0, BondMidpoints = 1, RingCenters = 2, ChargeSites = 3 };

// One source's per-frame point + identity, generic over cloud kind. Generalizes
// the one-off AnisoBond (FrameSpatialIndex.h) to all clouds. A SourceRef is the
// handle a selector returns; the attacher resolves it to geometry/identity.
struct SourceRef {
    CloudKind kind = CloudKind::Atoms;
    int32_t   cloud_index = -1;   // row in the cloud's per-frame point list
    int32_t   entity_index = -1;  // atom idx / bond idx / ring idx / charge-site atom idx
};

// Per-cloud, per-frame KD-tree. Generalizes FrameSpatialIndex (bond midpoints)
// to {atoms, bond-midpoints, ring-centers, charge-sites}. nanoflann L2_Simple,
// dim 3, std::size_t index (verified shape from FrameSpatialIndex.h).
class CloudTree {
public:
    // cutoff REQUIRED, recorded (substrate_conventions: no default cutoffs).
    std::vector<SourceRef> within(const Vec3& query, double cutoff) const;
    const Vec3& pointAt(std::size_t cloud_index) const;
    std::size_t size() const;
};

// The resident index set, built day 1 (SURFACE_DESIGN: "default is day-1").
struct ResidentIndexes {
    // per-cloud per-frame trees: trees[cloud][frame]
    std::array<std::vector<CloudTree>, 4> trees;

    // typed-atom index (residue, locant[, branch, element]) -> atom. Generalizes
    // QtResidue's N/CA/C/O/H/HA/CB cache to ALL locants. Powers atomOf(). The
    // identity-from-chemistry index — no name-scan, no positional guess.
    // (NEW — does NOT exist in the one-off; see ISSUE I1.)
    int32_t atomOf(int32_t residue_index, model::Locant,
                   model::BranchAddress = {}, model::Element = model::Element::Unknown) const;

    // own-ring / own-atom set per aromatic-H (for is_self_or_bonded). The one-off
    // rebuilt this per run (RingCurrentNeighborhood.cpp ownRingsByAtom); day-1 it.
    const std::unordered_set<int32_t>& ownRings(std::size_t atom) const;
    const std::unordered_set<int32_t>& ownAtoms(std::size_t atom) const;

    // slot-list backend is the H5 ring_membership_per_atom itself — no build,
    // named as the slots() source backend (SURFACE_DESIGN item 4).

    // charge-site cloud membership: which atoms are charge sites for a given
    // ChargeSource (drives CloudKind::ChargeSites). See ISSUE C1/C3.
    const std::vector<int32_t>& chargeSites() const;
};

// Built once after RunLoader::Load. NOT during traversal.
ResidentIndexes BuildResidentIndexes(const RunData& run);
```

### a.3 — primitive access verbs (Layer 1; pure; return views)

`SURFACE_DESIGN.md` Layer 1. These are the little APIs. In the **iterated-
closure** shape they are mostly *closure factories*: `near(cloud,cutoff)`
returns a configured `(atom,frame)->[SourceRef]` closure. The free-function
forms below are the underlying primitives the factories close over.

```cpp
// AtomState — pos + kernel handles at (atom,frame). Returned by at().
struct AtomState {
    Vec3            pos;
    SphericalTensor bare_kernel;   // this relationship's cross-check kernel (bs/mc/...)
    bool            kernel_present = false;
};

// FrameSpan — [frame-w .. frame] for Δ / rate-of-change. Atom-major storage
// makes this a contiguous span (SURFACE_DESIGN item 5).
struct FrameSpan { std::size_t atom; std::size_t frame_lo; std::size_t frame_hi; };

AtomState at      (const Body&, std::size_t atom, std::size_t frame);
FrameSpan window  (const Body&, std::size_t atom, std::size_t frame, std::size_t w);
std::vector<SourceRef> neighbors(const Body&, std::size_t atom, std::size_t frame,
                                 CloudKind, double cutoff);   // cutoff REQUIRED
// value(...) is Catalog::value (a.1) — the uniform catalog read.

// Geometry primitives the attachers close over (verified names from
// ConformationGeometry.h + the one-off):
Vec3        pos      (const Body&, std::size_t atom, std::size_t frame);  // conf.atomPosition
model::RingGeometry ringGeom(const Body&, std::size_t ring, std::size_t frame);  // RingGeometryAt
std::vector<int> ringsOf(const Body&, std::size_t atom);  // ownAromaticRing / memberships
```

### a.4 — combinator types (Layer 2; the little-API closures)

`SURFACE_DESIGN.md` Layer 2, plus `SourceClassifier` and `Attacher` named
explicitly per the deliverable. These are the closure *types*; the factories
(a.5) produce values of these types with config + body baked in (currying).

```cpp
// All take a const Body& (the resident body) + the (atom,frame) coordinate; the
// config is captured by the factory. std::function at the cool levels;
// templating the innermost Attacher is reserved for a profiling bite
// (SURFACE_DESIGN "Performance notes").

using Stratum        = std::function<std::vector<std::size_t>(const Body&)>;
using LocalFrameFn   = std::function<LocalFrame(const Body&, std::size_t atom, std::size_t frame)>;
using SourceSelector = std::function<std::vector<SourceRef>(const Body&, std::size_t atom, std::size_t frame)>;
//   cutoff/cloud are CAPTURED by the factory (near(cloud,cutoff)); the closure
//   signature is cutoff-free — the cutoff is baked in + recorded once.
using WindowFn       = std::function<FrameSpan(const Body&, std::size_t atom, std::size_t frame)>;  // w captured

// Classifier: validity + a typing key for per-key reduction. Generalizes the
// one-off's is_self_or_bonded + ring_type / bond_category.
struct SourceClass { bool valid = true; int key = -1; };
using SourceClassifier = std::function<SourceClass(const Body&, std::size_t atom,
                                                   std::size_t frame, const SourceRef&)>;

// Attacher: produce the per-source values, expressed in the LocalFrame. The
// frame closure is CAPTURED (disp_in(frameF)). Writes into a SourceSlot — the
// one-off's emitted per-source primitive (RediscoverTypes.h). One Attacher may
// fill geometry, another identity, another the kernel cross-check; a relationship
// supplies a list (SURFACE_DESIGN: attachers[]).
using Attacher = std::function<void(const Body&, std::size_t atom, std::size_t frame,
                                    const SourceRef&, const LocalFrame&, SourceSlot& out)>;

// Reducer: STATE-CARRYING fold over the per-source values for one (atom,frame).
// Produces the aggregated row's payload alongside the un-summed rows (the two
// row kinds). The classifier is captured so per-key sums work (sumByKey).
struct ReducedRow {
    int                 n_sources = 0;
    int                 n_sources_valid = 0;
    double              sum_dipolar_all = 0.0;        // l=2 dipolar pool (ring/mc)
    double              sum_dipolar_valid = 0.0;
    std::vector<double> per_key_sums;                 // per ring-type / bond-category / ...
    Vec3                pooled_vec  = Vec3::Zero();    // l=1 pool (charge dipole, E-field)
    Mat3                pooled_t2   = Mat3::Zero();    // l=2 tensor pool (quadrupole, EFG)
    double              cutoff_A    = std::numeric_limits<double>::quiet_NaN();  // recorded
};
using Reducer = std::function<ReducedRow(const Body&, std::size_t atom, std::size_t frame,
                                         const std::vector<SourceSlot>&,
                                         const SourceClassifier&)>;

// TargetFn: the shared DFT target (raw 3x3 + library-basis T0/T1/T2 + local-frame
// total). Generalizes ExtractionSupport::BuildTarget (verified).
using TargetFn = std::function<DftTarget(const Body&, std::size_t atom,
                                         std::size_t orig_frame, const LocalFrame&)>;
```

### a.5 — primitive-verb closure FACTORIES (currying)

The composition mechanism (`SURFACE_DESIGN.md` "Iterated closures"). Each
factory captures the config; the returned closure captures the config and,
at call time, the body. `near(cloud,cut)` returns the configured
`(atom,frame)->[SourceRef]`. This is the heart — losing it returns us to the
one-off's hand-inlined loop.

```cpp
// ── source selectors (two backends, unified) ──
SourceSelector near (CloudKind cloud, double cutoff);          // KD-tree radius query
SourceSelector slots(ArrayId ring_nbhd /* = RingNeighbourhood */); // precomputed-slot backend
SourceSelector self ();                                        // degenerate: returns [self] (ISSUE E1)

// ── local frames ──
LocalFrameFn aromaticHFrame();   // ringsOf -> ringGeom -> canon-normal -> atomOf(res,Γ|Δ) -> BuildAromaticHFrame
LocalFrameFn hnFrame();          // N/H/CA/Cprev -> BuildHNFrame
LocalFrameFn labFrame();         // identity (no per-atom frame; T0/|T2| only)

// ── attachers (geometry / identity / kernel) ──
Attacher dispIn      (const LocalFrameFn&);   // disp_local = toLocal(source_pt - pos)
Attacher dipolarForm ();                      // r, cosθ, (3cos²θ−1)/r³ about a source axis
Attacher ringNormalIn(const LocalFrameFn&);   // source_normal_local (the l=2 dipole axis)
Attacher bondAxisIn  (const LocalFrameFn&);   // bond_axis_local (unit B−A in frame)
Attacher ringIdentity();                      // QtRing virtuals (intensity, N count, …)
Attacher bondIdentity();                      // QtBond category / order / endpoint elems
Attacher chargeWeight(ArrayId charge_source); // q_i for charge multipoles (ISSUE C-family)
Attacher kernelOf    (ArrayId kernel);        // value(kernel,atom,frame) — cross-check

// ── classifiers ──
SourceClassifier selfBonded(/* uses idx.ownRings/ownAtoms */);  // is_self_or_bonded
SourceClassifier ringType();                                    // key = ring_type_index
SourceClassifier bondCategory();                                // key = bond_category
SourceClassifier always();                                      // valid, key = 0

// classifier composition: ⊕ ANDs validity, picks the right key.
SourceClassifier operator+(SourceClassifier, SourceClassifier);

// ── reducers (state-carrying folds) ──
Reducer sumDipolarByKey(const SourceClassifier&);   // ring/mc: per-key l=2 dipolar sums
Reducer poolVecByKey   (const SourceClassifier&);   // l=1: Σ q_i (r_i − r_atom)  (dipole / E-field)
Reducer poolT2ByKey    (const SourceClassifier&);   // l=2: Σ traceless 3rr−r²I    (quadrupole / EFG)
Reducer deltaOver      (const WindowFn&, Reducer);  // Δ over the window (rate-of-change)
Reducer passThrough    ();                          // un-summed only (rows already emitted)

// ── target ──
TargetFn dftTarget();   // BuildTarget: raw 3×3 + DecomposeLibrary(T0/T1/T2) + toLocal
```

### a.6 — the `Relationship` bundle (named, defined in code)

`SURFACE_DESIGN.md` Layer 2. A relationship is a **named bundle** — direct,
named closures composed in code. NOT discovered, NOT registry-loaded, NOT a
plugin. The nine named bundles in Part (c) each construct one of these.

```cpp
struct Relationship {
    QString                       name;        // "ring_current", "efg", …
    Stratum                       stratum;
    LocalFrameFn                  frame_fn;
    std::vector<SourceSelector>   selectors;   // heterogeneous sources = more selectors
    std::vector<Attacher>         attachers;
    SourceClassifier              classifier;
    std::optional<WindowFn>       window_fn;    // Δ / rate-of-change (optional)
    Reducer                       reducer;
    TargetFn                      target_fn;
    ArrayId                       cross_check_kernel = ArrayId::KernelBs;  // bare-kernel
    QString                       l_level;      // "l2" | "l1" | "l1+l0" | "feature" — target path
};

// The named bundles, defined in code (NOT a registry). One free function per
// relationship, returning the fully-curried bundle. main() selects by name.
namespace bundles {
    Relationship ring_current();
    Relationship mcconnell();
    Relationship buckingham_efield();
    Relationship efg();
    Relationship charge_dipole(ChargeSource);       // charge_source is a PARAMETER
    Relationship charge_quadrupole(ChargeSource);   // same parameter
    Relationship larsen_hbond();
    Relationship charge_response_gradient();
    Relationship aimnet2_embedding();
}

// ChargeSource — REQUIRED parameter, NO default (substrate_conventions).
enum class ChargeSource : int { FF14SB, AIMNet2, MOPAC };
```

### a.7 — the engine (one pure, order-free loop) + sink

`SURFACE_DESIGN.md` Layer 3. Map over the (atom,frame) index space, inner
map over sources, state-carrying fold. Order is free (referential
transparency) — schedule is a cache choice. Writes the two row kinds through
the existing `RecordSink` (verified: `WriteSourceRows` + `WriteAggregatedRow`,
`QSaveFile`, two CSVs per case).

```cpp
// Build the schema for a relationship's two row kinds (generalizes the per-
// extraction schema() in the one-off; column set is the union the bundle fills).
FeatureSchema SchemaFor(const Relationship&);

// THE ENGINE. Pure transforms over the immutable Body; the fold is the only
// state. Iterated closures: outer map over index_space, inner map over sources.
//
//   for (atom, frame) in index_space:                # order FREE → choose for cache
//       lf   = frame_fn(body, atom, frame)
//       src  = flatten( sel(body, atom, frame) for sel in selectors )   # cutoff captured
//       slots= [ fill(attachers, body, atom, frame, s, lf) for s in src ]
//       row  = reducer(body, atom, frame, slots, classifier)
//       sink.WriteSourceRows(rec(slots, lf, target_fn(...)))            # un-summed
//       sink.WriteAggregatedRow(rec, row)                              # aggregated
std::size_t RunEngine(const Body&, const Relationship&, RecordSink&);

// Index space = (stratum atom, dft frame row). dftRows() from FrameMap (verified).
struct IndexSpace { std::vector<std::size_t> atoms; std::vector<std::size_t> dft_rows; };
IndexSpace BuildIndexSpace(const Body&, const Relationship&);
```

### a.8 — the batch-file `main`

`SURFACE_DESIGN.md` "Python interface (RESOLVED — it's a batch file)".
`QCoreApplication`, `StructuredLogger::Install()` first, `QCommandLineParser`,
synchronous load → traverse → commit → rc. Generalizes the one-off
`main_extract.cpp` (verified) from two hard-coded extractions to a
JSON-scenario-driven named-bundle selection.

```cpp
// h5reader-rediscover --scenario spec.json --out <dir>   > <dir>/run.log 2>&1
//   (or the legacy --run <calcset> --case ring|mc|all path, kept working).
// rc classes (kept from codex, SURFACE_DESIGN): 0 ok, 2 bad-args, 1 load,
// 3 sink-open, 4 commit. (130 cancel reserved for the future control plane.)
int main(int argc, char** argv);   // declared; body is the stub's a.8
```

---

## Part (b) — the JSON scenario-spec + manifest

`SURFACE_DESIGN.md`: IN = a JSON scenario spec selecting a *named*
relationship + params; OUT = atomic CSV substrate + `manifest.json`.
Declarative format ⇒ flagged for discuss-before-finalize
(`feedback_data_format_discuss_first`) — see ISSUE D1. This is the stub
shape, not a finalized schema.

### b.1 — scenario spec (IN)

```json
{
  "scenario_schema_version": "0.1.0",
  "relationship": "ring_current",
  "run": "/shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft",
  "out": "/tmp/rediscover-out-ring",
  "params": {
    "cutoff_A": 8.0,
    "charge_source": null,
    "window_w": 0,
    "strata_filter": null
  },
  "emit": { "sources_csv": true, "aggregated_csv": true }
}
```

- `relationship` is the **named bundle** key (a.6) — selected, never loaded.
- `cutoff_A` REQUIRED for any relationship whose selector is `near()`; for
  `slots()`-backed relationships (ring_current) it is the recorded H5
  `ring_current_spatial_cutoff_A` provenance, not a query cutoff (ISSUE I2).
- `charge_source` REQUIRED (`"ff14sb" | "aimnet2" | "mopac"`) for the charge
  multipole relationships, `null` otherwise. No default (substrate_conventions).
- `window_w` 0 = instantaneous (the only validated mode); >0 selects the Δ
  reducer (a.5 `deltaOver`).

### b.2 — manifest (OUT, "the end bit")

Written via `QSaveFile` + `QJsonDocument::toJson(Indented)`. Atomic-on-commit
⇒ output-exists-iff-success. Records every convention in effect (the
provenance sidecar from substrate_conventions). Success = exit code + manifest
present.

```json
{
  "substrate_version": "0.1.0",
  "run_id": "01234567-89ab-cdef-0123-456789abcdef",
  "relationship": "ring_current",
  "rc": 0,
  "outputs": {
    "sources_csv": "ring_current_sources.csv",
    "aggregated_csv": "ring_current_aggregated.csv"
  },
  "counts": { "cases": 20500, "source_rows": 141000, "aggregated_rows": 20500 },
  "params": { "cutoff_A": 8.0, "charge_source": null, "window_w": 0 },
  "dft_frame_alignment": {
    "n_frames": 500, "mean_angle_deg": 8.9e-5, "max_angle_deg": 2.4e-4,
    "t2_components": "FRAME-ALIGNED (comparable as emitted)"
  },
  "conventions": {
    "spherical_harmonics": { "basis_ordering_ref": "src/Types.cpp::SphericalTensor::Decompose",
                             "normalization": "isometric_real_sh", "l_max": 2 },
    "dipolar_form": { "convention": "3cos2_minus_1", "reference": "GEOMETRIC_KERNEL_CATALOGUE.md" },
    "local_frames": { "ring_current": "ring_normal_x_anchored_CG_or_CD2",
                      "mcconnell": "HN_amide_plane_via_C_prev" },
    "charge_source": null,
    "multipoles": { "quadrupole": "traceless", "origin": "target_atom" },
    "cutoffs": { "policy": "required_no_defaults" },
    "charge_units": "elementary_charge_e", "distance_units": "angstrom",
    "angle_units_internal": "radians"
  },
  "build_metadata": { "git_commit": "<hash>", "build_preset": "linux-gcc" }
}
```

---

## Part (c) — a worked example for EVERY one of the nine items

Each is a bundle of curried closures (stratum, frame_fn, selector(s),
attacher(s), classifier, reducer, target). The DFT target (raw 3×3 +
library-basis T0/T1/T2 + local-frame) is the **shared target** for all nine —
`dftTarget()` everywhere. `charge_source` is a parameter (items 5/6), not a
separate item.

Verified-vs-inferred: items 1 and 2 are the one-off, re-expressed as curried
bundles (verified against `RingCurrentNeighborhood.cpp` /
`McConnellNeighborhood.cpp`). Items 3–9 are new bundles whose source arrays I
verified exist in `QtTrajectoryH5.h` (3,4,8,9) or established as gaps
(5,6,7 — see issues); the closure composition is inferred from the language.

### 1. ring_current — aromatic C–H; sources = rings via slots(); T0+T2; cross-check = bs

```cpp
Relationship bundles::ring_current() {
    return Relationship{
      .name      = "ring_current",
      .stratum   = atomsWhere([](const model::QtAtom& a){ return a.IsAromaticRingHydrogen(); }),
      .frame_fn  = aromaticHFrame(),                 // canon-normal + typed CG/CD2 anchor
      .selectors = { slots(ArrayId::RingNeighbourhood) },  // frame-0 membership snapshot
      .attachers = { dispIn(aromaticHFrame()),       // disp_local (ring center − H)
                     dipolarForm(),                  // r, cosθ=z/r, (3cos²θ−1)/r³
                     ringNormalIn(aromaticHFrame()), // source_normal_local (l=2 axis)
                     ringIdentity(),                 // QtRing virtuals (intensity, N, fused…)
                     kernelOf(ArrayId::KernelBs) },  // bare-kernel cross-check
      .classifier= selfBonded() + ringType(),        // drop self/fused; key=ring_type_index
      .window_fn = std::nullopt,
      .reducer   = sumDipolarByKey(selfBonded() + ringType()),  // sum_all vs sum_valid + per-type
      .target_fn = dftTarget(),
      .cross_check_kernel = ArrayId::KernelBs,
      .l_level   = "l2",
    };
}
```

Oracle (FINDINGS.md): leave-atoms-out k≈21 ppm·Å³, R²=0.62 coupled; T2 R²=0.44;
|T2| r=0.75; basis 4.9e-8. This bundle must reproduce it (the validation gate).

### 2. mcconnell — backbone HN; sources = bonds via near(BondMidpoints,cutoff); T0+T2; cross-check = mc

```cpp
Relationship bundles::mcconnell() {
    return Relationship{
      .name      = "mcconnell",
      .stratum   = atomsWhere([](const model::QtAtom& a){ return a.IsBackboneAmideHydrogen(); }),
      .frame_fn  = hnFrame(),
      .selectors = { near(CloudKind::BondMidpoints, /*cutoff_A=*/8.0) },  // captured + recorded
      .attachers = { dispIn(hnFrame()),              // disp_local (midpoint − HN)
                     dipolarForm(),                  // cosθ about BOND AXIS, (3cos²θ−1)/r³
                     bondAxisIn(hnFrame()),          // bond_axis_local (unit B−A)
                     bondIdentity(),                 // category / order / endpoint elems
                     kernelOf(ArrayId::KernelMc) },
      .classifier= always() + bondCategory(),        // no self-ring concept; key=category
      .window_fn = std::nullopt,
      .reducer   = sumDipolarByKey(bondCategory()),   // per-category sums; valid==all
      .target_fn = dftTarget(),
      .cross_check_kernel = ArrayId::KernelMc,
      .l_level   = "l2",
    };
}
```

Oracle (FINDINGS.md): scalar form recovers (h vs (3cos²θ−1)/r³ r=0.92); producer-
kernel reconstruction R²≈0.55 (the known gap — likely fuller anisotropy; ISSUE M1).

### 3. buckingham_efield — atoms; source = APBS E-field (l=1 vector); T0/T1

The E-field is a **per-atom** l=1 vector already at the target atom
(`apbs_efield_time_series`, verified `QtVec3TimeSeries`). It is not a
source-cloud sum — it is a degenerate selector returning [self], the field
read at the atom. The l=1 target path is T1 (see ISSUE L1 — does the substrate
need an l=1 target decomposition path?).

```cpp
Relationship bundles::buckingham_efield() {
    return Relationship{
      .name      = "buckingham_efield",
      .stratum   = atomsWhere([](const model::QtAtom&){ return true; }),  // all atoms
      .frame_fn  = labFrame(),                       // E-field is already in the lab/MD frame
      .selectors = { self() },                       // degenerate: [self] (ISSUE E1)
      .attachers = { /* read the field as the source value */
                     [](const Body& b, std::size_t a, std::size_t f, const SourceRef&,
                        const LocalFrame&, SourceSlot& s){
                          const Vec3 E = b.catalog.valueVec3(b, ArrayId::ApbsEfield, a, f);
                          s.disp_local = E;           // l=1 vector source value
                     } },
      .classifier= always(),
      .window_fn = std::nullopt,
      .reducer   = poolVecByKey(always()),            // trivial pool of one l=1 vector
      .target_fn = dftTarget(),                       // T0 + T1 the relevant levels
      .cross_check_kernel = ArrayId::KernelBs,        // no bare Buckingham kernel TS — ISSUE K1
      .l_level   = "l1+l0",
    };
}
```

### 4. efg — atoms; source = APBS EFG (T2); T2

The EFG is a **per-atom T2** already at the target atom
(`apbs_efg_time_series/t2`, verified `QtT2TimeSeries`, 5 components, library
m=−2..+2 order per its `irrep_layout`). Like item 3, a self-selector reading
a per-atom field, but the level is l=2. This is `PLAN.md`'s named "natural
next relationship."

```cpp
Relationship bundles::efg() {
    return Relationship{
      .name      = "efg",
      .stratum   = atomsWhere([](const model::QtAtom&){ return true; }),
      .frame_fn  = labFrame(),                       // EFG T2 emitted in the H5/MD frame
      .selectors = { self() },                       // degenerate [self]
      .attachers = { [](const Body& b, std::size_t a, std::size_t f, const SourceRef&,
                        const LocalFrame&, SourceSlot& s){
                          const std::array<double,5> t2 =
                              b.catalog.valueT2(b, ArrayId::ApbsEfg, a, f);  // source T2 value
                          // store the 5-vec on the slot (T2 carried, never collapsed)
                          for (int k=0;k<5;++k) s.disp_local[k%3] = t2[k];   // see ISSUE T1 (slot shape)
                     } },
      .classifier= always(),
      .window_fn = std::nullopt,
      .reducer   = poolT2ByKey(always()),             // single source T2 → pooled_t2
      .target_fn = dftTarget(),                       // T2 the relevant level
      .cross_check_kernel = ArrayId::ApbsEfg,         // the EFG IS its own cross-check ref
      .l_level   = "l2",
    };
}
```

### 5. charge_dipole — atoms; local Σ q_i (r_i − r_atom); charge_source PARAMETER; T0/T1

Genuine source-cloud sum: charge sites within a cutoff, weighted by `q_i`,
displacement from the target. `charge_source` is the knob (one bundle, a
parameter — per the conventions doc). Multipole conventions from
substrate_conventions: origin at target atom, target excluded, traceless not
relevant for the dipole, cutoff required.

```cpp
Relationship bundles::charge_dipole(ChargeSource cs) {
    const ArrayId q = (cs == ChargeSource::FF14SB)  ? ArrayId::Ff14sbCharge
                    : (cs == ChargeSource::AIMNet2) ? ArrayId::Aimnet2Charge
                                                    : ArrayId::MopacCharge;   // C1/C2/C3
    return Relationship{
      .name      = QStringLiteral("charge_dipole_%1").arg(chargeSourceTag(cs)),
      .stratum   = atomsWhere([](const model::QtAtom&){ return true; }),
      .frame_fn  = labFrame(),                        // μ reported atom-local OR lab (ISSUE F1)
      .selectors = { near(CloudKind::ChargeSites, /*cutoff_A=*/6.0) },  // required, recorded
      .attachers = { dispIn(labFrame()),              // (r_i − r_atom)
                     chargeWeight(q) },               // q_i for this source
      .classifier= excludeSelf(),                     // target atom excluded from the sum
      .window_fn = std::nullopt,
      .reducer   = poolVecByKey(excludeSelf()),       // μ = Σ q_i (r_i − r_atom)  → pooled_vec
      .target_fn = dftTarget(),                       // T0 + T1
      .cross_check_kernel = ArrayId::KernelBs,         // no bare charge-dipole kernel TS — ISSUE K1
      .l_level   = "l1+l0",
    };
}
```

### 6. charge_quadrupole — atoms; traceless local Σ q_i (3 r r − r² I); same charge_source; T2

Same charge-site source cloud + same `charge_source` parameter as item 5; the
reducer pools the traceless rank-2 form instead of the vector. Fixture from
substrate_conventions: +1e at r=1Å along z ⇒ Q_zz=+2, Q_xx=Q_yy=−1.

```cpp
Relationship bundles::charge_quadrupole(ChargeSource cs) {
    const ArrayId q = chargeArrayFor(cs);             // FF14SB | AIMNet2 | MOPAC
    return Relationship{
      .name      = QStringLiteral("charge_quadrupole_%1").arg(chargeSourceTag(cs)),
      .stratum   = atomsWhere([](const model::QtAtom&){ return true; }),
      .frame_fn  = labFrame(),
      .selectors = { near(CloudKind::ChargeSites, /*cutoff_A=*/6.0) },
      .attachers = { dispIn(labFrame()), chargeWeight(q) },
      .classifier= excludeSelf(),
      .window_fn = std::nullopt,
      .reducer   = poolT2ByKey(excludeSelf()),         // Q_ab = Σ q_i [3 d_a d_b − δ_ab |d|²], traceless
      .target_fn = dftTarget(),                        // T2
      .cross_check_kernel = ArrayId::ApbsEfg,          // EFG is the field-gradient cross-check
      .l_level   = "l2",
    };
}
```

### 7. larsen_hbond — exchangeable H; donor/acceptor geometry; T0+T2

Stratum = exchangeable H (amide HN, hydroxyl, etc.). `QtAtom::isExchangeable`
exists (verified). Sources = H-bond acceptors found geometrically (the Larsen
spatial criterion). The H5 carries Larsen H-bond shielding TS
(`larsenHBond1pHBShielding` etc., verified) as the cross-check, and a per-atom
`larsenHBondCount` scalar. The detection method is a REQUIRED parameter
(substrate_conventions `HBondDetection`) — see ISSUE H1 (stratum + detection).

```cpp
Relationship bundles::larsen_hbond() {
    return Relationship{
      .name      = "larsen_hbond",
      .stratum   = atomsWhere([](const model::QtAtom& a){ return a.isExchangeable; }),  // ISSUE H1
      .frame_fn  = hnFrame(),    // amide HN frame for HN donors; hydroxyl needs its own (ISSUE H2)
      .selectors = { near(CloudKind::Atoms, /*cutoff_A=*/3.5) },   // Larsen primary H...A
      .attachers = { dispIn(hnFrame()),               // donor→acceptor geometry, local
                     dipolarForm(),                   // r(H...A), angle(D-H...A)
                     // acceptor identity + Larsen class (ClassifyAcceptor): NEW attacher (ISSUE H3)
                     [](const Body& b, std::size_t a, std::size_t f, const SourceRef& s,
                        const LocalFrame&, SourceSlot& out){ /* classify acceptor type */ },
                     kernelOf(ArrayId::KernelBs) },    // see ISSUE K2 (which larsen TS is the kernel)
      .classifier= hbondAcceptorValid(),               // only true acceptors (N/O lone pairs)
      .window_fn = std::nullopt,
      .reducer   = sumDipolarByKey(hbondAcceptorValid()),  // T0 + T2 over acceptors
      .target_fn = dftTarget(),
      .cross_check_kernel = ArrayId::KernelBs,          // placeholder — ISSUE K2
      .l_level   = "l2",
    };
}
```

### 8. charge_response_gradient — atoms; source = AIMNet2 CRG (scalar+vector); T0+T2

The AIMNet2 charge-response gradient is **per-atom** scalar + vector
(`aimnet2ChargeResponseGradient()`, verified
`QtAimnet2ChargeResponseGradientTimeSeries`: `.scalar` (N,T), `.vec` (N,T,3)).
A self-selector reading both per-atom channels. This is CRG,
`d(sum(q^2))/dr`, not a polarizability or true alpha. It is an l=0
(scalar) + l=1 (vector) feature; against the DFT T2 it is the modulating
feature, not a source sum.

```cpp
Relationship bundles::charge_response_gradient() {
    return Relationship{
      .name      = "charge_response_gradient",
      .stratum   = atomsWhere([](const model::QtAtom&){ return true; }),
      .frame_fn  = labFrame(),
      .selectors = { self() },
      .attachers = { [](const Body& b, std::size_t a, std::size_t f, const SourceRef&,
                        const LocalFrame&, SourceSlot& s){
                          s.r          = b.catalog.value(b, ArrayId::Aimnet2ChargeRespScalar, a, f);
                          s.disp_local = b.catalog.valueVec3(b, ArrayId::Aimnet2ChargeRespVector, a, f);
                     } },
      .classifier= always(),
      .window_fn = std::nullopt,
      .reducer   = passThrough(),                      // one self-source; scalar+vec carried
      .target_fn = dftTarget(),                        // T0 + T2 (scalar response → T0; vec → T1)
      .cross_check_kernel = ArrayId::KernelBs,          // no dedicated kernel TS — ISSUE K1
      .l_level   = "l1+l0",
    };
}
```

### 9. aimnet2_embedding — atoms; per-atom 256-d learned embedding as feature → DFT σ

The 256-d per-atom AIMNet2 embedding (`aimnet2Embedding()`, verified
`QtEmbeddingTimeSeries`, float32). This is a **per-atom feature**, NOT a
source-neighborhood sum. It composes as a degenerate `self()` selector
returning [self] and an attacher that reads the 256-d span. This is awkward in
a source-sum language — surfaced as ISSUE E1/E2.

```cpp
Relationship bundles::aimnet2_embedding() {
    return Relationship{
      .name      = "aimnet2_embedding",
      .stratum   = atomsWhere([](const model::QtAtom&){ return true; }),
      .frame_fn  = labFrame(),
      .selectors = { self() },                         // degenerate: [self] (ISSUE E1)
      .attachers = { [](const Body& b, std::size_t a, std::size_t f, const SourceRef&,
                        const LocalFrame&, SourceSlot& s){
                          std::size_t n = 0;
                          const float* emb = b.catalog.valueEmbedding(b, ArrayId::Aimnet2Embedding,
                                                                       a, f, n);  // n == 256
                          // The 256-d vector does NOT fit SourceSlot's fixed fields — ISSUE E2.
                     } },
      .classifier= always(),
      .window_fn = std::nullopt,
      .reducer   = passThrough(),                      // feature row, no sum
      .target_fn = dftTarget(),                        // σ (T0 + T2) the target
      .cross_check_kernel = ArrayId::KernelBs,          // no kernel cross-check for a learned feature
      .l_level   = "feature",
    };
}
```

---

## Part (d) — SURFACED ISSUES (findings to fix, not deferrals)

Each is a finding. None is "out of scope" or "later"; each names what to do.

### Missing catalog entries

- **C1 — FF14SB charges are NOT in the reader model.** Verified: the H5 has
  no FF14SB charge dataset, and `QtAtom`/`QtTopology` carry only `formalCharge`
  (int8), not partial charges. The conventions doc says
  `traj.charge_ff14sb(atom)` reads "from prmtop." `RunData` loads from the
  `.LGS` calcset; whether the prmtop/`topol.top` partial charges are reachable
  through `QtProteinLoader` is unconfirmed (I read the header, not the .cpp).
  FIX: add a `Ff14sbCharge` catalog entry with `residence = StaticPrmtop`,
  backed by a new per-atom static charge load at `RunLoader::Load` time. If the
  prmtop is not in the calcset, this entry's `present()` returns false and
  charge_dipole(FF14SB)/charge_quadrupole(FF14SB) must fail-loud at validate,
  not emit zeros.
- **C2 — MOPAC charges are NOT present in this calcset.** Verified: no MOPAC
  charge per-atom TS on `QtTrajectoryH5` (only `mopacChargeWelford` rollup,
  `mopac_*_shielding`). The prompt says stub it anyway. FIX: `MopacCharge`
  catalog entry with `residence = Absent` until the per-frame MOPAC charge
  dataset lands; `Catalog::has(MopacCharge)` returns false; the charge_dipole/
  quadrupole(MOPAC) bundles validate-fail with a clear message rather than
  producing a silent zero field. Do NOT fabricate a fallback chain
  (substrate_conventions explicitly forbids `charge()` auto-pick).
- **C3 — "charge sites" cloud is undefined for each ChargeSource.** The
  `CloudKind::ChargeSites` tree (a.2) needs a membership rule: which atoms are
  charge sites, and with what charge. For FF14SB every atom carries a partial
  charge (so the cloud is all atoms). For AIMNet2 likewise (per-atom
  `aimnet2Charge`). For MOPAC, absent. FIX: `ResidentIndexes::chargeSites()`
  must be parameterised by `ChargeSource` (the stub currently shows it
  unparameterised — that is the bug). The KD-tree over charge sites is then
  per-(source, frame); for FF14SB/AIMNet2 it coincides with the atoms tree, so
  it can alias `CloudKind::Atoms` + a per-source weight rather than a 4th tree.

### Missing day-1 indexes

- **I1 — the typed-atom `(residue, locant[, branch, element]) -> atom`
  index does not exist.** Verified: the one-off does an O(ring) `Locant` scan
  inside `typedRingAnchor` (RingCurrentNeighborhood.cpp) and uses
  `QtResidue`'s N/CA/C cache for HN. SURFACE_DESIGN names this the day-1
  `atomOf` index that generalizes that cache to ALL locants — the
  identity-from-chemistry index. FIX: build `ResidentIndexes::atomOf` at load
  (a.2). Without it, `aromaticHFrame()` and `hnFrame()` re-derive anchors per
  call. This is also the guard against the IUPAC-positional trap
  (`feedback_identity_from_chemistry_not_position`): the index is keyed on
  typed `Locant`, never on name string or ring-walk position.
- **I2 — per-cloud KD-trees are lazy in the one-off, day-1 in the design.**
  Verified: `FrameSpatialIndex` is built per-frame inside the McConnell loop
  (`McConnellNeighborhood.cpp:139`, one tree per frame, bond-midpoints only).
  SURFACE_DESIGN requires all four clouds × all frames built up front so
  scoping never pays build cost mid-traversal. FIX: `BuildResidentIndexes`
  builds `trees[cloud][frame]` for {Atoms, BondMidpoints, RingCenters,
  ChargeSites} at load. ~751 frames × 4 clouds × ~850 pts is tens of KB/tree —
  cheap. RingCenters is new (the one-off reads ring geometry on demand via
  `RingGeometryAt`); it needs the per-frame ring-center point list.
- **I3 — the `slots()` selector's "query cutoff" is a category error.**
  The ring-current selector reads the H5 `ring_neighbourhood` slots (frozen
  frame-0 membership, `ring_current_spatial_cutoff_A = 15.0` recorded on the
  buffer, verified). It is NOT a runtime radius query — there is no cutoff to
  pass. The scenario spec's `cutoff_A` for ring_current is provenance (the
  producer's 15 Å), not a knob. FIX: `slots()` records the buffer's
  `ring_current_spatial_cutoff_A` into the manifest and the spec's `cutoff_A`
  for ring_current is validated == that value (or null), never used as a query.

### Undecided conventions

- **D1 — the JSON scenario-spec schema is a declarative format ⇒ discuss
  before finalizing.** `feedback_data_format_discuss_first` + SURFACE_DESIGN's
  own note. Part (b) is a stub shape, not a settled schema. Open questions:
  flat `params` vs per-relationship typed param blocks; how `charge_source`
  null-vs-required is encoded; whether one spec can request multiple
  relationships (a batch within a batch) or strictly one. Surface to Jessica
  before writing the parser.
- **F1 — charge-multipole reporting frame.** substrate_conventions says
  dipole/quadrupole are "always atom-local," origin at the target atom. The
  stub bundles for items 5/6 use `labFrame()` for the displacement but the
  multipole itself is origin-at-atom; whether the pooled vector/tensor is then
  rotated into a per-atom local frame (for equivariance) or left lab is
  undecided. For charged-residue targets there is no canonical aromatic/HN
  frame. FIX: decide per-atom-class frame vs lab for charge multipoles; record
  in the manifest `local_frames`.
- **F2 — local frame for "all atoms" strata (items 3,4,5,6,8,9).** The
  one-off only has `aromaticHFrame` and `hnFrame` (verified `FrameVariant`:
  Invalid, HN_Standard, HN_NTerminus, AromaticHRing). Items 3–9 run on all
  atoms or exchangeable H, which have no aromatic/HN frame. The stub uses
  `labFrame()` (identity), which means their T2 components are in the MD/lab
  frame, not a rotation-stable per-atom frame. T0 and |T2| are safe; T2
  *components* for the all-atom relationships are lab-frame. FIX: either accept
  lab-frame T2 components for these (and flag in the manifest, as the existing
  T2-frame-alignment check already does) or build the conventions doc's Cα /
  CO / HA frames (`LocalFrameBasis` has only HN + aromatic-H today). This is the
  same class as STATE.md caveat 2, generalized.

### Items the language can't cleanly express

- **E1 — the degenerate `self()` selector is a smell.** Items 3,4,8,9 are
  per-atom fields/features, not source sums. Modeling them as a selector that
  returns `[self]` works mechanically but inverts the language's intent (a
  selector is "scope outward to sources"). It is honest but awkward. FIX
  (small): keep `self()` as a named, documented degenerate selector — it
  composes cleanly and the engine needs no special case (which is the point).
  Note in the manifest that these relationships are `selector: self` so a
  downstream reader knows the "source set" is the atom itself. The awkwardness
  is acceptable; the alternative (a separate per-atom-feature code path) breaks
  "one engine, one loop."
- **E2 — `SourceSlot` cannot hold a 256-d embedding (item 9).** Verified:
  `SourceSlot` (RediscoverTypes.h) is a fixed struct of scalars + a few Vec3.
  The 256-d float32 embedding does not fit. The stub's attacher reads the span
  but has nowhere to put it. FIX: the embedding relationship needs either (a) a
  variable-width column block in the CSV schema (256 named `emb_0..emb_255`
  columns appended to the source row — straightforward, just wide), or (b) a
  side Parquet/NPY emitted alongside the CSV keyed by (atom,frame). The CSV-wide
  option is simplest and keeps the one-sink discipline; flag the row width
  (~280 columns) to Jessica. This is the embedding's surfaced awkwardness the
  prompt asked for: it is a feature, not a source, and the source-slot carrier
  was built for sources.
- **T1 — `SourceSlot` has no natural slot for a 5-component T2 source value
  (item 4 EFG).** Verified: `SourceSlot.disp_local` is a Vec3; the EFG source
  value is 5 components. The stub crammed it into `disp_local` indices, which is
  wrong. FIX: add a `std::array<double,5> source_t2` field to `SourceSlot` (or
  reuse `pooled_t2` on the reduced row and skip per-source for self-selector
  relationships). Since items 4/6 are the only T2-source bundles and 4 is a
  self-selector, the cleanest fix is: self-selector T2 relationships write the
  T2 only on the aggregated row (`ReducedRow.pooled_t2`) and emit no per-source
  T2 — the per-source row is degenerate anyway (one self source).
- **L1 — l=1 target decomposition path (items 3,5,8).** The DFT target is a
  rank-2 shielding tensor decomposed to T0/T1/T2. The Buckingham E-field and the
  charge dipole are l=1 *sources/features*; the question is whether the shielding
  T1 (the antisymmetric pseudovector part, present in `SphericalTensor.T1`) is
  the right target slice for an l=1 relationship, or whether l=1 sources only
  ever modulate T0 (isotropic) and T2 (the shielding has no physical l=1 the way
  the field does). Verified ambiguity: Types.h flags TWO T1 conventions
  (Cartesian pseudovector on the NPY path vs m-basis on the H5 path) and a
  parity inconsistency (catalog 1e vs SDK 1o). FIX: decide whether l=1
  relationships target T0 only, or T0 + the shielding-T1, and resolve the T1
  basis/parity question first (it is a latent bug for any T1 comparison). Until
  resolved, l=1 bundles should target T0 + |T2| invariants and flag T1
  unverified.
- **K1 — several relationships have no bare-kernel cross-check TS.** Verified
  H5 inventory: there are bare-kernel time-series for ring (`bs`), McConnell
  (`mc`), and the EFG (`apbs_efg`) — items 1, 2, 4 have a real cross-check. But
  Buckingham E-field (3), charge dipole (5), charge quadrupole (6),
  charge-response-gradient (8), and the embedding (9) have NO producer
  bare-kernel TS to cross-check against. The stub sets `cross_check_kernel =
  KernelBs` as a placeholder, which is wrong (it would compare an E-field
  relationship to the ring-current kernel). FIX: make `cross_check_kernel`
  optional (`std::optional<ArrayId>`); for relationships without a producer
  kernel, emit no bare-kernel columns and rely on the DFT target alone. The
  one-off always had a kernel (bs/mc); the generalized language must not assume
  one exists.
- **K2 — which Larsen TS is the H-bond cross-check (item 7)?** Verified: the
  H5 has four Larsen H-bond shielding TS (`larsen_1pHB`, `larsen_1pHaB`,
  `larsen_2pHB`, `larsen_2pHaB` — primary/secondary × HB/HaB). It is not obvious
  which (or which combination) is the bare kernel for the larsen_hbond
  relationship. FIX: read `src/LarsenHBondShieldingResult.cpp` to map the four
  TS to the donor classes, then set the cross-check per H-bond class (this may
  need a per-class kernel selection, not one `ArrayId`).
- **H1 — the H-bond stratum + detection method are undecided (item 7).**
  `QtAtom::isExchangeable` (verified) gives a candidate stratum, but
  substrate_conventions requires an explicit `HBondDetection` parameter
  (GeometricLarsen vs DSSP_KabschSander) with no default, and the two produce
  different H-bond sets. The stub hard-codes a geometric cutoff (3.5 Å) and the
  `isExchangeable` stratum. FIX: (a) confirm `isExchangeable` covers amide HN +
  hydroxyl + the intended donors; (b) add `detection` as a required scenario
  param; (c) the acceptor cloud is a KD-tree over N/O acceptors, and the Larsen
  class comes from `ClassifyAcceptor` — needs a new attacher (H3 below).
- **H2 — hydroxyl/other-donor local frames don't exist.** `LocalFrameBasis`
  has HN and aromatic-H only. Hydroxyl OH, thiol SH, and the ARG/LYS terminal
  donors in the H-bond stratum have no typed frame. FIX: add their frames to
  `LocalFrameBasis` (the conventions doc specifies HA / Cα / CO frames but not
  hydroxyl) or run those donors with `labFrame()` + a T2-unverified flag.
- **H3 — no acceptor-classification attacher exists.** The Larsen relationship
  needs `ClassifyAcceptor` (acceptor type → Larsen Primary1/Primary2/Secondary1
  /Secondary2) as an `Attacher`. Verified: nothing in the one-off or the reader
  model exposes this (it lives in the library's
  `LarsenHBondShieldingResult.cpp`, which the reader does not link). FIX: port
  the classification rule into a rediscover-side attacher, or read the Larsen
  class off the H5 if the producer recorded it per H-bond (check the Larsen TS
  side-data).

### Inconsistencies found between SURFACE_DESIGN.md and the one-off code

- **X1 — `near()` cutoff: required at call site vs captured by the factory.**
  SURFACE_DESIGN Layer 1 shows `neighbors(atom,frame,cloud,cutoff)` with cutoff
  at the call site, AND the iterated-closure section shows `near(cloud,cutoff)`
  capturing the cutoff so the closure is cutoff-free. Both appear. The one-off
  passes cutoff at the call site (`index.Within(atomPos, cutoff_A)`, verified).
  These are reconcilable (the primitive takes it; the factory captures it) and
  the stub does exactly that (a.3 primitive + a.5 factory), but it is worth
  noting the design text shows both forms — not a contradiction, a layering the
  stub makes explicit. No fix beyond this note.
- **X2 — `RunData` immutability vs a bolted-on index set.** SURFACE_DESIGN
  says "loaded once, resident, immutable" and lists indexes "built day 1." The
  one-off `RunData` (verified) has no index field; the McConnell tree is built
  lazily in the loop. The stub adds `ResidentIndexes` as a *separate* const
  carrier (a.0 `Body`) rather than mutating `RunData`, preserving immutability.
  This is a design choice the stub makes that SURFACE_DESIGN left implicit —
  flag it so the build doesn't instead grow a mutable cache on `RunData`.
- **X3 — the schema column union vs per-extraction schemas.** The one-off has
  two hand-written `schema()` methods (RingCurrent, McConnell) whose source
  columns are nearly identical supersets (verified: both list the full ring +
  bond column set). SURFACE_DESIGN's `SchemaFor(Relationship)` implies one
  generated schema. For nine relationships the union of all columns (ring + bond
  + charge + EFG-5 + E-field-3 + CRG scalar+vec + 256 embedding) is very wide if
  every row carries every column. FIX: `SchemaFor` should emit only the columns
  the bundle's attachers actually fill (a per-bundle column set), not the union
  — otherwise the embedding's 256 columns leak into the ring-current CSV. The
  one-off's superset approach does not scale to nine heterogeneous bundles.
- **X4 — manifest SH convention stub copied a stale value.** The
  substrate_conventions manifest example has
  `"spherical_harmonics": { "convention": "scipy_real", ... }` while its own
  prose (and `SphericalBasis.h`, verified) says the basis is the LIBRARY order
  (`src/Types.cpp`), isometric, NOT SciPy. The stub manifest in Part (b) uses
  `basis_ordering_ref: src/Types.cpp` and drops the `scipy_real` string. Flag:
  the conventions-doc manifest example is internally inconsistent (it says
  scipy_real in JSON but library-order in prose) — fix the conventions doc's
  example, not just the stub.

---

## Validation gate (carried from SURFACE_DESIGN)

The composed **ring_current** bundle (item 1) must reproduce the one-off
ORACLE: leave-atoms-out k≈21 ppm·Å³, R²=0.62 (coupled); equivariant T2
R²=0.44; |T2| r=0.75; basis check 4.9e-8; DFT frame rotation ~1e-4°. Match ⇒
faithful rebuild of the one-off as a curried bundle. The other eight bundles
have no oracle yet (they are new cells); their gate is that they run, emit
well-formed two-kind CSVs + manifest, and the DFT-frame-alignment check passes
(or flags T2 components lab-frame).
