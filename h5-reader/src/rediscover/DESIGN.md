# Rediscover — classical-shielding study substrate

**Status:** design, for Codex review. Branch `h5-reader-pysr-spike`.

Experimental, one-shot branch — no integration target, the viewer is
untouched, no fleet wiring. Judge on doability, physics, and coherence,
not merge risk.

## What this is

A headless CLI (`h5reader-rediscover`, `QCoreApplication` only, links
`h5reader_core` + vendored nanoflann) that loads one 1P9J calcset into an
all-frames-resident `RunData` and emits, per (target atom, frame), a
**method-agnostic substrate** for studying classical shielding:

- the **local source neighbourhood** — the set of sources (aromatic rings;
  anisotropic bonds) with their relative geometry as **vectors** + derived
  scalars + source identity, expressed in a recorded **local frame**;
- the **full DFT shielding tensor** target — raw 3×3 (total/dia/para) +
  decomposition (library basis) + isotropic T0;
- the producer's **bare classical kernels** as fixed-basis cross-checks.

It does **not** decide the fitter. The downstream method — ridge on the
kernels, scalar symbolic regression, equivariant/tensor SR, or an
equivariant sum-pooling network — is a separate, open choice. The
substrate serves all of them.

## Why method-agnostic

The problem is intrinsically **permutation-invariant** over a variable
source set, **rotation-equivariant**, **tensor-valued** (T2 is the thesis),
and a **local sum**: `σ_i = Σ_sources f(x_source)`. That is the shape an
equivariant sum-pooling model is built for; flat scalar regression is the
constrained special case. So the substrate emits the **un-summed source
set with displacement vectors and the target tensor**, which is exactly
what a sum-pooling/equivariant fitter consumes *and* what ridge / scalar
SR can be reduced from. Baking the sum or the T0-scalar into the
extraction would lock out the methods that fit this problem most
naturally. We keep it open — the equivariant / sum-pooling path is
explicitly NOT foreclosed (the Stage-1 MLP rejection says nothing about
it), and `σ_iso` (T0) is emitted alongside the full tensor so a
scalar-vs-equivariant comparison is itself a result.

## The per-(atom, frame) record (the substrate contract)

- **Identity / frame:** atom idx, residue, atom name, element, stratum
  label, H5 row, `original_index`, `time_ps`.
- **Source set (padded slots, un-summed):** per source —
  - identity/type (ring: type + `LiteratureIntensity` + nitrogen_count +
    jb_offset + aromaticity + ring_size + fused; bond: category + order +
    endpoint elements);
  - **relative displacement vector** to the source (in the atom's local
    frame) — not just scalars;
  - derived scalars `r`, `cosθ`, and for rings `z`, `rho`, `in_plane_angle`.
- **Local frame:** the recorded basis `{z, x, y}` for this atom, so the
  source vectors and the target tensor live in one rotation-stable frame.
- **Target:** DFT σ — **raw 3×3** (total/dia/para) + decomposition
  (library basis) + **σ_iso = T0** (the basis-free scalar). Total is the
  physical target; dia/para are diagnostics.
- **Cross-check:** the producer's per-atom bare kernels (`bs`/`mc`
  T0 + T2) — fixed-basis features for ridge and for checking the form.

**Two row kinds per case** (resolves the summed-law regression): the
per-source neighbourhood rows above (un-summed, for a sum-pooling /
equivariant fitter), plus one aggregated per-(atom,frame) row carrying
`Σ_sources (3cos²θ−1)/r³` + per-ring-type / per-bond-category summed
features against the same DFT target (the well-posed input for a scalar /
ridge fitter).

## Basis — both parts are needed

**(A) T2 component consistency.** The H5 kernels are decomposed in the
library order `[xy, yz, zz, xz, xx−yy]` (`src/Types.cpp:53`); the reader's
`SphericalDecomposition.cpp` uses a different order. Add
`rediscover::SphericalBasis::DecomposeLibrary(Mat3) -> SphericalTensor` in
the library order + isometric normalization, applied to the DFT **raw**
3×3 so DFT-T2 and kernel-T2 share a basis. Fixture (corrected): for
`3ẑẑᵀ−I = diag(−1,−1,2)`, `Szz=2` ⇒ `T2[2] = √(3/2)·2 = √6` (NOT √(3/2),
which is `(3ẑẑᵀ−I)/2`); and assert the dipolar component pattern.
**Caveat:** T2 *component* comparison is valid only if the ORCA tensors are
in the same Cartesian frame as the H5 kernels — the loader does no
rotation check, so record the frame and treat cross-frame T2 comparison as
unverified until checked. T0 (iso) is rotation-invariant and unaffected.

**(B) Per-atom local frame.** Typed frames per atom class (HN amide-plane
frame; aromatic-H ring-normal frame) to express the source displacement
vectors and the target tensor in a consistent, rotation-stable basis —
this is what makes the substrate equivariant-ready (a flat r/cosθ table
throws the vector information an equivariant fitter needs).

## The two extractions (concrete)

- `RediscoveryExtraction` — thin interface: `name()`, `schema() ->
  FeatureSchema`, `extract(const RunData&, RecordSink&)`. Not a framework
  (no scheduler/deps/phases); the driver runs a `vector<unique_ptr<…>>` in
  a plain loop. The `QtRing`-hierarchy idiom, kept thin.
- `RingCurrentNeighborhood` — stratum: aromatic ring-facing H. Reads the
  per-frame ring neighbourhood from the H5 (`ringNeighbourhood` TS:
  distance/rho/z/in_plane_angle; `ring_membership_per_atom` for ring
  identity → `QtRing` physics). `cosθ = z/r`. **Ring membership is the
  frame-0 snapshot** (the H5 freezes it); stated as such, not "every ring
  that ever enters cutoff." Bare `bsShielding` cross-check; DFT target.
- `McConnellNeighborhood` — stratum: backbone amide HN. Discovers
  anisotropic bonds per frame via a `FrameSpatialIndex` (nanoflann over
  bond midpoints). Per bond: **bond axis** (not "C=O axis") and
  `(3cos²θ−1)/r³` form — the `Δχ` lives in the parameter, the kernel is
  `/r³`. Categories: PeptideCO, PeptideCN, sidechain-CO, aromatic, …
  Bare `mcShielding` cross-check; DFT target.

## Class model

- `RunData` — immutable carrier: `unique_ptr<QtProtein>`,
  `unique_ptr<QtTrajectoryH5>`, `DftFrameSet`, `FrameMap`,
  `LocalFrameBasis`. Optional `PairwiseGeometry`.
- `RunLoader` — `Load(calcsetDir) -> RunData`: `QtProteinLoader::Load`,
  `QtTrajectoryH5`, walk `CalcsetManifest.dft.frames[]` via
  `DftShieldingLoader`, build `FrameMap`.
- `FrameMap` — H5 row → `original_index` → DFT frame. **Enforces**
  `frame_index_basis` matches the expected basis and
  `frameIndices().size() == frameCount()`. `dftRows()` = rows with a DFT
  target. (Reuse `DftShieldingStore::metaByOriginal_` mapping.)
- `DftFrameSet` — `vector<DftShieldingFrame>` keyed by original index;
  `target(atom, origIdx) -> optional<DftTarget>`.
- `SphericalBasis` — `DecomposeLibrary` (library T2 order) + the fixture.
- `LocalFrameBasis` — per-atom-class typed frame (HN, aromatic-H ring).
- `FrameSpatialIndex` — per-frame nanoflann KD-tree (bond midpoints),
  built lazily; bond discovery for McConnell.
- `PairwiseGeometry` — **optional / experimental** (`--pairwise`), not the
  core substrate. Full `[N,N,T]` distance + vector when requested.
- `DftTarget`, `SourceSlot`, `NeighborhoodRecord` — the typed primitives.
- `FeatureSchema`/`FeatureColumn`, `RecordSink` (CSV with padded source
  slots via `QSaveFile`; the schema names every column + unit).
- `main_extract.cpp` — `QCommandLineParser` (`--run`, `--out`, `--case`,
  `--pairwise`); load; case loop; write.

## Reuse + the only edits to existing reader code

Reuse: `QtProtein`/`QtTopology`/`QtAtom`/`QtRing`/`QtBond`; `QtTrajectoryH5`
(positions, `bsShielding`/`mcShielding` TS, `ringNeighbourhood` TS +
`ring_membership_per_atom`, frame meta); `ConformationGeometry`
(`RingGeometryAt`, `Distance`, `AngleDegrees`); `io::QtProteinLoader`,
`io::DftShieldingLoader`, `io::CalcsetManifest`; `StructuredLogger`,
`QSaveFile`, `QCommandLineParser`.

Additive edits (append-only, viewer-safe):
1. **`DftAtomShielding` gains `Mat3 total_raw/dia_raw/para_raw`;
   `OrcaShieldingParser` keeps the raw 3×3 it already parses** — mandatory
   (the basis fix needs it).
2. `QtAtom::isAromaticRingHydrogen()` — new predicate.
   `QtAtom::IsBackboneAmideHydrogen()` **already exists** — use it.
3. `QtRing` accessors are `TypeIndex/LiteratureIntensity/NitrogenCount/
   JohnsonBoveyLobeOffset/Aromaticity/RingSizeValue/IsFused` — use the real
   names.

## Coherence

One extractor, shared vocabulary (`RunData` + local frame + geometry
primitives + the record). Critique surface: the `FeatureSchema`, the thin
`RediscoveryExtraction` interface, the basis fixture. The fitter is
deliberately not in this design — the substrate is built to outlive the
choice of fitter.

## Discipline

Plain-C++ data classes (no `QObject`); `StructuredLogger` UDP for progress;
`QSaveFile` output; `QCommandLineParser` args. New CMake target
`h5reader_extract` links `h5reader_core` + nanoflann; GUI binary untouched.
