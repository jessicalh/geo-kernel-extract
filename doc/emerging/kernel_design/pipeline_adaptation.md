# Pipeline adaptation for the irrep kernel feature set

Status: **design, not implementation.** Neutral. Deprecate-and-add. Geometry
layer untouched. This document answers one question: when the kernel-design
emissaries (`ring.md`, `mcconnell.md`, `charge_efg.md`) converge on each kernel
emitted as **clean equivariant spherical-tensor features (`0e` / `1o` / `2e`),
raw in the molecular frame, with uncertain physical scales carried as parallel
channels and left for the fit to weigh** — what must the EMIT and the READER
*adapt* to carry it, and where (if anywhere) is that more than a field-set
registration?

The short version, established by reading the code: **the substrate is already
shaped for this.** The catalog (SDK and C++) already carries per-array irrep
strings, parity, frame, and per-source/per-channel NPYs; the rediscover
functional interface already reads tensor- and irrep-valued fields uniformly
(`valueT2`, `valueTensor`, `valueVec3`); per-relationship schemas already declare
their own columns and sidecar NPYs. The new feature set is, structurally, **more
entries of the kinds already present** — registration on both sides, plus the
new per-kernel maths inside the producer's emit. The geometry layer (spatial
indexes, KD-trees, frames, displacement verbs) does not move.

---

## 0. What "the new feature shape" is, in this pipeline's terms

From the three design docs, the converged feature per kernel is:

- **ring** → symmetric rank-2 contribution decomposed to **`1x0e + 1x2e`**
  (isotropic ring-current shift + the CSA anisotropy), distributed-source field
  the production form, raw molecular frame.
- **mcconnell** → the clean symmetric-traceless propagator-on-Δχ as the primary
  **`2e`** (5 numbers), per source category, with Δχ left as a calibratable
  coefficient (geometry emitted unscaled in Å⁻³); `0e`/`1o` only where a real
  trace/asymmetry exists.
- **charge/efg** → **`1o`** (E-field) **⊕ `2e`** (EFG) together, both orders,
  in the molecular frame, with **parallel channels** for charge source
  (ff14sb / aimnet2 / mopac) and screening model (vacuum-Coulomb / APBS-PB /
  explicit-water), plus optional `0e` invariants.

Two cross-cutting properties the pipeline must carry:

1. **Raw molecular-frame, equivariant.** No imposed per-atom local frame; the
   geometry is emitted unscaled and equivariance handles rotation. The substrate
   already records *which* frame each tensor lives in (per-relationship frame
   metadata; `SURFACE_DESIGN.md` decision #7) and the all-atom emit is already
   molecular/lab-frame by construction.
2. **Parallel channels, fit-weighed.** Multiple charge sources and screening
   models are emitted side by side as separate channels, not collapsed. The
   substrate already does this: `charge_source` is a parameter axis (a sweep, not
   N relationships), and the all-atom carrier already emits ff14sb / aimnet2 /
   mopac / APBS / water as distinct arrays.

So nothing in the feature shape is foreign to the substrate. The adaptation is
naming the new per-kernel irrep blocks and registering them.

---

## 1. EMIT side (the easy part — one-way write)

The emit surface has two registration ledgers and one place where real maths is
added. The write is one-way (`feedback_build_inmemory_export_dont_relitigate`,
project CLAUDE.md): nothing downstream couples back into the extractor, so adding
emitted arrays is cheap and safe.

### 1a. The two ledgers a new array registers in

- **SDK catalog** — `python/nmr_extract/_catalog.py`. The single source of truth:
  one `ArraySpec` per NPY. It **already carries the irrep apparatus** the new
  feature set needs: `irreps` (`"0e + 1e + 2e"`, `"1o"`, `"1x2e"`, `"2e"`,
  `"mixed"`, `"256x0e"` are all already in use), `parity`, `tensor_rank`,
  `native_axis`, `units`, `sign_convention`, `mechanism`, `is_feature`. Adding a
  new irrep feature = adding `ArraySpec` lines. No schema change to `ArraySpec`
  itself is required — every column the design docs ask for already exists.
- **C++ resident catalog** — `h5-reader/src/rediscover/Catalog.h` (`ArrayId`
  enum + `ArraySpec`). This is what the functional interface reads *through*. It
  already has `ArrayRank { Scalar, Vec3, T2_5, Tensor9, Embedding256, RingNbhd4 }`
  and the matching `value / valueVec3 / valueT2 / valueTensor / valueEmbedding`
  accessors — i.e. it **already addresses `1o` (Vec3), `2e` (T2_5), and full
  `0e+1e+2e` (Tensor9) irrep-valued fields uniformly.** A new resident array that
  a relationship needs to *read* (e.g. a new charge-source field) is a new
  `ArrayId` + spec entry; its rank is one of the existing `ArrayRank`s. No new
  accessor verb is needed for `0e`/`1o`/`2e` — they exist.

### 1b. The two emit paths, and which one the new features use

- **Per-frame NPY / per-pose NPY** (`FrameNpyEmitter` → `ConformationResult::
  WriteAllFeatures`; modes 1–5). Each producer `*Result::WriteFeatures` writes
  its NPYs; the SDK catalog mirrors them. A new producer kernel array is a new
  `WriteFeatures` line + a new `ArraySpec`. **This is the canonical-extractor
  surface and is frozen during reader/feature work** (project CLAUDE.md: "the
  extractor is not to be modified during viewer or reader feature work"). So the
  new irrep features are **not** added here in a reader-feature session; they
  land when the producer kernels themselves are re-cut (see 1c), scheduled by
  the user, and the catalog grows in the same commit.
- **Rediscover substrate NPYs** (`h5-reader/src/rediscover/*Sink`; the
  per-relationship `<case>_sources.csv` / `<case>_aggregated.csv` + sidecar
  NPYs). This is where the equivariant feature substrate actually lives, and
  where the new irrep blocks are emitted for the e3nn fit. The carrier is
  **already** the right shape:
  - **`SourceSlot` already holds the literature-scaled, frame-correct per-source
    kernels and the parallel-channel handles** (`RediscoverTypes.h`):
    `ring_jb_unit_kernel` + `ring_jb_kernel` (unit-current and fixed-literature
    ring kernels, `SphericalTensor` in the target frame), `bond_mc_lit_kernel`
    (Δχ-scaled McConnell PCS tensor), `source_normal_local` (the ring dipole
    axis the `2e` is oriented by), and the charge fields (`source_q_e`,
    `source_charge_source`). The `0e`/`1o`/`2e` content is *already a
    SphericalTensor per source*.
  - **`FeatureSchema` already lets each relationship declare its own columns +
    sidecar NPYs** (`RecordSink.h`; `SURFACE_DESIGN.md` "Output carrier"):
    `SchemaFor(relationship)` builds the per-relationship schema; wide/array
    payloads (the 5-component T2, per-source matrices, embeddings) go to
    relationship-named NPYs, one `ArraySpec` each. The existing
    `broad_backbone_aggregated_ring_literature_kernel_T2` /
    `_bond_literature_kernel_T2` / `_charge_literature_kernel_T2` and the
    `all_atom_equivariant_apbs_efield` (`1o`) / `_apbs_efg_T2` (`2e`) /
    `_aimnet2_*` channels are exactly the template.

### 1c. What is genuinely new on emit, and where it goes

The only *new* emit work is the per-kernel maths the design docs name — and it
lives **inside the producer kernels / the rediscover attachers**, not in the
emit plumbing:

- **ring**: designate the distributed-source (Johnson–Bovey / Haigh–Mallion) +
  Boyd–Skrynnikov full-tensor construction as the production tensor and emit its
  explicit `1x0e + 1x2e` split as the canonical e3nn feature. The producer
  already computes BS/HM/RingChi tensors and already decomposes to T0/T1/T2; the
  delta is *which object is the canonical input* and emitting the `0e`/`2e` split
  as a named feature (per `ring.md` §5).
- **mcconnell**: emit the clean symmetric-traceless propagator-on-Δχ as the
  primary `2e` per source category, geometry unscaled (Δχ a calibration layer).
  The producer already emits `mc_category_T2` (per-bond-category `2e`) and the
  rediscover side already carries `bond_mc_lit_kernel`; the delta is making the
  clean `2e` the primary object rather than the mixed three-term `M`
  (per `mcconnell.md` §4).
- **charge/efg**: emit `1o` field ⊕ `2e` EFG together per charge-source and
  per-screening channel. Every piece already exists as an array
  (`apbs_E`/`apbs_efg`, `coulomb_E`/`coulomb_efg_*`, `water_efield`/`water_efg`,
  `aimnet2_efg*`, the rediscover `*_field_local` / `*_charge_literature_kernel_T2`);
  the delta is presenting the **pair `1o ⊕ 2e`** as the feature with the parallel
  channels kept distinct (per `charge_efg.md` §3).

**Deprecate-and-add framing:** none of the existing arrays are removed. The
mixed-rank `mc_shielding` (full `M`), the per-type `T0`/`T2`, the existing
`coulomb_*`/`apbs_*`/`water_*` channels all stay. The new canonical irrep
features are **added** as new `ArraySpec` entries (e.g. a ring
`*_ring_current_0e` + `*_ring_current_2e`, a McConnell `*_mcconnell_2e`, a
charge `*_field_1o` + `*_efg_2e` per channel). The fitter reads the superset and
chooses; the old fields remain for cross-check and continuity.

**Net on emit:** registration is two ledger edits per array (already-shaped
`ArraySpec`). The real work is the per-kernel construction inside the producer /
attachers, gated by the user (producer-freeze rule) — that is a *physics* change
the design docs scope, not a plumbing change.

---

## 2. READER side (config/registration vs real code)

The reader has three layers that touch fields. The new feature set lands as
**registration** in all three; **no geometry-layer code moves**, and the one
place a *new H5 time-series* would force hand-written code is explicitly *not* on
the path the equivariant substrate uses.

### 2a. The NPY load boundary — pure registration (generated)

`h5-reader/src/io/QtFieldCatalog.gen.h` is **generated** from `_catalog.py` by
`h5-reader/scripts/gen_field_catalog.py` (committed, stdlib-only, parses the
catalog with `ast`). The reader resolves a filename stem to a typed `FieldSpec`
exactly once (`FindFieldByStem`) and dispatches on `FieldKind`/`FieldGroup`/
`NativeAxis`; `irreps`/`parity`/`units`/`mechanism` ride along as interop
metadata. **New catalog entry → re-run the generator → the reader's NPY boundary
knows the field.** This is registration, not code: the generator is the adapter,
and it already emits the `irreps` string verbatim. No hand-edit of the boundary,
no per-field accessor to write.

### 2b. The rediscover functional interface — already carries irrep-valued fields

This is the core of the question "does the functional interface already carry
tensor/irrep-valued fields?" — **yes, it does, today.**

- The **resident Catalog** (`Catalog.h`) addresses every rank uniformly:
  `value` (`0e` scalar), `valueVec3` (`1o`/`1e` vector), `valueT2` (the
  5-component `2e`), `valueTensor` (full `0e+1e+2e` `SphericalTensor`),
  `valueEmbedding` (`256x0e`). The verbs (`Verbs.h`) expose exactly these:
  `value` / `valueVec3` / `valueTensor` / `valueT2` are Layer-1 primitives. An
  irrep-valued field is read with the matching verb; **no new verb is needed**
  for `0e`/`1o`/`2e`.
- A **Relationship** (`Relationship.h`) is a named bundle of curried closures;
  its `attachers` write `SourceSlot` fields and its `schema` (`FeatureSchema`)
  declares the emitted columns + sidecar NPYs. `SourceSlot` already carries
  per-source `SphericalTensor` kernels and the parallel-channel handles (§1b).
  Adding a new kernel feature to a relationship = adding an attacher that fills
  the already-present `SphericalTensor`/`Vec3` slot fields (or one new slot
  field of an existing type) and declaring its column/NPY in the schema. This is
  **new closures + a schema declaration in `ComposedRelationships.cpp`**, the
  same code shape every existing relationship already uses — direct named
  lambdas, the allowed "minimal/clarifying abstraction"
  (`feedback_functional_api_minimal_clarifying_abstraction`), not new grammar.

So: the interface **already carries** tensor/irrep-valued fields. The "real code"
on this layer is **per-relationship attacher + schema additions** to emit the new
canonical `0e`/`1o`/`2e` blocks and parallel channels — bounded, local, and of a
kind the engine already iterates. It is *not* a change to the engine loop, the
verbs, the catalog accessors, or the index/frame machinery.

If a new charge-source or screening field needs to be *read* by an attacher and
is not yet a resident array, that is **one new `ArrayId` + `ArraySpec` line** in
`Catalog.h` (rank picked from the existing `ArrayRank` set) — registration, with
the read served by the existing `valueVec3`/`valueT2` accessor.

### 2c. The trajectory H5 boundary — the one hand-coded surface (and why it's off-path)

`h5-reader/src/io/QtTrajectoryH5.h` is the **one** field surface that is *not*
catalog-generated: it has a hand-written `unique_ptr` slot + accessor per TR
group (`bsShielding()`, `apbsEfg()`, …). A genuinely **new per-frame H5
time-series** would need a new slot + accessor here — real code, the single place
the adaptation is more than registration.

**But the equivariant feature substrate does not travel this way.** The
rediscover substrate is built by the reader from the resident body and emitted as
**NPYs + per-relationship CSV** (`SURFACE_DESIGN.md`: "the substrate emits NPYs";
Python consumes the substrate, never `trajectory.h5`). The new irrep features are
substrate NPYs, read through the resident `Catalog`, not new H5 time-series. So
**`QtTrajectoryH5` does not need to change** to carry the new feature set —
*unless* the user separately decides a new kernel's per-frame tensor should also
become a first-class animated H5 time-series for the viewer (a viewer-display
decision, distinct from the e3nn feature substrate, and a producer/fileformat
change the user schedules). For the feature set as scoped by the three design
docs, this boundary is untouched.

### 2d. Geometry layer — explicitly untouched

The spatial indexes (`SpatialIndexSet`, per-cloud per-frame KD-trees), the
`RingGeometryCache`, `TypedAtomIndex`, `TemporalIndex`, `LocalFrameBasis`, and
the `displacement` verb are the geometry layer. The new feature set changes
*what tensor is attached at a source*, not *how sources are found or framed*. The
raw-molecular-frame stance (no imposed per-atom frame; equivariance handles
rotation; `feedback_frames_from_physics_not_tests`) is the frame the substrate
already records and the all-atom emit already uses. **Nothing in the geometry
layer moves.**

---

## 3. The minimal defensible change — summary table

| Surface | Adaptation | Kind |
|---|---|---|
| SDK catalog `_catalog.py` | new `ArraySpec` per irrep block / channel (irrep apparatus already present) | **registration** |
| C++ resident `Catalog.h` | new `ArrayId` + `ArraySpec` only for any new array a relationship *reads*; rank from existing set | **registration** |
| Reader NPY boundary `QtFieldCatalog.gen.h` | re-run `gen_field_catalog.py` | **regeneration (registration)** |
| Functional interface verbs / engine | none — `valueT2`/`valueTensor`/`valueVec3` already carry `0e`/`1o`/`2e` | **none** |
| Rediscover relationships `ComposedRelationships.cpp` | new attacher(s) filling existing `SphericalTensor`/`Vec3` slot fields + schema column/NPY declaration | **real code (local, bounded; the design intent)** |
| Producer kernels (`ring`/`mcconnell`/`charge`) | construct & emit the canonical `0e`/`1o`/`2e` object per the design docs | **real code (physics; user-scheduled, producer-freeze)** |
| Geometry layer (indexes, KD-trees, frames) | none | **untouched** |
| Trajectory H5 boundary `QtTrajectoryH5.h` | none for the feature substrate (NPY path); only if a kernel becomes a new animated H5 TR (separate viewer decision) | **none on-path** |

**Deprecate-and-add throughout:** existing `*_shielding` / per-type / mixed-`M`
/ per-channel arrays stay; the canonical irrep features are added; the fitter
reads the superset.

**The single place it is more than a field-set registration:** the
**rediscover-side per-relationship attachers + schema** (`ComposedRelationships.cpp`
/ `FeatureSchema`), where the new canonical `0e`/`1o`/`2e` blocks and parallel
channels are *constructed and declared* for emit — and, upstream of that and
user-gated, the **producer kernel maths** the three design docs specify. Both are
the design's actual content (the physics object and its equivariant emit), not
plumbing; everything else — the catalog ledgers, the generated reader boundary,
the functional read verbs — is registration the substrate was already built to
absorb.
