# Handoff to codex — rediscover functional API (compile + oracle gate)

Branch `h5-reader-pysr-spike` (NEVER MERGE). Authored by the functional-API
agent on 2026-06-01. **The compiler/ninja was NOT reachable from this agent's
Bash** (sandbox denied `cmake --build`, even `which cmake`). So the code below
is authored + rigorously self-reviewed in place of a build; codex takes the
compile-fix-run loop + the oracle gate. This answers the open sandbox-knob
question: agent Bash here cannot invoke the compiler.

## What was authored (the functional API — SURFACE_DESIGN Layers 1–3)

New files (all under `src/rediscover/`, added to `h5reader_extract` in
`CMakeLists.txt`):

- **`Verbs.h` / `Verbs.cpp`** — Layer-1 primitive access verbs (pure, return
  views; THIN over the existing spine — they expose the indexes, don't rebuild
  them): `pos`, `at`, `window`, `near` / `nearPoint`, `value` / `valueVec3` /
  `valueTensor` / `valueT2` / `present`, `slots` (H5 ring-neighbourhood
  frozen-membership backend), `ringGeom`, `atomOf` / `selectAll` (wrap
  `TypedAtomIndex`), `ringsOf` / `ownRingAtoms` / `ownAromaticRing`.
- **`Relationship.h` / `Relationship.cpp`** — Layer-2 combinators AS CURRIED
  CLOSURES (`Stratum`, `LocalFrameFn`, `SourceSelector`, `Attacher`,
  `ClassifierPrep` + `Classifier` (SourceClassifier), `SourceFilter`, `Reducer`,
  `TargetFn`, `BareKernelFn`) bound into a `Relationship` named bundle. Curried
  verb builders: `atomsWhere(pred)`, `slotsBackend()`, `nearBackend(cloud,
  cutoff)` — each returns a closure capturing config; the body arrives at
  iteration.
- **`RelationshipEngine.h` / `RelationshipEngine.cpp`** — Layer-3: one pure,
  order-free loop `RunRelationship(rel, body, sink)`, exactly the SURFACE_DESIGN
  snippet (frame_fn → flatten selectors → attach per source [→ filter →
  classify] → reduce → sink.fold). Schedule: DFT-row-outer × stratum-atom-inner
  (matches the oracle's emission order so the gate file-compare is exact).
- **`ComposedRelationships.h` / `ComposedRelationships.cpp`** —
  `MakeRingCurrentRelationship()` and `MakeMcConnellRelationship(cutoff)`:
  ring_current and mcconnell re-expressed as composed `Relationship`s built from
  the verbs + curried closures, reproducing the procedural cells.
- **`analysis/oracle_parity.py`** — the gate runner: runs the extractor with
  `--engine composed` and `--engine procedural` into two dirs and diffs every
  CSV (numeric tol) + sidecar NPY. Exit 0 ⇒ byte-faithful.

Edited:
- **`main_extract.cpp`** — added `--engine {composed|procedural}` (default
  composed). ring/mc now run through the composed engine by default; the
  procedural cells run under `--engine procedural` (the reference oracle).
  charge_dipole unchanged (procedural; out of scope for the API).
- **`CMakeLists.txt`** — the 8 new files added to `h5reader_extract`.

The procedural cells (`RingCurrentNeighborhood.*`, `McConnellNeighborhood.*`)
are UNCHANGED and kept as the reference oracle.

## Untouched (discipline confirmed)
- The spine (`Catalog`, `TemporalIndex`, `TypedAtomIndex`, `SpatialIndexSet`,
  `RingGeometryCache`, `ChargeStore`, `ResidentIndexes`, `AnalysisBody`,
  `RecordSink`, `ExtractionSupport`, `RunData`, `SphericalBasis`,
  `LocalFrameBasis`) — verbs wrap it, nothing duplicated/regenerated.
- The GUI reader, the `nmr_shielding` library (`src/` at repo root), `RecordSink`
  serialization. Topology is reused from the resident indexes, never re-parsed.

## The build command (codex runs this)
```
cmake --build /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-gcc \
      --target h5reader_extract h5reader_rediscover_tests
ctest --test-dir /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-gcc \
      -R h5reader_rediscover_tests
```
`h5reader_rediscover_tests` is the SphericalBasis fixture only (unchanged —
should still pass). The functional API has no standalone unit test target (a
real test needs the full Body + a calcset on disk = the oracle gate below); I
did NOT expand the test target's link surface without a compiler to verify it.

## The oracle gate (codex runs this after a green build)
1. Byte-parity composed vs procedural (the API-is-faithful proof):
```
python src/rediscover/analysis/oracle_parity.py \
  --bin  build/linux-gcc/h5reader_extract \
  --run  /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft \
  --work /tmp/rediscover-parity --case all --mc-cutoff 8.0
```
Expect: every CSV + NPY identical (exit 0). Divergence = a composition bug in
ComposedRelationships.cpp (NOT a new result) — diff prints the column/row.

2. The physics oracle numbers (composed output, venv
`src/rediscover/analysis/venv`):
```
build/linux-gcc/h5reader_extract --run <calcset> --out /tmp/rdc-composed --case all
# ring: leave-atoms-out k≈21 ppm·Å³, coupled within-atom R²≈0.62
python src/rediscover/analysis/credibility2_instantaneous.py /tmp/rdc-composed
# equivariant T2 R²≈0.44, |T2| r≈0.75, basis ~4.9e-8  (torch is cu130 — put
#   nvidia/cu13/lib in LD_LIBRARY_PATH first or it segfaults)
python src/rediscover/analysis/equiv_t2.py /tmp/rdc-composed
# mcconnell: scalar R²≈0.55, kernel readout r≈0.918
python src/rediscover/analysis/sumpool_mcconnell.py /tmp/rdc-composed
```
The driver log also prints the DFT frame check (rotation ~1e-4°).

## Self-review (in place of the build) — traced, believed correct

**Types / signatures (all copied verbatim from the procedural cells or the
spine headers, so correct by construction):**
- `verbs::*` call exactly the spine methods the cells call: `valueVec3(Positions)`
  = `conf.atomPosition(frame, atom)`; `spatial.near(cloud, frame, query, cutoff)`;
  `ringGeometry.at(ring, frame)`; `typedAtoms.selectUnique(scope, sel, &err)`;
  `topo.ringMembershipsForAtom/ringMembershipAt/ringAt`; the ring-neighbourhood
  TS `ringNeighbourhood()/n_slots/ringIndexAt/at`. Arg ORDER verified against
  each header.
- `Relationship.h` forward-declares `struct QtAtom` (it IS a struct, QtAtom.h:33)
  for the `atomsWhere` signature; `RawSource` holds a complete `verbs::RingSlot`
  (Verbs.h is included). All closure types are `std::function` (the cool levels,
  per SURFACE_DESIGN); attachers are plain function pointers bound into them.
- `BuildTarget(body.run, atom, orig, frame)` — `body.run` is `const RunData&`,
  matches `BuildTarget(const RunData&, …)`.
- `main_extract` `RunnableCase` holds schema by value + a `std::function` runner
  capturing a `shared_ptr` (keeps the relationship/cell alive); calls const
  methods. `#include <functional>` added.

**Closure captures:** `nearBackend(cloud, cutoff)` captures cloud+cutoff by
value; `slotsBackend()` captures nothing (H5 supplies membership); `atomsWhere`
moves the predicate in; the mc reducer/filter lambdas capture `cutoff_A` by
value. The body is captured by NO closure — it is passed to each closure by the
engine as `const Body&` (the SURFACE_DESIGN "C++ reality": closures capture
config; the engine threads the body). [Design note: SURFACE_DESIGN says closures
"capture const Body&"; I pass it as a parameter instead so a Relationship is
body-independent and reusable across runs — same computation, the engine does
the body-threading. Flagged as a deliberate, faithful reading, not a procedural
revert.]

**Engine loop:** order-free per referential transparency; the chosen schedule
matches the cells' (dftRows outer, stratum inner) so emitted-row order is
identical — required for the byte-parity gate. frame_fn → flatten selectors →
(attach, filter, classify) per source → reduce → WriteSourceRows +
WriteAggregatedRow. Identical sink calls to the cells.

**Faithfulness deltas handled (each preserves byte-parity):**
1. ring `is_self_or_bonded`: the cell sets it inline before pushing the slot. In
   the composed path the attacher leaves it default; a `Classifier` (run before
   the row is written, with a per-atom `ClassifierPrep` building the own-ring /
   own-atom sets) stamps it — so it lands on the per-source row identically. The
   reducer then reads it off the slot.
2. mcconnell geometric rejects (degenerate axis / r≤1e-6 / own-H): the cell
   `continue`s (no row). The composed attacher marks the slot dipolar=NaN and a
   `SourceFilter` drops it before it becomes a row — so no extra NaN source row,
   and `n_sources == n_sources_valid` as in the cell.
3. ring H5-slot rejects (non-finite / distance≤1e-6): `verbs::slots` applies the
   cell's exact filter, so those slots never become RawSources (no row, not
   summed) — identical.
4. mcconnell `residueIndex<0`: the cell `continue`s the whole atom (no case).
   residueIndex is static, so I folded it into the stratum predicate
   (`IsBackboneAmideHydrogen() && residueIndex>=0`) → zero cases for such an
   atom, identical. (No backbone HN lacks a residue in practice.)
5. Source ORDER: ring = same H5 slot walk; mc = same `spatial.near` →
   nanoflann radiusSearch (deterministic, same tree+query) → same order.

**Risks to watch on first build (where a typo would surface):**
- `model::kAromaticRingTypeCount` is in Types.h (namespace model) = 8 — used as
  `model::kAromaticRingTypeCount`. The 8 ring-type / 4 bond-category per_type
  vector sizes must match the schema column counts (they do: 8
  `sum_dipolar_ringtype_*`, 4 `sum_dipolar_<cat>`).
- `RingTypeIndex::TrpBenzene`, `Locant::Gamma/Delta`, `Element::C`,
  `BondCategory::{PeptideCO,PeptideCN,SidechainCO,Aromatic}` — all used verbatim
  from the cells (exist).
- Eigen `Vec3` ops (`.norm()`, `.normalized()`, `.dot()`, `operator-`, scalar
  `*`) used exactly as the cells use them.

If the byte-parity gate diverges, the divergence is localized to one of the 5
deltas above — start there.
