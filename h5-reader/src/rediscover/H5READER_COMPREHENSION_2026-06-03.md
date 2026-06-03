# h5-reader comprehension snapshot — 2026-06-03

Session-start orientation. Jessica reordered the resume so that *understanding
the h5-reader substrate* (the CLI surface, the C++ functional programming, the
fast indexing, and above all the typed/data-rich protein model) came before the
process gate. This doc records what I read and what I now understand, so a
future session can re-acquire the same grip without re-reading the headers cold.

Not a spec, not a decision record. A comprehension capture. The governing
state-of-work doc is `STATE.md` (the 2026-06-03 handoff block); the architecture
contract is `NODE_STORE_CONTRACT_2026-06-02.md`.

## What I read (to ground each claim below)

- Rules / state: `h5-reader/CLAUDE.md`, `STATE.md:1-130` (the 2026-06-03 handoff),
  `FRESH_LOOK_2026-06-03.md`, `ALLATOM_FIT_SPEC_2026-06-03.md`.
- Typed model: `model/QtProtein.h`, `model/QtTopology.h`, `model/QtAtom.h`,
  `model/QtRing.h`, `model/QtFrame.h` (+ the `model/` and `io/` directory maps).
- Functional API: `rediscover/RediscoverTypes.h`, `rediscover/Verbs.h`,
  `rediscover/Relationship.h`, `rediscover/RelationshipEngine.h`,
  `rediscover/AnalysisBody.h`.
- Indexing: `rediscover/ResidentIndexes.h`, `rediscover/SpatialIndexSet.h`,
  `rediscover/TypedAtomIndex.h`, `rediscover/TemporalIndex.h`,
  `rediscover/FrameSpatialIndex.h`.
- CLI: `rediscover/main_extract.cpp` (case dispatch + the `per_atom_substrate`
  carrier at `:509-551`).
- Memories: `feedback_token_economy_codex_codes`, `project_between_calculator_network`,
  `feedback_python_complexity_is_the_primary_issue` (read this session) plus the
  named arc/vocabulary memories.

## The typed, data-rich protein model (the spine)

`QtProtein` is **identity + topology only** — non-copyable/non-movable (back-
pointers from conformations must never dangle); owns `QtAtom[]`, the parallel
`QtAtomNames[]` projection, `QtResidue[]`, and a `QtTopology`. Holds **no
geometry** — the wall is enforced in the type, not just the docs.

- **`QtAtom` IS the typed mechanical-identity tuple** — `(Element, Locant,
  BranchAddress, DiastereotopicIndex, BackboneRole)` plus the chemistry substrate
  (planar group, polar-H kind, ring-position labels, pseudoatom kind, **ff-atom-
  type as a typed enum, not a string**). It answers questions about itself:
  `IsBackboneNitrogen()`, `IsAromaticRingHydrogen()` (typed off HA/H4/H5 ff-types
  — the ring-current stratum), `IsAnyAlphaHydrogen()` (handles GLY HA2/HA3 per
  Markley). **No name strings on QtAtom** — amber/iupac/bmrb names live on
  `QtAtomNames` behind an explicit projection boundary (`QtProtein::atomNames(i)`).
  The reader is identity-driven, never label-driven. IUPAC-trap discipline is in
  the type itself.
- **`QtRing` is a class hierarchy, not enum+table.** Nine types (8 aromatic +
  `QtProPyrrolidineRing` saturated), each carrying its physics as virtuals —
  `LiteratureIntensity()` (−12.0 PHE, −11.28 TYR, −19.2 TRP9 indole perimeter,
  0.0 PRO), `JohnsonBoveyLobeOffset()`, `NitrogenCount()`, `Aromaticity()`,
  `RingSizeValue()`. `CreateQtRing()` is the *one* place ordinal→class is
  decided; all downstream code is ring-type-agnostic (`ring->LiteratureIntensity()`).
- **`QtTopology`** wraps bonds + rings + ring_membership and builds reverse-index
  caches at construction (CSR-equivalent: bonds-by-atom, ring-memberships-by-atom,
  memberships-by-ring in canonical-walk order). Fused TRP bridgeheads map to 2–3
  rings correctly.
- **Data-rich, per-frame**: `QtFrame` is a lightweight value view
  (`conformation* + tIndex`) into `TrajectoryConformation`'s per-TR buffers. It
  exposes every emitted physics channel per atom — `bs/hm/rsShielding`,
  `mcShielding`, `pqShielding`, `dispShielding`, `coulombShielding`/`apbsEfg`/
  `apbsEfield`, `aimnet2Shielding`, hbond/sasa/water/charges — each as a
  **`SphericalTensor` or `Vec3`, T2 preserved**, default-zero (= "no contribution")
  on an absent TR. The live model carries the full tensor surface, not scalars.

## The fast indexing — the resident spine

`Body{ const RunData& run, const ResidentIndexes& idx, const Catalog& catalog }`
is built **day-one, held immutable**; every closure captures it once.
`ResidentIndexes` =

- **`TypedAtomIndex`** — scoped, set-returning lookup over `QtAtom` *identity* via
  partial typed predicates; `selectUnique()` **fails loud on 0 or N matches**. Not
  name, not position — collision-safe, IUPAC-trap-safe.
- **`SpatialIndexSet`** — nanoflann KD trees, per cloud-kind × per-frame (5 clouds:
  Atoms / BondMidpoints / RingCenters / ChargeSites / AllBondMidpoints).
  `near`/`range` radius queries (cutoff REQUIRED at call site — no hidden cutoffs)
  return `SourceRef{cloud, cloud_index, entity_index}` — **provenance straight back
  to the typed entity.** (`FrameSpatialIndex` is the lazy per-frame
  anisotropic-bond-midpoint tree for McConnell.)
- **`RingGeometryCache`** — canonical normal-flipped ring geometry, day-one.
- **`TemporalIndex`** — no-copy frame-window arithmetic over atom-major storage
  (the Δ/rate combinator folds a contiguous span).

**This is the pair-index-as-pointer.** `SourceRef`/`SourceSlot` carry indices, not
copied pairs; queries are lazy radius lookups against resident trees. No resident
68 GB pair-dump — materialization is a transient, drop-old verb.

## The C++ functional programming — three layers

- **Layer 1 (Verbs, `Verbs.h`)** — pure thin wrappers that EXPOSE the resident
  spine and never regenerate it: `pos`, `at`, `window`, `value`/`valueVec3`/
  `valueTensor`/`valueT2`/`present`, `near`/`nearPoint`, `ringSlots` (frozen
  frame-0 membership backend), `ringGeom`, `atomOf`/`selectAll`, `heavyParent`,
  `ringsOf`/`ownRingAtoms`/`ownAromaticRing`, `displacement` (the ONE definition
  of the near-field (3cos²θ−1)/r³ geometry).
- **Layer 2 (Relationship, `Relationship.h`)** — curried closures bound into a
  named bundle defined IN CODE (`ComposedRelationships.cpp`), not a plugin
  registry: `Stratum`, `LocalFrameFn`, `SourceSelector[]`, `Attacher[]`,
  `ClassifierPrep`/`Classifier`, `SourceFilter`, `Reducer`, `TargetFn`,
  `BareKernelFn`. Building a relationship is partial application — config baked
  in, `(body, atom, frame)` arrive at iteration. Virtual dispatch on QtRing/QtAtom
  is the allowed polymorphism, inside attachers.
- **Layer 3 (Engine, `RelationshipEngine.h`)** — `RunTraversal` is ONE pure,
  order-free loop over the (atom, frame) index space. The carrier (`PerRecordSink`
  closure) is the ONLY thing that varies between drivers — the #29 unification:
  ring/mc/charge sum-sink, the broad reducer, and the per-atom-substrate carrier
  are all just different closures over the SAME walk. Referential transparency
  means the DFT-row-outer × stratum-atom-inner schedule is a CACHE choice, kept so
  emitted rows land byte-identical for the oracle gate.

The carrier types in `RediscoverTypes.h`: `DftTarget` (raw 3×3 total/dia/para +
library-basis decomposition + `total_local` in the per-atom local frame, with a
`local_frame_valid` flag — the T2 cross-frame caveat recorded, not assumed),
`SourceSlot` (ring/bond/charge superset, geometry in the target's local frame +
the source ring normal so the l=2 tensor can be reconstructed), `NeighborhoodRecord`
(identity + frame basis + un-summed sources + producer bare-kernel cross-check +
the DFT target).

## The CLI use

`h5reader_extract --run <calcset_dir_or_.LGS> --out <dir> --case <which>`, with
`which ∈ {ring_current, mcconnell, charge_dipole, broad_backbone,
all_atom_equivariant, per_atom_substrate, efg, buckingham_efield,
aimnet2_features, ring|mc|all}`. `per_atom_substrate` is explicitly **a carrier on
`RunTraversal`'s typed overload, NOT a RecordSink/WriteSource path → no
`*_sources.csv`** (pointer-index, not pair-dump). It writes
`per_atom_substrate_rows.csv` + NPY sidecars + a manifest stamping
`row_id == frame_slot * n_atoms + atom_index` and `raw_lab_frame`. The reader
(`main_reader.cpp`) is the Qt/VTK consumer of the same model; it never links the
library, never writes H5.

## The goal this rides toward

Recover the expected *classical* relationships from DFT on 1P9J; grade them by
statistical position (not "law"). The between-calculator network IS the
calibration + the calculator-inclusion gate. The typed model is the spine —
Python only fits and reads; it derives no physics value. Complexity is conserved,
not erased; the functional interface is the disciplined container that keeps
Python complexity knocked back.

## Where it stands + the loop (pointer; full detail in STATE.md)

v1 per-atom substrate committed + gated at `00ec168` (558,360 rows = 846 × 660).
Remaining arc (no new steps): all-atoms joint fit + partition → between-calculator
network → equations table (`equations/<mechanism>/`, pre-registered) → grade by
statistical position. Each step: ask-spec → vet-spec → execute → drift-assessment
+ postmortem; human-in-loop each loop; codex grinds, the lead vets/judges, the
lead's context is irreplaceable. Three fresh-look caveats to address: re-derive the
MOPAC-field r-values before citing; the all-atoms spec has LANDED (review it);
the combined-score stats are the old 54-atom prior and may move on the 846-atom
fit. Next action: review `ALLATOM_FIT_SPEC_2026-06-03.md` with the lead.
