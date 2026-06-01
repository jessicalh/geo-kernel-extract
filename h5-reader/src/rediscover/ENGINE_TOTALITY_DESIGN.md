# Engine totality design: fold/sink carrier

Status: proposal for lead + user review. No implementation in this pass.

Context: the broad-backbone run exposed that `RunRelationship` already has the
right closure traversal shape, but its tail is hardwired to the ring/McConnell
scalar-sum carrier:

```text
SourceSet -> AggregateResult -> RecordSink::WriteSourceRows()
                         \-> RecordSink::WriteAggregatedRow()
```

That forced `RunBroadBackbone`, a sibling runner with the same loop and a
different reducer/sink. The fix should cut that duplication by parameterizing the
loop tail, not by growing a scenario framework. Relationships stay direct, named
bundles defined in code.

## 1. Shape enumeration

The table below covers the agreed 9 relationship types from `SURFACE_DESIGN.md`
plus the current broad case. "Inferred" marks items where the current code or
stub notes fill a gap not fully settled by `SURFACE_DESIGN.md`.

| Relationship | `relationship_kind` | Source kinds | Reducer-output shape | Row kinds | Target | Cross-check kernel | Special carrier need |
|---|---|---|---|---|---|---|---|
| `ring_current` | `source_sum` | Aromatic rings from H5 `ring_neighbourhood` slots; source identity from `QtRing` | Scalar dipolar sums: `sum_dipolar_all`, producer-valid sum, per-ring-type sums, counts | Per-source + aggregated | DFT tensor for aromatic ring H; local aromatic-H frame | Present: producer `bs` kernel | Legacy byte-parity carrier; per-source ring geometry; target/bare-kernel T2 sidecar NPYs |
| `mcconnell` | `source_sum` | Anisotropic bonds from KD `bond-midpoints` | Scalar dipolar sums and per-bond-category sums | Per-source + aggregated | DFT tensor for backbone HN; HN frame | Present: producer `mc` kernel | Legacy byte-parity carrier; bond-axis local vector; target/bare-kernel T2 sidecar NPYs |
| `buckingham_efield` | `per_atom_feature` | APBS E-field at the target atom; no outward source cloud | Vec3 l=1 feature, plus scalar projections if declared by the relationship | Per-atom | DFT tensor for atoms; lab/MD frame unless a local-frame convention is chosen | Absent. APBS E-field is source data, not a producer bare shielding kernel | Vec3 payload should be NPY or compact per-atom columns; T1 target use remains flagged unverified |
| `efg` | `per_atom_feature` | APBS EFG T2 at the target atom | T2-5vec feature | Per-atom | DFT T2 for atoms in lab/MD frame | Absent per `SURFACE_DESIGN.md`: APBS EFG is the source, not a cross-check | T2 source payload sidecar NPY; no degenerate per-source CSV |
| `charge_dipole` | `source_sum` | Charge-site KD cloud; `charge_source` in `ff14sb`, `aimnet2`, `mopac` | Vec3-l1 dipole `mu = sum q_i (r_i - r_atom)`; optional norm/projections | Per-source + aggregated, but source emission must be policy-controlled | DFT tensor for atoms; lab-frame equivariant path or declared local frame | Absent | Charge source is required and recorded; high-cardinality source rows need compact NPY or explicit opt-out/opt-in policy |
| `charge_quadrupole` | `source_sum` | Charge-site KD cloud; same `charge_source` axis as `charge_dipole` | T2-5vec / traceless rank-2 charge quadrupole | Per-source + aggregated, but source emission must be policy-controlled | DFT T2 for atoms; lab-frame equivariant path or declared local frame | Absent | T2 aggregate sidecar NPY; charge source rows have same scale problem as `charge_dipole` |
| `larsen_hbond` | `source_sum` | H-bond acceptor/donor geometry around exchangeable H; detection method required | Scalar + T2 features by Larsen class or acceptor class (inferred) | Per-source + aggregated | DFT tensor for exchangeable H; HN frame for amide H, unresolved frames for other donors | Unresolved. Stub notes four Larsen H-bond time series exist, but mapping to one relationship kernel is not settled | Needs reviewer call on detection mode, donor stratum, non-HN frames, acceptor classification, and cross-check mapping |
| `charge_response_gradient` | `per_atom_feature` | AIMNet2 charge-response-gradient scalar + Vec3 at the target atom | Scalar + Vec3 feature; not polarizability | Per-atom | DFT T0/T2; T1 if/when target convention is resolved | Absent | Pair of per-atom payloads; manifest must name it CRG, not polarizability |
| `aimnet2_embedding` | `per_atom_feature` | AIMNet2 per-atom 256-d embedding | 256-d learned feature vector | Per-atom | DFT tensor for atoms | Absent | Embedding must be sidecar NPY keyed by per-atom row; no 256-column CSV and no `SourceSlot` stuffing |
| `broad_backbone` | `source_sum` | Heterogeneous KD sources: ring centers, bond midpoints, charge sites | Heterogeneous per-mechanism aggregate: ring scalar sum, bond scalar sum, charge field Vec3/`field_z`/norm plus `mu` | Per-source + aggregated | DFT tensor for backbone N/CA/C/O/HN/HA; backbone-class local frames | Absent | Target-once carrier already required; per-charge source CSV reaches 1.7-5.5 GB and needs compact NPY or optional source emission |

Two conclusions from the table are load-bearing:

1. `source_sum` and `per_atom_feature` are both first-class. The per-atom cases
   should not be modeled as fake `self()` source clouds in the output carrier.
2. Reducer output is not one scalar sum. The known closed set includes scalar
   sums, Vec3, T2-5vec, scalar+Vec3, 256-d embedding, and heterogeneous
   per-mechanism records.

## 2. Proposed minimal abstraction

Keep one engine loop and parameterize only the part that is not traversal:

```text
for (dft_row, atom) in selected index space:
    frame  = frame_fn(body, atom, row)
    base   = identity + frame + target + optional cross-check

    if relationship_kind == source_sum:
        sources = flatten(selectors(body, atom, row))
        slots   = attach/classify/filter each source
        record  = reducer(base, SourceSet(slots))
    else if relationship_kind == per_atom_feature:
        value   = feature_reader(body, atom, row, frame)
        record  = reducer(base, PerAtomValue(value))

    sink.fold(record)
```

This is still the existing Layer-3 closure traversal. The change is that the loop
does not know `AggregateResult`, `WriteSourceRows`, or
`WriteAggregatedRow`. It only constructs the common per-case base and delegates
the relationship-specific tail:

- `RelationshipTraversal`: name, kind, schema metadata, stratum, frame function,
  target function, optional cross-check, and either a source pipeline or a
  per-atom feature reader.
- `Reducer`: a direct named function/lambda that maps `SourceSet` or
  `PerAtomValue` to a typed record.
- `Sink`: a direct named object/lambda with `fold(record)`, responsible for CSV
  rows, NPY sidecars, row ids, source-row policy, and manifest-facing counts.

This can be implemented as a templated free function or an equivalent direct
function with a closed record variant. The design does not require abstract base
classes, registration, factories, config-string dispatch, or a plugin kit.

### Record shapes to support

The engine needs no open-ended "any payload" type. The total known set is small:

| Record shape | Relationships |
|---|---|
| `ScalarSourceSumRecord` | `ring_current`, `mcconnell`; possibly scalar portions of `larsen_hbond` |
| `Vec3SourceSumRecord` | `charge_dipole` |
| `T2SourceSumRecord` | `charge_quadrupole` |
| `Vec3FeatureRecord` | `buckingham_efield` |
| `T2FeatureRecord` | `efg` |
| `ScalarVec3FeatureRecord` | `charge_response_gradient` |
| `EmbeddingFeatureRecord` | `aimnet2_embedding` |
| `HbondSourceSumRecord` | `larsen_hbond`, once detection/classes are settled |
| `HeterogeneousMechanismRecord` | `broad_backbone` |

These are not a generic tensor framework. They are the minimal typed carriers
needed by the known relationships.

### Mapping every relationship onto the loop

| Relationship | Traversal input | Reducer | Sink |
|---|---|---|---|
| `ring_current` | Current aromatic-H stratum, aromatic frame, `slotsBackend`, ring attacher/classifier, `bs` cross-check | Current `ringReducer` -> `ScalarSourceSumRecord` | Legacy scalar-sum sink adapter first, preserving byte parity |
| `mcconnell` | Current HN stratum, HN frame, bond KD selector, bond attacher/filter, `mc` cross-check | Current `mcReducer` -> `ScalarSourceSumRecord` | Legacy scalar-sum sink adapter first, preserving byte parity |
| `buckingham_efield` | All-atom stratum, lab frame, per-atom `ApbsEfield` reader | Pass-through/projection reducer -> `Vec3FeatureRecord` | Per-atom feature sink; Vec3 as sidecar NPY or compact columns |
| `efg` | All-atom stratum, lab frame, per-atom `ApbsEfg` reader | Pass-through reducer -> `T2FeatureRecord` | Per-atom feature sink; T2 sidecar NPY |
| `charge_dipole` | Charge KD selector with required `charge_source` and cutoff | Sum `q_i * disp` -> `Vec3SourceSumRecord` | Source-sum Vec3 sink with source emission policy |
| `charge_quadrupole` | Charge KD selector with same charge-source axis | Sum traceless `q_i * (3 rr - r2 I)` -> `T2SourceSumRecord` | Source-sum T2 sink with source emission policy |
| `larsen_hbond` | Exchangeable-H stratum, detection-specific acceptor selector, donor frame | Sum/classify acceptor geometry -> `HbondSourceSumRecord` | H-bond source-sum sink once class schema and cross-check are decided |
| `charge_response_gradient` | All-atom stratum, lab frame, AIMNet2 CRG scalar+Vec3 reader | Pass-through reducer -> `ScalarVec3FeatureRecord` | Per-atom feature sink; scalar column + Vec3 sidecar/columns |
| `aimnet2_embedding` | All-atom stratum, lab frame, AIMNet2 embedding reader | Pass-through reducer -> `EmbeddingFeatureRecord` | Per-atom feature sink; 256-d sidecar NPY keyed by row id |
| `broad_backbone` | Current broad stratum, backbone frames, heterogeneous KD selectors and attachers | Current broad reducer -> `HeterogeneousMechanismRecord` | Broad sink behavior, but driven by the shared loop |

Nothing extra:

- No per-relationship runner.
- No registry or factory.
- No generic schema language beyond each sink declaring its CSV columns and NPY
  payloads.
- No fake `self()` selector for per-atom features in the output model.

Nothing missing:

- Scalar source sums still cover ring/McConnell and preserve the oracle path.
- Heterogeneous source sums cover broad.
- Source-sum Vec3 and T2 cover charge multipoles.
- Per-atom Vec3/T2/scalar+Vec3/embedding cover APBS, CRG, and AIMNet2 embedding.
- Optional cross-check is attached to the base record only when real producer
  kernel data exists.

### Where the abstraction is genuinely minimal

The loop has exactly one reason to vary: the value produced after traversal and
where that value is written. Reducer + sink are therefore the minimal parameters.
The `relationship_kind` branch is also necessary because per-atom features are
not source sets and should not be expressed as source-row artifacts.

Things deliberately cut:

- A plugin/registry/factory framework.
- One superset CSV or one superset `SourceSlot` that tries to hold embedding,
  T2, Vec3, ring, bond, and charge fields at once.
- A general tensor algebra DSL.
- New schedulers, parallel runtime, cancellation protocol, or config framework.
- Automatically generated relationships from JSON.

## 3. Per-source scale

Broad confirmed the target-repeat fix was necessary but not sufficient. Even
with target columns removed from source rows, broad source CSVs reached
1.7-5.5 GB across the 6/10/12 A sweep; charge alone produced 8.1M rows at 6 A.

Proposal:

1. Every relationship always emits the per-case aggregate/per-atom row and its
   declared NPY payloads.
2. Per-source emission becomes an explicit sink policy, recorded in the manifest:
   `none`, `compact_npy`, or `debug_csv`.
3. `compact_npy` is the default for high-cardinality charge source rows. Use
   columnar or structured NPY payloads keyed by `row_id` and `source_row_id`;
   keep CSV for human-readable low-cardinality/debug use only.
4. `debug_csv` has a documented source-row budget. The budget is a fail-loud
   preflight or early-run guard, never truncation. If the estimate/count exceeds
   the budget, the run fails with a message naming `compact_npy` or `none`.
5. Ring/McConnell keep their legacy CSV mode during the oracle migration. Do not
   change their byte-parity output in the same step that changes the engine.

This gives the fitter access to source rows when needed without silently writing
multi-GB CSVs or silently dropping data. It also keeps the broad/charge cases
from turning the source CSV into the de facto storage format for array-shaped
scientific data.

## 4. Migration plan

### Phase 0: review this design

No C++ changes until lead + user sign-off on the fold/sink boundary and the
source-row policy.

### Phase 1: factor the loop without changing outputs

Introduce the shared fold loop behind the current behavior:

- Keep current `Relationship` closure traversal semantics.
- Wrap ring/McConnell with a scalar reducer and a legacy `RecordSink` adapter.
- The adapter writes exactly the same source CSV, aggregated CSV, and T2 sidecar
  NPYs as today.
- Re-run the existing composed-vs-procedural oracle gate. Byte parity must remain
  exact.

Containment: the procedural cells stay untouched as the oracle, and the first
change is only a loop-tail refactor with identical sinks.

### Phase 2: fold broad into the shared loop

Route current `broad_backbone` through the shared fold loop using its existing
source pipeline, broad reducer, and broad sink behavior.

Gate:

- Counts and schema stable against current broad output.
- Broad source rows still have no `dft_*` target columns.
- Frame classes remain 0% invalid.
- Per-atom-type sigma_iso diagnostics stay stable within the existing
  correlate-not-match tolerance.

Containment: do not change broad carrier format in the same step as removing
`RunBroadBackbone`. First prove the loop unification.

### Phase 3: migrate `charge_dipole`

Route the current procedural charge-dipole case through the shared loop as a
`Vec3SourceSumRecord`. Preserve current output first, then apply the approved
source-row policy in a separate schema-versioned step.

Gate:

- Aggregated `mu` values stable.
- Source count stable.
- No charge-source fallback; absent source remains a ValidateScenario error.

### Phase 4: implement the per-atom feature branch

Bring up `efg` first because it exercises `per_atom_feature` and T2 sidecar
payloads. Then add `buckingham_efield`, `charge_response_gradient`, and
`aimnet2_embedding`.

Gate:

- Fail-loud when required arrays are absent.
- Per-atom rows align with DFT row order.
- Array sidecars have exact documented shapes.
- Manifest records frame, units, source array, and unverified T1/T2 flags where
  applicable.

### Phase 5: finish remaining source-sum shapes

Add `charge_quadrupole` after `charge_dipole` proves charge-source handling.
Add `larsen_hbond` only after reviewer decisions on detection, donor frames,
acceptor classes, and optional cross-check mapping.

## 5. Alternatives and trade-offs

### Alternative A: keep one runner per carrier

This is the current broad outcome: `RunRelationship` for scalar sums and
`RunBroadBackbone` for heterogeneous broad. It is easy locally, but it duplicates
the traversal loop and guarantees another runner for per-atom T2, embedding, and
charge multipoles. Rejected because it is exactly the multiplicity without
necessity the broad case exposed.

### Alternative B: one giant `NeighborhoodRecord` plus one giant sink

This would add fields for Vec3, T2, 256-d embedding, CRG, hbond classes, broad
field, and every possible source identity, then teach one sink to switch on
relationship name. It avoids function parameterization but creates a superset
carrier, wide schemas, and dead fields in most rows. Rejected because it repeats
the one-off carrier mistake at larger scale.

### Alternative C: relationship registry/factory with polymorphic sinks

This would make each relationship an object registered under a string and loaded
through a factory. It is flexible, but flexibility is not the need here. The
known relationship set is finite and named in code. Rejected by project rule and
by unnecessary machinery.

### Trade-off accepted

Typed record shapes mean adding a genuinely new carrier shape later may require
one new record and one new sink. That is acceptable: it is an explicit addition
to a closed scientific surface, not a hidden plugin system. The current known
set is covered without inventing carriers for unknown future cases.

## 6. Open questions for review

1. Source-row policy for charge-heavy cases: should default be `compact_npy` or
   `none`? What row/byte budget should make `debug_csv` fail loud?
2. Ring/McConnell carrier after migration: keep the legacy target-repeating
   source CSV indefinitely for oracle stability, or introduce a schema-v2
   target-once carrier after byte-parity gates are retired?
3. L1 target handling: do `buckingham_efield`, `charge_dipole`, and CRG target
   shielding T1, T0 only, or T0 plus T2 invariants until the T1 basis/parity
   issue is resolved?
4. All-atom frame convention: for APBS/charge/embedding relationships, is lab
   frame sufficient as the primary emitted frame, or should generic per-atom
   local frames be built before those cases become canonical?
5. `larsen_hbond`: which detection mode is canonical, which donors are in the
   stratum, what frames do non-HN donors use, and do the four existing Larsen
   time series map to a real cross-check payload?
6. For `aimnet2_embedding`, should the 256-d NPY be emitted by default when the
   relationship is requested, or require an explicit "large feature" opt-in in
   the same spirit as the Python loader's optional-large policy?

