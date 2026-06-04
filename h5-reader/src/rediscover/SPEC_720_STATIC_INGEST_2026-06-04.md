# SPEC 720 Static Ingest - 2026-06-04

Status: DESIGN ONLY, pending lead plan-vet. Nothing here is built by this pass.

Deliverable scope: one design document for the 720-WT static-pose ingest. No code, build, extraction, 720 run, ORCA run, fit, substrate edit, or git action is part of this pass. The brief makes this a spec-authoring pass only (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:5`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:110`).

## 0. Ground Truth Frame

**BUILT:** The current rediscover arc treats 1P9J as a within-axis instrument and explicitly defers between-axis probability to the 720-WT statics pilot. `NEXT_SESSION_PROMPT.md` says the 720-WT statics pilot is the B-path and the "ONLY between instrument" seed of the unified stats engine (`NEXT_SESSION_PROMPT.md:42`, `NEXT_SESSION_PROMPT.md:46`). `NOW.md` states that 1P9J's true between axis is null/retracted and the 720-WT owns between/transferability (`NOW.md:71`, `NOW.md:88`). `STATE.md` likewise records that 1P9J has no clean between axis and that between defers entirely to 720-WT (`STATE.md:69`, `STATE.md:75`).

**GIVEN:** The canonical input is the curated 720 mutant set, not the stale 725 duplicate-containing set in `/shared/2026Thesis/consolidated/`. The exact path is TBD and must not block the design (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:29`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:42`).

**GIVEN:** The curated 720 proteins are static and DFT-bearing. The upstream producer action is an `nmr_extract` re-run over each static pose with MOPAC ON, while reusing the existing r2SCAN DFT shielding targets; there is no new ORCA/DFT calculation. Per-protein producer output exists, keyed by protein ID: topology, positions, APBS/AIMNet2/MOPAC fields, kernel shadows, and raw DFT dia/para/total tensors as the WT absolute shielding target (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:35`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:42`).

**DESIGNED:** The static ingest consumes that per-protein producer output. It does not invoke, re-run, or re-implement `nmr_extract`; the MOPAC re-run is upstream producer work, and the extractor remains sacred and unmodified (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:35`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:40`; `GUIDANCE.md:17`, `GUIDANCE.md:20`).

## 1. The Foreclosure: Law, Not Advice

This is the headline constraint. The design must make the following two failure modes structurally impossible, not merely discouraged.

**LAW 1: no terabyte of Python.** The static ingest must not emit bulk pairwise, per-source, per-protein-model, or occupancy materializations for Python. The old pairwise/occupancy dump pattern is explicitly the anti-pattern; the 720 version would be a terabyte-scale failure (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:44`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:49`). The only default output is a lean keyed per-(protein, atom) digest under 15 GB total. All reductions, source association, geometry, dominance, bins, and between-axis accumulations happen in C++.

**LAW 2: no second protein model in Python, not even an aggregate one.** Python never owns proteins, atoms, residues, rings, bonds, KD-trees, source clouds, neighbour finding, ring/bond discovery, spatial joins, or an aggregate model over the 720. The typed C++ model and resident indexes own those semantics. Python's only role is fitting over the lean emitted substrate (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:50`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:60`).

**Rationalization guard:** If a future task seems to require Python to "just rebuild enough model/spatial structure to do X", that is a design failure signal. The correct move is to add a C++ emitted scalar/sidecar, add a C++ named reduction, or extend the C++ accumulator. This applies to both memory branches. In branch (b), "buffers of the things we care about" are C++ buffers, not a Python aggregate (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:56`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:60`).

**Existing architecture supports this law:** `GUIDANCE.md` says emit work lives on the C++ model side and Python only consumes the substrate and fits (`GUIDANCE.md:22`, `GUIDANCE.md:34`). `PARTITION_FILTER_DESIGN.md` says Python receives emitted scalars/bins/sidecars/manifests only, while C++ owns source association, geometry, exclusions, reductions, isolation, bin IDs, and case selection (`PARTITION_FILTER_DESIGN.md:16`, `PARTITION_FILTER_DESIGN.md:23`). `NODE_STORE_CONTRACT_2026-06-02.md` names the law directly: no maths model 2 in Python; C++ owns model semantics (`NODE_STORE_CONTRACT_2026-06-02.md:24`, `NODE_STORE_CONTRACT_2026-06-02.md:39`).

## 2. Existing Built Substrate Pieces To Reuse

**BUILT: typed run spine.** The current `RunData` object owns the typed protein, conformation, manifest, DFT frames, and frame map (`RunData.h:124`, `RunData.h:141`). `DftFrameSet` is keyed by original frame index and returns per-atom DFT shielding by `(atom, originalIndex)` (`RunData.h:35`, `RunData.h:78`). `FrameMap` maps H5 row to original frame index and exposes sorted rows with DFT targets (`RunData.h:99`, `RunData.h:122`).

**BUILT but trajectory-only today:** `RunLoader::Load` currently requires a trajectory calcset and rejects non-trajectory manifests (`RunData.cpp:76`, `RunData.cpp:87`). It loads the protein/conformation through `QtProteinLoader::LoadRunPath`, loads FF14SB charges from topology, validates PBC-wrapped bond lengths, walks `.LGS` DFT frames through `DftShieldingLoader::LoadAndValidate`, and builds a `FrameMap` (`RunData.cpp:76`, `RunData.cpp:146`). Static ingest must generalize this load shape rather than copy its logic into a separate Python or ad hoc reader path.

**BUILT: single-pose model seam.** The reader already has `SingleConformation`: one pose, no H5, positions from the snapshot `Pos` column, `frameCount()==1`, `timePicoseconds()==0.0`, and no fake one-frame trajectory (`../model/SingleConformation.h:1`, `../model/SingleConformation.h:10`; `../model/SingleConformation.cpp:18`, `../model/SingleConformation.cpp:30`). The base `Conformation` API already supplies `frameCount`, `originalFrameIndex`, and `atomPosition` as the shared seam for trajectories and single poses (`../model/Conformation.h:54`, `../model/Conformation.h:80`).

**BUILT: per-protein single-pose load.** `QtProteinLoader::LoadPose` exists for `SinglePose`: it loads sidecar topology, then loads flat per-atom calculator NPYs from the pose directory with `FrameNpyLoader`, then constructs `SingleConformation` (`../io/QtProteinLoader.h:64`, `../io/QtProteinLoader.h:70`; `../io/QtProteinLoader.cpp:148`, `../io/QtProteinLoader.cpp:208`).

**CAVEATED:** `FrameNpyLoader` is intentionally lenient today: it enumerates `*.npy`, skips unrecognized or malformed arrays, and returns a partial snapshot unless no recognized NPYs load (`../io/FrameNpyLoader.h:17`, `../io/FrameNpyLoader.h:21`; `../io/FrameNpyLoader.cpp:46`, `../io/FrameNpyLoader.cpp:121`). That is acceptable for the viewer/current reader, but not acceptable for the 720 static ingest. The 720 static path must use a strict expected-field loader over the same NPY primitive readers: resolve documented per-protein expected paths exactly, validate required fields, and log-and-stop on any missing expected file or required array. It must not reuse the `FrameNpyLoader::LoadSnapshotDir` glob/enumerate behavior for required static fields.

**BUILT: topology sidecar.** `QtTopologySidecar` loads `extraction_manifest.json`, `atoms_category_info.npy`, `residues.npy`, `bonds.npy`, `rings.npy`, and `ring_membership.npy` into typed `QtAtom`, `QtResidue`, `QtBond`, `QtRing`, and `QtRingMembership` structures (`../io/QtTopologySidecar.h:1`, `../io/QtTopologySidecar.h:21`). The current implementation uses canonical sidecar names (`../io/QtTopologySidecar.cpp:374`, `../io/QtTopologySidecar.cpp:384`) and populates the topology and ring atom order from those typed rows (`../io/QtTopologySidecar.cpp:542`, `../io/QtTopologySidecar.cpp:580`).

**BUILT: NPY primitive.** `QtNpyReader` validates NPY magic/version, dtype descriptor, shape, row count, and byte count for structured arrays (`../io/QtNpyReader.h:41`, `../io/QtNpyReader.h:60`; `../io/QtNpyReader.h:131`, `../io/QtNpyReader.h:203`) and widens rank-1/rank-2 numeric NPYs to doubles while rejecting unsupported rank/dtypes and structured dtypes through the numeric path (`../io/QtNpyReader.cpp:296`, `../io/QtNpyReader.cpp:441`).

**BUILT: DFT target loader.** `.LGS` carries `dft.frames[].meta_json`; each frame points to a DFT job meta JSON (`../io/CalcsetManifest.h:39`, `../io/CalcsetManifest.h:72`; `../io/CalcsetManifest.cpp:459`, `../io/CalcsetManifest.cpp:526`). `DftShieldingLoader` reads `files.out_primary` strictly from that meta JSON, parses the ORCA output, validates atom count and dia+para identity, and returns a `DftShieldingFrame` (`../io/DftShieldingLoader.h:16`, `../io/DftShieldingLoader.h:27`; `../io/DftShieldingLoader.cpp:45`, `../io/DftShieldingLoader.cpp:110`). `OrcaShieldingParser` stores total/dia/para raw 3x3 matrices and ORCA-input coordinates for frame checks (`../io/OrcaShieldingParser.cpp:99`, `../io/OrcaShieldingParser.cpp:122`).

**BUILT: resident indexes.** `ResidentIndexes` packages `TypedAtomIndex`, `SpatialIndexSet`, `RingGeometryCache`, and `TemporalIndex` (`ResidentIndexes.h:14`, `ResidentIndexes.h:21`). `BuildResidentIndexes` builds typed atom lookup, spatial source clouds, ring geometry cache, and temporal index from `RunData` (`ResidentIndexes.cpp:7`, `ResidentIndexes.cpp:13`). `SpatialIndexSet` builds per-frame KD trees for atoms, bond midpoints, ring centers, charge sites, and all bond midpoints (`SpatialIndexSet.h:21`, `SpatialIndexSet.h:85`; `SpatialIndexSet.cpp:75`, `SpatialIndexSet.cpp:139`). `RingGeometryCache` stores canonical per-frame ring centers/normals (`RingGeometryCache.h:1`, `RingGeometryCache.h:27`; `RingGeometryCache.cpp:26`, `RingGeometryCache.cpp:40`). `TypedAtomIndex` gives typed, set-returning atom lookup without name/position heuristics (`TypedAtomIndex.h:1`, `TypedAtomIndex.h:47`).

**BUILT: C++ verb/traversal contract.** `Verbs` expose position, KD-neighbour query, catalog reads, ring geometry, typed atom lookup, heavy-parent logic, and shared displacement geometry over the resident `Body`; they do not rebuild state (`Verbs.h:1`, `Verbs.h:30`). `RelationshipEngine::RunTraversal` already walks the `(atom, frame)` index space and lets carriers own record/reduction/sink shape (`RelationshipEngine.h:63`, `RelationshipEngine.h:91`; `RelationshipEngine.cpp:19`, `RelationshipEngine.cpp:79`).

**BUILT: current per-atom aggregate substrate.** `PerAtomSubstrate` is a lean all-atom aggregate emitter that writes one row per `(DFT frame slot, atom)` and deliberately does not emit a default source CSV (`PerAtomSubstrate.h:1`, `PerAtomSubstrate.h:5`). It has bounded query defaults, not resident pair dumps (`PerAtomSubstrate.h:25`, `PerAtomSubstrate.h:39`). It computes pair contributions, reductions, driver magnitudes, isolation, dominance, bins, support counts, sidecars, and row alignment in C++ (`PerAtomSubstrate.cpp:895`, `PerAtomSubstrate.cpp:915`; `PerAtomSubstrate.cpp:1987`, `PerAtomSubstrate.cpp:2104`; `PerAtomSubstrate.cpp:3300`, `PerAtomSubstrate.cpp:3445`).

## 3. Static Run Shape

**DESIGNED:** Add a `StaticPoseConformation` ingest path in rediscover. It should sit at the same model level as existing `SingleConformation`, not as a Python reader and not as a fake trajectory. It can reuse `SingleConformation` internally if that remains sufficient, but rediscover should name the semantic run shape explicitly because the 720 cohort has a different contract: one WT absolute static pose per protein, keyed by protein ID, with per-protein producer output and DFT target.

Minimum static run carrier:

```text
StaticPoseRunData
  protein_id
  curated_set_id
  source_producer_manifest
  unique_ptr<QtProtein>
  unique_ptr<Conformation>   # SingleConformation or StaticPoseConformation
  StaticFrameMap             # exactly one frame slot: static_frame_index=0
  DftFrameSet or StaticDftTargetSet
  Catalog or StaticCatalog
  ResidentIndexes
```

`RunData` can be generalized if the code owner prefers a single struct, but the static design must remove the trajectory-only assumptions:

- No requirement that `manifest.kind == Trajectory`; current code rejects non-trajectory at `RunData.cpp:86` and must not remain the static path gate.
- No `h5()` requirement for static arrays; current `Catalog` reads many time-series arrays through `run.h5()` (`Catalog.cpp:76`, `Catalog.cpp:210`), while a static pose's producer arrays live in the snapshot. Static catalog reads must work from the one resident snapshot.
- No `time_ps` meaning beyond provenance; static rows should carry `static_frame_index=0` and optional `source_frame_index/source_time_ps` if the producer records a source trajectory origin.

**DESIGNED:** Static pose is effectively one frame. Existing per-(atom, frame) machinery carries over as per-(protein, atom) by setting the frame axis length to 1 and retaining the same source-discovery/reduction semantics. The brief requires this exact shape (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:62`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:73`).

## 4. Input Path Discipline

**LAW:** No file discovery in the design. The curated set path is TBD. Execution receives a manifest or manifest-root filled in later by the lead. The static ingest resolves every expected per-protein file by the documented protein-ID convention and stops on the first missing required file. It must not glob for proteins, try alternate names, or infer co-location (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:29`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:42`; `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:115`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:121`).

**LAW:** Required NPYs are resolved by exact path, never by enumeration. The static loader constructs each required NPY path from the lead-vetted field-to-filename convention and opens that path directly. It must not call `entryList`, glob `*.npy`, enumerate the pose directory, or use `FrameNpyLoader::LoadSnapshotDir` for required static fields. A missing required NPY path is log-and-stop, not a partial snapshot or NaN-filled success.

**DESIGNED:** The execution-time input should be a curated cohort manifest, not a directory scan:

```json
{
  "cohort_id": "curated_720_wt_static",
  "schema_version": 1,
  "proteins": [
    {
      "protein_id": "TBD_PROTEIN_ID",
      "producer_lgs": "resolved/by/lead/later/TBD_PROTEIN_ID.LGS"
    }
  ]
}
```

The spec is path-agnostic. It does not target the stale 725 duplicate set. It does not assume NPYs are co-located. The cohort manifest is an execution fill, not part of this docs-only pass.

**DESIGNED:** At execution, the loader should accept exact `.LGS` paths from the cohort manifest. Although `CalcsetManifest::Load` can resolve a directory by listing a single `*.LGS` (`../io/CalcsetManifest.h:137`, `../io/CalcsetManifest.h:154`; `../io/CalcsetManifest.cpp:147`, `../io/CalcsetManifest.cpp:178`), the 720 static ingest should call it with exact `.LGS` paths constructed from the cohort manifest/protein ID. That avoids directory enumeration in the 720 path.

**DESIGNED:** Required per-protein checks:

- Manifest exists and `protein_id` equals cohort entry.
- Active run kind is `single_pose` or a new static pose kind accepted by the reader.
- `single_pose.pose_dir` and `single_pose.extraction_manifest` exist (`../io/CalcsetManifest.cpp:384`, `../io/CalcsetManifest.cpp:422`).
- Topology sidecar required files exist and decode: `atoms_category_info.npy`, `residues.npy`, plus bonds/rings/ring membership required for rediscover static work even if the viewer currently treats some as soft fail.
- Position NPY exists and has atom-axis rows equal to `QtProtein::atomCount`.
- Required kernel shadow NPYs exist for the selected static schema and have expected axis lengths.
- DFT raw target exists, is keyed by protein ID/static frame, and validates against topology atom count.

No missing required field becomes NaN silently in the 720 static ingest. Optional fields can be absent only if the schema marks them optional and the run manifest records absent support.

## 5. Kernel Sourcing, Resolved

**RESOLVED:** Kernel sourcing is not open. The ingest consumes the given per-protein producer output from the 720 static-pose `nmr_extract` re-run with MOPAC ON. It reuses the existing r2SCAN DFT shielding target, does not re-implement `nmr_extract`, does not recompute kernels reader-side, and does not generate new ORCA/DFT calculations (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:75`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:79`).

**DESIGNED: source table.**

| Static payload | Source | Existing reader grounding | Static ingest rule |
|---|---|---|---|
| Protein identity | Cohort manifest + per-protein `.LGS` | `.LGS` holds dataset/protein identity (`../io/CalcsetManifest.h:81`, `../io/CalcsetManifest.h:132`) | Must match cohort entry; mismatch stops. |
| Topology | `extraction_manifest.json`, `atoms_category_info.npy`, `residues.npy`, `bonds.npy`, `rings.npy`, `ring_membership.npy` | `QtTopologySidecar` decodes these into typed model objects (`../io/QtTopologySidecar.h:1`, `../io/QtTopologySidecar.h:21`) | Required for static rediscover; no soft topology loss. |
| Positions | Static producer NPY `pos.npy` / catalog `FieldKind::Pos` | Field catalog marks `pos` atom-axis 3-vector (`../io/QtFieldCatalog.gen.h:191`, `../io/QtFieldCatalog.gen.h:193`); `SingleConformation` reads `Pos` for positions (`../model/SingleConformation.cpp:18`, `../model/SingleConformation.cpp:26`) | Required, one static frame. |
| Ring shadows | Biot-Savart, HM, ring-chi, pi-quadrupole, dispersion, ring geometry/counts | Field catalog lists ring-current and ring-related fields (`../io/QtFieldCatalog.gen.h:196`, `../io/QtFieldCatalog.gen.h:212`) | Read from producer NPY surface; no recompute except C++ neighbourhood geometry/reduction. |
| Bond/McConnell shadows | McConnell and MOPAC-McConnell fields | Field catalog lists McConnell fields (`../io/QtFieldCatalog.gen.h:213`, `../io/QtFieldCatalog.gen.h:215`) and MOPAC Mc fields (`../io/QtFieldCatalog.gen.h:249`, `../io/QtFieldCatalog.gen.h:251`) | Read producer arrays; C++ computes neighbourhood source geometry from topology/positions. |
| Electrostatic/APBS/AIMNet2/MOPAC shadows | APBS E/EFG, AIMNet2 charge/CRG/embedding/EFG, MOPAC charges/scalars/bond orders/Coulomb/Mc | Catalog has these groups (`../io/QtFieldCatalog.gen.h:240`, `../io/QtFieldCatalog.gen.h:253`; `../io/QtFieldCatalog.gen.h:95`, `../io/QtFieldCatalog.gen.h:101`) | Required/optional by lead-vetted schema; MOPAC ON means MOPAC arrays expected. |
| DFT target | Existing raw r2SCAN ORCA dia/para/total tensors, WT absolute sigma | DFT loader reads `.out` via meta JSON (`../io/DftShieldingLoader.cpp:45`, `../io/DftShieldingLoader.cpp:77`); parser stores raw matrices (`../io/OrcaShieldingParser.cpp:99`, `../io/OrcaShieldingParser.cpp:110`) | Required target. Reuse existing DFT; no new ORCA. Use raw total as primary target; dia/para diagnostic splits. |

**DESIGNED:** Static target should use absolute WT sigma, not deltas. If `orca_total/orca_diamagnetic/orca_paramagnetic` NPYs are present in the producer surface, they are provenance/cross-check shadows. The canonical target for fitting should still be the existing raw dia/para/total tensors loaded into the typed DFT target path unless the lead explicitly vets the NPY target arrays as the canonical source. The design does not invent a new schema for that choice; it surfaces the fork for lead vet. Static T2 carries the same frame-alignment caveat as any raw tensor comparison: the manifest should record the single-pose DFT/producer coordinate alignment diagnostic, while T0 remains rotation-invariant.

**DESIGNED NEGATIVE:** Do not add mutation-delta fields as a default static substrate family. The 720 static design is WT absolute sigma, not a mutation-pair delta model. Delta surfaces in the producer catalog may be real, but adding them here would change the scientific question instead of rescuing a relationship inside the chosen static-ingest design (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:51`).

**DESIGNED BUILD REQUIREMENT - T2 orientation guard:** For between-protein maths, `T0` is rotation-invariant and is poolable across proteins. Raw lab-frame `T2` components are not. The build must not run a between-protein component fit over lab-frame `T2` from arbitrarily oriented static proteins. Before any such component fit, the C++ emit must provide one of: stable local-frame target/feature T2 for atom classes with valid frames, with `frame_valid`, `frame_variant`, and `frame_anchor_atom_index`; equivariant/node-store consumption for the affected family; or invariant summaries such as `|T2|`, tensor norms, and tensor inner products. Python must not rotate tensors, rebuild local frames, or derive substitute T2-frame features (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:11`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:29`).

## 6. Static Model And Index Build

**DESIGNED:** For each protein:

1. Load manifest by exact `.LGS` path from the cohort manifest.
2. Load topology into `QtProtein` via `QtTopologySidecar`.
3. Strict-load static snapshot arrays into `QtConformationSnapshot` by exact required NPY paths.
4. Construct `StaticPoseConformation` or `SingleConformation` with one frame.
5. Load raw DFT target into a one-frame target set keyed by static frame 0/protein ID.
6. Build `Catalog`/`StaticCatalog`.
7. Build `ResidentIndexes`.
8. Traverse all target atoms for one frame and compute static row reductions.

**DESIGNED:** `ResidentIndexes` is still per-protein for branch (b), and per-protein plus optional cross-protein metadata for branch (a). Neighbour discovery remains local to a protein unless the lead's between-axis maths explicitly needs cross-protein normalization/reduction over row summaries. The 720 proteins are not one physical spatial system; there is no cross-protein KD tree for geometry.

**DESIGNED:** Source neighbourhoods are C++ model queries:

- Ring centers from `RingGeometryCache` over topology/positions.
- Bond midpoint clouds from `SpatialIndexSet`.
- Charge sites from typed atom positions/charges.
- Typed strata from `QtAtom`/`TypedAtomIndex`.
- Displacements, `r`, `inv_r3`, `cos_theta`, and dipolar forms from `verbs::displacement`.

Existing code grounds this shape: `SpatialIndexSet` builds the source clouds per frame (`SpatialIndexSet.cpp:83`, `SpatialIndexSet.cpp:138`), `verbs::near` queries KD trees (`Verbs.cpp:34`, `Verbs.cpp:42`), and `verbs::displacement` is the single shared near-field geometry definition (`Verbs.h:134`, `Verbs.h:149`; `Verbs.cpp:102`, `Verbs.cpp:121`).

## 7. Accumulation And Memory Strategy

The memory strategy is OPEN and must follow the between-axis statistics maths. Do not bake either branch. The implementation should put a common C++ interface in front of both:

```text
StaticCohortAccumulator
  begin_cohort(cohort_manifest)
  observe_protein(StaticPoseRunData, ResidentIndexes, StaticCatalog)
  observe_row(StaticRowDigest)
  finalize()
  emit(StaticEmitWriter)
```

The output contract is the same for both branches. Only when rows become final differs.

### Branch (a): all-720-resident, emit once

**OPEN:** Hold the entire curated set in C++ resident model objects and indexes, then emit the lean digest at the end. This enables between-protein statistics that need all proteins co-resident, such as global centering, cohort-level stratified thresholds, or statistics whose definition depends on the full row distribution.

**C++ memory shape:** `vector<StaticPoseRunData>`, `vector<ResidentIndexes>`, and C++ accumulators over row summaries. Python sees nothing until the final lean emit.

**Pros:** Maximal flexibility for not-yet-settled between-axis maths. It can compute cohort-level ranks/quantiles/normalizations exactly in one pass over resident rows.

**Cost:** Higher RAM. The brief says this is acceptable: 128 GB, swap OK (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:80`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:88`). Choosing all-resident "to be safe" is still a between-axis maths commitment if it bakes global centering, cohort quantiles, or population-level reductions; it is not the default until the lead chooses the estimator.

### Branch (b): sequential, C++ running buffers

**OPEN:** Load one protein, build typed model/indexes, update only the reduced quantities the between-axis maths requires, release the protein, and emit at the end.

**C++ memory shape:** one `StaticPoseRunData` at a time, plus C++ buffers such as Welford cells, stratified sufficient statistics, row-key arrays, optional compact row digests, and lead-vetted reducers. These are not Python data structures and not a second protein model.

**Pros:** Lower RAM. This is better if the between-axis method can be expressed as streaming sufficient statistics plus a final pass/emit.

**Cost:** Some statistics may need two passes or compact C++ row spooling. If the maths needs all row vectors co-resident, branch (b) may not be valid. That decision belongs to the lead/statistics discussion, not this spec.

### Non-foreclosing recommendation

**RECOMMENDED DESIGN:** Build one `StaticCohortAccumulator` API with two storage backends:

- `AllResidentStaticAccumulator`
- `StreamingStaticAccumulator`

Both consume the same `StaticRowDigest` produced from the typed C++ model. Both produce the same `StaticEmitWriter` output. The CLI/config records `memory_strategy=all_resident|streaming`, but the default remains unset until the lead chooses after the between-axis maths is settled.

This structure honors the foreclosure: branch (a) keeps proteins/indexes in C++; branch (b) keeps only C++ reduced buffers. Neither branch hands Python a protein model, pairwise dump, or aggregate model.

## 8. Static Emit Contract

**DESIGNED:** The static substrate is per `(protein_id, atom_index)` with exactly one static frame per protein. It should be row-aligned NPY sidecars plus a compact CSV or JSONL identity table, matching the current row-alignment discipline. The current per-atom substrate records row alignment as `row_id == frame_slot * n_atoms + atom_index; sidecar row i == rows.csv row_id i` (`PerAtomSubstrate.cpp:3300`, `PerAtomSubstrate.cpp:3308`). Static should use:

```text
row_id == cumulative_atom_offset + atom_index_within_protein
join key == (protein_id, atom_index)
sidecar row i == static_rows.csv row_id i
```

**DESIGNED:** The lean static menu is review-folded and lead-vettable as requirements pending build. It should not blindly copy every current trajectory sidecar. Current `PerAtomSubstrate` writes many target/feature sidecars (`PerAtomSubstrate.cpp:2640`, `PerAtomSubstrate.cpp:2692`) and lists them in `PerAtomSubstrateSidecars` (`PerAtomSubstrate.cpp:3709`, `PerAtomSubstrate.cpp:3734`). For 720 static, the schema should be smaller and explicit, with the required grouped sidecars below and exact names/formats left for lead vet.

### Carry Over From Per-(atom, frame)

Carry these families, subject to lead vet:

- Identity: atom index, residue index/number, atom labels, element, atom role, stratum, topology ordinals, BMRB/IUPAC label projection. Existing row identity fields are emitted today (`PerAtomSubstrate.cpp:2870`, `PerAtomSubstrate.cpp:2890`).
- DFT target: total raw/decomposed T0/T2 as primary; dia/para T0/T2 as diagnostics. Current target build stores total/dia/para raw and decompositions (`ExtractionSupport.cpp:45`, `ExtractionSupport.cpp:68`), and current sidecars include total/dia/para T0/T1/T2 (`PerAtomSubstrate.cpp:2657`, `PerAtomSubstrate.cpp:2674`).
- Local frame fields for target atom classes that use them. Existing identity columns include frame axes/validity (`ExtractionSupport.cpp:70`, `ExtractionSupport.cpp:89`).
- Classical/reduced source features computed C++-side: ring JB, charge q/r^3, McConnell, MOPAC Coulomb shielding, MOPAC Mc shielding, APBS E/EFG, AIMNet2 charge/CRG/embedding presence, hbond/pq/disp/HM/ringchi/water/SASA/EEQ families as lead-vetted. Current direct features read these through `Catalog`/snapshot fields (`PerAtomSubstrate.cpp:1435`, `PerAtomSubstrate.cpp:1680`). The review-folded required structure is listed in the grouped-sidecar section below.
- C++ conditioning and isolation: nearest/gap counts, dominant fractions, source counts, self/same-residue exclusions, thin-support flags. Existing reductions and dominance live C++-side (`PerAtomSubstrate.cpp:957`, `PerAtomSubstrate.cpp:1023`; `PerAtomSubstrate.cpp:2040`, `PerAtomSubstrate.cpp:2104`).
- Separable row-aligned AIMNet2 256-d embedding sidecar with `embedding_present` and support flags. Current code supports a separable 256-d embedding sidecar (`PerAtomSubstrate.cpp:2689`, `PerAtomSubstrate.cpp:2692`), and the node-store contract already flags embedding as separable from the classical substrate (`NODE_STORE_CONTRACT_2026-06-02.md:150`, `NODE_STORE_CONTRACT_2026-06-02.md:158`).

### Required Review-Folded Grouped Sidecars

**DESIGNED:** The adversarial review identifies relationship families that the lean menu must not collapse away (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:23`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:61`). These are C++-emitted scalar or grouped sidecars only: never Python-derived features, never a Python protein model, and never a source/pair dump.

- Ring valid/self split: emit grouped `ring_jb_valid_T0`, `ring_jb_valid_T2`, `ring_jb_self_or_bonded_T0`, and `ring_jb_self_or_bonded_T2`, plus counts, using the existing `ringSourceIsSelfOrBonded` rule. This preserves through-space ring current separately from own-ring/bonded aromatic shielding (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:33`).
- Per-type/per-category mechanism structure: emit grouped `ring_{bs,hm,pq,disp}_per_type_T2`, `mc_category_T2`, and `mopac_mc_category_T2`, plus category counts/support. Do not collapse these to one total per mechanism, and do not expand them to per-source tables (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:35`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:37`).
- Separated electrostatic slabs: emit APBS `E` and `EFG_T2`; MOPAC Coulomb `E`; raw EFG backbone/aromatic and shielding T2; AIMNet2 charge/CRG/EFG splits; FF14SB, MOPAC, and EEQ charge scalars; and per-source presence/support flags. These slabs are not interchangeable and must not be collapsed into one field magnitude (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:39`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:41`).
- Compact hbond/solvation/exposure sidecar: emit geometric `hbond_T2`, Larsen per-class T2, water `E`/`EFG` total and first-shell terms, DSSP hbond energy/count/secondary-structure flags, SASA, SASA normal, and hydration-shell scalars (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:43`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:45`).
- AIMNet2 electronic-environment sidecar: keep `static_optional_aimnet2_embedding.npy` as a separable f32 256-d sidecar with `embedding_present` and support flags. Do not derive embedding features in Python (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:47`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:49`).

### Drop As Trajectory-only

Drop these from the default static substrate:

- `h5_row`, trajectory `frame_slot`, `original_frame_index`, and `time_ps` as primary axes. Replace with static provenance fields. The 1P9J trajectory fields may remain only in oracle-parity diagnostics.
- Welford/frame-modulation sidecars such as `per_atom_substrate_driver_modulation_by_atom.npy`. Current code computes modulation as SD across frames (`PerAtomSubstrate.cpp:2842`, `PerAtomSubstrate.cpp:2855`), which is trajectory-axis information and not meaningful for one static pose.
- Temporal windows and frame-ablation fields. `TemporalIndex` can exist in `ResidentIndexes`, but static frame count is one, so no between-frame modulation is emitted.
- Default named pair query outputs. Current pair queries are bounded and query-oriented (`PerAtomSubstrate.cpp:3479`, `PerAtomSubstrate.cpp:3639`), but 720 static must default them off. Any future pair/materialization output must be a named, capped, C++-produced diagnostic, never the fitter substrate.
- Any per-source flat table as a default emit. Source-level arrays are allowed only as small named-query diagnostics for an explicit gate, not as 720 fitter input. The cap rule is top-k sources over a small named atom/protein subset; never all-pairs x all-atoms, never all proteins by default, and never an input to the fitter.

### New Static Fields

Add these static provenance/key fields:

- `protein_id` and optional numeric `protein_ord`.
- `cohort_id`, `curated_set_id`, `producer_run_id`, `producer_manifest_id`.
- `static_frame_index` fixed to 0.
- `source_lgs_path_hash` or manifest digest, not necessarily full filesystem paths in the fitter substrate.
- `dft_source`: e.g. `raw_orca_out` or lead-vetted `producer_orca_npy`.
- `mopac_on` and required-kernel presence flags.
- `atom_count_in_protein`, `row_offset`, `row_id`.
- `support_thin_flag` per element/stratum/category, derived from cohort support counts.

**CAVEATED:** The final column names/formats need lead vet before code. The brief explicitly says to surface the format for lead vet and not invent a final schema (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:96`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:103`).

## 9. Manifest And Support Contract

**DESIGNED:** The output directory should contain:

```text
static_rows.csv or static_rows.jsonl
static_target_T0.npy
static_target_T2.npy
static_target_dia_T0.npy
static_target_dia_T2.npy
static_target_para_T0.npy
static_target_para_T2.npy
static_features_classical.npy
static_features_ring_split.npy
static_features_ring_per_type_T2.npy
static_features_mc_category_T2.npy
static_features_electrostatic_slabs.npy
static_features_hbond_solvation_exposure.npy
static_features_conditioning.npy
static_features_dominance.npy
static_optional_aimnet2_embedding.npy  # separable f32 256-d sidecar with embedding_present/support flags
static_aimnet2_embedding_support.npy
static_column_specs.json
static_manifest.json
static_protein_manifest.csv/jsonl
```

**DESIGNED:** `static_manifest.json` must record:

- cohort ID, row count, protein count, total atom count;
- memory strategy used;
- exact schema version and column spec digest;
- row alignment contract;
- all cutoffs and exclusion rules;
- DFT target source and frame-alignment status, including the single-pose Kabsch/coordinate diagnostic when T2 is compared across producer and DFT frames;
- per-feature support counts;
- absent optional slabs;
- disk budget calculation and pass/fail status;
- oracle-parity gate result if this is a gate run.

The current manifest pattern records row alignment, atom/frame counts, cutoffs, DFT frame alignment, sidecars, support rows, and absent slabs (`PerAtomSubstrate.cpp:3300`, `PerAtomSubstrate.cpp:3445`). Static should reuse that manifest style with protein-axis fields.

**DESIGNED:** `static_protein_manifest` must have one row per protein:

```text
protein_id,status,n_atoms,n_residues,n_rings,n_bonds,
dft_present_atoms,missing_required_field,missing_optional_fields,
thin_strata_json,source_manifest_digest,notes
```

Status values:

- `ok`: all required files/fields present and atom counts align.
- `skipped_optional_missing`: optional fields absent but row emitted with support flags.
- `failed_missing_required`: no substrate rows emitted for the cohort run; log-and-stop.
- `failed_validation`: atom count/schema/dtype/target validation failed; log-and-stop.

**DESIGNED:** Thin support flags must be C++-emitted. Do not ask Python to group proteins/atoms and determine thin element strata by reconstructing topology. Current support counts are C++ stats emitted into manifests (`PerAtomSubstrate.cpp:3736`, `PerAtomSubstrate.cpp:3767`) and the partition docs require support-flagging without overclaiming thin categories (`STATE.md:202`, `STATE.md:210`).

## 10. Disk Budget

The disk budget is the same under branch (a) and branch (b), because both emit the same lean substrate. Only C++ RAM strategy differs.

Let:

```text
P = 720 proteins
A_avg = average emitted atoms per protein
R = total rows = sum_p atom_count_p = P * A_avg
C64 = float64 columns across all numeric sidecars
C32 = float32 columns across separable embedding sidecar
Cint = int32/bin/provenance numeric columns
```

Recommended lead-vet lean defaults:

```text
C64 <= 220      # ceiling for a real static subset, not permission to carry every trajectory family
C32 = 256       # AIMNet2 embedding is separable and carries embedding_present/support flags
Cint <= 48      # bins/provenance ordinals if sidecarized
identity/provenance CSV/JSONL <= 256 bytes/row compressed or <= 512 bytes/row uncompressed
```

Budget with `A_avg = 1,000`:

```text
R = 720,000
float64: 720,000 * 220 * 8 = 1.27 GB
int32:   720,000 * 48  * 4 = 0.14 GB
identity uncompressed worst: 720,000 * 512 = 0.37 GB
subtotal before embedding sidecar: about 1.8 GB plus manifests
embedding sidecar: 720,000 * 256 * 4 = 0.74 GB
total with embedding sidecar: about 2.6 GB plus manifests
```

Budget with a conservative `A_avg = 2,000`:

```text
R = 1,440,000
float64: 1,440,000 * 220 * 8 = 2.53 GB
int32:   1,440,000 * 48  * 4 = 0.28 GB
identity uncompressed worst: 1,440,000 * 512 = 0.74 GB
subtotal before embedding sidecar: about 3.6 GB plus manifests
embedding sidecar: 1,440,000 * 256 * 4 = 1.47 GB
total with embedding sidecar: about 5.1 GB plus manifests
```

Budget with an intentionally high `A_avg = 4,000`:

```text
R = 2,880,000
float64: 2,880,000 * 220 * 8 = 5.07 GB
int32:   2,880,000 * 48  * 4 = 0.55 GB
identity uncompressed worst: 2,880,000 * 512 = 1.47 GB
subtotal before embedding sidecar: about 7.1 GB plus manifests
embedding sidecar: 2,880,000 * 256 * 4 = 2.95 GB
total with embedding sidecar: about 10.1 GB plus manifests
```

**DESIGNED gate:** Before writing rows, the static emitter computes the projected output size from manifest atom counts and enabled sidecars. It refuses the run if projected total output exceeds 15 GB. The current per-atom emitter has an explicit size gate for appended flat payload (`PerAtomSubstrate.cpp:3777`, `PerAtomSubstrate.cpp:3790`); static should use the same fail-loud pattern, but the static gate must sum all enabled sidecars, identity rows, manifests, and the separable AIMNet2 embedding sidecar, not only an append slab or column subset.

**FORECLOSED:** No per-source/pairwise default emit. A pairwise table would multiply rows by source count and violates the <15 GB substrate law. `NODE_STORE_CONTRACT_2026-06-02.md` allows materialization only as a named transient tool, never the default interface (`NODE_STORE_CONTRACT_2026-06-02.md:95`, `NODE_STORE_CONTRACT_2026-06-02.md:102`). The 720 static ingest should not create even a temporary Python-readable pairwise dump for fitting.

**CAVEATED:** The exact curated 720 atom counts are not grounded in this docs-only pass because the path is TBD. The design therefore includes a preflight size gate and reports exact `R` before execution.

## 11. Oracle-Parity Gate

**DESIGNED:** Before trusting any 720 output, run a known 1P9J WT static pose through the static path and compare it to the corresponding 1P9J trajectory-frame substrate.

### Gate Input

- One 1P9J trajectory frame with DFT target and producer per-frame NPYs.
- A static-pose producer output generated from the exact extracted positions of a named trajectory frame, the same topology, and the same existing raw ORCA target.
- MOPAC ON and required static fields present.

### Join

Compare after joining:

```text
static key:     (protein_id, atom_index)
trajectory key: (protein_id, selected_h5_row, atom_index)
```

The static row ID need not equal the trajectory row ID. The gate ignores row-order differences after the join, but static production should still be deterministic.

### Must Match

Exact/integer equality:

- atom count;
- atom indices;
- residue indices/numbers;
- element ordinals;
- role/stratum ordinals;
- topology-derived bond/ring IDs and source IDs;
- source counts and presence flags;
- optional field support flags when the same producer fields are present.

Floating equality:

- input positions: f64 round-trip / bit-identical where both paths read the same coordinate bytes, otherwise `abs <= 1e-12 Angstrom`; the `1e-9 Angstrom` budget is for derived geometry, not the oracle input positions;
- local frame axes: `abs <= 1e-10` per component, `frame_valid` identical;
- ring centers/normals and bond midpoint/axis geometry: `abs <= 1e-9 Angstrom` for positions, `abs <= 1e-10` for unit vectors/dimensionless terms;
- neighbour-derived displacement, `r`, `inv_r3`, `cos_theta`, `dipolar`: `abs <= 1e-9` for Angstrom-scaled values and `abs <= 1e-10` for dimensionless values unless f32-origin producer arrays force a looser threshold;
- copied f64 producer kernels after reader widening: `abs <= 1e-12` or bit-identical if both paths read the same NPY bytes;
- f32-origin embeddings/fields widened to double: `abs <= 5e-7`;
- DFT raw total/dia/para tensors and decomposed total/dia/para T0/T2: `abs <= 1e-9 ppm` when both paths parse the same raw ORCA text; if one path is lead-vetted NPY rounded output, `abs <= 1e-4 ppm` and dia+para identity must remain within the existing DFT loader tolerance of 0.1 ppm (`../io/DftShieldingLoader.cpp:39`, `../io/DftShieldingLoader.cpp:104`).

Feature-family parity:

- V1 oracle parity gates ring-current and McConnell/MOPAC-Mc reduced families only, because that is the current comparator coverage template.
- Charge, MOPAC Coulomb, APBS, AIMNet2, hbond/pq/disp/HM/ringchi/water/SASA/EEQ parity are deferred until the C++ parity report is explicitly extended for those families.
- Conditioning/isolation/dominance fields attached to those v1 ring/McConnell families.
- Column spec metadata for units/irreps/mechanism where the same v1 column exists.

### May Differ

Allowed differences:

- `row_id`, because static uses cumulative protein atom offsets while trajectory uses frame slot times atom count.
- `h5_row`, `frame_slot`, `original_frame_index`, and `time_ps`; static should replace them with `static_frame_index=0` plus optional source-frame provenance.
- trajectory-only Welford/frame-modulation columns, because static has no time axis.
- run-kind/provenance fields such as `run_kind`, `cohort_id`, `producer_manifest_digest`, `source_lgs_path_hash`.
- optional fields absent in one producer surface, but only if the gate is explicitly testing optional absence. The normal oracle gate should use a fully populated static output and fail on missing required fields.

### Gate Result

The gate produces:

- pass/fail summary;
- row count and atom count;
- max absolute difference per v1 numeric family, and explicit "not gated in v1" labels for deferred families;
- first mismatching key/column if any;
- manifest of allowed differences.

No Python spatial/model reconstruction is allowed in the gate. The comparator may be Python only over the lean emitted rows/sidecars. It must not open `trajectory.h5`, ORCA outputs, topology sidecars, or producer NPY directories. Current discipline already says analysis scripts must not open `trajectory.h5`, ORCA outputs, old run dirs, or per-source dumps to reconstruct values (`PARTITION_FILTER_DESIGN.md:16`, `PARTITION_FILTER_DESIGN.md:23`), and the brief repeats "never open trajectory.h5 in Python" (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:116`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:119`).

## 12. Execution Gates For The Future Build

This docs pass runs none of these gates. They are designed gates for a future implementation loop.

**Preflight gate:**

- Read cohort manifest by exact path.
- Validate exactly 720 unique protein IDs.
- Validate no ID duplicates.
- Validate every `.LGS` path exists and matches protein ID.
- Validate projected disk size < 15 GB.
- Validate memory strategy explicitly selected or lead-defaulted.

**Per-protein load gate:**

- Strict topology/snapshot/target load.
- Atom counts match across topology, positions, required atom-axis fields, and DFT target.
- Required MOPAC/APBS/AIMNet2/static fields present per schema.
- No file discovery, no fallback file names, no pose-directory NPY enumeration.

**Emit gate:**

- Rows emitted equals sum atom counts.
- Sidecar shapes match row count and column specs.
- Manifest support counts match emitted presence flags.
- No default source/pairwise materialization exists in output.

**Oracle gate:**

- 1P9J static-vs-trajectory ring+McConnell v1 equivalence passes before 720 numbers are accepted.

**Python gate:**

- Fitter consumes only `static_rows.*`, `static_*.npy`, `static_column_specs.json`, and `static_manifest.json`.
- Fitter does not import reader model code, open `.LGS`, open `trajectory.h5`, open ORCA outputs, open topology sidecars, scan producer directories, rotate T2 tensors, rebuild local frames, or derive substitute grouped sidecars.
- Frozen `get_C` remains the basis bridge on the Python/e3nn side; no Python recomputation of DFT tensor decomposition or source geometry.

## 13. Open / Undecided Running List

**OPEN - memory strategy.** All-resident vs sequential C++ buffers remains maths-dependent and must not be selected in this spec (`CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:80`, `CODEX_BRIEF_SPEC_720_INGEST_2026-06-04.md:88`).

**OPEN - final static schema.** The lean column menu above is no longer a bare proposal: after the adversarial review, the T2 orientation guard and grouped sidecar families in section 8 are lead-vetted requirements pending build (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:11`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:61`). What remains open is exact file/column naming, dtype/packing, compression, T1 diagnostics, and any additional optional diagnostics beyond those required families. Python-derived rotations/features and source/pair dumps remain foreclosed.

**OPEN - T2 orientation representation path under a DESIGNED guard.** The build requirement is fixed: do not pool lab-frame `T2` components across arbitrarily oriented static proteins. The open choice is which C++ path the schema uses before any between-protein component fit: stable local-frame T2 where frames are valid, equivariant/node-store consumption, or invariant summaries such as `|T2|`, norms, and tensor inner products. Python must not rotate tensors or rebuild frames (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:11`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:29`).

**OPEN - grouped sidecar packing/names under required content.** Required content is the ring valid/self split; ring per-type T2 and Mc/MOPAC-Mc category T2; separated APBS/MOPAC/AIMNet2/FF14SB/EEQ electrostatic slabs; compact hbond/solvation/exposure slabs; and the separable AIMNet2 256-d embedding sidecar with support flags. What remains open is packing and exact column names, not whether these relationship families are carried (`REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:33`, `REVIEW_codex_adversarial_SPEC_720_2026-06-04.md:49`).

**OPEN - canonical DFT target source if both raw ORCA and producer ORCA NPYs exist.** The design recommends existing raw ORCA via `DftShieldingLoader` as the canonical target and NPY target arrays as cross-check/provenance unless the lead vets otherwise. No new ORCA is part of this design.

**OPEN - exact cohort manifest shape.** The input path is TBD. The design requires a manifest of exact protein IDs and `.LGS` paths, but the final JSON keys can follow the lead's execution convention.

**OPEN - between-axis statistics method.** The accumulator must expose sufficient hooks for the method, but this spec does not choose centering, normalization, null model, or stratified thresholds.

**OPEN - full-menu oracle parity coverage.** V1 gates ring+McConnell/MOPAC-Mc only. Full-menu parity requires an explicit C++ parity-report extension before it can be treated as a gate.

**CAVEATED - static T2 frame alignment.** The one-pose manifest must carry the DFT/producer coordinate alignment diagnostic before T2 comparisons are trusted beyond T0/rotation-invariant checks. Coordinate alignment alone does not authorize pooled lab-frame T2 component fitting across proteins; the T2 orientation guard above still applies.

**OPEN - static output compression.** Row identity can be CSV/JSONL or compact binary/Arrow-like if the lead wants. The default spec is row-aligned CSV plus NPY sidecars because that matches the existing substrate style and keeps Python as a tabular fitter.

**CAVEATED - exact disk/RAM counts.** Actual atom counts for curated 720 were not grounded because the path is TBD and was intentionally not hunted. The design includes preflight arithmetic and fail-loud gates.

**CAVEATED - current loader leniency.** `FrameNpyLoader` currently enumerates and permits partial snapshots; the future static implementation must not reuse that leniency or its pose-dir glob for 720 required fields (`../io/FrameNpyLoader.h:17`, `../io/FrameNpyLoader.h:21`).

**CAVEATED - no claim of built behavior.** `StaticPoseConformation`, strict static snapshot loading, cohort accumulator, static manifest, and 1P9J static oracle comparator are designed here but not implemented in this pass.

## 14. Final Design In One Loop

The 720 static ingest is a C++ rediscover cohort path:

```text
curated cohort manifest
  -> exact per-protein .LGS paths
  -> strict per-protein producer output load
  -> typed QtProtein + one-frame StaticPoseConformation
  -> Catalog/StaticCatalog + ResidentIndexes
  -> C++ source neighbourhoods/reductions/accumulators
  -> lean row-aligned per-(protein, atom) substrate under 15 GB
  -> Python fitter only
```

It consumes the producer output as given, keeps the extractor sacred, keeps the memory strategy open, and makes both forbidden Python paths structurally unavailable.
