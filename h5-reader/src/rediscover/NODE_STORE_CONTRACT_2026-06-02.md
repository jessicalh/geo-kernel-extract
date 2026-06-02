# Rediscover — Node-Store Contract + Pair-Index + Query Surface (the build target)

Status: **DESIGN, agreed 2026-06-02** (Jessica). This is the thing we build to.
It supersedes the loose `FIX_TO_VISION_2026-06-02.md` recovery sketch and the
"spine-reduces / Python-thins" framing as the *governing architecture*. No code moves
until this is the shape we're building. Discipline unchanged: model-is-the-spine,
no Python physics (except labeled fixtures), equivariant = real e3nn + frozen
`get_C()`, no factories/ABCs, **T2 is the thesis target and is preserved end-to-end.**

## The one idea

There is **one resident representation**: the typed model's **node store** (small) +
the spatial index. The **pairwise view is kept — as an index of pointers, not copies.**
A flat materialized expansion (the old ~68 GB) is a **transient tool generated on
demand for the one job that needs it, used, and dropped** — never resident, never the
default, never re-emitted as a matter of course. *Tool in the pocket, not a tyre around
the neck.* (Jessica, 2026-06-02.)

The 68 GB dump was wrong because it (a) **copied** every source's payload into every
(target, source) pair row, (b) pre-computed `disp/cosθ/dipolar/sum_*` that e3nn and the
C++ reducers compute themselves, and (c) was the *only* interface, forcing Python to
re-model. All three die here.

## 0. The functional API is the contract

The **functional API itself is the contract and the description**. Code-as-doc here means
the C++ verbs/relationships own association, reduction, selection, and the design-matrix
queries. They emit **named, described products**: a reducer answer, a node-store tensor,
a pair-index query, an e3nn-ready materialization, or a small exploration view.

A dump is allowed only as a **named-query output** — for example a compressed straight-into-e3nn
tensor packet, or a small view for exploration. It is never a raw blob that Python reads and
remodels into the protein again.

**THE LAW:** no maths model 2 in Python. C++ owns model semantics. Python owns the edges:
e3nn/PySR/ridge fitting, CV, plots, and the frozen `get_C()` basis map. The pybind11/live
chewer is **deferred and delivery-agnostic**; this API-as-contract decouples the model from
whether the next consumer is a binding, a compact file, or a temporary run-framework output.

## 1. The node store (the e3nn contract)

Per frame, a typed **heterogeneous point cloud**. Payloads live **once**, on the node.

| node | count (1P9J) | position | features (e3nn irreps) |
|---|---|---|---|
| **atom** | 846 | ✓ (`radius_graph`) | type ord → `0e` embed; charge q ×3 sources `0e`; E-field `1o`; EFG `2e` (MOPAC-Coulomb-derived); AIMNet2 embedding 256×`0e` (optional/separable) |
| **ring center** | 16 | ✓ | ring-type ord `0e`; **ring normal `1e`** (axial, sign-even); intensity `0e` |
| **bond midpoint** | 862 | ✓ | category ord `0e`; **bond axis as `2e`** (b̂⊗b̂ traceless — the C6 sign-safe even form, never `1o`); Δχ, order `0e` |

**Label** (per target atom): DFT shielding **`2e`** (5-dim, via the frozen `get_C()` map
from library `[xy,yz,zz,xz,xx−yy]`); `0e` T0 optional. T2 preserved, never collapsed.

**Edges:** built by the dataloader (`radius_graph` at the recorded cutoff); edge attr =
`spherical_harmonics(rel_pos)`; tensor products do the equivariant message passing +
aggregation. **Charge and field sit on atom nodes** (point charges/fields are at atom
centres) — not separate nodes.

**NOT in the node store** (e3nn / the C++ reducer compute these): `disp/r/cosθ/dipolar`,
`sum_dipolar*`/`field_z`/`field_mag`/`q_over_r3`, duplicated targets, bare-kernel scratch.

## 2. The pair index — the "deep pairwise" in the functional API

The pairwise view (Stage-1, de-circ, codex, the law-hunter all use it) **stays**, as a
**pointer adjacency**, promoted to a first-class queryable structure in the engine:

- A pair is `(target_idx, source_idx, source_kind)` — pointers into the node store. **No
  payload copied.** Dereference for the source's payload (single resident copy).
- Per-pair derived quantities (disp, cosθ, the literature kernel tensor) are **computed
  lazily on dereference** by the Layer-1 verbs — never stored.
- Self/bonded / near-field status is a pointer attribute (the producer filter rule),
  not a recomputed Python proxy (kills the `SELF_R=3.0` hack).
- This adjacency **is** e3nn's `edge_index` **and** the de-circ pairwise — one structure,
  all consumers. `RunTraversal` already generates these pairs; we expose them as a
  queryable `NeighborIndex` instead of consuming them to emit copies. (Concrete data
  structure, not a new abstraction.)

## 3. The query surface (delivery deferred, vectorized when live)

Analysis rides the functional API. The eventual live binding should be **vectorized** —
arrays per query, **never a Python call per pair** (or we trade disk bloat for call overhead
across ~30–60M pairs). The binding itself is deferred; the API contract below is not:

- `pairs(stratum | atom-set | frame-range) -> numpy arrays` (target id, source id/kind,
  lazily-computed disp/kernel) — the Stage-1 / de-circ pairwise dataframe, as a **small
  transient query result**, generated from the pointer index on demand.
- `reduced(relationship, stratum) -> aggregate kernel` — the C++ reducer answer (the
  thing Python was re-summing).
- `hunter(law) -> least-confounded cases` (#55) — the model answers the isolation /
  dominance / support query; the hunter does **not** scan raw rows.

This is the `project_reader_as_platform` chewer, but the chewer is not required to define
the contract. The Stage-1-style flat table doesn't disappear — it becomes a named, small
query result for the specific analysis.

## 4. Materialization is a tool in the pocket

If a job genuinely needs a **flat expanded table** (a full-pairwise global fit; a one-shot
e3nn tensor dump), it is **materialized on demand under the run framework**
(`tools/rediscover_run.py`, drop-old / auto-clean), **used, and dropped**. Even ~60 GB for
a few minutes is fine. What is forbidden is keeping it resident, reading it as the default
interface, or re-emitting it casually. Materialization is a verb, not a noun.

## 5. Sizes (real, 1P9J, 660 labeled frames)

- **Node store (resident):** AIMNet2 embedding 846×256×4B×660 ≈ **571 MB**; everything else
  (positions, atom/ring/bond features, T2 labels) ≈ **~80 MB**. Total **≈ 650 MB**, or
  **≈ 80 MB** if the embedding is stored separately/optional. (vs 68 GB.)
- **Pair index:** ~30–60M pointer pairs × ~9 B ≈ **0.3–0.5 GB** if held, or **~0** lazy
  (regenerated from the per-frame KD-tree).
- **Transient full flat expansion:** ~60 GB — **momentary, on demand, dropped.**

≈100× smaller resident, with the heavy expansion available as a tool when wanted.

## 6. Consumers

- **e3nn** — trains on the node store; dataloader builds edges via `radius_graph`. *This is
  the run Jessica is curious to see; it is NOT hostage to the full chewer — the node store
  can be materialized once for the first run.*
- **de-circ / law fits (ridge/PySR)** — query `pairs()` / `reduced()` for the stratum/cases
  they fit; small transient matrices.
- **law-hunter #55** — queries the pair index for least-confounded-in-motion cases
  (select on geometry + input modulation, never on the DFT target).
- **stats** — query the model; no flat-file reads, no Python re-modeling.

## 7. What changes (migration shape — gated, disk-aware; NOT a line plan yet)

1. Emit/expose the **node store** (one producer; positions + typed node features + T2
   labels), not the pair expansion.
2. Promote `RunTraversal`'s pairs to a queryable **pair index** (pointer adjacency); this
   also folds the remaining engine bypasses onto one traversal.
3. Define the **typed relationship index over all record shapes**: source-sum carriers,
   the richer all-atom raw-source carrier, and no-source per-target feature carriers.
4. **Vectorized query surface** over the API when delivery lands (pybind11/live chewer is
   deferred + agnostic; compact named materializations are allowed meanwhile).
5. Migrate the analysis scripts to **query** (delete the `np.add.at` reassembly, the
   re-derived strata, the `SELF_R` proxy, the scratch-column reads).
6. Run **e3nn** on the node store — the result we're after.

Gates throughout: oracle parity is an explicit command, not a bare `--case all` shortcut:
`python src/rediscover/analysis/oracle_parity.py --bin <build/linux-gcc/h5reader_extract> --run <explicit single .LGS path> --work <tmp> --case all --mc-cutoff 8.0`.
The run directory has multiple `.LGS` files, so pass the desired `.LGS` path directly.
`--case all` currently covers **ring + mc parity only**; it does **not** cover broad,
all-atom, efg, buckingham, or aimnet2. Follow with scoped
`QT_QPA_PLATFORM=offscreen ctest --test-dir build/linux-gcc -R h5reader_rediscover`.
Any materialization goes through the run framework (drop-old).

## 8. Open decisions (Jessica's calls)

- **Sequencing to the first e3nn run:** materialize the node store once for the e3nn run
  *before* the full chewer (fast path to the result), with the chewer following for the
  query consumers? Or chewer-first?
- **AIMNet2 embedding:** in the node store (≈571 MB) or a separable sidecar loaded only
  when used?
- **Frames:** 660 DFT-labeled for the e3nn contract (the dense 1501-frame MOPAC trajectory
  feeds the unlabeled/static features).

## Salvage from the day (unchanged)
Physics + fixes are KEEPERS (MOPAC leg, Welford `a703701`, provenance `cd8d710`, the T2
labels, the typed payloads). The throwaway is the **shape** (pair-expansion, scratch, the
half-built #29) — rewritten here as the node store + pair index. Not a day rollback.
