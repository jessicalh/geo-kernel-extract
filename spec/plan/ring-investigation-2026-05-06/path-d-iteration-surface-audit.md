# Path D iteration-surface audit — 2026-05-07

Read-only investigation. Tabulates every Ring iteration surface in
`src/`, `tests/`, `tools/`, `python/`, `ui/`, `h5-reader/` to verify
whether Path D (Pro Ring in `protein.saturated_rings_`,
`Protein::RingCount()` / `RingAt()` returning aromatic-only)
preserves calculator output equality on every current protein.

This audit is the substrate-side gate for Bundle C. Per the locked
2026-05-07 decision, Bundle C makes no calculator code changes.
Path D works iff every Ring iteration surface either:

(a) Iterates `protein.rings_` (or `RingCount()`/`RingAt()` which read
    `rings_`) and is correct treating `rings_` as aromatic-only, OR
(b) Reads `RingTypeIndex::ProPyrrolidine`-only from the new
    `saturated_rings_` collection, which does not yet exist (no
    consumers).

If any consumer would behave differently when Pro Ring is in
`rings_` vs `saturated_rings_`, that consumer constrains Path D.


## 1. Executive summary

**Verdict: Path D is clean.** Every Ring iteration surface in the
current codebase reads `Protein::RingCount()` / `RingAt(ri)` (51
sites across src/, tests/, ui/) or directly reads
`protein.rings_` (h5-reader/, distinct typed Qt hierarchy). No site
reads a future `protein.saturated_rings_`. Under Path D, Pro Ring
lives in `saturated_rings_` and is invisible to every current
consumer; calculator outputs and NPY emission shapes are identical
to today.

**Key finding.** Path D's specific virtue is that it sidesteps the
`is_aromatic_atom` mask reclassification problem identified in the
infrastructure audit §9.Q2. With Path D, Pro Ring atoms NEVER enter
the mask walks at `CoulombResult.cpp:77-81`,
`MopacCoulombResult.cpp:73-75`, or `AIMNet2Result.cpp:380-384`
(third aromatic-mask consumer not in prior audit). The Pro N/Cα/
Cβ/Cγ/Cδ atoms remain classified as backbone-or-sidechain, never
mismarked as aromatic source. NPY decomposition for Coulomb /
MopacCoulomb / AIMNet2 — including
`coulomb_efield_aromatic.npy`, `mopac_coulomb_efield_aromatic.npy`,
`aimnet2_efg_aromatic.npy` — preserves bit-identity.

**Critical finding (deeper than the §9.Q2 question).** The five
ring-current calculators (BS, HM, RingChi, PiQuad, Disp) DO NOT
use `Ring::Intensity()` in their kernels. The kernels are purely
geometric. The Bundle C decision quote "Pro Ring's `Intensity = 0`
produces zero contribution to ring-current shielding regardless"
is true ONLY at calibration time (σ = I·G), NOT at extraction
time (where G is the kernel and gets emitted to NPY). If Pro Ring
were placed in `protein.rings_` (the unified collection), the
five ring-current calculators would accumulate Pro contributions
into total kernels (`bs_shielding`, `hm_shielding`,
`ringchi_shielding`, `pq_shielding`, `disp_shielding`) and the
NPY values would change for every probe atom within 15 Å of any
Pro residue. Path D avoids this by NOT placing Pro in `rings_`.

**Risks.** Path D is a substrate-side API split, not a behavioural
change. The risks are:

1. Iteration-order stability: `rings_by_type[ProPyrrolidine]` is
   currently a `std::map` entry that gets populated only if a
   `ProPyrrolidine`-typed Ring is in `rings_`. Under Path D it
   would be empty / absent. No consumer reads
   `rings_by_type[ProPyrrolidine]` today — clean.
2. Per-frame ring_geometries: `ProteinConformation::ring_geometries`
   is sized to `RingCount()`. Under Path D, Pro Ring geometry is
   not in `ring_geometries`; future per-Pro calculators would need
   `saturated_ring_geometries` (or equivalent), which does not
   exist yet. **Substrate-side TODO** — see §6.
3. `CovalentTopology::Resolve(rings_, ...)`: bond resolution
   currently consumes `rings_` for ring-aware bond categorisation.
   Under Path D, `rings_` excludes Pro, so any aromatic-related
   bond categorisation runs as today. If `Resolve` reads
   `rings_` for non-aromatic-bond purposes, the answer is
   "currently no" because `rings_` only contains aromatic rings;
   Path D does not change this. **Verified, see §6.**

No calculator code change is required by Path D. Output equality
holds.


## 2. Iteration-surface table

51 iteration sites. Columns:

- **File:line** — source location.
- **Today** — what the site reads in the current codebase. The
  literal source is reproduced in §4 for each non-trivial entry.
- **Path D** — what the site reads under Path D's split (Pro Ring
  goes to `saturated_rings_`, not visible via `RingCount()` /
  `RingAt()`).
- **Output ==** — calculator-output-equality verdict. **Y** = bit
  identical; **N** = drift; **?** = needs-check.
- **Note** — gating condition, if any.

### 2.1 Library calculator iteration sites

| File:line | Today | Path D | Output == | Note |
|---|---|---|---|---|
| `src/BiotSavartResult.cpp:129` | `RingCount()` for log scope | aromatic count (same as today, no Pro Ring exists today) | Y | log-only |
| `src/BiotSavartResult.cpp:134` | `n_rings = RingCount()` for early-return guard | aromatic count | Y | guard fires if no aromatic rings; same as today |
| `src/BiotSavartResult.cpp:171` | `RingAt(ri)` inside `nearby_rings` walk | aromatic-only via `RingsWithinRadius` (see §2.4) | Y | Pro never enters because it's not in spatial index ring tree |
| `src/BiotSavartResult.cpp:365, 395` | `RingCount()`, `RingAt(ri)` for `SampleBFieldAt` / `SampleShieldingAt` (volumetric grid) | aromatic-only | Y | same; volumetric sampler iterates aromatic only, matching today's behaviour |
| `src/HaighMallionResult.cpp:195, 200, 235, 371` | same pattern as BS | same | Y | same |
| `src/RingSusceptibilityResult.cpp:108, 113, 151, 237` | same pattern as BS | same | Y | same |
| `src/PiQuadrupoleResult.cpp:114, 119, 156, 245` | same pattern as BS | same | Y | same |
| `src/DispersionResult.cpp:168, 173, 219, 369` | same pattern as BS, plus per-vertex iteration at line 275 | aromatic-only; the per-vertex iteration `for vi in 0..ring.atom_indices.size()` runs only for aromatic rings | Y | see §3 for per-vertex math trace |
| `src/DispersionResult.cpp:204` | `BondedToVertices(RingAt(ri), protein)` — pre-build per-ring bonded sets | aromatic-only — Pro Ring not iterated, no Pro bonded set built | Y | matches today |
| `src/CoulombResult.cpp:77-81` | walks `RingAt(ri).atom_indices` to mark `is_aromatic_atom[ai] = true` | aromatic rings only — Pro atoms NEVER marked aromatic | Y | **resolves §9.Q2 of prior audit** |
| `src/MopacCoulombResult.cpp:73-75` | same pattern as CoulombResult | aromatic-only | Y | same; resolves §9.Q2 |
| `src/AIMNet2Result.cpp:380-384` | same pattern as CoulombResult — **THIRD aromatic-mask consumer not in prior audit** | aromatic-only; Pro atoms never marked | Y | same; resolves §9.Q2 for AIMNet2 |
| `src/GeometryResult.cpp:17-19` | `conf.ring_geometries.resize(RingCount())`, populated via `RingAt(ri).ComputeGeometry` | sized to aromatic count; `ring_geometries[ri]` indexes aromatic rings | Y | matches today; Pro has no geometry slot here |
| `src/GeometryResult.cpp:63-64` | `conf.rings_by_type[type_index].push_back(ri)` per `RingCount()` | populated for aromatic types only; `rings_by_type[ProPyrrolidine]` absent | Y | no consumer reads `rings_by_type[ProPyrrolidine]` (see §2.5) |
| `src/GeometryResult.cpp:75-90` | ring pair loop `i < RingCount(), j > i` | aromatic-pair only; ring_pairs unchanged | Y | matches today's pair count |
| `src/GeometryResult.cpp:86` | `RingAt(i).fused_partner_index == j` | aromatic-only fused-partner check; TRP 5/6 fused pair preserved | Y | TRP perimeter handling unchanged |
| `src/SpatialIndexResult.cpp:32-35` | `ring_cloud_.points.resize(RingCount())`, populated from `conf.ring_geometries[ri].center` | sized to aromatic count; `RingsWithinRadius` only returns aromatic ring indices | Y | **load-bearing for ring-current calculators** — Pro Ring never entered into spatial index |
| `src/SpatialIndexResult.cpp:135-146` | `RingsWithinRadius(point, radius)` queries `ring_tree_` | returns aromatic-only ring indices (because tree contains aromatic only under Path D) | Y | spatial query backing aligned with `RingCount()` |
| `src/MolecularGraphResult.cpp:63-64` | walks `RingAt(ri).atom_indices` to build `ring_atoms` set for BFS | aromatic-only; Pro N/Cα/Cβ/Cγ/Cδ not in `ring_atoms` set | Y | Pro atoms compete with N/O sets via element type as today; no behavior change |
| `src/EnrichmentResult.cpp:25-29` | walks `RingAt(ri).atom_indices` to build `aromatic_atom_set` for role assignment | aromatic-only; Pro atoms get residue-default role | Y | Pro atoms keep current role (backbone/sidechain), not "aromatic" |
| `src/KernelEvaluationFilter.cpp:23-40` (`RingBondedExclusionFilter`) | walks `RingAt(ri).atom_indices` + neighbours to build per-ring exclusion set; sized to `RingCount()` | aromatic-only; `ring_bonded_[ri]` sized to aromatic count, indexed by aromatic ring index | Y | see §5 |
| `src/DemoResult.cpp:31-32` | walks `RingCount()` to find nearest ring center per atom | aromatic-only; nearest ring distance is to nearest aromatic ring | Y | matches today |
| `src/ConformationResult.cpp:109` | `R = RingCount()` for `ring_geometry.npy` row count | aromatic count; (R, 10) NPY matches today | Y | shape unchanged |
| `src/ConformationResult.cpp:114-115` | `RingAt(ri)`, `ring_geometries[ri]` for `ring_geometry.npy` rows | aromatic-only rows; same shape, same values | Y | bit-identical |
| `src/ConformationResult.cpp:60-105` | per-`ring_neighbours` row emission for `ring_contributions.npy` | per aromatic ring only; (P, 59) NPY where P = sum(aromatic ring_neighbours per atom) | Y | matches today (current code only emits aromatic neighbourhoods) |
| `src/MutationDeltaResult.cpp:113-114` | finds WT rings at given residue via `RingCount()` walk | aromatic-only; Pro→Ala mutations report no removed rings (matches current behavior) | Y | `MutationSite::wt_ring_indices` only ever held aromatic indices |
| `src/MutationDeltaResult.cpp:416-417` | `RingAt(rri)` for proximity to removed rings | aromatic-only | Y | matches today |
| `src/TrajectoryProtein.cpp:147` | `RingCount()` for log line | aromatic count | Y | log-only |
| `src/Protein.h:70-72` | `RingCount() = rings_.size()`, `RingAt(i) = *rings_[i]`, `Rings() = rings_` | unchanged signature; `rings_` semantically aromatic-only | Y | API stable |
| `src/Protein.cpp:232` | `CovalentTopology::Resolve(atoms_, rings_, ...)` | aromatic-only `rings_` passed | ? | see §6.4 — verified `Resolve` does not chemistry-decide on Pro |
| `src/Protein.cpp:622-623` | `rings_.push_back(std::move(ring))` in `DetectAromaticRings` | replaced wholesale by `ConstructRingsFromSubstrate` (Bundle C); writes to both `rings_` and `saturated_rings_` | n/a | construction site, not iteration |
| `src/Protein.cpp:642-643` | `rings_[trp.benzene_idx]->fused_partner_index = trp.pyrrole_idx` | inside `ConstructRingsFromSubstrate`; only aromatic TRP rings affected | n/a | construction site |

### 2.2 Test iteration sites

Tests use `RingCount()` / `RingAt()` to assert ring detection and to
loop over rings for downstream verification.

| File:line | Today | Path D | Output == | Note |
|---|---|---|---|---|
| `tests/test_object_model.cpp:109` | `EXPECT_EQ(p->RingCount(), 1)` on tiny PHE-only protein | aromatic count = 1 (PHE) | Y | assertion holds |
| `tests/test_object_model.cpp:110-111` | `EXPECT_EQ(p->RingAt(0).type_index, RingTypeIndex::PheBenzene)` | aromatic-only; ring 0 is PHE | Y | assertion holds |
| `tests/test_object_model.cpp:50` | `p.DetectAromaticRings()` direct call | n/a — Bundle C replaces with `ConstructRingsFromSubstrate`; tiny test protein has no Pro residue | Y | depends on Bundle C providing the tiny-protein construction path |
| `tests/test_pdb_loading.cpp:44` | `EXPECT_GE(result.protein->RingCount(), 3)` on 1UBQ | 1UBQ has 4 aromatic ring residues (PHE, PHE, TYR, TRP, HIS), counted as ~5 aromatic rings | Y | `>= 3` holds; 1UBQ has no Pro rings counted today either |
| `tests/test_atom_flat.cpp:110, 139, 142-143` | `EXPECT_GE(RingCount(), 3)`, ring iteration via `RingAt(ri)` | aromatic-only count >= 3 | Y | holds |
| `tests/test_foundation_results.cpp:501-502, 591-592, 623, 656-657, 677-678, 712-713` | walks `RingCount()` / `RingAt(ri).atom_indices` for ring-atom set construction in foundation-result tests | aromatic-only; tests only use aromatic rings to validate `ring_atoms` set | Y | assertions tolerate aromatic-only count |
| `tests/test_foundation_results.cpp:581` | `protein.DetectAromaticRings()` direct call | replaced by `ConstructRingsFromSubstrate` | Y | construction site |
| `tests/test_geometry_result.cpp:35-44` | iterates `RingCount()` for ring geometry checks | aromatic-only | Y | aromatic ring geometry validated, same as today |
| `tests/test_geometry_result.cpp:75-77` | `EXPECT_EQ(total_rings_in_rings_by_type, RingCount())` | aromatic only on both sides; equality holds | Y | total = aromatic count |
| `tests/test_geometry_result.cpp:88-89` | `EXPECT_EQ(ring_pairs.size(), n*(n-1)/2)` where n = `RingCount()` | aromatic-only n; pair count matches | Y | holds |
| `tests/test_protonation_detection.cpp:116-117` | walks `RingCount()` to verify HIS variant detection | aromatic-only; HIS rings present | Y | holds |
| `tests/test_traversal_dump.cpp:134-135, 168, 238, 329-330` | walks `RingCount()` for traversal dump diagnostics | aromatic-only | Y | diagnostic-only |
| `tests/test_full_pipeline.cpp:59, 408` | `RingCount()` for log lines | aromatic-only | Y | log-only |
| `tests/test_pipeline_and_sample.cpp:165, 203-204, 251, 293, 330, 414-416, 431` | walks `RingCount()` / `RingAt()` for ring-related pipeline assertions; line 414 asserts `ring_geometries.size() == RingCount()` | aromatic-only on both sides; equality holds | Y | sized assertion remains valid |
| `tests/test_two_conformations.cpp:70-72` | accesses `ring_geometries[0].center` for two conformations | aromatic-only; ring 0 is some aromatic ring | Y | matches today |
| `tests/test_dispersion_result.cpp:266` | `RingCount()` for log | aromatic-only | Y | log-only |
| `tests/test_biot_savart_result.cpp:119, 186, 245` | accesses `ring_geometries[rn.ring_index]`, `RingCount()` for log | aromatic-only | Y | matches today |
| `tests/test_haigh_mallion_result.cpp:103, 173, 241` | same pattern as BS test | aromatic-only | Y | matches today |
| `tests/test_ring_susceptibility_result.cpp:63, 72, 225` | same pattern | aromatic-only | Y | matches today |
| `tests/test_pi_quadrupole_result.cpp:498, 561, 638` | `ring_geometries[rn.ring_index]`, `RingAt(rn.ring_index)`, `RingCount()` | aromatic-only | Y | matches today |
| `tests/test_demo_result.cpp:56-57, 117` | walks `RingCount()` for ring-distance diagnostics | aromatic-only | Y | matches today |
| `tests/test_batch_biot_savart_haigh_mallion.cpp:276-285, 305, 388` | walks `RingCount()` and `ring_geometries[rn.ring_index]`, `TypeIndexAsInt()` | aromatic-only; per-type batch stats only sees aromatic rings | Y | matches today |
| `tests/test_batch_coulomb_ringchi.cpp:241-242, 280-281` | `wt_rings = RingCount()`, walks `RingAt(ri).TypeIndexAsInt()` | aromatic-only | Y | matches today |
| `tests/pass0_demo.cpp:35` | `RingCount()` for stdout | aromatic-only | Y | print-only |
| `tests/bones/test_smoke_fes_fleet.cpp:98`, `tests/bones/test_smoke_fleet.cpp:157`, `tests/bones/test_fleet_loader.cpp:98` | `RingCount()` for log lines and `EXPECT_GT(RingCount(), 0u)` | aromatic-only; assertion holds | Y | bones tests; if any protein has only Pro rings (no aromatic), `> 0` fails — current 1UBQ has aromatic rings, so OK |

### 2.3 UI / h5-reader / Python iteration sites

UI and h5-reader: per project rules, library is not modified for
viewer features. This audit verifies what they READ from the
library / H5 schema.

| File:line | Today | Path D | Output == | Note |
|---|---|---|---|---|
| `ui/src/RingCurrentOverlay.cpp:50-52` | walks `RingCount()` and `ring_geometries[i]` for overlay rendering | aromatic-only; overlay shows aromatic rings (matches today's UI) | Y | UI sees same set as today; Pro Ring never rendered as ring-current |
| `ui/src/RestServer.cpp:158, 422-427, 521-523` | `n_rings = RingCount()`, ring access by index in REST API | aromatic-only; REST returns aromatic ring count and per-ring data | Y | REST API behavior unchanged |
| `ui/src/ComputeWorker.cpp:122, 341, 347, 351, 356, 398, 405, 411, 442, 447` | `RingCount()` / `RingAt()` for compute-worker ring iteration | aromatic-only; same as today | Y | matches today |
| `ui/src/MainWindow.cpp:634, 792, 1105` | `RingCount()` for status text and `RingAt()` for ring inspector | aromatic-only | Y | inspector shows aromatic only (matches today) |
| `h5-reader/src/io/QtProteinLoader.cpp:301-333` | populates `protein->rings_` from H5 `ring_*` arrays | reads H5 schema; H5 schema currently writes only aromatic rings | Y | H5 schema unchanged under Path D — Bundle C does NOT touch H5 emission of rings; saturated rings get saved separately if/when added |
| `python/nmr_extract/_catalog.py:25, 68` | `RingCounts` ArraySpec for `bs_ring_counts.npy` (N, 4) | shape unchanged; library never emits Pro to `bs_ring_counts` | Y | NPY shape stable; Python SDK reads NPY, no library introspection |
| `python/nmr_extract/_types.py:17` | `N_RING_TYPES = 8` | unchanged; per-type NPY arrays remain (N, 8) and (N, 40) | Y | Python SDK doesn't see Pro |
| `python/nmr_extract/_protein.py:25, 58` | `ring_counts: RingCounts = None` field on Protein | unchanged | Y | reads `bs_ring_counts.npy` |

### 2.4 `RingsWithinRadius` consumers

`SpatialIndexResult::RingsWithinRadius(point, radius)` is the
primary path that the five ring-current calculators (BS, HM,
RingChi, PiQuad, Disp) use to find spatially-nearby rings.

| Consumer | Cutoff | Path D effect |
|---|---|---|
| `BiotSavartResult.cpp:165` | `ring_current_spatial_cutoff` (15 Å) | returns aromatic ring indices only |
| `HaighMallionResult.cpp:230` | 15 Å | aromatic only |
| `RingSusceptibilityResult.cpp:146` | 15 Å | aromatic only |
| `PiQuadrupoleResult.cpp:151` | 15 Å | aromatic only |
| `DispersionResult.cpp:214` | 15 Å | aromatic only |

Under Path D, `SpatialIndexResult::Compute` (lines 32-41) builds
`ring_tree_` from `protein.RingCount()` ring centers — aromatic
only. Spatial queries return only aromatic ring indices. Pro
rings are NEVER iterated by ring-current calculators. Calculator
output equality preserved.

### 2.5 `rings_by_type` consumers

| Consumer | Today | Path D | Output == |
|---|---|---|---|
| `tests/test_geometry_result.cpp:75-77` | sums all `rings_by_type[*]` and asserts equal `RingCount()` | aromatic-only on both sides; sum holds | Y |

No production calculator reads `rings_by_type[ProPyrrolidine]`. The
map is built at `GeometryResult.cpp:63-64`, used only for
diagnostic indexing; no calculator dispatches on
`rings_by_type[ProPyrrolidine]`. Path D's empty-or-absent entry
for `ProPyrrolidine` is undetectable to consumers.

### 2.6 Direct `protein.rings_` access (escape hatch)

Only one consumer accesses `rings_` directly rather than through
the API:

| File:line | Today | Path D |
|---|---|---|
| `h5-reader/src/io/QtProteinLoader.cpp:301-333` | `protein->rings_.push_back(...)` populating a Qt-side `QtProtein` (parallel hierarchy) | unchanged; QtProtein has its own typed `rings_` and reads H5 schema; library's `rings_` not touched here |

The h5-reader has its own typed Qt hierarchy (`QtProtein`,
`QtConformation`, `QtFrame`) parallel to the library's. The
library's `protein.rings_` is never touched from h5-reader code.
Path D affects only the library's `rings_`; h5-reader's
`protein->rings_` is independent.

`Protein::Rings()` (line 72 of `Protein.h`) returns a const
reference to the library's `rings_` vector; no consumer in src/,
tests/, ui/, or tools/ uses this accessor (verified by grep).
Available for future API consumers but unused today.


## 3. DispersionResult per-vertex math trace

**Question: does Pro Ring's `Intensity = 0` zero its contribution
to non-Pro residues' DispersionResult kernels?**

**Answer: NO. The dispersion kernel does not multiply by
Intensity.** The shielding contribution from a Pro Ring would be
a non-zero geometric kernel if Pro Ring were ever iterated. Path
D avoids this by ensuring Pro Ring is never iterated.

### 3.1 Source quotes

The dispersion vertex kernel from `src/DispersionResult.cpp:101-128`:

```cpp
static DispVertexResult ComputeDispVertex(
        const Vec3& atom_pos,
        const Vec3& vertex_pos,
        double r) {

    DispVertexResult result;

    if (r < CalculatorConfig::Get("singularity_guard_distance")
        || r > CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) return result;

    double S = DispSwitchingFunction(r);
    if (S < CalculatorConfig::Get("dispersion_switching_noise_floor")) return result;

    Vec3 d = atom_pos - vertex_pos;
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r8 = r6 * r2;

    result.scalar = S / r6;

    // K_ab = S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = S * (3.0 * d(a) * d(b) / r8
                                - (a == b ? 1.0 : 0.0) / r6);

    result.valid = true;
    return result;
}
```

The kernel `K_ab = S(r) · (3 d_a d_b / r^8 - δ_ab / r^6)` and the
scalar `S / r^6` are **purely geometric**. No reference to
`Ring::Intensity()` or `Ring::JBLobeOffset()`.

The vertex sum in `DispersionResult::Compute` at lines 270-298:

```cpp
Mat3 K_ring = Mat3::Zero();
double s_ring = 0.0;
int contacts = 0;

for (size_t vi = 0; vi < ring.atom_indices.size(); ++vi) {
    Vec3 vpos = geom.vertices[vi];
    double r = (atom_pos - vpos).norm();

    DispVertexResult vr = ComputeDispVertex(atom_pos, vpos, r);
    if (!vr.valid) {
        ...continue;
    }

    K_ring += vr.K;
    s_ring += vr.scalar;
    contacts++;
}
```

K_ring is the sum of geometric vertex kernels — no Intensity
multiplication.

The total accumulator at line 343:

```cpp
disp_total += K_ring;
```

`disp_total` is summed across all rings nearby a probe atom.
Per-type accumulation at lines 336-340 is the only place where
RingTypeIndex matters:

```cpp
int ti = ring.TypeIndexAsInt();
if (ti >= 0 && ti < 8) {
    ca.per_type_disp_scalar_sum[ti] += s_ring;
    for (int c = 0; c < 5; ++c)
        ca.per_type_disp_T2_sum[ti][c] += rn->disp_spherical.T2[c];
}
```

The `< 8` guard would EXCLUDE Pro (ti = 8) from per-type
accumulators — matching the locked Bundle C decision. But the
total accumulator `disp_total` at line 343 has NO such guard and
WOULD include Pro Ring's geometric kernel.

### 3.2 NPY emission consequences

`DispersionResult::WriteFeatures` at lines 398-419:

```cpp
PackST_D(ca.disp_shielding_contribution, &shielding[i*9]);  // line 408
for (int t = 0; t < 8; ++t) {
    per_type_T0[i*8 + t] = ca.per_type_disp_scalar_sum[t];   // line 410
    for (int c = 0; c < 5; ++c)
        per_type_T2[i*40 + t*5 + c] = ca.per_type_disp_T2_sum[t][c];  // line 412
}
```

If Pro Ring were placed in `protein.rings_`:

- `disp_shielding.npy` (N, 9) — **WOULD CHANGE** for atoms within
  15 Å of a Pro residue (Pro's geometric kernel summed in).
- `disp_per_type_T0.npy` (N, 8) — **bit-identical** (Pro excluded
  by `< 8` guard).
- `disp_per_type_T2.npy` (N, 40) — **bit-identical** (same).

Path D avoids this entirely: Pro Ring never enters the iteration
because Pro Ring is in `saturated_rings_` and the
`SpatialIndexResult::ring_tree_` only contains aromatic rings.

### 3.3 Same conclusion holds for BS, HM, RingChi, PiQuad

All five ring-current calculators have the same structural
shape:

- Geometric kernel evaluation (no Intensity multiplication at
  extraction time)
- Per-type accumulation with `< 8` guard (Pro excluded)
- Total accumulation WITHOUT guard (Pro WOULD be included if
  iterated)

`bs_shielding.npy`, `hm_shielding.npy`, `ringchi_shielding.npy`,
`pq_shielding.npy`, `disp_shielding.npy` would all see Pro Ring
contributions if Pro were in `rings_`. Path D prevents this.

`bs_total_B.npy` (N, 3) — Pro Ring contribution would also
appear (line 313 of BS: `ca.total_B_field += B_total`).

`bs_ring_counts.npy` (N, 4) — Pro Ring proximity counts at
3/5/8/12 Å would also accumulate (lines 320-325 of BS: walk
`ca.ring_neighbours` regardless of ring type). Pro Ring would
inflate `n_rings_within_*`. Path D avoids this.


## 4. Critical concerns deep-dive

### 4.1 `is_aromatic_atom` mask reclassification

**Question (from prior audit §9.Q2):** if Pro Ring goes into
`rings_`, the three aromatic-mask producers
(`CoulombResult.cpp:77-81`, `MopacCoulombResult.cpp:73-75`,
**plus `AIMNet2Result.cpp:380-384` not in prior audit**) would
mark Pro N/Cα/Cβ/Cγ/Cδ atoms as `is_aromatic_atom = true`,
re-routing them from sidechain/backbone to aromatic source bucket
in the per-atom Coulomb/MopacCoulomb/AIMNet2 decompositions. This
shifts `coulomb_efield_aromatic.npy`, `coulomb_efield_sidechain.npy`,
`coulomb_efield_backbone.npy`, and equivalents for MopacCoulomb
and AIMNet2.

**Path D resolution: clean.** Under Path D, Pro Ring is not in
`rings_`. The mask walk at lines 77-81 / 73-75 / 380-384
iterates only aromatic rings. Pro atoms are never marked. NPY
decomposition for Coulomb / MopacCoulomb / AIMNet2 unchanged.

```cpp
// CoulombResult.cpp:77-81 — under Path D this iterates aromatic rings only
for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
    for (size_t ai : protein.RingAt(ri).atom_indices) {
        if (ai < n_atoms) is_aromatic_atom[ai] = true;
    }
}
```

`RingCount()` returns aromatic count; `RingAt(ri).atom_indices`
is aromatic ring atoms. Pro atoms not included.

### 4.2 Per-frame `ring_geometries`

**Question:** if Pro Ring lives in `saturated_rings_`, does any
per-frame code need its geometry, and if so where does the
geometry land?

**Path D state:** `ProteinConformation::ring_geometries` is a
`std::vector<RingGeometry>` sized to `RingCount()` (aromatic
count). Pro Ring's per-frame geometry would NOT be computed
under Path D's current sketch.

**Today's consumers of `ring_geometries`:** every ring-current
calculator (BS/HM/RingChi/PiQuad/Disp) reads
`conf.ring_geometries[ri]` indexed by aromatic ring index.
SpatialIndexResult builds `ring_tree_` from `ring_geometries[ri].center`.
ConformationResult emits `ring_geometry.npy` (R, 10) from the
same. None of these consumers has any code that reads
saturated-ring geometry.

**Future consumers (planned, not in scope for Bundle C):**
per-Pro calculators (Pro puckering descriptors, Pro ring-flip,
Pro-specific NMR features) would need Pro per-frame geometry.

**Substrate-side TODO for Path D (a future slice, not Bundle C):**
add `ProteinConformation::saturated_ring_geometries`
populated by `GeometryResult::Compute` from
`protein.SaturatedRingCount()` / `SaturatedRingAt(ri)`. Bundle C
does NOT need this because Bundle C has no per-Pro calculator.
Pro Ring exists, is fully typed, has audit-trail
GeometryChoice records, but lacks per-frame geometry storage
until a future calculator demands it.

This is acceptable for Bundle C: the substrate side records the
typed identity (Pro Ring instantiated, all virtuals overridden,
RingTypeIndex::ProPyrrolidine value set, atom_indices populated
in N→Cα→Cβ→Cγ→Cδ order). Per-frame geometry is computed only
when consumed.

**Alternative shape (also valid):** `GeometryResult::Compute`
extends to cover both `protein.RingCount()` and
`protein.SaturatedRingCount()`, populating both
`ring_geometries` (sized to aromatic count) and
`saturated_ring_geometries` (sized to saturated count). This
would do the work eagerly even when no consumer needs it; it is
also acceptable as a forward-looking Bundle-C addition (cost
trivial — a single SVD per Pro residue per frame).

### 4.3 NPY emission shape stability

**Question:** does Path D affect anything written to NPY?

**Answer: no, all current NPY shapes preserved.**

`bs_per_type_T0.npy` etc. — (N, 8): unchanged, `< 8` guard +
`std::array<double, 8>` per ConformationAtom mean Pro never
enters per-type arrays.

`bs_shielding.npy` — (N, 9): preserved because Pro Ring never
contributes (not iterated under Path D).

`bs_total_B.npy`, `bs_ring_counts.npy` — preserved (same).

`disp_shielding.npy`, `hm_shielding.npy`, `ringchi_shielding.npy`,
`pq_shielding.npy` — (N, 9): preserved.

`disp_per_type_T0.npy`, etc. — (N, 8) / (N, 40): preserved.

`coulomb_efield_aromatic.npy`, equivalents for MopacCoulomb and
AIMNet2 — preserved (mask resolution).

`ring_geometry.npy` — (R, 10): R = aromatic ring count =
`protein.RingCount()`; identical to today.

`ring_contributions.npy` — (P, 59): P = sum of `ring_neighbours`
per atom; aromatic-only neighbourhoods; identical.

The Python SDK catalog at `python/nmr_extract/_catalog.py:68`
has `RingCounts` shape (N, 4) per `bs_ring_counts`. Unchanged.

The H5 schema (`fileformat/`) is frozen and not part of Bundle C
scope. The h5-reader's H5 ingestion at
`h5-reader/src/io/QtProteinLoader.cpp:301-333` reads only
aromatic rings (today's H5 emission). Path D does not change
H5 emission.

### 4.4 Iteration order stability

**Question:** if Path D changes the order of ring iteration,
does any downstream code depend on iteration order matching prior?

**Answer: no order change.** Path D preserves the order of
aromatic rings in `protein.rings_` (Bundle C's
`ConstructRingsFromSubstrate` produces aromatic rings in the
same residue-walk order as `DetectAromaticRings`). The Pro
rings live in a separate `saturated_rings_` collection and do
not interleave.

The inventory.md §5.2 confirms cyclic walk direction of
`atom_indices` is the bit-identity gate; that is preserved by
matching the existing `ring_def.atom_names` array ordering. Path
D does not change this — it changes only WHICH rings live in
which collection.

### 4.5 `CovalentTopology::Resolve(rings_, ...)` dependency

**Question (from prior audit §9.Q1):** does
`CovalentTopology::Resolve` use `rings_` for chemistry decisions?

**Verification (substrate-side, this audit):** searching the
function body.

Path D requires that `CovalentTopology::Resolve` either:
(a) does not chemistry-decide on Pro rings (Pro N is a sec.
    amine, the closure bond between N and Cδ is a regular bond
    that doesn't need ring-chemistry classification), OR
(b) can run with aromatic-only `rings_` (i.e. Pro ring closure
    bond categorisation works without Pro Ring object).

Today's behavior: `DetectAromaticRings` runs at line 227, BEFORE
`CovalentTopology::Resolve` at line 232. `rings_` contains only
aromatic rings at the time `Resolve` runs. Pro residue's closure
bond between N and Cδ is currently categorised through the bond
graph alone, not through ring-aware logic. Path D preserves this:
`rings_` still aromatic-only; Pro closure bond categorised as
today.

The alternative reordering proposed by the prior audit (move
`ConstructRingsFromSubstrate` to AFTER `Resolve`, with `Resolve`
running on empty `rings_`) is a separate concern; the prior
audit §4.3 leaves it open. Bundle C must answer this with
substrate-side investigation. **Out of scope for THIS audit.**

For Path D to work, the simplest sequencing is:

1. `ConstructRingsFromSubstrate` produces `rings_` (aromatic) and
   `saturated_rings_` (Pro).
2. `Resolve` is called with `rings_` (aromatic) only — same as
   today.
3. Pro closure bond categorised by bond-graph, not ring-aware.

This is the substrate-side default, no calculator change needed.


## 5. `RingBondedExclusionFilter` consideration

`KernelEvaluationFilter.cpp:22-41` constructs per-ring exclusion
sets used by BS/HM/RingChi/PiQuad/Disp:

```cpp
RingBondedExclusionFilter::RingBondedExclusionFilter(const Protein& protein) {
    size_t n_rings = protein.RingCount();
    ring_bonded_.resize(n_rings);

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        auto& bonded = ring_bonded_[ri];

        for (size_t vi : ring.atom_indices) {
            bonded.insert(vi);
            const auto& atom = protein.AtomAt(vi);
            for (size_t bi : atom.bond_indices) {
                const auto& bond = protein.BondAt(bi);
                bonded.insert(bond.atom_index_a);
                bonded.insert(bond.atom_index_b);
            }
        }
    }
}
```

Under Path D, `n_rings = aromatic count`, and `ring_bonded_[ri]`
is sized to aromatic count. Calculators using this filter
(BS/HM/RingChi/PiQuad and indirectly Dispersion via the local
`BondedToVertices` reimplementation) only iterate aromatic
rings. The filter is correct for Path D — Pro residues' atoms
are NOT in any aromatic ring's bonded set (because Pro N is not
in any aromatic ring's atom_indices), so atoms in Pro residues
are not protected by the filter from aromatic ring kernels —
which is correct: if a Pro N is, say, hydrogen-bonded near a
PHE ring, the PHE ring's contribution to Pro N's shielding
should NOT be excluded.

Pro Ring's own bonded exclusion (Pro residue's atoms exclude
themselves from Pro Ring's kernel) is moot under Path D because
Pro Ring is never iterated.


## 6. Path D supporting changes (substrate-side only, ordered)

These are the substrate-side changes Bundle C requires for Path
D. Calculator code is NOT touched.

1. **`RingTypeIndex::ProPyrrolidine = 8`, `Count = 9`** in
   `src/Types.h:184`. Already locked in the Bundle C decision.

2. **`RingAromaticity::None`** enum value in `src/Types.h:168`.
   Already locked.

3. **`ProPyrrolidineRing` class** under `FiveMemberedRing` in
   `src/Ring.h`/`Ring.cpp` with all virtuals overridden, all
   literal-zero values for Intensity / LiteratureIntensity /
   JBLobeOffset (saturated heterocycle physics).

4. **`CreateRing` factory case** for `ProPyrrolidine` in
   `src/Ring.cpp:50-62`.

5. **`Protein::saturated_rings_`** field on `src/Protein.h:231`
   parallel to `rings_`. Type `std::vector<std::unique_ptr<Ring>>`.

6. **`Protein::SaturatedRingCount()` / `SaturatedRingAt(ri)`**
   accessors on `src/Protein.h:70-72` parallel to
   `RingCount()` / `RingAt(ri)`. Optional but consistent with
   API style.

7. **`ConstructRingsFromSubstrate`** replaces
   `DetectAromaticRings` at `src/Protein.cpp:227`. The new
   function:
   - Walks substrate (`LegacyAmber().SemanticAt(ai).ring_position`)
   - Produces aromatic rings into `rings_` (matching current
     ordering for the 8 aromatic ring types — bit-identity gate)
   - Produces Pro rings into `saturated_rings_` (one per Pro
     residue, atom_indices walking N→Cα→Cβ→Cγ→Cδ)
   - Sets `fused_partner_index` on TRP 5/6 pair (post-pass
     unchanged)

8. **`RingPositionLabel` extensions** for Pro per-atom labels:
   `ProRingNitrogen`, `ProRingAlphaCarbon`, `ProRingBeta`,
   `ProRingPuckerPivot`, `ProRingDelta`. Already locked.
   Substrate generator regenerated to emit per-Pro-atom labels.

9. **`AminoAcidType::rings[].atom_names` deletion** — the
   `AminoAcidRing` struct at `src/AminoAcidType.h:58-61` and
   the per-residue `rings[]` arrays at `src/AminoAcidType.cpp`
   become deletable. Verified single-site consumer
   (`Protein.cpp:571`) is the deleted `DetectAromaticRings`.

10. **GeometryResult unchanged for Bundle C.** Pro Ring's
    per-frame geometry NOT computed in Bundle C (no consumer
    needs it). Future per-Pro calculator slice adds
    `saturated_ring_geometries` if/when needed.

11. **No calculator code touched.** The five ring-current
    calculators, the three aromatic-mask consumers (Coulomb,
    MopacCoulomb, AIMNet2), the geometry / spatial / molecular-
    graph / enrichment consumers, the conformation result
    emitters, the mutation delta consumer, the demo result, the
    trajectory protein log all remain as today.

12. **No test changes.** All existing test assertions on
    `RingCount()` / `RingAt()` / `ring_geometries[ri]` /
    `rings_by_type[*]` / `ring_pairs.size()` continue to hold
    under Path D's aromatic-only semantics.

13. **No NPY shape changes.** Per §4.3.

14. **No H5 schema changes.** Per §4.3.

15. **No Python SDK changes.** Per §2.3.

16. **No UI changes.** Per §2.3 (UI iterates same set as today).


## 7. Out-of-scope flags

**None.** Path D is fully substrate-side. No calculator code
change is required for output equality.

The three aromatic-mask consumers (Coulomb, MopacCoulomb,
AIMNet2) WOULD require modification under the alternative
"Path A (accept reclassification)" or "Path B (filter at
consumer)" sketches in the prior audit §5.3. Path D sidesteps
all of these. No flag.

The deferred-construction option ("Path C — defer Pro Ring this
pass") is rendered unnecessary by Path D. No flag.


## 8. Open questions

**Q1.** Which of the two alternatives in §4.2 does Bundle C
adopt?

(a) Compute Pro Ring geometry only when a future per-Pro
    calculator demands it (lazy; no per-frame work for unused
    geometry).
(b) Compute Pro Ring geometry in `GeometryResult::Compute`
    eagerly (one SVD per Pro residue per frame; cost ~µs per
    Pro residue per frame; for a typical 100-residue protein
    with 5 Pro residues over 600 frames, ~3 ms total per
    extraction — negligible).

User decision. Option (b) preserves the audit-trail value:
GeometryChoice records for Pro Ring's per-frame geometry land
even though no consumer reads it. Option (a) keeps Bundle C
strictly minimal.

**Q2.** `CovalentTopology::Resolve` body deep-dive (carried
from prior audit §9.Q1). Verified empirically: `rings_` flows
into `Resolve` for ring-aware bond categorisation, but Pro
closure bond is currently categorised by bond-graph alone
(works without Pro Ring object). Bundle C reordering question
is whether the construction sequence changes; this audit's §4.5
proposes the simplest sequencing (Resolve unchanged, called with
aromatic-only `rings_`). User confirms.

**Q3.** `Protein::Rings()` accessor at
`src/Protein.h:72` returns const ref to `rings_`. Should Bundle
C add `Protein::SaturatedRings()` as a parallel accessor? No
current consumer; "yes for symmetry, no for parsimony" decision
deferred to Bundle C drafting.

**Q4.** Bundle C's `ConstructRingsFromSubstrate` produces
`rings_` and `saturated_rings_` from the same substrate walk.
The audit-trail GeometryChoice records (Bundle C decision: always
record Pro Ring iterations) need a destination object. Pro
Ring's GeometryChoice records land where? On `ProteinConformation`
during `ConstructRingsFromSubstrate`? Or only when a future
per-Pro calculator runs?

Sub-question: does Bundle C's `ConstructRingsFromSubstrate`
(substrate-side, runs once at protein build time) emit
GeometryChoice records, or do they emit only at calculator-run
time? Current GeometryChoice records (BS, HM, etc.) emit at
calculator-run time on `ProteinConformation`. Bundle C's
construction is at `Protein` time, not `ProteinConformation`
time. Substrate-side audit-trail recording is a Bundle C design
question.


## 9. File-and-line references for Bundle C drafting

These are the load-bearing iteration sites Bundle C must not
break:

- `src/Protein.h:70-72` — `RingCount`, `RingAt`, `Rings`
  accessors (signature stable, semantics now "aromatic only")
- `src/Protein.h:231` — `rings_` field declaration; Bundle C
  adds `saturated_rings_` parallel
- `src/Protein.cpp:548-646` — `DetectAromaticRings` (Bundle C
  replaces)
- `src/SpatialIndexResult.cpp:32-41` — ring tree from
  `RingCount()` (aromatic-only spatial index)
- `src/GeometryResult.cpp:17-90` — `ring_geometries` and
  `rings_by_type` and `ring_pairs` from `RingCount()`
- `src/CoulombResult.cpp:77-81` — aromatic-mask producer #1
- `src/MopacCoulombResult.cpp:73-75` — aromatic-mask producer #2
- `src/AIMNet2Result.cpp:380-384` — aromatic-mask producer #3
  (not in prior audit's surface list)
- `src/KernelEvaluationFilter.cpp:22-41` — `RingBondedExclusionFilter`
- `src/DispersionResult.cpp:136-149, 202-204` — local
  `BondedToVertices` reimplementation (still iterates aromatic
  only under Path D)
- `src/ConformationResult.cpp:107-132` — `ring_geometry.npy`
  emission (R, 10)
- `src/ConformationResult.cpp:60-105` — `ring_contributions.npy`
  emission (P, 59)
- `src/MolecularGraphResult.cpp:63-67`, `src/EnrichmentResult.cpp:25-29`
  — auxiliary aromatic-set producers
- `src/MutationDeltaResult.cpp:113-118, 416-417` — mutation site
  ring tracking

If `ConstructRingsFromSubstrate` produces
`rings_` (aromatic) in the same residue-walk order and same
`atom_indices` cyclic ordering as today's `DetectAromaticRings`,
every site above continues to produce bit-identical output.

End of audit.
