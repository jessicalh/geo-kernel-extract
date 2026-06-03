# Rediscover emit-surface audit — A/B/C/D buckets (2026-06-02)

Read-only audit on branch `h5-reader-pysr-spike`. **Nothing changed.** All
decisions reserved for the lead. This is the concrete form of the
"#29 model-placement / C++↔Python boundary" question: which emitted columns are
model-blessed physics predictions (the `jb_T*` ideal), which are legitimate raw
inputs, which are geometric-intermediate reassembly-scratch (the `sum_dipolar_*`
smell), and which are redundant.

## The lens, restated

The discipline (`feedback_model_is_spine`, `feedback_no_python_physics_except_
labeled_integrity_test`, the per-cell workflow in `REDISCOVERY_MAP.md` step 2):
the typed C++ model holds topology, typed rings/bonds, geometry, charges; the
spine asks IT the physics question and emits the literature-scaled PREDICTION;
Python only correlates. RING is the exemplar done right — per-ring-type
`LiteratureIntensity()` on the typed `QtRing` → emitted `jb_T*` (literature-scaled
ppm) → `ring_literature_decirc.py` only correlates it vs DFT. The counter-example
is the per-category geometric scalar sum that a Python script weights into physics.

**Headline finding.** The smell is real and it is *structural*: every
`SourceSum` relationship emits BOTH a model-blessed kernel (good) AND a stack of
per-category / per-type geometric scalar sums (`sum_dipolar_*`, `ring_sum_dipolar`,
`bond_sum_dipolar`, `field_z`, `field_mag`) that exist only so a Python OLS can
weight raw geometry into a coefficient. The codebase already KNOWS this — the
consumers and `REDISCOVERY_MAP.md` step 2 explicitly say the scalar sums are
"bare," "untyped and unweighted," "not directly comparable," and that the
intensity-weighted aggregate "is a C++ reducer to add, not a Python re-sum."
The C-bucket below is that backlog made into a column-by-column list.

A second, orthogonal finding: there is real **emit duplication** (bucket D) — the
ring path emits its T2 target three times (CSV columns + a sources NPY + an
aggregated NPY, ×{lab, local}); the broad path and the per-atom-feature paths each
re-derive their own NPY writer and re-emit the DFT target on every sink.

---

## Scope: what actually emits (2026-06-02)

Eight cases reach a sink. Three relationships are fail-loud stubs (no emit).

| case | run fn | sink | file |
|---|---|---|---|
| `ring_current` | `RingCurrentNeighborhood::extract` | `RecordSink` | RingCurrentNeighborhood.cpp |
| `mcconnell` | `McConnellNeighborhood::extract` | `RecordSink` | McConnellNeighborhood.cpp |
| `charge_dipole` | `ChargeDipoleNeighborhood::extract` | `RecordSink` (+`WriteChargeDipoleAggregatedRow`) | ChargeDipoleNeighborhood.cpp |
| `broad_backbone` | `RunBroadBackbone` | `BroadBackboneSink` | BroadBackbone.cpp / BroadBackboneSink.cpp |
| `efg` | `RunEfgPerAtomFeature` | `EfgFeatureSink` | EfgFeature.cpp / EfgFeatureSink.cpp |
| `buckingham_efield` | `RunBuckinghamEfieldPerAtomFeature` | `BuckinghamEfieldSink` | BuckinghamEfield.cpp / BuckinghamEfieldSink.cpp |
| `aimnet2_features` | `RunAimnet2PerAtomFeature` | `Aimnet2FeatureSink` | Aimnet2Feature.cpp / Aimnet2FeatureSink.cpp |
| (stubs) `charge_quadrupole`, `larsen_hbond`, `aimnet2_embedding` | — | fail-loud, no emit | main_extract.cpp:65-95 |

There are **two emit families** (a known coupling debt, `BroadBackbone.h:13-21`,
`#29`): the `RecordSink` family (ring/mc/charge_dipole, shared schema via
`FeatureSchema`) and four bespoke per-case sinks (Broad/Efg/Buckingham/Aimnet2),
each with its own hand-rolled `writeNpyF64` and its own CSV header string. The
duplication in bucket D below is largely a symptom of that split.

The Python SDK (`python/nmr_extract/_catalog.py:196+`) registers **only the NPY
sidecars** (`rediscover_*_target_T2`, `_target_local_T2`, `_bare_kernel_T2`); the
CSV columns are read directly by the analysis scripts, not through the SDK. So
"cut a CSV column" is a producer+consumer change, not an SDK-contract change.

---

## (A) MODEL-BLESSED PHYSICS PREDICTION — KEEP

A literature-scaled prediction a typed object computes (the `jb_T*` ideal). These
are the destination the whole surface should converge on.

| quantity | what computes it | file:line | why it's blessed |
|---|---|---|---|
| `jb_T0`, `jb_T2_local_0..4` (ring source) | `ScaleSphericalTensor(jb_unit, QtRing::LiteratureIntensity())` | RingCurrentNeighborhood.cpp:271-272; columns RingCurrentNeighborhood.cpp:87-92 | per-ring-type literature intensity from the typed `QtRing` virtual, applied C++-side; ppm. THE exemplar. |
| `jb_unit_T0`, `jb_unit_T2_local_0..4` (ring source) | `JohnsonBoveySourceUnitKernelLocal(...)` | RingCurrentNeighborhood.cpp:268-270; columns RingCurrentNeighborhood.cpp:81-86 | the unit-current (per-nA/T) source kernel in the target local frame — the bare form the literature constant scales. Diagnostic half of the blessed pair. |
| `literature_kernel_T2_0..4` (broad, summed) | `sumKernelT2(ring,bond,charge)` | BroadBackbone.cpp:475-477; BroadBackboneSink.cpp:219 | total fixed-coefficient T2 prediction over all mechanisms, the de-circularising target. |
| `ring_literature_kernel_T2_0..4` (broad) | `h5KernelT2Local(KernelBs)` | BroadBackbone.cpp:468-469; sink BroadBackboneSink.cpp:220 | producer BS ring-current kernel rotated into the typed backbone frame, C++-side. |
| `bond_literature_kernel_T2_0..4` (broad) | `h5KernelT2Local(KernelMc)` w/ `bondKernelT2FromSources` fallback | BroadBackbone.cpp:470-473; sink BroadBackboneSink.cpp:221 | McConnell kernel in the local frame; fallback rebuilds the dipolar tensor from the typed bond sources. |
| `charge_literature_kernel_T2_0..4` (broad) | `chargeKernelT2FromSources` (FF14SB Coulomb EFG, traceless) | BroadBackbone.cpp:474; def BroadBackbone.cpp:107-127; sink BroadBackboneSink.cpp:222 | fixed-form EFG-like T2 built C++-side from the attached charge sources. |

Consumers that treat these as the answer (correlate-only): `ring_literature_decirc.py`
(`jb_T*` / `jb_unit_T*`, with a self-audit that `jb_unit_T0 * ring_intensity == jb_T0`,
ring_literature_decirc.py:644-655), `backbone_literature_kernel_t2_corr.py` (all four
broad `*_literature_kernel_T2`), `johnson_bovey_region_recovery.py`. This is the
right shape — Python groups by `(atom,frame[,ring])` and correlates, no physics.

**Caveat (note, not a decision): the bare-kernel columns are an over-broad A.**
`bare_T0`, `bare_T1_0..2`, `bare_T2_0..4` (`BareKernelColumns()`,
ExtractionSupport.cpp:91-99) are the producer's per-atom BS/MC kernel read straight
from H5 — a legitimate cross-check (closer to bucket B "raw input the model
exposes"). They are emitted on BOTH the sources CSV and the aggregated CSV for
ring + mc (RecordSink.cpp:131-135, 224, 248) = 9 columns × 2 row kinds × 2 files,
plus the `_bare_kernel_T2` NPY ×2. See bucket D for the duplication.

---

## (B) LEGITIMATE RAW INPUT the model exposes — KEEP (with justification)

Positions, geometry (r, cosθ), the source's orientation vector, identity, and the
DFT target. The fitter and the de-circularising step genuinely need these; they are
not physics Python reassembles.

**Identity / provenance (every row, both families).** `atom_index`, `residue_index`,
`residue_number`, `amino_acid_ord`, `element_ord`, `atom_name`, `stratum`/
`backbone_frame_class`, `h5_row`, `original_index`, `time_ps`
(`IdentityColumns()`, ExtractionSupport.cpp:70-89; `FillIdentity`,
ExtractionSupport.cpp:15-43). KEEP — join keys + per-stratum grouping; `atom_name`
is display-only per `feedback_ui_data`.

**Local frame (every row).** `frame_z/x/y_{x,y,z}` (9), `frame_variant`,
`frame_valid`, `frame_anchor_atom_index` (ExtractionSupport.cpp:82-87). KEEP — the
basis the source vectors and `total_local` live in; an equivariant fit needs the
frame recorded (`SURFACE_DESIGN.md` #7). The anchor index is the IUPAC-trap-safe
typed provenance (`selectUnique`, RingCurrentNeighborhood.cpp:177-195).

**Per-source geometry & orientation.** Ring: `disp_local_{x,y,z}`, `r`, `cos_theta`,
`ring_z`, `ring_rho`, `ring_in_plane_angle`, `source_normal_local_{x,y,z}`
(RecordSink.cpp:137-152). Bond: `disp_local_*`, `r`, `cos_theta_bond_axis`,
`bond_axis_local_{x,y,z}` (RecordSink.cpp:154-163). Charge: `disp_local_*`, `r`,
`source_q_e` (RecordSink.cpp:165-173). KEEP — these are exactly the per-source
inputs the equivariant per-source-type radial fitter consumes
(`equiv_t2_backbone_e3nn.py`); the orientation vectors are a deliberate emit
(`REDISCOVERY_MAP.md` step 2, "#32"). NOTE: `dipolar_3cos2m1_over_r3` is the
borderline case — see (C); it is a *pre-combined* geometric quantity, not a raw input.

**Source identity (typed virtuals, raw).** Ring: `ring_type_index`, `ring_intensity`
(=`LiteratureIntensity`), `ring_nitrogen_count`, `ring_jb_offset`,
`ring_aromaticity_ord`, `ring_size`, `ring_fused`, `is_self_or_bonded`, `ring_index`
(RecordSink.cpp:141-145). Bond: `bond_category/order/elem_a/elem_b/index/atom_a/atom_b`
(RecordSink.cpp:155-159). Charge: `source_atom_index/residue_*/amino_acid/element/
atom_name` (RecordSink.cpp:166-170). KEEP — typed-object answers (the `QtRing`/`QtBond`
virtuals), used as stratifiers/keys, not reassembled. `ring_intensity` is the literature
constant ITSELF, exposed raw so the fitter can recover it un-circularly
(`ring_literature_decirc.py` recovers γ_bare vs the Pople value).

**DFT target (every row, the shared target).** `dft_present`, `dft_total_raw_*` (9),
`dft_dia_raw_*` (9), `dft_para_raw_*` (9), `dft_sigma_iso`, `dft_total_T1_0..2`,
`dft_total_T2_0..4`, `dft_dia_iso`, `dft_para_iso`, `dft_total_local_*` (9),
`dft_local_frame_valid` (`TargetColumns()`, ExtractionSupport.cpp:164-185;
`BuildTarget` ExtractionSupport.cpp:45-68). KEEP — the quantum reference; the
decomposition `DecomposeLibrary` lands DFT-T2 in the same basis as the kernel T2
(not a Python reprojection). NOTE: `dia_raw`/`para_raw`/`dia_iso`/`para_iso` are
labelled "diagnostics" (RediscoverTypes.h:27); none of the read scripts consume them
(grep: no hits) — borderline D (carried-but-unread), but cheap and physically
honest to keep. T1 is flagged `unverified` (`SURFACE_DESIGN.md` #4) — keep + flag.

**Field/feature raw inputs (per-atom-feature cases).** EFG: `efg_feature_lab_T2`,
`efg_feature_T2` (APBS EFG read from H5, rotated to local; EfgFeature.cpp:212-217),
`efg_T2_magnitude`, `efg_units`. Buckingham: `efield_local_{x,y,z}`, `E_proj`,
`E_mag`, `efield_units` (BuckinghamEfield.cpp:191-198). AIMNet2: `aimnet2_charge`,
`crg_scalar`, `crg_vector_local/lab`, `embedding` (Aimnet2Feature.cpp:227-247).
KEEP — these ARE the source data for those mechanisms (APBS field is the source for
#4, not a cross-check — `SURFACE_DESIGN.md` carrier note); the only transform is a
rigid frame rotation, done C++-side in the typed frame. AIMNet2 CRG/embedding are
learned reps carried as-is (correctly labelled "not_polarizability",
"learnable_ceiling_feature_not_physical_law", Aimnet2FeatureSink.cpp:180-183).

---

## (C) GEOMETRIC INTERMEDIATE FOR PYTHON REASSEMBLY — CANDIDATE TO CUT/CONSOLIDATE

Scalar sums / per-category pieces that exist only so a Python script can weight or
combine them into physics — the `sum_dipolar_*` smell. For each: the reassembling
consumer + the model-blessed emit (in the McConnell-Δχ → `jb_T*` image) that should
replace it.

### C1 — ring per-ring-type dipolar sums `sum_dipolar_ringtype_0..7`
- **Emit:** RingCurrentNeighborhood.cpp:107-109 (8 columns, aggregated CSV);
  computed RingCurrentNeighborhood.cpp:280-283 (unweighted `Σ(3cos²θ−1)/r³` per type).
- **Reassembled by:** `look03_coefficient.py` (per-type geometry split,
  look03_coefficient.py:5-7), `diag_differencing.py`.
- **Why it's a smell:** the per-type split exists so Python can apply a per-type
  weight (the ring intensity) — i.e. reassemble the intensity-weighted ring kernel
  outside the model. `look03_coefficient.py:86-87` says so explicitly: *"intensity-
  weighted aggregate intentionally omitted (it is a C++ reducer to add, not a Python
  re-sum of `ring_intensity*dipolar`)."*
- **Model-blessed replacement:** the per-ring-type literature-scaled `jb_T*` ALREADY
  EXISTS on the source rows (bucket A). The aggregated row should carry the C++
  per-type *intensity-weighted* `jb` sum (the reducer `look03` is waiting for), not
  the raw per-type geometry. Recommend: add an aggregated `jb_T0` / `jb_T2_local_*`
  (sum over through-space sources) computed in the C++ reducer; retire
  `sum_dipolar_ringtype_*`.

### C2 — ring scalar dipolar sums `sum_dipolar_all`, `sum_dipolar_producer_valid`
- **Emit:** RingCurrentNeighborhood.cpp:105-106 (aggregated CSV); RecordSink.cpp:244-247.
- **Reassembled by:** `look01_ring_triangulate.py` (look01:11), `look02_self_vs_
  throughspace.py`, `look03_coefficient.py:31` (`K_valid = sum_dipolar_producer_valid`),
  `credibility_check.py`, `credibility2_instantaneous.py`,
  `static_environment_calibration.py` (via broad), `variance_decomposition.py`,
  `differencing_system_id.py`, `johnson_bovey_region_recovery.py`,
  `aimnet2_ceiling_ensemble.py`.
- **Why it's a smell:** the unweighted geometric sum is the Python OLS x-variable to
  fit a coefficient against DFT — the model is NOT asked the ring-current question;
  Python is. (`look03` regresses DFT σ_iso on this sum to "recover the constant.")
- **Model-blessed replacement:** the aggregated `jb_T0` (C++ literature-scaled,
  intensity baked in) is the well-posed scalar; `sum_dipolar_*` becomes redundant
  once the aggregated `jb` sum exists. NOTE for the lead: these are the most
  widely-consumed columns; cutting them is the biggest consumer-rewrite (≈10 scripts).
  A staged path is "add the aggregated `jb` sum, migrate the scripts, then retire."

### C3 — McConnell per-category dipolar sums `sum_dipolar_peptide_co/_peptide_cn/_sidechain_co/_aromatic`
- **Emit:** McConnellNeighborhood.cpp:101-104 (`kCatCols`, 4 columns aggregated CSV);
  computed McConnellNeighborhood.cpp:190-191.
- **Reassembled by:** `look03_coefficient.py`, `credibility2_instantaneous.py`.
- **Why it's a smell:** THE named counter-example in the brief. Per-bond-category raw
  geometric sums whose only purpose is for a Python script to multiply each by a
  per-category Δχ and add — reassembling the McConnell tensor outside the model. The
  typed `QtBond` already knows its category; the Δχ weighting belongs on the typed
  object, not a Python coefficient table.
- **Model-blessed replacement:** a per-category McConnell *Δχ-scaled* T2 kernel,
  computed C++-side from the typed bond (mirror `bondKernelT2FromSources`,
  BroadBackbone.cpp:79-105, which already builds the McConnell dipolar tensor from
  typed bond sources). Emit `mc_T0`/`mc_T2_local_*` per category (or summed) the way
  ring emits `jb_T*`; retire `sum_dipolar_<category>`. This is the exact
  McConnell-Δχ → `jb_T*` image the brief asks for.

### C4 — McConnell scalar dipolar sums `sum_dipolar_all`, `sum_dipolar_producer_valid`
- **Emit:** McConnellNeighborhood.cpp:99-100. (For mc, valid==all — the two are
  identical by construction, McConnellNeighborhood.cpp:196-197.)
- **Reassembled by:** same scalar-fit consumers as C2.
- **Replacement:** the C++ McConnell kernel scalar (C3's `mc_T0`). NOTE: emitting
  both `_all` and `_producer_valid` when they are provably equal is also a D-duplication.

### C5 — broad-backbone per-mechanism dipolar sums `ring_sum_dipolar`, `bond_sum_dipolar`
- **Emit:** BroadBackboneSink.cpp:211-212 (CSV); reduced in `ReduceBroadBackboneSources`,
  BroadBackbone.h:99-111.
- **Reassembled by:** `static_environment_calibration.py::load_broad_scalars`
  (static_environment_calibration.py:421-453), `variance_decomposition.py:392-393`,
  `differencing_system_id.py:67-73`, `aimnet2_ceiling_ensemble.py`.
- **Why it's a smell:** the consumer itself documents that these are not physics:
  *"The broad column is a bare ring sum, so k≈21 ppm·Å³ is not a direct fixed-
  coefficient comparison"* and *"the broad scalar column is untyped and unweighted,
  so bond Δχ is not directly comparable"* (static_environment_calibration.py:428-451).
  The interpretable comparison in the SAME script uses the C++ `*_literature_kernel_T2`
  sidecars (`calibrated-to-physics`, static_environment_calibration.py:499-523).
- **Model-blessed replacement:** ALREADY EXISTS — `ring_literature_kernel_T2`,
  `bond_literature_kernel_T2` (bucket A). The raw `*_sum_dipolar` scalars are the
  half-built version of those kernels. Recommend retiring them once the σ_iso scalar
  path is served by a C++ scalar projection (e.g. T0 of the literature kernel) rather
  than the raw geometric sum.

### C6 — broad-backbone Coulomb field scalars `field_local_{x,y,z}`, `field_z`, `field_mag`, `mu_local_{x,y,z}`
- **Emit:** BroadBackboneSink.cpp:214-218 (CSV) + `field_local` NPY
  BroadBackboneSink.cpp:243-245; reduced BroadBackbone.h:96-128 (`E = Σ q(r_atom−r_i)/r³`
  in C++).
- **Reassembled by:** `static_environment_calibration.py:455-461` (fits
  `A·field_z + B·field_mag²`), `variance_decomposition.py:394`,
  `buckingham_efield_t0.py`, `differencing_system_id.py:80-81`,
  `aimnet2_ceiling_ensemble.py`, `rediscover_capstone_charts.py`.
- **Why it's a smell (partial):** the FIELD itself is a legitimate C++-computed
  physical quantity (Coulomb field of typed charges) — closer to B than C. But
  `field_z`, `field_mag`, and the `field_mag²` Python derives are *pre-projected
  scalar features chosen so a Buckingham OLS can fit `A·E + B·E²`*. The polynomial
  combination is the Python-side physics (Buckingham expansion). The consumer flags
  it form-only with "no clean literature scalar coefficient"
  (static_environment_calibration.py:460).
- **Model-blessed replacement:** the `charge_literature_kernel_T2` (the C++ EFG-like
  T2, bucket A) already serves the T2 side. For the σ_iso/l=1 Buckingham side the
  proper emit is the C++ Buckingham *shielding* prediction (field·polarizability,
  scaled by a literature A) — but that coefficient is genuinely open
  (`REDISCOVERY_MAP.md`: buckingham stub ❌). RECOMMEND TO LEAD: keep `field_local`
  (B-grade raw field, equivariant input) and `mu_local` (the dipole, also raw);
  drop the pre-derived scalar reductions `field_z`/`field_mag` (Python can project
  the raw field onto ẑ and take a norm trivially, and *should*, since those are not
  physics — they are `.z()` and `.norm()` of an emitted vector, BroadBackbone.h:126-127).
  This is the cleanest C-cut: the scalars are literally vector accessors.

### C7 — `dipolar_3cos2m1_over_r3` (per-source, ring + mc + broad)
- **Emit:** RingCurrentNeighborhood.cpp:64, McConnellNeighborhood.cpp:77,
  BroadBackboneSink.cpp source header; computed in each extract.
- **Reassembled by:** `look02_self_vs_throughspace.py:23,47-50` (sums it),
  `look01`, `look03`, the per-type/per-category sums (C1/C3) are sums OF this.
- **Why it's borderline:** `(3cos²θ−1)/r³` is the bare dipolar angular factor — a
  *pre-combined* geometric scalar. It is one algebraic step from the raw `r`, `cosθ`
  the model also emits (bucket B). It is the atom from which all the `sum_dipolar_*`
  reassembly is built.
- **Recommendation:** LOW priority. It is cheap and is a useful per-source diagnostic,
  but note for the lead that it is redundant with (`r`, `cos_theta`) and is the seed
  of the C1–C5 reassembly. If the per-type/per-category `sum_dipolar_*` are retired,
  the per-source `dipolar` loses its only reassembly role and could follow.

---

## (D) REDUNDANT / DUPLICATE — CONSOLIDATE

Same quantity under multiple names/scales, or columns no analysis script reads.

### D1 — DFT target T2 emitted up to 3× per row kind, ×2 frames
- The ring/mc/charge_dipole `RecordSink` emits the DFT target T2 as (a) CSV columns
  `dft_total_T2_*` (RecordSink.cpp:125), (b) a `_sources_target_T2.npy`, (c) an
  `_aggregated_target_T2.npy` (RecordSink.cpp:181-186), and the local-frame variant
  the same way → `target_T2` + `target_local_T2`, sources + aggregated = **4 NPYs
  + the CSV columns** per case. The aggregated NPY is the CSV column re-serialized.
- **Consolidate:** the lead's own carrier decision (`SURFACE_DESIGN.md`: "CSV carries
  scalars, NPYs carry arrays") implies the T2 should live in ONE place. The CSV
  `dft_total_T2_*` columns duplicate the `target_T2.npy`. NOTE: the broad sink already
  did the right thing (target ONCE, on the aggregated row + NPY, BroadBackboneSink.cpp:1-25)
  — but ring/mc still repeat the full ~50-column target on EVERY source row
  (RecordSink.cpp:225) AND in the NPY. The 828 MB charge_dipole bloat
  (BroadBackboneSink.cpp:2-5) is this exact duplication; the broad fix was applied
  there but NOT back-ported to `RecordSink`.

### D2 — bare-kernel emitted on both row kinds + as NPY
- `bare_T0/T1_*/T2_*` (12 cols) on the sources CSV (RecordSink.cpp:224) AND the
  aggregated CSV (RecordSink.cpp:248) AND as `_sources_bare_kernel_T2.npy` +
  `_aggregated_bare_kernel_T2.npy` (RecordSink.cpp:187-192). The aggregated value is
  identical to the per-source value (same `(atom,frame)` H5 read). Only the sources
  copy is needed.

### D3 — five hand-rolled identical `writeNpyF64` implementations
- RecordSink.cpp:53-93, BroadBackboneSink.cpp:45-85, EfgFeatureSink.cpp:27-67,
  BuckinghamEfieldSink.cpp:35-75, Aimnet2FeatureSink.cpp:32-75 (the last templated).
  Byte-for-byte the same NPY-v1 little-endian writer. **Consolidate** to one shared
  helper (a single `rediscover::WriteNpy`); this is pure DRY, no behavior change.

### D4 — five hand-rolled identical CSV identity/frame/target header blocks
- The identity+frame+DFT-target column block is re-spelled as a literal string in
  every bespoke sink (BroadBackboneSink.cpp:99-130, EfgFeatureSink.cpp:127-132,
  BuckinghamEfieldSink.cpp:77-89, Aimnet2FeatureSink.cpp:83-104) AND built
  programmatically in `IdentityColumns()`/`TargetColumns()` (ExtractionSupport.cpp:70-89,
  164-185) for the `RecordSink` family. Two sources of truth for the same schema; a
  drift risk (the bespoke sinks already diverge: some include `dia/para`, some don't;
  EFG omits the raw 3×3 and T1). **Consolidate** the shared blocks onto the
  `ExtractionSupport` helpers and have the bespoke sinks reuse them.

### D5 — `frame_variant` emitted twice per row
- Both `RecordSink` (`writeIdentity`, RecordSink.cpp:107) and every bespoke sink
  (e.g. BroadBackboneSink.cpp:204 then again :209; EfgFeatureSink.cpp:153 then :156)
  write `frame_variant` / `backbone_frame_class` AND a second `frame_variant` column.
  Two columns, same value. Drop one.

### D6 — carried-but-unread DFT diagnostics
- `dft_dia_raw_*` (9), `dft_para_raw_*` (9), `dft_dia_iso`, `dft_para_iso`
  (ExtractionSupport.cpp:172-179): no read script consumes them (grep: zero hits).
  Physically honest to keep (dia/para decomposition is real), but they are 20
  unread columns on every ring/mc/charge_dipole row. NOTE for the lead — keep or
  move to an opt-in diagnostic NPY.

### D7 — `n_sources_valid` / `sum_dipolar_producer_valid` == their `_all` twin (mc)
- For McConnell there is no self/bonded concept, so `n_sources == n_sources_valid`
  and `sum_dipolar_all == sum_dipolar_producer_valid` always
  (McConnellNeighborhood.cpp:196-197). Two pairs of identical columns. (For ring they
  genuinely differ.)

---

## Prioritized cut / consolidate / replace list (recommend; do not change)

Ordered by physics-integrity payoff first, then cheap DRY wins.

1. **C3 (McConnell per-category sums) → C++ per-category McConnell-Δχ kernel.**
   Highest-integrity: it is the brief's named counter-example, and the C++ machinery
   already exists (`bondKernelT2FromSources`). Emitting `mc_T0`/`mc_T2_local_*` per
   category puts the Δχ on the typed `QtBond`, completing the McConnell-Δχ → `jb_T*`
   image. Replaces C3 + C4.

2. **C1 (ring per-type sums) → aggregated C++ `jb` sum.** The source-row `jb_T*`
   already exist; add the intensity-weighted aggregated `jb` the reducer-comment in
   `look03_coefficient.py:86-87` is explicitly waiting for. Replaces C1; lets C2 retire.

3. **C6 field scalars → drop `field_z`/`field_mag`, keep raw `field_local`/`mu_local`.**
   Cheapest C-cut with real meaning: the scalars are literally `.z()` and `.norm()`
   of an emitted vector (BroadBackbone.h:126-127); Python projecting an emitted vector
   is fine and is not reassembly. Keep the vector (equivariant raw input, B-grade).

4. **C5 (broad ring/bond sum_dipolar) → retire in favor of the existing
   `*_literature_kernel_T2`.** The blessed kernels already exist and are what the
   consumer calls `calibrated-to-physics`; the raw sums are the form-only half the
   same script disclaims. Provide a C++ T0 projection for the σ_iso scalar path so
   the scalars aren't needed.

5. **C2 (ring scalar sums) → retire after steps 1-2 land + scripts migrate.** Largest
   consumer surface (~10 scripts); stage it: add the aggregated `jb` sum, migrate,
   then cut. Do NOT cut before the replacement is consumed.

6. **D1 + D2 (target/bare-kernel duplication) → back-port the broad sink's
   target-once fix to `RecordSink`.** The fix exists (BroadBackboneSink) but ring/mc
   still repeat the ~50-col target per source row + duplicate the NPY. Disk + clarity.

7. **D3 + D4 + D5 (DRY) → one shared `WriteNpy` + shared header/column builders;
   drop the doubled `frame_variant`.** Pure cleanup, no behavior change; also removes
   the schema-drift risk in D4. Naturally folds into the `#29` engine/sink unification
   already on the backlog.

8. **C7 (`dipolar`), D6 (dia/para), D7 (mc duplicate pairs) → low priority.** Cheap,
   honest, or trivially redundant; revisit only if the C1–C5 reassembly retires (C7
   then loses its role) or disk pressure bites (D6).

---

## Cross-cutting observations for the lead

- **The smell is converging, not spreading.** The newest path (broad_backbone) emits
  the blessed `*_literature_kernel_T2` AND still carries the raw `*_sum_dipolar`
  scalars side-by-side — the surface is mid-transition from C to A, exactly as
  `REDISCOVERY_MAP.md` step 2 describes the per-cell workflow. The C-bucket is the
  un-finished tail of that transition, not new debt.
- **`charge_quadrupole` / `larsen_hbond` / `aimnet2_embedding` are fail-loud stubs
  (main_extract.cpp:65-95) — they emit nothing.** When implemented, design them
  C++-blessed from the start (a typed quadrupole/H-bond T2 kernel), so they never
  acquire a `sum_*` scratch column. This audit is the moment to set that rule.
- **No IUPAC/topology archaeology touched.** The typed-atom anchoring uses the
  collision-safe `selectUnique(scope, TypedAtomSelector)` (RingCurrentNeighborhood.cpp:177-195,
  ChargeDipoleNeighborhood.cpp:32-66) — the IUPAC-trap-safe path, not positional
  indexing. Noted and moved on per instruction.
- **SDK impact of any C-cut is small:** the SDK registers only NPYs, not CSV columns
  (python/nmr_extract/_catalog.py:196+). Cutting `sum_dipolar_*` (CSV-only) needs no
  SDK change; cutting/adding an NPY needs one `ArraySpec` per file.
