# All-Atom Equivariant Emit — Pre-MOPAC Review

Scope: READ-ONLY review of commit `2ecbe59` ("Add all-atom equivariant rediscover
emit") on branch `h5-reader-pysr-spike`. This review gates the MOPAC-add (#51): it
assesses the all-atom emit foundation + the functional API before MOPAC welds more
weight on. Verdicts are reserved for the lead; what follows is evidence + a concrete
prescription. No git archaeology was performed. IUPAC noted and skipped.

Files read in full: `AllAtomEquivariant.{h,cpp}`, `AllAtomEquivariantSink.{h,cpp}`,
`SpatialIndexSet.{h,cpp}`, `main_extract.cpp` diff, `Catalog.{h,cpp}`, `Verbs.{h,cpp}`,
`Relationship.h`, `ExtractionSupport.cpp` (`BuildTarget`), `LocalFrameBasis.h`
(`LocalFrame`), the new test, `_catalog.py` new block, `STATE.md` top,
`analysis/PATTERNS.md`, `STAGE1_MUTANT_LANDSCAPE.md` (EFG/Coulomb/MOPAC sections),
producer `Mopac*TrajectoryResult` + `MopacCoulombResult.h`, `QtTrajectoryH5.{h,cpp}`
MOPAC surface.

---

## 1. C++ model discipline — VERDICT: CLEAN (typed spine respected; no string-dispatch on identity)

The all-atom generalization rides the typed object model. No identity is reconstructed
from names or category strings; every category/identity decision reads a typed property.

- **Typed dispatch, not string-dispatch on identity.** `makeRingSource`
  (`AllAtomEquivariant.cpp:192`) reads `ring.TypeIndex()`, `ring.TypeIndexAsInt()`,
  `ring.Aromaticity()`, `ring.RingSizeValue()`, `ring.NitrogenCount()`, `ring.IsFused()`,
  `ring.LiteratureIntensity()`, `ring.JohnsonBoveyLobeOffset()` — the `QtRing` virtuals.
  `makeBondSource` (`:235`) reads `b.category`, `b.order`, `b.atomIndexA/B` off the typed
  `QtBond`. The `bondCategoryName()`/`ringTypeName()` switches (`:64`, `:78`) are
  enum→label *string formatting for output columns*, not identity dispatch — the
  category was already determined by the typed object upstream. This is the allowed
  polymorphism per `Relationship.h:20-24`. No `if (name == "PHE")` anywhere.

- **Self/bonded is chemistry, not position-string.** `ringContainsAtom` tests the
  target's *heavy parent* (`heavyParent`, `:93`) against `ring.atomIndices`
  (`:230-231`) — membership via the typed ring, the same heavy-parent rule the verbs
  use (`Verbs.cpp:96`). Bond self/bonded (`:276-279`) is endpoint-identity OR a
  recorded near-field ratio. Both are typed/geometric, not name scans.

- **Frame discipline is genuinely raw lab-frame.** `noLocalFrame` is a
  default-constructed `LocalFrame` with `is_valid=false` (`LocalFrameBasis.h:69`).
  `BuildTarget` (`ExtractionSupport.cpp:63`) gates `TensorToLocal` on `frame.is_valid`,
  so with `noLocalFrame` the DFT `total_raw` is emitted unrotated in the ORCA/H5 frame.
  The "no imposed per-atom local frame" claim holds in code, and the test asserts the
  negative (`!contains("target_local")`, `!contains("disp_local")`,
  `tests/...:95,102`). This is the corrected-substrate intent, faithfully realized.

- **KD cloud is the right extension, not a hack.** `AllBondMidpoints`
  (`SpatialIndexSet.h:26`) is a sibling cloud built in the same day-1 ctor loop as the
  anisotropic `BondMidpoints`; the diff (`SpatialIndexSet.cpp:92-112`) refactors the
  single bond loop to populate both clouds from one pass, so `BondMidpoints` (the
  legacy McConnell-subset callers) is byte-for-byte preserved while `AllBondMidpoints`
  enumerates every category. `entity_index` carries the topology bond index, so the
  consumer re-reads the typed `QtBond`. Clean.

One minor note, not a discipline breach: `heavyParent` is duplicated in both
`Verbs.cpp` (anon-ns, `:96`) and `AllAtomEquivariant.cpp` (anon-ns, `:93`). Trivial,
but it is exactly the kind of copy `PATTERNS.md` rule 8 ("debt in the exemplar
propagates") warns about; fold it onto the verb surface (`verbs::heavyParent`) so the
all-atom case and the engine share one definition. See §2 prescription.

---

## 2. Functional API coherence — DIAGNOSIS: the all-atom case is the SECOND sibling that bypasses the engine; this is the moment to unify, not after MOPAC

### What the emit actually is

`RunAllAtomEquivariantEmit` (`AllAtomEquivariant.cpp:363`) is a **hand-written
double loop** (`for row : dftRows()` × `for atom : atomCount()`) that calls the
Layer-1 verbs directly (`verbs::near`, `verbs::pos`, `verbs::ringGeom`,
`body.catalog.*`) and a stack of `makeXSource` free functions, writing straight into
`AllAtomEquivariantSink`. It does NOT go through `Relationship` / the iterated-closure
engine (`RelationshipEngine::RunRelationship`). It is the same architectural shape as
`RunBroadBackbone` — a sibling runner with its own carrier sink — which
`BroadBackbone.h:15-21` and `PATTERNS.md` rule 1 (`#29`) already flag as the latent
narrowness: the engine's `RunRelationship` is coupled to the scalar-sum `RecordSink`,
so any carrier that isn't scalar-sum gets its own driver.

So there are now **three** drivers over the same (atom, frame, sources) traversal:

1. `RelationshipEngine::RunRelationship` — ring/mc/charge_dipole, scalar-sum `RecordSink`.
2. `RunBroadBackbone` — heterogeneous two-kind carrier → `BroadBackboneSink`.
3. `RunAllAtomEquivariantEmit` — all-atom raw-source carrier → `AllAtomEquivariantSink`.

### Is it spaghetti? Not yet — but it is at the inflection point

The all-atom loop is, in isolation, *clean and readable*: it reuses the verbs and the
day-1 spine, it does not regenerate topology, it does not open a parallel H5 path, and
the per-mechanism blocks are flat and obvious. It is good one-off code. The problem is
structural, not local: the substrate now has **three independent traversals of the
same index space, each re-implementing "for each DFT frame, for each atom in the
stratum, find near sources, attach typed fields"**, and each with its own
present/finite gating, its own stats struct, and its own sink. The shared logic
(stratum iteration, frame→original-index mapping, near-cloud queries, present-checks)
is copy-evolved across them, and the `makeXSource` helpers in the all-atom file
duplicate geometry already expressed as attachers in the ring/mc relationships
(disp/r/inv_r3/cos_theta/dipolar, ring-normal orientation, bond-axis orientation).

That is the accretion `PATTERNS.md` rule 1 names: *"don't let it sprawl into one-off
walks or sibling runners (the broad case forced one — #29 unifies them)."* The
all-atom case is the second forced sibling. MOPAC-add (#51) will push attachers into
all three (a MOPAC charge cloud, MOPAC EFG/Coulomb-shielding source rows), so without
unification MOPAC triples the surface that drifts.

### Concrete prescription (do this BEFORE MOPAC)

The goal is not a plugin framework (that violates `feedback_no_abstractions`); it is to
make the **carrier** the only thing that varies, and the traversal the one thing that
doesn't. Minimal-and-clarifying, designed across all three shapes at once
(`feedback_functional_api_minimal_clarifying_abstraction`):

1. **Generalize the engine's sink seam, not its closures.** `RunRelationship` already
   owns the canonical traversal (stratum → frame_fn → selectors flattened → attachers →
   classifier → filter). The coupling is only that its terminal step calls
   `RecordSink::WriteAggregatedRow`. Introduce one carrier interface the engine writes
   the *attached `SourceSlot` set* into per (atom, frame) — `RecordSink` (scalar-sum
   reducer), `BroadBackboneSink` (broad reducer), and `AllAtomEquivariantSink` (raw
   per-source rows, no reduce) each implement it. The reducer becomes a property of the
   carrier, not of the engine. This is the `#29` unification, now with three concrete
   carriers to design against instead of two — which is *better*, because the all-atom
   "emit raw sources, don't reduce" carrier is the cleanest possible proof the seam
   isn't reducer-specific.

2. **Express the all-atom case as a `Relationship` with `reducer = none`.** Its stratum
   is `atomsWhere([](const QtAtom&){ return true; })` (all atoms); its `selectors` is
   the LIST it already uses inline (`nearBackend(RingCenters, ring_cutoff)`,
   `nearBackend(AllBondMidpoints, bond_cutoff)`, `nearBackend(ChargeSites,…)`,
   `nearBackend(Atoms,…)` for AIMNet2 charge); its attachers are the `makeXSource`
   bodies refactored into `Attacher` closures that write into the `SourceSlot`. The
   per-target APBS/AIMNet2-feature "source rows" (the synthetic
   `makeVectorFeatureSource`/`makeTensorFeatureSource` rows, `:317`, `:339`) become a
   target-feature attacher, not a source. Then the all-atom case is the SAME four-line
   loop the trajectory-scope philosophy demands, and MOPAC adds one selector + one
   attacher to ONE place.

3. **Fold the duplicated geometry + `heavyParent` onto the verb/attacher surface.**
   disp/r/inv_r3/cos_theta/dipolar is computed three times now (ring attacher, bond
   attacher, all-atom `makeXSource`). Make it a Layer-1 verb
   (`verbs::displacement(body, target, sourcePoint, frame) -> {disp,r,inv_r3,…}`) or a
   shared attacher, so the McConnell self/bonded cutoff fix and the near-field ratio
   live in one definition. Same for `heavyParent`.

If unification is too large to land before MOPAC under the current capacity constraint
(`feedback_resource_constraint`), the honest fallback is: **do NOT add a fourth sibling
runner for MOPAC.** Add MOPAC as selectors+attachers to the all-atom case as it stands
(it is the right carrier for the law-discovery substrate), and book the engine-seam
unification as the explicit `#29` debt with these three carriers as its design inputs.
What must not happen is a `RunMopacAllAtomEmit` fourth driver — that is the spaghetti
the lead should refuse.

---

## 3. MOPAC-on-the-floor inventory — the field/EFG leg IS on the floor

### The Stage-1 mandate (why this section is load-bearing)

`STAGE1_MUTANT_LANDSCAPE.md` is explicit (`:139-140`, `:359-360`, `:633`, `STATE.md:21`):
in the backbone multi-source block, **`efg_coulomb` and `mopac_efg` have *moderate*
single-R²; `efg_apbs` is *tiny* (median ~0.0007–0.0025, zero drop-delta counts).** The
field/EFG signal lives in the Coulomb/MOPAC leg, not APBS. The all-atom emit currently
emits APBS EFG/E-field (the tiny leg) and emits **nothing** from the MOPAC/Coulomb leg.
So the emit captures the field channel that Stage-1 found least load-bearing and leaves
the moderate one on the floor. This is the priority the MOPAC-add must fix.

### What the emit captures today

Target-axis features actually written (`AllAtomEquivariantSink.h:38-50`,
`Sink.cpp` columns): `apbs_efield` (V/Å), `apbs_efg_T2` (V/Å²), `aimnet2_charge`,
`aimnet2_charge_response_gradient` (+scalar), `aimnet2_embedding` (256-d). Source-axis:
ring rows, all-bond-category rows, FF14SB charge-site rows, AIMNet2 charge-site rows
(`Atoms` cloud). **No MOPAC anything.**

### Every `mopac_*` in the dense substrate, and where it lands

The dense `.LGS` is `1p9j-calibration-dense-mopac-live-orca` — a FullFat
(`--mopac`) extraction, so the MOPAC family exists. But it splits into two storage
classes, and only one is even readable by the reader's trajectory H5 path:

**A. In the trajectory H5 as per-frame TRs — readable by the reader, NOT yet in the rediscover Catalog:**

| MOPAC output | What it is | Reader access | In rediscover Catalog? | On the floor? |
|---|---|---|---|---|
| `mopac_coulomb_shielding_time_series` | **T2 shielding from MOPAC-charge Coulomb EFG** (`MopacCoulombResult.h`: V/Å² EFG → shielding contribution). THIS is the moderate Stage-1 `efg_coulomb`/`mopac_efg` leg, as a shielding T2. | `QtTrajectoryH5::mopacCoulombShielding()` (`QtTrajectoryH5.h:89`, `QtT2TimeSeries`) | **NO** — no `ArrayId` exists | **YES — the priority leg** |
| `mopac_mc_shielding_time_series` | MOPAC-charge McConnell bond-anisotropy shielding | `mopacMcShielding()` (`:90`) | NO | YES |
| `mopac_vs_ff14sb_reconciliation` | per-atom signed cos(MOPAC-EFG-T2, FF14SB-EFG-T2), the FullFat charge-source probe | `mopacVsFf14sbReconciliation()` (`:91`) | NO | YES (diagnostic — see note) |
| `mopac_charge_welford` | per-atom mean/var of MOPAC Mulliken charge over frames | `mopacChargeWelford()` (`:150`) | NO | Partial (rollup, not per-frame) |
| `mopac_bond_order_welford` | per-atom mean/var Wiberg bond order | `mopacBondOrderWelford()` (`:151`) | NO | Partial (rollup) |

**B. Per-atom NPY surface only (SDK `_catalog.py:355-381`) — NOT written to the trajectory H5 as per-frame TRs, so NOT reachable on this dense-trajectory substrate:**

`mopac_charges`, `mopac_scalars`, `mopac_bond_orders`, `mopac_global`,
`mopac_coulomb_shielding`(9), `mopac_coulomb_E`(3), **`mopac_coulomb_efg_backbone`(5),
`mopac_coulomb_efg_aromatic`(5)**, `mopac_coulomb_scalars`, `mopac_mc_shielding`(9),
`mopac_mc_category_T2`(25), `mopac_mc_scalars`.

The crucial subtlety for the MOPAC-add: the *raw MOPAC Coulomb EFG tensor*
(`mopac_coulomb_efg_backbone/aromatic`, the literal V/Å² EFG) is a **per-atom-NPY** output,
not a trajectory TR. On the dense-trajectory substrate the only MOPAC field/EFG signal
physically present in the H5 is `mopac_coulomb_shielding_time_series` — the EFG already
contracted to a shielding T2. That contracted T2 *is* the moderate Stage-1 leg
(`mopac_efg` was assessed via its shielding/T2 footprint), so capturing it is the right
move — but the MOPAC-add brief should state plainly that it is capturing the
**MOPAC-Coulomb-EFG-derived shielding T2**, not the raw EFG tensor, because the raw EFG
is not on this substrate. If the law-discovery model wants the raw MOPAC EFG tensor (to
fit EFG→shielding rather than ingest a pre-contracted shielding), that requires a
producer-side change to write a `mopac_coulomb_efg_time_series` TR — which is a producer
modification, off-limits during this work and a separate scheduled item. Flag it; do not
silently substitute the shielding T2 and call it "the EFG."

### Note on the reconciliation probe

`mopac_vs_ff14sb_reconciliation` is the FullFat probe's *output* (a similarity scalar),
not a feature. It is legitimately diagnostic, not training material. It need not be in
the law-discovery feature substrate, but it IS the on-disk evidence that the MOPAC and
FF14SB EFG T2 vectors diverge — worth surfacing as a per-target provenance/QC column so
a downstream analysis can condition on "did MOPAC and FF14SB agree here," rather than
leaving it unread.

### MOPAC-add concrete requirements

1. **Capture the priority leg:** add `ArrayId::MopacCoulombShielding` (T2) wired to
   `mopacCoulombShielding()` in `Catalog.cpp`, emitted as a per-target source/feature
   row in the all-atom case (mirror the `apbs_efg` block, `:507`). This is the moderate
   Stage-1 field/EFG signal and must not stay on the floor.
2. **Capture `mopac_mc_shielding` (T2)** similarly — the MOPAC-McConnell variant Stage-1
   tracks alongside classical McConnell.
3. **MOPAC charge cloud:** add `ArrayId::MopacCharge` *actually wired* (today it is a
   stub: `Catalog.cpp:96` residence `Absent`, `present()` returns `false` at `:144`,
   `value()` returns `0.0` at `:190`). On the dense MOPAC substrate the per-frame MOPAC
   charge is the natural third `charge_source` alongside ff14sb/aimnet2 in the existing
   charge-site loop (`:463`, `:480`) — but verify whether per-frame MOPAC charge is in
   the trajectory H5 as a TR, or only as the Welford rollup. If only the rollup exists,
   the per-frame `mopac` charge-site rows are NOT available here and the brief must say so
   (capture the Welford mean as a static charge source, or note the gap).
4. **Do not double the 57 GB** (STATE.md NEXT): MOPAC rows extend the existing
   sources/targets carrier; re-emit-with-drop-old, not a parallel substrate.

---

## 4. Provenance completeness — the emit added irreps/units/sign/rank/parity/mechanism; FOUR Stage-1-mandated provenance fields are still missing or thin

`STAGE1_MUTANT_LANDSCAPE.md`'s mining warns: preserve exact source/family,
tensor-component, target, normalization, **atom-role**, **named-atom**, and **SUPPORT
(nonzero/rank)** metadata, or downstream recreates the confounds. The new ArraySpecs
(`_catalog.py:242-259`) added irreps/units/sign_convention/tensor_rank/parity/mechanism
— good, and consistent with the existing rediscover specs. What the mining demands that
is still missing:

- **Atom-role / named-atom: present in the CSV, ABSENT from the NPY sidecar provenance.**
  The CSV targets carry `atom_name`, `element_ord`, `amino_acid_ord`, `residue_*`, and
  the sources carry `source_atom_name`, `source_element_ord`, `category` — that part is
  good and is the Stage-1 named-atom/atom-role requirement met *in the CSV*. But the NPY
  sidecars (target_T2, target_raw, apbs_*, aimnet2_*) are bare value arrays whose only
  identity is *row order* aligned to the targets CSV. There is no `native_axis` stronger
  than `rediscover_target_row` and no atom-role/stratum column carried INTO the sidecar
  metadata. A consumer who loads `all_atom_equivariant_target_T2.npy` alone cannot
  recover atom-role without joining the CSV by row. That join is implicit and unstated.
  **Prescription:** document the row-alignment contract explicitly (sidecar row i ==
  targets.csv row_id i), and consider an `atom_role`/`stratum` and `element_ord` column
  echoed into a companion identity NPY, so the per-stratum requirement (`PATTERNS.md`
  rule 7) is mechanizable without the CSV.

- **Support / nonzero / rank metadata: MISSING.** Stage-1's central warning is that
  apparent EFG/field signal collapses under blocking and is *support-sensitive*
  (thin-N). The emit writes presence flags per row (`apbs_efg_present`,
  `aimnet2_charge_present`, etc.) — good, that is per-row support. But there is no
  aggregate support/rank provenance: no per-feature nonzero count, no per-stratum
  effective-N, no rank/condition of the T2 payload. The `STATE.md` ring finding
  ("effective signal-N ~2.2 atoms even in the 41-atom pool") is exactly the confound
  that support metadata exists to prevent downstream from missing. **Prescription:** the
  stats struct (`AllAtomEquivariantStats`) already counts ring/bond/charge rows and
  per-type/per-category rows — emit those counts into the manifest as per-feature
  support, and add a per-stratum present-count so a fitter can report effective-N
  without re-deriving it (which would be a Python recompute, forbidden by rule 5).

- **Normalization-variant: MISSING and this is a known Stage-1 confounder.** Stage-1
  found EFG visibility is "tensor-specific and normalisation-sensitive" (`:31`), and
  within-protein-Z vs raw changes which mechanism reads strong (`:499-500`). The emit
  writes raw lab-frame values with units, which is correct and primary — but it records
  no normalization-variant tag. Since the substrate is raw (good — normalization is a
  fitter choice), the requirement is lighter: **document in the ArraySpec/manifest that
  these are RAW un-normalized lab-frame values**, so a downstream Z-scoring is an
  explicit recorded transform, not an ambiguity. A `normalization="raw"` note on the
  feature specs would close it.

- **Minor: `mechanism` divergence.** The new
  `all_atom_equivariant_aimnet2_charge_response_gradient` uses `mechanism="aimnet2"`
  (`:255`) while the existing per-atom `aimnet2_charge_response_gradient` uses
  `mechanism="charges"` (`:470`). Both use `irreps="1o"`/`parity="odd"` consistently
  (so no irrep regression — the `1e`+`parity=odd` E-field entries are a pre-existing
  codebase convention, not introduced here). Pick one `mechanism` grouping for the
  AIMNet2-CRG family so the `mechanism`-keyed R analysis (`PATTERNS.md` rule 7 /
  `stage1_bmrb_dimension_independence.R`) groups them together.

The provenance the MOPAC-add must carry forward: when MOPAC rows land, they need the
same source/family (`charge_source="mopac"`, `mechanism`), the EFG-vs-shielding
distinction made explicit (§3 — it's the contracted shielding T2, not raw EFG), and the
`mopac_vs_ff14sb_reconciliation` similarity surfaced as the charge-source-divergence QC
column.

---

## Bottom line for the lead

- **Model discipline: clean.** Typed-spine, virtual dispatch on `QtRing`/`QtBond`, raw
  lab-frame target verified, `AllBondMidpoints` is a proper day-1 sibling cloud. No
  string-dispatch on identity, no model bypass, no parallel H5. Ship-quality on this axis.

- **Functional API: at the inflection point, not yet spaghetti — straighten it FIRST.**
  The all-atom case is the SECOND forced sibling runner (after broad_backbone) that
  bypasses the engine because the engine is coupled to the scalar-sum sink. MOPAC will
  push attachers into all three drivers. The prescription is the `#29` engine-seam
  unification (carrier owns the reducer; all-atom is the reducer-free carrier that proves
  the seam) PLUS folding the triplicated disp/geometry + `heavyParent` onto the verb
  surface. **If capacity won't allow unification before MOPAC, the hard floor is: add
  MOPAC to the existing all-atom case, never a fourth `RunMopac*Emit` driver.**

- **MOPAC on the floor: the moderate field/EFG leg is entirely unread.** The emit
  captures APBS EFG (Stage-1-tiny) and none of MOPAC. The reader already exposes
  `mopacCoulombShielding()` (the MOPAC-Coulomb-EFG-derived T2 shielding = the moderate
  Stage-1 leg) and `mopacMcShielding()`, but the rediscover Catalog has zero MOPAC
  ArrayIds (the lone `MopacCharge` is a hard-coded stub). The raw MOPAC EFG *tensor* is
  per-atom-NPY only and NOT on this trajectory substrate — capture the shielding T2 and
  label it honestly; do not call it "the EFG."

- **Provenance: irreps/units/sign/rank/parity/mechanism are in; atom-role-in-sidecar,
  support/effective-N, and normalization-variant are the gaps** Stage-1 explicitly warns
  recreate the confounds. None block the emit, but the MOPAC-add should fill them as it
  extends, and the row-alignment (sidecar row == CSV row_id) contract should be written
  down.
