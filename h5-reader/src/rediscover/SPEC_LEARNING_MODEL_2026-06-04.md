# Stage-3 Learning Model Spec - 2026-06-04

Status: **DESIGNED, not built**; companion references trued 2026-06-04. This document specifies the Stage-3
"ring-toss" equivariant-conditioned GNN. It is the architecture target for a
later build after the chewer and 720-WT static ingest exist. It does not
authorize training, extraction, or code changes.

## Grounding

Load-bearing local references:

- The rediscover substrate is method-agnostic because the problem is a
  permutation-invariant source set, rotation-equivariant, T2-valued local sum;
  the source set plus displacement vectors and target tensor are the natural
  equivariant sum-pooling input (`DESIGN.md:30-42`, `GUIDANCE.md:57-64`).
- Physics stays in the typed C++ model. Python consumes emitted arrays and fits;
  missing quantities are fixed by extending C++ emit/query surfaces, not by
  recomputing geometry or kernels in Python (`GUIDANCE.md:24-33`,
  `GUIDANCE.md:54-56`, `PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md:10-22`).
- The current all-atom equivariant foundation emits raw lab/molecular-frame
  geometry and targets with **no imposed per-atom local frame**; T2 is stored as
  five library-basis components and converted with frozen `get_C()`
  (`STATE.md:292-296`, `PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md:54-62`).
- The e3nn exemplar is the model template: e3nn spherical harmonics /
  `1o x 1o -> 2e` tensor products, scatter-add source pooling, C++-emitted
  inputs/targets only, and no Python projection fallback
  (`analysis/equiv_t2_e3nn.py:1-16`, `analysis/equiv_t2_e3nn.py:18-31`,
  `analysis/equiv_t2_e3nn.py:104-141`).
- The clean protocol is blocked/purged held-out frames, train-only centering and
  normalization (`analysis/e3nn_protocol.py:7-11`,
  `analysis/e3nn_protocol.py:14-30`, `analysis/e3nn_protocol.py:49-80`).
- The physics-architecture review independently re-derived the shared angular,
  per-type radial picture from the single `D_ab` object
  (`PHYSICS_ARCHITECTURE_UNIFICATION_2026-06-04.md:10-17`).

The companion specs now exist: `SPEC_720_STATIC_INGEST_2026-06-04.md` and
`SPEC_INTERP_1P9J_SCOPING_2026-06-04.md`. The bounded interpolation build also
landed in `BUILD_INTERP_RESULT_2026-06-04.md`; treat it as the v0/down-payment
result, not the full Stage-3 chewer.

## Non-Negotiable Boundary

This model is **e3nn**. There is no hand-rolled spherical-harmonic, Wigner-D,
tensor-product, or second equivariant model in the training path. The only
basis conversion in Python is the frozen constant `change_of_basis.get_C()`;
that is a convention map on emitted 5-vectors, not a tensor projection or
physics recomputation.

T2 is sacred. The primary target and prediction are the total shielding T2
five-vector. `sigma_iso` / T0 is emitted and trained/reported as a scalar
comparison head, never as a replacement for T2.

There is **no imposed per-atom local frame** for Stage 3. Inputs are raw
molecular-frame vectors/tensors from the C++ model; e3nn handles rotation.
Physical source axes are allowed because they are part of the source geometry
itself, for example ring normals, bond axes, hydrogen-bond partner axes, and
field vectors. They are not canonical frames attached to the target atom.

There is **no protein or spatial model in Python**, not even a secondary
aggregate one. Python does not open `trajectory.h5`, build atoms/residues,
construct KD-trees, discover rings/bonds/H-bonds, compute distances, decide
self/bonded exclusions, or assemble kernels. The chewer is the dependency that
binds Python to the live C++ model through pybind11 and returns vectorized
source/target tensors. If a feature is missing, the fix is a C++ emit/query
extension.

No terabyte Python dump. The default training substrate is lean row-aligned
arrays plus vectorized source-query tensors. Full pair/source materialization is
not the default product; the existing per-atom emitter explicitly writes one
row per `(DFT frame slot, atom)` and never emits a default source CSV
(`PerAtomSubstrate.h:1-5`).

## Thesis Role

Stages 1-2 are physics explanation: the question is whether named classical
kernels carry signal, whether coefficients are determinable, and whether the
one-object `D_ab` picture survives adversarial reading. R2 is diagnostic there,
not the thesis by itself.

Stage 3 is prediction. R2 **is** the metric, because the goal changes from
"does this kernel explain the DFT modulation?" to "can an equivariant,
T2-valued model predict held-out shielding from geometry and allowed
conditioners?" The earlier Stage-1 MLP rejection does not bind this path: that
was not an e3nn, permutation-invariant, source-sum, T2-valued model over raw
geometry. It was a different model class with the wrong inductive bias for the
source-set tensor problem. The substrate was deliberately kept open for this
path (`DESIGN.md:39-42`, `GUIDANCE.md:57-64`).

## Built / Designed / Open

**BUILT**

- All-atom raw equivariant emit for 1P9J: 846 atoms x 660 DFT rows, raw
  molecular-frame displacement/orientation vectors, no imposed per-atom frame,
  DFT raw/T2/T0 targets, and ArraySpec-style irreps/units/sign metadata
  (`STATE.md:292-296`).
- Per-atom aggregate substrate and sidecars:
  `per_atom_substrate_target_T2.npy`, `per_atom_substrate_target_T0.npy`,
  target dia/para/T1 diagnostics, classical features, ring paths, method paths,
  H-bond conditioning, generic conditioning, dominance, bins, driver modulation,
  and optional AIMNet2 embedding (`PerAtomSubstrate.cpp:3710-3733`).
- e3nn exemplars for T2 ring and EFG paths, frozen basis conversion, and the
  clean blocked/purged protocol.
- Stage-2 explanatory results: charge/ring/unified within-axis signals stand;
  true 1P9J between-axis recovery is null, so 720-WT owns transferability
  (`STATE.md:63-75`).

**DESIGNED HERE**

- The conditioned e3nn source-sum predictor over the variable source set.
- The source-query/chewer boundary that supplies per-source tensors to GPU while
  keeping spatial work in C++.
- The three-ground evaluation arc: train on DFT, transfer on 720-WT statics,
  shoot at BMRB experimental shifts in the dark.

**OPEN**

- The exact chewer API and whether training consumes a persisted tensor file or
  a pybind11-produced tensor cache. Either is acceptable if source selection and
  geometry remain C++-owned and vectorized.
- The exact source-type taxonomy and radial-network sharing groups for v1.
- The companion 720 ingest and 1P9J interpolation specs exist; interpolation v0
  has a result record, while the 720 ingest remains design-only.
- `OF3 embedding + OF3` and `boosts` are named in the brief but not grounded by
  any current column/file name found in this directory. They must be named by
  the lead or surfaced by a C++/producer spec before implementation.
- BMRB atom mapping, averaging policy, and uncertainty model.

**CAVEATED**

- BMRB is not clean validation for this MD run. Experimental shifts are ensemble
  averages; the short 1P9J trajectory samples a narrow geometry ensemble.
- AIMNet2 did not rescue current all-atom T2 fits after alpha-selection, but
  this is not a blanket ban on using AIMNet2 as a prediction conditioner. The
  stage changes from explanation to prediction.
- Field and H-bond standalone explanatory results are weak/null in current cuts,
  but their features can still condition a predictor.

## Input Object

Each training case is a target row:

```text
target case i = (protein_id, atom_index, frame_slot or static_pose_id)
target y_i = total DFT T2 in library basis, converted to e3nn 2e by get_C()
optional scalar target y_i_T0 = sigma_iso
source set S_i = variable set of C++-selected source objects around target i
```

For each source `s in S_i`, the chewer returns vectorized arrays:

```text
source_kind, source_type, source_category_ord
target_row_id / group id
disp = source_position - target_position        # 1o vector, molecular frame
r, inv_r3, cutoff/provenance flags
optional physical axis = ring normal / bond axis / H-bond axis / field vector
optional source scalar/tensor payloads emitted by the C++ model
```

This mirrors the existing `PairContribution` shape: mechanism/source identity,
`disp`, `r`, `inv_r3`, `cos_theta`, kernel T0/T2, contribution, and provenance
flags (`PerAtomSubstrate.h:54-77`). For Stage 3 the source tensor is not a
default CSV dump; it is a vectorized chewer/query result or a bounded persisted
training tensor produced from that same C++ path.

## Feature Pile

Prediction may use "understood" and "not-yet-understood" features. That is a
deliberate Stage-3 difference from Stage 1-2: the target is held-out predictive
R2, not explanation by a named kernel. The feature pile is still constrained by
the spine boundary: every feature must come from emitted C++/producer arrays or
from chewer queries over the live C++ model.

| Feature family | Current carrier / column anchor | Status for Stage 3 |
|---|---|---|
| DFT target | `per_atom_substrate_target_T2.npy`; T0 companion in `per_atom_substrate_target_T0.npy`; library T2 order converted by `get_C()` (`PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md:145-152`, `ALLATOM_FIT_SPEC_2026-06-03.md:221-234`) | **BUILT** for 1P9J; primary target is total T2 |
| Target identity | rows CSV fields: atom/residue/name/type/role/stratum/provenance (`PerAtomSubstrate.cpp:2870-2889`) | **BUILT**; categorical conditioners, not a Python protein model |
| Raw source geometry | chewer/query over C++ model; current all-atom foundation has raw `disp_*` / orientation vectors and no imposed frame (`STATE.md:292-296`) | **DESIGNED dependency** for GPU-scale source sets |
| Nearest geometry/counts | `nearest_ring_r`, `nearest_charge_r`, `nearest_bond_midpoint_r`, `nearest_heavy_atom_r`, ring/charge/bond count columns (`PerAtomSubstrate.cpp:2247-2254`) | **BUILT** aggregate conditioners |
| Dihedrals / local conformation | `dssp_chi_0..11`, `omega_actual`, `pyramidalization`, `dssp_ss8_*`, `ring_geometry_*` in `per_atom_substrate_features_hbond_conditioning.npy` (`PerAtomSubstrate.cpp:2163-2189`) | **BUILT** as conditioners; use only as emitted scalars/vectors |
| Motion / driver exercise | `sd_*_by_atom` in `per_atom_substrate_driver_modulation_by_atom.npy`; motion/quiet also appears in CaseHunter selection metadata (`PerAtomSubstrate.cpp:2277-2283`, `PerAtomSubstrate.cpp:2538-2548`, `STATE.md:82-84`) | **BUILT** as input-side conditioner; no target leakage |
| Ring-current shadows | `ring_jb_T0/T2`, source-level JB path; ring path sidecar has `bs_*`, `hm_*`, `ringchi_*`, `pq_*`, `disp_*` (`PerAtomSubstrate.cpp:2116-2132`, `PerAtomSubstrate.cpp:2192-2199`) | **BUILT/PARTIAL**; use JB as primary ring shadow, alternatives caveated by convention |
| Charge / field shadows | `charge_q_over_r3_T2_*`, `ff14sb_field_*`, `apbs_E_*`, `apbs_efg_T2_*`, `mopac_coulomb_shielding_T2_*` (`PerAtomSubstrate.cpp:2197-2216`, `PER_ATOM_SUBSTRATE_SPEC_2026-06-02.md:164-171`) | **BUILT/PARTIAL**; MOPAC field is important in the combine |
| McConnell / bond anisotropy | `mc_lit_T0_valid`, `mc_lit_T2_valid_*`, `mc_category_T2_*`, `mc_scalars_*`, `mopac_mc_*`, `mopac_bond_orders_*` (`PerAtomSubstrate.cpp:2135-2144`, `PerAtomSubstrate.cpp:2200-2209`) | **BUILT/PARTIAL**; standalone explanation weak, predictor may condition on it |
| H-bond geometry and Larsen | `hbond_T2_*`, nearest H-bond geometry, counts, donor/acceptor flags, `larsen_hbond_1pHB/2pHB/1pHaB/2pHaB_T2_*`, water term, DSSP backups (`PerAtomSubstrate.cpp:2163-2183`, `PerAtomSubstrate.cpp:2220-2226`) | **BUILT/CAVEATED**; Larsen ppm features must not be double-counted in explanation, but prediction may use them with labels |
| AIMNet2 | `aimnet2_charge`, `aimnet2_crg_scalar/x/y/z`, `aimnet2_efg*`, optional `per_atom_substrate_aimnet2_embedding.npy` (`PerAtomSubstrate.cpp:2151-2153`, `PerAtomSubstrate.cpp:2217-2219`, `PerAtomSubstrate.cpp:3232-3235`) | **BUILT/CAVEATED**; allowed for prediction despite weak current T2 lift |
| MOPAC | `mopac_coulomb_shielding_T2_*`, `mopac_mc_shielding_T2_*`, `mopac_coulomb_E_*`, `mopac_coulomb_efg_*`, `mopac_scalars_*`, `mopac_welford_mean_charge` (`PerAtomSubstrate.cpp:2141-2150`, `PerAtomSubstrate.cpp:2204-2209`, `PerAtomSubstrate.cpp:2883-2884`) | **BUILT**; distinguish MOPAC-Coulomb shielding from raw EFG |
| Other calculator outputs | water EFG/E-field, EEQ coordination, pi-quadrupole, dispersion, SASA, ring alternatives (`PerAtomSubstrate.cpp:2154-2159`, `PerAtomSubstrate.cpp:2227-2244`) | **BUILT** as optional conditioners/features |
| Dominance/isolation | `gap_to_2nd_*`, `dominant_fraction_*`, dominance bins (`PerAtomSubstrate.cpp:2089-2103`, `PerAtomSubstrate.cpp:2264-2275`, `PerAtomSubstrate.cpp:2313-2318`) | **BUILT**; primary for clean per-mechanism v1 scoping |
| OF3 embedding + OF3 | no local carrier found by name | **OPEN**; do not implement until lead identifies producer/source columns |
| Boosts | no local carrier found by name | **OPEN**; do not invent; require lead naming and C++/producer grounding |

## Architecture

The model is an e3nn-class equivariant source-sum network:

```text
for each target row i:
    for each source s in S_i:
        angular_2e_s = e3nn angular map(raw vectors / physical axes)
        radial_weights_s = per-type radial MLP(invariant source features,
                                               target/source conditioners)
        contribution_2e_s = radial_weights_s * angular_2e_s
    pred_T2_i = scatter_add_s(contribution_2e_s grouped by target i)
    pred_T0_i = optional scalar source-sum / target-conditioned scalar head
```

**Angular law.** The angular block is e3nn, not hand code. The required core is
`1o x 1o -> 2e` `o3.TensorProduct` or `FullyConnectedTensorProduct`, plus
e3nn spherical harmonics where a pure displacement or physical axis basis is
needed. This follows the exemplar, where e3nn computes `Y2(r_hat)`,
`Y2(n_hat)`, the `r_hat x n_hat -> 2e` cross term, and scatter-add pooling
(`analysis/equiv_t2_e3nn.py:13-16`, `analysis/equiv_t2_e3nn.py:119-141`).

**Shared angular, per-type radial.** The model embodies the one-kernel
`D_ab` insight for the `D_ab`-family mechanisms: charge/field, McConnell/bond
anisotropy, and geometric H-bond shadows differ mainly by source
axes/lists/intensities, not by inventing unrelated angular laws. The shared
angular tensor-product basis applies to that family. Ring-current enters through
its own e3nn angular maps, already present in the exemplar: `Y2(r_hat)`,
`Y2(n_hat)`, and the cross term, with JB as the primary ring shadow and
convention-caveated alternatives.
Radial/conditioner networks are per source type or per mechanism family:

```text
radial family examples:
ring type / ring calculator family
charge source family: FF14SB / AIMNet2 / MOPAC-Welford / static charge
bond category: peptide CO, peptide CN, sidechain, aromatic, ...
field/EFG source family: APBS, MOPAC-Coulomb, water, AIMNet2 EFG
H-bond/Larsen class
```

The radial MLP consumes invariant scalar quantities only: distance, cutoff flags,
source category, source scalar values, target/source type embeddings, and
allowed aggregate conditioners. It must not consume vector conditioners, target
values, residuals, or held-out-derived statistics. The failure modes are
specific: vector-fed radial gates break equivariance, and all-rows or
held-out-derived preprocessing leaks across the split. The exemplar keeps this
structural by feeding radial scalars such as `[r, ring_intensity]` and by using
train-only centering/normalization through `center_mask` (`analysis/equiv_t2_e3nn.py:115-118`,
`analysis/equiv_t2_e3nn.py:132-152`).

**Conditioning.** The feature pile conditions radial gates and optional scalar
heads. It is not a second spatial model. For example, AIMNet2 embedding, MOPAC
bond-order scalars, DSSP/dihedral features, dominance fractions, or driver
modulation may alter the radial response for a source family, but they do not
replace the C++ source set or geometry.

**Pooling.** Pooling is scatter-add over source rows keyed by target row. This
gives permutation invariance over the variable source set and preserves
equivariance because every per-source contribution is a proper e3nn irrep.

**Output.** Primary output is `1x2e` total shielding T2, five components in the
e3nn basis during training and mapped back to library basis for reporting if
needed. The scalar companion is `1x0e` `sigma_iso` / T0. Report scalar results
alongside the tensor model as comparison, never as a replacement.

**No local-frame branch.** There is no alternative architecture that first
rotates each target into a learned or hand-authored local target frame. Older
local-frame sidecars can remain compatibility/audit artifacts, but Stage 3 uses
raw molecular-frame vectors/tensors and e3nn equivariance.

## Three Grounds

| Ground | Role | Protocol | Caveat |
|---|---|---|---|
| DFT | **TRAIN**. This is the last stable ground: ORCA/r2SCAN shielding tensors, total T2 primary, T0 companion. | Held-out DFT splits: blocked/purged within frames for 1P9J, leave-atoms-out where meaningful, leave-protein/static-pose splits once 720-WT exists. | Do not train on dia/para split targets as the thesis target; total T2 is the stable target. |
| 720-WT | **TRANSFER**. This is the statics transferability/between-axis instrument. 1P9J true between-axis is null, so 720-WT owns this question (`STATE.md:63-75`, `NEXT_SESSION_PROMPT.md:43-51`). | Train on DFT source ground, test transfer on curated 720-WT static poses once `SPEC_720_STATIC_INGEST_2026-06-04.md` is built. Use protein/static-pose held-outs; no new ORCA in this spec. | Companion spec exists; exact input path and memory strategy remain open until ingest lands. |
| BMRB | **SHOOT IN THE DARK**. Experimental shifts are the thesis-facing dark target. | Freeze the DFT-trained model and compare predicted shifts to BMRB mappings without training or model selection on BMRB. | Not clean validation: BMRB is an ensemble average, and the short MD ensemble is incomplete. Interpolate-to-validate against held-out DFT: yes. Train on BMRB: no. |

## Protocol

Primary score is held-out R2 on the five-component total-T2 prediction. Also
report `|T2|` correlation and T0 R2, but do not use either as the primary T2
claim. Every table must distinguish tensor R2, invariant magnitude correlation,
and scalar T0 R2.

Splits:

- 1P9J within-axis: blocked/purged frame split using the clean e3nn protocol;
  train-only centering and normalization.
- 1P9J leave-atoms-out: report only where source support is not too thin; do
  not call it transferability after the true-LOAO audit.
- 720-WT: leave-protein, leave-mutant-family, or lead-approved static-pose
  splits once the ingest spec exists.
- BMRB: no training split because it is not a training ground; use only as a
  frozen-model dark shot.

Preprocessing rules:

- `get_C()` must be checked as orthogonal below `1e-12`; every T2 block used by
  the model has exactly five components and is transformed consistently
  (`ALLATOM_FIT_SPEC_2026-06-03.md:474-479`).
- Standardization, imputation, PCA, embedding compression, atom centering, and
  alpha/model selection are fit on training rows only.
- No partition, feature-selection rule, source-ranking rule, or favorable-case
  shortlist may use DFT target values, held-out residuals, or fitted
  coefficients (`ALLATOM_FIT_SPEC_2026-06-03.md:368-381`,
  `ALLATOM_FIT_SPEC_2026-06-03.md:490-496`).
- Every missing mechanism is represented by present flags plus NaN/masked
  payloads, not dropped rows.
- Thin support is flagged, not hidden.

## MSc-Realistic Scope

The full conditioned predictor is the forward arc, not the first build chunk.
The deliverable-sized version is:

1. Keep the per-mechanism dominance-clean equivariant fits as the explainable
   bridge from Stage 2 to Stage 3.
2. Build the 1P9J interpolation v0 as the runnable down-payment on this model:
   e3nn, T2-valued, within-axis, honest "direction not destination" graph.
3. Build the conditioned predictor incrementally by feature tier and mechanism
   family, with ablations. Start with a small, defensible source-family set
   rather than one grand all-calculator network.
4. Add 720-WT transfer only after the static ingest exists and passes its oracle
   parity gate.

The `.tex` "one network" picture is a unifying architecture, not the MSc build
plan. A single all-calculator net over every emitted feature would be too wide
and too easy to over-read on current data. The spec therefore recommends staged
feature tiers:

```text
v0: ring or D_ab-family e3nn interpolation on 1P9J held-out DFT
v1: per-mechanism dominance-clean e3nn fits
v2: conditioned predictor with classical + MOPAC + selected geometry conditioners
v3: AIMNet2 / Larsen / OF3 / boosts only after explicit ablation and lead vet
v4: 720-WT transfer
v5: frozen DFT-trained BMRB dark shot
```

Minimum viable thesis tier: v1 per-mechanism dominance-clean e3nn fits plus the
v2 DFT-ground conditioned predictor. v3-v5 are upside after ablation and lead
vet, not commitments required before the thesis result stands.

## Dependencies

The Stage-3 build is blocked on:

- **Chewer / pybind11 live C++ model binding.** Must expose vectorized source
  sets, target rows, emitted features, and metadata to Python/GPU without a
  Python protein model or per-pair Python loop.
- **720-WT static ingest.** Needed for transferability and the between-axis
  baseline. Cross-reference target:
  `SPEC_720_STATIC_INGEST_2026-06-04.md` (design present; build pending).
- **1P9J interpolation spec.** The v0/down-payment should align with this
  architecture. Cross-reference target:
  `SPEC_INTERP_1P9J_SCOPING_2026-06-04.md`; the built v0 result is recorded in
  `BUILD_INTERP_RESULT_2026-06-04.md`.
- **Feature naming decisions.** OF3 and boosts need concrete producer/C++ names.
- **BMRB mapping.** Needs atom-name/proton-equivalence mapping and an explicit
  ensemble caveat.

## Lead Open Items

1. Confirm the v1 source-family list: ring + charge + MOPAC-field + McConnell
   only, or include H-bond/Larsen/AIMNet2 from the first conditioned run.
2. Define OF3 and boosts concretely: source files, columns, units, irreps, and
   whether they are source features, target conditioners, or scalar heads.
3. Decide chewer delivery for training: persistent tensor cache produced by C++
   query vs live pybind11 tensor batches. Either way, no Python spatial model.
4. Decide radial sharing taxonomy: by mechanism, by source type, by atom
   stratum, or hierarchical source-type x target-type.
5. Decide BMRB comparison policy: per-atom mapping, proton averaging,
   uncertainty/replicate handling, and what graph is allowed to say.
6. Confirm whether any T1 diagnostic remains in the Stage-3 reporting. The
   architecture permits vector channels, but the thesis target remains total T2.

## Recommended Build Order

1. Finish and vet `SPEC_INTERP_1P9J_SCOPING_2026-06-04.md`; use it as the v0
   graph-producing down-payment.
2. Specify the chewer tensor API in C++ terms: target rows, source query, feature
   blocks, shapes, present masks, and ArraySpec metadata.
3. Implement the smallest e3nn conditioned source-sum model by extending the
   existing exemplar. Do not re-author the model class from scratch.
4. Run held-out DFT interpolation on 1P9J with tensor R2 as the primary score.
5. Add feature tiers and mechanism-family ablations.
6. After `SPEC_720_STATIC_INGEST_2026-06-04.md` exists and the ingest passes its
   oracle gate, run 720-WT transfer.
7. Only after DFT/720 behavior is understood, freeze and shoot at BMRB in the
   dark.
