# 13-TR build plan (Stage 2 trajectory queue closeout) — 2026-05-20

**Status**: live planning doc. Locked decisions below are commitments; the strawmen are the agreed *starting* point for codex annealing.

**Scope**: complete the trajectory-scope calculator queue with exactly 13 `*TrajectoryResult` classes. Retire all other pending calculator entries per the explicit triage at the end of this doc. After these 13 land + the 4 Python analysis scripts ship, the calculator-side of Stage 2 is closed and we move to chapter writing.

**Sequencing principle**: minimise the number of working sessions; within each session, ship the canonical TR of each pattern first so later TRs are pure clones.

---

## The 13 TRs

| # | TR | Source ConformationResult | Lifecycle | Cadence | Emission shape | Clone template |
|---|---|---|---|---|---|---|
| 1 | `AIMNet2EmbeddingTimeSeriesTrajectoryResult` | `AIMNet2Result` | FO | every frame | (N, T, 256) float32 + Blosc2 | `AIMNet2ChargeTimeSeriesTrajectoryResult` |
| 2 | `AIMNet2PolarisabilityTimeSeriesTrajectoryResult` | `AIMNet2PolarisabilityResult` (Vec3 + scalar) | FO | every frame | (N, T, 3) double + (N, T) double | `HydrationGeometryTimeSeriesTrajectoryResult` |
| 3 | `AIMNet2PolarisabilityWelfordTrajectoryResult` | same as #2 | AV | every frame | 4 Welford channels per atom (3 Vec3 components + 1 scalar) | `HydrationGeometryWelfordTrajectoryResult` |
| 4 | `ApbsEfgTimeSeriesTrajectoryResult` | `ApbsFieldResult.apbs_efg_spherical` (T2 SphericalTensor) | FO | every frame | (N, T, 5) double | `BsShieldingTimeSeriesTrajectoryResult` |
| 5 | `MopacChargeWelfordTrajectoryResult` | `MopacResult` (charges, scalar per atom) | AV | sparse (every 20 ps) | 1 Welford channel per atom | `HBondCountWelfordTrajectoryResult` + sparse gate from `TripeptideBackboneResidualVecTimeSeriesTrajectoryResult` |
| 6 | `MopacBondOrderWelfordTrajectoryResult` | `MopacResult` (bond orders) | AV | sparse | 1 Welford channel per bond (uses bonds.npy from topology sidecar) | clone of #5 with bond axis |
| 7 | `MopacCoulombShieldingTimeSeriesTrajectoryResult` | `MopacCoulombResult` (T2 SphericalTensor) | FO | sparse | (N, T, 5) double | combines #4 (T2 TS) + #5 (sparse gate) |
| 8 | `MopacMcConnellShieldingTimeSeriesTrajectoryResult` | `MopacMcConnellResult` (T2 SphericalTensor) | FO | sparse | (N, T, 5) double | clone of #7 |
| 9 | `MopacVsFf14SbReconciliationTrajectoryResult` | `MopacCoulombResult` + `CoulombResult` (FF14SB) | FO | sparse (both attached) | (N, T) double — `|cos(MOPAC_T2, FF14SB_T2)|` | new cross-source pattern; only TR of this shape so it is its own canonical |
| 10 | `RingNeighbourhoodTrajectoryStats` | `RingGeometryResult` + per-shielding-kernel internal ring loops | FO | every frame | Pattern A — (N, R_per_atom_max, K) rich struct vector; see strawman below | new canonical (Pattern A) |
| 11 | `RmsdTrackingTrajectoryResult` | positions only (no ConformationResult) | AV | every frame | (T,) double | new canonical (reference-frame scalar) |
| 12 | `RmsdSpikeSelectionTrajectoryResult` | `RmsdTrackingTrajectoryResult` (TR-on-TR read) | AV | every frame | SelectionBag entries + (T,) uint8 spike mask | `ChiRotamerSelectionTrajectoryResult` |
| 13 | `DftPoseCoordinatorTrajectoryResult` | `RmsdSpike` + `ChiRotamerSelection` SelectionBags | FO at Finalize | every frame (passive) + once at Finalize | deduplicated SelectionBag under DftPoseCoordinator kind | new canonical (cross-TR SelectionBag aggregator) |

### Use-case split

- **1P9J deep-dive narrative**: TRs 5-9, 11-13. The MOPAC family validates against DFT at 20-ps cadence; RmsdTracking + Spike + DftPoseCoordinator feed the ORCA pose-selection pipeline.
- **600+ fleet narrative**: TRs 1-4, 10. AIMNet2 polarisability replaces MOPAC at fleet scale; APBS EFG TS pairs with the already-shipped APBS E-field TS; RingNeighbourhood gives cross-protein ring-current decomposition.

---

## Locked decisions

These are commitments. Changes require explicit user sign-off, not a codex finding alone.

1. **RMSD reference frame** = trajectory frame 0. Per Wingens 2003 the starting-conformer choice matters for 1P9J's hinge-bending; frame 0 IS the starting conformer.
2. **RMSD atom selection** = backbone heavy atoms (N, CA, C, O). Captures hinge motion; excludes solvent + cap H noise.
3. **RmsdSpike threshold** = 2.0 Å vs baseline (frame 0) + 0.5 Å vs rolling 100-frame mean; cooldown = 100 frames (~2 ns at 20 ps stride).
4. **DftPoseCoordinator read sources** = `RmsdSpike` + `ChiRotamerSelection` SelectionBags. **Dedup key** = `(residue_index, frame_index // ns_bucket_size)` with `ns_bucket_size = 50` frames (~1 ns).
5. **AIMNet2 EmbeddingTS storage** = float32 + Blosc2-bitshuffle compression (per `feedback_embedding_float32` + `project_compression_plan`). Default-on. Group attr `optional_large=true` so SDK consumers can skip.
6. **AIMNet2 PolarisabilityWelford channels** = 4 (3 Vec3 components + 1 scalar). Don't recompute scalar as derived from Vec3; emit both.
7. **Pattern A layout for RingNeighbourhood**: strawman below; codex annealing in Session 3.

---

## Pattern A strawman for RingNeighbourhood (open — Session 3 codex pass)

Per-atom-per-ring rich struct vectors. **Strawman** (NOT locked):

```
/trajectory/ring_neighbourhood/
  ring_membership_per_atom: (N, R_per_atom_max) int32
    — ring_index per ring slot; sentinel = -1 for unfilled slots within cutoff
  distance:                  (N, T, R_per_atom_max) double — Å
  normal_projection:         (N, T, R_per_atom_max) double — Å (signed)
  in_plane_angle:            (N, T, R_per_atom_max) double — radians
  bs_contribution:           (N, T, R_per_atom_max, 5) double — T2 SphericalTensor, ppm
  hm_contribution:           (N, T, R_per_atom_max, 5) double — T2 SphericalTensor, ppm
  mc_contribution:           (N, T, R_per_atom_max, 5) double — T2 SphericalTensor, ppm
  per-frame metadata (frame_indices, frame_times)
  attrs: ring_cutoff_angstrom, k_max_rings_per_atom, source_attached_policy
```

`R_per_atom_max` is the maximum number of rings within 15 Å cutoff observed across the trajectory + a safety margin. Sentinel `-1` ring_index for unfilled slots.

**Open design questions for codex**:
- Flat fixed-R layout (above) vs ragged variable-length per-atom lists. Storage cost on 600-protein fleet matters.
- Are bs/hm/mc contributions worth the 3× storage cost, or should we recompute from geometry at analysis time?
- Should ring_membership_per_atom be static (frame-invariant since topology is static) or per-frame (allows for ring perception evolution)? Static recommended.

---

## Session sequencing

### Session 1: Pre-flight + AIMNet2 trio (~3-4 hrs) — **LANDED 2026-05-20**

**Pre-flight** (30 min): write this doc, lock decisions 1-6. ✓ done by reading this.

**TRs landed**:
- TR 1: `AIMNet2EmbeddingTimeSeriesTrajectoryResult` (canonical AIMNet2 always-attached TS — 256-dim float32, optional_large=true, per-atom hyperslab)
- TR 2: `AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult` (canonical Vec3+scalar TS — renamed from "Polarisability" per science S2)
- TR 3: `AIMNet2ChargeResponseGradientWelfordTrajectoryResult` (canonical full canonical Welford row — mean/std/m2/min/max/min_frame/max_frame per channel)

**Codex annealing pass** complete (R0 + R1 rounds — 7+5 findings landed across two fix bundles).

**Commits** (master, 510681f..02cd4df, 9 commits, all green):
- `510681f` AIMNet2EmbeddingTimeSeriesTrajectoryResult + 13-TR session plan
- `3e953ec` AIMNet2PolarisabilityTS + adversarial-review fix bundle
- `fac1e46` AIMNet2 TR pair: codex R0 followups
- `4ac6ccd` AIMNet2 fleet trio S1 complete: NaN-fill gate, TR #3 Welford
- `f339586` science+math review followups (WelfordFinalize, physics disclaimers)
- `8297535` S7 + S9: charge-conservation pre-flight script + stale docstring fix
- `58594f5` Rename AIMNet2Polarisability → AIMNet2ChargeResponseGradient
- `094ddb7` OBJECT_MODEL: fix stale polarisability dataset names
- `02cd4df` AIMNet2 CRG codex R1: F1+F2+F3+F4+F5+F6+F7 fix bundle

**Verification at close**:
- 6/6 C++ tests PASS on the trio test binary
- Integration1P9J test (real AIMNet2 backward through Trajectory::Run): 846/846 atoms populated, max|grad_x| ≈ 0.4 e²/Å
- 135/135 Python SDK tests PASS

**S2 entry conditions met**: templates for TR4 (T2 EFG TS, canonical Vec3+scalar TS) ready; the Welford-with-full-canonical-row pattern documented in `AIMNet2ChargeResponseGradientWelfordTrajectoryResult`; F1 non-finite guard pattern + ForceAttachResultForTesting helper landed.

### Session 2: ApbsEfgTS + MOPAC family (~4-5 hrs)

Canonical-first ordering — each TR clones the previous:

- TR 4: `ApbsEfgTimeSeriesTrajectoryResult` — canonical T2 EFG TS
- TR 5: `MopacChargeWelfordTrajectoryResult` — canonical sparse Welford scalar
- TR 6: `MopacBondOrderWelfordTrajectoryResult` — clone with bond axis
- TR 7: `MopacCoulombShieldingTimeSeriesTrajectoryResult` — canonical sparse T2 TS
- TR 8: `MopacMcConnellShieldingTimeSeriesTrajectoryResult` — clone of TR 7
- TR 9: `MopacVsFf14SbReconciliationTrajectoryResult` — cross-tensor |cos|

**Codex annealing pass** at session close.

### Session 3: RingNeighbourhood (~3-4 hrs)

- Step 1 (30 min): read existing ring-related code (RingGeometryResult, BS/HM ring loops in BiotSavart + HaighMallion) to enumerate per-(atom, ring) data already computed inside the shielding kernels. The TR captures this decomposition, doesn't re-derive.
- Step 2 (60 min): **codex Pattern A design pass**. User pastes a codex-bound prompt. Layout locks here.
- Step 3 (90 min): implement to locked layout.
- Step 4 (30 min): literature-anchored test (1P9J Phe33 ring → Hα of Cα-stem residue at ~5 Å in the X-ray).

**Codex annealing pass** at session close.

### Session 4: Rmsd trio + DftPoseCoordinator + 1P9J production (~3 hrs)

- TR 11: `RmsdTrackingTrajectoryResult` — canonical reference-frame scalar
- TR 12: `RmsdSpikeSelectionTrajectoryResult` — clone of `ChiRotamerSelection` with RMSD-threshold push
- TR 13: `DftPoseCoordinatorTrajectoryResult` — canonical cross-TR SelectionBag aggregator
- Production validation: full 1P9J run with all 13 TRs; confirm 13 H5 groups land + tests pass + Python SDK reads all groups

**Codex annealing pass** at session close.

### Session 5: Codex final annealing + Python chapter analyses (~3-4 hrs)

- Multiple codex annealing rounds on the full 13-TR bundle. **Codex review at this stage is an annealing set, not single-shot** — multiple rounds expected, each prompt narrower than the last.
- 4 Python analysis scripts in `learn/`:
  - `block_avg_convergence.py` — Grossfield-Zuckerman 2009 block-averaged SEM. Gates every other Stage 2 claim.
  - `sigma_lipari_szabo.py` — per-atom σ-S² + per-residue groupby + cross-report against BMRB 5801 NH-S².
  - `sigma_essential_dynamics.py` — SVD on per-frame σ tensor. Tells us if signal is low-rank.
  - `markwick_overlay_1p9j.py` — σ-residual ↔ BMRB 5801 Rex spatial overlap. **The 1P9J chapter-section deliverable.**

---

## Codex annealing loop policy

- **Cadence**: at every session close, plus an extended round during Session 5.
- **Form**: I draft a codex-bound prompt; user pastes into their codex session; codex investigates with full repo+tooling access; user pastes findings back; I bundle fixes; we iterate.
- **Annealing means multiple rounds per session, not one-shot**. Each round's prompt is narrower than the last (broad audit → specific concern → specific concern).
- **Final pass (Session 5) expects 3-5 rounds** to drive the bundle to commit-readiness.

---

## Markwick deliverable

**Goal**: a chapter section about 1P9J showing spatial overlap (or lack thereof) between our σ-residual-per-residue and BMRB 5801 Rex-per-residue. This is **the** 1P9J chapter deliverable; codex called it "something to say about 1P9J" and the user agreed it's high-value.

**Inputs (all already extracted or trivially pre-stageable)**:
- σ residual per residue — derived from the 7+ shielding TS NPYs (BS, HM, Mc, PiQuad, RingSusc, Dispersion, HBond) via groupby on `residue_index_per_atom`
- BMRB 5801 Rex per residue — already in `references/incoming/` for the calibration set; pre-stage as a per-protein CSV/NPY (30 min one-off)

**Output**: `learn/markwick_overlay_1p9j/figure.pdf` + 1-2 paragraphs of chapter prose.

**Python cost**: ~80 LOC. Zero protein-model rebuilding in Python (the residue_index_per_atom broadcast handles binding).

---

## Retired calculators (with rationale)

After the 13 TRs land, the following pending entries in `spec/PLANNED_CALCULATORS_*` retire. Update those docs in Session 5 with retire markers + rationale.

| Candidate | Reason for retirement |
|---|---|
| `SigmaMSMResult` | Markov state model needs ms-timescale sampling. 15 ns relaxation MD has ~3-4 rotamer transitions maximum. |
| `SigmaMemoryKernelResult` | Mori-Zwanzig formulation needs long-memory regime. 15 ns insufficient. |
| `PseudocontactShiftResult` | No paramagnetic substrate in scope. 1P9J is diamagnetic; fleet is too. |
| `CrossCorrelatedRelaxationResult` | Backbone CSA-dipolar cross-correlation needs ns-ms; 15 ns timescale-marginal. Methyl version interesting but not in current scope. |
| `GreenKuboSpectralDensityResult` | Spectral view redundant with `SigmaEssentialDynamics` PCA at our timescale. Revisit only if ED reveals residual frequency-domain structure. |
| `SigmaTuckerDecompositionResult` | Higher-order ED. Defer behind PCA. Conditional trigger: PCA shows multi-axis structure that 2D SVD can't unmix. |
| `Lie-group GP regression` | Stage 3 modeling work, not a calculator. |
| `SigmaLipariSzaboResult` (as C++ TR) | Ships as Python analysis script in Session 5, not as a TR. |
| `SigmaEssentialDynamicsResult` (as C++ TR) | Same. Python analysis. |
| `BlockAveragedConvergenceResult` (as C++ TR) | Same. Python analysis. |
| Aromatic-H geometry sanity check | Diagnostic only. mdtraj one-liner at write-up. |
| `Per-secondary-structure CSA stratification` (as C++ TR) | Groupby on existing DSSP TS NPYs. Python at chapter draft. |

Retire markers go into the source doc as `**RETIRED 2026-05-20**: <reason>` blocks, preserving the original entry for future archeology but flagging it as out-of-scope.

---

## Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Pattern A turns into half-day design | Medium | Half session | Codex annealing pass scheduled (Session 3 step 2). Strawman above is the starting point. |
| AIMNet2 EmbeddingTS too large (~3.6 GB/protein uncompressed, fleet ~2 TB) | High (re: size); Low (re: solving) | Storage budget | Blosc2-bitshuffle (already planned). Verify compression ratio ≥3× in Session 1. If lower, sub-sample frames. |
| MopacBondOrderWelford bond-axis storage may need new sidecar | Low | One extra commit | Topology sidecar already emits `bonds.npy` — verify Welford-against-bond-axis is supported in existing pattern. |
| MopacVsFf14SbReconciliation source-pairing edge cases | Medium | Wrong output on edge frames | Gate strictly on BOTH `MopacCoulombResult` AND `CoulombResult` attached; emit NaN for asymmetric attachment. |
| DftPoseCoordinator dedup ns-bucket may be too coarse | Low | Fewer poses than desired | ns_bucket_size = 50 frames is the strawman; revisit during validation if pose count is too low / too high. |
| Codex finds object-model violation late in Session 5 | Medium | One extra round | Session-close annealing catches most; Session 5 multi-round annealing catches the rest. |
| Emergent calculator surfaces during a session despite this plan | LOW (pre-flight clean) | Half session | Stop, plan, then add. The pre-flight on 2026-05-20 confirmed all sources exist; this should not happen. |

---

## What "done" looks like

End of Session 5:
- 13 TRs landed on master with all tests passing + sibling-TR regression clean
- `python/nmr_extract/_trajectory.py` exposes 13 new SDK groups + `JCouplingTimeSeriesGroup`-style dataclasses for each
- `OBJECT_MODEL.md` catalog rows updated; pending markers (⏳) flipped to landed (✓)
- `references/ANNOTATED_BIBLIOGRAPHY.md` updated where new citations land (Pattern A literature, Grossfield-Zuckerman, Lipari-Szabo if not already there)
- 4 Python analysis scripts under `learn/` produce figures + summary stats
- Markwick chapter-section figure + prose draft ready for thesis
- Retired calculator entries in `spec/PLANNED_CALCULATORS_*` carry **RETIRED 2026-05-20** markers
- Codex annealing pass on the final bundle returns no HIGH findings

Stage 2 calculator queue is closed. Forward work is calibration + chapter writing.
