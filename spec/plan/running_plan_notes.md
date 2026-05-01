# Running plan notes

**Status:** Live working notes. Refine as work lands. Captures the picture as understood at the end of 2026-04-30 conversation. Pairs with `thesis-discussion-part1-2026-04-30.md` (MSc framing + calculator pass items) and `md-rerun-685-discussion-priors-2026-04-30.md` (MD/MDP priors).

**Latest update 2026-04-30 PM:** §"Adversarial review findings" added at the bottom — substantial corrections to the assistant's earlier reasoning chain on cadence/storage/replica shape. Read that section together with the rest. Several earlier conclusions in this document are now superseded; specifically the "1 × 15 ns at slow=50 ps" landing was a phantom-storage-ceiling-induced miscalibration and the existing 1 × 25 ns × 40 ps single-tier protocol stands.

---

## The four chunks

1. **Chunk 1 — nmr-extract calculator pass.** IN PROGRESS. Steady coherent commits. AMBER charge slice cleanup + IUPAC integration + calculator improvements (ring-normal stability, single-loop default audit, H-bond geometry, methyl-specific kernel, EF coefficient verification, per-element CSA validation tests, differential calculator, S² aggregator, block-averaged convergence). Calculator pass items detailed in `thesis-discussion-part1-2026-04-30.md` §5.

2. **Chunk 2 — OF3 + tripeptide setup.** STARTING. New extraction tree ready. Existing OF3-CS pipeline at `gotham:/home/jessicalh/of3-cs-project/` to be made systematic. Tripeptide assembler exists but missing the long-range corrections Larsen 2015 ProCS15 included (see §First Model below). IUPAC stratification plumbed in.

**Two-tier protein set structure (both active, do not collapse):**

```
2400 proteins        RefDB pool                              raw availability

1200-2000 proteins   curated RefDB + BMRB                    OpenFold3 structure
    "broad set"      (BMRB experimental shifts attached)     + embedding tier

685 proteins         subset of broad set                     + MD trajectory
    "deep set"       (calculator-extraction binding)         + per-frame kernels

10 proteins × 26 fr  sub-subset of deep set                  + DFT shielding
    "anchor set"     (existing 260-DFT batch)                  full T2 tensor

2 proteins × 1.5 µs  separate dense validation runs          + DFT at ~40 ps
    "deep validation"  (1Z9B + 1P9J, in flight)               + published S²/T1/T2
```

The broad set (1200-2000) is the GNN training scope. The deep set (685) is where MD ensemble features live. The anchor set is where T2 ground truth is. The deep validation set is where high-density per-frame DFT plus published relaxation observables anchor trajectory-feature claims.

3. **Chunk 3 — MD runs + ingest.** LATER, but MD options must be set NOW because MD takes time both to run and to extract. Forward bets section below.

4. **Chunk 4 — Final model feed + play.** LATER. Combines outputs from 1+2+3, trains incrementally as data arrives.

---

## Calculator inventory — complete to my knowledge as of 2026-04-30

User rule: **all calculators proposed are kept on principle.** This is the comprehensive inventory pulled from existing nmr-extract code, the gotham `training_matrix.md` channels, Track A/B/C agent reports, and earlier conversation. If a calculator/feature is mentioned anywhere upstream, it's here. List-of-next items are marked.

### Existing in nmr-extract (per memory + project knowledge)

- BiotSavart ring-current kernel (T2 tensor output)
- Haigh-Mallion ring-current kernel (T2 tensor output)
- McConnell-style kernel
- Buckingham-style kernel
- EFG (electric field gradient) kernel
- Bond anisotropy kernel
- ApbsField (electric field from APBS PB solver)
- AIMNet2 group (per memory `project_aimnet2_contract_20260426.md`: required in production output)
- MOPAC charges + bond orders (per memory `project_mopac_role.md`: QM charges + bond orders for every conformation)
- MutationDeltaResult / DftCompare (Stage 1 mutant calibration; per memory completeness debt: needs WT original + mutant compared + diamagnetic + paramagnetic, six tensors + three deltas)
- AtomTopology / IupacAtomName / AtomReference (IUPAC topology layer; per memory landed 2026-04-26)
- Tripeptide assembler (Chunk 2 — exists but missing Larsen-style long-range corrections)

### Training matrix channels (gotham `of3-cs-project/docs/training_matrix.md`)

Channels A-H + D1-D4 from the existing systematic-decomposition design:

- **A** OF3 si_trunk transformer (DONE — current baseline)
- **A'** OF3 si_trunk MLP (failed baseline; confirms si_trunk is relational not per-residue)
- **B** Tripeptide tensors at minimum
- **C** rSCAN DFT at minimum
- **D** MACE-DFT at minimum (validation: MACE reproduces DFT)
- **E** APBS at minimum
- **F** MD-derived features
- **G** CPMHD-derived features (protonation switching MD, 10 ns)
- **H** CSPRED features
- **D1** MACE-DFT × MD ensemble, Boltzmann-averaged
- **D2** MACE-DFT × CPMHD ensemble, averaged
- **D3** MACE-DFT × MD + CPMHD combined
- **D4** MACE-DFT at 2-3 APBS minima, population-weighted (no MD needed)

### Chunk 1 calculator pass items (Track C cache + literature → kernel additions / improvements)

1. **Ring-normal stability** — two-normal averaging instead of single 3-atom normal (Sahakyan-Vendruscolo 2013). Affects BiotSavart + Haigh-Mallion kernels. Audit existing implementation; small change.
2. **Single-loop vs two-loop ring current** — audit BiotSavart loop separation default against Agarwal 1977 [10]-paracyclophane benchmark (single-loop S=0 fits experiment with r=0.9928 vs two-loop S=1.28 Å r=0.8883).
3. **H-bond angle θ as primary feature** — angle correlates with shifts; distance only weakly (Yi-McDermott 2024). New per-frame H-bond geometry calculator: donor, acceptor, d(donor-acceptor), θ(donor-H-acceptor).
4. **Buckingham σ^(1)·E coefficient verification** — sanity-check our EF kernel against analytical values (Buckingham 1960 / Case 1995): −3.0 ppm/au nucleic acids, −3.4 ppm/au proteins.
5. **Per-element CSA validation tests** — unit test cases asserting tensor magnitude reaches literature values at known geometries (Boyd-Skrynnikov 2002: 16.6 ppm CSA for N-H⋯π RC; Sahakyan-Vendruscolo 2013: up to 80 ppm ¹⁷O via EF).
6. **EF dominance audit for ¹³C / ¹⁵N / ¹⁷O** — verify our calibration weighting matches Sahakyan-Vendruscolo 2013 finding that EF dominates RC for heavy nuclei (R=0.702 vs 0.257 for ¹³C).
7. **CH3Shift methyl decomposition kernel** — new methyl-specific kernel mirroring Sahakyan 2011's δ = δ_rc^rot + Δδ_dih + Δδ_ring + Δδ_ma + Δδ_EF + Δδ_dist decomposition. Side-chain coverage gap currently.
8. **Static a^MD vs averaged a^Xray calibration discipline** — methods-chapter requirement; document explicitly in `learn/CLAUDE.md` and code comments. Per-frame-fitted vs frame-averaged-fitted `a` coefficients are NOT equivalent (Li-Brüschweiler 2012).
9. **Block-averaged convergence calculator** — `BlockAveragedConvergenceResult` TR (already in `PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` Idea 1). Methods-chapter requirement.
10. **Differential kernel calculator** — Δkernel(frame) = kernel(frame) − kernel(reference). Reference state choice still open (see §forward bets).
11. **Lipari-Szabo S² calculator** — backbone N-H autocorrelation; methyl C-H autocorrelation. Already in `PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` as Idea 2 `SigmaLipariSzaboResult`.
12. **Backbone N-H unit vector calculator (in protein-aligned frame)** — feeds S² and Lesovoy/Orekhov 2025-style direct relaxation back-calculation. Trivial: positions only.

### Additional candidates from Track C cache (not in Part 1 §5 but raised by sources)

13. **Differential Biot-Savart per-volume-element shielding density Σ_αβ(r)** (Pelloni-Ligabue-Lazzeretti 2004) — re-formulates ring-current integral as per-volume contribution; could become a finer kernel decomposition.
14. **Ring-current methods comparison** — HM, JB, Pople-PD-new directly comparable per Moyna 1998 (HM r=0.854, JB r=0.846, PD-new r=0.846 on 11406 ¹H shifts). Implementation choice + benchmarking calculator.
15. **CH3Shift Δδ_dist phenomenological distance term** — the active-region/neutral-region polynomial in Sahakyan 2011's decomposition; specific form to copy.
16. **Larsen ProCS15 explicit corrections** — Δσ_BB^(i±1), Δσ_HαB, Δσ_sc as named correction terms; fits into the tripeptide assembler architecture (Chunk 2).
17. **Per-residue secondary-structure-stratified analysis** — Karp 2014 disordered-region calculator; per-residue DSSP class as feature.
18. **Per-residue SASA** — Gadanecz 2026 / cache references; cheap geometric.

### Additional candidates from Track B literature

19. **LEGOLAS as per-frame predictor channel** (Darrows 2025) — only included if integrated into the extractor pipeline, NOT as external comparison. Outputs 6 backbone scalar shifts per residue per frame.
20. **TSENN-style equivariant tensor architecture** (Hsu 2026) — model-side, informs network design for T2 prediction.
21. **XPaiNN-style equivariant tensor architecture** (Wang 2026) — same.
22. **Max-entropy reweighting against BMRB** (Yu 2026 PNAS) — model-side validation/reweighting strategy; uses pre-computed CS predictions per frame.
23. **Hydration shell calculator** — 2 closest waters per polar atom (2026 Chem Rev recommendation; 3.5 Å cutoff per fragment-QM/MM literature). Per-frame indices + positions.
24. **Per-residue dihedral KDE** — Gadanecz 2026 EDMD-style as FEATURE (not as refinement target). Per-residue φ/ψ kernel-density estimates over the trajectory.

### Per-frame trajectory features feeding multiple calculators

25. φ/ψ/ω, χ1-χ4 dihedrals per residue per frame
26. DSSP class per residue per frame (IDPForge 2026 lineage)
27. Per-residue SASA per frame
28. Hydrogen bond network table per frame (donor, acceptor, d, θ — feeds calculator #3)
29. Solvent positions (for hydration shell calculator + per-water EF computation)
30. Box vectors per frame
31. Energy components in EDR (per-step)
32. Velocities per frame (extra data switch — opens NMRforMD-style direct relaxation analyses)
33. RMSF per residue per frame
34. RMSD from reference per frame
35. Radius of gyration per frame
36. Contact map per frame (pairwise distances)
37. Per-residue cross-correlation (dynamical cross-correlation matrix)

### Trajectory-derived aggregators

38. J(ω) at standard Larmor frequencies (600, 700, 850 MHz) — Lesovoy/Orekhov 2025; downstream of N-H autocorrelation
39. Spectral density mapping for relaxation observables
40. Per-atom kernel time-series moments (mean, variance, skew, kurtosis)
41. Autocorrelation tau_e (internal motion timescale per atom)
42. Per-residue dihedral autocorrelation
43. Differential time-series summaries (mean, variance, percentiles of Δkernel)
44. Block-averaged convergence statistics per kernel per atom
45. Existing TR aggregators (BondLengthStats, BsWelford, etc. per memory `project_tr_tests_start.md`)

### Methodology / discipline items (not calculators per se but load-bearing)

46. Diverging-errors discipline: DFT correction schemes do NOT transfer to ML CS predictors (Solid State NMR 140, 2025). Calibration must be empirical against BMRB, not against DFT-corrected references.
47. Cavender et al. 2025 LiveCoMS benchmarking review — methodology context for benchmark dataset selection.
48. RDC + RCSA experimental observable framing (Liu et al. 2018 *Nature Protocols*) — argues T2 preservation matters.
49. ENSEMBLE assignment / matching rule for pairing predicted CS to experimental BMRB shifts — already implied by IUPAC stratification.

### List-of-next (kept but timeline-deferred)

- Crankshaft cadence-aliasing diagnostic — sub-ps burst on 2-protein subset (Yi-McDermott 2024 motivated). Ruled out as fancy+novel for the main re-run; appropriate as list-of-next.
- ShiftML3-style ensemble equivariant transformer (Cordova 2026 ChemRxiv).
- Full Lipari-Szabo / Wingens-style relaxation analysis (R1, R1ρ, hetNOE fitting).
- Replica-averaged metainference (Papaleo 2018) — if we ever pull experimental restraints into MD.
- AMP-BMS/MM 2026 anchored sub-trajectory for DFT-quality MD on a small subset.
- Müntener-style high-throughput de novo BMRB expansion as calibration data evolves.
- ¹⁷O kernel investigations (Sahakyan-Vendruscolo 2013 shows up to 80 ppm sensitivity to E-fields).

### Cache hygiene flagged

- `references-text/bratholm-2017-procs15-qm-chemical-shielding-refinement-text-{1-7}.txt` files contain Venetos 2023 silicate text, not Bratholm 2017. PDF on disk; needs re-ingestion via standard script. Not blocking.

---

## First model = Chunks 1 + 2, static

The first working model does NOT need MD. It falls together when Chunks 1 and 2 land.

```
Tripeptide assembly  →  σ_BB per atom              (Chunk 2 enhancement)
nmr-extract calcs    →  Δσ_ring, Δσ_HB, Δσ_EF,    (Chunk 1 product)
                        Δσ_methyl, ...
OF3 si_trunk         →  per-residue context        (Chunk 2 product, existing)
                                ↓
                        small network
                                ↓
                T2 primary, isotropic projection,
                    with fallback isotropic head
```

Larsen 2015 ProCS15 decomposition is the assembly architecture:

```
σ^i  =  σ_BB(tripeptide context)        ← tripeptide assembler
      + Δσ_BB^(i±1)                     ← next-nearest backbone
      + Δσ_ring                         ← BiotSavart / HaighMallion (Chunk 1)
      + Δσ_HB                           ← H-bond calculator (Chunk 1)
      + Δσ_HαB                          ← Hα bonding correction
      + Δσ_sc, Δσ_EF, ...               ← side chain, electric field (Chunk 1)
```

Comparison baseline: existing `of3-cs-project` transformer (OF3 si_trunk only):
HA r=0.752, H r=0.542, N r=0.757, CA r=0.943, CB r=0.969, C r=0.688

The first model is what Chunks 1 + 2 compose into. It is the early deliverable that lands well before the 685-protein MD re-run completes.

---

## Chunk 2 specific items

- Add ring current to tripeptide assembler (matches Larsen 2015 Δσ_ring contribution).
- Add the full Larsen ProCS15 correction set: Δσ_ring, Δσ_HB, Δσ_HαB, Δσ_sc, Δσ_BB^(i±1) — likely sourced FROM the Chunk 1 calculators.
- IUPAC stratification plumbed into training data preparation.
- Set 75 expansion: 738 → ~1200 (subset of 2400 RefDB curated against BMRB).
- Pair representations (zij_trunk) wired into the model alongside si_trunk.
- Per-pose embedding extraction — single static structure pose per protein from OF3.
- Training infrastructure capable of incremental retraining as Chunk 3 trajectories complete (peel-off-as-we-go).
- Existing 11 model variants in `of3-cs-project/models/` represent prior architectural exploration; all kept as comparison points.
- Per-atom-type heads (per Stage 1 atom-type stratification reveal: sidechain N R²=0.887 vs backbone N R²=0.387).
- BMRB cross-reference plumbing (curated set already curated against BMRB; extraction needs to land per-atom shift records cleanly).

---

## Chunk 3 forward bets — must be set NOW because MD takes time

Hard rule: the MDP / output / extraction options we set now must align with what the Chunk 1 calculators are producing. Mismatch is very costly. Switches we can re-run if we screw up; the *thinking* about what to capture and how, we cannot revisit.

### Locked

- **Cadence: 10 ps** for all output. Crankshaft aliasing acknowledged as known limitation.
- **Force field: ff14SB / TIP3P** (status quo, has DFT calibration anchor).
- **No fancy + novel** (compressed-x-grps mixed cadence, OPC water swap, energygrps decomposition all OFF).

### Green-lit (extra data switches, low risk)

- `nstvout = 5000` — velocities at 10 ps. Opens NMRforMD-style relaxation analyses on list-of-next at zero MD-time cost.
- `nstenergy` denser. Finer kinetic/potential time series.
- H5MD output as an extra alongside XTC (if mdrun flag enables it cleanly). To-do item; not a dependency.

### OPEN forward bets — settle before MD starts

These need landing in this conversation or shortly after, because MD options must align with calculator work happening NOW.

1. **Replica shape on the 685.** Options (more configurations exist; these are common shapes):
   - `1 × 15 ns`: full autocorrelation horizon, no replica variance
   - `2 × 7.5 ns`: minimal replica variance, ~1.5 ns per-replica horizon
   - `3 × 5 ns`: replica variance + first-order S² (Papaleo lineage), ~1 ns per-replica horizon
   - `5 × 3 ns`: more independent samples, ~600 ps per-replica horizon
   - `n × m ns` for any (n × m) ≈ 15 ns
   - `1 × longer` (e.g., 1 × 20 ns or 1 × 30 ns) if budget allows
   The choice depends on what the model wants from the trajectory: replica variance, autocorrelation tau range, both, or other statistics. Cost of switching after run: re-run.

2. **Reference state for the differential calculator.** Options (more exist; these are the named candidates):
   - First post-equilibration MD frame
   - Minimum-energy frame from trajectory
   - Trajectory-mean kernel value (reference-free in spirit)
   - Frame-mean position (compute kernel on average geometry — synthetic pose)
   - OF3-predicted structure ("exact-to-bottle" reference)
   - PDB-deposited X-ray structure (where available)
   - NMR-deposited structure from BMRB (where available; multi-model ensembles)
   - Highest-population cluster centroid from MD trajectory clustering
   - Multi-reference scheme: different references for different atom types or residue classes
   - User may also want NO single reference; differential could be against multiple references and the model consumes all of them
   Cost of switching after run: re-extract from trajectory (cheap if positions kept).

3. **Solvent in XTC output.** Options:
   - Keep full system in XTC (enables hydration shell calculator + per-water EF + solvent dipole analyses)
   - Filter to protein-only at write time (smaller files)
   - Both via separate XTC writes (keep both, filtered subset is cheap convenience)
   Storage isn't binding. Cost of switching: re-run if filtered, free if kept.

4. **Aggregation list per trajectory.** What the extractor produces from each completed trajectory. The full inventory of aggregators is in §Calculator inventory items 38-45. The forward bet is which subset is computed at extraction time vs computed lazily on stored per-frame data. Cost of switching: re-extract (cheap if features are cheap; expensive if per-frame data was discarded).

5. **DFT pose strategy across the 685.** Existing 260-DFT batch is on 10 proteins × 26 frames (FES-forced poses; flawed by bad chain extractions per user). User's framing 2026-04-30: DFT is Stage 1 anchor / extractor-validation, not GNN training. Options for Stage 1 polish (more exist; user named these):
   - Re-run existing 260 with proper relaxation on OF3-exact-to-bottle structures
   - Run 100 curated proteins at minimum (proteins selected to be ones rSCAN converges on)
   - Try 685 at minimum (rSCAN won't converge on full variety)
   - Hold at existing 260 + 2-deep-runs DFT (~750 each at 40 ps)
   - Mix of above
   Cost of switching: re-run DFT (expensive but bounded; each ~hours).

6. **Reference-pose for ensemble averaging effect comparison.** Static-pose calculation kernel(representative pose) vs MD-ensemble kernel-mean(trajectory). The "static pose" choice options overlap with #2 reference state options. May be the same as #2 or different per chapter narrative.

7. **Per-frame DFT additions for differential channel anchor.** Existing 260 + 2-deep-runs ~1500 from the dense validation set. Forward bet: do we add per-frame DFT on a slice of the 685 (e.g., 50 selected proteins × 10 frames each = 500 more) to grow the differential anchor? Or stick with what we have?

8. **What the model's primary output channels are.** Options:
   - T2 tensor only (isotropic as projection)
   - T2 tensor + isotropic separate head
   - T2 + isotropic + S² head
   - T2 + isotropic + S² + Δprediction head
   - All-in-one tensor + scalar + relaxation (multi-task)
   Determines downstream architecture. Recoverable but shapes Chunk 4.

### 2 deep runs at 1.5 µs (existing dense validation set)

- 1Z9B + 1P9J in flight per memory entry `project_calibration_runs_underway_20260429.md` (1Z9B 2 µs scan1, 1P9J 5 µs scan2; ETA ~2026-05-05).
- DFT at ~40 ps cadence per memory entry `project_dense_validation_2proteins_20260424.md` → ~375 frames per protein.
- Role: validation/test for trajectory-feature claims. Demonstrates extract pipeline at high temporal density.

---

## What we are NOT deciding now

- Network architecture details (Chunk 4).
- Specific tuning hyperparameters (Chunk 4).
- Whether the GNN has a unique novel hypothesis (retired — see `thesis-discussion-part1-2026-04-30.md` §2 framing).
- Calibration architecture detail beyond static a^MD discipline (already in priors).

---

## Cross-chunk parallelism (loose)

```
Chunk 1 (nmr-extract)   ─────────────●──────────►   feeds Chunk 3 ingest
Chunk 2 (OF3 + tripep)  ─────●──────►              feeds first model + Chunk 4
Chunk 3 (MD + ingest)              ───MD running────●ingestion─►   feeds Chunk 4
Chunk 4 (model)                                         ●─partial─►─full→
```

First model lands at Chunks 1+2 confluence (much earlier than Chunk 3 completion). Ensemble features land incrementally as Chunk 3 ingestion peels off completed trajectories.

---

## Document hygiene

Live notes. Refine as the picture changes or work lands. When decisions in §"OPEN forward bets" lock, they migrate into a plan-of-record document. Until then, this is the running picture.

Cross-references:
- `thesis-discussion-part1-2026-04-30.md` — MSc framing, calculator pass items, raw findings.
- `md-rerun-685-discussion-priors-2026-04-30.md` — MD/MDP priors and assumptions.
- `gotham:/home/jessicalh/of3-cs-project/STATUS.md` — existing OF3-CS pipeline state.
- `gotham:/home/jessicalh/of3-cs-project/docs/training_matrix.md` — existing systematic-decomposition design.

---

## Adversarial review findings (2026-04-30 PM)

An adversarial agent reviewed the chain of reasoning leading to the proposed MD shape and extraction cadence for the 685-protein run. The review is substantial; it corrected several misframings and surfaced one urgent finding that has nothing to do with cadence. Captured here as the live position.

### Things that were wrong in the assistant's prior reasoning

1. **"40 ps was a previous-LLM claim" — misframing.** 40 ps is the **existing extraction protocol**. It's set in `spec/ANALYSIS_TRAJECTORY_2026-04-14.md` line 85 and traced to `nstxout-compressed = 10000` in `prod.mdp` (every 20 ps in MD) with stride 2 in extraction. 200 already-completed extractions run at 40 ps. The previous LLM didn't invent 40 ps; the user set it. The assistant should have read the existing config before treating it as suspect.

2. **The 10 TB storage ceiling was phantom.** Actual measurement on the existing schema:
   - 685 proteins × current uncompressed = ~1 TB total.
   - With realistic compression: ~400 GB total.
   - With compression + drop `aimnet2_embedding`: ~280 GB total.
   - The "few GB per protein" target the user named is not just hit — it's overshot by an order of magnitude on the side of headroom.
   The earlier conversation operated against a constraint that was 25-50× looser than framed.

3. **Compression ratio was optimistically stated as ~2.5×.** Actual measured aggregate: 2.39× on the full schema, 2.88× excluding `aimnet2_embedding`. Dense field tensors carrying physics signal compress at only 1.08-1.14× (numbers like `total_B_field`, `disp_shielding`, `mc_shielding`, `hbond_shielding`). Sparse arrays (which are all-zero — see urgent finding below) appear to compress at 1000× but that's sparsity, not compression.

4. **The N_eff "≥ 5 at fit τ_max" threshold was hand-waved and not the standard criterion.** The actual literature criterion for Lipari-Szabo fitting is **T ≥ 5-10 × τ_max** (Halford 2007, Genheden 2014), not N_eff at τ_max. The assistant mixed Sokal-style integrated-tau N_eff (a global quantity) with Lipari-Szabo fit stability (a local-to-τ quantity). Different math.

5. **"Specific maths at 10 ps" was wrong.** No specific physics in the current calculator catalog requires 10 ps cadence. Methyl rotation would (CH₃ rotation is on the 1-5 ps scale per Sahakyan 2011), but the existing schema does NOT have per-methyl orientation as a stored field, so 10 ps cadence wouldn't help unless the calculator catalog were extended. With the current calculators, 40 ps spans backbone N-H τ_e (50 ps to several ns) with margin and is the defensible mid-point. The assistant invented a constraint that wasn't there.

6. **Two-cadence on disk has hidden gotchas the assistant didn't surface:** every consumer (ui/, h5-reader/, python SDK, learn/) would have to learn two time bases; joint analyses across cadences require interpolation/decimation that introduces artifacts; the current H5 schema is single-cadence and a multi-cadence schema change touches 60+ catalog entries in `python/nmr_extract/_catalog.py`. Better path: **single fast cadence on disk + on-the-fly fitted parameters** (S², τ_e) computed in-memory at extraction and written as scalar `TrajectoryResult` fields. That's already the planned design (`SigmaLipariSzaboResult`).

7. **"Drop `aimnet2_embedding` based on uncertain utility."** The decision is probably right but skipped the prerequisite step — a one-day Ridge ablation on the existing 200 to confirm the 256-dim embedding doesn't add R² when downstream charges/multipoles are present. Skipping the test is fine for storage; not fine as scientific posture for a thesis. Add the ablation as a Chunk 2 prerequisite. AIMNet2 contract reconciliation (memory entry `project_aimnet2_contract_20260426.md` requires AIMNet2Group in production output) needs a check: does that mean the embedding specifically, or does AIMNet2-derived charges/multipoles in `charges/` and elsewhere satisfy the contract?

### The all-zero arrays — corrected diagnosis

**Initial alarm**: ~12 H5 datasets all-zero across all 5 proteins surveyed (1B1V_4292, 1CBH_192, 1HA9_5028, 1G26_4656, 1DV0_4757).

**After investigating the source**: most of these are **intentional skips by design**, not bugs. Only ~5 fields are genuine schema-without-producer orphans.

The corrected breakdown:

**INTENTIONAL EMPTY (skip_coulomb = true in `PerFrameExtractionSet`)**

CoulombResult is gated on `HasResult<ChargeAssignmentResult>() && !opts.skip_coulomb`. The `PerFrameExtractionSet` config (RunConfiguration.cpp lines 87-100, the production canonical for the 685-protein fleet) explicitly sets `skip_coulomb = true` because **vanilla Coulomb is O(N²) per frame and prohibitive on solvated boxes ~150K atoms × 625 frames × 685 proteins**. APBS supersedes — solves Poisson-Boltzmann for solvent-aware E-field per frame, much cheaper at large N.

So the following are reserved-but-empty BY DESIGN in the per-frame trajectory extraction:

```
efg/E_total, E_backbone, E_sidechain, E_aromatic, E_solvent
efg/E_magnitude, E_bond_proj, E_backbone_frac
efg/coulomb_total, coulomb_backbone, coulomb_aromatic, coulomb_shielding
```

These ARE populated when CoulombResult runs (Use Case A — static-pose / sparse-frame extraction). The 685-fleet uses Use Case B (`PerFrameExtractionSet`) which intentionally skips them.

**TRUE ORPHANS — schema slots with no producer in src/**

```
ring_current/G_iso_exp_sum                  declared in ConformationAtom.h
ring_current/G_T2_exp_sum                   read by AnalysisWriter +
ring_current/G_iso_var_8A                     fileformat/analysis_file
ring_current/mean_ring_distance             but no calculator code in src/
ring_current/nearest_ring_atom_distance       sets them to anything non-default
```

These are genuine orphan slots. Either implement the producer (Chunk 1 calculator pass decision: are the exponential-sum aggregates worth computing?) or remove from schema. Decision needed.

**POPULATED AND WORKING (the actual data the 685-fleet emits per frame)**

```
ring_current/bs_*, hm_*, rs_*  (per-type T0/T2 tensors + total_B_field)
ring_current/n_rings_3A/5A/8A/12A
efg/apbs_efield, apbs_efg              (per-frame E-field — the source)
efg/aimnet2_*                           (alternative E-field via AIMNet2 charges)
charges/aimnet2_charge, eeq_charge, eeq_cn
bond_aniso/* (McConnell)
dispersion/*, quadrupole/*, water/*, hbond/*, sasa/*, dihedrals/*, dssp/*
```

### Run state correction (2026-04-30 PM)

Earlier framing assumed 485 of the 685-fleet were "still running." User correction: **the 485 are cancelled, not in flight.** They were stopped because of bad chain extractions in the structure preparation. The next 685 run will use OF3-exact-to-bottle structures (matching the planned DFT re-run reference state). So:

- 200 completed extractions exist (Apr 14-15, 2026; Use Case B `PerFrameExtractionSet` output).
- 485 were cancelled before completion. No partial data to preserve.
- The 685-fleet starts fresh whenever the re-run is dispatched, with OF3-exact-to-bottle structures.

Implications:
- The earlier "don't change schema mid-flight for the 485" advice is moot — there is no mid-flight, the runs are stopped.
- Schema / cadence / compression / dtype / orphan-field decisions can be made cleanly without inhomogeneous-fleet concerns.
- The 200 already-completed extractions may themselves be flawed (same prep issue as the 485 likely). Whether to keep them, drop them, or re-run them is a separate decision.
- The next dispatch of 685 runs after Chunk 1+2 land + structures regenerated.

### Per use case (the actual extraction process structure)

```
USE CASE A — Static-pose / sparse-frame extraction
  RunConfiguration::FullFatFrameExtraction (or single-pose runs)
    skip_mopac    = false
    skip_apbs     = false
    skip_coulomb  = false       ← CoulombResult runs
    Output:       efg/E_*, efg/coulomb_*, mopac_* all populated
    Used for:     720-protein × 446K-atom Stage 1 mutant calibration

USE CASE B — Per-frame trajectory extraction (685 fleet target; the existing 200 had bad chain extractions, the 485 were cancelled before completion)
  RunConfiguration::PerFrameExtractionSet
    skip_mopac    = true        ← FullFat-only
    skip_apbs     = false       ← APBS every frame
    skip_coulomb  = true        ← APBS supersedes at N > 1000 atoms
    skip_dssp     = false
    stride        = 2           ← 25 ns × 1250 MD frames → 625 extracted at 40 ps
    Mandatory:    AIMNet2 every frame
    Output:       efg/apbs_*, efg/aimnet2_*, ring_current/bs_*/hm_*, etc.
                  (NOT efg/coulomb_*, NOT efg/E_*, NOT mopac_*)
    Per-frame E-field: efg/apbs_efield (solvent-aware via Poisson-Boltzmann)
    Per-frame EFG:     efg/apbs_efg

USE CASE C — DFT-anchor / scan-for-DFT-points
  RunConfiguration::ScanForDftPointSet
    skip_mopac    = true
    skip_apbs     = true
    skip_coulomb  = true
    skip_dssp     = false       ← dihedral bins need φ/ψ/χ
    Cheap subset for finding good DFT-input poses
    Attached:    BsWelfordTrajectoryResult
```

### Lesson from the corrected diagnosis

The earlier "calculator silently failing" framing was wrong. The adversarial agent saw all-zero arrays and assumed bug; the actual answer was `skip_coulomb=true` by design in PerFrameExtractionSet. **Adversarial review didn't read RunConfiguration.cpp** to find the run-shape choice. Lesson: when seeing schema-empty fields, check the run-config gates before declaring bugs. The check belongs in the same pass as reading the writer code.

The actual investigation finding is much narrower than first reported: only the 5 `ring_current/G_*_exp_sum` / `mean_ring_distance` / `nearest_ring_atom_distance` / `G_iso_var_8A` fields are true orphans needing a Chunk 1 decision. Everything else is working as designed.

### Other things missing from the conversation

8. **The 26 per-ns NPY-snapshot directories alongside each H5 are duplicating data at lower precision.** The data flow between H5 (per-frame at 40 ps) and the NPY directories needs to be named explicitly. If the H5 is authoritative, the NPY directories may be redundant; if the NPY directories are inputs to DFT scheduling and the H5 is downstream, the bifurcation needs documentation.

9. **`BlockAveragedConvergenceResult` (Idea 1 in `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`) is the actual diagnostic** that turns "626 frames" into "N_eff effective independent samples per atom per kernel." Until this calculator runs across the fleet, no defensible claim about "trajectory averaging gives independent samples" can be made. The assistant has been doing the cadence math without pointing at the user's own existing plan to compute exactly this.

10. **Schema versioning has not been addressed.** Cadence change / field drop / compression filter / dtype change all break consumers (ui/, h5-reader/, python SDK, learn/). A version-bump strategy is required before any schema change ships. This is a project-management gotcha that has burned the user before per memory entries.

11. **The 2-protein dense validation set (1Z9B + 1P9J at µs scales) IS the relaxation-observable bed.** The 685 fleet is **coverage**, not the deep validation per memory entry `project_dense_validation_2proteins_20260424.md`. Some of the "preserve every kernel time series feature for 685 proteins" anxiety relaxes if framing is correct: 685 is for breadth / model training; 2-deep is where relaxation observables live.

### Independent reconstruction of the 40 ps cadence justification

Constructed from the literature (Track C cache + the agent's own reading), arguing both for and against:

```
FOR 40 ps:
- Robustelli-Stafford-Palmer 2012 used 4.5 ps (overkill for ML-feature extraction)
- Li-Brüschweiler 2012 PPM: 100 ps cadence (defensible upper edge)
- Karp 2014 (IDPs): 5 ps full / 5 ns thermodynamic
- Yi-McDermott 2024: 100 fs short-window / 50 ns thermodynamic
- For backbone N-H τ_e ≈ 100 ps need cadence ≤ 30 ps; for τ_e = 200 ps need ≤ 70 ps
  → 40 ps spans both with margin
- For methyl rotation (τ ~ 5-10 ps), 40 ps captures rotamer occupancy averaged over
  many rotation periods (NOT direct rotation autocorrelation, which would need ≤ 2 ps)
- Storage at 40 ps: 0.30 GB/protein after compression and embedding drop, 200 GB total
  for 685 — well within 10 TB ceiling

AGAINST 40 ps:
- For backbone N-H S² fitting with τ_e ≈ 1 ns (slow folded limit), 40 ps gives 25
  points across τ_e — overkill
- For DFT-anchor frame extraction at 1 ns cadence, only every 25th frame matters
- For multimode MSM analysis (planned exploratory calculator), 40 ps is finer than
  needed (MSM lag times ≥ 1 ns typically)

VERDICT: 40 ps is a defensible mid-point. Not the only choice; 50-100 ps would also
be defensible. The "40 ps" inheritance from existing protocol is not suspect — it's
reasonable.
```

### Direct ACF measurement on existing data (1B1V_4292)

```
QUANTITY                                    τ_int (frames)   τ_int (ps)   N_eff (over 25 ns)
bs_iso (median high-var atom)               2-5              80-200       60-150
bs_T2_per_type (high-variance directions)   5-13             200-530      25-65
bs_T2_per_type (worst-case slow component)  13               530          23
bs_iso ACF: drops 1.0 (lag-0) → 0.15 (lag-1)
  → isotropic ring-current decorrelates within one frame at 40 ps
  → "ring-current temporal pattern is sacred" applies to T2 directional
     components (real τ_int 200-530 ps), NOT to the isotropic
     component (already past decorrelation timescale at 40 ps)
```

### Updated position on MD shape and extraction cadence

```
REPLICA SHAPE              1 × 25 ns (existing protocol; was the 200's shape)
                           Or whatever new shape is chosen for the next 685
                           re-run on OF3-exact-to-bottle structures.
                           No mid-flight constraint anymore — the 485 are
                           cancelled. Schema/cadence/replica decisions can
                           be made cleanly for the next dispatch.

ON-DISK CADENCE            40 ps single-tier
                           Existing protocol. Defensible mid-point for backbone
                           CS work. No specific physics in current calculator
                           catalog requires finer.

FAST CADENCE PROCESSING    In-memory at extraction time only. Lipari-Szabo fits
                           run on positions read at MD-write cadence (10-20 ps);
                           write only S² + τ_e per residue per atom + variance
                           summary. SigmaLipariSzaboResult-style.

SCHEMA CHANGES             Apply gzip-6 + shuffle compression (skip on
                           already-compression-resistant arrays).
                           Drop aimnet2_embedding/aim AFTER one-day Ridge
                           ablation on existing 200 confirms no R² loss.
                           Convert positions to float32 (precision floor for
                           atomic positions in MD is ~10⁻³ nm, well within
                           float32 precision).
                           Keep float64 on physics tensor outputs (ring_current,
                           efg, bond_aniso, etc.).

CHUNK 1 DECISION ITEM      The orphan ring_current accessory fields (no
                           producer in src/): G_iso_exp_sum, G_T2_exp_sum,
                           G_iso_var_8A, mean_ring_distance,
                           nearest_ring_atom_distance. Either implement
                           the producer (compute the exponential-sum
                           aggregates) or remove from schema. No deadline
                           pressure (485 cancelled), but should be settled
                           before the next 685 dispatch so the schema is
                           internally consistent. NOT a calculator bug;
                           schema-without-producer.

TOTAL STORAGE PROJECTION   685 × ~0.30 GB compressed = ~200 GB (with cleanups).
                           ~2-3% of the phantom 10 TB ceiling.

THESIS METHODS CAVEAT      At 25 ns we do not claim convergence in
                           Li-Brüschweiler's ensemble-averaged sense (their
                           floor is 100-200 ns). We claim per-frame kernel
                           features as MD-sampled snapshots. Per-frame
                           regression interpretation is defensible. Per-protein
                           average vs BMRB average comparison would be a
                           reviewer issue we'd need to name in the methods
                           chapter and discussion.
```

### What needs to land before scaling

In rough order:

1. **Investigate all-zero arrays.** Identify the calculator(s) emitting zeros. Decide: fix or remove from schema. Highest priority — pre-dates everything else.
2. **Schema versioning strategy.** Plan the version-bump for compression / dtype / field-drop changes; coordinate consumers (ui/, h5-reader/, python SDK, learn/) before any change ships.
3. **One-day Ridge ablation** on existing 200 with-and-without aimnet2_embedding to confirm safe to drop.
4. **AIMNet2 contract check** — what specifically is required by `project_aimnet2_contract_20260426.md` memory.
5. **Reconcile the H5/NPY data flow** — what's authoritative, what's input to DFT, what's redundant.
6. **`BlockAveragedConvergenceResult` calculator** — Idea 1 in PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md; the diagnostic that defends per-frame independence claims.
7. **Re-extraction strategy** — if compression and embedding-drop are applied, the 200 already-complete extractions need re-extraction from existing XTC. Bounded work; the XTC trajectories are authoritative.
