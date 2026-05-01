# MD re-run, 685 proteins — discussion priors

**Date:** 2026-04-30
**Status:** Discussion starter, not a plan-of-record. The MD/MDP/replica/cadence choices fall out of the thesis claim, not the other way around. The next step is a conversation that turns what we know into a thesis claim.

This document captures:

1. The MD-related context as it actually stands at the start of the conversation.
2. The assumptions Claude made during the walk-through that should be examined and discussed, not adopted.
3. The literature and local-cache findings that are live and bear on what the re-run should ask for — *not* "filed for later."
4. The questions we have not actually answered yet.

Three agent passes ran on 2026-04-30 and feed this document:

- **Track A** — comprehensive MDP audit against GROMACS 2026.0 reference. Full report at `/tmp/claude-1000/-shared-2026Thesis/1cf7a2bf-49d5-4964-96de-4162ec51dee0/tool-results/toolu_01H1y7mWY11H8aH8zftLMpZk.json` (~50 KB; saved because too large for inline display).
- **Track B** — Feb-Apr 2026 NMR-MD methodology literature scan, web-based.
- **Track C** — local references-text cache scan with Track A and Track B findings as context. Returned 2026-04-30 mid-walk-through.

---

## Section 1 — What is actually decided

Items the user has stated explicitly during the conversation, in their own words. These are not under discussion in this document; they are the constraints the discussion sits inside.

- **1.5 µs total per protein is OUT.** ("1.5us is right out.")
- **10-15 replicas of 10 ns each is doable with MPS.** ("With MPS 10 to 15 10ns is doable.")
  - "Doable" is what was said. Whether 10-15 × 10 ns is the right shape is part of what we have not decided.
- **No novel-to-us techniques on load-bearing one-shot runs.** ("The trouble with a student project on limited time with specialist techniques like moving things at different weights is there is no room for a full screwup with new tech.") Saved as memory entry `feedback_no_novel_techniques_on_oneshot.md`.
- **The 685-protein re-run is one-shot.** ("We have 1 re-run of all 685. Beyond that the re-run budget is infinite.")
- **Storage is not a binding constraint.** ("Well a big drive is purchasable.")
- **Wallclock is real.** Existing 15 ns runs on 1Z9B + 1P9J each took ~2 hours.
- **Re-run beyond this one is infinite cost.**
- **Calculator cadence determines model input volume.** ("The nature of our calculators limits cadence. We create a time series of 1.5 gb for 15ns or something. That feeds the model.")

---

## Section 2 — Assumptions Claude made during the walk-through that should be examined

These were presented as locked or recommended; they are not actually decided. Each is paired with the reasoning that *should have* preceded a decision but did not.

### 2.1 Frame cadence at 10 ps for everything

**What was said:** "Lock cadence at 10 ps for everything, no compressed-x-grps trickery."

**What was not actually established:**

- What does the calibration model do with the cadence? Per-frame ridge fitting? Frame-averaged kernel-DFT correlation? Something else?
- Is there a calculator for which 10 ps is too dense (oversampled, autocorrelated frames inflate effective sample size)?
- Is there a calculator for which 10 ps is too sparse (e.g., crankshaft motion at 0.23-0.89 ps state switching, see §4.5)?
- The cache (Li-Brüschweiler 2010, Robustelli 2012, Karp 2014) uses 1-10 ps for shift work. 10 ps is the upper edge — defensible but the upper edge.

**What we still need to know:**

- The model architecture, in enough detail to say "this cadence is right for this consumption pattern."
- Whether picosecond-aliased mechanisms (crankshaft) are tolerable noise or a known disagreement source we have to argue for in the methods chapter.

### 2.2 ff14SB / TIP3P as the force field

**What was said:** Treated as status quo throughout. "Keep ff14SB / TIP3P (the project's existing choice); the literature does not punish this for CS work."

**What was not actually established:**

- We have not asked whether the calibration corpus has any cross-FF requirement (e.g., does any DFT comparison require a different water model?).
- The re-run is downstream of fleet calibration; ff14SB is the per-pose force field used for the MD that produced the existing 260-DFT calibration set on 10 proteins. Switching now decouples the MD from that calibration anchor — but is that the right call given GROMACS 2026 has native ff19SB / OPC ports?
- TIP3P vs OPC for explicit-solvent shielding: Sahakyan-Vendruscolo 2013 says solvent E-field effects dominate for ¹³C and ¹⁵N; the question of *which* water model captures those E-fields most physically is unsettled in the local cache.

**What we still need to know:**

- Whether the existing DFT calibration on 10 proteins constrains the FF choice for the re-run.
- Whether the thesis methods chapter needs to defend the FF choice with literature, in which case "well-trodden" is enough; or with calibration evidence, in which case a single-protein FF comparison run might be needed.

### 2.3 10-15 × 10 ns replica strategy

**What was said:** "10-15 replicas × 10 ns. Locked." This was Claude treating user's "doable" as "decided."

**What was not actually established:**

- The cache's older lineage (Robustelli 2012, Yi-McDermott 2024, Li-Brüschweiler 2010, de Gortari 2010) prefers **single long trajectories** for free MD. Replicas in the cache appear when restraints are applied (Robustelli 2010 CS-MD with 21-temperature REMD; Papaleo 2018 with 4 replicas under CS+NOE restraints).
- The 10-20-replicas Lai-Brooks 2024 prescription is **S²-convergence-focused**. If the thesis does not claim S², the prescription does not bind.
- Papaleo 2018 found **4 replicas** optimal for ensemble refinement, not 10-20.
- 100-200 ns total simulation time can be built as: 1 × 100 ns; 4 × 25 ns; 10 × 10 ns; 15 × 10 ns. Each gives a different statistical surface.
  - 1 × 100 ns: best for chemical-shift ensemble averaging; full backbone autocorrelation tau range; no inter-replica variance.
  - 4 × 25 ns: Papaleo-style; modest replica variance; each replica has full autocorrelation range up to ~5 ns.
  - 10-15 × 10 ns: many independent samples; per-replica autocorrelation horizon is ~2 ns; replica variance is statistically meaningful; S² becomes possible.
- The choice depends on what the thesis is *defending*.

**What we still need to know:**

- What the thesis claims and which trajectory shape supports the strongest defense of that claim.
- Whether replica-variance bars on per-residue per-kernel quantities are part of what we want to report.
- Whether autocorrelation horizon (which limits what kernel time-series statistics can be computed) matters more than replica count for our purposes.

### 2.4 Velocities-on at 10 ps, forces-off

**What was said:** "Recommendation: `nstvout = 5000` (10 ps, matching positions). Forces I'd lean off."

**What was not actually established:**

- Velocities are needed for: kinetic-energy decomposition per group; some autocorrelation analyses; Lesovoy-style direct relaxation back-calculation. **Are any of these in the analysis plan?**
- "No paper argues for forces" was Claude filling silence, not a reason. **Forces are useful for force-field-consistency diagnostics** — comparing ff14SB-evaluated forces against expected analytical forces given charges and Lennard-Jones parameters. We may want this as a sanity check on the topology generation. We may not.
- The cache has nothing strong arguing for forces, but the cache is also not the field-current state.

**What we still need to know:**

- Whether the analysis plan uses velocities in any concrete calculator.
- Whether forces are part of any quality-control diagnostic we want for the re-run.

### 2.5 Other MDP knobs not surfaced

The Track A audit covers ~80 MDP options. Several were not even mentioned in the walk-through and warrant discussion:

- **`nstcheckpoint`** — current default is 15 minutes wallclock. For 10 ns runs that complete in ~2 hours under MPS, default is fine. But: **denser checkpointing enables finer-grained restart points if a run needs to be partially recovered.** Not a binding question; minor.
- **`nstenergy` / `nstcalcenergy`** — energy output cadence in the .edr. Currently coupled to position cadence. If we want denser energy reporting (for finer-grained kinetic / potential energy series), this is decouplable. Discussed nowhere yet.
- **`compressed-x-precision`** — XTC compression precision. Default 1000 (3 decimal Å). Trade-off: precision loss vs file size. Storage isn't binding, but we never asked whether the calibration is sensitive to 0.001 Å precision loss.
- **`refcoord-scaling`** — under pressure coupling, what reference coordinates rescale. Affects post-equilibration drift. Defaults vary by GROMACS version; worth checking.
- **`comm-mode`** — center-of-mass motion removal. Currently `linear` (default). Some MD-NMR work argues `angular` for finer rotational averaging analysis. Not investigated.
- **PME / cutoff parameters** — `rcoulomb = 1.0`, `rvdw = 1.0`, `pme-order = 4` are GROMACS defaults. Whether any of these should be tightened for protein NMR work is unknown.
- **Long-range LJ corrections (`DispCorr`)** — `EnerPres` recovers density better; affects pressure coupling stability. Not investigated.
- **Constraint algorithm (`LINCS` vs `SHAKE`)** — currently LINCS at default order. Robustelli 2010 used SHAKE on **all bond lengths AND angles** because CS-restrained MD demanded it. We are not CS-restrained; LINCS at default is presumably fine, but we never confirmed.

---

## Section 3 — What the three agent passes found, raw

### 3.1 Track A — MDP audit (full report on disk)

The full per-option matrix lives in the JSON file referenced above. Pulled here are findings that are load-bearing for the discussion, not summarized as recommendations:

- **`energygrps` is NOT supported on GPUs in mdrun.** Documented at `mdp-options.rst:474–477` ("not supported on GPUs"). Setting `energygrps = Protein SOL Ion` forces the short-range non-bonded calculation onto the CPU pathway with measurable wallclock penalty across all 685 × 15 (= 10,275) replica simulations. The penalty is enough that the agent surfaced it as the *single biggest finding*. 
- **System sizes verified:** 1P9J_5801 = ~73 K atoms (54 residues + 23933 SOL + 28 NA + 27 CL); 1Z9B_6577 = ~200 K atoms (135 residues + 66168 SOL + 284 K + 280 CL).
- **Run parameters verified from `dispatch_record.json`:** `duration_ns = 15.0`, `nsteps = 7,500,000`, `dt = 0.002 ps`, `nstxout-compressed = 5000` ⇒ frame cadence = 10 ps ⇒ 1500 frames per 15 ns trajectory.
- Recommended additions in the audit (these are agent recommendations, not decisions): velocities at modest cadence; H5MD output format; denser checkpointing; per-residue secondary structure; H-bond geometry; hydration shell descriptors.

The audit produced a **comprehensive option matrix** (every MDP option, current value, default, alternatives, mdrun runtime cost, scientific value added). It is on disk; we have not actually read through it together. Reading it together is part of the discussion.

### 3.2 Track B — Feb-Apr 2026 literature scan

Items the field-survey agent found in the 2026 publication window:

| Reference | What it argues |
|---|---|
| Piroozi et al. 2026 *Cell Rep Phys Sci* | Review framing equivariant GNNs and DFT+ML hybrids as state-of-the-art for ML-CS prediction |
| Akram et al. 2026 *Chem Rev* | Review of computational NMR; endorses MD-CS pipelines for flexible/disordered proteins |
| Gadanecz et al. 2026 *J Chem Inf Model* | EDMD: MD-derived φ/ψ distributions as Boltzmann-inversion potential for NMR-restraint refinement; 5 proteins, 13-19 kDa |
| Yu et al. 2026 *PNAS* | Max-entropy reweighting of CG-MD + ML CS predictors (UCBShift) against BMRB shifts for IDP ensembles |
| Cordova et al. 2026 ChemRxiv | ShiftML3 + UMA as drop-in DFT replacement for NMR crystallography (solid-state) |
| Müntener et al. 2026 bioRxiv | High-throughput automated NMR; 384 de novo proteins → 239 high-quality 2D NMR spectra |
| IDPForge 2026 bioRxiv | Transformer + diffusion model for IDP ensembles, validated against BMRB |
| AMP-BMS/MM 2026 ChemRxiv | DFT-quality protein MD up to ~100 ns at thousands of atoms |
| LEGOLAS 2025 *JCTC* | GPU-fast backbone CS predictor for MD trajectories; backbone-only |
| Lai & Brooks III 2024 *J Phys Chem B* | 10-20 replicas of tens of ns prescribed for S² convergence; ff14SB > CHARMM36m for NH order parameters |
| Lesovoy & Orekhov 2025 IJMS | Back-calculates R1, NOE, η_xy from N-H autocorrelations at 500 ps cadence |
| *Solid State NMR* 140, 2025 ("Diverging errors") | DFT and ML CS errors are distinct; correction schemes do NOT transfer between them |
| GROMACS 2026.0 release | Native AMBER ff14SB/ff19SB ports; OPC/OPC3 water models; NN/MM electrostatic embedding; H5MD trajectory output |

### 3.3 Track C — local references cache scan

The local cache has the **Vendruscolo-Cambridge through-line** that the field-survey agent could not see in the Feb-Apr 2026 window:

```
Cavalli 2007 CHESHIRE                CS → structure prediction (MC)
    ↓
Robustelli 2009                      same, refined MC
    ↓
Robustelli-Kohlhoff-Cavalli-         CS → MD as restraints (CS-MD with CamShift,
  Vendruscolo 2010 Structure           soft-square + tanh tail penalty,
                                       21-temperature replica exchange,
                                       SHAKE on all bonds AND angles,
                                       AMBER03 + Generalized Born, 2 fs,
                                       12 proteins, 11/12 folded to 0.8-2.2 Å)
    ↓
Robustelli 2012 JACS                 INTERPRETATION pivot — MD-averaged CS
                                       improves on X-ray for ALL backbone
                                       atoms in RNH (E. coli + T. thermophilus)
                                       100 ns, 99SB-ILDN, 4.5 ps cadence
    ↓
Papaleo 2018 PeerJ                   replica-averaged metainference on NCBD
                                       (Bayesian ensemble refinement,
                                       CS + NOE jointly, CHARMM22*,
                                       4 replicas optimal,
                                       28 × 50 ns separate runs for S²_relaxation)
    ↓
Gadanecz 2026 JCIM                   PHILOSOPHICAL INVERSION — instead of
                                       CS pulling MD, MD-derived φ/ψ
                                       distributions push NMR-restraint
                                       refinement. KDE-fitted Boltzmann-
                                       inversion potentials.
```

The cache also surfaces:

- **Li-Brüschweiler 2010 *J Phys Chem Lett*:** Plotted shift-RMSD vs trajectory length on calbindin D9k + ubiquitin. **Convergence sets in at 100-200 ns.** Argues "MD ~100 ns or longer should be most useful for this task."
- **Li-Brüschweiler 2012 *J Biomol NMR*:** PPM predictor. The load-bearing equation:
  ```
  δ_exp = Σ a_j^(k) · ⟨f_j^(k)(r_n)⟩ + δ_0^(k)
  ```
  Argues that **the fitted `a` parameters obtained from average X-ray structures are NOT equivalent to those fitted from MD ensembles**, even when ⟨r_n⟩_MD = r_Xray, because dynamic averaging gets absorbed into a_j^Xray.
- **de Gortari et al. 2010 *JACS*:** First-principles MD + GIPAW DFT averaging on the MLF tripeptide. **Fast bond-length oscillations (~20 fs period, 0.2 Å amplitude) produce instantaneous shift fluctuations of 20 ppm but average within 200 fs** — directly justifies our ps-spaced snapshots without zero-point corrections.
- **Yi-McDermott 2024 *JPCL*:** AF-QM/MM on DHFR-trimethoprim, 1000 ns, **20 snapshots at 50 ns spacing for QM/MM, plus a 100 ps segment at 100 fs cadence for torsion analysis**. Crankshaft motion in Ile60-Ile61: ψ_i and φ_{i+1} anticorrelated (r = -0.89), **state-switching residence times 0.23-0.89 ps, 25 ppm 15N excursion**. **Hydrogen bond ANGLE θ correlated strongly with ψ; H-bond DISTANCE only weakly.**
- **Sahakyan-Vendruscolo 2013 *JPCB*:** Direct analysis of Pople, Waugh-Fessenden-Johnson-Bovey, Haigh-Mallion ring-current contributions and electric-field effects on RNA base shifts. **Methodologically sharp finding: ring normals computed from three ring atoms are unstable to out-of-plane atomic fluctuations in MD. Their fix: average two independently computed normals.** Also: for ¹³C, EF dominates over RC (R = 0.702 vs 0.257); for ¹⁷O, EF sensitivity reaches 80 ppm.
- **Markwick et al. 2010 *JACS*:** AMD on IκBα. 15N RMSD 2.89 → 1.84 ppm with averaging (36% reduction). **15N is the largest-magnitude reporter of dynamics effects** (corroborated by Robustelli 2012, Yi-McDermott 2024).
- **Sahakyan et al. 2011 *JBNMR*:** CH3Shift method for methyl side chains. δ = δ_rc^rot + Δδ_dih + Δδ_ring + Δδ_ma + Δδ_EF + Δδ_dist. **Direct ancestor of our kernel-decomposition approach** — same physics families.
- **Buckingham 1960 *Can J Chem*:** Theoretical foundation of σ^(1)·E linear electric-field shielding. **Δσ ≈ −2×10⁻¹² E_z − 10⁻¹⁸ E².** At E ≈ 7 Å from a unit charge, the linear term is ~0.2 ppm — about 20× the quadratic term. Sets the analytical coefficient (-3.0 ppm/au for nucleic acids per Case 1995, -3.4 ppm/au for proteins).
- **Boyd-Skrynnikov 2002 *JACS*:** Closed-form σ_xz formula extending Johnson-Bovey σ_zz to the **full rank-2 ring-current tensor**. At G42 of fibronectin type-2 module, ring-current-dominant N-H⋯π hydrogen bond gives σ_rc CSA contribution of **16.6 ppm**.
- **Agarwal et al. 1977 *Can J Chem*:** [10]-paracyclophane benchmark. **Single-loop S=0 fits experiment best (r = 0.9928) vs two-loop (S = 1.28 Å gives r = 0.8883).** Direct flag for our `BiotSavartResult` loop-separation default.
- **Cache hygiene flag:** `references-text/bratholm-2017-procs15-qm-chemical-shielding-refinement-text-{1-7}.txt` and the corresponding summary file actually contain Venetos 2023 silicate text, not Bratholm 2017. PDF is on disk; needs re-ingestion.

---

## Section 4 — Live findings that bear on the re-run discussion (not "filed for later")

These are findings that Claude initially called "extraction-time decisions to defer" but which are load-bearing for shaping the re-run ask. They are open items for the discussion.

### 4.1 H-bond angle θ correlates with ψ and shifts; distance does not (Yi-McDermott 2024)

Currently, our extraction probably captures H-bond *distance* d(donor-acceptor). The Yi-McDermott finding says **angle is the stronger correlator** for backbone amide chemical shift dependence on backbone dihedral.

**Implications for the re-run discussion:**

- Per-frame extraction needs to compute θ(donor-H-acceptor). Positions are sufficient at our cadence; the MDP gives us what we need.
- **But** if H-bond geometry is not currently in any kernel, this is a *new* feature to add to the calibration model. That is a model-architecture decision, not an extraction nicety.
- We should look at what the existing 8 kernels actually consume and where H-bond geometry would slot in.

**Question:** Does any existing kernel use H-bond geometry as input? If not, is the right move to add a new kernel or to extend an existing one?

### 4.2 Ring-normal stability audit (Sahakyan-Vendruscolo 2013)

Three-atom ring-normal definitions are unstable to MD out-of-plane fluctuations. Their fix: compute two independent ring normals (e.g., atoms 1-2-3 and atoms 4-5-6 of a six-membered ring) and average.

**Implications for the re-run discussion:**

- Our `BiotSavartResult` and `HaighMallionResult` may have an MD-fluctuation bug. If the bug is real, we ship broken kernel time series into the re-run.
- The audit is **pre-re-run** — we want to know features are correct before we lock the MDP and run 685 × 10-15 trajectories that depend on the kernels.
- Audit shape: pick a 2-protein subset; compute ring-normal-1 and ring-normal-2 for each aromatic ring across all frames; histogram Δθ between them. If distributions are tight (Δθ < 1°), single-normal is fine. If distributions are wide (Δθ > 5°), two-normal averaging is required.

**Question:** Should this audit happen before the re-run, or alongside? Is it a 1-day task or a multi-day task?

### 4.3 Static a^MD vs averaged a^Xray calibration discipline (Li-Brüschweiler 2012)

The kernel coefficients fitted from per-frame DFT (a^MD) and from frame-averaged DFT (a^Xray) are not the same object. Mixing them gives wrong predictions.

**Implications for the re-run discussion:**

- This is the **first-order calibration architecture question.** What does our existing ridge regression in `learn/` do? Per-frame fitting? Frame-averaged fitting? Both?
- The 260-DFT calibration set we already have (10 proteins × 26 frames at 1 ns intervals) supports either path: per-frame (260 fits) or frame-averaged (10 fits).
- **The MD-side question that depends on this:** If we fit per-frame, the re-run needs to deliver per-frame kernel time series; cadence and replica count shape the per-frame statistics. If we fit frame-averaged, the re-run needs to deliver well-converged frame averages; total simulation time is the dominant variable.

**Question:** What is the calibration architecture? This is pre-MDP-lock for the re-run.

### 4.4 Convergence horizon: 100-200 ns total (Li-Brüschweiler 2010)

Cache empirically demonstrates shift-RMSD convergence sets in at 100-200 ns total simulation per protein.

**Implications for the re-run discussion:**

- This is the cleanest argument for total simulation time per protein, and it is **independent of replica strategy.** 100-200 ns can be 1×100 ns, 4×25 ns, 10×10 ns, or 15×10 ns.
- The replica-vs-single choice is downstream of: do we want autocorrelation horizon (single long), variance bars (replicas), or both?

**Question:** What does the analysis plan need from the trajectory shape — autocorrelation tau range, replica variance, per-trajectory ensemble averages, or a combination?

### 4.5 Crankshaft motion aliasing at 10 ps (Yi-McDermott 2024)

Picosecond-timescale crankshaft motion (ψ_i / φ_{i+1} anticorrelated, residence times 0.23-0.89 ps) produces 25 ppm 15N excursions. At 10 ps cadence we alias this motion.

**Implications for the re-run discussion:**

- 15N is the channel where the cache (Markwick, Robustelli, Yi-McDermott) consistently shows the largest dynamics signal.
- A 25 ppm aliased mechanism is **not noise.** If our 15N kernel-DFT correlation falls short by something like that magnitude on backbone amides, this is the candidate explanation.
- "Accept as finding" in the methods chapter is one option. "Run dense bursts" was ruled out (novel-to-us). A third option: live with aliasing as the baseline limitation but add a per-protein dense-burst diagnostic on a 2-protein subset to characterize what the aliasing costs us.

**Question:** Is the aliasing limitation tolerable for the thesis claim, or does it argue for a small dense-burst diagnostic alongside the main re-run?

### 4.6 Block-averaged convergence reporting (Markwick, de Gortari, Robustelli)

The cache's older lineage uses block averaging to defend convergence claims. Without it, our averaging claims do not survive review.

**Implications for the re-run discussion:**

- This is an **MD-design-affecting requirement**, not an extraction nicety.
- Single-trajectory block averaging: divide one trajectory into N blocks, compute kernel averages per block, report block-to-block variance. Standard MD analysis.
- Multi-replica block averaging: each replica is its own block; inter-replica variance is the convergence diagnostic.
- The two give different statistical surfaces.

**Question:** Which block-averaging story does the methods chapter want to tell? This shapes whether we want one trajectory long enough to block-divide, or many trajectories short enough that each is a block.

### 4.7 Per-element CSA expectations (Boyd-Skrynnikov 2002, Sahakyan-Vendruscolo 2013)

- Boyd-Skrynnikov: 16.6 ppm CSA for an N-H⋯π hydrogen bond from RC alone (rank-2 tensor σ_xz).
- Sahakyan-Vendruscolo: up to 80 ppm for ¹⁷O (mostly EF).

**Implications for the re-run discussion:**

- Validation targets for the T2 outputs from our kernels.
- If our σ_xz tensor never reaches these magnitudes at appropriate geometries, we are underestimating; if it routinely exceeds, we have a bug.
- This is post-extraction validation, not an MD-design choice. But it is **part of the methods chapter narrative**: "we validated kernel T2 outputs against literature CSA expectations at known geometries."

**Question:** Is this validation already done? If not, when does it happen relative to the re-run?

### 4.8 EDMD lineage as thesis context

The user found Gadanecz 2026 EDMD fascinating. The local cache supplies the lineage (Section 3.3 above).

**Implications for the re-run discussion:**

- The thesis methods chapter has a natural narrative arc through this lineage: CS-driven MD → MD-driven CS interpretation → ensemble refinement → MD-derived structural priors.
- Where does our work sit in this? We are not doing CS-driven MD (we run free MD). We are not doing structural refinement (we have crystal/NMR structures). We are doing **kernel-decomposed shielding prediction calibrated against DFT, with MD ensembles providing the conformation distribution that gets averaged over.**
- **Closest cache analog: Sahakyan 2011 CH3Shift** — direct ancestor of our kernel-decomposition approach; same physics families, polynomial functions of distances.

**Question:** Does the thesis position itself as kernel-decomposition (Sahakyan lineage, with MD as input) or as MD-ensemble-prediction (Vendruscolo lineage, with kernels as components)? The framing affects how the methods chapter argues for the MD setup.

---

## Section 5 — Questions we have not actually answered

These are the open questions whose answers determine the MD/MDP/replica/cadence choices, ordered roughly from most-foundational to most-immediate:

1. **What does the thesis actually claim?** Static shielding prediction via geometric kernels? Trajectory-averaged shielding prediction with explicit dynamics signal? Relaxation observables? All? The phrasing of the claim determines what observables the calibration model needs.

2. **What does the calibration model consume?** Per-frame kernel→shielding (static a^MD path) or frame-averaged kernel→shielding (a^Xray-like path) or both? See §4.3.

3. **Are relaxation observables (S², T1, T2, NOE, η_xy) part of the thesis story?** If yes, replica count matters per Lai-Brooks; if no, the cache argues for single long.

4. **What does the analysis plan need from the trajectory shape?** Autocorrelation horizon (single long), replica variance (many short), block-averaged convergence (either, but different stories), or some combination?

5. **Do the existing 8 kernels need any new inputs that the re-run must provide?** Specifically: H-bond angle θ; ring-normal stability inputs; per-polar-atom solvent positions for E-field calculation.

6. **Is the ring-normal stability audit pre-re-run or alongside?** If we ship broken `BiotSavartResult` / `HaighMallionResult` into the re-run, we waste the run.

7. **What is the calibration architecture today, and what should it be?** Reading `learn/CLAUDE.md` and the existing ridge regression code is part of answering this.

8. **Force-field choice — is ff14SB/TIP3P locked because of the existing 260-DFT calibration anchor, or is it a clean choice for the re-run?** GROMACS 2026's native ff19SB/OPC port is the alternative; we have not investigated whether it makes a defensible difference.

9. **Crankshaft aliasing at 10 ps — accept-as-finding or alongside-diagnostic?** §4.5.

10. **Velocities and forces — what does the analysis plan use them for?** §2.4.

11. **Replica directory layout, manifest schema, gen-seed schema** — downstream of the replica decision, but worth surfacing so we don't surprise ourselves later.

---

## Section 6 — What it would mean to land each open question

Brief notes on what answering each question above would require — not to answer them here, but so the discussion knows the cost.

- Q1 (thesis claim): user-led articulation; ~30 min discussion; needs `spec/MATHS_GOALS.md`, `spec/CONSTITUTION.md`, `learn/CLAUDE.md` re-read.
- Q2 (calibration consumption): same as Q1; concrete answer in `learn/` ridge regression code.
- Q3 (relaxation observables): user-led; cache supports either path with caveats per §2.3.
- Q4 (trajectory shape): falls out of Q1-Q3.
- Q5 (kernel inputs): code audit; Claude reads existing kernels + `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`; ~1 hour.
- Q6 (ring-normal audit): scope investigation; 1-2 hour task to characterize before deciding when to run the audit itself.
- Q7 (calibration architecture): code reading + discussion; ~1-2 hours.
- Q8 (FF choice): half-day investigation including a single-protein FF comparison run if needed.
- Q9 (crankshaft): user-led; if alongside-diagnostic, ~10 ps cadence on a 2-protein subset for 1 ns is sufficient.
- Q10 (velocities/forces): user-led; concrete answer in any pending or planned calculator.
- Q11 (layout/manifest): downstream; ~30 min once Q3-Q4 land.

---

## Section 7 — Re-run sequencing implication, if useful for the discussion

Not a plan; just a reminder of what depends on what:

1. Q1-Q3 (thesis claim, calibration consumption, relaxation scope) → drive Q4 (trajectory shape) → drive replica count and length → drive total wallclock budget.
2. Q5-Q6 (kernel inputs, ring-normal audit) → must land before MDP lock, because if a kernel is wrong, we waste the re-run.
3. Q7 (calibration architecture) → drives whether the re-run produces frame-by-frame data or aggregates; determines the H5 schema for the re-run.
4. Q8 (FF) → relatively independent; can be locked late.
5. Q9-Q10 (crankshaft, vel/forces) → MDP-level; lock late.
6. Q11 (layout) → lock last.

---

## Document hygiene

This document is `discussion priors`, not `plan-of-record`. When decisions land, they should be captured in a separate plan-of-record document with this one cross-referenced as the conversation that produced them. Until then, this is **not** a source of truth for what to do. It is a source of truth for **what we know and what we have not yet decided**.

If the document accumulates beyond what the discussion can productively use, sections should be moved to issue tickets or a separate scratch tree, not deleted.
