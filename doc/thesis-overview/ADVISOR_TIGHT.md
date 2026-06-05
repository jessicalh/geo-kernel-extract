# Investigation of Classical Frame Calculations in an NMR Learning Model
### Tight advisor overview — ERECTOR SET (scaffold for the ~1 h build)

> Shape = the précis, less stuff. Stage 1 reports the genuinely interesting ridge + the
> distribution-in-space, stats described properly. Stages 2–3 trimmed to essentials.
> `[source]` tags mark where each beam's content is staged. Numbers below are vetted
> (the in-protein truth-vet PASSed them) unless tagged `[confirm]`. Voice is yours; this is girders.

## The question  *(2 sentences, trimmed)*
- Can the NMR shielding **tensor** of a nucleus be accounted for, from **protein geometry alone**, by classical physics kernels (ring current, EFG-from-charge, bond anisotropy, H-bond) calibrated against DFT?
- Aim is **physics explanation, not a predictor**: *which* geometric mechanisms carry the signal, and *how much*. T2 (the angular part) is kept end-to-end.

---

## Stage 1 — Static mutant shielding  *(the reported centerpiece)*

**Direct question:** is the kernel set rich enough to carry the signal across many proteins — and *where in the geometry* does the signal live?

**Dataset (one line):** 720 WT+ALA AlphaFold pairs, mean pLDDT ≥ 95 (AFDB-confirmed, 226 taxa, ~87% archaeal); target = WT−ALA static DFT Δσ per matched atom. *(The same 720 DFT serves twice: as the WT−ALA **delta** here, and later as **absolute** shielding to train Model 1.)* `[orcas-720, alphafold-accessions]`

### A. The ridge — what it captures
- per-element, per-atom-type ridge over the **126-kernel T2 layout**.
- **R² = 0.818** on the 110-protein fair set (69K atoms); **holds at 0.718** across all 720 / 446K atoms — the drop is cross-protein generalisation, not a physics failure. `[LEARN_STAGE1]`
- backbone full-model R² (total |T2|): **0.775 (O) → 0.902 (HA)**; total-iso 0.374 (O) → 0.950 (HA). `[in-protein-stats]`
- all-atom medians: total-iso **0.852**, total |T2| **0.812** (242 strata). `[in-protein-stats]`
- "nitrogen is hard" is a **pooling artifact**: backbone-N ≈ 0.39 vs sidechain-N ≈ 0.89 (2nd-easiest type). `[LEARN_STAGE1]`
- **effective angular dimensions differ by element: H ≈ 20, C ≈ 6, N ≈ 3, O ≈ 12** — the spread *is* the story. `[LEARN_STAGE1]`

### B. Distribution in space — where the signal sits
- local conformation (φ/ψ): all-local T0 CA/N/O = 0.649 / 0.460 / 0.429; O is **ψ-led** (0.304); all-local T2 CA/N/O = 0.547 / 0.411 / 0.318. `[in-protein-stats: dihedral]`
- local-geometry / "geometry-cloud" probes: signal organises by **local geometric configuration**, selected by geometry not atom-name. `[positional summary; topology-workup]`
- backbone φ/ψ angle-contributor maps: signed z-statistic patterns (direction, not magnitude). `[in-protein-stats: OPUS packet]`
- **contributors** (who carries it): ring-current scalar leads **T0** (biot_savart 0.734, haigh_mallion 0.721, ring_susceptibility 0.723); q/r³ + EFG leads **T2** (coulomb 0.763, mopac_coulomb 0.756); ring current also carries T2 (0.662). Rings explicit. `[in-protein-stats: contributors]`

### C. Stats — properly described
- ridge: β = (XᵀX + λI)⁻¹ XᵀY; R² = 1 − SSE/SST; per-protein z-normalisation (skip groups < 3). `[stats-tests]`
- the test battery: drop-one **mechanism independence**, protein-**blocked CV**, **pairwise independence** (top-20 + named-atom), **weak-but-independent-signal**, **T0/T2 & dia/para split**, **atom-role stratification**, **geometry-cloud probe**, **torsion-pair path**. `[stats-tests]`
- R² is a **diagnostic** that kernels carry signal — never the optimised metric (model ≠ physics).

---

## Stage 2 — The laws  *(précis shape, trimmed)*

- **one classical object** behind the calculators: σ_ab(i) ≈ Σ_s I_s · D_ab(r_is), with D_ab = ∂_a∂_b(1/ρ) = (3ρ̂ρ̂ − δ)/ρ³. *[fig: Model-1 architecture-is-the-equation]*
- within-axis on **1P9J**: charge q/r³ (0.278) and ring (0.275) recover; the **unified D_ab-sum recovers the through-space tensor (0.432)** — a *real* combine, carried by MOPAC-field + McConnell, charge ≈ 0. `[in-protein-stats: Stage 2]`
- honest framing: 1P9J is the **within** instrument; the between-axis is deferred to the **720-WT statics pilot**.
- grading: **statistical position, not literature-match**. Correlate, do not match.

---

## Stage 3 — The learning model  *(forward arc — the part they'll push)*

- **Ethos flip:** the question changes here — from *which mechanisms carry the signal* (explanation; R² a diagnostic) to *do they help predict the measured shift* (utility; **R² is the metric**). The explanation goal isn't dropped — Stage 3 tests whether the explained physics **earns its place**.
- **Model 1 — small shielding emulator.** e3nn equivariant; **trains on the DFT** (720 statics + 1P9J frames). Stages 1–2 set *which features and what structure* it carries (kernel shadows, the D_ab combine), not its training data. Predicts σ where there is no DFT. Prototype already runs: held-out 1P9J **T2 R² ≈ 0.48**. *[fig: Model-1; capability graph]*
- **Model 2 — big shift predictor (the wild west).** Eats Model 1 + the extraction firehose on the ~656 no-DFT fleet, trained vs **BMRB** shifts. **Firewall:** the DFT-trained Model 1 is frozen; the noisy experiment only ever trains Model 2 — it never corrupts the physics. *[fig: Model-2 wild-west]*
- **The payoff — ablation.** Toggle structure quality / atom pLDDT, Model 1 in/out, and the calculators in/out → read what each contributes to the real observable. This is where *"does any of this help"* gets answered.
- arc: **distill → transfer → test** (1P9J is the calibration case — it carries both DFT and BMRB).

---

## The arc, in one line
**signal → equation calibrations (each law + coefficient + confidence) → ensemble → equivariant transferability pilot.**

---

## Figures / panels on hand  *(for selection in the hour)*
- `fig_overall_flow` — **the whole arc, tying it together** (Mermaid) — likely the opening figure
- `fig_shielding_model` — Model 1, architecture-is-the-equation (TikZ)
- `fig_shift_model` — Model 2, the wild west (Mermaid)
- `interp_1p9j_advisor_graph` — capability graph (held-out T2 R² ≈ 0.48)
- dataset panels: `taxa_map`, `plddt_distribution`, size histograms

> **Open beams for the hour:** title; which Stage-1 numbers make the cut vs. sit in back-pocket; figure selection (which 2–3 actually go on the desk); how hard to lean on Stage 3 (the advisor will); the one-line framing of the archaeal-set + pLDDT≥95 selection.
