# Markwick et al. 2010 — deferred deep-dive notes

**Status:** Deferred 2026-04-27 ~03:44, with submission due 14:00 same day. Decision: write a minimal-risk Markwick annotation now; defer the deeper version until time permits. If the deeper version doesn't get written, the Markwick depth flows into the literature-review outline section instead.

**Outline integration (Jessica's framing):** *"describe Markwick's approach to short ensemble and consider applicability and loss in adopted approach"* — this captures both the precedent (short trajectories with sophisticated method scaffolding) and the gap-to-Jessica's-kit (longer plain MD without that scaffolding) that Markwick raises but doesn't directly answer.

---

## The conditioning insight (the load-bearing observation)

The advisor cites Markwick 2010 as evidence that short MD trajectories can support chemical-shift prediction work. Read carefully, the paper actually demonstrates that *short trajectories work when they are paired with substantial methodological scaffolding* — twenty acceleration-level parameter sweeps, free-energy reweighting via MM/PBSA, and ensemble construction at each acceleration level. Without that scaffolding the short-trajectory claim does not transfer.

This matters because Jessica's kit uses 25 ns of plain MD per protein on 720 proteins, *without* the AMD machinery, parameter sweeps, or per-protein method-tuning that Markwick relied on. The advisor's reading of Markwick as a "short MD works" precedent compresses out the scaffolding-conditioning that makes the original result valid.

The legitimate intellectual move (which Jessica's outline could carry) is to articulate that the trade-off is real but defensible: Markwick optimises depth on one protein for absolute experimental match; Jessica optimises coverage across many proteins for kernel-vs-DFT correlation. The epistemic targets differ. Markwick's scaffolding is needed for absolute experimental match; Jessica's plain-MD-at-scale is sufficient for kernel signal extraction (because DFT calibration, not experimental matching, provides the absolute scale).

---

## Direct quotes worth retaining for later use

**Authors' framing of why ensemble averaging matters:**
> "the experimentally measured chemical shift for a given nucleus therefore reports on a free-energy-weighted average over all conformational substates explored by the protein up to the millisecond time scale (the so-called chemical shift coalescence limit). In light of this, there is a fundamental limitation on how accurately one can interpret or predict chemical shift data using a 'single-copy' or static-structure representation of the system."

**Headline result with full breakdown:**
> "a cumulative RMSD of 5.92 ppm, which represents a 20% improvement relative to the standard MD simulations and a 28% improvement relative to the static X-ray crystal structure."

**Cross-validation against an independent NMR observable (the strongest validation in the paper):**
> "The optimal reproduction of the chemical shift data was obtained from trajectory-averaged molecular ensembles generated using the acceleration parameters [E_b(dih) − V(dih) = 600 kcal/mol, R(dih) = 120 kcal/mol]. Interestingly, this coincides with the acceleration level that best reproduced the experimental RDC data for this system, as both chemical shift data and RDCs report on an ensemble and time average over the millisecond range."

**The averaging-is-essential claim:**
> "The extensive averaging procedure is very important: no individual molecular ensemble obtained at the optimal acceleration level gave a chemical shift prediction as good as the trajectory-averaged result."

**Per-residue mapping to dynamic regions:**
> "The improvement in the predicted chemical shift data is directly correlated to those regions of the protein that exhibit substantially enhanced conformational space sampling."
> "Significant improvement in the predicted ¹⁵N chemical shift data clearly coincides with those regions of the protein that exhibit backbone dynamics on longer time scales."

**¹⁵N as the dominant channel:**
> "The most notable improvement in the chemical shift results was obtained for the ¹⁵N nuclei."

**Closing scope statement:**
> "for proteins with a heterogeneous distribution of dynamics over a hierarchy of time scales, chemical shift prediction from free-energy-weighted AMD ensembles provides substantial improvement."

---

## Headline numbers (verified verbatim from paper)

**Per-channel RMSD (ppm) against experimental shifts on IκBα:**

| Treatment | HN | N | Cα | Cβ | C′ | Cumulative |
|---|---|---|---|---|---|---|
| X-ray (PDB 1NFI) | 0.62 | 2.89 | 1.45 | 1.31 | 1.92 | 8.19 |
| 5 × 10 ns plain MD averaged | 0.57 | 2.60 | 1.31 | 1.18 | 1.78 | 7.44 |
| Optimal AMD ensemble | 0.44 | 1.84 | 1.01 | 1.01 | 1.62 | 5.92 |

Improvement: 28% vs X-ray, 20% vs plain MD. ¹⁵N alone improved 36% (2.89 → 1.84).

**Optimal AMD parameters:**
- E_b(dih) − V(dih) = 600 kcal/mol
- α(dih) = 120 kcal/mol

Same parameters were optimal for RDC reproduction in Markwick et al. 2009 (J. Biomol. NMR 45, 17 — ref 12b in this paper).

---

## Two draft annotations (proposed but not adopted)

### Path A — extended (~180 words, covers more findings + the conditioning critique)

> Evidence that ensemble-averaged chemical shift prediction substantially outperforms static-structure prediction in proteins with heterogeneous dynamics, with the gain coming from accelerated-sampling methodology rather than from the trajectory length itself.
>
> This study predicts backbone chemical shifts for the ankyrin-repeat protein IκBα using the empirical SHIFTX predictor on each frame of an accelerated molecular dynamics ensemble, averaging across free-energy-weighted conformers. Cumulative RMSD against experimental shifts fell from 8.19 ppm for the X-ray crystal structure to 5.92 ppm for the optimal accelerated trajectory — a 28% improvement — with the largest single-channel gain on backbone nitrogen (2.89 → 1.84 ppm) and per-residue improvement matching regions of independently measured millisecond exchange. The same acceleration parameters also gave optimal agreement with residual dipolar couplings in earlier work, providing independent cross-validation; ubiquitin, a less dynamically heterogeneous protein, showed correspondingly smaller improvement. The authors argue that ensemble averaging is essential — no single frame from the optimal acceleration level matched the trajectory-averaged result — though the demonstration depends on substantial methodological scaffolding including twenty acceleration-level parameter sweeps, free-energy reweighting, and per-level ensemble construction. Relevant to deciding the simulation-time budget when the methodological scaffolding is not available.

### Path B — focused (~150 words, foregrounds the conditioning critique)

> Evidence that ensemble-averaged chemical shift prediction substantially outperforms static-structure prediction, with the gain scaling with the system's dynamic heterogeneity.
>
> This study computes backbone chemical shifts for the ankyrin-repeat protein IκBα using the empirical SHIFTX predictor on each frame of an accelerated molecular dynamics trajectory, averaging across the resulting free-energy-weighted ensemble. Cumulative RMSD against experimental shifts fell from 8.19 ppm for the X-ray crystal structure to 5.92 ppm for the optimal accelerated trajectory — a 28% improvement — with the largest single-channel gain on backbone nitrogen and the spatial pattern of improvement matching independently measured regions of millisecond exchange. The "short trajectory" result depends on substantial methodological scaffolding — twenty acceleration-level parameter sweeps, free-energy reweighting, and per-level ensemble construction — making this a sophisticated-method-conditioned demonstration rather than a license for short plain MD alone. Relevant to choosing the simulation-time budget when that scaffolding is not available.

---

## What Jessica would need to internalize before writing the deeper version

1. **AMD specifics**: how the dual-boost dihedral biasing actually works (Hamelberg-McCammon 2004 J. Chem. Phys. 120, 11919), why E_b(dih) and α(dih) are the two adjustable parameters, what changes when you sweep them.
2. **MM/PBSA reweighting**: how implicit-solvent Poisson-Boltzmann calculations let you reweight AMD frames to recover unbiased free-energy weights.
3. **Ensemble construction**: how the 20% surviving (post-pruning) frames are clustered and seeded into classical MD for free-energy statistics.
4. **The per-protein parameter sweep**: why each protein needs its own optimal acceleration level (and why Jessica's 720-protein scope can't afford this).
5. **The relationship between AMD-effective sampling time and physical wall-clock simulation time**: how 10 million AMD steps maps to ms-effective sampling.

Without internalizing at least 1-3, the deeper annotation reads as paraphrase rather than understanding. Jessica's instinct (defer rather than fake) is the right one.

---

## What ended up in the live entry

A minimal-risk version (TBD — to be written when returning to entries-working.txt). The conditioning critique flows into the literature-review outline as the bullet *"describe Markwick's approach to short ensemble and consider applicability and loss in adopted approach"*.
