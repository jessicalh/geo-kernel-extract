# Larsen et al. 2015 — ProCS15 deep dive

**Citation.** Larsen, A.S., Bratholm, L.A., Christensen, A.S., Channir, M. & Jensen, J.H. (2015) "ProCS15: a DFT-based chemical shift predictor for backbone and Cβ atoms in proteins." *PeerJ* 3, e1344. DOI: 10.7717/peerj.1344. **Open Access (CC-BY 4.0).** 19 pages. Submitted 25 Aug 2015, accepted 1 Oct 2015, published 20 Oct 2015. Authors at U. Copenhagen (Pharmacy + Chemistry) and U. Wisconsin–Madison.

**Disclosures.** "The authors declare there are no competing interests." Funded by the Lundbeck Foundation and the Danish e-Infrastructure Cooperation; the funders had no role in study design or analysis.

**Code & data availability (verbatim).** "ProCS15 is freely available at github.com/jensengroup/procs15 and all structures and DFT calculations, including the full NMR shielding tensors, are available at erda.dk/public/archives/YXJjaGl2ZS1TYk40VXo=/published-archive.html."

---

## Abstract verbatim

> "We present ProCS15: a program that computes the isotropic chemical shielding values of backbone and Cβ atoms given a protein structure in less than a second. ProCS15 is based on around 2.35 million OPBE/6-31G(d,p)//PM6 calculations on tripeptides and small structural models of hydrogen-bonding. The ProCS15-predicted chemical shielding values are compared to experimentally measured chemical shifts for Ubiquitin and the third IgG-binding domain of Protein G through linear regression and yield RMSD values of up to 2.2, 0.7, and 4.8 ppm for carbon, hydrogen, and nitrogen atoms. These RMSD values are very similar to corresponding RMSD values computed using OPBE/6-31G(d,p) for the entire structure for each proteins. These maximum RMSD values can be reduced by using NMR-derived structural ensembles of Ubiquitin. For example, for the largest ensemble the largest RMSD values are 1.7, 0.5, and 3.5 ppm for carbon, hydrogen, and nitrogen. The corresponding RMSD values predicted by several empirical chemical shift predictors range between 0.7–1.1, 0.2–0.4, and 1.8–2.8 ppm for carbon, hydrogen, and nitrogen atoms, respectively."

---

## System + method

ProCS15 is a lookup-and-additivity scheme on top of pre-computed DFT shielding values. The shift on atom *i* is

> δ_i = b − a σ_i (Eq. 1)

with a and b empirically fitted to experimental shifts by linear regression per atom type. The shielding σ_i is built additively from terms each computed once and stored:

> σ_i = σ_BB^i + Δσ_BB^{i−1} + Δσ_BB^{i+1} + Δσ_HB^i + Δσ_HαB^i + Δσ_RC^i + Δσ_w^i (Eq. 2)

In plain words: the shielding for an atom in residue i is the shielding it would have in an isolated **Ac-AXA-NMe** (alanine-X-alanine, capped) tripeptide for that residue's backbone and side-chain dihedrals, plus a small change for the side chain of the residue *before* it, plus a change for the side chain *after* it, plus changes from hydrogen bonding to the amide H, the carbonyl O, the Hα, and (for amide protons only) explicit water. Ring current is added for protons. The summation **assumes these effects are additive** — that's the price of the speed.

The pre-computed pieces:

- **Backbone scans.** For each amino acid X, an Ac-AXA-NMe tripeptide is built by FragBuilder, optimised at PM6 with backbone and side-chain dihedrals frozen, and a GIAO NMR calculation is run at OPBE/6-31G(d,p) with CPCM continuum solvation (ε = 78). Backbone φ/ψ are scanned on a 20° grid. For amino acids with many side-chain dihedrals (χ_1, χ_2, …) BASILISK is used to sample 1000 side-chain conformations per φ/ψ pair. Flanking alanines are held at φ = −140°, ψ = 120° (β-sheet-like). "In total the ProCS15 backbone terms are based on ∼2.35 million DFT calculations." Failed optimisations are filled by interpolation. The grid loaded into RAM is ∼25 GB.
- **Hydrogen-bond scans.** Δσ_HB and Δσ_HαB come from scanning a hydrogen-bonded N-methylacetamide (and similar Hα model systems) over r, θ, ρ. Bond length 1.5–3.0 Å in 0.125 Å steps; θ from 90–180° in 10° steps; ρ over the full circle. Computed at OPBE/6-31G(d,p)//PM6.
- **Ring current.** Simple point-dipole model Δσ_RC = i·B·(1 − 3cos²θ)/r³, cutoff 8 Å, ring intensities from Christensen, Sauer & Jensen 2011.
- **Water.** Δσ_w = 2.07 ppm added to amide H not hydrogen-bonded inside the protein, calibrated against an N-methylacetamide-water DFT calculation.

The benchmark proteins are **ubiquitin** (1UBQ X-ray + five NMR-derived ensembles 1D3Z, 2K39, 1XQQ, 2LJ5, 2KOX, with 10–640 structures) and **GB3** (the third IgG-binding domain of Protein G, 2OED). Geometries are pre-optimised either with PM6-D3H+ in PCM or with the CHARMM22/CMAP force field with GB/SA solvation. Experimental shifts are from BMRB IDs 17769 (ubiquitin) and 18531 (GB3), both at pH 6.5. Sequence-corrected random coil values are subtracted before correlation to suppress trivial sequence-driven r-inflation.

---

## Headline numbers (verified verbatim)

**ProCS15 vs full OPBE/6-31G(d,p)//PM6-D3H+ DFT on the same proteins (Table 1, additivity check):**

| | Cα | Cβ | C′ | Hα | HN | N |
|---|---|---|---|---|---|---|
| Ubiquitin | 1.7 (0.77) | 2.5 (0.69) | 2.1 (0.72) | 0.6 (0.82) | 0.7 (0.85) | 4.4 (0.74) |
| GB3 | 1.6 (0.84) | 2.3 (0.60) | 2.3 (0.65) | 0.7 (0.82) | 0.8 (0.82) | 4.5 (0.78) |

RMSD in ppm, r in parentheses, all after linear regression.

**ProCS15 vs experiment, single CHARMM22/CMAP-optimised structure (Table 3):**

| Method | Cα | Cβ | C′ | Hα | HN | N |
|---|---|---|---|---|---|---|
| Full DFT (PM6-D3H+ structure) | 2.1 (0.62) | 2.8 (0.56) | 1.8 (0.85) | 0.4 (0.83) | 0.6 (0.81) | 4.0 (0.80) |
| ProCS15 | 1.7 (0.70) | 2.0 (0.50) | 1.7 (0.81) | 0.4 (0.77) | 0.6 (0.72) | 4.0 (0.79) |
| CamShift | 1.1 (0.85) | 1.3 (0.71) | 1.0 (0.81) | 0.3 (0.73) | 0.5 (0.69) | 3.0 (0.63) |
| Sparta+ | 0.7 (0.93) | 1.1 (0.82) | 0.8 (0.88) | 0.2 (0.86) | 0.4 (0.72) | 2.2 (0.81) |
| ShiftX2 | 0.5 (0.97) | 0.7 (0.91) | 0.5 (0.96) | 0.1 (0.97) | 0.3 (0.91) | 1.8 (0.88) |

**ProCS15 with NMR-derived ensembles for ubiquitin (Table 4):**

ProCS15 RMSD shrinks monotonically as the ensemble grows: Cα 1.7 → 1.0, Cβ 2.0 → 1.7, C′ 1.7 → 1.6, Hα 0.3 → 0.2, HN 0.6 → 0.5, N 3.7 → 3.5 (going 1UBQ X-ray → 2KOX 640-structure ensemble). ShiftX2 over the same series goes the other way for several channels: Cβ 0.4 → 0.6, C′ 0.4 → 0.7, N 1.3 → 1.8.

---

## Cross-validation / what the additivity check buys

The non-trivial validation in this paper is **ProCS15 vs full DFT on the same protein** (Table 1). Both use the same level of theory — OPBE/6-31G(d,p) — so any RMSD between them is *the cost of replacing a full quantum calculation with an additive lookup over tripeptide fragments*. The largest channel-wise gap is N (4.4 ppm against 4.0 ppm full-DFT for ubiquitin). For the proton channels the ProCS15-DFT gap is 0.0–0.1 ppm. So the additivity assumption costs little where the underlying DFT method is itself precise, and the method's main limitation is the underlying functional/basis set, not the fragment decomposition. The authors compare ProCS15 vs full-DFT vs experiment and conclude:

> "the accuracies of the chemical shifts computed using ProCS15 (based on linear regression of the chemical shifts, cf. Eq. (1)) are very similar to the corresponding DFT calculations using single Ubiquitin and GB3 structures."

---

## Authors' interpretation

The paper is dry-voiced. The substantive claims are restrained:

> "ProCS15 reproduces the chemical shielding values computed using PCM/OPBE/6-31G(d,p)//PM6-D3H+ for Ubiquitin and GB3 with RMSD values (after linear regression) of up to 2.5 ppm for carbon atoms, 0.8 ppm for hydrogen atoms, and 4.5 ppm for nitrogen."

> "ProCS15 can reproduce experimental chemical shifts with an overall accuracy that is similar to full DFT chemical shielding calculations for Ubiquitin and GB3."

The empirical-predictor concession is direct:

> "Comparison of ProCS15 to the empirical methods (CamShift through ShiftX2) generally show considerably lower RMSD of the empirical predictions for all atoms types … The r values are also considerably higher for the empirical methods than for ProCS15 for Cα and, especially, Cβ."

The load-bearing argument for using DFT-based predictors *despite* this is sensitivity to structural change:

> "it appears that the use of ensemble structures does not lead to a significant increase in accuracy compared to using a single structure for any of the empirical methods, in contrast to ProCS15 and CheShift-2."

> "the empirical NMR prediction methods tend to be significantly less sensitive to changes in protein structure compared to DFT-based chemical shift predictors or chemical shifts computed using QM methods."

The closing forward-looking line scopes future work to refinement:

> "we are now planning similar refinement studies using all backbone atoms and Cβ chemical shifts."

The argument is **not** "DFT predicts shifts more accurately than empirical methods" — it's "DFT responds to structural change in ways empirical methods don't, which makes DFT more useful for structure refinement even though raw RMSDs are larger."

---

## Terms-of-art (defined inline for re-use)

- **OPBE/6-31G(d,p)//PM6** — single-point GIAO chemical-shielding calculation at the OPBE density functional with the 6-31G(d,p) Pople basis set, on geometries optimised by the semi-empirical PM6 method. The "//" separates the level of theory used for the property (left) from the level used for the geometry (right).
- **GIAO** — gauge-including atomic orbitals, the standard formalism for making computed shieldings independent of the gauge origin of the magnetic vector potential.
- **PM6 / PM6-D3H+** — semi-empirical Hamiltonian for fast geometry optimisation; the D3H+ variant adds dispersion (D3) and improved hydrogen-bonding (H+) corrections.
- **CPCM / PCM** — implicit solvation models that treat the solvent as a polarisable continuum; CPCM is the conductor-like variant. ProCS15 uses ε = 78 (water).
- **CHARMM22/CMAP + GB/SA** — molecular-mechanics force field plus a generalised-Born / surface-area implicit solvation; used as an alternative pre-optimisation route.
- **AXA tripeptide / Ac-AXA-NMe** — capped-alanine tripeptide with central residue X, the unit of the lookup table.
- **FragBuilder / BASILISK / PHAISTOS** — the Jensen-group Python and Markov-chain stack used to construct AXA conformations and host the runtime ProCS15 module.
- **Random-coil correction** — sequence-specific reference shifts (Tamiola, Acar & Mulder 2010) subtracted before correlation analysis to remove the dominant sequence component and avoid inflated r values.
- **Ring-current dipole** — closed-form approximation Δσ ∝ (1 − 3cos²θ)/r³ for the long-range shielding contribution from aromatic rings.

---

## What's relevant to your kit (one paragraph, for the entry / outline, not for the deep-dive section)

ProCS15 is the closest published precedent to a per-residue DFT-grounded shielding predictor and is methodologically Stage-1-shaped (geometric features → DFT shielding via additivity, with empirical regression to experiment as the last step). The architecture differs from yours in three ways: (i) ProCS15 is **isotropic**, collapsing the shielding tensor to its trace; your work preserves T0/T1/T2. (ii) ProCS15 is **per-residue context-blind** — neighbour effects are limited to immediate i±1 side chains, with backbone dihedrals on flanking alanines fixed at canonical β-sheet values. (iii) ProCS15 calibrates by linear regression of σ → δ per atom type; you calibrate kernel-vs-DFT correlation. The ensemble-improves-RMSD-for-DFT-but-not-for-empirical-methods finding (Table 4) is directly relevant to your Stage 2 framing: averaging across structurally realistic ensembles benefits physics-grounded predictors and degrades empirical ones, which is the thesis you're making about per-frame-tensor-averaged kernels. The 2.35 M DFT shielding tensors (open archive on ERDA) is the same dataset you have ~1.9 M of from a prior collaboration.
