# Master Chart: The T2 Shielding Perturbation — What the Tool Sees

2026-04-13.  720 proteins, 446,006 atoms.  All numbers final.

---

```
THE T2 SHIELDING PERTURBATION FROM AROMATIC RING DELETION
(WT - ALA mechanical mutant, DFT reference)
│
│  The DFT delta decomposes into diamagnetic + paramagnetic
│  contributions that nearly cancel: ~7 + ~7 → ~1-2 ppm.
│  Para/dia variance ratio: H=1.02  C=1.06  N=1.07  O=1.20
│  (Saito 2010: heavier atoms → more paramagnetic)
│
├── HYDROGEN  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
│   R²=0.856 (norm), LPOCV=0.844, 20 dims, NO nonlinear signal
│   Diamagnetic-dominated.  Linear.  The tool sees H clearly.
│   │
│   ├── Ring current [R²=0.782, cos→target=0.70]
│   │   THE dominant dimension.  17 kernels, 8-9 independent views.
│   │   Physics: pi-electron circulation → magnetic dipole field
│   │   Equation: G_ab = -n_b B_a  (Pople 1956, Johnson-Bovey 1958)
│   │   Distance: r⁻³.  Measured slope: -3.04 (theory: -3)
│   │   BS≡HM for same ring (cos=1.00, Case 1995 extended to T2)
│   │   PHE↔TRP_pyrrole cos=0.36 (different ring geometry = independent)
│   │   Ring type diversity provides most of the 20 dimensions
│   │
│   ├── EFG [R²=0.701 MOPAC / 0.660 ff14SB, cos→target=0.93]
│   │   Highest angular alignment with DFT target of any group.
│   │   Physics: partial charges → electric field gradient tensor
│   │   Equation: V_ab = Σ q_j K_ab(d_j)  (Buckingham 1960)
│   │   Pure T2 (traceless symmetric, Gauss's law)
│   │   ff14SB↔MOPAC cos=0.58 (same eq, different charges)
│   │   Enters forward selection FIRST (0.70 in one kernel)
│   │
│   ├── Dispersion [R²=0.387, cos→target=0.85]
│   │   Highest cos→target but spatially sparse (r⁻⁶ cutoff at 5Å)
│   │   Physics: London van der Waals from ring vertices
│   │   Equation: C6(3d_ad_b/r⁸ - δ_ab/r⁶)  (London 1937)
│   │   DispChi = disp_scalar × RingSusc T2 (proximity-gated angular)
│   │
│   ├── Bond aniso [R²=0.107]
│   │   Minor for H.  Hydrogen doesn't feel bond anisotropy directly
│   │   (bonded to only one atom with short bond).
│   │
│   └── Controls: PQ=0.030, solvation=0.021, HBond~0
│
│
├── CARBON  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
│   R²=0.484 (norm), LPOCV=0.456, 6 dims, +0.128 nonlinear signal
│   Paramagnetic contribution growing.  EFG-dominated.
│   Charge polarisation gap: +0.197 (MOPAC vs geo-only)
│   RF recovers +0.128 without MOPAC → gating partially bridges gap
│   │
│   ├── EFG [R²=0.380 MOPAC / 0.225 ff14SB]
│   │   Dominant.  The charge redistribution from ring removal changes
│   │   the EFG angular pattern — MOPAC captures this, ff14SB doesn't.
│   │   MOPAC-ff14SB gap = +0.155: the charge-polarisation dimension.
│   │   This is specifically a carbon phenomenon in this data.
│   │   AIMNet2 EFG: cos=0.34 (orthogonal), R²=0.07 (wrong projection)
│   │
│   ├── Ring current [R²=0.121]
│   │   Minor.  Carbon shielding is paramagnetic-term-dominated
│   │   (Sahakyan & Vendruscolo 2013) — ring current magnetic field
│   │   affects the diamagnetic term which is secondary for C.
│   │
│   ├── Dispersion [R²=0.115]
│   │   Enters forward selection rank 4-5.  Short-range angular
│   │   structure captures ring-proximity effects for C.
│   │
│   ├── Bond aniso [R²=0.025, MC_aromatic enters rank 3]
│   │   The removed aromatic bonds contribute.  MC_aromatic_total is
│   │   the only bond kernel that matters (backbone unchanged).
│   │
│   └── Nonlinear: RF=0.581 vs ridge=0.453 (+0.128)
│       Kernel interactions carry signal.  Which kernel matters depends
│       on the specific distance/angle geometry.  Gating helps because
│       the EFG angular pattern changes with distance to the ring in a
│       way that linear weighting cannot capture.
│       Paramagnetic term: 1/ΔE denominators couple geometry nonlinearly
│       (Ramsey 1950 — to be grounded at write time)
│
│
├── NITROGEN  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
│   R²=0.267 (norm), LPOCV=0.172, 3 blurred dims, +0.169 nonlinear
│   Paramagnetic-dominated.  Multi-mechanism.  The tool sees N poorly.
│   Largest overfit gap (+0.074): weights are protein-specific.
│   │
│   ├── Ring current [R²=0.124]
│   │   Enters rank 3.  Collective effect from 17 kernels spanning
│   │   geometric space — no single ring current kernel dominates.
│   │   Isotropic ring current on ¹⁵N is negligible (<0.6%, Han 2011)
│   │   but anisotropic (T2) effect is measurable.
│   │
│   ├── EFG [R²=0.064]
│   │   Weak.  ff14SB and MOPAC give same R² (0.07 vs 0.06).
│   │   Charge polarisation doesn't help N — the electronic
│   │   perturbation at N is not captured by point-charge EFG.
│   │   N lone pair and paramagnetic dominance mean the electrostatic
│   │   perturbation has a different effect than for H or C.
│   │
│   ├── Quadrupole [R²=0.063, enters rank 2]
│   │   More important for N than any other element.  PQ_total enters
│   │   forward selection at rank 2 (+0.043).  The r⁻⁵ quadrupole
│   │   from the ring's pi-cloud provides angular information at the
│   │   distances where N atoms sit (Stone 2013, measured slope -5.05).
│   │   PQ↔BS cos=0.45 (random) — genuinely independent dimension.
│   │
│   ├── Dispersion [R²=0.068 after norm, 0.015 raw]
│   │   Normalisation reveals dispersion for N (4.5× increase).
│   │   Enters 5 of last 10 forward selection steps.
│   │
│   ├── MOPAC bond [R²=0.033, enters rank 4]
│   │   MopacMC_total enters forward selection at rank 4.
│   │   Peptide bond partial double-bond character (Wiberg 1.3-1.5)
│   │   couples N electronically to the mutation site.
│   │   This is the only element where bond-order weighting matters.
│   │
│   ├── Bond aniso [R²=0.029, enters rank 8]
│   │   MC_aromatic_total — same as C but weaker.
│   │
│   ├── ¹⁵N CSA context (Hall & Fushman 2006):
│   │   Site-specific ¹⁵N CSA varies by 120 ppm in protein GB3.
│   │   Phe52 and Trp43 are explicit aromatic-adjacent outliers.
│   │   They measured the variation but could not explain it.
│   │   Our decomposition provides the mechanistic explanation.
│   │
│   └── Nonlinear: RF=0.316 vs ridge=0.148 (+0.169)
│       LARGEST nonlinear signal of any element.
│       The 3 blurred PCA dimensions are a linear projection of a
│       higher-dimensional nonlinear manifold.  RF accesses the
│       interaction effects: BS×PQ, EFG×MC, distance×angle combinations
│       that produce different shielding perturbations than their sum.
│       Grounding: paramagnetic shielding term has product terms
│       (geometry₁ × geometry₂ / ΔE) — Ramsey 1950.
│       ¹⁵N is paramagnetic-dominated (Saito 2010), making these
│       interaction effects dominant.
│
│
├── OXYGEN  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
│   R²=0.304 (norm), LPOCV=0.213, 12 dims, +0.013 nonlinear (marginal)
│   Most paramagnetic (para/dia ratio 1.20).  Dispersion-dominated
│   after normalisation.  Linear is nearly sufficient despite para dom.
│   │
│   ├── Dispersion [R²=0.234 norm / 0.058 raw — 4× from normalisation]
│   │   THE headline finding for oxygen.  After normalisation strips
│   │   the r⁻⁶ magnitude scale, dispersion becomes the dominant group
│   │   (0.234 vs ring current 0.195 vs EFG 0.158).
│   │   8 kernels with within-group cos range [0.33, 0.97] — rich
│   │   angular diversity from ring-type-specific vertex patterns.
│   │   DispChi for carbonyl O captures the ring's angular effect on
│   │   the C=O bond that raw magnitude hid.
│   │   4 consecutive dispersion steps in normalised forward selection.
│   │
│   ├── Ring current [R²=0.195]
│   │   Second group.  O atoms (mainly carbonyl) are further from
│   │   rings on average but still see the magnetic dipole field.
│   │
│   ├── EFG [R²=0.158 ff14SB / 0.150 MOPAC]
│   │   Third.  MOPAC does NOT help for O (slightly worse than ff14SB).
│   │   The charge-polarisation dimension is irrelevant for oxygen.
│   │
│   ├── Bond aniso [R²=0.041]
│   │   MC_aromatic_total enters at rank 13.  The dominant O bond
│   │   (C=O, largest Δχ) doesn't change in a mechanical mutant —
│   │   consistent with zero McConnell contribution from unchanged bonds.
│   │
│   ├── Quadrupole [R²=0.038, enters rank 7-8]
│   │   Modest.  Pi-quadrupole from rings at the distances where O sits.
│   │
│   ├── Near-field complexity: 8 predictive dims at 0-4 Å
│   │   Multiple short-range mechanisms compete: dispersion, ring current,
│   │   quadrupole all contribute in the near field.  This is where
│   │   the functional forms of the kernels diverge from the true
│   │   electron distribution.
│   │
│   └── Nonlinear: RF=0.215 vs ridge=0.203 (+0.013)
│       Marginal.  Despite being the most paramagnetic element (ratio
│       1.20), oxygen's nonlinear signal is small.  Possible reason:
│       oxygen has fewer competing mechanisms than nitrogen (dispersion
│       dominates after normalisation, so there are fewer interaction
│       terms).  The paramagnetic nonlinearity may be present but
│       masked by the strong dispersion linear signal.
│
│
├── CROSS-ELEMENT STRUCTURE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
│
│   Angular independence (group-to-group |cos|) is ELEMENT-INVARIANT:
│   Ring current ↔ EFG:      0.44 (near-random, all elements)
│   Ring current ↔ Bond:     0.40 (independent, all elements)
│   Ring current ↔ PQ:       0.45 (random — different multipolar orders)
│   EFG ↔ Dispersion:        0.53 (moderate — charge-C6 correlation)
│   ff14SB EFG ↔ MOPAC EFG:  0.58-0.67 (same eq, different charges)
│
│   The GEOMETRY of the kernel space is fixed.
│   The PROJECTION onto the DFT target is element-specific.
│
│   Normalisation is a physics operation:
│   - Does NOT change cosines (angular structure is intrinsic)
│   - DOES change R² (strips magnitude, reveals angular signal)
│   - H loses (magnitude was signal), N/O gain (magnitude was noise)
│   - Dispersion for O: 0.058 → 0.234 (4× — the biggest single effect)
│
│   Nonlinear signal follows paramagnetic ordering:
│   N (+0.169) > C (+0.128) > O (+0.013) > H (+0.002)
│   Grounding: Ramsey 1950 paramagnetic term, Saito 2010 element dep.
│
│
├── PROVEN ZEROS ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
│
│   Scalar features that add nothing to the kernel T2 prediction:
│   MOPAC valence, bond order, molecular dipole, SASA, DSSP SS8,
│   H-bond energy, AIMNet2 charge (scalar), MOPAC charge (scalar)
│   — all < 0.001.  The kernels encode everything these know.
│
│   Many are correctly zero BECAUSE of mechanical mutant design:
│   backbone unchanged → DSSP/HBond/SASA don't change.
│   For conformational variation (Stage 2) these may matter.
│
│   H-bond total: R² = 0.002 (negative control — correctly zero)
│   DeltaAPBS solvation: R² = 0.005 (2% of direct EFG — negligible)
│
│
├── DIA/PARA DECOMPOSITION ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
│
│   Parsed from raw orca .out files, validated to 0.003 ppm vs C++.
│   WT dia and para perturbations nearly cancel:
│     H: 6.44 + 6.49 → 0.78 ppm total
│     C: 7.42 + 7.67 → 1.54
│     N: 6.96 + 7.27 → 1.58
│     O: 6.55 + 7.15 → 2.13
│
│   Kernels predict the NET (R²=0.70 H) not the channels (<0.05).
│   Classical geometric kernels compute total interaction, not the
│   QM decomposition into dia/para.  Physically correct.
│
│
└── WHAT THE TOOL CANNOT SEE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    Per-protein ceiling: R²=0.81 within a protein, R²=0.35 across.
    The gap is global electrostatic environment (protein shape, bulk
    dielectric, crystal packing) that local geometric kernels cannot
    capture.  Motivates ensemble conformational sampling (Stage 2).

    Near-field accuracy: BS/HM wire-loop and surface-integral models
    approximate the true pi-electron current density.  At r < 4 Å
    the models diverge from QM.  Per-ring kernels partially address
    this but the functional form is fixed.

    HIE (imidazole): worst ring type.  Self-fit R² = 0.062.  The
    symmetric circular-loop BS approximation fails for the asymmetric
    5-membered ring with two nitrogen positions.

    Nitrogen's 3 blurred dimensions: the tool sees nitrogen in the
    corner of its eye.  The paramagnetic multi-mechanism response
    exceeds what geometric kernels can resolve linearly.  RF shows
    the information IS there (+0.169) but requires nonlinear extraction.
```

---

## Citation map

Every physics assertion maps to a citation or a measurement:

| Assertion | Source |
|-----------|--------|
| Ring current = magnetic dipole | Pople 1956, Johnson & Bovey 1958 |
| Surface integral model | Haigh & Mallion 1979 |
| BS-HM calibrated equivalence | Case 1995, extended to T2 here |
| EFG = Buckingham effect | Buckingham 1960 |
| Bond anisotropy | McConnell 1957 |
| Quadrupole from T-tensor | Stone 2013, Ch. 3 |
| Dispersion | London 1937 |
| H = diamagnetic dominated | Saito et al. 2010 |
| C = EFG dominated (isotropic) | Sahakyan & Vendruscolo 2013 |
| ¹⁵N CSA variability, aromatic outliers | Hall & Fushman 2006 |
| ¹⁵N isotropic ring current negligible | Han et al. 2011 |
| Paramagnetic term nonlinearity | Ramsey 1950 |
| BS slope -3.04, PQ slope -5.05 | Measured, this work |
| BS-HM ratio 1.009 ± 0.037 | Measured, this work |
| Per-element T2 decomposition | Novel, this work |
| Nonlinear signal ∝ paramagnetic dominance | Novel, this work |
| Dispersion dominant for O after normalisation | Novel, this work |
| Dia/para cancellation in T2 delta | Novel, this work |
