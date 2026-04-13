# New Array Assessment — Stage1Results Extraction

2026-04-13.  What the new arrays in Stage1Results tell us.

---

## Arrays new in Stage1Results (vs AzimuthalExtraction)

| Array | Shape | What |
|-------|-------|------|
| aimnet2_aim.npy | (N, 256) | AIMNet2 learned embedding |
| aimnet2_charges.npy | (N,) | AIMNet2 partial charges |
| aimnet2_efg_aromatic.npy | (N, 9) | EFG from AIMNet2 charges, aromatic decomp |
| aimnet2_efg_backbone.npy | (N, 9) | EFG from AIMNet2 charges, backbone decomp |
| aimnet2_efg.npy | (N, 9) | EFG from AIMNet2 charges, total |
| atom_sasa.npy | (N,) | Shrake-Rupley solvent-accessible surface area |
| dssp_chi.npy | (N, 12) | chi1-4 cos/sin/exists |
| dssp_hbond_energy.npy | (N, 4) | H-bond energies |
| dssp_ss8.npy | (N, 8) | 8-class secondary structure one-hot |
| orca_diamagnetic.npy | (N, 9) | WT absolute DFT diamagnetic shielding |
| orca_paramagnetic.npy | (N, 9) | WT absolute DFT paramagnetic shielding |

Arrays orca_total.npy (WT DFT total) was already present.
No arrays were removed from the old extraction.

## Every new scalar adds zero

Tested on 200 proteins.  Each feature alone as ridge predictor of
delta T2, and interacted with EFG_aro T2:

| Feature | Alone R² (all elem) | Best EFG interaction |
|---------|-------------------|-----------------------|
| SASA | < 0.001 | +0.007 for EFG+SASA on O |
| SS8 | 0.001-0.002 | +0.037 for EFG+SS8 on N |
| H-bond energy | < 0.001 | +0.003 |
| AIMNet2 charge | < 0.001 | ~ 0 |
| MOPAC charge | < 0.001 | ~ 0 |

The SS8+EFG interaction for N is the only thing above noise.
Consistent with nitrogen's backbone-angle sensitivity (Yao 2010) —
secondary structure captures N atom distribution, not a mutation
effect.

Conclusion: the geometric T2 kernels already encode everything
these scalars know, in the context of mechanical mutants where
backbone is unchanged.

## AIMNet2 EFG is orthogonal to Coulomb EFG

| Comparison | H | C | N | O |
|-----------|---|---|---|---|
| ff14SB ↔ MOPAC cos | 0.989 | 0.987 | 0.991 | 0.991 |
| ff14SB ↔ AIMNet2 cos | 0.338 | 0.338 | 0.337 | 0.336 |
| MOPAC ↔ AIMNet2 cos | 0.338 | 0.338 | 0.338 | 0.335 |

ff14SB and MOPAC produce nearly identical aromatic EFG tensors
(cosine 0.99).  AIMNet2 is orthogonal to both (cosine ~0.34, near
random in 5D = 0.45).

Despite being angular orthogonal, AIMNet2 EFG gets R²=0.53 for H
(meaningful).  But combined with ff14SB+MOPAC (all three together),
R² = 0.77 — same as MOPAC alone (0.77).  AIMNet2 adds nothing
beyond MOPAC for the mutation delta.

AIMNet2's value is for Stage 2 (ensemble): it runs at 0.17s/frame
vs MOPAC at ~10 min/frame.  The EFG from AIMNet2 charges captures
different angular structure that may matter for conformational
variation.

## Charge comparison: AIMNet2 vs MOPAC

| Element | Correlation | AIM mean | MOP mean |
|---------|------------|----------|----------|
| H | 0.892 | +0.122 | +0.200 |
| C | 0.980 | +0.043 | +0.017 |
| N | 0.247 | -0.356 | -0.516 |
| O | 0.530 | -0.441 | -0.624 |

AIMNet2 and MOPAC charges agree well for H and C but diverge
sharply for N (r=0.25) and O (r=0.53).  These are the elements
where shielding is hardest to model.  Stage 2/3 finding.

## SASA distribution

| Element | Mean Å² | Std | Frac > 0 |
|---------|---------|-----|----------|
| H | 7.96 | 8.07 | 72% |
| C | 1.17 | 2.18 | 34% |
| N | 1.56 | 3.40 | 31% |
| O | 13.94 | 13.92 | 72% |

Most C and N atoms are buried (SASA = 0).  H and O are mostly
exposed.  SASA as a scalar doesn't predict delta T2, but for
ensemble work it captures solvent exposure variation across frames.

## Dia/para decomposition

See separate note: dia_para_decomposition.md.
Validated to 0.003 ppm.  Key finding: dia and para perturbations
cancel (~7+7→1 ppm).  Kernels see the net, not the channels.
