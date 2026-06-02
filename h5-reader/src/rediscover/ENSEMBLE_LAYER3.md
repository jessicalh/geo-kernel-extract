# Ensemble Layer 3

Joint per-stratum ridge on the corrected emitted substrate. This is a within-protein signal-capture diagnostic, not a population prediction claim.

## Lead

- Physics-only between-atom LOAO R2 is positive in 8/12 target-stratum rows. Best physics-only row: sigma_iso HN R2=0.5332.
- Physics + AIMNet2 is positive in 9/12 rows. Best ceiling row: sigma_iso HN R2=0.6450.
- McConnell's joint drop-one result: 5 threshold-positive rows across the 12 between-led tests. This is the fair joint test; standalone McConnell remains a weaker read.
- Verdict bucket for the layer-3 physics ensemble: form-recovered-scale-fitted where between LOAO R2 is positive and stable; otherwise can't-make-it-work-for-now.

## Joint LOAO R2

Headline is between-atom LOAO R2. Within-frame R2 is a train/test frame diagnostic after train-only per-atom centering. Static sigma_iso is between-led.

| target | stratum | n_atoms_between_physics | variance_share_between_physics | between_LOAO_R2_physics | between_LOAO_R2_jackknife_se_physics | between_LOAO_R2_plus_aimnet2 | between_LOAO_R2_jackknife_se_plus_aimnet2 | aimnet2_gap_between_R2 | within_frame_R2_physics | within_frame_R2_plus_aimnet2 | within_N_eff_median_physics | median_lag1_rho_physics | thin_flag_physics | physics_verdict | plus_aimnet2_verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T2 | N | 54 | 0.4647 | 0.0249 | 0.1528 | 0.0739 | 0.2413 | 0.0489 | 0.5525 | 0.6245 | 2.2200e+04 | 0.0770 |  | can't-make-it-work-for-now | can't-make-it-work-for-now |
| T2 | CA | 54 | 0.5469 | 0.2015 | 0.1203 | 0.1701 | 0.1400 | -0.0313 | 0.2959 | 0.4211 | 2.0800e+04 | 0.0903 |  | form-recovered-scale-fitted | form-recovered-scale-fitted |
| T2 | C | 54 | 0.3330 | -0.4012 | 0.2308 | -0.2251 | 0.3179 | 0.1761 | 0.5826 | 0.6971 | 2.3612e+04 | 0.0547 |  | can't-make-it-work-for-now | can't-make-it-work-for-now |
| T2 | O | 54 | 0.3431 | -0.1184 | 0.0622 | -0.2856 | 0.1838 | -0.1673 | 0.5722 | 0.6342 | 2.2482e+04 | 0.0537 |  | can't-make-it-work-for-now | can't-make-it-work-for-now |
| T2 | HN | 52 | 0.5824 | 4.9275e-04 | 0.0778 | 0.0610 | 0.0823 | 0.0605 | 0.3192 | 0.4089 | 1.7877e+04 | 0.1575 |  | can't-make-it-work-for-now | can't-make-it-work-for-now |
| T2 | HA | 58 | 0.7077 | 0.0852 | 0.0790 | 0.0218 | 0.0956 | -0.0634 | 0.3498 | 0.4159 | 1.6975e+04 | 0.2332 |  | can't-make-it-work-for-now | can't-make-it-work-for-now |
| sigma_iso | N | 54 | 0.7707 | -0.4841 | 0.9210 | 0.3833 | 0.6254 | 0.8674 | 0.1664 | 0.3770 | 2.3273e+04 | 0.0500 |  | can't-make-it-work-for-now | form-recovered-scale-fitted |
| sigma_iso | CA | 54 | 0.5425 | 0.0180 | 0.2408 | 0.5253 | 0.1450 | 0.5073 | 0.1436 | 0.5895 | 2.4355e+04 | 0.0338 |  | can't-make-it-work-for-now | form-recovered-scale-fitted |
| sigma_iso | C | 54 | 0.0811 | -0.5034 | 0.2661 | -0.7291 | 0.3267 | -0.2256 | 0.2123 | 0.6460 | 2.5823e+04 | -0.0011 |  | can't-make-it-work-for-now | can't-make-it-work-for-now |
| sigma_iso | O | 54 | 0.2349 | 0.1046 | 0.1509 | 0.2944 | 0.1648 | 0.1899 | 0.2450 | 0.4982 | 2.3850e+04 | 0.0358 |  | form-recovered-scale-fitted | form-recovered-scale-fitted |
| sigma_iso | HN | 52 | 0.4378 | 0.5332 | 0.1408 | 0.6450 | 0.0933 | 0.1118 | 0.4326 | 0.7062 | 2.0788e+04 | 0.0920 |  | form-recovered-scale-fitted | form-recovered-scale-fitted |
| sigma_iso | HA | 58 | 0.5264 | 0.1429 | 0.2528 | 0.0743 | 0.2847 | -0.0686 | 0.2764 | 0.3935 | 2.0980e+04 | 0.1491 |  | form-recovered-scale-fitted | can't-make-it-work-for-now |

## Physics Contribution Shares

Shares are drop-one positive-delta shares within the physics-only joint fit. Each cell is `share (dR2 delta)`, so negative or redundant blocks can have zero share.

| target | stratum | ring | charge | mcconnell | field_buckingham | efg |
| --- | --- | --- | --- | --- | --- | --- |
| T2 | N | 0.0000 (dR2 -0.0850) | 0.0000 (dR2 -0.0298) | 0.9904 (dR2 0.1875) | 0.0096 (dR2 0.0018) | 0.0000 (dR2 -0.0094) |
| T2 | CA | 0.0070 (dR2 0.0024) | 0.0175 (dR2 0.0061) | 0.9755 (dR2 0.3406) | 0.0000 (dR2 -0.0042) | 0.0000 (dR2 -0.0270) |
| T2 | C | 0.0000 (dR2 -0.0149) | 0.0000 (dR2 -0.1747) | 0.0000 (dR2 -0.0474) | 1.0000 (dR2 0.1182) | 0.0000 (dR2 -0.1576) |
| T2 | O | 0.0048 (dR2 9.5341e-04) | 0.0282 (dR2 0.0056) | 0.5773 (dR2 0.1145) | 0.3897 (dR2 0.0773) | 0.0000 (dR2 -0.1001) |
| T2 | HN | 0.0000 (dR2 -0.0213) | 0.5426 (dR2 0.0402) | 0.0000 (dR2 -0.0524) | 0.1357 (dR2 0.0101) | 0.3217 (dR2 0.0238) |
| T2 | HA | 0.4666 (dR2 0.0228) | 0.5334 (dR2 0.0261) | 0.0000 (dR2 -0.0686) | 0.0000 (dR2 -0.0508) | 0.0000 (dR2 -0.0335) |
| sigma_iso | N | 0.0000 (dR2 -0.1108) | 0.0000 (dR2 -0.0192) | 0.0000 (dR2 -0.2156) | 0.0000 (dR2 -0.0827) | 0.0000 (dR2 -0.1036) |
| sigma_iso | CA | 0.3140 (dR2 0.1528) | 0.0461 (dR2 0.0224) | 0.0000 (dR2 -0.0114) | 0.0000 (dR2 -0.1064) | 0.6399 (dR2 0.3113) |
| sigma_iso | C | 0.0000 (dR2 -0.2060) | 1.0000 (dR2 0.0071) | 0.0000 (dR2 -0.2844) | 0.0000 (dR2 -0.0418) | 0.0000 (dR2 -0.0759) |
| sigma_iso | O | 0.0000 (dR2 -0.0600) | 0.0000 (dR2 -0.0413) | 0.5403 (dR2 0.3003) | 0.4597 (dR2 0.2555) | 0.0000 (dR2 -0.0880) |
| sigma_iso | HN | 0.0000 (dR2 -0.0358) | 0.0664 (dR2 0.0193) | 0.0000 (dR2 -0.0351) | 0.2642 (dR2 0.0768) | 0.6694 (dR2 0.1946) |
| sigma_iso | HA | 0.0000 (dR2 -0.2576) | 0.0000 (dR2 -0.0642) | 0.6056 (dR2 0.3406) | 0.0000 (dR2 -0.0229) | 0.3944 (dR2 0.2218) |

## McConnell Verdict

Drop-one values ask whether valid-source per-category McConnell improves the already-joint physics fit. The form-recovered-scale-fitted bucket requires positive dR2 above 0.02 and above the delete-atom jackknife SE.

| target | stratum | mcconnell_delta_between_R2 | mcconnell_delta_between_R2_se | mcconnell_share_between | mcconnell_verdict_bucket |
| --- | --- | --- | --- | --- | --- |
| T2 | N | 0.1875 | 0.1237 | 0.9904 | form-recovered-scale-fitted |
| T2 | CA | 0.3406 | 0.0776 | 0.9755 | form-recovered-scale-fitted |
| T2 | C | -0.0474 | 0.2022 | 0.0000 | can't-make-it-work-for-now |
| T2 | O | 0.1145 | 0.0742 | 0.5773 | form-recovered-scale-fitted |
| T2 | HN | -0.0524 | 0.0578 | 0.0000 | can't-make-it-work-for-now |
| T2 | HA | -0.0686 | 0.0456 | 0.0000 | can't-make-it-work-for-now |
| sigma_iso | N | -0.2156 | 0.1784 | 0.0000 | can't-make-it-work-for-now |
| sigma_iso | CA | -0.0114 | 0.1171 | 0.0000 | can't-make-it-work-for-now |
| sigma_iso | C | -0.2844 | 0.2306 | 0.0000 | can't-make-it-work-for-now |
| sigma_iso | O | 0.3003 | 0.1633 | 0.5403 | form-recovered-scale-fitted |
| sigma_iso | HN | -0.0351 | 0.0760 | 0.0000 | can't-make-it-work-for-now |
| sigma_iso | HA | 0.3406 | 0.2121 | 0.6056 | form-recovered-scale-fitted |

## Honest Bottom Line

The recovered kernels jointly capture real backbone signal only where the between-atom LOAO rows are positive with tolerable jackknife uncertainty. Ring, charge, field/Buckingham, and EFG compete in the joint basis; McConnell gets its fair valid-source per-category test as a component rather than as a standalone total-target fit. AIMNet2 remains a ceiling: embedding PCs, AIMNet2 charge, and CRG are learnable upper-bound features, not recovered laws.

Ring and charge remain visible but are not universal joint drop-one leaders in this one-protein layer-3 fit. Charge is clearest on HA/HN T2-style rows, ring appears in HA T2 and CA sigma_iso, while McConnell, field/Buckingham, and EFG carry several between-led deltas once all kernels compete in the same ridge.

This is n=1 protein evidence. The report is per stratum, correlate-not-match, and does not make population claims.

## Self-Audit

- Substrate: `/tmp/rediscover-corrected-backbone-snapshot-1p9j`.
- Broad bond cutoff unique values: `[10.0]`.
- McConnell near-field ratio unique values: `[0.5]`.
- Valid McConnell category counts: `{'PeptideCO': 2134696, 'PeptideCN': 2126586, 'SidechainCO': 683055}`.
- Python read emitted CSV/NPY features only; it did not open `trajectory.h5`, run ORCA, call the C++ emitter, or recompute kernels.
- Frozen basis audit `|C.T C - I|max`: 1.110223e-16.
- Ridge alpha: 10.0000; embedding PCs: 32 of 256.
- Artifacts: `/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis/ensemble_layer3`.
