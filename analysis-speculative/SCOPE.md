# Scope and Assumptions

## The thesis deliverable (non-negotiable)

Two things at the end:

1. **Analytical physics findings** from ensemble kernel data vs
   experimental shifts. The calibration proved the instrument works
   (28 realities, R^2=0.818 on mutants). The ensemble analysis shows
   the instrument sees real physics in nature. Even weak signals
   (4% of residual variance) are findings if they replicate known
   physics (r^-3 falloff, element specificity, ring-type dependence)
   at 685-protein scale. Emergent patterns (clusters, covariances,
   conditional structure) that nobody has seen are the prize.

2. **A working equivariant GNN** that predicts NMR shielding tensors
   from structure, using the classical kernels as warm-start features.
   The analytical findings and the GNN findings must corroborate each
   other. The GNN learns what the kernels can't reach. The analytical
   work shows what the kernels DO reach.

Neither alone is a thesis. Together they are.

## Hardware (what we actually have)

| Resource | Spec | Constraint |
|---|---|---|
| Primary GPU | 4x RTX 5090 | GNN training, batch extraction |
| CPU RAM | 128 GB | Limits in-memory analysis to ~100 GB working set |
| Disk | Sufficient for fleet | 685 x 500 MB = ~340 GB extracted features |
| R machine | ~2 TB RAM, short loan | Publication figures, massive cross-joins |
| Time | 3 months from 2026-04-10 | July 2026 deadline |

## What we compute ourselves (no lookup tables, no hunting old data)

The extraction engine IS the instrument. If we need a quantity, we
compute it from geometry during extraction. Anything that takes
>~2s/frame is too slow for ensemble processing.

Available at ~1.3s/frame (geometry-only):
- 4 ring current calculators (BS, HM, PQ, Disp) with full T0/T1/T2
- Ring susceptibility
- McConnell bond anisotropy (5 categories)
- Coulomb E-field and EFG (topology charges, backbone/aromatic split)
- H-bond geometry
- DSSP backbone angles
- Per-ring sparse contributions (59 columns per atom-ring pair)

NOT available without MOPAC (~10 min/frame):
- MOPAC charges, bond orders, orbital populations

MOPAC charges are NOT needed. The calibration showed topology EFG_aro
and MopacEFG_aro have angular correlation cos 0.9296 vs 0.9323 with
the DFT target (difference: 0.003). The EFG tensor direction is
determined by the geometric T-tensor, not by charge magnitudes. BS,
HM, and EFG all score R^2 ~0.77 for H independently -- three views
of the same ring geometry. See EFG_GEOMETRY_DOMINANCE.md.

## SPECULATIVE assumptions (will test, may not hold)

1. Ensemble-averaged kernel T0 correlates with experimental shift
   residuals after residue-type centering. ASSUMPTION: the baseline
   subtraction is clean enough that geometric perturbations emerge.
   May fail if: residue-type means are too noisy, atom matching is
   too lossy, or the geometric contribution is below noise floor.

2. Per-ring contributions across frames show conformational gating.
   ASSUMPTION: ring-atom distances vary enough across 5005 frames
   that some rings gate on/off. May fail if: the proteins are too
   rigid, or the 6-pose sampling misses the interesting motions.

3. Kernel covariance across frames reveals physical coupling.
   ASSUMPTION: BS and EFG co-vary because the same structural motion
   drives both. May fail if: the kernels are effectively independent
   (different structural modes dominate each), or the covariance is
   trivially explained by distance.

4. Clusters in kernel-ensemble profiles correspond to interpretable
   physics. ASSUMPTION: atoms with similar kernel distributions share
   structural features. May fail if: the distributions are too noisy,
   or the clusters are just element/residue-type.

5. The GNN learns residuals that the analytical work predicts.
   ASSUMPTION: the analytical findings (which kernels matter where)
   become feature-importance priors for the GNN. May fail if: the
   GNN finds completely different structure, or the kernel features
   are redundant with learned message-passing features.

6. BS, HM, and EFG divergence across conformations reveals
   near-field non-dipolar effects. ASSUMPTION: the three models
   agree in the far-field (calibration-verified) but diverge in
   specific near-field geometries. May fail if: they co-vary
   perfectly even in the near-field, or the divergence is below
   noise floor at ensemble scale.

## What we do NOT assume

- That R^2 will be high. We are looking for physics, not fit quality.
- That the GNN will beat classical methods. It might just corroborate.
- That all 685 proteins will be usable. Quality filters will shrink N.
- That the 6 initial-samples poses are sufficient. They're the control
  group. The full trajectories are the experiment.
- That we know what we'll find. The emergent layer is real.
