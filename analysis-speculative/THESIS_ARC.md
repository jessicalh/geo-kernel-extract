# Thesis Arc

Three steps: build the instrument, calibrate it, use it.

## 1. The instrument (what we built)

The C++ extraction engine. 10 calculators. Full T0/T1/T2
decomposition per kernel. Per-ring sparse contributions. The
spec (INDEX.md), the object model (ARCHITECTURE.md), the
constitution (CONSTITUTION.md). 287 passing tests.

The VTK visualisation that draws the tensors in 3D. The actual
picture of what a ring current shielding tensor looks like at an
atom near a phenylalanine ring. That illustration IS the thesis
in one image -- geometry producing a measurable physical quantity.

## 2. The calibration (what we learned about the geometric model)

Mechanical mutants (WT vs ALA) isolate the ring contribution.
DFT provides the reference tensor. 110 proteins, 69K atoms,
55 geometric kernels.

What we found:
- Per-element ridge R^2 = 0.818 (weighted). The geometric model
  explains 82% of the ring-deletion shielding perturbation.
- 28 validated realities grounding each kernel in NMR literature
  (distance falloff, element specificity, angular dependence,
  model convergence).
- BS, HM, and EFG are three views of the same geometry (R^2 ~0.77
  each for H, angular correlation cos ~0.93). The signal is
  geometric, not electronic.
- Per-element physics matches theory: H dominated by ring current
  (0.91), C by EFG (0.43), N multi-mechanism, O mixed.
- Ridge beats MLP. The coefficients are physical constants, not
  learned weights.

What we learned about the model's limits:
- HIE (imidazole, asymmetric 5-ring) is worst-modelled (R^2 = 0.06)
  -- the circular-loop approximation fails.
- Dispersion T2 doesn't decompose per-ring cleanly (many-body
  vertex sum).
- Valence, bond order, dipole add zero. The geometry encodes
  everything these scalars know.

## 3. The ensemble investigation (what we see in nature)

685 proteins. Conformational ensembles from well-tempered
metadynamics. Geometry-only extraction at 1.3s/frame. Experimental
RefDB shifts as ground truth.

What we're looking for:
- Do the kernel T0 values correlate with experimental shift
  residuals? Even 4% of variance is a finding at this scale.
- Does the distance dependence replicate (r^-3 for ring current)?
- Do BS, HM, and EFG co-vary across conformations (geometric
  equivalence holds under dynamics)?
- Where do they diverge (near-field, asymmetric rings)?
- Per-ring conformational gating -- which ring-atom interactions
  are stable vs dynamic?
- Emergent structure in the kernel-ensemble space that we didn't
  predict?

The analytical findings validate the instrument on real data. The
calibration validated it on DFT. The ensemble validates it on
experiment.

## + The GNN (what the geometry can't reach)

An equivariant GNN that takes the classical kernel outputs as
warm-start features and learns the residual through message passing.
The GNN's job is NOT to replace the geometric model. It's to learn
what the geometric model systematically misses: backbone conformation
effects, multi-ring interactions, charge redistribution, cooperative
H-bonding.

The connection: the analytical work tells the GNN which kernels carry
signal (feature importance), at what range (graph cutoff), through
which couplings (tensor product channels). The GNN's learned residual
should be interpretable against the analytical decomposition. If the
GNN improves prediction for N atoms by learning a backbone-torsion
correction, that's consistent with the calibration finding that N is
multi-mechanism.

The thesis claim: the geometric decomposition provides both a
complete analytical framework and an effective warm-start for
machine learning. The GNN's improvements are interpretable because
the baseline is interpretable.
