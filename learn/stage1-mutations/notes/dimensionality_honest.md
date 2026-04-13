# Dimensionality — Honest Account

2026-04-13.  300 proteins from Stage1Results.  No pre-committed
narrative.

---

## The raw eigenspectrum is universal

Every element shows the same raw structure: 5 eigenvalues at
16-24% each (one L=2 source, ~85% total), sharp drop to ~2%
(everything else is noise-level in raw space).

| Element | Top 5 cum% | PC6 | Interpretation |
|---------|-----------|-----|---------------|
| H | 84.6% | 3.3% | One dominant T2 source |
| C | 89.2% | 2.2% | Same |
| N | 87.0% | 2.7% | Same |
| O | 89.2% | 2.3% | Same |

This is the "where are the rings" information.  Every kernel knows
it.  It's shared, not predictive.

## Normalization reveals element-specific complexity

After per-protein normalization strips the magnitude structure,
each element reveals its own predictive dimensionality:

| Element | Norm plateau | Near-field dims | Far-field dims | R² |
|---------|-------------|----------------|---------------|------|
| H | 20 | 3 | 3 | 0.928 |
| C | 6 | 3 | 3 | 0.562 |
| N | 3 | 4 | 3 | 0.380 |
| O | 12 | 8 | 3 | 0.382 |

### Hydrogen: 20 dimensions, all contributing

The normalized eigenspectrum for H is nearly flat (3.1%, 2.9%,
2.7%...).  No clean blocks.  Ridge weight spreads across 96 PCs
for 90% coverage.  The prediction lives in the fine angular
disagreements between ALL kernel families.  20 PCA dimensions
contribute to out-of-sample R² before plateauing.

This is why H works so well: the tool sees hydrogen in 20
independent angular dimensions, each carrying a small slice of
the total R² = 0.928.

### Carbon: 6 dimensions

Tighter than H.  Dominated by EFG (Dimension 2) + charge
polarization.  The 6 predictive PCs likely correspond to the EFG
T2 (5 components) plus a scale/polarization factor.  Carbon needs
fewer angular dimensions because one mechanism dominates.

### Nitrogen: 3 dimensions, blurred

The simplest normalized space despite being the "hardest" element.
No single mechanism dominates — 5 families each contribute 0.015-
0.080 in forward selection.  The 3 predictive dimensions are not
clean physics channels; they're blurred mixtures of all mechanisms.
Near-field goes to 4 dimensions (the multi-mechanism complexity
concentrates close to the ring).

This is why N has the lowest R² (0.380): the tool sees nitrogen
through only 3 blurred dimensions.  The other mechanisms are there
but the kernels can't separate them.

### Oxygen: 12 dimensions, near-field complex

Second most complex after H.  Dispersion drives the near-field
(8 dimensions at 0-4 Å).  The paramagnetic dominance (para/dia
variance ratio 1.20) means oxygen responds to multiple perturbation
channels with different angular patterns.

## What this means

The geometric kernels are not a "3-dimensional" tool.  They see:
- H: richly, in 20 angular dimensions
- C: moderately, in 6 dimensions (EFG-dominated)
- N: poorly, in 3 blurred dimensions
- O: moderately, in 12 dimensions (dispersion + EFG)

The element-dependent dimensionality IS the physics.  Each element's
electronic structure determines which geometric perturbations it
responds to and how many independent angular patterns it can
resolve.  The tool sees some things clearly and some things in the
corner of its eye.

## The "3 dimensions" claim is wrong

The raw kernel space has 1 dominant direction (not 3).  The
normalized space has element-specific dimensionality (3 to 20).
The correspondence between PCA dimensions and physics groups is
unverified and probably element-dependent.

The thesis should present the element-specific dimensionality as
the finding, not flatten it to "3 dimensions."

## Distance dependence (stable across elements)

Far-field (>8 Å): universally 3 predictive dimensions.  This IS
the multipolar regime where all kernels converge to their
asymptotic angular patterns.

Near-field (<4 Å): element-dependent complexity.  H stays at 3
(ring current dominates at all distances).  O goes to 8 (multiple
short-range mechanisms compete).  This is where the functional
forms of the kernels (wire-loop vs surface integral vs quadrupole)
diverge from each other and from the true electron distribution.
