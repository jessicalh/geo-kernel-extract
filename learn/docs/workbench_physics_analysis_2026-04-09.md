# Physics Analysis of Workbench Results — 2026-04-09

Agent analysis of the R analytical workbench run (99 proteins, 127 kernels).

## Why MopacEFG_aro captures 75% of the ridge signal

The dominant L=2 perturbation from deleting an aromatic ring is
**electrostatic** (charge redistribution), not **magnetic** (ring current).

The EFG is a direct electrostatic tensor — it depends only on the charge
distribution and geometry, not on a model of the ring current. The ring
current mechanisms (BS, HM) compute the magnetic effect of delocalized pi
electrons via a model current loop. They are correct in form but their
amplitude depends on the ring current intensity parameter. The EFG,
computed from PM7 charges already calibrated against QM, does not need
a ring current intensity.

All ring current kernels (BS, HM for PHE/TYR/TRP variants) point in
nearly the same L=2 direction at each atom — they all encode "there is
a ring over there with this orientation." They differ in amplitude and
distance dependence but not in angular structure. The EFG encodes the
charge quadrupole of the ring, a different angular pattern.

The ~3 effective L=2 directions are:
1. EFG from aromatic charges (dominant)
2. Ring current dipolar (BS/HM/ring susceptibility, all collinear)
3. Residual from backbone EFG or bond anisotropy

## Why φ-modulated kernels fail but φ scalars work

The φ-modulated **kernels** pre-compute cos(φ) × BS_T2 — a fixed L=2
tensor weighted by a single global coefficient. The **scalar** cos_φ,
as an interaction with the top-10 kernel base, creates cos(φ) × w_k ×
kernel_k for every kernel k — the MLP can learn a **different**
φ-dependence for each kernel in each chemical environment.

ring_proximity (+0.305) nearly doubles the variance explained because
it converts 10 global coefficients into 10 × 61 = 610 spatially
modulated coefficients. The ridge can now say "weight BS more when
above the ring, weight EFG more when in the plane."

## Why per-ring (72 kernels) barely beats ring-type (32 kernels)

Per-ring kernels are computed as **residuals** from the per-type mean.
For single-ring-per-type proteins (common case), the residual is zero.
The 0.022 marginal gain comes from multi-ring proteins and per-ring
distance ordering. Summing over rings of the same type loses almost
nothing in the linear regime.

## Why bond_cat and total groups are useless

The WT-ALA delta zeroes most bond anisotropy:
- Backbone doesn't change (mechanical mutant)
- Sidechain changes only at the mutated residue (short-range)
- Coulomb total cancels — what remains is the charge difference,
  which EFG_aro already captures directly
- HBond barely changes (ring atoms rarely in H-bonds)

## Five hypotheses to test analytically

### A: Residual concentrates at the magic-angle region

Compute per-atom residuals from full ridge. Bin by (ρ, z) relative to
nearest ring. If residual concentrates near ρ~2-4Å, z~2-4Å (magic
angle cone), the kernel functional forms break down there.

**Test:** Add radial basis function kernels: exp(-(r-r₀)²/σ²) × existing
kernel, for a grid of r₀ values. Forward selection identifies if any
radial-localized kernel enters early.

### B: Per-protein variation is ring current intensity variation

Compute per-protein ridge coefficients for top-3 kernels. Plot against
protein-level descriptors (n_aromatics, size, dipole). If EFG coefficient
varies smoothly with dipole magnitude, a single scalar interaction
(protein_dipole × EFG_weight) could recover much of the 0.38→0.81 gap.

### C: HIE's problem is charge confounding, not weak ring current

Run forward selection on HIE-only atoms separately. If APBS_EFG or
DeltaAPBS_EFG ranks higher for HIE than for PHE, the point-charge EFG
is inadequate for the partially charged imidazole and the solvation
correction matters specifically for HIE.

### D: Ring-ring interactions are missing

Add per-ring distance-to-nearest-other-ring as a scalar. The shielding
perturbation from removing ring A depends on whether ring B is nearby
(screening/reinforcement). Ring-ring distances are computable from
ring_geometry.npy. Test as scalar interaction.

### E: Decompose ring_proximity to find the dominant coordinate

**THE FIRST TEST TO RUN.** Instead of all 61 ring_proximity scalars,
test subsets:
- (a) just 1/dist per ring (6 scalars)
- (b) just z + rho per ring (12 scalars)
- (c) just cos_φ + sin_φ per ring (12 scalars)
- (d) just disp_scalar + disp_contacts per ring (12 scalars)

The answer tells us whether the path forward is:
- Radial basis functions (if 1/dist dominates)
- Angular basis functions (if z/rho dominate)
- Azimuthal structure (if cos_φ/sin_φ dominate)

This costs nothing — subset analysis of existing export data.
