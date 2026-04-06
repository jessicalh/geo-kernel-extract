# The Argument

## What the system does

Eight classical calculators compute geometric kernels — the spatial
shape of each shielding mechanism (ring current, bond anisotropy,
electric field gradient, etc.) at every atom in a protein. Each kernel
is a rank-2 tensor that transforms correctly under rotation. The
kernels are physics; the parameters that scale them are learned.

## What DFT provides

Quantum-mechanical DFT calculations on 720 wild-type/alanine mutant
pairs provide a reference tensor at every atom. The WT-ALA difference
isolates the contribution of each aromatic ring. This is a meter
provided by nature — it measures what the ring actually does to
shielding, in full tensor detail.

## What training does

Training against the DFT meter learns two things:

1. **Global parameters** — the ring current intensity, bond
   susceptibility anisotropy, Buckingham coupling coefficients.
   These replace literature values from 1979 with values derived
   from modern DFT on 720 proteins.

2. **Environment-dependent corrections** — the optimal ring current
   contribution at an atom depends not just on distance to the ring
   but on what other fields are present. The Coulomb field screens
   the ring effect. Nearby bonds compete for the T2 angular
   structure. An equivariant model learns these interactions from
   the local environment description that the calculators already
   provide.

## What flows back

The trained model does not leave the system. It feeds back into the
C++ extractor as a ParameterCorrectionResult — a conformation result
like any other, with declared dependencies on the calculator results
it corrects. The next extraction pass produces better features because
the calculators now account for environment.

This is not a black box replacing physics. It is a calibrated
correction to physics calculators, bounded by the geometric kernels
they compute, measured against quantum mechanics.

## What goes up

Everything. All eight raw kernels, the corrected kernels, the trained
residual, the environmental descriptors. The upper NMR prediction
model receives physically motivated, equivariant tensor features with
known provenance. It does not need to rediscover ring currents from
scratch — it gets a calibrated ring current calculator.

## The result

A more accurate Biot-Savart shielding tensor.

By feeding back 720 mutant pairs against DFT, tuning filters
against geometry, and learning which environmental interactions
matter, we arrive at a ring current tensor that better describes
nature than the one we started with.

The improvement is small. Ring currents are a small effect. But
it is a measured improvement to a physical quantity — the angular
structure of magnetic shielding from aromatic ring currents in
proteins — derived from data, not assumed from literature.

This is analytical bioscience: better measurement of a natural
phenomenon, using the tools available.


## Validation: every constant justified by quantum mechanics

The classical calculators contain constants: distance cutoffs,
near-field filter extents, switching function boundaries. Each
is a claim about where a physical approximation is valid. The
720 DFT mutation pairs provide a concrete signal against which
every one of these constants can be tested.

### Empirical range validation (no re-extraction needed)

The WT-ALA delta isolates the contribution of one aromatic ring.
Each matched atom has a known position (ρ, z) in ring-local
cylindrical coordinates and a known delta DFT T2. Plotting
delta |T2| vs distance from the removed ring traces the
empirical falloff curve of the ring current effect.

This curve answers directly:
- At what distance does the ring current fall below DFT noise?
  That is the empirical cutoff, not a literature assumption.
- Does the falloff follow 1/r³ (dipolar) or 1/r⁵ (quadrupolar)?
  That validates which multipole order dominates.
- Is the falloff isotropic? Plotting on the (ρ, z) half-plane
  shows whether atoms above the ring and in the ring plane see
  different effective ranges. If the contours are not circular,
  the cutoff should be shaped, not spherical.

This analysis requires no new computation. The data is already
being extracted. It comes free from the mutation delta analysis.

### Filter parameter sweeps (re-extraction, feasible)

The DipolarNearFieldFilter rejects atoms inside the source
distribution where the multipole expansion diverges. Its extent
parameter determines where "inside" ends. The current value
(0.5 × source extent) is physically motivated but not optimised.

Sweeping this parameter across 10 values (0.3× to 1.5× source
extent), re-extracting features at each value, and training the
same model produces a residual-vs-filter-parameter curve. The
minimum of this curve is the empirically optimal filter boundary,
validated against quantum mechanics on 720 proteins.

Cost: 10 extraction passes × 1448 conformations × 3 seconds
≈ 12 hours of compute. 100 model trainings × 2 minutes each
≈ 3 hours of GPU time. Feasible on available hardware (DGX
Spark + 4× RTX 5090).

The same approach applies to:
- Ring calculator distance cutoff (currently 15 Å)
- McConnell bond cutoff (currently 10 Å)
- Dispersion switching function boundaries (4.3–5.0 Å)
- Sequential exclusion minimum separation (currently 2 residues)

Each sweep produces a curve. Each curve has a minimum. Each
minimum is a constant justified not by literature convention
but by 720 quantum-mechanical reference calculations.

### Shaped cutoffs from the (ρ, z) plane

The most interesting possibility: the data may show that
optimal cutoffs are not simple distance thresholds. An atom
5 Å directly above a ring plane is deep in the ring current
shielding cone. An atom 5 Å in the ring plane is in the
deshielding region where the effect changes sign. These are
not the same physics, and they should not have the same cutoff.

If the delta T2 contours on the (ρ, z) half-plane are
non-circular, a shaped filter — one that accounts for position
relative to the ring, not just distance from it — would be
more physically accurate. The equivariant model (Level B/C)
learns this shape implicitly through position-dependent weights.
The validation makes it explicit and interpretable.

### What this means for the thesis

Every constant in the calculator pipeline becomes empirically
validated against quantum mechanics. The committee does not
hear "we used the cutoff from the literature" or "the AI found
values that worked." They hear: "we swept each parameter
across its physically plausible range, trained against 720 DFT
reference calculations, and determined the optimal value from
the residual minimum. Here are the curves."

This is not parameter fitting. It is systematic validation of
physical approximations against a quantum-mechanical reference,
using the mutation delta as a controlled experiment and the
T2 angular residual as the diagnostic.
