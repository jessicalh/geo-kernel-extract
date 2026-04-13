# Polarisability Roadmap

2026-04-13.  Status: design, pre-implementation.

---

## Problem

Charge polarisation is a real physical dimension in the shielding
prediction.  Stage 1 dimension inventory shows it contributes +0.197
R-squared for carbon (the dominant element-specific gap between
geometry-only and charge-aware models).  The charge redistribution
from removing an aromatic ring changes the EFG angular pattern in
ways geometry alone cannot predict.

Three existing charge sources and their limitations:

| Source     | Speed (per frame) | Signal quality | Problem |
|------------|--------------------|----------------|---------|
| MOPAC PM7  | ~45 s (889 atoms)  | Strong (correct projection for mutation delta) | Too slow for per-frame trajectory extraction |
| AIMNet2    | 0.17 s (479 atoms) | Weak (orthogonal to correct projection, cos 0.34 with Coulomb EFG) | Wrong projection, but faint signal may improve for ensemble (non-mutation) case |
| ff14SB     | free (from topology)| Static (no geometry response) | Not polarisable |

Goal: obtain affordable geometry-responsive charges at MD timescale
(thousands of frames, 685 proteins) that capture the charge-
polarisation dimension accessible to MOPAC.

---

## Approaches (all additive -- no existing code modified)

### 1. Surface normal from SASA (SasaResult extension)

**What:** Compute per-atom outward surface normal from the exposed
Fibonacci sphere test points already evaluated in SasaResult::Compute.
The average direction of non-occluded test points IS the local surface
normal.  Store as `Vec3 surface_normal` on ConformationAtom, write
`sasa_normal.npy` (N, 3).

**Why:** The current HydrationShellResult uses atom-to-protein-COM as
a proxy for surface normal.  This is wrong for non-spherical proteins
(clefts, concave surfaces, multi-domain).  Every downstream half-shell
asymmetry and dipole alignment value is computed against this crude
axis.  A proper normal is needed before the water polarisation features
can be trusted.

**Effort:** Small.  The test points are already computed and looped
over.  Accumulate their average direction, normalise, store.

**Existing code impact:** SasaResult gains a new output field.  SASA
values unchanged.  No existing calculator modified.

### 2. HydrationGeometryResult (new calculator)

**What:** New ConformationResult.  Depends on SasaResult.  Reads the
SASA-derived surface normal and first-shell water positions from
SolventEnvironment.  Writes `water_polarization.npy` (N, 10):

| Column | Quantity | Physics |
|--------|----------|---------|
| 0-2 | Net water dipole vector (Vec3) | Direction + coherence of local electrostatic asymmetry |
| 3-5 | Surface normal (Vec3, from SASA) | Reference frame for interpreting water orientation |
| 6 | Half-shell asymmetry (SASA-normal) | Exposed vs buried -- dielectric environment type |
| 7 | Dipole alignment (SASA-normal) | Radial ordering against proper surface |
| 8 | Dipole coherence (|sum d_i| / n) | How ordered vs random the first shell is |
| 9 | First-shell count | Statistical weight -- 2 waters is noise, 12 is signal |

**Why:** The existing HydrationShellResult stores only scalar
projections (mean cos) against the crude COM normal.  The net water
dipole VECTOR is T1 (L=1), a new irrep channel the calibration has
never seen.  The dipole coherence tells you whether the local
electrostatic environment is ordered (ice-like, bridging water) or
disordered (bulk-exposed).  These are direct proxies for the local
dielectric response.

**Design rule:** Existing HydrationShellResult stays untouched.  The
COM-based values continue to exist under the same names.  The new
calculator writes distinct NPY files.

**Effort:** Moderate.  New .h/.cpp, new SDK catalog entry, new
`WaterPolarizationGroup` wrapper.

### 3. Water-embedded AIMNet2

**What:** Include first-shell water atoms (O + H) in the AIMNet2 input
tensor alongside the protein atoms.  The coordinate tensor, element
tensor, and neighbour matrix grow by ~500-2000 atoms (depending on
protein size and hydration).  The neural network sees the water
electrostatic environment.  The resulting Hirshfeld charges on the
protein atoms become solvent-aware.

Report only protein-atom charges (discard water-atom outputs).  Write
as new NPY files alongside existing AIMNet2 output -- add, don't
modify.

**Why:** AIMNet2's protein-only charges are orthogonal to the correct
charge-polarisation projection (cos 0.34 with Coulomb EFG).  This may
be because the neural network, seeing only protein atoms, has no
information about the solvent environment that drives polarisation.
Including the water atoms gives the network the same information that
QM/MM embedding gives MOPAC.  GROMACS 2026's NNPot interface confirms
this is the direction the field is moving.

**Cost:** Current AIMNet2 is 0.17 s for 479 protein atoms.  With
~1000-2000 additional water atoms, estimated ~0.3-0.5 s/frame.
Affordable for per-frame trajectory extraction.

**Existing code impact:** AIMNet2Result untouched.  New calculator
(e.g. `AIMNet2WaterResult`) takes both ProteinConformation and
SolventEnvironment.

**Effort:** Moderate-high.  Need to:
- Build expanded coordinate/element tensors including water atoms
- Expand neighbour matrix to include protein-water cross-neighbours
- Map protein atom indices through the expanded tensor
- Validate charges against protein-only AIMNet2 on a known system
- Test on largest fleet protein first

### 4. EEQ calculator (Caldeweyher 2019)

**What:** Extended Electronegativity Equilibration.  Geometry-dependent
charges from minimising E(q) = sum chi_i q_i + 1/2 sum eta_i q_i^2 +
sum_ij q_i q_j gamma_ij(R) subject to sum q = Q_total.  One N x N
linear solve per frame.  Milliseconds.

Reference: Caldeweyher, Ehlert, Hansen, Neugebauer, Spicher,
Bannwarth & Grimme, J. Chem. Phys. 150, 154122 (2019).
DOI: 10.1063/1.5090222.  Used inside GFN-FF and the D4 dispersion
model.

**Why:** EEQ provides geometry-responsive charges without any external
binary or neural network.  Pure C++, no CUDA dependency.  The charges
respond to conformational change through the Coulomb kernel gamma(R),
which encodes distance-dependent screening.  This is a fundamentally
different approximation to polarisation than AIMNet2 (physics-based
vs learned), making it a genuine independent comparison rather than a
weaker version of the same thing.

**Effort:** Moderate.  The algorithm is a constrained linear solve
(Lagrange multiplier for charge neutrality).  Needs element-dependent
chi (electronegativity), eta (hardness), and r_cov (covalent radius)
parameters from the D4 parameter set.  Pure C++ implementation,
no external dependencies beyond Eigen (already in project).

**Existing code impact:** None.  New `EeqResult` ConformationResult.

### 5. E-field temporal variance (already computed)

**What:** Per-atom variance of water E-field magnitude across
trajectory frames.  Atoms in ordered hydration shells (structural
water, ice-like packing) have low variance.  Bulk-exposed atoms
have high variance.

**Why:** Direct proxy for local dielectric response.  An atom whose
water E-field fluctuates wildly is in a high-dielectric environment
(bulk solvent).  An atom with stable E-field is in an ordered
environment (buried, or structural water bridge).  This variance
IS the dielectric constant, sampled empirically.

**Status:** Already computed.  `GromacsProteinAtom.water_emag` is a
Welford accumulator updated every frame.  The variance is available
after the trajectory scan as `water_emag.Variance()`.  Currently
written to `atom_catalog.csv`.

**Remaining work:** Route the per-atom E-field variance into the
`water_polarization.npy` block at extraction time so the SDK consumer
gets it alongside the other polarisation proxy features.  This may
mean the column count grows from 10 to 11, or the variance is written
as a separate scalar array.  Design decision deferred to
implementation.

---

## What we considered and rejected

### Drude-FF MD

Run the 685-protein fleet with CHARMM Drude-2019 polarisable force
field.  Polarisation would be free in every frame (the Drude particle
displacement x charge IS the induced dipole).  2-5x MD cost.

**Rejected:** GROMACS Drude integrator has no GPU support.  685
proteins x 20 ns on CPU would take months on 4 machines.  Not
feasible for the fleet timeline (starting 2026-04-15).

### xTB / GFN2-xTB

Tight-binding DFT.  Geometry-dependent Mulliken charges, ~100x cheaper
than MOPAC.

**Rejected:** Prior experience: xTB segfaulted on >450 atoms with 7%
success rate on our 725-protein dataset (replaced by MOPAC, see
`learn/bones/SESSION_MOPAC_INTEGRATION.md`).  O(N^3) diagonalisation
makes 5000-atom proteins prohibitive even if stable.  Not worth
revisiting for per-frame work.

### Water-embedded MOPAC

Include first-shell water point charges in the MOPAC QM/MM
calculation.  Standard approach, SolventEnvironment already identifies
the waters.

**Rejected for per-frame:** 45 s/frame is not viable at fleet scale.
Could be used for selected-frame validation against water-embedded
AIMNet2 to verify the projection is correct.

### Perturbative charge correction

First-order correction: delta_q = alpha x div(E) from existing
Coulomb field.  Microseconds.

**Rejected:** Tried with AIMNet2 and found consistency across frames
in later MD is uncontrollable.  The perturbative approximation breaks
down for large conformational changes where the linear response
assumption fails.

### Thole induced dipole model

Apply tabulated atomic polarisabilities to the water + protein
E-field.  Iterate to self-consistency.

**Not rejected but deferred:** Conceptually attractive (the water
E-field is already computed, just multiply by alpha).  But this
produces induced dipoles, not charges, and the EFG correction from
induced dipoles requires a different code path than the existing
charge -> Coulomb -> EFG pipeline.  May revisit after water-embedded
AIMNet2 and EEQ results clarify whether a third approach is needed.

---

## Implementation order

1. **SasaResult surface normal** -- small, unblocks (2)
2. **HydrationGeometryResult** -- new calculator, writes
   `water_polarization.npy`, depends on (1)
3. **Water-embedded AIMNet2** -- independent of (1)-(2), but wants
   the surface normal for validation.  Test on largest fleet protein
   first.
4. **EEQ calculator** -- fully independent.  Pure C++.  Can proceed
   in parallel with (3).
5. **E-field variance routing** -- after (2) defines the NPY format

Items 1-2 should complete before the fleet re-run starts (Wednesday
2026-04-15) so the extraction pipeline can produce the new arrays.
Items 3-4 can proceed in parallel, no deadline dependency.

---

## Validation plan

All new charge sources must be tested against MOPAC on the 723-protein
mutation set where we have ground truth:

- Compute EEQ charges and water-embedded AIMNet2 charges on the same
  proteins
- Feed through existing Coulomb EFG pipeline
- Compare direction cosine with MOPAC Coulomb EFG (the correct
  projection)
- Report per-element R-squared contribution vs the known +0.197
  carbon gap

The validation question: does EEQ/water-AIMNet2 access the same
charge-polarisation dimension as MOPAC, or a different one?  If
different, it may be complementary (additional dimension) rather than
a substitute.

---

## SDK packaging

New SDK group: `WaterPolarizationGroup`.  One coherent block for
"everything that tells the model about charge polarisation from the
water environment."

```python
protein = nmr_extract.load(path)
pol = protein.water_polarization   # WaterPolarizationGroup

pol.dipole_vector     # (N, 3) net first-shell water dipole
pol.surface_normal    # (N, 3) SASA-derived outward normal
pol.asymmetry         # (N,)   half-shell asymmetry (SASA normal)
pol.dipole_alignment  # (N,)   cos(dipole, SASA normal)
pol.coherence         # (N,)   |sum d_i| / n
pol.shell_count       # (N,)   first-shell water count
```

The surface normal is duplicated from SASA into this group
intentionally -- the consumer should not have to reconstruct the
relationship between SASA geometry and water orientation.

Existing `HydrationGroup` and `WaterFieldGroup` unchanged.

---

## Dependencies on fleet re-run

The fleet re-run (starting 2026-04-15, 685 proteins, single-walker
20 ns metadynamics, 4 machines, 3-4 weeks) produces trajectories only.
No extraction features are baked in.  Extraction runs after.

This means:
- All new calculators can be developed and tested after the fleet run
  starts
- The surface normal and HydrationGeometryResult are extraction-time
  changes, not MD-time
- Water-embedded AIMNet2 reads from the same full-system .xtc that
  the existing AIMNet2 reads -- no MD protocol change needed
- EEQ reads protein geometry only -- no trajectory dependency at all

The only hard constraint: `compressed-x-grps = System` in the .mdp
(include water in the .xtc).  This was already fixed for the re-run.
