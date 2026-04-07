# Calculator Audit: Environment Tests and Implicit Assumptions

**Date: 2026-04-02**

This document catalogues every place where a calculator tests the
environment — distance checks, atom classification, threshold
decisions, proxy relationships — and records what the test is doing,
whether the physics reason is documented, and whether it's a proxy
for something we could obtain from structure or calculation.

The calculators were written by an agent that 95% of the time knew
what it was doing. The problem is not that the decisions are wrong —
it's that they're not explained, so we can't tell the 5% from the 95%.
This audit makes every decision visible so that multiple models
(human, Opus, Gemini) can review them before we change anything.

---

## Classification of environment tests

Each entry is classified as one of:

- **PHYSICS**: first-principles criterion with documented reason.
  Leave alone unless the physics is wrong.
- **NUMERICAL**: guards against floating-point pathology. Not physics,
  but necessary. Should be documented with what it prevents.
- **PROXY**: uses one quantity (usually distance) as a stand-in for
  another (usually topology or electronic structure). Flag for
  review — could we check the real thing?
- **VERNIER**: an arbitrary choice that affects every downstream result.
  Must be documented, should be tunable.
- **UNDOCUMENTED**: a number or test with no stated reason. Highest
  priority for review.

---

## Global (PhysicalConstants.h)

### MIN_DISTANCE = 0.1 A
- **Class**: NUMERICAL
- **What it does**: singularity guard on all 1/r^n terms
- **Physics reason**: no two non-bonded atoms are ever 0.1A apart in
  a physical protein. This prevents division by ~zero from producing
  Inf/NaN.
- **Documented**: yes (comment says "singularity cutoff")
- **Concern**: none. Correct and necessary.

### NO_DATA_SENTINEL = 99.0
- **Class**: NUMERICAL
- **What it does**: initialises "nearest distance" trackers to a value
  larger than any real distance in a protein (~40A for large proteins)
- **Physics reason**: none — it's a sentinel value
- **Documented**: yes
- **Concern**: fragile if proteins exceed 99A diameter (they don't).
  Could use std::numeric_limits<double>::max() but 99.0 is fine.

### NEAR_ZERO_NORM = 1e-10
- **Class**: NUMERICAL
- **What it does**: threshold below which a vector is considered zero
  for normalisation purposes
- **Documented**: yes
- **Concern**: none

### RING_CALC_CUTOFF = 15.0 A
- **Class**: PHYSICS
- **What it does**: spatial search radius for ring-based calculators
  (BS, HM, RingSusceptibility, PiQuadrupole, Dispersion)
- **Physics reason**: ring current 1/r^3 at 15A is ~3e-4 A^-3.
  Constitution specifies 15A as the distance where ring current
  effects become negligible.
- **Documented**: yes (Constitution reference)
- **Concern**: no switching function at 15A, but the signal at the
  boundary is tiny. Dispersion has its own inner cutoff with a taper.
  For BS/HM/RingSusceptibility/PiQuadrupole the 1/r^3 or 1/r^5
  decay makes the boundary contribution negligible.

### EXP_DECAY_LENGTH = 4.0 A
- **Class**: VERNIER
- **What it does**: characteristic length for exponential decay
  weighting in MutationDeltaResult
- **Physics reason**: "the characteristic scale" (comment)
- **Documented**: minimally. No reference or derivation.
- **Concern**: this is a model vernier. It weights the importance of
  nearby vs far atoms in mutation delta analysis. The value 4A is
  plausible (roughly the range of ring current effects at 10% of
  maximum) but not derived from anything specific.

### PACKING_RADIUS = 8.0 A
- **Class**: VERNIER
- **What it does**: radius for heavy atom count (local packing density)
- **Physics reason**: not stated
- **Documented**: comment says "for heavy atom count" but not why 8A
- **Concern**: another model vernier. 8A is roughly twice the first
  solvation shell. Plausible but undocumented.

### RING_COUNT_SHELL_1..4 = 3, 5, 8, 12 A
- **Class**: VERNIER
- **What it does**: concentric shells for ring proximity counting
  (n_rings_within_3A, etc.)
- **Physics reason**: not stated individually. The set of shells
  provides a radial profile of ring environment.
- **Documented**: no physics reason for these specific radii
- **Concern**: these become features. The choice of shell radii
  is a model vernier.

### HBOND_COUNT_RADIUS = 3.5 A
- **Class**: VERNIER
- **What it does**: radius for counting H-bonds near each atom
- **Physics reason**: 3.5A is approximately the maximum N...O
  distance for a hydrogen bond. So "H-bonds within 3.5A" means
  "H-bonds whose midpoint is within 3.5A" — which is roughly
  "H-bonds where this atom is between the donor and acceptor."
- **Documented**: no
- **Concern**: the 3.5A is reasonable but should state what it
  represents physically.

### HBOND_MAX_DIST = 50.0 A
- **Class**: NUMERICAL
- **What it does**: maximum N...O distance for a resolved H-bond
- **Physics reason**: effectively no cutoff. 1/r^3 at 50A is ~8e-6.
  The real cutoff is that DSSP only identifies H-bonds up to ~3.5A
  N...O distance. The 50A limit is just a safety net.
- **Documented**: no
- **Concern**: misleading. The actual physics cutoff is DSSP's
  Kabsch-Sander criterion (~3.5A). 50A is a no-op guard.

### SEQUENTIAL_EXCLUSION_THRESHOLD = 2
- **Class**: PHYSICS
- **What it does**: minimum sequence separation for H-bond dipolar
  calculation
- **Physics reason**: through-bond coupling dominates through-space
  dipolar at |i-j| < 2
- **Documented**: yes (in SequentialExclusionFilter description)
- **Concern**: only used in HBondResult during H-bond resolution,
  not as a SequentialExclusionFilter in the filter set. The
  exclusion is applied before the filter framework sees it.

### APBS_SANITY_LIMIT = 100.0 V/A
- **Class**: NUMERICAL
- **What it does**: clamps E-field magnitude
- **Physics reason**: APBS grid interpolation can produce artifacts
  at grid edges or near charges. 100 V/A is well beyond physical
  E-fields in proteins (Case 1995: 1-10 V/A at backbone amide H).
- **Documented**: comment says "V/Angstrom" but not why 100
- **Concern**: also used in CoulombResult (see below).

### COULOMB_KE = 14.3996 V*A
- **Class**: PHYSICS
- **What it does**: Coulomb's constant in {e, A, eV} units
- **Documented**: yes, with derivation in comment and APPLIED_MATHS_FIXES
- **Concern**: none. This is a universal constant, correctly derived.

### KT_OVER_E_298K = 0.025693 V
- **Class**: PHYSICS
- **What it does**: converts APBS kT/e units to Volts at 298.15 K
- **Documented**: yes
- **Concern**: hardcoded temperature (298.15 K). If APBS is run at
  a different temperature, this constant would be wrong. Currently
  APBS temperature is also hardcoded to 298.15 K, so they're
  consistent.

---

## McConnellResult

### MCCONNELL_CUTOFF_A = 10.0 A
- **Class**: VERNIER
- **What it does**: spatial search radius for bond midpoints
- **Physics reason**: 1/r^3 at 10A is 0.001 A^-3. Beyond this the
  bond anisotropy contribution is < 0.1% of a nearby (2A) bond.
- **Documented**: comment says "Cutoff radius for bond midpoint
  search (Angstroms)" but no physics justification for 10A vs 15A
- **Concern**: no switching function. Hard cutoff creates a
  discontinuity for atoms near the boundary on MD trajectories.
  The signal at 10A is ~10x larger than at 15A, making this a
  sharper discontinuity than the ring cutoff.

### SelfSourceFilter
- **Class**: PHYSICS
- **What it does**: excludes bond endpoints from their own bond's
  dipolar field
- **Physics reason**: field undefined at source
- **Documented**: yes, via KernelEvaluationFilter framework
- **Concern**: none

### DipolarNearFieldFilter (source_extent = bond_length)
- **Class**: PHYSICS
- **What it does**: excludes atoms within half a bond length of
  the bond midpoint
- **Physics reason**: multipole expansion invalid inside source
- **Documented**: yes
- **Concern**: bond_length ~1.2-1.5A, so threshold ~0.6-0.75A.
  Rarely triggers beyond SelfSourceFilter. Honest but low impact.

### Category classification (PeptideCO, PeptideCN, etc.)
- **Class**: PHYSICS
- **What it does**: sorts bond contributions by bond type
- **Physics reason**: different bond types have different magnetic
  anisotropies (Delta_chi). CO vs CN vs aromatic have different
  magnitudes and orientations.
- **Documented**: yes (BondCategory enum, GEOMETRIC_KERNEL_CATALOGUE)
- **Concern**: none. Classification uses typed BondCategory, not strings.

### Nearest CO/CN tracking
- **Class**: PHYSICS (dominant contribution)
- **What it does**: tracks the nearest PeptideCO and PeptideCN bonds
  and their kernel values
- **Physics reason**: the nearest peptide CO bond is typically the
  dominant McConnell contribution (backbone amide H is ~2A from
  its own peptide CO). This is a key feature.
- **Documented**: implicitly
- **Concern**: none

### Traceless projection on category T2 totals (line 252-261)
- **Class**: NUMERICAL
- **What it does**: extracts symmetric part of accumulated M,
  removes trace, then decomposes for T2 features
- **Physics reason**: the symmetric traceless part IS the T2 content.
  Floating-point accumulation can break tracelessness.
- **Documented**: yes (inline comment)
- **Concern**: none. Correct.

---

## CoulombResult

### Self-exclusion: inline `if (i == j) continue`
- **Class**: PHYSICS
- **What it does**: excludes atom from its own Coulomb field
- **Physics reason**: same as SelfSourceFilter
- **Documented**: no (no filter framework, no logging)
- **Concern**: should use KernelFilterSet for consistency. Currently
  the only calculator that bypasses the framework.

### Charge threshold: `if (std::abs(q_j) < 1e-15) continue`
- **Class**: NUMERICAL
- **What it does**: skips atoms with zero or negligible charge
- **Documented**: no
- **Concern**: minor optimisation, physically harmless. The 1e-15
  threshold is effectively zero (charges are typically 0.01-1.0 e).

### Backbone classification via Residue backbone cache
- **Class**: PHYSICS
- **What it does**: classifies source atoms as backbone (N, CA, C, O,
  H, HA, CB) or sidechain using the residue's pre-built backbone
  index cache
- **Physics reason**: backbone and sidechain produce different angular
  patterns of electric field. Decomposition allows the model to weight
  them differently.
- **Documented**: yes (inline comment says "from residue backbone cache")
- **Concern**: CB is classified as backbone here. This is a structural
  choice — CB is technically the first sidechain carbon, but it's
  part of the backbone framework. This is a **PROXY**: the real
  question is "which atoms' charges contribute to the backbone
  electrostatic environment?" CB's charge influences the backbone
  E-field at amide H, so including it with backbone is defensible.
  But it should be documented as a choice.

### Aromatic classification via Ring atom membership
- **Class**: PHYSICS
- **What it does**: classifies atoms that belong to any ring as
  "aromatic" for E-field decomposition
- **Documented**: yes
- **Concern**: ring membership is topological (from Protein), not a
  distance proxy. Correct.

### Primary bond direction for E_bond_proj
- **Class**: PROXY
- **What it does**: for each atom, finds a "primary bond direction"
  to project E-field along. H atoms use the parent-H bond. Heavy
  atoms use the first bond in their bond list.
- **Physics reason**: Buckingham's E_z requires the E-field projected
  along the bond axis. For amide H, this is well-defined (N-H bond).
  For heavy atoms, "first bond" is arbitrary.
- **Documented**: inline comment, but the arbitrariness of heavy-atom
  choice is not discussed
- **Concern**: **this is a PROXY for the local chemical environment**.
  The Buckingham effect physically depends on E projected along the
  local symmetry axis of the atom's electron cloud. For sp2 carbons
  this might be the p-orbital axis (perpendicular to the plane), not
  any single bond. For sp3, there's no single preferred direction.
  The current implementation uses "first bond" which is effectively
  random for heavy atoms. The model can learn to ignore this for
  non-H atoms, but it's a known approximation.

### E-field clamping with APBS_SANITY_LIMIT
- **Class**: UNDOCUMENTED (in Coulomb context)
- **What it does**: clamps E_total magnitude to 100 V/A, scales all
  decomposition components uniformly
- **Physics reason**: originally for APBS grid artifacts. In Coulomb,
  large E-fields at short distances are physical (not artifacts).
- **Documented**: reuses APBS constant without separate justification
- **Concern**: **PROXY**. Clamping physical Coulomb E-fields treats
  extreme values as artifacts when they're actually real. At
  MIN_DISTANCE = 0.1A, E ≈ ke*q/r^2 ≈ 14.4*1/0.01 = 1440 V/A,
  well above the 100 V/A clamp. The SelfSourceFilter prevents i==j,
  but very close non-bonded pairs could still produce large fields.
  The clamp affects only the E-field, not the EFG (inconsistent).
  This should either be: (a) justified as a Winsorisation for
  numerical stability with documented reasoning, or (b) replaced
  by a switching function that tapers E-field contributions at
  short range.

### Solvent contribution: APBS - vacuum
- **Class**: PHYSICS
- **What it does**: E_solvent = E_apbs - E_coulomb_total
- **Physics reason**: the solvated (APBS) field includes screening
  by the protein dielectric and solvent. Subtracting vacuum Coulomb
  isolates the solvation effect.
- **Documented**: yes
- **Concern**: only computed if ApbsFieldResult is present. Correct
  use of HasResult<T>() pattern.

### backbone_frac: projection of E_backbone along E_total direction
- **Class**: PHYSICS (derived feature)
- **What it does**: measures how much of the total E-field comes from
  backbone charges (signed projection, not ratio)
- **Physics reason**: a signed projection is stable near cancellation
  (unlike |E_bb|/|E_total| which diverges when total is small)
- **Documented**: yes (inline comment explains stability choice)
- **Concern**: none. Good numerical thinking.

---

## RingSusceptibilityResult

### DipolarNearFieldFilter (source_extent = ring diameter)
- **Class**: PHYSICS
- **What it does**: excludes atoms inside the ring diameter
- **Documented**: yes, with comment about model validity
- **Concern**: ring atoms themselves are at radius ~1.4A from center,
  ring diameter ~2.8A, so threshold = 1.4A. Ring atoms are AT the
  threshold boundary. Some may pass, some may not, depending on
  exact geometry. This is a **boundary case** worth documenting.

### No self-exclusion for ring atoms
- **Class**: missing (potential issue)
- **What it does**: does NOT explicitly exclude ring atoms from their
  own ring's susceptibility field
- **Physics reason**: ring atoms are inside the ring's charge
  distribution. The dipolar model from the ring center is invalid
  for the ring's own atoms.
- **Documented**: no
- **Concern**: DipolarNearFieldFilter probably catches these (ring
  atoms are at radius ≈ ring radius, threshold = ring radius), but
  it's geometry-dependent. For non-planar rings (TRP pyrrole) some
  atoms could be further from center than others. A topological check
  (is this atom a vertex of this ring?) would be unambiguous. Compare
  to Dispersion's BondedToVertices.

### RingNeighbourhood creation and reuse
- **Class**: NUMERICAL (data structure)
- **What it does**: searches existing ring_neighbours by ring_index
  to avoid duplicates. If another calculator (BS) already created the
  RingNeighbourhood, reuses it and adds chi fields.
- **Documented**: inline comment
- **Concern**: linear search over ring_neighbours. Not a performance
  issue (typically < 10 neighbours) but the pattern is repeated in
  every ring calculator.

---

## HBondResult

### DSSP as H-bond source
- **Class**: PHYSICS
- **What it does**: uses DSSP Kabsch-Sander criterion to identify
  backbone H-bonds
- **Physics reason**: DSSP's energy criterion identifies geometrically
  valid H-bonds from backbone N-H...O=C pairs
- **Documented**: yes
- **Concern**: backbone only. Sidechain H-bonds (Ser-OH, Tyr-OH,
  Lys-NH3, etc.) are not included. This is a **known limitation**:
  DSSP does not identify sidechain H-bonds. For a complete H-bond
  calculator, we'd need a separate sidechain H-bond identification
  step using distance + angle criteria from the bond graph.

### H-bond resolution: donor N → acceptor O
- **Class**: PHYSICS
- **What it does**: resolves DSSP residue-level H-bond to atom
  positions using backbone index cache (res.N, res.O)
- **Documented**: yes (inline comments)
- **Concern**: correct use of typed topology, no string matching

### Deduplication by (donor_N, acceptor_O) pair
- **Class**: NUMERICAL
- **What it does**: prevents double-counting when DSSP reports the
  same H-bond from both residues' perspectives
- **Documented**: yes
- **Concern**: none. Correct.

### Sequential exclusion: `seq_sep < SEQUENTIAL_EXCLUSION_THRESHOLD`
- **Class**: PHYSICS
- **What it does**: excludes H-bonds between adjacent residues
  (|i-j| < 2)
- **Physics reason**: through-bond coupling dominates at small
  sequence separation. The dipolar model doesn't capture this.
- **Documented**: yes (references SequentialExclusionFilter)
- **Concern**: applied during H-bond resolution, not via the
  KernelFilterSet. The SequentialExclusionFilter exists but isn't
  used here — the exclusion happens at H-bond identification time
  (line 155), not at kernel evaluation time. Functionally equivalent,
  but the filter framework doesn't see it, so rejection counts
  don't include sequential exclusions.

### SelfSourceFilter + DipolarNearFieldFilter
- **Class**: PHYSICS
- **What it does**: excludes donor N and acceptor O from their own
  H-bond's field; excludes atoms inside the N...O source distribution
- **Documented**: yes
- **Concern**: source_extent = hb.distance (N...O distance, ~2.8A).
  Threshold = ~1.4A. The donor H atom is ~1.0A from the N, so
  ~0.4A from the midpoint — well inside the threshold. The acceptor
  C is bonded to O and can be ~1.2A from the midpoint. Both are
  correctly excluded. Good.

### Nearest H-bond tracking
- **Class**: PHYSICS (feature)
- **What it does**: tracks nearest H-bond midpoint to each atom
- **Concern**: recomputes the kernel for the nearest H-bond after
  the main loop (line 309). Minor inefficiency but correct.

### `hbond_is_donor` / `hbond_is_acceptor` detection (line 319-322)
- **Class**: PHYSICS
- **What it does**: checks whether this atom IS a donor N or acceptor O
  in any H-bond
- **Concern**: O(N × H) scan. Could be pre-built into a set. Not a
  correctness issue.

---

## BiotSavartResult

### Wire segment numerical guards (WireSegmentField)
- `lenA < 1e-25` and `lenB < 1e-25`: **NUMERICAL**. Field point at
  segment endpoint. In SI metres, so 1e-25 m ≈ 1e-15 A. Effectively
  zero distance.
- `crossSq < 1e-70`: **NUMERICAL**. Field point on the wire axis
  (cross product is zero). 1e-70 m^2 ≈ 1e-50 A^2.
- **Documented**: yes (inline comments from old code)
- **Concern**: none. These are deep numerical guards in SI that
  correspond to physically impossible situations.

### Johnson-Bovey: unit current I = 1.0 nA
- **Class**: PHYSICS
- **What it does**: computes B-field with unit current to produce a
  geometric kernel independent of intensity
- **Documented**: yes (inline comment line 182)
- **Concern**: none. Correct separation of kernel from parameter.

### G_ab = -n_b * B_a * PPM_FACTOR sign convention
- **Class**: PHYSICS
- **What it does**: the minus sign ensures sigma = I * G has the
  correct physical sign
- **Documented**: yes, extensively. Verified analytically.
- **Concern**: none.

### DipolarNearFieldFilter (source_extent = ring diameter)
- **Class**: PHYSICS
- **Same concern as RingSusceptibilityResult**: ring atoms at the
  boundary.

### Per-type accumulation: `if (ti >= 0 && ti < 8)`
- **Class**: NUMERICAL
- **What it does**: bounds check before storing per-type T0/T2 sums
- **Concern**: magic constant 8. Should reference enum size.

### Ring proximity counting (n_rings_within_3A, 5A, 8A, 12A)
- **Class**: VERNIER (feature design)
- **What it does**: counts rings in concentric shells
- **Lives in**: BiotSavartResult only
- **Concern**: general-purpose feature computed in a specific
  calculator. If BS isn't run, these are zero.

### `total_B_field += B_total` accumulation pattern
- **Class**: NUMERICAL
- **What it does**: uses += not =
- **Concern**: singleton guarantee prevents double-run. The +=
  would allow another calculator to contribute, but none does.
  Not a bug but inconsistent with other calculators' = pattern.

---

## HaighMallionResult

### Adaptive subdivision thresholds
- `HM_SUBDIV_THRESHOLD_L1 = 2.0 A`: **NUMERICAL**
- `HM_SUBDIV_THRESHOLD_L2 = 1.0 A`: **NUMERICAL**
- **What they do**: when any triangle vertex is within 2.0A (level 1)
  or 1.0A (level 2) of the field point, subdivide the triangle for
  better quadrature accuracy
- **Physics reason**: the 1/r^3 integrand varies rapidly when the
  field point is close to the surface. Standard quadrature misses
  this variation.
- **Documented**: yes (inline comments, APPLIED_MATHS_FIXES)
- **Concern**: the thresholds 2.0A and 1.0A come from the old project.
  They are numerical accuracy choices, not physics. The old code used
  these values and they were validated empirically. Could be verified
  by convergence testing (halve the thresholds, check if results change).

### Triangle area check: `if (triArea < 1e-20) return`
- **Class**: NUMERICAL
- **What it does**: skips degenerate triangles (zero area)
- **Concern**: none. Correct guard.

### Quadrature point distance check: `if (rhoMag < MIN_DISTANCE)`
- **Class**: NUMERICAL
- **What it does**: skips quadrature points where the field point is
  at the surface point (singularity)
- **Concern**: none. Uses global MIN_DISTANCE.

### Boyd-Skrynnikov construction: G = -n ⊗ (H·n)
- **Class**: PHYSICS
- **What it does**: constructs the full rank-1 shielding kernel from
  the raw surface integral H and the ring normal n
- **Documented**: yes (reference: Boyd & Skrynnikov JACS 2002)
- **Concern**: none. Verified analytically and numerically.

---

## PiQuadrupoleResult

### Stone T-tensor EFG formula
- **Class**: PHYSICS
- **What it does**: computes the EFG from a point axial quadrupole at
  the ring center
- **Documented**: yes (Stone Ch. 3 reference, derivation in catalogue,
  corrected from original)
- **Concern**: none. Tracelessness verified at machine precision.

### Scalar: (3 cos^2 theta - 1) / r^4
- **Class**: PHYSICS
- **What it does**: Buckingham A-term scalar for the quadrupole potential
- **Documented**: yes
- **Concern**: none

### DipolarNearFieldFilter (same as other ring calcs)
- **Class**: PHYSICS
- **Same concern**: ring atom boundary case

---

## DispersionResult

### DISP_VERTEX_R_CUT = 5.0 A
- **Class**: VERNIER (well-documented)
- **What it does**: maximum vertex-atom distance for dispersion contact
- **Physics reason**: C6/r^6 at 5A is 0.03% of value at 2A. Truncation
  error < 0.1% of total sum.
- **Documented**: yes, extensively. Physics reason stated with numbers.

### DISP_VERTEX_R_SWITCH = 4.3 A
- **Class**: VERNIER (well-documented)
- **What it does**: onset of smooth taper
- **Physics reason**: 0.7A taper width is comparable to MD fluctuations
  (~0.5A RMS), ensuring smooth features across frames
- **Documented**: yes. Reference: Brooks et al. 1983.

### CHARMM switching function
- **Class**: PHYSICS
- **What it does**: C^1 continuous taper from full strength to zero
- **Documented**: yes, with full formula and reference
- **Concern**: none. This is the gold standard for smoothness.

### BondedToVertices exclusion
- **Class**: PHYSICS (topology check)
- **What it does**: excludes atoms bonded to any ring vertex from
  dispersion evaluation for that ring
- **Physics reason**: through-space 1/r^6 kernel does not model
  through-bond electronic coupling
- **Documented**: yes
- **Concern**: none. This is the RIGHT way to do it — topology,
  not distance proxy. Other ring calculators should be reviewed
  for whether they need this too.

### DipolarNearFieldFilter (ring level)
- **Class**: PHYSICS
- **Same concern**: boundary case for ring atoms. But Dispersion's
  BondedToVertices already catches ring atoms and their bonded
  neighbours, so the filter is a belt-and-suspenders here.

---

## ApbsFieldResult

### VdW radius fallback (lines 128-136)
- **Class**: UNDOCUMENTED
- **What it does**: if vdw_radius < 0.5, substitutes element-specific
  defaults
- **Values**: C=1.70, N=1.63, O=1.48, S=1.78, H=1.00, default=1.50
- **Documented**: comment says "approximate ff14SB typical values"
  but no reference
- **Concern**: these define the molecular surface for the PB solve.
  They directly affect the APBS E-field at every atom. The trigger
  condition (r < 0.5) is a **PROXY** for "radius not set by
  ChargeAssignmentResult." Since ChargeAssignment is a declared
  dependency, this should either never trigger (in which case
  document that and add a diagnostic log if it does) or be
  documented as a known fallback with referenced values.

### Grid sizing: extent + 40A padding, coarse +30A
- **Class**: VERNIER
- **What it does**: determines APBS grid dimensions
- **Physics reason**: 20A padding per side for boundary condition
  accuracy. Coarse grid extends further for the focusing method.
- **Documented**: comment says "matches APBS mg-auto convention"
- **Concern**: standard APBS practice. No issue.

### Grid dimension: 161^3
- **Class**: VERNIER
- **What it does**: determines grid resolution (~0.3-0.5A spacing)
- **Documented**: comment says "targeting ~0.3-0.5A spacing"
- **Concern**: 161 is the APBS recommended odd dimension. Standard.

### Protein dielectric: pdie = 4.0
- **Class**: VERNIER
- **What it does**: sets the protein interior dielectric constant
- **Physics reason**: commonly used value. Published range: 2-20.
  4.0 is typical for implicit solvent PB calculations.
- **Documented**: no reference
- **Concern**: **model vernier**. Affects every APBS field. Should
  cite a reference. Consider making tunable.

### Solvent dielectric: sdie = 78.54
- **Class**: PHYSICS
- **What it does**: sets the bulk solvent dielectric
- **Physics reason**: dielectric constant of water at 25C
- **Documented**: comment says "(water, 25C)"
- **Concern**: standard value. Correct.

### Temperature: 298.15 K
- **Class**: PHYSICS
- **What it does**: Debye-Huckel ionic screening parameter
- **Documented**: yes
- **Concern**: coupled with KT_OVER_E_298K. Must stay consistent.

### Ionic strength: 0.15 M
- **Class**: VERNIER
- **What it does**: physiological ionic strength
- **Physics reason**: 150 mM NaCl is standard physiological conditions
- **Documented**: comment says "(physiological)"
- **Concern**: standard value. For non-physiological conditions
  (e.g., NMR buffer) this would need adjustment.

### E-field from grid: central differences
- **Class**: NUMERICAL
- **What it does**: E = -grad(phi) via central differences on the
  potential grid, EFG via central differences of E
- **Documented**: yes
- **Concern**: this means E uses 2 grid points (±h), EFG uses 4
  grid points (two E evaluations, each using ±h). The effective
  resolution for EFG is 4× the grid spacing. At 0.5A spacing,
  EFG resolution is 2A. This limits the accuracy of APBS EFG
  for atoms close together. This is a known limitation of grid-based
  methods and is correct.

### Traceless projection on EFG
- **Class**: PHYSICS + NUMERICAL
- **What it does**: removes the self-potential artifact from the grid
  EFG. The APBS potential includes each atom's own Coulomb field,
  whose Laplacian is a delta function at the atom site. The grid
  discretises this as a large finite trace.
- **Documented**: yes (excellent inline comment, lines 93-101)
- **Concern**: none. Correct and well-explained.

### NaN/Inf sanitisation
- **Class**: NUMERICAL
- **What it does**: zeros out pathological values from grid interpolation
  at atom positions near grid edges or at grid singularities
- **Documented**: yes
- **Concern**: none

### E-field clamping (APBS_SANITY_LIMIT)
- **Class**: NUMERICAL (in APBS context)
- **What it does**: clamps extreme E-fields from grid artifacts
- **Physics reason**: in APBS context, large fields indicate grid
  artifacts (atom too close to grid edge, insufficient resolution).
  This is correct for APBS.
- **Documented**: partially
- **Concern**: when reused in CoulombResult, the context changes
  (see CoulombResult section above).

### Unit conversion: kT/(e*A) to V/A
- **Class**: PHYSICS
- **What it does**: converts APBS native units to V/A for comparison
  with CoulombResult
- **Documented**: yes (inline comment + PhysicalConstants.h)
- **Concern**: none. Correct.

---

## Cross-calculator concerns

### Ring atom self-exclusion (BS, HM, RingSusceptibility, PiQuadrupole)

Four ring calculators use DipolarNearFieldFilter with source_extent =
ring diameter. Ring atoms are at distance ≈ ring radius from center,
and the filter threshold is source_extent/2 = ring radius. So ring
atoms are AT the threshold boundary.

Dispersion uses BondedToVertices (topology check) which unambiguously
excludes ring atoms and their immediate bonded neighbours.

**Question for review**: should BS/HM/RingSusceptibility/PiQuadrupole
also use a topological exclusion for ring atoms? The DipolarNearFieldFilter
is probably sufficient (ring atoms are inside the source distribution
by definition), but the boundary case means some ring atoms might
pass or fail depending on ring planarity.

### McConnell vs Dispersion: cutoff treatment

McConnell uses a hard 10A cutoff with no taper. Dispersion uses a
5A cutoff with a CHARMM switching function (onset 4.3A). Both model
through-space effects that decay with distance.

**Question for review**: should McConnell get a switching function?
At 10A the signal is ~0.001 A^-3, which is 10x the signal at
Dispersion's 5A cutoff. On an MD trajectory with 0.5A fluctuations,
an atom at 9.8A on one frame and 10.2A on the next gets a
discontinuous feature.

### HBOND_MAX_DIST = 50A vs actual DSSP range

The 50A is a no-op — DSSP N...O distances are always < 3.5A. The
HBondResult loop evaluates every identified H-bond against every atom
with no spatial cutoff on the evaluation side. The 1/r^3 decay handles
range naturally, but for very large proteins with many H-bonds, a
spatial index query would be more efficient.

### Filter framework usage across calculators

| Calculator         | KernelFilterSet? | Filters used                        |
|--------------------|------------------|-------------------------------------|
| McConnellResult    | YES              | SelfSource + DipolarNearField       |
| CoulombResult      | NO               | inline i!=j only                    |
| RingSusceptibility | YES              | DipolarNearField                    |
| HBondResult        | YES              | SelfSource + DipolarNearField       |
| BiotSavartResult   | YES              | DipolarNearField                    |
| HaighMallionResult | YES              | DipolarNearField                    |
| PiQuadrupoleResult | YES              | DipolarNearField                    |
| DispersionResult   | YES              | DipolarNearField + BondedToVertices |

CoulombResult is the only outlier.

### Magic constant 8 (ring type count)

Used in: BS, HM, PQ, Disp per-type accumulation.
Should reference enum size, not magic number.

---

## Summary: what needs multi-model review before changes

1. **McConnell switching function**: should we add one? What onset/cutoff?
2. **CoulombResult filter framework**: straightforward consistency fix
3. **CoulombResult E-field clamping**: is APBS_SANITY_LIMIT appropriate
   for vacuum Coulomb? Should it be a separate constant?
4. **Ring atom topological exclusion**: should BS/HM/RingSusceptibility/PQ
   use BondedToVertices like Dispersion?
5. **ApbsFieldResult VdW fallback**: does it ever trigger? If so, are
   the values correct?
6. **ApbsFieldResult pdie = 4.0**: should this be tunable? Reference?
7. **Primary bond direction in CoulombResult**: the heavy-atom "first bond"
   choice is arbitrary. Document or improve?
8. **EXP_DECAY_LENGTH, PACKING_RADIUS, RING_COUNT_SHELL radii**: model
   verniers without documented derivation

Each of these should be reviewed by at least two models before
any code change. The 95% that's correct should not be damaged by
fixing the 5% that's undocumented.
