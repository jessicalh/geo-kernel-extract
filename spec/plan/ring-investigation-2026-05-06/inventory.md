# Ring-investigation inventory — 2026-05-06

This document is the inventory phase output of the Bundle C ring
investigation. Its job is to enumerate every site in the runtime
library and the substrate generator that touches ring data, classify
each site, and feed dispatch decisions for the per-calculator phase.
Per the framing memo
(`spec/plan/ring-investigation-2026-05-06/README.md`), Bundle C plans
must be empirical: the per-site claims here cite specific files and
line ranges, not speculation.

The investigation runs read-only against repository HEAD on
2026-05-06; no code is modified.


## §1 — Coverage methodology

### Grep commands run

The prior-agent's coverage was reproduced verbatim. These three
commands form the load-bearing audit surface:

```bash
grep -rni "ring" src/ --include='*.h' --include='*.cpp' | wc -l
# Expected: 2908

grep -rln -i "ring" src/ --include='*.h' --include='*.cpp' | wc -l
# Expected: 160 files

grep -rnE "RingTypeIndex|RingNeighbourhood|RingSystemKind|RingPositionLabel|RingMembership|RingPosition\b" \
    src/ --include='*.h' --include='*.cpp' | wc -l
# Expected: 608 typed-grep matches

grep -rlE "RingTypeIndex|RingNeighbourhood|RingSystemKind|RingPositionLabel|RingMembership|RingPosition\b" \
    src/ --include='*.h' --include='*.cpp' | wc -l
# Expected: 22 files

grep -rni "ring" tools/ | wc -l
# Expected: 425 (420 in build_semantic_tables.cpp + ancillary)

grep -ni "ring" src/generated/LegacyAmberSemanticTables.cpp | wc -l
# Expected: 507
```

### Match counts

| Surface                                        | Matches | Files |
| ---------------------------------------------- | ------- | ----- |
| `src/` case-insensitive `ring`                  | 2908    | 160   |
| `src/` typed-grep (the 6 typed enum/struct names) | 608     | 22    |
| `tools/` case-insensitive `ring`                | 425     | 4     |
| `src/generated/LegacyAmberSemanticTables.cpp`   | 507     | 1     |

Verified independently 2026-05-06 by this agent. The user can re-run
any of the commands above to re-confirm.

### Generated-vs-runtime breakdown

`src/generated/LegacyAmberSemanticTables.cpp` is the substrate-table
generator's emitted output. Its 507 ring matches are typed enum
literals (`RingSystemKind::Benzene_Phe`, `RingPositionLabel::Ipso`,
etc.) inside the per-residue
`AtomSemanticTable` constexpr arrays, plus their headers and
comments. That file is consumed at runtime via `LookupBy` /
`LookupCap` / `ApplyCapDelta` (`src/LegacyAmberTopology.cpp`); no
calculator reads it directly.

Subtracting the generator-emitted file leaves 2908 - 507 = 2401
matches in 159 hand-written runtime files that touch ring data via
includes, struct fields, calculator code, comments, and tests.

The `tools/topology/build_semantic_tables.cpp` file (420 of the
425 tool-tree ring matches) is the substrate generator: it reads
CCD via cifpp, runs RDKit, and emits the typed-enum table. Per the
header comment in `src/SemanticEnums.h:14-23`, the runtime library
does NOT link RDKit; the string barrier is enforced at link time.
The generator is out of scope for Bundle C.


## §2 — Substantive ring-current calculator entries

This section enumerates the substantive consumers of Ring data:
calculators whose physics computation depends on Ring chemistry.
Each entry answers Q1-Q8 of the per-report template at inventory
depth (the per-calculator phase deepens these).

### 2.1 BiotSavartResult

**File and line range.**
`src/BiotSavartResult.cpp:1-503`. The Compute method spans lines
124-345; SampleBFieldAt/SampleShieldingAt span 356-418;
WriteFeatures spans 439-500.

**Physics, one sentence.** BiotSavart computes the magnetic field at
a probe atom from a circulating ring current via the Johnson-Bovey
two-loop line integral with current loops at ±jb_lobe_offset from
the ring plane (Johnson & Bovey 1958; Case 1995, J. Biomol. NMR 6,
341-346); shielding contribution σ = I × G with sign convention
σ_ab = -dB_a^sec/dB_{0,b} preserved in
`G_ab = -n_b · B_a · PPM_FACTOR` (verified analytically per
PATTERNS.md §22: Phe at literature I = -12.0 nA/T, probe 3 Å axial
above ring face, G_T0 = -0.116 → σ = (-12)(-0.116) = +1.40 ppm
shielded; sign convention catches early development bug not caught
by compilation or unit tests).

**Ring API surface consumed.**
- `Protein::RingCount()` (line 134), `Protein::RingAt(ri)` (line
  171, line 365, line 395) — base Ring access by index.
- `Ring::JBLobeOffset()` (lines 181, 186, 223, 380, 406) — typed
  virtual; reads `CalculatorConfig::Get` keyed per ring type.
- `Ring::Intensity()` (line 185) — recorded as a GeometryChoice
  AddNumber, NOT used in the kernel computation. The kernel
  evaluates with unit current `I = 1.0 nA` (line 223); intensity
  multiplication is calibration-time downstream.
- `Ring::TypeIndexAsInt()` (line 302) — for indexing into
  `ConformationAtom::per_type_G_T0_sum` (8-array) and
  `per_type_G_T2_sum` (8×5-array).
- `Ring::type_index` (line 248) — RingTypeIndex stored on each new
  RingNeighbourhood entry.
- `Ring::atom_indices` — read indirectly via the
  `RingBondedExclusionFilter` constructor at
  `src/KernelEvaluationFilter.cpp:22` (the filter walks
  `Protein::RingAt(ri).atom_indices` to build per-ring
  exclusion sets).
- `RingGeometry::vertices`, `::center`, `::normal`, `::radius` —
  read via `conf.ring_geometries[ri]` (line 172). The geometry is
  populated by `GeometryResult::Compute` at
  `src/GeometryResult.cpp:17-19` calling `Ring::ComputeGeometry`.
- `ConformationAtom::ring_neighbours` (lines 239-244, 276-277,
  320-325) — per-atom storage; this calculator owns
  `RingNeighbourhood::G_tensor`, `G_spherical`, `B_field`,
  `B_cylindrical`, `cos_phi`, `sin_phi`, plus the cylindrical
  coordinates `rho`, `z`, `theta`. (HM, PiQuad, RingChi, Disp
  populate the same struct's other tensor fields; they reuse the
  same RingNeighbourhood by `ring_index` matching.)
- `ConformationAtom::n_rings_within_3A`, `_5A`, `_8A`, `_12A` —
  proximity counts (lines 320-325).
- Spatial query: `SpatialIndexResult::RingsWithinRadius` (line
  165) with `ring_current_spatial_cutoff = 15.0 Å`.

**String discipline check.** None. No reads of
`atom.pdb_atom_name`, no `aatype.rings[].atom_names` access, no
`residue.three_letter` substring operations. Verified by inspection
of all 503 lines.

**`RingTypeIndex::Count` interaction.** **HARDCODED `8`.**
- Line 303: `if (ti >= 0 && ti < 8)` — guard before per-type
  accumulation.
- Line 455: `std::vector<double> data(N * 8)` — `bs_per_type_T0`
  output sized as `(N, 8)`.
- Line 457: `for (int t = 0; t < 8; ++t)` — loop bound.
- Line 459: `NpyWriter::WriteFloat64(... N, 8)` — explicit shape.
- Line 465: `std::vector<double> data(N * 40)` —
  `bs_per_type_T2` output sized as `(N, 40)` (8 ring types × 5
  T2 components).
- Line 467: `for (int t = 0; t < 8; ++t)`.
- Line 469: `data[i*40 + t*5 + c] = ...`.
- Line 470: `NpyWriter::WriteFloat64(... N, 40)`.

The `RingTypeIndex::Count = 8` constant is defined at
`src/Types.h:184`, but BiotSavart uses the literal `8` rather than
the enum constant. **This contradicts the prior agent's claim that
"no hardcoded `[8]` literals were found";** the literals are
present and they live on `ConformationAtom` as well
(`per_type_G_T0_sum: std::array<double, 8>`,
`per_type_G_T2_sum: std::array<std::array<double, 5>, 8>` at
`src/ConformationAtom.h:131-132`). Bundle C extending RingTypeIndex
to 9 (adding `ProPyrrolidine`) makes these sites Bundle C migration
targets.

**Recommended classification:** SUBSTANTIVE.

**Per-calculator report needed?** Yes.
Filename: `BiotSavartResult.md`.

**Planned-calculator implication.** BiotSavartResult is the
foundation for the planned `B4. Distributed ring current`
(NMR_EXTRACT_DESIDERATA §B; cross-referenced in
`spec/plan/planned-calculator-substrate-audit-2026-05-06.md` lines
264-272), which decomposes ring-current contributions per ring
vertex (per-`RingPositionLabel.position` typing). Today's per-ring
totals pool across vertices. Bundle C does not change BS; but
`ConstructRingsFromSubstrate`'s preservation of `atom_indices`
ordering (matching the substrate's `RingPositionLabel`-indexed
walk) is the gate for a future per-position-stratified BS variant.
B4 also exploits the `RingMembership.n_heteroatoms` field for
HID/HIE asymmetric per-vertex models; that field is reachable today
via `LegacyAmber().SemanticAt(ai).ring_position.primary.n_heteroatoms`
and does not require new substrate work.

**Pro pyrrolidine impact at inventory depth.** A
ProPyrrolidineRing object with `Intensity() = 0`, `Aromaticity() =
None`, `RingSizeValue() = 5`, `NitrogenCount() = 1` and
`atom_indices` walking N → Cα → Cβ → Cγ → Cδ enters the BS
iteration at line 170. The Johnson-Bovey kernel at line 221
evaluates with unit current (intensity is downstream), so a Pro
Ring produces a non-zero raw kernel `G` for nearby probe atoms
(within 15 Å). `Ring::TypeIndexAsInt() = 8` (assuming Pro adopts
index 8) at line 302: the current literal-`8` guard
`if (ti >= 0 && ti < 8)` would EXCLUDE the Pro contribution from
the per-type accumulators on `ConformationAtom`. The total
accumulator `G_total` (line 298) and `bs_shielding_contribution`
(line 317) WOULD include the Pro raw kernel, but with downstream
calibration the total comes from per-type×intensity sums; if
intensity is forced to zero in calibration for Pro, the
contribution vanishes. The verbatim `Builder.Record(...)` lambda
on lines 182-196 emits `AddRing(... Source/Included)`,
`AddNumber("intensity", 0.0, "nA")`, `AddNumber("lobe_offset", ...,
"A")` for Pro rings — documentary value.


### 2.2 HaighMallionResult

**File and line range.**
`src/HaighMallionResult.cpp:1-432`. Compute spans 190-357;
SampleShieldingAt spans 365-393; WriteFeatures spans 407-429.

**Physics, one sentence.** HaighMallion computes the dipolar surface
integral `H_ab = ∫ (3 ρ_a ρ_b / ρ⁵ - δ_ab / ρ³) dS` over a
triangulated ring face using 7-point Gaussian quadrature
(Stroud T2:5-1 / Dunavant degree-5 per `src/HaighMallionResult.cpp:42-43`)
with adaptive subdivision (level 1 below 2.0 Å vertex distance,
level 2 below 1.0 Å, max depth 2 → 7 → 28 → 112 quadrature points
per fan triangle); the effective field is `V = H · n` and the
shielding kernel is `G_ab = -n_b · V_a` (rank-1, same sign
convention as BS); BS-vs-HM cosine similarity is 0.999 (PATTERNS.md
§22 — parallel physics, different mathematical approximation).

**Ring API surface consumed.**
- `Protein::RingCount()` (line 200), `Protein::RingAt(ri)`
  (lines 235, 365 implicit via geom).
- `RingGeometry::vertices`, `::center`, `::normal`, `::radius` —
  via `conf.ring_geometries[ri]` (line 236).
- `Ring::TypeIndexAsInt()` (line 336) for the per-type 8-array
  accumulators.
- `Ring::type_index` (line 307).
- `RingBondedExclusionFilter` (line 216) — same filter as BS.
- Spatial query: `SpatialIndexResult::RingsWithinRadius` (line
  230) with cutoff 15.0 Å.
- Per-atom storage: `ConformationAtom::ring_neighbours[*]`'s HM
  fields `hm_H_tensor`, `hm_H_spherical`, `hm_B_field`,
  `hm_G_tensor`, `hm_G_spherical` (lines 326-330);
  `hm_shielding_contribution` (line 347);
  `per_type_hm_T0_sum[8]`, `per_type_hm_T2_sum[8][5]` (lines
  338-340).

**String discipline check.** None. Verified by inspection of all
432 lines.

**`RingTypeIndex::Count` interaction.** **HARDCODED `8`.**
- Line 337: `if (ti >= 0 && ti < 8)` guard.
- Line 413: `std::vector<double> per_type_T0(N * 8)`.
- Line 414: `std::vector<double> per_type_T2(N * 40)`.
- Line 419: `for (int t = 0; t < 8; ++t)`.
- Line 426: `NpyWriter::WriteFloat64(... N, 8)` for
  `hm_per_type_T0`.
- Line 427: `NpyWriter::WriteFloat64(... N, 40)` for
  `hm_per_type_T2`.

Same hardcoding pattern as BiotSavart.

**Recommended classification:** SUBSTANTIVE.

**Per-calculator report needed?** Yes.
Filename: `HaighMallionResult.md`.

**Planned-calculator implication.** HM shares the planned-calculator
surface with BS for the ring-current physics path
(`spec/plan/planned-calculator-substrate-audit-2026-05-06.md` §B4
"Distributed ring current" lines 264-272 — same per-vertex
decomposition is meaningful for both); `B14-B15` NICS probe
evaluation (planned-audit lines 185-196) consumes ring centroids,
which `RingGeometry::center` (Bundle C preserves) supplies. The
BS-vs-HM cosine 0.999 (PATTERNS.md §22) is what
`MoynaRingCurrentComparisonResult` (planned, planned-audit lines
434-435) will validate per-method per-ring. HM at probe atoms
inside the source distribution was the source of the 13% HM excess
artifact corrected 2026-04-02 (per README.md §6); the
`RingBondedExclusionFilter` topological exclusion is what makes
this physically sane.

**Pro pyrrolidine impact at inventory depth.** Pro Ring entering HM
at line 234: `SurfaceIntegral` (line 282) computes a non-zero raw
H from the Pro ring atoms' triangulation; downstream G is non-zero;
`G_total` accumulates non-zero. As with BS, the per-type 8-guard
EXCLUDES Pro from per-type accumulators if Pro = index 8;
`hm_shielding_contribution` (line 347) includes it. With
intensity-zero in calibration Pro vanishes from final shielding;
the raw G value is the documentary record of the iteration. The
`Builder.Record` lambda at lines 241-246 emits "surface integral"
for every ring including Pro.


### 2.3 RingSusceptibilityResult

**File and line range.**
`src/RingSusceptibilityResult.cpp:1-270`. Compute spans 103-228;
SampleShieldingAt spans 231-250; WriteFeatures spans 259-267.

**Physics, one sentence.** RingSusceptibility computes the
point-magnetic-dipole-at-ring-center kernel
`M_ab/r³ = [9 cos θ d̂_a n_b - 3 n_a n_b - (3 d̂_a d̂_b - δ_ab)] / r³`
(per RingSusceptibilityResult.cpp:84-91), the McConnell
derivation specialised to a magnetic ring with `b̂ → n̂` (ring
normal); the kernel decomposes into T0 (from `-3 n̂⊗n̂`),
T1 (from `9 cosθ d̂⊗n̂`), T2 (from `-(3 d̂⊗d̂ - I)`), and the
scalar `f = (3cos²θ - 1)/r³` is the Pople ring susceptibility
factor (Pople 1956); reference for the McConnell derivation:
McConnell 1957, J. Chem. Phys. 27, 226-228.

**Ring API surface consumed.**
- `Protein::RingCount()` (line 113), `Protein::RingAt(ri)` (line
  151).
- `RingGeometry::vertices`, `::center`, `::normal`, `::radius` —
  via `conf.ring_geometries[ri]` (line 152).
- `Ring::type_index` (line 189) — for RingNeighbourhood typing.
- `RingBondedExclusionFilter` (line 132).
- Spatial query: `SpatialIndexResult::RingsWithinRadius` (line
  146).
- `ConformationAtom::ring_neighbours[*]`'s
  `chi_tensor`, `chi_spherical`, `chi_scalar` fields (lines
  208-210); `ringchi_shielding_contribution` (line 218).

**String discipline check.** None.

**`RingTypeIndex::Count` interaction.** None directly. This
calculator does NOT emit per-type T0/T2 NPY arrays. It outputs only
the total `ringchi_shielding.npy` shape `(N, 9)` (line 265). No
hardcoded `8` array sizes.

**Recommended classification:** SUBSTANTIVE.

**Per-calculator report needed?** Yes.
Filename: `RingSusceptibilityResult.md`.

**Planned-calculator implication.** Per the substrate audit
(§A.2 lines 174-181), `BulkSusceptibilityAccumulator` (planned A3)
sums per-atom McConnell + RingSusceptibility + HBond contributions —
RingSusceptibility's per-atom M tensor IS one of the three sums.
`A4. NICSProbeEvaluator` (planned-audit lines 183-196) likewise
consumes ring susceptibility at probe positions; the sampling
machinery (`SampleShieldingAt`) in this calculator is the inroad.
The fact that RingSusceptibility doesn't emit per-type NPY today
makes it the cleanest target if the per-type stratification is
later judged unnecessary; conversely it's the noisiest if later
work wants HID/HIE asymmetry per ring type.

**Pro pyrrolidine impact at inventory depth.** Pro Ring entering at
line 150 with Aromaticity = None: `ComputeRingChiKernel` at line
54-96 produces a non-zero `M_over_r3` (the kernel is purely
geometric, not intensity-weighted in this calculator either).
`M_total` at line 213 accumulates non-zero. `chi_scalar` field on
RingNeighbourhood (line 210) records the Pople f factor for the Pro
geometry. There is no per-type 8-guard here, so Pro contribution is
NOT excluded; it lands in `ringchi_shielding_contribution` directly.
Downstream calibration must handle Pro's intensity-zero physics
either by upweighting RingTypeIndex::ProPyrrolidine to zero in the
fit, or by gating in the calibration pre-processing.


### 2.4 PiQuadrupoleResult

**File and line range.**
`src/PiQuadrupoleResult.cpp:1-292`. Compute spans 109-236;
SampleShieldingAt spans 239-258; WriteFeatures spans 267-289.

**Physics, one sentence.** PiQuadrupole computes the EFG geometric
kernel from a point axial quadrupole at the ring center, derived
via Stone's T-tensor formalism `V_ab = -(Θ/2) T_abcd n_c n_d`
(Stone 2013, "The Theory of Intermolecular Forces" Ch. 3) yielding
`G_ab = 105 dn² d_a d_b / r⁹ - 30 dn (n_a d_b + n_b d_a)/r⁷ -
15 d_a d_b / r⁷ + 6 n_a n_b / r⁵ + δ_ab (3/r⁵ - 15 dn²/r⁷)`
(symmetric, traceless by Laplace; verified numerically); the
scalar Buckingham A-term is `(3 cos² θ - 1)/r⁴` (Buckingham 1959).

**Ring API surface consumed.**
- `Protein::RingCount()` (line 119), `Protein::RingAt(ri)`
  (line 156).
- `RingGeometry::vertices`, `::center`, `::normal` — via
  `conf.ring_geometries[ri]` (line 157).
- `Ring::type_index` (line 192).
- `Ring::TypeIndexAsInt()` (line 213).
- `RingBondedExclusionFilter` (line 138).
- Spatial query: `SpatialIndexResult::RingsWithinRadius` (line
  151).
- `ConformationAtom::ring_neighbours[*]`'s `quad_tensor`,
  `quad_spherical`, `quad_scalar` (lines 208-210);
  `piquad_shielding_contribution` (line 226);
  `per_type_pq_scalar_sum[8]`, `per_type_pq_T2_sum[8][5]`
  (lines 215-217).

**String discipline check.** None.

**`RingTypeIndex::Count` interaction.** **HARDCODED `8`.**
- Line 214: `if (ti >= 0 && ti < 8)`.
- Line 272: `std::vector<double> per_type_T0(N * 8)`.
- Line 273: `std::vector<double> per_type_T2(N * 40)`.
- Line 278: `for (int t = 0; t < 8; ++t)`.
- Line 286: `NpyWriter::WriteFloat64(... N, 8)`.
- Line 287: `NpyWriter::WriteFloat64(... N, 40)`.

Same pattern. Same Bundle C migration target.

**Recommended classification:** SUBSTANTIVE.

**Per-calculator report needed?** Yes.
Filename: `PiQuadrupoleResult.md`.

**Planned-calculator implication.** The 1/r⁹ leading term means
PiQuadrupole is the most near-field-sensitive ring kernel —
exactly the regime where Pro pyrrolidine atoms are most relevant
(Pro ring atoms 5 Å from a probe contribute ~10× more than at 8 Å
versus 1/r³ kernels' ~3.4× change). Per the planned-calculator
audit, no specific planned calculator names PiQuadrupole, but the
audit's §B "Variations on existing calculators" includes "B2.
Smooth cutoffs on K_ab kernels" (lines 258-259) and "B3.
Multipole-expanded Coulomb" (lines 261-262); PiQuadrupole is the
existing example of a higher-multipole expansion.

**Pro pyrrolidine impact at inventory depth.** Pro Ring entering
at line 155: kernel is non-zero (geometric, not intensity-weighted).
The 1/r⁹ scaling makes Pro contributions decay rapidly but they're
non-zero for short-range probes (typical Pro Ring – probe Cα at
~3 Å gives `f ≈ 0.33/r⁴`). The per-type 8-guard at line 214
EXCLUDES Pro from `per_type_pq_*_sum`; `piquad_shielding_contribution`
includes it. Filter outcome: `RingBondedExclusionFilter` rejects
all Pro residue's own atoms; `DipolarNearFieldFilter` rejects when
distance < ring_diameter/2 (≈ 1.2 Å for a 2.4 Å pyrrolidine
diameter). The `Builder.Record` lambda at lines 170-178 fires for
filter exclusions only; passing Pro evaluations produce no
GeometryChoice record (the calculator's recording surface is
narrower than BS's).


### 2.5 DispersionResult

**File and line range.**
`src/DispersionResult.cpp:1-423`. Compute spans 163-360;
SampleShieldingAt spans 363-389; WriteFeatures spans 398-420.

**Physics, one sentence.** Dispersion computes the per-vertex London
dispersion kernel
`K_ab = S(r) (3 d_a d_b / r⁸ - δ_ab / r⁶)` with `r = |r_atom -
r_vertex|` and `S(r)` the CHARMM smooth switching function
(Brooks et al. 1983, J. Comput. Chem. 4, 187 — `S(r) = 1` for
`r ≤ R_switch = 4.3 Å`; cubic taper to zero at `R_cut = 5.0 Å`,
C¹ continuous at both boundaries); per-vertex summation with the
through-bond ring-bonded exclusion (a topological check using
`Protein::AtomAt(vi).bond_indices`, NOT a distance proxy);
traceless per vertex (`Tr(K) = S(r)·(3|d|²/r⁸ - 3/r⁶) = 0`).

**Ring API surface consumed.**
- `Protein::RingCount()` (line 173).
- `Protein::RingAt(ri)` (lines 204, 219, 369).
- `Ring::atom_indices` (lines 139, 275) — read directly to compute
  vertex sums and to build per-ring bonded exclusion sets.
- `Ring::type_index` (line 313).
- `Ring::TypeIndexAsInt()` (line 336).
- `RingGeometry::vertices`, `::center`, `::normal`, `::radius` —
  via `conf.ring_geometries[ri]` (line 220).
- Spatial query: `SpatialIndexResult::RingsWithinRadius` (line
  214).
- Filter set: `DipolarNearFieldFilter` only (line 189) — note: NO
  `RingBondedExclusionFilter`. Dispersion implements its own
  through-bond exclusion via the `BondedToVertices` helper at
  lines 136-149 plus the per-ring `ring_bonded` set at line 247.
  This is the prior agent's "RingBondedExclusionFilter applied"
  claim corrected: Dispersion has equivalent topology but
  reimplemented locally rather than using the shared filter
  class.
- `ConformationAtom::ring_neighbours[*]`'s
  `disp_tensor`, `disp_spherical`, `disp_scalar`, `disp_contacts`
  (lines 330-333); `disp_shielding_contribution` (line 348);
  `per_type_disp_scalar_sum[8]`, `per_type_disp_T2_sum[8][5]`
  (lines 338-340).

**String discipline check.** None.

**`RingTypeIndex::Count` interaction.** **HARDCODED `8`.**
- Line 337: `if (ti >= 0 && ti < 8)`.
- Line 403: `std::vector<double> per_type_T0(N * 8)`.
- Line 404: `std::vector<double> per_type_T2(N * 40)`.
- Line 409: `for (int t = 0; t < 8; ++t)`.
- Line 417: `NpyWriter::WriteFloat64(... N, 8)`.
- Line 418: `NpyWriter::WriteFloat64(... N, 40)`.

Same pattern.

**Recommended classification:** SUBSTANTIVE. Distinct from the other
ring-current calculators in two ways: (a) it iterates per ring
vertex inside the per-ring loop (line 275:
`for (size_t vi = 0; vi < ring.atom_indices.size(); ++vi)`),
making it the only calculator that reads `ring.atom_indices` for
geometric computation rather than only for filter set
construction; (b) it reimplements ring-bonded exclusion at
file scope (lines 136-149) rather than using
`RingBondedExclusionFilter`. PATTERNS.md anti-pattern: distinct
implementations of the same topological invariant.

**Per-calculator report needed?** Yes.
Filename: `DispersionResult.md`.

**Planned-calculator implication.** Per the substrate audit, "B7.
C8 dispersion term" (lines 283-284) is a planned extension of the
1/r⁶ to 1/r⁸ leading term; that future calculator would inherit
the per-vertex summation pattern and hence the hardcoded `[8]`
arrays + the BondedToVertices helper. Pro pyrrolidine geometry
matters more for Dispersion than for the ring-current calculators
because dispersion is intensity-INDEPENDENT — `Intensity = 0` does
not zero out the dispersion contribution. Pro Ring atoms ARE real
dispersion sources at ~5 Å, with nonzero K kernels weighted by S(r).

**Pro pyrrolidine impact at inventory depth.** Pro Ring entering at
line 218 with non-zero K per vertex. The per-vertex inner loop at
line 275 iterates the 5 Pro ring atoms (N, Cα, Cβ, Cγ, Cδ);
`ComputeDispVertex` at line 101 produces non-zero K for each
vertex within the 5.0 Å cutoff. **This is the calculator where Pro
contributes most non-trivially** — non-zero T0 + T2 from a
saturated heterocycle is real chemistry (Pro's electron cloud
participates in dispersion regardless of aromaticity). The per-type
8-guard at line 337 EXCLUDES Pro from `per_type_disp_*_sum`;
`disp_shielding_contribution` (line 348) includes it. Through-bond
exclusion at line 247 (`if (ring_bonded[ri].count(ai))`) skips
the entire ring iteration when the probe atom is itself bonded to
any Pro vertex — this protects the Pro residue's own backbone N,
Cα, etc. from spurious self-dispersion. The `Builder.Record`
lambdas at lines 233, 249, 262, 284 emit "near-field exclusion",
"through-bond exclusion", "dispersion taper", "switching noise
floor" respectively for Pro rings; non-trivial documentary record.


### 2.6 Trajectory-scope BS aggregators (4 results, bundled)

The four trajectory-scope ConformationResult subclasses that consume
BiotSavart per-frame output are NOT direct ring consumers; they read
`ConformationAtom::bs_shielding_contribution` (T0 / T2 magnitude)
that BiotSavart wrote each frame.

**Files.**
- `src/BsWelfordTrajectoryResult.{h,cpp}` (running mean / variance
  / min / max of T0, |T2|, T0 deltas; AV-pattern exemplar per
  PATTERNS.md §14).
- `src/BsAnomalousAtomMarkerTrajectoryResult.{h,cpp}` (z-score
  outlier detection vs BsWelford's running distribution; Phase 4
  cross-Result-read exemplar).
- `src/BsShieldingTimeSeriesTrajectoryResult.{h,cpp}` (per-atom
  per-frame SphericalTensor time series; canonical e3nn-irrep
  emission shape for SphericalTensor TRs).
- `src/BsT0AutocorrelationTrajectoryResult.{h,cpp}` (per-atom
  biased ACF on T0 across N_LAGS = 120 frames).

**Physics, one sentence each.**
- `BsWelfordTrajectoryResult` accumulates Welford's online
  variance estimator (Welford 1962) on T0, |T2|, and T0 frame-to-
  frame deltas across the trajectory's frames.
- `BsAnomalousAtomMarkerTrajectoryResult` flags per-atom per-frame
  events when `|z| > 2.0` after `MIN_BURN_IN_FRAMES = 20`,
  emitting two event kinds (`BsAnomalyHighT0`,
  `BsAnomalyLowT0`) into the per-atom event bag (Pattern C per
  PATTERNS.md).
- `BsShieldingTimeSeriesTrajectoryResult` archives the full
  SphericalTensor time series for downstream ACF / spectral
  density / S² extraction; pins parity `0e+1o+2e` and irrep layout
  `T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2`
  in H5 attributes.
- `BsT0AutocorrelationTrajectoryResult` computes the biased ACF
  estimator `ρ(k) = C(k)/C(0)` with
  `C(k) = (1/N) Σ_t (x_t - ⟨x⟩)(x_{t+k} - ⟨x⟩)`,
  guaranteeing `|ρ(k)| ≤ 1` (per the file header at lines 17-23).

**Ring API surface consumed.** None directly. Each reads
`conf.AtomAt(i).bs_shielding_contribution` (a SphericalTensor on
`ConformationAtom`); none calls `Protein::RingAt` or `Ring::*`
virtuals; none reads `ring_neighbours`. Verified by grep:
`grep -nE "ring|Ring" src/BsWelford*.cpp src/BsShielding*.cpp
src/BsT0*.cpp src/BsAnomalous*.cpp` returns only string matches in
log messages and comments, no executable ring code.

**String discipline check.** None.

**`RingTypeIndex::Count` interaction.** None.

**`AminoAcidType::rings` interaction.** None.

**Recommended classification:** SUBSTANTIVE (downstream); but only
because BS's T0/|T2| values they consume are themselves ring-
chemistry-derived. They do not need ring-API surgery; if BS is
correct, they are correct.

**Per-calculator report needed?** Yes, **bundled**.
Filename: `Trajectory-BS-Aggregators.md` (covers all four).

**Planned-calculator implication.** Per the planned-calculator
audit §1.3 (TS1-TS8) and §A.5 (A7-A9, A10), the trajectory-scope
ACF / spectral density / Mori-Zwanzig kernel / Lipari-Szabo / S²
results all have the same shape: read per-frame output of an
underlying ConformationResult, accumulate, finalise. The four BS
aggregators are the worked examples for this pattern; each
planned addition (TS1 BlockAveraged convergence, TS2 Lipari-Szabo,
TS5 CCR rates, TS7 memory kernel) clones the shape. None of those
would be affected by Bundle C — they consume scalar/tensor fields
on ConformationAtom whose semantics are unchanged by adding a Pro
ring.

**Pro pyrrolidine impact at inventory depth.** Indirect: if BS
produces a non-zero `bs_shielding_contribution` from a Pro Ring
on some atoms (because BS's per-RingTypeIndex 8-guard excludes the
Pro contribution from `per_type_*` accumulators but the total
`bs_shielding_contribution` includes it), the trajectory aggregators
ingest the (slightly-Pro-perturbed) values and propagate them
downstream. With Pro `Intensity = 0` from calibration, the actual
shielding contribution is zero and the aggregators see no change.
With raw geometric kernels (pre-calibration), they see the Pro's
geometric perturbation. This is the cleanest expression of the
pre-calibration / post-calibration boundary — kernels are what the
project outputs, shielding is what calibration produces from
kernels.


## §3 — Infrastructure entries

This section enumerates infrastructure: the runtime substrate that
ring-current calculators sit atop. Three named infrastructure
surfaces; one of them is the migration target Bundle C replaces.

### 3.1 DetectAromaticRings (the migration target)

**File and line range.**
`src/Protein.cpp:548-646` — the function being replaced by Bundle
C's `ConstructRingsFromSubstrate`.

**Purpose.** Walks `Protein::residues_`, for each residue with
`AminoAcidType::is_aromatic = true` reads `aatype.rings[]` (from
the static table at `src/AminoAcidType.cpp`), matches each
ring_def's `atom_names` const-char-pointer array against a per-
residue `pdb_atom_name → atom_index` map, creates the Ring
subclass via `CreateRing(effective_type)`, and pushes into
`rings_`. Sets `fused_partner_index` for TRP after both
sub-rings land.

**Ring API surface produced.** Builds `Protein::rings_` —
the `std::vector<std::unique_ptr<Ring>>` that all other ring
calculators consume via `RingAt(ri)`.

**String violation lines (verbatim).**

```cpp
// src/Protein.cpp:566-569
        // Build a name->atom_index map for this residue
        std::map<std::string, size_t> name_to_idx;
        for (size_t ai : res.atom_indices) {
            name_to_idx[atoms_[ai]->pdb_atom_name] = ai;
        }
```

```cpp
// src/Protein.cpp:571-582
        for (const auto& ring_def : aatype.rings) {
            // Check all ring atoms are present
            std::vector<size_t> atom_indices;
            bool all_present = true;
            for (const char* aname : ring_def.atom_names) {
                auto it = name_to_idx.find(aname);
                if (it == name_to_idx.end()) {
                    all_present = false;
                    break;
                }
                atom_indices.push_back(it->second);
            }
```

```cpp
// src/Protein.cpp:601-613 — HIS-tautomer fallback
                } else {
                    // PDB LOADING BOUNDARY fallback: protonation not yet
                    // detected (e.g., crystal structure with no H). Check
                    // hydrogen atom names to determine tautomer.
                    bool has_HD1 = name_to_idx.find("HD1") != name_to_idx.end();
                    bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();
                    if (has_HD1 && has_HE2) {
                        effective_type = RingTypeIndex::HisImidazole;
                    } else if (has_HD1) {
                        effective_type = RingTypeIndex::HidImidazole;
                    } else if (has_HE2) {
                        effective_type = RingTypeIndex::HieImidazole;
                    }
                }
```

**Typed translation.**

The string operations at line 568 read `atoms_[ai]->pdb_atom_name`
to build a map; in typed terms the operation is: for each atom in
the residue, query
`LegacyAmber().SemanticAt(atom_index).ring_position.primary.ring`;
group atoms by `(residue_index, RingSystemKind)` (the typed
substrate already records membership per atom).

The string operations at line 575 match
`ring_def.atom_names[k]` (a `const char*` like `"CG"`, `"CD1"`,
etc.) against the residue map keys; in typed terms, each
`ring_def.atom_names[k]` corresponds to a unique
`(RingPositionLabel, BranchAddress)` pair within the
`RingSystemKind` for the residue. Substrate-driven construction
walks the residue's atoms with
`LegacyAmber().SemanticAt(ai).ring_position.primary.ring` ==
target RingSystemKind, sorts by canonical
`RingPositionLabel` walk-order
(`Ipso → Ortho1 → Ortho2 → Meta1 → Meta2 → Para` for Phe; the
analogous walk for each ring per `SemanticEnums.h:528-589`), and
emits `atom_indices` in that order.

The HIS tautomer fallback at lines 604-613 reads
`name_to_idx.find("HD1")` and `name_to_idx.find("HE2")` to
discriminate HID / HIE / HIP. In typed terms, the substrate has
`RingPosition.primary.position == Heteroatom_NH` set on the
nitrogen carrying the proton; the H-bearing nitrogen's identity
follows from
`LegacyAmber().SemanticAt(N_atom).ring_position.primary.position ==
Heteroatom_NH` for the relevant N (Nδ1 for HID, Nε2 for HIE, both
for HIP). The protonation-detection result already carries
`protonation_variant_index`; the fallback path at lines 600-614
fires only when that index is `-1` (substrate not yet populated)
— a load-boundary case that becomes structurally impossible after
Bundle C lands at `FinalizeConstruction` order discipline.

**Ring API consumers in the substrate path.**

Verified that
`AminoAcidType::rings` is read in EXACTLY ONE site outside this
function:

```bash
grep -rn "rings\[\]\|aatype.rings" src/ --include='*.h' --include='*.cpp'
# Returns only: src/Protein.cpp:571
```

`AminoAcidRing::atom_names` is similarly only read at
`src/Protein.cpp:575`. Bundle C deletes both, plus the
`AminoAcidRing` struct entirely (`src/AminoAcidType.h:58-61`).

**Recommended classification:** INFRASTRUCTURE → MIGRATION TARGET.

**Per-calculator report needed?** Yes, with the Ring class hierarchy
bundled.
Filename: `Infrastructure-Ring-Detection-and-Class-Hierarchy.md`.


### 3.2 Ring class hierarchy

**File and line range.**
`src/Ring.h:1-205`. `src/Ring.cpp:1-65`.

**Purpose.** Defines the typed Ring class hierarchy, the
`RingGeometry` struct (per-conformation geometry), the
`RingAccumulated` struct (per-frame mutated state), and the
`CreateRing(RingTypeIndex)` factory.

**Class hierarchy.**

```text
Ring (abstract base)
├── SixMemberedRing (RingSizeValue = 6)
│   ├── PheBenzeneRing      I = -12.0,   JBLobeOffset 0.64 Å, Aromaticity::Full
│   ├── TyrPhenolRing       I = -11.28,  JBLobeOffset 0.64 Å, Aromaticity::Full
│   └── TrpBenzeneRing      I = -12.48,  JBLobeOffset 0.64 Å, Aromaticity::Full
├── FiveMemberedRing (RingSizeValue = 5)
│   ├── TrpPyrroleRing      I = -6.72,   JBLobeOffset 0.52 Å, Aromaticity::Reduced
│   ├── HisImidazoleRing    I = -5.16,   JBLobeOffset 0.50 Å, Aromaticity::Weak
│   ├── HidImidazoleRing    I = -5.16,   JBLobeOffset 0.50 Å, Aromaticity::Weak
│   └── HieImidazoleRing    I = -5.16,   JBLobeOffset 0.50 Å, Aromaticity::Weak
└── FusedRing
    └── IndolePerimeterRing I = -19.2,   JBLobeOffset 0.60 Å, Aromaticity::Full,
                            RingSizeValue = 9
```

**Virtual surface (`Ring.h:62-68`).**
- `Intensity() const` — calibration-time value via
  `CalculatorConfig::Get(<key>)`.
- `LiteratureIntensity() const` — hardcoded literature reference.
- `JBLobeOffset() const` — per-type lobe offset for Johnson-Bovey.
- `NitrogenCount() const`.
- `Aromaticity() const → RingAromaticity` (Full, Reduced, Weak;
  defined in `Types.h`).
- `RingSizeValue() const`.
- `TypeName() const → const char*` — short name for diagnostics
  ("PHE", "TYR", "TRP6", "TRP5", "TRP9", "HIS", "HID", "HIE").

**Non-virtual.**
- `IsFused() const` — returns true when `fused_partner_index !=
  SIZE_MAX`.
- `TypeIndexAsInt() const` — returns
  `static_cast<int>(type_index)`; the integer index into per-type
  arrays.
- `ComputeGeometry(positions) → RingGeometry` (`Ring.cpp:6-47`)
  — SVD-based normal, centroid, vertex collection, mean-radius;
  cross-product orientation fix to keep normal sign consistent
  with the right-hand rule on the first three vertices in
  `atom_indices` order.

**Identity fields (`Ring.h:53-57`).**
- `std::vector<size_t> atom_indices` — indices into
  `Protein::atoms_`.
- `RingTypeIndex type_index`.
- `size_t parent_residue_index`.
- `int parent_residue_number`.
- `size_t fused_partner_index = SIZE_MAX` — for TRP fused
  pyrrole/benzene Rings; SIZE_MAX otherwise.

**Pro pyrrolidine extension scope.** Adding `ProPyrrolidineRing
: public FiveMemberedRing` with overrides `Intensity() = 0`,
`LiteratureIntensity() = 0` (Joule & Mills 2010 ch. 7 — saturated
heterocycles do not carry circulating π current),
`JBLobeOffset() = 0`, `NitrogenCount() = 1`,
`Aromaticity() = RingAromaticity::None`, `RingSizeValue() = 5`,
`TypeName() = "PRO"`. Plus extending `RingTypeIndex` enum
(`Types.h:175-185`) with `ProPyrrolidine = 8` and incrementing
`Count = 9`. Plus a switch case in `CreateRing` factory
(`Ring.cpp:50-62`).

**Recommended classification:** INFRASTRUCTURE.

**Per-calculator report needed?** Yes; bundled with
DetectAromaticRings.
Filename: `Infrastructure-Ring-Detection-and-Class-Hierarchy.md`.


### 3.3 KernelEvaluationFilter / RingBondedExclusionFilter

**File and line range.**
`src/KernelEvaluationFilter.h:1-466` (header — defines the filter
ABC + concrete filters + `KernelFilterSet` aggregator).
`src/KernelEvaluationFilter.cpp:1-50+` — `RingBondedExclusionFilter`
constructor at line 22 walks the Protein topology to build per-ring
exclusion sets.

**Purpose.** Encapsulates physics-grounded filter rejection criteria
for kernel evaluations. Five concrete filters: `MinDistanceFilter`,
`DipolarNearFieldFilter`, `SelfSourceFilter`,
`SequentialExclusionFilter`, `RingBondedExclusionFilter`.
`KernelFilterSet` aggregates them; calculators add filters at
construction time and call `AcceptAll(ctx)` per evaluation.

**RingBondedExclusionFilter physics.** Per
`KernelEvaluationFilter.h:293-315`:
> An atom that is a vertex of a ring, or covalently bonded to a
> ring vertex, is in the through-bond regime for that ring's
> field. The through-space multipole expansion (dipolar,
> quadrupole, surface integral, dispersion) does not model
> through-bond electronic coupling.

The filter's constructor (line 22 of the .cpp) walks each ring's
`atom_indices` and each vertex's `bond_indices` to build per-ring
exclusion sets.

**Consumers (filter constructions).**
- `BiotSavartResult.cpp:151`.
- `HaighMallionResult.cpp:216`.
- `PiQuadrupoleResult.cpp:138`.
- `RingSusceptibilityResult.cpp:132`.

`DispersionResult.cpp` does NOT use this filter; it implements
equivalent topology locally at `DispersionResult.cpp:136-149`
(`BondedToVertices` helper) and stores per-ring sets in
`std::vector<std::set<size_t>> ring_bonded` at line 202.

**Recommended classification:** INFRASTRUCTURE.

**Per-calculator report needed?** Yes.
Filename: `Infrastructure-KernelFilterSet.md`.

**Bundle C implication.** No direct change. RingBondedExclusion
walks `Ring::atom_indices` regardless of how those indices were
obtained (string-matching today, substrate-driven after Bundle C);
bit-identity for non-Pro rings is the gate. The Pro
ProPyrrolidineRing simply gains an exclusion set entry — every Pro
residue's own backbone atoms (N, Cα, C, O, H, HA) and any
Cβ/Cγ/Cδ neighbours flow through the same logic.


## §4 — Incidental references

50+ files contain peripheral ring references that do not require a
per-calculator report. Clustering by file. Each row lists the file
and a one-line classification.

| File                                        | Line(s)                  | Classification                                                                                       |
| ------------------------------------------- | ------------------------ | ---------------------------------------------------------------------------------------------------- |
| `src/AIMNet2Result.{h,cpp}`                 | comments only            | Mentions ring chemistry in physics commentary; no executable ring code                                |
| `src/AmberChargeResolver.{h,cpp}`           | comments + variant       | "ring" appears in HID/HIE/HIP variant comments; no ring iteration                                    |
| `src/AmberLeapInput.{h,cpp}`                | comments                 | Ring-aware in tleap input templates; no Ring object access                                           |
| `src/AmberPreparedChargeSource.{h,cpp}`     | comments                 | Mentions ring atoms in charge resolution comments                                                    |
| `src/AminoAcidType.cpp`                     | 186-188, 201             | TRP `rings[]` table entries (Benzene/Pyrrole/Perimeter atom-name lists); TYR rings entry. Bundle C deletes |
| `src/AminoAcidType.h`                       | 58-61                    | `AminoAcidRing` struct definition (`type_index + atom_names`). Bundle C deletes                       |
| `src/ApbsFieldResult.{h,cpp}`               | comments                 | Ring atoms mentioned in residue iteration commentary                                                  |
| `src/Atom.{h,cpp}`                          | comments                 | Bond-graph code mentions ring; no Ring object access                                                  |
| `src/AtomEvent.h`                           | comments                 | Event-kind taxonomy mentions ring-flip events                                                         |
| `src/BondedEnergyResult.{h,cpp}`            | comments                 | Ring chemistry mentioned in bonded-energy commentary                                                  |
| `src/BondLengthStatsTrajectoryResult.{h,cpp}` | comments               | Trajectory bond-length stats; no ring access                                                          |
| `src/BuildResult.h`                         | comments                 | Ring-related diagnostic categories in build commentary                                                 |
| `src/CalculatorConfig.cpp`                  | 28-37, 51, 61-65         | Per-ring-type intensity + lobe offset config keys + ring_current_spatial_cutoff + ring_proximity_shell_{1,2,3,4} keys |
| `src/CalculatorConfig.h`                    | comments                 | Documents the keys                                                                                    |
| `src/ChargeAssignmentResult.{h,cpp}`        | comments                 | Ring-aware variant comments                                                                           |
| `src/ChargeSource.{h,cpp}`                  | comments                 | Ring-aware variant resolution                                                                         |
| `src/ChiRotamerSelectionTrajectoryResult.{h,cpp}` | header comment      | Mentions Pro ring rotamer selection                                                                   |
| `src/ConformationAtom.h`                    | 31-101, 121-122, 127-129, 131-134, 244-249 | `RingNeighbourhood` struct definition; `ring_neighbours` vector; per-type 8-arrays for BS/HM/PiQuad/Disp; ring proximity counts |
| `src/ConformationResult.{h,cpp}`            | comments                 | Ring kernels named in dependency commentary                                                           |
| `src/CoulombResult.cpp`                     | 39, 57, 77-81            | **Reads `Protein::RingCount()` + `Protein::RingAt(ri).atom_indices` to build per-atom `is_aromatic_atom[N]` mask for source classification (backbone / aromatic / sidechain).** Indirect ring consumer; not for ring-current physics |
| `src/CoulombResult.h`                       | comments                 | Documents the aromatic stratification                                                                 |
| `src/CovalentTopology.{h,cpp}`              | comments                 | Ring-related disulfide / aromatic mentions                                                            |
| `src/DemoResult.{h,cpp}`                    | comments                 | Demo result mentions ring type                                                                        |
| `src/DsspResult.{h,cpp}`                    | comments                 | DSSP secondary-structure mentions ring                                                                |
| `src/EeqResult.{h,cpp}`                     | comments                 | EEQ charge result mentions aromatic                                                                   |
| `src/EnrichmentResult.{h,cpp}`              | comments + ring class    | Enrichment categorises atoms; mentions ring class membership                                          |
| `src/ForceFieldChargeTable.{h,cpp}`         | comments                 | FF charge table includes per-residue ring atoms                                                       |
| `src/FullSystemReader.{h,cpp}`              | comments                 | Reader mentions ring-related TPR data                                                                 |
| `src/GeometryChoice.h`                      | comments                 | GeometryChoice taxonomy mentions ring                                                                 |
| `src/GeometryResult.{h,cpp}`                | h:22, cpp:17-19, 63-89, 96-97 | `RingGeometryAt(ri)` accessor; `Compute()` populates `conf.ring_geometries` via `Ring::ComputeGeometry`; populates `rings_by_type[type_index]`; populates ring-pair geometry |
| `src/GromacsEnergyResult.{h,cpp}`           | comments                 | TRR-energy mentions aromatic energetics                                                               |
| `src/GromacsFrameHandler.{h,cpp}`           | comments                 | Frame loader, no Ring access                                                                          |
| `src/GromacsFramePullResult.h`              | comments                 | Per-frame catch-all; no Ring                                                                          |
| `src/GromacsToAmberReadbackBlock.{h,cpp}`   | comments                 | Compiler-trace; no Ring access                                                                        |
| `src/HBondResult.{h,cpp}`                   | comments                 | H-bond physics mentions ring near-field                                                               |
| `src/HydrationGeometryResult.{h,cpp}`       | comments                 | Hydration around aromatic ring atoms                                                                  |
| `src/HydrationShellResult.{h,cpp}`          | comments                 | Hydration shell for ring atoms                                                                        |
| `src/JobSpec.{h,cpp}`                       | comments                 | Job spec mentions ring-current calculations                                                           |
| `src/KamlProtonator.{h,cpp}`                | comments                 | KaML protonator handles ring N-H tautomers                                                            |
| `src/LegacyAmberTopology.{h,cpp}`           | comments + substrate     | The runtime carrier of typed `ring_position` per atom; substrate composition                          |
| `src/McConnellResult.{h,cpp}`               | comments only            | Bond-anisotropy (NOT ring); mentions ring susceptibility as related physics                            |
| `src/MolecularGraphResult.{h,cpp}`          | comments                 | Graph result mentions ring detection                                                                  |
| `src/MopacCoulombResult.cpp`                | 53, 73-75                | **Same pattern as CoulombResult: builds per-atom aromatic mask from `Protein::RingAt(ri).atom_indices`.** Indirect consumer |
| `src/MopacCoulombResult.h`                  | comments                 |                                                                                                       |
| `src/MopacMcConnellResult.{h,cpp}`          | comments                 | Bond-anisotropy with MOPAC bond orders                                                                |
| `src/MopacResult.{h,cpp}`                   | comments                 | MOPAC subprocess wrapper; no Ring access                                                              |
| `src/MutationDeltaResult.{h,cpp}`           | comments                 | WT-Ala delta; mentions per-ring delta in commentary                                                   |
| `src/NamingRegistry.{h,cpp}`                | comments                 | Ring-aware FF name resolution                                                                         |
| `src/nmr_extract.cpp`                       | comments + dispatch      | Main dispatcher; ring calculators mentioned in registration                                           |
| `src/NpyWriter.h`                           | comments                 | NPY emitter; no Ring                                                                                  |
| `src/OperationLog.{h,cpp}`                  | log categories           | Log categories include `LogCalcBiotSavart`, `LogCalcHaighMal`, etc.                                   |
| `src/OperationRunner.{h,cpp}`               | comments                 | Operation orchestration                                                                               |
| `src/OrcaRunLoader.{h,cpp}`                 | comments                 | ORCA loader; no Ring                                                                                  |
| `src/OrcaShieldingResult.{h,cpp}`           | comments                 | DFT shielding loader; no Ring                                                                         |
| `src/PdbFileReader.{h,cpp}`                 | comments                 | PDB reader; no Ring                                                                                   |
| `src/PhysicalConstants.h`                   | comments                 | `BIOT_SAVART_PREFACTOR`, `NANOAMPERES_TO_AMPERES`, `PPM_FACTOR` documented for ring physics            |
| `src/PositionsTimeSeriesTrajectoryResult.{h,cpp}` | comments           | Position TR; no Ring                                                                                  |
| `src/PropkaProtonator.{h,cpp}`              | comments                 | Protonation; mentions HIS-ring tautomers                                                              |
| `src/Protein.h`                             | 14, 67-72, 203, 231      | `#include "Ring.h"`; `RingCount()`, `RingAt(ri)`, `Rings()` accessors; `DetectAromaticRings()` decl; `rings_` vector |
| `src/Protein.cpp`                           | 548-646                  | The migration target (§3.1); other ring mentions are comments + `RingCount()` returns                  |
| `src/ProteinBuildContext.h`                 | comments                 | Build context; no Ring access                                                                         |
| `src/ProteinConformation.h`                 | 112                      | `std::vector<RingGeometry> ring_geometries` per-conformation field (mutated by GeometryResult)        |
| `src/ProteinConformation.cpp`               | comments                 | Conformation construction                                                                             |
| `src/ProteinTopology.h`                     | comments                 | Topology ABC                                                                                          |
| `src/ProtonationDetectionResult.{h,cpp}`    | comments + variant       | HIS / HIP detection                                                                                   |
| `src/ProtonationState.{h,cpp}`              | comments                 | Variant state; no Ring                                                                                |
| `src/Protonator.h`                          | comments                 | ABC                                                                                                   |
| `src/RecordBag.h`                           | comments                 | Record bag template                                                                                   |
| `src/ReduceProtonation.{h,cpp}`             | comments                 | Reduce wrapper                                                                                        |
| `src/Residue.h`                             | comments                 | Residue identity                                                                                      |
| `src/RuntimeEnvironment.{h,cpp}`            | comments                 | RuntimeEnvironment; no Ring                                                                           |
| `src/RunConfiguration.{h,cpp}`              | comments + dispatch      | Per-frame config; ring calculators registered                                                         |
| `src/SasaResult.{h,cpp}`                    | comments                 | SASA; aromatic-vs-aliphatic mentions                                                                  |
| `src/SelectionRecord.h`                     | comments                 | Selection record                                                                                      |
| `src/SemanticEnums.h`                       | 484-632, 799             | The typed ring substrate vocabulary itself (`RingSystemKind`, `RingPositionLabel`, `RingMembership`, `RingPosition` definitions); reachable from runtime calculators via `LegacyAmber().SemanticAt(ai).ring_position`. NOT a calculator consumer; the substrate authority |
| `src/Session.{h,cpp}`                       | comments                 | Session orchestration                                                                                 |
| `src/SpatialIndexResult.cpp`                | 32-40, 135-141           | Builds KDTree over `conf.ring_geometries[ri].center`; `RingsWithinRadius(point, radius)` returns ring indices; ring-current calculators consume                  |
| `src/SpatialIndexResult.h`                  | 69                       | `RingsWithinRadius` decl                                                                              |
| `src/Trajectory.{h,cpp}`                    | comments                 | Trajectory orchestration                                                                              |
| `src/TrajectoryAtom.h`                      | comments                 | TrajectoryAtom field group; no direct ring                                                            |
| `src/TrajectoryProtein.{h,cpp}`             | comments                 | Trajectory-scope protein                                                                              |
| `src/TrajectoryResult.{h,cpp}`              | comments                 | TR ABC                                                                                                |
| `src/Types.{h,cpp}`                         | h:172-200, 184           | `RingTypeIndex` enum + `Count = 8` constant + `RingTypeName()` switch; `RingAromaticity` enum         |
| `src/WaterFieldResult.{h,cpp}`              | comments                 | Water field; no Ring                                                                                  |

Notable indirect consumers (CoulombResult, MopacCoulombResult): both
read `Protein::RingAt(ri).atom_indices` to build per-atom
classification masks, NOT for ring-current physics. These are
SUBSTANTIVE consumers of the Ring vector's identity surface but do
NOT need a per-calculator ring-investigation report — their ring
usage is a single-line membership check
(`is_aromatic_atom[ai] = true`), unaffected by Pro ring addition
(Pro has `is_aromatic = false` in `AminoAcidType`, so the source
class is unchanged; Pro ring atoms are still classified as
backbone/sidechain, never "aromatic source"). Documented here for
completeness.


## §5 — Test fixture list

60 test files contain ring references; all are test code consuming
the runtime ring API for assertions or fixture setup.

| File                                        | Classification                                                                                       |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `tests/test_ring_hierarchy.cpp`             | Direct ring-class hierarchy tests (subclass virtuals)                                                |
| `tests/test_biot_savart_result.cpp`         | BS calculator tests; ring assertions                                                                 |
| `tests/test_haigh_mallion_result.cpp`       | HM calculator tests                                                                                  |
| `tests/test_ring_susceptibility_result.cpp` | RingChi calculator tests                                                                             |
| `tests/test_pi_quadrupole_result.cpp`       | PiQuad calculator tests                                                                              |
| `tests/test_dispersion_result.cpp`          | Dispersion calculator tests                                                                          |
| `tests/test_batch_biot_savart_haigh_mallion.cpp` | Batch BS+HM tests across 100+ proteins                                                          |
| `tests/test_batch_coulomb_ringchi.cpp`      | Batch Coulomb+RingChi tests                                                                          |
| `tests/test_batch_piquad_disp.cpp`          | Batch PiQuad+Disp tests                                                                              |
| `tests/test_batch_mcconnell.cpp`            | Batch MC tests; mentions ring                                                                        |
| `tests/test_object_model.cpp`               | Ring class object-model conformance tests                                                            |
| `tests/test_amino_acid.cpp`                 | AminoAcidType including `rings[]` membership tests                                                   |
| `tests/test_geometry_result.cpp`            | RingGeometry tests                                                                                   |
| `tests/test_calculator_config.cpp`          | Per-ring intensity / lobe-offset config keys                                                         |
| `tests/test_string_barrier.cpp`             | String-barrier tests; verifies Ring class doesn't include cifpp/RDKit                                |
| `tests/test_traversal_dump.cpp`             | Traversal dump including ring data                                                                   |
| `tests/test_full_pipeline.cpp`              | End-to-end pipeline; loads rings                                                                     |
| `tests/test_pipeline_and_sample.cpp`        | Pipeline + sampling; rings present                                                                   |
| `tests/test_pdb_loading.cpp`                | PDB loader; ring detection                                                                           |
| `tests/test_amber_streaming.cpp`            | AMBER streaming; ring construction                                                                   |
| `tests/test_amber_trajectory.cpp`           | AMBER trajectory; per-frame rings                                                                    |
| `tests/test_two_conformations.cpp`          | Two-conformation comparison; ring identity preserved                                                 |
| `tests/test_protonation_pipeline.cpp`       | HIS variant detection affects ring type                                                              |
| `tests/test_protonation_detection.cpp`      | HIS detection ring tautomer assignment                                                               |
| `tests/test_main.cpp`                       | Test runner; no ring code                                                                            |
| `tests/test_calculation_runner.cpp`         | Calculation runner; runs ring calculators                                                            |
| `tests/test_write_features.cpp`             | NPY writers including per-type ring arrays                                                           |
| `tests/test_atom_flat.cpp`                  | Flat atom record; no ring                                                                            |
| `tests/test_smoke.cpp`                      | Smoke test; loads ring fixtures                                                                      |
| `tests/test_runtime_environment.cpp`        | RuntimeEnvironment tests                                                                             |
| `tests/test_dssp_result.cpp`                | DSSP test; aromatic SS                                                                               |
| `tests/test_naming_registry.cpp`            | FF-name registry; ring-bearing residues                                                              |
| `tests/test_amber_charge_resolver.cpp`      | AMBER charge resolution                                                                              |
| `tests/test_amber_prepared_charge_source.cpp` | Prepared AMBER charge                                                                              |
| `tests/test_amber_leap_input.cpp`           | tleap input including rings                                                                          |
| `tests/test_apbs_ff14sb.cpp`                | APBS field test                                                                                      |
| `tests/test_coulomb_result.cpp`             | Coulomb test; aromatic stratification                                                                |
| `tests/test_mcconnell_result.cpp`           | MC test; mentions ring                                                                               |
| `tests/test_hbond_result.cpp`               | H-bond test                                                                                          |
| `tests/test_demo_result.cpp`                | Demo test                                                                                            |
| `tests/test_foundation_results.cpp`         | Foundation Results including GeometryResult ring init                                                |
| `tests/test_job_spec.cpp`                   | Job spec test                                                                                        |
| `tests/test_mutation_delta.cpp`             | WT-Ala delta with rings                                                                              |
| `tests/test_spherical_tensor.cpp`           | Tensor decomposition tests; mentions T2 from ring kernels                                            |
| `tests/test_protonation_pipeline.cpp`       | (duplicate row-name; covered above)                                                                  |
| `tests/pass0_demo.cpp`                      | Pass 0 demo                                                                                          |
| `tests/TestEnvironment.{h,cpp}`             | Test fixture infrastructure                                                                          |
| `tests/BlessCompare.{h,cpp}`                | Bless-compare; tolerates ring-output drift                                                           |
| `tests/bones/test_job_spec_fleet.cpp`       | Bones-fleet test                                                                                     |

Bundle C verification: bless-compare on the standard fixture set
(per `tests/golden/blessed/`) is the gating non-Pro-ring bit-identity
check; new Pro Ring adds new NPY columns or alters existing column
shapes (`per_type_*` arrays going from 8 → 9 elements). The bless
policy in `tests/golden/blessed/bless_policy.toml` plus the drift
table in `tests/golden/blessed/BLESS_NOTES.md` document the
acceptable drift envelope; Pro Ring landing requires explicit
bless re-set on the affected per-type arrays.


## §6 — Per-calculator-report-needed list

The per-calculator phase produces these reports at
`spec/plan/ring-investigation-2026-05-06/<filename>`:

```text
BiotSavartResult.md
HaighMallionResult.md
RingSusceptibilityResult.md
PiQuadrupoleResult.md
DispersionResult.md
Infrastructure-Ring-Detection-and-Class-Hierarchy.md   (DetectAromaticRings + Ring class hierarchy)
Infrastructure-KernelFilterSet.md                       (RingBondedExclusionFilter + KernelFilterSet)
Trajectory-BS-Aggregators.md                            (BsWelford, BsShieldingTimeSeries, BsT0Autocorrelation, BsAnomalousAtomMarker)
```

8 reports total. The trajectory-scope BS aggregators are bundled
because they consume only `bs_shielding_contribution` (a
ConformationAtom field), not Ring directly — their analyses are
parallel and benefit from being in one document.


## §7 — Coverage diagnostic

Verification command (re-run by user):

```bash
cd /shared/2026Thesis/nmr-shielding/

grep -rni "ring" src/ --include='*.h' --include='*.cpp' | wc -l
# Expected: 2908 matches across 160 files

grep -rln -i "ring" src/ --include='*.h' --include='*.cpp' | wc -l
# Expected: 160 files

grep -rnE "RingTypeIndex|RingNeighbourhood|RingSystemKind|RingPositionLabel|RingMembership|RingPosition\b" \
    src/ --include='*.h' --include='*.cpp' | wc -l
# Expected: 608 matches in 22 files

grep -rni "ring" tools/ --include='*.cpp' --include='*.md' --include='*.txt' | wc -l
# Expected: 425 matches in 4 files (build_semantic_tables.cpp dominates)

grep -ni "ring" src/generated/LegacyAmberSemanticTables.cpp | wc -l
# Expected: 507 substrate-table emissions
```

Reproduced cleanly. The 2908/160/608/22 numbers from the prior
agent's analysis are independently confirmed.

**Coverage gap note.** Files that include `Ring.h` via header
declarations but never reference ring data in their .cpp do not
appear in the case-insensitive grep above (the include is in the
.h, the string `ring` appears nowhere in the .cpp). The 22-file
typed-grep is the more precise surface; it identifies the calculator
+ infrastructure files that actively typecheck against ring
substrate. The full 160-file count includes commentary-only and
indirect mentions covered in §4's table.


## §8 — Open questions for user input

These three questions remain open from the prior agent's analysis;
they are NOT proposed answers. The user calls them.

### Q1. Pro pyrrolidine atom ordering

Two candidate orderings for the new `ProPyrrolidineRing`'s
`atom_indices`:

(a) **Residue-walk order** (per README §C constraint): `N → Cα →
Cβ → Cγ → Cδ`. Stable across proteins; supports puckering
descriptors via the dihedral N-Cα-Cβ-Cγ-Cδ.

(b) **Substrate's `RingPositionLabel` walk-order**. Today's
`RingPositionLabel` for Pro is `Saturated` for all 5 atoms (per
`src/SemanticEnums.h:587-588`); the substrate does not currently
disambiguate Pro positions further. Adopting (b) would require
extending `RingPositionLabel` with Pro-specific labels, or relying
on `Locant` (Pro's Cα is `Locant::None`, Cβ/γ/δ are
`Beta/Gamma/Delta`).

The user picks. Bundle C drafting depends on the choice:
- (a) lets ConstructRingsFromSubstrate stable-sort Pro atoms by
  (BackboneRole, Locant) without substrate extension.
- (b) requires substrate extension (new RingPositionLabel values
  or a new mechanism) and a generator re-run.

### Q2. TRP perimeter representation

Two candidate implementations for `IndolePerimeterRing` (the
9-atom perimeter):

(a) **Synthesise at runtime** in
`ConstructRingsFromSubstrate`: union of `Indole_Trp_5` and
`Indole_Trp_6` atoms minus shared edge, with cyclic walk defined
locally.

(b) **Extend substrate** with `RingSystemKind::Indole_Trp_9`. The
generator emits perimeter atoms with secondary
`RingPosition.position` tagged with perimeter labels.

The user picks. Per the README's planned-calculator implications
(§C), (b) is more reusable for future per-position perimeter
calculators (e.g. shielding contribution from atoms at Trp
perimeter ortho1 vs meta2 across the trajectory). (a) is simpler
but pushes perimeter synthesis logic into runtime.

### Q3. GeometryChoice records for Pro Intensity = 0 rings

Should calculators that iterate Pro Ring (BS, HM, RingChi, PiQuad,
Disp) emit GeometryChoice records when the contribution is
mathematically zero? Two options:

(a) **Record** every Pro ring iteration with
`AddNumber("intensity", 0.0, "nA")`, `AddRing(... Source/Included)`
— documentary value; future calculators can ask "what was
considered? what was rejected?" by reading the choice log.

(b) **Skip recording** Pro ring iterations to keep the choice log
tidy.

The user picks. Per README §10, the project's discipline is
record-everything; (a) aligns with PATTERNS.md "When implementing
a calculator". (b) optimises for log size at the cost of audit
coverage.


## §9 — String discipline finding

The investigation confirms exactly ONE non-loader site reads atom-
name strings for ring chemistry inference: **DetectAromaticRings**
(`src/Protein.cpp:548-646`).

### The single-site claim, verified

The grep:
```bash
grep -rn "rings\[\]\|aatype.rings" src/ --include='*.h' --include='*.cpp'
```
returns ONE match: `src/Protein.cpp:571`. Similarly,
```bash
grep -rn "ring_def\|atom_names" src/ --include='*.h' --include='*.cpp' \
  | grep -v "AmberPreparedChargeSource\|OrcaRunLoader\|LegacyAmberTopology\|TrajectoryProtein"
```
returns matches only in `src/Protein.cpp:571,575` (the migration
target) and `src/AminoAcidType.h:60` (the struct definition).

The remaining `atom_names` matches in
`AmberPreparedChargeSource.cpp`, `OrcaRunLoader.cpp`,
`LegacyAmberTopology.cpp`, and `TrajectoryProtein.cpp` are
unrelated string surfaces (PRMTOP atom-name vectors for charge
resolution, ORCA atom labels, LegacyAmberTopology composition
loop, H5 emission for visualisation). Not Bundle C migration
targets.

### The five calculator string-discipline checks

Per inspection of all 503 + 432 + 270 + 292 + 423 = **1920 lines** of
calculator code (`BiotSavartResult.cpp`, `HaighMallionResult.cpp`,
`RingSusceptibilityResult.cpp`, `PiQuadrupoleResult.cpp`,
`DispersionResult.cpp`), no calculator reads
`atom.pdb_atom_name`. Verified:
```bash
grep -rn "atom.pdb_atom_name\|pdb_atom_name ==" \
  src/BiotSavartResult.cpp src/HaighMallionResult.cpp \
  src/PiQuadrupoleResult.cpp src/RingSusceptibilityResult.cpp \
  src/DispersionResult.cpp
# Returns no matches.
```

### Bundle C clearance target

`ConstructRingsFromSubstrate` reads typed `ring_position` directly
(via `LegacyAmber().SemanticAt(ai).ring_position.primary.ring`).
`DetectAromaticRings` is deleted. The const-char-pointer arrays in
`AminoAcidType::rings[].atom_names` (and the `AminoAcidRing`
struct definition itself at `AminoAcidType.h:58-61`) become
removable — no other consumer exists in the runtime library.

The HIS-tautomer fallback at `Protein.cpp:600-614` similarly
becomes structurally inaccessible after Bundle C lands at
`FinalizeConstruction` order discipline (the substrate has the
typed `Heteroatom_NH` per-N record before
ConstructRingsFromSubstrate runs; the protonation-detection-result
+ substrate composition both fire before the ring-construction
pass).

The framing from README §11 (fail-loud is scientific discipline,
not engineering preference): if a residue's substrate has no
`ring_position` set on atoms expected to be in a ring,
ConstructRingsFromSubstrate must fail-loud (FATAL + abort) rather
than silently degrade. Today's DetectAromaticRings silently
continues if `ring_def.atom_names` doesn't match (`if
(!all_present) continue;` at line 583); Bundle C's substrate-
driven version must enforce the stronger invariant.
