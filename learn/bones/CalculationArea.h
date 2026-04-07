#pragma once
//
// CalculationArea: every named geometric, topological, or parametric
// choice point in the NMR shielding tensor calculation pipeline.
//
// The 10 calculators make ~40 decisions using hardcoded constants:
// "is this source close enough?", "is this atom inside the source
// distribution?", "what is this ring type's current intensity?"
// Each decision is now a named, typed, documented object.
//
// The ABC carries identity and documentation. Subtypes add
// pattern-specific query interfaces. There is at most one
// inheritance layer between the ABC and concrete types.
//
// First pass: every area reproduces the current code exactly.
// The value is provenance and visibility, not behavior change.
//
// Full catalogue: spec/CALCULATION_AREA_CATALOGUE.csv
// Kernel physics: GEOMETRIC_KERNEL_CATALOGUE.md
//

#include "Types.h"
#include <climits>
#include <cstddef>
#include <set>
#include <vector>

namespace nmr {

class Protein;


// ============================================================================
// DomainKind: what kind of decision this area makes.
//
// Used by the UI to select visualization (draw a sphere, highlight
// atoms, show a scalar field) and by the calibration pipeline to
// group areas by pattern.
// ============================================================================

enum class DomainKind {
    Spatial,            // fixed-radius inclusion/exclusion
    Shell,              // cumulative bin edge for counting
    SourceRelative,     // boundary scaled by source extent
    Topological,        // bond-graph membership check
    Sequence,           // residue sequence separation
    Switching,          // smooth onset/cutoff pair
    Decay,              // exponential decay length
    RingMagnitude,      // ring current intensity (nA)
    RingGeometry,       // ring geometric parameter (lobe height)
    Numerical,          // adaptive refinement or clamp
    ValueThreshold,     // floor below which contributions are skipped
    Sentinel,           // marker for missing data
    WholeDomain,        // documents absence of cutoff
};


// ============================================================================
// CalculationArea (ABC)
//
// Every choice point in the calculator pipeline. The ABC carries
// identity and classification only. Subtypes own their query
// interfaces — there is no lowest-common-denominator API.
//
// This is NOT an evaluation interface. CalculationAreas do not
// Accept() or Reject() — KernelEvaluationFilter does that.
// CalculationAreas are the CONSTANTS that filters and calculators
// consume. The area says "the ring horizon is 15 Angstroms."
// The filter says "reject if beyond it."
// ============================================================================

class CalculationArea {
public:
    virtual ~CalculationArea() = default;

    // Short CamelCase identifier matching the catalogue name.
    // Example: "RingHorizon", "MultipoleInnerBoundary".
    virtual const char* Name() const = 0;

    // One-sentence physics reason this choice point exists and
    // why the value was chosen. Not a generic label — the actual
    // justification. Example: "Ring current B-field decays as
    // 1/r^3; at 15A contribution is <0.03% of the 3A value."
    virtual const char* Description() const = 0;

    // What kind of decision this is.
    virtual DomainKind Domain() const = 0;
};


// ============================================================================
// RadialThreshold
//
// A fixed-distance sphere around a source point. The evaluation at
// an atom either passes or fails based on distance compared to the
// threshold. Covers outer horizons (RingHorizon 15A, BondAnisotropy-
// Horizon 10A), inner guards (SingularityGuard 0.1A), proximity
// shells, and matching tolerances.
//
// The Sense distinguishes inclusion from exclusion:
//   IncludeWithin: atom must be CLOSER than Radius() (outer horizon)
//   ExcludeWithin: atom must be FARTHER than Radius() (inner guard)
//
// The center of the sphere comes from the evaluation context (source
// center, bond midpoint, ring center) — it is NOT stored on the area.
// The area holds only the radius and the comparison sense.
//
// UI visualization: draw a sphere of Radius() at the source center,
// colored by Sense (green for include, red for exclude).
//
// Catalogue instances (14):
//   RingHorizon           15.0 A  IncludeWithin  5 ring calculators
//   BondAnisotropyHorizon 10.0 A  IncludeWithin  McConnell
//   MopacBondAnisHorizon  10.0 A  IncludeWithin  MopacMcConnell (same value)
//   SingularityGuard       0.1 A  ExcludeWithin  all calculators
//   HBondMaxReach         50.0 A  IncludeWithin  HBond
//   HBondProximityShell    3.5 A  IncludeWithin  HBond (count feature)
//   PackingShell           8.0 A  IncludeWithin  feature (heavy atom count)
//   MutationMatchTol       0.5 A  IncludeWithin  MutationDelta (atom matching)
//   SpatialIndexHorizon   15.0 A  IncludeWithin  infrastructure (must >= max)
//
// ShellBoundary instances (4) use IncludeWithin but are cumulative
// bin edges for ring counting, not evaluation gates:
//   RingShellStrong        3.0 A
//   RingShellModerate      5.0 A
//   RingShellFarField      8.0 A
//   RingShellBackground   12.0 A
// ============================================================================

enum class Sense {
    IncludeWithin,   // accept if distance < radius (outer horizon)
    ExcludeWithin,   // accept if distance > radius (inner guard)
};

class RadialThreshold : public CalculationArea {
public:
    RadialThreshold(const char* name, const char* description,
                    DomainKind domain, double radius_angstroms, Sense sense);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The threshold distance in Angstroms.
    double Radius() const;

    // Whether atoms inside the radius are included or excluded.
    Sense ThresholdSense() const;

private:
    const char* name_;
    const char* description_;
    DomainKind domain_;
    double radius_;
    Sense sense_;
};


// ============================================================================
// SourceRelativeExclusion
//
// Inner boundary whose radius is a FRACTION of the source's spatial
// extent, not a fixed distance. The source extent is a measured
// property of each source (ring diameter, bond length, N...O distance),
// not a tuning parameter. The factor IS the tuning parameter.
//
// Physics: the multipole expansion (dipolar, quadrupole, surface
// integral) diverges when the field point is inside the source
// charge/current distribution. The expansion center is at the
// source midpoint; the distribution extends to ~half the source
// extent. Evaluating the kernel closer than this gives large
// numbers that are physically meaningless.
//
// Demonstrated effect (HBondResult, 466 proteins):
//   With filter:  max |T2| = 0.78 A^-3
//   Without:      max |T2| = 1908 A^-3
//
// The factor (0.5) is the choice point. It means "reject if closer
// than half the source extent." A smaller factor admits more atoms
// near the source boundary; a larger factor is more conservative.
// Jackson Ch. 4 (multipole convergence) justifies 0.5 as the
// natural boundary of a source distribution.
//
// Currently implemented as DipolarNearFieldFilter in
// KernelEvaluationFilter.h:140. Used by 6 calculators:
// BiotSavart, HaighMallion, McConnell, RingSusceptibility,
// PiQuadrupole, HBond.
//
// Catalogue instance (1):
//   MultipoleInnerBoundary  factor=0.5  6 calculators
// ============================================================================

class SourceRelativeExclusion : public CalculationArea {
public:
    SourceRelativeExclusion(const char* name, const char* description,
                            double factor);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The fraction of source extent that defines the exclusion
    // boundary. Accept if distance > factor * source_extent.
    double Factor() const;

private:
    const char* name_;
    const char* description_;
    double factor_;
};


// ============================================================================
// RingBondedExclusion
//
// Topological exclusion of atoms that are vertices of a ring or
// covalently bonded to a ring vertex. These atoms are in the
// through-bond regime for that ring's field — the through-space
// multipole expansion (dipolar, quadrupole, surface integral,
// dispersion) does not model through-bond electronic coupling.
// Evaluating the kernel at these atoms gives a large number that
// is physically meaningless.
//
// This is a TOPOLOGICAL check, not a distance proxy. It uses the
// bond graph to determine membership. Compare to the DipolarNear-
// FieldFilter which catches most of the same atoms by distance,
// but has a boundary case for ring atoms that sit exactly at the
// filter threshold.
//
// Construction: takes a const reference to the Protein whose
// topology defines which atoms are ring vertices and what they
// are bonded to. At construction, walks Protein -> Ring ->
// atom_indices -> bond graph to build a per-ring exclusion set.
//
// The TRP additivity proof depends on this filter: without it,
// atoms on the shared CD2-CE2 edge are evaluated by both TRP5
// and TRP6 surface integrals, producing a 13% excess in the HM
// ratio. With the filter, T0(TRP5)+T0(TRP6)/T0(TRP9) = 1.000
// exactly (verified on 720 proteins, 87,234 atom-TRP triplets).
//
// Currently implemented as RingBondedExclusionFilter in
// KernelEvaluationFilter.h:279. Used by 5 ring calculators:
// BiotSavart, HaighMallion, RingSusceptibility, PiQuadrupole,
// Dispersion.
//
// DispersionResult.cpp reimplements this walk independently
// (DispersionBondedExclusion in the catalogue). A future pass
// should unify them to use this CalculationArea.
//
// Catalogue instances (2, same pattern):
//   RingBondedExclusion       5 ring calculators
//   DispersionBondedExclusion Dispersion (reimplemented)
// ============================================================================

class RingBondedExclusion : public CalculationArea {
public:
    // Walks the Protein's topology once to build per-ring
    // exclusion sets. Each set contains the ring's vertices
    // and all atoms bonded to them.
    explicit RingBondedExclusion(const char* name, const char* description,
                                 const Protein& protein);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The set of atom indices excluded for a given ring.
    // Returns empty set if ring_index is out of range.
    const std::set<size_t>& ExcludedAtoms(size_t ring_index) const;

    // How many rings have exclusion sets.
    size_t RingCount() const;

private:
    const char* name_;
    const char* description_;
    std::vector<std::set<size_t>> ring_excluded_;
};


// ============================================================================
// SelfSourceExclusion
//
// Topological exclusion of an atom when it IS the source. The
// dipolar field of a source is undefined at the source itself —
// an atom that is an endpoint of a bond cannot be a field point
// for that bond's McConnell kernel.
//
// This is a type-based check: it compares atom indices, not
// distances. The physics justification is that the source atom's
// own electron cloud generates the field — it does not experience
// its own dipolar field as an external perturbation.
//
// Currently implemented as SelfSourceFilter in
// KernelEvaluationFilter.h:178. Used by 5 calculators:
// McConnell, Coulomb, HBond, MopacMcConnell, MopacCoulomb.
//
// Catalogue instance (1):
//   SelfSourceExclusion  5 calculators
// ============================================================================

class SelfSourceExclusion : public CalculationArea {
public:
    SelfSourceExclusion(const char* name, const char* description);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // No query method. SelfSourceExclusion holds no data —
    // it is pure identity. The comparison logic (atom_index
    // != source_atom_a, != source_atom_b) lives on
    // SelfSourceFilter in KernelEvaluationFilter.h.

private:
    const char* name_;
    const char* description_;
};


// ============================================================================
// ShellBoundary
//
// Cumulative bin edge for counting sources within distance shells.
// Unlike RadialThreshold (which gates evaluation), a ShellBoundary
// is a histogram boundary: "how many rings are within 3A? within 5A?"
// Each shell is cumulative — the 5A count includes rings also counted
// at 3A.
//
// The shell values define the ring proximity feature space for the
// learning model. The model sees ring density at multiple scales.
// The specific radii (3/5/8/12 A) are empirical, chosen to span
// from strong van der Waals contact (~2x vdW to ring center) through
// far-field onset (T0 sign agreement >99.3% beyond 8A).
//
// Currently applied in BiotSavartResult.cpp (lines 270-275) as
// distance_to_center <= RING_COUNT_SHELL_N comparisons.
//
// Catalogue instances (4):
//   RingShellStrong      3.0 A  ~2x vdW contact to ring center
//   RingShellModerate    5.0 A  intermediate ring current zone
//   RingShellFarField    8.0 A  BS/HM T0 sign >99.3% beyond here
//   RingShellBackground 12.0 A  negligible ring current for most atoms
// ============================================================================

class ShellBoundary : public CalculationArea {
public:
    ShellBoundary(const char* name, const char* description,
                  double bin_edge_angstroms);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The shell radius in Angstroms. Cumulative: count includes
    // all sources within this distance.
    double BinEdge() const;

private:
    const char* name_;
    const char* description_;
    double bin_edge_;
};


// ============================================================================
// SequenceGate
//
// Minimum residue sequence separation for a through-space kernel to
// be applicable. Atoms within a few residues share through-bond
// electronic coupling (inductive, hyperconjugation) that the
// through-space dipolar model does not capture. The dipolar kernel
// gives a geometric prediction at these distances, but the actual
// shielding perturbation is dominated by through-bond effects with
// a different angular structure.
//
// Currently implemented as SequentialExclusionFilter in
// KernelEvaluationFilter.h:221. Default min_separation = 2.
//
// Catalogue instance (1):
//   HBondSequenceGate  min_separation=2  HBond calculator
//   Reference: Cornilescu & Bax (2000)
// ============================================================================

class SequenceGate : public CalculationArea {
public:
    SequenceGate(const char* name, const char* description,
                 int min_separation);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // Minimum residue sequence separation. Accept if
    // |residue_i - residue_j| >= MinSeparation().
    int MinSeparation() const;

private:
    const char* name_;
    const char* description_;
    int min_separation_;
};


// ============================================================================
// SwitchingFunction
//
// A smooth onset/cutoff pair implementing a CHARMM-style switching
// function. Between Onset() and Cutoff(), contributions are
// continuously attenuated from 1.0 to 0.0. Below onset: full weight.
// Above cutoff: zero.
//
// Physics: the taper width (cutoff - onset) should approximate the
// thermal fluctuation of the distance (MD RMS ~0.7A). A hard cutoff
// creates discontinuities when atoms cross the boundary between MD
// frames. The switching function smooths this.
//
// Form (Brooks et al. 1983):
//   S(r) = (R_cut^2 - r^2)^2 * (R_cut^2 + 2r^2 - 3R_on^2)
//          / (R_cut^2 - R_on^2)^3
//
// Currently implemented in DispersionResult.cpp (lines 60-73)
// as DispSwitchingFunction with R_SWITCH=4.3A, R_CUT=5.0A.
//
// Catalogue instance (1):
//   DispersionTaper  onset=4.3A cutoff=5.0A  Dispersion calculator
//   Reference: Brooks et al. (1983) form; onset value empirical
// ============================================================================

class SwitchingFunction : public CalculationArea {
public:
    SwitchingFunction(const char* name, const char* description,
                      double onset_angstroms, double cutoff_angstroms);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // Distance where attenuation begins (full weight below this).
    double Onset() const;

    // Distance where contribution reaches zero.
    double Cutoff() const;

private:
    const char* name_;
    const char* description_;
    double onset_;
    double cutoff_;
};


// ============================================================================
// DecayFunction
//
// Exponential decay length for distance-weighted features.
// Weight = exp(-x / DecayLength()), where x is either a spatial
// distance (Angstroms) or a graph distance (bond hops).
//
// These are feature-engineering parameters, not evaluation gates.
// They weight contributions by proximity without hard cutoffs.
//
// Catalogue instances (2):
//   RingProximityDecay   tau=4.0 A     exp(-distance/4.0)
//     MutationDelta calculator. 4A ~ 2 backbone bond lengths.
//     Reference: Xu & Case (2001)?; likely empirical.
//   GraphDistanceDecay   tau=4.0 hops  exp(-hops/4.0)
//     MolecularGraph calculator. 4 hops ~ 1 residue.
//     Reference: empirical.
// ============================================================================

class DecayFunction : public CalculationArea {
public:
    DecayFunction(const char* name, const char* description,
                  double decay_length, const char* unit);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The characteristic decay length.
    double DecayLength() const;

    // "A" for Angstroms, "hops" for bond graph distance.
    const char* Unit() const;

private:
    const char* name_;
    const char* description_;
    double decay_length_;
    const char* unit_;
};


// ============================================================================
// RingCurrent
//
// Diamagnetic ring current intensity for one ring type. This is the
// I in sigma = I * G, where G is the intensity-independent geometric
// kernel from the Johnson-Bovey double-loop model.
//
// The geometric kernel G is evaluated at unit current (I = 1 nA).
// The intensity scales the kernel output uniformly — it affects T0
// magnitude but not T2 angular pattern. Whether I stays as a
// literature constant, becomes a learnable weight, or gets replaced
// by a calculated estimate from surrounding CalculationArea values
// is an open question for the calibration campaign.
//
// Currently defined as Intensity() virtual methods on Ring subclasses
// in Ring.h. Future: constexpr RING_TYPE_TABLE entries (see
// spec/RING_OBJECT_MODEL_PROPOSAL.md).
//
// Catalogue instances (6):
//   PheRingCurrent  -12.00 nA  Case (1995) JBNMR 6:341
//   TyrRingCurrent  -11.28 nA  Case (1995); Osapay & Case (1991)
//   Trp6RingCurrent -12.48 nA  Case (1995); Giessner-Prettre (1987)
//   Trp5RingCurrent  -6.72 nA  Case (1995)
//   Trp9RingCurrent -19.20 nA  = Trp5 + Trp6 exactly (verified 720 proteins)
//   HisRingCurrent   -5.16 nA  Case (1995); HIE/HID/HIP identical
// ============================================================================

class RingCurrent : public CalculationArea {
public:
    RingCurrent(const char* name, const char* description,
                RingTypeIndex ring_type, double intensity_nanoamperes);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // Which ring type this current applies to.
    RingTypeIndex RingType() const;

    // The ring current in nanoamperes (negative = diamagnetic).
    double Intensity() const;

private:
    const char* name_;
    const char* description_;
    RingTypeIndex ring_type_;
    double intensity_;
};


// ============================================================================
// LobeOffset
//
// Johnson-Bovey pi-cloud height for one ring type. This is the d in
// the JB double-loop model: two current loops at +/- d from the ring
// plane, each carrying I/2.
//
// The lobe offset directly controls the T2 angular pattern of the
// Biot-Savart kernel: larger d concentrates the B-field into sharper
// lobes above and below the ring. This is a GEOMETRIC parameter —
// it changes the shape of the kernel, not just its magnitude.
//
// If the T2 residual is systematically wrong for a ring type, the
// lobe offset is the first parameter to examine.
//
// Currently defined as JBLobeOffset() virtual methods on Ring
// subclasses in Ring.h.
//
// Catalogue instances (4):
//   SixMemberedLobeHeight 0.64 A  PHE/TYR/TRP6
//     Johnson & Bovey (1958); Haigh & Mallion (1979)
//   Trp5LobeHeight        0.52 A  TRP pyrrole
//     Giessner-Prettre & Pullman (1987)
//   Trp9LobeHeight        0.60 A  TRP indole perimeter
//     UNKNOWN — likely interpolated (weakest reference in catalogue)
//   HisLobeHeight         0.50 A  HIS imidazole (all protonation states)
//     Osapay & Case (1991)
// ============================================================================

class LobeOffset : public CalculationArea {
public:
    LobeOffset(const char* name, const char* description,
               RingTypeIndex ring_type, double offset_angstroms);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // Which ring type this offset applies to.
    RingTypeIndex RingType() const;

    // The pi-cloud height in Angstroms above/below ring plane.
    double Offset() const;

private:
    const char* name_;
    const char* description_;
    RingTypeIndex ring_type_;
    double offset_;
};


// ============================================================================
// NumericalAccuracy
//
// Proximity threshold that triggers adaptive refinement in a
// numerical integration. Not a physics decision — a numerical
// analysis decision about when the default quadrature is
// insufficiently accurate.
//
// Currently used in HaighMallionResult.cpp for the surface integral:
// the 7-point Gaussian quadrature is subdivided when the field point
// is close to the ring surface.
//
// Catalogue instances (2):
//   HmRefineNear   2.0 A  Subdivide once (L0 -> L1)
//   HmRefineClose  1.0 A  Subdivide twice (L1 -> L2)
//   Reference: standard adaptive quadrature practice
// ============================================================================

class NumericalAccuracy : public CalculationArea {
public:
    NumericalAccuracy(const char* name, const char* description,
                      double threshold_angstroms);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // Distance at which refinement triggers.
    double Threshold() const;

private:
    const char* name_;
    const char* description_;
    double threshold_;
};


// ============================================================================
// ValueClamp
//
// Ceiling on a computed value to prevent extreme outliers from
// propagating through the pipeline. Not physics — protective
// numerical hygiene against clash geometries.
//
// Catalogue instance (1):
//   EFieldSanityClamp  100.0 V/A  Coulomb/MopacCoulomb/APBS
//     Applied when E-field magnitude exceeds limit.
//     Reference: empirical
// ============================================================================

class ValueClamp : public CalculationArea {
public:
    ValueClamp(const char* name, const char* description,
               double limit, const char* unit);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The maximum allowed value.
    double Limit() const;

    // Unit of the clamped quantity (e.g. "V/A").
    const char* Unit() const;

private:
    const char* name_;
    const char* description_;
    double limit_;
    const char* unit_;
};


// ============================================================================
// ValueGate
//
// Floor below which a source contribution is skipped as numerically
// insignificant. These are noise floors, not physics decisions —
// the source exists but its contribution is below the precision of
// the calculation.
//
// Catalogue instances (3):
//   MopacBondOrderFloor  1e-6   dimensionless
//     Skip bonds with negligible electron sharing.
//     MopacMcConnellResult.cpp:147
//   ChargeNoiseFloor     1e-15  elementary charges
//     Skip atoms with negligible partial charge.
//     CoulombResult.cpp:142, MopacCoulombResult.cpp:126
//   NearZeroNorm         1e-10  Angstroms
//     Skip direction vectors too short to normalise safely.
//     PhysicalConstants.h (NEAR_ZERO_NORM), used in 7+ files.
// ============================================================================

class ValueGate : public CalculationArea {
public:
    ValueGate(const char* name, const char* description,
              double floor, const char* unit);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The minimum value below which contributions are skipped.
    double Floor() const;

    // Unit of the gated quantity.
    const char* Unit() const;

private:
    const char* name_;
    const char* description_;
    double floor_;
    const char* unit_;
};


// ============================================================================
// SentinelValue
//
// A marker value meaning "no data." Any real computed value replaces
// it. Used to initialise nearest-distance fields before the first
// real measurement.
//
// Catalogue instance (1):
//   NoDataMarker  99.0 A
//     Used for: nearest_CO_dist, nearest_CN_dist, hbond_nearest_dist,
//     demo_nearest_ring_distance. Convention, not physics.
// ============================================================================

class SentinelValue : public CalculationArea {
public:
    SentinelValue(const char* name, const char* description, double value);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

    // The sentinel value.
    double Value() const;

private:
    const char* name_;
    const char* description_;
    double value_;
};


// ============================================================================
// WholeDomain
//
// Documents the deliberate ABSENCE of a spatial cutoff. The Coulomb
// sum runs over all N atoms with no truncation because the 1/r
// potential has no convergent truncation radius — unlike 1/r^3
// (dipolar) or 1/r^6 (dispersion) which decay fast enough to
// truncate.
//
// This is a CalculationArea because "no cutoff" is a design decision
// that needs documentation and provenance, just like "cutoff at 15A."
//
// Catalogue instance (1):
//   CoulombWholeProtein  domain=infinite  Coulomb/MopacCoulomb
//   Reference: first principles (Coulomb 1/r decay)
// ============================================================================

class WholeDomain : public CalculationArea {
public:
    WholeDomain(const char* name, const char* description);

    const char* Name() const override;
    const char* Description() const override;
    DomainKind Domain() const override;

private:
    const char* name_;
    const char* description_;
};


// SpatialIndexHorizon (15.0 A) is a RadialThreshold instance, not its
// own type. It must be >= the maximum calculator cutoff (currently
// RingHorizon at 15A). Its Description() documents this dependency.


}  // namespace nmr
