# GeometryChoice Agent Brief — Batch 2

## What you are doing

Adding GeometryChoice recording to calculator Compute() methods. This is
runtime documentation for UI visualization. You are NOT changing physics.
.npy output must be binary identical.

## Rules — absolute, no exceptions

1. **Everything inside the lambda.** The `choices.Record(...)` lambda is the
   containment boundary. No computation, no detection, no variable assignment
   outside the lambda that exists only to feed the lambda. If you need to know
   something, read it from objects already in scope inside the lambda.

2. **Use the filter set.** Every exclusion goes through `filters.AcceptAll(ctx)`.
   MinDistanceFilter is now in the filter set — do NOT write independent
   `if (distance < MIN_DISTANCE) continue;` checks. Add MinDistanceFilter to
   the filter set instead. After AcceptAll returns false, call
   `filters.LastRejectorName()` inside the lambda to identify which filter.

3. **Record exclusions, not inclusions.** Don't record every atom that passes.
   Record the ones that get rejected, and why.

4. **Per-source facts are per-source.** Ring current intensity, lobe offset,
   etc. are properties of the ring, not the atom-ring pair. Record once per
   ring using `std::set<size_t> recorded_rings` to deduplicate.

5. **Per-atom aggregates are per-atom.** Shell counts, nearest distances — 
   record after the inner loop, one choice per atom.

6. **Triggered events are interesting.** E-field clamp, adaptive refinement —
   record only when they fire.

7. **Charge noise floor and switching function noise stay independent.** These
   check source properties (charge magnitude, switching function value), not
   geometry. They are not in the filter set. Record them only if they fire,
   inside a lambda at the point of the `continue`.

## GeometryChoice API

```cpp
#pragma once
//
// GeometryChoice: runtime record of one geometric decision made by a calculator.
//
// A GeometryChoice is a bag of model objects — ConformationAtom*, Ring*, Bond*,
// and named numbers — recording what a calculator decided, which entities were
// involved, and whether each was included or excluded.
//
// The conformation owns a flat vector<GeometryChoice>. Calculators populate it
// during Compute() via the GeometryChoiceBuilder factory, which enforces that
// every choice is populated within a lambda (keeping the recording code visually
// separate from the physics).
//
// The UI walks the list, follows pointers back into the live model, and draws.
//
// Entity roles:
//   Source  — the ring/bond/atom generating the field
//   Target  — the atom being shielded or evaluated
//   Context — referenced but not directly source or target (e.g. bond endpoints)
//
// Outcome:
//   Included     — entity passed the gate
//   Excluded     — entity failed the gate
//   Triggered    — a guard or clamp fired (singularity, sanity, refinement)
//   NotTriggered — a guard was checked but did not fire
//

#include "Types.h"
#include <string>
#include <vector>
#include <functional>

namespace nmr {

// Forward declarations — these are the live model objects.
class Ring;
class ConformationAtom;
struct Bond;
class ProteinConformation;


// ============================================================================
// Enums
// ============================================================================

enum class EntityRole   { Source, Target, Context };
enum class EntityOutcome { Included, Excluded, Triggered, NotTriggered };


// ============================================================================
// GeometryEntity — one entry in the bag.
//
// Exactly one of {atom, ring, bond} is non-null per entry.
// Named numbers use the NamedNumber struct instead.
// ============================================================================

struct GeometryEntity {
    const ConformationAtom* atom = nullptr;
    const Ring*             ring = nullptr;
    const Bond*             bond = nullptr;
    size_t                  atom_index = SIZE_MAX;  // into protein atom list

    EntityRole    role    = EntityRole::Target;
    EntityOutcome outcome = EntityOutcome::Included;

    // Optional: which filter rejected this entity (empty if not rejected)
    std::string   filter_name;
};


// ============================================================================
// NamedNumber — a numeric value with a name and unit.
//
// "horizon", 15.0, "A"
// "intensity", -12.0, "nA"
// "distance", 7.3, "A"
// ============================================================================

struct NamedNumber {
    std::string name;
    double      value = 0.0;
    std::string unit;
};


// ============================================================================
// GeometryChoice — the bag itself.
// ============================================================================

class GeometryChoice {
    friend class GeometryChoiceBuilder;
    friend void AddAtom(GeometryChoice&, const ConformationAtom*, size_t,
                        EntityRole, EntityOutcome, const char*);
    friend void AddRing(GeometryChoice&, const Ring*,
                        EntityRole, EntityOutcome, const char*);
    friend void AddBond(GeometryChoice&, const Bond*,
                        EntityRole, EntityOutcome, const char*);
    friend void AddNumber(GeometryChoice&, const char*, double, const char*);
    friend void SetSampler(GeometryChoice&, std::function<SphericalTensor(Vec3)>);

public:
    // Read-only access for the UI
    const std::string&              Label()      const { return label_; }
    CalculatorId                    Calculator() const { return calculator_; }
    size_t                          GroupKey()   const { return group_key_; }
    const std::vector<GeometryEntity>& Entities() const { return entities_; }
    const std::vector<NamedNumber>&    Numbers()  const { return numbers_; }

    // Optional field sampler: evaluates this choice's physics at any 3D point.
    // Captures the source geometry (ring vertices, bond midpoint, etc.) so
    // the UI can draw field lines, isosurfaces, or probe values interactively.
    // Returns SphericalTensor at the given point. Null if not applicable.
    bool HasSampler() const { return sampler_ != nullptr; }
    SphericalTensor SampleAt(Vec3 point) const {
        return sampler_ ? sampler_(point) : SphericalTensor{};
    }

private:
    GeometryChoice() = default;

    std::string label_;
    CalculatorId calculator_ = CalculatorId::BiotSavart;
    size_t group_key_ = 0;
    std::vector<GeometryEntity> entities_;
    std::vector<NamedNumber>    numbers_;
    std::function<SphericalTensor(Vec3)> sampler_;
};


// ============================================================================
// GeometryChoiceBuilder — factory that enforces population in a lambda.
//
// Usage in a calculator's Compute():
//
//   GeometryChoiceBuilder choices(conf);
//
//   choices.Record(CalculatorId::BiotSavart, group, "ring horizon",
//       [&](GeometryChoice& gc) {
//           AddRing(gc, ring, EntityRole::Source, EntityOutcome::Included);
//           AddAtom(gc, &conf.AtomAt(ai), ai, EntityRole::Target, EntityOutcome::Included);
//           AddNumber(gc, "horizon", 15.0, "A");
//           AddNumber(gc, "distance", dist, "A");
//       });
//
// The lambda body is the ONLY place entities and numbers are added.
// This keeps recording code visually offset from physics code.
// ============================================================================

class GeometryChoiceBuilder {
public:
    explicit GeometryChoiceBuilder(ProteinConformation& conf);

    // Record one geometric decision. The populate lambda fills the choice,
    // then it goes straight onto the conformation's list. No commit step.
    void Record(CalculatorId calculator,
                size_t group_key,
                const char* label,
                std::function<void(GeometryChoice&)> populate);

private:
    ProteinConformation& conf_;
};


// ============================================================================
// Convenience methods on GeometryChoice for use inside the populate lambda.
// ============================================================================

inline void AddAtom(GeometryChoice& gc,
                    const ConformationAtom* atom,
                    size_t atom_index,
                    EntityRole role,
                    EntityOutcome outcome,
                    const char* filter = nullptr)
{
    GeometryEntity e;
    e.atom = atom;
    e.atom_index = atom_index;
    e.role = role;
    e.outcome = outcome;
    if (filter) e.filter_name = filter;
    gc.entities_.push_back(std::move(e));
}

inline void AddRing(GeometryChoice& gc,
                    const Ring* ring,
                    EntityRole role,
                    EntityOutcome outcome,
                    const char* filter = nullptr)
{
    GeometryEntity e;
    e.ring = ring;
    e.role = role;
    e.outcome = outcome;
    if (filter) e.filter_name = filter;
    gc.entities_.push_back(std::move(e));
}

inline void AddBond(GeometryChoice& gc,
                    const Bond* bond,
                    EntityRole role,
                    EntityOutcome outcome,
                    const char* filter = nullptr)
{
    GeometryEntity e;
    e.bond = bond;
    e.role = role;
    e.outcome = outcome;
    if (filter) e.filter_name = filter;
    gc.entities_.push_back(std::move(e));
}

inline void AddNumber(GeometryChoice& gc,
                      const char* name,
                      double value,
                      const char* unit)
{
    gc.numbers_.push_back({name, value, unit});
}

inline void SetSampler(GeometryChoice& gc,
                       std::function<SphericalTensor(Vec3)> sampler)
{
    gc.sampler_ = std::move(sampler);
}


}  // namespace nmr
```

## KernelEvaluationFilter (with MinDistanceFilter and LastRejectorName)

```cpp
#pragma once
//
// KernelEvaluationFilter: ABC and concrete filters for deciding whether
// a geometric kernel should be evaluated at a given atom-source pair.
//
// Every calculator that evaluates a geometric kernel (dipolar, Coulomb,
// surface integral) holds a KernelFilterSet. Before computing the kernel
// at a field point, the calculator builds a KernelEvaluationContext from
// the evaluation geometry and passes it to the filter set. All filters
// must accept for the kernel to be evaluated.
//
// Filters operate on first-principles physics criteria:
//   - DipolarNearFieldFilter: multipole expansion invalid inside source
//   - SelfSourceFilter: field undefined at the source itself
//   - SequentialExclusionFilter: through-bond coupling not modelled
//
// Type-based filters (SelfSourceFilter) are permitted when the physics
// reason is explicit: an atom that IS the source cannot be a field point
// for that source's dipolar field.
//
// Filter sets are configurable from TOML. The test for a filter is:
// it works but not as well without it. A filter that changes nothing
// is dead code. A filter without which the calculator crashes is
// papering over a bug.
//

#include "PhysicalConstants.h"
#include <vector>
#include <set>
#include <memory>
#include <string>
#include <cstddef>

namespace nmr {

class Protein;
class ProteinConformation;


// ============================================================================
// KernelEvaluationContext: the geometric facts about one potential
// kernel evaluation. Assembled by the calculator from already-computed
// distances and topology. No new spatial computation happens here.
// ============================================================================

struct KernelEvaluationContext {
    // Distance from field point to source expansion center (Angstroms).
    // Already computed as part of the kernel setup.
    double distance = 0.0;

    // Spatial extent of the source distribution (Angstroms).
    // Bond: bond length. H-bond: N...O distance. Ring: ring diameter.
    // Zero if the source is a point (single atom in Coulomb sum).
    double source_extent = 0.0;

    // Sequence separation between field-point residue and source residue.
    // -1 if not applicable or not known.
    int sequence_separation = -1;

    // Field-point atom index.
    size_t atom_index = SIZE_MAX;

    // Source atom indices. For a bond: the two endpoints. For an H-bond:
    // donor N and acceptor O. For a ring: SIZE_MAX (no single atom).
    // For a Coulomb pair: source_atom_a = the source atom, source_atom_b = SIZE_MAX.
    size_t source_atom_a = SIZE_MAX;
    size_t source_atom_b = SIZE_MAX;

    // Source ring index in the Protein's ring list. SIZE_MAX if the
    // source is not a ring (bonds, H-bonds, Coulomb pairs).
    // Used by RingBondedExclusionFilter to check whether the field-point
    // atom is a vertex of this ring or bonded to one.
    size_t source_ring_index = SIZE_MAX;

    // Optional: the conformation that owns the evaluation geometry.
    // Available for future physics-aware filters that need to query
    // positions, bond directions, or already-computed results at
    // neighbouring atoms. Do NOT use this for simple distance/index
    // checks — use the fields above. This exists so that the mechanism
    // for richer filters stays within the object model rather than
    // requiring a parallel data path.
    //
    // This is safe to use: calculators run in dependency order and
    // all declared dependencies are attached before Compute() runs.
    // Do NOT add order-checking machinery to this library. The
    // dependency graph already guarantees safety.
    const ProteinConformation* conf = nullptr;
};


// ============================================================================
// ABC: one named physics criterion for kernel evaluation.
// ============================================================================

class KernelEvaluationFilter {
public:
    virtual ~KernelEvaluationFilter() = default;

    // Should this kernel evaluation proceed?
    virtual bool Accept(const KernelEvaluationContext& ctx) const = 0;

    // Short identifier for logging and TOML configuration.
    virtual const char* Name() const = 0;

    // Physics reason this filter exists.
    virtual const char* Description() const = 0;

    // Human-readable rejection reason with specific values from ctx.
    // Called only when Accept() returns false and logging is enabled.
    // Default: name + distance. Override for filter-specific detail.
    virtual std::string RejectReason(const KernelEvaluationContext& ctx) const {
        return std::string(Name()) + ": atom=" + std::to_string(ctx.atom_index)
            + " dist=" + std::to_string(ctx.distance);
    }
};


// ============================================================================
// DipolarNearFieldFilter
//
// Physics: the dipolar kernel K_ab = (3 d̂_a d̂_b - δ_ab) / r³ models
// the source as a point at the expansion center. When the field point
// is closer to the center than the source radius (half the source
// extent), the multipole expansion does not converge. The kernel value
// is large but physically meaningless — the field point is inside the
// charge/current distribution.
//
// Criterion: distance > source_extent / 2.
// Source extent is a measured property of each source, not a global
// tuning parameter.
//
// Demonstrated effect (HBondResult, 466 proteins):
//   With filter:  max |T2| = 0.78 A^-3
//   Without:      max |T2| = 1908 A^-3  (atoms 0.2A from H-bond midpoint)
// ============================================================================

// ============================================================================
// MinDistanceFilter
//
// Physics: at distances below 0.1 A the field point is essentially on top
// of the source. The 1/r^n kernels diverge numerically. The kernel helpers
// (ComputeBondKernel, AccumulateTensor, etc.) also guard internally, but
// this filter prevents the evaluation from starting at all and provides
// a named rejection for GeometryChoice recording.
//
// Criterion: distance >= MIN_DISTANCE (0.1 A).
// ============================================================================

class MinDistanceFilter : public KernelEvaluationFilter {
public:
    bool Accept(const KernelEvaluationContext& ctx) const override {
        return ctx.distance >= MIN_DISTANCE;
    }

    const char* Name() const override {
        return "MinDistanceFilter";
    }

    const char* Description() const override {
        return "Singularity guard: 1/r^n kernels diverge below 0.1 A";
    }

    std::string RejectReason(const KernelEvaluationContext& ctx) const override {
        return std::string("MinDistanceFilter: atom=")
            + std::to_string(ctx.atom_index)
            + " distance=" + std::to_string(ctx.distance)
            + "A < " + std::to_string(MIN_DISTANCE) + "A";
    }
};


class DipolarNearFieldFilter : public KernelEvaluationFilter {
public:
    bool Accept(const KernelEvaluationContext& ctx) const override {
        if (ctx.source_extent <= 0.0) return true;  // point source, no near field
        return ctx.distance > ctx.source_extent * 0.5;
    }

    const char* Name() const override {
        return "DipolarNearFieldFilter";
    }

    const char* Description() const override {
        return "Multipole expansion invalid when field point is within "
               "half the source extent (inside the source distribution)";
    }

    std::string RejectReason(const KernelEvaluationContext& ctx) const override {
        return std::string("DipolarNearFieldFilter: atom=")
            + std::to_string(ctx.atom_index)
            + " distance=" + std::to_string(ctx.distance)
            + "A < threshold=" + std::to_string(ctx.source_extent * 0.5)
            + "A (source_extent=" + std::to_string(ctx.source_extent) + "A)";
    }
};


// ============================================================================
// SelfSourceFilter
//
// Physics: the dipolar field of a source is undefined at the source
// itself. An atom that IS an endpoint of a bond cannot be a field point
// for that bond's McConnell kernel. An atom that IS the donor N or
// acceptor O cannot be a field point for that H-bond's kernel.
//
// This is type-based: it checks atom indices, not distances. The physics
// justification is that the source atom's own electron cloud generates
// the field — it does not experience its own dipolar field as an
// external perturbation.
//
// Criterion: atom_index != source_atom_a AND atom_index != source_atom_b.
// ============================================================================

class SelfSourceFilter : public KernelEvaluationFilter {
public:
    bool Accept(const KernelEvaluationContext& ctx) const override {
        if (ctx.atom_index == ctx.source_atom_a) return false;
        if (ctx.source_atom_b != SIZE_MAX && ctx.atom_index == ctx.source_atom_b)
            return false;
        return true;
    }

    const char* Name() const override {
        return "SelfSourceFilter";
    }

    const char* Description() const override {
        return "Dipolar field undefined at the source atom itself "
               "(atom cannot be field point for its own source)";
    }

    std::string RejectReason(const KernelEvaluationContext& ctx) const override {
        size_t match = (ctx.atom_index == ctx.source_atom_a)
            ? ctx.source_atom_a : ctx.source_atom_b;
        return std::string("SelfSourceFilter: atom=")
            + std::to_string(ctx.atom_index)
            + " IS source_atom=" + std::to_string(match);
    }
};


// ============================================================================
// SequentialExclusionFilter
//
// Physics: atoms within a few residues in sequence share through-bond
// electronic coupling (inductive, hyperconjugation) that the through-
// space dipolar model does not capture. The dipolar kernel gives a
// geometric prediction at these distances, but the actual shielding
// perturbation is dominated by through-bond effects that have a
// different angular structure.
//
// Criterion: sequence_separation >= min_separation.
// The default min_separation (2) excludes same-residue and immediate
// neighbour interactions. Configurable from TOML.
// ============================================================================

class SequentialExclusionFilter : public KernelEvaluationFilter {
public:
    explicit SequentialExclusionFilter(int min_separation = 2)
        : min_separation_(min_separation) {}

    bool Accept(const KernelEvaluationContext& ctx) const override {
        if (ctx.sequence_separation < 0) return true;  // unknown, allow
        return ctx.sequence_separation >= min_separation_;
    }

    const char* Name() const override {
        return "SequentialExclusionFilter";
    }

    const char* Description() const override {
        return "Through-bond coupling dominates through-space dipolar "
               "at small sequence separations (dipolar model inapplicable)";
    }

    std::string RejectReason(const KernelEvaluationContext& ctx) const override {
        return std::string("SequentialExclusionFilter: atom=")
            + std::to_string(ctx.atom_index)
            + " seq_sep=" + std::to_string(ctx.sequence_separation)
            + " < min=" + std::to_string(min_separation_);
    }

    int MinSeparation() const { return min_separation_; }

private:
    int min_separation_;
};


// ============================================================================
// RingBondedExclusionFilter
//
// Physics: an atom that is a vertex of a ring, or covalently bonded to a
// ring vertex, is in the through-bond regime for that ring's field. The
// through-space multipole expansion (dipolar, quadrupole, surface integral,
// dispersion) does not model through-bond electronic coupling. Evaluating
// the kernel at these atoms gives a large number that is physically
// meaningless — the atom is part of the source distribution.
//
// This is a topological check, not a distance proxy. It uses the bond
// graph to determine membership. Compare to DipolarNearFieldFilter which
// catches most of the same atoms by distance, but has a boundary case
// for ring atoms that sit exactly at the filter threshold.
//
// Construction: takes a const reference to the Protein whose topology
// defines which atoms are ring vertices and what they are bonded to.
// At construction, walks Protein → Ring → atom_indices → bond_indices
// to build a per-ring exclusion set. This is the one-way path from
// the formal protein definition down through typed topology.
//
// Criterion: atom_index NOT in excluded set for source_ring_index.
// Requires source_ring_index to be set on the evaluation context.
// ============================================================================

class RingBondedExclusionFilter : public KernelEvaluationFilter {
public:
    // Walks the Protein's topology once to build per-ring exclusion sets.
    // Each set contains the ring's vertices and all atoms bonded to them.
    explicit RingBondedExclusionFilter(const Protein& protein);

    bool Accept(const KernelEvaluationContext& ctx) const override {
        if (ctx.source_ring_index == SIZE_MAX) return true;  // not a ring source
        if (ctx.source_ring_index >= ring_bonded_.size()) return true;
        return ring_bonded_[ctx.source_ring_index].count(ctx.atom_index) == 0;
    }

    const char* Name() const override {
        return "RingBondedExclusionFilter";
    }

    const char* Description() const override {
        return "Through-space kernel invalid for ring vertices and their "
               "bonded neighbours (through-bond regime, not through-space)";
    }

    std::string RejectReason(const KernelEvaluationContext& ctx) const override {
        return std::string("RingBondedExclusionFilter: atom=")
            + std::to_string(ctx.atom_index)
            + " is bonded to ring " + std::to_string(ctx.source_ring_index);
    }

private:
    std::vector<std::set<size_t>> ring_bonded_;
};


// ============================================================================
// KernelFilterSet: a calculator's ordered list of filters.
// All must accept for the kernel to be evaluated.
// ============================================================================

class KernelFilterSet {
public:
    void Add(std::unique_ptr<KernelEvaluationFilter> filter) {
        filters_.push_back(std::move(filter));
        rejection_counts_.push_back(0);
    }

    // Returns true if all filters accept. Tracks per-filter rejection counts
    // and logs each rejection with the evaluation context.
    bool AcceptAll(const KernelEvaluationContext& ctx) {
        for (size_t i = 0; i < filters_.size(); ++i) {
            if (!filters_[i]->Accept(ctx)) {
                rejection_counts_[i]++;
                last_rejector_ = filters_[i]->Name();
                if (log_rejections_) {
                    rejections_.push_back({
                        filters_[i]->Name(),
                        ctx.atom_index,
                        ctx.source_atom_a,
                        ctx.source_atom_b,
                        ctx.distance,
                        ctx.source_extent
                    });
                    rejection_reasons_.push_back(
                        filters_[i]->RejectReason(ctx));
                }
                return false;
            }
        }
        return true;
    }

    // Enable/disable per-rejection logging. Off by default for batch
    // performance. Enable for single-protein diagnostics.
    void SetLogRejections(bool enable) { log_rejections_ = enable; }

    // Number of filters in the set.
    size_t Size() const { return filters_.size(); }

    // Name of the filter that rejected the most recent AcceptAll() call.
    // Only meaningful after AcceptAll() returns false.
    const char* LastRejectorName() const { return last_rejector_; }

    // Access for logging/reporting.
    const KernelEvaluationFilter& FilterAt(size_t i) const {
        return *filters_[i];
    }

    // How many evaluations this filter rejected.
    int RejectionCount(size_t i) const {
        return rejection_counts_[i];
    }

    // Total rejections across all filters.
    int TotalRejections() const {
        int total = 0;
        for (int c : rejection_counts_) total += c;
        return total;
    }

    // Report all filter names (for startup logging).
    std::string Describe() const {
        std::string s;
        for (size_t i = 0; i < filters_.size(); ++i) {
            if (i > 0) s += ", ";
            s += filters_[i]->Name();
        }
        return s;
    }

    // Report rejection counts per filter (for completion logging).
    std::string ReportRejections() const {
        std::string s;
        for (size_t i = 0; i < filters_.size(); ++i) {
            if (i > 0) s += ", ";
            s += filters_[i]->Name();
            s += "=" + std::to_string(rejection_counts_[i]);
        }
        return s;
    }

    // One logged rejection: what was rejected and why.
    struct RejectionRecord {
        const char* filter_name;
        size_t atom_index;
        size_t source_atom_a;
        size_t source_atom_b;
        double distance;
        double source_extent;
    };

    const std::vector<RejectionRecord>& Rejections() const {
        return rejections_;
    }

    // Human-readable rejection reasons (what, where, which).
    // Only populated when log_rejections_ is enabled.
    const std::vector<std::string>& RejectionReasons() const {
        return rejection_reasons_;
    }

private:
    std::vector<std::unique_ptr<KernelEvaluationFilter>> filters_;
    std::vector<int> rejection_counts_;
    std::vector<RejectionRecord> rejections_;
    std::vector<std::string> rejection_reasons_;
    const char* last_rejector_ = nullptr;
    bool log_rejections_ = false;
};


}  // namespace nmr
```

## PhysicalConstants

```cpp
#pragma once
//
// Physical constants for NMR field calculations.
//
// Only universal constants here -- no model-specific parameters.
// Model parameters (lobe offsets, ring current intensities) belong
// on ring type classes and calculator parameter structs.
//

#include <cmath>

namespace nmr {

// Mathematical
constexpr double PI = 3.14159265358979323846;

// Electromagnetic (SI)
constexpr double VACUUM_PERMEABILITY = 1.25663706212e-6;   // T*m/A (mu_0)

// Unit conversions
constexpr double ANGSTROMS_TO_METRES = 1.0e-10;
constexpr double NANOAMPERES_TO_AMPERES = 1.0e-9;
constexpr double PPM_FACTOR = 1.0e6;

// Electrostatics in {e, Angstrom, eV} units
// ke = e / (4 pi epsilon_0) = 14.3996 eV*A / e = 14.3996 V*A
// Converts E from e/A^2 (raw Coulomb sum) to V/A (physical E-field).
constexpr double COULOMB_KE = 14.3996;

// Thermal voltage at 298.15 K: kT/e = k_B * T / e = 0.025693 V
// Converts APBS potential/field from kT/e units to Volts.
constexpr double KT_OVER_E_298K = 0.025693;

// Biot-Savart prefactor: mu_0/(4*pi) in SI units (T*m/A)
constexpr double BIOT_SAVART_PREFACTOR = VACUUM_PERMEABILITY / (4.0 * PI);

// Constitution: numerical thresholds
constexpr double MIN_DISTANCE = 0.1;            // Angstroms -- singularity cutoff
constexpr double NO_DATA_SENTINEL = 99.0;       // sentinel for missing data
constexpr double NEAR_ZERO_NORM = 1e-10;        // near-zero vector norm
constexpr double NEAR_ZERO_FIELD = 1e-15;       // near-zero field magnitude

// Constitution: spatial shells for ring counting
constexpr double RING_COUNT_SHELL_1 = 3.0;      // Angstroms
constexpr double RING_COUNT_SHELL_2 = 5.0;
constexpr double RING_COUNT_SHELL_3 = 8.0;
constexpr double RING_COUNT_SHELL_4 = 12.0;

// Constitution: calculation cutoffs
constexpr double RING_CALC_CUTOFF = 15.0;       // Angstroms -- ring current cutoff
constexpr double EXP_DECAY_LENGTH = 4.0;        // Angstroms
constexpr double PACKING_RADIUS = 8.0;          // Angstroms -- for heavy atom count

// Constitution: H-bond thresholds
constexpr double HBOND_COUNT_RADIUS = 3.5;      // Angstroms
constexpr double HBOND_MAX_DIST = 50.0;         // Angstroms
constexpr double APBS_SANITY_LIMIT = 100.0;     // V/Angstrom

// Note: dispersion distance thresholds (R_SWITCH=4.3A, R_CUT=5.0A) are
// defined in DispersionResult.cpp with full physics documentation.
// They are NOT global constants because they are specific to the
// dispersion switching function (CHARMM form, Brooks et al. 1983).

// Constitution: sequence exclusion
constexpr int SEQUENTIAL_EXCLUSION_THRESHOLD = 2;

}  // namespace nmr
```

## Ring.h

```cpp
#pragma once
//
// Ring type class hierarchy.
//
// Ring types ARE classes with physics properties baked in.
// Each type provides const properties derived from its identity.
// Calculator code is ring-type-agnostic: ring.Intensity(), ring.JBLobeOffset().
//
// 8 types in 3 size categories:
//   SixMemberedRing: PheBenzeneRing, TyrPhenolRing, TrpBenzeneRing
//   FiveMemberedRing: TrpPyrroleRing, HisImidazoleRing, HidImidazoleRing, HieImidazoleRing
//   FusedRing: IndolePerimeterRing (TRP 9-atom)
//

#include "Types.h"
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <memory>

namespace nmr {

// ============================================================================
// Ring::Geometry -- conformation-dependent, computed by GeometryResult
// ============================================================================

struct RingGeometry {
    Vec3              center = Vec3::Zero();
    Vec3              normal = Vec3::Zero();
    double            radius = 0.0;
    std::vector<Vec3> vertices;
};

// ============================================================================
// Accumulated ring properties (set by ConformationResult post-pass updates)
// ============================================================================

struct RingAccumulated {
    Vec3 total_B_at_center = Vec3::Zero();
    double intensity_used = 0.0;
    double total_G_T0_diagnostic = 0.0;
    std::map<size_t, Vec3> mutual_B_from;
};

// ============================================================================
// Ring (base class)
// ============================================================================

class Ring {
public:
    // Structural identity (topology, set at construction)
    std::vector<size_t> atom_indices;
    RingTypeIndex       type_index = RingTypeIndex::PheBenzene;
    size_t              parent_residue_index = 0;
    int                 parent_residue_number = 0;
    size_t              fused_partner_index = SIZE_MAX;

    virtual ~Ring() = default;

    // Virtual const properties (overridden by each type class)
    virtual double Intensity() const = 0;
    virtual double LiteratureIntensity() const = 0;
    virtual double JBLobeOffset() const = 0;
    virtual int NitrogenCount() const = 0;
    virtual RingAromaticity Aromaticity() const = 0;
    virtual int RingSizeValue() const = 0;
    virtual const char* TypeName() const = 0;

    // Non-virtual queries
    bool IsFused() const { return fused_partner_index != SIZE_MAX; }
    int TypeIndexAsInt() const { return static_cast<int>(type_index); }

    // Compute geometry from positions (SVD normal)
    RingGeometry ComputeGeometry(const std::vector<Vec3>& positions) const;

    // Accumulated properties (set during extraction passes)
    RingAccumulated accumulated;
};


// ============================================================================
// Six-membered rings
// ============================================================================

class SixMemberedRing : public Ring {
public:
    int RingSizeValue() const override { return 6; }
};

class PheBenzeneRing : public SixMemberedRing {
public:
    PheBenzeneRing() { type_index = RingTypeIndex::PheBenzene; }
    double Intensity() const override { return -12.0; }
    double LiteratureIntensity() const override { return -12.0; }
    double JBLobeOffset() const override { return 0.64; }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "PHE"; }
};

class TyrPhenolRing : public SixMemberedRing {
public:
    TyrPhenolRing() { type_index = RingTypeIndex::TyrPhenol; }
    double Intensity() const override { return -11.28; }
    double LiteratureIntensity() const override { return -11.28; }
    double JBLobeOffset() const override { return 0.64; }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TYR"; }
};

class TrpBenzeneRing : public SixMemberedRing {
public:
    TrpBenzeneRing() { type_index = RingTypeIndex::TrpBenzene; }
    double Intensity() const override { return -12.48; }
    double LiteratureIntensity() const override { return -12.48; }
    double JBLobeOffset() const override { return 0.64; }
    int NitrogenCount() const override { return 0; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    const char* TypeName() const override { return "TRP6"; }
};


// ============================================================================
// Five-membered rings
// ============================================================================

class FiveMemberedRing : public Ring {
public:
    int RingSizeValue() const override { return 5; }
};

class TrpPyrroleRing : public FiveMemberedRing {
public:
    TrpPyrroleRing() { type_index = RingTypeIndex::TrpPyrrole; }
    double Intensity() const override { return -6.72; }
    double LiteratureIntensity() const override { return -6.72; }
    double JBLobeOffset() const override { return 0.52; }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Reduced; }
    const char* TypeName() const override { return "TRP5"; }
};

class HisImidazoleRing : public FiveMemberedRing {
public:
    HisImidazoleRing() { type_index = RingTypeIndex::HisImidazole; }
    double Intensity() const override { return -5.16; }
    double LiteratureIntensity() const override { return -5.16; }
    double JBLobeOffset() const override { return 0.50; }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIS"; }
};

class HidImidazoleRing : public FiveMemberedRing {
public:
    HidImidazoleRing() { type_index = RingTypeIndex::HidImidazole; }
    double Intensity() const override { return -5.16; }
    double LiteratureIntensity() const override { return -5.16; }
    double JBLobeOffset() const override { return 0.50; }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HID"; }
};

class HieImidazoleRing : public FiveMemberedRing {
public:
    HieImidazoleRing() { type_index = RingTypeIndex::HieImidazole; }
    double Intensity() const override { return -5.16; }
    double LiteratureIntensity() const override { return -5.16; }
    double JBLobeOffset() const override { return 0.50; }
    int NitrogenCount() const override { return 2; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Weak; }
    const char* TypeName() const override { return "HIE"; }
};


// ============================================================================
// Fused ring (TRP 9-atom indole perimeter)
// ============================================================================

class FusedRing : public Ring {};

class IndolePerimeterRing : public FusedRing {
public:
    IndolePerimeterRing() { type_index = RingTypeIndex::TrpPerimeter; }
    double Intensity() const override { return -19.2; }
    double LiteratureIntensity() const override { return -19.2; }
    double JBLobeOffset() const override { return 0.60; }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::Full; }
    int RingSizeValue() const override { return 9; }
    const char* TypeName() const override { return "TRP9"; }
};


// ============================================================================
// Factory: create a Ring subclass from a RingTypeIndex
// ============================================================================

std::unique_ptr<Ring> CreateRing(RingTypeIndex type);

}  // namespace nmr
```

## EXAMPLE: BiotSavartResult.cpp (already instrumented — study this pattern)

```cpp
#include "BiotSavartResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <set>

namespace nmr {


std::vector<std::type_index> BiotSavartResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// Wire segment B-field (Biot-Savart law).
//
// All computation in SI: positions in metres, current in amperes, B in Tesla.
//
//   B = (mu_0/4pi) * I * (dl x dA) / |dl x dA|^2 * (dl.dA/|dA| - dl.dB/|dB|)
//
// where:
//   dl = b - a  (wire segment vector)
//   dA = r - a  (field point from segment start)
//   dB = r - b  (field point from segment end)
//
// Numerical guards (in SI, so thresholds are small):
//   lenA, lenB < 1e-25 m:  field point at segment endpoint
//   crossSq < 1e-70 m^2:   field point on the wire axis
// ============================================================================

static Vec3 WireSegmentField(
        const Vec3& a_m, const Vec3& b_m,
        double I_A, const Vec3& r_m) {

    Vec3 dl_m = b_m - a_m;
    Vec3 dA_m = r_m - a_m;
    Vec3 dB_m = r_m - b_m;

    double lenA = dA_m.norm();
    double lenB = dB_m.norm();
    if (lenA < 1e-25 || lenB < 1e-25) return Vec3::Zero();

    Vec3 cross = dl_m.cross(dA_m);
    double crossSq = cross.squaredNorm();
    if (crossSq < 1e-70) return Vec3::Zero();

    // mu_0/(4*pi) = 1e-7 T*m/A (exact in SI)
    constexpr double mu0_4pi = 1e-7;
    double factor = mu0_4pi * I_A / crossSq;
    double dotTerm = dl_m.dot(dA_m) / lenA - dl_m.dot(dB_m) / lenB;

    return factor * dotTerm * cross;  // Tesla
}


// ============================================================================
// Johnson-Bovey double-loop model.
//
// Two current loops at +/- lobe_offset from the ring plane along the normal.
// Each loop carries half the total current (I/2). The total B-field is the
// sum over all wire segments of both loops.
//
// Input: vertex positions in Angstroms, current in nanoamperes.
// Converts to SI at the boundary, computes in pure SI, returns B in Tesla.
// ============================================================================

static Vec3 JohnsonBoveyField(
        const std::vector<Vec3>& vertices,
        const Vec3& normal,
        double lobe_offset_ang,
        double current_nanoamperes,
        const Vec3& point_ang) {

    int n = static_cast<int>(vertices.size());
    if (n < 3) return Vec3::Zero();

    // Unit conversion at the boundary: Angstroms -> metres, nA -> A.
    // After this block, all computation is pure SI.
    Vec3 offset_ang = normal * lobe_offset_ang;
    double halfI_A = 0.5 * current_nanoamperes * NANOAMPERES_TO_AMPERES;

    Vec3 B = Vec3::Zero();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;

        // Upper loop (z = +d)
        Vec3 a_upper = (vertices[i] + offset_ang) * ANGSTROMS_TO_METRES;
        Vec3 b_upper = (vertices[j] + offset_ang) * ANGSTROMS_TO_METRES;

        // Lower loop (z = -d)
        Vec3 a_lower = (vertices[i] - offset_ang) * ANGSTROMS_TO_METRES;
        Vec3 b_lower = (vertices[j] - offset_ang) * ANGSTROMS_TO_METRES;

        Vec3 r_m = point_ang * ANGSTROMS_TO_METRES;

        B += WireSegmentField(a_upper, b_upper, halfI_A, r_m);
        B += WireSegmentField(a_lower, b_lower, halfI_A, r_m);
    }

    return B;  // Tesla
}


// ============================================================================
// BiotSavartResult::Compute
//
// For each atom, find rings within RING_CALC_CUTOFF (15A). For each ring:
//   1. Compute B-field from JB double-loop model (unit current, I=1 nA)
//   2. Construct geometric kernel G_ab = n_b * B_a * PPM_FACTOR
//   3. Store G, SphericalTensor(G), B, cylindrical coords on RingNeighbourhood
//   4. Accumulate per-type T0 and T2 sums on ConformationAtom
// ============================================================================

std::unique_ptr<BiotSavartResult> BiotSavartResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("BiotSavartResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " rings=" + std::to_string(conf.ProteinRef().RingCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_rings = protein.RingCount();

    auto result_ptr = std::make_unique<BiotSavartResult>();
    result_ptr->conf_ = &conf;

    if (n_rings == 0) {
        OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
            "no rings — nothing to compute");
        return result_ptr;
    }

    // Filter set: DipolarNearFieldFilter with source_extent = ring diameter,
    // plus RingBondedExclusionFilter for topological exclusion of ring
    // vertices and their bonded neighbours.
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
        "filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);
    std::set<size_t> recorded_rings;

    int total_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, RING_CALC_CUTOFF);

        Mat3 G_total = Mat3::Zero();
        Vec3 B_total = Vec3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            if (geom.vertices.size() < 3) continue;

            // ---- GeometryChoice: ring current ----
            if (recorded_rings.find(ri) == recorded_rings.end()) {
                recorded_rings.insert(ri);
                auto verts_copy = geom.vertices;
                Vec3 normal_copy = geom.normal;
                double lobe_copy = ring.JBLobeOffset();
                choices.Record(CalculatorId::BiotSavart, ri, "ring current",
                    [&ring, verts_copy, normal_copy, lobe_copy](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddNumber(gc, "intensity", ring.Intensity(), "nA");
                        AddNumber(gc, "lobe_offset", ring.JBLobeOffset(), "A");
                        SetSampler(gc, [verts_copy, normal_copy, lobe_copy](Vec3 pt) -> SphericalTensor {
                            Vec3 B = JohnsonBoveyField(verts_copy, normal_copy, lobe_copy, 1.0, pt);
                            Mat3 G;
                            for (int a = 0; a < 3; ++a)
                                for (int b = 0; b < 3; ++b)
                                    G(a, b) = -normal_copy(b) * B(a) * PPM_FACTOR;
                            return SphericalTensor::Decompose(G);
                        });
                    });
            }

            double distance = (atom_pos - geom.center).norm();

            // Apply filter: source extent = ring diameter (2 * radius)
            KernelEvaluationContext ctx;
            ctx.distance = distance;
            ctx.source_extent = 2.0 * geom.radius;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) {
                // ---- GeometryChoice: near-field exclusion ----
                choices.Record(CalculatorId::BiotSavart, ri, "near-field exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                continue;
            }

            // B-field from JB model with unit current (1.0 nA).
            // The geometric kernel is independent of intensity.
            Vec3 B = JohnsonBoveyField(
                geom.vertices, geom.normal,
                ring.JBLobeOffset(), 1.0, atom_pos);

            // Geometric kernel: G_ab = -n_b * B_a * PPM_FACTOR
            // The minus sign comes from the shielding tensor definition:
            //   sigma_ab = -dB_a^sec / dB_{0,b}
            // With this convention, sigma = I * G gives the correct sign
            // using literature I values (I < 0 for diamagnetic rings).
            // Verified: I=-12, G_T0=-0.116 at (0,0,3A) above PHE
            //   -> sigma = (-12)(-0.116) = +1.40 ppm (shielded). Correct.
            Mat3 G;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    G(a, b) = -geom.normal(b) * B(a) * PPM_FACTOR;

            // Find or create RingNeighbourhood for this ring
            RingNeighbourhood* rn = nullptr;
            for (auto& existing : ca.ring_neighbours) {
                if (existing.ring_index == ri) {
                    rn = &existing;
                    break;
                }
            }
            if (!rn) {
                RingNeighbourhood new_rn;
                new_rn.ring_index = ri;
                new_rn.ring_type = ring.type_index;
                new_rn.distance_to_center = distance;

                Vec3 d = atom_pos - geom.center;
                new_rn.direction_to_center = d.normalized();

                // Cylindrical coordinates in ring frame
                double z = d.dot(geom.normal);
                Vec3 d_plane = d - z * geom.normal;
                double rho = d_plane.norm();
                double theta = std::atan2(d_plane.norm(), std::abs(z));
                new_rn.z = z;
                new_rn.rho = rho;
                new_rn.theta = theta;

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            // Store BS results on RingNeighbourhood
            rn->G_tensor = G;
            rn->G_spherical = SphericalTensor::Decompose(G);
            rn->B_field = B;

            // B-field in cylindrical coordinates (ring frame)
            Vec3 d = atom_pos - geom.center;
            double z_coord = d.dot(geom.normal);
            Vec3 d_plane = d - z_coord * geom.normal;
            double rho_mag = d_plane.norm();
            Vec3 rho_hat = Vec3::Zero();
            if (rho_mag > NEAR_ZERO_NORM) rho_hat = d_plane / rho_mag;
            rn->B_cylindrical = Vec3(
                B.dot(rho_hat),        // B_rho
                0.0,                   // B_phi (zero by axial symmetry choice)
                B.dot(geom.normal));   // B_z

            // Accumulate totals
            G_total += G;
            B_total += B;

            // Per-type T0 and T2 sums
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < 8) {
                ca.per_type_G_T0_sum[ti] += rn->G_spherical.T0;
                for (int c = 0; c < 5; ++c)
                    ca.per_type_G_T2_sum[ti][c] += rn->G_spherical.T2[c];
            }

            total_pairs++;
        }

        // Store accumulated totals on ConformationAtom
        ca.total_B_field += B_total;
        ca.total_G_tensor += G_total;
        ca.total_G_spherical = SphericalTensor::Decompose(
            ca.total_G_tensor);
        ca.bs_shielding_contribution = SphericalTensor::Decompose(G_total);

        // Ring proximity counts
        for (const auto& rn : ca.ring_neighbours) {
            if (rn.distance_to_center <= RING_COUNT_SHELL_1) ca.n_rings_within_3A++;
            if (rn.distance_to_center <= RING_COUNT_SHELL_2) ca.n_rings_within_5A++;
            if (rn.distance_to_center <= RING_COUNT_SHELL_3) ca.n_rings_within_8A++;
            if (rn.distance_to_center <= RING_COUNT_SHELL_4) ca.n_rings_within_12A++;
        }

        // ---- GeometryChoice: ring shells ----
        choices.Record(CalculatorId::BiotSavart, ai, "ring shells",
            [&ca, ai](GeometryChoice& gc) {
                AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Included);
                AddNumber(gc, "n_within_3A", static_cast<double>(ca.n_rings_within_3A), "count");
                AddNumber(gc, "n_within_5A", static_cast<double>(ca.n_rings_within_5A), "count");
                AddNumber(gc, "n_within_8A", static_cast<double>(ca.n_rings_within_8A), "count");
                AddNumber(gc, "n_within_12A", static_cast<double>(ca.n_rings_within_12A), "count");
            });
    }

    OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


// ============================================================================
// SampleBFieldAt / SampleShieldingAt: evaluate at arbitrary 3D points.
//
// Same physics as Compute(), but for a single point. No filter set
// (grid points are not atoms — no self-exclusion or bonded-exclusion
// applies). DipolarNearFieldFilter still applied for physical validity.
// ============================================================================

Vec3 BiotSavartResult::SampleBFieldAt(Vec3 point) const {
    if (!conf_) return Vec3::Zero();

    const Protein& protein = conf_->ProteinRef();
    const size_t n_rings = protein.RingCount();

    Vec3 B_total = Vec3::Zero();

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        double distance = (point - geom.center).norm();
        if (distance < MIN_DISTANCE) continue;

        // DipolarNearFieldFilter: multipole invalid inside source
        if (distance < geom.radius) continue;

        // Distance cutoff
        if (distance > RING_CALC_CUTOFF) continue;

        B_total += JohnsonBoveyField(
            geom.vertices, geom.normal,
            ring.JBLobeOffset(), 1.0, point);
    }

    return B_total;
}

SphericalTensor BiotSavartResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    const size_t n_rings = protein.RingCount();

    Mat3 G_total = Mat3::Zero();

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        double distance = (point - geom.center).norm();
        if (distance < MIN_DISTANCE) continue;
        if (distance < geom.radius) continue;
        if (distance > RING_CALC_CUTOFF) continue;

        Vec3 B = JohnsonBoveyField(
            geom.vertices, geom.normal,
            ring.JBLobeOffset(), 1.0, point);

        // G_ab = -n_b * B_a * PPM_FACTOR
        Mat3 G;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                G(a, b) = -geom.normal(b) * B(a) * PPM_FACTOR;

        G_total += G;
    }

    return SphericalTensor::Decompose(G_total);
}


// ============================================================================
// WriteFeatures: export what Compute() wrote on ConformationAtom.
//
// Every field this method reads was written by Compute() above. The arrays
// match Compute()'s accumulation: shielding contribution (the full
// SphericalTensor of the summed G over all rings), per-type T0 and T2
// sums (8 ring types, matching the RingTypeIndex enum), ring proximity
// counts, and the total B-field vector.
//
// Pack order for SphericalTensor: [T0, T1[0..2], T2[0..4]] = 9 doubles.
// ============================================================================

static void PackST(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int BiotSavartResult::WriteFeatures(const ProteinConformation& conf,
                                     const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // bs_shielding: (N, 9) — the full SphericalTensor sum over all rings
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i)
            PackST(conf.AtomAt(i).bs_shielding_contribution, &data[i*9]);
        NpyWriter::WriteFloat64(output_dir + "/bs_shielding.npy", data.data(), N, 9);
        written++;
    }

    // bs_per_type_T0: (N, 8) — isotropic kernel per ring type
    {
        std::vector<double> data(N * 8);
        for (size_t i = 0; i < N; ++i)
            for (int t = 0; t < 8; ++t)
                data[i*8 + t] = conf.AtomAt(i).per_type_G_T0_sum[t];
        NpyWriter::WriteFloat64(output_dir + "/bs_per_type_T0.npy", data.data(), N, 8);
        written++;
    }

    // bs_per_type_T2: (N, 40) — T2[5] per ring type[8]
    {
        std::vector<double> data(N * 40);
        for (size_t i = 0; i < N; ++i)
            for (int t = 0; t < 8; ++t)
                for (int c = 0; c < 5; ++c)
                    data[i*40 + t*5 + c] = conf.AtomAt(i).per_type_G_T2_sum[t][c];
        NpyWriter::WriteFloat64(output_dir + "/bs_per_type_T2.npy", data.data(), N, 40);
        written++;
    }

    // bs_total_B: (N, 3) — total B-field vector at each atom
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& B = conf.AtomAt(i).total_B_field;
            data[i*3+0] = B.x(); data[i*3+1] = B.y(); data[i*3+2] = B.z();
        }
        NpyWriter::WriteFloat64(output_dir + "/bs_total_B.npy", data.data(), N, 3);
        written++;
    }

    // bs_ring_counts: (N, 4) — proximity counts at 3/5/8/12 A shells
    {
        std::vector<double> data(N * 4);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            data[i*4+0] = static_cast<double>(ca.n_rings_within_3A);
            data[i*4+1] = static_cast<double>(ca.n_rings_within_5A);
            data[i*4+2] = static_cast<double>(ca.n_rings_within_8A);
            data[i*4+3] = static_cast<double>(ca.n_rings_within_12A);
        }
        NpyWriter::WriteFloat64(output_dir + "/bs_ring_counts.npy", data.data(), N, 4);
        written++;
    }

    return written;
}

}  // namespace nmr
```

## EXAMPLE: McConnellResult.cpp (already instrumented — study this pattern)

```cpp
#include "McConnellResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <algorithm>

namespace nmr {


std::vector<std::type_index> McConnellResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// The full McConnell shielding tensor from one bond at one atom.
//
// M_ab = 9 cos_theta d_hat_a b_hat_b
//      - 3 b_hat_a b_hat_b
//      - (3 d_hat_a d_hat_b - delta_ab)
//
// Returns M_ab / r^3 (Angstrom^-3).
//
// Also computes:
//   K_ab = (3 d_hat_a d_hat_b - delta_ab) / r^3  (symmetric traceless)
//   f    = (3 cos^2 theta - 1) / r^3              (McConnell scalar)
//
// Three terms in M:
//   Term 1: 9 cos_theta d_hat ⊗ b_hat     — asymmetric, gives T1
//   Term 2: -3 b_hat ⊗ b_hat              — symmetric, gives T0
//   Term 3: -(3 d_hat ⊗ d_hat - I)        — symmetric traceless, gives T2
// ============================================================================

struct BondKernelResult {
    Mat3 M_over_r3 = Mat3::Zero();   // full McConnell tensor (asymmetric)
    Mat3 K = Mat3::Zero();           // symmetric traceless dipolar kernel
    double f = 0.0;                  // McConnell scalar
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();   // unit vector from bond midpoint to atom
};


static BondKernelResult ComputeBondKernel(
        const Vec3& atom_pos,
        const Vec3& bond_midpoint,
        const Vec3& bond_direction) {

    BondKernelResult result;

    Vec3 d = atom_pos - bond_midpoint;
    double r = d.norm();

    if (r < MIN_DISTANCE) return result;

    result.distance = r;

    double r3 = r * r * r;
    Vec3 d_hat = d / r;
    result.direction = d_hat;

    double cos_theta = d_hat.dot(bond_direction);

    // McConnell scalar: (3 cos^2 theta - 1) / r^3
    result.f = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    // Symmetric traceless dipolar kernel K_ab
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = (3.0 * d_hat(a) * d_hat(b)
                              - (a == b ? 1.0 : 0.0)) / r3;

    // Full McConnell tensor M_ab / r^3
    //   = [9 cos_theta d_hat_a b_hat_b - 3 b_hat_a b_hat_b
    //      - (3 d_hat_a d_hat_b - delta_ab)] / r^3
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.M_over_r3(a, b) =
                (9.0 * cos_theta * d_hat(a) * bond_direction(b)
                 - 3.0 * bond_direction(a) * bond_direction(b)
                 - (3.0 * d_hat(a) * d_hat(b) - (a == b ? 1.0 : 0.0)))
                / r3;
        }
    }

    return result;
}


// ============================================================================
// McConnellResult::Compute
// ============================================================================

std::unique_ptr<McConnellResult> McConnellResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("McConnellResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " bonds=" + std::to_string(conf.ProteinRef().BondCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_bonds = protein.BondCount();

    auto result_ptr = std::make_unique<McConnellResult>();
    result_ptr->conf_ = &conf;

    // Filter set: SelfSourceFilter (atom is bond endpoint) +
    // DipolarNearFieldFilter (source extent = bond length).
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcMcConnell, "McConnellResult::Compute",
        "filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);

    int total_pairs = 0;
    int filtered_out = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        // Find nearby bonds via spatial index
        auto nearby_bonds = spatial.BondsWithinRadius(atom_pos, MCCONNELL_CUTOFF_A);

        // Per-category accumulators
        double co_sum = 0.0, cn_sum = 0.0, sidechain_sum = 0.0, aromatic_sum = 0.0;
        Mat3 M_backbone_total = Mat3::Zero();
        Mat3 M_sidechain_total = Mat3::Zero();
        Mat3 M_aromatic_total = Mat3::Zero();
        Mat3 M_total = Mat3::Zero();

        // Nearest CO and CN tracking
        double best_co_dist = NO_DATA_SENTINEL;
        double best_cn_dist = NO_DATA_SENTINEL;
        double best_co_f = 0.0;
        Vec3 best_co_midpoint = Vec3::Zero();
        Vec3 best_co_direction = Vec3::Zero();
        BondKernelResult best_co_kernel;
        BondKernelResult best_cn_kernel;

        for (size_t bi : nearby_bonds) {
            const Bond& bond = protein.BondAt(bi);
            Vec3 midpoint = conf.bond_midpoints[bi];
            Vec3 direction = conf.bond_directions[bi];

            BondKernelResult kernel = ComputeBondKernel(atom_pos, midpoint, direction);

            // Build evaluation context and apply filter set
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = conf.bond_lengths[bi];
            ctx.atom_index = ai;
            ctx.source_atom_a = bond.atom_index_a;
            ctx.source_atom_b = bond.atom_index_b;
            if (!filters.AcceptAll(ctx)) {
                // ---- GeometryChoice: filter exclusion ----
                choices.Record(CalculatorId::McConnell, bi, "filter exclusion",
                    [&](GeometryChoice& gc) {
                        AddBond(gc, &bond, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                filtered_out++; continue;
            }

            // Store in BondNeighbourhood
            BondNeighbourhood bn;
            bn.bond_index = bi;
            bn.bond_category = bond.category;
            bn.distance_to_midpoint = kernel.distance;
            bn.direction_to_midpoint = kernel.direction;
            bn.dipolar_tensor = kernel.K;
            bn.dipolar_spherical = SphericalTensor::Decompose(kernel.K);
            bn.mcconnell_scalar = kernel.f;
            ca.bond_neighbours.push_back(bn);

            // Accumulate full McConnell tensor by category
            M_total += kernel.M_over_r3;

            switch (bond.category) {
                case BondCategory::PeptideCO:
                    co_sum += kernel.f;
                    M_backbone_total += kernel.M_over_r3;
                    if (kernel.distance < best_co_dist) {
                        best_co_dist = kernel.distance;
                        best_co_f = kernel.f;
                        best_co_midpoint = midpoint;
                        best_co_direction = kernel.direction;
                        best_co_kernel = kernel;
                    }
                    break;

                case BondCategory::PeptideCN:
                    cn_sum += kernel.f;
                    M_backbone_total += kernel.M_over_r3;
                    if (kernel.distance < best_cn_dist) {
                        best_cn_dist = kernel.distance;
                        best_cn_kernel = kernel;
                    }
                    break;

                case BondCategory::BackboneOther:
                    M_backbone_total += kernel.M_over_r3;
                    break;

                case BondCategory::SidechainCO:
                    sidechain_sum += kernel.f;
                    M_sidechain_total += kernel.M_over_r3;
                    break;

                case BondCategory::Aromatic:
                    aromatic_sum += kernel.f;
                    M_aromatic_total += kernel.M_over_r3;
                    break;

                case BondCategory::SidechainOther:
                    M_sidechain_total += kernel.M_over_r3;
                    break;

                default:
                    break;
            }

            total_pairs++;
        }

        // Store per-atom totals
        ca.mcconnell_co_sum = co_sum;
        ca.mcconnell_cn_sum = cn_sum;
        ca.mcconnell_sidechain_sum = sidechain_sum;
        ca.mcconnell_aromatic_sum = aromatic_sum;

        // Nearest CO
        ca.mcconnell_co_nearest = best_co_f;
        ca.nearest_CO_midpoint = best_co_midpoint;
        ca.nearest_CO_dist = best_co_dist;
        ca.nearest_CN_dist = best_cn_dist;

        if (best_co_dist < NO_DATA_SENTINEL) {
            double dir_norm = best_co_direction.norm();
            ca.dir_nearest_CO = (dir_norm > NEAR_ZERO_NORM)
                ? Vec3(best_co_direction / dir_norm) : Vec3::Zero();
            ca.T2_CO_nearest = SphericalTensor::Decompose(best_co_kernel.K);
        }
        if (best_cn_dist < NO_DATA_SENTINEL) {
            ca.T2_CN_nearest = SphericalTensor::Decompose(best_cn_kernel.K);
        }

        // Category T2 totals (from symmetric dipolar kernel sums)
        // Apply traceless projection to fix floating-point drift
        auto project_traceless = [](Mat3& m) {
            double trace = m.trace();
            m -= (trace / 3.0) * Mat3::Identity();
        };

        // Extract symmetric part for T2 features
        Mat3 K_backbone = 0.5 * (M_backbone_total + M_backbone_total.transpose());
        K_backbone -= (K_backbone.trace() / 3.0) * Mat3::Identity();
        ca.T2_backbone_total = SphericalTensor::Decompose(K_backbone);

        Mat3 K_sidechain = 0.5 * (M_sidechain_total + M_sidechain_total.transpose());
        K_sidechain -= (K_sidechain.trace() / 3.0) * Mat3::Identity();
        ca.T2_sidechain_total = SphericalTensor::Decompose(K_sidechain);

        Mat3 K_aromatic = 0.5 * (M_aromatic_total + M_aromatic_total.transpose());
        K_aromatic -= (K_aromatic.trace() / 3.0) * Mat3::Identity();
        ca.T2_aromatic_total = SphericalTensor::Decompose(K_aromatic);

        // Full McConnell shielding contribution (from the full M tensor sum)
        ca.mc_shielding_contribution = SphericalTensor::Decompose(M_total);

        // ---- GeometryChoice: bond anisotropy ----
        choices.Record(CalculatorId::McConnell, ai, "bond anisotropy",
            [&ca, ai, best_co_dist, best_cn_dist](GeometryChoice& gc) {
                AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Included);
                AddNumber(gc, "nearest_CO_dist", best_co_dist, "A");
                AddNumber(gc, "nearest_CN_dist", best_cn_dist, "A");
            });
    }

    OperationLog::Info(LogCalcMcConnell, "McConnellResult::Compute",
        "atom_bond_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " bonds=" + std::to_string(n_bonds));

    return result_ptr;
}


// ============================================================================
// Query methods
// ============================================================================

double McConnellResult::CategorySum(size_t atom_index, BondCategory cat) const {
    const auto& ca = conf_->AtomAt(atom_index);
    switch (cat) {
        case BondCategory::PeptideCO: return ca.mcconnell_co_sum;
        case BondCategory::PeptideCN: return ca.mcconnell_cn_sum;
        case BondCategory::SidechainCO: return ca.mcconnell_sidechain_sum;
        case BondCategory::Aromatic: return ca.mcconnell_aromatic_sum;
        default: return 0.0;
    }
}

double McConnellResult::NearestCOContribution(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mcconnell_co_nearest;
}


// ============================================================================
// SampleShieldingAt: evaluate McConnell kernel at arbitrary 3D point.
// Sums over all bonds within MCCONNELL_CUTOFF_A. No atom-specific filters.
// ============================================================================

SphericalTensor McConnellResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 M_total = Mat3::Zero();

    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        auto kernel = ComputeBondKernel(
            point, conf_->bond_midpoints[bi], conf_->bond_directions[bi]);
        if (kernel.distance < MIN_DISTANCE) continue;
        if (kernel.distance > MCCONNELL_CUTOFF_A) continue;

        // DipolarNearFieldFilter: skip if inside the bond
        double bond_len = conf_->bond_lengths[bi];
        if (kernel.distance < 0.5 * bond_len) continue;

        M_total += kernel.M_over_r3;
    }

    return SphericalTensor::Decompose(M_total);
}


// ============================================================================
// WriteFeatures: mc_shielding (9), mc_category_T2 (5 categories × 5 T2),
// mc_scalars (CO/CN/sidechain/aromatic sums, nearest distances).
// ============================================================================

static void PackST_MC(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int McConnellResult::WriteFeatures(const ProteinConformation& conf,
                                    const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> cat_T2(N * 25);
    std::vector<double> scalars(N * 6);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_MC(ca.mc_shielding_contribution, &shielding[i*9]);

        const SphericalTensor* cats[5] = {
            &ca.T2_backbone_total, &ca.T2_sidechain_total,
            &ca.T2_aromatic_total, &ca.T2_CO_nearest, &ca.T2_CN_nearest
        };
        for (int c = 0; c < 5; ++c)
            for (int m = 0; m < 5; ++m)
                cat_T2[i*25 + c*5 + m] = cats[c]->T2[m];

        scalars[i*6+0] = ca.mcconnell_co_sum;
        scalars[i*6+1] = ca.mcconnell_cn_sum;
        scalars[i*6+2] = ca.mcconnell_sidechain_sum;
        scalars[i*6+3] = ca.mcconnell_aromatic_sum;
        scalars[i*6+4] = ca.nearest_CO_dist;
        scalars[i*6+5] = ca.nearest_CN_dist;
    }

    NpyWriter::WriteFloat64(output_dir + "/mc_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mc_category_T2.npy", cat_T2.data(), N, 25);
    NpyWriter::WriteFloat64(output_dir + "/mc_scalars.npy", scalars.data(), N, 6);
    return 3;
}

}  // namespace nmr
```

## EXAMPLE: CoulombResult.cpp (already instrumented — study this pattern)

```cpp
#include "CoulombResult.h"
#include "Protein.h"
#include "ChargeAssignmentResult.h"
#include "SpatialIndexResult.h"
#include "ApbsFieldResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <vector>

namespace nmr {


std::vector<std::type_index> CoulombResult::Dependencies() const {
    return {
        std::type_index(typeid(ChargeAssignmentResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}


// ============================================================================
// CoulombResult::Compute
//
// E_a(i) = ke * sum_{j!=i} q_j * (r_i - r_j)_a / |r_i - r_j|^3
//
// V_ab(i) = ke * sum_{j!=i} q_j * [3 (r_i-r_j)_a (r_i-r_j)_b / |r_i-r_j|^5
//                                   - delta_ab / |r_i-r_j|^3]
//
// ke = 14.3996 V*A  (Coulomb's constant in {e, A, eV} units)
//
// Decomposed by source atom classification:
//   backbone:  N, CA, C, O, H, HA, CB (from residue backbone cache)
//   aromatic:  atoms that are members of any Ring
//   sidechain: everything else
// ============================================================================

std::unique_ptr<CoulombResult> CoulombResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("CoulombResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t n_atoms = conf.AtomCount();

    auto result_ptr = std::make_unique<CoulombResult>();
    result_ptr->conf_ = &conf;

    // ------------------------------------------------------------------
    // Pre-build atom classification vectors from topology (no EnrichmentResult
    // dependency — uses Residue backbone cache and Ring atom membership).
    // ------------------------------------------------------------------

    std::vector<bool> is_backbone(n_atoms, false);
    std::vector<bool> is_aromatic_atom(n_atoms, false);

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        auto mark_bb = [&](size_t idx) {
            if (idx != Residue::NONE && idx < n_atoms) is_backbone[idx] = true;
        };
        mark_bb(res.N);
        mark_bb(res.CA);
        mark_bb(res.C);
        mark_bb(res.O);
        mark_bb(res.H);
        mark_bb(res.HA);
        mark_bb(res.CB);
    }

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            if (ai < n_atoms) is_aromatic_atom[ai] = true;
        }
    }

    // ------------------------------------------------------------------
    // For E_bond_proj: find each atom's primary bond direction.
    // H atoms: direction from parent heavy atom to H.
    // Heavy atoms: direction of first bond (arbitrary but consistent).
    // ------------------------------------------------------------------

    std::vector<Vec3> primary_bond_dir(n_atoms, Vec3::Zero());
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        if (atom.element == Element::H && atom.parent_atom_index != SIZE_MAX) {
            // H atom: bond direction from parent to H
            Vec3 d = conf.PositionAt(ai) - conf.PositionAt(atom.parent_atom_index);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        } else if (!atom.bond_indices.empty()) {
            // Heavy atom: first bond direction
            const Bond& b = protein.BondAt(atom.bond_indices[0]);
            size_t other = (b.atom_index_a == ai) ? b.atom_index_b : b.atom_index_a;
            Vec3 d = conf.PositionAt(other) - conf.PositionAt(ai);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        }
    }

    // ------------------------------------------------------------------
    // Main N^2 Coulomb sum
    // ------------------------------------------------------------------

    // Filter set: SelfSourceFilter (field undefined at source itself).
    // Coulomb is a point-source sum — no DipolarNearFieldFilter needed.
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());

    GeometryChoiceBuilder choices(conf);

    int aromatic_source_count = 0;
    for (size_t j = 0; j < n_atoms; ++j)
        if (is_aromatic_atom[j]) aromatic_source_count++;

    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos_i = conf.PositionAt(i);

        Vec3 E_total = Vec3::Zero();
        Vec3 E_backbone = Vec3::Zero();
        Vec3 E_sidechain = Vec3::Zero();
        Vec3 E_aromatic = Vec3::Zero();

        Mat3 EFG_total = Mat3::Zero();
        Mat3 EFG_backbone = Mat3::Zero();
        Mat3 EFG_sidechain = Mat3::Zero();
        Mat3 EFG_aromatic = Mat3::Zero();

        int n_sidechain_aromatic_sources = 0;

        for (size_t j = 0; j < n_atoms; ++j) {
            // Self-exclusion via filter framework (not inline check)
            KernelEvaluationContext ctx;
            ctx.atom_index = i;
            ctx.source_atom_a = j;
            ctx.distance = (pos_i - conf.PositionAt(j)).norm();
            if (!filters.AcceptAll(ctx)) continue;

            double q_j = conf.AtomAt(j).partial_charge;
            if (std::abs(q_j) < 1e-15) continue;

            Vec3 r = pos_i - conf.PositionAt(j);
            double r_mag = r.norm();

            double r3 = r_mag * r_mag * r_mag;
            double r5 = r3 * r_mag * r_mag;

            // E_a = q_j * r_a / r^3
            Vec3 E_j = q_j * r / r3;

            // V_ab = q_j * (3 r_a r_b / r^5 - delta_ab / r^3)
            Mat3 V_j = q_j * (3.0 * r * r.transpose() / r5
                              - Mat3::Identity() / r3);

            E_total += E_j;
            EFG_total += V_j;

            // Classify source atom
            if (is_aromatic_atom[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
                // Count aromatic sidechain sources near this atom
                // (atoms on aromatic residues that are sidechain)
                if (!is_backbone[j]) n_sidechain_aromatic_sources++;
            } else if (is_backbone[j]) {
                E_backbone += E_j;
                EFG_backbone += V_j;
            } else {
                E_sidechain += E_j;
                EFG_sidechain += V_j;
            }
        }

        // Apply Coulomb constant: convert from e/A^2 to V/A
        E_total     *= COULOMB_KE;
        E_backbone  *= COULOMB_KE;
        E_sidechain *= COULOMB_KE;
        E_aromatic  *= COULOMB_KE;
        EFG_total     *= COULOMB_KE;
        EFG_backbone  *= COULOMB_KE;
        EFG_sidechain *= COULOMB_KE;
        EFG_aromatic  *= COULOMB_KE;

        // Traceless projection on all EFG matrices.
        // Each individual term is traceless (Gauss's law), but floating-point
        // accumulation breaks this. Project to enforce the physics.
        auto project_traceless = [](Mat3& m) {
            m -= (m.trace() / 3.0) * Mat3::Identity();
        };
        project_traceless(EFG_total);
        project_traceless(EFG_backbone);
        project_traceless(EFG_sidechain);
        project_traceless(EFG_aromatic);

        // Sanitise NaN/Inf
        auto sanitise_vec = [](Vec3& v) {
            for (int d = 0; d < 3; ++d)
                if (std::isnan(v(d)) || std::isinf(v(d))) { v = Vec3::Zero(); return; }
        };
        auto sanitise_mat = [](Mat3& m) {
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    if (std::isnan(m(a,b)) || std::isinf(m(a,b))) m(a,b) = 0.0;
        };
        sanitise_vec(E_total);
        sanitise_vec(E_backbone);
        sanitise_vec(E_sidechain);
        sanitise_vec(E_aromatic);
        sanitise_mat(EFG_total);
        sanitise_mat(EFG_backbone);
        sanitise_mat(EFG_sidechain);
        sanitise_mat(EFG_aromatic);

        // Clamp extreme E-field magnitudes
        double E_mag = E_total.norm();
        if (E_mag > APBS_SANITY_LIMIT) {
            double scale = APBS_SANITY_LIMIT / E_mag;

            // ---- GeometryChoice: E-field clamp ----
            choices.Record(CalculatorId::Coulomb, i, "E-field clamp",
                [&conf, i, E_mag, scale](GeometryChoice& gc) {
                    AddAtom(gc, &conf.AtomAt(i), i, EntityRole::Target, EntityOutcome::Triggered);
                    AddNumber(gc, "actual_E_magnitude", E_mag, "V/A");
                    AddNumber(gc, "scale_factor", scale, "");
                });

            E_total     *= scale;
            E_backbone  *= scale;
            E_sidechain *= scale;
            E_aromatic  *= scale;
        }

        // ------------------------------------------------------------------
        // Store on ConformationAtom
        // ------------------------------------------------------------------
        auto& ca = conf.MutableAtomAt(i);

        ca.coulomb_E_total     = E_total;
        ca.coulomb_E_backbone  = E_backbone;
        ca.coulomb_E_sidechain = E_sidechain;
        ca.coulomb_E_aromatic  = E_aromatic;

        ca.coulomb_EFG_total   = EFG_total;
        ca.coulomb_EFG_total_spherical = SphericalTensor::Decompose(EFG_total);

        ca.coulomb_EFG_backbone = EFG_backbone;
        ca.coulomb_EFG_backbone_spherical = SphericalTensor::Decompose(EFG_backbone);

        ca.coulomb_EFG_aromatic = EFG_aromatic;
        ca.coulomb_EFG_aromatic_spherical = SphericalTensor::Decompose(EFG_aromatic);

        // Derived scalars
        ca.coulomb_E_magnitude = E_total.norm();

        // E projected along primary bond direction (for Buckingham E_z)
        ca.coulomb_E_bond_proj = E_total.dot(primary_bond_dir[i]);

        // Backbone projection: component of E_backbone along E_total direction.
        // Positive = backbone field aligned with total; negative = opposed.
        // Bounded by |E_backbone|. Stable near cancellation (unlike |bb|/|total|).
        if (ca.coulomb_E_magnitude > NEAR_ZERO_NORM) {
            Vec3 E_hat = E_total / ca.coulomb_E_magnitude;
            ca.coulomb_E_backbone_frac = E_backbone.dot(E_hat);
        } else {
            ca.coulomb_E_backbone_frac = 0.0;
        }

        // Aromatic E-field derived scalars
        ca.aromatic_E_magnitude = E_aromatic.norm();
        ca.aromatic_E_bond_proj = E_aromatic.dot(primary_bond_dir[i]);
        ca.aromatic_n_sidechain_atoms = n_sidechain_aromatic_sources;

        // Solvent contribution: APBS (solvated) minus vacuum Coulomb.
        // Only meaningful if ApbsFieldResult is present.
        // Both are in V/A, so this is a proper subtraction.
        if (conf.HasResult<ApbsFieldResult>()) {
            ca.coulomb_E_solvent = ca.apbs_efield - E_total;
            ca.coulomb_EFG_solvent = ca.apbs_efg - EFG_total;
        }

        // Shielding contribution: SphericalTensor of the total EFG.
        //
        // WARNING: This is ONLY the T2 geometric kernel. T0 is zero here
        // (EFG is traceless by Gauss's law). The actual T0 Coulomb shielding
        // comes from Buckingham's A*E_z + B*E_z^2, which requires element-
        // dependent parameters and a bond direction — it is NOT a pure
        // geometric kernel like McConnell's M/r^3.
        //
        // Unlike McConnell (where χ·K tensor contraction produces an
        // asymmetric tensor with non-zero T0+T1+T2 from geometry alone),
        // Coulomb has two SEPARATE geometric kernels:
        //   E_a  (rank-1) → T0 shielding via Buckingham A,B parameters
        //   V_ab (rank-2, symmetric, traceless) → T2 shielding via γ
        // There is no single "full tensor" that unifies them.
        //
        // This field stores Decompose(EFG) only. When Buckingham parameters
        // are available (from ParameterCorrectionResult), the full shielding
        // contribution including T0 from E should be computed and stored.
        ca.coulomb_shielding_contribution = SphericalTensor::Decompose(EFG_total);
    }

    OperationLog::Info(LogCalcOther, "CoulombResult::Compute",
        "atoms=" + std::to_string(n_atoms) +
        " aromatic_sources=" + std::to_string(aromatic_source_count) +
        " rejected={" + filters.ReportRejections() + "}");

    return result_ptr;
}


// ============================================================================
// Query methods
// ============================================================================

Vec3 CoulombResult::EFieldAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).coulomb_E_total;
}

Mat3 CoulombResult::EFGAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).coulomb_EFG_total;
}

SphericalTensor CoulombResult::EFGSphericalAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).coulomb_EFG_total_spherical;
}


Vec3 CoulombResult::SampleEFieldAt(Vec3 point) const {
    if (!conf_) return Vec3::Zero();

    const Protein& protein = conf_->ProteinRef();
    Vec3 E = Vec3::Zero();

    for (size_t j = 0; j < conf_->AtomCount(); ++j) {
        double q = conf_->AtomAt(j).partial_charge;
        if (std::abs(q) < 1e-15) continue;

        Vec3 d = point - conf_->PositionAt(j);
        double r = d.norm();
        if (r < MIN_DISTANCE) continue;

        double r3 = r * r * r;
        E += q * d / r3;
    }

    return E * COULOMB_KE;  // V/A
}


// ============================================================================
// WriteFeatures: coulomb_shielding (9), E-field (3), EFG decompositions,
// scalar features (magnitude, bond projection, backbone fraction).
// ============================================================================

static void PackST_C(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int CoulombResult::WriteFeatures(const ProteinConformation& conf,
                                  const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> efield(N * 3);
    std::vector<double> efg_bb(N * 9);
    std::vector<double> efg_aro(N * 9);
    std::vector<double> scalars(N * 4);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_C(ca.coulomb_shielding_contribution, &shielding[i*9]);

        efield[i*3+0] = ca.coulomb_E_total.x();
        efield[i*3+1] = ca.coulomb_E_total.y();
        efield[i*3+2] = ca.coulomb_E_total.z();

        PackST_C(ca.coulomb_EFG_backbone_spherical, &efg_bb[i*9]);
        PackST_C(ca.coulomb_EFG_aromatic_spherical, &efg_aro[i*9]);

        scalars[i*4+0] = ca.coulomb_E_magnitude;
        scalars[i*4+1] = ca.coulomb_E_bond_proj;
        scalars[i*4+2] = ca.coulomb_E_backbone_frac;
        scalars[i*4+3] = ca.aromatic_E_magnitude;
    }

    NpyWriter::WriteFloat64(output_dir + "/coulomb_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_E.npy", efield.data(), N, 3);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_efg_backbone.npy", efg_bb.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_efg_aromatic.npy", efg_aro.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_scalars.npy", scalars.data(), N, 4);
    return 5;
}

}  // namespace nmr
```

## EXAMPLE: HaighMallionResult.cpp (already instrumented — study this pattern)

```cpp
#include "HaighMallionResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <array>
#include <set>

namespace nmr {


std::vector<std::type_index> HaighMallionResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// 7-point Gaussian quadrature on a triangle (Stroud T2:5-1 / Dunavant degree-5).
//
// Barycentric coordinates (lambda_0, lambda_1, lambda_2) and weights.
// Three orbits: centroid (1 point), near-vertex (3 points), near-edge (3 points).
// Weights sum to 1.0 (normalised to reference triangle area = 1/2).
// ============================================================================

struct TriQuadPoint {
    double lambda[3];
    double weight;
};

static const std::array<TriQuadPoint, 7>& Gauss7() {
    static const double sqrt15 = std::sqrt(15.0);
    static const double a1 = (6.0 - sqrt15) / 21.0;          // ~0.1013
    static const double a2 = (6.0 + sqrt15) / 21.0;          // ~0.4701
    static const double w0 = 9.0 / 40.0;                      //  0.225
    static const double w1 = (155.0 - sqrt15) / 1200.0;       // ~0.1259
    static const double w2 = (155.0 + sqrt15) / 1200.0;       // ~0.1324

    static const std::array<TriQuadPoint, 7> pts = {{
        // Centroid
        {{ 1.0/3.0, 1.0/3.0, 1.0/3.0 }, w0},
        // Orbit 1 — near vertices
        {{ a1, a1, 1.0 - 2.0*a1 }, w1},
        {{ a1, 1.0 - 2.0*a1, a1 }, w1},
        {{ 1.0 - 2.0*a1, a1, a1 }, w1},
        // Orbit 2 — near edge midpoints
        {{ a2, a2, 1.0 - 2.0*a2 }, w2},
        {{ a2, 1.0 - 2.0*a2, a2 }, w2},
        {{ 1.0 - 2.0*a2, a2, a2 }, w2},
    }};
    return pts;
}


// ============================================================================
// Accumulate the dipolar kernel integral over one triangle.
//
// H_ab += integral_triangle [ 3 rho_a rho_b / rho^5 - delta_ab / rho^3 ] dS
//
// where rho = r - r_s (field point minus surface point).
//
// Uses 7-point Gaussian quadrature in barycentric coordinates.
// Triangle area computed from cross product of two edges.
// ============================================================================

static void AccumulateTensor(
        const Vec3& v0, const Vec3& v1, const Vec3& v2,
        const Vec3& r,
        const std::array<TriQuadPoint, 7>& qpts,
        Mat3& H) {

    double triArea = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    if (triArea < 1e-20) return;

    for (const auto& qp : qpts) {
        // Surface point in barycentric coordinates
        Vec3 rS = qp.lambda[0] * v0 + qp.lambda[1] * v1 + qp.lambda[2] * v2;
        Vec3 rho = r - rS;
        double rhoMag = rho.norm();
        if (rhoMag < MIN_DISTANCE) continue;

        double rho3 = rhoMag * rhoMag * rhoMag;
        double rho5 = rho3 * rhoMag * rhoMag;

        // K_ab = 3 rho_a rho_b / rho^5 - delta_ab / rho^3
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                H(a, b) += qp.weight * triArea *
                           (3.0 * rho(a) * rho(b) / rho5
                            - (a == b ? 1.0 : 0.0) / rho3);
    }
}


// ============================================================================
// Adaptive subdivision: when the field point is close to a triangle,
// subdivide into 4 sub-triangles at edge midpoints for better accuracy.
//
// Level 0 -> 1: if any vertex within 2.0 A of field point
// Level 1 -> 2: if any vertex within 1.0 A of field point
// Max depth: 2 (7 -> 28 -> 112 quadrature points per fan triangle)
// ============================================================================

constexpr double HM_SUBDIV_THRESHOLD_L1 = 2.0;  // Angstroms
constexpr double HM_SUBDIV_THRESHOLD_L2 = 1.0;  // Angstroms

static bool NeedsSubdivision(
        const Vec3& v0, const Vec3& v1, const Vec3& v2,
        const Vec3& r, double threshold) {
    return (r - v0).norm() < threshold
        || (r - v1).norm() < threshold
        || (r - v2).norm() < threshold;
}

static void AccumulateAdaptive(
        const Vec3& v0, const Vec3& v1, const Vec3& v2,
        const Vec3& r,
        const std::array<TriQuadPoint, 7>& qpts,
        Mat3& H, int level) {

    bool subdivide = false;
    if (level == 0)
        subdivide = NeedsSubdivision(v0, v1, v2, r, HM_SUBDIV_THRESHOLD_L1);
    else if (level == 1)
        subdivide = NeedsSubdivision(v0, v1, v2, r, HM_SUBDIV_THRESHOLD_L2);

    if (subdivide && level < 2) {
        Vec3 m01 = 0.5 * (v0 + v1);
        Vec3 m12 = 0.5 * (v1 + v2);
        Vec3 m02 = 0.5 * (v0 + v2);
        AccumulateAdaptive(v0,  m01, m02, r, qpts, H, level + 1);
        AccumulateAdaptive(m01, v1,  m12, r, qpts, H, level + 1);
        AccumulateAdaptive(m02, m12, v2,  r, qpts, H, level + 1);
        AccumulateAdaptive(m01, m12, m02, r, qpts, H, level + 1);
    } else {
        AccumulateTensor(v0, v1, v2, r, qpts, H);
    }
}


// ============================================================================
// Compute the HM surface integral for one ring at one atom position.
//
// Fan triangulation: n triangles from ring centroid to consecutive vertex pairs.
// Returns H_ab (symmetric, traceless, units Angstrom^-1).
// ============================================================================

static Mat3 SurfaceIntegral(
        const Vec3& point,
        const RingGeometry& geom) {

    const auto& verts = geom.vertices;
    int nv = static_cast<int>(verts.size());
    if (nv < 3) return Mat3::Zero();

    const auto& qpts = Gauss7();
    Mat3 H = Mat3::Zero();

    for (int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        AccumulateAdaptive(geom.center, verts[i], verts[j],
                           point, qpts, H, 0);
    }

    return H;
}


// ============================================================================
// HaighMallionResult::Compute
//
// For each atom, find rings within RING_CALC_CUTOFF. For each ring:
//   1. Compute H_ab = surface integral of dipolar kernel (symmetric, traceless)
//   2. Compute V = H . n (effective B-field from magnetised surface)
//   3. Construct G_ab = n_b * V_a (full shielding kernel, rank-1)
//   4. Store both H and G on RingNeighbourhood
//   5. Accumulate per-type T0 and T2 sums on ConformationAtom
// ============================================================================

std::unique_ptr<HaighMallionResult> HaighMallionResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("HaighMallionResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " rings=" + std::to_string(conf.ProteinRef().RingCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_rings = protein.RingCount();

    auto result_ptr = std::make_unique<HaighMallionResult>();
    result_ptr->conf_ = &conf;

    if (n_rings == 0) {
        OperationLog::Info(LogCalcHaighMal, "HaighMallionResult::Compute",
            "no rings — nothing to compute");
        return result_ptr;
    }

    // Filter set: DipolarNearFieldFilter with source_extent = ring diameter,
    // plus RingBondedExclusionFilter for topological exclusion.
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcHaighMal, "HaighMallionResult::Compute",
        "filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);
    std::set<size_t> recorded_rings;

    int total_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, RING_CALC_CUTOFF);

        Mat3 G_total = Mat3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            // ---- GeometryChoice: surface integral ----
            if (recorded_rings.find(ri) == recorded_rings.end()) {
                recorded_rings.insert(ri);
                choices.Record(CalculatorId::HaighMallion, ri, "surface integral",
                    [&ring](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        // TODO: HM sampler needs AccumulateTensor refactor
                    });
            }

            double distance = (atom_pos - geom.center).norm();

            // Apply filter
            KernelEvaluationContext ctx;
            ctx.distance = distance;
            ctx.source_extent = 2.0 * geom.radius;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) {
                // ---- GeometryChoice: near-field exclusion ----
                choices.Record(CalculatorId::HaighMallion, ri, "near-field exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                continue;
            }

            // ---- GeometryChoice: adaptive refinement ----
            if (distance < HM_SUBDIV_THRESHOLD_L1) {
                choices.Record(CalculatorId::HaighMallion, ri, "adaptive refinement",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Triggered);
                        AddNumber(gc, "distance", distance, "A");
                        AddNumber(gc, "L1_threshold", HM_SUBDIV_THRESHOLD_L1, "A");
                        AddNumber(gc, "L2_threshold", HM_SUBDIV_THRESHOLD_L2, "A");
                    });
            }

            // Step 1: Raw surface integral H_ab (symmetric, traceless, A^-1)
            Mat3 H = SurfaceIntegral(atom_pos, geom);

            // Step 2: Effective B-field V = H . n
            Vec3 V = H * geom.normal;

            // Step 3: Full shielding kernel G_ab = -n_b * V_a (rank-1)
            // Minus sign from sigma_ab = -dB_a^sec / dB_{0,b}.
            // Same convention as BiotSavartResult: sigma = I * G gives
            // correct sign with literature I (negative for diamagnetic).
            Mat3 G;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    G(a, b) = -geom.normal(b) * V(a);

            // Find or create RingNeighbourhood for this ring
            RingNeighbourhood* rn = nullptr;
            for (auto& existing : ca.ring_neighbours) {
                if (existing.ring_index == ri) {
                    rn = &existing;
                    break;
                }
            }
            if (!rn) {
                RingNeighbourhood new_rn;
                new_rn.ring_index = ri;
                new_rn.ring_type = ring.type_index;
                new_rn.distance_to_center = distance;

                Vec3 d = atom_pos - geom.center;
                new_rn.direction_to_center = d.normalized();

                double z = d.dot(geom.normal);
                Vec3 d_plane = d - z * geom.normal;
                double rho = d_plane.norm();
                double theta = std::atan2(d_plane.norm(), std::abs(z));
                new_rn.z = z;
                new_rn.rho = rho;
                new_rn.theta = theta;

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            // Store HM results on RingNeighbourhood
            rn->hm_tensor = H;                                    // raw integral (symmetric, traceless)
            rn->hm_spherical = SphericalTensor::Decompose(H);     // should be pure T2
            rn->hm_B_field = V;                                   // effective B-field

            // Accumulate the full shielding kernel G
            G_total += G;

            // Per-type T0 and T2 sums (from the full shielding kernel G)
            SphericalTensor G_st = SphericalTensor::Decompose(G);
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < 8) {
                ca.per_type_hm_T0_sum[ti] += G_st.T0;
                for (int c = 0; c < 5; ++c)
                    ca.per_type_hm_T2_sum[ti][c] += G_st.T2[c];
            }

            total_pairs++;
        }

        // Store accumulated HM shielding contribution (from full kernel G)
        ca.hm_shielding_contribution = SphericalTensor::Decompose(G_total);
    }

    OperationLog::Info(LogCalcHaighMal, "HaighMallionResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


// ============================================================================
// SampleShieldingAt: evaluate HM kernel at arbitrary 3D point.
// Same physics as Compute(). No atom-specific filters (grid points).
// ============================================================================

SphericalTensor HaighMallionResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 G_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        double distance = (point - geom.center).norm();
        if (distance < MIN_DISTANCE) continue;
        if (distance < geom.radius) continue;
        if (distance > RING_CALC_CUTOFF) continue;

        Mat3 H = SurfaceIntegral(point, geom);
        Vec3 V = H * geom.normal;

        // G_ab = -n_b * V_a
        Mat3 G;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                G(a, b) = -geom.normal(b) * V(a);

        G_total += G;
    }

    return SphericalTensor::Decompose(G_total);
}


// ============================================================================
// WriteFeatures: hm_shielding (9), per-type T0 (8), per-type T2 (40).
// Mirrors BiotSavart layout — same ring-type decomposition, different kernel.
// ============================================================================

static void PackST_HM(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int HaighMallionResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * 8);
    std::vector<double> per_type_T2(N * 40);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_HM(ca.hm_shielding_contribution, &shielding[i*9]);
        for (int t = 0; t < 8; ++t) {
            per_type_T0[i*8 + t] = ca.per_type_hm_T0_sum[t];
            for (int c = 0; c < 5; ++c)
                per_type_T2[i*40 + t*5 + c] = ca.per_type_hm_T2_sum[t][c];
        }
    }
    NpyWriter::WriteFloat64(output_dir + "/hm_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/hm_per_type_T0.npy", per_type_T0.data(), N, 8);
    NpyWriter::WriteFloat64(output_dir + "/hm_per_type_T2.npy", per_type_T2.data(), N, 40);
    return 3;
}

}  // namespace nmr
```

## TARGET: src/HBondResult.cpp (instrument this file)

```cpp
#include "HBondResult.h"
#include "Protein.h"
#include "DsspResult.h"
#include "SpatialIndexResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <set>

namespace nmr {


std::vector<std::type_index> HBondResult::Dependencies() const {
    return {
        std::type_index(typeid(DsspResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}


// ============================================================================
// An identified H-bond from DSSP, resolved to atom positions.
//
// DSSP identifies backbone H-bonds by the Kabsch-Sander energy criterion
// between residue pairs. Each H-bond has:
//   - donor residue (N-H donates)
//   - acceptor residue (C=O accepts)
//
// We resolve this to atoms:
//   - donor_N: the backbone N of the donor residue
//   - acceptor_O: the backbone O of the acceptor residue
//   - h_hat: unit direction from donor N to acceptor O
//   - midpoint: midpoint of N...O (the source point for the dipolar field)
//   - distance: |N...O| distance
//   - is_backbone: true (all DSSP H-bonds are backbone)
// ============================================================================

struct ResolvedHBond {
    size_t donor_N = SIZE_MAX;
    size_t acceptor_O = SIZE_MAX;
    size_t donor_residue = SIZE_MAX;
    size_t acceptor_residue = SIZE_MAX;
    Vec3 midpoint = Vec3::Zero();
    Vec3 h_hat = Vec3::Zero();        // donor N → acceptor O direction
    double distance = 0.0;            // N...O distance
    int sequence_separation = 0;
};


// ============================================================================
// The full H-bond dipolar tensor from one H-bond at one atom.
//
// Same derivation as McConnell with b_hat → h_hat:
//
//   M_ab = 9 cosθ d̂_a h_b  -  3 h_a h_b  -  (3 d̂_a d̂_b - δ_ab)
//
// Returns M_ab / r³ (Angstrom⁻³).
// ============================================================================

struct HBondKernelResult {
    Mat3 M_over_r3 = Mat3::Zero();
    double f = 0.0;
    double distance = 0.0;
};


static HBondKernelResult ComputeHBondKernel(
        const Vec3& atom_pos,
        const Vec3& hbond_midpoint,
        const Vec3& h_hat) {

    HBondKernelResult result;

    Vec3 d = atom_pos - hbond_midpoint;
    double r = d.norm();

    if (r < MIN_DISTANCE) return result;

    result.distance = r;
    double r3 = r * r * r;
    Vec3 d_hat = d / r;
    double cos_theta = d_hat.dot(h_hat);

    result.f = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.M_over_r3(a, b) =
                (9.0 * cos_theta * d_hat(a) * h_hat(b)
                 - 3.0 * h_hat(a) * h_hat(b)
                 - (3.0 * d_hat(a) * d_hat(b) - (a == b ? 1.0 : 0.0)))
                / r3;
        }
    }

    return result;
}


// ============================================================================
// HBondResult::Compute
// ============================================================================

std::unique_ptr<HBondResult> HBondResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("HBondResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& dssp = conf.Result<DsspResult>();
    // SpatialIndexResult is a declared dependency (ensures it's computed
    // before us) but we iterate H-bonds directly, not via spatial search.
    (void)conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_residues = protein.ResidueCount();

    auto result_ptr = std::make_unique<HBondResult>();
    result_ptr->conf_ = &conf;

    // ------------------------------------------------------------------
    // Step 1: Resolve DSSP H-bond partners to atom positions.
    //
    // DSSP provides up to 2 acceptor partners and 2 donor partners per
    // residue. Each partner is a residue index. We resolve to backbone
    // N (donor) and O (acceptor) atoms.
    //
    // Skip H-bonds where:
    //   - partner residue index is SIZE_MAX (no partner)
    //   - backbone N or O atoms are not present
    //   - sequence separation < SEQUENTIAL_EXCLUSION_THRESHOLD (too close)
    //   - N...O distance > HBOND_MAX_DIST
    // ------------------------------------------------------------------

    std::vector<ResolvedHBond> hbonds;

    // Use a set to deduplicate: (donor_N, acceptor_O) pairs
    std::set<std::pair<size_t, size_t>> seen;

    for (size_t ri = 0; ri < n_residues; ++ri) {
        const auto& dr = dssp.AllResidues()[ri];
        const Residue& res = protein.ResidueAt(ri);

        // This residue's N-H donates to acceptor residues
        for (int bi = 0; bi < 2; ++bi) {
            size_t acc_ri = dr.acceptors[bi].residue_index;
            if (acc_ri == SIZE_MAX || acc_ri >= n_residues) continue;

            const Residue& acc_res = protein.ResidueAt(acc_ri);
            if (res.N == Residue::NONE || acc_res.O == Residue::NONE) continue;

            int seq_sep = std::abs(static_cast<int>(ri) - static_cast<int>(acc_ri));
            if (seq_sep < SEQUENTIAL_EXCLUSION_THRESHOLD) continue;

            auto key = std::make_pair(res.N, acc_res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 N_pos = conf.PositionAt(res.N);
            Vec3 O_pos = conf.PositionAt(acc_res.O);
            Vec3 d = O_pos - N_pos;
            double dist = d.norm();

            if (dist < MIN_DISTANCE || dist > HBOND_MAX_DIST) continue;

            ResolvedHBond hb;
            hb.donor_N = res.N;
            hb.acceptor_O = acc_res.O;
            hb.donor_residue = ri;
            hb.acceptor_residue = acc_ri;
            hb.midpoint = 0.5 * (N_pos + O_pos);
            hb.h_hat = d / dist;
            hb.distance = dist;
            hb.sequence_separation = seq_sep;
            hbonds.push_back(hb);
        }

        // This residue's C=O accepts from donor residues
        for (int bi = 0; bi < 2; ++bi) {
            size_t don_ri = dr.donors[bi].residue_index;
            if (don_ri == SIZE_MAX || don_ri >= n_residues) continue;

            const Residue& don_res = protein.ResidueAt(don_ri);
            if (don_res.N == Residue::NONE || res.O == Residue::NONE) continue;

            int seq_sep = std::abs(static_cast<int>(ri) - static_cast<int>(don_ri));
            if (seq_sep < SEQUENTIAL_EXCLUSION_THRESHOLD) continue;

            auto key = std::make_pair(don_res.N, res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 N_pos = conf.PositionAt(don_res.N);
            Vec3 O_pos = conf.PositionAt(res.O);
            Vec3 d = O_pos - N_pos;
            double dist = d.norm();

            if (dist < MIN_DISTANCE || dist > HBOND_MAX_DIST) continue;

            ResolvedHBond hb;
            hb.donor_N = don_res.N;
            hb.acceptor_O = res.O;
            hb.donor_residue = don_ri;
            hb.acceptor_residue = ri;
            hb.midpoint = 0.5 * (N_pos + O_pos);
            hb.h_hat = d / dist;
            hb.distance = dist;
            hb.sequence_separation = seq_sep;
            hbonds.push_back(hb);
        }
    }

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "resolved " + std::to_string(hbonds.size()) +
        " unique backbone H-bonds from DSSP");

    // Store resolved H-bond geometry for SampleAt grid queries
    for (const auto& hb : hbonds) {
        result_ptr->hbond_midpoints_.push_back(hb.midpoint);
        result_ptr->hbond_directions_.push_back(hb.h_hat);
        result_ptr->hbond_distances_.push_back(hb.distance);
    }

    if (hbonds.empty()) {
        return result_ptr;
    }

    // ------------------------------------------------------------------
    // Step 2: Build the filter set for H-bond kernel evaluations.
    //
    // SelfSourceFilter: atom cannot be a field point for an H-bond
    //   where it is the donor N or acceptor O.
    // DipolarNearFieldFilter: point-source model invalid when field
    //   point is inside the N...O source distribution.
    // ------------------------------------------------------------------

    KernelFilterSet filters;
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "filter set: " + filters.Describe());

    // ------------------------------------------------------------------
    // Step 3: For each atom, compute the dipolar tensor from all H-bonds
    // that pass the filter set. 1/r³ decay handles range naturally.
    // ------------------------------------------------------------------

    int total_pairs = 0;
    int filtered_out = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        // Residue index of this atom (for sequence separation)
        size_t ai_res = protein.AtomAt(ai).residue_index;

        Mat3 M_total = Mat3::Zero();
        double nearest_dist = NO_DATA_SENTINEL;
        size_t nearest_hb_idx = SIZE_MAX;
        int count_3_5 = 0;

        for (size_t hi = 0; hi < hbonds.size(); ++hi) {
            const auto& hb = hbonds[hi];

            HBondKernelResult kernel = ComputeHBondKernel(
                atom_pos, hb.midpoint, hb.h_hat);

            if (kernel.distance < MIN_DISTANCE) continue;

            // Build evaluation context from already-computed geometry
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = hb.distance;  // N...O distance
            ctx.atom_index = ai;
            ctx.source_atom_a = hb.donor_N;
            ctx.source_atom_b = hb.acceptor_O;

            // Sequence separation: min distance to either endpoint residue
            int sep_don = std::abs(static_cast<int>(ai_res)
                                 - static_cast<int>(hb.donor_residue));
            int sep_acc = std::abs(static_cast<int>(ai_res)
                                 - static_cast<int>(hb.acceptor_residue));
            ctx.sequence_separation = std::min(sep_don, sep_acc);

            if (!filters.AcceptAll(ctx)) {
                filtered_out++;
                continue;
            }

            // Count H-bonds within 3.5A of this atom
            if (kernel.distance < HBOND_COUNT_RADIUS) count_3_5++;

            // Track nearest (among accepted evaluations only)
            if (kernel.distance < nearest_dist) {
                nearest_dist = kernel.distance;
                nearest_hb_idx = hi;
            }

            // Accumulate tensor from all H-bonds (1/r³ decay handles range)
            M_total += kernel.M_over_r3;
            total_pairs++;
        }

        ca.hbond_count_within_3_5A = count_3_5;

        if (nearest_hb_idx != SIZE_MAX) {
            const auto& nearest_hb = hbonds[nearest_hb_idx];
            ca.hbond_nearest_dist = nearest_dist;
            ca.hbond_nearest_dir = (atom_pos - nearest_hb.midpoint).normalized();
            ca.hbond_is_backbone = true;  // all DSSP H-bonds are backbone

            HBondKernelResult nearest_kernel = ComputeHBondKernel(
                atom_pos, nearest_hb.midpoint, nearest_hb.h_hat);

            ca.hbond_nearest_tensor = nearest_kernel.M_over_r3;
            ca.hbond_nearest_spherical = SphericalTensor::Decompose(
                nearest_kernel.M_over_r3);
            ca.hbond_inv_d3 = 1.0 / (nearest_dist * nearest_dist * nearest_dist);
        }

        // Check if this atom is a donor or acceptor
        for (const auto& hb : hbonds) {
            if (ai == hb.donor_N) ca.hbond_is_donor = true;
            if (ai == hb.acceptor_O) ca.hbond_is_acceptor = true;
        }

        ca.hbond_shielding_contribution = SphericalTensor::Decompose(M_total);
    }

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "atom_hbond_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " hbonds=" + std::to_string(hbonds.size()) +
        " atoms=" + std::to_string(n_atoms));

    return result_ptr;
}


SphericalTensor HBondResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_ || hbond_midpoints_.empty()) return SphericalTensor{};

    Mat3 M_total = Mat3::Zero();

    for (size_t hi = 0; hi < hbond_midpoints_.size(); ++hi) {
        auto kernel = ComputeHBondKernel(
            point, hbond_midpoints_[hi], hbond_directions_[hi]);
        if (kernel.distance < MIN_DISTANCE) continue;

        // DipolarNearFieldFilter: skip if inside the N...O distribution
        if (kernel.distance < 0.5 * hbond_distances_[hi]) continue;

        M_total += kernel.M_over_r3;
    }

    return SphericalTensor::Decompose(M_total);
}


static void PackST_HB(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int HBondResult::WriteFeatures(const ProteinConformation& conf,
                                const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> scalars(N * 3);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_HB(ca.hbond_shielding_contribution, &shielding[i*9]);
        scalars[i*3+0] = ca.hbond_nearest_dist;
        scalars[i*3+1] = ca.hbond_inv_d3;
        scalars[i*3+2] = static_cast<double>(ca.hbond_count_within_3_5A);
    }

    NpyWriter::WriteFloat64(output_dir + "/hbond_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/hbond_scalars.npy", scalars.data(), N, 3);
    return 2;
}

}  // namespace nmr
```

## TARGET: src/PiQuadrupoleResult.cpp (instrument this file)

```cpp
#include "PiQuadrupoleResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>

namespace nmr {


std::vector<std::type_index> PiQuadrupoleResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// EFG geometric kernel from a point axial quadrupole at the ring center.
//
// Stone T-tensor derivation: V_ab = -(Theta/2) * T_abcd n_c n_d
// We compute G_ab = T_abcd n_c n_d (the -Theta/2 prefactor goes into
// the learnable parameter Q_type).
//
//   G_ab = 105 dn^2 d_a d_b / r^9
//        - 30 dn (n_a d_b + n_b d_a) / r^7
//        - 15 d_a d_b / r^7
//        + 6 n_a n_b / r^5
//        + delta_ab (3/r^5 - 15 dn^2/r^7)
//
// Also: scalar = (3 cos^2 theta - 1) / r^4   (Buckingham A-term)
//
// Properties:
//   - Symmetric: yes (every term is symmetric in a,b)
//   - Traceless: yes (Laplace; verified numerically in tests)
//   - Units: G in A^-5, scalar in A^-4
// ============================================================================

struct PiQuadKernelResult {
    Mat3 G = Mat3::Zero();     // EFG geometric kernel (traceless, symmetric)
    double scalar = 0.0;       // (3 cos^2 theta - 1) / r^4
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();
};


static PiQuadKernelResult ComputePiQuadKernel(
        const Vec3& atom_pos,
        const Vec3& ring_center,
        const Vec3& ring_normal) {

    PiQuadKernelResult result;

    Vec3 d = atom_pos - ring_center;
    double r = d.norm();

    if (r < MIN_DISTANCE) return result;

    result.distance = r;
    result.direction = d / r;

    double r2 = r * r;
    double r5 = r2 * r2 * r;
    double r7 = r5 * r2;
    double r9 = r7 * r2;

    double dn = d.dot(ring_normal);       // height above ring plane
    double dn2 = dn * dn;
    double cos_theta = dn / r;

    // Scalar: (3 cos^2 theta - 1) / r^4
    result.scalar = (3.0 * cos_theta * cos_theta - 1.0) / (r2 * r2);

    // EFG tensor: G_ab = T_abcd n_c n_d (Stone Ch. 3)
    //
    //   105 dn^2 d_a d_b / r^9
    //   - 30 dn (n_a d_b + n_b d_a) / r^7
    //   - 15 d_a d_b / r^7
    //   + 6 n_a n_b / r^5
    //   + delta_ab (3/r^5 - 15 dn^2/r^7)
    double diag_term = 3.0 / r5 - 15.0 * dn2 / r7;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.G(a, b) =
                105.0 * dn2 * d(a) * d(b) / r9
                - 30.0 * dn * (ring_normal(a) * d(b) + ring_normal(b) * d(a)) / r7
                - 15.0 * d(a) * d(b) / r7
                + 6.0 * ring_normal(a) * ring_normal(b) / r5
                + (a == b ? diag_term : 0.0);
        }
    }

    return result;
}


// ============================================================================
// PiQuadrupoleResult::Compute
// ============================================================================

std::unique_ptr<PiQuadrupoleResult> PiQuadrupoleResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("PiQuadrupoleResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " rings=" + std::to_string(conf.ProteinRef().RingCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_rings = protein.RingCount();

    auto result_ptr = std::make_unique<PiQuadrupoleResult>();
    result_ptr->conf_ = &conf;

    if (n_rings == 0) {
        OperationLog::Info(LogCalcOther, "PiQuadrupoleResult::Compute",
            "no rings — nothing to compute");
        return result_ptr;
    }

    // Filter set: DipolarNearFieldFilter with source_extent = ring diameter,
    // plus RingBondedExclusionFilter for topological exclusion. The
    // quadrupole approximation is less accurate than the dipolar one
    // at close range (higher multipole → larger convergence radius),
    // making the topology check especially important here.
    KernelFilterSet filters;
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcOther, "PiQuadrupoleResult::Compute",
        "filter set: " + filters.Describe());

    int total_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, RING_CALC_CUTOFF);

        Mat3 G_total = Mat3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            PiQuadKernelResult kernel = ComputePiQuadKernel(
                atom_pos, geom.center, geom.normal);

            if (kernel.distance < MIN_DISTANCE) continue;

            // Apply filter: source extent = ring diameter
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = 2.0 * geom.radius;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) continue;

            // Find or create RingNeighbourhood
            RingNeighbourhood* rn = nullptr;
            for (auto& existing : ca.ring_neighbours) {
                if (existing.ring_index == ri) {
                    rn = &existing;
                    break;
                }
            }
            if (!rn) {
                RingNeighbourhood new_rn;
                new_rn.ring_index = ri;
                new_rn.ring_type = ring.type_index;
                new_rn.distance_to_center = kernel.distance;
                new_rn.direction_to_center = kernel.direction;

                Vec3 d_vec = atom_pos - geom.center;
                double z = d_vec.dot(geom.normal);
                Vec3 d_plane = d_vec - z * geom.normal;
                new_rn.z = z;
                new_rn.rho = d_plane.norm();
                new_rn.theta = std::atan2(d_plane.norm(), std::abs(z));

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            // Store quadrupole kernel on RingNeighbourhood
            rn->quad_tensor = kernel.G;
            rn->quad_spherical = SphericalTensor::Decompose(kernel.G);
            rn->quad_scalar = kernel.scalar;

            // Per-type accumulation
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < 8) {
                ca.per_type_pq_scalar_sum[ti] += kernel.scalar;
                for (int c = 0; c < 5; ++c)
                    ca.per_type_pq_T2_sum[ti][c] += rn->quad_spherical.T2[c];
            }

            // Accumulate EFG tensor
            G_total += kernel.G;
            total_pairs++;
        }

        // Store accumulated shielding contribution (pure T2)
        ca.piquad_shielding_contribution = SphericalTensor::Decompose(G_total);
    }

    OperationLog::Info(LogCalcOther, "PiQuadrupoleResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


SphericalTensor PiQuadrupoleResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 G_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];

        double distance = (point - geom.center).norm();
        if (distance < MIN_DISTANCE) continue;
        if (distance < geom.radius) continue;
        if (distance > RING_CALC_CUTOFF) continue;

        auto kernel = ComputePiQuadKernel(point, geom.center, geom.normal);
        G_total += kernel.G;
    }

    return SphericalTensor::Decompose(G_total);
}


static void PackST_PQ(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int PiQuadrupoleResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * 8);
    std::vector<double> per_type_T2(N * 40);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_PQ(ca.piquad_shielding_contribution, &shielding[i*9]);
        for (int t = 0; t < 8; ++t) {
            per_type_T0[i*8 + t] = ca.per_type_pq_scalar_sum[t];
            for (int c = 0; c < 5; ++c)
                per_type_T2[i*40 + t*5 + c] = ca.per_type_pq_T2_sum[t][c];
        }
    }

    NpyWriter::WriteFloat64(output_dir + "/pq_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/pq_per_type_T0.npy", per_type_T0.data(), N, 8);
    NpyWriter::WriteFloat64(output_dir + "/pq_per_type_T2.npy", per_type_T2.data(), N, 40);
    return 3;
}

}  // namespace nmr
```

## TARGET: src/RingSusceptibilityResult.cpp (instrument this file)

```cpp
#include "RingSusceptibilityResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>

namespace nmr {


std::vector<std::type_index> RingSusceptibilityResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// The full ring susceptibility shielding tensor from one ring at one atom.
//
// Same derivation as McConnell (GEOMETRIC_KERNEL_CATALOGUE.md) with
// b_hat → n_hat (ring normal):
//
//   M_ab = 9 cosθ d̂_a n_b  -  3 n_a n_b  -  (3 d̂_a d̂_b - δ_ab)
//
// Returns M_ab / r³ (Angstrom⁻³).
//
// Also computes:
//   K_ab = (3 d̂_a d̂_b - δ_ab) / r³   (symmetric traceless dipolar kernel)
//   f    = (3 cos²θ - 1) / r³           (ring susceptibility scalar)
//
// Three terms in M:
//   Term 1: 9 cosθ d̂ ⊗ n̂      — asymmetric, gives T1
//   Term 2: -3 n̂ ⊗ n̂           — symmetric, gives T0
//   Term 3: -(3 d̂ ⊗ d̂ - I)    — symmetric traceless, gives T2
// ============================================================================

struct RingChiKernelResult {
    Mat3 M_over_r3 = Mat3::Zero();   // full tensor (asymmetric)
    Mat3 K = Mat3::Zero();           // symmetric traceless dipolar kernel
    double f = 0.0;                  // susceptibility scalar
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();   // unit vector from ring center to atom
};


static RingChiKernelResult ComputeRingChiKernel(
        const Vec3& atom_pos,
        const Vec3& ring_center,
        const Vec3& ring_normal) {

    RingChiKernelResult result;

    Vec3 d = atom_pos - ring_center;
    double r = d.norm();

    if (r < MIN_DISTANCE) return result;

    result.distance = r;

    double r3 = r * r * r;
    Vec3 d_hat = d / r;
    result.direction = d_hat;

    double cos_theta = d_hat.dot(ring_normal);

    // Ring susceptibility scalar: (3 cos²θ - 1) / r³
    result.f = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    // Symmetric traceless dipolar kernel K_ab
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = (3.0 * d_hat(a) * d_hat(b)
                              - (a == b ? 1.0 : 0.0)) / r3;

    // Full tensor M_ab / r³
    //   = [9 cosθ d̂_a n_b - 3 n_a n_b - (3 d̂_a d̂_b - δ_ab)] / r³
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.M_over_r3(a, b) =
                (9.0 * cos_theta * d_hat(a) * ring_normal(b)
                 - 3.0 * ring_normal(a) * ring_normal(b)
                 - (3.0 * d_hat(a) * d_hat(b) - (a == b ? 1.0 : 0.0)))
                / r3;
        }
    }

    return result;
}


// ============================================================================
// RingSusceptibilityResult::Compute
// ============================================================================

std::unique_ptr<RingSusceptibilityResult> RingSusceptibilityResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("RingSusceptibilityResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " rings=" + std::to_string(conf.ProteinRef().RingCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_rings = protein.RingCount();

    auto result_ptr = std::make_unique<RingSusceptibilityResult>();
    result_ptr->conf_ = &conf;

    if (n_rings == 0) {
        OperationLog::Info(LogCalcOther, "RingSusceptibilityResult::Compute",
            "no rings — nothing to compute");
        return result_ptr;
    }

    // Filter set: DipolarNearFieldFilter with source_extent = ring diameter,
    // plus RingBondedExclusionFilter for topological exclusion of ring
    // vertices and their bonded neighbours. The distance filter catches
    // most ring atoms by geometry, but the topology check is unambiguous
    // at the boundary (ring atoms sit at exactly the distance threshold).
    KernelFilterSet filters;
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcOther, "RingSusceptibilityResult::Compute",
        "filter set: " + filters.Describe());

    int total_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        // Find nearby rings via spatial index
        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, RING_CALC_CUTOFF);

        Mat3 M_total = Mat3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            RingChiKernelResult kernel = ComputeRingChiKernel(
                atom_pos, geom.center, geom.normal);

            if (kernel.distance < MIN_DISTANCE) continue;

            // Apply filter: source extent = ring diameter (2 * radius)
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = 2.0 * geom.radius;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) continue;

            // Find or create RingNeighbourhood entry for this ring.
            // The ring_neighbours vector may already have an entry from
            // another calculator (BiotSavart). Find by ring index.
            RingNeighbourhood* rn = nullptr;
            for (auto& existing : ca.ring_neighbours) {
                if (existing.ring_index == ri) {
                    rn = &existing;
                    break;
                }
            }
            if (!rn) {
                RingNeighbourhood new_rn;
                new_rn.ring_index = ri;
                new_rn.ring_type = ring.type_index;
                new_rn.distance_to_center = kernel.distance;
                new_rn.direction_to_center = kernel.direction;

                // Cylindrical coordinates in ring frame
                Vec3 d = atom_pos - geom.center;
                double z = d.dot(geom.normal);
                Vec3 d_plane = d - z * geom.normal;
                double rho = d_plane.norm();
                double theta = std::atan2(d_plane.norm(), std::abs(z));
                new_rn.z = z;
                new_rn.rho = rho;
                new_rn.theta = theta;

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            // Store ring susceptibility kernel on this RingNeighbourhood
            rn->chi_tensor = kernel.M_over_r3;
            rn->chi_spherical = SphericalTensor::Decompose(kernel.M_over_r3);
            rn->chi_scalar = kernel.f;

            // Accumulate full tensor
            M_total += kernel.M_over_r3;
            total_pairs++;
        }

        // Store accumulated shielding contribution
        ca.ringchi_shielding_contribution = SphericalTensor::Decompose(M_total);
    }

    OperationLog::Info(LogCalcOther, "RingSusceptibilityResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


SphericalTensor RingSusceptibilityResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 M_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];

        double distance = (point - geom.center).norm();
        if (distance < MIN_DISTANCE) continue;
        if (distance < geom.radius) continue;
        if (distance > RING_CALC_CUTOFF) continue;

        auto kernel = ComputeRingChiKernel(point, geom.center, geom.normal);
        M_total += kernel.M_over_r3;
    }

    return SphericalTensor::Decompose(M_total);
}


static void PackST_RS(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int RingSusceptibilityResult::WriteFeatures(const ProteinConformation& conf,
                                             const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    std::vector<double> shielding(N * 9);
    for (size_t i = 0; i < N; ++i)
        PackST_RS(conf.AtomAt(i).ringchi_shielding_contribution, &shielding[i*9]);
    NpyWriter::WriteFloat64(output_dir + "/ringchi_shielding.npy", shielding.data(), N, 9);
    return 1;
}

}  // namespace nmr
```

## TARGET: src/DispersionResult.cpp (instrument this file)

```cpp
#include "DispersionResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <set>

namespace nmr {


std::vector<std::type_index> DispersionResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// Smooth switching function for dispersion cutoff (CHARMM functional form).
//
// S(r) = 1                                                for r <= R_switch
// S(r) = (Rc²-r²)²(Rc²+2r²-3Rs²) / (Rc²-Rs²)³         for R_switch < r < R_cut
// S(r) = 0                                                for r >= R_cut
//
// C¹ continuous at both boundaries. S(R_switch) = 1, S(R_cut) = 0,
// S'(R_switch) = 0, S'(R_cut) = 0.
//
// Physics reason for smooth taper: the 1/r^6 interaction does not
// physically stop at any distance. The switching function tapers it
// to zero over a finite range so that features vary smoothly with
// atomic position. Essential when the same kernel is evaluated on
// 100+ MD frames where atom positions fluctuate by ~0.5-1A.
//
// Reference: Brooks et al., J. Comput. Chem. 4, 187 (1983) — the
// CHARMM switching function for non-bonded interactions.
//
// R_switch = 4.3 A: onset of taper. Below this, full strength.
// R_cut = 5.0 A: zero beyond this. At R_cut, 1/r^6 = 6.4e-5 A^-6,
//   which is 0.03% of the value at 2A (the typical nearest non-bonded
//   contact). Truncation error from stopping here is < 0.1%.
// ============================================================================

// Dispersion range limits (Angstroms).
// R_CUT: at 5A, C6/r^6 = C6/15625. For the total sum over ~6 vertices
// within range, the contribution beyond 5A is < 0.1% of the total.
// This is a numerical precision choice: we are truncating a convergent
// sum, not asserting that dispersion stops at 5A.
//
// R_SWITCH: onset of the smooth taper. Set at 4.3A (0.7A before R_CUT)
// to give a gentle taper over the range where 1/r^6 is already small.
// The taper width (0.7A) is comparable to atomic position fluctuations
// in MD (~0.5A RMS), ensuring no atom ever jumps across the entire
// taper in one frame.
constexpr double DISP_VERTEX_R_CUT = 5.0;      // Angstroms
constexpr double DISP_VERTEX_R_SWITCH = 4.3;    // Angstroms

static double DispSwitchingFunction(double r) {
    if (r <= DISP_VERTEX_R_SWITCH) return 1.0;
    if (r >= DISP_VERTEX_R_CUT) return 0.0;

    double rc2 = DISP_VERTEX_R_CUT * DISP_VERTEX_R_CUT;
    double rs2 = DISP_VERTEX_R_SWITCH * DISP_VERTEX_R_SWITCH;
    double r2 = r * r;
    double num = (rc2 - r2) * (rc2 - r2) * (rc2 + 2.0 * r2 - 3.0 * rs2);
    double den = (rc2 - rs2) * (rc2 - rs2) * (rc2 - rs2);
    return num / den;
}


// ============================================================================
// London dispersion kernel from one ring vertex at one atom.
//
// Per vertex, with unit C6 = 1:
//
//   K_ab = S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)    (Angstrom^-6)
//   scalar = S(r) / r^6                                    (Angstrom^-6)
//
// where d = r_atom - r_vertex, r = |d|, S(r) is the switching function.
//
// The tensor is traceless per vertex:
//   Tr(K) = S(r) * (3|d|^2/r^8 - 3/r^6) = S(r) * 0 = 0.
//
// Vertex exclusion: atoms covalently bonded to a ring vertex are
// excluded because the through-space 1/r^6 kernel does not model
// through-bond electronic coupling. This is checked via the protein's
// bond connectivity, not a distance heuristic.
// ============================================================================

struct DispVertexResult {
    Mat3 K = Mat3::Zero();
    double scalar = 0.0;
    bool valid = false;
};


static DispVertexResult ComputeDispVertex(
        const Vec3& atom_pos,
        const Vec3& vertex_pos,
        double r) {

    DispVertexResult result;

    if (r < MIN_DISTANCE || r > DISP_VERTEX_R_CUT) return result;

    double S = DispSwitchingFunction(r);
    if (S < 1e-15) return result;  // below switching threshold

    Vec3 d = atom_pos - vertex_pos;
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r8 = r6 * r2;

    result.scalar = S / r6;

    // K_ab = S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = S * (3.0 * d(a) * d(b) / r8
                                - (a == b ? 1.0 : 0.0) / r6);

    result.valid = true;
    return result;
}


// ============================================================================
// Build the set of atoms bonded to any vertex of a ring.
// Used to exclude through-bond pairs from the through-space 1/r^6 kernel.
// ============================================================================

static std::set<size_t> BondedToVertices(
        const Ring& ring, const Protein& protein) {
    std::set<size_t> bonded;
    for (size_t vi : ring.atom_indices) {
        bonded.insert(vi);  // the vertex itself
        const auto& atom = protein.AtomAt(vi);
        for (size_t bi : atom.bond_indices) {
            const auto& bond = protein.BondAt(bi);
            bonded.insert(bond.atom_index_a);
            bonded.insert(bond.atom_index_b);
        }
    }
    return bonded;
}


// ============================================================================
// DispersionResult::Compute
//
// For each atom, find rings within RING_CALC_CUTOFF (15A coarse spatial
// query). Apply DipolarNearFieldFilter at ring level (same as all other
// ring calculators). For each passing ring, sum the dispersion kernel
// over vertices, excluding atoms bonded to any vertex (through-bond
// exclusion). Smooth switching function tapers the 1/r^6 contribution
// to zero at DISP_VERTEX_R_CUT.
// ============================================================================

std::unique_ptr<DispersionResult> DispersionResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("DispersionResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " rings=" + std::to_string(conf.ProteinRef().RingCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_rings = protein.RingCount();

    auto result_ptr = std::make_unique<DispersionResult>();
    result_ptr->conf_ = &conf;

    if (n_rings == 0) {
        OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
            "no rings — nothing to compute");
        return result_ptr;
    }

    // Ring-level filter: DipolarNearFieldFilter with source_extent = ring
    // diameter. Same physics as all other ring calculators: the per-vertex
    // summation is a discrete approximation that breaks down when the
    // field point is inside the ring.
    KernelFilterSet filters;
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
        "filter set: " + filters.Describe() +
        " | vertex range: [MIN_DISTANCE=" + std::to_string(MIN_DISTANCE) +
        ", R_CUT=" + std::to_string(DISP_VERTEX_R_CUT) +
        "] A, switch onset=" + std::to_string(DISP_VERTEX_R_SWITCH) + " A" +
        " | through-bond vertex exclusion: yes");

    // Pre-build bonded-to-vertex sets for each ring (once, not per atom).
    std::vector<std::set<size_t>> ring_bonded(n_rings);
    for (size_t ri = 0; ri < n_rings; ++ri)
        ring_bonded[ri] = BondedToVertices(protein.RingAt(ri), protein);

    int total_pairs = 0;
    int total_contacts = 0;
    int bonded_exclusions = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, RING_CALC_CUTOFF);

        Mat3 disp_total = Mat3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            if (geom.vertices.empty()) continue;

            double dist_to_center = (atom_pos - geom.center).norm();

            // Ring-level filter
            KernelEvaluationContext ctx;
            ctx.distance = dist_to_center;
            ctx.source_extent = 2.0 * geom.radius;  // ring diameter (A)
            ctx.atom_index = ai;
            if (!filters.AcceptAll(ctx)) continue;

            // Through-bond exclusion: skip this ring entirely if the
            // field atom is bonded to any vertex (it's part of or
            // immediately adjacent to the ring).
            if (ring_bonded[ri].count(ai)) {
                bonded_exclusions++;
                continue;
            }

            // Sum dispersion kernel over ring vertices
            Mat3 K_ring = Mat3::Zero();
            double s_ring = 0.0;
            int contacts = 0;

            for (size_t vi = 0; vi < ring.atom_indices.size(); ++vi) {
                Vec3 vpos = geom.vertices[vi];
                double r = (atom_pos - vpos).norm();

                DispVertexResult vr = ComputeDispVertex(atom_pos, vpos, r);
                if (!vr.valid) continue;

                K_ring += vr.K;
                s_ring += vr.scalar;
                contacts++;
            }

            if (contacts == 0) continue;

            // Find or create RingNeighbourhood
            RingNeighbourhood* rn = nullptr;
            for (auto& existing : ca.ring_neighbours) {
                if (existing.ring_index == ri) {
                    rn = &existing;
                    break;
                }
            }
            if (!rn) {
                RingNeighbourhood new_rn;
                new_rn.ring_index = ri;
                new_rn.ring_type = ring.type_index;
                new_rn.distance_to_center = dist_to_center;
                Vec3 d = atom_pos - geom.center;
                if (d.norm() > NEAR_ZERO_NORM)
                    new_rn.direction_to_center = d.normalized();

                double z = d.dot(geom.normal);
                Vec3 d_plane = d - z * geom.normal;
                new_rn.z = z;
                new_rn.rho = d_plane.norm();
                new_rn.theta = std::atan2(d_plane.norm(), std::abs(z));

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            // Store dispersion results on RingNeighbourhood
            rn->disp_tensor = K_ring;
            rn->disp_spherical = SphericalTensor::Decompose(K_ring);
            rn->disp_scalar = s_ring;
            rn->disp_contacts = contacts;

            // Per-type accumulation
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < 8) {
                ca.per_type_disp_scalar_sum[ti] += s_ring;
                for (int c = 0; c < 5; ++c)
                    ca.per_type_disp_T2_sum[ti][c] += rn->disp_spherical.T2[c];
            }

            disp_total += K_ring;
            total_contacts += contacts;
            total_pairs++;
        }

        ca.disp_shielding_contribution = SphericalTensor::Decompose(disp_total);
    }

    OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " vertex_contacts=" + std::to_string(total_contacts) +
        " bonded_exclusions=" + std::to_string(bonded_exclusions) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


SphericalTensor DispersionResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 K_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];

        // Ring-level distance check
        double ring_dist = (point - geom.center).norm();
        if (ring_dist < MIN_DISTANCE) continue;
        if (ring_dist < geom.radius) continue;
        if (ring_dist > RING_CALC_CUTOFF) continue;

        // Sum over ring vertices
        for (const auto& vertex : geom.vertices) {
            double r = (point - vertex).norm();
            if (r < MIN_DISTANCE || r > DISP_VERTEX_R_CUT) continue;

            auto vr = ComputeDispVertex(point, vertex, r);
            if (vr.valid) K_total += vr.K;
        }
    }

    return SphericalTensor::Decompose(K_total);
}


static void PackST_D(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int DispersionResult::WriteFeatures(const ProteinConformation& conf,
                                     const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * 8);
    std::vector<double> per_type_T2(N * 40);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_D(ca.disp_shielding_contribution, &shielding[i*9]);
        for (int t = 0; t < 8; ++t) {
            per_type_T0[i*8 + t] = ca.per_type_disp_scalar_sum[t];
            for (int c = 0; c < 5; ++c)
                per_type_T2[i*40 + t*5 + c] = ca.per_type_disp_T2_sum[t][c];
        }
    }

    NpyWriter::WriteFloat64(output_dir + "/disp_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/disp_per_type_T0.npy", per_type_T0.data(), N, 8);
    NpyWriter::WriteFloat64(output_dir + "/disp_per_type_T2.npy", per_type_T2.data(), N, 40);
    return 3;
}

}  // namespace nmr
```

## TARGET: src/MopacCoulombResult.cpp (instrument this file)

```cpp
#include "MopacCoulombResult.h"
#include "Protein.h"
#include "MopacResult.h"
#include "SpatialIndexResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <vector>

namespace nmr {


std::vector<std::type_index> MopacCoulombResult::Dependencies() const {
    return {
        std::type_index(typeid(MopacResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}


// ============================================================================
// MopacCoulombResult::Compute
//
// Same kernel as CoulombResult — same dipolar EFG, same Coulomb constant,
// same decomposition — but reading mopac_charge (PM7 Mulliken, conformation-
// dependent) instead of partial_charge (ff14SB, fixed per atom type).
//
// E_a(i) = ke * sum_{j!=i} q_mopac_j * (r_i - r_j)_a / |r_i - r_j|^3
// V_ab(i) = ke * sum_{j!=i} q_mopac_j * [3(r_i-r_j)_a(r_i-r_j)_b/|r_i-r_j|^5
//                                         - delta_ab / |r_i-r_j|^3]
// ============================================================================

std::unique_ptr<MopacCoulombResult> MopacCoulombResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("MopacCoulombResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t n_atoms = conf.AtomCount();

    auto result_ptr = std::make_unique<MopacCoulombResult>();
    result_ptr->conf_ = &conf;

    // ------------------------------------------------------------------
    // Atom classification: backbone, aromatic, sidechain.
    // Same topology walk as CoulombResult — from Residue backbone cache
    // and Ring atom membership. No EnrichmentResult dependency.
    // ------------------------------------------------------------------

    std::vector<bool> is_backbone(n_atoms, false);
    std::vector<bool> is_aromatic_atom(n_atoms, false);

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        auto mark_bb = [&](size_t idx) {
            if (idx != Residue::NONE && idx < n_atoms) is_backbone[idx] = true;
        };
        mark_bb(res.N);
        mark_bb(res.CA);
        mark_bb(res.C);
        mark_bb(res.O);
        mark_bb(res.H);
        mark_bb(res.HA);
        mark_bb(res.CB);
    }

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            if (ai < n_atoms) is_aromatic_atom[ai] = true;
        }
    }

    // ------------------------------------------------------------------
    // Primary bond direction for E_bond_proj (same as CoulombResult).
    // ------------------------------------------------------------------

    std::vector<Vec3> primary_bond_dir(n_atoms, Vec3::Zero());
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        if (atom.element == Element::H && atom.parent_atom_index != SIZE_MAX) {
            Vec3 d = conf.PositionAt(ai) - conf.PositionAt(atom.parent_atom_index);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        } else if (!atom.bond_indices.empty()) {
            const Bond& b = protein.BondAt(atom.bond_indices[0]);
            size_t other = (b.atom_index_a == ai) ? b.atom_index_b : b.atom_index_a;
            Vec3 d = conf.PositionAt(other) - conf.PositionAt(ai);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        }
    }

    // ------------------------------------------------------------------
    // Main N^2 Coulomb sum with MOPAC charges
    // ------------------------------------------------------------------

    KernelFilterSet filters;
    filters.Add(std::make_unique<SelfSourceFilter>());

    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos_i = conf.PositionAt(i);

        Vec3 E_total = Vec3::Zero();
        Vec3 E_backbone = Vec3::Zero();
        Vec3 E_sidechain = Vec3::Zero();
        Vec3 E_aromatic = Vec3::Zero();

        Mat3 EFG_total = Mat3::Zero();
        Mat3 EFG_backbone = Mat3::Zero();
        Mat3 EFG_sidechain = Mat3::Zero();
        Mat3 EFG_aromatic = Mat3::Zero();

        for (size_t j = 0; j < n_atoms; ++j) {
            KernelEvaluationContext ctx;
            ctx.atom_index = i;
            ctx.source_atom_a = j;
            ctx.distance = (pos_i - conf.PositionAt(j)).norm();
            if (!filters.AcceptAll(ctx)) continue;

            // MOPAC QM charge instead of ff14SB fixed charge
            double q_j = conf.AtomAt(j).mopac_charge;
            if (std::abs(q_j) < 1e-15) continue;

            Vec3 r = pos_i - conf.PositionAt(j);
            double r_mag = r.norm();

            if (r_mag < MIN_DISTANCE) continue;

            double r3 = r_mag * r_mag * r_mag;
            double r5 = r3 * r_mag * r_mag;

            // E_a = q_j * r_a / r^3
            Vec3 E_j = q_j * r / r3;

            // V_ab = q_j * (3 r_a r_b / r^5 - delta_ab / r^3)
            Mat3 V_j = q_j * (3.0 * r * r.transpose() / r5
                              - Mat3::Identity() / r3);

            E_total += E_j;
            EFG_total += V_j;

            if (is_aromatic_atom[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
            } else if (is_backbone[j]) {
                E_backbone += E_j;
                EFG_backbone += V_j;
            } else {
                E_sidechain += E_j;
                EFG_sidechain += V_j;
            }
        }

        // Apply Coulomb constant: convert from e/A^2 to V/A
        E_total     *= COULOMB_KE;
        E_backbone  *= COULOMB_KE;
        E_sidechain *= COULOMB_KE;
        E_aromatic  *= COULOMB_KE;
        EFG_total     *= COULOMB_KE;
        EFG_backbone  *= COULOMB_KE;
        EFG_sidechain *= COULOMB_KE;
        EFG_aromatic  *= COULOMB_KE;

        // Traceless projection: each term is traceless by Gauss's law,
        // but floating-point accumulation breaks this.
        auto project_traceless = [](Mat3& m) {
            m -= (m.trace() / 3.0) * Mat3::Identity();
        };
        project_traceless(EFG_total);
        project_traceless(EFG_backbone);
        project_traceless(EFG_sidechain);
        project_traceless(EFG_aromatic);

        // Sanitise NaN/Inf
        auto sanitise_vec = [](Vec3& v) {
            for (int d = 0; d < 3; ++d)
                if (std::isnan(v(d)) || std::isinf(v(d))) { v = Vec3::Zero(); return; }
        };
        auto sanitise_mat = [](Mat3& m) {
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    if (std::isnan(m(a,b)) || std::isinf(m(a,b))) m(a,b) = 0.0;
        };
        sanitise_vec(E_total);
        sanitise_vec(E_backbone);
        sanitise_vec(E_sidechain);
        sanitise_vec(E_aromatic);
        sanitise_mat(EFG_total);
        sanitise_mat(EFG_backbone);
        sanitise_mat(EFG_sidechain);
        sanitise_mat(EFG_aromatic);

        // Clamp extreme E-field magnitudes
        double E_mag = E_total.norm();
        if (E_mag > APBS_SANITY_LIMIT) {
            double scale = APBS_SANITY_LIMIT / E_mag;
            E_total     *= scale;
            E_backbone  *= scale;
            E_sidechain *= scale;
            E_aromatic  *= scale;
        }

        // ------------------------------------------------------------------
        // Store on ConformationAtom
        // ------------------------------------------------------------------
        auto& ca = conf.MutableAtomAt(i);

        ca.mopac_coulomb_E_total     = E_total;
        ca.mopac_coulomb_E_backbone  = E_backbone;
        ca.mopac_coulomb_E_sidechain = E_sidechain;
        ca.mopac_coulomb_E_aromatic  = E_aromatic;

        ca.mopac_coulomb_EFG_total   = EFG_total;
        ca.mopac_coulomb_EFG_total_spherical = SphericalTensor::Decompose(EFG_total);

        ca.mopac_coulomb_EFG_backbone = EFG_backbone;
        ca.mopac_coulomb_EFG_backbone_spherical = SphericalTensor::Decompose(EFG_backbone);

        ca.mopac_coulomb_EFG_aromatic = EFG_aromatic;
        ca.mopac_coulomb_EFG_aromatic_spherical = SphericalTensor::Decompose(EFG_aromatic);

        ca.mopac_coulomb_E_magnitude = E_total.norm();

        ca.mopac_coulomb_E_bond_proj = E_total.dot(primary_bond_dir[i]);

        if (ca.mopac_coulomb_E_magnitude > NEAR_ZERO_NORM) {
            Vec3 E_hat = E_total / ca.mopac_coulomb_E_magnitude;
            ca.mopac_coulomb_E_backbone_frac = E_backbone.dot(E_hat);
        } else {
            ca.mopac_coulomb_E_backbone_frac = 0.0;
        }

        // Shielding contribution: the total EFG SphericalTensor.
        // Pure T2 (EFG is traceless). gamma converts this to shielding.
        ca.mopac_coulomb_shielding_contribution =
            SphericalTensor::Decompose(EFG_total);
    }

    OperationLog::Info(LogCalcOther, "MopacCoulombResult::Compute",
        "atoms=" + std::to_string(n_atoms) +
        " rejected={" + filters.ReportRejections() + "}");

    return result_ptr;
}


// ============================================================================
// Query methods
// ============================================================================

Vec3 MopacCoulombResult::EFieldAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mopac_coulomb_E_total;
}

Mat3 MopacCoulombResult::EFGAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mopac_coulomb_EFG_total;
}

SphericalTensor MopacCoulombResult::EFGSphericalAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mopac_coulomb_EFG_total_spherical;
}


// ============================================================================
// WriteFeatures: mopac_coulomb_shielding (9), E-field (3),
// EFG decompositions, scalar features.
// ============================================================================

static void PackST_MCC(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int MopacCoulombResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> efield(N * 3);
    std::vector<double> efg_bb(N * 9);
    std::vector<double> efg_aro(N * 9);
    std::vector<double> scalars(N * 4);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_MCC(ca.mopac_coulomb_shielding_contribution, &shielding[i*9]);

        efield[i*3+0] = ca.mopac_coulomb_E_total.x();
        efield[i*3+1] = ca.mopac_coulomb_E_total.y();
        efield[i*3+2] = ca.mopac_coulomb_E_total.z();

        PackST_MCC(ca.mopac_coulomb_EFG_backbone_spherical, &efg_bb[i*9]);
        PackST_MCC(ca.mopac_coulomb_EFG_aromatic_spherical, &efg_aro[i*9]);

        scalars[i*4+0] = ca.mopac_coulomb_E_magnitude;
        scalars[i*4+1] = ca.mopac_coulomb_E_bond_proj;
        scalars[i*4+2] = ca.mopac_coulomb_E_backbone_frac;
        scalars[i*4+3] = ca.mopac_coulomb_E_aromatic.norm();
    }

    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_shielding.npy",
                            shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_E.npy",
                            efield.data(), N, 3);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg_backbone.npy",
                            efg_bb.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg_aromatic.npy",
                            efg_aro.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_scalars.npy",
                            scalars.data(), N, 4);
    return 5;
}

}  // namespace nmr
```

## TARGET: src/MopacMcConnellResult.cpp (instrument this file)

```cpp
#include "MopacMcConnellResult.h"
#include "Protein.h"
#include "MopacResult.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <algorithm>

namespace nmr {


std::vector<std::type_index> MopacMcConnellResult::Dependencies() const {
    return {
        std::type_index(typeid(MopacResult)),
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// Bond kernel computation (same physics as McConnellResult).
//
// M_ab = 9 cos_theta d_hat_a b_hat_b
//      - 3 b_hat_a b_hat_b
//      - (3 d_hat_a d_hat_b - delta_ab)
//
// Returns M_ab / r^3 (Angstrom^-3), plus symmetric traceless K and scalar f.
// ============================================================================

struct MopacBondKernelResult {
    Mat3 M_over_r3 = Mat3::Zero();
    Mat3 K = Mat3::Zero();
    double f = 0.0;
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();
};


static MopacBondKernelResult ComputeBondKernel(
        const Vec3& atom_pos,
        const Vec3& bond_midpoint,
        const Vec3& bond_direction) {

    MopacBondKernelResult result;

    Vec3 d = atom_pos - bond_midpoint;
    double r = d.norm();

    if (r < MIN_DISTANCE) return result;

    result.distance = r;

    double r3 = r * r * r;
    Vec3 d_hat = d / r;
    result.direction = d_hat;

    double cos_theta = d_hat.dot(bond_direction);

    // McConnell scalar: (3 cos^2 theta - 1) / r^3
    result.f = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    // Symmetric traceless dipolar kernel K_ab
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = (3.0 * d_hat(a) * d_hat(b)
                              - (a == b ? 1.0 : 0.0)) / r3;

    // Full McConnell tensor M_ab / r^3
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.M_over_r3(a, b) =
                (9.0 * cos_theta * d_hat(a) * bond_direction(b)
                 - 3.0 * bond_direction(a) * bond_direction(b)
                 - (3.0 * d_hat(a) * d_hat(b) - (a == b ? 1.0 : 0.0)))
                / r3;
        }
    }

    return result;
}


// ============================================================================
// MopacMcConnellResult::Compute
//
// Same loop as McConnellResult but each bond's contribution is weighted
// by its MOPAC Wiberg bond order. Bonds with no MOPAC order (0.0)
// contribute nothing — they are electronically insignificant.
// ============================================================================

std::unique_ptr<MopacMcConnellResult> MopacMcConnellResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("MopacMcConnellResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " bonds=" + std::to_string(conf.ProteinRef().BondCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const auto& mopac = conf.Result<MopacResult>();
    const size_t n_atoms = conf.AtomCount();

    auto result_ptr = std::make_unique<MopacMcConnellResult>();
    result_ptr->conf_ = &conf;

    KernelFilterSet filters;
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    int total_pairs = 0;
    int filtered_out = 0;
    int zero_bo_skipped = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_bonds = spatial.BondsWithinRadius(atom_pos, MOPAC_MCCONNELL_CUTOFF_A);

        // Per-category accumulators (bond-order-weighted)
        double co_sum = 0.0, cn_sum = 0.0, sidechain_sum = 0.0, aromatic_sum = 0.0;
        Mat3 M_backbone_total = Mat3::Zero();
        Mat3 M_sidechain_total = Mat3::Zero();
        Mat3 M_aromatic_total = Mat3::Zero();
        Mat3 M_total = Mat3::Zero();

        // Nearest CO and CN tracking (bond-order-weighted scalar)
        double best_co_dist = NO_DATA_SENTINEL;
        double best_cn_dist = NO_DATA_SENTINEL;
        double best_co_f_weighted = 0.0;
        MopacBondKernelResult best_co_kernel;
        MopacBondKernelResult best_cn_kernel;
        double best_co_bo = 0.0;
        double best_cn_bo = 0.0;

        for (size_t bi : nearby_bonds) {
            const Bond& bond = protein.BondAt(bi);

            // MOPAC Wiberg bond order for this topology bond
            double bo = mopac.TopologyBondOrder(bi);
            if (bo < 1e-6) { zero_bo_skipped++; continue; }

            Vec3 midpoint = conf.bond_midpoints[bi];
            Vec3 direction = conf.bond_directions[bi];

            MopacBondKernelResult kernel = ComputeBondKernel(atom_pos, midpoint, direction);
            if (kernel.distance < MIN_DISTANCE) continue;

            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = conf.bond_lengths[bi];
            ctx.atom_index = ai;
            ctx.source_atom_a = bond.atom_index_a;
            ctx.source_atom_b = bond.atom_index_b;
            if (!filters.AcceptAll(ctx)) { filtered_out++; continue; }

            // Bond-order-weighted accumulation
            Mat3 weighted_M = bo * kernel.M_over_r3;
            double weighted_f = bo * kernel.f;

            M_total += weighted_M;

            switch (bond.category) {
                case BondCategory::PeptideCO:
                    co_sum += weighted_f;
                    M_backbone_total += weighted_M;
                    if (kernel.distance < best_co_dist) {
                        best_co_dist = kernel.distance;
                        best_co_f_weighted = weighted_f;
                        best_co_kernel = kernel;
                        best_co_bo = bo;
                    }
                    break;

                case BondCategory::PeptideCN:
                    cn_sum += weighted_f;
                    M_backbone_total += weighted_M;
                    if (kernel.distance < best_cn_dist) {
                        best_cn_dist = kernel.distance;
                        best_cn_kernel = kernel;
                        best_cn_bo = bo;
                    }
                    break;

                case BondCategory::BackboneOther:
                    M_backbone_total += weighted_M;
                    break;

                case BondCategory::SidechainCO:
                    sidechain_sum += weighted_f;
                    M_sidechain_total += weighted_M;
                    break;

                case BondCategory::Aromatic:
                    aromatic_sum += weighted_f;
                    M_aromatic_total += weighted_M;
                    break;

                case BondCategory::SidechainOther:
                    M_sidechain_total += weighted_M;
                    break;

                default:
                    break;
            }

            total_pairs++;
        }

        // Store per-atom totals
        ca.mopac_mc_co_sum = co_sum;
        ca.mopac_mc_cn_sum = cn_sum;
        ca.mopac_mc_sidechain_sum = sidechain_sum;
        ca.mopac_mc_aromatic_sum = aromatic_sum;

        ca.mopac_mc_co_nearest = best_co_f_weighted;
        ca.mopac_mc_nearest_CO_dist = best_co_dist;
        ca.mopac_mc_nearest_CN_dist = best_cn_dist;

        if (best_co_dist < NO_DATA_SENTINEL) {
            ca.mopac_mc_T2_CO_nearest =
                SphericalTensor::Decompose(best_co_bo * best_co_kernel.K);
        }
        if (best_cn_dist < NO_DATA_SENTINEL) {
            ca.mopac_mc_T2_CN_nearest =
                SphericalTensor::Decompose(best_cn_bo * best_cn_kernel.K);
        }

        // Category T2 totals — extract symmetric traceless part
        Mat3 K_backbone = 0.5 * (M_backbone_total + M_backbone_total.transpose());
        K_backbone -= (K_backbone.trace() / 3.0) * Mat3::Identity();
        ca.mopac_mc_T2_backbone_total = SphericalTensor::Decompose(K_backbone);

        Mat3 K_sidechain = 0.5 * (M_sidechain_total + M_sidechain_total.transpose());
        K_sidechain -= (K_sidechain.trace() / 3.0) * Mat3::Identity();
        ca.mopac_mc_T2_sidechain_total = SphericalTensor::Decompose(K_sidechain);

        Mat3 K_aromatic = 0.5 * (M_aromatic_total + M_aromatic_total.transpose());
        K_aromatic -= (K_aromatic.trace() / 3.0) * Mat3::Identity();
        ca.mopac_mc_T2_aromatic_total = SphericalTensor::Decompose(K_aromatic);

        // Full shielding contribution
        ca.mopac_mc_shielding_contribution = SphericalTensor::Decompose(M_total);
    }

    OperationLog::Info(LogCalcMcConnell, "MopacMcConnellResult::Compute",
        "atom_bond_pairs=" + std::to_string(total_pairs) +
        " zero_bo_skipped=" + std::to_string(zero_bo_skipped) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms));

    return result_ptr;
}


// ============================================================================
// WriteFeatures: mopac_mc_shielding (9),
// mopac_mc_category_T2 (5 categories × 5 T2),
// mopac_mc_scalars (weighted CO/CN/sidechain/aromatic sums, nearest dists).
// ============================================================================

static void PackST_MMC(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int MopacMcConnellResult::WriteFeatures(const ProteinConformation& conf,
                                         const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> cat_T2(N * 25);
    std::vector<double> scalars(N * 6);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_MMC(ca.mopac_mc_shielding_contribution, &shielding[i*9]);

        const SphericalTensor* cats[5] = {
            &ca.mopac_mc_T2_backbone_total, &ca.mopac_mc_T2_sidechain_total,
            &ca.mopac_mc_T2_aromatic_total, &ca.mopac_mc_T2_CO_nearest,
            &ca.mopac_mc_T2_CN_nearest
        };
        for (int c = 0; c < 5; ++c)
            for (int m = 0; m < 5; ++m)
                cat_T2[i*25 + c*5 + m] = cats[c]->T2[m];

        scalars[i*6+0] = ca.mopac_mc_co_sum;
        scalars[i*6+1] = ca.mopac_mc_cn_sum;
        scalars[i*6+2] = ca.mopac_mc_sidechain_sum;
        scalars[i*6+3] = ca.mopac_mc_aromatic_sum;
        scalars[i*6+4] = ca.mopac_mc_nearest_CO_dist;
        scalars[i*6+5] = ca.mopac_mc_nearest_CN_dist;
    }

    NpyWriter::WriteFloat64(output_dir + "/mopac_mc_shielding.npy",
                            shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_mc_category_T2.npy",
                            cat_T2.data(), N, 25);
    NpyWriter::WriteFloat64(output_dir + "/mopac_mc_scalars.npy",
                            scalars.data(), N, 6);
    return 3;
}

}  // namespace nmr
```
