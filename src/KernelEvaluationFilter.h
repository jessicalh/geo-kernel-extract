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
