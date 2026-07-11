#include "RingSusceptibilityResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "RingNeighbourGeometry.h"
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
// The clean ring-susceptibility scalar from one ring at one atom.
//
// Same derivation as McConnell (GEOMETRIC_KERNEL_CATALOGUE.md) with
// b_hat → n_hat (ring normal):
//
//   M_ab = 9 cosθ d̂_a n_b  -  3 n_a n_b  -  (3 d̂_a d̂_b - δ_ab)
//
// The former M_ab / r³ tensor family is intentionally no longer computed
// or emitted; this result now preserves only the clean scalar rescue.
// ============================================================================

namespace ring_susceptibility_detail {

KernelResult ComputeKernel(
        const Vec3& atom_pos,
        const Vec3& ring_center,
        const Vec3& ring_normal) {

    KernelResult result;

    // ring→atom displacement
    Vec3 ring_to_atom = atom_pos - ring_center;
    double r = ring_to_atom.norm();
    result.distance = r;

    if (!std::isfinite(r) ||
        r <= CalculatorConfig::Get("singularity_guard_distance")) {
        return result;
    }

    double r3 = r * r * r;
    Vec3 d_hat = ring_to_atom / r;
    result.direction = d_hat;

    double cos_theta = d_hat.dot(ring_normal);

    // Ring susceptibility scalar: (3 cos²θ - 1) / r³
    result.scalar = (3.0 * cos_theta * cos_theta - 1.0) / r3;
    result.valid = true;

    return result;
}

}  // namespace ring_susceptibility_detail


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

    // Point-center singularity is handled by the kernel itself.  The only
    // generic ring filter is the typed through-bond topology exclusion.
    KernelFilterSet filters;
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcOther, "RingSusceptibilityResult::Compute",
        "filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);

    int total_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& atom = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        // nearby rings (spatial index)
        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, CalculatorConfig::Get("ring_current_spatial_cutoff"));

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            ring_susceptibility_detail::KernelResult kernel =
                ring_susceptibility_detail::ComputeKernel(
                atom_pos, geom.center, geom.normal);

            if (!kernel.valid) {
                choices.Record(CalculatorId::RingSusceptibility, ri,
                    "center singularity exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source,
                                EntityOutcome::Included);
                        AddAtom(gc, &atom, ai, EntityRole::Target,
                                EntityOutcome::Excluded,
                                "center_singularity_guard");
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "distance_guard",
                                  CalculatorConfig::Get(
                                      "singularity_guard_distance"), "A");
                    });
                continue;
            }

            // Point-center kernels have no ring-diameter source extent.
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = 0.0;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) {
                // record exclusion
                choices.Record(CalculatorId::RingSusceptibility, ri, "filter exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &atom, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                continue;
            }

            RingNeighbourhood& neighbour_ref = EnsureRingNeighbourGeometry(
                atom, geom, ring, ri, atom_pos);
            RingNeighbourhood* neighbour = &neighbour_ref;

            // Store only the scalar rescue. The manufactured rank-2
            // ring-chi tensor family was removed.
            neighbour->chi_scalar = kernel.scalar;

            total_pairs++;
        }
    }

    OperationLog::Info(LogCalcOther, "RingSusceptibilityResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


int RingSusceptibilityResult::WriteFeatures(const ProteinConformation& conf,
                                             const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    size_t P = 0;
    for (size_t ai = 0; ai < N; ++ai) {
        P += conf.AtomAt(ai).ring_neighbours.size();
    }

    std::vector<double> sparse(P, 0.0);
    std::vector<double> per_type(
        N * static_cast<size_t>(kAromaticRingTypeCount), 0.0);

    size_t row = 0;
    for (size_t ai = 0; ai < N; ++ai) {
        for (const RingNeighbourhood& neighbour :
             conf.AtomAt(ai).ring_neighbours) {
            sparse[row++] = neighbour.chi_scalar;
            const int ti = static_cast<int>(neighbour.ring_type);
            if (ti >= 0 && ti < kAromaticRingTypeCount) {
                per_type[ai * kAromaticRingTypeCount +
                         static_cast<size_t>(ti)] += neighbour.chi_scalar;
            }
        }
    }

    NpyWriter::WriteFloat64(
        output_dir + "/ringchi_scalar.npy", sparse.data(), P);
    NpyWriter::WriteFloat64(
        output_dir + "/ringchi_per_type_T0.npy", per_type.data(), N,
        kAromaticRingTypeCount);
    return 2;
}

}  // namespace nmr
