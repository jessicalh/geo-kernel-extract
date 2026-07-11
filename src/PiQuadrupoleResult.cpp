#include "PiQuadrupoleResult.h"
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


std::vector<std::type_index> PiQuadrupoleResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// Clean Buckingham A-term scalar from a point axial quadrupole at the
// ring center. The former geometric EFG tensor G_ab is intentionally no
// longer computed or emitted because its physical prefactor was deferred
// to a learnable parameter.
// ============================================================================

namespace pi_quadrupole_detail {

KernelResult ComputeKernel(
        const Vec3& atom_pos,
        const Vec3& ring_center,
        const Vec3& ring_normal) {

    KernelResult result;

    Vec3 d = atom_pos - ring_center;
    double r = d.norm();
    result.distance = r;

    if (!std::isfinite(r) ||
        r <= CalculatorConfig::Get("singularity_guard_distance")) {
        return result;
    }

    result.direction = d / r;

    double r2 = r * r;
    double d_dot_n = d.dot(ring_normal);       // height above ring plane
    double cos_theta = d_dot_n / r;

    // Buckingham A-term kernel: (3 cos^2 theta - 1) / r^4  (feeds pq_per_type_T0)
    result.scalar = (3.0 * cos_theta * cos_theta - 1.0) / (r2 * r2);
    result.valid = true;

    return result;
}

}  // namespace pi_quadrupole_detail


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

    // Point-center singularity is handled by the exact kernel distance;
    // topology separately excludes ring vertices/bonded neighbours.  A ring
    // diameter is not a valid exclusion sphere for this point kernel.
    KernelFilterSet filters;
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcOther, "PiQuadrupoleResult::Compute",
        "filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);

    int accepted_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& conf_atom = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, CalculatorConfig::Get("ring_current_spatial_cutoff"));

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            pi_quadrupole_detail::KernelResult kernel =
                pi_quadrupole_detail::ComputeKernel(
                atom_pos, geom.center, geom.normal);

            if (!kernel.valid) {
                choices.Record(CalculatorId::PiQuadrupole, ri,
                    "center singularity exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source,
                                EntityOutcome::Included);
                        AddAtom(gc, &conf_atom, ai, EntityRole::Target,
                                EntityOutcome::Excluded,
                                "center_singularity_guard");
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "distance_guard",
                                  CalculatorConfig::Get(
                                      "singularity_guard_distance"), "A");
                    });
                continue;
            }

            // Topological exclusion only; point-center kernels have no
            // source diameter in the filter context.
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = 0.0;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) {
                // record the exclusion
                choices.Record(CalculatorId::PiQuadrupole, ri, "filter exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &conf_atom, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                continue;
            }

            RingNeighbourhood& rn_ref = EnsureRingNeighbourGeometry(
                conf_atom, geom, ring, ri, atom_pos);
            RingNeighbourhood* rn = &rn_ref;

            // Store only the scalar rescue.
            rn->quad_scalar = kernel.scalar;

            // per-aromatic-ring-type feature bookkeeping (not part of the EFG sum)
            int ring_type_index = ring.TypeIndexAsInt();
            if (ring_type_index >= 0 && ring_type_index < kAromaticRingTypeCount) {
                conf_atom.per_type_pq_scalar_sum[ring_type_index] += kernel.scalar;
            }

            accepted_pairs++;
        }
    }

    OperationLog::Info(LogCalcOther, "PiQuadrupoleResult::Compute",
        "atom_ring_pairs=" + std::to_string(accepted_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


int PiQuadrupoleResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> per_type_T0(N * 8);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        for (int t = 0; t < 8; ++t) {
            per_type_T0[i*8 + t] = ca.per_type_pq_scalar_sum[t];
        }
    }

    NpyWriter::WriteFloat64(output_dir + "/pq_per_type_T0.npy", per_type_T0.data(), N, 8);
    NpyWriter::WriteFloat64(
        output_dir + "/piquad_axial_scalar_per_type_T0.npy",
        per_type_T0.data(), N, 8);
    return 2;
}

}  // namespace nmr
