#include "PiQuadrupoleResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
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

struct PiQuadKernelResult {
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

    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;

    result.distance = r;
    result.direction = d / r;

    double r2 = r * r;
    double d_dot_n = d.dot(ring_normal);       // height above ring plane
    double cos_theta = d_dot_n / r;

    // Buckingham A-term kernel: (3 cos^2 theta - 1) / r^4  (feeds pq_per_type_T0)
    result.scalar = (3.0 * cos_theta * cos_theta - 1.0) / (r2 * r2);

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
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
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

            PiQuadKernelResult kernel = ComputePiQuadKernel(
                atom_pos, geom.center, geom.normal);

            // Apply filter: source extent = ring diameter
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = 2.0 * geom.radius;
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

            // shared ring-neighbour geometry (consumed by other ring calculators)
            RingNeighbourhood* rn = nullptr;
            for (auto& existing : conf_atom.ring_neighbours) {
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

                // ring-frame cylindrical coordinates (z along normal, rho in plane)
                Vec3 d_vec = atom_pos - geom.center;
                double z = d_vec.dot(geom.normal);
                Vec3 d_plane = d_vec - z * geom.normal;
                new_rn.z = z;
                new_rn.rho = d_plane.norm();
                new_rn.theta = std::atan2(d_plane.norm(), std::abs(z));  // theta folded to [0, pi/2] via |z| — quadrupole is symmetric across the ring plane

                conf_atom.ring_neighbours.push_back(new_rn);
                rn = &conf_atom.ring_neighbours.back();
            }

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
    return 1;
}

}  // namespace nmr
