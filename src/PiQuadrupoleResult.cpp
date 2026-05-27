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
// EFG geometric kernel from a point axial quadrupole at the ring center.
//
// Stone T-tensor derivation: V_ab = -(Theta/2) * T_abcd n_c n_d
// We compute G_ab = T_abcd n_c n_d (the -Theta/2 prefactor goes into
// the learnable parameter Q_type).
//
//   G_ab = 105 (d·n)^2 d_a d_b / r^9
//        - 30 (d·n) (n_a d_b + n_b d_a) / r^7
//        - 15 d_a d_b / r^7
//        + 6 n_a n_b / r^5
//        + delta_ab (3/r^5 - 15 (d·n)^2/r^7)
//
// Also: scalar = (3 cos^2 theta - 1) / r^4   (Buckingham A-term)
//
// See header for the Stone T-tensor derivation and properties.
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

    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;

    result.distance = r;
    result.direction = d / r;

    double r2 = r * r;
    double r5 = r2 * r2 * r;
    double r7 = r5 * r2;
    double r9 = r7 * r2;

    double d_dot_n = d.dot(ring_normal);       // height above ring plane
    double d_dot_n_sq = d_dot_n * d_dot_n;
    double cos_theta = d_dot_n / r;

    // Buckingham A-term kernel: (3 cos^2 theta - 1) / r^4  (feeds pq_per_type_T0)
    result.scalar = (3.0 * cos_theta * cos_theta - 1.0) / (r2 * r2);

    // EFG tensor G_ab = T_abcd n_c n_d (Stone Ch. 3), five terms in order
    // ((d·n) = d_dot_n, the height above the ring plane):
    //   105 (d·n)^2 d_a d_b / r^9
    //   - 30 (d·n) (n_a d_b + n_b d_a) / r^7
    //   - 15 d_a d_b / r^7
    //   + 6 n_a n_b / r^5
    //   + delta_ab (3/r^5 - 15 (d·n)^2/r^7)
    double delta_term = 3.0 / r5 - 15.0 * d_dot_n_sq / r7;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.G(a, b) =
                  105.0 * d_dot_n_sq * d(a) * d(b) / r9
                - 30.0 * d_dot_n * (ring_normal(a) * d(b) + ring_normal(b) * d(a)) / r7
                - 15.0 * d(a) * d(b) / r7
                + 6.0 * ring_normal(a) * ring_normal(b) / r5
                + (a == b ? delta_term : 0.0);
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

        Mat3 G_total = Mat3::Zero();

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

            // Store quadrupole kernel on RingNeighbourhood
            rn->quad_tensor = kernel.G;
            rn->quad_spherical = SphericalTensor::Decompose(kernel.G);
            rn->quad_scalar = kernel.scalar;

            // per-aromatic-ring-type feature bookkeeping (not part of the EFG sum)
            int ring_type_index = ring.TypeIndexAsInt();
            if (ring_type_index >= 0 && ring_type_index < kAromaticRingTypeCount) {
                conf_atom.per_type_pq_scalar_sum[ring_type_index] += kernel.scalar;
                for (int c = 0; c < 5; ++c)  // 5 T2 components
                    conf_atom.per_type_pq_T2_sum[ring_type_index][c] += rn->quad_spherical.T2[c];
            }

            // Accumulate EFG tensor
            G_total += kernel.G;
            accepted_pairs++;
        }

        // decompose the total EFG tensor (T0=T1=0; pure T2)
        conf_atom.piquad_shielding_contribution = SphericalTensor::Decompose(G_total);
    }

    OperationLog::Info(LogCalcOther, "PiQuadrupoleResult::Compute",
        "atom_ring_pairs=" + std::to_string(accepted_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


SphericalTensor PiQuadrupoleResult::SampleKernelAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 G_total = Mat3::Zero();

    // Same physics as Compute(); no atom-specific filters apply to a bare grid
    // point — only the geometric guards (singularity / inside-radius / cutoff).
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];

        double distance = (point - geom.center).norm();
        if (distance < CalculatorConfig::Get("singularity_guard_distance")) continue;
        if (distance < geom.radius) continue;
        if (distance > CalculatorConfig::Get("ring_current_spatial_cutoff")) continue;

        auto kernel = ComputePiQuadKernel(point, geom.center, geom.normal);
        G_total += kernel.G;
    }

    return SphericalTensor::Decompose(G_total);
}


int PiQuadrupoleResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    // NPY ABI shapes: 9 = packed SphericalTensor (T0, T1x3, T2x5);
    // 8 = aromatic ring types; 40 = 8x5 (per-type x T2).
    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * 8);
    std::vector<double> per_type_T2(N * 40);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        ca.piquad_shielding_contribution.PackFull9(&shielding[i*9]);
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
