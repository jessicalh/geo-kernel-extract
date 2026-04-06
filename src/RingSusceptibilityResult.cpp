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
