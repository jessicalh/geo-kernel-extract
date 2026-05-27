#include "RingSusceptibilityResult.h"
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


std::vector<std::type_index> RingSusceptibilityResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// The full ring-susceptibility tensor kernel from one ring at one atom.
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
//   Term 1: 9 cosθ d̂ ⊗ n̂      — generally asymmetric; can contribute T1
//   Term 2: -3 n̂ ⊗ n̂           — symmetric; contributes T0 and T2
//   Term 3: -(3 d̂ ⊗ d̂ - I)    — symmetric traceless; contributes T2
// ============================================================================

struct RingChiKernelResult {
    Mat3 full_tensor_over_r3 = Mat3::Zero();  // full tensor M_ab / r³ (asymmetric)
    Mat3 dipolar_kernel = Mat3::Zero();       // symmetric traceless dipolar kernel K_ab
    double scalar_kernel = 0.0;               // ring susceptibility scalar f
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();   // unit vector from ring center to atom
};


static RingChiKernelResult ComputeRingChiKernel(
        const Vec3& atom_pos,
        const Vec3& ring_center,
        const Vec3& ring_normal) {

    RingChiKernelResult result;

    // ring→atom displacement
    Vec3 ring_to_atom = atom_pos - ring_center;
    double r = ring_to_atom.norm();

    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;

    result.distance = r;

    double r3 = r * r * r;
    Vec3 d_hat = ring_to_atom / r;
    result.direction = d_hat;

    double cos_theta = d_hat.dot(ring_normal);

    // Ring susceptibility scalar: (3 cos²θ - 1) / r³
    result.scalar_kernel = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    // Symmetric traceless dipolar kernel K_ab
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.dipolar_kernel(a, b) = (3.0 * d_hat(a) * d_hat(b)
                              - (a == b ? 1.0 : 0.0)) / r3;  // δ_ab

    // Full tensor M_ab / r³
    //   = [9 cosθ d̂_a n_b - 3 n_a n_b - (3 d̂_a d̂_b - δ_ab)] / r³
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.full_tensor_over_r3(a, b) =
                (9.0 * cos_theta * d_hat(a) * ring_normal(b)
                 - 3.0 * ring_normal(a) * ring_normal(b)
                 - (3.0 * d_hat(a) * d_hat(b) - (a == b ? 1.0 : 0.0)))  // δ_ab
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
    // close ring atoms by geometry, and the topology check excludes ring
    // vertices and their bonded neighbours directly.
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
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

        Mat3 M_total = Mat3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            RingChiKernelResult kernel = ComputeRingChiKernel(
                atom_pos, geom.center, geom.normal);

            // filter pair: source extent = ring diameter (2 * radius)
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = 2.0 * geom.radius;
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

            // Attach per-ring feature record to this atom (not kernel math).
            // The ring_neighbours vector may already have an entry from
            // another calculator (BiotSavart). Find by ring index.
            RingNeighbourhood* neighbour = nullptr;
            for (auto& existing : atom.ring_neighbours) {
                if (existing.ring_index == ri) {
                    neighbour = &existing;
                    break;
                }
            }
            if (!neighbour) {
                RingNeighbourhood new_neighbour;
                new_neighbour.ring_index = ri;
                new_neighbour.ring_type = ring.type_index;
                new_neighbour.distance_to_center = kernel.distance;
                // contract: unit vector pointing center→atom (see RingChiKernelResult::direction)
                new_neighbour.direction_to_center = kernel.direction;

                // Ring-frame coordinates: z along normal, rho in-plane,
                // theta = polar angle from normal (hemisphere-folded).
                Vec3 atom_offset = atom_pos - geom.center;
                double z = atom_offset.dot(geom.normal);
                Vec3 d_plane = atom_offset - z * geom.normal;
                double rho = d_plane.norm();
                double theta = std::atan2(rho, std::abs(z));
                new_neighbour.z = z;
                new_neighbour.rho = rho;
                new_neighbour.theta = theta;

                atom.ring_neighbours.push_back(new_neighbour);
                neighbour = &atom.ring_neighbours.back();
            }

            // store ring kernel
            // per-neighbour record (overwritten per ring); per-atom sum is M_total below
            neighbour->chi_tensor = kernel.full_tensor_over_r3;
            neighbour->chi_spherical = SphericalTensor::Decompose(kernel.full_tensor_over_r3);
            neighbour->chi_scalar = kernel.scalar_kernel;

            // accumulate atom tensor
            M_total += kernel.full_tensor_over_r3;
            total_pairs++;
        }

        // decompose atom total
        atom.ringchi_shielding_contribution = SphericalTensor::Decompose(M_total);
    }

    OperationLog::Info(LogCalcOther, "RingSusceptibilityResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


SphericalTensor RingSusceptibilityResult::SampleKernelAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 M_total = Mat3::Zero();

    // Grid-point sampling: same kernel as Compute, evaluated at an arbitrary
    // point with no per-atom filter set.
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];

        double distance = (point - geom.center).norm();
        if (distance < CalculatorConfig::Get("singularity_guard_distance")) continue;  // singularity
        if (distance < geom.radius) continue;                                          // inside-ring
        if (distance > CalculatorConfig::Get("ring_current_spatial_cutoff")) continue; // cutoff

        auto kernel = ComputeRingChiKernel(point, geom.center, geom.normal);
        M_total += kernel.full_tensor_over_r3;
    }

    return SphericalTensor::Decompose(M_total);
}


int RingSusceptibilityResult::WriteFeatures(const ProteinConformation& conf,
                                             const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    std::vector<double> shielding(N * 9);
    for (size_t i = 0; i < N; ++i)
        conf.AtomAt(i).ringchi_shielding_contribution.PackFull9(&shielding[i*9]);
    NpyWriter::WriteFloat64(output_dir + "/ringchi_shielding.npy", shielding.data(), N, 9);
    return 1;
}

}  // namespace nmr
