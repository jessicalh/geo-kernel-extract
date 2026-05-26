#include "HaighMallionResult.h"
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
#include <array>
#include <set>

namespace nmr {


std::vector<std::type_index> HaighMallionResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}

// Barycentric coordinates (lambda_0, lambda_1, lambda_2) and weights.
// Three orbits: centroid (1 point), near-vertex (3 points), near-edge (3 points).
// Weights sum to 1.0; the physical triangle area enters separately via
// triArea in AccumulateTriangleIntegral.

struct TriQuadPoint {
    double lambda[3];
    double weight;
};

// 7-point symmetric Gaussian quadrature on the reference triangle.
static const std::array<TriQuadPoint, 7>& TriangleGauss7Rule() {
    static const double sqrt15 = std::sqrt(15.0);
    static const double a1 = (6.0 - sqrt15) / 21.0;          // ~0.1013
    static const double a2 = (6.0 + sqrt15) / 21.0;          // ~0.4701
    static const double w0 = 9.0 / 40.0;                      //  0.225
    static const double w1 = (155.0 - sqrt15) / 1200.0;       // ~0.1259
    static const double w2 = (155.0 + sqrt15) / 1200.0;       // ~0.1324

    static const std::array<TriQuadPoint, 7> pts = {{
        {{ 1.0/3.0, 1.0/3.0, 1.0/3.0 }, w0},
        {{ a1, a1, 1.0 - 2.0*a1 }, w1},
        {{ a1, 1.0 - 2.0*a1, a1 }, w1},
        {{ 1.0 - 2.0*a1, a1, a1 }, w1},
        {{ a2, a2, 1.0 - 2.0*a2 }, w2},
        {{ a2, 1.0 - 2.0*a2, a2 }, w2},
        {{ 1.0 - 2.0*a2, a2, a2 }, w2},
    }};
    return pts;
}

// Accumulate the dipolar kernel integral over one triangle.
// H_ab += integral_triangle [ 3 rho_a rho_b / rho^5 - delta_ab / rho^3 ] dS
// where rho = r - r_s (field point minus surface point).

static void AccumulateTriangleIntegral(
        const Vec3& v0, const Vec3& v1, const Vec3& v2,
        const Vec3& field_point,
        const std::array<TriQuadPoint, 7>& qpts,
        Mat3& H) {

    double triArea = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    if (triArea < CalculatorConfig::Get("haigh_mallion_triangle_area_guard")) return;

    for (const auto& qp : qpts) {
        Vec3 surface_point = qp.lambda[0] * v0 + qp.lambda[1] * v1 + qp.lambda[2] * v2;
        Vec3 sep = field_point - surface_point;   // rho = r - r_s
        double rhoMag = sep.norm();
        if (rhoMag < CalculatorConfig::Get("singularity_guard_distance")) continue;

        double sep3 = rhoMag * rhoMag * rhoMag;
        double sep5 = sep3 * rhoMag * rhoMag;

        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                H(a, b) += qp.weight * triArea *
                           (3.0 * sep(a) * sep(b) / sep5
                            - (a == b ? 1.0 : 0.0) / sep3);
    }
}

// Adaptive subdivision: when the field point is close to a triangle vertex,
// subdivide into four sub-triangles at edge midpoints. Thresholds come from
// CalculatorConfig; max depth is 2.

static bool NeedsSubdivision(
        const Vec3& v0, const Vec3& v1, const Vec3& v2,
        const Vec3& field_point, double threshold) {
    return (field_point - v0).norm() < threshold
        || (field_point - v1).norm() < threshold
        || (field_point - v2).norm() < threshold;
}

static void AccumulateAdaptiveTriangleIntegral(
        const Vec3& v0, const Vec3& v1, const Vec3& v2,
        const Vec3& field_point,
        const std::array<TriQuadPoint, 7>& qpts,
        Mat3& H, int level) {

    bool subdivide = false;
    if (level == 0)
        subdivide = NeedsSubdivision(v0, v1, v2, field_point, CalculatorConfig::Get("haigh_mallion_subdivision_threshold_l1"));
    else if (level == 1)
        subdivide = NeedsSubdivision(v0, v1, v2, field_point, CalculatorConfig::Get("haigh_mallion_subdivision_threshold_l2"));

    if (subdivide && level < 2) {
        Vec3 m01 = 0.5 * (v0 + v1);
        Vec3 m12 = 0.5 * (v1 + v2);
        Vec3 m02 = 0.5 * (v0 + v2);
        AccumulateAdaptiveTriangleIntegral(v0,  m01, m02, field_point, qpts, H, level + 1);
        AccumulateAdaptiveTriangleIntegral(m01, v1,  m12, field_point, qpts, H, level + 1);
        AccumulateAdaptiveTriangleIntegral(m02, m12, v2,  field_point, qpts, H, level + 1);
        AccumulateAdaptiveTriangleIntegral(m01, m12, m02, field_point, qpts, H, level + 1);
    } else {
        AccumulateTriangleIntegral(v0, v1, v2, field_point, qpts, H);
    }
}

// Compute the HM surface integral for one ring at one atom position.
// Fan triangulation: n triangles from ring centroid to consecutive vertex pairs.
// Returns H_ab (symmetric, traceless, units Angstrom^-1).

static Mat3 ComputeSurfaceIntegralH(
        const Vec3& point,
        const RingGeometry& geom) {

    const auto& verts = geom.vertices;
    int nv = static_cast<int>(verts.size());
    if (nv < 3) return Mat3::Zero();

    const auto& qpts = TriangleGauss7Rule();
    Mat3 H = Mat3::Zero();

    for (int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        AccumulateAdaptiveTriangleIntegral(geom.center, verts[i], verts[j],
                                           point, qpts, H, 0);
    }

    return H;
}

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

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, CalculatorConfig::Get("ring_current_spatial_cutoff"));

        Mat3 G_total = Mat3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            // Record this ring as a source once.
            if (recorded_rings.find(ri) == recorded_rings.end()) {
                recorded_rings.insert(ri);
                choices.Record(CalculatorId::HaighMallion, ri, "surface integral",
                    [&ring](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                    });
            }

            double distance = (atom_pos - geom.center).norm();

            KernelEvaluationContext ctx;
            ctx.distance = distance;
            ctx.source_extent = 2.0 * geom.radius;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) {
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

            if (distance < CalculatorConfig::Get("haigh_mallion_subdivision_threshold_l1")) {
                choices.Record(CalculatorId::HaighMallion, ri, "adaptive refinement",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Triggered);
                        AddNumber(gc, "distance", distance, "A");
                        AddNumber(gc, "L1_threshold", CalculatorConfig::Get("haigh_mallion_subdivision_threshold_l1"), "A");
                        AddNumber(gc, "L2_threshold", CalculatorConfig::Get("haigh_mallion_subdivision_threshold_l2"), "A");
                    });
            }

            Mat3 H = ComputeSurfaceIntegralH(atom_pos, geom);

            Vec3 effective_field = H * geom.normal;

            // G_ab = -n_b * V_a.
            Mat3 G;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    G(a, b) = -geom.normal(b) * effective_field(a);

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

                Vec3 atom_from_center = atom_pos - geom.center;
                new_rn.direction_to_center = atom_from_center.normalized();

                double z = atom_from_center.dot(geom.normal);
                Vec3 d_plane = atom_from_center - z * geom.normal;
                double in_plane_radius = d_plane.norm();
                double theta = std::atan2(in_plane_radius, std::abs(z));
                new_rn.z = z;
                new_rn.rho = in_plane_radius;
                new_rn.theta = theta;

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            rn->hm_H_tensor = H;                                  // raw integral (symmetric, traceless, pure T2)
            rn->hm_H_spherical = SphericalTensor::Decompose(H);   // Decompose(H): T0 ~= 0, T1 ~= 0 by construction (FP quadrature)
            rn->hm_B_field = effective_field;                     // effective B-field V = H . n
            rn->hm_G_tensor = G;                                  // full shielding kernel (rank-1)
            rn->hm_G_spherical = SphericalTensor::Decompose(G);   // Decompose(G): T0, T1, T2

            G_total += G;

            int ti = ring.TypeIndexAsInt();
            // Hard-coded output width: aromatic ring types occupy indices 0..7.
            if (ti >= 0 && ti < 8) {
                ca.per_type_hm_T0_sum[ti] += rn->hm_G_spherical.T0;
                for (int c = 0; c < 5; ++c)
                    ca.per_type_hm_T2_sum[ti][c] += rn->hm_G_spherical.T2[c];
            }

            total_pairs++;
        }

        ca.hm_shielding_contribution = SphericalTensor::Decompose(G_total);
    }

    OperationLog::Info(LogCalcHaighMal, "HaighMallionResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}

// Same kernel as Compute(); grid-point exclusions differ (no bonded/atom
// filters; skips points within `radius` of the ring center, 3D distance).

SphericalTensor HaighMallionResult::SampleKernelAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 G_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        double distance = (point - geom.center).norm();
        if (distance < CalculatorConfig::Get("singularity_guard_distance")) continue;
        if (distance < geom.radius) continue;
        if (distance > CalculatorConfig::Get("ring_current_spatial_cutoff")) continue;

        Mat3 H = ComputeSurfaceIntegralH(point, geom);
        Vec3 effective_field = H * geom.normal;

        Mat3 G;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                G(a, b) = -geom.normal(b) * effective_field(a);

        G_total += G;
    }

    return SphericalTensor::Decompose(G_total);
}

// WriteFeatures: hm_shielding (9), per-type T0 (8), per-type T2 (40).

int HaighMallionResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * 8);
    std::vector<double> per_type_T2(N * 40);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        ca.hm_shielding_contribution.PackFull9(&shielding[i*9]);
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
