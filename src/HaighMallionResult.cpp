#include "HaighMallionResult.h"
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

#include <algorithm>
#include <cmath>
#include <array>
#include <limits>
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
// Weights sum to 1.0; the physical triangle area enters separately via
// triArea in AccumulateTriangleIntegral.
// ============================================================================

struct TriQuadPoint {
    double lambda[3];
    double weight;
};

// 7-point symmetric Gaussian quadrature on the reference triangle.
// Stroud, A.H. Approximate Calculation of Multiple Integrals,
// Prentice-Hall (1971), rule T2:5-1 (degree 5).
// Also: Dunavant, D.A. Int. J. Numer. Meth. Engng. 21, 1129-1148 (1985).
static const std::array<TriQuadPoint, 7>& TriangleGauss7Rule() {
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

static void AccumulateTriangleIntegral(
        const Vec3& v0, const Vec3& v1, const Vec3& v2,
        const Vec3& field_point,
        const std::array<TriQuadPoint, 7>& qpts,
        Mat3& H) {

    double triArea = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    if (triArea < CalculatorConfig::Get("haigh_mallion_triangle_area_guard")) return;

    for (const auto& qp : qpts) {
        // Surface point in barycentric coordinates
        Vec3 surface_point = qp.lambda[0] * v0 + qp.lambda[1] * v1 + qp.lambda[2] * v2;
        Vec3 sep = field_point - surface_point;   // rho = r - r_s
        double rhoMag = sep.norm();
        if (rhoMag < CalculatorConfig::Get("singularity_guard_distance")) continue;

        double sep3 = rhoMag * rhoMag * rhoMag;
        double sep5 = sep3 * rhoMag * rhoMag;

        // K_ab = 3 rho_a rho_b / rho^5 - delta_ab / rho^3
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                H(a, b) += qp.weight * triArea *
                           (3.0 * sep(a) * sep(b) / sep5
                            - (a == b ? 1.0 : 0.0) / sep3);
    }
}


// ============================================================================
// Adaptive subdivision: when the field point is close to a triangle vertex,
// subdivide into 4 sub-triangles at edge midpoints for better accuracy.
//
// Level 0 -> 1: if any vertex is within the configured L1 threshold
// Level 1 -> 2: if any vertex is within the configured L2 threshold
// Max depth: 2 (7 -> 28 -> 112 quadrature points per fan triangle)
// ============================================================================

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


// ============================================================================
// Compute the HM surface integral for one ring at one atom position.
//
// Fan triangulation: n triangles from ring centroid to consecutive vertex pairs.
// Returns H_ab (symmetric, traceless, units Angstrom^-1).
// ============================================================================

namespace haigh_mallion_detail {

namespace {

double PointToSegmentDistance(const Vec3& point,
                              const Vec3& segment_a,
                              const Vec3& segment_b) {
    const Vec3 segment = segment_b - segment_a;
    const double length_sq = segment.squaredNorm();
    if (length_sq <= std::numeric_limits<double>::epsilon()) {
        return (point - segment_a).norm();
    }
    const double t = std::clamp(
        (point - segment_a).dot(segment) / length_sq, 0.0, 1.0);
    return (point - (segment_a + t * segment)).norm();
}

}  // namespace


double PointToTriangleDistance(const Vec3& point,
                               const Vec3& triangle_a,
                               const Vec3& triangle_b,
                               const Vec3& triangle_c) {
    const Vec3 ab = triangle_b - triangle_a;
    const Vec3 ac = triangle_c - triangle_a;

    // Degenerate fan triangles do not define a surface.  Their segment
    // boundary still gives a finite, conservative distance.
    if (ab.cross(ac).squaredNorm() <=
        std::numeric_limits<double>::epsilon()) {
        return std::min({
            PointToSegmentDistance(point, triangle_a, triangle_b),
            PointToSegmentDistance(point, triangle_b, triangle_c),
            PointToSegmentDistance(point, triangle_c, triangle_a)});
    }

    const Vec3 ap = point - triangle_a;
    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return ap.norm();

    const Vec3 bp = point - triangle_b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return bp.norm();

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return (point - (triangle_a + v * ab)).norm();
    }

    const Vec3 cp = point - triangle_c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) return cp.norm();

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return (point - (triangle_a + w * ac)).norm();
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const Vec3 bc = triangle_c - triangle_b;
        const double w = (d4 - d3) /
                         ((d4 - d3) + (d5 - d6));
        return (point - (triangle_b + w * bc)).norm();
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    const Vec3 closest = triangle_a + v * ab + w * ac;
    return (point - closest).norm();
}


double MinimumDistanceToFanSurface(const Vec3& point,
                                   const RingGeometry& geom) {
    if (geom.vertices.size() < 3) {
        return std::numeric_limits<double>::infinity();
    }

    double minimum = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < geom.vertices.size(); ++i) {
        const size_t j = (i + 1) % geom.vertices.size();
        minimum = std::min(minimum, PointToTriangleDistance(
            point, geom.center, geom.vertices[i], geom.vertices[j]));
    }
    return minimum;
}


Mat3 ComputeSurfaceIntegralH(
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


bool ContributesToCanonicalTotals(RingTypeIndex ring_type) {
    return ring_type != RingTypeIndex::TrpPerimeter;
}

}  // namespace haigh_mallion_detail


// ============================================================================
// HaighMallionResult::Compute
//
// For each atom, find rings within CalculatorConfig::Get("ring_current_spatial_cutoff"). For each ring:
//   1. Compute H_ab = surface integral of dipolar kernel (symmetric, traceless)
//   2. Compute V = H . n (effective B-field from magnetised surface)
//   3. Construct G_ab = -n_b * V_a (full shielding kernel, rank-1)
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

    // Topological exclusion is separate from the finite-surface singularity
    // guard applied against the actual fan triangles below.
    KernelFilterSet filters;
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

            // record this ring as a source (once)
            if (recorded_rings.find(ri) == recorded_rings.end()) {
                recorded_rings.insert(ri);
                choices.Record(CalculatorId::HaighMallion, ri, "surface integral",
                    [&ring](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                    });
            }

            double distance = (atom_pos - geom.center).norm();

            const double minimum_surface_distance =
                haigh_mallion_detail::MinimumDistanceToFanSurface(
                    atom_pos, geom);
            const double surface_guard =
                CalculatorConfig::Get("singularity_guard_distance");
            if (minimum_surface_distance <= surface_guard) {
                choices.Record(CalculatorId::HaighMallion, ri,
                    "surface singularity exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source,
                                EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target,
                                EntityOutcome::Excluded,
                                "fan_surface_distance_guard");
                        AddNumber(gc, "minimum_surface_distance",
                                  minimum_surface_distance, "A");
                        AddNumber(gc, "surface_distance_guard",
                                  surface_guard, "A");
                    });
                continue;
            }

            // Through-bond topology filter.  The surface kernel has no
            // center-distance source extent.
            KernelEvaluationContext ctx;
            ctx.distance = distance;
            ctx.source_extent = 0.0;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) {
                // ---- GeometryChoice: near-field exclusion ----
                choices.Record(CalculatorId::HaighMallion, ri, "topology exclusion",
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

            // --- physics ---
            // Step 1: Raw surface integral H_ab (symmetric, traceless, A^-1)
            Mat3 H = haigh_mallion_detail::ComputeSurfaceIntegralH(
                atom_pos, geom);

            // Step 2: Effective B-field V = H . n
            Vec3 effective_field = H * geom.normal;

            // Step 3: Full shielding kernel G_ab = -n_b * V_a (rank-1)
            // Minus sign from sigma_ab = -dB_a^sec / dB_{0,b}.
            // Same convention as BiotSavartResult: sigma = I * G gives
            // correct sign with literature I (negative for diamagnetic).
            // outer product -> rank-1 shielding kernel
            Mat3 G;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    G(a, b) = -geom.normal(b) * effective_field(a);

            RingNeighbourhood& rn_ref = EnsureRingNeighbourGeometry(
                ca, geom, ring, ri, atom_pos);
            RingNeighbourhood* rn = &rn_ref;

            // Store HM results on RingNeighbourhood
            // raw surface integral H
            rn->hm_H_tensor = H;                                  // raw integral (symmetric, traceless, pure T2)
            rn->hm_H_spherical = SphericalTensor::Decompose(H);   // Decompose(H): T0 ~= 0, T1 ~= 0 by construction (FP quadrature)
            // rank-1 shielding kernel G
            rn->hm_B_field = effective_field;                     // effective B-field V = H . n
            rn->hm_G_tensor = G;                                  // full shielding kernel (rank-1)
            rn->hm_G_spherical = SphericalTensor::Decompose(G);   // Decompose(G): T0, T1, T2

            if (haigh_mallion_detail::ContributesToCanonicalTotals(
                    ring.type_index)) {
                G_total += G;

                const int ti = ring.TypeIndexAsInt();
                if (ti >= 0 && ti < kAromaticRingTypeCount) {
                    ca.per_type_hm_T0_sum[ti] += rn->hm_G_spherical.T0;
                    for (int c = 0; c < 3; ++c)
                        ca.per_type_hm_T1_sum[ti][c] += rn->hm_G_spherical.T1[c];
                    for (int c = 0; c < 5; ++c)
                        ca.per_type_hm_T2_sum[ti][c] += rn->hm_G_spherical.T2[c];
                }
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
// SampleKernelAt: evaluate HM kernel at arbitrary 3D point.
// Same kernel as Compute(); grid-point exclusions differ (no bonded/atom
// filters; skips points within `radius` of the ring center, 3D distance).
// ============================================================================

SphericalTensor HaighMallionResult::SampleKernelAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 G_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const Ring& ring = protein.RingAt(ri);
        if (!haigh_mallion_detail::ContributesToCanonicalTotals(
                ring.type_index)) continue;
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        double distance = (point - geom.center).norm();
        if (distance > CalculatorConfig::Get("ring_current_spatial_cutoff")) continue;
        if (haigh_mallion_detail::MinimumDistanceToFanSurface(point, geom)
            <= CalculatorConfig::Get("singularity_guard_distance")) continue;

        Mat3 H = haigh_mallion_detail::ComputeSurfaceIntegralH(point, geom);
        Vec3 effective_field = H * geom.normal;

        // outer product -> rank-1 shielding kernel: G_ab = -n_b * V_a
        Mat3 G;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                G(a, b) = -geom.normal(b) * effective_field(a);

        G_total += G;
    }

    return SphericalTensor::Decompose(G_total);
}


// ============================================================================
// WriteFeatures: hm_shielding (9), per-type T0 (8), T1 (24), T2 (40).
// Mirrors BiotSavart layout — same ring-type decomposition, different kernel.
// ============================================================================

int HaighMallionResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * 8);
    std::vector<double> per_type_T1(N * 24);
    std::vector<double> per_type_T2(N * 40);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        ca.hm_shielding_contribution.PackFull9(&shielding[i*9]);
        for (int t = 0; t < 8; ++t) {
            per_type_T0[i*8 + t] = ca.per_type_hm_T0_sum[t];
            for (int c = 0; c < 3; ++c)
                per_type_T1[i*24 + t*3 + c] = ca.per_type_hm_T1_sum[t][c];
            for (int c = 0; c < 5; ++c)
                per_type_T2[i*40 + t*5 + c] = ca.per_type_hm_T2_sum[t][c];
        }
    }
    NpyWriter::WriteFloat64(output_dir + "/hm_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/hm_per_type_T0.npy", per_type_T0.data(), N, 8);
    NpyWriter::WriteFloat64(output_dir + "/hm_per_type_T1.npy", per_type_T1.data(), N, 24);
    NpyWriter::WriteFloat64(output_dir + "/hm_per_type_T2.npy", per_type_T2.data(), N, 40);
    return 4;
}

}  // namespace nmr
