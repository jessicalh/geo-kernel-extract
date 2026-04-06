#include "HaighMallionResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <array>

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
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcHaighMal, "HaighMallionResult::Compute",
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

            double distance = (atom_pos - geom.center).norm();
            if (distance < MIN_DISTANCE) continue;

            // Apply filter
            KernelEvaluationContext ctx;
            ctx.distance = distance;
            ctx.source_extent = 2.0 * geom.radius;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) continue;

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
