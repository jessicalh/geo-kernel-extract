#include "QtHaighMallionCalc.h"

#include "QtPhysicalConstants.h"

#include <array>
#include <cmath>

namespace h5reader::calculators {

using model::Mat3;
using model::SphericalTensor;
using model::Vec3;

namespace {

// ---------------------------------------------------------------------------
// 7-point symmetric Gaussian quadrature on the reference triangle.
// Stroud T2:5-1 / Dunavant degree 5. Ported verbatim from
// HaighMallionResult.cpp::Gauss7.
// ---------------------------------------------------------------------------
struct TriQuadPoint {
    double lambda[3];
    double weight;
};

const std::array<TriQuadPoint, 7>& Gauss7() {
    static const double sqrt15 = std::sqrt(15.0);
    static const double a1 = (6.0 - sqrt15) / 21.0;           // ~0.1013
    static const double a2 = (6.0 + sqrt15) / 21.0;           // ~0.4701
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

// ---------------------------------------------------------------------------
// Accumulate the dipolar kernel integral over one triangle.
//   H_ab += int_triangle [ 3 rho_a rho_b / rho^5 - delta_ab / rho^3 ] dS
// where rho = r - r_s.
// ---------------------------------------------------------------------------
void AccumulateTensor(const Vec3& v0, const Vec3& v1, const Vec3& v2,
                       const Vec3& r,
                       const std::array<TriQuadPoint, 7>& qpts,
                       Mat3& H) {
    const double triArea = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    if (triArea < HM_TRIANGLE_AREA_GUARD) return;

    for (const auto& qp : qpts) {
        const Vec3 rS = qp.lambda[0] * v0 + qp.lambda[1] * v1 + qp.lambda[2] * v2;
        const Vec3 rho = r - rS;
        const double rhoMag = rho.norm();
        if (rhoMag < SINGULARITY_GUARD_DISTANCE) continue;

        const double rho3 = rhoMag * rhoMag * rhoMag;
        const double rho5 = rho3 * rhoMag * rhoMag;

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                H(a, b) += qp.weight * triArea *
                           (3.0 * rho(a) * rho(b) / rho5
                            - (a == b ? 1.0 : 0.0) / rho3);
            }
        }
    }
}

// Adaptive subdivision — if the field point is close to a triangle,
// split into 4 sub-triangles at edge midpoints. Max depth 2 (7 → 28 →
// 112 quadrature points per fan triangle).
bool NeedsSubdivision(const Vec3& v0, const Vec3& v1, const Vec3& v2,
                       const Vec3& r, double threshold) {
    return (r - v0).norm() < threshold
        || (r - v1).norm() < threshold
        || (r - v2).norm() < threshold;
}

void AccumulateAdaptive(const Vec3& v0, const Vec3& v1, const Vec3& v2,
                         const Vec3& r,
                         const std::array<TriQuadPoint, 7>& qpts,
                         Mat3& H, int level) {
    bool subdivide = false;
    if      (level == 0) subdivide = NeedsSubdivision(v0, v1, v2, r, HM_SUBDIV_THRESHOLD_L1);
    else if (level == 1) subdivide = NeedsSubdivision(v0, v1, v2, r, HM_SUBDIV_THRESHOLD_L2);

    if (subdivide && level < 2) {
        const Vec3 m01 = 0.5 * (v0 + v1);
        const Vec3 m12 = 0.5 * (v1 + v2);
        const Vec3 m02 = 0.5 * (v0 + v2);
        AccumulateAdaptive(v0,  m01, m02, r, qpts, H, level + 1);
        AccumulateAdaptive(m01, v1,  m12, r, qpts, H, level + 1);
        AccumulateAdaptive(m02, m12, v2,  r, qpts, H, level + 1);
        AccumulateAdaptive(m01, m12, m02, r, qpts, H, level + 1);
    } else {
        AccumulateTensor(v0, v1, v2, r, qpts, H);
    }
}

// Fan triangulation: n triangles from ring centroid to consecutive
// vertex pairs. Returns H_ab (symmetric traceless, units Å⁻¹).
Mat3 SurfaceIntegral(const Vec3& point,
                      const model::RingGeometry& geo,
                      const std::vector<Vec3>& vertices) {
    const int nv = static_cast<int>(vertices.size());
    if (nv < 3) return Mat3::Zero();
    const auto& qpts = Gauss7();
    Mat3 H = Mat3::Zero();
    for (int i = 0; i < nv; ++i) {
        const int j = (i + 1) % nv;
        AccumulateAdaptive(geo.center, vertices[i], vertices[j],
                            point, qpts, H, 0);
    }
    return H;
}

bool PointInValidRange(const Vec3& p, const model::RingGeometry& geo) {
    const double d = (p - geo.center).norm();
    return d >= SINGULARITY_GUARD_DISTANCE
        && d >= geo.radius
        && d <= RING_CURRENT_CUTOFF;
}

SphericalTensor DecomposeLocal(const Mat3& M) {
    SphericalTensor st;
    st.T0 = (M(0,0) + M(1,1) + M(2,2)) / 3.0;
    st.T1[0] = 0.5 * (M(0,1) - M(1,0));
    st.T1[1] = 0.5 * (M(0,2) - M(2,0));
    st.T1[2] = 0.5 * (M(1,2) - M(2,1));
    const double sxx = M(0,0) - st.T0;
    const double syy = M(1,1) - st.T0;
    const double sxy = 0.5 * (M(0,1) + M(1,0));
    const double sxz = 0.5 * (M(0,2) + M(2,0));
    const double syz = 0.5 * (M(1,2) + M(2,1));
    st.T2[0] = sxx - syy;
    st.T2[1] = 2.0 * sxy;
    st.T2[2] = 2.0 * sxz;
    st.T2[3] = 2.0 * syz;
    st.T2[4] = 2.0 * (M(2,2) - st.T0) - sxx - syy;
    return st;
}

}  // namespace

SphericalTensor EvaluateShielding(const Vec3&                   pointAng,
                                   const model::RingGeometry&   geo,
                                   const std::vector<Vec3>&     verticesAng,
                                   double                        intensityNA) {
    if (!PointInValidRange(pointAng, geo)) return {};
    if (verticesAng.size() < 3)             return {};

    const Mat3 H = SurfaceIntegral(pointAng, geo, verticesAng);
    const Vec3 V = H * geo.normal;

    // G_ab = -n_b * V_a, scaled by intensity (library applies the
    // ring current intensity to convert H units to shielding ppm).
    Mat3 G;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            G(a, b) = -geo.normal(b) * V(a) * intensityNA;
        }
    }
    return DecomposeLocal(G);
}

}  // namespace h5reader::calculators
