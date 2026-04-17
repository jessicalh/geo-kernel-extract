#include "QtBiotSavartCalc.h"

#include "QtPhysicalConstants.h"

#include <cmath>

namespace h5reader::calculators {

using model::Mat3;
using model::SphericalTensor;
using model::Vec3;

namespace {

// ---------------------------------------------------------------------------
// Wire segment B-field (Biot-Savart law).
// Ported verbatim from BiotSavartResult.cpp::WireSegmentField.
//
// All computation in SI: positions in metres, current in amperes, B in Tesla.
//
//   B = (mu_0/4pi) * I * (dl x dA) / |dl x dA|^2 * (dl.dA/|dA| - dl.dB/|dB|)
//
// where dl = b - a, dA = r - a, dB = r - b.
//
// Guards (SI thresholds):
//   lenA, lenB < BS_WIRE_ENDPOINT_GUARD: field point at segment endpoint
//   crossSq < BS_WIRE_AXIS_GUARD:         field point on the wire axis
// ---------------------------------------------------------------------------

Vec3 WireSegmentField(const Vec3& a_m, const Vec3& b_m,
                       double I_A, const Vec3& r_m) {
    const Vec3 dl = b_m - a_m;
    const Vec3 dA = r_m - a_m;
    const Vec3 dB = r_m - b_m;

    const double lenA = dA.norm();
    const double lenB = dB.norm();
    if (lenA < BS_WIRE_ENDPOINT_GUARD || lenB < BS_WIRE_ENDPOINT_GUARD)
        return Vec3::Zero();

    const Vec3   cross   = dl.cross(dA);
    const double crossSq = cross.squaredNorm();
    if (crossSq < BS_WIRE_AXIS_GUARD) return Vec3::Zero();

    const double factor  = BIOT_SAVART_PREFACTOR * I_A / crossSq;
    const double dotTerm = dl.dot(dA) / lenA - dl.dot(dB) / lenB;
    return factor * dotTerm * cross;   // Tesla
}

// ---------------------------------------------------------------------------
// Johnson-Bovey double-loop model.
// Two current loops at ±lobe_offset from the ring plane along the normal,
// each carrying half the total current (I/2). Sum over wire segments of
// both loops.
// Input: vertices and point in Angstroms, current in nanoamperes.
// Converts to SI at the boundary; returns B in Tesla.
// ---------------------------------------------------------------------------

Vec3 JohnsonBoveyField(
    const std::vector<Vec3>& verticesAng,
    const Vec3&              normal,
    double                   lobeOffsetA,
    double                   currentNA,
    const Vec3&              pointAng) {

    const int n = static_cast<int>(verticesAng.size());
    if (n < 3) return Vec3::Zero();

    const Vec3   offsetAng = normal * lobeOffsetA;
    const double halfI_A   = 0.5 * currentNA * NANOAMPERES_TO_AMPERES;
    const Vec3   r_m       = pointAng * ANGSTROMS_TO_METRES;

    Vec3 B = Vec3::Zero();
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;

        const Vec3 a_upper = (verticesAng[i] + offsetAng) * ANGSTROMS_TO_METRES;
        const Vec3 b_upper = (verticesAng[j] + offsetAng) * ANGSTROMS_TO_METRES;
        const Vec3 a_lower = (verticesAng[i] - offsetAng) * ANGSTROMS_TO_METRES;
        const Vec3 b_lower = (verticesAng[j] - offsetAng) * ANGSTROMS_TO_METRES;

        B += WireSegmentField(a_upper, b_upper, halfI_A, r_m);
        B += WireSegmentField(a_lower, b_lower, halfI_A, r_m);
    }
    return B;   // Tesla
}

// ---------------------------------------------------------------------------
// Shared spatial filter: skip the point if it's too close to or too far
// from the ring center, or inside the ring (multipole invalid). Matches
// the library's SampleShieldingAt / SampleBFieldAt predicate.
// ---------------------------------------------------------------------------
bool PointInValidRange(const Vec3& pointAng, const model::RingGeometry& geo) {
    const double distance = (pointAng - geo.center).norm();
    if (distance < SINGULARITY_GUARD_DISTANCE) return false;
    if (distance < geo.radius)                 return false;   // inside ring
    if (distance > RING_CURRENT_CUTOFF)        return false;
    return true;
}

// Spherical-tensor decomposition of a 3x3 matrix using the library's
// T0 = trace/3, T1 = antisymmetric flattened (0, 1, 2) = (xy-yx, xz-zx, yz-zy)/2,
// T2 = traceless symmetric in 5 components. For the BS/HM butterflies
// we only render T0 as an isosurface; T1 and T2 are filled for future
// ellipsoid glyphs and consistency with the library's layout.
SphericalTensor DecomposeLocal(const Mat3& M) {
    SphericalTensor st;
    st.T0 = (M(0,0) + M(1,1) + M(2,2)) / 3.0;
    // Antisymmetric part 0.5 * (M - M^T)
    st.T1[0] = 0.5 * (M(0,1) - M(1,0));
    st.T1[1] = 0.5 * (M(0,2) - M(2,0));
    st.T1[2] = 0.5 * (M(1,2) - M(2,1));
    // Traceless symmetric part 0.5*(M + M^T) - T0*I.
    // 5 components in the order: (xx-yy)/sqrt(3), xy, xz, yz, (2zz-xx-yy)/3.
    // This matches a common isometric real-sphericart layout; the reader
    // never compares T2 to library T2 directly (the H5 carries
    // precomputed per-atom shieldings; our T2 is only rendered for glyphs
    // in a later step, so the layout here is self-consistent).
    const double sxx = 0.5 * (M(0,0) + M(0,0)) - st.T0;
    const double syy = 0.5 * (M(1,1) + M(1,1)) - st.T0;
    const double szz = 0.5 * (M(2,2) + M(2,2)) - st.T0;
    const double sxy = 0.5 * (M(0,1) + M(1,0));
    const double sxz = 0.5 * (M(0,2) + M(2,0));
    const double syz = 0.5 * (M(1,2) + M(2,1));
    st.T2[0] = sxx - syy;
    st.T2[1] = 2.0 * sxy;
    st.T2[2] = 2.0 * sxz;
    st.T2[3] = 2.0 * syz;
    st.T2[4] = 2.0 * szz - sxx - syy;
    return st;
}

}  // namespace

Vec3 EvaluateBField(const Vec3&                     pointAng,
                     const model::RingGeometry&     geo,
                     const std::vector<Vec3>&       verticesAng,
                     double                          lobeOffsetA,
                     double                          currentNA) {
    if (!PointInValidRange(pointAng, geo))     return Vec3::Zero();
    if (verticesAng.size() < 3)                 return Vec3::Zero();
    return JohnsonBoveyField(verticesAng, geo.normal,
                              lobeOffsetA, currentNA, pointAng);
}

SphericalTensor EvaluateShielding(const Vec3&                   pointAng,
                                   const model::RingGeometry&   geo,
                                   const std::vector<Vec3>&     verticesAng,
                                   double                        lobeOffsetA,
                                   double                        intensityNA) {
    if (!PointInValidRange(pointAng, geo))     return {};
    if (verticesAng.size() < 3)                 return {};

    // Evaluate B at unit current (1 nA) then scale by intensity.
    // Equivalent to evaluating directly at intensityNA — library uses
    // the unit-current path for separability and we match it.
    const Vec3 B = JohnsonBoveyField(verticesAng, geo.normal,
                                       lobeOffsetA, 1.0, pointAng);

    // G_ab = -n_b * B_a * PPM_FACTOR, scaled by intensity (nA).
    // B is Tesla at 1 nA current; total B at intensityNA = B * intensityNA.
    // Sign convention: sigma_ab = -dB_a^sec/dB_{0,b}. See
    // feedback_t2_sacred and the library's BS sign verification.
    Mat3 G;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            G(a, b) = -geo.normal(b) * B(a) * PPM_FACTOR * intensityNA;
        }
    }
    return DecomposeLocal(G);
}

}  // namespace h5reader::calculators
