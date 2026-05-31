#include "LocalFrameBasis.h"

#include <cmath>

namespace h5reader::rediscover {

namespace {

// Gram-Schmidt: in-plane component of `ref` perpendicular to unit `z`,
// normalised. Returns false (and leaves `xOut` untouched) if the residual is
// degenerate (ref ∥ z).
bool inPlaneAxis(const Vec3& ref, const Vec3& z, Vec3& xOut) {
    const Vec3 perp = ref - z * ref.dot(z);
    const double n = perp.norm();
    if (!(n > 1e-9)) return false;
    xOut = perp / n;
    return true;
}

}  // namespace

LocalFrame BuildHNFrame(const Vec3& nPos, const Vec3& hPos, const Vec3& caPos,
                        const Vec3& cPrevPos, bool c_prev_valid) {
    LocalFrame f;

    // z = unit(H − N): the N–H bond direction.
    const Vec3 nh = hPos - nPos;
    const double nhNorm = nh.norm();
    if (!(nhNorm > 1e-9)) return f;  // is_valid stays false
    f.z = nh / nhNorm;

    // x: in-plane component of the amide-plane reference, ⟂ z.
    const Vec3 ref = c_prev_valid ? (cPrevPos - nPos) : (caPos - nPos);
    if (!inPlaneAxis(ref, f.z, f.x)) return f;

    f.y = f.z.cross(f.x).normalized();
    f.x = f.y.cross(f.z).normalized();  // re-orthogonalise (right-handed)
    f.is_valid = true;
    f.variant = c_prev_valid ? FrameVariant::HN_Standard : FrameVariant::HN_NTerminus;
    return f;
}

LocalFrame BuildAromaticHFrame(const Vec3& ringCenter, const Vec3& ringNormal,
                               const Vec3& anchorPos) {
    LocalFrame f;

    const double nNorm = ringNormal.norm();
    if (!(nNorm > 1e-9)) return f;
    f.z = ringNormal / nNorm;

    // x: centroid → anchor, projected into the ring plane (⟂ z).
    const Vec3 ref = anchorPos - ringCenter;
    if (!inPlaneAxis(ref, f.z, f.x)) return f;

    f.y = f.z.cross(f.x).normalized();
    f.x = f.y.cross(f.z).normalized();
    f.is_valid = true;
    f.variant = FrameVariant::AromaticHRing;
    return f;
}

}  // namespace h5reader::rediscover
