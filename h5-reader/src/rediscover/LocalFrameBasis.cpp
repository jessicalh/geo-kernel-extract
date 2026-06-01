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

namespace {

// Shared finish: given a unit z and an in-plane reference, Gram-Schmidt x,
// right-handed y, re-orthogonalise x. Sets is_valid + the variant on success.
// Leaves is_valid==false (no NaN) when the reference is parallel to z.
LocalFrame finishFrame(const Vec3& zUnit, const Vec3& ref, FrameVariant variant) {
    LocalFrame f;
    f.z = zUnit;
    if (!inPlaneAxis(ref, f.z, f.x)) return f;  // is_valid stays false
    f.y = f.z.cross(f.x).normalized();
    f.x = f.y.cross(f.z).normalized();
    f.is_valid = true;
    f.variant = variant;
    return f;
}

}  // namespace

LocalFrame BuildBackboneNFrame(const Vec3& nPos, const Vec3& caPos,
                               const Vec3& cRefPos, bool c_prev_valid) {
    // z = unit(CA − N): the N→Cα bond (the dominant covalent axis at N).
    const Vec3 nca = caPos - nPos;
    const double n = nca.norm();
    if (!(n > 1e-9)) return LocalFrame{};  // coincident N/CA — invalid, no NaN
    // x: in-plane component of the peptide reference (N→C_prev interior, else
    // N→C_own) ⟂ z. {N, CA, ref} spans the peptide plane N sits in.
    const Vec3 ref = cRefPos - nPos;
    return finishFrame(nca / n, ref,
                       c_prev_valid ? FrameVariant::BackboneN
                                    : FrameVariant::BackboneN_NTerminus);
}

LocalFrame BuildBackboneCaFrame(const Vec3& caPos, const Vec3& nPos, const Vec3& cPos) {
    // z = unit bisector of (Cα→N) and (Cα→C): the out-of-backbone direction.
    const Vec3 toN = nPos - caPos;
    const Vec3 toC = cPos - caPos;
    const double nN = toN.norm();
    const double nC = toC.norm();
    if (!(nN > 1e-9) || !(nC > 1e-9)) return LocalFrame{};  // coincident — invalid
    const Vec3 bisector = (toN / nN) + (toC / nC);
    const double nb = bisector.norm();
    if (!(nb > 1e-9)) return LocalFrame{};  // N–Cα–C collinear (anti-parallel)
    // x: in-plane component of (Cα→N) ⟂ z.
    return finishFrame(bisector / nb, toN, FrameVariant::BackboneCA);
}

LocalFrame BuildBackboneCarbonylCFrame(const Vec3& cPos, const Vec3& oPos, const Vec3& caPos) {
    // z = unit(O − C): the carbonyl bond direction (the McConnell reference).
    const Vec3 co = oPos - cPos;
    const double n = co.norm();
    if (!(n > 1e-9)) return LocalFrame{};  // coincident C/O — invalid
    // x: in-plane component of (CA − C) ⟂ z (in the peptide plane).
    const Vec3 ref = caPos - cPos;
    return finishFrame(co / n, ref, FrameVariant::BackboneCarbonylC);
}

LocalFrame BuildBackboneCarbonylOFrame(const Vec3& oPos, const Vec3& cPos, const Vec3& caPos) {
    // z = unit(C − O): O→C along the carbonyl bond.
    const Vec3 oc = cPos - oPos;
    const double n = oc.norm();
    if (!(n > 1e-9)) return LocalFrame{};  // coincident O/C — invalid
    // x: in-plane component of (CA − C) ⟂ z (the same peptide plane).
    const Vec3 ref = caPos - cPos;
    return finishFrame(oc / n, ref, FrameVariant::BackboneCarbonylO);
}

LocalFrame BuildBackboneHaFrame(const Vec3& haPos, const Vec3& caPos, const Vec3& nPos) {
    // z = unit(HA − Cα): the Cα→HA chirality direction.
    const Vec3 cah = haPos - caPos;
    const double n = cah.norm();
    if (!(n > 1e-9)) return LocalFrame{};  // coincident Cα/HA — invalid
    // x: in-plane component of (N − Cα) ⟂ z (Cα→N).
    const Vec3 ref = nPos - caPos;
    return finishFrame(cah / n, ref, FrameVariant::BackboneHA);
}

}  // namespace h5reader::rediscover
