// MolecularFrame -- per-frame molecular (chemical-group) reference frame and
// the lab->molecular tensor projection, for the viewer's standalone "look at a
// DFT that hasn't had a rediscover pass yet" path.
//
// This is the GEOMETRY CORE: pure Eigen + the shared model Vec3/Mat3 types, no
// Qt, no VTK, no QtProtein/Conformation coupling -- so it is headless-testable
// from tests/scene_math_tests and could later be shared with the extractor's
// AnalysisAtom in place of its inline copy. It takes a resolved MolFrameSpec
// (anchor atom indices already chosen) plus position / ring-geometry callbacks,
// and returns the per-frame orthonormal axes.
//
// The frame-BUILDING math (frameFromXAndPlane / XAndZ / XAndPlaneLocked) is
// grabbed faithfully from AnalysisAtom.cpp (the rediscover engine on the
// nmr-shielding-analysisatom worktree, ~lines 751-816 / 5166-5223) so the axes
// match the extractor. What is DELIBERATELY NOT ported is its way of CHOOSING
// the anchors: AnalysisAtom dispatches on IUPAC name strings ("C","OD1",...);
// the viewer's QtAtom is already typed (BackboneRole, PlanarGroupKind,
// parentAtomIndex), so the selector that fills a MolFrameSpec lives in a thin
// typed adapter (see SelectMolecularFrameSpec, alongside the viewer model) and
// dispatches off those enums -- saner, protein-agnostic, no string fragility,
// same frames. The extractor's selection way is not the viewer's; the geometry
// is shared.
//
// MolComp ordering is the engine's {xx, yy, xy, xz, yz, zz} (AnalysisAtom.cpp:
// 1421) -- it MUST match or every molecular-component relationship silently
// corrupts.

#pragma once

#include "Types.h"  // Vec3, Mat3 (Eigen)

#include <Eigen/Dense>

#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <utility>

namespace h5reader::model {

enum class MolecularFrameKind : std::uint8_t {
    None = 0,
    BackboneCarbonyl,
    BackboneAmideN,
    AromaticRingLocal,
    SaturatedRingLocal,
    MetSd,
    SidechainCarboxylate,
    SidechainGuanidinium,
    SidechainCarboxamide,
};

inline const char* MolecularFrameKindName(MolecularFrameKind k) {
    switch (k) {
    case MolecularFrameKind::None: return "none";
    case MolecularFrameKind::BackboneCarbonyl: return "backbone_carbonyl";
    case MolecularFrameKind::BackboneAmideN: return "backbone_amide_n";
    case MolecularFrameKind::AromaticRingLocal: return "aromatic_ring_local";
    case MolecularFrameKind::SaturatedRingLocal: return "saturated_ring_local";
    case MolecularFrameKind::MetSd: return "met_sd";
    case MolecularFrameKind::SidechainCarboxylate: return "sidechain_carboxylate";
    case MolecularFrameKind::SidechainGuanidinium: return "sidechain_guanidinium";
    case MolecularFrameKind::SidechainCarboxamide: return "sidechain_carboxamide";
    }
    return "none";
}

// Molecular-frame symmetric-tensor component order. MUST equal the engine's
// enum (AnalysisAtom.cpp:1421) -- a reorder silently corrupts every component
// relationship the viewer would show against the extractor's reveals.
enum class MolComp : std::uint8_t { xx = 0, yy, xy, xz, yz, zz };

// A resolved frame request: which kind, and the anchor ATOM INDICES the typed
// selector chose. Unused anchors stay -1. Geometry only -- no chemistry here.
struct MolFrameSpec {
    MolecularFrameKind kind = MolecularFrameKind::None;
    std::int32_t origin = -1;        // frame origin atom (the framed centre)
    std::int32_t xAnchor = -1;       // sets the x director
    std::int32_t planeAnchor = -1;   // defines the x-y plane (or 2nd ref)
    std::int32_t secondAnchor = -1;  // carboxylate: the 2nd carboxyl O
    std::int32_t ring = -1;          // aromatic ring index (ring-local kind)
    std::int32_t heavy = -1;         // aromatic: heavy atom whose radial sets x
};

// Cross-frame sign/handedness continuity for the "locked" kinds, so a glyph's
// axes do not flip frame-to-frame on a tiny twitch. Thread one instance per
// (atom) across the trajectory; the locked builders update it in place.
struct MolFrameContinuity {
    std::optional<Vec3> x;
    std::optional<Vec3> z;
};

// ---- geometry helpers (grabbed faithfully from AnalysisAtom.cpp:751-816) ----

inline std::optional<Vec3> normalizeFrameVec(const Vec3& v, double eps = 1e-12) {
    const double norm = v.norm();
    if (!std::isfinite(norm) || !(norm > eps)) return std::nullopt;
    return v / norm;
}

inline std::optional<Mat3> frameFromColumns(const Vec3& x, const Vec3& y, const Vec3& z) {
    Mat3 frame;
    frame.col(0) = x;
    frame.col(1) = y;
    frame.col(2) = z;
    if (!frame.allFinite()) return std::nullopt;
    return frame;
}

// x director from xVec; y in the (xVec, planeVec) plane; z = x cross y.
inline std::optional<Mat3> frameFromXAndPlane(const Vec3& xVec, const Vec3& planeVec) {
    const auto xOpt = normalizeFrameVec(xVec);
    if (!xOpt) return std::nullopt;
    const Vec3& x = *xOpt;
    const Vec3 plane = planeVec - planeVec.dot(x) * x;
    const auto y0Opt = normalizeFrameVec(plane);
    if (!y0Opt) return std::nullopt;
    const auto zOpt = normalizeFrameVec(x.cross(*y0Opt));
    if (!zOpt) return std::nullopt;
    const auto yOpt = normalizeFrameVec(zOpt->cross(x));
    if (!yOpt) return std::nullopt;
    return frameFromColumns(x, *yOpt, *zOpt);
}

// z director from zVec (e.g. a ring normal); x from xVec projected into the
// plane; y = z cross x. (Aromatic-ring-local frame.)
inline std::optional<Mat3> frameFromXAndZ(const Vec3& xVec, const Vec3& zVec) {
    const auto zOpt = normalizeFrameVec(zVec);
    if (!zOpt) return std::nullopt;
    const Vec3& z = *zOpt;
    const Vec3 xProj = xVec - xVec.dot(z) * z;
    const auto xOpt = normalizeFrameVec(xProj);
    if (!xOpt) return std::nullopt;
    const auto yOpt = normalizeFrameVec(z.cross(*xOpt));
    if (!yOpt) return std::nullopt;
    return frameFromColumns(*xOpt, *yOpt, z);
}

// Like frameFromXAndPlane, but flips x / z to stay continuous with the previous
// frame's directors (prevX / prevZ), then writes them back. Used by the kinds
// whose anchors do not pin an absolute handedness frame-to-frame.
inline std::optional<Mat3> frameFromXAndPlaneLocked(const Vec3& xVec,
                                                    const Vec3& planeVec,
                                                    std::optional<Vec3>& prevX,
                                                    std::optional<Vec3>& prevZ) {
    auto xOpt = normalizeFrameVec(xVec);
    if (!xOpt) return std::nullopt;
    Vec3 x = *xOpt;
    const Vec3 plane = planeVec - planeVec.dot(x) * x;
    const auto y0Opt = normalizeFrameVec(plane);
    if (!y0Opt) return std::nullopt;
    auto zOpt = normalizeFrameVec(x.cross(*y0Opt));
    if (!zOpt) return std::nullopt;
    Vec3 z = *zOpt;

    if (prevX && prevX->dot(x) < 0.0) x *= -1.0;
    if (prevZ && prevZ->dot(z) < 0.0) z *= -1.0;

    const auto yOpt = normalizeFrameVec(z.cross(x));
    if (!yOpt) return std::nullopt;
    const auto frame = frameFromColumns(x, *yOpt, z);
    if (!frame) return std::nullopt;
    prevX = x;
    prevZ = z;
    return frame;
}

// ---- per-frame axes from a resolved spec (pure; positions via callbacks) ----

// posOf(atomIndex)            -> lab position at the frame of interest, or null.
// ringCenterNormalOf(ringIdx) -> (center, normal) at that frame, or null.
// cont                        -> cross-frame continuity for the locked kinds.
using PositionLookup = std::function<std::optional<Vec3>(std::int32_t)>;
using RingFrameLookup = std::function<std::optional<std::pair<Vec3, Vec3>>(std::int32_t)>;

inline std::optional<Mat3> BuildMolecularFrameAxes(const MolFrameSpec& spec,
                                                   const PositionLookup& posOf,
                                                   const RingFrameLookup& ringCenterNormalOf,
                                                   MolFrameContinuity& cont) {
    switch (spec.kind) {
    case MolecularFrameKind::None:
        return std::nullopt;

    case MolecularFrameKind::BackboneCarbonyl:
    case MolecularFrameKind::BackboneAmideN: {
        const auto origin = posOf(spec.origin);
        const auto xAnchor = posOf(spec.xAnchor);
        const auto planeAnchor = posOf(spec.planeAnchor);
        if (origin && xAnchor && planeAnchor)
            return frameFromXAndPlane(*xAnchor - *origin, *planeAnchor - *origin);
        return std::nullopt;
    }

    case MolecularFrameKind::MetSd:
    case MolecularFrameKind::SidechainGuanidinium:
    case MolecularFrameKind::SidechainCarboxamide: {
        const auto origin = posOf(spec.origin);
        const auto xAnchor = posOf(spec.xAnchor);
        const auto planeAnchor = posOf(spec.planeAnchor);
        if (origin && xAnchor && planeAnchor)
            return frameFromXAndPlaneLocked(*xAnchor - *origin, *planeAnchor - *origin,
                                            cont.x, cont.z);
        return std::nullopt;
    }

    case MolecularFrameKind::AromaticRingLocal:
    case MolecularFrameKind::SaturatedRingLocal: {
        if (spec.ring < 0 || spec.heavy < 0) return std::nullopt;
        const auto heavy = posOf(spec.heavy);
        const auto ring = ringCenterNormalOf(spec.ring);
        if (heavy && ring)
            return frameFromXAndZ(*heavy - ring->first, ring->second);
        return std::nullopt;
    }

    case MolecularFrameKind::SidechainCarboxylate: {
        // Bisector frame: x along (O1 + O2 bisector) from C, plane toward O1.
        const auto c = posOf(spec.origin);
        const auto o1 = posOf(spec.xAnchor);
        const auto o2 = posOf(spec.secondAnchor);
        if (c && o1 && o2) {
            const auto d1 = normalizeFrameVec(*o1 - *c);
            const auto d2 = normalizeFrameVec(*o2 - *c);
            if (d1 && d2)
                return frameFromXAndPlaneLocked(*d1 + *d2, *o1 - *c, cont.x, cont.z);
        }
        return std::nullopt;
    }
    }
    return std::nullopt;
}

// ---- lab -> molecular projection + component readout ----

// Rotate a lab-frame (rank-2) tensor into the molecular frame: axes^T * T * axes
// (axes columns are the molecular x, y, z directors in lab coordinates).
inline Mat3 ProjectToMolecularFrame(const Mat3& molecularAxes, const Mat3& labTensor) {
    return molecularAxes.transpose() * labTensor * molecularAxes;
}

// Six unique symmetric components in MolComp order {xx, yy, xy, xz, yz, zz}.
// The tensor is symmetrised first (the antisymmetric T1 part is not a CSA
// component); this matches the engine's symmetricComponents.
inline std::array<double, 6> symmetricComponents(const Mat3& m) {
    return {m(0, 0),
            m(1, 1),
            0.5 * (m(0, 1) + m(1, 0)),
            0.5 * (m(0, 2) + m(2, 0)),
            0.5 * (m(1, 2) + m(2, 1)),
            m(2, 2)};
}

inline double component(const Mat3& m, MolComp c) {
    return symmetricComponents(m)[static_cast<std::size_t>(c)];
}

}  // namespace h5reader::model
