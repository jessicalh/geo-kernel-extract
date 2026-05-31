// LocalFrameBasis — per-atom-class right-handed local frames, so the source
// displacement vectors AND the target tensor land in one rotation-stable basis
// (what an equivariant fitter needs; a flat r/cosθ table throws the vector
// information away).
//
// Conventions are pinned by spec/substrate_conventions_2026-05-30.md
// ("Local frames per atom class"). All frames are right-handed: y = z × x.
// `is_valid == false` means the frame could not be constructed for this atom
// at this frame (an edge case below); callers must check before using the
// frame to express a vector.
//
//   HN frame (backbone amide hydrogen):
//     z = unit(H − N)                              (the N–H bond direction)
//     x = in-plane component of (C_prev − N) ⟂ z, normalised
//         (N-terminus fallback: use (Cα − N))
//     y = z × x
//
//   Aromatic-H frame (ring-facing aromatic hydrogen):
//     z = unit ring normal (fixed traversal order, "same direction every
//         frame" — fit from the ring vertex positions)
//     x = unit(anchor − centroid), anchor = chemistry-typed atom (CG / CD2
//         per the conventions doc), projected ⟂ z and normalised
//     y = z × x
//
// These are plain math over Vec3 (no QObject, no model coupling); the
// extraction code supplies the atom positions it already has from the
// Conformation / H5 and the ring geometry from ConformationGeometry.

#pragma once

#include "../model/Types.h"  // Vec3, Mat3

namespace h5reader::rediscover {

using model::Mat3;
using model::Vec3;

// Which atom-class construction produced a frame (recorded per record so the
// downstream fitter knows which convention is in force). Mirrors the
// FrameVariant idea in the conventions doc, trimmed to what the two
// extractions actually emit.
enum class FrameVariant : int {
    Invalid = 0,
    HN_Standard,     // interior residue, C_prev available
    HN_NTerminus,    // no C_prev: in-plane reference is (Cα − N)
    AromaticHRing,   // ring-normal frame, anchored on the typed ring atom
};

struct LocalFrame {
    Vec3 z = Vec3::UnitZ();
    Vec3 x = Vec3::UnitX();
    Vec3 y = Vec3::UnitY();
    bool is_valid = false;
    FrameVariant variant = FrameVariant::Invalid;

    // Express a lab-frame displacement vector in this frame's {x,y,z} basis.
    // (v·x, v·y, v·z). Identity when the frame is the lab axes.
    Vec3 ToLocal(const Vec3& v_lab) const { return Vec3(v_lab.dot(x), v_lab.dot(y), v_lab.dot(z)); }

    // Rotate a lab-frame rank-2 tensor T into the local frame: Rᵀ T R, with
    // R = [x y z] (columns). Used to express the DFT target tensor in the
    // same rotation-stable frame as the source vectors.
    Mat3 TensorToLocal(const Mat3& t_lab) const {
        Mat3 r;
        r.col(0) = x;
        r.col(1) = y;
        r.col(2) = z;
        return r.transpose() * t_lab * r;
    }
};

// HN amide-plane frame from the backbone atom positions of this residue.
// c_prev_valid==false selects the N-terminus fallback (in-plane reference
// from Cα instead of C_prev). Returns is_valid==false on degenerate geometry
// (coincident N/H, or an in-plane reference parallel to z).
LocalFrame BuildHNFrame(const Vec3& nPos, const Vec3& hPos, const Vec3& caPos,
                        const Vec3& cPrevPos, bool c_prev_valid);

// Aromatic ring-normal frame for an H facing ring `ringCenter` with unit
// normal `ringNormal` and a chemistry-typed anchor atom position. The anchor
// fixes the x-axis azimuth so the frame is stable across tautomers/frames.
// Returns is_valid==false if the normal is degenerate or anchor ∥ normal.
LocalFrame BuildAromaticHFrame(const Vec3& ringCenter, const Vec3& ringNormal,
                               const Vec3& anchorPos);

}  // namespace h5reader::rediscover
