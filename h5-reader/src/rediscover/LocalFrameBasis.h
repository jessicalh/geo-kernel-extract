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
//   Backbone frames (N / Cα / C(=O) / HA — broad-backbone, ADDED 2026-06-01):
//   the conventions doc (substrate_conventions_2026-05-30.md, "Local frames
//   per atom class") specifies HA / Cα / C=O; N is DEFINED here analogously
//   from N's typed backbone neighbours (see BuildBackboneNFrame below). All
//   four are pure Vec3 math, anchored on TYPED backbone atoms the caller looks
//   up collision-safe (selectUnique / the QtResidue backbone cache) — never a
//   positional index (the IUPAC-revert trap). The existing HN / aromatic-H
//   builders are UNTOUCHED so the ring/mc oracle byte-parity holds.
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
    // ── Backbone frames (broad-backbone, ADDED 2026-06-01) ──
    BackboneN,           // amide N frame (z along N→CA; in-plane ref N→C_prev)
    BackboneN_NTerminus, // N-terminus: no C_prev, in-plane ref is N→C(own)
    BackboneCA,          // Cα frame (z = bisector of Cα→N and Cα→C; x along Cα→N)
    BackboneCarbonylC,   // backbone carbonyl C frame (z along C→O; in-plane ref C→CA)
    BackboneCarbonylO,   // backbone carbonyl O frame (z along O→C; in-plane ref C→CA)
    BackboneHA,          // Hα frame (z along Cα→HA; x along Cα→N)
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

// ── Backbone frames (broad-backbone, ADDED 2026-06-01) ─────────────────────
// All four take TYPED backbone-atom positions (looked up collision-safe by the
// caller: QtResidue's N/CA/C/O/HA cache or selectUnique on a typed locant —
// never a positional index). All right-handed (y = z × x), all return
// is_valid==false on degenerate geometry (coincident/collinear anchors) with
// NO NaN poisoning. Conventions per substrate_conventions_2026-05-30.md
// ("Local frames per atom class") for Cα / C=O / HA; N is DEFINED here.

// Backbone amide-N frame (DEFINED here — the conventions doc gives HA/Cα/CO,
// not N). Convention: the N atom's bonded backbone neighbours fix the frame.
//   z = unit(CA − N)                   — the N→Cα bond, the dominant local axis
//   x = in-plane component of the peptide reference ⟂ z, normalised:
//         interior residue: ref = (C_prev − N) (the preceding carbonyl C —
//                           this is the peptide-bond partner, so {N, CA, C_prev}
//                           spans the amide plane the N sits in);
//         N-terminus (no C_prev): ref = (C_own − N) (the residue's own carbonyl
//                           C; {N, CA, C_own} still spans a backbone plane).
//   y = z × x
// Chosen so z is the strongest covalent bond direction at N and x lies in the
// peptide plane N is part of — directly analogous to the HN frame (z = N→H,
// x in the amide plane) but with the N→Cα bond as z since N has no H in PRO.
// `c_prev_valid==false` selects the N-terminus fallback.
LocalFrame BuildBackboneNFrame(const Vec3& nPos, const Vec3& caPos,
                               const Vec3& cRefPos, bool c_prev_valid);

// Cα frame (conventions doc "Cα frame"):
//   z = unit bisector of (Cα→N) and (Cα→C), pointing away from the backbone
//   x = in-plane component of (Cα→N) ⟂ z, normalised
//   y = z × x
LocalFrame BuildBackboneCaFrame(const Vec3& caPos, const Vec3& nPos, const Vec3& cPos);

// Backbone carbonyl-C frame (conventions doc "C=O carbonyl frame", referenced
// from the carbonyl carbon):
//   z = unit(O − C)                    — the carbonyl bond direction
//   x = in-plane component of (CA − C) ⟂ z, normalised (the in-peptide-plane
//       reference; CA is the residue's own α-carbon, always present)
//   y = z × x
// NOTE: this z-axis IS the McConnell kernel's reference direction.
LocalFrame BuildBackboneCarbonylCFrame(const Vec3& cPos, const Vec3& oPos, const Vec3& caPos);

// Backbone carbonyl-O frame (the O sits on the carbonyl bond; same plane as the
// C frame, referenced from O so the O atom's own displacement vectors are
// expressed relative to it):
//   z = unit(C − O)                    — O→C, the carbonyl bond direction
//   x = in-plane component of (CA − C) ⟂ z, normalised
//   y = z × x
LocalFrame BuildBackboneCarbonylOFrame(const Vec3& oPos, const Vec3& cPos, const Vec3& caPos);

// Hα frame (conventions doc "HA / Cα chirality frame"):
//   z = unit(HA − Cα)                  — the Cα→HA chirality direction
//   x = in-plane component of (N − Cα) ⟂ z, normalised (Cα→N)
//   y = z × x
LocalFrame BuildBackboneHaFrame(const Vec3& haPos, const Vec3& caPos, const Vec3& nPos);

}  // namespace h5reader::rediscover
