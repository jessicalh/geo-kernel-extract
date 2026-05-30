// CameraMode — typed sum for the active camera-lock target.
//
// Per spec/viewport_pipeline_2026-05-30.md §2.3.1: each variant carries
// the atom indices that define the lock and a built-in orientation
// policy default. The composer reads positions for those atoms each
// frame and writes the absolute camera state.
//
// Kind taxonomy:
//   Free      — no lock; the camera state is owned by accumulated user
//               gestures.
//   Atom      — 1 atom; focal point follows that atom each frame.
//   Bond      — 2 atoms; focal point at midpoint, sight axis = (b - a).
//   Dihedral  — 4 atoms; focal point at midpoint(b, c), sight down (c - b)
//               (Newman projection).
//   Plane     — 3 atoms; focal point at centroid, sight perpendicular to
//               the plane normal (sign-continuity guarded).
//   Subset    — N >= 3 atoms; Kabsch fit of the subset against its
//               reference positions; the camera follows the stabilised
//               local frame. The bridge to the position-side
//               TransformedConformation's FitSubset mode.
//
// Lock release semantics — agent decision per the implementation prompt:
//   * Plane locks release on selection change (matches existing behaviour;
//     the plane lock toolbar action already wires this).
//   * Atom / Bond / Dihedral / Subset locks STAY on selection change. The
//     user committed to a focus that survives small selection edits.
//
// User-gesture deltas reset on every mode change (cleanest mental model —
// each lock acquisition is a fresh start; the previous mode's azimuth
// does not bleed into the new mode's natural orientation).

#pragma once

#include "../diagnostics/ThreadGuard.h"

#include <cstddef>
#include <vector>

namespace h5reader::app {

struct CameraMode {
    enum class Kind {
        Free,
        Atom,
        Bond,
        Dihedral,
        Plane,
        Subset,
    };

    Kind kind = Kind::Free;
    std::vector<std::size_t> atoms;  // empty iff Kind::Free; size implies kind

    bool operator==(const CameraMode& other) const {
        return kind == other.kind && atoms == other.atoms;
    }
    bool operator!=(const CameraMode& other) const { return !(*this == other); }
};

// Expected atom counts per kind. Returns 0 for Free, SIZE_MAX for Subset
// (variable). Callers building a mode by hand should check against these.
inline std::size_t ExpectedAtomCount(CameraMode::Kind k) {
    switch (k) {
        case CameraMode::Kind::Free:     return 0;
        case CameraMode::Kind::Atom:     return 1;
        case CameraMode::Kind::Bond:     return 2;
        case CameraMode::Kind::Dihedral: return 4;
        case CameraMode::Kind::Plane:    return 3;
        case CameraMode::Kind::Subset:   return static_cast<std::size_t>(-1);
    }
    return 0;
}

inline const char* NameFor(CameraMode::Kind k) {
    switch (k) {
        case CameraMode::Kind::Free:     return "free";
        case CameraMode::Kind::Atom:     return "atom";
        case CameraMode::Kind::Bond:     return "bond";
        case CameraMode::Kind::Dihedral: return "dihedral";
        case CameraMode::Kind::Plane:    return "plane";
        case CameraMode::Kind::Subset:   return "subset";
    }
    return "?";
}

// Typed constructors. Each rejects mismatched cardinality by returning
// a Free mode (defensive default — the REST handler bounds-checks atoms
// before constructing, so this only fires on programmer error).
inline CameraMode FreeMode() { return {}; }

inline CameraMode AtomMode(std::size_t a) {
    CameraMode m;
    m.kind  = CameraMode::Kind::Atom;
    m.atoms = {a};
    return m;
}

inline CameraMode BondMode(std::size_t a, std::size_t b) {
    CameraMode m;
    m.kind  = CameraMode::Kind::Bond;
    m.atoms = {a, b};
    return m;
}

inline CameraMode DihedralMode(std::size_t a, std::size_t b,
                                std::size_t c, std::size_t d) {
    CameraMode m;
    m.kind  = CameraMode::Kind::Dihedral;
    m.atoms = {a, b, c, d};
    return m;
}

inline CameraMode PlaneMode(std::size_t a, std::size_t b, std::size_t c) {
    CameraMode m;
    m.kind  = CameraMode::Kind::Plane;
    m.atoms = {a, b, c};
    return m;
}

inline CameraMode SubsetMode(std::vector<std::size_t> atoms) {
    CameraMode m;
    m.kind  = CameraMode::Kind::Subset;
    m.atoms = std::move(atoms);
    return m;
}

}  // namespace h5reader::app
