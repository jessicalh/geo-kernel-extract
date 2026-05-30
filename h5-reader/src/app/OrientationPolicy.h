// OrientationPolicy — orientation override on top of a CameraMode's
// natural per-frame fit.
//
// Most camera modes (Atom, Bond, Dihedral, Plane) have a "natural"
// orientation implied by their geometry. The OrientationPolicy parameter
// to CameraComposer::setMode is the optional override. Kind::Default
// uses the natural pairing defined by the mode itself; the explicit
// kinds force one specific orientation.
//
// Natural pairings (per spec/viewport_pipeline_2026-05-30.md §2.3.3,
// agent decision §4-d in the implementation prompt):
//   CameraMode::Free     -> Free
//   CameraMode::Atom     -> Free
//   CameraMode::Bond     -> PerpendicularToBond
//   CameraMode::Dihedral -> DownAxis(atoms[1], atoms[2])
//   CameraMode::Plane    -> PerpendicularToPlane
//   CameraMode::Subset   -> Free
//
// SingleConformation handling: regardless of policy, a one-frame
// conformation runs the composer once at setFrame(0) and the per-frame
// fit reduces to that single sample. No code branch on conformation
// kind — frameCount() == 1 is just a degenerate (and cheap) case of the
// general path.

#pragma once

#include <array>
#include <cstddef>

namespace h5reader::app {

struct OrientationPolicy {
    enum class Kind {
        Default,               // use the CameraMode's natural default
        Free,                  // V/D inherited from gesture state (no override)
        PerpendicularToBond,   // ViewUp parallel to bond axis
        DownAxis,              // sight axis = axisAtoms[0] -> axisAtoms[1]
        PerpendicularToPlane,  // sight axis = plane normal
    };

    Kind kind = Kind::Default;
    // For DownAxis: which two atoms define the sight axis. The composer
    // expands axisAtoms[0..1] against the active CameraMode's atom list.
    std::array<std::size_t, 2> axisAtoms{0, 0};
};

inline const char* NameFor(OrientationPolicy::Kind k) {
    switch (k) {
        case OrientationPolicy::Kind::Default:              return "default";
        case OrientationPolicy::Kind::Free:                 return "free";
        case OrientationPolicy::Kind::PerpendicularToBond:  return "perp_bond";
        case OrientationPolicy::Kind::DownAxis:             return "down_axis";
        case OrientationPolicy::Kind::PerpendicularToPlane: return "perp_plane";
    }
    return "?";
}

inline OrientationPolicy DefaultPolicy() {
    OrientationPolicy p;
    p.kind = OrientationPolicy::Kind::Default;
    return p;
}

inline OrientationPolicy FreePolicy() {
    OrientationPolicy p;
    p.kind = OrientationPolicy::Kind::Free;
    return p;
}

inline OrientationPolicy DownAxisPolicy(std::size_t a, std::size_t b) {
    OrientationPolicy p;
    p.kind = OrientationPolicy::Kind::DownAxis;
    p.axisAtoms = {a, b};
    return p;
}

inline OrientationPolicy PerpToBondPolicy() {
    OrientationPolicy p;
    p.kind = OrientationPolicy::Kind::PerpendicularToBond;
    return p;
}

inline OrientationPolicy PerpToPlanePolicy() {
    OrientationPolicy p;
    p.kind = OrientationPolicy::Kind::PerpendicularToPlane;
    return p;
}

}  // namespace h5reader::app
