// PlaneFrameMath — pure-function plane geometry helpers extracted from
// MoleculeScene so they can be unit-tested without a VTK renderer or a
// loaded conformation.
//
// Two operations:
//
//   computePlaneFrame(positions) → optional<PlaneFrame>
//     Build an orthonormal {origin, x, y, z} frame from three atom
//     positions. Origin is the centroid; x is along (positions[1] -
//     positions[0]) normalised; z is the normalised cross product
//     (a→b) × (a→c); y is z × x. Returns nullopt for any degenerate
//     input (zero-length edge, collinear, coincident).
//
//   chooseContinuousNormal(previous, current) → Vec3
//     Frame-to-frame normal-sign continuity guard. If
//     dot(current, previous) < 0 the natural cross product flipped
//     between frames (geometry crossed a near-degenerate configuration
//     such as a ring flip or the third atom crossing the line through
//     the first two); return -current so the consumer's locked normal
//     stays on the same side it was on the previous frame. Otherwise
//     return current unchanged.
//
// Both are header-only; the only dependency is h5reader::model::Vec3
// (Eigen typedef from QtAtom.h). No VTK, no Qt, no I/O — by design.

#pragma once

#include "../model/QtAtom.h"

#include <array>
#include <cmath>
#include <optional>

namespace h5reader::math {

struct PlaneFrame {
    model::Vec3 origin = model::Vec3::Zero();
    model::Vec3 x = model::Vec3::Zero();
    model::Vec3 y = model::Vec3::Zero();
    model::Vec3 z = model::Vec3::Zero();
};

inline std::optional<PlaneFrame> computePlaneFrame(
    const std::array<model::Vec3, 3>& positions) {
    constexpr double kEdgeMin = 1e-6;

    const model::Vec3& a = positions[0];
    const model::Vec3& b = positions[1];
    const model::Vec3& c = positions[2];

    model::Vec3 x = b - a;
    if (x.norm() < kEdgeMin)
        return std::nullopt;
    x.normalize();

    model::Vec3 z = (b - a).cross(c - a);
    if (z.norm() < kEdgeMin)
        return std::nullopt;
    z.normalize();

    model::Vec3 y = z.cross(x);
    if (y.norm() < kEdgeMin)
        return std::nullopt;
    y.normalize();

    PlaneFrame frame;
    frame.origin = (a + b + c) / 3.0;
    frame.x = x;
    frame.y = y;
    frame.z = z;
    return frame;
}

inline model::Vec3 chooseContinuousNormal(const model::Vec3& previous,
                                          const model::Vec3& current) {
    return current.dot(previous) < 0.0 ? model::Vec3(-current) : current;
}

}  // namespace h5reader::math
