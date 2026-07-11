#include "RingNeighbourGeometry.h"

#include "CalculatorConfig.h"

#include <cmath>
#include <limits>

namespace nmr {

RingNeighbourhood& EnsureRingNeighbourGeometry(
        ConformationAtom& atom,
        const RingGeometry& geom,
        const Ring& ring,
        size_t ring_index,
        const Vec3& atom_pos) {
    RingNeighbourhood* row = nullptr;
    for (auto& existing : atom.ring_neighbours) {
        if (existing.ring_index == ring_index) {
            row = &existing;
            break;
        }
    }
    if (row == nullptr) {
        atom.ring_neighbours.emplace_back();
        row = &atom.ring_neighbours.back();
    }

    const double norm_guard =
        CalculatorConfig::Get("near_zero_vector_norm_threshold");
    const double nan = std::numeric_limits<double>::quiet_NaN();

    const Vec3 d = atom_pos - geom.center;
    const double distance = d.norm();

    row->ring_index = ring_index;
    row->ring_type = ring.type_index;
    row->distance_to_center = distance;
    row->direction_to_center =
        (std::isfinite(distance) && distance > 0.0)
            ? Vec3(d / distance) : Vec3::Zero();

    // Shared cylindrical geometry follows the stored ring normal.  Valid
    // RingGeometry normals are unit vectors; keeping these definitions in
    // one place also makes a malformed/degenerate geometry deterministic.
    row->z = d.dot(geom.normal);
    const Vec3 d_plane = d - row->z * geom.normal;
    row->rho = d_plane.norm();
    row->theta = std::atan2(row->rho, std::abs(row->z));

    // Deterministic in-plane gauge: vertex zero is fixed by RingTopology's
    // canonical cyclic walk.  Do not choose a fallback vertex when the gauge
    // is degenerate; a NaN azimuth exposes the invalid frame.
    const double normal_norm = geom.normal.norm();
    if (!std::isfinite(normal_norm) || normal_norm <= norm_guard ||
        geom.vertices.empty()) {
        row->cos_phi = nan;
        row->sin_phi = nan;
        return *row;
    }

    const Vec3 z_axis = geom.normal / normal_norm;
    const Vec3 reference = geom.vertices[0] - geom.center;
    const Vec3 reference_plane =
        reference - reference.dot(z_axis) * z_axis;
    const double reference_norm = reference_plane.norm();
    if (!std::isfinite(reference_norm) || reference_norm <= norm_guard) {
        row->cos_phi = nan;
        row->sin_phi = nan;
        return *row;
    }

    const Vec3 azimuth_plane = d - d.dot(z_axis) * z_axis;
    const double azimuth_radius = azimuth_plane.norm();
    if (!std::isfinite(azimuth_radius)) {
        row->cos_phi = nan;
        row->sin_phi = nan;
    } else if (azimuth_radius <= norm_guard) {
        // Azimuth is conventionally fixed on the ring axis.
        row->cos_phi = 1.0;
        row->sin_phi = 0.0;
    } else {
        const Vec3 d_hat = azimuth_plane / azimuth_radius;
        const Vec3 x_axis = reference_plane / reference_norm;
        row->cos_phi = d_hat.dot(x_axis);
        row->sin_phi = d_hat.cross(x_axis).dot(z_axis);
    }

    return *row;
}

}  // namespace nmr
