// ConformationGeometry — shared geometry helpers over a Conformation.
//
// Ring vertices + the least-squares plane fit used to be private to
// QtFrame (trajectory-only). They are pure functions of positions +
// topology, so they belong over the Conformation base: a SingleConformation
// (one pose, no H5) and a TrajectoryConformation share the same ring
// geometry the same way. The trajectory.h5 format does not store per-frame
// ring geometry — it is fit on demand from the ring's vertex positions,
// matching nmr::Ring::ComputeGeometry's least-squares plane fit.
//
// `frame` indexes the Conformation's frame axis (0 for a single pose); the
// positions come from Conformation::atomPosition, so both run shapes work.

#pragma once

#include "QtRing.h"  // RingGeometry
#include "Types.h"

#include <cstddef>
#include <vector>

namespace h5reader::model {

class Conformation;

// Vertices of aromatic ring `ringIdx` at `frame`, in canonical-walk order
// (Conformation::atomPosition over the ring's atomIndices). Empty if the
// ring index is out of range or the protein spine is absent.
std::vector<Vec3> RingVertices(const Conformation& conf, std::size_t ringIdx, std::size_t frame);

// Least-squares plane fit: center = centroid, normal = smallest singular
// vector of the centered vertex matrix (SVD), radius = mean spoke length.
// Returns a zero geometry (radius 0) for < 3 vertices.
RingGeometry FitRingGeometry(const std::vector<Vec3>& verts);

// Convenience: RingVertices + FitRingGeometry at (ringIdx, frame).
RingGeometry RingGeometryAt(const Conformation& conf, std::size_t ringIdx, std::size_t frame);

// ---------------------------------------------------------------------------
// Geometry of an ordered atom selection — the killer app's measurements.
//
// Distance/Angle/Dihedral are pure functions of positions; Measure() reads
// them through Conformation::atomPosition, so a single pose and a trajectory
// frame are measured identically. AngleDegrees takes the MIDDLE atom as the
// vertex. DihedralDegrees uses the signed Blondel-Karplus atan2 convention
// (Blondel & Karplus 1996, J. Comput. Chem. 17(9):1132), range (-180, 180],
// matching the library's omega_actual / chi / pucker sign so a
// measured-from-positions dihedral can be cross-checked against the extracted
// field. Degenerate input (coincident/collinear atoms -> undefined direction)
// yields NaN; Measure() reports that as valid == false.
// ---------------------------------------------------------------------------

// GeometryKind is owned by AtomSelection.h (it answers "what does this
// cardinality measure"); forward-declared here as an opaque scoped enum
// (underlying type int, matching AtomSelection's definition) so this geometry
// header stays free of the Qt model machinery. The .cpp includes
// AtomSelection.h for the enumerator values.
enum class GeometryKind;

double Distance(const Vec3& a, const Vec3& b);                                      // Å
double AngleDegrees(const Vec3& a, const Vec3& b, const Vec3& c);                   // vertex = b; [0, 180]
double DihedralDegrees(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d); // (-180, 180]

// The measured value of a k-tuple selection: its kind (from cardinality), the
// value (Å for Distance; degrees for Angle/Dihedral), and whether it is
// well-defined. kind value-initialises to GeometryKind::None (ordinal 0).
struct GeometryMeasurement {
    GeometryKind kind  = {};
    double       value = 0.0;
    bool         valid = false;
};

// Measure the ordered selection at `frame`: 2 atoms -> Distance, 3 -> Angle
// (vertex = middle), 4 -> Dihedral; any other count -> {None, 0, invalid}.
// Bounds-checks the frame and every atom index against the protein.
GeometryMeasurement Measure(const Conformation& conf, std::size_t frame,
                            const std::vector<std::size_t>& atoms);

}  // namespace h5reader::model
