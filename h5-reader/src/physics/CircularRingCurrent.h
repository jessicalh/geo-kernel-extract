#pragma once

#include "../model/Types.h"

#include <optional>
#include <vector>

namespace h5reader::physics {

struct CircularRingParameters {
    double radiusA = 0.0;
    double lobeOffsetA = 0.0;
    double currentNanoamperePerTesla = 0.0;
};

struct CircularRingPlane {
    model::Vec3 center = model::Vec3::Zero();
    model::Vec3 normal = model::Vec3::Zero();
    double planeRmsA = 0.0;
};

// Least-squares plane used by the frozen Candidate-A calculation. This is
// intentionally separate from the Reader's winding-stable display geometry.
std::optional<CircularRingPlane> FitCircularRingPlane(const std::vector<model::Vec3>& vertices);

// Frozen circular-loop model selected by the small-molecule ring-current
// comparison. Unsupported ring types return no parameters.
std::optional<CircularRingParameters> CandidateACircularParameters(model::RingTypeIndex type, int protonationVariant = -1);

// Exact field of two coaxial circular loops at +/- lobeOffsetA, each carrying
// half the declared current. Positions are Angstrom, current is nA/T, and the
// returned field is Tesla per tesla of applied field.
std::optional<model::Vec3> EvaluateCircularDoubleLoopField(const model::Vec3& pointA,
                                                           const CircularRingPlane& plane,
                                                           const CircularRingParameters& parameters);

// G_ab = -B_a n_b * 1e6, decomposed in the reader's library-compatible basis.
std::optional<model::SphericalTensor>
EvaluateCircularShielding(const model::Vec3& pointA, const CircularRingPlane& plane, const CircularRingParameters& parameters);

}  // namespace h5reader::physics
