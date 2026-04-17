// QtHaighMallionCalc — closed-form Haigh-Mallion kernel evaluator.
//
// Ports the 7-point Gauss quadrature (Stroud T2:5-1) on a fan
// triangulation of the ring, with adaptive subdivision up to level 2,
// from nmr-shielding/src/HaighMallionResult.cpp::SurfaceIntegral and
// SampleShieldingAt. Same sign convention: V = H * normal,
// G_ab = -n_b * V_a, scaled by ring-type intensity.
//
// Free functions — thread-safe by construction. See QtBiotSavartCalc.h
// for the parallel design notes.

#pragma once

#include "../model/QtRing.h"
#include "../model/Types.h"

#include <vector>

namespace h5reader::calculators {

// Evaluate the HM geometric shielding kernel G at a point, decomposed
// into a SphericalTensor. intensityNA is the ring's current intensity
// (literature or TOML-calibrated), in nanoamperes-equivalent scaling;
// the HM surface integral gives H in Å⁻¹ which becomes shielding via
// the intensity scale factor.
model::SphericalTensor EvaluateShielding(
    const model::Vec3&                 pointAng,
    const model::RingGeometry&         geo,
    const std::vector<model::Vec3>&    verticesAng,
    double                             intensityNA);

}  // namespace h5reader::calculators
