// QtBiotSavartCalc — closed-form Biot-Savart kernel evaluator.
//
// Ports the Johnson-Bovey two-loop wire integral from
// nmr-shielding/src/BiotSavartResult.cpp::JohnsonBoveyField +
// WireSegmentField + SampleShieldingAt. Same sign convention:
// G_ab = -n_b * B_a * PPM_FACTOR.
//
// Free functions on const inputs — no shared state, no Qt, no VTK.
// Runnable on any thread. Today the overlay calls these on the GUI
// thread inside setFrame(); if profiling later shows the cost hurts
// playback at high ring count, a worker thread can call them with
// zero changes here.

#pragma once

#include "../model/QtRing.h"   // RingGeometry
#include "../model/Types.h"    // Vec3, SphericalTensor

#include <vector>

namespace h5reader::calculators {

// Evaluate the B-field in Tesla at a point, from the Johnson-Bovey
// two-loop model for one ring. Vertices are the ring atom positions
// in Angstroms (the ring polygon); lobeOffsetA is the JB offset of
// each loop from the ring plane along the normal (ring-type-dependent);
// currentNA is the ring current in nanoamperes (ring-type intensity).
//
// Returns Vec3::Zero() if the point is inside the ring (distance <
// ring radius) or too close (distance < SINGULARITY_GUARD_DISTANCE)
// or too far (distance > RING_CURRENT_CUTOFF).
model::Vec3 EvaluateBField(
    const model::Vec3&                 pointAng,
    const model::RingGeometry&         geo,
    const std::vector<model::Vec3>&    verticesAng,
    double                             lobeOffsetA,
    double                             currentNA);

// Evaluate the geometric shielding kernel G at a point, decomposed
// into a SphericalTensor (T0 isotropic, T1 antisymmetric, T2 traceless
// symmetric). Implements:
//   G_ab = -n_b * B_a * PPM_FACTOR
// with the same spatial guards as EvaluateBField.
model::SphericalTensor EvaluateShielding(
    const model::Vec3&                 pointAng,
    const model::RingGeometry&         geo,
    const std::vector<model::Vec3>&    verticesAng,
    double                             lobeOffsetA,
    double                             intensityNA);

}  // namespace h5reader::calculators
