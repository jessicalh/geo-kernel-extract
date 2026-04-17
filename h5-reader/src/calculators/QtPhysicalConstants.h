// QtPhysicalConstants.h — constants the reader's kernel evaluators need.
//
// These are the subset of nmr-shielding's PhysicalConstants.h + TOML
// numerical guards that the volumetric BS/HM isosurface pipeline
// depends on. Values are copied verbatim so the reader reproduces the
// library's sign conventions and numerical behaviour bit-for-bit at
// open-space grid points.
//
// If any of these drift from the library, butterflies rendered here
// won't match what BiotSavartResult::SampleShieldingAt produces on
// the same geometry. Keep synchronised when library values change.

#pragma once

#include <cmath>

namespace h5reader::calculators {

// ---- Unit conversions (SI boundary) ----------------------------------
inline constexpr double ANGSTROMS_TO_METRES   = 1.0e-10;
inline constexpr double NANOAMPERES_TO_AMPERES = 1.0e-9;
inline constexpr double PPM_FACTOR            = 1.0e6;

// ---- Biot-Savart ----------------------------------------------------
// Pre-2019 SI: mu_0/(4*pi) is exactly 1e-7 T*m/A. DO NOT substitute
// the 2019 CODATA-derived value — it differs in the 10th digit and
// breaks binary reproduction of the library's BS output.
inline constexpr double BIOT_SAVART_PREFACTOR = 1e-7;

// Numerical guards from calculator_params.toml:
//   biot_savart_wire_endpoint_guard = 1e-25   (SI metres — field point at segment endpoint)
//   biot_savart_wire_axis_guard     = 1e-70   (SI metres² — field point on the wire axis)
inline constexpr double BS_WIRE_ENDPOINT_GUARD = 1e-25;
inline constexpr double BS_WIRE_AXIS_GUARD     = 1e-70;

// ---- Haigh-Mallion --------------------------------------------------
//   haigh_mallion_triangle_area_guard = 1e-10  (Å² — skip degenerate triangles)
//   haigh_mallion_subdivision_threshold_l1 = 2.0  (Å — subdivide at level 0→1)
//   haigh_mallion_subdivision_threshold_l2 = 1.0  (Å — subdivide at level 1→2)
inline constexpr double HM_TRIANGLE_AREA_GUARD = 1e-10;
inline constexpr double HM_SUBDIV_THRESHOLD_L1 = 2.0;
inline constexpr double HM_SUBDIV_THRESHOLD_L2 = 1.0;

// ---- Shared spatial guards -------------------------------------------
//   singularity_guard_distance      = 0.1 Å   (field-point proximity cutoff)
//   ring_current_spatial_cutoff     = 15.0 Å  (BS/HM decay length)
inline constexpr double SINGULARITY_GUARD_DISTANCE = 0.1;
inline constexpr double RING_CURRENT_CUTOFF        = 15.0;

// ---- Grid-sampling parameters (viewer choice, not physics) -----------
// Grid extent from ring center, in Angstroms. The library viewer uses
// 7 Å for the T0 field grid and 6 Å for the butterfly B-field grid.
// We use 7 for the isosurface grid.
inline constexpr double FIELD_GRID_EXTENT_A = 7.0;

// Grid resolution per axis. 20³ = 8000 evaluations per ring per frame.
// At typical BS cost ~1 μs/eval, that's ~8 ms per ring per frame —
// well under the 200 ms budget for 5 fps playback with a handful of rings.
inline constexpr int FIELD_GRID_DIM = 20;

}  // namespace h5reader::calculators
