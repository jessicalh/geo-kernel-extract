#pragma once
//
// TrajectoryAtom: per-atom trajectory-scope data store. See
// OBJECT_MODEL.md (trajectory-scope) and PATTERNS.md §13 + Lesson 25.
// Private constructor; only TrajectoryProtein constructs via friend.
//
// Per-Welford state lives in named substructs (one per writer TR)
// each holding WelfordMoments instances per channel. See
// spec/plan/welford-data-shape-design-2026-05-17.md for the design
// rationale. The substruct shape replaced ~136 loose fields of the
// original pattern; differentiated structure compresses visual
// surface ~9x without changing memory layout.
//
// Phase 2a (2026-05-17): refactor only — existing channels in
// substructs, behavior unchanged. Phase 2b lands the per-component
// T1 / T2, drift / abs / rms delta variants, and schema provenance
// per the design doc.
//

#include "AtomEvent.h"
#include "RecordBag.h"
#include "TrajectoryMoments.h"  // WelfordMoments

#include <array>
#include <cstddef>

namespace nmr {

// Per-Welford state struct definitions. Each one bundles the channels
// owned by a single Welford TR:
//   - one WelfordMoments per channel
//   - n_frames + delta_n at struct level (shared across channels)
//   - Welford state is read by the owning TR's Compute / Finalize /
//     WriteH5Group and (for BS) by BsAnomalousAtomMarker cross-result
//     read at PATTERNS §17 marker discipline.

// Written by BsWelfordTrajectoryResult.
// Source: ConformationAtom::bs_shielding_contribution (SphericalTensor,
// units = ppm·T/nA per OBJECT_MODEL drift table).
//
// Phase 2b expansion (2026-05-17): per-component T1[3] + T2[5] preserve
// tensor orientation that |T2| amplitude rollup discards; per-channel
// frame-to-frame delta variants (signed = drift, |Δ| = abs-fluctuation,
// Δ² → sqrt = RMS fluctuation) distinguish drift from dynamics that
// the existing t0_delta telescope hides. Per PATTERNS Lesson 25.
struct BsWelfordState {
    // T0 isotropic scalar
    WelfordMoments t0;

    // T1 antisymmetric, rank-1 (G = -n⊗B has T1 = ½(n×B)); m = -1, 0, +1
    std::array<WelfordMoments, 3> t1;

    // T2 symmetric traceless; m = -2, -1, 0, +1, +2
    std::array<WelfordMoments, 5> t2;

    // |T2| Frobenius L2 amplitude (scalar summary; pairs with per-component t2)
    WelfordMoments t2magnitude;

    // Frame-to-frame T0 delta variants — three signals distinguished:
    // - t0_delta: signed Δ. mean = (x_N - x_0)/(N-1) telescopes — drift.
    // - t0_abs_delta: |Δ|. mean captures total per-frame path length.
    // - t0_delta_squared: Δ². mean → sqrt at Finalize = RMS fluctuation.
    WelfordMoments t0_delta;
    WelfordMoments t0_abs_delta;
    WelfordMoments t0_delta_squared;
    double         t0_rms_delta = 0.0;   // sqrt(t0_delta_squared.mean), Finalize-derived

    // Shared denominators
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by HmWelfordTrajectoryResult.
// Source: hm_shielding_contribution (Å⁻¹, rank-1 same as BS but no
// PPM_FACTOR multiplier per OBJECT_MODEL drift table).
// Phase 2b expansion: identical shape to BsWelfordState.
struct HmWelfordState {
    WelfordMoments t0;
    std::array<WelfordMoments, 3> t1;          // m = -1, 0, +1
    std::array<WelfordMoments, 5> t2;          // m = -2..+2
    WelfordMoments t2magnitude;
    WelfordMoments t0_delta;                   // signed Δ
    WelfordMoments t0_abs_delta;               // |Δ|
    WelfordMoments t0_delta_squared;           // Δ²
    double         t0_rms_delta = 0.0;         // sqrt(<Δ²>), Finalize-derived
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by McConnellWelfordTrajectoryResult.
// Source: mc_shielding_contribution (Å⁻³, full asymmetric non-traceless
// three-term McConnell form per PATTERNS Lesson 19).
// Phase 2b expansion: T1 SKIPPED — McConnell-form has T1 = 0 by construction.
struct McConnellWelfordState {
    WelfordMoments t0;
    // NO t1 — McConnell-form has T1 ≡ 0 (PATTERNS Lesson 19)
    std::array<WelfordMoments, 5> t2;          // m = -2..+2
    WelfordMoments t2magnitude;
    WelfordMoments t0_delta;
    WelfordMoments t0_abs_delta;
    WelfordMoments t0_delta_squared;
    double         t0_rms_delta = 0.0;
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by EeqWelfordTrajectoryResult.
// Source: eeq_charge (double, elementary_charge).
// Phase 2b expansion: scalar source — no T1/T2 channels.
struct EeqWelfordState {
    WelfordMoments charge;
    WelfordMoments charge_delta;               // signed Δ
    WelfordMoments charge_abs_delta;           // |Δ|
    WelfordMoments charge_delta_squared;       // Δ²
    double         charge_rms_delta = 0.0;     // sqrt(<Δ²>), Finalize-derived
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by SasaWelfordTrajectoryResult.
// Source: atom_sasa (double, Å²).
struct SasaWelfordState {
    WelfordMoments sasa;
    WelfordMoments sasa_delta;
    WelfordMoments sasa_abs_delta;
    WelfordMoments sasa_delta_squared;
    double         sasa_rms_delta = 0.0;
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by HBondCountWelfordTrajectoryResult.
// Source: hbond_count_within_3_5A (int, promoted to double).
// Phase 2b expansion: scalar source + occupancy_fraction companion
// (Welford on indicator: count > 0 ? 1.0 : 0.0) captures binary
// participation as a separate physics question from ⟨N⟩ expected count.
struct HBondCountWelfordState {
    WelfordMoments count;
    WelfordMoments count_delta;
    WelfordMoments count_abs_delta;
    WelfordMoments count_delta_squared;
    double         count_rms_delta = 0.0;
    WelfordMoments occupancy_fraction;         // Welford on (count > 0 ? 1.0 : 0.0)
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};


class TrajectoryProtein;

class TrajectoryAtom {
    friend class TrajectoryProtein;
public:
    // Per-Welford state substructs (one writer per substruct).
    BsWelfordState         bs_welford;
    HmWelfordState         hm_welford;
    McConnellWelfordState  mc_welford;
    EeqWelfordState        eeq_welford;
    SasaWelfordState       sasa_welford;
    HBondCountWelfordState hbond_count_welford;

    // Pattern C — per-atom event bag.
    // Push via events.Push({emitter, kind, frame, time, metadata}) from
    // any TrajectoryResult::Compute or ::Finalize that emits per-atom
    // events at this atom. Queries are the standard RecordBag verbs.
    // Used by BsAnomalousAtomMarker (BsAnomalyHighT0 / BsAnomalyLowT0).
    RecordBag<AtomEvent> events;

private:
    explicit TrajectoryAtom() = default;
};

}  // namespace nmr
