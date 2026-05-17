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
struct BsWelfordState {
    WelfordMoments t0;             // BiotSavart isotropic shielding scalar
    WelfordMoments t2magnitude;    // |T2| Frobenius amplitude
    WelfordMoments t0_delta;       // frame-to-frame T0 difference (signed)
    std::size_t    n_frames = 0;   // primary denominator (T0 + |T2|)
    std::size_t    delta_n  = 0;   // delta-sample count (= n_frames - 1)
};

// Written by HmWelfordTrajectoryResult.
// Source: hm_shielding_contribution (Å⁻¹, rank-1 same as BS but no
// PPM_FACTOR multiplier per OBJECT_MODEL drift table).
struct HmWelfordState {
    WelfordMoments t0;
    WelfordMoments t2magnitude;
    WelfordMoments t0_delta;
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by McConnellWelfordTrajectoryResult.
// Source: mc_shielding_contribution (Å⁻³, full asymmetric non-traceless
// three-term McConnell form per PATTERNS Lesson 19).
struct McConnellWelfordState {
    WelfordMoments t0;
    WelfordMoments t2magnitude;
    WelfordMoments t0_delta;
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by EeqWelfordTrajectoryResult.
// Source: eeq_charge (double, elementary_charge).
struct EeqWelfordState {
    WelfordMoments charge;
    WelfordMoments charge_delta;
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by SasaWelfordTrajectoryResult.
// Source: atom_sasa (double, Å²).
struct SasaWelfordState {
    WelfordMoments sasa;
    WelfordMoments sasa_delta;
    std::size_t    n_frames = 0;
    std::size_t    delta_n  = 0;
};

// Written by HBondCountWelfordTrajectoryResult.
// Source: hbond_count_within_3_5A (int, promoted to double).
struct HBondCountWelfordState {
    WelfordMoments count;
    WelfordMoments count_delta;
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
