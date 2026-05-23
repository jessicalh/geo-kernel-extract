// QtWelfordBuffers.h — atom-axis Welford rollup buffers.
//
// Welford rollup TRs emit per-atom statistics across all frames (no T
// axis on the data — T is collapsed into mean/std/m2/min/max + min/max
// frame indices). Per design §11 the rollup is rich: per-component
// Welford on T1[3] and T2[5], plus delta variants (drift, |Δ|, Δ², dx/dt)
// and per-frame-count metadata.
//
// Storage strategy: one struct per Welford TR family, exposing typed
// rollup accessors. The shielding family (BS, HM, MC) shares the same
// shape via QtShieldingWelford; scalar-source Welfords (SASA, EEQ,
// HBondCount, MopacCharge, WaterField scalar) share QtScalarWelford.
// Specialised shapes (BondOrder Welford, Hydration Welford, AIMNet2
// charge-response Welford) get their own struct.

#pragma once

#include "Types.h"

#include <QString>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::model {

// ──────────────────────────────────────────────────────────────────
// QtWelfordMoments — scalar Welford moments + min/max + frame indices.
//
// Mirrors the per-channel rollup the writer emits: mean, std (sample
// ddof=1), m2 (second moment), min, max, min_frame, max_frame.
// ──────────────────────────────────────────────────────────────────

struct QtWelfordMoments {
    double mean = 0.0;
    double std = 0.0;
    double m2 = 0.0;
    double min = 0.0;
    double max = 0.0;
    uint64_t min_frame = 0;
    uint64_t max_frame = 0;
};


// ──────────────────────────────────────────────────────────────────
// QtShieldingWelford — per-atom rollup for shielding TRs (BS, HM, MC).
//
// Source dataset family in fixture (~75 datasets per TR group):
//   - t0_mean / t0_std / t0_m2 / t0_min / t0_max / t0_min_frame / t0_max_frame
//   - t1_mean (N,3) / t1_std (N,3) / t1_m2 (N,3) / ... (same 7 stats × 3 components)
//   - t2_mean (N,5) / t2_std (N,5) / ... (same 7 stats × 5 components)
//   - t2magnitude_mean / _std / _m2 / _min / _max / _min_frame / _max_frame
//   - t0_delta_mean / _std / _m2 / ... (signed Δ; 7 stats)
//   - t0_abs_delta_mean / ... (|Δ|; 7 stats)
//   - t0_delta_squared_mean / ... (Δ²; 7 stats)
//   - t0_dxdt_mean / ... (Δ/Δt; 7 stats)
//   - t0_rms_delta (derived at Finalize from t0_delta_squared.mean)
//   - n_frames_per_atom (N,) uint64
//   - delta_n_per_atom (N,) uint64
//   - dxdt_n_per_atom (N,) uint64
// ──────────────────────────────────────────────────────────────────

struct QtShieldingWelford {
    std::size_t n_atoms = 0;

    // T0 scalar Welford
    std::vector<QtWelfordMoments> t0;  // (N,)

    // T1 per-component (Cartesian or m-basis per irrep_layout_t1 attr)
    std::vector<std::array<QtWelfordMoments, 3>> t1;  // (N,)

    // T2 per-component (m=-2..+2 per irrep_layout_t2 attr)
    std::vector<std::array<QtWelfordMoments, 5>> t2;  // (N,)

    // |T2| Frobenius amplitude Welford
    std::vector<QtWelfordMoments> t2magnitude;  // (N,)

    // Delta variants on T0
    std::vector<QtWelfordMoments> t0_delta;          // signed Δ
    std::vector<QtWelfordMoments> t0_abs_delta;      // |Δ|
    std::vector<QtWelfordMoments> t0_delta_squared;  // Δ²
    std::vector<QtWelfordMoments> t0_dxdt;           // Δ/Δt
    std::vector<double> t0_rms_delta;                // derived at Finalize

    // Frame-count metadata
    std::vector<uint64_t> n_frames_per_atom;
    std::vector<uint64_t> delta_n_per_atom;
    std::vector<uint64_t> dxdt_n_per_atom;

    // Attrs — display
    QString units;
    QString result_name;
    QString irrep_layout_t1;  // "v_x,v_y,v_z"
    QString irrep_layout_t2;  // "m-2,m-1,m0,m+1,m+2"
    int ddof = 1;             // sample stddev divisor
    double mean_dt_ps = 20.0;
    std::array<uint64_t, 2> frame_index_range = {0, 0};
};


// ──────────────────────────────────────────────────────────────────
// QtScalarWelford — per-atom rollup for scalar-source Welfords.
// Used by SASA, EEQ charge, H-bond count, MOPAC charge, water field
// scalar variants. Schema mirrors QtShieldingWelford's t0-only subset
// (no T1, T2, or magnitude — scalar source has only the scalar channel).
//
// The writer emits dataset names like `<channel>_mean` (canonical) and
// some older fixtures used `mean` (legacy alias). The loader accepts
// either at the H5 boundary.
// ──────────────────────────────────────────────────────────────────

struct QtScalarWelford {
    std::size_t n_atoms = 0;

    std::vector<QtWelfordMoments> value;  // (N,)
    std::vector<QtWelfordMoments> delta;  // signed Δ
    std::vector<QtWelfordMoments> abs_delta;
    std::vector<QtWelfordMoments> delta_squared;
    std::vector<QtWelfordMoments> dxdt;
    std::vector<double> rms_delta;  // derived at Finalize

    std::vector<uint64_t> n_frames_per_atom;
    std::vector<uint64_t> delta_n_per_atom;
    std::vector<uint64_t> dxdt_n_per_atom;

    QString units;
    QString result_name;
    QString channel_name;  // e.g. "sasa", "charge", "count"
    int ddof = 1;
    double mean_dt_ps = 20.0;
};


// ──────────────────────────────────────────────────────────────────
// QtVec3Welford — per-atom Vec3 Welford (per-component rollup).
// Used by water_field_welford and aimnet2_charge_response_gradient_welford
// (when source is a vector).
// ──────────────────────────────────────────────────────────────────

struct QtVec3Welford {
    std::size_t n_atoms = 0;

    std::vector<std::array<QtWelfordMoments, 3>> components;  // (N,)
    std::vector<QtWelfordMoments> magnitude;                  // (N,)

    std::vector<uint64_t> n_frames_per_atom;

    QString units;
    QString result_name;
    QString irrep_layout;  // "x,y,z"
    int ddof = 1;
};


// ──────────────────────────────────────────────────────────────────
// QtBondOrderWelford — per-bond rollup for MOPAC Wiberg bond orders.
// Bond-axis (not atom-axis). Conditional-attach per MOPAC cadence.
// ──────────────────────────────────────────────────────────────────

struct QtBondOrderWelford {
    std::size_t n_bonds = 0;

    std::vector<QtWelfordMoments> bond_order;  // (B,)
    std::vector<uint64_t> n_frames_per_bond;

    QString units;
    QString result_name;
    int ddof = 1;
};


// ──────────────────────────────────────────────────────────────────
// QtHydrationWelford — per-atom hydration Welford (multi-channel).
// Used by hydration_shell_welford and hydration_geometry_welford.
// Channels: shell counts, distances, angles, etc. — variable per TR;
// loader populates whichever datasets are present.
// ──────────────────────────────────────────────────────────────────

struct QtHydrationChannel {
    QString name;                           // e.g. "first_shell_count"
    std::vector<QtWelfordMoments> moments;  // (N,)
    QString units;
};

struct QtHydrationWelford {
    std::size_t n_atoms = 0;

    std::vector<QtHydrationChannel> channels;
    std::vector<uint64_t> n_frames_per_atom;

    QString result_name;
    int ddof = 1;
};


// ──────────────────────────────────────────────────────────────────
// QtAutocorrelation — per-atom autocorrelation at fixed lag grid.
// Used by bs_t0_autocorrelation (BS T0 autocorr, ~120 lags).
// ──────────────────────────────────────────────────────────────────

struct QtAutocorrelation {
    std::size_t n_atoms = 0;
    std::size_t n_lags = 0;

    std::vector<double> rho;           // (N*n_lags,) row-major
    std::vector<double> lag_times_ps;  // (n_lags,)

    QString units;
    QString result_name;

    double at(std::size_t atomIdx, std::size_t lag) const {
        if (atomIdx >= n_atoms || lag >= n_lags)
            return 0.0;
        return rho[atomIdx * n_lags + lag];
    }
};


}  // namespace h5reader::model
