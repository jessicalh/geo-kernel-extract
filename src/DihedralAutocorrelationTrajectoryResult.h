#pragma once
//
// DihedralAutocorrelationTrajectoryResult: per-residue circular
// autocorrelation of the backbone phi/psi and sidechain chi torsions --
// the torsional decorrelation timescale ("the clicking"). Complements
// DihedralBinTransitionTrajectoryResult, which counts rotamer transitions
// but carries no timescale, and DihedralTimeSeriesTrajectoryResult, which
// emits the raw angle series.
//
// Circular autocorrelation (TrajectorySpectral.h CircularAcfAccumulator):
//   C(k) = < cos( theta(t+k) - theta(t) ) >, the proper periodic form
// (the naive < theta(t) theta(0) > is wrong across the +-pi branch cut).
// C(0) = 1; the long-time value is the circular order < cos >^2 + < sin >^2.
//
// Angles, per residue: phi = C(i-1)-N-CA-C, psi = N-CA-C-N(i+1), and
// chi[0..3] from Residue.chi[k]. Same IUPAC dihedral and bond-graph
// backbone-adjacency conventions as DihedralTimeSeriesTrajectoryResult,
// whose Dihedral() helper and Predecessor/Successor walk are cloned here
// (PATTERNS.md 17). Omega is omitted -- a near-rigid trans bond carries
// little decorrelation signal; DihedralTimeSeries keeps the omega series.
//
// Deferred (noted in the design doc and group attrs), to land carefully
// rather than rushed: a torsional power spectrum, and a rotamer-state
// survival ACF (jump time) on top of the transition counts.
//
// Lifecycle: FO. A circular-ACF accumulator per (residue, angle); Finalize
// produces C(k) and a 1/e decorrelation time. Structurally undefined
// angles (terminus phi, residue without chi[k]) emit NaN curves and NaN
// times -- isfinite() distinguishes. A non-finite per-frame dihedral
// (degenerate geometry, vanishingly rare for bonded atoms) holds the last
// finite value so the lag-to-time mapping stays intact.
//
// Emission /trajectory/dihedral_autocorrelation/:
//   phi_acf, psi_acf            (R, L)    float64  C(k) in [-1, 1]
//   chi_acf                     (R, 4, L) float64
//   phi_corr_time_ps, psi_corr_time_ps (R,)    float64  1/e decorrelation time
//   chi_corr_time_ps            (R, 4)    float64
//   phi_defined, psi_defined    (R,)   uint8 ; chi_defined (R, 4) uint8
//   residue_index_per_atom      (N,)   int32   atom-axis broadcast
//   lag_frames (L,) uint64 ; lag_times_ps (L,) float64
//   attrs: result_name, n_residues, n_atoms, n_frames, finalized,
//          sample_interval_ps, estimator, correlation_time_definition,
//          angle_convention, deferred
//

#include "TrajectoryResult.h"
#include "TrajectorySpectral.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class DihedralAutocorrelationTrajectoryResult : public TrajectoryResult {
public:
    // Lag count from the CalculatorConfig parameter `dynamics_n_lags`
    // (default 120), read at Create into n_lags_.

    std::string Name() const override {
        return "DihedralAutocorrelationTrajectoryResult";
    }

    // No declared dependency: positions and Residue.chi[k] are present
    // after tp.Seed (PATTERNS.md 15), like DihedralTimeSeries.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<DihedralAutocorrelationTrajectoryResult>
    Create(const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumFrames() const { return n_frames_; }

private:
    // One circular-ACF accumulator per residue per angle.
    std::vector<CircularAcfAccumulator> phi_acc_;       // R
    std::vector<CircularAcfAccumulator> psi_acc_;       // R
    std::vector<CircularAcfAccumulator> chi_acc_;       // R * 4

    // Static "structurally defined" masks (R / R*4), set at Create.
    std::vector<std::uint8_t> phi_defined_;
    std::vector<std::uint8_t> psi_defined_;
    std::vector<std::uint8_t> chi_defined_;             // R * 4

    // Hold state for the rare non-finite frame (carry last finite value).
    std::vector<double> phi_last_, psi_last_;           // R
    std::vector<double> chi_last_;                      // R * 4
    std::vector<std::uint8_t> phi_has_last_, psi_has_last_;
    std::vector<std::uint8_t> chi_has_last_;

    std::vector<std::int32_t> residue_index_per_atom_;  // N

    // Finalized, result-owned.
    std::vector<double> phi_acf_, psi_acf_;             // R * L
    std::vector<double> chi_acf_;                       // R * 4 * L
    std::vector<double> phi_corr_time_, psi_corr_time_; // R
    std::vector<double> chi_corr_time_;                 // R * 4

    std::vector<double> frame_times_;                   // per frame, for dt
    std::size_t n_lags_     = 120;  // CalculatorConfig dynamics_n_lags
    std::size_t n_residues_ = 0;
    std::size_t n_atoms_    = 0;
    std::size_t n_frames_   = 0;
    double sample_interval_ps_ = 0.0;
    bool finalized_ = false;
};

}  // namespace nmr
