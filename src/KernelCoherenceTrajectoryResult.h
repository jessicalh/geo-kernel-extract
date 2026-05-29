#pragma once
//
// KernelCoherenceTrajectoryResult: per-atom, zero-lag correlation matrix
// between the geometric shielding-kernel channels -- "which kernels move
// together at this atom." The relationships lens over the same coupled
// system KernelDynamics autocorrelates (IDENTITY_AND_DYNAMICS_ROLLUP 13.5,
// kernel-kernel covariance as signal). The chain, rings, disulfides, and
// packing impose correlations on atom motion; this surfaces them at each
// atom across the kernels.
//
// v1 is the zero-lag Pearson matrix r_xy = cov_xy / sqrt(var_x var_y).
// The lagged (lead/lag) cross-correlation -- which kernel anticipates
// which -- is a documented follow-up: it needs a second bounded
// accumulator and is better landed carefully than rushed.
//
// Same 13 channels as KernelDynamicsTrajectoryResult, cloned here rather
// than shared (PATTERNS.md 17: each TR self-contained, testable in
// isolation, safe to delete). FO lifecycle: per atom, running sums of
// each channel, its square, and each channel-pair product accumulate each
// frame; Finalize forms the matrix.
//
// Emission /trajectory/kernel_coherence/:
//   correlation_matrix (N, C, C) float64  symmetric, diagonal 1.0; NaN in
//                                          a row/column whose channel is
//                                          constant over the run (var ~ 0)
//   channel_names (C,) string ; channel_units (C,) string
//   attrs: result_name, n_atoms, n_channels, n_frames, finalized,
//          statistic="pearson_zero_lag", lagged_cross_correlation="deferred"
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class KernelCoherenceTrajectoryResult : public TrajectoryResult {
public:
    // Cloned from KernelDynamicsTrajectoryResult::Channel (same kernels,
    // same emission order) -- the two TRs stay independent.
    enum class Channel : std::uint8_t {
        BsT0, BsAbsT2, HmT0, HmAbsT2, McT0, McAbsT2,
        RingChiT0, RingChiAbsT2, HBondT0, HBondAbsT2,
        PiQuadAbsT2, DispAbsT2, ApbsAbsT2
    };
    static constexpr std::size_t N_CHANNELS = 13;

    std::string Name() const override {
        return "KernelCoherenceTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<KernelCoherenceTrajectoryResult> Create(
        const TrajectoryProtein& tp);

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
    // Running sums, flat per atom. sum_xy_ holds the upper triangle
    // including the diagonal, in (a <= b) row-major order (N_PAIRS per
    // atom). The full symmetric matrix is reconstructed at Finalize.
    static constexpr std::size_t N_PAIRS = N_CHANNELS * (N_CHANNELS + 1) / 2;
    std::vector<double> sum_x_;   // N * C
    std::vector<double> sum_xx_;  // N * C
    std::vector<double> sum_xy_;  // N * N_PAIRS

    std::vector<double> correlation_matrix_;  // N * C * C (finalized)
    std::size_t n_atoms_  = 0;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr
