#pragma once
//
// KernelDynamicsTrajectoryResult: per-atom autocorrelation and power
// spectrum of the classical geometric shielding kernels, across the whole
// protein. A protein under MD integrates a coupled differential equation;
// each kernel scalar at each atom is a coordinate of that system that
// rings in time. This TR records the ringing -- the autocorrelation (how
// long the coordinate remembers itself) and its Fourier partner the power
// spectrum (at what frequencies it oscillates). It is an observability
// instrument over the kernels, the move IDENTITY_AND_DYNAMICS_ROLLUP 13.5
// places before any modelling choice.
//
// Channels (per atom, per frame): the per-atom *_shielding_contribution
// SphericalTensors reduced two honest ways -- T0 where the kernel carries
// a physical isotropic part, |T2| = SphericalTensor::T2Magnitude() for
// all. piquad / disp / apbs(EFG) are structurally traceless (T0 == 0 by
// Laplace / symmetry), so only their |T2| channel is emitted.
//
// Lifecycle: FO. Per (atom, channel) an exact bounded-memory biased
// autocorrelation (TrajectorySpectral.h) accumulates each frame; Finalize
// turns each into rho(k), a Parzen power spectrum S(f) (guaranteed >= 0),
// and three scalar reductions. The rho(k) and S(f) curves are the product;
// the scalars are captions that point back at the curves.
//
// Storage: the finalized arrays are result-owned vectors written directly
// in WriteH5Group (the DihedralTimeSeries idiom). NOT DenseBuffer -- that
// is keyed one buffer per TR type, so the several double datasets this TR
// emits cannot share it without overwriting.
//
// Emission /trajectory/kernel_dynamics/:
//   acf                       (N, C, L) float64  rho(k) in [-1, 1]
//   power_spectrum            (N, C, F) float64  S(f) >= 0, Parzen PSD
//   decay_time_ps             (N, C)    float64  dt * sum rho to first zero crossing
//   peak_freq_per_ps          (N, C)    float64  argmax frequency (excl. DC bin)
//   spectral_centroid_per_ps  (N, C)    float64  sum(f S) / sum(S)
//   channel_names             (C,) string ; channel_units (C,) string
//   lag_frames (L,) uint64 ; lag_times_ps (L,) float64 ; frequencies_per_ps (F,) float64
//   attrs: result_name, n_atoms, n_channels, n_lags, n_freq, n_frames,
//          finalized, sample_interval_ps, estimator="biased",
//          mean_convention="full_range", window="parzen", spectrum_units
//
// A constant (or near-constant) signal has C(0) ~ 0: rho and spectrum are
// emitted as 0 and the reductions as NaN (the honest "no oscillation"
// answer), distinguishable downstream via isfinite() on the reductions.
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

class KernelDynamicsTrajectoryResult : public TrajectoryResult {
public:
    // The lag count is the CalculatorConfig parameter `dynamics_n_lags`
    // (default 120), read at Create into n_lags_; the one-sided
    // frequency-bin count is n_freq_ = n_lags_ + 1 (f_m = m/(2 L dt),
    // m = 0..L spans DC..Nyquist). The reduced channels are fixed below.

    // The reduced scalar channels, in emission order. T0 channels carry
    // the kernel's isotropic part; AbsT2 channels carry the angular
    // amplitude |T2|. Pure-T2 kernels (piquad/disp/apbs) have no T0 channel.
    enum class Channel : std::uint8_t {
        BsT0, BsAbsT2, HmT0, HmAbsT2, McT0, McAbsT2,
        RingChiT0, RingChiAbsT2, HBondT0, HBondAbsT2,
        PiQuadAbsT2, DispAbsT2, ApbsAbsT2
    };
    static constexpr std::size_t N_CHANNELS = 13;

    std::string Name() const override {
        return "KernelDynamicsTrajectoryResult";
    }

    // The source ConformationResults whose per-atom fields are read. All
    // are required by PerFrameExtractionSet; declaring them lets Phase 4
    // reject a configuration that cannot feed this instrument.
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<KernelDynamicsTrajectoryResult> Create(
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
    // One biased-ACF accumulator per (atom, channel), flat index
    // atom * N_CHANNELS + channel. Released after Finalize consumes it.
    std::vector<BiasedAcfAccumulator> accumulators_;
    std::size_t n_atoms_ = 0;
    std::size_t n_lags_ = 120;   // CalculatorConfig dynamics_n_lags, set at Create
    std::size_t n_freq_ = 121;   // n_lags_ + 1

    // Per-Compute frame times, for the median sample interval at Finalize
    // (self-recorded, not read from traj, so manual callers work -- same
    // rationale as BsT0Autocorrelation).
    std::vector<double> frame_times_;

    // Finalized, result-owned. acf_/spectrum_ flat (atom, channel, lag|freq);
    // the three reductions flat (atom, channel).
    std::vector<double> acf_;                       // N * C * n_lags_
    std::vector<double> spectrum_;                  // N * C * n_freq_
    std::vector<double> decay_time_ps_;             // N * C
    std::vector<double> peak_freq_per_ps_;          // N * C
    std::vector<double> spectral_centroid_per_ps_;  // N * C

    std::size_t n_frames_ = 0;
    double sample_interval_ps_ = 0.0;
    bool finalized_ = false;
};

}  // namespace nmr
