#include "KernelDynamicsTrajectoryResult.h"

#include "ApbsFieldResult.h"
#include "BiotSavartResult.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DispersionResult.h"
#include "HBondResult.h"
#include "HaighMallionResult.h"
#include "McConnellResult.h"
#include "OperationLog.h"
#include "PiQuadrupoleResult.h"
#include "ProteinConformation.h"
#include "RingSusceptibilityResult.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <typeinfo>

namespace nmr {

namespace {

using Channel = KernelDynamicsTrajectoryResult::Channel;

// The scalar each channel reads from a ConformationAtom. Typed field
// access on the per-atom shielding-contribution SphericalTensors -- T0 is
// the isotropic part, T2Magnitude() the angular amplitude |T2|. No string
// dispatch: the enum is the channel identity.
double ChannelValue(Channel c, const ConformationAtom& a) {
    switch (c) {
        case Channel::BsT0:         return a.bs_shielding_contribution.T0;
        case Channel::BsAbsT2:      return a.bs_shielding_contribution.T2Magnitude();
        case Channel::HmT0:         return a.hm_shielding_contribution.T0;
        case Channel::HmAbsT2:      return a.hm_shielding_contribution.T2Magnitude();
        case Channel::McT0:         return a.mc_shielding_contribution.T0;
        case Channel::McAbsT2:      return a.mc_shielding_contribution.T2Magnitude();
        case Channel::RingChiT0:    return a.ringchi_shielding_contribution.T0;
        case Channel::RingChiAbsT2: return a.ringchi_shielding_contribution.T2Magnitude();
        case Channel::HBondT0:      return a.hbond_shielding_contribution.T0;
        case Channel::HBondAbsT2:   return a.hbond_shielding_contribution.T2Magnitude();
        case Channel::PiQuadAbsT2:  return a.piquad_shielding_contribution.T2Magnitude();
        case Channel::DispAbsT2:    return a.disp_shielding_contribution.T2Magnitude();
        case Channel::ApbsAbsT2:    return a.apbs_efg_spherical.T2Magnitude();
    }
    return 0.0;
}

const char* ChannelName(Channel c) {
    switch (c) {
        case Channel::BsT0:         return "bs_T0";
        case Channel::BsAbsT2:      return "bs_absT2";
        case Channel::HmT0:         return "hm_T0";
        case Channel::HmAbsT2:      return "hm_absT2";
        case Channel::McT0:         return "mc_T0";
        case Channel::McAbsT2:      return "mc_absT2";
        case Channel::RingChiT0:    return "ringchi_T0";
        case Channel::RingChiAbsT2: return "ringchi_absT2";
        case Channel::HBondT0:      return "hbond_T0";
        case Channel::HBondAbsT2:   return "hbond_absT2";
        case Channel::PiQuadAbsT2:  return "piquad_absT2";
        case Channel::DispAbsT2:    return "disp_absT2";
        case Channel::ApbsAbsT2:    return "apbs_absT2";
    }
    return "unknown";
}

// Natural kernel units (bare decomposed kernel, not ppm), matching the
// SDK catalog unit strings for these fields. T0 and |T2| of one kernel
// share the kernel's units.
const char* ChannelUnits(Channel c) {
    switch (c) {
        case Channel::BsT0:
        case Channel::BsAbsT2:      return "ppm_T_per_nA";
        case Channel::HmT0:
        case Channel::HmAbsT2:      return "Angstrom^-1";
        case Channel::McT0:
        case Channel::McAbsT2:
        case Channel::RingChiT0:
        case Channel::RingChiAbsT2:
        case Channel::HBondT0:
        case Channel::HBondAbsT2:   return "Angstrom^-3";
        case Channel::PiQuadAbsT2:  return "Angstrom^-5";
        case Channel::DispAbsT2:    return "Angstrom^-6";
        case Channel::ApbsAbsT2:    return "V/A^2";
    }
    return "";
}

}  // namespace


std::vector<std::type_index>
KernelDynamicsTrajectoryResult::Dependencies() const {
    return {
        std::type_index(typeid(BiotSavartResult)),
        std::type_index(typeid(HaighMallionResult)),
        std::type_index(typeid(McConnellResult)),
        std::type_index(typeid(RingSusceptibilityResult)),
        std::type_index(typeid(PiQuadrupoleResult)),
        std::type_index(typeid(DispersionResult)),
        std::type_index(typeid(HBondResult)),
        std::type_index(typeid(ApbsFieldResult)),
    };
}


std::unique_ptr<KernelDynamicsTrajectoryResult>
KernelDynamicsTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<KernelDynamicsTrajectoryResult>();
    // Clamp to >= 1: a misconfigured dynamics_n_lags <= 0 must not size the
    // accumulators to zero (Finalize would then read an empty covariance
    // vector). codex review 2026-05-29.
    const double n_lags_raw = CalculatorConfig::Get("dynamics_n_lags");
    const std::size_t n_lags =
        (n_lags_raw >= 1.0) ? static_cast<std::size_t>(n_lags_raw) : 1;
    r->n_lags_ = n_lags;
    r->n_freq_ = n_lags + 1;
    r->n_atoms_ = tp.AtomCount();
    r->accumulators_.assign(r->n_atoms_ * N_CHANNELS,
                            BiasedAcfAccumulator(n_lags));
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Push each channel's scalar at each atom into its own bounded-memory
// biased-ACF accumulator. Record the frame time so Finalize can derive
// the sample interval without reading traj.FrameTimes() (unfilled for
// manually-orchestrated callers -- same rationale as BsT0Autocorrelation).

void KernelDynamicsTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)frame_idx;
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        const ConformationAtom& a = conf.AtomAt(i);
        for (std::size_t ch = 0; ch < N_CHANNELS; ++ch) {
            accumulators_[i * N_CHANNELS + ch].Push(
                ChannelValue(static_cast<Channel>(ch), a));
        }
    }
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// For each (atom, channel): C(k) from the accumulator, then rho(k) =
// C(k)/C(0), the Parzen power spectrum from C(k), and the three scalar
// reductions. Constant signal (C(0) ~ 0) or < 2 frames -> curves left 0,
// reductions NaN. The per-(atom,channel) accumulators are released after
// the sweep.

void KernelDynamicsTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                              Trajectory& traj) {
    (void)traj;
    if (finalized_) return;
    n_atoms_ = tp.AtomCount();
    const std::size_t N = n_atoms_;
    const std::size_t C = N_CHANNELS;
    const double nan_val = std::nan("");

    // Median consecutive dt from our own frame_times_ (0 if < 2 frames).
    sample_interval_ps_ = 0.0;
    if (frame_times_.size() >= 2) {
        std::vector<double> deltas;
        deltas.reserve(frame_times_.size() - 1);
        for (std::size_t t = 1; t < frame_times_.size(); ++t) {
            deltas.push_back(frame_times_[t] - frame_times_[t - 1]);
        }
        std::sort(deltas.begin(), deltas.end());
        sample_interval_ps_ = deltas[deltas.size() / 2];
    }

    acf_.assign(N * C * n_lags_, 0.0);
    spectrum_.assign(N * C * n_freq_, 0.0);
    decay_time_ps_.assign(N * C, nan_val);
    peak_freq_per_ps_.assign(N * C, nan_val);
    spectral_centroid_per_ps_.assign(N * C, nan_val);

    const std::vector<double> freqs =
        SpectrumFrequenciesPerPs(n_lags_, sample_interval_ps_, n_freq_);
    const bool dt_ok = sample_interval_ps_ > 0.0;

    for (std::size_t flat = 0; flat < N * C; ++flat) {
        const std::vector<double> cov = accumulators_[flat].Finalize();
        const double c0 = cov[0];
        // No oscillation to report: < 2 frames, a constant signal (c0 ~ 0),
        // or a non-finite c0. The last shouldn't occur -- the source
        // calculators sanitise NaN/Inf (PATTERNS "absent, not faked") -- but
        // the explicit isfinite() keeps a NaN channel from masquerading as a
        // clean constant via NaN-compare fall-through (codex review 2026-05-29).
        if (n_frames_ < 2 || !std::isfinite(c0) || !(c0 > 1e-15)) {
            continue;  // curves stay 0, reductions stay NaN
        }

        double* acf_slice = &acf_[flat * n_lags_];
        for (std::size_t k = 0; k < n_lags_; ++k) acf_slice[k] = cov[k] / c0;

        const std::vector<double> S =
            ParzenPowerSpectrum(cov, sample_interval_ps_, n_freq_);
        double* spec_slice = &spectrum_[flat * n_freq_];
        for (std::size_t m = 0; m < n_freq_; ++m) spec_slice[m] = S[m];

        if (dt_ok) {
            // Correlation time: dt * sum of rho up to the first zero
            // crossing (rho[0] = 1). Full window if no crossing -> a
            // lower bound (documented in the dataset note).
            double tau = 0.0;
            for (std::size_t k = 0; k < n_lags_; ++k) {
                if (acf_slice[k] <= 0.0) break;
                tau += acf_slice[k];
            }
            decay_time_ps_[flat] = tau * sample_interval_ps_;

            std::size_t argmax = 1;
            double best = -std::numeric_limits<double>::infinity();
            for (std::size_t m = 1; m < n_freq_; ++m) {
                if (S[m] > best) { best = S[m]; argmax = m; }
            }
            peak_freq_per_ps_[flat] = freqs[argmax];

            double num = 0.0, den = 0.0;
            for (std::size_t m = 0; m < n_freq_; ++m) {
                num += freqs[m] * S[m];
                den += S[m];
            }
            spectral_centroid_per_ps_[flat] =
                (den > 0.0) ? (num / den) : nan_val;
        }
    }

    std::vector<BiasedAcfAccumulator>().swap(accumulators_);
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "KernelDynamicsTrajectoryResult::Finalize",
        "ACF + Parzen power spectrum across " + std::to_string(N) +
        " atoms x " + std::to_string(C) + " channels, " +
        std::to_string(n_frames_) + " frames, " +
        std::to_string(n_lags_) + " lags; sample interval ~ " +
        std::to_string(sample_interval_ps_) + " ps");
}


// ── WriteH5Group ─────────────────────────────────────────────────

void KernelDynamicsTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    if (!finalized_) {
        OperationLog::Warn("KernelDynamicsTrajectoryResult::WriteH5Group",
                           "Finalize not called; nothing to write");
        return;
    }
    const std::size_t N = n_atoms_;
    const std::size_t C = N_CHANNELS;

    auto grp = file.createGroup("/trajectory/kernel_dynamics");
    grp.createAttribute("result_name",        Name());
    grp.createAttribute("n_atoms",            N);
    grp.createAttribute("n_channels",         C);
    grp.createAttribute("n_lags",             static_cast<std::size_t>(n_lags_));
    grp.createAttribute("n_freq",             static_cast<std::size_t>(n_freq_));
    grp.createAttribute("n_frames",           n_frames_);
    grp.createAttribute("finalized",          finalized_);
    grp.createAttribute("sample_interval_ps", sample_interval_ps_);
    grp.createAttribute("estimator",          std::string("biased"));
    grp.createAttribute("mean_convention",    std::string("full_range"));
    grp.createAttribute("window",             std::string("parzen"));
    grp.createAttribute("spectrum_sidedness", std::string("one_sided_0_to_nyquist"));
    grp.createAttribute("spectrum_units",     std::string("channel_units^2 * ps"));
    grp.createAttribute("constant_signal_policy", std::string(
        "C(0) < 1e-15 (constant signal) or < 2 frames: acf and "
        "power_spectrum are 0, the three reductions NaN. Use "
        "isfinite(decay_time_ps) to distinguish 'no oscillation' from a "
        "real measurement."));

    std::vector<std::string> names(C), units(C);
    for (std::size_t ch = 0; ch < C; ++ch) {
        names[ch] = ChannelName(static_cast<Channel>(ch));
        units[ch] = ChannelUnits(static_cast<Channel>(ch));
    }
    grp.createDataSet("channel_names", names);
    grp.createDataSet("channel_units", units);

    {
        std::vector<std::size_t> dims = {N, C, static_cast<std::size_t>(n_lags_)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>("acf", space);
        ds.write_raw(acf_.data());
        ds.createAttribute("units", std::string("dimensionless"));
        ds.createAttribute("layout", std::string("(atom, channel, lag)"));
    }
    {
        std::vector<std::size_t> dims = {N, C, static_cast<std::size_t>(n_freq_)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>("power_spectrum", space);
        ds.write_raw(spectrum_.data());
        ds.createAttribute("units", std::string("channel_units^2 * ps"));
        ds.createAttribute("layout", std::string("(atom, channel, frequency)"));
    }

    auto write_reduction = [&](const std::string& name,
                               const std::vector<double>& v,
                               const std::string& note) {
        std::vector<std::size_t> dims = {N, C};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(v.data());
        ds.createAttribute("layout", std::string("(atom, channel)"));
        ds.createAttribute("note", note);
    };
    write_reduction("decay_time_ps", decay_time_ps_,
        "dt * sum of acf from lag 0 to the first zero crossing; full "
        "window (a lower bound) when acf stays positive; NaN for "
        "constant signals. Units ps.");
    write_reduction("peak_freq_per_ps", peak_freq_per_ps_,
        "frequency of the largest power_spectrum bin excluding DC. "
        "Units 1/ps.");
    write_reduction("spectral_centroid_per_ps", spectral_centroid_per_ps_,
        "sum(f * S) / sum(S) over all bins. Units 1/ps.");

    std::vector<std::uint64_t> lag_frames(n_lags_);
    std::vector<double> lag_times(n_lags_);
    for (std::size_t k = 0; k < n_lags_; ++k) {
        lag_frames[k] = k;
        lag_times[k] = static_cast<double>(k) * sample_interval_ps_;
    }
    grp.createDataSet("lag_frames", lag_frames);
    grp.createDataSet("lag_times_ps", lag_times);
    grp.createDataSet("frequencies_per_ps",
                      SpectrumFrequenciesPerPs(n_lags_, sample_interval_ps_, n_freq_));
}

}  // namespace nmr
