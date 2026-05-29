#include "KernelCoherenceTrajectoryResult.h"

#include "ApbsFieldResult.h"
#include "BiotSavartResult.h"
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

#include <cmath>
#include <typeinfo>

namespace nmr {

namespace {

using Channel = KernelCoherenceTrajectoryResult::Channel;

// Channel value/name/units, cloned from KernelDynamicsTrajectoryResult so
// the two TRs stay independent (PATTERNS.md 17). Same 13 kernel scalars.
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
KernelCoherenceTrajectoryResult::Dependencies() const {
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


std::unique_ptr<KernelCoherenceTrajectoryResult>
KernelCoherenceTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<KernelCoherenceTrajectoryResult>();
    r->n_atoms_ = tp.AtomCount();
    r->sum_x_.assign(r->n_atoms_ * N_CHANNELS, 0.0);
    r->sum_xx_.assign(r->n_atoms_ * N_CHANNELS, 0.0);
    r->sum_xy_.assign(r->n_atoms_ * N_PAIRS, 0.0);
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Per atom: read the 13 channel scalars once, then accumulate sum, sum of
// squares, and the upper-triangle pairwise products. The pair index walks
// (a <= b) in the same row-major order Finalize reads back.

void KernelCoherenceTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)frame_idx; (void)time_ps;
    const std::size_t N = tp.AtomCount();
    const std::size_t C = N_CHANNELS;
    double v[N_CHANNELS];
    for (std::size_t i = 0; i < N; ++i) {
        const ConformationAtom& a = conf.AtomAt(i);
        for (std::size_t c = 0; c < C; ++c) {
            v[c] = ChannelValue(static_cast<Channel>(c), a);
        }
        double* sx  = &sum_x_[i * C];
        double* sxx = &sum_xx_[i * C];
        double* sxy = &sum_xy_[i * N_PAIRS];
        std::size_t p = 0;
        for (std::size_t a_ch = 0; a_ch < C; ++a_ch) {
            sx[a_ch]  += v[a_ch];
            sxx[a_ch] += v[a_ch] * v[a_ch];
            for (std::size_t b_ch = a_ch; b_ch < C; ++b_ch) {
                sxy[p++] += v[a_ch] * v[b_ch];
            }
        }
    }
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// r_xy = cov_xy / sqrt(var_x var_y), with cov_xy = <xy> - <x><y> and
// var_x = <x^2> - <x>^2 (the maximum-likelihood / biased moments, T in
// the denominator throughout so the ratio is unbiased). A channel with
// var ~ 0 (constant over the run) has undefined correlation: its whole
// row and column are NaN; the diagonal is 1.0 for non-constant channels.

void KernelCoherenceTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                               Trajectory& traj) {
    (void)traj;
    if (finalized_) return;
    n_atoms_ = tp.AtomCount();
    const std::size_t N = n_atoms_;
    const std::size_t C = N_CHANNELS;
    const double nan_val = std::nan("");
    const double T = static_cast<double>(n_frames_);

    correlation_matrix_.assign(N * C * C, nan_val);
    if (n_frames_ < 2) { finalized_ = true; return; }

    std::vector<double> mean(C), var(C);
    for (std::size_t i = 0; i < N; ++i) {
        const double* sx  = &sum_x_[i * C];
        const double* sxx = &sum_xx_[i * C];
        const double* sxy = &sum_xy_[i * N_PAIRS];
        for (std::size_t c = 0; c < C; ++c) {
            mean[c] = sx[c] / T;
            var[c]  = sxx[c] / T - mean[c] * mean[c];
        }
        double* mat = &correlation_matrix_[i * C * C];
        std::size_t p = 0;
        for (std::size_t a = 0; a < C; ++a) {
            for (std::size_t b = a; b < C; ++b) {
                const double cov = sxy[p++] / T - mean[a] * mean[b];
                double r = nan_val;
                if (var[a] > 1e-15 && var[b] > 1e-15) {
                    r = cov / std::sqrt(var[a] * var[b]);
                }
                mat[a * C + b] = r;
                mat[b * C + a] = r;  // symmetric
            }
        }
    }

    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "KernelCoherenceTrajectoryResult::Finalize",
        "zero-lag kernel correlation matrix across " + std::to_string(N) +
        " atoms x " + std::to_string(C) + " channels, " +
        std::to_string(n_frames_) + " frames");
}


// ── WriteH5Group ─────────────────────────────────────────────────

void KernelCoherenceTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    if (!finalized_) {
        OperationLog::Warn("KernelCoherenceTrajectoryResult::WriteH5Group",
                           "Finalize not called; nothing to write");
        return;
    }
    const std::size_t N = n_atoms_;
    const std::size_t C = N_CHANNELS;

    auto grp = file.createGroup("/trajectory/kernel_coherence");
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_channels",  C);
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("statistic",   std::string("pearson_zero_lag"));
    grp.createAttribute("lagged_cross_correlation", std::string("deferred"));
    grp.createAttribute("constant_channel_policy", std::string(
        "a channel constant over the run (var ~ 0) has undefined "
        "correlation: its row and column are NaN. Diagonal is 1.0 for "
        "non-constant channels. Use isfinite() to mask."));

    std::vector<std::string> names(C), units(C);
    for (std::size_t c = 0; c < C; ++c) {
        names[c] = ChannelName(static_cast<Channel>(c));
        units[c] = ChannelUnits(static_cast<Channel>(c));
    }
    grp.createDataSet("channel_names", names);
    grp.createDataSet("channel_units", units);

    std::vector<std::size_t> dims = {N, C, C};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("correlation_matrix", space);
    ds.write_raw(correlation_matrix_.data());
    ds.createAttribute("units", std::string("dimensionless"));
    ds.createAttribute("layout", std::string("(atom, channel, channel)"));
}

}  // namespace nmr
