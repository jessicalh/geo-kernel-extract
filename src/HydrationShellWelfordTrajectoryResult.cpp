#include "HydrationShellWelfordTrajectoryResult.h"
#include "HydrationShellResult.h"
#include "ConformationAtom.h"
#include "ProteinConformation.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <functional>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index>
HydrationShellWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(HydrationShellResult)) };
}


std::unique_ptr<HydrationShellWelfordTrajectoryResult>
HydrationShellWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<HydrationShellWelfordTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    r->prev_half_shell_asymmetry_.assign(N, 0.0);
    r->prev_mean_water_dipole_cos_.assign(N, 0.0);
    r->prev_nearest_ion_distance_.assign(N, 0.0);
    r->prev_nearest_ion_charge_.assign(N, 0.0);
    r->prev_valid_.assign(N, false);
    r->prev_ion_finite_.assign(N, false);
    r->prev_time_.assign(N, 0.0);
    return r;
}


void HydrationShellWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<HydrationShellResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);
    if (!source_attached) {
        for (std::size_t i = 0; i < N; ++i) {
            prev_valid_[i] = false;
            prev_ion_finite_[i] = false;
        }
        ++n_frames_;
        return;
    }

    for (std::size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        HydrationShellWelfordState& w = tp.MutableAtomAt(i).hydration_shell_welford;
        const std::size_t n_new = w.n_frames + 1;

        const double hsa = a.half_shell_asymmetry;
        const double dpc = a.mean_water_dipole_cos;
        const double nid = a.nearest_ion_distance;   // may be +inf
        const double nic = a.nearest_ion_charge;
        const bool   ion_finite = std::isfinite(nid);
        const double ion_indicator = ion_finite ? 1.0 : 0.0;

        // Unconditional channels (3 + presence indicator)
        WelfordUpdate(w.half_shell_asymmetry,  hsa,           n_new, frame_idx);
        WelfordUpdate(w.mean_water_dipole_cos, dpc,           n_new, frame_idx);
        WelfordUpdate(w.nearest_ion_charge,    nic,           n_new, frame_idx);
        WelfordUpdate(w.ion_present_fraction,  ion_indicator, n_new, frame_idx);

        // Conditional: nearest_ion_distance only when finite. Per-atom
        // counter n_ion_present is the Finalize divisor.
        if (ion_finite) {
            const std::size_t nip_new = w.n_ion_present + 1;
            WelfordUpdate(w.nearest_ion_distance, nid, nip_new, frame_idx);
            w.n_ion_present = nip_new;
        }

        w.n_frames = n_new;

        if (prev_valid_[i]) {
            const std::size_t dn_new = w.delta_n + 1;
            // Unconditional deltas (3 channels)
            auto upd = [&](double curr, double prev,
                           WelfordMoments& d_, WelfordMoments& ad,
                           WelfordMoments& sd) {
                const double delta = curr - prev;
                WelfordUpdate(d_, delta,           dn_new, frame_idx);
                WelfordUpdate(ad, std::abs(delta), dn_new, frame_idx);
                WelfordUpdate(sd, delta * delta,   dn_new, frame_idx);
            };
            upd(hsa, prev_half_shell_asymmetry_[i],
                w.half_shell_asymmetry_delta,
                w.half_shell_asymmetry_abs_delta,
                w.half_shell_asymmetry_delta_squared);
            upd(dpc, prev_mean_water_dipole_cos_[i],
                w.mean_water_dipole_cos_delta,
                w.mean_water_dipole_cos_abs_delta,
                w.mean_water_dipole_cos_delta_squared);
            upd(nic, prev_nearest_ion_charge_[i],
                w.nearest_ion_charge_delta,
                w.nearest_ion_charge_abs_delta,
                w.nearest_ion_charge_delta_squared);
            w.delta_n = dn_new;

            // Conditional delta on nearest_ion_distance: BOTH frames must
            // have finite distance. Per-atom counter n_ion_delta.
            if (ion_finite && prev_ion_finite_[i]) {
                const std::size_t nid_new = w.n_ion_delta + 1;
                const double delta = nid - prev_nearest_ion_distance_[i];
                WelfordUpdate(w.nearest_ion_distance_delta,
                              delta, nid_new, frame_idx);
                WelfordUpdate(w.nearest_ion_distance_abs_delta,
                              std::abs(delta), nid_new, frame_idx);
                WelfordUpdate(w.nearest_ion_distance_delta_squared,
                              delta * delta, nid_new, frame_idx);
                w.n_ion_delta = nid_new;
            }

            constexpr double MIN_DT_PS = 1e-12;
            const double dt = time_ps - prev_time_[i];
            if (std::abs(dt) > MIN_DT_PS) {
                const std::size_t dxn = w.dxdt_n + 1;
                WelfordUpdate(w.half_shell_asymmetry_dxdt,
                              (hsa - prev_half_shell_asymmetry_[i]) / dt, dxn, frame_idx);
                WelfordUpdate(w.mean_water_dipole_cos_dxdt,
                              (dpc - prev_mean_water_dipole_cos_[i]) / dt, dxn, frame_idx);
                WelfordUpdate(w.nearest_ion_charge_dxdt,
                              (nic - prev_nearest_ion_charge_[i]) / dt, dxn, frame_idx);
                w.dxdt_n = dxn;

                // Conditional dxdt on nearest_ion_distance: same finite-pair
                // gate AND dt gate. Per-atom counter n_ion_dxdt.
                if (ion_finite && prev_ion_finite_[i]) {
                    const std::size_t nid_dxn = w.n_ion_dxdt + 1;
                    WelfordUpdate(w.nearest_ion_distance_dxdt,
                                  (nid - prev_nearest_ion_distance_[i]) / dt,
                                  nid_dxn, frame_idx);
                    w.n_ion_dxdt = nid_dxn;
                }
            }
        }
        prev_half_shell_asymmetry_[i]  = hsa;
        prev_mean_water_dipole_cos_[i] = dpc;
        // Only update prev_nearest_ion_distance_ when finite — otherwise
        // the cached value stays at the last-known-finite, ready for
        // the next attached+finite pair. prev_ion_finite_ records the
        // CURRENT frame's finiteness for the next frame's pair test.
        if (ion_finite) prev_nearest_ion_distance_[i] = nid;
        prev_nearest_ion_charge_[i]    = nic;
        prev_ion_finite_[i]            = ion_finite;
        prev_time_[i]                  = time_ps;
        prev_valid_[i]                 = true;
    }
    ++n_frames_;
}


void HydrationShellWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        HydrationShellWelfordState& w = tp.MutableAtomAt(i).hydration_shell_welford;

        WelfordFinalize(w.half_shell_asymmetry,  w.n_frames);
        WelfordFinalize(w.mean_water_dipole_cos, w.n_frames);
        // R6: nearest_ion_distance uses its own per-atom counter — the
        // Welford only accumulated frames with finite ion-in-cutoff. Atoms
        // where n_ion_present == 0 finalize to NaN, consistent with
        // "no samples observed" semantics.
        WelfordFinalize(w.nearest_ion_distance,  w.n_ion_present);
        WelfordFinalize(w.nearest_ion_charge,    w.n_frames);
        WelfordFinalize(w.ion_present_fraction,  w.n_frames);

        auto fin = [&](WelfordMoments& d_, WelfordMoments& ad,
                       WelfordMoments& sd, WelfordMoments& dx,
                       double& rms,
                       std::size_t n_d, std::size_t n_dx) {
            WelfordFinalize(d_, n_d);
            WelfordFinalize(ad, n_d);
            WelfordFinalize(sd, n_d);
            WelfordFinalize(dx, n_dx);
            rms = (n_d == 0) ? std::nan("") : std::sqrt(sd.mean);
        };
        fin(w.half_shell_asymmetry_delta, w.half_shell_asymmetry_abs_delta,
            w.half_shell_asymmetry_delta_squared, w.half_shell_asymmetry_dxdt,
            w.half_shell_asymmetry_rms_delta,
            w.delta_n, w.dxdt_n);
        fin(w.mean_water_dipole_cos_delta, w.mean_water_dipole_cos_abs_delta,
            w.mean_water_dipole_cos_delta_squared, w.mean_water_dipole_cos_dxdt,
            w.mean_water_dipole_cos_rms_delta,
            w.delta_n, w.dxdt_n);
        // R6: nearest_ion_distance delta variants used their own per-atom
        // pair counter (BOTH prev and curr finite). dxdt variant gated
        // additionally on dt > MIN_DT_PS.
        fin(w.nearest_ion_distance_delta, w.nearest_ion_distance_abs_delta,
            w.nearest_ion_distance_delta_squared, w.nearest_ion_distance_dxdt,
            w.nearest_ion_distance_rms_delta,
            w.n_ion_delta, w.n_ion_dxdt);
        fin(w.nearest_ion_charge_delta, w.nearest_ion_charge_abs_delta,
            w.nearest_ion_charge_delta_squared, w.nearest_ion_charge_dxdt,
            w.nearest_ion_charge_rms_delta,
            w.delta_n, w.dxdt_n);
    }

    // mean_dt_ps + frame_index_range from source-attached subset only.
    const auto& times = traj.FrameTimes();
    std::vector<double> attached_times;
    attached_times.reserve(source_attached_per_frame_.size());
    for (std::size_t f = 0;
         f < source_attached_per_frame_.size() && f < times.size(); ++f) {
        if (source_attached_per_frame_[f]) attached_times.push_back(times[f]);
    }
    if (attached_times.size() >= 2) {
        mean_dt_ps_ = (attached_times.back() - attached_times.front()) /
                      static_cast<double>(attached_times.size() - 1);
    }
    const auto& fidx = traj.FrameIndices();
    std::vector<std::size_t> attached_fidx;
    attached_fidx.reserve(source_attached_per_frame_.size());
    for (std::size_t f = 0;
         f < source_attached_per_frame_.size() && f < fidx.size(); ++f) {
        if (source_attached_per_frame_[f]) attached_fidx.push_back(fidx[f]);
    }
    if (!attached_fidx.empty()) {
        frame_index_range_ = {attached_fidx.front(), attached_fidx.back()};
    }

    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "HydrationShellWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(N) + " atoms");
}


void HydrationShellWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();

    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "HydrationShellWelfordTrajectoryResult::WriteH5Group",
            "HydrationShellResult attached in 0/" +
            std::to_string(source_attached_per_frame_.size()) +
            " frames; skipping /trajectory/hydration_shell_welford/.");
        return;
    }

    auto grp = file.createGroup("/trajectory/hydration_shell_welford");

    grp.createAttribute("result_name",           Name());
    grp.createAttribute("n_frames",              n_frames_);
    grp.createAttribute("source_attached_count", source_attached_count);
    grp.createAttribute("finalized",             finalized_);
    grp.createAttribute("ddof",                  static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",            mean_dt_ps_);
    grp.createAttribute("frame_index_range",     frame_index_range_);
    grp.createAttribute("reference_frame",       std::string("COM"));
    // R6 codex 2026-05-18: changed from naive inf-poisoned Welford to
    // conditional Welford + presence-of-ion indicator. New semantics:
    grp.createAttribute("nearest_ion_distance_welford_policy",
        std::string("finite_only: only frames with finite "
                    "nearest_ion_distance contribute to Welford. "
                    "Per-atom counter `n_ion_present_per_atom` is the "
                    "divisor. Atoms that never see an ion in cutoff "
                    "finalize to NaN (n=0). Delta variants gated on "
                    "BOTH prev and curr frames being finite "
                    "(n_ion_delta_per_atom); dxdt likewise plus dt>0 "
                    "(n_ion_dxdt_per_atom). For probability of any ion "
                    "in cutoff, use `ion_present_fraction_mean ∈ [0,1]`."));

    auto get_w = [&tp](std::size_t i) -> const HydrationShellWelfordState& {
        return tp.AtomAt(i).hydration_shell_welford;
    };

    auto emit_1d = [&](const std::string& prefix,
                       const std::string& base_units,
                       const std::string& m2_units,
                       std::function<const WelfordMoments&(std::size_t)> get) {
        std::vector<double> mean(N), m2(N), std_(N), min_(N), max_(N);
        std::vector<std::size_t> minf(N), maxf(N);
        for (std::size_t i = 0; i < N; ++i) {
            const WelfordMoments& w = get(i);
            mean[i] = w.mean; m2[i] = w.m2; std_[i] = w.std;
            min_[i] = w.min; max_[i] = w.max;
            minf[i] = w.min_frame; maxf[i] = w.max_frame;
        }
        auto put = [&](const std::string& n_, const std::vector<double>& v,
                       const std::string& u) {
            auto ds = grp.createDataSet(n_, v); ds.createAttribute("units", u);
        };
        put(prefix + "_mean", mean, base_units);
        put(prefix + "_m2",   m2,   m2_units);
        put(prefix + "_std",  std_, base_units);
        put(prefix + "_min",  min_, base_units);
        put(prefix + "_max",  max_, base_units);
        auto dminf = grp.createDataSet(prefix + "_min_frame", minf);
        dminf.createAttribute("units", std::string("frame_index"));
        auto dmaxf = grp.createDataSet(prefix + "_max_frame", maxf);
        dmaxf.createAttribute("units", std::string("frame_index"));
    };

    struct Bundle {
        const char* base;
        const std::string& base_units;
        const std::string& m2_units;
        const std::string& sq_units;
        const std::string& sq_m2_units;
        const std::string& rate_units;
        const std::string& rate_m2_units;
        std::function<const WelfordMoments&(const HydrationShellWelfordState&)> val;
        std::function<const WelfordMoments&(const HydrationShellWelfordState&)> d;
        std::function<const WelfordMoments&(const HydrationShellWelfordState&)> ad;
        std::function<const WelfordMoments&(const HydrationShellWelfordState&)> sd;
        std::function<const WelfordMoments&(const HydrationShellWelfordState&)> dx;
        std::function<double(const HydrationShellWelfordState&)> rms;
    };
    const std::string kFrac    = "fraction";
    const std::string kFracSq  = "fraction^2";
    const std::string kFrac4   = "fraction^4";
    const std::string kFracRate = "fraction/ps";
    const std::string kFracRateSq = "fraction^2/ps^2";
    const std::string kCos     = "cos_angle";
    const std::string kCosSq   = "cos_angle^2";
    const std::string kCos4    = "cos_angle^4";
    const std::string kCosRate = "cos_angle/ps";
    const std::string kCosRateSq = "cos_angle^2/ps^2";
    const std::string kAng     = "Angstrom";
    const std::string kAngSq   = "Angstrom^2";
    const std::string kAng4    = "Angstrom^4";
    const std::string kAngRate = "Angstrom/ps";
    const std::string kAngRateSq = "Angstrom^2/ps^2";
    const std::string kCharge  = "e";
    const std::string kChargeSq= "e^2";
    const std::string kCharge4 = "e^4";
    const std::string kChargeRate = "e/ps";
    const std::string kChargeRateSq = "e^2/ps^2";

    const std::vector<Bundle> bundles = {
        { "half_shell_asymmetry", kFrac, kFracSq, kFracSq, kFrac4, kFracRate, kFracRateSq,
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_abs_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_delta_squared; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_dxdt; },
          [](const HydrationShellWelfordState& w) -> double { return w.half_shell_asymmetry_rms_delta; } },
        { "mean_water_dipole_cos", kCos, kCosSq, kCosSq, kCos4, kCosRate, kCosRateSq,
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.mean_water_dipole_cos; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.mean_water_dipole_cos_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.mean_water_dipole_cos_abs_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.mean_water_dipole_cos_delta_squared; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.mean_water_dipole_cos_dxdt; },
          [](const HydrationShellWelfordState& w) -> double { return w.mean_water_dipole_cos_rms_delta; } },
        { "nearest_ion_distance", kAng, kAngSq, kAngSq, kAng4, kAngRate, kAngRateSq,
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_distance; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_distance_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_distance_abs_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_distance_delta_squared; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_distance_dxdt; },
          [](const HydrationShellWelfordState& w) -> double { return w.nearest_ion_distance_rms_delta; } },
        { "nearest_ion_charge", kCharge, kChargeSq, kChargeSq, kCharge4, kChargeRate, kChargeRateSq,
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_charge; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_charge_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_charge_abs_delta; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_charge_delta_squared; },
          [](const HydrationShellWelfordState& w) -> const WelfordMoments& { return w.nearest_ion_charge_dxdt; },
          [](const HydrationShellWelfordState& w) -> double { return w.nearest_ion_charge_rms_delta; } },
    };

    for (const auto& b : bundles) {
        emit_1d(b.base,                                b.base_units, b.m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.val(get_w(i)); });
        emit_1d(std::string(b.base) + "_delta",        b.base_units, b.m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.d(get_w(i)); });
        emit_1d(std::string(b.base) + "_abs_delta",    b.base_units, b.m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.ad(get_w(i)); });
        emit_1d(std::string(b.base) + "_delta_squared", b.sq_units,  b.sq_m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.sd(get_w(i)); });
        emit_1d(std::string(b.base) + "_dxdt",         b.rate_units, b.rate_m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.dx(get_w(i)); });

        std::vector<double> rms(N);
        for (std::size_t i = 0; i < N; ++i) rms[i] = b.rms(get_w(i));
        auto rms_ds = grp.createDataSet(std::string(b.base) + "_rms_delta", rms);
        rms_ds.createAttribute("units", b.base_units);
    }

    // R6: ion_present_fraction — indicator-Welford on whether a finite
    // nearest_ion_distance is observed this frame. Mean ∈ [0, 1] is
    // interpretable as Pr(ion in cutoff per attached frame). Emitted
    // outside the bundle loop because it has no delta variants.
    {
        const std::string kProb = "probability";
        const std::string kProbSq = "probability^2";
        emit_1d("ion_present_fraction", kProb, kProbSq,
                [get_w](std::size_t i) -> const WelfordMoments& {
                    return get_w(i).ion_present_fraction;
                });
    }

    std::vector<std::size_t> n_frames(N), delta_n(N), dxdt_n(N);
    std::vector<std::size_t> n_ion_present(N), n_ion_delta(N), n_ion_dxdt(N);
    for (std::size_t i = 0; i < N; ++i) {
        const HydrationShellWelfordState& w = get_w(i);
        n_frames[i]      = w.n_frames;
        delta_n[i]       = w.delta_n;
        dxdt_n[i]        = w.dxdt_n;
        n_ion_present[i] = w.n_ion_present;
        n_ion_delta[i]   = w.n_ion_delta;
        n_ion_dxdt[i]    = w.n_ion_dxdt;
    }
    grp.createDataSet("n_frames_per_atom", n_frames)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("delta_n_per_atom",  delta_n)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("dxdt_n_per_atom",   dxdt_n)
       .createAttribute("units", std::string("frame_count"));
    // R6: per-atom counters specific to the conditional nearest_ion_distance
    // Welford. n_ion_present is the base-Welford divisor; n_ion_delta and
    // n_ion_dxdt are the divisors for the delta and dxdt variants of
    // nearest_ion_distance ONLY (other channels still use delta_n/dxdt_n).
    grp.createDataSet("n_ion_present_per_atom", n_ion_present)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("n_ion_delta_per_atom",   n_ion_delta)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("n_ion_dxdt_per_atom",    n_ion_dxdt)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_)
       .createAttribute("units", std::string("dimensionless"));
}

}  // namespace nmr
