#include "HydrationGeometryWelfordTrajectoryResult.h"
#include "HydrationGeometryResult.h"
#include "ConformationAtom.h"
#include "ProteinConformation.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"
#include "Types.h"
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
HydrationGeometryWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(HydrationGeometryResult)) };
}


std::unique_ptr<HydrationGeometryWelfordTrajectoryResult>
HydrationGeometryWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<HydrationGeometryWelfordTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    r->prev_half_shell_asymmetry_.assign(N, 0.0);
    r->prev_dipole_alignment_.assign(N, 0.0);
    r->prev_dipole_coherence_.assign(N, 0.0);
    r->prev_shell_count_.assign(N, 0.0);
    r->prev_valid_.assign(N, false);
    r->prev_time_.assign(N, 0.0);
    return r;
}


void HydrationGeometryWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<HydrationGeometryResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);
    if (!source_attached) {
        for (std::size_t i = 0; i < N; ++i) prev_valid_[i] = false;
        ++n_frames_;
        return;
    }

    for (std::size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        HydrationGeometryWelfordState& w =
            tp.MutableAtomAt(i).hydration_geometry_welford;
        const std::size_t n_new = w.n_frames + 1;

        // Vec3 channels: per-component + magnitude (dipole only)
        const Vec3& d = a.water_dipole_vector;
        const Vec3& n = a.water_surface_normal;
        const double d_mag = d.norm();
        WelfordUpdate(w.dipole_vector[0], d.x(), n_new, frame_idx);
        WelfordUpdate(w.dipole_vector[1], d.y(), n_new, frame_idx);
        WelfordUpdate(w.dipole_vector[2], d.z(), n_new, frame_idx);
        WelfordUpdate(w.dipole_magnitude, d_mag, n_new, frame_idx);
        WelfordUpdate(w.surface_normal[0], n.x(), n_new, frame_idx);
        WelfordUpdate(w.surface_normal[1], n.y(), n_new, frame_idx);
        WelfordUpdate(w.surface_normal[2], n.z(), n_new, frame_idx);

        // Scalar channels
        const double hsa  = a.sasa_half_shell_asymmetry;
        const double dal  = a.sasa_dipole_alignment;
        const double dco  = a.sasa_dipole_coherence;
        const double sct  = static_cast<double>(a.sasa_first_shell_count);
        WelfordUpdate(w.half_shell_asymmetry, hsa, n_new, frame_idx);
        WelfordUpdate(w.dipole_alignment,     dal, n_new, frame_idx);
        WelfordUpdate(w.dipole_coherence,     dco, n_new, frame_idx);
        WelfordUpdate(w.shell_count,          sct, n_new, frame_idx);

        w.n_frames = n_new;

        // Delta variants on the 4 primary scalars
        if (prev_valid_[i]) {
            const std::size_t dn_new = w.delta_n + 1;
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
            upd(dal, prev_dipole_alignment_[i],
                w.dipole_alignment_delta,
                w.dipole_alignment_abs_delta,
                w.dipole_alignment_delta_squared);
            upd(dco, prev_dipole_coherence_[i],
                w.dipole_coherence_delta,
                w.dipole_coherence_abs_delta,
                w.dipole_coherence_delta_squared);
            upd(sct, prev_shell_count_[i],
                w.shell_count_delta,
                w.shell_count_abs_delta,
                w.shell_count_delta_squared);
            w.delta_n = dn_new;

            constexpr double MIN_DT_PS = 1e-12;
            const double dt = time_ps - prev_time_[i];
            if (std::abs(dt) > MIN_DT_PS) {
                const std::size_t dxn = w.dxdt_n + 1;
                WelfordUpdate(w.half_shell_asymmetry_dxdt,
                              (hsa - prev_half_shell_asymmetry_[i]) / dt, dxn, frame_idx);
                WelfordUpdate(w.dipole_alignment_dxdt,
                              (dal - prev_dipole_alignment_[i]) / dt, dxn, frame_idx);
                WelfordUpdate(w.dipole_coherence_dxdt,
                              (dco - prev_dipole_coherence_[i]) / dt, dxn, frame_idx);
                WelfordUpdate(w.shell_count_dxdt,
                              (sct - prev_shell_count_[i]) / dt, dxn, frame_idx);
                w.dxdt_n = dxn;
            }
        }
        prev_half_shell_asymmetry_[i] = hsa;
        prev_dipole_alignment_[i]     = dal;
        prev_dipole_coherence_[i]     = dco;
        prev_shell_count_[i]          = sct;
        prev_time_[i]                 = time_ps;
        prev_valid_[i]                = true;
    }
    ++n_frames_;
}


void HydrationGeometryWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                       Trajectory& traj) {
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        HydrationGeometryWelfordState& w =
            tp.MutableAtomAt(i).hydration_geometry_welford;

        for (std::size_t k = 0; k < 3; ++k) WelfordFinalize(w.dipole_vector[k],  w.n_frames);
        WelfordFinalize(w.dipole_magnitude, w.n_frames);
        for (std::size_t k = 0; k < 3; ++k) WelfordFinalize(w.surface_normal[k], w.n_frames);
        WelfordFinalize(w.half_shell_asymmetry, w.n_frames);
        WelfordFinalize(w.dipole_alignment,     w.n_frames);
        WelfordFinalize(w.dipole_coherence,     w.n_frames);
        WelfordFinalize(w.shell_count,          w.n_frames);

        auto fin = [&](WelfordMoments& d_, WelfordMoments& ad,
                       WelfordMoments& sd, WelfordMoments& dx,
                       double& rms) {
            WelfordFinalize(d_, w.delta_n);
            WelfordFinalize(ad, w.delta_n);
            WelfordFinalize(sd, w.delta_n);
            WelfordFinalize(dx, w.dxdt_n);
            rms = (w.delta_n == 0) ? std::nan("") : std::sqrt(sd.mean);
        };
        fin(w.half_shell_asymmetry_delta, w.half_shell_asymmetry_abs_delta,
            w.half_shell_asymmetry_delta_squared, w.half_shell_asymmetry_dxdt,
            w.half_shell_asymmetry_rms_delta);
        fin(w.dipole_alignment_delta, w.dipole_alignment_abs_delta,
            w.dipole_alignment_delta_squared, w.dipole_alignment_dxdt,
            w.dipole_alignment_rms_delta);
        fin(w.dipole_coherence_delta, w.dipole_coherence_abs_delta,
            w.dipole_coherence_delta_squared, w.dipole_coherence_dxdt,
            w.dipole_coherence_rms_delta);
        fin(w.shell_count_delta, w.shell_count_abs_delta,
            w.shell_count_delta_squared, w.shell_count_dxdt,
            w.shell_count_rms_delta);
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
        "HydrationGeometryWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(N) + " atoms");
}


void HydrationGeometryWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();

    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "HydrationGeometryWelfordTrajectoryResult::WriteH5Group",
            "HydrationGeometryResult attached in 0/" +
            std::to_string(source_attached_per_frame_.size()) +
            " frames; skipping /trajectory/hydration_geometry_welford/.");
        return;
    }

    auto grp = file.createGroup("/trajectory/hydration_geometry_welford");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_frames",               n_frames_);
    grp.createAttribute("source_attached_count",  source_attached_count);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("ddof",                   static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",             mean_dt_ps_);
    grp.createAttribute("frame_index_range",      frame_index_range_);
    grp.createAttribute("irrep_layout_dipole",    std::string("v_x,v_y,v_z"));
    grp.createAttribute("irrep_layout_normal",    std::string("v_x,v_y,v_z"));
    grp.createAttribute("reference_frame",        std::string("SASA_normal"));
    grp.createAttribute("dipole_alignment_zero_sentinel",
        std::string("0.0 when |dipole_sum| < NEAR_ZERO_NORM "
                    "(HydrationGeometryResult.cpp:106-108) — bimodal: "
                    "real zero alignment vs zero-magnitude denominator"));
    grp.createAttribute("polarisation_signal_channels",
        std::string("dipole_alignment,dipole_coherence,half_shell_asymmetry"));
    grp.createAttribute("units",                  std::string("mixed_see_per_dataset"));

    auto get_w = [&tp](std::size_t i) -> const HydrationGeometryWelfordState& {
        return tp.AtomAt(i).hydration_geometry_welford;
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

    // Vec3 channels: per-component
    static const std::array<const char*, 3> kXyz = {"x", "y", "z"};
    // Dipole units: e·Å (raw H_charge·displacement sum from Water::Dipole(),
    // SolventEnvironment.h:32-38). Convert to Debye via 1 e·Å = 4.80320 D
    // if needed downstream. R5 codex 2026-05-18.
    const std::string kDipole = "e_Angstrom";
    const std::string kDipoleSq = "e^2_Angstrom^2";
    const std::string kUnit = "dimensionless";
    for (std::size_t k = 0; k < 3; ++k) {
        emit_1d(std::string("dipole_vector_") + kXyz[k], kDipole, kDipoleSq,
                [get_w, k](std::size_t i) -> const WelfordMoments& {
                    return get_w(i).dipole_vector[k];
                });
        emit_1d(std::string("surface_normal_") + kXyz[k], kUnit, kUnit,
                [get_w, k](std::size_t i) -> const WelfordMoments& {
                    return get_w(i).surface_normal[k];
                });
    }
    emit_1d("dipole_magnitude", kDipole, kDipoleSq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).dipole_magnitude;
            });

    // Scalar polarisation channels
    const std::string kFrac = "fraction";
    const std::string kFracSq = "fraction^2";
    const std::string kCos = "cos_angle";
    const std::string kCosSq = "cos_angle^2";
    // dipole_coherence: source formula |Σ d_i|/n_shell is e·Å, NOT a
    // [0,1] dimensionless order parameter. Same fix as the TS sibling,
    // R6 codex 2026-05-18.
    const std::string kOrd = "e_Angstrom";
    const std::string kOrdSq = "e^2_Angstrom^2";
    emit_1d("half_shell_asymmetry", kFrac, kFracSq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).half_shell_asymmetry;
            });
    emit_1d("dipole_alignment", kCos, kCosSq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).dipole_alignment;
            });
    emit_1d("dipole_coherence", kOrd, kOrdSq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).dipole_coherence;
            });
    emit_1d("shell_count", kUnit, kUnit,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).shell_count;
            });

    // Delta variants on the 4 primary scalars
    struct DeltaBundle {
        const char* base;
        const std::string& base_units;
        const std::string& m2_units;
        const std::string& sq_units;
        const std::string& sq_m2_units;
        const std::string& rate_units;
        const std::string& rate_m2_units;
        std::function<const WelfordMoments&(const HydrationGeometryWelfordState&)> d;
        std::function<const WelfordMoments&(const HydrationGeometryWelfordState&)> ad;
        std::function<const WelfordMoments&(const HydrationGeometryWelfordState&)> sd;
        std::function<const WelfordMoments&(const HydrationGeometryWelfordState&)> dx;
        std::function<double(const HydrationGeometryWelfordState&)> rms;
    };
    const std::string kFrac4   = "fraction^4";
    const std::string kCos4    = "cos_angle^4";
    const std::string kOrd4    = "e^4_Angstrom^4";
    const std::string kFracRate = "fraction/ps";
    const std::string kFracRateSq = "fraction^2/ps^2";
    const std::string kCosRate = "cos_angle/ps";
    const std::string kCosRateSq = "cos_angle^2/ps^2";
    const std::string kOrdRate = "e_Angstrom/ps";
    const std::string kOrdRateSq = "e^2_Angstrom^2/ps^2";
    const std::string kCountRate = "count/ps";
    const std::string kCountRateSq = "count^2/ps^2";
    const std::vector<DeltaBundle> bundles = {
        { "half_shell_asymmetry", kFrac, kFracSq, kFracSq, kFrac4, kFracRate, kFracRateSq,
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_abs_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_delta_squared; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.half_shell_asymmetry_dxdt; },
          [](const HydrationGeometryWelfordState& w) -> double { return w.half_shell_asymmetry_rms_delta; } },
        { "dipole_alignment", kCos, kCosSq, kCosSq, kCos4, kCosRate, kCosRateSq,
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_alignment_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_alignment_abs_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_alignment_delta_squared; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_alignment_dxdt; },
          [](const HydrationGeometryWelfordState& w) -> double { return w.dipole_alignment_rms_delta; } },
        { "dipole_coherence", kOrd, kOrdSq, kOrdSq, kOrd4, kOrdRate, kOrdRateSq,
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_coherence_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_coherence_abs_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_coherence_delta_squared; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.dipole_coherence_dxdt; },
          [](const HydrationGeometryWelfordState& w) -> double { return w.dipole_coherence_rms_delta; } },
        { "shell_count", kUnit, kUnit, kUnit, kUnit, kCountRate, kCountRateSq,
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.shell_count_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.shell_count_abs_delta; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.shell_count_delta_squared; },
          [](const HydrationGeometryWelfordState& w) -> const WelfordMoments& { return w.shell_count_dxdt; },
          [](const HydrationGeometryWelfordState& w) -> double { return w.shell_count_rms_delta; } },
    };
    for (const auto& b : bundles) {
        emit_1d(std::string(b.base) + "_delta",         b.base_units, b.m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.d(get_w(i)); });
        emit_1d(std::string(b.base) + "_abs_delta",     b.base_units, b.m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.ad(get_w(i)); });
        emit_1d(std::string(b.base) + "_delta_squared", b.sq_units,   b.sq_m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.sd(get_w(i)); });
        emit_1d(std::string(b.base) + "_dxdt",          b.rate_units, b.rate_m2_units,
                [get_w, &b](std::size_t i) -> const WelfordMoments& { return b.dx(get_w(i)); });

        std::vector<double> rms(N);
        for (std::size_t i = 0; i < N; ++i) rms[i] = b.rms(get_w(i));
        auto rms_ds = grp.createDataSet(std::string(b.base) + "_rms_delta", rms);
        rms_ds.createAttribute("units", b.base_units);
    }

    // Provenance counters
    std::vector<std::size_t> n_frames(N), delta_n(N), dxdt_n(N);
    for (std::size_t i = 0; i < N; ++i) {
        const HydrationGeometryWelfordState& w = get_w(i);
        n_frames[i] = w.n_frames; delta_n[i] = w.delta_n; dxdt_n[i] = w.dxdt_n;
    }
    grp.createDataSet("n_frames_per_atom", n_frames)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("delta_n_per_atom",  delta_n)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("dxdt_n_per_atom",   dxdt_n)
       .createAttribute("units", std::string("frame_count"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_)
       .createAttribute("units", std::string("dimensionless"));
}

}  // namespace nmr
