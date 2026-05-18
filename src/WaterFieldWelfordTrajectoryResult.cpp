#include "WaterFieldWelfordTrajectoryResult.h"
#include "WaterFieldResult.h"
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
WaterFieldWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(WaterFieldResult)) };
}


std::unique_ptr<WaterFieldWelfordTrajectoryResult>
WaterFieldWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<WaterFieldWelfordTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    r->prev_efield_mag_.assign(N, 0.0);
    r->prev_n_first_.assign(N, 0.0);
    r->prev_n_second_.assign(N, 0.0);
    r->prev_valid_.assign(N, false);
    r->prev_time_.assign(N, 0.0);
    return r;
}


// ── Compute ──────────────────────────────────────────────────────

void WaterFieldWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        WaterFieldWelfordState& w = tp.MutableAtomAt(i).water_field_welford;
        const std::size_t n_new = w.n_frames + 1;

        // E-field (Vec3) per-component + magnitude
        const Vec3& E       = a.water_efield;
        const Vec3& Efirst  = a.water_efield_first;
        const double Emag      = E.norm();
        const double Efirstmag = Efirst.norm();
        WelfordUpdate(w.efield[0], E.x(), n_new, frame_idx);
        WelfordUpdate(w.efield[1], E.y(), n_new, frame_idx);
        WelfordUpdate(w.efield[2], E.z(), n_new, frame_idx);
        WelfordUpdate(w.efield_magnitude, Emag, n_new, frame_idx);
        WelfordUpdate(w.efield_first[0], Efirst.x(), n_new, frame_idx);
        WelfordUpdate(w.efield_first[1], Efirst.y(), n_new, frame_idx);
        WelfordUpdate(w.efield_first[2], Efirst.z(), n_new, frame_idx);
        WelfordUpdate(w.efield_first_magnitude, Efirstmag, n_new, frame_idx);

        // EFG SphericalTensor (total + first-shell). T0 channel intentionally
        // absent — WaterFieldResult traceless-projects V before Decompose,
        // so T0 = 0 by construction. T1[3] + T2[5] + |T2| carry the signal.
        const SphericalTensor& efg  = a.water_efg_spherical;
        const SphericalTensor& efgf = a.water_efg_first_spherical;
        WelfordUpdate(w.efg_t2magnitude,        efg.T2Magnitude(), n_new, frame_idx);
        for (std::size_t k = 0; k < 3; ++k)
            WelfordUpdate(w.efg_t1[k],          efg.T1[k],        n_new, frame_idx);
        for (std::size_t k = 0; k < 5; ++k)
            WelfordUpdate(w.efg_t2[k],          efg.T2[k],        n_new, frame_idx);
        WelfordUpdate(w.efg_first_t2magnitude,  efgf.T2Magnitude(),n_new, frame_idx);
        for (std::size_t k = 0; k < 3; ++k)
            WelfordUpdate(w.efg_first_t1[k],    efgf.T1[k],       n_new, frame_idx);
        for (std::size_t k = 0; k < 5; ++k)
            WelfordUpdate(w.efg_first_t2[k],    efgf.T2[k],       n_new, frame_idx);

        // Shell occupancy counts (int → double for Welford)
        const double n_first_d  = static_cast<double>(a.water_n_first);
        const double n_second_d = static_cast<double>(a.water_n_second);
        WelfordUpdate(w.n_first,  n_first_d,  n_new, frame_idx);
        WelfordUpdate(w.n_second, n_second_d, n_new, frame_idx);

        w.n_frames = n_new;

        // Delta variants on the 4 primary scalars
        if (prev_valid_[i]) {
            const std::size_t dn_new = w.delta_n + 1;
            auto upd_deltas = [&](double curr, double prev,
                                  WelfordMoments& d, WelfordMoments& ad,
                                  WelfordMoments& sd) {
                const double delta = curr - prev;
                WelfordUpdate(d,  delta,                dn_new, frame_idx);
                WelfordUpdate(ad, std::abs(delta),      dn_new, frame_idx);
                WelfordUpdate(sd, delta * delta,        dn_new, frame_idx);
            };
            upd_deltas(Emag,        prev_efield_mag_[i],
                       w.efield_magnitude_delta,
                       w.efield_magnitude_abs_delta,
                       w.efield_magnitude_delta_squared);
            upd_deltas(n_first_d,   prev_n_first_[i],
                       w.n_first_delta, w.n_first_abs_delta, w.n_first_delta_squared);
            upd_deltas(n_second_d,  prev_n_second_[i],
                       w.n_second_delta, w.n_second_abs_delta, w.n_second_delta_squared);
            w.delta_n = dn_new;

            constexpr double MIN_DT_PS = 1e-12;
            const double dt = time_ps - prev_time_[i];
            if (std::abs(dt) > MIN_DT_PS) {
                const std::size_t dxn = w.dxdt_n + 1;
                WelfordUpdate(w.efield_magnitude_dxdt, (Emag       - prev_efield_mag_[i]) / dt, dxn, frame_idx);
                WelfordUpdate(w.n_first_dxdt,          (n_first_d  - prev_n_first_[i])    / dt, dxn, frame_idx);
                WelfordUpdate(w.n_second_dxdt,         (n_second_d - prev_n_second_[i])   / dt, dxn, frame_idx);
                w.dxdt_n = dxn;
            }
        }
        prev_efield_mag_[i] = Emag;
        prev_n_first_[i]    = n_first_d;
        prev_n_second_[i]   = n_second_d;
        prev_time_[i]       = time_ps;
        prev_valid_[i]      = true;
    }
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void WaterFieldWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                 Trajectory& traj) {
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        WaterFieldWelfordState& w = tp.MutableAtomAt(i).water_field_welford;

        for (std::size_t k = 0; k < 3; ++k) WelfordFinalize(w.efield[k],       w.n_frames);
        WelfordFinalize(w.efield_magnitude,       w.n_frames);
        for (std::size_t k = 0; k < 3; ++k) WelfordFinalize(w.efield_first[k], w.n_frames);
        WelfordFinalize(w.efield_first_magnitude, w.n_frames);

        for (std::size_t k = 0; k < 3; ++k) WelfordFinalize(w.efg_t1[k], w.n_frames);
        for (std::size_t k = 0; k < 5; ++k) WelfordFinalize(w.efg_t2[k], w.n_frames);
        WelfordFinalize(w.efg_t2magnitude,  w.n_frames);
        for (std::size_t k = 0; k < 3; ++k) WelfordFinalize(w.efg_first_t1[k], w.n_frames);
        for (std::size_t k = 0; k < 5; ++k) WelfordFinalize(w.efg_first_t2[k], w.n_frames);
        WelfordFinalize(w.efg_first_t2magnitude,  w.n_frames);

        WelfordFinalize(w.n_first,  w.n_frames);
        WelfordFinalize(w.n_second, w.n_frames);

        // Delta-variant finalize
        auto fin_deltas = [&](WelfordMoments& d, WelfordMoments& ad,
                              WelfordMoments& sd, WelfordMoments& dx,
                              double& rms) {
            WelfordFinalize(d,  w.delta_n);
            WelfordFinalize(ad, w.delta_n);
            WelfordFinalize(sd, w.delta_n);
            WelfordFinalize(dx, w.dxdt_n);
            rms = (w.delta_n == 0) ? std::nan("") : std::sqrt(sd.mean);
        };
        fin_deltas(w.efield_magnitude_delta, w.efield_magnitude_abs_delta,
                   w.efield_magnitude_delta_squared, w.efield_magnitude_dxdt,
                   w.efield_magnitude_rms_delta);
        fin_deltas(w.n_first_delta,          w.n_first_abs_delta,
                   w.n_first_delta_squared,  w.n_first_dxdt,
                   w.n_first_rms_delta);
        fin_deltas(w.n_second_delta,         w.n_second_abs_delta,
                   w.n_second_delta_squared, w.n_second_dxdt,
                   w.n_second_rms_delta);
    }

    const auto& times = traj.FrameTimes();
    if (times.size() >= 2)
        mean_dt_ps_ = (times.back() - times.front()) /
                      static_cast<double>(times.size() - 1);
    const auto& fidx = traj.FrameIndices();
    if (!fidx.empty()) frame_index_range_ = {fidx.front(), fidx.back()};

    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "WaterFieldWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(N) + " atoms");
}


// ── WriteH5Group ─────────────────────────────────────────────────

void WaterFieldWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    auto grp = file.createGroup("/trajectory/water_field_welford");

    grp.createAttribute("result_name",        Name());
    grp.createAttribute("n_frames",           n_frames_);
    grp.createAttribute("finalized",          finalized_);
    grp.createAttribute("ddof",               static_cast<int>(1));
    grp.createAttribute("mean_dt_ps",         mean_dt_ps_);
    grp.createAttribute("frame_index_range",  frame_index_range_);
    grp.createAttribute("irrep_layout_t1",    std::string("v_x,v_y,v_z"));
    grp.createAttribute("irrep_layout_t2",    std::string("m-2,m-1,m0,m+1,m+2"));
    // Group-level `units` describes the primary value channel (E-field).
    // Per-dataset `units` attributes are authoritative — the group holds
    // V/Å (E-field), V/Å² (EFG), and dimensionless (counts) datasets.
    grp.createAttribute("units",              std::string("V/Angstrom"));
    // EFG T0 emission absent — structurally zero (WaterFieldResult
    // traceless-projects V before Decompose). Same flag as TimeSeries TR.
    grp.createAttribute("efg_t0_structural_zero", true);

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
        auto put = [&](const std::string& n, const std::vector<double>& v,
                       const std::string& u) {
            auto ds = grp.createDataSet(n, v); ds.createAttribute("units", u);
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

    auto emit_2d = [&](const std::string& prefix, std::size_t K,
                       const std::string& base_units,
                       const std::string& m2_units,
                       std::function<const WelfordMoments&(std::size_t, std::size_t)> get) {
        std::vector<double> mean(N * K), m2(N * K), std_(N * K), min_(N * K), max_(N * K);
        std::vector<std::size_t> minf(N * K), maxf(N * K);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t k = 0; k < K; ++k) {
                const WelfordMoments& w = get(i, k);
                mean[i*K+k] = w.mean; m2[i*K+k] = w.m2; std_[i*K+k] = w.std;
                min_[i*K+k] = w.min;  max_[i*K+k] = w.max;
                minf[i*K+k] = w.min_frame; maxf[i*K+k] = w.max_frame;
            }
        }
        HighFive::DataSpace space({N, K});
        auto put = [&](const std::string& n, const std::vector<double>& v,
                       const std::string& u) {
            auto ds = grp.createDataSet<double>(n, space);
            ds.write_raw(v.data()); ds.createAttribute("units", u);
        };
        put(prefix + "_mean", mean, base_units);
        put(prefix + "_m2",   m2,   m2_units);
        put(prefix + "_std",  std_, base_units);
        put(prefix + "_min",  min_, base_units);
        put(prefix + "_max",  max_, base_units);
        auto dminf = grp.createDataSet<std::size_t>(prefix + "_min_frame", space);
        dminf.write_raw(minf.data());
        dminf.createAttribute("units", std::string("frame_index"));
        auto dmaxf = grp.createDataSet<std::size_t>(prefix + "_max_frame", space);
        dmaxf.write_raw(maxf.data());
        dmaxf.createAttribute("units", std::string("frame_index"));
    };

    const std::string kE   = "V/Angstrom";
    const std::string kEsq = "V^2/Angstrom^2";
    const std::string kEFG   = "V/Angstrom^2";
    const std::string kEFGsq = "V^2/Angstrom^4";
    const std::string kCount   = "dimensionless";
    const std::string kCountSq = "dimensionless";
    const std::string kErate    = "V/Angstrom/ps";
    const std::string kErateSq  = "V^2/Angstrom^2/ps^2";
    const std::string kEFGrate   = "V/Angstrom^2/ps";
    const std::string kEFGrateSq = "V^2/Angstrom^4/ps^2";
    const std::string kCountRate   = "count/ps";
    const std::string kCountRateSq = "count^2/ps^2";

    // Helpers: capture tp by reference so lambdas can read each per-atom
    // WaterFieldWelfordState through `tp.AtomAt(i).water_field_welford`.
    auto get_w = [&tp](std::size_t i) -> const WaterFieldWelfordState& {
        return tp.AtomAt(i).water_field_welford;
    };

    // E-field per-component (x/y/z) + magnitude  + first-shell variants
    static const std::array<const char*, 3> kXyz = {"x", "y", "z"};
    for (std::size_t k = 0; k < 3; ++k) {
        emit_1d(std::string("efield_") + kXyz[k], kE, kEsq,
                [get_w, k](std::size_t i) -> const WelfordMoments& {
                    return get_w(i).efield[k];
                });
        emit_1d(std::string("efield_first_") + kXyz[k], kE, kEsq,
                [get_w, k](std::size_t i) -> const WelfordMoments& {
                    return get_w(i).efield_first[k];
                });
    }
    emit_1d("efield_magnitude",       kE, kEsq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).efield_magnitude;
            });
    emit_1d("efield_first_magnitude", kE, kEsq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).efield_first_magnitude;
            });

    // EFG: T0 scalar + per-component T1[3] / T2[5] + |T2|, both total and first-shell
    // EFG T0 emission intentionally absent — structurally zero
    // (see WaterFieldWelfordState class comment).
    emit_1d("efg_t2magnitude", kEFG, kEFGsq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).efg_t2magnitude;
            });
    emit_2d("efg_t1", 3, kEFG, kEFGsq,
            [get_w](std::size_t i, std::size_t k) -> const WelfordMoments& {
                return get_w(i).efg_t1[k];
            });
    emit_2d("efg_t2", 5, kEFG, kEFGsq,
            [get_w](std::size_t i, std::size_t k) -> const WelfordMoments& {
                return get_w(i).efg_t2[k];
            });
    emit_1d("efg_first_t2magnitude", kEFG, kEFGsq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).efg_first_t2magnitude;
            });
    emit_2d("efg_first_t1", 3, kEFG, kEFGsq,
            [get_w](std::size_t i, std::size_t k) -> const WelfordMoments& {
                return get_w(i).efg_first_t1[k];
            });
    emit_2d("efg_first_t2", 5, kEFG, kEFGsq,
            [get_w](std::size_t i, std::size_t k) -> const WelfordMoments& {
                return get_w(i).efg_first_t2[k];
            });

    // Shell-occupancy counts
    emit_1d("n_first",  kCount, kCountSq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).n_first;
            });
    emit_1d("n_second", kCount, kCountSq,
            [get_w](std::size_t i) -> const WelfordMoments& {
                return get_w(i).n_second;
            });

    // Delta variants on the 3 primary scalars (efield_magnitude, n_first,
    // n_second). efg_t0 deltas omitted — channel structurally zero.
    //
    // Each bundle carries:
    //   base       — H5 dataset name prefix
    //   base_units — value unit (mean / std / min / max / abs_delta)
    //   m2_units   — sum-of-squared-deviations unit (= value²)
    //   sq_units   — *_delta_squared *value* unit (= value², the squared Δ)
    //   sq_m2_units — *_delta_squared *m2* unit (= value⁴, Welford M2 of Δ²)
    //   rate_units / rate_m2_units — dxdt / dxdt-m2 units
    //   d/ad/sd/dx/rms — channel getters
    //
    // Codex 2026-05-18: the `sq_m2_units` field replaces a fragile
    // string-dispatch lookup ("is base == 'efg_t0'?" etc.) that fell
    // through to "dimensionless" for the count channels — semantically
    // OK by luck but a typed lookup chain instead of an explicit unit.
    struct DeltaBundle {
        const char* base;
        const std::string& base_units;
        const std::string& m2_units;
        const std::string& sq_units;
        const std::string& sq_m2_units;
        const std::string& rate_units;
        const std::string& rate_m2_units;
        std::function<const WelfordMoments&(const WaterFieldWelfordState&)> d;
        std::function<const WelfordMoments&(const WaterFieldWelfordState&)> ad;
        std::function<const WelfordMoments&(const WaterFieldWelfordState&)> sd;
        std::function<const WelfordMoments&(const WaterFieldWelfordState&)> dx;
        std::function<double(const WaterFieldWelfordState&)> rms;
    };
    const std::string kEsq4    = "V^4/Angstrom^4";    // (V/Å)² -> m2 of squared-Δ
    const std::string kCount4  = "dimensionless";     // count⁴ is still dimensionless
    const std::vector<DeltaBundle> bundles = {
        { "efield_magnitude", kE,     kEsq,     kEsq,     kEsq4,    kErate,     kErateSq,
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.efield_magnitude_delta; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.efield_magnitude_abs_delta; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.efield_magnitude_delta_squared; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.efield_magnitude_dxdt; },
          [](const WaterFieldWelfordState& w) -> double { return w.efield_magnitude_rms_delta; } },
        { "n_first",          kCount, kCountSq, kCountSq, kCount4,  kCountRate, kCountRateSq,
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_first_delta; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_first_abs_delta; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_first_delta_squared; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_first_dxdt; },
          [](const WaterFieldWelfordState& w) -> double { return w.n_first_rms_delta; } },
        { "n_second",         kCount, kCountSq, kCountSq, kCount4,  kCountRate, kCountRateSq,
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_second_delta; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_second_abs_delta; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_second_delta_squared; },
          [](const WaterFieldWelfordState& w) -> const WelfordMoments& { return w.n_second_dxdt; },
          [](const WaterFieldWelfordState& w) -> double { return w.n_second_rms_delta; } },
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
        const WaterFieldWelfordState& w = get_w(i);
        n_frames[i] = w.n_frames;
        delta_n[i]  = w.delta_n;
        dxdt_n[i]   = w.dxdt_n;
    }
    auto nf_ds = grp.createDataSet("n_frames_per_atom", n_frames);
    nf_ds.createAttribute("units", std::string("frame_count"));
    auto dn_ds = grp.createDataSet("delta_n_per_atom",  delta_n);
    dn_ds.createAttribute("units", std::string("frame_count"));
    auto dxdtn_ds = grp.createDataSet("dxdt_n_per_atom", dxdt_n);
    dxdtn_ds.createAttribute("units", std::string("frame_count"));
}

}  // namespace nmr
