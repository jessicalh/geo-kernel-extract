#include "McConnellWelfordTrajectoryResult.h"
#include "McConnellResult.h"
#include "TrajectoryProtein.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Types.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <limits>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index> McConnellWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(McConnellResult)) };
}


std::unique_ptr<McConnellWelfordTrajectoryResult>
McConnellWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<McConnellWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_t0_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


void McConnellWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        const double t0    = ca.mc_shielding_contribution.T0;
        const double t2mag = ca.mc_shielding_contribution.T2Magnitude();

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        McConnellWelfordState& w = ta.mc_welford;
        const size_t n_new = w.n_frames + 1;

        WelfordUpdate(w.t0,          t0,    n_new, frame_idx);
        WelfordUpdate(w.t2magnitude, t2mag, n_new, frame_idx);
        w.n_frames = n_new;

        if (prev_valid_[i]) {
            const double delta_x = t0 - prev_t0_[i];
            const size_t dn_new  = w.delta_n + 1;
            WelfordUpdate(w.t0_delta, delta_x, dn_new, frame_idx);
            w.delta_n = dn_new;
        }
        prev_t0_[i] = t0;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


void McConnellWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                Trajectory& traj) {
    (void)traj;
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        McConnellWelfordState& w = tp.MutableAtomAt(i).mc_welford;
        WelfordFinalize(w.t0,          w.n_frames);
        WelfordFinalize(w.t2magnitude, w.n_frames);
        WelfordFinalize(w.t0_delta,    w.delta_n);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "McConnellWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms");
}


int McConnellWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 11;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const McConnellWelfordState& w = tp.AtomAt(i).mc_welford;
        data[i * K + 0]  = w.t0.mean;
        data[i * K + 1]  = w.t0.std;
        data[i * K + 2]  = w.t0.min;
        data[i * K + 3]  = w.t0.max;
        data[i * K + 4]  = w.t2magnitude.mean;
        data[i * K + 5]  = w.t2magnitude.std;
        data[i * K + 6]  = w.t2magnitude.min;
        data[i * K + 7]  = w.t2magnitude.max;
        data[i * K + 8]  = w.t0_delta.mean;
        data[i * K + 9]  = w.t0_delta.std;
        data[i * K + 10] = static_cast<double>(w.n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/mc_welford.npy",
                            data.data(), N, K);
    return 1;
}


void McConnellWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();

    std::vector<double> t0_mean(N), t0_std(N), t0_min(N), t0_max(N);
    std::vector<double> t2mag_mean(N), t2mag_std(N), t2mag_min(N), t2mag_max(N);
    std::vector<double> t0_delta_mean(N), t0_delta_std(N);
    std::vector<size_t> n_frames(N);

    for (size_t i = 0; i < N; ++i) {
        const McConnellWelfordState& w = tp.AtomAt(i).mc_welford;
        t0_mean[i]       = w.t0.mean;
        t0_std[i]        = w.t0.std;
        t0_min[i]        = w.t0.min;
        t0_max[i]        = w.t0.max;
        t2mag_mean[i]    = w.t2magnitude.mean;
        t2mag_std[i]     = w.t2magnitude.std;
        t2mag_min[i]     = w.t2magnitude.min;
        t2mag_max[i]     = w.t2magnitude.max;
        t0_delta_mean[i] = w.t0_delta.mean;
        t0_delta_std[i]  = w.t0_delta.std;
        n_frames[i]      = w.n_frames;
    }

    auto grp = file.createGroup("/trajectory/mc_welford");
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("units",       std::string("Angstrom^-3"));

    grp.createDataSet("t0_mean",       t0_mean);
    grp.createDataSet("t0_std",        t0_std);
    grp.createDataSet("t0_min",        t0_min);
    grp.createDataSet("t0_max",        t0_max);
    grp.createDataSet("t2mag_mean",    t2mag_mean);
    grp.createDataSet("t2mag_std",     t2mag_std);
    grp.createDataSet("t2mag_min",     t2mag_min);
    grp.createDataSet("t2mag_max",     t2mag_max);
    grp.createDataSet("t0_delta_mean", t0_delta_mean);
    grp.createDataSet("t0_delta_std",  t0_delta_std);
    grp.createDataSet("n_frames_per_atom", n_frames);
}

}  // namespace nmr
