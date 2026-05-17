#include "HmWelfordTrajectoryResult.h"
#include "HaighMallionResult.h"
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


std::vector<std::type_index> HmWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(HaighMallionResult)) };
}


std::unique_ptr<HmWelfordTrajectoryResult>
HmWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<HmWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_t0_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


void HmWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        const double t0    = ca.hm_shielding_contribution.T0;
        const double t2mag = ca.hm_shielding_contribution.T2Magnitude();

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        const size_t n_new = ta.hm_n_frames + 1;

        MomentsUpdate(ta.hm_t0_mean, ta.hm_t0_m2, t0, n_new);
        MinMaxUpdate(ta.hm_t0_min, ta.hm_t0_max,
                     ta.hm_t0_min_frame, ta.hm_t0_max_frame,
                     t0, frame_idx);

        MomentsUpdate(ta.hm_t2mag_mean, ta.hm_t2mag_m2, t2mag, n_new);
        MinMaxUpdate(ta.hm_t2mag_min, ta.hm_t2mag_max,
                     ta.hm_t2mag_min_frame, ta.hm_t2mag_max_frame,
                     t2mag, frame_idx);

        ta.hm_n_frames = n_new;

        if (prev_valid_[i]) {
            const double delta_x = t0 - prev_t0_[i];
            const size_t dn_new  = ta.hm_t0_delta_n + 1;
            MomentsUpdate(ta.hm_t0_delta_mean, ta.hm_t0_delta_m2,
                          delta_x, dn_new);
            MinMaxUpdate(ta.hm_t0_delta_min, ta.hm_t0_delta_max, delta_x);
            ta.hm_t0_delta_n = dn_new;
        }
        prev_t0_[i] = t0;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


void HmWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                        Trajectory& traj) {
    (void)traj;
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        ta.hm_t0_std       = MomentsStd(ta.hm_t0_m2,       ta.hm_n_frames);
        ta.hm_t2mag_std    = MomentsStd(ta.hm_t2mag_m2,    ta.hm_n_frames);
        ta.hm_t0_delta_std = MomentsStd(ta.hm_t0_delta_m2, ta.hm_t0_delta_n);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "HmWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms");
}


int HmWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 11;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        data[i * K + 0]  = ta.hm_t0_mean;
        data[i * K + 1]  = ta.hm_t0_std;
        data[i * K + 2]  = ta.hm_t0_min;
        data[i * K + 3]  = ta.hm_t0_max;
        data[i * K + 4]  = ta.hm_t2mag_mean;
        data[i * K + 5]  = ta.hm_t2mag_std;
        data[i * K + 6]  = ta.hm_t2mag_min;
        data[i * K + 7]  = ta.hm_t2mag_max;
        data[i * K + 8]  = ta.hm_t0_delta_mean;
        data[i * K + 9]  = ta.hm_t0_delta_std;
        data[i * K + 10] = static_cast<double>(ta.hm_n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/hm_welford.npy",
                            data.data(), N, K);
    return 1;
}


void HmWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();

    std::vector<double> t0_mean(N), t0_std(N), t0_min(N), t0_max(N);
    std::vector<double> t2mag_mean(N), t2mag_std(N), t2mag_min(N), t2mag_max(N);
    std::vector<double> t0_delta_mean(N), t0_delta_std(N);
    std::vector<size_t> n_frames(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        t0_mean[i]       = ta.hm_t0_mean;
        t0_std[i]        = ta.hm_t0_std;
        t0_min[i]        = ta.hm_t0_min;
        t0_max[i]        = ta.hm_t0_max;
        t2mag_mean[i]    = ta.hm_t2mag_mean;
        t2mag_std[i]     = ta.hm_t2mag_std;
        t2mag_min[i]     = ta.hm_t2mag_min;
        t2mag_max[i]     = ta.hm_t2mag_max;
        t0_delta_mean[i] = ta.hm_t0_delta_mean;
        t0_delta_std[i]  = ta.hm_t0_delta_std;
        n_frames[i]      = ta.hm_n_frames;
    }

    auto grp = file.createGroup("/trajectory/hm_welford");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("units",       std::string("Angstrom^-1"));

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
