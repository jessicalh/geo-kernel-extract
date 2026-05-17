#include "EeqWelfordTrajectoryResult.h"
#include "EeqResult.h"
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


std::vector<std::type_index> EeqWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(EeqResult)) };
}


std::unique_ptr<EeqWelfordTrajectoryResult>
EeqWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<EeqWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_value_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


void EeqWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const double q = conf.AtomAt(i).eeq_charge;

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        const size_t n_new = ta.eeq_charge_n_frames + 1;

        MomentsUpdate(ta.eeq_charge_mean, ta.eeq_charge_m2, q, n_new);
        MinMaxUpdate(ta.eeq_charge_min, ta.eeq_charge_max,
                     ta.eeq_charge_min_frame, ta.eeq_charge_max_frame,
                     q, frame_idx);

        ta.eeq_charge_n_frames = n_new;

        if (prev_valid_[i]) {
            const double delta_x = q - prev_value_[i];
            const size_t dn_new  = ta.eeq_charge_delta_n + 1;
            MomentsUpdate(ta.eeq_charge_delta_mean, ta.eeq_charge_delta_m2,
                          delta_x, dn_new);
            MinMaxUpdate(ta.eeq_charge_delta_min, ta.eeq_charge_delta_max, delta_x);
            ta.eeq_charge_delta_n = dn_new;
        }
        prev_value_[i] = q;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


void EeqWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                         Trajectory& traj) {
    (void)traj;
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        ta.eeq_charge_std       = MomentsStd(ta.eeq_charge_m2,       ta.eeq_charge_n_frames);
        ta.eeq_charge_delta_std = MomentsStd(ta.eeq_charge_delta_m2, ta.eeq_charge_delta_n);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "EeqWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms");
}


int EeqWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 7;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        data[i * K + 0] = ta.eeq_charge_mean;
        data[i * K + 1] = ta.eeq_charge_std;
        data[i * K + 2] = ta.eeq_charge_min;
        data[i * K + 3] = ta.eeq_charge_max;
        data[i * K + 4] = ta.eeq_charge_delta_mean;
        data[i * K + 5] = ta.eeq_charge_delta_std;
        data[i * K + 6] = static_cast<double>(ta.eeq_charge_n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/eeq_welford.npy",
                            data.data(), N, K);
    return 1;
}


void EeqWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();

    std::vector<double> mean(N), std_(N), min_(N), max_(N);
    std::vector<double> delta_mean(N), delta_std(N);
    std::vector<size_t> n_frames(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        mean[i]       = ta.eeq_charge_mean;
        std_[i]       = ta.eeq_charge_std;
        min_[i]       = ta.eeq_charge_min;
        max_[i]       = ta.eeq_charge_max;
        delta_mean[i] = ta.eeq_charge_delta_mean;
        delta_std[i]  = ta.eeq_charge_delta_std;
        n_frames[i]   = ta.eeq_charge_n_frames;
    }

    auto grp = file.createGroup("/trajectory/eeq_welford");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("units",       std::string("elementary_charge"));

    grp.createDataSet("mean",       mean);
    grp.createDataSet("std",        std_);
    grp.createDataSet("min",        min_);
    grp.createDataSet("max",        max_);
    grp.createDataSet("delta_mean", delta_mean);
    grp.createDataSet("delta_std",  delta_std);
    grp.createDataSet("n_frames_per_atom", n_frames);
}

}  // namespace nmr
