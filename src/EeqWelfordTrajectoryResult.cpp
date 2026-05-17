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
        EeqWelfordState& w = ta.eeq_welford;
        const size_t n_new = w.n_frames + 1;

        WelfordUpdate(w.charge, q, n_new, frame_idx);
        w.n_frames = n_new;

        if (prev_valid_[i]) {
            const double delta_x = q - prev_value_[i];
            const size_t dn_new  = w.delta_n + 1;
            WelfordUpdate(w.charge_delta, delta_x, dn_new, frame_idx);
            w.delta_n = dn_new;
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
        EeqWelfordState& w = tp.MutableAtomAt(i).eeq_welford;
        WelfordFinalize(w.charge,       w.n_frames);
        WelfordFinalize(w.charge_delta, w.delta_n);
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
        const EeqWelfordState& w = tp.AtomAt(i).eeq_welford;
        data[i * K + 0] = w.charge.mean;
        data[i * K + 1] = w.charge.std;
        data[i * K + 2] = w.charge.min;
        data[i * K + 3] = w.charge.max;
        data[i * K + 4] = w.charge_delta.mean;
        data[i * K + 5] = w.charge_delta.std;
        data[i * K + 6] = static_cast<double>(w.n_frames);
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
        const EeqWelfordState& w = tp.AtomAt(i).eeq_welford;
        mean[i]       = w.charge.mean;
        std_[i]       = w.charge.std;
        min_[i]       = w.charge.min;
        max_[i]       = w.charge.max;
        delta_mean[i] = w.charge_delta.mean;
        delta_std[i]  = w.charge_delta.std;
        n_frames[i]   = w.n_frames;
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
