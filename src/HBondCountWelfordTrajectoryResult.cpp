#include "HBondCountWelfordTrajectoryResult.h"
#include "HBondResult.h"
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


std::vector<std::type_index> HBondCountWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(HBondResult)) };
}


std::unique_ptr<HBondCountWelfordTrajectoryResult>
HBondCountWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<HBondCountWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_value_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


void HBondCountWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const double c = static_cast<double>(conf.AtomAt(i).hbond_count_within_3_5A);

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        HBondCountWelfordState& w = ta.hbond_count_welford;
        const size_t n_new = w.n_frames + 1;

        WelfordUpdate(w.count, c, n_new, frame_idx);
        w.n_frames = n_new;

        if (prev_valid_[i]) {
            const double delta_x = c - prev_value_[i];
            const size_t dn_new  = w.delta_n + 1;
            WelfordUpdate(w.count_delta, delta_x, dn_new, frame_idx);
            w.delta_n = dn_new;
        }
        prev_value_[i] = c;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


void HBondCountWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                 Trajectory& traj) {
    (void)traj;
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        HBondCountWelfordState& w = tp.MutableAtomAt(i).hbond_count_welford;
        WelfordFinalize(w.count,       w.n_frames);
        WelfordFinalize(w.count_delta, w.delta_n);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "HBondCountWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms");
}


int HBondCountWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 7;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const HBondCountWelfordState& w = tp.AtomAt(i).hbond_count_welford;
        data[i * K + 0] = w.count.mean;
        data[i * K + 1] = w.count.std;
        data[i * K + 2] = w.count.min;
        data[i * K + 3] = w.count.max;
        data[i * K + 4] = w.count_delta.mean;
        data[i * K + 5] = w.count_delta.std;
        data[i * K + 6] = static_cast<double>(w.n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/hbond_count_welford.npy",
                            data.data(), N, K);
    return 1;
}


void HBondCountWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();

    std::vector<double> mean(N), std_(N), min_(N), max_(N);
    std::vector<double> delta_mean(N), delta_std(N);
    std::vector<size_t> n_frames(N);

    for (size_t i = 0; i < N; ++i) {
        const HBondCountWelfordState& w = tp.AtomAt(i).hbond_count_welford;
        mean[i]       = w.count.mean;
        std_[i]       = w.count.std;
        min_[i]       = w.count.min;
        max_[i]       = w.count.max;
        delta_mean[i] = w.count_delta.mean;
        delta_std[i]  = w.count_delta.std;
        n_frames[i]   = w.n_frames;
    }

    auto grp = file.createGroup("/trajectory/hbond_count_welford");
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("units",       std::string("pairs"));
    grp.createAttribute("source_radius_A", 3.5);

    grp.createDataSet("mean",       mean);
    grp.createDataSet("std",        std_);
    grp.createDataSet("min",        min_);
    grp.createDataSet("max",        max_);
    grp.createDataSet("delta_mean", delta_mean);
    grp.createDataSet("delta_std",  delta_std);
    grp.createDataSet("n_frames_per_atom", n_frames);
}

}  // namespace nmr
