#include "SasaWelfordTrajectoryResult.h"
#include "SasaResult.h"
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


std::vector<std::type_index> SasaWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(SasaResult)) };
}


std::unique_ptr<SasaWelfordTrajectoryResult>
SasaWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<SasaWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_value_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


void SasaWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        const double s = conf.AtomAt(i).atom_sasa;

        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        const size_t n_new = ta.sasa_n_frames + 1;

        MomentsUpdate(ta.sasa_mean, ta.sasa_m2, s, n_new);
        MinMaxUpdate(ta.sasa_min, ta.sasa_max,
                     ta.sasa_min_frame, ta.sasa_max_frame,
                     s, frame_idx);

        ta.sasa_n_frames = n_new;

        if (prev_valid_[i]) {
            const double delta_x = s - prev_value_[i];
            const size_t dn_new  = ta.sasa_delta_n + 1;
            MomentsUpdate(ta.sasa_delta_mean, ta.sasa_delta_m2,
                          delta_x, dn_new);
            MinMaxUpdate(ta.sasa_delta_min, ta.sasa_delta_max, delta_x);
            ta.sasa_delta_n = dn_new;
        }
        prev_value_[i] = s;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


void SasaWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                          Trajectory& traj) {
    (void)traj;
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        ta.sasa_std       = MomentsStd(ta.sasa_m2,       ta.sasa_n_frames);
        ta.sasa_delta_std = MomentsStd(ta.sasa_delta_m2, ta.sasa_delta_n);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "SasaWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms");
}


int SasaWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 7;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        data[i * K + 0] = ta.sasa_mean;
        data[i * K + 1] = ta.sasa_std;
        data[i * K + 2] = ta.sasa_min;
        data[i * K + 3] = ta.sasa_max;
        data[i * K + 4] = ta.sasa_delta_mean;
        data[i * K + 5] = ta.sasa_delta_std;
        data[i * K + 6] = static_cast<double>(ta.sasa_n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/sasa_welford.npy",
                            data.data(), N, K);
    return 1;
}


void SasaWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();

    std::vector<double> mean(N), std_(N), min_(N), max_(N);
    std::vector<double> delta_mean(N), delta_std(N);
    std::vector<size_t> n_frames(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        mean[i]       = ta.sasa_mean;
        std_[i]       = ta.sasa_std;
        min_[i]       = ta.sasa_min;
        max_[i]       = ta.sasa_max;
        delta_mean[i] = ta.sasa_delta_mean;
        delta_std[i]  = ta.sasa_delta_std;
        n_frames[i]   = ta.sasa_n_frames;
    }

    auto grp = file.createGroup("/trajectory/sasa_welford");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("units",       std::string("Angstrom^2"));

    grp.createDataSet("mean",       mean);
    grp.createDataSet("std",        std_);
    grp.createDataSet("min",        min_);
    grp.createDataSet("max",        max_);
    grp.createDataSet("delta_mean", delta_mean);
    grp.createDataSet("delta_std",  delta_std);
    grp.createDataSet("n_frames_per_atom", n_frames);
}

}  // namespace nmr
