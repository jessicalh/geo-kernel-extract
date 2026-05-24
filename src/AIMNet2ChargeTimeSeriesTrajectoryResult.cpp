#include "AIMNet2ChargeTimeSeriesTrajectoryResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <typeinfo>

namespace nmr {


std::unique_ptr<AIMNet2ChargeTimeSeriesTrajectoryResult>
AIMNet2ChargeTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<AIMNet2ChargeTimeSeriesTrajectoryResult>();
    r->per_atom_charge_.assign(tp.AtomCount(), std::vector<double>{});
    return r;
}


void AIMNet2ChargeTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_charge_[i].push_back(conf.AtomAt(i).aimnet2_charge);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void AIMNet2ChargeTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                       Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<double>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_charge_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<double>().swap(per_atom_charge_[i]);
        ++atoms_written;
    }
    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<double>(
        std::move(buffer),
        std::type_index(typeid(AIMNet2ChargeTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "AIMNet2ChargeTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) + " frames) AIMNet2 charges to tp dense buffer");
}


void AIMNet2ChargeTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<double>(std::type_index(
            typeid(AIMNet2ChargeTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "AIMNet2ChargeTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/aimnet2_charge_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("irrep_layout", std::string("T0"));
    grp.createAttribute("parity",       std::string("0e"));
    grp.createAttribute("units",        std::string("elementary_charge"));

    std::vector<double> flat(N * T);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            flat[i * T + t] = buffer->At(i, t);
        }
    }
    std::vector<std::size_t> dims = {N, T};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("charge", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);
}

}  // namespace nmr
