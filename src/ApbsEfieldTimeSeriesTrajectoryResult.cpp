#include "ApbsEfieldTimeSeriesTrajectoryResult.h"
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


std::unique_ptr<ApbsEfieldTimeSeriesTrajectoryResult>
ApbsEfieldTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<ApbsEfieldTimeSeriesTrajectoryResult>();
    r->per_atom_efield_.assign(tp.AtomCount(), std::vector<Vec3>{});
    return r;
}


void ApbsEfieldTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_efield_[i].push_back(conf.AtomAt(i).apbs_efield);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void ApbsEfieldTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<Vec3>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_efield_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<Vec3>().swap(per_atom_efield_[i]);
        ++atoms_written;
    }
    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<Vec3>(
        std::move(buffer),
        std::type_index(typeid(ApbsEfieldTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "ApbsEfieldTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) + " frames) APBS E-field Vec3 to tp dense buffer");
}


void ApbsEfieldTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<Vec3>(std::type_index(
            typeid(ApbsEfieldTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "ApbsEfieldTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/apbs_efield_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("irrep_layout",  std::string("x,y,z"));
    grp.createAttribute("normalization", std::string("cartesian"));
    grp.createAttribute("parity",        std::string("1o"));
    grp.createAttribute("units",         std::string("V/Angstrom"));

    // (N, T, 3) via explicit .x()/.y()/.z() access.
    std::vector<double> flat(N * T * 3);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const Vec3& v = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 3;
            flat[base + 0] = v.x();
            flat[base + 1] = v.y();
            flat[base + 2] = v.z();
        }
    }
    std::vector<std::size_t> dims = {N, T, std::size_t(3)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);
}

}  // namespace nmr
