#include "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult.h"
#include "TripeptideBackboneShieldingResult.h"
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


std::unique_ptr<TripeptideBackboneMethodTagTimeSeriesTrajectoryResult>
TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        TripeptideBackboneMethodTagTimeSeriesTrajectoryResult>();
    r->per_atom_tag_.assign(tp.AtomCount(), std::vector<std::uint8_t>{});
    return r;
}


void TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_tag_[i].push_back(
            conf.AtomAt(i).tripeptide_bb_method_tag);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// Bounds-check Finalize idempotency.
void TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer =
        std::make_unique<DenseBuffer<std::uint8_t>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_tag_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<std::uint8_t>().swap(per_atom_tag_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<std::uint8_t>(
        std::move(buffer),
        std::type_index(typeid(
            TripeptideBackboneMethodTagTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) uint8 method-tag time-series to tp dense buffer");
}


// (N, T) flat uint8 emission. No inner axis: method_tag is a scalar
// categorical per (atom, frame). Downstream consumers read the legend
// attribute to interpret the integer values.

void TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<std::uint8_t>(std::type_index(typeid(
            TripeptideBackboneMethodTagTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/tripeptide_bb_method_tag_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    // Categorical metadata: explicit legend in the attrs so a downstream
    // reader doesn't need to consult the cpp for tag interpretation.
    grp.createAttribute("units",  std::string("categorical"));
    grp.createAttribute("dtype",  std::string("uint8"));
    grp.createAttribute("legend",
        std::string("0=no_match, 1=opbe_6-31g(d,p)_gaussian, "
                    "2=pbe_6-31g(d,p)_orca_ser"));

    // Flat (N, T) uint8 — already contiguous in DenseBuffer's storage,
    // but we still build a fresh vector via At() to avoid coupling to
    // the buffer's internal layout (same discipline the Vec3 /
    // SphericalTensor TRs follow).
    std::vector<std::uint8_t> flat(N * T);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            flat[i * T + t] = buffer->At(i, t);
        }
    }

    std::vector<std::size_t> dims = {N, T};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<std::uint8_t>("method_tag", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);
}

}  // namespace nmr
