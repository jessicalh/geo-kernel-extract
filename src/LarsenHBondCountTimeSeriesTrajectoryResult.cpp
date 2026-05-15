#include "LarsenHBondCountTimeSeriesTrajectoryResult.h"
#include "LarsenHBondShieldingResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <limits>
#include <typeinfo>

namespace nmr {


std::unique_ptr<LarsenHBondCountTimeSeriesTrajectoryResult>
LarsenHBondCountTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        LarsenHBondCountTimeSeriesTrajectoryResult>();
    r->per_atom_count_.assign(tp.AtomCount(), std::vector<int>{});
    return r;
}


void LarsenHBondCountTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    // "Absent, not faked" provenance: record whether the source
    // calculator (LarsenHBondShieldingResult) attached this frame.
    // Int-typed TR — no NaN-fill of the data (0 is the existing
    // "no pairs" sentinel) but the mask still gates H5 group
    // emission so the group is absent rather than all-zero when the
    // source calc never attached.
    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<LarsenHBondShieldingResult>();
    source_present_per_frame_.push_back(source_attached ? 1u : 0u);

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_count_[i].push_back(
            conf.AtomAt(i).larsen_hbond_n_pairs);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void LarsenHBondCountTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<int>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_count_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<int>().swap(per_atom_count_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<int>(
        std::move(buffer),
        std::type_index(typeid(
            LarsenHBondCountTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "LarsenHBondCountTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) int H-bond-count time-series to tp dense buffer");
}


void LarsenHBondCountTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<int>(std::type_index(typeid(
            LarsenHBondCountTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "LarsenHBondCountTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    // "Absent, not faked" — if the source ConformationResult was not
    // attached in any frame, skip emission. Group existence ⇒ source
    // ran in ≥1 frame. Downstream readers MUST tolerate group absence
    // for conditionally-attached-source TRs.
    std::size_t source_present_count = 0;
    for (auto v : source_present_per_frame_)
        if (v) ++source_present_count;
    if (source_present_count == 0) {
        OperationLog::Warn(
            "LarsenHBondCountTimeSeriesTrajectoryResult::WriteH5Group",
            "LarsenHBondShieldingResult was not attached in any of "
            + std::to_string(source_present_per_frame_.size()) +
            " frames; skipping /trajectory/larsen_hbond_count_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/larsen_hbond_count_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("units", std::string("pairs"));
    grp.createAttribute("dtype", std::string("int32"));
    grp.createAttribute("description",
        std::string("Per-atom H-bond pair count: sum over 1°HB + 2°HB + "
                    "1°HαB + 2°HαB contributions this frame "
                    "(LarsenContribDispatch / Larsen Table 2)."));

    std::vector<int> flat(N * T);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            flat[i * T + t] = buffer->At(i, t);
        }
    }

    std::vector<std::size_t> dims = {N, T};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<int>("count", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);

    // Provenance mask: per-frame source-attached flags.
    grp.createDataSet("source_attached_per_frame", source_present_per_frame_);
}

}  // namespace nmr
