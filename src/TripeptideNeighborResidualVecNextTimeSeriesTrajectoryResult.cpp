#include "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborShieldingResult.h"
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


std::unique_ptr<TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult>
TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult>();
    r->per_atom_residual_.assign(tp.AtomCount(), std::vector<Vec3>{});
    return r;
}


void TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    // "Absent, not faked" provenance: record whether the source
    // calculator (TripeptideNeighborShieldingResult) attached this
    // frame. When absent, the in-memory field is zero-default — we
    // capture that here but NaN-fill at WriteH5Group time so
    // downstream readers can distinguish "no measurement" from "real
    // measurement = 0."
    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<TripeptideNeighborShieldingResult>();
    source_present_per_frame_.push_back(source_attached ? 1u : 0u);

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_residual_[i].push_back(
            conf.AtomAt(i).tripeptide_neighbor_residual_vec_next);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<Vec3>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_residual_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<Vec3>().swap(per_atom_residual_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<Vec3>(
        std::move(buffer),
        std::type_index(typeid(
            TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) Vec3 next-residual time-series to tp dense buffer");
}


void TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<Vec3>(std::type_index(typeid(
            TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::WriteH5Group",
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
            "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::WriteH5Group",
            "TripeptideNeighborShieldingResult was not attached in any of "
            + std::to_string(source_present_per_frame_.size()) +
            " frames; skipping /trajectory/tripeptide_neighbor_residual_vec_next_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/tripeptide_neighbor_residual_vec_next_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("irrep_layout",  std::string("x,y,z"));
    grp.createAttribute("normalization", std::string("cartesian"));
    grp.createAttribute("parity",        std::string("1o"));
    grp.createAttribute("units",         std::string("angstrom"));

    // Flat (N, T, 3) via explicit component access. NaN values from
    // "no i+1 neighbour contribution this frame" propagate through
    // unchanged; also NaN-fill rows where source wasn't attached
    // (absent-not-faked).
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> flat(N * T * 3);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::size_t base = (i * T + t) * 3;
            if (t >= source_present_per_frame_.size()
                || source_present_per_frame_[t] == 0) {
                for (std::size_t k = 0; k < 3; ++k) flat[base + k] = kNaN;
                continue;
            }
            const Vec3& v = buffer->At(i, t);
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

    // Provenance mask: per-frame source-attached flags.
    grp.createDataSet("source_attached_per_frame", source_present_per_frame_);
}

}  // namespace nmr
