#include "PositionsTimeSeriesTrajectoryResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <typeinfo>

namespace nmr {


std::unique_ptr<PositionsTimeSeriesTrajectoryResult>
PositionsTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<PositionsTimeSeriesTrajectoryResult>();
    r->per_atom_positions_.assign(tp.AtomCount(), std::vector<Vec3>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Append this frame's position to each atom's growing buffer. Record
// frame_idx + time_ps in the frame-axis metadata so WriteH5Group can
// emit the index/time datasets alongside the xyz dataset without
// reaching into Trajectory's frame record (which is the authoritative
// Run-level record but may be strided; this Result records the
// frames it actually saw).

void PositionsTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_positions_[i].push_back(conf.PositionAt(i));
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Flatten per-atom growing buffers into an atom-major
// DenseBuffer<Vec3>, transfer ownership to TrajectoryProtein. The
// per-atom buffers are cleared after copy — post-Finalize the data
// lives on tp's dense-buffer pool.

void PositionsTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                   Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<Vec3>>(N, n_frames_);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_positions_[i];
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        // Release memory on the growing side; the dense buffer owns it now.
        std::vector<Vec3>().swap(per_atom_positions_[i]);
    }

    tp.AdoptDenseBuffer<Vec3>(
        std::move(buffer),
        std::type_index(typeid(PositionsTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "PositionsTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) + " frames) positions to tp dense buffer");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// /trajectory/positions/
//   xyz              (N, T, 3) float64  — atom-major positions in Å
//   frame_indices    (T,)      uint64   — original XTC indices
//   frame_times      (T,)      float64  — simulation times in ps
//   result_name      attr      string
//   n_frames, n_atoms, finalized   attrs
//
// The DenseBuffer<Vec3> is contiguous memory with Vec3 = 3 sequential
// doubles (Eigen Matrix<double,3,1> with no padding), so the raw
// buffer IS the (N*T*3) double array in atom-major layout. We write
// it to HighFive with an explicit 3D DataSpace.

void PositionsTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<Vec3>(std::type_index(
            typeid(PositionsTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "PositionsTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/positions");
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    // Emit xyz with explicit 3D shape. DenseBuffer<Vec3>::RawData()
    // returns Vec3* pointing at atom-major contiguous (N*T) Vec3s;
    // reinterpret as (N*T*3) doubles.
    const double* raw = reinterpret_cast<const double*>(buffer->RawData());
    std::vector<std::size_t> dims = {N, T, std::size_t(3)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(raw);

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);
}

}  // namespace nmr
