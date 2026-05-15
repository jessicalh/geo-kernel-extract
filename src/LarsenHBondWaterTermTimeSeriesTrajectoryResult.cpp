#include "LarsenHBondWaterTermTimeSeriesTrajectoryResult.h"
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


std::unique_ptr<LarsenHBondWaterTermTimeSeriesTrajectoryResult>
LarsenHBondWaterTermTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        LarsenHBondWaterTermTimeSeriesTrajectoryResult>();
    r->per_atom_water_term_.assign(tp.AtomCount(), std::vector<double>{});
    return r;
}


void LarsenHBondWaterTermTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    // "Absent, not faked" provenance: record whether the source
    // calculator (LarsenHBondShieldingResult) attached this frame.
    // When absent, the in-memory field is zero-default — we capture
    // that here but NaN-fill at WriteH5Group time so downstream
    // readers can distinguish "no measurement" from "real
    // measurement = 0."
    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<LarsenHBondShieldingResult>();
    source_present_per_frame_.push_back(source_attached ? 1u : 0u);

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_water_term_[i].push_back(
            conf.AtomAt(i).larsen_hbond_water_term);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// Bounds-check Finalize idempotency.
void LarsenHBondWaterTermTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<double>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_water_term_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<double>().swap(per_atom_water_term_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<double>(
        std::move(buffer),
        std::type_index(typeid(
            LarsenHBondWaterTermTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "LarsenHBondWaterTermTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) double water-term time-series to tp dense buffer");
}


// (N, T) flat double emission. Scalar L0 contribution per atom per
// frame. No inner axis. Same flat (N, T) shape as
// BsT0AutocorrelationTrajectoryResult but with a time axis instead of
// lag axis.

void LarsenHBondWaterTermTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<double>(std::type_index(typeid(
            LarsenHBondWaterTermTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "LarsenHBondWaterTermTimeSeriesTrajectoryResult::WriteH5Group",
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
            "LarsenHBondWaterTermTimeSeriesTrajectoryResult::WriteH5Group",
            "LarsenHBondShieldingResult was not attached in any of "
            + std::to_string(source_present_per_frame_.size()) +
            " frames; skipping /trajectory/larsen_hbond_water_term_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/larsen_hbond_water_term_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    // L0 scalar metadata for e3nn-consumable downstream.
    grp.createAttribute("irrep_layout", std::string("T0"));
    grp.createAttribute("parity",       std::string("0e"));
    grp.createAttribute("units",        std::string("ppm"));
    grp.createAttribute("description",
        std::string("Larsen Δσ_w water term — 2.07 ppm isotropic on amide H "
                    "atoms with zero H-bond pair contributions this frame; "
                    "zero on all other atoms (and on HN with ≥1 pair)."));

    // Flat (N, T) double. Same component-access discipline as the
    // Vec3/SphericalTensor TRs even though the cell type is already
    // a scalar — we go through DenseBuffer::At() so the layout
    // contract stays explicit at the emission boundary. NaN-fill rows
    // where source wasn't attached so downstream readers distinguish
    // "no measurement" from "measurement was zero."
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> flat(N * T);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            if (t >= source_present_per_frame_.size()
                || source_present_per_frame_[t] == 0) {
                flat[i * T + t] = kNaN;
                continue;
            }
            flat[i * T + t] = buffer->At(i, t);
        }
    }

    std::vector<std::size_t> dims = {N, T};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("water_term", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);

    // Provenance mask: per-frame source-attached flags.
    grp.createDataSet("source_attached_per_frame", source_present_per_frame_);
}

}  // namespace nmr
