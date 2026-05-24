#include "LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult.h"
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

std::unique_ptr<LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult>
LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_shielding_.assign(tp.AtomCount(),
                                  std::vector<SphericalTensor>{});
    return r;
}

void LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::Compute(
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
        per_atom_shielding_[i].push_back(
            conf.AtomAt(i).larsen_hbond_2pHaB_spherical);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}

void LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer =
        std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_shielding_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_shielding_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(
            LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) 2pHaB shielding time-series to tp dense buffer");
}

void LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(typeid(
            LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::WriteH5Group",
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
            "LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "LarsenHBondShieldingResult was not attached in any of "
            + std::to_string(source_present_per_frame_.size()) +
            " frames; skipping /trajectory/larsen_hbond_2pHaB_shielding_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/larsen_hbond_2pHaB_shielding_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("irrep_layout",
        std::string("T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("0e+1o+2e"));
    grp.createAttribute("units",         std::string("ppm"));

    // Flat (N, T, 9). NaN-fill rows where source wasn't attached —
    // readers use isfinite/isnan to distinguish "no measurement"
    // from "measurement was zero."
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> flat(N * T * 9);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::size_t base = (i * T + t) * 9;
            if (t >= source_present_per_frame_.size()
                || source_present_per_frame_[t] == 0) {
                for (std::size_t k = 0; k < 9; ++k) flat[base + k] = kNaN;
                continue;
            }
            const SphericalTensor& st = buffer->At(i, t);
            flat[base + 0] = st.T0;
            flat[base + 1] = st.T1[0];
            flat[base + 2] = st.T1[1];
            flat[base + 3] = st.T1[2];
            flat[base + 4] = st.T2[0];
            flat[base + 5] = st.T2[1];
            flat[base + 6] = st.T2[2];
            flat[base + 7] = st.T2[3];
            flat[base + 8] = st.T2[4];
        }
    }

    std::vector<std::size_t> dims = {N, T, std::size_t(9)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);

    // Provenance mask: per-frame source-attached flags.
    grp.createDataSet("source_attached_per_frame", source_present_per_frame_);
}

}  // namespace nmr
