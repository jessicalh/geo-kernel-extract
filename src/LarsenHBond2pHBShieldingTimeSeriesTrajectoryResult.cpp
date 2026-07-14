#include "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult.h"
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
#include <cstdint>
#include <typeinfo>

namespace nmr {

std::unique_ptr<LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult>
LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_shielding_.assign(tp.AtomCount(),
                                  std::vector<SphericalTensor>{});
    return r;
}

void LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    // "Absent, not faked" provenance: record whether the source
    // calculator (LarsenHBondShieldingResult) is present for this frame
    // through actual attachment or the test-only override.
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
            conf.AtomAt(i).larsen_hbond_2pHB_spherical);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}

void LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::Finalize(
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
            LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) 2pHB shielding time-series to tp dense buffer");
}

void LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(typeid(
            LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    // "Absent, not faked" — if no frame had the source-present flag,
    // skip emission. Group existence ⇒ source ran in ≥1 frame, or a
    // synthetic test forced presence. Downstream readers MUST tolerate
    // group absence for conditionally-attached-source TRs.
    std::size_t source_present_count = 0;
    for (auto v : source_present_per_frame_)
        if (v) ++source_present_count;
    if (source_present_count == 0) {
        OperationLog::Warn(
            "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "LarsenHBondShieldingResult was not attached in any of "
            + std::to_string(source_present_per_frame_.size()) +
            " frames; skipping /trajectory/larsen_hbond_2pHB_shielding_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/larsen_hbond_2pHB_shielding_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("source_attached_count",
        static_cast<std::uint64_t>(source_present_count));
    grp.createAttribute("source_attached_policy",
        std::string("conditional_larsen_grid_source"));
    grp.createAttribute("atom_axis", std::string("protein_atom_index"));
    grp.createAttribute("frame_axis", std::string("trajectory_frame_row"));

    grp.createAttribute("irrep_layout",
        std::string("PackFull9: [T0, T1_cartesian_xyz, T2_real_tesseral_m-2..m+2]"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("mixed"));
    grp.createAttribute("coordinate_frame",
        std::string("conformation_cartesian_xyz"));
    grp.createAttribute("tensor_basis",
        std::string("project_native_full9_spherical_tensor_v1"));
    grp.createAttribute("tensor_component_order", std::string(
        "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("tensor_frame",
        std::string("conformation_cartesian_xyz"));
    grp.createAttribute("tensor_t1_semantics", std::string(
        "Cartesian Levi-Civita dual x,y,z (not real-Y1m): "
        "a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2); "
        "axial a'=det(R) R a; generically nonzero"));
    grp.createAttribute("tensor_t1_structural_zero", false);
    grp.createAttribute("tensor_structural_zero_components",
        std::string("none"));
    grp.createAttribute("e3nn_export", std::string(
        "explicit project-basis to e3nn conversion required before use"));
    grp.createAttribute("normalization_scope", std::string(
        "xyz tensor payload: T2 uses isometric real-tesseral "
        "normalization; T1 uses the tensor_t1_semantics convention"));
    grp.createAttribute("transformation", std::string(
        "even_rank2 under proper rotations: T'=R T R^T; signed-rho DFT-grid "
        "lookup is chirality-conditioned and has no improper-transform contract"));
    grp.createAttribute("units",         std::string("ppm"));

    // Flat (N, T, 9). NaN-fill rows where the source-present flag is
    // 0 — readers use isfinite/isnan to distinguish "no measurement"
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
            st.PackFull9(&flat[base]);
        }
    }

    std::vector<std::size_t> dims = {N, T, std::size_t(9)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);

    // Provenance mask: per-frame source-present flags.
    grp.createDataSet("source_attached_per_frame", source_present_per_frame_);
}

}  // namespace nmr
