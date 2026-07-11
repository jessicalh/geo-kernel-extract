#include "MopacMcConnellShieldingTimeSeriesTrajectoryResult.h"

#include "McConnellResult.h"
#include "MopacMcConnellResult.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "ProteinConformation.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <limits>
#include <typeinfo>

namespace nmr {

std::unique_ptr<MopacMcConnellShieldingTimeSeriesTrajectoryResult>
MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        MopacMcConnellShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_.assign(tp.AtomCount(), std::vector<SphericalTensor>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Same sparse gate as TR7. The BO-channel D(r)Qhat response is a full
// rank-2 tensor, so we emit the full 9-component SphericalTensor at
// WriteH5Group time.

void MopacMcConnellShieldingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_present = conf.HasResult<MopacMcConnellResult>();
    SphericalTensor nan_tensor{};
    if (!source_present) {
        const double nan_d = std::numeric_limits<double>::quiet_NaN();
        nan_tensor.T0 = nan_d;
        nan_tensor.T1 = {nan_d, nan_d, nan_d};
        nan_tensor.T2 = {nan_d, nan_d, nan_d, nan_d, nan_d};
    }

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_[i].push_back(
            source_present ? conf.AtomAt(i).mopac_mc_shielding_contribution
                           : nan_tensor);
    }
    if (source_present) ++source_attached_count_;
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void MopacMcConnellShieldingTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_[i]);
        ++atoms_written;
    }
    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(
            MopacMcConnellShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "MopacMcConnellShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) + " frames) Mopac McConnell shielding "
        "(T0+T1+T2 = 9 components, source not traceless) to tp dense "
        "buffer (MOPAC ran on " +
        std::to_string(source_attached_count_) + " frames)");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// 9-component (N, T, 9) emission — source is NOT traceless, preserve
// all of T0+T1+T2 per user 2026-05-21 "if not traceless write both".
// Matches FF14SB McConnellShieldingTimeSeries layout exactly.

void MopacMcConnellShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    if (source_attached_count_ == 0) {
        OperationLog::Warn(
            "MopacMcConnellShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "MopacMcConnellResult attached on 0/" +
            std::to_string(n_frames_) + " frames; skipping "
            "/trajectory/mopac_mc_shielding_time_series/ per "
            "canonical 'absent, not faked' discipline.");
        return;
    }

    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(typeid(
            MopacMcConnellShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "MopacMcConnellShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/mopac_mc_shielding_time_series");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("source_attached_count",  source_attached_count_);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("tensor_basis",
        std::string(kMcConnellPackFull9TensorBasis));
    grp.createAttribute("tensor_component_order",
        std::string(kMcConnellPackFull9ComponentOrder));
    grp.createAttribute("tensor_frame",
        std::string(kMcConnellPackFull9TensorFrame));
    grp.createAttribute("tensor_parity", std::string("even"));
    grp.createAttribute("e3nn_export",
        std::string(kMcConnellPackFull9E3nnExport));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("units",         std::string("Angstrom^-3"));
    grp.createAttribute("source", std::string(
        "MopacMcConnellResult.mopac_mc_shielding_contribution "
        "(SphericalTensor; MOPAC Wiberg-weighted BO channel from the "
        "unified McConnellResult source model A=D(r)Qhat, emitted "
        "unscaled in Angstrom^-3. Delta-chi, sign, and unit conversion "
        "are learned downstream. T0 is trace(DQhat)/3, equal to the PCS "
        "scalar n^T Qhat n/r^3 for each traceless unit source shape. "
        "T1 is the even antisymmetric pseudovector from non-commuting D "
        "and Qhat. T2 is the symmetric-traceless branch in the same "
        "SphericalTensor convention. Aromatic McConnell is zeroed "
        "unconditionally because BS/HM always compute aromatic "
        "ring-current contributions; missing or below-floor bond order "
        "zeros only this BO channel, not the fixed channel."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- MopacMcConnellResult attaches sparsely per the "
        "Mopac cadence (OperationRunner.cpp:185, TimedAttach not "
        "RequireConformationResult). Compute's HasResult<MopacMcConnellResult>() "
        "gate emits NaN-fill + mask=0 on absent frames per canonical "
        "'absent, not faked'. WriteH5Group skips the entire group "
        "when source_attached_count==0."));

    // (N, T, 9) — full SphericalTensor (T0 + T1[3] + T2[5]).
    // Atom-major: [atom_0_frame_0_(T0,T1_x,T1_y,T1_z,T2_m-2..+2),
    //              atom_0_frame_1_..., atom_1_frame_0_..., ...].
    std::vector<double> flat(N * T * 9);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const SphericalTensor& st = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 9;
            st.PackFull9(&flat[base]);
        }
    }
    std::vector<std::size_t> dims = {N, T, std::size_t(9)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices",            frame_indices_);
    grp.createDataSet("frame_times",              frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);
}

}  // namespace nmr
