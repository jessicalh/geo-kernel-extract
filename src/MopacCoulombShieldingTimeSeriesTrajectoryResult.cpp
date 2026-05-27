#include "MopacCoulombShieldingTimeSeriesTrajectoryResult.h"

#include "MopacCoulombResult.h"
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

std::unique_ptr<MopacCoulombShieldingTimeSeriesTrajectoryResult>
MopacCoulombShieldingTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        MopacCoulombShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_.assign(tp.AtomCount(), std::vector<SphericalTensor>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Sparse cadence gate identical to TR5/TR6: HasResult<MopacCoulombResult>
// skip + NaN-fill on absent frames. Source field is T2-only per the
// MopacCoulombResult source comment ("Pure T2 (EFG is traceless)").

void MopacCoulombShieldingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_present = conf.HasResult<MopacCoulombResult>();
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
            source_present ? conf.AtomAt(i).mopac_coulomb_shielding_contribution
                           : nan_tensor);
    }
    if (source_present) ++source_attached_count_;
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void MopacCoulombShieldingTimeSeriesTrajectoryResult::Finalize(
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
            MopacCoulombShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "MopacCoulombShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) + " frames) Mopac Coulomb shielding "
        "T2 SphericalTensor to tp dense buffer (MOPAC ran on " +
        std::to_string(source_attached_count_) + " frames)");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// T2-only (N, T, 5) emission — source field is genuinely T2. Group
// skipped entirely when source never attached (canonical "absent, not
// faked").

void MopacCoulombShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    if (source_attached_count_ == 0) {
        OperationLog::Warn(
            "MopacCoulombShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "MopacCoulombResult attached on 0/" +
            std::to_string(n_frames_) + " frames; skipping "
            "/trajectory/mopac_coulomb_shielding_time_series/ per "
            "canonical 'absent, not faked' discipline.");
        return;
    }

    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(typeid(
            MopacCoulombShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "MopacCoulombShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/mopac_coulomb_shielding_time_series");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("source_attached_count",  source_attached_count_);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("irrep_layout",
        std::string("T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("2e"));
    grp.createAttribute("units",         std::string("V/Å^2"));
    grp.createAttribute("source", std::string(
        "MopacCoulombResult.mopac_coulomb_shielding_contribution "
        "(SphericalTensor; T2-only per source comment 'Pure T2 (EFG "
        "is traceless). gamma converts this to shielding.' at "
        "MopacCoulombResult.cpp:252-254). The stored value is the "
        "EFG_total kernel in V/Å² (Hessian of φ from MOPAC Mulliken "
        "charges, traceless-projected at MopacCoulombResult.cpp:175-178). "
        "NO γ multiplication is applied at extraction — despite the "
        "historical field name 'shielding_contribution', the stored "
        "quantity is the bare EFG kernel. Shielding is recovered at "
        "calibration time by multiplying by per-element γ."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- MopacCoulombResult attaches sparsely per the "
        "Mopac cadence (OperationRunner.cpp:183, TimedAttach not "
        "RequireConformationResult). Compute's HasResult<MopacCoulombResult>() "
        "gate emits NaN-fill + mask=0 on absent frames per canonical "
        "'absent, not faked'. WriteH5Group skips the entire group "
        "when source_attached_count==0."));

    // (N, T, 5) — T2 components only. Atom-major flat layout.
    std::vector<double> flat(N * T * 5);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const SphericalTensor& st = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 5;
            st.PackT2(&flat[base]);
        }
    }
    std::vector<std::size_t> dims = {N, T, std::size_t(5)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("t2", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices",            frame_indices_);
    grp.createDataSet("frame_times",              frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);
}

}  // namespace nmr
