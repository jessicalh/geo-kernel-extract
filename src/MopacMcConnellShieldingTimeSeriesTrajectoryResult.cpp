#include "MopacMcConnellShieldingTimeSeriesTrajectoryResult.h"

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
// Same sparse gate as TR7. Source field is NOT traceless (M_total
// has both symmetric and antisymmetric parts) — we emit the full
// 9-component SphericalTensor at WriteH5Group time.

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
    grp.createAttribute("irrep_layout", std::string(
        "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("0e+1o+2e"));
    grp.createAttribute("units",         std::string("Angstrom^-3"));
    grp.createAttribute("source", std::string(
        "MopacMcConnellResult.mopac_mc_shielding_contribution "
        "(SphericalTensor; bond-order-weighted bo·M/r³ kernel, NO "
        "Δχ × γ multiplication at extraction — bare kernel in Å^-3). "
        "M_total is NOT traceless: T0 = trace(M)/3 = bond-order-"
        "weighted McConnell f-scalar sum over ALL bond categories "
        "that enter the per-atom loop (PeptideCO, PeptideCN, "
        "BackboneOther, SidechainCO, SidechainOther, Aromatic — see "
        "MopacMcConnellResult.cpp:185-227). Note: the named "
        "ca.mopac_mc_*_sum scalars on ConformationAtom cover only "
        "{co_sum, cn_sum, sidechain_sum (SidechainCO only), "
        "aromatic_sum} — BackboneOther and SidechainOther bonds "
        "contribute to T0 here but NOT to any named per-category "
        "sum; Disulfide bonds (BondCategory::Disulfide) hit the "
        "switch default and contribute to neither (pre-existing "
        "issue across both Mopac and FF14SB McConnell, "
        "MopacMcConnellResult.cpp:185-227). T1 = antisymmetric "
        "McConnell pseudovector arising from the asymmetric "
        "9·cos_θ·d̂_a·b̂_b cross-coupling term (the -3·b̂_a·b̂_b "
        "and -3·d̂_a·d̂_b outer-product terms are symmetric in (a,b) "
        "and contribute only to T0+T2). T1 is a real geometric "
        "quantity, not numerical noise. T2 = symmetric traceless "
        "McConnell tensor (canonical bond-anisotropy contribution). "
        "All 9 components emitted per user direction 2026-05-21 "
        "'if not traceless write both' + the methods-accumulate "
        "principle. Per-category T2 totals T2_backbone / "
        "T2_sidechain / T2_aromatic ARE explicitly symmetrized at "
        "MopacMcConnellResult.cpp:252-262 but live on separate "
        "ConformationAtom fields, not in this TR."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- MopacMcConnellResult attaches sparsely per the "
        "Mopac cadence (OperationRunner.cpp:185, TimedAttach not "
        "RequireConformationResult). Compute's HasResult<MopacMcConnellResult>() "
        "gate emits NaN-fill + mask=0 on absent frames per canonical "
        "'absent, not faked'. WriteH5Group skips the entire group "
        "when source_attached_count==0."));

    // (N, T, 9) — full SphericalTensor (T0 + T1[3] + T2[5]).
    // Atom-major: [atom_0_frame_0_(T0,T1_x,T1_y,T1_z,T2_-2,T2_-1,T2_0,T2_+1,T2_+2),
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
