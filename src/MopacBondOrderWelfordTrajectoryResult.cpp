#include "MopacBondOrderWelfordTrajectoryResult.h"

#include "MopacResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cstddef>
#include <string>

namespace nmr {

std::unique_ptr<MopacBondOrderWelfordTrajectoryResult>
MopacBondOrderWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<MopacBondOrderWelfordTrajectoryResult>();
    // Trajectory::Run Phase 2 (Seed) precedes Phase 3 (factory
    // invocation), so the Protein is finalised here — BondCount()
    // returns the real value, not zero.
    const std::size_t B = tp.ProteinRef().BondCount();
    r->per_bond_.assign(B, WelfordMoments{});
    r->per_bond_n_.assign(B, 0);
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// MopacResult.TopologyBondOrders() is parallel to protein.Bonds(),
// so per_bond_[bi] aligns with protein.BondAt(bi). MopacResult zeros
// out bonds it didn't report (NOT NaN), so we accumulate every bond
// with the same n_per_bond as the source_attached_count.

void MopacBondOrderWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_present = conf.HasResult<MopacResult>();
    if (source_present) {
        const auto& mopac = conf.Result<MopacResult>();
        const auto& bond_orders = mopac.TopologyBondOrders();
        const std::size_t B = std::min(bond_orders.size(), per_bond_.size());
        for (std::size_t bi = 0; bi < B; ++bi) {
            const double bo = bond_orders[bi];
            const std::size_t n_new = per_bond_n_[bi] + 1;
            WelfordUpdate(per_bond_[bi], bo, n_new, frame_idx);
            per_bond_n_[bi] = n_new;
        }
        ++source_attached_count_;
    }
    // Sparse cadence is normal — no per-absent-frame logging.

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void MopacBondOrderWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                     Trajectory& traj) {
    (void)tp; (void)traj;
    if (source_attached_count_ > 0) {
        const std::size_t B = per_bond_.size();
        for (std::size_t bi = 0; bi < B; ++bi) {
            WelfordFinalize(per_bond_[bi], per_bond_n_[bi]);
        }
    }
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "MopacBondOrderWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames; MOPAC ran on " +
        std::to_string(source_attached_count_) + " frames; " +
        std::to_string(per_bond_.size()) + " bonds");
}


// ── WriteH5Group ─────────────────────────────────────────────────

void MopacBondOrderWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp, HighFive::File& file) const {
    (void)tp;
    const std::size_t B = per_bond_.size();
    const std::size_t T = n_frames_;

    if (source_attached_count_ == 0) {
        OperationLog::Warn(
            "MopacBondOrderWelfordTrajectoryResult::WriteH5Group",
            "MOPAC ran on 0/" + std::to_string(T) +
            " frames; skipping /trajectory/mopac_bond_order_welford/ "
            "per canonical 'absent, not faked' discipline.");
        return;
    }

    auto grp = file.createGroup("/trajectory/mopac_bond_order_welford");

    grp.createAttribute("result_name",           Name());
    grp.createAttribute("n_bonds",               B);
    grp.createAttribute("n_frames",              T);
    grp.createAttribute("source_attached_count", source_attached_count_);
    grp.createAttribute("finalized",             finalized_);
    grp.createAttribute("units",                 std::string("dimensionless"));
    grp.createAttribute("bond_axis",             std::string("bonds.npy"));
    grp.createAttribute("source", std::string(
        "MopacResult.TopologyBondOrders() (Wiberg bond orders from "
        "PM7+MOZYME, parallel to protein.Bonds() == bonds.npy axis). "
        "Per-bond Welford rollup; emits canonical row mean + std + "
        "m2 + min/max + min_frame/max_frame + n_per_bond. "
        "Minimum-viable v0 — no delta variants. Bonds with order "
        "near zero (no MOPAC entry) accumulate as 0.0 per frame, NOT "
        "as NaN — MopacResult.TopologyBondOrders() returns 0.0 for "
        "missing entries."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- MopacResult attaches sparsely per the Mopac "
        "cadence (OperationRunner.cpp:142, TimedAttach not "
        "RequireConformationResult). Compute's HasResult<MopacResult>() "
        "gate skips per-bond Welford update + records mask=0 on absent "
        "frames. WriteH5Group skips the entire group when "
        "source_attached_count==0."));

    std::vector<double>        mean(B), std_(B), m2(B), min_(B), max_(B);
    std::vector<std::uint64_t> min_frame(B), max_frame(B), n_per_bond(B);
    for (std::size_t bi = 0; bi < B; ++bi) {
        const auto& w = per_bond_[bi];
        mean[bi]      = w.mean;
        std_[bi]      = w.std;
        m2[bi]        = w.m2;
        min_[bi]      = w.min;
        max_[bi]      = w.max;
        min_frame[bi] = static_cast<std::uint64_t>(w.min_frame);
        max_frame[bi] = static_cast<std::uint64_t>(w.max_frame);
        n_per_bond[bi] = static_cast<std::uint64_t>(per_bond_n_[bi]);
    }

    auto with_units = [](HighFive::DataSet ds, const std::string& u) {
        ds.createAttribute("units", u);
    };
    // Bond-order is dimensionless; m2 is dimensionless squared
    // (still dimensionless, but using the squared label per the
    // canonical convention to keep readers aware of the m2 vs base
    // scaling — codex H3 2026-05-20).
    with_units(grp.createDataSet("order_mean",      mean),      std::string("dimensionless"));
    with_units(grp.createDataSet("order_std",       std_),      std::string("dimensionless"));
    with_units(grp.createDataSet("order_m2",        m2),        std::string("dimensionless^2"));
    with_units(grp.createDataSet("order_min",       min_),      std::string("dimensionless"));
    with_units(grp.createDataSet("order_max",       max_),      std::string("dimensionless"));
    with_units(grp.createDataSet("order_min_frame", min_frame), std::string("frame_index"));
    with_units(grp.createDataSet("order_max_frame", max_frame), std::string("frame_index"));
    with_units(grp.createDataSet("n_per_bond",      n_per_bond), std::string("frame_count"));

    grp.createDataSet("frame_indices",            frame_indices_);
    grp.createDataSet("frame_times",              frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);

    OperationLog::Info(LogCalcOther,
        "MopacBondOrderWelfordTrajectoryResult::WriteH5Group",
        "wrote /trajectory/mopac_bond_order_welford with " +
        std::to_string(B) + " bonds (" +
        std::to_string(source_attached_count_) + "/" + std::to_string(T) +
        " MOPAC-attached frames)");
}

}  // namespace nmr
