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
    r->per_bond_n_present_.assign(B, 0);
    r->per_bond_present_.assign(B, WelfordMoments{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// MopacResult.TopologyBondOrders() is parallel to protein.Bonds(),
// so per_bond_[bi] aligns with protein.BondAt(bi). MopacResult sets
// 0.0 (exact) for bonds it didn't report — the canonical "no
// observation" sentinel for a continuous QM observable. Sentinel-
// aware Welford per `feedback_conditional_welford_for_sentinels`:
//
//   - per_bond_ accumulates ONLY on frames where bo != 0.0
//     (per-bond divisor = per_bond_n_present_[bi]).
//   - per_bond_present_ indicator-Welford accumulates EVERY MOPAC-
//     attached frame: 1.0 if bo != 0.0, else 0.0
//     (per-bond divisor = source_attached_count_).
//
// Bond axis structural invariant: bond_orders.size() ==
// per_bond_.size() == protein.BondCount(). MopacResult sizes its
// topology_bond_orders_ vector from protein.BondCount() at
// extractor-side Compute (MopacResult.cpp). If they disagree at
// runtime, the protein's bond axis has changed across the lifecycle
// — fail loud rather than silently truncate.

void MopacBondOrderWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    // `did_update` is true ONLY when this frame produced a valid
    // Welford update (source attached AND bond axis invariant held).
    // Invariant-failure (size mismatch) is treated as "frame absent"
    // for accounting — mask=0 + source_attached_count_ unchanged —
    // so the H5 per-bond divisor stays consistent with the actual
    // sample count (codex 2026-05-21 finding 3).
    bool did_update = false;
    if (conf.HasResult<MopacResult>()) {
        const auto& mopac = conf.Result<MopacResult>();
        const auto& bond_orders = mopac.TopologyBondOrders();
        if (bond_orders.size() != per_bond_.size()) {
            OperationLog::Error(
                "MopacBondOrderWelfordTrajectoryResult::Compute",
                "bond axis size mismatch: MopacResult.TopologyBondOrders() "
                "= " + std::to_string(bond_orders.size()) +
                " vs TR per_bond_ = " + std::to_string(per_bond_.size()) +
                " — protein.BondCount() invariant violated; treating "
                "frame as source-absent (mask=0, count unchanged) so "
                "per-bond divisor stays consistent.");
        } else {
            const std::size_t B = per_bond_.size();
            const std::size_t n_total_new = source_attached_count_ + 1;
            for (std::size_t bi = 0; bi < B; ++bi) {
                const double bo = bond_orders[bi];
                // Indicator Welford — every MOPAC-attached frame
                // contributes a 0 or 1 sample per bond.
                const double indicator = (bo != 0.0) ? 1.0 : 0.0;
                WelfordUpdate(per_bond_present_[bi], indicator,
                              n_total_new, frame_idx);
                // Order Welford — only when bond was reported (bo != 0).
                if (bo != 0.0) {
                    const std::size_t n_pres_new = per_bond_n_present_[bi] + 1;
                    WelfordUpdate(per_bond_[bi], bo, n_pres_new, frame_idx);
                    per_bond_n_present_[bi] = n_pres_new;
                }
            }
            ++source_attached_count_;
            did_update = true;
        }
    }
    // Sparse cadence is normal — no per-absent-frame logging.

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(did_update ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void MopacBondOrderWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                     Trajectory& traj) {
    (void)tp; (void)traj;
    if (source_attached_count_ > 0) {
        const std::size_t B = per_bond_.size();
        for (std::size_t bi = 0; bi < B; ++bi) {
            // Order Welford finalized against per-bond present count
            // (NaN std when n_present == 0, i.e. MOPAC never reported
            // this bond on any frame).
            WelfordFinalize(per_bond_[bi], per_bond_n_present_[bi]);
            // Presence-fraction Welford finalized against protein-wide
            // attached count (every MOPAC-attached frame contributed an
            // indicator sample per bond).
            WelfordFinalize(per_bond_present_[bi], source_attached_count_);
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
        "Per-bond Welford rollup with sentinel-aware accumulation: "
        "the order Welford updates ONLY on frames where the bond was "
        "reported (bo != 0.0); the order_present_fraction "
        "indicator-Welford updates EVERY MOPAC-attached frame with a "
        "1.0/0.0 indicator. Mirrors the HydrationShellWelford "
        "ion_present_fraction pattern per "
        "feedback_conditional_welford_for_sentinels (R6 codex "
        "2026-05-18) — naive accumulation of the 0.0 sentinel "
        "biases the order mean toward 0 for intermittently-reported "
        "bonds (typical for MOZYME-merged sidechain interior bonds)."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- MopacResult attaches sparsely per the Mopac "
        "cadence (OperationRunner.cpp:138-148, Attach gated by "
        "!opts.skip_mopac AND non-null Compute return; NOT "
        "RequireConformationResult). Compute's HasResult<MopacResult>() "
        "gate skips per-bond Welford update + records mask=0 on absent "
        "frames. WriteH5Group skips the entire group when "
        "source_attached_count==0."));

    std::vector<double>        mean(B), std_(B), m2(B), min_(B), max_(B);
    std::vector<std::uint64_t> min_frame(B), max_frame(B), n_per_bond(B);
    std::vector<double>        pmean(B), pstd(B), pm2(B), pmin(B), pmax(B);
    std::vector<std::uint64_t> pmin_frame(B), pmax_frame(B), n_total_per_bond(B);
    for (std::size_t bi = 0; bi < B; ++bi) {
        const auto& w = per_bond_[bi];
        mean[bi]      = w.mean;
        std_[bi]      = w.std;
        m2[bi]        = w.m2;
        min_[bi]      = w.min;
        max_[bi]      = w.max;
        min_frame[bi] = static_cast<std::uint64_t>(w.min_frame);
        max_frame[bi] = static_cast<std::uint64_t>(w.max_frame);
        n_per_bond[bi] = static_cast<std::uint64_t>(per_bond_n_present_[bi]);

        const auto& p = per_bond_present_[bi];
        pmean[bi]      = p.mean;
        pstd[bi]       = p.std;
        pm2[bi]        = p.m2;
        pmin[bi]       = p.min;
        pmax[bi]       = p.max;
        pmin_frame[bi] = static_cast<std::uint64_t>(p.min_frame);
        pmax_frame[bi] = static_cast<std::uint64_t>(p.max_frame);
        n_total_per_bond[bi] = static_cast<std::uint64_t>(source_attached_count_);
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

    // Presence-fraction indicator Welford (mean ∈ [0, 1] = Pr(MOPAC
    // reports bond)). Per `feedback_conditional_welford_for_sentinels`.
    with_units(grp.createDataSet("order_present_fraction_mean",      pmean),      std::string("dimensionless"));
    with_units(grp.createDataSet("order_present_fraction_std",       pstd),       std::string("dimensionless"));
    with_units(grp.createDataSet("order_present_fraction_m2",        pm2),        std::string("dimensionless"));
    with_units(grp.createDataSet("order_present_fraction_min",       pmin),       std::string("dimensionless"));
    with_units(grp.createDataSet("order_present_fraction_max",       pmax),       std::string("dimensionless"));
    with_units(grp.createDataSet("order_present_fraction_min_frame", pmin_frame), std::string("frame_index"));
    with_units(grp.createDataSet("order_present_fraction_max_frame", pmax_frame), std::string("frame_index"));
    with_units(grp.createDataSet("n_total_per_bond",                 n_total_per_bond), std::string("frame_count"));

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
